//! ENA accession modeling, file-report resolution, and retrying HTTP byte-stream readers.

use std::{
    fmt,
    io::{self, Read as _},
    ops::Range,
    str::FromStr,
    thread,
    time::Duration,
};

use color_eyre::eyre::{Result, WrapErr, bail, ensure, eyre};
use md5::{Digest as _, Md5};
use reqwest::{
    StatusCode,
    blocking::{Client, Response},
    header::{
        ACCEPT_ENCODING, CONTENT_ENCODING, CONTENT_LENGTH, CONTENT_RANGE, CONTENT_TYPE,
        HeaderValue, RANGE,
    },
};
use url::Url;

use crate::record::MateSide;

const ENA_FILEREPORT_BASE_URL: &str = "https://www.ebi.ac.uk/ena/portal/api/filereport";

#[derive(Clone, Debug, Eq, PartialEq, Hash)]
/// ENA run accession constrained to the run-level prefixes accepted by this tool.
pub struct Accession(String);

impl Accession {
    /// Construct a validated ENA run accession.
    ///
    /// # Errors
    ///
    /// Returns an error when the value does not start with `SRR`, `ERR`, or `DRR` followed by
    /// digits.
    pub fn new(value: impl Into<String>) -> Result<Self> {
        let value = value.into();
        validate_run_accession(&value)?;
        Ok(Self(value))
    }

    /// Borrow the accession as a string slice.
    pub fn as_str(&self) -> &str {
        &self.0
    }
}

impl fmt::Display for Accession {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

impl FromStr for Accession {
    type Err = color_eyre::Report;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        Self::new(value)
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
/// One resolved ENA compressed FASTQ with catalogue integrity metadata.
pub(crate) struct EnaFastq {
    accession: Accession,
    mate: Option<MateSide>,
    url: Url,
    expected_bytes: u64,
    expected_md5: [u8; 16],
}

#[derive(Clone, Debug, Eq, PartialEq)]
/// A supported single- or paired-end ENA input after metadata validation.
pub(crate) enum EnaInput {
    /// One single-end compressed FASTQ.
    Single(EnaFastq),
    /// Two distinct compressed FASTQ mates.
    Paired { left: EnaFastq, right: EnaFastq },
}

impl EnaInput {
    fn from_filereport(accession: &Accession, body: &str) -> Result<Self> {
        let mut lines = body.lines();
        let header = lines
            .next()
            .ok_or_else(|| eyre!("ENA filereport response was empty for accession {accession}"))?;
        let row = lines.next().ok_or_else(|| {
            eyre!("ENA filereport did not return a data row for accession {accession}")
        })?;
        ensure!(
            lines.next().is_none(),
            "ENA filereport returned more than one run row for accession {accession}"
        );

        let names = header.split('\t').collect::<Vec<_>>();
        let values = row.split('\t').collect::<Vec<_>>();
        ensure!(
            names.len() == values.len(),
            "ENA filereport row shape did not match its header for accession {accession}"
        );

        let returned_accession = required_filereport_field(&names, &values, "run_accession")?;
        let layout = required_filereport_field(&names, &values, "library_layout")?;
        let urls = required_filereport_field(&names, &values, "fastq_ftp")?;
        let byte_counts = required_filereport_field(&names, &values, "fastq_bytes")?;
        let md5s = required_filereport_field(&names, &values, "fastq_md5")?;

        ensure!(
            returned_accession == accession.as_str(),
            "ENA filereport returned accession {returned_accession} while resolving {accession}"
        );

        let urls = urls.split(';').collect::<Vec<_>>();
        let byte_counts = byte_counts.split(';').collect::<Vec<_>>();
        let md5s = md5s.split(';').collect::<Vec<_>>();
        ensure!(
            urls.len() == byte_counts.len() && urls.len() == md5s.len(),
            "ENA FASTQ URL, byte-count, and MD5 cardinalities differ for accession {accession}"
        );

        let fastqs = urls
            .into_iter()
            .zip(byte_counts)
            .zip(md5s)
            .map(|((url, bytes), md5)| {
                let expected_bytes = bytes
                    .parse()
                    .wrap_err_with(|| format!("ENA fastq_bytes value was not numeric: {bytes}"))?;
                Ok(EnaFastq {
                    accession: accession.clone(),
                    mate: None,
                    url: parse_ena_fastq_url(url)?,
                    expected_bytes,
                    expected_md5: parse_md5_bytes(md5)?,
                })
            })
            .collect::<Result<Vec<_>>>()?;

        match layout {
            "SINGLE" => {
                let [single]: [EnaFastq; 1] = fastqs.try_into().map_err(|files: Vec<_>| {
                    eyre!(
                        "ENA SINGLE layout returned {} FASTQ files for accession {accession}",
                        files.len()
                    )
                })?;
                Ok(Self::Single(single))
            }
            "PAIRED" => {
                let [mut left, mut right]: [EnaFastq; 2] =
                    fastqs.try_into().map_err(|files: Vec<_>| {
                        eyre!(
                            "ENA PAIRED layout returned {} FASTQ files for accession {accession}",
                            files.len()
                        )
                    })?;
                ensure!(
                    left.url != right.url,
                    "ENA paired FASTQ URLs must be distinct for accession {accession}"
                );
                left.mate = Some(MateSide::Left);
                right.mate = Some(MateSide::Right);
                Ok(Self::Paired { left, right })
            }
            other => bail!("ENA returned unsupported library_layout: {other}"),
        }
    }
}

fn required_filereport_field<'a>(
    names: &[&str],
    values: &[&'a str],
    required: &str,
) -> Result<&'a str> {
    let value = names
        .iter()
        .position(|name| *name == required)
        .and_then(|index| values.get(index))
        .copied()
        .ok_or_else(|| eyre!("ENA filereport response did not include {required}"))?;
    ensure!(
        !value.is_empty(),
        "ENA filereport field {required} was empty"
    );
    Ok(value)
}

fn parse_ena_fastq_url(value: &str) -> Result<Url> {
    ensure!(!value.is_empty(), "ENA FASTQ URL must not be empty");
    let url = match Url::parse(value) {
        Ok(url) => url,
        Err(url::ParseError::RelativeUrlWithoutBase) => Url::parse(&format!("https://{value}"))?,
        Err(error) => return Err(error.into()),
    };
    ensure!(url.scheme() == "https", "ENA FASTQ URL must use https");
    ensure!(
        url.path().ends_with(".fastq.gz"),
        "ENA FASTQ URL path must end with .fastq.gz"
    );
    Ok(url)
}

fn parse_md5_bytes(value: &str) -> Result<[u8; 16]> {
    ensure!(
        value.len() == 32,
        "ENA fastq_md5 value must contain 32 hexadecimal characters"
    );
    let mut digest = [0_u8; 16];
    for (output, pair) in digest.iter_mut().zip(value.as_bytes().chunks_exact(2)) {
        let pair = std::str::from_utf8(pair)?;
        *output = u8::from_str_radix(pair, 16)
            .wrap_err_with(|| format!("ENA fastq_md5 value was not hexadecimal: {value}"))?;
    }
    Ok(digest)
}

#[derive(Clone, Debug)]
/// Blocking ENA client for metadata lookup and resumable HTTP stream construction.
pub struct EnaClient {
    http: Client,
    max_retries: u32,
    initial_backoff: Duration,
    max_backoff: Duration,
}

impl Default for EnaClient {
    fn default() -> Self {
        Self::new().expect("default ENA client configuration should be valid")
    }
}

impl EnaClient {
    /// Construct an ENA client with a default retry budget and backoff policy.
    ///
    /// # Errors
    ///
    /// Returns an error when the underlying blocking HTTP client cannot be constructed.
    pub fn new() -> Result<Self> {
        let http = Client::builder().build()?;
        Ok(Self {
            http,
            max_retries: 4,
            initial_backoff: Duration::from_millis(250),
            max_backoff: Duration::from_secs(10),
        })
    }

    /// Resolve an ENA accession into a validated single- or paired-end input.
    ///
    /// # Errors
    ///
    /// Returns an error when the filereport request fails, returns malformed data, or describes
    /// an unsupported FASTQ arrangement.
    pub(crate) fn resolve(&self, accession: &Accession) -> Result<EnaInput> {
        let response = self
            .http
            .get(ENA_FILEREPORT_BASE_URL)
            .query(&[
                ("accession", accession.as_str()),
                ("result", "read_run"),
                (
                    "fields",
                    "run_accession,library_layout,fastq_ftp,fastq_bytes,fastq_md5",
                ),
            ])
            .send()
            .wrap_err_with(|| {
                format!(
                    "ENA filereport request failed before a response was received for accession {accession}\n\
                     url: {ENA_FILEREPORT_BASE_URL}\n\
                     help: check network access to ebi.ac.uk from this runtime"
                )
            })?;

        let body = response_text_or_error(response)?;
        EnaInput::from_filereport(accession, &body).wrap_err_with(|| {
            format!(
                "ENA filereport could not be resolved for accession {accession}\n\
                 help: nuclease requires one single-end FASTQ or two distinct paired-end FASTQs with aligned byte-count and MD5 metadata"
            )
        })
    }

    /// Open and preflight a resolved ENA FASTQ stream.
    ///
    /// # Errors
    ///
    /// Returns an error when the initial response cannot be opened or does not match the
    /// catalogue metadata and gzip shape.
    pub(crate) fn stream(&self, fastq: EnaFastq) -> Result<EnaStream> {
        Ok(EnaStream::open(self.clone(), fastq)?)
    }
}

/// Blocking ENA FASTQ stream with ranged reconnects and compressed-object verification.
pub(crate) struct EnaStream {
    client: EnaClient,
    fastq: EnaFastq,
    body: Option<ActiveResponse>,
    compressed_offset: u64,
    digest: Md5,
    prefix: [u8; 2],
    prefix_delivered: usize,
    retries_remaining: u32,
    current_backoff: Duration,
    integrity_verified: bool,
    terminal_failure: bool,
    saw_eof: bool,
}

impl EnaStream {
    fn open(client: EnaClient, fastq: EnaFastq) -> io::Result<Self> {
        let current_backoff = client.initial_backoff;
        let retries_remaining = client.max_retries;
        let mut stream = Self {
            client,
            fastq,
            body: None,
            compressed_offset: 0,
            digest: Md5::new(),
            prefix: [0; 2],
            prefix_delivered: 0,
            retries_remaining,
            current_backoff,
            integrity_verified: false,
            terminal_failure: false,
            saw_eof: false,
        };

        loop {
            if let Err(error) = stream.ensure_connected() {
                stream.retry_or_return(error)?;
                continue;
            }

            let Some(body) = stream.body.as_mut() else {
                return Err(io::Error::other(
                    "missing ENA response body during preflight",
                ));
            };
            match body.read_exact(&mut stream.prefix) {
                Ok(()) => break,
                Err(error) => stream.retry_or_return(error)?,
            }
        }

        if stream.prefix != [0x1f, 0x8b] {
            return Err(invalid_data(format!(
                "ENA FASTQ for {} did not begin with gzip magic: {:02x} {:02x}",
                stream.fastq.accession, stream.prefix[0], stream.prefix[1]
            )));
        }

        Ok(stream)
    }

    fn ensure_connected(&mut self) -> io::Result<()> {
        if self.body.is_some() {
            return Ok(());
        }

        let mut request = self
            .client
            .http
            .get(self.fastq.url.clone())
            .header(ACCEPT_ENCODING, "identity");
        if self.compressed_offset > 0 {
            request = request.header(RANGE, format!("bytes={}-", self.compressed_offset));
        }

        let response = request.send().map_err(io::Error::other)?;
        let body = if self.compressed_offset == 0 {
            ActiveResponse::initial(response, &self.fastq)?
        } else {
            ActiveResponse::resumed(response, &self.fastq, self.compressed_offset)?
        };
        self.body = Some(body);
        Ok(())
    }

    fn retry_or_return(&mut self, error: io::Error) -> io::Result<()> {
        if !Self::should_retry_io_error(&error) {
            return Err(error);
        }
        if self.retries_remaining == 0 {
            return Err(io::Error::other("ENA stream retry budget exhausted"));
        }

        self.retries_remaining -= 1;
        self.body = None;
        thread::sleep(self.current_backoff);
        self.current_backoff = self
            .current_backoff
            .saturating_mul(2)
            .min(self.client.max_backoff);
        Ok(())
    }

    fn should_retry_io_error(error: &io::Error) -> bool {
        matches!(
            error.kind(),
            io::ErrorKind::BrokenPipe
                | io::ErrorKind::ConnectionAborted
                | io::ErrorKind::ConnectionReset
                | io::ErrorKind::Interrupted
                | io::ErrorKind::Other
                | io::ErrorKind::TimedOut
                | io::ErrorKind::UnexpectedEof
        )
    }

    fn reset_backoff_after_success(&mut self) {
        self.retries_remaining = self.client.max_retries;
        self.current_backoff = self.client.initial_backoff;
    }

    fn account(&mut self, bytes: &[u8]) -> io::Result<()> {
        let byte_count = u64::try_from(bytes.len())
            .map_err(|_| invalid_data("ENA read length did not fit in u64"))?;
        let next_offset = self
            .compressed_offset
            .checked_add(byte_count)
            .ok_or_else(|| invalid_data("ENA compressed byte offset overflowed"))?;
        if next_offset > self.fastq.expected_bytes {
            self.terminal_failure = true;
            return Err(invalid_data(format!(
                "ENA stream for {} exceeded catalogue size of {} bytes",
                self.fastq.accession, self.fastq.expected_bytes
            )));
        }

        self.digest.update(bytes);
        if next_offset == self.fastq.expected_bytes {
            let observed = self.digest.clone().finalize();
            if observed.as_slice() != self.fastq.expected_md5 {
                self.terminal_failure = true;
                return Err(invalid_data(format!(
                    "ENA compressed MD5 mismatch for {} after {next_offset} bytes",
                    self.fastq.accession
                )));
            }
            self.integrity_verified = true;
        }

        self.compressed_offset = next_offset;
        self.reset_backoff_after_success();
        Ok(())
    }
}

impl io::Read for EnaStream {
    fn read(&mut self, output: &mut [u8]) -> io::Result<usize> {
        if output.is_empty() {
            return Ok(0);
        }
        if self.saw_eof {
            return Ok(0);
        }
        if self.terminal_failure {
            return Err(invalid_data(
                "ENA stream previously failed integrity validation",
            ));
        }

        if self.prefix_delivered < self.prefix.len() {
            let pending = &self.prefix[self.prefix_delivered..];
            let count = pending.len().min(output.len());
            output[..count].copy_from_slice(&pending[..count]);
            self.account(&output[..count])?;
            self.prefix_delivered += count;
            return Ok(count);
        }

        loop {
            if let Err(error) = self.ensure_connected() {
                self.retry_or_return(error)?;
                continue;
            }

            let Some(body) = self.body.as_mut() else {
                return Err(io::Error::other("missing ENA response body"));
            };
            match body.read(output) {
                Ok(0) if self.integrity_verified => {
                    self.saw_eof = true;
                    return Ok(0);
                }
                Ok(0) => self.retry_or_return(io::Error::from(io::ErrorKind::UnexpectedEof))?,
                Ok(count) => {
                    self.account(&output[..count])?;
                    return Ok(count);
                }
                Err(error) => self.retry_or_return(error)?,
            }
        }
    }
}

#[derive(Debug)]
struct ActiveResponse {
    response: Response,
    promised: Range<u64>,
    bytes_read: u64,
}

impl ActiveResponse {
    fn initial(response: Response, fastq: &EnaFastq) -> io::Result<Self> {
        if response.status() != StatusCode::OK {
            return Err(io::Error::other(format!(
                "initial ENA stream request for {} failed with status {}",
                fastq.accession,
                response.status()
            )));
        }
        validate_identity_content_encoding(&response)?;

        let content_length = response
            .headers()
            .get(CONTENT_LENGTH)
            .ok_or_else(|| {
                invalid_data("initial ENA stream response did not include Content-Length")
            })
            .and_then(|value| parse_u64_header(value, "Content-Length"))?;
        if content_length != fastq.expected_bytes {
            return Err(invalid_data(format!(
                "ENA Content-Length was {content_length}, catalogue expected {} for {}",
                fastq.expected_bytes, fastq.accession
            )));
        }

        if response.url().path().ends_with('/') {
            return Err(invalid_data(format!(
                "ENA FASTQ for {} redirected to a directory URL: {}",
                fastq.accession,
                response.url()
            )));
        }
        if response
            .headers()
            .get(CONTENT_TYPE)
            .and_then(|value| value.to_str().ok())
            .and_then(|value| value.split(';').next())
            .is_some_and(|media_type| media_type.trim().eq_ignore_ascii_case("text/html"))
        {
            return Err(invalid_data(format!(
                "ENA FASTQ response for {} was HTML",
                fastq.accession
            )));
        }

        Ok(Self::new(response, 0..fastq.expected_bytes))
    }

    fn resumed(response: Response, fastq: &EnaFastq, expected_offset: u64) -> io::Result<Self> {
        if response.status() != StatusCode::PARTIAL_CONTENT {
            let message = format!(
                "ranged ENA stream request for {} failed with status {}",
                fastq.accession,
                response.status()
            );
            return if response.status().is_success() {
                Err(invalid_data(message))
            } else {
                Err(io::Error::other(message))
            };
        }
        validate_identity_content_encoding(&response)?;

        let content_range = response
            .headers()
            .get(CONTENT_RANGE)
            .ok_or_else(|| {
                invalid_data("ranged ENA stream response did not include Content-Range")
            })?
            .to_str()
            .map_err(|_| invalid_data("Content-Range header was not valid UTF-8"))?;
        let parsed = parse_content_range(content_range)?;
        if parsed.start != expected_offset {
            return Err(invalid_data(format!(
                "ranged ENA stream resumed at byte {}, expected byte {expected_offset}",
                parsed.start
            )));
        }
        if parsed.end < parsed.start || parsed.end.checked_add(1) != Some(parsed.total) {
            return Err(invalid_data(format!(
                "ranged ENA stream returned partial suffix Content-Range {content_range}"
            )));
        }
        if parsed.total != fastq.expected_bytes {
            return Err(invalid_data(format!(
                "ranged ENA stream total was {}, catalogue expected {} for {}",
                parsed.total, fastq.expected_bytes, fastq.accession
            )));
        }

        let content_length = response
            .headers()
            .get(CONTENT_LENGTH)
            .ok_or_else(|| {
                invalid_data("ranged ENA stream response did not include Content-Length")
            })
            .and_then(|value| parse_u64_header(value, "Content-Length"))?;
        let expected_length = parsed.end - parsed.start + 1;
        if content_length != expected_length {
            return Err(invalid_data(format!(
                "ranged ENA stream Content-Length was {content_length}, expected {expected_length} from Content-Range {content_range}"
            )));
        }

        Ok(Self::new(response, parsed.start..parsed.end + 1))
    }

    fn new(response: Response, promised: Range<u64>) -> Self {
        Self {
            response,
            promised,
            bytes_read: 0,
        }
    }

    fn promised_len(&self) -> u64 {
        self.promised.end - self.promised.start
    }

    fn bytes_remaining(&self) -> u64 {
        self.promised_len() - self.bytes_read
    }
}

impl io::Read for ActiveResponse {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if buf.is_empty() {
            return Ok(0);
        }

        let remaining = self.bytes_remaining();
        if remaining == 0 {
            return Ok(0);
        }

        let max_len = match usize::try_from(remaining) {
            Ok(remaining) => remaining.min(buf.len()),
            Err(_) => buf.len(),
        };
        let bytes_read = io::Read::read(&mut self.response, &mut buf[..max_len])?;
        if bytes_read == 0 {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                format!(
                    "ENA stream response ended early at byte {} of promised response span {}-{}",
                    self.promised.start + self.bytes_read,
                    self.promised.start,
                    self.promised.end,
                ),
            ));
        }

        let bytes_read_u64 = u64::try_from(bytes_read).map_err(|_| {
            invalid_data(format!(
                "ENA stream response returned a read larger than u64 can represent: {bytes_read} bytes"
            ))
        })?;
        if bytes_read_u64 > remaining {
            return Err(invalid_data(format!(
                "ENA stream response over-read promised response span {}-{}: read {} bytes with only {remaining} bytes remaining",
                self.promised.start, self.promised.end, bytes_read_u64
            )));
        }

        self.bytes_read += bytes_read_u64;
        Ok(bytes_read)
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct ContentRange {
    start: u64,
    end: u64,
    total: u64,
}

fn parse_content_range(value: &str) -> io::Result<ContentRange> {
    let Some(range) = value.strip_prefix("bytes ") else {
        return Err(invalid_data(format!(
            "Content-Range did not start with 'bytes ': {value}"
        )));
    };

    let Some((range_part, total_part)) = range.split_once('/') else {
        return Err(invalid_data(format!(
            "Content-Range did not include total size: {value}"
        )));
    };
    let Some((start, end)) = range_part.split_once('-') else {
        return Err(invalid_data(format!(
            "Content-Range did not include byte range: {value}"
        )));
    };

    Ok(ContentRange {
        start: start
            .parse()
            .map_err(|_| invalid_data(format!("Content-Range start was not numeric: {value}")))?,
        end: end
            .parse()
            .map_err(|_| invalid_data(format!("Content-Range end was not numeric: {value}")))?,
        total: total_part
            .parse()
            .map_err(|_| invalid_data(format!("Content-Range total was not numeric: {value}")))?,
    })
}

fn parse_u64_header(value: &HeaderValue, name: &str) -> io::Result<u64> {
    value
        .to_str()
        .map_err(|_| invalid_data(format!("{name} was not valid UTF-8")))?
        .parse()
        .map_err(|_| invalid_data(format!("{name} was not numeric")))
}

fn validate_identity_content_encoding(response: &Response) -> io::Result<()> {
    if let Some(value) = response.headers().get(CONTENT_ENCODING)
        && value
            .to_str()
            .map_or(true, |value| !value.eq_ignore_ascii_case("identity"))
    {
        return Err(invalid_data(format!(
            "ENA stream response used unexpected Content-Encoding {}; expected identity",
            value.to_str().unwrap_or("<non-UTF-8>")
        )));
    }

    Ok(())
}

fn invalid_data(message: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message.into())
}

fn response_text_or_error(response: Response) -> Result<String> {
    let status = response.status();
    ensure!(
        status == StatusCode::OK,
        "ENA filereport request failed with HTTP status {status}\n\
         help: this is a metadata lookup failure before FASTQ streaming starts"
    );
    Ok(response.text()?)
}

fn validate_run_accession(value: &str) -> Result<()> {
    let bytes = value.as_bytes();

    ensure!(
        bytes.len() >= 4,
        "ENA run accession must look like SRR12345, ERR12345, or DRR12345"
    );
    ensure!(
        matches!(&bytes[..3], b"SRR" | b"ERR" | b"DRR"),
        "ENA run accession must start with SRR, ERR, or DRR"
    );
    ensure!(
        bytes[3..].iter().all(u8::is_ascii_digit),
        "ENA run accession suffix must be numeric"
    );

    Ok(())
}

#[cfg(test)]
mod tests {
    use std::{
        io::{self, Read as _, Write as _},
        net::{SocketAddr, TcpListener, TcpStream},
        thread::{self, JoinHandle},
        time::Duration,
    };

    use color_eyre::eyre::{Result, bail, ensure, eyre};
    use reqwest::{StatusCode, blocking::Client};
    use url::Url;

    use crate::record::MateSide;

    use super::{Accession, EnaClient, EnaFastq, EnaInput};

    const TEST_FASTQ_BYTES: &[u8] = b"\x1f\x8babcdefghijklmnopqrstuvwxyz";
    const TEST_FASTQ_MD5: [u8; 16] = [
        0x45, 0x74, 0x41, 0xa3, 0x09, 0xc8, 0xfa, 0xc3, 0xf5, 0xc0, 0x59, 0xf2, 0x9b, 0x5f, 0x35,
        0x60,
    ];

    fn test_fastq(server: &TestServer) -> Result<EnaFastq> {
        Ok(EnaFastq {
            accession: Accession::new("SRR35939766")?,
            mate: None,
            url: Url::parse(&format!("http://{}{}", server.address, FASTQ_TARGET))?,
            expected_bytes: TEST_FASTQ_BYTES.len() as u64,
            expected_md5: TEST_FASTQ_MD5,
        })
    }

    fn test_client(max_retries: u32) -> Result<EnaClient> {
        let mut client = EnaClient::new()?;
        client.max_retries = max_retries;
        client.initial_backoff = Duration::from_millis(1);
        client.max_backoff = Duration::from_millis(2);
        Ok(client)
    }

    fn stream_error_after_partial(second_exchange: Exchange) -> Result<io::Error> {
        let server = TestServer::spawn(vec![
            Exchange::promised_suffix(TEST_FASTQ_BYTES, 0, 10),
            second_exchange,
        ])?;
        let fastq = test_fastq(&server)?;
        let mut stream = test_client(1)?.stream(fastq)?;
        let mut output = Vec::new();
        let error = stream
            .read_to_end(&mut output)
            .expect_err("invalid ranged response should fail");
        server.finish()?;
        Ok(error)
    }

    #[test]
    fn accession_accepts_valid_run_ids() -> Result<()> {
        assert_eq!(Accession::new("SRR35939766")?.as_str(), "SRR35939766");
        assert_eq!(Accession::new("ERR123456")?.as_str(), "ERR123456");
        assert_eq!(Accession::new("DRR987654")?.as_str(), "DRR987654");
        Ok(())
    }

    #[test]
    fn accession_rejects_invalid_run_ids() {
        assert!(Accession::new("PRJNA1247874").is_err());
        assert!(Accession::new("SRRabc").is_err());
        assert!(Accession::new("XYZ123").is_err());
    }

    #[test]
    fn filereport_resolves_valid_single_end_input() -> Result<()> {
        let accession = Accession::new("SRR35939766")?;
        let body = concat!(
            "run_accession\tlibrary_layout\tfastq_ftp\tfastq_bytes\tfastq_md5\n",
            "SRR35939766\tSINGLE\t",
            "ftp.sra.ebi.ac.uk/vol1/fastq/SRR35939766.fastq.gz\t",
            "26\td41d8cd98f00b204e9800998ecf8427e\n"
        );

        let input = EnaInput::from_filereport(&accession, body)?;
        let EnaInput::Single(fastq) = input else {
            bail!("expected single-end ENA input");
        };

        assert_eq!(fastq.accession, accession);
        assert_eq!(fastq.mate, None);
        assert_eq!(fastq.expected_bytes, 26);
        assert_eq!(
            fastq.expected_md5,
            [
                0xd4, 0x1d, 0x8c, 0xd9, 0x8f, 0x00, 0xb2, 0x04, 0xe9, 0x80, 0x09, 0x98, 0xec, 0xf8,
                0x42, 0x7e,
            ]
        );
        assert_eq!(
            fastq.url.as_str(),
            "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR35939766.fastq.gz"
        );
        Ok(())
    }

    #[test]
    fn filereport_resolves_valid_paired_end_input() -> Result<()> {
        let accession = Accession::new("SRR35939766")?;
        let body = concat!(
            "run_accession\tlibrary_layout\tfastq_ftp\tfastq_bytes\tfastq_md5\n",
            "SRR35939766\tPAIRED\t",
            "ftp.sra.ebi.ac.uk/vol1/fastq/SRR35939766_1.fastq.gz;",
            "ftp.sra.ebi.ac.uk/vol1/fastq/SRR35939766_2.fastq.gz\t",
            "17;19\t",
            "d41d8cd98f00b204e9800998ecf8427e;",
            "0cc175b9c0f1b6a831c399e269772661\n"
        );

        let input = EnaInput::from_filereport(&accession, body)?;
        let EnaInput::Paired { left, right } = input else {
            bail!("expected paired-end ENA input");
        };

        assert_eq!(left.mate, Some(MateSide::Left));
        assert_eq!(right.mate, Some(MateSide::Right));
        assert_eq!(left.expected_bytes, 17);
        assert_eq!(right.expected_bytes, 19);
        assert_eq!(
            left.url.as_str(),
            "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR35939766_1.fastq.gz"
        );
        assert_eq!(
            right.url.as_str(),
            "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR35939766_2.fastq.gz"
        );
        assert_eq!(
            left.expected_md5,
            [
                0xd4, 0x1d, 0x8c, 0xd9, 0x8f, 0x00, 0xb2, 0x04, 0xe9, 0x80, 0x09, 0x98, 0xec, 0xf8,
                0x42, 0x7e,
            ]
        );
        assert_eq!(
            right.expected_md5,
            [
                0x0c, 0xc1, 0x75, 0xb9, 0xc0, 0xf1, 0xb6, 0xa8, 0x31, 0xc3, 0x99, 0xe2, 0x69, 0x77,
                0x26, 0x61,
            ]
        );
        Ok(())
    }

    #[test]
    fn ena_stream_replays_gzip_prefix_and_verifies_md5() -> Result<()> {
        let server = TestServer::spawn(vec![
            Exchange::initial(TEST_FASTQ_BYTES)
                .with_header("Content-Type", "application/octet-stream"),
        ])?;
        let fastq = test_fastq(&server)?;
        let client = test_client(1)?;

        let mut stream = client.stream(fastq)?;
        let mut output = Vec::new();
        let read_result = stream.read_to_end(&mut output);
        server.finish()?;
        read_result?;

        assert_eq!(output, TEST_FASTQ_BYTES);
        Ok(())
    }

    #[test]
    fn ena_stream_retries_when_initial_body_ends_before_gzip_prefix() -> Result<()> {
        let server = TestServer::spawn(vec![
            Exchange::promised_suffix(TEST_FASTQ_BYTES, 0, 1),
            Exchange::initial(TEST_FASTQ_BYTES),
        ])?;
        let fastq = test_fastq(&server)?;

        let mut stream = test_client(1)?.stream(fastq)?;
        let mut output = Vec::new();
        let read_result = stream.read_to_end(&mut output);
        server.finish()?;
        read_result?;

        assert_eq!(output, TEST_FASTQ_BYTES);
        Ok(())
    }

    #[test]
    fn filereport_rejects_malformed_required_fields() -> Result<()> {
        let accession = Accession::new("SRR35939766")?;
        let cases = [
            (
                concat!(
                    "run_accession\tlibrary_layout\tfastq_ftp\tfastq_bytes\tfastq_md5\n",
                    "SRR35939766\tSINGLE\tftp.sra.ebi.ac.uk/a.fastq.gz\tnot-a-number\t",
                    "d41d8cd98f00b204e9800998ecf8427e\n"
                ),
                "not numeric",
            ),
            (
                concat!(
                    "run_accession\tlibrary_layout\tfastq_ftp\tfastq_bytes\tfastq_md5\n",
                    "SRR35939766\tSINGLE\tftp://example.test/a.fastq.gz\t12\t",
                    "d41d8cd98f00b204e9800998ecf8427e\n"
                ),
                "non-HTTPS URL scheme",
            ),
            (
                concat!(
                    "run_accession\tlibrary_layout\tfastq_ftp\tfastq_bytes\tfastq_md5\n",
                    "SRR35939766\tSINGLE\tftp.sra.ebi.ac.uk/a.fastq.gz\t12\tbad-md5\n"
                ),
                "32 hexadecimal",
            ),
            (
                concat!(
                    "run_accession\tlibrary_layout\tfastq_ftp\tfastq_md5\n",
                    "SRR35939766\tSINGLE\tftp.sra.ebi.ac.uk/a.fastq.gz\t",
                    "d41d8cd98f00b204e9800998ecf8427e\n"
                ),
                "did not include fastq_bytes",
            ),
            (
                concat!(
                    "run_accession\tlibrary_layout\tfastq_ftp\tfastq_bytes\tfastq_md5\n",
                    "SRR35939766\tSINGLE\tftp.sra.ebi.ac.uk/a.fastq.gz\t12\t\n"
                ),
                "fastq_md5 was empty",
            ),
            (
                concat!(
                    "run_accession\tlibrary_layout\tfastq_ftp\tfastq_bytes\tfastq_md5\n",
                    "SRR35939766\tSINGLE\tftp.sra.ebi.ac.uk/a.fastq.gz\t12\t",
                    "zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz\n"
                ),
                "was not hexadecimal",
            ),
        ];

        for (body, case) in cases {
            assert!(
                EnaInput::from_filereport(&accession, body).is_err(),
                "accepted malformed filereport case: {case}"
            );
        }
        Ok(())
    }

    #[test]
    fn filereport_rejects_inconsistent_run_metadata() -> Result<()> {
        let accession = Accession::new("SRR35939766")?;
        let cases = [
            (
                concat!(
                    "run_accession\tlibrary_layout\tfastq_ftp\tfastq_bytes\tfastq_md5\n",
                    "SRR35939766\tPAIRED\tftp.sra.ebi.ac.uk/a.fastq.gz;",
                    "ftp.sra.ebi.ac.uk/b.fastq.gz\t12\t",
                    "d41d8cd98f00b204e9800998ecf8427e;0cc175b9c0f1b6a831c399e269772661\n"
                ),
                "cardinalities differ",
            ),
            (
                concat!(
                    "run_accession\tlibrary_layout\tfastq_ftp\tfastq_bytes\tfastq_md5\n",
                    "SRR35939766\tSINGLE\tftp.sra.ebi.ac.uk/a.fastq.gz;",
                    "ftp.sra.ebi.ac.uk/b.fastq.gz\t12;13\t",
                    "d41d8cd98f00b204e9800998ecf8427e;0cc175b9c0f1b6a831c399e269772661\n"
                ),
                "SINGLE layout returned 2 FASTQ files",
            ),
            (
                concat!(
                    "run_accession\tlibrary_layout\tfastq_ftp\tfastq_bytes\tfastq_md5\n",
                    "SRR35939766\tPAIRED\tftp.sra.ebi.ac.uk/a.fastq.gz;",
                    "ftp.sra.ebi.ac.uk/a.fastq.gz\t12;12\t",
                    "d41d8cd98f00b204e9800998ecf8427e;d41d8cd98f00b204e9800998ecf8427e\n"
                ),
                "must be distinct",
            ),
            (
                concat!(
                    "run_accession\tlibrary_layout\tfastq_ftp\tfastq_bytes\tfastq_md5\n",
                    "ERR123456\tSINGLE\tftp.sra.ebi.ac.uk/a.fastq.gz\t12\t",
                    "d41d8cd98f00b204e9800998ecf8427e\n"
                ),
                "while resolving SRR35939766",
            ),
            (
                concat!(
                    "run_accession\tlibrary_layout\tfastq_ftp\tfastq_bytes\tfastq_md5\n",
                    "SRR35939766\tSINGLE\tftp.sra.ebi.ac.uk/a.fastq.gz\t12\t",
                    "d41d8cd98f00b204e9800998ecf8427e\n",
                    "SRR35939766\tSINGLE\tftp.sra.ebi.ac.uk/b.fastq.gz\t13\t",
                    "0cc175b9c0f1b6a831c399e269772661\n"
                ),
                "more than one run row",
            ),
        ];

        for (body, case) in cases {
            assert!(
                EnaInput::from_filereport(&accession, body).is_err(),
                "accepted malformed filereport case: {case}"
            );
        }
        Ok(())
    }

    #[test]
    fn ena_stream_rejects_catalogue_content_length_mismatch_during_open() -> Result<()> {
        let server = TestServer::spawn(vec![
            Exchange::initial(TEST_FASTQ_BYTES)
                .with_header("Content-Length", (TEST_FASTQ_BYTES.len() - 1).to_string()),
        ])?;
        let fastq = test_fastq(&server)?;

        let Err(error) = test_client(1)?.stream(fastq) else {
            bail!("catalogue size mismatch should fail during stream opening");
        };
        server.finish()?;

        assert_eq!(
            error.downcast_ref::<io::Error>().map(io::Error::kind),
            Some(io::ErrorKind::InvalidData)
        );
        assert!(error.to_string().contains("catalogue expected"));
        Ok(())
    }

    #[test]
    fn ena_stream_rejects_missing_content_length_during_open() -> Result<()> {
        let server = TestServer::spawn(vec![
            Exchange::initial(TEST_FASTQ_BYTES).without_header("Content-Length"),
        ])?;
        let fastq = test_fastq(&server)?;

        let Err(error) = test_client(3)?.stream(fastq) else {
            bail!("missing Content-Length should fail during stream opening");
        };
        server.finish()?;

        assert!(error.to_string().contains("Content-Length"));
        Ok(())
    }

    #[test]
    fn ena_stream_rejects_html_content_type_during_open() -> Result<()> {
        let server = TestServer::spawn(vec![
            Exchange::initial(TEST_FASTQ_BYTES)
                .with_header("Content-Type", "text/html; charset=utf-8"),
        ])?;
        let fastq = test_fastq(&server)?;

        let Err(error) = test_client(1)?.stream(fastq) else {
            bail!("HTML should fail during stream opening");
        };
        server.finish()?;

        assert!(error.to_string().contains("was HTML"));
        Ok(())
    }

    #[test]
    fn ena_stream_rejects_transformed_content_encoding_during_open() -> Result<()> {
        let server = TestServer::spawn(vec![
            Exchange::initial(TEST_FASTQ_BYTES).with_header("Content-Encoding", "gzip"),
        ])?;
        let fastq = test_fastq(&server)?;

        let Err(error) = test_client(3)?.stream(fastq) else {
            bail!("transformed Content-Encoding should fail during stream opening");
        };
        server.finish()?;

        assert!(error.to_string().contains("Content-Encoding"));
        Ok(())
    }

    #[test]
    fn ena_stream_rejects_wrong_gzip_magic_during_open() -> Result<()> {
        let payload = b"NOabcdefghijklmnopqrstuvwxyz";
        assert_eq!(payload.len(), TEST_FASTQ_BYTES.len());
        let server = TestServer::spawn(vec![Exchange::initial(payload)])?;
        let fastq = test_fastq(&server)?;

        let Err(error) = test_client(1)?.stream(fastq) else {
            bail!("wrong gzip magic should fail during stream opening");
        };
        server.finish()?;

        assert!(error.to_string().contains("gzip magic"));
        Ok(())
    }

    #[test]
    fn ena_stream_rejects_compressed_md5_mismatch() -> Result<()> {
        let server = TestServer::spawn(vec![Exchange::initial(TEST_FASTQ_BYTES)])?;
        let mut fastq = test_fastq(&server)?;
        fastq.expected_md5 = [0; 16];
        let mut stream = test_client(1)?.stream(fastq)?;

        let mut output = Vec::new();
        let error = stream
            .read_to_end(&mut output)
            .expect_err("compressed MD5 mismatch should fail");
        server.finish()?;

        assert_eq!(error.kind(), io::ErrorKind::InvalidData);
        assert!(error.to_string().contains("MD5 mismatch"));
        Ok(())
    }

    #[test]
    fn ena_stream_resumes_after_premature_eof() -> Result<()> {
        let server = TestServer::spawn(vec![
            Exchange::promised_suffix(TEST_FASTQ_BYTES, 0, 10),
            Exchange::ranged(TEST_FASTQ_BYTES, 10),
        ])?;
        let fastq = test_fastq(&server)?;
        let mut stream = test_client(2)?.stream(fastq)?;

        let mut output = Vec::new();
        let read_result = stream.read_to_end(&mut output);
        server.finish()?;
        read_result?;

        assert_eq!(output, TEST_FASTQ_BYTES);
        Ok(())
    }

    #[test]
    fn ena_stream_resumes_across_multiple_premature_eofs() -> Result<()> {
        let server = TestServer::spawn(vec![
            Exchange::promised_suffix(TEST_FASTQ_BYTES, 0, 8),
            Exchange::promised_suffix(TEST_FASTQ_BYTES, 8, 18),
            Exchange::ranged(TEST_FASTQ_BYTES, 18),
        ])?;
        let fastq = test_fastq(&server)?;
        let mut stream = test_client(3)?.stream(fastq)?;

        let mut output = Vec::new();
        let read_result = stream.read_to_end(&mut output);
        server.finish()?;
        read_result?;

        assert_eq!(output, TEST_FASTQ_BYTES);
        Ok(())
    }

    #[test]
    fn ena_stream_reports_exhausted_retry_budget() -> Result<()> {
        let error =
            stream_error_after_partial(Exchange::promised_suffix(TEST_FASTQ_BYTES, 10, 10))?;

        assert!(error.to_string().contains("retry budget exhausted"));
        Ok(())
    }

    #[test]
    fn ena_stream_rejects_invalid_ranged_responses() -> Result<()> {
        let cases = [
            (
                Exchange::ranged(TEST_FASTQ_BYTES, 10).with_status(StatusCode::OK),
                "failed with status",
            ),
            (
                Exchange::ranged(TEST_FASTQ_BYTES, 10).without_header("Content-Range"),
                "did not include Content-Range",
            ),
            (
                Exchange::ranged(TEST_FASTQ_BYTES, 10).with_header(
                    "Content-Range",
                    format!("bytes 11-27/{}", TEST_FASTQ_BYTES.len()),
                ),
                "resumed at byte",
            ),
            (
                Exchange::ranged(TEST_FASTQ_BYTES, 10).with_header(
                    "Content-Range",
                    format!("bytes 10-20/{}", TEST_FASTQ_BYTES.len()),
                ),
                "partial suffix",
            ),
            (
                Exchange::ranged(TEST_FASTQ_BYTES, 10).with_header("Content-Length", "17"),
                "Content-Length",
            ),
            (
                Exchange::ranged(TEST_FASTQ_BYTES, 10).with_header(
                    "Content-Range",
                    format!("bytes 10-28/{}", TEST_FASTQ_BYTES.len() + 1),
                ),
                "catalogue expected",
            ),
        ];

        for (exchange, expected_message) in cases {
            let error = stream_error_after_partial(exchange)?;
            assert!(
                error.to_string().contains(expected_message),
                "expected `{expected_message}` in `{error}`"
            );
        }
        Ok(())
    }

    #[test]
    fn ena_stream_rejects_directory_redirect_during_open() -> Result<()> {
        let redirect = Exchange {
            expected_target: FASTQ_TARGET.to_owned(),
            expected_range_offset: None,
            status: StatusCode::FOUND,
            headers: vec![
                ("Location".to_owned(), "/".to_owned()),
                ("Content-Length".to_owned(), "0".to_owned()),
            ],
            body: Vec::new(),
        };
        let mut directory = Exchange::initial(TEST_FASTQ_BYTES);
        directory.expected_target = "/".to_owned();
        let server = TestServer::spawn(vec![redirect, directory])?;
        let fastq = test_fastq(&server)?;

        let Err(error) = test_client(1)?.stream(fastq) else {
            bail!("directory redirect should fail during stream opening");
        };
        server.finish()?;

        assert!(error.to_string().contains("directory URL"));
        Ok(())
    }

    #[test]
    fn scripted_server_reports_request_range_mismatches() -> Result<()> {
        let server = TestServer::spawn(vec![Exchange::ranged(TEST_FASTQ_BYTES, 10)])?;
        let url = server.fastq_url()?;
        let client = Client::builder().build()?;

        let request_result = client
            .get(url.as_str())
            .header("Accept-Encoding", "identity")
            .send();
        assert!(request_result.is_err());

        let error = server
            .finish()
            .expect_err("request mismatch should reach the owning test");
        assert!(error.to_string().contains("expected `bytes=10-`"));
        assert!(error.to_string().contains("observed `<missing>`"));
        Ok(())
    }

    const FASTQ_TARGET: &str = "/reads.fastq.gz";
    const SHUTDOWN_TARGET: &str = "/__nuclease_test_shutdown";
    const MAX_REQUEST_HEADER_BYTES: usize = 16 * 1024;

    struct TestServer {
        address: SocketAddr,
        thread: Option<JoinHandle<Result<()>>>,
    }

    impl TestServer {
        fn spawn(exchanges: Vec<Exchange>) -> Result<Self> {
            let listener = TcpListener::bind("127.0.0.1:0")?;
            let address = listener.local_addr()?;
            let thread = thread::spawn(move || serve_exchanges(&listener, exchanges));

            Ok(Self {
                address,
                thread: Some(thread),
            })
        }

        fn fastq_url(&self) -> Result<Url> {
            Ok(Url::parse(&format!(
                "http://{}{}",
                self.address, FASTQ_TARGET
            ))?)
        }

        fn finish(mut self) -> Result<()> {
            let shutdown_result = request_shutdown(self.address);
            let server_result = self
                .thread
                .take()
                .ok_or_else(|| eyre!("test server thread was already joined"))?
                .join()
                .map_err(|_| eyre!("test server thread panicked"))?;

            server_result?;
            shutdown_result?;
            Ok(())
        }
    }

    impl Drop for TestServer {
        fn drop(&mut self) {
            if let Some(thread) = self.thread.take() {
                let _ = request_shutdown(self.address);
                let _ = thread.join();
            }
        }
    }

    struct Exchange {
        expected_target: String,
        expected_range_offset: Option<u64>,
        status: StatusCode,
        headers: Vec<(String, String)>,
        body: Vec<u8>,
    }

    impl Exchange {
        fn initial(payload: &[u8]) -> Self {
            Self::promised_suffix(payload, 0, payload.len())
        }

        fn ranged(payload: &[u8], offset: usize) -> Self {
            Self::promised_suffix(payload, offset, payload.len())
        }

        fn promised_suffix(payload: &[u8], offset: usize, delivered_end: usize) -> Self {
            let status = if offset == 0 {
                StatusCode::OK
            } else {
                StatusCode::PARTIAL_CONTENT
            };
            let mut headers = vec![
                (
                    "Content-Length".to_owned(),
                    (payload.len() - offset).to_string(),
                ),
                ("Accept-Ranges".to_owned(), "bytes".to_owned()),
                ("Content-Type".to_owned(), "application/x-gzip".to_owned()),
            ];
            if offset > 0 {
                headers.push((
                    "Content-Range".to_owned(),
                    format!("bytes {}-{}/{}", offset, payload.len() - 1, payload.len()),
                ));
            }

            Self {
                expected_target: FASTQ_TARGET.to_owned(),
                expected_range_offset: (offset > 0).then_some(offset as u64),
                status,
                headers,
                body: payload[offset..delivered_end].to_vec(),
            }
        }

        fn with_status(mut self, status: StatusCode) -> Self {
            self.status = status;
            self
        }

        fn with_header(mut self, name: &str, value: impl Into<String>) -> Self {
            let value = value.into();
            if let Some((_, existing)) = self
                .headers
                .iter_mut()
                .find(|(header, _)| header.eq_ignore_ascii_case(name))
            {
                *existing = value;
            } else {
                self.headers.push((name.to_owned(), value));
            }
            self
        }

        fn without_header(mut self, name: &str) -> Self {
            self.headers
                .retain(|(header, _)| !header.eq_ignore_ascii_case(name));
            self
        }

        fn serve(self, stream: &mut TcpStream, request: &[u8]) -> Result<()> {
            self.verify_request(request)?;
            self.write_response(stream)
        }

        fn verify_request(&self, request: &[u8]) -> Result<()> {
            let request = std::str::from_utf8(request)?;
            let mut lines = request.split("\r\n");
            let request_line = lines
                .next()
                .ok_or_else(|| eyre!("HTTP request did not include a request line"))?;
            ensure!(
                request_line == format!("GET {} HTTP/1.1", self.expected_target),
                "unexpected HTTP request line: expected `GET {} HTTP/1.1`, observed `{request_line}`",
                self.expected_target
            );

            let mut accept_encoding = None;
            let mut range = None;
            for line in lines.take_while(|line| !line.is_empty()) {
                let Some((name, value)) = line.split_once(':') else {
                    bail!("malformed HTTP request header: {line}");
                };
                let value = value.trim();
                if name.eq_ignore_ascii_case("Accept-Encoding") {
                    ensure!(
                        accept_encoding.replace(value).is_none(),
                        "HTTP request repeated Accept-Encoding"
                    );
                } else if name.eq_ignore_ascii_case("Range") {
                    ensure!(
                        range.replace(value).is_none(),
                        "HTTP request repeated Range"
                    );
                }
            }

            ensure!(
                accept_encoding.is_some_and(|value| value.eq_ignore_ascii_case("identity")),
                "HTTP request did not use Accept-Encoding: identity"
            );
            match self.expected_range_offset {
                Some(offset) => ensure!(
                    range == Some(format!("bytes={offset}-").as_str()),
                    "unexpected Range header: expected `bytes={offset}-`, observed `{}`",
                    range.unwrap_or("<missing>")
                ),
                None => ensure!(
                    range.is_none(),
                    "unexpected Range header on initial request: {}",
                    range.unwrap_or("<missing>")
                ),
            }

            Ok(())
        }

        fn write_response(&self, stream: &mut TcpStream) -> Result<()> {
            let reason = self.status.canonical_reason().unwrap_or("Unknown");
            let mut response = format!("HTTP/1.1 {} {reason}\r\n", self.status.as_u16());
            for (name, value) in &self.headers {
                response.push_str(name);
                response.push_str(": ");
                response.push_str(value);
                response.push_str("\r\n");
            }
            response.push_str("Connection: close\r\n\r\n");

            stream.write_all(response.as_bytes())?;
            stream.write_all(&self.body)?;
            stream.flush()?;
            Ok(())
        }
    }

    fn serve_exchanges(listener: &TcpListener, exchanges: Vec<Exchange>) -> Result<()> {
        for (index, exchange) in exchanges.into_iter().enumerate() {
            let (mut stream, _) = listener.accept()?;
            let request = read_http_request(&mut stream)?;
            ensure!(
                !is_shutdown_request(&request),
                "test server was shut down before scripted exchange {} arrived",
                index + 1
            );
            exchange.serve(&mut stream, &request)?;
        }

        let (mut stream, _) = listener.accept()?;
        let request = read_http_request(&mut stream)?;
        ensure!(
            is_shutdown_request(&request),
            "test server received an unexpected request after its script was exhausted: {}",
            String::from_utf8_lossy(&request)
                .split("\r\n")
                .next()
                .unwrap_or("<missing request line>")
        );
        Ok(())
    }

    fn read_http_request(stream: &mut TcpStream) -> Result<Vec<u8>> {
        stream.set_read_timeout(Some(Duration::from_secs(5)))?;
        let mut request = Vec::with_capacity(1024);
        let mut chunk = [0_u8; 512];

        loop {
            let bytes_read = stream.read(&mut chunk)?;
            ensure!(
                bytes_read > 0,
                "HTTP client closed before completing request headers"
            );
            request.extend_from_slice(&chunk[..bytes_read]);
            ensure!(
                request.len() <= MAX_REQUEST_HEADER_BYTES,
                "HTTP request headers exceeded {MAX_REQUEST_HEADER_BYTES} bytes"
            );
            if request.windows(4).any(|window| window == b"\r\n\r\n") {
                return Ok(request);
            }
        }
    }

    fn is_shutdown_request(request: &[u8]) -> bool {
        request.starts_with(format!("GET {SHUTDOWN_TARGET} HTTP/1.1\r\n").as_bytes())
    }

    fn request_shutdown(address: SocketAddr) -> Result<()> {
        let mut stream = TcpStream::connect(address)?;
        write!(
            stream,
            "GET {SHUTDOWN_TARGET} HTTP/1.1\r\nHost: {address}\r\nConnection: close\r\n\r\n"
        )?;
        stream.flush()?;
        Ok(())
    }
}
