//! ENA accession modeling, file-report resolution, and retrying HTTP byte-stream readers.

use std::{fmt, io, ops::Range, str::FromStr, thread, time::Duration};

use color_eyre::eyre::{Result, WrapErr, bail, ensure, eyre};
use reqwest::{
    StatusCode,
    blocking::{Client, Response},
    header::{
        ACCEPT_ENCODING, CONTENT_ENCODING, CONTENT_LENGTH, CONTENT_RANGE, HeaderValue, RANGE,
    },
};
use url::Url;

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
/// Validated HTTPS URL for a gzipped FASTQ resource.
pub struct FastqUrl(Url);

impl FastqUrl {
    /// Construct a validated HTTPS FASTQ URL ending in `.fastq.gz`.
    ///
    /// # Errors
    ///
    /// Returns an error when the URL is empty, invalid, not HTTPS, or does not point at a
    /// gzipped FASTQ path.
    pub fn new(value: impl Into<String>) -> Result<Self> {
        let value = value.into();
        if value.is_empty() {
            bail!("FASTQ URL must not be empty");
        }

        let url = Url::parse(&value)?;
        ensure!(url.scheme() == "https", "FASTQ URL must use https");
        ensure!(
            url.path().ends_with(".fastq.gz"),
            "FASTQ URL path must end with .fastq.gz"
        );

        Ok(Self(url))
    }

    /// Borrow the normalized URL as a string slice.
    pub fn as_str(&self) -> &str {
        self.0.as_str()
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
/// Distinct paired FASTQ URLs representing mate 1 and mate 2.
pub struct PairedFastqUrls {
    r1: FastqUrl,
    r2: FastqUrl,
}

impl PairedFastqUrls {
    /// Construct paired FASTQ URLs while enforcing that the two mates differ.
    ///
    /// # Errors
    ///
    /// Returns an error when the mate URLs are identical.
    pub fn new(r1: FastqUrl, r2: FastqUrl) -> Result<Self> {
        ensure!(r1 != r2, "paired FASTQ URLs must be distinct");
        Ok(Self { r1, r2 })
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
/// FASTQ URLs returned by ENA, grouped by inferred layout.
pub enum FastqUrlsByLayout {
    /// Single-end FASTQ layout.
    Single(FastqUrl),
    /// Paired-end FASTQ layout.
    Paired(PairedFastqUrls),
}

#[derive(Debug)]
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

    /// Resolve an ENA accession into validated FASTQ URLs grouped by layout.
    ///
    /// # Errors
    ///
    /// Returns an error when the filereport request fails, returns malformed data, or yields an
    /// unsupported FASTQ URL layout.
    pub fn lookup_fastq_urls(&self, accession: &Accession) -> Result<FastqUrlsByLayout> {
        let response = self
            .http
            .get(ENA_FILEREPORT_BASE_URL)
            .query(&[
                ("accession", accession.as_str()),
                ("result", "read_run"),
                ("fields", "run_accession,fastq_ftp,library_layout"),
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
        let fastq_ftp_field = extract_fastq_ftp_field(accession, &body)?;
        parse_fastq_urls_by_layout(&fastq_ftp_field).wrap_err_with(|| {
            format!(
                "ENA filereport fastq_ftp field could not be interpreted for accession {accession}\n\
                 fastq_ftp: {fastq_ftp_field}\n\
                 help: nuclease currently supports one FASTQ URL for single-end runs or two FASTQ URLs for paired-end runs"
            )
        })
    }

    /// Open a retrying byte stream for one FASTQ URL.
    pub fn open_retrying_stream(&self, url: FastqUrl) -> RetryingHttpRead {
        RetryingHttpRead::new(
            self.http.clone(),
            url,
            self.max_retries,
            self.initial_backoff,
            self.max_backoff,
        )
    }

    /// Open retrying byte streams for paired FASTQ URLs.
    pub fn open_retrying_paired_streams(
        &self,
        urls: PairedFastqUrls,
    ) -> (RetryingHttpRead, RetryingHttpRead) {
        let r1 = self.open_retrying_stream(urls.r1);
        let r2 = self.open_retrying_stream(urls.r2);
        (r1, r2)
    }
}

/// Blocking `Read` implementation that resumes ENA HTTP downloads with ranged reconnects.
pub struct RetryingHttpRead {
    http: Client,
    url: FastqUrl,
    body: Option<ActiveResponse>,
    next_byte_offset: u64,
    expected_total_bytes: Option<u64>,
    max_retries: u32,
    retries_remaining: u32,
    initial_backoff: Duration,
    current_backoff: Duration,
    max_backoff: Duration,
    saw_eof: bool,
}

impl RetryingHttpRead {
    /// Construct a retrying ENA reader with an explicit retry budget and backoff policy.
    pub fn new(
        http: Client,
        url: FastqUrl,
        max_retries: u32,
        initial_backoff: Duration,
        max_backoff: Duration,
    ) -> Self {
        Self {
            http,
            url,
            body: None,
            next_byte_offset: 0,
            expected_total_bytes: None,
            max_retries,
            retries_remaining: max_retries,
            initial_backoff,
            current_backoff: initial_backoff,
            max_backoff,
            saw_eof: false,
        }
    }

    fn ensure_connected(&mut self) -> io::Result<()> {
        if self.body.is_some() {
            return Ok(());
        }

        let response = if self.next_byte_offset == 0 {
            self.open_response(None)?
        } else {
            self.open_response(Some(self.next_byte_offset))?
        };
        self.body = Some(response);
        Ok(())
    }

    fn open_response(&mut self, offset: Option<u64>) -> io::Result<ActiveResponse> {
        let mut request = self
            .http
            .get(self.url.as_str())
            .header(ACCEPT_ENCODING, "identity");
        if let Some(offset) = offset {
            request = request.header(RANGE, format!("bytes={offset}-"));
        }

        let response = request.send().map_err(io_error_from_reqwest)?;
        let span = match offset {
            None => self.validate_initial_response(&response)?,
            Some(expected_offset) => self.validate_ranged_response(&response, expected_offset)?,
        };
        Ok(ActiveResponse::new(response, span))
    }

    fn validate_initial_response(&mut self, response: &Response) -> io::Result<Range<u64>> {
        if response.status() != StatusCode::OK {
            return Err(io::Error::other(format!(
                "initial ENA stream request failed with status {}",
                response.status()
            )));
        }
        validate_identity_content_encoding(response)?;

        let Some(content_length) = response.headers().get(CONTENT_LENGTH) else {
            return Err(invalid_data(
                "initial ENA stream response did not include Content-Length",
            ));
        };
        let total = parse_u64_header(content_length, "Content-Length")?;
        self.remember_expected_total(total)?;

        Ok(0..total)
    }

    fn validate_ranged_response(
        &mut self,
        response: &Response,
        expected_offset: u64,
    ) -> io::Result<Range<u64>> {
        if response.status() != StatusCode::PARTIAL_CONTENT {
            let message = format!(
                "ranged ENA stream request failed with status {}",
                response.status()
            );
            return if response.status().is_success() {
                Err(invalid_data(message))
            } else {
                Err(io::Error::other(message))
            };
        }
        validate_identity_content_encoding(response)?;

        let Some(content_range) = response.headers().get(CONTENT_RANGE) else {
            return Err(invalid_data(
                "ranged ENA stream response did not include Content-Range",
            ));
        };
        let content_range = content_range
            .to_str()
            .map_err(|_| invalid_data("Content-Range header was not valid UTF-8"))?;
        let parsed = parse_content_range(content_range)?;
        if parsed.start != expected_offset {
            return Err(invalid_data(format!(
                "ranged ENA stream resumed at byte {}, expected byte {}",
                parsed.start, expected_offset
            )));
        }
        if parsed.end < parsed.start {
            return Err(invalid_data(format!(
                "ranged ENA stream returned invalid Content-Range: {content_range}"
            )));
        }
        if parsed.end.checked_add(1) != Some(parsed.total) {
            return Err(invalid_data(format!(
                "ranged ENA stream returned partial suffix Content-Range {content_range}; expected byte range to end at {}",
                parsed.total.saturating_sub(1)
            )));
        }

        let Some(content_length) = response.headers().get(CONTENT_LENGTH) else {
            return Err(invalid_data(
                "ranged ENA stream response did not include Content-Length",
            ));
        };
        let observed_len = parse_u64_header(content_length, "Content-Length")?;
        let expected_len = parsed.end - parsed.start + 1;
        if observed_len != expected_len {
            return Err(invalid_data(format!(
                "ranged ENA stream Content-Length was {observed_len}, expected {expected_len} from Content-Range {content_range}"
            )));
        }

        self.remember_expected_total(parsed.total)?;

        Ok(parsed.start..parsed.end + 1)
    }

    fn remember_expected_total(&mut self, total: u64) -> io::Result<()> {
        match self.expected_total_bytes {
            Some(expected) if expected != total => Err(invalid_data(format!(
                "ENA stream total size changed across retries: first saw {expected} bytes, then saw {total} bytes"
            ))),
            Some(_) => Ok(()),
            None => {
                self.expected_total_bytes = Some(total);
                Ok(())
            }
        }
    }

    fn has_reached_expected_eof(&self) -> bool {
        match self.expected_total_bytes {
            Some(total) => self.next_byte_offset == total,
            None => true,
        }
    }

    fn retry_or_return(&mut self, error: io::Error) -> io::Result<()> {
        if !Self::should_retry_io_error(&error) {
            return Err(error);
        }

        self.consume_retry_budget()?;
        self.drop_connection();
        thread::sleep(self.current_backoff);
        self.advance_backoff();
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

    fn consume_retry_budget(&mut self) -> io::Result<()> {
        if self.retries_remaining == 0 {
            return Err(io::Error::other("ENA stream retry budget exhausted"));
        }

        self.retries_remaining -= 1;
        Ok(())
    }

    fn reset_backoff_after_success(&mut self) {
        self.retries_remaining = self.max_retries;
        self.current_backoff = self.initial_backoff;
    }

    fn drop_connection(&mut self) {
        self.body = None;
    }

    fn advance_backoff(&mut self) {
        self.current_backoff = self.current_backoff.saturating_mul(2).min(self.max_backoff);
    }
}

#[derive(Debug)]
struct ActiveResponse {
    response: Response,
    promised: Range<u64>,
    bytes_read: u64,
}

impl ActiveResponse {
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

fn extract_fastq_ftp_field(accession: &Accession, body: &str) -> Result<String> {
    let mut lines = body.lines();
    let Some(header) = lines.next() else {
        bail!(
            "ENA filereport response was empty for accession {accession}\n\
             help: confirm the accession is public and points to a run, not a study or sample"
        );
    };
    let Some(row) = lines.next() else {
        bail!(
            "ENA filereport did not return a data row for accession {accession}\n\
             help: ENA may not have public FASTQ files for this run, or the accession may not be run-level"
        );
    };

    let header_fields = header.split('\t').collect::<Vec<_>>();
    let row_fields = row.split('\t').collect::<Vec<_>>();
    ensure!(
        header_fields.len() == row_fields.len(),
        "ENA filereport row shape did not match header for accession {accession}\n\
         header_fields: {} row_fields: {}\n\
         help: this usually indicates an unexpected ENA API response shape",
        header_fields.len(),
        row_fields.len(),
    );

    let mut run_accession = None;
    let mut fastq_ftp = None;
    let mut library_layout = None;

    for (name, value) in header_fields.iter().zip(row_fields.iter()) {
        match *name {
            "run_accession" => run_accession = Some(*value),
            "fastq_ftp" => fastq_ftp = Some(*value),
            "library_layout" => library_layout = Some(*value),
            _ => {}
        }
    }

    ensure!(
        run_accession == Some(accession.as_str()),
        "ENA filereport returned an unexpected run accession while resolving {accession}\n\
         returned: {}\n\
         help: retry the request; if it repeats, ENA may have returned a stale or malformed row",
        run_accession.unwrap_or("<missing>"),
    );
    ensure!(
        library_layout.is_some(),
        "ENA filereport response did not include library_layout"
    );

    fastq_ftp
        .map(ToOwned::to_owned)
        .ok_or_else(|| {
            eyre!(
                "ENA filereport response did not include fastq_ftp for accession {accession}\n\
                 help: nuclease needs ENA-hosted FASTQ URLs; this run may only expose submitted BAM/CRAM or other file types"
            )
        })
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

/// Parse ENA `fastq_ftp` fields into validated FASTQ URLs grouped by layout.
///
/// # Errors
///
/// Returns an error when the field is empty, contains invalid FASTQ URLs, or yields an
/// unsupported number of FASTQ URLs.
pub fn parse_fastq_urls_by_layout(fastq_ftp_field: &str) -> Result<FastqUrlsByLayout> {
    let urls = fastq_ftp_field
        .split(';')
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .map(|value| {
            if value.starts_with("https://") {
                FastqUrl::new(value)
            } else {
                FastqUrl::new(format!("https://{value}"))
            }
        })
        .collect::<Result<Vec<_>>>()?;

    match urls.as_slice() {
        [single] => Ok(FastqUrlsByLayout::Single(single.clone())),
        [r1, r2] => Ok(FastqUrlsByLayout::Paired(PairedFastqUrls::new(
            r1.clone(),
            r2.clone(),
        )?)),
        [] => bail!(
            "ENA fastq_ftp field did not contain any FASTQ URLs\n\
             help: this run may not have generated FASTQ files available in ENA"
        ),
        _ => bail!(
            "ENA fastq_ftp field contained an unsupported number of FASTQ URLs\n\
             observed_url_count: {}\n\
             help: nuclease currently supports single-end runs with 1 URL and paired-end runs with 2 URLs; choose a specific run with a simpler layout or download/stage the desired FASTQs locally",
            urls.len(),
        ),
    }
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

impl io::Read for RetryingHttpRead {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if buf.is_empty() {
            return Ok(0);
        }

        if self.saw_eof {
            return Ok(0);
        }

        loop {
            if let Err(error) = self.ensure_connected() {
                self.retry_or_return(error)?;
                continue;
            }

            let Some(body) = self.body.as_mut() else {
                return Err(io::Error::other("missing ENA response body"));
            };

            match body.read(buf) {
                Ok(0) => {
                    if self.has_reached_expected_eof() {
                        self.saw_eof = true;
                        return Ok(0);
                    }

                    self.retry_or_return(io::Error::new(
                        io::ErrorKind::UnexpectedEof,
                        format!(
                            "ENA stream ended early at byte {} of expected {}",
                            self.next_byte_offset,
                            self.expected_total_bytes
                                .map_or_else(|| "unknown".to_owned(), |total| total.to_string())
                        ),
                    ))?;
                }
                Ok(n) => {
                    self.next_byte_offset += n as u64;
                    self.reset_backoff_after_success();
                    return Ok(n);
                }
                Err(error) => {
                    self.retry_or_return(error)?;
                }
            }
        }
    }
}

fn io_error_from_reqwest(error: reqwest::Error) -> io::Error {
    io::Error::other(error)
}

#[cfg(test)]
mod tests {
    use std::{
        io::{Read as _, Write as _},
        net::{SocketAddr, TcpListener, TcpStream},
        thread::{self, JoinHandle},
        time::Duration,
    };

    use color_eyre::eyre::{Result, bail, ensure, eyre};
    use reqwest::{StatusCode, blocking::Client};
    use url::Url;

    use super::{
        Accession, EnaClient, FastqUrl, FastqUrlsByLayout, PairedFastqUrls, RetryingHttpRead,
        extract_fastq_ftp_field, parse_content_range, parse_fastq_urls_by_layout,
    };

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
    fn paired_fastq_urls_keep_mates_distinct() -> Result<()> {
        let paired = PairedFastqUrls::new(
            FastqUrl::new("https://example.test/read_1.fastq.gz")?,
            FastqUrl::new("https://example.test/read_2.fastq.gz")?,
        )?;

        assert_eq!(paired.r1.as_str(), "https://example.test/read_1.fastq.gz");
        assert_eq!(paired.r2.as_str(), "https://example.test/read_2.fastq.gz");
        Ok(())
    }

    #[test]
    fn fastq_url_rejects_non_https_and_non_fastq_gz_urls() {
        assert!(FastqUrl::new("http://example.test/read_1.fastq.gz").is_err());
        assert!(FastqUrl::new("https://example.test/read_1.fastq").is_err());
        assert!(FastqUrl::new("https://example.test/").is_err());
    }

    #[test]
    fn paired_fastq_urls_reject_identical_mates() -> Result<()> {
        let mate = FastqUrl::new("https://example.test/read_1.fastq.gz")?;
        assert!(PairedFastqUrls::new(mate.clone(), mate).is_err());
        Ok(())
    }

    #[test]
    fn fastq_urls_by_layout_can_represent_paired_urls() -> Result<()> {
        let urls = FastqUrlsByLayout::Paired(PairedFastqUrls::new(
            FastqUrl::new("https://example.test/read_1.fastq.gz")?,
            FastqUrl::new("https://example.test/read_2.fastq.gz")?,
        )?);

        assert!(matches!(urls, FastqUrlsByLayout::Paired(_)));
        Ok(())
    }

    #[test]
    fn parse_fastq_urls_by_layout_accepts_single_url() -> Result<()> {
        let urls = parse_fastq_urls_by_layout("ftp.sra.ebi.ac.uk/vol1/fastq/SRR1.fastq.gz")?;

        assert!(matches!(urls, FastqUrlsByLayout::Single(_)));
        Ok(())
    }

    #[test]
    fn parse_fastq_urls_by_layout_accepts_paired_urls() -> Result<()> {
        let urls = parse_fastq_urls_by_layout(
            "ftp.sra.ebi.ac.uk/vol1/fastq/SRR1_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/SRR1_2.fastq.gz",
        )?;

        assert!(matches!(urls, FastqUrlsByLayout::Paired(_)));
        Ok(())
    }

    #[test]
    fn parse_fastq_urls_by_layout_rejects_more_than_two_urls() {
        let result = parse_fastq_urls_by_layout(
            "ftp.sra.ebi.ac.uk/a.fastq.gz;ftp.sra.ebi.ac.uk/b.fastq.gz;ftp.sra.ebi.ac.uk/c.fastq.gz",
        );

        assert!(result.is_err());
    }

    #[test]
    fn extract_fastq_ftp_field_reads_ena_tsv_row() -> Result<()> {
        let accession = Accession::new("SRR35939766")?;
        let body = concat!(
            "run_accession\tfastq_ftp\tlibrary_layout\n",
            "SRR35939766\t",
            "ftp.sra.ebi.ac.uk/vol1/fastq/SRR359/066/SRR35939766/SRR35939766_1.fastq.gz;",
            "ftp.sra.ebi.ac.uk/vol1/fastq/SRR359/066/SRR35939766/SRR35939766_2.fastq.gz\t",
            "PAIRED\n"
        );

        let fastq_ftp = extract_fastq_ftp_field(&accession, body)?;
        assert!(fastq_ftp.contains("SRR35939766_1.fastq.gz"));
        assert!(fastq_ftp.contains("SRR35939766_2.fastq.gz"));
        Ok(())
    }

    #[test]
    fn ena_client_is_constructible() -> Result<()> {
        let _client = EnaClient::new()?;
        Ok(())
    }

    #[test]
    fn retrying_http_read_streams_initial_response_without_range() -> Result<()> {
        let payload = b"abcdefghijklmnopqrstuvwxyz".to_vec();
        let server = TestServer::spawn(vec![Exchange::initial(&payload)])?;
        let url = server.fastq_url()?;
        let client = Client::builder().build()?;
        let mut reader = RetryingHttpRead::new(
            client,
            url,
            1,
            Duration::from_millis(1),
            Duration::from_millis(2),
        );

        let mut output = Vec::new();
        let read_result = reader.read_to_end(&mut output);
        server.finish()?;
        read_result?;

        assert_eq!(output, payload);
        Ok(())
    }

    #[test]
    fn retrying_http_read_requests_identity_encoding() -> Result<()> {
        let payload = b"abcdefghijklmnopqrstuvwxyz".to_vec();
        let server = TestServer::spawn(vec![Exchange::initial(&payload)])?;
        let url = server.fastq_url()?;
        let client = Client::builder().build()?;
        let mut reader = RetryingHttpRead::new(
            client,
            url,
            1,
            Duration::from_millis(1),
            Duration::from_millis(2),
        );

        let mut output = Vec::new();
        let read_result = reader.read_to_end(&mut output);
        server.finish()?;
        read_result?;

        assert_eq!(output, payload);
        Ok(())
    }

    #[test]
    fn retrying_http_read_resumes_after_premature_eof() -> Result<()> {
        let payload = b"abcdefghijklmnopqrstuvwxyz".to_vec();
        let server = TestServer::spawn(vec![
            Exchange::promised_suffix(&payload, 0, 10),
            Exchange::ranged(&payload, 10),
        ])?;
        let url = server.fastq_url()?;
        let client = Client::builder().build()?;
        let mut reader = RetryingHttpRead::new(
            client,
            url,
            2,
            Duration::from_millis(1),
            Duration::from_millis(2),
        );

        let mut output = Vec::new();
        let read_result = reader.read_to_end(&mut output);
        server.finish()?;
        read_result?;

        assert_eq!(output, payload);
        Ok(())
    }

    #[test]
    fn retrying_http_read_resumes_across_multiple_premature_eofs() -> Result<()> {
        let payload = b"abcdefghijklmnopqrstuvwxyz".to_vec();
        let server = TestServer::spawn(vec![
            Exchange::promised_suffix(&payload, 0, 8),
            Exchange::promised_suffix(&payload, 8, 18),
            Exchange::ranged(&payload, 18),
        ])?;
        let url = server.fastq_url()?;
        let client = Client::builder().build()?;
        let mut reader = RetryingHttpRead::new(
            client,
            url,
            3,
            Duration::from_millis(1),
            Duration::from_millis(2),
        );

        let mut output = Vec::new();
        let read_result = reader.read_to_end(&mut output);
        server.finish()?;
        read_result?;

        assert_eq!(output, payload);
        Ok(())
    }

    #[test]
    fn retrying_http_read_errors_when_premature_eof_exhausts_retry_budget() -> Result<()> {
        let payload = b"abcdefghijklmnopqrstuvwxyz".to_vec();
        let server = TestServer::spawn(vec![
            Exchange::promised_suffix(&payload, 0, 10),
            Exchange::promised_suffix(&payload, 10, 10),
        ])?;
        let url = server.fastq_url()?;
        let client = Client::builder().build()?;
        let mut reader = RetryingHttpRead::new(
            client,
            url,
            1,
            Duration::from_millis(1),
            Duration::from_millis(2),
        );

        let mut output = Vec::new();
        let error = reader
            .read_to_end(&mut output)
            .expect_err("exhausted retry budget should fail");
        server.finish()?;

        assert!(error.to_string().contains("retry budget exhausted"));
        Ok(())
    }

    #[test]
    fn retrying_http_read_resumes_from_nonzero_offset_with_range() -> Result<()> {
        let payload = b"abcdefghijklmnopqrstuvwxyz".to_vec();
        let expected_offset = 10_usize;
        let server = TestServer::spawn(vec![Exchange::ranged(&payload, expected_offset)])?;
        let url = server.fastq_url()?;
        let client = Client::builder().build()?;
        let mut reader = RetryingHttpRead::new(
            client,
            url,
            1,
            Duration::from_millis(1),
            Duration::from_millis(2),
        );
        reader.next_byte_offset = expected_offset as u64;

        let mut output = Vec::new();
        let read_result = reader.read_to_end(&mut output);
        server.finish()?;
        read_result?;

        assert_eq!(output, payload[expected_offset..]);
        Ok(())
    }

    #[test]
    fn retrying_http_read_rejects_bad_ranged_status() -> Result<()> {
        let payload = b"abcdefghijklmnopqrstuvwxyz".to_vec();
        let expected_offset = 10_usize;
        let server = TestServer::spawn(vec![
            Exchange::ranged(&payload, expected_offset).with_status(StatusCode::OK),
        ])?;
        let url = server.fastq_url()?;
        let client = Client::builder().build()?;
        let mut reader = RetryingHttpRead::new(
            client,
            url,
            1,
            Duration::from_millis(1),
            Duration::from_millis(2),
        );

        let error = reader
            .open_response(Some(expected_offset as u64))
            .expect_err("bad ranged status should fail");
        server.finish()?;
        assert!(
            error
                .to_string()
                .contains("ranged ENA stream request failed with status")
        );
        Ok(())
    }

    #[test]
    fn retrying_http_read_rejects_missing_content_range() -> Result<()> {
        let payload = b"abcdefghijklmnopqrstuvwxyz".to_vec();
        let expected_offset = 10_usize;
        let server = TestServer::spawn(vec![
            Exchange::ranged(&payload, expected_offset).without_header("Content-Range"),
        ])?;
        let url = server.fastq_url()?;
        let client = Client::builder().build()?;
        let mut reader = RetryingHttpRead::new(
            client,
            url,
            1,
            Duration::from_millis(1),
            Duration::from_millis(2),
        );

        let error = reader
            .open_response(Some(expected_offset as u64))
            .expect_err("missing content-range should fail");
        server.finish()?;
        assert!(
            error
                .to_string()
                .contains("ranged ENA stream response did not include Content-Range")
        );
        Ok(())
    }

    #[test]
    fn content_range_parser_requires_start_end_and_total() -> Result<()> {
        let parsed = parse_content_range("bytes 10-25/26")?;

        assert_eq!(parsed.start, 10);
        assert_eq!(parsed.end, 25);
        assert_eq!(parsed.total, 26);
        assert!(parse_content_range("items 10-25/26").is_err());
        assert!(parse_content_range("bytes 10-25/*").is_err());
        Ok(())
    }

    #[test]
    fn retrying_http_read_rejects_wrong_ranged_start() -> Result<()> {
        let payload = b"abcdefghijklmnopqrstuvwxyz".to_vec();
        let expected_offset = 10_usize;
        let server = TestServer::spawn(vec![
            Exchange::ranged(&payload, expected_offset).with_header(
                "Content-Range",
                format!(
                    "bytes {}-{}/{}",
                    expected_offset + 1,
                    payload.len() - 1,
                    payload.len()
                ),
            ),
        ])?;
        let url = server.fastq_url()?;
        let client = Client::builder().build()?;
        let mut reader = RetryingHttpRead::new(
            client,
            url,
            1,
            Duration::from_millis(1),
            Duration::from_millis(2),
        );

        let error = reader
            .open_response(Some(expected_offset as u64))
            .expect_err("wrong ranged start should fail");
        server.finish()?;
        assert!(error.to_string().contains("resumed at byte"));
        Ok(())
    }

    #[test]
    fn retrying_http_read_rejects_ranged_partial_suffix() -> Result<()> {
        let payload = b"abcdefghijklmnopqrstuvwxyz".to_vec();
        let expected_offset = 10_usize;
        let partial_end = expected_offset + 4;
        let server = TestServer::spawn(vec![
            Exchange::promised_suffix(&payload, expected_offset, partial_end)
                .with_header(
                    "Content-Length",
                    (partial_end - expected_offset).to_string(),
                )
                .with_header(
                    "Content-Range",
                    format!(
                        "bytes {}-{}/{}",
                        expected_offset,
                        partial_end - 1,
                        payload.len()
                    ),
                ),
        ])?;
        let url = server.fastq_url()?;
        let client = Client::builder().build()?;
        let mut reader = RetryingHttpRead::new(
            client,
            url,
            1,
            Duration::from_millis(1),
            Duration::from_millis(2),
        );

        let error = reader
            .open_response(Some(expected_offset as u64))
            .expect_err("partial ranged suffix should fail");
        server.finish()?;
        assert_eq!(error.kind(), std::io::ErrorKind::InvalidData);
        assert!(error.to_string().contains("partial suffix"));
        Ok(())
    }

    #[test]
    fn retrying_http_read_rejects_ranged_content_length_mismatch() -> Result<()> {
        let payload = b"abcdefghijklmnopqrstuvwxyz".to_vec();
        let expected_offset = 10_usize;
        let server = TestServer::spawn(vec![
            Exchange::ranged(&payload, expected_offset).with_header(
                "Content-Length",
                (payload.len() - expected_offset - 1).to_string(),
            ),
        ])?;
        let url = server.fastq_url()?;
        let client = Client::builder().build()?;
        let mut reader = RetryingHttpRead::new(
            client,
            url,
            1,
            Duration::from_millis(1),
            Duration::from_millis(2),
        );

        let error = reader
            .open_response(Some(expected_offset as u64))
            .expect_err("ranged content-length mismatch should fail");
        server.finish()?;
        assert_eq!(error.kind(), std::io::ErrorKind::InvalidData);
        assert!(error.to_string().contains("Content-Length"));
        Ok(())
    }

    #[test]
    fn retrying_http_read_rejects_unexpected_content_encoding() -> Result<()> {
        let payload = b"abcdefghijklmnopqrstuvwxyz".to_vec();
        let server = TestServer::spawn(vec![
            Exchange::initial(&payload).with_header("Content-Encoding", "gzip"),
        ])?;
        let url = server.fastq_url()?;
        let client = Client::builder().build()?;
        let mut reader = RetryingHttpRead::new(
            client,
            url,
            1,
            Duration::from_millis(1),
            Duration::from_millis(2),
        );

        let mut output = Vec::new();
        let error = reader
            .read_to_end(&mut output)
            .expect_err("unexpected content encoding should fail");
        server.finish()?;

        assert_eq!(error.kind(), std::io::ErrorKind::InvalidData);
        assert!(error.to_string().contains("Content-Encoding"));
        Ok(())
    }

    #[test]
    fn retrying_http_read_does_not_retry_invalid_response_headers() -> Result<()> {
        let payload = b"abcdefghijklmnopqrstuvwxyz".to_vec();
        let server = TestServer::spawn(vec![
            Exchange::initial(&payload).with_header("Content-Encoding", "gzip"),
        ])?;
        let url = server.fastq_url()?;
        let client = Client::builder().build()?;
        let mut reader = RetryingHttpRead::new(
            client,
            url,
            3,
            Duration::from_millis(1),
            Duration::from_millis(2),
        );

        let mut output = Vec::new();
        let error = reader
            .read_to_end(&mut output)
            .expect_err("invalid response headers should fail without retrying");
        server.finish()?;

        assert_eq!(error.kind(), std::io::ErrorKind::InvalidData);
        Ok(())
    }

    #[test]
    fn retrying_http_read_rejects_changed_total_size() -> Result<()> {
        let payload = b"abcdefghijklmnopqrstuvwxyz".to_vec();
        let expected_offset = 10_usize;
        let server = TestServer::spawn(vec![Exchange::ranged(&payload, expected_offset)])?;
        let url = server.fastq_url()?;
        let client = Client::builder().build()?;
        let mut reader = RetryingHttpRead::new(
            client,
            url,
            1,
            Duration::from_millis(1),
            Duration::from_millis(2),
        );
        reader.expected_total_bytes = Some(payload.len() as u64 + 1);

        let error = reader
            .open_response(Some(expected_offset as u64))
            .expect_err("changed total size should fail");
        server.finish()?;
        assert!(error.to_string().contains("total size changed"));
        Ok(())
    }

    #[test]
    fn scripted_server_reports_request_range_mismatches() -> Result<()> {
        let payload = b"abcdefghijklmnopqrstuvwxyz".to_vec();
        let server = TestServer::spawn(vec![Exchange::ranged(&payload, 10)])?;
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

        fn fastq_url(&self) -> Result<FastqUrl> {
            Ok(FastqUrl(Url::parse(&format!(
                "http://{}{}",
                self.address, FASTQ_TARGET
            ))?))
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
