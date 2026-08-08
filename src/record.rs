//! Core record and accounting types shared across ingress, parsing, and output layers.

use std::{fmt, marker::PhantomData, path::Path};

use serde::Serialize;

use crate::error::{HeaderExcerpt, InternalError, PairConstructionError, Result};

const RECORD_DEBUG_HEADER_BYTES: usize = 96;
const RECORD_DEBUG_CONTENT_BYTES: usize = 32;

/// Borrowed view of a biological sequence record.
///
/// Implementors expose header, sequence, and optional quality bytes without forcing ownership.
pub trait SequenceRecordRef {
    /// Return the raw record header bytes without the leading `@` or `>` marker.
    fn header(&self) -> &[u8];

    /// Return the raw sequence bytes for the record.
    fn sequence(&self) -> &[u8];

    /// Return quality bytes when the record format carries them.
    fn quality(&self) -> Option<&[u8]>;
}

/// Upstream source information attached to one borrowed record.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum InputSource<'a> {
    /// Record streamed directly from ENA.
    Ena { accession: &'a str },
    /// Record loaded from one local FASTQ file.
    LocalSingle { input: &'a Path },
    /// Record loaded from one interleaved local paired FASTQ file.
    LocalInterleavedPaired { input: &'a Path },
    /// Record loaded from one of two local paired FASTQ files.
    LocalPaired { input1: &'a Path, input2: &'a Path },
}

/// Mate identity for records originating from paired-end ingress.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum MateSide {
    /// First mate in a pair.
    Left,
    /// Second mate in a pair.
    Right,
}

impl fmt::Display for MateSide {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Left => formatter.write_str("left"),
            Self::Right => formatter.write_str("right"),
        }
    }
}

/// Canonical borrowed record view used by preprocessing plans.
#[derive(Clone, Copy)]
pub struct RecordView<'a> {
    header: &'a [u8],
    sequence: &'a [u8],
    quality: &'a [u8],
    source: Option<InputSource<'a>>,
}

mod debug {
    use std::fmt::{self, Write as _};

    pub(super) struct BytesPreview<'bytes> {
        bytes: &'bytes [u8],
        limit: usize,
    }

    impl<'bytes> BytesPreview<'bytes> {
        pub(super) const fn new(bytes: &'bytes [u8], limit: usize) -> Self {
            Self { bytes, limit }
        }
    }

    impl fmt::Debug for BytesPreview<'_> {
        fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
            formatter.write_char('"')?;
            for byte in self.bytes.iter().take(self.limit) {
                for escaped in std::ascii::escape_default(*byte) {
                    formatter.write_char(char::from(escaped))?;
                }
            }
            if self.bytes.len() > self.limit {
                formatter.write_str("…")?;
            }
            formatter.write_char('"')
        }
    }
}

impl fmt::Debug for RecordView<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("RecordView")
            .field(
                "header",
                &debug::BytesPreview::new(self.header, RECORD_DEBUG_HEADER_BYTES),
            )
            .field(
                "sequence",
                &debug::BytesPreview::new(self.sequence, RECORD_DEBUG_CONTENT_BYTES),
            )
            .field("sequence_len", &self.sequence.len())
            .field(
                "quality",
                &debug::BytesPreview::new(self.quality, RECORD_DEBUG_CONTENT_BYTES),
            )
            .field("quality_len", &self.quality.len())
            .field("source", &self.source)
            .finish()
    }
}

/// Positional claim that a record is the left mate of a pair.
#[derive(Clone, Copy, Debug)]
pub(crate) struct LeftMate;

/// Positional claim that a record is the right mate of a pair.
#[derive(Clone, Copy, Debug)]
pub(crate) struct RightMate;

/// Borrowed record carrying one typed positional mate claim.
#[repr(transparent)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct MateRecord<'record, Mate> {
    record: RecordView<'record>,
    _mate: PhantomData<Mate>,
}

impl<'record> From<RecordView<'record>> for MateRecord<'record, LeftMate> {
    fn from(record: RecordView<'record>) -> Self {
        Self {
            record,
            _mate: PhantomData,
        }
    }
}

impl<'record> From<RecordView<'record>> for MateRecord<'record, RightMate> {
    fn from(record: RecordView<'record>) -> Self {
        Self {
            record,
            _mate: PhantomData,
        }
    }
}

impl<'record, Mate> MateRecord<'record, Mate> {
    /// Return the claimed record without its positional marker.
    pub(crate) fn into_record(self) -> RecordView<'record> {
        self.record
    }
}

/// Two borrowed records whose mate relationship was established at ingress.
#[derive(Clone, Copy, Debug)]
pub(crate) struct RecordPair<'record> {
    left: RecordView<'record>,
    right: RecordView<'record>,
}

impl<'record> RecordPair<'record> {
    /// Validate two positional mate claims and establish their pair relationship.
    ///
    /// Successful construction performs no heap allocation. Diagnostic evidence for a rejected
    /// pair is copied into bounded inline storage.
    #[allow(
        clippy::result_large_err,
        reason = "rare pair failures retain bounded inline evidence without heap allocation"
    )]
    pub(crate) fn try_new(
        left: MateRecord<'record, LeftMate>,
        right: MateRecord<'record, RightMate>,
    ) -> std::result::Result<Self, PairConstructionError> {
        let left = left.into_record();
        let right = right.into_record();

        let left_sequence_len = left.sequence().len();
        let left_quality_len = left.quality().len();
        let right_sequence_len = right.sequence().len();
        let right_quality_len = right.quality().len();

        match (
            left_sequence_len != left_quality_len,
            right_sequence_len != right_quality_len,
        ) {
            (true, true) => {
                return Err(PairConstructionError::BothRecordLengths {
                    left_header: HeaderExcerpt::capture(left.header()),
                    left_sequence_len,
                    left_quality_len,
                    right_header: HeaderExcerpt::capture(right.header()),
                    right_sequence_len,
                    right_quality_len,
                });
            }
            (true, false) => {
                return Err(PairConstructionError::LeftRecordLength {
                    header: HeaderExcerpt::capture(left.header()),
                    sequence_len: left_sequence_len,
                    quality_len: left_quality_len,
                });
            }
            (false, true) => {
                return Err(PairConstructionError::RightRecordLength {
                    header: HeaderExcerpt::capture(right.header()),
                    sequence_len: right_sequence_len,
                    quality_len: right_quality_len,
                });
            }
            (false, false) => {}
        }

        let left_token = left
            .header()
            .split(u8::is_ascii_whitespace)
            .next()
            .unwrap_or(left.header());
        let right_token = right
            .header()
            .split(u8::is_ascii_whitespace)
            .next()
            .unwrap_or(right.header());

        match (left_token.ends_with(b"/2"), right_token.ends_with(b"/1")) {
            (true, true) => {
                return Err(PairConstructionError::BothHeadersContradictPositions {
                    left_header: HeaderExcerpt::capture(left.header()),
                    right_header: HeaderExcerpt::capture(right.header()),
                });
            }
            (true, false) => {
                return Err(PairConstructionError::LeftHeaderClaimsRight {
                    header: HeaderExcerpt::capture(left.header()),
                });
            }
            (false, true) => {
                return Err(PairConstructionError::RightHeaderClaimsLeft {
                    header: HeaderExcerpt::capture(right.header()),
                });
            }
            (false, false) => {}
        }

        if left.pair_key() != right.pair_key() {
            return Err(PairConstructionError::IdentifierMismatch {
                left_header: HeaderExcerpt::capture(left.header()),
                right_header: HeaderExcerpt::capture(right.header()),
            });
        }

        Ok(Self { left, right })
    }

    /// Return the left record.
    pub(crate) fn left(self) -> RecordView<'record> {
        self.left
    }

    /// Return the right record.
    pub(crate) fn right(self) -> RecordView<'record> {
        self.right
    }

    /// Replace both records while preserving their established ingress relationship.
    #[allow(
        clippy::unused_self,
        reason = "consuming an established pair proves replacement preserves its ingress relationship"
    )]
    pub(crate) const fn with_records(
        self,
        left: RecordView<'record>,
        right: RecordView<'record>,
    ) -> Self {
        Self { left, right }
    }
}

impl SequenceRecordRef for RecordView<'_> {
    fn header(&self) -> &[u8] {
        self.header
    }

    fn sequence(&self) -> &[u8] {
        self.sequence
    }

    fn quality(&self) -> Option<&[u8]> {
        Some(self.quality)
    }
}

impl<'a> RecordView<'a> {
    /// Construct a new borrowed FASTQ record view with no attached provenance.
    pub fn new(header: &'a [u8], sequence: &'a [u8], quality: &'a [u8]) -> Self {
        Self {
            header,
            sequence,
            quality,
            source: None,
        }
    }

    /// Attach the upstream input source to this record view.
    pub fn with_source(mut self, source: InputSource<'a>) -> Self {
        self.source = Some(source);
        self
    }

    /// Inherit the upstream input source from another record.
    pub fn inherit_source(mut self, record: Self) -> Self {
        self.source = record.source;
        self
    }

    /// Return the raw record header bytes without the leading `@` marker.
    pub fn header(&self) -> &'a [u8] {
        self.header
    }

    /// Return the sequence bytes.
    pub fn sequence(&self) -> &'a [u8] {
        self.sequence
    }

    /// Return the quality bytes.
    pub fn quality(&self) -> &'a [u8] {
        self.quality
    }

    /// Return the upstream input source when available.
    pub fn source(&self) -> Option<InputSource<'a>> {
        self.source
    }

    /// Return a new record view with updated sequence and quality slices while preserving metadata.
    ///
    /// # Errors
    ///
    /// Returns an error when the provided sequence and quality slices are not the same length.
    pub fn with_sequence_and_quality(self, sequence: &'a [u8], quality: &'a [u8]) -> Result<Self> {
        if sequence.len() != quality.len() {
            return Err(InternalError::ReplacementLength {
                header: String::from_utf8_lossy(self.header).into_owned(),
                sequence_len: sequence.len(),
                quality_len: quality.len(),
            }
            .into());
        }

        Ok(Self {
            header: self.header,
            sequence,
            quality,
            source: self.source,
        })
    }

    pub(crate) fn pair_key(&self) -> &'a [u8] {
        let first_token = self
            .header()
            .split(u8::is_ascii_whitespace)
            .next()
            .unwrap_or(self.header());
        match first_token {
            [prefix @ .., b'/', b'1' | b'2'] => prefix,
            _ => first_token,
        }
    }

    pub(crate) fn source_display(&self) -> String {
        match self.source() {
            Some(InputSource::Ena { accession }) => format!("ena:{accession}"),
            Some(InputSource::LocalSingle { input }) => format!("local:{}", input.display()),
            Some(InputSource::LocalInterleavedPaired { input }) => {
                format!("local-interleaved:{}", input.display())
            }
            Some(InputSource::LocalPaired { input1, input2 }) => {
                format!("local-paired:{}|{}", input1.display(), input2.display())
            }
            None => "unknown".to_owned(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{InputSource, LeftMate, MateRecord, RecordPair, RecordView, RightMate};
    use crate::error::PairConstructionError;

    #[test]
    fn pair_key_strips_slash_mate_suffix() {
        let left = RecordView::new(b"read123/1", b"A", b"I");
        let right = RecordView::new(b"read123/2", b"A", b"I");
        let bare = RecordView::new(b"read123", b"A", b"I");

        assert_eq!(left.pair_key(), b"read123");
        assert_eq!(right.pair_key(), b"read123");
        assert_eq!(bare.pair_key(), b"read123");
    }

    fn pair_error(
        left_header: &'static [u8],
        left_sequence: &'static [u8],
        left_quality: &'static [u8],
        right_header: &'static [u8],
        right_sequence: &'static [u8],
        right_quality: &'static [u8],
    ) -> PairConstructionError {
        RecordPair::try_new(
            RecordView::new(left_header, left_sequence, left_quality).into(),
            RecordView::new(right_header, right_sequence, right_quality).into(),
        )
        .expect_err("fixture should fail pair construction")
    }

    #[test]
    fn pair_construction_rejects_left_record_length() {
        assert!(matches!(
            pair_error(b"read/1", b"AA", b"I", b"read/2", b"TT", b"II"),
            PairConstructionError::LeftRecordLength { .. }
        ));
    }

    #[test]
    fn pair_construction_rejects_right_record_length() {
        assert!(matches!(
            pair_error(b"read/1", b"AA", b"II", b"read/2", b"TT", b"I"),
            PairConstructionError::RightRecordLength { .. }
        ));
    }

    #[test]
    fn pair_construction_rejects_both_record_lengths() {
        assert!(matches!(
            pair_error(b"read/1", b"AA", b"I", b"read/2", b"TT", b"I"),
            PairConstructionError::BothRecordLengths { .. }
        ));
    }

    #[test]
    fn pair_construction_rejects_left_header_claiming_right() {
        assert!(matches!(
            pair_error(b"read/2", b"A", b"I", b"read/2", b"T", b"I"),
            PairConstructionError::LeftHeaderClaimsRight { .. }
        ));
    }

    #[test]
    fn pair_construction_rejects_right_header_claiming_left() {
        assert!(matches!(
            pair_error(b"read/1", b"A", b"I", b"read/1", b"T", b"I"),
            PairConstructionError::RightHeaderClaimsLeft { .. }
        ));
    }

    #[test]
    fn pair_construction_rejects_both_headers_contradicting_positions() {
        assert!(matches!(
            pair_error(b"read/2", b"A", b"I", b"read/1", b"T", b"I"),
            PairConstructionError::BothHeadersContradictPositions { .. }
        ));
    }

    #[test]
    fn pair_construction_rejects_identifier_mismatch() {
        assert!(matches!(
            pair_error(b"left/1", b"A", b"I", b"right/2", b"T", b"I"),
            PairConstructionError::IdentifierMismatch { .. }
        ));
    }

    #[test]
    fn pair_construction_validates_lengths_before_headers_and_identifiers() {
        assert!(matches!(
            pair_error(b"left/2", b"AA", b"I", b"right/1", b"TT", b"II"),
            PairConstructionError::LeftRecordLength { .. }
        ));
    }

    #[test]
    fn pair_construction_validates_headers_before_identifiers() {
        assert!(matches!(
            pair_error(b"left/2", b"A", b"I", b"right/2", b"T", b"I"),
            PairConstructionError::LeftHeaderClaimsRight { .. }
        ));
    }

    #[test]
    fn typed_mate_claims_add_no_storage_or_drop_requirement() {
        assert_eq!(
            std::mem::size_of::<MateRecord<'_, LeftMate>>(),
            std::mem::size_of::<RecordView<'_>>()
        );
        assert_eq!(
            std::mem::size_of::<MateRecord<'_, RightMate>>(),
            std::mem::size_of::<RecordView<'_>>()
        );
        assert!(!std::mem::needs_drop::<MateRecord<'_, LeftMate>>());
        assert!(!std::mem::needs_drop::<MateRecord<'_, RightMate>>());
    }

    #[test]
    fn with_sequence_and_quality_preserves_source() {
        let record =
            RecordView::new(b"read123/1", b"ACGT", b"IIII").with_source(InputSource::Ena {
                accession: "ERR000002",
            });

        let rewritten = record
            .with_sequence_and_quality(b"AC", b"II")
            .expect("replacement slices should be valid");

        assert_eq!(rewritten.source(), record.source());
    }
}
