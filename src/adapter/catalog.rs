//! Structured metadata for curated adapter sequence catalogs.

use clap::ValueEnum;

use super::catalogs;

/// Adapter trimming presets exposed by the CLI.
#[derive(Clone, Copy, Debug, Eq, PartialEq, ValueEnum)]
pub enum AdapterPreset {
    /// Do not trim adapters.
    None,
    /// Trim Illumina `TruSeq` adapter overlap.
    #[value(name = "illumina-truseq", alias = "illumina-tru-seq")]
    IlluminaTruSeq,
    /// Trim Illumina Nextera/tagmentation adapter overlap.
    #[value(name = "illumina-nextera", alias = "illumina-tagmentation")]
    IlluminaNextera,
    /// Trim Illumina `TruSeq` small RNA adapter overlap.
    #[value(name = "illumina-truseq-small-rna", alias = "illumina-small-rna")]
    IlluminaTruSeqSmallRna,
    /// Trim Illumina stranded RNA ligation adapter overlap.
    #[value(name = "illumina-stranded-rna-ligation")]
    IlluminaStrandedRnaLigation,
    /// Trim Illumina `ScriptSeq`/`TruSeq` methylation adapter overlap.
    #[value(name = "illumina-scriptseq-truseq-methylation")]
    IlluminaScriptSeqTruSeqMethylation,
    /// Trim MGI/BGI/DNBSEQ short-read adapter overlap.
    #[value(name = "mgi-dnbseq", alias = "mgi-bgi", alias = "bgi-dnbseq")]
    MgiDnbSeq,
    /// Trim `QIAseq` miRNA adapter overlap.
    #[value(name = "qiaseq-mirna", alias = "qiaseq")]
    QiaSeqMirna,
}

impl AdapterPreset {
    /// Return the curated adapter catalog for presets that trim adapters.
    pub(crate) const fn catalog(self) -> Option<&'static AdapterCatalog> {
        match self {
            Self::None => None,
            Self::IlluminaTruSeq => Some(&catalogs::illumina::ILLUMINA_TRUSEQ),
            Self::IlluminaNextera => Some(&catalogs::illumina::ILLUMINA_NEXTERA),
            Self::IlluminaTruSeqSmallRna => Some(&catalogs::illumina::ILLUMINA_TRUSEQ_SMALL_RNA),
            Self::IlluminaStrandedRnaLigation => {
                Some(&catalogs::illumina::ILLUMINA_STRANDED_RNA_LIGATION)
            }
            Self::IlluminaScriptSeqTruSeqMethylation => {
                Some(&catalogs::illumina::ILLUMINA_SCRIPTSEQ_TRUSEQ_METHYLATION)
            }
            Self::MgiDnbSeq => Some(&catalogs::mgi::MGI_DNBSEQ),
            Self::QiaSeqMirna => Some(&catalogs::qiaseq::QIASEQ_MIRNA),
        }
    }

    #[cfg(test)]
    const CATALOG_PRESETS: &'static [Self] = &[
        Self::IlluminaTruSeq,
        Self::IlluminaNextera,
        Self::IlluminaTruSeqSmallRna,
        Self::IlluminaStrandedRnaLigation,
        Self::IlluminaScriptSeqTruSeqMethylation,
        Self::MgiDnbSeq,
        Self::QiaSeqMirna,
    ];
}

/// Curated sequence catalog selected by one CLI preset.
#[derive(Debug, Eq, PartialEq)]
pub(crate) struct AdapterCatalog {
    /// Stable machine-readable catalog ID.
    pub(crate) id: &'static str,
    /// Human-readable catalog name.
    pub(crate) name: &'static str,
    /// Curated entries in this catalog.
    pub(crate) entries: &'static [AdapterEntry],
}

/// Uppercase DNA or IUPAC ambiguity-code sequence curated from a catalog source.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct DnaSequence(&'static [u8]);

impl DnaSequence {
    /// Build a static DNA sequence, failing const evaluation for invalid bases.
    pub(crate) const fn from_iupac_ascii(sequence: &'static [u8]) -> Self {
        assert!(
            !sequence.is_empty(),
            "catalog DNA sequence must not be empty"
        );

        let mut index = 0;
        while index < sequence.len() {
            match sequence[index] {
                b'A' | b'C' | b'G' | b'T' | b'U' | b'R' | b'Y' | b'S' | b'W' | b'K' | b'M'
                | b'B' | b'D' | b'H' | b'V' | b'N' => {}
                _ => panic!("catalog DNA sequence must use uppercase IUPAC bases"),
            }
            index += 1;
        }

        Self(sequence)
    }

    /// Return sequence bytes for matching.
    pub(crate) const fn as_bytes(self) -> &'static [u8] {
        self.0
    }
}

impl AdapterCatalog {
    /// Return entries that the current 3' suffix-overlap trimmer may consume.
    pub(crate) fn current_suffix_trim_entries(
        &self,
    ) -> impl Iterator<Item = &'static AdapterEntry> {
        self.entries.iter().filter(|entry| {
            entry.role == SequenceRole::Adapter && entry.geometry == MatchGeometry::ThreePrimeSuffix
        })
    }
}

/// One curated adapter, primer, barcode, flank, or marker sequence.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct AdapterEntry {
    /// Stable machine-readable entry ID.
    pub(crate) id: &'static str,
    /// Human-readable entry name.
    pub(crate) name: &'static str,
    /// Sequence bytes as published by the source.
    pub(crate) sequence: DnaSequence,
    /// Biological or workflow role of this sequence.
    pub(crate) role: SequenceRole,
    /// Sequencing platform or vendor ecosystem.
    pub(crate) platform: SequencingPlatform,
    /// Kit or protocol family this entry belongs to.
    pub(crate) kit: &'static KitRef,
    /// Matching geometry required to interpret this sequence correctly.
    pub(crate) geometry: MatchGeometry,
    /// Source used for this curated entry.
    pub(crate) source: &'static SourceRef,
}

#[allow(dead_code)]
impl AdapterEntry {
    /// Build an adapter entry that the current 3' suffix-overlap trimmer may consume.
    pub(crate) const fn three_prime_adapter(
        id: &'static str,
        name: &'static str,
        sequence: DnaSequence,
        platform: SequencingPlatform,
        kit: &'static KitRef,
        source: &'static SourceRef,
    ) -> Self {
        Self::new(
            id,
            name,
            sequence,
            SequenceRole::Adapter,
            platform,
            kit,
            MatchGeometry::ThreePrimeSuffix,
            source,
        )
    }

    /// Build an adapter entry that needs both-end matching before it can be interpreted.
    pub(crate) const fn both_end_adapter(
        id: &'static str,
        name: &'static str,
        sequence: DnaSequence,
        platform: SequencingPlatform,
        kit: &'static KitRef,
        source: &'static SourceRef,
    ) -> Self {
        Self::new(
            id,
            name,
            sequence,
            SequenceRole::Adapter,
            platform,
            kit,
            MatchGeometry::BothEnds,
            source,
        )
    }

    /// Build an adapter entry retained only as curated metadata for now.
    pub(crate) const fn catalog_adapter(
        id: &'static str,
        name: &'static str,
        sequence: DnaSequence,
        platform: SequencingPlatform,
        kit: &'static KitRef,
        source: &'static SourceRef,
    ) -> Self {
        Self::new(
            id,
            name,
            sequence,
            SequenceRole::Adapter,
            platform,
            kit,
            MatchGeometry::CatalogOnly,
            source,
        )
    }

    /// Build an adapter entry that marks an internal sequence requiring future logic.
    pub(crate) const fn internal_adapter(
        id: &'static str,
        name: &'static str,
        sequence: DnaSequence,
        platform: SequencingPlatform,
        kit: &'static KitRef,
        source: &'static SourceRef,
    ) -> Self {
        Self::new(
            id,
            name,
            sequence,
            SequenceRole::Adapter,
            platform,
            kit,
            MatchGeometry::InternalMarker,
            source,
        )
    }

    /// Build a 5' primer entry retained for future primer-aware processing.
    pub(crate) const fn five_prime_primer(
        id: &'static str,
        name: &'static str,
        sequence: DnaSequence,
        platform: SequencingPlatform,
        kit: &'static KitRef,
        source: &'static SourceRef,
    ) -> Self {
        Self::new(
            id,
            name,
            sequence,
            SequenceRole::Primer,
            platform,
            kit,
            MatchGeometry::FivePrimePrefix,
            source,
        )
    }

    /// Build a 3' primer entry; these remain inert for the current adapter-only trimmer.
    pub(crate) const fn three_prime_primer(
        id: &'static str,
        name: &'static str,
        sequence: DnaSequence,
        platform: SequencingPlatform,
        kit: &'static KitRef,
        source: &'static SourceRef,
    ) -> Self {
        Self::new(
            id,
            name,
            sequence,
            SequenceRole::Primer,
            platform,
            kit,
            MatchGeometry::ThreePrimeSuffix,
            source,
        )
    }

    /// Build a catalog-only barcode entry.
    pub(crate) const fn catalog_barcode(
        id: &'static str,
        name: &'static str,
        sequence: DnaSequence,
        platform: SequencingPlatform,
        kit: &'static KitRef,
        source: &'static SourceRef,
    ) -> Self {
        Self::new(
            id,
            name,
            sequence,
            SequenceRole::Barcode,
            platform,
            kit,
            MatchGeometry::CatalogOnly,
            source,
        )
    }

    /// Build a catalog-only flank entry.
    pub(crate) const fn catalog_flank(
        id: &'static str,
        name: &'static str,
        sequence: DnaSequence,
        platform: SequencingPlatform,
        kit: &'static KitRef,
        source: &'static SourceRef,
    ) -> Self {
        Self::new(
            id,
            name,
            sequence,
            SequenceRole::Flank,
            platform,
            kit,
            MatchGeometry::CatalogOnly,
            source,
        )
    }

    #[expect(
        clippy::too_many_arguments,
        reason = "private funnel keeps role-specific catalog constructors consistent"
    )]
    const fn new(
        id: &'static str,
        name: &'static str,
        sequence: DnaSequence,
        role: SequenceRole,
        platform: SequencingPlatform,
        kit: &'static KitRef,
        geometry: MatchGeometry,
        source: &'static SourceRef,
    ) -> Self {
        Self {
            id,
            name,
            sequence,
            role,
            platform,
            kit,
            geometry,
            source,
        }
    }
}

/// Biological or workflow role of one catalog sequence.
#[allow(dead_code)]
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum SequenceRole {
    /// Library adapter sequence.
    Adapter,
    /// Primer sequence.
    Primer,
    /// Sample barcode sequence.
    Barcode,
    /// Sequence flanking a barcode or other marker.
    Flank,
}

/// Sequencing platform or vendor ecosystem.
#[allow(dead_code)]
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum SequencingPlatform {
    /// Illumina sequencing chemistry.
    Illumina,
    /// MGI/BGI/DNBSEQ short-read sequencing chemistry.
    MgiBgi,
    /// `QIAseq` library-preparation chemistry.
    QiaSeq,
    /// `PacBio` long-read sequencing chemistry.
    PacBio,
    /// Oxford Nanopore Technologies sequencing chemistry.
    OxfordNanopore,
}

/// Matching geometry required by a sequence entry.
#[allow(dead_code)]
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum MatchGeometry {
    /// Read 3' suffix overlaps the catalog sequence prefix.
    ThreePrimeSuffix,
    /// Catalog sequence is expected at the read 5' prefix.
    FivePrimePrefix,
    /// Catalog sequence may appear at either read end.
    BothEnds,
    /// Catalog sequence marks an internal region requiring future split/reject logic.
    InternalMarker,
    /// Cataloged metadata that the current trimmer must not consume.
    CatalogOnly,
}

/// Kit or protocol family for a catalog entry.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct KitRef {
    /// Vendor or project name.
    pub(crate) vendor: &'static str,
    /// Kit or protocol name.
    pub(crate) name: &'static str,
    /// Kit/protocol version when relevant.
    pub(crate) version: Option<&'static str>,
}

/// Provenance for a catalog entry.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct SourceRef {
    /// Short source label.
    pub(crate) label: &'static str,
    /// Source URL when available.
    pub(crate) url: Option<&'static str>,
    /// Source version/date when available.
    pub(crate) version: Option<&'static str>,
    /// Additional provenance notes.
    pub(crate) note: Option<&'static str>,
}

#[cfg(test)]
mod tests {
    use super::{AdapterPreset, DnaSequence, MatchGeometry, SequenceRole, catalogs};

    #[test]
    fn all_non_none_presets_have_catalog_entries() {
        for preset in AdapterPreset::CATALOG_PRESETS {
            let catalog = preset.catalog().expect("preset should have a catalog");
            assert!(
                !catalog.entries.is_empty(),
                "{} should not be empty",
                catalog.id
            );
        }
    }

    #[test]
    fn none_preset_has_no_catalog() {
        assert!(AdapterPreset::None.catalog().is_none());
    }

    #[test]
    fn dna_sequence_accepts_uppercase_iupac_bases() {
        assert_eq!(
            DnaSequence::from_iupac_ascii(b"ACGTURYSWKMBDHVN").as_bytes(),
            b"ACGTURYSWKMBDHVN"
        );
    }

    #[test]
    #[should_panic(expected = "catalog DNA sequence must use uppercase IUPAC bases")]
    fn dna_sequence_rejects_non_iupac_bases() {
        DnaSequence::from_iupac_ascii(b"ACGTX");
    }

    #[test]
    fn current_suffix_trim_entries_are_only_supported_adapters() {
        let catalog = &catalogs::ont::ONT_NATIVE_BARCODING_V14;

        for entry in catalog.current_suffix_trim_entries() {
            assert_eq!(entry.role, SequenceRole::Adapter);
            assert_eq!(entry.geometry, MatchGeometry::ThreePrimeSuffix);
        }
    }

    #[test]
    fn ont_native_barcodes_are_not_current_suffix_trim_targets() {
        let catalog = &catalogs::ont::ONT_NATIVE_BARCODING_V14;
        assert!(
            catalog
                .entries
                .iter()
                .any(|entry| entry.role == SequenceRole::Barcode)
        );
        assert!(
            catalog
                .current_suffix_trim_entries()
                .all(|entry| entry.role != SequenceRole::Barcode)
        );
    }
}
