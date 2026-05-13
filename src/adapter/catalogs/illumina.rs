//! Illumina adapter catalog entries.

use crate::adapter::catalog::{
    AdapterCatalog, AdapterEntry, DnaSequence, KitRef, SequencingPlatform, SourceRef,
};

const ILLUMINA_SOURCE: SourceRef = SourceRef {
    label: "Illumina adapter sequence documentation",
    url: Some(
        "https://knowledge.illumina.com/library-preparation/general/library-preparation-general-reference_material-list/000001314",
    ),
    version: None,
    note: Some("Official Illumina support documentation for adapter sequences."),
};

const TRUSEQ_KIT: KitRef = KitRef {
    vendor: "Illumina",
    name: "TruSeq",
    version: None,
};

const NEXTERA_KIT: KitRef = KitRef {
    vendor: "Illumina",
    name: "Nextera/tagmentation",
    version: None,
};

const SMALL_RNA_KIT: KitRef = KitRef {
    vendor: "Illumina",
    name: "TruSeq Small RNA",
    version: None,
};

const STRANDED_RNA_KIT: KitRef = KitRef {
    vendor: "Illumina",
    name: "Stranded RNA ligation",
    version: None,
};

const SCRIPTSEQ_METHYLATION_KIT: KitRef = KitRef {
    vendor: "Illumina",
    name: "ScriptSeq/TruSeq methylation",
    version: None,
};

pub(crate) static ILLUMINA_TRUSEQ: AdapterCatalog = AdapterCatalog {
    id: "illumina-truseq",
    name: "Illumina TruSeq",
    entries: &[
        AdapterEntry::three_prime_adapter(
            "illumina-truseq-r1",
            "Illumina TruSeq Read 1 adapter",
            DnaSequence::from_iupac_ascii(b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"),
            SequencingPlatform::Illumina,
            &TRUSEQ_KIT,
            &ILLUMINA_SOURCE,
        ),
        AdapterEntry::three_prime_adapter(
            "illumina-truseq-r2",
            "Illumina TruSeq Read 2 adapter",
            DnaSequence::from_iupac_ascii(b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"),
            SequencingPlatform::Illumina,
            &TRUSEQ_KIT,
            &ILLUMINA_SOURCE,
        ),
    ],
};

pub(crate) static ILLUMINA_NEXTERA: AdapterCatalog = AdapterCatalog {
    id: "illumina-nextera",
    name: "Illumina Nextera/tagmentation",
    entries: &[AdapterEntry::three_prime_adapter(
        "illumina-nextera-transposase",
        "Illumina Nextera transposase adapter",
        DnaSequence::from_iupac_ascii(b"CTGTCTCTTATACACATCT"),
        SequencingPlatform::Illumina,
        &NEXTERA_KIT,
        &ILLUMINA_SOURCE,
    )],
};

pub(crate) static ILLUMINA_TRUSEQ_SMALL_RNA: AdapterCatalog = AdapterCatalog {
    id: "illumina-truseq-small-rna",
    name: "Illumina TruSeq Small RNA",
    entries: &[AdapterEntry::three_prime_adapter(
        "illumina-truseq-small-rna-adapter",
        "Illumina TruSeq Small RNA adapter",
        DnaSequence::from_iupac_ascii(b"TGGAATTCTCGGGTGCCAAGG"),
        SequencingPlatform::Illumina,
        &SMALL_RNA_KIT,
        &ILLUMINA_SOURCE,
    )],
};

pub(crate) static ILLUMINA_STRANDED_RNA_LIGATION: AdapterCatalog = AdapterCatalog {
    id: "illumina-stranded-rna-ligation",
    name: "Illumina stranded RNA ligation",
    entries: &[AdapterEntry::three_prime_adapter(
        "illumina-stranded-rna-ligation-adapter",
        "Illumina stranded RNA ligation adapter",
        DnaSequence::from_iupac_ascii(b"ACTGTCTCTTATACACATCT"),
        SequencingPlatform::Illumina,
        &STRANDED_RNA_KIT,
        &ILLUMINA_SOURCE,
    )],
};

pub(crate) static ILLUMINA_SCRIPTSEQ_TRUSEQ_METHYLATION: AdapterCatalog = AdapterCatalog {
    id: "illumina-scriptseq-truseq-methylation",
    name: "Illumina ScriptSeq/TruSeq methylation",
    entries: &[
        AdapterEntry::three_prime_adapter(
            "illumina-scriptseq-truseq-methylation-r1",
            "Illumina ScriptSeq/TruSeq methylation Read 1 adapter",
            DnaSequence::from_iupac_ascii(b"AGATCGGAAGAGCACACGTCTGAAC"),
            SequencingPlatform::Illumina,
            &SCRIPTSEQ_METHYLATION_KIT,
            &ILLUMINA_SOURCE,
        ),
        AdapterEntry::three_prime_adapter(
            "illumina-scriptseq-truseq-methylation-r2",
            "Illumina ScriptSeq/TruSeq methylation Read 2 adapter",
            DnaSequence::from_iupac_ascii(b"AGATCGGAAGAGCGTCGTGTAGGGA"),
            SequencingPlatform::Illumina,
            &SCRIPTSEQ_METHYLATION_KIT,
            &ILLUMINA_SOURCE,
        ),
    ],
};
