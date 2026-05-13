//! `QIAseq` adapter catalog entries.

use crate::adapter::catalog::{
    AdapterCatalog, AdapterEntry, DnaSequence, KitRef, SequencingPlatform, SourceRef,
};

const FASTP_SOURCE: SourceRef = SourceRef {
    label: "fastp MIT adapter catalog",
    url: Some("https://github.com/OpenGene/fastp/blob/master/src/knownadapters.h"),
    version: None,
    note: Some("Official QIAGEN corroboration is still desirable for this entry."),
};

const QIASEQ_MIRNA_KIT: KitRef = KitRef {
    vendor: "QIAGEN",
    name: "QIAseq miRNA",
    version: None,
};

pub(crate) static QIASEQ_MIRNA: AdapterCatalog = AdapterCatalog {
    id: "qiaseq-mirna",
    name: "QIAseq miRNA",
    entries: &[AdapterEntry::three_prime_adapter(
        "qiaseq-mirna-adapter",
        "QIAseq miRNA adapter",
        DnaSequence::from_iupac_ascii(b"AACTGTAGGCACCATCAAT"),
        SequencingPlatform::QiaSeq,
        &QIASEQ_MIRNA_KIT,
        &FASTP_SOURCE,
    )],
};
