//! MGI/BGI/DNBSEQ adapter catalog entries.

use crate::adapter::catalog::{
    AdapterCatalog, AdapterEntry, DnaSequence, KitRef, SequencingPlatform, SourceRef,
};

const FASTP_SOURCE: SourceRef = SourceRef {
    label: "fastp MIT adapter catalog",
    url: Some("https://github.com/OpenGene/fastp/blob/master/src/knownadapters.h"),
    version: None,
    note: Some("Vendor corroboration is still desirable for these MGI/BGI entries."),
};

const MGI_DNBSEQ_KIT: KitRef = KitRef {
    vendor: "MGI/BGI",
    name: "DNBSEQ adapter catalog",
    version: None,
};

pub(crate) static MGI_DNBSEQ: AdapterCatalog = AdapterCatalog {
    id: "mgi-dnbseq",
    name: "MGI/BGI/DNBSEQ",
    entries: &[
        AdapterEntry::three_prime_adapter(
            "mgi-dnbseq-forward-adapter",
            "MGI/BGI forward adapter",
            DnaSequence::from_iupac_ascii(b"AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA"),
            SequencingPlatform::MgiBgi,
            &MGI_DNBSEQ_KIT,
            &FASTP_SOURCE,
        ),
        AdapterEntry::three_prime_adapter(
            "mgi-dnbseq-reverse-adapter",
            "MGI/BGI reverse adapter",
            DnaSequence::from_iupac_ascii(b"AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG"),
            SequencingPlatform::MgiBgi,
            &MGI_DNBSEQ_KIT,
            &FASTP_SOURCE,
        ),
    ],
};
