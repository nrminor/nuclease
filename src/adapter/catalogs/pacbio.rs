//! `PacBio` primer catalog entries.
#![allow(dead_code)]

use crate::adapter::catalog::{
    AdapterCatalog, AdapterEntry, DnaSequence, KitRef, SequencingPlatform, SourceRef,
};

const PACBIO_ISOSEQ_SOURCE: SourceRef = SourceRef {
    label: "PacBio Iso-Seq documentation",
    url: Some("https://isoseq.how/clustering/cli-workflow.html"),
    version: None,
    note: Some("Primer catalogs are not a replacement for lima/isoseq primer classification."),
};

const PACBIO_ISOSEQ_KIT: KitRef = KitRef {
    vendor: "PacBio",
    name: "Iso-Seq",
    version: None,
};

const PACBIO_MAS_ISOSEQ_KIT: KitRef = KitRef {
    vendor: "PacBio",
    name: "MAS-Seq/IsoSeqX",
    version: None,
};

pub(crate) static PACBIO_ISOSEQ: AdapterCatalog = AdapterCatalog {
    id: "pacbio-isoseq",
    name: "PacBio Iso-Seq primers",
    entries: &[
        AdapterEntry::five_prime_primer(
            "pacbio-isoseq-neb-5p-primer",
            "PacBio Iso-Seq NEB 5' primer",
            DnaSequence::from_iupac_ascii(b"GCAATGAAGTCGCAGGGTTGGG"),
            SequencingPlatform::PacBio,
            &PACBIO_ISOSEQ_KIT,
            &PACBIO_ISOSEQ_SOURCE,
        ),
        AdapterEntry::five_prime_primer(
            "pacbio-isoseq-clontech-5p-primer",
            "PacBio Iso-Seq Clontech 5' primer",
            DnaSequence::from_iupac_ascii(b"AAGCAGTGGTATCAACGCAGAGTACATGGGG"),
            SequencingPlatform::PacBio,
            &PACBIO_ISOSEQ_KIT,
            &PACBIO_ISOSEQ_SOURCE,
        ),
        AdapterEntry::three_prime_primer(
            "pacbio-isoseq-neb-clontech-3p-primer",
            "PacBio Iso-Seq NEB/Clontech 3' primer",
            DnaSequence::from_iupac_ascii(b"GTACTCTGCGTTGATACCACTGCTT"),
            SequencingPlatform::PacBio,
            &PACBIO_ISOSEQ_KIT,
            &PACBIO_ISOSEQ_SOURCE,
        ),
    ],
};

pub(crate) static PACBIO_MAS_ISOSEQ: AdapterCatalog = AdapterCatalog {
    id: "pacbio-mas-isoseq",
    name: "PacBio MAS-Seq/IsoSeqX primers",
    entries: &[
        mas_5p(
            "pacbio-mas-isoseq-bc01-5p-primer",
            "PacBio MAS-Seq BC01 5' primer",
            DnaSequence::from_iupac_ascii(b"CTACACGACGCTCTTCCGATCTACTACACGCAATGAAGTCGCAGGGTTGGG"),
        ),
        mas_5p(
            "pacbio-mas-isoseq-bc02-5p-primer",
            "PacBio MAS-Seq BC02 5' primer",
            DnaSequence::from_iupac_ascii(b"CTACACGACGCTCTTCCGATCTACTAGTAGCAATGAAGTCGCAGGGTTGGG"),
        ),
        mas_5p(
            "pacbio-mas-isoseq-bc03-5p-primer",
            "PacBio MAS-Seq BC03 5' primer",
            DnaSequence::from_iupac_ascii(b"CTACACGACGCTCTTCCGATCTAGTGTACGCAATGAAGTCGCAGGGTTGGG"),
        ),
        mas_5p(
            "pacbio-mas-isoseq-bc04-5p-primer",
            "PacBio MAS-Seq BC04 5' primer",
            DnaSequence::from_iupac_ascii(b"CTACACGACGCTCTTCCGATCTATCACTAGCAATGAAGTCGCAGGGTTGGG"),
        ),
        mas_5p(
            "pacbio-mas-isoseq-bc05-5p-primer",
            "PacBio MAS-Seq BC05 5' primer",
            DnaSequence::from_iupac_ascii(b"CTACACGACGCTCTTCCGATCTCAGCTGTGCAATGAAGTCGCAGGGTTGGG"),
        ),
        mas_5p(
            "pacbio-mas-isoseq-bc06-5p-primer",
            "PacBio MAS-Seq BC06 5' primer",
            DnaSequence::from_iupac_ascii(b"CTACACGACGCTCTTCCGATCTCAGTCACGCAATGAAGTCGCAGGGTTGGG"),
        ),
        mas_5p(
            "pacbio-mas-isoseq-bc07-5p-primer",
            "PacBio MAS-Seq BC07 5' primer",
            DnaSequence::from_iupac_ascii(b"CTACACGACGCTCTTCCGATCTCATGTATGCAATGAAGTCGCAGGGTTGGG"),
        ),
        mas_5p(
            "pacbio-mas-isoseq-bc08-5p-primer",
            "PacBio MAS-Seq BC08 5' primer",
            DnaSequence::from_iupac_ascii(b"CTACACGACGCTCTTCCGATCTCGTATGTGCAATGAAGTCGCAGGGTTGGG"),
        ),
        mas_5p(
            "pacbio-mas-isoseq-bc09-5p-primer",
            "PacBio MAS-Seq BC09 5' primer",
            DnaSequence::from_iupac_ascii(b"CTACACGACGCTCTTCCGATCTGACATGTGCAATGAAGTCGCAGGGTTGGG"),
        ),
        mas_5p(
            "pacbio-mas-isoseq-bc10-5p-primer",
            "PacBio MAS-Seq BC10 5' primer",
            DnaSequence::from_iupac_ascii(b"CTACACGACGCTCTTCCGATCTGAGTCTAGCAATGAAGTCGCAGGGTTGGG"),
        ),
        mas_5p(
            "pacbio-mas-isoseq-bc11-5p-primer",
            "PacBio MAS-Seq BC11 5' primer",
            DnaSequence::from_iupac_ascii(b"CTACACGACGCTCTTCCGATCTGTAGATAGCAATGAAGTCGCAGGGTTGGG"),
        ),
        mas_5p(
            "pacbio-mas-isoseq-bc12-5p-primer",
            "PacBio MAS-Seq BC12 5' primer",
            DnaSequence::from_iupac_ascii(b"CTACACGACGCTCTTCCGATCTGTATGACGCAATGAAGTCGCAGGGTTGGG"),
        ),
        AdapterEntry::three_prime_primer(
            "pacbio-mas-isoseq-3p-primer",
            "PacBio MAS-Seq/IsoSeqX 3' primer",
            DnaSequence::from_iupac_ascii(b"AAGCAGTGGTATCAACGCAGAGTAC"),
            SequencingPlatform::PacBio,
            &PACBIO_MAS_ISOSEQ_KIT,
            &PACBIO_ISOSEQ_SOURCE,
        ),
    ],
};

const fn mas_5p(id: &'static str, name: &'static str, sequence: DnaSequence) -> AdapterEntry {
    AdapterEntry::five_prime_primer(
        id,
        name,
        sequence,
        SequencingPlatform::PacBio,
        &PACBIO_MAS_ISOSEQ_KIT,
        &PACBIO_ISOSEQ_SOURCE,
    )
}
