//! Oxford Nanopore Technologies adapter and barcode catalog entries.
#![allow(dead_code)]

use crate::adapter::catalog::{
    AdapterCatalog, AdapterEntry, DnaSequence, KitRef, SequencingPlatform, SourceRef,
};

const ONT_CHEMISTRY_SOURCE: SourceRef = SourceRef {
    label: "Oxford Nanopore Technologies Chemistry Technical Document",
    url: Some("https://nanoporetech.com/document/chemistry-technical-document#barcode-sequences"),
    version: Some("CHTD_500_v1_revAT_17Feb2026"),
    note: Some(
        "Official ONT documentation; entries are cataloged for provenance, not full ONT preprocessing.",
    ),
};

const ONT_LIGATION_V14_KIT: KitRef = KitRef {
    vendor: "Oxford Nanopore Technologies",
    name: "Ligation/Native adapter",
    version: Some("V14"),
};

const ONT_RAPID_V14_KIT: KitRef = KitRef {
    vendor: "Oxford Nanopore Technologies",
    name: "Rapid adapter",
    version: Some("V14"),
};

const ONT_NATIVE_BARCODING_V14_KIT: KitRef = KitRef {
    vendor: "Oxford Nanopore Technologies",
    name: "Native Barcoding Kit",
    version: Some("V14"),
};

const ONT_CDNA_V14_KIT: KitRef = KitRef {
    vendor: "Oxford Nanopore Technologies",
    name: "cDNA/direct RNA adapters and primers",
    version: Some("V14"),
};

pub(crate) static ONT_LIGATION_V14: AdapterCatalog = AdapterCatalog {
    id: "ont-ligation-v14",
    name: "ONT ligation/native adapters V14",
    entries: &[
        AdapterEntry::both_end_adapter(
            "ont-v14-ligation-adapter-top",
            "ONT V14 ligation adapter top strand",
            DnaSequence::from_iupac_ascii(b"TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCT"),
            SequencingPlatform::OxfordNanopore,
            &ONT_LIGATION_V14_KIT,
            &ONT_CHEMISTRY_SOURCE,
        ),
        AdapterEntry::catalog_adapter(
            "ont-v14-ligation-adapter-bottom",
            "ONT V14 ligation adapter bottom strand",
            DnaSequence::from_iupac_ascii(b"GCAATACGTAACTGAACGAAGTACAGG"),
            SequencingPlatform::OxfordNanopore,
            &ONT_LIGATION_V14_KIT,
            &ONT_CHEMISTRY_SOURCE,
        ),
        AdapterEntry::catalog_adapter(
            "ont-v14-native-adapter-bottom",
            "ONT V14 native adapter bottom strand",
            DnaSequence::from_iupac_ascii(b"ACGTAACTGAACGAAGTACAGG"),
            SequencingPlatform::OxfordNanopore,
            &ONT_LIGATION_V14_KIT,
            &ONT_CHEMISTRY_SOURCE,
        ),
    ],
};

pub(crate) static ONT_RAPID_V14: AdapterCatalog = AdapterCatalog {
    id: "ont-rapid-v14",
    name: "ONT rapid adapter V14",
    entries: &[
        AdapterEntry::both_end_adapter(
            "ont-v14-rapid-adapter",
            "ONT V14 rapid adapter",
            DnaSequence::from_iupac_ascii(b"TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCT"),
            SequencingPlatform::OxfordNanopore,
            &ONT_RAPID_V14_KIT,
            &ONT_CHEMISTRY_SOURCE,
        ),
        AdapterEntry::catalog_flank(
            "ont-v14-rapid-barcode-flank",
            "ONT V14 rapid barcode adapter flank",
            DnaSequence::from_iupac_ascii(b"GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA"),
            SequencingPlatform::OxfordNanopore,
            &ONT_RAPID_V14_KIT,
            &ONT_CHEMISTRY_SOURCE,
        ),
    ],
};

pub(crate) static ONT_CDNA_V14: AdapterCatalog = AdapterCatalog {
    id: "ont-cdna-v14",
    name: "ONT cDNA/direct RNA adapters and primers V14",
    entries: &[
        AdapterEntry::five_prime_primer(
            "ont-v14-rt-primer",
            "ONT V14 RT primer",
            DnaSequence::from_iupac_ascii(b"CTTGCCTGTCGCTCTATCTTCAGAGGAG"),
            SequencingPlatform::OxfordNanopore,
            &ONT_CDNA_V14_KIT,
            &ONT_CHEMISTRY_SOURCE,
        ),
        AdapterEntry::internal_adapter(
            "ont-v14-cdna-rt-adapter",
            "ONT V14 cDNA RT adapter",
            DnaSequence::from_iupac_ascii(b"CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAG"),
            SequencingPlatform::OxfordNanopore,
            &ONT_CDNA_V14_KIT,
            &ONT_CHEMISTRY_SOURCE,
        ),
        AdapterEntry::five_prime_primer(
            "ont-v14-cdna-primer-forward",
            "ONT V14 cDNA primer forward",
            DnaSequence::from_iupac_ascii(
                b"ATCGCCTACCGTGACAAGAAAGTTGTCGGTGTCTTTGTGACTTGCCTGTCGCTCTATCTTC",
            ),
            SequencingPlatform::OxfordNanopore,
            &ONT_CDNA_V14_KIT,
            &ONT_CHEMISTRY_SOURCE,
        ),
        AdapterEntry::five_prime_primer(
            "ont-v14-cdna-primer-reverse",
            "ONT V14 cDNA primer reverse",
            DnaSequence::from_iupac_ascii(
                b"ATCGCCTACCGTGACAAGAAAGTTGTCGGTGTCTTTGTGTTTCTGTTGGTGCTGATATTGC",
            ),
            SequencingPlatform::OxfordNanopore,
            &ONT_CDNA_V14_KIT,
            &ONT_CHEMISTRY_SOURCE,
        ),
    ],
};

pub(crate) static ONT_NATIVE_BARCODING_V14: AdapterCatalog = AdapterCatalog {
    id: "ont-native-barcoding-v14",
    name: "ONT native barcoding V14",
    entries: &[
        native_flank(
            "ont-v14-native-forward-flank",
            "ONT V14 native barcode forward flank",
            "AAGGTTAA",
        ),
        native_flank(
            "ont-v14-native-forward-flank-3p",
            "ONT V14 native barcode forward 3' flank",
            "CAGCACCT",
        ),
        native_flank(
            "ont-v14-native-reverse-flank",
            "ONT V14 native barcode reverse flank",
            "GGTGCTG",
        ),
        native_flank(
            "ont-v14-native-reverse-flank-3p",
            "ONT V14 native barcode reverse 3' flank",
            "TTAACCTTAGCAAT",
        ),
        native_barcode(
            "ont-v14-native-barcode-01",
            "ONT Native Barcode 01",
            "AAGAAAGTTGTCGGTGTCTTTGTG",
        ),
        native_barcode(
            "ont-v14-native-barcode-02",
            "ONT Native Barcode 02",
            "TCGATTCCGTTTGTAGTCGTCTGT",
        ),
        native_barcode(
            "ont-v14-native-barcode-03",
            "ONT Native Barcode 03",
            "GAGTCTTGTGTCCCAGTTACCAGG",
        ),
        native_barcode(
            "ont-v14-native-barcode-04",
            "ONT Native Barcode 04",
            "TTCGGATTCTATCGTGTTTCCCTA",
        ),
        native_barcode(
            "ont-v14-native-barcode-05",
            "ONT Native Barcode 05",
            "CTTGTCCAGGGTTTGTGTAACCTT",
        ),
        native_barcode(
            "ont-v14-native-barcode-06",
            "ONT Native Barcode 06",
            "TTCTCGCAAAGGCAGAAAGTAGTC",
        ),
        native_barcode(
            "ont-v14-native-barcode-07",
            "ONT Native Barcode 07",
            "GTGTTACCGTGGGAATGAATCCTT",
        ),
        native_barcode(
            "ont-v14-native-barcode-08",
            "ONT Native Barcode 08",
            "TTCAGGGAACAAACCAAGTTACGT",
        ),
        native_barcode(
            "ont-v14-native-barcode-09",
            "ONT Native Barcode 09",
            "AACTAGGCACAGCGAGTCTTGGTT",
        ),
        native_barcode(
            "ont-v14-native-barcode-10",
            "ONT Native Barcode 10",
            "AAGCGTTGAAACCTTTGTCCTCTC",
        ),
        native_barcode(
            "ont-v14-native-barcode-11",
            "ONT Native Barcode 11",
            "GTTTCATCTATCGGAGGGAATGGA",
        ),
        native_barcode(
            "ont-v14-native-barcode-12",
            "ONT Native Barcode 12",
            "CAGGTAGAAAGAAGCAGAATCGGA",
        ),
        native_barcode(
            "ont-v14-native-barcode-13",
            "ONT Native Barcode 13",
            "AGAACGACTTCCATACTCGTGTGA",
        ),
        native_barcode(
            "ont-v14-native-barcode-14",
            "ONT Native Barcode 14",
            "AACGAGTCTCTTGGGACCCATAGA",
        ),
        native_barcode(
            "ont-v14-native-barcode-15",
            "ONT Native Barcode 15",
            "AGGTCTACCTCGCTAACACCACTG",
        ),
        native_barcode(
            "ont-v14-native-barcode-16",
            "ONT Native Barcode 16",
            "CGTCAACTGACAGTGGTTCGTACT",
        ),
        native_barcode(
            "ont-v14-native-barcode-17",
            "ONT Native Barcode 17",
            "ACCCTCCAGGAAAGTACCTCTGAT",
        ),
        native_barcode(
            "ont-v14-native-barcode-18",
            "ONT Native Barcode 18",
            "CCAAACCCAACAACCTAGATAGGC",
        ),
        native_barcode(
            "ont-v14-native-barcode-19",
            "ONT Native Barcode 19",
            "GTTCCTCGTGCAGTGTCAAGAGAT",
        ),
        native_barcode(
            "ont-v14-native-barcode-20",
            "ONT Native Barcode 20",
            "TTGCGTCCTGTTACGAGAACTCAT",
        ),
        native_barcode(
            "ont-v14-native-barcode-21",
            "ONT Native Barcode 21",
            "GAGCCTCTCATTGTCCGTTCTCTA",
        ),
        native_barcode(
            "ont-v14-native-barcode-22",
            "ONT Native Barcode 22",
            "ACCACTGCCATGTATCAAAGTACG",
        ),
        native_barcode(
            "ont-v14-native-barcode-23",
            "ONT Native Barcode 23",
            "CTTACTACCCAGTGAACCTCCTCG",
        ),
        native_barcode(
            "ont-v14-native-barcode-24",
            "ONT Native Barcode 24",
            "GCATAGTTCTGCATGATGGGTTAG",
        ),
    ],
};

const fn native_barcode(
    id: &'static str,
    name: &'static str,
    sequence: &'static str,
) -> AdapterEntry {
    AdapterEntry::catalog_barcode(
        id,
        name,
        DnaSequence::from_iupac_ascii(sequence.as_bytes()),
        SequencingPlatform::OxfordNanopore,
        &ONT_NATIVE_BARCODING_V14_KIT,
        &ONT_CHEMISTRY_SOURCE,
    )
}

const fn native_flank(
    id: &'static str,
    name: &'static str,
    sequence: &'static str,
) -> AdapterEntry {
    AdapterEntry::catalog_flank(
        id,
        name,
        DnaSequence::from_iupac_ascii(sequence.as_bytes()),
        SequencingPlatform::OxfordNanopore,
        &ONT_NATIVE_BARCODING_V14_KIT,
        &ONT_CHEMISTRY_SOURCE,
    )
}
