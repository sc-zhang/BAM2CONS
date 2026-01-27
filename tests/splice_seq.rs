// tests/splice_seq.rs

use bam2cons::cli::{LowConfFallback, Mode};
use bam2cons::consensus::{ConsensusParams, PosEvidence, WindowEvidence};
use bam2cons::splice::build_spliced_cds_seq;
use bam2cons::targets::{CdsBlock, SplicedTarget};

fn make_params_ref_mode() -> ConsensusParams {
    ConsensusParams {
        mode: Mode::Ref,
        lowconf_fallback: LowConfFallback::Ref,
        min_mapq: 0,
        min_baseq: 0,
        min_depth: 0,
        min_af: 0.0,
        min_iupac_af: 0.0,
        ignore_duplicates: false,

        enable_indel: false,
        max_indel_len: 50,
        min_indel_depth: 0,
        min_indel_af: 0.0,

        indel_end_dist: 0,
        indel_flank: 0,
        require_indel_strand_balance: false,

        indel_top2_ratio: 3.0,

        indel_mismatch_window: 0,
        indel_min_mismatch_bases: 0,
        indel_max_mismatch_rate: 1.0,
    }
}

fn make_window(ref_seq: &str) -> WindowEvidence {
    let ref_bytes = ref_seq.as_bytes().to_vec();
    let evidence = vec![PosEvidence::default(); ref_bytes.len()];
    WindowEvidence {
        chrom: "chr1".to_string(),
        start: 0,
        end: ref_bytes.len() as i64,
        ref_seq: ref_bytes,
        evidence,
    }
}

#[test]
fn splice_plus_orders_blocks_and_trims_phase() {
    // ref: 0..16
    let win = make_window("ACGTACGTACGTACGT");
    let params = make_params_ref_mode();

    // blocks intentionally unsorted; strand '+'
    // transcript-order should be start ascending => (0..4 phase=1) then (4..8 phase=0)
    // first block seq = "ACGT" trimmed phase 1 => "CGT"
    // second block seq = "ACGT"
    // result "CGTACGT"
    let st = SplicedTarget {
        id: 0,
        chrom: "chr1".to_string(),
        strand: Some('+'),
        gene_id: "g1".to_string(),
        transcript_id: "t1".to_string(),
        name: "g1|t1".to_string(),
        blocks: vec![
            CdsBlock {
                start: 4,
                end: 8,
                phase: 0,
            },
            CdsBlock {
                start: 0,
                end: 4,
                phase: 1,
            },
        ],
    };

    let (seq, stats) = build_spliced_cds_seq(&win, &st, &params, true).unwrap();
    assert_eq!(std::str::from_utf8(&seq).unwrap(), "CGTACGT");
    assert_eq!(stats.ref_len, 8);
    assert_eq!(stats.out_len, 7);
}

#[test]
fn splice_minus_respect_strand_changes_output_and_trims_phase() {
    // Use a non-palindromic reference to see revcomp difference clearly
    // ref: "AAAACCCCGGGGTTTT"
    // blocks: (0..4 "AAAA"), (8..12 "GGGG"), strand '-'
    // transcript order (descending start): 8..12 first (phase=2), then 0..4
    //
    // respect_strand=true:
    //   first block "GGGG" -> revcomp "CCCC" -> trim2 => "CC"
    //   second block "AAAA" -> revcomp "TTTT"
    //   => "CCTTTT"
    //
    // respect_strand=false:
    //   first block "GGGG" trim2 => "GG"
    //   second block "AAAA"
    //   => "GGAAAA"
    let win = make_window("AAAACCCCGGGGTTTT");
    let params = make_params_ref_mode();

    let st = SplicedTarget {
        id: 1,
        chrom: "chr1".to_string(),
        strand: Some('-'),
        gene_id: "g2".to_string(),
        transcript_id: "t9".to_string(),
        name: "g2|t9".to_string(),
        blocks: vec![
            CdsBlock {
                start: 0,
                end: 4,
                phase: 0,
            },
            CdsBlock {
                start: 8,
                end: 12,
                phase: 2,
            },
        ],
    };

    let (seq_rc, stats_rc) = build_spliced_cds_seq(&win, &st, &params, true).unwrap();
    assert_eq!(std::str::from_utf8(&seq_rc).unwrap(), "CCTTTT");
    assert_eq!(stats_rc.ref_len, 8);
    assert_eq!(stats_rc.out_len, 6);

    let (seq_no, stats_no) = build_spliced_cds_seq(&win, &st, &params, false).unwrap();
    assert_eq!(std::str::from_utf8(&seq_no).unwrap(), "GGAAAA");
    assert_eq!(stats_no.ref_len, 8);
    assert_eq!(stats_no.out_len, 6);
}
