// src/splice.rs

use anyhow::Result;

use crate::consensus::{
    consensus_from_window_with_stats, ConsensusParams, TargetStats, WindowEvidence,
};
use crate::targets::SplicedTarget;
use crate::utils::revcomp_in_place;

/// Build one transcript CDS sequence by concatenating block consensuses.
/// - Orders blocks by strand ( + ascending; - descending )
/// - If respect_strand && strand=='-': reverse-complements each block sequence before concatenation
/// - Applies phase trim (0/1/2) to the FIRST block in transcript order
pub fn build_spliced_cds_seq(
    win_ev: &WindowEvidence,
    st: &SplicedTarget,
    params: &ConsensusParams,
    respect_strand: bool,
) -> Result<(Vec<u8>, TargetStats)> {
    let mut blocks = st.blocks.clone();
    let strand = st.strand.unwrap_or('+');

    if strand == '-' {
        blocks.sort_by(|a, b| b.start.cmp(&a.start)); // descending
    } else {
        blocks.sort_by(|a, b| a.start.cmp(&b.start)); // ascending
    }

    let ref_len_sum: i64 = blocks.iter().map(|b| b.end - b.start).sum();

    let mut merged = TargetStats {
        target_id: st.id,
        name: st.name.clone(),
        chrom: st.chrom.clone(),
        start: blocks.first().map(|b| b.start).unwrap_or(0),
        end: blocks.last().map(|b| b.end).unwrap_or(0),
        ref_len: ref_len_sum,
        ..Default::default()
    };

    let mut out: Vec<u8> = Vec::with_capacity(ref_len_sum.max(0) as usize);

    for (bi, b) in blocks.iter().enumerate() {
        let (mut seg, seg_stats) =
            consensus_from_window_with_stats(win_ev, b.start, b.end, &st.name, st.id, params)?;

        // output transcript orientation (optional)
        if respect_strand && strand == '-' {
            revcomp_in_place(&mut seg);
        }

        // phase trim only at transcript 5' start (first block after strand ordering)
        if bi == 0 && b.phase > 0 {
            let k = b.phase as usize;
            if seg.len() > k {
                seg.drain(0..k);
            } else {
                seg.clear();
            }
        }

        // merge stats (sum per-block; out_len overwritten at end)
        merged.n_bases += seg_stats.n_bases;
        merged.ref_fallback_bases += seg_stats.ref_fallback_bases;
        merged.low_depth_bases += seg_stats.low_depth_bases;
        merged.low_af_bases += seg_stats.low_af_bases;

        merged.applied_insertions += seg_stats.applied_insertions;
        merged.applied_insertion_bp += seg_stats.applied_insertion_bp;
        merged.applied_deletions += seg_stats.applied_deletions;
        merged.applied_deletion_bp += seg_stats.applied_deletion_bp;

        merged.skipped_indel_lowconf += seg_stats.skipped_indel_lowconf;
        merged.skipped_indel_conflict += seg_stats.skipped_indel_conflict;
        merged.skipped_indel_topratio += seg_stats.skipped_indel_topratio;
        merged.skipped_indel_long += seg_stats.skipped_indel_long;

        merged.indel_rejected_end += seg_stats.indel_rejected_end;
        merged.indel_rejected_flankq += seg_stats.indel_rejected_flankq;
        merged.indel_rejected_nonacgt += seg_stats.indel_rejected_nonacgt;
        merged.indel_rejected_mismatch += seg_stats.indel_rejected_mismatch;
        merged.indel_rejected_mismatch_lowcomp += seg_stats.indel_rejected_mismatch_lowcomp;

        out.extend_from_slice(&seg);
    }

    merged.out_len = out.len() as i64;
    Ok((out, merged))
}
