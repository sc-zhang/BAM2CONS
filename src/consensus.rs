use crate::cli::{LowConfFallback, Mode};
use crate::utils::{base_to_index, index_to_base, is_acgt, iupac_code, passes_flank_qual};
use anyhow::{anyhow, Context, Result};
use rust_htslib::bam::{
    self,
    pileup::{Indel, Pileup},
    Read,
};
use rust_htslib::faidx;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::Path;

#[derive(Clone, Debug)]
pub struct ConsensusParams {
    pub mode: Mode,
    pub lowconf_fallback: LowConfFallback,
    pub min_mapq: u8,
    pub min_baseq: u8,
    pub min_depth: u32,
    pub min_af: f32,
    pub min_iupac_af: f32,
    pub ignore_duplicates: bool,

    pub enable_indel: bool,
    pub max_indel_len: usize,
    pub min_indel_depth: u32,
    pub min_indel_af: f32,

    pub indel_end_dist: usize,
    pub indel_flank: usize,
    pub require_indel_strand_balance: bool,

    pub indel_top2_ratio: f32,

    pub indel_mismatch_window: usize,
    pub indel_min_mismatch_bases: usize,
    pub indel_max_mismatch_rate: f32,
}

#[derive(Clone, Debug, Default)]
pub struct PosEvidence {
    pub base_counts: [u32; 4],
    pub base_depth: u32,

    pub ins: HashMap<Vec<u8>, u32>,
    pub ins_plus: HashMap<Vec<u8>, u32>,
    pub ins_minus: HashMap<Vec<u8>, u32>,
    pub ins_depth: u32,

    pub del: HashMap<usize, u32>,
    pub del_plus: HashMap<usize, u32>,
    pub del_minus: HashMap<usize, u32>,
    pub del_depth: u32,

    // diagnostics for JSON
    pub indel_rejected_end: u32,
    pub indel_rejected_flankq: u32,
    pub indel_rejected_long: u32,
    pub indel_rejected_nonacgt: u32,
    pub indel_rejected_mismatch: u32,
    pub indel_rejected_mismatch_lowcomp: u32,
}

#[derive(Clone, Debug)]
pub struct WindowEvidence {
    pub chrom: String,
    pub start: i64,
    pub end: i64,
    pub ref_seq: Vec<u8>,
    pub evidence: Vec<PosEvidence>,
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct TargetStats {
    pub target_id: usize,
    pub name: String,
    pub chrom: String,
    pub start: i64,
    pub end: i64,
    pub ref_len: i64,

    pub out_len: i64,
    pub n_bases: i64,
    pub ref_fallback_bases: i64,

    pub low_depth_bases: i64,
    pub low_af_bases: i64,

    pub applied_insertions: i64,
    pub applied_insertion_bp: i64,
    pub applied_deletions: i64,
    pub applied_deletion_bp: i64,

    pub skipped_indel_lowconf: i64,
    pub skipped_indel_conflict: i64,
    pub skipped_indel_topratio: i64,

    pub skipped_indel_long: i64,

    pub indel_rejected_end: i64,
    pub indel_rejected_flankq: i64,
    pub indel_rejected_nonacgt: i64,
    pub indel_rejected_mismatch: i64,
    pub indel_rejected_mismatch_lowcomp: i64,
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct SampleStats {
    pub sample: String,
    pub bam_path: String,
    pub params: serde_json::Value,

    pub targets: Vec<TargetStats>,

    pub total_targets: i64,
    pub total_ref_bp: i64,
    pub total_out_bp: i64,
    pub total_n_bp: i64,

    pub total_applied_ins: i64,
    pub total_applied_del: i64,
}

pub fn collect_window_evidence(
    fa: &faidx::Reader,
    bam_path: &Path,
    chrom: &str,
    start: i64,
    end: i64,
    params: &ConsensusParams,
) -> Result<WindowEvidence> {
    if end <= start {
        return Err(anyhow!("invalid window {}:{}-{}", chrom, start, end));
    }

    let ref_seq = fa
        .fetch_seq_string(chrom, start as usize, (end - 1) as usize)
        .with_context(|| format!("faidx fetch {}:{}-{}", chrom, start, end))?
        .into_bytes();

    let mut evidence = vec![PosEvidence::default(); ref_seq.len()];

    if matches!(params.mode, Mode::Ref) {
        return Ok(WindowEvidence {
            chrom: chrom.to_string(),
            start,
            end,
            ref_seq,
            evidence,
        });
    }

    let mut reader = bam::IndexedReader::from_path(bam_path)
        .with_context(|| format!("open indexed BAM/CRAM: {}", bam_path.display()))?;

    let header = reader.header().to_owned();
    let tid = header.tid(chrom.as_bytes()).ok_or_else(|| {
        anyhow!(
            "chrom '{}' not found in BAM header: {}",
            chrom,
            bam_path.display()
        )
    })?;

    reader.fetch((tid, start, end)).with_context(|| {
        format!(
            "bam fetch {}:{}-{} for {}",
            chrom,
            start,
            end,
            bam_path.display()
        )
    })?;

    for p in reader.pileup() {
        let pile: Pileup = p?;
        let pos = pile.pos() as i64;
        if pos < start || pos >= end {
            continue;
        }
        let idx = (pos - start) as usize;
        if idx >= evidence.len() {
            continue;
        }

        for aln in pile.alignments() {
            if aln.is_refskip() {
                continue;
            }
            let rec = aln.record();

            // ignore bad/secondary/supplementary/failed QC
            if rec.is_secondary() || rec.is_supplementary() || rec.is_quality_check_failed() {
                continue;
            }
            if params.ignore_duplicates && rec.is_duplicate() {
                continue;
            }
            if rec.mapq() < params.min_mapq {
                continue;
            }

            // base
            if !aln.is_del() {
                if let Some(qpos) = aln.qpos() {
                    let bq = rec.qual()[qpos] as u8;
                    if bq >= params.min_baseq {
                        let base = rec.seq().as_bytes()[qpos].to_ascii_uppercase();
                        if let Some(i) = base_to_index(base) {
                            evidence[idx].base_counts[i] += 1;
                            evidence[idx].base_depth += 1;
                        }
                    }
                }
            }

            // indel after this base
            if !params.enable_indel {
                continue;
            }
            let indel = aln.indel();
            if matches!(indel, Indel::None) {
                continue;
            }

            let qpos = match aln.qpos() {
                Some(q) => q,
                None => continue,
            };

            // read-end distance
            let read_len = rec.seq_len();
            if qpos < params.indel_end_dist
                || (read_len.saturating_sub(qpos + 1)) < params.indel_end_dist
            {
                evidence[idx].indel_rejected_end += 1;
                continue;
            }

            // flank quality around anchor
            if !passes_flank_qual(rec.qual(), qpos, params.indel_flank, params.min_baseq) {
                evidence[idx].indel_rejected_flankq += 1;
                continue;
            }

            // mismatch-rate heuristic around anchor
            match passes_mismatch_rate(
                &ref_seq,
                idx,
                &rec,
                qpos,
                params.indel_mismatch_window,
                params.indel_min_mismatch_bases,
                params.indel_max_mismatch_rate,
                params.min_baseq,
            ) {
                Ok(true) => {}
                Ok(false) => {
                    evidence[idx].indel_rejected_mismatch += 1;
                    continue;
                }
                Err(MismatchErr::TooFewComparable) => {
                    evidence[idx].indel_rejected_mismatch_lowcomp += 1;
                    continue;
                }
            }

            let is_rev = rec.is_reverse();

            match indel {
                Indel::None => {}
                Indel::Ins(ilen_i) => {
                    let ilen = ilen_i as usize;
                    if ilen > params.max_indel_len {
                        evidence[idx].indel_rejected_long += 1;
                        continue;
                    }
                    if let Some(ins_seq) = extract_insertion_bases(&aln, ilen) {
                        if ins_seq.iter().any(|&b| !is_acgt(b)) {
                            evidence[idx].indel_rejected_nonacgt += 1;
                            continue;
                        }
                        *evidence[idx].ins.entry(ins_seq.clone()).or_insert(0) += 1;
                        if !is_rev {
                            *evidence[idx].ins_plus.entry(ins_seq).or_insert(0) += 1;
                        } else {
                            *evidence[idx].ins_minus.entry(ins_seq).or_insert(0) += 1;
                        }
                        evidence[idx].ins_depth += 1;
                    }
                }
                Indel::Del(dlen_i) => {
                    let dlen = dlen_i as usize;
                    if dlen > params.max_indel_len {
                        evidence[idx].indel_rejected_long += 1;
                        continue;
                    }
                    *evidence[idx].del.entry(dlen).or_insert(0) += 1;
                    if !is_rev {
                        *evidence[idx].del_plus.entry(dlen).or_insert(0) += 1;
                    } else {
                        *evidence[idx].del_minus.entry(dlen).or_insert(0) += 1;
                    }
                    evidence[idx].del_depth += 1;
                }
            }
        }
    }

    Ok(WindowEvidence {
        chrom: chrom.to_string(),
        start,
        end,
        ref_seq,
        evidence,
    })
}

fn extract_insertion_bases(aln: &bam::pileup::Alignment, ilen: usize) -> Option<Vec<u8>> {
    let rec = aln.record();
    let qpos = aln.qpos()?;
    let seq = rec.seq().as_bytes();
    let start = qpos + 1;
    let end = start + ilen;
    if end <= seq.len() {
        let mut v = Vec::with_capacity(ilen);
        for &b in &seq[start..end] {
            v.push(b.to_ascii_uppercase());
        }
        Some(v)
    } else {
        None
    }
}

enum MismatchErr {
    TooFewComparable,
}

/// Heuristic mismatch filter around indel anchor (no CIGAR walk; conservative).
/// Compares read bases at qpos±d with ref bases at idx±d (within window),
/// counts only where baseQ>=min_baseq and ref is A/C/G/T.
/// - If comparable < min_compared: Err(TooFewComparable)
/// - Else returns mismatch_rate <= max_rate
fn passes_mismatch_rate(
    ref_seq: &[u8],
    idx: usize,
    rec: &bam::Record,
    qpos: usize,
    win: usize,
    min_compared: usize,
    max_rate: f32,
    min_baseq: u8,
) -> Result<bool, MismatchErr> {
    let read = rec.seq().as_bytes();
    let qual = rec.qual();

    let mut compared = 0usize;
    let mut mism = 0usize;

    for d in 1..=win {
        // left
        if idx >= d && qpos >= d {
            let rb = read[qpos - d].to_ascii_uppercase();
            let rq = qual[qpos - d] as u8;
            let fb = ref_seq[idx - d].to_ascii_uppercase();
            if rq >= min_baseq && is_acgt(fb) && is_acgt(rb) {
                compared += 1;
                if rb != fb {
                    mism += 1;
                }
            }
        }
        // right
        if idx + d < ref_seq.len() && qpos + d < read.len() {
            let rb = read[qpos + d].to_ascii_uppercase();
            let rq = qual[qpos + d] as u8;
            let fb = ref_seq[idx + d].to_ascii_uppercase();
            if rq >= min_baseq && is_acgt(fb) && is_acgt(rb) {
                compared += 1;
                if rb != fb {
                    mism += 1;
                }
            }
        }
    }

    if compared < min_compared {
        return Err(MismatchErr::TooFewComparable);
    }

    let rate = mism as f32 / compared as f32;
    Ok(rate <= max_rate)
}

pub fn consensus_from_window_with_stats(
    win: &WindowEvidence,
    target_start: i64,
    target_end: i64,
    name: &str,
    target_id: usize,
    params: &ConsensusParams,
) -> Result<(Vec<u8>, TargetStats)> {
    if target_end <= target_start {
        return Ok((Vec::new(), TargetStats::default()));
    }
    if target_start < win.start || target_end > win.end {
        return Err(anyhow!(
            "target {}-{} not inside window {}-{}",
            target_start,
            target_end,
            win.start,
            win.end
        ));
    }

    let s = (target_start - win.start) as usize;
    let e = (target_end - win.start) as usize;

    let mut stats = TargetStats {
        target_id,
        name: name.to_string(),
        chrom: win.chrom.clone(),
        start: target_start,
        end: target_end,
        ref_len: (target_end - target_start),
        ..Default::default()
    };

    let mut out: Vec<u8> = Vec::with_capacity(e - s);

    let mut i = s;
    while i < e {
        let ref_base = win.ref_seq[i].to_ascii_uppercase();
        let ev = &win.evidence[i];

        // diagnostics aggregate
        stats.indel_rejected_end += ev.indel_rejected_end as i64;
        stats.indel_rejected_flankq += ev.indel_rejected_flankq as i64;
        stats.indel_rejected_nonacgt += ev.indel_rejected_nonacgt as i64;
        stats.indel_rejected_mismatch += ev.indel_rejected_mismatch as i64;
        stats.indel_rejected_mismatch_lowcomp += ev.indel_rejected_mismatch_lowcomp as i64;
        stats.skipped_indel_long += ev.indel_rejected_long as i64;

        // base
        let (base_call, base_reason) =
            call_base_with_reason(ev.base_counts, ev.base_depth, ref_base, params);
        if base_call == b'N' {
            stats.n_bases += 1;
        }
        match base_reason {
            BaseReason::LowDepth => stats.low_depth_bases += 1,
            BaseReason::LowAf => stats.low_af_bases += 1,
            BaseReason::FallbackRef => stats.ref_fallback_bases += 1,
            BaseReason::Ok => {}
        }
        out.push(base_call);

        // indel after this base
        if params.enable_indel && !matches!(params.mode, Mode::Ref) {
            let ins_best = best_insertion_with_second(&ev.ins, &ev.ins_plus, &ev.ins_minus);
            let del_best = best_deletion_with_second(&ev.del, &ev.del_plus, &ev.del_minus);

            // Apply top1/top2 ratio filter (per-type)
            let ins_ratio_ok = ins_best.top2_count == 0
                || (ins_best.top1_count as f32)
                    >= (ins_best.top2_count as f32) * params.indel_top2_ratio;
            let del_ratio_ok = del_best.top2_count == 0
                || (del_best.top1_count as f32)
                    >= (del_best.top2_count as f32) * params.indel_top2_ratio;

            let ins_af = if ev.ins_depth > 0 {
                ins_best.top1_count as f32 / ev.ins_depth as f32
            } else {
                0.0
            };
            let del_af = if ev.del_depth > 0 {
                del_best.top1_count as f32 / ev.del_depth as f32
            } else {
                0.0
            };

            let ins_strand_ok = if params.require_indel_strand_balance && ins_best.top1_count > 0 {
                ins_best.top1_plus > 0 && ins_best.top1_minus > 0
            } else {
                true
            };
            let del_strand_ok = if params.require_indel_strand_balance && del_best.top1_count > 0 {
                del_best.top1_plus > 0 && del_best.top1_minus > 0
            } else {
                true
            };

            let ins_ok = ins_best.top1_count > 0
                && ev.ins_depth >= params.min_indel_depth
                && ins_af >= params.min_indel_af
                && ins_strand_ok
                && ins_ratio_ok;

            let del_ok = del_best.top1_len > 0
                && ev.del_depth >= params.min_indel_depth
                && del_af >= params.min_indel_af
                && del_strand_ok
                && del_ratio_ok;

            // Track ratio failures (when there is signal but ratio gate blocks)
            if ins_best.top1_count > 0 && !ins_ratio_ok {
                stats.skipped_indel_topratio += 1;
            }
            if del_best.top1_count > 0 && !del_ratio_ok {
                stats.skipped_indel_topratio += 1;
            }

            // If both pass, pick stronger; if neither pass do nothing
            if ins_ok || del_ok {
                let pick_insertion = match (ins_ok, del_ok) {
                    (true, false) => true,
                    (false, true) => false,
                    (true, true) => {
                        if (ins_af - del_af).abs() > 1e-6 {
                            ins_af > del_af
                        } else {
                            ins_best.top1_count >= del_best.top1_count
                        }
                    }
                    _ => false,
                };

                if pick_insertion {
                    if let Some(seq) = ins_best.top1_seq {
                        stats.applied_insertions += 1;
                        stats.applied_insertion_bp += seq.len() as i64;
                        out.extend_from_slice(&seq);
                    }
                } else {
                    // deletion
                    stats.applied_deletions += 1;
                    stats.applied_deletion_bp += del_best.top1_len as i64;
                    i = i.saturating_add(1 + del_best.top1_len);
                    continue;
                }
            } else {
                // if there is some signal but below AF/depth/strand/ratio thresholds
                if (!ev.ins.is_empty() || !ev.del.is_empty())
                    && (ins_best.top1_count > 0 || del_best.top1_count > 0)
                {
                    stats.skipped_indel_lowconf += 1;
                }
                // If both had strong-looking signal but conflict (rare here because we pick max if both pass),
                // you can treat as conflict; keep counter in case you extend logic later.
                if ins_best.top1_count > 0 && del_best.top1_count > 0 && !(ins_ok || del_ok) {
                    stats.skipped_indel_conflict += 1;
                }
            }
        }

        i += 1;
    }

    stats.out_len = out.len() as i64;
    Ok((out, stats))
}

#[derive(Copy, Clone, Debug)]
enum BaseReason {
    Ok,
    LowDepth,
    LowAf,
    FallbackRef,
}

fn call_base_with_reason(
    counts: [u32; 4],
    depth: u32,
    ref_base: u8,
    params: &ConsensusParams,
) -> (u8, BaseReason) {
    if matches!(params.mode, Mode::Ref) {
        return (ref_base, BaseReason::Ok);
    }

    if depth < params.min_depth {
        return match params.lowconf_fallback {
            LowConfFallback::Ref => (ref_base, BaseReason::FallbackRef),
            LowConfFallback::N => (b'N', BaseReason::LowDepth),
        };
    }

    let mut pairs: Vec<(u8, u32)> = (0..4).map(|i| (index_to_base(i), counts[i])).collect();
    pairs.sort_by(|a, b| b.1.cmp(&a.1));
    let (b1, c1) = pairs[0];
    let (b2, c2) = pairs[1];
    let f1 = c1 as f32 / depth as f32;
    let f2 = c2 as f32 / depth as f32;

    match params.mode {
        Mode::Majority => (b1, BaseReason::Ok),
        Mode::Strict => {
            if f1 >= params.min_af {
                (b1, BaseReason::Ok)
            } else {
                match params.lowconf_fallback {
                    LowConfFallback::Ref => (ref_base, BaseReason::FallbackRef),
                    LowConfFallback::N => (b'N', BaseReason::LowAf),
                }
            }
        }
        Mode::Iupac => {
            if f1 >= params.min_af {
                (b1, BaseReason::Ok)
            } else if f1 >= params.min_iupac_af && f2 >= params.min_iupac_af {
                (iupac_code(b1, b2).unwrap_or(b'N'), BaseReason::Ok)
            } else {
                match params.lowconf_fallback {
                    LowConfFallback::Ref => (ref_base, BaseReason::FallbackRef),
                    LowConfFallback::N => (b'N', BaseReason::LowAf),
                }
            }
        }
        Mode::Ref => (ref_base, BaseReason::Ok),
    }
}

struct BestIns {
    top1_seq: Option<Vec<u8>>,
    top1_count: u32,
    top2_count: u32,
    top1_plus: u32,
    top1_minus: u32,
}

fn best_insertion_with_second(
    all: &HashMap<Vec<u8>, u32>,
    plus: &HashMap<Vec<u8>, u32>,
    minus: &HashMap<Vec<u8>, u32>,
) -> BestIns {
    let mut v: Vec<(Vec<u8>, u32)> = all.iter().map(|(k, &c)| (k.clone(), c)).collect();
    v.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
    let top1 = v.get(0).cloned();
    let top2 = v.get(1).cloned();
    if let Some((seq, cnt)) = top1 {
        let p = *plus.get(&seq).unwrap_or(&0);
        let m = *minus.get(&seq).unwrap_or(&0);
        BestIns {
            top1_seq: Some(seq),
            top1_count: cnt,
            top2_count: top2.map(|x| x.1).unwrap_or(0),
            top1_plus: p,
            top1_minus: m,
        }
    } else {
        BestIns {
            top1_seq: None,
            top1_count: 0,
            top2_count: 0,
            top1_plus: 0,
            top1_minus: 0,
        }
    }
}

struct BestDel {
    top1_len: usize,
    top1_count: u32,
    top2_count: u32,
    top1_plus: u32,
    top1_minus: u32,
}

fn best_deletion_with_second(
    all: &HashMap<usize, u32>,
    plus: &HashMap<usize, u32>,
    minus: &HashMap<usize, u32>,
) -> BestDel {
    let mut v: Vec<(usize, u32)> = all.iter().map(|(&k, &c)| (k, c)).collect();
    v.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
    let top1 = v.get(0).cloned();
    let top2 = v.get(1).cloned();
    if let Some((len, cnt)) = top1 {
        let p = *plus.get(&len).unwrap_or(&0);
        let m = *minus.get(&len).unwrap_or(&0);
        BestDel {
            top1_len: len,
            top1_count: cnt,
            top2_count: top2.map(|x| x.1).unwrap_or(0),
            top1_plus: p,
            top1_minus: m,
        }
    } else {
        BestDel {
            top1_len: 0,
            top1_count: 0,
            top2_count: 0,
            top1_plus: 0,
            top1_minus: 0,
        }
    }
}
