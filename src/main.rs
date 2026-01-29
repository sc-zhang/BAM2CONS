// src/main.rs

mod cli;
mod consensus;
mod fasta_writer;
mod splice;
mod targets;
mod utils;

use anyhow::{Context, Result, bail};
use clap::Parser;
use cli::Cli;
use consensus::{
    ConsensusParams, SampleStats, collect_window_evidence, consensus_from_window_with_stats,
};
use fasta_writer::FastaWriter;
use rayon::prelude::*;
use rust_htslib::faidx;
use serde_json::json;
use std::collections::BTreeMap;
use std::ffi::OsStr;
use std::fs;
use std::path::{Path, PathBuf};
use tracing::{info, warn};
use walkdir::WalkDir;

use targets::{Target, read_bed, read_gff_like};
// Spliced CDS mode types / parser (you need these in src/targets.rs)
use splice::build_spliced_cds_seq;
use targets::{SplicedTarget, read_spliced_cds};
use utils::revcomp_in_place;

#[derive(Clone, Debug)]
struct Window {
    chrom: String,
    start: i64,
    end: i64,
    target_ids: Vec<usize>, // ids in either Target or SplicedTarget vector
}

fn main() -> Result<()> {
    tracing_subscriber::fmt().with_target(false).init();
    let cli = Cli::parse();

    fs::create_dir_all(&cli.out_dir)
        .with_context(|| format!("create out_dir: {}", cli.out_dir.display()))?;

    rayon::ThreadPoolBuilder::new()
        .num_threads(cli.jobs.max(1))
        .build_global()
        .ok();

    let bam_files = scan_bam_dir(&cli)?;
    if bam_files.is_empty() {
        bail!("No BAM/CRAM files found in {}", cli.bam_dir.display());
    }
    info!("Found {} alignment files", bam_files.len());

    let params = ConsensusParams {
        mode: cli.mode,
        lowconf_fallback: cli.lowconf_fallback,
        min_mapq: cli.min_mapq,
        min_baseq: cli.min_baseq,
        min_depth: cli.min_depth,
        min_af: cli.min_af,
        min_iupac_af: cli.min_iupac_af,
        ignore_duplicates: cli.ignore_duplicates,

        enable_indel: cli.enable_indel,
        max_indel_len: cli.max_indel_len,
        min_indel_depth: cli.min_indel_depth,
        min_indel_af: cli.min_indel_af,

        indel_end_dist: cli.indel_end_dist,
        indel_flank: cli.indel_flank,
        require_indel_strand_balance: cli.require_indel_strand_balance,

        indel_top2_ratio: cli.indel_top2_ratio,

        indel_mismatch_window: cli.indel_mismatch_window,
        indel_min_mismatch_bases: cli.indel_min_mismatch_bases,
        indel_max_mismatch_rate: cli.indel_max_mismatch_rate,
    };

    let params_json = json!({
        "mode": format!("{:?}", cli.mode),
        "lowconf_fallback": format!("{:?}", cli.lowconf_fallback),
        "min_mapq": cli.min_mapq,
        "min_baseq": cli.min_baseq,
        "min_depth": cli.min_depth,
        "min_af": cli.min_af,
        "min_iupac_af": cli.min_iupac_af,
        "ignore_duplicates": cli.ignore_duplicates,
        "enable_indel": cli.enable_indel,
        "max_indel_len": cli.max_indel_len,
        "min_indel_depth": cli.min_indel_depth,
        "min_indel_af": cli.min_indel_af,
        "indel_end_dist": cli.indel_end_dist,
        "indel_flank": cli.indel_flank,
        "require_indel_strand_balance": cli.require_indel_strand_balance,
        "indel_top2_ratio": cli.indel_top2_ratio,
        "indel_mismatch_window": cli.indel_mismatch_window,
        "indel_min_mismatch_bases": cli.indel_min_mismatch_bases,
        "indel_max_mismatch_rate": cli.indel_max_mismatch_rate,
        "batch_by_contig": cli.batch_by_contig,
        "merge_gap": cli.merge_gap,
        "max_window": cli.max_window,
        "respect_strand": cli.respect_strand,
        "write_json": cli.write_json,
        "keep_longest": cli.keep_longest,
        "append_transcript": cli.append_transcript
    });

    // CDS splicing mode: only for GFF/GTF with feature=CDS
    let cds_splice_mode = (cli.gff.is_some() || cli.gtf.is_some()) && cli.feature == "CDS";

    if cds_splice_mode {
        // -------- Spliced CDS path --------
        let (spliced, windows) = load_spliced_targets_and_windows(&cli)?;
        info!(
            "Loaded {} spliced CDS targets (keep_longest={}, append_transcript={})",
            spliced.len(),
            cli.keep_longest,
            cli.append_transcript
        );
        info!("Built {} fetch windows (spliced CDS)", windows.len());

        let spliced_by_id: BTreeMap<usize, &SplicedTarget> =
            spliced.iter().map(|t| (t.id, t)).collect();

        bam_files.par_iter().try_for_each(|bam_path| -> Result<()> {
            let sample = sample_name_from_path(bam_path);
            let out_fa = cli.out_dir.join(format!("{}.fa", sample));
            let out_json = cli.out_dir.join(format!("{}.json", sample));

            if cli.skip_existing && out_fa.exists() && (!cli.write_json || out_json.exists()) {
                info!("Skip existing: {}", sample);
                return Ok(());
            }

            let fa = faidx::Reader::from_path(&cli.ref_fa)
                .with_context(|| format!("open faidx: {}", cli.ref_fa.display()))?;

            info!("Processing sample: {} ({})", sample, bam_path.display());
            let mut fw = FastaWriter::create(&out_fa, cli.wrap)?;

            let mut sample_stats = SampleStats {
                sample: sample.clone(),
                bam_path: bam_path.display().to_string(),
                params: params_json.clone(),
                ..Default::default()
            };

            for w in windows.iter() {
                let win_ev = collect_window_evidence(&fa, bam_path, &w.chrom, w.start, w.end, &params)
                    .with_context(|| format!("collect evidence {} {}:{}-{}", sample, w.chrom, w.start, w.end))?;

                for tid in w.target_ids.iter() {
                    let st = *spliced_by_id.get(tid).expect("spliced target id exists");

                    let (seq, tstats) = build_spliced_cds_seq(&win_ev, st, &params, cli.respect_strand)
                        .with_context(|| format!("spliced CDS consensus {} {}", sample, st.name))?;

                    let header = format!(
                        "{}|{}|{}|mode={:?}|fallback={:?}|indel={}|min_depth={}|min_af={:.3}|min_indel_af={:.3}|top2_ratio={:.2}|mm_rate<= {:.2}",
                        sample,
                        st.name,
                        format_args!("{}:{}-{}", st.chrom, w.start, w.end),
                        cli.mode,
                        cli.lowconf_fallback,
                        cli.enable_indel,
                        cli.min_depth,
                        cli.min_af,
                        cli.min_indel_af,
                        cli.indel_top2_ratio,
                        cli.indel_max_mismatch_rate
                    );

                    fw.write_record(&header, &seq)?;

                    sample_stats.total_targets += 1;
                    sample_stats.total_ref_bp += tstats.ref_len;
                    sample_stats.total_out_bp += tstats.out_len;
                    sample_stats.total_n_bp += tstats.n_bases;
                    sample_stats.total_applied_ins += tstats.applied_insertions;
                    sample_stats.total_applied_del += tstats.applied_deletions;
                    sample_stats.targets.push(tstats);
                }
            }

            fw.flush()?;
            info!("Wrote {}", out_fa.display());

            if cli.write_json {
                let s = serde_json::to_string_pretty(&sample_stats)?;
                fs::write(&out_json, s)?;
                info!("Wrote {}", out_json.display());
            }

            Ok(())
        })?;
    } else {
        // -------- Simple interval path (BED or non-CDS feature) --------
        let mut targets = load_targets_intervals(&cli)?;
        targets.sort_by(|a, b| {
            (a.chrom.as_str(), a.start, a.end).cmp(&(b.chrom.as_str(), b.start, b.end))
        });
        info!("Loaded {} targets", targets.len());

        let windows = if cli.batch_by_contig {
            build_windows_intervals(&targets, cli.merge_gap, cli.max_window)
        } else {
            targets
                .iter()
                .map(|t| Window {
                    chrom: t.chrom.clone(),
                    start: t.start,
                    end: t.end,
                    target_ids: vec![t.id],
                })
                .collect()
        };
        info!(
            "Built {} fetch windows (batch_by_contig={})",
            windows.len(),
            cli.batch_by_contig
        );

        let target_by_id: BTreeMap<usize, &Target> = targets.iter().map(|t| (t.id, t)).collect();

        bam_files.par_iter().try_for_each(|bam_path| -> Result<()> {
            let sample = sample_name_from_path(bam_path);
            let out_fa = cli.out_dir.join(format!("{}.fa", sample));
            let out_json = cli.out_dir.join(format!("{}.json", sample));

            if cli.skip_existing && out_fa.exists() && (!cli.write_json || out_json.exists()) {
                info!("Skip existing: {}", sample);
                return Ok(());
            }

            let fa = faidx::Reader::from_path(&cli.ref_fa)
                .with_context(|| format!("open faidx: {}", cli.ref_fa.display()))?;

            info!("Processing sample: {} ({})", sample, bam_path.display());
            let mut fw = FastaWriter::create(&out_fa, cli.wrap)?;

            let mut sample_stats = SampleStats {
                sample: sample.clone(),
                bam_path: bam_path.display().to_string(),
                params: params_json.clone(),
                ..Default::default()
            };

            for w in windows.iter() {
                let win_ev = collect_window_evidence(&fa, bam_path, &w.chrom, w.start, w.end, &params)
                    .with_context(|| format!("collect evidence {} {}:{}-{}", sample, w.chrom, w.start, w.end))?;

                for tid in w.target_ids.iter() {
                    let t = *target_by_id.get(tid).expect("target id exists");
                    let (mut seq, tstats) = consensus_from_window_with_stats(&win_ev, t.start, t.end, &t.name, t.id, &params)
                        .with_context(|| format!("consensus {} {}:{}-{}", sample, t.chrom, t.start, t.end))?;

                    if cli.respect_strand && let Some('-') = t.strand {
                        revcomp_in_place(&mut seq);
                    }

                    let header = format!(
                        "{}|{}|{}:{}-{}|mode={:?}|fallback={:?}|indel={}|min_depth={}|min_af={:.3}|min_indel_af={:.3}|top2_ratio={:.2}|mm_rate<= {:.2}",
                        sample,
                        t.name,
                        t.chrom,
                        t.start,
                        t.end,
                        cli.mode,
                        cli.lowconf_fallback,
                        cli.enable_indel,
                        cli.min_depth,
                        cli.min_af,
                        cli.min_indel_af,
                        cli.indel_top2_ratio,
                        cli.indel_max_mismatch_rate
                    );

                    fw.write_record(&header, &seq)?;

                    sample_stats.total_targets += 1;
                    sample_stats.total_ref_bp += tstats.ref_len;
                    sample_stats.total_out_bp += tstats.out_len;
                    sample_stats.total_n_bp += tstats.n_bases;
                    sample_stats.total_applied_ins += tstats.applied_insertions;
                    sample_stats.total_applied_del += tstats.applied_deletions;
                    sample_stats.targets.push(tstats);
                }
            }

            fw.flush()?;
            info!("Wrote {}", out_fa.display());

            if cli.write_json {
                let s = serde_json::to_string_pretty(&sample_stats)?;
                fs::write(&out_json, s)?;
                info!("Wrote {}", out_json.display());
            }

            Ok(())
        })?;
    }

    Ok(())
}

fn load_targets_intervals(cli: &Cli) -> Result<Vec<Target>> {
    if let Some(bed) = &cli.bed {
        return read_bed(bed);
    }
    if let Some(gff) = &cli.gff {
        return read_gff_like(gff, &cli.feature, false);
    }
    if let Some(gtf) = &cli.gtf {
        return read_gff_like(gtf, &cli.feature, true);
    }
    bail!("You must provide one of --bed / --gff / --gtf");
}

fn load_spliced_targets_and_windows(cli: &Cli) -> Result<(Vec<SplicedTarget>, Vec<Window>)> {
    let is_gtf = cli.gtf.is_some();
    let path = if let Some(p) = &cli.gtf {
        p
    } else if let Some(p) = &cli.gff {
        p
    } else {
        bail!("CDS splicing mode requires --gff or --gtf");
    };

    // Parse spliced CDS targets (you implement this in targets.rs as discussed)
    let mut spliced = read_spliced_cds(path, is_gtf, cli.keep_longest, cli.append_transcript)?;
    // stable order
    spliced.sort_by(|a, b| {
        (
            a.chrom.as_str(),
            a.gene_id.as_str(),
            a.transcript_id.as_str(),
        )
            .cmp(&(
                b.chrom.as_str(),
                b.gene_id.as_str(),
                b.transcript_id.as_str(),
            ))
    });

    let windows = if cli.batch_by_contig {
        build_windows_spliced_bbox(&spliced, cli.merge_gap, cli.max_window)
    } else {
        spliced
            .iter()
            .map(|st| {
                let (s, e) = spliced_bbox(st);
                Window {
                    chrom: st.chrom.clone(),
                    start: s,
                    end: e,
                    target_ids: vec![st.id],
                }
            })
            .collect()
    };

    Ok((spliced, windows))
}

/// For spliced targets, build windows using each transcript's bounding box (min start, max end),
/// then optionally merge nearby transcript boxes (same contig) into bigger windows.
/// This guarantees every transcript's blocks are contained within at least one window that processes it.
fn build_windows_spliced_bbox(
    spliced: &[SplicedTarget],
    merge_gap: i64,
    max_window: i64,
) -> Vec<Window> {
    #[derive(Clone)]
    struct Boxed {
        chrom: String,
        start: i64,
        end: i64,
        id: usize,
    }

    let mut boxes: Vec<Boxed> = spliced
        .iter()
        .map(|st| {
            let (s, e) = spliced_bbox(st);
            Boxed {
                chrom: st.chrom.clone(),
                start: s,
                end: e,
                id: st.id,
            }
        })
        .collect();

    boxes.sort_by(|a, b| {
        (a.chrom.as_str(), a.start, a.end).cmp(&(b.chrom.as_str(), b.start, b.end))
    });

    let mut out = Vec::new();
    let mut i = 0usize;

    while i < boxes.len() {
        let chrom = boxes[i].chrom.clone();
        let mut j = i;
        while j < boxes.len() && boxes[j].chrom == chrom {
            j += 1;
        }
        let slice = &boxes[i..j];

        let mut cur_start = slice[0].start;
        let mut cur_end = slice[0].end;
        let mut cur_ids = vec![slice[0].id];

        for b in slice.iter().skip(1) {
            let gap = b.start - cur_end;
            let new_end = cur_end.max(b.end);
            let would_len = new_end - cur_start;

            let mergeable = gap <= merge_gap && would_len <= max_window;
            if mergeable {
                cur_end = new_end;
                cur_ids.push(b.id);
            } else {
                out.push(Window {
                    chrom: chrom.clone(),
                    start: cur_start,
                    end: cur_end,
                    target_ids: cur_ids,
                });
                cur_start = b.start;
                cur_end = b.end;
                cur_ids = vec![b.id];
            }
        }

        out.push(Window {
            chrom: chrom.clone(),
            start: cur_start,
            end: cur_end,
            target_ids: cur_ids,
        });
        i = j;
    }

    out
}

fn spliced_bbox(st: &SplicedTarget) -> (i64, i64) {
    let mut min_s = i64::MAX;
    let mut max_e = i64::MIN;
    for b in &st.blocks {
        if b.start < min_s {
            min_s = b.start;
        }
        if b.end > max_e {
            max_e = b.end;
        }
    }
    if min_s == i64::MAX || max_e == i64::MIN {
        (0, 0)
    } else {
        (min_s, max_e)
    }
}

fn build_windows_intervals(targets: &[Target], merge_gap: i64, max_window: i64) -> Vec<Window> {
    let mut out = Vec::new();

    let mut i = 0usize;
    while i < targets.len() {
        let chrom = targets[i].chrom.clone();
        let mut j = i;
        while j < targets.len() && targets[j].chrom == chrom {
            j += 1;
        }
        let slice = &targets[i..j];

        let mut cur_start = slice[0].start;
        let mut cur_end = slice[0].end;
        let mut cur_ids = vec![slice[0].id];

        for t in slice.iter().skip(1) {
            let gap = t.start - cur_end;
            let new_end = cur_end.max(t.end);
            let would_len = new_end - cur_start;

            let mergeable = gap <= merge_gap && would_len <= max_window;
            if mergeable {
                cur_end = new_end;
                cur_ids.push(t.id);
            } else {
                out.push(Window {
                    chrom: chrom.clone(),
                    start: cur_start,
                    end: cur_end,
                    target_ids: cur_ids,
                });
                cur_start = t.start;
                cur_end = t.end;
                cur_ids = vec![t.id];
            }
        }

        out.push(Window {
            chrom: chrom.clone(),
            start: cur_start,
            end: cur_end,
            target_ids: cur_ids,
        });
        i = j;
    }

    out
}

fn scan_bam_dir(cli: &Cli) -> Result<Vec<PathBuf>> {
    let mut out = Vec::new();
    for entry in WalkDir::new(&cli.bam_dir)
        .follow_links(true)
        .into_iter()
        .filter_map(|e| e.ok())
    {
        if !entry.file_type().is_file() {
            continue;
        }
        let path = entry.path();
        let ext = path
            .extension()
            .and_then(OsStr::to_str)
            .unwrap_or("")
            .to_ascii_lowercase();

        let is_bam = ext == "bam";
        let is_cram = ext == "cram";

        if (is_bam && !cli.include_bam) || (is_cram && !cli.include_cram) {
            continue;
        }
        if !is_bam && !is_cram {
            continue;
        }

        // Best-effort index warnings
        if is_bam {
            let bai = path.with_extension("bam.bai");
            let bai2 = path.with_extension("bai");
            let csi = path.with_extension("csi");
            if !(bai.exists() || bai2.exists() || csi.exists()) {
                warn!(
                    "Index not found for BAM (expected .bai/.csi): {}",
                    path.display()
                );
            }
        } else {
            let csi = path.with_extension("csi");
            if !csi.exists() {
                warn!(
                    "Index not found for CRAM (expected .csi): {}",
                    path.display()
                );
            }
        }

        out.push(path.to_path_buf());
    }
    out.sort();
    Ok(out)
}

fn sample_name_from_path(p: &Path) -> String {
    let file = p.file_name().and_then(|s| s.to_str()).unwrap_or("sample");
    file.strip_suffix(".bam")
        .or_else(|| file.strip_suffix(".cram"))
        .unwrap_or(file)
        .to_string()
}
