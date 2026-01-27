use clap::{Parser, ValueEnum};
use std::path::PathBuf;

#[derive(Copy, Clone, Debug, ValueEnum)]
pub enum Mode {
    Ref,
    Majority,
    Strict,
    Iupac,
}

#[derive(Copy, Clone, Debug, ValueEnum)]
pub enum LowConfFallback {
    Ref,
    N,
}

#[derive(Parser, Debug)]
#[command(
    name = "pileup2fa",
    version,
    about = "Consensus FASTA from BAM/CRAM over BED/GFF/GTF targets (no VCF). Supports indels + contig batching + JSON stats."
)]
pub struct Cli {
    /// Reference FASTA (must have .fai)
    #[arg(long)]
    pub ref_fa: PathBuf,

    /// Directory containing BAM/CRAM files (each should have index .bai/.csi)
    #[arg(long)]
    pub bam_dir: PathBuf,

    /// Output directory (one FASTA + one JSON per sample)
    #[arg(long)]
    pub out_dir: PathBuf,

    /// BED targets file (0-based, half-open). Mutually exclusive with --gff/--gtf.
    #[arg(long, conflicts_with_all = ["gff", "gtf"])]
    pub bed: Option<PathBuf>,

    /// GFF3 targets file (1-based inclusive). Mutually exclusive with --bed/--gtf.
    #[arg(long, conflicts_with_all = ["bed", "gtf"])]
    pub gff: Option<PathBuf>,

    /// GTF targets file (1-based inclusive). Mutually exclusive with --bed/--gff.
    #[arg(long, conflicts_with_all = ["bed", "gff"])]
    pub gtf: Option<PathBuf>,

    /// For GFF/GTF: which feature type to extract as targets (e.g. gene, CDS, exon)
    #[arg(long, default_value = "gene")]
    pub feature: String,

    /// For GFF/GTF + feature=CDS: keep only the longest transcript per gene
    #[arg(long, default_value_t = true)]
    pub keep_longest: bool,

    /// When --keep-longest=false, append transcript id to FASTA record name
    #[arg(long, default_value_t = true)]
    pub append_transcript: bool,

    /// Consensus mode
    #[arg(long, value_enum, default_value = "majority")]
    pub mode: Mode,

    /// What to output at low-confidence base positions (depth/AF insufficient)
    /// ref: output reference base (more complete); N: output N (more conservative)
    #[arg(long, value_enum, default_value = "ref")]
    pub lowconf_fallback: LowConfFallback,

    /// Minimum mapping quality
    #[arg(long, default_value_t = 20)]
    pub min_mapq: u8,

    /// Minimum base quality
    #[arg(long, default_value_t = 20)]
    pub min_baseq: u8,

    /// Minimum depth to emit a confident base
    #[arg(long, default_value_t = 10)]
    pub min_depth: u32,

    /// Minimum allele fraction for strict mode (and for choosing single base confidently)
    #[arg(long, default_value_t = 0.7)]
    pub min_af: f32,

    /// In IUPAC mode: minimum allele fraction to consider a base present (for 2-allele mixture)
    #[arg(long, default_value_t = 0.2)]
    pub min_iupac_af: f32,

    /// Ignore PCR duplicates
    #[arg(long, default_value_t = true)]
    pub ignore_duplicates: bool,

    /// Enable applying small indels from pileup (changes sequence length)
    #[arg(long, default_value_t = true)]
    pub enable_indel: bool,

    /// Maximum indel length to apply (bp)
    #[arg(long, default_value_t = 50)]
    pub max_indel_len: usize,

    /// Minimum supporting reads for indel to be applied
    #[arg(long, default_value_t = 6)]
    pub min_indel_depth: u32,

    /// Minimum allele fraction for indel to be applied
    #[arg(long, default_value_t = 0.7)]
    pub min_indel_af: f32,

    /// Reject indels that occur within this distance (bp) to read ends
    #[arg(long, default_value_t = 10)]
    pub indel_end_dist: usize,

    /// Indel flank size (bp) for base-quality check around indel anchor
    #[arg(long, default_value_t = 5)]
    pub indel_flank: usize,

    /// Require indel support from both strands (more reliable, may miss true indel in biased data)
    #[arg(long, default_value_t = false)]
    pub require_indel_strand_balance: bool,

    /// Indel "top1 vs top2" support ratio. top1 must be >= top2 * ratio (top2=0 passes).
    #[arg(long, default_value_t = 3.0)]
    pub indel_top2_ratio: f32,

    /// Indel-support read mismatch check window (bp) around indel anchor (heuristic, no realignment)
    #[arg(long, default_value_t = 10)]
    pub indel_mismatch_window: usize,

    /// Minimum number of comparable bases in mismatch window; otherwise reject that read's indel support
    #[arg(long, default_value_t = 6)]
    pub indel_min_mismatch_bases: usize,

    /// Maximum mismatch rate allowed for indel-support read in mismatch window
    #[arg(long, default_value_t = 0.20)]
    pub indel_max_mismatch_rate: f32,

    /// If set, output reverse-complement for targets on '-' strand (BED col 6 or GFF/GTF strand)
    #[arg(long, default_value_t = false)]
    pub respect_strand: bool,

    /// Wrap FASTA sequence at N columns per line
    #[arg(long, default_value_t = 60)]
    pub wrap: usize,

    /// Number of samples processed in parallel
    #[arg(long, default_value_t = 4)]
    pub jobs: usize,

    /// Include CRAM files in bam_dir
    #[arg(long, default_value_t = true)]
    pub include_cram: bool,

    /// Include BAM files in bam_dir
    #[arg(long, default_value_t = true)]
    pub include_bam: bool,

    /// Skip samples if output already exists
    #[arg(long, default_value_t = false)]
    pub skip_existing: bool,

    /// Merge nearby targets into larger fetch windows (per contig) to reduce random IO
    #[arg(long, default_value_t = true)]
    pub batch_by_contig: bool,

    /// When batching: max gap between neighboring targets to merge into one window
    #[arg(long, default_value_t = 10_000)]
    pub merge_gap: i64,

    /// When batching: max window size (bp)
    #[arg(long, default_value_t = 5_000_000)]
    pub max_window: i64,

    /// Write per-sample JSON stats
    #[arg(long, default_value_t = true)]
    pub write_json: bool,
}
