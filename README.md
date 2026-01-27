# bam2cons

Generate **per-sample consensus FASTA** for **target regions** directly from **BAM/CRAM pileup** — **without VCF/SNP calling**.  
Supports **small indels**, **contig-batched processing**, **multi-sample parallelism**, **JSON QC/stats**, and **spliced CDS reconstruction** from GFF/GTF.

> Design goal: apply **indels conservatively** (avoid false length changes), while keeping sequences **as complete as possible** (low-confidence bases fall back to REF by default).

---

## Features

- **No VCF / no SNP calling**: pileup evidence → consensus.
- Inputs:
  - Reference **FASTA** (requires `.fai`)
  - Targets via **BED** (0-based half-open) or **GFF3/GTF**
  - A directory of **BAM/CRAM** files (indexed)
- Output:
  - One `sample.fa` per BAM/CRAM
  - One `sample.json` QC/stats per sample (optional)
- **Indel support** (small ins/del, changes sequence length) with anti-false-positive heuristics:
  - read-end distance filter
  - flank base-quality filter
  - local mismatch-rate filter around indel anchor (heuristic, no realignment)
  - top1/top2 support ratio filter (reject ambiguous indel shapes)
  - optional strand-balance requirement
- **Batch-by-contig windows**: merge nearby targets and scan pileup once per window (fewer random fetches).
- **Parallel**: process multiple samples concurrently.
- **Spliced CDS mode (GFF/GTF + CDS)**:
  - **One transcript → one spliced CDS** (concatenate CDS blocks)
  - Default: **keep only the longest transcript per gene**
  - Optional: output all transcripts per gene (append transcript ID in record name)
  - If **no transcript/mRNA info exists**, treat subsequent CDS as belonging to the same gene (ordered-file assumption)

---

## Install / Build

```bash
cargo build --release
./target/release/bam2cons --help
```

---

## Quick Start

### BED targets + BAM directory

```bash
./target/release/bam2cons \
  --ref-fa ref.fa \
  --bed targets.bed \
  --bam-dir ./bams \
  --out-dir ./out \
  --mode strict \
  --lowconf-fallback ref \
  --min-mapq 30 --min-baseq 25 \
  --min-depth 8 --min-af 0.75 \
  --enable-indel \
  --min-indel-depth 8 --min-indel-af 0.85 \
  --indel-top2-ratio 3.0 \
  --indel-mismatch-window 10 --indel-min-mismatch-bases 6 --indel-max-mismatch-rate 0.20 \
  --batch-by-contig --merge-gap 10000 --max-window 5000000 \
  --write-json \
  --jobs 8
```

---

## Spliced CDS Reconstruction (GFF/GTF)

### Default: keep the longest transcript per gene

```bash
./target/release/bam2cons \
  --ref-fa ref.fa \
  --gtf annot.gtf \
  --feature CDS \
  --bam-dir ./bams \
  --out-dir ./out \
  --keep-longest \
  --write-json
```

### Output all transcripts (append transcript ID)

```bash
./target/release/bam2cons \
  --ref-fa ref.fa \
  --gff annot.gff3 \
  --feature CDS \
  --bam-dir ./bams \
  --out-dir ./out \
  --keep-longest false \
  --append-transcript true \
  --write-json
```

---

## File Formats

### Reference FASTA
Must have `.fai`:

```bash
samtools faidx ref.fa
```

### Targets

#### BED (BED3 / BED6)
- **0-based**, half-open intervals: `[start, end)`
- `--bed` selects BED as target source

Supported columns:

1. **chrom**: contig/chromosome name (must match reference FASTA)
2. **start**: 0-based start (inclusive)
3. **end**: 0-based end (exclusive)

Optional (BED6):

4. **name**: target name  
   - Used in FASTA record naming and JSON stats.
   - If absent, bam2cons auto-generates a name.

5. **score**: numeric score (often 0–1000 by convention)  
   - **Ignored by bam2cons** (no filtering/scoring).

6. **strand**: `+` or `-`  
   - Only used when `--respect-strand` is enabled.
   - If `strand == '-'`, output is **reverse-complemented** for that target.

Example (BED6):

```text
chr1    1000    2000    geneA   0   +
chr2    500     900     geneB   0   -
```

#### GFF3 / GTF
- **1-based inclusive** coordinates
- `--feature` selects which feature becomes targets (e.g., `gene`, `exon`, `CDS`)
- **CDS splicing mode** is enabled when `--feature CDS` is used with `--gff/--gtf`

---

## Output

### FASTA
- `--wrap N` controls FASTA line wrapping (format only; does not affect consensus).
- `--respect-strand` reverse-complements `-` strand outputs (BED col6 / GFF/GTF strand).

### JSON QC / Stats
If enabled, each `<sample>.json` contains:
- run parameters
- per-target statistics
- totals

Useful per-target counters:
- `ref_len`, `out_len`
- `n_bases`
- `ref_fallback_bases`
- `low_depth_bases`, `low_af_bases`
- `applied_insertions`, `applied_insertion_bp`
- `applied_deletions`, `applied_deletion_bp`
- indel diagnostics:
  - `indel_rejected_end`
  - `indel_rejected_flankq`
  - `indel_rejected_mismatch`
  - `indel_rejected_mismatch_lowcomp`
  - `indel_rejected_nonacgt`
  - `skipped_indel_topratio`
  - `skipped_indel_lowconf`

---

# Command Line Parameters (Full Reference)

> Below documents **all parameters** as implemented in `Cli` (`src/cli.rs`).

---

## Input / Output

### `--ref-fa <PATH>`
Reference FASTA (requires `.fai`). Used for:
- fetching reference bases
- generating consensus in `ref` mode
- fallback when low confidence (if `--lowconf-fallback ref`)

### `--bam-dir <DIR>`
Directory containing BAM/CRAM files. The tool scans recursively.

### `--out-dir <DIR>`
Output directory. Writes:
- `<sample>.fa`
- `<sample>.json` (if enabled)

### `--bed <PATH>`
BED targets (0-based, half-open). Mutually exclusive with `--gff/--gtf`.

### `--gff <PATH>`
GFF3 targets. Mutually exclusive with `--bed/--gtf`.

### `--gtf <PATH>`
GTF targets. Mutually exclusive with `--bed/--gff`.

### `--feature <STRING>` (default: `gene`)
For GFF/GTF: which feature type to extract as targets (e.g. `gene`, `exon`, `CDS`).

- If `--feature CDS` with `--gff/--gtf`, bam2cons enters **spliced CDS mode** (see below).

---

## Consensus Calling Mode

### `--mode <ref|majority|strict|iupac>` (default: `majority`)
Controls how the per-position consensus base is chosen.

- `ref`: output reference bases only (fastest; ignores BAM evidence).
- `majority`: choose the most frequent base among A/C/G/T (after filters).
- `strict`: choose top base only if its allele fraction ≥ `--min-af`; otherwise fallback.
- `iupac`: if top base AF < `--min-af` but both top1 & top2 AF ≥ `--min-iupac-af`, output IUPAC code; otherwise fallback.

### `--lowconf-fallback <ref|n>` (default: `ref`)
When confidence is insufficient (low depth / low AF), output:
- `ref`: reference base (more complete; fewer Ns)
- `n`: output `N` (more conservative)

---

## Read Filtering

### `--min-mapq <INT>` (default: `20`)
Minimum mapping quality for reads to contribute.

### `--min-baseq <INT>` (default: `20`)
Minimum base quality for bases to contribute.

### `--min-depth <INT>` (default: `10`)
Minimum depth (count of accepted A/C/G/T observations) required for a confident call.

- In `strict/iupac`, depth is checked first; if below => fallback.
- In `majority`, low depth still triggers fallback (per current implementation).

### `--min-af <FLOAT>` (default: `0.7`)
Minimum allele fraction required in `strict` (and also used as “high-confidence” threshold in `iupac`).

### `--min-iupac-af <FLOAT>` (default: `0.2`)
In `iupac` mode: minimum AF required for the second allele to be considered present and generate an IUPAC mixture call.

### `--ignore-duplicates <true|false>` (default: `true`)
If true, PCR duplicates do not contribute.

---

## Indel Calling / Application

> Indels change output sequence length. bam2cons applies indels **conservatively** to avoid false length changes.

### `--enable-indel <true|false>` (default: `true`)
Turn on/off applying small insertions/deletions inferred from pileup.

### `--max-indel-len <INT>` (default: `50`)
Maximum indel length to apply. Larger indels are rejected.

### `--min-indel-depth <INT>` (default: `6`)
Minimum number of reads (after filters) supporting an indel to be considered.

### `--min-indel-af <FLOAT>` (default: `0.7`)
Minimum indel allele fraction among indel-supporting reads required to apply.

### `--indel-end-dist <INT>` (default: `10`)
Reject indel support from reads if the indel anchor is too close to read ends.
This reduces alignment artifacts near ends.

### `--indel-flank <INT>` (default: `5`)
For reads supporting an indel: require that base qualities in `[qpos-flank, qpos+flank]` are all ≥ `--min-baseq`.

### `--require-indel-strand-balance <true|false>` (default: `false`)
If true, indel must have support on both strands (plus and minus) to be applied.
- More reliable
- May lose real indels under strong strand bias

### `--indel-top2-ratio <FLOAT>` (default: `3.0`)
Indel “main-shape consistency” filter.

For each position, for insertions and deletions separately:
- find top1 (most supported) indel and top2 (second)
- require `top1 >= top2 * ratio` (top2=0 passes)
If not, indel is considered ambiguous and will be skipped.

### `--indel-mismatch-window <INT>` (default: `10`)
Indel-support read mismatch heuristic: window size (bp) around indel anchor to sample mismatches vs reference.

### `--indel-min-mismatch-bases <INT>` (default: `6`)
Minimum number of comparable bases within the mismatch window (A/C/G/T, baseQ≥min) required.
If fewer, that read’s indel support is rejected (low comparability, often near clipping/poor regions).

### `--indel-max-mismatch-rate <FLOAT>` (default: `0.20`)
Maximum mismatch rate allowed for a read to count as indel-supporting.
Higher mismatch rate often indicates misalignment → false indel.

---

## Strand Handling

### `--respect-strand <true|false>` (default: `false`)
If true:
- For BED: uses column 6 strand if present
- For GFF/GTF: uses strand column
- If strand is `-`, output reverse-complemented sequence for that target/transcript

Notes:
- For **spliced CDS mode**, strand ordering is applied at the block level and then optional reverse-complementing produces transcript-oriented output.

---

## FASTA Formatting

### `--wrap <INT>` (default: `60`)
FASTA line width. Sequence is wrapped every N characters.
- Format-only (no effect on consensus)
- To effectively disable wrapping: set a very large value (e.g. `--wrap 1000000000`)

---

## Multi-sample Parallelism

### `--jobs <INT>` (default: `4`)
Number of samples processed in parallel (thread pool size).
Tune based on CPU + storage throughput.

---

## File Discovery

### `--include-cram <true|false>` (default: `true`)
Include `.cram` files when scanning `--bam-dir`.

### `--include-bam <true|false>` (default: `true`)
Include `.bam` files when scanning `--bam-dir`.

### `--skip-existing <true|false>` (default: `false`)
If true, skip a sample if its output FASTA exists (and JSON exists if `--write-json` is enabled).

---

## Contig-batched Processing (Windows)

### `--batch-by-contig <true|false>` (default: `true`)
If true, merge nearby targets into larger windows per contig, fetch pileup once per window, then generate each target from that window’s evidence.
- Faster on many small targets
- Reduces random IO

### `--merge-gap <INT>` (default: `10000`)
When batching, merge two adjacent targets into one window if the gap between them ≤ merge-gap.

### `--max-window <INT>` (default: `5000000`)
When batching, do not create windows larger than this many bp (guards memory/time).

---

## JSON Output

### `--write-json <true|false>` (default: `true`)
If true, write `<sample>.json` with per-target and total QC stats.

---

## Spliced CDS Mode (GFF/GTF + `--feature CDS` only)

> In this mode, bam2cons reconstructs **spliced CDS per transcript** by concatenating CDS blocks.

### `--keep-longest <true|false>` (default: `true`)
If true:
- for each gene, keep only the transcript with the maximum total CDS length.

If false:
- output all transcripts for each gene.

### `--append-transcript <true|false>` (default: `true`)
Only meaningful when `--keep-longest false`.

If true:
- FASTA record name becomes `gene_id|transcript_id`.

If false:
- record name may remain gene-only (not recommended when outputting multiple transcripts).

Fallback behavior if transcript/mRNA identifiers are missing:
- If CDS records do not contain transcript identifiers (e.g., missing `transcript_id` / `Parent`),
  bam2cons assumes subsequent CDS belong to the same gene (file order dependent) and assigns a default transcript (e.g., `tx1`).

---

## Practical Recipes

### “Complete output” (minimize Ns)
- `--lowconf-fallback ref`
- `--mode strict`
- moderate thresholds:
  - `--min-depth 6–10`
  - `--min-af 0.7–0.8`

### “Avoid false indels” baseline
- `--min-indel-af 0.85`
- `--min-indel-depth 8`
- `--indel-top2-ratio 3.0`
- `--indel-max-mismatch-rate 0.20`

If true indels are being missed, relax in this order:
1) `--indel-max-mismatch-rate 0.25`
2) `--indel-top2-ratio 2.0`
3) then consider lowering `--min-indel-af`

---

## Known Limitations

- Indels are pileup-based and **no local realignment** is performed.
- Very complex/repetitive regions may conservatively skip indels.
- Long indels / SVs beyond `--max-indel-len` are ignored by design.
- If transcript IDs are missing, the “subsequent CDS belong to the same gene” fallback depends on **annotation file order**.
