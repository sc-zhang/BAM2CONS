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
  - One `sample.fa` per input BAM/CRAM
  - One `sample.json` QC/stats per sample
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

### 1) BED targets + directory of BAMs

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
  --max-indel-len 50 \
  --indel-end-dist 12 --indel-flank 6 \
  --indel-top2-ratio 3.0 \
  --indel-mismatch-window 10 --indel-min-mismatch-bases 6 --indel-max-mismatch-rate 0.20 \
  --batch-by-contig --merge-gap 10000 --max-window 5000000 \
  --write-json \
  --jobs 8
```

Outputs:
- `out/<sample>.fa`
- `out/<sample>.json`

---

## Spliced CDS Reconstruction (GFF/GTF)

### Recommended: per gene keep the longest transcript (default)

Use `--feature CDS` with a GFF3 or GTF:

```bash
./target/release/bam2cons \
  --ref-fa ref.fa \
  --gtf annot.gtf \
  --feature CDS \
  --bam-dir ./bams \
  --out-dir ./out \
  --mode strict \
  --lowconf-fallback ref \
  --keep-longest \
  --write-json \
  --jobs 8
```

**Behavior**
- Reads all `CDS` records.
- Groups CDS blocks by **gene → transcript**.
- Reconstructs **spliced CDS** per transcript by concatenating blocks in transcript order:
  - `+` strand: ascending coordinates
  - `-` strand: descending coordinates
- Applies **phase/frame trimming** (`0/1/2`) at the **first CDS block** of the transcript.
- With `--keep-longest` (default), outputs **one CDS per gene** (the transcript with maximal total CDS length).

### Output all transcripts per gene (append transcript ID)

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

Record naming:
- keep-longest=true: `gene_id`
- keep-longest=false + append-transcript=true: `gene_id|transcript_id`

> If `--append-transcript false`, record name may be gene-only (not recommended if you output multiple transcripts).

### If transcript/mRNA info is missing

If CDS entries do not provide transcript identifiers (e.g., missing `transcript_id`/`Parent`):
- bam2cons assumes **subsequent CDS belong to the same gene** (file order matters)
- They are assigned to a default transcript (e.g., `tx1`) under that gene.

---

## File Formats

### Reference FASTA

Must have `.fai`:

```bash
samtools faidx ref.fa
```

### Targets

#### BED
- 0-based, half-open `[start, end)`
- Minimal: `chrom  start  end`
- Optional: `name` (col 4), `strand` (col 6)

```text
chr1    1000    2000    geneA   0   +
chr2    500     900     geneB   0   -
```

#### GFF3 / GTF
- 1-based inclusive coordinates
- `--feature` selects which feature becomes targets (e.g., `gene`, `exon`, `CDS`)
- **CDS splicing mode** is enabled when `--feature CDS` is used with `--gff/--gtf`

---

## Output

### FASTA

Each record header includes:
- sample
- target name (or spliced transcript name)
- coordinates / run parameters (mode, thresholds, indel filters)

Notes:
- `--wrap N` controls FASTA line wrapping (format only; does not affect consensus).
- `--respect-strand` will reverse-complement outputs for `-` strand targets (BED col6 / GFF/GTF strand).

### JSON QC / Stats

Each `<sample>.json` includes:
- run parameters
- per-target statistics
- totals

Selected per-target fields:
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

Quick interpretation:
- High `ref_fallback_bases` / `low_depth_bases`: insufficient coverage → consider adjusting `--min-depth` or filtering targets.
- High `indel_rejected_mismatch` / `skipped_indel_topratio`: ambiguous/repetitive region → indels are conservatively skipped (length preserved).

---

## Parameter Guide (Practical)

### Keep output “complete” (fewer Ns)
- `--lowconf-fallback ref` (default)
- In strict mode:
  - `--min-depth 6–10`
  - `--min-af 0.7–0.8`

### Avoid false indels (recommended baseline)
- `--min-indel-af 0.85`
- `--min-indel-depth 8`
- `--indel-top2-ratio 3.0`
- `--indel-max-mismatch-rate 0.20`
- `--indel-end-dist 10–12`
- `--max-indel-len 50`

If true indels are being missed:
- First relax:
  - `--indel-max-mismatch-rate 0.25`
  - `--indel-top2-ratio 2.0`
- Then consider relaxing:
  - `--min-indel-af` (e.g. 0.85 → 0.75)

---

## Known Limitations

- Indels are pileup-based and **no local realignment** is performed.
- Very complex/repetitive regions may conservatively skip indels.
- Long indels / SVs beyond `--max-indel-len` are ignored by design.
- If transcript IDs are missing, the “subsequent CDS belong to the same gene” fallback depends on **annotation file order**.
