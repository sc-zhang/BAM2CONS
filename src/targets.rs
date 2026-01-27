use anyhow::{anyhow, bail, Context, Result};
use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

#[derive(Clone, Debug)]
pub struct Target {
    pub id: usize,
    pub chrom: String,
    pub start: i64, // 0-based, inclusive
    pub end: i64,   // 0-based, exclusive
    pub name: String,
    pub strand: Option<char>,
}

#[derive(Clone, Debug)]
pub struct CdsBlock {
    pub start: i64, // 0-based inclusive
    pub end: i64,   // 0-based exclusive
    pub phase: u8,  // 0/1/2 (GFF phase or GTF frame); default 0
}

#[derive(Clone, Debug)]
pub struct SplicedTarget {
    pub id: usize,
    pub chrom: String,
    pub strand: Option<char>,
    pub gene_id: String,
    pub transcript_id: String,
    pub name: String,          // for FASTA header display
    pub blocks: Vec<CdsBlock>, // in genome coordinates; we will order later
}

pub fn read_bed(path: &Path) -> Result<Vec<Target>> {
    let f = File::open(path).with_context(|| format!("open BED: {}", path.display()))?;
    let r = BufReader::new(f);

    let mut out = Vec::new();
    for (lineno, line) in r.lines().enumerate() {
        let line = line?;
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 3 {
            bail!("BED parse error at line {}: need >=3 columns", lineno + 1);
        }
        let chrom = cols[0].to_string();
        let start: i64 = cols[1]
            .parse()
            .with_context(|| format!("BED start at line {}", lineno + 1))?;
        let end: i64 = cols[2]
            .parse()
            .with_context(|| format!("BED end at line {}", lineno + 1))?;
        if end <= start {
            return Err(anyhow!(
                "BED invalid interval at line {}: start>=end",
                lineno + 1
            ));
        }
        let name = if cols.len() >= 4 && !cols[3].is_empty() {
            cols[3].to_string()
        } else {
            format!("bed:{}:{}-{}", chrom, start, end)
        };
        let strand = if cols.len() >= 6 {
            match cols[5].chars().next() {
                Some('+') => Some('+'),
                Some('-') => Some('-'),
                _ => None,
            }
        } else {
            None
        };

        out.push(Target {
            id: out.len(),
            chrom,
            start,
            end,
            name,
            strand,
        });
    }
    Ok(out)
}

/// Lightweight GFF3/GTF parser, used only to define target intervals (NOT to make individual annotation).
pub fn read_gff_like(path: &Path, feature: &str, is_gtf: bool) -> Result<Vec<Target>> {
    let f = File::open(path).with_context(|| {
        format!(
            "open {}: {}",
            if is_gtf { "GTF" } else { "GFF" },
            path.display()
        )
    })?;
    let r = BufReader::new(f);

    let mut out = Vec::new();
    for (lineno, line) in r.lines().enumerate() {
        let line = line?;
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 9 {
            bail!(
                "{} parse error at line {}: need 9 columns",
                if is_gtf { "GTF" } else { "GFF" },
                lineno + 1
            );
        }
        let chrom = cols[0].to_string();
        let ftype = cols[2];
        if ftype != feature {
            continue;
        }

        let start1: i64 = cols[3]
            .parse()
            .with_context(|| format!("start at line {}", lineno + 1))?;
        let end1: i64 = cols[4]
            .parse()
            .with_context(|| format!("end at line {}", lineno + 1))?;
        if end1 < start1 {
            bail!("invalid interval at line {}: end<start", lineno + 1);
        }
        let start = start1 - 1;
        let end = end1;

        let strand = match cols[6].chars().next() {
            Some('+') => Some('+'),
            Some('-') => Some('-'),
            _ => None,
        };

        let attrs = cols[8];
        let name = extract_name(attrs, is_gtf)
            .unwrap_or_else(|| format!("{}:{}:{}-{}", feature, chrom, start, end));

        out.push(Target {
            id: out.len(),
            chrom,
            start,
            end,
            name,
            strand,
        });
    }

    if out.is_empty() {
        return Err(anyhow!(
            "No targets extracted from {} with feature='{}'. Check --feature or file format.",
            path.display(),
            feature
        ));
    }
    Ok(out)
}

fn extract_name(attrs: &str, is_gtf: bool) -> Option<String> {
    if is_gtf {
        let keys = ["transcript_id", "gene_id", "gene_name", "ID", "Name"];
        for k in keys {
            if let Some(v) = gtf_get(attrs, k) {
                return Some(v);
            }
        }
        None
    } else {
        let keys = ["ID", "Name", "gene_id", "transcript_id"];
        for k in keys {
            if let Some(v) = gff_get(attrs, k) {
                return Some(v);
            }
        }
        None
    }
}

fn gff_get(attrs: &str, key: &str) -> Option<String> {
    for part in attrs.split(';') {
        let part = part.trim();
        if part.is_empty() {
            continue;
        }
        let mut it = part.splitn(2, '=');
        let k = it.next()?.trim();
        let v = it.next()?.trim();
        if k == key && !v.is_empty() {
            return Some(v.to_string());
        }
    }
    None
}

fn gtf_get(attrs: &str, key: &str) -> Option<String> {
    for part in attrs.split(';') {
        let part = part.trim();
        if part.is_empty() {
            continue;
        }
        if part.starts_with(key) {
            let rest = part[key.len()..].trim();
            let rest = rest.strip_prefix('"').unwrap_or(rest);
            let val = rest.split('"').next().unwrap_or(rest).trim();
            if !val.is_empty() {
                return Some(val.to_string());
            }
        }
    }
    None
}

fn parse_phase(col8: &str) -> u8 {
    match col8.trim() {
        "0" => 0,
        "1" => 1,
        "2" => 2,
        _ => 0,
    }
}

/// Very small helper: fetch attribute for both GFF3 (k=v) and GTF (k "v")
fn attr_get(attrs: &str, key: &str, is_gtf: bool) -> Option<String> {
    if is_gtf {
        for part in attrs.split(';') {
            let part = part.trim();
            if part.is_empty() {
                continue;
            }
            if part.starts_with(key) {
                let rest = part[key.len()..].trim();
                let rest = rest.strip_prefix('"').unwrap_or(rest);
                let val = rest.split('"').next().unwrap_or(rest).trim();
                if !val.is_empty() {
                    return Some(val.to_string());
                }
            }
        }
        None
    } else {
        for part in attrs.split(';') {
            let part = part.trim();
            if part.is_empty() {
                continue;
            }
            let mut it = part.splitn(2, '=');
            let k = it.next()?.trim();
            let v = it.next()?.trim();
            if k == key && !v.is_empty() {
                return Some(v.to_string());
            }
        }
        None
    }
}

/// For GFF3: Parent often points to transcript/mRNA; gene may be Parent of mRNA.
/// We use Parent as transcript key when present; gene key fallback to gene_id/ID/Name/current gene.
pub fn read_spliced_cds(
    path: &Path,
    is_gtf: bool,
    keep_longest: bool,
    append_transcript: bool,
) -> Result<Vec<SplicedTarget>> {
    let f = File::open(path).with_context(|| format!("open {}", path.display()))?;
    let r = BufReader::new(f);

    // key: (gene_id, transcript_id) -> blocks + metadata
    #[derive(Default)]
    struct Acc {
        chrom: String,
        strand: Option<char>,
        blocks: Vec<CdsBlock>,
    }
    let mut acc: HashMap<(String, String), Acc> = HashMap::new();

    // "current gene" fallback when transcript info missing and file is ordered by gene
    let mut current_gene: Option<String> = None;
    let mut gene_counter: usize = 0;
    let mut fallback_tx_counter: BTreeMap<String, usize> = BTreeMap::new(); // per gene

    for line in r.lines() {
        let line = line?;
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 9 {
            continue;
        }

        let chrom = cols[0].to_string();
        let ftype = cols[2];
        let start1: i64 = cols[3].parse().ok().unwrap_or(0);
        let end1: i64 = cols[4].parse().ok().unwrap_or(0);
        if end1 < start1 {
            continue;
        }
        let start = start1 - 1;
        let end = end1;
        let strand = match cols[6].chars().next() {
            Some('+') => Some('+'),
            Some('-') => Some('-'),
            _ => None,
        };
        let phase = parse_phase(cols[7]);
        let attrs = cols[8];

        // Track current gene anchor (for “no transcript/mRNA => subsequent CDS belong to same gene”)
        if ftype == "gene" {
            // Prefer explicit gene identifiers if present
            let gid = attr_get(attrs, "gene_id", is_gtf)
                .or_else(|| attr_get(attrs, "ID", is_gtf))
                .or_else(|| attr_get(attrs, "Name", is_gtf))
                .unwrap_or_else(|| {
                    gene_counter += 1;
                    format!("gene{}", gene_counter)
                });
            current_gene = Some(gid);
            continue;
        }

        if ftype != "CDS" {
            continue;
        }

        // Gene id
        let gene_id = attr_get(attrs, "gene_id", is_gtf)
            .or_else(|| attr_get(attrs, "gene", is_gtf))
            .or_else(|| attr_get(attrs, "Name", is_gtf)) // sometimes abused
            .or_else(|| current_gene.clone())
            .unwrap_or_else(|| {
                // very last resort
                gene_counter += 1;
                format!("gene{}", gene_counter)
            });

        // Transcript id (preferred)
        let mut transcript_id = attr_get(attrs, "transcript_id", is_gtf)
            .or_else(|| attr_get(attrs, "Parent", is_gtf)) // GFF3 CDS Parent=transcript/mRNA
            .or_else(|| attr_get(attrs, "transcript", is_gtf));

        if transcript_id.is_none() {
            // fallback transcript for this gene (ordered-file assumption)
            let c = fallback_tx_counter.entry(gene_id.clone()).or_insert(0);
            if *c == 0 {
                *c = 1;
            }
            transcript_id = Some(format!("tx{}", *c));
        }

        let tid = transcript_id.unwrap();

        let key = (gene_id.clone(), tid.clone());
        let entry = acc.entry(key).or_insert_with(|| Acc {
            chrom: chrom.clone(),
            strand,
            blocks: Vec::new(),
        });

        // ensure consistent chrom/strand per transcript (best effort)
        entry.chrom = chrom.clone();
        entry.strand = strand.or(entry.strand);

        entry.blocks.push(CdsBlock { start, end, phase });
    }

    // Convert to SplicedTarget list and optionally keep-longest per gene
    let mut per_gene: BTreeMap<String, Vec<SplicedTarget>> = BTreeMap::new();
    let mut next_id: usize = 0;

    for ((gene_id, transcript_id), a) in acc {
        let total_len: i64 = a.blocks.iter().map(|b| b.end - b.start).sum();
        let name = if keep_longest {
            gene_id.clone()
        } else if append_transcript {
            format!("{}|{}", gene_id, transcript_id)
        } else {
            gene_id.clone()
        };

        let st = SplicedTarget {
            id: next_id,
            chrom: a.chrom,
            strand: a.strand,
            gene_id: gene_id.clone(),
            transcript_id: transcript_id.clone(),
            name,
            blocks: a.blocks,
        };
        next_id += 1;
        per_gene.entry(gene_id).or_default().push(st);

        // you can store total_len somewhere; easiest is recompute when selecting longest
        let _ = total_len;
    }

    let mut out = Vec::new();
    if keep_longest {
        for (_gid, mut v) in per_gene {
            // pick transcript with max total CDS length; tie-breaker: lexicographic transcript_id
            v.sort_by(|a, b| {
                let la: i64 = a.blocks.iter().map(|x| x.end - x.start).sum();
                let lb: i64 = b.blocks.iter().map(|x| x.end - x.start).sum();
                lb.cmp(&la)
                    .then_with(|| a.transcript_id.cmp(&b.transcript_id))
            });
            if let Some(first) = v.into_iter().next() {
                out.push(first);
            }
        }
    } else {
        for (_gid, mut v) in per_gene {
            // stable order by transcript_id
            v.sort_by(|a, b| a.transcript_id.cmp(&b.transcript_id));
            out.extend(v);
        }
    }

    Ok(out)
}
