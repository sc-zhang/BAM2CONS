pub fn base_to_index(b: u8) -> Option<usize> {
    match b {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

pub fn index_to_base(i: usize) -> u8 {
    match i {
        0 => b'A',
        1 => b'C',
        2 => b'G',
        3 => b'T',
        _ => b'N',
    }
}

pub fn iupac_code(a: u8, b: u8) -> Option<u8> {
    let (x, y) = if a <= b { (a, b) } else { (b, a) };
    match (x, y) {
        (b'A', b'G') => Some(b'R'),
        (b'C', b'T') => Some(b'Y'),
        (b'C', b'G') => Some(b'S'),
        (b'A', b'T') => Some(b'W'),
        (b'G', b'T') => Some(b'K'),
        (b'A', b'C') => Some(b'M'),
        _ => None,
    }
}

pub fn is_acgt(b: u8) -> bool {
    matches!(b, b'A' | b'C' | b'G' | b'T')
}

/// Check that qualities in [qpos-flank, qpos+flank] (inclusive) are all >= minq.
pub fn passes_flank_qual(qual: &[u8], qpos: usize, flank: usize, minq: u8) -> bool {
    if qual.is_empty() {
        return false;
    }
    let start = qpos.saturating_sub(flank);
    let end = (qpos + flank).min(qual.len().saturating_sub(1));
    for &q in &qual[start..=end] {
        if q < minq {
            return false;
        }
    }
    true
}

pub fn revcomp_in_place(seq: &mut [u8]) {
    seq.reverse();
    for b in seq.iter_mut() {
        *b = match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            b'R' => b'Y',
            b'Y' => b'R',
            b'S' => b'S',
            b'W' => b'W',
            b'K' => b'M',
            b'M' => b'K',
            b'N' => b'N',
            other => other,
        };
    }
}
