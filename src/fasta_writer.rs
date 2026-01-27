use anyhow::{Context, Result};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

pub struct FastaWriter {
    w: BufWriter<File>,
    wrap: usize,
}

impl FastaWriter {
    pub fn create(path: &Path, wrap: usize) -> Result<Self> {
        let f = File::create(path).with_context(|| format!("create {}", path.display()))?;
        Ok(Self { w: BufWriter::new(f), wrap: wrap.max(1) })
    }

    pub fn write_record(&mut self, header: &str, seq: &[u8]) -> Result<()> {
        writeln!(self.w, ">{}", header)?;
        for chunk in seq.chunks(self.wrap) {
            self.w.write_all(chunk)?;
            writeln!(self.w)?;
        }
        Ok(())
    }

    pub fn flush(&mut self) -> Result<()> {
        self.w.flush()?;
        Ok(())
    }
}
