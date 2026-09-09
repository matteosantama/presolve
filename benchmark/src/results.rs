use crate::{Kind, Result};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fs::{self, File, OpenOptions};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

#[derive(Debug, Serialize, Deserialize)]
pub struct Metadata {
    pub version: u32,
    pub name: String,
    pub kind: Kind,
    pub trials: usize,
    pub created_unix_seconds: u64,
    pub machine: String,
    pub settings: String,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Counts {
    pub variables: usize,
    pub linear_rows: usize,
    pub conic_rows: usize,
    pub a_nonzeros: usize,
    pub g_nonzeros: usize,
    pub p_nonzeros: usize,
}

impl From<presolve::result::Size> for Counts {
    fn from(s: presolve::result::Size) -> Self {
        Self {
            variables: s.variables,
            linear_rows: s.linear_rows,
            conic_rows: s.conic_rows,
            a_nonzeros: s.a_nonzeros,
            g_nonzeros: s.g_nonzeros,
            p_nonzeros: s.p_nonzeros,
        }
    }
}

impl Counts {
    pub fn values(&self) -> [usize; 6] {
        [
            self.variables,
            self.linear_rows,
            self.conic_rows,
            self.a_nonzeros,
            self.g_nonzeros,
            self.p_nonzeros,
        ]
    }
}

#[derive(Debug, Serialize, Deserialize)]
pub struct Measurement {
    // None for size runs: those runs intentionally collect no timings.
    pub elapsed_ns: Option<u64>,
    pub outcome: String,
    pub before: Counts,
    pub after: Option<Counts>,
    pub time_limit_reached: bool,
}

impl Measurement {
    pub fn same_result(&self, other: &Self) -> bool {
        self.outcome == other.outcome && self.before == other.before && self.after == other.after
    }
}

#[derive(Debug, Serialize, Deserialize)]
pub struct Case {
    pub input_hash: String,
    pub measurements: Vec<Measurement>,
    pub errors: Vec<String>,
}

#[derive(Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case")]
enum Record {
    Header(Metadata),
    Case { id: String, result: Case },
    Complete,
}

pub struct Run {
    pub metadata: Metadata,
    pub cases: BTreeMap<String, Case>,
    pub complete: bool,
}

pub struct RunWriter {
    pub path: PathBuf,
    writer: BufWriter<File>,
}

pub fn path(root: &Path, kind: Kind, name: &str) -> Result<PathBuf> {
    if name.is_empty()
        || name == "."
        || name == ".."
        || !name
            .bytes()
            .all(|c| c.is_ascii_alphanumeric() || b"-_.".contains(&c))
    {
        return Err("run names may contain only letters, numbers, '-', '_' and '.'".into());
    }
    Ok(root.join(kind.as_str()).join(format!("{name}.jsonl")))
}

impl RunWriter {
    pub fn create(root: &Path, metadata: Metadata) -> Result<Self> {
        let path = path(root, metadata.kind, &metadata.name)?;
        fs::create_dir_all(path.parent().unwrap())?;
        let file = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&path)
            .map_err(|e| {
                format!(
                    "cannot create {} (existing runs are never overwritten): {e}",
                    path.display()
                )
            })?;
        let mut out = Self {
            path,
            writer: BufWriter::new(file),
        };
        out.write(&Record::Header(metadata))?;
        Ok(out)
    }

    fn write(&mut self, record: &Record) -> Result<()> {
        serde_json::to_writer(&mut self.writer, record)?;
        writeln!(self.writer)?;
        self.writer.flush()?;
        Ok(())
    }

    pub fn case(&mut self, id: String, result: Case) -> Result<()> {
        self.write(&Record::Case { id, result })
    }

    pub fn finish(&mut self) -> Result<()> {
        self.write(&Record::Complete)
    }
}

pub fn load(root: &Path, kind: Kind, name: &str) -> Result<Run> {
    let path = path(root, kind, name)?;
    let file = File::open(&path).map_err(|e| format!("{}: {e}", path.display()))?;
    let mut records = BufReader::new(file).lines();
    let first = records.next().ok_or("empty run file")??;
    let Record::Header(metadata) = serde_json::from_str(&first)? else {
        return Err("run file must begin with metadata".into());
    };
    if metadata.version != 1
        || metadata.kind != kind
        || metadata.name != name
        || metadata.trials == 0
    {
        return Err("unsupported or inconsistent run metadata".into());
    }
    let mut run = Run {
        metadata,
        cases: BTreeMap::new(),
        complete: false,
    };
    for line in records {
        if run.complete {
            return Err("unexpected data after the completed run".into());
        }
        match serde_json::from_str(&line?)? {
            Record::Case { id, result } => {
                if run.cases.insert(id, result).is_some() {
                    return Err("duplicate case in run file".into());
                }
            }
            Record::Complete => run.complete = true,
            Record::Header(_) => return Err("duplicate run metadata".into()),
        }
    }
    Ok(run)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn run_files_preserve_partial_results_and_cannot_be_overwritten() {
        let root = std::env::temp_dir().join(format!(
            "presolve-results-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let metadata = || Metadata {
            version: 1,
            name: "example".into(),
            kind: Kind::Size,
            trials: 1,
            created_unix_seconds: 0,
            machine: "test".into(),
            settings: "test".into(),
        };
        let mut writer = RunWriter::create(&root, metadata()).unwrap();
        assert!(RunWriter::create(&root, metadata()).is_err());
        writer
            .case(
                "netlib/AFIRO/all".into(),
                Case {
                    input_hash: "hash".into(),
                    measurements: vec![],
                    errors: vec!["worker failed".into()],
                },
            )
            .unwrap();
        let partial = load(&root, Kind::Size, "example").unwrap();
        assert!(!partial.complete);
        assert_eq!(partial.cases.len(), 1);
        writer.finish().unwrap();
        assert!(load(&root, Kind::Size, "example").unwrap().complete);
        assert!(path(&root, Kind::Size, "../escape").is_err());
        drop(writer);
        fs::remove_dir_all(root).unwrap();
    }
}
