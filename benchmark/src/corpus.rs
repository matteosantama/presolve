//! Benchmark instances on disk: each subdirectory of the data directory is a
//! family, and each MPS file in it is an instance.

use std::io;
use std::path::{Path, PathBuf};

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct Instance {
    pub family: String,
    pub name: String,
    pub path: PathBuf,
}
impl Instance {
    /// `family/name`, the key that identifies an instance across snapshots.
    pub fn id(&self) -> String {
        format!("{}/{}", self.family, self.name)
    }
}

/// Instances under `data`, sorted by family then name. An empty `families`
/// selects every family; naming a missing family is an error. Files other
/// than `.mps` and `.mps.gz` are ignored, and two files with the same
/// instance name in one family are an error.
pub fn discover(data: &Path, families: &[String]) -> io::Result<Vec<Instance>> {
    let mut selected: Vec<String> = families.to_vec();
    if selected.is_empty() {
        for entry in std::fs::read_dir(data)? {
            let entry = entry?;
            if entry.file_type()?.is_dir() {
                selected.push(entry.file_name().to_string_lossy().into_owned());
            }
        }
    }
    selected.sort();
    selected.dedup();
    let mut instances = Vec::new();
    for family in selected {
        let dir = data.join(&family);
        let entries = std::fs::read_dir(&dir).map_err(|e| {
            io::Error::new(
                e.kind(),
                format!("family {family} ({}): {e}", dir.display()),
            )
        })?;
        for entry in entries {
            let path = entry?.path();
            let Some(file) = path.file_name().and_then(|f| f.to_str()) else {
                continue;
            };
            let Some(name) = file
                .strip_suffix(".mps.gz")
                .or_else(|| file.strip_suffix(".mps"))
            else {
                continue;
            };
            instances.push(Instance {
                family: family.clone(),
                name: name.to_string(),
                path: path.clone(),
            });
        }
    }
    instances.sort();
    if let Some(pair) = instances
        .windows(2)
        .find(|p| (&p[0].family, &p[0].name) == (&p[1].family, &p[1].name))
    {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "{} and {} are the same instance",
                pair[0].path.display(),
                pair[1].path.display()
            ),
        ));
    }
    Ok(instances)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn scratch(name: &str) -> PathBuf {
        let dir =
            std::env::temp_dir().join(format!("benchmark-corpus-{name}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        dir
    }

    #[test]
    fn families_are_directories_and_instances_are_mps_files() {
        let data = scratch("discover");
        for (family, file) in [
            ("b", "y.mps.gz"),
            ("b", "x.mps"),
            ("b", "PILOT.JA.mps.gz"),
            ("a", "z.mps.gz"),
            ("a", ".DS_Store"),
            ("a", "notes.txt"),
        ] {
            std::fs::create_dir_all(data.join(family)).unwrap();
            std::fs::write(data.join(family).join(file), "").unwrap();
        }
        std::fs::write(data.join("README"), "").unwrap();
        let ids = |families: &[&str]| -> Vec<String> {
            let families: Vec<String> = families.iter().map(|f| f.to_string()).collect();
            discover(&data, &families)
                .unwrap()
                .iter()
                .map(Instance::id)
                .collect()
        };
        assert_eq!(ids(&[]), ["a/z", "b/PILOT.JA", "b/x", "b/y"]);
        assert_eq!(ids(&["b", "a", "b"]), ["a/z", "b/PILOT.JA", "b/x", "b/y"]);
        let missing = discover(&data, &["c".to_string()]).unwrap_err();
        assert!(missing.to_string().starts_with("family c"), "{missing}");
        std::fs::write(data.join("a").join("z.mps"), "").unwrap();
        let clash = discover(&data, &[]).unwrap_err();
        assert!(
            clash.to_string().ends_with("are the same instance"),
            "{clash}"
        );
        std::fs::remove_dir_all(&data).unwrap();
    }
}
