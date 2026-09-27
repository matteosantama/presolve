//! Benchmark instances on disk: each subdirectory of the data directory is a
//! family, and each MPS file in it is an instance. Every file stays in the
//! repository; runs select a subset with per-family size limits.

use std::io;
use std::path::{Path, PathBuf};

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct Instance {
    pub family: String,
    pub name: String,
    pub path: PathBuf,
    /// Size of the file on disk, compressed when the file is.
    pub bytes: u64,
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
            let bytes = std::fs::metadata(&path)?.len();
            instances.push(Instance {
                family: family.clone(),
                name: name.to_string(),
                path: path.clone(),
                bytes,
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

/// A size limit for one family, in whole MiB of file size on disk.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SizeLimit {
    pub family: String,
    pub mib: u64,
}
impl std::str::FromStr for SizeLimit {
    type Err = String;
    /// Parse `FAMILY=MIB`, as in `miplib=4`.
    fn from_str(s: &str) -> Result<Self, String> {
        let (family, mib) = s
            .split_once('=')
            .ok_or_else(|| format!("expected FAMILY=MIB, got {s:?}"))?;
        let mib = mib
            .parse()
            .map_err(|_| format!("expected a whole number of MiB in {s:?}"))?;
        if family.is_empty() {
            return Err(format!("expected FAMILY=MIB, got {s:?}"));
        }
        Ok(Self {
            family: family.to_string(),
            mib,
        })
    }
}

/// Split `instances` into those within their family's limit and those above
/// it. A family without a limit keeps every instance. A limit naming a family
/// with no discovered instances is an error, since it is most likely a typo.
pub fn within_limits(
    instances: Vec<Instance>,
    limits: &[SizeLimit],
) -> Result<(Vec<Instance>, Vec<Instance>), String> {
    for limit in limits {
        if !instances.iter().any(|i| i.family == limit.family) {
            return Err(format!("size limit for unknown family {}", limit.family));
        }
    }
    Ok(instances.into_iter().partition(|instance| {
        limits
            .iter()
            .filter(|l| l.family == instance.family)
            .all(|l| instance.bytes <= l.mib << 20)
    }))
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
        let sized = discover(&data, &[]).unwrap();
        assert!(sized.iter().all(|i| i.bytes == 0));
        std::fs::write(data.join("a").join("z.mps"), "").unwrap();
        let clash = discover(&data, &[]).unwrap_err();
        assert!(
            clash.to_string().ends_with("are the same instance"),
            "{clash}"
        );
        std::fs::remove_dir_all(&data).unwrap();
    }

    fn instance(family: &str, name: &str, bytes: u64) -> Instance {
        Instance {
            family: family.into(),
            name: name.into(),
            path: PathBuf::from(format!("{family}/{name}.mps.gz")),
            bytes,
        }
    }

    #[test]
    fn size_limits_apply_per_family() {
        let mib = 1 << 20;
        let instances = vec![
            instance("big", "small", 4 * mib),
            instance("big", "large", 4 * mib + 1),
            instance("other", "large", 100 * mib),
        ];
        let limit: SizeLimit = "big=4".parse().unwrap();
        assert_eq!(
            limit,
            SizeLimit {
                family: "big".into(),
                mib: 4
            }
        );
        let (kept, excluded) = within_limits(instances.clone(), &[limit]).unwrap();
        let ids = |v: &[Instance]| v.iter().map(Instance::id).collect::<Vec<_>>();
        assert_eq!(ids(&kept), ["big/small", "other/large"]);
        assert_eq!(ids(&excluded), ["big/large"]);
        let (kept, excluded) = within_limits(instances.clone(), &[]).unwrap();
        assert_eq!((kept.len(), excluded.len()), (3, 0));
        let typo = within_limits(instances, &["bgi=4".parse().unwrap()]).unwrap_err();
        assert_eq!(typo, "size limit for unknown family bgi");
        for bad in ["big", "=4", "big=four", "big=-1"] {
            assert!(bad.parse::<SizeLimit>().is_err(), "{bad}");
        }
    }
}
