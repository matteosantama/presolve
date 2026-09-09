use crate::Result;
use flate2::read::GzDecoder;
use mps::Parser;
use mps::types::{BoundType, ObjectiveSense, RowType, WideLine};
use presolve::{
    matrix::CscMatrix,
    problem::{Bounds, Constraint, ProblemData},
};
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::Read;
use std::path::{Path, PathBuf};

pub const SUITES: [&str; 2] = ["netlib", "maros-meszaros"];

pub fn problems(suite: &str) -> Result<Vec<PathBuf>> {
    let directory = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("data")
        .join(suite);
    let mut paths = fs::read_dir(directory)?
        .map(|entry| entry.map(|e| e.path()))
        .collect::<std::io::Result<Vec<_>>>()?;
    paths.retain(|p| {
        p.file_name()
            .unwrap()
            .to_string_lossy()
            .ends_with(".mps.gz")
    });
    paths.sort();
    Ok(paths)
}

pub fn read(path: &Path) -> Result<ProblemData> {
    let mut text = String::new();
    GzDecoder::new(File::open(path)?).read_to_string(&mut text)?;
    parse(&text).map_err(|e| format!("{}: {e}", path.display()).into())
}

// Use the first named RHS, RANGES and BOUNDS sets, as in the corpus's solver runs.
fn first_set<'a>(lines: &'a [WideLine<'a, f64>]) -> impl Iterator<Item = (&'a str, f64)> {
    let name = lines.first().map(|line| line.name);
    lines
        .iter()
        .filter(move |line| Some(line.name) == name)
        .flat_map(|line| std::iter::once(&line.first_pair).chain(line.second_pair.as_ref()))
        .map(|pair| (pair.row_name, pair.value))
}

fn parse(text: &str) -> Result<ProblemData> {
    let source = Parser::<f64>::parse(text).map_err(|e| e.to_string())?;
    let objective = source
        .objective_name
        .or_else(|| {
            source
                .rows
                .iter()
                .find(|r| r.row_type == RowType::Nr)
                .map(|r| r.row_name)
        })
        .ok_or("missing objective row")?;
    let rows: Vec<_> = source
        .rows
        .iter()
        .filter(|r| r.row_type != RowType::Nr)
        .collect();
    let row_index: HashMap<_, _> = rows
        .iter()
        .enumerate()
        .map(|(i, r)| (r.row_name, i))
        .collect();
    let mut columns = HashMap::new();
    for line in &source.columns {
        let next = columns.len();
        columns.entry(line.name).or_insert(next);
    }
    let n = columns.len();
    let sign = if source.objective_sense == Some(ObjectiveSense::Max) {
        -1.
    } else {
        1.
    };
    let mut c = vec![0.; n];
    let (mut ai, mut aj, mut av) = (Vec::new(), Vec::new(), Vec::new());
    for line in &source.columns {
        let j = columns[line.name];
        for pair in std::iter::once(&line.first_pair).chain(line.second_pair.as_ref()) {
            if pair.row_name == objective {
                c[j] += sign * pair.value;
            } else if let Some(&i) = row_index.get(pair.row_name) {
                ai.push(i);
                aj.push(j);
                av.push(pair.value);
            }
        }
    }
    let mut rhs = vec![0.; rows.len()];
    let mut objective_constant = 0.;
    for (name, value) in first_set(source.rhs.as_deref().unwrap_or_default()) {
        if name == objective {
            objective_constant = -sign * value;
        } else if let Some(&i) = row_index.get(name) {
            rhs[i] = value;
        }
    }
    let mut bounds: Vec<_> = rows
        .iter()
        .zip(&rhs)
        .map(|(row, &b)| match row.row_type {
            RowType::Eq => Bounds::fixed(b),
            RowType::Leq => Bounds {
                upper: b,
                ..Bounds::FREE
            },
            RowType::Geq => Bounds {
                lower: b,
                ..Bounds::FREE
            },
            RowType::Nr => unreachable!(),
        })
        .collect();
    for (name, value) in first_set(source.ranges.as_deref().unwrap_or_default()) {
        if let Some(&i) = row_index.get(name) {
            match rows[i].row_type {
                RowType::Leq => bounds[i].lower = rhs[i] - value.abs(),
                RowType::Geq => bounds[i].upper = rhs[i] + value.abs(),
                RowType::Eq => {
                    bounds[i].lower = rhs[i] + value.min(0.);
                    bounds[i].upper = rhs[i] + value.max(0.);
                }
                RowType::Nr => unreachable!(),
            }
        }
    }
    let mut variable_bounds = vec![
        Bounds {
            lower: 0.,
            upper: f64::INFINITY
        };
        n
    ];
    let lines = source.bounds.as_deref().unwrap_or_default();
    let first = lines.first().map(|b| b.bound_name);
    for bound in lines.iter().filter(|b| Some(b.bound_name) == first) {
        let b = &mut variable_bounds[columns[bound.column_name]];
        match bound.bound_type {
            BoundType::Lo => b.lower = bound.value.ok_or("LO needs a value")?,
            BoundType::Up => b.upper = bound.value.ok_or("UP needs a value")?,
            BoundType::Fx => *b = Bounds::fixed(bound.value.ok_or("FX needs a value")?),
            BoundType::Fr => *b = Bounds::FREE,
            BoundType::Mi => {
                *b = Bounds {
                    upper: 0.,
                    ..Bounds::FREE
                }
            }
            BoundType::Pl => {
                *b = Bounds {
                    lower: 0.,
                    ..Bounds::FREE
                }
            }
            _ => return Err("benchmark loader supports continuous variables only".into()),
        }
    }
    // QPS files may store either triangle or both. Mirrored entries represent
    // the same Hessian coefficient and must not be added twice.
    let p = if let Some(terms) = &source.quadratic_objective {
        let mut pairs: HashMap<(usize, usize), [Option<f64>; 2]> = HashMap::new();
        for term in terms {
            let i = columns[term.var1];
            let j = columns[term.var2];
            let pair = pairs.entry((i.min(j), i.max(j))).or_default();
            *pair[usize::from(i > j)].get_or_insert(0.) += sign * term.coefficient;
        }
        let (mut pi, mut pj, mut pv) = (Vec::new(), Vec::new(), Vec::new());
        for ((i, j), pair) in pairs {
            let value = match pair {
                [Some(a), Some(b)] => {
                    if (a - b).abs() > 1e-12 * (1. + a.abs().max(b.abs())) {
                        return Err("inconsistent symmetric Hessian entries".into());
                    }
                    0.5 * a + 0.5 * b
                }
                [Some(a), None] | [None, Some(a)] => a,
                [None, None] => unreachable!(),
            };
            pi.push(i);
            pj.push(j);
            pv.push(value);
        }
        Some(CscMatrix::from_triplets(n, n, pi, pj, pv)?)
    } else {
        None
    };
    Ok(ProblemData {
        p,
        c,
        objective_constant,
        a: CscMatrix::from_triplets(rows.len(), n, ai, aj, av)?,
        rows: bounds.into_iter().map(Constraint::Linear).collect(),
        variable_bounds,
        cones: vec![],
    })
}
