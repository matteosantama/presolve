//! MPS reader for the benchmark corpus, with QUADOBJ quadratic objectives.
//!
//! Free format is tried first. A file that fails to parse is retried in
//! fixed format, whose column positions allow spaces inside names; the error
//! from the attempt that got further is reported. Supported sections are
//! NAME, ROWS, COLUMNS, RHS, RANGES, BOUNDS, QUADOBJ, and ENDATA.
//!
//! Conventions, matching common solvers:
//! - the first `N` row is the objective and later `N` rows are ignored;
//! - an RHS entry on the objective row is the negated objective constant;
//! - a range `R` makes an `L` row `[rhs - |R|, rhs]`, a `G` row
//!   `[rhs, rhs + |R|]`, and an `E` row `[rhs, rhs + R]` or `[rhs + R, rhs]`
//!   by the sign of `R`;
//! - bounds default to `[0, inf)`; `UP` with a negative value on a column
//!   without an explicit lower bound sets the lower bound to `-inf`;
//! - QUADOBJ lists one triangle of the symmetric `Q` in `0.5 xᵀ Q x`;
//! - set names in RHS, RANGES, and BOUNDS may be omitted, and every set is
//!   read as one;
//! - integrality markers are skipped, giving the continuous relaxation.
//!
//! Duplicate matrix or QUADOBJ entries are errors: in this corpus they only
//! arise from misreading a file's format.

use presolve::Problem;
use presolve::matrix::{CscMatrix, MatrixError};
use presolve::problem::{Bounds, Constraint};
use std::collections::HashMap;
use std::fmt;
use std::io::Read;
use std::path::Path;

#[derive(Debug)]
pub enum Error {
    Io(std::io::Error),
    /// `line` is one-based; zero means the error is not tied to a line.
    Parse {
        line: usize,
        message: String,
    },
    Matrix(MatrixError),
}
impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(e) => write!(f, "{e}"),
            Self::Parse { line: 0, message } => write!(f, "{message}"),
            Self::Parse { line, message } => write!(f, "line {line}: {message}"),
            Self::Matrix(e) => write!(f, "{e}"),
        }
    }
}
impl std::error::Error for Error {}
impl From<std::io::Error> for Error {
    fn from(e: std::io::Error) -> Self {
        Self::Io(e)
    }
}
impl From<MatrixError> for Error {
    fn from(e: MatrixError) -> Self {
        Self::Matrix(e)
    }
}

/// Read an MPS file, decompressing it first when the name ends in `.gz`.
pub fn read(path: &Path) -> Result<Problem, Error> {
    let file = std::fs::File::open(path)?;
    let mut text = String::new();
    if path.extension().is_some_and(|e| e == "gz") {
        flate2::read::MultiGzDecoder::new(file).read_to_string(&mut text)?;
    } else {
        std::io::BufReader::new(file).read_to_string(&mut text)?;
    }
    parse(&text)
}

/// Parse MPS text, trying free format and then fixed format.
pub fn parse(text: &str) -> Result<Problem, Error> {
    match parse_as(text, Format::Free) {
        Ok(problem) => Ok(problem),
        Err(free) => parse_as(text, Format::Fixed).map_err(|fixed| {
            if fixed.progress() > free.progress() {
                fixed
            } else {
                free
            }
        }),
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Format {
    Free,
    Fixed,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Section {
    None,
    Rows,
    Columns,
    Rhs,
    Ranges,
    Bounds,
    QuadObj,
    End,
}

#[derive(Clone, Copy)]
enum RowKind {
    Equal,
    Less,
    Greater,
}

#[derive(Clone, Copy)]
enum Row {
    Objective,
    Free,
    Constraint(usize),
}

impl Error {
    /// How far the failed attempt read. Errors not tied to a line are found
    /// after the last line.
    fn progress(&self) -> usize {
        match self {
            Self::Parse { line, .. } if *line > 0 => *line,
            _ => usize::MAX,
        }
    }
}

fn parse_as(text: &str, format: Format) -> Result<Problem, Error> {
    let mut builder = Builder::default();
    let mut section = Section::None;
    let mut tokens = Vec::new();
    for (index, line) in text.lines().enumerate() {
        let fail = |message: String| Error::Parse {
            line: index + 1,
            message,
        };
        if line.starts_with('*') || line.trim().is_empty() {
            continue;
        }
        if !line.starts_with(char::is_whitespace) {
            let keyword = line.split_whitespace().next().unwrap_or_default();
            section = match keyword {
                "NAME" => Section::None,
                "ROWS" => Section::Rows,
                "COLUMNS" => Section::Columns,
                "RHS" => Section::Rhs,
                "RANGES" => Section::Ranges,
                "BOUNDS" => Section::Bounds,
                "QUADOBJ" => Section::QuadObj,
                "ENDATA" => Section::End,
                other => return Err(fail(format!("unsupported section {other}"))),
            };
            if section == Section::End {
                break;
            }
            continue;
        }
        tokens.clear();
        match format {
            Format::Free => tokens.extend(line.split_whitespace()),
            Format::Fixed => fixed_fields(line, &mut tokens).map_err(fail)?,
        }
        builder.line(section, format, &tokens).map_err(fail)?;
    }
    if section != Section::End {
        return Err(Error::Parse {
            line: 0,
            message: "missing ENDATA".into(),
        });
    }
    builder.finish()
}

/// Split a fixed-format line into its six fields, trimmed. Empty fields are
/// kept so field positions stay meaningful.
fn fixed_fields<'a>(line: &'a str, fields: &mut Vec<&'a str>) -> Result<(), String> {
    const SPANS: [(usize, usize); 6] = [
        (1, 3),
        (4, 12),
        (14, 22),
        (24, 36),
        (39, 47),
        (49, usize::MAX),
    ];
    for (start, end) in SPANS {
        let end = end.min(line.len());
        let field = if start >= end {
            ""
        } else {
            line.get(start..end)
                .ok_or_else(|| "fixed-format field splits a character".to_string())?
        };
        fields.push(field.trim());
    }
    Ok(())
}

#[derive(Default)]
struct Builder<'a> {
    rows: HashMap<&'a str, Row>,
    row_names: Vec<&'a str>,
    kinds: Vec<RowKind>,
    has_objective: bool,
    columns: HashMap<&'a str, usize>,
    column_names: Vec<&'a str>,
    c: Vec<f64>,
    objective_set: Vec<bool>,
    a: Vec<(usize, usize, f64)>,
    rhs: Vec<f64>,
    ranges: Vec<Option<f64>>,
    constant: f64,
    bounds: Vec<Bounds>,
    lower_set: Vec<bool>,
    q: Vec<(usize, usize, f64)>,
}

impl<'a> Builder<'a> {
    fn line(&mut self, section: Section, format: Format, t: &[&'a str]) -> Result<(), String> {
        match (section, format) {
            (Section::Rows, Format::Free) => match t {
                [kind, name] => self.row(kind, name),
                _ => Err(shape("ROWS", t)),
            },
            (Section::Rows, Format::Fixed) => self.row(t[0], t[1]),
            (Section::Columns, Format::Free) => match t {
                [_, "'MARKER'", ..] => Ok(()),
                [column, row, value] => self.entry(column, &[(row, value)]),
                [column, r1, v1, r2, v2] => self.entry(column, &[(r1, v1), (r2, v2)]),
                _ => Err(shape("COLUMNS", t)),
            },
            (Section::Columns, Format::Fixed) => {
                if t[2] == "'MARKER'" {
                    Ok(())
                } else if t[4].is_empty() {
                    self.entry(t[1], &[(t[2], t[3])])
                } else {
                    self.entry(t[1], &[(t[2], t[3]), (t[4], t[5])])
                }
            }
            (Section::Rhs | Section::Ranges, _) => {
                let pairs: &[&str] = match format {
                    Format::Free if t.len() % 2 == 1 => &t[1..],
                    Format::Free => t,
                    Format::Fixed if t[4].is_empty() => &t[2..4],
                    Format::Fixed => &t[2..6],
                };
                if !matches!(pairs.len(), 2 | 4) {
                    return Err(shape("RHS or RANGES", t));
                }
                for pair in pairs.chunks(2) {
                    if section == Section::Rhs {
                        self.rhs(pair[0], pair[1])?;
                    } else {
                        self.range(pair[0], pair[1])?;
                    }
                }
                Ok(())
            }
            (Section::Bounds, Format::Free) => {
                let valued = match t.first() {
                    Some(&("UP" | "LO" | "FX")) => true,
                    Some(&("FR" | "MI" | "PL")) => false,
                    Some(other) => return Err(format!("unsupported bound type {other:?}")),
                    None => return Err(shape("BOUNDS", t)),
                };
                match (valued, t) {
                    (true, [kind, _, column, value] | [kind, column, value]) => {
                        self.bound(kind, column, Some(value))
                    }
                    (false, [kind, _, column] | [kind, column]) => self.bound(kind, column, None),
                    _ => Err(shape("BOUNDS", t)),
                }
            }
            (Section::Bounds, Format::Fixed) => {
                let value = (!t[3].is_empty()).then_some(t[3]);
                self.bound(t[0], t[2], value)
            }
            (Section::QuadObj, Format::Free) => match t {
                [first, second, value] => self.quadratic(first, second, value),
                _ => Err(shape("QUADOBJ", t)),
            },
            (Section::QuadObj, Format::Fixed) => self.quadratic(t[1], t[2], t[3]),
            (Section::None | Section::End, _) => Err("data line outside a section".into()),
        }
    }

    fn row(&mut self, kind: &str, name: &'a str) -> Result<(), String> {
        if name.is_empty() {
            return Err("row without a name".into());
        }
        let row = match kind {
            "N" if self.has_objective => Row::Free,
            "N" => {
                self.has_objective = true;
                Row::Objective
            }
            "E" | "L" | "G" => {
                self.kinds.push(match kind {
                    "E" => RowKind::Equal,
                    "L" => RowKind::Less,
                    _ => RowKind::Greater,
                });
                self.row_names.push(name);
                self.rhs.push(0.0);
                self.ranges.push(None);
                Row::Constraint(self.kinds.len() - 1)
            }
            other => return Err(format!("unknown row type {other:?}")),
        };
        if self.rows.insert(name, row).is_some() {
            return Err(format!("row {name} declared twice"));
        }
        Ok(())
    }

    fn find_row(&self, name: &str) -> Result<Row, String> {
        self.rows
            .get(name)
            .copied()
            .ok_or_else(|| format!("unknown row {name:?}"))
    }

    fn find_column(&self, name: &str) -> Result<usize, String> {
        self.columns
            .get(name)
            .copied()
            .ok_or_else(|| format!("unknown column {name:?}"))
    }

    fn entry(&mut self, column: &'a str, pairs: &[(&str, &str)]) -> Result<(), String> {
        if column.is_empty() {
            return Err("column without a name".into());
        }
        let n = self.columns.len();
        let j = *self.columns.entry(column).or_insert(n);
        if j == n {
            self.column_names.push(column);
            self.c.push(0.0);
            self.objective_set.push(false);
            self.bounds.push(Bounds {
                lower: 0.0,
                upper: f64::INFINITY,
            });
            self.lower_set.push(false);
        }
        for &(row, value) in pairs {
            let value = number(value)?;
            match self.find_row(row)? {
                Row::Objective => {
                    if std::mem::replace(&mut self.objective_set[j], true) {
                        return Err(format!("objective entry for column {column} given twice"));
                    }
                    self.c[j] = value;
                }
                Row::Free => {}
                Row::Constraint(i) => self.a.push((i, j, value)),
            }
        }
        Ok(())
    }

    fn rhs(&mut self, row: &str, value: &str) -> Result<(), String> {
        let value = number(value)?;
        match self.find_row(row)? {
            Row::Objective => self.constant = -value,
            Row::Free => {}
            Row::Constraint(i) => self.rhs[i] = value,
        }
        Ok(())
    }

    fn range(&mut self, row: &str, value: &str) -> Result<(), String> {
        let value = number(value)?;
        if let Row::Constraint(i) = self.find_row(row)? {
            self.ranges[i] = Some(value);
        }
        Ok(())
    }

    fn bound(&mut self, kind: &str, column: &str, value: Option<&str>) -> Result<(), String> {
        let j = self.find_column(column)?;
        let value = value.map(number).transpose()?;
        let b = &mut self.bounds[j];
        match (kind, value) {
            ("UP", Some(v)) => {
                if v < 0.0 && b.lower == 0.0 && !self.lower_set[j] {
                    b.lower = f64::NEG_INFINITY;
                }
                b.upper = v;
            }
            ("LO", Some(v)) => {
                b.lower = v;
                self.lower_set[j] = true;
            }
            ("FX", Some(v)) => {
                *b = Bounds::fixed(v);
                self.lower_set[j] = true;
            }
            ("FR", None) => {
                *b = Bounds::FREE;
                self.lower_set[j] = true;
            }
            ("MI", None) => {
                b.lower = f64::NEG_INFINITY;
                self.lower_set[j] = true;
            }
            ("PL", None) => b.upper = f64::INFINITY,
            ("UP" | "LO" | "FX", None) => return Err(format!("{kind} bound without a value")),
            ("FR" | "MI" | "PL", Some(_)) => return Err(format!("{kind} bound with a value")),
            (other, _) => return Err(format!("unsupported bound type {other:?}")),
        }
        Ok(())
    }

    fn quadratic(&mut self, first: &str, second: &str, value: &str) -> Result<(), String> {
        let (i, j) = (self.find_column(first)?, self.find_column(second)?);
        self.q.push((i.min(j), i.max(j), number(value)?));
        Ok(())
    }

    fn finish(self) -> Result<Problem, Error> {
        if !self.has_objective {
            return Err(Error::Parse {
                line: 0,
                message: "no objective row".into(),
            });
        }
        let m = self.kinds.len();
        let n = self.c.len();
        let rows = (0..m)
            .map(|i| {
                let rhs = self.rhs[i];
                let r = self.ranges[i].unwrap_or(f64::NAN);
                let (lower, upper) = match (self.kinds[i], self.ranges[i].is_some()) {
                    (RowKind::Equal, false) => (rhs, rhs),
                    (RowKind::Less, false) => (f64::NEG_INFINITY, rhs),
                    (RowKind::Greater, false) => (rhs, f64::INFINITY),
                    (RowKind::Less, true) => (rhs - r.abs(), rhs),
                    (RowKind::Greater, true) => (rhs, rhs + r.abs()),
                    (RowKind::Equal, true) if r >= 0.0 => (rhs, rhs + r),
                    (RowKind::Equal, true) => (rhs + r, rhs),
                };
                Constraint::Linear(Bounds { lower, upper })
            })
            .collect();
        let duplicate = |(row, column): (usize, usize), what: &str, rows: &[&str]| Error::Parse {
            line: 0,
            message: format!(
                "{what} entry ({}, {}) given twice",
                rows[row], self.column_names[column]
            ),
        };
        let a = csc(m, n, self.a).map_err(|e| match e {
            Built::Duplicate(at) => duplicate(at, "matrix", &self.row_names),
            Built::Matrix(e) => Error::Matrix(e),
        })?;
        let p = if self.q.is_empty() {
            None
        } else {
            Some(csc(n, n, self.q).map_err(|e| match e {
                Built::Duplicate(at) => duplicate(at, "QUADOBJ", &self.column_names),
                Built::Matrix(e) => Error::Matrix(e),
            })?)
        };
        Ok(Problem {
            p,
            c: self.c,
            c0: self.constant,
            a,
            rows,
            variable_bounds: self.bounds,
            cones: vec![],
        })
    }
}

enum Built {
    Duplicate((usize, usize)),
    Matrix(MatrixError),
}

/// Compressed columns from (row, column, value) triplets, dropping zeros
/// and rejecting repeated positions.
fn csc(
    rows: usize,
    columns: usize,
    mut triplets: Vec<(usize, usize, f64)>,
) -> Result<CscMatrix, Built> {
    triplets.sort_unstable_by_key(|&(i, j, _)| (j, i));
    if let Some(w) = triplets
        .windows(2)
        .find(|w| (w[0].0, w[0].1) == (w[1].0, w[1].1))
    {
        return Err(Built::Duplicate((w[0].0, w[0].1)));
    }
    let mut pointers = vec![0; columns + 1];
    let mut indices = Vec::with_capacity(triplets.len());
    let mut values = Vec::with_capacity(triplets.len());
    for (i, j, v) in triplets {
        if v != 0.0 {
            pointers[j + 1] += 1;
            indices.push(i);
            values.push(v);
        }
    }
    for j in 0..columns {
        pointers[j + 1] += pointers[j];
    }
    CscMatrix::new(rows, columns, pointers, indices, values).map_err(Built::Matrix)
}

fn number(token: &str) -> Result<f64, String> {
    token
        .parse::<f64>()
        .ok()
        .filter(|v| !v.is_nan())
        .ok_or_else(|| format!("invalid number {token:?}"))
}

fn shape(section: &str, tokens: &[&str]) -> String {
    format!("unexpected {section} line with {} fields", tokens.len())
}

#[cfg(test)]
mod tests {
    use super::*;

    const INF: f64 = f64::INFINITY;

    fn row_bounds(problem: &Problem) -> Vec<(f64, f64)> {
        problem
            .rows
            .iter()
            .map(|row| match row {
                Constraint::Linear(b) => (b.lower, b.upper),
                Constraint::Cone { .. } => unreachable!(),
            })
            .collect()
    }

    fn column(matrix: &CscMatrix, j: usize) -> Vec<(usize, f64)> {
        matrix.as_ref().column(j).collect()
    }

    fn bounds(problem: &Problem) -> Vec<(f64, f64)> {
        problem
            .variable_bounds
            .iter()
            .map(|b| (b.lower, b.upper))
            .collect()
    }

    #[test]
    fn free_format_linear_program() {
        let text = "\
* comment
NAME          TOY extra words
ROWS
 N  COST
 E  BAL
 L  CAP
 G  DEM
 N  SPARE
 E  NEG
COLUMNS
    MARKER    'MARKER'  'INTORG'
    X         COST      1.5   BAL   1.
    X         CAP       2.    SPARE 9.

    MARKER    'MARKER'  'INTEND'
    Y         BAL       -1    DEM   1e0
    Y         NEG       .5
    Z         COST      -2.
RHS
    RHS       COST      4.    BAL   3.
    RHS       CAP       10.   DEM   1.
    RHS       NEG       6.
RANGES
    RNG       BAL       2.    CAP   -4.
    RNG       DEM       5.    NEG   -1.
    RNG       SPARE     7.
BOUNDS
 UP BND       X         8.
 MI BND       Y
 UP BND       Y         3.
 FR BND       Z
ENDATA
";
        let problem = parse(text).unwrap();
        assert_eq!(problem.c, [1.5, 0.0, -2.0]);
        assert_eq!(problem.c0, -4.0);
        assert!(problem.p.is_none());
        assert_eq!(
            row_bounds(&problem),
            [(3.0, 5.0), (6.0, 10.0), (1.0, 6.0), (5.0, 6.0)]
        );
        assert_eq!(column(&problem.a, 0), [(0, 1.0), (1, 2.0)]);
        assert_eq!(column(&problem.a, 1), [(0, -1.0), (2, 1.0), (3, 0.5)]);
        assert_eq!(column(&problem.a, 2), []);
        assert_eq!(bounds(&problem), [(0.0, 8.0), (-INF, 3.0), (-INF, INF)]);
    }

    #[test]
    fn quadratic_objective_is_stored_as_an_upper_triangle() {
        let text = "\
NAME QP
ROWS
 N obj
 L r
COLUMNS
 x r 1
 y r 1
 z r 0
RHS
 rhs r 1
QUADOBJ
 x x 2
 y x -1
 y y 4
 z x 0
ENDATA
";
        let problem = parse(text).unwrap();
        let p = problem.p.as_ref().unwrap();
        assert_eq!((p.rows(), p.columns()), (3, 3));
        assert_eq!(column(p, 0), [(0, 2.0)]);
        assert_eq!(column(p, 1), [(0, -1.0), (1, 4.0)]);
        assert_eq!(column(p, 2), []);
        // An explicit zero coefficient is dropped from A.
        assert_eq!(column(&problem.a, 2), []);
    }

    #[test]
    fn omitted_set_names_are_accepted() {
        let text = "\
NAME
ROWS
 N obj
 L a
 L b
COLUMNS
 x a 1 b 1
 y a 1
RHS
 a 4 b 5
 obj 1
BOUNDS
 UP x 2
 LO y -1
 PL y
ENDATA
";
        let problem = parse(text).unwrap();
        assert_eq!(row_bounds(&problem), [(-INF, 4.0), (-INF, 5.0)]);
        assert_eq!(problem.c0, -1.0);
        assert_eq!(bounds(&problem), [(0.0, 2.0), (-1.0, INF)]);
    }

    #[test]
    fn negative_upper_bound_without_lower_frees_the_lower_side() {
        let text = "\
NAME
ROWS
 N obj
COLUMNS
 x obj 1
 y obj 1
BOUNDS
 UP b x -1
 LO b y 0
 UP b y -1
ENDATA
";
        let problem = parse(text).unwrap();
        assert_eq!(bounds(&problem), [(-INF, -1.0), (0.0, -1.0)]);
    }

    /// Lay out one fixed-format data line at the standard field positions.
    fn fixed(fields: [&str; 6]) -> String {
        let [f1, f2, f3, f4, f5, f6] = fields;
        format!(" {f1:<2} {f2:<8}  {f3:<8}  {f4:>12}   {f5:<8}  {f6:>12}")
            .trim_end()
            .to_string()
    }

    #[test]
    fn fixed_format_names_may_contain_spaces() {
        let lines = [
            "NAME          FIXED  (VARIANT)".to_string(),
            "ROWS".into(),
            fixed(["N", "COST", "", "", "", ""]),
            fixed(["E", "ROW 1", "", "", "", ""]),
            fixed(["L", "ROW 2", "", "", "", ""]),
            "COLUMNS".into(),
            fixed(["", "COL A", "COST", "1.", "ROW 1", "2."]),
            fixed(["", "COL A", "ROW 2", "3.", "", ""]),
            fixed(["", "COL B", "ROW 1", "-1.", "", ""]),
            "RHS".into(),
            fixed(["", "RHS 1", "ROW 1", "7.", "ROW 2", "8."]),
            "RANGES".into(),
            fixed(["", "RNG 1", "ROW 2", "3.", "", ""]),
            "BOUNDS".into(),
            fixed(["UP", "BND 1", "COL A", "5.", "", ""]),
            fixed(["FR", "BND 1", "COL B", "", "", ""]),
            "ENDATA".into(),
        ];
        let problem = parse(&lines.join("\n")).unwrap();
        assert_eq!(problem.c, [1.0, 0.0]);
        assert_eq!(row_bounds(&problem), [(7.0, 7.0), (5.0, 8.0)]);
        assert_eq!(column(&problem.a, 0), [(0, 2.0), (1, 3.0)]);
        assert_eq!(column(&problem.a, 1), [(0, -1.0)]);
        assert_eq!(bounds(&problem), [(0.0, 5.0), (-INF, INF)]);
    }

    fn error(text: &str) -> String {
        parse(text).unwrap_err().to_string()
    }

    #[test]
    fn errors_name_the_offending_line() {
        let head = "NAME\nROWS\n N obj\n L r\nCOLUMNS\n";
        assert_eq!(
            error(&format!("{head} x q 1\nENDATA\n")),
            "line 6: unknown row \"q\""
        );
        assert_eq!(
            error(&format!("{head} x r one\nENDATA\n")),
            "line 6: invalid number \"one\""
        );
        assert_eq!(
            error(&format!("{head} x r 1\nOBJSENSE\n MAX\nENDATA\n")),
            "line 7: unsupported section OBJSENSE"
        );
        assert_eq!(
            error(&format!("{head} x r 1\nBOUNDS\n BV b x\nENDATA\n")),
            "line 8: unsupported bound type \"BV\""
        );
        assert_eq!(error(&format!("{head} x r 1\n")), "missing ENDATA");
        assert_eq!(
            error(&format!("{head} x r 1\n x r 2\nENDATA\n")),
            "matrix entry (r, x) given twice"
        );
        assert_eq!(
            error(&format!("{head} x r 1\nQUADOBJ\n x x 1\n x x 1\nENDATA\n")),
            "QUADOBJ entry (x, x) given twice"
        );
        assert_eq!(
            error("NAME\nROWS\n L r\nCOLUMNS\n x r 1\nENDATA\n"),
            "no objective row"
        );
    }
}
