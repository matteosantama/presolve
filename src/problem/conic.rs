//! Optional conversion of ranged rows and bounds to Ax+s=b.
use crate::{
    matrix::CscMatrix,
    postsolve::{PrimalCertificate, Solution},
    problem::{Bounds, Cone, Constraint, Problem, row_indices},
};

/// CSC conic form. P is still upper triangular; the objective constant is kept.
#[derive(Clone, Debug)]
pub struct ConicData {
    pub p: Option<CscMatrix>,
    pub c: Vec<f64>,
    pub c0: f64,
    pub a: CscMatrix,
    pub b: Vec<f64>,
    pub cones: Vec<Cone>,
}
#[derive(Debug)]
pub struct ConicExport {
    pub problem: ConicData,
    pub map: ConicMap,
}
#[derive(Clone, Copy, Debug)]
enum Output {
    Row {
        index: usize,
        scale: f64,
        equality: bool,
    },
    Bound {
        column: usize,
        scale: f64,
        equality: bool,
    },
    Cone(usize),
}
/// Converts conic-form multipliers to native linear, bound and conic
/// multipliers. Compose this with `Postsolve` to recover original coordinates.
#[derive(Debug)]
pub struct ConicMap {
    outputs: Vec<Output>,
    variables: usize,
    linear: usize,
    conic: usize,
}
impl ConicMap {
    pub fn variable_count(&self) -> usize {
        self.variables
    }
    pub fn linear_row_count(&self) -> usize {
        self.linear
    }
    pub fn conic_row_count(&self) -> usize {
        self.conic
    }
    pub fn output_dimension(&self) -> usize {
        self.outputs.len()
    }

    pub fn conic_offset(&self) -> usize {
        self.outputs.len() - self.conic
    }
    pub fn dual(&self, y: &[f64]) -> PrimalCertificate {
        let mut out = PrimalCertificate {
            y: vec![0.; self.linear],
            z: vec![0.; self.variables],
            conic_dual: vec![0.; self.conic],
        };
        self.dual_into(y, &mut out.y, &mut out.z, &mut out.conic_dual);
        out
    }
    /// Allocation-free conversion, also used for infeasibility multipliers.
    pub fn dual_into(&self, y: &[f64], linear: &mut [f64], bounds: &mut [f64], conic: &mut [f64]) {
        linear.fill(0.);
        bounds.fill(0.);
        conic.fill(0.);
        for (&out, &v) in self.outputs.iter().zip(y) {
            match out {
                Output::Row { index, scale, .. } => linear[index] += scale * v,
                Output::Bound { column, scale, .. } => bounds[column] += scale * v,
                Output::Cone(i) => conic[i] = v,
            }
        }
    }
    /// Map native warm-start multipliers to individual constraint sides.
    pub fn warm_dual(&self, point: &Solution) -> Vec<f64> {
        self.outputs
            .iter()
            .map(|out| match *out {
                Output::Row {
                    index,
                    scale,
                    equality,
                } => {
                    let y = scale * point.y[index];
                    if equality { y } else { y.max(0.) }
                }
                Output::Bound {
                    column,
                    scale,
                    equality,
                } => {
                    let y = scale * point.z[column];
                    if equality { y } else { y.max(0.) }
                }
                Output::Cone(i) => point.conic_dual[i],
            })
            .collect()
    }
}
fn scalar(cones: &mut Vec<Cone>, equality: bool) {
    match (cones.last_mut(), equality) {
        (Some(Cone::Zero(n)), true) | (Some(Cone::Nonnegative(n)), false) => *n += 1,
        (_, true) => cones.push(Cone::Zero(1)),
        (_, false) => cones.push(Cone::Nonnegative(1)),
    }
}
impl Problem {
    /// Expand finite ranged-row and variable-bound sides directly into CSC.
    /// The returned map is necessary for dual and warm-start recovery.
    pub fn into_conic(self) -> ConicExport {
        let n = self.variable_count();
        let (linear_rows, conic_rows) = row_indices(&self.rows);
        let bounds_of = |i: usize| match self.rows[linear_rows[i]] {
            Constraint::Linear(b) => b,
            Constraint::Cone { .. } => unreachable!(),
        };
        let rhs_of = |i: usize| match self.rows[conic_rows[i]] {
            Constraint::Cone { rhs, .. } => rhs,
            Constraint::Linear(_) => unreachable!(),
        };
        let mut outputs = Vec::new();
        let mut b = Vec::new();
        let mut cones = Vec::new();
        let mut append = |bounds: Bounds, row: Option<usize>, column: usize| {
            let mut side = |rhs, scale, equality| {
                scalar(&mut cones, equality);
                b.push(rhs);
                outputs.push(match row {
                    Some(index) => Output::Row {
                        index,
                        scale,
                        equality,
                    },
                    None => Output::Bound {
                        column,
                        scale,
                        equality,
                    },
                });
            };
            if bounds.equality() {
                side(bounds.lower, -1., true);
            } else {
                if bounds.lower.is_finite() {
                    side(-bounds.lower, 1., false);
                }
                if bounds.upper.is_finite() {
                    side(bounds.upper, -1., false);
                }
            }
        };
        for i in 0..linear_rows.len() {
            append(bounds_of(i), Some(i), 0);
        }
        for j in 0..n {
            append(self.variable_bounds(j), None, j);
        }
        cones.extend_from_slice(&self.cones);
        for i in 0..conic_rows.len() {
            outputs.push(Output::Cone(i));
            b.push(rhs_of(i));
        }
        let map = ConicMap {
            outputs,
            variables: n,
            linear: linear_rows.len(),
            conic: conic_rows.len(),
        };
        let outputs = &map.outputs;
        let identity = outputs.len() == self.row_count()
            && outputs.iter().enumerate().all(|(i, out)| match *out {
                Output::Row { index, scale, .. } => linear_rows[index] == i && scale == -1.,
                Output::Cone(index) => conic_rows[index] == i,
                Output::Bound { .. } => false,
            });
        if identity {
            let Problem { p, c, c0, a, .. } = self;
            return ConicExport {
                problem: ConicData {
                    p,
                    c,
                    c0,
                    a,
                    b,
                    cones,
                },
                map,
            };
        }
        let a = {
            let matrix = &self.a;
            // Scatter each original column through at most two row sides.
            // Sorting is local to a column, never a scan over every (row,col).
            let mut row_map = vec![[None; 2]; self.row_count()];
            let mut bound_map = vec![[None; 2]; n];
            for (at, out) in outputs.iter().enumerate() {
                let (slots, scale) = match *out {
                    Output::Row { index, scale, .. } => (&mut row_map[linear_rows[index]], -scale),
                    Output::Cone(i) => (&mut row_map[conic_rows[i]], 1.),
                    Output::Bound { column, scale, .. } => (&mut bound_map[column], -scale),
                };
                let slot = usize::from(slots[0].is_some());
                slots[slot] = Some((at, scale));
            }
            let mut pointers = Vec::with_capacity(n + 1);
            let mut ri = Vec::new();
            let mut values = Vec::new();
            let mut column = Vec::new();
            pointers.push(0);
            for (j, bounds) in bound_map.iter().enumerate() {
                column.clear();
                for (i, v) in matrix.as_ref().column(j) {
                    for &(row, scale) in row_map[i].iter().flatten() {
                        column.push((row, scale * v));
                    }
                }
                column.extend(bounds.iter().flatten().copied());
                column.sort_unstable_by_key(|&(i, _)| i);
                for &(i, v) in &column {
                    ri.push(i);
                    values.push(v);
                }
                pointers.push(values.len());
            }
            CscMatrix::from_parts(outputs.len(), n, pointers, ri, values)
        };
        let Problem { p, c, c0, .. } = self;
        ConicExport {
            problem: ConicData {
                p,
                c,
                c0,
                a,
                b,
                cones,
            },
            map,
        }
    }
}
