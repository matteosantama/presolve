//! Public solution recovery and mappings between original, working, and reduced coordinates.
mod solution;

pub use solution::{CertificateRef, PrimalCertificate, Solution, SolutionMut, SolutionRef};

use crate::model::tape::{Point, Recovery, RecoveryTape};
use std::sync::Arc;

/// Working indices remain stable through rules; auxiliary variables are appended.
/// These maps lift compact reduced coordinates into that working index space.
/// Original variables occupy its initial prefix, while input row maps select the
/// caller's original linear rows and cone coordinates.
#[derive(Debug)]
pub(crate) struct Coordinates {
    pub(crate) compact_to_stable_columns: Arc<Vec<usize>>,
    pub(crate) compact_to_stable_linear_rows: Vec<usize>,
    pub(crate) compact_to_stable_conic_rows: Vec<usize>,
}

/// Reusable scratch space, allocated once for a particular postsolve map.
pub struct Workspace {
    point: Point,
    slacks: Vec<f64>,
}

/// Original-coordinate solution borrowed from a recovery workspace.
///
/// Primal variables and bound multipliers are contiguous slices. Linear and
/// conic rows are gathered lazily in their original block order, with the same
/// signs as [`Solution`]. No result buffers are allocated or copied. The view
/// must be released before the workspace can be reused.
pub struct RecoveredSolution<'a> {
    point: &'a Point,
    slacks: &'a [f64],
    original: OriginalMap<'a>,
}

impl<'a> RecoveredSolution<'a> {
    pub fn x(&self) -> &'a [f64] {
        &self.point.x[..self.original.columns]
    }

    pub fn z(&self) -> &'a [f64] {
        &self.point.z[..self.original.columns]
    }

    pub fn y(&self) -> impl ExactSizeIterator<Item = f64> + 'a {
        let y = &self.point.y;
        self.original.linear.iter().map(move |&i| y[i])
    }

    pub fn conic_dual(&self) -> impl ExactSizeIterator<Item = f64> + 'a {
        let y = &self.point.y;
        self.original.conic.iter().map(move |&i| -y[i])
    }

    pub fn conic_slack(&self) -> impl ExactSizeIterator<Item = f64> + 'a {
        let slacks = self.slacks;
        self.original.conic.iter().map(move |&i| slacks[i])
    }
}

#[derive(Debug)]
pub struct Postsolve {
    pub(crate) tape: RecoveryTape,
    pub(crate) coordinates: Coordinates,
    pub(crate) input_linear: Vec<usize>,
    pub(crate) input_conic: Vec<usize>,
    pub(crate) total_columns: usize,
    pub(crate) total_rows: usize,
    pub(crate) original_columns: usize,
    pub(crate) direct_slacks: Vec<(usize, usize)>,
}
/// The caller's original layout: a column prefix of the working point plus
/// the input positions of its linear rows and cone coordinates. A conic
/// multiplier is stored negated in the working point.
pub(crate) struct OriginalMap<'a> {
    pub columns: usize,
    pub linear: &'a [usize],
    pub conic: &'a [usize],
}
impl OriginalMap<'_> {
    /// Original coordinates of a working point, taking its buffers.
    pub(crate) fn gather(&self, mut p: Point, conic_slack: Vec<f64>) -> Solution {
        p.x.truncate(self.columns);
        p.z.truncate(self.columns);
        let conic_dual = self.conic.iter().map(|&i| -p.y[i]).collect();
        let y = if self.conic.is_empty() {
            p.y.truncate(self.linear.len());
            p.y
        } else {
            self.linear.iter().map(|&i| p.y[i]).collect()
        };
        Solution {
            x: p.x,
            z: p.z,
            y,
            conic_dual,
            conic_slack,
        }
    }
    /// Original multipliers of a working point: linear rows, bounds, cones.
    pub(crate) fn gather_dual(&self, p: &Point) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        (
            self.linear.iter().map(|&i| p.y[i]).collect(),
            p.z[..self.columns].to_vec(),
            self.conic.iter().map(|&i| -p.y[i]).collect(),
        )
    }
    fn gather_into(&self, p: &Point, slacks: &[f64], out: SolutionMut<'_>) {
        out.x.copy_from_slice(&p.x[..self.columns]);
        out.z.copy_from_slice(&p.z[..self.columns]);
        for (&i, v) in self.linear.iter().zip(out.y) {
            *v = p.y[i];
        }
        for (&i, v) in self.conic.iter().zip(out.conic_dual) {
            *v = -p.y[i];
        }
        for (&i, v) in self.conic.iter().zip(out.conic_slack) {
            *v = slacks[i];
        }
    }
    fn scatter(&self, p: &mut Point, s: SolutionRef<'_>) {
        p.x[..self.columns].copy_from_slice(s.x);
        p.z[..self.columns].copy_from_slice(s.z);
        for (&i, &v) in self.linear.iter().zip(s.y) {
            p.y[i] = v;
        }
        for (&i, &v) in self.conic.iter().zip(s.conic_dual) {
            p.y[i] = -v;
        }
    }
}
impl Coordinates {
    /// Compact reduced coordinates into the stable working point.
    fn scatter(&self, p: &mut Point, s: SolutionRef<'_>) {
        for (&j, &v) in self.compact_to_stable_columns.iter().zip(s.x) {
            p.x[j] = v;
        }
        for (&j, &v) in self.compact_to_stable_columns.iter().zip(s.z) {
            p.z[j] = v;
        }
        for (&i, &v) in self.compact_to_stable_linear_rows.iter().zip(s.y) {
            p.y[i] = v;
        }
        for (&i, &v) in self.compact_to_stable_conic_rows.iter().zip(s.conic_dual) {
            p.y[i] = -v;
        }
    }
    /// Compact coordinates of a working point, as `(x, z, y, conic_dual)`.
    fn gather(&self, p: &Point) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
        let pick = |map: &[usize], from: &[f64]| map.iter().map(|&k| from[k]).collect();
        (
            pick(&self.compact_to_stable_columns, &p.x),
            pick(&self.compact_to_stable_columns, &p.z),
            pick(&self.compact_to_stable_linear_rows, &p.y),
            self.compact_to_stable_conic_rows
                .iter()
                .map(|&i| -p.y[i])
                .collect(),
        )
    }
}
impl Postsolve {
    /// Positions of the original linear rows in the input problem's row ordering.
    pub fn original_linear_rows(&self) -> &[usize] {
        &self.input_linear
    }
    /// Positions of the original cone coordinates, in cone-block order.
    pub fn original_conic_rows(&self) -> &[usize] {
        &self.input_conic
    }
    fn original(&self) -> OriginalMap<'_> {
        OriginalMap {
            columns: self.original_columns,
            linear: &self.input_linear,
            conic: &self.input_conic,
        }
    }
    fn lift(&self, s: SolutionRef<'_>) -> Point {
        let mut p = Point::zeros(self.total_columns, self.total_rows);
        self.coordinates.scatter(&mut p, s);
        p
    }
    pub fn workspace(&self) -> Workspace {
        Workspace {
            point: Point::zeros(self.total_columns, self.total_rows),
            slacks: vec![
                0.;
                if self.input_conic.is_empty() {
                    0
                } else {
                    self.total_rows
                }
            ],
        }
    }
    pub fn solution(&self) -> Solution {
        Solution {
            x: vec![0.; self.original_columns],
            z: vec![0.; self.original_columns],
            y: vec![0.; self.input_linear.len()],
            conic_dual: vec![0.; self.input_conic.len()],
            conic_slack: vec![0.; self.input_conic.len()],
        }
    }
    /// Original/reduced coordinate pairs whose slacks survive without a
    /// coordinate transformation. Use `reduce_warm_start` for the full map.
    pub fn surviving_conic_coordinates(&self) -> &[(usize, usize)] {
        &self.direct_slacks
    }
    pub fn reduced_conic_dimension(&self) -> usize {
        self.coordinates.compact_to_stable_conic_rows.len()
    }
    /// Recover into reusable caller buffers without allocating. Dimensions and
    /// workspace compatibility are the caller's responsibility.
    pub fn recover_into(&self, point: SolutionRef<'_>, out: SolutionMut<'_>, work: &mut Workspace) {
        let recovered = self.recover_borrowed(point, work);
        recovered
            .original
            .gather_into(recovered.point, recovered.slacks, out);
    }

    /// Recover without copying the result out of the workspace. The workspace
    /// must have been created by this postsolve map; input dimensions must match
    /// the reduced problem, just as for [`Self::recover_into`].
    pub fn recover_borrowed<'a>(
        &'a self,
        point: SolutionRef<'_>,
        work: &'a mut Workspace,
    ) -> RecoveredSolution<'a> {
        let p = &mut work.point;
        p.x.fill(0.);
        p.y.fill(0.);
        p.z.fill(0.);
        self.coordinates.scatter(p, point);
        work.slacks.fill(0.);
        for (&i, &v) in self
            .coordinates
            .compact_to_stable_conic_rows
            .iter()
            .zip(point.conic_slack)
        {
            work.slacks[i] = v;
        }
        self.tape
            .recover_with_slacks(p, Recovery::Solution, &mut work.slacks);
        RecoveredSolution {
            point: p,
            slacks: &work.slacks,
            original: self.original(),
        }
    }
    /// Allocating convenience wrapper around `recover_into`.
    pub fn recover_solution(&self, point: SolutionRef<'_>) -> Solution {
        let mut out = self.solution();
        self.recover_into(point, out.as_mut(), &mut self.workspace());
        out
    }
    pub fn recover_primal_ray(&self, x: &[f64]) -> Vec<f64> {
        let mut p = self.lift(SolutionRef {
            x,
            y: &[],
            z: &[],
            conic_dual: &[],
            conic_slack: &[],
        });
        self.tape.recover(&mut p, Recovery::DualInfeasibility);
        p.x.truncate(self.original_columns);
        p.x
    }
    pub fn recover_primal_certificate(&self, c: CertificateRef<'_>) -> PrimalCertificate {
        let mut p = self.lift(SolutionRef {
            x: &[],
            y: c.y,
            z: c.z,
            conic_dual: c.conic_dual,
            conic_slack: &[],
        });
        self.tape.recover(&mut p, Recovery::PrimalInfeasibility);
        let (y, z, conic_dual) = self.original().gather_dual(&p);
        PrimalCertificate { y, z, conic_dual }
    }
    pub fn reduce_warm_start(&self, point: SolutionRef<'_>) -> Solution {
        let mut p = Point::zeros(self.total_columns, self.total_rows);
        self.original().scatter(&mut p, point);
        self.tape.reduce_point(&mut p);
        let (x, z, y, conic_dual) = self.coordinates.gather(&p);
        let conic_slack =
            if !self.input_conic.is_empty() && self.tape.transforms_conic_coordinates() {
                let mut slacks = vec![0.0; self.total_rows];
                for (&row, &value) in self.input_conic.iter().zip(point.conic_slack) {
                    slacks[row] = value;
                }
                self.tape.reduce_slacks(&mut slacks);
                self.coordinates
                    .compact_to_stable_conic_rows
                    .iter()
                    .map(|&row| slacks[row])
                    .collect()
            } else {
                self.direct_slacks
                    .iter()
                    .map(|&(original, _)| point.conic_slack[original])
                    .collect()
            };
        Solution {
            x,
            z,
            y,
            conic_dual,
            conic_slack,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn recovery_preserves_coordinate_order_and_dual_signs_with_reused_buffers() {
        let postsolve = Postsolve {
            tape: RecoveryTape::default(),
            coordinates: Coordinates {
                compact_to_stable_columns: Arc::new(vec![0, 2]),
                compact_to_stable_linear_rows: vec![2, 0],
                compact_to_stable_conic_rows: vec![3, 1],
            },
            input_linear: vec![0, 2],
            input_conic: vec![1, 3],
            total_columns: 3,
            total_rows: 4,
            original_columns: 3,
            direct_slacks: vec![(1, 0), (0, 1)],
        };
        let point = SolutionRef {
            x: &[10., 30.],
            y: &[3., 1.],
            z: &[4., 6.],
            conic_dual: &[7., 8.],
            conic_slack: &[9., 10.],
        };
        let mut work = postsolve.workspace();
        work.point.x.fill(99.);
        work.point.y.fill(99.);
        work.point.z.fill(99.);
        work.slacks.fill(99.);
        let mut out = postsolve.solution();
        postsolve.recover_into(point, out.as_mut(), &mut work);
        assert_eq!(out.x, [10., 0., 30.]);
        assert_eq!(out.y, [1., 3.]);
        assert_eq!(out.z, [4., 0., 6.]);
        assert_eq!(out.conic_dual, [8., 7.]);
        assert_eq!(out.conic_slack, [10., 9.]);
        let certificate = postsolve.recover_primal_certificate(CertificateRef {
            y: point.y,
            z: point.z,
            conic_dual: point.conic_dual,
        });
        assert_eq!(certificate.y, out.y);
        assert_eq!(certificate.z, out.z);
        assert_eq!(certificate.conic_dual, out.conic_dual);
        assert_eq!(postsolve.recover_primal_ray(point.x), out.x);
        let reduced = postsolve.reduce_warm_start(out.as_ref());
        assert_eq!(reduced.x, point.x);
        assert_eq!(reduced.y, point.y);
        assert_eq!(reduced.z, point.z);
        assert_eq!(reduced.conic_dual, point.conic_dual);
        assert_eq!(reduced.conic_slack, point.conic_slack);
    }
}
