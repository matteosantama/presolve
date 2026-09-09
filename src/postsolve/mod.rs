//! Public solution recovery and mappings between original, working, and reduced coordinates.
mod solution;
pub(crate) mod tape;

pub use solution::{CertificateRef, PrimalCertificate, Solution, SolutionMut, SolutionRef};

use std::sync::Arc;
use tape::{Point, Recovery, RecoveryTape};

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
pub(crate) fn original_point(
    mut p: Point,
    n: usize,
    linear: &[usize],
    conic: &[usize],
    slack: Vec<f64>,
) -> Solution {
    p.x.truncate(n);
    p.z.truncate(n);
    let conic_dual = conic.iter().map(|&i| -p.y[i]).collect();
    let y = if conic.is_empty() {
        p.y.truncate(linear.len());
        p.y
    } else {
        linear.iter().map(|&i| p.y[i]).collect()
    };
    Solution {
        x: p.x,
        z: p.z,
        y,
        conic_dual,
        conic_slack: slack,
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
    fn scatter(&self, p: &mut Point, x: &[f64], y: &[f64], z: &[f64], w: &[f64]) {
        for (&j, &v) in self.coordinates.compact_to_stable_columns.iter().zip(x) {
            p.x[j] = v;
        }
        for (&j, &v) in self.coordinates.compact_to_stable_columns.iter().zip(z) {
            p.z[j] = v;
        }
        for (&i, &v) in self.coordinates.compact_to_stable_linear_rows.iter().zip(y) {
            p.y[i] = v;
        }
        for (&i, &v) in self.coordinates.compact_to_stable_conic_rows.iter().zip(w) {
            p.y[i] = -v;
        }
    }
    fn lift(&self, x: &[f64], y: &[f64], z: &[f64], w: &[f64]) -> Point {
        let mut p = Point::zeros(self.total_columns, self.total_rows);
        self.scatter(&mut p, x, y, z, w);
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
    pub fn surviving_conic_coordinates(&self) -> &[(usize, usize)] {
        &self.direct_slacks
    }
    pub fn reduced_conic_dimension(&self) -> usize {
        self.coordinates.compact_to_stable_conic_rows.len()
    }
    /// Recover into reusable caller buffers without allocating. Dimensions and
    /// workspace compatibility are the caller's responsibility.
    pub fn recover_into(&self, point: SolutionRef<'_>, out: SolutionMut<'_>, work: &mut Workspace) {
        let p = &mut work.point;
        p.x.fill(0.);
        p.y.fill(0.);
        p.z.fill(0.);
        self.scatter(p, point.x, point.y, point.z, point.conic_dual);
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
        out.x.copy_from_slice(&p.x[..self.original_columns]);
        out.z.copy_from_slice(&p.z[..self.original_columns]);
        for (&i, v) in self.input_linear.iter().zip(out.y) {
            *v = p.y[i];
        }
        for (&i, v) in self.input_conic.iter().zip(out.conic_dual) {
            *v = -p.y[i];
        }
        for (&i, v) in self.input_conic.iter().zip(out.conic_slack) {
            *v = work.slacks[i];
        }
    }
    /// Allocating convenience wrapper around `recover_into`.
    pub fn recover_solution(&self, point: SolutionRef<'_>) -> Solution {
        let mut out = self.solution();
        self.recover_into(point, out.as_mut(), &mut self.workspace());
        out
    }
    pub fn recover_primal_ray(&self, x: &[f64]) -> Vec<f64> {
        let mut p = self.lift(x, &[], &[], &[]);
        self.tape.recover(&mut p, Recovery::DualInfeasibility);
        p.x.truncate(self.original_columns);
        p.x
    }
    pub fn recover_primal_certificate(&self, c: CertificateRef<'_>) -> PrimalCertificate {
        let mut p = self.lift(&[], c.y, c.z, c.conic_dual);
        self.tape.recover(&mut p, Recovery::PrimalInfeasibility);
        let p = original_point(
            p,
            self.original_columns,
            &self.input_linear,
            &self.input_conic,
            Vec::new(),
        );
        PrimalCertificate {
            y: p.y,
            z: p.z,
            conic_dual: p.conic_dual,
        }
    }
    pub fn reduce_warm_start(&self, point: SolutionRef<'_>) -> Solution {
        let mut p = Point::zeros(self.total_columns, self.total_rows);
        p.x[..self.original_columns].copy_from_slice(point.x);
        p.z[..self.original_columns].copy_from_slice(point.z);
        for (&i, &v) in self.input_linear.iter().zip(point.y) {
            p.y[i] = v;
        }
        for (&i, &v) in self.input_conic.iter().zip(point.conic_dual) {
            p.y[i] = -v;
        }
        self.tape.reduce_point(&mut p);
        Solution {
            x: self
                .coordinates
                .compact_to_stable_columns
                .iter()
                .map(|&j| p.x[j])
                .collect(),
            z: self
                .coordinates
                .compact_to_stable_columns
                .iter()
                .map(|&j| p.z[j])
                .collect(),
            y: self
                .coordinates
                .compact_to_stable_linear_rows
                .iter()
                .map(|&i| p.y[i])
                .collect(),
            conic_dual: self
                .coordinates
                .compact_to_stable_conic_rows
                .iter()
                .map(|&i| -p.y[i])
                .collect(),
            conic_slack: self
                .direct_slacks
                .iter()
                .map(|&(original, _)| point.conic_slack[original])
                .collect(),
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
