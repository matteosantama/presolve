//! Standard cone geometry. All dimensions and parameters are caller-validated.

/// Cone coordinates use the following fixed conventions. PSD coordinates are
/// upper-triangular, column-major `svec`, with off-diagonals scaled by sqrt(2).
/// This makes their Euclidean inner product the matrix trace inner product.
#[derive(Clone, Copy, Debug, PartialEq)]
#[non_exhaustive]
pub enum Cone {
    Zero(usize),
    Nonnegative(usize),
    /// (t,z) with t >= ||z||₂.
    SecondOrder(usize),
    /// (u,v,z) with u,v >= 0 and 2uv >= ||z||₂².
    RotatedSecondOrder(usize),
    PositiveSemidefinite {
        order: usize,
    },
    /// Closure of {(x,y,z): y > 0, y exp(x/y) <= z}.
    Exponential,
    /// x,y >= 0 and x^alpha y^(1-alpha) >= |z|, 0 < alpha < 1.
    Power {
        alpha: f64,
    },
    /// Geometry and coordinate convention belong to the caller. Its identifier
    /// is preserved; independent blocks may use the same identifier.
    Opaque {
        id: usize,
        dimension: usize,
    },
}
impl Cone {
    pub fn dimension(self) -> usize {
        match self {
            Self::Zero(n)
            | Self::Nonnegative(n)
            | Self::SecondOrder(n)
            | Self::RotatedSecondOrder(n) => n,
            Self::PositiveSemidefinite { order } => order * (order + 1) / 2,
            Self::Exponential | Self::Power { .. } => 3,
            Self::Opaque { dimension, .. } => dimension,
        }
    }
}

pub(crate) enum Membership {
    Inside,
    Outside(Vec<f64>),
    Unknown,
}
impl Cone {
    pub(crate) fn classify(self, x: &[f64], tolerance: f64) -> Membership {
        use Membership::*;
        let separate = |w: Vec<f64>| {
            let dot: f64 = w.iter().zip(x).map(|(a, b)| a * b).sum();
            let scale: f64 = w.iter().zip(x).map(|(a, b)| (a * b).abs()).sum();
            if dot.is_finite() && dot < -tolerance * (1. + scale) {
                Outside(w)
            } else {
                Unknown
            }
        };
        let coordinate = |i: usize, sign: f64| {
            let mut w = vec![0.; x.len()];
            w[i] = sign;
            separate(w)
        };
        let psd_witness = |mut v: Vec<f64>| {
            let scale = v.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
            if !scale.is_finite() || scale == 0. {
                return Unknown;
            }
            for v in &mut v {
                *v /= scale;
            }
            let mut w = Vec::with_capacity(x.len());
            for j in 0..v.len() {
                for i in 0..=j {
                    w.push(v[i] * v[j] * if i == j { 1. } else { std::f64::consts::SQRT_2 });
                }
            }
            separate(w)
        };
        match self {
            Self::Opaque { .. } => Unknown,
            Self::Zero(_) => {
                if let Some(i) = x.iter().position(|&v| v != 0.) {
                    coordinate(i, -x[i].signum())
                } else {
                    Inside
                }
            }
            Self::Nonnegative(_) => {
                if let Some(i) = x.iter().position(|&v| v < 0.) {
                    coordinate(i, 1.)
                } else {
                    Inside
                }
            }
            Self::SecondOrder(_) => {
                let norm = x[1..].iter().fold(0.0_f64, |n, &v| n.hypot(v));
                if x[0] >= norm {
                    Inside
                } else if norm == 0. {
                    coordinate(0, 1.)
                } else {
                    separate(
                        std::iter::once(1.)
                            .chain(x[1..].iter().map(|v| -v / norm))
                            .collect(),
                    )
                }
            }
            Self::RotatedSecondOrder(_) => {
                let q = std::f64::consts::FRAC_1_SQRT_2;
                let mut v = x.to_vec();
                v[0] = q * x[0] + q * x[1];
                v[1] = q * x[0] - q * x[1];
                if !v[..2].iter().all(|v| v.is_finite()) {
                    return Unknown;
                }
                match Self::SecondOrder(x.len()).classify(&v, tolerance) {
                    Outside(mut w) => {
                        let a = w[0];
                        let b = w[1];
                        w[0] = q * a + q * b;
                        w[1] = q * a - q * b;
                        separate(w)
                    }
                    other => other,
                }
            }
            Self::Exponential => {
                let (a, b, c) = (x[0], x[1], x[2]);
                if b < 0. {
                    return coordinate(1, 1.);
                }
                if c < 0. {
                    return coordinate(2, 1.);
                }
                if b == 0. && a <= 0. {
                    return Inside;
                }
                if b > 0. && c > 0. && a / b <= c.ln() - b.ln() {
                    return Inside;
                }
                let r = if b > 0. {
                    (a / b).clamp(-700., 700.)
                } else {
                    (c.max(f64::MIN_POSITIVE).ln() - a.ln() + 1.).max(0.)
                };
                if r >= 0. {
                    separate(vec![-1., r - 1., (-r).exp()])
                } else {
                    let e = r.exp();
                    separate(vec![-e, e * (r - 1.), 1.])
                }
            }
            Self::Power { alpha: a } => {
                if x[0] < 0. {
                    return coordinate(0, 1.);
                }
                if x[1] < 0. {
                    return coordinate(1, 1.);
                }
                if x[2] == 0. {
                    return Inside;
                }
                let logz = x[2].abs().ln();
                let mut lx = x[0].ln();
                let mut ly = x[1].ln();
                if a * lx + (1. - a) * ly >= logz {
                    return Inside;
                }
                if !lx.is_finite() && !ly.is_finite() {
                    lx = logz - 1.;
                    ly = logz - 1.;
                } else if !lx.is_finite() {
                    lx = (logz - 1. - (1. - a) * ly) / a;
                } else if !ly.is_finite() {
                    ly = (logz - 1. - a * lx) / (1. - a);
                }
                let lg = a * lx + (1. - a) * ly;
                let w0 = a.ln() + lg - lx;
                let w1 = (1. - a).ln() + lg - ly;
                let scale = w0.max(w1).max(0.);
                separate(vec![
                    (w0 - scale).exp(),
                    (w1 - scale).exp(),
                    -x[2].signum() * (-scale).exp(),
                ])
            }
            Self::PositiveSemidefinite { order: n } => {
                let mut a = vec![0.; n * n];
                let mut k = 0;
                for j in 0..n {
                    for i in 0..=j {
                        let v = x[k]
                            * if i == j {
                                1.
                            } else {
                                std::f64::consts::FRAC_1_SQRT_2
                            };
                        a[i * n + j] = v;
                        a[j * n + i] = v;
                        k += 1;
                    }
                }
                let mut l = vec![0.; n * n];
                let mut d = vec![0.; n];
                for j in 0..n {
                    l[j * n + j] = 1.;
                    d[j] = a[j * n + j]
                        - (0..j)
                            .map(|k| l[j * n + k] * l[j * n + k] * d[k])
                            .sum::<f64>();
                    if d[j] < 0. {
                        let mut v = vec![0.; n];
                        v[j] = 1.;
                        for i in (0..j).rev() {
                            v[i] = -(i + 1..=j).map(|k| l[k * n + i] * v[k]).sum::<f64>();
                        }
                        return psd_witness(v);
                    }
                    for i in j + 1..n {
                        let residual = a[i * n + j]
                            - (0..j)
                                .map(|k| l[i * n + k] * l[j * n + k] * d[k])
                                .sum::<f64>();
                        if d[j] == 0. {
                            if residual != 0. {
                                let diagonal = a[i * n + i]
                                    - (0..j)
                                        .map(|k| l[i * n + k] * l[i * n + k] * d[k])
                                        .sum::<f64>();
                                let lambda = diagonal * 0.5 - (diagonal * 0.5).hypot(residual);
                                let mut v = vec![0.; n];
                                v[j] = residual;
                                v[i] = lambda;
                                for k in (0..j).rev() {
                                    v[k] = -(k + 1..n).map(|r| l[r * n + k] * v[r]).sum::<f64>();
                                }
                                return psd_witness(v);
                            }
                        } else {
                            l[i * n + j] = residual / d[j];
                        }
                        if !l[i * n + j].is_finite() {
                            return Unknown;
                        }
                    }
                }
                Inside
            }
        }
    }
}
