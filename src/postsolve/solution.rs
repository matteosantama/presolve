//! Owned and borrowed solution data and infeasibility multipliers.
#[derive(Clone, Debug)]
pub struct Solution {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub z: Vec<f64>,
    pub conic_dual: Vec<f64>,
    pub conic_slack: Vec<f64>,
}
#[derive(Clone, Copy, Debug)]
pub struct SolutionRef<'a> {
    pub x: &'a [f64],
    pub y: &'a [f64],
    pub z: &'a [f64],
    pub conic_dual: &'a [f64],
    pub conic_slack: &'a [f64],
}
/// Caller-owned destination buffers. Sizes must match the postsolve dimensions.
pub struct SolutionMut<'a> {
    pub x: &'a mut [f64],
    pub y: &'a mut [f64],
    pub z: &'a mut [f64],
    pub conic_dual: &'a mut [f64],
    pub conic_slack: &'a mut [f64],
}
impl Solution {
    pub fn as_mut(&mut self) -> SolutionMut<'_> {
        SolutionMut {
            x: &mut self.x,
            y: &mut self.y,
            z: &mut self.z,
            conic_dual: &mut self.conic_dual,
            conic_slack: &mut self.conic_slack,
        }
    }
    pub fn as_ref(&self) -> SolutionRef<'_> {
        SolutionRef {
            x: &self.x,
            y: &self.y,
            z: &self.z,
            conic_dual: &self.conic_dual,
            conic_slack: &self.conic_slack,
        }
    }
}
/// Original-coordinate Farkas multipliers. Stationarity is
/// `Aᵀ y + z - Gᵀ conic_dual = 0`. Positive linear/bound multipliers
/// select lower sides; negative multipliers select upper sides. The selected
/// bound-weighted sum minus `hᵀ conic_dual` is positive for a contradiction.
#[derive(Clone, Debug)]
pub struct PrimalCertificate {
    pub y: Vec<f64>,
    pub z: Vec<f64>,
    pub conic_dual: Vec<f64>,
}
#[derive(Clone, Copy, Debug)]
pub struct CertificateRef<'a> {
    pub y: &'a [f64],
    pub z: &'a [f64],
    pub conic_dual: &'a [f64],
}
impl PrimalCertificate {
    pub fn as_ref(&self) -> CertificateRef<'_> {
        CertificateRef {
            y: &self.y,
            z: &self.z,
            conic_dual: &self.conic_dual,
        }
    }
}
