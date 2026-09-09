/// Native row and variable bounds, including unbounded sides.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Bounds {
    pub lower: f64,
    pub upper: f64,
}

impl Bounds {
    pub const FREE: Self = Self {
        lower: f64::NEG_INFINITY,
        upper: f64::INFINITY,
    };

    pub fn fixed(value: f64) -> Self {
        Self {
            lower: value,
            upper: value,
        }
    }

    pub fn equality(self) -> bool {
        self.lower.is_finite() && self.lower == self.upper
    }

    pub fn contains(self, x: f64) -> bool {
        self.lower <= x && x <= self.upper
    }

    pub fn recession(self) -> Self {
        Self {
            lower: if self.lower.is_finite() {
                0.0
            } else {
                f64::NEG_INFINITY
            },
            upper: if self.upper.is_finite() {
                0.0
            } else {
                f64::INFINITY
            },
        }
    }
}
