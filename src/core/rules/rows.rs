// SPDX-License-Identifier: Apache-2.0
// Modified for this library; copyright and attribution notices are in NOTICE.
//! Eliminate empty rows and turn singleton rows into variable bounds.

use crate::{
    core::model::{Model, RowDomain},
    postsolve::tape::{Certificate, Point, Recovery, Side},
};

impl Model {
    pub fn empty_rows(&mut self) -> Result<(), Certificate> {
        while let Some(i) = self.queues.empty_rows.pop() {
            let RowDomain::Linear(b) = self.rows[i] else {
                continue;
            };
            if !self.a.row(i).is_empty() {
                continue;
            }
            if self.separated(b.lower, 0.0) {
                return Err(self.row_certificate(i, 1.0));
            }
            if self.separated(0.0, b.upper) {
                return Err(self.row_certificate(i, -1.0));
            }
            // Apply the feasibility tolerance to empty rows. Retaining roundoff such as
            // 0 <= -7e-18 lets the embedding magnify it into an infeasibility
            // certificate. RecoveryTape still exposes that tiny original residual.
            self.delete_row(i);
        }
        Ok(())
    }

    pub fn singleton_rows(&mut self) -> Result<(), Certificate> {
        while let Some(i) = self.queues.singleton_rows.pop() {
            let Some(equation) = self.equation(i) else {
                continue;
            };
            let &[(j, a)] = equation.entries.as_slice() else {
                continue;
            };
            let b = equation.bounds;
            let (lower, upper) = if a > 0.0 {
                (b.lower / a, b.upper / a)
            } else {
                (b.upper / a, b.lower / a)
            };
            self.implied_bound(j, Side::Lower, lower, |_| equation.clone(), false)?;
            self.implied_bound(j, Side::Upper, upper, |_| equation, false)?;
            // An arithmetic overflow or a small inconsistent interval must not
            // make a row disappear without its bounds reaching the model.
            if self.bounds[j].lower >= lower && self.bounds[j].upper <= upper {
                self.delete_row(i);
            }
        }
        Ok(())
    }

    pub(super) fn row_certificate(&self, row: usize, multiplier: f64) -> Certificate {
        let mut point = Point::zeros(self.bounds.len(), self.rows.len());
        point.y[row] = multiplier;
        for (j, a) in self.a.row(row) {
            point.z[j] = -a * multiplier;
        }
        Certificate {
            mode: Recovery::PrimalInfeasibility,
            point,
        }
    }
}
