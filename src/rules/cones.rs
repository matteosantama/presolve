//! Incremental, exact structural cone simplification.
use crate::problem::{Bounds, Cone};
use crate::{
    model::tape::{Certificate, Rule, SocDirection},
    model::{Model, RowDomain},
};

impl Model {
    /// Free unit-coordinate tails with one common external linear row can
    /// be replaced by its normalized direction. Orthogonal directions have
    /// no effect outside the cone and can be zero in the recovered solution.
    fn simplify_soc_tail(&mut self, rows: &[usize]) -> Vec<usize> {
        let mut kept = Vec::with_capacity(rows.len());
        kept.push(rows[0]);
        let mut candidates = Vec::new();
        for (at, &row) in rows[1..].iter().enumerate() {
            if at != 0
                && at.is_multiple_of(256)
                && self
                    .deadline
                    .is_some_and(|d| std::time::Instant::now() >= d)
            {
                kept.extend_from_slice(&rows[at + 1..]);
                return kept;
            }
            let entries = self.a.row(row);
            if entries.is_empty() && self.cone_rhs(row) == 0.0 {
                self.remove_cone_row(row);
                continue;
            }
            kept.push(row);
            if entries.len() != 1 {
                continue;
            }
            let (column, coefficient) = entries.iter().next().unwrap();
            if coefficient.abs() != 1.0
                || self.a.column(column).len() != 2
                || self.bounds[column] != Bounds::FREE
                || self.objective.c[column] != 0.0
                || !self.objective.p.column(column).is_empty()
                || self.cone_rhs(row) != 0.0
            {
                continue;
            }
            let (external, a) = self
                .a
                .column(column)
                .iter()
                .find(|&(i, _)| i != row)
                .unwrap();
            if matches!(self.rows[external], RowDomain::Linear(_)) {
                candidates.push((external, column, row, a, -coefficient));
            }
        }
        if candidates.len() < 2 {
            return kept;
        }
        let mut aggregated = false;
        candidates.sort_unstable_by_key(|&(external, column, ..)| (external, column));
        let mut start = 0;
        while start < candidates.len() {
            if self
                .deadline
                .is_some_and(|d| std::time::Instant::now() >= d)
            {
                break;
            }
            let external = candidates[start].0;
            let mut end = start + 1;
            while end < candidates.len() && candidates[end].0 == external {
                end += 1;
            }
            let group = &candidates[start..end];
            start = end;
            if group.len() < 2 {
                continue;
            }
            let norm = group
                .iter()
                .fold(0.0_f64, |norm, &(_, _, _, a, _)| norm.hypot(a));
            // Subnormal norms may have a large relative rounding error,
            // which would no longer give a unit direction in the cone.
            if !norm.is_normal() || group.iter().any(|&(_, _, _, a, _)| a / norm == 0.0) {
                continue;
            }
            let column = group[0].1;
            let row = group[0].2;
            let members: Vec<_> = group
                .iter()
                .map(|&(_, column, row, a, sign)| SocDirection {
                    column,
                    row,
                    weight: a / norm,
                    sign,
                })
                .collect();
            // Members are ordered by column, so remove all but the retained
            // direction in one exterior-row merge, without a dense marker.
            let mut member = 0;
            let entries: Vec<_> = self
                .a
                .row(external)
                .iter()
                .filter_map(|(j, a)| {
                    if member < members.len() && members[member].column == j {
                        member += 1;
                        (j == column).then_some((j, norm))
                    } else {
                        Some((j, a))
                    }
                })
                .collect();
            if self
                .deadline
                .is_some_and(|d| std::time::Instant::now() >= d)
            {
                break;
            }
            self.replace_row(external, &entries, self.rows[external]);
            self.replace_row(row, &[(column, -1.0)], self.rows[row]);
            for member in &members[1..] {
                self.clear_row(member.row);
                self.alive[member.column] = false;
            }
            aggregated = true;
            self.postsolve.rules.push(Rule::SocAggregated {
                column,
                row,
                members,
            });
        }
        if aggregated {
            kept.retain(|&i| self.rows[i] != RowDomain::Deleted);
        }
        kept
    }

    fn cone_rhs(&self, i: usize) -> f64 {
        match self.rows[i] {
            RowDomain::Cone { rhs, .. } => rhs,
            _ => unreachable!(),
        }
    }
    fn zero_coordinate(&self, i: usize) -> bool {
        self.a.row(i).is_empty() && self.cone_rhs(i) == 0.0
    }
    fn remove_cone_row(&mut self, i: usize) {
        let rhs = self.cone_rhs(i);
        if rhs != 0. || !self.a.row(i).is_empty() {
            self.postsolve.rules.push(Rule::ConeSlack {
                row: i,
                rhs,
                entries: self.a.row(i).to_vec(),
            });
        }
        self.delete_row(i);
    }
    fn linear_cone_row(&mut self, i: usize, equality: bool) {
        let rhs = self.cone_rhs(i);
        let entries = self.a.row(i).to_vec();
        self.replace_row(
            i,
            &entries,
            RowDomain::Linear(if equality {
                Bounds::fixed(rhs)
            } else {
                Bounds {
                    lower: f64::NEG_INFINITY,
                    upper: rhs,
                }
            }),
        );
        self.postsolve.rules.push(Rule::ConeSlack {
            row: i,
            rhs,
            entries,
        });
    }
    pub fn simplify_cones(&mut self) -> Result<(), Certificate> {
        let mut visited = 0_usize;
        while let Some(block) = self.changed_cones.pop() {
            if visited.is_multiple_of(256)
                && self
                    .deadline
                    .is_some_and(|d| std::time::Instant::now() >= d)
            {
                self.changed_cones.push(block);
                break;
            }
            visited += 1;
            if block >= self.cones.len() {
                continue;
            }
            let rows: Vec<_> = self.cone_rows[block]
                .iter()
                .copied()
                .filter(|&i| matches!(self.rows[i], RowDomain::Cone { block: b, .. } if b == block))
                .collect();
            if rows.is_empty() {
                continue;
            }
            let cone = self.cones[block];
            if matches!(cone, Cone::Opaque { .. }) {
                continue;
            }
            // A separating dual vector is also an exact certificate in the
            // constant block's current coordinates, recovered by the tape.
            if rows.iter().all(|&i| self.a.row(i).is_empty()) {
                let rhs: Vec<_> = rows.iter().map(|&i| self.cone_rhs(i)).collect();
                match cone.classify(&rhs, self.settings.numerics.feasibility) {
                    crate::problem::Membership::Inside => {
                        for i in rows {
                            self.remove_cone_row(i);
                        }
                        continue;
                    }
                    crate::problem::Membership::Outside(dual) => {
                        return Err(self
                            .primal_certificate(rows.iter().zip(dual).map(|(&i, v)| (i, -v)), []));
                    }
                    crate::problem::Membership::Unknown => (),
                }
            }
            match cone {
                Cone::Zero(_) | Cone::Nonnegative(_) => {
                    for i in rows {
                        self.linear_cone_row(i, matches!(cone, Cone::Zero(_)));
                    }
                }
                Cone::SecondOrder(_) => {
                    let head = rows[0];
                    if self.a.row(head).is_empty()
                        && self.cone_rhs(head) < -self.settings.numerics.feasibility
                    {
                        return Err(self.primal_certificate([(head, -1.)], []));
                    }
                    if self.zero_coordinate(head) {
                        self.clear_row(head);
                        for &i in &rows[1..] {
                            self.linear_cone_row(i, true);
                        }
                        self.postsolve.rules.push(Rule::SocFace {
                            head,
                            tail: rows[1..].to_vec(),
                        });
                        continue;
                    }
                    let kept = self.simplify_soc_tail(&rows);
                    self.cones[block] = Cone::SecondOrder(kept.len());
                    self.cone_rows[block] = kept.clone();
                    if kept.len() == 1 {
                        self.linear_cone_row(head, false);
                    } else if kept.len() == 2 {
                        let tail = kept[1];
                        let h0 = self.cone_rhs(head);
                        let h1 = self.cone_rhs(tail);
                        let r0 = self.a.row(head).to_vec();
                        let r1 = self.a.row(tail).to_vec();
                        let combine = |sign: f64| {
                            let mut entries = std::collections::BTreeMap::new();
                            for &(j, v) in &r0 {
                                entries.insert(j, v);
                            }
                            for &(j, v) in &r1 {
                                *entries.entry(j).or_insert(0.) += sign * v;
                            }
                            entries
                                .into_iter()
                                .filter(|&(_, v)| v != 0.)
                                .collect::<Vec<_>>()
                        };
                        let plus = combine(1.);
                        let minus = combine(-1.);
                        if !plus.iter().chain(&minus).all(|&(_, v)| v.is_finite())
                            || !(h0 + h1).is_finite()
                            || !(h0 - h1).is_finite()
                        {
                            continue;
                        }
                        self.replace_row(
                            head,
                            &plus,
                            RowDomain::Linear(Bounds {
                                lower: f64::NEG_INFINITY,
                                upper: h0 + h1,
                            }),
                        );
                        self.replace_row(
                            tail,
                            &minus,
                            RowDomain::Linear(Bounds {
                                lower: f64::NEG_INFINITY,
                                upper: h0 - h1,
                            }),
                        );
                        self.postsolve.rules.push(Rule::SocToLinear { head, tail });
                        self.postsolve.rules.push(Rule::ConeSlack {
                            row: head,
                            rhs: h0,
                            entries: r0,
                        });
                        self.postsolve.rules.push(Rule::ConeSlack {
                            row: tail,
                            rhs: h1,
                            entries: r1,
                        });
                    }
                }
                Cone::RotatedSecondOrder(_) => {
                    let mut kept = rows[..2].to_vec();
                    for &i in &rows[2..] {
                        if self.zero_coordinate(i) {
                            self.remove_cone_row(i);
                        } else {
                            kept.push(i);
                        }
                    }
                    self.cones[block] = Cone::RotatedSecondOrder(kept.len());
                    self.cone_rows[block] = kept.clone();
                    if kept.len() == 2 {
                        for i in kept {
                            self.linear_cone_row(i, false);
                        }
                    }
                }
                Cone::PositiveSemidefinite { order } => {
                    // At the zero face, free equality multipliers on all off-
                    // diagonals always have a finite PSD completion on the
                    // identically-zero diagonals (Gershgorin). Partial faces
                    // need not admit a finite dual lift, so are not used here.
                    if (0..order).all(|j| self.zero_coordinate(rows[j * (j + 1) / 2 + j])) {
                        let mut k = 0;
                        for j in 0..order {
                            for i in 0..=j {
                                if i == j {
                                    self.clear_row(rows[k]);
                                } else {
                                    self.linear_cone_row(rows[k], true);
                                }
                                k += 1;
                            }
                        }
                        self.postsolve.rules.push(Rule::PsdZeroFace { rows, order });
                        continue;
                    }
                    // Connected components of structurally nonzero off-diagonal
                    // expressions identify independent principal blocks.
                    let mut parent: Vec<_> = (0..order).collect();
                    fn root(p: &[usize], mut i: usize) -> usize {
                        while p[i] != i {
                            i = p[i];
                        }
                        i
                    }
                    let mut k = 0;
                    for j in 0..order {
                        for i in 0..=j {
                            if i != j && !self.zero_coordinate(rows[k]) {
                                let a = root(&parent, i);
                                let b = root(&parent, j);
                                parent[b] = a;
                            }
                            k += 1;
                        }
                    }
                    let mut components = std::collections::BTreeMap::<usize, Vec<usize>>::new();
                    for i in 0..order {
                        components.entry(root(&parent, i)).or_default().push(i);
                    }
                    if components.len() == 1 && order != 1 {
                        continue;
                    }
                    let mut kept = vec![false; rows.len()];
                    for component in components.values() {
                        if component.len() == 1 {
                            let j = component[0];
                            let k = j * (j + 1) / 2 + j;
                            kept[k] = true;
                            self.linear_cone_row(rows[k], false);
                        } else {
                            let new_block = self.cones.len();
                            let mut new_rows = Vec::new();
                            for (at, &j) in component.iter().enumerate() {
                                for &i in &component[..=at] {
                                    let k = j * (j + 1) / 2 + i;
                                    kept[k] = true;
                                    let row = rows[k];
                                    let rhs = self.cone_rhs(row);
                                    self.rows[row] = RowDomain::Cone {
                                        rhs,
                                        block: new_block,
                                    };
                                    new_rows.push(row);
                                }
                            }
                            self.cones.push(Cone::PositiveSemidefinite {
                                order: component.len(),
                            });
                            self.cone_rows.push(new_rows);
                            self.changed_cones.push(new_block);
                        }
                    }
                    for (k, &i) in rows.iter().enumerate() {
                        if !kept[k] {
                            self.remove_cone_row(i);
                        }
                    }
                    self.cone_rows[block].clear();
                    self.revision += 1;
                }
                _ => (),
            }
        }
        Ok(())
    }
}
