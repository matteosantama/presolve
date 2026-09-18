//! Batch stable LP equality trees before falling back to individual pivots.
//! Eliminated variables must be free; one bounded representative is allowed.
//! Every incident row is rebuilt once, and the tape retains the original
//! equations and incidences for primal, dual, ray, and warm-start recovery.

use crate::{
    matrix::sparse::Entries,
    model::{
        Model, RowDomain, shifted,
        tape::{DoubletonStep, Rule},
    },
    problem::Bounds,
};
use std::time::Instant;

const MIN_EDGES: usize = 8;
const NONE: usize = usize::MAX;

#[derive(Clone, Copy)]
struct Edge {
    row: usize,
    columns: [usize; 2],
    values: [f64; 2],
    rhs: f64,
}

impl Model {
    pub(super) fn batch_doubleton_chains(&mut self, round: &[usize]) {
        // Avoid graph setup on QPs, conic models, and sparse late frontiers.
        // These cases retain the ordinary transactional doubleton path.
        if round.len() < MIN_EDGES
            || round.len().saturating_mul(16) < self.alive.len()
            || self.objective.p.nnz() != 0
            || !self.cones.is_empty()
            // The normal pipeline drains cheap bound cleanup first. With
            // singleton rows disabled, retain scalar substitutions instead of
            // building a graph for variables whose bounds remain in rows.
            || !self.queues.singleton_rows.is_empty()
        {
            return;
        }
        let mut edges = Vec::new();
        for (position, &row) in round.iter().enumerate() {
            if position.is_multiple_of(256) && self.deadline.is_some_and(|d| Instant::now() >= d) {
                return;
            }
            let RowDomain::Linear(bounds) = self.rows[row] else {
                continue;
            };
            if !bounds.equality() || self.a.row(row).len() != 2 {
                continue;
            }
            let mut entries = self.a.row(row).iter();
            let (j, a) = entries.next().unwrap();
            let (k, b) = entries.next().unwrap();
            // Unit-magnitude ratios neither amplify nor underflow along a
            // long chain. General ratios stay on the existing pivot path.
            if a.abs() != b.abs()
                || !a.is_finite()
                || !bounds.lower.is_finite()
                || (self.bounds[j] != Bounds::FREE && self.bounds[k] != Bounds::FREE)
            {
                continue;
            }
            edges.push(Edge {
                row,
                columns: [j, k],
                values: [a, b],
                rhs: bounds.lower,
            });
        }
        if edges.len() < MIN_EDGES {
            return;
        }
        let mut head = vec![NONE; self.alive.len()];
        let mut next = vec![NONE; edges.len() * 2];
        for (e, edge) in edges.iter().enumerate() {
            for side in 0..2 {
                let j = edge.columns[side];
                next[2 * e + side] = head[j];
                head[j] = 2 * e + side;
            }
        }
        let mut seen = vec![false; self.alive.len()];
        let mut mapping = Vec::new();
        let mut component = Vec::new();
        let mut rows = Vec::new();
        for edge in &edges {
            let start = edge.columns[0];
            if seen[start] {
                continue;
            }
            if self.deadline.is_some_and(|d| Instant::now() >= d) {
                break;
            }
            component.clear();
            rows.clear();
            component.push(start);
            seen[start] = true;
            let mut cursor = 0;
            let mut incidences = 0;
            while cursor < component.len() {
                if cursor.is_multiple_of(1024) && self.deadline.is_some_and(|d| Instant::now() >= d)
                {
                    return;
                }
                let j = component[cursor];
                let mut at = head[j];
                while at != NONE {
                    incidences += 1;
                    let edge = &edges[at / 2];
                    let other = edge.columns[1 - at % 2];
                    if !seen[other] {
                        seen[other] = true;
                        component.push(other);
                        rows.push(edge.row);
                    }
                    at = next[at];
                }
                cursor += 1;
            }
            // Cycles require consistency handling; ordinary substitutions
            // expose them as singleton/empty rows and retain their proofs.
            if rows.len() < MIN_EDGES || incidences / 2 != rows.len() {
                continue;
            }
            let mut bounded = component
                .iter()
                .copied()
                .filter(|&j| self.bounds[j] != Bounds::FREE);
            let first = bounded.next();
            if bounded.next().is_some() {
                continue;
            }
            let root = first.unwrap_or_else(|| {
                *component
                    .iter()
                    .max_by_key(|&&j| (self.a.column(j).len(), j))
                    .unwrap()
            });
            rows.sort_unstable();
            if mapping.is_empty() {
                mapping.resize(self.alive.len(), None);
            }
            // Parent-before-child traversal establishes the affine map to
            // the representative, but stores original equations on the tape.
            let mut order = vec![(root, NONE)];
            let mut steps = Vec::with_capacity(rows.len());
            mapping[root] = Some((0.0, 1.0));
            let mut valid = true;
            let mut cursor = 0;
            while cursor < order.len() {
                if cursor.is_multiple_of(1024) && self.deadline.is_some_and(|d| Instant::now() >= d)
                {
                    return;
                }
                let (parent, parent_edge) = order[cursor];
                let (offset, slope) = mapping[parent].unwrap();
                let mut at = head[parent];
                while at != NONE {
                    let e = at / 2;
                    if e != parent_edge {
                        let edge = edges[e];
                        let side = 1 - at % 2;
                        let column = edge.columns[side];
                        let pivot = edge.values[side];
                        let other = edge.values[1 - side];
                        let ratio = -other / pivot;
                        let child_offset = edge.rhs / pivot + ratio * offset;
                        let child_slope = ratio * slope;
                        if !child_offset.is_finite() || !child_slope.is_finite() {
                            valid = false;
                            break;
                        }
                        mapping[column] = Some((child_offset, child_slope));
                        order.push((column, e));
                        steps.push(DoubletonStep {
                            column,
                            parent,
                            row: edge.row,
                            pivot,
                            other,
                            rhs: edge.rhs,
                            objective: self.objective.c[column],
                            entries: self
                                .a
                                .column(column)
                                .iter()
                                .filter(|&(i, _)| i != edge.row)
                                .collect(),
                        });
                    }
                    at = next[at];
                }
                if !valid {
                    break;
                }
                cursor += 1;
            }
            // The representative remains an ordinary variable during row
            // merging; its existing coefficient is the initial accumulator.
            mapping[root] = None;
            if valid {
                self.commit_doubleton_chain(root, &rows, steps, &mapping);
            }
            for &j in &component {
                mapping[j] = None;
            }
        }
    }

    fn commit_doubleton_chain(
        &mut self,
        root: usize,
        chain_rows: &[usize],
        steps: Vec<DoubletonStep>,
        mapping: &[Option<(f64, f64)>],
    ) -> bool {
        let mut affected: Vec<_> = steps
            .iter()
            .flat_map(|s| s.entries.iter().map(|&(i, _)| i))
            .collect();
        affected.sort_unstable();
        affected.dedup();
        let mut updates = Vec::with_capacity(affected.len() + chain_rows.len());
        let mut fill = 0;
        for i in affected {
            if chain_rows.binary_search(&i).is_ok() {
                continue;
            }
            if self.deadline.is_some_and(|d| Instant::now() >= d) {
                return false;
            }
            let old_root = self.a.get(i, root);
            let mut root_value = old_root;
            let mut shift = 0.0;
            let mut entries: Entries = Vec::with_capacity(self.a.row(i).len() + 1);
            for (j, a) in self.a.row(i) {
                if let Some((offset, slope)) = mapping[j] {
                    let change = a * slope;
                    let value = root_value + change;
                    if !value.is_finite()
                        || (value != 0.0
                            && value.abs() < 1e-10 * root_value.abs().max(change.abs()))
                    {
                        return false;
                    }
                    root_value = value;
                    shift += a * offset;
                    if !shift.is_finite() {
                        return false;
                    }
                } else if j != root {
                    entries.push((j, a));
                }
            }
            if root_value != 0.0 {
                fill += usize::from(old_root == 0.0);
                // A conservative aggregate cap never grants more fill than a
                // single allowed pivot. Rejected batches fall back unchanged.
                if fill > self.settings.substitution_fill {
                    return false;
                }
                let at = entries.partition_point(|&(j, _)| j < root);
                entries.insert(at, (root, root_value));
            }
            let Some(domain) = shifted(self.rows[i], shift) else {
                return false;
            };
            updates.push((i, entries, domain));
        }
        let mut cost = self.objective.c[root];
        let mut constant = self.objective.constant;
        for step in &steps {
            let (offset, slope) = mapping[step.column].unwrap();
            cost += step.objective * slope;
            constant += step.objective * offset;
            if !cost.is_finite() || !constant.is_finite() {
                return false;
            }
        }
        if self.deadline.is_some_and(|d| Instant::now() >= d) {
            return false;
        }
        for &i in chain_rows {
            updates.push((i, Vec::new(), RowDomain::Deleted));
        }
        for step in &steps {
            self.alive[step.column] = false;
            self.objective.c[step.column] = 0.0;
        }
        self.objective.c[root] = cost;
        self.objective.constant = constant;
        self.replace_rows_batch(updates);
        self.postsolve.rules.push(Rule::DoubletonChain { steps });
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        matrix::{linked::LinkedMatrix, sparse::SymmetricMatrix},
        model::{
            objective::Objective,
            tape::{Point, Recovery},
        },
    };

    fn fixture(sign: f64) -> (Model, Point) {
        fixture_tree(sign, false)
    }

    fn fixture_tree(sign: f64, branched: bool) -> (Model, Point) {
        let n = 18;
        let m = 18;
        let x: Vec<_> = (0..n).map(|j| j as f64).collect();
        let mut coefficients = vec![Vec::new(); n];
        let mut rows = Vec::new();
        for j in 1..17 {
            let parent = if branched { (j - 1) / 2 } else { j - 1 };
            coefficients[parent].push((j - 1, -sign));
            coefficients[j].push((j - 1, 1.0));
            rows.push(RowDomain::Linear(Bounds::fixed(x[j] - sign * x[parent])));
        }
        for (j, column) in coefficients.iter_mut().enumerate() {
            column.push((16, if j == 17 { 3.0 } else { 1.0 }));
            column.push((17, if j % 2 == 0 { 2.0 } else { 1.0 }));
        }
        let activity = |i| {
            coefficients
                .iter()
                .zip(&x)
                .map(|(a, x)| {
                    a.iter()
                        .find(|&&(r, _)| r == i)
                        .map_or(0.0, |&(_, a)| a * x)
                })
                .sum()
        };
        rows.push(RowDomain::Linear(Bounds {
            lower: activity(16),
            upper: f64::INFINITY,
        }));
        rows.push(RowDomain::Linear(Bounds {
            lower: f64::NEG_INFINITY,
            upper: activity(17),
        }));
        let y: Vec<_> = (0..m)
            .map(|i| if i == 17 { -1.0 } else { (i % 3 + 1) as f64 })
            .collect();
        let mut z = vec![0.0; n];
        z[0] = 2.0;
        let c = coefficients
            .iter()
            .enumerate()
            .map(|(j, a)| a.iter().map(|&(i, a)| a * y[i]).sum::<f64>() + z[j])
            .collect();
        let a = LinkedMatrix::from_columns(m, n, |j| coefficients[j].iter().copied());
        let mut bounds = vec![Bounds::FREE; n];
        bounds[0] = Bounds {
            lower: 0.0,
            upper: 100.0,
        };
        let model = Model::from_parts(
            a,
            Objective {
                p: SymmetricMatrix::from_upper_columns(n, |_| std::iter::empty()),
                c,
                constant: 3.0,
                scratch: Default::default(),
            },
            rows,
            bounds,
        );
        (model, Point { x, y, z })
    }

    fn stationary(model: &Model, point: &Point, objective: bool) {
        for j in 0..model.alive.len() {
            if !model.alive[j] {
                continue;
            }
            let residual = if objective { model.objective.c[j] } else { 0.0 }
                - model
                    .a
                    .column(j)
                    .iter()
                    .map(|(i, a)| a * point.y[i])
                    .sum::<f64>()
                - point.z[j];
            assert!(residual.abs() < 1e-9, "column {j}: {residual}");
        }
    }

    #[test]
    fn batch_rewrites_rows_once_and_preserves_objective_primal_dual_and_warm_start() {
        for (sign, branched) in [(-1.0, false), (1.0, false), (-1.0, true), (1.0, true)] {
            let (mut model, original) = fixture_tree(sign, branched);
            let (source, _) = fixture_tree(sign, branched);
            let value = 3.0
                + model
                    .objective
                    .c
                    .iter()
                    .zip(&original.x)
                    .map(|(c, x)| c * x)
                    .sum::<f64>();
            model.batch_doubleton_chains(&(0..16).collect::<Vec<_>>());
            assert!(
                matches!(model.postsolve.rules.as_slice(),[Rule::DoubletonChain{steps}] if steps.len()==16)
            );
            assert_eq!(model.alive.iter().filter(|&&a| a).count(), 2);
            assert_eq!(model.a.nnz(), 4);
            let mut warm = original.clone();
            model.postsolve.reduce_point(&mut warm);
            stationary(&model, &warm, true);
            assert_eq!(
                model.objective.constant
                    + model
                        .objective
                        .c
                        .iter()
                        .zip(&warm.x)
                        .map(|(c, x)| c * x)
                        .sum::<f64>(),
                value
            );
            let mut recovered = warm;
            model.postsolve.recover(&mut recovered, Recovery::Solution);
            assert_eq!(recovered.x, original.x);
            assert_eq!(recovered.y, original.y);
            assert_eq!(recovered.z, original.z);
            stationary(&source, &recovered, true);
            // Rays omit all offsets in the chain, including alternating signs.
            let mut ray = Point::zeros(18, 18);
            ray.x[0] = 1.0;
            ray.x[17] = 2.0;
            model
                .postsolve
                .recover(&mut ray, Recovery::DualInfeasibility);
            for i in 0..16 {
                assert_eq!(
                    source
                        .a
                        .row(i)
                        .iter()
                        .map(|(j, a)| a * ray.x[j])
                        .sum::<f64>(),
                    0.0
                );
            }
            // Certificate recovery must omit the original objective costs.
            let mut cert = Point::zeros(18, 18);
            cert.y[16] = 1.0;
            cert.y[17] = -2.0;
            for &j in &[0, 17] {
                cert.z[j] = -model
                    .a
                    .column(j)
                    .iter()
                    .map(|(i, a)| a * cert.y[i])
                    .sum::<f64>();
            }
            model
                .postsolve
                .recover(&mut cert, Recovery::PrimalInfeasibility);
            stationary(&source, &cert, false);
            let rhs = |m: &Model, p: &Point| {
                m.rows
                    .iter()
                    .enumerate()
                    .map(|(i, r)| match r {
                        RowDomain::Linear(b) => {
                            if p.y[i] >= 0.0 {
                                b.lower * p.y[i]
                            } else {
                                b.upper * p.y[i]
                            }
                        }
                        _ => 0.0,
                    })
                    .filter(|x| x.is_finite())
                    .sum::<f64>()
            };
            assert!((rhs(&source, &cert) - rhs(&model, &cert)).abs() < 1e-9);
        }
    }

    #[test]
    fn unsupported_bounds_curvature_cycles_fill_and_overflow_leave_batch_transactional() {
        for variant in 0..7 {
            let (mut model, _) = fixture(1.0);
            match variant {
                0 => {
                    model.bounds[8] = Bounds {
                        lower: 0.0,
                        upper: 100.0,
                    }
                }
                1 => model.settings.substitution_fill = 0,
                2 => model.rows[3] = RowDomain::Linear(Bounds::fixed(f64::MAX)),
                3 => model.deadline = Some(Instant::now()),
                4 => {
                    model.objective.p =
                        SymmetricMatrix::from_upper_columns(18, |j| std::iter::once((j, 1.0)))
                }
                5 | 6 => {}
                _ => unreachable!(),
            }
            // Overflow needs two positive translations along the chain.
            if variant == 2 {
                model.rows[4] = RowDomain::Linear(Bounds::fixed(f64::MAX));
            }
            if variant == 1 {
                // Force a new representative entry in an external row.
                let entries: Entries = model.a.row(16).iter().filter(|&(j, _)| j != 0).collect();
                model.replace_rows_batch(vec![(16, entries, model.rows[16])]);
            }
            if variant == 5 {
                model.replace_rows_batch(vec![(
                    16,
                    vec![(0, 1.0), (16, -1.0)],
                    RowDomain::Linear(Bounds::fixed(-16.0)),
                )]);
            }
            if variant == 6 {
                model.replace_rows_batch(vec![(
                    16,
                    vec![(0, 1.0), (1, -1.0 + 1e-12)],
                    model.rows[16],
                )]);
            }
            let revision = model.revision;
            let c = model.objective.c.clone();
            model.batch_doubleton_chains(
                &(0..if variant == 5 { 17 } else { 16 }).collect::<Vec<_>>(),
            );
            assert_eq!(model.revision, revision, "variant {variant}");
            assert_eq!(model.objective.c, c);
            assert!(model.postsolve.rules.is_empty());
            assert!(model.alive.iter().all(|&x| x));
        }
    }
    #[test]
    fn infeasibility_after_a_batch_lifts_a_valid_bound_certificate() {
        let n = 17;
        let mut columns = vec![Vec::new(); n];
        for j in 1..n {
            columns[j - 1].push((j - 1, -1.0));
            columns[j].push((j - 1, 1.0));
        }
        for column in &mut columns {
            column.push((n - 1, 1.0));
        }
        let mut rows = vec![RowDomain::Linear(Bounds::fixed(0.0)); n];
        rows[n - 1] = RowDomain::Linear(Bounds {
            lower: f64::NEG_INFINITY,
            upper: -1.0,
        });
        let mut bounds = vec![Bounds::FREE; n];
        bounds[0] = Bounds {
            lower: 0.0,
            upper: f64::INFINITY,
        };
        let mut model = Model::from_parts(
            LinkedMatrix::from_columns(n, n, |j| columns[j].iter().copied()),
            Objective {
                p: SymmetricMatrix::from_upper_columns(n, |_| std::iter::empty()),
                c: vec![0.0; n],
                constant: 0.0,
                scratch: Default::default(),
            },
            rows,
            bounds,
        );
        model.batch_doubleton_chains(&(0..n - 1).collect::<Vec<_>>());
        assert!(matches!(
            model.postsolve.rules.first(),
            Some(Rule::DoubletonChain { .. })
        ));
        let mut certificate = model
            .singleton_rows()
            .expect_err("chain implies nonnegative sum");
        model
            .postsolve
            .recover(&mut certificate.point, certificate.mode);
        let point = certificate.point;
        assert!(point.y[n - 1] < 0.0);
        assert!(point.z[0] > 0.0);
        for (j, column) in columns.iter().enumerate() {
            let residual = column.iter().map(|&(i, a)| a * point.y[i]).sum::<f64>() + point.z[j];
            assert!(residual.abs() < 1e-12);
            if j != 0 {
                assert_eq!(point.z[j], 0.0);
            }
        }
        assert!(-point.y[n - 1] > 0.0);
    }
    #[test]
    fn pending_bound_cleanup_skips_graph_discovery_and_keeps_scalar_fallback() {
        let n = 17;
        let mut columns = vec![Vec::new(); n];
        for j in 1..n {
            columns[j - 1].push((j - 1, -1.0));
            columns[j].push((j - 1, 1.0));
        }
        for (j, column) in columns.iter_mut().enumerate() {
            column.push((n - 1 + j, 1.0));
        }
        let mut rows = vec![RowDomain::Linear(Bounds::fixed(0.0)); n - 1];
        rows.extend(vec![
            RowDomain::Linear(Bounds {
                lower: f64::NEG_INFINITY,
                upper: 100.0
            });
            n
        ]);
        for limit in [n - 2, n - 1] {
            let mut bounds = vec![Bounds::FREE; n];
            bounds[0] = Bounds {
                lower: 0.0,
                upper: f64::INFINITY,
            };
            let mut model = Model::from_parts(
                LinkedMatrix::from_columns(2 * n - 1, n, |j| columns[j].iter().copied()),
                Objective {
                    p: SymmetricMatrix::from_upper_columns(n, |_| std::iter::empty()),
                    c: vec![0.0; n],
                    constant: 0.0,
                    scratch: Default::default(),
                },
                rows.clone(),
                bounds,
            );
            model.settings.substitution_fill = limit;
            model.batch_doubleton_chains(&(0..n - 1).collect::<Vec<_>>());
            assert_eq!(model.revision, 0);
            assert!(model.postsolve.rules.is_empty());
            assert!(model.alive.iter().all(|&alive| alive));
            model.doubleton_equalities();
            assert_eq!(model.alive.iter().filter(|&&alive| alive).count(), 1);
            assert!(
                model
                    .postsolve
                    .rules
                    .iter()
                    .all(|r| !matches!(r, Rule::DoubletonChain { .. }))
            );
        }
    }
}
