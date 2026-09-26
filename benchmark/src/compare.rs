//! Field-by-field comparison of two snapshots.
//!
//! A difference in any field outside `work` is a *change*: a different
//! status, outcome, size, or reduction count. Differences under `work` are
//! informational. Reduction counts are stored only when nonzero, so a count
//! missing on one side is zero. Any other field present on only one side is
//! schema drift between the two binaries and is reported, not compared.

use serde_json::{Map, Value};
use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write;

pub type Snapshot = BTreeMap<String, Value>;

#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
pub enum Format {
    Text,
    Markdown,
}

/// One differing leaf of one instance, as a dotted path.
#[derive(Debug, PartialEq)]
pub struct Delta {
    pub instance: String,
    pub path: String,
    pub base: Value,
    pub head: Value,
}
impl Delta {
    pub fn is_work(&self) -> bool {
        self.path.starts_with("work.")
    }
}

#[derive(Debug, Default)]
pub struct Comparison {
    pub instances: usize,
    pub only_base: Vec<String>,
    pub only_head: Vec<String>,
    pub deltas: Vec<Delta>,
    /// Field path to the number of instances where only one side has it.
    pub fields_only_base: BTreeMap<String, usize>,
    pub fields_only_head: BTreeMap<String, usize>,
    /// Instances whose record on either side reached the time limit.
    pub time_limited: Vec<String>,
    base: Totals,
    head: Totals,
}

/// Corpus sums of one side, over the instances both sides share.
#[derive(Debug, Default)]
struct Totals {
    reductions: BTreeMap<String, i128>,
    after: BTreeMap<String, i128>,
    outcomes: BTreeMap<String, i128>,
    work: BTreeMap<String, i128>,
}

impl Comparison {
    pub fn new(base: &Snapshot, head: &Snapshot) -> Self {
        let mut c = Self {
            only_base: base
                .keys()
                .filter(|k| !head.contains_key(*k))
                .cloned()
                .collect(),
            only_head: head
                .keys()
                .filter(|k| !base.contains_key(*k))
                .cloned()
                .collect(),
            ..Self::default()
        };
        let empty = Value::Object(Map::new());
        for (instance, b) in base {
            let Some(h) = head.get(instance) else {
                continue;
            };
            c.instances += 1;
            if [b, h]
                .iter()
                .any(|v| v["time_limit_reached"] == Value::Bool(true))
            {
                c.time_limited.push(instance.clone());
            }
            let sizes = b["after"].is_object() && h["after"].is_object();
            c.base.add(b, sizes);
            c.head.add(h, sizes);
            if b.get("status") != h.get("status") {
                c.deltas.push(Delta {
                    instance: instance.clone(),
                    path: "status".into(),
                    base: status_summary(b),
                    head: status_summary(h),
                });
                continue;
            }
            let mut path = Vec::new();
            c.diff(instance, &mut path, Some(b), Some(h), false, &empty);
        }
        c
    }

    fn diff<'a>(
        &mut self,
        instance: &str,
        path: &mut Vec<&'a str>,
        base: Option<&'a Value>,
        head: Option<&'a Value>,
        zero_default: bool,
        empty: &'a Value,
    ) {
        let zero = || Value::from(0);
        match (base, head) {
            (Some(Value::Object(b)), Some(Value::Object(h))) => {
                let keys: BTreeSet<&'a String> = b.keys().chain(h.keys()).collect();
                for key in keys {
                    let zero_default = zero_default || (path.is_empty() && key == "reductions");
                    path.push(key);
                    self.diff(instance, path, b.get(key), h.get(key), zero_default, empty);
                    path.pop();
                }
            }
            (None, Some(v)) | (Some(v), None) if zero_default => {
                let missing = if v.is_object() { empty } else { &Value::Null };
                let (b, h) = if base.is_none() {
                    (missing, v)
                } else {
                    (v, missing)
                };
                if v.is_object() {
                    self.diff(instance, path, Some(b), Some(h), true, empty);
                } else {
                    let fill = |x: &Value| if x.is_null() { zero() } else { x.clone() };
                    self.push(instance, path, fill(b), fill(h));
                }
            }
            (None, Some(_)) => *self.fields_only_head.entry(path.join(".")).or_default() += 1,
            (Some(_), None) => *self.fields_only_base.entry(path.join(".")).or_default() += 1,
            (Some(b), Some(h)) => {
                if b != h {
                    self.push(instance, path, b.clone(), h.clone());
                }
            }
            (None, None) => {}
        }
    }

    fn push(&mut self, instance: &str, path: &[&str], base: Value, head: Value) {
        self.deltas.push(Delta {
            instance: instance.to_string(),
            path: path.join("."),
            base,
            head,
        });
    }

    /// Instances whose status, outcome, size, or reductions differ.
    pub fn changed_instances(&self) -> BTreeSet<&str> {
        self.deltas
            .iter()
            .filter(|d| !d.is_work())
            .map(|d| d.instance.as_str())
            .collect()
    }

    /// Whether the two snapshots differ outside the work counters.
    pub fn has_changes(&self) -> bool {
        !self.only_base.is_empty()
            || !self.only_head.is_empty()
            || !self.changed_instances().is_empty()
    }

    pub fn render(&self, format: Format, max_rows: usize) -> String {
        let mut out = String::new();
        let md = format == Format::Markdown;
        let changed = self.changed_instances();
        let work_only: BTreeSet<&str> = self
            .deltas
            .iter()
            .filter(|d| d.is_work() && !changed.contains(d.instance.as_str()))
            .map(|d| d.instance.as_str())
            .collect();
        let verdict = if self.has_changes() {
            "Reductions changed"
        } else {
            "No reduction changes"
        };
        if md {
            let _ = writeln!(out, "## {verdict}\n");
        } else {
            let _ = writeln!(out, "{verdict}.");
        }
        let _ = writeln!(
            out,
            "Compared {} instances: {} with status, outcome, size, or reduction changes, {} with work-only changes.",
            self.instances,
            changed.len(),
            work_only.len()
        );
        for (label, list) in [
            ("Only in base", &self.only_base),
            ("Only in head", &self.only_head),
        ] {
            if !list.is_empty() {
                let _ = writeln!(out, "{label}: {}.", list.join(", "));
            }
        }
        if !self.time_limited.is_empty() {
            let _ = writeln!(
                out,
                "Warning: {} instances reached the time limit and are machine dependent: {}.",
                self.time_limited.len(),
                self.time_limited.join(", ")
            );
        }
        for (label, fields) in [
            ("Fields only in base", &self.fields_only_base),
            ("Fields only in head", &self.fields_only_head),
        ] {
            if !fields.is_empty() {
                let list: Vec<String> = fields.iter().map(|(f, n)| format!("{f} ({n})")).collect();
                let _ = writeln!(out, "{label}, not compared: {}.", list.join(", "));
            }
        }
        let section = |out: &mut String, title: &str| {
            if md {
                let _ = writeln!(out, "\n### {title}\n");
            } else {
                let _ = writeln!(out, "\n{title}");
            }
        };
        let totals = |base: &BTreeMap<String, i128>, head: &BTreeMap<String, i128>, all: bool| {
            let keys: BTreeSet<&String> = base.keys().chain(head.keys()).collect();
            keys.into_iter()
                .filter_map(|k| {
                    let (b, h) = (
                        base.get(k).copied().unwrap_or(0),
                        head.get(k).copied().unwrap_or(0),
                    );
                    (all || b != h).then(|| {
                        vec![
                            k.clone(),
                            b.to_string(),
                            h.to_string(),
                            signed(h - b),
                            percent(b, h),
                        ]
                    })
                })
                .collect::<Vec<_>>()
        };
        let columns = ["", "base", "head", "delta", "change"];
        let rows = totals(&self.base.outcomes, &self.head.outcomes, false);
        if !rows.is_empty() {
            section(&mut out, "Outcomes");
            out += &table(md, &with_first(&columns, "outcome"), &rows);
        }
        section(
            &mut out,
            "Size after presolve, summed over instances presolved on both sides",
        );
        out += &table(
            md,
            &with_first(&columns, "field"),
            &totals(&self.base.after, &self.head.after, true),
        );
        let rows = totals(&self.base.reductions, &self.head.reductions, false);
        if !rows.is_empty() {
            section(
                &mut out,
                "Reductions by rule and kind, summed over the corpus",
            );
            out += &table(md, &with_first(&columns, "rule.kind"), &rows);
        }
        let instance_rows = |work: bool| {
            self.deltas
                .iter()
                .filter(|d| d.is_work() == work)
                .map(|d| {
                    let delta = match (d.base.as_i64(), d.head.as_i64()) {
                        (Some(b), Some(h)) => signed(i128::from(h) - i128::from(b)),
                        _ => String::new(),
                    };
                    vec![
                        d.instance.clone(),
                        d.path.clone(),
                        compact(&d.base),
                        compact(&d.head),
                        delta,
                    ]
                })
                .collect::<Vec<_>>()
        };
        let rows = instance_rows(false);
        if !rows.is_empty() {
            section(&mut out, "Changed instances");
            out += &capped(
                md,
                &["instance", "field", "base", "head", "delta"],
                rows,
                max_rows,
            );
        }
        let rows = totals(&self.base.work, &self.head.work, false);
        if !rows.is_empty() {
            section(
                &mut out,
                "Work counters, informational, summed over the corpus",
            );
            out += &table(md, &with_first(&columns, "field"), &rows);
        }
        out
    }
}

impl Totals {
    /// Count `record`, and its size after presolve when `sizes` is set.
    fn add(&mut self, record: &Value, sizes: bool) {
        let label = match (record["status"].as_str(), record["outcome"].as_str()) {
            (Some("presolved"), Some(outcome)) => outcome.to_string(),
            (Some(status), _) => status.to_string(),
            _ => "unknown".into(),
        };
        *self.outcomes.entry(label).or_default() += 1;
        sum_leaves("", &record["reductions"], &mut self.reductions);
        if sizes {
            sum_leaves("", &record["after"], &mut self.after);
        }
        sum_leaves("", &record["work"], &mut self.work);
    }
}

/// Add every numeric leaf under `value` into `sums`, keyed by dotted path.
fn sum_leaves(prefix: &str, value: &Value, sums: &mut BTreeMap<String, i128>) {
    match value {
        Value::Object(map) => {
            for (key, v) in map {
                let path = if prefix.is_empty() {
                    key.clone()
                } else {
                    format!("{prefix}.{key}")
                };
                sum_leaves(&path, v, sums);
            }
        }
        Value::Number(n) => {
            if let Some(n) = n.as_i64() {
                *sums.entry(prefix.to_string()).or_default() += i128::from(n);
            }
        }
        _ => {}
    }
}

fn status_summary(record: &Value) -> Value {
    match (
        record["status"].as_str(),
        record["message"].as_str(),
        record["outcome"].as_str(),
    ) {
        (Some(status), Some(message), _) => Value::from(format!("{status}: {message}")),
        (Some(status), None, Some(outcome)) => Value::from(format!("{status}: {outcome}")),
        _ => record["status"].clone(),
    }
}

fn compact(value: &Value) -> String {
    match value {
        Value::String(s) => s.clone(),
        other => other.to_string(),
    }
}

fn signed(delta: i128) -> String {
    if delta > 0 {
        format!("+{delta}")
    } else {
        delta.to_string()
    }
}

fn percent(base: i128, head: i128) -> String {
    if base == head {
        String::new()
    } else if base == 0 {
        "new".into()
    } else {
        format!("{:+.2}%", (head - base) as f64 * 100.0 / base as f64)
    }
}

fn with_first<'a>(columns: &[&'a str], first: &'a str) -> Vec<&'a str> {
    let mut columns = columns.to_vec();
    columns[0] = first;
    columns
}

fn capped(md: bool, headers: &[&str], mut rows: Vec<Vec<String>>, max_rows: usize) -> String {
    let hidden = rows.len().saturating_sub(max_rows);
    rows.truncate(max_rows);
    let mut out = table(md, headers, &rows);
    if hidden > 0 {
        let _ = writeln!(out, "\n{hidden} more rows not shown.");
    }
    out
}

/// A table with the first column left aligned and the rest right aligned.
fn table(md: bool, headers: &[&str], rows: &[Vec<String>]) -> String {
    let mut out = String::new();
    if md {
        let escape = |s: &str| s.replace('|', "\\|");
        let _ = writeln!(out, "| {} |", headers.join(" | "));
        let align: Vec<&str> = (0..headers.len())
            .map(|i| if i == 0 { ":--" } else { "--:" })
            .collect();
        let _ = writeln!(out, "| {} |", align.join(" | "));
        for row in rows {
            let cells: Vec<String> = row.iter().map(|c| escape(c)).collect();
            let _ = writeln!(out, "| {} |", cells.join(" | "));
        }
        return out;
    }
    let mut widths: Vec<usize> = headers.iter().map(|h| h.chars().count()).collect();
    for row in rows {
        for (w, cell) in widths.iter_mut().zip(row) {
            *w = (*w).max(cell.chars().count());
        }
    }
    let line = |cells: Vec<&str>| {
        let padded: Vec<String> = cells
            .iter()
            .zip(&widths)
            .enumerate()
            .map(|(i, (c, &w))| {
                if i == 0 {
                    format!("{c:<w$}")
                } else {
                    format!("{c:>w$}")
                }
            })
            .collect();
        padded.join("  ").trim_end().to_string()
    };
    let _ = writeln!(out, "{}", line(headers.to_vec()));
    for row in rows {
        let _ = writeln!(out, "{}", line(row.iter().map(String::as_str).collect()));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    fn snapshot(records: &[(&str, Value)]) -> Snapshot {
        records
            .iter()
            .map(|(k, v)| (k.to_string(), v.clone()))
            .collect()
    }

    fn presolved(variables: u64, reductions: Value, phases: u64) -> Value {
        json!({
            "status": "presolved",
            "outcome": "reduced",
            "time_limit_reached": false,
            "after": {"variables": variables, "a_nonzeros": 10},
            "reductions": reductions,
            "work": {"fast_phases": phases},
        })
    }

    #[test]
    fn identical_snapshots_have_no_changes() {
        let s = snapshot(&[(
            "a/x",
            presolved(3, json!({"fixed_variables": {"fixed": 2}}), 1),
        )]);
        let c = Comparison::new(&s, &s);
        assert!(!c.has_changes());
        assert!(c.deltas.is_empty());
        assert!(
            c.render(Format::Text, 10)
                .starts_with("No reduction changes.")
        );
    }

    #[test]
    fn missing_reduction_counts_are_zero() {
        let base = snapshot(&[(
            "a/x",
            presolved(3, json!({"fixed_variables": {"fixed": 2}}), 1),
        )]);
        let head = snapshot(&[(
            "a/x",
            presolved(
                2,
                json!({"fixed_variables": {"fixed": 2, "substituted": 1}, "cones": {"cone_slack": 4}}),
                1,
            ),
        )]);
        let c = Comparison::new(&base, &head);
        let paths: Vec<(&str, &Value, &Value)> = c
            .deltas
            .iter()
            .map(|d| (d.path.as_str(), &d.base, &d.head))
            .collect();
        assert_eq!(
            paths,
            [
                ("after.variables", &json!(3), &json!(2)),
                ("reductions.cones.cone_slack", &json!(0), &json!(4)),
                (
                    "reductions.fixed_variables.substituted",
                    &json!(0),
                    &json!(1)
                ),
            ]
        );
        assert!(c.fields_only_head.is_empty());
        assert!(c.has_changes());
        // And in the other direction.
        let c = Comparison::new(&head, &base);
        assert_eq!(c.deltas[1].path, "reductions.cones.cone_slack");
        assert_eq!(
            (&c.deltas[1].base, &c.deltas[1].head),
            (&json!(4), &json!(0))
        );
    }

    #[test]
    fn work_only_differences_are_informational() {
        let base = snapshot(&[("a/x", presolved(3, json!({}), 1))]);
        let head = snapshot(&[("a/x", presolved(3, json!({}), 2))]);
        let c = Comparison::new(&base, &head);
        assert!(!c.has_changes());
        assert_eq!(c.deltas.len(), 1);
        assert!(c.deltas[0].is_work());
        let text = c.render(Format::Text, 10);
        assert!(
            text.contains(
                "0 with status, outcome, size, or reduction changes, 1 with work-only changes"
            ),
            "{text}"
        );
        assert!(
            text.contains("work.fast_phases") || text.contains("fast_phases"),
            "{text}"
        );
    }

    #[test]
    fn fields_on_one_side_are_schema_drift() {
        let mut head_record = presolved(3, json!({}), 1);
        head_record["work"]["new_counter"] = json!(7);
        head_record["extra"] = json!(true);
        let base = snapshot(&[
            ("a/x", presolved(3, json!({}), 1)),
            ("a/y", presolved(1, json!({}), 1)),
        ]);
        let head = snapshot(&[("a/x", head_record), ("b/z", presolved(1, json!({}), 1))]);
        let c = Comparison::new(&base, &head);
        assert!(c.deltas.is_empty());
        assert_eq!(
            c.fields_only_head.keys().collect::<Vec<_>>(),
            ["extra", "work.new_counter"]
        );
        assert_eq!(
            (c.only_base.as_slice(), c.only_head.as_slice()),
            (&["a/y".to_string()][..], &["b/z".to_string()][..])
        );
        // A different instance set is a change even with identical records.
        assert!(c.has_changes());
    }

    #[test]
    fn a_status_change_is_one_delta() {
        let base = snapshot(&[(
            "a/x",
            presolved(3, json!({"fixed_variables": {"fixed": 2}}), 1),
        )]);
        let head = snapshot(&[(
            "a/x",
            json!({"status": "panic", "message": "index out of bounds"}),
        )]);
        let c = Comparison::new(&base, &head);
        assert_eq!(c.deltas.len(), 1);
        assert_eq!(c.deltas[0].path, "status");
        assert_eq!(c.deltas[0].base, json!("presolved: reduced"));
        assert_eq!(c.deltas[0].head, json!("panic: index out of bounds"));
        assert!(c.fields_only_base.is_empty());
        let md = c.render(Format::Markdown, 10);
        assert!(md.starts_with("## Reductions changed\n"), "{md}");
        assert!(md.contains("| reduced | 1 | 0 | -1 | -100.00% |"), "{md}");
        assert!(md.contains("| panic | 0 | 1 | +1 | new |"), "{md}");
    }

    #[test]
    fn a_certificate_replacing_a_size_is_a_leaf_change() {
        let base = snapshot(&[("a/x", presolved(3, json!({}), 1))]);
        let mut infeasible = presolved(3, json!({}), 1);
        infeasible["outcome"] = json!("infeasible");
        infeasible["after"] = Value::Null;
        let head = snapshot(&[("a/x", infeasible)]);
        let c = Comparison::new(&base, &head);
        let paths: Vec<&str> = c.deltas.iter().map(|d| d.path.as_str()).collect();
        assert_eq!(paths, ["after", "outcome"]);
        assert!(c.fields_only_base.is_empty());
        // The instance has no size on one side, so neither side sums it.
        assert!(c.base.after.is_empty() && c.head.after.is_empty());
    }

    #[test]
    fn long_instance_tables_are_capped() {
        let base: Snapshot = (0..5)
            .map(|i| (format!("a/{i}"), presolved(3, json!({}), 1)))
            .collect();
        let head: Snapshot = (0..5)
            .map(|i| (format!("a/{i}"), presolved(4, json!({}), 1)))
            .collect();
        let text = Comparison::new(&base, &head).render(Format::Text, 2);
        assert!(text.contains("3 more rows not shown."), "{text}");
    }
}
