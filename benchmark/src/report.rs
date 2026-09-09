use crate::Kind;
use crate::results::{Case, Measurement, Run};
use std::collections::BTreeSet;
use std::fmt::Write;
use std::io::IsTerminal;
use tabled::{builder::Builder, settings::Style};

fn paint(text: impl std::fmt::Display, code: u8, color: bool) -> String {
    if color {
        format!("\x1b[{code}m{text}\x1b[0m")
    } else {
        text.to_string()
    }
}

fn table(headers: &[&str], rows: Vec<Vec<String>>) -> String {
    let mut builder = Builder::default();
    builder.push_record(headers.iter().copied());
    for row in rows {
        builder.push_record(row);
    }
    builder.build().with(Style::rounded()).to_string()
}

fn checked(case: &Case, trials: usize, kind: Kind) -> std::result::Result<&Measurement, String> {
    if !case.errors.is_empty() {
        return Err(case.errors.join("; "));
    }
    if case.measurements.len() != trials || trials == 0 {
        return Err(format!(
            "expected {trials} trials, found {}",
            case.measurements.len()
        ));
    }
    let first = &case.measurements[0];
    for m in &case.measurements {
        if m.time_limit_reached {
            return Err("presolve time limit reached".into());
        }
        if !m.same_result(first) {
            return Err("outcome or model size differs between trials".into());
        }
        if m.elapsed_ns.is_some() != (kind == Kind::Time) {
            return Err("measurement does not match run type".into());
        }
    }
    Ok(first)
}

fn percentile(sorted: &[f64], fraction: f64) -> f64 {
    let position = fraction * (sorted.len() - 1) as f64;
    let lower = position.floor() as usize;
    let upper = position.ceil() as usize;
    sorted[lower] + (sorted[upper] - sorted[lower]) * position.fract()
}

fn duration(ns: f64) -> String {
    if ns >= 1e9 {
        format!("{:.3}s", ns / 1e9)
    } else if ns >= 1e6 {
        format!("{:.3}ms", ns / 1e6)
    } else if ns >= 1e3 {
        format!("{:.3}µs", ns / 1e3)
    } else {
        format!("{ns:.1}ns")
    }
}

fn timing(case: &Case) -> (f64, String) {
    let mut values: Vec<_> = case
        .measurements
        .iter()
        .map(|m| m.elapsed_ns.unwrap() as f64)
        .collect();
    values.sort_by(f64::total_cmp);
    let median = percentile(&values, 0.5);
    (
        median,
        format!(
            "{} [{}–{}] n={}",
            duration(median),
            duration(percentile(&values, 0.1)),
            duration(percentile(&values, 0.9)),
            values.len()
        ),
    )
}

fn percent(before: f64, after: f64) -> String {
    if before == 0. {
        "n/a (zero baseline)".into()
    } else {
        format!("{:+.2}%", 100. * (after / before - 1.))
    }
}

fn render(old: &Run, new: &Run, color: bool) -> (String, bool) {
    let mut out = format!(
        "{}: {} → {}\n",
        old.metadata.kind.as_str(),
        old.metadata.name,
        new.metadata.name
    );
    let mut issues = Vec::new();
    if !old.complete {
        issues.push(format!("{} is incomplete", old.metadata.name));
    }
    if !new.complete {
        issues.push(format!("{} is incomplete", new.metadata.name));
    }
    if old.metadata.settings != new.metadata.settings {
        issues.push(format!(
            "settings differ:\n  {}: {}\n  {}: {}",
            old.metadata.name, old.metadata.settings, new.metadata.name, new.metadata.settings
        ));
    }
    if old.metadata.kind == Kind::Time && old.metadata.machine != new.metadata.machine {
        issues.push(format!(
            "machines differ: {} versus {}",
            old.metadata.machine, new.metadata.machine
        ));
    }
    if old.metadata.kind != new.metadata.kind {
        issues.push("run types differ".into());
    }
    // Global incompatibilities make numeric comparisons misleading.
    if !issues.is_empty() {
        for issue in issues {
            writeln!(out, "{}", paint(issue, 33, color)).unwrap();
        }
        return (out, true);
    }
    let ids: BTreeSet<_> = old.cases.keys().chain(new.cases.keys()).collect();
    if ids.is_empty() {
        return (format!("{out}No cases to compare.\n"), true);
    }
    let mut rows = Vec::new();
    let mut matched = 0;
    let mut changed = 0;
    let mut time_totals = std::collections::BTreeMap::<&str, (usize, f64, f64)>::new();
    for id in ids {
        let (Some(a), Some(b)) = (old.cases.get(id), new.cases.get(id)) else {
            let missing = if old.cases.contains_key(id) {
                &new.metadata.name
            } else {
                &old.metadata.name
            };
            issues.push(format!("{id}: missing from {missing}"));
            continue;
        };
        if a.input_hash != b.input_hash {
            issues.push(format!("{id}: input file changed"));
            continue;
        }
        let (ma, mb) = match (
            checked(a, old.metadata.trials, old.metadata.kind),
            checked(b, new.metadata.trials, new.metadata.kind),
        ) {
            (Ok(ma), Ok(mb)) => (ma, mb),
            (ra, rb) => {
                if let Err(e) = ra {
                    issues.push(format!("{id} ({}): {e}", old.metadata.name));
                }
                if let Err(e) = rb {
                    issues.push(format!("{id} ({}): {e}", new.metadata.name));
                }
                continue;
            }
        };
        if ma.before != mb.before {
            issues.push(format!("{id}: input dimensions/nonzeros differ"));
            continue;
        }
        matched += 1;
        let outcome_changed = ma.outcome != mb.outcome;
        if old.metadata.kind == Kind::Time {
            let (ta, sa) = timing(a);
            let (tb, sb) = timing(b);
            let note = if outcome_changed {
                format!("{} → {}", ma.outcome, mb.outcome)
            } else if ma.after != mb.after {
                "output size changed".into()
            } else {
                String::new()
            };
            rows.push(vec![
                id.clone(),
                sa,
                sb,
                paint(
                    percent(ta, tb),
                    if tb < ta {
                        32
                    } else if tb > ta {
                        31
                    } else {
                        90
                    },
                    color,
                ),
                paint(note, 33, color),
            ]);
            // A different outcome should not look like a like-for-like speedup.
            if !outcome_changed {
                let rule = id.rsplit('/').next().unwrap();
                let total = time_totals.entry(rule).or_default();
                total.0 += 1;
                total.1 += ta;
                total.2 += tb;
            }
        } else {
            if ma.same_result(mb) {
                continue;
            }
            changed += 1;
            if outcome_changed {
                rows.push(vec![
                    id.clone(),
                    "outcome".into(),
                    ma.outcome.clone(),
                    mb.outcome.clone(),
                    paint("review", 33, color),
                ]);
            }
            match (&ma.after, &mb.after) {
                (Some(sa), Some(sb)) => {
                    for ((metric, a), b) in [
                        "variables",
                        "linear rows",
                        "conic rows",
                        "A nnz",
                        "G nnz",
                        "P nnz",
                    ]
                    .into_iter()
                    .zip(sa.values())
                    .zip(sb.values())
                    {
                        if a != b {
                            let delta = b as i128 - a as i128;
                            rows.push(vec![
                                id.clone(),
                                metric.into(),
                                a.to_string(),
                                b.to_string(),
                                paint(format!("{delta:+}"), if delta < 0 { 32 } else { 31 }, color),
                            ]);
                        }
                    }
                }
                (None, Some(_)) | (Some(_), None) => {
                    rows.push(vec![
                        id.clone(),
                        "final sizes".into(),
                        if ma.after.is_some() {
                            "available"
                        } else {
                            "unavailable"
                        }
                        .into(),
                        if mb.after.is_some() {
                            "available"
                        } else {
                            "unavailable"
                        }
                        .into(),
                        paint("review", 33, color),
                    ]);
                }
                (None, None) => (),
            }
        }
    }
    if !rows.is_empty() {
        let headers = if old.metadata.kind == Kind::Time {
            [
                "Case",
                "First: median [p10–p90]",
                "Second: median [p10–p90]",
                "Change",
                "Note",
            ]
        } else {
            [
                "Case",
                "Metric",
                "First: remaining",
                "Second: remaining",
                "Difference",
            ]
        };
        writeln!(out, "{}", table(&headers, rows)).unwrap();
    }
    if old.metadata.kind == Kind::Time {
        let totals = time_totals
            .into_iter()
            .map(|(rule, (n, a, b))| {
                vec![
                    rule.into(),
                    n.to_string(),
                    duration(a),
                    duration(b),
                    paint(
                        percent(a, b),
                        if b < a {
                            32
                        } else if b > a {
                            31
                        } else {
                            90
                        },
                        color,
                    ),
                ]
            })
            .collect();
        writeln!(
            out,
            "Sum of per-problem medians, grouped by rule (matching outcomes only):\n{}",
            table(&["Rule", "Cases", "First", "Second", "Change"], totals)
        )
        .unwrap();
        writeln!(out, "{matched} matched cases. Ranges show trial variability, not confidence intervals. Changes are descriptive, not significance tests.").unwrap();
    } else {
        writeln!(out, "{matched} matched cases; {changed} changed. Negative differences mean smaller models. P nnz counts both symmetric halves.").unwrap();
    }
    for issue in &issues {
        writeln!(out, "{}", paint(issue, 33, color)).unwrap();
    }
    writeln!(out, "{} unmatched or invalid cases.", issues.len()).unwrap();
    (out, !issues.is_empty())
}

pub fn compare(old: &Run, new: &Run) -> bool {
    let color = std::env::var_os("NO_COLOR").is_none()
        && (std::io::stdout().is_terminal()
            || std::env::var("CARGO_TERM_COLOR").as_deref() == Ok("always"));
    let (text, invalid) = render(old, new, color);
    print!("{text}");
    invalid
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::results::{Counts, Metadata};
    use std::collections::BTreeMap;

    fn run(kind: Kind, times: &[u64]) -> Run {
        let counts = Counts {
            variables: 3,
            linear_rows: 2,
            conic_rows: 0,
            a_nonzeros: 4,
            g_nonzeros: 0,
            p_nonzeros: 0,
        };
        Run {
            metadata: Metadata {
                version: 1,
                name: "test".into(),
                kind,
                trials: times.len(),
                created_unix_seconds: 0,
                machine: "test".into(),
                settings: "same".into(),
            },
            cases: BTreeMap::from([(
                "netlib/AFIRO/all".into(),
                Case {
                    input_hash: "same".into(),
                    errors: vec![],
                    measurements: times
                        .iter()
                        .map(|&t| Measurement {
                            elapsed_ns: (kind == Kind::Time).then_some(t),
                            outcome: "reduced".into(),
                            before: counts.clone(),
                            after: Some(counts.clone()),
                            time_limit_reached: false,
                        })
                        .collect(),
                },
            )]),
            complete: true,
        }
    }

    #[test]
    fn timing_uses_independent_trials_and_reports_variability() {
        let a = run(Kind::Time, &[10, 20, 30]);
        let b = run(Kind::Time, &[5, 10, 15]);
        let (text, invalid) = render(&a, &b, false);
        assert!(!invalid);
        assert!(text.contains("-50.00%"));
        assert!(text.contains("20.0ns [12.0ns–28.0ns] n=3"));
        assert!(!text.contains('\x1b'));
    }

    #[test]
    fn size_keeps_mixed_changes_and_missing_sizes_visible() {
        let a = run(Kind::Size, &[0]);
        let mut b = run(Kind::Size, &[0]);
        let m = &mut b.cases.values_mut().next().unwrap().measurements[0];
        m.after.as_mut().unwrap().variables -= 1;
        m.after.as_mut().unwrap().p_nonzeros += 2;
        let (text, invalid) = render(&a, &b, true);
        assert!(!invalid);
        assert!(text.contains("\x1b[32m-1"));
        assert!(text.contains("\x1b[31m+2"));
        let m = &mut b.cases.values_mut().next().unwrap().measurements[0];
        m.after = None;
        m.outcome = "infeasible".into();
        let (text, _) = render(&a, &b, false);
        assert!(text.contains("unavailable"));
        assert!(text.contains("infeasible"));
    }

    #[test]
    fn incompatible_inputs_and_failed_or_missing_trials_are_flagged() {
        let a = run(Kind::Time, &[10, 20]);
        let mut b = run(Kind::Time, &[10, 20]);
        b.cases.values_mut().next().unwrap().input_hash = "changed".into();
        assert!(render(&a, &b, false).1);
        b.cases.values_mut().next().unwrap().input_hash = "same".into();
        b.cases.values_mut().next().unwrap().measurements.pop();
        assert!(render(&a, &b, false).1);
        b.cases.clear();
        assert!(render(&a, &b, false).1);
        b.metadata.settings = "different".into();
        assert!(render(&a, &b, false).0.contains("settings differ"));
    }
}
