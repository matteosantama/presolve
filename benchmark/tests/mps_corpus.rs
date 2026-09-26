//! The reader against the benchmark corpus.
use benchmark::{corpus, mps};
use std::path::{Path, PathBuf};

fn data() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("data")
}

/// Rows excluding the objective, columns, and nonzeros including the
/// objective row, as listed in the Netlib LP index.
#[test]
fn netlib_dimensions_match_the_published_index() {
    for (name, rows, columns, nonzeros) in [
        ("AFIRO", 27, 32, 88),
        ("BLEND", 74, 83, 521),
        // Blank RHS and bound set names.
        ("SIERRA", 1227, 2036, 9252),
        // Fixed format with spaces inside names.
        ("FORPLAN", 161, 421, 4916),
    ] {
        let problem = mps::read(&data().join(format!("netlib/{name}.mps.gz"))).unwrap();
        let objective = problem.c.iter().filter(|&&c| c != 0.0).count();
        assert_eq!(
            (problem.row_count(), problem.variable_count()),
            (rows, columns),
            "{name}"
        );
        assert_eq!(
            problem.a.values().len() + objective,
            nonzeros,
            "{name} nonzeros"
        );
    }
}

/// QFORPLAN is FORPLAN with a quadratic objective, in the same fixed format.
#[test]
fn quadratic_fixed_format_variant_shares_its_linear_part() {
    let lp = mps::read(&data().join("netlib/FORPLAN.mps.gz")).unwrap();
    let qp = mps::read(&data().join("maros-meszaros/QFORPLAN.mps.gz")).unwrap();
    assert_eq!(qp.a, lp.a);
    assert_eq!(qp.rows, lp.rows);
    assert!(qp.p.is_some_and(|p| !p.values().is_empty()));
}

/// Every discovered instance parses. Slow: run with `--ignored`.
#[test]
#[ignore]
fn every_corpus_file_parses() {
    let instances = corpus::discover(&data(), &[]).unwrap();
    assert!(!instances.is_empty());
    let failures: Vec<String> = instances
        .iter()
        .filter_map(|instance| {
            mps::read(&instance.path)
                .err()
                .map(|e| format!("{}: {e}", instance.id()))
        })
        .collect();
    assert!(failures.is_empty(), "{}", failures.join("\n"));
    let mut families: Vec<&str> = instances.iter().map(|i| i.family.as_str()).collect();
    families.dedup();
    eprintln!("parsed {} instances in {families:?}", instances.len());
}
