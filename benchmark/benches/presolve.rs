//! Instruction counts of presolve under Callgrind on a fixed subset of the
//! corpus, for both profiles. Loading happens in setup and dropping the result
//! after the benchmark function, so only the presolve call is measured.
//!
//! The subset covers every rule and reduction kind the corpus triggers under
//! either profile, with at least one instance from each family and a spread
//! of sizes, and presolves in about 1 s natively for both profiles together.

use benchmark::mps;
use benchmark::run::Profile;
use gungraun::{
    Callgrind, LibraryBenchmarkConfig, library_benchmark, library_benchmark_group, main,
};
use presolve::result::PresolveResult;
use presolve::{Presolver, Problem};
use std::hint::black_box;
use std::path::Path;

fn load(instance: &str, profile: Profile) -> (Presolver, Problem) {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("data")
        .join(format!("{instance}.mps.gz"));
    let problem = mps::read(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    (Presolver::new(profile.settings()).unwrap(), problem)
}

fn load_default(instance: &str) -> (Presolver, Problem) {
    load(instance, Profile::Default)
}

fn load_aggressive(instance: &str) -> (Presolver, Problem) {
    load(instance, Profile::Aggressive)
}

// Both functions list the same instances; keep the lists identical.
#[library_benchmark(setup = load_default)]
#[bench::ken_11("kennington/ken-11")]
#[bench::cre_b("kennington/cre-b")]
#[bench::boyd1("maros-meszaros/BOYD1")]
#[bench::cvxqp1_l("maros-meszaros/CVXQP1_L")]
#[bench::q25fv47("maros-meszaros/Q25FV47")]
#[bench::academictimetablesmall("miplib/academictimetablesmall")]
#[bench::brazil3("miplib/brazil3")]
#[bench::irish_electricity("miplib/irish-electricity")]
#[bench::d2q06c("netlib/D2Q06C")]
#[bench::greenbea("netlib/GREENBEA")]
#[bench::qplib_8785("qplib/QPLIB_8785")]
#[bench::qplib_8790("qplib/QPLIB_8790")]
#[bench::qplib_8845("qplib/QPLIB_8845")]
fn default((presolver, problem): (Presolver, Problem)) -> PresolveResult {
    black_box(presolver.presolve(black_box(problem)))
}

#[library_benchmark(setup = load_aggressive)]
#[bench::ken_11("kennington/ken-11")]
#[bench::cre_b("kennington/cre-b")]
#[bench::boyd1("maros-meszaros/BOYD1")]
#[bench::cvxqp1_l("maros-meszaros/CVXQP1_L")]
#[bench::q25fv47("maros-meszaros/Q25FV47")]
#[bench::academictimetablesmall("miplib/academictimetablesmall")]
#[bench::brazil3("miplib/brazil3")]
#[bench::irish_electricity("miplib/irish-electricity")]
#[bench::d2q06c("netlib/D2Q06C")]
#[bench::greenbea("netlib/GREENBEA")]
#[bench::qplib_8785("qplib/QPLIB_8785")]
#[bench::qplib_8790("qplib/QPLIB_8790")]
#[bench::qplib_8845("qplib/QPLIB_8845")]
fn aggressive((presolver, problem): (Presolver, Problem)) -> PresolveResult {
    black_box(presolver.presolve(black_box(problem)))
}

library_benchmark_group!(name = profiles, benchmarks = [default, aggressive]);

main!(
    config = LibraryBenchmarkConfig::default().tool(Callgrind::with_args(["--cache-sim=no"])),
    library_benchmark_groups = profiles
);
