use presolve::{
    Outcome, Presolver, Settings,
    matrix::CscMatrix,
    postsolve::Solution,
    problem::{Bounds, Constraint, Problem},
    settings::Rules,
};
#[test]
fn final_scan_merges_rows_exposed_by_dual_fixing() {
    for (parallel, final_scan) in [(false, true), (true, false), (true, true)] {
        let p = Problem {
            p: None,
            c: vec![0., 0., 1.],
            c0: 0.,
            a: CscMatrix::from_triplets(
                2,
                3,
                vec![0, 0, 1, 1, 1],
                vec![0, 1, 0, 1, 2],
                vec![1.; 5],
            )
            .unwrap(),
            rows: vec![
                Constraint::Linear(Bounds {
                    lower: f64::NEG_INFINITY,
                    upper: 4.
                });
                2
            ],
            variable_bounds: vec![
                Bounds::FREE,
                Bounds::FREE,
                Bounds {
                    lower: 0.,
                    upper: f64::INFINITY,
                },
            ],
            cones: vec![],
        };
        let s = Settings {
            final_parallel_scan: final_scan,
            rules: Rules {
                parallel_rows: parallel,
                dual_propagation: true,
                ..Rules::none()
            },
            ..Settings::default()
        };
        let r = Presolver::new(s).unwrap().presolve(p);
        let size = r.stats.after.unwrap();
        assert_eq!(size.variables, 2);
        assert_eq!(size.linear_rows, if parallel && final_scan { 1 } else { 2 });
        let Outcome::Reduced(r) = r.outcome else {
            panic!("expected reduction")
        };
        let point = Solution {
            x: vec![0.; 2],
            y: vec![0.; size.linear_rows],
            z: vec![0.; 2],
            conic_dual: vec![],
            conic_slack: vec![],
        };
        let lifted = r.postsolve.recover_solution(point.as_ref());
        assert_eq!(lifted.x, vec![0.; 3]);
        assert_eq!(lifted.y, vec![0.; 2]);
        assert_eq!(lifted.z, vec![0., 0., 1.]);
    }
}
