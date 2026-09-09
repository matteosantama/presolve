//! Serial and parallel execution of independent discovery work.
use rayon::{ThreadPool, ThreadPoolBuildError, ThreadPoolBuilder, prelude::*};

#[derive(Debug)]
pub(crate) enum Executor {
    Serial,
    Parallel(ThreadPool),
}

impl Executor {
    pub fn new(threads: usize) -> Result<Self, ThreadPoolBuildError> {
        if threads == 1 {
            return Ok(Self::Serial);
        }
        // Rayon owns automatic selection, including RAYON_NUM_THREADS.
        let pool = ThreadPoolBuilder::new().num_threads(threads).build()?;
        Ok(if pool.current_num_threads() == 1 {
            Self::Serial
        } else {
            Self::Parallel(pool)
        })
    }

    pub fn threads(&self) -> usize {
        match self {
            Self::Serial => 1,
            Self::Parallel(pool) => pool.current_num_threads(),
        }
    }

    /// Discover independent candidates and sort them in a reproducible order.
    /// Callers decide whether their workload is large enough for parallel work.
    /// The serial branch uses ordinary iterators, even inside another Rayon pool.
    pub fn filter_map_sorted<T: Ord + Send>(
        &self,
        count: usize,
        enough_work: bool,
        map: impl Fn(usize) -> Option<T> + Sync,
    ) -> Vec<T> {
        if enough_work {
            if let Self::Parallel(pool) = self {
                return pool.install(|| {
                    let mut values: Vec<_> = (0..count).into_par_iter().filter_map(&map).collect();
                    values.par_sort_unstable();
                    values
                });
            }
        }
        let mut values: Vec<_> = (0..count).filter_map(map).collect();
        values.sort_unstable();
        values
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn serial_and_small_work_stay_on_the_caller_inside_another_pool() {
        let parallel = Executor::new(3).unwrap();
        let caller_pool = ThreadPoolBuilder::new().num_threads(2).build().unwrap();
        caller_pool.install(|| {
            let caller = std::thread::current().id();
            for (executor, enough_work) in [(&Executor::Serial, true), (&parallel, false)] {
                let values = executor.filter_map_sorted(2048, enough_work, |i| {
                    assert_eq!(std::thread::current().id(), caller);
                    (i % 2 == 0).then_some(2048 - i)
                });
                assert_eq!(values, (1..=1024).map(|i| 2 * i).collect::<Vec<_>>());
            }
        });
    }

    #[test]
    fn executors_keep_their_own_threads_across_repeated_calls() {
        let two = Executor::new(2).unwrap();
        let four = Executor::new(4).unwrap();
        let ids = |executor: &Executor| match executor {
            Executor::Parallel(pool) => pool.broadcast(|_| std::thread::current().id()),
            Executor::Serial => unreachable!(),
        };
        let before = ids(&two);
        for executor in [&two, &four, &two] {
            let values = executor.filter_map_sorted(4096, true, |i| {
                assert_eq!(rayon::current_num_threads(), executor.threads());
                Some(4095 - i)
            });
            assert_eq!(values, (0..4096).collect::<Vec<_>>());
        }
        assert_eq!(before, ids(&two));
    }
}
