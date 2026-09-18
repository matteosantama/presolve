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
    /// `key` must order values like their `Ord` prefix does, so that a radix
    /// pass on it followed by a comparison within equal keys gives the same
    /// order as sorting the values directly.
    pub fn filter_map_sorted<T: Ord + Copy + Send>(
        &self,
        count: usize,
        enough_work: bool,
        map: impl Fn(usize) -> Option<T> + Sync,
        key: impl Fn(&T) -> u64,
    ) -> Vec<T> {
        if enough_work && let Self::Parallel(pool) = self {
            return pool.install(|| {
                let mut values: Vec<_> = (0..count).into_par_iter().filter_map(&map).collect();
                values.par_sort_unstable();
                values
            });
        }
        let mut values: Vec<_> = (0..count).filter_map(map).collect();
        sort_by_radix_key(&mut values, key);
        values
    }
}

/// Order distinct values by a 64-bit key prefix with four counting passes,
/// then finish each run of equal keys by comparison. Fingerprint keys are
/// spread over the whole range, so this beats a comparison sort on models
/// with many rows or columns. Each pass sweeps 65536 buckets, so inputs
/// below the threshold use the library sort instead.
const RADIX_THRESHOLD: usize = 32_768;
fn sort_by_radix_key<T: Ord + Copy>(values: &mut Vec<T>, key: impl Fn(&T) -> u64) {
    let n = values.len();
    if n < RADIX_THRESHOLD {
        values.sort_unstable();
        return;
    }
    let mut keyed: Vec<(u64, u32)> = values
        .iter()
        .enumerate()
        .map(|(at, value)| (key(value), at as u32))
        .collect();
    let mut scratch = vec![(0u64, 0u32); n];
    for shift in (0..64).step_by(16) {
        let digit = |k: u64| ((k >> shift) & 0xffff) as usize;
        let mut starts = vec![0usize; 1 << 16];
        for &(k, _) in &keyed {
            starts[digit(k)] += 1;
        }
        // Skip a pass whose digit is constant, which is common in the high
        // half of a key built from a single 32-bit hash.
        if starts.contains(&n) {
            continue;
        }
        let mut total = 0;
        for start in &mut starts {
            let count = *start;
            *start = total;
            total += count;
        }
        for &entry in &keyed {
            let slot = &mut starts[digit(entry.0)];
            scratch[*slot] = entry;
            *slot += 1;
        }
        std::mem::swap(&mut keyed, &mut scratch);
    }
    let mut sorted: Vec<T> = keyed.iter().map(|&(_, at)| values[at as usize]).collect();
    let mut start = 0;
    while start < n {
        let mut end = start + 1;
        while end < n && keyed[end].0 == keyed[start].0 {
            end += 1;
        }
        sorted[start..end].sort_unstable();
        start = end;
    }
    *values = sorted;
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
                let values = executor.filter_map_sorted(
                    2048,
                    enough_work,
                    |i| {
                        assert_eq!(std::thread::current().id(), caller);
                        (i % 2 == 0).then_some(2048 - i)
                    },
                    |&v| v as u64,
                );
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
            let values = executor.filter_map_sorted(
                4096,
                true,
                |i| {
                    assert_eq!(rayon::current_num_threads(), executor.threads());
                    Some(4095 - i)
                },
                |&v| v as u64,
            );
            assert_eq!(values, (0..4096).collect::<Vec<_>>());
        }
        assert_eq!(before, ids(&two));
    }

    #[test]
    fn radix_ordering_matches_comparison_ordering_with_duplicate_prefixes() {
        type Candidate = ((u32, u32), Option<(u32, u32)>, usize);
        // Keys collide on the radix prefix; only the full value breaks ties.
        let mut random = 0x9e37_79b9_7f4a_7c15u64;
        let values: Vec<Candidate> = (0..2 * RADIX_THRESHOLD)
            .map(|i| {
                random ^= random << 13;
                random ^= random >> 7;
                random ^= random << 17;
                let hash = if i % 7 == 0 {
                    (3, 5)
                } else {
                    ((random >> 32) as u32 % 4, random as u32 % 512)
                };
                let extra = (i % 3 != 0).then_some((i as u32 % 5, 0));
                (hash, extra, i)
            })
            .collect();
        let mut expected = values.clone();
        expected.sort_unstable();
        for count in [RADIX_THRESHOLD - 1, RADIX_THRESHOLD, values.len()] {
            let sorted = Executor::Serial.filter_map_sorted(
                count,
                false,
                |i| Some(values[i]),
                |&(hash, _, _)| (u64::from(hash.0) << 32) | u64::from(hash.1),
            );
            let mut reference = values[..count].to_vec();
            reference.sort_unstable();
            assert_eq!(sorted, reference);
        }
        assert_eq!(expected.len(), values.len());
    }
}
