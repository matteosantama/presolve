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

/// Order distinct values by a 64-bit key prefix, then finish each run of
/// equal keys by comparison. Fingerprint keys are usually spread over the
/// high bits but not always: on ±1 matrices the coefficient word is nearly
/// constant, and on problems with few distinct supports only the low word
/// discriminates. An MSD radix with a constant-digit skip handles the first
/// shape; the second falls back to four LSD passes, and already sorted or
/// small inputs use the library sort.
const COMPARISON_THRESHOLD: usize = 4096;
fn sort_by_radix_key<T: Ord + Copy>(values: &mut Vec<T>, key: impl Fn(&T) -> u64) {
    let n = values.len();
    if n < COMPARISON_THRESHOLD {
        values.sort_unstable();
        return;
    }
    if values.is_sorted() {
        return;
    }
    // Sort records directly, avoiding a separate key/index array and final gather.
    let mut scratch = values.clone();
    let mut seen = vec![false; 1 << 16];
    let mut distinct = 0;
    for value in values.iter() {
        let digit = ((key(value) >> 32) & 0xffff) as usize;
        if !seen[digit] {
            seen[digit] = true;
            distinct += 1;
        }
    }
    if distinct < n / 64 {
        lsd_sort(values, &mut scratch, &key);
    } else {
        msd_sort(values, &mut scratch, 0, &key);
    }
    let mut start = 0;
    while start < n {
        let mut end = start + 1;
        while end < n && key(&values[end]) == key(&values[start]) {
            end += 1;
        }
        values[start..end].sort_unstable();
        start = end;
    }
}

/// Most-significant-digit radix sort by key. The digit width follows the
/// slice length so buckets stay a few elements deep; a digit that is constant
/// over the slice is consumed without a scatter.
fn msd_sort<T: Ord + Copy>(
    values: &mut [T],
    scratch: &mut [T],
    consumed: u32,
    key: &impl Fn(&T) -> u64,
) {
    let n = values.len();
    if n <= 64 || consumed >= 64 {
        values.sort_unstable();
        return;
    }
    let bits = ((n / 4).ilog2()).clamp(4, 16).min(64 - consumed);
    let shift = 64 - consumed - bits;
    let mask = (1usize << bits) - 1;
    let digit = |k: u64| ((k >> shift) as usize) & mask;
    let mut starts = vec![0u32; (1 << bits) + 1];
    for entry in values.iter() {
        starts[digit(key(entry)) + 1] += 1;
    }
    if starts[1..].contains(&(n as u32)) {
        return msd_sort(values, scratch, consumed + bits, key);
    }
    for b in 0..(1 << bits) {
        starts[b + 1] += starts[b];
    }
    let mut next = starts.clone();
    for &entry in values.iter() {
        let slot = &mut next[digit(key(&entry))];
        scratch[*slot as usize] = entry;
        *slot += 1;
    }
    values.copy_from_slice(&scratch[..n]);
    for b in 0..(1 << bits) {
        let (lo, hi) = (starts[b] as usize, starts[b + 1] as usize);
        if hi - lo > 1 {
            msd_sort(
                &mut values[lo..hi],
                &mut scratch[lo..hi],
                consumed + bits,
                key,
            );
        }
    }
}

/// Four 16-bit counting passes from the low word up, skipping constant digits.
fn lsd_sort<T: Copy>(keyed: &mut Vec<T>, scratch: &mut Vec<T>, key: &impl Fn(&T) -> u64) {
    let n = keyed.len();
    for shift in (0..64).step_by(16) {
        let digit = |k: u64| ((k >> shift) & 0xffff) as usize;
        let mut starts = vec![0usize; 1 << 16];
        for entry in keyed.iter() {
            starts[digit(key(entry))] += 1;
        }
        if starts.contains(&n) {
            continue;
        }
        let mut total = 0;
        for start in &mut starts {
            let count = *start;
            *start = total;
            total += count;
        }
        for &entry in keyed.iter() {
            let slot = &mut starts[digit(key(&entry))];
            scratch[*slot] = entry;
            *slot += 1;
        }
        std::mem::swap(keyed, scratch);
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
        let values: Vec<Candidate> = (0..16 * COMPARISON_THRESHOLD)
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
        for count in [COMPARISON_THRESHOLD - 1, COMPARISON_THRESHOLD, values.len()] {
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
    }

    #[test]
    fn radix_sort_handles_sorted_skewed_and_spread_keys() {
        let n = 50_000usize;
        let mut random = 0x2545_f491_4f6c_dd1du64;
        let mut next = move || {
            random ^= random << 13;
            random ^= random >> 7;
            random ^= random << 17;
            random
        };
        // Already sorted keys, keys with few distinct high words (only the low
        // word discriminates), keys with a near-constant low word, and spread
        // keys, each with duplicates.
        let sorted: Vec<(u64, usize)> = (0..n).map(|i| ((i as u64 / 3) << 20, i)).collect();
        let skewed: Vec<(u64, usize)> = (0..n)
            .map(|i| (((next() % 8) << 32) | (next() & 0xffff_ffff), i))
            .collect();
        let constant_low: Vec<(u64, usize)> = (0..n).map(|i| ((next() << 32) | 7, i)).collect();
        let spread: Vec<(u64, usize)> = (0..n).map(|i| (next() % 1000 * 977, i)).collect();
        for values in [sorted, skewed, constant_low, spread] {
            let sorted =
                Executor::Serial.filter_map_sorted(n, false, |i| Some(values[i]), |&(k, _)| k);
            let mut reference = values.clone();
            reference.sort_unstable();
            assert_eq!(sorted, reference);
        }
    }
}
