//! Per-thread heap allocation counts.
//!
//! [`Counting`] wraps the system allocator and counts on the calling thread,
//! so a measurement around serial code is exact and independent of what
//! other threads do. A binary opts in with
//! `#[global_allocator] static A: Counting = Counting;`; without it every
//! measurement is zero.

use serde::Serialize;
use std::alloc::{GlobalAlloc, Layout, System};
use std::cell::Cell;

/// Allocation counts on one thread over a measured region.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize)]
pub struct Allocations {
    /// Calls to `alloc`, `alloc_zeroed`, and `realloc`.
    pub count: u64,
    /// Bytes requested by those calls, counting a `realloc` at its new size.
    pub bytes: u64,
    /// The maximum of bytes allocated minus bytes freed since the region
    /// began. Freeing memory allocated before the region counts against it.
    pub peak_live_bytes: u64,
}

#[derive(Clone, Copy)]
struct Counters {
    count: u64,
    bytes: u64,
    live: i64,
    peak: i64,
}
const ZERO: Counters = Counters {
    count: 0,
    bytes: 0,
    live: 0,
    peak: 0,
};

thread_local! {
    // A const-initialized `Cell` of plain data needs neither allocation nor
    // a destructor, so the allocator can use it at any point in a thread's life.
    static COUNTERS: Cell<Counters> = const { Cell::new(ZERO) };
}

fn record(allocated: usize, freed: usize, calls: u64) {
    let _ = COUNTERS.try_with(|cell| {
        let mut c = cell.get();
        c.count += calls;
        c.bytes += allocated as u64;
        c.live += allocated as i64 - freed as i64;
        c.peak = c.peak.max(c.live);
        cell.set(c);
    });
}

/// The system allocator, counting on the calling thread.
pub struct Counting;

unsafe impl GlobalAlloc for Counting {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        let ptr = unsafe { System.alloc(layout) };
        if !ptr.is_null() {
            record(layout.size(), 0, 1);
        }
        ptr
    }

    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        let ptr = unsafe { System.alloc_zeroed(layout) };
        if !ptr.is_null() {
            record(layout.size(), 0, 1);
        }
        ptr
    }

    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        let new = unsafe { System.realloc(ptr, layout, new_size) };
        if !new.is_null() {
            record(new_size, layout.size(), 1);
        }
        new
    }

    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        unsafe { System.dealloc(ptr, layout) };
        record(0, layout.size(), 0);
    }
}

/// Run `f` and return its result with the allocations it made on this
/// thread. Measurements do not nest.
pub fn measure<R>(f: impl FnOnce() -> R) -> (R, Allocations) {
    COUNTERS.with(|c| c.set(ZERO));
    let result = f();
    let c = COUNTERS.with(Cell::get);
    let allocations = Allocations {
        count: c.count,
        bytes: c.bytes,
        peak_live_bytes: c.peak as u64,
    };
    (result, allocations)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[global_allocator]
    static ALLOCATOR: Counting = Counting;

    #[test]
    fn counts_allocations_bytes_and_peak() {
        let ((), a) = measure(|| {
            let mut v: Vec<u8> = Vec::with_capacity(100);
            v.extend_from_slice(&[1; 100]);
            let w = vec![0u64; 50];
            drop(v);
            let mut s = Vec::<u8>::with_capacity(10);
            s.reserve_exact(30);
            std::hint::black_box((&w, &s));
        });
        // 100 + 400 live at the peak; the realloc to 30 comes after 100 is freed.
        assert_eq!(
            a,
            Allocations {
                count: 4,
                bytes: 100 + 400 + 10 + 30,
                peak_live_bytes: 500,
            }
        );
    }

    #[test]
    fn freeing_earlier_memory_does_not_raise_the_peak() {
        let before = vec![0u8; 1000];
        let (n, a) = measure(|| {
            drop(before);
            let v = vec![1u8; 600];
            v.len()
        });
        assert_eq!(n, 600);
        assert_eq!(
            a,
            Allocations {
                count: 1,
                bytes: 600,
                peak_live_bytes: 0,
            }
        );
    }

    #[test]
    fn other_threads_are_not_counted() {
        let ((), a) = measure(|| {
            std::thread::scope(|s| {
                s.spawn(|| std::hint::black_box(vec![0u8; 1 << 20]));
            });
        });
        assert!(a.bytes < 1 << 20, "{a:?}");
    }
}
