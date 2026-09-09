//! Rule families organized by the model property they act on.
//!
//! Each family checks applicability and uses `Model` mutations to update sparse
//! storage, work queues, and recovery records together. The phase order and
//! cleanup loop live in `core::schedule`.

mod bounds;
mod cones;
mod dual_fixing;
mod parallel;
mod rows;
mod sparsification;
mod substitution;
mod variables;
