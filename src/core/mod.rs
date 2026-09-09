// SPDX-License-Identifier: Apache-2.0
// Modified for this library; copyright and attribution notices are in NOTICE.
//! Working model and rule scheduling for convex quadratic and conic presolve.
//!
//! Sparse storage and model mutations support semantic rule families in
//! `rules`; `schedule` orders their execution; `crate::postsolve` reverses their edits.
//! Native dual signs satisfy:
//! `P x + c = A^T y + z`, with positive multipliers on lower bounds.

pub(crate) mod activity;
pub(crate) mod model;
pub(crate) mod objective;
pub(crate) mod queues;
mod rules;
pub(crate) mod schedule;
