# TODO

Finished jobs are removed from this file: their rules are in `docs/refinement.md` and `docs/differences.md`, and
their measurements in `docs/testing.md`. Jobs that need the user's approval say so: don't start them before it is
given.

## Job 18: Benchmark with a prior taxonomy and score the curator workload

**Status: not started. Confirm with the user before starting.**

### Problem

- `scripts/benchmark.sh` runs curatomatic without `--prior-taxonomy`, so it measures the behaviour without downward
  refinement (Jobs 13 and 15 are not measured).
- It compares against the stored Python outputs, whose r226 runs had no prior, so Python is not given the same input
  as Rust.
- The curator workload, the user's metric of success, is scored by `wip/scripts/workload.py` and
  `wip/scripts/summary.py`, which are not committed (`wip/` is git-ignored). The Python metadata with the prior is
  made by `wip/scripts/make_py_meta.py`, also not committed.

The benchmark was run by hand on 2026-11-23 (`docs/testing.md`, "Curator workload at `3dc483a`", steps in "Running
it"). This job scripts it.

### Tasks

1. Take the prior files as arguments (the official GTDB releases joined across ar53 and bac120: r220 for the r226
   datasets, r226 for the r232 datasets, see `docs/testing.md`, "Prior taxonomy") and pass them with
   `--prior-taxonomy`.
2. Move `make_py_meta.py` into `scripts/`, and re-run Python (`main`) with the same prior instead of linking the
   stored outputs (r232_bac120 takes about 40 min).
3. Move the workload scorer into `scripts/`, run it in `benchmark.sh`, and document its metrics in
   `docs/testing.md`.
4. Check that `benchmark.sh` gives the numbers in `docs/testing.md`, "Curator workload at `3dc483a`".

## Open questions (ask the user before starting)

- Moving down (Job 13) adds work as well as saving it: the operations added against the input went from 764 at
  `32d3350` to 1,897 with Jobs 10 to 13 (Python adds 1,814), and the taxon specificity in r226_bac120 dropped from
  99.95 to 99.56. Should the added operations be looked at?
- The effect of Job 10 was only measured without a prior (+3 operations, r226_bac120 only). It was not measured
  separately with Job 13.
