# xbart-fork-workers: forked workers and a user cluster for xbart

agent: sonnet
rng: neutral (no draw moves; every unit already seeds from its own stream)
budget: ~250 lines (R ~80, tests ~120, man and NEWS ~50)

Status: PLANNED 2026-09-28

## Goal

`xbart` can run its (replication, fold) units on forked workers as well as
on separate R sessions, chooses between them by a reasonable default, and
accepts a cluster the caller already started (dec-B131).

## Context

- Today every `n.threads > 1` call starts `parallel::makeCluster(numChunks)`
  socket workers in `xbart`'s `runUnits` closure and stops them on exit. Each
  worker loads R and dbarts and receives a copy of `spec`: about 80 MB per
  worker beyond the data, plus startup time (measured 2026-09-28, dec-A05).
- Precedent: `boot::boot(parallel = c("no", "multicore", "snow"), ncpus, cl)`
  with its default read from `getOption("boot.parallel")`; `future` refuses
  multicore on Windows and in RStudio.
- Seeds: a unit's sampler seeds come from the control's seed slot and its
  split from `splitSeeds`, both drawn in the calling process, so results do
  not depend on the worker kind.

## Surface

- `parallel = getOption("dbarts.parallel", "auto")`, one of `"auto"`,
  `"fork"`, `"socket"` (match.arg semantics against those three).
- `cl = NULL`: a cluster from the parallel package (inherits `"cluster"`).
- `"auto"` forks unless the platform is Windows, the session runs under
  RStudio (`RSTUDIO` environment variable `"1"`) or Positron (`POSITRON`
  `"1"`), or `.Platform$GUI == "AQUA"` (R.app); otherwise sockets.
- `"fork"` on Windows is an error naming `"socket"`.
- A supplied `cl` is used whatever `parallel` says; units are split into
  `min(n.threads, length(cl))` chunks; `xbart` never stops a cluster it did
  not start. The workers must be able to load dbarts; a worker error
  propagates as it does today.
- At one chunk nothing changes: units run in the calling session.

## Steps

1. Add the two formals after `n.threads`; validate `cl` (NULL or inherits
   "cluster", non-empty) and `parallel` early with the other argument checks.
2. In `runUnits`, dispatch: `cl` -> `parallel::clusterMap` on it; fork ->
   `parallel::mclapply` over the chunk list with `mc.cores = numChunks`,
   `mc.preschedule = FALSE`, re-raising the first `try-error` result as an
   error in the caller with its message; socket -> today's path.
3. The verbose line names the worker kind.
4. man/xbart.Rd: document both arguments and the option in `\arguments`,
   including the fork hazard in plain words (a hang or crash when the calling
   session has already used a library that is not fork-safe, such as a
   multithreaded BLAS; set `parallel = "socket"` or
   `options(dbarts.parallel = "socket")`), and move the per-worker memory
   sentence under `parallel` since it applies to sockets only.
5. inst/NEWS.Rd 1.0-0: one item.
6. Tests (new file inst/tinytest/test-xbart-workers.R): fork, socket and a
   user cluster give identical results to `n.threads = 1` on a small k-fold
   call with a seed (fork arm skipped on Windows); `"fork"` on Windows
   errors (only runnable there; guard); a supplied `cl` is left running
   after the call; invalid `parallel` and `cl` are refused; an error raised
   in a forked worker (a loss function that stops) reaches the caller with
   its message; the option sets the default.

## Verification

- Full tinytest; lint gates per CLAUDE.local.md; R CMD check --as-cran
  on the built tarball (no new NOTE for the option or parallel usage).
- Equivalence is not affected (neutral); no exact gate reads xbart workers.

## Agent-made calls

The option name `dbarts.parallel`; the auto rule's GUI list; `cl` taking
precedence over `parallel` and the `min(n.threads, length(cl))` chunk
count; an error, not a silent fallback, for `"fork"` on Windows.
