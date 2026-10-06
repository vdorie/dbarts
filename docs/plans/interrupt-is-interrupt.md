# interrupt-is-interrupt: an interrupted fit raises R's interrupt

Status: PLANNED (dec-B256).

agent: sonnet implementer, one; opus reviewer.
rng: NEUTRAL. No draw changes; only what is raised after a cancelled run.
window: pre-release (dec-B256).
budget: ~250 lines (bridge ~60, R ~20, tinytest ~90, manual and header comments ~50, NEWS none). Plans have
run 1.5-2x low.

## Goal

A run stopped by the user raises R's interrupt condition, as base R's long computations do, on every entry
that runs the sampler: nothing is printed as an error, `try()` and error handlers do not catch it, and
`tryCatch(interrupt = )` sees it. The sampler is left as it is left today: usable, at a valid point of its
chain, the interrupted call returning nothing.

## Context

- The bridge polls for an interrupt through one function
  ([`bartcore_bridge::userInterrupted`](../../src/R_interface_bartcore.cpp)), which runs
  `R_CheckUserInterrupt` under `R_ToplevelExec` so the worker threads can be joined before control leaves.
  That consumes the pending interrupt. Every run entry then raises the error "sampler run interrupted"
  (`Rf_error`): the three R run entries in src/R_interface_bartcore.cpp and the flat
  [`dbarts_sampler_run`](../../src/C_interface.cpp).
- Measured at the branch, a loop of three fits in `try()` interrupted during the first: the error is caught
  by `try()` and the other two fits run to their end. 0.9-34 could not be stopped until its run finished
  and then raised a true interrupt.
- A stop asked for by a draw or sweep callback is a different thing and keeps its own error or normal
  return.
- inst/tinytest/test-monotone.R pins the error text on both routes through the test hook that arms the
  poll; inst/include/dbarts/dbarts.h documents that the flat run "RAISES" the error.
- The sampler stops with a pass over the trees partly applied where the poll sits inside one (a slow
  monotone count, large logistic counts, a latent refresh) and is not rolled back (dec-B256).

## Constraints

- No draw changes: the seeded snapshot files and every equivalence scenario are identical.
- Only entry points documented as API are called from C (the R-devel check in CI reports non-API calls).
  One way that needs none: after the workers have joined and nothing of the library's is live, evaluate an
  unexported R function that signals a condition of class `c("interrupt", "condition")` with
  `signalCondition` and, if no handler takes it, invokes the `abort` restart. That is what an unhandled
  interrupt does. A direct call into R's own interrupt entry is acceptable only if it is API on every
  platform CI builds.
- The jump out happens where the error is raised today, after the same cleanup; a callback's stop and an
  engine error keep their paths and their precedence over an interrupt.
- The flat C header's signatures do not change; its comments say what a caller now sees.
- An entry that cannot be interrupted at all today is listed in the landing note and queued in TODO, not
  built here.

## Steps

1. Inventory: every place the package turns a cancelled run into an error, and every long-running entry
   (`bart`, `bartBT`, `dbarts` and `$run`, `$sampleTreesFromPrior`, `xbart`, `rbart_vi`, `predict` with a
   latent refresh, `pdbart`, the flat run) with whether an interrupt reaches it.
2. One bridge function that raises the interrupt, called in place of each "sampler run interrupted" error.
3. tinytest, in process, through the existing hook that arms the poll: the condition has class `interrupt`
   and not `error`; `try()` around an interrupted run does not return (an enclosing
   `tryCatch(interrupt = )` receives it); calling handlers and `on.exit` run; the sampler runs again
   afterwards; a callback's stop is still what it was. The same on the flat route through the compiled
   consumer. The two pins of the error text are rewritten.
4. A subprocess check kept out of the CRAN tests (at home only): a real SIGINT sent to a loop of fits in
   `try()` stops the loop.
5. Manual: `$run` and the embedding page say what an interrupt does and leaves behind, including the
   partly applied pass; monotone.Rd's sentence and the flat header's comment are brought in line.

## Verification

- `R CMD INSTALL --preclean` into the slice's own library; full tinytest; `tests/cpp` builds and passes.
- The four seeded snapshot files pass unchanged on a reference build; the equivalence compares are
  identical.
- stan4bart's suite on a fresh install against this build (it calls the flat run).
- `lintr::lint_package()`, `air format --check .`, the three tools/ checks and
  `Rscript benchmarks/R/mutation-battery.R verify-anchors` clean; `R CMD check --as-cran` shows no new
  note.
