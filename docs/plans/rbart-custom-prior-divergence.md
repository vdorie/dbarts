# rbart-custom-prior-divergence

agent: opus
rng: neutral while investigating; the fix changes draws on the
     custom-prior path only (bug fix; the in-core path is the oracle)
budget: reproduction + diagnosis memo ~50 lines; fix + test ~150

## Goal

The rbart_vi custom-prior R loop and the in-core built-in-prior path
draw from the same posterior again: the 15-25x tau divergence at
n >~ 300 found during rbart-loop-profile (2026-07-07) is explained,
fixed, and regression-tested.

## Context

- Evidence (rbart-loop-profile.md, Measurements): with the custom
  prior verbatim rbart.priors$cauchy (so the R loop runs but the model
  is identical to prior = cauchy), tau settles ~15-25x above the
  in-core path's and the true value once n exceeds a few hundred;
  below that the paths agree. Both arms share creation semantics.
- The R loop: R/rbart.R rbart_vi_run (posteriorClosure/evalEnv from
  the prior builder around R/rbart.R:442-483); per-iteration blocks
  are the ranef draw, offset refresh, and tau slice sample (with an
  optim mode search). rel.scale = sd(y) at R/rbart.R:716.
- The in-core path: the decorator installed when the prior symbol
  matches rbart.priors (R/rbart.R:266-290, rel.scale at :272), running
  the same blocks in C.
- Suspects, in likelihood order: a scale mismatch between R-loop
  quantities and the engine's internal response scaling (the rewrite
  moved y-scaling in-core; original-vs-internal sigma or offset units
  in the ranef/tau conditionals); the slice sampler or its mode search
  operating on the wrong scale; weights/latents handling drift between
  the two implementations.

## Constraints

- Diagnose before fixing: a short memo (root cause, the mismatched
  quantity, the proposed fix) goes to VD before product code changes
  if the fix would change documented behavior; a plain bug fix with
  the in-core path as oracle proceeds directly.
- The in-core path must remain bitwise-untouched (it is the
  equivalence-gated oracle).
- Out of scope: the loop's bridge-overhead fix (rides capi-callbacks);
  any prior interface change.

## Steps

1. Minimal seeded reproduction: same data, both paths, tau traces at
   n in {200, 1e3, 1e4}; confirm the divergence and its onset.
2. Fixed-state conditional comparison: freeze the forest state, run
   one ranef/tau Gibbs block through each path, compare the
   conditionals directly (the in-core decorator computes the same
   quantities in C - instrument or replicate them in R).
3. Root cause + fix on the R-loop side only.
4. Regression test: a tinytest asserting statistical agreement of tau
   posteriors between paths (means within MC error on a small seeded
   problem); check whether any equivalence scenario exercises the
   custom-prior loop (none is expected - the baseline covers built-in
   priors) and report coverage.

## Verification

- The step-1 script shows divergence before, agreement after.
- Full tinytest green; equivalence exact 18/18 (in-core path
  untouched); the new path-agreement test in the suite.
