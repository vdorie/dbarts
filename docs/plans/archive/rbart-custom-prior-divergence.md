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
  the prior builder around [[R/rbart.R:442-483@95fa3899]]); per-iteration blocks
  are the ranef draw, offset refresh, and tau slice sample (with an
  optim mode search). rel.scale = sd(y) at [[R/rbart.R:716@95fa3899]].
- The in-core path: the decorator installed when the prior symbol
  matches rbart.priors ([[R/rbart.R:266-290@95fa3899]], rel.scale at [[R/rbart.R:272@95fa3899]]), running
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

## Findings (2026-07-07)

Root cause: R/rbart.R rbart_vi_fit, the first random-effects draw
conditioned on treeFit.train initialized from sampler$predict() on the
prior-sampled forest (former lines 772 and 800). A prior fit has ~zero
internal leaf values, so on the original scale it returns the response
midpoint (range/2 + min, the gaussian centering shift), not a fitted
value near the response mean. The reproduction only fired with normal
predictors: the Friedman f's 20*(x3-0.5)^2 term makes y strongly
right-skewed, so mean(y) and midpoint(y) diverge (e.g. 25.8 vs 143.1,
sd(y)=38.7). With uniform predictors (mean ~= midpoint) both paths agree
at every size, which is why the earlier small-n checks looked fine.

Mismatched quantity: the initial treeFit.train. resid = y.st -
treeFit.train inherited the ~constant (mean - midpoint) offset; the
first ranef block dumped it into every intercept, and thereafter it is a
stable wrong fixed point - the forest re-fits y - b (train ~= y) so
resid ~= b reproduces the same b, and the diffuse cauchy tau prior
(scale = 2.5*sd(y)) never pulls b back. tau settles near sd(y) instead
of the truth (~1). The in-core oracle (GroupedResponse) is immune
because it sweeps the forest before its first drawGroupEffects, so the
mean lands in the forest fit, not the intercepts.

Fix (R-loop side only, a plain bug fix - the interface contract is
unchanged): initialize treeFit.train from one fitted sweep
(sampler$run(0L, 1L)) instead of predict(), at both the n.burn>0 and
n.burn==0 init sites; the latter drops keepTrees for that throwaway
sweep. R/rbart.R +16 -7; new inst/tinytest/test-rbart-custom-prior.R
(skewed 8-group problem, asserts custom-loop tau ~= in-core tau).

Evidence (repro2.R, normal predictors, n.burn=n.samples=100, seed 1),
mean tau, custom loop (R) vs in-core (C):

  n=1e3 groups=10   before: tauR=94.8 tauC=1.13 (84x)
                    after:  tauR=0.95 tauC=1.13 (0.84x)
  n=1e4 groups=100  before: tauR=100.3 tauC=1.11 (91x)
                    after:  tauR=1.13  tauC=1.11 (1.02x)

All six {200,1e3,1e4} x {10,100} cells: pre-fix ratios 1.4-91;
post-fix 0.84-1.62 (independent RNG streams, so MC noise, not bias).
Equivalence exact 18/18 - no equivalence scenario touches rbart_vi.
