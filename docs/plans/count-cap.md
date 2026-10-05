# count-cap: cap logistic count weights, a fixed negative-binomial shape and multinomial trials at one million; interrupt long sweeps

Status: IMPLEMENTED 2026-10-04 on wt/count-cap, under dec-B202 in [decisions.md](../decisions.md) as the
maintainer extended it to multinomial trials; not yet landed. Logistic and nbinom at cc9a6a8a, multinomial
at 6334bffb, the first review's corrections at 58dc34e2 (the poll armed only around a sweep's own draws, a
linear leaf's cached U'WU dropped on a stop, one refusal phrasing). A sampled shape cannot reach the cap, its grid stopping at 50. An interrupted refresh
that drew the shape is put back whole, since its rows left at the previous omega would pair a draw at the old
shape with the new one. An interrupted multinomial glue draw leaves no latent in the chain's state: omega and
the margins are per-sweep scratch, written before each read, so what stands is the tree updates of the
categories before the stop.

agent: opus implementer, one; opus reviewer.
rng: NEUTRAL. Accepted inputs below the cap take the same path and the interrupt poll draws nothing; only
inputs above the cap, which no run could finish, are refused.
window: pre-release. All three families are new in 1.0-0.
budget: ~450 lines (R ~40, C++ ~120, tinytest ~150, tests/cpp ~80, manual ~20, records ~40). Plans have run
1.5-2x low.

## Goal

A logistic count weight, a fixed negative-binomial shape or a multinomial count row totalling more than one
million trials is refused by name on every surface that takes it, as a negative-binomial response above one
million already is. A user interrupt during
a sweep's Polya-Gamma draws stops the run within about a tenth of a second and leaves the sampler valid.

## Context

- The exact Polya-Gamma draw costs one draw per trial: about 160 ns per trial on the reference machine, so 100
  rows at a million trials take 16 seconds a sweep and a row of a billion 160 seconds, during which no
  interrupt is honoured (dec-B202 records the timings and the ruling).
- The response's cap is [`maximumCount`](../../src/bartcore/model.hpp), applied at creation in R and in
  the bridge at every mutation. The weight checks are [`refuseNonCountWeights`](../../R/spec.R) and the
  bridge's [`enforceBinaryWeightPolicy`](../../src/R_interface_bartcore.cpp); the shape enters through the
  `nbinom()` constructor and the shape resolution in R, and through a state install.
- A multinomial row of n_i trials costs n_i draws per category per sweep in
  [`MultinomialForestCombiner::drawForestGlue`](../../src/bartcore/combiner.hpp), interleaved with the
  category forests' tree updates. Count rows enter at creation (bart, dbarts) and through
  [`bartcore_setCounts`](../../src/R_interface_bartcore.cpp); xbart, `$setData` and state installs take
  none, and the posterior predictive draw is a binomial per category, costing nothing per trial.
- The run loop polls for interrupts between sweeps, throttled to about 100 ms, on the main thread
  ([`Sampler`](../../src/bartcore/sampler.hpp)); the monotone order count already aborts inside a move on an
  interrupt and leaves the chain valid, which is the precedent.

## Constraints

- One cap, one constant, shared by the response, the weights, the shape and the multinomial trials; messages name the argument and
  the cap.
- Out of scope: a sampled shape whose prior reaches past the cap (report whether it can, do not change the
  prior); posterior predictive trial counts, whose binomial draw costs nothing per trial.
- The engine stays R-agnostic: the interrupt reaches the draw loop through the run's existing poll, not an R
  call; threads other than the main one do not poll R. The poll is armed only around the sweep's own draws:
  the same draws reached from a host hook mid-run have nothing to catch the stop.
- Out of scope: an interrupt inside the draws made outside a run - `$setWeights`, `$setResponse`, `$setData`,
  `$setState`'s weight re-derivation and `dbartsDrawLatents`. The ruling names the in-sweep interrupt.
- After an interrupted sweep every latent is either this sweep's draw or the previous one, both valid states
  of the chain; say so in the code and test that the sampler runs and restores afterwards.
- No NEWS: the families are new in 1.0-0; the help pages for the logistic weights, `nbinom()` and the
  multinomial counts state the cap.

## Steps

1. Cap logistic count weights and a fixed shape at the shared constant on every surface (creation in R,
   `$setWeights`, `$setData`, the bridge, a state install carrying a shape). Tinytest per surface at the cap
   (accepted) and just above it (refused by name). tests/cpp for the engine-side install check.
2. Poll for an interrupt inside the logistic and negative-binomial latent refresh at the run loop's
   throttle; abort the sweep on an interrupt as the monotone count does. tests/cpp: an interrupt injected
   mid-refresh leaves every latent finite and the chain runs on. Tinytest where an interrupt can be
   simulated in-process; otherwise tests/cpp only.
3. Help pages: the cap and the reason (time per sweep grows with the total trials).
4. Cap a multinomial count row's trial total at the shared constant at creation and `$setCounts`, and poll
   inside the multinomial glue draws, stating what an interrupted one leaves. Tinytest at the cap and above
   it, and an in-process interrupt; tests/cpp: a stop in a category leaves the earlier categories updated and
   the rest untouched, and the chain runs on and restores.

## Verification

- `R CMD INSTALL --preclean` into the slice's own library; full tinytest; `tests/cpp`; the four seeded
  snapshot files on the reference build; equivalence compare in statistical mode, all scenarios identical.
- Timing check: a 100-row logistic fit at a million trials per row, and a multinomial fit at a million
  trials per row, still run, and an interrupt sent during the first sweep returns within a second.
- `lintr::lint_package()`, `air format --check .`, the three tools/ checks, and `R CMD check --as-cran`.
