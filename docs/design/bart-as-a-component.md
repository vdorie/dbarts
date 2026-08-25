# BART as a component: the internal contract

Two documents cover the same subject at different layers.
`vignettes/dbarts-as-a-component.Rmd` is the USER-FACING one: worked
composition recipes, which channel an outer block writes and which fit it
reads back. This is the INTERNAL one: what the engine guarantees a driver
loop, which mutations are legal between sweeps and under what predicate, and
what engine state a mutation does not carry. A recipe there is correct only
because a guarantee here holds; a guarantee here is worth stating only because
a recipe there depends on it.

The division of labour it serves is `r-c-division.md`'s principle: R addresses
the conditionals, C++ addresses the integrand. The family this document's
multi-forest clauses describe is `multiplier-combiner.md`'s; its "Surfaces"
and "What this family does not do" sections are the authority on what that
family reaches, and are not restated here.

## 1. The per-sweep driver loop is bitwise identical to a batched run

A host that calls `run(0L, 1L)` in its own loop gets, sweep for sweep, the
draws one `run(0L, N)` would have produced. This is exact, not statistical.
It holds because a sampler built with `control@seed` set carries its own
generator per chain, independent of R's stream, so the sweep boundary is not
a place where randomness enters or is re-seeded. It holds for a single-forest
sampler and for a multi-forest one, on the location channel, on sigma, and on
the per-forest fit and amplitude channels.

Two things it does not survive, both by construction: a thinning that differs
between the two routes (the loop advances by one recorded sample per call, so
`n.thin` must match), and a mutation applied between sweeps, which is the
whole point of driving the loop and is what section 2 governs.

The price is a per-call fixed cost, not a per-sweep one. See "Measured".

## 2. Which mutations are legal between sweeps

Single-forest samplers admit every conduit; the predicates below are all
multi-forest, and all live in the bridge's shared header so the R bridge and
the flat C API cannot state different rules
(`src/R_interface_bartcore_common.hpp`).

`refuseMultiForestMutation` (`src/R_interface_bartcore.cpp:2613`) fires on a
bare `numForests >= 2` and covers the whole-object conduits, which would
rebuild or reprice forest 0 alone: `bartcore_setData`, `bartcore_setModel`,
and the flat `dbarts_sampler_setTestOffset`.

`refuseMultiForestResponseMutation` (`:2636`) is the one multi-forest family
that is opt-in rather than refused. It passes a single-forest sampler
unconditionally, then asks the coupling whether it can express a response
swap at all - `Chain::supportsResponseMutation` (`chain.hpp:1068`), which is
the combiner's own answer and nothing else. `AmplitudeForestCombiner` returns
true (`combiner.hpp:1059`); the base `ForestCombiner` and the multinomial
coupling return false, the latter because its response is an n x K count matrix
that no flat conduit can carry, so its refusal names `bartcore_setCounts`
instead. There is no `family_ == gaussian` conjunct: it was removed when
`setResponse` began passing `combinedFits()` rather than forest 0's bare
totals, so a latent family refreshes its latents against the combined location
and the reason for the conjunct went with it. What remains is the scale clause:
a response or offset swap is admitted only at `updateScale == FALSE`, which
pins the response transform the per-forest leaf calibrations are stated
against. The test is `updateScale != FALSE`, so NA refuses too; the R5 methods
default the argument to FALSE. The weight conduit has no scale to pin and skips
the clause.

`refuseUndefinedTestFits` (`src/R_interface_bartcore.cpp:2867`) closes the
test surface, gated on `numForests >= 2 && !testFitsAreDefined` rather than
on the forest count, so a coupling whose test blend IS defined passes
through. It guards `setTestPredictor`, `setTestOffset`,
`setTestPredictorAndOffset` and `predict` on the bridge, and
`dbarts_sampler_setTestPredictors` and `dbarts_sampler_predict` flat.
Because the test predictors can never be installed, the whole surface is
closed rather than partly closed. That is a SAMPLER-level closure: the
combined location at new rows is still available at the fit level, where
`predict.bart` blends the per-forest replay with the stored glue and the
caller's own bases (R/generics.R, `predictBlend`).

R5 raises its own wording ahead of the bridge for a sampler carrying
amplitudes (`refuseAmplitudeMutation`, `R/bartcore.R:36`), on
`setResponse`/`setOffset` at `updateScale = TRUE`, on `setData`, `setModel`
and `setCalibration`. The probe is `!is.null(data@bases)` - a capability test,
deliberately not a forest count, because a K-forest multinomial carries
several forests and no amplitudes.

The predictor surface carries NO multi-forest guard, at any entry, and this is
a decision rather than an omission: the transactional update revalidates every
forest and the variance forest and rolls the whole change back if any leaf of
any tree would empty, and the per-observation session's cell guard caches
every forest pruned to the trees the column can move. A guard named for the
transactional path existed and was retired when those paths learned to loop
the forests. So `setPredictor`, whole-matrix and column-granular and
per-observation, is open at every sampler shape.

## 3. What engine state does not carry, and who reinstalls it

The engine RETAINS the pointer it is handed for the raw conditioning vectors
- response, offset, test offset, case weights, per-forest weights - which the
caller owns and must keep alive; `inst/include/dbarts/dbarts.h:53-66` states
which setters retain and which borrow for the call alone. It holds no
predictor matrix at all: predictors quantize into owned integer codes and the
raw is borrowed for the build or re-quantize call only.

What it DOES own is the derived state a mutation cannot regenerate: the
quantized cut grid, the pinned response transform (gaussian's offset-adjusted
min and max), the tree structure, and the rng position. A response, offset,
weight or sigma swap leaves all four untouched. `updateScale = TRUE` rebuilds
the transform and re-anchors sigma and the variance prior through it, which is
exactly why it is refused wherever a per-forest, variance-forest or grouped
calibration was stated against the old one.

The serialized state carries the derived quantities and NOT the raw ones. It
holds forests, sigma, the transform's endpoints, latents, cut points, DART
state, rng position and the amplitudes; it holds no y, weights, offset,
per-forest weight or active-row mask. So a re-created sampler must be given
them again, and who does that depends on the layer:

- The R5 sampler mirrors each raw vector onto its `data` object as the
  mutation lands (`data@y`, `data@offset`, `data@weights`, `data@x`,
  `data@bases`), so re-creation re-supplies them by construction, and mirrors
  the per-forest weight on an R5 field that `getPointer` and `setState`
  re-apply afterwards (`reapplyForestWeights`, `R/dbarts.R:1805`). There is no
  treatment slot: a Bayesian causal forest's z rides `data@bases` as forest
  2's basis, and moves only through `$setForestBasis`.
- Two holes remain, both known. A per-forest weight is not part of the state,
  so a pipeline that discards the R5 holder and installs a DONOR's state into
  a fresh engine starts with no per-forest weight whatever the donor had, and
  the two stored states compare equal while the fits diverge. This is a
  contract item, decided in `bcf.md:164-174` and pinned for the same-holder
  round trip in `inst/tinytest/test-forest-weights.R`. An active-row mask is
  mirrored nowhere at all and is lost on re-creation.
- A flat C consumer owns every one of them, having no R5 layer to mirror
  through.

## 4. What this document does not commit

Per-forest decomposition of TEST fits: `testFitsAreDefined()` is false under
the BCF coupling, so there is nothing to decompose. Mid-sweep hooks: inside
the tree loop `totalFits` is stale until rebuilt after it
(`src/bartcore/chain.hpp:1456`), so a host reading a fit there would read a
partially updated one; the refusal is on invariant grounds before it is on
threading grounds. Per-forest saved-tree replay. Cross-host bitwise
reproducibility - within-host across any SIMD dispatch only.

## 5. The sweep-boundary hooks that exist

Two ship, and neither is a mid-sweep hook.

`dbarts_sampler_setCallback` (`inst/include/dbarts/dbarts.h:792`) takes
`(userData, sampler, chainIndex, sweepIndex, isBurnIn)` and returns 0 to stop
the run early. It fires at the top of each iteration, unthrottled by thinning,
before sigma enters the sweep. It is refused when `numThreads > 1 &&
numChains > 1`, at registration and again at run: a callback requires chains
to run inline, and inline multi-chain runs them sequentially, so the hook sees
chain c finish before chain c+1 starts.

`bartcore_runWithCallback` (`src/R_interface_bartcore.cpp:4485`) is the
internal single-chain R hook behind `rbart_vi`'s Gibbs loop. It refuses more
than one chain outright, hands the closure one argument - the 0-based sweep
index - and carries no `GetRNGstate`/`PutRNGstate` bracket by design: the
chain's generator never touches R's stream while the closure draws from it, so
R owns `.Random.seed` throughout. The closure is evaluated under `R_tryEval`,
so an error becomes a cooperative stop re-raised in R rather than a longjmp
through C++ frames.

## 6. Graduation debt on the flat surface

A prototype written against R5 and then ported to a `LinkingTo: dbarts`
consumer meets four capabilities the flat surface does not carry: no
per-observation predictor update and no joint session across samplers (the
flat surface's whole predictor channel is `dbarts_sampler_setPredictor` and
`dbarts_sampler_updatePredictor`); no `setCutPoints` and no `setData`; no
`predictVariance`; and no forest-indexed `predict` - `dbarts_sampler_predict`
takes no forest index, and per-forest fits are IN-SAMPLE only, through
`dbarts_sampler_getForestFits`. All four are recorded doors elsewhere; naming
them together here is what makes them on-ramp debt rather than trivia.

## Measured (2026-08-19, arm64 macOS, one host)

Bitwise identity, section 1. A 30-tree two-chain gaussian sampler at n = 200:
10 calls to `run(0L, 1L)` reproduce `run(0L, 10L)` with `identical()` on both
the location and the sigma channel. The same at
`forests = list(forest(), forest(basis = ~ factor(z)))`, on the location,
per-forest fit and amplitude channels. Negative half: a loop at `n.thin = 3`
reproduces a batched run at `n.thin = 3` and differs from one at `n.thin = 1`,
so the pin discriminates. The suite carries the multi-forest half at
`inst/tinytest/test-bcf-mutation-pins.R:145-185` and
`inst/tinytest/test-bcf-reporting.R:84-105`.

Cost of the loop, section 1. Twenty alternating rounds of 500 sweeps, 200
trees, n = 500, p = 5, per-round minimum: batched 0.150-0.152 s, looped
0.165-0.170 s, so the loop costs 1.09-1.13x a batched run of the same length,
equivalently runs at about 0.9x its throughput. Same machine, ratios of
alternated arms; the host was loaded, which the alternation and the
per-round minimum are there to reject.

Cost of an R composition, section 1. K single-forest samplers driven a sweep
at a time against ONE batched engine sampler at matched total trees (200),
same data, same alternation: 1.9-2.1x at K = 2, 2.8x at K = 4, 4.6-5.2x at
K = 8. The tax grows roughly linearly in K, because each sampler carries its
own per-sweep fixed cost - a response install, a weight install, a result
allocation - over the full n. Quoting it without that denominator is what
makes a composition price meaningless.

The refusal matrix, section 2, driven entry by entry through R5 on a
two-forest sampler: `setResponse`/`setOffset` allowed at
`updateScale = FALSE` and refused at TRUE; `setWeights`, `setSigma`,
`setForestWeights`, `setForestBasis`, `getForestAmplitudes` allowed;
`setData`, `setModel`, `setCalibration`, `setTestPredictor`, `predict`
refused; `setPredictor` allowed whole-matrix, column-granular and
per-observation. A single-forest control arm takes `updateScale = TRUE`,
`setData` and `setTestPredictor`.

The rebuild divergence, section 3. A 10-tree treatment forest at n = 200
carrying an excluding per-forest weight, run 20 sweeps, state stored; the
donor continued 5 sweeps against a FRESH engine given that state and run 5
sweeps. The two stored states are `identical()`; the treatment forest's fits
differ by 0.057 against a fit range of 0.035, so the divergence exceeds the
quantity itself.

Semantic continuity, the claim the Gibbs-sampler vignette's multiplier section
states to users. A two-forest
sampler with both amplitudes pinned - glue exactly (a, b0, b1) = (1, 0, 1),
so the model is `y = mu + z tau` - against the same model composed from two
single-forest samplers in R, n = 400, p = 4, sigma fixed, 8 chains of 500
burn-in and 1000 draws each, tree counts and tree prior stated identically on
both routes. With the remedies (y centred at its midrange, k fixed and each
forest's leaf prior set from the engine's own calibration) the posterior mean
CATE agrees to 0.022 rms and 0.067 pointwise on a CATE of 1.96 - 3.4%
pointwise, against a chain-to-chain Monte Carlo standard error of 0.015, so
the rms agreement is inside Monte Carlo error. Without them the difference is
0.053 rms and 0.233 pointwise, 11.9% of the CATE and several times the Monte
Carlo error: observable, which is what makes the documented price honest
rather than reassuring. Both routes recover the truth about equally well
(RMSE against the true CATE 0.226 engine, 0.222 remedied, 0.229 bare), which
is the point - these are two priors that agree closely, not a right one and a
wrong one.

## Status

LANDED, 2026-08-19.
