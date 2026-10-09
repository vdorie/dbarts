# Per-chain responses and offsets in one sampler

Status: FOR DECISION, 2026-10-09 (dec-B419). Nothing is built. Measured on stan4bart's bartcore branch
and the bartcore tip.

## 1. The question

A dbarts sampler holds one response, one offset and one response transform (the multiplier and shift
that map the response less the offset to [-0.5, 0.5]) for all of its chains, as 0.9-34's did. stan4bart
and `rbart_vi` therefore run one sampler per chain. Each chain's offset is the parametric part of that
chain's larger Gibbs sampler. During warm-up each chain re-derives its transform from its own response
less offset (`updateScale = TRUE`), then freezes it, so the chains end on different transforms. The
maintainer asked whether responses and offsets, and so transforms, should be held per chain, so that one
sampler holds every chain of a larger sampler. The maintainer also noted that one shared transform would seem to
need the chains' threads to block and agree.

The note covers what that means for the engine (section 2), three rules for the transform (section 3),
a measurement of what warm-up re-derivation buys (section 4) and a recommendation (section 5).

## 2. What moves from shared to per chain

The engine is closer to per chain than the surface. Each chain owns a
[`ResponseModel`](../../src/bartcore/model.hpp), and the Gaussian one
([`GaussianResponse`](../../src/bartcore/model.hpp)) already holds its own working response (n doubles),
its own pointers to the response and offset, its own transform and its own sigma prior. Everything a
chain reads its transform from goes through its own response model: the fits, the saved-draw restatement
([`Chain::restateSavedDraws`](../../src/bartcore/chain.hpp)), predictions and the per-forest calibration
report ([`Sampler::forestCalibration`](../../src/bartcore/sampler.hpp) already takes a chain). The
sharing is enforced in three places only:

- [`Sampler::setResponse`](../../src/bartcore/sampler.hpp), [`Sampler::setOffset`](../../src/bartcore/sampler.hpp)
  and [`Sampler::setSigma`](../../src/bartcore/sampler.hpp) fan one vector, or one sigma, out to every
  chain.
- The sampler's one transform record (`anchorMin_`, `anchorMax_`) is read from chain 0 at creation and
  at every re-anchor ([`Sampler::recordAnchor`](../../src/bartcore/sampler.hpp)). An install is checked
  against it, and the R object records it as the model attribute `response.range`
  ([`recordAnchor`](../../R/dbarts.R), [`applyAnchor`](../../R/dbarts.R)).
- The bridge owns one copy of each vector, `ownedResponse` and `ownedOffset` on the holder
  ([`dbarts_sampler_t`](../../src/R_interface_bartcore_common.hpp)), which every chain borrows.
  [`dbarts_sampler_setResponse`](../../inst/include/dbarts/dbarts.h) and
  [`dbarts_sampler_setOffset`](../../inst/include/dbarts/dbarts.h) copy into them.

Under the division of [state-not-model.md](state-not-model.md):

| quantity | kind | today | per chain |
|---|---|---|---|
| response, offset | data | one per sampler | one per chain: a larger sampler's state, handed to BART as data |
| working response | scratch | per chain | per chain, unchanged |
| latents | state | per chain | per chain, unchanged ([`dbarts_sampler_getLatents`](../../inst/include/dbarts/dbarts.h) already returns n x chains) |
| sigma written by the host | state where drawn, model where fixed | one value fanned out | one per chain |
| response transform | model, and the units of the state | one per sampler | the open question (section 3) |
| weights, censoring status, predictors, test offset | data | one per sampler | unchanged |

Memory. The response and offset add 2 n (C - 1) doubles for C chains, in the bridge's owned buffers.
Every chain already carries an n-double working response and an n-double residual, so this is small
beside them. Under section 3's option (a) the transform record adds two doubles per chain.

The C API. Three new entries, each taking a chain index: a response setter, an offset setter and a
sigma setter. The existing entries keep fanning out to every chain, so nothing a consumer calls today
changes. It is a minor version bump and a new API hash.

The R surface. `$setResponse`, `$setOffset` and `$setSigma` gain a `chain` argument, where NULL means
every chain as today. A copy and a reload must carry the per-chain vectors. Otherwise a re-created
sampler would revert every chain to the data's one response, so the R object records them beside the
data. The data object itself stays one response.

What carries over unchanged. Chains already run in parallel, each with its own generator. A re-anchor
already rebuilds each chain's working response, restates its saved draws and re-anchors a variance
forest (`Chain::setResponse`, `Chain::setOffset`). The work is the chain index on the setters, the
per-chain buffers and the transform rule.

Threads. [`Sampler::run`](../../src/bartcore/sampler.hpp) starts its worker threads at each call and
joins them before it returns. It also advances the saved-draw cursors for every chain together. A host
that drives a larger sampler calls `run` for one sweep at a time, so every sweep already ends with every
chain stopped. A transform agreed between sweeps blocks nothing that is not already blocked. Threads
would only have to meet if each host thread drove its own chain through a per-chain run entry. That
entry does not exist, and the shared cursors argue against adding it.

## 3. Three rules for the transform

(a) Per chain, as stan4bart's chains freeze it today. Each chain re-derives its own transform and keeps
it. The chains then target slightly different priors: the leaf prior's centre is the transform's shift
and its spread is proportional to the range, so pooled draws are a mixture of posteriors. It breaks
state-not-model's "every chain of a sampler runs under one prior", and the model record becomes per
chain. That means `response.range` per chain, a leaf prior and `k.scale` reported per chain, and a
reload that re-creates each chain at its own transform. Cost over one transform: about 250 lines, and
a sampler whose model is no longer one thing.

(b) One barrier at the end of warm-up. Chains re-derive during warm-up and then agree one transform,
for example the mean over chains of each chain's (min, max). In one sampler this is a re-anchor whose
range reads every chain's response less offset. It is called between sweeps, while every chain is
stopped (section 2), so no thread blocks beyond the per-sweep join. Re-derivation during warm-up still
needs per-chain transforms until the barrier, which is (a)'s machinery, about 30 lines more. Keeping
one shared transform throughout instead, with every re-anchor reading all chains, needs none of it.
Across separate samplers, which is stan4bart's process-per-chain design today, the barrier would need
communication between processes and is not practical.

(c) Fixed at creation from the response less an initial parametric fit. There is one transform, never
re-derived, and no engine work beyond the per-chain setters. The host must hand in the initial fit as
its parametric part will carry it. stan4bart's parametric model has no intercept, because BART carries
the level, so the initial fit must leave its intercept out. Today stan4bart creates its sampler with
the transform taken from the response less the whole `lmer` fit, intercept included. That centres the
BART prior near 0 where the BART component sits near 14 (section 4).

## 4. What warm-up re-derivation buys

Setup. stan4bart's documented example model: Friedman's function, `X4` as a fixed effect, `(1 + X4 | g.1)
+ (1 | g.2)` with 5 and 8 groups, and BART on the remaining nine columns. Continuous and binary
responses, n = 100, 500 and 2000, three seeds each, 4 chains of 1000 warm-up and 1000 draws, 75 trees,
two cores. Four arms per continuous fit, on a scratch build of stan4bart with two switches:

- A, today: re-derive during warm-up, then freeze per chain.
- A2: arm A at another seed, the run-to-run yardstick.
- C0: never re-derive; the transform is creation's, from the response less the whole initial fit.
- C1: never re-derive; the transform is from the response less the initial fit less its mean, which
  stands in for the fit without its intercept. The offset itself is today's.

Fitted quantities are compared to arm A as |difference in posterior mean| / posterior sd, row by row.

Binary fits. The probit response carries no transform (its state records (0, 0)), and re-derivation is
a no-op. Arms A and C0 gave bitwise identical draws in all six fits. Re-derivation buys nothing there.

Continuous fits, the transform's width (which sets the leaf prior's spread). The chains of one fit froze
at widths within 20%, 6% and 2% of each other at n = 100, 500 and 2000. Their mean was within 13%, 5%
and 2% of the creation-time width. Re-derivation hardly moves the spread.

Continuous fits, the transform's centre. The chains of one fit froze at centres up to 0.78, 0.63 and
0.34 widths apart at the three sizes. A chain's centre tracks the level its BART component carries
(correlation 0.93 over the 36 chains). That level trades off against the random intercepts, and it
mixes slowly in every arm: its between-chain R-hat reached 14.4 with re-derivation and 13.5 without. So
the frozen per-chain centres record where each chain's level was when warm-up ended. They do not
estimate a quantity better than creation does. C1's creation centre lay within 0.37 widths of the
chains' mean centre, and within 0.2 in six of nine fits. C0's lay 0.56 to 1.18 widths below it.

Fitted quantities. The median over rows, averaged over seeds, at n = 100 / 500 / 2000:

| quantity | A2 (rerun) | C0 | C1 |
|---|---|---|---|
| expected value | 0.05 / 0.20 / 0.30 | 0.10 / 0.22 / 0.30 | 0.06 / 0.20 / 0.31 |
| BART component | 0.34 / 0.26 / 0.89 | 0.83 / 0.44 / 0.50 | 0.22 / 0.36 / 0.51 |
| random effects | 0.39 / 0.25 / 1.04 | 1.01 / 0.46 / 0.64 | 0.26 / 0.34 / 0.62 |
| residual sd | 0.17 / 0.18 / 0.17 | 0.31 / 0.22 / 0.17 | 0.17 / 0.09 / 0.29 |

The run-to-run yardstick grows with n because the posterior narrows faster than these chains mix. C1 is
indistinguishable from a rerun on every quantity at every size. C0 is too, except at n = 100, where it
moves the split between the BART component and the random effects by 0.8 to 1.0 posterior sd, more
than twice the rerun's difference, and the expected value by a median 0.10 sd against the rerun's 0.05.

Answer. Warm-up re-derivation buys nothing measurable on stan4bart fits, continuous or binary, against a
transform fixed at creation from the response less the initial fit without its intercept. With the
intercept left in, the fixed transform costs a visible shift in the BART / random-effect split on small
data. Only this one model shape was measured; `rbart_vi` has the same shape (random intercepts, BART
carrying the level). The scripts and outputs are in `scratch/pcr/`, untracked.

## 5. Recommendation

Rule (c): one transform per sampler, fixed at creation from the response less the host's initial
parametric fit as the host will carry it. A host that still wants to re-derive gets (b) with no thread
cost: on a sampler whose chains hold different responses or offsets, `updateScale = TRUE` re-derives one
transform from every chain's response less offset (the mean over chains of each chain's range). It is
called between sweeps, while every chain is stopped. Not (a): it gives one sampler several priors, and
the measurement shows nothing for it to buy.

Per-chain responses and offsets with that rule:

- Engine: per-chain setters on [`Sampler`](../../src/bartcore/sampler.hpp), and a pooled re-anchor
  built on [`Chain::moveScale`](../../src/bartcore/chain.hpp) and the re-anchor path's restatement,
  about 140 lines. Range readers: [`GaussianResponse`](../../src/bartcore/model.hpp), the count centre
  of [`NBResponse`](../../src/bartcore/model.hpp), and [`AFTResponse`](../../src/bartcore/model.hpp)
  through its Gaussian part. tests/cpp about 200.
- Bridge and C API: per-chain owned buffers, three dbarts.h entries with their contracts, the version
  and hash, about 200; the C API test about 120.
- R: the `chain` argument, and the copy and reload record of the per-chain vectors, about 120;
  tinytest about 300; help and docs about 180.
- About 1260 lines planned, 2500 to 2900 at the overrun recent engine slices have run.
- It fits after response-scale-rows, which rewrites the same range readers and the creation of the
  transform, so the queue reads k-internal, setState, response-scale-rows, then this. It adds dbarts.h
  entries and R arguments, so if it is built it goes before the release candidate.

What a consumer would need to use it is separate and larger. stan4bart would hold every chain's
parametric sampler in one process, run their steps on its own threads or in turn, call `run` once per
sweep for all chains, and seed its chains through the one sampler. Today's `cores =` would then mean
threads, not processes. That is several hundred lines in stan4bart, not sized here.

Rule (c) is also available to stan4bart today, without per-chain responses: create the sampler from the
initial fit less its intercept, and stop re-deriving during warm-up, about 10 lines. Every chain then
shares one transform, so a restore into one sampler created at that transform is exact without
conversion. That would make dec-B419's per-chain restore, about 120 lines, unnecessary, and it removes
the warm-up re-anchor that the k-internal plan's posterior change for a named sd goes through. It
changes stan4bart's continuous posterior by no more than a rerun does (C1 above). The maintainer's
call: dec-B419 stands until then.
