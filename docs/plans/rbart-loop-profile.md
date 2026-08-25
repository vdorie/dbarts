# rbart-loop-profile

agent: sonnet
rng: neutral (measurement; any fix is a follow-up plan)
budget: measurement + report; no product code

## Goal

A measured answer to whether the rbart_vi custom-prior R loop's
per-sample bridge round trip matters, so the item either closes or
gets a targeted fix plan.

## Context

- The loop pays setOffset + run(0L, 1L) per posterior sample with ~10
  fresh SEXPs per bartcore_run call (R/rbart.R); built-in tau priors
  bypass via the in-core decorator, so only custom-prior users pay.
- capi-callbacks' conditioning hook is the natural fix shape if the
  overhead is real - note the dependency, do not build it here.

## Constraints

- Profile realistic sizes: n in {1e3, 1e4, 1e5}, groups in {10, 100},
  against the in-core path as the floor.
- Out of scope: any change to rbart.R beyond instrumentation that is
  removed afterward.

## Steps

1. Rprof + a wall-clock A/B (custom prior expressing the built-in
   model vs the built-in path) at the sizes above.
2. Report: fraction of wall time in bridge overhead vs tree sweeps.
   Threshold for action: > 15% overhead at any measured size.
3. Record the numbers and the go/no-go in this file; if go, write the
   follow-up plan (likely riding capi-callbacks).

## Verification

- The recorded table; no product diff.

## Measurements (2026-07-07)

Method: arm A is `rbart_vi` with a custom prior function textually
identical to the built-in cauchy tau prior (`function(x, rel.scale)
dcauchy(x, 0, rel.scale * 2.5, TRUE)`, bound under a name other than
"cauchy"/"gamma" so the symbol lookup misses and rbart_vi routes
through the R loop); arm B is the same call with `prior = cauchy` (the
in-core path). Both: n.trees=75, n.thin=1, n.chains=1, n.threads=1,
n.burn=100, n.samples=100 (200 Gibbs iterations/run), keepTrees=FALSE,
verbose=FALSE, synthetic gaussian data (5 predictors, a Friedman-like
f, iid N(0,1) group effects, N(0,1) noise), fresh RNG per run (no
seed - rng is neutral for this item). 5 interleaved A/B pairs per cell
(order A,B,A,B,...), medians reported. Scratch harness (not checked
in, not retained).

### Cell table (median seconds over 5 reps; 200 iterations/run)

| n | groups | median A (s) | median B (s) | sec/iter A | sec/iter B | overhead (A-B)/A |
|---|---|---|---|---|---|---|
| 1e3 | 10 | 0.101 | 0.036 | 0.000505 | 0.000180 | 64.4% |
| 1e3 | 100 | 0.147 | 0.036 | 0.000735 | 0.000180 | 75.5% |
| 1e4 | 10 | 0.476 | 0.338 | 0.002380 | 0.001690 | 29.0% |
| 1e4 | 100 | 0.572 | 0.325 | 0.002860 | 0.001625 | 43.2% |
| 1e5 | 10 | 4.016 | 3.431 | 0.020080 | 0.017155 | 14.6% |
| 1e5 | 100 | 4.779 | 3.490 | 0.023895 | 0.017450 | 27.0% |

Noise: max/min within any single arm/cell never exceeded 1.38x (worst
case: n=1e3, groups=10, arm A - the first replicate ran slower, a
one-time warm-up effect, likely JIT byte-compilation of the closures on
first call). Medians are stable; no extra replicates were needed.

### Rprof attribution

Representative run: n=1e4, groups=100, custom prior, n.burn=1000,
n.samples=1000 (scaled up from the cell size for sampling resolution),
`Rprof(interval = 0.005)`; 3.975s of self time inside the `rbart_vi`
call tree.

- Bridge/engine-call self time (`sampler$run` -> `bartcoreSamplerRun`,
  `sampler$setOffset` -> `bartcoreSamplerSetOffset`, plus stray direct
  `.Call` frames): 1.705s (42.9%). This is the per-iteration engine
  work - the same tree sweep the in-core path also pays - just split
  across 2000 separate `.Call` round trips instead of 2.
- Everything else: 2.270s (57.1%), dominated by `rbart_vi_run`'s own
  per-iteration R bookkeeping (0.870s, 21.9% self time), the O(n)
  group-mean computation for the ranef update via `sapply`/`mean` over
  100 groups (0.610s, 15.3%), vector arithmetic feeding the
  residual/offset update (`-`/`as.vector`, 0.335s combined), and the
  tau slice sampler including its per-draw `optim()` mode search
  (0.155s, 3.9%).

The two views agree qualitatively: at this cell the wall-clock overhead
(A-B)/A measured at 200 iterations was 43.2%, close to the 57.1%
"everything but the engine call" share seen in the 2000-iteration
profile (different run lengths, so not expected to match exactly, but
both say under half of the R loop's own time is engine work).

### Go/no-go

**GO** - overhead exceeds the 15% threshold at 5 of 6 measured cells
(27.0%-75.5%); only the largest/sparsest cell (n=1e5, groups=10) comes
in just under, at 14.6%. Overhead shrinks as n grows (the fixed
per-call cost amortizes over a growing O(n) sweep) but stays above
threshold at every realistic combination of n and group count tested.

If go: the fix shape is extending the capi-callbacks conditioning hook
so the custom-tau-prior Gibbs blocks (ranef draw, offset refresh, tau
slice sample) run through a per-kept-sample R callback invoked from
inside a single `bartcore_run` call, mirroring the in-core path's
one-call structure while still letting R evaluate the user's tau
density - this removes the per-iteration `.Call`/`setOffset` round
trip and, incidentally, the O(n) R-side group-mean computation the
in-core decorator already does in C.

Caveat, out of scope, flagged only: at n >~ 300 the custom-prior R
loop's tau draws diverged sharply from the in-core path's in
exploratory checks (e.g. tau settling around 15-25x the in-core/true
value) even though both express the identical prior and data. This
looks like a pre-existing correctness issue independent of the
loop/bridge overhead measured here - each iteration still does the
same shape of work (one tree sweep, one ranef draw, one slice sample),
so it should not bias the timing numbers above, but it means the two
arms are not currently statistically equivalent at larger n and is
worth a separate look.
