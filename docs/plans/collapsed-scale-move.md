# collapsed-scale-move: the collapsed rescaling move for a drawn k under logistic and probit

Status: PLANNED 2026-10-09. Not yet critiqued; nothing is built. The logistic half is ruled (dec-B404:
built before 1.0-0 if a drawn k is kept, which [binary-k-fixed-or-drawn.md](../design/binary-k-fixed-or-drawn.md)
found). The probit half is the survey's recommendation (dec-B415) and is held for the maintainer
([Open calls](#open-calls), 1); negative binomial is recommended out (2).

agent: opus implementer, one; a blind critique of this plan first; one opus reviewer who runs the mutants
under [Mutants](#mutants).
rng: POSTERIOR-CHANGING by the README's rule (a default changes: every single-forest logistic fit with a
drawn k, and with the probit half every such probit fit, takes a new kernel); the stationary distribution is
unchanged by construction and the exact gate is what shows it. NEUTRAL, bit for bit, everywhere else
([Scope](#scope)).
window: before 1.0-0 (dec-B404); after response-scale-rows in the engine order of dec-A191; engine slices
stay serial.
budget: planned ~1240 lines, ~1860 to ~2480 at the usual overrun; the probit half ~145 of them. Stops:
[Budget and stops](#budget-and-stops).

## Summary for the maintainer

Logistic fits draw k by default, and their k mixes very slowly near separation: tens of thousands of
iterations per independent draw on the worst simulated datasets, and unconverged at 64,000 draws at the
sizes users fit. The probit step built on 2026-10-09 cannot help logistic. This plan builds the move the
survey measured: every leaf multiplied by one factor and k divided by it, the factor drawn from the
binary likelihood with the latent variables averaged out, then the latent variables drawn again. It is
exact, it has no setting, and it runs on about one iteration in four. On the survey's census it took
logistic's 99th percentile of k's autocorrelation time from about 23,000 iterations to 115.

The same code serves probit, beside the built step. There it cut the remaining slow tail about five
times (99th percentile 879 to 160) for 12 to 15 percent more time per iteration, measured on a loaded
laptop and to be confirmed on the x86 box. That half is yours to take or leave; it is a separate commit,
and leaving it out costs nothing to the logistic half.

## Goal

On a single constant-leaf logistic forest with a drawn k (and, with the probit half, a probit one), each
sampling sweep takes, with probability 1/4, a move that multiplies every occupied leaf by a factor r and
divides k by it, r drawn by one slice step from its conditional with the latent variables integrated out,
then redraws the latents at the new fit. Logistic's k reaches the survey's tail (99th percentile of tau(k)
about 115 sweeps on the census) and passes a frozen-tree exact gate and start-from-truth; a logistic-k SBC
arm becomes affordable and joins the matrix. No R, bridge, C API or state-format change.

## Context

- The survey that measured it: [k-tail-samplers.md](../design/k-tail-samplers.md), its "collapsed
  rescaling move" (candidate 4) and its Recommendation. Its prototype is scratch/survey/prototype.patch
  (an older engine; read, not applied), its drivers and results scratch/survey and scratch/survey/out.
- The k study: [binary-k-fixed-or-drawn.md](../design/binary-k-fixed-or-drawn.md) (a drawn k is kept;
  logistic's k unconverged at every length tried; probit's at n = 500 and 75 trees too).
- Logistic and negative binomial mixing: [k-mixing-pg-families.md](../design/k-mixing-pg-families.md).
- The built step: [probit-k-scale-move.md](../design/probit-k-scale-move.md) and its
  plan [probit-k-scale-move.md](probit-k-scale-move.md); in code [`Chain::drawForestRescale`](../../src/bartcore/chain.hpp),
  [`Chain::forestRescaleApplies`](../../src/bartcore/chain.hpp),
  [`sliceFromZero`](../../src/bartcore/chain.hpp), the cache bound
  [`totalFitsScaleBound`](../../src/bartcore/chain.hpp) with
  [`rebuildForestTotalFits`](../../src/bartcore/chain.hpp), the test-only toggle
  [`probitRescaleForest`](../../src/bartcore/chain.hpp) in
  [`SamplerOptions`](../../src/bartcore/chain.hpp), and its tests/cpp block
  [`runForestRescaleTests`](../../tests/cpp/test_ensemble.cpp).
- The sweep: [One sweep](../architecture.md#one-sweep); the cache rule:
  [Forests and combiners](../architecture.md#forests-and-combiners).
- The latents: [`ProbitResponse::refreshLatents`](../../src/bartcore/model.hpp) (truncated normals,
  inactive rows skipped), [`LogisticResponse::refreshLatents`](../../src/bartcore/model.hpp)
  (PG(w_i, psi) with integer count weights w_i, inactive rows skipped); the per-row likelihoods
  [`ProbitResponse::computeLogLikelihood`](../../src/bartcore/model.hpp) and
  [`LogisticResponse::computeLogLikelihood`](../../src/bartcore/model.hpp).
- Rulings: dec-B404 (the logistic move, built if a drawn k is kept), dec-B415 (the survey), dec-B371 and
  dec-B394 (the probit step and its switch), dec-B397 (no switch; an engine toggle for tests only),
  dec-A190 (the probit-k SBC arm at thin 700 "for now"; "We'll revisit everything k-related as we gather
  the evidence."), dec-A191 (the engine order), dec-B412 (negative binomial's prior predictive study).

## The move

One forest, constant leaves, k ~ s chi(nu), each occupied leaf mu ~ N(0, (c / k)^2). With the latents
integrated out the target over (trees, leaves, k) is the prior times the binary likelihood. The group
r > 0 acts as (mu, k) -> (r mu, k / r) on the M occupied leaves and k; eta = k mu / c and the trees are
invariant, and so is every fit up to the factor. With f_i the forest's cached fit
([`Forest`](../../src/bartcore/combiner.hpp)'s `totalFits`, the fit without the offset), o_i the
offset, the sums over active rows, C = k^2 / (2 s^2) and v = log r, the generalized Gibbs step of Liu and
Sabatti (2000) on Haar measure dr / r gives the Jacobian r^(M - 1) against the leaf prior's r^-M and k's
prior's r^-(nu - 1) exp(-C / r^2):

    log p(v) = -nu v - C e^(-2v) + sum_i l_i(e^v f_i + o_i)

    probit:    l_i(eta) = log Phi(s_i eta),                    s_i = 2 y_i - 1
    logistic:  l_i(eta) = w_i (y_i eta - log(1 + e^eta)),      w_i the count weight (1 unweighted)

M drops out with the leaf prior, as in the built step. Empty leaves are fixed at zero and stay there; a
leaf whose members are all inactive is drawn from its prior and scales like any other (the leaf walk of
the built step, so the set scaled is the set the k draw counts). The offset is data and is not scaled.

How the latents leave: under probit the integral of N(z; eta, 1) over s z > 0 is Phi(s eta); under
logistic the Polya-Gamma identity (Polson, Scott and Windle 2013) integrates omega ~ PG(w, 0) out of
exp(kappa eta - omega eta^2 / 2) 2^-w, kappa = w (y - 1/2), to exactly the binomial likelihood above. No
approximation enters; the conditional is evaluated with the likelihood itself.

The draw. The conditional has no closed form and log p is not concave in v in general (a row with
s_i f_i > 0 contributes a locally convex term), so neither an exact inverse nor adaptive rejection
applies. v is moved by one slice step from v = 0 through the built [`sliceFromZero`](../../src/bartcore/chain.hpp)
(Neal 2003: stepping out with the step limit split at random, then shrinkage), at the survey's width 0.5
and the built limit 1000. Each evaluation is one pass over the rows; a step takes about 6 to 10 at the
widths the census and the large-n cells showed (posterior sd of log k 0.1 to 0.9). Then the latents are
redrawn at the new fit by the family's own refresh: truncated normals under probit, PG(w_i, e^v f_i + o_i)
under logistic.

Exactness. Write x for everything but the latents z. If (x, z) is a posterior draw then x is a draw from
the latent-free posterior; the slice step is a reversible kernel for the one-dimensional conditional of
the orbit coordinate given the orbit-invariant ones (trees, eta), so the moved x' is again such a draw;
drawing z' from its exact conditional given x' restores the joint. Every decline below reads only
orbit-invariant quantities (the family, the forest's shape, whether k is drawn, finite and positive,
whether a tree map is stale, whether the prior's scale is infinite), and the coin that schedules the move
is independent of the state, so mixing the move with the identity keeps the posterior. Under C = 0 (an
infinite prior scale) log p grows without bound as v -> -infinity (k -> infinity, the fit -> 0, the
likelihood -> its value at the offset), so the conditional is improper and the move declines; the package
already documents k's posterior as improper there (test-binaryResponse-hyperprior.R). For C > 0 and
nu > 0 the conditional is proper: -C e^(-2v) dominates as v -> -infinity, and -nu v with a likelihood
term at most 0 as v -> +infinity.

What the prototype got wrong or left out, so the build does not copy it: probit declined under a mask;
logistic ignored the count weights and the mask; nothing declined under an infinite scale (the slice
wandered an improper conditional); the schedule was an environment variable with a counter that
setResponse reset; and the signs were read from the latents or the working response, not from y.

## Where it runs

In [`runSweeps`](../../src/bartcore/chain.hpp), after the built step and before the tree loop: level step,
built step, then the collapsed move. Every recorded channel (training and test fits, saved trees, k) is
written after it, as for the built step, and the tree loop reads the redrawn latents' working response and,
under logistic, their precisions as weights.

Eligible when all hold, read before any generator call: the engine toggle is on; one forest, no combiner,
no variance forest; the plain constant leaf; family logistic (and probit, with the probit half); the forest
draws k, k finite and positive; no tree's map stale; the prior's scale finite (C > 0). Otherwise nothing is
drawn. An eligible sweep draws one uniform; below 1/4 it takes the move: one pass for log p(0), the slice
step's passes, the leaf walk, k, the cache, then one latent refresh.

The schedule is a coin, not a counter ([Settled in planning](#settled-in-planning), 1): on average every
fourth sweep, as the survey measured, with no state of its own.

With the built step (probit half): both act on one forest in one sweep. They share the leaf scaling and the
cache rule through one helper extracted from [`Chain::drawForestRescale`](../../src/bartcore/chain.hpp):
occupied leaves times the factor, k divided by it, `totalFits` multiplied in place while the forest's
accumulated factor (`totalFitsScale`) stays within [1/2, 2] and re-derived from the leaves past that. Both
moves multiply into the one accumulated factor, so the gap the cache can carry stays bounded as today.
The built step then scales its latents as it does now; the collapsed move redraws them, so on a sweep that
takes both the built step's latent scaling is overwritten, while its leaves and k stand. Extracting the
helper is a neutral commit of its own: every equivalence scenario bitwise.

## State, chains and threads

- No new chain state. The coin reads the chain's generator, whose state rides every stored state, so
  `setState`, `copy()` and a `saveRDS` reload continue exactly where an uninterrupted chain would, the move
  taken on the same sweeps; a warm start is a new fit with its own generator and carries nothing for it.
  No [`ChainStateData`](../../src/bartcore/combiner.hpp) field, no state element in the bridge, no format
  version ([`stateFormatVersion`](../../src/R_interface_bartcore.cpp)).
- setResponse, setOffset, setWeights, setActiveRows and setData between runs need nothing: the move reads
  y, the offset, the counts and the mask from the response model at the sweep it runs.
- Chains: each chain decides and draws on its own generator; two chains take the move on different sweeps.
  Threads: the move runs inside its chain's sweep with no shared scratch (it needs none: the evaluator is a
  fused loop over the chain's own arrays), so n.threads 1 and 2 draw bitwise the same.
- The likelihood sum is a plain scalar loop in row order in both the shipped and the reference build. It
  must not be vectorized or reassociated: the sum's order is part of the draw, and the equivalence
  baselines are reference-build values the shipped build reproduces bitwise today.
- A cancel cannot land inside the move: the refresh-cancel poll is armed only around the sweep's own
  latent refresh ([`armRefreshCancel`](../../src/bartcore/chain.hpp)), and the move's refresh runs outside
  it, so a cancelled run never leaves the leaves scaled and the latents stale.

## Scope

Reached through, where k is drawn (the default for binary families,
[`drawsLeafKByDefault`](../../R/spec.R)) with a finite scale:
- logistic: `bart`, `dbarts` and `xbart` on a 0/1 response with family logistic, count-weighted fits and
  fits under a mask or an offset included.
- probit, with the probit half: `bart`, `dbarts` and `xbart` on a 0/1 response, `hazard` (one forest), the
  zero part of `hurdle.lognormal`, `rbart_vi`'s binary fit (through its offset), and `bartBT` given
  `k = chi(df, s)` with s finite. The same fits the built step reaches, less those under an infinite scale.

Bitwise unchanged, consuming no draw: every other family (gaussian, Student-t, aft, ordinal, multinomial,
negative binomial); a fixed k (bartBT's default, the SBC's fixed-k logistic arm, logistic-reference.R);
an infinite scale, chi(df, Inf) (equivalence's chik, logistic and wtlogistic; test-binaryResponse-hyperprior.R's
chi(1, Inf) and chi(1000, Inf) fits); BCF and every multi-forest or amplitude-coupled fit, the latent BCF
families included; linear, gp and monotone leaves; grow-from-root's sweeps (an initializer); a sweep whose
tree maps are stale. Without the probit half, every probit fit is bitwise unchanged too.

**The fixed-k case.** The predicate requires a drawn k (the forest's `updateK`), so a fixed k declines
before the coin and draws exactly as today; so does a k that a host installs fixed between runs, from the
next sweep. A k that runs to infinity (an infinite scale with a large df) declines on finiteness.

**Negative binomial: not in this slice** ([Open calls](#open-calls), 2). The move applies in principle
(the likelihood is cheap, the latents are Polya-Gamma), but: it was not measured, and dec-B404's option
was to measure first with a stop if it does not help; its k is not the slow part (99th percentile of
tau(k) 3,253 sweeps where measured, against logistic's 27,505; "For negative binomial alone, k is not the
priority"); its default prior is under its own study (dec-B412), which may change what k it draws; and
its latent redraw is a PG(y_i + r) pass whose cost grows with the counts, and its refresh redraws the
shape too, so the move's cost there is unknown. The code leaves it a door: an evaluator in
[`NBResponse`](../../src/bartcore/model.hpp) and the family in the predicate, plus a census arm and an
exact-gate arm with an exact negative binomial frozen-tree posterior.

## User surface and API

None. No `dbartsControl` argument, no bridge change, nothing in dbarts.h (no ABI event), no state block,
no help page beyond NEWS. The engine keeps a toggle, `SamplerOptions::collapsedRescaleForest = true`,
commented as [`probitRescaleForest`](../../src/bartcore/chain.hpp) is: not a setting, false only in the
engine's own tests. This extends dec-B397's rule to the new move ([Settled in planning](#settled-in-planning), 2).

## Changes

Engine ([chain.hpp](../../src/bartcore/chain.hpp), [model.hpp](../../src/bartcore/model.hpp); names are
proposals, cited by symbol once landed):
- Constants beside the built step's: `collapsedRescaleProbability = 0.25`,
  `collapsedRescaleSliceWidth = 0.5`; the step limit reuses `rescaleSliceStepLimit`.
- `SamplerOptions::collapsedRescaleForest`, as above.
- The helper extracted from `drawForestRescale` (leaves, k, cache), used by both moves; the built step
  keeps its latent scaling and its order of draws.
- `Chain::collapsedRescaleApplies()`: the predicate, no generator. `Chain::maybeDrawCollapsedRescale()`:
  predicate, coin, move; returns whether it moved. `Chain::drawCollapsedRescale(double* factorOut =
  nullptr)`: the move without the coin, factorOut receiving r (1 on a decline), the out-parameter
  tests/cpp reads, as `drawForestRescale`'s `alphaOut`.
- [`ResponseModel`](../../src/bartcore/model.hpp) gains
  `virtual double logLikelihoodAtScaledFits(const double* fits, double factor) const`: the sum over active
  rows, in row order, of the row's log-likelihood at factor * f_i + o_i, by the same per-row expression as
  that family's `computeLogLikelihood`, so the two agree bit for bit. The base returns NaN, unreachable
  past the predicate (as `scaleLatents`' base refuses). Overrides in `ProbitResponse` (the probit half)
  and `LogisticResponse` (count weights as multipliers). The redraw is the family's existing
  `refreshLatents` at the combined fits.
- After the move `runSweeps` re-reads its working-response and weight pointers, as it does after the
  sweep's own refresh (they are stable today; re-reading keeps the invariant local).
- No facade virtual, nothing in dbarts.h, no bridge, no state-format change; `--preclean` on every engine
  commit all the same (model.hpp's virtual table changes).

Scripts:
- [probit-k-scale-exact.R](../../benchmarks/R/probit-k-scale-exact.R) gains a link argument and logistic
  arms ([The exact gate](#the-exact-gate)); the name stays, so the workflow list and its cites stand.
- [probit-k-mixing.R](../../benchmarks/R/probit-k-mixing.R): a fresh sampler per dataset, seeded by the
  dataset's index, in [`runCensus`](../../benchmarks/R/probit-k-mixing.R) (today one sampler serves every
  dataset and the previous chain advances its generator, so per-dataset results cannot be compared across
  builds that draw differently; the survey's census.R made this fix and reproduced the published census,
  median 73 against 71), and a link argument for `census` and `truth`. Its section in
  benchmarks/README.md follows.
- [sbc.R](../../benchmarks/R/sbc.R): [`sbcConfigProbitK`](../../benchmarks/R/sbc.R) generalized to a
  family argument, a `logistic-k` arm, [`sbcMatrixConfigs`](../../benchmarks/R/sbc.R) and the matrix
  count M (84 to 91); [sbc.yaml](../../.github/workflows/sbc.yaml) gains the job; benchmarks/README.md
  its line.
- [equivalence.R](../../benchmarks/R/equivalence.R) gains one scenario, `logistick`: the sampler API,
  family logistic, the default chi(1.5, 2), count weights 1 to 3 and an offset, so the move's w and o
  terms have a bitwise anchor (every logistic scenario today uses chi(1.5, Inf), which the move declines).

Docs at landing: docs/design/collapsed-scale-move.md (the derivation and exactness above, the schedule,
the gates' landed numbers, the x86 cost) and its INDEX row; k-tail-samplers.md's and
k-mixing-pg-families.md's Status lines point to it; [One sweep](../architecture.md#one-sweep) step 2 gains
the move; inst/NEWS.Rd's tree-moves item gains a sentence (users' view, against 0.9-34: binary fits that
draw k take a further step on about one iteration in four that rescales the leaf values and k together,
drawn with the latent variables averaged out, so k and the fit's overall scale mix far faster near
separation; the posterior is unchanged), parse-gated per the README's checklist; TODO k-mixing-pg-families
closed; the MANIFEST row; a ledger entry for the calls made; this plan's Status and Landing note.

## Tests

tests/cpp, in [test_ensemble.cpp](../../tests/cpp/test_ensemble.cpp) beside
[`runForestRescaleTests`](../../tests/cpp/test_ensemble.cpp), through
[`TestPeer`](../../tests/cpp/test_peer.hpp). Fixtures: a burned-in drawn-k sampler per link with an offset,
a mask, two or more trees and a stranded nonzero empty leaf; the logistic one also with counts 1 to 3.
Reference constants (nu, C, the per-row terms) come from the raw state: y, the offset, the mask, the
counts, the leaves re-summed in tree order and the k prior, never the move's own code.
1. Evaluator: `logLikelihoodAtScaledFits(f, a)` equals the row-order sum of `computeLogLikelihood` at a f
   over active rows, bitwise, at several a, per link.
2. Conditional, in invariance form (a slice step from a fixed point is not a draw from the target): per
   replication draw v0 from the frozen state's target by inverse-CDF quadrature, move the frozen state to
   the orbit point e^v0 (leaves and k only; the move does not read the latents), take one
   `drawCollapsedRescale`, read v1 from the leaves and from k, which must agree. KS of the v1 against the
   target, D < 1.95 / sqrt(N) (p = 1e-3), N = 2e5, or 5e4 if a link's run passes 10 s. Per link.
3. Redraw: from a fixed state, many moves; under probit each active row's latent sits on y's side and its
   PIT under N(e^v f_i + o_i, 1) truncated to that side is uniform (KS on a few rows, pooled); under
   logistic each sampled omega_i standardized by PG(w_i, psi)'s mean w_i tanh(psi / 2) / (2 psi) and
   variance at the post-move psi averages to zero (|z| < 4.5). Inactive rows' latents unchanged.
4. Mapping, one move with r from the out-parameter: k' r = k; mu' = r mu on every occupied leaf of every
   tree, the empty leaf untouched; `totalFits` equal to the tree-order re-sum bitwise when the bound forced
   a rebuild and within 1e-14 relative otherwise, `totalFitsScale` multiplied by r; under probit the latents
   are not r z (redrawn).
5. Inertness: toggle off, fixed k, an infinite scale, k infinite, a stale map, gaussian, negative binomial,
   ordinal, a two-forest fit and a linear leaf consume no generator draw (the next uniform unchanged); an
   eligible sweep whose coin says no consumes exactly one.
6. The coin: `maybeDrawCollapsedRescale` called 4000 times on an eligible frozen state moves on a share
   within 4.5 binomial sd of 1/4.
7. Every chain: on a two-chain sampler each chain's k differs from a toggle-off twin after 50 sweeps.
8. With the built step (probit half): one sweep taking both leaves `totalFitsScale` the product of the two
   factors or 1 after a rebuild, and the debug cache check
   ([`assertForestCachesMatchLeaves`](../../src/bartcore/chain.hpp)) silent.

tinytest, a new file inst/tinytest/test-collapsed-rescale.R, both links unless said:
- Per chain, on any number of threads: a two-chain drawn-k fit draws identically at n.threads 1 and 2.
- Ahead of every recorded channel: a keepTrees fit's saved trees replayed reproduce each kept draw's
  training and test fits within 1e-10.
- Finite under an offset, a mask installed mid-run and (logistic) count weights, at a finite scale; a
  chi(1.5, Inf) logistic fit runs.
- State: a sampler restored by setState from a store and a reload of the same store continue identically.
- [test-composition-sequences.R](../../inst/tinytest/test-composition-sequences.R) gains `probitK` and
  `logisticK` arms in [`continuationArms`](../../inst/tinytest/test-composition-sequences.R): restore after 3
  sweeps (not a multiple of 4) and the move sequence matches the uninterrupted run exactly. This is the
  test a schedule counter left out of the state would fail.
- [test-probit-rescale-forest.R](../../inst/tinytest/test-probit-rescale-forest.R)'s header says its probit
  fits take both moves (probit half).

## Gates

On the slice tip against its own library, independently of the implementer
([RNG classes and their gates](README.md#rng-classes-and-their-gates), posterior-changing).

### The exact gate

[probit-k-scale-exact.R](../../benchmarks/R/probit-k-scale-exact.R), extended. One frozen tree, k ~
chi(1.5, 2), n = 150, the exact posterior of k and the leaf means by one quadrature a leaf
([`exactPosterior`](../../benchmarks/R/probit-k-scale-exact.R)), now per link: probit's l_i as today,
logistic's with plogis and the counts as powers, the leaf scale per link (3 for one probit tree, pi sqrt(3)
for logistic) read from the sampler rather than written in. Arms: probit pure, mixed, offset, mask (now with
both moves, probit half); logistic pure, mixed, offset, mask and counts (weights 1 to 3). Same statistics
and |z| bound 4.5. Mixing bounds:
- probit: full 0.01 and quick [`QUICK_MIXING_BOUND`](../../benchmarks/R/probit-k-scale-exact.R) 0.014
  unchanged. The built step alone sits at 0.0064 (arm64) and 0.0099 (x86-64), so no bound between it and
  the two moves together is safe across hosts; a probit half that silently declines is caught by tests/cpp
  6 and 7 and the census, not here.
- logistic: quick and full bounds set at landing between the move's value and the value without it,
  measured on a scratch build whose toggle default is flipped (never committed), and recorded in the
  header, so a logistic move that silently declines fails per push.
Run as `Rscript benchmarks/R/probit-k-scale-exact.R quick` (both links) in
[exact-gates.yaml](../../.github/workflows/exact-gates.yaml)'s existing list; quick mode grows from about
30 s to an estimated 2 to 3 min (logistic's sweeps cost more), inside the job's 40 min. The survey's
prototype passed quick with the probit move every sweep at worst |z| 2.9, the pure arm's error 0.0031;
it never ran a logistic arm.

### Start-from-truth

The mixing script's `truth` action per link: an exact posterior draw run on, log k at lags 300 to 3000;
|z| <= 3 at R = 2000. Survey: probit +0.57, logistic +1.23.

### The mixing gate (pre-registered)

The census per link, 100 datasets of the probit-k arm's generator (150 rows, 3 predictors, 50 trees, k from
chi(1.5, 2); logistic the same with a logistic response, as the survey), each on its own sampler, 30,000
sweeps of burn-in and 200,000 recorded at thin 10. tau in sweeps; the bar is 1.5 times the survey's every-
fourth-sweep figures, allowing for the coin (an estimated 1.0 to 1.35 times a fixed schedule where the move
dominates, less in the tail where other directions do) and the noise of the slowest one or two datasets:

| link | survey: tau(k) median / 90th / 99th / worst | tau(fit) worst | bar: tau(k) median, 99th | bar: tau(fit) |
|---|---|---|---|---|
| probit, both moves | 35 / 116 / 160 / 166 | 114 (99th 112) | 55, 240 | 99th 170 |
| logistic | 30 / 61 / 115 / 127 | 59 | 45, 175 | worst 90 |

Passing the bar is the landing's mixing evidence; missing it is stop 2. The run is in sweeps, so a laptop
under load serves. The census's per-dataset lag-thin autocorrelations also set the logistic-k arm's thin
and record what the probit-k arm could run at.

### The SBC arms

- logistic-k, new: R = 600, 99 draws, thin the smallest multiple of 50 at which every census dataset keeps
  50 or more effective draws of 99 by both of the probit-k arm's estimates (the lag-thin autocorrelation and
  k's tau; expected 100 to 200 from the survey's worst 127), burn at least nine times the census's slowest
  tau(k) (expected about 2000). Timeout 3 times its measured laptop time. It flags k on its ecdf band or a
  chi-square p below the band's alpha, as probit-k does. Without the move no such arm could run: logistic's
  k needed tens of thousands of sweeps a draw.
- probit-k: settings unchanged in this slice (R = 600, 99 draws, thin 700, 10,000 burn; dec-A190). With the
  probit half its draws move and it must pass as an ordinary SBC on the confirmation run; the census's
  reading of what thin would keep 50 effective draws (expected near 200) is recorded for the k revisit
  dec-A190 names, not applied here.
- Every other arm draws bitwise (fixed k or another family). This slice edits sbc.yaml, so its push starts
  the full matrix: that run is the confirmation (stop 7).

### Equivalence

On the reference build against the MANIFEST's current baseline at the slice's base (equivalence-734441f1
today; the slices ahead of this one re-record it first):
- With the probit half: every scenario but four identical draws (same RNG stream). The movers are exactly
  maskprobit, hazard, hurdle and bart2probit, the single-forest probit fits with a drawn k at a finite
  scale. chik takes the built step only (chi(1.5, Inf)) and stays bitwise, which also checks that the
  helper extraction is neutral; logistic and wtlogistic (chi(1.5, Inf)) and nbinom stay bitwise.
- Without it: every existing scenario identical.
- `logistick` is new: recorded fresh, reported skipped by the old baseline.
- bcf-equivalence 15 of 15 and multinomial-equivalence 11 of 11 bitwise (two forests decline).
- Re-record as equivalence-<sha>: the movers and `logistick` recorded fresh and merged into a copy of the
  base's baseline in scenario order, a `--bitwise --strict-coverage` compare against the merged file from a
  second `--preclean` reference install, the MANIFEST row. ORACLE (P17): the exact gate's quick and full
  passes on both links and the start-from-truth statistic. The z sizes the move and is not the oracle; the
  built step's re-record met |z| above 4 from its old baseline's short, unconverged chains in k, and the
  same reading may recur.

### Snapshots

The four reference-build snapshot files pass unchanged: test-reproducibility-binaryResponse.R fits through
bartBT at a fixed k = 4.5, the xbart file over a fixed k grid, the other two gaussian. A snapshot that
moves is a stop.

### Cost on the x86 box

The survey's costs came from a laptop at load 20 to 40, and logistic's was never measured. On the idle x86
box (`scratch/x86-box start`, stopped after), the built step's method: per-sweep time of the slice's build
over the same build with the toggle default off (a second private library from a scratch patch, never
committed), both restored from one stored state before every timed segment, ten interleaved repetitions in
alternating order, one chain, one thread; median ratio and range. Cells: 150 rows and 50 trees, then n in
{500, 2000, 5000, 20000} by trees in {75, 200}; probit (both moves against the built step alone) and
logistic. The ratio includes the coin, so it is the average cost; the cost of one move is (ratio - 1) / 0.25.
Laptop reference: every fourth sweep 1.15 at 150 rows and 50 trees, 1.12 at 5000 and 200 (probit).

If a cell passes its stop, the first remedy is a faster log Phi for the move's evaluator (erfc with the
asymptotic tail, as the prototype's), accepted only if it matches the family's `computeLogLikelihood` within
1e-14 relative over the rows' range; tests/cpp 1 then compares within that tolerance.

### The rest

tests/cpp plain and under address and undefined-behavior sanitizers; the R-loaded ASAN path over the new
and touched test files; the full tinytest suite (a seeded expectation on a drawn-k binary fit that moves is
re-pinned with its reason, a statistical one must pass as is); every exact-gates.yaml gate in quick mode and
this one in full (hazard-reduction and hurdle-reduction move both sides together and must stay bitwise
equal; logistic-reference and bcf-latent-exact fix k or use two forests and draw bitwise); `R CMD check
--as-cran` from a clean tarball; lintr, air, rc-codoc, win-drift, doc-freshness; bench-sampler compare
(maintainer): with the probit half
["run-binary-n1000-p10-t75"](../../benchmarks/R/bench-sampler.R) takes the move and is expected past 1.05
by its own cost, so its MANIFEST row gains that note or the maintainer re-records the speed baseline.
Consumers: stan4bart's, bartCause's, treatSens's and bairrtt's suites against the build; a seeded
expectation that moves on a drawn-k binary fit is re-recorded on that package's lockstep branch, the
orchestrator pushing.

## Mutants

Each must fail a named test or gate:
1. The k prior's power without the Jacobian (-(nu - 1) v), or with the leaf prior's M left in: tests/cpp 2,
   the exact gate.
2. The C term dropped or with the wrong sign of v: tests/cpp 2.
3. The offset scaled with the fit (e^v (f + o)): tests/cpp 1 and 2, the gate's offset arms.
4. Inactive rows in the sum: tests/cpp 1 and 2, the mask arms.
5. Logistic's count weight dropped: tests/cpp 1 and 2, the counts arm.
6. Probit's sign flipped for failures: tests/cpp 1 and 2.
7. The latents not redrawn, redrawn at the pre-move fit, or (probit) scaled by r instead: tests/cpp 3 and
   4, then the gate and start-from-truth.
8. k not divided by r: tests/cpp 4 (the gate is blind to this, as it was for the built step).
9. Only tree 0's leaves scaled, or the empty leaf scaled: tests/cpp 4.
10. `totalFits` multiplied without accumulating `totalFitsScale`, or never rebuilt: tests/cpp 4 and 8, the
    debug cache check.
11. The coin drawn before an ineligible forest declines, or the C = 0 guard removed: tests/cpp 5, chik,
    logistic and wtlogistic bitwise.
12. The coin at 1 or at 0: tests/cpp 6; at 0 also the logistic mixing bound and the census.
13. The move placed after the tree loop: the keepTrees replay.
14. Chain 0 only, or scratch shared across chains: tests/cpp 7, the thread-count tinytest.
15. The helper extraction reordering the built step's draws: chik's bitwise compare.
16. A schedule counter in place of the coin, not stored: the continuation arms.

Read, not run: the likelihood sum stays a scalar loop in row order, and the move's refresh sits outside the
armed cancel window.

## Budget and stops

Planned lines by part (the probit half's share in parentheses), with the usual overrun of 1.5 to 2 times:

| part | planned | at 1.5x to 2x |
|---|---|---|
| engine (chain.hpp, model.hpp) | 170 (25) | 255 to 340 |
| tests/cpp (test_ensemble.cpp, test_peer.hpp) | 380 (75) | 570 to 760 |
| tinytest | 120 (15) | 180 to 240 |
| exact gate script | 110 (0) | 165 to 220 |
| mixing script | 70 (0) | 105 to 140 |
| SBC harness, sbc.yaml, README line | 70 (0) | 105 to 140 |
| equivalence scenario | 30 (0) | 45 to 60 |
| design note, architecture, NEWS, README sections, TODO, indexes | 250 (30) | 375 to 500 |
| records (MANIFEST row, ledger, landing note) | 40 (0) | 60 to 80 |
| total | ~1240 (~145) | ~1860 to ~2480 |

The re-recorded baseline is data. Compute: two censuses (about an hour each at one core), two truth runs
(about 5 and 7 min), the gate in full (a few minutes a link), a local logistic-k SBC smoke at R = 100, and
the x86 cost cells (about two hours of the box).

Stop and report, without working around, when (1 to 6 pre-registered here):
1. The exact gate fails with the move (any |z| > 4.5, or a pure arm's error above its full bound) after one
   fix.
2. The census misses its bar on either link.
3. Start-from-truth passes |z| 3 at R = 2000 on either link.
4. On the x86 box, after the log Phi remedy, any probit cell's ratio passes 1.20 or any logistic cell's
   passes 1.25.
5. An equivalence scenario other than the expected movers moves, or a snapshot moves.
6. The diff passes 2480 lines, or one part passes twice its planned line.
7. After landing: the confirmation matrix flags k on probit-k or logistic-k (an end bin above 45 or a
   chi-square p below 0.01); it is investigated before anything else k-related lands.

## Order of work

1. Branch off origin/bartcore at a green engine tip carrying response-scale-rows.
2. Extract the leaf, k and cache helper from `drawForestRescale`; tests/cpp; the reference-build partition
   all bitwise; commit (neutral).
3. The logistic half: constants, toggle, predicate, coin, move, `LogisticResponse`'s evaluator; tests/cpp
   1 to 7 for logistic, plain and sanitized; commit.
4. The probit half, as its own commit: `ProbitResponse`'s evaluator, the family in the predicate, tests/cpp
   for probit and 8. If the maintainer rules it out, this commit is dropped before review.
5. tinytest; `--preclean`; full suite; commit.
6. The exact gate (both links, the logistic bounds from a flipped scratch build), quick and full; commit.
7. The mixing script's fix and link argument; census and truth on both links; commit (outputs stay in
   scratch).
8. The logistic-k arm, sbc.yaml, the README line; a local R = 100 smoke; commit.
9. `logistick`; the reference build with `--preclean`: the partition, the re-record and its MANIFEST row;
   the four snapshot files; commit.
10. The design note, architecture, NEWS, README sections, TODO; commit.
11. The x86 cost cells when the box is free; their table into the design note.
12. Consumers' suites.

## Queue slot and shared files

dec-A191's order: default-rule-per-gap (landing), then the k-internal slice (`k-internal-parameter.md`, on
branch wt/k-internal-plan), then setstate-force-update, then
[response-scale-rows.md](response-scale-rows.md), then this. It needs the k-internal slice in first, as
dec-A191 says: the move is written once against k's final meaning. That slice leaves the probit step
"Unchanged" and keeps a binary family's k against the fixed 3 or pi sqrt(3); an sd-spelled binary prior
becomes chi(df, k.scale / s), which this move reads like any other chi prior.

Files this slice shares, to rebase onto each:
- default-rule-per-gap: model.hpp (its tree-prior hunks, away from the response classes); the current
  equivalence baseline and MANIFEST (it re-records; this slice merges onto whatever is current); its
  exact-gates.yaml list edits (this slice changes no list entry); inst/NEWS.Rd; TODO; docs/plans/INDEX.md.
- k-internal: chain.hpp (its install, translation, setModel and `scaleDrawnK` hunks, away from
  `runSweeps`' step block and `drawForestRescale`); tests/cpp (its test_state.cpp and test_sampler.cpp;
  test_peer.hpp if it adds accessors); docs/architecture.md; inst/NEWS.Rd; TODO; docs/plans/INDEX.md.
- setstate-force-update: chain.hpp (`Chain::setState` and the state validity checks); equivalence.R (its
  scenarios beside `logistick`) and the baseline lineage; inst/NEWS.Rd; docs/plans/INDEX.md. This slice's
  continuation arms in test-composition-sequences.R run through its setState.
- response-scale-rows: chain.hpp and model.hpp (creation-time scale, binary families bitwise there);
  exact-gates.yaml; MANIFEST; inst/NEWS.Rd; docs/plans/INDEX.md.

Conflicts are expected in the records files only. After this slice: the binary hyperprior rerun of
dec-B390 (binary-k-default-revisit) runs on chains carrying both moves.

## Settled in planning

1. **A coin, not a counter.** The survey took the move on every fourth sweep by a counter. A counter is
   chain state: to keep a restored chain on the uninterrupted one's draws it must ride the stored state (a
   [`ChainStateData`](../../src/bartcore/combiner.hpp) field beside
   [`dartNumUpdatesSkipped`](../../src/bartcore/combiner.hpp), its element written and read at the
   bridge's two state parse sites, the sampler's state copy, a rule for what resets it), about 80 lines in
   files the k-internal and setState slices rewrite. A coin of probability 1/4 has the same average cost,
   rides the generator's state, and keeps the posterior trivially. Its cost is mixing where the move
   dominates: a fixed schedule decorrelates a little faster, by an estimated 1.0 to 1.35 times in tau. The
   census bar allows for that; if the coin misses it, the counter is the fallback, inside the budget's
   overrun.
2. **No switch.** dec-B397's rule for the built step ("the move runs wherever it applies ... with no
   R-level control; the engine keeps an internal toggle for its own tests") is extended to this move. Its
   premise was a nearly free step (0 to 3 percent of a sweep); this one costs 12 to 15 percent on the
   laptop, so the x86 table goes to the maintainer with the landing, and stop 4 reopens the question.
3. **Placement ahead of the trees, after the built step.** Placing the move just before the sweep's own
   latent refresh would save its separate redraw (under logistic a PG pass, possibly the larger part of
   the move's cost there), but every recorded channel would then be written before it and have to be rescaled in
   place mid-sweep (the test fits' accumulator and the saved trees' leaves), breaking the built step's
   "every recorded channel is written after it". Not taken; a door if the x86 cost presses.
4. **One slice step**, through the built `sliceFromZero` at the survey's width 0.5. Neither exact nor
   adaptive rejection is available (above). A two-evaluation Metropolis step would cost fewer passes but
   was not measured on the census and needs a step size that cannot follow the posterior's width (0.1 to
   0.9 across n) without depending on the state.
5. **A logistic equivalence scenario** (`logistick`): without it no baseline anchors the logistic move.
6. **The evaluator in the response model**, so the chain stays family-agnostic and the move's likelihood
   is the family's own reported one.
7. **The probit-k arm's thin is not changed here**; the census records what it could be (dec-A190).

## Open calls

1. **The probit half.** Recommended: include it, the survey's recommendation, unless the x86 cost at
   n >= 2000 passes 1.15 where the move gains least.
   - Including it: the census's remaining probit tail falls about five times (tau(k) 99th 879 to 160,
     worst 959 to 166; the fitted values' 99th 1003 to 112); on the ten slowest datasets the fitted values'
     worst falls from 1050 to 45. At larger n it gained up to two times or nothing (one dataset a cell). It
     costs 12 to 15 percent of a sweep on every default probit binary fit (laptop; x86 to confirm), about
     145 planned lines (220 to 290 at the overrun), four re-recorded equivalence scenarios, bench-sampler's
     binary row past its 1.05, and probit expectations re-recorded in consumers. It would let the probit-k
     SBC arm drop from thin 700 toward about 200, a third of its runner time, by a later call. Whether it
     also converges probit's k at n = 500 and 75 trees, where the k study found split-Rhat 1.12 to 1.21 at
     64,000 draws, is not measured.
   - Leaving it out: probit keeps the built step alone, at today's speed and draws; the slow tenth of small,
     near-separated datasets still needs several hundred sweeps per independent draw of k and of the fitted
     values. The logistic half is unaffected: the probit commit carries only the evaluator, the family in the
     predicate and its tests.
2. **Negative binomial.** dec-B404's option described the move as "covering negative binomial too".
   Recommended: not in this slice, for the reasons under [Scope](#scope); revisit after dec-B412's prior
   study, measured first as dec-B404's option said. Including it now: about 250 more planned lines (the
   evaluator, an exact negative binomial frozen-tree posterior for the gate, a census arm, the nbinom
   equivalence scenario re-recorded) and an unmeasured cost whose latent redraw grows with the counts.
