# collapsed-scale-move: the collapsed rescaling move for a drawn k under logistic and probit

Status: FOR DECISION 2026-10-09; revised from its blind critique (scratch/cmp/critique.md) and from the
maintainer's notes on the revision; rechecked 2026-10-11 against d5b637ee, bartcore with slice B of
response-scale-rows ([Recheck against bartcore](#recheck-against-bartcore)). Nothing is built. The
logistic half is RULED: dec-B404 builds the logistic move before 1.0-0 if a drawn k is kept, and dec-B420
keeps a drawn k under chi(1.5, 2) for both links
([binary-k-fixed-or-drawn.md](../design/binary-k-fixed-or-drawn.md)). Three calls stay open
([Open calls](#open-calls)): the probit half, which waits on the large-n reading, a switch, and negative
binomial. Steps 1 to 6 of [Order of work](#order-of-work) need none of the three and can start on d5b637ee.
The recheck changed no rule and opened no call.

agent: opus implementer, one; one opus reviewer who runs the mutants under [Mutants](#mutants). The blind
critique ran and is answered here.
rng: POSTERIOR-CHANGING by the README's rule (a default changes: every single-forest logistic fit with a
drawn k, and with the probit half every such probit fit, takes a new kernel); the stationary distribution is
unchanged by construction and the exact gate is what shows it. NEUTRAL, bit for bit, everywhere else
([Scope](#scope)).
window: before 1.0-0 (dec-B404). Everything dec-A191 queued ahead of it has landed, slice B of
response-scale-rows last; engine slices stay serial. The probit and switch calls wait on the early
measurements (step 7 of [Order of work](#order-of-work)), the probit call above all on their large cells.
budget: planned ~1445 lines, ~2170 to ~2890 at the overrun the last three engine slices showed; the probit
half ~195 of them. Stops: [Budget and stops](#budget-and-stops).

## Summary for the maintainer

Logistic fits draw k by default, and their k mixes very slowly near separation: tens of thousands of
iterations per independent draw on the worst simulated datasets, and unconverged at 64,000 draws at the
sizes users fit. The probit step built on 2026-10-09 cannot help logistic. This plan builds the move the
survey measured: every leaf multiplied by one factor and k divided by it, the factor drawn from the
binary likelihood with the latent variables averaged out, then the latent variables drawn again. It is
exact and runs on about one iteration in four. On the survey's small census it took logistic's 99th
percentile of k's autocorrelation time from about 23,000 iterations to 115. You ruled that binary fits
keep a drawn k (dec-B420), so this half is built.

What is still yours to decide:
1. The probit half, beside the built step: on the survey's small datasets it cut the remaining slow tail
   about five times (99th percentile 879 to 160) for 12 to 15 percent more time per iteration, measured
   on a loaded laptop with the same fast likelihood routine the build will use. The survey's largest
   test, 5,000 rows, gained up to two times or nothing. Before you rule, the plan measures fits of 10,000
   and 50,000 rows with 200 trees: what the move costs there and what it buys. If it buys nothing there,
   it can be switched off above a size, or taken less often at large sizes.
2. Whether the move has an off switch. The built step has none, but it costs 0 to 3 percent; this one
   costs several times that.
3. Negative binomial: recommended to wait for the redesigned study of its default, which you placed after
   the merge to main (dec-B432).

## Goal

On a single constant-leaf logistic forest with a drawn k (and, with the probit half, a probit one), each
sampling sweep takes, with probability 1/4, a move that multiplies every occupied leaf by a factor r and
divides k by it, r drawn by one slice step from its conditional with the latent variables integrated out,
then redraws the latents at the new fit. Logistic's k reaches the survey's tail (99th percentile of tau(k)
about 115 sweeps on the census) and passes a frozen-tree exact gate and start-from-truth; a logistic-k SBC
arm becomes affordable and joins the matrix. No R, bridge, C API or state-format change unless the
maintainer asks for a switch (open call 3). The likelihood the move evaluates runs over the rows in the
likelihood and no others ([The move](#the-move)).

## Context

- The survey that measured it: [k-tail-samplers.md](../design/k-tail-samplers.md), its "collapsed
  rescaling move" (candidate 4) and its Recommendation. Its prototype is scratch/survey/prototype.patch
  (the engine of 2026-10-09, before the k-internal slice; read, not applied), its drivers and results
  scratch/survey and scratch/survey/out. Its every-fourth-sweep censuses ran 2,000 sweeps of burn-in.
- The k study: [binary-k-fixed-or-drawn.md](../design/binary-k-fixed-or-drawn.md) (fixed k is not good
  enough; logistic's k unconverged at every length tried, split-Rhat 1.41 to 1.64 at the probe length;
  probit's at n = 500 and 75 trees too, 1.12 to 1.21). Ruled on as dec-B420: a drawn k under
  chi(1.5, 2) is kept for both links.
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
  [`ProbitResponse::computeLogLikelihood`](../../src/bartcore/model.hpp) (through `Rf_pnorm5` on the log
  scale) and [`LogisticResponse::computeLogLikelihood`](../../src/bartcore/model.hpp) (through
  [`logOnePlusExp`](../../src/bartcore/model.hpp)).
- Rulings: dec-B404 (the logistic move, built if a drawn k is kept), dec-B420 (a drawn k kept under
  chi(1.5, 2) for both links, so the logistic move is built), dec-B415 (the survey), dec-B371 and
  dec-B394 (the probit step and its switch), dec-B397 (no switch for the built step, ruled on its 0 to 3
  percent cost; an engine toggle for tests only), dec-A190 (the probit-k SBC arm at thin 700 "for now";
  "We'll revisit everything k-related as we gather the evidence."), dec-A191 (the engine order), dec-B412
  (negative binomial's prior predictive study), dec-B431 (that study read: the negative binomial default
  stays k drawn under chi(1.5, 2), the study judged flawed and the default improvable), dec-B432 (a
  redesigned study before 1.0-0, after the merge to main), dec-B390 (the binary hyperprior rerun).
- What the move is written against, landed since this plan was first drawn: k as the chain's parameter
  against a scale the family fixes ([k-internal-parameter.md](k-internal-parameter.md), its
  [Landing note](k-internal-parameter.md#landing-note); dec-B416 to dec-B418, dec-A196), an install that
  puts a state in as stored ([setstate-force-update.md](setstate-force-update.md)), `run()$k` recorded
  exactly when k is drawn in that run (dec-A202), and the rows in the likelihood
  ([response-scale-rows.md](response-scale-rows.md)).

## The move

One forest, constant leaves, k ~ s chi(nu), each occupied leaf mu ~ N(0, (c / k)^2). k is the chain's
parameter against a scale the family fixes: for a binary family c sqrt(m) is 3 (probit) or pi sqrt(3)
(logistic) whatever the prior's spelling, and no re-anchor moves it, so nu and s are read off the forest's
k prior as it stands (`kHyperprior` on [`Forest`](../../src/bartcore/combiner.hpp)), an sd-spelled prior
included, which the engine holds as chi(nu, `k.scale` / s). With the latents
integrated out the target over (trees, leaves, k) is the prior times the binary likelihood. The group
r > 0 acts as (mu, k) -> (r mu, k / r) on the M occupied leaves and k; eta = k mu / c and the trees are
invariant, and so is every fit up to the factor. With f_i the forest's cached fit
([`Forest`](../../src/bartcore/combiner.hpp)'s `totalFits`, the fit without the offset), o_i the
offset, the sums over the rows in the likelihood, C = k^2 / (2 s^2) and v = log r, the generalized Gibbs
step of Liu and Sabatti (2000) on Haar measure dr / r gives the Jacobian r^(M - 1) against the leaf prior's
r^-M and k's prior's r^-(nu - 1) exp(-C / r^2):

    log p(v) = -nu v - C e^(-2v) + sum_i l_i(e^v f_i + o_i)

    probit:    l_i(eta) = log Phi(s_i eta),                    s_i = 2 y_i - 1
    logistic:  l_i(eta) = w_i (y_i eta - log(1 + e^eta)),      w_i the count weight (1 unweighted)

M drops out with the leaf prior, as in the built step. A check with no data: under the prior eta is
independent of k, so k' = k e^-v must be s chi(nu), which gives exactly e^(-nu v) exp(-C e^(-2v)) (the
critique's cross-check). Empty leaves are fixed at zero and stay there; a leaf whose members are all
inactive is drawn from its prior and scales like any other (the leaf walk of the built step, so the set
scaled is the set the k draw counts). The offset is data and is not scaled. Probit carries no weights
(the family refuses them), so its sum has none. The term C e^(-2v) is evaluated as 0.5 (k e^-v / s)^2,
so it cannot underflow to zero at a small k and leave the slice an improper target.

Which rows. A row is in the likelihood when its weight is positive and no mask switches it off, the
sense [response-scale-rows.md](response-scale-rows.md) fixed. For these two families that is the response
model's own active-row test and nothing else, the test the latent refresh and `computeLogLikelihood`
already apply: probit carries no weight channel, and a 0/1 weight vector given to it installs as the mask
([`isMaskedWeightFamily`](../../R/spec.R)); a logistic weight is a count, refused unless a positive integer
at creation and at `setWeights` ([`refuseNonCountWeights`](../../R/spec.R)), so no logistic row sits at
weight 0, and a count enters as the multiplier w_i. A row out contributes no term to the sum, keeps its
leaf's occupancy (so its leaf is scaled with the rest) and keeps its latent as it stood. The sum is never
taken over every row while a mask is installed: [Mutants](#mutants), 5. Neither family derives anything
from the response, so no response range, center or kept scale enters the move, and a `setResponse` or
`setOffset` on these families always answers `accepted`
([`ResponseUpdateResult`](../../src/bartcore/sampler.hpp)).

How the latents leave: under probit the integral of N(z; eta, 1) over s z > 0 is Phi(s eta); under
logistic the Polya-Gamma identity (Polson, Scott and Windle 2013) integrates omega ~ PG(w, 0) out of
exp(kappa eta - omega eta^2 / 2) 2^-w, kappa = w (y - 1/2), to exactly the binomial likelihood above. No
approximation enters; the conditional is evaluated with the likelihood itself.

The probit evaluator. log Phi is evaluated from the start by the erfc form with an asymptotic tail, the
form the survey's costs were measured with, not by `Rf_pnorm5`, which the critique's microbenchmark put at
24.8 ns a row against the erfc form's 16.1 (1.5 times slower) and which would make the slice's 6 to 10
passes, most of the move's cost, half again as dear. With x = s_i eta_i:
- x >= 0: log1p(-0.5 erfc(x / sqrt 2)), so a probability near 1 keeps its relative accuracy (the
  prototype's log(0.5 erfc(-x / sqrt 2)) loses it there: at x = 8 it is several percent off);
- -37 < x < 0: log(0.5 erfc(-x / sqrt 2)), finite to x = -37 (erfc underflows near -38);
- x <= -37: -x^2 / 2 - log(sqrt(2 pi)) - log(-x) + log1p(t), t the asymptotic series
  -1/x^2 + 3/x^4 - 15/x^6 + ... carried until a term falls below 1e-17 (seven terms at x = -37; the
  prototype's three leave an error near 1e-8 in log Phi there, far outside the tolerance below).
It must agree with `Rf_pnorm5(x, 0, 1, 1, 1)` within 1e-14 relative on a grid of x from -60 to 40 and at
each branch's edges (tests/cpp 1). `ProbitResponse::computeLogLikelihood`, which a fit reports, keeps
`Rf_pnorm5`; the move's evaluator and it then differ in the last bits, which leaves the move exact to
floating point and no reported value changed. Logistic needs no such form: its per-row term through
[`logOnePlusExp`](../../src/bartcore/model.hpp) is already cheap. Both stay finite at large r.
`std::erfc` is new to the draw path (exp, log, log1p, expm1, lgamma and pow are on it today and the
recorded baselines reproduce bit for bit on CI's macOS runner): whether the runner's erfc returns the
laptop's bits is not known until the re-recorded probit scenarios run there, and is stop 9 of
[Budget and stops](#budget-and-stops).

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
the orbit coordinate given the orbit-invariant ones (trees, eta, y, the offset, the mask), never reading
z, so the moved x' is again such a draw; drawing z' from its exact conditional given x' before anything
else reads z restores the joint (van Dyk and Park's condition for a partially collapsed step). Every
decline below reads only orbit-invariant quantities (the family, the forest's shape, whether k is drawn,
finite and positive, whether a tree map is stale, whether the prior's scale is infinite), and the coin
that schedules the move is independent of the state, so mixing the move with the identity keeps the
posterior. Under an infinite prior scale (C identically 0) log p grows without bound as v -> -infinity
(k -> infinity, the fit -> 0, the likelihood -> its value at the offset), so the conditional is improper
and the move declines; the package already documents k's posterior as improper there
(test-binaryResponse-hyperprior.R). For a finite scale and nu > 0 the conditional is proper:
-C e^(-2v) dominates as v -> -infinity, and -nu v with a likelihood term at most 0 as v -> +infinity.
With every row inactive the likelihood term vanishes and the conditional is k's prior along the orbit,
still proper: the move runs (unlike the built step, which declines there at R = 0), redraws no latent, and
moves k and the leaves together as the prior allows. A mask with every row out installs from R today
(`setActiveRows(rep(0, n))`, ran), so this case is reachable. One decline falls after the coin: a
log p(0) that is not finite, which only a state of zero prior density gives (k / s overflowing under a
scale near the smallest double, which `chi()` accepts); it costs no exactness.

Inside a larger sampler. The move is a kernel for (leaves, k, latents) given everything else the chain
holds, so it is valid as part of the BART block wherever the family's own latent refresh is: the code
that runs the BART sampler inside a larger sampler must reach the binary response only through y, the
offset, the counts and the mask, which the move reads at the sweep it runs.

What the prototype got wrong or left out, so the build does not copy it: probit declined under a mask;
logistic ignored the count weights and the mask; nothing declined under an infinite scale (the slice
wandered an improper conditional); the schedule was an environment variable with a counter that
setResponse reset; the signs were read from the latents or the working response, not from y; and its
erfc-based log Phi lost accuracy near Phi = 1 and in its three-term tail (fixed above).

## Where it runs

In [`runSweeps`](../../src/bartcore/chain.hpp), after the built step and before the tree loop: level step,
built step, then the collapsed move, the survey's measured order. Every recorded channel (training and
test fits, saved trees, k) is written after it, as for the built step, and the tree loop reads the
redrawn latents' working response and, under logistic, their precisions as weights.

Eligible when all hold, read before any generator call: the engine toggle is on; one forest, no combiner,
no variance forest; the plain constant leaf; family logistic (and probit, with the probit half); the forest
draws k, k finite and positive; no tree's map stale; the k prior's scale finite (tested on the scale
itself, never on C, which depends on k). Otherwise nothing is drawn. An eligible sweep draws one uniform;
below 1/4 it takes the move: one pass for log p(0), the slice step's passes, the leaf walk, k, the cache,
then one latent refresh.

The schedule is a coin, not a counter ([Settled in planning](#settled-in-planning), 1): on average every
fourth sweep, as the survey measured, with no state of its own.

With the built step (probit half): both act on one forest in one sweep. They share the leaf scaling and the
cache rule through one helper extracted from [`Chain::drawForestRescale`](../../src/bartcore/chain.hpp):
occupied leaves times the factor, k divided by it, `totalFits` multiplied in place while the forest's
accumulated factor (`totalFitsScale`) stays within [1/2, 2] and re-derived from the leaves past that. Both
moves multiply into the one accumulated factor, so the gap the cache can carry stays bounded as today.
The built step then scales its latents as it does now; the collapsed move redraws them, so on a sweep that
takes both the built step's latent scaling is overwritten, while its leaves and k stand. Extracting the
helper is a neutral commit of its own: every equivalence scenario bitwise. `drawForestRescale` has the
shape the extraction assumes, unchanged by the slices landed since: after the slice step it walks the
leaves, divides k, calls `scaleLatents`, then scales or re-derives the cache, and draws nothing after the
slice step, so the helper (leaves, k, cache) with the latent scaling on either side of it draws bitwise.
[`Chain::forestRescaleApplies`](../../src/bartcore/chain.hpp) holds the conditions both predicates share
(one forest, no combiner, no variance forest, the plain constant leaf, k drawn, finite and positive, no
stale map); the same commit splits them from its family test so the two moves cannot disagree on them.

## State, chains and threads

- No new chain state. The coin reads the chain's generator, whose state rides every stored state
  ([`ChainStateData`](../../src/bartcore/combiner.hpp)'s `rngState`), so after `setState`, `copy()` or a
  `saveRDS` reload the coin falls on the sweeps it would have in the uninterrupted chain. A restore does
  not continue that chain bit for bit and never did: an install rebuilds the fits from the flattened trees,
  so a copy, a reload and a sampler's own restore agree with each other bit for bit and with the
  uninterrupted chain to rounding (the k-internal slice's
  [Landing note](k-internal-parameter.md#landing-note)). Ran on the base, drawn-k probit and logistic,
  restored after 3 sweeps and after 30: the tree moves identical over the 5, 40 and 200 sweeps that
  followed, the fits and k within 2e-14 relative; a copy and a reload identical. The move adds nothing to
  that: an install leaves the forest's accumulated factor at 1 and the cache re-derived, and the move's
  own rebuild past the bound is the same re-sum. A warm start is a new fit with its own generator and
  carries nothing for it.
- What a stored state must carry for the move, and does today: each forest's leaves and its k (written
  only where k is drawn), the latents, the generator. The move redraws the latents it reads past, and
  reads no other state. The mask is not part of a state and R puts it back after an install
  ([`reapplyActiveRows`](../../R/dbarts.R)); the counts and the offset are data. No state field, no state
  element in the bridge, no format version ([`stateFormatVersion`](../../src/R_interface_bartcore.cpp)).
  An install goes in as stored (dec-B418) and a binary family's mapping is fixed, so neither form of
  `setState` has anything to convert or report for the move.
- setResponse, setOffset, setWeights, setActiveRows and setData between runs need nothing: the move reads
  y, the offset, the counts and the mask from the response model at the sweep it runs. `setLeafPrior`
  between runs, which can turn a drawn k into a fixed one and back, is read at the sweep too (the forest's
  `updateK`); `run()$k` is there exactly when k is drawn in that run (dec-A202), which is when the move
  can run.
- Chains: each chain decides and draws on its own generator; two chains take the move on different sweeps.
  Threads: the move runs inside its chain's sweep with no shared scratch (it needs none: the evaluator is a
  fused loop over the chain's own arrays), so n.threads 1 and 2 draw bitwise the same.
- The likelihood sum is a plain scalar loop in row order in both the shipped and the reference build. It
  must not be vectorized or reassociated: the sum's order is part of the draw, and the equivalence
  baselines are reference-build values the shipped build reproduces bitwise today.
- A cancel cannot land inside the move: the refresh-cancel poll is armed only around the sweep's own
  latent refresh ([`armRefreshCancel`](../../src/bartcore/chain.hpp)), and the move's refresh runs outside
  it, so a canceled run never leaves the leaves scaled and the latents stale.

## Scope

Reached through, where k is drawn (the default for binary families,
[`drawsLeafKByDefault`](../../R/spec.R)) with a finite scale:
- logistic: `bart`, `dbarts` and `xbart` on a 0/1 response with family logistic (`xbart` where its k is
  left at the default or given as `chi()`), count-weighted fits and fits under a mask or an offset
  included, and `hazard.logistic`, the hazard family's third spelling, which is a logistic fit on the
  person-period rows (one forest; its subject weights are counts, so a subject at weight 0 is refused).
- probit, with the probit half: `bart`, `dbarts` and `xbart` on a 0/1 response, a 0/1 weight vector
  included (it is the mask), `hazard` and `hazard.probit` (one forest), the zero part of
  `hurdle.lognormal`, `rbart_vi`'s binary fit (through its offset), `bartBT` given `k = chi(df, s)` with
  s finite, and through the flat C entries a binary stan4bart fit, which is probit, draws k unless its
  `bart_args` fix it, and carries its parametric part as the offset (read in stan4bart, not run). The same
  fits the built step reaches, less those under an infinite scale.

Each ran on the base: `bart` and `hazard.logistic` under logistic and `bartBT` with `k = chi(1.5, 2)`
return a k that moves; `xbart` on a binary response runs under both links with k left out and with
`chi(1.5, 2)` (that its default there is the drawn k was read, not run); a probit fit under 0/1 weights
mirrors the mask and holds no weights; a logistic sampler takes new counts by `setWeights`; an sd-spelled
logistic prior draws k against a `k.scale` of pi sqrt(3).

Bitwise unchanged, consuming no draw: every other family (gaussian, Student-t, aft, ordinal, multinomial,
negative binomial); a fixed k (bartBT's default, the SBC's fixed-k logistic arm, logistic-reference.R,
`xbart` over a numeric k grid); an infinite scale, chi(df, Inf), which `sd = invchi(df, 0)` also names
(equivalence's chik, logistic and wtlogistic; test-binaryResponse-hyperprior.R's chi(1, Inf) and
chi(1000, Inf) fits); BCF and every multi-forest or amplitude-coupled fit, the latent BCF families
included; linear, gp and monotone leaves; grow-from-root's sweeps (an initializer); a sweep whose tree
maps are stale. Without the probit half, every probit fit is bitwise unchanged too.

Why each is out. Several forests: the move as derived scales one forest's leaves against its own k and
reads that forest's cached fit as the whole location; under a combiner the location is a blend of
forests with amplitude draws between them, which is not derived and not measured, as for the built step.
Multinomial is several forests by construction. Ordinal and multinomial draw no k by default
([`drawsLeafKByDefault`](../../R/spec.R)), a drawn k's mixing there was never measured, and each would
need its own likelihood in the evaluator. A weighted probit does not exist (the family has no weight
channel). Counts under logistic, a mask and an offset are all in.

**The fixed-k case.** The predicate requires a drawn k (the forest's `updateK`), so a fixed k declines
before the coin and draws exactly as today; so does a k that `setLeafPrior` fixes between runs, from the
next sweep. A k that runs to infinity (an infinite scale with a large df) declines on finiteness, and so
does a k of exactly 0, which a prior scale near the smallest double can draw.

**Negative binomial: not in this slice** ([Open calls](#open-calls), 4). The move applies in principle
(the likelihood is cheap, the latents are Polya-Gamma), but: it was not measured, and dec-B404's option
was to measure first with a stop if it does not help; its k is not the slow part (99th percentile of
tau(k) 3,253 sweeps where measured, against logistic's 27,505; "For negative binomial alone, k is not the
priority"); its default prior was studied (dec-B412) and kept at chi(1.5, 2) with the study judged flawed
and the default improvable (dec-B431), and a redesigned study runs before 1.0-0, after the merge to main
(dec-B432), which may still change what k it draws; and
its latent redraw is a PG(y_i + r) pass whose cost grows with the counts, and its refresh redraws the
shape too, so the move's cost there is unknown. The code leaves it a door: an evaluator in
[`NBResponse`](../../src/bartcore/model.hpp) and the family in the predicate, plus a census arm and an
exact-gate arm with an exact negative binomial frozen-tree posterior. It keeps a TODO entry
([Changes](#changes), docs).

## User surface and API

As planned, none: no `dbartsControl` argument, no bridge change, nothing in dbarts.h (no ABI event), no
state block, no help page beyond NEWS. The engine keeps a toggle, `SamplerOptions::collapsedRescaleForest
= true`, commented as [`probitRescaleForest`](../../src/bartcore/chain.hpp) is: not a setting, false only
in the engine's own tests. Whether to ship a switch is open call 3; if the maintainer asks for one, the R
surface follows the built step's withdrawn switch (a logical `dbartsControl` argument fixed at creation,
the bridge's control parse, help, the surface tests), about 120 lines more. Nothing a fit reports changes
meaning: `fit$k`, `run()$k` and `getK` read the same k against the same `k.scale`.

## Changes

Engine ([chain.hpp](../../src/bartcore/chain.hpp), [model.hpp](../../src/bartcore/model.hpp); names are
proposals, cited by symbol once landed):
- Constants beside the built step's: `collapsedRescaleProbability = 0.25`,
  `collapsedRescaleSliceWidth = 0.5`; the step limit reuses `rescaleSliceStepLimit`.
- `SamplerOptions::collapsedRescaleForest`, as above.
- The helper extracted from `drawForestRescale` (leaves, k, cache), used by both moves; the built step
  keeps its latent scaling and its order of draws. The conditions the two predicates share, split out of
  `forestRescaleApplies`.
- `Chain::collapsedRescaleApplies()`: the predicate, no generator. `Chain::maybeDrawCollapsedRescale()`:
  predicate, coin, move; returns whether it moved. `Chain::drawCollapsedRescale(double* factorOut =
  nullptr)`: the move without the coin, factorOut receiving r (1 on a decline), the out-parameter
  tests/cpp reads, as `drawForestRescale`'s `alphaOut`.
- [`ResponseModel`](../../src/bartcore/model.hpp) gains
  `virtual double logLikelihoodAtScaledFits(const double* fits, double factor) const`: the sum over the
  rows in the likelihood, in row order, of the row's log-likelihood at factor * f_i + o_i.
  `LogisticResponse`'s uses the per-row expression of its `computeLogLikelihood`, count weights as
  multipliers. The two agree to rounding and not bit for bit: at R's flags on arm64 the compiler fuses the
  factor's multiply into the offset's add (ran), which a caller that scales the fits first cannot
  reproduce, so tests/cpp 1 compares within 1e-14 relative under both links. `ProbitResponse`'s (the
  probit half) uses the erfc form of [The move](#the-move), within 1e-14
  relative of its `computeLogLikelihood`, which keeps `Rf_pnorm5`. The base returns NaN, unreachable past
  the predicate (as `scaleLatents`' base refuses). The redraw is the family's existing `refreshLatents`
  at the combined fits, which for one forest are its `totalFits` after the helper has scaled them.
- After the move `runSweeps` re-reads its working-response and weight pointers, as it does after the
  sweep's own refresh (they are stable today; re-reading keeps the invariant local).
- No facade virtual, nothing in dbarts.h, no bridge, no state-format change; `--preclean` on every engine
  commit all the same (model.hpp's virtual table changes).

Scripts:
- [probit-k-scale-exact.R](../../benchmarks/R/probit-k-scale-exact.R) gains a link argument and logistic
  arms ([The exact gate](#the-exact-gate)); the name stays, so the workflow list and its cites stand. Its
  written-in leaf scale of 3 becomes the sampler's own, `getLeafPrior()$k.scale` on the one-tree fit (3
  under probit, pi sqrt(3) under logistic; ran). The counts arm gives the gate's sampler its counts by
  `setWeights` between arms, which a logistic sampler takes (ran), or builds its own sampler.
- [probit-k-mixing.R](../../benchmarks/R/probit-k-mixing.R): a fresh sampler per dataset, seeded by the
  dataset's index, in [`runCensus`](../../benchmarks/R/probit-k-mixing.R) (today one sampler serves every
  dataset and the previous chain advances its generator, so per-dataset results cannot be compared across
  builds that draw differently; the survey's census.R made this fix and reproduced the published census,
  median 73 against 71), and a link argument for `census` and `truth`. Its section in
  benchmarks/README.md follows.
- [sbc.R](../../benchmarks/R/sbc.R): [`sbcConfigProbitK`](../../benchmarks/R/sbc.R) generalized to a
  family argument, a `logistic-k` arm, [`sbcMatrixConfigs`](../../benchmarks/R/sbc.R) and the matrix
  count M (84 to 91, [`sbcMatrixFunctionals`](../../benchmarks/R/sbc.R));
  [sbc.yaml](../../.github/workflows/sbc.yaml) gains the job and its header's count; benchmarks/README.md
  its line.
- [equivalence.R](../../benchmarks/R/equivalence.R) gains one scenario, `logistick`, its 56th: the sampler
  API, family logistic, the default chi(1.5, 2), count weights 1 to 3 and an offset, so the move's w and o
  terms have a bitwise anchor (every logistic scenario today uses chi(1.5, Inf), which the move declines).
- [binary-hyperprior.R](../../benchmarks/R/binary-hyperprior.R) is used as it stands for the
  pre-registered cells ([Early measurements](#early-measurements)); no change.

Docs at landing: docs/design/collapsed-scale-move.md (the derivation and exactness above, the schedule,
the gates' landed numbers, the early measurements, the x86 cost) and its INDEX row; k-tail-samplers.md's
and k-mixing-pg-families.md's Status lines point to it; [One sweep](../architecture.md#one-sweep) step 2
gains the move; inst/NEWS.Rd's tree-moves item gains a sentence (users' view, against 0.9-34: binary fits
that draw k, `hazard.logistic` among them, take a further step on about one iteration in four that rescales
the leaf values and k together, drawn with the latent variables averaged out, so k and the fit's overall
scale mix far faster near separation; the posterior is unchanged), parse-gated per the README's
checklist; TODO's k-mixing-pg-families item narrowed to negative binomial (measure k's mixing under the
move after dec-B432's redesigned study of its default, and build its arm if it helps), not closed; the
MANIFEST row; the feature matrix's line that names the current baseline; a ledger entry for the calls
made; this plan's Status and Landing note.

## Tests

tests/cpp, in [test_ensemble.cpp](../../tests/cpp/test_ensemble.cpp) beside
[`runForestRescaleTests`](../../tests/cpp/test_ensemble.cpp), through
[`TestPeer`](../../tests/cpp/test_peer.hpp). Fixtures: a burned-in drawn-k sampler per link with an offset,
a mask, two or more trees and a stranded nonzero empty leaf; the logistic one also with counts 1 to 3.
Reference constants (nu, C, the per-row terms) come from the raw state: y, the offset, the mask, the
counts, the leaves re-summed in tree order and the k prior, never the move's own code. The block reuses
what the built step's block holds (its fixture builder, the quadrature target, the KS distance, the
decline-without-a-draw check and the peer's accessors); the builder sets the new toggle off for the built
step's own fixtures, so `runForestRescaleTests` draws as it does today. Every generator a test creates is
destroyed: the k-internal slice's leak-sanitizer leg found four that were not.
1. Evaluator: `logLikelihoodAtScaledFits(f, a)` equals the row-order sum of `computeLogLikelihood` at a f
   over the rows in at several a, within 1e-14 relative under both links (not bitwise:
   [Changes](#changes)); and
   probit's per-row log Phi against `Rf_pnorm5(x, 0, 1, 1, 1)` within 1e-14 relative on a grid of x from
   -60 to 40 and on each side of x = 0 and x = -37. Again after `setResponse` with a new y and after `setData`
   with new rows, so the evaluator reads the response model's current y (the prototype's sign-source bug).
   Under the fixture's mask the sum has exactly the rows in: with one row out given a fit of 1e6 against
   its response, the value does not move.
2. Conditional, in invariance form (a slice step from a fixed point is not a draw from the target): per
   replication draw v0 from the frozen state's target by inverse-CDF quadrature, move the frozen state to
   the orbit point e^v0 (leaves and k only; the move does not read the latents), take one
   `drawCollapsedRescale`, read v1 from the leaves and from k, which must agree. KS of the v1 against the
   target, D < 1.95 / sqrt(N) (p = 1e-3), N = 2e5, or 5e4 if a link's run passes 10 s. Per link, and on a
   near-separated fixture too. The quadrature grid is set per fixture so the target's mass outside it is
   below 1e-12: the right tail decays only as e^(-nu v) where the data are separated, so the built test's
   [-6, 6] does not serve (at nu = 1.5 the grid reaches about v = 20).
3. Redraw: from a fixed state, many moves; under probit each active row's latent sits on y's side and its
   PIT under N(e^v f_i + o_i, 1) truncated to that side is uniform (KS on a few rows, pooled); under
   logistic each sampled omega_i standardized by PG(w_i, psi)'s mean w_i tanh(psi / 2) / (2 psi) (its
   limit w_i / 4 at psi = 0, which a row with zero fit and offset exercises) and variance at the post-move
   psi averages to zero (|z| < 4.5). Inactive rows' latents unchanged.
4. Mapping, one move with r from the out-parameter: k' r = k; mu' = r mu on every occupied leaf of every
   tree, the empty leaf untouched; `totalFits` equal to the tree-order re-sum bitwise when the bound forced
   a rebuild and within 1e-14 relative otherwise, `totalFitsScale` multiplied by r; under probit the latents
   are not r z (redrawn).
5. Inertness: toggle off, fixed k, an infinite scale, k infinite, k zero, a stale map, gaussian, negative
   binomial, ordinal, a two-forest fit and a linear leaf consume no generator draw (the next uniform
   unchanged); an eligible sweep whose coin says no consumes exactly one. A tiny k (1e-170, whose square
   underflows) under a finite scale is eligible and its move finite. A scale of 1e-200 with k at 1, where
   log p(0) is not finite, is eligible and its move declines with no draw of its own.
6. The coin: `maybeDrawCollapsedRescale` called 4000 times on an eligible frozen state moves on a share
   within 4.5 binomial sd of 1/4.
7. Every chain: on a two-chain sampler each chain's k differs from a toggle-off twin after 50 sweeps.
8. With the built step (probit half): one sweep taking both leaves `totalFitsScale` the product of the two
   factors or 1 after a rebuild, and the debug cache check
   ([`assertForestCachesMatchLeaves`](../../src/bartcore/chain.hpp)) silent.
9. Every row inactive: the move runs, its v1 passes test 2's KS against the target with no likelihood
   term (k's prior along the orbit, N = 5e4), and no latent changes.

tinytest, a new file inst/tinytest/test-collapsed-rescale.R, both links unless said:
- Per chain, on any number of threads: a two-chain drawn-k fit draws identically at n.threads 1 and 2.
- Ahead of every recorded channel: a keepTrees fit's saved trees replayed reproduce each kept draw's
  training and test fits within 1e-10.
- Finite under an offset, a mask installed mid-run, a mask with every row out and (logistic) count
  weights, at a finite scale; a probit fit given 0/1 weights draws as the same fit given them as a mask
  (identical today, ran); a chi(1.5, Inf) logistic fit runs; a `hazard.logistic` fit with a drawn k runs
  finite.
- State: a sampler restored by setState from a store and a reload of the same store continue identically,
  as does a copy.
- [test-composition-sequences.R](../../inst/tinytest/test-composition-sequences.R) gains `probitK` and
  `logisticK` arms in [`continuationArms`](../../inst/tinytest/test-composition-sequences.R): restore after 3
  sweeps (not a multiple of 4), then the arm's own assertions against the uninterrupted run, the tree moves
  identical and the fits within 1e-9 of their scale. Both pass those today, before the move (ran on a
  binary response of the file's size and settings).
  The file runs every arm 5 sweeps past the restore, over which a coin of 1/4 takes no move about one
  time in four; these two arms run 24 past it (one in a thousand), their lengths made arguments of
  [`continuationArm`](../../inst/tinytest/test-composition-sequences.R). This is the test a schedule
  counter left out of the state would fail.
- [test-probit-rescale-forest.R](../../inst/tinytest/test-probit-rescale-forest.R)'s header says its probit
  fits take both moves (probit half).

## Early measurements

Run before review (step 7 of [Order of work](#order-of-work)), on the engine of steps 2 to 6 and a scratch
copy whose toggle default is flipped off (a second private library from an uncommitted patch). They feed
the maintainer's calls; review does not start until they are in.

1. **Cost on the x86 box.** The survey's costs came from a laptop at load 20 to 40, and logistic's was
   never measured. The survey timed the same erfc form of log Phi the build uses, so its probit figures
   stand as the expectation: every fourth sweep 1.15 at 150 rows and 50 trees and 1.12 at 5000 rows and
   200 trees, 12 to 15 percent. The method is the built step's: per-sweep time of the slice's build over
   the flipped copy, both restored from one stored state before every timed segment, ten interleaved
   repetitions in alternating order, one chain, one thread; median ratio and range, on the idle box
   (started for the run and stopped after). Cells: 150 rows and 50 trees; n in {500, 2000, 5000, 20000}
   by trees in {75, 200}; and the large cells of measurement 4 (10,000 and 50,000 rows, 200 trees, both
   signals), timed in segments of 100 sweeps; probit (both moves against the built step alone) and
   logistic. The ratio includes the coin, so it is the average cost; one move costs (ratio - 1) / 0.25.
   The table goes to the maintainer with open calls 2 and 3; the probit ruling waits on it.
2. **The sizes the k study ran** (pre-registered here). [binary-hyperprior.R](../../benchmarks/R/binary-hyperprior.R)
   at the k study's probe length (`BINARY_HYPERPRIOR_CHAINS=8`, `BINARY_HYPERPRIOR_BURN=8000`,
   `BINARY_HYPERPRIOR_DRAWS=8000`), the incumbent arm only (`BINARY_HYPERPRIOR_ARMS="chi(1.5, 2)"`), four
   repetitions, with the move and on the flipped copy:
   - logistic (`BINARY_HYPERPRIOR_LINK=logistic`): the nine near-separable cells, `sim:separable:500`,
     where the study read k's split-Rhat at 1.55 to 1.64;
   - probit: the same nine and the three strong-signal cells with 20 predictors, `sim:strong:500:20`,
     where the study read 1.12 to 1.27 with the built step alone.
   Reported per link and build: the share of fits with k's split-Rhat above 1.05, k's median and 10th
   percentile effective draws, and the held-out probability's split-Rhat share. The reading, fixed now:
   the move fixes what the study found if, with it, k's split-Rhat is at most 1.05 in at least three
   quarters of the fits and k's median effective draws reach 200, September's reopen condition; it
   helps if the share above 1.05 halves; otherwise it does not reach these sizes. Not a stop: the result
   goes to the maintainer beside the open calls, and the probit reading is the evidence for open call 2's
   n = 500 question. About 1.5 hours of two laptop cores.
3. **The logistic gate's mixing bound on both machines.** The gate's logistic pure arm in quick mode, three
   seeds, with the move and on the flipped copy, on the arm64 laptop and on the x86 box (exact-gates runs
   on x86-64 Linux; the probit bound needed both, 0.0064 against 0.0099). Separated when the largest
   with-move error over machines and seeds, times 1.5, is below the smallest error without it; the quick bound
   is then set at the geometric mean of those two and the full bound likewise from full runs. Not
   separated is stop 8.
4. **Large fits** (pre-registered here). The worry is a per-sweep tax on every default probit fit that
   buys nothing at large n: the survey's largest test, 5000 rows and 200 trees with one dataset a cell,
   gained up to two times or nothing, and n = 500 or 2000 are not large fits.
   - Data: the survey's largen.R generator, 10 uniform predictors, eta = a (2 sin(pi x1 x2) +
     2 (x3 - 0.5)^2 + x4 - 1), y = 1{eta + e > 0} with e from the link's own noise (normal for probit,
     logistic for logistic), a set per link so the sign of eta matches y on about 90 percent of rows
     (near-separated) or about 76 percent (moderate), the survey's two signals; the realized share is
     reported.
   - Cells: n in {10,000, 50,000} by the two signals by the two links, 200 trees, the default chi(1.5, 2),
     one chain, one thread; two datasets a cell (seeds 1 and 2); each with the move and on the flipped
     copy (probit: both moves against the built step alone; logistic: the move against none). A 20,000-row
     cell is added for a link only if the reading below finds 10,000 paying and 50,000 not, to place the
     boundary.
   - Length: 5,000 sweeps of burn-in and 40,000 recorded at thin 10, in chunks of 500 kept draws so that
     no run holds more than 500 draws of 50,000 fitted values; per chunk k and the mean training fit are
     kept. A cell whose tau(k) without the move passes 1,000 sweeps is rerun at 120,000 recorded.
   - Measured per cell, as geometric means over its two datasets: the gain, tau without the move over tau
     with it (Geyer's estimate, as the mixing script's), for k and for the mean fit; the cost ratio from
     measurement 1 at the same cell; and the net, gain over cost ratio, effective draws per second with
     the move over without.
   - Reading, fixed now: a cell buys nothing when the gain is below 1.25 for k and below 1.1 for the mean
     fit while its cost ratio is above 1.03; it pays when either net is at least 1.
   - Decision rule for the probit half (open call 2): every large cell pays, include it as planned; the
     large cells buy nothing while the census and the survey's 5000-row cells gained, include it with one
     of the two large-n options of open call 2, the threshold at the smallest measured size whose cells
     buy nothing at both signals; near-separated cells pay and moderate ones buy nothing, the frequency
     option is the better fit (the size cannot see the signal). The same reading is reported for
     logistic: its move is ruled, but where it buys nothing at large n the same option is put to the
     maintainer for it too.
   - Compute: about 6 ms a sweep at 10,000 rows and 30 ms at 50,000 under probit (the survey's 2.55 ms at
     5000 rows and 200 trees, scaled), half again under logistic; 45,000 sweeps a run make 4.5 to 34
     minutes; 8 cells, 2 datasets and 2 builds come to about 9 core-hours, 4.5 hours on two laptop cores
     (tau is in sweeps, so the laptop's load does not bias it). Their x86 timing adds about 30 minutes of
     the box to measurement 1.

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
  the two moves together is safe across machines; a probit half that silently declines is caught by tests/cpp
  6 and 7 and the census, not here.
- logistic: from early measurement 3, recorded in the header with both machines' values, so a logistic move
  that silently declines fails per push.
Run as `Rscript benchmarks/R/probit-k-scale-exact.R quick` (both links) in
[exact-gates.yaml](../../.github/workflows/exact-gates.yaml)'s existing list, which stays at 42 commands.
Timed on the base, before the move, on the gate's own fixture on the laptop at load 4 to 7: 1e6 sweeps
take 5 to 9 s under probit, 27 to 29 s under logistic and 52 to 58 s under logistic with counts 1 to 3.
The five logistic arms are then about 3 minutes in quick mode and 12 in full before the move's own cost,
which the coin's share of a slice step and a second Polya-Gamma pass puts near a third more; timed again
at step 6. Quick stays inside the job's 40 minutes with the other gates (if the logistic arms
pass 15 minutes on the runner, their quick sweeps halve and the bound is re-read); full runs in the weekly
job, whose per-gate limit is raised for it if needed. The survey's prototype passed quick with the probit
move every sweep at worst |z| 2.9, the pure arm's error 0.0031; it never ran a logistic arm.

### Start-from-truth

The mixing script's `truth` action per link: an exact posterior draw run on, log k at lags 300 to 3000;
|z| <= 3 at R = 2000. Survey: probit +0.57, logistic +1.23.

### The mixing gate (pre-registered)

The census per link, 100 datasets of the probit-k arm's generator (150 rows, 3 predictors, 50 trees, k from
chi(1.5, 2); logistic the same with a logistic response, as the survey), each on its own sampler, 30,000
sweeps of burn-in (the survey's every-fourth runs used 2,000) and 200,000 recorded at thin 10. tau in
sweeps; the bar is 1.5 times the survey's every-fourth-sweep figures:

| link | survey: tau(k) median / 90th / 99th / worst | tau(fit) worst | bar: tau(k) median, 99th | bar: tau(fit) |
|---|---|---|---|---|
| probit, both moves | 35 / 116 / 160 / 166 | 114 (99th 112) | 55, 240 | 99th 170 |
| logistic | 30 / 61 / 115 / 127 | 59 | 45, 175 | worst 90 |

Why 1.5. The survey's schedule was fixed; this one is a coin. If one move nearly decorrelates k (lag
correlation a after a move) and nothing else does, a fixed fourth-sweep schedule gives tau about
4 (1 + a) / (1 - a) and the coin (7 + a) / (1 - a), a ratio of (7 + a) / (4 (1 + a)): 1.75 at a = 0, 1 at
a = 1. The census says the move is not alone: taking it every sweep rather than every fourth cut k's
median only from 35 to 24 and its 99th percentile from 160 to 124, so the other directions dominate and
a sits well above 0, where the ratio is an estimated 1.3 or less; 1.5 adds the noise of the slowest one or two
datasets. Missing the bar is stop 2, and a coin that misses it where a counter would pass returns to the
maintainer ([Settled in planning](#settled-in-planning), 1).

The run is in sweeps, so a laptop under load serves; each census took the survey 19 to 28 minutes at
2,000 burn-in, so about 30 to 40 minutes at 30,000. The census's per-dataset lag-thin autocorrelations also
set the logistic-k arm's thin and record what the probit-k arm could run at.

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

On the reference build against the MANIFEST's current baseline at the slice's base: equivalence-308c45be,
55 scenarios, of which five fit a single probit forest with a drawn k (chik, maskprobit, hazard, hurdle,
bart2probit) and none a logistic forest with a drawn k at a finite scale.
- After the helper's commit (step 2), and after the logistic half (step 3): 55 of 55 identical draws
  (same RNG stream). The five probit scenarios run the helper every sweep, so they are what shows the
  extraction neutral; no recorded scenario takes the logistic move.
- With the probit half: 51 of 55 identical. The movers are exactly maskprobit, hazard, hurdle and
  bart2probit, the single-forest probit fits with a drawn k at a finite scale. chik takes the built step
  only (chi(1.5, Inf)) and stays bitwise; logistic and wtlogistic (chi(1.5, Inf)) and nbinom stay bitwise.
- Without it at the landing: every existing scenario identical, and the re-record only adds a scenario.
- `logistick` is new: recorded fresh, reported skipped by the old baseline.
- bcf-equivalence-1b7d730c 15 of 15 and multinomial-equivalence-80b1c8d4 11 of 11 bitwise (two forests
  decline); both keep their recordings.
- Re-record as equivalence-<sha>, 56 scenarios: the movers and `logistick` recorded fresh and merged into a
  copy of the base's baseline in scenario order, a `--bitwise --strict-coverage` compare against the merged
  file from a second `--preclean` reference install, the MANIFEST row, and the new name written where the
  old one is pinned: cpp-tests.yaml, equivalence.yaml, mutation-battery.R and the feature matrix. ORACLE
  (P17): the exact gate's quick and full passes on both links and the start-from-truth statistic. The z
  sizes the move and is not the oracle; the built step's re-record met |z| above 4 from its old baseline's
  short, unconverged chains in k, and the same reading may recur. The bitwise compare then runs on CI's
  macOS runner: the probit movers reproducing there is stop 9.

### Snapshots

The four reference-build snapshot files pass unchanged (27 results): test-reproducibility-binaryResponse.R
fits through bartBT at a fixed k = 4.5, the xbart file over a fixed k grid, the other two gaussian. A
snapshot that moves is a stop.

### The rest

tests/cpp plain and under address and undefined-behavior sanitizers, and on the x86 box under the leak
sanitizer too; the R-loaded ASAN path over the new and touched test files; the full tinytest suite (a
seeded expectation on a drawn-k binary fit that moves is re-pinned with its reason, a statistical one must
pass as is; 42 of the suite's files fit a logistic family or name chi(1.5, 2), and how many pin a seeded
value on such a fit is not counted); in tests/cpp only the built step's block builds a binary forest with
a drawn k (the eight other `updateK = true` sites are gaussian, read), so no other fixture there takes
the move; every exact-gates.yaml gate in
quick mode, 42 commands, and this one in full (hazard-reduction and hurdle-reduction move both sides
together and must stay bitwise equal, hazard-reduction on both links; logistic-reference and
bcf-latent-exact fix k or use two forests and draw bitwise); `R CMD check --as-cran` from a clean tarball;
lintr, air, rc-codoc, win-drift, doc-freshness; bench-sampler compare (maintainer): with the probit half
["run-binary-n1000-p10-t75"](../../benchmarks/R/bench-sampler.R) takes the move and is expected past 1.05
by its own cost, so its MANIFEST row gains that note or the maintainer re-records the speed baseline.
Consumers: stan4bart's, bartCause's, treatSens's and bairrtt's suites against the build; a seeded
expectation that moves on a drawn-k binary fit is re-recorded on that package's lockstep branch, the
orchestrator pushing. The logistic half reaches none of stan4bart's fits (its binary fits are probit);
the probit half reaches every binary stan4bart fit that draws k, one coin a BART sweep.

## Mutants

Each must fail a named test or gate:
1. The density in r taken against dr instead of dr / r, giving -(nu - 1) v: tests/cpp 2, the exact gate.
2. The Jacobian dropped, giving -(M + nu - 1) v: tests/cpp 2, the exact gate.
3. The C term dropped, with the wrong sign of v, or computed as C e^(-2v) from a C that underflows: tests/cpp
   2 and 5 (the tiny-k fixture).
4. The offset scaled with the fit (e^v (f + o)): tests/cpp 1 and 2, the gate's offset arms.
5. Inactive rows in the sum: tests/cpp 1 and 2, the mask arms.
6. Logistic's count weight dropped: tests/cpp 1 and 2, the counts arm.
7. Probit's sign flipped for failures, or logistic's y flipped (y eta replaced by (1 - y) eta): tests/cpp 1
   and 2, the gate.
8. The evaluator reading a y cached at creation: tests/cpp 1 after setResponse.
9. The latents not redrawn, redrawn at the pre-move fit, or (probit) scaled by r instead: tests/cpp 3 and
   4, then the gate and start-from-truth.
10. k not divided by r: tests/cpp 4 (the gate is blind to this, as it was for the built step).
11. Only tree 0's leaves scaled, or the empty leaf scaled: tests/cpp 4.
12. `totalFits` multiplied without accumulating `totalFitsScale`, or never rebuilt: tests/cpp 4 and 8, the
    debug cache check.
13. The coin drawn before an ineligible forest declines, the infinite-scale guard removed, or the guard
    tested on C: tests/cpp 5, chik, logistic and wtlogistic bitwise.
14. The coin at 1 or at 0: tests/cpp 6; at 0 also the logistic mixing bound and the census.
15. The move placed after the tree loop: the keepTrees replay.
16. Chain 0 only, or scratch shared across chains: tests/cpp 7, the thread-count tinytest.
17. The helper extraction reordering the built step's draws, or the predicate's split dropping a shared
    condition: the five probit scenarios' bitwise compare at step 2 (chik's once the probit half is in),
    and the built step's own inertness checks in tests/cpp.
18. A schedule counter in place of the coin, not stored: the continuation arms.
19. A decline added for every row inactive: tests/cpp 9.
20. The erfc form's x >= 0 branch as log(0.5 erfc(-x / sqrt 2)), or its tail cut to the prototype's three
    terms: tests/cpp 1's grid.
21. A log p(0) that is not finite handed to the slice step: tests/cpp 5's last case.

Read, not run: the likelihood sum stays a scalar loop in row order, and the move's refresh sits outside the
armed cancel window.

## Budget and stops

Planned lines by part (the probit half's share in parentheses). The overrun is what the three engine
slices landed since the plan was drawn showed against their own plans: the built step 0.9 times (its plan
file left out), the k-internal slice 1.7, slice B of response-scale-rows 1.9; their tests alone 1.5, 2.0
and 2.6 times. The tests/cpp line is raised for that: the built step's block, four tests on one link,
landed at 570 lines with its peer accessors, and this one is nine tests on two links sharing that
block's helpers.

| part | planned | at 1.5x to 2x |
|---|---|---|
| engine (chain.hpp, model.hpp) | 205 (45) | 310 to 410 |
| tests/cpp (test_ensemble.cpp, test_peer.hpp) | 520 (110) | 780 to 1040 |
| tinytest | 140 (15) | 210 to 280 |
| exact gate script | 110 (0) | 165 to 220 |
| mixing script | 70 (0) | 105 to 140 |
| SBC harness, sbc.yaml, README line | 70 (0) | 105 to 140 |
| equivalence scenario and the baseline's pins | 35 (0) | 50 to 70 |
| design note, architecture, NEWS, README sections, TODO, indexes | 255 (25) | 380 to 510 |
| records (MANIFEST row, ledger, landing note) | 40 (0) | 60 to 80 |
| total | ~1445 (~195) | ~2170 to ~2890 |

A switch (open call 3) adds about 120 planned lines; a large-n option for the probit half (open call 2)
about 40 to 45. The re-recorded baseline is data. Compute: two
censuses (30 to 40 minutes each), two truth runs (about 5 and 7 minutes), the gate in full (about a
quarter hour for logistic on the laptop, by the timing under [The exact gate](#the-exact-gate)), a local
logistic-k SBC smoke at R = 100, the pre-registered n = 500 cells (about 1.5 hours
of two cores), the large cells (about 4.5 hours of two cores), and the x86 cells (cost, the large cells'
timing and the gate bound, about three and a half hours of the box).

Stop and report, without working around, when (1 to 6 and 8 to 10 pre-registered here):
1. The exact gate fails with the move (any |z| > 4.5, or a pure arm's error above its full bound) after one
   fix.
2. The census misses its bar on either link. If a counter in place of the coin would be the fix, that is a
   state-format change (a [`ChainStateData`](../../src/bartcore/combiner.hpp) field beside
   [`dartNumUpdatesSkipped`](../../src/bartcore/combiner.hpp), the bridge's two state parse sites, the
   sampler's state copy, a reset rule; about 80 lines), outside this plan: it returns to the maintainer.
3. Start-from-truth passes |z| 3 at R = 2000 on either link.
4. On the x86 box any probit cell's ratio passes 1.20 or any logistic cell's
   passes 1.25.
5. An equivalence scenario other than the expected movers moves, or a snapshot moves.
6. The diff passes 2890 lines, or one part passes twice its planned line.
7. After landing: the confirmation matrix flags k on probit-k or logistic-k (an end bin above 45 or a
   chi-square p below 0.01); it is investigated before anything else k-related lands.
8. The logistic gate's with-move and without-move errors are not separated on both machines (early
   measurement 3): the gate would then carry no per-push check that the logistic move runs, which is the
   maintainer's to accept.
9. With the probit half, the re-recorded probit scenarios do not reproduce bit for bit on CI's macOS
   runner (cpp-tests) while the rest of the file does: the runner's `std::erfc` differs from the
   laptop's. The remedies are the evaluator through `Rf_pnorm5` (R's own code, the same bits everywhere,
   about half again the slice's passes) or an erfc carried in the package; either changes the cost the
   probit call is ruled on, so it returns with open call 2. It can first be read when the slice's commits
   run on the runner. The logistic half is not exposed: its evaluator calls only exp and log1p, which the
   recorded draw paths already call.
10. An expectation of the suite, or of a consumer's, that pins no seeded value fails with the move and
    passes without it.

## Order of work

Steps 1 to 6 wait on no ruling: the logistic half is ruled, and the probit half is built as a commit the
maintainer's call can drop. Steps 7 on are as before; review does not start until step 7's results are in.

1. Branch off origin/bartcore at d5b637ee, the tip that carries slice B of response-scale-rows, or at a
   later green engine tip; a private library, `R_LIBS` naming it on every R call.
2. Extract the leaf, k and cache helper from `drawForestRescale` and split the conditions the two
   predicates share out of `forestRescaleApplies`; tests/cpp; on a `--preclean` reference install the
   partition against equivalence-308c45be 55 of 55 identical, bcf-equivalence 15 of 15, multinomial
   11 of 11; commit (neutral).
3. The logistic half: constants, toggle, predicate, coin, move, `LogisticResponse`'s evaluator; tests/cpp
   for logistic, plain and sanitized; the partition still 55 of 55; commit.
4. The probit half, as its own commit: `ProbitResponse`'s evaluator in the erfc form, the family in the
   predicate, tests/cpp for probit and 8; the partition 51 of 55 with the four movers of
   [Equivalence](#equivalence). It stays on the branch until open call 2 is ruled, and is
   dropped before review if the maintainer rules it out. Nothing after it in steps 5 and 6 may depend on
   it: the tinytest file's probit cases and the gate's probit arms pass with and without this commit.
5. tinytest; `--preclean`; full suite; commit.
6. The exact gate's logistic arms, timed; commit.
7. Early measurements 1 to 4, on the step-6 engine and its flipped copy. Results to the maintainer with
   open calls 2 and 3; review waits for them, and the box's availability sets when this step ends.
8. The mixing script's fix and link argument; census and truth on both links; commit (outputs stay in
   scratch).
9. The logistic-k arm, sbc.yaml, the README line; a local R = 100 smoke; commit.
10. After the rulings: the switch if asked for (open call 3); the probit commit dropped, or a large-n
    option added to it (open call 2), its tests/cpp case (a decline without a draw above the threshold, or
    the coin's share at a large fixture) with it.
11. `logistick`; the reference build with `--preclean`: the partition, the re-record, its MANIFEST row and
    the baseline's name where it is pinned; the four snapshot files; commit.
12. The design note, architecture, NEWS, README sections, TODO; commit.
13. Review with mutants; consumers' suites.

## Queue slot and shared files

Everything dec-A191 queued ahead of this slice has landed: gp-copy-continuation, the per-gap cut rule,
the k-internal slice, both parts of setstate-force-update, and slice B of
[response-scale-rows.md](response-scale-rows.md) at d5b637ee. The move is therefore written once against
k's final meaning, which was the reason for the order: the k-internal slice left the built step unchanged
and keeps a binary family's k against the fixed 3 or pi sqrt(3), and an sd-spelled binary prior is
chi(df, `k.scale` / s) in the engine, which the move reads like any other chi prior.

Engine slices stay serial. A slice that lands between this plan's recheck and the slice's landing is
rebased onto; what this slice touches, so that the overlap can be read off the other's diff:
- chain.hpp: `runSweeps`' step block between the level step and the tree loop, `drawForestRescale`,
  `forestRescaleApplies`, the constants beside `rescaleSliceWidth`, one field of `SamplerOptions`.
- model.hpp: one virtual on [`ResponseModel`](../../src/bartcore/model.hpp) and its overrides in
  `LogisticResponse` and `ProbitResponse`. The response setters, the scale and the rows-in members that
  slice B of response-scale-rows added are not touched.
- tests/cpp: test_ensemble.cpp's rescale block and test_peer.hpp's accessors beside it.
- benchmarks: probit-k-scale-exact.R, probit-k-mixing.R, sbc.R, equivalence.R (one scenario appended),
  the baseline and the MANIFEST (this slice merges onto whatever baseline is current at its landing and
  re-records last), the two workflows and mutation-battery.R that pin the baseline's name, sbc.yaml.
- inst/tinytest: one new file, two arms in test-composition-sequences.R, a header line in
  test-probit-rescale-forest.R.
- Records: inst/NEWS.Rd, TODO, docs/architecture.md, the feature matrix, docs/plans/INDEX.md and
  docs/design/INDEX.md.
- A switch (open call 3) would add R/dbarts.R, R/A_class.R, the bridge's control parse and man pages.

Conflicts are expected in the records files and the baseline only. After this slice: the binary
hyperprior rerun of dec-B390 (binary-k-default-revisit) runs on chains carrying the move.

## Settled in planning

1. **A coin, not a counter.** The survey took the move on every fourth sweep by a counter. A counter is
   chain state: to keep a restored chain on the uninterrupted one's draws it must ride the stored state,
   a state-format change this plan otherwise avoids (stop 2 says what it takes). A coin of probability 1/4
   has the same average cost, rides the generator's state, and keeps the posterior trivially. Its cost is
   mixing where the move dominates, by the ratio derived under the mixing gate (an estimated 1.3 or less
   at the census's reading, at most 1.75). If the coin misses the census bar, the counter does not follow as
   overrun: it returns to the maintainer.
2. **Placement ahead of the trees, after the built step.** Placing the move just before the sweep's own
   latent refresh would save its separate redraw (under logistic a PG pass, possibly the larger part of
   the move's cost there), but every recorded channel would then be written before it and have to be
   rescaled in place mid-sweep (the test fits' accumulator and the saved trees' leaves), breaking the built
   step's "every recorded channel is written after it". Not taken; a door if the x86 cost presses.
   Taking the collapsed move first and the built step after would let the built step act on fresh latents
   on the sweeps that take both, rather than have its latent scaling discarded; it departs from the order
   the survey measured, so it is not taken here either.
3. **One slice step**, through the built `sliceFromZero` at the survey's width 0.5. Neither exact nor
   adaptive rejection is available (above). A two-evaluation Metropolis step would cost fewer passes but
   was not measured on the census and needs a step size that cannot follow the posterior's width (0.1 to
   0.9 across n) without depending on the state.
4. **A logistic equivalence scenario** (`logistick`): without it no baseline anchors the logistic move.
5. **The evaluator in the response model**, so the chain stays family-agnostic and the move's likelihood
   is the family's own reported one (for probit, through the faster erfc form, within 1e-14 of it).
6. **The probit-k arm's thin is not changed here**; the census records what it could be (dec-A190).

## Open calls

1. **The logistic half: RULED.** dec-B404 builds the move before 1.0-0 "if a drawn k is kept"; dec-B420
   keeps a drawn k under chi(1.5, 2) for both links ("Keep k modeled and its current prior."), so this half
   is built. Its schedule at large n is not part of either ruling: if early measurement 4 finds the
   logistic move buying nothing there, the option of open call 2 is put to the maintainer for it as well.
2. **The probit half.** The recommendation follows early measurement 4's reading at large n, with the
   x86 cost: include it as planned if every large cell pays; include it with option (a) below if the
   large cells buy nothing while the small ones gained; with option (b) if near-separated large cells pay
   and moderate ones do not. Leave it out only if it buys nothing at the census size either, which the
   survey's numbers make unlikely.
   - Including it: the census's remaining probit tail falls about five times (tau(k) 99th 879 to 160,
     worst 959 to 166; the fitted values' 99th 1003 to 112); on the ten slowest datasets the fitted values'
     worst falls from 1050 to 45. At 5000 rows it gained up to two times or nothing (one dataset a cell);
     10,000 and 50,000 rows are measurement 4. Cost on every default probit binary fit: 12 to 15 percent of
     a sweep on the laptop, with the erfc form the build uses; x86 to confirm. A binary stan4bart fit
     that draws k pays it too, on each of its BART sweeps. About 195 planned lines
     (290 to 390 at the overrun), four re-recorded equivalence scenarios, bench-sampler's binary row past
     its 1.05, probit expectations re-recorded in consumers, and stop 9's exposure to the runner's erfc.
     It would let the probit-k SBC arm drop from
     thin 700 toward about 200, a third of its runner time, by a later call. Whether it converges probit's
     k at n = 500 and 75 trees (split-Rhat 1.12 to 1.21 at 64,000 draws in the k study) is early
     measurement 2.
   - Option (a), the probit move off above n0: the predicate also requires the fit's row count (the
     data's n, which the move does not change, so the decline is exact) to be at most n0, a named engine
     constant with no R setting, n0 the smallest measured size whose large cells buy nothing at both
     signals (10,000, or 20,000 if that cell is run). Above n0 a probit fit draws exactly as with the
     built step alone and pays nothing. About 45 planned lines (engine 10, tests/cpp 20, docs 15), 70 to
     90 at the overrun. Fits either side of n0 take different kernels, both exact; a setData across n0
     switches kernels mid-chain, exact either way. Recommended over (b) when large cells buy nothing at
     both signals, since (b) would still tax them.
   - Option (b), a lower frequency at large n: the coin's probability 1/4 up to n0 and 1/4 n0 / n above
     it, exact because it depends on the data alone; the tax above n0 falls in proportion to n0 / n and
     some of the gain stays. About 40 planned lines (engine 10, tests/cpp 20 for the coin's share at a
     large fixture, docs 10), 60 to 80 at the overrun. Recommended when near-separated large cells pay and
     moderate ones do not.
   - Leaving it out: probit keeps the built step alone, at today's speed and draws; the slow tenth of
     small, near-separated datasets still needs several hundred sweeps per independent draw of k and of
     the fitted values. The logistic half is unaffected: the probit commit carries only the evaluator, the
     family in the predicate and its tests.
3. **A switch.** dec-B397 ruled out a switch for the built step on a measured 0 to 3 percent ("If it's
   basically always beneficial we can consider making it the only option."); this move costs more, so the
   ruling does not cover it.
   - No switch (recommended, unless measurement 4 finds large cells buying nothing and neither large-n
     option of open call 2 is taken):
     one kernel to test and document; every fit gets the faster-mixing k; draws change throughout 1.0
     anyway; a shipped switch is hard to remove. Cost: every eligible fit pays the move, logistic an
     unmeasured share, probit the survey's 12 to 15 percent, for gains at large n that measurement 4
     reads; with option (a) or (b) of open call 2 the large-n tax is gone or cut without a switch.
   - A logical `dbartsControl` switch, default on, fixed at creation: lets a large-n user drop the cost
     where the move gains little; about 120 more lines (R, bridge, help, surface tests).
4. **Negative binomial.** dec-B404's option described the move as "covering negative binomial too".
   Recommended: not in this slice, for the reasons under [Scope](#scope); TODO's k-mixing-pg-families item
   is narrowed to it, to be measured as dec-B404's option said (measure first, stop if it does not help).
   What changed since the plan was drawn: dec-B412's prior study ran and the default stays k drawn under
   chi(1.5, 2) (dec-B431), so the k the move would act on is today's; the maintainer judged that study
   flawed and the default improvable, and a redesigned study runs before 1.0-0, after the merge to main
   (dec-B432), so the default may still move.
   - Not in this slice, measured after dec-B432's study (recommended): the measurement is then made once,
     on the default that ships. It falls after the merge, on an engine that already carries the move for
     the binary links; if it helps, its arm is a slice of its own before 1.0-0.
   - Measured now, at step 7, on the default as it stands: an evaluator for
     [`NBResponse`](../../src/bartcore/model.hpp) and the family in the predicate on a scratch build, not
     committed, and a census with and without it. It would say before the merge whether the move helps
     negative binomial's k at all. Its compute is not known (the latent draw's cost grows with the counts,
     and the default prior draws datasets whose counts run very large), and it is repeated if dec-B432's
     study changes the default.
   - Including it now: about 250 more planned lines (the evaluator, an exact negative binomial frozen-tree
     posterior for the gate, a census arm, the nbinom equivalence scenario re-recorded) and an unmeasured
     cost whose latent redraw grows with the counts.

## Recheck against bartcore

On 2026-10-11, against d5b637ee (bartcore with slice B of response-scale-rows): the notes of what landed
since the plan was drawn were read (the k-internal slice's
[Landing note](k-internal-parameter.md#landing-note), both parts of
[setstate-force-update.md](setstate-force-update.md), the state's record of missable columns, the sparse
starting sigma, gp-copy-continuation, dec-A202, [response-scale-rows.md](response-scale-rows.md) and the
re-recorded baselines), then every symbol, file, test, gate, count and baseline the plan names was read
at that commit, and the runs below were made on a shipped build of it. No rule changed and no call
opened; the text above is written to that code throughout.

What did not change:
- The derivation, the exactness argument, the coin, the placement, the slice step and every rule under
  [Settled in planning](#settled-in-planning).
- The code the move is built from. `drawForestRescale`, `forestRescaleApplies`, `sliceFromZero`, the
  cache bound and its rebuild, the test-only toggle and the built step's tests/cpp block are as the plan
  read them; the k-internal slice left the built step alone, and the forest's `updateK`, `k` and k prior
  keep their names and meanings for a binary family. The latent refreshes, the per-row likelihoods and
  the cancel window are unchanged.
- The stored state: it carries the leaves, a drawn k, the latents and the generator, which is all the
  move needs. No slice since added or removed anything the move reads.
- The equivalence movers (maskprobit, hazard, hurdle and bart2probit with the probit half; none without
  it), the snapshot files, the SBC matrix (M 84, the probit-k arm at thin 700), the mixing script (one
  sampler still serves every dataset), the exact gate script and the study script's settings.

What changed, and where the plan now says so:
- The base and the queue: every slice ahead has landed, so [Order of work](#order-of-work) names its base
  and [Queue slot and shared files](#queue-slot-and-shared-files) lists what this slice touches in place
  of what it waited on.
- Baselines and counts: equivalence-308c45be with 55 scenarios (56 with `logistick`), 42 exact-gate
  commands, the baseline's name pinned in two workflows, the mutation battery and the feature matrix.
- Rows in the likelihood: [The move](#the-move) says which rows the sum runs over and why, for these two
  families, that is the mask alone; tests/cpp 1 pins it, and the tinytest file gains the every-row-out
  and the 0/1-weight cases.
- A restore: it puts the coin on the same sweeps and agrees with the uninterrupted chain to rounding,
  not bit for bit ([State, chains and threads](#state-chains-and-threads)); the continuation arms assert
  what the file's own arms assert, and run long enough for the coin to have fallen.
- The logistic evaluator agrees with the reported log-likelihood to rounding, not bit for bit, the
  compiler fusing a multiply into an add ([Changes](#changes)); tests/cpp 1 compares within 1e-14.
- `std::erfc` on the draw path is a new exposure of the bitwise compare on CI's runner, for the probit
  half only: stop 9.
- Scope: `xbart` with k left at its default, a probit fit under 0/1 weights, `hazard.logistic`'s subject
  weights, an sd-spelled prior and a binary stan4bart fit are placed; what is out says why.
- Negative binomial: dec-B431 and dec-B432 are in [Scope](#scope) and open call 4, which gains the option
  of measuring now.
- The gate's logistic arms are timed on the base, and the budget is recomputed from the three engine
  slices landed since: ~1445 planned lines, stop at 2890.
- Two declines the text lacked: a k of exactly 0 before the coin, a log p(0) that is not finite after it.

Runs, all on a shipped build of d5b637ee: a one-tree fit's `k.scale` under each link; a logistic
sampler's stored state; new counts, a mask, a mask with every row out and an offset on a logistic sampler;
0/1 weights on a probit, a logistic and a `hazard.logistic` fit; `bartBT` with a drawn k, `bart` and
`xbart` under logistic, `hazard.logistic`, an sd-spelled logistic prior; the degenerate chi priors R
accepts; a restore, a copy and a reload of drawn-k probit and logistic samplers against the uninterrupted
chain at three lengths; a probit fit under 0/1 weights against the same mask set by hand; the gate's
fixture timed under probit, logistic and logistic with counts; R's compiler flags on a multiply-add.

Not verified: nothing of the slice is built, so no test, mutant, gate or measurement of the move was run,
and every cost and mixing figure is the survey's, from its prototype on the engine of 2026-10-09. The
prototype was not applied to today's engine. Whether the runner's erfc matches the laptop's is not known.
How many seeded expectations in the suite and in the consumers sit on a drawn-k logistic fit was not
counted. stan4bart was read for its family and its `bart_args`, not built or run; bartCause, treatSens
and bairrtt were not read again. No reference build, sanitizer, x86 run or `R CMD check` was made. The
slices in flight beside this one on 2026-10-11 were not read for overlap.
