# binary-kforest-prior-default

Status: ARC COMPLETE 2026-08-15 (S0 4dbf2dbc, S1 0faeb416, S2 e623fbf3;
  landing notes at EOF). Specced 2026-08-15, two blind critique rounds discharged
  (0 BLOCKER, 6 MAJOR, 13 MINOR, all adopted; see "Critique adjudication").
  Owner forks 1 and 2 SETTLED by VD 2026-08-15 (see "Decisions"). Design
  evidence is GIT-IGNORED session material under
  `.claude/binary-kforest-prior/` (the verified census, and both critique
  rounds), cited the way M4's plan cites `.claude/m4-basis-design/`; every
  load-bearing number is carried HERE rather than pointed at, since this file
  is the only tracked record.
agent: S0 opus (a mandated oracle; the fixture is the whole deliverable and a
  unit-valued factor vacates it). S1 opus (engine + flat C + bridge, three
  layers of one struct). S2 opus (a shipped prior default; the R change is four
  lines and the argument is the slice). Serialized: one implementer, each slice
  lands before the next starts.
rng: S0 DRAW-NEUTRAL (tests only, no `src/` delta). S1 DRAW-NEUTRAL by
  construction - it adds read paths and stores a value the constructor already
  computes and discards; the leaf-scale expression is untouched, the path takes
  no rng, and the trio must come back bitwise on all three baselines. Any
  divergence at S1 is a LEAK, not a re-record: ABORT. S2 DRAW-LAW-CHANGING on
  exactly two named configuration classes - every latent-family K-forest fit
  carrying a basis-free forest, and every K >= 3 fit in every family - and
  BITWISE everywhere else. **No baseline re-records**, because
  `benchmarks/R/bcf-equivalence.R` drives gaussian K = 2 exclusively and both
  moved defaults are identity there (`sqrt(2/2) = 1`, gaussian keeps 2).
window: none. Everything is pre-release and the sister packages update in
  lockstep. S1 moves `sizeof(dbarts_forest_calibration)`, so a `LinkingTo`
  consumer recompiles; that is a migration cost to enumerate, not a design
  input.
budget: S0 ~120-190 dense-equivalent tests (stop 260). S1 engine 20-34 (stop
  48), bridge 10-20 (stop 30), flat C 12-24 (stop 36), R 2-8, docs 50-85, tests
  100-160 (stop 215) - revised up twice from the first draft's four-field
  shape: for the fifth field, the `test-calibration-midchain.R` pin and the new
  `consumer.c` leg, and again for Fork 3b's three truthfulness rules at the two
  state-install sites and falsifier F8's four arms. S2 R 16-32 (stop 46; Fork 6c's two-part K = 1 fix joins),
  docs 60-110 (stop 150), tests 75-135 (stop 185). Dense-equivalent throughout,
  per the arc's standing convention (`multiplier-combiner.md`, "Budget units");
  never report raw against these. Recorded deviation from the ticket:
  `TODO` prices the observability leg at "R 8-14 dense, bridge 0-4", which
  priced a shape (Fork 3a) that does not survive - see Discrepancies 7.

## Goal

The shipped prior default for BINARY (probit and logistic) K-forest models
becomes one a defaulting user can defend, and one they can READ. Three changes
and one pin:

1. the basis-free forest's half-Cauchy median (`aPriorScale`, `forest(sd =)`
   without a basis) becomes FAMILY-AWARE - 2 under gaussian, 1 under probit and
   logistic;
2. the basis forest's node-scale factor default becomes K-AWARE - `sqrt(2/K)`
   rather than the literal 1 - which HOLDS the fixed-variance channel's induced
   index prior at its two-forest value in the ALL-BASIS shape and CAPS it there
   in every other shape, where today it grows without bound in K (leg (b)
   carries both statements and their numbers; the difference is load-bearing
   and no sentence in this document, in `man/forest.Rd` or in `inst/NEWS.Rd`
   may state the stronger form);
3. the five map quantities that are in force but unreadable (`v_f`, the
   half-Cauchy amplitude scale, the node-scale factor and divisor, and the
   basis row norm `c_f`) join the calibration reader at all three layers;
4. the `nodeScaleFactor * anchor / (divisor * rowNorm)` PRODUCT gets its first
   pin, with every factor discriminating, and `tests/cpp` gets its first
   calibration-map coverage of any kind.

The anchor itself does not move. M4.4 settled Option L (`latentScaleAnchor`,
`src/bartcore/chain.hpp`: probit 1.0, logistic `pi/sqrt(3)`), and nothing here
re-opens it - every number below is stated AT that anchor.

## The framing constraint, stated first

The mixing hypothesis this slice was ticketed on is DEAD. M4.4's arm E
(`.claude/m4-basis-design/arm-e-2026-08-14.md`, GIT-IGNORED session evidence,
so its numbers are carried in the M4.4 landing note at the EOF of
`docs/plans/multiforest-extension-surface.md` and repeated here) measured the
native K-forest probit sampler against an R-composed arm B and an independent
reference: AGREES at max |z| 1.46 and 1.73 on 12 functionals, power
precondition met at |z| = 772.49, IACT at or below arm B's on 8 of 12
functionals at mean ratio 0.98. The pinned sigma costs nothing measurable.

So the entire case here is PRIOR COVERAGE: what prior the shipped default puts
on the latent index, what a dbarts user reasonably expects it to be, and
whether the gap is defensible. Arm E is not evidence for this slice in either
direction - it measured a posterior under a DGP the prior easily covers, and
its own LIMITS section records that `update.amplitude = FALSE` on both forests
in every arm, so the half-Cauchy channel - leg (a)'s entire subject - was
pinned OFF throughout.

## The comparator: binding choice

Two comparators are in circulation and the TODO quotes both in one paragraph.
This design uses each where labelled and never mixes them:

- **PRIMARY, for coverage statements: the SHIPPED single-forest binary
  default.** `resolveNodeHyperprior` (`R/model.R`) defaults a binary response to
  `chi(1.5, 2.0)`, so the index prior is `3 s / k` with `k = 2 chi(1.5)`
  (analytic median 1.906131, mean 2.092099, 95% (0.229024, 5.010313); the
  simulation reproduces them at 1.9049 / 2.092 / (0.229, 5.013)). Its induced
  prior on the index: robust sd 1.5847 s, `P(p < 0.01 or p > 0.99)` = 0.2387,
  `P(|eta| > 6)` = 0.0673. This is what a user actually gets from `dbarts()` or
  `bart2()` on a binary y today, so it is what "reasonably expect" means. It is
  also the denominator M4.4's own documentation mandate names (0.238).
- **SECONDARY, for analytic budget statements: the fixed `k = 2` arm.** Index
  prior sd exactly `1.5 s`, tail mass 0.1213, `P(|eta| > 6)` = 0.0001. Used
  wherever a closed form is wanted, because the K-forest map PINS `k = 1` and
  forbids the hyperprior (`resolveNodeHyperprior`'s `multiForest` arm; the
  bridge backstops it in `refuseUnsupportedBCFComposition`, and
  `Chain::buildBCFForest` then sets `forest.k = 1.0` regardless), so the
  K-forest prior is normal-tailed and has no hyperprior analogue. This is the
  denominator behind the TODO's "1.978x" and behind the "0.989x" identity.

The two models are not nested and the ratio changes by a factor of ~2 between
them. Every ratio below names its denominator.

## Problem statement, on prior coverage

### What is in force

Per forest, with `s = latentScaleAnchor(family)`, `m_f` trees, `c_f` the median
nonzero basis row norm, `sd_f` = `forest(sd =)`, `v_f` =
`amplitude.prior.variance` (`Chain`'s BCF constructor, `src/bartcore/chain.hpp`;
`struct ForestSpec`, `src/bartcore/combiner.hpp`):

    leaf.scale_f = sd_f * s / (0.674 * c_f * sqrt(m_f))          [k pinned to 1]
    nodeScale_f  = sd_f * s / (0.674 * c_f)   = prior.scale_f reported by
                                                Chain::forestCalibration

A basis forest contributes `m_f(i) g_f(x_i)` with `m_f(i) ~ N(0, v_f)` at unit
row norm, so its index sd is `nodeScale_f sqrt(v_f)` = `1.049120 sd_f s` at the
defaults. The forest total `g_f` is EXACTLY normal at any x - `m_f` independent
leaves each `N(0, (nodeScale_f/sqrt(m_f))^2)` - so every SD statement here is
exact rather than approximate. The DISTRIBUTION of a contribution is not: it is
a product of two independent normals, which is leptokurtic, and that is why
every tail figure below is simulated rather than read off a normal (a normal at
sd 1.48368 would put 0.1169 outside `(0.01, 0.99)`; the product puts 0.1089).
A basis-FREE forest contributes `a g_0(x_i)` with
`a ~ Cauchy(0, aPriorScale)` (`ForestAmplitudePrior` plus
`BCFForestCombiner::drawShippedGlue`'s `IG(1/2, scale^2/2)` mixture, both
`src/bartcore/combiner.hpp`) and has NO finite sd.

Defaults, all three spellings verified by symbol:
`BCFSpec::aPriorScale = 2.0` (`combiner.hpp`, read only through
`expandForestSpecs` when `spec.forests` is empty);
`forestParams`' `if (withBasis) declared(spec$sd, 1) else 1` at position 4 and
`if (withBasis) 0 else declared(spec$sd, 2)` at position 7 (`R/model.R`);
`bartcoreBCFSampler(sd.control = 2, sd.moderate = 1, ...)` (`R/bartcore.R`).
Two further spellings the census does not enumerate and which a K law must NOT
touch: `BCFSpec::bPriorVariance = 0.5` and `BCFSpec::sdModerate = 1.0`, both the
K = 2 values of the same quantities.

### Leg (a). The half-Cauchy median, under a pinned sigma

VERDICT: **CHANGE, under latent families only.**

Re-derived (seed 20260814, N = 2e6; see Appendix), shipped K = 2 shape
(forest 1 basis-free, forest 2 fixed-variance), probit, at the Option L anchor:

| arm | med abs eta | robust sd | eta 2.5% | eta 97.5% | P(p<0.01 or >0.99) | P(abs eta > 6) |
|---|---|---|---|---|---|---|
| K = 2 shipped, `aPriorScale = 2` | 1.5206 | 2.2544 | -20.29 | 20.35 | **0.3764** | **0.1643** |
| PRIMARY comparator, 1 forest, `k ~ chi(1.5,2)` | 1.0689 | 1.5847 | -7.36 | 7.39 | **0.2387** | **0.0673** |
| SECONDARY comparator, 1 forest, `k = 2` | 1.0118 | 1.5001 | -2.94 | 2.94 | 0.1213 | 0.0001 |

Logistic is the same story at its own anchor: shipped 0.3527 against the
primary comparator's 0.2155.

The honest reading: the K = 2 default puts 58 percent more prior mass outside
`(0.01, 0.99)` than the prior dbarts ships for single-forest binary BART, and
2.4x as much mass within 1e-9 of a boundary. The 2.5/97.5 index quantiles are
+/- 20 where the comparator's are +/- 7.4. The source is the PROGNOSTIC channel
alone: `a mu` has median 1.1749 and 95th percentile 20.26, re-derived.

Why this is a defect and not a shrinkage prior doing its job. The Cauchy
amplitude is a legitimate adaptive-magnitude prior and its heavy tail is the
point; the comparator is heavy-tailed too (a small drawn `k`). What is wrong is
the SCALE, and specifically the transplanted UNIT. Under gaussian, `s` is
`scaledResponseSd()` - the sd of the range-scaled RESPONSE, which already
contains the signal - and sigma is drawn, so a prognostic total at a median
2 sd(y) is both a statement about total variation and correctable by the
sampler. Under probit and logistic, `s` is the LINK'S OWN ERROR sd, sigma is
pinned (`Chain`'s `sigmaIsFixed_` under a non-gaussian, non-aft family), and
`2 s` is a statement that the median prognostic signal is twice the noise -
before any moderator forest is added. Nothing in the sampler absorbs it, and
dbarts's own single-forest binary convention (`defaultNodeScale`: 3.0 for
probit, `pi*sqrt(3)` for logistic - exactly 3 anchor units, divided by a k with
median 1.9) exists precisely to keep the index near (-3, 3).

The candidate values, re-derived on the shipped K = 2 probit shape:

| `aPriorScale` | med abs eta | robust sd | P(p<0.01 or >0.99) | P(abs eta > 6) |
|---|---|---|---|---|
| 2 (today) | 1.5172 | 2.2494 | 0.3756 | 0.1638 |
| 1.5 | 1.2555 | 1.8615 | 0.3175 | 0.1265 |
| **1** | **0.9941** | **1.4739** | **0.2468** | **0.0867** |
| 0.674 | 0.8151 | 1.2083 | 0.1909 | 0.0598 |
| 0.5 | 0.7155 | 1.0607 | 0.1577 | 0.0445 |

**1 reproduces the PRIMARY comparator's coverage to within 3.4 percent on tail
mass** (0.2468 against 0.2387) and is slightly tighter on scale (1.4739 against
1.5847). The same number works in BOTH latent families - logistic at
`aPriorScale = 1` gives tail 0.2255 against the comparator's 0.2155, robust sd
2.6730 against 2.8730 - because `defaultNodeScale` and Option L's anchor are
the same multiple of the latent sd in both.

It also reads correctly: `forest(sd =)`'s default becomes 1 for EVERY forest
under a latent family, one anchor unit, which is the sentence
`man/forest.Rd`'s `sd` item wants to say and today cannot.

### Leg (b). The sqrt(K) dispersion law

VERDICT: **CHANGE, all families, as a DEFAULT and not as map algebra.**

The map has no per-K renormalization: `ForestSpec` carries factor and divisor
per forest, `expandForestSpecs` never reads a count, and the BCF constructor's
calibration loop is a plain per-f loop with no `forestSpecs.size()` in the
expression. So the all-basis index disperses as

    sd(eta) = sqrt(0.5)/0.674 * sqrt(K) * s = 1.049120 sqrt(K) s

exponent exactly 1/2 by construction. Re-derived: `sqrt(0.5)/0.674 =
1.049119853`, giving 1.48368 at K = 2, 2.96736 at K = 8, and
`2.96736/1.5 = 1.97824` against the SECONDARY comparator.

The coverage consequence, re-derived on the all-basis probit shape (the
decisive arm, because `aPriorScale` does not appear in it):

| K | index sd | `P(p<0.01 or >0.99)`, current | same, under `sqrt(2/K)` |
|---|---|---|---|
| 2 | 1.48368 | 0.1089 | 0.1089 (identical) |
| 4 | 2.09824 | 0.2293 | 0.1118 |
| 8 | 2.96736 | 0.4001 | 0.1137 |

The current default's tail mass nearly QUADRUPLES from K = 2 to K = 8 while the
law holds it flat to within 0.005. That is the whole argument: a user who adds
a fourth moderating forest to a three-forest model did not ask for a more
diffuse prior on p, and under a pinned sigma nothing tells them they got one.

**Is sqrt(K) a defect at all?** Not on its own terms - adding an additive
component adds prior variance to the index, and that is what additive priors
do; the per-forest anchor is a statement about ONE forest and the caller owns
the sum. That position is defensible for the KNOB and indefensible for the
DEFAULT: a default is what the caller gets when they have said nothing about
the sum. So the knob keeps its per-forest reading (a declared `sd = 1` means
exactly what it means today, at every K) and only the DEFAULT acquires the
budget reading.

The law: **default `nodeScaleFactor` = `sqrt(2/K)`**, K the forest count.
Re-derived, `sqrt(2/K)` = 1.414214, 1, 0.816497, 0.707107, 0.5, 0.353553 at
K = 1, 2, 3, 4, 8, 16, and the all-basis index sd is 1.48368 s at every one of
them. Its three properties:

- **K = 2 is the fixed point.** `sqrt(2/2) = 1`, so every shipped two-forest
  configuration - bcf, every `tests/cpp` fixture, every `bcf-equivalence`
  scenario, `test-bcf-family.R`'s anchor arm - is BITWISE unchanged.
- **It states the budget in the comparator's own units - IN THE ALL-BASIS
  SHAPE.** There, and only there, the index prior is 1.48368 s at every K,
  which is 0.98912x the SECONDARY comparator's `1.5 s`: the default splits the
  classic binary node-scale budget among the K forests instead of handing each
  forest a full copy of it. The 0.989x identity does NOT hold for the shipped
  one-basis-free shape; see the cap paragraph below, which is the only correct
  statement about that shape and is the one the user-facing sentences take.
- **It is defined at K = 1** (1.414214, same index total), which is what the
  law owed and what `test-bcf-family.R`'s K-written index assertion needs.

**The cap, stated exactly, because it is what the shipped shape actually
gets.** On the shipped shape - one basis-free forest plus K-1 basis forests,
which is bcf's shape, the only one `bartcoreBCFSampler` can build, and the
default whenever forest 1 declares no basis - the fixed-variance channel under
the law is

    1.049120 sqrt(K-1) sqrt(2/K) s = sqrt((K-1)/K) s / 0.674

| K | 2 | 3 | 4 | 8 | limit |
|---|---|---|---|---|---|
| channel sd | 1.0491199 | 1.2114193 | 1.2849042 | 1.3878551 | 1.4836795 |
| x the `k = 2` budget `1.5 s` | 0.699 | 0.808 | 0.857 | 0.925 | 0.989 |

So on the SHIPPED shape the law CAPS rather than pins: the channel rises
monotonically toward 0.989x and never reaches it, against today's unbounded
`1.049120 sqrt(K-1) s`, which is 1.850x the budget at K = 8 and passes 2x at
K = 10 (re-derived from the exact `sqrt(0.5)/0.674`:
`u * sqrt(7) / 1.5 = 1.850473`, `u * 3 / 1.5 = 2.098240`). The reason it caps
rather than pins is
structural: the Cauchy channel has no finite variance to enter a
root-sum-square budget with, so a K-forest shape carrying one cannot have its
total pinned by any per-forest law. That is the same fact that killed index
normalization, and every user-facing sentence this slice writes must state the
cap and not the pin.

### Leg (c). The observability gap

VERDICT: **CHANGE, at all three layers.**

Verified: `bartcore_getCalibration` (`src/R_interface_bartcore.cpp`) emits
exactly seven columns - `prior.scale`, `prior.sd`, `prior.mean`, `k`,
`k.has.hyperprior`, `response.scale`, `response.shift` - plus a `leaf.model`
attribute. `$getForestAmplitudes` (`R/dbarts.R:1450`) returns DRAWS.
`ForestCalibration` (`src/bartcore/chain.hpp`) carries the same seven fields.
`v_f` reaches no reader; `c_f` is computed inside `Chain::basisRowNorm` at
construction and again inside `Chain::setForestBasis`, and STORED NOWHERE.

The census's "cannot reconstruct at all" is now too strong -
`man/forest.Rd`'s `amplitude.prior.variance` item ships the induced-index
formula, the `v_f` default and the `$getCalibration` source of `s_f`, and
`data@bases` is a documented public slot. What survives, and is the actual
defect: `v_f` in force, the factor and divisor in force, and `c_f` in force are
all unreadable, and `c_f` is genuinely unobservable rather than merely
undocumented - after a `$setForestBasis` swap the user would have to
re-implement `Chain::basisRowNorm`'s median-of-nonzero-row-norms rule,
including its R-style even-count averaging, against `data@bases[[f]]`.

This leg is what makes the other two checkable by the person who gets them
wrong, and after leg (b) it is what makes a K-aware default legible: the
`node.scale.factor` column prints `sqrt(2/K)` and the user sees the budget.

### Leg (d). The unpinned product

VERDICT: **CHANGE (tests only), and it lands FIRST.**

Verified exhaustively: every `forest(sd = )` occurrence in `inst/tinytest` is
under a gaussian response (`test-bcf-creation.R:51`, `:103`, `:321`, `:400`,
`:713`, `:836`; `test-forest-basis-r5.R:199`), and every latent-family
K-forest fixture omits `sd` - `test-bcf-family.R`'s `basisSampler` passes no
`forests =` at all, and `test-forest-basis-r5.R`'s `probitForests()` uses bare
`forest()`. So `nodeScaleFactor` is 1 in every anchor assertion the arc has,
and `test-bcf-family.R:109-113` pins `1.0 * s / (0.674 * c)` with the `1.0`
written as a literal. `tests/cpp` has NO calibration-map coverage of any kind:
no `0.674`, no `basisRowNorm`, no `latentScaleAnchor`, no non-unit
`nodeScaleFactor` anywhere in `tests/cpp/*.cpp`; its BCF fixtures set
`forests[f].leaf.scale` directly.

This is the arc's own recorded hazard - "a pin fixture must give every factor
in the pinned expression a DISCRIMINATING value; unit values silently vacate
pins" (`multiplier-combiner.md`, M4.0's required fix and M4.3's Arm 5(iii) dead
pin) - sitting on the exact expression this slice moves the default of.

### What does NOT change, and why

- **`amplitude.prior.variance`, default 0.5.** It is the PARTNER of the 0.674
  divisor, not a free scale: with `b0, b1 ~ N(0, 1/2)` the contrast
  `b1 - b0 ~ N(0, 1)` has half-normal median 0.674, so dividing by 0.674 is
  exactly what places the effect `(b1 - b0) tau` at `sd_f` anchor units.
  Move `v` and that identity - which `docs/design/bcf.md` and
  `multiplier-combiner.md` both state - stops being true. A K law routed
  through `v` would ALSO shrink the reported amplitudes themselves
  (`$getForestAmplitudes`) and the `(b1 - b0)` effect reading. Fork 2 records
  this as the rejected placement.
- **The 0.674 divisor**, for the same reason.
- **The Option L anchor.** Settled at M4.4; the L/C ratio is exactly 3
  regardless of this default, and the decisive all-basis arm does not involve
  `aPriorScale`.
- **The `k = 1` pin and the hyperprior refusal.** Out of scope; they are what
  make the K-forest prior normal-tailed and the comparator not nested.
- **`calibrationMapName`'s "two-forest calibration map"**
  (`src/R_interface_bartcore.cpp`), wrong for K != 2 and reachable through
  `$setCalibration`'s refusal. Ticketed to `bcf-naming-generalization`; this
  slice edits adjacent prose and must not widen into it.

## Decision forks

### Fork 1. The `aPriorScale` default: value, and family-awareness

Alternatives:

- **(1a) Leave at 2 everywhere; document the diffuseness.** Cost: zero code,
  zero re-record. Buys the M4.4 documentation debt discharged and nothing else.
  Rejected: the mandated sentence would then read "the prior we ship is 1.6x
  more diffuse than the one we ship, and here is the knob" - a default nobody
  would choose, shipped into the 1.0-0 window that locks it.
- **(1b) Move to 1 in EVERY family.** Cost: re-records
  `bcf-equivalence-8b047f8b` (12 scenarios), moves the gaussian bcf default
  away from Hahn/Murray/Carvalho's `use_muscale` convention that
  `docs/design/bcf.md:272-274` states dbarts reproduces, and turns
  `test-bcf-creation.R`'s public-vs-internal draw-for-draw oracle red until
  both spellings move together. Rejected: the gaussian case has a DRAWN sigma
  and a response-sd anchor, so the coverage argument does not reach it, and
  paying a baseline re-record plus a fidelity regression for a case with no
  defect is unpriced.
- **(1c) FAMILY-AWARE: 2 under gaussian, 1 under probit and logistic.**
  RECOMMENDED. Cost: a `defaultAmplitudePriorScale(family)` helper in
  `R/model.R` beside `defaultNodeScale`, a `family` argument on `forestParams`
  (`R/spec.R`'s call site already has `family` in scope), and
  `bartcoreBCFSampler`'s `sd.control` defaulting to NULL and resolving after
  the effective family is known (`R/bartcore.R` already resolves
  `family = NULL` off `sampler$model@family`). `bcf-equivalence` is gaussian
  (`y <- mu + z * tau + rnorm(n, sd = 0.2)`, no `family =` on any of its 12
  `bartcoreBCFSampler` calls), so it stays BITWISE and does not re-record.
- **(1d) A value chosen to match the comparator exactly** (solve for tail mass
  0.2387; lands near 0.96). Rejected: a two-significant-figure default derived
  from a Monte Carlo target is unreadable and unquotable, and 1 is already
  within 3.4 percent while carrying a statable meaning ("one anchor unit").

Sub-fork: **does `BCFSpec::aPriorScale = 2.0` follow?** RECOMMENDED: NO, and
say why in its comment. `createBCFSampler` has no flat-C entry point (verified:
`inst/include/dbarts/dbarts.h` carries no BCF creation symbol; the only callers
are `src/R_interface_bartcore.cpp`, which fills `spec.forests`
unconditionally, and `tests/cpp/test_model.cpp`). So that member initializer is
a FIXTURE default reachable by no consumer, and the two-forest spelling it
belongs to is gaussian bcf's, where 2 is correct. Making it family-aware would
put a family branch inside `expandForestSpecs`, which is the one adapter the
`bcf as the K = 2 instance` contract requires to stay a pure re-spelling.

### Fork 2. The K law: whether, where, and on what count

Alternatives:

- **(2a) Do nothing; document the exponent.** Cost: zero. Already half-done -
  `multiplier-combiner.md`'s "Two further facts" paragraph and contract
  sentence (4) carry the exponent and the endpoints. Rejected on the coverage
  table above: 0.1089 -> 0.4001 tail mass from K = 2 to K = 8 is a default
  drifting off the comparator by a factor of 3.3 with no signal to the user.
- **(2b) K-aware `nodeScaleFactor` default, `sqrt(2/K)`.** RECOMMENDED.
  One line in `forestParams`: `if (withBasis) declared(spec$sd, sqrt(2 /
  numForests)) else 1`. The knob's contract is untouched (a declared `sd`
  overrides), the map stays K-free (so the plan's rejection of index
  normalization stands undisturbed), the value is directly readable off
  `prior.scale` today and off `node.scale.factor` after S1, and K = 2 is the
  fixed point so nothing shipped moves.
- **(2c) K-aware `amplitude.prior.variance` default, `1/K`.** Reaches the same
  index total. Rejected: it breaks the 0.674-divisor identity (see "What does
  not change"), it shrinks the reported amplitudes rather than the leaf scales,
  and `v_f` is the one map quantity NO getter reports - routing the law through
  the invisible knob is the opposite of leg (c).
- **(2d) Normalize on the count of BASIS forests `Kb` rather than on K.**
  Rejected with a number: at bcf's shape `Kb = 1`, so any law normalized to 1
  at `Kb = 1` leaves bcf alone but sends the ALL-BASIS K = 2 shape to
  `1.049120 s` - which moves M4.4's item-25 anchor gate, the one arm the whole
  Option L verdict was argued on (0.989x the comparator). Normalizing on K
  leaves both K = 2 shapes exactly as they are.
- **(2e) Apply the law only under latent families. CLOSED 2026-08-15: VD
  declined 2e; the law ships all-families.** The tradeoff prose below stands as
  the record of what was weighed, not as an open question. REJECTED, but on a
  re-argued basis - the first draft's reason ("a family-dependent default for
  `forest(sd =)` cannot be stated in one Rd sentence") is refuted by Fork 1c,
  which is exactly such a default, in the same slice and the same
  `man/forest.Rd` item. After 1c that Rd sentence is family-branched anyway, so
  2e's marginal documentation cost is ONE CLAUSE, not a sentence. The reason
  that survives is a gaussian coverage argument stated in gaussian's own
  units, and it has to be stated because the `P(p)` table is defined only for a
  binary link and does not reach gaussian:

  **The MAGNITUDE argument is withdrawn for gaussian; it does not
  discriminate.** The first revision argued that the all-basis default puts the
  prior sd of `E[y|x]` at `1.049120 sqrt(K)` times `sd(y)` - 1.48 / 2.10 / 2.97
  at K = 2 / 4 / 8, which is right, `s = scaledResponseSd()` being `sd(y)` in
  internal units - and called 2.97 indefensible. Applying the SAME criterion to
  what dbarts ALREADY ships for a single-forest gaussian fit (`node.scale` 0.5,
  `k` 2, so `0.25 / s` in `sd(y)` units) refutes it. Re-derived, 2000
  replicates:

  | y | n = 50 | 200 | 1000 | 5000 |
  |---|---|---|---|---|
  | normal | 1.128 | 1.378 | 1.621 | 1.835 |
  | uniform | 0.836 | 0.858 | 0.865 | 0.866 |
  | **t3** | 1.446 | **2.288** | **3.940** | **6.795** |

  On t3 data at n = 1000 the SHIPPED single-forest gaussian default already
  sits at 3.94 `sd(y)`, above the 2.97 the first revision called indefensible,
  and it grows without bound in n through range-scaling. So no magnitude
  threshold separates the K = 8 default from this package's own conventions,
  and the first revision reached its band only by excluding the one refuting
  arm with the adjective "normal-ish". Withdrawn.

  **What survives is a COMPOSITIONALITY argument, and it is the one the law
  actually embodies.** The single-forest default is one ensemble's whole
  budget, however large; the `sqrt(K)` growth compounds that budget ACROSS
  ensembles, so the prior on `E[y|x]` depends on how the user DECOMPOSED the
  mean function rather than on what they said about it. Writing the same model
  as two forests or as four changes the prior on the total by `sqrt(2)` while
  changing nothing the modeller intended. A default should be invariant to a
  re-expression the modeller regards as cosmetic, and `sqrt(2/K)` is exactly
  that invariance; the magnitude at K = 2 - bcf's own settled convention, from
  Hahn, Murray and Carvalho's prognostic-dominant causal model - is preserved
  rather than judged. This argument is family-free, needs no threshold and no
  data-dependent comparator, and it is why the recommendation survives the
  withdrawal above.

  **Under a binary link the magnitude argument DOES discriminate, and both
  arguments apply.** There `s` is the link's own error sd, fixed by the family
  rather than by the data, so the comparator does not move with `n` or with the
  tails of `y`: 0.4001 tail mass at K = 8 against the shipped single-forest
  binary default's 0.2387, with no arm of the comparator anywhere near it.

  **The migration cost, stated as a cost and not as an argument.** Every
  gaussian fit with K >= 3 draws differently after S2. No equivalence baseline
  and no tinytest snapshot pins those draws (`bcf-equivalence` is K = 2;
  `equivalence` and `multinomial-equivalence` never build a combiner; the
  suite's three K = 3 fixtures assert reconstruction identities and ratios).
  That absence is a MIGRATION fact - it means the change is cheap to land, not
  that it is right to make - and it must never be written as a reason for the
  change. A user wanting the previous model declares `forest(sd = 1)`, which
  falsifier F3 pins bitwise.

  **The fallback, priced, and the tradeoff as the owner should see it.** 2e is
  a live fork, not a settled one. Taking it costs one clause in the `sd` Rd
  item (1c family-branches that sentence anyway), leaves gaussian K >= 3
  bitwise, and leaves the gaussian default's prior on `E[y|x]` dependent on the
  user's decomposition. Declining it moves shipped gaussian K >= 3 draws on a
  compositionality argument rather than on a coverage number, in a family where
  a drawn sigma does correct for the excess in the posterior. The asymmetry
  that decides my recommendation is the release window: pre-1.0-0 a default
  move costs a re-run and nothing else, and after it the surface is locked, so
  the branch that can be revisited cheaply is the one that ships the
  invariance now. An owner who weighs "do not touch a working gaussian
  default" above that reaches 2e, and nothing in this document should read as
  if that were a mistake.

Sub-fork: **does the basis-FREE forest's default also become K-aware?**
RECOMMENDED: NO. A Cauchy channel has no finite variance to enter the budget
with; dividing its median by `sqrt(K)` would be an arbitrary scaling of a
quantity the budget cannot contain. The consequence is stated rather than
hidden: on the shipped shape the induced tail mass still drifts, re-derived at
0.2458 / 0.2734 / 0.2894 for K = 2 / 4 / 8 under the recommended pair, against
0.3757 / 0.4479 / 0.5496 today.

### Fork 3. Closing the observability gap: which surfaces grow what

Alternatives:

- **(3a) R-only, off the control attribute.** `attr(control,
  "bartcore.bcf")$params[[f]][4:6]` already carries factor, divisor and `v_f`
  R-side, so a reader could be written with zero engine work (the costing
  study's 8-14 dense estimate assumed this). Rejected: it creates a SECOND
  source of truth for a quantity the engine re-derives -
  `Chain::setForestBasis` recomputes `c_f` and rewrites the leaf scale from it,
  so an R-side reader would report the CREATION-time row norm after a
  mid-sample swap. `Chain::forestCalibration`'s own docstring calls it "the
  authoritative reader of what is IN FORCE"; a sibling reader that can go stale
  against it is worse than no reader.
- **(3b) FIVE new fields on `ForestCalibration`, surfaced at all three
  layers, every one a SPEC ECHO except the row norm.** RECOMMENDED. Engine:
  `ForestCalibration` gains `amplitudePriorVariance`, `amplitudePriorScale`,
  `nodeScaleFactor`, `nodeScaleDivisor`, `basisRowNorm`; `Chain` gains three
  small vectors written in the BCF constructor straight from `forestSpecs`
  (the two amplitude-prior values) and from the value the constructor and
  `setForestBasis` already compute (the row norm). Bridge: five columns on
  `bartcore_getCalibration`. Flat C: five `double*` fields appended to
  `dbarts_forest_calibration` BELOW the documented "1.0-0 field boundary",
  five `offsetof` static_asserts plus the `sizeof` assert updated, five `FILL`
  lines.

  **A spec echo is a LIE after a state install, and the fix is a truthfulness
  bit, not a smaller promise.** Verified by opening both paths:
  `installForest` (`src/bartcore/chain.hpp:3235`) and the restore path
  (`:3323`) each execute `if (fs.leafScale > 0.0) forest.leaf.scale =
  fs.leafScale`, by design - `ForestStateData::leafScale`'s own doc block says
  a continuation restoring `k` without it "would pair a donor's trees with the
  destination's calibration" - and each then calls `combiner_->restoreGlue`
  (`:3246`, `:3355`), which overwrites `glue_.prior[f].variance` from
  `state.amplitudeVariances[f]` (`src/bartcore/combiner.hpp:1044-1045`).
  Neither `$setState` (`R/dbarts.R:1580`, `refuseHostMutation` only) nor
  `$installTrees` (`R/dbarts.R:1612`) carries a BCF refusal, `glueIsValid`
  checks widths and totals only, and `inst/tinytest/test-forest-basis-r5.R:110-136`
  already drives `setState` on a K-forest sampler. So a pure echo would report
  the RECIPIENT's decomposition beside a `prior.scale` derived from the
  DONOR's - the map identity false inside one returned matrix, which is
  exactly the failure Fork 3a was rejected for. Three rules make every column
  true on every path:

  1. **The amplitude columns follow the state.** Both install sites update
     `amplitudePriorVariances_[f]` from `state.amplitudeVariances[f]`, under
     the same guard `restoreGlue` uses to reach that loop (`state.hasBCF`,
     non-empty `state.amplitudeWidths`, and `combiner_->glueIsValid(state)` -
     an EXISTING base virtual, `combiner.hpp:693`, already called from
     `chain.hpp:2809`, so no new virtual and Fork 3b's whole reason survives).
     `amplitudePriorScale` is not serialized at all, so it needs nothing.
  2. **The two map-decomposition columns carry a truthfulness bit.**
     `Chain` gains `nodeScaleIsMapDerived_`, one flag per forest, set false at
     either install site when `fs.leafScale > 0.0 && fs.leafScale !=
     forest.leaf.scale` - compared BEFORE the assignment, so a self-restore
     (bitwise-identical scale) keeps its columns and only a genuinely foreign
     calibration loses them - and set true again by `$setForestBasis`, which is
     the call that re-imposes the map from the stored factor and divisor
     (`chain.hpp:934-945`). `node.scale.factor` and `node.scale.divisor` report
     NaN while the flag is false. The stored values are NOT cleared: they are
     what `setForestBasis` re-derives from.
  3. **`basis.row.norm` needs no rule.** Bases are not state, and
     `setForestBasis` already rewrites it.

  The result is that the identity is available exactly when it is true, which
  is the same NaN discipline the off-map samplers already get, extended to the
  one on-map case where the decomposition becomes genuinely unknown. Rejected
  alternatives, each on a stated ground: NARROWING THE Rd CONTRACT ("the
  identity does not survive a restore") ships a reader whose numbers are wrong
  on a live sampler and asks the user to remember when - a documented lie is
  what Fork 3a was rejected for; RE-DERIVING a factor from the installed leaf
  scale (`factor = leafScale * sqrt(m) * divisor * rowNorm / anchor`) keeps the
  identity but INVENTS a decomposition, since the state carries only the scale
  and the donor's own divisor and row norm are unknowable; REFUSING `$setState`
  and `$installTrees` on K-forest samplers deletes shipped, tested and
  deliberately designed behaviour (the restore-then-widen semantics M4.3 built,
  `multiplier-combiner.md`, "Persistence") to make a getter tidy.

  **The engine reads its own SPEC, never the combiner.** The first draft had
  `v_f` come from `glue_.prior[f].variance`, which would need a new virtual on
  the `ForestCombiner<L, ResidT>` base (`Chain` holds it by `unique_ptr`) -
  and that perturbs the vtable of every combiner, which can shift inlining and
  therefore FMA formation, in a codebase that already documents
  contraction-dependent bitwise divergence (`BCFSpec::generalAmplitudeDraw`'s
  comment) and sets no `-ffp-contract` in `src/Makevars.in`. Echoing the spec
  instead adds no virtual, no combiner member and no call into the combiner at
  all, which is what makes S1's draw-neutrality argument structural rather
  than empirical. The cost is that the LIVE scale-mixture auxiliary stays
  unreadable - correctly, because it is STATE (it serializes as `aVariance`)
  and not calibration; a state reader is a separate question and a door. R: the `$getCalibration` docstring and two Rd sections.
- **(3c) (3b) plus an `$getIndexPrior()` reader returning `sd(eta_i)`
  quantiles.** Rejected for this slice, recorded as a door. With (3b) every
  INPUT to the formula is readable and the formula itself already ships in
  `man/forest.Rd`; a derived reader would have to decide what to report when a
  basis-free forest makes the sd undefined, which is a modelling statement
  belonging in prose. Ship a six-line worked example in the Rd instead.

Sub-decisions inside (3b), all recommended:

- **NaN, not a neutral value, off the map.** The five fields are non-NaN
  exactly when the forest has a calibration-map entry (`f <
  nodeScaleFactors_.size()`, sized only by the BCF constructor). A single-forest
  or multinomial sampler reports NaN, which is a positive signal ("this
  sampler's scale is not map-derived") rather than a plausible 1.0 a caller
  would multiply by.
- **The two amplitude columns are EXCLUSIVE, and that is what makes their names
  honest.** `ForestAmplitudePrior` carries `{variance, halfCauchyScale}` and a
  forest is one kind or the other (`forest.ridge = amplitudePriorScale > 0.0`,
  `applyBCFSpec`). So a fixed-variance forest reports
  `amplitude.prior.variance` = its prior variance and
  `amplitude.prior.scale` = NaN, and a scale-mixture forest reports the
  reverse. This is the resolution of a real naming collision: reporting a
  moving auxiliary under the name `amplitude.prior.variance` would use the name
  of an argument `resolveForests` REFUSES on exactly the forest it would be
  reporting for. Each column now reports a prior the caller may set, or NaN.
  It also closes leg (a)'s own observability: after S2 the user can READ that
  their basis-free forest's half-Cauchy median is 1 rather than 2.
- **No `node.scale.anchor` column.** `s` is recoverable exactly as
  `prior.scale * divisor * rowNorm / factor`, and the Rd states that identity
  WITH its condition: it holds whenever `node.scale.factor` is not `NA`, which
  by rule 2 above is exactly when the calibration in force is the map's.
  Recorded as a door, since under gaussian `s = scaledResponseSd()` is
  data-dependent and the derivation is the only route to it.
- **No `DBARTS_C_API_MINOR` bump and no hash re-bake.** `dbarts_apiHash()` is
  FNV-1a over the stringized `DBARTS_C_API_LIST` (function signatures), which
  does not move; the header's own append rule conditions the MINOR bump on
  "after 1.0-0" and we are before it. `test-capi.R:59`'s literal hash pin
  stays. `sizeof(dbarts_forest_calibration)` DOES move, so a `LinkingTo`
  consumer recompiles - the sixth of the recorded consumer gotchas, and the
  reason `DBARTS_HAS_FIELD` exists.

### Fork 4. What the product-pinning test asserts

Alternatives: (i) `prior.scale_f == sd_f * s / (0.674 * c_f)` with all four
factors discriminating; (ii) the same after a `$setForestBasis` swap; (iii) the
induced index at non-default `sd_f`.

RECOMMENDED: **(i) and (ii) both, in `inst/tinytest/test-bcf-family.R` and
`inst/tinytest/test-forest-basis-r5.R` respectively, plus a `tests/cpp` twin;
(iii) is subsumed** - the index assertion is a function of the `prior.scale`
values (i) already pins, so it adds a second reading of the same numbers rather
than a second constraint. Full design in "Test design" below.

The load-bearing constraint the alternatives hide: **the probit arm cannot
discriminate a DROPPED anchor**, because `s = 1` there and
`sd_f * 1 / (0.674 c_f)` equals `sd_f / (0.674 c_f)`. The pin is therefore
meaningless unless it runs under LOGISTIC (`s = 1.813799`) as well. This is the
same shape as M4.4's recorded near-miss (the naive anchor's logistic arm missing
by only 0.865 percent) and it is why the fixture cannot be probit-only.

### Fork 5. M4.4's undischarged documentation mandate

The mandate, verified verbatim at `docs/plans/multiforest-extension-surface.md`
item 10: M4.4 "must NOT move it", and "must DOCUMENT that the shipped binary
K-forest prior is more diffuse than the shipped single-forest one
(`P(p < 0.01 or p > 0.99)` = 0.376 versus 0.238) and point at the deferred
slice."

Verified UNDISCHARGED: `0.376`, `0.238` and any "more diffuse" statement appear
nowhere in `docs/`, `man/`, `inst/NEWS.Rd`, `R/`, `src/` or `vignettes/` outside
the plan's own item 10. (Two unrelated hits: a CATE interval in
`grow-from-root-default.md` and a table cell in `composition-mixing-probe.md`.)

Alternatives: **(5a) discharge it as written** (document today's diffuseness,
point at this slice) - now self-defeating, since this slice removes the thing
the sentence describes; **(5b) discharge it in its SUCCESSOR form** -
RECOMMENDED - the same comparison, at the NEW defaults, written into
`man/forest.Rd` and `multiplier-combiner.md`, plus the landing note recording
that M4.4's obligation was met here and how the numbers moved; **(5c) leave
it** - rejected, it would leave the tree in the same state twice, which is the
census's own phrasing and correct.

The successor sentence carries, in this form and no stronger (M1's correction;
the pin/cap distinction is the part a reader will over-read): at the shipped
defaults a K = 2 binary K-forest prior puts `P(p < 0.01 or p > 0.99)` = 0.2468
(probit) against 0.2387 for the single-forest binary default, where before this
slice it was 0.3764; the fixed-variance channel's induced index prior is HELD
at 1.4837 s - 0.989x the `k = 2` index budget - at every K when every forest
carries a basis, and BOUNDED BY it, rising monotonically from 0.699x, when one
forest does not.

### Fork 6. K = 1 through the K-forest path

**Re-traced by opening both routes, and the two K = 1 routes do DIFFERENT
things - neither of which is the single behaviour the first draft assumed.**

- **The data route.** `dbartsData(x, y, bases = list(b))` - `validateForestBases`
  (`R/data.R:635-668`) imposes no length floor, `R/spec.R`'s
  `if (!is.null(data@bases))` block computes `numForests <- length(data@bases)`
  = 1 with no floor of its own, `forestParams` returns a length-1 list, and
  `applyBCFSpec` (`src/R_interface_bartcore.cpp`) errors with "bcf parameters
  must be a list of at least two per-forest parameter vectors". Internal
  vocabulary, as the census records.
- **The forests route.** `dbarts(x, y, forests = list(forest(basis = ~z)))`
  does NOT reach that error. `forestBasisDeclarations` (`R/model.R:689-692`)
  returns NULL on `length(forests) < 2L`, so `R/dbarts.R:542-561` never sets
  `dataCall$bases`, `data@bases` stays NULL, `R/spec.R`'s bcf block never runs,
  and `resolveForests` returns a length-1 spec whose only surviving effect is
  `firstForest$n.trees`/`interactions`/`blocks` (`R/spec.R:217-226`). **The
  declared basis is SILENTLY DROPPED and an ordinary single-forest model is
  fit.** Same on `dbartsSpec` (`R/spec.R:653-664`, the same helper).

That second finding is not in the census, not in the TODO and not in the
critique that raised this fork; it is a silent-wrong-model defect on a public
route, and it is strictly worse than the error the fork was opened about.

Alternatives: **(6a) out of scope entirely** - rejected, because this slice's
own law has to state a K = 1 value (1.414214), and because one of the two
routes silently fits a different model than the user wrote; **(6b) make K = 1
reachable** - deferred FOR THIS SLICE, on GATE grounds rather than on cost (see
the corrected pricing below): a new shipped configuration needs acceptance
evidence of its own, and it has nameable enabling value (a single
varying-coefficient forest with no prognostic partner is VCBART's shape,
recorded as D4 in `docs/design/model-space-survey.md`), so it earns its own
ticket rather than a refusal that forecloses it; **(6c) RECOMMENDED, in the
corrected two-part form:**

1. **Drop `forestBasisDeclarations`' `length(forests) < 2L` floor to
   `< 1L`**, so a length-1 `forests` declaration carrying a basis reaches
   `data@bases` like every other one. Verified safe: a length-1 list declaring
   NO basis still expands to `list(NULL)`, still trips the
   `any(!vapply(expanded, is.null, ...))` gate, and still falls through to
   `resolveForests`' existing refusal - so `test-bcf-creation.R:713`
   (`forests = list(forest(sd = 1.5))` -> "single-forest 'forests' has none")
   stays green. **Three behaviours move, not one, and all three are
   enumerated** (none is pinned by a fixture - no test anywhere declares a
   length-1 `forests` carrying a basis): (i) the silent drop becomes the
   designed refusal, which is the point; (ii) `forests = list(forest(basis =
   ..., sd = ...))` today answers "'sd' configures the amplitudes ... a
   single-forest 'forests' has none" because `hasBasis` is `logical(0)`, and
   afterwards `hasBasis` is TRUE so `R/model.R:849`'s branch is skipped and the
   new K = 1 refusal answers instead - a message change, not an acceptance
   change; (iii) a length-1 `forests` declaring a basis OVER a data object that
   already carries K >= 2 bases: `R/spec.R:405-410` replaces `data@bases` with
   the declaration, so the count drops to 1 and the new refusal fires on a call
   that named two bases. Today that call fits K = 2 with the declaration
   silently dropped, so both behaviours are wrong; the refusal text must
   therefore be written against the RESOLVED basis count and name where that
   count came from, rather than saying "K = 1" flatly at a user who wrote two
   bases.
2. **One designed refusal in `R/spec.R`'s `!is.null(data@bases)` block**, where
   `numForests <- length(data@bases)` is already computed and where the family
   gate already lives - the site BOTH routes now pass through, verified by
   opening both traces above. It names what a one-forest amplitude model is
   missing (a second ensemble for the amplitudes to distinguish it from) and
   points at the two-forest spelling and at the K = 1 ticket.

Corrected pricing for 6b, since the first draft's deferral cited engine work
that does not exist. Verified live: `BCFForestCombiner::drawGlue`
(`src/bartcore/combiner.hpp`) dispatches on
`forests.size() == 2 && !generalAmplitudeDraw_ && shippedShape()` and sends
everything else to the K-general `drawAmplitudes`, so K = 1 takes the general
path with no new code; `synthesizeIndicatorBasis` is already gated on
`numForests_ > 1`; `shippedShape()` requires `glue_.basis.size() == 2` and is
simply false at K = 1; `expandForestSpecs`' two-forest synthesis is unreachable
from R (`applyBCFSpec` fills `spec.forests` unconditionally); and
`createBCFSampler` (`src/bartcore/facade.hpp`) imposes no forest-count floor at
all. The only floor in the tree is `applyBCFSpec`'s `Rf_xlength < 2`. So 6b is
plausibly a bridge-plus-R change with tests, NOT engine work - which makes the
deferral a scheduling and gate decision, and the ticket must say so rather than
inherit a cost claim that does not survive.

### Fork 7. Does the internal `bartcoreBCFSampler` route follow?

Alternatives: **(7a) yes, family-aware in both routes** - RECOMMENDED; **(7b)
no, pin the value explicitly in the creation oracle and let the routes
diverge** - rejected, since `test-bcf-creation.R`'s public-vs-internal
draw-for-draw oracle exists precisely to hold the two routes on one set of
defaults, and a divergence would make the gate harness (`bcf-equivalence`,
`bcf-exact*`) run a model the public surface cannot express.

Because `bcf-equivalence` is gaussian and gaussian keeps `sd.control = 2`,
(7a) costs no re-record. `sd.moderate` does NOT become K-aware in the internal
route: it is a two-forest spelling and `sqrt(2/2) = 1`.

## Decisions

The two OWNER-level forks are settled. Decided by VD 2026-08-15, on an
in-session walkthrough of the fork list after both critique rounds.

1. **Fork 1: ADOPTED 1c.** The `aPriorScale` default becomes FAMILY-AWARE - 2
   under gaussian, 1 under probit and logistic. The sub-fork lands as
   recommended: `BCFSpec::aPriorScale = 2.0` (`src/bartcore/combiner.hpp`) is
   NOT made family-aware, staying the gaussian-only fixture default that no
   consumer reaches.
2. **Fork 2: ADOPTED 2b AT FULL SCOPE.** The `sqrt(2/K)` default on
   `nodeScaleFactor` applies in ALL families, gaussian included. **2e is
   CLOSED** - declined - so gaussian K >= 3 draws move, on the
   compositionality argument recorded in Fork 2e rather than on a magnitude
   one. The sub-fork lands as recommended: the basis-FREE forest's default is
   not K-aware, a Cauchy channel having no finite variance to enter the budget
   with.

Forks 3 through 7 proceed on this design's own recommendations under the
standing grant, no separate owner call owed: the five-column calibration reader
with Fork 3b's three truthfulness rules, the product pin under both latent
families plus the first `tests/cpp` calibration coverage, M4.4's documentation
debt discharged in its cap-not-pin successor form, Fork 6c's two-part K = 1 fix
with reachability ticketed separately, and the internal `bartcoreBCFSampler`
route following the public one.

## RNG / draw-law classification

Per proposed change:

| change | class | reach |
|---|---|---|
| S0 product pins, `tests/cpp` calibration coverage | DRAW-NEUTRAL | no `src/` delta at all |
| S1 stored `basisRowNorms_` + 4 `ForestCalibration` fields | DRAW-NEUTRAL | read-only additions; no rng in the path |
| S1 bridge columns, flat fields, Rd | DRAW-NEUTRAL | reporting only |
| S2 `aPriorScale` 2 -> 1 under probit/logistic | DRAW-LAW-CHANGING | every latent K-forest fit with a basis-free forest |
| S2 default `nodeScaleFactor` -> `sqrt(2/K)` | DRAW-LAW-CHANGING | every K >= 3 fit, ALL families |
| S2 K = 1 designed refusal | DRAW-NEUTRAL | an error path that today errors elsewhere |

S1's neutrality is argued from MECHANISM, not from the name of the surface (the
M4.4 lesson). Two conditions the implementer must hold, and the reviewer must
check:

1. The leaf-scale expression must keep its exact written shape. Compute
   `double c = basisRowNorm(...); basisRowNorms_[f] = c;` and then pass
   `nodeScaleFactor * s / (nodeScaleDivisor * c)` - the same association, the
   same operand order. Hoisting the subexpression differently, or reusing a
   stored value where the source recomputes, is where a bitwise change would
   enter.
2. `Chain::setForestBasis` must write the stored norm from the SAME call it
   already makes, not from a second call.

Expected verdicts, per baseline, at each slice:

- **S0**: trio NOT RUN (zero engine delta; the `zero-weight-exactness` S0
  precedent).
- **S1**: `equivalence-8b047f8b` 37/37 identical draws at 37 compared / 0
  skipped under `--strict-coverage`; `bcf-equivalence-8b047f8b` 12/12 identical
  on every channel; `multinomial-equivalence-1027be5` 10/10 identical on every
  channel. ANY divergence is a leak: ABORT, do not re-record.
- **S2**: the SAME three verdicts, all still bitwise. This is the slice's
  sharpest gate and its most counter-intuitive claim, so state the mechanism in
  the landing note: `bcf-equivalence` drives `bartcoreBCFSampler` exclusively,
  all 12 scenarios are gaussian and two-forest, gaussian keeps `sd.control = 2`
  and `sqrt(2/2) = 1`, so both moved defaults are the identity there. A
  divergence on `bcf-equivalence` means the family gate or the K gate leaked;
  a divergence on either of the other two means the change escaped the combiner
  entirely. ABORT in both cases. **No baseline is re-recorded anywhere in this
  arc**, and no MANIFEST entry is owed.
- The mandated BCF exact oracles (`bcf-exact.R`, `-restricted`, `-weak`) are
  gaussian K = 2 and therefore unmoved by construction; run `bcf-exact.R quick`
  at S2 as a cheap leak detector, and record that it is a confirmation and not
  a re-validation.

## Shipped-surface deltas (all user-visible)

- `man/forest.Rd`, `sd` item: "The defaults are 2 for a forest with no basis and
  1 for one with" is FALSIFIED by S2 and must state the family-aware and
  K-aware defaults plus the declared-value override.
- `man/forest.Rd`, `amplitude.prior.variance` item: unchanged default (0.5),
  gains the budget sentence and the worked `sd(eta_i)` example.
- `man/dbartsSampler-class.Rd`: `getCalibration`'s value section (five new
  columns, the NaN rule, the exclusivity of the two amplitude columns, the
  anchor-recovery identity) and the forest-index paragraph.
- `inst/tinytest/test-calibration-midchain.R` and `inst/tinytest/capi/consumer.c`:
  not user-facing, listed here because both are MANDATORY same-commit edits at
  S1 and one of them fails S1's own gate if missed (S1 items 5 and 6).
- `inst/NEWS.Rd`: one UPGRADING bullet (results move for the two named
  configuration classes) plus one bullet per slice under the ordinary
  headings.
- `inst/include/dbarts/dbarts.h`: `dbarts_forest_calibration` gains five fields
  below the 1.0-0 boundary; `structSize` moves; the Doxygen for each field.
- `docs/design/multiplier-combiner.md`: the calibration-map section's contract
  sentence (4) and its "Two further facts" paragraph (the sqrt(K) statement
  becomes the budget statement), the "Surfaces" list, and a landing note per
  slice.
- `docs/design/bcf.md`: the two "median 2 sd(y)" statements
  (`:94`, `:272-274`) gain the latent-family qualifier; `sd.control` /
  `sd.moderate` prose likewise.
- `docs/design/nameable-calibration.md`: the calibration-reader rows.
- `docs/design/feature-matrix.md`: per-landing check, and it is NOT vacuous at
  both slices. At S2 a default is not a capability, so no cell moves - record
  the check. At S1 the `nameable calibration` capability's own footnote
  `[f16]` (body from `:368` at tip `ffb9959c`) enumerates what
  `$getCalibration` reports and which tests carry it, so its PROSE moves while
  its cells do not. Do not bump the Status stamp; `feature-matrix-anchor-refresh`
  still owns the full pass, and `5277cd54` re-anchored the file wholesale, so
  every citation into it must be re-derived by opening it.
- `TODO`: the entry closes WITH a recorded deviation - the `sd(eta_i)`
  quantile reader, the fourth deliverable of the entry's observability leg, is
  deferred as a door (Fork 3c) rather than shipped, and the entry must say so
  rather than reading as fully discharged.
  `binary-kforest-k1-reachability` opens.

## Slice decomposition

Order is deliberate: pin the product BEFORE anything moves, make the map
readable BEFORE the default changes, then change the default into a tree where
both hold.

### S0. The product pin and the first `tests/cpp` calibration coverage

No `src/` change. Content is "Test design" below.

Gates: `tests/cpp` from clean (`cd tests/cpp && make && ./test_bartcore`); full
tinytest against a private-library install of the UNCHANGED tree; `air format
--check`. No `--preclean` (no header edit), no trio, no ASAN leg, no R CMD
check (no R/Rd delta).

Budget: tests ~120-190 dense, stop 260. The band is wide because the arc has
under-priced test work against a mandated oracle twice (M4.4 348 against a
155-260 band).

### S1. The observability columns

1. `src/bartcore/chain.hpp`: `ForestCalibration` gains
   `amplitudePriorVariance`, `amplitudePriorScale`, `nodeScaleFactor`,
   `nodeScaleDivisor`, `basisRowNorm`, each defaulting to a quiet NaN; `Chain`
   gains `amplitudePriorVariances_` and `amplitudePriorScales_` (written in the
   BCF constructor straight from `forestSpecs`, NaN on whichever branch the
   forest does not take) and `basisRowNorms_` (written in the constructor and
   in `setForestBasis` from the value each already computes);
   `Chain::forestCalibration` fills the five when
   `f < nodeScaleFactors_.size()` and leaves them NaN otherwise. NO virtual is
   added to `ForestCombiner` and nothing calls into the combiner except the
   EXISTING `glueIsValid` - see Fork 3b for why that is the load-bearing
   choice. Plus Fork 3b's three truthfulness rules, which are engine work at
   the two state-install sites (`chain.hpp:3235` and `:3323`): the
   `nodeScaleIsMapDerived_` flag written there and cleared/restored as
   specified, and the `amplitudePriorVariances_` update from
   `state.amplitudeVariances` under `restoreGlue`'s own guard. Both sites are
   INSTALL paths, not draw paths, and the added work is one comparison, one
   flag write and one copy - no rng, nothing feeding an accumulator, so S1's
   neutrality argument is unchanged.
2. `src/R_interface_bartcore.cpp`: `bartcore_getCalibration` grows five
   columns, named `amplitude.prior.variance`, `amplitude.prior.scale`,
   `node.scale.factor`, `node.scale.divisor`, `basis.row.norm` - 12 in all.
3. `src/C_interface.cpp` and `inst/include/dbarts/dbarts.h`: five fields, five
   `offsetof` asserts, the `sizeof` assert, five `FILL` lines.
4. `R/dbarts.R`: `$getCalibration`'s docstring.
5. **`inst/tinytest/test-calibration-midchain.R:41-58`**: the shipped pins
   `expect_equal(dim(calibration), c(2L, 7L))` and the exact-set
   `expect_identical(colnames(calibration), c(...))` go RED under S1's own
   full-tinytest gate and MUST move in the same commit - to `c(2L, 12L)` and
   the twelve-name vector, keeping the EXACT-SET form rather than relaxing it
   to a subset, since that pin is what makes an accidental column reordering
   visible.
6. **`inst/tinytest/capi/consumer.c` plus `inst/tinytest/test-capi.R`**: a NEW
   partial-`structSize` leg. The existing one (`consumer.c:807`) sets
   `structSize = offsetof(dbarts_forest_calibration, leafModel)`, which
   simulates a pre-`leafModel` caller, not a pre-S1 one; falsifier F6 needs
   `offsetof(..., amplitudePriorVariance)` - the S1 boundary - with the
   `test-capi.R` assertion that the five new buffers are left untouched and the
   original eight are filled.
7. **`docs/design/feature-matrix.md` `[f16]`**, the `nameable calibration`
   footnote, whose body begins at `:368` at tip `ffb9959c` (re-derived at the
   new tip by opening the file, since `5277cd54` re-anchored it wholesale). It
   enumerates what `$getCalibration` reports and which tests carry it -
   including `test-calibration-midchain.R`, which item 5 moves - so its PROSE
   is a required S1 edit. No cell moves; the capability stays `S`.
8. Rd and design docs per "Shipped-surface deltas".

Gates: `R CMD INSTALL --preclean` into a PRIVATE library (headers change);
delete `benchmarks/kernels` binaries; `tests/cpp` from clean, plain AND under
ASAN+UBSAN (`ASAN_OPTIONS=detect_container_overflow=0`) - S1 is
engine-reachable and adds a vector; full tinytest; `test-capi.R` including the
flat-C leg and the unchanged `dbarts_apiHash()` literal; the trio BITWISE with
the ABORT clause above; full `R CMD check` from a clean copy staged outside the
tree (R/Rd touched); `air format --check`; `lintr` on touched R files;
`NEWS.Rd` parses.

Budget: engine 20-34 (stop 48), bridge 10-20 (stop 30), flat C 12-24 (stop 36),
R 2-8, docs 50-85, tests 100-160 (stop 215).

### S2. The defaults, the documentation debt, and the K = 1 refusal

1. `R/model.R`: `defaultAmplitudePriorScale(family)` beside `defaultNodeScale`
   (gaussian 2, probit and logistic 1; no other family reaches it, since
   `R/spec.R`'s family gate refuses every other one at the door);
   `forestParams(specs, hasBasis, family)` with the two new defaults.
2. `R/spec.R`: pass `family` at the one call site.
3. `R/bartcore.R`: `sd.control = NULL`, with the effective family RESOLVED
   FIRST. Two traps, both live and both fatal to a literal reading of "resolve
   after the family is known": the family is currently computed at the BOTTOM
   of `bartcoreBCFSampler` (`model <- sampler$model; if (!is.null(family))
   model@family <- family`) while `bcfParams` is built at the TOP, and a NULL
   left in `as.double(c(0, 0, 0, 1, 1, 1, sd.control, update.a))` COLLAPSES the
   vector to length 7 and surfaces as the bridge's "must be a length-8 numeric
   vector per forest". So the implementer MOVES the family resolution above
   `bcfParams`, and the helper ERRORS on an unknown family rather than
   returning `switch()`'s invisible NULL - `bartcoreBCFSampler` has no R-side
   family gate at all (its refusal is `refusedBCFFamilyReason` at the bridge),
   so the R helper cannot borrow the public path's gate.
4. NO C twin for the new helper, and say why in its comment: unlike
   `defaultNodeScale`, nothing backstops it - `applyBCFSpec` always receives
   explicit per-forest params, so there is no route on which a default could be
   silently taken engine-side.
5. `R/model.R` `forestBasisDeclarations` and `R/spec.R`'s `data@bases` block:
   Fork 6c's two-part K = 1 fix (drop the `< 2L` floor, add the designed
   refusal at the site both routes reach).
6. Docs, Rd, NEWS per "Shipped-surface deltas", including Fork 5's successor
   sentence in its cap-not-pin form.

Gates: `R CMD INSTALL` into a private library (no header edit, so no
`--preclean` is required - but re-run S1's `--preclean` install if the two land
close together); `tests/cpp` (unchanged, as a leak check); full tinytest; the
trio with the S2 verdicts above; `bcf-exact.R quick`; full `R CMD check`; `air
format --check`; `lintr`; `NEWS.Rd` parses; `pkgdown::check_pkgdown` (no new
topic is added, so no `_pkgdown.yml` entry is owed).

Budget: R 16-32 (stop 46), docs 60-110 (stop 150), tests 75-135 (stop 185).

## Test design: the product pin

### The oracle, stated exactly

For a forest f of a K-forest sampler under a latent family:

    $getCalibration(f)[, "prior.scale"] == sd_f * s / (0.674 * c_f)

with `s = 1` (probit) or `pi/sqrt(3) = 1.813799364` (logistic), `c_f` the
median of the NONZERO row norms of `data@bases[[f]]`, and `sd_f` the value
declared by `forest(sd = )`. Exact, at tolerance 1e-12: under a latent family
`response_->fitScale()` is 1 and `k` is 1, so `Chain::forestCalibration`'s
`priorScale = leaf.scale * fitScale * sqrt(numTrees)` cancels the map's
`sqrt(m)` exactly and the reported number IS the map's node scale.

### The fixture, and why each value is what it is

    bases   = list(3 * ones, 5 * zBasis)      # c = (3, 5)
    forests = list(forest(sd = 2.5), forest(sd = 0.4))

Every factor discriminates and no two are equal: `sd` in {2.5, 0.4}, `s` in
{1, 1.813799}, `c` in {3, 5}, divisor 0.674. Resulting oracle values,
re-derived:

| forest | probit | logistic |
|---|---|---|
| 1 | 1.236399604352 | 2.242580816313 |
| 2 | 0.118694362018 | 0.215287758366 |

The mutations each assertion catches, and the one it cannot:

- `sd_f` dropped (anchor used alone): forest 1 reads 0.494715 instead of
  1.236400. CAUGHT.
- forest 1's `sd` used for forest 2: 0.741840 instead of 0.118694. CAUGHT.
- divisor dropped: 0.833333 instead of 1.236400. CAUGHT.
- row norm dropped: 3.709199 instead of 1.236400. CAUGHT.
- forest 1's `c` used for forest 2: 0.197824 instead of 0.118694. CAUGHT.
- **anchor dropped**: UNCAUGHT under probit, where `s = 1`. Caught under
  logistic (1.236400 against 2.242581, a 45 percent miss). **The pin is
  therefore invalid unless it runs under both latent families** - write that
  reason into the fixture comment, because it is exactly the class of
  near-miss M4.4 recorded (the naive anchor's logistic arm missing by 0.865
  percent).
- divisor and row norm interchanged: UNCAUGHT and harmless - both sit in the
  denominator and only the product is claimed.

### Where each cell lives

1. `inst/tinytest/test-bcf-family.R`, beside the existing anchor block (it
   already carries `medianRowNorm`, the per-family loop and the non-unit-norm
   `scaledBases` cell): the fixture above, both families, at 1e-12. This
   EXTENDS the `1.0 * s / (0.674 * ...)` assertion at `:109-113` whose literal
   `1.0` is the hole - keep that cell too, since it pins the DEFAULT while the
   new one pins the DECLARED value.
2. `inst/tinytest/test-forest-basis-r5.R`, in the `probitForests()` block
   (`:335-372`): the same fixture, then `$setForestBasis(2, 7 * zBasis)`, then
   `prior.scale_2 == 0.4 * s / (0.674 * 7)` = 0.084781687156 under probit.
   This is the only test of `Chain::setForestBasis`'s leaf-scale re-derivation
   and today it runs at `nodeScaleFactor = 1`, so it cannot see the stored
   factor being lost on a swap.
3. `tests/cpp`, a new `testBCFCalibrationMap` (natural home:
   `tests/cpp/test_sampler.cpp`, beside the BCF combiner cases; `test_state.cpp`
   already carries the "BCF's two forests calibrate differently" leaf-scale
   block, and `test_sampler.cpp:5982-6009` a BCF calibration READER block that
   asserts positivity, `k == 1` and the tree-count factor but pins no map
   value, as the nearest neighbours). Build through `createBCFSampler`
   (`src/bartcore/facade.hpp`) with `spec.forests` of length THREE - factors
   {2.5, 0.4, 1.75}, divisors {0.674, 0.8, 0.25}, bases with row norms
   {3, 5, 6} - under `ResponseFamily::probit` and again under `logistic`, and
   assert `sampler->forestCalibration(0, f).priorScale` against the literal
   expression at 1e-12. K = 3 rather than 2 so the fixture also covers a forest
   the two-forest spelling cannot describe.

   **The fixture's own admissibility check, which belongs in its comment
   because this is the third time these values have been wrong.** All nine
   values distinct and none 1; the three DENOMINATORS `divisor * rowNorm` =
   2.022, 4.0, 1.5 distinct, none 1, and none equal to any of the nine values;
   the resulting `prior.scale` under probit = 1.236400, 0.100000, 1.166667 and
   under logistic = 2.242581, 0.181380, 2.116266, each distinct from the other
   two, from `s`, and from its own forest's factor, divisor and row norm. The
   first draft gave forest 3 divisor 1.0 and row norm 1 (both factors
   invisible); the first revision gave it 0.5 and 2, whose PRODUCT is 1.0, so
   the denominator as a whole was still invisible there AND its probit
   `prior.scale` came out at exactly 1.75, its own `nodeScaleFactor` - which
   after S1 adds a `node.scale.factor` column would hide a bridge that filled
   one buffer from the other. Two rounds of the arc's own vacated-pin hazard
   inside the fixture built to close it; the check above is what stops a third.
4. `tests/cpp`, the row-norm CONVENTION, which no test anywhere exercises:
   a basis with an EVEN number of nonzero rows (so the average-of-two-central
   -order-statistics branch runs) and at least one all-zero row (so the
   exclusion rule runs), plus an all-zero basis asserting the 1.0 fallback.
   `Chain::basisRowNorm` is private, so assert it THROUGH `forestCalibration`.
5. `inst/tinytest/test-bcf-creation.R`: transport pins that
   `params[[f]][4]` carries the declared `sd` for a basis forest and
   `params[[1]][7]` the declared `sd` for a basis-free one - extending
   `test-bcf-creation.R:400` (`forest(sd = 3.5)` into `params[[1]][7]`) and
   `test-forest-basis-r5.R:196-205` (the `params[[f]][6:7]` ridge derivation),
   both of which are gaussian - under a latent family for the first time.

### What S2 adds on top

- The K-invariance assertion. `test-bcf-family.R:118-125` currently asserts
  `sqrt(sum(scales^2 * 0.5)) == 1.04912 * sqrt(K) * s` at K = 2 and 3. Under
  the law it becomes `== 1.04912 * sqrt(2) * s` - a CONSTANT - at K = 2, 3 and
  4, which is strictly stronger and is the design statement itself.
- The default-transport pins: `params[[1]][7]` is 2 under gaussian and 1 under
  probit and logistic; `params[[f]][4]` is `sqrt(2/K)` for a basis forest at
  K = 2, 3, 4 (1, 0.816497, 0.707107) and unchanged when `sd` is declared.
  **The GAUSSIAN arm of this pin is mandatory and its reason goes in the
  comment**: after S2 a basis-free forest's params vector under a latent family
  is `c(50, 0.25, 3, 1, 1, 1, 1, 1)` - positions 4 through 8 all carry the
  literal 1 - so a latent-only fixture cannot see a transport confusion among
  slots 4-8 at all. Only the gaussian arm, where slot 7 is 2, discriminates
  position. This is the same "unit values silently vacate pins" hazard the
  slice exists to close, arriving through the slice's own new default.
- The internal-route agreement: `bartcoreBCFSampler` on a probit host and the
  public two-forest probit call report the same `prior.scale` AND the same
  amplitude prior scale.
- The K = 1 refusal, by message.
- The `scaledBases` cell carrying `4 * zBasis` is LOAD-BEARING (M4.4's landing
  note: delete it and the dropped-divisor mutation goes green). Neither slice
  may tidy it away; say so in a comment beside it.

## Falsifiers (pre-registered)

- **F0 (S1, STOP).** Any non-bitwise scenario on any of the three baselines at
  S1 means a read-only change moved a draw: stop, do not re-record, and
  diagnose along BOTH admissible branches, not just the first. (i) The
  ARITHMETIC branch: an expression was re-associated or a stored value reused
  where the source recomputes - re-derive the leaf-scale expression
  character by character. (ii) The CODEGEN branch: the expression is provably
  identical and the object code is not, which this codebase has already met
  once (`BCFSpec::generalAmplitudeDraw`'s comment records four accumulators
  agreeing in exact arithmetic and differing in where FMAs form, over 21
  variants). Confirm (ii) by rebuilding the PRE-slice tree with the new header
  present but unreferenced, and treat a divergence that survives that as a
  finding about the build, not about the design. A falsifier admitting only
  diagnosis (i) would misclassify (ii) as a design defect. Fork 3b's spec-echo
  choice exists to keep branch (ii) empty: no virtual is added, so no vtable
  moves.
- **F1 (S2, STOP).** Any non-bitwise scenario on `bcf-equivalence` at S2 means
  the family gate or the K gate leaked into gaussian K = 2. Stop.
- **F2 (S2).** A gaussian K = 2 fit built through `dbarts(..., forests = ...)`
  before and after S2 is BITWISE identical, on the public route as well as the
  internal one. The trio only covers the internal one.
- **F3 (S2).** A latent-family fit that DECLARES `forest(sd = 1)` on every
  basis forest and `forest(sd = 2)` on the basis-free one is BITWISE identical
  before and after S2 at every K. This is the composability claim: the law is a
  default, not map algebra, and a caller who states the old values gets the old
  model.
- **F4 (S1).** `$getCalibration` on a single-forest sampler and on a
  multinomial sampler returns NaN in exactly the five new columns and is
  otherwise unchanged, including the `leaf.model` attribute. On a K-forest
  sampler the two amplitude columns are EXCLUSIVE per forest: exactly one of
  `amplitude.prior.variance` and `amplitude.prior.scale` is finite on each,
  and which one it is agrees with `params[[f]][6:7]`.
- **F5 (S1).** The reported `basis.row.norm` changes after `$setForestBasis`
  and equals the R-side `median(norms[norms > 0])` of the new basis, including
  a basis with an even nonzero-row count. A stored-and-stale value is the
  failure this catches.
- **F6 (S1).** A flat-C caller with a PRE-S1 `structSize` - set to
  `offsetof(dbarts_forest_calibration, amplitudePriorVariance)`, the S1
  boundary, in a NEW `inst/tinytest/capi/consumer.c` leg beside the existing
  pre-`leafModel` one at `:807` - still reads the EIGHT original fields and
  leaves the five new buffers untouched (poison them and assert they survive,
  as the existing leg does). `DBARTS_HAS_FIELD`'s contract, asserted rather
  than assumed.
- **F8 (S1), the truthfulness falsifier Fork 3b's rules owe.** Three arms, all
  on K-forest samplers, all reachable today. (a) FOREIGN CALIBRATION: with
  `A <- dbarts(..., forests = list(forest(), forest(basis = ~z, sd = 2)))` and
  `B` the same at `sd = 0.5`, `B$setState(A$state)` is ACCEPTED (same widths,
  same tree counts - neither `glueIsValid` nor the forest gate looks at
  calibration), and afterwards `B$getCalibration(2)` reports `NA` in
  `node.scale.factor` and `node.scale.divisor` while `prior.scale` reflects the
  donor. Assert the NA and assert the identity is NOT computable, so nobody
  re-derives a wrong `s`. (b) SELF-RESTORE: store, run, restore - every column
  survives non-NA and the identity still holds, because the installed scale is
  bitwise the one in force. (c) RE-IMPOSITION: after (a), a `$setForestBasis`
  on forest 2 restores both columns and the identity, since the map is
  re-derived from the stored factor and divisor. Also assert (d): after (a) the
  `amplitude.prior.variance` column equals the DONOR's, matching what
  `drawAmplitudes` will use. Without (a) and (d) the rules are unpinned and
  the reader is the stale one Fork 3a was rejected for.
- **F7 (S0/S2).** The product pin FAILS under each of the five mutations
  enumerated in "Test design" and PASSES as shipped. Substantiate by mutation
  build, as M4.4 substantiated item 25; report the assertion count moved by
  each, and record any mutation caught only by the logistic arm.

## Edge cases the tests must name

A forest whose basis has an all-zero row (excluded from the median, not counted
small). An all-zero basis (`basisRowNorm` returns 1.0; the forest contributes
nothing and its reported scale is still finite). A basis-free forest under the
law (factor stays the literal 1; the law touches only the `withBasis` branch).
`forest(sd = )` declared on a basis-free forest under a latent family (it is
the half-Cauchy median, now defaulting to 1, and a declared value overrides in
both families). K = 1 (refused by name, with the message naming the two-forest
spelling). A `dbartsData(bases = )` object with NO `forests =` declaration at
K >= 3 (every forest at the new default - this is the route that most silently
inherits the law). `$getCalibration` on a K-forest sampler between sweeps: the
two amplitude columns are EXCLUSIVE per forest and the reported variance does
NOT move with the scale mixture's auxiliary, because the reader carries the
PRIOR and the auxiliary is state (reachable through `$state`'s `aVariance`,
recorded as a door). A state install that brings a foreign calibration
(`$setState`, `$installTrees`), and the `$setForestBasis` that re-imposes the
map afterwards - F8's three arms, and the one edge case where a column
legitimately becomes NA on a K-forest sampler.

## Doors held open (recorded, not scheduled)

- **K = 1 as a shipped configuration** (Fork 6b), with its D4 enabling value
  named. New ticket.
- **`$getIndexPrior()`**, a reader returning `sd(eta_i)` quantiles (Fork 3c).
- **A `node.scale.anchor` column**, needed only because gaussian's anchor is
  data-dependent; recoverable by identity until then.
- **`bcf-naming-generalization`**, which owns `calibrationMapName` and every
  `BCF*` spelling this slice writes prose around.
- **A family-aware `amplitude.prior.variance`**: not proposed, and the reason
  is recorded so it is not rediscovered - `v` is the 0.674 divisor's partner,
  not a free scale.

## NEWS bullets (inst/NEWS.Rd, one per slice, same commit)

- S1: `$getCalibration` reports the five calibration-map quantities that were
  in force but unreadable - the forest's amplitude prior variance or
  half-Cauchy scale, its node scale factor and divisor, and the basis row norm
  the map divides out - so the induced prior on a multi-forest model's combined
  location can be reconstructed from the sampler rather than from the call. The
  flat C `dbarts_forest_calibration` gains the matching fields.
- S2: the prior defaults for a multi-forest model change. A forest with no
  basis takes a half-Cauchy median of 1 rather than 2 under `"probit"` and
  `"logistic"` (unchanged at 2 under a gaussian response, where a drawn sigma
  absorbs the difference), and a forest with a basis takes `sqrt(2/K)` rather
  than 1 for K forests, so that declaring more forests no longer widens the
  prior on the combined location without bound: with every forest carrying a
  basis it is now the same at every K, and otherwise it is bounded by that
  value instead of growing. Declare `forest(sd = )` for the previous values;
  results move for binary multi-forest fits and for any fit with three or more
  forests.

## Discrepancies (against the census)

Recorded where this design disagrees with the verified census
(`.claude/binary-kforest-prior/census.md`, GIT-IGNORED session evidence; every
fact it carries that this plan relies on is restated here) or adds to it. The
census's symbol citations all held; these are corrections of detail and of
consequence.

1. **`dbarts_forest_calibration` carries EIGHT field pointers, not seven.**
   The census (leg (c), and section 4.3) says "the same seven pointers plus
   `structSize`". Verified: `priorScale`, `priorSd`, `priorMean`, `k`,
   `responseScale`, `responseShift`, `kHasHyperprior`, `leafModel`, each with
   its own `offsetof` static_assert in `src/C_interface.cpp`, and a `sizeof`
   assert of `sizeof(size_t) + 8 * sizeof(double*)`. The R getter's SEVEN
   columns plus a `leaf.model` attribute is right; the flat twin carries the
   leaf model as a field. This is load-bearing for Fork 3's cost - the append
   touches nine asserts, not eight.
2. **The flat struct is DESIGNED for appends.** The census reads the flat twin
   as "the same gap at the LinkingTo surface", which understates it: the header
   carries an explicit "1.0-0 field boundary: appends go below, never above"
   comment, `DBARTS_HAS_FIELD` gates every fill on the caller's `structSize`,
   and the size assert's own message tells an appending author what to do. So
   Fork 3's flat leg is a versioned append, not a break, and the `apiHash` does
   not move.
3. **The census's leg (a) site list is complete for the "2" and incomplete for
   the map's other defaults.** `BCFSpec::bPriorVariance = 0.5` and
   `BCFSpec::sdModerate = 1.0` (`src/bartcore/combiner.hpp`) are two further
   engine spellings of shipped defaults. Neither moves here - both are the
   K = 2 values - but a K law defined on the wrong count would have hit
   `sdModerate`, which is why Fork 2d is priced.
4. **"A change to `forestParams` alone turns the creation oracle red" is true
   only under gaussian.** The census says so implicitly (its section 5 marks
   every `sd` site gaussian) but the CONSEQUENCE is not drawn: it is exactly
   what makes Fork 1c cheap and Fork 1b expensive, and it is why this slice
   re-records nothing.
5. **The census's open question 1 (which comparator) is answered by M4.4's own
   mandate.** The undischarged sentence names 0.238 - the hyperprior arm - so
   the primary comparator is settled by the obligation this slice inherits, not
   chosen freshly here. The fixed-k arm survives as the analytic denominator.
6. **The census's K = 1 statement is right about one route and wrong about the
   other.** It says `dbarts(x, y, forests = list(forest(basis = ~z)))` "reaches
   a bridge-level error message rather than a designed refusal". Re-traced by
   opening `forestBasisDeclarations` (`R/model.R:689-692`, floor
   `length(forests) < 2L`), `R/dbarts.R:542-561` and `R/spec.R:206-226`: that
   call reaches NO error. The basis never reaches `data@bases`, the bcf block
   never runs, and an ordinary single-forest model is fit with the declared
   basis SILENTLY DROPPED. The bridge error belongs to the DATA route
   (`dbartsData(x, y, bases = list(b))`). Fork 6 is re-scoped around both.
7. **The observability leg's own budget escalates more than 10x over the
   TODO's estimate, and the estimate was for a design that does not survive.**
   `TODO` prices it "R 8-14 dense, bridge 0-4"; S1 here is roughly 145-270
   dense-equivalent across five layers. The cause is Fork 3a's refutation: the
   ticket's figure assumed an R-side reader off
   `attr(control, "bartcore.bcf")$params`, which cannot report the row norm in
   force because `Chain::setForestBasis` re-derives it. The escalation is a
   consequence of correctness, not of scope creep, but it is a discrepancy
   against the ticket and is recorded as one.
8. **Every number in the census that this design places load on was
   re-derived and reproduces**: `sqrt(0.5)/0.674 = 1.049119853`; 1.48368 at
   K = 2 and 2.96736 at K = 8; 1.97824 against the fixed-k comparator; the
   shipped K = 2 probit shape at med 1.5206 / robust sd 2.2544 / tail 0.3764 /
   `P(|eta| > 6)` 0.1643 (census: 1.519 / 2.252 / 0.3759 / 0.1643); the
   single-forest hyperprior arm at 1.0689 / 1.5847 / 0.2387 (census: 1.069 /
   1.584 / 0.2383); `k = 2 chi(1.5)` median 1.9049, mean 2.092, 95% (0.229,
   5.013) against analytic 1.906131 / 2.092099 / (0.229024, 5.010313);
   `a mu` median 1.1749 and 95th percentile 20.26 (census: 1.175 and
   20.17). No number failed to reproduce.

## Appendix: reproducibility

Every simulated figure in this document: R 4.6.x, base only, `set.seed(20260814)`
once at the top, N = 2e6 draws per arm, arms drawn in the order they appear.
The generative model is exact rather than approximate, and that is why no
sampler is involved: a forest of `m` trees with leaf prior
`N(0, (nodeScale/sqrt(m))^2)` has total `N(0, nodeScale^2)` at any x, so an arm
is a product of two draws per forest -

    basis forest f:      rnorm(N, 0, sqrt(v_f)) * rnorm(N, 0, sd_f * s / (0.674 * c_f))
    basis-free forest:   rcauchy(N, 0, aPriorScale) * rnorm(N, 0, s)
    single-forest ref:   rnorm(N, 0, 3 * s / k),  k = 2 * sqrt(rchisq(N, 1.5))

- summed over forests, with `p = pnorm(eta)` under probit and `plogis(eta)`
under logistic. "Robust sd" is the interquartile range over `2 * qnorm(0.75)`,
used because every arm carrying a half-Cauchy has infinite variance. The
literal `0.674` is the code's, not `qnorm(0.75) = 0.674489750`. The gaussian
range-to-sd figures in Fork 2e are the mean over 200 replicates of
`0.5 / sd((y - min(y))/(max(y) - min(y)) - 0.5)`, which is the ratio the
engine's `scaledResponseSd()` produces.

Independently re-checked at N = 4e6 on a different seed by the blind critique:
no figure failed to reproduce, and the two closed forms it added (analytic
`k` quantiles, and `eta = (3/(2 sqrt(1.5))) t_{1.5}` for the primary
comparator) are folded in above.

## Critique adjudication

Blind critique at `.claude/binary-kforest-prior/critique.md` (GIT-IGNORED
session evidence; both rounds' findings and their adjudications are recorded in
full below, so nothing load-bearing lives only there), read against tree
`55cc1756`; this document is now at `ffb9959c` (three docs-only commits; no
code moved, so every code symbol above is unaffected and the one
`feature-matrix.md` citation was re-derived by opening the file). Verdict
counts: 0 BLOCKER, 4 MAJOR, 8 MINOR. **All twelve ADOPTED; none contested.**
Two of them found more when re-traced than the finding claimed, and those
extensions are recorded here because they change what the slice does.

**M1 (pin/cap overclaim) - ADOPTED, in three places plus a fourth the critique
did not list.** The Goal item 2, Fork 5's successor sentence and the S2 NEWS
bullet now state HELD for the all-basis shape and BOUNDED for every other, and
leg (b) carries the exact expression `sqrt((K-1)/K) s / 0.674` with the four
values and their ratios to the `k = 2` budget (0.699 / 0.808 / 0.857 / 0.925,
limit 0.989). The fourth site was Fork 2's own "states the budget in the
comparator's own units" bullet, which asserted the same identity unscoped.
Also folded in: today's shipped-shape channel is 1.850x the budget at K = 8
and passes 2x at K = 10, which is the number that makes the cap worth having
and which the first draft never computed.

**M2 (Fork 2e's rejection refuted by Fork 1c) - ADOPTED; the reasoning is
replaced and the RECOMMENDATION STANDS.** The Rd-sentence argument is
withdrawn: 1c family-branches that item anyway, so 2e costs one clause. The
replacement is a gaussian coverage argument in gaussian's own units -
`s = scaledResponseSd()` is `sd(y)` internally, so the all-basis default puts
the prior sd of `E[y|x]` at exactly `1.049120 sqrt(K)` times `sd(y)`, i.e.
2.97 `sd(y)` at K = 8, and a prior asserting the mean varies three times as
much as the response is not a defensible default; "a drawn sigma absorbs it" is
a posterior-correction claim, not a defence of the prior. The single-forest
gaussian comparator is included as CONTEXT with its data dependence measured
(`0.5/s` = 2.25 / 2.73 / 3.24 / 3.66 for normal y at n = 50 / 200 / 1000 /
5000, 1.73 uniform, 4.80 t3) rather than as the argument, because it is not
scale-free. The gate-absence sentence is struck as an argument and restated as
a migration fact, per the standing rule that gate absence never argues for a
change. 2e is now presented to the owner as a LIVE fallback with its cost
priced at one Rd clause, since the gaussian argument is the whole of what
stands between the two branches.

**M3 (refusal site unreachable from one route) - ADOPTED, and the finding's own
trace is CORRECTED in the adopting direction.** The critique is right that
`resolveForests` returns at `R/model.R:811` on `forests = NULL` and so cannot
see the `dbartsData(bases = )` route. Re-tracing the OTHER route by opening
`forestBasisDeclarations` (`R/model.R:689-692`) shows it does not reach the
bridge error either: its `length(forests) < 2L` floor means
`dbarts(x, y, forests = list(forest(basis = ~z)))` never puts the basis on the
data object, never enters the bcf block, and fits an ordinary single-forest
model with the declared basis SILENTLY DROPPED - a silent-wrong-model defect
on a public route, worse than the error the fork was about, and in neither the
census nor the critique. Fork 6c is therefore two parts: drop the floor to
`< 1L` (verified not to disturb `test-bcf-creation.R:713`'s existing refusal,
since a length-1 list declaring no basis still expands to `list(NULL)` and
still falls through), then one designed refusal in `R/spec.R`'s
`!is.null(data@bases)` block where `numForests` is computed - the site both
routes now reach, verified by opening both traces.

**M4 (two missing S1 edit sites) - ADOPTED in full.**
`inst/tinytest/test-calibration-midchain.R:41-58`'s `c(2L, 7L)` and exact
`colnames` pins are now S1 item 5, moving to `c(2L, 12L)` in the same commit
and keeping the exact-set form; `inst/tinytest/capi/consumer.c` plus
`test-capi.R` are S1 item 6, with F6 restated to name the S1 boundary
(`offsetof(..., amplitudePriorVariance)`) rather than the existing
pre-`leafModel` leg at `:807`. Both are listed in "Shipped-surface deltas" as
mandatory same-commit edits.

**m1 (latent basis-free params all-ones) - ADOPTED.** The gaussian arm of the
transport pin is now MANDATORY with its reason stated: after S2 a latent
basis-free vector is `c(50, 0.25, 3, 1, 1, 1, 1, 1)`, so positions 4-8 are
mutually indiscriminable and only the gaussian arm (slot 7 = 2) can see a
transport confusion. The slice's own new default reintroduces the arc's
vacated-pin hazard, which is worth the sentence.

**m2 (the `tests/cpp` third forest is itself vacuous) - ADOPTED.** Divisors
become {0.674, 0.674, 0.5} and row norms {3, 5, 2}, so all nine fixture values
are non-unit and no two within a row are equal. The critique is right that the
first draft claimed "unequal per-forest values throughout" while giving the
K = 3 forest divisor 1.0 and row norm 1.

**m3 (`sd.control = NULL` ordering and NULL-collapse traps) - ADOPTED.** S2
step 3 now says to MOVE the family resolution above `bcfParams` and states the
failure mode explicitly (a NULL inside `as.double(c(...))` yields a length-7
vector and the bridge's length-8 error), and requires the helper to ERROR on an
unknown family rather than return `switch()`'s invisible NULL - the critique is
right that `bartcoreBCFSampler` has no R-side family gate to borrow.

**m4 (Fork 6b priced on unverified engine cost) - ADOPTED.** Verified live:
`drawGlue` dispatches K != 2 to the general path, `synthesizeIndicatorBasis` is
already gated on `numForests_ > 1`, `shippedShape()` is simply false at K = 1,
`expandForestSpecs`' synthesis is unreachable from R, and `createBCFSampler`
imposes no count floor - the only floor is `applyBCFSpec`'s. The deferral now
rests on GATE grounds (a new shipped configuration owes acceptance evidence),
and the ticket must carry the corrected pricing rather than the false one.

**m5 (`amplitude.prior.variance` names a refused knob) - ADOPTED, and it
changed the design.** The reader now carries TWO exclusive amplitude columns,
`amplitude.prior.variance` and `amplitude.prior.scale`, each reporting a prior
the caller may set or NaN, mirroring `ForestAmplitudePrior`'s own
`{variance, halfCauchyScale}` pair. Combined with the vtable risk the critique
raised in its "could not refute" section, this replaced the live-auxiliary read
with a SPEC ECHO: no new virtual on `ForestCombiner`, no call into the combiner,
and S1's draw-neutrality becomes structural instead of empirical. The live
scale-mixture auxiliary is state, not calibration, and is now a named door. Net
effect on the slice: five columns rather than four, and one fewer risk.

**m6 (feature-matrix justification covers only S2) - ADOPTED.** The
per-landing check is now stated per slice: vacuous at S2, prose-moving at S1
through `[f16]`, whose body was re-derived at tip `ffb9959c` (`:368`) rather
than carried over from `55cc1756`.

**m7 (numeric slips) - ADOPTED in full, none changing a verdict.** `sqrt(2/16)`
0.353554 -> 0.353553; the K = 8 cap 1.38786 -> 1.3878547; "within 3 percent"
-> 3.4 percent; the `k` quantiles now carry the critique's closed forms
(1.906131 / 2.092099 / (0.229024, 5.010313)) beside the simulated ones; and
leg (a)'s exactness sentence is split - the SD statements are exact, the
DISTRIBUTION of a contribution is a leptokurtic product of normals, which is
why the tail figures are simulated (0.1089 against a normal's 0.1169 at the
same sd).

**m8 (TODO closes with a deliverable deferred; budget escalation) - ADOPTED.**
The `TODO` line now closes WITH a recorded deviation naming the deferred
`sd(eta_i)` reader, and Discrepancies gains entry 7, the >10x escalation over
the ticket's "R 8-14 dense, bridge 0-4" and its cause (that estimate priced
Fork 3a, which the critique itself could not refute the refutation of).

**The critique's residual risk inside "could not refute" - ADOPTED as two
changes.** The vtable/inlining hazard is removed at the source by Fork 3b's
spec echo, and F0 now admits two diagnostic branches (arithmetic
re-association, and codegen) with a stated procedure for separating them,
rather than the single diagnosis that would have misclassified a codegen
artifact as a design defect.

**Nothing contested.** Each finding was re-verified by opening the cited target
before adoption, and the two that were incomplete (M3's trace, m5's remedy)
were extended rather than argued with.

### Round 2

Second blind pass over the revised content only, at tip `ffb9959c`: 0 BLOCKER,
2 MAJOR, 5 MINOR. **All seven ADOPTED; none contested.** Fork 6c survived both
of its claims (including the new silent-drop defect, reproduced independently
and checked against all ten length-1 `forests` fixtures in the suite and
against every consumer), M4's 7 + 5 = 12 arithmetic survived at every site,
M1's four sites survived, and Fork 3b's draw-neutrality survived as structural.

**R2-M1 (the spec echo is not truthful after `$setState` or `$installTrees`) -
ADOPTED; this is the round's real design change.** Verified by opening both
paths: `installForest` (`chain.hpp:3235`) and the restore path (`:3323`) each
install a donor `leaf.scale` BY DESIGN - `ForestStateData::leafScale`'s doc
block states the rationale - and each then calls `restoreGlue`, which
overwrites `glue_.prior[f].variance` (`combiner.hpp:1044-1045`); neither
mutator carries a BCF refusal (`R/dbarts.R:1580`, `:1612`) and
`test-forest-basis-r5.R:110-136` already drives `setState` on a K-forest
sampler. A pure echo would therefore print the recipient's decomposition beside
a donor-derived `prior.scale` - Fork 3a's own failure mode, one call deeper.
The fix is Fork 3b's three rules: the amplitude columns FOLLOW the state
(updated at both install sites under `restoreGlue`'s own guard, using only the
EXISTING `glueIsValid` virtual, so nothing about the no-new-virtual argument
moves); the two map-decomposition columns carry a `nodeScaleIsMapDerived_` bit,
cleared only when an install brings a scale that differs BITWISE from the one
in force - so a self-restore keeps its columns - and restored by
`$setForestBasis`; `basis.row.norm` needs no rule. The Rd identity is stated
with its condition ("whenever `node.scale.factor` is not `NA`") rather than
narrowed by prose, and falsifier F8 pins all four arms. The three rejected
alternatives are recorded in Fork 3b with their grounds: narrowing the
contract (a documented lie, which is what Fork 3a was rejected for), deriving a
factor from the installed scale (invents a decomposition the state cannot
carry), and refusing the two mutators on K-forest samplers (deletes shipped,
tested, deliberately designed restore-then-widen behaviour to tidy a getter).

**R2-M2 (the gaussian argument permits but does not support) - ADOPTED; the
argument is replaced and the recommendation stands.** The counterexample
reproduces: re-derived at 2000 replicates, the SHIPPED single-forest gaussian
default puts the prior sd of `E[y|x]` at 3.940 `sd(y)` for t3 data at n = 1000
and 6.795 at n = 5000 - above the 2.97 the first revision called indefensible -
so no magnitude threshold separates the K = 8 default from this package's own
conventions, and the first revision reached its band only by excluding that arm
with an adjective. The magnitude argument is WITHDRAWN for gaussian and the
table is now IN Fork 2e, so the owner sees the counterexample. What replaces it
is the critique's own suggested discriminator, which is also the argument the
law embodies: COMPOSITIONALITY - the single-forest default is one ensemble's
whole budget however large, while `sqrt(K)` compounds it across ensembles, so
today the prior on `E[y|x]` depends on how the user decomposed the mean rather
than on what they said about it, and writing one model as two forests or four
changes it by `sqrt(2)`. Family-free, threshold-free, comparator-free. The
magnitude argument is RETAINED for the binary families, where `s` is fixed by
the link rather than by the data and the comparator does not move. Fork 2's
scope is therefore UNCHANGED at all families, on a different and better basis,
with 2e presented as a live fork whose tradeoff is now stated in both
directions and whose tiebreak is named as the release window rather than as
evidence.

**R2-m1 (the revised fixture restores a unit denominator) - ADOPTED.** `0.5 * 2
= 1.0` made forest 3's denominator invisible as a whole, and its probit
`prior.scale` came out at exactly 1.75, its own `nodeScaleFactor` - which after
S1's `node.scale.factor` column would hide a cross-buffer wiring error.
Divisors are now {0.674, 0.8, 0.25} against norms {3, 5, 6}, giving denominators
2.022 / 4.0 / 1.5, and the fixture carries an explicit admissibility check in
its comment (all nine values distinct and non-unit; the three denominators
distinct, non-unit, and distinct from the nine; every `prior.scale` distinct
from its own factor, divisor, norm and from `s`) so a fourth round of this is
not needed.

**R2-m2 (the adjudication claimed a fix it did not make) - ADOPTED.** The wider
claim was true of the fault and false of the remedy: `{0.674, 0.674, 0.5}` did
leave two divisors equal. The new `{0.674, 0.8, 0.25}` makes the wider claim
true, so both sentences now agree.

**R2-m3 (a stale edge case now contradicting Fork 3b and F4) - ADOPTED.** The
"the reported variance MOVES" line is replaced by the exclusivity case, the
door for the live auxiliary (`$state`'s `aVariance`), and F8's restore arms -
the one place a K-forest column legitimately reads `NA`.

**R2-m4 ("the only behaviour that moves" too strong) - ADOPTED.** Fork 6c
fix 1 now enumerates three moves, not one, with the two the critique found: the
message change on `forests = list(forest(basis =, sd =))`, and the length-1
declaration over a data object already carrying K >= 2 bases, where the count
drops to 1 and the refusal would otherwise say "K = 1" to a user who wrote two
bases. The refusal text is specified against the RESOLVED count, naming where
that count came from. Neither is pinned by a fixture, so nothing goes red.

**R2-m5 (residual digit slips) - ADOPTED.** The cap table is now 1.0491199 /
1.2114193 / 1.2849042 / 1.3878551 (re-derived), and `u * sqrt(7) / 1.5` is
1.850473 rather than 1.85042. The limit 1.4836795 and every ratio were already
right.

**FROZEN.** No further self-initiated edits to this document.

## Landing note, S0 (appended 2026-08-15)

LANDED as 4dbf2dbc, tests only, four files: test-bcf-family.R cell (d)
(the product at a declared forest(sd =) under BOTH latent families),
test-forest-basis-r5.R (the same fixture across a setForestBasis swap -
the sole assertion that sees a stored factor lost on an install),
test-bcf-creation.R (ragged transport under a latent family), and
tests/cpp testBCFCalibrationMap (K = 3, nine all-distinct non-unit
values, plus the row-norm convention arms), the first calibration-map
coverage tests/cpp has carried. Budget ~161 dense-equivalent against
the 120-190 band.

Battery on the slice's own private lib: preclean install clean;
tests/cpp all ok from make clean and again under ASAN+UBSAN
(detect_container_overflow=0); full tinytest 4740 results 0 failures;
equivalence trio BITWISE on all three baselines (equivalence-8b047f8b
37/37, bcf-equivalence-8b047f8b 12/12, multinomial-equivalence-1027be5
10/10 - no max-|z| line anywhere); air and lintr clean on the three
touched R files; R CMD check --as-cran from a clean staged tarball
Status OK 0/0/0. Independently re-verified by the orchestrator:
tests/cpp re-run from make clean, diff-vs-plan review, oracle
arithmetic re-derived by hand.

Mutation proof (fresh mutated preclean install per arm; every arm run):
factor dropped, factor misattributed, divisor dropped, row norm
dropped, row norm misattributed, anchor dropped, stored-factor-lost-
on-swap - each moved its named assertions; the anchor arm moved the
LOGISTIC pins only, probit blind exactly as the plan states; the swap
arm moved only the new r5 assertion. Nothing stayed green under its
own falsifier.

Deviations, both recorded: (a) the plan's logistic oracle for the
tests/cpp third forest reads 2.116266; the correct value is 2.1160993
(1.75/(0.25*6) * pi/sqrt(3)). The shipped test asserts the EXPRESSION,
so no test value depends on the slip; the plan body stays frozen and
this note is the correction of record. (b) Test-design cells 3 and 4
ship as one testBCFCalibrationMap with two arms, and the row-norm
convention arm runs under both families rather than probit alone -
strictly more coverage, ~10 fewer lines.

## Landing note, S1 (appended 2026-08-15)

LANDED as 0faeb416. Five layers as specced. Engine: ForestCalibration gains the
five map fields (quiet NaN by default), Chain gains
amplitudePriorVariances_/amplitudePriorScales_/basisRowNorms_ and the
nodeScaleIsMapDerived_ flag, forestCalibration fills the five under the
f < nodeScaleFactors_.size() guard. Bridge: 12 columns. Flat C: five
appended double* fields below the 1.0-0 boundary, five offsetof
asserts, the sizeof assert at 13 pointers, five FILL lines - apiHash
unmoved, sizeof moved. R: the $getCalibration docstring and
bartcoreForestCalibration's. Both mandatory same-commit edits landed:
test-calibration-midchain's dim/colnames pins moved to c(2L, 12L) and
the twelve-name exact set (two further shape pins in the same file, the
BCF c(2L, 7L) at :375 and the multinomial c(1L, 7L) at :391, moved with
them - not in the plan's list, found by the gate), and consumer.c's new
PRE-APPEND leg.

Fork 3b's three truthfulness rules are engine work at the two install
sites, as specced: noteInstalledLeafScale compares BEFORE the
assignment (so a self-restore keeps its columns) and
adoptInstalledAmplitudePriors runs under restoreGlue's own guard
through the EXISTING glueIsValid virtual. One refinement the plan does
not state and the exclusivity sub-decision requires: the amplitude
update writes only forests whose stored variance is non-NaN, i.e. the
FIXED-VARIANCE ones. A scale-mixture forest's serialized
amplitudeVariance is the live inverse-gamma auxiliary, not a prior, so
copying it would have reported a moving quantity under a prior's name
and broken the two columns' exclusivity - the collision m5 was adopted
to close.

DRAW-NEUTRALITY, both written conditions honored and cited: (1) the
leaf-scale expression keeps its exact written shape - the constructor
computes `double c = basisRowNorm(...); basisRowNorms_[f] = c;` and
passes `nodeScaleFactor * s / (nodeScaleDivisor * c)`, the same
association and operand order, the call merely NAMED rather than
recomputed; (2) setForestBasis writes the stored norm from the SAME
call it already makes, not a second one. No virtual joined
ForestCombiner, so F0's codegen branch stays empty by construction.

Battery on the slice's own private lib: preclean install clean;
tests/cpp all ok from make clean (242 ok lines) and again under
ASAN+UBSAN (detect_container_overflow=0), no sanitizer diagnostic;
full tinytest 4826 results 0 failures; equivalence trio BITWISE on all
three baselines (equivalence-8b047f8b 37 compared / 0 skipped under
--strict-coverage, every scenario "identical draws (same RNG stream)";
bcf-equivalence-8b047f8b 12/12 identical on every channel;
multinomial-equivalence-1027be5 10/10 identical on every channel - no
max-|z| line anywhere); air format --check clean; lintr clean on every
touched R file; R CMD check --as-cran from a clean staged tarball
outside the tree Status OK 0/0/0. No Rd topic is new, so no
_pkgdown.yml entry was owed.

Mutation proof (fresh mutated --preclean install per arm; the tree
restored and byte-verified after each): F8(a) foreign calibration, flag
never cleared -> test-forest-basis-r5 :453, :454, :456-461 move (the
two NaN columns and the non-computable anchor); F8(b) cleared on ANY
install, the bitwise comparison dropped -> :468, :469 (forest 1's
columns surviving a foreign install) and :508, :509, :511-518 (the
self-restore arm) move; F8(c) setForestBasis's re-imposition removed ->
:484, :485, :487-494 move; F8(d) adoptInstalledAmplitudePriors made a
no-op -> :475 moves, the donor's 0.125 read back as the recipient's
0.5. F6 in two arms, since its two halves fail differently: (i)
DBARTS_HAS_FIELD's >= weakened to > inside the calibration FILL ->
test-capi :963-965, :1001 and the flat-C BCF leg at :1123-1125 move;
(ii) the size guard dropped from the five appended FILLs -> the
poisoned pointers are dereferenced and R aborts (SIGSEGV, exit 139) at
the omitting-caller leg. No arm stayed green.

Weak-gate audit: tests/cpp's empty-calibration field enumeration
(test_sampler.cpp) was extended with the five NaN checks, consumer.c's
calibration leg with the five buffers and the BCF map read, and
test-capi's partial-caller comparison with all eight originals. NOT
extended, with the reason: test-calibration-midchain's
$setCalibration-isolation comparison enumerates four columns on a
SINGLE-FOREST sampler, where all five new columns are always NaN, so
adding them would be a silent no-op; the NaN is pinned directly at the
shape block instead. No R-level helper compares whole calibration
matrices.

Budget, both readings stated because the arc's two data points and its
written rule disagree. RAW NET added lines: engine 97, bridge 11, flat
C 39, R 2, docs 105, tests 413. Under the standing WRITTEN convention
(raw ~ 2x dense, multiplier-combiner.md "Budget units"): engine 49
(band 20-34, stop 48), bridge 6, flat C 20, R 1, docs 53, tests 207
(stop 215) - only the engine is over its stop, at 1.02x. Under the
ratio the S0 and M4.4 notes actually practiced (dense ~ 0.8 x raw):
engine 78, tests 330, both over their stops. The engine overrun is
Doxygen: the added block ran at 1.24 comment lines per code line
against chain.hpp's own 0.51, and was trimmed to 0.87 before landing;
the 55 code lines it documents are irreducible given five fields, three
vectors, a flag and two install-site helpers. The tests figure is the
mandated falsifier set - F4, F5, F6 and F8's four arms - which the
100-160 band predates in its costing.

## Landing note, S2 (appended 2026-08-15)

LANDED as e623fbf3. Both defaults, M4.4's documentation debt in cap-not-pin form,
and Fork 6c's K = 1 fix. R: defaultAmplitudePriorScale(family) beside
defaultNodeScale, forestParams(specs, hasBasis, family) taking
sqrt(2/length(specs)) on the withBasis branch and the family helper on
the other, R/spec.R passing family at its one call site,
bartcoreBCFSampler's sd.control defaulting NULL with the model's family
resolved ABOVE bcfParams (m3's two traps both honored). Engine:
BCFSpec::aPriorScale stays 2.0 and says why in its comment - a fixture
default no consumer reaches, and expandForestSpecs stays a pure
re-spelling. No C twin for the new helper, with the reason in its
comment.

ONE DEVIATION, and it is a correction to Fork 6's trace rather than a
choice. The plan says the data route's K = 1 error is applyBCFSpec's
"at least two per-forest parameter vectors"; re-traced by running it,
`dbartsData(x, y, bases = list(b))` never reaches the bridge - the
dbartsData S4 VALIDITY function refuses first with "'bases' must be null
or name at least two forests" (R/A_class.R), which validateForestBases
(the symbol the plan opened) does not carry. That floor also swallows
the forests route once forestBasisDeclarations' floor drops, since
dbarts() writes the declaration through dbartsData(bases = ), so the
designed refusal would have been reachable from dbartsSpec ALONE - which
assigns data@bases directly and so bypasses validity - and gate 11's
"both creation routes" would have failed. The floor therefore relaxes to
one forest (a container may carry one basis; the refusal is a modelling
statement, not a well-formedness one), and both routes now reach the
single designed refusal in resolveSamplerSpec. Mutation arm H pins that
this is load-bearing rather than tidying. Nothing pinned the old message.
The refusal names the RESOLVED count and its source, per R2-m4: "the
data object carries 1" where the count rode the data object, "this
call's 'basis' declarations resolve to 1" where a declaration replaced
it - so the length-1-over-a-two-basis-object case is not told it wrote
one basis. The K = 1 ticket is named in the code comment rather than in
the message; an R error naming an internal ticket would be worse
vocabulary than the two-forest spelling it points at instead.

One design detail the plan does not state and m3's trap requires:
defaultAmplitudePriorScale is TOTAL over the package's six-family
vocabulary (gaussian and aft 2, the four latent-scale families 1),
mirroring defaultNodeScale's own shape, rather than defined on the three
the multi-forest path builds. It still ERRORS on an unknown family, as
specced; but a helper defined on three would have fired BEFORE the
bridge on bartcoreBCFSampler(host, z, family = "aft") and replaced
"BCF does not support an AFT" - a deliberately placed backstop pin
(test-bcf-family.R) - with a worse message from the wrong layer.

DRAW LAW, verified by building both trees and comparing draws rather
than by argument. BITWISE across S2: gaussian K = 2 through the public
forests route (F2), the same through bartcoreBCFSampler (F2, internal),
and latent fits DECLARING forest(sd = 2) basis-free plus forest(sd = 1)
per basis forest at K = 2 probit, K = 3 probit and K = 3 logistic (F3 -
the composability claim, that the law is a default and a caller who
states the old values gets the old model). MOVED, as the two named
classes: probit K = 2 default, logistic K = 2 default, gaussian K = 3
default. Nothing else.

Battery on the slice's own private lib: preclean install clean;
tests/cpp all ok from make clean (242 ok lines, unchanged from S1 - the
cpp fixtures declare their scales explicitly and none moved) and again
under ASAN+UBSAN (detect_container_overflow=0), no diagnostic; full
tinytest 4872 results 0 failures; equivalence trio BITWISE on all three
baselines (equivalence-8b047f8b 37 compared / 0 skipped under
--strict-coverage, every scenario "identical draws (same RNG stream)";
bcf-equivalence-8b047f8b 12/12 "identical (all 7|8 channels)";
multinomial-equivalence-1027be5 10/10 identical - no max-|z| line
anywhere, and no baseline re-recorded in this arc); bcf-exact.R quick
OK as a confirmation, not a re-validation; air format --check clean;
lintr clean on all seven touched R files; R CMD check --as-cran from a
clean staged tarball outside the tree Status OK 0/0/0; NEWS.Rd parses
(229 entries); pkgdown::check_pkgdown no problems, no new Rd topic.

ASSERTIONS MOVED, all three the plan predicts and no others - the full
suite went from 4826 results to 4872 with exactly these red:
(1) test-bcf-family.R cell (c), the induced-index law, probit arm:
1.04912 * sqrt(K) * s -> 1.04912 * sqrt(2) * s, a CONSTANT, and the loop
gains K = 4 ("What S2 adds on top", first bullet); the value was
1.81712914 and is 1.48367953. (2) the same cell's logistic arm,
3.29590768 -> 2.69109698. (3) test-forest-basis-r5.R, F8(d)'s donor arm:
mapColumn(recipient, 1, "amplitude.prior.scale") 2 -> 1, that fixture
being probit with an undeclared basis-free forest. Values RE-DERIVED
from the law rather than regenerated from a run. Every other latent
K = 2 fixture either declares its sd or asserts an identity.

ADDED, per "What S2 adds on top": the K-invariance loop at K = 2, 3, 4
(a fourth unit-norm basis joins); the default-transport pins in their
own slots, with the MANDATORY gaussian arm and m1's reason in its
comment (slot 7 is 2 there and 1 under both latent families, and a
latent-only fixture cannot discriminate positions 4-8, which are then
all the literal 1); the K-aware slot-4 pins on BOTH shapes at K = 2, 3,
4, the shipped shape's being what refutes Fork 2d - K - 1 basis forests
take sqrt(2/K), not sqrt(2/(K-1)); the declared-sd override at K = 3 on
both shapes; the internal-route agreement, prior.scale AND amplitude
prior scale, plus its gaussian arm at 2; and the K = 1 refusal by
message on all three surfaces with the count's source discriminated.
The scaledBases cell keeps its LOAD-BEARING comment untouched.

Mutation proof, fresh mutated --preclean install per arm into a separate
library, the tree byte-verified restored after each:
(A) defaultAmplitudePriorScale not family-aware, 2 everywhere -> 6 moved
(the shipped-shape reader pin under both families, both latent slot-7
transport pins, the internal-route pin, and r5's donor arm);
(B) the K law dropped, factor default back to the literal 1 -> 16 moved
(both families' K-invariance loops, and every slot-4 pin);
(C) the law applied to the BASIS-FREE forest too -> 7 moved (the slot-4
pins' factors[1] == 1 arm, and the declared-shape pin);
(D) normalized on the BASIS count Kb rather than K, i.e. Fork 2d -> 15
moved across three files, including test-bcf-creation's K = 2 transport
pins and r5's staleness cell;
(E) forestBasisDeclarations' floor back at < 2L, restoring the silent
drop -> 4 moved (the forests-route refusal, the declaration-over-data
case, and the basis+sd message change);
(F) the designed refusal in R/spec.R removed -> 6 moved, i.e. every
K = 1 assertion on all three surfaces;
(G) bartcoreBCFSampler back to the literal sd.control = 2 -> 2 moved
(the internal-route agreement, both halves);
(H) the dbartsData validity floor left at two forests -> 4 moved, which
is what makes the deviation above load-bearing rather than tidying.
No arm stayed green.

Budget, both readings, as S1's note established. RAW NET added lines:
R 115, engine 9 (a comment block and one declaration), docs 124, tests
204. Under the standing WRITTEN convention (dense ~ raw/2): R 58 (band
16-32, stop 46), docs 62 (60-110), tests 102 (75-135) - only R over.
Counting added NON-COMMENT non-blank lines instead: R 46 (exactly the
stop), docs 153, tests 139. R is the leg that runs long either way, and
the cause is enumerable: 18 of the 46 are the designed refusal's own
message body as air formats it, and 12 are the family helper's switch -
both mandated shapes, and the "four lines" the agent line prices was the
defaults alone. feature-matrix.md: CHECKED, no cell moves - a default is
not a capability, and no footnote enumerates one (the S1 pass on [f16]
covers what $getCalibration reports).
