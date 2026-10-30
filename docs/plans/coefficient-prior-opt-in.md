# coefficient-prior-opt-in: a heavy-tailed prior on a forest's drawn coefficient

Status: LANDED 2026-10-30 (dec-A245), built in two slices with one review and a text-only fix round; no call
was held (dec-B509 to dec-B511 ruled the open ones). Rulings: dec-B262 (the opt-in, its user and its
measurement), dec-B271 (one argument, `amplitude`, beside `fixed()`; `student(3)` the form shown), dec-B509 (a
forest with no basis held by default, bcf's prior the opt-in `student(1)`), dec-B510 (`student()` is df 1, a df
finite and at least 1) and dec-B511 (the live change its own item), with dec-B18, dec-B253, dec-B254, dec-B261,
dec-B266, dec-B290 to dec-B292, dec-B354, dec-B503, dec-B505 and dec-B508.

agent: one opus implementer for both slices; one opus review over both that runs the mutants below, the same
reviewer checking the fix round unless the fixes change engine behavior. Tier "Changes draws" (a prior and the
engine's coefficient draw; [Process by risk](README.md#process-by-risk)). The blind critique that tier asks for
has run.
rng: NEUTRAL, bit for bit, for every call main accepts, by structure: a forest that states no df takes today's
lines, verbatim, behind a branch on whether a df is stated. A forest with no basis stating `student(1)` at the
median its default has is main's drawn default restated and draws what it draws, bit for bit. Any other
`amplitude = student(df)` is a new prior, posterior-defining for those models: gated by new exact-posterior arms
and a design-note amendment. No re-record.
window: after the merge to main, before 1.0-0 (dec-B262: "Sure, opt-in for 1.0. Post main-merge."). Additive:
`amplitude = student(3)` is "could not find function" today. It lands before forests-shared-total's slice F, the
default flip of dec-B509, and before bartCause's `bcf()` states `student(1)` ([Sequencing](#sequencing)).
budget: ~1130 lines in the tree. Slice 1 ~330 (engine ~85, support library ~3, bridge ~35, tests/cpp ~205).
Slice 2 ~740 (R ~140, help ~70, NEWS ~5, tinytest ~280, exact arms in forest-prior-exact.R ~190, design notes
~55). Records ~40. Out of the tree, ~320 lines of scratch A/B script.

## Summary for the maintainer

You ruled that a caller may give a forest's drawn coefficient a heavy-tailed prior, written
`forest(..., amplitude = student(3))` beside `fixed()`, with the default unchanged. This plan builds it.

With `student(df)` the coefficient of a forest with a basis is a Student-t with `df` degrees of freedom instead
of a normal, scaled so that its median size is exactly the normal's. The forest's `sd` keeps its meaning, the
prior median of the forest's contribution, and only the tails move. A coefficient ten times its median size has
prior probability 6 percent under `student(1)` (the half-Cauchy of the measurement you were shown), 0.5 percent
under `student(3)`, and less than one in ten billion under the normal. That tail is what lets the coefficient
follow a multiplier whose scale drifts during a run. The sampler draws it exactly, by the same two-step update
the forest with no basis already uses; nothing is approximated. Default fits do not move, bit for bit. How well
a heavy-tailed coefficient mixes in this sampler is not measured; your ruling's numbers came from a prototype.

A note on wording: the earlier ruling's "nothing stated is the normal prior" describes a forest with a basis,
the multiplied forest your opt-in ruling is about. The coefficient of a forest with no basis has been
half-Cauchy since the causal model was matched to bcf, which is `student(1)` in this spelling; probed today, its
draws match that, and a numeric column's coefficient a standard normal.

Your three rulings on this plan's questions are built in. `student()` means `student(1)`, the half-Cauchy, and a
stated df must be finite and at least 1 (dec-B510). It is taken on a forest with no basis, where it turns on a
drawn coefficient: by your ruling that such a forest is held by default (dec-B509), it is how bcf's prognostic
prior is written, and bartCause's `bcf()` will state it. Changing a coefficient's prior on a live sampler is its
own item, after this one (dec-B511). This plan lands before the flip to a held default, so that bartCause can
state its prior first; until the flip, `student(1)` on a forest with no basis restates today's default exactly.

The size is about 1130 lines with tests and gates; the ruling's 300 counted the code alone. About 14 agent hours,
and plans here run 1.5 to 2 times low.

## Goal

`forest(amplitude = student(df))` gives a multiplied forest's drawn coefficient a Student-t prior with `df`
degrees of freedom whose median size is the normal default's, and gives a forest with no basis a drawn coefficient
under a t whose median is its sd (dec-B509), at every door that takes several forests; the
reader, `print`, `summary` and a fit's `leaf.prior` say so; a copy, a reload and a state install carry it; and
every call main accepts draws what it drew.

## What exists today

Run on main at 177fd104 (shipped build, arm64 macOS, private library scratch/libs/cpo-plan; probes
scratch/cpo-plan/p1.R to p4.R, untracked). Data: 100 to 200 rows, gaussian `y`, `w ~ N(0, 3)` or `dose ~ N(5, 2)`,
a two-level `z`. The critic reproduced every premise below on its own install (scratch/crit-cpo).

- No `student()` reaches `amplitude`. In a formula term and in a `forests = list()` at `dbarts()`:
  `could not find function "student"` (p1, p4). `amplitude = "student"`: "a forest's 'amplitude' must be fixed(),
  which holds its coefficient, or left out, which draws it" (p4;
  [`validateForestAmplitude`](../../R/model.R)). `names(dbartsPriors)` is cgm, dart, normal, linear, gp, chisq,
  fixed, chi, invchi; `dbartsFamilies` has `student`, the response family (p1;
  [`student`](../../R/family.R), `df = NULL` meaning "estimate"). `family = student(3)` fits a Student-t
  response (p3).
- The name is reserved where a prior is not stated: `basis = ~ student(dose)` is refused, "a prior is not stated
  on a basis term. Give the forest's size as 'sd' and its coefficient's prior as 'amplitude'" (p3;
  [`BASIS_RESERVED_CALLS`](../../R/forestBasis.R)).
- What the engine draws ([`forestMapEntry`](../../src/bartcore/combiner.hpp),
  [`ForestAmplitudePrior`](../../src/bartcore/combiner.hpp)). A drawn forest with a basis: N(0, 1) on a numeric
  column's coefficient, N(0, 1/2) on each of a factor's two, a fixed variance; its leaf scale is its sd over
  0.674. A drawn forest with no basis: N(0, v) with v a live inverse-gamma auxiliary refreshed after each draw
  (`drawAmplitudes`, IG((1 + q)/2, (c^2 + |a|^2)/2)), marginally Cauchy(0, c), c the sd as a multiple of the
  unit; its leaf scale is the unit. Measured under a flat likelihood (sigma held at 1e6, one tree a forest,
  60000 sweeps, p2): the no-basis coefficient's |a| median 2.02 and 90th percentile 12.96 against a Cauchy(0, 2)'s
  2 and 12.63, a Kolmogorov-Smirnov p of 0.79; the numeric coefficient 0.673 and 1.645 against N(0, 1)'s 0.674
  and 1.645, p 0.32; the factor's difference 0.675 and 1.646, the two coefficients uncorrelated (0.00), sd 0.707
  each. The recent critique's "normal with an inverse-gamma variance, marginally Cauchy" is the forest with no
  basis alone.
- The state already carries one variance a forest and a restore installs it where the forest's prior is a
  mixture ([`AmplitudeForestCombiner`](../../src/bartcore/combiner.hpp)'s `restoreGlue`); the shipped C header
  has no forest-creation entry ([`DBARTS_C_API_HASH`](../../inst/include/dbarts/dbarts.h)).
- The record. The control's `bartcore.forests` holds `params` (four numbers a forest, the fourth 0 where held),
  `sd` (NA where none stated) and `heldAmplitude` (p4: `75 0.95 2 1`, `50 0.25 3 0`, `NA NA`, `FALSE TRUE`);
  a `copy()` keeps the hold (p4). The bridge reads the statements in
  [`applyForestStatements`](../../src/R_interface_bartcore.cpp) through
  [`forestNumbers`](../../src/R_interface_bartcore.cpp), which returns nothing for an absent entry.
- The reader ([`reportLeafPrior`](../../R/dbarts.R)): `amplitude` is `fixed()` on a held forest and absent on a
  drawn one; `prior.sd.of` is "forest (prior median)" for a drawn forest with no basis, "forest, per unit of basis
  (prior median)" for a drawn numeric one, "level difference" for a held factor; `k.scale` 1.93 (the unit) on the
  first, 0.893 on the numeric (its sd 0.602 over 0.674) (p1). The bridge's raw row carries
  `amplitude.prior.variance` and `amplitude.prior.scale`, which no R code reads.
- Print and summary ([`printedForest`](../../R/forestPrint.R)): the short block has no word on the coefficient;
  the full block ends a held forest's line "; coefficient held", decided by `!is.null(entry[["amplitude"]])`
  (p3: "3 z: factor basis, 2 levels; sd 1.553 between the two levels (default); coefficient held"). A fit carries
  the reader's entries in `fit$leaf.prior` and the hold in `fit$fixed$amplitude` (p3).
- A live sampler: `$setLeafPrior(forests = list(forest(), forest(amplitude = fixed())))` is "$setLeafPrior
  restates only a forest's 'sd': 'amplitude' is fixed at creation" (p3;
  [`resolveForestSpreads`](../../R/dbarts.R)).
- The family door resolves families and priors in one vocabulary, `c(dbartsFamilies, dbartsPriors)`, whose names
  its comment calls disjoint ([`resolveFamily`](../../R/family.R)). Where two entries share a name the FIRST is
  bound: `vocabularyEnv` assigns `vocabulary[[name]]` (the critic's probe).
- Other prior sites refuse a foreign prior by name: `leaf.prior = fixed(1)` is "'leaf.prior' must be a leaf prior
  specification; see ?dbartsPriors", and the same for `tree.prior`.
- [`forestMapEntry`](../../src/bartcore/combiner.hpp) lives in combiner.hpp, which does not include the support
  library's Rmath wrapper; the engine's model header does, for `Rf_dt`.
- Consumers. bartCause's `bcf()` passes `amplitude = NULL` or `fixed()` on its two forests, and states the
  prognostic forest's sd only where it holds it under gaussian; under dec-B509 it states `student(1)` and its sd
  ([Sequencing](#sequencing)). stan4bart passes one `forest(n.trees = )` and no amplitude and does not change.
- Units (run on main's code at 1ffd6f92, 2026-10-29): a forest with no basis beside a factor forest stating `sd`
  equal to its default draws `identical()` to the default under probit (`sd = 1`) and logistic (`sd = pi /
  sqrt(3)`); under gaussian R's `2 * sd(y)` and the engine's 2U differ in the last place (1.2166873523600565
  against 1.2166873523600563) and the draws differ. `qt(0.75, 1)` is exactly 1.

## What 0.9-34 did

Nothing here. Probed in scratch/libs/cran: no `forest`, no `student` in the namespace, no `forests` argument and
no multiplied forest; its `cauchy` is `rbart_vi`'s random-effect prior. Nothing to mirror.

## Design

The plan is written for dec-B510 (df 1 by default, finite and at least 1), dec-B509 (taken on a forest with no
basis) and dec-B511 (at creation only).

### The prior, exactly

A multiplied forest's coefficient block a has q coordinates: 1 for a numeric column, 2 for a factor. The default
puts a ~ N(0, kappa I), kappa 1 for a numeric column and 1/2 for a factor's levels, and gives the forest's own fit
the leaf scale sd / 0.674.

`amplitude = student(nu)`:

- Prior: a | v ~ N(0, v I_q), v ~ IG(nu / 2, nu c^2 / 2), so a is multivariate t with nu degrees of freedom and
  scale c, one v for the block.
- The median match: c = sqrt(kappa) z_.75 / t_nu,.75, z_.75 = 0.6744897501960817 the normal's 0.75 quantile and
  t_nu,.75 the t's. A numeric coefficient's |a| then has median sqrt(kappa) z_.75, the normal's to the last
  place; a factor's difference b1 - b0 is t_nu with scale sqrt(2) c and median z_.75, the default's
  N(0, 1) difference's. At nu = 1, c = 0.6745 (numeric); at nu = 3, 0.8818.
- The forest's leaf scale and `sd` do not move: the node scale stays sd / 0.674, so the contribution's median is
  the default's exactly and only the tails differ. A coefficient beyond ten times its median has prior mass
  0.063 at nu = 1, 0.0046 at nu = 3, 1.5e-11 under the normal; the 0.90 over 0.50 quantile ratio of its size is
  6.31, 3.08 and 2.44.
- The update, a two-block Gibbs scan with no Metropolis step and acceptance 1: a | v and the data, the Gaussian
  full conditional `drawForestAmplitude` draws today, reading the live v as its prior variance as it already does
  for the forest with no basis; then v | a ~ IG((nu + q) / 2, (nu c^2 + |a|^2) / 2), conjugate, reading no data.
  One gamma draw a stated forest a sweep. No other move reads the coefficient's prior.
- v starts at kappa, the variance the map already gives the block, so the first draw of a is the default's.
- nu is 1 where none is given and otherwise finite and at least 1 (dec-B510): below about 0.01 c underflows
  (t_nu,.75 is 6e28 at 0.01 and infinite at 1e-6). Infinite nu, the normal, is refused by the same rule.

### On a forest with no basis (dec-B509)

Today such a forest beside a forest with a basis draws a scalar coefficient a | v ~ N(0, v), v the live
inverse-gamma auxiliary, marginally Cauchy(0, c) with c its median, its leaf scale at the unit and its size
carried by a. `amplitude = student(nu)` keeps that arrangement and takes nu: v ~ IG(nu / 2, nu c^2 / 2),
c = m / t_nu,.75 with m the median. At nu = 1, t_1,.75 is exactly 1 and the draw takes today's unstated lines
(the map's `mixtureDf` left NaN, the statement's df still reported), so `student(1)` is today's prior operand for
operand. v starts at 1, as today.

- The median m: the stated `sd`, or the prior median of the size the forest would otherwise have. On main until
  forests-shared-total's slice F, that is the default's own median, 2U under gaussian and U under probit and
  logistic, so `student(1)` with no sd restates the default and draws `identical()` to it. After slice F the
  forest is held by default and m is the prior median of its held size: `k.scale` / k at a held k, `k.scale`
  over the chi prior's median k at a drawn one, a named sd as given (forests-shared-total, Slice F engine item 3:
  1.574 under probit at `bart()`'s chi(1.5, 2)).
- A drawn k: a forest stating `student()` has no k of its own; its size is in the coefficient, and a `chi()` at
  the fitting function enters only through its median. The probit rescaling step does not apply to it (it needs a
  drawn k), as it does not to today's drawn default; under logistic nothing changes.
- Where it is taken: in a model of several forests, on a forest with no basis that states sd or that shares the
  total with no other (forests-shared-total D3, which refuses it among several sharers). On a lone forest it is
  refused with `fixed()`'s refusal, reworded (R item 5): the single-forest chain has no coefficient.

Why that is the whole proof of correctness: the first block is today's draw, unchanged and already gated at a
fixed variance by benchmarks/R/bcf-exact.R's mode 2b and the numeric arms of benchmarks/R/forest-prior-exact.R;
the second reads no data, so the flat-likelihood prior arm below gates it exactly, and it tells one shared v from
one a coordinate (P(|b1 - b0| <= median) is 0.456 under the latter against 0.5). Data arms with the t in place
add the joint check. A posterior functional of bcf-exact.R's 2b design cannot tell the normal, this t, independent
t's or a Cauchy apart (the critic's quadrature: E[(b1 - b0) tau] equal to four places, E[mu] within 0.013 of a
0.07 tolerance), so that script gains no mode.

### Engine (slice 1)

R-agnostic: no R name, record or object. `--preclean` on every commit.

1. [`ForestStatement`](../../src/bartcore/combiner.hpp) gains `amplitudeDf`; NaN means none is stated. It rides
   into `ForestSpec` and the chain's map record, so every restatement
   ([`Chain::setForestMapSd`](../../src/bartcore/chain.hpp)) re-resolves with it.
2. [`ForestAmplitudePrior`](../../src/bartcore/combiner.hpp): `halfCauchyScale` becomes `mixtureScale` (a
   positive one making the variance the live auxiliary, as now) and gains `mixtureDf`, NaN by default, the
   stated nu. The rename touches about twelve sites and tests/cpp; no expression changes.
3. [`forestMapEntry`](../../src/bartcore/combiner.hpp): a drawn forest of levels or numbers that states a df keeps
   its node scale, its variance kappa, and takes `mixtureDf` nu and `mixtureScale` c from the median match; a
   drawn forest with no basis that states a df keeps today's arrangement and takes nu and c = m / t_nu,.75, nu = 1
   leaving `mixtureDf` NaN ([On a forest with no basis (dec-B509)](#on-a-forest-with-no-basis-dec-b509)); a held
   forest ignores a df. Where none is stated every output is today's, operand for operand. A stated df below 1 or
   not finite throws `std::invalid_argument` (the backstop; R refuses each first). t_nu,.75 comes from the support library: `ext_quantileOfStudentT`, one
   line in src/include/external/stats.h over Rmath's `Rf_qt`, beside
   [`ext_quantileOfNormal`](../../src/include/external/stats.h), and combiner.hpp includes that header.
4. `drawAmplitudes`: where `mixtureDf` is NaN, today's two lines of shape and rate, verbatim; where stated, the
   nu expressions. Neutrality is then structure, not an argument about `0.5 * 1.0`.
5. [`ForestCalibration`](../../src/bartcore/chain.hpp) gains `amplitudeDf`, the stated df, NaN where none;
   `amplitudePriorScale` reports the mixture scale on any mixing forest.
6. No facade virtual, no state block, no `stateFormatVersion` move: the state's per-forest variance is v on a
   mixing forest, installed where the recipient mixes, as now.

### Bridge (slice 1)

- [`applyForestStatements`](../../src/R_interface_bartcore.cpp) reads `amplitude.df`, one number or NA a
  forest, through [`forestNumbers`](../../src/R_interface_bartcore.cpp); NA and an absent entry state nothing.
- [`bartcore_getLeafPrior`](../../src/R_interface_bartcore.cpp) gains the column `amplitude.df`.
- No registration's arity changes; the shipped header does not move.

### R (slice 2)

1. A constructor `student(df = 1)`, exported only as `dbartsPriors$student` (an internal name of its own, the
   family's being `student`), class `dbartsStudentPrior`, one slot, valid for a single finite number at least 1
   (dec-B510). Refused by name otherwise: "a coefficient's student(df) takes one finite df of 1 or more;
   student() is student(1), the half-Cauchy".
2. Vocabularies. `student` joins `fixed` in [`forestConstructors`](../../R/model.R), the `forests` entry of
   [`FOREST_ARGUMENT_VOCABULARIES`](../../R/model.R), [`processHit`](../../R/formulaTerms.R)'s vocabulary and the
   hint list of `forceCallerCode`. [`resolveFamily`](../../R/family.R) keeps the families first and drops the
   duplicate explicitly, `c(dbartsFamilies, dbartsPriors[setdiff(names(dbartsPriors), names(dbartsFamilies))])`,
   so `student` there is the family whatever the order; its comment says the two lists share that name.
3. [`validateForestAmplitude`](../../R/model.R) takes `fixed()` or a `dbartsStudentPrior` on any forest of a
   model of several; a bare `student`, as a bare `fixed` is `fixed()`, is `student()`. Refused by name: a family
   object ("a forest's 'amplitude' is fixed() or student(df); dbartsFamilies$student is the response's family")
   and anything else, naming both; `student(Inf)` by item 1's rule.
4. [`forestStatements`](../../R/model.R) adds `amplitude.df` only where some forest states one, so a default
   fit's record and every fit object made without it are as today. `params` and `heldAmplitude` do not change.
5. The one-forest refusal reads "which has no coefficient to give a prior", covering both values. (Whether
   `fixed()`, the default restated after the flip, is taken on a lone forest is forests-shared-total's open call
   4; `student()` stays refused there.)
6. [`reportLeafPrior`](../../R/dbarts.R): `amplitude` is `student(df)` where stated, `fixed()` where held, absent
   otherwise; a fit's `leaf.prior` follows. `prior.sd.of` and `k.scale` do not change.
7. [`printedForest`](../../R/forestPrint.R): the full block's test becomes "is `fixed()`", and a stated forest's
   line ends "; coefficient student(3)". The short block stays silent, as for a hold.
8. `$setLeafPrior(forests = )` keeps refusing `amplitude` as fixed at creation until coefficient-prior-live is
   built (dec-B511). Nothing here should make that switch harder: v starts at the normal's variance kappa on a
   multiplied forest (Calls made 4), where dec-B511 starts a t's variance on a live switch, and the df rides the
   statement and the map record, which a restatement re-resolves.

### Help, NEWS and docs

- man/forest.Rd `amplitude`: the default by kind (normal with a basis; half-Cauchy with none until
  forests-shared-total's slice F, held after), `student(df)` and the median rule on each kind, bcf's prognostic
  prior written `amplitude = student(1), sd = 2 * sd(y)`, what a tail buys and when; man/dbartsPriors.Rd: a
  `student(df)` item, saying that a bare `student()` means df = 1 here and an estimated df as a response family
  (dec-B510); man/dbartsSampler-class.Rd: `getLeafPrior`'s
  `amplitude` entry.
- NEWS, in 1.0-0's forest item: one clause.
- Amended, not new notes: [1. The rule](../design/forest-multiplier.md#1-the-rule) gains the drawn
  `student(df)` row and section 7's heavy-tail bullet goes;
  [The calibration map, general in K](../design/multiplier-combiner.md#the-calibration-map-general-in-k) gains
  the mixture's nu and the median match.

## Interaction with the three pending plans

forest-after-merge (8adc5d1a): dec-B508 removed its long form, `forest(leaf.prior = )`, and its first slice, so
`normal` is not added to the forest vocabularies and `forest(leaf.prior = student(3))` is R's unused-argument
error with no refusal of its own. Its remaining slice's `restateForestBasisScale` re-resolves
[`forestMapEntry`](../../src/bartcore/combiner.hpp) with "the forest's own default statement": it must keep the
statement's `amplitudeDf` and default only the sd. Under Calls made 1 the coefficient's prior does not depend on
the column's scale, so its "no amplitude prior moves" then holds for a `student()` forest. Whichever lands second
pins: a `student(3)` numeric forest swapped with `updateBasisScale = TRUE` keeps the reader's `amplitude`, and its
draws equal the creation twin's (a restatement that drops the df fails it). No shared facade edit.

forests-shared-total (its branch, amended with this one for dec-B509): this plan lands first. That plan carries the default flip
(its slice F), the median of a `student()` forest with no basis and no sd after the flip (its slice F engine item
3), and the refusal of `student()` on a sharer among several (its D3: no shared coefficient is built). Shared code:
[`forestMapEntry`](../../src/bartcore/combiner.hpp) (its share and the kind-none hold; this item's nu and c),
[`validateForestAmplitude`](../../R/model.R) and benchmarks/R/forest-prior-exact.R; it rebases onto this item.
This item's pins of `student(1)` against the pre-flip default are re-pinned by its slice F.

setmodel-parameters-not-structure (8342bc31):

- Slices 1 to 3 refuse or never reach a sampler with amplitudes; the coefficient's prior rides the forests'
  record on the control, which no model comparison reads. Its family checks meet this item's vocabulary change:
  whichever lands second pins `$setModel` with a model of `family = student(5)` resolved as the family.
- Slice 4 opens `$setModel` on such samplers for the tree and residual priors and never routes through the
  leaf-scale write. Whichever lands second pins: on a `student(3)` forest, `$setModel` with a new tree prior
  leaves the reader's `amplitude`, the coefficient and v as they were, and the next draws equal a copy's. Its
  `$setTreePrior(forests = )`, evaluated in the widened vocabulary, refuses `forest(amplitude = student(3))` as
  any knob but `base` and `power`.
- coefficient-prior-live (dec-B511), after its slice 4 and this item, adds the live change through
  `$setForest(amplitude = )` (dec-B518); the coefficient's prior is not the model's, so `$setModel` never carries it.

### Sequencing

The default flip must not land before `student(1)` on a forest with no basis exists, or bartCause's `bcf()` loses
its prior. Order of landing:

1. This item, both slices, on dbarts's main: `student(1)` on a forest with no basis is main's drawn default
   restated.
2. bartCause's `bcf()`, on its own branch and main: where it draws the prognostic coefficient it states
   `amplitude = student(1)` and `sd` 2 sd(y) under gaussian or the link's unit under probit and logistic (1, pi /
   sqrt(3)); where it holds it, the link's unit under probit and logistic too (the held binary default moves to
   the drawn k at the flip). Its draws stay bit for bit under probit and logistic and differ in rounding under
   gaussian (the unit's last place); its prior does not move. About 20 lines and a test, ~1.5 agent-hours.
3. forests-shared-total, slices B, 1, 2 and F, merged once, its slice F built after a rebase onto 1 and merged
   after 2.

## Slices

### Slice 1: engine and bridge

Steps: engine items 1 to 6 and the bridge; tests/cpp below. Nothing in R reaches the new branch.

RNG: neutral, every draw of every call. Gates: tests/cpp plain and under address and undefined sanitizers; the
equivalence trio bitwise; the full tinytest suite.

Cost: ~330 lines; 2.5 agent hours to build.

### Slice 2: R surface, help and the exact arms

Steps: R items 1 to 8, help, NEWS, the design-note amendments, tinytest, the exact arms.

RNG: neutral for every call main accepts; posterior-defining for a forest stating `student()`. Gates: the whole
list below on the slice's tip.

Cost: ~740 lines in the tree and ~320 of scratch A/B; 4.5 agent hours to build.

Both: one review over both slices, 2 hours; the gate battery 3 hours, most of it compute; the fix round 1.5;
records 0.5. About 14 agent hours; plans here run 1.5 to 2 times low, so budget 21 to 28. The x86-64 leg before
landing is the coordinator's (GCC contracts multiply-adds by default).

## Tests and gates

tests/cpp (slice 1):

- The map: no df stated gives today's entry bitwise for every kind, held and drawn; `student(1)` and `student(3)`
  on a numeric and a factor forest leave the node scale and kappa bitwise and give c against literals
  (0.6744897501960817, 0.8818100603, and sqrt(1/2) times each); a held forest ignores a df; NaN states none; a df
  of 0, a negative, 0.5 and Inf throw. A drawn forest with no basis: df 1 bitwise today's entry, `mixtureDf` NaN
  and the reported df 1; df 3 keeps the node scale at the unit and gives c = m / t_3,.75.
- The draw: a K = 3 fixture stating no df draws bitwise what it drew before the slice (a recorded vector); with a
  student forest, v after each sweep against IG((nu + q) / 2, (nu c^2 + |a|^2) / 2) read through `TestPeer` on a
  seeded generator; v starts at kappa; at nu = 1, 10000 sweeps leave every a and v finite; a forest with no
  basis at df 1 draws bitwise what the unstated K = 3 fixture drew.
- State: serialize and restore carry v on a student forest; a normal donor's variance installs as v; a student
  donor's v does not reach a normal recipient.
- The calibration reports the stated df and the mixture scale; NaN where none is stated.

tinytest, a new test-forest-amplitude-student.R, small n:

- Every door: a formula term, `forests = list()`, a data object's `bases` with `forests`, and `dbartsSpec()`;
  gaussian, probit and logistic; a numeric and a factor basis.
- The reader's `amplitude` `identical()` to `student(3)`; `fit$leaf.prior`; the full block's line and the short
  block's silence; `fixed()` still printed as held; a default fit's forests record without `amplitude.df`.
- `student()` and a bare `student` read back as `student(1)` (dec-B510).
- Refusals by message, counted with `withCallingHandlers`: `student(0)`, `student(0.5)`, `student(NA)`,
  `student(NaN)`, `student(c(1, 3))`, `student(Inf)`, `dbartsFamilies$student(3)`, `amplitude` on one forest,
  `$setLeafPrior(forests = list(forest(amplitude = student(3))))` (dec-B511).
- A forest with no basis (dec-B509): beside a factor forest, `student(1)` with no sd, and with `sd = 1` under
  probit, draws `identical()` to the default; `student(3)` there runs, the reader's `amplitude` `student(3)` and
  its `sd` the default's median.
- `family = student(3)` still the response family, at `bart()`, `dbarts()` and `dbartsSpec()`;
  `leaf.prior = student(3)` refused as a foreign prior; `dbartsPriors$student(3)` outside the argument builds the
  prior.
- `copy()`, and `storeState`, `saveRDS`, `readRDS` and `setState` after a run: the next draws `identical()` to
  the live sampler's; a state installed from a default sampler into a student one runs.
- A coarse, seeded prior check in process: under a flat likelihood (sigma held at 1e6), 4000 sweeps, the numeric
  coefficient's median size within 0.1 of 0.674, and its 0.90 over 0.50 ratio above 2.8 under `student(3)` and
  below 2.7 under the default (the critic's 300 seeds: lowest 2.89, highest 2.57).

Exact posterior, new arms of benchmarks/R/forest-prior-exact.R in quick mode, run by exact-gates.yaml, which
already lists the script:

- A prior arm: sigma held at 1e6 so the likelihood is flat to 1e-9; one single-tree forest with no basis (held),
  one numeric and one factor, drawn; and a second model whose forest with no basis is drawn under `student(1)` and
  `student(3)`, P(|a| <= q_p) against the t's closed form at its median; P(|a| <= q_p) at p = 0.25, 0.5 and 0.9 against the closed form for
  `student(1)`, `student(3)` and the default, and the same for the factor's b1 - b0; the script's verdict, a gap
  past 0.0005 plus 4 standard errors over seeds.
- Numeric arms: the existing drawn numeric arm at `student(1)` and `student(3)`, its oracle integrating the
  coefficient on a midpoint grid of the prior's quantile scale, a = c qt(u, nu); functionals unchanged.
- A drawn-factor arm, at the default and at `student(3)`: the oracle a two-dimensional quadrature over (b0, b1),
  the cell values integrated in closed form given them; functionals the posterior mean per class and
  P(|b1 - b0| <= t) at t = 0.5, 1 and 2. It closes the script's "Not seen here: a drawn coefficient on a factor".
- The grid rule: every oracle's functionals agree to 1e-4 between its grid and one twice as wide (twice as fine
  in u on the quantile scale); a prior-mass cutoff is not the criterion (grids to 7 and to 14 agree to 1e-4).

Mutants, each to fail a named expectation:

- c without sqrt(kappa) on a factor (the prior arm's factor difference; the map literal); z_.75 replaced by 0.674
  (the map literal); nu dropped from the rate (the nu = 3 prior arm; the numeric arm, P(|a| < 0.5) 0.909 against
  0.831); the shape's (1 + q) kept at nu = 3 (the prior arm; the numeric arm, 0.759); v per coordinate on a factor
  (the prior arm's P(|b1 - b0| <= median), 0.456).
- The node scale moved with the prior (the map's bitwise node scale; the reader's `k.scale`).
- The unstated path rewritten in nu (the bitwise draw fixture; the bcf equivalence compare; the A/B).
- v not restored on a student forest (the tests/cpp state case; the tinytest continuation).
- `amplitude.df` not written to the record (the copy test); written on every fit (the default record test); not
  read by the bridge (the reader test).
- `student` missing from the formula term's vocabulary (the formula door); the family door taking the prior
  (`family = student(3)` pin).
- Print's held test left as `!is.null` (the full block's student line).
- On a forest with no basis, c not divided by t_nu,.75 at nu = 3 (the no-basis prior arm; the map literal); the
  multiplied forest's arrangement taken there, node scale over 0.674 (the bitwise df 1 pin).

Gates, on the branch's own libraries:

- `R CMD INSTALL --preclean`, plain and `--configure-args=--enable-reference-build`, `buildInfo()$mode` checked.
- tests/cpp plain and under address and undefined sanitizers; the R-loaded sanitizer run over the new test file
  (`Rf_qt` and the new update are reachable from R).
- A neutrality A/B (scratch, ~300 lines): a 177fd104 library against the branch's, on the reference build and on
  the shipped build, two and three chains, every recorded channel compared with `identical()`: bartCause's shape
  drawn and held under gaussian, probit and logistic, and its drawn shape again with `student(1)` stated on the
  forest with no basis, against the 177fd104 library's default (identical); three forests (none, a factor, a numeric column) with and
  without case weights and with zero weights; the data object's `bases` door and `dbartsSpec()`; a held forest
  with no basis; `$setLeafPrior(forests = )` restating the no-basis forest's sd (the renamed scale's setter) and a
  multiplied one's; `$setForestBasis`; a warm start; each continued after `copy()` and after `saveRDS` and
  `readRDS`. A difference is a stop.
- The equivalence trio on the reference build, `compare --bitwise --strict-coverage`: 55 of 55
  (equivalence-32418b6e), 19 of 19 (bcf-equivalence-20eed7ed) and 11 of 11 (multinomial-equivalence-80b1c8d4),
  identical, no `max |z|` line; the four snapshot files unchanged.
- exact-gates.yaml's list in quick mode, the new arms among them, every existing arm passing as before.
- The full tinytest suite in one process; `lintr::lint_package()`, `air format --check .`,
  `tools/check-rc-codoc.R`, `tools/check-win-drift.R`, `tools/check-doc-freshness.R`, each on its own exit
  status; `inst/NEWS.Rd` parses; `R CMD check --as-cran` from a tarball of a clean copy staged outside the tree.
- No `bench-sampler.R` compare: the change runs once a forest a sweep, never per node or row, and the default's
  lines are unchanged.
- Before landing, the x86-64 leg (the coordinator's): tests/cpp, the suite and the A/B there.

## Cost

- One more prior to explain: the default differs by kind (normal with a basis; half-Cauchy without until
  forests-shared-total's flip, held after), which the help now says outright.
- The heaviest tail is the default df (dec-B510's cost), likely the slowest to mix along the ridge with the
  leaves; dec-B509's check compared a held prognostic forest with it, not mixing within it.
- A heavy-tailed coefficient mixes along the ridge it shares with the forest's leaves more slowly than a normal
  one; the opt-in's user accepts that for a scale that drifts. Not measured here.
- `student` names two things, the response family at `family =` and the coefficient's prior at `amplitude =`,
  each resolved by where it is written, and a bare `student()` reads two ways: df 1 here, an estimated df there
  (dec-B510).

## Calls made

1. The median match keeps the node scale: the t replaces the normal on the coefficient alone, scaled to the
   normal's median, so the forest's `sd`, leaf scale and `k.scale` do not move, and a swap or a restated sd
   touches no coefficient prior. Alternative: the forest with no basis's arrangement, the leaf scale at the unit
   and the sd carried by the t's median; the same prior on the contribution to 0.07 percent, but switching the
   prior would move the leaf scale and every restatement would have to restate the mixture.
2. The median matched is the normal's own, z_.75, not the 0.674 the node scale divides by, so the default and the
   opt-in differ in the tails alone. Alternative: 0.674, which makes the contribution's median the stated sd
   exactly and the two priors' medians 0.07 percent apart.
3. One v for a factor's pair: the pair is bivariate t and the difference a t with the matched median.
   Alternative: one v a coordinate, whose difference has no closed-form median (0.456 of the mass inside the
   target) and a different tail. Low stakes; worth the maintainer's nod.
4. v starts at kappa, the variance the map gives the block today, and where dec-B511 starts a t's variance on a
   live switch. Alternative: the prior median of v, which starts a stated forest off the default's first draw.
5. The constructor lives in `dbartsPriors`, the home of `fixed()`, under an internal name of its own; the family
   door keeps the families first and drops the duplicate name explicitly. Alternative: taking the family object
   `student(3)` as a coefficient prior, a response distribution whose missing df means "estimate it".
6. RULED by dec-B510: `student(Inf)` and any df below 1 are refused by name.
7. The df rides the forests' record as `amplitude.df`, written only where some forest states one, and the
   engine's statement as `amplitudeDf`, NaN for none; `params` keeps its four numbers. Alternative: a code in the
   fourth number, which every reader of the hold would have to decode.
8. The reader reports the statement: `student(df)` where stated, absent at the default, though the default now
   differs by kind. Alternative: a label of the prior in force on every forest, which changes every reader of
   `amplitude`, print's among them. A public-surface choice worth the maintainer's nod.
9. The full printed block names the prior; the short block stays silent, as for a hold.
10. t_nu,.75 is computed in the engine, through the support library's Rmath wrapper, since the engine holds the
    map and every default. Alternative: R hands c over, against the architecture's rule that R derives none of
    the map.
11. `halfCauchyScale` is renamed `mixtureScale`, a name that stays true at any df.
12. The unstated path keeps today's shape and rate lines verbatim behind a branch on a stated df, so neutrality
    is structural, including under a compiler that contracts multiply-adds. Alternative: one expression in nu,
    bitwise at nu = 1.0 only by argument.
13. The exact test is new arms of forest-prior-exact.R: a flat-likelihood prior arm, numeric arms, and a
    drawn-factor arm. bcf-exact.R gains nothing, its functionals having no power over the prior's form.
14. The design notes are amended, not joined by a new one: the rule gains a row.
15. stan4bart is not changed. bartCause's `bcf()` states `student(1)` and its sd as [Sequencing](#sequencing)
    says, a commit of its own between this item and the flip (dec-B509).
16. On a forest with no basis, `student()` with no sd takes as median the prior median of the size the forest
    would otherwise have: today's default median until the flip, the held size's after. Alternatives: bcf's 2U or
    U always, which keeps bartCause's gaussian draws bit for bit with no sd stated but gives `student(1)` a size
    unrelated to the model's leaf prior, against dec-B509's "its sd the prior median"; or refusing it without an
    sd, a refusal at the binary default, where the drawn k leaves no held size.
17. On a forest with no basis the t keeps that forest's arrangement (leaf scale at the unit, the size in the
    coefficient's median) and `student(1)` takes today's lines, so it is the default operand for operand.
    Alternative: the multiplied forest's median match, which moves the leaf scale and every draw.
18. `student()` is refused on a lone forest, whose single-forest chain has no coefficient. Alternative: route a
    lone forest stating it through the K-forest chain at K = 1, a path nothing builds or tests today.
19. A bare `student` is `student()`, as a bare `fixed` is `fixed()`.

## Open calls

1. RULED, dec-B510 (VD 2026-10-29: "Use your revised recommendation."). `amplitude = student()` is a t with 1
   degree of freedom, the half-Cauchy; `student(df)` states a finite df of at least 1, anything else refused by
   name; the help says a bare `student()` means df = 1 here and an estimated df as a response family. The
   alternatives, df required at least 1 or any positive, were not taken.
2. SETTLED by dec-B509 (VD 2026-10-29: "OK, use your recommendation for this decision."). A forest with no basis
   is held by default in every model, and `student(df)` is accepted there as the opt-in that turns on a drawn
   coefficient, `student(1)` with `sd` 2 sd(y) being bcf's prognostic prior, which bartCause's `bcf()` states. The
   refusal this plan had recommended is gone; [On a forest with no basis (dec-B509)](#on-a-forest-with-no-basis-dec-b509)
   and Calls made 16 to 18 are the calls under it.
3. RULED, dec-B511 (VD 2026-10-29: "if (a) is not too much work we can support it, even if I don't see why."). On a
   live sampler a coefficient may change between held and drawn, between the normal and `student(df)`, and a t's
   df, at the next sweep; built as its own 1.0-0 item, coefficient-prior-live, after setmodel-parameters' slice 4.
   This plan stays at creation (R item 8).

## Not verified

- Nothing is built; the engine steps are from reading the code, not a prototype.
- The mixing of a `student()` coefficient is not measured; the tail's value to its user is dec-B262's prototype
  measurement, not repeated here.
- The drawn-factor arm's run time: the two-dimensional quadrature is new (the critic timed the script at 27 s in
  quick mode on two cores, about 0.8 s an arm, before it).
- The numeric arm's power against a normal drawn in place of the t (0.814 against 0.831 at four values) at 32
  seeds' standard error; the prior arm carries that mutant regardless.
- The pending plans were read at 8342bc31 and 8adc5d1a, and forests-shared-total on its branch, amended with this one;
  their calls may move.
- A held prognostic forest against the drawn default was checked under gaussian (dec-B509's check: not clearly
  worse for the effects at equal time); the mixing of `student(df)` at df above 1 on a forest with no basis is not
  measured.

## Landing note

Landed 2026-10-30 (dec-A245), rebased onto forest-after-merge's landing: `forest(amplitude = student(df))` as
both slices describe, with its help, NEWS clause, design notes and exact-posterior arms. Neutral for every call
main accepts, structurally and by the A/B on both builds; `student(1)` on a forest with no basis is main's
default bit for bit. x86-64 leg and CI green on the branch tip. The live change stays its own item,
coefficient-prior-live, through `$setForest(amplitude = )` (dec-B518).
