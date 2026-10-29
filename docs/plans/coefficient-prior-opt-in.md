# coefficient-prior-opt-in: a heavy-tailed prior on a forest's drawn coefficient

Status: PLANNED 2026-10-29 on 177fd104; amended the same day after its blind critique (build after amendments;
scratch/crit-cpo/report.md) and the coordinator's rulings on it (scratch/crit-cpo/my-rulings.md). Three calls held
for the maintainer ([Open calls](#open-calls)); slice 2 waits on the first two. Rulings: dec-B262 (the opt-in,
its user and its measurement) and dec-B271 (one argument, `amplitude`, beside `fixed()`; `student(3)` the form
shown), with dec-B18, dec-B253, dec-B254, dec-B261, dec-B266, dec-B290 to dec-B292, dec-B354, dec-B503 and
dec-B505.

agent: one opus implementer for both slices; one opus review over both that runs the mutants below, the same
reviewer checking the fix round unless the fixes change engine behavior. Tier "Changes draws" (a prior and the
engine's coefficient draw; [Process by risk](README.md#process-by-risk)). The blind critique that tier asks for
has run.
rng: NEUTRAL, bit for bit, for every call main accepts, by structure: a forest that states no df takes today's
lines, verbatim, behind a branch on whether a df is stated. A forest stating `amplitude = student(df)` is a new
prior, posterior-defining for those models: gated by new exact-posterior arms and a design-note amendment. No
re-record.
window: after the merge to main, before 1.0-0 (dec-B262: "Sure, opt-in for 1.0. Post main-merge."). Additive:
`amplitude = student(3)` is "could not find function" today.
budget: ~1000 lines in the tree. Slice 1 ~300 (engine ~75, support library ~3, bridge ~35, tests/cpp ~185).
Slice 2 ~640 (R ~130, help ~60, NEWS ~5, tinytest ~240, exact arms in forest-prior-exact.R ~160, design notes
~45). Records ~40. Out of the tree, ~300 lines of scratch A/B script.

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

Three questions for you:

1. What `student()` takes. I recommend that the degrees of freedom be required and at least 1: `student(1)` is
   the half-Cauchy you were shown measured, `student(3)` the form shown with the ruling, and a tail heavier than
   the Cauchy's was neither measured nor asked for (below about 0.01 the scale is too small for the computer to
   hold). The response family `student()` reads a missing df as "estimate it", which a coefficient's prior does
   not do; requiring it keeps that reading free. The alternatives are any positive df, or a default of 1.
2. Forests with no basis. Their coefficient is already half-Cauchy. Two views. Mine was to accept `student(df)`
   there too, one rule for every drawn coefficient. The critic's was to refuse it: your ruling's user and its
   measurement are about a multiplier; accepting it adds a second spelling of that forest's default and a rule
   for forests that share one total; and accepting it later breaks nothing. I now recommend refusing it, for the
   critic's last reason above all.
3. Changing it on a live sampler. My reading, which you would confirm, is that your `setModel` ruling makes a
   coefficient's prior, and whether it is held, parameters a live sampler should accept. Neither can be changed
   today, the hold included. I recommend building both together as their own 1.0-0 item and keeping this one to
   creation. The alternative is the prior's change here and the hold later, about 250 more lines.

The size is about 1000 lines with tests and gates; the ruling's 300 counted the code alone. About 13 agent hours,
and plans here run 1.5 to 2 times low.

## Goal

`forest(amplitude = student(df))` gives a multiplied forest's drawn coefficient a Student-t prior with `df`
degrees of freedom whose median size is the normal default's, at every door that takes several forests; the
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
- Consumers. bartCause's `bcf()` passes `amplitude = NULL` or `fixed()` on its two forests; stan4bart passes one
  `forest(n.trees = )` and no amplitude; neither changes.

## What 0.9-34 did

Nothing here. Probed in scratch/libs/cran: no `forest`, no `student` in the namespace, no `forests` argument and
no multiplied forest; its `cauchy` is `rbart_vi`'s random-effect prior. Nothing to mirror.

## Design

The plan is written for open calls 1 (a) and 2 (b): df required and at least 1, and refused on a forest with no
basis. What (a) of call 2 would add is listed with it.

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
- nu at least 1 (open call 1): below about 0.01 c underflows (t_nu,.75 is 6e28 at 0.01 and infinite at 1e-6).
  Infinite nu is the normal, refused as a second spelling of the default (Calls made 6).

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
   its node scale, its variance kappa, and takes `mixtureDf` nu and `mixtureScale` c from the median match; a held
   forest ignores a df. Where none is stated every output is today's, operand for operand. A df stated on a forest
   with no basis, and a stated df below 1 or infinite, throws `std::invalid_argument`
   (the backstop; R refuses each first). t_nu,.75 comes from the support library: `ext_quantileOfStudentT`, one
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

1. A constructor `student(df)`, exported only as `dbartsPriors$student` (an internal name of its own, the
   family's being `student`), class `dbartsStudentPrior`, one slot, valid for a single finite number at least 1.
   With no df: "student() on a coefficient takes its degrees of freedom, 1 or more: student(1) is a half-Cauchy,
   student(3) a lighter tail".
2. Vocabularies. `student` joins `fixed` in [`forestConstructors`](../../R/model.R), the `forests` entry of
   [`FOREST_ARGUMENT_VOCABULARIES`](../../R/model.R), [`processHit`](../../R/formulaTerms.R)'s vocabulary and the
   hint list of `forceCallerCode`. [`resolveFamily`](../../R/family.R) keeps the families first and drops the
   duplicate explicitly, `c(dbartsFamilies, dbartsPriors[setdiff(names(dbartsPriors), names(dbartsFamilies))])`,
   so `student` there is the family whatever the order; its comment says the two lists share that name.
3. [`validateForestAmplitude`](../../R/model.R) takes `fixed()` or a `dbartsStudentPrior`; a bare `student`, as
   a bare `fixed` is `fixed()`, is refused with item 1's message. Refused by name: a family object ("a forest's
   'amplitude' is fixed() or student(df); dbartsFamilies$student is the response's family") and anything else,
   naming both. `student(Inf)`, worded for a multiplied forest, whose default is the normal: "a t with infinite
   degrees of freedom is the normal, this forest's default: leave 'amplitude' out". A forest with no basis
   refuses any `student()` first, by kind: "a forest with no basis draws its coefficient under a half-Cauchy
   prior, its sd the median; 'amplitude' takes fixed() there". It runs at fit time, where the basis is known.
4. [`forestStatements`](../../R/model.R) adds `amplitude.df` only where some forest states one, so a default
   fit's record and every fit object made without it are as today. `params` and `heldAmplitude` do not change.
5. The one-forest refusal reads "which has no coefficient to give a prior", covering both values.
6. [`reportLeafPrior`](../../R/dbarts.R): `amplitude` is `student(df)` where stated, `fixed()` where held, absent
   otherwise; a fit's `leaf.prior` follows. `prior.sd.of` and `k.scale` do not change.
7. [`printedForest`](../../R/forestPrint.R): the full block's test becomes "is `fixed()`", and a stated forest's
   line ends "; coefficient student(3)". The short block stays silent, as for a hold.
8. `$setLeafPrior(forests = )` keeps refusing `amplitude` as fixed at creation (open call 3).

### Help, NEWS and docs

- man/forest.Rd `amplitude`: the default by kind (normal with a basis, half-Cauchy with none), `student(df)` on a
  multiplied forest and the median rule, what a tail buys and when; man/dbartsPriors.Rd: a `student(df)` item,
  naming the family of the same name in `dbartsFamilies`; man/dbartsSampler-class.Rd: `getLeafPrior`'s
  `amplitude` entry.
- NEWS, in 1.0-0's forest item: one clause.
- Amended, not new notes: [1. The rule](../design/forest-multiplier.md#1-the-rule) gains the drawn
  `student(df)` row and section 7's heavy-tail bullet goes;
  [The calibration map, general in K](../design/multiplier-combiner.md#the-calibration-map-general-in-k) gains
  the mixture's nu and the median match.

## Interaction with the three pending plans

forest-after-merge (8adc5d1a; its tip 8c189a18 does not move the plan):

- Slice 1, the long form, adds `normal` to the same three vocabularies this item adds `student` to: a textual
  conflict, both names kept. With both in, `forest(leaf.prior = student(3))` resolves and must be refused,
  pointing to `amplitude`: "a coefficient's prior is 'amplitude', as amplitude = student(3); 'leaf.prior' is the
  forest's size, normal(sd = )". Whichever lands second adds the refusal and its pin.
- Slice 2's `restateForestBasisScale` re-resolves [`forestMapEntry`](../../src/bartcore/combiner.hpp) with "the
  forest's own default statement": it must keep the statement's `amplitudeDf` and default only the sd. Under
  Calls made 1 the coefficient's prior does not depend on the column's scale, so its "no amplitude prior moves"
  then holds for a `student()` forest. Whichever lands second pins: a `student(3)` numeric forest swapped with
  `updateBasisScale = TRUE` keeps the reader's `amplitude`, and its draws equal the creation twin's (a
  restatement that drops the df fails it).
- No shared facade edit: this item adds no virtual.

forests-shared-total (e48fdcf3, dec-B505: an equal part of the total per forest, each sharer its own leaf
scale):

- Under open call 2 (b), the plan's case, no sharer states `student()`: a forest with no basis refuses it by kind,
  and the sharers' coefficient keeps its half-Cauchy. Its D10 rules on `fixed()` are untouched.
- Under (a), the sharers' one coefficient takes one statement for the group: every sharer states the same
  `student(df)` or none, a mix refused by name as D10 refuses a mixed hold. Checked against dec-B505: the shared
  draw's design row is the sum of the sharers' fits, each at its own leaf scale, and that is the data block; the
  coefficient's prior reads no leaf scale. Its c is the group's median over t_nu,.75 and v | a reads only a, so a
  per-sharer scale changes neither, and its engine item 2 refreshes one v with the first sharer's nu unchanged.
- Shared files: [`forestMapEntry`](../../src/bartcore/combiner.hpp) (its share on a sharer's node scale, this
  item's nu and c on a multiplied forest's coefficient; disjoint branches), `drawAmplitudes` (its group draw, this
  item's stated branch) and benchmarks/R/forest-prior-exact.R (each adds arms). Textual conflicts; the second to
  land rebases and, under (b), pins `student(3)` refused on a sharer by the kind message. This item no longer
  edits bcf-exact.R.
- Engine features merge serially; either may land first.

setmodel-parameters-not-structure (8342bc31):

- Slices 1 to 3 refuse or never reach a sampler with amplitudes; the coefficient's prior rides the forests'
  record on the control, which no model comparison reads. Its family checks meet this item's vocabulary change:
  whichever lands second pins `$setModel` with a model of `family = student(5)` resolved as the family.
- Slice 4 opens `$setModel` on such samplers for the tree and residual priors and never routes through the
  leaf-scale write. Whichever lands second pins: on a `student(3)` forest, `$setModel` with a new tree prior
  leaves the reader's `amplitude`, the coefficient and v as they were, and the next draws equal a copy's. Its
  `$setTreePrior(forests = )`, evaluated in the widened vocabulary, refuses `forest(amplitude = student(3))` as
  any knob but `base` and `power`.
- Open call 3 (b), or the later item it recommends, would add a second leaf-side door beside its tree-prior one;
  the coefficient's prior is not the model's, so `$setModel` never carries it.

## Slices

### Slice 1: engine and bridge

Steps: engine items 1 to 6 and the bridge; tests/cpp below. Nothing in R reaches the new branch.

RNG: neutral, every draw of every call. Gates: tests/cpp plain and under address and undefined sanitizers; the
equivalence trio bitwise; the full tinytest suite.

Cost: ~300 lines; 2.5 agent hours to build.

### Slice 2: R surface, help and the exact arms

Steps: R items 1 to 8, help, NEWS, the design-note amendments, tinytest, the exact arms. Waits on open calls 1 and
2 (the plan is written for 1 (a) and 2 (b)).

RNG: neutral for every call main accepts; posterior-defining for a forest stating `student()`. Gates: the whole
list below on the slice's tip.

Cost: ~640 lines in the tree and ~300 of scratch A/B; 4 agent hours to build.

Both: one review over both slices, 2 hours; the gate battery 3 hours, most of it compute; the fix round 1.5;
records 0.5. About 13 agent hours; plans here run 1.5 to 2 times low, so budget 20 to 26. The x86-64 leg before
landing is the coordinator's (GCC contracts multiply-adds by default).

## Tests and gates

tests/cpp (slice 1):

- The map: no df stated gives today's entry bitwise for every kind, held and drawn; `student(1)` and `student(3)`
  on a numeric and a factor forest leave the node scale and kappa bitwise and give c against literals
  (0.6744897501960817, 0.8818100603, and sqrt(1/2) times each); a held forest ignores a df; NaN states none; a df
  of 0, a negative, 0.5 and Inf throw, and so does any df on a forest with no basis.
- The draw: a K = 3 fixture stating no df draws bitwise what it drew before the slice (a recorded vector); with a
  student forest, v after each sweep against IG((nu + q) / 2, (nu c^2 + |a|^2) / 2) read through `TestPeer` on a
  seeded generator; v starts at kappa; at nu = 1, 10000 sweeps leave every a and v finite.
- State: serialize and restore carry v on a student forest; a normal donor's variance installs as v; a student
  donor's v does not reach a normal recipient.
- The calibration reports the stated df and the mixture scale; NaN where none is stated.

tinytest, a new test-forest-amplitude-student.R, small n:

- Every door: a formula term, `forests = list()`, a data object's `bases` with `forests`, and `dbartsSpec()`;
  gaussian, probit and logistic; a numeric and a factor basis.
- The reader's `amplitude` `identical()` to `student(3)`; `fit$leaf.prior`; the full block's line and the short
  block's silence; `fixed()` still printed as held; a default fit's forests record without `amplitude.df`.
- Refusals by message, counted with `withCallingHandlers`: `student()`, a bare `student`, `student(0)`,
  `student(0.5)`, `student(NA)`, `student(c(1, 3))`, `student(Inf)`, `dbartsFamilies$student(3)`, `student(3)` on
  a forest with no basis (open call 2), `amplitude` on one forest,
  `$setLeafPrior(forests = list(forest(amplitude = student(3))))`.
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
  one numeric and one factor, drawn; P(|a| <= q_p) at p = 0.25, 0.5 and 0.9 against the closed form for
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
- The kind refusal removed (the no-basis refusal test; the engine backstop then throws under the generic text).

Gates, on the branch's own libraries:

- `R CMD INSTALL --preclean`, plain and `--configure-args=--enable-reference-build`, `buildInfo()$mode` checked.
- tests/cpp plain and under address and undefined sanitizers; the R-loaded sanitizer run over the new test file
  (`Rf_qt` and the new update are reachable from R).
- A neutrality A/B (scratch, ~300 lines): a 177fd104 library against the branch's, on the reference build and on
  the shipped build, two and three chains, every recorded channel compared with `identical()`: bartCause's shape
  drawn and held under gaussian, probit and logistic; three forests (none, a factor, a numeric column) with and
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

- One more prior to explain: the default differs by kind (normal with a basis, half-Cauchy without), which the
  help now says outright.
- A heavy-tailed coefficient mixes along the ridge it shares with the forest's leaves more slowly than a normal
  one; the opt-in's user accepts that for a scale that drifts. Not measured here.
- `student` names two things, the response family at `family =` and the coefficient's prior at `amplitude =`,
  each resolved by where it is written.

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
4. v starts at kappa, the variance the map gives the block today. Alternative: the prior median of v, which
   starts a stated forest off the default's first draw.
5. The constructor lives in `dbartsPriors`, the home of `fixed()`, under an internal name of its own; the family
   door keeps the families first and drops the duplicate name explicitly. Alternative: taking the family object
   `student(3)` as a coefficient prior, a response distribution whose missing df means "estimate it".
6. `student(Inf)` is refused, worded for a multiplied forest; a forest with no basis refuses every `student()` by
   kind before it. Alternative: accept Inf as the normal, a second spelling of the default.
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
15. bartCause and stan4bart are not changed; a `student()` opt-in on bartCause's `bcf()` is its own addition.

## Open calls

1. What `student()` takes on a coefficient. The family `student(df = NULL)` reads a missing df as "estimate it";
   the prior constructors `chisq()` and `chi()` have defaults, `invchi()` requires its `scale`. A df near 0 cannot
   be held: c underflows below about 0.01.
   - (a) Recommended. `df` required, a finite number at least 1: the half-Cauchy dec-B262 measured is the
     heaviest tail, `student(3)` as dec-B271 showed it. Cost: one more thing to type; no tail heavier than the
     Cauchy; a later drawn df can take `student()` with nothing given, as the family does.
   - (b) `df` required, any positive finite number, refused only where c is not a positive normal double, with a
     finite-draw test at the smallest df taken. Cost: tails nobody measured, and a guard at an arbitrary-looking
     machine limit.
   - (c) `df = 1` by default, the measured half-Cauchy, at least 1. Cost: `student()` means half-Cauchy here and
     "estimate" in the family; a drawn df later needs another spelling.
2. `student(df)` on a forest with no basis, whose coefficient is half-Cauchy today (dec-B18).
   - (a) The planner's first view: accepted, one rule for every drawn coefficient, a t at the median the
     forest's sd sets, `student(1)` exactly today's prior. Cost: `student(1)` restates that forest's default, a
     second spelling of it; a lighter tail nobody asked for; forests-shared-total gains a rule (every sharer
     states the same df or none); a no-basis prior arm and the engine's kind-none branch, ~60 lines; `student(Inf)`
     then needs a message for a forest whose default is not the normal.
   - (b) Recommended, the critic's view, which I now share: refused there by kind, "a forest with no basis draws
     its coefficient under a half-Cauchy prior". dec-B262's user redraws a multiplier and its measurement is of
     multiplied forests; accepting it later is an addition that changes no call, where taking it now cannot be
     undone. Cost: no lighter tail for a causal model's prognostic forest until asked for.
3. Changing the coefficient's prior on a live sampler. That the coefficient's prior, and its hold, are parameters
   under dec-B254 ("parameters yes, structure no", the hand between drawn and held a parameter change) is this
   plan's reading, for the maintainer to confirm; if so, dec-B254 asks for it "for 1.0", and
   `$setLeafPrior(forests = )` refuses both today.
   - (a) Recommended. This item stays at creation; a new TODO item for 1.0-0 builds the hold and the prior on a
     live sampler together, through `$setLeafPrior(forests = )`, after setmodel-parameters' slice 4. Cost: one more
     1.0-0 item, ~400 lines.
   - (b) The prior's change here (normal, `student(df)`, a new df), the hold later. Cost: ~250 more lines, an
     engine setter and a facade virtual, and the two halves of one ruling built apart.
   - (c) Neither is a parameter: both stay fixed at creation, as built. Cost: none in code; dec-B254 read narrowly.

## Not verified

- Nothing is built; the engine steps are from reading the code, not a prototype.
- The mixing of a `student()` coefficient is not measured; the tail's value to its user is dec-B262's prototype
  measurement, not repeated here.
- The drawn-factor arm's run time: the two-dimensional quadrature is new (the critic timed the script at 27 s in
  quick mode on two cores, about 0.8 s an arm, before it).
- The numeric arm's power against a normal drawn in place of the t (0.814 against 0.831 at four values) at 32
  seeds' standard error; the prior arm carries that mutant regardless.
- The pending plans were read at e48fdcf3, 8342bc31 and 8adc5d1a, before their builds; their calls may move.
