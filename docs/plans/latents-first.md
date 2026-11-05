# latents-first: a sweep draws its augmentation variables before the trees

Status: PLANNED 2026-11-03 on main 730dd0adf (dec-B586); blind critique 2026-11-05 on main 7937d252f, BUILD AFTER
AMENDMENTS; amended 2026-11-05 ([Critique and amendments](#critique-and-amendments)). Not built.

agent: one opus implementer each for slices 2a and 2b (engine) and one sonnet implementer for slices 1 and 3 (a
gate, then the documentation). Review by tier ([Process by risk](README.md#process-by-risk)): slice 1 one sonnet
review (gate-only); slices 2a and 2b one opus review each that runs that slice's mutants below, the same
reviewer checking the fix round unless the fixes change engine behavior; slice 3 one sonnet review against the
engine as built. Slice 4 is the orchestrator's.
rng: POSTERIOR-CHANGING overall. For a plain fit of a family with augmentation variables the draws move and the
posterior does not (shifting); for gaussian, heteroscedastic gaussian, multinomial and aft without censored rows
the draws stay bit for bit. The posterior changes where a caller changes the sampler between sweeps and its
outer step integrates the latents out (now exact, off before), and on a logistic sampler handed a new offset
(a defect fixed). Gates: [Gates](#gates).
window: before 1.0-0. Build order: coefficient-prior-live, then latents-first, then frequency-weights, then
host-neutral-engine S2a onward ([Order and overlaps](#order-and-overlaps)).
budget: ~3100 lines in the tree: slice 1 ~400 (a new exact gate ~340, benchmarks README ~60); slice 2a ~650
(engine ~80 added and ~40 removed, tests/cpp ~330, tinytest ~180, workflow ~10); slice 2b ~1250 (engine ~70
added and ~380 removed, bridge ~140 removed, tests/cpp ~150 added and ~280 removed, tinytest ~200 added and
~430 rewritten or removed, a benchmark header ~20, one snapshot file and three equivalence baselines
re-recorded); slice 3 ~750 (help in the docstrings and Rd ~340, vignette ~50, header ~30, embedding page ~70,
architecture and design notes ~180, a new design note ~70, NEWS ~15); records ~50. About 56 agent-hours in all
([Slices](#slices)).

## Summary for the maintainer

Every family that augments its likelihood - the probit and ordinal latents, the Polya-Gamma precisions of
logistic, binomial and negative binomial, the Student-t scale mixing, the aft imputed times - draws those values
once a sweep. Today it draws them after the trees. So when the code that runs the BART sampler inside a larger
sampler changes something between sweeps (an offset, a predictor, sigma), the next sweep's trees fit values
drawn before the change. That is right only when the larger sampler's own step reads those values; when its
step averages them out, as treatSens's treatment model, bairrtt's ability step and the vignette's fourth recipe
do, the answer is off. Drawing them first, right after the point where the larger sampler makes its changes, is
right either way and lets every setter stop drawing. It costs nothing extra a sweep for a one-forest fit; a
two-forest probit or logistic fit blends its forests' fits once more a sweep, a cost the build measures against
today's.

What changes for a user:

- Every family with augmentation variables draws differently at the same seed, once. Gaussian, multinomial,
  heteroscedastic gaussian and aft without censored rows draw exactly as before (checked).
- No setter draws any more: setResponse, setActiveRows, setWeights and setState store what they are given and
  the next sweep draws the rest. getLatents after a run returns the values the last sweep's trees were drawn
  from; between a setter and the next run, the values from before the setter, except that swapping in a new data
  object starts them as a new sampler would, and a new aft response leaves a censored row at its new bound.
- A Student-t setResponse keeps the drawn degrees of freedom and scales, where today it resets them.
- Swapping in a data object keeps the drawn ordinal cutpoints, Student-t degrees of freedom and negative binomial
  shape (dec-B589), always: the ruling's restart for a changed number of categories cannot arise, since a
  sampler already refuses data whose categories differ from those it was created with.
- A saved state no longer records which rows were at weight zero, nor fingerprints of the weights and censoring
  it was stored under. Their only use was redrawing at a restore, which nothing needs once each sweep draws
  first. The aim of the 2026-10-06 ruling that added the record (dec-B277), a restore under other weights that
  is exact, stays and gets stronger: setWeights and setState now undo each other exactly in either order.
- On a small probit model whose larger sampler averages the latent out, today's order is off by about 15
  standard errors; drawing first is exact (checked on the real engine). A logistic model averaged the same way is
  very likely off too; the first slice measures each case on today's package and on the new order.
- Two defects found while planning are fixed with it: a logistic sampler handed a new offset kept fitting the
  old one for a sweep (it reads the old offset through a buffer the new one has already overwritten), and a mask
  installed from the per-sweep hook was ignored by that sweep.
- Documentation states the rule in each place a user or package author meets it: getLatents and every setter,
  the help page and vignette on running the sampler inside a larger one, the C header.

0.9-34 drew probit's latent after the trees, redrew it when given a new response and kept it across setOffset,
setPredictor and setWeights, so for probit this is the same behavior made correct, not a regression, except that
setResponse no longer draws; NEWS says so. The other families are new in 1.0-0.

The work comes in five pieces: the new exactness check, the reordering itself, the setters' draws and the state
records removed, the documentation, and the consumer runs and records; about 56 agent-hours.

## Goal

After this lands, each sweep of every sampler draws its augmentation variables (with the ordinal cutpoints and
the negative-binomial shape drawn beside them) after the per-sweep hook and before anything reads them, so a
change made between sweeps or inside the hook never leaves the trees a stale latent. No setter draws. The
sampler is exact whether an outer step conditions on the latents or integrates them out, and the contract is
written where users and package authors read it.

## Context

- Ruling: dec-B586 (direction, cost, documentation places, Student-t nu) and dec-B589 (setData keeps the family
  parameters); related dec-B583 (accept-or-keep Polya-Gamma latent), dec-B580, dec-B581 (draw callback);
  superseded in mechanism here, their goal kept: dec-A160 and dec-B277 (the Student-t redraws under weights and
  the state's zero-weight record). Evidence: scratch/stc-chk/report.md and its probes; the critique
  scratch/ltf-crit/report.md and its probes.
- The sweep: [`Chain::runSweeps`](../../src/bartcore/chain.hpp), with [`drawLevelShift`](../../src/bartcore/chain.hpp),
  [`drawForestRescale`](../../src/bartcore/chain.hpp), [`combinedFits`](../../src/bartcore/chain.hpp),
  [`armRefreshCancel`](../../src/bartcore/chain.hpp), [`sweepVarianceForest`](../../src/bartcore/chain.hpp) and
  [`storeSample`](../../src/bartcore/chain.hpp); the initializer [`growForestFromRoot`](../../src/bartcore/chain.hpp).
- The families: [`ResponseModel`](../../src/bartcore/model.hpp) and its `refreshLatents`;
  [`ProbitResponse`](../../src/bartcore/model.hpp), [`OrdinalResponse`](../../src/bartcore/model.hpp),
  [`LogisticResponse`](../../src/bartcore/model.hpp), [`NBResponse`](../../src/bartcore/model.hpp),
  [`TResponse`](../../src/bartcore/model.hpp), [`AFTResponse`](../../src/bartcore/model.hpp),
  [`MultinomialResponse`](../../src/bartcore/model.hpp), [`GaussianResponse`](../../src/bartcore/model.hpp).
- Copy-on-set: the bridge copies a new offset, response or weights into the buffer it already owns
  ([`adoptVector`](../../src/R_interface_bartcore.cpp)), so a family's old and new pointers are the same memory.
- Embedded use: [bart-as-a-component.md](../design/bart-as-a-component.md) and the help page
  [dbarts-embedding.Rd](../../man/dbarts-embedding.Rd); the vignette's recipes 4 and 6
  ([dbarts-as-a-component.Rmd](../../vignettes/dbarts-as-a-component.Rmd)).
- Exact-gate template: [`mask-redraw-exact.R`](../../benchmarks/R/mask-redraw-exact.R) (pointer-level loop on the
  real engine against an enumerated one-tree posterior).
- Branches: coefficient-prior-live (built, landing records on its branch, not merged), frequency-weights (plan
  516d02e40), host-neutral-engine (plan dd24d7670; S1 landed), approximate-draws (plan); consumers stan4bart,
  treatSens, bartCause and bairrtt.

## Design

### The question

Where in a sweep are the augmentation variables drawn, so that a change made by the code driving the sampler
between sweeps reaches the next sweep's trees?

### Sweep order

Today, per sweep: the per-sweep hook; mean weights under a variance forest; the level shift; the probit
rescaling step; the trees, forest by forest (a multi-forest combiner draws its per-forest glue, and the
multinomial its latents, inside this loop); the total fits; **the augmentation draw** (the ordinal cutpoints
then z, the negative-binomial shape then omega, the Student-t lambda then nu, all inside `refreshLatents`); the
working response and weights pointers read again and the vector-leaf statistics invalidated where the weights
vary per sweep; sigma; the combiner's glue; the variance forest; k and DART; the record.

After, per sweep:

1. The per-sweep hook (the caller's conditioning point; any setter may run here).
2. **The augmentation draw**, against the combined fits as the previous sweep left them (trees, glue, k), with
   its cancel point armed exactly as today; inside it the same pairs in the same order (cutpoints then z, shape
   then omega, lambda then nu). The combined fits are formed here only for a family that draws (the existing
   test `latents() != nullptr`): a single forest's are its totals, a bare pointer; a coupling's blend
   ([`ForestCombiner::combinedFits`](../../src/bartcore/combiner.hpp)) is an O(n) pass per forest into the
   combiner's scratch, and the multinomial's softmax is not formed here at all (its family draws nothing at
   step 2).
3. The working response and working weights pointers read from the family, every sweep, then the vector-leaf
   statistics invalidated where the family's weights vary per sweep.
4. Mean weights under a variance forest (they replace the weights pointer for the trees).
5. The level shift (fit-preserving; reads no latent).
6. The probit rescaling step (reads z, the offset and the fits; scales z, the leaves and k).
7. The trees and total fits; the per-forest glue and multinomial latents inside the loop, unchanged.
8. Sigma; the combiner's glue and post-combine move; the variance forest; k and DART. These read the combined
   fits formed again after the trees, as today; never step 2's pointer, which on a coupling holds the
   pre-tree blend.
9. The record and the draw callback.

Reading the pointers every sweep (step 3) replaces today's read once per run and after the end-of-sweep draw;
that is what makes a setter in the hook take effect in the same sweep (probe below).

| family | augmentation drawn at step 2 | before | after |
|---|---|---|---|
| probit, hazard.probit, hurdle's zero part | z per active row | after trees | first; rescale reads it |
| ordinal | cutpoints (marginal MH, z integrated), then z | after trees | first, pair kept |
| logistic, binomial, hazard.logistic | omega ~ PG(count, eta) per active row | after trees | first |
| nbinom | shape r (grid, omega integrated), then omega | after trees | first, pair kept |
| Student-t residuals | lambda every row, then nu | after trees | first |
| aft | censored log-times per active row | after trees | first; reads the variance surface the previous sweep drew, as today |
| BCF probit, BCF logistic | z or omega against a mu + b_z tau | after trees | first, against the previous sweep's glue |
| multinomial | none at step 2 (per-forest omega in step 7) | unchanged | unchanged, bit for bit |
| gaussian, heteroscedastic gaussian, hurdle's positive part | none | - | bit for bit |
| variance forest (gaussian or aft) | none of its own | - | its sweep stays in step 8 |

The level shift and the rescaling step could sit on either side of the mean weights; they stay in today's
relative order. The draw goes before the rescaling step, not after it: placed after, the rescaling reads a z
drawn before the caller's change (probe 3 shows that placement off by as much as drawing last).

### growForestFromRoot

The grow-from-root initializer runs its own sweep (trees, then latents, sigma, glue, k, DART). Decision: it
takes the same order, drawing the latents at the start of each grow sweep. Today its last latent draw is never
read by anything but getLatents, since the run that follows draws first; drawn first, the grown trees fit latents
drawn at the fits they start from, and there is one order to describe.

### Creation-time latents and the first sweep

Creation keeps its deterministic starts (probit and ordinal z at the sign of y, logistic omega at w/4, the
negative-binomial omega at its PG(y + r, 0) mean, Student-t lambda at 1 and nu at its grid's middle, aft censored
times at their bounds); no draw is added at creation. They are now read only by getLatents and storeState before
the first run. The first sweep draws the latents from the starting fits (zero, or trees from the prior or a warm
start) before its trees, where today its trees fit the deterministic start.

### Store, restore, copy, saveRDS

A stored latent is the value the stored trees were drawn from, a valid joint draw with them. Under the new order
no next sweep reads the previous latent (each family's conditional draws from the fits, sigma and its family
parameters, and the two collapsed draws integrate the latent out), so a restored or copied sampler continues the
original bit for bit as before, and the stored latents serve getLatents after a restore. The state format does
not change, except that the Student-t record of zero-weight rows and the weights and survival digests are
retired (below); a development-build state carrying them loads and the attributes are ignored (only shipped
versions are owed compatibility; no migration).

### getLatents

After a run: the values the last sweep's trees were drawn from. Between a setter and the next run: the values
from before the setter, with two exceptions: setData, a new data object whose rows may differ in number, starts
the latents as creation does (as today); and aft setResponse, where a censored row holds its new bound (data)
until the next sweep draws it. setState installs the stored values (an aft row the state has as an event and this
sampler censors reads its bound, as today). In the per-sweep hook: the previous sweep's values. In the draw
callback: this sweep's, the ones its trees were drawn from (today the callback saw the values the next sweep's
trees would fit).

A new fixed negative-binomial shape through setModel
([`NBResponse::installShape`](../../src/bartcore/model.hpp)) today restarts omega at its cold start (y + r)/4.
Decision: it keeps omega and rebuilds the working response and composite under the new shape. Nothing reads omega
before the next sweep draws it from the new shape, so the draws are the same either way; keeping it removes a
third exception from the rule above and a cold start no reader needs.

### Setters

Each setter stores what it is given, rebuilds what is deterministic, and draws nothing:

| setter | today | after |
|---|---|---|
| setResponse, probit / ordinal / logistic / nbinom | redraws every active row | stores y; ordinal recomputes its scales and keeps the cutpoints; nbinom keeps r and recomputes its kernel; the working response rebuilt from the latents in force |
| setResponse, aft | redraws censored rows | installs the times; censored rows at their new bounds |
| setResponse, Student-t | lambda reset to 1, an estimated nu to its grid's middle | lambda and nu kept |
| setOffset, logistic | shifts the working response by old minus new offset, reading the old one through the overwritten buffer | rebuilds the working response from omega and the new offset |
| setOffset, others; setPredictor; setSigma; setForestBasis | no draw | unchanged |
| setActiveRows | redraws the rows it switches back on (Student-t: lambda) | stores the mask and keeps today's deterministic rebuilds: ordinal recomputes its scales, nbinom its count kernel, logistic and nbinom their composite, Student-t its membership marks and composite, aft hands the mask to its gaussian part |
| setWeights, logistic | inactive rows cold-started, then every active row redrawn | stores the counts; the working response rebuilt from omega |
| setWeights, Student-t | redraws rows whose weight turns positive | stores the weights; membership marks and the composite recomposed |
| setSurvivalStatus | no draw | unchanged |
| setData | cold-starts the latents and restarts the cutpoints, nu and the nbinom shape | latents start as at creation (rows may differ in number); cutpoints, nu and shape kept, always (dec-B589, [Ruled](#ruled)) |
| setModel, a new fixed nbinom shape | omega restarted at (y + r)/4 | omega kept; working response and composite rebuilt under the new shape |
| setState reconciliation | logistic omega and Student-t lambda redrawn on a weights-digest mismatch (Student-t by the stored zero-weight record); aft censored rows redrawn on a survival-digest mismatch | none; the next sweep draws them |
| sampleTreesFromPrior | (comment: follow with setResponse for a true restart) | the next sweep draws the latents from the new trees |

Retired with the redraws: `redrawEntering`, `reapplyWeights`, `zeroWeightRows`, `reapplySurvivalStatus`,
`weightsDigest` and `survivalDigest` (the model members, the chain, facade and sampler methods and virtuals,
their tests and the tests/cpp facade spies), the bridge's `weights.zero`, `weights.digest` and `survival.digest`
state attributes, the digest encoders, and the setState reconciliation that reads them. The digests' only readers
are that reconciliation's two mismatch flags. The Student-t record existed because a setter that drew needed to
know which rows were out (dec-B277); with every latent drawn at the next sweep there is nothing to reconcile. A
development-build state carrying any of the three attributes loads and the attribute is ignored.
`findReactivatedRows` stays only if another caller needs it.

### Draw callback and the C header

The draw callback (dec-B581) fires at the record, unchanged; `dbarts_draw` carries no latents. The header's latent
reader [`dbarts_sampler_getLatents`](../../inst/include/dbarts/dbarts.h) gains the sentence of the getLatents
section above. The run's interrupt paragraph changes: a stop inside the latent draws, which now open a sweep,
leaves no tree of that sweep drawn and each latent this sweep's draw or the previous one (a sampled
negative-binomial shape's refresh put back whole), a valid state to run on from. The flat C API has no sweep
hook, so a C caller makes its changes between runs; nothing else in the header moves, and no signature changes.

### Cost

The draw moves; no draw is added. On a single forest nothing else is added either, the combined fits being the
forest's totals. A two-forest probit or logistic fit forms its blend twice a sweep where it formed it once (one
O(n) pass per forest), and a multinomial or gaussian-family fit forms nothing extra (step 2 forms the fits only
for a family that draws). Saved: the setter draws, which a caller swapping the response every sweep pays today
(about 22 percent of a sweep for probit and 50 percent for logistic at n = 5,000 and 50 trees,
scratch/stc-chk/report.md), and logistic setWeights' redraw. Slice 2a measures it: a same-machine A/B of main
against the branch, alternating arms, per-round minima, loadavg stamped, at most 2 cores, on probit, logistic,
BCF probit and multinomial at n = 5,000 and 50 trees, plain fits. Pass: each ratio within the 3 percent
single-run grain, re-run before a lone flag counts; a BCF ratio above it is reported with the blend's measured
share, not treated as a failure.

### RNG

Draw-changing for every family with augmentation variables, with or without setters. Bit for bit for gaussian,
heteroscedastic gaussian, multinomial, aft with no censored row (probe 6); BCF gaussian and hurdle's positive
part are predicted to be, having no latent. Posterior-changing for embedded use whose outer step integrates the
latents out and for logistic setOffset; shifting for a plain fit.

### Accept-or-keep latents (approximate-draws)

Not built here. The constraint the approximate-draws plan must meet: an accept-or-keep Polya-Gamma latent is
exact for an outer step that integrates it out only where it is a fresh draw from its current conditional, so
every row whose conditional a setter or a restore moved since the last sweep (every row for setOffset, setSigma,
setForestBasis, setResponse, setData, setModel's shape and setState; the moved rows for setPredictor; the
changed rows for setWeights and setActiveRows; rows changed inside the per-sweep hook included) is drawn
exactly at step 2 of the next sweep, and the mark is cleared only after that draw. With this plan the setters stay draw-free; they set marks.

### Consumers

| consumer | outer step | effect | to re-record |
|---|---|---|---|
| stan4bart | conditions on z (contract A): setOffset, run(0, 1), getLatents | draws move, exact before and after | nothing pinned: its binary test checks a tolerance, its seed tests compare fits with themselves |
| rbart_vi (dbarts) | conditions on z | draws move; init reads the creation-time z, unchanged | nothing pinned moved under the prototype |
| treatSens, trt.model = bart() | integrates z out: setOffset then run | now exact; draws move | its tests' pinned values, if any move (not run, below) |
| treatSens outcome model | gaussian setResponse | bit for bit | nothing |
| bartCause | BCF and probit fits through run; grouped through stan4bart | draws move | the probit snapshots in its regression tests |
| bairrtt | integrates z out: its ability step scores the probit assignment model with pnorm, installs through updatePredictorPerObservationJointly, then run(0, 1) | off on main, now exact; draws move | nothing pinned (read): its tests check shapes, dimensions and tolerances |

### 0.9-34 and NEWS

0.9-34 drew probit's latents after the trees, redrew them at setResponse and kept them at setOffset, setPredictor
and setWeights; it had no other augmenting family. NEWS, USER-VISIBLE CHANGES of 1.0-0:

> A sampler draws its augmentation variables (the probit and ordinal latents, the Polya-Gamma precisions, the
> Student-t scales and aft's imputed times) at the start of each sweep, before its trees, and no setter draws
> them. A larger sampler that changes the offset, predictors or sigma between sweeps is exact whether its own
> step reads `$getLatents()` or integrates the latents out; in 0.9-34 the second was not, the next sweep's
> trees fitting latents drawn before the change. `$setResponse` on a probit sampler no longer redraws the
> latents, as it did in 0.9-34; `$getLatents()` returns the earlier values until the next sweep draws them.

## Documentation

Shipped text states the contract; history stays in this plan and the ledger.

- [`dbartsSampler$getLatents`](../../man/dbartsSampler-class.Rd) (R/dbarts.R docstring and the Rd value section):
  "The values the last sweep's trees were drawn from. Each sweep draws them first, after the per-sweep hook and
  before its trees, so between a setter and the next run they are the values from before the setter (setData
  starts them as a new sampler does; an aft setResponse leaves a censored row at its new bound) and the next
  sweep draws them afresh."
- setResponse, setOffset, setPredictor, setSigma, setForestBasis, setWeights, setActiveRows, setData, setModel
  (a new negative-binomial shape), setState (docstrings and Rd): one shared sentence, "This draws nothing; the
  next sweep draws the augmentation variables against the change." Removed: the Student-t redraw under weights,
  the logistic count-swap redraw and inactive cold start, the reactivation redraw and the joint-draw advice under
  the mask, the aft status-change redraw, the setState reconciliation against the destination's weights and
  censoring status (the paragraph that says a restore is reconciled against the destination's own) and the
  weights.zero record; the undo paragraph
  becomes "a weight change is undone by setWeights and setState in either order, bit for bit". The Student-t
  setResponse text says nu and the scales are kept. setData's text says it keeps the trees, the ordinal
  cutpoints, the Student-t degrees of freedom and the negative-binomial shape, and starts the latents as a new
  sampler does.
- The predict-then-setPredictor advice: adds that the accepted predictors reach the next sweep's latents, so a
  Metropolis step that scores proposals with the latents integrated out is exact.
- [dbarts-embedding.Rd](../../man/dbarts-embedding.Rd): a new subsection "When the augmentation variables are
  drawn": the order, the two kinds of outer step (one that reads getLatents and one that integrates them out) and
  that both are exact; "Which fit to read back" says what getLatents returns between sweeps.
- Vignette recipe 4: one paragraph, that scoring with the latent integrated out (pnorm) is exact because the next
  sweep draws sB's latents under the accepted theta. Recipe 6: one sentence, that z and f read after a run are a
  joint draw, z being the value f was drawn from; and its closing sentence's hard-coded numbers ("about 0.29,
  roughly seventeen oracle standard errors"), recorded from a probit run whose draws move, re-measured on the
  built engine and replaced.
- [`mask-redraw-exact.R`](../../benchmarks/R/mask-redraw-exact.R)'s header (slice 2b, with the change it
  describes): the condition that a row switched back in has its latent redrawn before the trees next move is now
  met by the next sweep drawing every latent first, not by the setter; the measured size of a sampler that keeps
  a stale latent stays as the gate's sizing.
- [`dbarts_sampler_getLatents`](../../inst/include/dbarts/dbarts.h) and the run's interrupt paragraph: as in
  [Draw callback and the C header](#draw-callback-and-the-c-header).
- [One sweep](../architecture.md#one-sweep): the step list rewritten in the new order; the RNG paragraph drops
  "the probit latent redraw on setResponse"; [Callbacks and errors](../architecture.md#callbacks-and-errors) says
  what the hook and the callback see.
- Design notes that state the old order or a setter draw: ordinal.md (draw order, setResponse),
  negative-binomial.md (r first, setResponse), robust-errors.md (setResponse cold init, zero-weight rows),
  active-rows-mask.md (reactivation redraws), aft-status-setter.md, aft-variance-forest.md, weighted-logistic.md,
  probit-k-scale-move.md (the rescale reads the fresh z), bart-as-a-component.md (the contract), survival.md,
  empty-leaf-veto.md where it names the order, and public-surface.md, weighted-logistic.md and aft-status-setter.md
  where they name the weights and survival digests or the weights.zero record; each says the new order where it
  describes it.
- A new design note, docs/design/latents-first.md, for the posterior-changing class: the question, the two outer
  steps, why the draw precedes the rescaling step, the gate, with a row in the design index.
- NEWS: the bullet above, and nothing under BUG FIXES (the two defects never reached a release).

## Gates

### New exact gate: outer-step-exact.R (slice 1)

benchmarks/R/outer-step-exact.R, on the real engine through the pointer-level entries as
[`mask-redraw-exact.R`](../../benchmarks/R/mask-redraw-exact.R) is. 32 rows in 4 cells of one ordinal predictor,
one tree, `n.cuts = 3`, the tree space enumerated under the CGM prior; an offset gamma * v with gamma on the grid
-1, -0.5, ..., 4 under a uniform prior, drawn every sweep by an outer step; the exact posterior of gamma and of
each cell's fit computed by enumeration over trees, quadrature over each leaf, and quadrature over k under its
prior. Arms:

- probit, k drawn (the rescaling step runs), outer step integrating z out (pnorm);
- probit, k drawn, outer step conditioning on z;
- logistic, outer step integrating omega out (plogis);
- logistic, outer step conditioning on omega.

Pass rule |z| <= 4.5 on every cell fit and every P(gamma) grid point whose exact mass is at least 1e-4 (8 of the 11;
the three left out hold under 3e-4 of the mass together, and a point visited that rarely has no normal batch error); quick mode 4 chains of 250,000 sweeps after 2,000 of
burn-in (the probe runs 20,000 sweeps a second). Which arms fail today is measured in slice 1, not assumed: the
slice runs every arm on a main build, on the prototype with the draw first, and on the prototype with the draw
after the rescaling step (mutant M2), and records each arm's largest |z| in the benchmarks README section. The
planning probes predict three failing arms on main, the first, third and fourth: probe 3's probit arm
integrating z out has P(gamma = 4) 0.142 against 0.112, about 15 standard errors at 400,000 sweeps; its logistic
arm integrating omega out, today's order, has the first cell's fit at -1.413 and -1.416 against -1.423 to -1.429
over seven draw-first runs and P(gamma = 4) at 0.3136 and 0.3148 against 0.3172 to 0.3189; the fourth arm is the
setOffset defect. The prototype is predicted to pass all four, and its M2 build to fail the first. The gate
checks where the draw sits (M1, M2) and nothing else: a stale working pointer, a missed invalidation and a split
pair (M4 to M6) leave its arms exact or out of its reach, and are the unit tests' job. Added to
.github/workflows/exact-gates.yaml in slice 2a.

### Tests (slice 2a: the order)

- tests/cpp, per augmenting family: snapshot the generator and the fits, run one sweep whose hook applies a
  setter (setOffset, setPredictor, setResponse, setActiveRows, setWeights where taken), and replay
  `refreshLatents` on a twin from the snapshot at the hook's state (taken after the setter returns); the sweep's
  latents equal the twin's bit for bit. Probit at a fixed k, so the rescaling step does not scale them; this is
  why the test cannot see M2.
- tests/cpp: under linear leaves on logistic, with leaves of at least the cache's minimum size (32 rows; smaller
  leaves rescan and never touch the cache), two or more sweeps draw bit for bit as a twin whose crossproduct
  cache is cleared before every tree update.
- tinytest, from the hook-vs-between probe: for probit, logistic, nbinom, Student-t and gaussian, setOffset and
  setActiveRows applied in the per-sweep hook (the internal `runWithCallback` entry) and between `run(0, 1)`
  calls give identical draws.
- tinytest: after each `run(0, 1)` on an ordinal sampler, every z lies in its category's interval under the
  thresholds that run recorded.
- The speed A/B of [Cost](#cost), recorded in the slice's landing note.

### Tests (slice 2b: the setters and the state)

- tests/cpp: every setter, setData, setModel's shape and setState under other weights and another censoring
  status consume no variate (generator state unchanged across the call).
- tests/cpp: after logistic setOffset (and setWeights) the working response equals the one rebuilt from omega,
  the response, the counts and the new offset, every row, inactive rows included.
- tinytest: a logistic sampler created with offset o2 and one created with o1 and handed o2 draw identically.
- tinytest: getLatents is unchanged across each setter on each family (aft setResponse: event rows the new times,
  censored rows the new bounds; setModel's new nbinom shape: omega kept); a Student-t setResponse keeps nu (the
  next run's resid.df draw starts from it) and the lambdas.
- tinytest: setData on an ordinal, a Student-t (resid.df drawn) and an nbinom (shape drawn) sampler keeps the
  cutpoints, nu and the shape (read from the state before and after the call), and an ordinal setData holding
  fewer categories than the sampler keeps all K - 1 cutpoints.
- tinytest: a weight change undone by setWeights and setState in either order continues bit for bit, on
  logistic and Student-t; on aft, a state installed under another censoring status draws nothing at the install.
- Retired or rewritten: test-state-zero-weight-rows.R; the redraw assertions in test-active-rows-reactivation.R,
  test-logistic-weight-swap.R and test-aft.R; the digest and reconciliation arms of test-state-weight-pairing.R and
  test-setstate-force-update.R; the files listed by `grep -li redraw inst/tinytest`; their tests/cpp twins in
  test_sampler.cpp and test_state.cpp, and test_facade.cpp's spies and digest checks for `weightsDigest`,
  `reapplyWeights`, `zeroWeightRows`, `survivalDigest` and `reapplySurvivalStatus`.

### Mutants the reviewer runs

Slice 2a's reviewer runs M1 to M6 and M10; slice 2b's runs M7 to M9 and M11.

| mutant | caught by |
|---|---|
| M1 draw after the trees (today) | outer-step-exact's failing arms as measured in slice 1; the replay test |
| M2 draw after the rescaling step | outer-step-exact's probit arm integrating z out, alone (the replay test runs at a fixed k, and both placements are valid scans of a plain fit, so no other gate sees it) |
| M3 draw before the per-sweep hook | hook-vs-between (setOffset arm); the replay test |
| M4 the vector-leaf invalidation after the draw deleted | the linear-leaf twin test |
| M5 ordinal cutpoints drawn at the end, z at the start (or nbinom r apart from omega) | the ordinal interval test; the replay test |
| M6 working pointers read once per run | hook-vs-between (setActiveRows arm; probe 1) |
| M7 logistic setOffset keeps the reshift | the tests/cpp working-response read alone |
| M8 a setter or setState keeps its redraw | the no-variate test; the getLatents-unchanged test |
| M9 Student-t setResponse cold-starts | the Student-t setResponse test |
| M10 one combined-fits call at the sweep top, its pointer reused by sigma, the glue and the variance forest | bcf-equivalence.R `--bitwise` on its gaussian scenarios against the baseline (a gaussian BCF's sigma reads the stale blend) |
| M11 setData restarts a kept family parameter | the setData tinytest |

Leaving the invalidation at the end of the sweep, the planned M4, is not a mutant: nothing fills the
crossproduct cache between that point and the next sweep's draw, so the next tree loop starts from an empty
cache either way and the draws are bitwise those of the correct order. The debug assertion at the end of each
sweep checks a constant leaf's total fits only, not this cache, so it is not a killer. M7 has no visible effect
end to end once the draw comes first: the draw rewrites every active row's working response before the trees
read it (the prototype is M7, and its created-versus-handed run is identical, probe 2); the created-versus-handed
test stays as a regression test against M1 with M7 together, which is main. M5 is a valid sampler for the fits
and cutpoints, so no exact gate on them sees it; only the joint read does.

### Battery and how to tell a legitimate move from a bug

The posterior-changing battery ([RNG classes and their gates](README.md#rng-classes-and-their-gates)), run at the
end of each of 2a and 2b: tests/cpp, the full tinytest suite, sanitizers on tests/cpp, R CMD check, every exact
gate in quick mode with the new one. The equivalence trio and test-reproducibility-binaryResponse.R are
re-recorded once, at the end of 2b (the snapshot regenerated by replaying the whole file,
tools/regenerate-snapshots.R); the reproducibility leg on the branch is red from 2a's tip until then, which the
merge-once rule allows. 2a records a scratch baseline of the trio at its tip, not committed, for 2b's partition
check. A move is legitimate when:

- at 2a, the moved set is exactly the augmenting families: `equivalence.R compare --bitwise` against main's
  baseline is identical on every gaussian, heteroscedastic gaussian, multinomial and uncensored-aft scenario,
  multinomial-equivalence.R identical throughout, bcf-equivalence.R identical on its gaussian scenarios;
- at 2b, nothing moves against 2a's scratch baseline but a scenario that calls a setter, setData or setState
  (plain fits draw nothing at a setter on either side);
- every moved scenario passes the statistical compare against main's baseline (plain fits are shifting);
- the other three snapshot files pass unchanged and only the binary one is regenerated;
- every exact gate passes, the new one included; the latent-family gates (ordinal, negbin, t, logistic-reference,
  aft, aft-hetero-pit, hazard, hurdle, probit-k-scale, bcf-latent, mask-redraw) are the oracle for the re-recorded
  baselines, as the MANIFEST's oracle rule asks.

A scenario that moves outside the augmenting families, or a latent family that stays bit for bit, is a bug.

### Consumer runs

revdep-smoke on the branch tip; stan4bart's tinytest, treatSens's, bartCause's and bairrtt's tests against the
built branch, re-recording bartCause's probit regression snapshots and any pinned treatSens value that moves.

## Calls made in planning

Recorded at landing in an A entry of the ledger:

1. growForestFromRoot takes the new order.
2. Every setter redraw goes (not only those the ruling names), and the Student-t zero-weight record, the weights
   and survival digests, their state attributes and setState's reconciliation retire with them. This follows from
   dec-B586 and is no new call, but it replaces the mechanism of two earlier entries: the A entry cites dec-A160
   (the Student-t redraws under weights and at a restore) and dec-B277 (the state's zero-weight record) as
   superseded in mechanism, their goal, an exact restore under other weights, kept and strengthened (setWeights
   and setState undo each other bit for bit in either order); the landing note to the maintainer says so in one
   line.
3. The working pointers are read every sweep, which fixes a mask installed from the per-sweep hook.
4. Logistic setOffset rebuilds its working response from omega and the new offset. The draw first already
   corrects the stale-offset sweep (it rewrites every active row's working response before the trees read it), so
   this is for the state between setOffset and the next sweep: the working response stays a function of omega,
   the response, the counts and the offset at every moment, an inactive row's included, for any reader between
   runs. A plain store of the offset was the alternative, the same size, leaving that state stale until the draw.
5. Student-t setResponse keeps lambda and nu. dec-B586 records the fix ("which keeping nu fixes"); keeping
   lambda as well is this plan's reading, harmless since the next sweep draws every lambda first; the landing's
   A entry records it as built.
6. Creation-time starts stay deterministic and unchanged.
7. aft setResponse leaves a censored row at its new bound until the next sweep (rather than keeping an old
   imputed time that may lie below the new bound).
8. setModel's new fixed nbinom shape keeps omega (the [getLatents](#getlatents) section says why); setData keeps
   starting the latents as creation does.
9. dec-B589 is built as "always keep": the ruling's restart when the number of categories changes cannot be
   reached, because setData refuses a response outside the categories the sampler was created with (the
   ordinal check in the bridge's response validation, "make a new sampler to change the category set"), and data
   holding fewer of them keep the sampler's K (critique probe ordinal-setdata.R). No restart branch is built; the
   A entry records the reading.
10. Step 2 forms the combined fits only for a family that draws, and sigma, the glue and the variance forest read
    them formed again after the trees.

## Ruled

### 1. Does setData keep the drawn ordinal cutpoints and Student-t degrees of freedom?

Ruled (a), dec-B589. The maintainer, 2026-11-03: "Yes, keep them." Built as "always keep" (call 9 of [Calls
made in planning](#calls-made-in-planning)): setData already refuses data whose categories differ from the
sampler's, so the restart clause of option (a) below cannot be reached.

Background: `$setData` replaces the whole data object. Today an ordinal sampler restarts its cutpoints at the
prior's default spacing and a Student-t sampler restarts nu at its grid's middle (and lambda at 1); the
negative-binomial shape restarts too. A caller who swaps data every sweep (an imputed predictor, say, written
through the data object rather than setPredictor) therefore forces those parameters back every sweep. Probe 4
(400 swaps of the same data on 150 rows): ordinal cutpoints average 0, 1.04, 2.13 against 0, 1.96, 3.81 without
swaps, under either order; nbinom's shape is unaffected once drawn first (3.240 with and without), its draw
integrating omega out and not reading the old shape; Student-t nu 8.45 against 8.65, inside the noise of that
probe. setResponse already keeps the cutpoints (ordinal.md). No ruling names setData.

Options:

- (a) Keep them, and the negative-binomial shape for one rule, restarting the cutpoints only when the number of
  categories changes: about 40 lines and 2.5 agent-hours in slice 2, draw-changing only for setData callers.
- (b) Keep today's restart and say so in setData's help: about 10 lines, 0.5 agent-hours; the swap loop stays
  biased.

Recommendation: (a). The next sweep draws the latents from the kept values first, so keeping them is a valid
state under any data, and a user who swaps data reads the same as one who swaps the response.

## Slices

| slice | tier | content | lines | agent-hours |
|---|---|---|---|---|
| 1 | gate-only (sonnet build, sonnet review) | outer-step-exact.R; benchmarks README section; each arm's largest \|z\| measured and recorded on main, on the prototype with the draw first and with it after the rescaling step (M2); the gate's reach stated (placement only) | ~400 | 4 build, 1.5 review, 0.5 fix |
| 2a | draws (opus build, opus review running M1 to M6 and M10) | sweep order, combined fits formed only where read, pointers every sweep, the invalidation after the draw, growForestFromRoot, the interrupt path; the replay, linear-leaf twin, hook-vs-between and ordinal interval tests; the gate in exact-gates.yaml; the speed A/B; battery, partition against main, a scratch baseline at the tip | ~650 | 8 build, 4 review, 2 fix |
| 2b | draws (opus build, opus review running M7 to M9 and M11) | setter redraws out, logistic setOffset and setWeights rebuilds, Student-t setResponse keeps, setData keeps (dec-B589), setModel's shape keeps omega, the zero-weight record, both digests and setState's reconciliation retired through model, chain, sampler, facade and bridge; their tests retired or rewritten; mask-redraw-exact.R's header; battery; snapshot and equivalence trio re-recorded once | ~1250 | 11 build, 5 review, 3 fix |
| 3 | shipped help (sonnet build, sonnet review against the engine) | every place under [Documentation](#documentation), recipe 6's numbers re-measured; codoc, doc-freshness | ~750 | 7 build, 2 review, 1 fix |
| 4 | orchestrator | consumer runs and re-records (bartCause, treatSens, stan4bart, bairrtt), revdep-smoke, records, landing | ~50 here, consumer repos apart | 4 |

About 53 agent-hours by slice, plus about 3 for rebasing onto coefficient-prior-live and re-gating what it
touches (probit-k-scale-exact among them): about 56, against the first plan's 34. Plans here have run 1.5 to 2
times low; the slices above are sized up from the first estimate, and the critique's 50 to 70 (34 at that rate)
brackets the total. Slices 2a, 2b and 3 land on the branch before it merges; docs land with code. Slices 2a and 2b
each run `R CMD INSTALL --preclean` (facade virtuals change in 2b; 2a's chain change reaches every family's
objects).

## Order and overlaps

- coefficient-prior-live (first): shares src/bartcore/chain.hpp (its hunks are live-coefficient installs, none in
  `runSweeps`), facade.hpp, sampler.hpp, R_interface_bartcore.cpp, R/dbarts.R, man/dbartsSampler-class.Rd,
  inst/NEWS.Rd, tests/cpp/test_facade.cpp and test_sampler.cpp, and touches probit-k-scale-exact.R. Lands first as
  the brief orders; this branch rebases onto it and reruns probit-k-scale-exact.
- latents-first (second).
- frequency-weights (third, for its setter slices): shares `ProbitResponse`, `OrdinalResponse`,
  `LogisticResponse`, `NBResponse` setWeights and mask composition in model.hpp, the setWeights, setData and
  setActiveRows help and the embedding page. Build it on this; three of its clauses go, and its slices get
  simpler:
  - its S4: a row the weight composition reactivates "redraws its latent as a mask change does now"; after this
    plan no mask change draws, so that clause and its tests go;
  - its S3: `$setCounts` on a binomial sampler "swaps both and redraws as a count swap does now"; after this plan
    a count swap stores the counts and the next sweep draws, so `$setCounts` does the same;
  - its design call 5 and plan mutant 5: the reading enters the weights digest so that a state installed under
    the other reading redraws; after this plan there is no digest and no reconciliation, so the call, the mutant
    and test-state-weight-pairing.R's new arm go (the next sweep draws under the reading in force).
- approximate-draws (plan, on its branch): its accept-or-keep route proposes "after the trees" and its setters
  draw "AT the setter" on a separate exact path. After this plan the route runs at step 2, and a row a setter or
  restore moved gets its exact draw there, at the next sweep, marked by the setter rather than drawn by it
  ([Accept-or-keep latents (approximate-draws)](#accept-or-keep-latents-approximate-draws)); its slice S2's
  separate setter path becomes that step-2 branch. Its plan is amended to say so before it builds.
- host-neutral-engine S2a onward (S1 landed; S2a bridge refusals as values, R_interface_bartcore.cpp throughout)
  overlaps only the state install region this plan trims (the digest and record reads go); S3 rewrites the same
  header sentences (the run's interrupt paragraph, the latent reader). Whichever lands second carries the other's
  text; this plan first, so S3 rewrites settled text.

Build order (the orchestrator's call): coefficient-prior-live, then latents-first (slices 1, 2a, 2b, 3, 4), then
frequency-weights, then host-neutral-engine S2a onward.

## Verification

    R CMD INSTALL --preclean .
    cd tests/cpp && make && ./test_bartcore
    Rscript -e 'tinytest::test_package("dbarts")'
    Rscript benchmarks/R/outer-step-exact.R quick        # OK on every arm
    for g in ordinal-exact negbin-exact t-exact logistic-reference aft-exact aft-hetero-pit hazard-exact \
      hurdle-exact probit-k-scale-exact bcf-latent-exact mask-redraw-exact; do Rscript benchmarks/R/$g.R quick; done
    # 2a: against main's baselines; 2b: against 2a's scratch baseline (identical but setter scenarios)
    Rscript benchmarks/R/equivalence.R compare <baseline> --bitwise   # fails exactly on augmenting scenarios (2a)
    Rscript benchmarks/R/equivalence.R compare <baseline>             # passes
    Rscript benchmarks/R/multinomial-equivalence.R compare <baseline> --bitwise # identical
    Rscript benchmarks/R/bcf-equivalence.R compare <baseline> --bitwise         # identical on gaussian scenarios
    Rscript tools/regenerate-snapshots.R    # end of 2b; only test-reproducibility-binaryResponse.R changes
    Rscript tools/check-rc-codoc.R . && Rscript tools/check-doc-freshness.R .

Slice 1 on its own: `outer-step-exact.R quick` against a main library and the prototype library at each
placement, its per-arm |z| recorded. Expected: everything passes on the branch; the mutants of [Mutants the
reviewer runs](#mutants-the-reviewer-runs) each fail their named test.

## Critique and amendments

A blind opus critique on main 7937d252f (scratch/ltf-crit/report.md, probes in scratch/ltf-crit/probes/) found
the design sound and every premise holding, and asked for ten amendments: BUILD AFTER AMENDMENTS. Each finding
was checked against the code at 0f6f344f5 before it was written in; all ten are applied.

### Amendments after the critique

| # | amendment | where it went |
|---|---|---|
| 1 | Mutant table: M7 caught only by a tests/cpp working-response read; M4 is "the invalidation deleted", the end of the sweep being no mutant; M2 caught by the gate's probit arm alone; M10 "sigma reuses the pre-tree fits" added | Mutants the reviewer runs; Tests (slice 2a), Tests (slice 2b) |
| 2 | The weights and survival digests, their methods, encoders, state attributes and tests retire with weights.zero; frequency-weights' call 5 goes | Setters; Tests (slice 2b); Order and overlaps |
| 3 | bairrtt in the consumer table and the consumer runs | Consumers; Consumer runs |
| 4 | dec-B589 built as "always keep", with the reason (setData refuses another category set); M11 and a setData test | Ruled; Calls made in planning, 9; Setters; Tests (slice 2b) |
| 5 | getLatents names the setData start; setModel's shape keeps omega (chosen over listing it: no reader needs the cold start, and the draws are the same); the setActiveRows rebuilds kept in the setter table | getLatents; Setters; Calls made in planning, 8 |
| 6 | Slice 1 measures and records each arm's largest \|z\| on main and on the prototype (and at the M2 placement), and states that the gate covers placement only | New exact gate; Slices |
| 7 | The landing's A entry cites dec-B277 and dec-A160 as superseded in mechanism, with one line to the maintainer | Calls made in planning, 2; Context |
| 8 | Combined fits formed only where read; M10; a same-machine speed A/B on probit, logistic, BCF probit and multinomial | Sweep order, steps 2 and 8; Cost; Calls made in planning, 10 |
| 9 | Vignette recipe 6's numbers re-measured; mask-redraw-exact.R's header; the approximate-draws and frequency-weights setCounts overlaps; the probit setResponse change from 0.9-34 in NEWS | Documentation; Order and overlaps; 0.9-34 and NEWS |
| 10 | Slice 2 split into 2a (the order) and 2b (the setters and the state), baselines re-recorded once at the end of 2b; slice 3 resized; about 56 agent-hours | Slices; Battery; the front block |

### Where this plan adds to the critique

- M4: the planned second killer, the debug cache assertion, checks only a constant leaf's total fits, so it
  cannot see either form of the mutant; the linear-leaf test must use leaves of at least the cache's minimum size
  and compare against a twin with the cache cleared.
- M10: the gaussian BCF compare sees the mutant only in the form "one combined-fits call at the sweep top for
  every family"; once step 2 forms the fits only for a family that draws, a gaussian BCF has no pre-tree blend to
  reuse.
- The setter table lacked setData and setModel's shape rows, and the documentation list lacked
  public-surface.md, which names the digests.
- No finding is declined.

### Open calls

None. The critique found no user-visible call the rulings leave open, and this pass confirms it: the setData
start and setModel's kept omega change only what getLatents shows before the next sweep, which draws the same
either way; the digests and record belong to development-build states only (dec-B461); "always keep" is dec-B589
read against a refusal the package already makes.

## Premises run

Prototype: a copy of main 730dd0adf with the draw in `runSweeps` movable by an environment variable (0 first, 1
after the rescaling step, 2 last as today) and the pointers read after the draw; setters unchanged. Built as a
private library; main built as another. Scripts and outputs in scratch/ltf-plan/probes/.

1. hook-vs-between.R: setOffset in the hook equals between runs on main for probit and logistic; setActiveRows in
   the hook differs from between runs on main from sweep 1 (probit 2.11, logistic 2.42, gaussian 0.094 max
   difference); identical on the prototype.
2. logistic-setoffset.R: created with o2 against created with o1 and handed o2: probit identical, logistic
   differs by 2.66 on main, identical on the prototype; nbinom differs by design (its count center is kept).
3. contractB-probit.R (probit and logistic, both outer steps, 400,000 sweeps, two seeds; logistic five): probit
   conditioning on z agrees under every placement; integrating z out agrees with it only with the draw first,
   after the rescaling step and last both off by about 0.03 in P(gamma = 4). Logistic integrating omega out with
   the draw first agrees with conditioning on it drawn first; conditioning on omega with today's order is off
   (P(gamma = 0.5) 0.0153 against 0.0172), which is probe 2's defect, not the order. The same output shows the
   logistic arm integrating omega out off under today's order as well (the critique's reading: first cell's fit
   -1.413 and -1.416 against -1.423 to -1.429, P(gamma = 4) 0.3136 and 0.3148 against 0.3172 to 0.3189), so
   slice 1 measures every arm rather than assuming which fail.
4. setdata-resets.R: as in ruled call 1.
5. Exact gates on the prototype, draw first, quick: ordinal, negbin, t, logistic-reference, aft, aft-hetero-pit,
   hazard, hurdle, probit-k-scale and bcf-latent all OK. Full tinytest suite: 40,650 tests, 0 failures (the
   snapshot files exit off the reference flag; no other test pins a latent family's draws).
6. neutral-partition.R: main and prototype identical for gaussian, heteroscedastic gaussian, multinomial and
   uncensored aft; different for probit, logistic, censored aft and Student-t.
7. Read: the sweep, each family's setters, the bridge's copy-on-set, growForestFromRoot, 0.9-34's bartFit.cpp
   (latents after the trees; setResponse redraws; setOffset keeps), stan4bart's init and loop, rbart_vi's loop,
   treatSens's treatment model update.

## Not verified

- tests/cpp under the prototype (not built; the prototype did not change setters).
- M3, M5, M8, M9, M10 and M11 are predicted from the code, not run. M7's invisibility end to end is shown by the
  prototype (probe 2); M4 as first planned being no mutant, and M2 having one killer, are read from the code.
- Which outer-step-exact arms fail on main: predicted (three), measured in slice 1.
- The speed A/B: not run; slice 2a.
- The consumer test suites against the prototype: the installed stan4bart, treatSens, bartCause and bairrtt are
  built against an older dbarts and were not rebuilt; bairrtt's pinned values were read, not run.
- Student-t nu under setData swaps: the probe's difference is inside its noise.
- BCF gaussian bit for bit, and the equivalence trio's partition, predicted from probe 6 and the code.
