# latents-first: a sweep draws its augmentation variables before the trees

Status: PLANNED 2026-11-03 on main 730dd0adf (dec-B586), not built. A blind critique follows.

agent: one opus implementer for slice 2 (engine) and one sonnet implementer for slices 1 and 3 (a gate, then
the documentation). Review by tier ([Process by risk](README.md#process-by-risk)): slice 1 one sonnet review
(gate-only); slice 2 one opus review that runs the mutants below, the same reviewer checking the fix round
unless the fixes change engine behavior; slice 3 one sonnet review against the engine as built. Slice 4 is the
orchestrator's.
rng: POSTERIOR-CHANGING overall. For a plain fit of a family with augmentation variables the draws move and the
posterior does not (shifting); for gaussian, heteroscedastic gaussian, multinomial and aft without censored rows
the draws stay bit for bit. The posterior changes where a caller changes the sampler between sweeps and its
outer step integrates the latents out (now exact, off before), and on a logistic sampler handed a new offset
(a defect fixed). Gates: [Gates](#gates).
window: before 1.0-0. After coefficient-prior-live lands on main; before frequency-weights' slices that touch
the setters and before host-neutral-engine's header rewrite ([Order and overlaps](#order-and-overlaps)).
budget: ~2400 lines in the tree: slice 1 ~380 (a new exact gate ~330, benchmarks README ~50); slice 2 ~1450
(engine ~110 added and ~330 removed, bridge ~70 removed, tests/cpp ~380 added and ~200 removed, tinytest ~260
added and ~400 rewritten or removed, one snapshot file and three equivalence baselines re-recorded); slice 3
~520 (help ~200, vignette ~40, header ~30, embedding page ~60, architecture and design notes ~150, a new design
note ~60, NEWS ~10); records ~40. About 34 agent-hours in all ([Slices](#slices)).

## Summary for the maintainer

Every family that augments its likelihood - the probit and ordinal latents, the Polya-Gamma precisions of
logistic, binomial and negative binomial, the Student-t scale mixing, the aft imputed times - draws those values
once a sweep. Today it draws them after the trees. So when the code that runs the BART sampler inside a larger
sampler changes something between sweeps (an offset, a predictor, sigma), the next sweep's trees fit values
drawn before the change. That is right only when the larger sampler's own step reads those values; when its
step averages them out, as treatSens's treatment model and the vignette's fourth recipe do, the answer is off.
Drawing them first, right after the point where the larger sampler makes its changes, is right either way,
costs nothing extra a sweep, and lets every setter stop drawing.

What changes for a user:

- Every family with augmentation variables draws differently at the same seed, once. Gaussian, multinomial,
  heteroscedastic gaussian and aft without censored rows draw exactly as before (checked).
- No setter draws any more: setResponse, setActiveRows, setWeights and setState's reconciliation store what they
  are given and the next sweep draws the rest. getLatents between a setter and the next run returns the values
  from before the setter; after a run, the values the last sweep's trees were drawn from.
- A Student-t setResponse keeps the drawn degrees of freedom and scales, where today it resets them.
- On a small probit model whose larger sampler averages the latent out, today's order is off by about 15
  standard errors; drawing first is exact (checked on the real engine).
- Two defects found while planning are fixed with it: a logistic sampler handed a new offset kept fitting the
  old one for a sweep (it reads the old offset through a buffer the new one has already overwritten), and a mask
  installed from the per-sweep hook was ignored by that sweep.
- Documentation states the rule in each place a user or package author meets it: getLatents and every setter,
  the help page and vignette on running the sampler inside a larger one, the C header.

One open call: whether swapping in a data object keeps the drawn ordinal cutpoints and Student-t degrees of
freedom when nothing about the categories changes. Today it restarts them, which biases a loop that swaps data
every sweep; the recommendation is to keep them.

0.9-34 drew probit's latent after the trees and kept it across setOffset, setPredictor and setWeights, so for
probit this is the same behavior made correct, not a regression; the other families are new in 1.0-0.

## Goal

After this lands, each sweep of every sampler draws its augmentation variables (with the ordinal cutpoints and
the negative-binomial shape drawn beside them) after the per-sweep hook and before anything reads them, so a
change made between sweeps or inside the hook never leaves the trees a stale latent. No setter draws. The
sampler is exact whether an outer step conditions on the latents or integrates them out, and the contract is
written where users and package authors read it.

## Context

- Ruling: dec-B586 (direction, cost, documentation places, Student-t nu); related dec-B583 (accept-or-keep
  Polya-Gamma latent), dec-B580, dec-B581 (draw callback). Evidence: scratch/stc-chk/report.md and its probes.
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
- Branches: coefficient-prior-live (under review), frequency-weights (plan 516d02e40), host-neutral-engine
  (plan dd24d7670), approximate-draws (plan).

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
   then omega, lambda then nu).
3. The working response and working weights pointers read from the family, every sweep, then the vector-leaf
   statistics invalidated where the family's weights vary per sweep.
4. Mean weights under a variance forest (they replace the weights pointer for the trees).
5. The level shift (fit-preserving; reads no latent).
6. The probit rescaling step (reads z, the offset and the fits; scales z, the leaves and k).
7. The trees and total fits; the per-forest glue and multinomial latents inside the loop, unchanged.
8. Sigma; the combiner's glue and post-combine move; the variance forest; k and DART.
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
not change, except that the Student-t record of zero-weight rows is retired (below); a development-build state
carrying it loads and the attribute is ignored (only shipped versions are owed compatibility; no migration).

### getLatents

After a run: the values the last sweep's trees were drawn from. Between a setter and the next run: the values
from before the setter, with one exception, aft setResponse, where a censored row holds its new bound (data)
until the next sweep draws it. In the per-sweep hook: the previous sweep's values. In the draw callback: this
sweep's, the ones its trees were drawn from (today the callback saw the values the next sweep's trees would fit).

### Setters

Each setter stores what it is given, rebuilds what is deterministic, and draws nothing:

| setter | today | after |
|---|---|---|
| setResponse, probit / ordinal / logistic / nbinom | redraws every active row | stores y; ordinal recomputes its scales and keeps the cutpoints; nbinom keeps r and recomputes its kernel; the working response rebuilt from the latents in force |
| setResponse, aft | redraws censored rows | installs the times; censored rows at their new bounds |
| setResponse, Student-t | lambda reset to 1, an estimated nu to its grid's middle | lambda and nu kept |
| setOffset, logistic | shifts the working response by old minus new offset, reading the old one through the overwritten buffer | rebuilds the working response from omega and the new offset |
| setOffset, others; setPredictor; setSigma; setForestBasis | no draw | unchanged |
| setActiveRows | redraws the rows it switches back on (Student-t: lambda) | stores the mask; Student-t recomputes its membership marks |
| setWeights, logistic | inactive rows cold-started, then every active row redrawn | stores the counts |
| setWeights, Student-t | redraws rows whose weight turns positive | stores the weights; membership marks and the composite recomposed |
| setSurvivalStatus | no draw | unchanged |
| setState reconciliation | Student-t rows redrawn by the stored zero-weight record; aft censored rows redrawn on a status change | none; the next sweep draws them |
| sampleTreesFromPrior | (comment: follow with setResponse for a true restart) | the next sweep draws the latents from the new trees |

Retired with the redraws: `redrawEntering`, `reapplyWeights`, `zeroWeightRows`, `reapplySurvivalStatus` (the
model members, the facade and sampler virtuals, their tests), the bridge's `weights.zero` state attribute and
the setState reconciliation that reads it. The Student-t record existed because a setter that drew needed to
know which rows were out (dec-B277); with every lambda drawn at the next sweep there is nothing to reconcile.
`findReactivatedRows` stays only if another caller needs it.

### Draw callback and the C header

The draw callback (dec-B581) fires at the record, unchanged; `dbarts_draw` carries no latents. The header's latent
reader [`dbarts_sampler_getLatents`](../../inst/include/dbarts/dbarts.h) gains the sentence of the getLatents
section above. The run's interrupt paragraph changes: a stop inside the latent draws, which now open a sweep,
leaves no tree of that sweep drawn and each latent this sweep's draw or the previous one (a sampled
negative-binomial shape's refresh put back whole), a valid state to run on from. The flat C API has no sweep
hook, so a C caller makes its changes between runs; nothing else in the header moves, and no signature changes.

### Cost

None added per sweep: the draw moves. Saved: the setter draws, which a caller swapping the response every sweep
pays today (about 22 percent of a sweep for probit and 50 percent for logistic at n = 5,000 and 50 trees,
scratch/stc-chk/report.md), and logistic setWeights' redraw.

### RNG

Draw-changing for every family with augmentation variables, with or without setters. Bit for bit for gaussian,
heteroscedastic gaussian, multinomial, aft with no censored row (probe 6); BCF gaussian and hurdle's positive
part are predicted to be, having no latent. Posterior-changing for embedded use whose outer step integrates the latents out and for logistic
setOffset; shifting for a plain fit.

### Accept-or-keep latents (approximate-draws)

Not built here. The constraint the approximate-draws plan must meet: an accept-or-keep Polya-Gamma latent is
exact for an outer step that integrates it out only where it is a fresh draw from its current conditional, so
every row whose conditional a setter or a restore moved since the last sweep (every row for setOffset, setSigma,
setForestBasis, setResponse and setState; the moved rows for setPredictor; the changed rows for setWeights and
setActiveRows; rows changed inside the per-sweep hook included) is drawn exactly at step 2 of the next sweep, and
the mark is cleared only after that draw. With this plan the setters stay draw-free; they set marks.

### Consumers

| consumer | outer step | effect | to re-record |
|---|---|---|---|
| stan4bart | conditions on z (contract A): setOffset, run(0, 1), getLatents | draws move, exact before and after | nothing pinned: its binary test checks a tolerance, its seed tests compare fits with themselves |
| rbart_vi (dbarts) | conditions on z | draws move; init reads the creation-time z, unchanged | nothing pinned moved under the prototype |
| treatSens, trt.model = bart() | integrates z out: setOffset then run | now exact; draws move | its tests' pinned values, if any move (not run, below) |
| treatSens outcome model | gaussian setResponse | bit for bit | nothing |
| bartCause | BCF and probit fits through run; grouped through stan4bart | draws move | the probit snapshots in its regression tests |

### 0.9-34 and NEWS

0.9-34 drew probit's latents after the trees, redrew them at setResponse and kept them at setOffset, setPredictor
and setWeights; it had no other augmenting family. NEWS, USER-VISIBLE CHANGES of 1.0-0:

> A sampler draws its augmentation variables (the probit and ordinal latents, the Polya-Gamma precisions, the
> Student-t scales and aft's imputed times) at the start of each sweep, before its trees, and no setter draws
> them. A larger sampler that changes the offset, predictors or sigma between sweeps is exact whether its own
> step reads `$getLatents()` or integrates the latents out; in 0.9-34 the second was not, the next sweep's
> trees fitting latents drawn before the change.

## Documentation

Shipped text states the contract; history stays in this plan and the ledger.

- [`dbartsSampler$getLatents`](../../man/dbartsSampler-class.Rd) (R/dbarts.R docstring and the Rd value section):
  "The values the last sweep's trees were drawn from. Each sweep draws them first, after the per-sweep hook and
  before its trees, so between a setter and the next run they are the values from before the setter (an aft
  setResponse leaves a censored row at its new bound) and the next sweep draws them afresh."
- setResponse, setOffset, setPredictor, setSigma, setForestBasis, setWeights, setActiveRows, setState (docstrings
  and Rd): one shared sentence, "This draws nothing; the next sweep draws the augmentation variables against the
  change." Removed: the Student-t redraw under weights, the logistic count-swap redraw and inactive cold start,
  the reactivation redraw and the joint-draw advice under the mask, the aft status-change redraw, the setState
  reconciliation and the weights.zero record; the undo paragraph becomes "a weight change is undone by setWeights
  and setState in either order, bit for bit". The Student-t setResponse text says nu and the scales are kept.
- The predict-then-setPredictor advice: adds that the accepted predictors reach the next sweep's latents, so a
  Metropolis step that scores proposals with the latents integrated out is exact.
- [dbarts-embedding.Rd](../../man/dbarts-embedding.Rd): a new subsection "When the augmentation variables are
  drawn": the order, the two kinds of outer step (one that reads getLatents and one that integrates them out) and
  that both are exact; "Which fit to read back" says what getLatents returns between sweeps.
- Vignette recipe 4: one paragraph, that scoring with the latent integrated out (pnorm) is exact because the next
  sweep draws sB's latents under the accepted theta. Recipe 6: one sentence, that z and f read after a run are a
  joint draw, z being the value f was drawn from.
- [`dbarts_sampler_getLatents`](../../inst/include/dbarts/dbarts.h) and the run's interrupt paragraph: as in
  [Draw callback and the C header](#draw-callback-and-the-c-header).
- [One sweep](../architecture.md#one-sweep): the step list rewritten in the new order; the RNG paragraph drops
  "the probit latent redraw on setResponse"; [Callbacks and errors](../architecture.md#callbacks-and-errors) says
  what the hook and the callback see.
- Design notes that state the old order or a setter draw: ordinal.md (draw order, setResponse),
  negative-binomial.md (r first, setResponse), robust-errors.md (setResponse cold init, zero-weight rows),
  active-rows-mask.md (reactivation redraws), aft-status-setter.md, aft-variance-forest.md, weighted-logistic.md,
  probit-k-scale-move.md (the rescale reads the fresh z), bart-as-a-component.md (the contract), survival.md,
  empty-leaf-veto.md where it names the order; each says the new order where it describes it.
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

Pass rule |z| <= 4.5 on every P(gamma) cell and cell fit; quick mode 4 chains of 250,000 sweeps after 2,000 of
burn-in (the probe runs 20,000 sweeps a second). On main the first and fourth arms fail (probe 3: P(gamma = 4)
0.142 against 0.112, about 15 standard errors at 400,000 sweeps); the prototype passes all four. Added to
.github/workflows/exact-gates.yaml in slice 2.

### Tests (slice 2)

- tests/cpp, per augmenting family: snapshot the generator and the fits, run one sweep whose hook applies a
  setter (setOffset, setPredictor, setResponse, setActiveRows, setWeights where taken), and replay
  `refreshLatents` on a twin from the snapshot at the hook's state; the sweep's latents equal the twin's bit for
  bit. Probit at a fixed k, so the rescaling step does not scale them.
- tests/cpp: every setter consumes no variate (generator state unchanged across the call).
- tests/cpp: under linear leaves on logistic, the leaf caches match the leaves after the draw (the invalidation).
- tinytest, from the hook-vs-between probe: for probit, logistic, nbinom, Student-t and gaussian, setOffset and
  setActiveRows applied in the per-sweep hook (the internal `runWithCallback` entry) and between `run(0, 1)`
  calls give identical draws.
- tinytest: a logistic sampler created with offset o2 and one created with o1 and handed o2 draw identically.
- tinytest: after each `run(0, 1)` on an ordinal sampler, every z lies in its category's interval under the
  thresholds that run recorded.
- tinytest: getLatents is unchanged across each setter on each family (aft setResponse: event rows the new times,
  censored rows the new bounds); a Student-t setResponse keeps nu (the next run's resid.df draw starts from it)
  and the lambdas.
- Retired or rewritten: test-state-zero-weight-rows.R and the redraw assertions in test-active-rows-reactivation.R,
  test-logistic-weight-swap.R, test-aft.R, test-state-weight-pairing.R, test-setstate-force-update.R and the files
  listed by `grep -li redraw inst/tinytest`; their tests/cpp twins in test_sampler.cpp, test_state.cpp and
  test_facade.cpp.

### Mutants the reviewer runs

| mutant | caught by |
|---|---|
| M1 draw after the trees (today) | outer-step-exact probit-integrating and logistic arms; the replay test |
| M2 draw after the rescaling step | outer-step-exact probit-integrating arm (probe 3); the replay test |
| M3 draw before the per-sweep hook | hook-vs-between (setOffset arm); the replay test |
| M4 vector-leaf invalidation left at the end of the sweep | the linear-leaf cache test; the debug cache assertion |
| M5 ordinal cutpoints drawn at the end, z at the start (or nbinom r apart from omega) | the ordinal interval test; the replay test |
| M6 working pointers read once per run | hook-vs-between (setActiveRows arm; probe 1) |
| M7 logistic setOffset keeps the reshift | the created-versus-handed offset test (probe 2); outer-step-exact logistic arms |
| M8 a setter keeps its redraw | the no-variate test; the getLatents-unchanged test |
| M9 Student-t setResponse cold-starts | the Student-t setResponse test |

M5 is a valid sampler for the fits and cutpoints, so no exact gate on them sees it; only the joint read does.

### Battery and how to tell a legitimate move from a bug

The posterior-changing battery ([RNG classes and their gates](README.md#rng-classes-and-their-gates)): tests/cpp,
the full tinytest suite, sanitizers on tests/cpp, R CMD check, every exact gate in quick mode with the new one,
the equivalence trio re-recorded, test-reproducibility-binaryResponse.R regenerated by replaying the whole file
(tools/regenerate-snapshots.R). A move is legitimate when:

- the moved set is exactly the augmenting families: `equivalence.R compare --bitwise` against the previous
  baseline is identical on every gaussian, heteroscedastic gaussian, multinomial and uncensored-aft scenario,
  multinomial-equivalence.R identical throughout, bcf-equivalence.R identical on its gaussian scenarios;
- every moved scenario passes the statistical compare against the previous baseline (plain fits are shifting);
- the other three snapshot files pass unchanged and only the binary one is regenerated;
- every exact gate passes, the new one included; the latent-family gates (ordinal, negbin, t, logistic-reference,
  aft, aft-hetero-pit, hazard, hurdle, probit-k-scale, bcf-latent, mask-redraw) are the oracle for the re-recorded
  baselines, as the MANIFEST's oracle rule asks.

A scenario that moves outside the augmenting families, or a latent family that stays bit for bit, is a bug.

### Consumer runs

revdep-smoke on the branch tip; stan4bart's tinytest, treatSens's and bartCause's tests against the built
branch, re-recording bartCause's probit regression snapshots and any pinned treatSens value that moves.

## Calls made in planning

Recorded at landing in an A entry of the ledger:

1. growForestFromRoot takes the new order.
2. Every setter redraw goes (not only those the ruling names), and the Student-t zero-weight record, its state
   attribute and setState's reconciliation retire with them.
3. The working pointers are read every sweep, which fixes a mask installed from the per-sweep hook.
4. Logistic setOffset rebuilds its working response from omega and the new offset, fixing the stale-offset sweep.
5. Student-t setResponse keeps lambda and nu. dec-B586 records the fix ("which keeping nu fixes"); keeping
   lambda as well is this plan's reading, harmless since the next sweep draws every lambda first; the landing's
   A entry records it as built.
6. Creation-time starts stay deterministic and unchanged.
7. aft setResponse leaves a censored row at its new bound until the next sweep (rather than keeping an old
   imputed time that may lie below the new bound).

## Open calls

### 1. Does setData keep the drawn ordinal cutpoints and Student-t degrees of freedom?

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
| 1 | gate-only (sonnet build, sonnet review) | outer-step-exact.R; benchmarks README section; shown failing on main, passing on the prototype | ~380 | 3 build, 1.5 review |
| 2 | draws (opus build, opus review with mutants) | sweep order, growForestFromRoot, pointers, setter redraws out, logistic setOffset, Student-t setResponse, the zero-weight record retired, open call 1 if (a); tests/cpp and tinytest; exact-gates list; snapshot and equivalence trio re-recorded | ~1450 | 13 build, 5 review, 2 fix |
| 3 | shipped help (sonnet build, sonnet review against the engine) | every place under [Documentation](#documentation); codoc, doc-freshness | ~520 | 5 build, 1.5 review |
| 4 | orchestrator | consumer runs and re-records (bartCause, treatSens), revdep-smoke, records, landing | ~40 here, consumer repos apart | 3 |

About 34 agent-hours. Slices 2 and 3 land on the branch before it merges; docs land with code. Slice 2 runs
`R CMD INSTALL --preclean` (facade virtuals change).

## Order and overlaps

- coefficient-prior-live (first): shares src/bartcore/chain.hpp (its hunks are live-coefficient installs, none in
  `runSweeps`), facade.hpp, sampler.hpp, R_interface_bartcore.cpp, R/dbarts.R, man/dbartsSampler-class.Rd,
  inst/NEWS.Rd, tests/cpp/test_facade.cpp and test_sampler.cpp, and touches probit-k-scale-exact.R. Lands first as
  the brief orders; this branch rebases onto it and reruns probit-k-scale-exact.
- latents-first (second).
- frequency-weights (third, for its setter slices): shares `ProbitResponse`, `OrdinalResponse`,
  `LogisticResponse`, `NBResponse` setWeights and mask composition in model.hpp, the setWeights, setData and
  setActiveRows help and the embedding page. Its slice 4 says a row the weight composition reactivates "redraws its
  latent as a mask change does now"; after this plan no mask change draws, so that clause and its tests go, and
  its slices get simpler. Build it on this.
- host-neutral-engine S2a (bridge refusals as values, R_interface_bartcore.cpp throughout) overlaps only the state
  install region this plan trims; S3 rewrites the same header sentences (the run's interrupt paragraph, the latent
  reader). Whichever lands second carries the other's text; this plan first, so S3 rewrites settled text.

Build order: coefficient-prior-live, latents-first slices 1 to 4, frequency-weights, host-neutral-engine.

## Verification

    R CMD INSTALL --preclean .
    cd tests/cpp && make && ./test_bartcore
    Rscript -e 'tinytest::test_package("dbarts")'
    Rscript benchmarks/R/outer-step-exact.R quick        # OK on every arm
    for g in ordinal-exact negbin-exact t-exact logistic-reference aft-exact aft-hetero-pit hazard-exact \
      hurdle-exact probit-k-scale-exact bcf-latent-exact mask-redraw-exact; do Rscript benchmarks/R/$g.R quick; done
    Rscript benchmarks/R/equivalence.R compare <previous baseline> --bitwise   # fails exactly on augmenting scenarios
    Rscript benchmarks/R/equivalence.R compare <previous baseline>             # passes
    Rscript benchmarks/R/multinomial-equivalence.R compare <baseline> --bitwise # identical
    Rscript tools/regenerate-snapshots.R    # only test-reproducibility-binaryResponse.R changes
    Rscript tools/check-rc-codoc.R . && Rscript tools/check-doc-freshness.R .

Expected: everything passes; the mutants of [Mutants the reviewer runs](#mutants-the-reviewer-runs) each fail
their named test.

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
   (P(gamma = 0.5) 0.0153 against 0.0172), which is probe 2's defect, not the order.
4. setdata-resets.R: as in open call 1.
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
- M3, M4, M5, M8 and M9 are predicted from the code, not run.
- The consumer test suites against the prototype: the installed stan4bart, treatSens and bartCause are built
  against an older dbarts and were not rebuilt.
- Student-t nu under setData swaps: the probe's difference is inside its noise.
- BCF gaussian bit for bit, and the equivalence trio's partition, predicted from probe 6 and the code.
