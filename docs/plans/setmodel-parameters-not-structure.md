# setmodel-parameters-not-structure: setModel changes a model's parameters and refuses its structure

Status: PLANNED 2026-10-29 on e2f37f6d (dec-B254); amended the same day after its blind critique (build after
amendments; scratch/crit-2plans/report.md). Two calls are open for the maintainer ([Open calls](#open-calls)):
the door for the tree priors of forests past the first, and whether an unstated variance forest follows
`$setModel`'s tree prior. Slice 4 waits on both and on forests-shared-total; slices 1 to 3 wait on neither.

agent: a sonnet implementer for slice 1 (R only) and one opus implementer for slices 2 to 4 (engine). Review by
tier ([Process by risk](README.md#process-by-risk)): slice 1 one sonnet review; slices 2, 3 and 4 one opus review
each that runs the mutants below, the same reviewer checking the fix round unless the fixes change engine
behavior. The blind opus critique of this plan has run.
rng: by slice. Slice 1 NEUTRAL, bit for bit: refusals and R checks only. Slice 2 POSTERIOR-CHANGING for a fit
that states a zero split probability, and bitwise neutral for every other fit. Slices 3 and 4 NEUTRAL for every
call main honors; they make honored what main accepts and ignores (a Student-t df, a negative-binomial shape, a
lengthscale), which has no main draw to compare against, and add what main refuses. Under open call 2's
recommendation one call main honors moves in slice 4: `setModel` changing the tree prior of a heteroscedastic
sampler whose variance forest stated none, posterior-changing for that call. How each class is gated:
[Tests and gates](#tests-and-gates).
window: before 1.0-0, after the merge to main (dec-B254: "Build it for 1.0. It can be post main-merged."). Slice 4
lands after forests-shared-total and is planned against that plan's design (its D7, D8 and D10).
budget: ~3450 lines in the tree: slice 1 ~630 (R ~240, help and NEWS ~60, tinytest ~330), slice 2 ~620 (engine
~50, bridge ~40, tests/cpp ~110, tinytest ~150, a new exact gate ~230, design note and NEWS ~40), slice 3 ~1150
(engine ~300, bridge ~160, R ~90, help ~50, tests/cpp ~250, tinytest ~300), slice 4 ~1000 (engine ~120, bridge
~170, R ~250, help ~70, tests/cpp ~100, tinytest ~290); records ~50. Slice 3's engine was ~130 before the
critique; four installs, each with its fixed and drawn handoffs, a check-first `Sampler` method and a facade
virtual, the DART grid rebuild and the kernel flush, do not fit there. Out of the tree, ~400 lines of scratch
A/B script.

## Summary for the maintainer

`setModel` is to change the numbers inside a model and refuse a model of another shape, by name, before
anything is stored (dec-B254). Today it does neither fully: a model of another family, a Student-t change, an
interaction limit, blocks and Gaussian-process settings are taken and quietly not used, so the live sampler and
its copy run different models; some numbers it should take, it cannot (the Student-t degrees of freedom, a fixed
negative-binomial shape, a Gaussian-process lengthscale, the tree prior of any forest but the first); and a DART
sampler refuses everything.

Of the five defects the item names, two are already fixed on main: the k column of `run()` after k is switched
between fixed and drawn (since 2026-10-11), and the move mixture reaching the variance forest (since
2026-10-11, kept by the per-forest move settings that landed today). The other three are as described: a zero
split probability on a column the trees use is taken and lingers for hundreds of sweeps, and a zero does not
exclude a column once the columns with positive probability run out along a path; editing `leaf.prior@k` alone
is taken and changes nothing; a residual prior on probit is stored and ignored.

The work is four slices. The first refuses every structural difference and the two record defects, R only. The
second makes a zero split probability mean the column is never split on, which changes results for fits that
state a zero, 0.9-34's included, and refuses a zero on a column in use. The third installs the numbers the model
carries; a number switched from fixed to drawn starts where it stands rather than where a new sampler would
start it. The first three come first. The fourth comes after the item that lets a model hold
several forests with no basis, and is planned against it: a tree prior given to `setModel` reaches every such
forest that did not state its own, and the size those forests share is changed through `setModel` as one forest's
is.

The fourth needs two answers from you. First, where a caller states the tree prior of the second forest or the
variance forest on a live sampler, since today no model object carries them; I recommend a `$setTreePrior()`
beside `$setLeafPrior()`. Second, whether a variance forest that never stated its own tree prior follows a tree
prior changed through `setModel`. Today it keeps the one it was created with, though at creation it took the
model's. I recommend it follow, as forest settings left unstated already follow a changed control; the cost is
that one call that works today changes what it does.

## Goal

`$setModel` refuses, by name and before anything is stored, a model that differs from the sampler's in any fact
dec-B254 calls structure, and installs every parameter it names on every chain, so that the live sampler, its
copy and its reload run the same model, and a model installed before any run draws what creating with it draws,
a parameter handed from fixed to drawn excepted by name (call 8).
A zero split probability excludes its column on every path, and the tree prior of every forest can be changed.

## What exists today

Run on e2f37f6d (shipped build, arm64 macOS; library scratch/libs/smp-plan, probes scratch/smp-plan, the
2026-10-05 survey's probes rerun from scratch/smp-plan/old). "Three routes" is the live sampler, `copy()` and
`storeState` plus `saveRDS` and `readRDS`.

- Where the model's facts live. Slots of the model: the tree prior (`cgm()` or `dart()`, with split
  probabilities and DART's settings), the leaf prior (kind, columns, a gp's lengthscale and `max.leaf.size`, and
  the raw `k` and `prior.sd` it was written with), the resolved leaf hyperprior and `prior.scale`, the residual
  prior, `leaf.scale` and `family`. Attributes of the model: `family.spec`, `resid.df`, `monotone`,
  `monotone.prior`, `interaction.max.order`, `interaction.forbidden`, `block.of.column`, `block.tree.counts`,
  `forest.columns`, `response.range`. Attributes of the control: `bartcore.forests` (every forest's tree count,
  base, power, sd, vars, labels and more), `bartcore.variance` (the variance forest's tree count, base, power,
  columns), `bartcore.shape`, `bartcore.n.categories`, `bartcore.hazard.periods`, `bartcore.survival` and
  `bartcore.moves`. dec-B239 ruled the five model facts off the control and onto the model (2026-10-05); that
  migration is not built and has no TODO item.
- How `setModel` works. [`dbartsSampler$setModel`](../../R/dbarts.R) refuses a multinomial sampler
  ([`refuseCountsMutation`](../../R/bartcore.R)), one whose forests carry amplitudes
  ([`refuseAmplitudeMutation`](../../R/bartcore.R)), a DART prior on either side, and a changed `forest.columns`;
  checks the split probabilities ([`refuseNoSplittableColumn`](../../R/model.R)); overwrites the new model's
  `family` and `response.range` with the sampler's; stores the model and calls
  [`bartcore_setModel`](../../src/R_interface_bartcore.cpp), restoring the old model on an error. The bridge
  refuses a monotone difference ([`monotoneModelChangeRefusal`](../../src/R_interface_bartcore.cpp)) and a leaf
  covariate difference, and passes base, power, split probabilities, the control's move mixture, the leaf scale,
  k or its hyperprior and the residual prior to [`Chain::setModel`](../../src/bartcore/chain.hpp), which writes
  them on the first forest, the mixture and the residual calibration on the variance forest too. Everything else
  the model carries is stored and never reaches the engine; a copy re-creates from it
  ([`recreatePointer`](../../R/dbarts.R)).
- Structure accepted and not honored, as the item says. A gaussian sampler given a probit sampler's model: taken,
  `family` left "gaussian" and `family.spec` stored as probit, sigma held at 1 and `k.scale` 18.35 on all three
  routes. Student-t: added to or removed from a gaussian sampler, taken, the live sampler unchanged and the copy
  and reload refused ("state is not consistent with this sampler"); df 3 to 30, 3 to estimated and estimated to 3,
  live keeps its own while copy and reload run the model's (live 3, copy 30). Interactions: `max.order = 1` added,
  removed or loosened, and a forbidden pair added, taken and ignored live; copy and reload run the stored model
  (live 3 columns in a tree, copy 1). Blocks added, removed or re-counted: likewise (block groups in a tree live
  2, copy 1). A gp lengthscale 1 to 0.02: live roughness 0.11, copy and reload 0.92; `max.leaf.size` 2: copy and
  reload fall back on 86% of evaluations, live on none; with kept draws a lengthscale change makes copy and
  reload fail. A negative binomial with shape 4 given a shape-10 model: taken, `family.spec` says 10, the engine
  and every route keep 4, the control's `bartcore.shape` holding 4.
- Refused today: DART on either side, its own model included ("changing a DART tree prior is not available on an
  existing sampler"), so neither power nor `alpha` can change; a sampler of several forests, its own model
  included; a multinomial sampler; the leaf kind and columns; monotone; `forest.columns`. A DART sampler takes
  `$setLeafPrior` and, since dec-A239, `$setControl`. Called directly, the bridge takes power, base, k and the
  residual prior under DART (k 5, leaves per tree 2.34 to 1.24, the split weights still drawn).
- The setters the item names do not exist. `$setLeafPrior(forests = list(forest(power = 3)))` is refused ("'power'
  is fixed at creation"). The variance forest's base and power are resolved into `bartcore.variance` at creation
  and stay there: after `setModel` with power 4 and base 0.5 on a heteroscedastic sampler, the record still says
  2 and 0.95 (TODO variance-forest-inherited-tree-prior). The engine has the Student-t and negative-binomial
  pieces below the facade ([`TResponse::restoreResidualDf`](../../src/bartcore/model.hpp),
  [`NBResponse::restoreShape`](../../src/bartcore/model.hpp)), the gp kernel only at creation and state install
  ([`GPGaussianLeaf`](../../src/bartcore/model.hpp), [`Chain::lengthscaleStateFeasible`](../../src/bartcore/chain.hpp)),
  and DART's grid only at creation ([`DartPrior`](../../src/bartcore/model.hpp)).
- Defect 1, the k column: fixed on main by 1e4dd72b (2026-10-11), which records k in `run()` whenever the first
  forest draws it at that run. Run: `setModel` from a fixed k of 2 to `chi(1.5, 2)` after 20 sweeps, live
  `run()$k` 1.63 1.68 1.56 ..., the copy draws too; from `chi()` to a fixed 3, `NULL` on both; created fixed,
  `NULL`. 0.9-34 gives the same: draws after fixed to drawn, `NULL` after drawn to fixed. The pin covers
  `$setLeafPrior` only (["sampler$setLeafPrior(normal(k = 3))"](../../inst/tinytest/test-run-k-after-setleafprior.R));
  `$setModel` is unpinned. The flat C API's [`dbarts_sampler_kIsSampled`](../../src/C_interface.cpp) reports the
  creation option, and the API has no prior setter, so stan4bart, which reads it once at creation, cannot see a
  switch.
- Defect 2, a zero split probability: present. `setModel` with `split.probs` of `c(0, 0.5, 0.5)` after 50 sweeps
  on trees splitting on the first column: taken; 35 such splits at the call, 1 at sweeps 200 and 400, none from
  600. The run-out half is pinned as current behavior
  (["zeroed <- runOut(tree.prior = cgm(split.probs = c(0.5, 0.5, 0)))"](../../inst/tinytest/test-single-forest-vars.R),
  182 of 182 passing on this build), and the numbers are in
  [Why not zeros in the split probabilities](../design/single-forest-vars.md#why-not-zeros-in-the-split-probabilities).
  The cause: where every column available at a node has probability zero,
  [`CGMTreePrior::drawSplitVariable`](../../src/bartcore/model.hpp) draws a cutoff of 0 and returns the first
  available column, and [`CGMTreePrior::splitVariableLogProbability`](../../src/bartcore/model.hpp) prices it as
  log(0 / 0).
- Defect 3, the mixture and the variance forest: fixed on main by dec-A205
  ([variance-forest-mixture-reinstall.md](variance-forest-mixture-reinstall.md), 2b0df1c0 to d02e61d6), and kept
  by dec-A239, whose move entry writes every forest. Run: a heteroscedastic sampler set to change-only moves by
  `$setControl` after 50 sweeps changes its mean and variance node counts in 0 and 0 of 100 sweeps on all three
  routes (90 and 89 under the default). The mixture has been the control's since then, so no model states one.
- Defect 4, raw leaf-prior slots: present. `m@leaf.prior@k <- 9` alone: taken, engine k 2, the stored hyperprior
  2, the stored spec 9. `m@leaf.prior@prior.sd <- 0.3` alone: taken, `k.scale` 3.06 and `prior.scale` NA as
  before, and `$getLeafPrior()` reports an sd of 0.3 that is not in force.
- Defect 5, a residual prior on a family with no residual scale: present. A probit sampler given `fixed(4)`, then
  `chisq(10, 0.5)`: both taken and stored, sigma 1 on all three routes.
- Other doors. `$setControl` routes a changed mixture through `$setModel(model)` unless
  [`setModelRefused`](../../R/bartcore.R) holds (amplitudes, counts, DART), then through the move entry.
  `$setLeafPrior` has its own entries. xbart's per-cell model ([`cellModel`](../../R/xbart.R)) calls the bridge
  directly, varying power, base and k. No consumer calls `$setModel` or `$setLeafPrior` (search of stan4bart,
  bartCause, treatSens and bairrtt).

## Design

### The line

From dec-B254. Structure names the model: the response family, the kind of leaf and its columns, the forests
with their tree counts and columns, monotone, interaction limits and blocks, the variance forest, the ordinal
category count, the hazard periods, and DART against the ordinary tree prior. Parameters are the numbers inside
it: power and base on any forest, the variance forest included; split probabilities; the leaf scale as k or sd;
the residual prior and a residual sd held fixed; the Student-t df; a negative-binomial shape held fixed; a gp
lengthscale; DART's own settings. Handing a parameter between drawn and held is a parameter change.

The rule every slice is held to: after `$setModel(B)` (or the slice 4 setter) the live sampler runs the model a
copy re-created from the stored records runs, and before any run it draws bit for bit what creating with B
draws. So every setter writes the record creation reads, on the model or, until dec-B239's migration, on the
control. One exemption, by name: a k, Student-t df, negative-binomial shape or DART concentration handed from
fixed to drawn starts at the value in force (call 8), where creation starts at the grid median or the prior's
start, so before a run it draws what creation draws from that start, not creation's own; `$setLeafPrior`
already documents the same for k ("on a sampler not yet run it is not the start creation under the new prior
would make"). A copy, and a reload after `storeState`, still draw what the live sampler draws.

Facts on the control that no model can state (the forests, the variance forest, the category count, the hazard
periods) need no refusal until the migration moves them; this item adds none for them.

### Slice 1: refuse structure, the record defects, DART's ordinary settings (R)

1. `modelStructureRefusal(sampler, newModel)` in R/bartcore.R, called first in `$setModel`, before any other
   check and before anything is stored. It compares, in order, and refuses the first difference by name with the
   sampler's value and the model's, ending "make a new sampler instead":
   - the family: the `family` slot, a model's "auto" resolved against the sampler's response as a fitting
     function resolves it; and the kind `family.spec` names (gaussian against student, among the rest);
   - the leaf kind and its columns, moved up from the bridge with a message that names the leaf model (the
     bridge's check stays as a backstop);
   - a gp leaf's `max.leaf.size` (call 3);
   - interactions (`interaction.max.order`, `interaction.forbidden`) and blocks (`block.of.column`,
     `block.tree.counts`), an absent attribute being a statement of none (dec-B248);
   - DART against the ordinary prior;
   - monotone and `forest.columns`: the R check of the columns moves here, the bridge's monotone check stays, each
     message as now.
2. Parameters the engine does not yet install are refused by name until their slice: a changed Student-t df,
   negative-binomial shape, gp lengthscale or DART setting ("cannot be changed on an existing sampler yet"). The
   branch merges once, after slice 4, so no release carries these.
3. DART. The blanket DART refusal goes. A DART sampler takes a model of the DART kind whose DART settings are its
   own (item 2 until slice 3), with power, base, k and the residual prior changed; the bridge installs them as it
   does now. `setModelRefused` is unchanged, so `$setControl` keeps its routes and dec-A239's gates hold.
4. Defect 4. A model whose raw leaf spec (`leaf.prior@k`, `leaf.prior@prior.sd`) does not resolve to its
   `leaf.hyperprior` and `prior.scale` is refused, naming `$setLeafPrior` as the way to restate the leaf prior.
   The resolution is the one `$setLeafPrior` applies, against the scale stored with the sampler, never one
   re-derived from the current response, so the sampler's own model passes at every point of its life:
   `$setModel(sampler$model)` after `$setResponse` (with and without `updateScale`) and `$setOffset`, and
   `$setControl`'s route through `$setModel(model)`, which must not newly refuse.
5. Defect 5. On a family with no residual scale (probit, logistic, ordinal, negative binomial, a hazard fit,
   multinomial) a residual prior other than the sampler's own is refused: "family 'probit' has no residual scale
   to give a prior".
6. Defect 1. Nothing to build; `$setModel`'s switch is pinned beside `$setLeafPrior`'s.
7. The `$setModel` docstring and man/dbartsSampler-class.Rd list what it takes and what it refuses.

### Slice 2: a zero split probability excludes its column (engine)

1. The rule. A column of split probability zero is never split on: it is outside the forest's allowed set, as a
   column `vars` leaves out is. A node whose remaining columns all have probability zero cannot split.
2. The engine reads it through the forest's column mask, so every site that asks which columns a node may split
   on (the prior's growth probability and draws, `growFromRoot`, the rule Gibbs and NOG moves, a state install)
   agrees without a change of its own: the mask is the `vars` set intersected with the positive-probability set,
   built at creation and rebuilt by `Chain::setModel` from the new probabilities, the block rows re-intersected.
   A forest whose probabilities are all positive keeps no mask, as now. `drawSplitVariable`'s run-out fallback
   goes: where no available column has positive probability it returns the no-column value, which no caller
   reaches through the mask, asserted in debug builds. Under DART nothing changes: its weights are drawn and
   never zero, and DART and `split.probs` are exclusive.
3. `$setModel` with a zero on a column that some chain's live trees split on is refused, naming the column, before
   any chain is written: the bridge asks the engine, which walks every chain's live trees. Kept draws are not
   checked; they are replayed for prediction, not priced.
4. What follows. A state, a copy or a warm-start donor splitting on a zero-probability column meets the mask:
   declined unforced, cut back forced, as for `vars`; a 0.9-34 state that splits on such a column is declined
   too. The reported split probabilities do not change.
5. Draws. A fit stating no zero keeps no mask and draws bit for bit as main. A fit stating a zero draws from the
   prior it states, which differs from 0.9-34's and main's where positive columns run out along a path.

### Slice 3: the parameters the model carries (engine)

1. Engine, each an install on every chain behind a [`SamplerBase`](../../src/bartcore/facade.hpp) virtual and a
   `Sampler` method that checks before it writes any chain:
   - Student-t df: a fixed value sets nu; fixed to drawn sets the estimate flag, nu starting at the value in force
     and drawn from its full conditional on the grid at the next sweep; drawn to fixed sets the stated value. The
     scale latents are redrawn given nu at the next sweep, as now.
   - Negative-binomial shape: the same three cases through `NBResponse::restoreShape`'s contract, the working
     response rebuilt.
   - gp lengthscale: per-column lengthscales installed as creation installs supplied ones, every chain's kernel
     cache flushed ([`GPGaussianLeaf`](../../src/bartcore/model.hpp)'s `invalidateKernels`) before the next
     sweep. A model that states no lengthscale states no change: the value in force stays, the one the data
     object's standardization records, and nothing is re-derived, so `$setModel(power)` after `$setPredictor` or
     `$setResponse` neither moves a call main honors nor newly refuses one (call 14). A stated change is
     refused by name on a sampler holding kept draws, which replay only under their own kernel, as the gp
     refusals of dec-B237 are.
   - DART: `a`, `b`, `rho` (NA the allowed column count), `alpha`, `update.alpha` and `update.delay` installed,
     the concentration grid rebuilt, the split weights and the count of skipped updates kept; `alpha` handed
     between drawn and held as nu is.
2. Bridge: `bartcore_setModel` reads `resid.df`, the shape `family.spec` states, the lengthscale and the DART
   slots, checks them all before the engine's `setModel`, and installs them after it, so a refusal changes
   nothing. A state keeps stating these only where drawn (dec-B196), so no state block changes.
3. R: slice 1's interim refusals go; a changed shape is also written to the control's `bartcore.shape`, which a
   copy reads; after a lengthscale call, stated or not, the record a copy reads the lengthscale from carries the
   value in force, so a model stating none after one that stated a value leaves the copy on that value.

### Slice 4: tree priors on every forest, and the shared size (engine; after forests-shared-total and open calls 1 and 2)

Planned against forests-shared-total (c0bb924d's plan): a model may hold several forests with no basis
("sharers"), each taking the fitting function's tree prior where it states none (its D7, dec-B274), and, where
no forest has a basis, the fitting function's leaf prior as the group's one size, k drawn or held (its D1 and
D8).

1. Engine: `Sampler::setTreePriors` takes a base and power per mean forest and optionally the variance forest's,
   checks the counts before writing any chain, and writes each forest's tree prior; nothing else moves. A second
   `Sampler` method installs a group's leaf prior on a K-forest chain of sharers with no forest with a basis:
   k fixed or drawn with its chi prior, or a named sd, written to every sharer through the chain's group k
   (forests-shared-total's engine items 4 and 5), a forest that stated sd untouched; fixed to drawn starts at
   the k in force (call 8). Each behind a facade virtual.
2. Bridge: an entry beside `bartcore_setMoveSettings` that resolves every forest's tree prior from the records
   and installs it, all checks first, and one for the group's leaf prior.
3. The records. `bartcore.forests` holds a forest with no basis's base and power as stated or NULL, resolved at
   the bridge to the model's tree prior where NULL, so every sharer that states none follows `$setModel`'s tree
   prior on every re-creation, as it took the fitting function's at creation. A forest with a basis that states
   none keeps its kind's default (dec-B274). `bartcore.variance` follows open call 2: under its recommendation it
   holds base and power as stated or NULL, resolved the same way, and `bartcore_setModel` installs the variance
   forest's through the same engine call, since `Chain::setModel` writes no variance tree prior.
4. The door, by open call 1. Recommended: `$setTreePrior(tree.prior = NULL, forests = NULL, variance = NULL)`,
   shaped as `$setLeafPrior`: `tree.prior` a `cgm()` or `dart()` of the sampler's kind, the model's tree prior,
   reaching every forest with no basis that states none; `forests` a list of `forest(base = , power = )` by
   label or position, a forest set there stating it from then on; `variance` a `varianceForest(base = ,
   power = )`; any other knob refused as fixed at creation; NULL leaves a part as it is.
5. `$setModel` on a sampler with amplitudes and on a multinomial one: the blanket refusals go. It takes a model
   that differs from the sampler's in the tree prior, installed on every forest with no basis that states none,
   and the residual prior. On a sampler of sharers with no forest with a basis the model's leaf prior is the
   group's (dec-B254 makes the leaf scale and its hand between drawn and held a parameter): `$setModel` installs
   it through item 1's group install, and `$setLeafPrior(leaf.prior = )` does the same. Elsewhere, a forest with a
   basis present or a multinomial sampler, a leaf prior difference is refused, naming `$setLeafPrior(forests = )`
   or `$setLeafPrior(normal(k = ))`. The bridge routes these samplers through the tree-prior, group and residual
   installs and never through `Chain::setModel`'s leaf-scale write, whose discarding of a calibrated scale was
   the reason for the refusal. `setModelRefused` keeps its meaning for `$setControl`'s route.
6. TODO variance-forest-inherited-tree-prior is closed by item 3 under either answer to open call 2: followed, or
   settable through the door.

The engine sees no record, no model object and no R name in any slice: the bridge converts.

## Fit and state

- A stored state holds chains; no block is added and `stateFormatVersion` does not move. The df, shape and DART
  concentration are stored only where drawn, as now.
- `copy()`, a reload and a dead-pointer `setState` re-create from the stored model and control, which every setter
  writes, so each runs what the live sampler runs.
- A fit or sampler saved by a development build is owed nothing (dec-B461); 0.9-34 objects carry none of the
  records touched.
- From slice 2, a state that splits on a column of split probability zero, a 0.9-34 one included, is declined by
  `setState` unforced and cut back forced, as one splitting outside `vars` is. Main installs it.

## Help and NEWS

- man/dbartsSampler-class.Rd: `setModel` says what it changes and what it refuses, each by name; slice 3 adds the
  df, shape, lengthscale and DART settings; slice 4 adds the forests and the door.
- man/dbarts.Rd and the `split.probs` help of `cgm()`: a zero excludes the column, `vars`-like, on every path.
- docs/design/single-forest-vars.md: "Left out" loses its first sentence and the zero rule gets a section, the
  design note slice 2's class needs.
- NEWS, 1.0-0: `setModel` refuses a model of another structure by name where 0.9-34 installed another family's
  model with sigma held at 1; it takes the df, shape, lengthscale, DART settings and every forest's tree prior; a
  zero split probability now excludes its column, which changes results of 0.9-34 fits that state one where the
  columns with positive probability run out along a path, and declines a stored state that splits on it.

## Tests and gates

Every slice: the twin oracle. A model installed before any run draws bit for bit what creation with it draws; a
model installed mid-run leaves the live sampler, its copy and its reload drawing the same, to the reassociation
of fits dec-A239 measured. And a bitwise A/B against an e2f37f6d library, one scratch script, two chains, every
channel compared with `identical()`, over every call main accepts and the slice does not change: `$setModel` with
power, base, split probabilities, k and sd in each spelling and the residual prior on gaussian, probit, aft,
Student-t, heteroscedastic, monotone, linear-leaf and gp samplers, before and after a run, then a copy and a
reload; `$setModel(sampler$model)` after `$setResponse` and `$setOffset`; `$setLeafPrior`; `$setControl` on each
route, the mixture change through `$setModel(model)` included, after a `$setResponse`; a k-fold `xbart()` with
two cells per tree count. The masked forests, whose code slice 2 changes, each at creation and through
`$setModel`, where the mask is rebuilt, the block rows re-intersected and each tree's mask pointer reset: a
`vars`-restricted forest, `blocks()`, the BCF moderator mask (through `$setModel` only from slice 4, main
refusing it on a sampler with amplitudes) and a variance forest with `vars`. From slice 3 on, a gp sampler's
`$setModel(power)` after `$setPredictor` and after `$setResponse`, with and without kept draws. A moved stream
is a stop.

Slice 1 (R surface): a new test-setmodel-structure.R. Each structural difference refused with its message, in
each direction where there are two, and the stored model and the next draws untouched after the refusal (twin);
the interim parameter refusals; DART taking power, base, k and the residual prior, bitwise as creation before a
run; defects 4 and 5 refused, a consistent raw spec taken; the k column after `$setModel` both ways; xbart's
cells unaffected. Pins that relied on acceptance are respelled. Gates: the files touched, the A/B, the four
snapshot files unchanged, lint, air, doc-freshness, codoc; CI the rest. Mutants: the family compare reading the
slot alone; an absent interaction attribute read as no difference; the leaf-spec check on `k` only; DART's
settings not compared.

Slice 2 (posterior-changing): tests/cpp: the mask from probabilities alone, with `vars` and with blocks, rebuilt
by `setModel` and dropped when every probability is positive; a node whose positive columns are used up has
growth probability zero; the refusal walks every chain and writes none. tinytest: the run-out pin turns (0
splits on the zero column, where `vars` gives 0); `split.probs = c(0.5, 0.5, 0)` draws bit for bit what the same
probabilities with `vars` naming the first two draw; the setModel refusal names the column, the sampler
untouched; a state declined unforced and cut back forced, and a forced warm start cut back. tests/cpp
`test_split_probs_runout` calls `CGMTreePrior::drawSplitVariable` directly on a node whose available columns all
have probability zero and expects the no-column value, so the fallback kept fails in a release build too, not
only a debug assert. Gates: tests/cpp plain and under ASan/UBSan, the full suite,
the A/B on fits stating no zero, the equivalence trio bitwise (55, 19 and 11 against equivalence-32418b6e,
bcf-equivalence-20eed7ed and multinomial-equivalence-80b1c8d4; no scenario states `split.probs`, so no
re-record), the four snapshot files unchanged, and the exact-posterior evidence: a new split-probs-exact.R in
exact-gates.yaml, a single tree over two 0/1 columns held to one cut and a third at probability zero, enumerated
with the split-variable prior over the available columns of positive probability, and a non-uniform positive
pair; and a `vars` arm, the third column left out by `vars` instead, whose enumeration sets the growth
probability to 0 explicitly at a node with no allowed column left, so the two arms agree only if the zero is
read as `vars` is; plus bcf-exact-restricted.R and mask-redraw-exact.R quick, on the mask path. The R-loaded
sanitizer run ([Gate hygiene](README.md#gate-hygiene)) over the slice's tinytest files. `R CMD check --as-cran`. No
bench compare: a fit stating no zero takes the same code. Mutants: the mask built at creation and not rebuilt by
`setModel`; the vars set not intersected; the refusal reading chain 0 alone; the fallback kept, caught by
`test_split_probs_runout`.

Slice 3 (engine, neutral): tests/cpp: each install on every chain, the count and value checks writing nothing on
a refusal, the DART grid rebuilt with weights and skipped count kept. tinytest: for each parameter, the twin
oracle, the three Student-t and shape cases each, fixed to drawn under its exemption (the start pinned at the
value in force); a lengthscale change refused with kept draws; a gp model stating no lengthscale after
`$setPredictor` and after `$setResponse`, taken with the lengthscale in force unchanged and kept draws not
refused; a stated lengthscale, then a model stating none, the copy and reload running the stated one; a
lengthscale change's next draws equal a sampler created with it on the same data; the control's shape
record after a change. Gates: tests/cpp plain and sanitized, the R-loaded sanitizer run over the slice's tinytest
files, full suite, the A/B with its gp arm, the equivalence trio bitwise,
snapshots unchanged, exact gates quick on t-exact.R, negbin-exact.R and heteroscedastic-exact.R as smoke on the
changed paths, `--preclean` on every commit (a facade virtual), check. Mutants: drawn to fixed not resetting the
flag; the shape installed without the working response rebuilt; the DART grid not rebuilt, or the weights reset;
the bridge installing before checking; the lengthscale installed with `invalidateKernels` skipped, a stale kernel
cache, caught by the next-draws pin against a fresh sampler; a model with no lengthscale re-deriving the
heuristic.

Slice 4 (engine, neutral but, under open call 2's recommendation, one call): tests/cpp: the per-forest install
and its count refusal on two chains; the group install writing every sharer and no forest that stated sd.
tinytest: each forest's and the variance forest's tree prior through the door, twin oracle and copies; on three
sharers, one stating base and power, `$setModel`'s tree prior reaching the other two and not it, and a copy
re-created equal; on sharers with no forest with a basis, `$setModel`'s leaf prior and `$setLeafPrior(leaf.prior
= )` setting the group's k fixed, drawn (from the k in force) and as an sd, every sharer equal, a forest that
stated sd unmoved; beside a forest with a basis and on a multinomial sampler, a leaf prior difference refused by
name; `$setModel` on amplitude and multinomial samplers taking the tree and residual priors; an unstated variance
forest following `$setModel`'s tree prior and a stated one not (under (b), neither follows); `$setControl`'s
routes bitwise as main; once forest-after-merge has landed (whichever lands second adds it),
`$setTreePrior(forests = )` given a bare `normal()` refused as a non-forest entry, and
`forest(leaf.prior = normal(sd = 2))` refused as `sd` is. The A/B adds bcf and three-forest samplers through
`$setModel` with an unchanged tree prior and excludes the one changed call (a heteroscedastic `$setModel` with a
new tree prior and an unstated variance forest), which the twin oracle and the exact arm below cover. Gates as
slice 3, plus forest-prior-exact.R (its sharer arms with the group's leaf prior installed by `$setModel` before
the run), bcf-exact.R and multinomial-exact.R quick, and, under open call 2's recommendation, a new arm of
heteroscedastic-exact.R's enumerated part: a sampler created under one tree prior and given another by
`$setModel` before the run, its variance-forest posterior held to the enumeration under the new one. Mutants: the
variance forest's NULL resolved at creation; the amplitude route reaching the leaf-scale write; the install
writing forest 1 alone; the tree prior written over a sharer that stated its own; the group's k written to the
first sharer only.

## Order of work

1. Slice 1 (R). 2. Slice 2 (engine, posterior-changing), independent of slice 1's code and gated alone. 3. Slice 3
(engine). Slices 1 to 3 come first and are independent of forests-shared-total and of forest-after-merge.
4. forests-shared-total lands (its own branch, merged into main). 5. Slice 4, rebased onto it, after open calls 1
and 2. Each slice lands green on the branch with its review; engine features merge into main serially, and this
branch merges once, after slice 4 and the records: a ledger entry for the calls made, this plan's Status and
Landing note, the INDEX row, TODO's entry removed, single-forest-vars-leftovers' first half and
variance-forest-inherited-tree-prior removed. The order is call 15.

With forest-after-merge (8adc5d1a), independent of this item: whichever is green first lands first. Shared file:
[`SamplerBase`](../../src/bartcore/facade.hpp), to which both add virtuals (this plan's slices 3 and 4), a
textual conflict, `--preclean` on each. Its long form folds into `sd`, so `$setTreePrior(forests = )` evaluates
in its widened vocabulary, where `normal` resolves: a bare `normal()` there reaches the door's refusal of a
non-forest entry rather than "could not find function", and `forest(leaf.prior = normal(sd = ))` is refused as
`sd` is. Its `updateBasisScale` leaves a re-derived prior that `$setModel(s$model)` must not undo, slice 4 never
taking the leaf-scale write on a sampler with amplitudes; whichever lands second pins update, `$setModel(s$model)`,
then the reader's `sd` and the next draws unchanged, and a copy equal. Its own plan names
[`forestMapEntry`](../../src/bartcore/combiner.hpp) as its conflict with forests-shared-total; this plan does
not touch it.

## Cost

- Four checks or installs where `$setModel` had one route, and in slice 4 a new method (under option a).
- A zero split probability changes 0.9-34's results on fits where the positive columns run out; a state, copy or
  donor splitting on such a column is now declined unforced.
- `$setModel` refuses models it took: another family's, a Student-t change, constraints, a raw leaf spec out of
  step, a residual prior on probit. None is in a consumer.
- Until dec-B239's migration, `$setModel` writes the control's shape record and the slice 4 door the forests'
  records, as `$setLeafPrior(forests = )` already does.
- Under open call 2's recommendation, a heteroscedastic `$setModel` with a new tree prior moves the variance
  forest's prior, which main leaves.
- The k, df, shape and concentration handoffs from fixed to drawn are not creation's start, so the twin oracle
  carries a named exemption.

## Calls made

1. A structural difference is checked in R against the stored model before anything is stored; the bridge keeps
   its checks as a backstop. Alternative: the bridge alone, which stores first and rolls back.
2. A model's family "auto" is resolved against the sampler's response as a fitting function would, then compared.
   Alternative: refuse "auto", which breaks a model built by `new()`.
3. A gp's `max.leaf.size` is structure: it decides which leaves the gp covers, and the ruling lists "Gaussian-process
   settings" beside the lengthscale. Alternative: a parameter with a setter.
4. Parameters not yet installable are refused by name in slices 1 to 3 rather than left accepted; the branch
   merges once. Alternative: leave them accepted until their slice.
5. A raw leaf spec out of step with the resolved slots is refused, naming `$setLeafPrior`. Alternative: resolve
   the model from its raw spec, which makes two slots authoritative.
6. A zero split probability on a column in use is refused, not collapsed: `$setModel` has no forced form, and
   dec-B283 has a set call reject what is plainly broken. Alternative: cut the splits back, a silent repair.
7. A zero is read through the forest's column mask, which every site already consults. Alternative: the prior's
   own functions, which `growFromRoot` and the rule moves would not see.
8. Fixed to drawn starts the k, df, shape or concentration at the value in force; drawn to fixed sets the stated
   one. These handoffs are exempt by name from the twin oracle's before-run clause, as `$setLeafPrior` already
   documents for k. Alternative: restart at the grid median, as creation does, which passes the oracle as
   written but is a jump the ruling does not ask for.
9. A lengthscale change is refused while kept draws are held, as dec-B237's gp refusals are. Alternative: drop
   the kept draws.
10. Moved to [Open calls](#open-calls) 2: whether an unstated variance forest follows `$setModel`'s tree prior.
11. `$setModel` on amplitude and multinomial samplers takes tree and residual priors, and the leaf prior where it
    is a group's (call 16); elsewhere it refuses a leaf prior, naming `$setLeafPrior`, whose writers keep the
    calibration. Alternative: keep refusing those samplers.
12. `setModelRefused` and `$setControl`'s routes do not change, so per-forest-move-settings' neutrality holds.
13. Defect 3 needs no work; defect 1 only a pin.
14. A gp model that states no lengthscale keeps the value in force; nothing is re-derived, so `$setModel(power)`
    after `$setPredictor` or `$setResponse` neither moves a call main honors nor newly refuses one, as call 8
    keeps a value in force. Alternative: the heuristic over the current data, as creation would take it, which
    moves those calls and, with kept draws, refuses them.
15. Order: slices 1 to 3 first, independent; then forests-shared-total; then slice 4, planned against it.
    Alternative: slice 4 before forests-shared-total, which would write the first forest's tree prior alone on
    samplers of sharers, against that plan's D7, and be redone.
16. On a sampler of forests with no basis and no forest with a basis, the model's leaf prior is the group's, and
    `$setModel` and `$setLeafPrior(leaf.prior = )` set it (dec-B254: the leaf scale and its hand between drawn
    and held are parameters). Alternative: refuse it there, as forests-shared-total leaves it, which gives the
    group's k no setter.
17. `$setModel`'s tree prior reaches every forest with no basis that states none, the record holding such a
    forest's base and power as stated or NULL and resolving NULL at the bridge, as that forest took the fitting
    function's at creation (forests-shared-total's D7, dec-B274). Alternative: the first forest alone, which
    gives one sharer a prior the others do not follow.

## Open calls

1. Where a caller states the tree prior of a forest past the first, or of the variance forest, on a live sampler.
   No model object carries them: they sit on the control until dec-B239's migration, which is not built.
   - (a) Recommended: `$setTreePrior(tree.prior = , forests = , variance = )`, beside `$setLeafPrior(leaf.prior =
     , forests = )` and shaped like it. Cost: one more method, and the model's tree prior reachable by two doors,
     as the leaf prior already is; with sharers both doors write every forest with no basis that states none, and
     a forest set through `forests` states it from then on and leaves that group's following. After the migration
     `$setModel` takes a model carrying every forest too.
   - (b) `$setModel(newModel, forests = , variance = )`. Cost: one door, but it takes a whole model and pieces of
     one at once; with sharers, the model's tree prior and a forest's own arrive in one call, and which wins on a
     forest named in both must be stated.
   - (c) Build the forests' part of dec-B239's migration first, so `$setModel(model)` carries every forest. Cost:
     about 30 sites moved first, a slice outside this item, after forests-shared-total; the migrated record must
     carry, per forest, whether its tree prior was stated, sharers taking the model's where not. The cleanest end
     state.
2. Whether a variance forest that states no tree prior follows a tree prior changed through `$setModel`. At
   creation it takes the model's; today `$setModel` leaves it at creation's (run: power 4 and base 0.5 leave the
   variance record at 2 and 0.95), and TODO had it as "decide".
   - (a) Recommended. It follows: `bartcore.variance` holds its base and power as stated or NULL, resolved at the
     bridge, as slice 4 does for every forest with no basis that states none (call 17) and as an unstated move
     setting follows `$setControl` (dec-B477). Cost: it changes the prior in force for a call main honors (a
     heteroscedastic `$setModel` with a new tree prior and an unstated variance forest), the posterior-changing
     class for that call: a new arm of heteroscedastic-exact.R's enumerated part, the tree prior changed by
     `$setModel` before the run and the variance forest's posterior held to the enumeration under the new one,
     the A/B excluding that call. No consumer calls `$setModel`.
   - (b) It keeps creation's value, as today, and changes only through the slice 4 door. Cost: an unstated
     variance forest follows the model's tree prior at creation and not after, unlike every other unstated
     setting; no main call moves and no exact arm is needed.

## Not verified

- The engine claims for slices 3 and 4 (the Student-t, shape, lengthscale and DART installs, the per-forest tree
  prior) are from reading the code, not prototypes.
- Whether a fit stating a zero draws bit for bit as main until its first run-out is not measured; only fits
  stating none are claimed neutral.
- Slice 4 is planned against forests-shared-total's plan, which is not built; the group install's engine size
  assumes that plan's group k as written.
- dec-B239's migration has no TODO item; this plan assumes it is not built first.
- 0.9-34 was rerun for the k column only; its behavior elsewhere is taken from the 2026-10-05 survey.
