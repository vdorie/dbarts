# setmodel-parameters-not-structure: setModel changes a model's parameters and refuses its structure

Status: PLANNED 2026-10-29 on e2f37f6d (dec-B254); amended the same day after its blind critique (build after
amendments; scratch/crit-2plans/report.md). Slices 1 to 3 built and reviewed LAND on e2f37f6d (tip abdcea31).
Slice 4 amended 2026-10-30 against main 0766f069 for dec-B518 (a forest changed through `$setForest`, the
sampler's doors keeping the fitting function's level), forest-door's slice 1 (the per-forest tree-prior install,
branch forest-door-s1, reviewed LAND), tree-prior-flat's landing (dec-A246) and forests-shared-total as built on
its branch; again the same day for the rulings dec-B519, dec-B522, dec-B523 and dec-B524 and forest-door as
rewritten for them. Open calls 1 and 2 are ruled. Slice 4 waits on forests-shared-total (slice F included),
control-migration's slice 1a and forest-door's slices 1 and 2.

agent: a sonnet implementer for slice 1 (R only) and one opus implementer for slices 2 to 4 (engine). Review by
tier ([Process by risk](README.md#process-by-risk)): slice 1 one sonnet review; slices 2, 3 and 4 one opus review
each that runs the mutants below, the same reviewer checking the fix round unless the fixes change engine
behavior. The blind opus critique of this plan has run.
rng: by slice. Slice 1 NEUTRAL, bit for bit: refusals and R checks only. Slice 2 POSTERIOR-CHANGING for a fit
that states a zero split probability, and bitwise neutral for every other fit. Slices 3 and 4 NEUTRAL for every
call main honors; they make honored what main accepts and ignores (a Student-t df, a negative-binomial shape, a
lengthscale), which has no main draw to compare against, and add what main refuses. How each class is gated:
[Tests and gates](#tests-and-gates).
window: before 1.0-0, after the merge to main (dec-B254: "Build it for 1.0. It can be post main-merged."). Slice 4
lands after forests-shared-total (slice F included), control-migration's slice 1a and forest-door's slices 1 and 2,
and before forest-door's slice 3; it is planned against forests-shared-total's D1, D7 and D8 and forest-door's
per-forest stated flags.
budget: ~3390 lines in the tree: slice 1 ~630 (R ~240, help and NEWS ~60, tinytest ~330), slice 2 ~620 (engine
~50, bridge ~40, tests/cpp ~110, tinytest ~150, a new exact gate ~230, design note and NEWS ~40), slice 3 ~1150
(engine ~300, bridge ~160, R ~90, help ~50, tests/cpp ~250, tinytest ~300), slice 4 ~940 (engine ~100, bridge
~130, R ~205, help ~60, tests/cpp ~110, tinytest ~320, the variance forest's help and pin ~15); records ~50. Slice 4's tree-prior install, ~44 lines of engine, is forest-door's slice 1. Slice 3's engine was ~130 before the
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
start it. The first three are built. The fourth comes after the items that let a model hold several forests
with no basis and that give each forest its own door, `$setForest` (your dec-B518): a tree depth given to
`setModel` or to a new `$setTreePrior(power, base, split.probs)` reaches every forest with no basis that did not
state its own, and the size forests with no basis share is changed through `setModel` or `$setLeafPrior` as one
forest's is. As you ruled, where only one forest could take the value it is written whether it stated its own or not,
the later write winning; where several could and all state their own the call is refused, naming `$setForest`. A
variance forest stays as it is through all of these, as it does today, and changes through `$setVarianceForest`.

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

The rule every slice is held to: after `$setModel(B)` (or slice 4's `$setTreePrior`) the live sampler runs the model a
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
   - the family: the `family` slot, a model's "auto" stating no change (call 2); and the kind `family.spec`
     names (gaussian against student, among the rest);
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
   - Negative-binomial shape: the same three cases; a new r in force restarts omega at its cold start under that
     r, the working response rebuilt, as creation leaves it, since the first sweep's trees read it. Keeping omega
     and rebuilding the working response alone, as a state install does, would not draw before a run what
     creation draws.
   - gp lengthscale: per-column lengthscales installed as creation installs supplied ones, every chain's kernel
     cache flushed ([`GPGaussianLeaf`](../../src/bartcore/model.hpp)'s `invalidateKernels`) before the next
     sweep. A model that states no lengthscale states no change: the value in force stays, the one the data
     object's standardization records, and nothing is re-derived, so `$setModel(power)` after `$setPredictor` or
     `$setResponse` neither moves a call main honors nor newly refuses one (call 14). A stated change is
     refused by name on a sampler holding kept draws, which replay only under their own kernel, as the gp
     refusals of dec-B237 are.
   - DART: `a`, `b`, `rho` (NA the allowed column count), `alpha`, `update.alpha` and `update.delay` installed,
     the concentration grid rebuilt, the split weights and the count of skipped updates kept; `alpha` handed
     between drawn and held as nu is, and a drawn `alpha` given a new stated start restarts there, as creation
     would, while the same start keeps the value in force.
2. Bridge: `bartcore_setModel` reads `resid.df`, the shape from the control's record, which R writes from the
   shape `family.spec` states, the lengthscale and the DART slots (one reader, shared with creation), checks them
   all before the engine's `setModel`, and installs them after it, the lengthscale, which the sampler's kept
   draws can refuse, before it, so a refusal changes nothing. A state keeps stating these only where drawn (dec-B196), so no state block changes.
3. R: slice 1's interim refusals go; a changed shape is also written to the control's `bartcore.shape`, which a
   copy reads; after a lengthscale call, stated or not, the record a copy reads the lengthscale from carries the
   value in force, so a model stating none after one that stated a value leaves the copy on that value.

### Slice 4: the fitting function's level on every forest, and the shared size (engine; after forests-shared-total and forest-door)

Planned against:

- forests-shared-total, as built on its branch (slices B, 1 and 2 at 839cbb67, slice 2 in review; slice F not
  built) and its plan: forests with no basis that state no sd share the fitting function's leaf prior as one
  total, held, with one k for the group (its D1 under dec-B509, D3); each takes the fitting function's power and
  base where it states none (D7); after its slice F a lone forest with no basis beside a forest with a basis is a
  sharer too, a held k of 2 kept on dec-B290's map entry (its engine item 4). "Sharers" below are the forests its
  `forestSharers` names.
- dec-B518 and the forest-door plan as rewritten for dec-B519 to dec-B525 (branch forest-door): a forest's own sd,
  power, base and move settings change through `$setForest(forest, ...)`, the variance forest's through
  `$setVarianceForest`. Its slice 2 keeps, on the model records control-migration's slice 1a creates, each forest's
  values in force and whether it states them: the record's `sd`, NA where none is stated, and `"power"` and `"base"`
  in its `stated`; forest 1's depth stays in the model's tree-prior slot, which the bridge reads. Nothing records the
  fitting function's level apart from the forests' (dec-B522), and no flag records a variance forest's statement
  (dec-B524). Slice 4 reads and writes those records and adds none.
- The reach of a fitting-function door (dec-B523). Its candidates are the forests with no basis. One candidate, a lone
  forest, or the one beside forests with a basis, or a multinomial model's one entry: it is written whether it states
  its own or not. Several: those stating none of their own, and where all state their own the door is refused, naming
  `$setForest`. None: refused per name. A value equal to the one in force is no change, on every sampler.
- tree-prior-flat (dec-A246): plain `power`, `base` and `split.probs`; its D5 (a fitting-function value that
  reaches no forest is refused per name), D6 (`split.probs` is the mean forest's; a non-uniform one is refused on
  several forests) and its call 13 (the live door refuses per name a statement that reaches no forest, naming the
  forests that stated it and `$setForest`).

1. Engine. The tree-prior install is forest-door's slice 1: `Sampler::setForestTreePrior(forestIndex, base,
   power)` behind a facade virtual, `varianceForestIndex` naming the variance forest, and the bridge entry
   `bartcore_setForestTreePrior(ptr, forest, base, power)`, forest -1 the variance forest. Slice 4 adds no
   tree-prior method; `setTreePriors` is dropped. Built on branch forest-door-s1 (75507119 to da875712), it is on
   main with forest-door before slice 4 starts (call 23). New: one `Sampler` method behind a
   facade virtual that installs the group's leaf prior on every chain: k held or drawn with its chi prior, or a
   named sd, written to every sharer as forests-shared-total's `installSharedLeafPrior` writes it at creation, a
   named sd translated by `translateNamedSd`'s sharer branch against the single-forest chain's node scale, never
   one sharer's leaf scale. Held to drawn starts at the k in force (call 8); drawn to drawn keeps the k in force
   under the new chi prior, as `Chain::setModel` does on one forest. Every sharer's leaf scale (each its own,
   dec-B505) and every forest that is not a sharer stay as they are. It checks before writing any chain that the
   chain has sharers. After slice F, where the lone sharer beside a forest with a basis sits on dec-B290's map
   entry at a held k of 2, the install moves it to the k path for any other statement and back for a held 2, so
   the twin with creation holds both ways. Under dec-B523 the same switch takes the one candidate where it states its
   own sd, a value written whether it states its own or not: its statement is cleared, its entry becomes a sharer's of
   share 1 and it takes the leaf prior given, so `$setLeafPrior(normal())` on `forest(x1, sd = 2) + forest(x2, basis =
   z)` under probit returns k to drawn. A forest of a group of several that has left stays out (dec-B525): it is not the
   one candidate. Where the one candidate states `student(df)` or `normal()` it has no k, and the leaf prior given
   writes its coefficient's median (dec-B540): `normal(sd = v)` v, as `$setForest(1, sd = v)` does; `normal(k = v)` the
   scale creation set over v; a `chi()` k that scale over the chi prior's median k, once, with no k drawn afterward.
   The df stays; a stated sd goes. Only the spelling's reading comes from creation: nothing is re-read from the
   data, since a door never recalibrates, and nothing is divided. The scale creation set is carried on the record,
   this forest's reported `k.scale` being its median. ~40 lines beside the pins.
2. Bridge. A group entry beside `bartcore_setForestTreePrior`, reading k, its drawn flag, the chi prior and a
   named sd from the model as creation's K-forest path reads them (forests-shared-total's bridge item 1), every
   check first. `bartcore_setModel` gains a route for samplers of several forests and multinomial ones that never
   calls `Chain::setModel`, whose leaf-scale write would discard a sharer's share and a calibrated scale: the
   residual prior through the residual half of `Chain::setModel`, factored out so both routes share it; each
   changed forest's depth through `bartcore_setForestTreePrior`; the group's leaf prior through item 1's entry.
   Every check comes before the first write, so a refusal changes nothing.
3. R, `$setModel` on a sampler of several forests and on a multinomial one. The blanket refusals
   ([`refuseAmplitudeMutation`](../../R/bartcore.R), [`refuseCountsMutation`](../../R/bartcore.R)) go;
   [`setModelRefused`](../../R/bartcore.R) keeps its meaning for `$setControl`'s route. After control-migration
   the new model carries every forest's record, so `$setModel` reads a difference by level:
   - The fitting function's level, installed, by the reach above: `power` and `base` per name, a multinomial
     sampler's category forests taking one value. Nothing holds the level but the forests that run it, so a level's
     value in force is that of the first forest it reaches, and a changed value is one on all of those. Where forest 1
     takes the level, the model's tree-prior slot is the level, so a hand edit `m@tree.prior@power <- 3` changes the
     level, as main takes it.
   - The leaf prior goes where `$setLeafPrior(leaf.prior = )` sends it (dec-B515), one writer behind both
     doors: on a sampler with sharers, to the group through item 1; on a multinomial sampler, to every category
     forest's k, as `$setLeafPrior(normal(k = ))` writes it today; on the one candidate stating its own sd, through
     item 1's switch, its record's `sd` going to NA; where several candidates all state sd, or every forest has a
     basis, the model holds the fixed default creation writes there, so that default restated is equal to the value in
     force and changes nothing, as creation takes it (forests-shared-total's call 21), and anything else is refused,
     naming `$setForest(f, sd = )`. `$setLeafPrior(NULL)` stays refused, as on main. A group is read from the engine
     (`samplerSharers`), not `sharingForests`, so that forest-door's slice 3 finds a group of one, or none, here with
     nothing to change.
   - The residual prior, installed as on one forest.
   - A forest's own statement, a stated sd, power or base that differs or a change in which it states, on a sampler
     with several candidates, and the variance forest's own: refused by name, naming `$setForest(label, ...)` or
     `$setVarianceForest(...)` (dec-B518: one forest through the door that names it). A value that differs only
     because the forest follows the level is not its own. The variance forest follows nothing: a changed mean depth
     leaves it (dec-B524).
   - `split.probs`: on several forests a uniform one is the default and changes nothing, and a non-uniform one
     or `dart()` is refused by D6's message.
   - Structure: slice 1's check with control-migration's slice 1b fields, as built.
   A level value equal to the one in force writes nothing (forest-door's call 4). One that differs is written on the
   record of every forest it reaches, and forest 1's goes into the model's tree-prior slot, so a copy and a reload
   re-create what the live sampler runs. It does not make those forests state it: they go on following.
4. R, the door. `$setTreePrior(power, base, split.probs)`, the fitting function's level only, a thin writer
   over item 3's route (the stored model with the level changed, through `$setModel`'s path), so the two doors
   cannot differ. A name left out changes nothing. A `split.probs` vector changes the values; `dart(...)` on a
   DART sampler changes DART's settings through slice 3's install, a setting left out taking creation's default;
   `1` restores equal probabilities; `dart()` on an ordinary sampler, or a vector on a DART one, is structure and
   refused (tree-prior-flat's live door). No `tree.prior`, `forests` or `variance` argument. Every formal defaults to
   NULL and a written NULL is a name left out (dec-B519, dec-B522): `$setTreePrior(power = NULL)` changes nothing.
   There is no spelling that clears a level or returns a forest to it. On a one-forest sampler whose forest term stated
   its own depth, `$setTreePrior(power = )` and `$setModel` with a changed depth write it, the later write winning
   (dec-B523), as `$setModel` does on main.
5. The variance forest (dec-B524): `$setModel` and `$setTreePrior` leave it, as main does; the help says so; ~15 lines,
   a pin and a help sentence, part of items 3 and 4 and no last step. TODO variance-forest-inherited-tree-prior closes:
   it changes through `$setVarianceForest`.

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
  df, shape, lengthscale and DART settings; slice 4 adds samplers of several forests and multinomial ones, what
  reaches the forests that follow the fitting function's level and what is refused naming `$setForest`, and the
  usage and value of `$setTreePrior`, its formals NULL by default (tools/check-rc-codoc.R covers it), and that the
  variance forest stays as it is.
- man/dbarts.Rd and the `split.probs` help of `cgm()`: a zero excludes the column, `vars`-like, on every path.
- docs/design/single-forest-vars.md: "Left out" loses its first sentence and the zero rule gets a section, the
  design note slice 2's class needs.
- NEWS, 1.0-0: `setModel` refuses a model of another structure by name where 0.9-34 installed another family's
  model with sigma held at 1; it takes the df, shape, lengthscale and DART settings, and on a sampler of several
  forests the fitting function's tree and leaf priors; `$setTreePrior` is new (one item with forest-door's
  `$setForest`, whichever lands second writing it); a zero split probability now excludes its column, which
  changes results of 0.9-34 fits that state one where the columns with positive probability run out along a
  path, and declines a stored state that splits on it.

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

Slice 4 (engine, neutral). tests/cpp: the group install
on two chains writing every sharer's k, drawn flag and chi prior and no other forest, sharers of unequal counts
(38 and 37) keeping their own leaf scales, a named sd translated against the group's node scale, held to drawn
starting at the k in force, the install refused on a chain with no sharers with nothing written, and an install
before any run drawing what creation with the statement draws; after slice F, the lone sharer moved between
dec-B290's map entry and the k path both ways, and the one candidate moved from a stated sd to the k path. The tree-prior install's pins are forest-door's slice 1
(`testForestTreePriorInstall`, `testForestTreePriorMatchesCreation`, test-forest-tree-prior-install.R).

tinytest, a new test-setmodel-forests.R, every twin on two chains:

- On three sharers, one stating power: `$setModel`'s and `$setTreePrior`'s power reaching the other two and not
  it, the twin oracle against creation with the fitting function's power and that forest's own, and a copy and a
  reload `identical()` to the live sampler; the per-name refusal where several forests with no basis all state
  their own, naming them and `$setForest`; `$setTreePrior` and `$setModel` drawing bitwise alike for one level change.
- Later write wins (dec-B523): on a lone forest and on the one candidate beside a forest with a basis, stating its own
  power or sd, `$setTreePrior`, `$setModel` and `$setLeafPrior` writing it, then `$setForest` writing it back, the
  records agreeing with the engine after each; `$setLeafPrior(normal())` on a probit `forest(sd = 2)` returning k to
  drawn, from the spread in force (call 8's exemption). A forest following the level still following after a level
  change. `$setTreePrior(power = NULL)` and a round trip of `$getForest` writing nothing.
- The one candidate that draws a coefficient (dec-B540): `$setLeafPrior` with `normal(sd = )`, `normal(k = )` and a
  `chi()` k writing the median creation's reading gives on the same data, the df kept, `getK` still reporting no k,
  and draws matching `$setForest(1, sd = )` at that median; a stated sd cleared; after `$setResponse` the median still
  read from creation's scale; a copy and a reload continue.
- The group's leaf prior through `$setModel` and through `$setLeafPrior(leaf.prior = )`: k held, drawn (from the
  k in force, call 8's exemption) and as an sd, every sharer's k equal and each its own leaf scale, a forest that
  stated sd unmoved; after slice F the same beside a forest with a basis, the lone sharer's held 2 restated
  bitwise inert and a held 3 drawing as creation with it; a multinomial sampler's level and k, its one entry written
  whether it stated its own or not; the default restated taken and anything else refused where several candidates
  all state sd or every forest has a basis.
- Refused by name with the sampler untouched: a forest's own difference naming `$setForest`, the variance
  forest's naming `$setVarianceForest`, a non-uniform `split.probs` on several forests; a uniform one inert.
- Inert: `$setModel(s$model)` on every kind, stated flags included, so a following `$setForest` does what it
  did; `$setControl`'s routes bitwise as main.
- Carried from items that landed first, slice 4 landing second: coefficient-prior-opt-in's (on a `student(3)`
  forest, `$setModel` with a new depth leaves the reader's `amplitude`, the coefficient and v as they were, the
  next draws equal to a copy's) and its family pin (`$setModel` with a model of `family = student(5)` resolved as
  the family, carried at the branch's rebase); forest-after-merge's (`$setForestBasis(updateBasisScale = TRUE)`,
  then `$setModel(s$model)`, the reader's sd and the next draws unchanged, a copy equal).
- The variance forest: `$setModel` and `$setTreePrior` with a new depth leaving it, stated or not, its draws those of
  a twin stating its old depth; `$setVarianceForest` then changing it.

The A/B adds bcf and three-forest samplers through `$setModel(s$model)` and an unchanged level, and a heteroscedastic
`$setModel` with a new depth, which moves nothing main does. Gates as slice 3, plus forest-prior-exact.R (its sharer arms with the group's leaf
prior installed by `$setModel` before the run), probit-k-scale-exact.R's shared arms, bcf-exact.R and
multinomial-exact.R quick. No heteroscedastic-exact arm: no main call moves.

Mutants: the route installing forest 1's depth alone; the depth written over a forest that stated its own; the
group's k written to the first sharer only; a named sd translated against the first sharer's leaf scale, or the
install resetting every sharer to one leaf scale (each fails the unequal-count pin); the route reaching
`Chain::setModel`'s leaf-scale write (the `student(3)` and forest-after-merge pins); a forest's own difference
installed rather than refused; the stated flag written on an equal call (the inert pin); `$setTreePrior` given
installs of its own that skip the group write (the agreement pin); after slice F, the lone sharer left on the map
entry at a held 3; the one candidate stating its own refused or left unwritten (the later-write pins); several
candidates all stating their own taken in silence (the refusal pin); a level change marking the forests it reaches as
stating it (the still-following pin); a NULL read as a clearing (the NULL pin); the variance forest moved by a mean
depth change (the variance pin).

## Order of work

1. Slice 1 (R). 2. Slice 2 (engine, posterior-changing), independent of slice 1's code and gated alone. 3. Slice 3
(engine). Slices 1 to 3 are built and reviewed LAND. 4. On main before slice 4 starts: forests-shared-total with
its slice F, control-migration (slice 1a's records on the model) and forest-door (its slice 1, the tree-prior
install, and its slice 2, the stated flags and `$setForest`). 5. Slice 4, the branch rebased onto that main by its
builder. 6. Forest-door's slice 3 follows, on its own, reading the group this slice's writers read. Each slice lands green on the branch with its review; engine features merge into main serially, and this branch merges once, after slice 4 and the records: a
ledger entry for the calls made, this plan's Status and Landing note, the INDEX row, TODO's entry removed,
single-forest-vars-leftovers' first half and variance-forest-inherited-tree-prior removed. The order is call 15.

Before slice 4 is built, re-read on main, forests-shared-total's slice 2 and fix round having landed:

- the refusal of `leaf.prior` on a live sampler with sharers in [`dbartsSampler`](../../R/dbarts.R)'s
  `setLeafPrior`, which item 3 replaces with the group install; the refusal of a sharer's sd is `$setForest`'s by then,
  until forest-door's slice 3 replaces it (dec-B525);
- the engine's sharer paths (`forestSharers`, `installSharedLeafPrior`, `translateNamedSd`'s sharer branch,
  `drawSharedK`, `ForestCalibration`'s share), and every reader of the chain options' leaf prior after
  construction (a response re-anchor, a state install, a warm start), which item 1's install must update or
  leave consistent;
- the bridge's K-forest creation hand-off of k, its drawn flag, the chi prior and a named sd (its bridge item 1),
  which item 2's entry reads alike;
- [`resolveSamplerSpec`](../../R/spec.R)'s resolution of the group's leaf prior (the binary default, the default
  restated where none shares) and [`forestParams`](../../R/model.R)'s depth for sharers, so `$setModel` resolves
  what creation resolves;
- [`reportLeafPrior`](../../R/dbarts.R) and `getK` on sharers (its D8), which the pins read;
- slice F as built: whether `forestSharers` names a lone sharer beside a forest with a basis and how the held-2
  map entry is kept, which item 1's path switch follows;
- forest-door's slice 2 as built: where the per-forest stated flags live and how the model's tree-prior slot and
  forest 1's record are kept in step.

forest-after-merge and coefficient-prior-opt-in landed (dec-A244, dec-A245): slice 4 lands second to both and
carries their pins (slice 4's tinytest). Shared file: [`SamplerBase`](../../src/bartcore/facade.hpp), to which
slices 3 and 4 and forest-door's slice 1 add virtuals, a textual conflict, `--preclean` on each. `forest()` takes
no `leaf.prior` (dec-B508). forest-after-merge's `updateBasisScale` leaves a re-derived prior that
`$setModel(s$model)` must not undo, slice 4 never taking the leaf-scale write on a sampler of several forests.
[`forestMapEntry`](../../src/bartcore/combiner.hpp) is forests-shared-total's and coefficient-prior-opt-in's;
this plan does not touch it.

## Cost

- Four checks or installs where `$setModel` had one route; in slice 4 a second route for samplers of several
  forests, one more sampler method (`$setTreePrior`) and one more engine install (the group's leaf prior).
- A zero split probability changes 0.9-34's results on fits where the positive columns run out; a state, copy or
  donor splitting on such a column is now declined unforced.
- `$setModel` refuses models it took: another family's, a Student-t change, constraints, a raw leaf spec out of
  step, a residual prior on probit. None is in a consumer.
- Until control-migration lands, slice 3's `$setModel` writes the control's shape record; slice 4 writes the
  model's records control-migration and forest-door create.
- `$setModel` on a sampler of several forests refuses a forest's own difference that a model object can carry,
  naming `$setForest`: two calls where one model would have done, by dec-B518's one door per forest.
- Slice 4 waits on three other items (forests-shared-total with slice F, control-migration, forest-door).
- Two doors can overwrite each other in silence; a forest that states its own and is the only candidate is written
  by the fitting function's door (dec-B523).
- The k, df, shape and concentration handoffs from fixed to drawn are not creation's start, so the twin oracle
  carries a named exemption.

## Calls made

1. A structural difference is checked in R against the stored model before anything is stored; the bridge keeps
   its checks as a backstop. Alternative: the bridge alone, which stores first and rolls back.
2. A model's family "auto" states no change: it is read as the sampler's own family and kind, as a model with no
   `family.spec` states no shape, and a model stating no kind at all (no `family.spec`, no `resid.df`) is stored
   with the sampler's, so a copy runs what the live sampler runs. No call main honors is newly refused.
   Alternatives: resolve "auto" against the sampler's response as a fitting function would, which refuses a
   hand-built model main takes on a gaussian fit of a 0/1 response and on a negative binomial; refuse "auto",
   which breaks a model built by `new()`.
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
10. Moved to [Open calls](#open-calls) 2, ruled: an unstated variance forest does not follow `$setModel`'s tree prior.
11. `$setModel` on samplers of several forests and multinomial ones takes the fitting function's level (tree,
    leaf and residual priors) and refuses a forest's own difference by name, naming `$setForest` or
    `$setVarianceForest` (call 19). Alternative: keep refusing those samplers, which leaves the group's k and the
    sharers' depth no setter.
12. `setModelRefused` and `$setControl`'s routes do not change, so per-forest-move-settings' neutrality holds.
13. Defect 3 needs no work; defect 1 only a pin.
14. A gp model that states no lengthscale keeps the value in force; nothing is re-derived, so `$setModel(power)`
    after `$setPredictor` or `$setResponse` neither moves a call main honors nor newly refuses one, as call 8
    keeps a value in force. Alternative: the heuristic over the current data, as creation would take it, which
    moves those calls and, with kept draws, refuses them.
15. Order: slices 1 to 3 first, independent; then forests-shared-total (slice F included), control-migration and
    forest-door; then slice 4, planned against them. Alternative: slice 4 before forests-shared-total, which would
    write the first forest's tree prior alone on samplers of sharers, against that plan's D7, and be redone.
16. On every sampler with sharers the model's leaf prior is the group's, and `$setModel` and `$setLeafPrior(leaf.prior
    = )` set it (dec-B254: the leaf scale and its hand between drawn and held are parameters; dec-B515: the
    fitting function's leaf prior reaches every forest that states no sd). After forests-shared-total's slice F
    that includes a lone sharer beside a forest with a basis (dec-B509), which this plan had refused. Alternative:
    refuse it there, which gives the group's k no setter.
17. The fitting function's power and base reach every forest with no basis that states none, read from
    forest-door's per-forest stated flags, as that forest took the fitting function's at creation
    (forests-shared-total's D7, dec-B274). Alternative: the first forest alone, which gives one sharer a prior the
    others do not follow. The plan's earlier record of NULL for an unstated base and power is dropped for
    forest-door's flags. Where one forest is the only candidate, dec-B523 writes it whether it states its own or not.
18. The per-forest tree-prior install is forest-door's slice 1, one forest per call with the variance forest at
    index -1; slice 4 builds none. Alternative: the planned `setTreePriors` over every forest at once, which
    duplicates the counting the bridge does before its per-forest calls.
19. A forest's own sd, power or base differing in the model given to `$setModel`, and the variance forest's own,
    are refused naming `$setForest` and `$setVarianceForest` (dec-B518: one forest through the door that names
    it). Alternative: install each forest's own from the model's records, a second per-forest door.
20. `$setModel`'s leaf prior and `$setLeafPrior(leaf.prior = )` share one writer on every sampler, the multinomial
    one included. Alternative: refuse a leaf prior difference on a multinomial `$setModel`, naming
    `$setLeafPrior`, as this plan had it: two answers to one statement.
21. `$setTreePrior` is a thin writer over `$setModel`'s route, as `$setControl` routes a changed mixture through
    `$setModel(model)`. Alternative: installs of its own, which can drift from `$setModel`'s.
22. A level value equal to the one in force writes nothing (forest-door's call 4); one that differs is written and
    marks no forest as stating it, the level having no record (dec-B522). Alternative: a forest written by a level
    door then states it, which ends the following the door exists for.
23. Slice 4 reads forest-door's stated flags and waits for its slice 2. Alternative: record the flags itself if
    built first, two writers for one record.
24. Where a level has no forest to be read from, none reaching it, a `$setModel` difference in the slot is a forest's
    own and refused. Alternative: a record of the level, which dec-B522 rules out.

## Open calls

1. RULED, dec-B518 (2026-10-29). Where a caller states the tree prior of a forest past the first, or of the
   variance forest, on a live sampler. Shown `$setTreePrior(tree.prior, forests, variance)` (recommended),
   `$setModel(newModel, forests = , variance = )` and the forests' migration first, the maintainer asked for
   forests as first-class members and, after a design pass, took a door keyed by label: `$setForest(forest, sd = ,
   power = , base = , ...)` and `$setVarianceForest(...)`, built by the forest-door item; the sampler's
   `$setLeafPrior(leaf.prior)` and `$setTreePrior(power, base, split.probs)` take only the fitting function's
   level. Slice 4 is amended to it: its door is that `$setTreePrior`, and `$setModel` refuses a forest's own
   change naming `$setForest` (calls 18 to 21).
2. RULED, dec-B524 (2026-10-30). Whether a variance forest that states no tree prior follows a live change of the
   depth it took at creation. Shown it following (recommended first) and keeping its value, the maintainer: "Will that
   confuse users? Personally, and maybe this is just because I tend to know how these things are coded, I wouldn't
   assume that the objects were tied together just because they shared something at creation. Feel free to push back
   and tell me if I'm wrong.", then, shown the second: "Yes, use your revised recommendation." `$setTreePrior`,
   `$setModel` and `$setForest(1, )` leave it and `$setVarianceForest` changes it; the help says so.

Ruled in forest-door, followed here: dec-B519 and dec-B522 (a written NULL at the door is a name left out; no return
spelling; no record of the level) and dec-B523 (the later write wins; the reach above).

## Not verified

- The engine claims for slices 3 and 4 were from reading the code; slice 3 is built, and slice 4's per-forest
  tree-prior install is built and reviewed as forest-door's slice 1. The group install is from reading
  forests-shared-total's slice 1 as built (3b0bab56), not a prototype.
- Whether a fit stating a zero draws bit for bit as main until its first run-out is not measured; only fits
  stating none are claimed neutral.
- Slice 4 is planned against forests-shared-total's slices 1 and 2 as built on its branch, slice 2 not yet
  reviewed, and against its slice F, not built; the lone sharer's path switch assumes slice F keeps dec-B290's
  map entry as its plan says.
- Slice 4 is planned against forest-door's slice 2, not built, and control-migration's slice 1a, built on its branch:
  how the tree-prior slot and forest 1's record are kept in step is from their plans.
- The one candidate's switch from a stated sd to the k path is from reading, and so is the median a `student(df)` or
  `normal()` forest takes there (dec-B540), slice F being unbuilt.
- 0.9-34 was rerun for the k column only; its behavior elsewhere is taken from the 2026-10-05 survey.
