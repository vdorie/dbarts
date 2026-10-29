# setmodel-parameters-not-structure: setModel changes a model's parameters and refuses its structure

Status: PLANNED 2026-10-29 on e2f37f6d (dec-B254). One call is open for the maintainer, the door for the tree
priors of forests past the first ([Open calls](#open-calls)); slice 4 waits on it, slices 1 to 3 do not.

agent: a sonnet implementer for slice 1 (R only) and one opus implementer for slices 2 to 4 (engine). Review by
tier ([Process by risk](README.md#process-by-risk)): slice 1 one sonnet review; slices 2, 3 and 4 one opus review
each that runs the mutants below, the same reviewer checking the fix round unless the fixes change engine
behavior. A blind opus critique of this plan first.
rng: by slice. Slice 1 NEUTRAL, bit for bit: refusals and R checks only. Slice 2 POSTERIOR-CHANGING for a fit
that states a zero split probability, and bitwise neutral for every other fit. Slices 3 and 4 NEUTRAL for every
call main honors; they make honored what main accepts and ignores (a Student-t df, a negative-binomial shape, a
lengthscale), which has no main draw to compare against, and add what main refuses. One call main honors moves
in slice 4: `setModel` changing the tree prior of a heteroscedastic sampler whose variance forest stated none
([Calls made](#calls-made) 10). How each class is gated: [Tests and gates](#tests-and-gates).
window: before 1.0-0, after the merge to main (dec-B254: "Build it for 1.0. It can be post main-merged."). Slice 4
writes the forests' record and is serial with forests-shared-total, whose plan is being written.
budget: ~2850 lines in the tree: slice 1 ~600 (R ~230, help and NEWS ~60, tinytest ~310), slice 2 ~550 (engine
~50, bridge ~40, tests/cpp ~90, tinytest ~150, a new exact gate ~180, design note and NEWS ~40), slice 3 ~850
(engine ~130, bridge ~160, R ~90, help ~50, tests/cpp ~150, tinytest ~270), slice 4 ~850 (engine ~60, bridge
~140, R ~220, help ~70, tests/cpp ~80, tinytest ~280); records ~50. Out of the tree, ~300 lines of scratch A/B
script.

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
carries. The fourth gives every forest's tree prior a setter, and needs one answer from you: where a caller states
the tree prior of the second forest or the variance forest on a live sampler, since today no model object
carries them. I recommend a `$setTreePrior()` beside `$setLeafPrior()`.

## Goal

`$setModel` refuses, by name and before anything is stored, a model that differs from the sampler's in any fact
dec-B254 calls structure, and installs every parameter it names on every chain, so that the live sampler, its
copy and its reload run the same model, and a model installed before any run draws what creating with it draws.
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
control.

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
   The resolution is the one `$setLeafPrior` applies.
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
   becomes unreachable and is asserted in debug builds. Under DART nothing changes: its weights are drawn and
   never zero, and DART and `split.probs` are exclusive.
3. `$setModel` with a zero on a column that some chain's live trees split on is refused, naming the column, before
   any chain is written: the bridge asks the engine, which walks every chain's live trees. Kept draws are not
   checked; they are replayed for prediction, not priced.
4. What follows. A state, a copy or a warm-start donor splitting on a zero-probability column meets the mask:
   declined unforced, cut back forced, as for `vars`. The reported split probabilities do not change.
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
   - gp lengthscale: per-column lengthscales installed as creation installs supplied ones; a model with no
     lengthscale takes the heuristic over the current data, as creation would. Refused by name on a sampler
     holding kept draws, which replay only under their own kernel, as the gp refusals of dec-B237 are.
   - DART: `a`, `b`, `rho` (NA the allowed column count), `alpha`, `update.alpha` and `update.delay` installed,
     the concentration grid rebuilt, the split weights and the count of skipped updates kept; `alpha` handed
     between drawn and held as nu is.
2. Bridge: `bartcore_setModel` reads `resid.df`, the shape `family.spec` states, the lengthscale and the DART
   slots, checks them all before the engine's `setModel`, and installs them after it, so a refusal changes
   nothing. A state keeps stating these only where drawn (dec-B196), so no state block changes.
3. R: slice 1's interim refusals go; a changed shape is also written to the control's `bartcore.shape`, which a
   copy reads.

### Slice 4: tree priors on every forest (engine; waits on open call 1)

1. Engine: `Sampler::setTreePriors` takes a base and power per mean forest and optionally the variance forest's,
   checks the counts before writing any chain, and writes each forest's tree prior; nothing else moves.
2. Bridge: an entry beside `bartcore_setMoveSettings` that resolves every forest's tree prior from the records
   and installs it, all checks first.
3. The records. `bartcore.forests` holds each forest's base and power as now; `bartcore.variance` holds them as
   stated or NULL, resolved at the bridge to the model's tree prior where NULL, so an unstated variance forest
   follows the model's tree prior on every re-creation, and `bartcore_setModel` installs it there through the
   same engine call, since `Chain::setModel` writes no variance tree prior (call 10). A forest with a basis
   that states none keeps its kind's default (dec-B274).
4. The door, by open call 1. Recommended: `$setTreePrior(tree.prior = NULL, forests = NULL, variance = NULL)`,
   shaped as `$setLeafPrior`: `tree.prior` a `cgm()` or `dart()` of the sampler's kind for the only or first
   forest, `forests` a list of `forest(base = , power = )` by label or position, `variance` a
   `varianceForest(base = , power = )`; any other knob refused as fixed at creation; NULL leaves a part as it is.
5. `$setModel` on a sampler with amplitudes and on a multinomial one: the blanket refusals go. It takes a model
   that differs from the sampler's in the first forest's tree prior and the residual prior, and refuses a leaf
   prior difference, naming `$setLeafPrior(forests = )` or `$setLeafPrior(normal(k = ))`; the bridge routes these
   samplers through the tree-prior and residual installs and never through `Chain::setModel`'s leaf-scale write,
   whose discarding of a calibrated scale was the reason for the refusal. `setModelRefused` keeps its meaning for
   `$setControl`'s route.
6. TODO variance-forest-inherited-tree-prior is closed by item 3.

The engine sees no record, no model object and no R name in any slice: the bridge converts.

## Fit and state

- A stored state holds chains; no block is added and `stateFormatVersion` does not move. The df, shape and DART
  concentration are stored only where drawn, as now.
- `copy()`, a reload and a dead-pointer `setState` re-create from the stored model and control, which every setter
  writes, so each runs what the live sampler runs.
- A fit or sampler saved by a development build is owed nothing (dec-B461); 0.9-34 objects carry none of the
  records touched.

## Help and NEWS

- man/dbartsSampler-class.Rd: `setModel` says what it changes and what it refuses, each by name; slice 3 adds the
  df, shape, lengthscale and DART settings; slice 4 adds the forests and the door.
- man/dbarts.Rd and the `split.probs` help of `cgm()`: a zero excludes the column, `vars`-like, on every path.
- docs/design/single-forest-vars.md: "Left out" loses its first sentence and the zero rule gets a section, the
  design note slice 2's class needs.
- NEWS, 1.0-0: `setModel` refuses a model of another structure by name where 0.9-34 installed another family's
  model with sigma held at 1; it takes the df, shape, lengthscale, DART settings and every forest's tree prior; a
  zero split probability now excludes its column, which changes results of 0.9-34 fits that state one where the
  columns with positive probability run out along a path.

## Tests and gates

Every slice: the twin oracle. A model installed before any run draws bit for bit what creation with it draws; a
model installed mid-run leaves the live sampler, its copy and its reload drawing the same, to the reassociation
of fits dec-A239 measured. And a bitwise A/B against an e2f37f6d library, one scratch script, two chains, every
channel compared with `identical()`, over every call main accepts and the slice does not change: `$setModel` with
power, base, split probabilities, k and sd in each spelling and the residual prior on gaussian, probit, aft,
Student-t, heteroscedastic, monotone and linear-leaf samplers, before and after a run, then a copy and a reload;
`$setLeafPrior`; `$setControl` on each route; a k-fold `xbart()` with two cells per tree count. A moved stream
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
untouched; a state and a forced warm start cut back. Gates: tests/cpp plain and under ASan/UBSan, the full suite,
the A/B on fits stating no zero, the equivalence trio bitwise (55, 19 and 11 against equivalence-32418b6e,
bcf-equivalence-20eed7ed and multinomial-equivalence-80b1c8d4; no scenario states `split.probs`, so no
re-record), the four snapshot files unchanged, and the exact-posterior evidence: a new split-probs-exact.R in
exact-gates.yaml, a single tree over two 0/1 columns held to one cut and a third at probability zero, enumerated
with the split-variable prior over the available columns of positive probability, and a non-uniform positive
pair; plus bcf-exact-restricted.R and mask-redraw-exact.R quick, on the mask path. `R CMD check --as-cran`. No
bench compare: a fit stating no zero takes the same code. Mutants: the mask built at creation and not rebuilt by
`setModel`; the vars set not intersected; the refusal reading chain 0 alone; the fallback kept.

Slice 3 (engine, neutral): tests/cpp: each install on every chain, the count and value checks writing nothing on
a refusal, the DART grid rebuilt with weights and skipped count kept. tinytest: for each parameter, the twin
oracle, the three Student-t and shape cases each; a lengthscale refused with kept draws; the control's shape
record after a change. Gates: tests/cpp plain and sanitized, full suite, the A/B, the equivalence trio bitwise,
snapshots unchanged, exact gates quick on t-exact.R, negbin-exact.R and heteroscedastic-exact.R as smoke on the
changed paths, `--preclean` on every commit (a facade virtual), check. Mutants: drawn to fixed not resetting the
flag; the shape installed without the working response rebuilt; the DART grid not rebuilt, or the weights reset;
the bridge installing before checking.

Slice 4 (engine, neutral but one call): tests/cpp: the per-forest install and its count refusal on two chains.
tinytest: each forest's and the variance forest's tree prior through the door, twin oracle and copies; an
unstated variance forest following `$setModel`'s tree prior and a stated one not; `$setModel` on amplitude and
multinomial samplers taking the tree and residual priors and refusing a leaf prior; `$setControl`'s routes
bitwise as main. The A/B excludes the one changed call (a heteroscedastic `$setModel` with a new tree prior and
an unstated variance forest) and covers it by the twin oracle. Gates as slice 3, plus forest-prior-exact.R,
bcf-exact.R and multinomial-exact.R quick. Mutants: the variance forest's NULL resolved at creation; the
amplitude route reaching the leaf-scale write; the install writing forest 1 alone.

## Order of work

1. Slice 1 (R). 2. Slice 2 (engine, posterior-changing), independent of slice 1's code and gated alone. 3. Slice 3
(engine). 4. Slice 4, after open call 1 and in turn with forests-shared-total. Each lands green on the branch with
its review; the branch merges into main once, after slice 4 and the records: a ledger entry for the calls made,
this plan's Status and Landing note, the INDEX row, TODO's entry removed, single-forest-vars-leftovers' first
half and variance-forest-inherited-tree-prior removed.

## Cost

- Four checks or installs where `$setModel` had one route, and in slice 4 a new method (under option a).
- A zero split probability changes 0.9-34's results on fits where the positive columns run out; a state, copy or
  donor splitting on such a column is now declined unforced.
- `$setModel` refuses models it took: another family's, a Student-t change, constraints, a raw leaf spec out of
  step, a residual prior on probit. None is in a consumer.
- Until dec-B239's migration, `$setModel` writes the control's shape record and the slice 4 door the forests'
  records, as `$setLeafPrior(forests = )` already does.

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
8. Fixed to drawn starts the df, shape or concentration at the value in force; drawn to fixed sets the stated one.
   Alternative: restart at the grid median, as creation does, a jump the ruling does not ask for.
9. A lengthscale change is refused while kept draws are held, as dec-B237's gp refusals are. Alternative: drop
   the kept draws.
10. An unstated variance forest follows the model's tree prior through `$setModel`, its default being the
    model's, as an unstated move setting follows `$setControl` (dec-B477). Alternative: keep creation's value, as
    today.
11. `$setModel` on amplitude and multinomial samplers takes tree and residual priors and refuses a leaf prior,
    naming `$setLeafPrior`, whose writers keep the calibration. Alternative: keep refusing those samplers.
12. `setModelRefused` and `$setControl`'s routes do not change, so per-forest-move-settings' neutrality holds.
13. Defect 3 needs no work; defect 1 only a pin.

## Open calls

1. Where a caller states the tree prior of a forest past the first, or of the variance forest, on a live sampler.
   No model object carries them: they sit on the control until dec-B239's migration, which is not built.
   - (a) Recommended: `$setTreePrior(tree.prior = , forests = , variance = )`, beside `$setLeafPrior(leaf.prior =
     , forests = )` and shaped like it. Cost: one more method, and the first forest's tree prior reachable by two
     doors, as the leaf prior already is; after the migration `$setModel` takes a model carrying every forest too.
   - (b) `$setModel(newModel, forests = , variance = )`. Cost: one door, but it takes a whole model and pieces of
     one at once.
   - (c) Build the forests' part of dec-B239's migration first, so `$setModel(model)` carries every forest. Cost:
     about 30 sites moved first, a slice outside this item, serial with forests-shared-total; the cleanest end
     state.

## Not verified

- The engine claims for slices 3 and 4 (the Student-t, shape, lengthscale and DART installs, the per-forest tree
  prior) are from reading the code, not prototypes.
- Whether a fit stating a zero draws bit for bit as main until its first run-out is not measured; only fits
  stating none are claimed neutral.
- forests-shared-total's plan is not written; slice 4's overlap with it is assumed, not read.
- dec-B239's migration has no TODO item; this plan assumes it is not built first.
- 0.9-34 was rerun for the k column only; its behavior elsewhere is taken from the 2026-10-05 survey.
