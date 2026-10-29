# per-forest-move-settings: forest() and varianceForest() take their own move settings

Status: PLANNED 2026-10-29, revised the same day from its blind critique (scratch/pfm-crit/critique.md, verdict
"revise first") and amended from its re-check (verdict "build" after three amendments, made;
[Revision from critique](#revision-from-critique)); not built. Implements dec-B477, which
confirms that call of dec-A205. One call is held for the maintainer ([Open calls](#open-calls)): the plan builds
without it, and the landing waits on the answer.

agent: one opus implementer for both slices; one opus review over both that runs the mutants below, the same
reviewer checking the fix round unless the fixes change engine behavior. Tier "Changes draws"
([Process by risk](README.md#process-by-risk)): slice 1 edits the engine's sweep and adds an install.
rng: NEUTRAL, bit for bit, for every call that states no move setting on a forest and calls `setControl` only
as it may today. That includes every non-default control: `treeShift` "always", "never" and "auto" on a frozen
mixture, perturb and rule_gibbs shares, any `n.perturb.cuts`, `setControl`'s install on the variance forest, and
every caller of `bartcore_setModel` (`$setModel`, `$setLeafPrior`, xbart's per-cell model install).
No recorded scenario or snapshot reaches those paths, so the gate for them is a bitwise A/B against a 7c264c37
build ([Tests and gates](#tests-and-gates)); expected 55, 19 and 11 identical against equivalence-32418b6e,
bcf-equivalence-20eed7ed and multinomial-equivalence-80b1c8d4 and the four snapshot files unchanged, no
re-record. A moved stream is a stop. New behavior, with nothing on main to compare against: a forest's own
setting, a changed `treeShift` through `setControl`, and `setControl`'s install on a sampler of several forests,
a multinomial one or a DART fit. Each runs a kernel built from the same posterior-preserving moves, and the
scratch exact variants below are the evidence the posterior holds.
window: after the merge to main (dec-B477), which is done. The ruling names no release; the arguments are
additive, so landing before or after 1.0-0 breaks nothing.
budget: ~1100 lines in the tree: engine ~60, bridge ~180, R ~155, help and NEWS ~75, tests/cpp ~150, tinytest
~460, docs and records ~35. Out of the tree, ~290 lines of scratch scripts: the A/B (~190) and three exact-gate
variants (~100).

## Summary for the maintainer

Today every forest of a sampler proposes tree changes under the one control object: its move mixture, its
perturb window and its tree-shift step. After this item a forest can state its own, written on the forest the
way its tree count is: `forest(proposal.probs = , n.perturb.cuts = , treeShift = )` and
`varianceForest(proposal.probs = , n.perturb.cuts = )`. A forest that states none takes the control's, so every
existing fit draws exactly what it drew.

`setControl` changes all three settings on every forest that states none of its own, and leaves a forest's own
as created. Two things it could not do before, it now does. It changes the tree-shift step, which until now was
fixed when the sampler was created. And it changes the settings on a model of several forests, on a
multinomial one and on a DART fit, which until now refused any change: the mixture was installed together with
the priors, and those samplers refuse a prior change, the first two because their leaf scales are calibrated
when they are created and a DART fit because its split weighting is set up then. Nothing about the mixture
itself needed the refusal. This also removes the cost of today's rule, that a mean forest
cannot be frozen alone on a heteroscedastic sampler: state the variance forest's mixture when the sampler is
created, then freeze the control's. A sampler created without one still cannot.

The shipped C header does not change. One question for you: the variance forest never takes the tree-shift step,
so a tree-shift argument on `varianceForest()` could do nothing. The ruling lists it for both constructors; I
recommend leaving it off `varianceForest()`. The rest is built without it, and adding it later in either other
form is a few lines.

## Goal

`forest()` takes `proposal.probs`, `n.perturb.cuts` and `treeShift`, and `varianceForest()` takes
`proposal.probs` and `n.perturb.cuts`; each one a forest states governs that forest, and each one it leaves
out takes the control's. `setControl` with a changed mixture, window or `treeShift` changes every forest that
states none of its own for that setting, on every kind of sampler. A live sampler, its copy and its reload run
the same kernel on every forest.

## What exists today

Read, and run on main at 7c264c37 (shipped build, arm64 macOS):

- The arguments. [`forest`](../../R/model.R) takes `vars`, `basis`, `sd`, `n.trees`, `base`, `power`,
  `amplitude`, `interactions` and `blocks`; [`varianceForest`](../../R/model.R) takes `vars`, `n.trees`,
  `base` and `power`. Run: both refuse `proposal.probs`, and `forest()` `n.perturb.cuts` and `treeShift`, as
  unused arguments. The control's mixture is filled by [`resolveProposalProbs`](../../R/model.R), which checks
  names only; its range, sum and birth checks live in the `dbartsControl` validity method in
  [A_class.R](../../R/A_class.R) (the critique's run: `c(birth_death = 0.9, change = 0.4)` resolves to a swap of
  -0.3 without error). The window is checked by [`perturbCutsRefusal`](../../R/A_class.R) and the shift by
  [`resolveTreeShift`](../../R/dbarts.R), kept in the slot `levelGibbs` as NA ("auto"), TRUE or FALSE. The
  setting is new in 1.0-0 and has never shipped (dec-B139).
- How the per-forest settings that exist default
  ([Where the numbers live](../design/forest-defaults-by-kind.md#where-the-numbers-live)). A forest's
  `n.trees`, `base` and `power` each take the default of its kind where not stated. The first forest's are
  written over the control's `n.trees` and the model's tree prior, where the bridge reads them
  ([`resolveSamplerSpec`](../../R/spec.R)); the other forests' ride `attr(control, "bartcore.forests")`
  ([`forestParams`](../../R/model.R)). The variance forest's ride `attr(control, "bartcore.variance")`: 40
  trees by default, and `base` and `power` the model's, resolved once at creation (run: `list(n.trees = 5,
  base = 0.95, power = 2, columns = NULL)`), so they do not follow `setModel` (TODO
  variance-forest-inherited-tree-prior). `setControl` refuses a changed `n.trees`; `$setLeafPrior(forests = )`
  restates a forest's `sd` and refuses any other knob as "fixed at creation"
  ([`resolveForestSpreads`](../../R/dbarts.R); run).
- The mixture and window at creation. The bridge reads the control's slots
  ([`parseProposalProbs`](../../src/R_interface_bartcore.cpp)) into
  [`SamplerOptions`](../../src/bartcore/chain.hpp). The single-forest chain copies them onto its one
  [`Forest`](../../src/bartcore/combiner.hpp) and [`buildVarianceForest`](../../src/bartcore/chain.hpp) onto
  the [`VarianceForest`](../../src/bartcore/chain.hpp); [`applyAmplitudeSpec`](../../src/R_interface_bartcore.cpp)
  writes the same values onto every forest's [`ForestStructureSpec`](../../src/bartcore/combiner.hpp), and
  [`buildMultinomialSampler`](../../src/R_interface_bartcore.cpp) onto the one
  [`MultinomialForestSpec`](../../src/bartcore/combiner.hpp) every category forest is built from. Every
  creation, a copy's, a reload's and a dead-pointer `setState`'s included, enters at
  [`bartcore_create`](../../src/R_interface_bartcore.cpp), which sends a multinomial sampler to
  [`createMultinomialDataHolder`](../../src/R_interface_bartcore.cpp), then
  [`createMultinomialCountsHolder`](../../src/R_interface_bartcore.cpp) and `buildMultinomialSampler`, and every
  other to [`createHolder`](../../src/R_interface_bartcore.cpp); xbart's folds go through
  [`bartcore_createFromHandle`](../../src/R_interface_bartcore.cpp). The engine holds the six probabilities and
  the width on each `Forest` and on the `VarianceForest`, and reads them there every sweep.
- The tree shift. One [`LevelGibbsMode`](../../src/bartcore/chain.hpp) on `SamplerOptions`, which
  [`Chain::run`](../../src/bartcore/chain.hpp) reads for every mean forest, "auto" deciding per forest and per
  sweep on that forest's mixture ([`structureIsFrozen`](../../src/bartcore/moves.hpp)). The level step loops
  over the mean forests alone and [`drawLevelShift`](../../src/bartcore/chain.hpp) declines a linear or gp leaf
  and a forest of one tree. Run: a heteroscedastic fit whose mean forest has one tree draws bitwise the same
  under "always" and "never" (40 variance trees, 100 sweeps), and with ten mean trees the two differ; the
  variance forest never takes the step. The enum is defined in chain.hpp after combiner.hpp is included, so
  `Forest` cannot name it today.
- `setControl`. [`dbartsSampler$setControl`](../../R/dbarts.R) refuses a changed `treeShift` (run:
  "changing 'treeShift'", pinned in test-level-fiber.R before this item), carries every `bartcore.*` attribute of
  the old control onto the new one, and on a changed mixture or window calls `$setModel(model)`;
  [`bartcore_setModel`](../../src/R_interface_bartcore.cpp) re-reads the control's mixture and
  [`Chain::setModel`](../../src/bartcore/chain.hpp) writes it to the first forest and the variance forest
  ([variance-forest-mixture-reinstall.md](variance-forest-mixture-reinstall.md)). Run: frozen by `setControl`
  after 50 sweeps, a heteroscedastic sampler's mean and variance structures both stand for 10 sweeps, and its
  copy carries `bartcore.variance`. `Chain::setModel` does more than install the mixture: on a model with a
  fixed sigma it puts sigma back to the fixed value, undoing a host's `setSigma`.
- Why some samplers refuse. A sampler whose forests carry amplitudes, a multinomial one and a DART fit refuse
  the change and roll the stored control back (run, all three; a DART fit with "changing a DART tree prior is
  not available"): `setControl` routes the mixture through
  [`dbartsSampler$setModel`](../../R/dbarts.R), which those samplers refuse
  ([`refuseAmplitudeMutation`](../../R/bartcore.R), [`refuseCountsMutation`](../../R/bartcore.R), the DART
  check, and behind the first two [`refuseMultiForestMutation`](../../src/R_interface_bartcore.cpp)). The
  prior install was refused on several forests because `Chain::setModel` writes the first forest's leaf scale
  from the prior alone and would discard a calibrated one (7560877d;
  [Mutation surface](../design/bcf.md#mutation-surface)), and on DART because its Dirichlet machinery is fixed
  at creation; the mixture install was refused with it (variance-forest-mixture-reinstall's calls). Nothing
  calibrated or fixed at creation depends on the mixture, the window or the shift mode, which the sweep reads
  per forest each time.
- Monotone. Creation and `setControl` rewrite a control mixture equal to the package default to birth/death
  only and refuse one that proposes other moves ([`monotoneProposalProbs`](../../R/model.R); run: the stored
  control holds `1 0 0 0 0 0.5`).
- The doors. `bart()`, `dbarts()` and `dbartsSpec()` take `forest()` terms and `forests` lists; a formula's
  term calls `forest()` with its own arguments ([`processHit`](../../R/formulaTerms.R)), so a new formal
  reaches both spellings. Run: a single forest with a variance forest is accepted as a term and as a one-entry
  list; a multinomial model takes a one-entry `forests` list, which sets every category forest (its `n.trees`
  of 4 was used), and refuses a formula term; `xbart()` and the hurdle family refuse a `forest()` term;
  `bartBT()` has none. `dbartsSpec(proposal.probs = )` writes the control's slot (run).
- A control taken from a fit loses its `bartcore.*` attributes at `dbarts()` and `dbartsSpec()`
  (`resolveSamplerSpec`) and is refused at `bart()` and `xbart()`
  ([`refuseFitStateControl`](../../R/dbarts.R); run).
- The verbose creation summary ([`printInitialSummary`](../../src/R_interface_bartcore.cpp)) prints the
  control's shift mode and mixture as the fit's.
- The shipped header. [`DBARTS_C_API_HASH`](../../inst/include/dbarts/dbarts.h) covers no creation and no
  mixture install; nothing here touches the header. stan4bart reads none of the three settings (search of its R
  code; the critique also checked its C side).

## Design

### R surface

1. `forest()` gains `proposal.probs = NULL, n.perturb.cuts = NULL, treeShift = NULL` after its last formal, and
   `varianceForest()` gains `proposal.probs = NULL, n.perturb.cuts = NULL` after its last; `treeShift` there is
   held ([Open calls](#open-calls)). NULL is "not stated", as for every other knob, and the constructors store
   the values unvalidated, as they do the others.
2. Validation at fit time. The validity method's mixture block is factored into `proposalProbsRefusal()` in
   A_class.R beside `perturbCutsRefusal`, and the validity method calls it. A forest's `proposal.probs` is
   filled by `resolveProposalProbs` and checked by `proposalProbsRefusal`, its `n.perturb.cuts` by
   `perturbCutsRefusal` and kept as a double, its `treeShift` resolved by `resolveTreeShift`. Each refusal names
   the forest by its label or position, in [`validateForestKnobs`](../../R/model.R)'s style ("forest 2's
   'proposal.probs' ..."), and the variance forest's as `varianceForest()`'s.
3. A positive `rule_gibbs` share stated on `varianceForest()` is refused: that move does nothing on the
   variance forest's leaves.
4. Monotone: a forest's stated mixture is refused if it proposes swap, change, perturb or rule_gibbs, and is
   never rewritten; the control's keeps its rule as built.
5. The record. `resolveSamplerSpec` attaches `attr(control, "bartcore.moves")`: `list(forests = , variance = )`,
   `forests` one entry per mean forest in creation order (a multinomial model's one entry covering every
   category forest), `variance` an entry or NULL. An entry is `list(proposal.probs = , n.perturb.cuts = ,
   levelGibbs = )`, each NULL where the forest states none; NA in `levelGibbs` is a stated "auto". It is
   attached only where some forest states a setting. The control's slots keep the control's values and never a
   forest's.
6. `setControl`: "levelGibbs" leaves the list of settings fixed at creation. A changed mixture, window or
   `treeShift` installs on every forest that states none, by one of two routes. A changed mixture or window, with
   or without a changed `treeShift`, goes through `$setModel(model)` on a sampler whose `$setModel` takes it, as
   today, so main's call there is unchanged. Everything else goes through the new move entry: any change on a
   sampler whose `$setModel` refuses - one carrying amplitudes or counts
   ([`samplerCarriesAmplitudes`](../../R/bartcore.R), [`samplerCarriesCounts`](../../R/bartcore.R)) or a DART
   tree prior - and a change to `treeShift` alone on any sampler, since `$setModel` would also put a fixed sigma
   back. One R predicate, true where `$setModel` refuses, picks the route. A refusal rolls the stored control
   back on either route. The record rides with the other `bartcore.*` attributes.
7. `format.dbartsVarianceForest` prints the two new knobs, `<control's>` where NULL.

### Bridge

1. One parser for a mixture and a width from a pair of SEXPs, which `parseProposalProbs` calls for the control's
   slots and the record reader for each entry: the same checks and messages, the cap clamp and `Inf` handled as
   now. The shift mode's slot read is shared the same way.
2. A record reader that gives, for every mean forest and for the variance forest, the mixture, width and shift
   mode in force: the entry's where stated (NA a stated "auto"), else the control's. Its `forests` holds one
   entry per mean forest, checked against the sampler's count, except on a multinomial sampler, where it holds
   exactly one entry, which covers every category forest, and `variance` is NULL. Any other count is an error.
   Resolving parses and checks everything, the monotone check on forest 1's resolved mixture
   ([`proposalMixIsRefused`](../../src/R_interface_bartcore.cpp)) included, and writes nothing.
3. The creation install, at two sites, which between them take every creation, copy, reload and dead-pointer
   `setState`: `createHolder`, for one forest or several, and `createMultinomialCountsHolder`, after
   `buildMultinomialSampler`, for a multinomial sampler. Where the control carries the record, the settings are
   resolved and checked, then installed through the engine's new move install before any run. Without the record
   nothing new runs, so creation is main's.
4. `bartcore_createFromHandle` refuses a control carrying the record; xbart cannot make one.
5. A new entry, `bartcore_setMoveSettings(ptr, control)`, registered beside `bartcore_setModel`: resolves and
   checks the settings from the control and its record, then installs them on every forest of every chain. A
   failed check writes nothing, so R's rollback of the stored control leaves the two agreeing. It is
   `setControl`'s route of R item 6.
6. `bartcore_setModel` resolves and checks the settings before the engine's `setModel`, and after it installs
   them through the same engine entry, which can no longer fail: a refused call changes nothing in the engine.
   So a stated setting survives a direct `$setModel`, `$setLeafPrior` and xbart's per-cell model install, and a
   changed shift mode reaches forest 1 on that route. With no record and an unchanged `treeShift`, it writes the
   values `setModel` just wrote.
7. `printInitialSummary` prints forest 1's resolved mixture and mode, and one line for each other forest, and
   the variance forest, whose values differ.
8. Comments brought to the new rule, among them `applyAmplitudeSpec`'s and
   [`refuseUnsupportedAmplitudeComposition`](../../src/R_interface_bartcore.cpp)'s notes that the mixture is the
   fit's and not a forest's.

The engine never sees the record or the control: it is handed each forest's values.

### Engine

1. `LevelGibbsMode` moves from chain.hpp to moves.hpp, beside `structureIsFrozen`, so combiner.hpp can name it.
   Beside it, a `MoveSettings` aggregate: the six probabilities, `perturbWidth` and a `LevelGibbsMode`, shipped
   defaults; the variance forest's ignores the mode.
2. `Forest::levelGibbs`, set from `SamplerOptions::levelGibbs` at the three places a `Forest` is built: the
   single-forest constructor, [`buildSpecifiedForest`](../../src/bartcore/chain.hpp) and
   [`buildMultinomialForest`](../../src/bartcore/chain.hpp).
3. `Chain::run`'s level step reads each forest's mode, at the same place in the sweep, ahead of
   `drawForestRescale` and the forest loop: forest f takes the step where its mode is on, or automatic with its
   mixture frozen, decided before any generator call.
4. The move install. `Sampler::setMoveSettings` takes one `MoveSettings` per mean forest and an optional one
   for the variance forest, checks the counts against the sampler's forests before writing any chain, then
   calls `Chain::setMoveSettings` on every chain, which writes each forest's seven move fields and mode and
   cannot fail, so no call leaves some chains written and others not.
   [`SamplerBase`](../../src/bartcore/facade.hpp) gains the virtual and `SamplerFacade` the forwarder.
   Installed before any run, it matches creating with the settings: creation reads neither the mixture nor the
   mode, and nothing does before a sweep.
5. `SamplerOptions`, [`ModelParameters`](../../src/bartcore/chain.hpp) and the forest specs do not change, nor
   does `Chain::setModel`.
6. Comments brought to the new rule: `SamplerOptions::levelGibbs` becomes the mode a forest is built with, and
   [`sweepVarianceForest`](../../src/bartcore/chain.hpp)'s line that the variance forest reads the mean
   forest's mixture goes.

A facade virtual changes, so every engine commit installs with `--preclean`.

## Fit and state

- A fit is its sampler, and the record rides the sampler's stored control, only where a forest stated a
  setting. A fit made without the new arguments carries exactly what it carries today.
- A stored state holds the chains, not settings: no block is added and `stateFormatVersion` does not move. A
  state installs across samplers whatever their move settings, as it does across controls today.
- `copy()`, a reload and a dead-pointer `setState` re-create the engine through `createHolder`, or
  `createMultinomialCountsHolder` for a multinomial sampler, from the stored control, record included, so each
  runs the live sampler's kernel on every forest, the stored control's `treeShift` after a `setControl`
  included.
- `$setLeafPrior(forests = )` refuses the three as fixed at creation through its existing knob loop; pinned.
- A control taken from such a fit: its record is dropped at `dbarts()` and `dbartsSpec()` and the control
  refused at `bart()` and `xbart()`, as for every `bartcore.*` attribute.
- A fit or sampler saved by a development build before this item carries no record and runs every forest on the
  control's, as it did; nothing is owed it beyond that (only shipped versions are).

## Help and NEWS

- man/forest.Rd: usage, and one item for the three. Each stated here governs this forest and each left out
  takes the control's; `setControl` changes a forest's setting only where it states none, and a forest's own
  is fixed for the sampler's life. A stated `proposal.probs` is a whole mixture filled against the package
  defaults as `dbartsControl()`'s is, not edits to the control's (an omitted `birth` is 0.5). A forest frozen
  by its own mixture is still grown by `n.grow.sweeps` and replaced by `warm.start`, as under the control's.
  `treeShift` acts on constant leaves of two or more trees, as the control's does. And `n =` or `p =` no longer
  abbreviate a single argument.
- man/varianceForest.Rd: usage, and one item for the two, the same rule; a positive `rule_gibbs` share is
  refused there.
- man/dbartsControl.Rd: the `proposal.probs`, `n.perturb.cuts` and `treeShift` items say each is the value of
  every forest that states none of its own; the `setControl` sentence "on a heteroscedastic sampler's variance
  forest as on its mean forest" becomes "on every forest that states none of its own, on every sampler, DART
  fits, models of several forests and multinomial ones included"; the
  `treeShift` item loses "Fixed when the sampler is created" and its "a mixture frozen with `$setModel`"
  becomes `$setControl`.
- man/dbartsSampler-class.Rd, the control item: `treeShift` leaves the list `setControl` cannot change, with
  its guard sentence, and "which `setModel` may freeze" becomes "which `setControl` may freeze, on every forest
  that states no mixture of its own".
- man/dbarts.Rd, the `forests` item: its list of what each forest carries gains the three, and its sentence
  that the control's mixture "is a property of the fit and reaches every forest" becomes every forest that
  states none of its own. man/dbartsSpec.Rd: its flat `proposal.probs` is the control's and does not reach a
  forest that states its own.
- docs/architecture.md: [One sweep](../architecture.md#one-sweep) step 2, "When the control's `treeShift`"
  becomes "When a forest's `treeShift`, its own or else the control's"; in [Tree moves](../architecture.md#tree-moves),
  the frozen mixture is read per forest, from each forest's own or the control's.
- NEWS, one sentence: in 1.0-0's "Tree moves" item if it lands before 1.0-0 ships, else under the next version:
  `forest()` and `varianceForest()` take `proposal.probs` and `n.perturb.cuts`, and `forest()` `treeShift`, of
  their own, defaulting to the control's, and `setControl` changes the forests that state none.

## Tests and gates

tests/cpp, in test_sampler.cpp:

- The move install writes each forest's fields and mode, on one forest with a variance forest, on two forests
  and on a multinomial chain, on every chain; a count that is not the sampler's is refused and leaves every
  chain as it was, on a sampler of two chains.
- Installed before any run, it draws bitwise what creation under the same settings draws, for each sampler kind.
- `Forest::levelGibbs` is the options' mode at each of the three build sites.
- The level-step identity, two forests, the second frozen: modes (off, on) draw bitwise what (auto, auto)
  draws, auto being off on a growing forest and on on a frozen one; (on, on) and (off, off) each differ.

tinytest, a new test-forest-move-settings.R, small n and several cuts per predictor:

- Cross-twins, bitwise in draws and trees: a forest stating X under a control holding Y draws what the same call
  with the control holding X draws, for every path - one forest, the heteroscedastic mean forest, the variance
  forest, each forest of two (the other stating X as well, so every forest is accounted for), the multinomial
  forests - and every setting, each run where it acts: the mixture a moving one with swap, perturb and
  rule_gibbs positive; the window at 3 and `Inf` with perturb positive; `treeShift` "always" and "never" on a
  growing forest and on a forest frozen by its own mixture.
- Stated "auto" is not unstated: `forest(treeShift = "auto")` on a forest frozen by its own mixture, under a
  control's "never", draws what a control's "auto" draws.
- Forests apart. A heteroscedastic sampler with the variance forest's moving mixture stated, grown, then frozen
  by `setControl`: the mean trees stand and the variance trees move (read off the stored state, as
  test-variance-forest-mixture.R reads them); the converse with the mean forest's stated. Two forests, the
  second frozen by its own: its trees stay stumps while the first grows.
- `setControl` reaches exactly the forests that state none, as a twin: `forest(X)` beside an unstated forest
  under control Y, then `setControl` to Y' before any run, draws what `forest(X)` beside `forest(Y')` draws;
  the mixture, the window and `treeShift` each, on one forest with a variance forest, on two forests and on a
  multinomial sampler. `setControl` to X before any run draws what creation under X draws on every kind, a DART
  fit included, `treeShift` included; on a DART fit after some sweeps it takes, where main refuses.
- A change to `treeShift` alone on a gaussian model with a fixed sigma leaves a sigma the host set with
  `setSigma`; with the mixture it is put back, as on main.
- No partial writes: a record of the wrong count handed straight to the move entry, and to `bartcore_setModel`
  with a changed tree prior, is refused, and the sampler then draws what a twin never handed it draws.
- The copy, and `saveRDS` then `readRDS` after `storeState`, of each sampler above run what the live sampler
  runs, before and after a `setControl`.
- `$setModel(sampler$model)` on a sampler with stated settings draws what a twin that does not call it draws;
  after `$setLeafPrior` with a new hyperprior, a forest stating X draws what its twin under a control holding X
  draws after the same call.
- Refusals: each setting's bad values, `c(birth_death = 0.9, change = 0.4)` and `birth = 2` among them, through
  a `forests` list and a formula term, the message naming the forest; through `varianceForest()`, a positive
  `rule_gibbs` share among them; the three at `$setLeafPrior(forests = )`; a monotone forest's stated mixture
  that proposes other moves, and one equal to the package default; `varianceForest(treeShift = )` as R's
  unused-argument error while the call is held.
- `dbarts()` handed such a fit's control runs every forest on the control's; a verbose creation prints each
  forest's values.
- `format()` of a `varianceForest()` shows the two knobs.

Pins updated: `forest()`'s structure and formals in
["names(formals(forest))"](../../inst/tinytest/test-forest-predictors.R); the printed `varianceForest()` in
test-argument-surface.R; the `treeShift` refusals in test-level-fiber.R become installs; the multinomial
refusal in test-perturb-width.R (`fixedPrior$setControl(refusedControl)` before this item) becomes an
install that takes.

Mutants, each to fail a named expectation: `Chain::run` reading the options' mode instead of the forest's; the
per-forest level call moved into the forest loop or after `drawForestRescale`; `buildSpecifiedForest` or the
single-forest constructor leaving the forest's mode at its default; the move install skipping the variance
forest, or its window; `bartcore_setModel` not installing the resolved settings, or checking them after the
engine's `setModel`; `createMultinomialCountsHolder` not installing the record; `createHolder` honoring a stated
mixture and not a stated window; the record reader taking a stated NA as unstated, or refusing a multinomial
sampler's one entry; the shared parser clamping `Inf` other than to the cap; the move entry writing forest 1
alone; `Sampler::setMoveSettings` checking each chain inside its write loop; R resolving an unstated setting to
the control's value in the record; `setControl` still refusing `treeShift`, routing a DART sampler through
`$setModel`, or routing a change to `treeShift` alone through it; the monotone rule rewriting a stated
mixture.

Gates, by the neutral class and the engine tier, on the branch's own libraries:

- `R CMD INSTALL --preclean`, plain and `--configure-args=--enable-reference-build`, `buildInfo()$mode`
  checked; tests/cpp plain and under address and undefined sanitizers.
- The A/B, on both slices: one scratch script run on the reference build against a library built from 7c264c37
  and against the branch's, every recorded channel compared with `identical()`, two chains, small n, several
  cuts per predictor. Samplers: one forest, heteroscedastic, two forests on a factor basis, multinomial, DART,
  and a linear-leaf and a monotone fit under "always". Each under the control's `treeShift` "always", "never",
  and "auto" with a control-frozen mixture; a moving mixture with swap, perturb and rule_gibbs positive under
  `n.perturb.cuts` 1, 3 and `Inf`; `setControl` freezing after some sweeps on one forest and on a
  heteroscedastic sampler (the dec-A205 install) and changing the window; and a copy and a `saveRDS`/`readRDS`
  continuation of each. Then every caller of `bartcore_setModel`, whose code changes on every call, under
  non-default settings ("always", perturb with `n.perturb.cuts` 3, and frozen): `$setModel` mid-run with a
  changed tree prior, `$setLeafPrior` with a new hyperprior, and a k-fold `xbart()` with two cells per tree
  count, so its per-cell model install runs. All identical; a difference is a stop.
- The full tinytest suite, one process, 0 failures.
- The equivalence trio on the reference build, `compare --bitwise --strict-coverage`: 55, 19 and 11 identical,
  0 skipped, no `max |z|` line. The four snapshot files unchanged on the reference build.
- Exact gates in quick mode on the changed paths: heteroscedastic-exact.R, aft-exact.R, aft-hetero-pit.R,
  bcf-exact.R and forest-prior-exact.R. Three scratch variants, the evidence that a forest's own kernel keeps
  the posterior: bcf-exact.R with its one ordinal predictor raised from 3 cells to 6, so a split has five cuts
  and windows 1 and 3 propose different ones, and the two single-tree forests under different moving mixtures
  (swap, perturb and rule_gibbs positive on one) and windows 3 and 1, its enumeration resized to match;
  forest-prior-exact.R, whose forests have several trees, with the forests under different moving mixtures and
  `treeShift` "always" on one and "never" on another; and heteroscedastic-exact.R part (b) with the mean and
  variance forests under different moving mixtures, where the window is vacuous (one binary predictor, one cut)
  and is not claimed. CI's exact-gates workflow runs the full list on the branch.
- `lintr::lint_package()`, `air format --check .`, `tools/check-rc-codoc.R`, `tools/check-win-drift.R`,
  `tools/check-doc-freshness.R`, each on its own exit status; `inst/NEWS.Rd` parses with its entry count;
  `pkgdown::check_pkgdown(".")`; `R CMD check --as-cran` from a tarball of a clean copy staged outside the
  tree. The R-loaded sanitizer run is CI's: no new numerics.

## Order of work

1. Slice 1, engine and bridge: the engine items, the bridge items (the new entry and its R wrapper included),
   the tests/cpp pins, the comments. Nothing writes the record yet and `setControl` still refuses what it
   refuses, so the sampler draws as before: verified by tests/cpp, the A/B, the full tinytest suite, the
   equivalence trio and the snapshots.
2. Slice 2, R surface: the arguments, validation and `proposalProbsRefusal()`, the record, the monotone rule,
   `setControl`'s routes and `treeShift`, `format`, help, NEWS, architecture.md, the new tinytest file and the
   pins. Then the review over both, the fix round, and the whole gate list, the A/B again on the final tip.
3. Records: the ledger entry for the calls made, dec-B139's Marked line noting it revised in part by dec-B477
   (the step is set per forest as well as on the control), this plan's Status and Landing note, the INDEX row,
   and TODO per-forest-move-settings removed with a note that the dec-A205 cost is closed for samplers created
   with a stated setting, an existing heteroscedastic sampler still unable to freeze its mean forest alone. The
   landing waits on the held call. One merge at the end.

## Cost

- Three more knobs on `forest()` and two on `varianceForest()`, a fourth `bartcore.*` record on a stated fit's
  control, and one more engine entry and facade virtual.
- A forest's own setting cannot be changed on a live sampler; the route is to leave it to the control.
- `setControl` keeps two routes, the prior install for a changed mixture or window on samplers that take it,
  and the move entry for the rest.
- The sweep reads one mode per forest per sweep in place of one per sampler; not measured, and not a hot path by
  the README's sense (no per-node or per-observation work moves).

## Calls made

1. Names: the control's own, `proposal.probs`, `n.perturb.cuts`, `treeShift`, so the forest's knob reads as the
   control's restated. Alternative: forest-style dotted `tree.shift`, against the control's spelling.
2. Each setting defaults on its own, as `n.trees`, `base` and `power` do: a forest stating a mixture and not a
   window takes the control's window, and a changed one. Alternative: stating any setting fixes all three.
3. Every forest that states none takes the control's value, whatever its kind. Alternative: kind defaults, as
   `n.trees` has; nothing measured asks for one.
4. A forest's own setting rides a record beside the control's slots, never written over them, unlike the first
   forest's `n.trees`: writing it over the slot would hand it to the variance forest and leave `setControl`
   unable to tell a forest's own from the control's. Alternative: the slot, as `n.trees`.
5. One new record, `bartcore.moves`, attached only when a forest states a setting. Alternative: fields on
   `bartcore.forests` and `bartcore.variance`, which a model of one forest does not carry.
6. Stated settings reach the engine through one move install, run after creation and before any run, so a call
   that states none creates exactly as main does. Alternative: new fields on `SamplerOptions`, the forest specs
   and `ModelParameters`, filled at every creation path, which moves every creation's code.
7. A forest's own setting is fixed for the sampler's life, as its tree count is; `$setLeafPrior(forests = )`
   already refuses it. Alternative: a per-forest setter, which the ruling does not ask for.
8. `setControl` changes `treeShift` on every forest that states none, as dec-B477's move settings include the
   step; a changed step is as safe mid-run as a changed mixture, both read each sweep. Alternative, from the
   first draft: keep refusing it.
9. `setControl` installs on a sampler of several forests, a multinomial one and a DART fit, through a move
   entry that writes no prior. The refusal there stood because the mixture rode the prior install, which those
   samplers refuse, for calibrated leaf scales and for DART's fixed Dirichlet machinery; nothing ties the
   mixture to either. Alternative: keep refusing, leaving TODO setmodel-parameters-not-structure to lift it.
10. `setControl` keeps its route through `$setModel(model)` for a changed mixture or window on samplers whose
    `$setModel` takes it, so its install there stays main's call for call, re-install of the model included;
    a change to `treeShift` alone, which main has no call to match, takes the move entry, since the re-install
    is not a no-op (it puts a fixed sigma back). Alternative: the move entry for every change on every sampler,
    which drops the re-install where main makes it.
11. A control stating a mixture beside a forest's own is not refused: the forest's governs it, as a forest's
    `n.trees` governs the control's at `dbarts()`, and the control's still reaches the forests that state none.
    Alternative: refuse it, as `bart()` refuses its own `n.trees` beside the plain forest's.
12. Monotone: a stated forest mixture is never rewritten, so one equal to the package default is refused like
    any other that proposes change; the control's keeps its rewrite, since its default was not chosen.
    Alternative: rewrite a stated default-equal mixture as the control's is, a silent change.
13. On a multinomial model the one `forests` entry's settings set every category forest, as its `n.trees` does.
    Alternative: refuse them there.
14. A positive `rule_gibbs` share stated on `varianceForest()` is refused; the control's share reaching the
    variance forest stays as built. Alternative: accept it, a stated no-op.
15. `forest(treeShift = )` is accepted on a linear or gp leaf, a no-op there as the control's is, because the
    argument acts on the constant leaves most forests have; on `varianceForest()` it could never act, which is
    the held call. Alternative: refuse it on those leaves.
16. A forest's refusals name the forest; `proposalProbsRefusal()` is shared with the control's validity method.
    Alternative: the control's messages, which name no forest.
17. xbart's handle path refuses a control carrying the record, which no door can hand it. Alternative: ignore
    the record there.
18. Two stale help sentences that name `$setModel` as the way to freeze a mixture are corrected to `$setControl`.
19. No `bench-sampler.R` compare: see [Cost](#cost).

## Open calls

1. Held for the maintainer: does `varianceForest()` take `treeShift`? The ruling lists the tree-shift setting
   for both constructors, but the variance forest never takes the step: its leaves multiply, and no shift is
   built for them (run: see [What exists today](#what-exists-today)). The plan builds (a); either
   other is a few lines added after the answer, and the landing waits on it.
   - (a) Recommended: `varianceForest()` takes `proposal.probs` and `n.perturb.cuts` only; a `treeShift` there
     is R's unused-argument error. Cost: the two constructors differ by one argument; should a step for the
     variance forest ever be built, the argument is added then without breaking anything.
   - (b) It takes `treeShift` and ignores it, as the control's is ignored there, so the two constructors match.
     Cost: an argument that can never do anything.
   - (c) It takes `treeShift` and refuses any value but "auto". Cost: an argument whose only legal value is its
     default.

## Not verified

- stan4bart was checked by a search of its R code; the critique read its C side.
- The sweep's cost is not measured (see [Cost](#cost)).
- What `$setModel(model)`'s re-install inside `setControl` changes beyond the mixture is not mapped in full;
  one effect is known (a fixed sigma put back). Call 10 keeps it only where main already makes it.

## Revision from critique

The critique (scratch/pfm-crit/critique.md, on 7b8fde3e) found the design sound and asked for a revision first.
The orchestrator's calls on it, as applied here:

- 1 and 2, blocking: the A/B against 7c264c37 and the cross-twins, both added, with the critique's mutants.
- 3: the premise of the old call 7 was wrong (the setting never shipped, dec-B139); `setControl` now installs a
  changed `treeShift` (call 8).
- 4: the refusal on several forests was traced to the prior install and nothing in the mixture; the install is
  built there (call 9), and the summary says so.
- 5 to 11: `proposalProbsRefusal()` with refusals that name the forest; the multi-cut and several-forest exact
  variants; the verbose summary; the monotone rule listed (call 12); `LevelGibbsMode` moved to moves.hpp; the
  handle path; the help pages, architecture.md, the dec-B139 mark and the wording of the TODO note.
- 12: the variance forest's `rule_gibbs` refused (call 14), the linear-leaf `treeShift` reasoned (call 15),
  partial matching noted in the help.
- The design moved from creation-time fields on `SamplerOptions`, the specs and `ModelParameters` to one move
  install after creation (call 6), which the new `setControl` routes needed anyway.
- Open call 1 stands, confirmed by the critique's run and again here, and is held for the maintainer.

The critique's re-check of 591587bb said build after three amendments, applied here with its other points:

- Every sampler whose `$setModel` refuses, DART among them (run on 7c264c37), takes the move entry; DART is in
  the A/B, the `setControl` twins and the help (R item 6, call 9).
- The multinomial install site is `createMultinomialCountsHolder`, not `createHolder`, and its record holds one
  entry for every category forest (Bridge 2 and 3).
- Every check runs before any write, in `bartcore_setModel`, the move entry and `Sampler::setMoveSettings`, so a
  refusal leaves the engine and the stored control agreeing (Bridge 2, 5 and 6, Engine 4).
- The A/B runs every caller of `bartcore_setModel` under non-default settings; the bcf-exact variant has six
  cells, so its window claim is not vacuous; a change to `treeShift` alone takes the move entry (call 10).
