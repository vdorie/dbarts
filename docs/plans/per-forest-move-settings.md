# per-forest-move-settings: forest() and varianceForest() take their own move settings

Status: PLANNED 2026-10-29, not built. Implements dec-B477, which confirms that call of dec-A205. One open call
([Open calls](#open-calls)); the plan builds its recommendation unless ruled otherwise. Engine tier: a blind
critique comes before the implementer.

agent: one opus implementer for both slices; one opus review over both that runs the mutants below, the same
reviewer checking the fix round unless the fixes change engine behavior. Tier "Changes draws"
([Process by risk](README.md#process-by-risk)): slice 1 edits the engine's sweep and its prior install.
rng: NEUTRAL, bit for bit, for every call that states no move setting on a forest, which is every call that
exists today. A forest that states one runs a kernel no build ran before, built from the same
posterior-preserving moves, so the posterior is unchanged ([Tests and gates](#tests-and-gates)). No recorded
scenario or snapshot can state one: expected 55, 19 and 11 identical against equivalence-32418b6e,
bcf-equivalence-20eed7ed and multinomial-equivalence-80b1c8d4, the four snapshot files unchanged, no
re-record. A moved stream is a stop.
window: after the merge to main (dec-B477), which is done. The ruling names no release; the arguments are
additive, so landing before or after 1.0-0 breaks nothing.
budget: ~800 lines: engine ~45, bridge ~110, R ~110, help and NEWS ~45, tests/cpp ~130, tinytest ~300, docs and
records ~60.

## Summary for the maintainer

Today every forest of a sampler proposes tree changes under the one control object: its move mixture, its
perturb window and its tree-shift step. After this item a forest can state its own, written on the forest the
way its tree count is: `forest(proposal.probs = , n.perturb.cuts = , treeShift = )` and
`varianceForest(proposal.probs = , n.perturb.cuts = )`. A forest that states none takes the control's, so every
existing fit draws exactly what it drew. `setControl` keeps changing the mixture and the window of every forest
that states none of its own, and leaves a forest's own settings as they were created. That removes the one
cost of today's rule, that the mean forest cannot be frozen alone on a heteroscedastic sampler: state the
variance forest's mixture when the sampler is created, then freeze the control's. The engine already keeps a mixture for each
forest; what changes is where each forest's values come from, and the tree-shift step moving from the sampler
to each forest. The shipped C header does not change.

One question for you: the variance forest never takes the tree-shift step, so a tree-shift argument on
`varianceForest()` could do nothing. I recommend leaving it off that one constructor.

## Goal

`forest()` takes `proposal.probs`, `n.perturb.cuts` and `treeShift`, and `varianceForest()` takes
`proposal.probs` and `n.perturb.cuts`; each one a forest states governs that forest, and each one it leaves
out takes the control's. `setControl` with a changed mixture or window changes every forest that states
none of its own for that setting. A live sampler, its copy and its reload run the same kernel on every forest.

## What exists today

Read, and run on main at 7c264c37 (shipped build, arm64 macOS):

- The arguments. [`forest`](../../R/model.R) takes `vars`, `basis`, `sd`, `n.trees`, `base`, `power`,
  `amplitude`, `interactions` and `blocks`; [`varianceForest`](../../R/model.R) takes `vars`, `n.trees`,
  `base` and `power`. Run: both refuse `proposal.probs`, and `forest()` `n.perturb.cuts` and `treeShift`, as
  unused arguments. The control's three are validated by [`resolveProposalProbs`](../../R/model.R),
  [`perturbCutsRefusal`](../../R/A_class.R) and [`resolveTreeShift`](../../R/dbarts.R); `treeShift` is kept
  in the slot `levelGibbs` as NA ("auto"), TRUE or FALSE.
- How the per-forest settings that exist default
  ([Where the numbers live](../design/forest-defaults-by-kind.md#where-the-numbers-live)). A forest's
  `n.trees`, `base` and `power` each take the default of its kind where not stated. The first forest's are
  written over the control's `n.trees` and the model's tree prior, where the bridge reads them
  ([`resolveSamplerSpec`](../../R/spec.R)); the other forests' ride `attr(control, "bartcore.forests")`
  ([`forestParams`](../../R/model.R)). The variance forest's ride `attr(control, "bartcore.variance")`: 40
  trees by default, and `base` and `power` the model's, resolved once at creation (run: `list(n.trees = 5,
  base = 0.95, power = 2, columns = NULL)`), so they do not follow `setModel` (TODO
  variance-forest-inherited-tree-prior). `setControl` refuses a changed `n.trees`; `setModel` changes the
  first forest's tree prior and refuses a sampler of several forests
  ([`refuseAmplitudeMutation`](../../R/bartcore.R), [`refuseCountsMutation`](../../R/bartcore.R),
  [`refuseMultiForestMutation`](../../src/R_interface_bartcore.cpp)); `$setLeafPrior(forests = )` restates a
  forest's `sd` and refuses any other knob as "fixed at creation" ([`resolveForestSpreads`](../../R/dbarts.R);
  run).
- The mixture and window at creation. The bridge reads the control's slots
  ([`parseProposalProbs`](../../src/R_interface_bartcore.cpp)) into
  [`SamplerOptions`](../../src/bartcore/chain.hpp). The single-forest chain copies them onto its one
  [`Forest`](../../src/bartcore/combiner.hpp) and [`buildVarianceForest`](../../src/bartcore/chain.hpp) onto
  the [`VarianceForest`](../../src/bartcore/chain.hpp); [`applyAmplitudeSpec`](../../src/R_interface_bartcore.cpp)
  writes the same values onto every forest's [`ForestStructureSpec`](../../src/bartcore/combiner.hpp), and
  [`buildMultinomialSampler`](../../src/R_interface_bartcore.cpp) onto the one
  [`MultinomialForestSpec`](../../src/bartcore/combiner.hpp) every category forest is built from. The engine
  already holds the six probabilities and the width on each `Forest` and on the `VarianceForest`, and reads
  them there.
- The tree shift. One [`LevelGibbsMode`](../../src/bartcore/chain.hpp) on `SamplerOptions`, which
  [`Chain::run`](../../src/bartcore/chain.hpp) reads for every mean forest, "auto" deciding per forest on that
  forest's mixture ([`structureIsFrozen`](../../src/bartcore/moves.hpp)). The level step loops over the mean
  forests alone and [`drawLevelShift`](../../src/bartcore/chain.hpp) declines a linear or gp leaf, so the
  control's `treeShift` is inert on the variance forest and on those leaves.
- `setControl`. [`dbartsSampler$setControl`](../../R/dbarts.R) refuses a changed `treeShift` (run:
  ["changing 'treeShift'"](../../inst/tinytest/test-level-fiber.R)), carries every `bartcore.*` attribute of
  the old control onto the new one, and on a changed mixture or window calls `$setModel(model)`;
  [`bartcore_setModel`](../../src/R_interface_bartcore.cpp) re-reads the control's mixture and
  [`Chain::setModel`](../../src/bartcore/chain.hpp) writes it to the first forest and the variance forest
  ([variance-forest-mixture-reinstall.md](variance-forest-mixture-reinstall.md)). Run: frozen by `setControl`
  after 50 sweeps, a heteroscedastic sampler's mean and variance structures both stand for 10 sweeps, and its
  copy carries `bartcore.variance`. A two-forest sampler and a multinomial one refuse the change and roll the
  stored control back (run, both).
- Monotone. Creation and `setControl` rewrite a defaulted control mixture to birth/death only and refuse one
  that proposes other moves ([`monotoneProposalProbs`](../../R/model.R); run: the stored control holds
  `1 0 0 0 0 0.5`).
- The doors. `bart()`, `dbarts()` and `dbartsSpec()` take `forest()` terms and `forests` lists; a formula's
  term calls `forest()` with its own arguments ([`processHit`](../../R/formulaTerms.R)), so a new formal
  reaches both spellings. Run: a single forest with a variance forest is accepted as a term and as a one-entry
  list; a multinomial model takes a one-entry `forests` list, which sets every category forest (its `n.trees`
  of 4 was used), and refuses a formula term; `xbart()` and the hurdle family refuse a `forest()` term;
  `bartBT()` has none. `dbartsSpec(proposal.probs = )` writes the control's slot (run).
- A control taken from a fit loses its `bartcore.*` attributes at `dbarts()` and `dbartsSpec()`
  (`resolveSamplerSpec`) and is refused at `bart()` and `xbart()`
  ([`refuseFitStateControl`](../../R/dbarts.R); run).
- The shipped header. [`DBARTS_C_API_HASH`](../../inst/include/dbarts/dbarts.h) covers no creation and no
  mixture install; nothing here touches the header.

## Design

### R surface

1. `forest()` gains `proposal.probs = NULL, n.perturb.cuts = NULL, treeShift = NULL` after its last formal, and
   `varianceForest()` gains `proposal.probs = NULL, n.perturb.cuts = NULL` after its last (open call 1). NULL is
   "not stated", as for every other knob. The constructors store them unvalidated, as they do the others.
2. Validation at fit time, by the control's own resolvers and with their messages: in
   [`validateForestKnobs`](../../R/model.R) for a forest, and in the variance block of `resolveSamplerSpec`
   for the variance forest. `proposal.probs` resolves to the canonical six (partial spellings filled as the
   control's are), `n.perturb.cuts` passes `perturbCutsRefusal` and is kept as a double, and `treeShift`
   resolves through `resolveTreeShift` to NA, TRUE or FALSE.
3. Monotone: a forest's own mixture passes `monotoneProposalProbs` at creation, as the control's does.
4. The record. `resolveSamplerSpec` attaches `attr(control, "bartcore.moves")`: `list(forests = , variance = )`,
   `forests` one entry per mean forest in creation order (a multinomial model's one entry covering every
   category forest), `variance` an entry or NULL; an entry is `list(proposal.probs = , n.perturb.cuts = ,
   levelGibbs = )`, each NULL where the forest does not state it (`levelGibbs` absent from the variance
   entry). It is attached only where some forest states a setting. The control's slots keep the control's
   values and never a forest's.
5. `setControl` needs no change: it carries the record with the other `bartcore.*` attributes, a changed
   mixture or window goes through `$setModel` as now, and the bridge gives it to each forest whose entry is
   NULL for that setting. A changed `treeShift` stays refused.
6. `format.dbartsVarianceForest` prints the two new knobs, `<control's>` where NULL.

### Bridge

1. One parser for a mixture and a width, from a pair of SEXPs, which `parseProposalProbs` calls for the
   control's slots and the record reader for each entry: the same checks and messages.
2. A record reader that gives, for forest f and for the variance forest, the mixture, width and level mode in
   force: the entry's where stated, else the control's.
3. Creation, one mean forest ([`createHolder`](../../src/R_interface_bartcore.cpp), after
   [`optionsFromParsed`](../../src/R_interface_bartcore.cpp)): forest 1's values onto `SamplerOptions`' mixture
   fields and `levelGibbs`, and the variance forest's onto the new `varianceMixture`, set whenever a variance
   forest exists. The monotone check after creation
   ([`proposalMixIsRefused`](../../src/R_interface_bartcore.cpp)) reads forest 1's mixture.
4. Creation, several forests (`applyAmplitudeSpec`): each forest's mixture onto its `ForestStructureSpec`, and
   its level mode onto the spec's new `levelGibbs`.
5. Creation, multinomial (`buildMultinomialSampler`): forest 1's values onto the one spec and onto
   `SamplerOptions::levelGibbs`.
6. `bartcore_setModel`: forest 1's mixture onto `ModelParameters`' mixture fields, the variance forest's onto
   `ModelParameters::varianceMixture`, and the monotone check on forest 1's. The multi-forest refusal stands.

The engine never sees the record or the control: it is handed each forest's values.

### Engine

1. A `ProposalMixture` aggregate in chain.hpp: the six probabilities and `perturbWidth`, shipped defaults.
2. `SamplerOptions::varianceMixture` and [`ModelParameters`](../../src/bartcore/chain.hpp)`::varianceMixture`,
   each `std::optional<ProposalMixture>`. `buildVarianceForest` and `Chain::setModel` write it to the variance
   forest where set, and the mean forest's values where not, which is today's engine behavior and what
   tests/cpp relies on ([`testVarianceForestMixtureInstall`](../../tests/cpp/test_sampler.cpp)).
3. `Forest::levelGibbs`, set at the three places a `Forest` is built: from `SamplerOptions::levelGibbs` in the
   single-forest constructor and in [`buildMultinomialForest`](../../src/bartcore/chain.hpp); in
   [`buildSpecifiedForest`](../../src/bartcore/chain.hpp) from the spec's new
   `std::optional<LevelGibbsMode> levelGibbs` where set, else the options'.
4. `Chain::run`'s level step reads each forest's mode: forest f takes the step where its mode is on, or
   automatic with its mixture frozen, decided before any generator call as now. `SamplerOptions::levelGibbs`
   becomes the mode of a forest built without one of its own.
5. Comments brought to the new rule: `SamplerOptions::levelGibbs`, `ModelParameters`, the
   [`sweepVarianceForest`](../../src/bartcore/chain.hpp) line that the variance forest reads the mean forest's
   mixture, and the bridge's notes in `applyAmplitudeSpec` and
   [`refuseUnsupportedAmplitudeComposition`](../../src/R_interface_bartcore.cpp) that the mixture is the fit's
   and not a forest's.

No virtual of the facade changes; `SamplerOptions` and `ModelParameters` change layout, so every engine
commit installs with `--preclean`.

## Fit and state

- A fit is its sampler, and the record rides the sampler's stored control, only where a forest stated a
  setting. A fit made without the new arguments carries exactly what it carries today.
- A stored state holds the chains, not settings: no block is added and `stateFormatVersion` does not move. A
  state installs across samplers whatever their move settings, as it does across controls today.
- `copy()`, a reload and a dead-pointer `setState` re-create the engine from the stored control, record
  included, so each runs the live sampler's kernel on every forest; the agreement dec-A205 restored is
  tested with a stated setting on each forest.
- `$setLeafPrior(forests = )` refuses the three as fixed at creation through its existing knob loop; pinned.
- A control taken from such a fit: its record is dropped at `dbarts()` and `dbartsSpec()` and the control
  refused at `bart()` and `xbart()`, as for every `bartcore.*` attribute.
- A fit or sampler saved by a development build before this item carries no record and runs every forest
  on the control's, as it did; nothing is owed it beyond that (only shipped versions are).

## Help and NEWS

- man/forest.Rd: usage, and one item for the three: each stated here governs this forest, each left out
  takes the control's, `setControl` changes a forest's setting only where the forest states none, and a
  forest's own is fixed for the sampler's life; `treeShift` is inert on linear and gp leaves, as the
  control's is.
- man/varianceForest.Rd: usage and one item for the two, the same rule; a `rule_gibbs` share is spent on no
  move there, as the control's already is
  ([Landing note](variance-forest-mixture-reinstall.md#landing-note)).
- man/dbartsControl.Rd: the `proposal.probs`, `n.perturb.cuts` and `treeShift` items say each is the value of
  every forest that states none of its own; the `setControl` sentence "on a heteroscedastic sampler's
  variance forest as on its mean forest" becomes "on every forest that states none of its own"; the
  `treeShift` item's "a mixture frozen with `$setModel`" becomes `$setControl`.
- man/dbartsSampler-class.Rd, the control item: "which `setModel` may freeze" becomes "which `setControl` may
  freeze, on every forest that states no mixture of its own".
- docs/architecture.md, [One sweep](../architecture.md#one-sweep) step 2: "When the control's `treeShift`"
  becomes "When a forest's `treeShift`, its own or else the control's".
- NEWS, one sentence: in 1.0-0's "Tree moves" item if it lands before 1.0-0 ships, else under the next
  version: `forest()` and `varianceForest()` take `proposal.probs` and `n.perturb.cuts`, and `forest()`
  `treeShift`, of their own, defaulting to the control's, and `setControl` changes the forests that state
  none.

## Tests and gates

tests/cpp, in test_sampler.cpp:

- A variance mixture handed at creation lands on the variance forest and the mean forest keeps the options';
  absent, the variance forest takes the mean's (the existing pin stays). The same through `setModel`.
- `Forest::levelGibbs` is the options' mode on a single-forest and a multinomial chain, and the spec's where
  set on a two-forest chain.
- The level-step identity, two forests, the second frozen: modes (off, on) draw bitwise what (auto, auto)
  draws, auto being off on a growing forest and on on a frozen one. (on, on) and (off, off) each differ from
  it.

tinytest, a new test-forest-move-settings.R, one chain, small n:

- Twins, bitwise in draws and trees: a forest stating the control's own values draws what the same call
  without them draws, for each setting, on one forest, on a heteroscedastic sampler (mean and variance), on
  two forests (each), and on a multinomial one.
- `forest(treeShift = "always")` under a control's "never" draws what a control's "always" draws. Two
  forests, the second frozen by its own mixture: the first `"never"` and the second `"always"` under a
  control's "never" draw what a control's "auto" draws.
- Kernels apart. A heteroscedastic sampler created with `varianceForest(proposal.probs = )` a moving mixture,
  grown, then frozen by `setControl`: the mean trees stand and the variance trees move (read off the stored
  state, as test-variance-forest-mixture.R reads them); the converse with the mean forest's stated. Two
  forests, the second frozen at creation: its trees stay stumps while the first grows.
- `setControl` reaches only forests without their own: the mixture and the window separately, a forest that
  states one and not the other taking the other.
- The copy, and `saveRDS` then `readRDS` after `storeState`, of each sampler above run what the live sampler
  runs.
- `$setModel(sampler$model)` on a sampler with stated settings draws what a twin that does not call it draws.
- Refusals by name: each setting's bad values, through `forest()` in a list and in a formula term, and
  through `varianceForest()`; the three at `$setLeafPrior(forests = )`; a monotone forest's own mixture that
  proposes other moves; `varianceForest(treeShift = )` as R's unused-argument error (open call 1).
- `dbarts()` handed such a fit's control runs every forest on the control's.
- `format()` of a `varianceForest()` shows the two knobs.

Pins updated: the structure and the formals of `forest()` in
["names(formals(forest))"](../../inst/tinytest/test-forest-predictors.R), and the printed `varianceForest()`
in test-argument-surface.R.

Mutants, each to fail a named expectation: `Chain::run` reading the options' mode instead of the forest's;
`buildSpecifiedForest` ignoring the spec's mode; the single-forest constructor leaving the forest's mode at
its default; `buildVarianceForest` and, separately, `Chain::setModel` ignoring `varianceMixture`; the bridge's
`setModel` ignoring forest 1's entry; the amplitude path reading forest 1's entry for every forest; R
resolving an unstated setting to the control's value in the record; the monotone rule skipped on a forest's
own mixture.

Gates, by the neutral class and the engine tier, on the branch's own libraries:

- `R CMD INSTALL --preclean`, plain and `--configure-args=--enable-reference-build`, `buildInfo()$mode`
  checked; tests/cpp plain and under address and undefined sanitizers.
- The full tinytest suite, one process, 0 failures.
- The equivalence trio on the reference build, `compare --bitwise --strict-coverage`: 55, 19 and 11 identical,
  0 skipped, no `max |z|` line. The main corpus' heteroscedastic scenarios cover the variance forest's
  creation path, the bcf corpus the per-forest specs, the multinomial corpus its spec.
- The four snapshot files unchanged on the reference build.
- Exact gates in quick mode that reach a changed path: heteroscedastic-exact.R, aft-exact.R,
  aft-hetero-pit.R, bcf-exact.R; and a scratch variant of heteroscedastic-exact.R part (b) that gives the mean
  and variance forests different moving mixtures and windows, which is the evidence that a forest's own
  kernel keeps the posterior. CI's exact-gates workflow runs the full list on the branch.
- `lintr::lint_package()`, `air format --check .`, `tools/check-rc-codoc.R`, `tools/check-win-drift.R`,
  `tools/check-doc-freshness.R`, each on its own exit status; `inst/NEWS.Rd` parses with its entry count;
  `R CMD check --as-cran` from a tarball of a clean copy staged outside the tree. The R-loaded sanitizer run
  is CI's: no new numerics.

## Order of work

1. Slice 1, engine and bridge: the engine items, the bridge items, the tests/cpp pins, the comments. Nothing
   writes the record yet, so the sampler draws as before: verified by tests/cpp, the full tinytest suite, the
   equivalence trio and the snapshots.
2. Slice 2, R surface: the arguments, validation, the record, the monotone rule, `format`, help, NEWS,
   architecture.md, the new tinytest file and the pins. Then the review over both, the fix round, and the
   whole gate list.
3. Records: the ledger entry for the calls made, this plan's Status and Landing note, TODO
   per-forest-move-settings removed (and the dec-A205 cost noted closed), the INDEX row. One merge at the end.

## Cost

- Two more lines of knobs on `forest()`, and a fourth `bartcore.*` record on a stated fit's control.
- A forest's own setting cannot be changed on a live sampler; the route is to leave it to the control.
- The sweep reads one byte per forest per sweep for the level step in place of one per sampler; not measured,
  and not a hot path by the README's sense (no per-node or per-observation work moves).

## Calls made

1. Names: the control's own, `proposal.probs`, `n.perturb.cuts`, `treeShift`, so the forest's knob reads as
   the control's restated. Alternative: forest-style dotted `tree.shift`, against the control's spelling.
2. Each setting defaults on its own, as `n.trees`, `base` and `power` do: a forest stating a mixture and not
   a window takes the control's window, and a changed one. Alternative: stating any setting fixes all three.
3. Every forest that states none takes the control's value, whatever its kind; no kind gets a mixture of its
   own. Alternative: kind defaults, as `n.trees` has; nothing measured asks for one.
4. A forest's own setting rides a record beside the control's slots, never written over them, unlike the first
   forest's `n.trees`. Writing the first forest's over the slot would hand it to the variance forest and
   leave `setControl` unable to tell a forest's own from the control's. Alternative: the slot, as `n.trees`.
5. One new record, `bartcore.moves`, read by one bridge reader, attached only when a forest states a setting.
   Alternative: fields on `bartcore.forests` and `bartcore.variance`, which a model of one forest does not
   carry, so it would need a third home anyway.
6. A forest's own setting is fixed for the sampler's life, as its tree count and its tree prior past the
   first forest are; `$setLeafPrior(forests = )` already refuses it. Alternative: a per-forest setter,
   which the ruling does not ask for.
7. `setControl` still refuses a changed `treeShift`: it has been fixed at creation since it shipped, and the
   ruling's clause covers what `setControl` changes. Alternative: make it changeable, for the forests that
   state none.
8. A control stating a mixture beside a forest's own is not refused: the forest's governs it, as a forest's
   `n.trees` governs the control's at `dbarts()`, and the control's still reaches the forests that state
   none. Alternative: refuse it, as `bart()` refuses its own `n.trees` beside the plain forest's.
9. Monotone: a forest's own mixture passes the creation rule; `setControl`'s rewrite and refusal of the
   control's mixture stay as built, even where the one forest states its own. Alternative: skip the
   control's check where no forest takes it.
10. A sampler of several forests and a multinomial one keep refusing a changed control mixture, as built
    (dec-B477 keeps the install as built); lifting it belongs with TODO setmodel-parameters-not-structure.
    Alternative: a mixture-only install on every sampler.
11. On a multinomial model the one `forests` entry's settings set every category forest, as its `n.trees`
    does. Alternative: refuse them there.
12. `varianceForest(proposal.probs = )` with a `rule_gibbs` share is accepted and the help says it is spent on
    no move, as the control's is today. Alternative: refuse a positive share there.
13. The engine takes optional variance mixture and per-spec mode, absent meaning today's engine behavior, and
    the bridge always hands the resolved values. Alternative: required fields, which would move the defaults
    of every tests/cpp caller.
14. Two stale help sentences that name `$setModel` as the way to freeze a mixture are corrected to
    `$setControl` as their items are edited.
15. No `bench-sampler.R` compare: see [Cost](#cost).

## Open calls

1. Does `varianceForest()` take `treeShift`? The ruling lists the tree-shift setting for both constructors,
   but the variance forest never takes the step: its leaves multiply, so it has no additive level to shift,
   and the control's `treeShift` is ignored there today.
   - (a) Recommended: `varianceForest()` takes `proposal.probs` and `n.perturb.cuts` only; a `treeShift` there
     is R's unused-argument error. Cost: the two constructors differ by one argument; should a step for the
     variance forest ever be built, the argument is added then without breaking anything.
   - (b) It takes `treeShift` and ignores it, as the control's is ignored there, so the two constructors
     match. Cost: an argument that can never do anything.
   - (c) It takes `treeShift` and refuses any value but "auto". Cost: an argument whose only legal value is its
     default.

## Not verified

- That the variance forest takes no level step is read from the code (the level loop covers the mean forests
  alone), not shown by a run: no channel reports the shift.
- stan4bart was checked by a search of its R code only: it reads none of the three settings and no
  `bartcore.*` record.
- The sweep's cost is not measured (see [Cost](#cost)).
- No run covers a sampler of several chains here; `Sampler::setModel` loops the chains, and dec-A205's
  pins cover the install on every chain.
