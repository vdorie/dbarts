# forest-after-merge: re-derive a default forest's sd from a swapped column, and the long form of a forest's sd

Status: PLANNED 2026-10-29 on 917f02c1; two calls held for the maintainer ([Open calls](#open-calls)), slice 2
waits on the first, slice 1 on neither. Rulings: dec-B276 (the window), dec-B265 and dec-B303 (updateBasisScale
and its rows), dec-B266 (the long form), with dec-B346, dec-B362, dec-B364, dec-B465, dec-B479 and dec-B488.

agent: slice 1 a sonnet implementer and one sonnet review (R surface, no draw change). Slice 2 an opus
implementer and one opus review that runs the mutants below, the same reviewer checking the fix round unless
the fixes change engine behavior ([Process by risk](README.md#process-by-risk): the engine and a prior). A
blind critique of this plan first, as that tier asks; it may be short, slice 2 being neutral for every call main
accepts and its new path held to a creation twin.
rng: slice 1 NEUTRAL, bit for bit: the long form folds into the `sd` a forest already states. Slice 2 NEUTRAL,
bit for bit, for every call main accepts (`updateBasisScale` exists on no main call, and FALSE takes today's
path); a call with `updateBasisScale = TRUE` restates a prior the caller asked to restate, held to a forest
created on the new column ([Tests and gates](#tests-and-gates)). No re-record.
window: after the merge to main, before 1.0-0 (dec-B276: "For deferring, just until it lands on main"). Additive:
both are R's unused-argument error today.
budget: ~760 lines in the tree. Slice 1 ~220 (R ~45, help ~30, NEWS ~5, tinytest ~140). Slice 2 ~540 (engine
~80, facade ~10, bridge ~45, R ~70, help ~30, tests/cpp ~130, tinytest ~150, design note and NEWS ~25).
Records ~40.

## Summary for the maintainer

Two small pieces of the forest surface you deferred until after the merge.

First, `setForestBasis(updateBasisScale = TRUE)`. A forest multiplied by a number takes its default size per
standard deviation of that column, read when the sampler is created and then held. Today a swap keeps that
number, so a caller whose first column was a placeholder samples under a prior set by the placeholder unless
they state an sd. With the argument, the swap also reads the new column's spread, over the rows in the
likelihood at that call, and restates that forest's prior from the next sweep, as `updateScale` does for the
response. The trees are left as they are. It is for burn-in. A forest whose sd was stated, or one multiplied
by a factor, has no column-derived default, and the argument is refused there by name. A copy or a reload
afterwards runs the new prior.

Second, `forest(x, basis = ~ dose, leaf.prior = normal(sd = 30))` is accepted as the long spelling of
`forest(x, basis = ~ dose, sd = 30)`. The two are one statement, and giving both is refused. A k, a prior on the
sd, or another leaf model is refused by name. Those belong on the fitting function.

Two questions for you, below. First, what happens when a swap with the argument meets a column that is
constant over the rows in. You ruled a refusal. Every later rule for a scale re-derived in burn-in keeps the
number in force with one warning, and I recommend that here too. Second, whether `normal()` with no sd should
put a stated forest back to its default. I recommend leaving that out for now.

Neither piece moves any draw of a call that works today.

## Goal

`$setForestBasis(forest, basis, updateState = NULL, updateBasisScale = FALSE)` re-derives, when TRUE, a numeric
forest's default sd from the column it installs, over the rows in the likelihood at the call, records the new
scale where every re-creation reads it, and leaves the trees and coefficients as they are. `forest()` takes
`leaf.prior = normal(sd = v)` wherever it is written, as the same statement as `sd = v`.

## What exists today

Run on main at 917f02c1, the shipped build, arm64 macOS, private library scratch/libs/fam-plan; the probes are
scratch/fam-plan/p1.R to p4.R, untracked. Data: 120 rows, a gaussian `y`, `dose ~ N(5, 2)`, one chain unless
stated.

- No `updateBasisScale`. `s$setForestBasis(2, 3 * dose, updateBasisScale = TRUE)`: "unused argument
  (updateBasisScale = TRUE)" (p1). [`dbartsSampler$setForestBasis`](../../R/dbarts.R) takes `forest`, `basis`
  and `updateState`.
- A swap keeps the scale. On `y ~ forest(x1) + forest(x2, basis = ~ dose)` the reader gives forest 2 `sd`
  0.4616 per unit of dose, `sd.stated` FALSE, `basis.scale` 2.0678 = `sd(dose)`; after
  `setForestBasis(2, 3 * dose)` the same 0.4616 and 2.0678, `data@basis.scale` unchanged, and a `copy()` the
  same (p1). The engine's swap ([`Chain::setForestBasis`](../../src/bartcore/chain.hpp)) installs values and
  nothing else (dec-B452).
- The rows the scale is read over. Created with rows 1 to 30 at weight zero: `basis.scale` 2.1170 =
  `sd(dose[31:120])`, against 2.0678 over every row; `setWeights(rep(1, 120))` leaves 2.1170 (p2). The engine
  takes them in [`Chain::rowsForBasisScale`](../../src/bartcore/chain.hpp) at creation only.
- The record. [`recordBasisScales`](../../R/dbarts.R) writes `data@basis.scale` for a forest whose entry is
  empty and never overwrites one; [`applyForestBases`](../../src/R_interface_bartcore.cpp) hands it to every
  later creation. The bridge's reader gives the scale, the map's sd and the stated sd per chain: two rows on a
  two-chain sampler, 2.0678, 0.4669, NaN on each (p3).
- Restating a forest. `$setLeafPrior(forests = list(forest(), forest(sd = 0.4)))` makes forest 2 stated at 0.4
  and keeps `basis.scale` 2.0678; `forest()` written back afterwards leaves it stated at 0.4 (p1). No route puts
  a stated forest back to its default. [`Chain::setForestMapSd`](../../src/bartcore/chain.hpp) is the write.
- Constant columns. Created on a column of ones with no sd: refused, "forest 2 ('one') states no sd and its basis
  is constant (1 on every row in the likelihood) ..." ([`basisScaleMessage`](../../src/R_interface_bartcore.cpp));
  with `sd = 1`: created. A swap to `rep(3, n)` and to `rep(0, n)`: both taken (p1; dec-B346 is built).
- A re-derived response scale on several forests: `s$setResponse(2 * y, updateScale = TRUE)` is refused,
  "setResponse(updateScale = TRUE) does not support a sampler that carries forest amplitudes" (p1; dec-B362).
  `setActiveRows` is taken on such a sampler and `setData` refused (p2).
- No long form. `leaf.prior = normal(sd = 2)` inside a formula's `forest()` term, inside `forests = list()`,
  and in `$setLeafPrior(forests = list(forest(), forest(leaf.prior = normal(sd = 2))))`: "unused argument
  (leaf.prior = normal(sd = 2))" (p2). `normal` does not resolve inside a forest argument
  ([`FOREST_ARGUMENT_VOCABULARIES`](../../R/model.R), [`processHit`](../../R/formulaTerms.R)): a bare
  `normal(sd = 2)` in `$setLeafPrior(forests = )` is "could not find function "normal"" (p1).
- What a forest's `sd` is. [`forest`](../../R/model.R) keeps `sd` as given and
  [`validateForestKnobs`](../../R/model.R) checks it at fit time; the record carries it to the engine as written
  ([`forestStatements`](../../R/model.R)). On a model of one forest `forest(x1, sd = 2)` builds the model
  `leaf.prior = normal(sd = 2)` builds, `identical()` (p3), and both stated is refused: "this model has one
  forest and states its size twice" (p2). At the fitting function on several forests, `normal(sd = 2)` is
  refused with a pointer to `forest(sd = )` and `normal(k = 2)` taken as nothing (p1).
- What the readers return. A single-forest sampler's `leaf.prior` is a `normal()` object: `normal(k = 2)` at the
  default, `normal(sd = 2)` when stated (p2, p3). A forest of several reports `forest()` or `forest(sd = v)`
  ([`reportLeafPrior`](../../R/dbarts.R)).
- The precedents. `updateScale` on `setResponse`, `setOffset`, `setWeights` and `setActiveRows` and
  `updateStandardization` re-derive over the rows in at the call and keep the number in force, with one warning
  per sampler ([`warnScaleKept`](../../R/bartcore.R); the bridge's `scaleKeptWarned` flag), where those rows
  hold fewer than two distinct response values (dec-B364) or one value of a leaf covariate (dec-B465);
  [`Sampler::rederiveStandardization`](../../src/bartcore/sampler.hpp) reads the rows through
  [`Sampler::leafRowsFrom`](../../src/bartcore/sampler.hpp).
- 0.9-34 (scratch/libs/cran, p4) has no `forest()`, no multiplier and no `$setForestBasis`; it has
  `setOffset(offset, updateScale = FALSE, updateState = NA)`, the idiom this follows. Nothing to mirror beyond it.
- No consumer calls `$setForestBasis` or writes `leaf.prior` inside a `forest()` (search of bartCause, stan4bart,
  treatSens and bairrtt R sources).

## Design

### updateBasisScale (slice 2)

The rule, from dec-B265 and dec-B303 with the later keep rulings:

1. `updateBasisScale = TRUE` on a numeric forest that states no sd: the swap installs the column, then the
   forest's scale is taken from that column over the rows in the likelihood at the call (positive weight, kept
   by the mask), its default sd is restated from it, and the new scale is recorded on the data object. Trees,
   coefficients, the unit, the anchor and every other forest are untouched; the new prior acts from the next
   sweep.
2. Too few rows to read: fewer than two rows in, or under gaussian rows in that hold fewer than two distinct
   values of the response. The column is installed and the scale in force is kept, with one warning per sampler
   (dec-B364, dec-B465; dec-B479 has a re-derivation keep where creation falls back).
3. A column constant, or constant to rounding, over those rows: OPEN CALL 1. The plan is written for its
   recommendation, the same keep and warning as item 2.
4. Refused by name, before anything is installed: a forest that states an sd ("forest 2 ('dose') states its sd,
   0.4, so it has no default to update; restate it with $setLeafPrior(forests = ) or leave updateBasisScale
   FALSE"); a forest multiplied by a factor ("... its default is one unit and reads no column"); a forest with no
   basis is refused today by the swap itself. Any value but a single TRUE or FALSE is refused.
5. FALSE, the default, is today's swap, bit for bit.

R. [`dbartsSampler$setForestBasis`](../../R/dbarts.R) gains `updateBasisScale = FALSE` after `updateState`
(Calls made 1). Before the bridge, it reads the forest's statement and kind from the reader's raw row
(`named.sd`, `multiplier`) and refuses item 4 by name. After a TRUE call it overwrites that forest's entry of
`data@basis.scale` with the engine's scale (a one-forest variant of [`recordBasisScales`](../../R/dbarts.R)
that writes even where an entry exists), restoring both mirrors on an error as the bases mirror is restored
now, then raises item 2's warning on the bridge's first-keep answer through a basis wording of
[`warnScaleKept`](../../R/bartcore.R): "fewer than two distinct values of the response are in the likelihood,
so the scale of forest 2 ('dose')'s basis is kept as it was (said once per sampler)", and under probit and
logistic the row-count lead the creation fallback uses. The docstring and man/dbartsSampler-class.Rd say: burn-in
only; the rows read are those in at the call, which may differ from the rows the response's scale was read at
creation (dec-B362); re-installing the column in force with TRUE is how a caller rescales after the weights or
the mask change (dec-B303).

Bridge. [`bartcore_setForestBasis`](../../src/R_interface_bartcore.cpp) takes a fourth argument, the flag, and
returns TRUE on the sampler's first keep (the `scaleKeptWarned` flag `setOffset` uses, which a sampler of
several forests never sets otherwise, `updateScale` being refused there); its registration's arity goes 3 to
4. A refusal from the engine stays "the sampler refuses this basis for the forest", R having named the
reachable ones.

Engine. R-agnostic, no R name or record:

- [`Sampler::setForestBasis`](../../src/bartcore/sampler.hpp) gains the flag and returns an outcome (refused,
  installed, installed with the scale kept). Under TRUE it decides once, from chain 0, before any chain moves:
  forest f is a map forest of kind numeric that states no sd (else refused); the rows in through
  [`Sampler::leafRowsFrom`](../../src/bartcore/sampler.hpp) over chain 0's weights and the mask; the gaussian
  distinct-values test the response's keep uses
  ([`ResponseModel::lacksDistinctRowsIn`](../../src/bartcore/model.hpp)); then
  [`basisColumnScale`](../../src/bartcore/combiner.hpp) over the incoming column and those rows, a
  `BasisScaleError` read as a keep (open call 1). Then the existing install on every chain, then the restatement
  on every chain.
- `Chain::restateForestBasisScale(f, scale)`, beside [`Chain::setForestMapSd`](../../src/bartcore/chain.hpp):
  sets the forest's `basisScale` and clears `basisScaleFromAllRows`, resolves the entry again through
  [`forestMapEntry`](../../src/bartcore/combiner.hpp) with the forest's own default statement, and resets its
  leaf scale through `mapLeafScale`. A numeric coefficient's prior is N(0, 1) whatever the scale, so no
  amplitude prior moves. False, writing nothing, off a numeric default map forest.
- [`SamplerBase`](../../src/bartcore/facade.hpp)'s `setForestBasis` virtual takes the flag and the outcome:
  `--preclean` on every commit. The shipped header has no such entry
  ([`DBARTS_C_API_HASH`](../../inst/include/dbarts/dbarts.h) does not move).
- No RNG is consumed. No state block or `stateFormatVersion` changes: a state holds no scale, and a reload
  reads the record.

Out of scope: `setData(updateBasisScale = )`, which dec-B265 names, waits for `setData` on several forests,
refused today ([What this leaves for later, and why each is an addition](forest-multiplier.md#what-this-leaves-for-later-and-why-each-is-an-addition)).

### The long form (slice 1)

R only. [`forest`](../../R/model.R) gains `leaf.prior = NULL` after `sd`, and folds it at the call: NULL is
absent; a `normal()` whose `prior.sd` is one number becomes `sd`, the leaf prior itself not kept, so every reader
of a forest specification ([`validateForestKnobs`](../../R/model.R), [`forestStatements`](../../R/model.R),
[`resolveForestSpreads`](../../R/dbarts.R), the single-forest rule) sees the one statement and none changes.
Refused at the call, by name:

- both given: "forest()'s 'sd' and 'leaf.prior' are one statement, the forest's size; give one";
- `normal(k = )`, `normal(sd = invchi())`, `normal()` with neither (open call 2), and any other class: "forest()'s
  'leaf.prior' states the forest's size as normal(sd = ); a k, a prior on the sd or another leaf model is the
  fitting function's 'leaf.prior'".

`normal` joins the vocabulary of the forest arguments: `forestConstructors` and the `forests` entry of
[`FOREST_ARGUMENT_VOCABULARIES`](../../R/model.R), and the vocabulary [`processHit`](../../R/formulaTerms.R)
evaluates a formula term's arguments in, so `leaf.prior = normal(sd = 30)` works in a formula, in `forests =
list()` and in `$setLeafPrior(forests = )`, as `fixed()` does. man/forest.Rd gains the `leaf.prior` item,
pointing to `sd`; the examples keep `sd = `, the spelling dec-B266 has the help show. NEWS, in 1.0-0's forest
item: one clause.

## Interaction with the two pending plans

setmodel-parameters-not-structure (96bf027c):

- Its refusal of a raw leaf spec out of step with the resolved slots (slice 1, item 4) reads the model's
  `leaf.prior`. The long form never reaches a model as such: on one forest it builds the model `sd = ` builds,
  pinned `identical()` here, and on several it is the forests' record's `sd`. No conflict.
- `updateBasisScale` writes the data object's scale, not the model or the control, and its slice 4 routes a
  sampler with amplitudes through the tree-prior and residual installs, never the leaf-scale write. So a
  `$setModel` after an update leaves the re-derived prior in force. Whichever lands second adds the pin:
  update, `$setModel(s$model)`, then the reader's `sd` and the next draws unchanged, and a copy equal.
- Its open door for per-forest tree priors (`$setTreePrior(forests = list(forest(base = , power = )))`) meets
  the long form: folded, `forest(leaf.prior = normal(sd = 2))` there is an `sd` and is refused as `sd` would be.
  The added `normal` in the forest vocabulary is harmless to it.
- Both add facade virtuals (its slices 3 and 4): a textual conflict in facade.hpp, `--preclean` on each.

forests-shared-total (f4b9256b):

- `forest(leaf.prior = normal(sd = v))` on a forest with no basis is that forest's own size and takes it out of
  the share, as `sd = v` does (its D4); folded, its sharer test reads it unchanged. A restated sd on a sharer
  through `$setLeafPrior(forests = list(forest(leaf.prior = normal(sd = v))))` meets its refusal as `sd` does.
  Whichever lands second pins both.
- Where no forest has a basis, its item 4 has the fitting function's `leaf.prior = normal(sd = )` state the
  group's total, while `normal(sd = )` inside a `forest()` states one forest's own size. That is the one place
  one constructor means two sizes by where it is written; its help and this item's must say so side by side.
  Its D8 gives sharers a k; `forest(leaf.prior = normal(k = ))` stays refused here, and the message names the
  fitting function, where that k lives.
- `updateBasisScale` touches only forests with a basis, which never share. Its open call 2 (refuse a constant
  basis column at creation) and this plan's open call 1 (keep on a swap) agree with dec-B346: creation strict,
  a swap lenient.
- Both call [`forestMapEntry`](../../src/bartcore/combiner.hpp): it adds the sharer count to the signature, and
  this item adds a call site, which it then extends.

Order of landing, recommended: this item first, then forests-shared-total, then setmodel-parameters' slice 4,
which its plan already puts after forests-shared-total; setmodel's slices 1 to 3 are independent of this item
and may land in any order with it. This item is the smallest, neutral and file-near both others, so landing it
first leaves each of them one pin to add and one call site to extend, rather than this item re-gating against
two moved engines.

## Slices

### Slice 1: the long form (R surface)

Steps: the formal and the fold in [`forest`](../../R/model.R); the vocabulary in three places; man/forest.Rd;
NEWS; tinytest in test-forest-arguments.R.

Tests: at every door (formula term, `forests` list, `$setLeafPrior(forests = )`, on a single forest and on
three) the long form gives `identical()` model, control record and train draws to `sd = `; on one forest also to
the fitting function's `leaf.prior = normal(sd = )`; the reader says `sd.stated` TRUE; each refusal by its
message, counted with `withCallingHandlers`; `leaf.prior = NULL` is absent.

Mutants: the fold dropping the value; both given taken; `normal(k = 3)` taken as no statement; `normal` left out
of the formula term's vocabulary.

Gates: the touched test files, the four reproducibility snapshot files unchanged, `lintr::lint_package()`,
`air format --check .`, `tools/check-rc-codoc.R`, `tools/check-win-drift.R`, `tools/check-doc-freshness.R`,
each on its own exit status; `inst/NEWS.Rd` parses. CI: the full suite, check and sanitizers. RNG: neutral.

Cost: ~220 lines; 1.5 agent hours to build, 1 to review and fix.

### Slice 2: updateBasisScale (engine, bridge, R)

Steps as in [updateBasisScale (slice 2)](#updatebasisscale-slice-2), then the help, NEWS and the design note:
[2. A numeric column's scale and its record](../design/forest-multiplier.md#2-a-numeric-columns-scale-and-its-record)
gains the row for a swap with the argument, and section 7's bullet goes.

Tests and gates are below. RNG: neutral for every main call; the new path held to a creation twin.

Cost: ~540 lines; 4 agent hours to build, 1.5 to review, 2 for the gate battery and the fix round.

Both: these estimates run 1.5 to 2 times low, as plans here do; budget ~11 to 15 agent hours in all.

## Tests and gates

tests/cpp (slice 2), a fixture of two forests, the second numeric, two chains:

- After a swap with the flag, every chain's leaf scale and calibration equal those of a chain created on the
  new column, bitwise; FALSE leaves every chain bitwise as today's swap.
- Rows: a weight of zero and a mask entry each take their rows out of the read, compared with the sd over the
  rows in by hand.
- Keeps: no row in, one row in, one response value over the rows in (gaussian), a constant column (open call
  1): the column installed, the scale in force, the outcome says kept.
- Refusals write nothing on any chain: a stated forest, a factor forest, a forest with no basis.
- Trees, coefficients and the other forest's calibration unchanged across the call.

tinytest, test-forest-basis-scale-update.R, small n:

- The twin. A sampler created on `dose`, swapped to `3 * dose` with TRUE before any run, against one created on
  `3 * dose`: the reader's `sd` and `basis.scale` `identical()`; the train draws of a run `identical()`.
- After a run: `sd` is the default over the new column's sd; `sd.stated` stays FALSE; the coefficient and the
  trees unchanged at the call; `data@basis.scale` the new scale.
- Rows: created with rows of weight zero, `setWeights` to all ones, the same column re-installed with TRUE:
  the scale over every row; under probit, through the mask.
- `copy()`, and `storeState`, `saveRDS` and `readRDS`, after an update: next draws `identical()` to the live
  sampler's.
- Each refusal by message, the sampler and both mirrors untouched after it (twin); the keep's warning once per
  sampler, counted with `withCallingHandlers`.
- A swap without the argument pinned as today (the scale kept).

Mutants (each must fail a named expectation): the read over every row; the record not written (the copy test);
the restatement on chain 0 alone (the two-chain tests/cpp case); the scale read from the column before the
install; a stated forest not refused; the keep not taken for one response value; the warning every call;
`sd.stated` flipped to TRUE.

Gates, on the slice's own libraries, `--preclean` on every commit:

- tests/cpp plain and under address and undefined sanitizers; the R-loaded sanitizer run over the new test file
  (`basisColumnScale` becomes reachable from a setter).
- The equivalence trio on the reference build, `compare --bitwise --strict-coverage`: 55 of 55
  (equivalence-32418b6e), 19 of 19 (bcf-equivalence-20eed7ed, whose numeric_swap_write covers a swap without the
  flag) and 11 of 11 (multinomial-equivalence-80b1c8d4), identical, no `max |z|` line; the four snapshot files
  unchanged.
- benchmarks/R/forest-prior-exact.R and bcf-exact.R at quick, as smoke on the touched map path.
- The full tinytest suite in one process; the lint chain and checks of slice 1; `R CMD check --as-cran` from a
  clean staged tarball.
- No `bench-sampler.R` compare: the code runs once per call, never per sweep.

## Calls made

1. `updateBasisScale` goes after `updateState`, as `setWeights` put its new arguments after the ones it had.
   Alternative: before it, as `setOffset(offset, updateScale, updateState)`; `$setForestBasis` is new in 1.0-0, so
   either is free, and after keeps every call written today.
2. The rows read are those in at the call, positive weight and kept by the mask, under every family (dec-B303),
   read as a leaf covariate's re-derivation reads them. Alternative: the rows the response's scale was read at
   creation, which dec-B303 does not say and which never follows a weight change.
3. Too few rows, or one response value over them, keeps the scale in force with one warning per sampler,
   sharing the flag the response's keep uses (dec-B364, dec-B465, dec-B479). Alternative: every row, as creation
   falls back, which jumps a mixture component's prior to the whole data's.
4. TRUE on a factor forest is refused by name, by dec-B265's reason for a stated one: no column-derived default to
   update. Alternative: a silent no-op, kinder to a loop over forests.
5. The decision is made once before any chain moves, and the swap and the restatement land together or not at
   all. Alternative: install, then restate, and roll back on a refusal.
6. R names the refusals from the reader's row; the engine refuses the same cases as a backstop. Alternative:
   an engine outcome per reason, worded by the bridge.
7. The new scale overwrites `data@basis.scale`, the record every re-creation reads. Alternative: a second
   record for the updated scale, which nothing else needs.
8. `setData(updateBasisScale = )` waits for `setData` on several forests.
9. The long form folds into `sd` at the call, so no reader changes and both pending plans read it for free.
   Alternative: keep `leaf.prior` on the specification and resolve it at fit time, which every reader of a
   forest would have to learn.
10. The long form takes `normal(sd = number)` alone, on one forest and on several (dec-B266: "the only leaf
    prior a forest can take is a normal with an sd"). Alternative: on a model of one forest, any leaf prior the
    fitting function takes, k and `linear()` among them, which would let a single-forest reader's `normal(k = 2)`
    be handed back, at the cost of two meanings by model shape.
11. Messages that name a forest's `sd` keep saying `sd` when the caller wrote the long form. Alternative: carry
    the spelling on the specification for messages alone.
12. A design-note amendment, not a new note: the rule of the scale does not change, only who may re-derive it.

## Open calls

1. A swap with `updateBasisScale = TRUE` whose column is constant over the rows in. dec-B265 (2026-10-06): "a
   constant column is refused wherever a default would have to be derived from it, at creation with no stated
   sd and on an update". Since then dec-B346 took a constant or zero column on a swap, dec-B364 kept the
   response's scale where a re-derivation finds one value, and dec-B465 kept a leaf covariate's numbers where a
   re-derivation finds one value, each refusing at creation.
   - (a) Refuse, as dec-B265 rules: the whole call, the column not installed, the message naming the constant.
     Cost: an outer sampler imputing a treatment that, during burn-in, makes every row in treated stops with an
     error mid-run, unless it catches it and calls again without the argument.
   - (b) Recommended. Install the column, keep the scale in force, one warning per sampler, as dec-B364 and
     dec-B465 do. Revises dec-B265 in part; creation keeps its refusal. Cost: a mistaken constant column is
     reported once, not stopped.
2. `forest(leaf.prior = normal())`, with no sd. The forest-multiplier plan listed "`normal()` for the default
   back" beside the long form; no ruling names it. A forest stated at creation holds no scale, so putting its
   default back means deriving one.
   - (a) Recommended. Refused for now, the message saying to leave `leaf.prior` out for the default. Both
     meanings stay additions. Cost: no way to return a stated forest to its default.
   - (b) "The default": at creation the same as no statement; through `$setLeafPrior(forests = )` it returns a
     stated forest to its default, from the recorded scale where one exists and from the column in force over
     the rows in otherwise. Cost: ~250 more lines (an engine restatement and its rows, as slice 2), and one
     more engine-tier review.
   - (c) No statement everywhere, as `forest()` is. Cost: closes the spelling to (b) later.

## Not verified

- The twin's train draws being `identical()` assumes creation draws nothing from the column but the scale; if
  it does, the twin falls back to the reader's numbers and the copy test, and the implementer reports it.
- The engine steps are from reading the code, not a prototype.
- That `scaleKeptWarned` is never set on a sampler of several forests is read from dec-B362's refusal, not probed
  on every setter.
- 0.9-34 was probed for the setters' signatures only.
- The pending plans were read at 96bf027c and f4b9256b, before their reviews; their calls may move.
