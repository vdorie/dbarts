# cross-family-state-install: a state goes only into a sampler of the family that stored it

Status: PLANNED.

agent: opus implementer, one (bridge, engine, tests); opus reviewer.
rng: by call sequence.
- NEUTRAL, bit for bit, for every sampler that installs no state, for every state installed into a sampler
  of the family that stored it (by `setState`, `copy` or a reload, with the new record or without it), and
  for every warm start.
- A changed result, not a generator matter: a state of another family is refused where today it is
  installed, and a Student-t, logistic or negative-binomial sampler refuses a latent block holding a value
  that is not positive and finite. A refusal draws nothing.
Proved by the bitwise gates below (no baseline scenario, exact gate or snapshot file installs a state across
families), by a seeded digest of twenty kinds of sampler on the base and slice builds, and by the new
tests' identity with a twin.
window: pre-release, before 1.0-0. Serial with [monotone-unforced-refusal.md](monotone-unforced-refusal.md)
and [leaf-conversions.md](leaf-conversions.md), which edit the same bridge file, chain.hpp and the sampler's
manual page in other functions and items; shares no file with
[factor-column-update-forms.md](factor-column-update-forms.md) or [written-surface.md](written-surface.md).
Recommended order: after monotone-unforced-refusal, before leaf-conversions or between its pushes.
budget: ~430 lines (bridge ~50, engine ~35, R none, tests/cpp ~70, tinytest ~230, manual, design notes,
comments and TODO ~45), upper figure 650. Plans have run 1.5-2x low; the scratch build under Context came to
79 lines of C++.

## Goal

A state says which response family stored it, and `setState`, `copy` and a reload refuse a state of another
family by name, before anything is touched, where today most such states install and some leave a sampler
whose fits are not finite or whose next sweep does not return. A state that carries no such record installs
as it does today, except that no sampler takes scales or Polya-Gamma variates that are not positive. A warm
start keeps taking the trees of any family.

## Context

All numbers were run on the tip's build (shipped mode): 60 rows, three predictors, 5 trees, one chain,
twenty kinds of sampler on the same predictors (the eight families; Student-t and negative binomial with
their scalar fixed and drawn; the hazard model under each link; a variance forest; two forests under
gaussian, probit and logistic; linear, gp and monotone leaves; probit with a monotone leaf). Each ordered
pair ran in a process of its own under a 20-second limit ("no return"): an install, then 23 sweeps.

- What in a state says which model stored it. Every state has class `bartcoreState` and the same seven
  attributes (cut points, two counters, two versions, two digests); a Student-t state adds `weights.zero`.
  Per chain, `latents` is present for every family but gaussian and multinomial, `thresholds` for ordinal,
  `resid.df` and `shape` where drawn. Nothing separates a probit state from a logistic or a hazard one.
- What the refusal keys on. [`Chain::stateIsValid`](../../src/bartcore/chain.hpp) judges counts, block
  shapes and each tree against the data. Of the family it asks only that a sampler without latents is
  given none, that a Student-t or negative-binomial sampler is given some, and for an ordinal sampler's
  thresholds and a negative-binomial one's shift. [`Chain::setState`](../../src/bartcore/chain.hpp) then
  copies the block in ([`TResponse::restoreLatents`](../../src/bartcore/model.hpp) and its siblings).
- The eight families against each other, `setState`, one forest, constant leaf (rows: the state; columns:
  the sampler; "refused" is `state is not consistent with this sampler`):

  | state | gaussian | Student-t | probit | logistic | ordinal | nbinom | aft | multinomial |
  |---|---|---|---|---|---|---|---|---|
  | gaussian | - | refused | runs | runs | refused | refused | runs | refused |
  | Student-t | refused | - | runs | runs | refused | runs | runs | refused |
  | probit | refused | not finite | - | no return | refused | refused | runs | refused |
  | logistic | refused | runs | runs | - | refused | refused | runs | refused |
  | ordinal | refused | not finite | runs | no return | - | refused | runs | refused |
  | nbinom | refused | runs | runs | runs | refused | - | runs | refused |
  | aft | refused | not finite | runs | no return | refused | no return | - | refused |
  | multinomial | refused | refused | refused | refused | refused | refused | refused | - |

  Of the 56 pairs 31 are refused, 18 install and run, 3 leave fits that are not finite and 4 a sweep that
  does not return. Among gaussian, Student-t, probit and logistic: 4, 6, 1 and 1.
- All twenty kinds, 380 pairs across kinds: 309 refused (196 in the words above, 68 `malformed cut points
  in bartcore state`, 32 `bartcore state is missing required block 'tree.params'`, 13 others), 56 run, 8
  not finite, 7 no return. Of the 328 pairs of different families 265 are refused, 48 run and those 15
  break, the hazard and the two-forest forms of probit into logistic as the plain pair does. After each
  of the 309 refusals the stored state is byte for byte the one before and the next three draws are an
  untouched twin's.
- How it breaks. Every breaking pair puts real-valued latents (probit, ordinal, aft; 25 to 58 rows
  negative here) where the sampler holds precisions: as Student-t scales, as low as -3.2, they give fits
  that are not finite within three sweeps; as Polya-Gamma variates the next sweep does not return. Of the
  48 that run, 30 leave the sampler holding another family's latent block as stored.
- `copy` and a reload install the object's `state` field through the same entry
  ([`copy`](../../R/dbarts.R), [`getPointer`](../../R/dbarts.R)), so they meet another family's state only
  when the field was assigned by hand. They refuse the same 309 pairs in the same words and break on the
  same 15 and one more. A refused `copy` leaves the original running as its twin; a reloaded object
  raises the error at its first use and at its second.
- The warm start reads no latents ([`readWarmStartState`](../../src/R_interface_bartcore.cpp)): it took
  134 of the 380 pairs, each ran 23 finite sweeps, and none held the state's latent block afterwards.
- Pairs of one family, 52 here. 44 differ in forests, leaf model, a variance forest or the hazard's rows
  and are refused for that. Of the eight that differ in a setting or a leaf constraint, six install
  whatever the state: Student-t with the degrees of freedom fixed and drawn, and negative binomial with
  the shape fixed and drawn, each way; a monotone leaf's state into the plain sampler, probit and
  gaussian. `setState` returns `TRUE`, the sampler's trees and latent block afterwards are byte for byte
  the state's, what it holds fixed stays its own, and 200 sweeps are finite. The two reverse installs go
  into a monotone sampler when the state's leaf values are in order.
- A sampler's own precision blocks are positive and finite in 52 of 52 corners (before any sweep, under a
  mask, at weight zero, with count weights, on a response of all zeros); the smallest value was 0.019.
- One existing pin installs across families:
  ["a gaussian state leaves a probit sampler's sigma, pinned at 1, where it is"](../../inst/tinytest/test-state-not-model.R)
  (the sampler it builds is logistic).
- The rule below, tried on a scratch copy of the tip: 372 of the 380 pairs are refused, the 328 of
  different families by name, and the 8 that install are the pairs of one family above; `copy` and a
  reload agree on every pair; nothing fails to return. With the record removed from each state 325 are
  refused, the 15 breaking pairs among them, and 55 run. Own-kind installs with and without the record
  give identical draws, 20 of 20; the digest under Verification is equal on the two builds; tests/cpp
  passes unchanged (349); the tinytest suite fails at the one pin only.
- Consumers (read only). stan4bart splices its per-chain samplers' chains into one state object,
  attributes kept, for a sampler of the same model; that installs on the scratch build. No other consumer
  stores a state.

## The rule

What is stored. Every state carries a top-level attribute `family`, one string: the response family the
sampler runs, as the `family` argument spells it - `gaussian`, `student`, `probit`, `logistic`, `ordinal`,
`nbinom`, `aft` or `multinomial`. A hazard sampler writes its link's family, the model it runs on its
expanded rows. The leaf model, the forests and a variance forest are not in it.

The install, by `setState`, `copy` and a reload alike:

1. The record is the sampler's family: installed by today's rules.
2. The record is another family's: refused with `state is not consistent with this sampler: its family is
   "probit" and the sampler's is "student"; to start one fit from another's trees use installTrees`. Asked
   after the class, the format version and the chain count, and before the digests, the cut points and the
   blocks are read, so past those three a state of another family is named as that whatever else about it
   differs.
3. The record is present and is not one string: `malformed family in bartcore state`.
4. No record, as on a state stored before this: installed by today's rules.
5. Whatever the record says, a Student-t, logistic or negative-binomial sampler refuses a latent block
   holding a value that is not finite and positive, with today's `state is not consistent with this
   sampler`. The other families' latents are real numbers and are not judged.
6. Rules 2 and 3 are raised before the engine is given the state, and rule 5 with the engine's other
   validity checks: after any of them the sampler, its stored state and its generators are as they were.

Every pair of different families is refused, the 18 of the table that run today included. Pairs of one
family that differ in a fixed value or a leaf constraint stay accepted, on the evidence under Context.
`copy` raises the refusal and returns nothing; a reloaded object raises it at each use until a state of
its own is assigned to the field, as today for any state it refuses. The warm start (`installTrees`,
`warm.start`) is not changed and does not read the record: it is the route across families.

## Constraints

- An install into the state's own family draws and returns what it does now, with the record or without.
- State format: an added top-level attribute; `formatVersion` and the oldest readable version stay 1,
  under the rule at [`stateFormatVersion`](../../src/R_interface_bartcore.cpp). Every state gains it, so
  none is byte for byte what it is now; the gates compare draws.
- No facade virtual changes ([`ResponseModel`](../../src/bartcore/model.hpp) gains one) and
  [`DBARTS_C_API_HASH`](../../inst/include/dbarts/dbarts.h) does not move; `--preclean` all the same.
- The message starts with today's words, so a caller matching them still matches. The lines the mutation
  battery anchors in the state reader stay as they are.
- No R code changes, and no NEWS entry: the bartcore state is new in 1.0-0.

## Steps

1. Engine. [`ResponseModel`](../../src/bartcore/model.hpp) gains a const reader saying whether a stored
   latent block is one the family can hold, true by default;
   [`TResponse`](../../src/bartcore/model.hpp), [`LogisticResponse`](../../src/bartcore/model.hpp) and
   [`NBResponse`](../../src/bartcore/model.hpp) answer "every value finite and positive".
   [`Chain::stateIsValid`](../../src/bartcore/chain.hpp) asks it where it checks the block's length.
2. Bridge. One static function names a sampler's family from its
   [`SamplerShape`](../../src/bartcore/facade.hpp): `multinomial` where `supportsCountsMutation`,
   `student` where `carriesResidualDf`, else by [`ResponseFamily`](../../src/bartcore/model.hpp), as
   [`dbarts_sampler_family`](../../src/C_interface.cpp) reads the first and the last.
   [`storeState`](../../src/R_interface_bartcore.cpp) writes it;
   [`setState`](../../src/R_interface_bartcore.cpp) reads it by name before the digests and applies rules
   2 and 3 through the error it already accumulates, printing the stored name to a bounded width.
3. tests/cpp, beside [`testStateValidation`](../../tests/cpp/test_state.cpp),
   [`testStateRoundTripLatents`](../../tests/cpp/test_state.cpp) and
   [`testStateRoundTripStudentT`](../../tests/cpp/test_state.cpp): a logistic, a Student-t and a
   negative-binomial state with one latent set to 0, -1, NaN and infinity is refused and the state read
   back is the one before; 1e-300 installs; a probit state with a negative latent installs. Fails today.
4. tinytest, a new file `test-state-family.R`: the eight families on one predictor matrix, a hazard
   sampler under each link, and two forests under probit and logistic.
   - What a state carries: the family's name for each of the eight, `probit` and `logistic` for the
     hazard states, the response family for two-forest and monotone states. Fails today: absent.
   - Each of the 56 ordered pairs, the hazard pair and the two-forest pair: an error naming both
     families, the stored state byte for byte the one before, the next three draws a twin's. Fails today
     in 25 of the 56 (installed) and by its message in 31. The draws are asked only of a sampler that
     refused, so the file returns on the tip, where 7 of these installs would not.
   - One family, another setting: the six pairs that install whatever the state, `TRUE`, trees and
     latents the state's, three sweeps finite; a plain state into a monotone sampler installs or is
     refused for its leaf values, never for its family. Holds today.
   - A state from before, the attribute removed: into its own family it returns and draws what the same
     state with the attribute does, each family, and a Student-t state into a logistic sampler installs
     and runs (rule 4; both hold today). A probit state into a Student-t and into a logistic sampler is
     refused in the plain words, sampler untouched. Fails today: installed.
   - Malformed records: a number, two strings, `NA_character_` and a raw vector give `malformed family`;
     an empty string and a name no family has are refused naming it. A latent block edited by hand: 0,
     -1, `NaN` and `Inf` refused for Student-t, logistic and negative binomial; a probit state edited the
     same way installs. The sampler untouched after each refusal.
   - `copy` with another family's state assigned to the field: the error by name, the original's next
     draws its twin's. A reload: the error at first use and at the second; its own state assigned back,
     it runs. Two chains spliced from two samplers of one family install; another family refuses them.
     `installTrees` takes a probit state into a Student-t sampler, with the attribute and without.
   In [test-state-not-model.R](../../inst/tinytest/test-state-not-model.R) the pin under Context installs
   its gaussian state with the attribute removed and says so; its expectation stands.
5. Mutations (Verification): apply, install with `--preclean`, run, report the counts, revert, `touch`.
6. Records.
   - Manual, [`dbartsSampler$setState`](../../man/dbartsSampler-class.Rd). The `newState` item opens "For
     `setState`, a state object previously produced by this sampler, or by another of the same response
     family over the same data". In the paragraph on restoring, after "a refused restore leaves the
     sampler exactly as it was": "A state records the response family of the sampler that stored it, and
     a state of another family is refused, naming both: its trees and latent variables mean something
     only under the family they were drawn in. To start a sampler from a fit of another family, use
     `installTrees`, which takes the donor's trees and, where this sampler draws them, its `sigma`, `k`
     and DART state, and none of its latent variables. A state that carries no such record, as one stored
     before the record existed, is judged by its contents. A Student-t, logistic or negative-binomial
     sampler also refuses a state whose latent variables are not all positive and finite." After the
     sentence on assigning the field directly: "`copy` and a reload install the field by the same rules,
     so an object whose field holds a state it would refuse cannot be copied, and after a reload raises
     that refusal at each use until a state of its own is assigned."
   - Design record: the table under
     [What a state carries](../design/state-not-model.md#what-a-state-carries) gains the row (`family`:
     which family's chain this is; model, by name; compared, never installed) and the note a dated
     paragraph with the rule and the table above. [public-surface.md](../design/public-surface.md) lists
     `family` with the added attributes and stops saying a state must hold a binary sampler's latents (a
     gaussian state installs into a probit sampler today, and after this when it carries no record).
   - TODO: the entry `cross-family-state-install` names this plan.

## Verification

- `R CMD INSTALL --preclean` into the slice's own library; `tests/cpp` passes, clean under ASan and UBSan;
  the full tinytest suite; the new file under ASan on the R-loaded path (the bridge prints an R string
  into a fixed buffer).
- On a reference build: the four `test-reproducibility-*.R` files pass unchanged, and the three compares
  are bitwise with `EQUIVALENCE_CORES=2`, every scenario reporting identical draws, counted per scenario
  with no `max |z|` line, at the counts the last landing recorded (55, 15 and 11). Nothing is re-recorded.
- Every gate `.github/workflows/exact-gates.yaml` lists, in `quick` mode, unchanged: a fit's state
  carries one attribute more, and no gate installs a state across families.
- One script on the base and slice builds digesting, for each of the twenty kinds, the draws of a run, a
  restore into a second sampler with the value `setState` returns, a copy and a reload, and the stored
  states with the new attribute removed. Equal.
- The 400 ordered pairs of the twenty kinds through `setState`, `copy` and a reload, each in a process of
  its own under a time limit: every one returns with finite fits.
- Mutations, each expected to fail the named test:
  - the writer omits the attribute, or names a Student-t sampler `gaussian` or a multinomial one by its
    engine family: tinytest "what a state carries" and those families' pairs;
  - the reader does not compare: the 56 pairs, the hazard and two-forest pairs, `copy` and the reload;
  - the reader compares after the engine's install: "byte for byte the one before" and the twin's draws;
  - the reader refuses a state with no record, or skips the type check: "a state from before",
    "malformed records";
  - the floor is removed, or admits zero: tests/cpp step 3, tinytest "a probit state into a Student-t
    and into a logistic sampler" and "a latent block edited by hand";
  - the floor is asked of every family: "a probit state edited the same way installs" and every probit,
    ordinal and aft restore in the suite.
- `lintr::lint_package()`, `air format --check .`, `Rscript tools/check-rc-codoc.R .`,
  `Rscript tools/check-win-drift.R .`, `Rscript tools/check-doc-freshness.R .`,
  `Rscript benchmarks/R/mutation-battery.R verify-anchors`, each on its own exit status; `R CMD build`
  with every vignette rebuilt and `R CMD check --as-cran` on the tarball from a clean copy (man/ changes).
- Not a hot-path change: nothing is added to a sweep. An install scans the latent block once more.

## Out of scope, and where it goes

- Other structure than the family. A DART sampler's state and an ordinary one's install into each other,
  and a monotone leaf's state into a plain sampler; linear and gp leaves, forest counts, a variance forest
  and the ordinal category count are refused by the shape of their blocks, in other words. dec-B254 calls
  all of these structure; no such pair broke, and none is this entry's. The coordinator has the figures.
- Latents that are not finite in a probit, ordinal or aft state: installed today and after.
- The sweep that does not return when a precision is not positive. After this no state can bring one.

## Calls made in planning

- A field, not a refusal alone. The blocks cannot tell probit latents from logistic ones, and the message
  is to name both families. One string at the top level is the smallest thing that does it. It is compared
  and never installed, as the two digests are, so a state still carries no model.
- Every family writes it, gaussian included, and every pair of different families is refused, the ones
  that run today included. The TODO entry's words are "a state whose latents are another family's"; a
  gaussian state has no latents, installs into a probit, logistic or aft sampler today, and one tinytest
  pins that. The rule follows dec-B254's line for `setModel` - a chain's trees and latents mean something
  only under the family they were drawn in - and leaves `installTrees` as the route. The narrower rule, a
  state with no latent block taken by any family, keeps that pin; it is the question in the coordinator's
  notes.
- The name is the user's word for the family: the engine calls a Student-t sampler gaussian and has no
  value for multinomial. A hazard sampler is named by its link, which is what tells its two forms apart.
  An absent record is "not known", never a mismatch, which keeps a state stored before this installing.
- The floor (rule 5) is in the slice. Without it a state with no record, or one edited by hand, can still
  leave the two breaking installs. It is about 25 lines of engine and can be struck without touching the
  rest. No sign check is made on probit, ordinal or aft latents: a state stored before a response change
  installs today and is re-drawn at the next sweep, and such a check would refuse it.
- Order. leaf-conversions edits the state writer and reader at the leaf calibration blocks and model.hpp's
  leaves, in two pushes three times this size: this goes before it or between its pushes, not beside it.
- The tip against the TODO entry's words: its counts hold for the four families it names; the breaking
  pairs are 15 of 380, not two; the coordinator's notes list the rest. Nothing was found done already.
- After review, 2026-10-07. The refusal by name ends "; to start one fit from another's trees use
  installTrees", so the message itself points the way. The stored name is printed up to its first byte
  outside printable ASCII and to 32 bytes, with dots where it was cut: a width counted in bytes cut a
  character in two and left a message that was not valid text. A warm start takes more than the trees - the
  donor's sigma, k and DART state where the sampler draws them - and the help's sentence, first written here
  as "takes the trees alone", says so; run on a gaussian donor into a Student-t sampler (sigma moves to the
  donor's), a drawn k and a DART donor, with the latents, the df and the generator left the sampler's own.
- A hazard sampler's refusal names its link's family, `probit` or `logistic`, and stays so: the expansion to
  person-period rows is made in R and the engine runs a probit or logistic sampler on them, so the bridge
  cannot say "hazard" without R code or a second record, and the record must stay the link's for the same
  states to be accepted. The help says so in one sentence.
- Known, inside the rule. Rule 4's "installed by today's rules" has the floor as its exception in more pairs
  than Context counts: of 784 ordered pairs of 28 kinds without the record, 13 that ran finite on the base
  build are refused by the floor - probit, ordinal and aft states into a Student-t sampler with zero-weight
  rows or a logistic sampler with count weights, where the weights reconciliation redrew the block right
  after the install. None is a restore of a state into its own family.
- Known. The record is one per state, so chains spliced by hand across families are judged by the first
  chain's record: a probit chain and a logistic chain under `probit` install into a two-chain probit
  sampler, and only the floor catches a real-valued chain under a precision family's record. stan4bart
  splices chains of one model only.
