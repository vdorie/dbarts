# cross-family-state-install: a sampler of precisions refuses latents that are not positive and finite

Status: LANDED 2026-10-07 (257ce57d to 96251ac7; dec-A174).

agent: opus implementer, one (bridge, engine, tests); opus reviewer.
rng: by call sequence.
- NEUTRAL, bit for bit, for every sampler that installs no state, for every state installed into a sampler
  of the family that stored it (by `setState`, `copy` or a reload), and
  for every warm start.
- A changed result, not a generator matter: a Student-t, logistic or negative-binomial sampler refuses a
  latent block holding a value that is not positive and finite. A refusal draws nothing.
Proved by the bitwise gates below (no baseline scenario, exact gate or snapshot file installs a state across
families), by a seeded digest of twenty kinds of sampler on the base and slice builds, and by the new
tests' identity with a twin.
window: pre-release, before 1.0-0. Serial with [monotone-unforced-refusal.md](monotone-unforced-refusal.md)
and [leaf-conversions.md](leaf-conversions.md), which edit the same bridge file, chain.hpp and the sampler's
manual page in other functions and items; shares no file with
[factor-column-update-forms.md](factor-column-update-forms.md) or [written-surface.md](written-surface.md).
Recommended order: after monotone-unforced-refusal, before leaf-conversions or between its pushes.
budget: ~380 lines (bridge none, engine ~35, R none, tests/cpp ~70, tinytest ~230, manual, design notes,
comments and TODO ~45), upper figure 650. Plans have run 1.5-2x low; the scratch build under Context came to
79 lines of C++.

## Goal

No sampler takes scales or Polya-Gamma variates that are not positive and finite: `setState`, `copy` and a
reload refuse such a state before anything is touched, where today a state of another family can bring them
and leave a sampler whose fits are not finite or whose next sweep does not return. A state does not say
which family stored it and is not asked (dec-B283): one whose blocks fit the sampler installs as it does
today. A warm start keeps taking the trees of any family.

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
- The rule below, on the build that has it: 325 of the 380 pairs are refused, the 15 breaking pairs among
  them, and 55 run, 47 of them pairs of different families; `copy` and a reload agree on every pair;
  nothing fails to return.
- Consumers (read only). stan4bart splices its per-chain samplers' chains into one state object,
  attributes kept, for a sampler of the same model; that installs on the scratch build. No other consumer
  stores a state.

## The rule

The install, by `setState`, `copy` and a reload alike:

1. A state does not say which family stored it and is not asked: one whose blocks fit the sampler is
   installed by today's rules, and what the sampler then holds of another family's state is not promised.
2. A Student-t, logistic or negative-binomial sampler refuses a latent block
   holding a value that is not finite and positive, with today's `state is not consistent with this
   sampler`. The other families' latents are real numbers and are not judged.
3. Rule 2 is raised with the engine's other validity checks: after it the sampler, its stored state and
   its generators are as they were.

The 18 pairs of the table that run today still run. Pairs of one
family that differ in a fixed value or a leaf constraint stay accepted, on the evidence under Context.
`copy` raises the refusal and returns nothing; a reloaded object raises it at each use until a state of
its own is assigned to the field, as today for any state it refuses. The warm start (`installTrees`,
`warm.start`) is not changed and reads no latents: it is the route across families.

## Constraints

- An install into the state's own family draws and returns what it does now.
- State format: unchanged, and a stored state is byte for byte what it is now.
- No facade virtual changes ([`ResponseModel`](../../src/bartcore/model.hpp) gains one) and
  [`DBARTS_C_API_HASH`](../../inst/include/dbarts/dbarts.h) does not move; `--preclean` all the same.
- The message is today's words, so a caller matching them still matches. The lines the mutation
  battery anchors in the state reader stay as they are.
- No R code changes, and no NEWS entry: the bartcore state is new in 1.0-0.

## Steps

1. Engine. [`ResponseModel`](../../src/bartcore/model.hpp) gains a const reader saying whether a stored
   latent block is one the family can hold, true by default;
   [`TResponse`](../../src/bartcore/model.hpp), [`LogisticResponse`](../../src/bartcore/model.hpp) and
   [`NBResponse`](../../src/bartcore/model.hpp) answer "every value finite and positive".
   [`Chain::stateIsValid`](../../src/bartcore/chain.hpp) asks it where it checks the block's length.
2. Bridge. Nothing: the floor is the engine's, and the state's writer and reader stay as they are.
3. tests/cpp, beside [`testStateValidation`](../../tests/cpp/test_state.cpp),
   [`testStateRoundTripLatents`](../../tests/cpp/test_state.cpp) and
   [`testStateRoundTripStudentT`](../../tests/cpp/test_state.cpp): a logistic, a Student-t and a
   negative-binomial state with one latent set to 0, -1, NaN and infinity is refused and the state read
   back is the one before; 1e-300 installs; a probit state with a negative latent installs. Fails today.
4. tinytest, a new file `test-state-family.R`: the eight families on one predictor matrix, a hazard
   sampler under each link, and two forests under probit and logistic.
   - A state names no family, and an attribute of that name on a state is not read.
   - One family, another setting: the six pairs that install whatever the state, `TRUE`, trees and
     latents the state's, three sweeps finite. Holds today.
   - Latent responses offered as precisions: each of the 15 breaking pairs under Context is refused in the
     plain words, the stored state byte for byte the one before and the next three draws a twin's, and
     installs once its latents are made positive. Fails today: installed. The draws are asked only of a
     sampler that refused, so the file returns on the tip. A Student-t state installs into a logistic sampler.
   - A latent block edited by hand: 0, -1, `NaN` and `Inf` refused for Student-t, logistic and negative
     binomial, on the last chain and on a row the mask has out too; a probit state edited the same way
     installs. The sampler untouched after each refusal.
   - `copy` with a probit state assigned to a Student-t sampler's field: the refusal, the original's next
     draws its twin's. A reload: the error at first use and at the second; its own state assigned back,
     it runs. Two chains spliced from two samplers of one family install. `installTrees` takes a probit
     state into a Student-t sampler.
   In [test-state-not-model.R](../../inst/tinytest/test-state-not-model.R) the pin under Context stands.
5. Mutations (Verification): apply, install with `--preclean`, run, report the counts, revert, `touch`.
6. Records.
   - Manual, [`dbartsSampler$setState`](../../man/dbartsSampler-class.Rd). In the paragraph on restoring,
     after "a refused restore leaves the sampler exactly as it was": "A Student-t, logistic or
     negative-binomial sampler also refuses a state whose latent variables are not all positive and
     finite. Beyond that a state is not checked against the response family of the sampler that stored
     it: one whose contents fit this sampler is installed, and what the sampler then holds of another
     family's state is not promised. To start a sampler from a fit of another family, use `installTrees`,
     which takes the donor's trees and, where this sampler draws them, its `sigma`, `k` and DART state,
     and none of its latent variables." After the sentence on assigning the field directly: "`copy` and a
     reload install the field by the same rules, so an object whose field holds a state it would refuse
     cannot be copied, and after a reload raises that refusal at each use until a state of its own is
     assigned."
   - Design record: [state-not-model.md](../design/state-not-model.md) gains the section
     [A state of another family](../design/state-not-model.md#a-state-of-another-family), with the rule
     and the table as it is under the rule. [public-surface.md](../design/public-surface.md) stops saying
     a state must hold a binary sampler's latents (a gaussian state installs into a probit sampler).
   - TODO: the entry `cross-family-state-install` names this plan.

## Verification

- `R CMD INSTALL --preclean` into the slice's own library; `tests/cpp` passes, clean under ASan and UBSan;
  the full tinytest suite.
- On a reference build: the four `test-reproducibility-*.R` files pass unchanged, and the three compares
  are bitwise with `EQUIVALENCE_CORES=2`, every scenario reporting identical draws, counted per scenario
  with no `max |z|` line, at the counts the last landing recorded (55, 15 and 11). Nothing is re-recorded.
- Every gate `.github/workflows/exact-gates.yaml` lists, in `quick` mode, unchanged: no gate installs a
  state across families.
- One script on the base and slice builds digesting, for each of the twenty kinds, the draws of a run, a
  restore into a second sampler with the value `setState` returns, a copy and a reload, and the stored
  states. Equal.
- The 400 ordered pairs of the twenty kinds through `setState`, `copy` and a reload, each in a process of
  its own under a time limit: every one returns with finite fits.
- Mutations, each expected to fail the named test:
  - the floor is removed, or admits zero: tests/cpp step 3, tinytest "latent responses offered as
    precisions" and "a latent block edited by hand";
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
- The sweep that does not return when a precision is not positive. After this no stored state brings one.
  A value edited by hand to the edge of what a double holds still can: a Polya-Gamma variate of 1e-320, or
  the largest finite double as a scale, passes the floor and breaks the sweep, within one family as across
  two.

## Calls made in planning

- A state is not asked its family (dec-B283). The blocks cannot tell probit latents from logistic ones,
  and nothing is added to a state to say which: one whose blocks fit the sampler installs, a gaussian
  state into a probit, logistic or aft sampler among them, which one tinytest pins. What the sampler then
  holds is not promised; `installTrees` is the route across families that promises something.
- The floor (rule 2) is about 25 lines of engine. It refuses what is plainly broken, by the stored values
  alone. No sign check is made on probit, ordinal or aft latents: a state stored before a response change
  installs today and is re-drawn at the next sweep, and such a check would refuse it.
- Order. leaf-conversions edits the state writer and reader at the leaf calibration blocks and model.hpp's
  leaves, in two pushes three times this size: this goes before it or between its pushes, not beside it.
- The tip against the TODO entry's words: its counts hold for the four families it names; the breaking
  pairs are 15 of 380, not two; the coordinator's notes list the rest. Nothing was found done already.
- After review, 2026-10-07. A warm start takes more than the trees - the donor's sigma, k and DART state
  where the sampler draws them - and the help's sentence, first written here as "takes the trees alone",
  says so; run on a gaussian donor into a Student-t sampler (sigma moves to the donor's), a drawn k and a
  DART donor, with the latents, the df and the generator left the sampler's own.
- Known, inside the rule. The floor refuses more pairs than the 15 that broke: of 784 ordered pairs of 28
  kinds, 13 that ran finite on the base build are refused by it - probit, ordinal and aft states into a
  Student-t sampler with zero-weight rows or a logistic sampler with count weights, where the weights
  reconciliation redrew the block right after the install. None is a restore of a state into its own
  family.
- Known. Chains spliced by hand across families are judged as any state is: a probit chain and a
  logistic chain install into a two-chain probit sampler, and a two-chain logistic sampler refuses them
  by the floor. stan4bart splices chains of one model only.

## Landing note

Landed 2026-10-07 as 257ce57d to 96251ac7 on bartcore, 9 commits, 661 lines added over 10 files against a
planned 430 to 650. One review told to refute: LAND, with nothing that blocked. On its own build it ran
every ordered pair of 48 kinds of sampler with and without the record (4608 installs): every pair of
different families was refused by name and nothing hung or left fits that were not finite; 363 installs
within one model across 55 routes were all taken; states written before the record installed exactly as
before (496 installs and 209 same-family pairs equal on the base and the slice); 2495 refusals each left
the sampler bit for bit as it was; and 729 warm starts were byte for byte the base's, 262 of them pairs
`setState` now refuses.

What the review changed. Four versions of the code that were wrong passed every test and each now fails
one: the record written only on a state of one chain, the floor skipped for the last of several chains,
the floor skipped for rows the mask has out (under which a masked logistic sampler given a zero at an
inactive row installed and its next sweep did not return), and the family compared by a prefix. The
refusal ends by naming `installTrees` as the way to start one fit from another's trees. The help said a
warm start takes the trees alone; it takes the donor's trees and, where the sampler draws them, its sigma,
k and DART state, and none of its latents, and says so. A name too long or not plain text in a hand-edited
record is printed cut, so the message is always valid text.

Known and left (dec-A174): a gaussian state, which holds no latents, is refused by every other family
with the rest; 13 of 784 pairs of a state without the record and a sampler of another family that ran on
the base build, all with weighted recipients, are refused by the floor; a state spliced by hand across
families is judged by its first chain's record; a hazard sampler's refusal names its family as probit or
logistic, the expansion to rows being R's.

Gates at landing, on a clean copy of the rebased tree in a library of its own (shipped mode), run in
series: tests/cpp 351 ok, 0 failed; the full tinytest suite 16654 results, 0 failed, 228 files; lintr no
lints; air, rc-codoc, win-drift, doc-freshness and the mutation battery's anchors clean; `R CMD build`
with every vignette rebuilt and `R CMD check --as-cran` with the Date NOTE alone. By the implementer at
its last commit: the four snapshot files on a reference build; the three bitwise compares at 55, 15 and
11 scenarios, all identical; the exact gates in `quick`, 28 of 28, and the monotone gates; the two test
files clean under ASan; 25 mutations, none surviving. stan4bart, which stores a sampler's state on its
fit, needs no edit: its fits are gaussian or probit and fits saved on the base build restore.

The family record and the refusal by family, which this note describes, were removed under dec-B283
([state-family-record-removal.md](state-family-record-removal.md)); the floor on precisions stands.
