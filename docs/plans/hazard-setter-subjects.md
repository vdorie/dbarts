# hazard-setter-subjects: a hazard fit's row setters count its subjects

Status: LANDED 2026-10-25 (33731694); implements dec-B490 and dec-B491, recorded as dec-A232 ([Built](#built),
[Fix round](#fix-round), [setPredictor (dec-B491)](#setpredictor-dec-b491), [Second fix round](#second-fix-round),
[Third fix round](#third-fix-round)).

agent: one opus implementer (engine flag, bridge, R text, tests, help, notes); one opus review that runs mutants.
rng: draw-changing only where a setter on a hazard fit re-derived a covariate at fewer than two subjects in and
now keeps it. NEUTRAL, bit for bit, for every other call; no recorded scenario has a hazard setter at one subject
in: expected 55, 19 and 11 identical against equivalence-32418b6e, bcf-equivalence-20eed7ed and
multinomial-equivalence-80b1c8d4, no re-record.
window: before the merge (TODO `hazard-setter-subjects`).
budget: ~100 lines of package code plus tests; built at 92 code lines added and 32 removed, comments aside;
dec-B491 34 and 14, the fix round 151 and 43 (30 of them the bridge's `using` list sorted).

## Goal

`setWeights` and `setActiveRows` with `updateStandardization = TRUE` on a hazard fit where fewer than two subjects
are in the likelihood after the call keep every linear or gp leaf covariate's center and scale, and a gp leaf's
default lengthscale, in force, the period column included, under one warning, once per sampler, that counts
subjects. Two or more subjects in re-derive each covariate by its own values over their period rows (dec-B465).

## Context

The ruling is dec-B490, revising that call of dec-A231; it extends dec-B489's subject count from creation and a
deriving `setData` to the row setters, beside dec-B465 (the setters keep a covariate with too little) and dec-B486
(one call's covariates kept together). The build it extends is [leaf-one-row-fallback.md](leaf-one-row-fallback.md)
([Fix round](leaf-one-row-fallback.md#fix-round)): R counts a hazard fit's subjects
([`likelihoodUnits`](../../R/bartcore.R), [`tooFewUnitsIn`](../../R/bartcore.R)) and tells the engine a flag at
creation (`SamplerOptions::tooFewUnitsIn`) and at `setData` (`DataReplacement::tooFewUnitsIn`), which reads every
row in place of the rows in. The setters' re-derivation is
[`Sampler::rederiveStandardization`](../../src/bartcore/sampler.hpp), reached from the bridge's
`rederiveAfterInstall` ([`bartcore_setWeights`](../../src/R_interface_bartcore.cpp),
[`bartcore_setActiveRows`](../../src/R_interface_bartcore.cpp)); it already keeps every covariate where fewer than
two rows are in, and the bridge reports a keep on the holder's one "warn now" flag (`firstKeep`). R words the keep
([`warnRowsKept`](../../R/bartcore.R), [`warnScaleKept`](../../R/bartcore.R)). A probit hazard fit's
`setWeights` passes its 0/1 weights to `setActiveRows`; a logistic one's weights are positive counts, so only the
mask in force leaves its rows out there. Rule 9 of [response-scale-rows.md](../design/response-scale-rows.md)
states the setters' rule.

## Steps

1. Engine: `Sampler::rederiveStandardization` takes the host's flag, as `setData` does, reading the rows as too
   few in (every marked covariate kept) where it is set; the facade virtual follows.
2. Bridge: the two setter entries take a fifth argument, TRUE or NULL as the create and `setData` entries take it,
   passed to the re-derivation.
3. R: the setters count the subjects in after the call, over the new weights or mask and the one in force, and pass
   the flag; the keep warning counts subjects where it was passed.
4. Tests, help (the setters' docstrings and the sampler help), rule 9.

## Verification

The probe (center 2.78 and scale 1.39 over subjects 4 and 5; subject 5 alone by either setter keeps them, one
warning counting subjects); a two-library compare of setter cases, moving only the targeted ones; the engine
tier's gates.

## Built

Built in 2eb23614. [`rowSetterUnits`](../../R/bartcore.R) counts a hazard fit's subjects in after the call through
[`likelihoodUnits`](../../R/bartcore.R), which gains a `weights` argument: `setActiveRows` over the weights in force
and the new mask, `setWeights` over the new weights and the mask in force. Where fewer than two of more than one
are in ([`tooFewUnitsIn`](../../R/bartcore.R)) and the call re-derives, the setter hands TRUE to a fifth argument
of [`bartcore_setActiveRows`](../../src/R_interface_bartcore.cpp) or
[`bartcore_setWeights`](../../src/R_interface_bartcore.cpp), as the create and `setData` entries take it, which
`rederiveAfterInstall` passes to [`Sampler::rederiveStandardization`](../../src/bartcore/sampler.hpp); there it
reads the rows as too few in, `LeafRows{nullptr, true}`, as `setData` does on `DataReplacement::tooFewUnitsIn`, so
every covariate keeps through the path no row in already takes. The bridge reports the keep on the holder's one
flag (`firstKeep`), so the warning is said once per sampler, shared with every other keep, and again by a copy or
a reload. [`warnScaleKept`](../../R/bartcore.R) words it from the count: "one subject is in the likelihood (1 of 60
subjects is in it), so the centers and scales of "x1" and "period" are kept as they were (said once per
sampler)". A probit hazard fit's `setWeights` passes its 0/1 weights to `setActiveRows`; a logistic one takes its
own door, where weights are positive counts and only the mask in force leaves a subject out. dbarts.h is
unchanged; the direct `.Call` sites in four test files and mask-redraw-exact.R take the fifth argument.

Probe (scratch, 60 subjects, breaks 0:5, a linear leaf on period, then on x1 and period; subjects 4 and 5 in: period
center 2.78 and scale 1.39). On 0aa94644 `setActiveRows` or `setWeights` to subject 5 alone with
`updateStandardization = TRUE` read period at 3 and 1.58, silently with period alone and under x1's one-value
warning with both. Built, both keep 2.78 and 1.39 under the one subject text; creation with subject 5 alone still
reads every row (2.33 and 1.25).

What moves, by a two-library compare against 0aa94644 over 113 setter cases (probit and logistic hazard fits,
linear on x1, on period and on both, gp on both, subject 5, subjects 4 and 5, every subject and no subject in, both
setters with `updateStandardization` TRUE and FALSE, the logistic `setWeights` door under a mask, a one-subject fit,
gaussian and probit fits at 0, 1, 2 and 50 rows in): 81 identical in record, warnings and draws. The 32 that move
are all fewer than two subjects in with `updateStandardization = TRUE`; the record and the draws move only in the
12 where the period column was read again (period, both and gp at subject 5), and in the other 20 (x1 alone at one
subject, every leaf at none) only the warning's text changes, to the subject count.

Calls made, each with its alternative:

- The engine is told, by a new argument of `Sampler::rederiveStandardization` and its facade virtual. Alternative:
  the bridge keeps every covariate itself without calling the engine, as it does where the response's scale is
  kept (dec-B486), with no engine change. Measured: that route moves the draws of all 24 linear-leaf cases at
  fewer than two subjects in, the 16 already kept today among them, since the engine's keep rebuilds a linear
  leaf's fits under the numbers held and the bridge's skip does not; told, they stay bit for bit what they were,
  and the setters take the trigger creation and `setData` take.
- The count runs only where the call re-derives. Alternative: at every setter call, which the bridge reads only for
  a re-derivation; it would cost a pass over the period column and, while setPredictor still took a period column
  out of the layout, refuse every later setter call on such rows (corrected in the [Fix round](#fix-round), which
  refuses that column where it is handed in).
- The warning counts subjects and names every leaf covariate, all kept. Alternative: the one-value text,
  "fewer than two distinct values of the leaf covariate", which is false of the period spread over one subject's
  periods.
- No subject in takes the subject text too, "no subject is in the likelihood (0 of 60 subjects are in it)", as
  creation's no-subject text does; the record and draws there are unchanged. Alternative: the one-value text kept
  at none, the engine keeping everything there already, so that a hazard fit's setters would count subjects at one
  and not at none.
- Where `updateScale = TRUE` keeps the response's scale in the same call, dec-B486's text names the response.
  Unreached on a hazard fit, whose probit and logistic links read no response scale. Alternative: the subject
  text there.
- Not extended: `setPredictor` and `updatePredictor` with `updateStandardization = TRUE` re-derive over the rows in
  at that call too, and on a hazard fit at one subject in still read the period column again over its rows
  (probe: 3 and 1.58, x1 kept); dec-B490 names the row setters.

Tests. test-leaf-standardization-rows.R, the "either count" pin restated: linear and gp on x1 and period, both
setters, subject 5 alone keeping the record, the engine's reader agreeing, finite sweeps; said once by either
setter and again by a copy; period alone, x1 alone and no subject in, each in its text; the logistic `setWeights`
door under a mask; two subjects in re-derived over their nine rows under probit and logistic, and with x1 at one
value over them x1 kept under its one-value text and period read again; a one-subject fit re-derived silently.
tests/cpp `testLeafStandardizationRows`: the flag keeps every covariate over rows in that spread one.

Mutants, each alone, installed with `--preclean`, failures in test-leaf-standardization-rows.R (670), tests/cpp
for the engine's: the setters' count dropped, 22; the threshold at three subjects, 7 (the two-subject pins and the
logistic door's); the flag not passed at the `setActiveRows` door, 18, and at the `setWeights` door, 2; the bridge's
`setActiveRows` entry ignoring it, 18, and its `setWeights` entry, 2; the engine ignoring it, 20, tests/cpp 1; the
facade not forwarding it, 20, tests/cpp 1; the warning counting rows, 12; the warning said on every keep, 8; the
mask in force ignored at the `setWeights` door, 2; the count without `tooFewUnitsIn`'s guards, 2 (the one-subject
fit's pin and a probit fit's one-row setter text). All twelve killed. Equivalent by construction and not run: the
`setWeights` door counting over the old weights (a logistic hazard fit's weights are positive, so only the mask
leaves a subject out there) and the count run where the call does not re-derive (the bridge then reads no flag).

Gates, on the package code of 2eb23614 from `--preclean` installs in a private library: tests/cpp 412 ok plain and
under ASAN and UBSAN, no diagnostic; `tinytest::test_package` at its default `ncpu` in one process, 268 files,
31042 results, none failed; the 43 commands of exact-gates.yaml at `quick` on the shipped build, all passing; on a
`--enable-reference-build` install the four snapshot files (14, 3, 7 and 3, none failed) and the bitwise compares,
55 of 55 identical against equivalence-32418b6e (`--strict-coverage`, 0 skipped), 19 of 19 against
bcf-equivalence-20eed7ed and 11 of 11 against multinomial-equivalence-80b1c8d4 on every channel, no "max |z|"
line, and both cross-host compares within tier 1; lint (0 lints), air, check-rc-codoc (48 methods),
check-win-drift, check-doc-freshness and verify-anchors clean; `R CMD check --as-cran --no-manual` from a tarball
built from a clean export, the Date field's NOTE alone. Not run: the R-loaded path under ASAN (no new numerics: the
flag sends a re-derivation down the path no row in already takes), stan4bart (it calls neither changed entry, both
the R bridge's), the x86-64 leg.

## Fix round

The review of dc6fc00c (one opus reader, FIX FIRST) ran an 848-case two-library compare (240 moved, all hazard
fits with `updateStandardization = TRUE` at fewer than two subjects in, the draws bit for bit after 30 sweeps where
a covariate was already kept) and 13 mutants: 8 killed, 2 equivalent, 3 surviving. Built in 9d9a83d8, on
wt/hazard-setpredictor beside dec-B491, as the orchestrator ruled each finding.

1. A fit of one subject was exempt everywhere: `tooFewUnitsIn`'s `n > 1` let subject 5 alone as a fit re-read its
   period over two of its five rows (1.5 and 0.71) silently at `setActiveRows`, say the one-value text with none
   in, and read two rows at a deriving `setData` and by hand. Narrowed to dec-A231's analogy, one row of a one-row
   store being every row in: a fit of one subject is exempt only where every row is in, nothing then being fallen
   back from. [`likelihoodUnits`](../../R/bartcore.R) reports whether every row is in, and
   [`tooFewUnitsIn`](../../R/bartcore.R) flags fewer than two subjects in unless it is, at every door that counts
   (creation, by hand too, a deriving `setData`, the row setters, `setPredictor`): with some rows out a one-subject
   fit now reads every row at creation and `setData`, "one subject is in the likelihood (1 of 1 subjects is in it),
   so the center and scale of "period" are taken from all 1", and keeps at a setter; with every row in it reads its
   rows, silently, as before.
2. `setPredictor` took a period column out of the layout, after which the setters' count refused with no prefix
   (0aa94644 re-derived there). Every door that writes a hazard fit's training period column now refuses one that
   could leave the layout, whatever `updateStandardization`, before anything moves, in the layout text behind the
   door's name: `setPredictor` whole (the dense and the sparse-bearing route), by column and row by row
   (`"partial"`), and `updatePredictorPerObservationJointly`. A row-by-row install can put in any subset of its
   rows, so there each row is taken at its value in force or its new one
   ([`refuseHazardLayout`](../../R/bartcore.R), [`refusePeriodUpdate`](../../R/bartcore.R)): subject 5's periods 1
   to 5 handed as 1, 2, 1, 2, 3, in the layout as a whole, are refused, since 1, 2, 3, 2 could go in; handed back
   as they are, they go in. The setters' count takes the caller's prefix. The plan's claim that counting at every
   call would change nothing is corrected above.
3. M13 (`n > 2L`, a two-subject fit exempt) survived: pinned, subjects 4 and 5 as a fit with a leaf on x2 and
   period, `setActiveRows` to the second subject keeping 2.78 and 1.39 under "1 of 2 subjects is in it".
4. M11 (the bridge keeping every covariate itself) survived: pinned bit for bit, two seeded fits over subjects 4, 5
   and 10 with a linear leaf on x1, 30 sweeps, then one `setActiveRows` to subject 5 with `updateStandardization`
   and the other the bridge entry without the flag, which keeps x1 for its one value: the next ten draws agree.

Notes, all taken: the comments of `likelihoodUnits` and `tooFewUnitsIn` state every door and that at a setter the
flag keeps; the header's parameter names take the source's stems (`tooFewUnits`, six entries); the count is
skipped where the leaf takes no covariate; rule 9 is rewrapped to the note's width; the bridge's `using` list is
sorted.

Calls made, each with its alternative:

- The exemption by whether every row is in (the orchestrator's). Alternative: no exemption, under which a fit of
  one subject with every row in would read every row with a warning at creation, the same rows it reads, and keep
  at a setter what it would read again unchanged.
- A row-by-row door checks every way the column in force and the new values can mix. Alternatives: the new values
  alone, which takes a column whose partial install leaves the layout; refusing any change to the period column
  row by row, which refuses a split of a subject that keeps it.
- The creation text says "taken from all 1" of a one-subject fit, as it says "all 60". Alternative: words of its
  own for one subject (taken in the [Second fix round](#second-fix-round)).
- The count at `setPredictor` reads the rows in force before the call, ahead of the period column's check.
  Corrected in the second fix round: this said the column was checked first and that a count over the new column
  gave the same answer once checked; the count ran first, and a regrouping kept in the layout gave another
  answer. Since the second fix round the column cannot change at that door, so the count is the same before or
  after the install.
- The setters keep their prefix though their count can no longer refuse, every door that writes the period column
  now checking it. Not checked: a test set's period column, which `setTestPredictor` writes, is laid out period by
  period across the test subjects, a contract of its own the subject count does not read.
- `setWeights` on a probit hazard fit passes its 0/1 weights to `setActiveRows`, whose prefix its count takes.

## setPredictor (dec-B491)

Built in 7869e3b0. `setPredictor(updateStandardization = TRUE)` on a hazard fit where fewer than two subjects are
in keeps every replaced leaf covariate's center and scale, and a gp leaf's default lengthscale, the period column
among them, under dec-B490's subject text, once per sampler on the flag the row setters' keeps share; two or more
subjects in re-derive each by its own values. [`bartcoreSamplerSetPredictor`](../../R/bartcore.R) counts through
`rowSetterUnits` over the weights and the mask in force, before anything moves, and hands the flag to both bridge
entries, [`bartcore_setPredictor`](../../src/R_interface_bartcore.cpp) (whole matrix, the dense and the
sparse-bearing route) and [`bartcore_updatePredictor`](../../src/R_interface_bartcore.cpp) (by column), as a sixth
and seventh argument, which pass it to `Sampler::rederiveStandardization`, the engine function dec-B490 changed; no
other engine function is touched, and dbarts.h is unchanged. The method warns through `warnScaleKept` with the
count. A hazard fit that `dbarts()` makes takes the dense route; one with a sparse column takes the other.

Probes (60 subjects, breaks 0:5, subjects 4 and 5 in, then the mask to subject 5): on dc6fc00c the whole matrix kept
x1 under the one-value text and read period at 3 and 1.58, and the period column alone read it at 3 and 1.58
silently; built, both keep 2.78 and 1.39 under the subject text.

What moves, by a two-library compare against dc6fc00c over 209 `setPredictor` cases (probit and logistic hazard
fits, linear on x1, on period and on both, gp on both, subject 5, subjects 4 and 5, every subject and none in, the
whole matrix, period alone, x1 alone, x2 alone and x1 with period, `updateStandardization` TRUE and FALSE, a
one-subject fit, gaussian and probit fits at 0, 1, 2 and 50 rows in): 153 identical in record, warnings and draws;
of the 56 that move, 18 move record and draws (a period column replaced and read again at subject 5, now kept) and
38 only the warning's text (x1 at one subject, any leaf at none). x2 alone, every non-hazard case and dec-B490's 113
row-setter cases are identical; after the fix round both compares are unchanged.

Calls made, each with its alternative:

- The count reuses the row setters' (`rowSetterUnits`), the one rule for every setter. Alternative: a count of its
  own in the predictor setters.
- The flag reaches the engine through both bridge entries' new argument. Alternative: R keeping the replaced
  covariates by passing `updateStandardization = FALSE`, which would bypass the engine's keep path and the
  holder's once-per-sampler flag, as the bridge-only route measured for dec-B490 does.
- The keep is only of the replaced covariates, as the re-derivation already marks them; a replaced column that is
  no leaf covariate keeps nothing and says nothing.

Tests (test-leaf-standardization-rows.R, 710 results then, 750 after the fix round): the whole matrix and the period
column alone, linear and gp on x1 and period, each keeping the record and the reader's numbers under its text, with
finite sweeps; a period-only leaf; x2 alone silent; no subject in; once per sampler with the row setters, either
first; two subjects in, the whole matrix keeping x1 for its one value and reading period again, and the period
column alone read again silently under probit and logistic; the sparse-bearing route. The direct `.Call` sites of
the two entries in four test files take the new argument.

Mutants on the final code (9d9a83d8), each alone, `--preclean`, failures in test-leaf-standardization-rows.R of 750
(tinytest lists up to seven failures a run): dec-B491's, the count dropped 20; the flag not passed at the dense
whole-matrix call 7, the sparse-bearing call 3, the column call 9; the whole-matrix entry ignoring it 10, the
column entry 9; the flag's reader taking NULL for TRUE 67; the warning without the count 9; the mask ignored by the
count 20; the flag on every re-derivation 67. The fix round's: the one-subject exemption dropped 10; the old
exemption restored 11; a two-subject fit exempt (M13) 13, the two-subject pin among them; the bridge keeping
itself (M11) 1, the bit-for-bit pin; the whole-matrix door unchecked 6, the column door 5, the row-by-row door 5,
the joint door 5; the row-by-row door, and the joint door, checking the new values alone 4 and 3; the first row
checked in one of the choices alone 2. All 21 killed. Before the fix round dec-B491's ten were run on 7869e3b0,
all killed. Equivalent and not run: the setters' prefix dropped and the count run for a leaf with no covariates,
neither reachable once every door checks the layout.

Gates, once, on the package code of 9d9a83d8 (the final commit adds this record), from `--preclean` installs in a
private library: tests/cpp 412 ok plain and under ASAN and UBSAN, no diagnostic; `tinytest::test_package` at its
default `ncpu` in one process, 268 files, 31122 results, none failed; the 43 commands of exact-gates.yaml at
`quick` on the shipped build, all passing; on a `--enable-reference-build` install the four snapshot files (14, 3,
7 and 3, none failed) and the bitwise compares, 55 of 55 identical against equivalence-32418b6e
(`--strict-coverage`, 0 skipped), 19 of 19 against bcf-equivalence-20eed7ed and 11 of 11 against
multinomial-equivalence-80b1c8d4 on every channel, no "max |z|" line, and both cross-host compares within tier 1;
lint (0 lints), air, check-rc-codoc (48 methods), check-win-drift, check-doc-freshness and verify-anchors clean;
`R CMD check --as-cran --no-manual` from a tarball built from a clean export, the Date field's NOTE alone. Not
run: the R-loaded path under ASAN (no new numerics), stan4bart (it calls none of the changed entries, all the R
bridge's), the x86-64 leg.

## Second fix round

The review of 9eac74df (the same opus reader, FIX FIRST) closed the first review's four findings, found dec-B491
right on every route and nothing moved that should not (a 1420-case compare, 78 moved, all hazard fits with
`updateStandardization = TRUE` at fewer than two subjects in or a one-subject fit with rows out), and ran 11
mutants, 6 killed. Built in 9630d190, as the orchestrator ruled each finding.

1. A period column handed twice passed the check, which read the first occurrence while the install writes the
   last (`setPredictor(cbind(p, rev(p)), c("period", "period"))`, dense and dgCMatrix), after which the setters'
   count refused; N03 (the column route reading the first column handed) and N16 (the joint door checking its
   first sampler) survived; and a regrouping kept in the layout (subjects 5 and 6 relabeled as one subject's
   periods 1 to 6) was taken, leaving `survivalProbabilities` 59 subjects and partial dependence an error. The
   predictor doors, `setPredictor` whole (both routes), by column and row by row and
   `updatePredictorPerObservationJointly`, now take a hazard fit's period column only unchanged, every occurrence
   of it checked, and refuse any other values before anything moves, behind the door's name: "the period column of
   a hazard fit is its expansion of each subject's time at risk, set with the response, so it can only be handed
   back unchanged; replace the data with setData to change it" ([`refusePeriodUpdate`](../../R/bartcore.R)). The
   check of every mix of values at the row-by-row doors goes; creation, `new()` and `setData` keep dec-A231's
   layout check ([`refuseHazardLayout`](../../R/bartcore.R), one vector again). With the column fixed at those
   doors, the subject count is the same before or after the install; the fix round's two sentences that said
   otherwise are corrected above.
2. A fit of one subject was worded "(1 of 1 subjects is in it)" and "taken from all 1", and at no subject in
   "(0 of 1 subjects are in it)", in the leaf's text ([`fewRowsInLead`](../../R/bartcore.R),
   [`warnRowsFallBack`](../../R/bartcore.R)) and the periods' text (`hazardGridTimes`): now "(1 of 1 subject is in
   it)", "(0 of 1 subject is in it)" and "taken from the fit's one subject".
3. The test set's period column, which `setTestPredictor` writes, stays unchecked: the reviewer found nothing that
   reads it after the fit. No help text.
4. N10 (a one-row column passing the layout check) is pinned on the check itself, a one-row store at period 2
   refused; built by hand it cannot be made cheaply, `dbartsData` stopping first on a one-row store's cut count.
   N11 (the last row out of the layout named) is killed by the existing pins on rows reversed, where several rows
   are out.

Calls made, each with its alternative:

- The predictor doors take the period column only unchanged (the orchestrator's, as the reviewer recommended).
  Alternative: keep a layout check there, every occurrence and every mix checked, under which a regrouping kept in
  the layout still goes in, renumbering subjects behind the response they were expanded from; a row's period says
  which interval its 0/1 response belongs to, so a changed history goes through `setData`, which checks it with
  the response.
- Equality is checked on the values, a missing value counting as a change. Alternative: within a tolerance, which
  a column of whole periods does not need.
- The one-subject wording changes only where n is one. Alternative: "every subject" for any n, a change to every
  hazard fit's text.

Tests (test-leaf-standardization-rows.R, 763 results): at every predictor door, whole matrix with and without
`updateStandardization`, by column, row by row and jointly, the period column reversed and regrouped in the layout
refused, the sampler left as it was (its predictors and the draws of a twin); the column twice with the second
changed, dense and dgCMatrix; the period after x1, changed and unchanged; the joint update with the hazard fit
second; the column handed back unchanged at every door, the joint one with the hazard fit second too; the sparse
route refusing a changed column and taking it unchanged; the one-subject texts at creation by hand, `setData`, the
setters and no subject in, and the periods' text of a one-subject fit at weight 0; the one-row layout check.

Mutants on 9630d190, each alone, `--preclean`, failures in test-leaf-standardization-rows.R of 763 (a door that
lets a changed column in stops the file, a later hand-back of the column in force then being refused): the
whole-matrix door unchecked, stops; the column door, 11; the row-by-row door, stops; the joint door, stops; the
first occurrence checked alone, stops; the column route reading the first column handed (N03), stops at the
unchanged period after x1; the joint door checking its first sampler (N16), stops; the layout check in place of
equality, 16; the singular noun dropped, 6; "the fit's one subject" dropped from the leaf's text, 2, and from the
periods' text, 1; a one-row column passing the layout check (N10), 1; the last row out named (N11), 3. All 13
killed. Both compares (209 `setPredictor` and 113 row-setter cases) are unchanged by this round.

Gates, once, on the package code of 9630d190 (the final commit adds this record), from `--preclean` installs in a
private library: tests/cpp 412 ok plain and under ASAN and UBSAN, no diagnostic; `tinytest::test_package` at its
default `ncpu` in one process, 268 files, 31135 results, none failed; the 43 commands of exact-gates.yaml at
`quick` on the shipped build, all passing; on a `--enable-reference-build` install the four snapshot files (14, 3,
7 and 3, none failed) and the bitwise compares, 55 of 55 identical against equivalence-32418b6e
(`--strict-coverage`, 0 skipped), 19 of 19 against bcf-equivalence-20eed7ed and 11 of 11 against
multinomial-equivalence-80b1c8d4 on every channel, no "max |z|" line, and both cross-host compares within tier 1;
lint (0 lints), air, check-rc-codoc (48 methods), check-win-drift, check-doc-freshness and verify-anchors clean;
`R CMD check --as-cran --no-manual` from a tarball built from a clean export, the Date field's NOTE alone. Not
run: the R-loaded path under ASAN (no new numerics), stan4bart, the x86-64 leg.

## Third fix round

The review of 3a2d9b3d (FIX FIRST for missing tests only, no defect in the package) closed the second review's three
blocking findings and the regrouping, found the hand-back accepted at every door tried and the texts and help
matching, and ran 7 mutants: P06 and P08 killed, Q01 equivalent, P01 with the same end state, P02, P03 and P07
surviving. A 1420-case compare was 1414 identical and 6 moved in the one-subject text alone. Built in f7e925ad, as
the orchestrator directed for each finding.

1. P02: `na.rm = TRUE` in the equality of [`refusePeriodUpdate`](../../R/bartcore.R) survived, a period column with
   one missing value (otherwise unchanged) refused on the slice and taken under the mutant. Pinned at every door,
   whole matrix with and without `updateStandardization`, by column, row by row and jointly, the sampler left as it
   was.
2. P03: a tolerance of 0.5 in the equality survived, `period + 0.3` taken. Pinned at the same doors.
3. P07: `if (subjects && n == 1L)` weakened to `if (n == 1L)` in [`warnRowsFallBack`](../../R/bartcore.R) survived,
   the "taken from the fit's one subject" wording leaking to a fit that is no hazard fit. Pinned: a gaussian fit of
   3 rows, `setData` of one row at weight 0 with `updateScale = TRUE`, says "taken from all 1".
4. P01 (a note): checking only the last occurrence of the column handed twice has the same end state, the install
   writing the last. Pinned all the same, the first changed and the last as it is, dense and dgCMatrix, since the
   rule is that every occurrence is checked.
5. A row-by-row update of the wrong number of values, through `setPredictor`'s `"partial"` route or
   `updatePredictorPerObservationJointly`, got the period refusal in place of "requires one value per observation":
   the equality read the vector past its end as missing, and a longer one by its first rows. `refusePeriodUpdate`
   now returns for a vector that is not the rows of each column named long, which the bridge's length check then
   refuses before anything moves; the whole-matrix and column doors check their length before the period, and a
   matrix is read as it was. Pinned for n - 1, n + 1 and 2n values, the column's own or reversed, both doors,
   nothing moved. A column update handed as one vector, column after column (`c(period, reversed)` for two
   columns), is still refused as a change: pinned, since a skip on any length but one column's takes it.
6. The joint door's refusal opened "updatePredictorPerObservationJointly: " where the bridge's errors there open
   "updatePredictorPerObservationJointly(): ": now the bridge's form. No help text quoted it.

Calls made, each with its alternative:

- The period check leaves a vector of the wrong length to the bridge's length check. Alternative: a length check in
  R ahead of the period check at the two row-by-row doors, a second copy of the bridge's words.

Tests (test-leaf-standardization-rows.R, 791 results): the above, 28 results over the 763.

Mutants on f7e925ad, each alone, `--preclean`, failures in test-leaf-standardization-rows.R of 791: P02 11, P03 9,
P07 1, P01 2 (the dense and the sparse pin); the period check ahead of the length check again 8; a skip on any vector
not of one column's length 7 listed, then the file stops, the changed column in force refusing the later hand-back;
the joint door's old prefix 5. All 7 killed. No compare was run: a call changed by this round is an error before
and after, in its words alone.

Gates, once, on the package code of f7e925ad (the final commit adds this record), from a `--preclean` install in a
private library: `tinytest::test_package` at its default `ncpu` in one process, 268 files, 31163 results, none
failed; lint (0 lints), air, check-rc-codoc (48 methods), check-win-drift and check-doc-freshness clean. Not run: no
C or C++ changed, so tests/cpp, the sanitizer legs, the exact gates, the snapshot files, the bitwise compares and
`R CMD check` stand as of the second fix round.

The review of e5c100cd found the package code correct and one survivor: R03, the `is.null(dim(x))` guard dropped from
the early return of `refusePeriodUpdate`, which lets a sparse-bearing fit's container, handed whole to `setPredictor`,
skip the check, its `length()` counting its fields. Pinned (795 results): the container with the period's dense column
reversed refused, the sampler left as it was, and the container handed back as it is taken; R03 fails 2 of 795.
