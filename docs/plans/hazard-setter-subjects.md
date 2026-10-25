# hazard-setter-subjects: a hazard fit's row setters count its subjects

Status: BUILT 2026-10-24 on wt/hazard-setter, off bartcore 0aa94644, not landed; implements dec-B490
([Built](#built)).

agent: one opus implementer (engine flag, bridge, R text, tests, help, notes); one opus review that runs mutants.
rng: draw-changing only where a setter on a hazard fit re-derived a covariate at fewer than two subjects in and
now keeps it. NEUTRAL, bit for bit, for every other call; no recorded scenario has a hazard setter at one subject
in: expected 55, 19 and 11 identical against equivalence-32418b6e, bcf-equivalence-20eed7ed and
multinomial-equivalence-80b1c8d4, no re-record.
window: before the merge (TODO `hazard-setter-subjects`).
budget: ~100 lines of package code plus tests; built at 92 code lines added and 32 removed, comments aside.

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
- The count runs only where the call re-derives. Alternative: at every setter call; the bridge reads the flag only
  for a re-derivation, so it would cost a pass over the period column and change nothing.
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
