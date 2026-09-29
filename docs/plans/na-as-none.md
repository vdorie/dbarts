# na-as-none: NULL for absent on the public surface

agent: sonnet
rng: neutral (defaults resolve to the same values; no draw moves)
budget: ~450 lines (R ~200, tests ~180, man and NEWS ~70)

Status: LANDED 2026-09-29 (0209919c, 97e49b43, cb799bb0, 12f46c05)

## Goal

Apply dec-B132 (public surfaces follow base R: NULL means absent, NA is a
missing value) to every public argument the 2026-09-28 sweep found using NA
for absent or default, per the rulings dec-B136 to dec-B140.

## Rules for every site

- A spelling 0.9-34 documented (NA on a 0.9-34 formal) keeps working for
  one release: warnOnce with a key per site, a message naming NULL and
  tombstoneExpiry, then the NULL behavior. Register it in R/tombstones.R's
  registry like the other tombstones, so the 1.1-0 sweep removes it.
- A spelling that never shipped (new in 1.0-0) is refused now, as a missing
  value, naming NULL.
- Slots may keep NA internally; only the public argument changes.
- Validation: each argument accepts only its documented values; anything
  else is refused with a message naming the argument.

## Sites

1. sigest (dec-B136): NULL default on bart, dbarts, dbartsSpec, xbart;
   explicit NA warns (0.9-34's bart, bart2 documented NA). The retired
   sigma spelling on dbarts and dbartsSpec defaults to NULL and forwards as
   sigest does. bartBT keeps sigest = NA (BayesTree's list). Fix NULL, which
   today fails "missing value where TRUE/FALSE needed".
2. seed (dec-B137): NULL stays the default; NA warns on bart, bartBT,
   dbarts, xbart, dbartsControl (0.9-34 defaulted to NA); refused now on
   dbartsSpec and dbartsValidateComposition. Drop the "NA is accepted ...
   silently" sentences from the five Rd pages. resolveSeedArg's NA path
   becomes the warning path.
3. updateState (dec-B138): NULL default on the twenty sampler methods and
   updatePredictorPerObservationJointly; NULL defers to control@updateState,
   TRUE/FALSE override, anything else refused; explicit NA warns (0.9-34's
   methods defaulted to NA). resolveUpdateState validates.
4. levelGibbs -> treeShift (dec-B139): dbartsControl(treeShift = c("auto",
   "always", "never")), stored for the bridge as today's tri-state (the
   bridge and engine keep their LevelGibbsMode; only the R surface and the
   slot name, if the slot is public in the Rd, change - keep the slot name
   if renaming it would touch the bridge, and say which). Remove levelGibbs
   from cgm() and dart(), their S4 slots and validity, and spec.R's copy to
   the control; fix spec.R's comment that calls it "the categorical-split
   level Gibbs step". Development-only: no tombstone. Update
   man/dbartsControl.Rd, dbartsPriors.Rd, docs/design/level-fibre.md's
   surface paragraph, benchmarks that pass levelGibbs, and tests.
5. nbinom(dispersion = NULL) (dec-B140): NULL estimates, a positive value
   fixes; NA refused (development-only).
6. run(numBurnIn = NULL, numSamples = NULL) (dec-B140): NULL takes the
   control's; explicit NA warns (0.9-34 documented "missing or NA").
7. dbartsControl(n.samples = NULL) (dec-B140): the slot keeps NA as "not
   set"; explicit NA warns.
8. normal/linear/gp sd and scale, dart rho and update.delay (dec-B140):
   NA refused (development-only).
9. Tombstones (dec-B140): xbart(sigma = ) warns and forwards to sigest
   (0.9-34 formal); $setTestPredictor, $setTestPredictorAndOffset and
   $setTestOffset accept updateState and ignore it with a once-per-session
   warning (0.9-34 formals; test data is not stored state).
10. gaussian(link = identity) unquoted (dec-B140): read a symbol link by
    name as stats::gaussian does.

## Sister packages (lockstep)

bartCause (branch dbarts-1.0) R/bcf.R passes dbartsControl(seed = NA) when
unseeded; stan4bart (branch bartcore) R/mvbart.R passes dbarts(seed =
NA_integer_). Each changes to NULL; their own seed formals stay. Report the
exact lines; the orchestrator commits them.

## Tests

One file, inst/tinytest/test-na-as-none.R: per site, the NULL default,
the NA warning or refusal, the validation refusal, and that the resolved
value equals today's for the default call. Existing tests asserting NA
behavior are updated. tombstone registry test stays green.

## Verification

Full tinytest (unwrapped); lint gates per CLAUDE.local.md; exact gates
quick; R CMD check --as-cran --no-manual; the sister packages' suites
against the change (bartCause dbarts-1.0, stan4bart bartcore) where they
run locally.

## Agent-made calls

Warning keys per site; the treeShift slot storage.

## Landing note (2026-09-29)

Landed as 0209919c (the ten sites, test-na-as-none.R, man and NEWS),
97e49b43 (the negative-binomial exact gate passes NULL for an estimated
dispersion), cb799bb0 (review fixes) and 12f46c05 (two wording slips).
Sister packages in lockstep: bartCause 3a90eb4 on dbarts-1.0 (an unset seed
reaches dbartsControl as NULL, and an explicit bartc(seed = NA) is kept out
of the forwarded call) and stan4bart b2dc25c on bartcore (an unseeded
mvbart component passes NULL).

The warn-or-refuse split follows the rule above checked against 0.9-34,
which corrects several premises in the rulings; dec-A118 records it. NA is
refused on dbarts's sigest and seed, xbart's sigest, dbartsControl's seed
and dbartsSpec's sigest; the test-data setters take no updateState. NaN is
refused by name wherever the change touched. treeShift is stored in the
kept levelGibbs slot, since the bridge reads that name.

Review found two warning leaks on default calls (the hurdle components
forwarded an NA seed; the front door rebuilt a control from NA slots),
NaN read as absent at several sites, and NEWS and Rd text describing
spellings 0.9-34 did not have; all fixed in cb799bb0.

Gates, macOS arm64: full tinytest 10074 pass, 0 fail; test-na-as-none.R
238; nine default fits and a seeded hurdle fit draw bitwise as on the
parent; exact gates quick 25 of 25; R CMD check --as-cran: the expected
NOTE only; air, lintr, rc-codoc, win-drift, doc-freshness clean. bartCause
774 tests and stan4bart 565 pass against the change with their fixes, with
no warnings from seeded or unseeded bartc and mvbart calls.
