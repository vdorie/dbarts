# r5-cleanup

agent: sonnet
rng: neutral
budget: ~200 lines

## Goal

Three R-surface leftovers resolved: the startThreads/stopThreads no-ops
and their purely internal call sites are gone; the
testUsesRegularOffset offset sync is explicit or documented; S4
validity errors stop leaking internal class names.

## Context

- No-ops: [[R/dbarts.R:1012-1018@a2802563]]; callers all internal
  ([[R/bart.R:469@a2802563]], [[R/bart.R:491@a2802563]], [[R/bart.R:494@a2802563]], [[R/bart.R:634@a2802563]], [[R/bart.R:656@a2802563]], [[R/bart.R:659@a2802563]]; [[R/partialDependence.R:17@a2802563]], [[R/partialDependence.R:208@a2802563]], [[R/partialDependence.R:439@a2802563]];
  [[R/rbart.R:762@a2802563]], [[R/rbart.R:826@a2802563]]). No revdep dependency - dbarts controls both ends.
  Caveat: check whether any revdep calls the methods on the sampler
  object directly (grep the revdep sources during the sweep); if one
  does, keep the methods and delete only the internal calls.
- Offset sync: testUsesRegularOffset flips implicitly - construction
  sets it ([[R/data.R:175-254@a2802563]]), setOffset re-propagates when TRUE
  ([[R/bartcore.R:175-240@a2802563]]), any explicit setTestOffset/
  setTestPredictorAndOffset breaks the link ([[R/dbarts.R:662-745@a2802563]]).
  Ported classic statefulness; users cannot see or restore the link.
- Validity prefix: "invalid class ... object:" prepended by
  validObject on every user-facing validation error.

## Constraints

- Behavior of the offset sync does not change in this pass; it gets a
  documented name. Making it restorable (a setter) is optional, ask in
  review.
- Out of scope: dbartsSampler method removal beyond the two no-ops;
  any bridge change.

## Steps

1. Delete startThreads/stopThreads and every internal call site.
2. Document the offset-sync rule in man/dbarts.Rd (sampler methods
   section) with one paragraph: when the link exists, what breaks it,
   that it never re-forms.
3. Validity messages: wrap stop() sites (or validity functions) so
   user-facing errors read as sentences without the S4 prefix;
   package-wide sweep of validity strings.
4. Regenerate any tests asserting exact error strings.

## Verification

- Full tinytest (error-message assertions updated).
- R CMD check clean; no exported-method removal warnings (the no-ops
  were R5 methods, not S4 generics - confirm NAMESPACE untouched).
- Equivalence exact (neutral).

## Landing note (2026-07-07, af6cb5c)

Landed: startThreads/stopThreads no-ops removed (methods, internal
call sites, man aliases; revdeps checked - neither stan4bart nor
bartCause calls them); the testUsesRegularOffset sync documented as a
Test offset synchronization subsection in man/dbarts.Rd's details
(the plan's named "sampler methods section" does not exist - closest
fit taken); S4 validity boilerplate stripped from user-facing errors
via newValidated/validateObject wrappers at the 9 construction sites
(nested new() validation propagates, so call sites suffice); one
expect_error pattern updated. Gates: tinytest 2468 ok, equivalence
exact 18/18, codoc clean, lint zero. +42/-36 lines.
