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

- No-ops: [R/dbarts.R:1012-1018](https://github.com/vdorie/dbarts/blob/0fcea39bfe19bacb136847b302759e68a3e8c64d/R/dbarts.R#L1012-L1018); callers all internal
  ([R/bart.R:469](https://github.com/vdorie/dbarts/blob/0fcea39bfe19bacb136847b302759e68a3e8c64d/R/bart.R#L469), [R/bart.R:491](https://github.com/vdorie/dbarts/blob/0fcea39bfe19bacb136847b302759e68a3e8c64d/R/bart.R#L491), [R/bart.R:494](https://github.com/vdorie/dbarts/blob/0fcea39bfe19bacb136847b302759e68a3e8c64d/R/bart.R#L494), [R/bart.R:634](https://github.com/vdorie/dbarts/blob/0fcea39bfe19bacb136847b302759e68a3e8c64d/R/bart.R#L634), [R/bart.R:656](https://github.com/vdorie/dbarts/blob/0fcea39bfe19bacb136847b302759e68a3e8c64d/R/bart.R#L656), [R/bart.R:659](https://github.com/vdorie/dbarts/blob/0fcea39bfe19bacb136847b302759e68a3e8c64d/R/bart.R#L659); [R/partialDependence.R:17](https://github.com/vdorie/dbarts/blob/0fcea39bfe19bacb136847b302759e68a3e8c64d/R/partialDependence.R#L17), [R/partialDependence.R:208](https://github.com/vdorie/dbarts/blob/0fcea39bfe19bacb136847b302759e68a3e8c64d/R/partialDependence.R#L208), [R/partialDependence.R:439](https://github.com/vdorie/dbarts/blob/0fcea39bfe19bacb136847b302759e68a3e8c64d/R/partialDependence.R#L439);
  [R/rbart.R:762](https://github.com/vdorie/dbarts/blob/0fcea39bfe19bacb136847b302759e68a3e8c64d/R/rbart.R#L762), [R/rbart.R:826](https://github.com/vdorie/dbarts/blob/0fcea39bfe19bacb136847b302759e68a3e8c64d/R/rbart.R#L826)). No revdep dependency - dbarts controls both ends.
  Caveat: check whether any revdep calls the methods on the sampler
  object directly (grep the revdep sources during the sweep); if one
  does, keep the methods and delete only the internal calls.
- Offset sync: testUsesRegularOffset flips implicitly - construction
  sets it ([R/data.R:175-254](https://github.com/vdorie/dbarts/blob/a2802563422608c15e77e364053a7399a78fab91/R/data.R#L175-L254)), setOffset re-propagates when TRUE
  ([R/bartcore.R:175-240](https://github.com/vdorie/dbarts/blob/a2802563422608c15e77e364053a7399a78fab91/R/bartcore.R#L175-L240)), any explicit setTestOffset/
  setTestPredictorAndOffset breaks the link ([R/dbarts.R:662-745](https://github.com/vdorie/dbarts/blob/a2802563422608c15e77e364053a7399a78fab91/R/dbarts.R#L662-L745)).
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
