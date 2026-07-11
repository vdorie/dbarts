# cran-readiness

agent: opus
rng: neutral - every fix must be behavior-neutral (equivalence
  identical); doc/example/build-flag changes only. Anything that
  cannot be fixed neutrally is REPORTED, not fixed.
budget: findings-driven; expect Rd/example/DESCRIPTION touch-ups plus
  compiled-code warning hygiene.

## Goal

Find the CRAN rejections early: a clean R CMD check --as-cran run on
the built tarball, a full ASAN+UBSAN pass over the compiled core, and
a triaged findings list - mechanical fixes landed, judgment calls
reported. Release itself waits (VD); this buys the lead time.

## Context

- R CMD check --as-cran has never run cleanly in the session harness
  (it breaks there); run it in a plain shell against a built tarball
  (R CMD build first). The season added many Rd examples
  (survivalProbabilities, samplePriorPredictive, growFromRoot,
  n.grow.sweeps) that have never seen the example-timing limits.
- CRAN's extra checks run gcc-ASAN/clang-UBSAN on Linux; this machine
  is arm64 macOS, where tests/cpp sanitizes directly and the R-loaded
  .so sanitizes with more ceremony. A sanitizer caught a real defect
  this season (the callback landing), so the pass has teeth.
- Sanitizer CI gating is a separate follow-up decision (TODO soft
  call), not this item. Manual pass only.

## Steps

1. R CMD build, then R CMD check --as-cran on the tarball in a clean
   shell. Triage every NOTE/WARNING/ERROR: fix the mechanical ones
   (Rd nits, cross-refs, example runtimes via smaller inputs or
   \donttest, DESCRIPTION fields, non-ASCII, stray files in the
   tarball via .Rbuildignore), report the judgment ones.
2. Compiled-code warning hygiene: build src/ with CRAN-like strict
   flags (-Wall -Wextra -pedantic) and clear what is clearable
   without behavior change.
3. Sanitizers: tests/cpp under -fsanitize=address,undefined
   (mandatory, full run); the R-loaded package under UBSAN and, if
   macOS tooling permits, ASAN over the full tinytest suite
   (best-effort; document exactly what ran and what remains for
   CRAN's Linux images).
4. Findings memo appended to this plan: what was fixed, what was
   found clean, what needs VD, what could not be exercised on macOS.

## Verification

- check --as-cran status recorded verbatim (target: no ERROR/WARNING,
  NOTEs enumerated and justified or fixed).
- Sanitized tests/cpp: zero reports. Sanitized R suite: zero reports
  on whatever subset ran.
- Behavior neutrality: full tinytest 2727+ / 0; equivalence vs
  equivalence-ac6ec2c.rds 22/22 IDENTICAL; tests/cpp clean rebuild.
