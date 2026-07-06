# repo-modernization

agent: sonnet
rng: neutral
budget: small, recurring

## Goal

The standing hygiene bundle, unchanged in intent from the pre-rewrite
TODO: CI pins current, sanitizer coverage widened, optional tooling
migrations decided when someone wants them.

## Context / Steps

1. CI action pins: recheck on each push; r-lib/actions stays @v2 (the
   latest major); bump when GitHub deprecates a runtime.
2. Sanitizers: add a gcc-asan/valgrind variant beside the existing
   workflow; gate the sanitizers workflow on pull_request once it has
   proven stable on pushes.
3. codecov: optional; wire only if a public badge is wanted (coverage
   measured locally on demand; 88 percent at last run).
4. roxygen2 / NEWS.md migrations: optional, taste-dependent, huge
   one-time diffs; do not start without VD asking (hand-written Rd is
   codoc-clean today; pkgdown renders NEWS.md natively).

## Constraints

- Each sub-item is independent; land separately.
- Out of scope: new CI jobs beyond the sanitizer variant
  (equivalence-ci owns the equivalence job).

## Verification

- Workflows green on push after each change; for sanitizers, a known
  clean run plus one deliberately injected error caught in a scratch
  branch (then reverted).
