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
   proven stable on pushes. DROPPED (2026-07-22): this is a push-based,
   no-PR repo (main is protected by the bartcore integration-branch-
   then-merge pattern), so a pull_request-only trigger would silence
   the sanitizer instead of gating it. Replaced by concurrency
   cancel-in-progress + push paths-ignore, which answers the cost
   concern that motivated this item without losing coverage.
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

## Landing note (2026-07-07, 56e47f6): sub-items 1-2

Pin audit: all 8 workflows checked against current majors; only
sanitizers.yaml lagged (checkout@v5 -> v7); no deprecated runtimes.
Sanitizer widening: the clang job became a fail-fast-off matrix over
the r-hub clang-asan/gcc-asan containers; a new valgrind job runs the
suite under R -d valgrind on the r-hub valgrind container (180-minute
timeout, ERROR SUMMARY grep gate), same install/test/grep shape as
the proven asan job. Triggers unchanged (push + dispatch; the
pull_request gate stays deferred until stable). Gates: YAML parses;
runtime behavior maintainer-deferred to the first push. Sub-item 3
(codecov badge) awaits VD's call; sub-item 4 remains dormant by
design.

Sub-item 3 resolved (VD, 2026-07-07): no codecov badge; coverage
stays local-on-demand. Closed as a documented no.

## Landing note (2026-07-22): CI cost/hygiene pass

Added top-level concurrency (group ${{ github.workflow }}-
${{ github.ref }}, cancel-in-progress) to check-standard.yaml,
cpp-tests.yaml, lint.yaml, and sanitizers.yaml - the push/pull_request
workflows that lacked one. pkgdown.yaml already serializes deploys via
a job-level concurrency group; left as-is. equivalence.yaml, rchk.yaml,
sbc.yaml, revdep-smoke.yaml, and valgrind.yaml are schedule/dispatch-
only (no push or pull_request trigger), so cancel-in-progress does not
apply to them.

push paths-ignore (docs/**, **.md, TODO, benchmarks/baselines/**) was
already present, as an equal-or-broader benchmarks/** ignore, on
check-standard.yaml, cpp-tests.yaml, sanitizers.yaml, and lint.yaml;
verified and left as-is, no duplicate narrower entry added.

Sub-item 2 (pull_request gate for sanitizers) is DROPPED, not
deferred - see Context / Steps above. The cost concern that motivated
it is now addressed by cancel-in-progress + paths-ignore instead,
without silencing the gate on a no-PR repo.
