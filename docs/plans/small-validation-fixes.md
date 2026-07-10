# small-validation-fixes

agent: sonnet
rng: item 1 (input-precision-validation) neutral - a pure R-boundary
     warning; no data reaching the engine changes. item 2
     (missing-incorporate-matrix-interface) shifting - the raw-matrix
     "incorporate" path now keeps NA-bearing rows it previously
     complete-cases-dropped; any consumer that was hitting the bug sees
     a real data-composition change (more rows), so its draws move,
     while the documented "incorporate" semantics themselves do not
     change.
budget: R/data.R ~55 lines; man/dbartsPriors.Rd + docs/design/gp-leaves.md
     ~25 lines; two new tinytest files + one existing-test fix ~150
     lines.

## Goal

Two review findings land together: (1) dbartsData warns when a response
is precision-degenerate before the engine silently fits a collapsed
target, plus the matching GP predict-vs-training-fit doc note; (2) the
raw matrix (x.train, y.train) interface actually honors
missing = "incorporate" instead of complete-cases-dropping NA rows
regardless of the argument.

## Item 1: input-precision-validation (review-6, seed 7)

Context: TODO entry input-precision-validation;
docs/plans/architecture-numerical-review.md L172-180 (finding + GP doc
note text); R/data.R response validation block (formerly ending at the
`anyNA(y)` check).

Fix: R/data.R, right after the unconditional `anyNA(y)` stop (now
L733-735), a single check warns (L737-766):

- range(y)/max(abs(y)) < 1e-10 (max(abs(y)) == 0 guarded, since an
  all-zero response is exactly, not approximately, constant - no
  precision loss to report). 1e-10 is ~1e6x double epsilon
  (~2.22e-16), so real precision loss trips it with several orders of
  magnitude of headroom before a merely low-variance response would.
  The warning message reports the distinct-value count for context.

Deliberately NO distinct-value-count check (a first draft had one at
3 <= distinct(y) <= 3, n >= 100; reviewer-cut). The cardinality check
is subsumed by the ratio check for every genuine double-precision
collapse, and what remains is only false positives. Ulp arithmetic:
near magnitude s, representable doubles are spaced ~2.2e-16 * s. For
rounding to collapse n >= 100 draws of a continuous variable to <= 3
distinct doubles, the data range can span at most a few ulps, so
range/scale is ~1e-15 - five orders of magnitude below the 1e-10
threshold, which therefore already fires. Conversely, any y that would
reach a cardinality check WITHOUT tripping the ratio (e.g. 1e6 +
{0,1,2}: ratio 2e-6) has its values separated by ~1e10 ulps -
perfectly distinguishable doubles that the pre-fit standardization
maps cleanly to [-0.5, 0.5] and the engine fits without degradation.
So a cardinality check could only ever fire on legitimately discrete
data (3-level counts/ordinal at n >= 100, e.g. sample(0:2, 200) with
no offset at all), where the "indistinguishable at double precision"
message would be factually wrong.

GP doc note (same TODO entry): predict() at a training row re-krigs the
jitter-free posterior mean via the cached kernel and drawn training
values, while the recorded training fit carries the conditioning
nugget (1e-6, model.hpp) on the diagonal; the two differ by roughly
2e-3 to 3e-3 at training rows only (MCMC always reads the recorded fit
directly). Documented in man/dbartsPriors.Rd's gp() entry and as a new
section in docs/design/gp-leaves.md (before ## Status).

## Item 2: missing-incorporate-matrix-interface (review-4, tier B)

Context: TODO entry missing-incorporate-matrix-interface; R/data.R
dense/data.frame branch (`is.numeric(formula) || is.data.frame(formula)
|| is.factor(formula)`), the `if (!xIsMixed) { ... }` block.

Decision: (a), make the raw-matrix path honor "incorporate" by keeping
NA rows, not (b) error. Reasoning: downstream of this branch, R/data.R
already has interface-agnostic NA handling shared by the formula,
sparse, and mixed-container paths (L733 onward) - `anyNA(y)`
unconditionally stops (L733-735), `xHasNA` is computed generically off
whatever `x` is (L797-812), and `missing == "error"` (L820) is the only
thing that turns xHasNA into a stop. The mixed-container branch already
skips row-dropping entirely and relies on that shared code (comment:
"a mixed container keeps its rows... validated below, like the
sparse-matrix branch above"). The dense branch's own "error" pre-check
(`if (missing == "error") { if (anyNA(y)) stop(...); if (anyNA(x))
stop(...) }`) already runs before the row-drop, so by the time
`y <- y[completeCases]` executes under "error" there is provably
nothing left to drop - the row-drop was live ONLY under "incorporate",
which is exactly backwards. Fix: `completeCases` is `complete.cases(x,
y)` under "error" (unchanged, dead-code-safe) and `rep_len(TRUE,
length(y))` under "incorporate" (R/data.R L611-616); the now-dead
"N row(s) dropped" warning branch is deleted. No hidden trap found: (a)
was cheap, and end-to-end testing (below) confirms the shared
machinery genuinely supports it.

## Verification

- tests/cpp: not applicable, R-only change.
- inst/tinytest/test-data-precision-warning.R (new, 11 checks): fires
  on the huge-offset/tiny-range case and on an exactly-constant
  response; silent on ordinary continuous data, on a moderate-offset
  3-value response (well-separated at double precision - the case
  that motivated cutting the cardinality check), on binary at any n,
  on a small-n categorical response (both kept as regression guards
  against reintroducing a cardinality check), and on an all-zero
  response (0/0 guard).
- inst/tinytest/test-data-missing-matrix.R (new, 16 checks): raw-matrix
  "incorporate" keeps NA rows with no spurious warning;
  missing = "error" still rejects NA-bearing predictors; a missing
  response is still rejected outright (never silently dropped); an
  end-to-end dbarts() fit through the matrix interface reproduces the
  identical formula-interface fit bit for bit
  (`expect_identical(samples.matrix$train, samples.formula$train)`)
  given the same seed - the two interfaces now agree.
- inst/tinytest/test-boundary-inputs.R: the pre-existing constant-
  response boundary case (y.const, L33-40) now correctly trips the new
  item-1 warning; wrapped in expect_warning rather than left to leak
  into the suite's warning log.

## Status (2026-07-10)

Review revision folded in before landing: the first draft's second
warning arm (distinct-value count) was cut per the subsumption
argument recorded under Item 1 above; the ratio check and the 0/0
guard are unchanged, and the message keeps the distinct-value count
as context. Item 2 approved as-is (man/dbarts.Rd already documents
"incorporate" with no matrix caveat, so the fix makes code match
docs). Gates below re-run after the revision from the worktree
(.claude/worktrees/small-validation-fixes):

- R CMD INSTALL . : clean.
- tinytest::test_package("dbarts") -> 2526 pass / 0 fail (2498 baseline
  + 27 new checks + 1 wrapped pre-existing warning).
- air format: clean on all touched files (R/data.R,
  inst/tinytest/test-data-precision-warning.R,
  inst/tinytest/test-data-missing-matrix.R,
  inst/tinytest/test-boundary-inputs.R).
- lintr: 0 lints on the same files.
- man/dbartsPriors.Rd: ASCII-clean (grep verified).

Equivalence vs benchmarks/baselines/equivalence-0e9ccca.rds: 20/21
IDENTICAL, 1 DIVERGED. The "missing" scenario (equivalence.R L218-234)
diverges: max |z| = 72.64, 27/37 summaries |z| > 3, 10 with disjoint
seed ranges. Adjudicated by review as the legitimate shifting class;
the anchor is re-recorded at landing. Not re-run after the item-1
revision: a warning-only change cannot alter any scenario (none has a
low-cardinality or precision-degenerate response), so the 20/21 +
missing-divergence result stands.

Root cause confirmed, not a regression: the "missing" scenario builds
`x` as a bare matrix with 8% NA injected per cell and
`missing = "incorporate"`, i.e. exactly the raw-matrix path item 2
fixes. Under the pre-fix code, `complete.cases(xm, ym)` dropped 286 of
500 rows (57.2%) despite "incorporate" being requested (reproduced
directly: `set.seed(5111L); xm <- matrix(runif(5000), 500); xm[sample
...] <- NA; sum(!complete.cases(xm))` == 286). The fix is the entire
point of item 2: those rows are now correctly retained and routed
through the MIA machinery, so the scenario fits a materially different
(more complete) dataset and its draws necessarily move. All 20 other
scenarios remain byte-identical (same RNG stream) - fix is precisely
scoped. This matches the repo's established "shifting" RNG-class
precedent (sigma-df-zero-weights re-recorded the anchor after a
similar legitimate single-scenario shift; see MANIFEST).

- 2026-07-10: LANDED as de67cbb (squash of wt/small-validation-fixes;
  reviewer stripped a process reference from the R/data.R comment).
  Reviewer re-ran the gates from a private library while the shared
  install stayed on the previous build: tinytest 2526/0; equivalence
  vs 0e9ccca reproduced 20/21 identical with the adjudicated missing
  shift (sigma.mean 1.48 -> 0.68, the fuller dataset fitting better).
  Anchor re-recorded as equivalence-de67cbb.rds (21 scenarios,
  self-compare clean); MANIFEST and the equivalence workflow updated.
