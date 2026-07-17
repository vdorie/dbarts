# data-review-remediation

agent: Opus (C1), Sonnet (C2)
rng: neutral throughout (C1's quantize change is bit-identical by
  argument AND gate; nothing touches a draw)
budget: C1 ~150 lines; C2 ~120 lines

## Goal

Fix the three verified bugs and the one first-tier performance finding
from the 2026-07-17 four-reviewer data-handling review (abstractions,
performance, C++ factorization, R surface). The larger refactors the
review proposed are recorded as TODO entries, not taken here.

## Findings addressed

- B1 (two reviewers independently): UpdateSessionImpl::commitObservation
  (sampler.hpp:1180-1188) re-inlines ColumnStore::setCell
  (data.hpp:1173-1180) minus its `hasMissing[j] = 1` line, and bridge
  validateColumnValues (R_interface_bartcore.cpp:1505-1512) checks
  categoricals only - an ordinal NA committed per-observation leaves the
  MIA gauge stale: descent routes naCode by a missing bit that
  dropStaleMissingDirections clears, while the NA-free partition kernel
  routes it always-right. Verified by direct read.
- B2: setColumnJournaled's CSC branch (data.hpp:1124-1130) is dead
  (refuseViewSampler blocks all mutation on CSC stores,
  R_interface_bartcore.cpp:1381-1387) and wrong if ever live (snapshots
  through codeOffsets[j], unused by rank columns). Verified.
- B3: bartcore_setResponse (R_interface_bartcore.cpp:2508-2525) checks
  type/length only - NA response reaches the sampler though
  construction rejects anyNA(y) (data.R:899-901); setWeights skips the
  R-side checks dbartsData's validity enforces (A_class.R:465-478)
  while setPredictor/setOffset validate. Verified.
- P1: codeFor's ordinal path (data.hpp:300-309) walks cuts linearly;
  a full predictor swap is O(n*C*p) ~ 5e9 comparisons at n=1e6, more
  than a sweep. The loop counts cuts strictly below value, exactly
  std::lower_bound over the same sorted doubles with the same
  comparison - the returned code is IDENTICAL for every input incl.
  ties (draws untouched by construction; gates confirm).

## Commits

C1 (Opus, C++): B1 - commitObservation routes through setCell (delete
  the re-inlined body; setCell gains its one production caller, or if
  the call shape does not fit, add the hasMissing line and delete
  setCell - implementer's choice, one cell-writer either way). B2 -
  delete the CSC branch; the doc comment states dense-only and cites
  the bridge refusal. P1 - lower_bound in codeFor's ordinal arm.
  tests/cpp: per-observation NA into a clean ordinal column sets
  hasMissing and routes consistently (session commit then a descent /
  partition agreement check); a codeFor tie-and-boundary check vs the
  old linear scan transcribed in-test.
C2 (Sonnet, R + bridge): B3 - setResponse rejects NA (R side, matching
  the setPredictor precedent in bartcore.R:111-121; keep the C length
  check); setWeights gains the length/non-negativity/NA checks the
  class validity states; setSigma/setOffset audited for the same gap
  and aligned. One shared R helper if it falls out naturally; do not
  build a framework. tinytest for each rejection.

## Constraints

- All gates bitwise: equivalence-ac6ec2c 22/22, bcf-99205ee,
  multinomial-8c2b5fc, suite 3050+ (grows with new tests). Any
  divergence = stop, report, never re-record.
- No dbarts.h / inst/include change. C2's R touches pass
  air format --check . and lintr.
- Serialize C1 then C2.
- Machine may carry other load; do NOT run bench-sampler (the P1
  speed claim is verified at the next quiet window by the orchestrator
  via the setPredictor arms).

## Verification

Per commit: R CMD INSTALL --preclean .; tests/cpp make clean && make &&
./test_bartcore; full tinytest; the three equivalence compares
(bitwise). C2 additionally: air format --check . exits 0.
