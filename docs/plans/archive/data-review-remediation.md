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
  ([[sampler.hpp:1180-1188@bd5279a1]]) re-inlines ColumnStore::setCell
  ([[data.hpp:1173-1180@bd5279a1]]) minus its `hasMissing[j] = 1` line, and bridge
  validateColumnValues ([[R_interface_bartcore.cpp:1505-1512@bd5279a1]]) checks
  categoricals only - an ordinal NA committed per-observation leaves the
  MIA gauge stale: descent routes naCode by a missing bit that
  dropStaleMissingDirections clears, while the NA-free partition kernel
  routes it always-right. Verified by direct read.
- B2: setColumnJournaled's CSC branch ([[data.hpp:1124-1130@bd5279a1]]) is dead
  (refuseViewSampler blocks all mutation on CSC stores,
  [[R_interface_bartcore.cpp:1381-1387@bd5279a1]]) and wrong if ever live (snapshots
  through codeOffsets[j], unused by rank columns). Verified.
- B3: bartcore_setResponse ([[R_interface_bartcore.cpp:2508-2525@bd5279a1]]) checks
  type/length only - NA response reaches the sampler though
  construction rejects anyNA(y) ([[data.R:899-901@bd5279a1]]); setWeights skips the
  R-side checks dbartsData's validity enforces ([[A_class.R:465-478@bd5279a1]])
  while setPredictor/setOffset validate. Verified.
- P1: codeFor's ordinal path ([[data.hpp:300-309@bd5279a1]]) walks cuts linearly;
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
  the setPredictor precedent in [[bartcore.R:111-121@bd5279a1]]; keep the C length
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

## Landings

C1 d9fc90d (2026-07-18). B1: commitObservation routes through setCell
(the dead re-inlined body deleted; setCell gains its one production
caller) - safe because the session's cuts are fixed, so setCell's
requantize reproduces the precommitted code, and the hasMissing mark
is restored; the fix is discriminating (the old body fails all three
new assertions). B2: the CSC branch deleted, comment states
dense-only citing refuseViewSampler. P1: codeFor's ordinal arm is
std::lower_bound - the old loop counted cuts below value, identical
result for every input including ties and NaN (returned early).
One fix round: the implementer's new tests perturbed the shared test
RNG history and broke the BCF grow-from-root snapshot downstream,
initially misdiagnosed as pre-existing - the reviewer reproduced
clean HEAD passing, and both tests now carry the local-generator +
rngState save/restore insulation convention (second occurrence;
the convention is now load-bearing for any new tests/cpp test that
builds a sampler). Gates: suite 3050/0, all three anchors bitwise,
component tests green with the snapshot at its recorded value;
reviewer re-ran the battery after a provenance reinstall. Diff 127.

C2 7fcc635 (2026-07-18). setResponse rejects NA R-side; setWeights
enforces the class-validity trio (length, NA, non-negative);
setSigma, which had NO validation on either side, gains length-1 /
non-NA / positive (Rf_asReal previously truncated a vector silently -
now an explicit error); setOffset gains the NA rejection
construction already enforced. BONUS FIND: setResponse mutated
data@y BEFORE the C length check, so a rejected call corrupted the
R-side slot with no rollback - reordered to install-on-success.
16 new tinytest expectations assert each rejection fires before any
sampler state changes. Gates: suite 3066/0, equivalence 22/22
bitwise, air format clean, no lints. Diff 91. ARC CLOSED; the
refactor-scale findings live in data-store-consolidation.md and
r-ingestion-cleanups.md.
