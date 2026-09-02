# data-ownership-3-mutation

agent: opus (the whole plan is R-surface + design/docs; one owner keeps the
  mutation contract and the data@x reconciliation coherent).
rng: neutral - every rewired path passes the SAME raw values to the SAME
  .Call entry points; codes, quantization order, the per-observation scan
  permutation, the installed mask, and every draw are byte-for-byte
  unchanged. Only who holds the R-side vector and whether an EXTRA copy is
  taken changes. Gate: equivalence 22/22 IDENTICAL vs equivalence-ac6ec2c.rds
  at EVERY commit; full tinytest from a preclean install; no snapshot regen.
window: plan 3 of docs/design/data-ownership.md (FROZEN 2026-07-11); builds
  on the plan-1 container (2e2b1c9) and plan-2 frame ingestion (6c90507).
  Plan 4 does views/sharing; plan 5 the sparse-categorical kernel.
budget: ~200-350 lines across ~9 files (R/mixedMatrix.R, R/bartcore.R,
  R/updatePredictorPerObservationJointly.R, R/generics.R, a mutation
  tinytest, man/*.Rd, the design landing note; NEWS deferred to the
  1.0-0 NEWS pass). NO engine (.hpp)
  change, NO bridge change, NO dbarts.h change - the hot mutation path is
  already correct (write-through died at plan 1); plan 3 removes the R-side
  interim copy and settles the recorded OPEN decision.
bench: CORRECTED AT COMMIT 1 - the original success criterion (recover
  setPredictor-accept-n1000-t75 to ~1.0 vs the 32fc7c8 baseline) was based
  on a wrong attribution. Measured: (a) the bench scenario is
  MATRIX-backed (genFriedman), which this plan deliberately leaves on
  copy-modify; (b) even container-backed, the interim helper was ~10 us of
  a ~260 us update - the ~50 us gap vs baseline is structural to the
  plan-1 write-through death (engine-side update into owned codes,
  R-side copy on matrix sources), on BOTH paths, and predates this plan.
  Revised success criteria: the helper is measurably cheaper in isolation
  (measured ~3x, ~10 -> ~3 us/call, container path), reference semantics
  verified at the address level (tracemem), run-* and setPredictor
  metrics DO NOT REGRESS, and the O(p)-spine allocation profile holds.
  The residual ~1.2-1.3x accept-vs-32fc7c8 gap is ACCEPTED as plan-1
  structural cost; attributing its engine-side share precisely is a
  recorded follow-up, not a plan-3 gate. Metrics carry 6-8% noise;
  re-run single flags before belief.

## Goal

Replace the plan-1/2 interim - the R layer copy-modifies data@x on every
accepted predictor mutation (assignIntoPredictorSource, [[mixedMatrix.R:273@f7efebf4]]) -
with reference-install: the accepted column vector is installed into the
container's per-column list BY REFERENCE, no extra copy on top of R's own
copy-on-write. Settle the design's one recorded OPEN item (the engine-owned
mutable-raw flag) and reconcile what data@x means. Draw-neutral throughout.

## Baseline (what holds at plan 2)

- The engine keeps NO predictor matrix. The only owned raw is the
  leaf-covariate gather ([[data.hpp:193@f7efebf4]] gatheredRawColumns, populated by
  setupGatheredColumns:536 for linear/GP leaves, or every column of a data
  handle). No `mutable` designation exists in the engine - the design's
  proposed mutable-raw flag was NEVER BUILT; plan 1 shipped the gather only.
- Off-hot-path raw consumers take raw AT CALL TIME from the live R slot:
  setCutPoints ([[bartcore.R:318@f7efebf4]]), setState ([[dbarts.R:1089@f7efebf4]], [[dbarts.R:1113@f7efebf4]]), getTrees
  saved-tree replay ([[dbarts.R:1243@f7efebf4]]), all via rawPredictorMatrix(data@x)
  ([[utility.R:371@f7efebf4]]). PROT_DATA is the create-time GC anchor for the flat-C
  path; the R5 path passes current data@x fresh, so a data@x reassignment
  never staleness-traps the engine.
- data@x for a frame fit is a dbartsMixedMatrix dense flavor: a per-column
  list ([[mixedMatrix.R:1-20@f7efebf4]]). Accepted mutations keep it CURRENT via
  assignIntoPredictorSource. extract(sampler, "predictors") reads it
  ([[generics.R:823@f7efebf4]]) - already container-authoritative, already current.
- The interim copy, per site:
  - full-column updatePredictor ([[bartcore.R:176@f7efebf4]]): assignIntoPredictorSource
    with rows = NULL does `matrix(as.double(values), ...)` (O(n) coercion) +
    `column[] <- values[,k]` (O(n) copy-on-modify of the shared column) +
    the list/S4 reassembly. Two O(n) touches the engine did not need. THIS
    is the setPredictor-accept regression.
  - per-observation partial ([[bartcore.R:94@f7efebf4]]) and jointly
    ([[updatePredictorPerObservationJointly.R:88@f7efebf4]]): the inherent O(n) merge
    (old-at-!installed, new-at-installed) - no EXTRA copy beyond it.
  - full-MATRIX setPredictor ([[bartcore.R:123-129@f7efebf4]]): already a pointer swap
    (`sampler$data@x <- x`), essentially reference-install already.
  - test predictor column update ([[bartcore.R:353@f7efebf4]]): copy-modifies the x.test
    matrix to assemble the engine's whole-matrix argument.

## THE decision: the engine-owned mutable-raw flag does NOT survive

RECOMMEND: kill it - do not build it. Reference-install makes the container
itself the current raw store, so the sole rationale the design left open (a
CoW-free home for updated columns) is gone: R's copy-on-write over a
per-column list is already O(spine) + O(one column) once the extra copy is
removed. Consequences, all favorable:

- data.hpp gather designations are UNCHANGED: gatheredRawColumns stays
  leaf-covariate-only (an actual per-fit standardization need, model.hpp,
  not a mutation need). No new engine designation, no new width, no new NA
  bookkeeping.
- Memory: zero new engine allocation. A mutable flag would have re-added a
  full owned double column per mutable predictor (the ~400 MB borrow the
  program set out to delete, re-entering by the side door). Killing it keeps
  the ~8x u8 win intact.
- setState / setCutPoints / getTrees replay paths are UNCHANGED: they source
  raw from the live R container (rawPredictorMatrix / extract), which
  reference-install keeps current. There is no second engine copy to keep in
  sync, so no divergence surface and no cross-grid ambiguity beyond what
  plan 1 already closed.
- No dbarts.h change, no stan4bart lockstep delta (getTrees already gained
  its trailing trainingData arg at plan 1; nothing else moves).

The mutable flag was a hedge against R-side CoW cost; plan 2 already cut that
cost from O(n x p) (matrix) to O(columns x n) (per-column list), and plan 3
cuts the accepted-column case to O(spine). The hedge is unneeded.

## data@x after plan 3 (the reconciliation)

The design says "data@x is a creation-time snapshot BY DEFINITION" - but that
sentence was a CONSEQUENCE of the mutable-flag model (engine holds current
raw, data@x frozen, extract reads the engine). With the flag killed and the
container authoritative, that model is superseded. RECONCILE to:

  data@x is the COLLECTED RAW SOURCE: the per-column values R mapped at
  creation, with an accepted mutation swapping the affected column(s) in
  by reference at the mutating call. It is NEVER a representation of the
  engine's quantized state, and the sampler NEVER writes to it during
  sampling - there is no per-iteration maintenance; the only writes are
  the mutation entry points recording what the engine accepted (for the
  per-observation partial case, the old/new merge at the installed mask).
  extract(sampler, "predictors") is the on-demand materializer of the
  historical numeric-code-matrix form over whatever vehicle @x holds
  (matrix, dgCMatrix, either container flavor) - constructed when called,
  never cached, no engine reroute (there is no engine raw to read).

This matches plan-2 reality (assignIntoPredictorSource already keeps data@x
current; extract already reads it). Plan 3 changes the install MECHANIC, not
the invariant. The design's snapshot phrasing is struck in the landing note.
(VD confirmed this contract 2026-07-13: mutation-collection by R reference
is fine; a maintained public view of internal quantized state would not be,
and this is not that.)

## Reference-install mechanic

The install is an R-surface operation; the engine is untouched. For an
accepted column j on the dense-flavor container:

- WHICH object mutates: the R5 sampler's data@x slot is reassigned to a
  container that shares every unchanged column vector BY REFERENCE and whose
  slot j points at the newly supplied vector. The prior container is
  untouched (copy-on-write preserved: a shallow copy holding it diverges -
  [[test-sampler-predictors.R:61-67@f7efebf4]] stays green, no assertion change).
- WHO allocates: the caller (R's `as.double(x)`, or the user's own vector).
  The helper allocates only the O(p) list spine of the shallow container
  duplicate. The engine allocates nothing (it re-quantizes into its existing
  codes buffer, as today). No O(n) column copy, no matrix() coercion.
- HOW R GC sees it: the prior slot-j vector loses its container reference on
  reassignment and is collected when no other binding holds it; shared
  columns stay alive through the new spine. PROT_DATA (flat-C anchor) is
  unaffected - the R5 raw source is always the live data@x passed fresh.
- The supported in-place re-install pattern: a perf caller keeps its own
  column vector, mutates it in place, and passes THAT vector; the helper
  installs it by reference. No aliasing hazard - R's copy-on-write makes it
  safe: the moment either the caller or the container mutates a shared
  vector, R copies, so the two never silently desync. Reference-install
  merely declines to add an engine/helper copy ON TOP of R's; it never
  weakens R's guarantee.

New helper (mixedMatrix.R), replacing assignIntoPredictorSource:

  installPredictorColumns(x, rows, columns, values)   # returns the container

- matrix x: `x[, columns] <- values` for a matrix source stays (a matrix has
  no columnar spine to share; the pointer-swap full-matrix path is already
  reference-install). Column/per-obs updates on a matrix source keep their
  CoW - matrices are the flat-C / hand-coded-matrix case, off the frame hot
  path.
- container x, rows = NULL (full column): install the passed vector by
  reference into `x$dense[[x$map[columns[k]]]]` - drop the `matrix(...)`
  coercion and the `column[] <-` copy. A factor source column decays to its
  code vector once (the engine already validated the values as existing
  codes), as today. This is the setPredictor-accept recovery.
- container x, rows given (partial): the merge is inherent for a fresh
  supplied vector (the new column differs from BOTH old and supplied at the
  rejected rows), so this stays O(n) for one column - same as the interim,
  no worse. The O(changed cells) goal is met only by the in-place
  re-install pattern (the caller's vector already carries the merged state,
  so install is O(spine)); document both.

## Column-subset, per-observation, and partial rollback

Preserved EXACTLY - reference-install touches only the data@x write-back
after the engine returns, never the engine transaction:

- The engine's per-observation session ([[sampler.hpp:912@f7efebf4]]
  updatePredictorPerObservation -> UpdateSessionImpl:1023,
  observationWouldRemainValid:1051, commitObservation:1072, finalize:1089)
  installs one column in random scan order and rolls back exactly the
  observations whose move would empty a leaf in any tree of any chain. The
  bridge ([[R_interface_bartcore.cpp:2296@f7efebf4]]) returns the `installed` LGLSXP mask
  unchanged; the R5 partial path returns it ([[bartcore.R:100@f7efebf4]]). Plan 3 keeps
  this return value byte-for-byte - callers (IRT-style Metropolis samplers)
  depend on it to know which observations moved so a rejected proposal rolls
  back only the affected rows.
- The R-side maintenance stays `installPredictorColumns(data@x, installed,
  column, x[installed])`: it starts from the OLD column and overwrites only
  the installed rows, so the non-installed rows keep the old values the
  engine kept (matching the engine codes). Semantics identical to the
  interim; only the copy discipline changes.
- The jointly path ([[facade.hpp:325@f7efebf4]] updatePredictorPerObservationJointly, one
  shared column across samplers, index-aligned) keeps its per-sampler
  write-back loop ([[updatePredictorPerObservationJointly.R:87@f7efebf4]]); each sampler
  reference-installs its own slice of the shared installed mask.

## Test-predictor (x.test) symmetry

Deliberate asymmetry, documented: the engine has NO per-column test entry
point - bartcore_setTestPredictor re-quantizes the WHOLE test matrix, so a
single-column x.test update must assemble the full matrix R-side to pass it
([[bartcore.R:353@f7efebf4]]). There is no reference-install win to take: x.test stays a
plain matrix, the column update stays a copy-modify (it is building the
engine's argument, not maintaining a live snapshot). No bench metric covers
it; leave it. (Containerizing x.test for a columnar per-column path is
plan-4/5 scope if it ever earns its keep.)

## extract() authority

No change beyond a comment/Rd clarification: extract(sampler, "predictors")
stays pure R over the container ([[generics.R:823@f7efebf4]]) and is already current for
every column, because the container IS the current raw store. The design's
plan-3 note "reroute extract to owned engine raw for mutable columns" is
DROPPED - there is no engine raw to reroute to, and the container is
authoritative. Rd gains one line: the returned codes reflect current values
after mutation.

## Commits

1. Reference-install helper + full-column and full-matrix rewire (R).
   installPredictorColumns replaces assignIntoPredictorSource; bartcore.R's
   updatePredictor column branch ([[generics.R:176@f7efebf4]]) and the full-matrix branch ([[generics.R:123@f7efebf4]])
   install by reference. Gate: equivalence 22/22 IDENTICAL vs
   equivalence-ac6ec2c.rds; full tinytest (test-sampler-predictors,
   test-sampler-setData, no snapshot regen); the revised bench criteria
   (header) - helper cheaper in isolation, reference semantics verified,
   nothing regresses. Abort: any draw, installed mask, or snapshot moves,
   or anything regresses. LANDED: 4283e3a - all gates held; the original
   recover-to-~1.0 criterion was corrected (see bench header): the
   bench's setPredictor scenario is matrix-backed and the accept gap is
   plan-1 structural, not the interim. Full-matrix setPredictor stayed a
   pointer swap (routing it through the helper would add cost; the
   plan's costs-nothing clause was not met).
2. Per-observation partial + jointly rewire (R). [[bartcore.R:94@f7efebf4]] and
   [[updatePredictorPerObservationJointly.R:87@f7efebf4]] install by reference; the
   installed-row-only merge is preserved. Gate: equivalence 22/22;
   test-sampler-setPredictorPerObservation,
   test-sampler-updatePredictorPerObservationJointly, test-mutate-then-
   serialize pass unregenerated (the installed mask is unchanged). Abort:
   the installed mask or a rolled-back row's value differs from the interim.
3. Comments + extract Rd + data@x contract (R + docs). Strike the "interim
   of design plans 1-2 / until plan 3" comments ([[bartcore.R:92@f7efebf4]], [[bartcore.R:172@f7efebf4]],
   [[updatePredictorPerObservationJointly.R:84@f7efebf4]], [[mixedMatrix.R:265@f7efebf4]]) for the
   reference-install final model; extract Rd + dbarts.Rd note data@x is the
   live current source; a small tinytest asserting extract == current after
   a column and a per-observation mutation, and that a shallow copy still
   diverges. Gate: checkRd/undoc/codoc clean; full tinytest. Abort: extract
   returns non-current or the shallow-copy divergence assertion fails.
4. Landing note (docs). The design landing note (plan-1/2 style) records:
   the mutable-flag KILL and its consequences, the data@x reconciliation
   (snapshot phrasing struck), extract left container-authoritative, no
   dbarts.h change / no stan4bart lockstep. The user-visible bits (the
   supported in-place re-install pattern; data@x reflects current values)
   are flagged for the 1.0-0 NEWS pass, not written here. Gate: R CMD
   check man; full tinytest. Abort: none (docs).

## Verification

- Full tinytest per commit from a preclean install (R side only; no header
  edit, so no --preclean strictly needed, but install fresh).
- Equivalence 22/22 IDENTICAL vs equivalence-ac6ec2c.rds at every commit -
  the frozen-default gate; any deviation is a defect, stop.
- The five mutation-focused files pass UNREGENERATED: test-sampler-
  predictors, test-sampler-setData, test-sampler-setPredictorPerObservation,
  test-sampler-updatePredictorPerObservationJointly, test-mutate-then-
  serialize. Reference-install is a copy-discipline change, not a value
  change - a forced regen anywhere means a value leaked.
- bench-sampler interleaved A/B on the quiet machine (>= 5 reps each build,
  medians) vs bench-sampler-32fc7c8.csv: setPredictor-accept within the
  6-8% noise band of baseline (ratio ~1.0); setPredictor-reject and every
  run-* metric unmoved. Single-run flags are re-run before belief.
- No tests/cpp change (no engine edit); no rchk change (no PROT edit).

## Open questions for VD

- RESOLVED (VD, 2026-07-13): the data@x contract. VD approved the
  reconciliation conditional on the model being mutation-COLLECTION (an R
  list of columns holding what was mapped, mutation swapping columns in
  and out by R reference) and NOT a maintained public view of the engine's
  quantized state - which is exactly the model (see the reconciliation
  section). The mutable-raw flag stays dead; extract stays the on-demand,
  constructed-when-called materializer. Direct data@x access remains
  supported; extract is a convenience/stability shim over the slot's
  varying vehicle, and un-exporting it pre-release remains cheap if it
  proves dead weight.

## Landing notes

Commit 1 = 4283e3a. installPredictorColumns replaces the copy-modify
interim on the full-column updatePredictor path: the accepted vector
installs into the container's per-column list by reference (one shared
coercion hoisted above the .Call), unmutated columns stay shared, and
only the O(p) spine allocates. Bench-attribution FINDING recorded here:
the original success criterion (recover setPredictor-accept-n1000-t75 to
~1.0 vs the 32fc7c8 baseline) was mis-attributed. Measured: the bench
scenario is matrix-backed (genFriedman), out of this plan's scope by
design; even container-backed, the interim was ~10 us of a ~260 us
update, so the ~50 us gap is structural to the plan-1 write-through
death (engine-side update cost, on both paths), not the interim removed
here. Correction landed df55a97. Neutral: equivalence 22/22 identical
draws, tinytest 2763 no regen.

Commit 2 = 0ff4842. The per-observation partial branch and the jointly
path route through installPredictorColumns, starting from the old
column and overwriting only the engine-installed rows, so non-installed
rows keep the values the engine kept. The installed-mask return value
verified byte-identical pre-vs-post on a forced partial rejection
(identical masks, identical kept and applied values, extract matching
data@x). assignIntoPredictorSource had no remaining callers and was
deleted. Neutral: equivalence 22/22 identical draws, tinytest 2763 no
regen.

Commit 3 = landed with this arc. Struck the "interim of design plans 1-2"
comments at the per-observation partial site (bartcore.R) and the
jointly path (updatePredictorPerObservationJointly.R) for the
reference-install final model - the full-column site already carried
that wording from commit 1. extract.dbartsSampler.Rd gained the data@x
contract line: each mutating call collects the accepted change into the
sampler's data object directly, nothing maintained per iteration, so the
returned codes always reflect current values. New tinytest coverage in
test-sampler-predictors.R: extract(sampler, "predictors") tracks
data@x's current materialization through both a full-column and a
per-observation partial mutation, and a shallow copy taken before either
mutation does not follow it (extends the existing full-column
shallow-copy divergence assertion rather than duplicating it). Gates: R
CMD INSTALL clean; checkRd/undoc/codoc clean on the touched Rd; full
tinytest 2771 (8 new) no regen; equivalence 22/22 identical vs
equivalence-ac6ec2c.rds.

Commit 4 = this commit (docs). This landing note, plus
docs/design/data-ownership.md: implementation-split item 3 marked
LANDED, and the second-round "OPEN AT PLAN 3" line amended in place to
record the resolution. Gate: full tinytest (confirmation; docs-only
commit, no other re-runs needed).

The mutable-raw flag is KILLED, never built. Reference-install makes the
container itself the current raw store, so the flag's sole rationale (a
CoW-free home for updated columns) is gone; building it now would only
re-add a full owned double column per mutable predictor - the ~400 MB
borrow this program set out to delete, re-entering by the side door.
data@x is the collected raw source per the approved contract (VD
confirmed 2026-07-13): the per-column values R mapped at creation, with
an accepted mutation swapping the affected column(s) in by reference at
the mutating call; it is never a representation of the engine's
quantized state, and the sampler never writes to it during sampling.
extract(sampler, "predictors") stays container-authoritative - pure R
over data@x, constructed when called, never cached, with no engine
reroute, because there is no engine raw left to reroute to. No dbarts.h
change, no stan4bart lockstep delta. NEWS is deferred to the 1.0-0 pass;
the user-visible bits are the supported in-place re-install pattern and
that extract reflects current values after mutation.

FOLLOW-UP RESOLVED (2026-07-14, phase-timer decomposition on both
revisions, env-gated diagnostic patch, never committed): the engine is
EXONERATED - phase for phase the tip's transaction is FASTER than
32fc7c8 (166.7 vs 171.7 us/update at n=1e3; 87 us faster at n=1e4,
because the write-through and oldValues memcpys the baseline paid were
O(n)). Every prior candidate cleared: validation ~0 on both, snapshot
0.16 vs 0.50 us, rewrite flat, bridge 0.08 us and byte-identical. The
whole gap is R-side maintenance the write-through era never performed:
per accepted update, ~30 us for the matrix-source O(n*p) subassign in
installPredictorColumns (container source: ~3 us) + ~14 us for the S4
data@x slot write-back, on top of ~44 us of R5 wrapper dispatch that
the BASELINE ALSO paid. Net new cost vs 32fc7c8: ~17 us/update
(container) to ~44 us (matrix) at n=1e3; at n=1e4 a container-backed
sampler already matches the baseline end to end (2799 vs 2733 us,
within noise) while matrix-backed is the worst cell (O(n*p) subassign
~195 us). The regression is a small-n fixed-overhead phenomenon plus
the matrix-source copy.

Recovery DECISION (VD, 2026-07-14): ACCEPT + DOCUMENT. The R-side
collection cost is the price of R-level consistency; a truly
performance-bound client drives the engine through the C ABI
(dbarts.h), which pays neither the R5 dispatch nor the collection and
supplies its own raw for replay (the documented plan-1 contract).
Documented in dbartsSampler-class.Rd (Mutation cost subsection:
frame input preferred for mutation-heavy loops; C interface for
microsecond-sensitive clients). The opt-out flag (b) is NOT built -
revisit only if a real workload demonstrates the need; transparent
lazy (c) stays rejected (direct data@x access is a supported contract
and would silently see stale values). Draw-neutrality was never at
risk (the write-back touches no codes, cuts, RNG, or fits).

