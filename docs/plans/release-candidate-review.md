# Release-candidate review program

Status: SPECCED 2026-08-17 (base 96ab54e0), amended per the blind
critique of the same date (verdict: execute with amendments; all
eleven applied - untracked record at
.claude/rc-review-critique-2026-08-17.md). The TODO entry
release-candidate-review is the charter; this file is the derived
slate and the program's record. Landing notes append at EOF.

## 1. Charter (restated)

An orchestrated review wave before declaring a release candidate for
serious human review, targeting the two defect families the per-slice
gate structure cannot see. Family 1, correctness beyond the gates'
reach: the slice batteries prove diffs neutral against agent-recorded
baselines, never that the baselines were right. Family 2,
accumulation: agent-written slices breed sediment, missed
factorization, surface inconsistency, dialect drift, and
over-complication. The charter's seed lists are a floor; this slate
was re-derived from two fresh censuses (2026-08-17, untracked working
notes at .claude/rc-review-census-gates-2026-08-17.md and
.claude/rc-review-census-accumulation-2026-08-17.md, both
citation-dense against 96ab54e0; the blind critique re-verified six
of six sampled load-bearing claims). Load-bearing ordering:
accumulation cleanups land before the coverage/mutation/coherence
passes; behavior-level oracles may run parallel to the code lanes,
subject to the freeze protocol below. RC declaration is VD-held.

Two standing facts every CI-touching item must respect:

- GitHub registers workflows from the DEFAULT branch only. Ten of the
  eleven workflow files exist only on bartcore, so they have never
  been registered: schedule, workflow_dispatch, everything. Any item
  that wants a CI leg (P2, P11's full-mode schedule, P14, P15, P17's
  mechanized half, equivalence-statistical-mode) is gated on a
  default-branch commit, which is VD-held (nothing is ever pushed to
  main by this program). Until then the fallback is local execution
  of the underlying scripts. main also has no benchmarks/ tree, so
  any shim landed there must be dispatch-only, schedule: stripped,
  or it installs a permanent red on the release branch. The API also
  still lists a fileless registered workflow (abi-contract, file
  deleted at 1f8a360e); enumeration tooling must not trip on it.
- Installed-package staleness: R CMD INSTALL without --preclean
  leaves a stale .so after header edits (facade.hpp virtuals
  bus-error). A build-stamp guard (stamp the build, assert the stamp
  moved) lands in wave 0 and every src-touching battery uses it.

## 2. What the censuses established (the derivation evidence)

Family 1 census, headline facts, each verified in-tree (counts as
corrected by the critique):

- Fourteen of the sixteen exact/balance math gates run single-tree
  (only heteroscedastic-exact at 20 and multinomial-exact at 50
  exceed one); the shipped default is 200 trees. The backfit loop -
  residual bookkeeping, totalFits accumulation, sigma conditional on
  the full ensemble - is gated by nothing at ensemble scale.
  gate-blindspot-audit.md reached this and it was never closed.
- Ten of the eleven CI workflow files have never executed anywhere
  (default-branch registration, above): rchk, valgrind, ALL SBC, the
  59-scenario statistical equivalence gate, revdep-smoke, and the
  full-mode/scheduled halves of exact-gates. equivalence.yaml's YAML
  itself is unvalidated.
- Three recorded defects are live and ungated: multinomial silently
  drops DART (chain.hpp hard-sets useDart = false while bart2
  threads it); Student-t pointwise log-likelihood evaluates the
  gaussian density - the packaged fit records family =
  fit$model@family, and resid.dist is orthogonal to family, so a t
  fit is indistinguishable from a gaussian one in the fit object and
  R/generics.R's dnorm branch runs (the fit carries neither
  resid.dist nor the estimated df); student/grouped + variance =
  constructs unrefused and untested (R/spec.R:369-377).
- The tinytest suite is ~2% math-vs-code (~80 of 4204 assertions),
  62% structural, 14% dbarts-vs-dbarts. 225 equality assertions pass
  when the compared field is NULL on both sides; the shared helper
  statesAgree() returns TRUE for two empty lists and carries 18
  assertions across 15 files. 405 assertions are
  expect_true(all(...)), which passes on an empty operand. 77 refusal
  assertions are pattern-free (50 of them in test-xbart-error.R).
  The helper's stated premise for collapsing to a boolean (tinytest
  ignores nested expectations) is FALSE - measured: expectations
  nested in a helper are collected and attributed to the call site -
  so per-field inlining is available.
- 12% of the suite can silently not run: test-capi.R (the whole
  shipped-ABI contract, 201 assertions) vanishes if R CMD SHLIB
  fails; test-simd.R is disabled by exactly the CPU-detection defect
  class it exists to catch; test-aft.R's survival-gated assertions
  never run (survival is not in Suggests); tests/tinytest.R passes
  with zero tests; the six sparse/mixed test files exit_file on
  Matrix.
- xbart is the least-verified surface: zero equivalence scenarios,
  no SBC, no exact gate, loss values never checked against a hand
  computation, weighted-vs-unweighted indistinguishable in its tests.
- Mutation testing is unrepeatable (the prior sweep's driver was
  never committed); no C/C++ coverage has ever been measured;
  vignettes are never built or executed in any CI leg; twelve
  benchmarks/R scripts belong to no workflow, and two of them are
  the only checks the shipped hazard and hurdle families have.

Family 2 census, headline facts:

- Three live cross-surface defects found by refusal-parity reading,
  all absences: bart() never attaches $weights, so ppd/loglik on a
  weighted bart() fit is silently unweighted (bart2 attaches them);
  bart() has no thinning guard where bart2 and rbart_vi have
  divergent copies; xbart(weights =, family = "probit") skips the
  binary weight refusal every other surface makes, because xbart
  bypasses the shared resolver.
- The dbartsMixedMatrix container parse is validated three times in
  the bridge with a shared predicate spelled differently, and has
  drifted: creation silently ignores malformed CSC metadata both
  mutation paths reject.
- Prose sediment survives the landed readability review and the
  dcc8262e comment sweep specifically where a comment's falsity is
  only visible from another file (R/bartcore.R asserts capabilities
  are unreachable that shipped the next day). Both prior passes read
  files in isolation; that is the nameable blind spot.
- inst/NEWS.Rd's unreleased 1.0-0 section is a 150-item slice log
  that narrates its own reversals and advertises a retired spelling.
  NEWS.Rd SHIPS (not in .Rbuildignore, not path-ignored), so NEWS
  edits are shipped-doc edits, not docs-only.
- The bridge carries ~13 hand-spelled chain ternaries, 8 forest-index
  copies, 10 present-tense "classic engine" comments, and duplicated
  guards across the language boundary (verbatim shared error strings,
  differing predicates for the same rule).
- Binary xbart's k axis is broken, not merely ambiguous: the default
  evaluates to a chi hyperprior object, kIsGrid comes out FALSE, the
  cross-validation axis collapses to ONE cell, and k is then sampled
  within it - the k cross-validation the function exists to perform
  does not happen, and the comment at R/xbart.R:311-312 asserts the
  opposite (Fork 1).

New review kinds this derivation adds to the seeds, both families:
cross-surface refusal parity (the sharpest findings are absences;
duplication-hunting cannot see them), gate-survival ranking (a
finding that survived a dedicated prior gate names a structural blind
spot), cross-language contract duplication, release-note
accumulation, arbiter-doc decay, comparison-to-a-vanished-baseline,
discharged-constraint audit, executable composition matrix,
mutation-sequence composition (distinct from the matrix: see wave 3),
refusal-path census, baseline-provenance audit, reduced-power-gate
audit, re-record discipline.

Census findings are tool-verified hypotheses; every slice below
re-verifies its own claims against the live tree before editing (the
standing discipline). Error-path resource leaks are a real class not
sliced here: they are covered by the valgrind leg (~900 expect_error
paths under memcheck) once its registration lands.

## 3. The slate

Slice letters (A-L, K2-K6) are the accumulation census's; P-numbers
are the gate census's; I-numbers are the re-verified interleavable
tickets (.claude/interleavables-verification-2026-08-17.md); FX are
fix slices this derivation adds. Sizes are dense non-comment lines.

### Wave 0 - rule, repairs, premise checks, live-defect fixes

Internal order matters and is part of the spec:

0a. The message-rule DRAFT (Fork 5's artifact, drafted first): a
    docs/design document fixing argument quoting, verb-per-refusal-
    kind, period policy, and the C-side caller-prefix policy. Wave
    0/1 slices write new messages against the DRAFT; VD signs the
    rule off before the wave-2 sweep executes, and the sweep absorbs
    any drafted-vs-final delta once. (Without a rule there is no
    "current house style" to write against - measured: six
    length-mismatch templates, four enum-rejection terms, four
    missing-argument verbs, three caller-prefix conventions.)
0b. Suite repairs and battery hardening, before the slices that lean
    on them: P3i statesAgree() repair (refuse zero-length and NULL
    fields, per-field attribution - inlining is available since the
    nesting premise measured false); P4 skip accounting and
    assertion floors (missing tinytest hard-fails under CI;
    test-capi.R compile failure fails under CI; test-simd.R exits
    when it would assert nothing; survival to Suggests or the dead
    assertions dropped; per-file floors for high-value files,
    sparse/mixed files included); the build-stamp guard (section 1).
0c. Premise checks, cheapest first - these test the program's
    founding claim before a dozen slices gate against the baselines:
    P1a ensemble-scale sum-invariance in tests/cpp (at m = 200,
    ensemble fit == from-scratch resum of per-tree fits; totalFits
    == sum of treeFits after every sweep); P10 baseline-provenance
    audit (classify every current baseline's re-record chain as
    oracle-proven / argument-adjudicated / expected-class-accepted;
    the third category is where a wrong baseline would have
    entered - reading only, no build). P17's WRITTEN half (a
    re-record must name its oracle in the MANIFEST row or commit
    body) also lands here, BEFORE the owed re-record; its CI
    mechanization waits on the default-branch gate.
0d. Live-defect fixes: I1 multinomial silent-drop refusals (five
    leaking arguments - dart, a DART tree.prior, split.probs,
    variance, and monotone, which is HALF-applied: proposal.probs
    reach the engine, directions do not, so draws move; SMALL, in
    flight at spec time). I2 hurdle variance doc fix (the
    per-observation sigma(x) promise is unreachable by design and
    lives in BOTH man/bart2.Rd and man/dbarts.Rd; TINY). I3
    getForestVariableCounts dimnames (rownames from
    colnames(data@x); the raw run()$varcount stays unnamed -
    separate, larger decision; SMALL). A bart() parity (attach
    $weights/$weights.test as bart2 does; thinning guard; zero-draw
    refusal; fixes a silent wrong answer; SMALL ~15). FX1-guard:
    record resid.dist on the packaged fit and REFUSE
    pointwiseLogLikelihood on a non-gaussian resid.dist, so a t fit
    stops silently returning gaussian values (R-only; the compute
    option is Fork 2 and its channel slice is post-K2). FX2
    student/grouped + variance = : Fork 6 decides refuse vs
    adjudicate; implemented per the fork's outcome (scheduled last
    in 0d to give the fork air).
0e. B xbart through the shared resolver (closes the binary-weight
    silent acceptance; the two doc-defended divergences survive by
    construction; MEDIUM ~60). Then B2 container-parse unification
    (share the PREDICATE - the three call sites keep their
    caller-specific messages and side outputs; creation tightens to
    the mutation paths' validation; MEDIUM ~70, src-touching).
    B2 runs AFTER P4 because its only real gate is the sparse/mixed
    tinytest files, all currently behind exit_file on Matrix: its
    gate is those files CONFIRMED RUN plus a malformed-metadata pin
    at creation. The bitwise trio is a NON-GATE for B2 (no harness
    constructs a dbartsMixedMatrix) and must not be recorded as
    evidence; a mixed-matrix scenario joins the owed equivalence
    re-record.
0f. P2 statistical legs, locally: run benchmarks/R/equivalence.R's
    statistical mode and the SBC suite as plain local Rscripts
    (that is all "the dormant workflows" can be until the
    default-branch gate opens); run the twelve orphaned
    benchmarks/R scripts once (hazard-reduction.R and
    hurdle-reduction.R first - the only checks those families
    have), then adjudicate each: gate / document as manual /
    delete. rchk and valgrind cannot run on this host (Linux
    containers) and stay gated on Fork 4. revdep-smoke by hand is
    subsumed by P12.

### Wave 1 - accumulation cleanups (Family 2 execution)

Recommended order within the wave: C, K3 (comment/prose, highest
gate-survival value) -> D, K2, K4, K5 (mechanical collapses) -> G, E,
K -> F (carries Fork 1) -> H, I, J (docs). New messages follow the
0a draft rule.

- C stale-prose sweep, cross-file-checked (4 R files + the NEWS
  rngSeed spot): the rule that catches this class - a comment
  asserting a capability's absence must be checked against the file
  that would provide it. SMALL. Note: the NEWS touch makes this a
  shipped-doc slice (gate class below), not byte-identical.
- K3 bridge/engine sediment sweep: ten "classic engine" comments,
  five comparison-to-a-vanished-baseline comments replaced by the
  invariants they stood in for, dbarts.h marker placement, the 2016
  rc/bounds.h TODO. SMALL, comment-only, byte-identity provable.
- D column-name resolution + pdbart prologue extraction. SMALL ~35,
  neutral.
- K2 run-result packaging compaction (one contiguous bridge region;
  ~170 lines -> ~95; makes the next result channel a one-line add -
  FX1's channel slice depends on this). MEDIUM ~75 removed, neutral,
  src-touching.
- K4 bridge indirection trim (the shared header itself is explicitly
  out of scope). SMALL ~50 removed, neutral.
- K5 getListElement shadow removal (the vendored rc_ copy adds the
  missing PROTECT) + bartcore_setModel dead parameter. SMALL,
  --preclean.
- G dead elements incl. n.burn[3] (a discharged-compatibility
  carve-out; narrows a public signature, so it is surface-bearing and
  lands pre-release by the utility rule). SMALL ~55.
- E withFixedSeed adoption at six sites (stops leaking the fixed seed
  into .Random.seed on the error path). SMALL ~45.
- K calibrationMapName K != 2 message (the user-visible sliver of the
  ledgered bcf-naming-generalization item, explicitly separable).
  SMALL ~4.
- F prior-ladder/collision-refusal factorization: refuseColliding
  takes the collision set as data so xbart's defended divergence
  survives; the binary-xbart k-axis repair is Fork 1's outcome.
  MEDIUM ~55.
- FX1-channel (post-K2, Fork 2 permitting): the per-draw df channel
  (the nbinom dispersion channel is the shape precedent) and the t
  marginal in pointwiseLogLikelihood. src-touching, --preclean,
  ASAN leg, dbarts.h addition.
- H NEWS 1.0-0 consolidation: rewrite the 150-item slice log as a
  delta from 0.9-31. Shipped-doc slice: fires full CI, R CMD check
  parses it; gate class below. It IS the release note.
- I core-generalization.md current-state extraction (the Family 2
  arbiter must answer "is this divergence intended" without a history
  read). Docs-only.
- J bcf-b-ridge exponent-rule promotion to
  docs/design/multiplier-combiner.md (the ledgered door-memo
  question, answered: promote). Docs-only.

### Wave 2 - message normalization

- L error-message sweep: VD-signed rule (0a draft -> Fork 5) applied
  repo-wide (~120 strings), updating every message-pinned test once.
  K6 (cross-language guard reconciliation: who owns each
  boundary-duplicated message, which predicates are defence-in-depth
  for the directly-reachable .Call) lands immediately after, written
  against the rule. Ordering rationale: wave-3 passes pin messages;
  pinning before the sweep would pin twice.

### Wave 3 - test-strength and measurement passes

Run after waves 0-2 so reshaping cannot invalidate them.

- P8 refusal-path census: every stop()/Rf_error site reached by a
  message-pinned assertion; unreached sites are the finding. Feeds
  P3ii-iv (expectFieldsEqual, expect_all, patterns for the 77 bare
  refusal assertions) as one combined pass, plus a lint gate against
  new instances of each vacuity shape.
- P16 xbart hardening: hand-computed loss values, 1-thread == n-thread
  pin, weights-change-the-answer pin, patterns for the 50 bare
  errors; the equivalence scenario rides the next re-record.
- P5 executable composition matrix: generate every (family x
  capability) cell feature-matrix.md marks S or ?, record
  construct/refuse/error, diff against the matrix, fail on
  disagreement. Doubles as the Rd-vs-behavior leg of coherence.
- P5b mutation-sequence composition (NOT derivable from the matrix,
  whose shape is one row per model x one column per capability):
  seeded from gate-blindspot-audit.md Half A's eight OPEN items -
  per-observation setPredictor x probit/logistic, installTrees warm
  start, setCutPoints, multi-chain posterior x {BCF, DART,
  linear/GP}, BCF x mutation axes - the between-sweep mutation
  surface that distinguishes this package. Includes the round-trip
  DEPTH question that audit left open: state round-trips must assert
  posterior CONTINUATION, not just serialization equality (P3i fixes
  vacuity, not depth).
- P6 repeatable package-scale mutation harness: committed driver
  (mutation list + the wave-0 build-stamp guard + battery runner),
  seeded from the recorded 16, extended to R/ (never mutated, no
  recompile cost). The durable asset of the program.
- P7 C/C++ branch-reach one-shot: gcov/llvm-cov over tests/cpp + the
  tinytest run, delivered as a never-hit-branch list (dead code ->
  Family 2; untested paths -> P6 targets). Explicitly not a coverage
  gate (repo-modernization closed that).
- P11 reduced-power-gate audit: re-run each gate's originating poison
  under quick mode; any quick gate that cannot kill its own poison is
  mis-tuned. The weekly mode=full schedule half is gated on the
  default-branch commit (section 1).

### Wave 4 - coherence and infrastructure guards

- Feature-matrix anchor resync over the FULL anchor namespace the
  matrix declares (RIB/MOD/CH/FAC/COM/MOV/SAM/CAPI, sampler.Rd,
  bart.Rd, and the R files), not just the R surfaces - waves 0-2
  rewrite the bridge in seven slices and man/ in three. Add the
  missing man/bart2.Rd alias (S14 created the file; the alias table
  predates it). Anchor passes run last, after content edits are
  final, and are verified by sampling for identical-delta offset
  signatures.
- P15 index/living-reference freshness gate (INDEX completeness +
  feature-matrix staleness ceiling); CI half gated on the
  default-branch commit.
- P13 .win drift guard (version literals vs DESCRIPTION; record the
  deliberately-absent macro set in config.h.win).
- P14 vignette execution leg: local run now; the CI leg is gated on
  the default-branch commit.
- equivalence-harness-statistical-mode (already ledgered): the
  cross-platform statistical compare the bitwise gate cannot give;
  its CI legs are gated on the default-branch commit.

### Wave 5 - the outer verifications

- Fresh adversarial whole-branch review (independent agents, refute
  posture, whole-tree scope; scheduled last so it reviews the settled
  surface).
- P12 consumer-bundle reality check: build/check bartCause,
  stan4bart, treatSens compat branches against the tip; Fork 3
  decides bairrtt; fix revdep-smoke to check compat branches rather
  than CRAN releases (checking CRAN releases against a breaking dev
  tree can only fail-as-designed or pass-vacuously).
- Then the RC call - VD's.

### Oracle lane - parallel, under the freeze protocol

FREEZE PROTOCOL: if any oracle fires against the tip, the cleanup
lane freezes until the defect's fix lands and the baselines
re-record; landed-slice neutrality evidence against a
known-wrong baseline is void and those slices re-gate. This is why
P1a and P10 run FIRST (wave 0c) - they are the two cheapest tests of
the program's premise and they precede the slices that gate on it.

- P1b conditional-exactness at scale (benchmarks/R/backfit-exact.R):
  freeze all but one tree, check its leaf draw against the closed
  form given the other m-1 as offset (test_model.cpp already codes
  the marginal independently). Unordered once P1a has reported.
- P9 Geweke marginal-conditional check, gaussian family at m = 5-10,
  using the sampler's own between-sweep setResponse and SBC's band
  machinery. Unordered once P1a has reported.
- SBC deepening for the five families with neither SBC nor a gated
  posterior oracle, on the sbc-family-tiers machinery.

### Riding items

- binary-kforest-k1-reachability: still decision-gated on acceptance
  evidence; the re-verification (2026-08-16) found four count-keyed
  bridge guards beyond the ticketed floor. An evidence memo may be
  commissioned during waves 1-3; the open/keep-refusing call is VD's.
- M0 component-contract doc stays deferred at discretion.
- The next baseline re-record carries: bart2 gaussian/probit/
  two-forest scenarios (bart2-argument-consolidation.md section 7
  preamble), an xbart scenario (P16), and a mixed-matrix scenario
  (B2) - and it happens under P17's oracle-naming rule, which is why
  that rule lands in wave 0c.

## 4. Forks for VD (batched; recommendations attached)

1. Binary xbart's k axis (surfaced by F): today the default
   evaluates to a chi hyperprior, the CV axis collapses to one cell,
   and k is SAMPLED inside it - the k cross-validation does not
   happen, and the in-file comment claims otherwise. Alternatives:
   (a) binary default becomes a genuine fixed-k grid (restores the
   function's purpose; mirrors the gaussian arm's semantics);
   (b) keep hyperprior-sampled k and drop k from the binary CV axes
   explicitly (documented, honest, smaller);
   (c) status quo + comment fix only (records the collapse as
   intended - hard to defend). Recommendation: (a), exact default
   grid values decided at slice time from the gaussian arm's shape.
2. Student-t pointwise log-likelihood, the compute half (the wave-0
   guard ships regardless): (a) leave it refused; (b) evaluate the
   gaussian density conditional on each draw's latent lambda_i (the
   sampler's working likelihood; needs the latents in the fit);
   (c) evaluate the t marginal integrating lambda (needs a per-draw
   df channel; the nbinom dispersion channel is the precedent).
   (b) and (c) give DIFFERENT loo/waic answers. Recommendation: (c)
   - the marginal is the observation-level likelihood loo/waic are
   defined on, and the channel cost is one scalar per draw landing
   post-K2 where channels are cheap.
3. bairrtt: fourth consumer, R-API-only, no compat branch, drifting.
   In the bundle (gets a compat branch) or declared out?
   Recommendation: out of the release bundle, with a one-line
   record; re-homing is its own decision.
4. The default-branch commit that registers CI legs (rchk, valgrind,
   SBC, statistical equivalence, revdep-smoke; also unlocks P11's
   schedule, P14's leg, P15's gate, P17's mechanization): the only
   viable shape is dispatch-only shims with schedule: stripped until
   bartcore merges (main lacks benchmarks/ entirely; schedules would
   run permanently red). This commit is VD's to make or delegate -
   the program never pushes to main. Until then: local statistical
   legs, no rchk/valgrind. Recommendation: make the commit early;
   rchk and valgrind are the two legs the charter names that nothing
   else can substitute for.
5. The message rule: 0a's draft goes to VD for sign-off; the wave-2
   sweep executes only the signed version. (Drafting is not the
   fork; the rule's content is.)
6. FX2, heteroscedastic student/grouped (variance = constructs
   unrefused and unadjudicated): (a) refuse now, door memo recorded
   - refusal-to-support is a compatible later move, support-now
   without adjudication risks shipping wrong draws; (b) adjudicate
   support now (its own verification arc: does the variance forest
   compose with the scale-mixture/grouped machinery correctly?).
   Recommendation: (a) refuse with the door named
   (heteroscedastic-robust and heteroscedastic-grouped regression
   are nameable capabilities; the door stays open).

### Resolutions (VD, 2026-08-17)

All six forks answered the day the plan landed:

1. Binary xbart k axis: MATCH GAUSSIAN - the binary default becomes a
   genuine fixed-k grid mirroring the gaussian arm's semantics
   (option a). Slice F implements; grid values derived from the
   gaussian arm's shape at slice time.
2. Student-t loglik: T MARGINAL (option c). The wave-0 guard ships
   first; the post-K2 channel slice adds the per-draw df channel
   (nbinom dispersion precedent) and evaluates the t marginal.
3. bairrtt: IN THE MAINTAINED BUNDLE. It is a candidate to merge
   into stan4bart, and until that is decided it stays compatible
   with the bartcore branch. P12's reality check covers four
   downstream packages: stan4bart, bartCause, treatSens, bairrtt.
4. Default-branch commit: DEFERRED TO THE COORDINATED MERGE. When
   the release candidate is ready and the four downstream packages
   are ready, all bartcore work merges into the mains simultaneously
   - workflow registration happens then, not before. Consequences:
   rchk and valgrind have NO leg for the duration of this program
   (recorded limitation; error-path leak coverage waits with them);
   P2 is local-only; P11's schedule half, P14's CI leg, P15's CI
   half, P17's mechanization, and equivalence-statistical-mode's CI
   legs all land at or after the merge. Their local halves proceed.
5. Message rule: FOLLOW BEST PRACTICES FROM HIGHLY REGARDED R
   PACKAGES, refined same day: VD dislikes tidyverse in general, so
   its guides carry NO authority weight - a tidyverse-originated
   practice is adopted only on its own scrutinized merits. Evidence
   precedence in the revised rule (docs/design/error-style.md): the
   measured practice of highly regarded base-style packages (base,
   stats, Matrix, survival, lme4, mgcv) wins where it clearly
   speaks; guide-originated adoptions each state their own merit
   case; otherwise the in-repo majority stands. The revised rule is
   the one slices write against and the wave-2 sweep executes; its
   delta from the draft is reported to VD before the sweep runs.
   SIGNED BY VD 2026-08-17: the wave-2 L sweep and K6 are unblocked;
   L queues behind the in-flight FX1-channel landing (file overlap).
6. FX2: REFUSE WITH A DOOR MEMO (option a), with a binding shape
   constraint: no interface friction for adding support
   post-release - the variance formal stays, the refusal is a
   validation error only, its message says the composition is not
   supported (never "cannot be"), and the door memo records that
   support arrives by adjudication plus refusal-removal with no new
   surface.

## 5. Gate classes for this program

- Behavior-affecting R slices: full battery (private-lib install,
  tinytest FAILURES == 0 with the wave-0 floors active, equivalence
  trio bitwise vs the canonical baselines WHERE THE TRIO EXERCISES
  THE TOUCHED CODE - a slice whose surface no harness constructs
  must say so and name its real gate instead, per B2), air --check,
  lintr on touched files, R CMD check from a clean-staged tarball,
  plus a discrimination proof: new refusals/pins shown failing on
  the base build.
- src-touching slices add --preclean under the build-stamp guard,
  tests/cpp from make clean, and the ASAN local leg + CI sanitizer
  watch when new code becomes reachable.
- Neutral cleanups: the bitwise trio is mandatory; comment-only
  changes prove built-package byte-identity.
- Shipped documentation (man/**, inst/NEWS.Rd): full battery minus
  equivalence, plus NEWS parse (the from-package_NEWS_Rd invocation)
  and R CMD check codoc. NEWS.Rd ships and fires full CI; it is
  never "docs-only".
- Non-shipped docs (docs/**, TODO): byte-identity of the built
  package to the last gated code commit (these fire no CI).
- Workflow/CI slices: the workflow's own green run is the gate,
  which presupposes registration (section 1).
- Oracle additions gate themselves by construction (a new oracle must
  pass on the tip and kill a planted poison).

## 6. Deliberately not proposed

- The bcf-naming-generalization rename itself: ledgered separately,
  engine-sized, priced there (only its user-visible K message, slice
  K, is taken).
- bart()/rbart_vi result-packaging merge: the array shapes make it a
  rewrite, not a factorization.
- A comment-style rule (## vs #) invented only to normalize toward:
  fails the enabling-value gate.
- A C/C++ coverage gate or badge: closed by repo-modernization; P7 is
  a one-shot measurement.
- Naming the raw run()$varcount channel (I3's larger sibling):
  separate decision, touches the C bridge and FB14's pins.
- New standalone slices for serialization format, platform-float
  divergence, long-chain memory growth, RNG-across-threads: probed
  by the critique, no plan-sized counterfactual found beyond what
  P5b covers; error-path leaks ride the valgrind leg once
  registered.

## Landing notes

### Wave 0c - ensemble-scale gaussian SBC, the premise oracle (2026-08-17)

The provenance audit's named highest-value action ran at the tip:
Rscript benchmarks/R/sbc.R gaussian 200 200 30 (R = 200
replications, L = 200 retained draws, thin = 30), against the
a39da5d9 build (d60057b1 is records-only, binary identical), with
SBC_FAIL_ON_FLAG=1. ALL SEVEN FUNCTIONALS PASS: avg.f chisq p 0.294
/ ks 0.442, f.star1-5 chisq p 0.255-0.849 / ks 0.051-0.688, sigma
chisq 0.962 / ks 0.998; every ecdf-diff (0.026-0.094) inside the
0.132 simultaneous band; the sigma side-check P(sigma < sigest) =
0.8998 vs target 0.90. Exit 0. Consequences: the freeze protocol
does not trigger; the cleanup lane is clear; the audit's residue R1
(the fused roll+suffstat pass, oracle-checked only at single-tree
scale since c8f661a) is RETIRED at this tip, and R2's cutover-window
re-anchor is refreshed - the un-retired remainder is unchanged
(gp/wtgp/quants/missing/sparse/setdata have no varying arm; the
xbart pin still has no oracle, P16's job). Log preserved untracked
at .claude/rc-review-sbc-gaussian-ensemble-2026-08-17.log.

### FX1-channel - per-draw Student-t df and the t marginal (c3af16a1, 2026-08-18)

Fork 2 executed. The engine's df dynamics were verified FIRST: nu
is sampled per sweep in grid mode (TResponse::refreshLatents draws
every lambda then nu from the two lambda statistics; fixed mode
holds the creation value), so the plan's per-draw shape is right.
The channel mirrors the nbinom dispersion channel at all four
layers - Results::residualDf with the storeSample write gated on
carriesResidualDf (a pure read of state the sweep settled, after
refreshLatents, so no RNG is consumed), the per-chain slice, the
SamplerShape bit, and the bridge's installChannel("resid.df") -
K2's one-line-add promise held. dbarts.h gains double* residualDf
appended above the documented 1.0-0 field boundary per that
boundary's own rule (no version constant moves; offset
static_assert added, size assert 10 -> 11; the signature token is
layout-invisible and unmoved); capi contract 232 -> 237 with a
capi_run_residual_df shim. pointwiseLogLikelihood's student branch
scores the marginal location-scale t (dt((y - ev)/sd, nu) -
log(sd), sd = sigma/sqrt(w), df recycled on sigma's axis; a
student fit predating the channel refuses by name). sampleFromPPD
stays refusing; DOOR: with $resid.df the t draw is sd * rt(n, df)
in the gaussian arm's shape - one branch, no new state. NEWS: the
student item extended AND the refusal item corrected in place (it
claimed loglik refuses student, which this slice makes false).
Budget overran (~194 dense vs ~105 cap), flagged by the
implementer and accepted on review: the overage is the capi C shim
the local idiom requires plus mandated-oracle assertions, each
mapped to a mandate clause. Battery twice (implementer, then
independent gate-runner): install --preclean; tests/cpp 69 checks
incl. 3 student state round-trips; ASAN/UBSAN zero diagnostics;
tinytest 5852/0 with capi 237 confirmed run; trio bitwise 37/37
12/12 10/10 with the student scenario among the 37; air/lintr
clean; NEWS 268 parses; R CMD check Status OK; discrimination -
base has no $resid.df and refuses loglik by name, slice matches an
independently written dt() computation exactly, ppd refuses on
both. Patch-id-verified rebase over the P1b landing. Prior-slice
CI: H six-green, which also covers J's concurrency-cancelled
exact-gates/sanitizers by containment (recorded here as promised);
P1b's benchmarks push went exact-gates green.

### P1b - backfit conditional-exactness oracle; PASS (ded60c1d, 2026-08-18)

benchmarks/R/backfit-exact.R, 345 lines, 11.4s full / 3s quick.
Design: nothing frozen by hand - the sweep's own ordering supplies
the conditioning (entering tree j the running residual holds trees
1..j-1 at this sweep and j+1..m at the last, both carried by
consecutive keepTrees records), so EVERY tree is checked, not one.
Each recorded leaf value is standardized by the re-derived
ConstantGaussianLeaf posterior (precision (k sqrt(m)/node.scale)^2
+ sum w_i/sigma^2) into the exact standard normal the engine
consumed; iid across leaves/trees/sweeps, so moment and KS tests
are exact, with the prior-informed small-leaf stratum reported
apart. VERDICT: PASS on every functional, both arms - A at
ensemble scale (200 trees x 1000 obs x 400 sweeps, 193,141 leaf
draws; routing counts exact, fit-sum 8.97e-14) and B under case
weights + offset (all |z| <= 1.51); a borderline arm-A KS p =
0.047 dissolved under four extra seeds pooled to 568,647 draws
(p = 0.66, kurtosis 3.00001, per-tree-position chi-square p 0.64).
Teeth measured by poisoning the REFERENCE: a stale residual reads
z(var) +434, a dropped sqrt(m) z(mean) -9.3 pooled / -19.3 on
prior-informed leaves. Not seen by design: the tree-structure MH
step (the exact-posterior gates' job), the sigma/k draws, non-
gaussian families, non-constant leaves. Freeze protocol DORMANT.
Gates: fresh-Rscript clean run; air clean; lintr clean; diff is
the one new file (benchmarks push - fires exact-gates only,
watched). Run log preserved at the session scratchpad
p1bimpl-run.log.

### H - NEWS 1.0-0 consolidation; WAVE 1 CLOSES (7623af19, 2026-08-17)

The unreleased section now reads as the delta from 0.9-31: 155
items -> 153 (149 kept as already-final surface descriptions, 6
trimmed of branch-time narration or stale spellings, 2
supersession chains merged to final state - rngSeed -> seed and
the retired variance companion knobs -> varianceForest(), 3
dropped as intra-branch churn by their own admission, 1 added: the
deferred B2 item, dbartsMixedMatrix reference metadata validated
at creation to the mutation standard, wording verified against
requireCscReferenceMeta). The accreted duplicate subsection blocks
merge to one UPGRADING / NEW FEATURES / C API / BUG FIXES set. The
rngSeed rename disclosure stays IN FULL as the correctly-worded
final-state record of a real 0.9-31 breaking change (the
consumer-gotcha item); the census's actual complaint - rngSeed
advertised as the current spelling - is dead. Item-by-item
classification preserved untracked at
.claude/../scratchpad/himpl-classification.md (session scratch);
prior-release sections verified byte-untouched by the orchestrator;
every kept "gains" item verified to name a live surface. Gates
twice (implementer, then independent gate-runner on the
patch-id-verified rebase): parse 268 entries; rngSeed confined to
the disclosure item plus the historical 0.9-19 section; R CMD
check Status OK from a clean-staged tarball; diff touches only
inst/NEWS.Rd.

WAVE 1 COMPLETE: C, K3, D, K2, K4, K5, G, E, K, F, I, J, H landed
plus the VD-directed CI hang-guard slice; every code commit
six-green (H's legs in flight at landing, watched). Open residue
riding forward: grouped-mixing ratio flag (oracle lane),
mixedMatrix predictorSourceColumn naming residual (wave-3 refusal
census), the wave-2 L sweep gated on VD signing the revised
message rule, and FX1-channel now UNBLOCKED (Fork 2 resolved, K2
landed - the per-draw df channel is a one-line add).

### J - exponent-rule promotion (76be5e5b, 2026-08-17)

The ledgered door-memo follow-up discharged: the b-move Jacobian
exponent rule's derivation and prototype evidence (the operational
c-form substitution, the GIG match giving the general p = (k - d)/2,
and the proto-b.R run rejecting the naive Jacobian at KS 1.6e-21)
move from docs/plans/bcf-b-ridge.md sections 2.2/2.3/5a into
docs/design/multiplier-combiner.md, "The exponent rule", placed
after "The ASIS ridge". The plan keeps its BCF-specific invariance
argument, edge cases, and routing narrative, with MOVED pointer
stubs at the three sections (the composition-mixing-probe stub
convention). combiner.hpp's ASIS-ridge comment now cites the
design doc for BOTH facts, and test_sampler.cpp's shipped comment
sheds its docs/plans reference for the same citation. Repo-wide
citation sweep verified by orchestrator sampling per the anchor
discipline: mixed deltas (-34 between the moved sections, -57
below them, both spot-verified line-exact), in-move cites
redirected to the design doc with origin notes, and the one
seemingly-off recompute traced to a pre-existing loose range whose
semantics the recompute preserved. Gates: R CMD INSTALL --preclean
clean (comment-only src edit compiles); tests/cpp from make clean
all pass; line counts reconcile (bcf-b-ridge 530 -> 473,
multiplier-combiner 604 -> 696, glue accounted). TODO's follow-up
sentence resolved in this records commit. Prior-slice CI: F
six-green.

### I - core-generalization.md current-state extraction (a55be8b5, 2026-08-17)

The arbiter doc's 533-line Phases DONE-log and the pre-rewrite
Starting-point constraints section move below the architecture, so
the body reads as current state: Current architecture -> Goals ->
Architecture -> Extensions -> Performance -> Validation -> Risks,
then the quarantined history. Destination is a same-file appendix
with the Phases heading renamed Landing notes - the convention
every sibling multi-phase design doc uses (forest-combiner, bcf,
multiplier-combiner, sparse-columns, grow-from-root) and
docs/README.md defines - not a docs/plans/ file, which no sibling
precedent supports. The move is a verified pure permutation: the
sorted-multiset diff of the file is exactly the one heading rename
(orchestrator re-ran it independently), 802 lines before and
after. All 56 repo-wide citations of the file resolve unchanged -
none anchor by line or literal section heading; the phase-number
citations and the one exact-quote "Wave 2 models" citation are
preserved verbatim inside the moved block, and
state-continuation.md's pre-existing "landing notes" phrasing is
now literally true. Docs-only: fires no CI; the built package is
byte-identical to the gated F tree. Prior-slice CI: F's six legs
in flight at landing, watched.

### F - prior-ladder factorization and the Fork 1 binary-k repair (621da478, 2026-08-17)

refuseColliding promoted to R/model.R taking (matchedCall,
objectName, shorthands) - the collision set stays caller data, so
bart2's tree.prior collides with power/base while xbart's
doc-defended grid-axes divergence survives untouched; xbart's
drifted "the tree prior" wording is gone. The dart/cgm ladder is
one resolveDartShorthand with buildDart/buildCgm constructor
closures, because the callers build in different currencies (bart2
an unevaluated call forwarded to dbarts(), xbart a live object);
both formerly-duplicated literals now have one source. A_class.R's
third phrasing stays: it is a setValidity object invariant, not an
argument-shorthand collision. FORK 1 EXECUTED: xbart's grid
default k is the fixed value 2 for every family - .kDefault
deleted (it reached control out of the caller's frame and
duplicated resolveNodeHyperprior, whose chi(1.5, 2) binary default
remains the bart family's single source, dissolving
chi-hyperprior-df.md's lockstep concern); the once-false fixed-k
comment is now true; man/xbart.Rd rewritten; one NEWS item (270
parse) also records the side effect that drop = FALSE binary
results carry the length-one k dimension gaussian always had.
Orchestrator review caught the new test using a length-3 n.burn
against the just-narrowed surface; amended before gating. Battery
twice (implementer, then independent gate-runner on the
patch-id-verified rebase): tinytest 5833/0; trio bitwise 37/37
12/12 10/10 (the no-bart-family-contamination gate); air clean;
lintr clean on all six files; R CMD check Status OK from a
clean-staged tarball; discrimination three ways - default binary
xbart carries dbartsChiHyperprior on base and fixed k = 2 on the
slice, bart2's binary default is chi on BOTH (the non-change
gate), and the slice's test file fails 3 on base / 0 on slice.
Prior-slice CI: K six-green.

### K - calibration-map refusal above K = 2 (ceffd276, 2026-08-17)

The bcf-naming-generalization ledger item's one user-visible leg,
landed alone as scoped. calibrationMapName (4 lines + comment)
deleted; its ternary hoisted into bartcore_setCalibration's single
Rf_error as a three-way on the coupling: softmax and K = 2 keep
byte-identical text, K >= 3 now names the multi-forest calibration
map - the same noun the R layer's refuseBCFMutation refusal uses,
so the two routes agree while staying independently worded. The
surviving "owns both halves of its calibration" clause is
K-independent by the engine's own vocabulary (chain.hpp: BCF and
multinomial "derive both halves from their own maps"; softmax
carried the clause at K >= 3 all along). Three pins in
test-calibration-midchain.R reach the engine message past the R
guard through the low-level handle: two-forest text at K = 2 via
bcf, multi-forest presence and two-forest absence on a
three-forest sampler. No NEWS: the wrong text and its fix are both
inside the unreleased cycle. Battery twice (implementer, then
independent gate-runner): tests/cpp clean from make clean;
tinytest 5830/0; trio bitwise 37/37 12/12 10/10; air clean; lintr
clean; R CMD check Status OK from a clean-staged tarball (standard
invocation, manual built); discrimination - base fails exactly the
two K >= 3 assertions, slice 125/125. TODO trimmed to the rename
remainder in this records commit. Prior-slice CI: E six-green.

### E - withFixedSeed adoption (1fd22c6f, 2026-08-17)

The six open-coded .Random.seed save/restore sites route through
the on.exit-based withFixedSeed (xbart's chunk seeds - the one
already-correct copy - plus bart2's hurdle component seeds and four
rbart sites: the per-thread seed draw, the single-threaded chain
loop via a local fitChains closure whose returned list reassigns
the caller's, rbart_vi_fit_bartcore's whole
sampler-creation-through-run region as a block passed by promise
through a local seed-conditional delegate so its bindings land in
the function frame, and packageRbartResults' recorded seed). Two
deliberate behavior changes, both the leak the census named: an
error inside the fixed-seed window now restores the caller's
.Random.seed instead of leaving the fixed state, and a call that
CREATED .Random.seed removes it on exit (the helper's documented
leave-no-trace contract; the old open-coded sites skipped restore
when no seed pre-existed). Success path draw-identical. One NEWS
item. Discriminating test at the single-threaded chain loop: a
custom prior function that stops routes rbart_vi through the R
callback fit path so the error surfaces inside the window;
.Random.seed asserted equal to its pre-call state. Battery twice
(implementer, then independent gate-runner on the
patch-id-verified rebase onto the CI-guard tip): tinytest 5827/0;
trio bitwise 37/37 12/12 10/10 (the rbart_vi scenarios are the
leak detector); air clean; lintr clean; NEWS parses at 269; R CMD
check OK from a clean-staged tarball; discrimination both ways -
base fails the restore assertion (1/6) and leaves a created
.Random.seed behind, slice passes 6/6 and removes it. Prior-slice
CI: ci-hang-guards six-green, pkgdown inside its normal window.

### CI hang guards - job timeouts everywhere (6000efb7, 2026-08-17)

VD-directed after two same-day hangs. Diagnosis from the cancelled
attempts' logs: pkgdown run 32086436327 (2h20m) and lint run
32096324117 (~40m) both stalled inside r-lib/actions/setup-r's
"Updating system package data" step - apt-get update announced
Get:5 noble-security InRelease from the azure.archive.ubuntu.com
mirror and never completed the fetch; a same-day non-hung pkgdown
run hit the identical mirror throttling on package installs and
took 76 minutes to self-resolve, corroborating the class (mirror
stall, not an R download). Guard: job-level timeout-minutes on
every job of the eight workflow files that lacked one (12
insertions; revdep-smoke, sbc, valgrind already carried reasoned
timeouts), sized ~3x normal wall time, minimum 15: R-CMD-check 80 /
NEON leg 15, cpp-tests 15, exact-gates 25, lint+format 15,
sanitizers 40, pkgdown 90 (sized above the observed 76m benign
throttled run, not 3x the 16m reference), and generous bounds for
the not-yet-registered equivalence (120/90/90) and rchk (120) so
they carry the guard when the default-branch merge registers them.
Push events read the pushed ref's workflow file, so the guard is
live on bartcore now. Gates: all eleven files parse
(yaml.safe_load); the landing push itself fires all six live
workflows as the in-vivo check. Patch-id-verified rebase from base
75d82925 onto eb8c292d. Prior-slice CI: G six-green, no hangs.

### G - dead elements (99f98043, 2026-08-17)

xbart's n.burn narrowed to two non-negative integers (the third was
never read; formal default, rep_len, the stale prose comment, the
man/xbart.Rd compat sentence and its example all updated). Posture
toward a longer vector is silent truncation, matched to the
measured sibling posture - n.test is unconditionally reduced to
n.test[1L] regardless of supplied length - and pinned by a surface
assertion that a length-3 call still runs and equals the length-2
call, plus a default-formals pin that discriminates against base.
One NEWS item. convertSamplesFromBartsToDbarts deleted with its
only exercises (its inverse keeps its many callers); validateObject
inlined at its sole caller; xbartRunChunk's family == "auto" branch
dropped after live verification that xbart forces the family before
the model is built. F-J5 NOT CUT: its census premise is FALSE by
live instrumentation - a plain dbarts(formula, data) call, control
omitted, reaches validateArgumentsInEnvironment's func-reading
branch because match.call drops unsupplied arguments; the formal
stays. Battery twice (implementer, then independent gate-runner
with base-vs-slice probes): tinytest 5825/0; trio bitwise 37/37
12/12 10/10; air clean; lintr clean on all six touched files; NEWS
parses at 268 entries; R CMD check OK from a clean-staged tarball;
base runs length-3 n.burn identically, slice lacks the deleted
helper. Prior-slice CI: K5 six-green - its lint job hung at
setup-r past 40 minutes and was cancelled and re-run to green, the
second same-day r-lib-action hang after K4's pkgdown; the
timeout-minutes guard slice was queued in response.

### K5 - vendored getListElement, setModel arity (71691eb7, 2026-08-17)

The bridge's local getListElement shadowed the in-scope
rc_getListElement and read the names attribute unprotected across
the lookup; the vendored copy is contract-identical (R_NilValue on
a missing names attribute or an absent name, first match wins) and
holds the PROTECT, so the shadow is deleted and its 100 call sites
routed to rc_getListElement (zero stragglers by grep; every other
hunk in the bridge TU is line-rewrap from the longer name).
bartcore_setModel never read its control argument: dropped from the
entry point, the header declaration, the call-table arity (4 -> 3),
both R callers (dbarts.R's setModel method including its
error-rewriting quote(), bartcore.R's wrapper), xbart's call site,
and three tinytest files; xbart's cellControl stays live at sampler
creation. Battery (independent gate-runner, own lib): --preclean
install fresh by Built stamp; tests/cpp from make clean; tinytest
5825/0; trio bitwise 37/37 12/12 10/10; air clean; lintr clean on
all six touched files; R CMD check OK from a clean-staged tarball.
Prior-slice CI: K4 six-green - its pkgdown job hung past 2h and was
cancelled and re-run to green (run 32086436327); K2's cancelled
R-CMD-check is covered by containment (K4's green check builds the
superset tree).

### K4 - bridge indirection trim (1e229a77, 2026-08-17)

Deletions and inlinings, net -29 (-54 pure code; 25 of the
additions are doc comments relocated out of the shared header):
refuseRequantizeWithoutSource inlined into its three callers with
its constraint folded into refuseMutationOnView's comment; the
three refusal names consolidated onto two predicates
(refusePredictorMutation on acceptsNewRawPredictors,
refuseMutationOnView on hasRequantizeSource);
rawTrainingColumn/rawParsedTestColumn folded into rawViewColumn and
the one-field DataHandle dropped for a raw ColumnStore*;
createBCFHolder and installForests moved out of the shared header
after grep-verifying zero C_interface.cpp consumers (installForests
needed a forward declaration in the bridge TU - the one
non-pure-deletion, recorded); rc/bounds.h's dead clang-11 disjunct
replaced by the honest C++20 predicate with the comment restated
(CXX20 in both Makevars, so the clause could never decide); K2's
single-caller installResult folded into installChannel. Battery:
tests/cpp clean and under ASAN/UBSAN (zero diagnostics); tinytest
5825/0 with test-capi.R's 232 contract assertions confirmed run;
trio bitwise 37/37 12/12 10/10; R CMD check OK from a built
tarball.

### K2 - run-result packaging compaction (447b8c81, 2026-08-17)

The bridge's packaging region (live :4188-4462, 161 lines of channel
assembly) becomes 87: thirteen hand-rolled chain-ternary allocations
collapse onto one allocChannel(type, {leading dims}) lambda (a
single-chain per-draw scalar stays a bare dim-less vector - the
byte layout every consumer relies on - by construction) with scalar
and widening wrappers beside it; the four recomputed prefix sums
and the trailing 19-line names block become a running slot counter
plus installChannel(name, value), names attaching as slots are
claimed - which also retires the old free-the-counts-before-the-
names-block OOM-longjmp ordering dance; seven verbatim forest-index
reads become forestIndexFrom (the census said eight - the getTrees
variant keeps its own out-of-range message deliberately). Net
-79 lines. The enabling value, demonstrated: a new channel is now
one installChannel line plus one term in the slot sum, versus the
four-slot prefix-sum audit plus a names-block edit (the post-K2
Student-t df channel rides this). Neutrality, the strongest stack
yet: 14/14 full-result identical() probes covering every channel
(gaussian 1/2-chain, nulls, chi-k, DART varprobs, variance with
and without test, ordinal cutpoints, nbinom dispersion, BCF
forestFits/glue/widened varcount, multinomial train/test, grouped
tau/ranef, zero-sample NULL); trio bitwise 37/37 12/12 10/10;
tests/cpp clean AND under ASAN/UBSAN; the instrumented R-package
ASAN leg ran all 14 probes clean (ULP-level value drift vs the
plain build is instrumented-codegen artifact; names/dims/types
identical); tinytest 5825/0; R CMD check --as-cran OK. K4 note
recorded: installResult now has exactly one caller.

### D - column-resolution and pdbart-prologue extractions (3158d96e, 2026-08-17)

Two pure extractions, net -8 lines. resolveColumnIndex(source,
column, what) collapses the three near-identical resolution blocks
in R/bartcore.R (setPredictor/setCutPoints/setTestPredictor), the
per-site not-found messages preserved via the what parameter; the
pdbart/pd2bart prologue collapse renames the existing
pdbart.getSampler to pdbart.prologue and folds the duplicated
comment into it. One correctly-refused non-extraction, recorded:
match.call()/parent.frame() are frame-sensitive and stay inline at
each call site - moving them into the helper would capture the
wrong frame. Neutrality: a per-site probe table with (ok, value,
message) triples identical base-vs-slice for every success AND
error path (bad name, out-of-range index, unnamed source, bad
x.train type); tinytest 5825/0; trio bitwise 37/37 12/12 10/10;
air clean; lintr - partialDependence.R zero, bartcore.R's seven
object_usage warnings proven pre-existing and identical on base;
R CMD check --as-cran OK from a built tarball.

### K3 - bridge/engine sediment sweep (28ead6e2, 2026-08-17)

Comment-only, 54/54 lines across six files, all five findings
re-derived live. The ten present-tense "classic engine" mentions in
the bridge are gone (grep returns zero): the two describing the
real getTrees data.frame format reworded keeping the substance, the
two-engine invariant dropped outright, and the undefined
"reference-engine" claim replaced by refuseBinaryWeightChange's
actual predicate. The five comparison-to-a-vanished-baseline
comments (combiner.hpp x4, data.hpp) are now DERIVED INVARIANTS a
reader can check against the adjacent code - e.g. "yields exactly
today's pointer" became "off an offset this returns forest k's own
totalFits buffer, so no train-side reader sees a copy and raw_
stays empty; the offset costs an offset-free run neither an
allocation nor an addition". bartcore.hpp's header prose refreshed
to the shipped family/leaf roster (verified against
ResponseFamily/LeafModelKind); dbarts.h's 1.0-0 field-boundary
marker moved below the fields it governs with its apology deleted,
declarations byte-untouched; the 2016 rc/bounds.h TODO deleted -
the tree-wide marker grep over src/ and inst/include/ now returns
ZERO. Two candidates adjudicated OUT with reasons that sharpen the
kind: model.hpp's "classic form" is textbook algebra, and nine
surviving "today" statements have their right-hand side IN the
tree (checkable current-set claims, not vanished baselines).
Neutrality: clang -E -P identity on six TUs base-vs-slice - the
full engine via bartcore.hpp, the bridge, combiner, data, the
shipped header as -x c (the LinkingTo view), bounds - all
byte-identical, zero __LINE__ shifts. Battery: tests/cpp clean;
tinytest 5825/0; trio bitwise 37/37 12/12 10/10; R CMD check OK
from a built tarball (raw-directory check trips on Authors@R -
staging artifact, build-then-check is the recorded procedure).

### C - the cross-file-checked stale-prose sweep opens wave 1 (5f2ced48, 2026-08-17)

Five sediment findings, each verified still-stale against the live
tree then corrected as prose only: R/bartcore.R's comment asserting
the logistic family and categorical ingestion are unreachable
capabilities (both shipped - bart2's public family argument, the
data constructor's factors = "categorical" default - the comment
now states the shipped reality); its "future multinomial creation
route" (shipped); rbart.R's "we're are" typo; the NEWS rngSeed
advertisement corrected to seed (that one line only - the wholesale
consolidation is H's); utility.R's "now that the rename lands"
tense. The enforced rule, stated for the record: a comment
asserting a capability's ABSENCE must be checked against the file
that would PROVIDE it - the class both the landed readability
review and the dcc8262e sweep missed by reading files in isolation.
No finding needed a code change. Neutrality: parse-tree identity
(keep.source = FALSE) on all three touched R files, proven by the
implementer and INDEPENDENTLY re-derived by the orchestrator; full
battery green (tinytest 5825/0, trio bitwise 37/37 12/12 10/10,
air/lintr clean, NEWS parses at 267, R CMD check --as-cran OK).
13 insertions / 12 deletions. Landing gated on B2's six-green
(sanitizers included), which arrived first.

### B2 - the shared container-metadata check; WAVE 0 CLOSES (0a1d56fe, 2026-08-17)

One file-local requireCscReferenceMeta helper in the bridge
validates a mixed container's sparseReference/sparseCategoryCount
(integer, one entry per CSC column) and borrows the pointers; all
three funnels call it - parseTestContainer, parseMutationSource,
parseData - each keeping its own message and side outputs (the
critique's amendment honored: the shared part is the predicate,
never one function). Creation, which silently IGNORED all six
malformation classes both mutation paths refused (two fields x
absent/wrong-type/wrong-length), now refuses at arrival: measured
on base, a malformed all-ordinal-CSC container CONSTRUCTED AND RAN
(finite, statistically fine - ordinal columns never consume the
metadata) while the identical object was refused at every mutation
entrance; a categorical-CSC malformed container was already caught
by a downstream presence check, now dead and deleted (its message
becomes the shared one; nothing pinned the old string). Gates
twice, src-slice battery: tests/cpp clean AND under ASAN/UBSAN
(zero diagnostics; the implementer additionally ran an
INSTRUMENTED R PACKAGE build - probes clean); tinytest 5825/0 with
all six sparse/mixed files confirmed run at exact per-file counts
(69/15/44/16/20/137); trio bitwise and recorded as a NON-GATE for
the parse (no harness constructs a dbartsMixedMatrix); R CMD check
--as-cran OK; discrimination 3/3 corruption classes (base silent
construct-and-run vs slice refusal), well-formed fixed-seed sparse
fit bitwise across builds. RESIDUAL recorded, not fixed: deleting
sparseReference and then calling setPredictor errors through a
pre-existing R-side short-circuit ("argument is of length zero",
R/mixedMatrix.R predictorSourceColumn) instead of the container
message - one malformation class, one path, naming-only; wave-3
refusal-census material. NEWS deliberately deferred to slice H's
1.0-0 consolidation (this item rides that rewrite).

WAVE 0 IS COMPLETE. Landed, in order: I1 multinomial refusals, I2
hurdle doc fix, I3 varcount dimnames, the SBC ensemble premise
oracle (PASS, freeze protocol dormant), P10 provenance audit +
MANIFEST corrections + P17 rule, suite repairs (statesAgree
per-field, skip accounting, build-freshness guard), the capi
Windows fix (the ABI contract's first-ever x64 Windows execution),
CI hardening (reduction gates + floors), P1a (the ensemble
sum-invariance oracle, census B1 closed), A (bart() weights
parity), FX1-guard (student loglik+ppd), FX2 (hetero refusals +
door), B (xbart shared helpers), B2 (this slice). Every code
commit six-green; the message rule adopted and merit-rebased; the
orphaned scripts adjudicated; the equivalence recipe executed for
the first time. Open wave-0 residue: the grouped-mixing.R ratio
discrepancy (oracle-lane verification item) and the mixedMatrix
naming residual above. Wave 1 (accumulation cleanups) begins in
the recommended order C, K3, D, K2, K4, K5, G, E, K, F (fork 1:
fixed-k grid matching gaussian), H, I, J.

### B - xbart through the shared validation helpers (ef7730a7, 2026-08-17)

The four drifts closed by extracting three helpers into R/spec.R -
isBinaryFamily, estimateStartingSigma (one shared "unable to obtain
a starting estimate of sigma; provide one instead" failure),
enforceWeightPolicy (the whole binary/ordinal/nbinom policy with
the all-ones courtesy) - called symmetrically from
resolveSamplerSpec and xbart(): the weaker binary predicate, the
bare lm() error, the absent weight policy, and the
override-only-when-unsupplied resid.prior are all gone. FULL
resolver routing was deliberately NOT taken: resolveSamplerSpec's
parsePriors path would run resolveNodeHyperprior on xbart's binary
k and silently move the k-grid - slice F's fenced decision - so the
helper extraction is the correct mechanism, not a compromise. The
census's "silent-acceptance hole" framing was CORRECTED by
measurement: the C bridge's backstop already refused every weighted
binary xbart call, in its own dialect and with NO all-ones
courtesy - so the live defects were a WRONG refusal of unit weights
(now legally run, the courtesy applied R-side) and a cross-language
dialect split (the R policy now fires first; unification is K6/L's;
both dialects quoted in the gate record) - census kind 9 observed
in the wild. Defended divergences pinned by construction and test:
the chisq resid.prior default, the tree-prior grid-axes-override
rule, and binary always-override to fixed(1) even against an
explicit user prior. Gates twice: tinytest 5820/0; trio bitwise
(37/37, 12/12, 10/10); air/lintr clean; NEWS 267; R CMD check
--as-cran OK; discrimination five-way clean incl. the fixed-seed
gaussian draw-identity pin (bitwise across builds) and the
substituted sigma-failure probe (p >= n silently NaNs rather than
erroring - measured - an Inf column used instead).

### FX2 - unadjudicated heteroscedastic compositions refused (35b2b92d, 2026-08-17)

Reproduction first, per the mandate, and it halved the scope:
resid.dist = student() alongside variance = CONSTRUCTED AND RAN
silently on dbarts() and bart2() (no warning, class as normal) -
now a construction-time validation error ("a variance forest does
not support Student-t residuals: the two are not yet shown to
compose", fork 6's shape: formal stays, "does not support", fires
in 17 ms before any sampling) at the single resolveSamplerSpec site
all three construction routes share; the grouped x variance
composition is NOT REACHABLE at all (rbart_vi declares no variance
formal and its dots are rejection-only) - documented, not
dead-refused. The door memo is docs/design/heteroscedastic.md
section 15: support arrives by adjudicating the weight-channel
divisor against the scale-mixture weights (or, for grouped, adding
the same formal plus the existing collision check) and dropping the
refusal - no new surface either way. One NEWS item; one sentence in
each of man/dbarts.Rd and man/bart2.Rd, worded unadjudicated-not-
by-design. Four test arms: two refusal pins, two alive-neighbor
pins (gaussian+variance and student-sans-variance both stay legal -
and the harness's hetforce and student scenarios, which exercise
exactly those two halves, stayed bitwise). Gates twice: tinytest
5814/0; trio bitwise (37/37, 12/12, 10/10); air/lintr clean; NEWS
266 entries; R CMD check --as-cran OK; discrimination five-way
clean including construction-time verification.

### FX1-guard - student fits stop scoring gaussian (4975c20b, 2026-08-17)

Every packaged bart/bart2 fit now records a resid.dist token
("gaussian"/"student", read from the model's resid.df attribute in
packageBartResults - the single choke point; both response branches,
so the field can never be a student-only vacuity). Consumers guard:
pointwiseLogLikelihood's gaussian branch and sampleFromPPD (shared
by predict/extract for bart AND rbart) refuse a present
non-"gaussian" token - the ppd half was the mandated check finding
its own defect: rnorm noise drawn unconditionally, the same silent
class, guarded in the same shape. An ABSENT field reads as gaussian
so pre-field serialized fits keep working. The t-marginal compute
stays post-K2 per fork 2; surfaces verified: dbarts()/dbartsSpec()
return samplers, not packaged fits (guards unreachable); rbart_vi
has no resid.dist formal. Discrimination (gate-runner's own probes):
base's student-fit loglik succeeded IDENTICAL to the hand gaussian
density and its ppd drew silently; slice refuses both by name;
gaussian fits bitwise-unaffected; the deleted-field legacy shape
runs. Gates twice: tinytest 5810/0; trio bitwise (37/37 incl. the
student scenario, 12/12, 10/10); air/lintr clean; NEWS 265 entries;
R CMD check OK (the gate-runner's first WARNING was its own
--no-build-vignettes flag, not the slice). Process note: the
gate-runner FAILED the slice on hygiene - a docs/plans reference on
an added test-file line the orchestrator's own strip had missed -
fixed by a one-comment amend (4975c20b vs the gated 2a8082c4;
comment-byte delta, gates 1-7 evidence carries, the touched file
re-run 39/39). Budget overage (~4x on R lines) flagged by the
implementer and accepted: two guard sites and rationale comments,
in-mandate scope.

### A - bart() weights parity and zero-draw guard (f3f6c80c, 2026-08-17)

bart()'s tail now attaches $weights/$weights.test exactly as bart2's
does, fixing a silent wrong answer: on a weighted bart() fit the
loglik/ppd machinery computed UNWEIGHTED quantities - the
gate-runner's probe showed the base build matching the unweighted
gaussian density exactly where the slice matches
dnorm(y, ev, sigma/sqrt(w)) bit-identically. A zero (or
thinned-to-zero) draw count now refuses as "'ndpost' must be a
positive integer" (bart2's guard mirrored in this surface's own
formal spelling) instead of faulting deep in the empty-array reshape
("dim(X) must have a positive length", reproduced verbatim on base).
One NEWS BUG FIXES item; no Rd edit (no surface documents the
internal $weights component - checked bart2.Rd/rbart.Rd too); no
feature-matrix cell records the gap - checked. Two observations
recorded for later waves, not fixed here: rbart_vi's thinning guard
diverges from bart2's (== 0L vs <= 0L, different message - slice
F/L material), and bart()/bart2() draws diverge at matched settings
(surface-coherence question for wave 4, no contract claims
otherwise). Test file test-bart-weights-parity.R: 6 assertions
incl. the exact weighted-density pin. Gates twice (implementer +
independent gate-runner): tinytest 5805/0; trio bitwise (37/37,
12/12, 10/10); air/lintr clean; NEWS 264 entries; R CMD check
--as-cran OK; discrimination three-way clean. Test dense count ran
over budget on assertion formatting density - flagged by the
implementer, accepted as precision rather than scope.

### P1a - the ensemble-scale sum-invariance oracle (6bfe46bc, 2026-08-17)

tests/cpp gains test_ensemble.cpp: 200 trees x 503 observations x 30
sweeps, asserting after EVERY sweep over EVERY observation that
totalFits sums the per-tree fits and that the rolled residual equals
the working response minus every tree but the last - two independent
gates on different code (totalFits is rebuilt as y - treeY +
fits_last, so a finalize reading the wrong tree breaks the first
alone), plus a non-degeneracy check. Real symbols re-derived from
the live engine (the plan's names were descriptions):
Chain::totalFits/treeFits/residualForTesting/
workingResponseForTesting; the sweep's only aggregate writer is
finalizeTotalFits, which never sums trees - exactly why the oracle
compares against a fresh tree-order sum. Tolerance 1e-11 absolute,
not bitwise (incremental roll vs fresh sum legitimately reassociate;
worst measured deviation 1.17e-15, four orders of headroom, nine
below one leaf value). Own seeded generator, shared-stream neutral
(full-run diff shows only the one added line), green standalone
under its own filter; ~30-50 ms, full suite 15.6 -> 16.1 s. POISON
PROOF, both directions with hash-verified reverts and post-revert
touch+rebuild: (A) a 1e-3 residual drift at tree 150 fails 60
checks compounding across sweeps; (B) finalizeTotalFits reading
last-1 at numTrees > 100 fails 59 totalFits checks while sweep 0's
residual check stays green (the independence claim demonstrated).
UNDER BOTH POISONS EVERY PRE-EXISTING GATE STAYED GREEN - the
census's B1 blind spot confirmed experimentally, now closed.
Orchestrator re-ran make clean + full suite + standalone filter
independently, reproducing the reported worst deviations exactly.
71 dense lines vs ~80 budget. Landed by rebase onto 0ddc233e
(docs-only).

### CI hardening - reduction gates and assertion floors (5e1ee0a1, 2026-08-17)

Three workflow files, package untouched. exact-gates.yaml gains
hazard-reduction.R and hurdle-reduction.R (the only checks those two
shipped families have; ~0.3s each; neither reads commandArgs, so the
mode argument passes harmlessly in quick and full). sanitizers.yaml's
result post-processing gains a 5200-result total floor (~10% headroom
under the measured 5799 baseline, derivation dated in the yaml) plus
per-file floors set at 60-80% of measured counts (capi 150/232,
bcf-r5 30/45, multinomial-surface 80/126, argument-surface 60/96,
bartcore 160/250; test-simd floored at 0 with the reason recorded -
MEASURED: an exit_file'd script yields a zero-length results object
indistinguishable in aggregate from a silently dropped file, so the
real SIMD floor lives on check-standard's NEON-forced arm64 leg,
raised to the exact >= 2). Judgment call ACCEPTED at review: no
duplicate floor inside the five-leg check matrix - those legs run
tests opaquely inside rcmdcheck, a floor there would need a second
install+suite per leg, and the same-push sanitizers floor already
catches the class; redundancy, not blindness. Local gates: YAML
valid under two parsers; floor expressions dry-run against the real
5799-result object and shown tripping on subsets and on a bumped
literal; both reduction scripts exercised exactly as invoked;
paths-ignore/concurrency blocks unchanged (exact-gates' deliberate
benchmarks/** non-ignore preserved). CI: SIX GREEN on its own push -
the new gates and floors executed live and passed.

### capi-winfix - the ABI contract test runs on x64 Windows (9a299af6, 2026-08-17)

The fix for the finding below: test-capi.R now writes a Makevars
file into its build dir (PKG_CPPFLAGS = -I"<include>") instead of
passing the flag through system2's env= (platform-variant on
Windows), and probes for dbarts/dbarts.h under the include dir
BEFORE compiling (CI: stop naming the path; off-CI: exit_file) so a
missing header and a broken flag are distinguishable forever after.
The CI hard-fail on compile failure stays. One file, +13/-1. Local
gates: tinytest 5799/0 with test-capi.R at 232/232 through the new
path; Makevars confirmed tempdir-only; air/lintr clean; R CMD check
--as-cran Status OK from a clean-staged tarball; negative probe
(bogus include dir under CI=1) fails naming the path. CI: SIX GREEN
- R-CMD-check's windows-latest leg passes with the hard-fail live,
which by construction means the 201-assertion shipped-ABI contract
executed on x64 Windows for the first time in the leg's history.

### Wave 0b first CI catch - the ABI contract test never ran on x64 Windows (2026-08-17)

f6c8979d's first CI outing turned P4's hard-fail into a finding: the
windows-latest (x64) R-CMD-check leg ERRORS in test-capi.R -
"consumer.c:10:10: fatal error: dbarts/dbarts.h: No such file or
directory" - while windows-11-arm compiles and runs all the C API
assertions. The mechanism: test-capi.R passes PKG_CPPFLAGS via
system2(env = ...), whose child-environment handling is
platform-variant on Windows; the -I flag evidently never reached the
x64 compiler, so the shipped-ABI contract test (201 assertions, the
whole dbarts.h surface stan4bart/treatSens link against) has been
silently exit_file-skipping on that platform for as long as the leg
has existed - green checks, zero contract coverage. Exactly the
silent-skip class P4 was built to expose, caught on its first run.
Fix (next slice, jumps the queue - every push reds this leg until it
lands): replace the env= passing with a Makevars file in the build
dir carrying PKG_CPPFLAGS = -I"<include>" (R CMD SHLIB reads it
natively; mechanism verified locally), plus a pre-compile
header-exists probe so future failures distinguish a missing header
from a broken flag. The five other legs (three Ubuntu, macOS, ARM
Windows) stayed green.

### Wave 0b - suite repairs: P3i + P4 in-repo halves + build guard (f6c8979d, 2026-08-17)

statesAgree() now asserts PER FIELD: inst/common/stateContinuation.R
emits an info-labeled expect_identical per compared field (plus the
sigma expect_equal at 1e-8), fails by name on zero-length input, and
keeps a silent boolean mode via expect = FALSE for the one
negative-sense caller; all 11 source() sites pass local = TRUE
(LOAD-BEARING - tinytest collects helper-nested expectations only
under a local source; the census's premise-false claim was half
right and the implementer measured the real contract before writing
code) and all 17 call sites (11 files - the census's 18/15 was
miscounted) drop their expect_true wrappers. THE NULL RESULT
MATTERS: with the vacuous boolean replaced by genuine comparison,
zero call sites fail - the vacuity was real but masked no live
state disagreement. Suite grows 5591 -> 5799, the +208 reconciled
per-file exactly by the gate-runner (586 -> 794 across the touched
files). Skip accounting: tests/tinytest.R hard-fails under CI when
tinytest is missing; test-capi.R's SHLIB failure hard-fails under CI
(compiler output included) and exit_files off-CI; test-simd.R
exit_files when no non-scalar SIMD is detected instead of passing at
zero assertions; survival added to Suggests (gate-runner correction:
3 of test-aft.R's 51 assertions were survival-gated, not 51 - they
run and pass). New tools/check-build-freshness.R guards the
stale-install hazard: compares the installed package's Built stamp
(not Packaged - unpopulated by the in-place install path, measured)
against source mtimes under R/, src/, inst/, excluding
configure-generated and compiled outputs (they legitimately postdate
the stamp on a fresh build); exercised both directions. The
workflow-side floors and exact-gates additions are the next slice.
Gates run twice (implementer + independent gate-runner, own libs):
tinytest 5799/0; equivalence trio bitwise (37/37 identical-draws
lines, 12/12, 10/10); air clean; lintr zero on five touched files;
R CMD check --as-cran Status OK, Suggests addition NOTE-free;
discrimination - zero-length and per-field attribution probes fail
on the slice build and pass vacuously/opaquely on base; CI=1 run of
tests/tinytest.R from the staged copy runs the full suite. No
trailers; 17-file stat reconciled to the four parts.

### I3 - getForestVariableCounts dimnames (cd918d9a, 2026-08-17)

The live R5 read $getForestVariableCounts(forest) now names its rows
from colnames(data@x) when the training matrix carries them (five
dense lines in the method body); the raw run()$varcount channel is
deliberately untouched (the recorded separate decision). One
man/dbartsSampler-class.Rd sentence, one NEWS BUG FIXES item (the
implementer matched the file's existing "$X now Y" phrasing), and a
test block in test-bcf-r5-surface.R whose named fixture pins BOTH
halves: rownames equal the colnames AND the numeric counts are
identical to the unnamed sampler's (names added, values untouched),
plus the colname-free no-op guard. Gates run twice (implementer,
then independent gate-runner in its own lib): tinytest 5591/0;
equivalence trio bitwise - 37/37 "identical draws" lines with zero
max-|z| lines on the main harness, 12/12 and 10/10 on the other two;
air clean; lintr zero on both R files; NEWS parses at 263 entries;
R CMD check --as-cran Status OK from a clean-staged tarball;
discrimination re-proof - rownames appear only on the slice build,
values bitwise-identical across builds. Process note recorded in the
runbook: the main harness's terminal "statistically
indistinguishable" line prints even on fully-bitwise runs and is
never evidence; the per-scenario identical-draws COUNT is. Landed by
rebase onto 51c8f43f (docs-only, gated build byte-identical). CI
fired into the GitHub outage - backfill pending with the rest.

### Wave 0f - orphaned scripts and the first equivalence-leg run (2026-08-17)

Twelve benchmarks/R scripts belong to no workflow; all ran (private
lib at the a39da5d9 build, <2 min each): 8 PASS, 3 BROKEN, 1 not
standalone (bartcore-shim.R, a sourced helper). hazard-reduction.R
and hurdle-reduction.R both PASS bitwise ("reduces bitwise to the
binary fit (both links)" / "to its two standalone component fits") -
they are the only checks those two shipped families have, run under
half a second each, and are RECOMMENDED for the exact-gates.yaml
gate list (queued behind the GitHub outage; a workflow edit fires
full CI and lands only when CI can actually run). BROKEN, all
permanently absent their inputs: change-fix-instrumentation.R and
change-fix-stage2.R read CSVs only a reverted-before-commit
instrumentation patch ever wrote; parallel-falsifiers.R exits 0 with
zero data (its BC_FALSIFIER_* env hooks are inert, self-declared).
Their disposition (document-as-manual vs delete) routes to the
wave-1 accumulation lane. FLAG for the oracle lane: grouped-mixing.R
runs clean but measures the forest-ranef confounding ratio at
3.3x (K = 3) / 1.9x (K = 10), well below its own header's cited
~25x / ~3x from the prior review - 6-8 seed averages, not obviously
noise; a small verification item should decide whether behavior
moved or the header's provenance differs. The statistical-equivalence
leg (equivalence.yaml's recipe, executed for the FIRST time
anywhere): all three jobs reproduced locally, all three baselines
bitwise-identical on every scenario (37/12/10) despite 21/4/32
intervening src-touching commits - the MANIFEST's current
classifications corroborated. The YAML parses cleanly; the
main-has-no-benchmarks red-schedule hazard is CONFIRMED for the
merge-time shim design. Full report untracked at
.claude/rc-review-orphan-legs-2026-08-17.md.

### I2 - hurdle variance doc fix (4c2cdb9a, 2026-08-17)

The hurdle docs promised natural-scale predictions consume the
positive part's per-observation sigma(x) under a heteroscedastic
variance surface - a state no bart2 call reaches (the occupancy
probit's variance refusal, live at R/spec.R:392, fires before either
component fits) and one the code documents as DELIBERATE
(hurdleSigmaVec's comment, R/generics.R:989-996). Doc-side fix, both
copies: man/bart2.Rd's hurdle \item and the independent copy in
man/dbarts.Rd now state the always-homoscedastic contract (one sigma
per draw, recycled), with the heteroscedastic positive part kept in
the recorded-limitation list as a follow-up, not a promise.
feature-matrix hurdle/variance cell ? -> R spec.R:392, [f34]
rewritten as resolved, the Gaps bullet updated. One NEWS BUG FIXES
item (implementer found live precedent for recording doc-only
corrections - the offset.test entry - overriding the no-NEWS
expectation, correctly). Gates: tinytest 5588/0; R CMD check OK from
a clean-staged tarball; built package differs from base only in
NEWS.Rd + the two Rd files; NEWS parses at 262 entries
(re-verified independently); no R code touched so no
equivalence/lintr legs. Orchestrator note: the interleavables
re-verification's spec.R:370 anchor was loose - the live stop is
:392; other feature-matrix rows citing :370 for the same shared
refusal are stale by the same delta and ride the wave-4 full-
namespace resync. Prior-slice CI: five green on a39da5d9 plus
exact-gates green on the build-identical d60057b1 (the records push
fired exact-gates - its paths-ignore deliberately omits
benchmarks/** - and cancel-in-progress killed a39da5d9's own run;
an infra-flake rerun, codeload 429s, is recorded in the run
history).

### I1 - multinomial silent-drop refusals (a39da5d9, 2026-08-17)

bart2(family = "multinomial") now refuses by name the five arguments
it previously accepted and silently dropped through the
gaussian-resolved host-sampler seam: dart (flag or object), a DART
tree.prior (resolved with the same prior vocabulary parsePriors
uses, so tree.prior = dart() is caught), split.probs, variance, and
monotone - the last was HALF-applied (its proposal-probability
rewrite reached the engine while the direction constraints did not,
so draws moved). One refusal block in R/bart.R's multinomial branch
ahead of the factor/count-matrix dispatch; one man/bart2.Rd
sentence; feature-matrix multinom DART and variance cells re-pointed
at the refusal (the variance cell's R FAC:833 was a recorded value
fault - dropped at the R seam, never refused at the factory) with
[f33] rewritten past-tense and its stale engine anchors corrected
during orchestrator review (CH:4878/4813 live); a
docs/design/multinomial.md limitation sentence; one NEWS item; seven
tinytest assertions via a multinomialRefuses helper, dart covered on
both entry shapes, variance probed at K = 3 (the K = 2 host
auto-resolves probit, whose own variance refusal masks the
multinomial one - recorded so a tidy-up cannot "simplify" the
fixture back). Gates, run twice (implementer, then independent
gate-runner on the rebased commit): tinytest 5588/0 FAILURES;
equivalence trio bitwise (37/0 strict, 12/12, 10/10); air clean;
lintr clean on both touched R files; NEWS parse 261 entries; R CMD
check OK from a clean-staged tarball; discrimination re-proof - all
seven probes error by name on the slice build and silently fit on
the base build (monotone did not fail differently on base: it
silently fit). Commit hygiene verified: no trailers, six files
exactly, no plan references on added lines.
