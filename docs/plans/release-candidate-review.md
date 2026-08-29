# Release-candidate review program

Status: SPECCED 2026-08-17 (base 96ab54e0), amended per the blind
critique of the same date (verdict: execute with amendments; all
eleven applied). The TODO entry
release-candidate-review is the charter; this file is the derived
slate and the program's record. Landing notes go newest-first under the Landing notes heading.

## 1. Charter (restated)

An orchestrated review wave before declaring a release candidate for
serious human review, targeting the two defect families the per-slice
gate structure cannot see. Family 1, correctness beyond the gates'
reach: the slice batteries prove diffs neutral against agent-recorded
baselines, never that the baselines were right. Family 2,
accumulation: agent-written slices breed sediment, missed
factorization, surface inconsistency, dialect drift, and
over-complication. The charter's seed lists are a floor; this slate
was re-derived from two fresh censuses (2026-08-17, both
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
- Five of the eleven CI workflow files have never executed anywhere
  (default-branch registration, above): rchk, valgrind, sbc,
  equivalence (the statistical equivalence gate) and revdep-smoke; the
  other six run on every bartcore push. equivalence.yaml's YAML parses
  but has never run.
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
tickets; FX are
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
- The baseline re-record carrying the xbart scenario (P16, f009eff8)
  and the bart2 gaussian/probit/two-forest and mixed-matrix (B2)
  scenarios (wave-4 batch, 4a42620a; equivalence-4a42620a.rds, 42
  scenarios) LANDED, under P17's oracle-naming rule as spec'd; see
  the landing notes below.

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

### Repo-state audit and clean before the human review (018225ca..HEAD, 2026-08-25)

Scope: VD's directive that present-facing files carry no edit history
(history lives in retired plan docs and git), that memory be cleaned and
fact-checked, and that nothing be missing before the human-review and RC
legs. Four read-only audits (forward-facing files, memory accuracy, repo
hygiene, RC completeness) drove the work; their reports were session
files and are summarized here.

Design-doc anchors. An independent sample audit of 018225ca's
re-derivation found 42 percent of anchors not landing (59 live sampled;
multinomial-mutation-arc's capability table 0/37); the deltas were
content-derived, so the failure was coverage. Every anchor in the 25
live docs/design files plus feature-matrix.md was re-derived by content
(60aea703 four docs, 89e1ecd3 four, 207d4ca0 six, 5a7afe57 eleven,
27c63859 feature-matrix; 9f538bbf the sample fixes), the prose the code
contradicted was reconciled (6e148281, 30 items; multinomial-mutation-arc
now LANDED with sections 1-4 scoped as the pre-arc tree), and a
fresh-seed certification sample (n = 105) landed 97.1 percent with no
clustered file (69cb6491). feature-matrix.md lost its refresh chronology
(b2e360c6, ba0e79ba; 1255 -> 1042 lines) and two landed fixes it still
listed as open gaps. After the code fixes below moved lines, a
diff-derived line map re-aligned 87 anchors (7bf3c600); anchors are
verified against the code tree at e654012d.

Forward-facing strip (39fa379a, 9a908bbd): TODO 583 -> 329 lines (nine
closed-narrative entries deleted, ten trimmed, three stale claims
corrected - the anchor resync had landed, the data@x dimnames drop was
fixed, the bcf/multinomial pins live in equivalence.yaml); INDEX.md
statuses re-read from each doc; docs/README.md and docs/plans/README.md
rules rewritten to describe practice; workflow comments stripped of
history. The TODO's claim that the five schedule-only workflows 404 on
manual dispatch was verified TRUE and kept. Baselines pruned 65 -> 11
files (01bd49a7; 18.6 MB), MANIFEST rows with them;
two dead .gitignore lines dropped.

Shipped fixes (654e4681, 3298da0e, ae5b91d8):
man/dbartsSampler-class.Rd documented wrong argument names for
setControl/setModel/setData and omitted storeState's ptr - R CMD check's
codoc cannot see reference-class methods, so tools/check-rc-codoc.R now
guards them in lint.yaml (mutation-proved); man/pdbart.Rd's example
re-wrapped in \donttest; six plot sites restore par() (CRAN policy; five
tinytest assertions rewritten to a sentinel); R/diagnostics.R named a
file that does not exist; sanitizers.yaml now installs survival and
posterior, so seven test files run under ASAN (7279 results, was 7246
there); NEWS.Rd's three self-contradicting pairs resolved and three
surfaces plus two warning classes introduced. The sanitizer step then
failed on ae5b91d8 from an apostrophe inside its single-quoted R block
(comment strip); cured at 674be0c4 - fire-and-cure, the red run stands.
e654012d reunified the forest-weights refusal text across the boundary
with the bridge file's line count held.

Record: the review's evidence is tracked at docs/plans/review-2026-08-24/
(658869ac) and the tour refreshed and tracked at
docs/plans/bartcore-review-tour.md (0b89ab8b, anchored at ae5b91d8);
session-path citations retargeted or stripped (4c018187).
tools/check-doc-freshness.R now checks every docs/design anchor at the
line (strict tier for identifiers beside the anchor, advisory otherwise),
commit-hash resolvability (skipped on shallow clones, cross-repo cites
exempt) and INDEX status labels (3ef2190f; ff7fee0a fixed what it found;
70 advisories remain, listed by shape in the commit). air excludes docs/
(3a148caf).

Memory: 156 claims checked, eleven false (the x86 bench box is reachable;
the baseline trio and counts; dbarts.h's last change; a phantom CI
workflow); all files rewritten forward-facing.

Doors left (unchanged unless noted): fit-time test-basis; logistic
weights never ride saved state; DESCRIPTION Date bumps at the RC; the 70
freshness advisories; plot.bart/plot.rbart still leak par() through
plotSigmaTrace (the six own-class sites restore; plot.bart/plot.rbart fixed fcbbc478).

### The docs/ citation strip and xbart's k-grid sort (ad4d801d, e1223266, 2026-08-25)

Citation strip (ad4d801d, three commits): the maintainer's ruling that shipped files cite
no docs/ path of any kind - the tarball strips docs/, so a path is a pointer to something
not shipped - executed as a READING pass over 316 sites in 81 files (R/ 13, src/ 18, inst/
49 incl. the tests and dbarts.h): every citation clause dropped, every constraint kept in
place. The independent fidelity read sampled ~180 sites (every src/bartcore and bridge
hunk in full) and found zero load-bearing statements lost; not one comment needed a rule
inlined from a design doc - each citation sat beside a self-contained statement. Also
cut: the bridge's state-format history narrative (the registry rule and the attributes
clause stay and agree with stateFormatVersion = 3); the 14 review-labelled test sections
renamed to what they pin; 13 bare design-doc filename mentions swept under the same
ruling; a WG14 standards URL that merely contains "docs/" restored. feature-matrix.md's
251 affected anchors re-derived by content (deltas 0..-7, non-uniform; a 15-anchor audit
plus a mechanical check of all 463 resolvable anchors: 446 byte-identical targets, 17 on
lines the strip itself edited). Census rule going forward: `grep -rn "docs/" R src
inst/include inst/tinytest inst/common man | grep -v http` is empty, and so is the bare
`*.md` grep. Bitwise 43/12/11; tinytest 7240/0. FOUND, NOT FIXED HERE: feature-matrix.md's
anchors had rotted BEFORE this pass (its legend is stamped ad4a131b; ~9 verified cells
name a construct their line no longer holds - the freshness check verifies presence, not
construct identity), and ~705 file:line anchors in the other design docs and plans are
verified by nothing; a whole-file content resync of the matrix, stamped at the tip, runs
next, ahead of the review-tour refresh.

xbart k-grid sort (e1223266, three commits): 0.9-34 sorted a numeric k grid decreasing
before the sweep and un-permuted the reported axis, so the loss array was invariant to the
order k is listed; the rewrite had dropped it, and a cell's loss shifts toward its sweep
predecessor through the warm start (measured t = -6.6 at short burn-in, 0.27 percent at
the shipped n.burn). Restored: the sort where k becomes a numeric grid (a sampled k is
untouched), the un-permute on the full array before the drop so the array and its dimnames
agree, one Rd sentence saying exactly what the sort buys - order-invariance, not per-cell
unbiasedness - and pins: identical() arrays for k = c(2, 8) and c(8, 2) after alignment
(red against the previous tip, 6 percent disagreement), caller-order dimnames, no k axis
for a hyperprior k. test-reproducibility-xbart.R's snapshot regenerated (its grid now
sweeps 8 before 4). THE PROGRAM'S THIRD BASELINE RE-RECORD: equivalence-736bfb05.rds
(recorded at the sort commit 736bfb05, the engine tip; reproduces itself 43/43 under
--strict-coverage; four places updated in the recording commit: equivalence.yaml, MANIFEST
with the P17 row, feature-matrix [f39], TODO). Partition against 5a3bc276: 42 of 43
scenarios bitwise, xbart alone moves (8 summaries, max |z| = 1.94; its grid is k = c(1, 3),
now swept 3 then 1). Adjudicated independently: the new order's c(1, 3) is BITWISE the old
build's c(3, 1) on all 40 seeds - a pure relabelling of the previous sweep - and a 40-seed
paired per-cell comparison puts the pooled difference at -0.019 (SE 0.014, 0.8 percent of
the loss), with the k = 1 cell's warm start shifted toward its k = 3 predecessor
(-0.044, t = -2.6) and the k = 3 cell's not (-0.006): the target is unmoved to
Monte-Carlo precision, and the warm-start bias is real at n.burn[2] = 5, which is why the
Rd disclaims unbiasedness. The MANIFEST row says "adjudication of the expected class"
verbatim, as the P17 rule requires. bcf 12/12 and multinomial 11/11 bitwise at the same
tip; tinytest 7246/0. The adjudicator also corrected the TODO, which had listed the
already-landed sigma cure as open.

### Threaded predict landed, and the constant-response sigma cure (e0e59097, b2df6522, 2026-08-25)

Sigma cure (e0e59097, three commits): estimateSigmaFromLinearModel now routes a non-finite
linear-model sigma (no residual degrees of freedom - n = 2 with a predictor, and n = 1,
which dbartsData does not refuse) through the marginal residual sd with a classed
dbartsSigmaFallbackWarning (the sparse branch's shape, which gained the shared class), and
then floors whatever value results at sqrt(.Machine$double.eps) times max(1, max |residual|)
- so a perfectly-fit constant response no longer installs whatever rounding noise the
host's BLAS kernel leaves (exactly 0 on SSE2-class kernels, the clang-asan red recorded
below). Pins in test-boundary-inputs.R and test-multinomial-r5-surface.R with
tolerance = 0 (all.equal's default tolerance had made the exact-floor pin vacuous against
the old 1e-15 value - caught by the implementer's own red-then-green); the n = 2 case, which
used to be an accidental "sigma estimate is NaN" refusal, now runs with sd(y). Bitwise
43/12/11, tinytest 7194/0, multinomial-mutation-arc.md's measured 1.52e-16 reworded.

Threaded predict (b2df6522, three commits: engine + surface, the design doc's landing
note, a reviewer comment fix): docs/design/threaded-predict.md is LANDED. predict's
n.threads is real: Sampler::predictColumns partitions the (chain, draw) slab list across a
std::thread fan-out mirroring run() (SIGINT blocked around the spawn on POSIX, per-worker
scratch, every body catches and the caller rethrows after the join, no cross-thread
reduction - so the replay is bitwise identical at every worker count, and the inline
below-cutoff path is the same body); numThreads = 0 means the sampler's own count floored
at 1; dbarts_sampler_predict gained size_t numThreads and both API-hash literals were
re-baked (DBARTS_C_API_MINOR NOT bumped: the header's own rule keeps the constants still
before the first release, and a parameter gained is breaking, not additive); the R5
predict/predictForests default to the fit's thread count, refuse a non-positive / NA /
non-numeric value by name, and the formal is appended LAST on all six generics (bartCause's
positional predict(fit, x, group.by, combineChains = FALSE) still binds);
test-generics-multithreaded.R is a real pin (44 assertions: bitwise identity across thread
counts on gaussian, binary, sparse, variance, multinomial, BCF glue and per-forest replay,
plus a deterministic partition report - resolved count, worker count, block contiguity -
through a test-only channel with a writable traversal cutoff so unit-sized fixtures
exercise the fan-out); tests/cpp covers counts 1/2/3/7/64, the setNumThreads(0) floor and
the inline path, and test_facade.cpp the arity crossing; test-capi.R's consumer.c passes 0
and had been SILENTLY SKIPPING (264 assertions) since the prototype changed - fixed. Three
mutation classes caught (a worker skipping its last slab; a partition-dependent tree-sum
order; the resolved count forced to 1). ASAN+UBSAN and TSAN on tests/cpp: 0 diagnostics.
Same-process scaling probe (n = 1000, 200 trees, 500 draws, 2000 rows, 5 alternating
rounds): 1.92x at 2 threads, 3.72x at 4 - consistent with the memo's serial-fraction
bound, still information rather than a gate. Consumer lockstep: stan4bart's one call site
passes 0 (its own bartcore branch, committed there, not pushed); bartCause needs no edit.
Reported, not fixed (mirrors run()): a throw from std::thread construction itself would
leave earlier workers unjoined; Rf_error after the join skips the frame's destructors; a
fractional n.threads truncates silently. 0.9-31's NEWS entry is corrected in place.

### Second whole-branch review: fix wave 3 - the judgement-gated slices (42b12ac7..7ad0bbea, 2026-08-25)

Scope: every slice the review left gated on a maintainer judgement (the nine groups
decided 2026-08-24/25) plus the three questions the report still owed the maintainer,
answered this session: the three engine divergences are FIXED (all three, no
re-record); shipped code cites NO docs/ path of any kind (docs/design included: the
tarball strips docs/, so a path is a pointer to something not shipped - the strip is a
separate reading pass, below); xbart's k-grid sort is RESTORED (0.9-34 parity,
order-invariance; its equivalence re-record lands last, after the threaded-predict
arc). Decision substance with probes: docs/plans/review-2026-08-24/decision-brief.md.

Landed, in order, each implemented in its own worktree off the then-tip, each
independently gate-run (--preclean private lib, tests/cpp full + filtered where src/
moved, full tinytest, trio 43/12/11 bitwise, air, lintr, R CMD check --as-cran from a
clean staged copy) before its push, CI six-green on every push (one cure below):
- 42b12ac7 + 33837351: every `default:` arm on a ResponseFamily switch deleted, the
  bridge's open-coded family chain made an exhaustive switch; families that cannot
  reach a site are listed as fall-through labels with the reason. The package build
  only warns (CRAN forbids -Werror in Makevars), so the hard gate is tests/cpp's
  Makefile (-Werror=switch -Werror=return-type) and a cpp-tests CI step writing
  `CXX20FLAGS += ...` to ~/.R/Makevars - CXXFLAGS is inert there because R passes
  CXXFLAGS='$(CXX20FLAGS)' on the make command line under CXX_STD = CXX20; verified by
  a temporary seventh enumerator failing the install. Bitwise.
- cd615a1e: the three divergences. Variance-leaf positivity is now one predicate at
  BOTH warm-start install arms (installTrees on a hand-edited state used to take a
  non-positive scale leaf that setState refused; fits blew four orders of magnitude);
  the amplitude combiner's shippedShape() reads the prior as well as the basis shape
  (a positive half-Cauchy scale past forest 0 forces the general sweep; every shipped
  route stays on the two-scalar path, bitwise on all three baselines; a tests/cpp case
  prints forest 1's variance moving against the frozen 0.5); recoverTreeParameters and
  applyNewData assert the single-forest precondition and say what a multi-forest lift
  needs (bases in the same call). The combiner.hpp static_assert residue is swept.
  Local ASAN+UBSAN tests/cpp leg: 0 diagnostics.
- 4e5bcebf (four commits): xbart oracles above the loss call -
  test-xbart-fold-oracle.R pins replication 1's fold rows reconstructed OUTSIDE xbart
  from the seed, disjointness and coverage, per-fold alignment of y.test against the
  test channel versus a permutation null, no leakage (held-out rmse at or above sd(y)
  with the in-sample fit far below), and axis placement (the k = 20 slice at sd(y),
  equal-length axes so a transposition fails); test-xbart-reproducibility.R pins
  chunk-1 prefix identity across n.threads = 1 vs 2 (the 1-vs-4 arm needs more workers
  than CRAN's core limit allows). Five planted mutations caught; the leakage band was
  re-calibrated by a 20-seed sweep after the Opus read found 10 percent false alarms
  at n.test = 12 (now n.test = 30, threshold 1.0 sd(y): clean 20/20, mutant 20/20).
- 2ef009a7 (three commits): the R surface. bart2()/rbart_vi() lose their
  rejection-only `...` (rejectUnknownDotsArgs, retiredDotsNames and the rngSeed sunset
  row deleted with it - a recorded loss; R's own "unused argument" is the wall);
  bart() refuses all ten bart2 tokens BY NAME echoing the token typed (twopart named
  twopart, never its alias), and dbarts() captures the typed family ahead of the fold
  and the hazard remap so every downstream refusal echoes it; ONE offset name
  everywhere with shape by family (multinomial requires n x K, a vector is refused as
  the softmax's null direction; every other family refuses a matrix), dbartsData's
  offset.category / offset.category.test formals are gone (slots and R5 method names
  keep theirs), and the false "same number of observations as 'y'" message for a
  Surv / n x 2 / n x K response is a by-kind refusal on the matrix, sparse AND formula
  branches; monotone loses "inc"/"dec" (the "0" token STAYS - a positional character
  vector needs it; the case fold is documented); one shared zero-sample refusal for
  bart2/xbart/rbart_vi with a dbartsControl.Rd sentence; $n.chains on every
  keepSampler fit, names(fit) no longer lists NULL components (forest.labels carried
  across the filter by hand), dbartsDrawLatents accepts its own default,
  xbart validates family before ingesting the response and n.threads' length,
  defaultNodeScale refuses an unknown family; dataSlotOrNULL deleted (bare slot reads).
  Reviewer fixes folded: dbarts(family = "twopart") echoed the resolved token; two
  offset messages named the removed arguments; a NEWS claim about "0" was false.
- fdcbabe5 (two commits): the own-class generics. bartMultinomial/bartOrdinal/
  bartNegbin/bartHurdle honour combineChains and ci.level and refuse forest /
  contribution / sample-on-fitted / type-on-residuals / vars-on-summary BY NAME (the
  reviewer closed 13 further swallows on fitted/residuals/predict); plotTree and
  survivalProbabilities refuse by name instead of "no applicable method";
  extract(type = "trees") on a keepSampler-only fit follows plotTree's disclosed
  fallback, silently; predict.bartMultinomial's formal is `offset`; each family gets
  plot (siblings of plot.bart: scalar-parameter trace left, observed-vs-fitted interval
  right; ordinal degrades at K = 2; multinomial branches to observed proportions on
  multi-trial rows; hurdle is 2x2 with the composed E[y|x]), extract(type = "loglik")
  (dnbinom at the per-draw dispersion; log of the stored Phi difference; dmultinom
  WITH its coefficient, unit = the observation row; hurdle log1p(-pi) at zeros and
  log(pi) + dlnorm at positives on the NATURAL scale with the Jacobian, no truncation
  because the positive part is lognormal - the Rd carries the four statements a porter
  needs), and as_draws_array/as_draws_df (scalar parameters only: dispersion;
  threshold[k] with the pinned first; meanProb[level]; occupancy./positive. prefixes).
  The survey behind the semantics (BART, brms, rstanarm, pscl, ordinal, nnet, VGAM;
  formulas verified numerically) is docs/plans/review-2026-08-24/generics-survey.md, the
  ruling generics-phase2-spec.md. The reviewer's own oracles agree to <= 3e-15 on all
  four families, combined and uncombined. Four planted formula mutations caught.
- 7ad0bbea: the 31 low-level bartcore* handle wrappers moved out of the namespace into
  inst/common/bartcoreHandle.R (same names, bodies the .Call on dbarts:::C_ symbols;
  three keep in-package copies and alias them), 559 call sites in 41 test/benchmark
  files prefix-stripped by an idempotent script and each file sources the helper;
  Tree::rightChildOf and Sampler::setCurrentSampleNum deleted; adoptPointer and
  reapplyForestWeights (3 live call sites, NOT dead) documented on the R5 class;
  core-generalization.md's "nothing dispatches per observation" amended to what
  facade.hpp does. Consumers below the R5 class use the shipped header; there is no
  R escape hatch. tinytest total unchanged by the move (7181 before and after).

Plan and ledger: docs/plans/review-2026-08-24/wave3-plan.md (Opus; every count
tool-verified; its 714-site figure counted the creators, the true scope is 559).
tinytest 7040 -> 7181. Baselines UNCHANGED; no re-record; dbarts.h untouched.

CI fire-and-cure: the clang-asan sanitizer job went RED on fdcbabe5 - a PRE-EXISTING
pin (dbarts() on a constant response, test-multinomial-r5-surface.R:91) hit the
bridge's "sigma estimate must be greater than 0" - and GREEN on a rerun of the same
commit and image. Diagnosed in the r-hub clang-asan container (OPENBLAS_CORETYPE
sweep, fresh processes, allocator offsets): the starting estimate
summary(lm(y ~ x))$sigma for a perfectly-fit constant y is rounding noise whose
value is a BLAS-KERNEL property - ~3.8e-16 on AVX2-class kernels, ~2.2e-16 Nehalem,
EXACTLY 0 on SSE2/SSE3-class kernels - and the two runs landed on different runner
hardware (the failing run's per-file times were uniformly 13 percent faster). Not a
commit regression (the five sigma-path files are identical across the two commits;
main carries the same unfloored estimate and the same > 0 gate). The fix - a relative
epsilon floor in estimateSigmaFromLinearModel, host-independent, with the existing
"indistinguishable response" warning kept as the user-facing signal - lands as its
own slice after this note (cured at e0e59097, the note above). The red run stays red in
history by design. A second
cure, same day: lint went RED on 7ad0bbea - 46 object_usage_linter hits in 10 test
files, because the handle names now come from a sourced helper that the linter cannot
see inside closures (per-file lintr on the two edited package files was clean; only
CI's tree-wide lint_package() showed it). Cured at ee7bb783 by a .lintr exclusion of
that one linter for the two helper-sourcing directories, listed per file at lint time
because the installed lintr drops per-linter scoping on a directory entry; an
.onLoad-side globalVariables hook was tried and rejected as lint-only machinery in
shipped code. Lesson: any slice that moves names or touches many R files runs
lint_package() as CI does, not per-file lintr.

Design settled alongside (not yet implemented): docs/design/threaded-predict.md -
predict's n.threads WIRED (maintainer ruling): per-(chain, draw) slab partition,
structural bitwise identity at every thread count, dbarts_sampler_predict gains
numThreads (stan4bart one line, bartCause none under the append-last rule), honest
Amdahl ceiling rather than a claimed speedup. Memo, revision 2 and the 35-finding
critique: docs/plans/review-2026-08-24/memos/threaded-predict-memo.md,
threaded-predict-memo-r2.md, and threaded-predict-critique.md.

Doors and residue recorded here, not in code: negbin rootogram panel; ordinal
log_diff_exp tail precision; a negbin burn-in dispersion channel (bart2 negbin runs
one run(n.burn, n.samples), so plot draws no burn-in segment); par(mfrow) is not
restored by any plot method (overbroad when written - the six own-class sites already restored, line 610; plot.bart/plot.rbart's half fixed fcbbc478);
as_draws_array.bartMultinomial(vars =) ignores a non-meanProb value (documented);
predict(bases =), predict(sample =), fitted(combineChains =) on the own-class fits
are still silent no-ops (fixed d48aef8a, refused by name; n.threads real since the
threaded-predict arc); dbartsData(counts =) with a formula reaches its own later
refusal; the `# ---- fix NN` section labels in tests (14 sites) and the 162
pre-existing docs/design paths in shipped files are the strip pass's scope;
setForestBasis(k, ~var) evaluating the formula in environment(basis) and plotTree's
dead padding branch stay out (no judgement names them).

Lessons: tinytest 1.4.3 has no `!` method on its result object - count failures with
sum(sapply(out, function(x) !isTRUE(x))); a leakage-style statistical pin needs a
multi-seed sweep at design time (the 5-seed check passed while 20 seeds showed 10
percent false alarms); a nondeterministic CI red is diagnosed by a same-commit rerun
BEFORE any bisection, and a green rerun does not close it - the flip is the defect;
strip scripts using process substitution need bash, not /bin/sh.

### Second whole-branch review: lenses, findings, fix waves 1-2 (b102e17c..07ad73e4, 2026-08-24)

Design (review-lenses-memo.md: eight agent-code failure modes, in-repo evidence; VD's
six decisions): cross-surface consistency first, since findings here have usually been
an ABSENCE across surfaces, not a bug inside one; run the matrix and gate ledger over
the WHOLE surface, not a diff; measure test strength by mutation testing, scoped first
to touched files then widened to a full-suite replant after a spot check found
survivors the narrower scope had hidden; four legs in parallel; reserve fill-vs-refuse
calls and all deletions for VD; verify-triage-one-fix-wave, draw-moving fixes last, a
calibration lane alongside.

Consistency legs: matrix.R reproduced the builder's 877-cell grid byte-identical;
matrix-review-entries.md added 655 cells (1532 total), swept 314 Rd claims (71 agree,
8 contradict); matrix-review-generics.md ran 4110 generic-execution cells over 40 fit
variants x 11 generics. reading-R/engine-list.md read R/ and the
engine/bridge/support libs for duplication and stale comments, candidates only.
gate-ledger.md inventoried every workflow and baseline (64-row MANIFEST, 11 workflows
on bartcore vs 1 on main, 5 declared statistical gates, zero CI runs);
gate-ledger-read.md re-ran it independently, confirmed most, corrected overstatements,
named misses it missed (the C API's own ABI-hash gate, the untracked evidence
substrate).

Mutation legs: A (85/84 scored, the 65 touched tinytest files); B (80/63 caught, 15
gaps, touched tests/cpp); C (132/129) and D (101/66, budget stop) split the 102
untouched files, each narrower than full-suite. A spot check found leg C/D
"survivors" a full-suite replant actually killed, so all 31 outstanding zero-killers
were replanted whole-suite: 16 survived (confirmed reach gaps), 15 killed (narrowed or
refuted).

Calibration lane: anchor-main.md built released 0.9-34 and ran the equivalence
harness's statistical mode against 1.0-0 over 20 scenarios (10 high-precision) - zero
unexplained disagreements; every |z| > 4 traces to a documented change (change-move
balance, zero-weight sigma df, rbart_vi re-anchoring, xbart fold leakage, binary-k
default). calibration-sbc.md re-ran sbc.R at n.trees = 200: 75/83 functionals pass,
all 3 Bonferroni flags pre-adjudicated (an nbinom mixing ridge, a grouped-gaussian
harness artifact); aft and heteroscedastic stay uncovered at ensemble scale -
leaf-scale-pin-memo.md found that fix is already-shipped surface.

consolidated-report.md re-ran every claim against its own build, replanted 63 more
mutations: CONFIRMED 8 BLOCKERS/44 MAJORS/35 MINORS; REFUTED 5; QUALIFIED 6; NARROWED
6; 1 fixed pre-review. Bins: 44 agent-fix entries plus 15 minor one-liners, 9
VD-judgement groups (the other 8 MAJORS, ~260 grid cells, 8 deletion/seam items), 11
defers.

Five most consequential defects, plainly: (1) extract(fit, type="trees", sample=1)
silently returns a DIFFERENT sample's trees, no warning; sample="train" crashes on an
internal NA. (2) survivalProbabilities() fails on every hazard fit whose predictors
carry column names - the normal spelling. (3) fitted()/residuals() on an rbart fit
SEGFAULT when the random-effect array has no dimnames, reachable by one assignment.
(4) dbartsData() falsely claims mismatched row counts for every two-column response;
four surfaces inherit it. (5) the facade, the one hop between the C API and the
engine, had no conformance test: 5 of 7 planted forwarding defects passed the whole
C++ suite.

Fix wave 1 (b102e17c..fe505ae3, five slices, merged battery PASS, CI six-green):
8042cc2c fixes B1 (extract's trees call refuses
sample/combineChains/forest/contribution by name) and B2 (survivalProbabilities names
the period column both branches); b657e8ae corrects eight Rd/comment claims the matrix
review falsified (bart.Rd shapes/defaults, dbartsSampler-class.Rd, three stale engine
comments, a phantom configure switch); 7318b266 hardens rchk.yaml, adds [main, master]
to sanitizers.yaml, scopes lint/pkgdown's PR trigger, deletes three unrunnable
harnesses, wires two drift checks into CI; 52d3b5ff pins ten R-side reach gaps, one
assertion each; fe505ae3 adds the facade conformance suite (test_facade.cpp, 59 rows,
one per SamplerBase virtual, through the base reference not Sampler<L>) plus direct
pins for B3-B7, tests/cpp 257 -> 266 cases. Gates: tinytest 6984/0; tests/cpp 266/0,
0 ASAN+UBSAN diagnostics; trio 43/12/11 bitwise, no re-record; --as-cran 1 NOTE.

Fix wave 2 (fe505ae3..07ad73e4, three slices, merged battery PASS): e35c8797 refuses
an unmatched or all-NA group index before .Call in fitted()/residuals() on rbart and
bounds-checks it C-side too, closing a SEGFAULT reachable by one assignment (new
BLOCKER B8); 66ac05b3 makes logistic-reference's 200-tree arm - the only true
ensemble-scale comparison against an independent implementation on the push path -
actually gate, and installs BART in CI so it is not skipped; 07ad73e4 pins 21 more
files of reach gaps and retires test-sampler-state-emptyLeafVeto.R's false regression
claim - the invariant lives in tests/cpp/test_moves.cpp instead, proven there by
planting the same mutation. Gates: tinytest 7040/0; tests/cpp 266/0, 0 ASAN+UBSAN
diagnostics; trio 43/12/11 bitwise, no re-record; --as-cran 1 NOTE.

What the gates and anchor establish (gate-ledger-read.md sec 4, faithfully restated):
the green checks are strong evidence of memory safety, portability, API stability, and
per-step engine math against independent brute-force oracles on every push, and of a
flat C API that cannot change a signature, ABI enum, or struct layout without failing
the build - and NO evidence about calibration at ensemble scale or behavior on any
other machine, since equivalence/SBC/rchk/valgrind/revdep-smoke have zero CI runs and
every baseline descends from bartcore's own output.

Pending. Nine VD-judgement groups, one per message: VD-A family-spelling parity; VD-B
the unknown-argument diagnostic; VD-C dbartsData's offset.category/survival spellings;
VD-D own-class generics' argument vocabulary; VD-E which generics an own-class fit
carries; VD-F undocumented option vocabularies; VD-G the four ResponseFamily default:
arms; VD-H deletions and seams; VD-I predict's n.threads formal. Gated: S3
(VD-A/VD-C), S7 (docs/ citation policy), S8 (M11's three draw-moving divergences, last
by design), deletions (VD-H). S10 (leaf-scale-pin-memo.md's shipped-surface pin for
the aft SBC arm and grouped-sigma artifact) is in progress. The xbart oracle slice
(fold/axis/thread probes, R-only, ~1 Sonnet session) is recommended pre-RC, not built.
The review-tour refresh and VD's human review follow.

Lessons. The spot-check rule paid off: replanting a handful of leg C/D zero-killers
whole-suite, instead of trusting the touched/half-suite score, found survivors the
narrower scope had missed - why the consolidator widened to all 31 and replanted 32
more. Assertions that cannot fail keep recurring: test-plot-generics.R:103-104's
expect_true(is.character(capture.output(...))) literally cannot fail, and print.bart's
family label was wrong under it. The facade is a sentinel-pattern gap like the
ordinal deep-grow segfault the last cycle found: every other test drove the engine
directly, so nothing proved the one hop safe until a conformance suite was built on
purpose to reach it. A refute-posture consolidator earns its keep: re-running three
of the lenses memo's own claims (rbart_vi thinning, n.threads length-2, SamplerBase's
virtual count) refuted all three - the discipline that kept this pass off phantoms.

### Residue burn-down and value-scan defect fixes (a0eaf348..044a9098, 2026-08-24)

Two independent censuses -
docs/plans/review-2026-08-24/memos/backlog-value-scan-2026-08-24.md and its blind
critique - re-verified every open TODO entry against live code and found four
undocketed defects (heteroscedastic loglik/PPD, the variance setState column mask,
dbartsData(bases=)/dbarts() alignment, no summary method for ordinal/nbinom/hurdle).
Three memo/critique arcs (weights-saved-state, multinomial-forest-replay,
fuzz-grow-mask) plus one independent adjudication of the ordinal draw-law change
round out nine slices: design memo -> blind critique -> implementer spec ->
orchestrator diff review -> (adjudication, for the draw-law change) -> cherry-picked
stack, five hand-merged append-point conflicts -> one merged-tree gate battery.

a0eaf348: getPointer/setState bound a re-created engine to selfEnv$pointer
BEFORE its state installed, so a refused install left a live unfitted engine that
the next run sampled from and storeState overwrote the fit with. Fix binds only
after install succeeds (test-sampler-state-format.R, 67 lines).

221ec7af: loglik, ppd, and summary all scored a heteroscedastic fit at the
scalar sigma - fixed, no posterior content - instead of s(x): measured summed loglik
-3592.1 reported vs -2031.7 correct. Fix scores/draws at s(x)/sqrt(w), reports
mean.s, and refuses ppd at test rows with no s.test (214-line test file).

c95a5e83: setState held only mean forests to installForests' column-mask rule;
reproduced, 32 forbidden variance splits installed, 12 still live 5 sweeps later.
Fix shares columnMaskStateFeasible across every forest (C++ and R pins: refused with
installTrees' own text, recipient left unchanged).

99317fec: multinomial/ordinal/nbinom generics used bare match.arg instead of
validateType, so type=response/link gave reason-free errors; multinomial's
type=bart/forest refusal now names itself on predict/fitted too, and
predict.bartMultinomial gains offset.category.test. Closes the
multinomial-forest-replay door as a DECISION - every candidate use collapses to
log(predict(...)) to 1e-13 or needs a declined identification; re-open trigger is a
level-sensitive K-latent composition.

b4b9119d + 044a9098: the ordinal scan dropped missing rows from
the split likelihood while the no-split term kept them, biasing toward not-splitting
under n.grow.sweeps; fixed by scoring each missing direction over the rows it seats.
Independently adjudicated to a max abs diff of 1.53e-16 against an R re-derivation
across five fixtures, chi-square rejecting both the dropped and halved-mass laws.
The adjudication also caught a sentinel gap the passing suite never saw - it
segfaults a deep-grow probe - closed in the amendment.

47cdb96a: dbartsData(bases=)'s formula path checked the post-subset row count,
dbarts()'s the full-data count; at equal counts the two silently accepted different
bases for one call - dbartsData now applies the dbarts() rule. 1583140b:
ordinal, nbinom and hurdle.lognormal fits fell through to summary.default's
Length/Class/Mode table; they now get real summaries.

d788bfef: OP_GROW was absent from five multi-forest and four single-forest
fuzz configs. Widening the mask alone is vacuous - the op ends in a
getState/setState round trip that heals corruption first - so the fix checks
invariants before that round trip. A planted defect is caught at seed 1 on all five
multi-forest configs, missed by the stock mask; a legal reordering control stays
green at 60 seeds.

da82ec23: logistic's PG latents are weight-shaped but no state carried the
weights that shaped them, so a restore paired stale latents with new counts for one
sweep. Fix: a weights.digest (byte hash, since moments collide on unequal weight
vectors) rides the state; setState re-derives latents through setWeights on a
mismatch, gated so a matched restore stays the identity. Deleting the write, the
compare, or digesting the wrong side each turns a distinct mutation pin red;
dbarts.h changed in prose only.

Gates, final merged-tree battery at 044a9098: tinytest 6946/0; tests/cpp 257/0 full,
sampler filter, and 200-seed fuzz green; ASAN+UBSAN full suite 0 diagnostics;
--as-cran 1 expected NOTE; pkgdown and doc-freshness clean; every slice trio-bitwise (43/12/11), no baseline re-recorded.

Decisions (VD, 2026-08-24): residue burn-down first, then a second adversarial
whole-branch review (discussed first), then the review-tour refresh. Fit-time test
basis DEFERRED POST-1.0 - bartCause's counterfactual arms are already served at fit
time, unlocking nothing new. group.by, survival entry=, sparse-extensions' two
halves, and rbart_vi's logistic token deferred post-1.0, window not lock-in.
Approximate Polya-Gamma DECLINED for 1.0; hurdle samplerOnly stays refused; the
twin-create deletion is STRUCK as a relitigation.

Doors left: the active-row mask still does not ride saved state (bounded,
documented, not the hybrid closed here); per-forest multinomial replay stays
refused, re-open trigger above; the four additive items above; fit-time test basis,
post-1.0; variance-forest-mutation-routing.md's other two doors (scale-leaf
staleness, prior-predictive) untouched.

Lessons: append-point conflicts concentrated in inst/NEWS.Rd, long Rd lines, and
R/diagnostics.R. The sentinel gap is a reminder that a passing suite is not
evidence of an assertion's reach - found only by building a probe the landed tests
had no reason to write. The value-scan pattern - scan, then a critique that REFUTES
rather than confirms - surfaced four silently-wrong exported channels no per-slice
gate had reason to look for.

### Pre-RC burn-down (2026-08-24)

The RC gate exists to prepare for the serious human review the charter
names in section 1, so before the coverage/mutation/coherence passes,
the content-heavy backlog went first. Four slices landed: 53525f4d,
the multinomial fuzz arm (docs/plans/archive/multinomial-counts-mutation.md);
139a1976, the predict-time blend
(docs/plans/archive/bcf-bartcause-relocation.md, 2026-08-24 note); d0701a6a,
the logistic weight channel (docs/design/r-c-division.md's "The
latent-family weight channel"; TODO's latent-family-weight-channel
entry); and 124259d0 + 41661523, the saved-tree store cursor fix plus
pdbart, recorded here (no plan file owns the pair).

The store's write cursor carries across `run()` calls; all six readers
- predict, predictPerForest, predictVariance, getTrees, printTrees,
the warm-start donor pool - walked slots 0..capacity-1, so a second
RECORDED run replayed draws rotated against the recorded-draw channels
(burn-in itself is never recorded;
docs/plans/review-2026-08-24/memos/tree-store-burnin-memo.md
has the corrected mechanism). Fix: draw i = slot (cursor + capacity -
filled + i) mod capacity, filled = min(recorded, capacity), oldest
first; a store with no recorded draws refuses. recordedDraws is now a
required state field (format 3, readable floor 3); numSavedSamples
reports recorded draws, not capacity; dbarts.h changed in prose only.
Option C (reset the cursor per run) was rejected - bart2's nbinom
records each kept draw in its own `run(0, 1)`, depending on the carry.
`pdbart(keeptrees=TRUE, nskip=0)` hit the same bug:
`bart(sampleronly=TRUE)` hands back an unrun sampler and pdbart's
keepTrees branch never ran the sampling phase, returning a flat
surface; fixed in the shared prologue, both surfaces now bitwise
identical at a seed. Price 253/335 raw (+ pdbart 6/63); mutations
M1/M2/M3 (11/3/2 fails) discriminate; gates tinytest 6814/0, trio
43/12/11 bitwise no re-record, tests/cpp + ASAN clean, --as-cran 1
expected NOTE.

Doors left: the fit-time test-basis channel (VD-held modelling
decision); weights never ride saved state (logistic); per-forest
off-sample replay on multinomial and OP_GROW outside the multi-forest
fuzz mask (fuzz-arm doors); docs/plans/bartcore-review-tour.md is
refreshed (anchored at ae5b91d8) and tracked - read it first.
Consumer facts: stan4bart is compatible with numSavedSamples' new
meaning (resets storage per iteration, reads the count after its
loop); bartCause's `bcf(keepTrees=TRUE)` was rotated before the fix
and is cured by it, no bartCause change; its `predict.bartBCF` error
text is now stale too - it claims dbarts lacks per-forest replay,
shipped at 63df524e - a bartCause-side follow-up.

### The 0.9-34 CRAN patch lands on main; bartcore rebases onto it (main b9d42948; bartcore d3acbdc1 + 0b424eeb, 2026-08-20)

CRAN broke 0.9-33 on six flavors: R-devel's RNGkind gained a
binom.kind argument, so the deparse-based kind extraction in
A_class.R's validity method left invalid R at install time, and Apple
clang 14 on Intel macOS accepts -std=gnu++20 while its libc++ lacks
std::bit_cast. Main took PR #81 (Brian King) by rebase-merge plus
five follow-ups: dead bit_cast-shim removal, release prep (Date,
NEWS), a per-field .Random.seed validity check in
external/randomBase.c replacing the seed0 > 11000 bound (R-devel
encodes binom.kind at the 100000s digit, pushing the default first
element to ~110403, so every non-native RNG creation warned and
clock-seeded), and an install-time hardening of the validity
fallback (it read .Random.seed before it exists and restored the
seed into a local). Gated: two independent --as-cran runs clean,
tinytest 1714/0, lintr delta zero, and an R-devel container leg
(baseline reproduces both CRAN failures; fixed tree installs;
1716/1716; guard discrimination shown in both directions), then all
five main CI legs green (one r-hub setup-r infrastructure rerun).
dbarts_0.9-34.tar.gz staged for submission. stan4bart needs no
release: all six of its CRAN errors trace to dbarts, and its CXX17
build routes around the bit_cast landmine.

bartcore then rebased onto the new main to restore
fast-forwardability: 1178 commits replayed, 8 conflict stops all
resolving to the bartcore side, author/subject md5 identical across
the range; range-diff shows 10 commits with modified diffs (the 8
stops plus two context shifts). d3acbdc1 ports what the rewrite
would otherwise lose: the last std::bit_cast in the registration
tables becomes a (DL_FUNC) cast (include <bit> dropped from the
bridge; the engine's remaining <bit> uses are popcount/countr_zero,
not affected), and the per-field seed validation lands against
bartcore's own trimmed sentinels. 0b424eeb splits the released
0.9-34 section back out of 1.0-0's NEWS (the four fix items plus
verbatim-duplicated feature items move to history; 1.0-0's section
now reads relative to 0.9-34). Full battery on the rebased tip:
tinytest 6463/0 locally and 6480/0 in the R-devel container (which
also shows unseeded nchain=2/nthread=2 runs bit-identical - the
0.9-x ambient-seed race is closed by the rewrite, as the UPGRADING
notes claim), tests/cpp full and sampler-filtered green, trio
bitwise 42/12/11 identical-draw counts with no max |z| lines,
--as-cran Status OK. HISTORY NOTE: every bartcore hash cited in
docs/plans and docs/design below this note predates the rebase; the
old history remains reachable at tag bartcore-pre-cran-rebase (old
tip fc2af5ce).

Residue, recorded not fixed: ext_rng_createDefault decodes
seed0 % 100 against R's kind numbering while bartcore's trimmed
ext_rng_algorithm_t renumbers (MERSENNE_TWISTER = 0); the path is
dead in-tree (the only caller passes useNative = true) and the new
guard now rejects rather than mis-mapping, so it waits for the next
randomBase.c edit.

Addendum, same day: VD redirected on rchk - CRAN's additional-issues
report gets preempted, not risked at the incoming inspection. Main
gained three more commits (tip cb290550): a defensive-PROTECT pass
over every function CRAN's rchk flags (REPROTECT_SLOT idiom in
R_interface_common.cpp mirroring bartcore's own earlier pass;
createMatrix, rbart_getFitted, the two randomBase.c getDefault
readers) plus a follow-up making the crossvalidation UNPROTECT
counts constant (rchk abandons functions with computed counts) and
protecting the custom-loss closure temporaries. Verified with the
kalibera/rchk container itself: baseline 48 [UP] findings, fixed
tree ZERO real findings (two "address taken ... results will be
incomplete" analyzer coverage notes on fitFinalizer/isValidPointer
remain by design - they are unfixable without re-keying the pointer
registry and CRAN's own rchk revision does not emit them). Bitwise
A/B across 27 draw channels identical; tinytest 1714/0; --as-cran
clean; all five main CI legs green. The submission tarball is
rebuilt from cb290550. bartcore re-rebased on top (1181 commits, 5
stops, log md5 preserved): the port audit found bartcore ALREADY
carried every protection (its pass landed first; main's is the
mirror), so the only tree delta vs the prior bartcore tip is the
one NEWS item. Hash map for this note's citations after the
replay: the port commit d3acbdc1 -> e75a2c79, the NEWS split
0b424eeb -> f04e8686; the records commit rides at a4e6033d. The
pre-rebase tag still anchors the original line.

### Cleanup wave 2: the amplitude family sheds its BCF spelling (8215a53c + be1d06e4, 2026-08-20)

The bcf-naming-generalization item executes under the full arc
discipline: Opus design memo (census, vocabulary, state-key
migration, baseline-impact analysis) ->
independent blind critique (verdict EXECUTE WITH AMENDMENTS; it
refuted the memo's error-text census as 4x undercounted and its
one-assertion version-bump claim as three, amended the refusal
message's wording off a factually-false forests = phrasing, and
CONFIRMED everything load-bearing: the bump is genuinely required
because a renamed OPTIONAL block would silently default, structSize
does not move, no baseline re-record, no consumer reach across all
four sister packages, "glue" wins the key name on shipped-vocabulary
precedent) -> two-commit execution -> independent gate-run.

8215a53c: 20 identifiers renamed across engine/bridge/R/tests
(AmplitudeSpec, ForestStructureSpec, AmplitudeState,
AmplitudeForestCombiner, createAmplitudeSampler, hasAmplitudes,
samplerCarriesAmplitudes, refuseAmplitudeMutation,
refuseUndefinedTestFits and kin), ~21 error texts rewritten to the
amplitude vocabulary, state key "bcf" -> "glue" with
stateFormatVersion AND minReadableStateFormatVersion 1 -> 2, control
attribute bartcore.forests. The v1 refusal was verified LIVE twice
(implementer and gate-runner, the latter through the public
forests = route): "state encoding version 1 ... predates the oldest
this dbarts (2) can read", fired before any block read. Three named
assertion moves in test-sampler-state-format.R plus one added
v1-refusal test. Zero old-name residue (three independent greps).
The KEEP class held: bcfGlue and bartcoreBCFSampler (semantics ARE
two-forest BCF), bartcore_createBCF/createBCFHolder, test file
names, harness files and baselines - NO baseline re-record, trio
bitwise 42/12/11 on both independent batteries, tinytest 6463/0,
R CMD check --as-cran OK twice.

be1d06e4: docs/design + architecture re-pointed (11 files);
multiplier-combiner.md's debt paragraph records the discharge AND
corrects its own two mispriced claims (structSize; the M4.3
re-encode analogy); check-doc-freshness green (503 anchors, 35
symbols). Residue, recorded: two never-emitted static_assert
strings still say "BCF is a constant-leaf model" (combiner.hpp,
chain.hpp) - swept opportunistically with the next edit there.

### The oracle lane closes: hazard and hurdle exact-posterior gates (0787b8e3, 2026-08-20)

The fork taken as recommended under VD's discretion grant.
benchmarks/R/hazard-exact.R and hurdle-exact.R gate the two families
whose only prior evidence was bitwise reduction - blind by
construction to transformation defects - against exact posteriors
derived independently of the package's expander and splitter, on the
original data scale.

hazard: 2 covariate cells x 3 periods, single tree, tree space
enumerated exactly (62 trees, prior mass sums to 1); the reference
collapses the subject-level discrete-hazard likelihood to per-cell
Bernoulli factors read straight off (x, time, status) - no
person-period expansion anywhere in the reference. PASS, max gaps
0.0008/0.0005 vs tol 0.004. hurdle: split re-derived from y alone,
occupancy by quadrature, positive factor conjugate on the engine's
internal rescaling, all three channels gated. PASS,
0.0002/0.0000/0.0002. Teeth measured on BOTH axes per family -
reference poisons AND transformation-semantics poisons (the
wrong-censoring-convention risk set fires at 0.69; the y > 1 split
fires at 0.20-0.32) - the second axis is exactly what the reduction
gates cannot see. INDEPENDENT ADJUDICATION (own re-derivations to 6
dp by different methods, own poisons, tolerance audit: margins
2.5-8x over measured seed spread, no detectable bias at a 12-seed
mean) CONFIRMED both derivations; its three comment/format amends
are applied, including honest-scope headers (the hazard survival
channel's moment refinement and hurdle's structure-prior terms are
exercised but not discriminated at these designs - stated, not
implied away). Wired into exact-gates.yaml's exact-posterior group;
first live CI run GREEN on this commit (all six workflows).
Design finding recorded in the script headers: a constant filler
column is NOT inert on this engine - a non-quantile grid gives a
degenerate column its full cut complement and availability is
index-based, so it shifts the tree prior; both gates keep every
column live. Every shipped family now carries SBC or a gated
posterior oracle, or both.

### Cleanup wave 1: perfect-fit muffle, filtered-run fixture, feature-matrix prose (6d82ac86 + e74380f7 + b1559457, 2026-08-19/20)

The RC soak is VD-paused for a CRAN patch on main; the window goes to
cleanup at orchestrator discretion, VD's grant covering the two open
forks.

b1559457 (the muffle fork, taken as recommended):
estimateSigmaFromLinearModel muffles base R's "essentially perfect
fit" condition at the source - only that condition, everything else
passes through - and test-boundary-inputs.R's constant-response pin
re-tightens to an exact count of one. Container-verified on the
configuration that warns (amd64 R 4.6.1 + OpenBLAS: 2 warnings
before, exactly 1 after; the fixed test file 13/13 there). CI's
ubuntu legs green on the landed commit confirm cross-platform.

e74380f7 (TODO test-bartcore-filtered-run-false-fail, discharged):
testBCFGrowForestFromRoot builds its fixture from a local xorshift
stream instead of the shared runif01() state, so its hardcoded
characteristic value (re-derived, -0.028618738206336595) is
bit-identical in filtered and full runs; `./test_bartcore sampler`
alone now passes. Zero downstream pinned values moved - the old code
already restored the shared state it consumed, so isolation shifts
nothing behind it.

6d82ac86 (TODO feature-matrix-prose-staleness, discharged): the
three verified findings fixed in place - test-capi.R's BCF reach
correctly attributed to forests = list(forest(basis = ...)); [f23]
quotes the multi-forest calibration-map string its cited test
actually pins (the "two-forest" literal survives at the separate
direct-engine-gate site by design); both historical back-references
re-derived to live locations.

Batch mechanics: feature-matrix landed alone (docs-only, fires
nothing); the muffle and fixture slices landed under the batch
clause with one merged-tree battery (tests/cpp full AND filtered
modes, tinytest 6462/0, trio bitwise 42/12/11, R CMD check --as-cran
Status OK), CI six-green on b1559457.

### Soak wave 3 addendum: the anchor sweep (662a7ef8, 2026-08-19)

Ledger item 15 discharges, run LAST after every content edit per the
serialization rule. The five known-stale sites re-derived by content
(data-store.md's transaction anchors, grow-from-root-default.md's
growFromRoot, model-space-survey.md's setData,
within-chain-threading.md's fan-out half, data-ownership.md's
citation of deleted const_cast code - rewritten to cite the file
without a line, the doc's own convention for deleted content). Drift
sweep over every file this session's landings moved: ~150 anchor
instances across 48 design docs checked, ~100 re-derived, ~50
confirmed unmoved; deltas non-uniform throughout, one content
move across files (data.R -> R/utility.R), and error-style.md's
"Violation" examples rewritten where the violations had since been
fixed in code. The sweeper self-caught three of its own
manual-counting errors on range endpoints; orchestrator audit
re-opened a sample and every sampled anchor landed on its named
content. Docs-only; fires no CI; the built package is byte-identical
to 49b8cac3's.

### Soak wave 3: cross-arch dead kernel families + error-style residue (ded1f54e + 49b8cac3, 2026-08-19)

Two slices land under the batch clause (each independently
implemented and gate-run; one merged-tree battery on the stack -
tests/cpp clean, tinytest 6463/0, trio bitwise 42/12/11, R CMD check
Status OK, air clean; rebase onto the records tip verified
patch-id-identical to the gated commits).

ded1f54e (ledger item 16, opened and closed same-day out of item 8's
caller map): the callerless kernel families are gone on every arch -
misc_addVectors, misc_addVectorsWithMultiplier, the two
*AlignedVectorsInPlace dispatch generics (whose per-arch entries were
already mere rebinds), misc_transposeMatrix with its _c/_sse2/_avx
implementations and static block-transpose helpers, and
misc_simd_alignment (its only readers were the aligned bindings; its
documented consumer, misc/memalign.h, was deleted in an earlier
pruning). misc_subtractVectors stays - the engine's response rescale
calls it. 314 deleted / 12 added; docs/design/kernel-vocabulary.md
corrected to match. x86 verification: both x86 files syntax-compiled
via clang --target=x86_64-apple-darwin against a doctored config
defining the x86 ISA macros (implementer and gate-runner each
re-derived this independently); full x86 compilation rides CI.

49b8cac3 (ledger item 7): the error-style residue clears. R6:
spec.R's probit-weights refusal collapses to refusal plus one remedy.
R14: setControl's per-slot refusal takes the rule's verb shape. R11:
the seven competing length-mismatch templates unify onto the rule's
own settled shape across 13 sites (augmentation.R's sprintf form
stays as the rule's named exception), with 9 message pins updated;
error-style.md's corpus and appendix counts re-derived against the
live tree (stop 549, Rf_error 348, ext_throwError 2 - the
gate-runner's independent recount matched exactly). The two stale
comments the M0 verification flagged are fixed (R/bartcore.R's
predictor-session comment; refuseMultiForestMutation's docstring no
longer lists the response conduits its sibling governs). Text-only
proven twice by parse-tree comparison with string literals
normalized (implementer, then gate-runner: every divergence confined
to string literals and comments).

### Soak wave 2: count-exact warning assertions + the M0 component contract (c4971bbb + cf4a290c, 2026-08-19)

Two slices land under the batch clause (each implemented and
independently gate-run on 3e9ef3aa; one merged-tree battery on
cf4a290c - tinytest 6462/0, trio bitwise 42/12/11, R CMD check
Status OK with the new vignette chunk building, air clean - gated
the push).

c4971bbb (ledger item 10): all 17 expect_warning sites converted to
count-exact assertions. tinytest's expect_warning(expr, pattern) was
verified to pass with an extra unrelated warning (either order), a
doubled match, or an interleaved message; no site asserted a count.
inst/common/captureWarnings.R (withCallingHandlers, muffle+record,
returns every warning condition in emission order) now backs each
site: expect_equal on the exact count plus expect_match on
conditionMessage (and class where checked), all unqualified.
Discrimination proven three ways: synthetic probes (old form TRUE /
new form counts 2 under injected extras), the implementer's two
real-site injections (test-bcf-creation, test-multipleAssignment:
"Expected '1', got '2'" where the old assertion still passed), and
the gate-runner's own re-plant at a THIRD site (test-rbart-weights
via R/data.R's train-only-weights warning). +24 assertions
(6438 -> 6462).

cf4a290c (M0, closing the multiforest-extension-surface arc): the
on-ramp contract ships. docs/design/bart-as-a-component.md is the
internal layer (driver-loop bitwise identity, the multi-forest
mutation-predicate matrix entry by entry, what state a mutation does
not carry vs who reinstalls it, the two sweep-boundary hooks, flat
graduation debt, a Measured section);
vignettes/gibbs_sampler_mixture_model.Rmd gains the K-forest
multiplier-composition section (runnable backfitting recipe with the
exact-zero snap, cost with its denominator, five honest losses, the
graduation remedies). EVERY measured claim re-derived on this tip -
and three of the older notes' figures did NOT survive: the R
composition tax is 1.9-2.1x at K = 2 rising roughly linearly to
4.6-5.2x at K = 8 against one matched-tree batched sampler (not
5.1-5.4x flat in K); the unremedied two-route CATE divergence is
11.9 percent pointwise (not ~9), falling to 3.4 percent pointwise /
rms-inside-Monte-Carlo-error with the remedies; the excluding-weight
rebuild divergence is 0.057 against a fit range of 0.035. The doc
ships the measured values with methodology. Anchor drift recorded as
it now stands: the transactional multi-forest guard is retired (the
predictor surface is open at every entry BY DESIGN),
Chain::supportsResponseMutation carries no family conjunct, and z
rides data@bases - there is no treatment slot. Bitwise driver-loop
identity was verified twice (implementer, then the gate-runner's own
independent probe, all identical() TRUE with the n.thin-mismatch
negative half discriminating). Residue for the error-style slice:
two stale comments found in existing files (R/bartcore.R:167;
src/R_interface_bartcore_common.hpp:105-108 predicate listing).

CI FIRE AND CURE: cf4a290c went red on the three ubuntu R-CMD-check
jobs and both sanitizer jobs - one site, test-boundary-inputs.R's
constant-response pin, "Expected '1', got '2'" - while macOS and
both Windows jobs passed. Diagnosis (container-reproduced on amd64
R 4.6.1 + OpenBLAS, the CI configuration; native arm64 Debian
r-base does NOT reproduce): the starting-sigma lm on a constant
response leaves residual variance within a factor of two of
summary.lm's 1e-30 essentially-perfect-fit threshold, and which
side the least-squares solve lands on is a BLAS/architecture
property - OpenBLAS falls below (second, leaked base-R warning),
Accelerate and reference BLAS stay above. The count-exact pin did
its job: the old pattern-only assertion was blind to this leak on
every platform. CURE 07b02f06: the site partitions by pattern -
exactly one degenerate-response warning, at most one perfect-fit
warning, total must equal the two matched counts - verified on
macOS (6463/0) and in the amd64 container pre/post. cf4a290c's two
red workflows stay red in history by design; 07b02f06 carries the
green. OPEN QUESTION for VD (not taken here): whether
estimateSigmaFromLinearModel should muffle base R's leaked
"essentially perfect fit" condition at the source - it is an
internal implementation detail whose text means nothing to a
dbarts user, and the degenerate-response signal is already carried
by dbarts's own warning; muffling would restore an exact count on
every platform.

### Soak wave 1: NEON dead kernels + fuzz reach; bairrtt compatible; oracle-lane scoping (4a6acdbc + 15e1174a, 2026-08-19)

Two ledger slices land under the batch clause (each implemented and
independently gate-run on 3e9ef3aa; one merged-tree battery on
15e1174a - tests/cpp clean full run, tinytest 6438/0, trio bitwise
42/12/11 - gated the push), plus two report-only outcomes.

4a6acdbc (ledger item 8): the seven unreachable NEON kernels are
gone. Four (addVectors, subtractVectors, addVectorsWithMultiplier,
addAlignedVectorsInPlaceWithMultiplier _neon) had no dispatch entry -
their generics are plain scalar functions, so nothing could reach
them; three (addAlignedVectorsInPlace, subtractAlignedVectorsInPlace,
transposeMatrix _neon) were bound in misc_simd_setSIMDInstructionSet
but their generics have zero callers anywhere in the tree. The NEON
branch now binds the aligned generics to the unaligned NEON kernels -
exactly what the AVX and SSE2 branches already do - and
transposeMatrix to the scalar kernel. 418 deleted / 3 added; both
batteries (implementer, then independent gate-runner from a fresh
--preclean build) fully green with the trio bitwise. FOLLOW-UP
(ledger item 16): the *AlignedVectorsInPlace and transposeMatrix
families are provably callerless on EVERY arch - generic, sse2, avx
and _c variants alike, plus misc_simd_alignment whose only consumers
were the aligned variants - a clean cross-arch dead-code slice when
taken.

15e1174a (ledger item 9): the whole-variance-forest exactness
assertion in tests/cpp/test_fuzz.cpp (formerly never executed - 0
fires at 3 and at 40 seeds) is reachable: a fourth heteroscedastic
shape confines the variance forest by column mask to columns 0-1, so
columns 2-3 are untouched by construction, and the site now fires 4
times (2 columns x 2 chains) at the default entry. A permanent
varianceForestsSkipped > 0 gate keeps the branch from going dark
silently. Teeth proven twice (implementer, then gate-runner's own
probe): inverting the comparison fails exactly 4 times at that check
string and nowhere else. Existing shapes 0-2 keep their seeds and
streams; fuzz-suite runtime delta below noise (0.281 -> 0.276 s/run).

bairrtt compat check (report-only, no commit): bairrtt main @ 6167423
against tip 3e9ef3aa in a private lib - dbarts and bairrtt both
install clean, tinytest 206/206 across all six files, every consumed
surface exercised live (dbarts(), dbartsControl, setCutPoints,
sampleTreesFromPrior, run, predict, setPredictor(forceUpdate=TRUE),
updatePredictorPerObservationJointly across the response/assignment
pair, extract). VERDICT compatible; no migration cost, nothing
resembling a dbarts regression. The coordinated-merge gate for
bairrtt is clear.

Oracle-lane scoping (report-only): of the five families the
SBC-deepening item targeted (aft, hazard, hurdle, heteroscedastic,
monotone), three - aft, heteroscedastic, monotone - have carried
gated exact-posterior oracles in exact-gates.yaml since e22ed536
(2026-07-24, predating the item's spec), and SBC is structurally
ill-posed for hazard/hurdle as designed (person-period/positive-
subset DGP depends on y, breaking exchangeability;
sbc-family-tiers.md). Their reduction gates prove family ==
target-family-on-expanded-data bitwise but are blind to expander
defects by construction. The lane's real remainder is a hazard and a
hurdle posterior oracle; the fork (accept reduction-plus-target-
oracle composition vs build small exact-conditional oracles on the
original data scale, recommendation build) is WITH VD. No SBC runs
are owed; nbinom's deepening was discharged 2026-08-18.

### rc-gate (e) valgrind: BCF spec leak, model-matrix OOB read, and the clean full-suite pass (7e42dc93 + 7be7a126, 2026-08-19)

The first full-suite local memcheck (native arm64 ubuntu:24.04
container, stock R 4.3.3, VALGRIND_OPTS matching valgrind.yaml - the
r-hub amd64 image crashes under Rosetta emulation on this host)
against db5f88ae found exactly two defects: zero invalid reads/writes,
zero uninitialised uses.

Defect 1, caught by memcheck directly: a definite leak, 2464 bytes in
7 blocks, at applyBCFSpec's forests.assign - bartcore::BCFSpec, a
stack local inside createHolder's unwindProtect closure BODY; cleanup
frees the closure's CAPTURES only, and an Rf_error longjmp skips body
locals, leaking the spec's vector once per bridge refusal raised after
the assign. Census matched: 7 expected-refusal BCF creations == 7
blocks (R refuses most bad compositions before C; only
hand-built-object backstops reach the hole). Fix (7e42dc93): the spec
joins the capture list beside the storage it borrows, both creation
routes (createBCFHolder preventive, its family door fires before the
assign). Targeted proof: a 3-refusal repro leaks 1056 bytes/3 blocks
unfixed, 0/0 fixed. Gates: tinytest 6427/0, tests/cpp full pass,
equivalence bitwise 42/12/11, CI green all six; no NEWS entry (pattern
new in 1.0-0, never shipped).

Defect 2, not caught by memcheck (reports nothing - R carves nodes
from its own arena) but by valgrind CHANGING THE HEAP: under
zero-filled pages the extra-factor-level expect_error in
test-bart-bart2.R stopped firing, halting the suite. Investigation found a live
out-of-bounds heap read on every R version since 802daf36 (2015),
pre-bartcore: on the indicators route (bart() x/y, bart2
factors=indicators, predict) the replay indexes TRAINING-sized
per-level counts by the TEST column's level count at four sites in
makeModelMatrixFromDataFrame.c; the documented extra-level refusal
fired only when the garbage read was positive. A masked accept fits
silently wrong (mean |yhat-y| 1.60 on phantom-level rows vs 0.435
elsewhere; unbounded - 5003 test levels emitted 4623 garbage columns),
the flip pure heap state even on R 4.6.1 (repeated calls gave ncol
4,4,4,4,6,6,6,6); categorical route unaffected (real level check),
fewer-levels factors index in bounds.

Fix (7be7a126), two layers since LinkingTo consumers reach the .Call
entry directly: refuseWiderTestColumns in validateXTest (factor,
character, matrix columns) BEFORE the replay, plus a shared
validateDropPatternLength at all four C read sites. Both refuse the
WIDER direction only - a training level declared but never observed
contributes no column, so a test frame lacking it still aligns and
predicts correctly (verified pre-fix); a != guard would have broken
that. Extra-level pin kept and tightened (distinct messages,
extra-level vs fewer-levels); +11 tests including a 5000-extra-level
arm and direct-.Call C-guard arms; mutation proofs per layer (R guard
dropped: C still refuses every route; C dropped: R refuses every
R-level route, direct .Call accepts - the hole C closes). Gates:
tinytest 6438/0, tests/cpp full pass, equivalence bitwise 42/12/11,
air/lintr clean, R CMD check --as-cran OK, NEWS parses at 281 (shipped
bug, gets the entry), CI green all six incl. Windows; tests/cpp does
not compile makeModelMatrixFromDataFrame.c, so local ASAN carries no
signal here by construction.

Final full-suite memcheck against 7be7a126: COMPLETE - "All ok, 6419
results (1h 47m)" (platform-skips account for the count vs 6438 on the
dev host), ERROR SUMMARY 0 errors, all lost bytes 0, zero invalid or
uninitialised accesses, applyBCFSpec absent.
rc-gate (e) is discharged in full (rchk clean at the merged tip per
the prior note; valgrind clean full-suite at the tip), and with it the
whole rc-gate slate (a)-(e). RC declaration remains VD-held.

### rchk gate: triage and defensive PROTECTs (69de27ac, 2026-08-19)

rc-gate (e)'s rchk half ran locally (kalibera/rchk docker under amd64
emulation, ~5 min/run; the CI workflow cannot register until the
branch reaches main). Initial verdict FAIL per the workflow's own gate
grep: 11 [UP], 0 [PB], clean maacheck, 2 analysis bailouts on
bartcore_bridge::setState.

Independent triage found ALL 11 false positives, each with a
mechanism: ordinalExpr/dimsExpr are attributes of rooted .Call
arguments (rchk classifies Rf_isInteger allocating because R's
isInteger calls inherits()); xExpr is an attribute of the
R_PreserveObject'd data spec, the rooting dbarts.h ships as the flat-C
contract; the seven installChannel flags miss that SET_VECTOR_ELT into
the PROTECTed result runs BEFORE Rf_mkChar - rchk does not model
rooting by container assignment. The setState bailouts are analyzer
capacity on a pure reader off the rooted stateExpr (balanced PROTECT
pairs, refusals defer through errorMessage so no longjmp crosses
them). All four sites confirmed empirically under gctorture(TRUE)
including the flat C API legs (capi_get_trees, capi_set_state restored
prediction bit-identical).

DECISION (orchestrator): apply a 10-line defensive PROTECT patch
anyway - CRAN publishes rchk output and standing false positives would
blind this gate's next real finding; the sites are one-time bridge
paths, so the sampling loop pays nothing. Landed as 69de27ac (+19/-3,
each hunk carrying a constraint comment naming the real rooting).

Gates on its own base: tinytest 6413/0, tests/cpp full pass,
equivalence trio bitwise 42/12/11 identical-draws lines. Merged-tree
gate stacked on 8dbc0ce9: tinytest 6427/0, tests/cpp 252 ok, main
equivalence 42/42, and a fresh rchk over the merged tree CLEAN -
"Analyzed 13502 functions, traversed 3049440 states", 0 [UP], 0 [PB],
empty maacheck, only the 5 expected bailouts.

GATE LESSON, durable: an empty [UP] grep PASSES on an OOM-killed
analyzer (observed: "Killed" mid-analysis left an empty .bcheck) - the
gate must also require the "Analyzed N functions" completion line and
the absence of "Killed".

The valgrind half of (e) also ran (native arm64 ubuntu:24.04 container
fallback, stock R 4.3.3, the r-hub amd64 image crashes under
emulation; VALGRIND_OPTS matching the workflow): zero invalid
reads/writes, zero uninitialised uses, but ONE definite leak (2464
bytes in 7 blocks) allocated at applyBCFSpec
(R_interface_bartcore.cpp:2173, spec.forests.assign) via
createHolder's unwindProtect lambda on the BCF creation path - fix
slice in flight; and on that R 4.3.3 the test-bart-bart2.R
extra-factor-level expect_error did not fire, halting tinytest partway
(in-contract divergence, DESCRIPTION floor is R >= 4.0.0;
investigation in flight; CI oldrel-1 R 4.4 is green). Full-suite
valgrind re-run on the fixed tip remains the gate.

### rc-gate (d): slot-sourced heteroscedastic warm starts pair the sample's own scale surface (8dbc0ce9, 2026-08-19)

VD decision 2026-08-19 fix-now (decision brief verified in advance
that no equivalence scenario or RNG-pinned test reaches hetero
slot-sourced installTrees, so the draw-moving-in-principle fix
re-records nothing). Design memo -> independent blind critique verdict
SOUND-WITH-AMENDMENTS -> implementation with amendments binding.

A1: the guard is stride equality - savedVarianceTrees.size() must
equal capacity * nvt (capacity from the mean side's saved/live ratio,
the same quantity that bounds slot) plus the base bound; a size-only
guard would let a state whose live variance block was replaced by a
shorter one slice ACROSS slot boundaries and install a cross-slot
mixture presenting the right count.

A2: a saved slot is exempt from the state-side occupancy pass but
becomes occupancy-checked where it turns live, so a donor that kept
sweeps then moved its rows can hold a slot the destination rightly
refuses - heteroscedastic.md's two contrary statements rewritten, the
varianceMismatch message extended with the scale-leaf positivity
clause.

A3's premise was CORRECTED in implementation: the enum value landed on
an existing line so no wholesale anchor shift occurred; moved anchors
re-derived by content with distinct deltas and one deliberately
unmoved.

Mechanics: readWarmStartState now parses the variance.saved.* blocks
(optional - adjudication belongs to installForests, which knows the
slot); slot >= 0 slices savedVarianceTrees and the pooled-categorical
mask channel at slot*nvt; the sliced trees are held to the scale-leaf
positivity law (hand-buildable surface multiplied by capacity, and a
rebuild scatters a leaf straight into a divisor); new
WarmStartResult::varianceSlotMismatch with its own bridge error.

Tests: new inst/tinytest/test-heteroscedastic-warm-start.R (s-surface
state pin at slot k with expected values computed R-side from the
donor's saved blocks, pooled-mask arm at slot 0 where the expectation
is a buffer prefix, refusal arm at k = capacity-1 - the only slot a
one-short buffer trips, live-sourced regression arm, occupancy-edge
arm), plus testVarianceWarmStartSlot in tests/cpp/test_state.cpp with
the stride arm. Five mutations all killed: live-copy revert fails 5 R
assertions, parse-discard refuses the positive arm, constant-slot-0
fails 2 C++ assertions, mask-slice-to-live fails the pooled rebuild,
stride-to-bound-only fails the C++ stride arm (A1's own proof).

Gates: tinytest 6427/0; tests/cpp + ASAN/UBSAN clean; equivalence trio
BITWISE 42/12/11 identical-draws lines - zero re-records, the brief's
prediction confirmed in practice; R CMD check --as-cran from a
clean-copy tarball Status OK; air/lintr clean; NEWS parses at 280.

RESIDUE, verified stale AT BASE and left untouched (the 95f6fd7f
full-namespace resync missed them): data-store.md:198/201/206/214
(constructs sit ~320 lines below the cited sampler.hpp lines),
grow-from-root-default.md:141, model-space-survey.md:71 (quoted text
absent from the cited file), within-chain-threading.md:204's :675-702
half, data-ownership.md:213 (names deleted code) - ledgered below.

### rc-gate wave 1: NEWS binary-k, sigest sparse surface, equal-rank-1 coverage (105f2bd6 + 3d5d2ed5 + b9c3f313, 2026-08-19)

rc-gate items (a), (b) and (c) landed as three commits; (b) and (c)
were built in parallel worktrees off db5f88ae/105f2bd6, stacked by
rebase, and gated by one shared merged-tree battery per the batch
clause.

(a) 105f2bd6: inst/NEWS.Rd entry for 5b6e4825's binary-k default
change (chi(1.5, Inf) -> chi(1.5, 2.0)) in the 1.0-0 NEW FEATURES
subsection; man/dbartsPriors.Rd verified already correct, untouched;
the entry's "prior median 1.9" claim reproduced numerically
(2*sqrt(qchisq(0.5, 1.5)) = 1.906). Gates: NEWS parse 277 entries;
clean-copy tarball check (vignette-flag warnings only, from
--no-build-vignettes).

(b) 3d5d2ed5: sigest's silent sd(y) fallback on sparse-backed
predictors now warns - a classed warning
(dbartsSparseSigmaFallbackWarning, dbartsWarning) at
estimateSigmaFromLinearModel's predictorSourceIsSparse branch,
fallback estimate unchanged. sparseVector/dgCMatrix formula columns
get sparseFactor's guarded refusal via a new shared
refuseSparseFormulaColumns helper (R/mixedMatrix.R), replacing the
sparseFactor-only inline loops in dbartsData's formula path and
extractMultinomialFormulaData (a census-found second duplicate);
refusal message widened to all sparse kinds ("sparse predictors must
be specified through the x/y interface"), sparseVector normalized off
its concrete subclass; man/dbartsData.Rd corrected (had described the
raw S4 error verbatim); two NEWS BUG FIXES items. test-data-mixed.R:
two message-content expect_error arms plus a warning-count pin (1
mixed, 0 dense control); anti-vacuity: 4 new-assertion failures
pre-change, green after. Gates: tinytest 6413/0, equivalence 42/42
with no max |z| line, R CMD check tarball Status OK, air clean, lintr
0 on all six touched files, NEWS parse 278.

(c) b9c3f313: tests/cpp only, zero production change -
testEqualRankOneComparison in test_moves.cpp (206 lines, local ext_rng
so filtered runs stay clean). Fixture: a one-cut binary split column
under an all-zero weight vector puts root, split branch and every leaf
at veto rank 1, so birth/death/change realize ONLY equal-rank-1
comparisons and each acceptance is a hand-derived prior x transition
constant (birth/death alpha 0.25/0.75 bitwise, change alpha 1.0),
invariant to response, sigma, k and leaf model; the vetoed branch
score is pinned 0.0 bitwise under LinearGaussianLeaf whose naive
vetoed-marginal sum is 8.882e-16; law cited to
docs/design/empty-leaf-veto.md. A dev-time counter (removed before
commit) proved execution: 7 equal-rank-1 resolutions (4 birth, 2
death, 1 change), zero unequal-rank realizations. Mutation proofs:
equal-nonzero-ranks-as-unequal fails 9 checks; -HUGE_VAL both sides
fails 10; dropping the vetoed-leaf skip in logLikelihoodForBranch
fails 3, all in the new test. Suite 251 ok (was 250); ASAN/UBSAN
clean. Discharges the rc-gate requirement that the flagship veto fix
carry no untested branch.

Batch gate on the merged tree: tinytest 6413/0, tests/cpp 251 ok, NEWS
parse 279, air clean, conflict-marker sweep clean. CI on 105f2bd6: a
live Ubuntu-mirror incident stalled setup-r on pkgdown, exact-gates
and the three ubuntu R-CMD-check jobs - each failed fast at the named
step's 12-minute bound exactly as the 887787fd hardening intends, and
bounded reruns cured all three (one to two cycles); macOS/Windows legs
stayed green throughout.

Also this session, for continuity: rc-gate (e)'s rchk leg RAN locally
(kalibera/rchk docker under amd64 emulation, ~5 min): verdict FAIL, 11
[UP] unprotected-variable findings all in the .Call bridge (a
bartcore_run result-assembly cluster, parseControl,
predictorsFromDataExpr, bartcore_setCounts) plus two analysis bailouts
on setState; a
triage-and-fix slice is in flight. The valgrind leg and the
VD-approved installTrees hetero pairing fix (decision 2026-08-19: fix
now) are in flight, own notes to come.

### Consumer lockstep rebuilds against the new ABI token (no dbarts commit, 2026-08-19)

The rebuild obligation the ABI-token fold created (ledger items 11
remainder and 12), discharged with ZERO consumer source changes -
the fold was rebuild-only for consumers, as designed. stan4bart
(bartcore @ 54e157b): --preclean rebuild clean; resolved-sigma
obligation VERIFIED SUPPLIED (stan4bart_fit.R sets data.bart@sigma
from the lmer/lm init before create, dbarts's spec layer fills only
when NA so the value survives, and init.cpp calls setSigma after
create); gates green - tinytest 531/0 under NOT_CRAN,
compare-posterior all tiers pass vs its MANIFEST baseline row,
recompute-exactness 3/3 at 7e-14. treatSens (dbarts-1.0 @ 1d7c697):
rebuild clean, its own load-time hash check passes, testthat 186/0,
--as-cran with only machine-flag noise. USER LIBRARY refreshed:
dbarts, stan4bart, treatSens, bartCause all rebuilt 2026-08-19
against the tip; stan4bart fit+predict smoke clean; bartCause's
suite 0 failures - the ledgered stan4bart-ABI-staleness test error
is CURED. One durable finding from the stale-binary probe: the
enforced handshake is NOT retroactive - a consumer binary compiled
BEFORE the enforcement never runs the check (a pre-migration
stan4bart loaded and fit silently, failing only inside reshaped
entries on dbarts's own range checks), so pre-enforcement binaries
must be assumed silently wrong and cured by rebuild; binaries built
at or after the enforcement fail loudly at create.

### Full-namespace anchor resync, pre-RC pass (95f6fd7f, 2026-08-19)

The wave-4 precedent (f953127c) repeated over the FINAL content
tree, closing ledger item 13 and the queue's last pre-RC step.
Whole-file walk of docs/design/feature-matrix.md: 535 anchors
checked (502 base citations + 29 comma-list continuations the
freshness checker's regex drops + 4 bare no-token continuations
found by hand) - 230 moved, 305 unchanged, 0 file-changed, 0
content mismatches; the two historical "the earlier X" citations
left per the file's own convention. Every anchor re-derived BY
CONTENT against the live tree from the 8e1e674c baseline blobs,
never by offset; deltas region-dependent, from +1 (TODO) and +3
(chain.hpp head) through +73..+119 (dbarts.h, the ABI-token
reshape) to +127 (chain.hpp tail - the hetero-state-carry drift
item 13 named) and +175 (test_sampler.cpp). Verification was
contested and therefore strong: concurrent passes cross-audited
the same file and surfaced and fixed real errors in each other
(range endpoints, duplicated values, a missed range), and an
independent sample re-derivation settled the one remaining
disagreement (the buildMultinomialSampler span starts at its doc
comment, 3238-3294, preserving the old span's length and
semantics) and confirmed the +73 CAPI delta against the true
8e1e674c baseline rather than the intervening tree. Gate:
check-doc-freshness OK (505 anchors, 37 symbols, 3 scenario-count
claims) on the slice base and again on the rebased tree. Docs-only
(+250/-197, one file); package byte-identical to ad4a131b's fully
gated build.

### Flat-C ABI token layout fold + enforced handshake (ad4a131b, 2026-08-19)

Closes the flat-C struct-append surface decision (VD: fold the
layouts, 2026-08-19) and the ledger's hash half of items 4/11.
DBARTS_C_API_HASH (now 0x6c9776ae1197e8f5) folds, beyond the
entry-point signatures: each ABI struct's size and per-field
name+offset in pointer units (one macro token feeds both, so name
and offset cannot drift; offsets are compiler-reported, never
hand-kept layout text), the two ABI enums' stringized enumerator
lists (the enum bodies are X-macro generated, so an added enumerator
cannot escape the fold), and the callback's hoisted parameter list.
Integers enter FNV-1a via shift-based byte extraction (endian-proof);
a private signature-only literal stays in C_interface.cpp so a
failed build names WHICH half moved. The stubs now ENFORCE the
handshake on first symbol resolution (Rf_error naming the cure),
replacing the opt-in consumer check; test-capi.R re-pins the token,
demotes the layout-blind literal to the must-not-equal list, and
adds a NEGATIVE test: a consumer compiled with a forced-wrong token
(the header's #ifndef guard is the only door) fails at create.
Design memo critiqued blind, verdict SOUND-WITH-AMENDMENTS; all
five mandatory amendments bound and implemented (field names not
field sizes; X-macro enums; pointer-unit normalization instead of a
64-bit guard; shift-based bytes; enforced + negatively-tested
handshake). Documented residue: same-width in-place field type
swaps only. Gates dual-run (implementer + independent runner):
flip-digit re-bake proof; mutation probes discriminate (struct
append, field swap, enum addition, callback change -> combined
assert only; signature change -> both asserts); tinytest 6409/0
(test-capi.R 254 incl. the negative test); equivalence bitwise
42/42, 12/12, 11/11 with zero max |z| lines; air and lintr clean; R
CMD check from an out-of-tree tarball Status OK. tests/cpp does not
link C_interface.cpp: the CI sanitizers leg on this commit is the
ASAN evidence for the stub path. CONSEQUENCE: both consumer
branches (stan4bart bartcore, treatSens dbarts-1.0) now REQUIRE a
--preclean rebuild at lockstep - the token moved and the stubs
enforce it.

### CI hardening - apt stall diagnosis and bounds (887787fd, 2026-08-19)

VD-prioritized mid-arc. The recurring "cancelled" CI runs were
TIMEOUTS misread as the concurrency quirk: job logs show apt-get
update/install (inside setup-r and setup-r-dependencies) emitting
zero output until the job limit, on three consecutive lint attempts,
two exact-gates attempts, pkgdown and three ubuntu R-CMD-check legs.
apt has no default network timeout, so a stalled connection to the
Ubuntu archive mirror waits forever. Not suite creep, not R-package
cache loss (cache hits resolved in seconds). Fix: an apt.conf.d
drop-in (Acquire::Retries 3, http/https Timeout 15s) written before
any apt call in the four exposed workflows; step-level
timeout-minutes on setup-r (12; 15 on the check matrix) and
setup-r-dependencies (15; 20 on the check matrix); job timeouts
resized from hang-sized to observed duration plus headroom (lint
15->20, exact-gates 25->30, pkgdown 90->30, check-standard 80->30).
No gate's checks changed; package tree byte-identical to 77349d29's
gated build. Validated the same afternoon against a live Canonical
incident: the drop-in skipped the dead azure.archive.ubuntu.com farm
in ~23s of retries, the archive.ubuntu.com fallback then trickled
past the 15s inactivity timeout, and the step bound killed it at
exactly 12 minutes with the step named - a silent 25-90 minute burn
became a fast attributed failure. Known residue: no native step
retry exists (reruns stay manual); setup-r's own R-download curl is
outside the apt conf; actions are pinned to moving major tags;
earlier genuine cancel-at-creation events (zero steps run) remain a
separate unexplained phenomenon. The apt-immune legs (cpp-tests,
sanitizers) stayed green throughout the incident.

### Arc-close comment sweep, residual pass (77349d29, 2026-08-19)

The reading-pass sweep owed at arc close, before the RC call. The
bulk sweep landed at dcc8262e (2026-08-17); this pass covers what
wave 5 added afterward plus what that sweep missed: 8 sites in 6
files (+19/-19, comments only). Every site keeps its rule text and
its docs/design citations; only provenance went: plan slice tags
dropped from two docs/design section citations (the files remain
cited, one codenamed heading unpacked to its descriptive name),
"slice-0" renamed to CSC code-validation in the bridge, the
test-capi polarity comment now pins the shipped polarity without
the landed-plans history, a dangling "Q4" citation into
docs/design/multinomial.md (the label exists only in the plan)
dropped with its self-contained rule kept, and a ticket name
removed from R/spec.R's kforest refusal comment. Tree-wide scans
after the pass: zero docs/plans references and zero slice
codenames remain in shipped files (R/, src/, inst/, man/).
Implementer battery and an INDEPENDENT gate-runner battery both
green: tinytest 6407/0, equivalence bitwise 42/42 + 12/12 + 11/11
with zero max |z| lines, tests/cpp clean, air format clean, lintr
clean on all five touched R files, R CMD check from an
out-of-tree tarball Status OK. The arc-close sweep obligation is
discharged; the feature-matrix anchor resync (content edits now
final) is the remaining pre-RC step.

### Wave-5b veto-law enforcement + main-baseline re-record (d15a2bfb + 9ba57750, 2026-08-19)

Fix-queue item 3, the queue's last slice - WAVE 5B CLOSES. Closes
the veto BLOCKER (empty-leaf weight law unenforced by setWeights/
setActiveRows: a grown forest froze bitwise with NaN ratios). The
design memo's 3-valued branch veto rank ships as specified with all
three critique amendments binding: Tree::leafVetoRank (2 member-
empty / 1 weight-vetoed / 0 admissible, degenerating to the old
count test with null weights), moves.hpp::logLikelihoodForBranch
returning BranchScore with the rank taken over fillBottom and only
rank-0 leaves summed (A3), and resolveVetoRank feeding the moves'
EXISTING single exp expressions - the worse rank substitutes
-HUGE_VAL, so rank-0 pairs reproduce the old operands bit for bit
(A1; no factoring, no re-association). Review verdict LAND with
every claim reproduced: operand provenance verified by code
reading; both mutations re-run (dropping the rank gate fails 9
tests + an assert; collapsing rank 2 into 1 fails 6 incl. the
masked state round trip - rank 2 is load-bearing); RNG discipline
confirmed (rejections consume the same draws); the monotone leaf's
weight-law early return safely dropped (feasibility sentinel kept,
the NaN path now explicitly rejected). The adjudication (A2) was
REPRODUCED with an independent counter build: draw shifts land
exactly where an equal-rank comparison is realized - maskprobit 2
realizations / maskordinal 1, all controls zero despite 54k-79k
vetoed scores each, bcf (masked included) and multinomial zero
across their harnesses - so ONLY the main baseline re-records; bcf
and multinomial stand. The new bd-balance veto arm is the oracle:
target derived by R-side enumeration over the restricted posterior,
never the engine's ratio; pre-fix the stranded start NEVER absorbs
in 20000 sweeps, post-fix it absorbs in 1, zero vetoed-partition
visits in 300000 draws, max |z| = 1.8; default arm unchanged.
Records commit 9ba57750: equivalence-d15a2bfb.rds (verified 42/42
strict-bitwise; old-baseline cross-check reproduces exactly 40/42
plus maskprobit 0.48 / maskordinal 0.65), MANIFEST supersession
(c7546233 historical, draw-shifting note with the realization
counts), equivalence.yaml, feature-matrix, and the follow-ups
ledger reconciled (all seven prior items verified still open;
seven residues appended: NEON dead-code kernels, fuzz-config
reach, count-blind warnings, the flat-C struct-append hash gap,
the stan4bart lockstep obligations, feature-matrix anchor resync,
monotone-leaf double-fill + equal-rank-1 coverage). Gates on the
rebased tree: main 40/42 with exactly the two adjudicated lines,
bcf 12/12, multinomial 11/11, tinytest 6407/0, tests/cpp + ASAN
clean, R CMD check --as-cran OK, check-doc-freshness OK, CI green.

### Wave-5b test-adequacy repairs (03b97db7, 2026-08-18)

Fix-queue item 7, test-only (production code untouched; two hand-
backs below). Closes the vacuity BLOCKER: expectTwinsAgree's
assertions were invisible to tinytest because namespace-qualified
tinytest::expect_* resolves past the runner's locally-masked
capturing wrappers to the raw namespace copy whose return value is
discarded; dropping the prefix registers all 13
(test-mutate-sparse-valued.R 16 -> 29 assertions; the original
sabotage now fails 13, restored green). Review swept the whole tree
for the hazard class: zero other namespace-qualified expect_* calls,
and every inst/tinytest source() passes local=TRUE, so no helper
escapes the masked env. SIMD gate: reach extended 3 -> 7 of 14 NEON
kernels - all seven that are live (review verified the other seven
dead by dispatch-table + caller audit: four have no dispatch entry,
three are bound but never called from src/; ledgered as production
dead code) - and, after the review caught reach-without-gating (the
planted tail off-by-one passed the extended test while 47 other
assertions failed), the fix round routed a genuine nonzero offset
through the kernel into compared channels; the same mutation now
fails test-simd.R directly (1/12), restored green. Flat-C family
coverage: logistic/ordinal/aft create+short-run cells added at the
SHIPPED ABI (dbarts_sampler_create -> resolveFamily via the
compiled test-capi.R consumer; a wrong family token demonstrably
refused), after the review found the first version tested only the
engine factory - that engine-level arm is kept for its independent
value, with rngState save/restore verified correct and matching
house idiom. Minors m2 (declined-row slice now genuinely declines,
anti-vacuity assertion, plus an unname() fix for a real base-R
subsetting artifact it exposed) and m5 (formal-existence checks
before NULL pins) fixed; m1 (fuzz-config reach) and m3 (count-blind
warning assertions, a 17-site tinytest semantic limit) handed back
to the ledger; m4 was the hetero slice's statesAgree extension,
landed there. Suite grows 6028 -> 6045 assertions on its base.
Gates on the rebased tree: trio bitwise 42/12/11, tinytest 6402/0
(test-capi.R 252 with both this slice's cells and the flat-C
guards), tests/cpp green, R CMD check --as-cran (core-limited)
Status OK, both mutation proofs re-verified post-rebase, CI green.

### Wave-5b flat-C guards + header truth (61ca56d2, 2026-08-18)

Fix-queue item 6, draw-neutral, no dbarts.h signature changes (zero
consumer migration; stan4bart/treatSens call none of the touched
entries today). Closes the three flat-C BLOCKERs:
setTestPredictors(NULL) now clears the test offset (the header's
promise made true - the pointer is the whole state); a test-row-count
change under a standing offset is refused with the R bridge's message
at the ONLY entry that can resize test rows, killing the OOB read;
and requireResolvedSigmaEstimate at the three creation paths that
hand sigmaEstimate to the engine (createHolder, createBCFHolder,
createFromHandle) closes the silent NA-sigma NaN-draw path - review
re-derived the family law (only gaussian/T/aft read sigmaEstimate;
the multinomial builder takes none; sigmaIsFixed resolves before the
gate) and confirmed the gating exact. Both header-truth MAJORs
rewritten: the ownership preamble now names the retaining setters
(setTestOffset included, which the original audit missed) with a
per-setter ownership sentence on every raw-array setter, and the
state-continuation promise is semantic-not-bitwise, matching
stateContinuation.R. Minors 6-12 all fixed (8 in code: the counts
rule scoped to logistic, the only law reading a weight; 11 verified
field-by-field and by -Wall -Wextra consumer compiles); INFO 13
(the version handshake cannot see struct appends) handed to the
ledger - the remedy re-bakes the API-hash token and is a surface
decision. Review verdict LAND-WITH-FIXES; the fix round corrected
four untrue text claims (ownership "each doc states it", the
every-sweep read claim - gaussian/aft bake buffers at set time, the
sigma rule restated family-keyed, a wrong refusePinnedSigmaChange
cross-reference) and put the BCF family refusal ahead of the sigma
gate. Mutation proofs on every guard (review re-ran two: the
row-count widening and the sigma gate both bite). Migration note:
a gaussian/aft spec with NA @sigma now hard-errors where it silently
drew NaN; in-tree R paths all resolve sigma; stan4bart's create path
owes a verification at its lockstep rebuild (ledgered). Sanitizer
coverage note: tests/cpp never links C_interface.cpp, so the local
ASAN leg cannot reach the guards; the CI sanitizers workflow runs
the full suite (test-capi.R floor 150) against the instrumented
package and covers them. Gates on the rebased tree: trio bitwise
42/12/11 (post-re-record baselines), tinytest 6369/0, ASAN leg
clean, R CMD check --as-cran OK, CI green.

### Wave-5b R-surface fixes + main-baseline re-record (c7546233 + 178d1491, 2026-08-18)

Fix-queue item 5. Closes the summary()/as_draws BLOCKER (silently
wrong on combined multi-chain fits for non-scalar vars) via the
diagnostics reshape, plus the Fork-1 chain-major collapse sweep
(combined multi-chain scalar fields now collapse chain-major to
match yhat.train/varcount/ranef; rbart.R verified to need no change -
all packaging routes through convertSamplesFromDbartsToBart and the
unmeasured-ranef draw uses per-chain raw tau; a naive dnorm loglik
now matches extract("loglik") exactly where the old order was 2.72
off), matchedCall[["prior"]] partial-match, keepsampler||keeptrees,
named factor-level mismatch, truthful missing message, and the
subset misreport. Review verdict LAND-WITH-FIXES; the fix round
closed the review's own regression find (scalarFields omitted
first.sigma/first.k/first.tau/resid.df - mis-split chain margins on
multi-chain fits; failed 4/4 pre-fix), corrected both Rd row-order
formulas to kept-draws-per-chain, made combineChains' 2-D branch
invert the new vector branch with a round-trip test, added a
chain-major order pin that fails on the interleaved base, and made
the factor-mismatch advice truthful. The NaN-sigma minor was handed
back (C-side; discharged by the flat-C guards slice). Every fix
carries a fail-on-pre-fix proof. THE ADJUDICATION: the main harness
fell to 41/42 with `bart2probit 37 summaries, max |z| = 0.00` - the
only scenario pairing adaptive k with pooled multi-chain output.
Independent per-channel verification (20 seeds x 37 channels):
every RNG-stream channel bitwise identical; only k.mean/k.sd moved,
<= 7 ULP, pure accumulation reassociation from the chain-major k
reorder. Decision: draw-neutral RE-RECORD over a permanent
fallback, because bart2probit is the sole tripwire for exactly the
configuration this work touches and |z| < 4 is far too loose there.
Records commit 178d1491: equivalence-c7546233.rds (42 scenarios,
verified 42/42 strict-bitwise; old-baseline cross-check reproduces
the adjudicated profile exactly), MANIFEST supersession row
(4a42620a historical, draw-neutral note), equivalence.yaml,
feature-matrix pointers; check-doc-freshness OK. Gates on the
rebased tree: trio 41/42+adjudicated / 12 / 11, tinytest 6361/0,
R CMD check --as-cran OK, air/lintr clean, CI green. Next owed:
the bartCause lockstep slice (sample-major -> chain-major sigma
reshape, consumer gotcha 2).

### Wave-5b hetero saved-variance state carry (73e14af4, 2026-08-18)

Fix-queue item 2. Closes the hetero BLOCKER: getState/setState kept
the variance forest's live trees but dropped its saved (keepTrees)
buffer, so predict on a sampler re-created from stored state replayed
the identity fill and reported s(x) == 0. The carry now round-trips
the saved variance trees and their mask channels, including the
>63-level mask case (an 80-level fixture stores and round-trips
nonzero live and saved mask bytes bitwise). The critique's A1 landed
STRONGER than specified - a size-equality gate at the restore door
subsumes the presence check and D2 - and its A5 guard, A6 comment,
and D3 identity fill are all in place; a pooled saved tree with no
mask channel is refused at both doors, so the empty-mask assign is
unreachable. statesAgree gained the saved-variance fields and the
extension is proven biting: review re-ran both mutation proofs
independently (dropping the chainFields loop flips the agree pin;
deleting the saved-tree assign fails 10 tinytest assertions on the
1.0 identity fill plus two tests/cpp arms). Review verdict LAND;
its no-op name-presence pin was strengthened to a non-NULL check at
the amend. Recorded residues: the memo's O4 test_fuzz mutator family
for saved variance leaves was dropped (test_state's four refusal
arms cover the ground more directly); O5's pins live in
test-heteroscedastic.R rather than test-sampler-state-format.R with
equal-or-better coverage; feature-matrix CH: line anchors have
drifted and wait for the next full resync (check-doc-freshness
range-checks pass, so its OK is not evidence there). Gates on the
rebased tree (independent lib): trio bitwise 42/12/11 against the
post-re-record baselines, tests/cpp green, tinytest 6338/0, ASAN
clean, R CMD check --as-cran from a clean-copy tarball Status OK,
NEWS parses via the check. CI green.

### Wave-5b multinomial prior fix + baseline re-record (4d9a3337 + 509b68cc, 2026-08-18)

Fix-queue item 4, the queue's only draw-shifting slice. Closes the
multinomial BLOCKER: the level-shift combiner divided the prior sd by
sqrt(m) twice. It now draws from the conditional the leaf prior
states - with per-leaf prior sd scale/k (ConstantGaussianLeaf::
drawFromPrior), shifting forest k's occupied leaves by c/m_k gives
prec = sum_k L_k/(m_k^2 s_k^2); the old form's precision was a factor
m too large (intercept-only case: tau^2/(K m) where the design doc
states N(0, tau^2/K)). docs/design/multinomial.md updated to state
the conditional. Verification at review was independent end to end:
the reviewer re-derived the conditional from the prior, hand-
recomputed the C++ oracle's constants (arm B solves mean and sd out
of the OUTPUT shift at two seeds, sharing no expression with
combiner.hpp), and re-proved discrimination by mutation (re-inserting
the old divisor fails 5 checks with exactly the MANIFEST's numbers).
The permanent prior-predictive arm 7 - the only arm reading the
non-identified level - shows 8.3x separation against a pre-fix
install; the k3countsswap scenario rider is DISCHARGED. Records
commit 509b68cc: multinomial-equivalence-4d9a3337.rds (11 scenarios),
MANIFEST current row naming both P17 oracles with 1027be5 demoted to
historical under a draw-shifting supersession note, equivalence.yaml
baseline line, feature-matrix [f39]=11, TODO trims - nothing else.
The baseline's embedded meta$rev is the pre-rebase sha cd5efd97;
MANIFEST discloses it (the two trees are byte-identical outside TODO
and this plan file). Gates re-run by the reviewer in an independent
lib: multinomial 11/11 identical, gaussian 42/42 --strict-coverage,
bcf 12/12, zero max-z lines anywhere; tests/cpp 246 ok; tinytest
6028 passed, FAILURES 0; multinomial-exact all 7 arms in band. CI
green on both pushes.

### Wave-5b docs-coherence slice (ea630c21, 2026-08-18)

Fix-queue item 8, docs-only. All six MAJOR doc-review findings fixed
plus the minors: makeind.Rd tells the truth about dbartsMixedMatrix
returns; a NEW man/dbartsMixedMatrix.Rd documents the container (its
four S3 methods, supported vs unsupported operations, the as.matrix
escape, the row-not-linear single-index quirk) with its _pkgdown.yml
entry; pdbart.Rd/pd2bart carry the real field lists (bartcall rename,
dropped fields); sparseFactor.Rd states the silent sigest fallback;
dbartsData.Rd scopes the sparse-ordinal formula-path failure; dead
docs/design pointers dropped from shipped help text (docs/ never
installs); getTrees chain-column conditionality; interactions.Rd
wording + bare-vector shorthand; NEWS gains the posterior-integration
entry. Gates: NEWS parses (274 entries), check_pkgdown green,
R CMD check --as-cran from a clean-copy tarball Status OK, tinytest
green, six CI workflows green. Behavioral follow-ups handed to the
ledger (TODO wave5-review-followups): sigest silent fallback should
warn or keep the lm fit; sparseVector/dgCMatrix formula columns need
the sparseFactor-style guard; NEWS owes an entry for 5b6e4825's
binary-k default change.

### Wave-5b bridge-guards slice (613d59b9, 2026-08-18)

Fix-queue item 1. Closes review BLOCKERs 4-5 of 10 (the two bridge
crashes) plus three MAJORs and five MINORs, all draw-neutral. An
ordinal column now carries >= 1 cut: the quantile grid's constant-
column case floors to one degenerate cut at the observed value
(exactly what the uniform grid places), closing both the
printCutoffs unsigned-wrap segfault and the unrestorable-state
$copy() failure; setCutPoints refuses an empty grid; the summary
skips zero-cut columns defensively. The model-matrix builder
validates factor codes: NA codes flow into the missing-data route
(indicator rows go NA, bart() reports the canonical missing-values
refusal, missing="incorporate" models it) and out-of-table codes are
refused by column name - no more wild writes. Warm-start zero-forest
donors refused; assignInPlace requires a non-empty index; the five
minors (unwindProtect on the basis buffer, bounded map negation,
interned attribute symbols, ordered lastIndex test, R_alloc for the
column-type array). Six new/rewritten tinytest files; mutation
proofs on both guard families (planted reverts caught). Gates: full
battery green, ASAN clean, equivalence trio BITWISE (42/12/10
identical-draws lines, zero max-z lines), six CI workflows green
incl. sanitizers. A consequence recorded in
test-sampler-degenerate-cuts.R: zero-cut stores are no longer
R-reachable, so the moves-level degenerate-root guard is now reached
via a single-level factor.


### Wave 5a - adversarial whole-branch review; NO COMMIT (review record, 2026-08-18)

Six independent refute-posture reviewers against 13a8867c (engine
numerics, bridge/memory safety, R surface, flat-C contract, test
adequacy, docs coherence), each with build/run probes in a private
lib from a pristine archive stage. Full evidence in
rc-review-artifacts-2026-08-18/w5{eng,bridge,rsurf,capi,test,docs}-
findings.md; every orchestrator re-verification (nine findings
re-run or code-walked independently) CONFIRMED the reviewer claim -
zero false positives observed.

VERDICT: not RC-clean. Ten BLOCKERs, ~19 MAJORs, ~25 MINORs.
BLOCKERs: hetero predict s(x)==0 after any state round-trip (saved
variance trees not carried by get/setState); empty-leaf veto weight
law unenforced by setWeights/setActiveRows (grown forest freezes
bitwise, NaN ratios); multinomial level-shift sd divided by sqrt(m)
twice (prior defect in the tree-structure marginal - the FIX SHIFTS
DRAWS, so it rides a multinomial baseline re-record);
printCutoffs+useQuantiles zero-cut segfault (bart/bart2/dbarts);
NA-factor wild write aborting R via makeModelMatrixFromDataFrame
and bart(df,y); summary()/as_draws_*() silently wrong on combined
multi-chain fits for non-scalar vars; flat-C setTestPredictors(NULL)
leaves a stale test offset (contra header); flat-C test-row-count
change with installed offset reads heap garbage (R bridge refuses
it, C path does not); flat-path gaussian NA data@sigma yields silent
all-NaN draws (consumer gotcha 4 formalized); the sparse-mutation
pin file's central assertion never registers with tinytest (13
call sites vacuous - proven by surviving sabotage).

Countervailing positive evidence, also part of the record: an
independently coded exact-posterior gate validated the MH tree
kernel to 2.1e-14 with large-sample stationarity at 3.8e-4-7.2e-4
across ordinal/categorical/missing probes; birth/death detailed
balance hand-derived clean; PROTECT balance walked whole-bridge;
the C API signature/hash discipline cannot drift; all 4352 R
expect sites and 52 S3 methods execute; examples and vignettes
reproduce. The core sampler is sound; the defects cluster in
state round-tripping, guard-less boundaries, and reporting layers.

FIX QUEUE (wave 5b), serialized: (1) bridge crash guards
(in flight at this writing); (2) hetero saved-variance state carry;
(3) veto weight-law enforcement at the mutators; (4) multinomial
level-shift fix + baseline re-record (multinomial baseline ONLY,
carrying the k3countsswap rider; the bart2/xbart/mixed-matrix
riders were already discharged at 4a42620a/f009eff8 - this note
originally said otherwise; P17 = the MANIFEST oracle-naming rule); (5) R-surface
fixes (diagnostics layout, prior.scale match.call partial-match,
keeptrees/keepsampler, bart() route gaps); (6) flat-C guards +
header doc truth; (7) test-adequacy repairs (expectTwinsAgree
registration, SIMD gate reach, flat-C family coverage); (8) docs
coherence sweep. Draw-shifting: item 4 ONLY; all others gate
bitwise on the canonical trio.


### Full-namespace anchor resync - WAVE 4 CLOSES (f953127c, 2026-08-18)

The whole-file walk the matrix's value audit recorded as owed: every
anchor in docs/design/feature-matrix.md re-located by symbol against
8e1e674c (the walked tree; the landing rebased onto 142a50b6, which
touches no anchored file, so the stamp stays truthful). MOV, four
byte-identical design docs, and tests/cpp/test_model.cpp carried
over by byte-identity; the bart2.Rd -> man/bart2.Rd alias was added
(man/bart2.Rd postdated the alias table); NO anchor changed file;
two [f15] test_model.cpp anchors stale at c05322a8 itself were
corrected; the two historical numbers stay as written; no cell value
re-adjudicated. Guard exits 0 at the new tip (502 anchors resolved).

Orchestrator verification (the mandated sampling, per the anchor
discipline): per-alias delta histograms show the healthy signature -
RIB 12+ distinct deltas in both signs, data.R correctly unmoved, MOD
uniformly -1 past its single small edit - and ~40 anchors were
opened by content across every alias class. The sampling caught two
defect families the implementer's own two "converged" passes both
missed, fixed before landing: COM:723 cited BCFForestCombiner's doc
comment rather than the declaration (true line 731), and the seven
loglik M-cell refusal anchors plus both ranges were one line short
(generics.R:120 -> 121, 120-125 -> 121-126; the file's own
convention cites the stop() refusal site, not the enclosing else,
while the S-cells of the same dispatch chain cite the branch-guard
lines). Process note for the runbook: the implementer spawned a
subagent that exceeded its remit and raced it on the branch;
convergence of two entangled passes is NOT independent verification
- the orchestrator sampling pass stays load-bearing.

### P12 consumer bundle + revdep-smoke repoint (142a50b6, 2026-08-18)

All four downstream packages (Fork 3 puts bairrtt in the bundle)
built and checked against tip 8e1e674c, each agent in its own
private lib from a pristine git-archive stage:

- stan4bart bartcore@54e157b: --preclean install with zero dbarts.h
  diagnostics, no C-API version mismatch on load, tinytest 531/531
  (it ships tinytest, not testthat), R CMD check OK - zero E/W/N.
- treatSens dbarts-1.0@1d7c697 (its dbarts-1.0 worktree): compile
  clean, testthat 186/186, a runtime smoke over the eight flat-C
  entries finite and correctly shaped, --as-cran matching the
  2026-07-22 baseline (pre-existing compilation-flags WARNING and
  CRAN-history NOTE only; an ERROR manifests only when NOT_CRAN=true
  is forced past a skip_on_cran guarding CRAN's 2-core cap).
- bairrtt main@6167423: all nine dbarts touchpoints inventoried and
  signature-stable at the tip, including the joint per-observation
  mutation path that is its distinguishing dependency; tinytest
  206/206; needs NOTHING - no compat branch, no fix-list.
- bartCause dbarts-1.0: the sweep's only real findings, neither a
  dbarts regression. (1) Three test-06-regression snapshots stale
  from 31e52644's intentional prior-init draw shift - regenerated
  in-suite against the tip and pushed to the compat branch
  (bartCause 1b37f91; suite FAIL 0, with an unchanged constant
  reproduced exactly as the fidelity control). (2) The
  stan4bart-dependent test's garbage-column ABI error traced to the
  user lib's July-built stan4bart 0.0-14, not to the tip: with the
  tip-built stan4bart on the lib path, test-11 runs 24/24 green.
  Standing flag: any installed stan4bart predating the dbarts.h
  reshapes needs a --preclean rebuild.

revdep-smoke.yaml now shallow-clones each consumer's compat branch
(stan4bart bartcore, bartCause dbarts-1.0, and treatSens dbarts-1.0
- ADDED, because the header's "stan4bart is the only
compiled-boundary consumer" rationale was stale; treatSens crosses
the boundary with eight entries) instead of fetching CRAN sources
that cannot pair with dev dbarts; the matrix branch refs flip to the
repos' mains at the lockstep merge. Validated by YAML parse and live
clone sanity of all three branch/URL pairs. One registration fact
learned: workflow_dispatch requires the workflow file on the DEFAULT
branch, so the smoke is fully dark - not even manually dispatchable
- until the merge registers it (the header's "workflow_dispatch
only" note was optimistic). CI on 142a50b6 itself: all six green.


### P15 + P14 - freshness gate and vignette leg; WAVE 4's content closes (60c92fa8, 2026-08-18)

P15 (local half): tools/check-doc-freshness.R - INDEX completeness
both directions for docs/design and docs/plans, every
feature-matrix anchor resolved against the live tree through the
matrix's own alias table (file exists, line in range, cited symbol
still present, with a Class::member fallback), and the [f39]
scenario counts recomputed from the harnesses. A dead-reference and
drifted-count detector by design, not a semantic reviewer. At HEAD
it found exactly ONE real drift across the living references:
bart2-argument-consolidation.md had no row in docs/plans/INDEX.md
(fixed in this records commit; the guard now exits 0 - 46 design
docs, 152 plan docs, 498 anchors, 37 symbols, 3 count claims).
Self-tested against an injected phantom row, dead-symbol anchor,
and dead-file anchor - all caught, no false positives. The
inventory feeding the resync pass is at
rc-review-artifacts-2026-08-18/p15-staleness-inventory.md; note its
signal is symbol-adjacent presence, so moved-but-alive line numbers
remain the resync's job. CI half default-branch-gated.

P14 (local half): NO-COMMIT - all three vignettes execute cleanly
everywhere exercised. R CMD build knits all three (47.6s), R CMD
check --as-cran re-builds vignette outputs OK, and direct renders
run 2.5-4.2s each with no warnings. CENSUS CORRECTED: the claim
"vignettes are never built or executed in any CI leg" is stale -
pkgdown.yaml has knitted all three as articles on every push since
the site build landed; the real, deliberate gap is
check-standard.yaml's --no-build-vignettes, whose flip is the CI
half, default-branch-gated. A local gotcha recorded for the
runbook: plain R CMD INSTALL on a source tree silently skips
vignettes (no inst/doc), so tarball-based checks are the only local
executions. Run log at
rc-review-artifacts-2026-08-18/p14-vignette-run.log.

With these, wave 4's content items are done; what remains of the
wave is the FULL-NAMESPACE anchor resync, deliberately last.

### Wave-4 batch - harness extension, statistical mode, P13 (4a42620a, 2026-08-18)

Three individually-gated slices stacked and pushed under the batch
clause (disjoint files; one merged-tree battery re-ran everything on
the stack).

HARNESS EXTENSION (the re-scoped rider, discharged): equivalence.R
gains four scenarios - bart2gauss (formula bart2 with non-default
weights + subset, unreachable from the x/y family scenarios),
bart2probit (offset/offset.test + two pooled chains), bart2twoforest
(the canonical zf:forest(x1 + x2) term sugar, K = 2
amplitude-coupled, recording forestFits/glue/forest-widened
varcount), and mixedmatrix (a dbartsMixedMatrix dense+sparse
container through the sampler API, train/test/varcount). Seeds and
sizes literal outside the guarded settings (settingsList()
unchanged); the implementer also caught and fixed a presence-only
yhat.train guard that would have leaked a new channel into every
bart()-routed existing scenario - the 38/38 neutrality compare is
the proof it did not. Recording at the landed sha becomes
equivalence-4a42620a.rds (42 scenarios), self-compared 42/42 strict
and cross-checked bitwise against the implementer's own-build
recording; four places bumped in this commit. Budget honestly over
(175 dense vs ~90 - four scenarios plus three fitters; the estimate
was low, the same class as P6's data-row underweighting).

STATISTICAL MODE (TODO ledger item, local half): bcf-equivalence.R
and multinomial-equivalence.R now fall through from a bitwise
mismatch to equivalence.R's shape - Welch z per summary over each
scenario's draws-axis statChannels (BCF sigma/train/varcount;
multinomial train/test/runVarcount) against baseline-recorded
summaries, |z| >= 4 fails; the record path writes the summaries
forward. Snapshot channels (live-sampler mu/tau/glue/varcount
reads, accept/install verdicts) have no distribution to test and
keep gating bitwise. LOUD degradation when a baseline predates
summaries - never a silent pass. Verified: this host's verdicts
byte-identical (12/12, 10/10, identical stdout); mutation-proven on
both harnesses (BCF glue-precision mutant max |z| 42.6, multinomial
omega mutant 12.9-61.7, both exiting through the statistical
verdict against summaries-carrying scratch baselines). The TODO
entry records the two-part residual honestly: the pinned baselines
predate summaries (re-record only valid from the recording host),
AND a re-record alone cannot make a cross-host run pass - the
snapshot channels mismatch across hosts by construction, so
cross-host needs a design decision (exempt-under-flag vs
draws-axis conversion) not taken here. CI legs stay
default-branch-gated.

P13 (.win drift guard, local half): tools/check-win-drift.R checks
the six checked-in Windows variants - version literals against
DESCRIPTION, and each config header's macro nameset against its .in
template under a RECORDED 38-entry expected-absent table (zero
UNKNOWN reasons; every gap traced to a use site). No real drift at
HEAD; the guard was self-tested against four injected drift classes
(version mismatch, new-check-missing, stale-now-defined,
stale-dropped-from-.in), all caught. CI half default-branch-gated.

Merged battery on the stack: all eight gates green from a
fresh build (tinytest 5991/0; the win-drift guard exit 0; 38/38
neutrality vs f009eff8; the 42/42 record self-compare AND a 42/42
bitwise cross-check against the implementer's own-build recording;
bcf 12/12 and multinomial 10/10 with the fallback code inert;
air/lintr with zero new lints, the two equivalence.R shifts
confirmed as moved-not-new; tarball builds). CI watched to
six-green at the sha before this records push.

### P5b - mutation-sequence composition arms; WAVE 3 CLOSES (61310e0e, 2026-08-18)

inst/tinytest/test-composition-sequences.R (82 assertions), seeded
from the gate-blindspot audit's uncovered combinations with every
cell's openness RE-DERIVED against the current suite first (the
Step-0 table with per-cell citations is in the implementer report;
covered-since cells were cited, not retested). What landed, each
with a draw-level oracle and a non-vacuity arm: per-observation
setPredictor on probit AND logistic fits (accepted and rolled-back
sessions, continuation checked against an uninterrupted control);
installTrees warm starts on multi-chain, linear-leaf, and
missing-data samplers plus the BCF refusal; setCutPoints composed
with BCF, leaf-covariate columns, and missing data (grid actually
in force, with the still-drawable count as the non-vacuity arm);
multi-chain-reproduces-single-chain draw-level oracles for BCF,
DART, linear and GP (converting shape tests into posterior ones);
BCF between-sweep mutation axes at the R5 surface; and the
round-trip DEPTH arms - fresh-sampler setState continuation on all
five model classes, asserting the move sequence EXACTLY and the
numbers to 1e-9 relative.

Two behaviors pinned that no test stated: (i) continuation past a
fresh-sampler restore is SEMANTIC, not bitwise - train diverges at
1e-15 on the first post-restore sweep from the canonical fit resum
while varcount stays bitwise on all five classes - which is the
SETTLED 2026-07-06 state-continuation decision, so the spec's
"non-bitwise is a finding" clause was superseded by the recorded
design, not weakened; (ii) setPredictor's per-observation path
advances the chain generator even on a complete no-op (it draws
the scan-order permutation before any row decision) where the
whole-column path does not - a Gibbs host reading the mask as an
MH accept mask inherits that advance; pinned, and worth an Rd
sentence when the mutation surface's docs are next touched
(ledgered).

The mutation evidence is the headline: reversing the order
createChainRngs seeds per-chain generators is INVISIBLE to the
entire pre-existing suite - 5909 results, all passing - and the
new multi-chain oracles kill it (8 failures, bcf/dart/linear/gp,
train and sigma pins each). Independently re-planted and confirmed
by the gate-runner: blindness half 5909/5909, kill half 8/82. The
rng-restore mutation is killed by both old and new arms; P6's
m21-m23 documented survivals remain survivals (their guards fire
only on hand-corrupted states no continuation arm constructs - as
P6's note predicted). Budget 231 dense vs ~220. Gates green twice
(implementer + independent gate-runner, fresh builds): tinytest
5991/0, the file 82/82 standalone, trio bitwise 38/38 / 12/12 /
10/10, air/lintr clean, R CMD check --as-cran with the core limit
Status: OK. Tests-only; no NEWS. WAVE 3 IS COMPLETE: P8, P5, P7,
P16, P6, P5b, P11 all landed.

### P11 - reduced-power-gate audit; two quick gates retuned (2b30ba95, 2026-08-18)

The audit ran as a measurement pass against 87bf44af using the P6
driver's poison list: of the 16 poisons, 15 have a quick-capable
killer (tests/cpp has no quick concept). Five kill under quick
(m01, m02, m04, m14, m15). Eight are N-A BY DESIGN: the equivalence
family's compare mode refuses a quick run against a full-recorded
baseline loudly - the settings-identity guard stops before touching
any scenario (verified live: exit 1, "baseline was recorded with
different settings") - so quick equivalence cannot silently
underpower, which was itself worth confirming. Exactly TWO were
MIS-TUNED, each confirmed killed at full power on the same mutated
build (no new gate gap): change-balance quick at nKept 100000 read
|z| = 3.6 on the poison-3 forward-correction drop against its own
|z| < 4 gate, and bcf-exact-weak quick's toleranceProb 0.05 sat
above the poison-13 a-glue-precision gap of 0.0314. Full table at
rc-review-artifacts-2026-08-18/p11-reduced-power.md.

The follow-through landed as 2b30ba95: nKept 100000 -> 125000
(z scales ~sqrt(n)) and quick toleranceProb 0.05 -> 0.03, comments
stating the kill constraint with the measured margins. Gated
measured, not assumed: clean-tip pass at the new settings with
comfortable margins (largest clean |z| 1.1 against the 4 gate; clean
gap 0.0141 against 0.03), and both mutated builds now KILLED under
quick (m03 |z| = 4.5 exit 1; m13 gap 0.0314 > 0.03 exit 1), builds
verified fresh through the build-freshness guard. air/lintr clean,
diff confined to the two gate scripts. P11's weekly full-mode CI
schedule half stays deferred on the default-branch gate, as
specified.

### P6 - the repeatable mutation harness lands (b46c0704, 2026-08-18)

The durable asset the July poison sweep never committed:
benchmarks/R/mutation-battery.R, a self-contained driver (modes
list / verify-anchors / run all|ids, --keep-going) holding 23
entries - the 16 poisons re-derived BY SYMBOL against the current
tree (16/16 still live; the BCF glue pair followed the relocation
into combiner.hpp), four R/-level mutations with tinytest killers,
and three SURVIVE_DOCUMENTED entries from P7's untested-path feed
(bartcore_bridge::setState currentSampleNum and dart-skip
negative-checks, readWarmStartState's k length-check) that document
the measured state-restore gap executably. Every mutation lands on
a fresh git-archive copy - the working tree is never touched -
builds --preclean into a per-run lib, passes the wave-0
build-freshness guard, then runs its DESIGNATED killers
(equivalence kills via EQUIVALENCE_SCENARIOS subsets); kill
verdicts come from completed killer exits with harness errors
distinguished and retried, and each entry's copy is discarded
before the next. Demonstrated run: 23 entries, 20 KILL_EXPECTED
killed, 3 SURVIVE_DOCUMENTED survived, 0 unexpected, ~23 min
(p6-run.log preserved in the artifacts dir).

Honest deviations and finds. Budget: 423.5 dense vs the spec's
~250 - the estimate underweighted the mutation list's data rows;
the implementer trimmed genuinely, stopped at the cap, and
over-reported rather than gamed. Two first-attempt mutations
survived their designated killers and were re-targeted with the
reason recorded in-file: poison 3's forward-correction drop must
target the single-sided ordinal/categorical branch change-balance's
MIXED arm actually exercises AND needs change-balance FULL - quick
mode's Monte Carlo error is too loose to kill it, a live
reduced-power datum carried into P11; two R/ candidates were
C++-backstopped or logically undetectable and were replaced by
purely R-level checks. The implementer also caught a real driver
bug mid-build (R's system2(env=) wants "NAME=value" strings, not
named vectors) that had equivalence kills passing vacuously via
subprocess crashes - fixed, and the verifier's audit confirms kills
now read the harness's own completed exit. Verification: my diff
review, then an independent cold-start pass - inventory 23,
verify-anchors clean, subset rerun (m10 equivalence-kill, m17
tinytest-kill, m21 documented-survival) all as designed, worktree
untouched before and after, air/lintr clean. Landed at b46c0704
after a clean rebase (patch byte-identical) with verify-anchors
re-run at the rebased HEAD - the anchors survive compose's
combiner.hpp additions. Benchmarks-only: no NEWS, no baseline
motion; CI legs that fire are exact-gates and lint.

### Compose slice - per-forest veto weights reach the prior initializer (1e020abb, 2026-08-18)

The residue VD routed at the fix landing, widened by its blind
critique, closed in one slice. Chain::sampleTreesFromPrior now
conditions each forest's rejection draw on the vector that forest's
OWN moves veto against: a new weights-only combiner virtual,
formForestVetoWeights (BCF: w * m^2 under the near-zero snap, so
creation glue b0 = 0 zeroes every control row; multinomial:
omega * activeRows against its null global vector), composed with
setForestWeights' per-forest channel inside the loop - NOT
formForestResponse, whose documented invariant owes an immediately
preceding drawForestGlue no initializer performs. The empty
conditioning event is settled UNCONDITIONALLY before the tree loop
by one O(n) scan: a forest whose composed vector carries no positive
entry takes bare roots with zero parameters (the unique structure no
later weight restore can strand a member-empty leaf in), which
removes the live fault the critique measured at the fix's tip
(setActiveRows all-zero + sampleTreesFromPrior exhausted the attempt
cap on a legal state); the cap stays as the backstop, its comment
rewritten - exhaustion now means the scan and the predicate
disagree. growTreeFromRoot's occupancy law moved from member count
to positive weight in the scans it reads (scanOrdinalCuts,
scanCategoricalPartitions), so $growFromRoot can no longer land the
illegal states the prior-init fix closed; the scans are reachable
ONLY from grow-from-root, so no MH stream moves. sbc.R's BCF theta0
generator now draws AND INSTALLS the glue (state round-trip) before
drawing prior trees - under the composed law the tree prior is
glue-dependent, and the old order drew theta0 from two inconsistent
priors.

NOT re-record class, confirmed twice: no recorded baseline scenario
calls sampleTreesFromPrior or growFromRoot, and the trio is bitwise
at the landed tree - 38/38 (xbart included) / 12/12 / 10/10.
Mutation evidence, one deviation from the spec's plan: the mandated
intersection-law mutation is PROVABLY INERT rather than undetected -
every combiner forms its per-forest weights by multiplying the
chain's, so composed support is a SUBSET of global support in every
shipped configuration and the intersection IS the composed set
(draw digests byte-identical over four fixtures); the implementer
substituted two mutations that discriminate for real regressions:
conditioning on the global vector alone (the pre-slice law) is
killed by the composed-law arm's tauUnreached = 0 against 670
forbidden leaves, and reverting the scan occupancy to member count
is killed by grownUnreached = 0 against 1839/1900. The >= scan
mutation and the forestWeights-only scan mutation are killed by the
attempt-cap fault arms. New tests: test-prior-init-composed-law.R
(default-glue BCF fixture - creation amplitudes (1, 0, 1) make the
defect the DEFAULT construction - per-forest oracle through the
internal reader, both forests' tree priors stated so the prior
cannot stand in for the law) plus tests/cpp extensions to
test_grow/test_scan.

Gates: implementer battery green, then ALL 12 re-run independently
on the rebase (1e020abb) by a separate gate-runner from its own
fresh build: tests/cpp 246/246, ASAN zero reports, tinytest 5909/0,
trio bitwise as above, geweke-mc FULL exit 0 (largest |z| 1.90),
backfit-exact PASS, bcf-exact quick OK, sbc BCF smoke end-to-end,
law arms 19/19 and 3/3 with their non-vacuity contrasts, air/lintr
clean, NEWS parses at 273, R CMD check --as-cran Status: OK. CI at the
sha: five green (sanitizers included - the engine-landing watch);
pkgdown stalled at 68 minutes in the known mirror-stall class,
cancelled and rerun live (attempt 2), concluding under the standing
timeout guards. The setForestWeights residue and the critique's four
findings are all DISCHARGED; the empty-leaf-veto design doc carries
the composed-law and bare-root-guard paragraphs.

### SBC-deepening - nbinom third ladder point flags; adjudicated MIXING, mechanism demonstrated (2026-08-18)

The owed family-tiers item, unfrozen by the prior-init fix and run
at the corrected theta0 law: sbc.R nbinom R = 200, L = 150, thin =
300 (10x spacing; 9009s wall). avg.mu PASSES cleanly (ecdf-diff
0.058 in the 0.128 band); r and agg.psi STILL FLAG (ecdf-diff
0.1519 / 0.1419, chisq p 0.000, agg.psi with a 3.6x bottom-rank
spike). The earlier "crosses into the band at 5x spacing" reading
rested on reduced-R ladder points with wider bands - at full R the
failure is spacing-invariant, which strained the recorded H-MIX
adjudication enough to route it, per the multinomial f_ik
precedent, to an exact-conditional derivation rather than more SBC.

The adjudication (independent Opus pass, everything tool-verified)
came back MIXING - the H-MIX reading CONFIRMED and upgraded from
hypothesis to demonstrated mechanism, with NO defect anywhere:
(1) the r conditional is a 13-point capped integer grid (NBResponse
/ NBDispersionPrior, model.hpp), and the code matches the
derivation term by term - kernel lgamma sums, collapsed
-log1pexp statistic, grid prior, current-not-stale conditioning
(refreshLatents draws r, then omega at the NEW r, then trees);
(2) the SBC generator is exact - sbc.R's grid and weights are
bit-for-bit the engine's, the likelihood transcription matches to
5e-14, and burn is 24000 SWEEPS at both thins, so the ladder was a
fair test; (3) negbin-exact.R PASSES on the tip (max grid gap
0.0100 vs tol 0.025) and it DOES gate the full grid r posterior -
at single-tree two-cell scale; (4) the mechanism: at n = 150 the
r-vs-psi ridge is effectively NON-ERGODIC - integrated
autocorrelation 2199-6359 sweeps for r (ESS ~7-21 in the thin=300
window, ~1-2 at thin=30), and two chains on one dataset held
DISJOINT basins (r in {10,12,15} vs r pinned at 50) for 100000
sweeps with identical avg.mu - zero crossings, so thinning cannot
manufacture transitions that never occur; the cold start at the
grid median makes the rank error one-sided exactly as observed
(large r0 -> rank at L, small r0 -> ties -> uniform); (5) the
CONTROL: the identical harness at n = 20, where tau(r) is 35-157,
passes EVERYTHING at R = 200 (r ecdf-diff 0.0360, all in the
0.0924 band, no end spikes). avg.mu = r e^psi is the identified
coordinate and mixes at lag ~1, which is why it is clean at every
spacing.

RECORDED CONSEQUENCES: the nbinom r/agg.psi flag is a KNOWN MIXING
LIMITATION at ensemble n, not a correctness bug; the honest SBC
gate for nbinom calibration is the small-n configuration. DOOR
(VD-timed, evidence in hand, not designed unprompted): the joint
ridge move the sweep lacks - propose r' from the grid with a
compensating psi shift -log(r'/r) on the ensemble, MH-accepted -
the same remedy family the tiers plan already named for this ridge
and ordinal's; re-record class if built. Freeze protocol NOT
invoked (pre-existing recorded flag, no baseline neutrality
evidence rests on these functionals, and the adjudication found no
defect) - orchestrator judgment, flagged for VD. Logs:
rc-review-artifacts-2026-08-18/sbc-nbinom-thin300.log; the
adjudicator's rank/trace artifacts under the session scratchpad
(nbderiv-*).

### P16 - xbart hardening; the loss plumbing gets its oracle (f009eff8, 2026-08-18)

The least-verified surface hardened. inst/tinytest/test-xbart-oracle.R
(16 assertions): a capturing loss records the exact (y.test,
testSamples, weights) triples xbart hands its loss, and rmse /
weighted rmse / log / mcr are transcribed BY HAND from those triples
and fold-averaged outside xbart, plus two pure-arithmetic fold-size
pins (12 rows -> 4,4,4 avg 4; 3,3,2,2,2 avg 2.4). Discrimination
proven twice: the implementer's planted mutations (fold-count
off-by-one kills 6/16; swapped log branches + reversed rmse rows
kills exactly the right 3/16 with mcr unaffected) and the
independent gate-runner's own re-plant (6/16). Weights pins: non-unit
weights move the loss (max abs draw gap 1.19); unit-vs-absent is
equal only to re-association (measured 4.4e-16 draws / 6.8e-16
loss), pinned at 1e-10 with an exact fold-partition identity.
n.samples = 0 now refused by name in R/xbart.R (census 5b closed;
the stale array-dimension pattern in test-xbart-error.R moved with
it; NEWS item, 272 entries parse). The xbart equivalence scenario
lands per the re-scope: fitViaXbart drives the k-fold loop over a
(n.trees x k) grid, loss array as the recorded channel, sizes
literal, settingsList() untouched.

SPEC DEVIATION, measured and stopped per instruction: the mandated
1-thread == n-thread pin CANNOT exist. xbart seeds per chunk
(chunkSeeds = sample.int over numChunks = min(n.threads, n.reps)),
so at n.reps > 1 results are thread-count-dependent BY DESIGN
(measured: rep 1 reproduces across 1/2/4 threads, rep 2 diverges);
the design is documented in R/xbart.R and already gated as an
inequality in test-xbart-reproducibility.R. The reproducibility
contract is (seed, n.threads) pairs. DOOR, not designed here:
per-replication seeding would make xbart thread-invariant at an
RNG-stream cost - a surface-behavior call for VD if wanted before
1.0-0. Also discharged by drift: the spec's "patterns for the 50
bare errors" was already landed by P8's repo-wide pinning.

Gates: implementer battery green (tinytest 5890/0; trio vs
31e52644/6e3b9fb8/1027be5 = 37 identical + xbart uncovered / 12/12
/ 10/10; scenario record+self-compare 38/38 strict; air/lintr clean
with only equivalence.R's two pre-existing warnings; R CMD check
--as-cran with _R_CHECK_LIMIT_CORES_=TRUE Status: OK) and re-run
INDEPENDENTLY by a separate gate-runner from its own fresh build,
ALL GREEN. Landed at f009eff8 (R-only; engine binary is 31e52644's).
This records commit installs equivalence-f009eff8.rds (the
gate-runner's fresh-build recording, self-compared 38/38) as
canonical with the four-place obligation; the xbart scenario's P17
oracle is the hand-computed loss oracle landed in the same commit.

### Prior-initializer fix LANDS (31e52644, 2026-08-18) + main-baseline re-record

VD's decisions recorded: LAND the prepared fix; the setForestWeights
residue COMPOSES with an empty-set guard (the all-zero-forest-weight
case adjudicated below, slice in flight); the P9 freeze scope as
applied during the fire is RATIFIED as recorded - no landed slice
re-gates. ccf4b687 rebased onto c98e5c14 with zero conflicts as
31e52644; the rebased patch is hunk-for-hunk byte-identical (only
inst/NEWS.Rd's base blob moved - intervening entries sit outside the
fix's hunks). The full battery re-ran INDEPENDENTLY on the rebase:
install, tests/cpp 246/246, ASAN clean, tinytest 5874/0 (the new
pin: 0/13225 forbidden leaves, non-vacuity arm 3062), trio partition
exactly as prepared (33/37 bitwise; movers grouped 2.26 /
grouped_aft 1.24 / hazard 2.30 / hurdle 3.01, all under 4; BCF
12/12, multinomial 10/10 bitwise), geweke-mc FULL exit 0 (largest
tree-shaped |z| 1.90 vs the pre-rebase run's 1.2 - the trio's
bitwise match shows the build's draws are identical, so the spread
is the oracle harness's own run-to-run variance, inside its bands),
backfit-exact PASS, bcf-exact quick OK (gaps <= 0.0009), air/lintr
clean, NEWS parses at 271, R CMD check --as-cran Status: OK. CI
six-green at the sha (pkgdown, lint, cpp-tests, exact-gates,
R-CMD-check, sanitizers).

Re-record at the final sha under P17's oracle-naming rule:
equivalence-31e52644.rds recorded from a fresh --preclean build,
self-compare 37/37 under --strict-coverage, bitwise-identical to the
pre-rebase candidate 37/37, partition vs 8b047f8b as above.
bcf-equivalence-6e3b9fb8 and multinomial-equivalence-1027be5 are NOT
re-recorded: both bitwise at this tip (12/12, 10/10), and the
MANIFEST convention (the 21fc29c precedent) keeps unmoved baselines
current with neutrality noted - the three-baseline obligation is
discharged as fresh record+compare evidence for all three with only
the moved harness renamed. Four places bumped in this commit:
equivalence.yaml main line, MANIFEST (new current row, 8b047f8b
superseded), TODO, feature-matrix f39. The riding scenario additions
do NOT ride this re-record: bart2 gaussian/probit/two-forest plus
the B2 mixed-matrix scenario become a dedicated harness-extension
slice (harness extensions self-record at their own sha, the 33f6fdc
S0 precedent; queued behind the compose slice), the xbart scenario
stays with P16 - TODO updated so nothing drops.

The compose adjudication ran ahead of its slice: an independent
blind critique (own code reading, own probes) CONFIRMED the residue
and WIDENED it. (1) The empty conditioning event is reachable at
THIS tip with no combiner at all - setActiveRows all-zero plus an
explicit sampleTreesFromPrior exhausts priorTreeDrawMaxAttempts and
throws, on a state test-active-rows-pins.R documents as
accepted-and-runs - so the guard must be UNCONDITIONAL, and the
attempt-cap comment's "unreachable" claim is wrong until the compose
slice lands. (2) BCF's creation-default amplitudes are {1, 0, 1}, so
the DEFAULT dbarts(forests=) build zeroes control rows in the
treatment forest's veto vector - measured 8/1263 forbidden leaves
over 50 prior draws. (3) growTreeFromRoot still enforces
non-emptiness by MEMBER COUNT while growForestFromRoot hands it the
composed vector, so $growFromRoot can land the illegal states this
fix closed. (4) sbc.R's BCF theta0 generator must draw AND install
glue before the prior trees once the tree prior is glue-dependent.
None of it blocks this landing (no test or consumer reaches the
fault; the primary law is strictly better) and none of it is
re-record class (no recorded baseline scenario calls
sampleTreesFromPrior or growFromRoot - they drive bartcoreRun from a
bare-root start). The all-zero adjudication: BARE ROOT with zero
parameters - the conditional law does not exist over an empty event,
the kernel freezes structure under an all-zero vector and recovers
on restore (both measured), and the bare root is the unique
structure legal under ANY later weight restore. Slice spec critiqued
and in flight. UNFROZEN by this landing: SBC-deepening, P16, P6.

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
xbart pin still has no oracle, P16's job).

### P7 - branch-reach one-shot; measurement only (2026-08-18)

Nothing committed to the package by design. llvm-cov over an
instrumented --preclean install plus instrumented tests/cpp, both
suites run under instrumentation (tests/cpp 246/246, tinytest
5871/0 - coverage did not perturb behavior), scoped to the
project's own sources. Headline: 11,902 regions 91.10% hit, 7,735
branches 83.36%; engine 93.34%/89.47%; the bridge is the laggard
at 87.63%/74.59% with R_interface_bartcore.cpp worst. All 411
never-hit function-level clusters classified against the P8
census: 372 UNTESTED-PATH (the P6 mutation-target feed - the two
largest clusters in the codebase are bartcore_bridge::setState and
readWarmStartState, a systemic state-restore gap beyond the
census's spot checks), 16 DEAD-CANDIDATE (incl. one non-census
find: bartcore::Tree::rightChildOf has zero callers in-tree - a
Family 2 removal candidate), 23 GUARD (expected-unexecuted refusal
paths). Full classified report at
rc-review-artifacts-2026-08-18/p7-branch-reach.md. Explicitly not
a coverage gate.

### P5 - executable composition matrix (7997e50c, 2026-08-18)

benchmarks/R/composition-matrix.R (740 lines): parses
feature-matrix.md's tables 1-4 programmatically (header/row-order
drift fails loudly), executes all 180 S and ? cells across the 13
families with minimal seeded fixtures (n.threads = 1L throughout),
classifies CONSTRUCTS / REFUSES / ERRORS, and exits nonzero on any
disagreement with the matrix's claim. First run: 175 confirmations,
ZERO disagreements, five ? resolutions - four match their footnotes'
expectations (the two silently-wrong pointwiseLoglik cells per
f19/f28, both grouped/hetero constructions per f30) and ONE
corrected stale doc prose: student + variance REFUSES at spec.R:423
(the wave-0 FX2 refusal) where footnote f30 still claimed
"constructs today with no refusal site". This records commit
updates the cell to R spec.R:423, rewrites f30 to the measured
split verdict, and fixes the student-section prose - the runner's
Rd-vs-behavior coherence role exercised on its first run. Gates:
runner end-to-end clean from a fresh Rscript; air/lintr clean; diff
was the one new file; tinytest unchanged 5871/0 after the run (no
fixture state leaked). Run log preserved at
rc-review-artifacts-2026-08-18/p5-run.log. Prior CI: P8 six-green.

### P8 - refusal census + anti-vacuity pass (5d06f641, 2026-08-18)

Wave 3 opens. Census over 912 sites (887 refusals, 25 warnings;
L's 941-row enumeration reconciled by content): 749 refusals
PINNED by a message-matching assertion, 138 unreached - 32 DEAD
(foreclosed by earlier checks or self-documented unreachable), 30
HARD (interrupt/alloc/raw-entry/rethrow), 79 GENUINE GAPS, the
real finding, concentrated in the bridge's entry validation and
R/generics.R's per-generic argument checks; unreached verdicts
verified by two independent read-only passes with 2 hand
corrections, including a 5b wrong-message finding (xbart's
n.samples = 0 path errors about 'x' dimensions). The 73 live bare
expect_error assertions (77 had
drifted) all gain literal patterns from their sites' captured
messages, plus 6 bare expect_warning the new gate flagged and 10
opportunistic pins at trivially-reachable gap sites. .lintr gains
bare_expect_pattern_linter (XPath, both named and positional
pattern absence) - repo-wide lint clean with zero nolint
exceptions. Strengthening proven twice on disjoint samples (5 + 4
broken messages -> patterned assertions fail, bare ones had
passed). Gate-runner's --as-cran check caught one
environment-fragile pin - xbart's auto n.threads trips CRAN's
core limit before the target validation - fixed with an explicit
n.threads = 1L and the as-cran-for-thread-spawning-tests lesson
recorded in the runbook. Battery: tinytest 5871/0 (with and
without _R_CHECK_LIMIT_CORES_); diff touches only inst/tinytest/
and .lintr; air clean; NEWS untouched parses 270; R CMD check
--as-cran Status OK. The 79 genuine gaps are the follow-up ledger
for the remaining wave-3 passes. Prior CI: K6 six-green (its
exact-gates hit the NEW 25-minute timeout on a mirror stall at
setup-r-dependencies - the guard's first live catch - and went
green on rerun).

### K6 - boundary guard reconciliation; WAVE 2 CLOSES (ee249420, 2026-08-18)

The predicate-parity policy lands as error-style.md's R8 addendum
(R8 had already set message ownership - R canonical, C
independently-worded backstop - and deferred predicates here): the
non-canonical side must accept and reject the identical input set;
the C backstop stays because the low-level bartcore* .Call entries
skip R validation entirely; an R predicate over a possibly-NA/NaN
value refuses it deterministically in ONE guard, never a separate
anyNA pre-check relying on evaluation order (NA_real_ >= 0.0 is NA
and if (NA) raises the wrong error). Census correction by content:
the case-weight NaN "gap" was stale - a separate NA guard already
refused it, now merged into one bridge-matching predicate
(any(is.na(w) | w < 0.0)); the REAL gap was setForestWeights
admitting Inf past the R guard to fail three frames deep in a
tryCatch rethrow - now !all(is.finite(w)) || any(w < 0.0) refuses
at the surface with the finiteness named. Audit of the other
boundary pairs (setActiveRows, multinomial counts/labels/
category-offset, survival-time): R-only or already
parity-conforming. Case-weight Inf is ACCEPTED on both routes -
parity holds; its semantic sense rides the wave-3 refusal census.
Eight pins; one NEWS item (270). Battery twice (implementer, then
independent gate-runner with its own route x value parity matrix
on both libs): tests/cpp clean; tinytest 5860/0; trio bitwise
37/37 12/12 10/10; air/lintr clean; R CMD check Status OK;
discrimination - base's Inf forest weight errors in doTryCatch
three frames deep, slice refuses at the surface frame. Prior CI:
L six-green.

WAVE 2 COMPLETE: L (941 messages, 111 rewritten under the signed
rule) + K6. Wave 3 opens with P8, pinning against settled text as
the ordering rationale prescribed.

### L - the error-message sweep under the signed rule (efdfa74a, 2026-08-18)

Wave 2 opens. Enumeration: 941 reachable messages (570 R stop/
warning, 344 C Rf_error/Rf_warning, 27 indirect C literal
definitions traced to their defining sites); 111 REWRITTEN (46 R,
65 C), 830 already conform - the census's ~120 estimate matched
the rewrite bucket, not the corpus. Rewrite classes: quoting,
lowercase-initial, no terminal period, at most one explanatory or
remedy clause (3-clause messages trimmed with the remedy kept),
consistent C-side caller-colon (41 hardcoded-literal sites - the
multinomial mutation channel predated the rule - plus 9 dynamic
prefixes), settled templates for the six recurring shapes, and
already-in-hand value enrichments where the rule's templates name
the expected quantity.
14 message-pinned assertions updated across 10 test files; one
NEWS item (269 parse). Battery twice (implementer, then
independent gate-runner whose LOAD-BEARING gate was text-only
verification: R side by parse-tree comparison with stop/warning
string literals normalized to placeholders - every divergence
confined inside those calls' argument lists, zero control-flow or
condition changes; C side by reading all 53 hunks - format-string
text and wrapping only): install --preclean; tests/cpp clean;
tinytest 5852/0; trio bitwise 37/37 12/12 10/10; air clean; lintr
clean; R CMD check --as-cran Status OK; 20-message style
spot-audit all conforming. Prior CI: FX1-channel and both oracle
pushes fully green (P1b's cancelled lint covered by FX1's
containment).

### Prior-initializer fix PREPARED, NOT LANDED (ccf4b687 on wt/prior-init-fix, 2026-08-18)

The adjudicated fix implemented and fully gated in a preserved
worktree (session-local, not retained); LANDING IS VD-HELD
because every prior-initializing fit's stream moves (re-record
class). Engine delta ~21 dense lines: Chain::sampleTreesFromPrior
regrows each tree whole until Tree::bottomNodesHaveWeight (new,
the veto's own leafHasNoWeight law, distinguished in its comment
from the restore path's membership law) accepts it - rejection,
not projection, whole-tree so the conditioning tilts the parent's
rule draw; priorTreeDrawMaxAttempts = 10000 faults loudly and is
unreachable unless the conditioning event is empty (acceptance
bounded below by P(bare root) = 1 - base). The mid-chain collapse
sites are repair, not prior draws, and keep the member count -
left alone. DECISIVE GATES: geweke-mc full exits 0 with every
tree-shaped z inside |1.2| (was |3.7-6.6|); backfit-exact PASS
byte-identical to its landing note; bcf-exact quick at the
documented gaps. Trio vs canonical: 33/37 bitwise, the four
movers are exactly the prior-initializing fits (grouped,
grouped_aft, hazard, hurdle; max |z| 3.01 on 1 of 85 summaries -
posteriors statistically unmoved); BCF 12/12 and multinomial
10/10 bitwise. Candidate baselines recorded and self-compared
bitwise 37/37 12/12 10/10, preserved with the run logs and sweep
classifications.
Snapshots regenerated (test-reproducibility-rbart.R via the tool;
test-rbart-loop-callback.R's 15 literals by hand - MANIFEST wants
this sha at landing); the calibration non-vacuity bar sat ON its
statistic at 3.0 and moves to 2.5 with the measurement recorded.
New pin test-prior-init-empty-leaves.R: 0/13,225 forbidden leaves
under a zero-weight half-space, non-vacuity arm 3,062. Full
battery otherwise green (tests/cpp 246, ASAN zero, tinytest
5855/0, air/lintr clean, NEWS 269 parses, R CMD check OK).
RESIDUE FOR VD, deliberately outside the adjudicated shape: the
predicate reads workingWeights, so per-forest setForestWeights
zeroes can still leave a from-prior leaf that forest's veto
forbids; composing it is one line but an all-zero forest weight
is a legal tested state whose conditioning event is EMPTY, so it
needs an empty-set-guard adjudication first. At landing the fix
rebases over the message sweep (NEWS append point conflicts,
trivial) and the re-record carries the four-place obligation.

### P9 - Geweke marginal-conditional oracle lands FIRING; diagnosis CONFIRMED (2278c929, 2026-08-18)

benchmarks/R/geweke-mc.R (644 lines, 2:57 full / 16s quick) lands
recording a REAL defect, so its full run exits 1 against this tip
by design. Design: iid prior draws with simulated y versus a
stationary chain alternating one sweep with a between-sweep
setResponse refresh - stationary from step 0, no burn-in
assumption; prior source is the engine's own machinery
(sampleTreesFromPrior + sampleNodeParametersFromPrior; sigma
transcribed via sbc.R's sbcSigmaDraw); 11 functionals with two
exact pivots; sbc.R band machinery; 2000 chains x 200 steps x 2
arms. Teeth measured: three poisons fire at |z| 8-23, and the
skip-refresh negative control correctly does not. Response side
CLEAN everywhere; the tree-shaped functionals FIRE.

DIAGNOSIS (independent adversarial adjudication, own code reading
+ own probe, CONFIRMED): sampleTreesFromPrior draws the
unrestricted CGM prior then PROJECTS via collapseEmptyNodesBelow,
while the moves price the CGM prior RESTRICTED and renormalized to
the empty-leaf-free set (logLikelihoodForBranch's leafHasNoWeight
veto, docs/design/empty-leaf-veto.md) - projection != conditioning.
Sharpest signature: P(single-leaf tree) = 0.04993 from the
initializer vs 0.05492 from the sweep against theory 0.0500 /
0.0552 (z = 17.61). The sweep realizes the DOCUMENTED deliberate
target; the initializer realizes a third, undocumented law - the
initializer is the deviant. KERNEL PURITY: Chain::run never calls
the initializer path; the mid-chain collapse sites are forced
data-mutation repair outside the kernel; converged posteriors are
unaffected. BLAST RADIUS: fit initialization + rbart_vi warmup +
n.grow.sweeps (init-only, benign); samplePriorPredictive and the
$sampleTreesFromPrior surface R AND C API (directly wrong as
prior-draw products); SBC theta0 (~1% assumption violation, below
SBC's resolution). Baselines and landed neutrality evidence NOT
impeached (build-vs-build determinism; law-level correctness never
enters). SECONDARY DEFECT, same path: collapse triggers on
numObservations() == 0 but the veto on positive WEIGHT, so under a
zero-weight half-space 1954/2000 from-prior forests hold a leaf
the moves forbid (0/2000) and a birth from it compares -HUGE_VAL
to -HUGE_VAL - a NaN ratio that rejects and never faults. FIX
SHAPE (adjudicated, NOT implemented here): per-tree rejection loop
on Tree::leafHasNoWeight - exactly p(T | no empty leaf) because
the event factorizes over the recursion - fixing both findings;
reject rate ~9% at n=30, ~1.3% at n=400; moves the initial
forest's RNG stream everywhere (snapshot + baseline re-record
class, posterior summaries unchanged). The fix slice is being
prepared GATED IN A WORKTREE with candidate baselines; LANDING IS
VD-HELD. Freeze scope stands as recorded in the fire note;
SBC-deepening remains frozen until the fix lands.

### P9 ORACLE FIRE - initializer/prior inconsistency (2026-08-18, ADJUDICATION IN FLIGHT)

The Geweke marginal-conditional oracle (benchmarks/R/geweke-mc.R,
built and run on wt/geweke-mc at da018e3f, NOT yet landed) FIRED on
the tree-shaped functionals while the response side and both exact
pivots stayed clean in both arms: n.leaves z = -4.41/-6.62,
leaf.depth -3.66/-5.79 against the band. The implementer's
adjudication: sampleTreesFromPrior does not realize the
forest-structure law the Metropolis moves price. An in-script
probe at sigma pinned to 1e4 (structure-move likelihood ratios ~1)
reads the gap directly - initializer expected leaves 19.650 vs
sweep-stationary 19.504 at m = 10 - and a scaling probe shows the
relative gap tracking the incidence of prior draws carrying a
leaf of <= 2 observations (-1.11% at n = 30 shrinking to z =
-0.12 at n = 1600). Mechanism HYPOTHESIS, unproven:
growSubtreeFromPrior grows into empty children by rule-based
availability and collapseEmptyNodes projects them out - a
projection of the unrestricted law, where the birth/death moves
price a restricted-and-renormalized one. An INDEPENDENT
adversarial adjudicator (own code reading + own standalone probe)
is running; P9's landing waits for its verdict.

FREEZE SCOPE APPLIED (orchestrator judgment under VD's AFK
discretion grant, FLAGGED for VD review): the fire impeaches the
initializer's distributional claim, not the sweep (response side +
exact pivots clean) and not any recorded baseline - every landed
neutrality evidence item is bitwise determinism against recorded
streams, and the exact-posterior gates run past burn-in where
initialization bias decays; therefore NO landed slice re-gates and
the message lane (L, K6 - string edits with no engine interaction)
CONTINUES. SBC-deepening is FROZEN: its theta0 draws come from the
impeached initializer, so tree-shaped SBC evidence generated now
would be against a known-slightly-wrong prior sampler. The FIX is
an engine change that moves draws (baseline re-record class) and
is VD-HELD: direction (initializer -> moves' law vs both -> the
documented prior) rides the adjudicator's report to VD.

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
classification preserved untracked (session scratch);
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
move from docs/plans/archive/bcf-b-ridge.md sections 2.2/2.3/5a into
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
and without test, ordinal thresholds, nbinom dispersion, BCF
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
merge-time shim design. Full report untracked.

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
