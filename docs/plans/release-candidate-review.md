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
