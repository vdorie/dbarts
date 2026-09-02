# multinomial-counts

agent: opus (a LIKELIHOOD arc: the count-matrix generalizes the softmax data
  model from single-trial labels to n_i-trial rows via PG(n_i, .) = sum of n_i
  PG(1, .) draws. The care is that the n_i = 1 reduction stays BITWISE on every
  existing multinomial path - the interleaved PG draw and the working-response
  construction both live on the per-forest hot hook, where a reordered draw is a
  silent posterior change - and that the new count path is correct against its
  own exact-posterior arm).
rng: NEW CAPABILITY, existing paths NEUTRAL. The n_i > 1 count path is a new
  likelihood; every EXISTING path stays bitwise: single-forest equivalence 22/22
  IDENTICAL vs equivalence-ac6ec2c.rds, BCF 5x6 IDENTICAL vs
  bcf-equivalence-99205ee.rds, and the multinomial fixture's TWO single-trial
  scenarios (all 5 channels) IDENTICAL vs multinomial-equivalence-5afb09a.rds.
  The new count path is gated by a NEW exact-posterior arm and a NEW fixture
  scenario; the 5afb09a baseline is NOT re-recorded for the single-trial
  scenarios (they reproduce inside the compare), only re-recorded to ADD the
  count scenario. A moved draw on any anchor means the n_i = 1 reduction is not
  byte-identical - stop and rethink the loop shape.
window: NONE. Internal engine + bartcore .Call + the bart2 matrix R surface;
  dbarts.h FROZEN (multinomial is internal, stan4bart never runs it). family =
  "multinomial" already ships (bb29d00); the count matrix is a new response FORM
  under it, no new family slot.
budget: ~520-620 lines. C1 ~380 (engine ~120: count-native combiner - the
  PG(n_i) summing hook, the (y_if - n_i/2) working response, MultinomialSpec
  counts/trials fields, the ctor; bridge ~90: a count entry point parsing the
  n x K matrix + deriving/validating trials + ownedCounts/ownedTrials; internal
  R ~30; exact-gate count arm ~60; fixture count scenario + re-record ~40;
  tests/cpp ~40). C2 ~140 (public matrix interface + generics/print + tinytest).
  C3 ~40 docs. Chiefly src/bartcore/combiner.hpp, chain.hpp,
  R_interface_bartcore.cpp, R/bartcore.R, R/bart.R, R/generics.R,
  benchmarks/R/multinomial-exact.R + multinomial-equivalence.R + a re-recorded
  baseline, benchmarks/baselines/MANIFEST, tests/cpp/test_sampler.cpp,
  inst/tinytest/test-multinomial-surface.R. Header edits -> --preclean; delete
  stale tests/cpp binaries.

## Goal

The shipped multinomial softmax model fits GROUPED counts: an observation is a
row of an n x K nonnegative integer count matrix with row sum n_i (trials), and
category k's one-vs-rest conditional is binomial(n_i, sigmoid(eta_ik)) with
success count y_ik, augmented by omega_ik ~ PG(n_i, eta_ik). The single-trial
label path (n_i = 1, one-hot rows) is the exact reduction and stays bitwise. The
capability is reachable via the matrix interface bart2(x.train, Y, family =
"multinomial") with Y an n x K count matrix. Closes the count-matrix follow-up
the C7 landing filed (docs/plans/archive/multinomial.md C7 scope-narrowings;
docs/plans/archive/multi-forest-models.md; TODO).

## Ordering (this arc lands BEFORE multinomial-formula; Q1)

The count matrix completes the multinomial DATA CONTRACT (factor label OR n x K
counts). Landing it first lets the formula arc design one ingestion pass against
both response forms (a factor LHS and a cbind(...) matrix LHS) instead of
shipping a factor-only formula surface and re-touching it when counts arrive.
This is also the riskier (engine, Opus) half; the formula arc is then a single
RNG-neutral surface commit. See multinomial-formula.md, which is gated on this
arc for its count branch. The counter-case (formula is the more visible
usability gap) is real but loose: formula can accept counts additively later, so
the churn saved by counts-first is a second ingestion touch, not a rewrite. Q1
records the fork; VD picks the routing order, neither plan is internally blocked.

## The augmentation reuses the shipped sampler, no new numerical code

- PG(n_i, psi) for INTEGER n_i is the sum of n_i iid PG(1, psi) draws; the only
  shipped sampler is ext_rng_simulatePolyaGamma(rng, psi) = PG(1, psi)
  ([[src/include/external/random.h:126@6a48351b]]; [[src/external/random.c:581@6a48351b]]). There is NO
  direct PG(n, .) sampler. LogisticResponse already sums a weight's worth of
  PG(1, psi) draws ([[model.hpp:2236-2242@6a48351b]], `reps = lround(weights_[i])`); the count
  combiner does the same per (i, f) with reps = n_i. Cost is O(n_i) draws per
  (i, f) per sweep - linear in TOTAL counts (the accepted weighted-logistic cost,
  weighted-binary.md), never in the per-observation kernels.
- NO r-update analog constrains this path. Negative-binomial's caveat (integer
  vs real PG shape; negative-binomial.md, TODO) is that its shape y_i + r
  contains a SAMPLED dispersion r that goes non-integer. Here the PG shape is
  n_i, the OBSERVED trial count - fixed integer DATA, never sampled - so the
  summing sampler is always exact and no real-shape gap opens. The only
  real-shape concern is fractional/weighted counts, which stay OUT OF SCOPE
  exactly as case weights were refused for single-trial (the same real-shape PG
  gap weighted-binary/negative-binomial carry).

## The data contract

- Y is an n x K matrix of nonnegative integers; trials n_i = sum_k y_ik.
- K = 2 is a binomial(n_i, p) per row; the softmax reduces to a two-forest
  logistic on counts (distributionally the weighted-logistic path, not
  bitwise - K forests vs one, different rng consumption).
- Single-trial reduction: a one-hot Y with every n_i = 1 IS the current label
  path. This is the neutrality anchor (below).
- n_i = 0 (empty row): PG(0, .) is a point mass at 0, so omega_ik = 0 breaks the
  (y_if - 0)/0 working response. REFUSE at ingestion (require n_i >= 1); Q2.

## Context (seams, read in code)

- The combiner reads the response ONLY in two functions; combinedFits (651-665),
  combinedTestFits (676-689), afterCombine (710-743), and the margin (748-757)
  are count-INDEPENDENT (they read forest fits, so the softmax output, test
  blend, and level-centering are unchanged):
  - drawForestGlue(f) ([[combiner.hpp:620-630@6a48351b]]): today one PG(1, .) draw per
    observation, `omega[i] = ext_rng_simulatePolyaGamma(rng, fFits[i] - margin)`.
    Becomes the LogisticResponse summing loop with reps = trials_[i].
  - formForestResponse(f) ([[combiner.hpp:635-647@6a48351b]]): today
    `yif = labels_[i] == f ? 1.0 : 0.0; forestResponse_[i] = (yif - 0.5)/omega
    + margins_[i]`. Becomes `yif = counts_[f*n + i]; (yif - trials_[i]*0.5)/omega
    + margins_[i]`.
- MultinomialSpec ([[combiner.hpp:181-191@6a48351b]]) carries `const int* labels`; the ctor
  (591-602) stores labels_, cold-starts omega_ at 1/4, sizes n x K scratch.
  Combiner members at 760-766.
- MultinomialResponse ([[model.hpp:2339-2382@6a48351b]]) is unchanged - it owns no labels and
  no counts (the combiner does), and its seams are already vestigial no-ops.
- Chain multinomial ctor ([[chain.hpp:413-434@6a48351b]]) builds K forests + the combiner
  from the spec; createMultinomialSampler ([[facade.hpp:545-551@6a48351b]]) is the single
  ConstantGaussianLeaf instantiation.
- createMultinomialHolder ([[R_interface_bartcore.cpp:1670-1748@6a48351b]]): parses labelsExpr
  (integer code per observation, validated 0..K-1 at 1699-1708), builds
  MultinomialSpec, moves the buffer into holder->ownedLabels (1744). The .Call is
  bartcore_createMultinomial (2012-2016). bartcoreMultinomialSampler
  ([[R/bartcore.R:559-583@6a48351b]]) is the internal wrapper.
- The fixture (benchmarks/R/multinomial-equivalence.R) drives the internal
  surface at two SINGLE-TRIAL scenarios (k3 covariate 104-132, k2 logistic
  135-151), recording 5 channels; the exact gate (multinomial-exact.R) has three
  train arms + one test arm over single-trial labels.

## Design

The combiner becomes COUNT-NATIVE: it holds counts_ (n x K int, category-major
to match omega_) and trials_ (n). The single-trial label entry re-expresses
labels as a one-hot counts matrix with trials_ == 1 at the bridge, so the two
draw-reading functions reduce EXACTLY at n_i = 1:

- drawForestGlue: `double omega = PG(rng, psi); for (long c = 1; c < trials_[i];
  ++c) omega += PG(rng, psi);` At trials_[i] == 1 the loop is empty -> one PG
  draw with the identical psi = fFits[i] - margin. Byte-identical.
- formForestResponse: `(yif - trials_[i]*0.5)/omega + margin`. At trials_[i] == 1,
  yif in {0,1}, trials*0.5 == 0.5 -> identical doubles.

This is the neutrality guarantee: the single-trial fixture scenarios reproduce
5afb09a bitwise because the reduced code is the same draw stream. UNIFY (one
count-native path) over a DUAL representation (keep labels_, branch on which is
set): unify is one code path, provably bitwise by the reduction; dual is
marginally safer-by-construction for the gate but a permanent branch on the hot
hook. The fixture is the proof either way - the established multinomial-arc
discipline (the seam is exercised end to end only through surface + fixture).

Bridge: keep bartcore_createMultinomial (labelsExpr) as the single-trial entry -
it converts labels to a one-hot counts buffer + trials == 1 and fills the spec's
counts/trials pointers (so its DRAWS are the reduction, proven by the fixture);
add bartcore_createMultinomialCounts(countsExpr, numCategoriesExpr) parsing an
n x K integer matrix, deriving trials as row sums, validating nonneg and n_i >= 1
and column count == K, storing holder->ownedCounts (n*K) and ownedTrials (n).
MultinomialSpec swaps `const int* labels` for `const int* counts` + `const int*
trials` (the label entry builds the one-hot). One new .Call entry -> one rchk
note.

## Commits

C1. The count-matrix likelihood, one gated commit (the arc's C4 discipline: the
   engine seam is unreachable from R without the bridge + internal surface, so it
   lands with the gate that exercises it). Sub-parts:

   (a) ENGINE (combiner.hpp): MultinomialSpec counts/trials fields; count-native
   ctor (build one-hot from labels in the bridge, so the ctor takes counts +
   trials); drawForestGlue PG(n_i) summing loop; formForestResponse
   (y_if - n_i/2) working response; members counts_/trials_ replacing labels_.
   combinedFits/combinedTestFits/afterCombine/margin UNCHANGED.

   (b) BRIDGE (R_interface_bartcore.cpp): the label entry converts to one-hot +
   trials == 1; a new bartcore_createMultinomialCounts parsing/validating the
   n x K matrix (nonneg integers, n_i >= 1, K columns), ownedCounts/ownedTrials.

   (c) INTERNAL R (R/bartcore.R): bartcoreMultinomialSampler unchanged (label
   path); add bartcoreMultinomialCountSampler(sampler, counts, K) routing to the
   new .Call, validating the matrix R-side (safe over fast in R).

   (d) EXACT-GATE ARM (benchmarks/R/multinomial-exact.R): an intercept-only
   K = 3 COUNT arm (n_i > 1, fixed counts), matching the posterior mean of the
   identified softmax probabilities to the same correlated-difference-Gaussian
   quadrature the label intercept arm uses (the likelihood is now multinomial(n_i,
   p) per row; the quadrature integrates that). Note the K = 2 count arm reduces
   to binomial(n_i, p) - state it, gate the K = 3 count arm. This is the only
   exact gate on the PG(n_i) summing and the (y - n_i/2) working response.

   (e) FIXTURE (benchmarks/R/multinomial-equivalence.R): add a k3-COUNTS scenario
   (n_i > 1) recording the 5 channels; re-record as
   multinomial-equivalence-<C1-hash>.rds, demote 5afb09a to historical, MANIFEST
   NEUTRALITY note (the two single-trial scenarios reproduce 5afb09a bitwise
   inside the compare - the n_i = 1 reduction is byte-identical; only the count
   scenario is new).

   (f) COMPONENT TESTS (tests/cpp/test_sampler.cpp, beside testMultinomial): a
   count fit's drawForestGlue omega positivity + the PG(n) MEAN moment (sum of n
   PG(1, psi) has mean n * PG(1, psi)-mean, weighted-logistic's moment style);
   the SINGLE-TRIAL REDUCTION (a one-hot counts fit reproduces the label fit's
   working-response/weights construction bit-for-bit at a fixed seed - the
   in-process neutrality proof); the K = 2 count == binomial reduction; a
   growForestFromRoot case for the count branch (the G4 lesson - the combiner
   branch is unreachable from R).

   Files above. Gate: all THREE existing anchors IDENTICAL (single-forest 22/22,
   BCF 5x6, the multinomial fixture's two single-trial scenarios all 5 channels)
   + the new count scenario identical to itself + the exact gate's count arm to
   MC error (existing arms unchanged) + tests/cpp from make clean + full tinytest
   (no regen - existing paths neutral) + air + rchk note (bartcore_create-
   MultinomialCounts is a new .Call entry). Size: XL. --preclean; delete
   tests/cpp binaries. Abort: any existing anchor or single-trial channel moves =
   the n_i = 1 reduction is not byte-identical.

C2. Public matrix interface (RNG-neutral surface). bart2(family = "multinomial")
   accepts a count matrix as the response (the y.train positional): detect Y with
   is.matrix(Y) && integer/numeric && ncol(Y) >= 2 in the multinomial branch
   ([[bart.R:472-499@6a48351b]], beside the factor path); levels from colnames(Y) (or
   seq_len(K); Q4); validate nonneg integers + n_i >= 1; route through a
   count analog of bart2Multinomial to bartcoreMultinomialCountSampler.
   packageMultinomialResults threads colnames(Y) as levels on the K margin; the
   fit stores the count matrix as $y. Generics: extract/fitted/predict return
   per-category probabilities exactly as the label fit (the softmax output is
   count-independent); fitted(type = "class") argmax stays defined (a modal-
   category convenience); print reports the count input. tinytest: a count fit
   reproduces the internal bartcoreMultinomialCountSampler channel bit-for-bit at
   a fixed seed (the C7 reproduction pattern); a one-hot count matrix reproduces
   the factor fit's probabilities. Files: R/bart.R, R/generics.R, inst/tinytest/
   test-multinomial-surface.R. Gate: full tinytest (grows) + R CMD check man +
   all three anchors + the multinomial fixture identical. Size: L. Abort: a
   public count fit diverges from the internal path on the same seed.

C3. Docs. docs/design/multinomial.md "The surface": the count-matrix data
   contract (n x K, trials n_i, the PG(n_i) = sum-of-PG(1) additivity, the n_i = 1
   reduction as the neutrality anchor, n_i = 0 refused, K = 2 == binomial, the
   NO-r-update-analog fact). Mark the item landed in
   docs/plans/archive/multi-forest-models.md and the TODO; this plan's landing notes.
   Files: docs/design/*, docs/plans/*, TODO. Gate: full tinytest; R CMD check man
   unaffected (no man/ topic added - the count matrix is a documented response
   form of the existing bart2 topic). Size: S.

## Verification

- Every commit: the three standing anchors bitwise (single-forest, BCF, the
  multinomial fixture's single-trial scenarios). A moved draw = the n_i = 1
  reduction touched a draw path, stop.
- C1 additionally: the exact gate's count arm to MC error, the new fixture count
  scenario identical to itself, and the tests/cpp single-trial reduction case.
- C2 additionally: the public==internal count-fit reproduction, and the one-hot
  count == factor probabilities.
- No bench-sampler: the count path adds work only on the multinomial per-forest
  hook (n_i PG draws vs 1) and only for n_i > 1; single-forest/BCF/single-trial
  multinomial draw the same stream. Note it; skip the quiet-window bench unless a
  reviewer disputes the reduction claim.
- dbarts.h unchanged -> no stan4bart lockstep; C1's new .Call entry earns "rchk
  on next scheduled run".

## Open questions for VD

ADOPTED AS WORKING DEFAULTS (orchestrator, 2026-07-17, VD afk): every
recommendation below - counts-first, refuse n_i = 0, colnames-or-index
naming. Each is branch-local and reversible (refusal can be relaxed
additively; naming moves dimnames only); VD may overrule any of them
on return, before release.

- Q1 (ORDERING: counts-first vs formula-first). RECOMMEND counts-first: the
  formula arc then designs one ingestion pass against the complete data contract
  (factor + count-matrix LHS), and the riskier engine half lands while it is
  fresh; formula collapses to a single RNG-neutral commit. AGAINST: the formula
  interface is the more visible usability gap the C7 landing named. COST of
  formula-first: a second formula-ingestion touch when counts land (additive, not
  a rewrite - the coupling is loose). Gates the orchestrator's routing order;
  neither plan is internally blocked. RECOMMEND counts-first.
- Q2 (n_i = 0 empty rows: refuse vs drop vs carry-zero-weight). RECOMMEND REFUSE
  at ingestion (require n_i >= 1), the simplest contract and "fast over safe in
  C/C++": a zero-trial row carries no information for any channel. AGAINST:
  carry-zero-weight is the statistically coherent reading (an uninformative row),
  but omega = 0 makes the working response NaN and the suffstat kernel would need
  a weight-0 special case on the hot path. Dropping silently is the third option
  but hides a likely user error. Gates C1's data-contract validation. RECOMMEND
  refuse.
- Q4 (count-matrix category naming). RECOMMEND colnames(Y) as the levels when
  present, else as.character(seq_len(K)); matches how a user labels the columns
  and mirrors the factor path's levels(y). AGAINST: none material; state the
  fallback so an unnamed matrix is not silently unlabeled. Gates C2's surface.
  RECOMMEND colnames-or-index.

## Landings

### C1 (2026-07-17, 2bd34db; baseline record 5c44819)

As specced, plus two implementer findings worth keeping. The combiner
is count-native (counts_/trials_, the PG(n_i) summing loop, the
(y - n_i/2) working response); the label bridge entry builds one-hot
counts + unit trials; bartcore_createMultinomialCounts is the new
.Call (rchk on next scheduled run); bartcoreMultinomialCountSampler
is the internal wrapper. Gates, independently re-run at landing:
single-forest 22/22 identical draws, BCF 5x6 identical, multinomial
single-trial scenarios identical vs 5afb09a (the byte-identical
reduction, proven end to end) AND vs the new baseline; exact-gate
count arm max gap 0.0000 (tol 0.008), other arms unchanged;
tests/cpp from clean incl. the bit-for-bit reduction case, the PG(n)
moment, K = 2 == binomial, and a count grow-from-root case; tinytest
2938/0 no regen. Implementer notes: (1) the fixture's count-scenario
seeds are inline literals kept OUT of the guarded seeds vector so
settingsList() stays identical to 5afb09a and the neutrality compare
actually runs; (2) the engine no longer has a label path to diff
against, so the in-process reduction test compares the count path at
unit trials to a hand-rolled replay of the old label arithmetic.
Baseline multinomial-equivalence-2bd34db.rds (3 scenarios x 5
channels) recorded in the follow-up commit per the 5afb09a/d66ece4
precedent; 5afb09a demoted with its neutrality trail.

### C2 (2026-07-17, 38b2482)

bart2(family = "multinomial") accepts the n x K count matrix response
beside the factor path: R-side validation (nonneg whole numbers, no
NA, row sums >= 1, K >= 2, each with its own refusal), levels from
colnames or the index fallback (Q4 as adopted), dispatch through
bart2MultinomialCounts mirroring the factor path's host-sampler
construction exactly. packageMultinomialResults takes levels
explicitly (internal signature only). Generics needed no changes -
none read the response's type. Gates: public==internal and
one-hot==factor reproduction both bit-for-bit in tinytest (2959/0);
fixture + both other anchors identical; R CMD check Rd checks OK; air
clean. In passing: bart.Rd's stale "no varcount" claim corrected (the
per-category channel has shipped since the varcounts arc).

### C3 (2026-07-17)

docs/design/multinomial.md's surface section rewritten around the
count contract (the reduction as the neutrality anchor, the
no-r-update-analog fact, n_i = 0 refused, K = 2 == binomial); TODO's
multi-forest entry updated (formula interface is the one remaining
follow-up). The Q1/Q2/Q4 defaults were adopted during VD's absence as
annotated under Open questions; none was contradicted by
implementation.
