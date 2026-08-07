# runsbcbcf-repair

status: LANDED 62caed0 2026-08-05, acceptance PASS (FIX-B per the
        Design + critique with all three blocking amendments; setModel
        guard precondition 6744aca; the R=200 acceptance adjudicated
        across thin=30/90/120 - see the Acceptance addendum)
agent: sonnet diagnosed; opus if an engine capability ships
rng: neutral - a creation-time calibration override defaults to the
     current derivation, so existing draws are bitwise-unchanged
budget: decision first; if fix, ~30-60 engine lines + harness adaptation

## Goal

`runSbcBCF` (benchmarks/R/sbc.R ~1101) runs again, or the reason it
cannot is recorded and the arm is retired from the harness surface.

## Diagnosis (2026-08-05)

- Repro: `Rscript benchmarks/R/sbc.R bcf 3 10 1` fails in
  .bcfSetResponse: "bartcore_setResponse: a multi-forest sampler fixes
  its data at creation; make a new sampler instead".
- Cause: 9030d93 added refuseMultiForestMutation
  (src/R_interface_bartcore.cpp ~1756), refusing
  setResponse/setData/setWeights at numForests >= 2 - correctly, since
  whole-data mutation rebuilds only forest 0. runSbcBCF structurally
  depends on the refused pattern: one sampler built from the fixed
  anchor yBuild, reused across reps via setResponse(updateScale=FALSE).
  The FALSE is load-bearing (bcf-sigma-residual.md): build-scale
  carryover keeps a, aVariance, b0, b1 and sigma on one prior across
  reps.
- Why rebuild-per-rep (the ordinal/nbinom/multinomial adaptation) is
  NOT valid here: BCF's per-forest leaf-scale calibration
  (scaledResponseSd, src/bartcore/chain.hpp ~645) derives from
  y - offset's range at construction (GaussianResponse::rescale,
  src/bartcore/model.hpp ~2800) and no API sets it independently of
  the fit target. A per-rep rebuild calibrates each rep's posterior to
  sd(y0) while theta0 came from the fixed anchor scale - prior and
  posterior no longer share one prior, the SBC precondition.
- Scope: BCF is the only arm where response-derived calibration +
  per-rep response swap + multi-forest collide. The sbc.yaml CI matrix
  (discrete-selfcheck, gaussian, ordinal, nbinom, t, multinom) is
  unaffected; runSbcBCF is manual-only.
- ABI: bcfParams is .Call-bridge surface (an 8-double vector in
  R_interface_bartcore.cpp), NOT in inst/include/dbarts/dbarts.h, so a
  creation-time override has no ABI checklist / lockstep cost.

## Decision (VD)

- FIX-A, calibration-scale override at creation: bcfParams grows one
  slot (or an explicit arg) pinning the leaf-scale calibration;
  default = derive from y as today (bitwise-neutral). runSbcBCF then
  rebuilds per rep with the anchor scale pinned. Restores a
  re-runnable whole-posterior gate for BCF. ~30-60 engine lines +
  harness adaptation; the modified .Call entry re-earns its rchk note.
- FIX-B, scale-pinned multi-forest setResponse: teach the mutation
  path to rebuild ALL forests' fit state when updateScale = FALSE.
  Surgery on the mutation transaction the 9030d93 guard exists to
  protect; strictly larger and riskier than FIX-A for the same
  harness-facing capability.
- DOCUMENT: retire runSbcBCF with a pointer here; the recorded BCF SBC
  acceptance runs (bcf-sigma-residual.md, bcf-ridge-interweaving.md)
  predate the guard and stand. BCF loses re-runnable SBC until an
  engine capability ships.

## Constraints

- The 9030d93 guard stays for the general case; only an explicit
  opt-in may bypass what it protects against.
- Existing families and default BCF draws bitwise-unchanged
  (equivalence trio + bcf-equivalence identical under FIX-A default).
- No shipped-header (dbarts.h) change.

## Steps

1. VD picks FIX-A / FIX-B / DOCUMENT.
2. If fix: implement, adapt runSbcBCF, tiny-config smoke, then a real
   `bcf` arm run at the recorded config; equivalence trio +
   bcf-equivalence bitwise; rchk note for the touched .Call entry.

## Verification

- Diagnosis verified by reproduction at HEAD and code reading of the
  guard, the calibration path, and the recorded updateScale=FALSE
  rationale; worktree left clean (no code change shipped with this
  record).

## Survey addendum (2026-08-05, multi-forest setResponse as a feature)

An independent survey (opus, read-only, working probes) of whether
FIX-B has value beyond the harness. Findings that AMEND the Decision
section above:

- Guard surface, verified: refuseMultiForestMutation has exactly
  three sites - setResponse (~3060), setData (~3087), setWeights
  (~3276) in R_interface_bartcore.cpp. setOffset (~3045) is UNGUARDED
  and, on BCF, CORRECT: Chain::setOffset touches no forest and
  BCFForestCombiner::formForestResponse re-derives every per-forest
  residual from y each sweep, caching nothing. Proven by a working
  two-outcome SUR-BCF probe (setOffset + setSigma + run(0,1) +
  inverse-Wishart), which answers the open question at
  docs/design/correlated-outcomes.md ~132-134: residual-conditioning
  embeddings (multilevel/SUR) ALREADY work through offset; only
  latent-RESPONSE schemes need the guarded path.
- FIX-A as stated above is INSUFFICIENT: sbcMakeBCF recovers
  fitScale/fitShift once, so pinning the leaf scale s alone leaves
  the response range map unpinned; a correct override must pin the
  whole calibration (scale AND shift).
- FIX-B cost, revised downward with a caveat: the "rebuild ALL
  forests' fit state" premise is vacuous for BCF (nothing per-forest
  is cached across sweeps), estimated ~35 engine lines / 60-120
  total, no ABI change, and it fixes runSbcBCF verbatim. This
  contradicts the diagnosis's "strictly larger and riskier" reading;
  if FIX-B is picked the arc still runs design + blind critique,
  which must resolve what Chain::setResponse itself must touch per
  chain (residual/totalFits bookkeeping, scale handling at
  updateScale=FALSE) before trusting the low estimate.
- Constituency, verified on-package: TobitBART (Imports: dbarts)
  hand-wires per-sweep setResponse/setSigma/setWeights across two
  dbarts samplers; Chen et al.'s AoAS replication builds five dbarts
  samplers 6000 times over for want of a response swap. Binary BCF
  demand is real in the literature but served on CRAN since 2026
  (stochtree outcome_model, flexBART). Censored, robust-t,
  missing-data, measurement-error and multinomial demand are
  UNASSESSED (the survey retracted a fabricated citation batch rather
  than launder it; only tool-verified claims kept).
- New doors, outside this item: setModel is UNGUARDED at
  numForests >= 2 and rewrites forest 0 only, silently discarding
  BCF's calibrated mu leaf scale - a correctness gap on the shipped
  surface wanting a guard (or a real multi-forest path); setOffset
  and setTreatment are silent no-ops on multinomial and should refuse
  or work.

## Design + critique (2026-08-05)

Process: design (opus, working probes built out-of-tree against
src/*.a) then an independent blind critique (opus) instructed to
refute, which rebuilt the probe from scratch, re-derived the state
walk, and ran the counterfactual (+43/-3 patch on a scratch clone).
Full documents lived in the session scratchpad (disposable); this
section is the durable record.

VERDICT on the cost discrepancy: both prior readings wrong, in
opposite directions. There is NOTHING to rebuild: at
updateScale=FALSE, GaussianResponse::setResponse (model.hpp ~2670)
re-maps yRescaled_ through the pinned (min_, range_) and touches no
forest, and BCFForestCombiner re-derives every per-forest residual
from y each sweep, caching nothing per-forest. Chain::setResponse is
ALREADY correct for BCF; FIX-B is ~40 lines of opt-in predicate
plumbing (combiner virtual supportsResponseMutation(), default FALSE
on the base so future multi-forest models stay refused until audited,
BCF override true; the bridge refuses at numForests >= 2 &&
(updateScale || !supports)). Decisive probe, independently reproduced
by the critique: setResponse(yNew, false) is BITWISE identical
(sigma, both forests' totalFits, glue a/b0/b1, 270 sweeps) to the
already-unguarded setOffset(yBuild - yNew, false) - the same re-map
with a different pointer, so the guard forbids on one entry what it
permits on another. fitScale/fitShift are pinned by construction on
the FALSE branch (rescale() is never called), the whole-calibration
pin FIX-A lacked; sigma, sigmaSqPrior, the per-forest leaf
calibration and the glue all carry over (the one-prior-across-reps
SBC precondition). runSbcBCF runs VERBATIM; zero harness edits; no
man/*.Rd change (the documented R5 surface cannot reach a multi-forest
sampler); rchk risk nil.

Critique BLOCKING amendments (all bind the implementation):
1. The flat C API hole is TWO holes, not one (dbarts_sampler_
   setResponse at C_interface.cpp ~194 AND dbarts_sampler_setWeights
   at ~204; both unreachable today - ResponseFamily has no
   multi-forest member), and refuseMultiForestMutation has internal
   linkage (anonymous namespace), so closing them needs the helper
   promoted into bartcore_bridge + declared in
   R_interface_bartcore_common.hpp, not a one-liner.
2. ORDER INVERSION: the setModel guard is a PRECONDITION of FIX-B.
   "No mutation path writes forest.leaf.scale" is false at HEAD -
   Chain::setModel (chain.hpp ~1122) writes it and is reachable from
   the unguarded bartcore_setModel. The guard ships FIRST; FIX-B's
   pin argument holds only once it exists.
3. The gaussian-response conjunct must be structural, not a comment:
   the chain-level permission is combiner_->supportsResponseMutation()
   && family == gaussian (a latent family would read
   forests_[0].totalFits as the combined location).
Non-blocking, for the implementer: the reused refusal message says
"fixes its data" - reword for the prior/opt-in case; stale comment at
model.hpp ~3383-3384; NA updateScale must refuse (asLogical NA would
take the permitted branch); multinomial setOffset is the same silent
no-op the design refuses for setResponse (door, not taken); setState/
installForests restore fitMin/fitMax with leaf.scale pinned (sibling
door, not taken); spec.sdModerate is discarded at construction
(chain.hpp ~636 its only read) but recoverable from the live scales
(s = mu.leaf.scale * sqrt(mu.numTrees), sdModerate =
tau.leaf.scale * sqrt(tau.numTrees) * 0.674 / s), so updateScale=TRUE
is implementable in principle; the governing refusal is d489986's -
recomputing s mid-run is a data-dependent prior refresh
(2026-08-06 correction); the tau-forest leaf scale
is untestable from tests/cpp without a new accessor (Chain::leaf() is
forest-0-only).

COUNTERFACTUAL (critique, patched scratch clone): the smoke
`Rscript benchmarks/R/sbc.R bcf 3 10 1` RUNS, 15/15 PASS; equivalence
27/27, bcf-equivalence 5/5, multinomial-equivalence 3/3 all identical
(default draws bitwise-neutral); tinytest 1 fail of 3474, exactly the
predicted refusal assertion at test-multi-forest-seam.R:75-78 (flips
to expect_silent in the implementation); tests/cpp pass; a short arm
`bcf 25 100 10` gave 13/15 with both flags explained (the harness's
documented sign-symmetric `a` functional, abs.a PASS; prog3 at ~1
bin-sd on R=25), so the R=200 recorded-config arm is the real
acceptance gate.

Confirmed refusals, each with evidence: updateScale=TRUE decalibrates
(fitScale 6.42 -> 64.16 with leaf.scale unmoved) and stays refused;
multinomial setResponse is an empty override (silent no-op) and stays
refused at both updateScale values; setData stays refused
(applyNewData hard-codes forest 0). setWeights is vacuously guarded
for BCF by the same argument as setResponse - a door recorded, not
taken.

## Guard landing (2026-08-05, 6744aca)

The FIX-B precondition shipped: refuseMultiForestMutation at
bartcore_setModel (after the SamplerBase binding, before the class
check, so the refusal fires regardless of the argument's class); the
helper promoted into bartcore_bridge with a declaration in
R_interface_bartcore_common.hpp (critique amendment 1); both flat
C API holes (dbarts_sampler_setResponse / _setWeights, unreachable
today) guarded defensively; shipped header unchanged. Two new seam
tests (BCF setModel refusal, single-forest setModel inert). Gates,
implementer + independent re-run from a fresh --preclean install:
tests/cpp from clean all pass; tinytest 3476/0 (+2); equivalence trio
bitwise identical (27/27, BCF 5 scenarios x 6 channels, multinomial
3 x 5 vs the ec2a3d0 baseline); air format clean; R CMD check
Status: OK with a zero pre/post delta. FIX-B's leaf-scale pin
argument now holds at HEAD.

## FIX-B landing (2026-08-05, 62caed0)

Implemented per the Design + critique with the three blocking
amendments: the flat C API holes and helper promotion had landed with
the guard precondition (6744aca); the chain forwarder carries the
structural gaussian conjunct (combiner opt-in AND family == gaussian);
NA updateScale refuses (the bridge permits only an explicit FALSE).
supportsResponseMutation defaults false on the ForestCombiner base
(future couplings refused until audited), BCF overrides true with its
two conditions stated; the bridge's two refusal messages both carry
"multi-forest" and describe the opt-in era accurately; the stale
multinomial comment at model.hpp ~3383 rewritten. Nine files,
+276/-18 (engine/bridge 84, tests 210).

Gates, implementer + independent re-run from a fresh --preclean
install: tests/cpp from clean 175 checks all pass, including the new
BCF block whose strongest arm reproduces setOffset(yBuild - yNew,
FALSE) bitwise (sigma, both forests' totalFits, glue; exact-delta
construction) and whose totals-vs-tree-sums check documents a
pre-existing swap-independent ~2.2e-8 tau-forest gap (the b_z 1e-9
multiplier floor in the residual cancellation, measured identical
with and without the swap); ASAN/UBSAN tests/cpp leg clean
(implementer); tinytest 3484/0 (+8: seam refusal flipped to success,
TRUE/NA/multinomial refusals added, BCF round-trip retargets with the
reported/internal map identical across the swap); equivalence trio
bitwise identical (27/27, BCF 5x6, multinomial 3x5 vs ec2a3d0); smoke
`sbc.R bcf 3 10 1` RUNS 15/15 PASS (5.6 s/rep); air clean; R CMD
check Status: OK, zero delta (implementer). runSbcBCF needed zero
edits. rchk: NOT dispatchable pre-merge (rchk.yaml is unregistered
until it lands on the default branch, the known limitation in its
header); the design assessed the entry's rchk risk as nil (no SEXP
allocation or PROTECT change, retain unmoved) and the release
procedure's as-cran rchk run covers it. The R=200 recorded-config
acceptance arm runs post-landing, verdict recorded below as an
addendum.

## Acceptance addendum (2026-08-05, R=200 arms at 62caed0)

Three R=200 L=200 runs on the landed build (5.9-7.0 s/rep; prior
moment checks PASS identically in all three: sigma coverage 0.8998 vs
0.90, glue Cauchy IQR 3.9899 vs 4.0, b sd 0.7065 vs 0.7071), spanning
thinning regimes thin=30/90/120 with burn pinned at 72000 absolute
sweeps; per-functional 5% band 0.0917 (the CLI band - the sbc.yaml
matrix's Bonferroni discipline would widen it to ~0.12).

- The two sign-ill-posed functionals a and b1.minus.b0 FLAG at every
  regime with their absolute counterparts PASS (abs.a 0.046-0.081,
  abs.diff 0.051-0.082) - exactly the recorded pre-guard baseline
  (bcf-sigma-residual.md sec 3, "FLAG by design").
- Every other cell: thin=30 put eff1/prog3/prog4 at 1.00-1.16x the
  band; the 3x ladder (thin=90) released prog4 but retained
  eff1/prog3 (the two runs share pinned seeds, so persistence there
  is partially shared randomness, not independent evidence); the
  recorded acceptance regime thin=120 RESOLVES eff1 (0.0592) and
  prog3 (0.0431) while prog2 (0.1119) and eff3 (0.1028) sit marginal
  instead. No well-posed cell is elevated across regimes - which 2-3
  cells sit at 1.0-1.2x the per-functional band reshuffles with the
  regime, the multiplicity signature of 13 correlated per-row cells
  at alpha = 0.05 - and every well-posed cell at every regime clears
  the Bonferroni band. The recorded healthy interim R=200 showed the
  same picture (one marginal cell, adjudicated PASS at R=400).

VERDICT: PASS. FIX-B restores a healthy re-runnable BCF SBC gate:
ranks uniform at the recorded regime, sign-symmetric flags matching
the recorded baseline. Scope note: a defect here could not have
implicated FIX-B in any case - the permitted swap is bitwise
setOffset(yBuild - yNew, FALSE), and BCF's draw stream is pinned
bitwise by bcf-equivalence-99205ee across every landing this run.
Door, not taken: the runSbcBCF CLI band is per-functional 5%;
adopting the matrix's Bonferroni band would stop wandering marginal
cells from reading as flags.

## setData door survey (2026-08-06, at f592572)

The setData refusal above was surveyed (design memo plus an
independent refuting critique, both against live code) after VD
reopened the door on enabling-models grounds: mutation on multi-forest
samplers can enable model classes ahead of any packaged consumer.
Verdict: the door stays OPEN, held at the joint x/y/z semantics with
n free; any taken form MUST carry z in the same call. Findings:

- Three holes, not one. G1 structural (the recorded one):
  recoverTreeParameters and applyNewData address forests_[0] only
  (chain.hpp ~1456, ~1472). G2 scale: GaussianResponse::setData
  rescales unconditionally (model.hpp ~2695) - the decalibration
  setResponse(TRUE) is refused for; a taken door needs the pinned
  shape, and restoreScale re-anchors the prior but not sigma
  (~2769/2785 vs ~2686-2698), so a naive getScale/setData/restoreScale
  sandwich double-anchors. G3 memory safety, previously unrecorded:
  the combiner indexes borrowed z (holder.ownedTreatment, sized at
  creation or setTreatment) over live numObservations
  (combiner.hpp ~482/501/525/662) and setData carries no z channel -
  an n-growing per-forest setData over-reads the heap, and a fixed-n
  one silently re-pairs the old z with new rows.
- The split remap is variable-preserving: mapOldCutPointsOntoNew
  reads rule.variableIndex and never writes it (tree.hpp ~1201,
  ~1244), so a tree's used-column set can only shrink and tau's
  columnMask containment is invariant under applyNewData. The state
  path's containment checks (chain.hpp ~2377) guard unconstrained
  donors, not the remap; a per-forest applyNewData needs neither.
- Already served on BCF at fixed n: forced setPredictor (the
  transactional refusal gates only forceUpdate != TRUE), setTreatment,
  setResponse/setOffset scale-pinned under the opt-in, setWeights,
  and setCutPoints - which has and needs no multi-forest guard and
  installs an arbitrary grid including a quantile shrink. The door's
  unique remainder: a row-count change, the nearest-value split
  remap, and joint-swap atomicity.
- Calibration fork to resolve at design time: pinned (construction's
  calibration persists across the swap - the accepted posture of
  every pinned mutation) vs recalibrating (recompute s: rejected by
  d489986 as a data-dependent prior refresh).
- Demand today is zero (no consumer calls setData; none creates a
  BCF sampler; stan4bart/bartCause/treatSens grepped at their compat
  heads). Recorded as fact, not as the gate - the enabling rationale
  governs.
- The motivating class - Bayesian sensitivity analysis with latent
  treatment, the treated set resampled per iteration - runs TODAY at
  fixed n: the glue's default (a = 1, b0 = 0, b1 = 1;
  combiner.hpp ~293) held fixed (updateA/updateB = false) is exactly
  y = mu + z tau; a b0 = 0 control row's multiplier is floored to
  1e-9 so its tau-suffstat weight is ~1e-18, effectively excluding
  it; and setTreatment swaps z mid-run. What that leaves on the
  table - exact zero-weight exclusion and tau paying O(n) per sweep
  regardless of treated count - is a narrower, efficiency-shaped
  door (per-forest row subsetting), distinct from whole-data setData.
