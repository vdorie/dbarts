# runsbcbcf-repair

status: OPEN (diagnosed 2026-08-05: engine-side, decision-gated; the
        TODO's "small, separate" reading is rescoped by the diagnosis)
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
