# tau-cauchy-exact-ig

agent: opus
rng: posterior-preserving but DRAW-CHANGING for grouped fits with the
  cauchy tau prior (rbart_vi's default): the grouped and grouped_aft
  equivalence scenarios shift (statistical verdict) and the anchor
  re-records at landing (orchestrator step). Every ungrouped scenario
  and the gamma-tau-prior path must stay BIT-IDENTICAL.
budget: ~60-120 lines engine + component tests + snapshot regen.

## Goal

Replace the Neal step-out/shrinkage slice sampler for the grouped
random-effect scale tau UNDER THE CAUCHY PRIOR ONLY with the
Makalic-Schmidt inverse-gamma scale-mixture two-block Gibbs draw: an
exact conjugate draw with a constant 2 RNG draws per sweep, deleting
the step-out-hang failure class that the 7c1c7c9 cap merely bounds.
Value is robustness, RNG-budget constancy, and code cleanliness -
NOT mixing (docs/plans/tau-slice-review.md measured slice == exact
conjugate in isolation; the real grouped mixing bottleneck is
forest-ranef confounding, filed separately).

## Context

- The review memo docs/plans/tau-slice-review.md section 4a carries
  the full derivation, quadrature-validated under this codebase's
  exact convention (half-Cauchy scale A = 2.5 * rel.scale /
  sigmaScale; tau is the SD of b_j ~ N(0, tau^2), J groups,
  SS = sum b_j^2). Implement from it; do not re-derive from the
  paper without checking the memo's parameterization notes.
- The M&S representation: tau^2 | aux ~ IG with shape (J + 1) / 2
  and the aux-dependent rate; aux | tau^2 ~ IG(1, ...). Both
  conditionals are exact inverse-gamma draws; the auxiliary is drawn
  fresh each sweep from its conditional given the current tau, so NO
  new persistent state field and NO state-format change.
- The GAMMA tau prior KEEPS the slice sampler (gamma-on-SD admits no
  exact draw; the memo's section 4b shows GIG does not apply). The
  slice machinery therefore stays; only the cauchy branch reroutes.
- The custom-prior/callback R path (rbart_vi_run posteriorClosure)
  is untouched by design.
- The step-out cap and tauSliceSteps remain live for the gamma
  branch; note in the code that the cauchy branch no longer consumes
  them.

## Constraints

- Bit-identity everywhere the change cannot reach: all ungrouped
  equivalence scenarios; grouped fits under the GAMMA prior (the
  equivalence "grouped" scenario uses gamma - CHECK which scenarios
  use which prior before predicting the verdict; only cauchy-prior
  grouped paths may shift).
- The grouped-aft reduction gate (rbart_vi aft all-uncensored ==
  gaussian on log T, bitwise) must keep holding - both arms use the
  same new draw, so it survives any correct implementation.
- RNG draw count per sweep documented at the draw site (2 inverse
  gamma draws = 2 ext_rng gamma consumptions; state the exact
  consumption so equivalence reasoning stays predictable).

## Steps

1. Engine: reroute the cauchy branch of the tau draw (chain.hpp) to
   the two-block IG draw per the memo; keep the slice for gamma.
   Component tests (tests/cpp): (a) the sampled tau posterior under
   a fixed b matches 1-D quadrature (promote the memo's validation
   into a checkNear gate, tolerance from MC error at the test's
   draw count); (b) draw-count constancy (a privately seeded rng
   advances by exactly 2 gamma draws per tau update); (c) the gamma
   branch is bit-identical to before (seeded draw comparison).
2. R-level: regenerate any tinytest snapshots that hardcode grouped
   cauchy-prior draws (whole-file replay, never single values);
   confirm the grouped-aft reduction gate still passes bitwise.
3. Docs: rbart.Rd / dbartsPriors.Rd only if they describe the
   sampling mechanism (check; likely no user-visible change).
   Landing note in the plan file.

## Verification

- tests/cpp all pass incl. the new quadrature gate (poison-proof it
  once: perturb the IG rate, watch it fail, restore).
- Full tinytest 0 fails after snapshot regen.
- Equivalence vs equivalence-ac6ec2c.rds: ungrouped scenarios
  IDENTICAL; grouped scenarios' verdict per their prior (gamma
  identical if untouched, cauchy statistical); overall statistical
  pass. Orchestrator re-records the anchor at landing.
- air + lintr on touched R; no dbarts.h change; no state-format
  change.

## Landing

Landed on wt/tau-exact-ig off bartcore 57bbdb5.

- Engine: GroupedResponse::refreshLatents (model.hpp) branches on
  priorKind_. The cauchy branch now calls drawTauCauchyExactIG (a new
  header-only helper); the gamma branch keeps the slice loop verbatim.
  The helper's conditionals as coded, with SS = sum b_j^2, J groups,
  A = priorScale_ = 2.5 * rel.scale / sigmaScale:
    xi    | tau   ~ IG(1,          1/tau^2 + 1/A^2)
    tau^2 | b, xi ~ IG((J + 1)/2,  0.5 SS + 1/xi)
  Each is one reciprocal ext_rng_simulateGamma draw; the primitive
  takes a SCALE, so an IG(shape, rate) is 1 / Gamma(shape, 1/rate). The
  update consumes EXACTLY two ext_rng_simulateGamma draws per sweep (a
  fixed count) and never persists xi - it is redrawn fresh against the
  current tau, so no state-format change. The cauchy branch no longer
  reads sliceSteps_ or the step-out cap; both stay live for gamma.

- IMPORTANT correction to the front matter's RNG prediction: BOTH
  equivalence grouped scenarios (grouped, grouped_aft) were recorded
  with prior = gamma ([[equivalence.R:556@c6601616]], [[equivalence.R:577@c6601616]]; gate-hardening chose gamma
  because the half-Cauchy tail can stall the slice sampler). NO
  equivalence scenario exercises the cauchy grouped path, so the whole
  suite stays BIT-IDENTICAL vs equivalence-ac6ec2c.rds (22/22 identical
  draws, overall OK) - there are no cauchy-grouped z-summaries to
  report and NO anchor re-record is needed. Only the tinytest
  grouped-cauchy tripwire shifted.

- Gates: R CMD INSTALL --preclean clean; tests/cpp all pass incl. the
  three new gates in testGroupedMath (exact-IG quadrature-moment match,
  two-gamma draw-count constancy, gamma-branch bit-identity) - poison
  (0.5 SS -> 1.0 SS in the tau^2 rate) drove the quadrature gate to
  FAIL, then restored; full tinytest 2727/2727 after regenerating
  inst/tinytest/test-reproducibility-rbart.R (its grouped default-prior
  ranef snapshot) by whole-file replay; equivalence 22/22 identical;
  air --check clean (lintr not installed here). rbart.Rd's tau-sampling
  sentence updated to note the exact cauchy draw (checkRd clean;
  pkgdown site is not version-controlled, nothing to regenerate).
