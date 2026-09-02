# gate-blindspot-audit

agent: one opus poison-sweeper (worktree, never lands) + one sonnet
       matrix reader (read-only); fable adjudicates
rng: neutral (the tree does not change; poisons live and die in a
     throwaway worktree)
budget: poison x gate sensitivity matrix + feature x gate coverage
        matrix, recorded here; gaps feed the SBC review

Review 3 of the retrospective program (retrospective-reviews.md).

## Goal

Measure what the gate battery can actually see. Two defects this
cycle were invisible to every gate until they were found by hand:
the change-move cross-type composition (caught only in meta review;
the MIXED-TYPE arm was added afterward) and the zero-weight sigma
df (self-consistency class: the code disagreed with its own
documented ignore-semantics). A third audit note - zero-weight rows
influence tree structure but not the numeric posterior - is exactly
the kind of thing no current gate would flag. The question is not
"do the gates pass" but "which deliberate breakages do they fail".

## Method

1. POISON SWEEP (opus, own worktree, reverts after every poison,
   NEVER pushes or lands): ~15 deliberate single-site kernel
   breakages spanning the audited surface - birth/death tree-prior
   depth term, node-selection count, change-move forward and
   reverse correction sides separately (the mixed-type arm must
   catch the composition poison), swap validity, sigma df
   (positive-weight count) and SSR weighting, chi-k shape and rate,
   chi hyperprior df, DART Dirichlet update, Polya-Gamma working
   weight, grouped-intercept precision, BCF glue precision and
   per-forest weight transform, calibration constant, linear-leaf
   ridge, GP nugget. For each: apply, rebuild, run the battery
   (tests/cpp; tinytest; equivalence in statistical mode;
   change-balance.R, logistic-reference.R, categorical-exact.R,
   bcf-exact.R as relevant), record WHICH gates fail and which stay
   green, revert, verify clean. A poison that every gate survives
   is a finding.
2. FEATURE x GATE MATRIX (sonnet, read-only): enumerate the feature
   axes (family x weights x offset x DART x grouped x BCF x
   missingness x keepTrees x mutation entry points x multi-chain)
   and map every gate (tests/cpp cases, tinytest files, the 18
   equivalence scenarios, the four exact-posterior gates) onto the
   combinations it exercises. Output the uncovered combinations,
   ranked by how much posterior-relevant code they route through.

Fable adjudicates: gate hardening items for poisons that survived,
SBC targeting from the uncovered combinations.

## Constraints

- One review at a time; these two agents are the review's entire
  fan-out. FOREGROUND only, no sub-agents.
- The poison worktree never lands anything; the sweep ends with
  git status clean and an explicit revert-verification.
- bench-sampler is NOT part of the battery (quiet-machine only).

## Verification

- Both matrices recorded here; surviving poisons become gate items;
  uncovered combos feed sbc-calibration's prioritization.

## Status

- 2026-07-08: started; sweeper and matrix reader dispatched.
- 2026-07-09: FEATURE x GATE MATRIX complete (sonnet reader; the
  poison sweep is still running - this is an interim record of the
  coverage half, to be finalized with the poison matrix). Uncovered
  combinations, ranked by posterior-relevant code routed through
  (each becomes an SBC-review target and/or a gate item):
  1. Linear/GP leaf x missing (NA) leaf-covariate values. The
     imputation "NAs enter at the standardized mean (zero)"
     ([[model.hpp:57@7496d571]], [[model.hpp:77@7496d571]], [[model.hpp:92@7496d571]], [[model.hpp:173-179@7496d571]]) is live deliberate code no gate
     anywhere executes, sitting directly behind the just-landed
     linear crossproduct cache and the GP kernel cache (both
     cache-invalidation sensitive). HIGHEST priority.
  2. BCF x {weights, missing, keepTrees, DART, grouped, linear/GP
     forests, logistic}. Every BCF gate uses one identical minimal
     composition (gaussian, constant leaves, no weights, single
     chain); the glue-precision, per-forest weight transform, and
     calibration code has never run under any other config. BCF is
     the narrowest-composed subsystem in the battery.
  3. Grouped intercepts x zero-weight rows - the SAME
     self-consistency class as the sigma-df defect, unguarded on
     the grouped-tau side (does the group-effect effective count
     exclude zero-weight members?).
  4. Per-observation setPredictor x probit/logistic - the finalize
     step re-derives leaf stats from the family working response
     (PG omega z / probit latent), a path never executed under a
     latent family. Directly relevant to the IRT/bairrtt joint-
     sampler use case.
  5. DART x linear/GP leaf; 6. missing x grouped, missing x BCF;
     7. wide/pooled categorical x {DART, grouped, linear/GP, BCF};
     8. installTrees warm start x {grouped, BCF, linear/GP, missing,
     categorical, multi-chain}; 9. setCutPoints x {categorical,
     missing, grouped, BCF, leaf-covariate cols}; 10. multi-chain
     posterior check x {BCF, DART, linear/GP}.
  Assertion-QUALITY gaps (paths executed but a wrong posterior
  would not be detected): test-dart-mixed-columns.R and
  test-rbart-options.R (DART: dims/simplex-sum only, no signal
  concentration); test-rbart-weighted-binary.R / -weights.R
  (grouped x weights: finite-only); test-bcf.R (blend recovered
  from the same sample = tautology; all real BCF burden on
  bcf-exact.R's single-tree 3-cell kernel); linear/GP probit
  compose blocks (finite-only); thread-invariance tests (bitwise
  across threads, not posterior correctness); change-balance.R
  CONTROL arm (|z|<30 loosened); most state round-trips assert
  serialization not posterior continuation. GOOD tight oracles
  (categorical-exact, logistic-reference pt1, bcf-exact) are all
  tiny single-tree problems - none gates realistic multi-tree scale.
- 2026-07-09: POISON SWEEP complete (16 single-site kernel
  breakages; opus, throwaway worktree, nothing landed; ended clean
  and green). Full detail in the session tmp poison-sweep-results.md.
  Columns: cpp = tests/cpp, tt = tinytest, equiv = equivalence
  statistical verdict, exact = the relevant exact-posterior gate;
  F = caught, P = passed-blind.
    1 bd depth off-by-one ([[model.hpp:1398@7496d571]])        cpp F tt F equiv F  exact change/cat F
    2 bd reverse count wrong tree (moves ~245)    cpp P tt F equiv weak(8/18) change PASS cat F
    3 change FORWARD corr dropped (moves:459)     cpp P tt F equiv P   change F(MIXED arm) cat P
    4 change REVERSE corr dropped (moves:474)     cpp P tt F equiv P   change F(MAIN arm)  cat P
    5 swap validity walk skipped (moves:634)      cpp F(crash) tt F    equiv P change CRASH cat P
    6 sigma df counts all rows (model:1756)       cpp F tt F equiv F(zeroweights)
    7 sigma SSR drops w_i (model:1751)            cpp F tt F equiv F(weighted)
    8 chi-k shape mislabel (model:1729)           cpp F tt P equiv P
    9 chi-k rate drops /leafScale^2 (model:1731)  cpp F tt F equiv F(chik/linear/gp/logistic)
    10 DART counts+1 (model:1680)                 cpp F tt F equiv F(dart)
    11 Polya-Gamma weight omega^2 (model:2180)    cpp F tt F equiv F(logistic) logistic-ref F
    12 grouped precision count (model:2334)       cpp F tt HANG equiv P (no exact gate)
    13 BCF glue drops 1/aVariance (chain:2049)    cpp P tt P equiv P  bcf-exact PASS
    14 BCF forest weight w*m not w*m^2 (chain:2022) cpp F tt F equiv P bcf-exact F
    15 linear ridge drops sigma^2 (model:304)     cpp F tt P equiv F(linear) (no exact gate)
    16 GP nugget drops /w_i (model:721)           cpp F tt P equiv P  (no exact gate)

  SURVIVORS / single-point-of-failure (the findings):
  - Poison 13 BCF a-glue PRIOR PRECISION (1/aVariance) - TRUE
    SURVIVOR, ZERO gates. Verified by orchestrator: [[chain.hpp:2049@7496d571]]
    correctly carries the term (block-6 CONFIRMED); at bcf-exact's
    data size the data precision swamps the prior so dropping it
    shifts E[a*mu] by 0.0001. Needs a prior-dominated bcf-exact mode
    or an a-posterior-VARIANCE check (mean-only is blind).
  - Poison 8 chi-k SHAPE - only one cpp unit test; tt+equiv blind
    (nu=1.5 shifts shape by only +0.5, invisible; the chik
    equivalence scenario df is too small to separate).
  - Poison 16 GP nugget /w_i - only cpp; both R gates blind to
    WEIGHTED-GP (equivalence gp scenario is unweighted; test-gp-
    leaves is finiteness-only and its non-unit-weight fit routes a
    different, unpoisoned path).
  - Poison 12 grouped precision - only cpp catches cleanly; tinytest
    "catch" is merely a CI HANG (extreme weights blow up the ranef
    and stall the tau slice sampler); equivalence fully blind (NO
    grouped/rbart scenario exists). The hang is itself a robustness
    finding -> numerical-robustness review.
  - Poison 15 linear ridge - no exact gate; cpp + the single linear
    equivalence scenario only; tt blind.
  - Poison 2 birth/death reverse count - cpp AND change-balance both
    blind (change-balance targets the change move); only categorical
    -exact, tt reproducibility, and a barely-significant equivalence
    catch it. No detailed-balance oracle for birth/death exists.
  - Confirmed load-bearing: poisons 3/4 (change one-sided
    corrections) are INVISIBLE to equivalence (|z|<=3.4 at 20 seeds
    x1000 draws); only change-balance's MIXED arm (fwd) and MAIN arm
    (rev) catch them - each individually necessary, vindicating the
    mixed-type arm added at the change-move landing.

  METHODOLOGY finding (dev/CI footgun): R CMD INSTALL WITHOUT
  --preclean reuses stale .o and leaves the installed .so SILENTLY
  UNPOISONED after a bartcore header edit - tinytest/equivalence
  falsely passed against a clean .so until --preclean was used.
  Vindicates the standing "--preclean after header edits" rule;
  tinytest/equivalence cannot be trusted after a header edit without
  it. Candidate: a CI build-stamp guard.

  DELIVERABLES: gate-hardening items filed as gate-hardening-1.0
  (TODO). SBC targets = the coverage-matrix uncovered combos above
  (linear/GP x missing highest), carried into review 4. Review 3
  COMPLETE.
