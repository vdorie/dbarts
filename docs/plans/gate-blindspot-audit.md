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
     (model.hpp:57,77,92,173-179) is live deliberate code no gate
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
