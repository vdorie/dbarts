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
