# multi-forest-models

agent: opus
rng: posterior-changing (new models)
budget: one design note + plan per model, on pickup

## Goal

The multi-forest model family queued behind the Forest split lands one
at a time, each with its own design note: multinomial, heteroscedastic
(HBART), hurdle/zero-inflated.

## Context

- All three were the design's motivating cases for multi-forest
  Samplers (docs/design/core-generalization.md:138-144); none can
  start before forest-split-bcf's refactor.
- Sketches the notes will expand:
  multinomial - K or K-1 forests with Polya-Gamma logistic
  (Held-Holmes style stick-breaking or symmetric formulation; the PG
  sampler ships);
  heteroscedastic - mean forest + variance forest with multiplicative
  positive leaves (Pratola et al.; needs a non-Gaussian leaf prior on
  the variance forest, exercising the leaf-model seam);
  hurdle - binary-occupancy forest + positive-part forest, sharing
  predictors through the data handle.

## Constraints

- Blocked on forest-split-bcf (the refactor half; BCF itself need not
  precede these).
- Each model needs an exact-posterior gate in the established
  single-tree enumeration style before it ships.
- Out of scope here: everything except keeping the queue and its
  blockers recorded.

## Steps

1. On pickup of any one: design note (likelihood, forest coupling,
   prior defaults, surface), VD review, fresh implementation plan.

## Verification

- Per model on landing: exact-posterior gate, component tests,
  recovery tests, bench on single-forest paths unchanged.
