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
  Samplers (docs/design/core-generalization.md:138-144); none could
  start before forest-split-bcf's refactor, and none could compose
  with a polymorphic coupling seam before forest-combiner's (below).
- Sketches the notes will expand:
  multinomial - LANDED 2026-07-15 (commit bb8855e,
  docs/design/multinomial.md): K symmetric forests coupled through a
  softmax likelihood with an interleaved one-vs-rest Polya-Gamma
  augmentation and a level-centering move;
  heteroscedastic - mean forest + variance forest with multiplicative
  positive leaves (Pratola et al.; needs a non-Gaussian leaf prior on
  the variance forest, exercising the leaf-model seam);
  hurdle - binary-occupancy forest + positive-part forest, sharing
  predictors through the data handle.

## Constraints

- forest-combiner LANDED 2026-07-14 (docs/design/forest-combiner.md):
  the multi-forest coupling seam is a polymorphic ForestCombiner<L>
  hierarchy, BCF its first instance, BCF-touching gates green
  throughout. That design note also names, honestly, what still
  RE-CARVES at the Chain level for these three models and is NOT
  discharged by the combiner refactor: the combined-fit output is one
  n-vector (multinomial's n x K needs signature changes in
  combinedFits/refreshLatents/drawSigma/results.trainingFits); Chain
  holds a single response_/sigma_ (heteroscedastic's variance forest
  needs either the combiner's unused weight-channel route or a
  per-observation sigma, still an open decision); Chain holds a
  single-leaf-type vector<Forest<L>> (hurdle's per-forest response
  families break that invariant); and the state wire format
  (ChainStateData's BCF-shaped a/aVariance/b0/b1) needs a format bump
  for any non-BCF combiner's glue. Each design note below plans
  against this list, not against a combiner that already generalizes
  past it.
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
