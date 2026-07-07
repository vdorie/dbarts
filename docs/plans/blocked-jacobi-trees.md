# blocked-jacobi-trees

agent: opus
rng: neutral at defaults (b = 1 is today's sampler); b > 1 is a new
     exact kernel, gated as posterior-changing
budget: memo + prototype + experiment; implementation planned on the
        experiment's outcome

## Goal

A measured decision on exact parallel tree updates via noise-splitting
augmentation, as a candidate mechanism for within-chain parallelism:
does the wall-clock win of updating b trees concurrently survive the
mixing cost of b-fold precision on structure moves?

## Context

- Construction: augment with per-tree pseudo-responses y_k ~
  N(g_k, sigma^2/b) for the b trees of a batch, constrained to sum to
  the batch residual. Marginally the model is unchanged (the
  convolution recovers N(sum g_k, sigma^2)), so Gibbs alternation is
  exact: (1) draw the decomposition - a Gaussian bridge, cheap;
  (2) update the b trees INDEPENDENTLY, each a single-tree update
  against its own y_k with noise variance sigma^2/b. b = 1 is
  sequential backfitting (Gauss-Seidel); b = m is full Jacobi.
- The price: structure moves within a batch see precision b/sigma^2 -
  a move changing fit by delta pays exp(-b delta^2 / 2 sigma^2), so
  birth/death grows conservative with b. ESS per second is the only
  honest metric.
- Sync profile vs the data-parallel alternative
  (within-chain-threading): m/b barriers per sweep (25 at b = 8,
  m = 200) instead of ~3m (~600), each worker doing a cache-coherent
  full single-tree update rather than a slice of a reduction. At
  b = m the same construction is the one BART formulation with enough
  concurrent work per step to interest a GPU.
- Binary/logistic families: the augmentation applies on the working
  (z, w) gaussian representation the tree stage already consumes;
  per-observation weights split the batch precision as w_i * c_k. The
  memo pins this down.

## Constraints

- Exactness is verified, not assumed: the exact-posterior gates
  (benchmarks/R/logistic-reference.R pattern) must pass at b in
  {2, 8}, and b = 1 must reproduce current draws bitwise.
- The experiment precedes any engine integration; a "b > 1 never
  wins" outcome closes this item and is recorded in
  within-chain-threading.md as well.
- Out of scope: DART/grouped interactions beyond noting them in the memo.

## Steps

1. Memo (docs/design/blocked-jacobi.md): the augmentation formalized
   including weights and working-response families, sigma and k
   updates on the recombined residual, DART split-count accounting
   across batch updates, and the bridge draw's cost.
2. Prototype single-threaded: batch update path behind an internal b
   parameter; gates - b = 1 bitwise vs current, exact-posterior at
   b in {2, 8}.
3. Experiment: ESS per second (sigma, k, train-fit summaries) at
   b in {1, 4, 8, 16}, n in {1e4, 1e5}, threaded across the batch;
   compare against within-chain-threading's data-parallel prototype
   on the same hardware.
4. Record the decision in this file and within-chain-threading.md;
   whichever mechanism wins gets the implementation plan.

## Verification

- b = 1 bitwise equivalence; exact-posterior gates at b > 1; the ESS
  table recorded here.
