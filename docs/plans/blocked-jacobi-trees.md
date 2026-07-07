# blocked-jacobi-trees

agent: opus
rng: neutral at defaults (b = 1 is today's sampler); b > 1 is a new
     exact kernel, gated as posterior-changing
budget: memo + prototype + experiment; implementation planned on the
        experiment's outcome

## Goal

A measured answer to VD's originating question (2026-07-07): can a
forest's trees jump INDEPENDENTLY - instead of sequentially and
conditionally - with a jump on a latent parameter preserving the
correct posterior, and does that admit GPU-level parallelism? The
memo surveys the construction family; the best member gets the
prototype and the ESS-per-second experiment. "No member is sound or
none wins" is a valid, recorded outcome.

## Context

- The family has at least two orderings. Latent-first: draw a latent
  decomposition, then update trees independently given it - the
  noise-splitting augmentation below is the worked-out member, exact
  by construction. Jump-then-correct: propose every tree's move in
  parallel against the current state, then restore the stationary
  distribution afterward (a joint acceptance, a delayed-acceptance
  stage, or a corrective move on an expanded space) - unworked; the
  memo assesses whether any such member is exact and how its mixing
  and sync profile compare.
- Noise-splitting construction: augment with per-tree pseudo-responses y_k ~
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
- Sync profile vs within-chain-threading's data parallelism: m/b
  barriers per sweep instead of ~3m, each worker a cache-coherent
  single-tree update rather than a reduction slice. At b = m this is
  the one BART formulation with enough concurrent work per step to
  interest a GPU (gpu-bart weighs it against other directions).
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

1. Memo (docs/design/blocked-jacobi.md): the construction family
   first - the latent-first augmentation formalized (weights and
   working-response families, sigma and k updates on the recombined
   residual, DART split-count accounting, the bridge draw's cost)
   AND the jump-then-correct orderings assessed for exactness and
   sync profile; the memo names the member the prototype implements.
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
