# Monotone "leaf" prior: a hybrid birth/death move that never runs a large count

Status: PROPOSED 2026-09-29; the upgrade dec-B149 names; not ruled.

Under monotone(prior = "leaf") every birth/death move needs Z_T0 / Z_T* = m theta
([Counting: algorithm](../plans/monotone-exact-birth-death.md#counting-algorithm)): T* is the finer
tree, c1 and c2 its two new leaves (c1 lower in the order), U the components of T*'s leaf order that hold them,
m = |U|, and theta the probability that c2 immediately follows c1 in a uniform linear extension of U. Counting
theta is exponential in the order's width; measured moves take up to 2.4 s. dec-B149 ships no limit and a
slow-count warning, and names this hybrid as the upgrade: small moves keep Metropolis-Hastings (MH) with the
exact count, large moves switch to an exact move that never counts.

## 1. Sources, read in full

- Goncalves, Latuszynski and Roberts, "Barker's algorithm for Bayesian inference with intractable
  likelihoods", arXiv 1709.07710v1 (the preprint of the Braz. J. Probab. Stat. 31 paper), sections 1-2. Barker's acceptance is
  a_B = pi(y) q(y, x) / (pi(x) q(x, y) + pi(y) q(y, x)) (eq. 3), with a_MH / 2 <= a_B <= a_MH. The two-coin
  algorithm outputs 1 with probability c1 p1 / (c1 p1 + c2 p2): draw C1 ~ Ber(c1 / (c1 + c2)); on 1 flip a
  p1-coin and output 1 on heads; on 0 flip a p2-coin and output 0 on heads; on tails repeat. It needs known
  constants c1, c2 > 0 and coins that are exactly Bernoulli(p1) and Bernoulli(p2), p in [0, 1]; nothing else is
  bounded. The number of loops is Geometric((c1 p1 + c2 p2) / (c1 + c2)).
- Latuszynski and Roberts, "CLTs and asymptotic variance of time-sampled Markov chains", arXiv 1102.2171
  (Methodol. Comput. Appl. Probab. 15, 2013), Corollary 1 and Theorem 4, which the paper above cites as its
  Proposition 1. Theorem 4(ii): if a CLT holds for f under MH, it holds under Barker with
  sigma2_MH <= sigma2_B <= 2 sigma2_MH + Var_pi(f). The proof applies Peskun's ordering to
  min(1, R) >= R / (1 + R) >= min(1/2, R/2), the last being the acceptance of the lazy chain
  (1/2) Id + (1/2) P_MH, whose variance Corollary 1 gives.
- Vats, Goncalves, Latuszynski and Roberts, "Efficient Bernoulli factory MCMC for intractable posteriors",
  arXiv 2004.07471v3, sections 2-4. Algorithm 2 restates the two-coin algorithm with local
  (per-pair) bounds allowed. Theorem 1: an acceptance pi(y) q(y, x) / (pi(x) q(x, y) + pi(y) q(y, x) + d(x, y))
  is pi-reversible if and only if d is symmetric. The portkey variant (Algorithm 3, Theorems 2-3) caps the
  expected loops at 1 / (1 - beta) for acceptance at most beta a_B.
- Huber, "Fast perfect sampling from linear extensions", Discrete Math. 306 (2006) 420-428, all sections. The
  Karzanov-Khachiyan chain (hold 1/2, else a uniform adjacent transposition when the order allows) has the
  uniform law on linear extensions. Section 3 builds a bounding chain on the positions' upper bounds, coupled to
  the chain non-Markovianly (Theorem 3). Theorem 5: it bounds a single state after
  (16 / pi^2) n^3 (ln n + ln(2 / eps) / 2) steps with probability at least 1 - eps. Section 4 gives coupling from
  the past for non-Markovian couplings (Theorem 8: exact output), where each doubling level's randomness must be
  reused when the forward pass reruns it. Lemma 10: with doubling, the expected running time is at most
  4.3 n^3 ln n steps. Conditions: any partial order on n elements, relabelled so the identity is an extension.

The plan's statements of these results check out: expected coin draws c / (1 + c theta) <= min(c, 1 / theta)
(from the loop count above with p2 = 1), Barker at least half of MH's acceptance, and Huber's O(n^3 log n).

## 2. The move

### Switch rule

A move is switched when W(U) > B, where W(U) is the sum over U's components C of D(C) |C|, D the number of
down-sets: the layered count's work in down-set x leaf units (the plan's unit). B is an engine constant, 2^22 by
default (section 4).

W depends only on T*'s order and the pair, and T* and the pair are the same whichever tree the chain holds, so
the rule is a function of the unordered pair {T0, T*}. It reads no cache, no leaf value and no clock. A
wall-clock or "is the current tree's count cached" rule would not be symmetric and is ruled out.

W is decided cheaply in most moves. A minimum chain partition of each component (bipartite matching, width w)
brackets D: 2^w <= D <= prod over the chains of (length + 1). The lower bound holds because subsets of a
maximum antichain generate distinct down-sets; both bounds held on all 15,794 measured moves.

- Upper bound times |C|, summed, at most B: MH.
- 2^w |C|, summed, over B: switched.
- Otherwise: run the layered count on T*'s side with a work budget B. If it finishes, the move is MH and uses
  that count; if it passes B, abandon it and switch. The abandoned work is at most B units.

In one component the count on T*'s side counts C* and C0, and D(C0) <= D(C*), so an MH move's count costs at
most about 2 B units. With the hybrid the plan's lazy-count case disappears: a birth inside a merged component
has W(C*) >= D(C0), so it is switched whenever that component is large, and never counts it.

### Barker acceptance

Let r1 be the rest of the birth's ratio (tree prior, proposal, and the touched leaves' integrated likelihood
with the frozen leaves fixed, as the plan's Context states it). The birth's ratio is R = r1 m theta; put
c = r1 m.

- Birth (T0 -> T*): accept with probability c theta / (1 + c theta). Loop: with probability c / (1 + c) flip
  the coin, and accept on heads or loop again on tails; otherwise reject. (In the two-coin notation: constant c
  with the theta coin, constant 1 with a coin that always shows heads.)
- Death (T* -> T0): accept with probability 1 / (1 + c theta), the same c computed on the pair. The coin sits
  on the reject side: with probability c / (1 + c) flip the coin, and reject on heads or loop on tails;
  otherwise accept.
- A birth whose pi(T*, same) = 0 (the empty-cone sentinel) has c = 0 and rejects without a coin. A death whose
  merged cone is empty, pi(T0, same) = 0, has c = infinity and rejects at once instead of looping.
- Free bound: theta <= 1 is the only bound known, and the two-coin algorithm already uses it. Its first loop
  ending without a coin is the plan's shortcut: a birth rejected, or a death accepted, with no coin, with
  probability 1 / (1 + c). The MH shortcut (draw u, decide from m and r1) applies to unswitched moves only, after the switch
  is decided, since the pair's acceptance rule must be fixed before u is read.

### The coin

Draw one uniform linear extension of U exactly and call heads when c2 immediately follows c1: P(heads) = theta
by the identity in the plan. Each component of U gets its own draw by Huber's coupling from the past (as
`huberDraw` in [monotone_ratio.cpp](../../benchmarks/kernels/monotone_ratio.cpp) implements it). With two
components (sizes a and b), riffle the two draws uniformly, choosing c1's component's positions as a uniform
a-subset of the m slots. This is exact because e(C1 + C2) = e(C1) e(C2) C(m, a), with each interleaving equally
likely.

After acceptance the move proceeds as the plan's MH move does: the exact pair redraw after a birth, the merged
leaf's draw after a death. The redraw is conditional on T and the frozen leaves, so it does not depend on which
acceptance rule was used.

## 3. Exactness

- Per pair. Integrating the touched leaves gives the marginal move on (T, same) with the proposal q unchanged.
  MH and Barker each satisfy pi(x) q(x, y) a(x, y) = pi(y) q(y, x) a(y, x) (Vats et al., Theorem 1 with d = 0).
  Detailed balance is checked pair by pair, and the switch assigns each unordered pair one rule for both
  directions, so the mixture over pairs is pi-reversible. The leaf redraw on acceptance is an exact
  conditional draw, as in the plan.
- The two-coin output. It is exactly Bernoulli(a_B) given exact coins (Goncalves et al., section 2).
  - c comes from r1, whose parts the MH move already computes.
  - The coin is exact: by Theorem 8, coupling from the past returns an exact draw from the Karzanov-Khachiyan
    chain's stationary law, which is uniform on extensions.
  - The riffle is exact (above).
  - theta > 0 because T0 is a valid tree (e(C0) >= 1), so the loop ends with probability 1.
- Support. a_B > 0 wherever a_MH > 0, so reachability and irreducibility are unchanged.
- Engine seams. The coin's randomness comes only from the chain's own stream. The switch rule reads nothing the
  chain's history could change.

## 4. Cost

Moves: 15,794 births and deaths in 13 probe fits of the plan's
[monotone-order-size.R](../../benchmarks/R/monotone-order-size.R) pairs mode. These are the plan's eight fits
(1, 5 and 10 trees) and five one-tree fits with 1-2 constrained and 1-3 free axes, the p 4 fit included.
Counts, thetas and Huber draws come from [monotone_ratio.cpp](../../benchmarks/kernels/monotone_ratio.cpp)'s
draw mode, on arm64 macOS under background load, so aggregate times are high; the rows below were retimed
alone.

- Coin validity. Draws matched the exact theta on 14,293 moves with theta < 1 (mean z 0.001, mean z^2 1.004).
  Full-extension chi-square on six random 6-7 element orders gave |z| <= 1.12.
- Coins per switched move. c / (1 + c theta) = P(the move ends at T*) / theta, and never more than 1 / theta.
  - A birth pays about its acceptance over theta.
  - A death pays about its rejection over theta: most deaths are rejected, so deaths pay close to 1 / theta.
  - On moves with W > 2^22, theta has quantiles 0.031 (min), 0.045 (10%), 0.083 (median), 0.204 (90%).
- Huber per draw: about 35-45 ns x n^3 ln n in the kernel, 3 ms at 29 elements and 5.5 ms at 36; 12 ms at 45.
  At 128 elements (the mask width) it would be about 0.4 s.

The slowest measured moves:

| fit | U | count on T*'s side | theta | ms per draw | Barker, at most 1 / theta draws |
|---|---|---|---|---|---|
| 1 tree, 1 constrained + 3 free | 36 leaves, D 3.2e6 | 2.38 s | 0.048 | 5.7 | 0.12 s |
| 1 tree, 1 constrained + 1 free | 45, D 1.9e6 | 0.95 s | 0.034 | 12 | 0.35 s |
| 5 trees, 2 constrained + 1 free | 29, D 1.2e6 | 0.89 s | 0.041 | 3.0 | 0.074 s |
| same fit, the death that merges into 3.6e8 | 24 + 30, D 4.6e3 + 1.0e5 | 0.053 s | 0.027 | 5.6 | 0.20 s (not switched) |

The threshold, over all 15,794 moves (count time at ~8 ns per unit on large orders):

| B (units) | moves switched | max unswitched count | decided by the bounds alone | budgeted count needed |
|---|---|---|---|---|
| 2^20 | 20% | 0.08 s | 66% | 34% |
| 2^22 | 11% | 0.23 s | 74% | 26% |
| 2^24 | 5.6% | 0.56 s | 79% | 21% |

At 2^22, the switched moves all sit in few-tree fits. Among the one-tree fits, 100% of moves switch in the
p 4 fit, 83% and 35% in the two 1-constrained + 1-free fits, and 13% across the two 2-constrained + 1-free fits.
The 5-tree, 2-constrained fit switches 15%. Nothing switches in the 10-tree fit, the other 5-tree fits, or the
remaining one-tree fits. The chain-partition bound is loose (median 5.6x, up to 1,500x over D), hence the
budgeted count. A tighter cheap bound would only cut wasted work, not change which moves switch.

- Mixing, switched moves only. Per pair a_MH / 2 <= a_B <= a_MH, so for the birth/death kernel alone
  Latuszynski and Roberts' Theorem 4(ii) bounds the loss: sigma2_MH <= sigma2_hybrid <= 2 sigma2_MH + Var(f).
  For the full sampler (leaf, sigma and other-tree updates between moves) this is not a verified theorem. The
  loss is small where R << 1 (most births: a_B ~ a_MH) and largest near R = 1 (1/2 against 1). In the
  enumerated designs below, with nearly every move switched, switched moves accepted 0.63-0.67 of what MH would
  have. On real fits it is unmeasured; the engine census can log a_B / a_MH per switched move.
- If theta were ever tiny, the portkey two-coin (Vats et al., Algorithm 3) caps the loops at a known acceptance
  cost. The measured minimum, 0.026, does not call for it.

## 5. RNG and reproducibility

- Two-coin: one uniform per loop from the chain's ext_rng, then the coin.
- Huber: each coupling-from-the-past level needs its randomness replayed on the forward pass. Draw a 64-bit
  seed per level from the chain's ext_rng and run the level from a small local generator (splitmix64 or
  xoshiro256**, about 15 lines; the engine has none today). Storing the level's steps instead costs a few MB at
  45 elements. The riffle draws the a-subset from ext_rng.
- A switched move consumes a random number of draws, so B, the base block length and the generator are part of
  the chain's stream. Given the seed and B, a fit reproduces exactly and does not depend on the thread count,
  since the streams are per chain. With B infinite (test hook) no move switches, and the engine matches the
  plan's MH-only move bit for bit. Unconstrained and "joint" fits are untouched.

## 6. What it removes and adds

Removes, under "leaf":

- The moves behind the slow-count warning. A move's count work is bounded by about 2 B units (~0.1-0.2 s at
  2^22) and its memory by the down-sets that fit under B, so the 1-3 s moves and the 424 MB-1 GB layers no
  longer occur.
- The lazy merged-component recount.
- The warning itself, if step 6's given-T prior draw also uses Huber above B. It is exact either way, needs no
  layers, and is off the MCMC path. The interrupt poll stays, moved between coins; the allocation-failure path
  becomes near-unreachable but stays.

Adds, before the plan's 1.5-2x:

- Engine, ~330 lines: the matching bounds (~50), a work budget on the layered count (~30), Huber's sampler with
  its local generator (~130), the riffle and coin (~30), and the two-coin loop in the `logNormalizerRatio` seam
  of `birthOrDeathMove` (~60). The seam then returns a decision instead of a ratio for switched moves. Plus a
  census tally of switched moves and coins (~30).
- Bridge ~30: the B test hook, like the slow-count threshold hook.
- tests/cpp ~220:
  - Huber uniformity by chi-square over all extensions of random 6-7 element orders and a 5-leaf N-plus-chain.
  - Riffle uniformity.
  - The two-coin probability against c theta / (1 + c theta) with a coin of known theta.
  - Switch symmetry: every birth and its reverse death decide the same rule, on random trees, with B at
    several values.
  - Budget determinism: abort at exactly B, independent of any cache.
- Gates ~70: an exact-gate arm running [monotone-exact-enumeration.R](../../benchmarks/R/monotone-exact-enumeration.R)'s
  engine mode with B forced to 0 through the hook, so every move with a component pair switches, and a
  tinytest that the hook changes draws but not the target on a tiny fit.
- Docs ~50.

Total ~700; expect up to ~1,400.

## 7. Brute-force check

[monotone-exact-enumeration.R](../../benchmarks/R/monotone-exact-enumeration.R) `hybrid` mode runs the plan's R
prototype with every move whose U holds 3 or more leaves switched. Each switched move takes the two-coin
acceptance, the coin an exact extension of U from an R port of Huber's sampler, per component and riffled. Pairs
of two leaves keep MH, so both kernels mix in one chain.

- Quick (300k draws, 4 chains), every root group passes:
  - c1: p 0.70, 0.77, 0.99;
  - c2: 0.44, 0.84, 0.67;
  - c3: 0.80, 0.96, 0.99;
  - cN (N-shaped orders, two-component pairs): 0.93.
- 132k-238k moves were switched per design, at 0.34-0.79 coins each.
- `hybrid-nocoin` (coin always heads, dropping theta) fails as it must: c3 p 3.5e-94 and 5.5e-116, cN 7e-77.
- An R chi-square of the port over the 5-leaf N-plus-chain's 7 extensions, 100k draws: p 0.39.

## 8. Open

- B's default. 2^22 keeps every unswitched count under ~0.25 s and switches nothing at 10 or more trees in the
  measured fits. The corrected engine's checkpoint fits (plan, Staging: the checkpoint) should confirm it, and
  give the switched share and a_B / a_MH per fit.
- Whether it goes in before release: the plan's checkpoint decides (dec-B149).
