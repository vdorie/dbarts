# nbinom-large-shape-pg: an exact Polya-Gamma draw whose cost does not grow with the shape

Status: RESEARCH, not scheduled (dec-B143 as revised 2026-09-29: for 1.0-0 the cost is documented; whether to build this or admit approximate draws is the post-release approximate-distributions arc)

agent: opus (engine numerics; one C file, one seam edit, tests, one gate arm)
rng: shifting (a negative-binomial row with y + r > 64 draws from a new exact sampler; every recorded snapshot,
baseline and pinned test has shape at most 64 and stays bitwise)
window: post-release, if at all (dec-B143 in [decisions.md](../decisions.md); TODO approximate-distributions)
budget: ~1000 lines (C ~550 with comments, seam ~20, tests/cpp ~220, R harness ~150, gate arm ~40, design
section ~40). Plan estimates have run 1.5-2x low: expect up to ~2000.

## Goal

PG(b, psi) for b > 64 is drawn exactly at a cost bounded in b, so an nbinom fit at r = 1e5 costs about what
r = 64 does. b <= 64 keeps today's integer sum of PG(1) draws bit for bit.

## Context

- Today: [`simulatePolyaGammaShape`](../../src/bartcore/model.hpp) sums round(b) draws of
  [`ext_rng_simulatePolyaGamma`](../../src/external/random.c) (Devroye). Callers: [`NBResponse::drawOmega`](../../src/bartcore/model.hpp)
  and the bridge's nbinom arm of [`drawAugmentation`](../../src/R_interface_bartcore.cpp). Logistic count
  weights ([`LogisticResponse::refreshLatents`](../../src/bartcore/model.hpp), the bridge's logistic arm) and
  multinomial trials ([`MultinomialForestCombiner::drawForestGlue`](../../src/bartcore/combiner.hpp)) sum
  PG(1) in their own loops and are not routed through the seam.
- Where large b comes from: estimated r is capped at 50 by [`NBDispersionPrior`](../../src/bartcore/model.hpp),
  so large shapes come from a large fixed dispersion or from large counts (b = y + r).
- Notation (Windle, Polson and Scott 2014, arXiv:1405.0506, "WPS"): PG(b, psi) = J*(b, z)/4, z = |psi|/2.
  J*(1, z) = sum_n g_n / d_n, g_n ~ Exp(1), d_n = pi^2 (n + 1/2)^2 / 2 + z^2 / 2 (WPS Fact 3). K(t) is the
  cumulant generating function of J*(1, z), xbar = X / b, t(xbar) solves K'(t) = xbar, phi(xbar) = K(t) - t xbar,
  sp_b(xbar) = sqrt(b / 2 pi) K''(t)^-1/2 exp(b phi) is the saddle-point density, f_b is the true density of xbar.

### What the primary source says (read 2026-09-29, not recalled)

- The saddle-point (SP) sampler is approximate. WPS section 5 is titled "An Approximate J*(b, z) Sampler"; its
  algorithm accepts when U k(X) <= sp_b(X), so it draws exactly from the normalized saddle-point density, not
  from PG; section 5.3: it "generates approximate J*(n, z) random variates". BayesLogit 2.4 (GPL >= 3) does the
  same (`PolyaGammaApproxSP::draw` tests against `sp_approx`), and its hybrid `rpg` uses a moment-matched normal
  above b = 170 and a truncated sum of gammas for most non-integer b (its "alternate" sampler is commented out,
  "Need to review"). dec-B14 bars all three.
- The "alternate" sampler (WPS section 4) is exact only if WPS Conjecture 8 holds, checked numerically for
  h in [1, 4]; its acceptance 1/c(h, z) falls as h grows, so large h is a sum of pieces in (1, 4]: linear in b.
- So no published exact PG sampler has cost bounded in b. The design below supplies the missing step.
- Measured in BayesLogit 2.4's SP sampler (20,000 draws per cell): proposals per draw 1.03-1.26, flat in b; but
  at b = 1e5 (any psi), b = 1e4 with psi >= 20 and b = 170 with psi = 80 every draw hits the 200-iteration cap and
  returns a biased value (mean 1.10 times the true mean): its envelope weights are formed with exp() of O(b)
  exponents and overflow. Licence: BayesLogit is GPL (>= 3), dbarts GPL (>= 2); reimplement from WPS, copy
  nothing (its `y_func` also computes `1/3` and `2/15` in integer arithmetic, i.e. as zero).

## Decision

Outcome (2026-09-29): not built for 1.0-0. The help page states the cost
linear in y + r; whether to build this construction or to admit approximate
draws under stated criteria is the post-release approximate-distributions
research arc (TODO), with real dispersion and non-integer logistic weights.


Question: dec-B143 asks for an exact sampler "such as the saddle-point rejection sampler". That sampler is
approximate. Does the ruling stand for the exact construction below, which reuses the SP envelope as a proposal
and accepts against the true density?

Recommendation: yes. Alternatives: (a) WPS SP or BayesLogit's hybrid as shipped - approximate, bars on dec-B14;
(b) the alternate sampler - exact only on a conjecture, still linear in b; (c) a Metropolis step on omega with
an SP proposal - the chain stays exact but the draw is not, and the bridge's replay has no previous omega;
(d) leave cost linear (the rejected alternatives of dec-B143). What would change it: the benchmark (Step 6)
showing one density inversion costs more than ~64 PG(1) draws, which moves the threshold, not the design.

## Design

Sampler for xbar ~ f_b, b > 64, any z >= 0 (then PG = b xbar / 4):

1. Envelope: WPS Proposition 17 as specified in section 5.2 (x_l = m = tanh(z)/z, x_c = 1.1 m, x_r = 1.2 m;
   tangent lines to eta = phi - delta at x_l and x_r). k(x) >= sp_b(x); left piece an inverse-Gaussian
   (mu = 1/sqrt(rho_l), lambda = b) kernel on (0, x_c], right piece a Gamma(b, b rho_r) kernel on (x_c, inf).
2. Stage 1 (cheap, as WPS): draw X ~ k, U1 ~ U(0, 1); continue only if U1 k(X) <= sp_b(X).
3. Stage 2 (new): U2 ~ U(0, 1); accept iff U2 C_b sp_b(X) <= f_b(X), else go to 1.
   Net acceptance f_b/(C_b k): X ~ f_b exactly, given f_b <= C_b sp_b <= C_b k.

Density bound (derived here). By inversion on any line Re s = t below d_0,
f_b(xbar) = (b / 2 pi) int exp(b [K(t + iy) - t xbar - i y xbar]) dy, and |exp(K(t + iy) - K(t))| =
prod_n (1 + y^2 / a_n^2)^-1/2 with a_n = d_n - t > 0. Since prod (1 + w_n) >= 1 + sum w_n and
K''(t) = sum_n 1 / a_n^2, the integrand is bounded by exp(b phi) (1 + y^2 K'')^(-b/2), whose integral is
K''^-1/2 sqrt(pi) Gamma((b - 1)/2) / Gamma(b/2). Hence, for every b > 1, xbar and z,

    f_b / sp_b <= C_b = sqrt(b/2) Gamma((b - 1)/2) / Gamma(b/2) ~ 1 + 3/(4b)   (1.0117 at b = 65).

The bound holds on any line t, so root-finding error in t costs acceptance, never exactness. Measured by
quadrature (b in 4..1e5, z in 0..40): f_b/sp_b lies in [1 - 1/(12b), 1] to within rounding, so stage 2
accepts about 1/C_b - 1/(12 b) of its arrivals: about one density evaluation per draw.

Exactness rests on three things: the bound above (proved); k >= sp_b (WPS Proposition 17 taking alpha_l,
alpha_r at x_c, which needs WPS Conjecture 15: K''/x^3 decreasing and K''/x^2 increasing in x); and deciding
stage 2 correctly. Conjecture 15 depends on neither b nor z (K'' as a function of xbar is z-free). Since x_c = 1.1 m
lies in (0, 1.1], all four monotonicities are used (both ratios, on both branches). In closed form, with
x = tanh(s)/s below 1 and tan(s)/s above, K''/x^3 = coth^2 s - s cosh s / sinh^3 s and
K''/x^2 = 1 + cot^2 s - cot(s)/s. All four were checked here on a 4000-point grid per branch, with no
violation; the fine-grid violations were rounding noise. Step 2 closes it by proof or by a tabulated rigorous
lower bound, so exactness hangs on no numerical conjecture. For stage 2, f_b is evaluated by the trapezoid rule on the vertical
line through the saddle point with a certified error bound (analytic in the strip |Im y| < a_0 = d_0 - t; the
same product bound bounds the strip integrals and the truncated tails). Accept if U2 C_b sp < f - err, reject
if U2 C_b sp > f + err, otherwise halve the step and retry (probability ~1e-30). Measured: 33 nodes (17
evaluations by conjugate symmetry) at spacing 0.5 sd over +-8 sd reach the 1e-11..1e-9 rounding floor for b
from 64 to 1e5, all z. At z = 0 the saddle point is v = 0, where the series branch applies.

Numerical hazards, each a named requirement:
- Saddle equation tan(sqrt v)/sqrt v = xbar (tanh for v < 0): monotone. Use a bracketed Newton in s = sqrt(|v|),
  with the Taylor branch 1 + v/3 + 2v^2/15 + 17v^3/315 for |v| < 1e-3 (floating-point constants). xbar -> 0
  gives v ~ -1/xbar^2; xbar -> inf gives s -> pi/2 (use s = pi/2 - 2/(pi xbar) as the start).
- Overflow: never form cosh^b, exp(b phi) or the envelope weights; carry logs, use log cosh z = z + log1p(e^-2z)
  - log 2, the inverse-Gaussian CDF as a log-sum-exp of the two log Phi terms (its e^(2 lambda / mu) factor
  overflows at lambda = b), and the gamma tail as log Rf_pchisq(2 b rho_r x_c, 2b, lower = 0, log = 1). This
  is the failure BayesLogit shows.
- Cancellation: b [K(t + iy) - K(t)] loses digits as b grows (1e-9 relative at b = 1e5 in the naive form).
  Form log(cos w / cos w0) as log1p(cos delta - 1 - tan(w0) sin delta), delta = 2iy/(w + w0), not as a
  difference of logs. Complex arithmetic by hand in real sin/cos/sinh/cosh (no C99 complex.h dependence).
- Proposal pieces: the right piece is a gamma tail beyond x_c, far past the mode for large b, so it needs an
  exponential-proposal tail sampler, never draw-and-reject. The left piece is a general inverse-Gaussian
  (Michael-Schucany-Haas) truncated at x_c; draw-and-reject is fine there, the mass above x_c being
  exponentially small. The static `simulateTruncatedInverseGaussian` is lambda = 1, t = 0.64 only; do not
  reuse it.
- Extreme arguments: psi = +-700, b = 1e9 must return finite positive draws (tests/cpp pins).

RNG consumption contract ([RNG architecture](../architecture.md#rng-architecture),
[Threading model](../architecture.md#threading-model)): draws come only from the chain's `ext_rng`, as
uniforms, normals and exponentials, in a count that depends on (b, psi), exactly like Devroye's. The skip rule
for inactive rows is unchanged. There is no lazily built table and no mutable static: constants are
`static const`, so draws stay bitwise identical at any thread count. The quadrature sum is scalar in fixed
order, so no dispatch level can flip an accept
([Reproducibility contract](../architecture.md#reproducibility-contract)).

Placement: a new host-agnostic C file, src/external/randomPolyaGamma.c, beside the Devroye sampler, exporting
`ext_rng_simulatePolyaGammaLargeShape(generator, b, psi)` through src/include/external/random.h. It uses only
the support library's existing math ([`Rf_pchisq`](../../src/include/external/stats.h), `Rf_pnorm5`, lgamma).
The threshold lives at the seam: [`simulatePolyaGammaShape`](../../src/bartcore/model.hpp) keeps the integer sum
for b <= 64 and calls the new primitive above that. Both nbinom callers go through the seam, so both are
covered. The seam's comment promising "an approximate primitive" for fractional b is rewritten (dec-B14).

Threshold 64: estimated r is at most 50, so an estimated-r fit with counts up to 14 never leaves today's stream.
The largest shape in any recorded artifact is 60 (the equivalence nbinom scenario: max y 10 plus r <= 50), and
the tinytest nbinom data are Poisson/NB with means <= 6. Measured cost of the old path: ~165 ns per PG(1) (from
3.3 s for 10 sweeps at n = 200, r = 1e4), so b = 64 costs ~10 us, against an estimated 2-4 us for the new path.
Step 6 measures both; if the new path is slower than the sum at b = 64 the threshold rises, and VD may lower it
at the cost of re-recording the nbinom baseline.

What it unblocks (recorded, not in scope): the construction is exact for real b > 1 (C_b is finite; acceptance
1/C_b degrades toward b = 1, C_2 = 1.77). Real b >= 2 splits exactly as an integer sum plus a remainder in [2, 3)
drawn by this sampler, though the quadrature needs many more nodes at small b (tails decay like
exp(-b sqrt|y|)). Real b <= 1 stays open: WPS gives no exact method there. So TODO negbin-real-dispersion is
unblocked for r >= 2 (every row then has b >= 2) and weighted-binary's real weights for w >= 2; r or w below 2
still waits. The same primitive would bound logistic count weights and multinomial trials above 64; that moves
their draws and needs its own ruling (dec-B143 names nbinom only).

## Constraints

- dec-B14: every draw exact; no normal, truncated-gamma or saddle-point-only fallback anywhere, even for
  "extreme" arguments. A failure is a thrown error, never a biased return. No iteration cap that returns.
- b <= 64 bitwise unchanged; the logistic, multinomial and bridge-logistic loops untouched.
- No new dependency; nothing from BayesLogit's source.
- Out of scope: real dispersion, real weights, routing logistic/multinomial, any R surface change.

## Steps

1. Primitive helpers in randomPolyaGamma.c: stable log cosh, the saddle solver, K, K'', phi, and the certified
   trapezoid density with its error bound. tests/cpp: the density agrees with the alternating series (WPS eq. 11)
   to 1e-8 where that series is still accurate (b <= 60 near the mode); it integrates to 1 with mean b m.
2. Close WPS Conjecture 15 (proof of the two one-variable inequalities, or a rigorous tabulated lower bound for
   alpha_l, alpha_r over x_c), recorded in the design section.
3. Envelope and proposal (log-space weights, truncated IG, gamma tail sampler) plus the two-stage loop;
   `ext_rng_simulatePolyaGammaLargeShape`.
4. Seam: threshold 64 in [`simulatePolyaGammaShape`](../../src/bartcore/model.hpp); rewrite its comment.
5. Tests: extend [`testNBPolyaGammaShapeMoments`](../../tests/cpp/test_model.cpp) (b = 1..64 still equals the
   sum bit for bit at seeded streams; b in {65, 100, 1e3, 1e4, 1e5, 1e6} x psi in {0, 1, 4, 20, 80}: mean
   b tanh(psi/2)/(2 psi) within 5 SE and variance b (sinh psi - psi)/(4 psi^3 cosh^2(psi/2)) within 3%);
   extreme-argument finiteness; sanitizers ([Gate hygiene](README.md#gate-hygiene)).
6. benchmarks/R/pg-large-shape.R (new): distributional check and benchmark table, below.
7. [`rFixed`](../../benchmarks/R/negbin-exact.R): add a fixed-r arm at r = 200, so every row runs the new path.
8. Design section in [negative-binomial.md](../design/negative-binomial.md): the construction, the bound, the
   threshold.

## Verification

- `cd tests/cpp && make && ./test_bartcore`: all pass, including step 1 and 5 pins.
- Distributional (Step 6): through the exported [`dbartsDrawLatents`](../../R/augmentation.R), which needs no new
  entry point. `family = "nbinom", y = 0, dispersion = b` draws the new sampler; `family = "logistic",
  weights = b` still draws the sum of b PG(1). For b in {65, 100, 250, 1000} x psi in {0, 0.5, 2, 8, 30},
  1e5 draws each, two-sample KS: no p below 0.05/20; moment z-scores for mean, variance and third cumulant
  within 4; at b in {1e4, 1e5, 1e6} moments only (the sum is too slow). Mutation check: replace f_b by
  sp_b in stage 2 (the WPS approximate sampler) and show the third-cumulant check fails at b = 65.
- Benchmark table (quiet machine; the 2026-09-29 host was at load 8-10, so no timings were taken for this
  plan): seconds per 1e5 draws, old vs new, b in {1, 10, 64, 65, 100, 1e3, 1e4, 1e5, 1e6} x psi in {0, 2, 30};
  and nbinom fits, n = 200, 10 sweeps, r in {50, 1e4, 1e5} (today 0.02, 3.3, 36 s). Pass: new-path time flat
  within 2x from b = 65 to 1e6; the r = 1e5 fit within 2x of the r = 64 fit.
- Fit level: `negbin-exact.R` including the r = 200 arm, within its tolerances (runs in
  [exact-gates.yaml](../../.github/workflows/exact-gates.yaml)).
- RNG class shifting ([RNG classes and their gates](README.md#rng-classes-and-their-gates)), but no recorded
  draw moves: the equivalence compare against the current MANIFEST baseline reports identical draws for every
  scenario, nbinom ([`fitViaNbinom`](../../benchmarks/R/equivalence.R)) included (max shape 60); the four
  test-reproducibility files pass unchanged on the reference build (none fits nbinom); full tinytest passes. No
  re-record, no z-mode compare needed; if any of these moves, the threshold is wrong - stop.
- bench-sampler.R compare (maintainer-run): the seam adds one branch per row on the hot path; no measurable
  change at b <= 64.
