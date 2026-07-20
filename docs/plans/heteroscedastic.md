# heteroscedastic

agent: opus (new leaf numerics, the two-leaf-type Chain change, the exact gate;
  C4's R surface rides along - routing per step).
rng: C1 neutral; C2 neutral for the homoscedastic default (posterior-changing only
  on the not-yet-reachable variance path); C3/C4 posterior-changing (variance
  forest becomes user-reachable).
budget: ~1950-3000 lines total across four commits (design section 10); C2 is the
  large structural one.

## Goal

A second, variance-modeling forest: y = f(x) + s(x) eps, s^2(x) a product of
scaled-inverse-chi-squared leaves (HBART, param A). The variance forest is
CONJUGATE - it reuses the existing conjugate move, adds no MoveStrategy, and does
NOT touch monotone's seams. `variance = ~x1 + x2` on dbarts()/bart2(). Homoscedastic
fits stay byte-identical at every commit.

## Context

docs/design/heteroscedastic.md is the spec (VD resolutions inline: param A, weight
channel, nullable variance forest; two riders in section 2). Anchors, re-verified
2026-07-19:
- Chain<L> holds vector<Forest<L>> with one shared leaf type (chain.hpp ~286); the
  variance leaf differs -> a nullable distinctly-typed member (design section 6).
- Conjugate move + per-leaf marginal reused across leaf types: moves.hpp ~61-78;
  metropolisJumpForTree chain.hpp:702. Additive roll rollTreeResidual chain.hpp:713
  is NOT reused (the variance roll is multiplicative, divide by s^2_{-j}).
- Weight channel: response_->workingWeights(), pulled at chain.hpp:783 (latent
  families own it -> gaussian-only refusal). LinearGaussianLeaf own-span suffstat
  accumulation (model.hpp:828-855) is the scale-leaf precedent.
- chi^-2 prior + draw: ChiSquaredScalePrior (model.hpp:2287-2307; posteriorScale =
  df*scale + weighted SSR at :2304); ext_rng_simulateChiSquared random.h:137.

## Constraints

- Homoscedastic (no variance forest) fits BYTE-IDENTICAL at every commit: the
  variance forest is a nullable member behind if(varianceForest_); the mean-forest
  hot path and RNG stream unperturbed. Equivalence trio bitwise + suite at every
  commit; bench-sampler neutral (quiet-machine grant).
- Variance suffstat carries user weights: sum_i w_i e_i^2 / s^2_{-j} (design section
  3). Construction-time refusal: variance forest requires ResponseFamily::gaussian,
  refuses probit/logistic/nbinom/ordinal (workingWeights collision).
- Rider (i): the roll component test must perturb the tree's OWN leaf and assert the
  suffstat does not move (perturbing another tree scales both sides). Rider (ii):
  the m'=2 gate's non-identified a_c d_c ridge may mix slowly - budget long runs; keep
  the quadrature reference finer than the sampler MCse (the monotone-gate lesson).
- No dbarts.h change expected (R-surface + bridge, like monotone); confirm at C4.
  v1 scope: constant scale leaves only; no linear/gp variance leaves.

## Steps

1. C1 - the scale leaf in isolation, opus. ScaleLeafModel concept +
   ConstantVarianceLeaf: own-span scale suffstat, chi^-2 marginal + draw, section-3.4
   calibration. NO Chain wiring. Gate: tests/cpp vs quadrature (marginal, draw
   moments); equivalence trio bitwise; suite. ~350-500 lines.
2. C2 - Chain integration + neutrality, opus (the large, risk-bearing commit). The
   nullable variance forest, the multiplicative roll, the variance sweep, the
   weight-channel coupling (two n-scratch) + gaussian-only refusal. Gate: equivalence
   trio bitwise (homoscedastic unchanged) + a self-identical heteroscedastic fixture +
   tests/cpp incl. the divisor guard (rider i). ~800-1300 lines.
3. C3 - reporting + gates, opus. s(x) train/test channels (NEW separately-typed
   reporting), predict for both forests, scale-leaf state (NEW FlatNode branch),
   benchmarks/R/heteroscedastic-exact.R: (a) m'=1 homoscedastic reduction, (b) the
   m'=2 CLOSING exact gate (two variance trees, one binary predictor, per-cell 2-D
   quadrature on identified s(x)), (c) SBC. Gate: all pass; equivalence trio bitwise.
   ~400-600 lines.
4. C4 - the R surface, opus. variance/n.trees.variance/.variance hypers on
   dbarts()/bart2(); bridge SamplerOptions fields + refusal message; fit class +
   generics; tinytest incl. the recovery test. Confirm no dbarts.h diff. air format.
   ~400-600 lines.

## Verification

Per commit: R CMD INSTALL --preclean; tests/cpp from make clean; equivalence trio
bitwise (equivalence-f494156, bcf-99205ee, multinomial-8c2b5fc); full tinytest. C3:
the exact/reduction gates pass (reference finer than MCse). C4: air format --check .;
a bench-sampler compare confirms the homoscedastic path is speed-neutral (grant).

## Landing

All four commits landed 2026-07-19/20, each independently battery-verified and its
CI sanitizer watched to green before the next (the sanitizer discipline the monotone
arc taught: the plain battery does not catch heap/UB, and R-reachability is where a
memory bug hides):

1. C1 3775437 - the conjugate scale leaf (ScaleLeafModel + ConstantVarianceLeaf):
   closed-form chi^-2 marginal, scaled-inverse-chi-squared draw, section-3.4
   calibration collapsing to the sigma prior at one tree; component-tested vs
   quadrature. No sampler instantiates it yet.
2. C2 dafe96b - Chain integration: a nullable distinctly-typed VarianceForest, the
   multiplicative residual roll (divide by s^2_{-j}, the divisor guard excluding the
   tree's own factor), the weight-channel coupling (w_i^mean = w_i / s^2, sigma fixed
   at 1), gaussian-only refusal, and the five move templates relaxed to
   MoveScorableLeafModel with byte-identical existing-leaf codegen.
3. C3 043d1ba - reporting (train/test s(x), a separately-typed channel), predict,
   and by-name state serialization with a strictly-positive scale-leaf validation.
4. C4 994ec7e - the R surface (variance = ~x / TRUE / selector; n.trees.variance;
   base/power.variance; s.train/s.test; predict returns the variance surface), the
   bridge route (control attribute -> options; append-only variance state block; NO
   dbarts.h change), and the two-part exact-posterior gate.

The conjugacy premise-correction (the TODO/forest-combiner.md "non-integrable,
non-conjugate" framing was wrong; the variance forest is conjugate and reuses the
existing move) was verified against HBART eq. 6-7 by two independent critiques and is
recorded in section 0 and 12; the TODO and forest-combiner.md were corrected at the
design landing.

Correctness: the m'=2 closing exact gate (benchmarks/R/heteroscedastic-exact.R)
matches the sampler to nested quadrature on the identified s(x) - m'=1 reduction gap
0.0010 (tol 0.05), m'=2 s^2 gap 0.0027 (tol 0.045), f gap 0.0012 (tol 0.015), gaps
shrinking with sample count (MC, not bias) - proving the multiplicative roll. Homo-
scedastic fits bitwise throughout (no equivalence re-record); suite 3304/0; every
commit's CI sanitizer (ASAN+UBSAN + valgrind) green, incl. C4's R-reachable paths.

Bench-sampler homoscedastic default-path speed-neutrality CONFIRMED 2026-07-20 (record
docs/plans/monotone-bart.md, one shared same-machine A/B: HEAD carries monotone +
heteroscedastic + the residuals batch, all post-b9d53c7). A naive compare against the
cold b9d53c7 CSV flagged four arms, but the identical b9d53c7 binary re-measured under
the same conditions drifted +7.4% geomean too (zero code change) - pure cross-session
machine drift. Apples-to-apples, HEAD vs a fresh b9d53c7 baseline was geomean 0.981,
every arm within the +/- 6% same-code noise. Draws already bitwise; no incidental
slowdown. v2 doors
(design section 11): non-constant (linear/gp) variance leaves; a k.variance /
variance-df hyper (C4 calibrates from the mean resid.prior); hurdle (the other
multi-forest model). macOS SIP blocks local R-under-ASAN - the CI sanitizer is the
authoritative R-reachable ASAN gate.
