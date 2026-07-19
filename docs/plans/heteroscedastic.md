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
