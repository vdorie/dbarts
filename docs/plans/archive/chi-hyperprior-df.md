# chi-hyperprior-df

agent: opus
rng: posterior-changing (defaults draw-neutral; manual chi() draws change)
budget: small (engine one-liner + R/doc relabels + one test)

## Goal

The `chi(degreesOfFreedom = nu)` hyperprior on the leaf precision `k` samples
an actual chi distribution with `nu` degrees of freedom, matching the
argument and the docs. The empirically-validated binary default is preserved
bit-for-bit under an honest name (`chi(1.5, Inf)`); manual `chi(nu)` users get
the `nu` they asked for.

## Derivation

Leaf model: given `k`, each of the `M` leaf parameters is normal with
precision proportional to `k^2`, contributing `k^M exp(-k^2 sumSq / (2 s_leaf^2))`.
Prior: `k ~ chi(nu, scale s)`, i.e. `k^2 ~ chi-squared(nu) = Gamma(nu/2, 1/2)`,
so the prior contributes shape `nu/2` to `u = k^2` and rate `1/(2 s^2)`. The
`M` normal leaves add `M/2` to the shape and `sumSq/(2 s_leaf^2)` to the rate.
Hence `u = k^2 | . ~ Gamma(shape = 0.5 (M + nu),
rate = 0.5 (sumSq/s_leaf^2 [+ 1/s^2]))`. The correct shape is `0.5 (M + nu)`.
The code used `0.5 (M + 2 nu - 1)` = the exact update for a `chi(2 nu - 1)`
prior: a dropped `|du/dk|` Jacobian in the original `k^2` derivation inflated
the prior's shape contribution from `nu/2` to `nu - 1/2`. The two agree only
at `nu = 1`. Rate terms (leaf, finite-scale `0.5/s^2`, infinite-scale limit)
were already correct and are untouched.

## Fix (two halves, VD-approved)

1. Engine: `ChiKHyperprior::draw` shape `0.5*(M + 2*nu - 1)` -> `0.5*(M + nu)`
   (src/bartcore/model.hpp).
2. Relabel the binary default `chi(1.25, Inf)` -> `chi(1.5, Inf)` at every
   construction and doc site. Because `2*1.25 - 1 = 1.5`, the corrected shape
   at `nu = 1.5` equals the old shape at `nu = 1.25`, so default binary fits
   draw bit-identically. Sites: R/bart.R `.kDefault`, R/model.R
   `resolveNodeHyperprior` and the `chi()` signature default, the dead C
   struct default `kDf` in src/R_interface_bartcore.cpp, and the illustrative
   `chi(1.25)` comments; docs man/bart.Rd, man/dbartsPriors.Rd; NEWS.

## Context

- src/bartcore/model.hpp ChiKHyperprior (shape line).
- [[R/bart.R:269@14bd6b52]] .kDefault; [[R/model.R:347@14bd6b52]] resolveNodeHyperprior, [[R/model.R:451@14bd6b52]] chi().
- correctness-audit.md Block 3 ADJUDICATED FINDING.

## Test

tests/cpp testChiKHyperprior rewritten: 1e5 draws at small `M`, large `nu`
(so the df term is a large share of the shape) pin the sample mean and
variance of `k^2` to the analytic `Gamma(0.5(M+nu), 1/rate)` moments, for
infinite and finite prior scale. This separates `chi(nu)` from the old
`chi(2 nu - 1)` (shape 7 vs 10.5 at the chosen inputs).

## Verification / Status

Gate 1 R CMD INSTALL --preclean .: PASS.
Gate 2 tests/cpp (make && ./test_bartcore): PASS ("all tests passed";
"ok: chi-k hyperprior").
Gate 3 tinytest::test_package: PASS (2473 tests, 0 failures).
Gate 4 equivalence compare equivalence-235bebc.rds: PASS statistically,
exit 0. 16/18 identical draws, including the binary-default `chik`
scenario (draw-neutral). The two manual chi(1.25) scenarios diverge
legitimately: `linear` max |z| = 2.04, `gp` max |z| = 2.60 (both < 4).
Baseline NOT re-recorded.
Gate 5 default draw-neutrality (seeded binary bart2, OLD clean install vs
NEW fixed install): PASS - yhat.train, sampled k, and varcount all
identical().
No tinytest snapshot regenerated: the non-default chi tests (chi(1.25),
chi(1.1)) assert only structure/distribution; chi(1.0) tests sit at the
draw-neutral fixed point.
