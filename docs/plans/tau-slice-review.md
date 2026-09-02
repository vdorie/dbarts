# Review: tau slice sampler in the in-core grouped sampler (rbart_vi engine path)

Repo dbarts @ bartcore 8d763c9. READ-ONLY review. Files cited file:line.

## VERDICT (summary; detail below)

KEEP the slice sampler as the correct, competently-implemented default, but
there is a clean, exact, drop-in REPLACEMENT for the *default cauchy* prior that
an implementer could land: the Makalic-Schmidt half-Cauchy inverse-gamma scale
mixture, a two-block conjugate Gibbs identical in spirit to the code the repo
already ships for the BCF `a` scalar (chain.hpp drawGlue). It removes the slice
sampler's failure modes (the step-out hang the season just capped), consumes a
bounded 2 RNG draws/sweep, and is the same machinery already trusted elsewhere.
Whether it *improves mixing* is the empirical question measured in section 3.
The gamma prior has NO exact conjugate/GIG draw under its current (gamma-on-sd)
parameterization, so slice must stay for it (or the prior is reconsidered - VD).

---

## 0. What is being sampled (parameterization, established precisely)

Model (design/grouped-random-effects.md, confirmed in code):

    y_i = f(x_i) + b_{g(i)} + offset_i + eps_i
    b_j ~ N(0, tau^2),  j = 1..J          (tau is the SD, NOT the variance)
    tau ~ built-in prior

- tau is the standard deviation of the group intercepts. Confirmed: R preStep
  uses prior precision 1/tau^2 ([[rbart.R:617@4a521760]], `1/state$tau^2`) and draws
  `rnorm(.., sqrt(post.var))`; C++ drawGroupEffects uses `1/(tau*tau)` as the
  prior precision ([[model.hpp:2592@4a521760]]). Reported draws de-scale by sigmaScale()
  ([[chain.hpp:2375@4a521760]]), i.e. original-scale SD.

- Built-in priors, on tau (the SD), original scale with rel.scale = sd(y)
  continuous / 0.5 binary ([[rbart.R:1-6@4a521760]], [[model.hpp:2507-2514@4a521760]],
  design doc:42-46):
    cauchy: dcauchy(tau; 0, 2.5*rel.scale)      -- half-Cauchy on the SD
    gamma:  dgamma(tau; shape 2.5, scale 2.5*rel.scale)  -- gamma on the SD
  The shape 2.5 and the 2.5 scale-multiplier are BOTH fixed constants. The
  cauchy is the default. Internal-scale conversion is one division of the scale
  by sigmaScale() (both are scale families): priorScale_ = 2.5*tauPriorScale/scale
  ([[model.hpp:2626@4a521760]]).

- Conditional target p(tau | b) prop prior(tau) * prod_j N(b_j;0,tau^2):
      logTauPosterior = -J log(tau) - 0.5*sum(b^2)/tau^2 + log p(tau)
  on (0, Inf) ([[model.hpp:2519-2525@4a521760]]). Sampled DIRECTLY on tau (no transform),
  so no Jacobian is needed or present - correct. The R loop's posteriorClosure
  is the identical expression ([[rbart.R:870-876@4a521760]]), so the two paths target the
  same conditional.

The R loop (custom-prior path) and the in-core path are DIFFERENT samplers of
the same target (statistically equivalent, not bit-identical, by design):
- R sliceSample ([[sliceSample.R:46-312@8ee812dc]]): adaptive - L-BFGS-B mode-find + numeric
  Hessian EACH call, width from a Gaussian curvature approx, works on the
  linear (exp) density scale, rejection-sampling fallback if the start density
  is tiny. Expensive (an optim per MCMC iteration); that R-loop cost is exactly
  what the in-core path was built to avoid.
- C++ sliceSampleOnce ([[model.hpp:2539-2563@4a521760]]): fixed width = priorScale_, log
  scale, Neal (2003) step-out + shrinkage, both step-outs capped at 1e4/side.

Group-effect draw also differs: R uses the UNWEIGHTED group mean ([[rbart.R:619@4a521760]]);
C++ uses working-weight sums ([[model.hpp:2588-2593@4a521760]]) so Polya-Gamma logistic
weights compose. For unweighted gaussian data the two coincide (design doc:27-30).

---

## 1. CORRECTNESS WALK

Target density ([[model.hpp:2519-2525@4a521760]]): CORRECT.
- Normal likelihood prod_j (2pi tau^2)^{-1/2} exp(-b_j^2/2tau^2) contributes
  -J log tau - 0.5 SS/tau^2 (SS = sum b^2). Matches exactly. Prior added on
  log scale. Support guard `tau<=0 || isinf(tau) -> -HUGE_VAL` ([[model.hpp:2522@4a521760]]).
- cauchy logTauPrior: -log(pi*scale) - log1p((tau/scale)^2) ([[model.hpp:2508-2510@4a521760]]).
  This is the FULL Cauchy normalizer; the missing factor 2 for the half-Cauchy
  is a constant, irrelevant to slice sampling. log1p is the right numeric choice.
- gamma logTauPrior: 1.5 log tau - tau/scale - 2.5 log scale - lgamma(2.5)
  ([[model.hpp:2512-2513@4a521760]]) = dgamma(tau; shape 2.5, scale), correct and normalized.

Edge cases:
- K=1 (single group): J=1, SS=b_1^2. Target prop p(tau)/tau * exp(-b_1^2/2tau^2).
  Proper under either prior for any b_1 != 0. Fine.
- Tiny K / empty groups: an empty group's b_j is drawn from its prior N(0,tau^2)
  in drawGroupEffects (weightScratch[j]=0 => precision=1/tau^2, mean=0;
  [[model.hpp:2591-2595@4a521760]]) - a legitimate latent, but it has NO mean reversion, so a
  mostly-empty grouping lets tau random-walk upward. This is the documented
  driver of the step-out hang (tau-slice-stepout-cap.md). It is a property
  of the *model/Gibbs kernel*, NOT a bug in the sampler; but it is exactly the
  regime a conjugate replacement (section 4) handles gracefully (a bounded IG
  draw cannot hang).
- All-b-zero / zero-variance: if SS -> 0 the target prop p(tau)/tau^J near 0 is
  non-integrable for J>=1, i.e. improper. b_j are continuous draws so SS>0 a.s.;
  never hit in practice, but it explains why a tiny SS drags tau toward 0 and can
  make step-outs long on the small side. No divide-by-zero: tau starts positive
  and the support guard rejects tau<=0.
- tau -> Inf guard: isinf(tau) scores -HUGE_VAL, so a step-out that reaches a
  huge finite tau still evaluates finite densities; HUGE_VAL upper bound is only
  a clamp, never evaluated as a point. Fine.

Slice mechanics vs Neal (2003) ([[model.hpp:2540-2563@4a521760]]):
- logHeight = logDensity(x) - Exp(1): correct (log of Uniform(0,f(x))).
- Initial bracket left=x-u*width, right=left+width, u~U(0,1): correct random
  placement.
- Step-out both sides with boundary clamp: correct; matches R getInterval
  ([[sliceSample.R:169-196@4a521760]]).
- Shrinkage: propose U(left,right), accept if in slice, else move the side the
  proposal fell on *relative to x* ([[model.hpp:2560@4a521760]]): correct Neal shrinkage.
- Shrinkage cap 1000 iters returning x ([[model.hpp:2556-2562@4a521760]]): a numeric safety
  valve; in exact arithmetic shrinkage always terminates, and returning the
  current point on pathology is a valid (identity) kernel step. Benign.

Step-out cap correctness (the 7c1c7c9 landing, tau-slice-stepout-cap.md):
- The cap is `steps-- > 0 &&` PREPENDED to each while's condition
  ([[model.hpp:2548@4a521760]], [[model.hpp:2552@4a521760]]), short-circuiting before the FP comparison, so any run
  where the cap does not engage evaluates a bit-identical FP sequence. The
  equivalence gate (21/21 identical draws) confirms it never engages in the
  suite. Bias-freeness claim when it DOES engage: "a capped bracket still
  contains the current point, so shrinkage samples correctly inside it"
  ([[model.hpp:2535-2536@4a521760]]).
  - NUANCE (worth recording, not a defect at the shipped cap): this is NOT
    exactly Neal's m-limited step-out. Neal (2003, Fig 3) splits a SINGLE budget
    m randomly across the two sides (J ~ U{0..m-1} left, m-1-J right) precisely
    so the interval's construction probability is symmetric between the current
    and proposed points, preserving detailed balance when the cap binds. Two
    INDEPENDENT per-side caps of 1e4 do not reproduce that randomization, so a
    bound-binding step *can* in principle bias the stationary distribution. In
    practice this is a non-issue: the cap only binds in the runaway regime
    (mode ~ 9e11, tau > ~1e8) where any finite sample is already meaningless and
    the alternative is an indefinite hang; healthy runs measure <100 expansions
    (tau-slice-stepout-cap.md). So the cap is the right safety valve, but
    if one wanted strict Neal correctness the fix is a single randomized budget,
    not two per-side counters. (A conjugate replacement moots this entirely.)

Verdict on correctness: the slice sampler targets the right 1-D conditional with
correct density, no missing Jacobian, sane support/overflow handling, and a
faithful Neal step-out+shrink. The only theoretical blemish is the independent
per-side cap vs Neal's randomized budget, which is immaterial at the shipped 1e4.

## 2. IMPLEMENTATION QUALITY

- Log-scale evaluation throughout the C++ path (good); log1p for the Cauchy tail
  (good); lgamma constant kept so component tests compare to R (design doc:220).
- RNG consumption per tau update, BEFORE the cap: 1 exponential (height) + 1
  uniform (bracket) + (K_left+K_right) uniforms during step-out + (#shrink)
  uniforms. UNBOUNDED before 7c1c7c9; now bounded by the cap but still
  DATA-DEPENDENT and not a fixed count - this is why every grouped fit's RNG
  stream is fragile to any change here (the equivalence gate is bit-exact only
  because the arithmetic is untouched). A conjugate replacement would make the
  per-sweep RNG count a FIXED 2 draws - a real robustness win.
- The fixed width = priorScale_ ([[model.hpp:2671@4a521760]]) is the prior scale, a
  reasonable a-priori guess but NOT adaptive to the posterior sd of tau. When J
  is large the posterior concentrates far tighter than the prior scale, so each
  update pays several shrink steps; when tau has walked large, several step-outs.
  Slice sampling is robust to width mis-tuning (correctness unaffected, only a
  few extra evals), so this is a minor efficiency point, not a defect. The R
  path's per-call mode-find is the opposite extreme (accurate width, huge cost).
- sliceSteps_ = n.thin (design doc:38-40,171-172): tau takes n.thin slice steps
  per update and keeps the last, reproducing the R loop's coupling of the slice
  step count to the thinning interval. Slightly surprising (couples an MCMC
  efficiency knob to a bookkeeping knob) but faithful to the R loop and
  documented.
- Numerical hygiene: no under/overflow traps beyond the guards above; densities
  are finite on (0, huge). The shrinkage 1000-cap and step-out 1e4-cap prevent
  hangs. Clean.

Verdict on implementation: competent, faithful to Neal and to the R loop, with
the documented caps. The one substantive quality gap is the unbounded/data-
dependent RNG consumption (now capped but not fixed), which a conjugate draw
removes.

## 3. MIXING (measured on the shared library)

Two experiments; scripts in this job's tmp dir (mix_grid.R, sampler_compare.R,
exact_check.R, attrib.R). ESS via coda::effectiveSize; IACT = N/ESS.

### 3a. Full rbart_vi tau chain, built-in cauchy, gaussian, n.thin=1, 1 chain,
n.samples=4000, n.burn=1500, n.trees=50, sigma_eps=1 (mix_grid.R):

    K    signal   true_tau  postmean   ESS/4000   IACT   lag1
    3    small(.2)   0.20     0.863       35      115.5  0.692
    10   small(.2)   0.20     0.211      110       36.5  0.775
    40   small(.2)   0.20     0.195      265       15.1  0.774
    3    large(2)    2.00     2.798      735        5.4   0.605
    10   large(2)    2.00     2.669     1335        3.0   0.240
    40   large(2)    2.00     2.669     1884        2.1   0.093

=> tau mixes FINE with a strong group signal (IACT 2-5). It mixes POORLY with a
weak signal, worst at small K (IACT 115 at K=3). So there IS a real mixing
problem, in exactly the weak-signal/small-K corner - the requester's instinct is
not baseless. The question is WHAT causes it.

### 3b. Isolate the SAMPLER (model held fixed). Pure Gaussian variance-component
Gibbs y_ij = b_j + eps, sigma known, group means directly observed; same b-draw
and #sweeps (6000), only the tau kernel varies (sampler_compare.R):

    K   gsd  | slice IACT | exactIG IACT | asis IACT
    3   0.2  |    6.2     |     5.9      |   2.2
    10  0.2  |   23.4     |    25.6      |   2.1
    40  0.2  |    8.7     |     8.7      |   2.1
    3   2.0  |    3.5     |     1.4      |   1.4
    10  2.0  |    1.6     |     1.1      |   1.2
    40  2.0  |    1.3     |     1.2      |   1.1

=> slice ~ exactIG everywhere (identical in the weak regime; exactIG modestly
better in the strong regime but both already ~1-3). CONCLUSION: replacing the
slice with an EXACT conjugate draw does essentially nothing for mixing in the
regime that matters. The slice is NOT the bottleneck.
=> asis (ASIS interweaving) collapses IACT to ~2 UNIFORMLY: an ~11x gain at K=10
weak, ~3x at K=3 weak. Reparameterization, not a better 1-D draw, is what buys
mixing on the tau-b funnel.

### 3c. Exact-posterior validation (exact_check.R). 1-D quadrature of the
marginal tau posterior (b analytically integrated: ybar_j ~ N(0, tau^2+sigma^2/n_j)):

    K   gsd  | EXACT m/sd/med    | slice           | exactIG         | asis
    3   0.2  | 0.200/0.262/0.127 | 0.197/.244/.130 | 0.211/.296/.131 | 0.198/.254/.127
    10  0.2  | 0.254/0.108/0.242 | 0.260/.108/.247 | 0.254/.106/.242 | 0.255/.109/.243
    3   2.0  | 3.076/1.664/2.661 | 3.063/1.64/2.66 | 3.089/1.72/2.66 | 3.088/1.75/2.65
    10  2.0  | 1.642/0.416/1.569 | 1.647/.420/1.57 | 1.639/.410/1.57 | 1.644/.416/1.57

=> all three kernels reproduce the exact quadrature posterior to MC error. The
slice is correct; my exactIG and asis derivations are correct (same stationary
law). This quadrature is the exact-posterior GATE a replacement must pass.

### 3d. ATTRIBUTION - what actually drives the full-model bottleneck (attrib.R).
Full rbart, weak signal, WITH the strong forest f vs WITHOUT f (f=0, pure
intercept + ranef):

    config                    IACT (full rbart tau)
    K=3  gsd=.2  with f          208.3
    K=3  gsd=.2  NO f              8.3     <- matches isolation (6.2)
    K=10 gsd=.2  with f           51.5
    K=10 gsd=.2  NO f             16.0     <- matches isolation (23.4)

=> Removing f drops IACT ~25x at K=3 and ~3x at K=10. The dominant full-model
bottleneck is FOREST-RANEF CONFOUNDING (f and the group intercepts competing to
explain group-level structure), NOT the tau-b funnel and NOT the slice sampler.
The tau-b funnel (the part ASIS fixes) is the SECONDARY, smaller contributor. The
without-f IACTs match the isolation experiment, confirming the decomposition.

MIXING VERDICT: the slice sampler is not the mixing bottleneck in any regime.
The real weak-signal mixing problem is mostly f-b confounding (out of scope of
"the tau sampler"; needs a forest-ranef joint/interweaving move - a separate,
larger research item), with the tau-b funnel a secondary effect that ASIS could
fix. A drop-in tau-sampler swap changes mixing by ~0.

## 4. REPLACEMENT CANDIDATES (derived, this codebase's exact conventions)

### (a) Half-Cauchy via inverse-gamma scale mixture (Makalic-Schmidt 2016) -- EXACT, drop-in for the CAUCHY prior
Prior (default): tau ~ C+(0, A), A = priorScale_ = 2.5*tauPriorScale/sigmaScale
([[model.hpp:2626@4a521760]]; the cauchy scale is exactly this priorScale_, [[model.hpp:2508-2510@4a521760]]).
Auxiliary representation (derivation checked against quadrature in 3c):
    tau^2 | xi ~ IG(1/2, 1/xi),   xi ~ IG(1/2, 1/A^2)   ==>  tau ~ C+(0,A)
Full conditionals (b_j ~ N(0,tau^2), SS = sum b^2, J groups):
    xi   | tau     ~ IG(1,       1/tau^2 + 1/A^2)
    tau^2| b, xi   ~ IG((J+1)/2, 0.5*SS + 1/xi)
Both are exact iid draws via 1/Gamma, EXACTLY the codebase's existing idiom for
the BCF `a` scalar ([[chain.hpp:2290@4a521760]]: `1.0/ext_rng_simulateGamma(rng, shape, 1/rate)`).
Engine sketch replacing [[model.hpp:2670-2672@4a521760]] (ext_rng_simulateGamma takes SCALE):
    double A = priorScale_;                       // cauchy scale, internal
    double xi = 1.0/ext_rng_simulateGamma(rng, 1.0,
                    1.0/(1.0/(tau_*tau_) + 1.0/(A*A)));
    double tau2 = 1.0/ext_rng_simulateGamma(rng, 0.5*(numGroups+1.0),
                    1.0/(0.5*sumSquaredEffects + 1.0/xi));
    tau_ = std::sqrt(tau2);
Properties: FIXED 2 RNG draws/sweep (vs the slice's data-dependent, now-capped
count); cannot hang (removes the entire step-out-cap failure class the season
just patched, tau-slice-stepout-cap.md); faster (no density evals). xi need NOT
be persisted in ChainStateData - it is conditionally independent of everything
given tau, so restore tau then redraw xi|tau; bitwise-continuation is preserved
with no state-format change. MIXING IMPACT: ~none (3b: exactIG ~ slice).

### (b) The gamma prior -- NO exact conjugate/GIG draw as parameterized
The code's gamma prior is dgamma(tau; shape 2.5, scale) -- gamma on the SD tau
([[model.hpp:2512@4a521760]], [[rbart.R:3-4@4a521760]]), NOT the variance, NOT the precision. Posterior:
    p(tau|b) prop tau^{1.5-J} exp(-tau/scale - 0.5*SS/tau^2).
This is NOT GIG: GIG(p,a,b) prop x^{p-1} exp(-0.5(a x + b/x)) needs a 1/x term;
here the likelihood gives 1/tau^2, and substituting u=tau^2 turns -tau/scale into
-sqrt(u)/scale (neither u nor 1/u). So the src/external GIG sampler
(ext_rng_simulateGeneralizedInverseGaussian) does NOT apply. No conjugacy either.
=> Keep the slice sampler for the gamma prior. (Two model-CHANGING alternatives,
VD only: gamma on the VARIANCE tau^2 makes the posterior GIG(2.5-J/2, chi=SS,
psi=2/scale) -> exact GIG draw; gamma on the PRECISION 1/tau^2 is conjugate
inverse-gamma. Either changes what rbart.priors$gamma means to users.)

### (c) ASIS / parameter-expansion interweaving on the (tau, b) ridge
Yu-Meng (2011) global interweaving, the bcf-ridge-interweaving.md precedent
applied to the ranef scale. After the centered draws (b|tau conjugate, then the
CP-sufficient tau|b draw), add a non-centered ANCILLARY draw: eta = b/tau (fixed,
eta_j ~ N(0,1)); tau is then a coefficient with half-Cauchy prior via
tau|v ~ N(0,v) trunc>0, v ~ IG(1/2, A^2/2) (the drawGlue mixture, [[chain.hpp:2285-2290@4a521760]]):
    v      | tau      ~ IG(1, (tau^2 + A^2)/2)
    tau    | eta,y,v  ~ N(m, 1/prec) truncated to (0,Inf),
        prec = sum_j n_j eta_j^2 / sigma^2 + 1/v,
        m    = (sum_j n_j ybar_j eta_j / sigma^2) / prec;   then b <- tau*eta.
Validated against quadrature (3c). MIXING IMPACT: large ON THE tau-b FUNNEL (3b:
IACT 23->2). BUT in the full engine the tau-b funnel is the SECONDARY mode (3d);
the dominant weak-signal bottleneck is forest-ranef confounding, which ASIS on
(tau,b) does NOT touch. So ASIS is a genuine but PARTIAL win whose full-model
payoff must be measured in a prototype before committing (my measurement is the
isolation, not the engine). A "b-ancillary" draw also needs the group means in
the ranef block; more engine surgery than (a).

### What BART-adjacent tools do
rstanarm/brms marginalize the ranef via HMC (not applicable to a Gibbs-in-C
BART). Classic BayesTree-era rbart used the same slice (this is its descendant).
For this exact 1-D half-Cauchy-scale conditional the literature's "smarter" move
is precisely the M&S IG mixture (a) for an exact draw and ASIS/PX (c) for the
funnel - both already precedented in THIS repo (drawGlue, bcf-ridge-interweaving).
Nothing more exotic is warranted for a 1-D conditional.

## 5. MIGRATION COST (any replacement is draw-changing for grouped fits)

- The 20 ungrouped scenarios must stay bit-identical (they do not touch this
  code); the grouped/grouped_aft scenarios' RNG streams CHANGE and their anchors
  re-record. Gate battery a replacement needs:
  1. R CMD INSTALL --preclean; tests/cpp rebuild + ./test_bartcore (testGroupedMath
     already checks priors/posterior/moments - extend it: the new kernel's moments
     vs the R quadrature in exact_check.R, and vs the slice's, on fixed rng).
  2. tinytest full suite; regenerate the grouped hardcoded snapshots by replaying
     whole files (last-ulp depends on process history, design doc:225-229).
  3. Equivalence compare: ungrouped 20 "identical draws"; grouped scenarios are a
     STATISTICAL-equivalence verdict (z-tests on tau/ranef/fit posterior
     means+intervals vs the slice baseline, as the grouped landing did:
     design doc:210-219), NOT bitwise.
  4. Exact-posterior check: the 1-group / 2-group quadrature in exact_check.R
     (marginal tau posterior by 1-D integration) - a replacement must match it.
  5. The custom-prior R loop ([[rbart.R:531-696@4a521760]]) is untouched and must keep working
     (a custom prior forcing the cauchy density is the cross-check the landing used).
- (a) exact-IG: smallest surface (swap [[model.hpp:2670-2672@4a521760]]; no state change; keep
  slice for gamma via the existing priorKind_ switch). Draw-changing => gates 1-5.
- (c) ASIS: larger (new ancillary draw needs group sums in the b/tau block, a
  truncated-normal helper, and it is draw-changing) + a PROTOTYPE mixing gate on
  the full engine to confirm the payoff survives forest-ranef confounding.

## VERDICT

KEEP, with one optional FIX and a redirection of the real mixing lever.

- CORRECTNESS: the slice sampler is correct - right conditional, no missing
  Jacobian (samples tau directly), sound support/overflow guards, faithful Neal
  step-out+shrink, and it reproduces the exact quadrature posterior (3c). One
  academic blemish: the step-out cap is two INDEPENDENT per-side counters, not
  Neal's single randomized budget, so it is not provably unbiased WHEN THE CAP
  BINDS - immaterial at the shipped 1e4 (binds only in a runaway regime where any
  sample is meaningless, and a conjugate draw would moot it).
- IMPLEMENTATION: competent; the only real gap is data-dependent (now-capped) RNG
  consumption, which a conjugate draw would fix to a constant 2 draws/sweep.
- MIXING (the requester's hypothesis): tested honestly both ways. The slice is
  NOT the bottleneck - an exact conjugate draw mixes identically (3b) and hits the
  same posterior (3c). The weak-signal mixing problem is real but is mostly
  FOREST-RANEF CONFOUNDING (3d: ~25x of the IACT at K=3), which no tau-sampler
  change addresses; the tau-b funnel that ASIS fixes is secondary. So the instinct
  "the slice is ripe for replacement to fix mixing" is, on the evidence, wrong; a
  drop-in swap changes mixing by ~0.

RECOMMENDATION
- Implementer-ready NOW (no model decision; draw-changing, needs gates 1-5):
  swap the CAUCHY slice for the Makalic-Schmidt exact IG two-block (4a). It is the
  same trusted idiom as drawGlue, deletes the step-out-hang failure class, fixes
  RNG to a constant count, and is faster - but it does NOT improve mixing and the
  hang is already capped, so its value is robustness/cleanliness, weighed against
  a real grouped-snapshot re-record. DEFENSIBLE BUT OPTIONAL; keep the slice for
  the gamma prior (no exact draw exists as parameterized, 4b).
- Needs VD SIGN-OFF (research, draw-changing): if mixing on grouped fits is a
  priority, the lever is NOT the tau draw. Rank: (1) a forest-ranef joint/
  interweaving move (biggest full-model payoff, the real bottleneck, separate
  item); (2) ASIS on (tau,b) per 4c/bcf-ridge-interweaving (validated ~11x on the
  tau-b funnel in isolation, but partial in the full engine - gate with a
  prototype). exact-IG (4a) is a prerequisite building block for (2).
- If nothing else: the strict-Neal single-randomized-budget step-out cap is a
  ~5-line correctness tidy, but it is cosmetic given where the cap binds.
