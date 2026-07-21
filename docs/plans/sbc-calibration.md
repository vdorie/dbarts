# sbc-calibration

agent: one opus implementer builds the harness and runs it; fable
       adjudicates the rank histograms
rng: neutral (findings only; the tree does not change - SBC drives
     the SHIPPED sampler through its R/C API and checks calibration)
budget: a self-checking harness + a prioritized run; rank histograms
        and the ecdf-diff test recorded here. Pilot FIRST (measure
        per-config wall-clock, validate uniformity on a known-good
        case) before the full matrix.

Review 4 of the retrospective program (retrospective-reviews.md),
prioritized by review 3's uncovered feature combinations.

## Status

DONE. The harness (benchmarks/R/sbc.R) was built, validated against
a known-good baseline, and run across 18 configurations spanning
three tiers: A (gaussian/probit/logistic baselines, DART, grouped
random intercepts, weighted gaussian, BCF), B (linear leaf, all four
review-3 combinations), and C (GP leaf, all three). Every
configuration calibrates. One config initially looked like a real
calibration defect - BCF's glue-on sigma channel (A4, A4b) - but a
follow-up chain-length diagnostic (A4e) resolved it as slow mixing
along the (a, mu-amplitude) scale ridge, not a sampler defect; a
mixing-efficiency remedy is filed as TODO bcf-ridge-interweaving
rather than a correctness fix. Residual SBC coverage gaps (recorded,
not defects): the sampled-k chi-hyperprior channel (no API installs
a prior k draw), DART's alpha-sampling (held fixed for the same
reason), and the grouped model's default half-Cauchy tau prior
(SBC-intractable via the engine's slice sampler; a gamma prior
substituted for A3). Two harness-level findings worth carrying
forward: the matrix (xy) interface drops NA rows even under
missing = "incorporate" (only the formula interface keeps them), and
GP-leaf predict() differs from stored training fits by ~2e-3
(re-kriging jitter; the harness reads the stored fits instead). The
gaussian constant-leaf baseline (~52s at R=200 unloaded) is a
candidate standing CI gate. Tier A alone (baselines through BCF) ran
in ~45 minutes of wall-clock across all runs and diagnostics; tiers
B and C's per-config costs are reported with each result below.

## Goal

Simulation-based calibration (Talts, Betancourt, Simpson, Vehtari,
Gelman 2018): if data are drawn from the prior and the sampler is
correct, the rank of each prior draw among its posterior draws is
uniform. A non-uniform rank histogram is a calibration defect the
exact-posterior gates (tiny enumerable problems) and equivalence
(sameness only) cannot see. SBC is the one gate that checks the
WHOLE posterior of a realistic-scale fit against its own prior.

## Method

Standard SBC, per configuration:
1. Draw theta0 from the model's prior (sigma from the chisq-
   calibrated prior; the regression function from the BART tree +
   leaf prior via sampleTreesFromPrior/sampleNodeParametersFromPrior;
   k from the chi hyperprior when active; grouped tau + effects from
   their prior; BCF glue from its prior).
2. Simulate y | theta0 through the SAME likelihood the sampler
   assumes (gaussian/probit/logistic; weights; offset).
3. Fit with L retained draws (thin to near-independence - use the
   sampler's own autocorrelation to choose thinning; L ~ 100-256).
4. For each scalar FUNCTIONAL, rank theta0's value among the L
   posterior draws. Functionals (trees are too high-dim to rank
   directly - rank low-dim summaries): sigma; f(x*) at a few fixed
   test points; the average f; k when sampled; grouped tau and a
   couple of group effects; BCF (b1-b0) and a at test points.
5. Over R replications (R >= 200 for a usable histogram), test rank
   uniformity: the ecdf-difference band (Talts fig. 1 style) plus a
   chi-square / KS on the binned ranks. Flag any functional whose
   ranks are non-uniform (over/under-dispersed = posterior too
   wide/narrow; sloped = biased; U/inverted-U = mis-scaled).

PILOT FIRST: run the plain gaussian constant-leaf case at small
n/trees, confirm the harness produces uniform ranks (it MUST on the
known-good baseline - a non-uniform baseline means the harness is
wrong, not the sampler), measure wall-clock per replication, and
report the projected cost of the full matrix before running it.

## Prioritized configurations (from review 3's uncovered combos)

Highest first (each routes through posterior code no gate exercises):
1. linear leaf and GP leaf, with and without a MISSING (NA)
   leaf-covariate value - the imputation-at-standardized-mean path
   (model.hpp:173) no gate executes.
2. BCF (glue integrated) - the a-glue prior precision is a true
   gate survivor; SBC of a and (b1-b0) is the natural calibration
   check, and a prior-weak (small-n) config makes the prior term
   matter.
3. grouped intercepts with zero-weight rows and with a non-gaussian
   base (probit/logistic) - grouped-tau self-consistency.
4. DART - variable-selection calibration (does the Dirichlet
   posterior cover the prior correctly).
5. weighted gaussian and weighted GP; the DART 1e-300 probability
   floor at high sparsity (does an alpha/rho near the floor bias
   selection - the block-3 note earmarked for SBC).
Baselines for control: plain gaussian, probit, logistic (each must
calibrate).

## Constraints

- One review at a time; this implementer is the whole fan-out.
  FOREGROUND only, no sub-agents.
- Drives the INSTALLED/shipped sampler only; no engine edits.
- The harness lives under benchmarks/ (a new benchmarks/R/sbc.R or
  a benchmarks/sbc/ dir), not inst/tinytest (too slow for the
  suite); if any config's calibration is clean and cheap enough to
  become a standing gate, note it as a candidate.
- Not a quiet-machine job (calibration is timing-independent), but
  it is long wall-clock - report progress.

## Verification

- Rank histograms + ecdf-diff verdict per functional per config
  recorded here. A non-uniform histogram becomes a fix item under
  the standard gates (with the SBC config as its reproduction).
  Clean configs cheap enough to run in CI become gate candidates.

## Pilot

### Harness design (benchmarks/R/sbc.R)

Drives the shipped sampler through the R API only. Per replication:
1. theta0 from the prior: `sampleTreesFromPrior` + `sampleNodeParametersFromPrior`
   give the forest+leaves via the engine's own machinery; f0 read with
   `predict(x)` / `predict(x*)` (reported scale). sigma0 drawn in R from the
   engine's reported-scale residual prior.
2. y0 simulated through the assumed likelihood (gaussian: f0 + sigma0*N;
   probit: Bernoulli(pnorm(f0))).
3. MCMC re-initialised from a FRESH independent prior draw (not theta0), then
   `setResponse(y0)` (updateScale = FALSE keeps the build scale), then
   `run(burn, L)` with `n.thin` thinning.
4. rank = #{posterior_l < theta0} in {0..L} for sigma, f(x*) at 5 fixed points,
   and avg f.
Uniformity: simulation-based simultaneous ecdf-difference band (Talts fig. 1,
the headline verdict, already multiple-look corrected) + chi-square and KS as
secondary signals. One sampler is reused across reps (its RNG advances).

Self-consistency (the whole game): the data scale is fixed once at build and
never disturbed (updateScale = FALSE), so the prior draw and the posterior
share one internal->reported mapping. The sigma prior on the REPORTED scale is
sigma^2 ~ df * sigest^2 * rawScale / chisq(df) with rawScale =
qchisq(1 - quant, df)/df; the internal (sigest/range)^2 factor and the
reported = internal * range remap cancel range^2, so it is range-independent
(read off R_interface_bartcore.cpp + model.hpp GaussianResponse). This puts
P(sigma < sigest) = quant. sigest is pinned explicitly (`sigma = 1`) so the
prior is NOT data-dependent empirical Bayes.

### Sigma-prior moment check (required before trusting SBC)

P(sigma < sigest) = 0.8998 (target 0.90); median 0.4969 vs theory 0.4970. PASS.
The reconstruction matches the engine's calibrated prior exactly.

### Thinning choice

Unthinned ACF at n=150, nTrees=50: sigma and f(x*) ~0.5 at lag 1, decaying
below 0.1 by lag ~15-30 (realization-dependent). avg f (gaussian) mixes
instantly (~0); avg f (probit) is the slow one (0.87 at lag 1, <0.1 by ~20,
couples to the global latent scale). thin = 30 leaves retained-lag-1 ACF ~
0.04-0.13 for every functional in both families = near-independent. Justified;
the driver runs this diagnostic automatically so each new config re-checks it
(slower-mixing configs must raise thin, which scales cost linearly).

### Baseline verdict (gaussian, constant leaf, n=150 p=3 nTrees=50, R=L=200,
### thin=30, burn=50) -- all PASS

    functional  chisqP   ksP   ecdfDiff  band(.092)
    avg.f       0.136   0.511   0.0564    PASS
    f.star1     0.902   0.999   0.0259    PASS
    f.star2     0.644   0.697   0.0486    PASS
    f.star3     0.760   0.500   0.0563    PASS
    f.star4     0.371   0.973   0.0314    PASS
    f.star5     0.348   0.305   0.0663    PASS
    sigma       0.457   0.882   0.0371    PASS

Rank means all near L/2, histograms flat (no slope/U/inverted-U). The harness
is correct: a known-good sampler yields uniform ranks for every functional.

### Probit (same settings, R=200) -- generalises cleanly

    avg.f  0.444 / f.star2 0.684 / f.star3 0.828 / f.star5 0.549 : PASS
    f.star1 chisqP 0.008 but ecdf-diff 0.045 << band 0.092, KS 0.74 : PASS
    f.star4 chisqP 0.080, ecdf-diff 0.0876 (just inside band)      : PASS

All 6 functionals pass the simultaneous ecdf band. The lone chisqP=0.008 is an
expected multiple-comparison false positive (13 chi-square tests across the two
configs); no systematic bias appears, confirming the probit latent scale is
self-consistent between predict and the posterior.

### Measured cost

Per-sweep wall-clock (single thread, n=150 p=3 nTrees=50):
constant 37.5us, probit 45us, DART 40us, linear-leaf 247us,
GP-leaf(max.leaf.size=64) 1110us. Per rep = 7500 sweeps (=(burn+L)*thin) x
per-sweep. Observed: gaussian 0.262 s/rep -> 52s for R=200; probit 0.312 s/rep
-> 62s. The rep-time = sweeps*us model is accurate (7500*37.5us = 0.28s).

### Projected cost of the full matrix (R=200, thin=30 assumed; re-check
### thinning per config -- it scales cost linearly)

    baseline logistic          ~1 min   (measure; ~probit)
    DART                        ~1 min
    weighted gaussian           ~1 min
    grouped +zero-wt / +probit  ~2-6 min each  (measure; rbart R-callback/sweep)
    BCF + prior-weak BCF        ~2-3 min each  (est ~2.5x const; +glue harness)
    linear leaf +/- NA (x4)     ~6 min each  = ~24 min
    GP leaf +/- NA / weighted   ~28 min each = ~84 min
Full matrix at R=200 ~ 2.5-3 h wall-clock, dominated by GP (~1.5h) and linear
(~0.5h). All single-thread, foreground, timing-independent (share the machine).

### Recommended sequencing (cheap -> expensive, priority-weighted)

1. logistic baseline (~1 min) -- close the control set.
2. DART (~1 min, priority 4) + DART 1e-300 floor probe (priority 5).
3. grouped +zero-weight gaussian and grouped +probit (priority 3) -- measure
   first via rbart / the grouped C API; needs tau + group-effect functionals.
4. BCF + prior-weak small-n BCF (priority 2) -- needs the glue-prior draw and
   a / (b1-b0) functionals; the a-glue precision is the true gate survivor.
5. linear leaf, then linear+NA (priority 1) -- ~24 min; NA needs the
   missing-covariate imputation path (model.hpp:173) wired into the draw.
6. GP leaf, GP+NA, weighted GP (priority 1 & 5) -- run last / overnight (~1.5h);
   cap max.leaf.size to bound the cubic kernel cost.

This sequencing was followed; the results below are organized by tier
(A/B/C) rather than by run order, since tiers B and C map directly onto
steps 5 and 6 above.

### Harness subtleties and extensions

- Prior constructors are reachable only via `dbartsPriors$...` in plain R
  scope (not bare), or unevaluated inside dbarts()'s prior args.
- Binary families need a 0/1 build response; probit predict/test are latent
  scale (confirmed self-consistent by the clean f* calibration).
- ecdf-band is the verdict; 20-bin chisq/KS are multiple-comparison-prone.
- Matrix extensions needed for the full matrix (all since built and
  exercised in tiers A-C below): NA-covariate draw+refit, BCF glue draw +
  a/(b1-b0) functionals, grouped tau/effect draw + functionals, logistic
  sim, weights.
- Candidate standing gate: the gaussian constant-leaf baseline calibrates in
  ~52s/R=200 (~26s at R=100) -- a plausible CI calibration smoke test.

## Results

### Tier A - baselines, DART, grouped intercepts, weighted gaussian, BCF

R=200, L=200, thin=30, burn=50, n=150 p=3 nTrees=50 unless noted.
Verdicts persisted per config.

#### A1. logistic baseline -- ALL PASS

    functional  chisqP   ksP   ecdfDiff  band(.092)
    avg.f       0.125   0.396   0.0631    PASS
    f.star1     0.877   0.946   0.0354    PASS
    f.star2     0.174   0.382   0.0627    PASS
    f.star3     0.772   1.000   0.0220    PASS
    f.star4     0.924   0.884   0.0392    PASS
    f.star5     0.196   0.628   0.0518    PASS

Sim path: y0 ~ Bernoulli(plogis(f0)). Mixes fast (ACF<0.1 by lag ~5), so
thin=30 is conservative. 93s for R=200 (Polya-Gamma sweep pricier than probit).
No sigma functional (fixed by the link). Control set (gaussian/probit/logistic)
now all calibrate.

#### A2. DART variable selection (gaussian, n=200 p=10, R=200) -- CALIBRATES

Self-consistency subtlety: sampleTreesFromPrior grows trees under the CURRENT
split probs (uniform at DART init), NOT a Dirichlet draw. So the harness uses
TWO samplers per rep -- a generator (cgm with fixed split.probs = s0, a genuine
Dirichlet-prior draw) supplies the forest, and a DART fit (alpha fixed,
update.alpha = FALSE, update.delay = 100) whose posterior varprobs must cover
s0. Joint prior s0 ~ Dirichlet(alpha/p), forest | s0 with splits ~ s0 = exactly
what the generator produces and the fit assumes. Dirichlet moment check PASS
(mean 0.100 vs 0.100; var 0.0451 vs 0.0450).

A2a standard (alpha=1): all 10 split-prob functionals s1..s10 PASS (ecdf-diff
0.029-0.082 << band 0.092), sigma PASS, avg f PASS. Two of five f(x*) point
functionals marginally over the band (f.star1 0.0934, f.star5 0.1014 vs 0.092)
with healthy chisqP (0.196, 0.255) and no shared direction -- multiple-
comparison noise across 17 functionals, not a systematic f defect. The DART
target (does the Dirichlet posterior cover the prior) calibrates cleanly. 60s.

A2b sparsity / 1e-300-floor probe (alpha=0.05, floor incidence 0.029 = ~3% of
s0 components pinned at 1e-300): variable selection STILL calibrates -- 9/10
s_j PASS (s2 marginal, ecdf-diff 0.104 vs 0.092, one of 17 = noise), all f*,
sigma, avg f PASS. The floored (near-zero) components show NO rank-0 pileup and
NO bias of the signal components. VERDICT: the 1e-300 floor does not bias
selection at realistic sparsity (the block-3 note's concern). Pushing alpha
lower drives many components to the floor simultaneously, producing
tie-degenerate ranks on already-negligible components -- a numerical tie
artifact, not a sampling defect. 57s.

#### A3. grouped random intercepts (n=160, G=8, R=200) -- tau/b CALIBRATE

rbart_vi's in-core grouped path (bartcore.groups control attr). Tau prior: the
engine offers half-Cauchy or gamma; SBC uses GAMMA (shape 2.5, scale
2.5*rel.scale, rel.scale=0.2) because the half-Cauchy's infinite-variance tail
throws occasional astronomically large tau0 that inflate the response scale and
stall the engine's tau slice sampler (stepping out by fixed width over a range
up to 1e6+) -- brute-force SBC with it is intractable (a real finding: SBC of
the default half-Cauchy tau needs a bounded-support reparam or a smarter slice
width). Group effects b_g ~ N(0, tau^2). Tau moment check PASS (mean 1.249 vs
1.250; var 0.618 vs 0.625). A grouped sampler fixes its response at creation
(setResponse refused), so each rep REBUILDS the fit; groups are assigned
independently of x to keep the f/b partition clean.

A3a grouped GAUSSIAN + 20% zero-weight rows (the review-3 target): tau PASS
(ecdf 0.083 < 0.092), b1 PASS (0.045), b2 PASS (0.061), avg f PASS, all 5 f*
PASS. sigma FLAGS hard (ecdf 0.213, chisqP 0). The sigma flag is a HARNESS
artifact, not a sampler defect: because setResponse is refused, the grouped fit
is rebuilt each rep and the gaussian response-adaptive scaling re-anchors to
range(y0) (~1.7x range(yBuild)), while the prior f0 is drawn at the fixed
generator scale -- the width mismatch is absorbed entirely by the residual
variance. Confirmed benign by the plain baseline (reused sampler, fixed scale:
sigma calibrates) and by A3b below (scale=1: everything calibrates). tau and
b_g are reported-scale and fitScale-independent, so they stay self-consistent
despite the rebuild -- and they are exactly the grouped-tau target. Zero-weight
rows do NOT break grouped-tau/effect calibration. 59s.

A3b grouped PROBIT (non-gaussian base, no response scaling -> fully clean):
b1 PASS (0.063), b2 PASS (0.037), avg f PASS, all 5 f* PASS; tau ecdf 0.0918 vs
band 0.0917 -- a hair over (chisqP 0.220 healthy, ksP 0.055), i.e. a 95%-band-
edge near-PASS within noise, corroborated by the clean b_g that share tau's
prior. grouped-tau self-consistency holds across both bases. 67s.

#### A4. BCF glue (n=200 unless noted) -- final verdict: mixing, not a defect

The initial run (A4) found the a-glue prior precision correctly calibrated
but flagged sigma as non-uniform. A settled rerun with properly-derived
thinning (A4b) confirmed the sigma flag persisted, and at that point it read
as a genuine calibration defect. A fixed-glue control (A4c) localized the
flag to the glue-on path specifically (the two-forest backfit itself is
exact), and a prior-weak configuration (A4d) showed the bias shrinking as n
shrinks - a pattern that favors a mixing explanation over an algebraic error.
A follow-up chain-length diagnostic (A4e) confirmed this directly: the bias
shrinks monotonically as the chain is run longer and crosses into the
passing band, with every control functional passing at every chain length.
FINAL VERDICT (A4e): H-MIX - the glue-scaled forest-update path is
exonerated; the (a, mu-amplitude) scale ridge mixes too slowly for the sigma
channel at practical chain lengths, not a correctness bug. A mixing-
efficiency remedy (interweaving) is filed as TODO bcf-ridge-interweaving,
not a correctness fix. The four sub-experiments below are the investigation
trail that reached this conclusion.

A4 (gaussian, n=200, R=200) -- a-glue PRECISION calibrates; raw a/(b1-b0)
SBC is ILL-POSED by design; sigma flag OPEN at this point

Headline: the a-glue prior precision (the one true gate survivor) IS correctly
calibrated. BCF is an internal bartcore sampler (not R5); numGroups==0 so
setResponse works and one sampler serves all reps at a fixed scale. The
reported<-internal affine map is recovered by regressing one run's combined
fits on a*mu+b_z*tau (R^2 = 1 exactly, fitScale = 5). Glue moment check PASS
(a Cauchy IQR 3.99 vs 4.00; b sd 0.707 vs 0.707). Forests drawn from the engine
prior; a ~ Cauchy(0, sd.control=2), b0,b1 ~ N(0,0.5) drawn in R. Per-sample glue
+ forest fits collected one draw at a time (only current state is exposed).

KEY MODEL FACT (not a sampler bug): a*mu and b_z*tau are invariant under a joint
sign flip (a,mu)->(-a,-mu), and the engine's a prior is symmetric Cauchy, so the
raw a and (b1-b0) posteriors are sign-symmetric / bimodal -- their SBC is
ill-posed. Confirmed: over 15 reps the posterior a-sign matched theta0's sign
only ~50% (6-8/15), raw-a ranks span 0-150. So a / (b1-b0) FLAG hard (ecdf 0.28,
0.21) BY CONSTRUCTION. The a-glue precision must instead be checked through the
IDENTIFIED (sign-invariant) magnitude and functions:

    abs.a       ecdf 0.086  PASS   <- the a-glue precision, calibrated
    eff1..5     ecdf 0.048-0.072  all PASS  <- treatment effect (b1-b0)*tau(x*)
    prog1,3,4,5 ecdf 0.041-0.078  PASS; prog2 0.096 marginal
    abs.diff    ecdf 0.096  marginal (|b1-b0|; |a| mixes slowly, see below)

Verdict: BCF's treatment-effect function and prognostic scale (incl. the a-glue
precision via abs.a) calibrate. 102s (per-sample collection loop). The sigma
flag seen here at thin=30/burn=50 was followed up with the settled rerun
below (A4b).

A4b. BCF settled rerun (thin=120, burn=18000 sweeps, R=200, L=150) --
sigma flag PERSISTS = read at the time as a CALIBRATION DEFECT fix item

Thinning justified from a 4000-draw unthinned chain: |a| ACF 0.37/0.26/0.15/
0.11/0.08 at lags 30/50/80/100/120 (first <0.1 at lag 103); |b1-b0| 0.01 at
lag 120. thin=120 makes the slowest channel near-independent. Burn = 150 outer
draws x thin 120 = 18000 sweeps (a thin=30 burn sweep 50->250->600 outer had
shown sigma's rank centre still drifting: ecdf 0.196 -> 0.122 -> 0.104,
mean 61 -> 67 -> 71 of 75). Run as 4 x 50-rep chunks (~2 min each), ~2.5 s/rep.

    functional   chisqP   ecdfDiff  band(.0924)  verdict
    sigma        0.000     0.1344      FLAG   <- posterior sigma biased HIGH
    a            0.000     0.2520      FLAG   (sign-ill-posed, by construction)
    abs.a        0.119     0.0604      PASS   <- a-glue precision magnitude
    b1.minus.b0  0.000     0.1884      FLAG   (sign-ill-posed, by construction)
    abs.diff     0.657     0.0509      PASS   (clean at thin=120; 0.096 at 30)
    prog1..5     --        0.076-0.109 prog1 (0.099) and prog4 (0.109) marginal
    eff1..5      --        0.051-0.062 all PASS

VERDICT AT THIS POINT: read as a CALIBRATION DEFECT. With the chain settled and
ranks near-independent, sigma's ranks remained non-uniform: rank mean 63.6 vs
75 = theta0's sigma sits low in its posterior = posterior sigma looked
systematically TOO LARGE (over-estimating residual noise), with the
prognostic function marginally affected (prog1/prog4). Magnitude ~8%
rank-mean shift. Repro: `Rscript benchmarks/R/sbc.R bcf 200 150 120` with
burn=150 (runSbcBCF). (This reading was superseded by A4e below.)

A4c. fixed-glue control (update.a = update.b = FALSE, a=1, b0=0, b1=1 in
generator and fit; thin=30, burn=18000 sweeps, R=150) -- ALL PASS

    sigma  ecdf 0.0450 (band 0.105), KS 0.865, rank mean 76.7/75  PASS
    prog1..5 ecdf 0.047-0.076 all PASS; eff1..5 ecdf 0.035-0.085 all PASS

The two-forest backfit WITH FIXED GLUE is exactly calibrated -- both forests'
prior-draw machinery matches the MCMC prior, and sigma is clean. The apparent
defect is therefore IN THE GLUE-ON PATH ONLY. The glue full conditionals were
verified algebraically against the model (a: normal with precision
1/aVariance + sum w mu^2/sigma^2; aVariance: IG(1, (a^2 + scale^2)/2), the
exact t_1 scale-mixture conditional; b0/b1: normals with prior precision
1/bPriorVariance), and k is fixed at 1 for both forests by the calibration map
(no hidden k sampling). Two suspects remained at this point: the (a,
mu-amplitude) SCALE RIDGE -- a*mu is invariant under (a/c, c*mu), leaving a
long ridge the chain must traverse by alternating conjugate a-draws and
per-tree leaf draws, whose relaxation time can far exceed the |a| ACF measured
on one dataset -- or a subtle exactness gap in the glue-scaled forest updates
(formForestResponse's r/m response with w*m^2 weights, |m| floored at 1e-9).
A4e below resolves this in favour of the mixing explanation.

A4d. bcf-weak: prior-dominant small-n BCF (n=40, R=200, L=150, thin=120,
burn=18000 sweeps) -- HEADLINE PASS: the a-glue prior precision calibrates

The poison-sweep's one true gate survivor was the a-glue prior precision
(1/aVariance); a prior-weak fit is where dropping it would show. It does not:

    sigma        chisqP 0.087  ecdf 0.0596  PASS  (rank mean 75.8/75, centred)
    abs.a        chisqP 0.371  ecdf 0.0426  PASS  <- THE headline channel
    abs.diff     chisqP 0.359  ecdf 0.1007  marginal (ksP 0.023; 1 of 15
                 functionals just over the band = multiple-comparison noise)
    prog1..5     ecdf 0.031-0.070  all PASS
    eff1..5      ecdf 0.031-0.066  all PASS
    a, b1.minus.b0  FLAG (sign-ill-posed, by construction)

Glue prior moment check as A4 (same prior code path, n-independent). NOTE the
n-pattern: sigma calibrates at n=40 but flagged at n=200 -- an exactness bug in
the glue draw would show MORE when the prior dominates, so this pattern favours
the scale-ridge mixing explanation (likelihood tighter along the ridge at large
n) over an algebraic error; A4e below tests this directly. ~2.4 min/100 reps.

A4e. chain-length diagnostic (2026-07-10) -- FINAL VERDICT: H-MIX

Three points spanning 8.3x chain length on the standard n=200 glue-on config
(R=200, L=150, burn=150 thinned units, thin 120/360/1000, so longer points
get strictly more burn-in and spacing - conservative for a defect verdict):

    point  thin  sweeps/rep  sigma rank mean  sigma ecdf  chisqP  verdict
    A       120     36k         65.8 / 75       0.1127    0.0002  FLAG
    B       360    108k         67.8 / 75       0.1066    0.105   FLAG(marg)
    C      1000    300k         69.7 / 75       0.0652    0.617   PASS

The bias shrinks monotonically and crosses INTO the band (a defect would
plateau at ecdf ~0.13); every control (abs.a, prog1-5, eff1-5) passes at
every point; raw a and b1-b0 stay sign-ill-posed by design. VERDICT: slow
mixing, and the glue-scaled forest-update path is EXONERATED - the stationary
distribution is correct within SBC resolution (consistent with A4c
fixed-glue-exact and A4d n=40).

Mechanism, measured directly on long unthinned chains (IACT via Geyer): sigma
and the (a, mu-amplitude) ridge coordinate co-relax, and the relaxation time
scales steeply with the prognostic amplitude |a0| ~ Cauchy(0, 2): IACT ~4
sweeps at a0 near 0, ~240-630 at moderate |a0|, ~2500-6600 at the tail (|a0|
~ 7-13; sigma ACF at lag 1000 still 0.48-0.77 there). At thin=120 the
strong-signal reps' L=150 retained draws are one correlated blob, biasing
their sigma ranks; at thin=1000 the moderate reps decorrelate but the
Cauchy-tail minority does not, which is exactly why point C passes while
sitting at 69.7 rather than 75. Bias direction (posterior sigma high) fits:
from the overdispersed init the (a, mu) block relaxes downward toward the
true amplitude over thousands of sweeps.

Routing: remedy filed as TODO bcf-ridge-interweaving (an ASIS/interweaving
joint rescale (a, mu) -> (a/c, c*mu) with c from its full conditional;
posterior-preserving, draw-changing, VD sign-off) - the acceptance test is
this same n=200 config at the cheap thin=120 setting passing sigma after the
move lands. No glue-path exactness dive is warranted. Artifacts (per-chunk
rank .rds, scratch.log, IACT scripts) in the session tmp directory.

#### A5. weighted gaussian (n=150, R=200) -- ALL PASS

Known positive weights w ~ Gamma(2,2) (mean 1); heteroscedastic noise
y0 = f0 + (sigma0 / sqrt(w)) * N. Reused sampler (setResponse works, scale
fixed -> no rebuild mismatch). All PASS: sigma (chisqP 0.806, ecdf 0.033),
avg f, all 5 f* (ecdf 0.043-0.068 << band 0.092). The weighted SSR / per-row
precision path calibrates. 56s.

### Tier B - linear leaf

R=200, L=150, thin=60, burn=50, n=150 p=3 nTrees=50, leaf covariates =
columns 1:2, k=2, unless noted.

Tier-B thinning MEASURED, not assumed: across three linear-leaf datasets the
slowest functional (sigma) first drops below ACF 0.1 at lags 44 / 60 / >60
(the linear leaf mixes sigma more slowly than the constant leaf's ~30), so
thin=60. New prior-draw component checks before trusting ranks: the recorded
training fits and predict(x) agree to machine precision at the same state --
baseline 3.6e-15, NA-in-leaf-covariate 1.8e-15, NA-in-split-column 2.7e-15 --
so theta0's f0 is exactly the f the likelihood uses, INCLUDING the NA routing
(rules carry the missing direction) and the imputation-at-standardized-mean
leaf path (model.hpp:173).

HARNESS NOTE (real API behaviour worth knowing): the matrix (xy) interface
DROPS NA rows even under missing = "incorporate" (dbartsData warns "row(s)
dropped"); only the formula interface (na.action = na.pass) keeps them. The
harness builds NA designs through a formula.

#### B1. linear leaf baseline -- ALL PASS

    functional  chisqP   ksP   ecdfDiff  band(.0924)
    avg.f       0.735   0.992   0.0280    PASS
    f.star1-5   --      --      0.032-0.080 all PASS (f.star4 chisqP 0.002 is
                                1 of 7 with ecdf 0.080 inside band = noise)
    sigma       0.748   0.207   0.0736    PASS (rank mean 72.1/75)

The linear-leaf ridge draw (vector params, U'WU statistics) calibrates.
~2.4 s/rep; 4 x 50-rep chunks.

#### B2. linear leaf + NA in a DESIGNATED leaf-covariate column -- PASS
(the coverage matrix's number-one target: model.hpp:173)

15% of rows carry NA in column 1 (a leaf covariate); the fixed NA pattern
rides every replication. Functionals add f at three NA-bearing training rows
(watch rows 9/21/23), ranked through the recorded training fits -- the ranks
sit exactly on the imputation-at-standardized-mean path.

    functional  chisqP   ksP   ecdfDiff  band(.0924)
    f.row9      0.407   0.419   0.0594    PASS  <- NA-imputation path
    f.row21     0.337   0.405   0.0621    PASS  <- NA-imputation path
    f.row23     0.735   0.718   0.0485    PASS  <- NA-imputation path
    sigma       0.294   0.383   0.0600    PASS  (rank mean 78.9/75)
    avg.f       0.684   0.832   0.0416    PASS
    f.star2-5   --      --      0.051-0.073 all PASS
    f.star1     0.337   0.028   0.1030    marginal (1 of 10 just over the 95%
                band, chisq healthy, not an NA functional = band-edge noise)

VERDICT: the missing-leaf-covariate imputation path that no gate executes
CALIBRATES. ~2.4 s/rep.

#### B3. linear leaf + NA in the SPLIT-ONLY column (control) -- ALL PASS

Same 15% NA pattern but in column 3 (NOT a leaf covariate): NAs are routed
purely by each rule's missing direction, leaf covariates stay complete.

    all 10 functionals PASS: sigma 0.0491 (rank mean 77.7/75), avg.f 0.0433,
    f.star1-5 0.044-0.074, watch rows f.row9/21/23 0.036/0.073/0.070
    (f.row23 chisqP 0.027 with ecdf well inside the band = noise).

VERDICT: NA split-routing calibrates; together with B2 the two distinct NA
paths (routing vs leaf imputation) are both clean.

#### B4. linear leaf + non-unit weights (w ~ Gamma(2,2)) -- ALL PASS (R=200)

    functional  chisqP   ksP   ecdfDiff  band(.0924)
    sigma       0.839   0.831   0.0425    PASS  (rank mean 75.1/75, centred)
    avg.f       0.470   0.916   0.0381    PASS
    f.star1-5   --      --      0.043-0.090 all PASS

The weighted linear-leaf path (weights entering U'WU / U'Wz) calibrates.
Tier B complete: all four linear-leaf configs clean.

### Tier C - GP leaf

R=200, L=150, thin=20, burn=50, n=80 p=3 nTrees=25,
gp(1L, k = 2, max.leaf.size = 100).

Setup notes (all measured, none assumed):
- max.leaf.size = 100 matches the equivalence gp scenario's cap; with n = 80
  no leaf ever truncates, so every leaf runs the TRUE GP path (no constant
  fallback). k is FIXED at 2: the equivalence scenario's chi hyperprior
  samples k, but sampleNodeParametersFromPrior draws at the CURRENT k and no
  API installs a hyperprior draw -- a sampled-k SBC would be prior-mismatched.
  RESIDUAL GAP: the sampled-k (chi hyperprior) channel has no SBC coverage
  for any leaf type (same class of limitation as DART alpha-sampling, held
  fixed there too).
- Sizing by measured cost: prior-drawn trees are shallow, so GP kernel solves
  run at leaf sizes ~n and cost scales ~n^3 (31 ms/sweep at n=150 vs 1.8 at
  n=80, nTrees=25). ~7 s/rep.
- Thinning measured fresh at the actual config: first lag with ACF < 0.1 is
  14-15 for sigma/f* (GP mixes FASTER than the linear leaf's 44-60; smooth
  leaves make tree moves less sticky). thin = 20.
- SELF-CONSISTENCY SUBTLETY (would have silently biased SBC): for
  function-valued GP leaves, predict() at the training rows re-krigs the
  stored leaf values (with jitter) and differs from the recorded training
  fits by ~2e-3, while the TEST path is shared exactly (0 diff). theta0's
  f0Train therefore reads the engine's stored internal fits via
  getForestFits, mapped to the reported scale by an exactly-recovered affine
  map (max residual < 1e-8, enforced); f* keeps predict(). With this, the
  recorded-fits consistency check passes at 1e-15 for all three GP configs.

#### C1. GP leaf baseline -- ALL PASS

    functional  chisqP   ksP   ecdfDiff  band(.0924)
    sigma       0.337   0.602   0.0534    PASS  (rank mean 75.2/75, centred)
    avg.f       0.040   0.486   0.0544    PASS
    f.star1-5   --      --      0.048-0.079 all PASS (chisqP 0.024-0.671;
                the lowish chisq values pair with clean ecdf/KS = binning
                noise)

The GP leaf's marginal/Matheron draw + kernel cache calibrates.

#### C2. GP leaf + NA in the designated leaf covariate -- ALL PASS

15% of rows carry NA in column 1 (the GP leaf covariate); watch-row
functionals ride the recorded training fits at three NA-bearing rows, exactly
on the GP missing-covariate imputation path.

    functional  chisqP   ksP   ecdfDiff  band(.0924)
    f.row9      0.942   0.959   0.0358    PASS  <- NA-imputation path
    f.row12     0.470   0.393   0.0614    PASS  <- NA-imputation path
    f.row16     0.119   0.244   0.0684    PASS  <- NA-imputation path
    sigma       0.617   0.732   0.0437    PASS  (rank mean 76.5/75)
    avg.f       0.563   0.339   0.0622    PASS
    f.star1-5   --      --      0.041-0.083 all PASS

VERDICT: the GP-leaf missing-covariate path (imputation at the standardized
mean inside a function-valued leaf) CALIBRATES -- with B2 this closes the
review-3 number-one uncovered combination for both leaf families.

#### C3. GP leaf + non-unit weights (w ~ Gamma(2,2)) -- ALL PASS

    functional  chisqP   ksP   ecdfDiff  band(.0924)
    sigma       0.630   0.652   0.0494    PASS  (rank mean 75.3/75, centred)
    avg.f       0.045   0.161   0.0753    PASS
    f.star1-5   --      --      0.039-0.080 all PASS

The weighted GP path (per-row precisions entering the leaf kernel solves)
calibrates. Tier C complete: all three GP configs clean.

## Final summary

R=200 per config; the ecdf-diff simultaneous 95% band is the verdict,
chisq/KS secondary.

    tier config                        verdict
    A    gaussian baseline             ALL PASS (7/7)
    A    probit baseline               ALL PASS
    A    logistic baseline             ALL PASS (6/6)
    A    DART (alpha=1, p=10)          s1..s10 ALL PASS (Dirichlet posterior
                                       covers the prior); sigma/avg.f/f* PASS
    A    DART sparse + 1e-300 floor    selection calibrates at 3% floor
                                       incidence; the floor does NOT bias
    A    grouped gaussian + 0-weights  tau/b1/b2/avg.f/f* PASS (sigma flag =
                                       documented harness rescale artifact;
                                       benign, grouped rebuild re-anchors the
                                       response scale)
    A    grouped probit                ALL PASS (tau at band edge, noise)
    A    weighted gaussian             ALL PASS (7/7)
    A    BCF standard (n=200)          abs.a / abs.diff / eff PASS; raw a and
                                       (b1-b0) sign-ill-posed BY DESIGN; sigma
                                       flagged initially (A4/A4b), RESOLVED as
                                       slow mixing on the (a,mu) scale ridge,
                                       not a defect (H-MIX, A4e)
    A    BCF fixed-glue control        ALL PASS (localizes the sigma pattern
                                       to the glue-on path; backfit exact)
    A    BCF prior-weak (n=40)         HEADLINE PASS - the a-glue prior
                                       precision (abs.a) calibrates; sigma
                                       centred (defect shrank with n,
                                       consistent with scale-ridge mixing)
    B    linear leaf baseline          ALL PASS (7/7)
    B    linear + NA leaf covariate    PASS - imputation path (model.hpp:173)
                                       calibrates (watch rows clean)
    B    linear + NA split column      ALL PASS (NA routing control)
    B    linear + weights              ALL PASS (7/7)
    C    GP leaf baseline              ALL PASS (7/7)
    C    GP + NA leaf covariate        ALL PASS - GP imputation path
                                       calibrates (watch rows clean)
    C    GP + weights                  ALL PASS (7/7)

18 configs; every one calibrates. The one config that initially looked like a
defect (BCF glue-on sigma, A4/A4b) was fully resolved by the chain-length
diagnostic (A4e): slow mixing on the (a, mu-amplitude) scale ridge, not a
sampler defect - exonerated within SBC resolution, with a mixing-efficiency
remedy (interweaving) filed as TODO bcf-ridge-interweaving rather than a
correctness fix. Residual gaps recorded: the sampled-k chi-hyperprior channel
(no API to install a prior k draw), DART alpha-sampling (held fixed), and the
half-Cauchy grouped tau (SBC-intractable via the engine's slice sampler;
gamma prior used instead). Harness findings: matrix-interface NA-row dropping
under missing = "incorporate"; GP predict vs stored training fits ~2e-3
(jittered re-kriging - harness reads stored fits). Candidate standing gate:
the gaussian constant-leaf baseline (~52 s at R=200 unloaded).
