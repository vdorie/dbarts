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

## Status

- 2026-07-09: plan authored; pilot dispatched.
