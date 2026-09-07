# bcf-latent-evidence

Status: PROPOSED, 2026-09-07
agent: opus for the oracle and the arms; the derivation check is a second, independent pass
rng: neutral - `benchmarks/` and `.github/` only, so every baseline replays
budget: one harness (~600 lines), ~160 lines in `sbc.R`, ~20 lines in `sbc.yaml` plus one word in
  `exact-gates.yaml`, three equivalence scenarios and one re-record

## Goal

The K-forest amplitude family ships under gaussian, probit and logistic
([`refusedAmplitudeFamilyReason`](../../src/R_interface_bartcore.cpp)); every calibration and exactness gate it has is
gaussian. The [Gaps](../design/feature-matrix.md#gaps) row for the latent sub-families names three, equivalence, SBC and
the active-rows mask. After this item a probit and a logistic BCF each carry a deterministic exact-posterior gate, a
recorded SBC verdict and a bitwise equivalence scenario, and the row narrows to the mask alone.

## Context

Design: [Exact-posterior gate](../design/bcf.md#exact-posterior-gate),
[Calibration (2026-07-07)](../design/bcf.md#calibration-2026-07-07),
[The calibration map, general in K](../design/multiplier-combiner.md#the-calibration-map-general-in-k). Instruments:
[`exactBCF`](../../benchmarks/R/bcf-exact.R), [`runSbcBCF`](../../benchmarks/R/sbc.R),
[`exactPredictive`](../../benchmarks/R/logistic-reference.R) for the binary leaf quadrature idiom,
[`batchMeanSE`](../../benchmarks/R/linear-exact.R). Adjudication:
[The chain-length ladders (the A4e protocol)](sbc-family-tiers.md#the-chain-length-ladders-the-a4e-protocol),
[The flags and their adjudication](review-2026-08-24/calibration-sbc.md#the-flags-and-their-adjudication).

Four things are new under a latent family and they are what has to be gated. The latent refresh runs against the
COMBINED location `a mu + b_z tau`, not one forest's fits ([`refreshLatents`](../../src/bartcore/model.hpp),
[`Chain`](../../src/bartcore/chain.hpp)); the amplitude draw consumes the WORKING response and weights, so under
logistic its precision is a Polya-Gamma variate redrawn every sweep
([`drawAmplitudes`](../../src/bartcore/combiner.hpp)). Sigma is pinned at exactly 1
([`sigmaIsFixed_`](../../src/bartcore/chain.hpp), [`ProbitResponse`](../../src/bartcore/model.hpp),
[`LogisticResponse`](../../src/bartcore/model.hpp)). And the map's anchor is the link's latent scale, 1 under probit and
`pi/sqrt(3)` under logistic ([`latentScaleAnchor`](../../src/bartcore/chain.hpp)), at a default half-Cauchy scale of 1,
not gaussian's 2 ([`defaultAmplitudePriorScale`](../../R/model.R)). The glue IS drawn under both links, the sweep's
[`drawGlue`](../../src/bartcore/combiner.hpp) call being unconditional.

## Decision 1 - the SBC arms

The gaussian arm today. [`sbcBCFGlueDraw`](../../benchmarks/R/sbc.R) draws `a ~ Cauchy(0, sd.control)` and `b0, b1 ~
N(0, bPriorVariance)`; [`sbcInstallBCFGlue`](../../benchmarks/R/sbc.R) installs them through the state BEFORE the
forests, the prior tree draw conditioning on each forest's own veto vector
([`formForestVetoWeights`](../../src/bartcore/combiner.hpp)); the engine's prior-draw entry points supply the forests;
sigma comes from the reported-scale scaled-inverse-chi-squared; `y0` is the affine-mapped index plus gaussian noise. The
fit re-inits from a second prior draw, `setResponse` swaps `y0` at the pinned scale, and draws come one at a time, the
glue and per-forest fits being current-state only. Fifteen functionals: `sigma`, raw `a` and `b1.minus.b0`, `abs.a` and
`abs.diff`, and `prog_j = a mu(x_j)`, `eff_j = (b1-b0) tau(x_j)` at five rows; the raw glue pair is ill-posed (the index
is sign-flip invariant) and is reported to show it. The latent arms change five things.

1. `sbcAddBCF(sbcConfig(family = "probit"|"logistic", n = 200L, nTest = 3L), sdControl = 1)` - `n` the gaussian arm's,
   `sdControl` being `sbcAddBCF`'s argument and not `sbcConfig`'s - with [`sbcMakeBCF`](../../benchmarks/R/sbc.R)
   passing `family` to [`bartcoreBCFSampler`](../../R/bartcore.R), whose `family` formal writes it into the model the
   bridge reads. The transform is then the identity, so the regressed map becomes a self-check: scale 1, shift 0, R2 1.
2. `y0 ~ Bernoulli(pnorm(index0))` / `Bernoulli(plogis(index0))`, no offset, no sigma drawn. That functional goes and
   the sigma moment check gives way to [`sbcCheckBCFGlue`](../../benchmarks/R/sbc.R) plus a latent consistency check:
   recorded combined train fits equal `a mu + b_z tau` to 1e-12.
3. Add `p_j`, the link at each evaluation row's index - the reported deliverable, bounded so its histogram reads
   cleanly, where neither `prog_j` nor `eff_j` alone is the index and the increasing link makes its RANKS the index's.
   At `nTest = 3` an arm scores thirteen: 4 glue, 3 prog, 3 eff, 3 p.
4. Chain length. The 72000-sweep burn is NOT the sigma channel's and pinning sigma does not discharge it: the recorded
   mechanism is the `(a, mu)` amplitude ridge co-relaxing with tree-structure mixing at a settle time scaling in
   `|a|/sigma`, with sigma the READOUT
   ([Burn-in under strong prognostic signal (2026-07-10)](../design/bcf.md#burn-in-under-strong-prognostic-signal-2026-07-10)).
   Pinning sigma leaves the ridge untouched and removes the functional that reported it, so the misfit sigma absorbed
   lands in the index instead - in `prog_j` and `p_j`, which these arms rank. The family default `sd.control = 1` halves
   the Cauchy tail rather than removing it: `P(|a| > 5)` is 0.126 and `P(|a| >= 40)` is 0.016, the stratum bcf.md
   records as a mixing limit no burn fixes. So no burn is pre-registered, and in particular not thin 120 over
   18000-36000 sweeps, the setting at which the gaussian arm flagged. Reprice the ladder: give
   [`sbcFamilyConfig`](../../benchmarks/R/sbc.R) and [`sbcFamilySpec`](../../benchmarks/R/sbc.R) BCF branches so
   [`sbcBurnLadder`](../../benchmarks/R/sbc.R) takes `burn-bcf-probit` and `burn-bcf-logistic`, and run 24 prior-drawn
   datasets - three expected over `|a| > 5`, at least one with probability 0.96, against 0.33 at three datasets - plus
   four at PRESCRIBED `|a| in {0.5, 2, 5, 10}`. Thin and burn come out of that, into a table keyed by ARM name:
   [`sbcBurnSweeps`](../../benchmarks/R/sbc.R)'s `config$family` key would collide with a plain probit arm. `R = 200`,
   `L = 150`.
5. One arm per link, not one arm with counts: the links are separate sampling code with different working-weight
   behaviour, and probit refuses weights outright. The trial-count channel - logistic weights are frequency weights,
   `omega ~ PG(w, psi)` ([The two interpretations](../design/weighted-logistic.md#the-two-interpretations)) - is better
   gated deterministically by the exact gate's aggregated arm than by a third expensive statistical arm.

Admission and CI. An arm is admitted only when its ladder shows, at the chosen thin and burn, every ranked functional's
ACF under 0.1 within one thinning interval and the block-mean z transient settled ON the worst `|a|` stratum the prior
gives non-negligible mass, and when the `R = 200` run then passes with its tail replications kept. If the `|a| >= 40`
stratum cannot be cleared at an affordable burn - the recorded expectation - the arm is reported as a finding naming
that stratum, not admitted. The matrix exclusion names the glue-on SIGMA channel
(["Not in the matrix: BCF"](../../.github/workflows/sbc.yaml)); pinning sigma removes that readout, not the ridge behind
it, so admission is earned from the ladder, not inherited. On admission: both arms with `SBC_EXPECTED_FLAGS:
a,b1.minus.b0`, the ill-posed pair waived by name as nbinom's ridge pair is
(["SBC_EXPECTED_FLAGS: r,agg.psi"](../../.github/workflows/sbc.yaml)), and
[`sbcMatrixFunctionals`](../../benchmarks/R/sbc.R) 30 to 56 beside [`sbcMatrixConfigs`](../../benchmarks/R/sbc.R), so
the band WIDENS for every existing arm - no rank moves, no recorded PASS at risk. The gaussian arm stays out.

Adjudicating a flag, in order. The ill-posed pair is waived; then the A4e ladder, three points over at least 8x chain
length at fixed `R`, monotone shrinkage in `ecdfDiff / band` being H-MIX and a plateau a defect candidate; then the
arm's controls, `fixedGlue = TRUE` isolating the backfit from the glue draw and `bcf-weak` at `n = 40` making the glue
prior dominate, where an exactness error shows MORE and a ridge LESS. On a plateau the exact gate can exonerate only the
single-tree, two-cell conditional it runs; a plateau surviving that points at ensemble scale or the continuous cut grid,
and the next instrument is an independent derivation of the flagged conditional.

## Decision 2 - the exact gate

Alternatives: (A) keep the tree enumeration and integrate the leaf parameters by adaptive Gauss-Hermite quadrature; (B)
a Monte-Carlo importance oracle over the same space; (C) a reduction to a covered family.

(C) does not exist. A probit BCF reduces to no shipped single-forest model: the forests carry different tree priors and
leaf scales, the treatment forest enters through a basis rather than a split, and a model is refused below two
forests. (B) is dominated as a gate, and adopted as a check. At block dimension 3 the prior-importance ESS fraction
is about 2e-3, so 1e7 draws give a posterior-MEAN error near 8e-4 - eighty times the bar below, and the order of the
sampler's own MC error. A Laplace-centred proposal fixes that and then carries (A)'s machinery plus noise; it is the
right cross-check on two configurations (Decision 3).

Recommendation: (A), on one ordinal predictor with `K = 2` cells, `n.cuts = 1`, 200 rows per cell and `z` balanced
within each - four `(cell, z)` groups of 100, `n = 400` - no offset, and leaf values `mu = (0.0, 0.3)`, `tau = (0.5,
0.8)`, so the index `mu_c + z tau_c` at the four groups is `(0.00, 0.50, 0.30, 1.10)` and `p` is `(0.500, 0.691, 0.618,
0.864)` under probit and `(0.500, 0.623, 0.574, 0.750)` under logistic, well inside the unit interval, which is also what
keeps `0 < s_g < n_g` at every group for the aggregated arm below. The contrast is 0.3 in BOTH forests so that more than
one configuration stays alive: at 100 rows per group the split-versus-stump Bayes factor is `exp(O(n))`, and at the
design's expected counts the four configurations' weights are `(0.002, 0.007, 0.899, 0.092)` under probit and `(0.045,
0.011, 0.886, 0.058)` under logistic. Over 200 simulated realizations of the design the largest weight runs 0.44-0.95
with a median near 0.85 while the runner-up carries 0.05-0.20 and never vanishes, so the guard is stated from below: the
script recomputes the weights at its realized counts and refuses to run if the SECOND-largest falls under 0.02. (A cap of
0.9 on the largest would refuse a third of realizations, and the expected counts themselves under probit.) At two cells each
forest's tree has AT MOST TWO LEAVES by construction: after the one split each child holds a single cell with an empty
cut interval, so the CGM growth probability there is exactly zero and [`enumerate`](../../benchmarks/R/bcf-exact.R)
yields two trees per forest, four joint configurations. Balanced `z` leaves every realizable leaf positively weighted,
so the veto never restricts the prior.

Dimension count. Conditional on a tree pair the leaf parameters couple only through the cells, cell `c` linking mu-leaf
`l(c)` to tau-leaf `m(c)`, so the integral factorizes over that bipartite graph's connected components: blocks of
dimension 2; 3; 3; and 2 + 2, the maximum THREE. Per block, Newton to the mode - the log integrand is strictly concave -
then a `d`-fold adaptive Gauss-Hermite product rule at the mode and curvature, at 12 nodes per dimension, where the
dim-3 block moves 2.6e-13 against a 200-node reference. That is 3,888 integrand evaluations per glue point over the four
configurations. (Three cells would raise the maximum block to 4 and that count 38-fold, to 148,752, unaffording the
free-glue modes.)

Modes and matched quantities, on the index scale, read through
[`bartcoreForestFits`](../../inst/common/bartcoreHandle.R),
[`bartcoreForestAmplitudes`](../../inst/common/bartcoreHandle.R):

    1   glue fixed (1, 0, 1)                E[mu_c], E[tau_c]
    2a  a ~ Cauchy(0, sd.control = 1)       E[a mu_c], E[tau_c]
    2b  b0, b1 ~ N(0, 0.5), a = 1           E[mu_c], E[(b1-b0) tau_c]

plus, in every mode, `E[F(eta_cz)]` at the four groups - the reported probability surface, ridge-invariant, and what the
binary gates already match. The sigma quadrature disappears, the initial sigma being exactly 1 with no draw to move it.

Glue axes. Both take the same substitution, and the truncated tensor grid the gaussian gate uses does not meet the bar.
`a = sd.control tan(t)` on an open trapezoid is spectrally accurate - 101 points give 8e-12 and 201 machine precision
against a 2001-point reference - and `b0, b1` take `sqrt(bPriorVariance) tan(t)` on the same rule: at 81 points per axis
`E[mu]` and `E[(b1-b0) tau]` sit at 2.9e-8 and 5.7e-9, where the truncated 45x45 on [-3.5, 3.5] is at 1.1e-3 and even
101x101 on [-4.5, 4.5] is at 3.8e-5 - above the bar, and a systematic offset the gate could not see.

Quadrature accuracy against MC error. Gate on batch-means standard errors at `zBound = 4`, the
[`batchMeanSE`](../../benchmarks/R/linear-exact.R) idiom, rather than the gaussian gate's absolute tolerance: the binary
leaf posteriors are 0.13 to 0.31 wide by channel (probit mu the narrowest, logistic tau the widest) against that gate's
0.0034, so an absolute bound carries no fixed meaning across the two - keeping 0.015 would leave a bound of about 12 se,
not 4. The standard error is per channel and per seed, the seeds pooled as the mean of per-seed means with
`sqrt(sum se_s^2) / S`, never batched across a chain seam. At 100000 kept draws thinned by 10 (25000 by 5 in quick mode)
over three seeds the standard error lands near 6e-4 on the narrowest channel and about 2.5x that on logistic tau, the
detection threshold near 2.4e-3 and 6e-3, so the quadrature must be good to about 1e-5. Mode 2b is the exception: `E[mu]`
couples to the slowly mixing `(b0 + b1) / (b1 - b0)` ratio, which the gaussian gate met with thin 200 and eight seeds, so
this gate does the same there and reports the batch-mean lag-1 autocorrelation beside every z; a 400-batch se over a
thin-10 chain would understate the error and fail a correct sampler. Nine configurations by up to eight matched
quantities is about 70 tests at `|z| <= 4`, a family-wise false-failure rate near 4e-3, accepted here in advance. The
refinement self-check therefore runs on the GLUE axis, where the error is: recompute the
dim-3 configuration at 121 nodes per `b` axis and 401 `t`-points in `a`, and refuse to run if any reported quantity
moves more than 1e-6. Refining the LEAF rule would certify nothing - at 1e-13 by 12 nodes, a doubling reports zero.

A fourth arm at no oracle cost: the logistic mode-1 arm run again on AGGREGATED data, two rows per group, `y = 1` with
count `s` and `y = 0` with count `n - s`. The oracle reads the group sufficient statistics alone, so it is the same
target bit for bit and both fits must meet it - the trial-count path, the one channel probit cannot reach.

A fifth arm for the tree prior's power. At two cells the only interior node is the root, where `base / (1 + depth)^power`
equals `base` for every power, so a poisoned power passes every mode and both links silently, and no other gate covers
it under a latent family (the gaussian gates cover the cut factor at `K = 3` and the variable factor at two predictors).
A `K = 3`, mode-1-only arm at each link - five trees per forest, 25 configurations, maximum block dimension 4 - costs
148,752 integrand evaluations at its single glue point, under one percent of the free-glue grids, and restores
depth-decay and cut-selection coverage. Left uncovered, by design: the basis forest's ridge (`ridgeB` ships off, and mode
2b holds the prognostic block) and the multiplier snap's near-tolerance boundary (mode 1 fires it on every control row,
modes 2a and 2b never).
Cost, estimated here and measured in the slice: the quadrature is about a minute per link at full settings, the sampler
side dominates, the per-sweep cost at `n = 400` with one tree per forest being the latent refresh, one truncated normal
or Polya-Gamma draw per row. Nine sampler configurations (three modes x two links, the aggregated arm and the two `K = 3`
arms) at three seeds put the full grid near 12 minutes and quick near 3, so quick is what CI runs inside exact-gates.yaml's shared
30-minute job while the full grid is workflow_dispatch. Harness about 600 lines in a new
`benchmarks/R/bcf-latent-exact.R`, against [`exactBCF`](../../benchmarks/R/bcf-exact.R)'s 446 for gaussian, closed-form and one link.

## Decision 3 - the independent derivation

Before either gate lands, a second pass derives four objects from the design notes and the engine, not the harness: the
per-configuration marginal likelihood (the block decomposition and its justification, and the integrand summing `s_g log
F(eta_g) + (n_g - s_g) log(1 - F(eta_g))` over groups at an index linear in the leaf parameters); the tree prior mass
(the CGM product with the forced-leaf term on an empty cut interval and the uniform cut-selection factor); the glue
prior under a latent family (`a` the half-Cauchy scale mixture [`drawAmplitudes`](../../src/bartcore/combiner.hpp)
samples, whose marginal must be the Cauchy the oracle integrates, at `sd.control` defaulting to 1 and NOT rescaled by
sigma, the b pair normal at `bPriorVariance`); and the leaf scales (mu normal at `s`, tau at `sd.moderate s / 0.674`,
for `s` of 1 and `pi/sqrt(3)`, `k` pinned at 1). Checked against three things, not one: the engine source; the shipped
calibration reader, which decomposes `prior.scale` as `node.scale.factor * anchor / (node.scale.divisor *
basis.row.norm)` and reports `amplitude.prior.scale`, read off a public-route latent BCF built with the matching
declaration, where a disagreement between the creation routes is itself the finding; and a prior Monte-Carlo estimate of
one block marginal and posterior mean at two configurations, agreeing to its own MC error - the decisive one, the only
check not sharing the quadrature's assumptions.

Outcome (2026-09-07). Derived blind from the engine and then compared: the enumeration, the block dimensions 2/3/3/2+2,
the integrand, the glue laws (`a` an inverse-gamma scale mixture whose marginal is Cauchy(0, 1), each `b` normal at 1/2,
sigma nowhere), the leaf scales and the pinned sigma all match the engine. `getCalibration` agrees across the
evaluated-basis route, the `~ factor(z)` route and [`bartcoreBCFSampler`](../../R/bartcore.R). A prior Monte Carlo of the
block marginals at 2e7 draws agrees with two independent quadratures to its own error, and the shipped sampler matches
the derived oracle at fixed glue to 7e-4 (probit) and 1.6e-3 (logistic). The pass corrected this plan's design line (the
treated cell-2 index is 1.10, not 0.80), its configuration weights and its weight guard, and added the `K = 3` arm; the
anchor readings and poison (ii)'s magnitudes below are as remeasured.

## Decision 4 - discrimination

Each poison is run once, by hand, and recorded here; none lands, and none edits the engine. Exact gate: (i) score the
probit fits against the LOGISTIC oracle and the logistic fits against the probit one; both must fail, or the gate is
measuring the two-forest structure and not the link. (ii) Run mode 2a with `update.a = FALSE` against the free-`a`
oracle: on one fitted design the gap in `E[a mu]` is 2.0e-2 under probit and 5.6e-2 under logistic, and `E[tau]` moves
6.9e-3 and 2.9e-2, every one above the detection threshold, so the poison is named on both quantities; the magnitudes are
realization-dependent. The free-`a` posterior has `E[a] = 0` exactly by the `(a, mu) -> (-a, -mu)` symmetry, which is why
`E[a mu]` and not `E[a]` is the mode-2a target. (iii) Build the oracle's leaf scales at the sample sd of the RANGE-SCALED
binary response, 0.482, in place of the link's latent anchor: the factors are 2.08 under probit and 3.76 under logistic
and every matched quantity must fail (gaps of 2.8e-2 and 8.8e-2 measured). The anchor has to be named, because the other
reading of "gaussian's anchor" - the engine's own, the sd of the cold-start working response - is 0.96 under probit and
1.98 under logistic, prior errors of -4% and +9% whose gaps (7e-4 and 1.9e-3) sit under the threshold and would not
land, the logistic one at 0.8x the threshold and realization-sensitive. (iv) Score the AGGREGATED logistic arm against
an oracle built at unit counts, two rows per group in place of the
true trial counts: it must fail, or that arm gates the row layout, not the counts.

SBC arms. (i) Simulate `y0` through `plogis` while fitting probit: `p_j` and `prog_j` must FLAG. (ii) Draw theta0's glue
at gaussian's `sd.control = 2` while the sampler runs at the family default 1, so every replication's `a` comes from the
wrong prior and `abs.a` and the `p_j` cells must FLAG. This replaces drawing the forests before the glue, which cannot
flag: that swap moves the tree law in about one replication of 200. (iii) Draw a sigma in the generator as noise the fit cannot model: every latent functional must FLAG. Only (ii)
is glue-targeted, so only (ii) carries the control that `fixedGlue = TRUE` stays clean; (i) and (iii) redden it too.

## Constraints

Out of scope, each named rather than absorbed:

- Active-rows-mask evidence for the latent sub-families. Different shape, not to be smuggled into either instrument: an
  inactive row consumes NO latent variate under either link, a STREAM property, so the evidence is a bitwise pair - a
  latent BCF fit on the retained rows against one on the full data with the complement masked - not a posterior gate.
  Roughly 80 lines; it is what remains of the Gaps row, and pairs with
  [Per family](../design/active-rows-mask.md#per-family).
- The engine. Nothing here needs an engine, `R/`, `inst/` or `man/` change and none is proposed; every channel is
  reachable from the installed package.

The three latent `bcf-equivalence` scenarios ARE in scope and cheap: the bcf compare loops over the BASELINE's scenario
names, so scenarios appended with literal seeds kept out of the guarded settings list (the precedent the
forced-`setPredictor` scenario set) perturb no stream and invalidate no pinned baseline; they go unchecked until a
re-record, which repins one hash in one file. Order: the exact gate FIRST - deterministic, cheaper to build and run, it
gates every push, and it forces the derivation the SBC generator then reuses.

## Steps

1. `benchmarks/R/bcf-latent-exact.R`: the design and its runner-up weight check, enumeration, block decomposition,
   adaptive GH at 12 nodes, both glue axes on the tan substitution, the glue-axis refinement check, three modes x two
   links plus the aggregated logistic arm and the two `K = 3` mode-1 arms, per-seed batch-means z at `zBound = 4` with
   mode 2b at thin 200 and eight seeds; and its word in `.github/workflows/exact-gates.yaml`.
2. The independent derivation (Decision 3, done) and the four exact-gate poisons (Decision 4), recorded here; then three
   `bcf-equivalence` scenarios and one re-record.
3. `benchmarks/R/sbc.R`: the two latent configs, the latent branch of the BCF maker and runner, the `p_j` functional,
   and the `sbcFamilyConfig` / `sbcFamilySpec` branches the burn ladder needs.
4. Run the repriced ladders (24 prior-drawn datasets plus the four prescribed `|a|` strata), record the sweeps and
   per-sweep costs, fix thin and burn; then the `R = 200` verdicts plus the fixed-glue and weak controls, an A4e ladder
   point for anything flagged, and the three SBC poisons.
5. Only if the admission criterion is met: the `.github/workflows/sbc.yaml` entries, `sbcMatrixConfigs` and
   `sbcMatrixFunctionals`, the gaussian arm's ranks confirmed to replay with only the band moving. Then the records: the
   feature matrix Gaps row, the design note's exact-gate section, the review tour's section 5 bullet.

## Verification

    Rscript benchmarks/R/bcf-latent-exact.R quick     # and without 'quick'
    Rscript benchmarks/R/bcf-equivalence.R compare benchmarks/baselines/<new>.rds
    Rscript benchmarks/R/sbc.R burn-bcf-probit 40000 24     # and burn-bcf-logistic
    Rscript benchmarks/R/sbc.R bcf-probit   200 150 <thin>  # and bcf-logistic
    Rscript benchmarks/R/sbc.R gaussian 100 200 30          # ranks byte-identical
    Rscript tools/check-doc-freshness.R .

Expected: every exact-gate quantity inside `|z| <= 4` at both links and in all four arms, the glue-axis refinement under
1e-6, the runner-up configuration at or above 0.02 posterior weight; every SBC functional PASS but the two ill-posed
glue ones.
