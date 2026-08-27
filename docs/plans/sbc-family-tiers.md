# sbc-family-tiers

agent: opus implementer; fable adjudicates the rank histograms
status: BUILT. See "Results" for the measured ladders and the R=200 verdicts.
rng: neutral - benchmarks/ and .github/ only; gaussian CI ranks replay exactly
budget: one sbc.R block per family, a discrete-rank helper, a CI matrix

## Goal

A recorded SBC verdict for every shipped family, or the reason it is ill-posed.

## Context

- sbc-calibration.md; sbc-ci-gate.md, whose admission rule step 4 replaces; the
  design notes. runSbcBCF errors at HEAD (9030d93): OUT of scope, fixed separately.

## Decision - scope

The fork: which families, at which prior. IN below at R=200, L=150, thin=30, burn
MEASURED (step 2). Changes it: a k=8 cost over budget makes nbinom manual-only; a
burn ladder against the 72000 floor re-opens burn for every arm.

    R=200 at the SHIPPED 72000-sweep burn, 38/48/56 us/sweep: 10-15m, not ~2m
    ordinal   gamma_2..K-1, eta, agg p  marginal-MH cutpoints     ~14m
    nbinom    r, avg mu, agg psi [k=8]  grid r cond., PG y+r         ?
    t         sigma, nu, avg f, agg f*  lambda mixture, nu grid   ~12m
    multinom  agg p_k(x*), raw f_ik     interleaved PG, centering    ?

- nbinom at a TIGHTENED k = 8 (psi sd = node.scale/k = pi*sqrt(3)/8 ~ 0.68, vs 2):
  simulatePolyaGammaShape loops sum(y_i + r) times/sweep, lognormal-tailed and
  unbounded under default-k draws (13.5x over six, tail reps hours), so it is not
  budgetable; a tightened prior still validates NB. RE-MEASURE the sweep cost at 8.
- multinom also ranks raw f_ik at a few (i, k) cells via bartcoreForestFits (the
  BCF accessor): the softmax flat direction is no prior symmetry, so f_ik SBC is
  well-posed and tests MultinomialForestCombiner::afterCombine's level-centering
  draw, whose precision (sum of invV*n) approximates the exact per-leaf
  conditional - the pre-registered suspect if f_ik alone flags.
- r and nu are DISCRETE, state-only: rank = #{draws < theta0} + Uniform{0..#ties}
  - theta0 takes its own iid tie-breaker beside the L draws, so the null is
  rankUniformity's uniform on {0..L}; collect by per-sample run(0,1)+storeState.

OUT, with the change that would enable each:
- aft: status is fixed at creation, so a status-varying DGP forces a rebuild that
  re-anchors range(y0). Enabler: a status setter.
- hazard: the at-risk person-period expansion makes the design - and so the tree
  prior (sampleTreesFromPrior's cut grid, empty-node collapse) - depend on y0,
  breaking exchangeability; second, it owns no sampling code (hazard-reduction.R
  gates it bitwise vs a hand-expanded fit).
- hurdle: no new sampling code; same design-depends-on-y0 defect (y > 0 subset).
- heteroscedastic: prior draws never reach varianceForest_; liftable R-side via
  setState (state carries variance.vars/values/sizes/flags). Deferred, not blocked.
- monotone: the enabler (a constrained prior draw) LANDED 173a710 the day this
  plan was written (monotone-prior-draw.md), so monotone is now liftable; first
  follow-on arm after the initial matrix, not in it.

## Constraints

R entry points of the installed package only, no engine edits. Out of scope: BCF
repair, zero-inflation, grouped families, xbart, Tier A-C, sampled-k/DART/tau,
and the critique's monotone bonus finding (a SEPARATE item, not this plan).

## Steps

1. Discrete tie-break helper + a sbcCheckSigmaPrior-style self-check: a synthetic
   discrete conjugate case with closed-form posterior must pass rankUniformity.
2. Measure a per-family burn ladder first - 72000 was a BCF-specific floor - and
   for ordinal pre-commit to the A4e chain-length ladder before reading any
   cutpoint/eta flag as a defect (ordinal.md section 9's cutpoint-shift ridge).
3. Per family in table order: prior draw + moment check, likelihood, functionals,
   thin/burn, R=200 run, verdict. ordinal and nbinom REBUILD per rep (safe: scales
   are 1, node.scale constant; it clears the kept gamma_/r_ that would break rank
   iid-ness); only t reuses a pinned sampler (setResponse cold-inits nu, lambda).
4. sbc.yaml matrix (fail-fast: false, timeout from step 2). Admission: the alpha
   Bonferroni'd to 0.05/M over the matrix's functional count (full-matrix pass
   ~0.95 on a fresh stream); drop gamma_1 (pinned at 0), rank gamma_2..K-1 with
   K >= 4, aggregate the (k, x*) cells. Seed pinned; risk: a draw-shifting commit
   or CI image re-rolls the stream, so the gate stays non-blocking.

## Verification

- `Rscript benchmarks/R/sbc.R gaussian 100 200 30` -> identical ranks.
- `... ordinal|nbinom|t|multinom 200 150 30`, one each -> every functional PASS.
- `... discrete-selfcheck` -> rankUniformity PASS on the conjugate case.

## Results

Machine: the maintainer's laptop (arm64, single thread), private library, one R
process at a time. Admission is the ecdf-diff band Bonferroni'd to
0.05/30 = 0.00167 over the matrix's functional count (gaussian 7 + ordinal 10 +
nbinom 3 + t 4 + multinomial 6); at R=200, L=150 that band is 0.1282.

### Step 1 - the discrete rank and its self-check

`sbcDiscreteRank`: rank = #{draws < theta0} broken by an iid Uniform(0, 1) tag
carried by every draw AND by theta0, so an atom contributes a Uniform{0..#ties}
increment and the total is uniform on {0..L}. `sbc.R discrete-selfcheck` ranks
r0 against a CLOSED-FORM conjugate posterior (r0 from the engine's own nbinom
grid prior, y ~ NB(r0, p) at known p, posterior propto prior x product dnbinom
over the same grid), so the L draws are exact and iid and any non-uniformity is
the rank rule's fault:

    grid priors: nbinom max cell diff 0.0018, mean 10.181 vs 10.183 -> PASS
                 t      max cell diff 0.0009, mean  9.411 vs  9.427 -> PASS
    ranks (R=400, L=200, 71.4% of draws tied on average):
                 chisqP 0.881  ksP 0.861  ecdfDiff 0.0266  band 0.0659  PASS

Negative control on the SAME case ranked WITHOUT the tie-break: ecdfDiff 0.5607
vs band 0.0659, chisqP 0 - the helper is load-bearing, not decorative.

FINDING (harness, and the reason the driver tie-breaks EVERY functional rather
than only the two grid parameters): a numerical atom needs the tie-break exactly
as a discrete parameter does. The ordinal top-category probability
p_K = mean_i (1 - Phi(gamma_K-1 - eta_i)) UNDERFLOWS to exactly 0 whenever the
prior draws the top cutpoint far out - over a quarter of replications at K = 4,
the empty-cell case ordinal.md section 9 names. A first ordinal R=200 run with
the plain rank flagged p4 (ecdfDiff 0.1319 vs 0.1282) with a pure rank-0 pileup
(33 of 200 in the first bin, flat elsewhere); instrumenting 60 replications
showed EVERY rank-0 replication had p4_0 == 0 exactly, zero observations in
category 4, and up to 110 of 150 posterior draws tied at 0. Same class as the
DART 1e-300 floor probe's tie-degenerate ranks: a ranking-rule defect, not a
sampler one. With the tie-break p4 lands at 0.0926 and p3 at 0.0413.

### Step 2 - the measured burn ladders

`sbc.R burn-<family> 40000 3` runs 40000 UNTHINNED sweeps from an independent
prior init on each of three prior-drawn datasets and reports, per functional,
each 4000-sweep block's mean as a z against the settled half plus the first ACF
lag under 0.1. 72000 was a BCF number and none of these arms needs it; two of
them are set by a LIKELIHOOD RIDGE rather than a transient.

    family       us/sweep   slowest functionals (first acf < 0.1, 3 datasets)
    ordinal        96       gamma2/gamma3 NA(>200) / 10-179 / 110-114; p2 NA/14/96
                            every eta functional 4-24; p1/p3/p4 5-16
    nbinom        277       r NA(>200) / 80 / NA; agg.psi the same, MIRRORED
                            avg.mu 1 (the identified quantity)
    t              93       sigma 37-41, nu 13-60, avg.f 2-4, agg.f.star 4-23
    multinomial   245       everything 1-9

The two slow arms are ridges, and the ladder names them: ordinal's free
cutpoints trade against the mean level (ordinal.md section 9's
f-vs-cutpoint-shift ridge - which the design explicitly leaves mitigated, not
eliminated), and nbinom's r trades against the psi level because only
mu = r exp(psi) is identified (r and agg.psi z-blocks mirror each other sign for
sign while avg.mu decorrelates at LAG 1). Burn taken as ordinal 36000, nbinom
24000, t 12000, multinomial 6000 absolute sweeps (`sbcBurnSweeps`).

Measured per-sweep cost at k = 8 for nbinom (the plan's open number): 277
us/sweep in the unthinned ladder, which includes one R round trip and the
functional evaluation per sweep; in situ at thin = 30 the arm runs 17.0 s/rep,
i.e. ~57 min for R=200 - budgetable, so nbinom stays in the matrix rather than
becoming manual-only.

### Step 3 - the R=200 verdicts (L=150, thin=30, band 0.1282)

t - ALL PASS, 173 s (0.87 s/rep)

    sigma       chisqP 0.671  ecdfDiff 0.0655  PASS
    nu          chisqP 0.028  ecdfDiff 0.0771  PASS   <- the grid df draw
    avg.f       chisqP 0.161  ecdfDiff 0.0644  PASS
    agg.f.star  chisqP 0.255  ecdfDiff 0.0739  PASS

The lambda scale-mixture, the conjugate sigma under the composite precisions,
and the grid nu full conditional all calibrate.

ordinal - 9 of 10 PASS, 422 s (2.11 s/rep)

    gamma2      chisqP 0.042  ecdfDiff 0.0536  PASS
    gamma3      chisqP 0.006  ecdfDiff 0.1378  FLAG   <- see the ladder below
    avg.eta     chisqP 0.910  ecdfDiff 0.0397  PASS
    eta.star1-3               ecdfDiff 0.0435-0.0619  PASS
    p1-p4                     ecdfDiff 0.0413-0.0926  PASS

nbinom - 1 of 3 PASS, 3092 s (15.5 s/rep)

    r           chisqP 0.000  ecdfDiff 0.1725  FLAG
    avg.mu      chisqP 0.617  ecdfDiff 0.0362  PASS   <- the identified quantity
    agg.psi     chisqP 0.000  ecdfDiff 0.1782  FLAG

r's histogram is a U with a heavy right tail (19 in the first bin, 39 in the
last, 2-12 through the middle) and agg.psi mirrors it. A U means the L retained
draws are too tightly clustered for the spread the prior implies - the
correlated-blob signature - and that is exactly what the burn ladder predicted
for this arm: r's ACF is still above 0.1 past lag 200, so at thin = 30 the 150
retained draws carry an effective sample size of a few dozen, while the
identified mu = r exp(psi) mixes at lag 1 and calibrates perfectly.

multinomial - ALL PASS, 457 s (2.29 s/rep)

    p1  p2  p3          ecdfDiff 0.0428 / 0.0397 / 0.0524  PASS
    f.1.1 f.2.2 f.3.3   ecdfDiff 0.1107 / 0.1144 / 0.1170  PASS, but see below

The aggregated softmax probabilities are clean. The three RAW per-forest f_ik
cells are not comfortable: all three sit just inside the matrix band and would
FLAG at a per-functional 5% band (0.0924), all three carry chisqP 0.000, and all
three have the SAME shape - a U (f.1.1: 20 in the first bin, 29 in the last, 4-8
through the middle), i.e. posterior draws too tightly clustered for the spread
the prior implies. That is the plan's pre-registered suspect,
MultinomialForestCombiner::afterCombine's level-centering draw, whose precision
(sum of invV*n) approximates the exact per-leaf conditional; the mixing
alternative is weak here because the ladder puts every f cell's ACF under 0.1 by
lag 9 at thin = 1, so thin = 30 draws are near-independent.

### The chain-length ladders (the A4e protocol)

Every flagged or near-flagged arm was re-run at 3x the spacing with at least as
much burn, the A4e discipline: if the bias SHRINKS the flag is mixing, if it
PLATEAUS it implicates the sampler. Run at REDUCED R=100 to stay inside the
day's budget, so the comparable quantity is ecdfDiff / band (both scale like
1/sqrt(R)); the full-R commands are recorded below for a later pass.

    arm / functional     point A (thin 30, R=200)   point B (3x-5x spacing)
    ordinal gamma3       0.1378 / 0.1282 = 1.08     0.0574 / 0.1872 = 0.31
    multinomial f.1.1    0.1107 / 0.1282 = 0.86     0.2013 / 0.1872 = 1.08
    multinomial f.2.2    0.1144 / 0.1282 = 0.89     0.0990 / 0.1872 = 0.53
    multinomial f.3.3    0.1170 / 0.1282 = 0.91     0.1383 / 0.1872 = 0.74
    nbinom r             0.1725 / 0.1282 = 1.35     0.1938 / 0.2405 = 0.81
    nbinom agg.psi       0.1782 / 0.1282 = 1.39     0.1874 / 0.2405 = 0.78
    nbinom avg.mu (ctl)  0.0362 / 0.1282 = 0.28     0.1994 / 0.2405 = 0.83

Point B is thin 90 for ordinal and multinomial (R=100) and thin 150 for nbinom
(R=60, the only reduction the budget forced past R=100).

ordinal gamma3: RESOLVED as mixing/noise, not a defect. Two independent things
say so. (1) The flag does not reproduce across streams: the FIRST ordinal R=200
run (the one that also carried the plain-rank p4 artifact, so its stream differs
from rep 1 onward) put gamma3 at ecdfDiff 0.0603, comfortably inside the band,
with the same sampler and the same configuration; gamma3's rank is atomless, so
the tie-break change cannot have altered how it is computed, only which stream
it was computed on. (2) Point B, at 3x the spacing and 3x the burn, drops it to
0.31 of the band with chisqP 0.549. This is the ridge ordinal.md section 9
declares mitigated-not-eliminated, showing up as advertised.

nbinom r and agg.psi: H-MIX, the r-vs-psi ridge, on weaker evidence than A4e's.
At 5x the spacing both flagged channels cross INTO the band (1.35 -> 0.81 and
1.39 -> 0.78, chisqP 0.000 -> 0.001/0.003), which is the shrinking-not-
plateauing signature. Honest caveat: point B ran at R=60, where the noise floor
is high enough that the control avg.mu also reads 0.83, so this is two points
with a wide second one, not A4e's three clean ones - it establishes "the bias
moves the right way under more spacing", not "the stationary distribution is
exactly right". The mechanism is not in doubt: only mu = r exp(psi) is
identified, avg.mu decorrelates at lag 1 and calibrates at 0.28 of the band at
R=200, and the ladder shows r still autocorrelated past lag 200.

multinomial f_ik: NOT resolved by the ladder - reported as a finding, not tuned.
The normalized statistic does not shrink (mean over the three cells 0.89 at
point A, 0.78 at point B, with f.1.1 rising to 1.08 and keeping chisqP 0.000 at
both points), and the mixing story was already weak because every f cell's ACF
clears 0.1 by lag 9 at thin = 1. The aggregate softmax probabilities are clean
at BOTH points (0.044-0.086), so whatever this is, it is confined to the raw
per-forest level and does not reach the reported deliverable. Pre-registered
suspect, unchanged: MultinomialForestCombiner::afterCombine's level-centering
draw, whose precision (sum of invV*n) approximates the exact per-leaf
conditional; a U shape is what an approximated precision that is too LARGE would
produce. Routed as a follow-up below rather than diagnosed here.

### Step 4 - the CI matrix

.github/workflows/sbc.yaml becomes a `fail-fast: false` matrix over
discrete-selfcheck + gaussian + the four arms, scheduled and workflow_dispatch
only (unchanged: a statistical test must not gate pull requests). Admission is
the Bonferroni'd band above, computed inside sbc.R (`sbcMatrixAlpha`) and
applied only to the matrix configs, so every previously recorded non-matrix
result still reads at the 5% band it was recorded at. Timeouts are ~3x the
measured laptop wall-clock. Seeds are pinned constants, so a healthy matrix
replays identically and a draw-shifting commit reshuffles the same fixed-seed
stream - the second reason the gate stays non-blocking.

### Recorded commands

    Rscript benchmarks/R/sbc.R discrete-selfcheck
    Rscript benchmarks/R/sbc.R burn-<family> 40000 3
    Rscript benchmarks/R/sbc.R ordinal|nbinom|t|multinom 200 150 30
    # the 5th positional argument is the burn in absolute sweeps, for the ladder
    Rscript benchmarks/R/sbc.R ordinal   200 150  90 108000   # ran at R=100
    Rscript benchmarks/R/sbc.R multinom  200 150  90  18000   # ran at R=100
    Rscript benchmarks/R/sbc.R nbinom    200 150 150  24000   # ran at R=60

The ladder points were reduced to fit the day (R=100 for ordinal and
multinomial, R=60 for nbinom, whose per-rep cost is 20 s even at point B); the
R=200 forms above are the recorded full-R commands for a later pass. Nothing
else was reduced - every verdict run in "Step 3" is R=200.

Neutrality (`Rscript benchmarks/R/sbc.R gaussian 100 200 30`, before vs after
the whole change): chisqP, ksP and ecdfDiff BYTE-IDENTICAL for all seven
functionals, so the ranks replay exactly. The only column that moves is the
band, 0.1323 -> 0.1872, which is the intended Bonferroni. rankUniformity's
faster tabulate-based band was separately checked to reproduce 0.1323 exactly
at alpha = 0.05, so nothing previously recorded is invalidated.

## Verdict

    family        verdict
    ordinal       9/10 PASS at R=200; gamma3 flagged in one stream and is
                  RESOLVED as the cutpoint-vs-mean-level ridge mixing slowly
                  (does not reproduce across streams; 0.31 of the band at 3x
                  the chain length). The cutpoint block, the latent eta and
                  the K category probabilities all calibrate.
    t             ALL PASS (4/4). The lambda scale mixture, the composite-
                  precision sigma draw and the grid nu conditional calibrate.
    nbinom        avg.mu (the identified mean) PASSES cleanly; r and agg.psi
                  flag at thin = 30 and cross into the band at 5x the spacing.
                  Read as the r-vs-psi ridge mixing slowly (H-MIX), on two
                  ladder points rather than three - a stronger statement wants
                  the recorded full-R third point.
    multinomial   aggregate p_k(x*) ALL PASS at both chain lengths; the raw
                  per-forest f_ik cells carry a persistent U that the ladder
                  does not shrink. OPEN, routed below.

## Follow-ups

- multinomial raw f_ik. CLOSED: the pre-registered suspect was the cause, and
  the fix landed as ec2a3d0; the closure is recorded in
  docs/plans/archive/multinomial-level-centering.md (Landing, :177-192, with the
  re-run `sbc.R multinom 200 150 30` showing every functional PASS).
  The one result this pass could not discharge. Evidence:
  three cells, same U shape, chisqP 0.000, statistic 0.86-0.91 of the matrix
  band at thin = 30 and not shrinking at thin = 90, against per-cell ACF under
  lag 10 - so a mixing explanation has to explain why 3x the spacing does not
  help. Pre-registered suspect: MultinomialForestCombiner::afterCombine's
  level-centering precision (sum of invV*n) as an approximation to the exact
  per-leaf conditional. The reported deliverable (the softmax probabilities) is
  unaffected. Next step is an exact-conditional derivation for the centering
  draw checked against the code, not another SBC run.
- The two ridges are mixing-efficiency items, not correctness items, and both
  already have a named remedy shape in their design notes: ordinal.md section 9
  asks for a joint (f-level, gamma) shift move (the centering-move analog) and
  the nbinom r-vs-psi ridge would take the same treatment - a joint rescale of
  (r, psi) - or the interweaving move bcf-ridge-interweaving already proposes
  for the structurally identical (a, mu) ridge.
- Third ladder points (thin ~300) at full R=200 for ordinal, nbinom and
  multinomial, if a stronger monotone-shrinkage statement is wanted than two
  points support.
- Still OUT, unchanged: aft (needs a status setter), hazard and hurdle (design
  depends on y0), heteroscedastic (liftable R-side via setState, deferred),
  monotone (now liftable after 173a710, the first follow-on arm).
