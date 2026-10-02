# Review 3, cross-implementation anchor: 0.9-34 vs 1.0-0 at 01dee4b4

Lens: anchor. A re-run of docs/plans/review-2026-08-24/anchor-main.md (then b102e17c) at tip
01dee4b4, against released 0.9-34. The question was whether any posterior has changed, in a model
0.9-34 can also fit, in a way no recorded decision explains.

Covered: all 17 of the record's scenarios plus 11 new ones (gaussian, probit, rbart_vi, xbart; front
doors bart/bart2, bartBT/bart, rbart_vi, xbart). A 10x-precision arm. Change-move-off controls for
every residual flag. A re-run of the in-tree release gate benchmarks/R/classic-compare.R (dec-B120).
Not covered: anything 0.9-34 cannot fit (families, multi-forest, DART, monotone, MIA, subset splits,
sparse), threading beyond one thread, the sampler mutation API beyond what classic-compare runs,
predict-from-saved-trees.

Verdict: **no unexplained posterior disagreement.** Every |z| > 4, disjoint range, KS or aggregate
anomaly traces to a recorded decision. In most cases a control removes it, or a fresh seed block
clears it. One expected difference from the 08-24 record (E3, rbart_vi response scaling) has gone,
as the rbart_vi port predicts. Two incidental findings are not posterior changes: anchor-01 (the
release gate's 0.9-34 side no longer records) and anchor-02 (`fn <- bart; fn(...)` fails).

## 1. Setup

- 0.9-34: `git archive main` (cb290550, 0.9-34 plus protection fixes, no sampler change). 1.0-0:
  `git archive HEAD` of the pinned tree (01dee4b4). Both built with `R CMD INSTALL --preclean -l` into
  scratchpad r3-anchor-libold / r3-anchor-libnew. Every R process ran with exactly one lib on R_LIBS.
  At most 4 R processes ran at once, each with n.threads 1.
- The original harness was gone, so it was rebuilt from the record's description:
  scratchpad/r3-anchor-out/anchor.R (one scenario, one engine, 20 seeds, one summary vector per seed)
  and compare.R / aggregate.R. Statistics are those of equivalence.R's statistical mode:
  - per-summary Welch z (|z| > 4 is a FAIL);
  - disjoint seed ranges;
  - Spearman and two-sample KS between the engines' per-observation posterior means;
  - Welch z on per-seed aggregates (mean |fitted|, mean posterior sd, sigma, k, tau).
- Summaries per fit: posterior mean and sd of each yhat.train/yhat.test, sigma mean and sd,
  per-variable inclusion proportion, and k mean/sd/median where k is sampled. rbart_vi adds per-group
  ranef mean and sd, and tau mean/sd/median/q10/q90. xbart: rmse per grid cell, averaged over reps,
  plus the overall mean.
- Data: Friedman, n = 400, p = 10, 100 test rows, noise sd 1. Probit: y ~ Bern(pnorm(scale(f))).
  Weak-signal binary cell: n = 200, effect 0.3. rbart: 20 groups of 20, ranef sd 1.5 (Friedman) or
  1.0 (symmetric response, range about 8). Data seeds are fixed per scenario; MCMC seeds are 1001-1020,
  and 1101-1120 for the "fresh block" re-checks.
- Main arm: ndpost 1000, nskip 1000, n.trees 75, 1 chain, n.thin 1. The 10x arm uses ndpost 10000
  and nskip 2000.

## 2. Prior and control mapping

The 08-24 record's mapping holds, spelled in 1.0-0's current vocabulary:

| 0.9-34 bart2 | 1.0-0 bart |
| --- | --- |
| sigdf 3, sigquant 0.9 | `family = gaussian(sigma = chisq(3, 0.9))` |
| power 2, base 0.95, split.probs | `tree.prior = cgm(2, 0.95[, split.probs])` |
| proposal.probs c(.5, .1, .4, birth .5) | `control = dbartsControl(proposal.probs = c(birth_death .5, swap .1, change .4, perturb 0, rule_gibbs 0, birth .5))` |
| k 2 | k 2 |

The proposal mixture now has to be pinned on the 1.0-0 side: its default drops swap (dec-B02).

- The chi mapping holds: 0.9-34 `chi(1.25, s)` corresponds to 1.0-0 `chi(1.5, s)` (dec-B04). It is
  confirmed again at s = 2 and s = Inf (probit) and, new here, at s = 2 for a continuous response.
- rbart_vi: prior cauchy and k = 2 explicit on both. 1.0-0 takes the mixture through `...` into
  dbartsControl. 0.9-34's rbart_vi takes no mixture, so the change-off control patches the default of
  0.9-34's `dbarts()` with assignInNamespace (section 5, rbart_sym).
- xbart: n.burn is c(200, 150, 150) on 0.9-34 and c(200, 150) on 1.0-0, and these cannot be matched
  (dec-B77). n.trees c(25, 75) x k c(1, 2, 4), 5-fold, 4 reps.
- "DEFAULTS" scenarios pin nothing but the draw counts, the chain count and threads = 1.

## 3. Results

summ is the number of summaries; n>4 and ndisj count summaries over |z| 4 and with disjoint seed
ranges. Main arm unless noted.

| scenario | what | summ | max abs z | n>4 | ndisj | verdict |
| --- | --- | --- | --- | --- | --- | --- |
| g_cont | bart/bart2 gaussian + test | 1012 | 3.13 | 0 | 0 | AGREE |
| g_mixed | coarse x6-10, uniform grid | 1012 | 3.95 | 0 | 0 | AGREE |
| g_weights | weights 0.5/1/2 | 1012 | 3.66 | 0 | 0 | AGREE |
| g_offset | offset + test offset | 1012 | 3.63 | 0 | 0 | AGREE |
| g_bart | bartBT vs 0.9-34 bart, pinned | 1012 | 2.97 | 0 | 0 | AGREE |
| b_probit_k2 | probit, k = 2 | 1010 | 3.43 | 0 | 0 | AGREE (E1 residual at 10x, sec 5) |
| b_probit_chi2 | chi(1.25,2) vs chi(1.5,2) | 1013 | 3.05 | 0 | 0 | AGREE |
| b_probit_chiinf | chi(1.25,Inf) vs chi(1.5,Inf) | 1013 | 2.99 | 0 | 0 | AGREE |
| b_probit_DEFAULTS | each engine's own k, strong signal | 1013 | 2.94 | 0 | 0 | AGREE (E5 inert here) |
| rbart | rbart_vi, Friedman | 857 | 3.44 | 0 | 0 | AGREE |
| rbart_sym | rbart_vi, symmetric y | 857 | 3.16 | 0 | 0 | AGREE (was E3) |
| g_quants | useQuantiles, coarse x6-10 | 1012 | 22.82 | 65 | 6 | DIFFER E1 |
| g_splitprobs | split.probs 8:8:4:4:2:1x5 | 1012 | 11.57 | 9 | 2 | DIFFER E1 |
| g_quantsplit | both | 1012 | 12.23 | 10 | 1 | DIFFER E1 |
| g_zeroweights | 80/400 w = 0 | 1012 | n/a | n/a | 0 | DIFFER E2 |
| xbart | 5-fold rmse grid | 7 | 27.84 | 6 | 6 | DIFFER E4 |
| b_weak_DEFAULTS | weak-signal binary, own k | 513 | 59.07 | 488 | 457 | DIFFER E5 |
| **new** g_DEFAULTS | bart vs bart2 at defaults, 4 chains (swap 0 vs 0.1) | 1012 | 3.94 | 0 | 0 | AGREE |
| **new** g_bartBT_DEFAULTS | bartBT vs 0.9-34 bart at defaults | 1012 | 3.53 | 0 | 0 | AGREE |
| **new** g_chi2 | gaussian chi(1.25,2) vs chi(1.5,2) | 1015 | 3.29 | 0 | 0 | AGREE |
| **new** b_offset | probit + offset/test offset | 1010 | 3.96 | 0 | 0 | AGREE |
| **new** g_factor | factor column, factors = "indicators" | 1016 | 3.16 | 0 | 0 | AGREE |
| **new** rbart_bin | binary rbart_vi, k = 2 | 855 | 3.58 | 0 | 0 | AGREE |
| **new** b_zeroweights | probit, 80/400 w = 0 | 1010 | n/a | n/a | 0 | DIFFER E6 |
| **new** g_zwsub2 | 1.0-0 w = 0 vs 0.9-34 rows dropped | 852 | 3.75 | 0 | 0 | AGREE |
| **new** b_zwsub2 | same, probit | 850 | 3.29 | 0 | 0 | AGREE |
| **new** b_zwsub | as b_zwsub2, grids differ | 850 | 5.52 | 1 | 0 | DIFFER by construction (sec 5) |
| (10x) g_cont / g_bart / b_probit_k2 | | | 3.29 / 3.66 / 3.98 | 0 / 0 / 0 | 0 | AGREE |
| (10x) rbart / rbart_sym | | | 3.95 / 4.38 | 0 / 2 | 0 | AGREE / E1 residual |
| (10x) g_DEFAULTS / g_zwsub2 | | | 4.22 / 4.35 | 3 / 1 | 0 | AGREE (fresh block: g_DEFAULTS 3.62, 0 over 4) |

On the matched rbart_bin and b_zwsub2 fits the posterior agrees with 0.9-34's. With all rows weighted
(rbart_bin) or the zero-weight rows dropped on the 0.9-34 side (b_zwsub2), 0.9-34 never reaches the
weighted-probit path whose output in b_zeroweights is all non-finite (E6).

**Z calibration, main arm.** The 16 AGREE scenarios pooled: n = 15,888, mean z -0.058, sd 1.002,
|z| > 2 in 4.59% (t38 nominal 5.27), |z| > 3 in 0.34% (0.48), none over 4 (4.5 expected under
independent t38; the summaries are correlated). No disjoint range anywhere.

- Spearman on per-observation posterior means: 0.99968-0.99991.
- KS D: 0.010-0.0175, p = 1.000.
- mean(1.0-0 minus 0.9-34) of fitted values: within ±0.0014 on the gaussian responses (mean 14.5)
  and the probit latents.
- The rbart_vi tree fits are the exception, at -0.027 (rbart, mean 15.1) and +0.008 (rbart_sym).
  Their aggregate z is 1.21 and -0.43, and at 10x they shrink to -0.005 and +0.003, so this is Monte
  Carlo error in how the intercepts and the trees split the fit.

**Z calibration, 10x arm.** Seven scenarios, n = 6,612: mean +0.036, sd 1.091, 6 summaries over 4
(1.9 expected). The fat shoulder and the six exceedances resolve as follows:

- The 3 in g_DEFAULTS clear on a fresh seed block (max 3.62).
- The 2 in rbart_sym are the E1 residual (section 5).
- The 1 in g_zwsub2 stands alone, and its aggregates are clean.

The 08-24 record also saw this shoulder (sd 1.088).

## 4. Explained differences

These are the 08-24 record's E1-E5, re-measured, plus E6.

**E1. Change move without detailed balance in 0.9-34.** [inst/NEWS.Rd:376](https://github.com/vdorie/dbarts/blob/01dee4b4f1a21565c91eceea391e666851662a0d/inst/NEWS.Rd#L376); dec-B03; the design note is
docs/design/change-move-balance.md.

- **g_quants.** 0.9-34 over-uses the coarse x6-x10: vprop 0.130/0.097/0.067/0.134/0.098 against
  1.0-0's 0.068/0.059/0.047/0.069/0.052. TV distance 0.233. Sigma 0.991 vs 0.845 (-14.8%, z 12.1).
- **g_splitprobs.** 0.9-34 over-uses the high-probability variables: x1-x2 0.317/0.327 against
  0.286/0.285, prescribed 0.25. TV 0.074.
- **New control: the change move off on both sides** (birth_death 0.9, swap 0.1, change 0).
  - g_quants, g_splitprobs and g_quantsplit all agree: max |z| 3.53 / 3.82 / 3.50, 0 over 4, no
    disjoint range, sigma z -0.70 / 0.72 / 0.02.
  - 1.0-0's vprop with the change move on equals its vprop with it off (g_quants x1-x4
    0.172/0.166/0.149/0.129 vs 0.171/0.173/0.144/0.125). The biased side is 0.9-34's.
- **Sigma moves in g_splitprobs here** (-3.4%, z 3.44), where the 08-24 record saw no shift (z 0.74).
  With the change move off the shift disappears, so it is E1 at this design's strength and no new
  cause.

**E2. Zero-weight rows in the sigma degrees of freedom.** [inst/NEWS.Rd:384](https://github.com/vdorie/dbarts/blob/01dee4b4f1a21565c91eceea391e666851662a0d/inst/NEWS.Rd#L384); dec-A11.

- g_zeroweights: 0.9-34 is non-finite on 3 of 20 seeds. On the 17 finite seeds its sigma is 0.375
  against 1.0-0's 0.786 (-52%, z -34) and its posterior sd is 0.370 against 0.630.
- vprop (finite on every seed) z 4.19.
- The anchor that settles it: 0.9-34 fitting only the positive-weight rows agrees with 1.0-0's
  weighted fit. In g_zwsub2 sigma is 0.711 vs 0.737 at the main arm and 0.6875 vs 0.6821 (z 0.82)
  at 10x, and every summary agrees.

**E3. rbart_vi response scaling: gone.** The 08-24 cause was the in-core grouped Gibbs. That code has
been removed (dec-A01), and rbart_vi is now 0.9-34's R loop on the new sampler (dec-B130, dec-A99),
including `setOffset(updateScale = TRUE)` during warmup.

- Main arm, rbart_sym: sigma 0.9077 vs 0.9094 (z -0.68), against the record's -1.78% (z 4.72).
- The small residual left at 10x is E1 (section 5).
- The expected-difference row "rbart_sym sigma -1.8%" in anchor-main.md section 7 is now wrong:
  rbart_sym must AGREE.

**E4. xbart fold leakage in 0.9-34.** [inst/NEWS.Rd:177](https://github.com/vdorie/dbarts/blob/01dee4b4f1a21565c91eceea391e666851662a0d/inst/NEWS.Rd#L177); dec-A05, dec-B77. Mean 5-fold rmse is 1.284 vs
1.429: 0.9-34 reports a 10.2% lower loss, z -27.8, with 6 of 7 summaries disjoint. Every cell is
lower on 0.9-34. The record measured 10.5%.

**E5. Binary k default.** 0.9-34 uses chi(1.25, Inf), 1.0-0 chi(1.5, 2) ([inst/NEWS.Rd:138](https://github.com/vdorie/dbarts/blob/01dee4b4f1a21565c91eceea391e666851662a0d/inst/NEWS.Rd#L138); dec-A07,
dec-A13).

- On b_weak_DEFAULTS, 0.9-34's sampled k has median 37,921 against 3.92. Its fit collapses to the
  intercept: mean |posterior mean latent| 0.0076 vs 0.232, mean posterior sd 0.0247 vs 0.378
  (aggregate z -34).
- The k channel itself reads only z 1.5-1.7, the degenerate case the harness warns about. The
  disjoint channel catches it: 457 summaries.
- The strong-signal b_probit_DEFAULTS agrees (k 1.91 vs 1.84, z 1.23).

**E6 (new). Weighted probit with zero weights.** 0.9-34 returns non-finite yhat on all 20 seeds of
b_zeroweights. 1.0-0 reads 0/1 weights as the row mask ([inst/NEWS.Rd:80](https://github.com/vdorie/dbarts/blob/01dee4b4f1a21565c91eceea391e666851662a0d/inst/NEWS.Rd#L80), "0.9-34 fit a weighted probit,
which was incorrect"; dec-B13).

- vprop, the only finite channel, agrees (max |z| 2.77).
- 1.0-0's masked fit agrees with 0.9-34 fitting the kept rows (b_zwsub2).
- 0.9-34's NaN itself is not itemized in NEWS. The same was true of the gaussian NaN in the 08-24
  record.

## 5. Flags traced, sub-threshold included

- **b_probit_k2, 10x: 1.0-0's mean posterior sd is larger by about 0.2-0.3%.**
  - Measured: +0.17% (aggregate z -3.34; test +0.24%, z -3.72), and on a fresh block +0.29%
    (z -5.08; test +0.32%, z -4.50). It has the same sign as the main arm (+0.11%).
  - This is the 08-24 record's sub-threshold observation 2, which then changed sign. It no longer
    does.
  - The control: with the change move off on both sides it vanishes, at -0.04% (z 0.55) and +0.06%
    (z -1.05) on two blocks. 1.0-0's change-on value (0.5211) equals both engines' change-off values
    (0.5206-0.5211), while 0.9-34's change-on value is lower (0.5196-0.5203).
  - Explained as E1: the change-move bias is not confined to unequal root cut counts. Valid-cut counts
    differ by variable at every node below the root, which is the residual change-move-balance.md
    concedes.
- **rbart_sym, 10x: 1.0-0's sigma is lower by 0.32-0.41%** (z 3.41, then 2.11 on a fresh block), and
  2 per-observation means sit over 4.
  - With the change move off on both sides (0.9-34 patched through its `dbarts()` default):
    -0.13% (z 1.07) and -0.05% (z 0.45). 0.9-34's sigma drops from 0.9112 to 0.9093/0.9088, and
    1.0-0's stays at 0.908.
  - Explained as E1, and not scaling (E3). This also accounts for the 08-24 record's sub-threshold
    observation 1, that 0.9-34's sigma was the larger in every gaussian or grouped arm.
- **g_bart, 10x: sigma -2.38% (z 3.21).** Cleared on a fresh block: -0.37%, z 0.40. Monte Carlo.
- **g_DEFAULTS.** The 10x arm had 3 summaries over |z| 4 and the main arm a test-sd aggregate at
  z -3.05. A fresh 10x block gives max 3.62, 0 over 4, aggregates within |z| 1.65. Monte Carlo. Swap
  off (dec-B02) leaves the posterior where it was, matching classic-compare's Table 2.
- **b_zwsub: |z| 5.52, then 4.71 / 4.13 on a fresh block, with the same rows recurring.**
  - Cause: the cut grid. It spans every row, zero-weight rows included, on both engines. Dropping the
    rows in 0.9-34 and masking them in 1.0-0 therefore gives different grids.
  - Control within 1.0-0 (rows dropped vs rows masked): the same signature, max 5.63 on the same rows
    (61, 180). 0.9-34-dropped vs 1.0-0-dropped: max 3.61, 0 over 4. b_zwsub2, whose zero rows avoid
    every column's extremes so the grids coincide, agrees.
  - A difference built into the scenario, not an engine change.

## 6. The in-tree release gate (classic-compare.R, dec-B120)

- **Compare at tip.** Recorded at 01dee4b4 (mixture=classic) and compared with the stored
  benchmarks/baselines/classic-compare-0.9-34.rds.
  - 22 of 26 scenarios are clean. Every per-scenario maximum matches docs/plans/classic-compare.md
    Table 1 to the second decimal: friedman 2.63, offset 4.20, mixedcuts 6.44, zeroweights 47.28.
    1.0-0's draws for those scenarios have not moved since the table was taken.
  - xbart is 53.17 against 48.94 in the table, with the same 7 disjoint cells. Still E4.
- **Re-recording the 0.9-34 side fails.** This is anchor-01.

## Findings

**anchor-01, MAJOR (release process; no user-facing effect). The 0.9-34 side of the classic-compare
release gate no longer runs.**

- Location: benchmarks/R/classic-compare.R, `makeSampler` (line 523) and `fitViaXbart` (line 666).
- Claim: since d89c9405 (2026-09-14) these pass `family = gaussian(sigma = chisq(3, 0.9))`
  unconditionally, so `R_LIBS=<0.9-34> Rscript classic-compare.R record` errors in 6 of 26
  scenarios. That breaks the reproduction recipe in docs/plans/classic-compare.md, which also claims
  the 0.9-34 recording "reproduces bitwise when it is [re-taken]". It also breaks
  CLASSIC_COMPARE_SEED_OFFSET re-checks: a different block refuses to compare against the stored
  baseline, so the 0.9-34 side must be re-recorded. The TODO's submission battery depends on that
  re-check when a scenario flags.
- Probe:

  ```
  cd <git archive of 01dee4b4>
  CLASSIC_COMPARE_SCENARIOS=offset,chik,gibbsloop,setpredictor,xbart,xbart1rep \
    CLASSIC_COMPARE_CORES=1 R_LIBS=<0.9-34 lib> \
    Rscript benchmarks/R/classic-compare.R record x.rds quick
  -> Error in dbarts(scn[["x"]], scn[["y"]], ...):
       unused argument (family = gaussian(sigma = chisq(3, 0.9)))
  ```

  The full record at CORES = 4 gives "fit failed for offset 1, ..., setpredictor 20, ...
  (120 function calls resulted in an error)".
- Why the gates missed it: the gate runs only at submission, and nothing runs it between releases
  (dec-B120). Compare mode reads only the stored 0.9-34 baseline. The spelling sweep edited the
  harness without re-recording the 0.9-34 side.
- Fix: spell `resid.prior = chisq(3, 0.9)` when `packageVersion("dbarts") < "1.0"`, as the harness
  already branches bart/bartBT. Then re-record the 0.9-34 side once and confirm it is bitwise equal to
  the stored rds.

**anchor-02, MINOR (incidental, outside this lens). Calling `bart` through an alias named `fn`
fails.**

- Location: R/utility.R, `redirectCall`, `originalFn <- eval(call[[1L]])`.
- Claim: `redirectCall` evaluates the call's head symbol in its own frame, where `fn` is its formal
  (the redirect target, `dbartsControl`). Because 1.0-0's `dbartsControl` gained `...`, the
  dots-forwarding branch then keeps every `bart` argument that `dbartsControl` lacks, and
  `dbartsControl` refuses them. 0.9-34's `dbartsControl` had no `...`, so the same call worked; this
  is an undocumented regression for this spelling.
- Probe (R_LIBS = the 1.0-0 lib):

  ```
  fn <- bart; fn(y ~ ., data = df, n.samples = 10L, n.burn = 10L, n.chains = 1L, verbose = FALSE)
  -> unused arguments 'formula', 'data', 'family', 'factors' passed to 'dbartsControl'
  ```

  Under 0.9-34, `fn <- bart2; fn(...)` gives "ok". The aliases h, f, fit, fun and g all work.
  `call <- bart` works but warns "formals(fun): argument is not a function". An alias local to a
  function fails on both releases ("object 'myfit' not found"), so that part is not a regression.
- Why the gates missed it: tests call the entry points by name.
- Fix: evaluate the head in the caller's frame (pass `sys.function(sys.parent())` or
  `parent.frame(2L)` from the entry point) rather than in `redirectCall`'s own frame.

## Checked and found correct

- Matched-prior agreement of bart/bart2 for gaussian (plain, weighted, offset, test set, coarse
  predictors, indicator-coded factor) and probit (fixed k, chi at s = 2 and s = Inf, offset).
- bartBT reproduces 0.9-34's bart, pinned and at defaults: the compatibility door keeps BayesTree's
  mixture and defaults.
- Front-door defaults: gaussian bart vs bart2, 4 chains. Swap off leaves the posterior where it was.
- The chi(nu) relabelling identity, for continuous and binary responses.
- The rbart_vi port: continuous (two designs) and binary, ranef, tau and sigma at 10x.
- Zero weights mean dropped rows in 1.0-0: gaussian and probit, against 0.9-34 on the subset.
- E1 lives entirely in 0.9-34: in every design tried, 1.0-0's posterior is invariant to the
  change-move share.
- classic-compare at tip: 1.0-0 draws unchanged since Table 1, the 22 clean scenarios still clean.

Artifacts: scratchpad r3-anchor-out/ (anchor.R, compare.R, aggregate.R, main/, x10/, x10s100/,
x10nc/, x10ncs100/, nc/, zwgrid*/, cc/; compare.txt and aggregate.txt in each).
