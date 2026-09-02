# Cross-implementation anchor: 0.9-34 (classic engine) vs 1.0-0 (bartcore)

The only independent implementation of dbarts' posterior in reach is the released 0.9-34 tree on
`main`. All three canonical equivalence baselines descend from bartcore's own output; this
comparison does not.

## 1. Setup

Both engines built from clean archives into private libs, neither touching the checkout: `git
archive main | tar -x` -> 0.9-34 and `git archive HEAD` -> 1.0-0 (b102e17c), each `R CMD INSTALL
--preclean -l <lib> .`, every run prefixed with its own `R_LIBS`. Method is the equivalence
harness's STATISTICAL mode, re-implemented to drive two installs: data fixed per scenario, 20 MCMC
seeds per engine, one summary vector per seed, per-summary across-seed Welch z
([[equivalence.R:1830-1841@b102e17c]]; |z| > 4 = FAIL) plus the parameter-free disjoint-seed-range complement
([[equivalence.R:1821-1829@b102e17c]]). Summaries per fit: per-observation posterior mean and sd of yhat.train/yhat.test,
sigma mean/sd, per-variable inclusion proportion, sampled k mean/sd/median, and for rbart_vi
per-group ranef mean/sd and tau mean/sd/median/q10/q90; plus, on top of the harness, Spearman and a
two-sample KS between the two engines' per-observation posterior means. Main arm ndpost 1000 / nskip
1000 / n.trees 75 / n.chains 1 / n.threads 1 / n.thin 1, a second at ndpost 10000 / nskip 2000 over
the matched-prior scenarios, dropping the detection floor ~3x. Harness, recordings and raw compare
output under a private scratch prefix's out/ (anchor.R, compare.R, *.rds,
*.txt).

## 2. Prior and control mapping

Matched by setting both sides explicitly (identical constructors on both trees):

- n.trees 75, power 2.0, base 0.95, sigdf 3.0, sigquant 0.90, n.cuts 100, useQuantiles FALSE,
  n.chains 1, n.threads 1, n.thin 1, combineChains TRUE.
- resid.prior chisq(3, 0.9); `chisq()` and `fixed()` are byte-identical constructors.
- proposal.probs c(birth_death .5, swap .1, change .4, birth .5) and split.probs uniform 1/p:
  spelled differently (1.0-0 promoted the first from NULL to a literal, writes the second `NULL`
  where 0.9-34 writes `1/num.vars`), resolving identically; an explicit non-uniform split.probs
  vector in the two splitprobs arms. k = 2.0 fixed for every gaussian arm and for b_probit_k2.
- BINARY k HYPERPRIOR, the one non-obvious mapping: classic draws k^2 with shape 0.5*(M + 2*nu - 1)
  (main:[[src/dbarts/parameterPrior.cpp:166@b102e17c]]) where bartcore uses 0.5*(M + nu)
  ([[src/bartcore/model.hpp:2521@b102e17c]], the dropped-Jacobian fix bcdcc07/14bd6b52), and the rate term
  including the finite-scale 0.5/s^2 is identical. So classic chi(nu,s) IS bartcore chi(2*nu - 1,
  s): b_probit_chi2 runs chi(1.25, 2) against chi(1.5, 2), b_probit_chiinf chi(1.25, Inf) against
  chi(1.5, Inf). Both hold, which is itself a check on that identity.
- rbart_vi: prior = cauchy on both, rel.scale = sd(y) on both, k = 2.0 explicit; `group.by` given as
  a factor COLUMN of the data, because 0.9-34 resolves an external group.by vector positionally and
  returns a per-observation ranef matrix instead of a per-group one, so an external vector is not a
  comparable surface there.
- Numeric matrices throughout, so 1.0-0's `factors`/`missing` defaults are never reached.
- Where the default differs by design, 0.9-34 is set to 1.0-0's value: binary k, 0.9-34 chi(1.25,
  Inf) == chi(1.5, Inf) against 1.0-0 chi(1.5, 2.0) (5b6e4825), held equal in b_probit_chi2/chiinf,
  with each engine's OWN default compared separately (E5).
- Not matchable, compared as-is and read as designed: xbart n.burn (three elements on 0.9-34, two on
  1.0-0, which never carries a chain across folds) run c(200,150,150) against c(200,150); xbart fold
  cut grids (per-fold vs shared); rbart_vi Gibbs (R loop vs in-core).

## 3. Scenario results

n>4 / ndisj are counts over the scenario's summaries; (10x) is the high-precision arm.

  scenario           entry / family          summ  max|z|  n>4  ndisj  verdict
  g_cont             bart2 gaussian+test     1012    4.44    1     0   AGREE
  g_mixed            bart2, coarse noise x   1012    3.50    0     0   AGREE
  g_weights          bart2 + weights         1012    4.96    2     0   AGREE
  g_offset           bart2 + offset/test     1012    3.91    0     0   AGREE
  g_bart             bart() shim             1012    3.38    0     0   AGREE
  b_probit_k2        bart2 probit, k=2       1010    3.58    0     0   AGREE
  b_probit_chi2      probit, chi(.,2)        1013    4.22    1     0   AGREE
  b_probit_chiinf    probit, chi(.,Inf)      1013    3.73    0     0   AGREE
  b_probit_DEFAULTS  each engine's own k     1013    3.42    0     0   AGREE (E5)
  rbart              rbart_vi, Friedman       857    3.45    0     0   AGREE
  (10x) g_cont/g_bart/b_probit_k2/rbart max|z| 3.64/3.81/5.52/3.33, n>4 0/0/2/0  AGREE
  g_quants           useQuantiles=TRUE       1012   25.77   46     9   DIFFER E1
  g_splitprobs       non-uniform split.probs 1012   12.93   11     2   DIFFER E1
  g_quantsplit       both                    1012   14.03   12     4   DIFFER E1
  g_zeroweights      bart2, 80/400 w = 0     1012     n/a  n/a     0   DIFFER E2
  rbart_sym          rbart_vi, symmetric y    857    6.06   11     0   DIFFER E3
  xbart              5-fold CV, rmse            7   54.75    7     7   DIFFER E4

Z CALIBRATION, the headline. Pooling every summary of the nine AGREE scenarios at ndpost 1000: n =
8953, mean z +0.044, sd z 0.975, |z| > 2 in 4.33% (nominal 4.55), |z| > 3 in 0.346% (0.27), 4
summaries over 4 against 2.53 expected under t_38, and no disjoint seed range anywhere. The 10x arm:
n = 3891, mean -0.047, sd 1.088, 2 over 4 against 1.10 expected (shoulder mildly fat, tail nominal).
Where the priors are matched, the z field is standard normal: the two engines sample the same
posterior. Rank/KS on the per-observation posterior means, every AGREE scenario: Spearman
0.99822-0.99996, KS D 0.010-0.040 at p = 1.000, mean(1.0-0 minus 0.9-34) between -0.0143 and +0.0173
on responses of mean 14.4 (probit latents -0.0024 to +0.0002).
DETECTION FLOOR at 10x: a per-observation posterior mean resolves to 0.185 (1.3% of the
across-observation sd), the aggregate mean fitted value to 0.009%, mean posterior sd to 2.6%, sigma
to 3.0-3.6%, rbart tau to 2.0%. Measured gaps: sigma -0.66% (g_cont, z 0.88), -1.37% (g_bart, z
1.51), -1.07% (rbart, z 1.26), tau -0.03% (z 0.05). A systematic sigma shift above ~3% is excluded
on the matched arms, a 1-2% one is not.

## 4. Explained differences

E1. CHANGE-MOVE DETAILED BALANCE. 0.9-34's changeRule.cpp accepts on the pure pi-ratio with no
proposal-density term, inherited from the CGM-lineage original; 1.0-0 adds the hybrid
correction ([[docs/design/change-move-balance.md:18-77@658869ac]], [[NEWS.Rd:1154-1163@658869ac]]; exact gate
benchmarks/R/change-balance.R, defect z +255 -> repair z -1.2,
[[docs/plans/archive/change-move-fix.md:200-207@658869ac]]). The doc's bias factor, [p_var(v')/p_var(v)] *
[|Valid(v)|/|Valid(v')|], is invisible at equal cut counts AND equal split probabilities; both
halves reproduce, each in its own scenario.
- g_quants (coarse columns get 2-4 quantile cuts, the rest ~100): 0.9-34 OVER-uses the
  low-cardinality variables - vprop x6-x10 0.109/0.086/0.146/0.082/0.115 against 1.0-0's
  0.053/0.056/0.068/0.070/0.065, x1-x5 correspondingly under-used. max |z| 25.8, 8 of 10 vprop cells
  with disjoint seed ranges, total-variation distance 0.2245. It reaches the fit: sigma 1.0758 vs
  0.9253 (-14.0%, z 15.9, disjoint), 11 of 400 obs means over |z| 4.
- g_splitprobs (all continuous, equal cut counts, split.probs 8:8:4:4:2:1:1:1:1:1): 0.9-34 instead
  over-uses the HIGH-probability variables - x1 0.3218 against the prescribed 0.25 (1.0-0: 0.2815),
  x9 0.0048 against 0.031 (1.0-0: 0.0132), the documented effectively-squared rule prior. max |z|
  12.9, TV 0.0740, and here sigma does NOT move (-0.73%, z 0.74): a tree-structure-only distortion.
  The fitted-value RANKING survives throughout - g_quants still gives Spearman 0.9987, KS D 0.0175
  at p = 1.000.

E2. ZERO-WEIGHT SIGMA DF. 0.9-34 counts all n rows in the sigma posterior df
(main:[[src/dbarts/parameterPrior.cpp:109@658869ac]]) though zero-weight rows contribute nothing to the sum of
squares; 1.0-0 counts positive-weight rows only (cf99a00, [[NEWS.Rd:1164-1172@658869ac]]). With 80 of 400 weights
exactly 0 and true sigma 1: 1.0-0 gives 1.072 (sd 0.078) on every seed, 0.9-34 gives 0.43-0.54 on 12
of 20 seeds and NON-FINITE yhat/sigma on the other 8, leaving 1002 of its 1012 summary columns
unusable. The documented feedback loop (deflated sigma -> trees absorb residual -> smaller S ->
smaller sigma) runs here to numerical breakdown; the NaN is not itemized in NEWS. Only vprop stayed
finite, and it too diverges (|z| 4.98).

E3. rbart_vi RESPONSE SCALING. 0.9-34 re-anchors the response transform during warmup with the
current random effects subtracted; 1.0-0's in-core Gibbs fixes it at creation from y alone, so b
never touches it - measured as "sigma ~1.4% lower in-core" at
[[docs/design/grouped-random-effects.md:212-219@658869ac]] ([[NEWS.Rd:521-537@658869ac]]). Reproduced with the right
conditional structure: the effect scales with ranef sd relative to response range. rbart_sym (range
~8, ranef sd 1.0) gives sigma 0.8256 vs 0.8109, -1.78%, z 4.72, and mean posterior sd of yhat
+3.13%, z -4.11 (11 of 857 over 4); rbart on Friedman (range ~30, ranef sd 1.5) gives -1.07%, z 1.26
at 10x, same sign, below floor. tau, ranef means and vprop agree in both (10x tau -0.03%, z 0.05):
the shift is in the leaf prior's effective scale, not the variance components.

E4. xbart FOLD LEAKAGE. 0.9-34 carries the chain across folds, so a fold's held-out rows were
training rows moments earlier; 1.0-0 starts fresh per split and shares one cut grid across folds
([[NEWS.Rd:618-638@658869ac]]). 5-fold rmse over 20 reps: 0.9-34 1.5537 against 1.0-0 1.7170, i.e. 0.9-34 reports
10.5% LOWER loss, z -54.8, all 7 summaries with disjoint seed ranges; 14.7% at n.trees = 25 against
6.4% at 75, the documented concentration at small ensembles. A discriminator separates leakage from
under-burning: raising burn-in to 2000 leaves 1.0-0 unmoved (1.7574 -> 1.7576, so its 200-sweep
fresh burn is adequate) while 0.9-34's loss RISES 1.562 -> 1.691 as the carried state washes out.
Most of the gap is leakage; the 1.691 vs 1.758 residual is carry-over plus the per-fold cut grid.

E5. BINARY k DEFAULT. 0.9-34 defaults binary k to chi(1.25, Inf) (== chi(1.5, Inf)), 1.0-0 to
chi(1.5, 2.0) (5b6e4825, docs/plans/archive/chi-default-research.md). On the strong-signal probit scenario
it is invisible (b_probit_DEFAULTS max |z| 3.42, 0 over 4, sampled k 2.196 vs 2.119), correctly,
since the research names weak-signal balanced cells as the regime. A dedicated
weak-signal cell (n = 200, base rate 0.5, effect 0.3) reproduces the documented failure: 0.9-34's
sampled k has median 14848 and max 519000 against 1.0-0's 3.61 and 7.01, and the runaway collapses
the fit - mean |posterior mean latent| 0.00046 vs 0.2591, mean posterior sd 0.0098 vs 0.3962 (z ~
-47), i.e. 0.9-34's default binary fit degenerates to the intercept. The k channel gives only z
1.2-1.7: the runaway inflates across-seed variance, the degeneracy [[equivalence.R:1815-1820@658869ac]] warns
about.

## 5. Unexplained disagreements

NONE. Every |z| > 4 failure is accounted for by a change recorded in NEWS.Rd or docs/design, in the
direction and under the conditions those documents state; every scenario whose priors could be
matched agrees with a standard-normal z field. Two sub-threshold observations, so a later pass need
not re-derive them:
1. 0.9-34's sigma is the larger in all four matched gaussian/grouped arms, by 0.66-1.78%. Only
  rbart_sym clears |z| 4 (E3); the other three sit at z 0.88-1.51 against a ~3% floor. Consistent
  with E3 where it applies and with the change-move residual the design doc concedes even at equal
  ROOT cut counts (CONTROL z -15.7 at exact-enumeration scale, docs/plans/archive/change-move-fix.md), but
  not separated from Monte Carlo error here.
2. b_probit_k2's mean per-observation posterior sd differs by 0.16-0.18% against a 0.21% floor,
  reaching z 3.02 at 10x, but changes SIGN between the two precision arms, so it is not systematic.

## 6. What this anchor cannot reach

0.9-34 has no `family` argument, so logistic, multinomial, ordinal, negbin, hurdle, aft, hazard and
heteroscedastic have no independent implementation here; nor do multi-forest/BCF, DART, monotone,
MIA missingness, categorical subset splits, blocks, bases, storage = "single", the mid-chain
mutation API, or any weighted binary fit (0.9-34 fit a weighted probit 1.0-0 refuses as incorrect).
Threading and multi-chain axes were held at one. Covered: gaussian (plain, weighted, zero-weighted,
offset, test set, both wrapper entries), probit (fixed k and both chi parameterizations), rbart_vi
grouped effects, xbart cross-validation - the families both engines ship, and no more.

## 7. Recommendation for benchmarks/baselines/MANIFEST

ADD an anchored row, but not as a re-recordable baseline and not in equivalence.R's format:

- Record the 0.9-34 side once. Unlike equivalence-5430fdb.rds (the one `historical-classic` row,
  whose engine is deleted and can never be re-run) this one CAN be regenerated from `git archive
  main` - the only oracle in the tree that is external to bartcore AND still executable.
- Give it its own role token, `external-anchor`, distinct from `historical-classic` (permanent,
  dead) and `current` (a live gate): not compared on push, never bitwise, statistical only, and
  one-sided - a FAIL is a question, not a regression, because 0.9-34 is known-wrong on E1/E2/E4 and
  known-different on E3/E5.
- The row must carry the EXPECTED-DIFFERENCE table, not a max |z|. Its value is naming which
  scenarios must AGREE (the ten rows above, at |z| > 4 with no disjoint range) and which must
  DISAGREE with a stated sign and size (g_quants vprop toward low-cardinality variables;
  g_splitprobs vprop toward high-probability ones; g_zeroweights sigma deflated to NaN; rbart_sym
  sigma -1.8%; xbart loss -10.5%; weak-signal binary k runaway). A change that made g_quants AGREE
  would mean the change-move fix had been undone - which no bartcore-descended baseline can detect,
  since they all inherit it.
- Cadence: schedule-only, beside the existing equivalence workflow, or on demand before a release;
  two installs and ~35 minutes of one core, no quiet machine needed. Keep the harness at
  `benchmarks/R/anchor-093.R` next to equivalence.R with sec 2's mapping in its header - above all
  the classic-chi(nu) == bartcore-chi(2nu-1) identity, the one mapping a later reader cannot
  rediscover from the R surface.

If only one thing is adopted, adopt the expected-difference table: the current baselines are
self-referential because nothing in the tree records what the released engine's posterior looked
like, and sections 3-4 are that record - cheap, because `main` still builds.
