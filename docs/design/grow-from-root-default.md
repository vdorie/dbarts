# grow-from-root-default: should the XBART warm start default on

Status: KILLED (measured), 2026-08-08. TODO `grow-from-root-default: RESOLVED
NO`. This is the durable record of the pre-registered study and its result;
working papers, the harness, and 98 checkpoints live in
an untracked session directory (gitignored, not preserved).
Pre-registration: `docs/plans/archive/grow-from-root-default-study.md` (its Results
section is the short form; this doc carries the full data and does not
contradict it - see the Provenance section for one numeric discrepancy
between the two that is reported, not resolved). Surface under study:
`docs/plans/archive/grow-from-root-warm-start.md` (LANDED 2026-07-10, opt-in
`n.grow.sweeps`). Validity floor: `docs/design/grow-from-root.md` section
3(a) - the posterior is invariant to the warm start, so "should it default"
is entirely a question of what a data-fitted init costs a sticky sampler,
never a question of correctness.

Summary: dbarts' four-move sampler (birth/death/swap/change) is not the
sampler either formal mixing-rate lower bound in the literature analyzes,
and the literature that IS relevant reports no mixing diagnostic and cannot
separate "warm start helped" from "many chains helped". The study built to
close that gap measured, rather than assumed, everything the literature
could not: it found the warm start stays overdispersed relative to the
target everywhere (M0 branch 1, not the anticipated branch 2), passes an
aggregate early-iteration benefit test and every aggregate non-inferiority
test in both size strata - and still gets KILLED in both, because a
per-cell harm check (not the aggregate) catches a plateau-posterior-mean
cost that concentrates in noise-heavy and large-n regimes and is invisible
once averaged away. `n.grow.sweeps` stays opt-in.

---

## 1. The question

VD's framing: BART's local-move MCMC sampler is empirically known to mix
slowly over tree space (Pratola arXiv 1312.1895; Hill, Linero and Murray,
Annu. Rev. Stat. Appl. 7:251-278 section 2.3.2). A warm start built from
XBART-style greedy-stochastic grow-from-root construction
(`n.grow.sweeps > 0`, landed opt-in in `grow-from-root-warm-start.md`)
replaces today's prior-drawn stump forest with a data-fitted one before the
exact Metropolis-Hastings sampler takes over. The posterior is unchanged
(warm start only, never a standalone sampler) - but if the sampler is sticky
enough that four independently-seeded chains never leave the neighborhood
of where they started, defaulting the warm start silently converts those
four chains from a convergence diagnostic (do independent starts agree?)
into four correlated replicates of one point estimate. The question is not
"does the posterior change" (it provably does not) but: **what does walking
a sticky posterior from a data-fitted init cost, in long-run sampling
behavior and in the diagnostic value of running multiple chains, relative
to today's prior-drawn stump default** - and is that cost small enough to
flip the default.

## 2. Literature, as corrected

Full detail: an untracked companion document
(gitignored working paper; this section compresses it and keeps every
`[verified: url]` tag). The corrected reading is materially different from
a first-pass reading of the same sources, and the corrections all widen the
gap between what is in print and VD's question rather than closing it.

**The two formal mixing lower bounds do not bind dbarts, and both name a
data-fitted init as a remedy, not a hazard.**

- Kim and Rockova, "On Mixing Rates for Bayesian CART", EJS 2025, 19(2)
  3041-3067 (DOI 10.1214/25-EJS2397)
  [verified: https://doi.org/10.1214/25-EJS2397]. The published version
  says **superpolynomial**, not the preprint's "exponential" - a
  misquotation the memo that started this arc made by citing the journal
  while quoting arXiv 2306.00126v1. The model is a SINGLE tree,
  one-dimensional, dyadic, with only grow-a-node / prune-two-sisters moves;
  Theorem 5.1 is the worst case over initializations, scoped to
  "initializations which need [to] grow through noise to catch" the target.
  The EJS version adds a referee-prompted paragraph absent from the
  preprint naming a dynamic-programming tree point estimate as "a good
  initialization for which faster mixing rate could be achieved" - the
  paper invoked as the stickiness premise names a data-fitted
  initialization as its own remedy.
- Ronen, Saarinen, Tan, Duncan, Yu, arXiv 2210.09352
  [verified: https://arxiv.org/abs/2210.09352]. Four simplifications, not
  two: single tree, Grow/Prune only ("thereby excluding the Change and
  Swap moves"), discrete features, fixed leaf-prior variance. Their own
  caveat: allowing Change or Swap breaks the bottleneck-set proof. Their
  recommendation list includes initializing "at an intelligent guess for
  the possible trees" - the second negative-result paper recommending a
  warm start.
- dbarts has all four moves:
  `enum class StepType { birth, death, swap, change };`
  (src/bartcore/moves.hpp:838). Neither lower bound analyzes a sampler
  dbarts runs, on a model (one tree, dyadic, two moves) dbarts fits.

**What does support the stickiness premise is empirical, narrower, and
about coverage, not point estimates.** Pratola arXiv 1312.1895
[verified: https://arxiv.org/abs/1312.1895]; Hill, Linero, Murray
[verified: https://doi.org/10.1146/annurev-statistics-031219-041110];
Carnegie, Stat. Sci. 34:90-93 (2019)
[verified: https://doi.org/10.1214/18-STS696] - whose modifications all
help coverage but where "use of TMLE is least beneficial for coverage and
results in considerably wider confidence intervals" while SATT bias stayed
"minimal, regardless of the BART formulation": coverage bought with
interval width is the first instance of a pattern that recurs through this
literature.

**The positive warm-start evidence is real, but one lab, and reports no
mixing diagnostic at all.** He and Hahn, arXiv 2002.03375
[verified: https://arxiv.org/abs/2002.03375]: XBART-initialized MCMC
"considerably improves credible interval coverage and reduces total
run-time". Their Table 4 (n=10000, nominal 95%, XBART/BART/WS-BART),
corrected in full:

| DGP | k=1 XBART/BART/WS | k=2 XBART/BART/WS |
|---|---|---|
| Linear | 0.78/0.77/0.99 | 0.50/0.83/0.98 |
| Max | 0.86/0.78/0.95 | 0.88/0.84/0.97 |
| Single Index | 0.77/0.73/0.87 | 0.81/0.83/0.91 |
| Trig+Poly | 0.90/0.74/0.96 | 0.90/0.82/0.96 |

A prior pass on this arc had misread the table's 0.50 cell as BART's
minimum (it is XBART's) and printed only 4 of 8 cells, hiding that Single
Index is the one WS cell materially below nominal (0.87, 0.91). Correctly
read: warm start is a real +8 to +22 point coverage gain, not a uniform
0.95-0.99. Three confounds, all verified in the same paper: (1) the
headline contrast is 25 warm chains x 100 iterations against ONE cold
chain of 1000+2500 - Hill-Linero-Murray's own "run more chains" prescription
confounded with initialization, and the cold arm got MORE iterations and
wall clock, not less; (2) WS intervals are uniformly 30-60% wider (the
Carnegie/TMLE pattern again); (3) zero occurrences of R-hat, ESS,
autocorrelation or trace plots anywhere in the paper. Krantsevich, He and
Hahn (XBCF), arXiv 2209.06998
[verified: https://arxiv.org/abs/2209.06998] is better controlled (a real
20000+20000 cold control that still underperforms warm-start BCF) but
shares both caveats, with the widest intervals of any arm (CATE 0.376 vs
BCF 0.262).

**The seed-independence inversion.** The memo's central structural claim
was that dbarts' per-chain independently-seeded grows reproduce He-Hahn's
design. Backwards: He-Hahn's 25 "independent" chains are 25 successive,
autocorrelated states of ONE 40-sweep XBART chain (15 discarded, 25
retained), distributed one-per-chain. stochtree does the same and
documents it (`forest_ind <- num_gfr - chain_num`,
`num_gfr >= num_chains` required). dbarts' k grow sweeps run
INDEPENDENTLY per chain, each on its own Mersenne Twister
(src/bartcore/sampler.hpp:1260-1290, `Chain::growForestFromRoot`) - a
structurally different design that, if anything, resembles Tan et al.'s
per-chain-varied greedy fit more than it resembles either published
positive design. Whether dbarts' design is more, equally, or less
dispersed than the published ones is unmeasured anywhere in print, and it
is the hinge: a warm start MORE dispersed than the posterior costs nothing
(Gelman-Rubin recommends overdispersed starts); one that concentrates the
chains costs exactly the diagnostic value in question. This is why M0
(below) is a hard premise gate.

**stochtree's default is 5, not 10, and the paper misreports its own
package.** Herren, Hahn, Murray, Carvalho, arXiv 2512.12051
[verified: https://arxiv.org/abs/2512.12051]: prose says `num_gfr = 10`
twice; the shipped default in `R/bart.R` (main and the CRAN 2026-06-30
manual, and at both the paper's v1 and v2 dates) is `num_gfr = 5`. Cite
the source, not the prose - the correction sharpens the point that the
only published default sits between dbarts' unmeasured R5 method default
of 2 and the study's grid, equally unmeasured.

**Tan et al., restated at the source's own strength.** Tan, Ronen,
Saarinen, Yu, arXiv 2406.19958v2
[verified: https://arxiv.org/abs/2406.19958]. Section 9: a fitted-XGBoost
init "can increase R-hat values and exaggerate their increasing trend...
the chains become stuck around different local maxima." This is the only
R-hat measurement of a data-fitted BART init in print, and it points the
direction VD suspects - but a prior pass overstated it as a flat claim.
The source hedges ("seems to"), states no effect size, and three
corrections narrow it further: (1) the init is NOT "essentially
deterministic" or "a greedy point" - the paper's own repository
(`yanshuotan/bart-comp-efficiency`, `runner.py:66-98`,
`consolidated_configs/initialization.yaml`) shows each of the 8 chains gets
a SEPARATELY, randomly tuned XGBoost ensemble (depth, learning rate,
subsample, colsample_bytree, `RandomizedSearchCV(random_state=seed)`
inside `run_chain_bart`) - so the raised R-hat is as consistent with
OVERDISPERSED starts as with stuck ones, and the paper does not
distinguish "stuck around different modes" from "started far apart, not
yet met"; (2) the paper's own companion Experiment 2 finds "increasing the
number of trees consistently dampens the increasing trend for R-hat", the
m=20 choice (vs dbarts' 75) was deliberately chosen "to emphasize the
initialization effect", so the m=20 experiment is not evidence at
ship-default ensemble sizes; (3) Theorem 7.2's O_P(1) mixing bound - the
theoretical result closest to rescuing either side - holds for a LAZIFIED
chain of a MODIFIED sampler with a fully marginalized MH filter, discrete
covariates, and a chained-proposal-length condition (one proposal must
traverse the whole state space); section 8 concedes the results only
address how iteration count scales with n, never per-iteration cost, and
that the modifications "may not improve the overall effectiveness of the
sampler". A one-shot initializer collects none of this.

**Tree-space stickiness is benign only for function-space functionals.**
Deshpande, flexBART, arXiv 2211.04459
[verified: https://arxiv.org/abs/2211.04459] appendix B3 and Tan et al.
section 1.3 both make the same point: BART's tree-structure non-identifiability
means failure to mix over TREE STRUCTURE need not mean failure to mix over
f(x) - "point estimates and reasonably well-calibrated uncertainty
intervals" survive despite structural stickiness, scoped explicitly to
"using only grow and prune moves". This benignness does NOT extend to
tree-structure functionals: variable-inclusion proportions, interaction
reporting, `plotTree`, anything in bartCause/treatSens that reads
structure. flexBART section 2 supplies the mechanism a more informative
proposal is plausibly counter-productive: an informative grow proposal's
prior/proposal density ratio is far from 1, deflating acceptance and
INHIBITING exploration - a statement about proposals, not initializers,
but the only published mechanism for the effect VD is asking about.

**Where the literature is silent, once every confound above is counted:**
nobody has measured R-hat, cross-chain dispersion, or mode recovery under
a STOCHASTIC grow-from-root warm start (the two papers using one report no
diagnostic; the one paper reporting R-hat used XGBoost); nobody has
separated "warm start helped" from "25 chains helped"; nobody has
published a tree-structure functional as a function of initialization;
nobody has published a sensitivity study over the grow-sweep count k
(XBART's 40/15 and stochtree's 5 are both asserted, not measured, and the
stochtree paper misreports its own value); no constructed-multimodality
warm-vs-cold benchmark exists. Two general results bound the question
without answering it: Gelman and Rubin, Statist. Sci. 7:457-472 (1992)
[verified: https://doi.org/10.1214/ss/1177011136] recommend overdispersed
starting points - cutting against a warm start only if it is
UNDER-dispersed, an empirical question about dbarts, not a fact about
data-fitted inits in general; Koehler, Lee and Vuong, arXiv 2411.09117
[verified: https://arxiv.org/abs/2411.09117] point the other way but more
weakly than first read - the result needs a kth-order spectral gap, a
Poincare inequality, and is scoped by its own framing sentence to "score
matching methods" over a continuous state space via Langevin diffusion or
Glauber dynamics; it needs samples FROM the stationary law, and
grow-from-root draws are not that.

**The honest answer, as the literature stood before this study:** nothing
in print answers VD's question, and every correction above widens the gap
rather than closing it. The two formal negative results both name a
data-fitted init as a remedy for the very bound they prove. The strongest
positive evidence shows real coverage gains but is one lab, confounds
chain-count with initialization, buys coverage partly with interval width,
and reports no mixing diagnostic. The single contrary result hedges,
reports no effect size, used a deliberately amplifying ensemble size, and
initialized from per-chain-varied overdispersed ensembles rather than a
common point. Above all, dbarts' design (independent per-chain grows)
matches neither published design (one shared sweep chain distributed to
many BART chains), so no result transfers by analogy in either direction.
Whether independently-seeded grow-from-root yields overdispersed, matched,
or under-dispersed starts relative to dbarts' own posterior was, before
this study, unmeasured anywhere - which is why the study's own premise
gate (M0, below) had to run before anything else in it.

## 3. The design that survived measurement

Full detail: an untracked companion document
(gitignored). Three scratch probes against the installed package replaced
several of the plan's asserted variances with measured ones and killed
three proposed remedies outright - the design below is what was left after
that filtering, not what was first proposed.

**What survived.** Paired warm-cold contrasts sharing a seed (data seed
`BASE_SEED + s`, sampler seed `s`, both arms of a pair share `s`, per the
`grouped-mixing.R` idiom) so every contrast is a matched-pair difference,
not two independent samples. Stratification on n (SMALL: n<=5000, LARGE:
n=50000) with INDEPENDENT verdicts, because a probe found the plateau cost
at +0.4% (n=2000) growing to +4.2% (n=20000) on the same metric - a pooled
verdict would average away the only regime where the cost appears. A hard
PREMISE GATE (M0) run first (~8 minutes) that can stop the whole study
before the expensive part runs, because seed-independence is not the same
claim as output dispersion and a 20-minute measurement can invalidate a
2.5-hour study. Power-first margins: every gated contrast needed BOTH a
scientific ceiling (the largest difference callable benign) and a
statistical floor of 4x the per-replicate SE at the planned R (~0.99 power
at zero effect) strictly inside that ceiling, with thresholds FROZEN at the
end of a Stage-0 pilot before Stage 1's confirmatory run - a metric whose
floor exceeds its ceiling leaves the gate and is reported descriptively
rather than adjusted after the fact. A per-cell KILL clause requiring BOTH
a Holm-corrected one-sided harm rejection AND a point estimate exceeding
its margin, in at least two cells, with a single flagged cell mandating a
fresh-seed re-run before the flag counts - because an earlier aggregate-only
design (the memo's) could not have caught a cost that concentrates in one
or two regimes and cancels in a pooled average.

**What was deleted by measurement, and why.**

- **R-hat left the gate entirely.** A probe found the per-replicate SD of
  (R-hat_warm - R-hat_cold) on the held-out-RMSE functional is 0.198 at
  n=2000 (12 reps) and 0.482 at n=20000 (6 reps), against absolute R-hats
  of 1.28 and 1.53 - neither arm converges on that functional at these
  chain lengths. A 0.05 margin at 4x SE needs roughly 250 replicates at
  n=2000 and roughly 1500 at n=20000; more chains does not rescue a
  scalar collapsing a p-dimensional disagreement. Replaced by C4 (see
  section 4), which measures the same concept - persistent cross-chain
  disagreement - by pooling p inclusion-proportion columns instead of one
  scalar; its probe SD/rep was 0.0011-0.0015 against a cold level of
  0.0075 and was already significant at 12 reps. R-hat survives only as
  UNGATED reported context (X2).
- **The m=75 ensemble self-averages any structural label.** At the ship
  default of 75 trees, an ensemble cannot express "the root is on x1 or
  x2" as a binary outcome - the label self-averages at a rate of
  1/sqrt(#splits) with no mixing required to produce that average, so no
  statistic recovers structural mode collapse at m=75 at any feasible
  replicate count. The structure family (S1m/S2m, below) therefore only
  runs at m=1, on its own scenario (S5), and the m=75 cells are gated on
  function-space cost plus C4 rather than on a powerless structural
  statistic.
- **The duplicated-column design proves nothing about the sampler.** The
  original structure probe used two bitwise-identical predictor columns
  (x1==x2) so the two "modes" (root on x1, root on x2) are exactly
  exchangeable by symmetry. But with x1 and x2 identical, a CHANGE move
  swapping the split variable between them has integrated-likelihood
  ratio EXACTLY 1 and prior ratio 1 - it is accepted with probability ~1.
  The two "modes" are trivially connected; there is no valley for
  anything to get stuck in, so a duplicated pair measures the estimator,
  not the sampler. Replaced by an XOR construction (S5:
  `4*XOR(x1>.5, x2>.5)`) whose root split on x1 or x2 ALONE yields zero
  fit improvement (both children have the same conditional mean), so the
  chain must pass through a likelihood-neutral intermediate state to
  switch roots, while the two modes remain exactly exchangeable with
  analytic answer 1/2. A confirming probe (m=1, n=2000, 4 chains x 300
  draws, cold init) found real separation: per-chain root-on-x1 fraction
  q = 0.057/0.440/0.100/0.000 across the 4 chains, root switches
  8/8/12/0, between-chain SD of q = 0.198 against a mixing null near
  0.05 - one chain locked onto x2 for 272 of 300 draws while the others
  moved 8-12 times. S6 (exact duplicates, `x1==x2`) is demoted to a NULL
  CONTROL, where its trivial connectivity is exactly the useful property:
  both arms must return pooled p1 within MC error of 1/2 with non-zero
  switch counts in every chain, or the estimators themselves are broken.

## 4. Results, in full

Executed against the pre-registration at commit `84d2594`. Read-only on the
package throughout: no engine edit, no default change, no baseline
re-record. Every number below is a function of the draw sequence; wall
clock enters no criterion. Sign convention throughout: **positive = warm
worse** for C1-C4, X1, S1m, S2m; B1 and B2 are benefit metrics, signed
positive = warm ahead / faster.

### 4.1 M0: the premise gate

`t=0` is the first KEPT draw of the `n.burn=0` trajectory (not the
pre-sweep initial forest - see the deviation note below). rho0 = the
between-chain SD of the t=0 test fit, divided by a reference cold
posterior SD measured separately per n (10000 burn + 10000 kept, 3 reps,
R=20 for the grid itself).

rho0 by n and warm-sweep count k (k=0 is cold):

| n | k=0 | k=1 | k=2 | k=4 | k=5 | k=8 | k=16 |
|---|---|---|---|---|---|---|---|
| n500 | 2.25 | 1.51 | 1.43 | 1.31 | 1.26 | 1.19 | 1.12 |
| n5000 | 5.40 | 2.54 | 2.25 | 1.87 | 1.75 | 1.58 | 1.38 |
| n50000 | 11.51 | 4.06 | 3.31 | 2.44 | 2.26 | 1.93 | 1.65 |

M0c, modal-root-variable share across chains (per tree, meaned):

| n | k=0 | k=1 | k=2 | k=4 | k=5 | k=8 | k=16 |
|---|---|---|---|---|---|---|---|
| n500 | 0.295 | 0.317 | 0.321 | 0.326 | 0.324 | 0.322 | 0.327 |
| n5000 | 0.297 | 0.347 | 0.350 | 0.351 | 0.353 | 0.356 | 0.363 |
| n50000 | 0.296 | 0.442 | 0.415 | 0.391 | 0.383 | 0.379 | 0.383 |

M0b, mean pairwise L1 between per-chain inclusion vectors at t=0 (reported,
not gated):

| n | k=0 | k=1 | k=2 | k=4 | k=5 | k=8 | k=16 |
|---|---|---|---|---|---|---|---|
| n500 | 0.299 | 0.284 | 0.281 | 0.262 | 0.276 | 0.259 | 0.266 |
| n5000 | 0.313 | 0.198 | 0.179 | 0.166 | 0.166 | 0.162 | 0.166 |
| n50000 | 0.307 | 0.088 | 0.091 | 0.111 | 0.110 | 0.113 | 0.120 |

Reference posterior SDs (rho0's denominator): n=500 0.9045, n=5000 0.3649,
n=50000 0.1692.

Branch rule (the plan branches on one rho0; M0 measures a k x n grid): take
the best (most dispersed) warm k per n, then branch on the worst such value
across n. Best warm k per n = 1/1/1, giving rho0 = 1.51/2.54/4.06; the
worst over n is rho0 = 1.511 with M0c = 0.442 there.

**M0 BRANCH = 1** - rho0 >= 1 everywhere: the warm start is at least as
dispersed as the target in every cell measured. Stage 1 runs as written,
with NO mandatory caveat (branch 2's caveat, and branch 3's stop, both did
not fire). This is not what the design's own probe expected (probe:
n=10000, cold 1.92 / k=2 0.81 / k=8 0.53 against a posterior SD near 0.35,
i.e. branch 2) - the fuller in-harness measurement landed one branch more
favorable than the probe that sized the study.

Supplementary, not the gated quantity: at the PRE-SWEEP initial forest the
cold arm's between-chain SD of the test fit is EXACTLY 0 at every n,
because bart2's cold path calls `sampleTreesFromPrior` only - tree
structure is drawn but every node parameter is left at 0, so all 8 chains
predict the same constant before the first sweep runs. Pre-sweep modal-root
share (seed 301): n=500 cold 0.29 / k=8 0.33; n=5000 cold 0.29 / k=8 0.36;
n=50000 cold 0.29 / k=8 0.38.

### 4.2 Frozen thresholds (end of Stage 0; not adjusted afterward)

Stage-0 pilot: R=8, warm k=4 (k* not yet selected), disjoint seed block
from Stage 1. Floor = 4 x SD_rep / sqrt(planned R); a floor above the
margin means the metric leaves the gate and is reported descriptively.

| stratum | metric | scale | SD/rep | 4xSE | margin | R | status |
|---|---|---|---|---|---|---|---|
| SMALL | C1 | rel | 0.02340 | 0.01911 | 0.03000 | 24 | GATED |
| SMALL | C2 | rel | 0.03998 | 0.03264 | 0.06000 | 24 | GATED |
| SMALL | C3 | abs | 0.00938 | 0.00766 | 0.02000 | 24 | GATED |
| SMALL | C4 | rel | 0.12750 | 0.10411 | 0.30000 | 24 | GATED |
| LARGE | C1 | rel | 0.02022 | 0.01808 | 0.03000 | 20 | GATED |
| LARGE | C2 | rel | 0.04762 | 0.04259 | 0.06000 | 20 | GATED |
| LARGE | C3 | abs | 0.01711 | 0.01530 | 0.02000 | 20 | GATED |
| LARGE | C4 | rel | 0.25795 | 0.23072 | 0.30000 | 20 | GATED |

S1m/S2m margins from their computable nulls (S5, m=1): cold-arm
root-indicator ESS 184.4 gives a perfect-mixing between-chain SD of 0.0368
(nominal at 2000 independent draws: 0.0112); perfect collapse is 0.5.
S1m margin = 0.25 x (0.5 - mixing null) = 0.1158, 4xSE = 0.1334 -> **leaves
the gate**. S2m margin = log(2) (a halving of the root-switch rate) =
0.6931, 4xSE = 0.9585 -> **leaves the gate**.

**S6 null control (plan Verification clause) FAILS.** Stage-0 pilot (R=8,
64 chain-replicates): pooled p1 cold 0.5292, warm 0.5271 (chain-level MC SE
0.0509) - the pooled fractions look fine - but zero-switch chains: cold
4/64, warm 57/64. Per the plan's own Verification clause this VOIDS the
S1m/S2m structure family; both are reported descriptively in section 4.4
below and enter no verdict.

### 4.3 Stage 1 Part A: the k grid

600 draws, 4 chains, R=8, B2 = first-hit iterate against each cell's
cold-derived `band_hi` (95th percentile of the cold arm's last-500-iterate
RMSE).

B2, mean first-hit iterate over reps x chains (censoring value 601):

| cell | band_hi | k=0 | k=1 | k=2 | k=4 | k=5 | k=8 | k=16 |
|---|---|---|---|---|---|---|---|---|
| S1-500 | 1.3841 | 114.6 | 87.0 | 97.8 | 55.2 | 65.2 | 82.2 | 46.8 |
| S1-5000 | 0.5891 | 136.3 | 172.9 | 197.2 | 221.9 | 178.9 | 151.3 | 98.0 |
| S1-50000 | 0.3727 | 129.0 | 186.1 | 194.6 | 138.1 | 109.9 | 48.5 | 5.8 |
| S2-5000 | 1.6195 | 63.9 | 192.5 | 136.7 | 91.4 | 90.3 | 84.5 | 61.5 |
| S3-5000 | 0.2920 | 64.0 | 226.7 | 160.8 | 153.6 | 169.7 | 121.6 | 143.4 |
| S4-5000 | 1.0176 | 123.1 | 33.7 | 5.6 | 1.0 | 1.0 | 1.0 | 1.0 |

Per-iterate test RMSE at t=100 (the B-vs-k curve):

| cell | k=0 | k=1 | k=2 | k=4 | k=5 | k=8 | k=16 |
|---|---|---|---|---|---|---|---|
| S1-500 | 1.2590 | 1.2217 | 1.2336 | 1.2293 | 1.2340 | 1.2589 | 1.2245 |
| S1-5000 | 0.5935 | 0.6178 | 0.5971 | 0.5888 | 0.5930 | 0.5766 | 0.5589 |
| S1-50000 | 0.4019 | 0.4164 | 0.3992 | 0.3681 | 0.3617 | 0.3343 | 0.3086 |
| S2-5000 | 1.4277 | 1.5238 | 1.4798 | 1.4799 | 1.4912 | 1.4538 | 1.4387 |
| S3-5000 | 0.2488 | 0.2739 | 0.2700 | 0.2595 | 0.2640 | 0.2609 | 0.2689 |
| S4-5000 | 1.0764 | 0.8232 | 0.7465 | 0.6891 | 0.6691 | 0.6279 | 0.6090 |

**k\* = 16 is a FALLBACK, and it is a grid artefact, not a measured
optimum.** The plan's selection rule (smallest k within 10% of the
per-cell warm-grid minimum, in EVERY cell) selected NOTHING: k=16 is the
warm-grid minimum in S1 (all three n), S2 and S4, while S3 bottoms out at
k=8 and reads 1.179x its k=16 value there - the pre-registered band is
empty, so a fallback (k minimizing the worst-cell relative excess over
that cell's own minimum) was applied at analysis time, landing on k=16 at
a worst-cell excess of 1.179. Two things this says: (1) no single k serves
every scenario in this grid; (2) k=16 is the LARGEST value in the plan's
grid and B2 is STILL falling at the edge (S1-50000 goes from 129 at cold
to 6 at k=16) - the grid does not bracket an optimum. Any future flip
proposal must re-run Part A on a grid extended past 16 before choosing a
number.

### 4.4 Stage 1 Part B: per-cell statistics, cold vs warm(k\*=16)

Raw C1-C3 (SMALL: R=24, LARGE: R=20, 8 chains):

| cell | stratum | R |
|---|---|---|
| S1-5000 | SMALL | 24 |
| S2-5000 | SMALL | 24 |
| S3-5000 | SMALL | 24 |
| S4-5000 | SMALL | 24 |
| S5-5000-m75 | SMALL | 24 |
| S1-50000 | LARGE | 20 |
| S3-50000 | LARGE | 20 |
| S7-5000 | SMALL | 24 |
| S5-5000-m1 | SMALL | 24 |
| S6-5000-m1 | SMALL | 24 |

C1, mean test RMSE over the last 500 iterates:

| cell | cold | warm |
|---|---|---|
| S1-5000 | 0.4905 | 0.4986 |
| S2-5000 | 1.4443 | 1.4456 |
| S3-5000 | 0.2560 | 0.2568 |
| S4-5000 | 0.5510 | 0.5544 |
| S5-5000-m75 | 0.3581 | 0.3488 |
| S1-50000 | 0.2485 | 0.2596 |
| S3-50000 | 0.1012 | 0.1014 |
| S7-5000 | 0.1580 | 0.1580 |
| S5-5000-m1 | 0.7078 | 0.6274 |
| S6-5000-m1 | 0.1699 | 0.1707 |

C2, test RMSE of the pooled posterior mean:

| cell | cold | warm |
|---|---|---|
| S1-5000 | 0.3052 | 0.3195 |
| S2-5000 | 0.8730 | 0.8820 |
| S3-5000 | 0.1542 | 0.1551 |
| S4-5000 | 0.3086 | 0.3424 |
| S5-5000-m75 | 0.2873 | 0.2808 |
| S1-50000 | 0.1517 | 0.1678 |
| S3-50000 | 0.0705 | 0.0702 |
| S7-5000 | 0.0930 | 0.0936 |
| S5-5000-m1 | 0.4917 | 0.4298 |
| S6-5000-m1 | 0.1039 | 0.1373 |

C3, 95% pointwise coverage of true mu (200 held-out points):

| cell | cold | warm |
|---|---|---|
| S1-5000 | 0.9715 | 0.9579 |
| S2-5000 | 0.9852 | 0.9858 |
| S3-5000 | 0.9858 | 0.9858 |
| S4-5000 | 0.9833 | 0.9662 |
| S5-5000-m75 | 0.9817 | 0.9817 |
| S1-50000 | 0.9640 | 0.9477 |
| S3-50000 | 0.9367 | 0.9390 |
| S7-5000 | 0.9910 | 0.9900 |
| S5-5000-m1 | 0.9800 | 0.9790 |
| S6-5000-m1 | 0.9256 | 0.6998 |

C4, between-chain SD of time-averaged inclusion, meaned over columns:

| cell | cold | warm |
|---|---|---|
| S1-5000 | 0.00674 | 0.00689 |
| S2-5000 | 0.00849 | 0.00915 |
| S3-5000 | 0.01000 | 0.01028 |
| S4-5000 | 0.00254 | 0.00271 |
| S5-5000-m75 | 0.01525 | 0.01393 |
| S1-50000 | 0.00616 | 0.00727 |
| S3-50000 | 0.00809 | 0.00895 |
| S7-5000 | 0.01778 | 0.01905 |
| S5-5000-m1 | 0.05634 | 0.05635 |
| S6-5000-m1 | 0.05242 | 0.06174 |

Contrasts with per-replicate SE, shown as `estimate (SE)`. C1, C2:

| cell | C1 | C2 |
|---|---|---|
| S1-5000 | 0.0168 (0.0049) | 0.0480 (0.0097) |
| S2-5000 | 0.0011 (0.0029) | 0.0102 (0.0048) |
| S3-5000 | 0.0032 (0.0024) | 0.0066 (0.0048) |
| S4-5000 | 0.0064 (0.0049) | 0.1110 (0.0127) |
| S5-5000-m75 | -0.0273 (0.0073) | -0.0255 (0.0103) |
| S1-50000 | 0.0447 (0.0059) | 0.1066 (0.0102) |
| S3-50000 | 0.0016 (0.0030) | -0.0038 (0.0053) |
| S7-5000 | -0.0001 (0.0070) | 0.0048 (0.0084) |
| S5-5000-m1 | -0.1088 (0.0374) | -0.1134 (0.0427) |
| S6-5000-m1 | 0.0046 (0.0094) | 0.3282 (0.0327) |

C3, C4:

| cell | C3 | C4 |
|---|---|---|
| S1-5000 | 0.0135 (0.0026) | 0.0301 (0.0270) |
| S2-5000 | -0.0006 (0.0015) | 0.1005 (0.0472) |
| S3-5000 | 0.0000 (0.0018) | 0.0387 (0.0298) |
| S4-5000 | 0.0171 (0.0026) | 0.0789 (0.0278) |
| S5-5000-m75 | 0.0000 (0.0005) | -0.0776 (0.0275) |
| S1-50000 | 0.0162 (0.0046) | 0.2050 (0.0509) |
| S3-50000 | -0.0022 (0.0029) | 0.1263 (0.0486) |
| S7-5000 | 0.0010 (0.0012) | 0.0871 (0.0315) |
| S5-5000-m1 | 0.0010 (0.0063) | 0.0215 (0.0414) |
| S6-5000-m1 | 0.2258 (0.0215) | 0.2632 (0.1000) |

X1 (interval width, the width confound), B2 (log2 first-hit ratio, benefit):

| cell | X1 | B2 |
|---|---|---|
| S1-5000 | -0.0025 (0.0044) | 0.3538 (0.1522) |
| S2-5000 | -0.0054 (0.0033) | -2.2700 (0.2988) |
| S3-5000 | 0.0019 (0.0030) | 0.8411 (0.2990) |
| S4-5000 | -0.0400 (0.0057) | -1.8166 (0.2366) |
| S5-5000-m75 | -0.0217 (0.0052) | -7.0646 (0.3327) |
| S1-50000 | 0.0123 (0.0056) | 0.6965 (0.0855) |
| S3-50000 | 0.0104 (0.0035) | 1.7572 (0.1063) |
| S7-5000 | -0.0009 (0.0070) | 0.7797 (0.1027) |
| S5-5000-m1 | -0.1668 (0.0538) | -2.6438 (0.5937) |
| S6-5000-m1 | -0.2836 (0.0200) | -7.7058 (0.5761) |

B1 (benefit, cold-minus-warm relative gap, positive = warm ahead) at t=1,
t=25:

| cell | B1@t1 | B1@t25 |
|---|---|---|
| S1-5000 | 0.7427 (0.0032) | 0.3736 (0.0070) |
| S2-5000 | 0.4100 (0.0065) | 0.0492 (0.0083) |
| S3-5000 | 0.5545 (0.0052) | -0.0701 (0.0105) |
| S4-5000 | 0.8087 (0.0021) | 0.6028 (0.0034) |
| S5-5000-m75 | 0.7997 (0.0089) | 0.6789 (0.0143) |
| S1-50000 | 0.8554 (0.0019) | 0.5839 (0.0044) |
| S3-50000 | 0.7784 (0.0032) | 0.0883 (0.0106) |
| S7-5000 | -0.7287 (0.0276) | -0.1810 (0.0232) |
| S5-5000-m1 | 0.2137 (0.0526) | 0.2163 (0.0528) |
| S6-5000-m1 | 0.8674 (0.0024) | 0.7856 (0.0057) |

B1 at t=100, t=250:

| cell | B1@t100 | B1@t250 |
|---|---|---|
| S1-5000 | 0.0547 (0.0072) | -0.0287 (0.0063) |
| S2-5000 | -0.0112 (0.0063) | -0.0054 (0.0060) |
| S3-5000 | -0.0469 (0.0088) | -0.0084 (0.0067) |
| S4-5000 | 0.4567 (0.0060) | 0.2571 (0.0078) |
| S5-5000-m75 | 0.4393 (0.0195) | 0.2463 (0.0191) |
| S1-50000 | 0.2322 (0.0062) | 0.0338 (0.0070) |
| S3-50000 | -0.1196 (0.0114) | -0.0804 (0.0077) |
| S7-5000 | -0.0068 (0.0124) | -0.0172 (0.0147) |
| S5-5000-m1 | 0.2235 (0.0518) | 0.2304 (0.0527) |
| S6-5000-m1 | 0.6392 (0.0088) | 0.4505 (0.0150) |

Mean first-hit iterate B, raw (censoring value 2001):

| cell | band_hi | B cold | B warm |
|---|---|---|---|
| S1-5000 | 0.5394 | 300.4 | 472.9 |
| S2-5000 | 1.6278 | 51.6 | 50.4 |
| S3-5000 | 0.2924 | 31.8 | 112.2 |
| S4-5000 | 0.6077 | 868.2 | 454.9 |
| S5-5000-m75 | 0.4941 | 396.7 | 59.4 |
| S1-50000 | 0.2704 | 780.7 | 1295.0 |
| S3-50000 | 0.1122 | 140.4 | 422.4 |
| S7-5000 | 0.1977 | 59.1 | 119.6 |
| S5-5000-m1 | 1.7511 | 635.3 | 443.5 |
| S6-5000-m1 | 0.1956 | 852.9 | 118.4 |

X2, UNGATED context: R-hat, means over replicates:

| cell | RhatR cold | RhatR warm | RhatSig cold | RhatSig warm |
|---|---|---|---|---|
| S1-5000 | 1.011 | 1.148 | 1.001 | 1.035 |
| S2-5000 | 1.026 | 1.042 | 1.006 | 1.008 |
| S3-5000 | 1.015 | 1.038 | 1.003 | 1.007 |
| S4-5000 | 1.008 | 1.245 | 1.003 | 1.116 |
| S5-5000-m75 | 1.021 | 1.315 | 1.005 | 1.034 |
| S1-50000 | 1.003 | 1.116 | 1.000 | 1.032 |
| S3-50000 | 1.007 | 1.037 | 1.000 | 1.001 |
| S7-5000 | 1.006 | 1.002 | 1.000 | 1.000 |
| S5-5000-m1 | 1.281 | 1.561 | 1.261 | 1.419 |
| S6-5000-m1 | 1.016 | 1.237 | 1.009 | 1.005 |

X2, IACT, means over replicates:

| cell | IACT-R cold | IACT-R warm | IACT-Sig cold | IACT-Sig warm |
|---|---|---|---|---|
| S1-5000 | 40.2 | 97.7 | 35.7 | 24.5 |
| S2-5000 | 23.2 | 22.7 | 12.5 | 2.9 |
| S3-5000 | 15.4 | 23.8 | 5.1 | 2.1 |
| S4-5000 | 108.0 | 96.5 | 82.2 | 20.8 |
| S5-5000-m75 | 95.9 | 38.9 | 70.6 | 2.5 |
| S1-50000 | 45.4 | 267.8 | 28.0 | 27.5 |
| S3-50000 | 14.6 | 95.5 | 5.9 | 1.0 |
| S7-5000 | 6.5 | 6.8 | 3.3 | 4.5 |
| S5-5000-m1 | 949.4 | 518.5 | 888.9 | 403.7 |
| S6-5000-m1 | 177.6 | 43.3 | 89.9 | 1.1 |

Post-hoc transparency check: the frozen floors were measured on the k=4
Stage-0 pilot; Stage 1 ran at k*=16. Realised SD/rep beside the frozen
value; NO THRESHOLD WAS CHANGED ON THE BASIS OF THIS TABLE.

| stratum | metric | frozen SD/rep | realised SD/rep | margin | OK? |
|---|---|---|---|---|---|
| SMALL | C1 | 0.02340 | 0.02350 | 0.03000 | yes |
| SMALL | C2 | 0.03998 | 0.04415 | 0.06000 | yes |
| SMALL | C3 | 0.00938 | 0.00955 | 0.02000 | yes |
| SMALL | C4 | 0.12750 | 0.16048 | 0.30000 | yes |
| LARGE | C1 | 0.02022 | 0.02105 | 0.03000 | yes |
| LARGE | C2 | 0.04762 | 0.03631 | 0.06000 | yes |
| LARGE | C3 | 0.01711 | 0.01718 | 0.02000 | yes |
| LARGE | C4 | 0.25795 | 0.22256 | 0.30000 | yes |

(4xSE at frozen and realised SD/rep both stayed under margin in every row;
omitted here for width, present in the working-paper `results.md`.)

### 4.5 Gate arithmetic - SMALL stratum

Cells: S1-5000, S2-5000, S3-5000, S4-5000, S5-5000-m75.

**(a) BENEFIT** - one aggregate one-sided test on B1@t100, alpha 0.05, warm
ahead:

| cell | B1@100 | SE | weight |
|---|---|---|---|
| S1-5000 | 0.05475 | 0.00725 | 19032.9 |
| S2-5000 | -0.01117 | 0.00631 | 25106.6 |
| S3-5000 | -0.04692 | 0.00876 | 13023.5 |
| S4-5000 | 0.45672 | 0.00602 | 27595.7 |
| S5-5000-m75 | 0.43927 | 0.01950 | 2628.8 |

aggregate = 0.15916, SE_agg = 0.00338, z = 47.05, one-sided p ~ 0 ->
**PASS**.

**(b) NON-INFERIORITY** - inverse-variance-weighted aggregate per metric,
H0: d>=delta, reject (conclude non-inferior) at z<-1.645:

| metric | d_agg | SE_agg | delta | z | status |
|---|---|---|---|---|---|
| C1 | 0.00290 | 0.00159 | 0.03000 | -17.07 | non-inferior |
| C2 | 0.01491 | 0.00296 | 0.06000 | -15.24 | non-inferior |
| C3 | 0.00091 | 0.00046 | 0.02000 | -41.61 | non-inferior |
| C4 | 0.02334 | 0.01340 | 0.30000 | -20.65 | non-inferior |

All four PASS.

**(c) KILL check** - one-sided harm tests over every gated C1-C4 x cell
contrast, Holm at family-wise 0.05, AND point estimate exceeding delta;
KILL needs both in >= 2 cells. Estimate/SE/z (only contrasts relevant to
the outcome shown; full 20-row table in the working paper):

| metric | cell | d | SE | z |
|---|---|---|---|
| C2 | S4-5000 | 0.11098 | 0.01267 | 8.76 |

| metric | cell | raw p | Holm p | >delta | KILL |
|---|---|---|---|---|---|
| C2 | S4-5000 | 9.86e-19 | 1.97e-17 | yes | **YES** |

Every other C1-C4 x cell contrast in this stratum is either not
Holm-significant or does not exceed its margin (full 20-row table in
`results.md`). Only S4-5000/C2 flags.

Fresh-seed re-run of S4-5000 (seeds 601-624, deviation D16 rule):

| metric | d | SE | z | delta | confirms? |
|---|---|---|---|---|---|
| C2 | 0.11579 | 0.01068 | 10.84 | 0.06000 | **YES** |

Re-run **CONFIRMS** - counts as the second cell, KILL clause met.

**SMALL verdict: KILL.** (a) benefit PASS; (b) non-inferiority PASS on
every gated metric; (c) no-cell-triggers-KILL FAILS.

### 4.6 Gate arithmetic - LARGE stratum

Cells: S1-50000, S3-50000.

**(a) BENEFIT:**

| cell | B1@100 | SE | weight |
|---|---|---|---|
| S1-50000 | 0.23222 | 0.00621 | 25952.9 |
| S3-50000 | -0.11960 | 0.01144 | 7637.0 |

aggregate = 0.15223, SE_agg = 0.00546, z = 27.90, one-sided p = 1.33e-171
-> **PASS**.

**(b) NON-INFERIORITY:**

| metric | d_agg | SE_agg | delta | z | status |
|---|---|---|---|---|---|
| C1 | 0.01048 | 0.00269 | 0.03000 | -7.26 | non-inferior |
| C2 | 0.01938 | 0.00467 | 0.06000 | -8.69 | non-inferior |
| C3 | 0.00290 | 0.00244 | 0.02000 | -7.02 | non-inferior |
| C4 | 0.16386 | 0.03515 | 0.30000 | -3.87 | non-inferior |

All four PASS.

**(c) KILL check** - flagged contrasts:

| metric | cell | d | SE | z |
|---|---|---|---|
| C1 | S1-50000 | 0.04473 | 0.00593 | 7.54 |
| C2 | S1-50000 | 0.10660 | 0.01021 | 10.44 |

| metric | cell | raw p | Holm p | >delta | KILL |
|---|---|---|---|---|---|
| C1 | S1-50000 | 2.38e-14 | 1.67e-13 | yes | **YES** |
| C2 | S1-50000 | 7.84e-26 | 6.28e-25 | yes | **YES** |

Both flags land on the SAME cell (S1-50000), which is why the plan's rule
still routes through a re-run rather than an immediate KILL: one CELL
flagged, even with two metrics.

Fresh-seed re-run of S1-50000 (seeds 601-620):

| metric | d | SE | z | delta | confirms? |
|---|---|---|---|---|---|
| C1 | 0.03450 | 0.00545 | 6.33 | 0.03000 | **YES** |
| C2 | 0.09873 | 0.01255 | 7.87 | 0.06000 | **YES** |

Re-run **CONFIRMS** - counts as the second cell, KILL clause met.

**LARGE verdict: KILL.** (a) benefit PASS; (b) non-inferiority PASS on
every gated metric; (c) no-cell-triggers-KILL FAILS.

### 4.7 S7: the descriptive categorical cell

Pre-registered prediction: a NULL. `grow.hpp`'s categorical scan skip
scales only the split-vs-no-split ODDS by A = prior mass on scannable
available columns (bounded in [1/p, 1]) against integrated-likelihood
ratios of exp(O(n)) - the grow.hpp arithmetic is real but the harm was
never established, so this cell tests it rather than assuming it.
Reopening trigger (pre-registered): S7's contrast landing beyond the COST
kill margin.

| metric | contrast | SE | margin | z | beyond margin? |
|---|---|---|---|---|---|
| C1 | -0.00007 | 0.00696 | 0.03000 | -0.01 | no |
| C2 | 0.00482 | 0.00842 | 0.06000 | 0.57 | no |
| C3 | 0.00104 | 0.00116 | 0.02000 | 0.89 | no |
| C4 | 0.08710 | 0.03147 | 0.30000 | 2.77 | no |

**S7 outcome: consistent with the null.** No C1-C4 contrast lands beyond
its COST kill margin. The documented-behavior disposition (D1(ii) in the
pre-registration) stands; the categorical renormalization is NOT
reopened. Also reported: B1@100 = -0.00684 (SE 0.01243), X1 = -0.00086
(SE 0.00703).

### 4.8 Structure family and the root-rotation stickiness datum

The S6 null control FAILED at Stage 0 (section 4.2), so per the plan's own
Verification clause the S1m/S2m estimators are VOID and neither enters a
verdict. Both are reported below because the failure is itself a finding
about the m=1 sampler, not about the estimators' arithmetic.

**S5-5000-m1** (R=24, 8 chains): pooled p1 (root on x1) cold 0.3707, warm
0.3008; between-chain SD of q: cold 0.3619, warm 0.3249; root switches per
chain: cold mean 123.5 (min 0, zero in 5/192), warm mean 80.8 (min 0,
zero in 53/192); S1m = -0.03700 (SE 0.03549, margin 0.1158); S2m = 1.34782
(SE 0.27328, margin 0.6931).

**S6-5000-m1** (R=24, 8 chains): pooled p1 cold 0.5143, warm 0.5241;
between-chain SD of q: cold 0.3924, warm 0.4600; root switches per chain:
cold mean 3.8 (min 0, zero in 10/192), warm mean 0.3 (min 0, **zero in
147/192**); S1m = 0.06757 (SE 0.01349); S2m = 1.22796 (SE 0.04275).

**The root-rotation stickiness datum.** The 147/192 zero-switch warm
chains at S6 (exact x1==x2 duplicates, where a change move swapping the
label has likelihood ratio exactly 1 and should switch freely) is not an
estimator bug: it is a structural fact about dbarts' change move.
`changeMove` (src/bartcore/moves.hpp:475-616) draws its target node
uniformly over ALL non-bottom nodes INCLUDING the root
(`tree.fillNotBottom(0, notBottom)` at moves.hpp:482, seeded from index 0;
uniform draw at moves.hpp:485-487) - the root is not structurally excluded
from proposal. But once proposed, the move keeps the ENTIRE descendant
skeleton fixed (every split rule below the changed node is untouched) and
only swaps the changed node's own rule, then reroutes every observation
through that unchanged skeleton (`Tree::refreshSubtree`, called at
moves.hpp:593) and rescores the FULL subtree's likelihood under the new
routing (`logLikelihoodForBranch`, moves.hpp:588 and :599, walking every
bottom descendant via `fillBottom` and vetoing with `-HUGE_VAL` if any
leaf empties, moves.hpp:79). Acceptance is `alpha = exp((belowY-belowX) +
(yLogL-xLogL) + correction)`, capped at 1 (moves.hpp:609-611). At the root
of a tree with any depth, "the subtree below the changed node" is the
entire tree: rerouting all n observations through a new top split while
every rule beneath it stays fixed makes it overwhelmingly likely some
downstream leaf empties or scores poorly, so acceptance collapses
combinatorially with depth even when swapping the ROOT's own label would
cost nothing in principle (as at S6). This pins the root variable once a
single tree grows past a stump - it is a mechanical property of the change
move's fixed-skeleton design, not a bug and not specific to grow-from-root:
the SAME mechanism would pin a cold-started tree's root too, but a cold
start reaches its first non-trivial root only after passing through
shallow states where the move CAN still act on it, while a grow-from-root
warm start starts already deep. Root mobility at m=1 is therefore mostly a
measurement of how often the tree collapses back to depth 1 (which the
grow-from-root arm essentially never does), not of general tree-space
mixing - and any future structure-metric design needs a statistic that
does not route through the change move's root sensitivity.

**Erratum (2026-08-09, tree-mixing synthesis).** The datum cannot
separate acceptance collapse from proposal rarity: the change move
proposes the root at only ~0.4 x 1/|notBottom| x 1/p per sweep
(~2e-3/sweep at the S6 shape), which alone predicts a switch count
matching the cold arm's recorded rate with no acceptance-collapse
contribution required. The fixed-skeleton acceptance mechanism above is
real as mechanism but is NOT what S6 measured; the corrected reading
and the probe that would separate the two live in
docs/design/tree-mixing-proposals.md (mode A).

The m=75 arms cannot detect structural mode collapse at all (the ensemble
self-averages the x1/x2 label at rate 1/sqrt(#splits) with no mixing
required); no m=75 cell is gated on a structural statistic; those cells
are gated on function-space cost plus C4, as pre-registered.

### 4.9 Invariance diagnostic

S1, n=5000, one cold and one warm(k=4) chain, 40000 draws, second half:
|z| on the sigma posterior-mean difference = 3.51 (threshold 3), held-out
fit gap mean|fhat_w - fhat_c| / posterior SD = 0.5693 (threshold 0.1).
**Disagrees** - the discriminator (4x length) runs.

At 160000 draws: |z| = 7.28, fit gap = 0.3904 -> **shrinks with length:
SLOW MIXING**, reported as a primary result, not a bug. Detail:

- Fit-gap ratio 0.3904/0.5693 = 0.686 (1/sqrt(4)=0.5 would be the ideal
  1/sqrt(length) rate) - the two single long chains are still converging
  toward each other, but slowly.
- **The sigma posterior-mean difference (warm - cold) FLIPS SIGN** between
  lengths: -0.00138 at 40000, +0.00188 at 160000. A systematic
  init-dependent bias cannot change sign with chain length; a
  slowly-mixing chain can. |z| grows only because coda's spectral MC SE
  shrinks with length faster than the (noise-driven) difference does - the
  known failure mode of a spectral ESS on a very sticky chain. The sign
  flip, not the growing |z|, is the informative part.
- Read against Stage 0's other finding (neither arm converges on the
  held-out RMSE functional at these lengths - absolute R-hat 1.28 at
  n=2000 and 1.53 at n=20000 for BOTH arms), this says the study compares
  two NON-CONVERGED samplers, which is exactly why the verdict rests on
  the plateau (a quantity both arms reach in the sampled window) and not
  on a convergence diagnostic (a quantity neither arm reaches).

## 5. Verdict and consequences

**KILL in both strata**, each confirmed by a mandated fresh-seed re-run.
`n.grow.sweeps` stays opt-in; the default forest-init path (prior-drawn
stump, `sampleTreesFromPrior`) is unchanged.

- The early-iteration benefit is real and passes its aggregate test in
  both strata (SMALL z=47.05, LARGE z=27.90 in `results.md`'s own gate
  arithmetic - see the Provenance note on a discrepancy with the
  pre-registration's summary prose for the SMALL number specifically).
  Every aggregate non-inferiority test on C1-C4 passes in both strata.
  Aggregation alone would have shipped GREEN.
- The per-cell KILL clause - not any aggregate - is what catches the real
  finding: **the plateau posterior-mean-RMSE cost concentrates in specific
  regimes and pooling averages it away.** SMALL: S4-5000 (Friedman p=50,
  45 of 50 columns pure noise) carries C2 = +11.10% against a 6% margin
  (Holm p=1.97e-17), re-run +11.58% (z=10.84). LARGE: S1-50000 carries
  C1 = +4.47% (margin 3%) and C2 = +10.66% (margin 6%), both
  Holm-significant; re-run C1 +3.45% (z=6.33), C2 +9.87% (z=7.87). Both
  re-runs CONFIRM.
- Warm starts stay TARGET-OVERDISPERSED everywhere measured (M0 branch 1,
  rho0 >= 1 in every grid cell) - the diagnostic-value concern that
  motivated the premise gate did not materialize; the cost that did
  materialize is a plateau accuracy cost in noise-heavy and large-n
  regimes, not a lost-diagnostic-value cost.
- **The categorical scan skip (D1 in the pre-registration) is CLOSED as
  documented behavior.** S7 matched its pre-registered null; the
  renormalization stays undone; the old prerequisite "(1) categorical scan
  support" is MOOT for this flip specifically (it stands only as ordinary
  opt-in-quality backlog, unscheduled). **Landing note (2026-08-12):** real
  categorical split support has since landed
  (docs/design/grow-from-root.md section 8; 995002ef..7f82f560) - a
  categorical column is now scanned and split at a prior-mass-commensurate
  weight rather than counted-but-skipped. The renormalization named above
  stays closed, superseded rather than revisited. This carries no verdict
  on the `n.grow.sweeps` default: the default question stays closed.
- **The change-move root-rotation limitation** (section 4.8) is recorded
  as a structural datum about dbarts' sampler, independent of this flip's
  outcome, for any future structure-metric design.
- **k\*=16 is a grid artefact**, not a validated recommendation for any
  future opt-in guidance: B2 was still falling at the largest grid value
  tested.
- Family consistency (four bart2 families refuse `n.grow.sweeps` via
  `checkFamilyUnsupportedArgs`, R/bart.R:618-646; `rbart_vi` has no
  `n.grow.sweeps` formal at all) is now an ordinary opt-in-surface
  question - extend the surface or document the gap - independent of this
  verdict, low priority, unscheduled.

## 6. Reopen clause

Reopening the DEFAULT question requires BOTH of the following, not either
alone:

1. **A changed design that raises the plateau cost's regime specificity
   without giving up dispersion or the early benefit** - candidates named
   in the pre-registration and not built or measured here: warm-starting
   only a SUBSET of chains (so at least one chain per fit still runs a
   from-prior, unconfounded check), or post-init TEMPERING (an annealed
   or partially-randomized handoff from the grow-from-root state into the
   exact sampler, softening the "already deep" start that section 4.8
   identifies as the mechanical reason the root gets pinned).
2. **A re-run of this battery on a k grid extended past 16**, since B2 was
   still falling at the grid's edge and k*=16 was never a measured
   optimum - any new default proposal needs Part A's B-vs-k curve to
   actually turn over before a k is chosen, in every scenario or with an
   honest per-scenario fallback rule agreed in advance (the plan's own
   10%-of-minimum rule returned empty on this grid and needs revision, not
   just a wider range, if it is to select anything at all).

Separately, and narrower: the categorical scan skip's own reopening
trigger (S7's contrast landing beyond the COST kill margin) is
pre-registered and did NOT fire here. It stays available for a future
study with a heavier categorical load (more scannable columns, a design
where the A-factor bites harder) without needing to reopen the
grow-from-root default question at all.

## 7. Provenance

```
tip           84d2594e24faa900f11b36e8c8c6bdfbaea44948
package       dbarts 1.0.0
R             R version 4.6.1 (2026-06-24), arm64 darwin, 10 cores
install       R CMD INSTALL --preclean -l <isolated library> .
n.threads     8 (pinned; ungated - wall clock enters no criterion)
BASE_SEED     20260808
compute       6933 s measured (1.93 h), summed from per-cell checkpoint
              timings; the pipeline was interrupted once and resumed, so
              this is a floor, not a ceiling, on true compute
```

**Process note.** An intermediate analysis pass over partial checkpoints
printed a GREEN verdict before the mandated fresh-seed re-runs for the
single flagged cell in each stratum had been run. The verdict logic was
corrected to a PENDING state pending those re-runs, and the premature
GREEN was never a valid result - both strata's final, re-run-confirmed
verdict is KILL (section 5). Recorded here because it is exactly the
failure mode the per-cell KILL clause with mandatory confirmation exists
to prevent, and because a partial run of this study, read incautiously,
would have shipped the wrong default.

**17 deviations were logged against the pre-registration** (full log in
the working-paper `results.md`); the consequential ones:

- `t=0` was defined as the first KEPT draw, not the pre-sweep forest (the
  pre-sweep cold state is degenerate by construction - section 4.1).
- M0's branch rule was extended from a single rho0 to a k x n grid rule
  (worst-over-n of the best-per-n warm k), since the plan branches on one
  number but M0 measures a grid.
- k\* selection needed a fallback rule invented at analysis time because
  the pre-registered 10%-of-minimum band was empty across cells (section
  4.3) - this is the single largest deviation and the reason k\*=16 is
  flagged as a grid artefact rather than a result.
- Frozen C1-C4 thresholds were measured on the k=4 Stage-0 pilot while
  Stage 1 ran at the later-selected k\*=16; the realised Stage-1 SD/rep is
  reported beside the frozen value as a transparency check (section 4.4)
  and no threshold was changed on the basis of it.
- The S6 null control failed its non-zero-switch clause, voiding the
  S1m/S2m structure family per the plan's own stated consequence (section
  4.2, 4.8) - not an estimator bug, but the root-rotation datum.
- The single-flagged-cell re-run's confirmation rule (any originally
  flagged metric again showing one-sided harm p<0.05 AND exceeding its
  margin, Holm dropped for a single cell) was fixed before any re-run data
  existed, in the direction more permissive to KILL.

**Numeric discrepancy, reported and not resolved.** The pre-registration's
committed Results section (`docs/plans/archive/grow-from-root-default-study.md`)
states the SMALL-stratum aggregate benefit as "+0.199, z=25.5"; the
gate-arithmetic table in the working paper `results.md` (section 4.5
above) computes aggregate=0.15916, z=47.05 from the same five cells' B1@100
estimates and SEs. The LARGE-stratum numbers match exactly between both
sources (+0.152/z=27.9 vs 0.15223/z=27.90). Neither this document nor its
author has reconciled the SMALL discrepancy; it may originate in the
premature-GREEN intermediate pass noted above computing from a partial or
differently-weighted set of cells. Flagged for the record, not corrected
here - `results.md`'s figure is what is transcribed into this document's
own tables per the instruction that every transcribed number match the
full working-paper record exactly.

**Where the work lives.** The harness (`common.R`, `cells.R`, `m0.R`,
`stage0.R`, `freeze.R`, `partA.R`, `partB.R`, `analyze.R`, plus
`smoke.R`/`diag-root.R`/`probe-shapes.R`/`verify-threads.R`/`rerun.R`/
`run-all.sh`), its 98 `.rds` checkpoints (~86 MB), and the working papers
`results.md`, `lit-conclusions.md` and `synthesis.md` all live under
an untracked session directory in this repository, which is
entirely gitignored (`.gitignore:27`). None of it is preserved in git
history at any commit. This design document is the only durable copy of
the study's data; see section 8 for reconstructing the harness if a future
study needs to re-run this battery.

## 8. Re-running

The harness is NOT copied into `benchmarks/R/` (see the task note this
document was written under, and the reasoning below), so re-running this
battery means reconstructing it from the pre-registration
(`docs/plans/archive/grow-from-root-default-study.md`) and this document, not
recovering a script. The harness as built was an 8-file, checkpoint-driven
pipeline (`common.R` sourced by every other file via a HARDCODED absolute
path; `m0.R` -> `stage0.R` -> `freeze.R` -> `partA.R` -> `partB.R` ->
`analyze.R`, each stage reading the previous stage's `.rds` checkpoint by
name) plus a manual re-run step (`rerun.R`) invoked by hand once the
per-cell KILL check named a flagged cell - not a design that survives
being dropped into a shared, path-independent `benchmarks/R/` script
without a rewrite. Reconstructing it means, in order:

1. Re-read the pre-registration in full (`docs/plans/
   grow-from-root-default-study.md`) for the metric definitions (M0a-c,
   B1, B2, C1-C4, S1m/S2m, X1/X2), the scenario constructors (S1-S7,
   section "Scenarios and run plan"), and the seed convention
   (`grouped-mixing.R`'s idiom: data `set.seed(BASE_SEED + s)`, sampler
   `seed = s`, shared `s` per paired arm).
2. Re-run Stage 0 (the calibration pilot, R=8 warm k=4, plus the M0 grid
   at R=20 and the invariance diagnostic) to get FRESH per-replicate SDs;
   do not reuse this document's frozen thresholds verbatim unless
   deliberately reproducing this exact study - a design change (section 6)
   changes what Stage 0 should measure.
3. FREEZE thresholds against the fresh Stage-0 SDs before looking at any
   Stage-1 contrast, following the floor-above-ceiling rule (section 3)
   exactly - this ordering is the load-bearing part of the design, not an
   implementation detail.
4. Re-run Part A's k grid EXTENDED PAST 16 (section 6) before selecting a
   k\* for any new proposal; the pre-registered 10%-of-minimum selection
   rule returned empty on the original grid and needs either a wider grid
   or a revised rule, decided before running, not after seeing results.
5. Run Part B (cold vs warm(k\*)) at the frozen thresholds, apply the
   per-cell KILL clause with Holm correction and MANDATORY fresh-seed
   re-run of any single flagged cell before treating the flag as
   confirmed - the premature-GREEN process note in section 7 is the
   concrete failure this step order prevents.
6. If reconstructing the structure family (S1m/S2m), account for the
   root-rotation stickiness datum (section 4.8) in the metric design
   itself - a statistic that routes through the change move's root
   sensitivity will not discriminate deep-tree mixing from
   shallow-collapse frequency, independent of the initializer under test.
