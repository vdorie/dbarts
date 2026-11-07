# Review 4, slice S-C: the prior predictive under the sampler's prior

Status: built, in review (branch review4-prior-predictive).

rng: neutral for `run()` and every entry that existed before: every complete run, and every state the
package writes, draws bit for bit as before; the equivalence baselines and the four snapshot files are
unchanged ([Gates](#gates)). One new engine entry, `Chain::sampleHyperparametersFromPrior`, draws from the
chains' own generators. `samplePriorPredictive`'s output changes wherever a hyperprior, DART or a
non-gaussian family is in force. The private sampler runs one chain: chain 1's stream is unchanged, so
"ev" draws as before wherever nothing new is drawn, but a multi-chain sampler at `seed = NULL` takes
fewer seeds from R's stream at creation, so its "ppd" noise and R's stream after the call move; under
a pinned seed, or on a one-chain sampler, a gaussian sampler with k fixed (the default) draws as
before, bit for bit.

Findings fixed: state-01 and state-V-02 of the fourth whole-repo review, as amended by its blind critique
(amendments 1, 2 and 4 to 10 binding as written, 11 as accepted, 12 decided: no floor on k, 13 rebase
notes), with amendment 3 replaced by dec-B608. Tests added for mutation-02 and mutation-V-02. Rulings
this rests on: dec-B361 (samplePriorPredictive draws under the prior the sampler runs under), dec-B599
(nbinom, ordinal and multinomial built), dec-B594 (predict's ppd shape without weights), dec-B608 (under
dart(), whole-forest rejection with a fallback), dec-B605 (the later split rule this draw carries over to).

## Fixes

### state-V-02: the hyperparameters are drawn, not held at their start

Before: each sample redrew the trees and the leaves at the k, df, shape, cutpoints and DART weights the
private sampler started at, so on a default probit, logistic or nbinom sampler (k ~ chi(1.5, 2)) the
prior on the latent was a normal at one k rather than its t mixture, and DART's weights were uniform.

Now: each sample first calls `Chain::sampleHyperparametersFromPrior` through
`bartcore_sampleHyperparametersFromPrior`, then the existing tree and leaf entries:

- k for each forest whose k carries a chi hyperprior, `ChiKHyperprior::draw` at no leaves (shape
  df / 2, rate 1 / (2 scale^2)). A forest with neither a k prior nor DART draws nothing and consumes
  no variate.
- Under DART, `DartPrior::drawFromPrior`: alpha from the grid's prior weights alone (a new
  `gridLogPrior_`, filled beside `gridConstant_`, whose expression and order are unchanged), then the
  weights from Dirichlet(alpha / p) over the allowed columns, every other column 0;
  `numUpdatesSkipped_` is left alone. `DartPrior::update` is not touched. How the weights and the
  trees are drawn together is dec-B608's, below.
- Refusals, before any draw: an amplitude-combined chain throws `std::invalid_argument` ("the prior
  draw of a model of several forests is not built"), R refusing first as before; a k prior with an
  infinite scale throws "a k prior with an infinite scale has no prior draws", the backstop behind R's
  own refusal: "samplePriorPredictive does not support an improper leaf prior (invchi() with scale 0, or
  chi() with an infinite scale): it has no prior draws".

The entry returns chain 0's values in force, drawn or fixed, as a list (forests.drawn, attempts, k,
split.probabilities, residual.df, shape, thresholds), NULL where a family has none.

Tests: `testHyperparametersFromPrior` (test_moves.cpp) - 20,000 prior k draws under chi(1.5, 2) against
`2 sqrt(qchisq(p, 1.5))` at five quantiles (worst z 1.82), k installed and nothing else drawn; a fixed k
draws nothing (generator streams compared) and reports k; the infinite-scale and amplitude refusals
leave the generator unmoved. `testDartDrawFromPrior`: at alpha 2 and p 5, 20,000 draws, s_1 against
Beta(alpha / p, alpha (1 - 1 / p)) by KS; alpha drawn, its frequencies over ten grid bins against the
Beta(a, b) kernel by chi-squared and E[s_1^2] against its grid average; under a mask the excluded
columns exactly 0, the rest Dirichlet over their count; the delay count unchanged. test_facade.cpp's
row for the new virtual. test-prior-predictive.R - on a default probit sampler (chi(1.5, 2), 20 trees),
(qnorm(ev[, 1]) - prior.mean) 2 sqrt(1.5) / k.scale against t(1.5) by KS, and under
normal(sd = invchi(3, 0.5)) (k ~ chi(3, k.scale / 0.5)) against t(3); both truncated at a latent of
+-7, where pnorm resolves the latent, over at least 2,000 values. A fixed-k gaussian sampler's "ev" is
identical to a hand loop of `sampleTreesFromPrior`, `sampleLeafParametersFromPrior` and `predict` under
one seed. Both improper spellings are refused in R's words.

### dec-B608: under dart(), the weights and the trees are one draw

The sampler never holds a tree with an empty leaf, so the prior it runs on the split weights and the
trees together is the stated Dirichlet restricted to forests in which every tree is held: the weights'
prior is the stated one times P(every tree held | s), which depends on the data. A per-tree retry at
fixed weights leaves that factor out.

`Chain::drawForestsWithSplitWeights` draws alpha and the weights (`DartPrior::drawFromPrior`), then each
tree once from the tree prior given them (`Chain::drawTreeFromPriorOnce`), and redraws the whole, weights
included, at the first tree refused, at most a given number of times. `drawTreeFromPriorOnce` is the
predicate `sampleTreesFromPrior` already rejected on, factored out of it with its RNG order unchanged:
the tree drawn by `growSubtreeFromPrior` from `CGMTreePrior`, refused for an empty leaf on membership or,
under the monotone "joint" prior, for leaves outside the cone. A later change to the tree prior's split
rule (dec-B605) reaches this draw through the same code. Each tree accepted inside an attempt is made a
chain state at once (`resetPriorTreeParameters`), and a refused one goes back to the bare root, so an
allocation failure part way finds every tree consistent and `settleAfterAllocationFailure` settles the
one in flight. On exhaustion the trees are drawn by `sampleTreesFromPrior`'s per-tree rejection under
the last weights, so the chain is a state the sampler holds, and the entry reports forests.drawn FALSE.

R (`drawPriorForests`) pools a budget of 50 whole draws per requested sample over the call; an
unexported, undocumented option, `dbarts.priorForestAttempts`, sets the call's total instead, so a test
reaches the fallback at any seed. While it
lasts each sample is the accepted joint draw, the entry having drawn the trees, so R skips
`sampleTreesFromPrior`. When it runs out, every sample of the call comes instead from the private
sampler's run with every row inactive (`setActiveRows` to 0), whose stationary distribution is that
prior: the control's `n.burn` sweeps, then one sample every `n.thin` sweeps, k, the family's
hyperparameters and the leaves drawn afresh for each through the entry at a budget of 0 (with no row
active the tree and weight moves do not read them, so this is a Gibbs block of the same chain). The call
then warns: "samplePriorPredictive: under dart() fewer than one in 50 independent draws of the split
weights and trees was one the sampler holds, so the draws come from the sampler run with every row
inactive and successive draws are dependent; see ?samplePriorPredictive". The help says which was used
and when. Acceptance measured on this build: 0.01 at 60 rows, 8 columns and 200 trees; 0.035 at 200 rows;
0.36 at 1000 rows.

Tests: test-prior-predictive.R, the small case - x1 binary, x2 on {0, 1, 2} (two equal cuts), four of the
six cells occupied, 20 trees, alpha fixed at 2 so s_1 is uniform under the stated prior. There
P(every tree held | s_1) has a closed form, a recursion over the node's open columns and cuts with
growth 0.95 / (1 + depth)^2, so the exact density of s_1 under whole-draw rejection is
P(held | s_1)^20 on a grid of 1,000 points. 2,000 entry draws, every one accepted, against its CDF by KS;
the same draws against the stated uniform, p below 1e-6 (the per-tree retry's distribution). Under one
seeded control, dart()'s "ev" differs from the same sampler's at uniform split.probs. At 60 rows, 8
columns and 200 trees with the budget forced to 0, the call warns and returns finite draws of the right
shape at n.samples 1, 2 and 5; on one column of two values, where every tree is held, a budget of 1
serves one sample silently and runs out on the second of two, the last, which warns and returns both.
Whether the fallback reaches its prior waits on the maintainer's calls on the budget and on the
fallback's burn-in and thinning ([Measured for the open calls](#measured-for-the-open-calls)).
tests/cpp: an accepted
joint draw leaves every tree held, totals and leaves zero, the reported weights the chain's own, and a
sweep runs from it; a budget of one at 60 rows, 8 columns and 200 trees runs out within 20 calls and
leaves the same; test_alloc.cpp fails each allocation of the joint draw at budgets of 1000 and 1 on the
constant, linear, gp and monotone leaves.

### state-01: each family's observation layer

Before: every family but the binary ones took the gaussian arm, so nbinom "ppd" held negative and
fractional counts, ordinal and multinomial categories came out as real numbers and Student-t noise was
gaussian.

Now (R/dbarts.R `samplePriorPredictive`, through the posterior path's helpers):

- Student-t: "ppd" adds sigma_i times t noise at the family's fixed df, or at the df the entry drew
  from `ResidualDfPrior`'s grid weights (installed through `restoreResidualDf`). The family is read off
  the entry's residual.df, since a Student-t sampler reports family "gaussian".
- nbinom: "ev" is exp(predict(x.test, offset.test)), which carries the center and the offset; "ppd" is
  `negbinPpd(ev, r_i)` at the fixed shape or the one the entry drew from `NBShapePrior`'s grid weights
  (`restoreShape`).
- ordinal: the entry draws the cutpoints, gamma_1 = 0 and each gap exp(N(0, 1.5^2)) at
  `OrdinalResponse::priorLogGapMean_` and `priorLogGapSd_` (now public static), installed through
  `restoreOrdinalThresholds`; one set per sample, shared by every row. "ev" is
  `ordinalCategoryProbabilities(predict, thresholds)`, n.samples x rows x K; "ppd" integer codes 1..K
  (`multinomialPpdFromProbs`). A seeded ordinal "ev" does not touch R's stream.
- multinomial: "ev" is predict's [, , 1], n.samples x rows x K; "ppd" integer codes on a sampler whose
  rows each hold one trial, a one-trial count vector (`multinomialCountPpdFromProbs`) otherwise; a
  sampler carrying a category offset refuses a NULL offset.test: "'offset.test' must be given on a
  sampler that carries a category offset: one row per row of 'x.test' and one column per category (an
  all-zero matrix for the offset-free surface)".
- Hazard, aft and the variance forest are as before.

Tests (test-prior-predictive.R, each under a seeded control with set.seed before the "ev" and "ppd"
calls, so ppd given ev is the observation layer alone; every KS on at least 2,000 values):

- Student-t: student(df = 3, sigma = fixed(1)), (ppd - ev) / sigma against t(3) over 30,000 values;
  student(sigma = fixed(1)), the same against the nine-point grid mixture, its gamma(2, 0.1) weights
  written out; the entry's df over 3,000 calls against those weights by chi-squared; student(df = 3)
  with sigma drawn, (ppd - ev) at one row against sigma t(3), its CDF averaged over 2,000 quantiles of
  sigma's scaled inverse chi-squared prior.
- nbinom: (log ev - prior.mean) 2 sqrt(1.5) / k.scale against t(1.5); nbinom(shape = 2) at k = 2, a
  randomized PIT of ppd against pnbinom(, mu = ev, size = 2); nbinom() at k = 2, the same against the
  13-point grid mixture, and the entry's shape over 3,000 calls against its weights by chi-squared.
- ordinal (K = 4, k = 2): one x row at offsets -16 to 2400 by 8, so each cutpoint reads back at a row
  where its cumulative probability lies in (1e-6, 1 - 1e-6). f = -qnorm(ev[, , 1]) less the offset, at
  each sample's row nearest 1/2, against N(prior.mean, (k.scale / 2)^2); gamma_s = qnorm(cumulative ev)
  + f + offset equal across rows to 1e-8, gamma_1 = 0 to 1e-8, the 4,000 log gaps against N(0, 1.5^2);
  ppd frequencies at the zero offset against the summed ev by chi-squared; a seeded "ev" identical
  across calls with `.Random.seed` unmoved.
- multinomial (K = 3, k = 2): rows of ev sum to 1 to 1e-12; log(ev_2 / ev_1) against
  N(0, 2 (k.scale / 2)^2); an offset.test matrix moves each log ratio by o_k - o_1 to 1e-10; ppd
  frequencies against the summed ev; a count-row sampler's ppd one-trial vectors; the category-offset
  refusal, and an all-zero offset.test served.

### mutation-02 and mutation-V-02: the nbinom count kernel under a mask

test-setters-draw-nothing.R, after the `recomputed` loop: a masked nbinom sampler handed a new response
(`setResponse(updateScale = FALSE)`, which holds the shift) draws what a twin handed it before the mask
came, both from the donor's stored state, the rows out shifted by 40 so a kernel that reads them moves
the draws (mutation-02); a mask installed then cleared draws what a never-masked twin draws
(mutation-V-02); and the mask moves the draws, the control. The oracle holds the shift because creation
derives it over every row and setResponse under updateScale over the rows in. Tests only.

## Gates

Run on arm64 macOS, at most 2 cores. Fix round 1 reran, after the rebase onto 61f610bae and on
`--preclean` installs of both libraries: tests/cpp plain (559 ok), the full suite, the trio, the
snapshot files, the extra seeds, the lint set and verify-anchors; the numbers below are that round's.
The sanitizer runs and R CMD check are the build round's.

- tests/cpp, shipped config: plain 559 ok, all passed; under `-O2 -g -fsanitize=address,undefined`
  (`detect_container_overflow=0`) 559 ok, 0 sanitizer reports. The new suites: the DART prior draw's s_1
  KS p 1.00 at a fixed alpha and 0.84 under a mask, alpha's grid chi-squared p 0.849, E[s_1^2] 0.1081
  against 0.1086; k's quantiles at worst z 1.82; the alloc harness's DART arms at budgets 1000 and 1 over
  every leaf type, 14 to 17 allocations each, all failing cleanly.
- Full tinytest suite, shipped, one process: 42053 results in 303 files, 0 failed.
  test-prior-predictive.R 81 results, test-setters-draw-nothing.R 121.
- Equivalence trio on the reference build, `--bitwise` (`--strict-coverage` on equivalence.R): 55, 19
  and 11 scenarios identical against equivalence-38da480c.rds, bcf-equivalence-38da480c.rds and
  multinomial-equivalence-80b1c8d4.rds, 0 skipped.
- The four snapshot files on the reference build: 7, 3, 14 and 3 results, 0 failed.
- Each KS and chi-squared test of test-prior-predictive.R at seeds 2611, 7193, 40417 and 91 besides its
  own: 16 tests at 4 seeds, every p > 1e-3, the smallest 0.0125 (the nbinom k test). The one test that
  must reject (the draws against the uniform, not the oracle, under dart()) rejects at every seed,
  p < 1e-133.
- Fix round 1's draws: the private sampler's one chain leaves "ev" identical to the four-chain code at
  every seed tried, gaussian and probit, `seed = NULL` and pinned; "ppd" and R's stream after the
  call are identical under a pinned seed or on a one-chain sampler and move on a multi-chain sampler
  at `seed = NULL`.
- dec-B608's small-case oracle: the accepted s_1 against the exact whole-draw rejection density, KS
  p 0.19; against the uniform, p < 2e-16.
- R-loaded ASan (`-fsanitize=address,undefined` in CFLAGS, CXXFLAGS and CXX20FLAGS, the bridge compile
  line checked, run through bin/exec/R with the runtime preloaded, mode shipped):
  test-prior-predictive.R 67, test-setters-draw-nothing.R 121, test-calibration-prior-draws.R 29,
  test-dart-mixed-columns.R 13, test-generics-posteriorPredictiveDistribution.R 204,
  test-prior-init-empty-leaves.R 5; 0 failed, 0 sanitizer reports.
- R CMD check --as-cran on the built tarball: 2 NOTEs. The incoming-feasibility NOTE is not this
  slice's (an expired certificate at two existing URLs, the Date field). The other, `rt` used without
  an import, is fixed since; codetools then reports nothing. Tests OK in 217 s.
- lintr 0 lints; air format --check clean; check-rc-codoc OK (53 methods); check-win-drift OK;
  check-doc-freshness OK.
- mutation-battery.R verify-anchors clean, m44 re-anchored on the refactored monotone draw; `run m44`
  killed by its monotone killer.

Mutants, each planted alone in a copy of the tree, all killed:

| mutant | gate | failures |
| --- | --- | --- |
| k not drawn | tests/cpp; test-prior-predictive.R | 1; 3 of 67 |
| k's rate 1 / scale^2, not 1 / (2 scale^2) | tests/cpp | 1 |
| k's shape n / 2 + nu, not (n + nu) / 2 | tests/cpp | 1 |
| DartPrior::drawFromPrior skipped | tests/cpp | 4 |
| ordinal gaps at sd 1 | test-prior-predictive.R | 1 of 67 |
| ordinal gamma_1 drawn, not 0 | test-prior-predictive.R | 1 of 67 |
| DART weights drawn once, only the trees retried | test-prior-predictive.R | 3 of 67 |
| nbinom draw at prob = ev (R) | test-prior-predictive.R | 3 |
| Student-t at a grid df when df is fixed (R) | test-prior-predictive.R | 1 |
| nbinom mean without exp (R) | test-prior-predictive.R | 3 |
| mutation-02: nbinom kernel reads masked rows | test-setters-draw-nothing.R | 1 of 121 |
| mutation-V-02: cleared mask left stale | test-setters-draw-nothing.R | 1 of 121 |
| Student-t df left at its start (C++) | test-prior-predictive.R | 1 of 81 |
| nbinom shape left at its start (C++) | test-prior-predictive.R | 1 of 81 |
| Student-t sigma of sample 1 for all (R) | test-prior-predictive.R | 1 of 81 |
| finding 1 restored: a NULL sample deletes its element (R) | test-prior-predictive.R | file errors |
| the forced budget ignored (R) | test-prior-predictive.R | 1 of 81 |

## Measured for the open calls

For the maintainer's calls on the fallback's budget (call 1) and its burn-in and thinning (call 2),
unchanged in the code. Eight uniform columns, a gaussian response, 200 trees at 60 rows (the collapse
case) and 500 at 200 rows (the larger), one chain, 1,000 samples per row of the table, the fallback
forced by a budget of 0. max s is the largest split weight; the contrast is f(x_a) - f(x_b) at rows
differing in column 1 alone (0.1 against 0.9, the rest 0.5). The fallback's figures come from 1,000
dependent draws, so at a lag-1 autocorrelation near 0.9 its contrast sd is itself noisy.

| case | draw | n.burn | n.thin | P(max s > 0.9) | contrast sd | lag-1 acf of max s | s per 1,000 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| collapse | exact | - | - | 0.951 | 0.487 | 0.04 | 2.9 |
| collapse | fallback | 200 | 1 | 0.003 | 0.350 | 0.99 | 0.4 |
| collapse | fallback | 200 | 10 | 0.729 | 0.234 | 0.99 | 2.3 |
| collapse | fallback | 200 | 50 | 0.965 | 0.065 | 0.68 | 8.6 |
| collapse | fallback | 1000 | 1 | 0.000 | 0.578 | 0.86 | 0.5 |
| collapse | fallback | 1000 | 10 | 0.968 | 0.049 | 0.90 | 2.8 |
| collapse | fallback | 1000 | 50 | 0.924 | 0.565 | 0.83 | 12.7 |
| collapse | fallback | 5000 | 1 | 1.000 | 0.102 | 0.88 | 0.9 |
| collapse | fallback | 5000 | 10 | 0.902 | 0.092 | 0.97 | 7.3 |
| collapse | fallback | 5000 | 50 | 0.928 | 0.094 | 0.81 | 36.5 |
| larger | exact | - | - | 0.891 | 0.321 | 0.00 | 140.9 |
| larger | fallback | 200 | 1 | 0.000 | 0.735 | 0.94 | 0.9 |
| larger | fallback | 200 | 10 | 0.244 | 0.300 | 0.96 | 6.5 |
| larger | fallback | 200 | 50 | 0.139 | 0.470 | 0.94 | 32.0 |
| larger | fallback | 1000 | 1 | 0.000 | 0.495 | 0.91 | 1.2 |
| larger | fallback | 1000 | 10 | 0.062 | 0.749 | 0.98 | 9.8 |
| larger | fallback | 1000 | 50 | 0.155 | 0.526 | 0.90 | 45.6 |
| larger | fallback | 5000 | 1 | 0.000 | 0.104 | 0.96 | 3.2 |
| larger | fallback | 5000 | 10 | 0.184 | 0.216 | 0.94 | 29.0 |
| larger | fallback | 5000 | 50 | 0.063 | 0.436 | 0.93 | 147.0 |

No setting in the grid brings the fallback's lag-1 autocorrelation below 0.68. Its contrast sd scatters
from 0.05 to 0.75, against the exact 0.49 and 0.32, as slow mixing predicts. At the larger case
P(max s > 0.9) stays at 0.25 or below, against the exact 0.89.

The exact draw's cost, up to 50 accepted draws or 150 s per case, on other data seeds:

| rows | trees | acceptance | s per 1,000 samples |
| --- | --- | --- | --- |
| 60 | 200 | 1.5e-2 | 2 |
| 200 | 500 | 4.9e-4 | 84 (141 in the table above) |
| 60 | 500 | 3.1e-4 | 59 |
| 60 | 1000 | 3.6e-6 | 5,772 |
| 30 | 500 | 1.4e-6 | 8,826 |

## Calls made

1. The fallback switches the whole call, not the sample, once a budget of 50 whole draws per requested
   sample, pooled over the call, runs out: below an acceptance of about 2 percent the call is served by
   the run. Switching per sample would mix independent and dependent draws in one result, and a
   per-sample cap would spend the budget unevenly. User-visible (the threshold); stated for the
   maintainer with the recommendation as built. Alternatives: a fixed total cap; an argument. The
   review judged the threshold wrong (the exact draw is cheap down to about 1e-4); open with the
   maintainer, unchanged in fix round 1.
2. The fallback's burn-in and spacing are the control's `n.burn` and `n.thin`. User-visible; stated for
   the maintainer with the recommendation as built. Alternatives: fixed counts; arguments. The review
   judged these wrong (the run at the prior mixes slowly at the control's defaults); open with the
   maintainer, unchanged in fix round 1.
3. The fallback warns once per call, beside the help's account. User-visible; stated for the maintainer
   with the recommendation as built. Alternatives: an attribute on the result; the help alone.
4. In the fallback k, the family's hyperparameters and the leaves are drawn afresh for each sample by
   the entry and the leaf entry, rather than read off the run: with no row active the tree and weight
   moves read none of them, so the draw is a Gibbs block of the run's chain and each sample's k, df,
   shape and cutpoints are independent of the others. Alternative: take the run's own values, which
   mix as slowly as the run.
5. Under DART the entry keeps the accepted trees and R skips `sampleTreesFromPrior`. Alternative: the
   entry returns the weights only and R redraws the trees, which is the per-tree retry dec-B608
   replaces.
6. On exhaustion inside the engine the trees are drawn by the per-tree rejection under the last
   weights, leaving a state the sampler holds, and forests.drawn is FALSE. Alternatives: restore the
   trees from before the call (a copy per call); leave the last refused draw (not a chain state).
7. A refused tree inside an attempt goes back to the bare root at once, before the attempt stops.
   Alternative: leave it until the next attempt redraws it, which a throw in between would find
   holding an empty leaf.
8. The entry returns the values in force, drawn or fixed, so R reads the fixed df and shape from it
   too. Alternative: NULL where nothing was drawn, R then reading the fixed values from the model.
9. The df and shape prior weights are read from a local `ResidualDfPrior` or `NBShapePrior`, whose
   constructors build them, through a new public `logPriorWeights()`. Alternative: an accessor to the
   response's own instance through `ResponseModel`.
10. A Student-t sampler is told by the entry's residual.df, since its model's family reads "gaussian".
    Alternative: a family field on the model.
11. A multinomial sampler draws category codes when every row holds one trial and one-trial count
    vectors otherwise, the reading of dec-B594's "fitted to one category per row" the stored data
    allows: a factor response is stored as one-hot counts, and nothing else records how it was given.
    Alternative: record the response's form at creation.
12. No floor on a drawn k (amendment 12, decided): a gamma draw that underflows to 0 at a df near 0.01
    gives non-finite draws in that share. Alternative: floor k at DBL_MIN in the prior draw.
13. The k draw passes the forest's leaf scale as `ChiKHyperprior::draw`'s scale argument rather than
    1.0, the plan's; at no leaves it multiplies 0, so the draw is the same either way.
14. Only the first forest's DART weights are reported; on a multinomial sampler each category forest
    draws its own. Alternative: a list of every forest's weights, which nothing reads.
15. The small-case DART oracle is the exact density, a closed-form recursion for P(held | s_1), rather
    than a simulated independent rejection compared by a two-sample KS: exact, no simulation error, and
    independent of the engine's code. Alternative: the simulated rejection.
16. The ordinal cutpoints are recovered through an offset grid at one x row, so each cutpoint reads
    back where its cumulative probability resolves in doubles; the latent is read at the row nearest
    one half. The k tests on probit truncate at a latent of +-7 for the same reason, against the
    truncated t. Alternatives: read the cutpoints off the entry directly (not the public result).
17. R's refusals run in the order: sampler class, several forests, improper leaf prior, type,
    n.samples, category offset, so a sampler no draw serves is named before its arguments.
    Alternative: the arguments first.
18. Fix round 1: the budget is forced by an unexported, undocumented option read per call, as
    `dbarts.sparseIndicators` is, and it sets the call's total rather than a per-sample count, so a
    test can run it out on the last sample. Alternatives: an internal argument of `drawPriorForests`,
    which a test of `samplePriorPredictive` cannot reach; a per-sample count, which cannot.
19. Fix round 1: a sample is held in a local before it enters the list. Alternative:
    `samples[i] <- list(...)`.
20. Fix round 1: the last-sample test runs on one column of two values, where every tree is held
    (each split leaves a row on either side and no node below has a cut), so a budget of 1 accepts
    the first sample and runs out on the second at every seed; a silent one-sample call checks that.
    Alternative: a seed that happens to run out there.
21. Fix round 1: the df and shape tests call the entry directly, since the public result reports
    neither. Alternative: infer them from "ppd", which the grid-mixture tests already do with too
    little power.
22. Fix round 1: under the restored finding-1 code the test file errors ("subscript out of bounds")
    rather than failing an expectation; an erroring file fails the suite and the check.
    Alternative: catch the error into an expectation.

Not in this slice: the dart() help's account of the prior the sampler runs (dec-B608's rule, first
sentence) lands with TODO dart-split-rule and dec-B605's change.

## Design-doc register

- docs/architecture.md, beside the empty-leaf veto: the prior draws condition on the same event;
  `Chain::drawTreeFromPriorOnce` is the shared predicate; `bartcore_sampleHyperparametersFromPrior` and
  the DART joint draw.
- man/samplePriorPredictive.Rd: the type item for each family, offset.test on a multinomial sampler,
  the hyperparameter draws and the improper refusal, the DART draw and its fallback, the seeding
  paragraph, the value's shapes.
- man/dbartsSampler-class.Rd: the see-also line names the hyperparameter draws; the "MEAN-forest entries
  by contract" passage still holds, the new entry being separate and bridge-only.
