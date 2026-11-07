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
dart(), whole-forest rejection), dec-B631 (with no cap and no fallback, replacing dec-B608's), dec-B605
(the later split rule this draw carries over to).

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
Under dec-B631 every engine call that holds nothing ends this way, and R discards that forest, since
the next call redraws the whole. That is wasted work but harmless: a bare root is held with
probability at least 0.05, so the per-tree rejection's 10,000-attempt throw cannot occur in practice.
It is left as it is.

R (`drawPriorForests`) draws each sample by whole-forest rejection with no cap and no fallback
(dec-B631): it calls the entry with a bounded count of attempts until one call accepts, the entry having
then drawn the trees, so R skips `sampleTreesFromPrior`. The count starts at 1, doubles while a call
takes under 0.05 s and halves while one takes over 0.2 s, up to 2^20. A call is bounded by 2^20
attempts, about 9 s at 9 microseconds an attempt, though in practice it lasts one sample's time: where
a sample takes under 0.05 s the count climbs to the cap and each call ends at its first held draw. A
call that accepts nothing has also drawn k and the family's hyperparameters, which the next
call redraws; neither bears on acceptance, so the accepted sample's are still prior draws. After each call
the bridge polls for an interrupt (`userInterrupted`, R_CheckUserInterrupt under R_ToplevelExec, as a
run's poll is) once the engine's objects are gone and the sampler holds every tree it carries, and
raises it there (`raiseInterrupt`); an interrupt leaves only the private sampler, which is dropped, and
a time limit set by `setTimeLimit` is raised the same way. Once a second has gone into the attempts (the
threshold over 60), each call projects the call's length: the time so far times n.samples over the
samples accepted, or with none accepted yet, the time so far times n.samples, a lower bound. The first projection over 60 seconds warns, once per call, printed when
raised (`immediate. = TRUE`) so that it arrives while the draw runs rather than after the call returns, as at 60 rows and 500 trees for 3,000 samples, after a second: "samplePriorPredictive: under
dart() 19 of 69,426 joint draws of the split weights and trees so far were ones the sampler holds
(acceptance 0.00027), so the call is projected to take about 3 minutes; it can be interrupted", or with
none accepted, as at 30 rows and 500 trees for 200 samples, "none of 98,303 joint draws ... so far was one the sampler
holds (acceptance none in 98,303), so the call will take at least about 4 minutes; it can be
interrupted". An unexported, undocumented option, `dbarts.priorPredictiveWarnSeconds`, sets the
threshold, so a test can force the warning at any seed. Acceptance measured on this build: 0.01 at 60
rows, 8 columns and 200 trees; 0.035 at 200 rows; 0.36 at 1000 rows.

Tests: test-prior-predictive.R, the small case - x1 binary, x2 on {0, 1, 2} (two equal cuts), four of the
six cells occupied, 20 trees, alpha fixed at 2 so s_1 is uniform under the stated prior. There
P(every tree held | s_1) has a closed form, a recursion over the node's open columns and cuts with
growth 0.95 / (1 + depth)^2, so the exact density of s_1 under whole-draw rejection is
P(held | s_1)^20 on a grid of 1,000 points. 2,000 entry draws, every one accepted, against its CDF by KS;
the same draws against the stated uniform, p below 1e-6 (the per-tree retry's distribution). Under one
seeded control, dart()'s "ev" differs from the same sampler's at uniform split.probs. At 60 rows, 8
columns and 200 trees, about one in two hundred held, 400 samples of a negative-binomial sampler are
each an accepted joint draw, over more than 8,000 attempts, and their k follows its chi(1.5, 2) prior
(KS) and their shape the grid's weights (chi-squared), so the hyperparameters are redrawn on every
engine call however large the batches grow. With the warning threshold forced to 0 seconds a call
warns once, in the words above. On the small case the default threshold stays silent. With the
interrupt poll armed to fire once, a call is interrupted (class "interrupt"), and a seeded call
afterward returns what it did before. Where no forest is ever held (two points on the diagonal, 200
trees, each tree held with probability about 0.6), a call armed to interrupt at its twelfth poll is
interrupted, so a call holding nothing never ends a sample; and at a threshold of 3 seconds the first
projection, the time so far times n.samples, warns in well under 1.5 s, "none of" and "at least
about". The fallback's tests went with the fallback; the last-sample case of fix round 1 went too,
since with no cap every sample is accepted and none can come back empty.
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

Run on arm64 macOS, at most 2 cores. Fix round 2 reran, after the rebase onto 1aaba0199 and on
`--preclean` installs of both libraries: tests/cpp plain and under ASan/UBSan, the full suite, the
trio, the snapshot files, the extra seeds, the lint set and verify-anchors; the numbers below are that
round's. The R-loaded ASan run is the build round's. Fix round 3 changed only R, the tests, the help and
these records, so it reran the full suite, the extra seeds, the lint set, verify-anchors and R CMD
check, whose numbers are given as that round's; the engine and bridge are unchanged, so tests/cpp, the
trio and the snapshot files stand at fix round 2's.

- tests/cpp, shipped config: plain 559 ok, all passed; under `-O2 -g -fsanitize=address,undefined`
  (`detect_container_overflow=0`) 559 ok, 0 sanitizer reports. The new suites: the DART prior draw's s_1
  KS p 1.00 at a fixed alpha and 0.84 under a mask, alpha's grid chi-squared p 0.849, E[s_1^2] 0.1081
  against 0.1086; k's quantiles at worst z 1.82; the alloc harness's DART arms at budgets 1000 and 1 over
  every leaf type, 14 to 17 allocations each, all failing cleanly.
- Full tinytest suite, shipped, one process: fix round 3, by `run_test_dir` (at home), 42150 results
  in 304 files, 0 failed; fix round 2, 42049 in 303. test-prior-predictive.R 83 results,
  test-setters-draw-nothing.R 121.
- Equivalence trio on the reference build, `--bitwise` (`--strict-coverage` on equivalence.R): 55, 19
  and 11 scenarios identical against equivalence-38da480c.rds, bcf-equivalence-38da480c.rds and
  multinomial-equivalence-80b1c8d4.rds, 0 skipped.
- The four snapshot files on the reference build: 7, 3, 14 and 3 results, 0 failed.
- Each KS and chi-squared test of test-prior-predictive.R at its own seed, 2611, and at 7193, 40417
  and 91 (fix round 3): 18 tests at 4 seeds, every p > 1e-3, the smallest 0.0125 (the nbinom k test);
  the DART block's k KS p 0.41 to 0.83 and shape chi-squared p 0.52 to 0.92. The one test that
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
- R CMD check --as-cran on the built tarball, fix round 3: 1 NOTE, the incoming-feasibility one, not
  this slice's (an expired certificate at https://www.rob-mcculloch.org in man/bartBT.Rd and
  man/pdbart.Rd, and the Date field over a month old); tests OK in 245 s. The build round's second
  NOTE, `rt` used without an import, was fixed then.
- The warning, printed when raised (fix round 3): at a threshold of 3 seconds, a 2,000-sample call at
  60 rows, 8 columns and 200 trees printed it at 0.07 s and returned at 8.2 s.
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
| finding 1 restored: a NULL sample deletes its element (R; fix round 1, since removed) | test-prior-predictive.R | file errors |
| the forced budget ignored (R; fix round 1, since removed) | test-prior-predictive.R | 1 of 81 |
| fix round 2: warns on every engine call (R) | test-prior-predictive.R | 1 of 77 |
| fix round 2: never warns (R) | test-prior-predictive.R | 2 of 77 |
| fix round 2: a cap of 50 attempts per sample, pooled (R) | test-prior-predictive.R | 1 of 77 |
| fix round 2: the bridge's interrupt poll removed (C++) | test-prior-predictive.R | 1 of 77 |
| H1: k, the df and the shape drawn only on calls of at most one attempt (C++) | test-prior-predictive.R | 2 of 83 |
| A1: a refused tree after the first kept (C++) | test-prior-predictive.R | 5 of 83 |
| B1: an empty call's per-tree forest kept once the batch is 64 or more (R) | test-prior-predictive.R | 1 of 83 |
| W1: the no-acceptance projection without the n.samples factor (R) | test-prior-predictive.R | 1 of 83 |

## Measured for the open calls

Ruled: dec-B631 settled calls 1 and 2 on this evidence (always draw exactly; no fallback). Measured
for the maintainer's calls on the fallback's budget (call 1) and its burn-in and thinning (call 2), on
the fix-round-1 code. Eight uniform columns, a gaussian response, 200 trees at 60 rows (the collapse
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

1. Settled by dec-B631: no cap. As built before the ruling, the fallback switched the whole call once a
   budget of 50 whole draws per requested sample, pooled over the call, ran out. Alternatives: a fixed
   total cap; an argument; a cap of about 10,000 attempts per sample.
2. Settled by dec-B631: no fallback, so no burn-in or thinning. As built before the ruling, the
   fallback's burn-in and spacing were the control's `n.burn` and `n.thin`. Alternatives: fixed
   counts; arguments.
3. Settled by dec-B631: the fallback's warning went with it; a warning of the projected time takes
   its place (call 23).
4. Gone with the fallback, which drew k, the family's hyperparameters and the leaves afresh for each
   sample rather than reading them off the run.
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
18. Gone with the cap (fix round 2): fix round 1 forced the budget by an unexported option.
19. Gone with the cap (fix round 2): fix round 1 held a sample in a local before it entered the list,
    since a NULL sample deleted its element; with no cap no sample is NULL, and the samples are built
    by `lapply`.
20. Gone with the cap (fix round 2): fix round 1's last-sample test, on a design where every tree is
    held.
21. Fix round 1: the df and shape tests call the entry directly, since the public result reports
    neither. Alternative: infer them from "ppd", which the grid-mixture tests already do with too
    little power.
22. Gone with the cap (fix round 2): the finding-1 mutant's erroring test file.
23. Fix round 2: the warning's rule. Projections start once a second has gone into the attempts (the
    threshold over 60, so a forced threshold of 0 projects at once), then follow every engine call; the
    first over the threshold warns, once per call. The projection is the time so far times n.samples
    over the samples accepted, a rate over the whole call, or with none accepted, the time so far times
    n.samples, stated as a lower bound with acceptance "none in" the attempts made, never a
    division by zero. Fix round 3: the warning is printed when raised (`immediate. = TRUE`), since R
    otherwise holds it until the top-level call ends. Alternatives: project from the first call (noisier, a slow first sample warning
    a call that is not long); a fixed number of attempts before the first projection, which scales
    with neither the attempt's cost nor the threshold; a Bayesian bound on acceptance with none
    accepted.
24. Fix round 2: the attempt loop runs in R over engine calls of a bounded count, sized to take 0.05
    to 0.2 s and capped at 2^20 attempts (about 9 s, in practice one sample's time), the bridge polling
    after each, rather than an uncapped loop in the engine with
    a poll inside: each call ends at a state the sampler holds, the engine stays R-agnostic and its
    entry and tests are unchanged, and an interrupt waits at most one call (12 ms measured). Alternatives: the engine
    loop with a poll callback, as a run takes; R's own interrupt check between calls, which comes only
    every thousand or so evaluations.
25. Fix round 2: the threshold is an unexported, undocumented option read per call, as
    `dbarts.sparseIndicators` is, default 60 seconds. Alternatives: an internal argument; a
    documented argument.
26. Fix round 2: forcing verbose off in the private control existed only for the fallback's run and is
    removed. Alternative: keep it, which nothing now needs.
27. Fix round 2: the cap test reads each sample's forests.drawn through `drawPriorForests`, which the
    public result does not report. Alternative: infer the cap from the split weights' concentration,
    which needs thousands of samples. Fix round 3: the same block reads k and the shape on a
    negative-binomial sampler, the hyperparameters' only test under DART.
28. Fix round 3: the batch and the no-acceptance projection are tested on a design where no forest is
    ever held, two points on the diagonal at 200 trees, so both tests are deterministic. The
    projection's test times the warning against a 3-second threshold. Alternatives: a low-acceptance
    design (each test then holds with a probability below 1); read the projected seconds off the
    message (its units change with the value).
29. Fix round 3: call 6's per-tree draw, which R discards on every call that holds nothing, is left in
    place. Alternative: leave the bare-root reset as the settled state and zero the totals, saving the
    work.

Not in this slice: the dart() help's account of the prior the sampler runs (dec-B608's rule, first
sentence) lands with TODO dart-split-rule and dec-B605's change.

## Design-doc register

- docs/architecture.md, beside the empty-leaf veto: the prior draws condition on the same event;
  `Chain::drawTreeFromPriorOnce` is the shared predicate; `bartcore_sampleHyperparametersFromPrior` and
  the DART joint draw.
- man/samplePriorPredictive.Rd: the type item for each family, offset.test on a multinomial sampler,
  the hyperparameter draws and the improper refusal, the exact DART draw, its time, warning and
  interrupt, the seeding
  paragraph, the value's shapes.
- man/dbartsSampler-class.Rd: the see-also line names the hyperparameter draws; the "MEAN-forest entries
  by contract" passage still holds, the new entry being separate and bridge-only.
