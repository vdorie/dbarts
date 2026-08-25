# Mutation testing, leg C: untouched tinytest files, first half (b102e17c)

SCOPE the 102 tinytest files NOT in leg A's table (`ls inst/tinytest/test-*.R` minus leg A's
65), sorted alphabetically, files 1..51 - test-aft.R through test-ppd-sigma-pairing.R (the
exact list is mutation-C-evidence/scope-files.txt). Leg D takes 52..102, beginning at
test-predict-sparse.R. Priority order inside the half followed the brief: the six
generics/plot files first (no xbart file falls in this half; all seven are leg D's).

METHOD leg A's driver, catalogue format and severity scale, with one change that bought the
coverage: each mutation is scored against ITS TARGET FILE ONLY (install 9-25s, run 5-20s,
~15s a cycle against ~75s for the 51-file sweep), and a survivor whose finding depends on
"nothing else sees it" is re-run WIDE against all 51. Leg A's mtime fix is carried over
verbatim (restore touches, src/ restores touch every TU); a pristine build reproduced the 51
files at 0/0 before the first mutation and the driver restores between entries. 132 mutations
planted, 129 scored, 3 ANCHOR-FAIL and re-anchored. 6 wide escalations.

STATIC 1105 expect_ sites over the 51 files (1738 assertions at run time; the gap is loops).
content pin 88, shape 121, value-vs-value comparison 384, expect_true(all()) 130, smoke /
finiteness / inherits 136, refusal 246. The distribution matches leg A's. No tautology found.
Two structural notes: test-mutate-then-serialize.R has ZERO expect_ calls of its own (its 111
assertions all come from inst/common/stateContinuation.R's statesAgree), and
test-data-compatibility.R / test-data-formula.R are 7-of-8 and 11-of-13 expect_inherits -
the same class-only shape leg A flagged as unfalsifiable in test-bart-formula.R.

## Per-file table

n = expect_ sites; pin/shp/cmp/all/smk/ref as leg A. P = mutations aimed here; M = mutations
that moved an assertion here (aimed or incidental). OK >= 2 moved, WEAK = 1, GAP = 0.

file                                      n pin shp cmp all smk ref   P   M  verdict
aft                                      46   2   6  15   6   2  15   3   3  OK
assignInPlace-bounds                      9   0   0   0   0   0   9   3   3  OK
augmentation                             59   1   0  32   4   3  19   5   5  OK
bart-weights-parity                       6   0   1   4   0   0   1   3   3  OK
bart2-grow-from-root                      6   0   0   3   1   0   2   3   2  OK
bcf-forest-channel                       23   1   9  11   0   0   2   3   3  OK
bcf-zero-multiplier                       6   1   1   0   3   1   0   3   2  OK
binaryResponse-hyperprior                12   0   2   7   2   1   0   3   3  OK
calibration-prior-draws                  18   0   0  13   5   0   0   5   4  OK
composition-sequences                    52   4   1  27  12   6   2   1   1  WEAK
control-errors                           36   0   0   0   0   0  36   3   3  OK
control-valuesAreUsed                    11   1   0   9   1   0   0   6   3  OK
convergence-diagnostics                  44   2  21  12   0   8   1   3   2  OK
dart-mixed-columns                       13   0   4   0   8   1   0   2   1  WEAK
data-categorical-declared                21   2   4   0   9   3   3   6   2  OK
data-categorical-wide                    21   1   4   4   9   3   0   6   3  OK
data-categorical                         28   1   6  11   4   5   1   6   2  OK
data-compatibility                        8   0   0   1   0   7   0   3   1  WEAK
data-formula                             13   0   0   2   0  11   0   6   2  OK
data-handle                              21   0   3   3   3   2  10   2   2  OK
data-missing-matrix                      16   1   0   9   0   4   2   2   2  OK
data-missing                             26   1   3   7   4   7   4   5   4  OK
data-testOffset                          40  13   5  21   0   0   1   6   2  OK
embedding-recipes                        32   4   0  22   3   2   1   2   2  OK
factor-response                          24   6   0   6   0   3   9   3   2  OK
family-offset                            16   2   0   7   7   0   0   5   1  WEAK
family                                   25   8   0   9   2   4   2   3   4  OK
generics-correctValues                   13   0   0  12   0   1   0   3   3  OK
generics-errors                           5   0   0   0   0   0   5   3   3  OK
generics-intervals                       17   0   5   8   2   0   2   4   5  OK
generics-multithreaded                    2   0   0   2   0   0   0   5   1  WEAK
generics-posteriorPredictiveDistribution  6   0   1   2   3   0   0   2   2  OK
generics-sequentialExecution              2   0   0   1   0   1   0   1   0  GAP
gp-leaves                                30   0   4   6  10   2   8   4   3  OK
grouped-swap                             26   0   0   8   2   8   8   2   1  WEAK
grow-from-root-categorical                9   0   1   3   1   4   0   1   1  WEAK
hazard                                   37   6   3  10   4   4  10   3   2  OK
heteroscedastic-mutation                 36   5   2   9   8  11   1   1   1  WEAK
hurdle-surface                           17   2   1   2   0   4   8   2   2  OK
interactions                             16   2   0   3   0   3   8   3   3  OK
linear-leaves                            25   0   3   6   6   2   8   2   2  OK
model-errors                             19   0   0   0   0   0  19   5   1  WEAK
model-priors                             39   9   7   6   1   7   9   4   2  OK
monotone                                 29   9   3   6   6   1   4   4   2  OK
multinomial-counts-mutation              44   0   0  18   2   5  19   5   2  OK
multithreaded                            11   0   3   6   0   2   0   5   1  WEAK
mutate-sparse-valued                     17   0   1   5   0   6   5   1   1  WEAK
mutate-then-serialize                     0   0   0   0   0   0   0   4   1  WEAK
pdbart                                   12   0   0   8   0   2   2   4   3  OK
pointwise-loglik                         46   4  11  20   1   2   8   3   3  OK
ppd-sigma-pairing                        15   0   6   8   1   0   0   2   3  OK

OK 38, WEAK 12, GAP 1. No file went unreached; every one of the 51 ran under every mutation
aimed at it. Rate: 1 GAP in 51 (2%) against leg A's 4-in-10 spot check and 8-in-55 in-scope,
so the first half of the untouched population is NOT worse than the changed set - but the
WEAK band is wide (12 of 51), and the findings below are drawn from it as much as from the GAP.

## Findings

BLOCKER = the mutation left all 51 files in this half green AND names a defect in shipped
behavior. MAJOR = an aimed mutation the target file cannot see. MINOR = reach gap with a
killer elsewhere, or an equivalent mutant worth a decision.

1. BLOCKER  R/dbarts.R:1079 and R/generics.R:294  `n.threads` is a SILENT NO-OP on predict.
   The R5 method is `predict = function(x.test, offset.test, n.threads = control@n.threads)`
   and its body ends `.Call(C_dbarts_bartcore_predict, ptr, x.test, offset.test)` - three
   arguments, no thread count (`DEF_FUNC("dbarts_bartcore_predict", bartcore_predict, 3)`,
   R_interface.cpp:224). predict.bart coerces `n.threads <- as.integer(n.threads)[1L]` and
   hands it to that method, which discards it. Replacing the coerced value with -99L (gm1)
   and pinning the R5 default to 1L (gm2) both leave ALL 51 FILES GREEN on the wide run.
   test-generics-multithreaded.R exists solely to gate this ("test that predict gives same
   result when single or multi-threaded", 2 assertions) and is therefore vacuous: it never
   runs multi-threaded, and is a byte-for-byte duplicate of
   test-generics-correctValues.R:16-31 with `n.threads = 2L` written in the call.
   SHOULD EXIST either wire the argument through (a 4-argument bartcore_predict, or
   `sampler$setNumThreads()` around the call) and then pin it in
   test-generics-multithreaded.R with a timing- or state-free discriminator, OR drop the
   argument from both signatures and delete the file. Do not leave a documented argument that
   the engine cannot see.  VD-judgement (the fix is a public-signature decision, and
   inst/include/dbarts/dbarts.h consumers may reach the same entry)
2. MAJOR  R/plot.R:219  plot.pdbart's `for (i in xind)` can be narrowed to `xind[1L]` -
   rendering only the first requested predictor - with all 51 files green (wide).
   test-pdbart.R:88-91 states the claim in a comment ("the plot method renders each requested
   predictor into a null device") and then asserts `expect_silent(plot(pdb1))`, which sees
   neither how many panels were drawn nor which. Same shape as leg A finding 11.
   SHOULD EXIST test-pdbart.R:90: count the panels, e.g.
   `pdf(NULL); plot(pdb1); expect_equal(length(recordPlot()[[1L]]), <k>)`, or simplest,
   `expect_error(plot(pdb1, xind = 99L))` plus a per-panel `expect_silent(plot(pdb1, xind = 2L))`.
   agent-fix
3. MAJOR  R/bart.R:1236-1237  `n.grow.sweeps`'s COUNT is unpinned: `sampler$growFromRoot(
   n.grow.sweeps, ...)` -> `growFromRoot(1L, ...)` leaves all 51 green (wide). This is leg A
   finding 1's twin one level up - leg A pinned the sampler method, this is bart2's forwarding
   of the user's number. test-bart2-grow-from-root.R does gate the FEATURE (:107's early-RMSE
   claim kills a mutation that skips growth entirely) but not the count.
   SHOULD EXIST test-bart2-grow-from-root.R after :107: at one seed,
   `expect_false(isTRUE(all.equal(bart2(..., n.grow.sweeps = 1L)$yhat.train,
   bart2(..., n.grow.sweeps = 3L)$yhat.train)))`.  agent-fix
4. MAJOR  R/data.R:536-548  the 3+-level factor refusal in the single-forest entry points
   (`if (!splitMultinomialMessage) stop(caller, " does not fit a ", K, "-level ", ...)`) can
   be deleted with all 51 green (wide). test-factor-response.R:85-103 is the file that exists
   to hold it (five expect_error calls naming "multinomial"), and every one of them still
   passes: the SECOND stop() a few lines down (the family == "auto" arm) fires instead and
   carries the same word. The guard for the explicit-family callers - xbart, rbart_vi - is
   the one that goes unheld.
   SHOULD EXIST test-factor-response.R: anchor the refusals that must come from the
   non-auto arm on their own text, e.g. `expect_error(xbart(x, y3, family = "probit"), "does
   not fit a 3-level factor")`.  agent-fix
5. MAJOR  R/bartcore.R:407-414  `setResponse(updateScale = TRUE)`'s refusal under an
   amplitude-carrying or grouped sampler can be deleted with all 51 green (wide).
   test-grouped-swap.R:64-74 is written to hold exactly this ("updateScale = TRUE is refused
   on both response-side conduits, and the refusal names the two quantities it is
   protecting", three expect_error calls with pattern "tau") - and they still pass, because
   the C bridge refuses independently with a message that also contains "tau". The R-layer
   guard, which is the one that names the calibration anchor, is unheld.
   SHOULD EXIST test-grouped-swap.R:67: pin the R text that distinguishes the two,
   `pattern = "leaf calibration stated against the anchor"`.  agent-fix
6. MAJOR  test-generics-sequentialExecution.R (2 assertions, the leg's only GAP)  making
   `$sampleTreesFromPrior()` a complete no-op (se1: drop the `.Call`) leaves both assertions
   green, while it IS caught elsewhere (test-calibration-prior-draws.R:171-172). The file's
   :42 compares a bart2 fit against a hand-driven sampler loop, and BOTH sides call the
   mutated function, so the equality is common-mode blind to anything the two paths share.
   Its second assertion, :65 `expect_inherits(sampler, "dbartsSampler")`, cannot fail short of
   an error - it is the "sequential runs don't overflow with fixed trees" claim, unpinned.
   SHOULD EXIST test-generics-sequentialExecution.R:65: after the 6 sweeps past n.samples = 5,
   assert what "no overflow" means - `expect_equal(dim(sampler$getTrees()),
   dim(<the 5-sample store>))`, or that the store's recorded-draw count is 5 not 6.
   agent-fix
7. MAJOR  R/data.R:769  `validateXYOffset` can return `offset = NULL` - silently discarding
   the user's offset on the whole x/y interface - and test-data-compatibility.R (8
   assertions, 7 of them `expect_inherits(..., "dbartsData")`) and test-data-formula.R
   (13, 11 of them the same) both stay green. Their subject IS offset/weights/subset
   ingestion; class membership is invariant to whether any of it arrived.
   SHOULD EXIST test-data-compatibility.R: replace the inherits chain with content, e.g.
   `expect_equal(dbartsData(x, y, offset = offset)@offset, offset)` and
   `expect_equal(dbartsData(x, y, subset = 1:10, weights = weights)@weights, weights[1:10])`.
   The same rewrite for test-data-formula.R:14-55.  agent-fix
8. MINOR  R/plotTree.R:9-35  decodeCategoricalSplits' PADDING BRANCH IS DEAD CODE. Skipping
   the whole decode (dca3) and padding with "R" instead of "L" (dcd2) BOTH leave
   test-data-categorical.R, -declared.R and -wide.R green - including
   test-data-categorical-declared.R:53's `nchar(directions) == 4L` on a factor whose 4th
   level is never observed. The C side already emits one character per DECLARED level, so
   `padding > 0` never holds on any path these files reach. The one kill dcd1 earns
   (dropping the declared-count lookup) is incidental: `max(x[, j])` returns NA on the
   NA-carrying wide fixture and the file dies with "missing value where TRUE/FALSE needed",
   not on the claim.
   VD-judgement (delete the padding loop and the factorLevels lookup with it, or find the
   path that still needs it and give that path a test).
9. MINOR  R/model.R:1554 (resolvePriorScale's `node.prior@prior.sd * node.hyperprior@k`)
   is unreached by this half: dropping the `* k` conversion leaves all of
   test-calibration-prior-draws.R green because every sampler there spells `scale =`, never
   `sd =`. The sd spelling's own refusal (a drawn k) is tested; its arithmetic is not.
   SHOULD EXIST test-calibration-prior-draws.R: one constant-leaf arm with
   `node.prior = normal(k = fixedK, sd = priorSd)` measured against the same band as :57.
   agent-fix
10. MINOR  R/generics.R:125 (pointwiseLogLikelihood's zero-weight NaN flag) can report -Inf
    instead of NaN with test-pointwise-loglik.R green - the file's zero-weight arm asserts
    only that the value is not finite. The comment says the channel flags it "rather than
    reporting the -Inf an infinite sd would give", which is precisely the distinction lost.
    SHOULD EXIST test-pointwise-loglik.R: `expect_true(all(is.nan(ll[, zeroWeightRow])))`.
    agent-fix
11. MINOR  R/bart.R:2709 (`control@n.burn <- control@n.burn %/% control@n.thin` on bart()'s
    path) can be deleted with test-control-valuesAreUsed.R green; its keepevery arm
    (:40-55) pins only `nrow(yhat.train) == n.sims %/% keepevery`, the SAMPLE side. The burn
    side of the same thinning is unpinned on both bart() and bart2() (mt1 is the bart2 twin,
    also green everywhere).
    SHOULD EXIST test-control-valuesAreUsed.R: assert the sampler's resolved control, e.g.
    `expect_equal(fit$fit$control@n.burn, nskip %/% keepevery)`.  agent-fix
12. MINOR  Six reach gaps confined to this half - an out-of-file assertion is the sole
    killer, so no new assertion is needed, only a note not to trim the killer: dcm3 (x/y
    weights installed in reverse row order) survives test-data-compatibility.R and
    test-data-formula.R and dies only at test-bart-weights-parity.R:22 and
    test-generics-posteriorPredictiveDistribution.R:90; mp3 (named split probabilities
    resolving to nothing) survives test-dart-mixed-columns.R and dies at
    test-model-priors.R:86; me2/ml1/gl1 (the node-prior designation refusals) survive
    test-model-errors.R and die in test-gp-leaves.R / test-linear-leaves.R; fo1
    (getFitsWithoutOffset shifted) survives test-family-offset.R and dies in
    test-augmentation.R and test-embedding-recipes.R.  defer
13. MINOR  Equivalent mutants worth recording, not test gaps: bz2 (formForestVetoWeights'
    near-zero snap, leg A finding 2) again moved nothing - test-bcf-zero-multiplier.R, the
    one untouched file that names the snap, does not see it either, so leg A's finding 2
    stands unchanged and this file is not the missing killer. mts1 (stateFormatVersion 3 ->
    4) is inert because the writer and the reader move together, which is what
    minReadableStateFormatVersion is for. mth1 (bart2 collapsing n.chains) never reached
    test-multithreaded.R, which is a bart() file, so it is unproven rather than survived.
    defer

## Real defects exposed incidentally

ONE, and it is finding 1: `n.threads` is accepted and documented on `predict.bart` and on the
sampler's `$predict` and reaches nothing. It is not a wrong answer - the engine predicts on
the thread count the sampler was BUILT with - so the visible symptom is only that
`predict(fit, x, n.threads = 8)` does not use 8 threads. Reproduction:
`fit <- bart(x, y, keeptrees = TRUE, verbose = FALSE); body(fit$fit$predict)` shows the
three-argument `.Call`; `dbarts:::C_dbarts_bartcore_predict` is registered with 3 arguments.

Nothing else failed on the pristine tree: the 51 files reproduced at 0/0 before the run, and
every red assertion in the 129 scored entries traced to its planted mutation. Two mutations
killed by CRASHING the R session rather than failing an assertion (ai3, removing
assignInPlace's empty-index guard; dh1, removing the data-handle view's mutation refusal) -
both are kills, and both confirm the guards are load-bearing against a segfault, which is
what those two files say they are for.

## Not reached

Nothing in the 51-file half went unrun. Left for a later leg: per-file depth beyond the ~2.5
aimed mutations the budget bought (the WEAK band would move under a third aimed probe each);
wide escalation of the remaining survivors (6 of ~20 were escalated - dcm2, cvu1, cd3, pl2,
fo2, fo3, gl2, mo2, hz2, dca2, dfm2, dto2, dto3, dm2, mc2, mth3, me3 were scored against
their target file only, so their findings are stated as "the target file cannot see it", not
as "nothing sees it"); and inst/common/stateContinuation.R as a mutation TARGET, which is
where ALL 111 of test-mutate-then-serialize.R's assertions actually live.
