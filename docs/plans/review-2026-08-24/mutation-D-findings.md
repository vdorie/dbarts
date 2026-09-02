# Mutation testing, leg D: the untouched tinytest files, second half (b102e17c)

SCOPE the 102 files leg A did not table (`ls inst/tinytest/test-*.R` minus its 65), sorted
alphabetically, entries 52-102 - 51 files, test-predict-sparse.R through test-zero-weights.R.
Leg C takes 1-51 (ending at ppd-sigma-pairing). Priority inside the half followed leg A's
recommendation: the xbart loss arms and the report/generic surfaces first.

METHOD leg A's harness reused verbatim - drive.py, runfails.R, classify.py and its
touch-on-restore mtime fix. Each mutation was applied to a `git archive HEAD` copy under the
scratch dir, installed into a private lib, and run against all 51 files; the checkout was
never touched. A pristine build reproduced the 51 at 1443 assertions / 0 failures before the
first mutation; a second, independent pristine tree reproduced the full 167-file suite at
6946 / 0.

BUDGET STOP the machine was shared with leg C's battery and three SBC arms for the second
half of the run, and throughput fell from ~25 s/mutation to ~5 min/mutation. The leg stopped
at the 5-hour mark with a queue of 30 prioritized mutations partly run; the files that queue
was aimed at and did not reach are named under "Not reached". The three full-suite replants
of zero-killer mutations were cancelled for the same reason - severities below are therefore
scoped to the 51 files plus, where stated, direct inspection of tests/cpp.

STATIC every expect_ site classified by reach (classify.py, leg A's classes folded into leg
A's six columns): 684 sites over the 51 files - content pin 51, shape 88, value-vs-value 200,
expect_true(all()) 63, finiteness/typecheck/silence smoke 136, refusal 146. NO tautology of
any kind: zero self-comparisons, zero `is.character(capture.output(print(x)))`-style
cannot-fail smoke (leg A's finding 11 pattern does not recur here), and the two files that
capture print output (validate-composition, sampler-degenerate-cuts) both grep the text.

## Per-file table

n = expect_ sites; pin/shp/cmp/all/smk/ref as in leg A. P = mutations aimed here that were
scored; M = mutations that moved an assertion here; Mn = the same excluding "global" killers
(a mutation that moves more than 10 of the 51 files proves nothing about one file's reach).
Verdict is on Mn: OK >= 2, WEAK = 1, GAP = 0.

file                                                    n pin shp cmp all smk ref   P   M  Mn  verdict
predict-sparse                                         20   0   1  12   0   2   5   3   1   0  GAP
prior-init-empty-leaves                                 3   1   0   0   0   2   0   2   3   2  OK
prior-predictive                                       20   0   7   5   3   2   3   3   3   2  OK
rbart-aft                                              33   3   4   7   9   4   6   0   5   4  OK
rbart-bartcore                                         16   2   5   1   3   3   2   1   6   5  OK
rbart-custom-prior                                      2   0   0   0   0   2   0   0   3   2  OK
rbart-error                                             6   0   0   1   0   0   5   2   3   2  OK
rbart-example                                          16   0  12   4   0   0   0   0   3   2  OK
rbart-generics                                         15   0   0  14   0   1   0   2   5   4  OK
rbart-loop-callback                                     6   0   0   6   0   0   0   0   4   3  OK
rbart-multithreaded                                     1   0   0   0   0   1   0   1   1   0  GAP
rbart-options                                          15   3   2   2   1   6   1   5   4   3  OK
rbart-performance                                       3   0   0   1   0   2   0   0   3   2  OK
rbart-reproducibility                                   7   0   0   6   0   1   0   0   1   0  GAP
rbart-weighted-binary                                  10   0   0   0   7   2   1   0   1   0  GAP
reproducibility-binaryResponse                          7   0   2   5   0   0   0   0   0   0  GAP
reproducibility-continuousResponse-multithreaded        3   0   0   3   0   0   0   1   1   0  GAP
reproducibility-continuousResponse-singleThreaded      14   0   3  11   0   0   0   0   1   0  GAP
reproducibility-rbart                                   1   0   0   1   0   0   0   0   3   2  OK
reproducibility-xbart                                   3   0   0   2   1   0   0   1   7   6  OK
rng                                                     7   0   0   6   0   1   0   1   2   1  WEAK
robust-errors                                          22   4   6   1   1   3   7   2   3   2  OK
sampler-degenerate-cuts                                15   1   3   1   7   1   2   1   2   1  WEAK
sampler-model                                           3   0   0   0   0   3   0   2   2   1  WEAK
sampler-offset                                         22   0   7  15   0   0   0   3   3   3  OK
sampler-prior                                           1   0   0   0   0   1   0   0   1   0  GAP
sampler-prior-midchain                                  6   0   0   6   0   0   0   1   2   1  WEAK
sampler-residuals                                       2   0   0   2   0   0   0   0   1   0  GAP
sampler-saveLoad                                        2   0   0   1   0   0   1   1   2   1  WEAK
sampler-setData                                         1   0   0   1   0   0   0   1   1   0  GAP
sampler-setPredictorPerObservation                     30   0   3   8   5  10   4   0   1   0  GAP
sampler-splitProbabilities                             17   0   6   3   3   0   5   3   3   2  OK
sampler-state-emptyLeafVeto                             4   0   0   2   0   2   0   1   0   0  GAP
sampler-updatePredictorPerObservationJointly           44   0   3  16   4  14   7   1   2   1  WEAK
sampler-updateState                                    11   0   0   2   0   9   0   1   3   2  OK
serialization-storeState                                5   0   0   3   0   1   1   1   4   3  OK
simd                                                    3   0   0   3   0   0   0   0   1   0  GAP
slice-sample                                           22  10   1   0   5   2   4   4   3   3  OK
spec                                                   30  15   4   4   0   2   5   2   4   3  OK
sum-to-one-tolerance                                    6   0   0   0   0   3   3   0   1   0  GAP
utility-chains                                         10   0   2   8   0   0   0   4   4   4  OK
validate-composition                                   55   7   3   6   3  23  13   4   4   4  OK
warm-start                                             14   0   0   4   3   2   5   1   1   0  GAP
weighted-logistic                                      10   1   1   0   2   3   3   0   0   0  GAP
xbart-error                                            52   1   0   1   0   0  50   2   5   4  OK
xbart-loss                                              8   0   2   2   0   4   0   4   3   2  OK
xbart-method                                           13   0   4   2   1   5   1   3   3   2  OK
xbart-model                                            37   0   4   6   0  15  12   2   6   5  OK
xbart-oracle                                           16   3   1  10   1   1   0   4   8   7  OK
xbart-reproducibility                                   7   0   2   2   2   1   0   1   6   5  OK
zero-weights                                            8   0   0   4   2   2   0   1   2   1  WEAK


OK 28, WEAK 7, GAP 16
mutations scored: 66


Most GAPs above are UNPROVEN rather than proven bad: 10 of the 16 took no aimed mutation at
all before the budget stop (see "Not reached"). The six that DID take an aimed probe and
survived it are predict-sparse, rbart-multithreaded, reproducibility-continuousResponse-
multithreaded, sampler-setData, sampler-state-emptyLeafVeto and warm-start - findings 1, 2,
4, 5, 6 and 8 below.

## Findings

Severity here is leg A's, adjusted for the budget stop: BLOCKER would require a full-suite
plus tests/cpp replant, which this leg could not afford, so nothing is filed as one.

1. MAJOR  [[src/bartcore/moves.hpp:120-123@b102e17c]]  resolveVetoRank's `-HUGE_VAL` can be replaced by
   the pre-fix finite `-1.0e7` with all 51 files green. test-sampler-state-emptyLeafVeto.R
   exists for exactly this ("Regression for the -1e7 -> -HUGE_VAL fix") and none of its 4
   assertions move: its fixture (n = 2000, resid.prior = fixed(1e-6), 50 trees) no longer
   drives a valid branch's score past -1e7, so the regime the file's header describes is not
   reached. [[tests/cpp/test_moves.cpp:1212-1218@b102e17c]] pins `currentLogL == -HUGE_VAL` directly, so
   the C++ gate does own the invariant - the R file's claim to be the regression does not
   hold.  SHOULD EXIST  either scale the fixture (smaller sigma / deeper trees) until the
   assertion at [[tests/cpp/test_moves.cpp:43-48@b102e17c]] actually separates, or drop the regression claim from the header and
   point it at tests/cpp.  VD-judgement

2. MAJOR  [[R/xbart.R:367-371@b102e17c]], [[test-xbart-method.R:114-132@b102e17c]]  the k-fold remainder distribution
   (`rep.int(c(1L, 0L), c(n %% k, k - n %% k))`) can be deleted - every fold takes
   floor(n / k) rows and the remainder rows are never held out - with test-xbart-method.R
   green. That includes its own section headed "test that k-fold subdivides data correctly
   when data do not divide evenly by k", whose single assertion is
   `expect_inherits(xval, "array")`. Only [[test-xbart-oracle.R:141@b102e17c]] catches it.
   SHOULD EXIST  [[test-xbart-method.R:132@b102e17c]], a fold-size oracle rather than a class check:
   `expect_equal(xbart(x, y, method = "k-fold", n.test = 5, n.reps = 1L,
   loss = function(y.test, s, w) length(y.test), ...), 24 / 5)`.  agent-fix

3. MAJOR  [[R/xbart.R:712@b102e17c]] and [[R/xbart.R:731@b102e17c]], test-xbart-loss.R  the file whose subject is the custom
   loss channel survives BOTH "the loss function's value is discarded"
   (`lossValues[cell, ] <- 0.0`) and "the k-fold losses are summed instead of averaged". Its
   8 assertions are inherits / dim / dimnames / !anyNA, twice; not one reads a loss value.
   Its n.threads = 1 and n.threads = 2 arms are likewise indistinguishable (finding 6).
   SHOULD EXIST  [[test-xbart-loss.R:41@b102e17c]], a constant-loss identity
   (`expect_equal(as.vector(xbart(..., loss = function(a, b, c) 42)), rep(42, ...))`) and a
   mad() value oracle mirroring [[test-xbart-oracle.R:63-70@b102e17c]].  agent-fix

4. MAJOR  [[R/dbarts.R:836@b102e17c]], [[test-prior-predictive.R:104-111@b102e17c]]  samplePriorPredictive's sigma
   prior draw (`sqrt(df * sigest^2 * rawScale / rchisq(n.samples, df))`) can be replaced by
   the point estimate `rep_len(sigest, n.samples)` with all 51 green: the file's (d) block
   asserts only `var(ppd.pinned) >= var(pinned1)`, which holds for any positive noise scale.
   The prior-predictive's whole point is that sigma carries its own prior uncertainty.
   SHOULD EXIST  [[test-prior-predictive.R:111@b102e17c]], at the pinned control seed
   `expect_true(sd(apply(ppd.pinned - pinned1, 1L, sd)) > 0)` - under the mutation every draw
   shares one sigma and that spread is 0.  agent-fix

5. MAJOR  [[R/utility.R:797-816@b102e17c]], test-predict-sparse.R  the sparse test column's level-recode
   arm has NO live coverage. Every passing sparse case in the file declares the training
   level set, so `identical(oldLevels, trainingLevels)` short-circuits before the recode: two
   mutations inside it - category codes off by one, and the reference level left on the
   training coding - leave all 51 green, and so does deleting its unseen-level refusal, whose
   error the engine raises independently. The file's headline claim ("bitwise against the
   dense twin") is therefore tested only on the arm that does no work.
   SHOULD EXIST  test-predict-sparse.R, a sparse test frame whose factor declares MORE levels
   than training but uses only training ones, predicted `expect_identical` against its dense
   twin.  agent-fix

6. MINOR  test-rbart-multithreaded.R (1 assertion, `expect_inherits(..., "rbart")`)  forcing
   dbartsControl's n.threads to 1L leaves all 51 green, test-reproducibility-
   continuousResponse-multithreaded.R included - its pinned draws are thread-count
   independent by design ([[test-rng.R:133-135@b102e17c]] says so). Nothing in the half can observe a lost
   thread, and this file cannot fail for any reason short of an error.
   SHOULD EXIST  fold it into test-rbart-options.R, or pin the reachability directly:
   `expect_equal(fit$fit[[1L]]$control@n.threads, 2L)`.  VD-judgement

7. MINOR  [[R/model.R:1614-1627@b102e17c]]  `student()`'s own `df <= 0.0` refusal is dead: deleting it
   leaves all 51 green because `newValidated` runs dbartsStudentDist validity
   ([[R/A_class.R:223@b102e17c]]), whose message also contains "positive finite". test-robust-errors.R
   [[R/A_class.R:28-31@b102e17c]] pins the S4 validity, not the constructor. A redundancy, not a coverage hole.
   defer

8. MINOR  [[R/bartcore.R:441@b102e17c]] (setData), test-sampler-setData.R  quantizing the replacement data
   onto a single cut per column leaves the file's ONE assertion green: it compares
   `sd(rowMeans(samples1$train) - y)` against `sd(rowMeans(samples2$train) - y)` at tol 1e-2,
   and both fits degrade together, so the comparison is self-cancelling. The file's stated
   subject - "setData yields a valid model when redefining cut points" - has no pin.
   SHOULD EXIST  test-sampler-setData.R, an absolute band on the post-setData fit
   (`expect_true(sd(rowMeans(samples2$train) - y) < 0.3)`) alongside the paired comparison.
   agent-fix

9. MINOR  [[R/model.R:378-380@b102e17c]]  the uniform-split-probability canonicalization to `numeric()` is
   unreached by test-sampler-splitProbabilities.R: with `if (FALSE)` all 51 stay green,
   because its three `length(splitProbabilities) == 0L` assertions ([[R/model.R:21-41@b102e17c]]) take the earlier
   scalar/NULL path instead.  SHOULD EXIST  the same assertion for an explicitly uniform
   length-p vector (`splitprobs = rep(1, ncol(x))`).  agent-fix

10. MINOR  [[R/rbart.R:1201@b102e17c]] vs [[R/rbart.R:1269@b102e17c]]  the two arms of rbart_vi's `$varprobs` assembly are
    separately unpinned: [[test-rbart-options.R:167-168@b102e17c]] exercises the DART report at
    n.chains = 1, so the mutation aimed at the multi-chain arm moved nothing. The
    single-chain replant was queued and did not run inside the budget.  defer

11. MINOR  [[R/sliceSample.R:185@b102e17c]]  dropping the upper-boundary clamp on the stepping-out
    interval is an EQUIVALENT mutant here (the beta target is 0 outside [0, 1], so every
    out-of-range proposal is shrunk away). Recorded so a later leg does not re-file it.
    defer

## Real defects exposed incidentally

One, and it is a robustness gap rather than a live defect. [[R/generics.R:2040-2046@b102e17c]]
(fitted.rbart, type = "ev") reads `ranefNames <- dimnames(object$ranef)[[length(...)]]` and
passes `match(object$group.by, ranefNames)` straight into `.Call(C_rbart_fitted, ...)`. When
that name vector is absent the match is all-NA and the C entry point indexes on NA_INTEGER:
the R session SEGFAULTS ("caught segfault ... cause 'invalid permissions'", traceback
`1: fitted.rbart(rbartFit)`) rather than raising an error. Shipped code always populates the
dimnames, so this is not reachable from the public surface today; it is reachable from any
consumer that builds an rbart-shaped list by hand, which the flat C API invites.
REPRODUCTION: in [[R/bart.R:40@b102e17c]], `colnames(res) <- dimnames(samples)[[1L]]` -> `NULL`; install;
`Rscript -e 'tinytest::run_test_file("inst/tinytest/test-rbart-performance.R")'`. Suggested
fix: bounds-check the group index in src/R_interface_rbart.cpp's rbart_fitted, or refuse an
NA match R-side before the .Call.  agent-fix

## Not reached

The 10 GAP files below took NO aimed mutation before the budget stop; their verdict is
"unproven", not "proven bad". Their aimed mutations are written and anchor-checked in
mutations.jsonl (ids in parentheses) and can be run as-is by a follow-up:
rbart-reproducibility (rr1), rbart-weighted-binary (wl3), reproducibility-binaryResponse and
-continuousResponse-singleThreaded (no aimed probe was written - any RNG-shifting mutation
serves), sampler-prior (pr1), sampler-residuals (sr1),
sampler-setPredictorPerObservation (po1, po2), simd (sm1), sum-to-one-tolerance (st1, st2),
weighted-logistic (wl1, wl2). Warm-start's ws1/ws2/wsr and the joint file's jt1/jt2 were
still running at the stop; their results, if any, are in run-results.txt.
Also cancelled: the full-suite (167-file) replants of the 12 zero-killer mutations, which
would have settled whether findings 1, 4, 5 and 8 are BLOCKERs rather than MAJORs.
