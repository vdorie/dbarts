# Mutation testing, leg A: the R test suite (b102e17c)

SCOPE 55 tinytest files changed since 03b97db7 (`git diff --name-only 03b97db7 HEAD --
inst/tinytest`, 56 paths, one of them capi/consumer.c) plus 10 untouched files sampled at
seed 20260824: bart-formula, bcf-reporting, data-mixed-mutation, data-sparse,
heteroscedastic, plot-generics, prior-init-composed-law, sampler-trees,
weighted-binary-ppd, xbart-weight.

METHOD 85 mutations planted, 84 scored, 1 wedged the runner (b04). Each was applied to a
`git archive HEAD` copy under the scratch dir, installed into a private lib, and run against
all 65 files; the checkout was never touched. A mutation that moved no assertion in any of
the 65 was replanted against the full 167-file suite, and if it still moved nothing and lived
in src/, against tests/cpp. Two independent pristine builds reproduced 167 files at 0/0.

HARNESS NOTE the first pass restored mutated files with mtime-preserving copies, so R CMD
INSTALL skipped recompilation and an R-only mutation inherited the previous C++ one. Caught
from an impossible kill pattern (an R/data.R mutation failing multinomial fits), fixed by
touching on restore, every affected entry re-run. benchmarks/R/mutation-battery.R is NOT
exposed: it archives and --precleans per entry.

STATIC every expect_ site classified by reach: 2917 assertions in the 55 in-scope files,
3363 over all 65 - content pin 233, shape 299, value-vs-value comparison 1042,
expect_true(all(...)) 250, finiteness/typecheck/silence smoke 392, refusal 634. THREE reach
classes the memo expected are NOT here. Tautologies: one candidate in 3363 sites (a
NULL-vs-NULL expect_equal), no self-comparison, no compare-to-same-call-path pair. RNG-locked
snapshot literals: all five drift-tripwire files are untouched and out of scope; the 11
high-precision literals in scope carry their derivation inline ([[test-bcf-family.R:150-174@658869ac]]).
Vacuous all(): an instrumented run shadowing base::all recorded ZERO zero-length operands
across all 167 files - the 250 sites are a latent hazard, not a live one.

## Per-file table

n = expect_ sites; pin/shp/cmp/all/smk/ref = content pin, shape, value comparison,
expect_true(all()), smoke, refusal; P = mutations aimed here; M = mutations that moved an
assertion here (aimed or incidental). OK >= 2 moved, WEAK = 1, GAP = 0.

file                                  n pin shp cmp all smk ref    P    M  verdict
active-rows-pins                     45   0   0  34   5   0   6    4    4  OK
argument-surface                    117  21  18  18   5  18  20    2    2  OK
bart-bart2                           32   4   1   3   0   4  14    1    0  GAP
bartcore                            182   8  29  61  30  28  24    1    3  OK
bcf-creation                        132  17  16  29   4   6  59    3    4  OK
bcf-family                           65  19   4  18   6   5  13    2    1  WEAK
bcf-mutation-pins                    20   0   0   4   3   7   6    3    0  GAP
bcf-r5-surface                       45   0   6  16   5   8  10    3    4  OK
bcf                                  74   1   5  30  22   6  10    3    1  WEAK
blocks                               25   0   1   4   0  10  10    1    1  WEAK
boundary-inputs                      13   0   2   0   9   1   1    1    1  WEAK
calibration-creation                 40   8   2  13   3   4  10    4    2  OK
calibration-midchain                110  11   5  48  16  10  20    4    4  OK
capi                                233  40  19  78  26  28  40    0    8  OK
data-errors                          17   0   0   0   0   1  16    4    2  OK
data-mixed                           75   5   8  29   3   9  13    6    3  OK
data-precision-warning               13   0   3   1   0   7   0    1    0  GAP
dispersion-channel                   66   1  28  15   9  13   0    1    1  WEAK
empty-leaf-veto-weights              12   0   1   1   8   2   0    3    1  WEAK
error-quality                         8   0   1   0   0   0   6    3    0  GAP
fits-without-offset                  44   0   1  27   0  14   2    0    2  OK
forest-basis-r5                      88  21  13  28   1  11  14    3    4  OK
forest-basis-subset                  15   0   1  10   0   0   4    3    2  OK
forest-weights-r5                    16   0   0   5   0   4   7    3    5  OK
forest-weights                       28   0   1   5   8   3  11    2    1  WEAK
formula-terms                        29   0   2   8   0  16   3    0    1  WEAK
grow-from-root                       18   0   1   4   5   5   3    3    0  GAP
heteroscedastic-channels             24   0   2  17   0   5   0    3    1  WEAK
heteroscedastic-warm-start           23   0   1   9   5   5   3    1    2  OK
host-shell-pins                      22   3   7   7   0   4   1    1    0  GAP
hurdle                               41   0   6  20   6   4   4    3    3  OK
logistic-weight-swap                 39   1   1  17   3   8   9    2    1  WEAK
makeModelMatrix                      91   8  20  49   3   1  10    1    1  WEAK
multi-forest-seam                    67   0   4  13  11  24  15    1    1  WEAK
multinomial-category-offset          50   0   0  22   4   7  17    1    2  OK
multinomial-generics                 48   1   7  18   2  11   9    1    2  OK
multinomial-r5-surface              117   4  12  32   1  15  52    2    2  OK
multinomial-surface                 121   1  10  59  16   7  28    0    1  WEAK
multinomial-test-offset              55   0   1  21   3  10  20    1    2  OK
multipleAssignment                   27  16   3   0   0   0   5    2    2  OK
nbinom                               86   5   6  42   7   9  16    3    3  OK
ordinal                              87   9   7  42   7   8  12    4    1  WEAK
pdbart-keeptrees                      6   0   2   2   2   0   0    3    1  WEAK
predict-blend                        49   0   2  33   0   1  13    7    8  OK
predict-forest                       31   0   7  16   2   0   6    3    2  OK
rbart-groupby                        26   0   5  14   1   1   0    1    0  GAP
rbart-weights                         3   0   1   0   0   0   0    1    0  GAP
sampler-bridge-errors                 7   0   0   0   0   1   6    3    3  OK
sampler-errors                       59   3   1  16   2   6  31    5    2  OK
sampler-predictors                   34   0   2  26   1   5   0    1    1  WEAK
sampler-state-format                 26   1   3  11   0   1   9    1    1  WEAK
sparse-factor                       137  20   9  39   5  14  38    4    5  OK
state-weight-pairing                 28   2   3  17   0   4   2    2    3  OK
summary-nondefault-families          17   0   7   2   0   8   0    3    2  OK
tree-store-order                     34   3   2   9   1  13   6    3    3  OK
bart-formula                          2   0   0   0   0   0   0    0    0  GAP
bcf-reporting                        46   2  11  22   4   7   0    0    3  OK
data-mixed-mutation                  15   0   2   7   0   3   2    0    2  OK
data-sparse                          44   0   5  16   1  10   7    2    1  WEAK
heteroscedastic                      29   0   3  13   2   4   5    1    0  GAP
plot-generics                        24   0   2   0   0  12  10    2    0  GAP
prior-init-composed-law              19   3   2   8   6   0   0    0    5  OK
sampler-trees                        32   0   5  14   7   2   4    1    2  OK
weighted-binary-ppd                  10   0   1   3   6   0   0    1    1  WEAK
xbart-weight                          5   0   1   0   0   2   2    2    0  GAP

OK 34, WEAK 19, GAP 12. No file went unreached: all 65 ran under all 85 mutations.

## Findings

BLOCKER = the mutation left all 167 tinytest files AND tests/cpp green, so the defect class
is invisible to the gates in shipped behavior.
1. BLOCKER  [[R/dbarts.R:1025@658869ac]]  `$growFromRoot(n.sweeps)` can ignore n.sweeps entirely
   (`.Call(C_dbarts_bartcore_growFromRoot, ptr, n.sweeps)` -> `, 1L`). 167 files green,
   tests/cpp green. test-grow-from-root.R spends its 18 assertions on finiteness, same-seed
   reproducibility and the R-level "positive integer" refusal - none of which move when the
   argument is discarded. A documented, user-facing parameter with zero behavioral coverage.
   SHOULD EXIST test-grow-from-root.R: at one seed, `expect_false(isTRUE(all.equal(
   growFromRootFit(1L), growFromRootFit(3L))))`.  agent-fix
2. BLOCKER  [[src/bartcore/combiner.hpp:906-913@658869ac]]  BCF `formForestVetoWeights` can drop its
   near-zero multiplier snap (`fabs(m) < zeroMultiplierTolerance ? 0.0 : w*m*m` -> `w*m*m`).
   167 files green, tests/cpp green. The snap is what makes the creation glue b0 = 0 leave
   every control row weightless in the treatment forest, so without it
   `$sampleTreesFromPrior()` may seat a treatment leaf holding only control rows - the exact
   invariant the compose slice (1e020abb) landed to hold. Narrowing the band to `m == 0.0`
   (d03) IS caught, but only by tests/cpp, never by R.
   SHOULD EXIST test-bcf-mutation-pins.R: after `$sampleTreesFromPrior()` on a b0 = 0 BCF,
   `expect_true(all(fit$getTrees(forest = 2L, current = TRUE, newdata = x[z == 1, ,
   drop = FALSE])$n[var == -1L] > 0L))`.  agent-fix
3. BLOCKER  [[src/bartcore/combiner.hpp:1822-1825@658869ac]]  multinomial `formForestVetoWeights` can
   ignore the active-row mask (`activeRows_.empty() ? omega[i] : omega[i]*activeRows_[i]`
   -> `omega[i]`). 167 files green, tests/cpp green. Same channel as finding 2, same
   silence: the mask reaches the sweep's precisions (c05 IS caught at
   [[test-active-rows-pins.R:470-471@658869ac]]) but not the prior tree draw's veto vector.
   SHOULD EXIST test-active-rows-pins.R: the multinomial mirror of finding 2's assertion,
   after `$setActiveRows(a)` then `$sampleTreesFromPrior()`.  agent-fix
4. MAJOR  [[src/bartcore/model.hpp:3415@658869ac]]  ordinal `updateCutpoints`'s loop bound
   (`s < numCategories_` -> `s + 1 < numCategories_`) freezes the last free cutpoint; at
   K = 3 that is every free cutpoint, so the ordinal cutpoint sampler is inert. All 167
   files green; only tests/cpp ([[test_model.cpp:4885@658869ac]]) kills it. test-ordinal.R's only value
   assertion on gamma_2 is [[test_model.cpp:356@658869ac]] `abs(mean(cutpoints[, 2L]) - 0.8) < 0.35`, and the frozen
   cold start (spacing 1.0) sits inside that band at |1 - 0.8| = 0.2.
   SHOULD EXIST test-ordinal.R after [[test_model.cpp:356@658869ac]]: `expect_true(sd(fitRec$cutpoints[, 2L]) > 0.02)`.
   agent-fix
5. MAJOR  [[src/bartcore/grow.hpp:186@658869ac]]  `if (growth <= 0.0) return;` -> `< 0.0` disables
   growTreeFromRoot's depth/availability veto. All 167 green; caught only by tests/cpp
   (test_grow.cpp "a growth-vetoed node draws nothing"). test-grow-from-root.R is 11/18
   smoke and never bounds a grown tree's depth.
   SHOULD EXIST test-grow-from-root.R: with a tree prior whose growth is 0 past depth 2,
   `expect_true(max(sampler$getTrees()$depth) <= 2L)` (or the equivalent walk).  agent-fix
6. MAJOR  [[R/data.R:1495@658869ac]], [[test-data-precision-warning.R:20-24@658869ac]], [[test-data-precision-warning.R:44-47@658869ac]]  the precision-degeneracy
   threshold is unpinned from below: `1e-10 -> 1e-30` leaves all 167 files green, because
   both warn-cases have `diff(range(y))` EXACTLY 0 (y.huge collapses to one double; y.const
   is constant), so only "range == 0" is tested and the constant the file exists to defend
   is never exercised. [[test-data-precision-warning.R:36@658869ac]] is the one tautology in the suite - it recomputes the code's own
   ratio and compares it to the code's own literal, testing nothing.
   SHOULD EXIST test-data-precision-warning.R: `y.near <- 1e15 + runif(n) * 10` (ratio ~1e-14,
   many distinct doubles) must warn "indistinguishable".  agent-fix
7. MAJOR  [[R/data.R:627-629@658869ac]]  the `family = "ordinal"` minimum-category refusal can be deleted
   with all 167 files green. Neither test-ordinal.R nor test-data-errors.R ever offers a
   one-level ordered response.
   SHOULD EXIST test-ordinal.R: `expect_error(bart2(x, ordered(rep("a", n)), family =
   "ordinal"), "at least 2 categories")`.  agent-fix
8. MAJOR  R/xbart.R (mcr arm, `sum(weights * misclassified) / sum(weights)`)  the
   misclassification loss can ignore case weights entirely with all 167 green. The rmse arm
   is covered (the same mutation there is killed by [[test-xbart-oracle.R:79-88@658869ac]]); mcr is not,
   and test-xbart-weight.R carries 5 assertions, none of them a weighted-loss oracle.
   SHOULD EXIST test-xbart-weight.R: an mcr weighted-loss oracle mirroring
   [[test-xbart-oracle.R:79-88@658869ac]], or at minimum `expect_false(isTRUE(all.equal(xbart(...,
   loss = "mcr", weights = w), xbart(..., loss = "mcr"))))` at a fixed seed.  agent-fix
9. MAJOR  [[R/rbart.R:270@658869ac]]  the "survival (aft) fits do not support 'weights'" refusal can be
   deleted with all 167 green. test-rbart-weights.R holds 3 assertions and never asks;
   test-rbart-aft.R does not either.
   SHOULD EXIST test-rbart-weights.R: `expect_error(rbart_vi(y ~ x, group.by = g, weights =
   w, family = "aft", ...), "do not support 'weights'")`.  agent-fix
10. MAJOR  [[R/dbarts.R:1929@658869ac]]  `getTrees`'s `current = TRUE` can be ignored (`useSaved <-
    control@keepTrees && !current` -> `control@keepTrees`), reading the saved store where the
    caller asked for the live trees; all 167 green. This is load-bearing:
    test-empty-leaf-veto-weights.R's entire structural oracle is getTrees(current = TRUE),
    and test-sampler-trees.R has 32 assertions over the same surface.
    SHOULD EXIST test-sampler-trees.R: after `run()` on a keepTrees sampler,
    `expect_false(identical(sampler$getTrees(current = TRUE), sampler$getTrees()))`.
    agent-fix
11. MAJOR  [[R/generics.R:2364-2369@658869ac]]  print.bart's whole synopsis body can print wrong values
    (family always "gaussian"; n.chains under the n.trees label) with all 167 green.
    [[test-plot-generics.R:103-104@658869ac]] are `expect_true(is.character(capture.output(print(x))))`,
    which CANNOT FAIL - capture.output always returns character. The file is 24 assertions:
    12 expect_silent, 10 refusals, 2 that cannot fail, zero content pins.
    SHOULD EXIST [[test-plot-generics.R:103@658869ac]]: `expect_true(any(grepl("family: probit",
    capture.output(print(fit.bin)), fixed = TRUE)))`.  agent-fix
12. MINOR  [[R/data.R:735@658869ac]]  validateForestBases's finiteness refusal is double-guarded at
    creation ([[A_class.R:624@658869ac]] refuses first, which is what [[test-bcf-creation.R:352-355@658869ac]]
    actually exercises) and untested on the predict path ([[generics.R:647@658869ac]]), so deleting it
    leaves all 167 green while `predict(fit, newdata, bases = <NaN>)` would reach the blend.
    SHOULD EXIST test-predict-blend.R: `expect_error(predict(fit, xNew, bases = list(NULL,
    cbind(rep(NaN, m), 0))), "must all be finite")`.  agent-fix
13. MINOR  [[R/data.R:381@658869ac]]  the "cannot find test offset" message text is unpinned (replacing it
    leaves all 167 green), while the sibling text at [[R/data.R:231@658869ac]] IS pinned three times
    ([[test-data-errors.R:90@658869ac]], [[test-sampler-errors.R:43@658869ac]], [[test-sampler-trees.R:274@658869ac]]).
    test-error-quality.R is 8 assertions and does not cover it.
    SHOULD EXIST test-error-quality.R: `expect_error(dbartsData(y ~ x, df, test.offset =
    nosuch), "cannot find test offset")`.  agent-fix
14. MINOR  [[src/bartcore/sampler.hpp:485@658869ac]]  `setCurrentSampleNum` is an equivalent mutant -
    discarding its argument changes nothing anywhere, because it has no caller in src/, R/,
    tests/ or inst/. Dead public accessor on the engine's shape surface, not a test gap.
    VD-judgement (delete, or wire it into the setState restore that sets
    state.currentSampleNum directly at [[R_interface_bartcore.cpp:6738@658869ac]]).
15. MINOR  test-bart-formula.R (2 assertions, both `expect_inherits(..., "bart")`) is
    class-only smoke: no mutation, aimed or incidental, can move it, and none did.
    VD-judgement (fold into test-data-formula.R, or give it one content pin).
16. MINOR  b04 (logistic PG latent ignoring the case weight, [[model.hpp:3554@658869ac]]) makes the suite
    HANG rather than fail - the 65-file run did not finish in 300s, where every other
    mutation finished in under 45s. The only signal a CI job gets is a timeout. Worth a
    separate look as robustness (an unbounded loop reachable from an inflated working
    response), not as a test gap.  VD-judgement
17. MINOR  Seven reach gaps confined to the CHANGED set - an untouched file is the SOLE
    killer, so no new assertion is needed, only a note not to trim that file: r12
    ppdNoiseScale's per-row weight ([[test-generics-posteriorPredictiveDistribution.R:90@658869ac]],
    [[test-ppd-sigma-pairing.R:152@658869ac]]); r22 extract(type="forest", contribution=TRUE) with `*`
    made `+` ([[test-bcf-forest-channel.R:90@658869ac]]); b06 the interaction subtree walk
    ([[test-interactions.R:96@658869ac]]); b07 bart2's grow-from-root falling back to the prior draw
    ([[test-bart2-grow-from-root.R:107@658869ac]], [[test-grow-from-root-categorical.R:75@658869ac]]); b15 rbart_vi's
    group.by type refusal ([[test-rbart-error.R:26@658869ac]]); d06 xbart's rmse weights
    ([[test-xbart-oracle.R:79@658869ac]]); d03 the BCF snap band (tests/cpp only).  defer

## Real defects exposed incidentally

None. No assertion failed on the pristine tree in any run, both independent pristine builds
reproduced the 167-file suite at 0/0, and every red assertion traced to its planted mutation.
The one behavior worth flagging is not a defect but finding 16's hang.

## Spot check: do the ten untouched files hold?

NO - the leg must widen. 6 of 10 came back OK/WEAK (bcf-reporting OK, data-mixed-mutation OK,
prior-init-composed-law OK, sampler-trees OK, data-sparse WEAK, weighted-binary-ppd WEAK) and
4 came back GAP. Two of the four failed an AIMED probe: plot-generics survived both print-
synopsis mutations (finding 11), xbart-weight survived both loss-weight mutations (finding 8).
bart-formula is structurally unfalsifiable (finding 15). heteroscedastic (29 assertions) took
only an incidental probe and did not move, so it is unproven rather than proven bad.

Rate, stated honestly: 8 of 55 in-scope files are GAP (14.5%) against 4 of 10 untouched
(40%), but in-scope files took ~1.5 aimed mutations each and spot files ~0.6, so the true gap
is at most a factor of two, not four. Either way the wave-3 licence - that the untouched
suite needs no re-examination - does not survive: the untouched population carries the same
reach class at no lower a rate. Recommend a second leg over the 112 untouched files, aimed
first at the print/plot/summary surface and the xbart loss arms.

## Not reached

Nothing in the 65-file set went unrun. Left for a later leg: inst/tinytest/capi/consumer.c as
a mutation TARGET (it served as a killer throughout - test-capi.R moved under 8 mutations -
but no mutation was aimed at the flat header itself); and per-file depth beyond the ~1.5
aimed mutations the budget bought.
