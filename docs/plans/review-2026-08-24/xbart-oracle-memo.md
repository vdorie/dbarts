# ORACLE MEMO - what an xbart() oracle is, and what is already covered

Read-only, tip b102e17c. Every claim is tool-verified; probe numbers come from a staged
`git archive HEAD` + `R CMD INSTALL --preclean` build (dbarts 1.0.0).

FIRST CORRECTION. A prior note, whose location could not be placed at any candidate
commit (unresolved: [[review-lenses-memo.md:142@b102e17c]] - no file of that name ever existed in this
repository's history), says "xbart still has zero equivalence
scenarios, no exact gate, no SBC arm". The first clause is FALSE at the tip: the
canonical baseline equivalence-5a3bc276.rds carries 43 scenarios and `xbart` is one of
them (verified by reading the .rds). The phrase "no oracle at any date" appears nowhere
in any scratch note or docs/. The accurate statement is narrower, and is section 1.

## 1. Decomposition and coverage

xbart's answer = split assembly x per-cell fit x loss x aggregation x chunking.

| piece | what it is | oracle today |
|---|---|---|
| split assembly | [[R/xbart.R:718@5a3bc276]] `sample.int(n)` per rep, cut into `foldSizes` ([[R/xbart.R:353@5a3bc276]]); train = complement | PARTIAL. test-xbart-oracle.R pins fold SIZES arithmetically (12 -> 4,4,4 and 3,3,2,2,2 -> 2.4) and pins a partition by VALUE (`sort(unlist(fold y)) == sort(y)`). WHICH rows, and that train is the complement: nothing. |
| per-cell fit | `bartcoreSamplerFromHandle` row-subset view ([[R/bartcore.R:689@5a3bc276]]), fresh per n.trees, warm across cells at n.burn[2] | PARTIAL. The engine's own gates cover the sampler. test-data-handle.R pins a FULL-row view bitwise == an ordinary sampler, and for a proper fold only shapes, finiteness, and one offset mean-within-1. Test-row ALIGNMENT and no-leakage: nothing. Warm-start-across-cells: nothing anywhere. |
| loss | rmse / weighted rmse / log / mcr / custom ([[R/xbart.R:564-611@5a3bc276]]) | FULL. f009eff8's test-xbart-oracle.R, 16 assertions: a capturing loss records the (y.test, testSamples, weights) triples and every built-in formula is transcribed BY HAND outside xbart. Mutation-proven twice (fold-count off-by-one kills 6/16; swapped log branches + reversed rmse rows kills 3/16). |
| aggregation | fold mean, then placement into (rep, n.trees, k, power, base, loss) by `linearIndex` ([[R/xbart.R:509-517@658869ac]]) | PARTIAL. The fold mean is pinned (a loss returning `length(y.test)` reports 4 and 2.4). Placement is SHAPE-ONLY: dims and dimnames are asserted in test-xbart-method.R / -loss.R; no test ties the slice LABELLED k = X to a fit that used k = X. Two equal-length axes could be transposed silently. |
| chunking | `numChunks = min(n.threads, n.reps)`, one seed per chunk ([[R/xbart.R:452-461@658869ac]]) | PARTIAL. Fixed (seed, n.threads) reproducibility is pinned; the DIFFERENCE across thread counts is pinned as an inequality (`any(xval.1 != xval.3)`). Nothing pins what threading must preserve. |
| drift channel | equivalence scenario `xbart`, recorded channel = the loss array | NOT AN ORACLE. It is a self-comparison: the recording and the re-run are the same code. It catches movement, not error. |

f009eff8's oracle proves exactly one thing: GIVEN the triples xbart hands its loss, the
arithmetic from there to the reported number is right. Every defect UPSTREAM of the loss
call - wrong test rows, a train set overlapping the test set, a cell landing in the wrong
array slot, a chunk clobbering shared state - passes all 16 assertions. That is the gap.

## 2. Candidate oracles for the uncovered pieces

(a) FOLD-ASSEMBLY. Feasible, and cheaper than the brief assumed - see probe A1.
    xbart does not expose fold membership, but it is RECONSTRUCTIBLE outside: with a
    seed and n.threads = 1, `set.seed(seed); s <- sample.int(.Machine$integer.max, 1)`
    then `set.seed(s); sample.int(n)` reproduces replication 1's permutation exactly,
    and a capturing loss reveals the rows it produced. Three arms:
      a1 exact rep-1 fold rows + disjointness + train == complement.
      a2 row alignment: per-fold cor(y.test, rowMeans(testSamples)) against a
         permutation null - catches a misordered gather in the view's test channel.
      a3 leakage: y independent of x; held-out rmse must sit at or above sd(y) while an
         in-sample fit at the same hyperparameters sits far below.
    Proves: the split is the split the design claims, the fits score the rows they held
    out, and they never trained on them. Cost: one tinytest file, R only, ~90 lines,
    seconds of runtime, no engine build. Sonnet.
    NOT feasible as the brief framed it: reproducing the reported loss bitwise from k
    separate bart() fits. A fold sampler is a VIEW over the full data's cut grid, while
    bart(x[trainRows, ]) re-quantizes on the subset's own range; and replications past
    the first draw their permutation after the engine has consumed R's stream (probe A1:
    rep 2 reconstructible from R alone = FALSE). A handle-level replay would reproduce
    it, but it re-executes xbart's own body - a change-detector, not an oracle.

(b) KNOWN-ANSWER / AXIS PLACEMENT. Two forms.
    b1 extreme cells: k = 20 shrinks the sum-of-trees to the training mean, so that
       slice's rmse must land at sd(y.test); n.trees = 1 must beat nothing. Deterministic
       enough to assert with a wide band, and it ties each labelled slice to the
       hyperparameter it claims. Cost: one arm, ~25 lines, ~20 s. Sonnet.
    b2 rank-the-truth (a nested design whose best cell is known, asserted to rank first
       with high probability): needs many reps for power and buys little over b1, since
       the failure mode here is placement, not statistics. Skip.

(c) THREADING. A 1-thread == n-thread pin CANNOT exist; P16 measured this and R/xbart.R
    documents it (seeds are per CHUNK, chunks change with n.threads). The achievable
    invariant, currently unpinned, is FIRST-CHUNK PREFIX IDENTITY: chunk 1's seed is
    `sample.int(m, numChunks)[1]`, which does not depend on numChunks, so the
    replications in chunk 1 must be bitwise identical at every thread count. That is the
    statement that the parallel path only re-seeds and corrupts nothing. Probe D confirms
    it holds. Cost: ~15 lines added to test-xbart-reproducibility.R beside the existing
    inequality. Sonnet.

(d) 0.9-34 CROSS-VERSION. `git show main:R/xbart.R` is 144 lines handing the whole loop to
    `.Call(C_dbarts_xbart)`, and src/crossvalidate.cpp holds an independent C++ fold
    permutation, fold assembly, and loss functors - a genuinely independent implementation
    of the same target, the strongest oracle available in principle. Statistical only, and
    the defaults have diverged (n.burn length 3 vs 2, node.scale, sigma vs sigest, the
    k-sort in (e)), so each must be matched by hand. Cost: two installs, a benchmarks/R
    harness, a quiet machine, an Opus adjudication. High. Defer unless (a)-(c) fire.

(e) SIDE FINDING, not an oracle. 0.9-34 sorts the k grid DECREASING before sweeping and
    un-permutes the result (`kOrder`/`kOrder.inv`), so its answer is invariant to the
    order a user lists k in. The rewrite drops that sort - cells sweep in `expand.grid`
    order, warm-starting off each other - and nothing in docs/, inst/NEWS.Rd, or the
    landing notes records the removal (grepped). Probe C measures the consequence.

## 3. Live probe (staged build, foreground)

A1 fold assembly, n = 24, 4 folds, seed 7, 1 thread. Hand-reconstructed rep-1 folds
   `2,5,7,16,18,19 | 3,4,6,11,12,13 | 8,9,10,14,17,23 | 1,15,20,21,22,24` are IDENTICAL
   to the rows the capturing loss saw; the four folds cover 1:24 exactly. Rep 2
   reconstructible from R's stream alone: FALSE (as expected - engine draws intervene).
A2 row alignment, same run: per-fold cor(y.test, rowMeans(testSamples)) =
   0.749 0.888 0.969 0.813 0.985 0.943 0.693 0.946; permutation null over 200 shuffles =
   -0.033 .. +0.034. Aligned, and a shuffle is detectable.
A3 leakage, n = 60, y ~ N(0,1) independent of x, k = 0.5, 50 trees: sd(y) = 0.9639,
   xbart held-out rmse = 1.2751 / 1.4449, in-sample rmse at the same hyperparameters =
   0.4306. No leakage; the discrimination is 3x.
B  axis placement, Friedman n = 100, sd(y) = 5.219, mean over 2 reps:
   k=0.5: n.trees 1 -> 4.099, 50 -> 2.752;  k=20: 1 -> 4.673, 50 -> 5.024.
   The k = 20 slice sits at sd(y) as the shrink-to-mean argument requires. Placement is
   correct today and is exactly what b1 would pin.
C  order dependence, same data, n.reps = 1 so both arms share the rep-1 folds:
   k = c(2, 8) -> 2.1346, 3.3444;  k = c(8, 2) -> 3.2739, 1.9543. The k = 2 cell moves
   0.1804 (8.4%) purely by where it sits in the sweep, k = 8 moves 0.0705. One draw, so
   this does not separate warm-start bias from stream noise - but it is the quantity
   0.9-34's k-sort existed to remove, and it is unmeasured on this branch.
D  threading, n.reps = 4, seed 0: 1 vs 2 threads, per-rep bitwise identity = T T F F;
   1 vs 4 = T F F F. Exactly the chunk-prefix structure. The invariant in (c) holds and
   is free to assert.

## 4. Recommendation

BUILD (a1+a2+a3) + (b1) + (c) as one pre-RC slice: R-only, no engine touch, one tinytest
file plus ~15 lines in test-xbart-reproducibility.R, ~1 Sonnet session with an Opus read
of the assertions. DEFER (d) - the only oracle that could catch a systematically-wrong
loop, but it costs a second install plus a hand-matched default table, and (a)-(c) cover
the same failure classes far cheaper. SEPARATELY, put (e) to VD as a surface question.

THE ONE FACT THAT DECIDES IT: f009eff8's oracle transcribes whatever triple xbart hands
its loss, so it is blind to every upstream defect - wrong rows, leakage, transposed axis,
chunk corruption - and probe A1 shows that upstream is exactly reconstructible from
outside at a fixed seed for replication 1. The gap is real and the fix is cheap.
