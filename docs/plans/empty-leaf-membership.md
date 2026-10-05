# empty-leaf-membership: a leaf is empty only when no row reaches it

Status: PLANNED.

agent: opus implementer, one; opus reviewer.
rng: POSTERIOR-CHANGING for a fit with a zero case weight, a zero per-forest weight, an active-row mask or a
zero-trial multinomial row: the set of trees the sampler may hold changes. NEUTRAL for every other fit, whose
draws are bit for bit unchanged.
window: pre-release (dec-B238).
budget: ~600 lines (C++ ~60, tests/cpp ~150, tinytest ~150, a tracked exact harness ~150, design notes, manual
and records ~90). Plans have run 1.5-2x low.

## Goal

The rule that rejects a tree move leaving a leaf empty judges emptiness by membership, as 0.9-34 did: a leaf
is empty only if no row at all reaches it. A leaf that holds only rows of zero weight, or only rows the mask
switches off, is legal; it contributes nothing to the likelihood and its value is drawn from the prior. The
prior over trees then does not depend on the weights or the mask, so a larger sampler that redraws the mask
every sweep samples the model it assumes.

## Context

- Today a leaf with members but no positive-weight member loses outright to any branch a likelihood term
  reaches ([`Tree::leafVetoRank`](../../src/bartcore/tree.hpp), rank 1; the cut scan's sentinels in
  scan.hpp read the weight on each side). That made a fixed zero weight the same as deleting the row, and made
  the set of allowed trees depend on the mask.
- Measured against the exact posterior of a two-part mixture whose membership is redrawn every sweep, ten
  rows, one and two trees: under today's rule the long-run membership probabilities are off by up to 0.066
  (about 320 standard errors) and the fitted means by a quarter of a posterior standard deviation, the trees
  too small; judging by membership, every quantity is within Monte Carlo error. At 300 rows and 50 trees the
  two rules differ by up to 0.044 in a membership probability, where a region is mostly switched off, and by
  under 0.005 in 290 of 300 rows.
- With no zero weight and no mask the two rules are the same numbers and the draws are identical.
- A leaf's draw at zero total weight is already a draw from its prior, and its integrated likelihood already
  0 on the log scale, for the conjugate leaves.
- What does not change: the outright rejection of a leaf no row reaches (dec-A12); the residual variance's
  degrees of freedom, which count positive-weight rows; every sufficient statistic; the merge of member-empty
  leaves after a data or cut change; a masked row's NaN pointwise log-likelihood.

## Constraints

- A fit that installs no zero weight, no per-forest zero weight, no mask and no zero-trial row draws exactly
  what it draws now: the seeded snapshot files and every equivalence scenario without one are unchanged.
- One rule on every path that decides whether a branch is legal: the moves, the cut scans (ordinal and
  categorical), grow-from-root, the per-forest weight composition of a multi-forest sampler, and the variance
  forest.
- A leaf holding only switched-off rows must behave on every leaf model: constant, linear, gp, the monotone
  leaves and the variance leaf. Where a leaf model's integrated likelihood or draw at zero weight is not the
  prior's today, make it so or say why it cannot be.
- An all-zeros mask still runs, every forest at its prior.
- No NEWS entry: this restores what 0.9-34 counted, and the mask is new in 1.0-0.

## Steps

1. The rule: emptiness by member count in [`Tree::leafVetoRank`](../../src/bartcore/tree.hpp) and the scan
   sentinels, and wherever else a zero weight sum is read as "empty". Remove what the middle rank needed and
   nothing else uses. tests/cpp: a move that leaves a leaf of only zero-weight rows is accepted on its
   likelihood, on each leaf model; such a leaf's value is a prior draw; a leaf no row reaches is still
   rejected; the scans agree with the moves on which candidate is legal; grow-from-root under a mask builds no
   member-empty leaf.
2. tinytest: the tests that pin the old rule are rewritten to the new one, each keeping what it protects (an
   all-zeros mask runs; a forest grown into rows the mask then removes keeps moving; zero-trial rows; binary
   0/1 weights; per-forest weights). A fit with zero weights no longer has to equal the fit on the remaining
   rows alone: where a test asserted that, it asserts instead what still holds (the sufficient statistics and
   the residual degrees of freedom) and that a leaf of only zero-weight rows can occur.
3. A tracked exact harness in benchmarks/R: the small mixture whose membership is redrawn every sweep, its
   exact posterior by enumeration, and the combined sampler's long-run membership probabilities and fitted
   means against it, with a quick mode sized for CI. It fails on the old rule. Add it to the exact gates.
4. The exact gates and balance gates that install zero weights or a mask are re-read: each either does not
   move or has its expectation restated under the new rule, with the reason.
5. The compare against 0.9-34 is not re-recorded. Re-run its zero-weight row and restate what the list of
   explained differences says about empty leaves: 1.0-0 again counts a zero-weight row as a member of its
   leaf, as 0.9-34 did, and still rejects outright where 0.9-34 charged a finite penalty.
6. Design notes and manual: the empty-leaf note's section on what counts as empty, the active-row mask note
   and the mask's seven rules, the grow-from-root note; the sampler help for `setWeights`, `setActiveRows` and
   `setForestWeights`, and `bart`'s `weights`: a zero-weight or masked row is in the design but not in the
   likelihood, so a fit with one is not the fit on the remaining rows alone, and a larger sampler that redraws
   the mask or the weights every sweep samples the joint model it assumes.

## Verification

- `R CMD INSTALL --preclean` into the slice's own library; full tinytest; `tests/cpp` builds and passes, clean
  under ASan and UBSan.
- The four seeded snapshot files pass unchanged on a reference build.
- The equivalence compares in statistical (z) mode against the current baselines: every scenario without a
  zero weight or a mask identical; any scenario with one is re-recorded as the manifest's earlier
  posterior-changing re-records were, its z statistics against the previous baseline reported.
- Every exact gate in `.github/workflows/exact-gates.yaml` in quick mode, the new one included.
- `lintr::lint_package()`, `air format --check .`, the three tools/ checks and
  `Rscript benchmarks/R/mutation-battery.R verify-anchors` clean.
