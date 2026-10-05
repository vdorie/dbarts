# empty-leaf-membership: a leaf is empty only when no row reaches it

Status: PLANNED.

agent: opus implementer, one; opus reviewer.
rng: POSTERIOR-CHANGING for a fit with a zero case weight, a zero per-forest weight, an active-row mask or a
zero-trial multinomial row, and for a multi-forest fit whose basis and amplitudes give some rows zero weight
in a forest although no weight is installed - a treatment forest whose amplitudes are held at (b0, b1) =
(0, 1), where every control row is weightless for the whole run: the set of trees the sampler may hold
changes, and a latent family redraws the latents of rows a mask switches back in. DRAW-SHIFTING in principle
for any other two-forest fit: its amplitudes start at (1, 0, 1), so its treatment forest is in the same
position for the prior draw and the first sweep, until b0 is first drawn. NEUTRAL for every other fit. Shown
bit for bit on the shipped build: 52 of 55 scenarios of the gaussian harness (all but zeroweights, maskprobit
and maskordinal), 13 of 15 of the BCF harness (all but masked and glue_toggle; the thirteen include every
drawn-amplitude scenario that installs no mask, none of which met the first-sweep case), and 11 of 11
multinomial; the four seeded snapshot files pass unchanged on the reference build.
window: pre-release (dec-B238).
budget: ~800 lines (C++ ~150, tests/cpp ~200, tinytest ~150, a tracked exact harness ~200, design notes, manual
and records ~100). Plans have run 1.5-2x low.

## Goal

The rule that rejects a tree move leaving a leaf empty judges emptiness by membership, as 0.9-34 did: a leaf
is empty only if no row at all reaches it. A leaf that holds only rows of zero weight, or only rows the mask
switches off, is legal; it contributes nothing to the likelihood and its value is drawn from the prior. The
prior over trees then does not depend on the weights or the mask, so a larger sampler that redraws the mask
every sweep samples the model it assumes.

## Context

- Today a leaf with members but no positive-weight member loses outright to any branch a likelihood term
  reaches (retired: [`Tree::leafVetoRank`](../../src/bartcore/tree.hpp), rank 1, which step 1
  removes; the cut scan's sentinels in scan.hpp read the weight on each side). That made a fixed zero weight the same as deleting the row, and made
  the set of allowed trees depend on the mask.
- Measured against the exact posterior of a two-part mixture whose membership is redrawn every sweep, ten
  rows, one and two trees: under today's rule the long-run membership probabilities are off by up to 0.066
  (about 320 standard errors) and the fitted means by a quarter of a posterior standard deviation, the trees
  too small; judging by membership, every quantity is within Monte Carlo error. At 300 rows and 50 trees the
  two rules differ by up to 0.044 in a membership probability, where a region is mostly switched off, and by
  under 0.005 in 290 of 300 rows.
- With no zero weight and no mask the two rules are the same numbers and the draws are identical, in a
  single-forest fit. Found in implementation: a multi-forest sampler composes each forest's precisions from
  its amplitudes, and a multiplier of zero is a zero weight in that forest that nobody installed. A
  two-forest fit starts at amplitudes (1, 0, 1), so every control row is weightless in the treatment forest
  until b0 is first drawn, and for good when (b0, b1) is held at (0, 1). The old rule conditioned that
  forest's prior draw and its moves on holding a treated row in every leaf; the new one does not. The BCF
  equivalence scenario with held amplitudes (glue_toggle) leaves its recorded stream; over 20 seeds no
  posterior summary of it moves beyond a standard error, the change being 0.05 treatment-forest leaves per
  sweep that hold control rows only. The drawn-amplitude scenarios were bitwise.
- The same holds, measured the same way, with an unordered factor (the old rule's worst case: membership off
  by 0.21, the new within noise), with zero case weights in place of the mask, with a mask that empties whole
  regions, and after grow-from-root under a mask.
- A probit fit needs one thing more. A row switched back in keeps the latent it had when it was switched off,
  which was drawn against an older fit, and with the mask redrawn every sweep that alone leaves the combined
  sampler off (membership by 0.004, the fit by 0.05 of a posterior standard deviation, its spread 5 to 7
  percent too narrow). Redrawing the latents after the mask is installed makes it exact. A logistic sampler
  already redraws its latents when its counts are swapped.
- What a fixed zero weight now costs against deleting the row, 300 rows and 50 trees, split points and scales
  made equal: at active rows the fitted mean moves by at most 0.02 of a posterior standard deviation with a
  random half masked and 0.05 with a region masked; the residual scale does not move; mixing is not worse.
- A leaf's draw at zero total weight is already a draw from its prior, and its integrated likelihood already
  0 on the log scale, for the conjugate leaves.
- What does not change: the outright rejection of a leaf no row reaches (dec-A12); the residual variance's
  degrees of freedom, which count positive-weight rows; every sufficient statistic; the merge of member-empty
  leaves after a data or cut change; a masked row's NaN pointwise log-likelihood.

## Constraints

- A single-forest fit that installs no zero weight, no mask and no zero-trial row draws exactly what it
  draws now: the seeded snapshot files and every equivalence scenario of the gaussian and multinomial
  harnesses without one are unchanged. A multi-forest fit does too unless a forest's own multiplier is zero
  at some rows (Context): held there, the fit is posterior-changing; zero only at creation, its draws can
  shift from the first sweep on with the posterior unchanged.
- One rule on every path that decides whether a branch is legal: the moves, the cut scans (ordinal and
  categorical), grow-from-root, the per-forest weight composition of a multi-forest sampler, and the variance
  forest.
- A leaf holding only switched-off rows must behave on every leaf model: constant, linear, gp, the monotone
  leaves and the variance leaf. Where a leaf model's integrated likelihood or draw at zero weight is not the
  prior's today, make it so or say why it cannot be.
- An all-zeros mask still runs, every forest at its prior.
- On every family that carries a per-row latent (probit, ordinal, logistic, negative binomial, aft's censored
  times, Student-t's scales, multinomial), a row the mask switches back in has its latent redrawn from its
  conditional given the current fit before the call returns, from the sampler's own generators, as a logistic
  count swap already does. Rows that stay active keep their latents, and a mask that reactivates no row draws
  nothing, so a sampler whose mask does not change consumes no extra variates. Found in implementation:
  multinomial holds no latent between sweeps - each category's Polya-Gamma column is drawn inside the sweep,
  immediately before that forest reads it - so it has nothing to redraw and redraws nothing.
- No NEWS entry: this restores what 0.9-34 counted, and the mask is new in 1.0-0.

## Steps

1. The rule: emptiness by member count in retired: [`Tree::leafVetoRank`](../../src/bartcore/tree.hpp), which
   [`Tree::leafIsEmpty`](../../src/bartcore/tree.hpp) replaces, and the scan
   sentinels, and wherever else a zero weight sum is read as "empty". Remove what the middle rank needed and
   nothing else uses. tests/cpp: a move that leaves a leaf of only zero-weight rows is accepted on its
   likelihood, on each leaf model; such a leaf's value is a prior draw; a leaf no row reaches is still
   rejected; the scans agree with the moves on which candidate is legal; grow-from-root under a mask builds no
   member-empty leaf.
2. Reactivated latents: the redraw above, per family. tests/cpp: a reactivated row's latent is drawn against
   the current fit, an always-active row's is untouched, and a mask equal to the one in force consumes no
   variate. The exact harness of step 4 runs a probit arm that fails without it.
3. tinytest: the tests that pin the old rule are rewritten to the new one, each keeping what it protects (an
   all-zeros mask runs; a forest grown into rows the mask then removes keeps moving; zero-trial rows; binary
   0/1 weights; per-forest weights). A fit with zero weights no longer has to equal the fit on the remaining
   rows alone: where a test asserted that, it asserts instead what still holds (the sufficient statistics and
   the residual degrees of freedom) and that a leaf of only zero-weight rows can occur.
4. A tracked exact harness in benchmarks/R: the small mixture whose membership is redrawn every sweep, its
   exact posterior by enumeration, and the combined sampler's long-run membership probabilities and fitted
   means against it, with a quick mode sized for CI. It fails on the old rule. Add it to the exact gates.
5. The exact gates and balance gates that install zero weights or a mask are re-read: each either does not
   move or has its expectation restated under the new rule, with the reason.
6. The compare against 0.9-34 is not re-recorded. Re-run its zero-weight row and restate what the list of
   explained differences says about empty leaves: 1.0-0 again counts a zero-weight row as a member of its
   leaf, as 0.9-34 did, and still rejects outright where 0.9-34 charged a finite penalty.
7. Design notes and manual: the empty-leaf note's section on what counts as empty, the active-row mask note
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
