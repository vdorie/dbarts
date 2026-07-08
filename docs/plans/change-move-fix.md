# change-move-fix

agent: opus (both stages)
rng: posterior-changing (the whole point: the current chain targets
     the wrong distribution; both candidate fixes change draws AND
     the stationary distribution to the correct one)
budget: stage 1 instrumentation-only (falsifier pattern, nothing
        lands but a driver + numbers); stage 2 dual prototype in a
        worktree behind a switch, one variant lands on VD's pick

## Goal

The change move satisfies detailed balance (correctness-audit block
2, confirmed z = +479 by benchmarks/R/change-balance.R), chosen
between two exact variants by measured mixing cost under typical
usage - VD's criterion: users will not know to lengthen chains to
compensate for silently rejected proposals.

## The variants (both exact)

(a) Propose-from-prior: draw the new rule from the prior (full
    ancestor interval / all 2^R - 2 patterns); a draw whose
    descendants become unsatisfiable is an automatic rejection
    (pi = 0), restoring the birth/death cancellation with no
    counting. Also removes the categorical 64-attempt asymmetry.
    Cost: no-op proposals at nodes with grandchildren.
(b) Hybrid: ordinals keep the restricted good-set proposal and the
    acceptance gains the counted |Valid_old|/|Valid_new| ratio
    (menus already enumerated by findGoodOrdinalRules);
    categoricals use (a) (counting descendant-valid patterns is
    exponential for wide masks - infeasible). Cost: two mechanisms.
VD leans (b) on the silent-mixing-cost argument; stage 1 can
discharge or confirm that concern before stage 2 decides.

## Steps

1. Instrumentation (no fix): on stock fits, log per change proposal
   the valid-menu fraction |Valid|/|SI| at the target node (ordinal:
   both already computed; categorical: the rejection loop's attempt
   count is a geometric draw of the fraction). Report the implied
   (a) no-op rate by config: Friedman n in {1e3, 1e4} x m in
   {75, 200} (defaults); a mixed continuous+categorical design; a
   deep single-tree config; a DART run. Instrumented engine reverts
   per the falsifier pattern; only the driver and numbers land.
2. Dual prototype in a worktree behind an internal switch: variants
   (a) and (b) plus the current (defective) kernel for reference.
   Gate first: BOTH variants pass benchmarks/R/change-balance.R
   (engine matches the exact arm; the current kernel's failure
   reproduces). Then ESS-per-second and per-sweep on sigma, k,
   train fits, and variable-inclusion traces, same config grid as
   step 1. One table.
3. VD picks; the winner lands as its own gated item
   (posterior-changing: snapshot regeneration, equivalence anchor
   re-record, all exact-posterior gates incl. change-balance.R as a
   permanent gate, design note, bench-sampler on a quiet machine).
   NEWS entry: seeded results change; variable-selection behavior
   at unequal cut counts corrected.

## Constraints

- Nothing lands from stage 2 until VD picks.
- change-balance.R is the correctness arbiter for both variants
  before any performance comparison.
- Out of scope: replacing change with different move types
  (rotate etc.) - a separate research item if ever.

## Verification

- Stage-1 table + stage-2 table recorded here; the landed variant's
  gates run under its own item.

## Stage 1 results

Driver: benchmarks/R/change-fix-instrumentation.R (raw numbers in
benchmarks/results/change-fix-instrumentation.csv). The engine's
changeMove was instrumented log-only (moves.hpp changefix_instr,
reverted before commit): per change proposal it recorded whether the
target node's two children are both leaves (Valid = SI exactly, zero
waste), and, for the deeper nodes, the variant-(a) no-op probability
1 - |Valid(v)|/|SI(v)| - ordinals exact from findGoodOrdinalRules good
count vs splitInterval size (the missing-direction coin doubles both
and cancels), categoricals from the rejection loop's attempt count
(1 - found/attempts, a truncated-geometric with mean 1/valid-fraction).
Move-type counts came from a per-sweep counter. Reads only values the
move already computed; the ext_rng stream is untouched. 500 burn +
1000 logged sweeps per config, single chain/thread, seed 20260708.

Move mix is stable across all configs: birth ~0.26, death ~0.24, swap
~0.10, change ~0.39 of moves (change is what this item touches).

Per config (change = fraction of moves that are change proposals;
leaf-share = change proposals whose children are both leaves, exactly
zero-waste; deep no-op = 1 - |Valid|/|SI| over the remaining, deeper
proposals; overall = mean no-op over ALL change proposals incl. the
zero-waste leaves; /sweep = expected no-op change proposals per sweep):

  config                         change  leaf   deep no-op          overall  /sweep
                                  /moves  share  mean  med   p90    /prop
  Friedman  n=1e3  m=75           0.375  0.784  0.051 0.00  0.100   0.0110   0.311
  Friedman  n=1e3  m=200          0.368  0.808  0.056 0.00  0.100   0.0108   0.795
  Friedman  n=1e4  m=75           0.391  0.658  0.055 0.00  0.150   0.0190   0.556
  Friedman  n=1e4  m=200          0.365  0.822  0.053 0.00  0.090   0.0095   0.692
  Mixed     n=1e4  m=75           0.394  0.729  0.066 0.00  0.220   0.0179   0.530
  Deep 1-tree n=1e3 power=1       0.384  0.365  0.183 0.00  0.737   0.1164   0.045
  Mixed+DART n=1e4  m=75          0.393  0.727  0.151 0.00  0.699   0.0410   1.211

Mixed configs also carry a categorical share (0.38 non-DART, 0.18 with
DART) whose pooled draw-level no-op (1 - successes/draws over all
categorical proposals) is 0.140 (Mixed) and 0.235 (Mixed+DART) -
folded into the numbers above.

Reading:

- Leaf-children traffic dominates default usage: 66-82% of change
  proposals hit nodes whose children are both leaves, where variant (a)
  wastes nothing. The exception is the deep single-tree stress (36%):
  deep trees put most change traffic at internal nodes.
- The deep-node no-op median is 0 in every config - most deeper nodes
  still have no same-variable descendant split, so even there variant
  (a) usually costs nothing. The waste lives in a heavy tail: p90 rises
  from 0.09-0.15 (Friedman) to 0.22 (Mixed), 0.70 (Mixed+DART), 0.74
  (deep tree).
- Overall variant-(a) no-op rate per change proposal is small under
  default shallow forests: ~1-2% (Friedman, Mixed), i.e. 0.3-0.8
  wasted change proposals per sweep. It climbs under the regimes the
  defect actually distorts: 4.1% (Mixed+DART, 1.2 no-ops/sweep) and
  11.6% per proposal for deep single trees (only 0.045/sweep there
  because ntree=1 yields <1 change proposal per sweep).
- DART roughly doubles the waste of the same mixed design (per-proposal
  0.018 -> 0.041; deep p90 0.22 -> 0.70; categorical pooled 0.14 ->
  0.24): concentrating split-variable draws onto a few variables raises
  same-variable descendant stacking, the source of every no-op.

Verdict on VD's silent-mixing-cost concern (variant (a) vs (b)): under
flat default fits the cost is negligible (<1 silently-rejected change
per sweep, ~1-2% of change proposals), so (a) is nearly free there. But
deep trees, mixed designs, and especially DART - the very regimes where
the defect's between-variable bias is worst - push the per-proposal
no-op rate to 4-12% with a fat tail (p90 up to ~0.7), where a subset of
change proposals become near-certain no-ops that a user would not know
to lengthen chains for. The concern is real but concentrated in the
deep/DART regime, not the default. This does not decide (a) vs (b); it
scopes the mixing cost stage 2 must measure with change-balance.R-gated
prototypes.

## Stage 2 results

Prototype: src/bartcore/moves.hpp changeMove behind a runtime switch
(DBARTS_CHANGE_KERNEL in {current, prior, hybrid}, read once when the C++
sampler is built and stored in MoveContext; all three kernels build
together). Default "current" is bit-identical: the full tinytest suite
(2473 results) and tests/cpp both pass unchanged, so unset behavior is
untouched. Env-gated change-move counters (DBARTS_CHANGE_STATS, a no-op
that never touches the RNG when unset) report acceptance/no-op rates.
Balance driver benchmarks/R/change-balance.R extended to loop the switch;
ESS driver benchmarks/R/change-fix-stage2.R reuses the stage-1
generators. Raw: benchmarks/results/change-fix-stage2.csv (+ .txt) and
change-balance-stage2.txt. ESS estimator: coda::effectiveSize throughout
(single-chain spectral); ESS/sweep = ESS / nSamples at n.thin=1.

Derivations. Both variants drop the changed node's own split-variable and
rule-prior factors; growth(node) and every factor above the node cancel
between T and T' as in birth/death, leaving the prior of the subtree
STRICTLY BELOW the node, B(T) = treeLogProbability(leftChild) +
treeLogProbability(leftChild+1), plus the branch likelihood.

- (a) prior. The proposal draws v' from p_var and the rule from the
  unrestricted node prior over the ancestor-only interval, so q(T'|T) =
  p_var(v') * p_ruleprior(newRule) is exactly the node's own prior factor
  in pi(T'); likewise q(T|T') is the node factor in pi(T). Substituting
  into alpha = [pi(T') L(T') q(T|T')] / [pi(T) L(T) q(T'|T)] cancels both
  node factors AND the variable prior identically, giving
    alpha = exp( B(T') - B(T) + yLogL - xLogL ).
  A draw whose descendant good set is empty makes pi(T') = 0: an automatic
  no-op that restores state. Ordinal: uniform over the whole splitInterval
  + missing coin, rejected iff the split falls outside findGoodOrdinalRules'
  range. Categorical: one unrestricted gauge draw (density 1/(2^R - 2)
  cancels the closed-form rule prior to sub-epsilon, as in birth/death),
  rejected iff a descendant loses its gauge. Removes the 64-try asymmetry.

- (b) hybrid. Ordinals keep the restricted good-set proposal, so
  q(T'|T) = p_var(v') / |Valid_T(v')| and q(T|T') = p_var(v) /
  |Valid_T'(v)| (missing coin symmetric). The variable prior and missing
  coin cancel against the node's rule prior; what survives multiplies the
  (a) form by |SI(v)|/|SI(v')| (ancestor-only interval sizes, from the
  node rule prior) times |Valid_T(v')|/|Valid_T'(v)| (good-set counts):
    alpha = exp( B(T')-B(T) + yLogL - xLogL
                 + log|SI(v)| - log|SI(v')|
                 + log|Valid_T(v')| - log|Valid_T'(v)| ).
  Since findGoodOrdinalRules and splitInterval both ignore the node's own
  rule, the reverse count re-enumerates the OLD variable on the changed
  tree and always contains the old rule (>= 1, never zero); same-variable
  and equal-cut changes give correction 1 (the defect's cancelling case).
  Categoricals in hybrid use mechanism (a) - counting descendant-valid
  gauge patterns is exponential for wide masks.

Balance gate (change-balance.R, 300k kept draws, P(root=x2 | root split)
vs the exact-enumeration arbiter; MCse ~ 5e-4):

  kernel   P(x2|split)  exact   z vs exact   verdict
  current    0.2978    0.0774    +255.0      FAIL - reproduces the defect
  prior      0.0770    0.0774     -0.7       MATCH
  hybrid     0.0768    0.0774     -1.2       MATCH
  control (equal 19-vs-19 cut counts): current z = -15.7 (the honest
  deep-node residual the audit noted), prior +1.2, hybrid -0.8 - all
  match, and the repairs also remove current's residual. Both repairs
  pass; the defect fails decisively. Gate cleared, ESS run below.

ESS (1000 burn + 3000 sweeps, n.thin=1, single chain/thread, seed
20260708, serial, otherwise-idle machine). acc = change-acceptance rate;
noop = realized silent-no-op rate (prior's is the quantity stage 1
predicted); S = ESS/sec, W = ESS/sweep, for sigma / mean over 6 fixed
train-fit coords / mean over per-variable varcount traces (k fixed in
every config, so not sampled):

  config             kernel   acc   noop  |  sigS  trnS  varS | sigW   trnW   varW
  Friedman 1e3 m75   current 0.109 0.0005 | 173.7  74.2 149.3 |0.033 0.014 0.029
                     prior   0.105 0.0123 | 100.2  86.7 179.4 |0.019 0.016 0.034
                     hybrid  0.103 0.0002 | 101.8  88.8 157.0 |0.020 0.017 0.030
  Friedman 1e3 m200  current 0.164 0.0005 |  37.5  65.2  46.4 |0.019 0.034 0.024
                     prior   0.156 0.0106 |  35.4  65.3  41.7 |0.018 0.032 0.021
                     hybrid  0.160 0.0003 |  20.6  58.0  45.2 |0.010 0.029 0.023
  Friedman 1e4 m75   current 0.024 0.0014 |  55.0  10.6  24.9 |0.114 0.022 0.052
                     prior   0.022 0.0177 |  23.8  14.7  26.3 |0.048 0.030 0.053
                     hybrid  0.026 0.0012 |   7.9   6.0  17.3 |0.016 0.012 0.036
  Friedman 1e4 m200  current 0.080 0.0006 |  30.8   3.1   4.1 |0.160 0.016 0.021
                     prior   0.083 0.0085 |  37.6   4.0   4.6 |0.192 0.020 0.024
                     hybrid  0.085 0.0005 |  35.0   4.2   5.5 |0.186 0.022 0.030
  Mixed 1e4 m75      current 0.020 0.0029 |  17.6   8.2  14.5 |0.038 0.018 0.031
                     prior   0.018 0.0187 |  18.7   6.7  15.1 |0.041 0.015 0.033
                     hybrid  0.018 0.0039 |   7.0   8.5  18.3 |0.015 0.019 0.040
  Deep 1-tree 1e3    current 0.071 0.0394 | 322  4519  1130  |0.002 0.028 0.007
                     prior   0.057 0.1099 | 467  1117  1044  |0.003 0.006 0.006
                     hybrid  0.075 0.0233 | 567   955  1628  |0.003 0.005 0.009
  Mixed+DART 1e4 m75 current 0.024 0.0042 |  27.9   6.2  34.8 |0.059 0.013 0.073
                     prior   0.019 0.0420 |  17.0   4.2   5.8 |0.035 0.009 0.012
                     hybrid  0.023 0.0041 |  35.1   6.5  15.5 |0.075 0.014 0.033

Reading. The realized prior no-op rates reproduce stage 1's forecast
(~1-2% on the flat Friedman/Mixed configs, 11% on the deep single tree,
4.2% under DART); current and hybrid stay below 0.4% everywhere except
the deep tree. Change-acceptance rates are indistinguishable across the
three kernels in every config, so the correction reweights which
proposals stand, not how often a change lands. On the flat default
configs the three kernels' ESS/sweep for train fits and varcount are
comparable, with no consistent ordering and sigma ESS too noisy at a
single 3000-draw chain to separate them. In exactly the two regimes stage
1 flagged, the cost surfaces: on the deep single tree prior and hybrid
lose train-fit ESS/sweep against current (0.028 -> 0.006 / 0.005), and
under DART prior's varcount ESS/sweep collapses (0.073 -> 0.012) while
hybrid - which spends no-ops only on categoricals - sits between (0.033).
Caveat: "current" targets the wrong stationary distribution, so its ESS
measures mixing on a different (incorrect) target and is only loosely
comparable to the repairs'; single-chain ESS at 3000 draws, sigma
especially, carries large variance. No recommendation - measurements
only.

## Status: landed (2026-07-08)

VD picked variant (b), the hybrid, from the stage-2 measurements. It is
now THE change move: the DBARTS_CHANGE_KERNEL switch, the "current" and
"prior" code paths, and the DBARTS_CHANGE_STATS counters are gone;
changeMove carries the hybrid acceptance unconditionally (moves.hpp).
Ordinals keep the restricted good-set proposal with the counted
log|SI(v)|/|SI(v')| + log|Valid_T(v')|/|Valid_T'(v)| correction;
categoricals propose from the node prior (drawCategoricalRuleFromPrior).
change-fix-stage2.R and its result files are kept unchanged as the
experiment record; change-balance.R is restored to a single-engine form
and is now the permanent exact-posterior gate. Design note:
docs/design/change-move-balance.md.

Gate results at the landing. change-balance.R: MAIN engine matches
exact at z = -1.2, CONTROL z = -0.8 (the defect's deep-node residual is
gone), gate PASS. logistic-reference.R OK (max gap 0.0003, tol 0.005);
categorical-exact.R OK (max gap 0.0008 - the propose-from-prior
categorical mechanism is exact); bcf-exact.R OK (mode 2b at its
documented looser tolerance). tests/cpp all pass. Full tinytest suite
green (2489 results) after regenerating the RNG-locked snapshots
(reproducibility singleThreaded/binaryResponse/xbart, rbart-loop-
callback) and re-tuning four seed-fragile statistical assertions
(binaryResponse-hyperprior cor/median-k bounds, sampler-splitProbabilities
run length, capi transactional-mutation column seed). equivalence.R vs
equivalence-bcdcc07.rds: fit-level summaries stay |z| < 4 across all 18
scenarios; the variable-proportion (vprop.*) and sigma summaries exceed
it exactly in the unequal-cut designs (splitprobs, quants, categorical,
dart) - the between-variable reweighting the fix intends. Baseline NOT
re-recorded here; the maintainer re-records after landing. bench-sampler
deferred to a quiet-machine grant.

Review catch before release: the landed correction branched on the NEW
variable's type only, but the proposal-density ratio composes per side
(forward = new variable's mechanism, reverse = old variable's), so both
MIXED cases were wrong - ordinal-to-categorical dropped the old
ordinal side's counted terms, and categorical-to-ordinal ran
splitInterval/findGoodOrdinalRules on a categorical column (meaningless
gauge-mask arithmetic feeding the logs). Fixed to the per-side form
  (v' ordinal ? log|Valid_T(v')| - log|SI(v')| : 0)
    + (v ordinal ? log|SI(v)| - log|Valid_T'(v)| : 0);
all-ordinal fits are bit-identical to the first landing build (verified
against saved oracles), categorical-fit draws shift. change-balance.R
gained a MIXED arm built to see exactly this defect - a 9-cut ordinal
against a 4-level factor aligned with its step blocks, so cross-type
changes at stacked nodes are likelihood-neutral and carry the root
balance: the one-sided kernel fails it at z = -8.0 (300k draws), the
per-side kernel passes at z = +0.7, and MAIN/CONTROL reproduce
bit-for-bit (z = -1.2 / -0.8). Re-run after the fix: categorical-exact,
logistic-reference, bcf-exact all OK unchanged; tests/cpp pass; tinytest
green at 2464 results (the count moved from 2489 because
test-data-categorical*'s per-sampled-root assertion loops see different
draws; no snapshot regeneration was needed and no all-ordinal file
changed); equivalence vs bcdcc07 identical verdicts except the
categorical scenario's offenders moved (max |z| 52.25 -> 53.99, still
vprop/sigma only).
