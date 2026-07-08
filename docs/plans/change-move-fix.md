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
