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
