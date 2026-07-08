# The change move's detailed balance (fix, 2026-07-08)

The tree change move proposes a new split variable and rule at an
internal node and accepts by a Metropolis-Hastings ratio. Since the
package's origin that ratio omitted the proposal-density term, so the
chain targeted a distribution that over-weights low-cardinality and
descendant-constrained split variables whenever the design's effective
cut counts are unequal. This note records the defect, the acceptance
that repairs it, the two proposal mechanisms it uses, the evidence
behind the choice, and the gate that pins it. The blow-by-blow numbers
live in the plan files (docs/plans/correctness-audit.md Status block 2;
docs/plans/change-move-fix.md stages 1 and 2); this note is the
standing account.

## The defect

The move draws the new variable v' from the split-variable prior
p_var, then draws the new rule NOT from that variable's rule prior but
uniformly over its descendant-valid good set (ordinal:
findGoodOrdinalRules' interval; categorical: reject-until-valid). The
old acceptance kept only the pi ratio,
exp(yLogPi + yLogL - xLogPi - xLogL), with no q(T|T')/q(T'|T) factor.
The node's rule prior in pi is normalized over the ancestor-only
interval, not the good set, so the proposal and prior densities cancel
ONLY for a same-variable redraw. A cross-variable change is
mis-weighted by [p_var(v')/p_var(v)] * [|Valid(v)|/|Valid(v')|]: the
stationary distribution carries an effectively SQUARED rule prior at
changed nodes (mass ~ (p_var/a_v)^2 rather than p_var/a_v, a_v = the
variable's available cut count), biased toward variables with few cuts
or heavily descendant-constrained good sets. It is invisible when every
variable has equal cut counts and trees are shallow - exactly the
regime the pre-existing exact-posterior gates covered, which is why it
went uncaught.

The defect is INHERITED, not a rewrite regression: the deleted classic
engine's changeRule.cpp computes the identical pure-pi-ratio acceptance
(git b354f3a~1), so dbarts carried it from its CGM-lineage origin.
Birth and death always targeted the correct pi, so the realized bias
sits between the correct and the fully-squared target rather than at
the latter.

## The hybrid acceptance

Both a node's own split-variable factor and its rule-prior factor, and
every prior factor at or above the node, cancel between T and T' as in
birth/death. What survives is the prior of the subtree STRICTLY BELOW
the changed node,

  B(T) = treeLogProbability(leftChild) + treeLogProbability(leftChild + 1),

plus the branch likelihood. The landed acceptance is

  alpha = exp( B(T') - B(T) + yLogL - xLogL + logProposalCorrection ).

The correction is the surviving proposal-density ratio q(T|T')/q(T'|T)
and COMPOSES PER SIDE: the forward density q(T'|T) uses the new
variable v''s proposal mechanism, the reverse density q(T|T') the old
variable v's, and each side either cancels its own node-prior factor
or leaves a counted term, independently of the other side:

  logProposalCorrection =
      (v' ordinal ? log|Valid_T(v')| - log|SI(v')| : 0)
    + (v  ordinal ? log|SI(v)| - log|Valid_T'(v)| : 0).

An ordinal side draws uniformly over the descendant-valid good set
(q = p_var / |Valid|) while the node's rule prior normalizes over the
ancestor-only split interval |SI|; the variable prior and the symmetric
missing-direction coin cancel within the side, leaving the counted
ratio. findGoodOrdinalRules and splitInterval both ignore the node's
OWN rule, so the reverse count re-enumerates the old variable v on the
CHANGED tree T' (equal to its value on T) and always contains the old
rule: |Valid_T'(v)| >= 1, never zero. A same-variable or equal-cut
ordinal change gives correction 1 - the case the old ratio already got
right.

## Categoricals: propose from the prior

Counting descendant-valid gauge patterns is exponential in the mask
width, so the counted correction is infeasible for categoricals.
Instead a categorical side draws one unrestricted gauge pattern
straight from the node prior over the reachable category set
(drawCategoricalRuleFromPrior). That prior density, 1/(2^R - 2),
cancels the node's closed-form rule prior exactly (to sub-epsilon, as
in birth/death), so the side contributes nothing to the correction. A
forward draw whose directions strand a same-variable descendant makes
pi(T') = 0 - an automatic no-op that restores state. This also removes
the old rejection sampler's variable-dependent 64-attempt abort
asymmetry.

The per-side composition matters exactly in the CROSS-TYPE cases. An
ordinal-to-categorical change keeps the old ordinal side's counted
ratio (log|SI(v)| - log|Valid_T'(v)| > 0 whenever a same-variable
descendant split constrains the node) even though the forward side
cancels; a categorical-to-ordinal change keeps only the forward ordinal
ratio, and the ordinal counters must never run on the old categorical
column (their interval arithmetic is meaningless for gauge masks). The
initial landing branched the whole correction on the new variable's
type alone - wrong in both mixed directions; review caught it before
release and the gate below gained a mixed-type arm that detects it.

The two mechanisms coexist because propose-from-prior is exact for both
column kinds but wastes proposals as no-ops at nodes with same-variable
grandchildren, whereas the counted correction is exact and no-op-free
for ordinals but intractable for wide categorical masks. The hybrid
confines the residual no-op cost to categoricals.

## Evidence and decision

Stage 1 (change-move-fix.md) instrumented stock fits: the
propose-from-prior no-op rate is negligible under flat default forests
(~1-2 percent of change proposals, < 1 per sweep) and concentrates in
the deep-tree and DART regimes (per-proposal 4-12 percent, p90 tails to
~0.7) - the same regimes whose between-variable bias the defect most
distorts. Stage 2 built both exact repairs behind a switch. Both passed
the balance gate and reproduced the defect's failure; change-acceptance
rates were indistinguishable across kernels (the correction reweights
which proposals stand, not how often one lands); ESS/sweep was
comparable on flat configs, with the pure-prior variant losing
varcount ESS under DART where the hybrid - spending no-ops only on
categoricals - sat between. VD picked the hybrid on 2026-07-08.

## The gate

benchmarks/R/change-balance.R is the permanent exact-posterior gate,
three scenarios of long single-tree fits against region-DP enumerations
of the true single-tree posterior. MAIN (x1 19 cuts vs x2 2 cuts)
compares root split frequencies against the exact posterior and the
squared-rule-prior wrong target; the engine must MATCH exact (|z| < 4)
and sit far from wrong. CONTROL (equal 19-vs-19 cut counts) must stay
clean. MIXED pits a 9-cut ordinal against a 4-level factor whose levels
exactly align with the ordinal's step blocks, so every stacked tree
shape has a likelihood-neutral cross-type partner and the x1-vs-x2 root
balance is carried almost purely by cross-type change acceptances -
the DP there enumerates both rule kinds (ordered gauge assignments,
2^|S| - 2 per reachable set S). The one-sided correction fails this arm
at z = -8.0 (300k draws) while the per-side kernel passes; at the
landing the engine matched exact at z = -1.2 (MAIN), z = -0.8 (CONTROL,
the residual the defect left behind now gone), and z = +0.7 (MIXED).
