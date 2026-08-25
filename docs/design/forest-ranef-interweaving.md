# Forest-ranef confounding: a joint / interweaving move (design)

Status: NO-GO now, recorded door (2026-07-20; blind critique SOUND WITH
CAVEATS - see section 9, authoritative, which corrects the section-2 magnitude
framing and the section-3.2 cost framing)

The in-engine grouped sampler (rbart_vi with a built-in tau prior, the
GroupedResponse path) fits y_i = f(x_i) + b_{g(i)} + eps_i with b_j ~ N(0,
tau^2). In the WEAK-signal / FEW-groups corner the group scale tau and the
intercepts b mix poorly. A prior review (docs/plans/tau-slice-review.md)
established by measurement that the tau draw is NOT the bottleneck (an exact
conjugate draw mixes identically -- since landed, 103a9ef / docs/plans/tau-
cauchy-exact-ig.md), and that the DOMINANT weak-signal bottleneck is FOREST-
RANEF CONFOUNDING: the mean forest f and the group intercepts b competing to
explain group-level structure. This note designs the remedy. It is design +
measurement only; nothing here is implemented. rng class of any future build:
posterior-preserving but DRAW-CHANGING for grouped fits (statistical, not
bitwise, re-record of the grouped equivalence scenarios); every ungrouped path
and the tau kernel stay BIT-IDENTICAL (the move touches only the grouped
forest/ranef blocks).

Landscape (general knowledge, not re-verified this pass): the confounding of a
flexible mean with a random intercept is the classic variance-component mixing
pathology. The textbook remedies are (i) the COLLAPSED / marginal sampler --
integrate the intercept out of the mean update (the gold standard for
random-intercept confounding); (ii) INTERWEAVING / ASIS (Yu-Meng 2011) and
parameter expansion (Liu-Wu 1999; Papaspiliopoulos-Roberts-Skold 2007), which
alternate a centered and a non-centered parameterization to travel the ridge.
BART-adjacent, stan4bart (Dorie, the dbarts LinkingTo consumer) sidesteps the
Gibbs confounding by alternating a dbarts f-draw with a Stan HMC block over
(b, Sigma) that marginalizes the ridge in one joint proposal -- the strongest
existing evidence that a JOINT (f, b) move is what pays, and the model this note
tries to bring in-engine. No maintained pure-Gibbs BART package fixes this;
BayesTree-era rbart (this code's ancestor) carried the same confounded blocked
Gibbs.

**Erratum (2026-08-09, tree-mixing synthesis).** "Marginalizes the ridge
in one joint proposal" overstates stan4bart: the HMC block marginalizes
only WITHIN the parametric block (b, Sigma); the (f, b) ridge itself is
still traversed by alternation - stan4bart runs one BART sweep per
parametric draw (stan4bart init.cpp:642). The strongest-evidence claim
should read: joint moves pay within the parametric block, and whether
composition helps the cross-block ridge is measured, not established -
see docs/design/tree-mixing-proposals.md (parametric absorption).

## 0. The load-bearing finding: there is NO cheap interweave; the fix is a collapse

The natural hope is that this mirrors the BCF a-ridge move (docs/plans/bcf-ridge-
interweaving.md), where the prognostic amplitude `a` and the leaf values `mu`
sit on an explicit shared ridge `a*mu` and a single scalar rescale `c` travels
it. That structure does NOT exist here, and recognizing why is the whole design:

- BCF's ridge is between two EXPLICIT parameters (`a` and the `mu` leaves) that
  both live in the model as drawable scalars; a rescale move is a one-line
  reparameterization with a closed-form GIG conditional.
- The f-b ridge is between a NONPARAMETRIC forest and a per-group intercept.
  The confounded quantity is each group's location, `fbar_j + b_j` where
  `fbar_j = mean_{i in j} f(x_i)`. But f is organized by LEAVES, not groups:
  there is no per-group knob in f to swap a constant into or out of. A move that
  "shifts a constant s_j from f into b_j" is not expressible in leaf space -- to
  change `fbar_j` you must redraw leaves, which IS the tree sweep. So the
  ASIS/PX family (which needs a shared scalar coordinate to reparameterize) has
  nothing to grip. The isolation prototype confirms this operationally: an
  ASIS (tau, b) move does not touch the f-b ridge and can even WORSEN it at
  small K (section 5).

The only move that attacks the f-b ridge directly is therefore a COLLAPSE:
draw f with b MARGINALIZED, then draw b | f -- a joint (f, b) block. That is
exact and strong, but marginalizing a per-group intercept induces within-group
correlation that COUPLES the forest's per-leaf sufficient statistics across the
group partition, landing the change on the hottest path in the engine. This
inverts the "cheap like the BCF ridge" expectation: the correct remedy is the
widest grouped engine surgery in the queue, and its payoff is confined to a
narrow corner (section 6). That tension is the go/no-go.

## 1. The model and its reduction to the engine sweep

The Gibbs blocks per raw sweep, as run() drives them:

1. `f | b` -- one tree sweep of the mean forest against the working response
   `z_i - b_{g(i)}` (GroupedResponse::workingResponse, model.hpp:4902-4906),
   backfit tree-by-tree (chain.hpp:1418-1530). This is where f sees the group
   intercepts subtracted, so f fits the residual-of-b.
2. `b | f` -- refreshLatents (chain.hpp:1535) calls GroupedResponse::
   refreshLatents (model.hpp:4746-4779), which draws b_j conjugately from the
   group means of `z_i - F_i` with F = f-only fits (drawGroupEffects,
   model.hpp:4668-4691; called at :4750-4753 with the combined = f-only fits).
3. `tau | b` -- exact Makalic-Schmidt cauchy draw (drawTauCauchyExactIG,
   model.hpp:4647-4657) or slice for the gamma prior.
4. `sigma | f, b` -- drawSigma on the shifted fits (chain.hpp:1545,
   model.hpp:4782-4786).

Blocks 1 and 2 are the ridge: f conditions on the current b, b conditions on the
current f, and neither integrates the other out. When x carries group-level
structure (the common applied case: groups correlate with covariates), f can
absorb part of each `fbar_j` and b absorbs the rest; the pair traverses the
`fbar_j + b_j = const` ridge one block at a time, slowly. The GroupedResponse
decorator (model.hpp:4706) presents only `(z, w)` to the forest, so it CANNOT
see the coupling -- the collapse cannot live in the decorator (section 3).

## 2. The measured bottleneck (HEAD, benchmarks/R/grouped-mixing.R)

Measured on the installed 1.0-0 (in-engine cauchy = exact-IG), tau chain,
n=1000, n.trees=50, 6000 kept draws, cauchy prior, IACT = N / coda::
effectiveSize averaged over 6-8 replicate seeds. The single-seed IACT of the
weak/small-K tau chain is heavy-tailed and very noisy, so ranges are reported.

Part A -- tau IACT over K x signal (all WITH the forest):

    K    signal      IACT    postmean (true)
    3    weak(0.2)   88.0    0.815 (0.2)
    10   weak(0.2)    7.4    0.291
    40   weak(0.2)   34.9    0.217
    3    strong(2.0) 25.4    4.586 (2.0)
    10   strong(2.0)  6.3    2.644
    40   strong(2.0)  2.1    2.190

tau mixes well with a strong group signal at moderate/many groups (IACT 2-6) and
poorly in the weak-signal corner, worst at small K (biased posterior mean, IACT
88) and elevated again at K=40 (weak groups barely identified from noise).

Part B -- the ATTRIBUTION: weak signal, the same fits WITH the strong forest f
vs WITHOUT it (f == 0, intercept-only control):

    K    IACT with f    IACT no f    ratio    IACT range (with f | no f)
    3    71.2           21.6         3.3      [3-310 | 11-42]
    10    9.1            4.9         1.9      [5-24  |  3-10]
    40   39.9           16.1         2.5      [9-91  |  5-40]

Removing f drops the tau IACT most at small K (~3.3x at K=3, ~2x at K=10), and
the WITH-f chain shows catastrophic single-seed excursions (IACT to 310) that
the intercept-only chain never produces (max 42). This CONFIRMS the review's
direction -- forest-ranef confounding is the dominant single lever in the weak-
signal corner, larger there than the tau kernel (~0) or the tau-b funnel -- and
UPDATES its magnitude: the mean multiplier is ~3x on HEAD, not the review's
slice-era ~25x single-DGP point estimate. The confounding is real and dominant-
in-the-tail; the precise multiplier is highly DGP- and estimator-dependent
because the small-K tau posterior is barely identified and heavy-tailed. (The
review's ~25x came from a lost throwaway DGP on the slice sampler; the qualitative
finding reproduced across every DGP tried, the magnitude did not.)

## 3. Primary remedy: the collapsed / joint (f, b) move

Marginalize b out of the mean-forest update, then redraw b | f -- a joint (f, b)
block that collapses the ridge.

### 3.1 The math

With b integrated out, the data covariance is block compound-symmetric: within
group j (n_j observations),

    Sigma_j = sigma^2 I_{n_j} + tau^2 1 1^T,
    Sigma_j^{-1} = (1/sigma^2)( I - (tau^2 / (sigma^2 + n_j tau^2)) 1 1^T )   (Woodbury).

So the forest sees a rank-1 per-group DOWNWEIGHT of the group-mean direction:
the group-mean component of each group's residual carries variance
`sigma^2/n_j + tau^2` (inflated by tau^2), the within-group component carries
`sigma^2`. This is exactly the information f should NOT double-count with b.
Drawing f from p(f | y, tau, sigma) under this covariance, then b_j | f from its
conjugate group posterior (unchanged, drawGroupEffects), is a valid blocked
Gibbs step drawing (f, b) jointly given tau -- the collapse. tau | b is
untouched (block 3), independent of f given b.

### 3.2 What it costs in the engine (the sharp edge)

The constant leaf's marginal and draw depend on a node ONLY through
`(sumWeights, sumWeightedResponse)` against a SCALAR `residualVariance`
(logIntegratedLikelihood model.hpp:165-183; drawFromPosterior model.hpp:186-198;
posteriorPrecision = sumWeights/residualVariance, model.hpp:174,190), and the
birth/death score sums that per-leaf marginal independently across leaves
(logLikelihoodForBranch moves.hpp:67; metropolisJumpForTree moves.hpp:844; the
cached node stat moves.hpp:280-281). The compound-symmetry covariance BREAKS
this: the Woodbury correction couples every leaf that a group touches, so a
node's `(sumWeights, sumWeightedResponse)` is no longer a sufficient statistic
-- the leaf marginal needs, per group overlapping the leaf, the group's total
residual sum (which spans OTHER leaves). Two honest implementation shapes:

- (P-exact) Thread a per-group residual-sum channel through the leaf marginal
  and draw so each leaf's score carries its `sum_{i in leaf ∩ group j}` terms
  and the per-group `tau^2/(sigma^2 + n_j tau^2)` correction. This is a new
  coupled leaf model (a `GroupCollapsedGaussianLeaf`) plus per-sweep per-group
  reductions maintained across the tree loop. It is EXACT and it changes the
  hottest path (every birth/death score, every leaf draw) -- the widest grouped
  engine change, comparable to heteroscedastic's C2 (heteroscedastic.md) but
  ON the constant-leaf path rather than a nullable second forest. It cannot live
  in the GroupedResponse decorator, which only sees `(z, w)`; it must reach into
  the leaf model, so grouping stops being a pure response decorator.

- (P-approx) A cheaper heuristic: pre-shrink only the GROUP-MEAN component of the
  forest's residual by the compound-symmetry factor before the forest sees it,
  leaving the within-group residual alone -- expressible as a per-group
  correction to the working response the decorator forms. This keeps the change
  in the decorator (no leaf-model surgery) but is NOT exact (a per-group
  location shrink is not the full rank-1 covariance; the forest's own splitting
  redistributes the shrunk mean across leaves imperfectly). It would need its own
  exact-posterior gate to prove it targets the right law, and may only recover
  part of the collapse's gain. Recorded as a door, not the recommendation --
  shipping an unvalidated approximate kernel is against the repo bar.

### 3.3 Draws changed / RNG

Draw-changing for grouped fits (the f-draw stream changes). Ungrouped paths and
the tau/sigma kernels are untouched -> bitwise-neutral. No new persistent state
(the collapse is a within-sweep change to how f is drawn; b, tau, sigma channels
are unchanged). Per-sweep cost: the per-group reductions are O(n) plus the leaf
score's per-group terms O(#leaves touched per group) -- a real constant-factor
tax on every tree move, paid by grouped fits only.

## 4. Secondary remedy: ASIS / parameter expansion on (tau, b)

The Yu-Meng interweave on the ranef SCALE, the docs/plans/bcf-ridge-
interweaving.md precedent applied to (tau, b), review section 4c. After the
centered draws (b | f conjugate, tau | b exact-IG -- the landed cauchy step is
its prerequisite building block), add a non-centered ANCILLARY draw: set
eta_j = b_j / tau (fixed), then redraw tau as a ridge coefficient with a
half-Cauchy prior via the drawGlue IG mixture, then b <- tau * eta. The
conditionals (validated against quadrature, review 3c):

    v   | tau      ~ IG(1, (tau^2 + A^2)/2)
    tau | eta, y,v ~ N(m, 1/prec) truncated > 0,
        prec = sum_j n_j eta_j^2 / sigma^2 + 1/v,
        m    = (sum_j n_j rbar_j eta_j / sigma^2) / prec;   then b <- tau eta,

with `rbar_j` the group mean of `y - f`, `A = priorScale_`. This is CHEAP and
decorator-local: a few O(J) reductions plus a truncated-normal draw, no leaf-
model change, and it builds directly on the landed exact-IG cauchy step.

CRITICAL LIMITATION: it attacks the tau-b FUNNEL, not the f-b confounding. The
review measured a large gain on the isolated funnel (IACT 23 -> 2, 3b) but that
funnel is the SECONDARY full-model mode; ASIS on (tau, b) leaves the dominant
f-b ridge untouched. The prototype (section 5) shows it can even WORSEN mixing
at small K, where the confounding dominates and the extra ancillary randomization
of the ridge coordinate is counterproductive -- exactly the failure the
bcf-ridge landing saw when it tried refreshing the conditioned auxiliary
(bcf-ridge-interweaving.md). So ASIS is a genuine but PARTIAL win whose
sign is regime-dependent; it is not a substitute for the collapse.

## 5. Isolation prototype and measured gain (benchmarks/R/forest-ranef-collapse-proto.R)

A tractable surrogate isolates the SAME beta-b ridge with a closed-form
collapse, the way the review prototyped ASIS in isolation: y_i = beta0 +
beta1 x_i + b_{g(i)} + eps_i, sigma known, x clustered by group so beta1 x
(the "forest") competes with b. Three tau chains, identical target, only the
mean+ranef update differs; tau drawn by the engine's exact cauchy step in every
arm. CENTERED = the engine's separate blocks; COLLAPSED = beta with b
marginalized (the primary, exact Woodbury), then b | beta; ASIS = centered plus
the (tau, b) interweave (the secondary). Median tau agrees across arms (same
target: correctness self-check); only IACT differs:

    K    CENTERED    ASIS     COLLAPSED    median tau (agree)
    3    56.1        114.6    9.3          0.223 / 0.254 / 0.239
    10   13.1        2.8      11.7         0.214 / 0.214 / 0.215

Reading: at K=3 (strong confounding) the COLLAPSE cuts IACT ~6x (56 -> 9) --
the primary remedy directly collapses the ridge -- while ASIS does NOT (56 ->
115, worse: it fixes a funnel that is not the bottleneck and perturbs the
dominant ridge). At K=10 the confounding is mild so the collapse is ~flat
(13 -> 12) while ASIS shows its funnel gain (13 -> 3). This is the note's thesis
in one table: the collapse is the lever for the confounding-dominated corner;
ASIS is a funnel fix that does not transfer to it. (Optimistic bound: the
surrogate's mean is LINEAR, so its collapse is exact and cheap; the real forest
collapse costs section 3.2 and its in-engine gain will be smaller than 6x
because the forest is not a 2-parameter mean.)

## 6. Recommendation: go/no-go

Honest weighing, given the measurements:

- The confounding is REAL and the dominant single lever in its corner, but the
  corner is NARROW: weak group signal AND few (or barely-identified many) groups.
  With a strong group signal, or moderate K, tau mixes fine (IACT 2-7, section 2).
- The mean multiplier on HEAD is ~3x (section 2), not the review's ~25x; the
  drama is in the heavy-tailed excursions. There is a trivial user workaround for
  the corner -- run a longer chain (IACT ~70 at K=3 still yields ESS ~85 per 6000
  and is embarrassingly parallel across n.chains).
- The exact collapse (P-exact) is the WIDEST grouped engine change: it moves
  grouping off the clean ResponseModel-decorator seam into a coupled leaf model
  on the hottest path (section 3.2), is draw-changing (grouped re-record), and
  carries the full exact-posterior gate burden. That is a large, VD-sign-off,
  from-scratch research arc for a narrow-corner ~3x with a longer-chain
  workaround.

RECOMMENDATION:

- NO-GO, for now, on the from-scratch exact collapse (P-exact). The payoff (a
  narrow-corner ~3x with a trivial workaround) does not justify moving grouping
  off the decorator seam onto the constant-leaf hot path and re-recording the
  grouped snapshots. Keep it a RECORDED DOOR: reopen with VD sign-off if a
  concrete weak-signal / few-groups application (e.g. a small-J multisite trial,
  an IRT-style grouping) needs calibrated tau intervals that a longer chain
  cannot cheaply buy. Estimated cost when it lands: comparable to
  heteroscedastic's C2, ~800-1300 lines engine (a coupled leaf model + per-group
  reductions through the tree loop) + the exact-posterior gate + grouped
  re-record; the biggest risk is protecting the ungrouped bitwise guarantee while
  the constant-leaf path grows a group-aware branch.

- CONDITIONAL-GO, only IF grouped mixing is prioritized, on the SECONDARY ASIS
  (tau, b) (section 4) as a bounded, decorator-local win that builds on the
  landed exact-IG -- BUT gated on a full-engine prototype FIRST, because the
  isolation shows its sign flips against the confounding (it helped the funnel at
  K=10 and hurt at K=3). Estimated cost: small (~60-120 lines engine + a
  truncated-normal helper + component/quadrature gates, the tau-cauchy-exact-ig
  shape), but its expected full-model payoff is modest and possibly negative in
  the very corner that motivates the work. If the prototype shows a net win
  across the K grid it is a cheap ship; if it flips sign at small K (as the
  isolation warns), DROP it.

- Default: DO NEITHER; document the corner and the longer-chain workaround. The
  tau kernel is already exact (103a9ef); the residual mixing issue is a
  known-narrow confounding corner, not a defect. This note is the record so the
  next agent does not re-derive that the cheap interweave does not exist.

## 7. Gates any engine build would need

Per review section 5 (the migration battery), unchanged in shape:

- Exact-posterior GATE (the load-bearing one): the 1-D / 2-group marginal-tau
  quadrature (b analytically integrated, review 3c; ybar_j ~ N(0, tau^2 +
  sigma^2/n_j)). A collapse or ASIS build must reproduce it to MC error -- this
  is what catches a self-consistently-wrong collapse (the class of bug the
  heteroscedastic m'=2 gate and the monotone stationary-bias gate guarded).
  Promote it into a tests/cpp checkNear as tau-cauchy-exact-ig did.
- tests/cpp: extend testGroupedMath (the new kernel's moments vs the R
  quadrature and vs the current draw on fixed rng); a poison check (perturb the
  Woodbury correction, watch the quadrature gate fail, restore).
- Full tinytest; regenerate the grouped hardcoded snapshots by WHOLE-FILE replay
  (last-ulp depends on process history, grouped-random-effects.md).
- equivalence.R compare: ungrouped scenarios IDENTICAL (bitwise); the grouped /
  grouped_aft scenarios re-record as a STATISTICAL verdict (z-tests on tau /
  ranef / fit posterior means+intervals vs the current baseline), the grouped
  landing's shape (grouped-random-effects.md). NOTE the grouped
  equivalence scenarios use the GAMMA prior (tau-cauchy-exact-ig.md);
  a cauchy-branch collapse would need a cauchy-grouped z-summary added, or the
  scenarios re-recorded under cauchy.
- The custom-prior R loop (rbart.R:868-1042) stays untouched and must keep
  working -- a custom prior forcing the cauchy density is the cross-check the
  grouped landing used.
- grouped-mixing.R (this arc's gate) re-run: the collapse must drop the Part B
  with-f IACT toward the no-f floor at K=3; ASIS must show its (possibly
  sign-flipped) funnel effect. This is the acceptance measurement.

## 8. Doors / out of scope

- The P-approx group-mean-shrink decorator variant (section 3.2): a cheaper
  non-exact heuristic; reopen only if a validated approximate kernel is wanted
  and the exact collapse is judged too invasive.
- Slope / non-intercept random effects, crossed or nested groupings: the collapse
  generalizes in principle (the marginal covariance stops being block-diagonal
  for crossed effects, breaking the per-group Woodbury), out of scope here.
- A stan4bart-style HMC block over (b, tau, sigma) alternated with the dbarts
  f-draw: the joint move done OUTSIDE the engine (the LinkingTo consumer already
  exists). If the in-engine collapse is judged not worth it, this is the honest
  place the joint move already lives for users who need it -- recorded as the
  external alternative to a from-scratch engine arc.
- dbarts.h exposure: none (grouping is rbart_vi-internal; universal precedent).

## 9. Critique resolutions (2026-07-20, authoritative)

An independent blind critique reproduced both scripts to the digit, verified the
engine internals, and ran a worst-case-regime sensitivity probe. Verdict: SOUND
WITH CAVEATS - the NO-GO-now / recorded-door recommendation is CORRECT and ROBUST
(it survives even a 9x reading; the corner is narrow and the longer-chain
workaround is honest). Three framing corrections, which supersede the sections
they touch:

- MAGNITUDE (supersedes section 2's "~3x" headline). The 3x is real but fragile
  and OVERSOLD. The with-f/no-f contrast structurally CANNOT reach the worst case:
  where the forest most tightly aliases group means (tight clustering), the
  intercept-only "no-f" control itself fits b through the predictor, so its
  between-group fitted variance jumps ~10x and the ratio COLLAPSES to ~1.8 - an
  estimand artifact, not low confounding. And the number is estimator-noise
  dominated (a 5-seed rerun of the same K=3 cell gave 7.5 vs the 8-seed 3.3). The
  honest statement is "several-fold, ~3-9x, estimator-unstable, worst case
  unmeasurable by this contrast"; no regime reproduced the review's 25x. The
  cleaner magnitude is the PROTOTYPE's direct centered-vs-collapsed 6x (section 5,
  immune to control contamination) - lead with that, not the section-2
  attribution. (Note lines ~58 mislabel S_WITHIN=0.8 as "where f and b compete
  hardest"; it is the largest aliasing at which the estimand still works.)

- COST/NOVELTY (softens section 3.2). "Widest surgery / off the clean
  ResponseModel-decorator seam / from-scratch" OVERSTATES novelty: a coupled-
  scoring seam ALREADY exists - ParamScoringLeafModel + logLikelihoodForBranch-
  WithParams (model.hpp:145-152,729), used live by MonotoneConstantGaussianLeaf
  in the sweep (chain.hpp:1479-1480). A GroupCollapsedGaussianLeaf would RIDE that
  path, not invent it. The CORE claim stands - the constant-leaf suffstat
  (sumWeights, sumWeightedResponse) is insufficient for the per-group Woodbury
  correction, which couples leaves across the group partition AND across trees
  (wider than monotone's within-tree frozen-neighbor coupling) - so it is real,
  substantial engine work (~800-1300 lines) on the hot path with the
  ungrouped-bitwise guarantee to protect, but not greenfield.

- PROTOTYPE BOUND (drops a section-5 claim). The assertion that the real-forest
  collapse gain is NECESSARILY below the prototype's 6x is unjustified: a
  nonparametric forest can alias per-group means MORE freely than a 2-parameter
  line (K free constants vs one trend), which argues for MORE confounding and a
  potentially LARGER gain; the only credible attenuator is BART's leaf-shrinkage
  prior, which the note never argues. So the real-forest gain could be above or
  below 6x - genuinely uncertain.

Confirmed in the note's favor (no change): there is no cheap ASIS/PX here (the
forest exposes no per-group knob - workingResponse = z - b_g, fbar_j is emergent
leaf-space, not a parameter; contrast BCF's explicit a-scalar), and ASIS's
small-K worsening (56->115) is a real mixing effect, not a bug (median tau agrees
across arms), with in-repo precedent (the reverted BCF aVariance-refresh,
bcf-ridge-interweaving.md). The exact-IG landing did NOT absorb the
slice-era 25x (the review's 3b showed slice ~ exactIG mixing), so the 25x->3x gap
is a DGP/estimator difference, as stated.
