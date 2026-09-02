# BCF b-ridge (treatment-scale) interweaving / rescale move: derivation

Note (added 2026-08-28): the recommendation below still holds - the
treatment-scale ridge is not implemented; `AmplitudeForestCombiner` in
src/bartcore/combiner.hpp defaults `ridgeB = false`. The code citations
below (`BCFState`, `bcf_`, `buildBCFForest`, the hardcoded `forests_[1]`)
predate the multi-forest generalization and no longer resolve; that
machinery is now `AmplitudeForestCombiner` and its `ridgeB` flag in
src/bartcore/combiner.hpp. Treat citations below as historical.

Read-only derivation memo for the `bcf-sigma-residual` follow-up. Successor to
`bcf-ridge-interweaving` (a-ridge, LANDED 9617c94): that move collapsed the
prognostic `(a, mu)` ridge but left the thin=120 sigma SBC flag standing in
BOTH pre- and post-move builds at R=400 (rank mean ~67-69/75). The prime
suspect named on the TODO is the treatment-side scale ridge. Repo dbarts @
bartcore 3bcb664. Format mirrors docs/plans/archive/bcf-ridge-interweaving.md
(the a-move memo): assumptions, algebra, result, checks. Every model fact is
cited file:line. Uncertainty is flagged inline with [FLAG].

The move: after the conjugate `b0, b1` draws in drawGlue (and the a-ridge
move that already ships), jointly rescale the treatment-scale coordinates

    (b0, b1, tau_1..tau_Lt) -> (b0/c, b1/c, c*tau_1..c*tau_Lt),   c > 0,

`c` from its exact full conditional. It travels the likelihood ridge
`b_{z(i)} * tau(x_i) = (b_{z(i)}/c) * (c*tau(x_i))`, collapsing the slow
`(b, tau-amplitude)` mode -- the exact analogue, on the treatment forest, of
the non-identifiability the a-move fixed on the prognostic forest.

---

## 1. Model facts used (all cited, chain.hpp @ 3bcb664)

### 1.1 Treatment (tau) forest leaf prior (constant Gaussian leaf)

The tau forest is `forests_[1]`, built by `buildBCFForest(spec.tau,
sdModerate * s / 0.674)` ([[chain.hpp:475@9cebb352]]). buildBCFForest sets `k = 1`,
`updateK = false`, `useDart = false`, and
`leaf.scale = nodeScale / sqrt(numTrees)` ([[chain.hpp:2079-2083@9cebb352]]). So with
`nodeScale_tau = sdModerate * s / 0.674` (s = scaledResponseSd,
[[chain.hpp:473@9cebb352]]; sdModerate default 1, kHalfNormalMedian = 0.674):

    tau_l ~ N(0, leafVar_tau)  iid across occupied leaves,
    leafVar_tau := (leaf.scale_tau / k)^2 = leaf.scale_tau^2   (k = 1)
                 = (sdModerate * s / 0.674)^2 / T_tau,   T_tau = # tau trees.

Define the tau-leaf precision `P_tau := (k / leaf.scale_tau)^2 =
1/leafVar_tau` (exactly the a-move's `(forest.k/forest.leaf.scale)^2`,
[[chain.hpp:558@9cebb352]], but on forests_[1]). Keep the general `(k/scale)^2` form in
code so the move stays correct if k ever varies; for shipped BCF k = 1.

Empty tau leaves are forced to value 0 and are NOT prior draws
([[chain.hpp:2024@9cebb352]], same rule as mu). So the free tau parameters are the
OCCUPIED bottom nodes only. Define

    Lt = # occupied tau bottom nodes (over all T_tau trees),
    Mt = sum over occupied tau leaves of (leaf value)^2.

The tau totals are `forests_[1].totalFits`; per-tree leaf-value fits live in
`forests_[1].treeFits` (n * T_tau). Recover a leaf's value exactly as the
a-move does: `treeFits[t][indices[node.begin]]` ([[chain.hpp:549@9cebb352]]).

### 1.2 The b prior -- a PLAIN Gaussian, NO auxiliary (contrast a)

drawGlue, [[chain.hpp:2162-2175@9cebb352]]. Unlike a (half-Cauchy via the inverse-gamma
auxiliary aVariance), the treatment coefficients have a fixed-variance normal
prior with NO scale-mixture auxiliary:

    b0 ~ N(0, bPriorVariance),  b1 ~ N(0, bPriorVariance),  independent,
    bPriorVariance = 0.5  (fixed; BCFState [[chain.hpp:313@9cebb352]], [[chain.hpp:327@9cebb352]]; spec 480).

They are drawn JOINTLY in drawGlue (one `if (bcf_->updateB)` block,
[[chain.hpp:2162@9cebb352]]): each is an independent conjugate normal from its own arm's
rows -- b0 from control rows (z=0), b1 from treated rows (z=1):

    b0 ~ N(n0/p0, 1/p0),  p0 = 1/bPriorVariance + sum_{z_i=0} w_i tau_i^2/sigma^2,
                          n0 = sum_{z_i=0} w_i tau_i (y_i - a mu_i)/sigma^2;
    b1 ~ N(n1/p1, 1/p1),  same over treated rows ([[chain.hpp:2163-2174@9cebb352]]).

Because bPriorVariance is a FIXED constant, there is NO auxiliary to condition
on and no aVariance-style one-sweep lag: the b-move's conditional is exact and
unconditional. This is strictly SIMPLER than the a-move -- the single biggest
structural difference between the two derivations.

### 1.3 Likelihood and the combining response

    y_i = a*mu(x_i) + b_{z_i}*tau(x_i) + eps_i,  eps ~ N(0, sigma^2/w_i)

([[chain.hpp:2127-2128@9cebb352]] combinedFits, 2136-2140 drawGlue; bcf.md). The mu forest
and a are NOT touched by the b-move.

### 1.4 Sweep placement

run() loop: forest loop (mu backfit, tau backfit) -> refreshLatents (769) ->
drawSigma (779, on the OLD combined fit) -> `drawGlue` (782) ->
`interweaveGlueRidge` [a-move] (783) -> storeSample. The b-move sits
immediately after the a-move (both after drawGlue, before storeSample). Under
keepTrees, `storeSavedTreeRecord` flattens each tau tree during a recorded
sweep BEFORE drawGlue -- same sharp edge as the a-move (section 4).

---

## 2. Derivation of the c full conditional

### 2.1 One-coordinate Gibbs update in a reparameterization (invariance proof)

Reparameterize the Lt+2 treatment-scale coordinates `(b0, b1, tau_{1:Lt})` by
their Lt+1 likelihood invariants plus one scale coordinate. Assume `b0 != 0`
(measure zero otherwise; handled in 3) and use the map

    r := b1/b0,   psi_l := b0 * tau_l  (l=1..Lt),   b0 := b0,

a bijection of `(b0, b1, tau_{1:Lt})`. Inverse: `b1 = r b0`,
`tau_l = psi_l/b0`. Jacobian of `(b0,b1,tau) -> (b0, r, psi)`:

    row b0   = [1,       0,   0..0]
    row b1   = [r,      b0,   0..0]
    row tau_l= [-psi_l/b0^2, 0, (1/b0) delta_{lm}]

block lower-triangular, so
`|det d(b0,b1,tau)/d(b0,r,psi)| = 1 * b0 * (1/b0)^Lt = |b0|^{1-Lt}`, i.e.

    db0 db1 dtau_{1:Lt} = |b0|^{1-Lt} db0 dr dpsi_{1:Lt}.

NOTE the extra `+1` in the exponent relative to the a-move's `|a|^{-L}`: the
second glue scalar `b1 = r b0` contributes one factor of `b0`. This is what
distinguishes the two-scalar move from the one-scalar a-move.

The Gaussian likelihood depends on `(b0, b1, tau)` only through the combined
fit; for the treatment block that means only through `b0 tau_l = psi_l` (over
control rows in leaf l) and `b1 tau_l = r psi_l` (over treated rows) -- both
functions of `(r, psi)` ALONE, constant along each orbit `{b0 varies, (r,psi)
fixed}`. The move `(b0,b1,tau)->(b0/c, b1/c, c tau)` keeps `r = b1/b0` and
`psi_l = b0 tau_l = (b0/c)(c tau_l)` fixed and changes only `b0 -> b0/c`. So
in `(b0, r, psi)` coordinates the move updates ONLY the coordinate `b0`. That
is a legitimate Gibbs sub-block update; it preserves the posterior iff `b0` is
drawn from its full conditional given `(r, psi, rest)`.

Full conditional of `b0` given `(r, psi, rest)`:

    q(b0 | r,psi,rest) prop pi(b0, b1=r b0, tau=psi/b0, rest) * |b0|^{1-Lt}   [Jacobian]
      prop [likelihood: fn of (r,psi) only, DROP]
           * exp(-b0^2/(2 sB2))                       [b0 prior, 1.2]
           * exp(-(r b0)^2/(2 sB2))                   [b1 prior]
           * prod_l N(psi_l/b0; 0, leafVar_tau)       [tau prior, 1.1]
           * |b0|^{1-Lt}.

Combine the two b priors: `exp(-(1+r^2) b0^2/(2 sB2))` and note
`(1+r^2) b0^2 = b0^2 + b1^2`. The tau factor is
`exp(-(sum psi_l^2)/(2 leafVar_tau b0^2))` with `sum psi_l^2 = S_psi` fixed on
the orbit (at the current point `S_psi = b0^2 Mt`). Hence

    q(b0) prop exp(-(b0^2+b1^2)/(2 sB2)) * exp(-S_psi/(2 leafVar_tau b0^2)) * |b0|^{1-Lt}.  (*)

Symmetric in `b0 -> -b0`; RESTRICT to the current sign (c > 0), so neither b0
nor b1 nor any tau leaf ever flips sign -- a valid Gibbs step conditioning on
the ancillary sign coordinate the move holds fixed. `(b1-b0)` stays
sign-ill-posed by design (A4b), as for the a-move. [checks 5a]

### 2.2 In terms of c (operational form) -- MOVED

Moved to docs/design/multiplier-combiner.md, "The exponent rule".

### 2.3 Result: c^2 is Generalized Inverse Gaussian -- MOVED

Moved to docs/design/multiplier-combiner.md, "The exponent rule" (section 4
below still uses the result: `gigP = (Lt-2)/2`, `gigB =
(b0^2+b1^2)/bPriorVariance`).

### 2.4 No aVariance analogue -- nothing to condition on or refresh

The b prior is a fixed-variance normal (1.2), so there is no scale-mixture
auxiliary. B is a deterministic function of the current (b0, b1) and the fixed
bPriorVariance. Hence the b-move has NONE of the a-move's aVariance subtleties:
no "condition on the live auxiliary", no optional refresh, no one-sweep lag.
The b-conditional is exact outright. [This removes the a-move's most delicate
correctness argument.]

---

## 3. Edge cases

- **b0 near 0** (part of `B -> 0` only if ALSO b1 near 0): `B = (b0^2+b1^2)/sB2`.
  Since b0, b1 have continuous priors and (all-treated aside) at least one arm
  has a likelihood term, `B > 0` a.s. If b0=0 exactly, `b0 <- 0/c = 0` stays 0
  (the scale group cannot cross 0), harmless; measure zero. With `B>0, A>0` the
  GIG is proper for EVERY real p, including negative -- so unlike the a-move
  there is no propriety worry at small Lt.
- **All-treated or all-control data**: if every z=1 (no control rows), b0 has
  NO likelihood term; drawGlue draws it from its prior `N(0, sB2)`
  (p0=1/sB2, n0=0, [[chain.hpp:2173@9cebb352]]). b0 then multiplies tau only over control
  rows -- of which there are none -- so it does NOT appear in the likelihood,
  and the move `b0 -> b0/c` is STILL likelihood-invariant (only `b1 tau`
  matters, `-> (b1/c)(c tau) = b1 tau`). The derivation never assumed both arms
  occupied: b0 remains a genuine `N(0,sB2)` parameter and its prior term
  `b0^2/sB2` belongs in B exactly as written. Including a prior-only b0 in B is
  correct (the joint is preserved by construction); it is a benign passenger.
  Prototype (5a) confirms combined-fit invariance in the all-treated case to
  1e-16. Symmetric for all-control.  ==> the move is UNCHANGED for one-arm data.
- **Empty / one-leaf tau forest**: `Lt = T_tau` for all-stumps. If every tau
  leaf were empty (`Lt=0, Mt=0`) skip (nothing to rescale). `Lt=1` gives
  `p=-1/2`; GIG(-1/2, A>0, B>0) is proper (an inverse-Gaussian relative), so
  the move is valid, but mirror the a-move and GUARD `Lt >= 2 && Mt > 0` --
  T_tau default is 50 so `Lt >= ~50`; the guard is theoretical and a skip only
  forgoes mixing, never corrupts state. [FLAG: the guard threshold `Lt>=2` is a
  convention copied from the a-move, not a propriety requirement here (B>0
  makes every Lt proper); `Lt>=1` would also be safe. Keep `Lt>=2` for
  symmetry with interweaveGlueRidge.]
- **1e-9 multiplier floor** (formForestResponse, [[chain.hpp:2114@9cebb352]]): floors
  `|b0|`, `|b1|` only when forming the tau backfit response `resid/m`, `w m^2`.
  It does NOT enter the move (which uses exact b0,b1,tau; no division by m). If
  the move drove `|b|<1e-9`, the next sweep's tau backfit floors it -- benign,
  consistent with existing behaviour. The `exp(-Mt c^2/(2 leafVar_tau))` term
  self-regulates against pathological shrinkage.
- **Weights / sigma**: the Gaussian likelihood is EXACTLY invariant on the
  orbit, so w, sigma, y, a, mu ALL drop out of (**). The c-conditional depends
  ONLY on `Mt, leafVar_tau, b0, b1, bPriorVariance`. Weighted and unweighted
  BCF use the identical move.
- **updateB gating**: the move CHANGES b0, b1; gate on `bcf_->updateB` (the
  drawGlue-B guard, [[chain.hpp:2162@9cebb352]]). If b is pinned the ridge is absent and the
  move must not run. `updateA` is irrelevant to the b-move (they act on
  disjoint blocks -- section 3a).

### 3a. Interaction with the LANDED a-move (order, commutation)

The a-move rescales `(a, forests_[0])`; its conditional reads `(M_mu,
leafVar_mu, a, aVariance)`. The b-move rescales `(b0, b1, forests_[1])`; its
conditional reads `(Mt, leafVar_tau, b0, b1, bPriorVariance)`. The two blocks
are DISJOINT: neither move touches any input of the other (the a-move does not
change b0/b1/tau; the b-move does not change a/mu/aVariance). Therefore the two
moves COMMUTE and are each an exact Gibbs update conditional on the rest. Run
BOTH every sweep, in either order, immediately after drawGlue:

    drawGlue(y, weights);
    interweaveGlueRidge(record, sampleNum);     // a-move (mu, a)  -- ships
    interweaveTreatmentRidge(record, sampleNum);// b-move (tau, b0,b1) -- NEW

No need to alternate or share machinery beyond the GIG RNG. Order-independence
is worth a cheap test (run both orders on a fixed seed; combined fit and the
joint law must agree). [FLAG: the ONLY shared state is the rng stream -- the
two GIG draws consume rng in call order, so a fixed-seed bitwise gate must fix
the order a<-b; this is a stream-ordering convention, not a correctness issue.]

---

## 4. State that must be rescaled for consistency

Enumerated against chain.hpp, exactly mirroring interweaveGlueRidge but on
`forests_[1]` and both b scalars. After drawing c (b0_0,b1_0 = pre-move):

Required:
- `bcf_->b0 *= 1.0/c;  bcf_->b1 *= 1.0/c`         ([[chain.hpp:324@9cebb352]])
- `forests_[1].treeFits *= c` over all `n*T_tau` entries -- MANDATORY: next
  sweep's residual roll reads each treeFits slab once before overwriting; a
  stale slab desyncs the roll. One `misc_scalarMultiplyVectorInPlace` over
  `forests_[1].treeFits.data()`, length `n*T_tau` ([[chain.hpp:572@9cebb352]] analogue).
- `forests_[1].totalFits *= c` over `n` (read by the roll and combinedFits).

Required when recorded AND test data present:
- `forests_[1].totalTestFits *= c`, `forests_[1].currTestFits *= c`
  (n_test each). [FLAG] testFits is diagnostically dead under BCF
  (bcf-testfits-guard) but scale for STATE self-consistency, cheap.

Required when recorded AND keepTrees on (THE SHARP EDGE, identical to a-move):
- this sweep's saved TAU tree slot was flattened with PRE-move leaf values; after
  the move multiply every FlatNode leaf value in the slot's T_tau saved tau
  trees by c ([[chain.hpp:588-597@9cebb352]] analogue, `node.variable == invalidVariable`),
  so `b_saved * tau_saved` keeps the identified product. Verify with a keepTrees
  BCF round-trip: predicted `b_z*tau` from saved trees must track the live fit.

Held FIXED (not rescaled):
- `bcf_->a`, `aVariance`, the mu forest (`forests_[0]` everything), `sigma_`,
  latents: untouched (invariant, and disjoint from the b block). There is NO
  auxiliary to refresh (contrast the a-move's aVariance note).

Computing Lt and Mt: iterate `forests_[1]`'s trees; for each occupied bottom
node take `treeFits[t][indices[node.begin]]`, accumulate `Mt += v*v; ++Lt`,
skipping empty nodes -- run UNCONDITIONALLY (not gated on updateK, which BCF
leaves false), exactly as interweaveGlueRidge does for the mu forest
([[chain.hpp:541-553@9cebb352]]). Cost O(Lt + total tau nodes).

Recording note: storeSample records `scale*(a mu + b_z tau)+shift`
([[chain.hpp:2195-2197@9cebb352]]) which is INVARIANT -> trainingFits unaffected.

The implementation is interweaveGlueRidge with forests_[0] -> forests_[1],
the a-block -> the (b0,b1)-block, aVariance -> (dropped), gigP (L-1)/2 ->
(Lt-2)/2, gigB a^2/aVariance -> (b0^2+b1^2)/bPriorVariance. Factor the shared
scaffolding into one templated helper taking the forest index, the glue
scalar(s), and (p,A,B).

---

## 5. Verification plan

### 5a. Pure-R prototype (RUN -- adversarial check on the algebra) -- PASSED, MOVED

Moved to docs/design/multiplier-combiner.md, "The exponent rule": the
`proto-b.R` run (N=20000, two parameterizations, three `Lt`), its DISCRIMINATION
block, and the KS 1.6e-21 rejection of the off-by-one exponent that CONFIRMS
`p=(Lt-2)/2`. PASSES stands; 5b-5d below read it as established.

### 5b. Existing gates (expected, by the same argument as the a-move)

- `bcf-exact.R`, `bcf-exact-weak.R` MUST STAY EXACT: they match posterior
  EXPECTATIONS (E[a mu], E[tau], E[(b1-b0)tau]) which a posterior-preserving
  move leaves unchanged. The three modes ([[bcf-exact.R:16-18@9cebb352]], 379-413) are the
  clean gate:
    mode 1  (updateA=F, updateB=F): BOTH moves gated off -- unchanged.
    mode 2a (updateA=T, updateB=F): a-move active, b-move gated off -- unchanged
            by the b-move; matches E[a mu], E[tau].
    mode 2b (updateA=F, updateB=T): b-move ACTIVE (the b-move's headline gate) --
            must keep E[mu] and E[(b1-b0)tau] exact. It will (posterior-
            preserving); mode 2b is exactly to the b-move what mode 2a is to the
            a-move. The landing note that mode 2b "mixes poorly / b trades off
            with tau scale" ([[bcf-exact.R:34@9cebb352]]) is the very ridge this move
            collapses -- expect mode 2b to match AND mix better.
- `equivalence.R`: no BCF scenario and the move is bcf_-gated -> every existing
  baseline stays BITWISE identical, NO re-record. Confirm with a compare
  against the current anchor.
- `inst/tinytest/test-bcf.R`: structural/finiteness + glue round-trip, no
  hardcoded draws -> PASS; add a move-active + keepTrees smoke pair mirroring
  the a-move's 7 checks.

### 5c. Acceptance test (the fix target)

The n=200 glue-on config at thin=120 (L=150, burn=150), pooled to R=400. The
a-move landing measured BOTH builds FLAG sigma there:

    build   R    sigma mean  ecdf    band    verdict
    new(a)  400   67.6/75     0.1003  0.0656   FLAG   (a-move only)

The ACCEPTANCE for the b-move: with BOTH the a- and b-ridge moves active,
sigma at thin=120, R>=400 must PASS -- rank mean ~= 75 (i.e. ~L/2 on the 0..L
scale), sigma ecdfDiff inside the simultaneous band, chisqP not significant.
Controls (abs.a, abs.diff, prog1-5, eff1-5) stay PASS. Raw a and (b1-b0) stay
sign-ill-posed by design. IF the b-move alone still does not clear sigma, the
localization below (5d) says whether a burn-transient fix is additionally
required.

### 5d. Localization controls (RUN -- see section 6)

Two discriminating SBC controls decide whether the b-ridge is actually the
residual's carrier before anyone implements the move:
  (a) fixed-b: updateA=TRUE (a-move active), updateB=FALSE (b pinned) -- if
      sigma CALIBRATES, the b-path is implicated (implement the b-move); if it
      still FLAGS, the b-ridge is exonerated.
  (b) burn-2x: standard full glue, burn doubled in absolute sweeps at fixed
      thin -- if the flag fades, the residual is a burn TRANSIENT (cheap
      remedy: longer default burn / better init), not another ridge.

---

## 6. Localization results

Both controls: n=200 glue-on config, thin=120, L=150, pooled 4x50 = R=200,
current shared library (a-move LANDED 9617c94, b-move NOT built). Table format
= A4e (rankMean out of L=150 so target 75; ecdf vs simultaneous band; chisqP;
verdict). Chunk .rds + run-log under .../tmp/b-ridge/.

Reference points (from the a-move landing, thin=120, sigma):
    full glue, a-move active, R=200:  sigma 66.2/75  ecdf 0.1013  band 0.0924  FLAG
    full glue, a-move active, R=400:  sigma 67.6/75  ecdf 0.1003  band 0.0656  FLAG

### Control (a) -- FIXED-b  (updateA=TRUE a-move active, updateB=FALSE b pinned)

    func      rankMean   ecdf    band    chisqP  verdict
    sigma       61.69   0.1359  0.0924   0.0002  FLAG
    prog1..5    72.9-75.9  <=0.065  0.0924  ...   PASS (all)
    eff1..5     69.8-77.9  <=0.066  0.0924  ...   PASS (all)
    (mapR2 = 1.00000; wall 385s)

VERDICT: pinning the treatment b-ridge does NOT calibrate sigma -- it FLAGS
HARDER (rank 61.7 vs the full config's 66.2; ecdf 0.136 vs 0.101, further out
of band). Direction: low rank mean => sig0 (truth) sits below most posterior
sigma draws => posterior sigma biased HIGH (residual over-estimated), i.e. the
forests/glue UNDER-explain the signal; with b frozen there is even less
capacity to fit the treatment arm, so more variance leaks to sigma. This
EXONERATES the b-ridge as the carrier of the residual sigma bias: the a-only
memo's "prime suspect" is WRONG. [Caveat: fixed-b (b0=0) is a slightly
different DGP -- the control arm carries no tau -- so the magnitude is not
directly comparable; the QUALITATIVE result (still flags, harder) is robust.]

### Control (b) -- BURN-2x  (full glue, burn 150->300 thinned units, same thin)

    func         rankMean   ecdf    band    chisqP  verdict
    sigma          72.05   0.0888  0.0924   0.0020  PASS   (<- was 66.2 / 0.1013 FLAG at burn=150)
    abs.a          72.42   0.0570  0.0924   0.1359  PASS
    abs.diff       75.81   0.0447  0.0924   0.5493  PASS
    a              75.30   0.2419  0.0924   0.0000  FLAG   (sign-ill-posed by design, A4b)
    b1.minus.b0    74.56   0.1484  0.0924   0.0000  FLAG   (sign-ill-posed by design, A4b)
    prog1..5       75.6-78.5  <=0.079  0.0924  ...    PASS (all)
    eff1..4        77.4-82.8  <=0.092  0.0924  ...    PASS (all)
    eff5           82.03   0.0926  0.0924   0.3593  FLAG   (a hair over band, chisqP=0.36 -> noise)
    (mapR2 = 1.00000; wall 1150s, contended with the 5-core grid)

VERDICT: doubling burn (150->300 thinned units = 18000->36000 absolute sweeps
at thin=120), holding thin fixed, moves sigma from FLAG (66.2/75, ecdf 0.1013)
to PASS (72.05/75, ecdf 0.0888, just inside the R=200 band) -- rank mean
climbs toward 75, the low bias shrinks. The residual sigma bias is a BURN
TRANSIENT, confirmed. (The a and b1-b0 FLAGs are the DESIGNED sign-ill-posed
raw glue -- their abs.* identified counterparts PASS. eff5 sits 0.0002 over the
band with chisqP 0.36 -- multiple-comparison noise across 10 spatial
functionals, not a real flag.) This also RE-READS A4e: its "PASS at thin=1000"
was really a pass at longer ABSOLUTE burn (burn=150 thinned * 1000 = 150k
sweeps), the confound the TODO named -- longer thin PASSED because it dragged
burn along, not because thinning per se was needed.

## 7. Routed recommendation

Both controls point the same way and AWAY from the b-move:

  * Control (a): pinning the treatment b-ridge does NOT calibrate sigma
    (worsens it) -> the b-ridge is NOT the residual's carrier.
  * Control (b): doubling burn alone calibrates sigma -> the residual is a
    BURN TRANSIENT (a fraction of reps, the Cauchy-tail strong-|a0| minority,
    have not reached stationarity in burn=150; the under-fit signal inflates
    sigma). Consistent with sigma IACT ~40 << thin=120 (draws near-independent,
    so the bias is a transient offset, not autocorrelation).

RECOMMENDATION -- ROUTE: FIX INITIALIZATION / BURN, do NOT implement the
b-move for the sigma flag.

  1. bcf-sigma-residual (the sigma calibration defect): remedy is a longer BCF
     burn-in default and/or an amplitude-aware initialization, NOT another
     ridge move. Cheapest: raise the BCF default burn-in. Better (targets the
     mechanism): initialize `a` (and the aVariance auxiliary) near the
     data-implied prognostic scale so the strong-|a0| reps do not spend the
     transient climbing from a=1; the a-ridge move then keeps them mixed.
     ACCEPTANCE for the fix: the n=200 glue-on config at thin=120, pooled to
     R>=400 (the tighter band ~0.0656), must show sigma ecdfDiff INSIDE the
     band and rank mean ~75, with abs.a / abs.diff / prog* / eff* still PASS.
     CAUTION: burn=300 cleared sigma at R=200 (band 0.0924, margin 0.0036) but
     ecdf 0.0888 would still EXCEED the R=400 band (~0.0656) -- so burn=300 is a
     demonstrated LOWER BOUND, not a proven sufficient default. The implementer
     must tune burn upward (or add the init fix, which should need less burn)
     until sigma ecdf clears the R=400 band; verify at R>=400.

  2. The b-move itself: KEEP the derivation on the shelf as a SEPARATE,
     independently-justified mixing improvement for the treatment ridge
     (exactly what the a-move is for the prognostic ridge; [[bcf-exact.R:34@9cebb352]] notes
     mode-2b "mixes poorly, b trades off with tau scale" -- that IS this
     ridge). If implemented, judge it on its OWN IACT payoff (e.g. |b1-b0| or
     tau-amplitude IACT on a strong-treatment-signal DGP, mirroring the a-move's
     |a| IACT 2.5x measurement), NOT on the sigma SBC flag, which it does not
     and need not fix. It is cheap (~4%/sweep, same two-scal structure) and the
     GIG machinery already ships; correctness ACCEPTANCE = bcf-exact mode-2b
     stays exact AND a keepTrees BCF round-trip tracks. Implementation template:
     interweaveGlueRidge with the substitutions in section 4 (forests_[1], the
     (b0,b1) block, p=(Lt-2)/2, B=(b0^2+b1^2)/bPriorVariance, no aVariance).

  BOTH is defensible (fix burn now; land the b-move later as a mixing win), but
  the two are DECOUPLED: the sigma defect is fixed by burn/init, and the b-move
  is a mixing nicety. Do not gate the b-move on sigma, and do not expect the
  b-move to substitute for the burn/init fix.

## 8. Certainty flags summary

- CERTAIN (algebra + prototype-verified): the GIG result in 2.3 with
  p=(Lt-2)/2 (discrimination-tested), the invariance proof in 2.1,
  likelihood/weight/sigma drop-out, sign never flips, all-treated/all-control
  invariance, the b-move commutes with the a-move (disjoint blocks), NO
  auxiliary so the conditional is exact outright.
- FLAG (implementation traps, not algebra):
  * keepTrees under BCF: saved TAU leaves flattened pre-move must be scaled by
    c (or reorder) -- the same sharp edge as the a-move.
  * forests_[1].treeFits MUST be scaled (roll reads it) -- silent corruption if
    missed.
  * gate on updateB; guard Lt>=2 && Mt>0 (convention; B>0 makes any Lt proper).
  * two GIG draws (a then b) consume rng in a fixed order -- pin it for bitwise
    gates.
  * bcf-exact mode-1 has updateB ON -- confirm the b-move keeps E[(b1-b0)tau]
    exact there.
- LOW-RISK RESIDUAL: none of the a-move's aVariance-lag class (no auxiliary).
