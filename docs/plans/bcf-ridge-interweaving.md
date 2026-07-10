# BCF ridge interweaving / rescale move: derivation

Read-only derivation memo for the `bcf-ridge-interweaving` implementer
(remedy for `bcf-glue-sigma-calibration`, H-MIX confirmed 2026-07-10;
sbc-calibration.md A4e). Repo dbarts @ bartcore cc154b4. Format mirrors
docs/plans/correctness-audit.md derivation blocks: assumptions, algebra,
result, checks. Every model fact is cited file:line. Uncertainty is
flagged inline with [FLAG].

The move: after the conjugate `a` draw and the per-tree prognostic leaf
draws, jointly rescale `(a, mu_1..mu_L) -> (a/c, c*mu_1..c*mu_L)`, `c > 0`,
`c` from its exact full conditional. It travels the likelihood ridge
`a*mu(x) = (a/c)*(c*mu(x))`, collapsing the slow (a, mu-amplitude) mode
whose IACT scales to thousands of sweeps at large `|a0|` (A4e).

---

## 1. Model facts used (all cited)

### 1.1 Prognostic leaf prior (constant Gaussian leaf)

model.hpp:95-138. Each occupied bottom node value:

    mu_l ~ N(0, (scale/k)^2)   iid across leaves,     (model.hpp:95,137-138)

with `scale = nodeScale/sqrt(numTrees)`. For the BCF prognostic forest
`buildBCFForest` sets `nodeScale = s = scaledResponseSd()`
(chain.hpp:474, 1958-1969), `k = 1`, `updateK = false`
(chain.hpp:1987,1985), so

    leafVar := (scale/k)^2 = s^2 / T_mu,   T_mu = # prognostic trees.

Define the leaf prior precision `P_leaf := (k/scale)^2 = 1/leafVar`
(model.hpp:112,129). Keep the general `(k/scale)^2` form in code so the
move stays correct if k ever varies; for shipped BCF `k = 1`.

The prognostic total is `mu(x_i) = sum_l mu_l 1[i in leaf l]`, stored as
`forests_[0].totalFits` (chain.hpp:2030,505-508); each tree's per-obs
leaf-value fits live in `forests_[0].treeFits` (chain.hpp:270,615).

Empty leaves are forced to value 0 and are NOT prior draws
(chain.hpp:1884-1895, "a forced-zero empty leaf is not a draw from the
k-scaled prior"). So the free leaf parameters are the OCCUPIED bottom
nodes only. Define

    L = # occupied prognostic bottom nodes (over all T_mu trees),
    M = sum over occupied prognostic leaves of (leaf value)^2.

These are exactly what `kSumSquaredParams`/`kNumLeaves` accumulate
(chain.hpp:1892-1894) -- but that accumulation is gated on
`forest.updateK`, which is false for BCF, so the implementer must
recompute L and M unconditionally (see section 4).

### 1.2 The a prior via the inverse-gamma auxiliary

drawGlue, chain.hpp:2042-2066; audit block 6, correctness-audit.md:345-350.
Scale-mixture representation of the half-Cauchy:

    aVariance ~ IG(1/2, aPriorScale^2/2)      (prior on the auxiliary)
    a | aVariance ~ N(0, aVariance)            (chain.hpp:2049 uses 1/aVariance as the a-prior precision)
    => marginally a ~ Cauchy(0, aPriorScale)   (Monte-Carlo verified, audit:348-350)

Conditional refresh after the a-draw (chain.hpp:2060-2065):

    aVariance | a ~ IG(1, (a^2 + aPriorScale^2)/2).

The conjugate a-draw itself (mu held fixed, chain.hpp:2048-2058):
`a ~ N(aNum/aPrec, 1/aPrec)`, `aPrec = 1/aVariance + sum_i w_i mu_i^2/sigma^2`,
`aNum = sum_i w_i mu_i (y_i - b_{z_i} tau_i)/sigma^2`.

### 1.3 Likelihood and the combining response

    y_i = a*mu(x_i) + b_{z_i}*tau(x_i) + eps_i,  eps ~ N(0, sigma^2/w_i)

(chain.hpp:2026-2036 combined fit; bcf.md:11). `b0,b1 ~ N(0, bPriorVariance)`
(default 1/2, chain.hpp:2068-2080). The tau forest and b0/b1 are NOT
touched by the move.

### 1.4 Sweep placement

run() loop chain.hpp:579-709: forest loop (mu backfit, then tau backfit)
-> refreshLatents (678) -> drawSigma (688, on the OLD combined fit) ->
`drawGlue(y, weights)` (690) -> k/DART (692-706) -> `storeSample` (708).
Inside the forest loop, `storeSavedTreeRecord` (651) flattens each mu tree
during a recorded sweep, BEFORE drawGlue. This ordering matters for
keepTrees (section 4).

The move sits at the END of drawGlue (after the a, aVariance, b0, b1 draws)
or immediately after it returns -- BEFORE storeSample.

---

## 2. Derivation of the c full conditional

### 2.1 The formally correct framing: a one-coordinate Gibbs update in a
reparameterization (proof of invariance)

Reparameterize the L+1 prognostic-scale coordinates `(a, mu_1..mu_L)` by
the L likelihood invariants and one scale coordinate. The map

    psi_l := a * mu_l   (l = 1..L),        a := a

is a bijection of `(a, mu_{1:L})` with `a != 0`: inverse `mu_l = psi_l/a`.
The Jacobian of `(a, mu_{1:L}) -> (a, psi_{1:L})`:

    d(a, mu)/d(a, psi):  row a = [1, 0..0];  row mu_l = [-psi_l/a^2, (1/a) delta_{lm}]
    |det| = |a|^{-L}    (lower block-triangular, det = 1 * det((1/a) I_L))

so Lebesgue measure transforms as `da dmu_{1:L} = |a|^{-L} da dpsi_{1:L}`.

Because the Gaussian likelihood depends on `(a, mu)` only through the
combined fit `a*mu(x_i) = sum_l psi_l 1[i in leaf l]`, it is a function of
`psi_{1:L}` alone -- CONSTANT along each orbit `{a varies, psi fixed}`.
The move `(a, mu) -> (a/c, c mu)` keeps `psi_l = a mu_l = (a/c)(c mu_l)`
fixed and changes only `a -> a/c`. So in `(a, psi)` coordinates the move
updates ONLY the coordinate `a`, holding `psi` and every other parameter
fixed. That is a legitimate Gibbs sub-block update; it preserves the
posterior iff `a` is drawn from its full conditional given `(psi, rest)`.

The full conditional of `a` given `psi` and everything else (condition on
the auxiliary aVariance -- see 2.3):

    q(a | psi, rest) ∝ pi(a, mu = psi/a, rest) * |a|^{-L}       [Jacobian]
                     ∝ [likelihood: fn of psi only, DROP]
                       * exp(-a^2/(2 aVariance))                [a-prior, 1.2]
                       * prod_l N(mu_l = psi_l/a; 0, leafVar)   [leaf prior, 1.1]
                       * |a|^{-L}.

With `sum_l psi_l^2 = a^2 * M` (M evaluated at the current point, invariant
along the orbit), the leaf-prior factor is
`exp(-(sum_l psi_l^2)/(2 leafVar a^2)) = exp(-a^2 M/(2 leafVar a^2))`... no:
keep psi fixed, `sum psi_l^2 = S_psi` constant, so it is
`exp(-S_psi/(2 leafVar a^2))`. Hence

    q(a | psi) ∝ exp(-a^2/(2 aVariance)) * exp(-S_psi/(2 leafVar a^2)) * |a|^{-L}.  (*)

This is symmetric in `a -> -a`. The full both-sign conditional would allow
a sign flip (which also flips mu, keeping a*mu invariant); since raw a is
sign-ill-posed anyway (A4/A4b), we RESTRICT to the current sign, i.e.
`c > 0`. Restricting to `sign(a)` fixed is conditioning on an extra
(ancillary) discrete coordinate the move holds constant -- still a valid
Gibbs step. Result: the sign of a NEVER flips through the move. [checks 5a]

### 2.2 In terms of c (the operational form)

Set `a = a0/c`, `c > 0`, `a0` = current a. Then `S_psi = a0^2 M`, and
`da = -(a0/c^2) dc`. Substituting into (*) and folding constant `|a0|`
powers into the normalizer:

    exp(-a^2/(2 aVariance))       -> exp(-a0^2/(2 aVariance c^2))
    exp(-S_psi/(2 leafVar a^2))   -> exp(-M c^2/(2 leafVar))
    |a|^{-L}                      -> |a0|^{-L} c^{L}
    |da/dc|                       -> |a0| c^{-2}

    q_c(c) ∝ c^{L-2} * exp( - M c^2 / (2 leafVar)  -  a0^2 / (2 aVariance c^2) ).   (**)

### 2.3 Result: c^2 is Generalized Inverse Gaussian

Substitute `v = c^2` (`dc = dv/(2 v^{1/2})`, `c^{L-2} dc = (1/2) v^{(L-3)/2} dv`):

    q_v(v) ∝ v^{(L-3)/2} * exp( -alpha v - beta/v ),
      alpha = M/(2 leafVar),   beta = a0^2/(2 aVariance).

Matching the GIG density `∝ v^{p-1} exp(-(A v + B/v)/2)`:

    -----------------------------------------------------------------
    v = c^2  ~  GIG( p = (L-1)/2,   A = M/leafVar,   B = a0^2/aVariance )
    -----------------------------------------------------------------
      p = (L-1)/2
      A = M / leafVar = M * (k/scale)^2   (= M * T_mu / s^2 for BCF, k=1)
      B = a0^2 / aVariance
    then  c = sqrt(v),  a <- a0/c,  mu_l <- c * mu_l.

GIG density convention: `f(v) ∝ v^{p-1} exp(-(A v + B/v)/2)`, `v>0`. Sample
with any standard GIG generator (Devroye 2014 ratio-of-uniforms / Hormann-
Leydold; `GIGrvg::rgig` in R, `Rf_*`-free C ports exist). Do NOT
special-case: a single GIG draw covers every regime below.

### 2.4 Condition on aVariance -- do NOT integrate it out

Conditioning on the live auxiliary `aVariance` gives the clean GIG above
and is EXACT, not an approximation: the augmented sampler already targets
the joint `(a, aVariance, mu, ...)` whose a-marginal is Cauchy (1.2), and
an exact Gibbs move on c given aVariance preserves that joint, hence the
Cauchy marginal. This is consistent with audit block 6, which verified all
glue conditionals conditional on aVariance.

Integrating aVariance out (marginal Cauchy prior `p(a) ∝ 1/(a^2 +
aPriorScale^2)`) gives instead

    q_c(c) ∝ c^{L} exp(-M c^2/(2 leafVar)) / (a0^2 + aPriorScale^2 c^2),

which is NOT a standard family (the denominator breaks conjugacy) and would
need slice/griddy sampling. [FLAG] Not tractable in closed form; REJECT the
integrated form. Use the GIG (condition on aVariance).

---

## 3. Edge cases

- **a0 near 0** (`B -> 0`): GIG(p, A, B->0) -> Gamma(shape p, rate A/2),
  i.e. `v = c^2 ~ Gamma((L-1)/2, rate = M/(2 leafVar))`. Proper iff
  `p = (L-1)/2 > 0`, i.e. `L >= 2`. Well-behaved; the standard GIG
  generator handles B=0. If `a0 == 0` exactly then `a <- 0/c = 0` stays 0
  (the move cannot restore a nonzero a -- the scale group cannot cross 0);
  it only rescales an unidentified mu. a0=0 is measure zero.
- **Empty prognostic forest / all-stumps**: L = T_mu (one occupied leaf per
  tree). Formula holds, `p = (T_mu-1)/2`. If every prognostic leaf were
  empty (mu ident. 0) then `L = 0, M = 0` -> skip the move (nothing to
  rescale, a unidentified). Never happens with T_mu=200 default.
- **L = 1 total** (only if T_mu=1 with a single occupied leaf): `p = 0`.
  GIG(0, A, B) is proper for `A>0, B>0`; if additionally `a0=0` (B=0) it is
  improper. BCF T_mu default is 200, so `L >= ~200`; theoretical only.
  [FLAG] Guard `L >= 2 && M > 0`; else skip the move that sweep (rare,
  harmless -- the missed move only forgoes mixing, never corrupts state).
- **1e-9 multiplier floor** (formForestResponse, chain.hpp:2020): floors
  `|a|` only when forming the mu backfit response `resid/m`, `w*m^2`. It
  does NOT enter the move (the move uses exact a, mu, aVariance; no division
  by m). If the move ever drove `|a| < 1e-9`, the NEXT sweep's mu backfit
  floors it -- consistent with existing behavior, benign (audit:353-355
  verified the floor benign for the sufficient statistics). The move draws
  from the exact conditional, which self-regulates against pathological
  shrinkage (the `exp(-M c^2/(2 leafVar))` term penalizes large c).
- **Weights**: the Gaussian likelihood (hence w and sigma) is EXACTLY
  invariant on the orbit, so w, sigma, y, tau, b0, b1 ALL drop out of (**).
  The c-conditional depends ONLY on `M, leafVar, a0, aVariance`. Weighted
  and unweighted BCF use the identical move. [This is the strongest
  self-consistency signal -- section 5a.]
- **updateA gating**: the move CHANGES a; if `bcf_->updateA == false` (a
  pinned, e.g. the bcf-exact mode-1 control a=1), the move must NOT run
  (it would break the pin, and with a fixed the ridge is absent). Gate the
  move on `bcf_->updateA` (chain.hpp:2048's own guard). updateB is
  irrelevant to the move.

---

## 4. State that must be rescaled for consistency

Enumerated against chain.hpp data structures. After drawing c (a0 = pre-move
`bcf_->a`):

Required:
- `bcf_->a *= 1.0/c`                        (chain.hpp:234,324)
- `forests_[0].treeFits *= c`  over all `n*T_mu` entries (chain.hpp:270).
  MANDATORY: the next sweep's residual roll reads each `treeFits[t]` slab
  once (chain.hpp:617-630, `oldFits`) before overwriting it; a stale
  unscaled slab desyncs the roll. One contiguous `misc_scalarMultiply` /
  scal over `forests_[0].treeFits.data()`, length `n*forest.numTrees`.
- `forests_[0].totalFits *= c`  over `n` (chain.hpp:271). Read by the roll
  (chain.hpp:619) and by combinedFits/storeSample.

Required when this sweep is recorded AND test data present:
- `forests_[0].totalTestFits *= c` and `forests_[0].currTestFits *= c`
  (n_test each, chain.hpp:271-272). [FLAG] Under BCF the testFits results
  channel is NaN'd / the test surface refused (bcf-testfits-guard, LANDED
  2026-07-08, TODO:180-189), so these arrays are diagnostically dead, but
  scale them anyway to keep STATE self-consistent (cheap, O(n_test)).

Required when this sweep is recorded AND keepTrees is on:
- The mu-forest SAVED tree records for this sweep's slot were flattened at
  chain.hpp:651 with the PRE-move leaf values. After the move, the reported
  `bcf_->a` is `a0/c` while the saved mu leaves are unscaled -> a stored *
  mu_saved = (a0/c)*mu0 desyncs the identified product. [FLAG -- the sharp
  edge.] FIX (choose one): (i) after the move, multiply every FlatNode value
  in this slot's T_mu saved mu trees by c (chain.hpp:1127-1156 slot layout;
  O(L)); or (ii) restructure so storeSavedTreeRecord for the mu forest runs
  AFTER drawGlue. Option (i) is the smaller change. Verify with a keepTrees
  BCF round-trip: predicted `a*mu` from saved trees must match the live
  combined fit.

Held FIXED (not rescaled):
- `bcf_->aVariance`: the move conditions on it; it stays. It was drawn
  `| a0` at chain.hpp:2065; leaving it conditioned on the pre-move a is a
  one-sweep lag, not an error (next sweep's drawGlue refreshes it `| a_new`).
  [Optional, not required] redraw `aVariance | a_new ~ IG(1, (a_new^2 +
  aPriorScale^2)/2)` right after the move to erase the lag -- one gamma draw.
- `bcf_->b0, bcf_->b1`, the tau forest (`forests_[1]` everything), `sigma_`,
  latents: untouched (tau term and likelihood invariant).

Computing L and M (the move's only new reduction): iterate
`forests_[0]`'s trees; for each occupied bottom node take its value
`treeFits[t][indices[node.begin]]` (exactly recoverParametersFromFits,
chain.hpp:1774-1788) and accumulate `M += v*v; L += 1`, skipping empty
nodes -- mirroring the k-accumulator (chain.hpp:1892-1894) but run
UNCONDITIONALLY (not gated on updateK). Cost O(L + total nodes).

Recording note: `storeSample` records `scale*(a*mu + b_z*tau)+shift`
(chain.hpp:2100-2103) which is INVARIANT -> trainingFits, and results.k /
variableCounts (mu-forest diagnostics) are unaffected by the move.

---

## 5. Verification plan

### 5a. Pure-R prototype (RUN -- adversarial check on the algebra)

`/Users/vdorie/.claude/jobs/7fe13675/tmp/proto.R`. Reduction used: since the
likelihood is constant along the orbit, the move preserves the posterior IFF
it preserves the PRIOR's along-orbit conditional. So draw `(a, mu)` from the
prior (`a ~ half-normal(sd sqrt(aVariance))` on the a>0 branch, `mu_l ~
N(0, leafVar)` iid), push through the move, and test the pushed sample still
has the prior law. A wrong GIG parameter shows immediately.

Two independent samplers -- a grid inverse-CDF on the operational form (**)
and a grid inverse-CDF on GIG(p=(L-1)/2, A=M/leafVar, B=a0^2/aVariance) via
v=c^2 -- both preserve the prior:

    L=3 (v0=0.4, aVar=1.5, N=60000):
      q_c:  KS a=0.159  mu1=0.410  M=0.168   (M ~ Gamma(L/2, rate 1/2v0))
      GIG:  KS a=0.978  mu1=0.395  M=0.067
    L=8 (v0=0.4, aVar=1.5, N=60000, seed 7):
      q_c:  KS a=0.782  mu1=0.261  M=0.703
      GIG:  KS a=0.493  mu1=0.696  M=0.401
    combined-fit invariance under the move: max|delta| = 1.8e-15  (~0)

All KS p-values non-significant (cannot reject prior preservation) at two L
values and two independent parameterizations; a>0 preserved (no sign flip);
combined fit invariant to machine precision. The GIG parameterization the
implementer will code (a standard rgig call) is the one validated. This is
the memo's adversarial self-check and it PASSES.

### 5b. Existing gates

- `benchmarks/R/bcf-exact.R`, `bcf-exact-weak.R` MUST STAY EXACT. They match
  posterior EXPECTATIONS to MC tolerance (E[a*mu], E[tau], E[(b1-b0)tau];
  bcf-exact.R:16-18,33-38), all of which the posterior-preserving move
  leaves unchanged -- it only mixes better. Mode 1 (a=1 fixed) has the move
  gated OFF (updateA=false). Mode 2a (a free) matches the IDENTIFIED E[a*mu]
  (not bare E[mu]), still exact. Expect PASS (possibly tighter, since mode 2a
  "mixes well" already, bcf-exact.R:34).
- `benchmarks/R/equivalence.R`: its scenarios (makeScenarios, lines 60-285:
  friedman/probit/weighted/splitprobs/chik/chains/setdata/wtoffset/quants/
  categorical/missing/dart/linear/gp/logistic) include NO BCF scenario, and
  the move only fires under `bcf_`. So every existing equivalence baseline
  stays BITWISE identical -- NO re-record needed. (Confirm with an
  equivalence compare against the current anchor; expect identical.) The
  TODO "anchor re-record" (TODO:174) bites only a BCF-specific snapshot: if a
  BCF scenario is ever added to equivalence it will shift and must be
  recorded fresh. `inst/tinytest/test-bcf.R` uses structural/finiteness
  checks and a save/restore glue round-trip (no hardcoded draw values), so it
  stays PASS.

### 5c. Acceptance test (the fix target)

sbc-calibration.md A4e, the n=200 glue-on config at the CHEAP thin=120
setting. Currently (A4e point A, sbc-calibration.md:615-618):

    thin=120: sigma rank mean 65.8/75, ecdf 0.1127 (band ~0.092), chisqP 0.0002  FLAG

After the move this must PASS: sigma ecdf inside the ~0.092 band and rank
mean ~= 75 (uniform). Repro: `Rscript benchmarks/R/sbc.R bcf 200 150 120`
(runSbcBCF, burn=150 thinned units). Secondary: on long unthinned chains the
strong-signal (|a0| ~ 7-13) sigma IACT should drop 1-2 orders of magnitude
(A4e measured ~2500-6600 sweeps pre-move; TODO:178). Controls (abs.a,
prog1-5, eff1-5) already pass and must remain passing; raw a and (b1-b0) stay
sign-ill-posed by design (A4b) -- the move does not and should not fix those.

---

## 6. Cost

Per-move arithmetic:
- L, M reduction: O(L + total mu nodes), ~a few thousand ops (L ~ few hundred
  occupied leaves at T_mu=200).
- one GIG draw: O(1), ~tens-hundreds of flops.
- state rescale: `forests_[0].treeFits` scal over `n*T_mu` (the dominant
  term), plus `totalFits` O(n), plus O(n_test) and O(L saved) when recording.

The `n*T_mu` fit-slab scal is the cost driver. At the A4 config (n=200,
T_mu=200) that is 40000 doubles = 320 KB streamed (read+write ~640 KB), on
the order of ONE extra backfit fit-write pass. A full BCF sweep streams that
same slab several times (residual roll + setNodeAverages + parameter/fit set,
plus the tau forest and glue). Honest estimate: the STRAIGHTFORWARD two-scal
move is on the order of a few percent (~5-15%) of a BCF sweep at these sizes
-- NOT the ~1% the plan implies. [FLAG: correct the ~1% claim.]

To recover ~<1%: FUSE the `c` multiply into the residual roll that already
streams the slab next sweep -- scale only `totalFits` (O(n)) at move time,
carry `c` as a pending factor, and apply it to `oldFits` inside the existing
roll (chain.hpp:617-630) which reads each `treeFits[t]` once before
overwriting it. That adds one multiply per element to a pass already paid
for, dropping the move's marginal cost to O(n + L) + one GIG draw. The
fusion is delicate (the incremental roll for t>=1 mixes an old, c-needing
`oldFits[t]` with an already-fresh `prevFits[t-1]`); recommend shipping the
simple two-scal version first (correctness), optimizing only if profiling
shows the move on the critical path.

Either way the per-sweep overhead is irrelevant next to the payoff: the move
lets BCF calibrate at thin=120 instead of needing thin >> 1000 for the
Cauchy-tail reps (A4e), a 1-2 order-of-magnitude reduction in sweeps-per-ESS.
Net compute for a target accuracy drops sharply.

---

## 7. Certainty flags summary

- CERTAIN (algebra + prototype-verified): the GIG result in 2.3, the
  invariance proof in 2.1, likelihood/weight/sigma drop-out, sign never
  flips, condition-on-aVariance is exact, integrate-out is intractable.
- FLAG (implementation traps, not algebra):
  * keepTrees under BCF: saved mu leaves flattened pre-move must be scaled
    by c (or reorder) -- section 4, the sharp edge.
  * treeFits MUST be scaled (roll reads it) -- easy to forget, silent
    corruption if missed.
  * gate on `updateA`; guard `L>=2 && M>0`.
  * the ~1% cost claim is optimistic for the simple implementation (~5-15%);
    <1% needs the fused-roll optimization.
- LOW-RISK RESIDUAL: aVariance one-sweep lag (benign; optional refresh);
  the 1e-9 floor never interacts with the move.

## Status

- 2026-07-10: derivation complete (read-only Opus deriver), prototype
  PASSED (prior-preservation at L=3 and L=8 via two independent
  samplers, KS non-significant, combined fit invariant to 1.8e-15).
  Result: c^2 ~ GIG((L-1)/2, M/leafVar, a0^2/aVariance), conditioned
  on the IG auxiliary (exact, audit-block-6 consistent). Two
  corrections to the earlier chat-level cost picture: equivalence.R
  has NO BCF scenario and the move is BCF-gated, so NO anchor
  re-record is needed; the naive treeFits rescale costs ~5-15% of a
  sweep, recovered to <1% by fusing the c multiply into the residual
  roll. Sharp edge: keepTrees flattens saved mu leaves BEFORE
  drawGlue (chain.hpp:651) - the saved slot must be rescaled or the
  ordering changed. VD approved the remedy in principle 2026-07-10;
  sequencing (orchestrator's call): after c-api-growth.
  NOT YET IMPLEMENTED.
