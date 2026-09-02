# nameable-calibration

agent: S0 opus (it freezes a header signature whose content is engine
  semantics). S1 opus (engine conversion, three refusal sites, creation
  surface). S2 opus (engine get/set, `facade.hpp` virtuals, bridge; the R5 and
  man-page half rides the slice). S3 opus, folded INTO dbarts-h-reshape S1.
  Serialized: one implementer, each slice lands before the next starts.
rng: NEUTRAL on every default path, EVERY slice, and MEASURED rather than
  asserted. The creation half is byte-identical while `prior.scale` is
  `NA_real_` - the family switch still resolves `node.scale` - and the
  mid-chain half is inert until called and bitwise-SKIPS a write that
  reproduces the current internal scale. DRAW-LAW-CHANGING only on the opt-in
  named path, which is the arc's point: naming a calibration names a different
  prior, and no default moves. Gate: the trio bitwise at every slice; a trio
  deviation is a bug in the slice, NEVER a re-record.
window: pre-release. S3 is an item INSIDE dbarts-h-reshape S1, the second and
  last re-bake of the window, so the flat spelling is decided there or never;
  S0 must therefore land before that slice starts. S1 and S2 are pure
  additions with no header and may land on either side of the reshape.
  Everything pre-release is breakable and the sister packages migrate in
  lockstep at the freeze.
budget: S0 ~50 (no code); S1 ~358; S2 ~360; S3 ~295 INSIDE dbarts-h-reshape
  S1's own budget; tests ~530 across S1 and S2, plus ~125 of test riding the
  reshape's own files inside S3. Total ~1593 in repo, of which ~530 is test.
  Out of repo: 0 required, 2 optional one-liners.
tip: every anchor below read against b2290614 in this worktree.

## Goal

An R program can NAME the leaf-prior calibration it composes against, in
response units, instead of inheriting one from the range of whatever vector
the sampler happened to be constructed on. After the arc: `prior.scale` at
creation (a `dbartsModel` slot plus a `bart2`/`rbart_vi`/`bart` formal),
`$getCalibration` / `$setCalibration` mid-chain over every chain, and two flat
C entries; and the getter is the only authoritative reader of what is in
force, so a re-anchoring channel stops being silent.

(artifacts untracked - memo, critique, synthesis,
recritique)

## The wall this closes

docs/design/r-c-division.md, "The wall that remains is calibration, and it is
measured": a structurally correct pure-R probit composition inherits its leaf
prior from the construction range, a 16x sweep of that accident moves the
implied leaf sd 16x and the posterior sd of `f` by 4.6x with no error and no
warning, and correctness is recovered only near one lucky range. The
composed model's BLOCKING ports to R; its PRIOR does not, unless it can be
named. The principle's clause 4 ("What that asks") states the obligation - R
should be able to name the calibration it composes against, and read/write
symmetry is the default - and the adoption slate entry for this item schedules
the design before the dbarts.h reshape's S1 because it may leave a header
footprint. This plan is that design.

## 1. What the calibration IS

For a forest with `m` trees, leaf scale `s` (internal units), precision
divisor `k`, and response transform `(fitScale, fitShift)`:

    prior.scale = fitScale * s * sqrt(m)         # response units, sd at k = 1
    prior.sd    = prior.scale / k                # response units, at current k
    prior.mean  = fitShift                       # response units

`prior.scale` is the identified, nameable object. Only the ratio `s/k` enters
any draw law - measured bitwise on every path, and the reason `k` is not in
the setter at all - so under a fixed `k` the pair `(prior.scale, k)` has one
degree of freedom and only `prior.scale` is identified; `k` matters when and
only when a hyperprior draws it, and then `prior.scale` is the constant of the
model while `prior.sd` moves every sweep. `prior.sd` is the ergonomic reading
and is exact whenever `k` is fixed. `prior.mean` is a read-only diagnostic
whose lever is the offset channel.

The engine already spells the anchor internally (`[[combiner.hpp:285-286@4c018187]]`, "tau
= nodeScale / k is the per-forest total-fit prior sd"), but that text is the
multinomial combiner's own spec and multinomial is a forest this surface
refuses, so it is a VOCABULARY precedent only. The load-bearing evidence is
the measurement.

NOT part of it, each recorded rather than deferred: the residual-variance
prior (its anchor is the named `dbarts(..., sigma =)`, already original-scale
end to end and already pinned per sweep by `setSigma`); the tree prior (scale
free); the variance forest's leaf prior (an inverse-gamma-shaped
`ConstantVarianceLeaf`, `[[model.hpp:240@4c018187]]`, which sits outside `forests_` and is
not addressable - a door, not coverage).

## 2. Units, per response family

The getter reports `response.scale = response_->fitScale()` and
`response.shift = response_->fitShift()` through one virtual call, so the
implementation carries no family switch; this table is documentation.

| family | units of prior.scale / prior.mean | response.scale | response.shift |
|---|---|---|---|
| gaussian | response | `range_` | `range_*0.5 + min_` |
| aft | log survival time, anchored on the OBSERVED log times at creation (censored rows sit at their bound) and never re-anchored as censored latents are imputed, so late-chain imputed times routinely leave `[min_, max_]` while `prior.scale` stays the log-T-units forest-total sd | `range_(log T)` | `range_*0.5 + min_` |
| Student-t | response (delegates) | `range_` | as gaussian |
| grouped | response (delegates); the tau prior is a separate object | `range_` | as base |
| probit, weighted binary | probit latent | 1 | 0 |
| ordinal | probit latent, relative to the PINNED first cutpoint `gamma_1 = 0` | 1 | 0 |
| logistic | logistic latent (log-odds) | 1 | 0 |
| nbinom | log-odds `psi`; `E[y] = r exp(psi)` with `r` redrawn per sweep, so a mean-scale reading needs the current dispersion | 1 | 0 |
| multinomial | softmax log-odds, per category forest; scales come from the combiner map (`pi*sqrt(3)/sqrt(2)`, `[[combiner.hpp:285-287@4c018187]]`) | 1 | 0 |
| BCF | response; forest 0 `range_ * s`, forest 1 `range_ * sdModerate * s / 0.674`, `k` fixed at 1 | `range_` | `range_*0.5 + min_` |

The BCF row is the SCALED working response's sd (`scaledResponseSd`,
`[[chain.hpp:4252@4c018187]]`) multiplied out to response units; the getter multiplies by
`fitScale()` uniformly, so this is a table fact, not a code branch. Degenerate
cases, stated because R-side arithmetic gets them wrong: a constant
offset-adjusted response takes `range_ = 1.0`, and `n = 0` takes
`(min_, max_, range_) = (0, 0, 1)` (the `rescale()` guards, `model.hpp`).

## 3. What prior.sd MEANS, per leaf model

One sentence governs the whole tag vocabulary, and both the help and the
header lead with it:

> `prior.scale` and `prior.sd` describe the LEAF-PARAMETER scale of the forest
> total. They equal the prior sd of `f(x)` at every x for the constant leaf
> only. For the other three leaf models the prior of `f(x)` is x-dependent and
> `prior.sd` bounds it in a leaf-model-specific direction, stated below.

Four location-scale leaves ship, all carrying `double scale`, all reachable
from R: `ConstantGaussianLeaf` (`[[model.hpp:154@4c018187]]`),
`MonotoneConstantGaussianLeaf` (`[[model.hpp:506@4c018187]]`, `priorSd = (constrained ? cInflation :
1.0) * scale / k` with `cInflation = sqrt(pi/(pi-1))` installed in the `Chain`
constructor, `[[chain.hpp:568-569@4c018187]]`), `LinearGaussianLeaf` (`[[chain.hpp:914@4c018187]]`,
`drawFromPrior` draws all q+1 coefficients iid `N(0, (scale/k)^2)`) and
`GPGaussianLeaf` (`[[chain.hpp:1276@4c018187]]`, prior draw `f = (scale/k) L_C eps`). Public
surfaces: `dbarts(..., monotone =)`, `dbartsPriors$linear()` / `$gp()`.

The getter carries `leaf.model` in `{constant, monotone, linear, gp}` and the
help states, per tag:

- **constant**: EXACT. `prior.sd` is the prior sd of `f(x)` at every x and
  `prior.mean` is its prior mean.
- **linear (q leaf covariates)**: `prior.sd` is a LOWER bound, attained at the
  standardized covariate origin; in general `sd(f(x)) = prior.sd *
  sqrt(1 + ||z(x)||^2)` with `z` the internally standardized leaf covariates
  (a missing value maps to `z_j = 0`; a constant or all-missing column keeps
  sd 1 and contributes 0). `prior.mean` is exact. The help gives the formula.
  MEASURED, 400 draws at four rows: predicted 1.000/1.414/1.732/3.000 against
  measured 0.977/1.331/1.646/2.946; reconfirmed at 4000 draws
  with a constant column, training missingness and a predict-time NA
  (predicted 1.204/1.769/1.019/1.012/2.782/1.919/3.519 against measured
  1.200/1.777/0.998/0.999/2.756/1.906/3.518).
- **gp**: `prior.sd` is an UPPER bound over x. It is attained at rows that
  reproduce a leaf member (the drawn function values at member rows ARE the
  parameters) and on over-cap leaves, which draw as constant leaves. At any
  other row - every test row, every `$predict` row - the draw is the
  conditional mean `c(x*)' C^-1 f`, whose prior variance is
  `prior.sd^2 * c(x)' C^-1 c(x) <= prior.sd^2` and decays to 0 as `x` leaves
  the leaf's data cloud, at which point every prior draw equals `prior.mean`
  exactly. `prior.mean` is exact everywhere and is the whole prior under
  extrapolation. MEASURED, 2000 draws, gp trained on `x1` in [0,1],
  `sd/prior.sd`: 0.989 / 1.003 at training rows, 1.005 at an interior row,
  then 0.974 / 0.723 / 0.0355 / ~1e-39 / 0.000 at `x1` = 1.25 / 1.5 / 2 / 5 /
  20. The kernel diagonal is `1.0 + nugget` (`[[model.hpp:1692@4c018187]]`) but that
  governs MEMBER rows only, and `max.leaf.size = 256` (the `gp()` default,
  `[[R/model.R:967@4c018187]]`) means an n = 200 fit never trips the over-cap constant-leaf
  fallback, so this is the DEFAULT gp configuration, not a corner.
- **monotone**: `prior.sd` is a LOWER bound in the interior. The per-leaf
  c-inflation targets a single leaf's post-truncation variance while the prior
  draw is per-tree from the monotone cone (`drawFromPriorForTree`,
  `[[chain.hpp:1555@4c018187]]`), a different law, so the realized sd runs 3% to 20% above
  `prior.sd`, largest in the middle of the constrained axis - real, not Monte
  Carlo (MEASURED 1.067-1.090 at z = 4.2-5.7 over 2000 draws; +2.9% to +20.4%
  at z up to 15.2 in a second run). `prior.mean` is NOT the prior mean of
  `f(x)` under an active constraint: the constrained marginal is skew with an
  x-dependent mean that tracks the constraint direction and spans at least
  +-1 and up to +-3.3 prior sds across the axis (MEASURED -1.05 to +1.16 at
  2000 draws; +-0.8 at 400; +-3.3 in a second run). The help says so and
  points at the monotone design note.

The setter is TOTAL over all four: it writes the one `scale` field each leaf
model carries, and that write is exactly well defined in every case. The
asymmetry between a total write and a bounded read is stated, not hidden -
principle clause 3's "named, measured and testable" difference.

## 4. Creation surface

- `dbartsModel` gains one slot, `prior.scale` (numeric, `NA_real_` default,
  validity `> 0` or NA). When finite it OVERRIDES the family switch; the
  engine sets `nodeScale = prior.scale / response_->fitScale()` at the one
  site that already sets the leaf scale (`[[chain.hpp:558-559@4c018187]]`), where the
  response object and its `rescale()` already exist - the grouped decoration
  above it delegates its transform, every family's `rescale()` ran in its own
  constructor, and the degenerate guards make the divisor never zero.
  `node.scale` stays the internal-unit primitive and is what the bridge reads
  when `prior.scale` is NA.
- **The same conversion runs in `setModel`.** The bridge's `setModel` reads
  `model.nodeScale` (`[[R_interface_bartcore.cpp:4469@4c018187]]`) and `Chain::setModel`
  writes `forest.leaf.scale = model.nodeScale / sqrt(m)`
  (`[[chain.hpp:1287-1288@4c018187]]`), so without this the round trip
  `$setModel(sampler$model)` - a documented no-op today - silently REVERTS a
  named calibration (MEASURED 1.5 -> 12.0, 8x, on a range-24 gaussian).
  `Chain::setModel` therefore applies `model.priorScale /
  response_->fitScale()` whenever `priorScale` is finite, exactly as creation
  does, RE-DERIVED against the CURRENT transform. Consequence for the one live
  in-repo caller: `xbart` re-uses a sampler across cells that share `n.trees`
  and calls `bartcoreSetModel(sampler, cellModel(cell), ...)`
  (`[[R/xbart.R:619@4c018187]]`, `cellModel` at `[[R/xbart.R:572@4c018187]]`), so carrying `prior.scale` on that
  model makes every cell - created or re-modelled - run the same named
  calibration and the CV loss surface stops depending on cell ordering.
  Alternatives rejected: writing the resolved `node.scale` back into the model
  SEXP (it freezes a transform-dependent number into a slot that means
  response units), and refusing `setModel` on a `prior.scale` sampler (it
  breaks xbart outright).
- **Authority rule.** `model@prior.scale` records the NAMED INTENT, in
  response units. The ENGINE holds what is IN FORCE, and `$getCalibration` is
  the only authoritative reader. The intent is applied at creation and at
  every `setModel`, and is never rewritten by the engine. A channel that
  re-anchors the transform (`setResponse` / `setOffset` with
  `updateScale = TRUE`, `setData`) moves what is in force without touching the
  intent - MEASURED, a named 1.5 is 0.1875 in force after an `updateScale`
  re-anchor - and the getter shows the move. Re-applying the intent is
  `$setCalibration(prior.scale = ...)` or `$setModel(sampler$model)`. The slot
  is documented as an intent record, the live value has a total reader, and O6
  pins the difference.
- R sugar on the node priors: `normal(k = 2, sd =, scale =)`, and the same two
  arguments on `linear()` / `gp()`. `sd` is converted with the resolved k
  (`prior.scale = sd * k`) and is REFUSED when k is a hyperprior, naming
  `scale =`. Exactly one of `sd` / `scale` may be given.
- New formal on `bart2()`, `rbart_vi()` and `bart()`: `prior.scale` (response
  units), forwarded into the spec. `node.scale` is not it. This is the formal
  bartCause forwards for free - its whitelist is DYNAMIC, computed against
  `formals(eval(bartCall[[1L]]))` union `formals(dbartsControl)`
  (bartCause's `R/responseFit.R` lines 213-216, bartCause's `R/bartc.R` line 22), and `node.prior` is a
  formal of `dbarts()` only - and it is what `xbart`'s own family switch must
  honor.
- **Refused at creation, by name, at three NEW sites.** BCF and multinomial
  build their forests in SEPARATE `Chain` constructors (`[[chain.hpp:661@4c018187]]`,
  `[[chain.hpp:701@4c018187]]`), so the single-forest conversion site is not on their path and a
  named `prior.scale` would be dropped in silence. The existing gates test
  `model@node.scale != 0.5` (`[[R/spec.R:400@4c018187]]`) and `model.nodeScale != 0.5`
  (`[[R_interface_bartcore.cpp:2277@4c018187]]`), neither of which fires on a
  `prior.scale`-named model, and multinomial has no node-scale gate at all
  (its creation path documents that the host node scale is deliberately not
  read, `[[R_interface_bartcore.cpp:3014-3015@4c018187]]`). So this design ADDS: (i) `"a named 'prior.scale'"` to
  the offender list in `R/spec.R`; (ii) `std::isfinite(model.priorScale)` to
  the BCF offender chain in the bridge; (iii) a first node-scale-class refusal
  on the multinomial creation path. Three sites in two languages, all of which
  must agree; O5 tests all three. Nothing here is inherited.

## 5. Mutation surface, mid-chain

| channel | prior.scale (response units) | leaf.scale (internal) | other |
|---|---|---|---|
| `setResponse(updateScale = FALSE)` (default) | unchanged | unchanged | the composition path |
| `setResponse(updateScale = TRUE)` | MOVES with the re-anchor | unchanged | the getter makes the re-anchor visible for the first time |
| `setOffset(FALSE)` | unchanged | unchanged | |
| `setOffset(TRUE)` | MOVES | unchanged | stan4bart's decaying warmup schedule becomes inspectable |
| `setWeights` | unchanged | unchanged | |
| `setSigma` | unchanged | unchanged | original-scale value |
| `setData` | ALWAYS re-anchors, so it moves | unchanged | structural replacement |
| `setModel` | RE-DERIVED from the model's `prior.scale` against the CURRENT transform when finite; otherwise replaced via `node.scale`, internal units | replaced | AND re-pins sigma for gaussian/aft with no variance forest (`[[chain.hpp:1301-1310@4c018187]]`, MEASURED 3.5 -> 1); AND re-calibrates a variance forest's scale leaf |
| `$setCalibration` (new) | WRITTEN, response units, every chain | derived | touches nothing else - no sigma, no tree prior, no DART |
| `storeState` / `setState` | adopted from the state | adopted (`[[chain.hpp:3055-3060@4c018187]]`) | fixes the stale-`model@node.scale` reader, since the getter reads the engine |
| warm start (`installTrees`) | ADOPTED from the donor | adopted (`installForest`, `[[chain.hpp:2951@4c018187]]`) | documented semantic, not a defect: a donor's trees were drawn under the donor's scale. Recipe: re-issue `$setCalibration` after `installTrees` to keep your own |
| xbart, k grid (`k` numeric) | held across cells, including the `setModel` branch | - | the grid sweeps `prior.sd = prior.scale / k` about a fixed anchor |
| xbart, k hyperprior (`k = chi(...)`) | held across cells | - | k is DRAWN every sweep in every cell, so the loss is computed under a moving shrinkage; the help says so, and a `prior.sd`-flavoured argument meets the sampled-k refusal here as everywhere |

## 6. The R surface

    cal <- sampler$getCalibration(forest = 1L)
    # numeric matrix, one ROW per chain, columns:
    #   prior.scale prior.sd prior.mean k k.has.hyperprior
    #   response.scale response.shift
    # plus attr(, "leaf.model") in {constant, monotone, linear, gp}

    sampler$setCalibration(prior.scale = 3.0, forest = 1L, updateState = NA)
    sampler$setCalibration(prior.sd = 1.5, forest = 1L)   # sugar; fixed k only

`updateState = NA` matches `setResponse` / `setOffset` / `setWeights` /
`setSigma`; without it the setter inherits the recorded save/load defect that
`inst/tinytest/test-sampler-state-format.R` pins ("a `setModel(node.scale)`
issued after the last `storeState()` no longer survives a save/load
re-creation").

Refusals, each with an honest reason in the message:

1. `prior.mean =` is refused; the message names the lever,
   `setOffset(rep(-prior.mean, n))`, and the recipe.
2. `prior.sd =` under a k hyperprior is refused; the message names BOTH
   remedies - write `prior.scale`, or pin k through `setModel`. The binary
   default IS the sampled path (`[[R/model.R:392@4c018187]]`, `k <- if (monotone ||
   !binary) 2.0 else chi(1.5, 2.0)`, so probit/ordinal/logistic/nbinom default
   to a hyperprior and gaussian to fixed 2.0), and a named sd under it drifts
   silently (MEASURED 2.97 -> 1.98 in five sweeps with `leaf.scale` pinned).
   `prior.scale` is honored EXACTLY under a sampled k, because the hyperprior
   scales it rather than replacing it. The engine's own hyperprior defaults
   are `degreesOfFreedom = 1.5`, `scale = HUGE_VAL` (`[[model.hpp:2448-2449@4c018187]]`).
3. BCF and multinomial forests are refused; the message names the calibration
   map that owns them.
4. Host-shell samplers are refused through the existing `refuseHostMutation`.
5. A non-finite or non-positive `prior.scale` is an ERROR, not a refusal.
6. DART samplers are NOT refused. The `setModel` DART guard wraps only
   `splitProbabilities` (`[[chain.hpp:1331-1342@4c018187]]`) while the R method refuses
   outright (`[[R/dbarts.R:1004-1011@4c018187]]`), so this is the concrete gain over the
   existing route, and it is MEASURED: a live DART sampler refuses `$setModel`
   with the documented message while accepting a `leaf.scale` write through
   `setState` (0.1 -> 0.3) and running 40 clean sweeps afterwards with a full
   `varcount`. O5 keeps it as a shipped test.
7. The variance forest is not addressable (it is outside `forests_`), stated
   as a door.

Exactness rule the implementer MUST honor: `setCalibration` derives
`leaf.scale = prior.scale / (fitScale * sqrt(m))`, compares it BITWISE against
the current value, and SKIPS the write when equal, so a get-then-set cannot
perturb the last bit and move the equivalence baseline.

## 7. Oracles

| oracle | shape | tolerance | why it discriminates |
|---|---|---|---|
| O1a | two composed probit arms, ranges 16x apart, CENTERED so the shift matches, same named `prior.scale`, 120 Albert-Chib sweeps, **`n.chains = 2` and every chain compared** | `max|diff| <= 1e-12` | today the arms differ ~4x in posterior sd of f; MEASURED 1.44e-14 on a correct implementation, and it catches a no-op setter (5.68), an internal-units setter (7.62) and a chain-1-only setter (2.33 on chain 2) |
| O1b | same, arms NOT centered, each applying the documented `setOffset(rep(-getCalibration()[["prior.mean"]], n))` recipe, same chain rule | same | without the recipe it measures 2.59; it also catches a getter reporting `prior.mean = 0` (2.59); it is the only test exercising getter, setter and recipe as one contract |
| O1 SCOPE | stated in the help AND the test file: O1a/O1b compare arms at the SAME `numTrees`, so they see only `fitScale`-dependent errors and are BLIND to a setter error that is a fixed function of `prior.scale` (MEASURED: a `sqrt(m)`-forgetting setter passes both at 2.6e-13 and 2.0e-14 while reporting 10.6066 for a requested 1.5) | - | naming the blind spot is what sends the reader to O4b |
| O2 | named composition vs engine-native probit: posterior sd of f and RMSE | band, single DGP | smoke check (0.4571 vs 0.4561 measured) |
| O3 | prior draws, per LEAF MODEL, each with a row exercising its bound: constant `= prior.sd`; linear `= prior.sd*sqrt(1+||z||^2)` including a constant column, training missingness and a predict-time NA; **gp on member rows `= prior.sd` AND on an extrapolation row `<< prior.sd`, decaying to 0**; **monotone sd `> prior.sd` in the interior AND mean sweeping several `prior.sd` with the constraint's sign** | RELATIVE, stated with the draw count (se of an sd estimate ~ `sd/sqrt(2R)`; 6% band at 1500 draws), and the gp/monotone rows assert INEQUALITIES, not equalities | it is the only test of the section-3 contract, and an equality-only version would pass an engine doing the opposite |
| O4a | get-then-set is bitwise inert | bitwise | protects the equivalence baseline |
| O4b | set-then-get FIDELITY: `setForestPriorScale(P)` then `forestCalibration()` returns `P` bitwise, on every chain | bitwise | the only member with power against a setter error that is a fixed function of `P`: the `sqrt(m)`-forgetting setter returns 10.60660 for `P = 1.5` |
| O4c | static `m` falsifier: a correct setter reports `prior.scale = 1.5` at `m = 50` and at `m = 200`; the `sqrt(m)`-forgetting setter reports 10.60660 and 21.21320 (ratio 2.0000) | exact | cheap, and it catches the `m` dependence arms at equal `m` cannot see. Two O1a arms at DIFFERENT `m` are not bitwise comparable even under a correct implementation, so this is the right shape |
| O4d | per-chain divergence after a divergent `setState` is REFUSED, not flattened | exact | the write-every-chain rule |
| O5 | refusal matrix: BCF and multinomial at CREATION (all three new sites) and mid-chain, host shell, `prior.sd` under a hyperprior (including through `xbart`), `prior.mean`, non-finite/non-positive, out-of-range forest, DART ACCEPTED, monotone/linear/gp ACCEPTED with the tag | exact | the contract's boundary |
| O6 | every row of the section-5 table, one assertion each, including `setModel`'s `prior.scale` re-derivation, `setModel`'s sigma re-pin, the warm-start adoption, the authority rule (intent slot vs in-force engine value across an `updateScale` re-anchor), `updateState` save/load survival, and both `xbart` k arms | exact | the channels the surface must not surprise on |
| O7 | full tinytest, equivalence trio bitwise, `tests/cpp` | bitwise | RNG neutrality when `prior.scale` is unnamed; MEASURED, not asserted |
| O8 | flat-C parity in reshape S1: both entries, `structSize` canary INCLUDING a caller that omits the trailing fields, `DBARTS_HAS_FIELD`, refusal matrix, bit-parity with the R getter | bitwise | the C contract |

The `1e-12` tolerance is PINNED, not deferred: MEASURED 1.443e-14 at 120
sweeps, 2.887e-14 at 400, 5.129e-14 at 1000 and 6.239e-14 at 2000, growing
like `sqrt(sweeps)`, so it holds with 16x margin at 16x the shipped sweep
count. Re-measure at the shipped sweep count before pinning the test file's
constant; do not loosen it.

## 8. The dbarts.h footprint, for reshape S1

Creation half: NONE, truthfully - the model crosses as SEXP
(`dbarts_sampler_create`, `[[dbarts.h:175-177@4c018187]]`) and the conversion is
engine-side, so a flat-C consumer reaches it with no header change.

Mid-chain half: one output POD, one enum, two X-list entries appended at the
END of `DBARTS_C_API_LIST`.

    typedef enum {
      DBARTS_LEAF_CONSTANT = 0,
      DBARTS_LEAF_MONOTONE = 1,
      DBARTS_LEAF_LINEAR   = 2,
      DBARTS_LEAF_GP       = 3
    } dbarts_leaf_model;

    /// Caller-owned output buffers for dbarts_sampler_forestCalibration, the
    /// dbarts_results contract: set structSize; EVERY member is a pointer and
    /// is filled only when both present-by-size and non-null, each over
    /// numChains; a zero structSize errors. Fields append below the marked
    /// boundary and never reorder. All quantities are in RESPONSE units (the
    /// family's latent units where the response is not rescaled).
    typedef struct dbarts_forest_calibration_t {
      size_t  structSize;      ///< caller sets to sizeof(dbarts_forest_calibration)
      double* priorScale;      ///< numChains; forest-total prior sd at k = 1
      double* priorSd;         ///< numChains; priorScale / k at the current k
      double* priorMean;       ///< numChains; prior mean of the forest total
      double* k;               ///< numChains
      double* responseScale;   ///< numChains; internal-to-response multiplier
      double* responseShift;   ///< numChains; internal-to-response offset
      int*    kHasHyperprior;  ///< numChains; THIS FOREST's own k law (not the
                               ///< sampler-wide dbarts_sampler_kIsSampled,
                               ///< which reads the sampler option and
                               ///< disagrees on BCF and multinomial)
      int*    leafModel;       ///< numChains; dbarts_leaf_model, qualifying
                               ///< priorSd and priorMean (see below)
      /* 1.0-0 field boundary: appends go below, never above. */
    } dbarts_forest_calibration;

    #define DBARTS_FOREST_CALIBRATION_INIT { sizeof(dbarts_forest_calibration) }

    X(int, dbarts_sampler_forestCalibration, \
      (const dbarts_sampler* sampler, size_t forest, \
       dbarts_forest_calibration* out), \
      (sampler, forest, out)) \
    X(int, dbarts_sampler_setForestPriorScale, \
      (dbarts_sampler* sampler, size_t forest, double priorScale), \
      (sampler, forest, priorScale))

**Why the two trailing members are `int*` and not `int`** (compiled and
measured, `cc -std=c99`, arm64 LP64). With two trailing `int`s, tail padding
makes `sizeof` identical with and without the last field - 56 through
`responseShift`, 64 with one `int`, 64 with two, `offsetof(leafModel) = 60` -
so a caller omitting `leafModel` sets `structSize = 64` and
`DBARTS_HAS_FIELD(..., leafModel)` returns TRUE for a field it does not carry,
while the exact-`sizeof` static assert cannot see a future sub-word append at
all because it lands in existing padding. That would be the first POD to break
the header's stated invariant (`dbarts_results`' own Doxygen: the library
fills only fields whose end offset falls within `structSize`). With pointers
the arithmetic is honest: 64 with one, 72 with two, and the omitting caller's
`HAS` is FALSE. The uniform pointer shape also gives the two members the "null
member skips" spelling the POD's Doxygen claims for every member, and the
per-chain shape the rest of the struct already has.

Prototype-view Doxygen (the `#else` branch) states, in this order: what the
getter fills and that `priorScale` is the quantity the setter writes; that
`priorSd` is `priorScale / k` per chain and moves every sweep under
`kHasHyperprior` while `priorScale` does not; the LEAF-PARAMETER sentence of
section 3 with equality for `DBARTS_LEAF_CONSTANT` only; then per tag -
`DBARTS_LEAF_LINEAR` a LOWER bound attained at the standardized covariate
origin, larger by `sqrt(1 + ||z(x)||^2)`; `DBARTS_LEAF_GP` an UPPER bound
attained at rows reproducing a leaf member and on over-cap leaves, elsewhere
`priorSd^2 c(x)' C^-1 c(x)` decaying to 0 as x leaves the leaf's data cloud,
where every draw equals `priorMean`; `DBARTS_LEAF_MONOTONE` a LOWER bound in
the interior (realized sd a few per cent to ~20% above it) whose `priorMean`
is NOT the prior mean of `f(x)` under an active constraint, that marginal
being skew with an x-dependent mean spanning several `priorSd` along the
constrained axis; and that `priorMean` is exact for the constant, linear and
gp leaves. Returns 1, or 0 without touching `out` when `forest` names no
forest; errors on a zero `structSize`.

The setter's Doxygen states: it restates forest `forest`'s leaf prior on EVERY
chain so the forest total's prior sd at `k = 1` is `priorScale`, in response
units; `k`, the response transform, sigma and the tree prior are untouched; it
takes effect on the NEXT sweep and never reinterprets leaf values already
drawn; a write reproducing the current internal scale bitwise is skipped, so a
read-then-write is inert; to move the prior MEAN, shift the reported fit with
`dbarts_sampler_setOffset`; the leaf model qualifies the write exactly as it
qualifies the read. TWO ERROR CHANNELS, stated because the header's global
contract raises on invalid arguments: a CAPABILITY answer is a RETURN VALUE -
0, touching nothing, when `forest` names no forest or a combiner owns this
forest's calibration (a two-forest or multinomial sampler) - while a MALFORMED
VALUE RAISES, namely a non-finite or non-positive `priorScale`. 1 = accepted.
The flat surface deliberately carries no `prior.sd` spelling, so it has no
sampled-k refusal; that sugar and its refusal are R-side only.

Conventions honored, each checked: `1 = accepted, 0 = refused` and
out-of-range returns 0 for an `int` entry (`dbarts_sampler_forestFits`,
`[[C_interface.cpp:489@4c018187]]`; reshape binding decision 5 reserves `Rf_error` for the
`size_t` probes, which have no error channel); NO chain index, filled over
chains (the `dbarts_results` and `forestFits` shape, and reshape S1's
`forestAmplitudes`); `Forest`-infix naming (`setForestWeights`,
`setForestBasis`, `numForestAmplitudes`); `structSize`-first POD with a marked
boundary, read through reshape S1 item 1's `DBARTS_HAS_FIELD` generalization
rather than a bespoke `HAS` macro, with exact-`offsetof` asserts beside the
`dbarts_results` ones (`[[C_interface.cpp:68-77@4c018187]]`) and the exact-`sizeof` assert
whose job is to force an appending author to update them (`[[C_interface.cpp:78-80@4c018187]]`); a zero
`structSize` errors (`[[C_interface.cpp:135@4c018187]]`); additive append at the END; every
member pointer-shaped so the presence test cannot lie; the entries borrow
nothing and retain nothing, so no ownership sentence is owed; no field name
collides with a shipped entry. `DBARTS_C_API_MAJOR` / `MINOR` do not move
(reshape binding decision 8); the hash re-bake is the acknowledgment.

Forward compatibility: `forest` is a parameter, so the general basis family
(multiforest-extension-surface M4) relaxes the refusal in its own guard body
and moves no header. Record that in the reshape landing note.

Cost of NOT carrying this in reshape S1: a post-release `DBARTS_C_API_MINOR`
bump plus a lockstep bump of stan4bart's floor - the cost that plan's resolved
question 2 refused to pay for `setForestWeights`.

## S0. Signature freeze. No code. ~50 lines.

The signature must be settled before the last re-bake; the CODE need not be.

1. Write the section-8 footprint into `docs/plans/archive/dbarts-h-reshape.md` S1 as a
   numbered item (POD, enum, two X entries, the per-leaf-model Doxygen, the
   two error channels, the pointer rationale) and into TODO.
2. Record the forward-compat note (the basis family relaxes a guard body, not
   the header) and the `setModel` re-derivation rule the signature assumes.
3. PRECONDITION verified live and re-verified at slice start:
   `dbarts_sampler_numForests` is in the X-list (`[[dbarts.h:264@4c018187]]`), so reshape
   S1's own start condition is met. The reshape plan's zero-`structSize`
   anchor is already correct at this tip, so S0 inherits no errata.

The reshape item stays AMENDABLE until reshape S1 starts: if S2's
implementation falsifies a signature choice, the plan is corrected rather than
frozen. LIMIT, stated: S0 freezes a SIGNATURE, and the `setModel`
re-derivation was not a signature question, which is why that rule is written
into the design text above and into S0's record rather than left to the
implementer.

rng: NEUTRAL (no code). Gates: none beyond review - the two documents carry
the signature verbatim and agree.

## S1. Creation half. ~358 lines.

1. Engine `SamplerOptions.priorScale` plus the conversion at `[[chain.hpp:558-559@4c018187]]`
   (~18).
2. The SAME conversion in `Chain::setModel` plus its bridge read (~15).
3. Bridge slot read and validation (~20).
4. The three new BCF/multinomial refusal sites, two languages (~35).
5. `A_class.R` slot, prototype and validity (~20).
6. `normal` / `linear` / `gp` gain `sd =` / `scale =` plus the hyperprior
   refusal (~55).
7. `spec.R` resolution (~40).
8. `prior.scale` formal on `bart2` / `rbart_vi` / `bart` plus forwarding (~50).
9. `xbart`: its own family switch, the `setModel`-branch path, the
   hyperprior-k help sentence (what the loss means under a moving shrinkage)
   and the `prior.sd` refusal there (~35).
10. Man pages, the authority-rule paragraph, a `docs/design/prior-defaults.md`
    paragraph (~70).

Tests riding this slice: O1a/O1b (`n.chains >= 2` arms), O2, O3, O5's creation
refusals, and O6's `setModel` and `xbart` rows.

rng: NEUTRAL on every default path - byte-identical while `prior.scale` is
`NA_real_`; the named path is opt-in and no default moves. Gates:
`R CMD INSTALL --preclean .` into a private library (headers move);
`cd tests/cpp && make clean && make && ./test_bartcore` plain AND the
ASAN+UBSAN leg (new engine code in the `Chain` constructor and `setModel`);
full `tinytest::test_package("dbarts")` with FAILURES == 0 and NO snapshot
regenerated; the trio BITWISE against the MANIFEST hashes current at slice
start (re-read the MANIFEST rather than trusting this line;
`equivalence-a825263.rds` -> 35/35 "identical draws (same RNG stream)" under
`--strict-coverage`, `bcf-equivalence-a825263.rds` -> 11 scenarios,
`multinomial-equivalence-1027be5.rds` -> 10, no max-|z| line anywhere);
`air format --check .` and `lintr::lint_package()`; full local `R CMD check`
on a clean-copy tarball (R/ and Rd move); a `_pkgdown.yml` entry if any new Rd
topic appears - the new formals ride existing topics, so expect none.
ABORT: any trio deviation; a forced snapshot regeneration.

## S2. Mid-chain get/set. ~360 lines.

1. `Chain::forestCalibration` / `Chain::setForestPriorScale`, the `facade.hpp`
   virtuals and `sampler.hpp` (~75).
2. Bridge entries, refusals, registration, NAMESPACE (~120).
3. R5 `$getCalibration` / `$setCalibration` with `updateState` (~100).
4. `man/dbartsSampler-class.Rd` and a new `docs/design/nameable-calibration.md`
   (~65).

Tests riding this slice: O4a-d (including set-then-get fidelity and the static
`m` falsifier), the rest of O5, the rest of O6, and the `tests/cpp` component
test.

Carried from the S1 landing (the S1 landing note has the detail), inside this
slice's budget: (i) decide the `NaN` refusal site - tighten the R-side checks
(`is.na(NaN)` is TRUE, so `validateNamedScale` passes it today) or record the
bridge as the site; (ii) the NEWS bullet dedupes the sampled-k refusal content
S1's bullet already carries; (iii) fix the stale coverage comment at
`[[test-calibration-creation.R:100-104@4c018187]]`; (iv) pin the HETEROSCEDASTIC
creation-time calibration - it already works ungated (the conversion runs
before the variance branch) but ships no test, so the feature matrix carries
`?` for it; a one-arm prior-draw pin settles the cell.

rng: NEUTRAL - the surface is inert until called, and a write reproducing the
current internal scale is bitwise-SKIPPED, so a get-then-set cannot move the
baseline. Gates: as S1 (`--preclean` is MANDATORY here - `facade.hpp` virtuals
move and stale objects bus-error), with the trio expected bitwise for the same
reason; plus the `saveRDS`/`readRDS` state round trip through `updateState`.
ABORT: any trio deviation; a get-then-set that is not bitwise inert.

## S3. Flat C. An item INSIDE dbarts-h-reshape S1. ~295 lines.

Not a separate slice and not a separate re-bake.

1. POD, enum, the two X entries and the per-leaf-model Doxygen (~80).
2. `C_interface.cpp` bodies plus exact-`offsetof` and exact-`sizeof` asserts
   (~90).
3. O8 in `inst/tinytest/capi/consumer.c` and `inst/tinytest/test-capi.R`,
   including an omitting-caller `structSize` canary (~125).

Carried from the S2 landing (detail in the S2 landing note): (a) the engine
`Chain::forestCalibration` getter lacks a forest-index bounds check while the
setter has one; the flat getter's frozen signature owes the "return 0, touch
nothing" capability answer, so the check goes in the engine or is repeated in
`C_interface.cpp` - decide and test here. (b) Cosmetic, fix if touching the
site anyway: `refuseBCFMutation` fires before argument validation, so a BCF
`setCalibration(prior.mean = 0)` reports the BCF refusal rather than the
`prior.mean` one.

rng: NEUTRAL (header and `C_interface` only). Gates: reshape S1's list as it
actually reads - `tests/cpp` plain AND ASAN from clean (that plan's S0 states
`ASAN_OPTIONS=detect_container_overflow=0`); `test-capi.R` as the load-bearing
gate; full tinytest; the trio bitwise, LABELLED A FORMALITY; `air format
--check .`; `dbarts.h` ASCII-clean, C99-clean and CXX17-clean. ABORT: a
re-signed X-list leaving `DBARTS_C_API_HASH` unchanged from the literal
recorded at slice start; any movement in `DBARTS_C_API_MAJOR` or
`DBARTS_C_API_MINOR`; any sparse-vs-dense divergence at anything but bitwise;
a missing `dbarts_sampler_numForests` in the live X-list at slice start.

## Slice order

S0 records the signature into the reshape plan first. S1 and S2 land next and
are pure additions with no header. S3 then executes as a numbered item inside
reshape S1, which is where its gates live. That order satisfies both
constraints - the signature settled before the last re-bake, the code free to
follow.

## Consumer migration, verified

- stan4bart `bartcore`: 0 required edits, rebuild only (already owed by
  reshape S1). It never calls `setResponse` and drives by `setOffset` +
  `setSigma`; `node.prior` already forwards through `bart_args` ->
  `dbartsSpec`, so the creation surface reaches stan4bart users with zero
  stan4bart edits. OPTIONAL one-liner: `mvbart.R` sets `node.prior` itself and
  can opt in.
- treatSens `dbarts-1.0`: 0 required edits from this arc. It builds its model
  in R (treatSens's `R/cibart.R` lines 83-87, `node.scale = if (binary) 3.0 else 0.5`) and
  passes the model SEXP through to `dbarts_sampler_create`, so the engine-side
  conversion serves it; it already owes an independent 3rd-argument edit at
  its four flat `setResponse` sites. OPTIONAL one-liner at treatSens's `R/cibart.R` lines 83-87
  closes its own calibration hazard.
- bartCause `dbarts-1.0`: 0 edits, and it forwards the new `bart2` formal for
  free (its whitelist is dynamic).
- bairrtt: R-API only, 0 edits.
- In-repo `inst/tinytest/capi/consumer.c`: the two entries (S3).

## NEWS bullets (inst/NEWS.Rd; the first rides S1, the rest S2. ~15 lines,
additive to the design budget - flagged and added at plan commit)

- A model's leaf-prior calibration can be named at creation, in response
  units: `prior.scale` (the forest-total prior sd at `k = 1`) as a
  `dbartsModel` slot and a `bart`/`bart2`/`rbart_vi` formal. Unset, nothing
  changes.
- A sampler reports the calibration in force per chain (`$getCalibration`:
  `prior.scale`, `prior.sd`, `prior.mean`, `k`, `k.has.hyperprior`, the
  transform and the leaf-model tag) and can rewrite its scale half mid-chain
  (`$setCalibration`). For non-constant leaf models the reported `prior.sd`
  is a stated bound, not an equality.
- Under a sampled `k`, `prior.scale` is honored and `prior.sd` is refused
  with both remedies named.

## Decision record

Every option below was ADOPTED AT ORCHESTRATOR DISCRETION under VD's
2026-08-12 proceed-at-discretion grant.

- **V1, the named primitive.** `prior.scale`, the forest-total prior sd at
  `k = 1`, with `prior.sd` as k-conditioned sugar. Rejected: `prior.sd` as the
  primitive (ill-posed on the binary default path); documented `node.scale`
  plus a reader (~350 lines cheaper, but every composition does the division
  itself and gets it wrong silently, which is the wall).
- **V2, the sampled-k rule.** Accept `prior.scale`, refuse `prior.sd` naming
  both remedies. Rejected: pinning k on write (a hyperprior is a model change);
  honoring at the current k with a warning (the drift is real and a warning
  does not stop it). `xbart` is not a fixed-k island, so the refusal has an
  xbart surface to meet and xbart's help owes the sentence.
- **V3, leaf-model scope.** Total setter plus a tagged, honest getter.
  Rejected: constant leaf only (refuses a well-defined write, leaves
  monotone/linear/gp users with no lever); reporting the x-dependent quantity
  properly (a larger design - `prior.sd` becomes a function of x and the
  setter's inverse stops being a scalar).
- **V4, dbarts.h footprint shape.** Output POD plus a scalar setter, every
  member a per-chain pointer. Rejected: an input POD (its own INIT idiom errors
  twice); shipping no flat-C half (forfeits the window).
- **V5, creation mechanism.** Engine-side conversion, a `prior.scale` model
  slot, and the `bart2`/`rbart_vi`/`bart` formal - including the `setModel`
  re-derivation, the authority rule and the three refusal sites. Rejected:
  R-side spec resolution (a hand-built model is refused or silently unserved);
  a `calibration =` argument object (+~120 lines and a new public S4 class
  every consumer's argument forwarding must tolerate).
- **V6, is `prior.mean` writable?** Read-only, with the offset recipe (the
  offset lever reproduces the built-centered arm to 1.97e-14 with the
  transform verifiably pinned). Writable at creation only was rejected (+~80
  lines and it must decide what the getter reports after a `setData`
  re-anchor); mid-chain is REFUSED outright (leaves are stored internally).
- **V7, multi-forest writer scope in v1.** Refuse BCF and multinomial by name,
  at the three new sites. Rejected: relaxing now (+~90 engine plus a BCF
  re-record, and the creation-side refusals still have to exist).
- **V8, getter shape in R.** A chains-by-fields numeric matrix with a
  `leaf.model` attribute. Rejected: a `chain =` argument returning a named
  vector; a list.
- **V9, setter's per-chain semantics.** Write every chain, and refuse a
  round-trip write when chains have diverged. Rejected: a `chain =` argument
  (+~40 lines and the C entry then breaks the no-chain-index convention). A
  single-chain O1a PASSES a chain-1-only setter, which is why the oracles pin
  `n.chains >= 2`.
- **V10, budget and slice count.** All four slices, ~1593. Rejected: S0 + S1
  (forfeits read/write symmetry, the `setState` staleness fix, the mid-chain
  lever, and the header window); S0 + S2 + S3 (leaves `bart2`/`bart` users
  with no calibration argument).

## Carried caveats, to verify at implementation time

1. **The monotone mean magnitude is configuration-dependent** (+-1.2 prior sds
   in one 2000-draw run, +-3.3 in another). The design states a BOUND, not a
   point value. Measure the shipped test's own configuration before pinning
   O3's monotone row.
2. **Three measurements are CARRIED, not re-run**: the anchor honored across
   nine family and decoration paths at 2000 draws (all within ~2 se); DART
   accepting a `leaf.scale` write through `setState` and running 40 clean
   sweeps; the `1e-12` tolerance out to 2000 sweeps (6.239e-14). Each closes a
   caveat in the design's own favour, so re-run each at implementation time -
   they are cheap, and O3/O5 re-derive two of the three as shipped tests.

## Records at landing

docs/design/r-c-division.md's mvbart sentence and its sized-budget line landed
with the records commit at this tip and need no further edit. What is still
owed: record in this arc's landing note that the 4.6x posterior-sd sweep cited
in "The wall that remains is calibration" was reproduced twice independently,
at 4.1x and 4.7x on different DGPs. Also bump
`docs/design/nameable-calibration.md`'s `Status:` line at landing, per
docs/plans/README.md.

## Landing notes

S0 - LANDED 4c866286, 2026-08-12. The frozen `dbarts_sampler_get/
setCalibration` signatures and the `prior.scale` slot shape were
written into `docs/plans/archive/dbarts-h-reshape.md` S1 (items 7+8, beside
the active-rows footprint) and `docs/plans/archive/c-api-growth.md`; no code.

S1 - LANDED c2a7e89b, 2026-08-13. All ten items shipped; the
conversion is a private `Chain::resolvedNodeScale` helper shared by
the constructor and `setModel` sites, returning `nodeScale` VERBATIM
(no arithmetic) when `priorScale` is non-finite - the default path is
bit-exact by construction, and the trio confirmed it. Oracles
O1a/O1b/O2/O3/O5/O6 each shown RED against two deliberately broken
preclean builds (creation conversion removed; setModel conversion +
refusals removed) before green. Carried caveats re-measured at
implementation time: the monotone shipped-configuration magnitude is
pinned as a BOUND (span 1.98 prior sd, sd ratio 1.095-1.136x, max <
1.25), never a point value; the anchor holds across all nine
family/decoration paths (0.73438, grouped 0.74210, vs 0.75 - within
2.1 percent); DART accepts a named calibration; the 1e-12 tolerance
carries ~20x margin at 16x the shipped sweep count. Budget: code
+319/-39 vs ~358 (under); tests 720 raw / 574 non-comment vs the ~320
oracle share - the independent review CORRECTED the implementer's raw
2.25x figure to ~1.55x dense-equivalent (293 lines are air-mandated
one-argument/closer formatting) and ruled ACCEPT AS IS: the nine-path
sweep is plan-directed (carried caveat 2a ships as a test) and is the
only coverage of the `fitScale` divisor across families. Process
note: the implementer finished-and-reported instead of stopping at
the raw 1.5x threshold; future budgets are stated in dense-equivalent
terms. Deviation, reviewer-accepted: O1b computes the offset shift in
R from the plan's own section-2 formula (`getCalibration` is S2's).
Implementation fact: the bridge slot read needs `RC_NA | RC_YES` plus
an EXPLICIT non-finite check - `+Inf` passes the GT-0 constraint
because `assertDoubleConstraint` returns early on NA and `Inf <= 0`
is false (precedent [[RIB:1150-1152@4c018187]], the sigma estimate). Gates:
implementer and independent reviewer batteries both fully green
(preclean private libs, tests/cpp plain + ASAN with zero reports,
tinytest 4288/0, trio identical on all three baselines under
--strict-coverage, air + lint_package clean, R CMD check Status OK
from a clean-copy tarball, pkgdown no new topic); x86 box leg green
(tests/cpp plain+ASAN, tinytest 4276/0 incl. both new test files,
equivalence.R statistical OK; gcc emits benign
-Wmissing-field-initializers in sampler.hpp - noted, non-functional);
CI green on the c2a7e89b push. Carried to S2: (i) `NaN` escapes the
R-side `validateNamedScale`/validity checks because `is.na(NaN)` is
TRUE in R, and errors only at the bridge ("named prior scale is NaN")
- the refusal-5 outcome holds via an unintended path; S2 decides
tighten-R-side vs record-the-bridge-as-the-site; (ii) S2's NEWS
bullet must DEDUPE the sampled-k refusal content the S1 bullet
already carries, not restate it; (iii) fix the stale comment at
[[test-calibration-creation.R:100-104@4c018187]], which understates the shipped
coverage (test-calibration-prior-draws.R pins the absolute prior sd
at m = 20, so a sqrt(m)-forgetting conversion fails outright).

S2 - LANDED d809b944 + 7da36dc3 (record correction), 2026-08-13. The
four items plus carried items i-iv, all reviewer-verified: (i) NaN now
refused R-side via `is.nan` in `validateNamedScale` and the validity
method, bridge check retained as depth; (ii) NEWS sampled-k content
deduped into the mid-chain bullet; (iii) the stale coverage comment
rewritten; (iv) heteroscedastic pinned as the tenth anchor-sweep arm
(absolute 0.75 target, 9 percent band - the unnamed path reports ~1.5,
RED by 100 percent), settling the feature-matrix `?` cell. Get-then-set
inertness confirmed bitwise by the implementer three ways AND by the
reviewer's independent 40-cell mid-chain sweep (named/unnamed x five m
x four response scalings x three chains; a 1-ulp-nudged write moves the
draws, so the harness is non-vacuous). O4b DEVIATION, accepted as
STRENGTHENED: the plan asked set-then-get BITWISE; shipped is < 4 eps
plus bitwise-across-chains. The general argument stands ((P/f)*f != P
for ~9.7-9.8 percent of positive pairs, independently measured twice),
but the originally recorded fixture evidence was FALSE - the shipped
fixture is bitwise-exact at n.trees 20/50/200 and rounds ONLY at
m = 25, P = 0.25 (the record had it inverted; root cause: a number
transplanted from a scratch script with a different fitScale).
Corrected in 7da36dc3, which also adds the m = 25/P = 0.25 arm so the
tolerance is exercised by a genuinely rounding cell; the pin is a
robustness choice - exactness is a property of the (P, f) pair, not of
the implementation, so a bitwise assertion would pin an accident.
Budget: code+docs 389 dense-equivalent vs ~360 (1.08x; the docs item
was 2.6x because ~65 covered a man page plus a whole new design doc);
tests 430 dense-equivalent vs the ~210 S2 residual = 2.05x, reviewer
ruled ACCEPT AS IS - the ~210 predates carried items i-iv and O6's
twelve-row channel table, and the reviewer's own count found zero
padding (~13 optional lines named, not taken). Implementation facts:
`leafScale` already rode `ForestStateData`, so no state field grew;
`statesAgree` in inst/common/stateContinuation.R extended to walk `k`
and `leaf.scale` anyway (discrimination verified); the R5 `forest`
argument is 1-BASED via `resolveForestIndex` per this plan's own
spelling, while `$getForestFits`/`$getForestVariableCounts` are
0-based - the mixed convention is TICKETED (TODO
r5-forest-indexing) for the reshape re-bake. NAMESPACE needed no edit
(useDynLib .registration covers the new entries). Gates: implementer
and reviewer batteries both green from scratch (preclean private
libs, tests/cpp plain + ASAN zero reports, tinytest 4400/0, trio
identical x3, air + lint_package clean, R CMD check Status OK
clean-copy, pkgdown clean, saveRDS/readRDS updateState round trip
asserted); x86 box leg green at d809b944 (tinytest 4388/0 incl. both
calibration files, equivalence.R statistical OK); CI on the push per
the records commit. Carried to S3 (inside reshape S1): (a) the engine
`Chain::forestCalibration` getter has NO bounds check on the forest
index (the setter does; harmless today, bridge is the only caller) -
the flat getter owes the frozen signature's "return 0, touch nothing"
capability answer, so put the check in the engine or repeat it in
C_interface.cpp; (b) cosmetic: on a BCF sampler `refuseBCFMutation`
fires before argument validation, so setCalibration(prior.mean = 0)
reports the BCF refusal rather than the prior.mean one.

S3 - LANDED at dbarts-h-reshape S1, ab3aa2fa, 2026-08-13 (implemented as
3a977b6d, amended during independent review). All three items shipped
exactly to this section's spec: the `dbarts_forest_calibration` POD, the
`dbarts_leaf_model` enum, the two X-list entries and the per-leaf-model
Doxygen (dbarts.h); the `C_interface.cpp` bodies
(`dbarts_sampler_forestCalibration` [[src/C_interface.cpp:856-887@c85d6db6]], `dbarts_sampler_setForestPriorScale`
[[src/C_interface.cpp:889-899@c85d6db6]]) plus the exact-`offsetof`/exact-`sizeof` asserts; the `consumer.c`
+ `test-capi.R` coverage, including the omitting-caller `structSize` canary
(`[[test-capi.R:955-961@4c018187]]`: a caller whose `structSize` stops below `leafModel`
gets that member skipped and poisoned, everything else still fills) and the
zero-`structSize` error canary (`[[test-capi.R:962-965@4c018187]]`). Both carried items
from the S2 landing note above shipped alongside: (a) the engine bounds check
on `Chain::forestCalibration` - `if (f >= forests_.size()) return
ForestCalibration{};` (`[[chain.hpp:985@4c018187]]`) - so the reader now answers an
out-of-range forest with a default-constructed calibration exactly as the
writer already refused it, pinned at the engine level
(`[[tests/cpp/test_sampler.cpp:4722-4733@4c018187]]`, `testForestCalibration`) and through
the flat getter's "return 0, touch nothing" contract (`[[test-capi.R:950-954@4c018187]]`);
(b) the `refuseBCFMutation` reorder in R5 `$setCalibration` - argument
validation (the `prior.mean`-not-writable and exactly-one-of
`prior.scale`/`prior.sd` checks) now runs BEFORE the BCF refusal
(`[[dbarts.R:1469-1489@4c018187]]`, comment "so a malformed call is answered on its own
terms rather than by the refusal that would follow a well-formed one"), so
`setCalibration(prior.mean = 0)` on a BCF sampler reports the `prior.mean`
refusal rather than the BCF one. No further items are owed forward from S2.
Budget, review verdict, deviations and gates are recorded at the
dbarts-h-reshape S1 landing note, not repeated here - this arc contributed
~295 of that slice's ~1310 re-priced total, per the budget line above.

THE CALIBRATION ARC IS COMPLETE. `prior.scale` is nameable at creation
(S1), readable and writable mid-chain over every chain (S2), and reachable
from the flat C API (S3) - the R-user-facing capability was already `S` in
the feature matrix before S3, since the flat-C gap never gated a
`dbartsSampler` caller; S3 closes the flat-C gap itself. Nothing here is
owed forward.
