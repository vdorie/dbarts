# one-forest-basis: open the lone basis-carrying forest

Status: NOT TAKEN 2026-09-29 (dec-A109: the one-forest basis model stays refused; the refusal points to the p + 1 forest form and to linear() leaves). Kept as the record of the option and of the engine findings.

agent: opus for the engine, bridge and R slice and for both evidence instruments; the exact gate's oracle gets a
second, independent derivation pass (sonnet is not enough for either)
rng: neutral on every configuration that exists today (the trio replays bitwise); the one-forest model is a NEW
configuration with no prior draws, so it carries the posterior-changing class's evidence (an exact-posterior gate, a
design note) plus an SBC measurement and equivalence scenarios, all gathered before landing
window: 1.0-0 (dec-A109, the maintainer: "You can open it for 1.0-0. It seems like a small addition.")
budget: ~1,600 lines, two slices landed in one push. S1 code: engine ~20 (5 code), bridge ~50, R ~120, man ~60,
NEWS ~8, docs/design ~80, tinytest ~300, tests/cpp ~80. S2 evidence: one exact-gate harness ~550, sbc.R ~180,
bcf-equivalence.R ~40 plus one re-record, two workflow words. Machine time ~4-5 hours, ~1 of it on the x86 box.

## Goal

A model of one forest carrying a basis, `y = shift + (B_i . a) f(x_i) + e`, fits under gaussian, probit and logistic
from `dbarts(forests = list(forest(basis = ...)))`, `dbartsSpec(forests = )`, and a data object carrying one basis
(`dbartsData(bases = list(b))`, which is also how `bart()` reaches it). Every reader, guard and fit-object surface
treats it as the amplitude-coupled model it is, and it lands with an exact-posterior gate, an SBC verdict and
equivalence scenarios for all three families.

## Context

- Ruling: dec-A109 reversed by the maintainer; the refusal it recorded is
  [`resolveSamplerSpec`](../../R/spec.R)'s `numForests < 2L` block, the one site both creation routes reach. The
  bridge keeps its own floor in [`applyAmplitudeSpec`](../../src/R_interface_bartcore.cpp) ("at least two
  per-forest parameter vectors").
- Record of the deferral:
  [Fork 6. K = 1 through the K-forest path](archive/binary-kforest-prior-default.md#fork-6-k--1-through-the-k-forest-path). Model
  and calibration: [The model](../design/multiplier-combiner.md#the-model),
  [The calibration map, general in K](../design/multiplier-combiner.md#the-calibration-map-general-in-k),
  [What this family does not do](../design/multiplier-combiner.md#what-this-family-does-not-do). Enabling-value
  claim: [D4. The general per-forest multiplier (basis/amplitude) channel - CLOSED (2026-08-13 to 2026-08-14)](../design/model-space-survey.md#d4-the-general-per-forest-multiplier-basisamplitude-channel---closed-2026-08-13-to-2026-08-14).
- Evidence template: [bcf-latent-evidence](bcf-latent-evidence.md#bcf-latent-evidence), whose
  [Decision 1 - the SBC arms](bcf-latent-evidence.md#decision-1---the-sbc-arms) and
  [Decision 2 - the exact gate](bcf-latent-evidence.md#decision-2---the-exact-gate) set how a new amplitude
  configuration earns acceptance: a deterministic exact gate per family, a measured SBC verdict adjudicated by
  [The chain-length ladders (the A4e protocol)](sbc-family-tiers.md#the-chain-length-ladders-the-a4e-protocol),
  equivalence scenarios, and poisons proving each instrument discriminates.

Four facts found while planning correct the record this ruling was made on.

1. **It is not VCBART's shape.** VCBART (Deshpande et al., Bayesian Analysis 2026; eq. 2 of arXiv:2003.06416v8) is
   `y = beta_0(z) + sum_j beta_j(z) x_j + e`, one ensemble per coefficient INCLUDING the intercept `beta_0`. So
   VCBART with p covariates is K = p + 1 here - a basis-free first forest plus one `forest(basis = ~ x_j)` each,
   amplitudes held with `update.amplitude = FALSE` - which already ships. The one-forest model is the rank-one
   product `(B_i . a) f(x)`: with basis `cbind(1, x)` the intercept and slope functions are forced proportional,
   `beta_0 = a_0 f`, `beta_1 = a_1 f`. Its exact VCBART reduction is one term with the intercept pinned: basis `x`,
   amplitude held, `y = shift + x f(z) + e`.
2. **The engine does not take K = 1 as it stands.** [`AmplitudeForestCombiner`](../../src/bartcore/combiner.hpp)
   rounds its forest count up to two, so a one-forest chain carries a phantom second forest in the amplitude
   layout: an all-ones basis and one extra amplitude, never drawn (the draw and combine loops run over the chain's
   own forests), but counted by `totalAmplitudes`, so the run's glue channel,
   [`bartcore_getForestAmplitudes`](../../src/R_interface_bartcore.cpp) and
   [`serializeGlue`](../../src/bartcore/combiner.hpp) would each carry q + 1 entries. The chain already clamps the
   one read this was known to break ([`Chain`](../../src/bartcore/chain.hpp)'s `numVariableCountForests`, pinned by
   [`testBCFLegacyVarcount`](../../tests/cpp/test_shape.cpp)). Fork 6's "bridge plus R, not engine work" missed it.
3. **Six bridge guards and four R gates key on a forest count of two or more**, where the property they mean is
   "carries amplitudes": [`isMultiForest`](../../src/R_interface_bartcore.cpp) (setData, setModel),
   [`refuseMultiForestWarmStart`](../../src/R_interface_bartcore.cpp),
   [`responseConduitIsFixed`](../../src/R_interface_bartcore.cpp) and the scale arm of
   [`refuseMultiForestResponseMutation`](../../src/R_interface_bartcore.cpp),
   [`testFitsAreUndefined`](../../src/R_interface_bartcore.cpp),
   [`refuseMultiForestTestOffset`](../../src/R_interface_bartcore.cpp); R's
   [`refuseMultiForestWarmStart`](../../R/bartcore.R), and in the fit object
   [`packageBartResults`](../../R/bart.R)'s `hasForestReporting` and forest naming,
   [`refuseDroppedForestChannel`](../../R/generics.R), [`fitAllowsKHyperprior`](../../R/generics.R). At one forest
   each would read the sampler as single-forest: `setModel` would reprice the forest off the host model, a donor
   warm start would pair one draw's trees with another's amplitudes, and the sampler's own `predict` would return
   the bare forest without its multiplier. The flat C API calls the same predicates, so it follows the fix; no
   header changes.
4. **Two fronts mishandle a bases-carrying data object at every forest count today.**
   [`xbart`](../../R/xbart.R) ignores `data@bases` and cross-validates a single-forest model (a K = 2 data object
   gives a result identical to the same data without bases, checked on the installed build);
   [`rbart_vi`](../../R/rbart.R) builds the coupled sampler and fails on a reader shape. Separately, a data object
   carrying `bases = list(NULL, NULL)` fits two basis-free forests, which the `forests =` route refuses ("forest 2
   needs a 'basis'"): the data route never runs [`resolveForests`](../../R/model.R).

## The one-forest model

Everything below is what the existing code already computes once the floors go; nothing new is designed.

- Mean: `y_i = shift + (B_i . a) f(x_i) + e_i` under gaussian; the index `offset_i + (B_i . a) f(x_i)` under probit
  and logistic. The shift is the response transform's, fixed at creation: the midrange of `y` net of the offset
  under gaussian, 0 under the latent families. Nothing learns a baseline. A row whose basis row is zero has its mean
  at the shift (gaussian) or `p = 0.5` plus offset (latent).
- Intercept: none is implied beyond the shift. A factor basis spans the constant (its level indicators sum to one
  on every row), so `forest(basis = ~ factor(z))` gives each level its own scale on one shared function. A numeric
  basis does not; `cbind(1, x)` is how a caller asks for one. This is open question Q2.
- Amplitudes: the forest carries a basis, so [`forestParams`](../../R/model.R) gives it the fixed-variance channel,
  `a ~ N(0, 0.5 I_q)` (`amplitude.prior.variance`), never the half-Cauchy (that channel is only for a basis-free
  forest). Amplitudes enter at 1.0 ([`rebuildAmplitudeLayout`](../../src/bartcore/combiner.hpp));
  `update.amplitude = FALSE` holds them there. `(a, f) -> (c a, f / c)` leaves the likelihood unchanged, so only the
  product is identified, as at bcf's treatment forest.
- Leaf calibration: k fixed at 1; the node scale is `sqrt(2/K) s / (0.674 c)` at K = 1, factor 1.414214, so the
  forest total's prior sd is 2.0982 s / c, s the family anchor ([`latentScaleAnchor`](../../src/bartcore/chain.hpp))
  and c the basis's median non-zero row norm. The index prior sd at a row of norm c is 1.4837 s, the all-basis value
  the K-aware law holds at every K. Its law is a product of two normals, peaked at zero: a probit model's
  `P(p < 0.01 or p > 0.99)` at such a row is 0.104 (logistic 0.087; Monte Carlo, 4e6 draws), against 0.239 for the
  single-forest binary default and 0.247 for the two-forest bcf shape. This is the shipped law's own K = 1 value,
  not a new default.
- Readers: `$getLeafPrior()` reports one forest - `amplitude.prior.variance` 0.5, `amplitude.prior.scale` NaN,
  `leaf.scale.factor` 1.414214, `leaf.scale.divisor` 0.674, `basis.row.norm` c; `$getForestAmplitudes()` a q x
  n.chains matrix; `$getForestFits()` the forest's internal-scale total, n x n.chains; `$predictForests()` n.new x 1
  x draws; `run()` carries `forestFits` and a q-row `glue`, which the bridge already emits on the coupling rather
  than the count. A `bart()` fit carries `n.forests = 1`, `forestFits` with a length-one trailing forest margin
  named `forest1`, `glue`, `bases`, and the single-forest `varcount` shape (its margin follows the count).
  `extract(type = "forest")` and `predict` recombine as at K >= 2; a bare `bases =` value positions itself on the
  one carrier.
- Tree readers: a `forests =` declaration already carries the forest column whatever its count (dec-A96, the
  `bartcore.forestsDeclared` attribute read by `hasForestColumn` in [`dbartsSampler`](../../R/dbarts.R)). A
  data-route K = 1 sampler has no declaration and one forest, so today's rule gives it no column (Q3).

## Decision

VD signs off on these before S1 starts. Each changes what users get.

- **Q1. Open it on the corrected premise?** (a) Open as ruled. Costs the budget above - engine, bridge guards, fit
  object, two front refusals, three evidence instruments - rather than a small addition. (b) Keep the refusal,
  rewrite its message to point VCBART users at the K = p + 1 spelling that ships, and correct D4 and the TODO.
  ~40 lines, no compute. Recommendation: (a) if VD sees users for the rank-one product (a dose- or group-scaled
  effect with one shared shape); otherwise (b). What changes it: a named consumer.
- **Q2. Intercept.** (a) Document only: state the model, say factor bases span the constant and numeric ones do
  not, show `cbind(1, x)`. (b) Refuse a basis with any all-zero row at K = 1, which catches `basis = z` but not
  `basis = x`. (c) Refuse a basis whose column space lacks the constant (one QR per creation and per
  `$setForestBasis`), which forbids a legitimate no-baseline model. (d) An engine intercept: a scalar location drawn
  each sweep, a new conditional needing its own exact gate - not 1.0-0 sized. (e) Append a ones column
  automatically, which silently changes the amplitude layout the caller reads. Recommendation: (a), the base-R
  precedent being `lm(y ~ 0 + x)`, which lets a caller drop the intercept.
- **Q3. Forest column on a data-route K = 1 sampler.** (a) Key the column on the model kind - several forests OR
  carries amplitudes - so both routes to the same model print the same table. (b) Keep dec-A96's literal rule
  (declared with `forests =`, or several forests), so the data route prints no column. Recommendation: (a); it is
  what dec-A96's "model kinds that have several forests" is reaching for.
- **Q4. The SBC bar.** (a) Landing needs a measured verdict with no plateau flag, as the latent bcf arms landed;
  matrix admission is recorded separately. (b) Landing needs admission to the SBC matrix. Recommendation: (a). (b)
  risks the chain-length finding the bcf arms hit, which here would block a model whose exact gate passes.
- **Q5. `bart()` reach.** (a) Through a data object only; a `forest()` formula term keeps its meaning of an
  ADDITIONAL forest. (b) New grammar, e.g. `y ~ 0 + z:forest(x1 + x2)`. Recommendation: (a).

Not a decision, and done in S1 unless VD objects: the two front refusals and the `list(NULL, NULL)` refusal of
Context fact 4. Each replaces a silent wrong model or a crash with an error naming the cause; the bases slot is new
in 1.0-0, so none needs a NEWS entry of its own.

## Constraints

- Every configuration that stays accepted draws bitwise as before: the combiner change only moves a count below two,
  and the re-keyed guards answer as before wherever the count is two or more. Any trio divergence is a leak: abort.
- No flat C header change: [`DBARTS_C_API_HASH`](../../inst/include/dbarts/dbarts.h) must not move. The flat API
  cannot create a multi-forest sampler; it drives one built from R and inherits the re-keyed predicates.
- Families: gaussian, probit, logistic only. aft, ordinal and nbinom stay refused at the existing gate, with no
  wording change past "multi-forest".
- A lone forest must carry a basis. A one-forest basis-free model (`bases = list(NULL)`) is a single BART with a
  redundant Cauchy scale, and stays refused, now by its own message.
- Out of scope: `xbart` support for amplitude models (its grid axes include `k`, which the coupling pins); a
  `bart()` formula spelling of K = 1 (Q5); an intercept conditional (Q2 d); the per-draw amplitude channel in flat C.
- Evidence is gathered against the S1 worktree's own library before anything lands (docs/plans/README.md,
  Landing). S1 and S2 land in one push, S1 first; if the decision rule below fails, neither lands and the finding is
  reported.

## Steps

S1 - code, in one worktree off `origin/bartcore`, with a private library.

1. Engine: [`AmplitudeForestCombiner`](../../src/bartcore/combiner.hpp)'s constructor floors its count at one, not
   two. Update the comments that state the round-up (the constructor's, `numVariableCountForests`'s, and
   [`Chain`](../../src/bartcore/chain.hpp)'s clamp, which stays and becomes inert). tests/cpp: a one-forest
   amplitude spec reports `numAmplitudes == q`, serializes one width, restores bitwise, and runs a sweep; extend
   [`testBCFLegacyVarcount`](../../tests/cpp/test_shape.cpp)'s leg (ii). Mutation: restore the round-up and confirm
   the new checks fail. `--preclean` (facade virtuals are untouched, but the header is).
2. Bridge: one predicate, "amplitude-coupled or several forests" (`numForests >= 2 || numAmplitudes > 0` off
   `SamplerShape`), replaces the count in the six guards of Context fact 3, keeping every message (a count still
   names "N forests" where it does). [`applyAmplitudeSpec`](../../src/R_interface_bartcore.cpp) accepts one
   parameter vector; its message says "at least one". Comments that say "numForests >= 2" as the meaning of
   multi-forest are rewritten to the capability.
3. R creation: [`resolveSamplerSpec`](../../R/spec.R)'s `numForests < 2L` refusal becomes two refusals at the same
   site: a lone forest without a basis ("a single forest needs a 'basis' to be an amplitude model; drop 'bases' for
   a single-forest model"), and any forest past the first without one on the data route, mirroring
   [`resolveForests`](../../R/model.R)'s message. Update the comments in
   [`forestBasisDeclarations`](../../R/model.R) and [`forestParams`](../../R/model.R) that name the K = 1 refusal.
   `unsupported`'s message "fit a single-forest model" stands.
4. R fronts: [`xbart`](../../R/xbart.R) and [`rbart_vi`](../../R/rbart.R) refuse a data object carrying bases, by
   name, beside `refuseCountsCarryingData`. Check by test that `bart()`, `dbarts()`, `dbartsSpec()`, `bartBT()` and
   `pdbart()` each fit K = 1 or refuse it by name (the last two through the test-fit guard).
5. R fit object and readers: [`packageBartResults`](../../R/bart.R) names forests and stores `forestFits`, `glue`
   and `bases` on the coupling, not on `numForests > 1L`, and stores `bases` whenever the data carries them, so a
   `keepFits = FALSE` fit still knows it was coupled. One R predicate on the fit (`bases` present) replaces the count
   in [`refuseDroppedForestChannel`](../../R/generics.R) and [`fitAllowsKHyperprior`](../../R/generics.R); R's
   [`refuseMultiForestWarmStart`](../../R/bartcore.R) takes the capability too (via
   [`samplerCarriesAmplitudes`](../../R/bartcore.R)). The `varcount` margin, [`fitSynopsis`](../../R/generics.R) and
   `extract(type = "varcount")` stay on the count. Tree readers per Q3.
6. Sweep: `git grep -n -i "at least two forests\|two-forest\|second forest\|past the first\|single-forest 'forests'\|more than one forest\|several forests"`
   over R/, src/, man/, vignettes/, inst/, docs/design/; rewrite every sentence that states or implies a two-forest
   floor. Known hits: man/dbarts.Rd (`forests`, "Options a two-forest model does not read"), man/forest.Rd
   (`basis`, `update.amplitude`'s "refused on a single-forest forests", Details' "Both forests' leaf scales" and
   "when a second forest is declared"), man/dbartsData.Rd (`bases`), man/xbart.Rd (the new refusal), man/bart.Rd
   (fit components), man/dbartsSampler-class.Rd (guards stated as "multi-forest"), the tree-reader sentences in
   man/bartBT.Rd and vignettes/working_with_saved_trees.Rmd if Q3 (a). Also settle man/forest.Rd's claim that
   `update.amplitude = FALSE` holds amplitudes "at their prior center": the engine enters them at 1.0; read a held
   fit's `$getForestAmplitudes()` and correct whichever is wrong.
7. tinytest: a new `inst/tinytest/test-bcf-one-forest.R` (the test-bcf*.R family). Per family: creation on the
   three routes and their agreement (`$getLeafPrior()` equal across routes); every reader's shape and value from
   The one-forest model; the reconstruction identity `train == shift + (B a) f` under gaussian and index equality
   under the latent families to 1e-12; the held-amplitude reduction (basis `x`, `update.amplitude = FALSE`, glue
   exactly 1 every draw); `storeState`/`setState` and save/load continuing bitwise; `$setForestBasis` width change
   remapping amplitudes; `n.grow.sweeps` running. Every re-keyed guard, each asserted by message: setData, setModel,
   donor warm start (R and bridge), `$predict`, `$setTestPredictors`, test offset, `setResponse`/`setOffset` with
   and without `updateScale`. The `bart()` fit: descriptors, `extract(type = "forest")` with and without
   `contribution`, `predict` with a bare `bases`, the `keepFits = FALSE` and `type = "k"` messages. Refusals: aft,
   ordinal, nbinom; `list(NULL)`; `list(NULL, NULL)`; xbart and rbart_vi at K = 1 and 2. Warnings counted with
   `withCallingHandlers`. Rewrite the pins in inst/tinytest/test-bcf-creation.R that expect
   ["needs at least two forests"](../../inst/tinytest/test-bcf-creation.R).
8. Docs, present tense, landing with the code: docs/design/multiplier-combiner.md gains "The one-forest instance"
   (the model section above, minus numbers already stated in "The calibration map, general in K") and loses the
   floor from its surfaces; docs/design/model-space-survey.md's D4 VCBART bullet states the intercept ensemble and
   the K = p + 1 spelling, and names the one-forest model as the rank-one product; docs/design/feature-matrix.md's
   bcf row and Gaps table (the xbart row gains the refusal). inst/NEWS.Rd 1.0-0: the multi-forest item says any
   number of forests from one, one sentence on the lone basis forest and on its intercept. Remove the TODO item
   `binary-kforest-k1-reachability`.

S2 - evidence, in the same worktree against S1's library. Run sizes and bars are fixed here, before any run.

9. Exact gate, `benchmarks/R/amplitude-one-exact.R`, built on [`exactBCF`](../../benchmarks/R/bcf-exact.R)'s
   enumeration and bcf-latent-exact.R's adaptive Gauss-Hermite. Design: one ordinal predictor, two cells, `n.cuts =
   1`, one tree, so the tree space is two trees and the leaf block at most two-dimensional; `z` balanced within
   cells; `n = 400`. Arms per family: basis `1` (q = 1) and `cbind(1, z)` (q = 2), each with amplitudes free and
   held - 12 arms. Given the amplitudes, gaussian leaves integrate in closed form and sigma by 1-D quadrature on
   the sampler's own prior; latent leaves by 12-node adaptive GH. Amplitude axes on the `sqrt(0.5) tan(t)` open
   trapezoid (201 points at q = 1, 81 per axis at q = 2), with a refinement check under 1e-6. Matched quantities:
   `E[(B_g . a) f_c]` per (cell, group), `E[||a||]` on free arms, `E[sigma]` under gaussian, `E[F(eta_g)]` under the
   latent families. Runner-up tree weight guard at 0.02. Chains: 100,000 kept draws thinned by 10 over three seeds
   (quick: 25,000 by 5, one seed); [`batchMeanSE`](../../benchmarks/R/linear-exact.R) per seed, `zBound = 4`. Add
   the script to `.github/workflows/exact-gates.yaml`'s gate list in quick mode.
10. Independent derivation: a second agent derives, blind from the engine and the design notes, the tree prior
    mass, the amplitude prior, the leaf scales at K = 1 and the per-configuration marginal, checks them against
    `$getLeafPrior()` read off each creation route, and checks one marginal and one posterior mean by prior Monte
    Carlo (2e7 draws) to its own error.
11. Exact-gate poisons, each run once in quick mode and recorded: (i) oracle leaf scale factor 1 instead of
    sqrt(2) - every free arm fails (names the K-aware default at K = 1); (ii) held-amplitude fits scored against the
    free oracle - fails on `E[(B . a) f]`; (iii) latent fits scored against the other link's oracle - fails. S1's
    engine mutation already covers the phantom amplitude.
12. SBC, in [`runSbcBCF`](../../benchmarks/R/sbc.R)'s and [`sbcBurnLadder`](../../benchmarks/R/sbc.R)'s machinery:
    arms `one-gaussian`, `one-probit`, `one-logistic`, basis `cbind(1, z)`, `n = 200`, 50 trees, keyed in
    [`sbcBurnSweeps`](../../benchmarks/R/sbc.R) by arm name. theta0 installs the drawn amplitudes through the state
    before the prior trees, as [`sbcInstallBCFGlue`](../../benchmarks/R/sbc.R) does. Functionals: `sigma`
    (gaussian); `norm.a`, `a1.over.a0`; raw `a0`, `a1` reported and waived as sign-ambiguous; `index_j` at three rows
    with each `z` arm present; `p_j` at those rows (latent). Ladder: 40,000 sweeps over 24 prior-drawn datasets;
    thin and burn from it. Verdict: `R = 200`, `L = 150`. Controls on any flag: held amplitudes; `n = 40`; an A4e
    point at 3x chain length, `R = 80`. Poisons: amplitude prior variance 1.0 in the generator (must flag `norm.a` or
    an `index_j`), and the link swap on `one-probit` (must flag `p_j`).
13. Equivalence: three scenarios appended to benchmarks/R/bcf-equivalence.R with literal seeds outside the guarded
    settings - gaussian with `cbind(1, dose)` (a non-orthogonal basis), probit with `~ factor(z)`, logistic with
    `cbind(1, dose)` - and one re-record; update benchmarks/baselines/MANIFEST and the baseline name in
    `.github/workflows/cpp-tests.yaml` and `exact-gates.yaml`.
14. Record: results under this plan's `## Results`, the exact-gate section in multiplier-combiner.md's new
    section, the feature-matrix evidence cell; if admitted per Q4 (b) or by its own criterion, the sbc.yaml matrix
    entries, else its exclusion note names the arms and why.

### Decision rule (fixed before any run)

Landing needs all of:

- L1. Exact gate, full mode: every matched quantity at `|z| <= 4` in all 12 arms (about 100 tests, family-wise
  false failure near 6e-3, accepted in advance); refinement under 1e-6; weight guard met. Quick mode passes too.
- L2. The equivalence trio replays every existing scenario bitwise (`--bitwise`, full scenario count, no
  "max |z|" line) on the reference build; the three new scenarios replay bitwise after the re-record.
- L3. Each exact-gate poison and each SBC poison reddens its named quantity. A poison that does not land means the
  instrument is not evidence: stop and redesign it, do not land.
- L4. SBC, per Q4: under (a), no functional's flag survives adjudication as a plateau - that is, every flag either
  clears under held amplitudes or shrinks in `ecdfDiff / band` at the A4e point. A plateau blocks landing and goes
  to an independent derivation of the flagged conditional.

Stop rule: a verdict settled at `R = 200` is not re-run at larger `R`; controls and A4e points run only on a flag.

### Where it runs, and compute

The local arm64 machine for everything except the x86 leg that an engine landing owes (benchmarks/README.md, "The
x86 leg"). Per-sweep costs measured today on the installed build, two-forest samplers as the upper bound: one tree
per forest at `n = 400`, 10.8 / 32.0 / 73.7 us (gaussian / probit / logistic); 50 + 50 trees at `n = 200`,
78.8 / 90.2 / 111.4 us. A one-forest sweep costs less; the numbers below use these.

| run | estimate |
|---|---|
| tests/cpp plain + ASAN; R-loaded ASAN on the new tinytest file | 20 min |
| full tinytest, twice (implementer, independent battery) | 30-40 min |
| equivalence trio `--bitwise` on the reference build, plus the re-record | 30-45 min |
| exact gate full (12 arms x 3 seeds x 1e6 sweeps: ~2 / 5 / 12 min by family) plus quadrature | 25-35 min |
| exact gate quick, three poisons, prior Monte Carlo | 15 min |
| SBC ladders (1.1e6 sweeps per family) | 5 min |
| SBC verdicts, three families, at a bcf-like 19,500 sweeps per replication | 15-25 min |
| SBC poisons, two runs at `R = 100` | 10 min |
| contingent: controls and A4e points (up to 3x chain length) | 0-90 min |
| x86 leg: tests/cpp plain + ASAN, full tinytest, trio statistical | 60 min |

About 3.5 hours without contingencies, 5 with them. The SBC verdicts may run as three parallel processes (they make
no timing claim); the ladders, which record per-sweep cost, run one at a time on a quiet machine.

## Verification

    cd tests/cpp && make && ./test_bartcore
    R CMD INSTALL --preclean -l <lib> .
    R_LIBS=<lib> Rscript -e 'tinytest::test_package("dbarts")'
    R_LIBS=<lib> Rscript -e 'tinytest::run_test_file("inst/tinytest/test-bcf-one-forest.R")'
    R_LIBS=<lib> Rscript benchmarks/R/equivalence.R compare <MANIFEST current> --bitwise      # and the two siblings
    R_LIBS=<lib> Rscript benchmarks/R/amplitude-one-exact.R quick                            # and without quick
    R_LIBS=<lib> Rscript benchmarks/R/sbc.R burn-one-probit 40000 24                          # and the other two
    R_LIBS=<lib> Rscript benchmarks/R/sbc.R one-probit 200 150 <thin> <burn>                  # and the other two
    Rscript -e 'lintr::lint_package()' && air format --check . && Rscript tools/check-rc-codoc.R . \
      && Rscript tools/check-win-drift.R . && Rscript tools/check-doc-freshness.R .
    tools/check-api-hash.sh                                                                  # hash unmoved

Plus the exact gates of `.github/workflows/exact-gates.yaml` in quick mode (a fit object's contents change), R CMD
check --as-cran from a clean staged tarball, the `inst/NEWS.Rd` parse gate, and `git grep` for the removed TODO
item. Expected: the decision rule's four clauses; tinytest failures 0; no API hash move.
