# bartcore <- main: a guided review tour

A reading guide for the `main...bartcore` diff, anchored at branch tip
**ae5b91d8**. Every stop names its files, its anchors, its paired design docs,
the evidence behind its claims, and the command that pulls up its diff.

    merge base: cb290550 (main tip; bartcore is 1273 commits ahead)
    scale:      680 files, +229,948 / -32,374
    docs:       +80,194 / -0  (205 files, ~35% of the diff, zero deletions)

Regenerate any number with `git diff --numstat main...bartcore -- <path> |
sort -rn`. The file count fell since the previous refresh (733 -> 680) because
a baseline prune deleted branch-added recordings.

**Citation convention**: prefer a symbol name to a bare line number
(`docs/plans/README.md`'s rule). Shipped code carries no `docs/` path at
all: a citation strip removed 316 such comments from 81 files, on the ruling
that a tarball install ships no `docs/` tree. And line numbers rot: every code
anchor in 29 of the 50 `docs/design` files was re-derived by content in six
passes, but `docs/plans` was not swept (3 of its 153 files moved), and
`tools/check-doc-freshness.R` checks anchor *identity* for
`docs/design/feature-matrix.md` alone. Every line number below was
grep-verified against ae5b91d8 at write time.

**One stamp to distrust.** `docs/design/feature-matrix.md:37` says its anchors
are verified by content "against the tree at 54dec2ab". That sha is not an
ancestor of this branch - it is a pre-rebase twin of `9f538bbf` (identical
patch-id; `git diff 54dec2ab 9f538bbf` is empty), so the anchors are sound but
`git show 54dec2ab` fails in a fresh clone.

## How to use this

Twelve stops, ordered by blast radius rather than diff size. Stops 1-9 are the
change; 10 is the evidence it works; 11 is sediment. Two are reshaped this
refresh: Stop 2 walks all thirteen response-model rows rather than listing
anchors, and Stop 9 walks the capability axes that cut across them.

    [ ] 0  Orientation                    docs, no diff
    [ ] 1  The engine swap                -8,718 / +21,085
    [ ] 2  model.hpp and the family rows  4,923
    [ ] 3  Moves and trees                869 + 2,158
    [ ] 4  The data layer                 1,857
    [ ] 5  Multi-forest                   2,085
    [ ] 6  The bridge                     -4,554 / +8,232
    [ ] 7  The shipped C ABI              -1,753 / +3,740
    [ ] 8  The R surface                  +22,656 / -3,227
    [ ] 9  The capability axes            cuts across 2, 5, 6, 8
    [ ] 10 The gate apparatus             ~88,800
    [ ] 11 Build sediment                 mostly deletion

Each stop closes with "Known open": doors and residue, not oversights. The
appendices index the arcs that landed most recently, with an entry anchor and
a gate for each, and the documents where decisions live.

---

## Stop 0 - Orientation

No diff reading. ~40 minutes.

1. `docs/architecture.md` (398 lines, current-state): Layering, Dispatch
   tiers, The facade, Model concepts, Tree moves, ColumnStore, The mutation
   surface and its transaction semantics, Tree storage forms, RNG
   architecture, Threading model, Gates.
2. `docs/design/feature-matrix.md` (1,032 lines). **The most valuable second
   read**: one row per response model, 26 capability columns across five
   tables, every `S`/`R` cell carrying an anchor, 49 footnotes, and a "Gaps"
   section (:924) collecting every `M` and `?` cell as a candidate work item.
   Stops 2 and 9 narrate its rows and columns; it is the authority where they
   disagree.
3. `docs/design/INDEX.md` (29 of 50 `LANDED`), `docs/plans/INDEX.md` (153
   plans), and `docs/README.md` on what each navigation surface is NOT for.
4. `inst/NEWS.Rd` (2,551 lines): 1.0-0 runs UPGRADING :5-77 (11 items), NEW
   FEATURES :78-955, C API :956-1023, BUG FIXES :1024-1995. Read a family's
   NEW FEATURES paragraph before its design doc - it states the claim in the
   terms a user would defend.
5. root `TODO` (329 lines, open work only): `second-review-followups` and
   `rc-gate` are the live status pointers, and the `release` section is the
   one ordered procedure.

**The one structural idea**: the leaf model is a compile-time template
parameter `L` (`accumulate`/`logIntegratedLikelihoodForNode` must inline),
forcing the `facade.hpp` type-erasure layer; the response family is a runtime
virtual chosen once per chain. Almost every "why is this file shaped like
this" question resolves to that sentence.

**RNG-effect gate classes**: neutral (tests/cpp + tinytest), shifting (+
regenerated snapshots, re-recorded equivalence baseline, statistical-z),
posterior-changing (+ exact-posterior gates + design note).

**Three renames that make greps lie**: `BCFForestCombiner` is now
`AmplitudeForestCombiner` (combiner.hpp:741), one of 20 identifiers a "BCF
sheds its spelling" cleanup renamed (state key `"bcf"` -> `"glue"`) - the old
name finds nothing, which does not mean BCF was removed;
`bartcore_createMultinomial(Counts)` are retired, multinomial going through
the same `bartcore_create` as every other family; and BCF's own R verb,
`bcf()`/`bartBCF`, lives in **bartCause**'s `dbarts-1.0` branch while dbarts
keeps the general engine it was built on (Stop 5).

**Two doc-vs-code divergences still open**: `ColumnStore` (data.hpp:501) uses
uniform 16-bit codes, no per-column width (`hot-layer-u8.md` is a measured
NO-GO), and it is still not arena-allocated - `docs/architecture.md` says both
in its own voice.

**One document that contradicts itself where a reviewer will land.**
`docs/plans/variance-forest-mutation-routing.md:415-426` opens its last
recorded door with "CLOSED: built at c95a5e83" and then, without deleting it,
continues with the original present-tense description of the hole
("`rebuildVarianceForest` has no `columnMaskSubtreeIsValid` backstop ... a
column-restricted variance forest restored through `setState` is not"). CLOSED
is the true half: `Chain::columnMaskStateFeasible` (chain.hpp:3418) is live
and both install arms gate on it (sampler.hpp:946, :1134).
`feature-matrix.md:1016-1022` records it correctly.

---

## Stop 1 - The engine swap

    git diff main...bartcore -- src/dbarts src/bartcore

`src/dbarts/` stays deleted (-8,718). `src/bartcore/` is 21,085 lines across
11 headers - growth from multi-forest generalization, nameable calibration,
threaded predict, and a response-family exhaustiveness sweep.

**Read order**: `bartcore.hpp` (24, umbrella); `facade.hpp` (906) -
`SamplerBase` :137-412 with **59 pure virtuals** (count it: `sed -n '137,412p'
src/bartcore/facade.hpp | grep -c '= 0;'`; matches
`tests/cpp/test_facade.cpp:4`'s own note - a higher count has picked up a
helper declared outside the class body), `SamplerFacade` :413, factories
`createConstantLeafSampler` :730, `createSampler` :787,
`createSamplerOverStore` :840, `createAmplitudeSampler` :868,
`createMultinomialSampler` :894; `sampler.hpp` (1,999) - `Sampler` :132,
`PredictorUpdateResult` :85, `PredictorUpdateSession` :114, `run` :298,
`predictColumns` :705 fanning out over `std::thread` workers through
`fanOutPredictSlabs` :602 exactly as `run` does, a test-only
`predictPartition` channel :51, the saved-tree ring readers `filledSavedDraws`
:515 / `savedSlotForDraw` :520, `installForests` :983; `chain.hpp` (5,512, the
biggest file) - `SamplerOptions` :40, `Results` :312, `VarianceForest` :404,
`Chain` :548, `setActiveRows` :1666, `interweaveGlueRidge` :1331,
`weightsDigest` :1611, `columnMaskStateFeasible` :3418, `buildVarianceForest`
:4176.

**What to look for**:

- RNG architecture and threading: read `docs/architecture.md`'s own sections
  (per-chain generator, worker `SIGINT` masking, `QueuedProgressSink`,
  callback-inline-only restriction) rather than any paraphrase.
- **The embedding contract**: `docs/design/bart-as-a-component.md`, the
  internal counterpart to `vignettes/dbarts-as-a-component.Rmd`. Its
  "Measured" section proves the per-sweep loop bitwise identical to a batched
  run and prices the tax of driving it that way (1.09-1.13x per call) and of
  composing K single-forest samplers in R rather than one batched multi-forest
  one (1.9x at K=2 up to 5.2x at K=8). Its mutation-legality table
  (`refuseMultiForestMutation`, `refuseMultiForestResponseMutation`,
  `refuseUndefinedTestFits`, R5's `refuseAmplitudeMutation`) is required
  reading before touching multi-forest mutation code.
- **The response-family switch is exhaustive.** Every `default:` arm on a
  `ResponseFamily` switch is deleted, so a 7th enumerator is a hard compile
  error (`-Werror=switch` in `tests/cpp/Makefile`, plus a cpp-tests CI step
  writing `CXX20FLAGS`, `CXXFLAGS` being inert under `CXX_STD = CXX20`) rather
  than a silent gaussian fallback. One missing arm was a confirmed-live
  defect: a negative nbinom count could underflow into a ~1.8e19 allocation.
- `recoverTreeParameters`/`applyNewData` used to hard-code `forests_[0]` where
  siblings looped every forest, so a whole-data replacement on a BCF or
  multinomial chain left forests 1..K-1 on the old grid; now asserted
  single-forest as a precondition.
- State restore is semantic, not bitwise (`state-continuation.md`).
  `stateFormatVersion` is 3 (`src/R_interface_bartcore.cpp:6380`), the
  readable floor also 3 since `recordedDraws` became required.

**Known open**: a per-forest weight is not part of saved state - a pipeline
installing a donor's state into a fresh engine starts with none, and the two
stored states compare `identical()` while fits diverge (0.057 fit difference
against a 0.035 range; decided in `bcf.md`, pinned only for the same-holder
round trip). An active-row mask is mirrored nowhere and lost on re-creation.
`within-chain-threading.md` stays CLOSED/NO-GO; `data-layout.md` stays
CLOSED-SHELVED, orphaned by block-fusion's own WONT-DO.

---

## Stop 2 - model.hpp and the thirteen family rows

    git show bartcore:src/bartcore/model.hpp

4,923 lines. **Leaf models**: `LeafModelCore` :40, `ConstantGaussianLeaf`
:155, `ConstantVarianceLeaf` :241, `MonotoneConstantGaussianLeaf` :507,
`LinearGaussianLeaf` :974, `GPGaussianLeaf` :1336. **Priors**: `CGMTreePrior`
:2117, `DartPrior` :2433, `ChiKHyperprior` :2507,
`ChiSquaredScalePrior::drawSigmaSqFromPosterior` :2540, `ResidualDfPrior`
:3994. **Responses**: `ResponseFamily` enum :2580 (six tokens),
`ResponseModel` :2618, `GaussianResponse` :2782, `ProbitResponse` :3057,
`OrdinalResponse` :3196, `LogisticResponse` :3517, `MultinomialResponse`
:3716, `AFTResponse` :3794, `TResponse` :4041, `NBResponse` :4338,
`GroupedResponse` :4706.

The enum's six tokens do not enumerate the shipped models. Thirteen rows reach
a user, and four (student, hazard, hurdle, the K-forest amplitude family)
arrive by an attribute, an ingestion rewrite, an R-level composition or a
combiner rather than by an enum value - which is why the exhaustive-switch
sweep in Stop 1 mattered, and why `feature-matrix.md`'s row list, not the
enum, is the inventory. `docs/plans/correctness-audit.md` (7 blocks, all
CONFIRMED) is the answer key for the math. Each row below names what to judge;
the matrix cell and its footnote carry the rest.

**gaussian** (:2782). The reference row: every construction surface, every
mutation channel, every composition column `S`; 21 equivalence scenarios, SBC
PASS 7/7, ~20 tinytest files, and `feature-matrix.md:931` records no gap at
all. A sibling row that behaves differently is the thing to question.

**student** (`TResponse` :4041). Not a family token: `resid.dist =
student(df)` (`R/model.R:1612`) resolves to a finite `resid.df` attribute at
`R/spec.R:397`, refused for a non-gaussian family at :390; the engine family
stays gaussian, so the gaussian row applies in tables 2-4. **Judge**: the
loglik now scores the *marginal* t density (`R/generics.R:118`), reading a
stored `resid.df` and erroring rather than guessing when it is absent - a
recorded live defect, fixed, but still a `?` cell because no test pins the
formula. `getLatents` returns scale-mixture *precisions*, not locations; `type
= "ppd"` is refused outright (:2798) since the noise would be gaussian;
composition with a variance forest is refused by name (`R/spec.R:521`, "not
yet shown to compose", reopenable). No `rbart_vi()`/`xbart()` reach, one
tinytest file, exact gate `t-exact.R`, and, alone among the heavy families, a
full CI SBC arm, PASS 4/4.

**probit** (:3057). The binary workhorse and the substrate under hazard; every
construction surface `S`, and `feature-matrix.md:941` records no gap.
**Judge**: the binary `k` default moved to `chi(1.5, 2)`, with `xbart()`
carved out at a fixed `k = 2` so a sampled `k` cannot collapse its own grid
axis (`R/xbart.R:222-224`, :328-331); `setWeights`/`setSigma` refused by name
(`src/R_interface_bartcore.cpp:2773`, :2898); `updateScale` *ignored* rather
than refused, the transform being the identity - a `-` cell, not `R`. And
there is **no probit arm in the CI SBC matrix** (`sbc.yaml` runs
discrete-selfcheck, gaussian, ordinal, nbinom, t, multinom), so in CI binary
calibration rests entirely on `logistic-reference.R` against
`BART::pbart`/`lbart`. Shipped weighted binary is logistic only; its PPD draw
for a weight-w row is `Binomial(w, p)`, not `w * Bernoulli(p)`.

**logistic** (:3517). Weights are observation counts, so `setWeights` is a
*model* change: MOD:3600 redraws the Polya-Gamma omegas against the new
counts. Saved state carries a byte-hash `weights.digest`
(`Chain::weightsDigest` chain.hpp:1611, virtual facade.hpp:159) and `setState`
re-derives the latents on a mismatch, the identity on a match; the weights
themselves still do not ride the state (`test-state-weight-pairing.R`). Open:
no `rbart_vi()` token, and no flat-C test coverage although the family is
reachable through `dbarts.h`.

**ordinal** (:3196). Cumulative probit with sampled cutpoints; `bart2` and
`dbarts()` only. **Judge**: the grow-from-root scan dropped missing rows from
the split likelihood while the no-split term kept them, biasing
`n.grow.sweeps` toward not splitting. `scanOrdinalCuts` (scan.hpp:105) now
doubles the entry layout when a node holds missing members and scores each
missing direction over the rows it seats; `growTreeFromRoot` (grow.hpp:179)
sizes `cutLogLikelihood` for the doubled layout unconditionally. Verified
against an independent R re-derivation to 1.53e-16; oracle
`tests/cpp/test_scan.cpp:337` `testMissingRouted`. Doors: no grouped-ordinal
surface, no warm start or `n.grow.sweeps` formals, an undocumented selecting
header attribute, a recorded `log_diff_exp` tail-precision residue.

**nbinom** (:4338). `r` restricted to positive integers; real `r` is a
recorded door. `setActiveRows` rebuilds the count-histogram kernel at every
mask change - the channel's one per-install cost. **Judge**: the calibration
lane found `r`/`agg.psi` a genuine, forest-size-invariant mixing ridge (two
chains hold disjoint `r` basins for 100,000 sweeps at identical `avg.mu`),
adjudicated H-MIX on two ladder points; the third is owed. The missing
rootogram panel is a rejected design alternative, and `plot` draws no burn-in
dispersion segment because `bart2`'s negbin runs one `run(n.burn, n.samples)`.

**multinomial** (:3716 plus `MultinomialForestCombiner` combiner.hpp:1574). K
forests behind a softmax; `bart2` builds the `dbartsSampler` directly, no host
shell, so `$fit` is the engine that ran. The mutation table reads `-`
throughout: eight methods refused by name, the public response channels being
`setCounts`, `setCategoryOffset`, `setCategoryTestOffset`. **Judge**: DART is
refused by name (`R/bart.R:901`) where it used to be dropped silently
(`forest.useDart = false`, chain.hpp:5219), alongside `split.probs`,
`monotone` and `variance`; per-forest row masking is a permanent model-level
refusal, the softmax margin being a log-sum-exp over the other K-1 forests. No
`getLatents`, no flat-C creation path, calibration refused. An ensemble-scale
SBC re-run moved the raw `f_ik` cells from undischargeable to PASS.

**aft** (:3794). Log-normal AFT, reducing exactly to Gaussian on uncensored
rows - a free correctness gate, and riAFTBART's model. `setSigma` open,
`setWeights` and `setData` refused. **Judge**: the censoring status is fixed
at creation, which is what keeps AFT out of the SBC matrix; the named enabler
is a status setter. No standalone equivalence scenario (`grouped_aft` carries
it), and `aft-exact.R` is not a MANIFEST entry.

**hazard**. Person-period ingestion sugar: `R/dbarts.R:481-541` expands the
design and remaps the token to `"probit"`/`"logistic"` at :538 before any
model is built, so the fit records `family = "probit"` and the probit row
applies. No engine code, no C-API token. **Judge**: the loglik channel is per
person-period row, not per subject - a real `loo`/`waic` trap, stated nowhere
in the code. The expander refuses a formula interface, a `subset` and a `test`
set by name; `survivalProbabilities()` used to die on hazard fits whose
predictors carried column names, and a fixture now pins it.

**hurdle**. An R-level composition with no sampler of its own - `dbarts()`
refuses it (`R/dbarts.R:447`) and `bart2Hurdle` glues an occupancy probit to a
lognormal positive part at report time. **Judge**: `hurdleLogLik`
(`R/generics.R:2045`) is `log1p(-pi)` at zeros and `log(pi) + dnorm(log y, f,
sigma, log = TRUE) - log(y)` at positives - natural-y scale with its Jacobian,
explicitly **not** the sum of the components' channels. A heteroscedastic
positive part is unreachable because `redirectCall` forwards `variance =` to
the occupancy component too; deliberate, evidenced by the Rd text, stated in
no comment in `R/`.

**bcf / K-forest amplitude family** (`AmplitudeForestCombiner`
combiner.hpp:741): its own stop, Stop 5.

**grouped** (:4706). The in-core replacement for `rbart_vi`'s R-level Gibbs
loop, which diverged 15-91x on tau. An `rbart_vi()`-only surface:
`bartcore.groups` is written at one site (`R/rbart.R:384`) and no other entry
point carries a `group.by` formal. `setWeights`/`setSigma` read off the *base*
family, so grouped probit refuses both while grouped aft takes `setSigma`;
`updateScale` is refused only under a base family with a data-derived
transform. **Judge**: grouped + variance forest constructs, unrefused and
unadjudicated - a `?` cell. `fitted()`/`residuals()` used to segfault on a
dimnames-less random-effect array; an unmatched or all-NA group index is now
refused before the `.Call`, with a C-side bounds check.

**heteroscedastic** (`ConstantVarianceLeaf` :241, `VarianceForest`
chain.hpp:404, built at :4176). A second forest over log s(x), reached by
`variance =` on `bart2`/`dbarts()` (or a `varianceForest()` object,
`R/model.R:1761`); refused on `rbart_vi()`, no `xbart()` reach. `$predict()`
returns a named list `(mean, variance)`. **Two things to judge, both recent.**
First, every report channel used to score at the scalar `sigma` - under this
parameterization a fixed unit residual times the response range, carrying no
posterior content - instead of at s(x); the measured error was a summed loglik
of -3592.1 against -2031.7 correct. `heteroscedasticScale` (`R/generics.R:32`)
now takes precedence in `pointwiseLogLikelihood` (:95) and, via
`ppdNoiseScale` (:2754), in `sampleFromPPD` (:2785), which refuses PPD
sampling where no variance surface was replayed; `summary` reports the
synthetic `mean.s` (`R/diagnostics.R:78`). Evidence
`test-heteroscedastic-channels.R`; the cell stays `?` because the tests pin
the routing, not the formula. Second, `setState` held only mean forests to the
column-mask rule - a reproduction installed 32 forbidden variance splits, 12
still live 5 sweeps later. `columnMaskStateFeasible` (chain.hpp:3418) now
carries a variance pass, `rebuildVarianceForest`'s backstop is at :4377, and
both install arms gate on it (sampler.hpp:946, :1134), each surfacing one
refusal by name (`src/R_interface_bartcore.cpp:7118`, :7433); evidence
`test-heteroscedastic-warm-start.R` and `tests/cpp/test_state.cpp`
`testColumnMaskContainment` :846 / `testVarianceWarmStart` :1184. **Not
covered at ensemble scale by SBC** - deferred not blocked, one of only two
families in that position (aft is the other). `heteroscedastic-exact.R` is the
only exact gate with arms above a single tree (2 and 20), both small-m'.

**Known open across the rows**: `weighted-binary.md` stays parked, post-1.0;
`MonotoneConstantGaussianLeaf` (:507) still forces birth/death-only proposals
and fixed `k = 2`; real `r` for nbinom, the ordinal `log_diff_exp` tail and
the negbin panels are recorded doors; grouped-ordinal and grouped-nbinom block
interleaving are unbuilt doors, not refusals.

---

## Stop 3 - Moves and trees

    git diff main...bartcore -- src/bartcore/moves.hpp src/bartcore/tree.hpp

`moves.hpp` (869): `MoveScratch` :23, `MoveContext` :34, `StepType` :838,
`metropolisJumpForTree` :844. `tree.hpp` (2,158): `InteractionConstraint` :77,
`Rule` :125, `FlatKind` :171, `FlatNode` :190, `Node` :217, `Tree` :257,
`columnMaskSubtreeIsValid` :534 - flat arena unchanged.

**What to look for**: change-move detailed balance and the empty-leaf veto at
`-HUGE_VAL` are both landed and worth reading as diffs with context; the veto
reach gap the second review found on `combiner.hpp`'s per-forest veto-weight
path (`formForestVetoWeights` :1825) is now covered from
`tests/cpp/test_moves.cpp`. `growTreeFromRoot`'s categorical scan landed since
the last refresh (`grow-from-root-categorical-scan.md`, ARC CLOSED), inverting
the v1 "categoricals never split here" contract; with `dart = TRUE`
categorical predictors now receive split counts during grow sweeps instead of
the structurally zero counts they always drew. `n.grow.sweeps` stays opt-in by
design, not omission - `grow-from-root-default.md` is a pre-registered study
that KILLED defaulting it on in both size strata, confirmed by fresh-seed
re-runs. Rule encoding (63-category inline cap, bit 63 doubling as the
missing-value flag) is unchanged.

**Known open**: `setpredictor-leafof-rebuild.md` is CLOSED, roll-only taken.
`tree-mixing-proposals.md` (COMPLETE survey) withdrew its own top pick after a
harm clause fired and ranks a same-variable "perturb" move first among
unimplemented candidates. `monotone-leaf-branch-fill` is a live TODO entry.

---

## Stop 4 - The data layer

    git diff main...bartcore -- src/bartcore/data.hpp

`data.hpp` (1,857). `ColumnStore` moved :215 -> :501 - real growth (typed
ingestion, sparse extensions). Read `docs/design/data-store.md` first; it is a
standing reference required before data-adjacent work. Anchors:
`SparseColumnData` :104, `CscColumnSlice` :121, `ColumnSourceKind` :128,
`ColumnSource` :144, `CodeBlock` :448, `ColumnStore` :501, `ScopedCutGrid`
:1830.

**What to look for**: ownership inversion, cut-point centralization, the
20%-density sparse crossover, and the mutation transaction semantics are
unchanged in shape and are now current-state fact in `docs/architecture.md`
("ColumnStore", "The mutation surface...") rather than tour paraphrase -
prefer that prose; it names the exact enum, `PredictorUpdateResult{accepted,
rolledBack, invalidCutPoints}`. The mutation fuzzer more than tripled:
`tests/cpp/test_fuzz.cpp` is 2,000 lines (was 597), carrying the widened
`OP_GROW` mask, the multi-forest arm, and
`OP_SET_COUNTS`/`OP_SET_CATEGORY_OFFSET`. `typed-ingestion.md` replaced
`SamplerOptions`' 8 ingestion fields and the 7-arg `setTestData` with one
borrowed `PredictorSource` view; no `dbarts.h` change.

**Known open**: `sparse-extensions.md` - in-place nonzero mutation LANDED; the
streaming range kernel and per-column u8 widths stay gated ("gate STANDS,
reopen with numbers"); `setData` on CSC/mixed and per-observation CSC mutation
are recorded typed-ingestion doors. `group-by-exposure.md` stays
RESEARCH-OPEN.

---

## Stop 5 - Multi-forest

    git diff main...bartcore -- src/bartcore/combiner.hpp

`combiner.hpp` (2,085) - more than doubled by `multiplier-combiner.md`'s
generalization from BCF's hardcoded two-forest glue to an arbitrary K-forest
basis/amplitude family. Anchors: `ForestStateData` :35, `ChainStateData` :67,
`MultinomialForestSpec` :382, `ForestCombiner` :525, `AmplitudeForestCombiner`
:741, `shippedShape()` :1478, `MultinomialForestCombiner` :1574.

**What to look for**: `docs/design/multiplier-combiner.md` (LANDED) is the
primary doc; BCF's `a*mu + b_z*tau` is its K=2 instance. Gaussian, probit and
logistic are the reachable families - under probit/logistic the forests
combine into the index with sigma pinned - and aft, ordinal and nbinom are
refused by name at three independent sites (`R/spec.R:639`,
`refusedAmplitudeFamilyReason` `src/R_interface_bartcore.cpp:2265`,
`createAmplitudeSampler` facade.hpp:868). **BCF's own R verb left this
repository** (`bcf-bartcause-relocation.md`, ARC CLOSED); dbarts kept three
R-surface guards on the construction seam and a per-draw per-forest varcount
channel that makes bartCause's contract literal. A BCF fit is an ordinary
`dbartsSampler` via `dbarts(forests = list(forest(vars = ...), forest(basis =
..., vars = ...)))`, and on `bart2` via a `forest()` formula term,
byte-identical to a `forests =` fit at the same seed. `shippedShape()` used to
test only basis width and canonicality, not the half-Cauchy flag, so two
admissible specs got two different models; it now reads the prior too.
`interweaveGlueRidge` (chain.hpp:1331) is landed; the treatment-scale
counterpart stays shelved (`bcf-b-ridge.md`, NO-GO). `bcf-sigma-residual.md`
is RESOLVED: the sigma burn-in transient is slow forest-structure mixing, not
a glue defect, and its `burn = ceiling(72000/thin)` is live in `sbc.R`.

**New since the last refresh**: per-forest replay at new rows
(`Chain::predictPerForestFromSavedSample`, surfaced as `predict(type =
"forest", forest =)` and `$predictForests`, R side `predictForest`
`R/generics.R:635`), and the *blend* that recombines it - `predict(type =
"ev"/"ppd"/"bart")` on an amplitude-coupled fit now forms `eta = shift + sum_k
(glue_k %*% t(B_k)) * F_k + offset`, link after, with a new `bases` formal
after `forest`. Both amplitude arms now require `keepTrees`, which `type =
"forest"` silently did without before. Evidence: `test-predict-blend.R` (390
lines), in-sample identity against `yhat.train` at 1e-12 over six fixture
configs.

**Known open**: `mu.blocks`, per-forest saved-tree replay on multinomial and
grouped BCF stay refused by name; the K-forest batched front door's spelling
is undecided. `updateScale` is refused under *every* family on a multi-forest
sampler, keyed on bases rather than family, "though its transform is the
identity and the re-anchoring the refusal guards against cannot occur" - a
stated divergence from the arc's own plan. Warm start and grow-from-root are
unrefused and *untested* at two forests. Two ticketed residues from the
relocation: `bartcause-subset-pscore` (two pre-existing bartCause bugs,
neither caused by `bcf()`) and `getforestvariablecounts-dimnames`. A bartCause
rebuild against two adjacent dbarts tips showed the same three failing
assertions at both, all `dbarts:::alignForestBasisToSubset` refusing a
partial-data `basis` - a bartCause-side follow-up, which is why a grep here
finds nothing. Calibration: BCF's sigma, flagged at ensemble size 50, is clean
at 200, and the evidence is gaussian-only - no equivalence scenario and no SBC
arm reaches a latent K-forest.

---

## Stop 6 - The bridge

    git diff main...bartcore -- src/R_interface_bartcore.cpp \
      src/R_interface_bartcore_common.hpp src/R_interface.cpp

`R_interface_bartcore.cpp` is 7,759 lines. Entry map: `bartcore_create` :3489,
`_createFromHandle` :3567, `_createBCF` :3739 (no `_createMultinomial*` -
retired); `_run` :4241, `_runWithCallback` :4485, `_growFromRoot` :4566;
`_setOffset` :4583, `_setResponse` :4603, `_setSigma` :4627, `_setData` :4634,
`_setWeights` :4869; `_setPredictor` :5053, `_updatePredictor` :5101,
`_setActiveRows` :4022; `_setForestBasis` :3929, `_getForestAmplitudes` :4052,
`_getForestVariableCounts` :4214; `_getCalibration` :4139, `_setCalibration`
:4194; `_storeState` :5618, `_setState` :5623, `_installForests` :5636;
`_predict` :5771, `_predictPerForest` :5872, `_getTrees` :5956, `_getLatents`
:6105. Guards: `refusedAmplitudeFamilyReason` :2265,
`refuseMultiForestMutation` :2613, `refuseMultiForestResponseMutation` :2636,
`refuseBinaryWeightChange` :2773, `refuseUndefinedTestFits` :2867,
`refusePinnedSigmaChange` :2898.

**What to look for**: **the facade had no conformance test - now it does.**
The second review's sharpest finding was that 5 of 7 planted forwarding
defects in `facade.hpp` passed the whole C++ suite, because every other test
drives `Sampler<L>` directly; `tests/cpp/test_facade.cpp` (1,202 lines) is the
answer, one row per `SamplerBase` virtual, driven through the base reference.
Then read `refusePinnedSigmaChange`'s own comment (:2885-2897) - the clearest
statement in the tree of why a guard is keyed on family rather than on an
internal flag.

**Known open**: PROTECT balance and leak safety rest on manual local runs, not
on a repeatable gate (Stop 10).

---

## Stop 7 - The shipped C ABI

    git diff main...bartcore -- inst/include src/C_interface.cpp \
      inst/tinytest/capi inst/tinytest/test-capi.R

`dbarts.h` is 1,158 lines (was 477 - nameable calibration's flat-C half, the
dbarts.h reshape, and threaded predict all landed here); `C_interface.cpp`
1,175. Constants: `DBARTS_C_API_MAJOR 1` (:103), `DBARTS_C_API_MINOR 0`
(:104), `DBARTS_C_API_HASH 0x66d33f1613892406ULL` (:142).

**What to look for**: **threaded predict** (`docs/design/threaded-predict.md`,
LANDED; section 11 is its landing note). `predict`'s `n.threads` formal
existed before but was inert; it is now wired, since dbarts is pre-1.0-0 and
the ABI is free to change. Partition axis: (chain, draw) slabs with no
cross-slab accumulation, so replay is bitwise identical at every thread count.
`dbarts_sampler_predict` gained `size_t numThreads` as its new last parameter
(:861); `numThreads = 0` means the sampler's own count floored at
1. **Version constants did not move** - the header's rule is that they do not
   move before the first release, and a signature that gains a parameter is
   breaking rather than additive, so only the two hash literals were re-baked.
   Consumer migration: stan4bart's one call site gains a `0`; bartCause needs
   zero lines, because the formal is appended **last** on all six generics
   specifically so its positional call stays correct. `test-capi.R` grew 379
   -> 1,759 lines; its compiled `consumer.c` had been passing `0` at the *old*
   arity and silently skipping 264 assertions since the prototype changed -
   caught and fixed by this landing. `bart-as-a-component.md` section 6
   ("Graduation debt") lists what the flat surface still lacks against R5: no
   per-observation predictor update or joint session, no `setCutPoints`/
   `setData`, no `predictVariance`, no forest-indexed `predict`.

**Known open**: `run()`'s own per-call thread override stays inert by design
(door D1 - `run`'s count also sizes `routeTestRows` and mutates state).
Interrupt polling during a long threaded predict is deferred (door D3), judged
safe because predict never calls `R_CheckUserInterrupt` at all. **Reported,
not fixed**: a fractional `n.threads` truncates silently
(`validatePredictThreads` `R/generics.R:234` - `is.numeric(2.7)` passes every
predicate and `as.integer` truncates); a throw from `std::thread` construction
would leave earlier workers unjoined, a hole `run()` shares and no D-number
records; `Rf_error` after the join skips frame destructors. The scaling-curve
gate the design calls for, `benchmarks/R/bench-predict.R`, **does not exist**
(grep-confirmed absent) - the only numbers on record are a same-process probe
in the landing note (1.92x at 2 threads, 3.72x at 4), explicitly "information
rather than a gate". Whether predict's R-level default stays the fit's own
thread count, and the tuned cutoff constant (`Sampler::predictParallelCutoff`,
1e7 traversals), are both stated open. Four selecting control attributes
(`bartcore.n.categories`, `bartcore.dispersion`, `bartcore.groups`,
`bartcore.variance`) are reachable through `dbarts_sampler_create` and
documented nowhere in the header.

---

## Stop 8 - The R surface

    git diff main...bartcore -- R NAMESPACE man

+22,656/-3,227 across 55 files. Biggest: `R/generics.R` 2,985, `R/bart.R`
2,896, `R/dbarts.R` 2,093, `R/model.R` 1,812, `R/data.R` 1,729, `R/rbart.R`
1,309, `R/bartcore.R` 1,107.

**Own-class families and their generics**: `bartMultinomial`/`bartOrdinal`/
`bartNegbin`/`bartHurdle` each get `extract`/`fitted`/`residuals`/`predict`/
`print`, and now **`plot`**, **`extract(type = "loglik")`**,
**`as_draws_array`/`as_draws_df`**, and real **`summary`** methods
(`R/diagnostics.R:302`, :313, :325, :390; `print.summary.bart` :429).
`combineChains` and `ci.level` are honoured on those fits, and 13 further
argument swallows were closed. The reviewer's independently coded loglik
oracles agreed to <= 3e-15 on all four families.

**Family vocabulary, offset shape, removed dots**: `bart()` refuses all ten
`bart2` family tokens **by name**, echoing the token typed
(`bartRedirectedFamilies` `R/bart.R:2621`, `refuseBartRedirectedFamily`
:2634). One `offset` argument everywhere, shape enforced by family
(multinomial requires n x K, every other family refuses a matrix);
`dbartsData`'s separate `offset.category`/`offset.category.test` R-facing
formals are gone, though the slots keep their names. `bart2()`/`rbart_vi()`
lose their rejection-only `...` entirely, so R's bare `unused argument` is now
the wall for those two entries. `checkFamilyUnsupportedArgs` (`R/bart.R:620`)
is the single site refusing `samplerOnly`/`warm.start`/`n.grow.sweeps` per
family.

**Other landed items**: `$n.chains` is set unconditionally on every fit that
keeps a sampler. The 31 low-level `bartcore*` handle wrappers moved out of the
namespace to `inst/common/bartcoreHandle.R`, 559 call sites across 41 files
prefix-stripped; that file ships in the tarball by design.
`man/dbartsSampler-class.Rd` documented three wrong RC formals
(`control`/`model`/`data` for `newControl`/`newModel`/`newData`, and
`storeState` missing `ptr`) - fixed, and `tools/check-rc-codoc.R` (259 lines,
wired into `lint.yaml`) now fails CI on a renamed RC formal, which
`tools::codoc()` structurally cannot see because it walks S3/S4 methods rather
than a `setRefClass` generator.

**Known open**: **`predict(bases =)`, `predict(sample =)` and
`fitted(combineChains =)` on the own-class fits are still silent no-ops** -
explicitly outside the judgement table, and distinct from the
`combineChains`/`ci.level` honouring that same pass DID add; do not assume the
generics work closed every swallow. `setForestBasis(k, ~var)` still evaluates
the formula in `environment(basis)` rather than against the sampler's data, so
`~x3` naming a predictor raises a raw `object 'x3' not found`. `plotTree`'s
dead padding branch stays out ("no judgement names it"). `dbartsData(counts
=)` with a formula still reaches its own later refusal. **`par(mfrow)`
restoration is now partial**: the five own-class plot methods, `plot.pd2bart`
and `plotTree` save and restore, but `plot.bart` and `plot.rbart` set `mfrow`
through the shared `plotSigmaTrace` (`R/plot.R:9-12`) and still leak it.

---

## Stop 9 - The capability axes

The axes that cut across every row of Stop 2. `feature-matrix.md`'s five
tables are the authority; this is the narrative for the ones needing a
judgement rather than a check.

**DART** (`DartPrior` model.hpp:2433). New in 1.0-0 - nothing on `main`
carries it: a Dirichlet prior over split-variable probabilities inducing
variable selection, with optional sampling of the concentration. Two
deliberately colliding spellings: the prior object `dart(power, base, a, b,
rho, alpha, update.alpha, update.delay)` (`R/model.R:1633`, exported through
`dbartsPriors`) passed as `tree.prior`, or the `dart = TRUE` shorthand on
`bart2`/`rbart_vi`/`xbart`; the ladder is `resolveTreePrior` (`R/bart.R:484`)
and supplying both is refused. `update.delay` defaults to half of `n.burn`,
the BART package's `startdart` convention (`docs/design/prior-defaults.md`).
Runs return per-sample probabilities as `varprobs` beside the variable counts.
**Judge**: DART is refused on the amplitude/multi-forest path
(`refuseUnsupportedAmplitudeComposition`, `src/R_interface_bartcore.cpp:2323`)
and by name on multinomial; the variance forest never takes it while the mean
forest does; and there is **no SBC arm** - sampled `k`, DART and tau are out
of scope. The one dedicated test is
`test-dart-mixed-columns.R` (76 lines), covering `varprobs` indexing across
pooled-categorical and NA-incorporate column mappings, otherwise unguarded.

**setSigma** (R5 `R/dbarts.R:1504`, bridge :4627, flat C
`dbarts_sampler_setSigma`). `refusePinnedSigmaChange`
(`src/R_interface_bartcore.cpp:2898`) refuses two cases: a variance forest
owns the residual scale, so there is no single sigma to set; and any family
outside `{gaussian, aft}` fixes the scale by definition. Keyed on family
rather than on the internal `sigmaIsFixed_` flag *deliberately* - a gaussian
sampler with `resid.prior = fixed()` also pins sigma, and driving it per sweep
is the supported outer-Gibbs idiom. Ladder pinned at
`test-sampler-errors.R:128-200`.

**Starting sigma**. `estimateSigmaFromLinearModel` (`R/utility.R:542`),
wrapped by `estimateStartingSigma` (`R/spec.R:24`); the user override is
`sigest`. Three recent fixes, innermost out: a non-finite linear-model sigma
(typically zero residual df) now routes to the marginal residual sd with a
classed `dbartsSigmaFallbackWarning` rather than erroring; the result is
floored at `sqrt(.Machine$double.eps) * max(1, max|residual|)` (:617), so a
perfectly reproduced response no longer installs host-BLAS rounding noise - it
used to land on exactly 0 on SSE2-class kernels, and the sampler refuses
non-positive sigma, so `dbarts()` succeeded or failed by hardware; and the
finiteness test is re-checked *after* the fallback (:618), because `sd()` of a
length-1 residual is `NA`. **Worth an eye**: the sparse branch returns at
:562, before the floor, so a CSC-backed design with a degenerate residual
still reaches the sampler's own refusal.

**setResponse(updateScale =)**. Default **FALSE** (`R/dbarts.R:1301`), locking
the creation-time transform - the safe embedded-Gibbs behavior, and the
resolution of the drifting-leaf-prior defect `range-scaling.md` records.
Latent families no-op the flag rather than refusing it (a `-` cell, not `R`);
it is refused on a multi-forest sampler under every family, and on grouped
only under a base family with a data-derived transform. **The flat C API
passes `updateScale = true`**, so `dbarts.h` consumers keep historical
re-anchoring bit-identically. Open: no equivalence scenario and no RNG-locked
test exercises mid-run `setResponse`.

**Zero-weight rows and `setActiveRows`** are two mechanisms, easy to conflate.
Zero *weights* are accepted, not refused, on the families that take weights,
and the sigma posterior uses `nu_0 + N_positive` rather than `nu_0 + N`
(`drawSigmaSqFromPosterior` model.hpp:2540; `test-zero-weights.R`). The mask
is the general mechanism (`R/dbarts.R:1398` -> bridge :4022 ->
`Chain::setActiveRows` chain.hpp:1666 -> facade.hpp:368, flat C
`dbarts_sampler_setActiveRows`), `ARC COMPLETE`, `supportsActiveRows() ==
true` on every family: an inactive row leaves every sufficient statistic,
every family-level parameter update and its own latent draw, but keeps its
leaf occupancy and its fitted value. Zero-weight rows report `NaN` from the
loglik channel, not `-Inf`, deliberately. Open: the masked-vs-compacted
equivalence oracle is piloted, not shipped, and must use *different* seeds -
same-seed runs are bitwise-locked and the assertion is vacuous; zero-weight
occupancy still does not match a compacted fit; two `-Inf` sites at zero user
weight remain unguarded in the engine's pointwise log-likelihood.

**Warm start** is three distinct things, and `warm.start` names only the
first. (1) Donor-forest install: R5 `installTrees(donor, samples)`
(`R/dbarts.R:1888`) and `bart2(warm.start =)`, bridge :5636, engine
`Sampler::installForests` (sampler.hpp:983). It seeds live trees, sigma and
`k` **only** - rng, latents, group effects and saved buffers stay fresh: a
starting position, not a continuation. Tree-count and DART mismatches are
refused by name; cross-grid donors remap under a scoped cut-grid swap. (2)
Grow-from-root: `growFromRoot(n.sweeps)` (`R/dbarts.R:1013`) /
`bart2(n.grow.sweeps =)`, refusing linear and GP leaves; its memo's NO-GO
stands - warm start only, never a standalone posterior sampler. Combining (1)
and (2) errors. (3) Prior draw: `sampleTreesFromPrior` (`R/dbarts.R:985`),
whose law was corrected to the CGM prior *conditioned* on the no-empty-leaf
set. Open: `rbart_vi()` carries neither formal - a surface gap, not an engine
one - and neither do ordinal, nbinom, multinomial, hazard or hurdle.

**Pointwise loglik**. Six methods carry an `extract(type = "loglik")` channel:
`extract.bart` (`R/generics.R:425`, covering gaussian, student,
heteroscedastic, probit, logistic and aft through the shared
`pointwiseLogLikelihood` :56), `.bartMultinomial` :966, `.bartOrdinal` :1297,
`.bartNegbin` :1576, `.bartHurdle` :2001, `.rbart` :2379. It is extract-only
and train-only by design; `predict`/`fitted` refuse it, and an unrecognized
family errors rather than being scored with a formula that does not fit. The
engine-side per-family `Results::logLikelihood` channel is deferred behind
`c-api-growth`.

**Nameable calibration** (`prior.scale`) and **interaction constraints** are
both `ARC COMPLETE`; `$getCalibration`'s docstring (`R/dbarts.R:1731`) is the
authoritative statement of what its eleven columns mean and when the two
`node.scale` columns go `NaN`.

---

## Stop 10 - The gate apparatus

    git diff --shortstat main...bartcore -- tests/cpp inst/tinytest benchmarks .github

Two adversarial whole-branch reviews and a calibration lane sit on top of the
mechanical apparatus, and their headline finding is that green checks prove
less than they look like they prove.

**tests/cpp**: 13 `test_*.cpp`, 26,424 lines (26,842 with the harness), led by
`test_model.cpp` 7,291, `test_sampler.cpp` 6,841, `test_fuzz.cpp` 2,000,
`test_state.cpp` 1,858, `test_moves.cpp` 1,759, `test_grow.cpp` 1,718,
`test_facade.cpp` 1,202 and `test_data.cpp` 1,177.

**inst/tinytest**: 169 files, 43,098 lines. Pin files worth knowing by name:
`test-capi.R` 1,759, `test-argument-surface.R` 923, `test-active-rows-pins.R`
572, `test-pointwise-loglik.R` 490, `test-predict-blend.R` 390,
`test-generics-multithreaded.R` 373 (`identical()` not tolerance, thread
counts 1-4, every family), `test-state-weight-pairing.R` 287,
`test-xbart-fold-oracle.R` 284 (folds reconstructed outside `xbart` from the
seed; its leakage band was recalibrated after a 20-seed sweep found a 10%
false-alarm rate), `test-tree-store-order.R` 235.

**benchmarks**: 35 harnesses (`sbc.R` 2,603, `equivalence.R` 1,913).
`benchmarks/baselines` was pruned to **11 files plus the MANIFEST**: four
`role: current` (`equivalence-736bfb05.rds` 43 scenarios,
`bcf-equivalence-6e3b9fb8.rds` 12, `multinomial-equivalence-4d9a3337.rds` 11,
`bench-sampler-ab1dc52.csv`), six `historical` still reachable by a workflow
or harness, and `equivalence-5430fdb.rds`, the unregenerable classic-engine
cutover evidence. **The P17 rule** (MANIFEST header): a re-record that changes
any recorded draw must name the ORACLE establishing the new values are right -
"'the deviation looks like the expected class' is an adjudication, not an
oracle." The three equivalence current rows honour it explicitly; the bcf row
states its oracle informally.

**CI** (11 workflow files; only `check-standard.yaml` exists on `main`). Six
are push-triggered on `bartcore` and do run: `check-standard`, `cpp-tests`,
`exact-gates` (quick), `lint`, `pkgdown`, `sanitizers`. **Five are `schedule`
+ `workflow_dispatch` only - `equivalence`, `rchk`, `sbc`, `valgrind`,
`revdep-smoke` - and have never run, not once, by any trigger, on any
branch**: GitHub binds both triggers to the default branch, which does not
carry these files, so they cannot even be hand-dispatched. Every clean rchk /
valgrind / equivalence / SBC claim rests on a **manual local run**: rchk
clean; valgrind clean ("All ok, 6419 results") after fixing one real leak and
one pre-bartcore OOB read in R's own `makeModelMatrixFromDataFrame.c`;
gaussian SBC PASS; the calibration run below. **Read
`docs/plans/release-candidate-review.md:34-38` and `:59-63` with a
correction**: they state as standing fact that "ten of the eleven workflow
files ... have never been registered". That was the plan's premise and it is
wrong for the five push-triggered ones - five, not ten, are dormant. The tip
you are reviewing is itself docs-only and fired no workflow, `lint.yaml`'s
`paths-ignore` covering `docs/**`, `TODO`, `**.md`.

**exact-gates** runs 20 scripts. Most are `n.trees = 1` by construction - that
is how an exact posterior is available at all - the exceptions being
`heteroscedastic-exact.R` (2 and 20), `multinomial-exact.R` (50 and 75),
`logistic-reference.R` (50 in quick, 200 in full), `hazard-reduction.R` (40)
and `hurdle-reduction.R` (30). `logistic-reference.R`'s ensemble arm **does**
gate CI: `exact-gates.yaml` installs `any::BART` for it and the loop exits
nonzero on failure. Full mode has never run in CI, only quick.

**The two reviews.** The first (2026-08-18, no commit, six independent
reviewers) found ~10 BLOCKERs / ~19 MAJORs / ~25 MINORs, all fixed in the wave
that followed; its own `geweke-mc.R` fired a real defect
(`sampleTreesFromPrior` projected the unrestricted CGM prior through
`collapseEmptyNodesBelow` where the moves price it restricted and renormalized
to the empty-leaf-free set). The second (2026-08-24) confirmed **8 BLOCKERS,
44 MAJORS, 35 MINORS; 5 refuted; 6 qualified; 6 narrowed**, after a
whole-suite replant of 31 "zero-killer" mutations found 16 real gaps the
per-leg file scope had hidden and 15 false alarms. All nine VD-judgement
groups are decided and landed. Its five most consequential defects, all fixed:
`extract(type = "trees", sample = 1)` silently returning a different sample's
trees; `survivalProbabilities()` dead on named-column hazard fits;
`fitted()`/`residuals()` on `rbart` segfaulting on a dimnames-less
random-effect array; `dbartsData()`'s false row-count message; the facade
conformance gap.

**Calibration lane**: `sbc.R` re-run at `n.trees = 200` - bracketing both
shipped defaults, 75 and the classic 200 - against the 50 every prior SBC run
used; Bonferroni band over 83 functionals. **11 families carry an
ensemble-scale verdict, 75/83 functionals PASS, 3 flag, all 3
pre-adjudicated** (nbinom's `r`/`agg.psi` ridge; a grouped-gaussian sigma
harness artifact confirmed three independent ways). Multinomial's raw `f_ik`
cells, undischargeable at m=50, are clean at 200; BCF's sigma clears. **Not
covered at ensemble scale**: aft and heteroscedastic, both with real sampling
code that does not reduce to a covered family; hazard and hurdle are
calibration-inherited by bitwise reduction, not directly SBC'd.

**Known open**: 11 of 35 `benchmarks/R` harnesses and two of the five `tools/`
scripts are wired to no workflow; five harnesses call themselves "repeatable"
or "permanent" gates in their own headers yet have exactly one recorded
(manual) run each - `mutation-battery.R`, `grouped-mixing.R` (whose
re-measured IACT numbers already diverge from the figures its header cites,
undetected because nothing re-runs it), `backfit-exact.R`, `geweke-mc.R`,
`composition-matrix.R`. All baselines are single-host (arm64 macOS); the x86
bench box is unavailable; and the bcf equivalence baseline predates the
`summaries` field, so the most expensive dormant gate will degrade loudly on
its first live cross-host run. Reproducibility is within-host bitwise across
SIMD dispatch only.

---

## Stop 11 - Build sediment

    git diff --stat main...bartcore -- configure configure.ac tools src/misc \
      src/external src/include src/Makevars.in

Skim, do not read. The `configure`/`configure.ac`/`tools/m4` deletions are
`autoconf-dead-code.md`'s. `tools/` holds five R scripts -
`check-doc-freshness.R` 450, `check-rc-codoc.R` 259 (new, wired into
`lint.yaml`), `regenerate-snapshots.R` 134, `check-win-drift.R` 128,
`check-build-freshness.R` 58 (inherently local-only) - and `tools/` is not
`.Rbuildignore`d, so all five ship in the tarball. `simd.c`'s
AVX2-misdetected-as-AVX fix (`x86-simd-plan.md`) is the one thing here worth
an actual look; `src/external/randomBase.c`'s per-field `.Random.seed`
validity check, inherited from the 0.9-34 CRAN patch this branch rebased onto,
is the other.

---

## Appendix: recent arcs, entry anchors, evidence

Everything landed since 2026-08-20, newest first: one place to start reading
and the gate that holds it. The arcs above the rule have **no landing note
anywhere** - their only record is the commit message body, which is unusually
full; read them with `git log a5d0d43b..ae5b91d8`.

- NEWS coherence (`ae5b91d8`): three self-contradicting 1.0-0 pairs removed
  (the binary-`k` default listing `xbart` after it was carved to fixed `k =
  2`; the dots described as rejection-only *and* as removed; a multinomial
  `offset.category*` formal that never existed). Entry `inst/NEWS.Rd:5-77`.
- Sanitizer deps (`3298da0e`, `sanitizers.yaml` installs `survival` and
  `posterior`, so seven `requireNamespace`-gated files stop silently skipping
  under ASAN/UBSAN); RC Rd, pdbart cost and `par()` (`654e4681`, entry
  `tools/check-rc-codoc.R` and `R/plot.R`, evidence the tool's own
  planted-rename check); baseline prune (`01bd49a7`, 65 -> 11 files, entry
  `benchmarks/baselines/MANIFEST`); backward-facing strip (`39fa379a`,
  `9a908bbd`, `TODO` 583 -> 329 lines).
- Anchor re-derivation (`018225ca`, `9f538bbf`, `60aea703`, `89e1ecd3`,
  `207d4ca0`, `5a7afe57`, `69cb6491`; `ffa13a05`, `27c63859`, `b2e360c6`,
  `ba0e79ba` for the matrix; `6e148281` for prose): every code anchor in 29
  `docs/design` files re-derived by content, after an audit measured a 42
  percent non-landing rate on a 67-anchor sample - coverage, not method, was
  the failure. `docs/plans` was not swept.

---

- **Citation strip** (`dfb6dc0a`, `c7efa857`, `ad4d801d`): 316 `docs/`-path
  comments removed from 81 files, each constraint restated in place. Entry
  `inst/include/dbarts/dbarts.h`; the census rule `grep -rn "docs/" R src
  inst/include inst/tinytest inst/common man | grep -v http` is empty.
- **xbart k-grid sort** (`736bfb05`, `e73ce948`): the rewrite had dropped
  0.9-34's decreasing sort and un-permute, so loss depended on the order `k`
  was listed. Entry `R/xbart.R`; evidence an order-invariance pin in
  `test-xbart-method.R` plus the third baseline re-record (42/43 bitwise,
  xbart alone moves, max |z| 1.94).
- **Threaded predict** (`c44fcbc5`, `e9be0c7a`, `b2df6522`): entry
  `Sampler::predictColumns` sampler.hpp:705 (fan-out :602); record
  `docs/design/threaded-predict.md` section 11.
- **Constant-response sigma** (`ff27b1ba`, `afb2b80a`, `e0e59097`): entry
  `estimateSigmaFromLinearModel` `R/utility.R:542`; evidence
  `test-boundary-inputs.R:25-82`, pinned with `tolerance = 0` after
  `all.equal`'s default had made the old 1e-15 pin vacuous.
- **Fix wave 3** (`42b12ac7..7ad0bbea`): six independently gated slices -
  exhaustive family switches; the three engine divergences (`cd615a1e`); xbart
  fold oracles; the R-surface vocabulary pass; the own-class generics pass;
  the handle-wrapper move. tinytest 7040 -> 7181, baselines unchanged,
  `dbarts.h` untouched.
- **Second review, waves 1-2** (`b102e17c..07ad73e4`): entry
  `tests/cpp/test_facade.cpp`; evidence `consolidated-report.md` and
  `gate-ledger.md` under `docs/plans/review-2026-08-24/`.
- **Heteroscedastic report channels** (`221ec7af`): entry
  `heteroscedasticScale` `R/generics.R:32`; evidence
  `test-heteroscedastic-channels.R`.
- **Variance column mask on setState** (`c95a5e83`): entry
  `columnMaskStateFeasible` chain.hpp:3418; evidence
  `test-heteroscedastic-warm-start.R`, `tests/cpp/test_state.cpp:846`.
- **Ordinal missing-row scan** (`b4b9119d`, `044a9098`): entry
  `scanOrdinalCuts` scan.hpp:105; evidence an independent R re-derivation to
  1.53e-16 and `tests/cpp/test_scan.cpp:337` `testMissingRouted`.
- **Saved-tree store cursor, `pdbart(keeptrees=)`** (`124259d0`, `41661523`):
  entry `savedSlotForDraw` sampler.hpp:520 and
  `pdbart.getAndInitializeSampler` `R/partialDependence.R:1`; evidence
  `test-tree-store-order.R`, `test-pdbart-keeptrees.R`.
- **`dbartsData(bases =)` alignment** (`47cdb96a`): entry
  `validateForestBases` `R/data.R:721`; evidence `test-forest-basis-subset.R`.
- **Summary methods** (`1583140b`): entry `summary.bartOrdinal`
  `R/diagnostics.R:302`; evidence `test-summary-nondefault-families.R`.
- **Type refusals, multinomial replay decision** (`99317fec`): entry
  `validateType` `R/generics.R:1833`. Per-forest off-sample replay on
  multinomial is a **decided refusal** - every candidate use collapses to
  `log(predict(...))` to 1e-13 or needs a declined identification; re-open
  trigger is a level-sensitive K-latent composition.
- **Weights digest on saved state** (`da82ec23`): entry `Chain::weightsDigest`
  chain.hpp:1611; evidence `test-state-weight-pairing.R`.
- **`OP_GROW` fuzz reach** (`d788bfef`), **multinomial fuzz arm**
  (`53525f4d`): entry the `tests/cpp/test_fuzz.cpp` op table. Widening the
  mask alone was vacuous - the op ends in a getState/setState round trip that
  heals corruption first - so the check runs before that round trip.
- **Predict-time blend** (`139a1976`), **per-forest replay** (`63df524e`):
  entry `predictForest` `R/generics.R:635` and
  `Chain::predictPerForestFromSavedSample`; evidence `test-predict-blend.R`,
  `test-predict-forest.R`.
- **Logistic weight channel** (`d0701a6a`), **engine binding order**
  (`a0eaf348`), **multinomial mutation arc S0-S4** (`5e586587`, `b2d1749f`,
  `5a3bc276`, `2619ac9e`): records in
  `docs/design/multinomial-mutation-arc.md`, `docs/design/r-c-division.md`.

## Appendix: gate commands

    R CMD INSTALL .             # --preclean after headers/Makevars.in/configure.ac
                                 # ALWAYS --preclean after facade.hpp virtuals
    cd tests/cpp && make && ./test_bartcore
    Rscript -e 'tinytest::test_package("dbarts")'

    # equivalence - NEVER run in CI; this is the only pre-landing check
    Rscript benchmarks/R/equivalence.R compare \
      benchmarks/baselines/equivalence-736bfb05.rds

    # exact-posterior gates (quick is what CI runs; full mode never has)
    Rscript benchmarks/R/bcf-exact.R quick

    # speed - never concurrent with other load, never in CI either
    Rscript benchmarks/R/bench-sampler.R compare \
      benchmarks/baselines/bench-sampler-ab1dc52.csv

## Appendix: where decisions live

- `docs/architecture.md` - current state, not history; prefer its prose to
  this tour's paraphrase wherever they overlap.
- `docs/design/feature-matrix.md` - the present-state matrix: 13 family rows,
  26 capability columns across five tables, 49 footnotes, and a Gaps section
  (:924) collecting every `M` and `?` cell as a candidate work item. The one
  document that surfaces every family and capability gap in one place.
- `docs/design/INDEX.md`, `docs/plans/INDEX.md` - complete manifests with
  status, NO-GO and closed items included.
- `docs/plans/README.md` - process, RNG gate classes, and the
  symbol-over-line-number citation rule.
- `docs/plans/release-candidate-review.md` (3,831 lines) - the master log. New
  entries are inserted right after the `## Landing notes` header (:541) rather
  than appended at EOF, so it reads newest-first. Sections 1-2 are the plan's
  premises at its writing, not current state; see Stop 10 on the dormancy
  claim.
- `docs/plans/review-2026-08-24/` - the second review's artifacts:
  `consolidated-report.md` (the ledger), `matrix-results.md` and
  `matrix-review-{entries,generics}.md`, `gate-ledger.md`,
  `mutation-{A,B,C,D}-findings.md`, `reading-{engine,R}-list.md`,
  `calibration-sbc.md`, `decision-brief.md`, `generics-survey.md`,
  `generics-phase2-spec.md`, `wave3-plan.md`, `xbart-oracle-memo.md`; the
  per-slice memos and their blind critiques are under `memos/`.
- root `TODO` - open work only, alphabetical, 329 lines.
  `second-review-followups` and `rc-gate` are the "what is left" pointers; the
  `release` section is the one ordered procedure.

Statuses that mean "read this before assuming a gap": NO-GO (measured and
rejected, with numbers), CLOSED/SUPERSEDED, RESEARCH-OPEN (gated on inability
to name enabling value, not on an absent caller), REFERENCE (standing doc, no
status line by design), DECIDED (a judgement group resolved it - the decision
is the record, do not relitigate).
