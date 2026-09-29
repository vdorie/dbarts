# multiforest-leaf-prior-writer: setLeafPrior on multi-forest samplers

Status: PLANNED 2026-09-29. The Decision below was settled by the orchestrator as an agent-made call (dec-A127): the reader reports forest(sd = ) on map forests, and creation refuses forest(sd = Inf) as the writer does. Blind critique next, then implementation.

agent: opus (engine, facade, bridge and R in one worktree; one writer)
rng: neutral (no draw moves for a sampler that never calls the writer; the constructors are untouched)
budget: ~700 lines: engine ~70, facade ~20, bridge ~45, R ~150, man/NEWS/docs ~65, tests ~350 (R ~220,
tests/cpp ~130); plans have run 1.5-2x low, so expect 1000-1400. The ruling's ~450 omits the reader
change and the identity tests. stan4bart and bartCause 0.

## Goal

`$setLeafPrior` accepts on multi-forest samplers what their creation accepts (dec-B142): on a
multinomial sampler, `normal(k = )` with a fixed k, applied to every category forest; on a sampler
that carries forest amplitudes, a forest's `sd` stated as at creation, `forests = list(forest(sd = ),
...)`. The single-forest contract holds: a write equal to what is in force is bitwise inert, takes
effect on the next sweep, reinterprets no drawn value, and survives re-creation.

## Decision (reader on a map forest)

On a forest whose scale the calibration map sets, `getLeafPrior()$leaf.prior` is today
`normal(sd = )` at the map's leaf spread in response units (dec-A124, not yet ruled). Neither creation
nor this writer takes that on such a forest, so the ruled round trip cannot hold as stated.

- Recommended: report creation's own term there, `forest(sd = )`: the half-Cauchy median
  (`amplitude.prior.scale`) on a forest with no basis, `leaf.scale.factor` on a basis forest.
  `prior.sd.of` reads `"amplitude scale"` or `"forest total"`. Nothing is lost: k is pinned at 1, so
  the response-unit spread is `anchor`. Writing every forest's `leaf.prior` back through
  `forests =` is then inert. Cost: ~25 R lines and the reader tests' map-forest expectations.
- Alternative: keep the reader and have the writer also take `normal(sd = )` per forest. That is a
  second spelling creation never takes, and on a forest with no basis it would set a leaf spread that
  creation cannot set, since that forest's `sd` is its amplitude's median.
- What would change it: a consumer that reads `leaf.prior`'s sd on a map forest. None exists in
  dbarts; bartCause's bcf reads `response.scale` and `response.shift` only (confirm at the sweep).

## Context

- Today [`setLeafPrior`](../../R/dbarts.R) refuses both samplers through
  [`refuseCountsMutation`](../../R/bartcore.R), [`refuseAmplitudeMutation`](../../R/bartcore.R), and
  the engine refuses through [`Chain::setForestPriorScale`](../../src/bartcore/chain.hpp) (false
  whenever a combiner exists). The single-forest contract is in the [Landing](leaf-prior-k-or-sd.md#landing)
  note of leaf-prior-k-or-sd.md. The reader's shape is in [Reader shape](leaf-prior-reader-shape.md#reader-shape).
- Multinomial: [`Chain::buildMultinomialForest`](../../src/bartcore/chain.hpp) builds every category
  forest from [`MultinomialSpec`](../../src/bartcore/combiner.hpp)'s constant anchor and the host k;
  the per-leaf sd `leaf.scale / k` is read live by the leaf draws and by
  [`MultinomialForestCombiner::afterCombine`](../../src/bartcore/combiner.hpp). Nothing caches k.
- Amplitude coupling: the map is described in
  [The calibration map, general in K](../design/multiplier-combiner.md#the-calibration-map-general-in-k).
  [`forestParams`](../../R/model.R) sends a forest's `sd` down one of two channels. With no basis, it
  is the half-Cauchy median ([`ForestAmplitudePrior`](../../src/bartcore/combiner.hpp)'s
  `halfCauchyScale`, echoed in `amplitudePriorScales_`), and the node scale stays at the anchor s.
  With a basis, it is `nodeScaleFactor` in `factor * s / (0.674 * c)`.
  [`Chain::setForestBasis`](../../src/bartcore/chain.hpp) already re-derives that leaf scale from the
  retained `nodeScaleAnchor_`, and re-imposes the map by setting `nodeScaleIsMapDerived_`.
- State: [`ForestStateData`](../../src/bartcore/combiner.hpp) carries each forest's k and leaf scale,
  not the half-Cauchy scale (the saved variance is the live auxiliary); see
  [`Chain::noteInstalledLeafScale`](../../src/bartcore/chain.hpp), [`Chain::adoptInstalledAmplitudePriors`](../../src/bartcore/chain.hpp).
  Creation reads the per-forest sd from `attr(control, "bartcore.forests")$params`
  ([`applyAmplitudeSpec`](../../src/R_interface_bartcore.cpp)); [`bartcoreSamplerSetResponse`](../../R/bartcore.R)
  already mirrors a mutation into a control attribute so re-creation reads it.

## Settled

1. Surface: `setLeafPrior(leaf.prior, forests = NULL, updateState = NULL)`, exactly one of the first
   two. Creation states per-forest spreads only through `forests = list(forest(sd = ))`, so the
   writer takes the same argument name, the same constructor (resolved by
   [`evalInForestVocabulary`](../../R/family.R)) and the same positional correspondence. One call
   can restate several forests, all validated before any engine write, and `getLeafPrior()`'s
   per-forest list maps onto it. Rejected: a `forest` index as on `setForestBasis`, which puts a
   `forest()` in `leaf.prior`, a slot creation never uses for one. A short list reaches the first
   forests, as at creation; an undeclared `sd` leaves its forest as it is (a write, not a creation).
   `setLeafPrior` never reached main, so inserting `forests` second breaks no released call.
2. Multinomial: nothing in the softmax map is recomputed. The anchor and the leaf scale stay; only
   each forest's k moves. The model field records the write (via [`restateLeafPrior`](../../R/dbarts.R)),
   so `getPointer` re-creation and `copy` build with it. A fixed k is not "drawn", so the rule that a
   drawn k keeps its value across an anchor change does not apply.
3. Amplitude coupling:
   - basis forest: the write sets `nodeScaleFactors_[f]`, re-derives the leaf scale with
     `setForestBasis`'s expression in its order, and re-imposes the map. Kept: s, the divisor 0.674,
     the row norm c and `amplitude.prior.variance`.
   - forest with no basis: the write sets the half-Cauchy scale, both in the combiner and in the echo
     the reader reads. The leaf scale, factor, divisor and the live variance auxiliary stay; the
     auxiliary is refreshed under the new scale after the next sweep's block draw.
   - both: R mirrors the new sd into `params[[f]]`, slot 4 or slot 7 (slot 7 > 0 marks the
     no-basis channel, the engine's own test), so `getPointer`, `setState` re-creation, `copy` and
     save/load build with it.
4. Interactions. `setState`: the installed state's k and leaf scale win, as for every other state
   install; a pre-write state marks a basis forest foreign (its sd reads NA) until `setForestBasis`
   or a write re-imposes the map; the half-Cauchy scale stays as written. `setForestBasis` after a
   write keeps the written factor. No re-anchoring channel exists to restate: every updateScale swap
   and `setData` is refused on both, the multinomial anchor is a constant, and
   [`reissueNamedLeafSd`](../../R/dbarts.R) stays a no-op. Neither admits a variance forest.
   `setCounts` and `setCategoryOffset` leave k alone.
5. Still refused (R-side, before the .Call; bridge backstops keep their wording):
   - multinomial, a named sd: "$setLeafPrior on a multinomial sampler takes normal(k = ) with a
     fixed k, as its creation does: the softmax calibration map sets every category forest's leaf
     scale, so a named 'sd' has nowhere to land"; a `chi()` law or k string: "... a 'k' hyperprior
     is not supported ..., at creation or after"; `forests =`: "its forests are its categories;
     normal(k = ) states every one".
   - amplitude sampler, `leaf.prior`: "the calibration map sets every forest's leaf scale; state a
     forest's spread as at creation, forests = list(forest(sd = ), ...)".
   - amplitude sampler, any other `forest()` knob: "$setLeafPrior changes only a forest's 'sd';
     '<knob>' is fixed at creation" (`basis`: "change it with $setForestBasis"); a list longer than
     the forest count, a non-`forest()` element, an sd not positive and finite; an NA sd gets the
     reader's NA message (dec-A124) ahead of [`validateForestKnobs`](../../R/model.R).
   - `forests =` on a single-forest sampler: refused by name, as [`resolveForests`](../../R/model.R)
     refuses `sd` there. `linear()`/`gp()`: the existing leaf-model refusal.
   The test is the capability [`samplerCarriesAmplitudes`](../../R/bartcore.R), never a forest
   count, so a one-forest sampler with a basis is covered.
6. RNG class, neutral: the constructors, the sweep and the state format are unchanged. The only
   engine refactor, `setForestBasis` calling the new shared leaf-scale helper, keeps its expression
   and its operand order.

## Constraints

- Neutral-class gates, plus the causal-forest and multinomial exact gates the ruling names. The
  equivalence baselines must stay bitwise.
- No flat C entry and no `dbarts.h` change, so there is no ABI event. The multiplier-combiner.md
  bullet saying a flat entry answers 0 is stale: no such entry exists. Step 6 rewrites that bullet.
- Facade virtuals change: every install is `--preclean`.
- Out of scope: writing `amplitude.prior.variance` or `update.amplitude`, a k hyperprior or named sd
  on multinomial, `setModel` on multi-forest samplers, xbart. Creation's `validateForestKnobs`
  accepts `sd = Inf`. That is a pre-existing gap: the writer refuses it in the bridge, and tightening
  creation is a one-line call for the maintainer.

## Steps

1. Engine, [`Chain`](../../src/bartcore/chain.hpp):
   - `setForestFixedK(f, k)`: false when f names no forest, the forest draws its k, or f has a map
     entry (k pinned at 1). It skips the write when `k == forest.k`.
   - `setForestMapSd(f, sd)`: false off a map forest. It skips the write when the value in force is
     the same double; on a basis forest, the map must also still be in force for the skip.
   - `mapLeafScale(f)`: the helper `setForestBasis` now calls.
   - `ForestCombiner::setAmplitudePriorScale(f, scale)`: a virtual, default false;
     `AmplitudeForestCombiner` writes it only on a scale-mixture block.
   - Fan-outs in [`Sampler`](../../src/bartcore/sampler.hpp), as for `setForestPriorScale`: every
     chain, each skipping independently.
2. Facade: two `SamplerBase` virtuals with their impl forwarding
   ([`SamplerBase`](../../src/bartcore/facade.hpp)), and the spy table
   [`FacadeVirtual`](../../tests/cpp/test_facade.cpp).
3. Bridge: `bartcore_setForestK(ptr, k)` for every forest, and `bartcore_setForestSd(ptr, forest,
   sd)`, each refusing a non-finite or non-positive value. They are registered in
   [`R_callMethods`](../../src/R_interface.cpp). The multinomial predicate is the same on every
   forest, so the first forest's refusal comes before any write.
   [`bartcore_setLeafPrior`](../../src/R_interface_bartcore.cpp)'s backstop message names both routes.
4. R: `setLeafPrior` dispatches by capability per Settled 1-5, stores state per `updateState`, and
   updates its docstring. [`reportLeafPrior`](../../R/dbarts.R) changes per the Decision.
5. Tests (below).
6. Docs: man/forest.Rd (a sentence on the writer); multiplier-combiner.md's
   [What this family does not do](../design/multiplier-combiner.md#what-this-family-does-not-do)
   bullet; the existing 1.0-0 NEWS entry for `$setLeafPrior` (setLeafPrior never reached main, so no
   new entry).

## Tests

New inst/tinytest/test-multiforest-leaf-prior-writer.R. Every draw comparison covers all chains and
every channel; warnings are counted.

- Round trip: twins with the same seed. One writes `getLeafPrior()`'s entries back: `normal(k = )`
  on a multinomial, and `forests = ` the per-forest `leaf.prior` list on a K = 2 gaussian causal
  forest, a probit one, and a K = 3 `forests =` sampler. Draws and `getLeafPrior()` stay identical.
- Identity (the dec-B122 pattern): A, created with P, runs, writes P' and stores state. B is created
  with P'. Both install A's state, and both runs are bitwise identical, with `getLeafPrior()`
  identical on every forest. Covered: multinomial k 2 -> 3, the basis forest's sd, the half-Cauchy
  median, and both at once. Installing the state on both sides removes the restore's
  last-ulp difference.
- Discrimination: a changed write moves draws. `getK()` reads the new k on every forest and chain,
  and `anchor` follows a basis forest's write.
- Interactions: `setForestBasis` after a write equals creation with the written sd followed by the
  same `setForestBasis`. A pre-write `setState` reads NA until a write re-imposes the map. `copy()`
  and a saveRDS/readRDS re-creation keep the write. A partial list, or `forest()` with no `sd`, is
  bitwise inert. `setResponse(updateScale = FALSE)` keeps the write.
- Every refusal in Settled 5, by message pattern. Also flip
  ["threeForests$setLeafPrior(normal(k = 2))"](../../inst/tinytest/test-forest-basis-r5.R) and
  ["sampler$setLeafPrior(normal(k = 3))"](../../inst/tinytest/test-multinomial-r5-surface.R), and
  update test-calibration-midchain.R's map-forest reader expectations.
- tests/cpp, beside [`testForestCalibration`](../../tests/cpp/test_sampler.cpp): the chain-level
  identity for both writers and both branches, bitwise in `forestCalibration` and the draws;
  equal-write inertness; each false return.
- Poisons, each reverted and touched, each failing its gate: drop the combiner write (half-Cauchy
  identity), reassociate the leaf-scale expression (bitwise leaf-prior identity), skip the control
  mirror (re-creation test).

## Verification

Against the slice's library (`R_LIBS=$LIB` on every R call), each gate on its own exit status:

```sh
R CMD INSTALL --preclean -l $LIB .
(cd tests/cpp && make && ./test_bartcore)
Rscript -e 'lintr::lint_package()' && air format --check .
Rscript tools/check-rc-codoc.R . && Rscript tools/check-win-drift.R . && Rscript tools/check-doc-freshness.R .
Rscript -e 'tinytest::test_package("dbarts")'
for g in bcf-exact bcf-exact-weak bcf-exact-restricted bcf-latent-exact multinomial-exact; do
  Rscript benchmarks/R/$g.R quick || echo "FAIL $g"; done
```

- Reference build: the three equivalence compares against MANIFEST's current baselines, `--bitwise`;
  count "identical draws (same RNG stream)" lines, 53/53, bcf 15/15, multinomial 11/11, no "max |z|"
  line. The four seeded-drift snapshot files pass unchanged.
- tests/cpp under ASan/UBSan, and the new tinytest file on the R-loaded ASan path (README, Gate
  hygiene).
- `R CMD check --as-cran` on a clean staged tarball; `inst/NEWS.Rd` parses, entry count unchanged.
- stan4bart and bartCause suites against `$LIB` pass (confirms the Decision's claim about bcf).

## Landing

Pending. Records: this note, the INDEX row, the TODO entry closing, dec-B142's record line.
