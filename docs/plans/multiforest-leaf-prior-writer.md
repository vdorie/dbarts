# multiforest-leaf-prior-writer: setLeafPrior on multi-forest samplers

Status: PLANNED 2026-09-29; revised after a blind critique. The reader question was settled by the
orchestrator as an agent-made call (dec-A127): the reader reports forest(sd = ) on map forests, and
creation refuses forest(sd = Inf) as the writer does. Implementation next.

agent: opus (engine, facade, bridge and R in one worktree; one writer)
rng: neutral on every baselined path. One unbaselined path changes, on purpose: a gaussian amplitude
sampler that is re-created after a response or offset swap now keeps its creation anchor (Settled 4).
budget: ~850 lines: engine ~85, facade ~20, bridge ~60, R ~185, man/NEWS/docs ~90, tests ~410 (R ~270,
tests/cpp ~140). Plans have run 1.5-2x low, so expect 1300-1700. That is ~2x the ruling's ~450: the
reader change, the anchor carry and the identity tests make the difference. stan4bart and bartCause 0.

## Goal

On multi-forest samplers, `$setLeafPrior` accepts what their creation accepts (dec-B142):

- a multinomial sampler: `normal(k = )` with a fixed k, applied to every category forest;
- a sampler that carries forest amplitudes: a forest's `sd` stated as at creation,
  `forests = list(forest(sd = ), ...)`.

The single-forest contract holds. A write equal to what is in force is bitwise inert. A write takes
effect on the next sweep, reinterprets no drawn value, and survives re-creation.

## Context

- Today [`setLeafPrior`](../../R/dbarts.R) refuses both samplers through
  [`refuseCountsMutation`](../../R/bartcore.R), [`refuseAmplitudeMutation`](../../R/bartcore.R). The
  engine refuses too: [`Chain::setForestPriorScale`](../../src/bartcore/chain.hpp) returns false
  whenever a combiner exists.
- Contracts: the single-forest one is in the [Landing](leaf-prior-k-or-sd.md#landing) note of
  leaf-prior-k-or-sd.md; the reader's is in [Reader shape](leaf-prior-reader-shape.md#reader-shape).
- Multinomial: [`Chain::buildMultinomialForest`](../../src/bartcore/chain.hpp) builds every category
  forest from [`MultinomialSpec`](../../src/bartcore/combiner.hpp)'s constant anchor and the host k.
  The per-leaf sd, `leaf.scale / k`, is read live each sweep by the leaf draws and by
  [`MultinomialForestCombiner::afterCombine`](../../src/bartcore/combiner.hpp). Nothing caches k.
  Creation accepts `normal(k = Inf)`, and a sampler built with it runs.
- Amplitude coupling: the map is set out in
  [The calibration map, general in K](../design/multiplier-combiner.md#the-calibration-map-general-in-k).
  The engine decides each forest's channel from its amplitude prior, never from its basis:
  - Scale-mixture forest (`amplitudePriorScale > 0`, params slot 7): `sd` is the half-Cauchy median,
    held in [`ForestAmplitudePrior`](../../src/bartcore/combiner.hpp)'s `halfCauchyScale` and echoed
    in `amplitudePriorScales_`. Its node scale stays at the anchor s.
  - Fixed-variance forest: `sd` is `nodeScaleFactor` in `factor * s / (0.674 * c)`.
  - [`forestParams`](../../R/model.R) picks the channel from whether a forest has a basis, but only
    at creation. After [`Chain::setForestBasis`](../../src/bartcore/chain.hpp) gives forest 1 a basis,
    forest 1 is still a scale-mixture forest.
  - `setForestBasis` re-derives a fixed-variance forest's leaf scale from the retained
    `nodeScaleAnchor_`, and re-imposes the map by setting `nodeScaleIsMapDerived_`.
- Anchor s: under gaussian it is the sample sd of the working response at construction
  (`latentScaleAnchor`). Re-creation recomputes it from the current `data@y`. So after
  `setResponse(updateScale = FALSE)`, a copy, a save/load or a `getPointer` re-creation builds on a
  different s (measured 1.4986 against the live 1.9505). The installed state's leaf scales then read
  as foreign, and a later `setForestBasis` re-imposes the map on the wrong s. This is a pre-existing
  defect; the writer would inherit it.
- State: [`ForestStateData`](../../src/bartcore/combiner.hpp) carries each forest's k and leaf scale,
  and a state install overwrites both. It does not carry the half-Cauchy scale; the saved variance is
  the live auxiliary. See [`Chain::noteInstalledLeafScale`](../../src/bartcore/chain.hpp),
  [`Chain::adoptInstalledAmplitudePriors`](../../src/bartcore/chain.hpp).
- Creation reads per-forest values from `attr(control, "bartcore.forests")`
  ([`applyAmplitudeSpec`](../../src/R_interface_bartcore.cpp)).
  [`bartcoreSamplerSetResponse`](../../R/bartcore.R) already mirrors a mutation into a control
  attribute. With `updateState = FALSE` no state is stored, and `copy` builds from control, model and
  data alone.

## Settled

1. Surface: `setLeafPrior(leaf.prior, forests = NULL, updateState = NULL)`, taking exactly one of the
   first two.
   - `forests =` mirrors creation: the same argument name, the same constructor (resolved by
     [`evalInForestVocabulary`](../../R/family.R)) and the same positions.
   - One call restates several forests, and everything is validated before any engine write.
   - `getLeafPrior()`'s per-forest list maps onto it.
   - Names on the list must equal the creation labels; different ones are refused.
   - A short list reaches only the first forests, as at creation. An undeclared `sd` leaves its
     forest as it is.
   - Rejected: a `forest` index argument, which would put a `forest()` in `leaf.prior`. R would
     partially match `forest =` to `forests =`, so the method refuses it by name.
   - `setLeafPrior` never reached main, so the new argument order breaks no released call.
2. Multinomial: nothing in the softmax map is recomputed. The anchor and the leaf scale stay; only
   each forest's k moves.
   - k is any positive value creation accepts, `Inf` included, so the reader's `normal(k = Inf)`
     writes back.
   - The model field records the write (via [`restateLeafPrior`](../../R/dbarts.R)), which is what a
     stateless re-creation reads.
3. Amplitude coupling, keyed on the channel, never on `data@bases`:
   - Fixed-variance forest: the write sets `nodeScaleFactors_[f]`, re-derives the leaf scale with
     `setForestBasis`'s expression in its order, and re-imposes the map. Kept: s, 0.674, c, and
     `amplitude.prior.variance`.
   - Scale-mixture forest: the write sets `halfCauchyScale` in the combiner and in the reader's echo.
     The leaf scale and the live auxiliary stay. The auxiliary is refreshed under the new scale after
     the next sweep's block draw.
   - Both: R mirrors the sd into `params[[f]]`, slot 7 when slot 7 > 0 and slot 4 otherwise. Every
     re-creation builds with it.
   - Accepted as no-op writes, because creation accepts them: `leaf.prior = normal()` and
     `normal(k = 2)`.
4. Anchor carry. At first creation R records s, read off the engine, as `bartcore.forests$anchor`.
   [`applyForestAttributes`](../../src/R_interface_bartcore.cpp) passes it to a new
   `AmplitudeSpec` field (implemented there rather than in `applyAmplitudeSpec`, which never sees
   the attribute list). When that field is finite, the
   constructor uses it instead of `latentScaleAnchor`.
   - It is the same double creation computed, so re-creation is bitwise.
   - A control without the field (a fit saved before this slice) behaves as today.
   - The engine exposes s through a new `ForestCalibration::mapAnchor`, carried as one more internal
     column of [`bartcore_getLeafPrior`](../../src/R_interface_bartcore.cpp). The reader does not
     report it.
5. Reader ([`reportLeafPrior`](../../R/dbarts.R), dec-A127):
   - On a map forest, `leaf.prior` is `forest(sd = )`, taken from `amplitude.prior.scale` on a
     scale-mixture forest and from `leaf.scale.factor` otherwise.
   - `prior.sd.of` there reads `"amplitude scale"` or `"forest total"`.
   - It is an S3 `dbartsForest`, so callers read `$sd`; `@k` and `@prior.sd` fail there.
   - After a foreign state install, or when chains disagree, it reads `forest(sd = NA)`.
   - A scale-mixture forest that `setForestBasis` gave a basis reports its median. That spec
     round-trips through `setLeafPrior`. A fresh `dbarts()` given the same bases would route the
     same `sd` to the fixed-variance channel, so it would build a different prior; the manual says
     so.
6. Interactions:
   - `setState` installs the state's k and leaf scale. A pre-write state marks a fixed-variance
     forest foreign until `setForestBasis` or a write re-imposes the map; the half-Cauchy scale stays
     as written.
   - `setForestBasis` after a write keeps the written factor.
   - `setResponse` and `setOffset` at `updateScale = FALSE` keep the live s, and with Settled 4 so
     does every re-creation. The updateScale swaps and `setData` stay refused on both samplers, and
     [`reissueNamedLeafSd`](../../R/dbarts.R) stays a no-op.
   - Neither sampler admits a variance forest. `setCounts` and `setCategoryOffset` leave k alone.
7. Still refused (R-side, before the .Call; bridge backstops keep their wording):
   - Multinomial:
     - a named sd: "$setLeafPrior on a multinomial sampler takes normal(k = ) with a fixed k, as its
       creation does: the softmax calibration map sets every category forest's leaf scale, so a
       named 'sd' has nowhere to land";
     - a `chi()` law or a k string: "... a 'k' hyperprior is not supported ..., at creation or after";
     - `forests =`: "its forests are its categories; normal(k = ) states every one".
   - Amplitude sampler:
     - any other `leaf.prior`: "the calibration map sets every forest's leaf scale; state a forest's
       spread as at creation, forests = list(forest(sd = ), ...)";
     - any other `forest()` knob: "'<knob>' is fixed at creation" (for `basis`: "change it with
       $setForestBasis");
     - a list too long, mismatched names, or a non-`forest()` element;
     - an sd that is not positive and finite. An NA sd gets the reader's NA message (dec-A124).
   - Creation now refuses `forest(sd = Inf)` too: [`validateForestKnobs`](../../R/model.R) adds
     `is.finite` for `sd`.
   - `forests =` on a single-forest sampler, as [`resolveForests`](../../R/model.R) refuses `sd`
     there.
   - The capability test is [`samplerCarriesAmplitudes`](../../R/bartcore.R), never a forest count.

## Constraints

- Gates: the neutral-class gates, plus the causal-forest and multinomial exact gates. The
  equivalence baselines must stay bitwise; that is the evidence that the constructors, including the
  anchor override when it is absent, moved nothing.
- No flat C entry and no `dbarts.h` change. The multiplier-combiner.md bullet saying a flat entry
  answers 0 is stale; Step 6 rewrites it.
- Facade virtuals and `ForestCalibration` change, so every install is `--preclean`.
- Out of scope:
  - writing `amplitude.prior.variance` or `update.amplitude`;
  - a k hyperprior or a named sd on multinomial;
  - `setModel` on multi-forest samplers;
  - xbart.

## Steps

1. Engine:
   - [`Chain`](../../src/bartcore/chain.hpp):
     - `setForestFixedK(f, k)`: false when f names no forest, the forest draws its k, or f has a map
       entry. It skips the write when `k == forest.k`.
     - `setForestMapSd(f, sd)`: false off a map forest. It skips the write when the value in force
       is the same double; a fixed-variance forest must also still be map-derived.
     - `mapLeafScale(f)`, shared with `setForestBasis`.
   - `ForestCombiner::setAmplitudePriorScale(f, scale)`: a virtual, default false. The amplitude
     combiner writes only a scale-mixture block.
   - The `AmplitudeSpec` anchor override and `ForestCalibration::mapAnchor`.
   - Fan-outs in [`Sampler`](../../src/bartcore/sampler.hpp), as for `setForestPriorScale`.
2. Facade: two [`SamplerBase`](../../src/bartcore/facade.hpp) virtuals, and the spy table
   [`FacadeVirtual`](../../tests/cpp/test_facade.cpp).
3. Bridge:
   - `bartcore_setForestK(ptr, k)`, for every forest, refusing k that is not positive; `Inf` is
     accepted. The multinomial predicate is the same on every forest, so a refusal comes before any
     write.
   - `bartcore_setForestSd(ptr, forest, sd)`, refusing sd that is not positive and finite.
   - `applyForestAttributes` reads `anchor`, and `bartcore_getLeafPrior` gains the `mapAnchor`
     column (named `map.anchor`). `bartcore_setForestK` also gates on the counts capability, since a
     single-forest sampler at a fixed k meets the engine predicate.
   - Registration in [`R_callMethods`](../../src/R_interface.cpp). `bartcore_setLeafPrior`'s
     backstop names both routes.
4. R:
   - `setLeafPrior` per Settled 1-3 and 7.
   - The anchor record in `initialize`, only when the attribute is absent.
   - `reportLeafPrior` per Settled 5.
   - `validateForestKnobs`.
   - The `setLeafPrior` and `getLeafPrior` docstrings; `getLeafPrior`'s currently says
     `normal(sd = )` on map forests.
5. Tests (below).
6. Docs:
   - man/forest.Rd: the writer, and the scale-mixture round-trip caveat.
   - multiplier-combiner.md: its
     [What this family does not do](../design/multiplier-combiner.md#what-this-family-does-not-do)
     bullet, and a sentence on the anchor carry.
   - leaf-prior-reader-shape.md's map-forest bullet.
   - The existing 1.0-0 NEWS entry for `$setLeafPrior`. No new entry: the anchor defect never
     reached main.

## Tests

New file inst/tinytest/test-multiforest-leaf-prior-writer.R. Every draw comparison covers all chains
and every channel, and warnings are counted. Cases: multinomial; a K = 2 gaussian causal forest; a
probit one; a K = 3 `forests =` sampler.

- Twin identity (the main oracle). Creation consumes R's RNG and a write does not, so same-seed twins
  are bitwise identical. A is created with P under seed S and writes P' at once. B is created with P'
  under seed S. Their runs and `getLeafPrior()` must be bitwise identical. Cases:
  - multinomial k 2 -> 3, and -> Inf;
  - a fixed-variance sd;
  - a scale-mixture median;
  - both at once.
- Mid-run half-Cauchy (that scale is not in the state). A is created with P, runs, writes the median
  P', and stores its state. B is created with P'. Both install A's state, and their runs must be
  bitwise identical.
- Round trip: twins, one writing every forest's `leaf.prior` back (`normal(k = )`, or the
  `forests =` list). Draws and the reader are unchanged. A no-op `normal()` is inert too.
- Discrimination: a changed write moves draws, `getK()` reads the new k, and `anchor` follows a
  fixed-variance write.
- Re-creation:
  - gaussian: `setResponse(updateScale = FALSE)`, then `copy()`, a saveRDS/readRDS and a
    `setForestBasis`. Every forest's factor stays non-NA, and `anchor` is bitwise the live
    sampler's. A write afterwards equals the same write on the live sampler.
  - With `updateState = FALSE`: a write followed by `copy()` keeps the k (multinomial), the factor
    and the median.
  - A pre-write `setState` reads `forest(sd = NA)` until a write re-imposes the map.
  - `setForestBasis(1, ...)` after a median write keeps the median and the channel.
- Refusals: every one in Settled 7, by message pattern, and creation's `forest(sd = Inf)`. Also flip
  ["threeForests$setLeafPrior(normal(k = 2))"](../../inst/tinytest/test-forest-basis-r5.R) to a
  no-op, flip ["sampler$setLeafPrior(normal(k = 3))"](../../inst/tinytest/test-multinomial-r5-surface.R)
  to an accepted write, and update test-calibration-midchain.R's map-forest reader expectations.
- tests/cpp, beside [`testForestCalibration`](../../tests/cpp/test_sampler.cpp): the chain-level twin
  identity for both writers and both channels; the anchor override reproducing a construction
  bitwise; equal-write inertness; each false return.
- Poisons. Each is reverted and the file touched, and each must fail the named assertion:
  - dropping the combiner write fails the mid-run half-Cauchy test and the median twin;
  - reassociating the leaf-scale expression fails the fixed-variance twin's bitwise `anchor` (choose
    an sd for which the reassociation moves bits, and show that it does);
  - skipping the control mirror fails the `updateState = FALSE` copy (factor and median);
  - skipping the model-field record fails the same copy's `getK()`;
  - dropping the anchor override fails the gaussian re-creation arm.

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

- On the reference build:
  - the three equivalence compares against MANIFEST's current baselines, with `--bitwise`;
  - count the "identical draws (same RNG stream)" lines: 53/53, bcf 15/15, multinomial 11/11, with
    no "max |z|" line;
  - the four seeded-drift snapshot files pass unchanged.
- The new file's twin-identity and gaussian re-creation arms pass on both the shipped and the
  reference build.
- tests/cpp under ASan/UBSan, and the new tinytest file on the R-loaded ASan path (README, Gate
  hygiene).
- `R CMD check --as-cran` on a clean staged tarball. `inst/NEWS.Rd` parses with its entry count
  unchanged.
- The stan4bart and bartCause suites pass against `$LIB`.

## Landing

Pending. Records: this note, the INDEX row, the TODO entry closing, dec-B142's record line.
