# adoption-slate

Status: IN PROGRESS 2026-08-15 (S1 LANDED 5bedf923, S2 LANDED da3c76f9,
  S3-S8 pending;
  landing notes at EOF). Two blind critique rounds and one
  fix-verification round discharged (5 BLOCKER, 15 MAJOR, 23 MINOR, all
  adopted). Owner fork F1 SETTLED by VD 2026-08-15; F2-F6 and the slice
  adoptions resolved under the delegated grant; see "Decisions". Design evidence is GIT-IGNORED session material
  (the tool-verified census, both
  critiques, the verification pass and the working draft); every load-bearing
  number is carried HERE, since this file is the only tracked record. The
  slate this arc adopts is `docs/design/r-c-division.md`, "Adoption slate".
agent: S1 opus (an engine accessor through the facade, plus a `storeSample`
  refactor). S2 opus (a new reporting channel). S3 opus (a correctness
  adjudication, then a bridge refusal). S4 opus (augmentation math and an RNG
  contract). S5 sonnet (R-only, one diagnostic port). S6 sonnet (docs and
  recipe tests). S7 sonnet (three refusals, no algebra). S8 opus (C API
  design). Serialized: one implementer, each slice lands before the next.
rng: NEUTRAL on every slice. Not one slice moves a draw, so the equivalence
  trio expects BITWISE on all three harnesses at every slice and NO baseline
  re-record is authorized anywhere in this arc; a deviation is a LEAK, never a
  re-record. Three acts carry that risk and are fenced by name: S1's
  `storeSample` refactor (draw-neutral by construction - the function performs
  no draws - under the operation-order constraint in S1); S2's `bart2Negbin`
  driver replacement (OUT OF SCOPE, a different rng class); S3's optional
  grouped-swap equivalence scenario (NOT taken, it would be a new stream).
window: S1 precedes S2 and S3 (the slate's own text schedules G3 and G10 after
  the getLatents slice, and S1 rewrites the `\value` block both extend). S4
  precedes S6 (recipe 1 consumes the helpers). S5 precedes S6 (recipe tests
  call the validator) but is INDEPENDENT of S4 - its oracle is an analytic
  conjugate step built from `rnorm`. S8 is last: it is the arc's only
  irreversible act and its content is downstream of S2 and S4. Order:
  S1, S2, S3, S4, S5, S6, S7, S8.
budget: **RAW ADDITIONS-ONLY (`git show --numstat`), NOT dense-equivalent; the
  1.5x stop applies to these raw numbers.** This is a DELIBERATE ARC-LEVEL
  DEVIATION from `multiplier-combiner.md`'s "Budget units" convention, taken
  because the repo records two conflicting conversions (raw ~ 2x dense as
  written there; dense ~ 0.8x raw as practised, per this arc's own reading of
  `binary-kforest-prior-default.md`'s S1 landing note) and a contested factor
  makes a dense figure unfalsifiable. Every extracted per-slice file, if this
  arc is ever split, repeats this marker verbatim. Divide by 1.25 to 2 for a
  dense reading. Arc: non-test ~2142, tests ~1738, total ~3880 raw additions.
  Per slice, non-test / tests: S1 259/345 (stop 389/518), S2 180/150,
  S3 145/223, S4 336/230, S5 411/200, S6 431/360, S7 135/105, S8 245/125.
  Core S1-S6 alone is 1762/1508. Each figure includes its own records.
  Budgets are sized to the MANDATED ORACLE, not the code delta: S3 is under 50
  lines of bridge code behind a simulation-based calibration arm, and that
  asymmetry is the point.

## Goal

The six open items of the r-c-division adoption slate ship, plus the residue
and the one ABI event they imply. When the arc is done: `getLatents` states
what it returns per family on both surfaces and a caller can read the combined
fit without the offset on any sampler that reports an additive location; the
negative-binomial dispersion is a first-class per-draw channel instead of a
state-block read; a grouped sampler accepts the response swap its model
already delegates correctly and refuses the one mode where its group effects
would silently change meaning; the augmentation draws the engine runs
internally are callable from R against R's own RNG stream; a host can run
simulation-based calibration over its OWN one-sweep step; six named, tested
recipes and an embedding page say so where an author is looking; three carried
refusals land; and one `dbarts.h` re-bake carries every header-touching item
the arc generates.

## Context

- The slate and its frame: `docs/design/r-c-division.md`, "Adoption slate"
  and "Defects found by this arc" items 4-6. All items pre-release; utility
  decides and price sizes a budget only.
- Landed context this arc must not contradict: the empty-leaf veto weight law
  (`docs/plans/archive/empty-leaf-veto.md`, "Landing"), and G1/G5/G11 at 33f6fdc
  (`docs/design/r-c-division.md`, "Defects").
- Gate vocabulary: `docs/plans/README.md`, "RNG classes and their gates" and
  "Implementation protocol"; its ABI review rule governs S8.
- The two blind critiques corrected two claims the census carried, both now
  load-bearing here: `dbarts_apiHash` hashes entry-point SIGNATURES and not
  struct layout, so a `dbarts_results` append is hash-neutral (which is why
  G3's R half costs no ABI event); and the CI path filters are not uniform
  (see "Verification").

## Decisions (settled 2026-08-15)

**F1 = (C).** `$getFitsWithoutOffset()`, argumentless, offset-free, and NO
with-offset twin. VD's principle, as stated: add surface where FUNCTIONALITY
is missing, not where recovery is convenient - the with-offset value is
recovered by adding the offset the caller installed. The accessor's whole
justification is therefore the multi-forest leg: on a Bayesian causal forest
there is no route from R to the combined offset-free fit by any existing
method (`$getForestFits` gives one forest's internal-scale totals, `$predict`
is refused by `refuseUndefinedTestFits`, `run()$train` carries the offset). The
single-forest identity is DOCUMENTED, not built. The surface division that
makes this coherent, and which S1 must state: `extract` generalizes
`fitted`/`predict` for R-convention users on FIT objects, so the R5 sampler
surface has one audience - model writers - and returns engine primitives.

**F2 = (C)** for the dispersion channel: run-result slot plus R5 getter, R side
only, no `dbarts.h` touch, no hash move, no consumer action. The flat half is
DEFERRED, not dropped: it lands in S8.

**F3 = (B).** Relax grouped `setResponse` at `updateScale = FALSE`; refuse at
`TRUE` on a two-way `gaussian || aft` family predicate; add the matching
`setOffset` guard; `setData` stays refused. The engine re-anchor of `b`, `tau`
and the tau prior scale stays a named door. The pre-registered falsifier
stands.

**F4.** The validator rewrites its driver and ports only the diagnostic
(`rankUniformity` and `sbcDiscreteRank`, with `.Random.seed` hygiene).

**F5.** `.Call` bridge to the existing C primitives, drawing R's stream.

**F6.** The embedding page folds into the recipes slice; the six recipes ship
as a vignette and the embedding page as an Rd topic.

**S7 and S8 adopted** as slices rather than riders, per the standing rule that
surface-bearing items with nameable value land pre-release.

**G4 and the host-shell read guards are TICKETED OUT** of this arc; see
"Tickets drafted out of the arc".

## Shipped-surface deltas (all user-visible)

1. `$getFitsWithoutOffset()`, a new `dbartsSampler` method (S1).
2. A fit-surface table and a boundary sentence in
   `man/dbartsSampler-class.Rd` (S1).
3. Corrected per-family `getLatents` semantics in the Rd, the R5 docstring and
   `dbarts.h`'s comment (S1).
4. A `dispersion` element in `run()`'s result list and a `$getDispersion()`
   method, nbinom only (S2).
5. `setResponse` and `setOffset` accepted on a grouped sampler at
   `updateScale = FALSE`, refused at `TRUE` for gaussian, Student-t and aft
   (S3).
6. `dbartsDrawLatents` and `dbartsWorkingResponse`, exported (S4).
7. `dbartsValidateComposition`, exported (S5).
8. A `dbarts-as-a-component` vignette and a `dbarts-embedding` Rd topic (S6).
9. A magnitude cap on nbinom responses; `hostFor` on bart2's ordinal and
   nbinom hosts; the flat `dbarts_sampler_setWeights` guard (S7).
10. `dbarts_results.dispersion`, `dbarts_sampler_dispersion`,
    `dbarts_drawLatents`, `dbarts_workingResponse`, and a re-baked
    `DBARTS_C_API_HASH` (S8).

## Slice decomposition

### S1. What `getLatents` returns, and the fit without the offset

Scope. `getLatents` is documented per family - a LOCATION for probit, ordinal
and aft; a PRECISION for logistic, nbinom and Student-t - on the Rd, the R5
docstring and `dbarts.h`, whose "(gaussian)" parenthetical is wrong at tip
because a gaussian-family sampler built with `resid.dist = student()` reaches
`TResponse`, which overrides `latents()`. `run()$train` CARRIES THE OFFSET, so
the natural `getLatents() - train` reading is biased (measured at 7-10 oracle
SEs; `docs/design/r-c-division.md` defect 5) and correct usage is incremental.
A new accessor returns the combined per-observation location on the response
scale without the offset, on any sampler reporting a single additive location.

The engine work, and why it is engine work. `Chain::combinedFits` is private,
`SamplerBase` declares `fitScale()` but no `fitShift()`, and its only fit
reader is `forestTotalFits`, which memcpys ONE forest's internal-scale totals.
There is no facade route to the combined, response-scaled, offset-free fit.
One new virtual, three signatures, ONE refusal site:

- `chain.hpp`, public, non-const: `bool fitsWithoutOffset(double* out)`,
  opening `if (numReportedLocations() > 1) return false;` then the scale/shift
  write. `Chain::numReportedLocations()` is already public.
- `sampler.hpp`: `bool fitsWithoutOffset(size_t chainNum, double* out)`
  forwarding to `chains_[chainNum]`.
- `facade.hpp`: the `SamplerBase` pure virtual plus the `SamplerFacade`
  override forwarding to `impl_`, on the `forestTotalFits` pattern.

The bridge does NOT re-test the shape; it maps `false` to the named error. A
second bridge-side test would make the engine branch unreachable from R and
from `tests/cpp` alike - untested dead code.

**`storeSample` is refactored onto it.** `fitScale`/`fitShift` composition
already appears three times in `chain.hpp` (`storeSample`'s training write,
`predictFromSavedSample`, `predictFromCurrentTrees`), so a new accessor that
copies the arithmetic is a fourth. Instead `storeSample`'s `L == 1` training
write calls `Chain::fitsWithoutOffset(dst)` and then adds the offset in place;
the multi-location branch is precisely the branch the accessor refuses, so the
shared path is exactly `L == 1`. Binding constraints: preserve the EXACT
operation order (`dst[i] = scale * src[i] + shift` over all n, THEN
`misc_addVectorsInPlace(offset, n, dst)`, never fused), or the slice leaves
its bitwise class; `storeSample` performs no draws, so the refactor is
draw-neutral by construction and the trio proves it; the two predict sites are
NOT folded in - they walk trees into a test-sized buffer rather than reading
`combinedFits`, so two of the three copies remain and this plan says so.

Constness is a CHOICE, not a compulsion: `unique_ptr` does not propagate
constness to its pointee, so a const `Sampler::fitsWithoutOffset` could call a
non-const `Chain` method exactly as the const `Sampler::forestTotalFits` does
today. Non-const is chosen because the read REFILLS the combiner's scratch
buffer, so a const signature would lie about observable state;
`forestTotalFits` can stay const only because it is a memcpy of a resident
vector.

Naming. `$getFitsWithoutOffset()`, no arguments. `$getFits(offset = )` is
rejected: it sits one letter from `$getForestFits`, which returns
INTERNAL-scale per-forest totals, so a misreading swaps both the scale and the
aggregation; and an `offset = TRUE` mode would duplicate `run()$train` under a
second name.

Multinomial, and the two R surfaces it forces. The accessor REFUSES on a
multi-location coupling. That refusal is not writable on the R5 surface: a
multinomial sampler is never a `dbartsSampler` - `bartcoreMultinomialSampler`
(R/bartcore.R) returns a bare `new.env(parent = emptyenv())` carrying `$ptr`,
`$x`, `$K` - while `man/dbarts.Rd` documents a `forests = ` sampler as "an
ordinary `dbartsSampler`". So: (i) S1 adds `bartcoreFitsWithoutOffset` to
`R/bartcore.R` alongside the R5 method, the two-surface shape d809b944 took;
(ii) the multinomial cell runs against that wrapper on a
`bartcoreMultinomialSampler` handle via `getFromNamespace`, the pattern
`benchmarks/R/sbc.R`'s `sbcMakeMultinomial` uses; (iii) the accessor refuses
on `hostFor` through a READ-SIDE sibling of `refuseHostMutation`, called from
the new method only - the one R5 object fronting a multinomial fit is the host
shell, a gaussian sampler whose `numReportedLocations` is 1, so without a
guard the new read answers from the placeholder. The four shipped readers stay
unguarded, deliberately and by ticket (`host-shell-read-guards`).

`$predict(data@x)` serves the multinomial read the refusal declines, but only
at `savedTreeCapacity == 0`: `predictFromSource` sets
`numSamples = capacity > 0 ? capacity : 1`, so under `keepTrees` it replays
saved samples rather than reporting the current state. The caveat rides both
the refusal message and the Rd.

F1's two documentation items. A fit-surface table in
`man/dbartsSampler-class.Rd` - METHOD / SCALE / OFFSET / CURRENT-OR-STORED
over `run()$train`, `$getForestFits`, `$predict`, `$getLatents`,
`$getFitsWithoutOffset` - and one boundary sentence stating that fit-OBJECT
accessors follow R modelling conventions while the R5 sampler methods are the
composition surface. **IMPLEMENTER INSTRUCTION, binding: open `extract.bart`
(R/generics.R) and read what its `type` arms actually do with the offset
before writing that sentence, and check `extract.dbartsSampler`, which takes
`type = "predictors"` and is not a fitted-value accessor at all. Do not assert
the fit-object side's offset behaviour from convention.**

Anchors. `ProbitResponse::latents`, `OrdinalResponse::latents`,
`LogisticResponse::latents`, `AFTResponse::latents`, `TResponse::latents`,
`NBResponse::latents`, `GroupedResponse::latents` (model.hpp; `GaussianResponse`
declares none and the base default returns null). `Chain::combinedFits`,
`Chain::storeSample`, `Chain::predictFromSavedSample`,
`Chain::predictFromCurrentTrees`, `Chain::forestTotalFits`,
`Chain::numReportedLocations` (chain.hpp). `SamplerBase::fitScale`,
`SamplerBase::forestTotalFits`, `SamplerFacade<L, ResidT>` and its by-value
`impl_`, `SamplerShape::numReportedLocations` (facade.hpp).
`Sampler::forestTotalFits` and the `chains_` member (sampler.hpp).
`bartcore_getLatents`, `bartcore_getForestFits`, `refuseBCFTestSurface`,
`predictFromSource` (R_interface_bartcore.cpp). `dbartsSampler`'s
`getLatents`, `getForestFits`, `getSigmas`, `refuseHostMutation` and the
`hostFor` field (R/dbarts.R). `bartcoreForestFits`, `bartcoreGetLatents`,
`bartcorePredict`, `bartcoreMultinomialSampler` (R/bartcore.R).
`man/dbartsSampler-class.Rd`'s `\value` paragraphs for `getLatents` and
`getForestFits`, and its `\item{active}` block. The comment block above
`dbarts_sampler_getLatents` (dbarts.h); comment-only, so no re-bake is owed.

Budget: **RAW ADDITIONS-ONLY; the 1.5x stop applies to these raw numbers.**
non-test ~259 (engine ~48 = chain.hpp ~26 with the refactor, sampler.hpp ~6,
facade.hpp ~13, ~3 net churn; bridge ~31 incl. one `DEF_FUNC` line; R ~36 = R5
~16 + wrapper ~8 + read guard ~12; Rd ~24 = accessor ~10 + table ~10 +
boundary sentence ~4; dbarts.h ~14; NEWS ~10; records ~96), stop 389; tests
~345 (tinytest ~290, `tests/cpp` ~55), stop 518. Comparable, additions-only,
whole commit: d809b944 "Read and rewrite the leaf calibration mid-chain",
non-test 568 (R 122, docs 170, NEWS 16, inst/common 9, man 25, bridge 91,
ENGINE 135 = chain 69 + facade 25 + model 20 + sampler 21), tests 652 - a
read/write PAIR over a five-field struct with state continuation, whose engine
135 is the honest ceiling for S1's ~48.

Oracle: TWO instruments, because the refactor blinds one of them.
(a) THE IDENTITY CELLS GATE PLUMBING, NOT ARITHMETIC.
`$getFitsWithoutOffset()` equals `run(0, 1)$train` minus the installed offset
at the same state - but after the refactor both sides evaluate the SAME
`Chain::fitsWithoutOffset`, so any defect inside it moves both sides and the
identity stays green. What these cells DO gate is the bridge entry, the R5
method, the `R/bartcore.R` wrapper, the per-chain indexing, the return shape,
the refusals, and the offset add - the one step the two sides do not share.
Thirteen cells: gaussian with a NON-UNIT response range and a non-null offset;
probit; logistic; ordinal; aft; Student-t; nbinom; BCF; grouped;
variance-forest; masked; multinomial (the refusal, via the internal wrapper);
the `hostFor` shell (the refusal).
(b) THE ARITHMETIC GATES ARE TWO, NEITHER SHARING THE SITE. A `tests/cpp` cell
in `test_sampler.cpp` pinning `Chain::fitsWithoutOffset` against a hand-built
`fitScale * combined + fitShift` - SINGLE-FOREST necessarily, since
`combinedFits` stays private and the public `forestTotalFits` equals it only
off a combiner, and THIS fixture is the one that must carry a non-unit
response range. Plus an R RECOMBINATION cell on a single-forest and a BCF
fixture, building the expected value through the INDEPENDENT read path -
`$getForestFits` (a memcpy that never touches the accessor), `$getForestAmplitudes`,
and `response.scale`/`response.shift` from `$getCalibration` - and comparing
it to the accessor. Single-forest expected:
`response.scale * getForestFits(1) + response.shift`. BCF expected:
`response.scale * (a * mu + b_z * tau) + response.shift`, the recombination
`refuseUndefinedTestFits`'s own message directs BCF consumers to.

Mutations, with the GATE and the GREEN set named per mutation.
1. Make the accessor add the offset too (the mutation is an added add: after
   the refactor no subtraction exists). GATE: every identity cell carrying a
   non-null offset. A green offset-bearing cell means that fixture is vacuous.
2. Return the internal scale. GATE: the `tests/cpp` cell and the single-forest
   recombination cell. ALL THIRTEEN identity cells stay GREEN.
3. Report forest 0's totals instead of the combiner blend. GATE: the BCF
   recombination cell ALONE - the identity cells move together and the
   `tests/cpp` cell is single-forest, where forest 0's totals ARE the combined
   fit. Outside this slice only the BCF equivalence baseline would notice, so
   that cell is load-bearing and must not be tidied away.
4. Delete the `numReportedLocations() > 1` early return - one site. GATE: the
   multinomial cell, and it is the only one: the buffer is location-major with
   stride n, so this is an IN-BOUNDS read of category 0's channel, a silent
   wrong value no memory tool can see.
4b. Delete the read-side `hostFor` guard. GATE: the shell cell.
5. Reverse the location/precision assignment in the Rd. Caught by NO test;
   the doc half is gated by reading.

Gates: `R CMD INSTALL . --preclean` - MANDATORY, and for the facade reason
(`CLAUDE.local.md`, "Recompile before testing": always after changing virtuals
in `facade.hpp`, or stale objects bus-error); delete the `benchmarks/kernels`
shim binaries before any re-measure. `tests/cpp` from clean, plain AND under
ASAN/UBSAN (new engine code becomes reachable). Full tinytest FAILURES == 0,
no snapshot re-pinned. Equivalence trio BITWISE, no re-record - and the trio
is this slice's leak detector for the refactor specifically. `air format
--check` and `lintr` on R/dbarts.R and R/bartcore.R. `pkgdown::check_pkgdown`
(the Rd changed; no new topic, so no `_pkgdown.yml` entry). Full `R CMD check`
from a clean-copy tarball. rchk on the next scheduled run.

### S2. The nbinom dispersion as a first-class per-draw channel

Scope. The engine reports the dispersion `r` once per kept draw, on the terms
it reports `sigma` and `tau`, and a caller reads the current `r` mid-sweep
without serializing state. `bart2Negbin`'s per-sweep state-read idiom is
SUPPLEMENTED, not replaced - replacing it changes the sweep structure and moves
the `nbinom` equivalence scenario, which is a different rng class and a
different slice.

F2 = (C): no `dbarts.h` touch. The precedent is `bartcore::Results::cutpoints`,
an engine reporting channel with no `dbarts_results` counterpart. The flat half
lands in S8.

Anchors. `NBResponse::dispersion`, `::carriesDispersion`, `::restoreDispersion`
(model.hpp) - the state already exists, nothing new is drawn.
`bartcore::Results` and `Chain::storeSample`'s `results.tau` write (chain.hpp),
the exact shape the dispersion write copies. `bartcore_run`
(R_interface_bartcore.cpp), which is NOT a simple allocate/name/assign triple:
it computes `numResultSlots = 8 + (hasCutpoints ? 1 : 0) + (hasVariance ? 2 :
0) + (hasForestReporting ? 2 : 0)` and four dependent indices
(`varianceTrainSlot`, `varianceTestSlot`, `forestFitsSlot`, `glueSlot`), each
an offset expression over the preceding conditional slots, plus the keyed
`SET_STRING_ELT` names - that arithmetic is the bulk of the bridge budget.
`bart2Negbin` and its `state[[chain]]$dispersion` read (R/bart.R); the
`"dispersion"` entry of the state slot-name table (R_interface_bartcore.cpp).

Budget: **RAW ADDITIONS-ONLY; the 1.5x stop applies to these raw numbers.**
non-test ~180 (engine ~10, bridge ~40, R ~26, Rd ~8, NEWS ~8, records ~88),
band 160-215, stop 270; tests ~150. Comparable, additions-only, whole commit:
0faeb416 "Report the K-forest calibration map", non-test 368 (R 7, docs 169,
NEWS 15, dbarts.h 29, man 18, C_interface 13, R_interface_bartcore 12,
chain.hpp 105), tests 443.

Oracle: the state read it replaces. `run(0, 1)$dispersion` at sweep s equals
`storeState()[[chain]]$dispersion` read immediately after that sweep, for every
chain, under BOTH fixed r (`dispersion > 0`) and grid-estimated r; the R5
getter equals the same value.

Mutations. (1) Write the dispersion one sample late - RED only on the
ESTIMATED-r arm, so that arm is mandatory and must carry a non-vacuity
measurement (assert r takes at least two distinct values across the recorded
draws). (2) Fill the slot with the grid median rather than the draw - RED on
the same arm. (3) Off-by-one in a downstream slot index after the insertion -
RED on the BCF and heteroscedastic run cells, which must therefore be included
even though they carry no dispersion. (4) Drop the null check so a gaussian run
writes an unallocated slot - caught by ASAN/valgrind only.

Gates: `--preclean` install; `tests/cpp` from clean; full tinytest
FAILURES == 0, no snapshot re-pinned; equivalence trio BITWISE, no re-record;
air + lintr on R/bart.R and R/dbarts.R; `pkgdown::check_pkgdown`; full
`R CMD check` from a clean-copy tarball; **local ASAN `tests/cpp` REQUIRED**
(new engine code reachable; mutations 3 and 4 are only visible there); rchk on
the next scheduled run.

### S3. The grouped response swap

Scope and law. On a sampler carrying grouped random intercepts, `setResponse`
is ACCEPTED at `updateScale = FALSE` and REFUSED at `TRUE` when the base family
carries a data-dependent response transform. `setOffset` gains the matching
refusal, closing a live silent decalibration. `setData` stays refused.

The predicate, a static shape fact:
`shape.numGroups > 0 && updateScale == TRUE && (shape.family ==
ResponseFamily::gaussian || shape.family == ResponseFamily::aft)`.
`ResponseFamily` has no `student` member and reports `gaussian` for a
Student-t sampler (the `supportsActiveRows` comment in `SamplerShape` records
exactly this), so Student-t is covered free. The closest house precedent is
`refusePinnedSigmaChange`, which already carries the identical
`family != gaussian && family != aft` two-way shape for the same reason;
`refuseVarianceForestScaleUpdate` is the condition-keying precedent. Groupable
families are gaussian, Student-t, probit, logistic and aft: grouping is applied
only in `createHolder`, which refuses ordinal and nbinom by name, and
`createMultinomialHolder`/`buildMultinomialSampler` never call
`applyGroupAttribute`. So a grouped probit or logistic sampler takes
`updateScale = TRUE` as the documented no-op it already is, and the Rd says so.

The correctness analysis, which is the slice's real work. `GroupedResponse`
holds `groupEffects_`, `tau_` and `priorScale_` on the BASE MODEL'S INTERNAL
SCALE (its class comment says so; the constructor divides by
`base_->sigmaScale()` once, and `Chain::storeSample` de-scales for reporting).
`GroupedResponse::setResponse` is `shiftFits` (building f + b) then
`base_->setResponse` then `rebuildWorking` (z_new minus b), touching neither b
nor tau; group indices are per-observation and fixed, so a same-length swap
keeps every row in its group. `Chain::run` sweeps the trees against
`workingResponse()` FIRST and calls `refreshLatents` after, which redraws the
base latents against f + b, then b conjugately, then tau - an ordinary
blocked-Gibbs order, so carrying b for one tree sweep is what a Gibbs block
does. It breaks in exactly one place: at `updateScale = TRUE` on a re-anchoring
base, `GaussianResponse::setResponse` recomputes `range_` through `rescale()`
and converts exactly `*sigmaInOut` and `sigmaSqPrior_.scale` - so b, tau and
the tau prior scale silently change meaning on the original scale while sigma
does not. That is defect 3's class, and `GroupedResponse::setOffset`'s own
comment already admits it for the offset channel, which is unrefused today.

Consumer reachability: `rbart_vi`'s in-core path returns the live grouped
sampler, assembled by its callers into `$fit`, and
`inst/tinytest/test-rbart-bartcore.R` reaches
`fit$fit[[1L]]$setResponse(y)` to assert today's refusal.

Anchors. `bartcore_setResponse`'s `shape.numGroups > 0` refusal (becomes
predicate-gated), `bartcore_setOffset` (gains the sibling),
`bartcore_setData` (UNCHANGED), `applyGroupAttribute`,
`createMultinomialHolder`, `refusePinnedSigmaChange`
(R_interface_bartcore.cpp). `GroupedResponse::setResponse`, `::setOffset`,
`::shiftFits`, `::rebuildWorking`, `::refreshLatents`,
`::restoreGroupEffects`, `drawGroupEffects`, `GaussianResponse::setResponse`,
`::setOffset`, `::rescale`, `ResponseFamily` (model.hpp).
`SamplerShape::numGroups`, `::family` (facade.hpp). `Chain::getState`,
`::setState` group arms (chain.hpp) and the `"ranef"`/`"tau"` state slot
names. `dbartsSampler$storeState` (which assigns the R5 field and returns
`invisible(NULL)`) and `$setState` (which `inherits(..., "bartcoreState")`
checks) (R/dbarts.R). `inst/tinytest/test-rbart-bartcore.R`'s grouped
`setResponse` refusal cell and its save/load cell. `runSbcGrouped`,
`sbcMakeGroupedFit` (whose control sets `updateState = FALSE`),
`sbcMakeGroupGenerator`, `sbcTauDraw`, `sbcDiscreteRank`, `rankUniformity`,
`sbcBurnSweeps` and the `grouped-gaussian`/`grouped-probit` CLI configs
(benchmarks/R/sbc.R); `.github/workflows/sbc.yaml`'s measured-wall-clock
header.

Budget: **RAW ADDITIONS-ONLY; the 1.5x stop applies to these raw numbers.**
non-test ~145 (bridge ~39, Rd ~8, NEWS ~8, records ~90), band 130-175, stop
220; tests ~223 (sbc.R swap branch ~40, `tests/cpp` ~60, tinytest ~123 of
which ~14 is the shipped-cell replacement). The measured sibling guard,
`refuseVarianceForestScaleUpdate` at 33f6fdc, is 45 raw bridge lines across
three files with ~18 on the test side; the whole 33f6fdc commit is non-test
245 / tests 90, so ~39 is 16% of a five-defect stop-loss slice's non-test cost
for one guard.

Oracle: SBC over the swap, as a FOUR-step replication. `runSbcGrouped` gains a
`swap` argument; at `swap = TRUE` the grouped fit is built ONCE before the
loop and a pristine state captured with `fit$storeState(); state0 <-
fit$state` - two lines, because `$storeState()` returns `invisible(NULL)` and
`sbcMakeGroupedFit` sets `updateState = FALSE` so nothing else populates the
field. Per replication: (1) theta0 from the priors through the separate
generator sampler, exactly as the rebuild arm does; (2) simulate y0 with b0
folded in; (3) RE-INITIALIZE - `fit$setState(state0)`, then
`sampleTreesFromPrior`, `sampleNodeParametersFromPrior` and `setSigma`, the
overdispersed init `sbcReplication` performs, then
`fit$setResponse(y0, updateScale = FALSE)`; (4) run and rank tau, b1, b2,
avg.f and sigma through `sbcDiscreteRank` and `rankUniformity`. The state
install is the only channel that returns b and tau to a y0-independent start -
there is no `$setTau` - and `Chain::getState`/`setState` do carry them, a round
trip `test-rbart-bartcore.R`'s save/load cell already exercises. A constant
start is legitimate: SBC needs independence from y0, not a prior draw, and the
one thing `state0` does not give - a fresh prior draw for tau - is stated in
the arm's comment.

Wall-clock, from a MEASUREMENT not a timeout. `sbc.yaml`'s header records the
maintainer's measured single-thread arms (ordinal 6.6 min, nbinom 46 min, t
2.6 min, multinomial 7.3 min at R = 200; gaussian 27 s at R = 100) and states
the matrix timeouts are ~3x those. Source arm ORDINAL: 36000 + 150 * 30 =
40500 sweeps/rep, 8.1M at R = 200, in 6.6 min = ~20k sweeps/s. The swap arm is
78000 sweeps/rep and 15.6M at R = 200, so ~13 minutes BEFORE the grouped
decorator's overhead, which is a MULTIPLIER to be measured by a pilot at
R = 10 and recorded, not folded in. The arm runs as a ONE-TIME LOCAL GATE and
is NOT added to `sbc.yaml`'s matrix - not on cost (the matrix already carries
nbinom at 46 minutes) but because it is a one-time adjudication whose mutated
companion cannot live in CI at all, and because `grouped-gaussian` and
`grouped-probit` already exist as CLI configs a maintainer runs on demand.

Additional oracles. A `tests/cpp` cell in `test_model.cpp` pinning the
delegation directly: after a same-length swap `groupEffects()` and
`groupTau()` are UNCHANGED and `workingResponse()` equals the base's new z
minus b elementwise. tinytest: the swap runs and MOVES the fit at
`updateScale = FALSE`; `TRUE` errors naming b and tau on grouped gaussian,
Student-t and aft; a grouped PROBIT sampler ACCEPTS `TRUE` as a no-op; the same
pair for `setOffset`; `setData` still refused; and a statistical agreement cell
(created-on-y2 versus created-on-y1-then-swapped, agreeing after a long burn -
never bitwise, the streams differ). **The shipped cell in
`test-rbart-bartcore.R` is REPLACED, not deleted**: `expect_silent` on the
default-`updateScale` swap plus `expect_error(..., updateScale = TRUE, pattern
= "tau")`; its `setData` and save/load neighbours are untouched.

Mutations. (1) Widen the predicate to allow `updateScale = TRUE` on gaussian -
the swap arm's tau and b ranks must go NON-UNIFORM. **Pre-registered
falsifier: if they stay uniform the analysis is wrong and G10 relaxes
unconditionally. Run the unmutated arm and record its verdict BEFORE the
mutation.** (2) Delete `rebuildWorking()` from `GroupedResponse::setResponse`
- the `tests/cpp` cell goes RED, and the mechanism is that `workingResponse_`
is a separate vector `base_->setResponse` never touches, so it is left at
z_old - b_old entirely. (3) Narrow the predicate to gaussian only - the
grouped-aft error cell goes RED. (4) Guard `setResponse` but not `setOffset` -
the `setOffset` error cell goes RED.

Draw-neutrality, on the decisive fact: a grep for `setOffset|setResponse` over
`benchmarks/R/equivalence.R`, `bcf-equivalence.R` and
`multinomial-equivalence.R` returns NOTHING - no harness scenario calls either
mutator on any sampler - so both the relaxation and the new `setOffset`
refusal are unreachable from every recorded scenario. Do NOT add a
grouped-swap equivalence scenario: it would be a new stream forcing a
re-record, and the SBC arm is the stronger gate.

Gates: `--preclean` install; `tests/cpp` from clean, plain AND **under
ASAN/UBSAN - REQUIRED**, because existing engine code becomes REACHABLE FROM R
for the first time (`docs/plans/archive/interaction-constraints.md`'s
new-reachable-code lesson); full tinytest FAILURES == 0; equivalence trio
BITWISE, no re-record; the SBC swap arm run once locally, unmutated then
mutated, both verdicts recorded; air + lintr on benchmarks/R/sbc.R and the
touched tinytest files; `pkgdown::check_pkgdown`; full `R CMD check` from a
clean-copy tarball; rchk on the next scheduled run.

### S4. Exported augmentation helpers

Scope. The augmentation draws the engine runs internally become callable from
R, against R's own stream, returning the same quantity `$getLatents` returns
for that family.

API, two exported names.

```
dbartsDrawLatents(family, fit, y, weights = NULL, offset = NULL, sigma = 1,
                  dispersion = NULL, thresholds = NULL, df = NULL)
dbartsWorkingResponse(family, latent, y, weights = NULL, offset = NULL,
                      dispersion = NULL)
```

`family` is one of "probit", "logistic", "ordinal", "aft", "nbinom",
"student", matched with `match.arg`, no default. **`fit` is the
per-observation location WITHOUT the offset** - exactly
`$getFitsWithoutOffset()`, or f(x) however the host computes it - and the
helper forms the linear predictor internally as `fit + offset`, mirroring
`ProbitResponse::refreshLatents`'s `totalFits[i] + offset_[i]` and
`LogisticResponse::refreshLatents`'s `psi = totalFits[i] + offset`;
`offset = NULL` means zero. The engine's `totalFits` is offset-free at every
call site, so this is the engine's own convention. `y` is validated with the
SAME rule `bartcore_bridge::validateResponseSupport` applies, which needs a
stated string-to-enum map: `ResponseFamily` has no `student` member, so
"student" maps to `gaussian`, whose arm constrains nothing, and "aft"
likewise. `weights` is LOGISTIC ONLY, entering as `lround(weights[i])`
independent PG(1, psi) draws summed - the identical reduction
`LogisticResponse::refreshLatents` performs - and refused by name for every
other family. `dispersion` is NBINOM ONLY, the current r, entering as
`PG(y[i] + r, psi)` through the integer-shape reduction
`bartcore::simulatePolyaGammaShape` provides and `NBResponse::drawOmega` uses.
`thresholds` is ordinal only; `sigma` aft and student; `df` student.
`dbartsDrawLatents` returns a numeric of length n with
`attr(result, "quantity")` of "location" or "precision".

Two functions and not one, because for logistic and nbinom the drawn latent is
a PRECISION and the quantity a host regresses on is the working response
(`w (y - 1/2) / omega - o`; `0.5 (y - r) / omega - o`, both read off the
engine's refresh sites). Every measured composition defect in the origin survey
was a host getting exactly that step wrong.

RNG contract. The helpers draw from R's stream via an `ext_rng` created with
`ext_rng_createDefault(TRUE)` - which installs `unif_rand` as
`EXT_RNG_ALGORITHM_USER_UNIFORM` and `norm_rand` as
`EXT_RNG_STANDARD_NORMAL_USER_NORM`, consuming no draw at creation - inside a
`GetRNGstate`/`PutRNGstate` pair. `set.seed()` reproduces them and they compose
with any other R draw. They never touch a sampler's chain RNG, which
`createChainRngs` seeds separately and which "never advances R's stream". Both
halves belong in the Rd.

First-caller obligations. `ext_rng_createDefault` has no shipped caller; its
one existing call site is `tests/cpp/rshim.cpp`, whose pattern is the model:
`GetRNGstate(); ext_rng* rng = ext_rng_createDefault(true); if (rng == nullptr)
{ PutRNGstate(); Rf_error(...); }`. S4 owes that NULL path, `ext_rng_destroy`
on every exit INCLUDING `Rf_error`'s longjmp (so the generator is owned by an
`unwindProtect` body or a scope guard, not a bare local), a PER-CALL rather
than cached generator (native mode holds no state of its own), and one Rd
sentence: `random.h` records native mode as "generally not thread safe", so
the helpers must not be called from a worker thread.

OUT OF SCOPE, explicitly: the in-sampler draw sites are NOT refactored to
share code with the helpers. Refactoring them is a draw-law risk and a
different rng class.

Rd and pkgdown: ONE topic, `man/dbartsAugmentation.Rd`, with an `\alias` for
each function, and a single `- dbartsAugmentation` entry under `_pkgdown.yml`'s
existing "The mutable sampler" section, whose `desc` already reads "for use
inside larger MCMC schemes".

Budget: **RAW ADDITIONS-ONLY; the 1.5x stop applies to these raw numbers.**
non-test ~336 (bridge ~110 with two `DEF_FUNC` lines, R ~80, Rd ~45, NEWS ~8,
`_pkgdown.yml` ~3, records ~90), band 300-390, stop 505; tests ~230. Against
33f6fdc's measured bridge 149 for two guard functions with five-way family
branching, ~110 is the honest band. The guide memo's ~240 was the whole item
including tests; the honest total is ~566 raw.

Oracle: three. (1) DISTRIBUTIONAL - truncated-normal mean and variance against
closed forms at several truncation points and both tails; `E[PG(1, psi)] =
tanh(psi/2) / (2 psi)` and `E[PG(b, psi)] = b` times that, including the
`psi -> 0` limit 1/4; each a z on the mean with a stated Monte Carlo error at a
fixed seed. (2) AGREEMENT - a pure-R Albert-Chib probit loop built ONLY from
the two helpers reproduces the engine's probit posterior on a fixed DGP to a
stated tolerance. (3) OFFSET DISCRIMINATION - `dbartsDrawLatents(fit = f,
offset = o)` and `dbartsDrawLatents(fit = f + o, offset = NULL)` agree in the
LATENT and DIFFER in `dbartsWorkingResponse`, so a caller who confuses the
conventions is told; the Rd carries the matching sentence.

Mutations. (1) Truncate on the wrong side - the moment test goes RED and the
Albert-Chib slope inverts. (2) Ignore the PG shape - RED only where weights
exceed 1, so a weighted-logistic cell is MANDATORY, and RED on nbinom only if
y + r > 1, so that fixture must not use y = 0, r = 1. (3) Drop
`GetRNGstate`/`PutRNGstate` - caught ONLY by the reproducibility test, which
must assert BOTH halves (two `set.seed(k)` calls agree, AND the stream
ADVANCES). (4) Make `dbartsWorkingResponse` the identity for logistic - RED on
the agreement arm only, which is why that arm exists.

Gates: `--preclean` install; `tests/cpp` from clean; full tinytest
FAILURES == 0; equivalence trio BITWISE, no re-record; air + lintr; NEW Rd
topic so a `_pkgdown.yml` entry is REQUIRED plus `pkgdown::check_pkgdown`;
full `R CMD check` from a clean-copy tarball; **local ASAN `tests/cpp` and a
valgrind pass over the new tinytest file REQUIRED** (new C code with
caller-sized buffers and a heap-owned rng becomes reachable); rchk on the next
scheduled run.

### S5. The composition validator

Scope. `dbartsValidateComposition` runs simulation-based calibration over a
HOST-SUPPLIED one-sweep step and reports, per scalar functional, whether the
host's composed kernel targets the posterior it claims. It catches the measured
posterior-mean-as-draws defect class; it is not a calibration fix.

API.

```
dbartsValidateComposition(drawPrior, simulate, init, step, functionals,
                          n.replications = 200L, n.draws = 200L,
                          n.thin = 30L, n.burn = 200L, alpha = 0.05,
                          seed = NULL)
```

Five properties an implementer must not guess at. (1) It ranks
`functionals(init(theta0, data))` - the state AT theta0 - among
`functionals(state_1..L)`, NOT a name-intersection between `drawPrior`'s
return and `functionals`' names; SBC's interesting functionals are DERIVED
(sbc.R's own headline is `avg.f = mean(f0Train)`), so `drawPrior` needs no name
contract at all and a mistyped functional name becomes impossible rather than
silently reported. (2) `init` MUST start the chain at theta0 and the Rd states
it as a requirement; `n.burn` defaults to 200 rather than 0 so a host whose
`init` is approximate is not flagged for its early draws, and the validator
errors if `functionals(init(...))` disagrees in length or names with
`functionals(state)` after one step. (3) `sbcDiscreteRank` ports and is used
for every functional - the tie artifact it prevents is a rank-computation
concern, and a host-supplied functional is MORE likely to tie than a dbarts
one. (4) The port SAVES AND RESTORES `.Random.seed`: `rankUniformity` calls
`set.seed(seed)` internally and never restores it, which is fine in a dev
harness and not in an exported function under "safe over fast in R". (5) The
Bonferroni rule is an ADAPTATION, not an inheritance: `sbcMatrixAlpha` corrects
over the SBC matrix's 30 functionals across five family arms, i.e. a whole CI
run, while the validator corrects over `length(functionals(state))` within one
call. Returns an object of class `dbartsCompositionValidation` carrying the
rank matrix, L, per-functional verdicts and the alpha used, with a `print`
method.

What ports and what does not. `rankUniformity` and `sbcDiscreteRank` port -
self-contained, no dbarts dependency, no prior-draw assumption. Everything
upstream does not: `sbcConfig`, `sbcMakeSampler`, `sbcReplication` and every
per-family generator draw theta0 from the SAMPLER'S OWN prior via
`sampleTreesFromPrior`, which is the self-consistency sbc.R's header calls
"the whole game". The validator inverts that; the host owns the prior. That is
why the recorded ~350-line budget, which assumed reuse, is wrong.

The diagnostic is the ecdf-difference statistic against a simulation-based
simultaneous band (Talts et al. 2018 fig. 1), already corrected for multiple
looks across the rank grid, with chi-square and a jittered KS as secondary -
the judgement `rankUniformity` already records and argues.

Blind spot, stated in the Rd: the validator cannot detect a host whose
`drawPrior` and `simulate` disagree with each other.

Budget: **RAW ADDITIONS-ONLY; the 1.5x stop applies to these raw numbers.**
non-test ~411 (R ~240 = driver ~90 + `rankUniformity` port with the seed guard
~75 + `sbcDiscreteRank` ~15 + print/format ~40 + validation ~20; Rd ~70;
NEWS ~8; `_pkgdown.yml` ~3; records ~90), band 370-470, stop 615; tests ~200.
TODO records ~350 for the whole item; the honest total is ~611 raw and the
stop applies to 411.

Oracle: a discriminating pair on an ANALYTIC step - a conjugate normal-normal
one-sweep step whose posterior is closed form, so no dbarts sampler runs in
tinytest. (1) KNOWN-GOOD, the exact conjugate draw: every functional PASSES at
a fixed seed. (2) KNOWN-BAD, the posterior MEAN substituted for the draw: must
FLAG. (3) OVER-DISPERSED, drawing at 2x the posterior sd: must FLAG, proving
the diagnostic catches both directions. A dbarts-driven arm goes in
`benchmarks/`, not tinytest.

Test-cell parameters, with the two justifications assigned to their right
owners. The cells pass `n.thin = 1L` because the conjugate step draws from the
exact posterior so consecutive states are INDEPENDENT, and `n.burn = 0L`
because the cells' `init` is exact at theta0 so the chain starts stationary.
(The API's own `n.burn = 200L` default serves hosts whose `init` is only
approximate.) That drops each cell to 200 x 200 = 40k `step` calls. The band
simulation is `max(nSim, ceiling(20 / alpha))` with `nSim` defaulting to 2000:
at the three discriminating cells' scale `alpha = 0.05/3` gives 1200 < 2000, so
nothing is pushed past the default - it bites only from M >= 6, i.e. at the
many-functional cell. Targets: under 5 s for each of the six ordinary cells
(under 30 s), under 15 s for the many-functional cell, **file under 45 s**,
measured and recorded; a cell that misses drops to `n.replications = 100L`
before anything else changes.

Mutations. (1) Rank with `<=` instead of `<` - RED only if the battery asserts
the mean rank against L/2, so that assertion is mandatory. (2) Drop the
Bonferroni - RED only on the many-functional cell. (3) Take the band from the
wrong quantile tail - the KNOWN-BAD cell goes GREEN, which is why that cell is
the whole gate on the diagnostic's direction. (4) Rank
`functionals(theta0-as-named-vector)` instead of `functionals(init(...))` - RED
only on a cell whose functional is DERIVED, so at least one must be. (5)
Remove the `.Random.seed` restore - RED on the dedicated seed cell.

Gates: install before tinytest; `tests/cpp` from clean; full tinytest
FAILURES == 0; equivalence trio BITWISE, no re-record, engine binary IDENTICAL
(the R-only precedent recorded in `benchmarks/baselines/MANIFEST`); air +
lintr; NEW Rd topic so a `_pkgdown.yml` entry is REQUIRED - a single
`- dbartsValidateComposition` under the existing **Diagnostics** section - plus
`pkgdown::check_pkgdown`; full `R CMD check` from a clean-copy tarball. No ASAN
leg owed.

### S6. Named recipes and the embedding page

Scope. Six named, tested compositions ship as a vignette
(`dbarts-as-a-component.Rmd`), cross-referenced from `?bart` and
`?dbartsSampler`, plus a `dbarts-embedding` Rd topic. The claim defended is
`docs/design/r-c-division.md`'s own: the capability gap is smaller than the
perception gap, and recipes and documentation are the competitive answer.

Recipes, each with a tinytest cell and, where the recipe has a one-sweep step,
a `dbartsValidateComposition` call: (1) AUGMENTATION - pure-R Albert-Chib
probit over a gaussian sampler using S4's helpers, tested against the engine's
probit posterior; (2) OFFSET BLOCK - an outer block writing `$setOffset`,
tested on both blocks' truth, stating the seeding contract; (3) K-FOREST - K
gaussian samplers with the sqrt(K) leaf-prior rescale through
`$setCalibration`, tested against a single-sampler fit; (4) LATENT COVARIATE -
`updatePredictorPerObservationJointly`, tested on the install mask as an MH
accept mask; (5) OUTER-OWNED SIGMA - `resid.prior = fixed()` plus `$setSigma`,
tested on the pin and the guard message; (6) INCREMENTAL USE -
`$getFitsWithoutOffset()` and `getLatents`, showing why
`getLatents() - run()$train` is wrong, tested on the slope bias. Recipe 6
closes defect 5's measured trap and is why S1 builds the accessor.

Carrier split, principled. The recipes are runnable narrative with output, so
they are a VIGNETTE, which `R CMD check` builds - an erroring recipe fails the
check outright. The embedding page is reference material wanting `\link{}`s and
a `?` entry point, so it is an Rd TOPIC, which REQUIRES a `_pkgdown.yml`
`reference:` entry - a single `- dbarts-embedding` under the existing
**Utilities** section, the only section not about fitting, diagnosing or
driving a sampler. `_pkgdown.yml` has no `articles:` key, so the vignette needs
no yaml edit.

Why each recipe ALSO needs a tinytest cell: **a vignette asserts nothing.** A
recipe that runs to completion and produces a WRONG answer - the defect class
this arc exists to prevent - builds cleanly and ships.

Per-recipe cost caps, because both shipped vignettes already run real samplers
under `VignetteBuilder: knitr, rmarkdown`: `n.samples = 100L`,
`n.burn = 100L`, `n.trees = 25L`, n = 200, K = 2, with a **60-second vignette
build target**, measured and recorded.

Budget: **RAW ADDITIONS-ONLY; the 1.5x stop applies to these raw numbers.**
non-test ~431 (vignette ~210, embedding Rd ~110, `\seealso` ~8,
`_pkgdown.yml` ~5, NEWS ~8, records ~90), band 390-490, stop 645; tests ~360,
about 60 per recipe. TODO records ~150 plus ~80, which priced
cross-references rather than six runnable recipes with tests.

Mutations. (1) Drop the sqrt(K) rescale from the K-forest recipe - its test
must go RED; **run this mutation as a non-vacuity measurement BEFORE accepting
the tolerance**, since the lever it defends was measured at 16x. (2) Use
`run()$train` instead of `$getFitsWithoutOffset()` in recipe 6 - the slope-bias
test goes RED at the measured 7-10 oracle SEs. (3) Substitute a posterior mean
for a draw in recipe 2's outer block - the recipe's own validator call must
FLAG, which is the end-to-end proof that S5 works on a real composition. (4)
Break a `\link{}` - caught by `R CMD check`, not tinytest.

Gates: install; `tests/cpp` from clean; full tinytest FAILURES == 0;
equivalence trio BITWISE, no re-record; air + lintr on every touched R file and
Rmd chunk; NEW Rd topic so a `_pkgdown.yml` entry is REQUIRED plus
`pkgdown::check_pkgdown`; full `R CMD check` from a clean-copy tarball WITH
vignettes built, and the build time recorded. No ASAN leg owed.

### S7. Residue

Scope. Three fixes carried unattached through the design and adopted at the
walkthrough, landed together because they share a gate battery and nothing
else. (1) An nbinom response MAGNITUDE CAP, refusing y above a stated bound at
creation and on both mutation conduits, on both surfaces - the deliberate
deviation left open at the G1 landing, whose mechanism is
`NBDispersionPrior::computeKernel` sizing a count histogram from
`lround(max y)`, so y = 1e9 hangs on unbounded allocation. (2) `hostFor` set on
BOTH bart2's ordinal and nbinom hosts, so their placeholder `$fit` refuses
mutation as the multinomial host already does; today the field is assigned at
exactly two sites, both multinomial. (3) `refuseBinaryWeightChange` promoted
out of the anonymous namespace of `R_interface_bartcore.cpp` into
`bartcore_bridge`, declared in `R_interface_bartcore_common.hpp` with the
Doxygen its siblings carry, and called from `dbarts_sampler_setWeights`, which
carries only `refuseMultiForestResponseMutation` today - so probit, ordinal,
aft and nbinom silently ignore a flat weight change and logistic reaches a
division by zero on a negative weight.

Anchors. `NBDispersionPrior::computeKernel` (model.hpp);
`bartcore_bridge::validateResponseSupport` (declared in
`R_interface_bartcore_common.hpp`), whose nbinom arm gains the bound;
`hostFor`, `refuseHostMutation` (R/dbarts.R) and the two multinomial
assignment sites (R/bart.R); `refuseBinaryWeightChange`
(R_interface_bartcore.cpp) and `dbarts_sampler_setWeights` (C_interface.cpp).
The promotion pattern is the one `validateResponseSupport` and
`refuseVarianceForestScaleUpdate` already follow; the measured sibling cost is
11 header lines plus a `using` plus the call site.

Budget: **RAW ADDITIONS-ONLY; the 1.5x stop applies to these raw numbers.**
non-test ~135 (cap 25, `hostFor` 14 for BOTH halves at ~7 each, `setWeights`
promotion 26, records ~70 shared), tests ~105 (35 + 30 + 40).

Mutations. Cap: raise the bound above the fixture's max y - the creation and
both mutation cells go GREEN, so the fixture must sit just over the bound.
`hostFor`: drop either assignment - that host's refusal cell goes RED, which
is why both halves need their own cell. `setWeights`: drop the flat call site
- the capi cell goes RED while the R-side cell stays GREEN, which is the whole
point of the promotion.

Gates: `--preclean` install; `tests/cpp` from clean; full tinytest
FAILURES == 0; equivalence trio BITWISE, no re-record (all three are refusals
on paths no harness scenario reaches); air + lintr on R/bart.R and R/dbarts.R;
`pkgdown::check_pkgdown`; full `R CMD check` from a clean-copy tarball; rchk on
the next scheduled run. No ASAN leg owed - no new engine code becomes
reachable, and the cap makes an existing path unreachable.

### S8. The one header re-bake

Scope. The arc spends exactly ONE `dbarts.h` re-bake and it carries everything
needing one, so no second window is owed before the freeze. Four items:

1. `dbarts_results` gains `double* dispersion`, appended at the marked 1.0-0
   field boundary, with the `offsetof`/`sizeof` lock in `src/C_interface.cpp`
   updated and one `FILL(dispersion, dispersion)` line added. Hash-neutral by
   itself, and `structSize` bounds the write so a consumer built against the
   smaller header is never written past.
2. `dbarts_sampler_dispersion`, the flat mid-sweep getter - the C-side analogue
   of S2's R5 getter, for a host composing a count block in C++.
3. F5's LinkingTo exposure of the augmentation surface.
4. The decision-8 carve-out comment. `dbarts.h` says in two places - the
   `dbarts_results` doc block and the "1.0-0 field boundary" comment - that an
   append bumps `DBARTS_C_API_MINOR`, while binding decision 8
   (`docs/plans/archive/dbarts-h-reshape.md`) pins MINOR at 0 through the window. Item
   1 makes the two collide, so this slice lands the one-sentence carve-out
   ("pre-1.0-0 appends extend the initial field set and move no version
   constant") in the same commit; without it the next author reads the header
   and bumps MINOR against the decision.

**Which entry points cross, decided: the WRAPPED forms `dbarts_drawLatents`
and `dbarts_workingResponse`, mirroring S4's two R helpers and drawing from
R's stream internally - NOT the raw `ext_rng` primitives.** Reasons in order:
exposing `ext_rng_simulatePolyaGamma` and the three truncated-normal entries
would make `ext_rng` itself an ABI type plus create/destroy/seed entries, a
five-to-six-entry surface with its own lifetime contract frozen at 1.0-0, to
expose four functions; every `LinkingTo` consumer is an R package, so R's RNG
is the stream it already owns and `GetRNGstate`/`PutRNGstate` the discipline it
already follows, and the wrapped forms need no rng handle at all; and the
wrapped forms carry the family dispatch, the support validation and the
working-response step, which is the half every measured composition defect got
wrong, so exposing the primitives alone would export the sharp edge and keep
the guard rail in R. The raw primitives stay internal, recorded as a door.

Anchors. `dbarts_results`, `DBARTS_RESULTS_HAS`, the "1.0-0 field boundary"
comment, `DBARTS_C_API_LIST`, `DBARTS_C_API_HASH`, `DBARTS_C_API_MINOR`
(dbarts.h). The `offsetof`/`sizeof` static_assert block, the `FILL` macro, and
`static_assert(dbarts_fnv1a(DBARTS_C_API_DECLS) == DBARTS_C_API_HASH)`
(C_interface.cpp) - **that last assert, not the offsetof block, is what the
re-bake acknowledges.** `bartcore::Results::dispersion` (added by S2) and
`NBResponse::dispersion` (model.hpp). The S4 augmentation cores, promoted from
the bridge's anonymous namespace into `bartcore_bridge` so both surfaces call
one implementation. `inst/tinytest/capi/consumer.c`, the only consumer-visible
hash check in the repo.

Budget: **RAW ADDITIONS-ONLY; the 1.5x stop applies to these raw numbers.**
non-test ~245 (dbarts.h ~54: the field and its doc, three X-list entries,
three prototype/Doxygen blocks, the carve-out - the re-baked literal is a
change, not an addition; `C_interface.cpp` ~57: the offsetof line, the `FILL`
line, three implementations, `using` lines;
`R_interface_bartcore_common.hpp` ~26 for the two promoted declarations;
`R_interface_bartcore.cpp` ~6 net for the move; NEWS ~10; records ~90), band
220-310, stop 368; tests ~125 (`capi/consumer.c` ~70 for four new legs,
`test-capi.R` ~55; the hash handshake needs no new cell, since `consumer.c`
already checks `dbarts_apiHash()` against the header constant, which is
exactly the assertion the re-bake must keep true). Comparable, additions-only:
ab3aa2fa (the reshape S1 re-bake) - dbarts.h 354, C_interface 455,
consumer.c 455, test-capi 644 - which RE-SIGNED four predictor entries onto a
new POD; S8 adds three entries and one field and re-signs nothing, so it is a
small fraction of that, and the comparable confirms the shape: a re-bake's
weight sits in the capi tests.

Lockstep cost, enumerated. An APPEND-ONLY X-list entry changes ZERO call sites
in any consumer - nothing existing is re-signed - so no consumer source
changes. What is owed: re-bake the `DBARTS_C_API_HASH` literal; rebuild
stan4bart (bartcore), treatSens (dbarts-1.0), bartCause (dbarts-1.0) and
bairrtt; and note that any consumer hard-checking
`dbarts_apiHash() == DBARTS_C_API_HASH` fails to LOAD until rebuilt, the
handshake working as designed. That rebuild is already committed: the release
checklist records "re-verify all four sister packages against the final
re-baked header once, at the freeze", so landing S8 before that verification
makes its marginal consumer cost zero. For contrast, the reshape's SIGNATURE
change cost five stan4bart call sites; S8 costs none.

Gates: `--preclean` install; `tests/cpp` from clean; full tinytest
FAILURES == 0; equivalence trio BITWISE, no re-record (three read entries and
one output slot; the augmentation entries call the same cores S4 ships
alongside, never inside, the draw sites); full `R CMD check` from a clean-copy
tarball; **the capi consumer leg is this slice's own gate and must be RUN, not
assumed**; rchk on the next scheduled run; and the ABI review step
`docs/plans/README.md` mandates - "any diff under inst/include/dbarts/ is an
ABI event ... build+test stan4bart against the change before landing" - which
is a HARD GATE on S8, not a follow-up.

## Stale-doc corrections

1. `man/dbartsSampler-class.Rd`, the `\item{active}` block. The sentence
   listing count-based accounting names "the empty-leaf veto", which is wrong
   at tip: the mask composes into the working weights, those become
   `ctx.weights`, and `Tree::leafHasNoWeight` therefore vetoes a leaf all of
   whose members are inactive. The correction is CONDITIONAL, not a flat
   reversal - `leafHasNoWeight` degenerates to the member count bit-for-bit
   when no weight vector is installed. Lands in S1, which already opens this
   Rd; it is a shipped-doc defect, not a plan-doc nit.
2. `docs/plans/latent-subset-mask.md`, rule 2 of "Semantics of inactive: the
   seven uniform rules", which still says "the veto sees it" although the same
   file's S1 landing note records the pin as dropped. **The correction is a
   LINE-COUNT-PRESERVING in-place rewrite plus an EOF errata line**, because
   seven documents cite line numbers into that file and any drift falsifies
   the citations below the edit. Rule 2 is five lines, widths 80/72/73/62/35;
   the replacement is five lines, widths 62/65/71/72/69, all within 80:

   ```
   2. The row STILL OCCUPIES its leaf for COUNT-based accounting:
      `numObservations()` counts it, the scan's `count` sees it, and
      `collapseEmptyNodes` TRIGGERS on member count regardless of weights.
      The empty-leaf VETO does NOT: since empty-leaf-veto-fix (21fc29c3) it
      counts POSITIVE-WEIGHT members, so an all-inactive leaf IS vetoed.
   ```

   Commit that block verbatim; do not redraft it. Lands in S3.
3. `docs/design/r-c-division.md`'s latent-subset-mask slate bullet, written in
   the present tense from inside its own landing commit ("SHIPPED, S0-S4 ...
   this slice"), where "this slice" has no referent for a later reader.
   Replace with the landing commit and drop the present tense. Lands in S1,
   which already opens that file.

## Doors held open (recorded, not scheduled)

- The `GroupedResponse` re-anchor: convert b, tau and the tau prior scale by
  the `sigmaScale` ratio inside `setResponse`/`setOffset`, making
  `updateScale = TRUE` legal rather than refused (F3 alternative (C), ~10
  engine additions). S3's refusal is exactly the mode this would open, so it
  can be taken later without breaking anything S3 ships.
- `fitShift()` on `SamplerBase`, absent while `fitScale()` is present - a real
  asymmetry left by S1's one-virtual choice.
- Promoting `Chain::combinedFits` to public, which would let S1's `tests/cpp`
  cell cover BCF; not taken, because the R recombination cell buys the same
  coverage for ~15 lines without widening the engine surface.
- The raw `ext_rng` primitives on `dbarts.h` (S8 exposes wrapped forms only).
- A grouped-swap equivalence scenario (a new stream; the SBC arm is stronger).

## Tickets drafted out of the arc

Two items are DECIDED-AND-TICKETED, priced into no slice: TODO
`latent-family-weight-channel` (census G4, adjudicated as a PARTIAL BUILD -
logistic built, the other four declined by identification, with the decline
recorded in `docs/design/r-c-division.md`) and TODO `host-shell-read-guards`
(whether the four shipped readers that answer from a `hostFor` placeholder
shell should refuse as the new accessor does). S1's Rd cross-references the
second by name so its asymmetry reads as a decision.

## Verification

Per slice, the gate battery in that slice's section. Arc-wide:

```sh
R CMD INSTALL . --preclean            # every slice; MANDATORY on S1 (facade virtuals)
cd tests/cpp && make && ./test_bartcore
```
```r
tinytest::test_package("dbarts")      # FAILURES == 0, no snapshot re-pinned
```
```sh
Rscript benchmarks/R/equivalence.R compare \
  benchmarks/baselines/equivalence-8b047f8b.rds
Rscript benchmarks/R/bcf-equivalence.R compare \
  benchmarks/baselines/bcf-equivalence-8b047f8b.rds
Rscript benchmarks/R/multinomial-equivalence.R compare \
  benchmarks/baselines/multinomial-equivalence-1027be5.rds
```
Expected: IDENTICAL on every scenario of every harness, every slice. **No
baseline re-record is authorized anywhere in this arc**; a deviation is a leak
and the slice ABORTS.

CI, tool-verified and NOT uniform - six workflows fire on push to `bartcore`
(check-standard, cpp-tests, exact-gates, lint, pkgdown, sanitizers; the other
five are schedule/dispatch only). `exact-gates` has no `benchmarks/**`
path-ignore and `lint` excludes only `benchmarks/baselines/**`, so **S3's
`benchmarks/R/sbc.R` edit fires both**; `pkgdown` has no `**.md` ignore because
README is an input. A `docs/**`-only or `TODO`-only push fires nothing, so the
two plan-doc corrections above fire nothing. Any `inst/include/`, `src/`,
`R/`, `man/`, `vignettes/`, `_pkgdown.yml`, `NAMESPACE` or `DESCRIPTION` touch
fires all six.

## Critique adjudication

Three adversarial passes were run against this plan's working draft, all
findings adopted: critique-1 (2 BLOCKER, 9 MAJOR, 12 MINOR), critique-2
(3 BLOCKER, 6 MAJOR, 11 MINOR) and a fix-verification pass (20/20 findings
landed, 2 residual items). The substantive corrections that survive into this
spec, recorded so a reader knows which claims were adversarially tested: S1 is
an ENGINE slice, not a bridge slice (no facade route to the combined fit
exists); the `storeSample` refactor makes the R identity oracle INSENSITIVE to
defects inside the accessor, so the arithmetic gates are the `tests/cpp` cell
and the recombination cells; the SBC swap arm needs the per-replication
re-initialization `sbcReplication` performs or its falsifier cannot
discriminate; the grouped refusal is family-keyed, not blanket; S4's `fit`
argument is OFFSET-FREE; and two comparable figures were corrected under
`git show --numstat` (0faeb416 tests 443, 33f6fdc non-test 245), which is what
the raw-currency choice exists to make possible. Full per-finding records are
in an untracked, gitignored revision-notes memo.

## Landing note, S1 (appended 2026-08-15)

Five layers as specced, one virtual and one refusal site. Engine:
`Chain::fitsWithoutOffset(double*)` public and non-const beside
`forestTotalFits`, opening on `numReportedLocations() > 1` and writing
`fitScale * combinedFits + fitShift`; `Sampler::fitsWithoutOffset(chainNum,
double*)`; the `SamplerBase` pure virtual plus the `SamplerFacade` override on
the `forestTotalFits` pattern (`SamplerFacade` is the only implementer, so the
new pure virtual breaks nothing). `storeSample`'s training write splits on
`numLocations == 1`: that branch calls the accessor and then adds the offset as
a SEPARATE `misc_addVectorsInPlace` pass, never fused; the multi-location
branch keeps the original loop verbatim. Bridge: `bartcore_getFitsWithoutOffset`
(one arg), mapping the engine's `false` to a named error, one prototype, one
`DEF_FUNC`; no second shape test. R: `$getFitsWithoutOffset()`,
`bartcoreFitsWithoutOffset`, and `refuseHostRead`, a read-side sibling of
`refuseHostMutation` called from the new method only. Shape: an
n.observations x n.chains matrix always, `$getForestFits`'s convention rather
than `$getLatents`'s single-chain vector.

Docs: the per-family `getLatents` semantics on all three surfaces (Rd `\value`,
R5 docstring, the comment above `dbarts_sampler_getLatents` - comment-only, no
hash move), read off `model.hpp`'s overrides rather than asserted (locations
`ProbitResponse::latents_`, `OrdinalResponse::latents_`, `AFTResponse::logT_`;
precisions `LogisticResponse::omega_`, `NBResponse::omega_`,
`TResponse::lambda_`; `GaussianResponse` and `MultinomialResponse` declare
none). The fit-surface table and the boundary sentence went into a new
`\subsection{Reading the fit}` of the Rd's `\details`. The BINDING INSTRUCTION
was discharged: `extract.bart` ([[R/generics.R:211@70b30c18]]) reads `object$yhat.train` /
`object$yhat.test` directly and NO `type` arm removes the offset - "bart"
returns those stored draws as they stand, "ev" maps them through
`probabilityFromLatents`, "ppd" samples from them, "loglik" evaluates against
them - so every arm is offset-INCLUSIVE; `extract.dbartsSampler`
([[R/generics.R:1449@70b30c18]]) takes `type = "predictors"` only and returns the design
matrix, not a fitted value. Both stale-doc corrections landed: the
`\item{active}` veto clause is now conditional (count-based accounting is
`numObservations`, the scan count and leaf collapsing; the veto counts
POSITIVE-WEIGHT members, degenerating to the member count only with no weight
vector installed - and a mask IS one), and `r-c-division.md`'s
latent-subset-mask bullet names 93afd635 instead of "this slice".

Budget, raw additions against 94747480: non-test 184 + this note, tests 362
(tinytest 290, `tests/cpp` 72). Stops 389 / 518.

Battery, all on the slice's private lib (`/tmp/s1-lib`, since `~/.Renviron`
overrides `R_LIBS_USER`): `--preclean` install clean; `tests/cpp` from
`make clean` 243 ok, all passed, and again under
`OPT="-O2 -g -fsanitize=address,undefined"` with
`ASAN_OPTIONS=detect_container_overflow=0`, no sanitizer diagnostic; full
tinytest 4919 results, 0 failures, no snapshot re-pinned; equivalence trio
BITWISE with no re-record (equivalence-8b047f8b 37 compared / 0 skipped, every
scenario "identical draws (same RNG stream)"; bcf-equivalence-8b047f8b every
channel identical on all 12 scenarios; multinomial-equivalence-1027be5 every
channel identical on all 10 - no max-|z| line anywhere), which is this slice's
leak detector for the `storeSample` refactor and reports none; `air format
--check` clean; `lintr::lint_package()` 0 lints; `pkgdown::check_pkgdown` no
problems (no new topic, so no `_pkgdown.yml` entry was owed); `R CMD check
--as-cran` from a tarball staged outside the tree, Status OK, 0/0/0.

Mutation proof, each applied, run, then reverted and byte-verified (the R
identity cells are 11 value comparisons plus 22 shape/null assertions; the two
refusal cells and the two recombination cells complete the 47).

- (1a) Added offset at an ACCESSOR-ONLY site (the R5 method): 13 assertions
  move - all 11 identity comparisons and both recombination cells. No fixture
  is vacuous.
- (1b) The same added offset inside `Chain::fitsWithoutOffset`, the SHARED
  site: every identity cell stays GREEN, because `storeSample` then
  double-adds and the identity holds. **This is the refactor's blind spot and
  the plan's mutation-1 wording does not distinguish the two sites.** It is
  caught anyway, by the arithmetic gates: both `tests/cpp` assertions and both
  recombination cells go RED (plus the shipped "sums of squared residuals
  descale" pin). The `tests/cpp` cell asserts the recorded training write
  equals the accessor plus the offset precisely to cover this.
- (2) Return the internal scale: all 13 identity cells GREEN as predicted; RED
  on the `tests/cpp` cell, both recombination cells, and 22 other `tests/cpp`
  pins (the training channel is engine-wide).
- (3) Report forest 0's totals instead of the combiner blend: EXACTLY ONE
  R assertion moves, the BCF recombination cell, plus the shipped
  "BCF training fits are the a*mu + b_z*tau blend" cell in `tests/cpp`.
  Identity cells, single-forest recombination and the new `tests/cpp` cell all
  GREEN, confirming the BCF cell is load-bearing and must not be tidied away.
- (4) Delete the `numReportedLocations() > 1` early return: exactly the two
  multinomial refusal assertions move; `tests/cpp` fully green under ASAN too,
  which is the point - the read is in-bounds and silent.
- (4b) Delete the read-side `hostFor` guard: exactly the shell cell moves.
- (5) Reverse the location/precision assignment in the Rd: caught by NO test,
  as specced; verified by reading `model.hpp`'s `latents()` overrides.

Deviations and carried obligations.

- The Rd first drew a `\link{student}`, which `R CMD check --as-cran` flagged
  as a missing link (`student` is unexported and has no topic); it reads
  `\code{student()} (see \code{\link{dbarts}})` instead.
- The identity helper cannot hold the expectations: `tinytest::expect_*` is
  declined for registration ("Remove the 'tinytest::' prefix"), while an
  UNQUALIFIED `expect_*` inside a helper is an `object_usage_linter` finding
  under `lintr::lint_package()`. The helper therefore gathers and the three
  assertions sit at top level. Note for later slices: the shipped
  `tinytest::expect_identical` in `[[test-mutate-sparse-valued.R:48@70b30c18]]` is
  unregistered for the same reason.
- The `tests/cpp` cell carries one assertion beyond the plan's - the recorded
  training write equals the accessor plus the offset - which is what makes
  mutation 1b visible.
- `docs/design/feature-matrix.md`'s file:line anchors into `R/dbarts.R` and
  `R_interface_bartcore.cpp` shifted and were NOT refreshed here, following the
  previous arc's pattern of one anchor-refresh commit at the END of the arc
  (cd27822f). Owed before the arc closes.

Landed: 5bedf923, pushed 2026-08-15. Independent gate-runner CONFIRM on its
own lib (243/243 tests/cpp plain and ASAN, tinytest 4919/0 failures, trio
bitwise 37/12/10, R CMD check --as-cran OK 0/0/0, provenance verified);
orchestrator diff review against the spec found no deviation beyond those
recorded above.

## Landing note, S2 (appended 2026-08-15)

Five layers, no `dbarts.h` touch (F2 = (C); the flat half stays S8's). Engine:
`Results::dispersion`, a numSamples channel on the `Results::cutpoints`
precedent (no `dbarts_results` counterpart, so no hash move); `Chain::dispersion`
and `::carriesDispersion`, one-line forwards to `response_->` beside
`numCutpoints`, and the `Sampler` pair; `Sampler::run`'s per-chain stride
`+ c * numSamples`, sigma's own; and `Chain::storeSample`'s write, guarded
`results.dispersion != nullptr && response_->carriesDispersion()`, placed after
the cutpoint write and copying the `results.tau` shape. Nothing new is drawn -
`NBResponse` already held r - so it is a pure read of state the sweep settled
on. `SamplerShape::carriesDispersion` plus a `SamplerBase::dispersion(chainNum)`
pure virtual and its `SamplerFacade` override carry it to the bridge
(`SamplerFacade` is still the only implementer).

Bridge, the bulk of the work: the slot inserts directly after the ordinal
cutpoint slot, so `numResultSlots` gains its term and `varianceTrainSlot`,
`varianceTestSlot`, `forestFitsSlot` and `glueSlot` each gain
`(hasDispersion ? 1 : 0)`, with the keyed `SET_STRING_ELT` names following. No
family carries both channels; the arithmetic composes regardless. Shape is
sigma's. Plus `bartcore_getDispersion` on `bartcore_getSigmas`'s pattern,
returning `R_NilValue` off a family carrying none, one prototype, one `DEF_FUNC`.

R: `$getDispersion()`, argumentless, docstringed, beside `getSigmas` as its
count analog, guarded by `refuseHostRead` on S1's precedent (a NEW read is
guarded; the shipped four stay unguarded by ticket `host-shell-read-guards`).
The non-nbinom answer is `NULL`, not an error, which is why the Rd states the
`!is.null` law and every assertion in the test file leads with a shape check.

**`bart2Negbin` WAS touched**, under the spec's conditional permission. Its
one-kept-sample-at-a-time driver STAYS - the sweep structure is untouched and
the `nbinom` equivalence scenario does not move - but its per-sweep read
switched from `bartcoreStoreState(bc)` to the run's own `r$dispersion` row.
Draw-neutral by construction (`bartcore_storeState` reaches `Sampler::getState`,
which draws nothing and mutates nothing, and `state` was read for the dispersion
and nothing else) and by measurement (the trio's `nbinom` scenario reports
identical draws). Simpler by a line, and it removes a FULL tree-state
serialization from every sweep of every `bart2(family = "nbinom")` fit, which is
the channel's whole point: the package now reads it the way its documentation
tells a host to. The stale comment above `bart2Negbin` ("the engine has no
per-sample dispersion output channel") was corrected in place.

Docs: the `run` `\value` paragraph gains the element with its shape, its two
r-modes and the absent-not-NULL law; a `\value` entry for `getDispersion`; the
alias and `\usage` line; one NEWS item.

Budget, raw additions against 7e6a4a3e: non-test 116 + this note (total 244,
band 160-215, stop 270), tests 293. Both overshoot, both for the same reason -
cell count. The test oracle is {run slot, getter, state} x {fixed r, estimated
r} x {1 chain, 2 chains}, twelve comparisons before the shape, NULL,
multi-sample, ordinal-neighbour, BCF, heteroscedastic, host-shell and end-to-end
cells, and `air`'s vertical call formatting inflates raw lines against a dense
reading; the note carries four mutations and five deviations. The spec's
"engine ~10" costed the `storeSample` write alone and not the shape flag, pure
virtual and facade override the R5 getter it also mandates. No `tests/cpp`
cell; the `thresholds` precedent carries none either.

Battery, all on the slice's private lib (`/tmp/s2-lib`, since `~/.Renviron`
overrides `R_LIBS_USER`): `--preclean` install clean; `tests/cpp` from
`make clean` 243 ok, all passed, and again under
`OPT="-O2 -g -fsanitize=address,undefined"` with
`ASAN_OPTIONS=detect_container_overflow=0`, 243 ok, no sanitizer diagnostic;
full tinytest 4985 results, 0 failures, no snapshot re-pinned (4919 at S1 plus
this file's 66, so no shipped count moved); equivalence trio BITWISE with no
re-record (equivalence-8b047f8b 37 compared / 0 skipped, every scenario
"identical draws (same RNG stream)" - the `nbinom` scenario among them, this
slice's leak detector for the `bart2Negbin` read swap, reporting none;
bcf-equivalence-8b047f8b identical on all 12; multinomial-equivalence-1027be5
identical on all 10; no max-|z| line anywhere); `air format --check` clean;
`lintr::lint_package()` 0 lints; `pkgdown::check_pkgdown` no problems (no new
topic, so no `_pkgdown.yml` entry was owed); `R CMD check --as-cran` from a
tarball staged outside the tree, Status OK, 0/0/0.

Mutation proof, each applied, run, reverted and byte-verified by `cmp`.

- (1) Write the dispersion one sample late (a cached previous r, lazily
  initialized so the first store is still correct): EXACTLY TWO assertions
  move, the estimated-r arm's `slot == state` comparisons at one and two
  chains. Every fixed-r cell stays green, as specced - a repeated constant
  cannot lag - which is why that arm is mandatory and carries the
  at-least-two-distinct-values measurement. The GETTER cells stay green too,
  and correctly: it reads live state and is a separate instrument. (A first
  attempt left the cache uninitialized, moving ten assertions including the
  fixed arm; that is a different defect, and the two above are the faithful count.)
- (2) Fill the slot with the shipped grid's median (8.0): 8 assertions move,
  all on the estimated arm (both `slot == state` cells, all three non-vacuity
  cells, the multi-sample grid and distinctness cells, and the `bart2`
  end-to-end `mu = r exp(psi)` cell), plus 3 in the shipped `test-nbinom.R`.
  The fixed arm stays GREEN by construction - its r IS 8, which is why that
  value was chosen for it.
- (3) Off-by-one in the downstream slot indices (`varianceTrainSlot` and
  `forestFitsSlot` each +1, `numResultSlots` left correct): the heteroscedastic
  run errors "attempt to set index 10/10 in SET_STRING_ELT" and the BCF run
  "attempt to set index 10/10 in SET_VECTOR_ELT", while the nbinom and gaussian
  runs stay green and correctly named. Exactly why two models carrying no
  dispersion at all are in this battery.
- (4) Drop both guards so a gaussian run writes through the unallocated
  channel: under `-fsanitize=address,undefined` UBSAN names it precisely -
  "runtime error: store to null pointer of type 'double'" at `[[chain.hpp:5054@70b30c18]]` -
  and ASAN follows with `SEGV on unknown address 0x0` after 35 ok lines.
  **Deviation from the plan's wording:** it is not sanitizer-ONLY. The plain
  binary also dies (exit 139) at the same write; the sanitizer leg adds the
  file, line and diagnosis rather than a bare signal. "Caught by ASAN only"
  would hold for an in-bounds-but-wrong write, not for a null one.

Deviations and carried obligations.

- `R CMD check --as-cran` flagged `\link{storeState}` as a missing link (it is
  a reference-class method alias, not a topic); the sentence reads
  `\code{storeState}`. Same shape as S1's `\link{student}` finding.
- `bart2Negbin` was touched, with the reasoning and proof above; the driver
  itself is unchanged.
- `refuseHostRead` on `$getDispersion()` is not in the spec's letter. It
  follows S1's precedent and will matter once S7 puts `hostFor` on bart2's
  nbinom host, where an unguarded read would return the placeholder shell's r
  as if it were the fit's.
- The multi-sample cell's "last recorded draw equals the post-run state" did
  NOT move under mutation 1 (r happened not to change on that final step). The
  sample stride is gated instead by the grid-membership and distinctness
  assertions over the whole slab, which do move; the per-sweep cells are the
  load-bearing ones.
- `docs/design/feature-matrix.md`'s file:line anchors into `R/dbarts.R`,
  `R/bart.R` and `R_interface_bartcore.cpp` shifted again and were NOT
  refreshed here, following S1 and the previous arc's one anchor-refresh commit
  at the END of the arc. Owed before the arc closes.

Landed: da3c76f9, pushed 2026-08-15. Independent gate-runner CONFIRM on its own
lib (243/243 tests/cpp plain and ASAN, tinytest 4985/0 failures, trio bitwise
37/12/10 incl. the nbinom scenario, R CMD check --as-cran OK 0/0/0, provenance
verified; weak-gate and mandated-cell audits clean). Orchestrator adjudication
of the test-budget overage (293 raw vs the ~150 estimate, past the nominal
1.5x): ACCEPTED - every cell traces to the spec's mandated oracle grid or its
mutation-3 index neighbors, and budgets size to the MANDATED ORACLE, not the
estimate; the estimate was low the same way the memo prices this arc corrected
were.

## Landing note, S3 (appended 2026-08-15)

One guard, stated once and called from four entries. `refuseGroupedScaleUpdate`
(declared beside its two siblings in `R_interface_bartcore_common.hpp`, defined
under `refuseVarianceForestScaleUpdate`) replaces `bartcore_setResponse`'s
blanket numGroups refusal and is added to `bartcore_setOffset`, closing the
unrefused decalibration `GroupedResponse::setOffset`'s own comment admitted.
The predicate is the plan's - grouped, re-anchoring base (gaussian, which is
Student-t's report, or aft) - with one keying deviation: ANYTHING BUT FALSE is
refused, `refuseVarianceForestScaleUpdate`'s condition-keying, not the literal
`== TRUE`. Reason: the two surfaces convert `updateScale` to the engine's bool
differently (the R bridge on `== TRUE`, the flat API on `!= 0`), and FALSE/0
is the only value both read as "pin it"; the sole behavioural delta is NA, now
refused. Grouped probit and logistic take TRUE as the documented no-op;
`setData` is untouched.

**The flat C API is guarded too**, a scope addition the plan did not carry:
`dbarts_sampler_create` reads the same `bartcore.groups` attribute, so a
flat-API sampler can be grouped, and both sibling scale guards already guard
both surfaces - guarding one layer would have left the identical defect open
one layer down. `dbarts_sampler_setResponse` and `_setOffset` call the same
function; no `dbarts.h` declaration moved, so no ABI or apiHash event.

`model.hpp`: three comments asserted the blanket refusal this slice replaces
and were corrected in place (the class comment, `::setOffset`'s,
`::setData`'s), 9 raw lines - mandated by truthfulness, not the plan.

The SBC swap arm. `runSbcGrouped(swap = TRUE)` builds ONE fit on
`config$yBuild`, captures state0 by the two-line `storeState()`/`state` idiom
the plan gives, and per replication installs state0, re-overdisperses
(`sampleTreesFromPrior`, `sampleNodeParametersFromPrior`, `setSigma` where the
family has one) and swaps y0 in at the pinned scale. CLI configs
`grouped-gaussian-swap`/`grouped-probit-swap`; NOT in sbc.yaml's matrix, per
the plan. The state install restarts the chain rng at the same point each
replication, so each rank is a deterministic function of that replication's
own iid (theta0, y0) draw - ranks stay iid across replications, recorded in
the arm's comment.

Wall-clock, measured: pilot at R = 10, swap arm 3.094 s/rep against 3.086
rebuild and 2.746 ungrouped - decorator multiplier 1.13x, swap-vs-rebuild
1.00x. The plan's 78000 sweeps/rep over-counts: the harness's burn default
reaches `run()` as burn-in sweeps, so a replication is 8400 sweeps and the
full arm ~11 minutes, not ~13+.

Verdicts, in the pre-registered order. UNMUTATED first, R = 200: all 10
functionals PASS, tau and b uniform - THE FALSIFIER DID NOT FIRE, the refusal
stands as designed. Mutation 1 (predicate widened so the arm re-anchors),
R = 200: 3 FLAGs against 0 - tau (ecdf 0.0975 over the 0.0917 band), sigma
(0.2412, chisq 0.000) and f.star1. **The plan's prediction half-landed**: it
said tau AND b go non-uniform; tau and sigma move, b does not, because the
tau prior scale is converted exactly once at construction while b's conjugate
conditional re-equilibrates against the new scale within the burn. The gate
discriminates cleanly regardless; the tool-verified-claims discipline catches
another plan-stated behavioural prediction.

Other mutations, applied, run, reverted, source verified clean: (2) drop
`rebuildWorking()` from `GroupedResponse::setResponse` - 60 RED, every one the
tests/cpp elementwise "new z minus b" assertion, nothing else, the private
mechanism the plan named; (3) predicate narrowed to gaussian - exactly the two
grouped-aft error cells RED; (4) `setResponse` guarded but not `setOffset` -
exactly the three `setOffset` error cells RED.

Tests. tests/cpp `testGroupedResponseSwap` runs on a PRIVATELY seeded rng so
the shared stream's downstream fixed expectations stay unmoved; it pins both
halves - the pinned swap leaving b, tau and sigma untouched and rebuilding the
working response as the reference base's new z minus b, exact equality, and
the engine-level updateScale defect itself (sigmaScale moves, sigma converts,
b and tau keep their numbers). tinytest `test-grouped-swap.R` (26 results):
the pinned swap moves the fit TOWARD the new response (the correlation flips
sign); refusals matched on "tau" AND "random intercepts b"; the probit no-op
is a BITWISE two-sampler comparison at a fixed rngSeed, both conduits; the
shipped `test-rbart-bartcore.R` cell replaced with its setData and save/load
neighbours untouched. **The agreement cell deviates from the plan's letter to
its own benefit**: the swapped sampler is created on a PERMUTATION of y2, not
on a different y1, so sd and hence the creation transform, leaf calibration
and tau prior scale are IDENTICAL across the pair and the two chains target
exactly the same posterior; created-on-y1 would have compared two subtly
different targets.

Docs: the latent-subset-mask.md rule-2 block landed verbatim, 5 lines for 5,
EOF errata only (numstat 11/5); the Rd's `updateScale` paragraph gains the
grouped clause and a "Grouped random effects" subsection states the law; one
NEWS item; the NEWS db builds.

Budget, raw additions against 05a2b73a: non-test 104 + this note (band
130-175, stop 220 - under; the plan's ~90 records allowance is this note);
tests 316 against ~223, stop 335. Overage ADJUDICATED ACCEPTED: the
independent gate-runner traced every cell to a mandated oracle or its direct
extension (the second error-pattern check, post-swap ranef/tau finiteness,
correlation-direction, pinned-offset-accept cells), none scope creep.

Battery, twice, separate libs (implementer /tmp/s3-lib, gate-runner
/tmp/s3-gate-lib): `--preclean` install; tests/cpp from `make clean` 244 ok
plain AND under `-fsanitize=address,undefined` with
`detect_container_overflow=0`, zero diagnostics; tinytest 5012 results, 0
failures (4985 at S2 + 26 this file + 1 net from the replaced cell's split);
trio BITWISE, no re-record (equivalence-8b047f8b 37/37, bcf-equivalence-
8b047f8b 12/12, multinomial-equivalence-1027be5 10/10, no max-|z| line
anywhere); air clean; lintr 0 on all three touched R files against the slice
lib; pkgdown no problems (no new topic); `R CMD check --as-cran` from a
git-archive tarball staged outside the tree, Status OK, 0/0/0.

Landed: eeedc07c, pushed 2026-08-15. Independent gate-runner CONFIRM (all nine
gates, four audits: oracle traceability, weak-gate, rule-2 numstat, no
docs/plans in code comments). feature-matrix.md's grouped-row cells updated in
the records commit; the arc-end scoped anchor refresh stays owed.

## Landing note, S4 (appended 2026-08-15)

The API as specced: `dbartsDrawLatents` and `dbartsWorkingResponse` exported,
`match.arg` over the six family names, `fit` offset-free per the engine's
`totalFits` convention, the latent's quantity attribute, the two DEF_FUNC
lines, one Rd topic `dbartsAugmentation` with per-function aliases, the
`_pkgdown.yml` entry under "The mutable sampler", one NEWS item. The six draw
laws and six working responses are RESTATED in the bridge, not shared with the
response models (the spec's OUT OF SCOPE ban; the models run against a
different generator), and `model.hpp` is absent from the diff entirely.

The RNG contract and the first-caller obligations, all four discharged and
independently audited: the per-call generator is init-captured as a
`unique_ptr` with an `ext_rng_destroy` deleter INTO the house `unwindProtect`
closure, whose R_UnwindProtect cleanup deletes the heap-held closure on
Rf_error's longjmp as well as the normal return - ownership by mechanism, not
by claim; the NULL-creation path errors after `PutRNGstate`; nothing is
cached; the Rd carries the thread-safety sentence, both halves of the
stream contract, and the offset-discrimination sentence.

Implementation deviations, all adopted: no second enum - "student" maps to
`ResponseFamily::gaussian`, which then NAMES the scale-mixture arm of both
switches; `dbartsWorkingResponse` skips `validateResponseSupport` for ordinal
only (y never enters that arm and there is no K to state support against);
`sigma`, carrying a default, is refused off aft/student on `!missing(sigma)`
rather than non-NULL; `weights` and `sigma` are refused off-family but
OPTIONAL within it, while `dispersion`/`thresholds`/`df` are two-way (their
family REQUIRES them - the law has no fallback); the R layer adds
positivity/wholeness/sortedness validation the spec implies but does not
itemize, safe-over-fast.

Oracles, measured (seed 5, N = 20000 per cell, threshold |z| < 4):
truncated-normal moments at 10 probit points both tails, 4 aft points, and -
added post-audit, see below - 2 ordinal points (two-sided category z = 2.553,
one-sided z = -0.291), max |z| = 3.19, MC errors 0.0018-0.0067, variances
within 2-5 percent; PG means at psi in {0, 0.5, 2, 5} for b = 1, 3, 5
including the psi -> 0 limits b/4 and b/24, max |z| = 1.61; the Student-t
Gamma mixer at df 3 and 10. Agreement: an Albert-Chib probit loop built ONLY
from the two helpers against the engine's posterior, cor 0.9967 / latent rms
0.057 / probability rms 0.016, plus a LOGISTIC arm (0.9963 / 0.081 / 0.012;
worst over three more seeds 0.9952) - beyond the spec's probit-only letter,
but mutation 4's kill ("identity for logistic - RED on the agreement arm
only") is unreachable without it, so the spec's own mutation matrix requires
the arm. Offset discrimination: latents `identical`, working responses
differ, on a location AND a precision family.

Mutations, each built to its own lib, run, reverted, revert verified: (1)
wrong-side truncation - 32 RED, all probit moment and agreement cells, the
Albert-Chib slope INVERTS (cor -0.993), logistic untouched; (2) PG shape
ignored - 16 RED, exactly the weighted-logistic (weights 3) and nbinom
(y = 2, r = 3 - the spec's y+r > 1 constraint honoured) cells, every shape-1
cell green; (3) Get/PutRNGstate dropped - 1 RED, the stream-ADVANCE half
only, the reproduce half stays green, which is why the cell asserts both;
(4) logistic working response as identity - 5 RED, the logistic agreement
cells plus the direct quotient identity; (5, added post-audit) ordinal
interior arm drawn at mean 0 instead of psi - 1 RED, the new interior moment
z moving +2.553 -> -8.305 while containment stays GREEN, which is precisely
the gap the added cell closes.

Post-audit amendment, one commit kept: the independent gate-runner flagged
the ordinal law as the one cell instrumented by CONTAINMENT alone - true on
any draw inside the interval - so a 16-line two-sided/one-sided moment cell
was added under oracle 1's mandate, proven discriminating by mutation 5, and
the commit AMENDED before push. Full tinytest and `R CMD check` re-ran on the
amended tree; the C side did not change, so the ASAN and trio legs stand from
the pre-amend run.

Gate corrections this slice surfaced, both now standing battery text: the
NEWS-parse invocation `tools:::.build_news_db("dbarts", ".")` is a SILENT
NO-OP ("." binds to lib.loc; it returns NULL without parsing) - the working
forms are `.build_news_db_from_package_NEWS_Rd("inst/NEWS.Rd")` or an
explicit `lib.loc`, gated on a non-NULL entry count (234 here). And the
spec's valgrind leg CANNOT run on this host or this branch - darwin/arm64
has no valgrind, and valgrind.yaml, like every schedule/dispatch workflow,
is unregistered off bartcore - so the memory-safety gate was the local ASAN
tests/cpp leg plus the per-push CI sanitizers workflow (full tinytest, this
file included, under two instrumented R builds), with the memcheck leg
riding the first nightly after bartcore reaches main. Recorded, not waived.

Budget, raw additions against f8e630fe, BOTH SIDES PAST THEIR STOPS -
adjudicated, not slipped: non-test 517 (bridge 218 with the hpp and
DEF_FUNC lines, R 164, Rd 122, NEWS 11, NAMESPACE 1, pkgdown 1) against the
spec's ~336-with-records band 300-390 stop 505; tests 418 against ~230 stop
345. Compaction was performed BEFORE the report (bridge 241 -> 211 dropping
a second enum and a spec struct, tests 512 -> 402, R 192 -> 164), no
assertion was dropped to fit. ADJUDICATED ACCEPTED on the gate-runner's
line-complete audit: the bridge inventories to exactly the mandated content
(two entries, family map, marshalling, six laws, six responses, RNG
plumbing - no cache, no refactor, no extra export), and all 57 pre-amend
test call sites trace to a mandated oracle, a mutation-kill arm the spec's
own matrix requires, or a validation the shipped R code performs. The ~110
bridge estimate was scaled from a five-way GUARD - which branches to a
message - while this branches to six SAMPLING LAWS; the same
estimate-vs-mandate correction as S2's, and the spec itself flagged its
figure ("the honest total is ~566 raw").

Battery, twice, separate libs (implementer /tmp/s4-lib, gate-runner
/tmp/s4-gate-lib, pre-amend; the amend re-ran tinytest/air/lintr/check):
`--preclean` install; tests/cpp from `make clean` 244 ok plain AND under
ASAN/UBSAN, zero diagnostics; tinytest 5131 results, 0 failures (5012 at S3
+ 119 the new file); trio BITWISE, no re-record (equivalence-8b047f8b
37/37, bcf-equivalence-8b047f8b 12/12, multinomial-equivalence-1027be5
10/10); air clean; lintr 0 on both touched R files; NEWS parses (234
entries, by the corrected invocation); `pkgdown::check_pkgdown` no problems
with the new topic's entry verified in place; `R CMD check --as-cran` from
a git-archive tarball staged outside the tree, Status OK, 0/0/0.

Landed: 890efd3d, pushed 2026-08-15. Independent gate-runner CONFIRM
(pre-amend tree; all nine gates, audits A-D including the budget audit
above). No feature-matrix cell bears on the helpers; none updated. The
arc-end scoped anchor refresh stays owed.

## Landing note, S5 (appended 2026-08-15)

R-ONLY, stated plainly: the commit touches NAMESPACE, `_pkgdown.yml`,
`inst/NEWS.Rd` and three new files (`R/validateComposition.R`,
`inst/tinytest/test-validate-composition.R`,
`man/dbartsValidateComposition.Rd`); `src/` is absent from the diff, so the
engine binary is e650ab8c's - the 33f6fdc precedent the MANIFEST records.

The five properties, each verified in code by the independent gate-runner
with line citations: (1) the ranked quantity is
`functionals(init(theta0, simulate(theta0)))`, theta0 passed opaquely and
never name-indexed; (2) `compositionFunctionals` errors on a length or name
disagreement, checked at "at init" AND "after a step" against the FIRST
replication's reference - later replications' inits are checked too, a
strengthening over the letter that catches a host whose init is stochastic
in shape; (3) `sbcDiscreteRank` ported and applied to every functional's
rank; (4) the `withFixedSeed` guard restores `.Random.seed` via `on.exit`
(error paths included) and REMOVES the seed it itself created when the
caller had none; (5) Bonferroni over this call's functional count. The
returned `dbartsCompositionValidation` carries the rank matrix, L, the
verdicts frame and both alphas, with a registered S3 print method.

**The spec's mutation-3 prediction was WRONG, and the correction is the
note's headline**: the band is `quantile(nullMax, 1 - alpha)` compared as
`observed <= band`, and `nullMax` - a max-absolute-ecdf-difference at
R = L = 200 - is bounded near 0.15-0.2 under the null, so NO quantile of
it, either tail, can reach the known-bad cell's 0.53 statistic. Taking the
wrong tail SHRINKS the band and reddens the calibrated cells instead. Both
directions are gated regardless: M3 as specced (wrong tail: 10 RED on
known-good/derived/many-functional) plus M3b, added at implementation (band
inflated 5x: 4 RED, three on the OVER-DISPERSED cell at 0.184 - the cell
actually sensitive to a loose band; known-bad at 0.53 is out of reach of
any null-shaped band). The gate-runner checked the algebra against sbc.R's
identical band formula, confirmed the prediction unreachable, noted the one
alternate reading (max -> min in the OBSERVED statistic) is a different bug
locus than "quantile tail", and recommended accepting M3b as the intent
realized. The tool-verified-claims discipline catches a third plan-stated
behavioural prediction this arc (after S3's M1-b and S2's mutation-4
wording).

Other mutations, each built to its own lib, run, reverted, revert
re-verified 55/55: (1) `<=` for `<` in the rank - 3 RED on the derived/atom
cell including the mandated mean-rank-vs-L/2 assertion, which also reddens
alone; (2) Bonferroni dropped - 4 RED, the many-functional cell only;
(4) name-intersection against drawPrior's return - the derived cell aborts
on an NA rank while the plain cell reproduces the clean run exactly, which
is the spec's point about derived functionals being unrankable by name;
(5) `.Random.seed` restore removed - 2 RED, the seed cell only (both its
halves: restore-when-present and remove-when-absent).

Oracle cells as specced, on the analytic conjugate normal-normal step, no
dbarts sampler in tinytest: known-good PASSES every functional; known-bad
(posterior mean as draw) FLAGS with `ecdf.diff > band` asserted, not just
"ran"; over-dispersed (2x posterior sd) FLAGS - both directions; the
derived-functional, many-functional (M = 8), seed and validation-error
cells complete the set. Cell parameters carry their spec-assigned owners:
`n.thin = 1L`/`n.burn = 0L` justified in the TEST FILE (exact-posterior
step, exact init), the 200/200 defaults justified in the Rd (approximate
hosts). Timings, measured: every ordinary cell under 0.15 s, the
many-functional cell 0.36 s, the FILE 1.1 s against its 45 s budget (1.2 s
inside R CMD check); the n.replications = 100L fallback was never needed.

Budget, raw additions against e650ab8c: non-test 452 (R 297, Rd 141, NEWS
8, NAMESPACE 3, pkgdown 3) - over the 280-380 band, UNDER the 525 stop;
tests 267 against ~200, under the 300 stop. The overshoot is the
11-argument Rd (141 vs ~70) and air's vertical formatting of the driver's
signature and validation ladder; compaction ran before the report
(R 316 -> 297), nothing dropped to fit. No stop crossed, no adjudication
owed; the gate-runner's cell audit still traced all 55 assertions to the
mandated set with one structural extra (unnamed-functional labeling,
traceable to the documented labeling behavior).

Battery, twice, separate libs (implementer /tmp/s5-lib, gate-runner
/tmp/s5-gate-lib): install; tests/cpp from `make clean` 244 ok (no-change
verification); tinytest 5186 results, 0 failures (5131 at S4 + 55); trio
BITWISE, no re-record (equivalence-8b047f8b 37/37, bcf-equivalence-8b047f8b
12/12, multinomial-equivalence-1027be5 10/10); air clean; lintr 0 on both
new R files; NEWS parses by the corrected invocation (235 entries);
`pkgdown::check_pkgdown` no problems, the entry under Diagnostics;
`R CMD check --as-cran` from a git-archive tarball staged outside the
tree, Status OK, 0/0/0. No ASAN leg owed (R-only), and the blind-spot
sentence - a drawPrior/simulate mutual inconsistency is undetectable -
ships in the Rd as mandated.

Landed: 37d9ec81, pushed 2026-08-15. Independent gate-runner CONFIRM (all
nine gates, audits A-D, M3-deviation algebra checked and accepted). No
feature-matrix cell bears on the validator; none updated. The arc-end
scoped anchor refresh stays owed.

## Landing note, S6 (appended 2026-08-15)

R/docs-ONLY: the commit touches `_pkgdown.yml`, `inst/NEWS.Rd`,
`man/bart.Rd`, `man/dbartsSampler-class.Rd` (the two `\seealso` legs) and
three new files (`vignettes/dbarts-as-a-component.Rmd`,
`man/dbarts-embedding.Rd`, `inst/tinytest/test-embedding-recipes.R`);
`src/` and `inst/include/` are absent from the diff, so the engine binary
is 8a15d0ad's. NAMESPACE is untouched, correctly - the page and the
vignette export nothing.

The six recipes ship as specced, each with its tinytest cell, all six
traced by the gate-runner with none under-asserting: (1) AUGMENTATION,
the Albert-Chib loop over S4's helpers against the engine's probit
posterior; (2) OFFSET BLOCK, both blocks' truth plus the seeding-contract
identity cell, AND the recipe's own `dbartsValidateComposition` run -
clean PASS, mean-substituted FLAG - which is the end-to-end proof S5
works on a real composition; (3) K-FOREST with the sqrt(K)
`$setCalibration` rescale; (4) LATENT COVARIATE, the install mask read as
an MH accept mask; (5) OUTER-OWNED SIGMA, the pin surviving a sweep and
the guard message; (6) INCREMENTAL USE, closing defect 5's measured trap -
the vignette pins `z - train = (z - f) - offset` numerically and the cell
gates the slope bias at 6 oracle SEs.

Non-vacuity, measured where the spec demanded it and once more where it
did not: M1's carrier MOVED as a result. The posterior-fit comparison at
n = 200 is likelihood-dominated and separates a dropped sqrt(K) rescale
by only ~1.2x, so the single-sampler fit cell ships green as specced but
the MUTATION rides the composed PRIOR sd (clean ratio 0.9706, mutated
1.3726 against sqrt(2) = 1.414; tolerance 0.15 set AFTER measurement,
margins 5.1x below and 2.5x over). `samplePriorPredictive` was rejected
as the instrument - it reads the model object, which `setCalibration`
never rewrites, so it is blind to the lever. Recipe 4 runs at n = 80
with `cgm(power = 0.5)` because at n = 200 the install mask goes
all-TRUE and the cell would be vacuous; the cell asserts its own
non-vacuity as a hard conjunction, so an all-TRUE mask FAILS rather than
silently passing.

Mutations: M1 above (2 RED); M2, `run()$train` for
`$getFitsWithoutOffset()` - 2 RED, clean -3.17 SE vs mutated -16.78 SE
against the 6 SE threshold, 13.6 SE separation (the plan predicted 7-10;
the direction and the gate hold, the magnitude was larger); M3, posterior
mean for a draw in recipe 2's outer block - the validator FLAGS
(ecdf.diff 0.0863 -> 0.3824 against band 0.2225), 2 RED; M4, a broken
`\link{}` - caught by `R CMD check`'s cross-reference leg only, invisible
to tinytest, run on a staged copy, never the worktree. M1-M3 live in the
test file's own recipe copies, so no mutated install was needed; every
revert verified 32/32 green.

One documentation-only gap, recorded not fixed: recipe 6's prose quotes
the MUTATED variant's numbers ("about 0.29, roughly seventeen oracle
standard errors"), which the gate-runner reproduced exactly (0.2856,
16.78 SE) but which no shipped test pins - the mutated variant never
runs in CI. A future RNG-stream change can drift those two hedged words;
the mechanism itself stays gated by the recipe's 6-SE cell.

Budget, raw additions against 8a15d0ad: non-test 474 (vignette 349,
embedding Rd 108, NEWS 9, `\seealso` 7, pkgdown 1) - over the 300-400
band, UNDER the 555 stop; tests 428 against ~360, under the 540 stop.
The overshoot is the vignette: the plan's ~210 priced six recipes but
not recipe 2's SBC harness (45 lines) that the spec itself mandates. No
stop crossed, no adjudication owed.

Battery, twice, separate libs (implementer /tmp/s6-lib, gate-runner
/tmp/s6-gate-lib): install; tests/cpp from `make clean` 244 ok
(no-change verification); tinytest 5218 results, 0 failures (5186 at S5
+ 32); trio BITWISE, no re-record (37/37, 12/12, 10/10); air clean;
lintr 0 (`.lintr` excludes vignettes/ at directory level, so the new Rmd
is handled exactly as the two shipped ones); NEWS parses (236 entries);
`pkgdown::check_pkgdown` no problems, `- dbarts-embedding` under
Utilities; `R CMD check --as-cran` from a git-archive tarball staged
outside the tree WITH vignettes, Status OK, 0/0/0. Vignette build 5.16 s
(implementer) / 3.4 s (gate-runner) against the 60 s target.

Landed: 55811082, pushed 2026-08-15. Independent gate-runner CONFIRM
(all nine gates; recipe traceability, statistical-prose audit against
the shipped Rd/landing-note record, tolerance recomputation on its own
bitwise-identical run, deviation audits all clean). No feature-matrix
cell bears on the recipes; none updated. The arc-end scoped anchor
refresh stays owed.

## Landing note, S7 (appended 2026-08-15)

**The bound is 1e6, derived, not chosen**: `NBDispersionPrior::computeKernel`
allocates its count histogram as `maxCount + 1` doubles - 8 bytes per unit of
the largest count - so y = 1e9 requested 8 GB where no R error can be raised,
while the bound pins the request at 8,000,008 bytes and the kernel rebuild at
1.3e7 multiply-adds, measured at 0.02 s. The bytes-at-bound arithmetic lives
in the guard comment; the bound is stated in the error message, three Rd
files where the support rule already appears (bart.Rd, dbarts.Rd,
dbartsSampler-class.Rd), `dbarts.h`'s setResponse comment and NEWS. The
Polya-Gamma O(y + r) draw cost (0.63 s at the bound) is recorded as a family
cost, deliberately NOT the refusal's reason. The comparison sits inside
`validateResponseSupport`'s nbinom arm, so one edit covers creation and every
y-swapping conduit on both surfaces - and, by the same channel, S4's
`dbartsDrawLatents`/`dbartsWorkingResponse`, which share the O(y) hang;
intended, no shipped cell moved. The gate-runner traced all six caller sites
at HEAD; the flat surface has NO `dbarts_sampler_setData` (verified), so its
conduit set is {creation, setResponse} and the spec's "both mutation
conduits" reads R-side.

`hostFor` now marks all FOUR bart2 placeholder hosts: the two new
assignments (ordinal, nbinom) mirror the multinomial pair exactly - same
field, same string shape, same placement after the bartcore sampler is
built - and the shipped `refuseHostMutation`/`refuseHostRead` wiring does
the rest, including `$getDispersion` on the count shell, whose placeholder
r is not the fit's.

`refuseBinaryWeightChange` left the R bridge's anonymous namespace for
`bartcore_bridge` - the gate-runner diffed the two bodies byte-identical
(only the comment gained an external-linkage note), confirmed exactly one
definition remains, and read `dbarts_sampler_setWeights` end to end: the
call lands after the multi-forest guard, before the install, and nothing
else changed. Signature deliberately unchanged (no caller label; the
messages name the FAMILY and four shipped cells match on that text).

The `dbarts.h` touch, the slice's sharpest review point: a three-line
comment edit on `dbarts_sampler_setResponse`'s Doxygen block, outside
`DBARTS_C_API_LIST` (which carries only type/name/param tuples, no prose),
so the FNV-1a hash CANNOT see it - and neutrality is GATED, not asserted:
test-capi.R's `0xcd88efcd67de55d7` text-hash cell and its
`matches.header` companion are untouched by the diff and green in both
batteries, and the install's own `static_assert` held. NO apiHash move, no
consumer-lockstep event.

Mutations, each to its own lib, reverted, verified: bound raised to 2e6 -
SIX cells go GREEN (R creation/$setResponse/$setData, flat
creation/setResponse, the refused-swap pin), proving the fixture (1000001,
with 1000000 accepted at-bound) sits just over; ordinal `hostFor` dropped -
3 RED all ordinal, nbinom 79/79 green; nbinom `hostFor` dropped - 4 RED
including `$getDispersion`, ordinal 80/80 green; flat call site dropped -
2 RED both capi-side, every R-side cell green, which is the promotion's
whole point. METHOD LESSON, now in the standing mechanics: restoring a
mutated `.cpp` by `mv` preserves the old mtime, so `make` keeps the mutated
object and the "clean" verification runs mutated code - `touch` after every
revert (one mutation run was re-run after this was caught; the recorded
table is the clean run).

Budget, raw additions against ed04b6f3: non-test 113 as `numstat` counts
it, of which 25 is the setWeights MOVE charged as additions (bodies
byte-identical, verified) - net-new 89 against the ~65 estimate, inside
the 55-95 band; NEWS is 22 for three user-visible refusals. Tests 112
against ~105, inside. No stop approached.

Battery, twice, separate libs (implementer /tmp/s7-lib, gate-runner
/tmp/s7-gate-lib): install; tests/cpp from `make clean` 244 ok; tinytest
5239 results, 0 failures (5218 at S6 + 21); trio BITWISE, no re-record
(37/37, 12/12, 10/10 - all three fixes are refusals no harness scenario
reaches); air clean; lintr 0 on all five touched R files; NEWS parses (239
entries); `pkgdown::check_pkgdown` no problems, no topic owed;
`R CMD check --as-cran` from a git-archive tarball staged outside the
tree, Status OK, 0/0/0. No ASAN leg owed - no new engine code becomes
reachable and the cap makes an existing path unreachable.

Landed: cedf4c34, pushed 2026-08-15. Independent gate-runner CONFIRM (all
eight gates, audits A-F all verified against source). This closes the
response-support-validation entry's whole residue: the G1 landing's
deliberate deviation (the cap), both of its open follow-ons (the flat
setWeights guard, the hostFor extension), with defect 5 closed at S1/S6.
The arc-end scoped anchor refresh stays owed; S8 is the arc's last slice.

## Landing note, S8 (appended 2026-08-15)

**The hash moved by design: 0xcd88efcd67de55d7 -> 0x85bd1ef04beb3848**, for
exactly three appended X-list entries - `dbarts_sampler_dispersion`,
`dbarts_drawLatents`, `dbarts_workingResponse` - with the 44 pre-existing
entries byte-identical (gate-runner diffed the list; nothing re-signed).
`DBARTS_C_API_MAJOR`/`MINOR` stay 1/0; the decision-8 carve-out ("a
pre-1.0-0 append extends the initial field set and moves no version
constant") landed at the two mandated `dbarts.h` sites PLUS
`C_interface.cpp`'s prose lock, which carried the same stale claim, PLUS -
added at the review round - `dbarts_forest_calibration`'s marker, which the
audit found unreconciled: its five S1-landed fields sit BELOW its marker,
`dispersion` sits ABOVE `dbarts_results`', and the marker now states that
both sides of the pre-1.0-0 line are the initial field set, so the next
author has one rule instead of two precedents.

Shipped: `dbarts_results.dispersion` as the struct's LAST field
(offsetof/sizeof lock updated, one FILL line; structSize bounds the write,
so a smaller-header consumer is never written past);
`dbarts_sampler_dispersion` with the R5 getter's channel-test semantics
(returns 0 without touching out off-family); the WRAPPED augmentation
forms, drawing R's stream per-call under the S4 discipline - the raw
`ext_rng` primitives stay internal, the recorded door. The S4 cores moved
to `bartcore_bridge` (five entities, verbatim bodies, single definitions -
audited) so both surfaces run one implementation.

**The ABI hard gate ran TWICE, independently**: implementer and gate-runner
each staged stan4bart (bartcore @ 7ce9a763) by `git archive` into scratch,
compiled against the NEW header (clean), loaded through its
`checkDbartsAPIVersion` hard-equality handshake (accepts the re-baked
token; the implementer also proved the check LIVE - the same binary
against the older user-lib dbarts fails at load), and ran its full suite:
531 results, 0 failures, 21/21 files, both runs. An append-only re-bake
costs zero consumer call sites, exactly as the plan priced; the
four-package re-verification stays where the release checklist put it, at
the freeze.

Flat-path validation, audited adversarially: nothing reaches an OOB access
or unbounded allocation (nbinom y is bounded by S7's cap through the shared
validator; ordinal k-indexing is branch-safe). Three silent-degradation
contracts the R surface validates and the flat surface CONTRACTS are now
STATED in the Doxygen, added at the review round after the audit found them
undisclosed: fractional dispersion is rounded into a different-shape
Polya-Gamma draw, not refused; unordered thresholds corrupt the draw
silently (the rejection loop degrades to its fallback); a non-finite y
under "aft"/"student" propagates NaN into out, per the same contract
`dbarts_sampler_setResponse` states. Known-weak cell, recorded, no action:
consumer.c's `disp$guarded` is TRUE-by-survival on the pre-existing
`DBARTS_RESULTS_HAS`/structSize machinery, which no S8 mutation falsifies.

Mutations: FILL line dropped - exactly 3 RED, all the dispersion-slot leg;
flat draw ignoring offset - exactly 2 RED, the flat-vs-R agreement cell and
the convention-discrimination cell (S4's oracle 3 mirrored); the appended
field misplaced mid-struct - COMPILE FAILURE at
`offsetof(dbarts_results, train)` ("24 == 16"), cascading through every
downstream offset, which is the lock doing its job. Each reverted, touched
(the S7 mtime lesson applied), re-verified 228/0.

Post-audit amendment, one commit kept: the two comment-only header fixes
above, proven comment-only by stripping comments (366 code lines before
and after, IDENTICAL) and hash-gated unmoved (the static_assert compiled
and test-capi.R's 0x85bd1ef04beb3848 cell passed on the amended tree);
install, full tinytest and a fresh clean-copy `R CMD check` re-ran green.

Budget, raw additions against d101ebc9: non-test 350 as numstat counts it,
70 of it the verbatim promotion MOVE - net-own 280 against the 278 stop,
the final +11 being the two ORCHESTRATOR-DIRECTED comment fixes.
ADJUDICATED ACCEPTED: shortening a stated ABI rule to fit a line count
would defeat the fix, and pre-amendment the slice sat at 269, inside.
Tests 187 against the 188 stop. The spec's own comparable (ab3aa2fa)
predicted the weight sits in the capi tests; it does (consumer.c 87,
test-capi.R 100).

Battery, twice, separate libs (implementer /tmp/s8-lib, gate-runner
/tmp/s8-gate-lib): `--preclean` install, zero warnings; tests/cpp from
`make clean` 244 ok plain AND under ASAN/UBSAN, zero diagnostics; tinytest
5256 results, 0 failures (5239 at S7 + 17), test-capi.R 228/0 with
consumer.c COMPILED AND RUN both times; trio BITWISE, no re-record (37/37,
12/12, 10/10 - the augmentation entries call the cores alongside, never
inside, the draw sites); air clean; lintr 0; NEWS parses (242 entries);
pkgdown no problems (flat entries carry no R surface, no topic owed);
`R CMD check --as-cran` from a git-archive tarball staged outside the
tree, Status OK, 0/0/0, re-run on the final amended commit.

Landed: 13393350, pushed 2026-08-15. Independent gate-runner CONFIRM (all
ten gates including its own stan4bart leg; audits A-F). ALL EIGHT SLICES
OF THE ARC ARE LANDED; the arc-end scoped feature-matrix anchor refresh is
the one remaining arc obligation.

## Arc closure (appended 2026-08-15)

The scoped anchor refresh landed at a118ca29: 532 anchors accounted (108
carried over from untouched files, 142 verified unmoved in touched files,
279 moved across 27 DISTINCT deltas - the healthy-sweep signature - and
one dead citation), stamp moved to c05322a8, the S3 header exception
collapsed, the pass recorded in the file's own header. The orchestrator's
sample verification opened seven moved anchors across five delta classes
(including a -16, a +36 and the +1112 relocation) and every one resolved
to its symbol. The one dead citation - the Gaps section still calling
grouped setResponse "an unbuilt door" at the guard S3 deleted - was a
VALUE fix outside the refresh's no-re-adjudication rule, made in the same
commit by the orchestrator: the gap is recorded CLOSED at S3, pointing at
the section-2 cell. CI: docs-only, fires nothing; gate evidence is
built-package byte-identity to 13393350.

THE ADOPTION-SLATE ARC IS COMPLETE: eight slices, every one landed with
an implementer battery plus an independent gate-runner battery on
separate libs, orchestrator diff review, and CI six-green; zero baseline
re-records anywhere, exactly as designed - the trio stayed bitwise at
every slice. Decision 6's slate is discharged in full. Residue living
elsewhere: the freeze-time four-package re-verification (release
checklist), the valgrind/memcheck and rchk legs riding the first
scheduled runs after bartcore reaches main, ticket
host-shell-read-guards, ticket latent-family-weight-channel, the raw
ext_rng door, the grouped re-anchor door (F3 alternative (C)), and the
doors this plan's "Doors held open" section records.
