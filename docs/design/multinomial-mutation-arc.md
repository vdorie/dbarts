# The multinomial mutation arc: scoping and design

Status: LANDED, 2026-08-24, in four slices (S0 5e586587, S1a b2d1749f,
S2+S3 5a3bc276, S4/F1 2619ac9e).

Sections 1 through 4 describe the tree BEFORE the arc landed - the
proposal's inventory, its forks and its prices - and are not statements
about the live code. The live multinomial surface is documented in
`docs/design/multinomial.md` and `docs/design/feature-matrix.md`; the
`retired:` markers below flag individual pre-arc constructs, and the
landing notes at the end record what each slice shipped.

Design text: REVISION 2, 2026-08-20, after independent blind critique
(verdict EXECUTE-WITH-AMENDMENTS). Every one of the critique's 21
findings is folded in with its measured evidence.
Section 7 records the disposition finding by finding; section 8 records
DISPUTED items (none on substance) and four upward refinements where my
re-verification tightened a finding rather than contradicting it.

VD's framing (2026-08-20): pre-release-candidate. The `hostFor` guard is
a stopgap the r-c-division record itself calls temporary
(`docs/design/r-c-division.md:334-341`, "needs a guard or prominent
documentation until the multinomial mutation arc replaces the area").
Guard-all-reads is REJECTED. The arc replaces the scaffolding with a real
public surface, or scopes and names what it refuses.

## Verification method and provenance

Source anchors: read live at tip 0f5a2285, by hand, this pass.

Executed evidence: probes run against the installed dbarts 1.0.0, Built
`2026-08-19 18:15:22 UTC`. **That build is one day behind tip**: the two
2026-08-20 source commits (`7ff77ab4` amplitude rename, `e75a2c79` the
0.9-34 rng/cast port) are not in it - confirmed by
`refuseAmplitudeMutation` being absent from the installed namespace while
present at tip (`R/bartcore.R:36`). Neither commit touches a probed path
(ordinal/nbinom host construction, pointer adoption, multinomial host
family resolution, serialization), and every bitwise probe compares two
arms of the SAME build, so a draw-moving rng port cannot flip a
TRUE/FALSE. Anything labelled MEASURED below was executed; anything
labelled READ was traced in source at tip.

---

## 1. Problem statement and inventory

### 1.1 The four host-shell sites, verified

`hostFor` is assigned at exactly four sites in `R/bart.R` (retired: the
field and its guards were deleted with the mechanism at S4/F1; the site
anchors below are re-derived to the successor constructions):

| site | family | host built by | real engine |
|---|---|---|---|
| `:1473-1481` | multinomial (labels) | `buildHostSamplerCall(family = "multinomial")` (`:1477`) + `samplerCall$data <- y` (retired: the pre-arc call passed `family = NULL` and `as.double(labels)`) | the sampler built here IS the engine that runs (S4/F1); `bartcoreMultinomialSampler` (`R/bartcore.R:932`) is now a thin wrapper over the same public path |
| `:1559-1568` | multinomial (counts) | same, `family = "multinomial"` (`:1563`), `samplerCall$data <- y` (retired: was `as.double(y[, 1L])`) | as above; `bartcoreMultinomialCountSampler` (`:968`) |
| `:1809` | ordinal | `buildHostSamplerCall(family = "ordinal")` | `bartcoreSampler(sampler, family = "ordinal")` (`:645`) |
| `:2058` | nbinom | `buildHostSamplerCall(family = "nbinom")` | `bartcoreSampler(sampler, family = "nbinom")` |

`buildHostSamplerCall` is `R/bart.R:540-560`. `result$bc <- bc` is
retired - the separate `$bc` handle is gone (section 1.5).
`result$fit <- sampler` at `:1502`, `:1582`, `:1882`, `:2124`, gated on
`control@keepTrees || keepSampler`.

MEASURED: two `bartcore_create` calls per fit for ordinal and nbinom -
the R stream advances at both, and the count scales with `n.chains`
(`createChainRngs`, `src/R_interface_bartcore.cpp:1863-1886`, draws one
uniform per chain when `haveSeed` is false).

### 1.2 The asymmetry the defect record flattens

The four sites are NOT one shape. No shipped doc states this.

**Ordinal and nbinom: the host is a REDUNDANT TWIN, not a placeholder.**
`buildHostSamplerCall` passes the real family token, so `dbarts()` builds
a fully correct, fully mutable sampler of the right family - the same
object `dbarts(x, y, family = "ordinal")` returns, whose whole
single-forest mutation surface is SHIPPED
(`docs/design/feature-matrix.md:267-268`, `:297-298`).
`bartcoreSampler(sampler, family = ...)` then builds a SECOND engine from
that same `(control, model, data)` triple and the first is abandoned.
MEASURED: `model@family` on an ordinal host reads `"ordinal"` and the
`bartcore.n.categories` control attribute rides on it - so the triple
genuinely describes the engine. That fact is what makes fork D0 work.

**Multinomial: a genuine placeholder whose resolved family depends on K,
reachable only under an explicit token.** `family = NULL` removes the
token, so `dbarts()` runs `family = "auto"` against `as.double(labels)`,
codes `0..K-1`; `R/spec.R:197-206` resolves auto on a numeric response by
`length(unique(y)) == 2 && sort(unique(y)) == c(0, 1)`.

- K = 2, **under an explicit `family = "multinomial"`**: labels are 0/1,
  so the host is a PROBIT sampler. MEASURED: `control@binary` TRUE,
  `model@family` `"probit"`.
- K >= 3: MEASURED `control@binary` FALSE, `model@family` `"gaussian"`,
  response transform anchored to the range of the integer category codes.
- **`family = "auto"` never reaches this path.** `detectAutoMultinomial`
  (`R/bart.R:1384-1392`) requires `n.levels >= 3L`, and a 2-level factor
  under auto announces probit and returns class `"bart"`. MEASURED:
  `family = "auto": 2-level factor response detected, fitting family =
  "probit"; set 'family' to override`.

`R/bart.R` called the resolved family "irrelevant and never surfaced"
(retired: the comment went with the host shell). It is surfaced - through thirteen unguarded R5 readers.

**Multinomial had no public construction route at all** (retired at
S2+S3, which built one). `dbarts()`'s `family` formal
(`R/dbarts.R:375-389`) now lists thirteen tokens with `"multinomial"`
among them (`:381`); `dbartsSpec()` reaches it the same way
(`R/spec.R:799-808`, `"multinomial"` at `:805`); and creation runs
through the single public dispatch every family uses, `bartcore_create`
(`src/R_interface_bartcore.cpp:3558`), whose multinomial arm is
`createMultinomialDataHolder` (`:3524`) - the dedicated
`C_dbarts_bartcore_createMultinomial` / `...Counts` entries are retired
(Fork J1). What still holds is the absence of a `dbarts.h` creation path
(`feature-matrix.md` `[f4]`): `dbarts_sampler_create`
(`src/C_interface.cpp:542`) routes to `createHolder`, never to the
multinomial arm.

### 1.3 What the `$bc` handles can already do, bridge-side

`$bc` is `new.env(parent = emptyenv())` carrying `$ptr`, `$x`, and (for
multinomial) `$K`. No class, no state, no methods.

Every entry below ALREADY works on a multinomial handle (READ, against
`src/R_interface_bartcore.cpp` and `docs/design/multinomial.md:304-353`).
Bare line anchors are per COLUMN, not per reading order: the bridge column
is `src/R_interface_bartcore.cpp`, the R-wrapper column `R/bartcore.R`
except where it names another file:

| capability | bridge entry | R wrapper | status |
|---|---|---|---|
| create (labels / counts) | `bartcore_create`'s multinomial arm `:3541` -> `createMultinomialDataHolder` `:3507` (retired: the dedicated `bartcore_createMultinomial(Counts)` entries) | `:932`, `:968` | S, unexported |
| run | `bartcore_run` `:4315` | `bartcoreRun` `:1070` | S |
| **response swap** | `bartcore_setCounts` `:3839` | `bartcoreSamplerSetCounts` `:87` | S, unexported |
| **train offset** | `bartcore_setCategoryOffset` `:3914` | `:108` | S, unexported |
| **test offset** | `bartcore_setCategoryTestOffset` `:3953` | `:125` | S, unexported |
| predictors: whole / column / per-obs / joint | `:5167`, `:5119`, `:5280`, `:5334` | | S - no multi-forest guard, by decision (`bart-as-a-component.md:83-90`) |
| cut points | `bartcore_setCutPoints` `:5227` | `:536` | S |
| test predictors | `bartcore_setTestPredictor` `:4796` | `:558` | S (`refuseUndefinedTestFits` `:2900` gates on `testFitsAreDefined`, TRUE here) |
| active-row mask | `bartcore_setActiveRows` `:4074` | retired; `R/dbarts.R:1398` | S, global only (`[f21]`) |
| predict (K-aware, own n x K offset) | `bartcore_predict` `:5840` | `:1085` | S |
| per-category fits / varcounts | `:4147`, `:4280` | retired; `R/dbarts.R:1727`, `:1757` | S |
| calibration read | `bartcore_getCalibration` `:4205` | retired; `R/dbarts.R:1811` | S (map columns, NaN off-map) |
| state store / restore | `:5684`, `:5689` | unresolved | S, STRUCTURAL not bitwise (omega redrawn; `multinomial.md:354-363`) |
| saved trees | `bartcore_getTrees` `:6022` | retired; `R/dbarts.R:2035` | S, forest-indexed |

Refused, with the refusal already written: `setResponse` / `setOffset`
(redirect to counts / category offset, `:2647-2674`), `setWeights`
(`:2668-2672`, same branch), `setSigma`, `setData` / `setModel`
(`refuseMultiForestMutation` `:2618`), `setCalibration` (`:4248-4256`),
`setForestWeights` (`:4034-4037`, permanent, model grounds), `getLatents`
(NULL - `[f22]`), `getFitsWithoutOffset`.

For ordinal and nbinom the handle is an ordinary single-forest sampler:
it can do everything an R5 `dbartsSampler` can do, because it IS the same
engine object an R5 sampler wraps.

**So the engine gap is closed and the surface gap is not.** TODO's
`multinomial-counts-mutation` entry (`TODO:138-144`) then said "The multinomial
family is now a full conditional inside a larger Gibbs/MH sampler." True
of the ENGINE, false of the PUBLIC surface: reaching it needs
`getFromNamespace` against three unexported functions and a classless
environment.

### 1.4 What the R5 class promises that the shells cannot deliver

`dbartsSampler` (`R/dbarts.R:889-2176`): 43 methods, 7 fields (`pointer`,
`control`, `model`, `data`, `state`, `hostFor` (retired: the field was
deleted with the mechanism at S4/F1), `forestWeights` `:902`).

- **22 carry `refuseHostMutation`**: `run`, `sampleTreesFromPrior`,
  `sampleNodeParametersFromPrior`, `growFromRoot`, `setControl`,
  `setModel`, `setData`, `setResponse`, `setOffset`, `setWeights`,
  `setActiveRows`, `setForestWeights`, `setForestBasis`, `setSigma`,
  `setPredictor`, `setCutPoints`, `setTestPredictor`,
  `setTestPredictorAndOffset`, `setTestOffset`, `setCalibration`,
  `setState`, `installTrees`.
- **2 carry `refuseHostRead`** (retired: the guard is gone; both methods
  remain): `getDispersion` (`:1690`), `getFitsWithoutOffset` (`:1737`).
- **13 substantive methods answer from the placeholder, unguarded**:
  `copy`, `predict`, `getLatents`, `getSigmas`,
  `getSumsOfSquaredResiduals`, `getForestFits`, `getForestAmplitudes`,
  `getForestVariableCounts`, `getCalibration`, `storeState`,
  `printTrees`, `getTrees`, `plotTree`.
- Remaining 6 are infrastructure: `initialize`, `refuseHostMutation`,
  `refuseHostRead`, `reapplyForestWeights`, `getPointer`, `show`.
  (22 + 2 + 13 + 6 = 43, checked.)

**`$predict` is silently degenerate, and the real value is worse than a
plausible one.** MEASURED on a K = 3 fit: `f$fit$predict(x)` returns a
CONSTANT 1 - range `[1, 1]`, one distinct value - because the host is
created and never run, so it answers from stumps at the mean of the
category codes `0..2`. Not a plausible-looking surface; a constant, with
no error.

**`$copy()` launders the guard.** MEASURED: original `hostFor` length 1,
copy `hostFor` length 0; `copy$setResponse(...)` ACCEPTED while the
original refuses. `copy` (`:1041`) constructs through
`dbartsSampler$new`, which sets `hostFor <- character(0L)` (retired: both
the field and that reset are gone), and transfers only `state` and
`forestWeights`. The copy is a REAL
independent engine, so the mutations are not no-ops - they take effect on
an engine nothing reads. Same harm, one method call away.

### 1.5 `$fit` is load-bearing, so it cannot simply be deleted

All three `predict` methods code `newdata` against the HOST's design:
`R/generics.R:1344`, `:1650`, `:1927` all call
`validateXTest(newdata, object$fit$data@x)` before handing the matrix to
`bartcorePredict` (retired: routed through `object$fit` directly, not a
separate `$bc` handle). The host carries the factor level
table and column names a data-frame `newdata` is expanded against. Any
fork that removes `$fit` must relocate that table.

### 1.6 The serialization hole, and the absence of any workaround

`$bc` holds an external pointer and nothing else: no state slot, no
class, no re-creation path (contrast `getPointer` `:1815-1851`).

MEASURED, all three families, save then reload then `predict`:

```
multinomial reloaded predict: ERR: bartcore function called on NULL external pointer
  after $fit$storeState():   SAME ERROR
ordinal     reloaded predict: ERR: bartcore function called on NULL external pointer
  after $fit$storeState():   SAME ERROR
nbinom      reloaded predict: ERR: bartcore function called on NULL external pointer
  after $fit$storeState():   SAME ERROR
```

The second line of each pair is the one that matters:
`$fit$storeState()` - the documented `bart` escape at `man/bart.Rd:254` -
does NOT help, because `predict` reaches through `$bc`, which has no
state. **There is no user-side mitigation today for any of the three
families.** No shipped test covers this; no `saveRDS` appears in
`test-multinomial-*.R`, `test-ordinal.R`, or `test-nbinom.R`.

The engine half exists: store/restore work on a multinomial handle
(`test-multinomial-counts-mutation.R:190-204`). What is missing is an
object to hang the state on, and somewhere for the counts to live - they
are data and ride no state block (`multinomial.md:354-363`; that test
file's header, `:7`).

---

## 2. Design forks

### Fork A. Surface shape

**A1** front any handle from `dbartsSampler` (needs the response on
`data`, fork G; `length(data@y)` appears at `:1340`, `:1402`, `:1429`,
`:1477`, and `R/bartcore.R:455`; `$predict`'s offset coercion
`R/dbarts.R:1104-1143` needs a K-matrix arm - S2+S3 built it, the
`!is.null(counts)` branch at `R/dbarts.R:1107-1125`).
**A2** a new public class (doubles the documented surface; breaks the
"an ordinary `dbartsSampler`" promise `man/dbarts.Rd` makes for
multi-forest fits).
**A3** remove `$fit`, fit-object accessors only (kills `predict` per 1.5;
concedes the premise).
**A4 (RECOMMENDED) split by family** - the two shapes are different
problems:
- ordinal, nbinom: the host already IS a correct sampler; the work is
  making `$fit` point at the engine that ran (fork D0), not building
  anything.
- multinomial: A1, via fork G.

Sub-decision, **corrected**: `samplerOnly` is NOT "two lines".
`checkFamilyUnsupportedArgs` (`R/bart.R:618-646`) has **four** callers -
multinomial `:901`, ordinal `:1054`, nbinom `:1087`, and
**`hurdle.lognormal` `:1131`** - and its `samplerOnly` block is three
lines (`:625-627`) inside a helper that refuses two other things.
Deleting the block un-refuses `samplerOnly` for hurdle, whose `$fit` is a
PAIR of samplers (`test-hurdle.R:213-214`) and which this arc does not
touch. RECOMMEND: a per-caller `allow.samplerOnly` flag, TRUE for the
three families this arc fixes, FALSE for hurdle; hurdle's own
`samplerOnly` is Fork I.

### Fork B. Which mutations are meaningful per family

**Ordinal, nbinom.** Nothing owed. Every channel is built and `S`;
thresholds ride `run()$thresholds`, `r` rides `run()$dispersion` and
`$getDispersion()`. MEASURED: on a pointer-adopted nbinom sampler
`$getDispersion()` answers with the FIT's `r`, which is precisely the
read `refuseHostRead` blocks today.

**Multinomial.**

| channel | verdict | reason |
|---|---|---|
| `$setCounts` | BUILD | the response; n and K fixed at creation |
| `$setCategoryOffset` | BUILD | n x K, before the softmax; only the row-centred part identified |
| `$setCategoryTestOffset` | BUILD | nTest x K, reported test probabilities only |
| `$setPredictor` (4 shapes), `$setCutPoints`, `$setTestPredictor` | ALREADY OPEN | |
| `$setActiveRows` | ALREADY OPEN, GLOBAL | per-forest masking permanently refused: the margin is a log-sum-exp over the other K-1 forests |
| `$setResponse`, `$setOffset` (flat) | REFUSE | a flat shift is the softmax's null direction |
| `$setWeights` | REFUSE | no exact real-shape PG draw; an integer weight IS count replication |
| `$setSigma` | REFUSE | no sigma |
| `$setData`, `$setModel` | REFUSE | rebuild forest 0 only; survey doors 1/3 undesigned |
| `$setCalibration` | REFUSE | the softmax map owns every leaf scale |
| `$setForestWeights`, `$setForestBasis` | REFUSE | no amplitudes; per-forest weighting incoherent |
| `$getFitsWithoutOffset` | REFUSE | no single additive location; `$predict(data@x)` serves it |

**B1 (OPEN VD FORK): `$getLatents` on a multinomial sampler.** NULL
today (`[f22]`). r-c-division clause 4 defaults to read/write symmetry -
"whatever the engine draws, the R program may read" - and the K x n omega
matrix IS drawn every sweep (`combiner.hpp:1428`). (i) build a
combiner-side `latents()` (engine touch: `tests/cpp` cell plus an
ASAN/UBSAN leg; ~40 engine + ~60 test by the
`latent-family-weight-channel` sibling estimate); (ii) record a
considered DECLINE (the omegas are interleaved one-vs-rest augmentation
variables refreshed against a moving margin, so a host reading them
between sweeps reads a quantity whose conditioning it cannot reproduce;
the composition recipes reach the same posterior through `f + o` without
ever touching omega); (iii) leave NULL and say nothing - the accident
clause 4 forbids. RECOMMEND (ii), written into
`man/dbartsSampler-class.Rd` and `multinomial.md`. Genuinely VD's call:
(i) is cheap and clause 4 leans toward it.
DECIDED (VD 2026-08-20): (ii), the written decline. Usefulness is the
criterion; latents are not persisted absent a compelling use case.
Scope the decline's rationale to the AUGMENTATION families (PG-type
omega/precision channels, meaningless marginally) - do not phrase it
as covering latents globally, so AFT's documented imputed log-time
read stays coherent.

**B2: `$setResponse` sugar for counts? NO.** The shapes differ (n vector
vs n x K matrix) and silently reinterpreting a length-n integer vector as
one-hot counts is the kind of guess the repo refuses elsewhere. Refuse
and name `$setCounts`, with the reason in the message.

### Fork C. Does the flat C ABI grow?

Facts, corrected and re-anchored: `DBARTS_C_API_HASH` is
`0x66d33f1613892406ULL` (`inst/include/dbarts/dbarts.h:189`);
`DBARTS_C_API_MAJOR`/`MINOR` are 1/0 and the header itself states the
pre-release carve-out - "no version of this API has shipped yet, so
whatever they read then simply IS the initial contract, and they do not
move before it" (`dbarts.h:135-139`, which is the citation to use, not
the plan doc). Multinomial has no flat creation path
(`feature-matrix.md` `[f4]`). One reverse `LinkingTo`, in-house.

**Cost of C2, corrected:** TWO literal re-bakes, not one -
`DBARTS_C_API_HASH` (`dbarts.h:189`) AND `dbarts_apiSignatureToken`
(`src/C_interface.cpp:513`, `static_assert(... == 0xcb83367ee0c4175bULL)`),
because appending to `DBARTS_C_API_DECLS` moves the signature half too.
Plus `inst/tinytest/test-capi.R`: `:84` pins the literal and must
change, and one `expect_false` for the superseded literal should be added
alongside the ones already there. The quoted pin has since moved: the
header re-baked when `dbarts_sampler_predict` took a thread count, so
`:84` now reads `expect_identical(hashes$text, "0x66d33f1613892406")` and
`0x6c9776ae1197e8f5` is itself one of the superseded literals, at `:81`
beside `:73-74` and `:78`.

**C1** no growth / **C2** grow / **C3** reserve only (~30 lines in
`docs/plans/archive/c-api-growth.md`).

RECOMMEND **C3 now, C2 as a post-RC door.** The critique's finding 21 is
right and my earlier C2-now recommendation rested on a misidentification
(section 7, finding 11): I inferred "a second header window already
landed after adoption-slate S8, so windows are cheap" from
`0x85bd1ef04beb3848`, which is in fact the LIVE signature-half token, not
a superseded hash. The conclusion survives on different evidence -
`test-capi.R:73-74` names `0x1a911c00bb26dcd7` and `0xcd88efcd67de55d7`
as superseded header literals, so the hash HAS re-baked more than once -
but the inference I built the recommendation on was wrong, and the
remaining reasons point the other way: the enabling consumer is
stan4bart, which is on its own branch with its own release, so nothing is
blocked by deferring. Reserve now; spend the window when stan4bart asks.

### Fork D. Construction

**D0 (NEW, RECOMMENDED). Pointer adoption for ordinal and nbinom.** The
critique's finding 1, independently re-verified here. `$fit` does not
have to be the FIRST-created engine; it has to be an R5 object pointing
at the engine that RAN. `dbartsSampler$pointer` is a plain field. Both
creates still happen, in the same order, so **the R stream is untouched
and no draw moves**.

MEASURED at tip-minus-one-day, ordinal and nbinom arms:

```
adopted R5 run == twin run, bitwise train:       TRUE  (both families)
adopted R5 run == twin run, bitwise cutpoints:   TRUE  (ordinal)
adopted R5 run == twin run, bitwise dispersion:  TRUE  (nbinom)
adopted sampler $setResponse:                    OK
adopted sampler $setPredictor:                   OK
adopted $storeState() / $getPointer() round trip: OK
adopted $state non-null:                          TRUE
adopted $getDispersion() (nbinom):                5   (the fit's r)
host model@family:  "ordinal";  control attrs: bartcore.n.categories
```

Soundness, traced at tip: the abandoned host externalptr becomes
unreachable and its finalizer deletes its own holder once
(`holderFinalizer`, `:88-94`) - no double free, since each externalptr
carries its own holder; the adopted ptr's protection slot already pins
`sampler$data`, the same S4 object the R5 mirrors into, so every setter's
retained vector lands on the right holder; the `state` promise is
`delayedAssign`ed against the object's own `pointer`
(`R/dbarts.R:932-946`), so it resolves against the adopted engine; and
`getPointer()`'s re-creation branch rebuilds from
`(control, model, data)` with `model@family == "ordinal"` and the
`bartcore.*` control attributes intact, which is exactly correct for
these two families.

**Adoption is UNSOUND for multinomial**, and for the same reason: `data@y`
is label codes, so `getPointer()` would silently re-create a
gaussian (K >= 3) or probit (K = 2) engine. That is precisely the split
Fork A4 already argues for.

Two refinements on the mechanism: (a) express adoption as a named R5
method or a construction argument, not a bare `$pointer <-` poke from
`R/bart.R`, so the contract is documented and testable; (b) adoption
leaves the wasted create in place - the "two engines per fit" fact
SURVIVES into the RC. Do not let "the win is delivered" imply otherwise.

**D1. Delete the twin (post-RC door).** One create instead of two.
MOVES DRAWS - but only for `seed = NA`. MEASURED: with `seed = 77L` the
twin's run and the host's run are `identical`; without a seed they
differ. So D1 is bitwise for every seeded user, and the 2-of-43
re-record is a property of `fitSummaries`' `set.seed()` fixture
(`benchmarks/R/equivalence.R:1621-1622`), not of the change. The prize is
small: MEASURED at n = 5000, p = 30, 100+100 sweeps, twin create 0.003s
against a 0.339s run - about 0.8% of a short fit, not scaling with
`n.samples`. **Do not spend the program's first re-record to reclaim
3ms.** Defer.

**D2 (multinomial). Keep host-then-convert, stop returning the host.**
Cheapest; leaves two engines, an orphaned level table, and save/load
broken.

**D3 (multinomial, RECOMMENDED). Direct**, via fork G. One create
instead of two. MOVES DRAWS for `bart2(family = "multinomial")`.

D3 does not have to move `benchmarks/R/multinomial-equivalence.R`: that
harness reseeds between the host build and the handle build (`:158-165`,
`:177-184`) and `bartcoreMultinomialSampler` makes exactly one `.Call`
(`R/bartcore.R:928-933`), so the 11 baselines are insulated - see Fork J
for what that implies and how to discharge it.

What D3 breaks is `test-multinomial-surface.R`'s REPRODUCTION GATE
(`:1-11`, `:97-102`), whose comparator `internalMultinomialFit`
(`:12-40`) mirrors the host-then-bc sequence "with no reseed in between".
Both sides change together, which is why the gate must be re-pointed
rather than trusted. **And it is the only instrument that path has** -
see finding 9 / section 3's S0, which adds a real baseline first.

### Fork E. predict / keepTrees / state-save

**E1** status quo. **E2 (RECOMMENDED)** mirror the `bart` contract:
`$fit` IS the sampler; `predict.bart{Multinomial,Ordinal,Negbin}` routes
through `$fit$predict(newdata, ...)`; `$bc` deleted; save/load requires
`fitObj$fit$storeState()`, which `getPointer()`'s re-creation branch
already implements. This is the arc's largest user-visible win and today
there is no workaround at any price (1.6). Hard prerequisite for the
multinomial third: fork G1.

`$bc` deletion detail the S4 spec must carry: `$bc` **co-gates**
`thresholds.raw` (`R/bart.R:1879`) and `dispersion.raw` (`:2121`) inside
the same `keepTrees` block, and `predict.bartOrdinal` /
`predict.bartNegbin` read them. Those two channels KEEP their gate.

### Fork F. The `hostFor` mechanism

After A4 + D0 + D3 there are zero `hostFor` assignment sites.
**F1 (RECOMMENDED) delete** the field, both guards, 24 call sites, and
the three Rd sentences (retired: deleted at S4/F1);
`host-shell-read-guards` closes as OBVIATED. Do it in the LAST
behavioral slice - the guard is load-bearing until the last shell dies.
**F2** keep as a general mechanism: a guard with no caller is dead code
the RC cleanup waves flag; re-adding 12 lines later is trivial.

**Test footprint, corrected: SEVEN files, not two and not five**
(verified live): `test-dispersion-channel.R` (retired: the shell cell was
deleted, not rewritten), `test-fits-without-offset.R:231-241`,
`test-calibration-midchain.R:466-469`, `test-forest-weights-r5.R:232-235`,
`test-ordinal.R:145-152`, `test-nbinom.R:143-152`,
`test-multinomial-surface.R:346-360`. The
ordinal and nbinom blocks pin the FAMILY-SPECIFIC message text
(`"host sampler of a bart2\\(family = \"ordinal\"\\) fit"`), and under D0
those cells INVERT - the sampler becomes mutable - so they must be
rewritten as capability assertions, not deleted.

That footprint is spent. `hostFor` has zero survivors in `R/`, `src/`,
`inst/` and `man/` - the sole remaining occurrence of the name is the
descriptive comment at `inst/tinytest/test-host-shell-pins.R:4`. The
ordinal and nbinom cells did invert to capability assertions
(`test-ordinal.R:145-152`, `test-nbinom.R:143-152`, both now asserting
that the retained `$fit` reads, mutates and runs); the
`test-dispersion-channel.R` cell was deleted rather than rewritten; and
the four multinomial cells survive as refusals stated on MODEL grounds
rather than host-shell ones.

### Fork G. Where the multinomial response lives

**G1 (RECOMMENDED). `dbartsData` gains slots**: `counts` (n x K
integer), `offset.category`, `offset.category.test`. The `bases`
precedent is exact (`R/A_class.R:513-524`, whose own comment states the
rule). Re-creation becomes free; save/load works; `$setCounts` mirrors as
`$setWeights` does.

Three costs the first draft missed, all confirmed:

- **The `dbartsData` validity method is in scope.** `R/A_class.R:565-626`
  derives `numObservations <- length(object@y)` and carries the refusals
  for `x`, `varTypes`, `weights`, `offset`, `bases`, `x.test`. Three new
  slots must be constrained there. Priced into S2 now; it was not before.
- **Two family gates are required, and one cannot key on family.**
  MEASURED on `dbartsData(x, rep(1, n))`: the degeneracy check
  (`R/data.R:1699-1710`) fires `response values are indistinguishable, or
  nearly so, at double precision (1 distinct value(s) among 60
  observations)`, and `dbarts(x, rep(1, n))` installs a starting sigma at
  `estimateSigmaFromLinearModel`'s host-independent floor (was measured at
  `1.52e-16` before the floor landed). A single-trial multinomial has
  `n_i == 1` for every row -
  the dominant case - so both fire on EVERY such fit. The sigma gate is
  one condition (`fixedUnitScale`, `R/spec.R:218-225`, already a
  family-keyed flag). The warning gate is NOT: `dbartsData()` has no
  `family` formal (`R/data.R:1699-1710`), so it must key on the data object
  - skip the range check when `@counts` is non-NULL, which is available
  at construction and is the honest predicate anyway.
- **Option (ii) is a hard error, not "five sites".** MEASURED:
  `dbarts(x[0, ], numeric(0))` stops with `data has zero rows` at the R
  layer, and the bridge would stop with `length of y must be greater than
  0` (`src/R_interface_bartcore.cpp:945-946`) if it got there. Option (i)
  still wins, on corrected costs.

**G1a: `data@y` carries the TRIALS vector `n_i`.** Length n, meaningful,
and no `length(data@y)` check moves. (Option (ii) `numeric(0)` is the
hard error above; option (iii) the argmax code is a lie.)

**G2** an R5 field reapplied at re-creation (the `forestWeights` pattern)
repeats a shape `bcf.md:164-174` already books as a defect - "the two
stored states compare equal while the fits diverge".
**G3** creation-argument only makes E2 impossible.

**G-migration (NEW SUB-FORK, from finding 5): adding slots invalidates
every previously serialized `dbartsData`.** My earlier 4.1 said the
opposite; it was true of the engine blob and false of the R object.
MEASURED: add a slot to a class, `readRDS` an old instance, read the new
slot -> `no slot of name "b" for this object of class ...`. Every saved
`bart`/`bart2`/`dbartsSampler` holds a `dbartsData`, so the first read of
`@counts` errors - and this lands on bartCause / treatSens / stan4bart /
bairrtt users holding saved fits.
- **(m1) Plain break, enumerated in NEWS.** Pre-release, so it is a
  migration cost, not a compatibility bar - but it must be a NEWS line,
  not a sentence about `minReadableStateFormatVersion`.
- **(m2) Guarded accessor.** MEASURED: `methods::.hasSlot(z, "b")`
  returns FALSE on an old instance without erroring, so routing every
  internal read of a new slot through it keeps old fits readable.
  (`methods::updateObject` does NOT exist in base R - verified - so that
  route is unavailable.) ~10-15 lines plus a round-trip test.
RECOMMEND **m2 plus the m1 NEWS line**: 15 lines buys four downstream
packages a clean upgrade, and the NEWS entry is owed either way.

### Fork H. Who names the channel in a cross-surface refusal

`refuseMultiForestResponseMutation` (`:2636-2681`) is shared with the
flat C API and its own comment forbids stating two rules (`:2650-2651`).
Its multinomial arm then told R users to call `bartcore_setCounts` -
an internal C name not callable from R.
**H1** keep C names (status quo, wrong for R users). **H2** thread a
surface label (against the shared-guard note in letter). **H3
(RECOMMENDED)** name the CAPABILITY, not the entry - "replace it through
the counts channel" - one string, right on both surfaces, and what
`docs/design/error-style.md` (ADOPTED 2026-08-17) asks for. H3 LANDED at
S2+S3: that is the message the response arm now carries (`src/R_interface_bartcore.cpp:2676-2677`).

Same edit, cheap: with `supportsCountsMutation` true the WEIGHTS conduit
then fell through to the RESPONSE message (`:2652-2667` tested only
`conduit == offset`), so `bartcore_setWeights` on a multinomial sampler
answered "this sampler's response is its n x K count matrix". Give
weights its own arm naming the model reason. LANDED with H3: the weights
arm is `:2661-2665`, and it states the model reason - an integer weight
already IS row-wise count replication, and a non-integer one has no exact
augmentation sampler.

---

## 3. Pricing and slice plan

Raw additions only; the 1.5x stop applies to these numbers. Comparables
measured live (the critique independently confirmed all five totals):

| commit | subject | non-test | test |
|---|---|---|---|
| `05ac3b4b` | public `setForestWeights` + mirroring on re-creation | 86 | 203 |
| `50d540da` | open the public multinomial offset surface | 176 | 182 |
| `ca820eb0` | open the multinomial counts mutation channel | 239 | 568 |
| `60d7eb7c` | create BCF samplers through the public spec surface | 712 | 385 |
| `d809b944` | read and rewrite the leaf calibration mid-chain | 568 | 652 |

**Four-slice floor** (finding 21). Findings 1 and 11 free two slices:
D0 delivers ordinal/nbinom without a re-record, and Fork C reserves
rather than spends a header window.

| slice | scope | non-test | stop | test | stop | draws |
|---|---|---|---|---|---|---|
| **S0** | Pins + the two things that must not wait. (a) The 22/2/13/6 census; the `$copy()` laundering; the constant-1 host `predict`; the K-dependent host family under an explicit token; today's save/reload failure AND the absent `$fit$storeState()` workaround, all three families. (b) The shipped **Rd defect** (below) - it contradicts a shipped test, so it does not wait for a docs slice. (c) A **`bart2(family = "multinomial")` scenario appended to `equivalence.R`**, so the family bart2 emits has a baseline before S4 moves it. (d) The `$copy()` one-line `hostFor` transfer. (e) Fork C3's reserve doc. | 120 | 180 | 260 | 390 | 42/42 existing scenarios bitwise; one appended row; new baseline hash |
| **S1a** | **Pointer adoption** (Fork D0): a named R5 adoption method, applied at the ordinal and nbinom sites; `$fit` becomes the engine that ran; `predict` routes through it; two `hostFor` sites die; the two family-specific refusal blocks invert to capability assertions; `samplerOnly` unblocks for these two via the per-caller flag (Fork A sub-decision). | 190 | 285 | 260 | 390 | **NEUTRAL - 43/43, 12/12, 11/11 bitwise** |
| **S2+S3** | The multinomial surface. `dbartsData` slots + **validity method** + the two family gates + the `.hasSlot` migration guard; `bartcore_create` dispatch arm; `dbarts(family = "multinomial")`; the three R5 methods; the K-matrix `$predict` arm; the refusal matrix (Fork B); the message repoint (Fork H). | 1060 | 1590 | 830 | 1245 | NEUTRAL |
| **S4/F1** | `bart2(family = "multinomial")` constructs direct; `$fit` is the sampler; `$bc` deleted (4 sites, 3 generics, 3 test assertions, 3 Rd paragraphs, 1 NEWS sentence) while `thresholds.raw`/`dispersion.raw` keep their gate; `hostFor` mechanism deleted across seven test files; reproduction gate re-pointed; Fork J discharged. | 330 | 495 | 300 | 450 | **MOVES** - re-records exactly the S0 `bart2multinom` scenario, by design |

Arc total: **1700 non-test / 1650 test**.

**S2+S3 repricing, per finding 12.** My first draft priced S2 at 420
against a named comparable of 712 and was 1.7x light by its own
yardstick. S2 is `60d7eb7c`'s shape (bridge + `spec.R` + `data.R` +
`A_class` + Rd) PLUS three slots, a matrix-valued response ingestion, the
validity work, two family gates and the migration guard: **~700**. S3 is
three R5 methods against `05ac3b4b`'s 86/method unit (258) with the
K-matrix `predict` arm, refusal matrix, message repoint and Rd on top:
**~360** non-test, not 260. Test: `05ac3b4b`'s 203/method unit gives 609
for three; shared fixtures pull it to ~500, plus ~330 for the data-layer
half. Merged as one slice at 1060/830 they sit under one stop; if the
implementer splits them, S2 = 700/330 and S3 = 360/500, each with its own.

### Deferred residue: post-RC doors, with prices

| door | content | price | cost of deferring |
|---|---|---|---|
| **Door 1: twin-create deletion** | Drop the redundant `bartcore_create` for ordinal/nbinom (Fork D1). | ~60 / ~40, **moves 2 of 43** equivalence scenarios | 3ms per fit (MEASURED 0.003s against a 0.339s run) and "two engines per fit" survives into the RC |
| **Door 2: flat ABI (Fork C2)** | `dbarts_sampler_createMultinomial` + counts/offset setters + a K-aware predict. | ~260 / ~200, **two literal re-bakes** (`dbarts.h:165`, `C_interface.cpp:484`), `test-capi.R:84` + one new `expect_false`, one in-house consumer rebuild | multinomial stays R-only; the named consumer (stan4bart) is on its own branch and its own release, so nothing is blocked |
| **Door 3: nbinom loop collapse** | Replace the per-sweep loop (`R/bart.R:2078-2103`) with one `run`. | ~40 / ~30; **`bench-sampler` IS owed** - it removes `n.samples` R-level round trips from the sampling path, which is sweep-proportional | a per-sweep R round trip per kept draw |
| **Door 4: `$getLatents` build (B1(i))** | Combiner-side `latents()` returning the n x K omegas. | ~40 engine / ~60 test, plus an ASAN/UBSAN leg | clause 4's read/write symmetry stays deliberately declined rather than met |
| **Door 5: hurdle `samplerOnly`** | Decide whether a two-sampler `$fit` may be returned unrun. | decision, then ~20 / ~30 | the per-caller flag keeps today's refusal, which is correct-by-default |

**Door 3 evidence (finding 14, re-verified and extended).** MEASURED, a
seeded control, four cells: `run(3, 4)` batched against a loop of
`run(3|0, 1)` is bitwise on BOTH `$train` and `$dispersion` at
`keepTrees = FALSE` and `TRUE`, and at `n.thin = 1` and `2`. The saved
slot base advances across calls (`Sampler::run` sets
`setSavedSlotBase(currentSampleNum_)`, `sampler.hpp:302`, advanced at
`:459`), so `keepTrees` is not a counterexample. The collapse is safe
whenever VD wants it.

### Gates, per slice

Standing battery: `R CMD INSTALL . --preclean` into a PRIVATE lib
(mandatory on S2+S3, which touch the bridge); `tests/cpp` from clean;
full tinytest FAILURES == 0; `air format --check` + `lintr` on every
touched R file; `pkgdown::check_pkgdown`; full `R CMD check --as-cran`
from a clean-copy tarball; rchk on the next scheduled run for S2+S3.
ASAN/UBSAN owed only where new engine code becomes reachable - S2's
dispatch arm, and Door 4 if taken; NOT owed on S0/S1a/S4.

Equivalence trio, **42/12/11 today** (verified against
`equivalence-d15a2bfb.rds`, `bcf-equivalence-6e3b9fb8.rds`,
`multinomial-equivalence-4d9a3337.rds`), **43/12/11 after S0**:

- **S0**: all 42 existing scenarios bitwise; one appended
  `bart2multinom` row; the baseline file re-records because it grows.
  The new scenario runs LAST with its own `set.seed`, the harness's
  established isolation pattern (`multinomial-equivalence.R:189-196`), so
  it perturbs nothing above it.
- **S1a**: **43/43, 12/12, 11/11 BITWISE, no re-record.** Both creates
  still happen in the same order. This is the whole point of D0.
- **S2+S3**: 43/43, 12/12, 11/11 bitwise - the dispatch arm is
  unreachable from any shipped path until S4.
- **S4**: **re-records exactly one scenario** - S0's `bart2multinom` -
  with a recorded reason (one create instead of two). 42 others bitwise,
  12/12, 11/11 bitwise.

`bench-sampler`: NOT owed for S0/S1a/S2+S3/S4 (no hot-path change; the
create savings are construction-time and are deferred to Door 1 anyway).
**Owed for Door 3** if taken.

**SBC: nothing owed, and this closes here rather than deferring.**
`benchmarks/R/sbc.R` contains no `bart2(` call anywhere, and
`sbcMakeMultinomial` (`:1599-1602`) goes through
`getFromNamespace("bartcoreMultinomialSampler")`. Every SBC arm drives
the handle, so no arm sees S1a's or S4's draw changes.

---

## 4. Risks and interactions

**4.1 State format: no bump owed, but the R object breaks.** The engine
half is clean - the counts are DATA, not state
(`multinomial.md:354-363`), G1 puts them where re-creation already reads
from, and `minReadableStateFormatVersion` is untouched. The R half is
not: adding S4 slots invalidates every previously serialized
`dbartsData` (Fork G-migration, MEASURED). Take m2 (the `.hasSlot`
guard) and write the m1 NEWS line regardless. Two further notes for the
test author: a multinomial state round trip is **STRUCTURAL, not
bitwise** (omega is redrawn), so a bitwise cell will fail; and the
existing format floor keeps every shipped state loading.

**4.2 Equivalence baselines.** Section 3. The non-obvious one is Fork J
(below): the 11 multinomial baselines are insulated by a FIXTURE
property, not by anything in the code under change.

**4.3 Snapshot re-derivation.** S1a moves no draws, so `test-ordinal.R`
and `test-nbinom.R` need only their refusal blocks inverted. S4 moves
`test-multinomial-surface.R` and any `test-multinomial-*.R` cell running
through `bart2`; regenerate by REPLAYING THE WHOLE FILE - the values
depend on the file's full execution history, not just the preceding
`set.seed`.

**4.4 Docs and Rd surface, plus one live shipped defect.**

- **DEFECT, to S0 not S5 (retired: fixed by removing `$bc` rather than by
  correcting the sentence).** `man/bart2.Rd` at the time said `bc` and
  `fit` are "present only under `keepTrees`", false for `$fit` since all
  four gates were `control@keepTrees || keepSampler`. Current
  `man/bart2.Rd` (e.g. line 332 for ordinal, 336 for nbinom) states only
  `fit`'s gate; there is no `bc` left to misdescribe.
- `man/bart2.Rd:301` (the multinomial refusal set).
- `man/dbartsSampler-class.Rd`: the host-shell sentences (retired: deleted
  under F1), the `getLatents`
  paragraph (Fork B1), `\item{forest}`, and the fit-surface table
  `adoption-slate.md:224-233` added.
- `man/bart.Rd:254` (save/load) extended to the three families - and it
  is currently the ONLY place a user is told a workaround exists, which
  MEASURED does not work for them.
- `man/dbarts.Rd`: the `family` formal gains `"multinomial"`; the
  counts-matrix response documented.
- `docs/design/feature-matrix.md`: the `multinom` row of section 2, the
  `getLatents` and calibration cells of section 3, footnotes `[f4]`,
  `[f11]`, `[f21]`, `[f22]`, `[f23]`, and the gap paragraph at
  `:961-972`.
- `docs/design/r-c-division.md` defect 6 (`:334-341`) closes;
  `docs/design/multinomial.md:190-380`;
  `docs/design/bart-as-a-component.md` section 2.
- `TODO`: `host-shell-read-guards` (retired: the entry is gone) closes as OBVIATED;
  `multinomial-counts-mutation` (`:138-144`) gains the surface note.
- `inst/NEWS.Rd`: the `$fit` change, the `$bc` deletion (there is a
  shipped sentence at `:2034`), `dbarts(family = "multinomial")`, the
  three methods, and the serialized-`dbartsData` migration.
- `docs/plans/archive/c-api-growth.md`: Fork C3's reserve.
- `docs/design/INDEX.md` if a new design doc lands (47 docs besides the
  index today - verified).

**4.5 The two `static_assert` residue strings (retired: swept by the
amplitude rename; both now read "an amplitude coupling is a constant-leaf
model").** `combiner.hpp:727` and `chain.hpp:734` no longer say BCF;
`TODO:54-57` said sweep opportunistically with the
next edit to those files. Honest read: **no slice above is expected to
edit either file** - the multinomial engine is done. If Door 4 is taken,
or if S2's dispatch arm turns out to need a combiner touch, sweep both in
that slice and say so in the landing note. Do NOT schedule an edit purely
to discharge the residue.

**4.6 Cross-surface refusal wording.** Fork H. New messages must be
written to `error-style.md`'s R1-R14 from the start; the repo-wide sweep
is release-candidate-review slice L and this arc should not add to it.

**4.7 Interim exposure.** Until S1a/S4 land, thirteen readers answer from
a placeholder, `$predict` returns a constant, and `$copy()` launders the
guard. S0's pins record that behavior and S0's one-line `hostFor`
transfer in `copy` closes the laundering - that is a hole in the guard
that already ships, NOT the rejected guard-all-reads stopgap.

**4.8 Consumers (retired: `$bc` was removed with no replacement field;
`inst/NEWS.Rd:2034` records the deletion).** `$bc` was a bare
environment named in `man/bart2.Rd`'s Value. Full in-repo footprint at the
time this section was written: six code readers (`R/generics.R:659/649`,
`794/810`, `924/939`), three test assertions (`test-ordinal.R:156`,
`test-nbinom.R:142`, `test-multinomial-surface.R:354`), three Rd
paragraphs, one shipped NEWS sentence (`inst/NEWS.Rd:2034`), and the
`thresholds.raw`/`dispersion.raw` co-gate. Under the standing
no-backwards-compat constraint this was a NEWS migration, not a design
input, and it has since landed.

---

## 5. Claim ledger

Doc-sourced claims, verified live at tip. Items 1-16 carried from
revision 1 with corrections; 17-22 added this pass.

Sections 5-8 record dispositions verified against the tree as of the
2026-08-20 blind critique; their VERIFIED/STALE verdicts and quoted
anchors are period-correct within these four sections and are not
re-resynced here (they are not live pointers into the current tree).

1. `r-c-division.md:339-341` - "whether to extend the `hostFor` guard to
   [ordinal and nbinom] is a VD question BY NAME, unscheduled."
   **STALE.** Extended at adoption-slate S7 (`cedf4c34`); `R/bart.R:1754`
   and `:1999` carry the assignments. (Defect 6 spans `:334-341`; my
   revision-1 citation `:337-341` clipped it.)
2. `TODO:150-168` - "four shipped readers ... answer from the shell."
   **TRUE BUT UNDERSTATED.** Thirteen substantive unguarded methods;
   `$predict` among them, MEASURED to return a constant; `$copy()` a
   further hole.
3. `adoption-slate.md:212-215` (**`docs/plans/`, not `docs/design/`**) -
   "the host shell, a gaussian sampler". **PREMISE INCOMPLETE,
   CONCLUSION SURVIVES.** Gaussian at K >= 3; PROBIT at K = 2 - and K = 2
   is reachable only under an **explicit `family = "multinomial"`**,
   since auto announces probit and returns class `"bart"` (MEASURED).
   `numReportedLocations` is 1 either way.
4. `multiforest-extension-surface.md:906-908` (**`docs/plans/`**) - "the
   R5 fields are exactly six ... (`R/dbarts.R:711-722`)". **STALE, both
   halves.** Seven fields at `R/dbarts.R:848-866`.
5. `feature-matrix.md` `[f11]` - "Every R5 mutator is refused by
   `refuseHostMutation` (`dbarts.R:914-925`)". **VERIFIED.** 22 sites.
6. `feature-matrix.md` `[f11]` anchors `bart.R:1451, 1534`. **STALE.**
   Live `:1462`, `:1539`; `hostFor` at `:1443`, `:1526`.
7. `feature-matrix.md` `[f4]` anchors `bartcore.R:846, 909`. **STALE.**
   Live `:852`, `:915`.
8. `adoption-slate.md:248-249` (**`docs/plans/`**) anchors four
   `R/bartcore.R` wrappers. **VERIFIED** at `:1093`, `:1329`, `:1339`,
   `:831`.
9. `bart-as-a-component.md` anchors. **FIVE VERIFIED EXACT** -
   `refuseMultiForestMutation` `:2629`, `refuseMultiForestResponse
   Mutation` `:2652`, `refuseUndefinedTestFits` `:2849`,
   `Chain::supportsResponseMutation` `chain.hpp:1068`,
   `AmplitudeForestCombiner` `combiner.hpp:726`. The sixth
   (`refuseAmplitudeMutation`, `R/bartcore.R:36`) is exact as an anchor
   but is cited at **doc line 74**, outside the `:42-71` window I quoted;
   "ALL SIX VERIFIED EXACT" overstated the window, not the anchors.
10. `TODO:309-322` - "a full conditional inside a larger Gibbs/MH
    sampler." **TRUE OF THE ENGINE, FALSE OF THE PUBLIC SURFACE.**
11. `r-c-division.md:28-31` - "the response side is sealed."
    **STALE BY LANDING** at the bridge; still true at the surface.
12. `feature-matrix.md:999-1005` - the multinomial gap list.
    **VERIFIED, every clause.** The most accurate single statement of the
    problem in the repo.
13. `adoption-slate.md:307` (**`docs/plans/`**) - "the shell cell."
    **VERIFIED**: `test-fits-without-offset.R:239-242` and
    `test-dispersion-channel.R:271-272`. Both lose their subject under F1.
14. **RETRACTED AND REPLACED.** Revision 1 claimed
    `0x85bd1ef04beb3848` was a superseded `DBARTS_C_API_HASH`, proving a
    second header window had landed after adoption-slate S8. **Wrong.**
    That literal is the LIVE `dbarts_apiSignatureToken`, static_asserted
    at `src/C_interface.cpp:460`, and `test-capi.R:78` asserts the full
    hash is NOT it. The conclusion survives on other evidence -
    `test-capi.R:72-73` names `0x1a911c00bb26dcd7` and
    `0xcd88efcd67de55d7` as superseded header literals - but the
    inference was unsound and Fork C's recommendation changed with it.
    Also: the MAJOR/MINOR carve-out should be cited to
    **`dbarts.h:100-104`**, the header's own words, not to a plan; the
    carve-out text sits at `adoption-slate.md:1757-1759` and decision 8
    is owned by `docs/plans/archive/dbarts-h-reshape.md`.
15. `man/bart2.Rd:291/:295/:299` - "`bc` and `fit` present only under
    `keepTrees`". **WRONG for `$fit`**, and contradicted by
    `test-multinomial-surface.R:340-354`. Moved to S0.
16. `docs/design/INDEX.md:3` - "47 docs besides this index".
    **VERIFIED** (48 files).
17. `benchmarks/R/equivalence.R` - **no multinomial scenario exists**
    (`grep -i multinom` returns nothing). So
    `bart2(family = "multinomial")` has no equivalence baseline at all,
    and revision 1's "S4 re-records NONE" was technically true and
    substantively hollow. S0 now adds one.
18. `benchmarks/R/sbc.R` - **no `bart2(` call anywhere**;
    `sbcMakeMultinomial:1360-1362` uses
    `getFromNamespace("bartcoreMultinomialSampler")`. No SBC re-run is
    owed for S1a or S4.
19. `R/bart.R:596-624` `checkFamilyUnsupportedArgs` has **four callers**:
    `:891` multinomial, `:1041` ordinal, `:1072` nbinom, `:1114`
    **hurdle.lognormal**. The `samplerOnly` block is `:604-606`.
    Revision 1's "two lines" was wrong and named three families where
    there are four.
20. `R/A_class.R:548-610` - the `dbartsData` **validity method**, which
    derives `numObservations <- length(object@y)` and carries every slot
    refusal. Absent from revision 1's Fork G scope and price; now in S2.
21. `R/data.R:1259-1270` degeneracy check and `R/spec.R:189-190` sigma
    estimate both fire on a single-trial multinomial response.
    **MEASURED**: the warning text, and a starting sigma at
    `estimateSigmaFromLinearModel`'s host-independent floor (was
    `sigma = 1.52e-16` before the floor landed).
    `dbartsData()` has no `family` formal (`R/data.R:735-746`), so the
    warning gate must key on `@counts`, not on family.
22. `combiner.hpp:728` / `chain.hpp:760` carry the two residue
    strings (retired: both now read "an amplitude coupling", not BCF).

---

## 6. Critique disposition, finding by finding

| # | disposition |
|---|---|
| 1 | ACCEPTED and independently re-verified; **Fork D0 added and recommended**, with my own bitwise probe (train + cutpoints + dispersion, mutability, state round trip) and a source-level soundness trace. Twin deletion demoted to Door 1. |
| 2 | ACCEPTED; the 3ms measurement is now the stated reason NOT to spend a re-record, and the "two engines per fit" framing is qualified as surviving into the RC under D0. |
| 3 | ACCEPTED; draw movement is now stated as conditional on `seed = NA`, with the seeded-identical probe and the `equivalence.R:1396` fixture attribution. |
| 4 | ACCEPTED in full; Fork G re-scoped (validity method into S2), both family gates named with the `dbartsData`-has-no-`family` complication, option (ii) corrected to a hard error. Ledger 20, 21. |
| 5 | ACCEPTED; 4.1 corrected, **Fork G-migration added** with a measured `.hasSlot` mitigation and a NEWS line. |
| 6 | ACCEPTED; four callers named, hurdle called out, per-caller flag recommended, Door 5 opened. Ledger 19. |
| 7 | ACCEPTED; section 1.2 and ledger 3 now say "under an explicit `family = "multinomial"`", with the auto announcement measured. |
| 8 | ACCEPTED; "plausible-looking numeric vector" replaced with the measured constant 1. |
| 9 | ACCEPTED; **a `bart2multinom` scenario is added in S0** so S4's move has a baseline; trio becomes 43/12/11 after S0. Ledger 17. |
| 10 | ACCEPTED; **Fork J** added (below) with a recommended fate and a fallback. |
| 11 | ACCEPTED; ledger 14 RETRACTED and replaced, MAJOR/MINOR re-cited to `dbarts.h:100-104`, C2's cost corrected to two re-bakes, **and Fork C's recommendation changed from C2-now to C3-now**. |
| 12 | ACCEPTED; S2 repriced 420 -> ~700 against its own named comparable, S3 260 -> ~360, test 380 -> ~830 for the merged slice. |
| 13 | ACCEPTED; `bench-sampler` is now marked OWED for Door 3 and not-owed elsewhere. |
| 14 | ACCEPTED and extended; I re-ran it across four cells (keepTrees x n.thin) on both channels. |
| 15 | ACCEPTED; the SBC question is resolved in place, not deferred. Ledger 18. |
| 16 | ACCEPTED and refined upward: **seven** files carry host-shell cells, not five; the ordinal/nbinom blocks pin family-specific text and INVERT under D0; `thresholds.raw`/`dispersion.raw` co-gating recorded in the S4 spec. |
| 17 | ACCEPTED; all anchor corrections applied (plans-vs-design paths, `:334-341`, `:309-322`, `:1757-1759`, `combiner.hpp:730-731`, and the `:42-71` window overstatement). |
| 18 | ACCEPTED; "silent no-op" replaced with "takes effect on an engine nothing reads", with the measured probe. |
| 19 | ACCEPTED; the Rd defect moves from S5 to S0. |
| 20 | ACCEPTED; the no-workaround line is now quoted from my own run, all three families. |
| 21 | ACCEPTED; **restructured to the four-slice floor**, with the deferred residue priced as Doors 1-5. |

### Fork J (NEW, from finding 10). The fate of `bartcore_createMultinomial`

After S4, `bart2` routes through `bartcore_create`'s dispatch arm, but
the 11 multinomial baselines drive
`dbarts:::bartcoreMultinomialSampler`. Leaving both alive means an
internal-only construction path frozen by a fixture - scaffolding in the
reviewed surface under a different name, which is exactly what VD
asked to remove.
- **(J1, RECOMMENDED)** Make `bartcoreMultinomialSampler` /
  `...CountSampler` thin wrappers that build a counts-carrying
  `dbartsData` and call `bartcore_create`, then **verify** the 11 stay
  bitwise (one create either way, same factory, same rng consumption).
  Retire `bartcore_createMultinomial` / `...Counts` in the same slice.
- **(J2)** If the 11 do not stay bitwise, re-record them with a recorded
  reason and retire the entries anyway.
- **(J3)** Keep the C entries as a documented fixture-only path.
  Rejected: it is the scaffolding, renamed.

---

## 7. DISPUTED

**None on substance.** Every finding I probed reproduced, including all
five whose evidence I re-ran independently (1, 3, 7, 8, 14, 18, 20, and
the S4-slot half of 5). I could not refute finding 1's soundness: I
traced the finalizer, the protection slot, the `state` promise's
evaluation environment, and `getPointer()`'s re-creation triple at tip,
and all four support adoption for ordinal and nbinom.

Four **upward refinements**, where re-verification tightened a finding
rather than contradicting it. Recording them so the record is exact:

1. **Finding 11(b)**, "the 'one text-hash cell edit' is likely zero edits
   plus one new `expect_false`." Refined: `test-capi.R:80` pins the
   literal with `expect_identical(hashes$text, "0x6c9776ae1197e8f5")`, so
   C2 edits that line AND adds an `expect_false`, on top of the two
   source re-bakes. Costlier than stated; same direction.
2. **Finding 4**, option (ii) "is a hard bridge error ...
   `R_interface_bartcore.cpp:937-953`." Refined: it dies one layer
   earlier - MEASURED `dbarts(x[0, ], numeric(0))` stops at the R layer
   with `data has zero rows` before the bridge is reached. Both are hard
   errors; same conclusion.
3. **Finding 16**, "F1's test cost is five files, not two." Refined:
   **seven** - the two shell cells plus `test-calibration-midchain.R`,
   `test-forest-weights-r5.R`, `test-ordinal.R`, `test-nbinom.R`,
   `test-multinomial-surface.R`.
4. **Finding 11**, `dbarts_apiSignatureToken` "static_asserted at
   `src/C_interface.cpp:462`." Refined: the `static_assert` statement
   begins at **`:460`**; `:462` is inside its message string.

One **provenance caveat on my own evidence**, stated rather than buried:
the probes ran against a build one day behind tip (see "Verification
method"). Every bitwise probe compares two arms of the same build, so a
draw-moving rng port cannot flip a result; every corresponding source
anchor was read at tip 0f5a2285.

---

## 8. Fork slate for VD

Ordered by consequence. Forks the critique effectively CLOSED are marked
as such, with the closing evidence, rather than re-asked.

**1. Fork D0/D1 - pointer adoption now, twin deletion later. CLOSED by
evidence; confirm only.** Revision 1 said ordinal/nbinom required
deleting a `bartcore_create` and therefore the program's first
equivalence re-record. That is false: assigning the R5's `pointer` field
to the engine that ran delivers the entire win - `$fit` mutable,
`storeState`/`getPointer` working, both `hostFor` sites dead,
`samplerOnly` unblocked - with both creates still happening in order, so
nothing moves. MEASURED bitwise on train, cutpoints and dispersion, with
mutation and state round trips OK. The deleted create was worth 0.003s
against a 0.339s run. ADOPT now (S1a); delete the twin post-RC (Door 1).

**2. Fork C - flat ABI. CLOSED, reversing revision 1.** I recommended
growing the ABI in-arc on the inference that header windows were already
cheap; that inference rested on misreading `0x85bd1ef04beb3848` as a
superseded hash when it is the live signature-half token. Corrected cost
is two literal re-bakes plus test edits, and the named consumer
(stan4bart) is on its own branch with its own release, so nothing is
blocked. RESERVE now (~30 lines in `c-api-growth.md`); build post-RC
(Door 2) when stan4bart asks.

**3. Fork G-migration - old serialized fits. OPEN, cheap, recommend m2.**
Adding `dbartsData` slots breaks every previously serialized fit on first
access to a new slot (MEASURED), and that lands on bartCause, treatSens,
stan4bart and bairrtt users. Options: plain break with a NEWS line (m1),
or route internal reads through `methods::.hasSlot`, which MEASURED
returns FALSE on an old instance without erroring (m2, ~15 lines).
`methods::updateObject` does not exist in base R, so that route is out.
RECOMMEND m2 plus the m1 NEWS line - fifteen lines buys four downstream
packages a clean upgrade.

**4. Fork B1 - `$getLatents` on a multinomial sampler. DECIDED (VD
2026-08-20): DECLINE, written down; see the B1 section for the decided
scope of the rationale.** Original framing kept below for the record.
It returns NULL today because the omegas live in the combiner, not the
response model. r-c-division clause 4 makes read/write symmetry the
default, and the K x n omegas are drawn every sweep, so the status quo is
the inherited accident clause 4 forbids. Build it (~40 engine + ~60 test
plus an ASAN leg), or record a considered decline on the grounds that an
interleaved one-vs-rest augmentation variable refreshed against a moving
margin is not a quantity a host can condition on reproducibly.
RECOMMEND the decline, written down - but this is the one fork where the
principle leans against my recommendation, so it is VD's.

**5. Fork J - the fate of `bartcore_createMultinomial`. OPEN, recommend
J1.** After S4 the two unexported wrappers and their C entries would
survive purely to keep 11 equivalence baselines bitwise: scaffolding
frozen by a fixture. Make the wrappers thin shims over the new
`bartcore_create` dispatch and VERIFY the 11 stay bitwise; retire the
dedicated C entries in the same slice. If they do not stay bitwise,
re-record with a reason and retire anyway. Do not keep a fixture-only
construction path.

**6. Fork I - hurdle `samplerOnly`. OPEN, small, recommend defer.**
`checkFamilyUnsupportedArgs` is shared by four callers, hurdle included,
whose `$fit` is a PAIR of samplers. A per-caller flag unblocks
`samplerOnly` for the three families this arc fixes and leaves hurdle
refused. Whether hurdle should ever return two unrun samplers is a
separate question (Door 5).

**7. Fork A / E / F / H - shape, semantics, guard, wording. CLOSED by
critique endorsement; no objection raised to any.** A4 split-by-family;
E2 mirror the `bart` save/load contract (there is no workaround today -
MEASURED, `$fit$storeState()` does not help for any of the three); F1
delete the `hostFor` mechanism in the last behavioral slice, across seven
test files not two; H3 name the capability rather than the C entry point,
plus a weights arm so `setWeights` stops answering about the response.

**8. Fork B2, G1a - settled by argument, no VD input sought.** No
`$setResponse` sugar for counts (the shapes differ and the guess is the
kind the repo refuses). `data@y` carries the trials vector (option (ii),
`numeric(0)`, is a hard error - MEASURED).

**Shape of the ask:** four slices before the RC - S0 (pins + the shipped
Rd defect + a `bart2multinom` baseline + the `copy` one-liner + the ABI
reserve), S1a (pointer adoption), S2+S3 (the multinomial surface), S4/F1
(direct construction and guard deletion) - at ~1700 non-test / ~1650 test
raw additions, with five priced doors deferred. Only S4 moves a draw, and
it moves exactly the one scenario S0 creates for it to move.

## Landing notes (appended per slice, oldest first)

**S0 LANDED 5e586587** (2026-08-20). All five scope items shipped in one
commit: the pins (inst/tinytest/test-host-shell-pins.R, 307 lines), the
man/bart2.Rd correction (three paragraphs now state $fit rides
keepTrees-or-keepSampler while $bc needs keepTrees), the bart2multinom
scenario appended last in benchmarks/R/equivalence.R with its own literal
seed, the $copy() hostFor transfer, and Fork C3's reserve section in
docs/plans/archive/c-api-growth.md (inst/include/dbarts/dbarts.h untouched).
Price: 5 non-test dense (budget 120) + 389 test-and-harness dense
(stop 390); docs/NEWS/MANIFEST/CI uncounted per convention.

One census correction against section 1.4: predict-replay (63df524e,
landed after this design was written) added `predictForests`, a
fourteenth substantive unguarded method, so the live census is
22/2/14/6 of 44 own methods, not 22/2/13/6 of 43. The pin reads each
method's own body rather than hardcoding the split, which is how the
drift was caught and why the next addition will be caught too. The
other S0-relevant anchors in section 1 were re-verified exact; the
test-capi.R literals in Fork C had drifted by one line and the reserve
doc carries live-verified numbers.

BASELINE RE-RECORD, the program's first: equivalence-00474606.rds, 43
scenarios, named for the commit's own parent (a single combined commit
cannot cite its own hash; S0 adds no src/ file, so the engine binary IS
00474606's). Neutrality partition: all 42 pre-existing scenarios
reproduce equivalence-d15a2bfb.rds bitwise in a non-strict compare
(bart2multinom skipped as not-in-baseline, zero max|z| lines); the new
baseline reproduces itself 43/43 under --strict-coverage;
bcf-equivalence-6e3b9fb8.rds 12/12 and multinomial-equivalence-4d9a3337
.rds 11/11 stay bitwise on every channel. Four-place obligation
discharged in the same commit (equivalence.yaml, MANIFEST, feature-
matrix [f39]; TODO cites no affected baseline name).

Gating: orchestrator diff review (LAND) plus an independent gate
battery: tinytest 6526/0; the full equivalence sequence above re-run;
tests/cpp full + sampler green; --as-cran Status OK; air/lintr/pkgdown/
checkRd/NEWS (292)/doc-freshness clean; the copy-fix pin re-verified to
bite (reverting the one-line fix fails exactly the five copy-laundering
assertions, clean restore passes 32/32). No ASAN owed (no engine code);
no bench.

One CI fire-and-cure: 5e586587's lint leg went red - lint_package()
cannot see tinytest's run-time attachment of expect_* inside a helper
function body, so the two save/reload assertions linted as undefined
globals, while every local pass had loaded tinytest first and was
clean. Cured comment-only at 288ae9e7 with the nolint annotation
test-formula-terms.R already uses; the other five workflows were green
on 5e586587 itself.

Next: S1a (pointer adoption), expected trio 43/43, 12/12, 11/11 bitwise
with NO re-record.

**S1a LANDED b2d1749f** (2026-08-23). Pointer adoption (Fork D0) at the
ordinal and nbinom sites: a documented dbartsSampler$adoptPointer method
(the selfEnv idiom; the caller-owns-soundness contract is in the
docstring) replaces both hostFor assignments, so $fit is the engine that
ran. predict.bartOrdinal/.bartNegbin replay through $fit$getPointer() -
the $bc keepTrees gate stays, so refusal behavior is unchanged while
save/reload-predict now round-trips (verified by inverting the S0 pin
into checkSaveReloadPredicts). samplerOnly unblocks for these two
families via checkFamilyUnsupportedArgs's per-caller allow.samplerOnly
flag; hurdle and multinomial stay refused, both pinned. Method census
is now 45, adoptPointer classed infrastructure in the pins file. The
refusal blocks in test-ordinal.R and test-nbinom.R and the S0 pins
inverted to capability assertions; the other five Fork F files pin
multinomial only (checked, untouched). Multinomial behavior unchanged
end to end.

Price: ~59 non-test dense (budget 190) + ~69 test dense (budget 260);
the adoption mechanism collapsed to nine lines of R/dbarts.R. One
deviation, accepted: predict builds a list(ptr = ...) handle for
bartcorePredict rather than delegating to $fit$predict, keeping
validation byte-identical (bartcorePredict reads only $ptr, and $K,
which is NULL either way for these families); full delegation is S4's
territory. One orchestrator amend pre-landing: the NEWS phrase "for
every family" tightened to "all three host-shell families" (gaussian
and binary fits always had the documented storeState workaround).

Gating: orchestrator diff review (LAND) plus an independent gate
battery: tinytest 6538/0; trio 43/43 (--strict-coverage, 43 compared /
0 skipped), 12/12, 11/11 identical-draws lines, zero max|z| lines -
NEUTRAL exactly as designed, NO re-record; tests/cpp full + sampler
green; TWO independent mutation probes (implementer's and the
gate-runner's own hostFor revert) fail the inverted assertions and
restore green (under mutation: test-ordinal.R 7 failures,
test-nbinom.R hard-errors at refuseHostRead, test-host-shell-pins.R
reproduces the NULL-external-pointer defect); --as-cran Status 1 NOTE
(days-since-update only); air/lintr (8 files)/pkgdown/NEWS (293
entries)/doc-freshness all clean. No ASAN owed (zero src/ files
touched); no bench.

Next: S2+S3 (multinomial surface); ASAN owed for the bridge dispatch
arm.

**S2+S3 LANDED 5a3bc276** (2026-08-24, one merged commit). The public
multinomial sampler surface: dbartsData gains counts/offset.category/
offset.category.test slots constrained in the validity method, data@y
derived as the trials vector so the two cannot disagree (Fork G1/G1a);
dbarts(x, y, family = "multinomial") accepts a counts matrix or a
factor/character/integer-code vector (one-hot expanded), matrix
interface only, and counts-carrying data resolves from family =
"auto"; the bridge creation arm builds from the spec triple so
re-creation mirroring is free - the slots ARE what creation reads, no
reapply path (the design's own G1 rationale, delivered); the three R5
methods setCounts/setCategoryOffset/setCategoryTestOffset ride the
existing internal channels; the K-matrix predict arm; the Fork B
refusal matrix with capability-named messages (H3 repoint executed;
the weights conduit has its own arm); B1's getLatents decline written
into dbartsSampler-class.Rd and multinomial.md, scoped to the
augmentation families with AFT's imputed log-time read preserved.
Host-shell census now 48 own methods / 25 host-mutation-guarded.
Beyond-spec, kept: nine R-side refusals for compositions the
multinomial factory silently drops (monotone, DART, split.probs,
linear/GP leaves, k hyperprior, named prior.scale, groups,
storage = "single", variance forest) and a bridge nullptr backstop
for the variance-forest case, which previously crashed latent.

Price: 576 non-test dense (budget 1060, stop 1590) + 428 test dense
(budget 830, stop 1245). Under budget because G1's mirroring-for-free
held exactly as designed.

Review cycle, three rounds: Opus implementer -> Opus adversarial
review (BLOCK: bart2()/xbart() accepted a counts-carrying dbartsData
through family = "auto" resolution and silently returned a
mis-reshaped single-forest fit - zero warnings; plus formula-LHS
silently discarded at the data layer, a dead .hasSlot guard inside
the validity method, creation-side bridge checks weaker than
mutation's, an NA-sentinel predict edge regressed, two codename
comments) -> fixes -> reviewer re-measurement (residual: bart() still
accepted the same object through the passthrough branch; bridge
refusal misreported K on the public route) -> fixes -> reviewer LAND
(all cures re-measured; empirical entry sweep: every exported
single-location packager refuses a counts-carrying dbartsData -
bart/bart2/xbart by the new gate, its mutation proof failing exactly
the 5 pins, rbart_vi by the groups gate now pinned, pdbart by its
type check). The lesson the pin carries: family = "auto" resolution
makes a new family reachable from EVERY dbartsData-accepting entry,
not just the ones the slice built.

Gating on the final tree (landed as 5a3bc276 after a
patch-id-identical rebase onto the docs-only records tip): two
independent legs. Gate battery PASS: tinytest 6686/0; trio 43/43
strict (43 compared / 0 skipped), 12/12, 11/11 identical-draws lines,
zero max|z| - the dispatch arm is unreachable from shipped paths
until S4, verified not asserted; tests/cpp full + sampler; ASAN/UBSAN
39 files zero reports (macOS mechanics: sanitizer flags via
R_MAKEVARS_USER, run through bin/exec/R with DYLD_INSERT_LIBRARIES -
the bin/R wrapper strips the insert under SIP); --as-cran 1 NOTE
(days-since-update); air/lint_package/pkgdown/NEWS (296)/
doc-freshness clean. Reviewer leg re-ran both equivalence harnesses
and eight multinomial-adjacent test files independently, same counts.
clang-tidy warnings in the bridge all blame to pre-slice commits.

Next: S4/F1 (bart2 multinomial constructs direct; $bc deleted;
hostFor mechanism deleted; reproduction gate re-pointed; Fork J
discharged; re-records exactly the S0 bart2multinom scenario).

**S4/F1 LANDED 2619ac9e** (2026-08-24) - the arc's final behavioral
slice. bart2(family = "multinomial") constructs its sampler directly
through the public dispatch arm: one bartcore_create per fit, $fit IS
the sampler that ran, save/reload-predict round-trips for all three
families now. $bc deleted everywhere (thresholds.raw/dispersion.raw
keep their keepTrees gate and readers); the hostFor mechanism -
field, refuseHostMutation, refuseHostRead, every call site, the Rd
sentences, the S0 copy() transfer - deleted with zero survivors
(host-shell-read-guards closed as OBVIATED in the TODO). Fork J1:
bartcoreMultinomialSampler/...CountSampler are thin wrappers over
the public path, the bartcore_createMultinomial/...Counts C entries
retired from bridge, header, and registration; the 11 multinomial
baselines verified bitwise through the wrappers. samplerOnly
unblocked for multinomial (hurdle stays refused); the
refuseCountsCarryingData gates stay (bart/bart2/xbart do not serve a
counts-carrying dbartsData; dbarts/dbartsSpec do). test-host-shell-
pins.R rewritten to capability assertions (census now 46 own methods
/ 41 substantive, name-read since the guards are gone); the
reproduction gate re-pointed with independent teeth: bart2's direct
construction vs a hand-built dbartsSpec + wrapper construction,
different one-hot code and different triple resolution converging on
the shared factory, with a verified-discriminating one-hot mutation.

THE RE-RECORD, the program's second: equivalence-5a3bc276.rds (43
scenarios, named for the commit's parent), moving EXACTLY the S0
bart2multinom scenario - partition vs the old baseline verbatim
"bart2multinom 87 summaries, max |z| = 2.06", 42/43 bitwise;
recorded reason: one create instead of two. Four-place obligation
in-commit: equivalence.yaml, MANIFEST (old baseline demoted to
historical, file kept per S0 precedent), feature-matrix [f39]; TODO
cites no affected name. bcf 12/12 and multinomial 11/11 untouched.

Price: 169 non-test dense (budget 330) + 142 test dense (budget
300). Review cycle: Sonnet implementer -> Opus adversarial review
(LAND-WITH-FIXES: dead categoriesAreDeclared branch + false comment
in the bridge; three stale 1.0-0 NEWS items describing $bc/hostFor
as live, added by S0/S1a after the design's 4.8 census; codename
comments; stale multinomial.md creation path; plus small message/
pin/alignment items - all fixed) -> re-measurement (BLOCK on one
content-free defect: the amend re-stamped the baseline's meta$rev to
the dangling pre-amend hash; repaired by restoring the byte-identical
correctly-stamped recording, orchestrator-applied) -> final LAND.
Gate battery PASS twice (once pre-landing; the final tree differs only
by that .Rbuildignore'd baseline metadata revert, package content
identical): tinytest 6672/0; the re-record sequence as above;
tests/cpp; mutation probe (cyclic one-hot shift caught by 6/126 +
9/11, restore-green confirmed); --as-cran 1 NOTE; air/lint_package/
NEWS (293)/doc-freshness/pkgdown clean. No ASAN owed (deletions
plus one R-level check). Landed as 2619ac9e after a
patch-id-identical rebase onto the records tip.

Process notes for the record: the implementer's own subagents
performed unauthorized git operations mid-slice (a commit and an
amend); content was independently re-verified end-to-end by both
verification legs before landing, and later rounds kept git in the
implementer's own hands. A gate-runner self-caught pointing R CMD
check --library at its live shared lib (mid-run clobber, false
failures) - isolate check libraries per run.

ARC BEHAVIORAL WORK COMPLETE. Remaining arc close-out: the
serialized stale-cite anchor sweep and residual codename/scaffolding
sweep across docs (the Queue's recorded list), the feature-matrix/
r-c-division/TODO closures already partially discharged in-slice,
then the exhaustive human review VD sequenced after the arc.
