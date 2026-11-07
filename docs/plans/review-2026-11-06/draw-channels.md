# Review 4, slice S-B: draw channels, the flat run and the cut-grid floor

Status: built, in review (branch review4-draw-channels).

rng: neutral. Every run draws bit for bit as before, with or without a draw callback and whatever the caller
stores; the equivalence baselines and the four snapshot files are unchanged ([Gates](#gates)). What changes is
what a draw callback sees, the scalar sigma a sampler with a variance forest reports through the C API, which
runs write `results.k`, the uniform grid over a span wider than `DBL_MAX`, and which cut grids are refused.

Findings fixed: bridge-02, bridge-03, bridge-04, bridge-V-02, ingest-04, ingest-02 and docs-10's share in the
shipped header, of the fourth whole-repo review, with the mutation-04 test addition; as amended by the blind
critique of the slice (amendments 1 to 7 and 9 to 13 binding, amendment 1's accessor and amendment 8 superseded
by dec-B615). Rulings this rests on: dec-B600 (a registered callback receives every channel the fit carries,
the log-likelihood included, paid only by a run with a callback), dec-B601 (sigma is NaN in every draw of a
sampler with a variance forest), dec-B604 (setCutPoints refuses a cut that is not finite; a saved state's grid
takes the same floor) and dec-B615 (no accessor here; the sampler facts move into one struct later, TODO
capi-sampler-info).

ABI: not an ABI event. `inst/include/dbarts/dbarts.h` changes in comments only; no declaration, enumerator or
layout moves, so `DBARTS_C_API_HASH` and `dbarts_apiSignatureToken` are unchanged and the header's
`static_assert` on the token still holds.

## Fixes

### bridge-02: a draw callback receives every channel the fit carries

Before: a callback saw a channel only where the caller stored it. The flat run with `results` NULL handed
null train, test and variance pointers, and the R run's draw log-likelihood was always null, though the header
promised every channel.

Now: `bartcore::Results` carries a per-draw stride for every draw-indexed array channel (`glueStride`,
`variableCountsStride`, `splitProbabilitiesStride`, `logLikelihoodStride` and `ordinalThresholdsStride` join
the five per-observation ones), read by `Sampler::run`'s per-chain slab base and copy, `storeSample`'s write
and `fillDraw`'s read. `Chain::carriedChannels` states once which channels the fit carries and its coupling
defines; `fillDraw` hands null by it, and `Sampler::provideObserverScratch`, called by `Sampler::run` only when
a draw callback is set, points each carried channel the caller left null at one draw of scratch per chain with
stride 0. The scratch is owned by the run's frame, so every exit frees it; no callback, nothing allocated. The
R run's `keepFits = FALSE` scratch for the per-observation channels stays the bridge's; neither bridge
allocates anything new (amendment 3, [Calls made](#calls-made) 1).

Tests: `testDrawCallbackCarriesEveryChannel` (test_sampler.cpp) - gaussian under DART, ordinal, a variance
forest and BCF, each on two chains and two threads over an empty `Results`: each carried channel non-null on
every call, every other null, each equal to a twin storing every channel at the same draw. test-capi.R - the
consumer's recording callback (capi/consumer.c) over five arms (DART, ordinal, variance, aft with a variance
forest, BCF) on four routes (flat with `results` NULL, flat with every field, R at `keepFits` TRUE and
FALSE): the channels seen match the fit, and the copy each takes at chain 2, draw 2 equals the stored array of
the route that stores it, the log-likelihood among them.

### bridge-03: `kIsSampled` follows the model as it now stands

Before: `Sampler::kIsSampled` returned the creation option, so after a `setLeafPrior` or `setModel` that gave
or took away a k hyperprior the flat run wrote `results.k` where k was fixed or left it where k was drawn,
while the R run computed its own answer.

Now: `Sampler::kIsSampled` computes it from chain 0's forest calibrations (forest 0 draws k, or a forest
sharing a finite prior total does), the loop `bartcore_run` held; `bartcore_run` reads `shape.kIsSampled`, so
both routes have one source.

Tests: test-capi.R - a fixed-k sampler given `normal(k = chi(1.5, 2))` through `setLeafPrior` and through
`setModel` reports 1, writes `results.k` with varying draws, and the R run returns k; given a fixed k again it
reports 0, leaves `results.k` at the caller's sentinel, and the R run returns no k (amendment 9).

### bridge-04 and bridge-V-02: sigma is NaN under a variance forest

Before: a sampler with a variance forest reported its pinned global sigma (1 times the response scale) in the
draw and in `results.sigma`, which is no residual scale of the model.

Now: `fillDraw` and `storeSample` write NaN wherever the chain has a variance forest. The R run already
returns no sigma there.

Tests: `testDrawCallbackCarriesEveryChannel`'s variance arm, NaN in the draw and the stored slab; test_fuzz's
liveness invariant takes NaN on a heteroscedastic sampler; test-capi.R - NaN in both on gaussian and aft with a
variance forest, and the pinned 1 on probit and on gaussian at `sigma = fixed(1)` (amendment 7's control).

### ingest-04: the uniform grid over a span wider than `DBL_MAX`

Before: a column with ends of opposite signs near `DBL_MAX` gave an infinite increment and infinite cut
points, so the column could not be split and its record could not be reloaded.

Now: `ColumnStore::fillCutsOverRange` keeps the increment loop wherever the span is finite, bit for bit, and
otherwise places each cut at `xMin (1 - t) + xMax t`, every term bounded by an end.

Tests: `testUniformGridOverflow` (test_data.cpp) - (-1.7e308, 1.7e308) and (-DBL_MAX, DBL_MAX) at 10 and 100
cuts finite, strictly increasing and inside the range; (0, 1) and (-3, 7) bitwise the old loop (amendment 11);
a +-1e308 dense and sparse column with zeros. test-cut-grid-distinct.R - 40 values from -1e308 to 1e308 at
n.cuts 5 on a dense column, a sparse one holding zeros and a setPredictor refresh: 5 finite cuts, the column
split, and the copy runs.

### ingest-02: no grid holds an infinite point

Before: `setCutPoints` took an infinite cut, and a state or warm-start donor whose grid ended at an infinity
was installed.

Now: `cutGridIsValid` requires finite ends, so with a strict increase every point is finite; the state's grid
check and the donor's take it, and `readCutPointsRecord`'s own finite check goes as implied.
`bartcore_setCutPoints` refuses a non-finite cut after the NA check: "$setCutPoints: 'cuts' must be finite; a
cut at -.Machine$double.xmax separates values at -Inf". The help line and a NEWS line (0.9-34 took one) are
added.

Tests: `testDegenerateGridRestores` (test_state.cpp) - {0.5, Inf}, {-Inf, 0.5}, {NaN} and {Inf} refused, a
state ending at Inf or starting at -Inf refused with the sampler unchanged. test-cut-grid-distinct.R - an
infinite end refused through `setState` and through a warm-start donor by the existing messages.
test-cut-points-undo.R - Inf and -Inf refused by the new message, and a cut at `-.Machine$double.xmax` taken,
the copy running.

### mutation-04: the repeated-cut refusal names every pair

The refusal already scans every adjacent pair; test-cut-grid-distinct.R gains a repeat at the first pair and
at the last, through a state and through a donor, so a scan that starts late or stops early is caught.

### docs-10: the header says caller

Each `host` in dbarts.h reads "caller".

## Gates

Run on arm64 macOS against the branch's own libraries (shipped and reference builds, each installed with
`--preclean`).

- tests/cpp: every suite passes, plain and under `-fsanitize=address,undefined`, 561 checks reported ok in
  each, no sanitizer report.
- Full tinytest, one process, shipped build: 42283 results in 303 files, 0 failed.
- Equivalence trio on the reference build, `--bitwise` (and `--strict-coverage` on equivalence.R): 55, 19
  and 11 scenarios compared, 0 skipped, every one identical draws.
- The four snapshot files on the reference build: 14, 3, 7 and 3 assertions, 0 failed (unchanged).
- R-loaded ASan and UBSan (the bridge and the flat C entry file compiled with
  `-fsanitize=address,undefined`, run through `bin/exec/R` with the runtime inserted): the four touched files
  (test-capi.R, whose compiled consumer runs every channel route, test-callback.R, test-cut-grid-distinct.R,
  test-cut-points-undo.R) and six neighbors (heteroscedastic, its channels, aft heteroscedastic, BCF
  log-likelihood, ordinal, degenerate cuts), 1617 results, 0 failed, no report.
- Consumer sweep, each consumer built from source against the branch's library and against main's
  (980bc50f9) into its own library: stan4bart 1860 results, 0 failed, on both; treatSens 491 expectations in
  103 tests, 1 failed, 10 skipped, on both - the same failure, in treatSens' treatment-model forwarding
  test, which expects `bart` to be unresolved where the sweep's driver had attached dbarts.
- R CMD check --as-cran on the built tarball: Status 1 NOTE, CRAN incoming feasibility only (an external
  site's expired certificate, cited from two help pages this slice does not touch).
- `lintr::lint_package()` 0 lints; `air format --check .`, check-rc-codoc.R (53 methods), check-win-drift.R,
  check-doc-freshness.R and `mutation-battery.R verify-anchors` all OK.
- Mutants, each built into a copy of the tree, the engine ones against tests/cpp and the rest installed and
  run against their test file:

  | mutant | killed by |
  |---|---|
  | the observer's scratch left at the channel's natural stride | a trap in `testDrawCallbackNeutrality` (the chains' one-draw blocks overrun) |
  | the observer's scratch skipped | `testDrawCallbackCarriesEveryChannel`'s presence check, 8 checks |
  | the log-likelihood left out of the scratch | the same presence check, 6 checks |
  | sigma NaN in the draw only | the variance arm's stored-twin check and test_fuzz's heteroscedastic liveness, 3 seeds |
  | sigma NaN in the stored slab only | the variance arm's stored-twin check |
  | sigma NaN keyed on a fixed sigma in place of the variance forest | 46 checks, the ordinal arm's first; test-capi.R's probit and `sigma = fixed(1)` arms are written for it |
  | `kIsSampled` returning the creation option | test-capi.R, 8 results, the setLeafPrior arm's first (tests/cpp does not reach it) |
  | the floor on the first cut only; on the last only | `testDegenerateGridRestores`' finite floor and its Inf-end (resp. -Inf-start) setState refusal |
  | setCutPoints' finite refusal removed | test-cut-points-undo.R, 8 results |
  | the old increment in the overflow branch | `testUniformGridOverflow`, 3 checks |
  | the repeat scan starting at the second pair; stopping before the last | test-cut-grid-distinct.R's first- and last-pair arms, 4 results each |

## Calls made

1. Amendment 3 is taken: the scratch is allocated once, in `Sampler::run`, from the predicates `fillDraw` reads
   (`Chain::carriedChannels`), not in the two bridges. Alternative: the plan's bridge-side scratch per
   amendments 2 and 4, two allocations to keep in step with the engine's predicates.
2. The ordinal thresholds' stride for a caller's buffer keeps `results.numOrdinalThresholds`; the scratch takes
   the response's count. Alternative: one count for both, which would move a caller's layout.
3. A state or donor whose grid has an infinite end is refused by the existing messages ("state is not
   consistent with this sampler", the donor's grid refusal), as the plan states. Alternative: a message naming
   the infinity, a new string at two doors for a hand-edited object.
4. The setCutPoints message names `-.Machine$double.xmax` as the cut that separates values at `-Inf`, per
   dec-B604. Alternative: refuse without the hint.
5. mutation-04 is tests only: the refusal already names every pair, so there is no code to fix.
6. Amendment 7's three bridge-02 mutants are equivalent under call 1: with no copy of the stride into a chain's
   `Results` on the scratch path (`provideObserverScratch` sets stride 0 directly), dropping a stride from the
   copy changes nothing there, and the R route has no log-likelihood claim inside `if (!keepFits)` to remove.
   Named in their place: the scratch left at the natural stride, and the scratch skipped
   ([Gates](#gates), mutants).
7. Amendment 9's setModel arms are included beside setLeafPrior's.
8. An aft arm with a variance forest is added to the channel and sigma tests, the second family that can carry
   one.
9. The NEWS line for the infinite cut is in UPGRADING, since it can stop 0.9-34 code. Alternative: BUG FIXES
   beside the NaN refusal.
10. docs-10 changes each `host` in dbarts.h to "caller"; every hit is about the C caller of the API.
    Alternative: the long form where the sentence is about running the sampler inside a larger one.

## Design-doc register

- inst/include/dbarts/dbarts.h (comments only): the PRESENCE AND ABSENCE paragraph says every carried channel
  reaches a callback, the log-likelihood included, and sigma is NaN under a variance forest; the results
  paragraph says k follows the model and the log-likelihood is computed for a callback; `dbarts_draw::sigma`
  and `dbarts_sampler_kIsSampled` gain their comments; "caller" for "host".
- docs/architecture.md, "Callbacks and errors": the per-draw callback sees every carried channel through
  `Sampler::run`'s scratch, the log-likelihood computed for it, sigma NaN under a variance forest.
- docs/design/per-draw-callbacks.md: the channel paragraph states presence and the NaN sigma; the stride
  paragraph says every draw-indexed channel has a stride and `Sampler::run` supplies the observer's scratch
  (`Chain::carriedChannels`).
- man/dbartsSampler-class.Rd: the run method's callback item (amendment 12); setCutPoints' cuts are finite.
- inst/NEWS.Rd, UPGRADING: setCutPoints refuses an infinite cut point, which 0.9-34 took.
