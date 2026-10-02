# nbinom-dispersion-name: the negative-binomial r is named shape

Status: LANDED 2026-10-02 (32ea44b2, 0a4680c6) under dec-B189, dec-B190 and dec-A145 in [decisions.md](../decisions.md).

agent: sonnet implementer, one, serialized; opus reviewer.
rng: NEUTRAL. A rename: no draw, no RNG call order and no default moves. Every equivalence scenario stays
identical and every exact gate passes unchanged.
window: pre-release, before the 1.0-0 merge. The old name never reached a release, so nothing is tombstoned.
budget: ~700 changed lines of code, tests and manual, ~150 of docs. Plans have run 1.5-2x low.

## Goal

The negative-binomial parameter r, the shape of the gamma distribution the Poisson rates are mixed over, is called
shape everywhere the package names it: the family argument, the sampler's reader, extract's type, the fit's and
run's components, the saved state, the shipped C header and the engine's own vocabulary. The word dispersion no
longer names r anywhere in the tree outside frozen records.

## Context

The name sits today on `nbinom(dispersion = )`, `getDispersion`, `extract(type = "dispersion")`, the fit's
`dispersion` and `dispersion.raw`, `run()$dispersion`, the saved state's `dispersion` block, two struct fields of
`inst/include/dbarts/dbarts.h`, and engine names such as `NBDispersionPrior`. The rulings and what each covers are
in the three ledger entries the Status line names. The model itself is described in
[negative-binomial.md](../design/negative-binomial.md).

## The mapping

| Today | After |
|---|---|
| `nbinom(dispersion = NULL)`, its settings entry and its messages | `nbinom(shape = NULL)` |
| `$getDispersion()` and its C entry point | `$getShape()` |
| `extract(type = "dispersion")`, summary's `vars = "dispersion"` | `"shape"` |
| the fit's `dispersion`, `dispersion.raw` | `shape`, `shape.raw` |
| `run()$dispersion` | `run()$shape` |
| the saved state's `dispersion` block | `shape` |
| the header's two `dispersion` fields | `shape` |
| engine and bridge identifiers (`NBDispersionPrior`, `setDispersion`, the state slot enum, locals) | the same with Shape or shape |
| `inst/tinytest/test-dispersion-channel.R` | `test-shape-channel.R`, by `git mv` |

`getShape` is an overloaded reader (dec-B190): it returns the shape parameter of whatever family the sampler runs,
one value per chain, and NULL on a family that has none. Only nbinom carries one today. Its documentation says so
in those terms rather than as a negative-binomial reader.

## Constraints

- The CONCEPT keeps its word. "overdispersion", "overdispersed counts" and a comparison to glm's dispersion stay
  as they are; only text that names the parameter r moves. Where a sentence said "the dispersion r", it says
  "the shape r".
- The manual states the direction once at the family entry and once at `getShape`'s value: the variance is
  mu + mu^2 / shape, a larger shape is closer to Poisson, and the same quantity is `size` in `rnbinom` and `theta`
  in MASS and mgcv.
- No tombstone and no alias for the old spellings. A top-level `dispersion =` was already an unused-argument
  error; the tests that pin a family-only argument being refused at top level move to `shape`.
- The saved-state block is renamed with no change to either state-format version constant: no serialized format
  has shipped, as the registry comment at `stateFormatVersion` says.
- A renamed header field moves the header's API hash. Re-bake `DBARTS_C_API_HASH` as the header's own comment
  describes; no version constant moves.
- No baseline is re-recorded. The equivalence harness keeps whatever key names the recorded baselines carry for
  this channel, reads the fit's new component into them, and says so in a comment.
- Windows `*.win` variants move with their autoconf counterparts.
- Records are not reworded. Everything under `docs/plans/` other than this file, `docs/decisions.md` and the
  review directories describe what was true when written: there, change only a cite the freshness check rejects,
  marking it retired and naming the new spelling. Documents that describe the present - the design docs' current
  sections, `docs/architecture.md`, the feature matrix, `benchmarks/README.md`, `TODO` item text - move to the new
  name. TODO item ids do not change.
- `inst/NEWS.Rd` names only the final spelling.
- Out of scope: a dispersion spelling for 1/r, a real-valued r, any change to the grid or its prior, and the
  consumer packages, none of which reads this channel.

## Steps

1. Engine, bridge and header: identifiers, the state block name, the two header fields and the hash, with
   `tests/cpp` following.
2. R: the family constructor and its validation, the reader, extract, summary, plot, diagnostics, the fit's
   components and the augmentation helpers.
3. Manual pages, `inst/NEWS.Rd`, tinytest, `benchmarks/R`, the exact-gates workflow comment.
4. Present-facing docs and the cites the freshness check rejects.
5. One commit for steps 1-3 and one for step 4.

## Verification

Against a private library, installed with `--preclean`:

- `cd tests/cpp && make && ./test_bartcore` passes.
- `tinytest::test_package("dbarts")` passes with no new warning.
- Every gate script `.github/workflows/exact-gates.yaml` lists passes in `quick` mode.
- The equivalence compare on a reference build reports identical draws for every scenario against the current
  baseline in `benchmarks/baselines/MANIFEST`, and the four seeded-drift snapshot files pass on that build.
- `lintr::lint_package()`, `air format --check .`, `tools/check-rc-codoc.R`, `tools/check-win-drift.R` and
  `tools/check-doc-freshness.R` pass, each on its own exit status.
- `git grep -i dispersion -- R src inst man tests benchmarks/R .github` lists only the concept, glm comparisons and
  the baseline key names; the report gives that list.

## Landing

LANDED 2026-10-02: the rename (32ea44b2) and the present-facing docs (0a4680c6).

- One departure from the mapping: the facade virtual is `shapeParameter(chainNum)`, `shape()` there already
  returning the sampler's `SamplerShape`. `DBARTS_C_API_HASH` is 0xa33182bf349aa60d.
- The review, an independent reader with its own builds, found the rename had also taken the concept's word in
  four places - the AFT variance forest's "dispersion of log survival time" in two manual pages and two test
  comments - and those were put back before landing. Nothing else in the diff was other than a substitution.
- A state carrying the old block name is refused by an nbinom sampler as inconsistent, not misread.
- Verification, by the reviewer against its own libraries: tests/cpp; tinytest 11621 tests, 0 failures, the
  warning table identical to the base commit's; the lint chain; `R CMD check --as-cran` with the Date NOTE only;
  the two negative-binomial exact gates in quick mode; a mutation of the recorded channel turning
  test-shape-channel.R red; stan4bart's suite built against the new header; on a reference build the equivalence
  compare identical on 55 of 55 scenarios, multinomial 11 of 11, BCF 15 of 15, and the four snapshot files.
- Left as they were: R locals named `disp`, and the baseline key names in the equivalence harness.
- Found on the way, not from this slice: one tests/cpp run in sixteen failed the monotone slow-count tally check
  (TODO, monotone-slow-count-test-intermittent).
