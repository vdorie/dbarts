# leaf-vocabulary: leaf where a name means a leaf

agent: sonnet (R surface, bridge strings, manual; no engine change)
rng: neutral (names only; the tombstone path fits bitwise what the new spelling fits)
budget: dbarts ~700 changed lines over ~110 files, nearly all one-token; stan4bart ~20 lines / 6 files; bartCause 1 line

## Goal

Every user-facing name that belongs only to terminal nodes says leaf:
- `leaf.prior` on every entry point;
- `$sampleLeafParametersFromPrior`;
- `$getLeafPrior` and `$setLeafPrior`;
- the reader's `leaf.scale.factor` and `leaf.scale.divisor` columns.

Two old names that shipped in 0.9-x stay until 1.1-0 as once-per-session warning tombstones: `node.prior` on
`dbarts()`, and `$sampleNodeParametersFromPrior`. No other rename gets a tombstone. stan4bart and
bartCause move in lockstep.

## Context

- Rulings in docs/decisions.md:
  - dec-B129: the rename.
  - dec-A02: a name 0.9-x shipped warns once, names its successor, uses the value, and expires at 1.1-0.
  - dec-B128: a name that existed only on the development branch gets no tombstone. It gets an unused-argument
    error.
- What 0.9-x shipped: `node.prior` on `dbarts()` only, and `$sampleNodeParametersFromPrior`. 0.9-x's `bart`,
  `bart2` and `xbart` took `k`. `node.prior` on `bart`/`xbart`, `getCalibration`/`setCalibration` and the
  `node.scale.*` columns are branch-only, so NEWS names them only in their new spelling.
- The `dbartsSpec` rule: a tombstone lives only on an entry point where 0.9-34 had the argument. `dbartsSpec` is
  new, so it takes none: `sigma`, `node.prior` and `resid.prior` are refused by name, each message naming its
  successor (`sigest`, `leaf.prior`, `family = gaussian(sigma = )`). The registry holds no entry with owner
  "dbartsSpec", and `dbartsSpec` carries `...` only so that
  [`refuseForeignFrontDoorArgs`](../../R/tombstones.R) can say so.
- Tombstone machinery:
  - Registry: [`dbartsTombstones`](../../R/tombstones.R).
  - Once-per-session flag: [`warnOnce`](../../R/utility.R) over [`onceWarnState`](../../R/utility.R).
  - Renamed-formal template, `resolveRenamedSigma`: the old formal stays before `...`. Supplying both spellings
    is an error ("'sigma' and 'sigest' name the same estimate on 'dbarts'; supply one"). Otherwise it warns once
    and the value is used. The door then rewrites `matchedCall` to the new name, so the stored call carries one
    spelling.
  - Method template: [`noOpThreadMethod`](../../R/tombstones.R).
  - inst/tinytest/test-tombstones.R checks every registry entry against NAMESPACE, the owner's formals, the
    generator's methods and the NEWS tombstone-list item (marker "Every name and argument kept reachable").
- Foreign-argument refusal: [`refuseForeignFrontDoorArgs`](../../R/tombstones.R) refuses any name in a door's
  `...` that is neither a tombstone nor accepted, with "unused argument 'x' passed to 'bart'". On `bart` it runs
  after `forwardToLegacyDoor` and `refuseLegacyPositionalCall`.
- Prior plumbing:
  - [`resolveSamplerSpec`](../../R/spec.R) turns the door's `matchedCall` into a call to
    [`parsePriors`](../../R/model.R), then fills defaults with [`setDefaultsFromFormals`](../../R/utility.R).
    `redirectCall` keeps only names that are formals of `parsePriors`, which has no `...`.
  - `bart` builds its prior in [`buildSamplerPriors`](../../R/bart.R) and forwards it through
    [`buildHostSamplerCall`](../../R/bart.R). [`bartBT`](../../R/bart.R) forwards through its `args` list.
    [`xbart`](../../R/xbart.R) reads the prior argument itself.
- Consumers of internals:
  - treatSens (dbarts-1.0) calls `parsePriors(..., node.prior = )` by name and reads `priors$node.prior`.
  - CRAN stan4bart 0.0-13 calls `parsePriors` positionally, reads `bart_priors$node.prior`, and builds
    `new("dbartsModel", ..., node.scale = )`.

  Keeping `parsePriors`' formal, its list names and the model slots keeps both working.
- Engine, state and header:
  - The saved-state blocks already say `leaf.scale` (see [`stateFormatVersion`](../../src/R_interface_bartcore.cpp)).
  - No engine identifier changes.
  - inst/include/dbarts/dbarts.h has no node-prior name. Its `dbarts_leaf_model` enum already says leaf, and one
    comment reads "an R-side calibration read". DBARTS_C_API_HASH does not fold comments.

## Constraints

- RNG neutral; bitwise equivalence against the MANIFEST baselines is the proof.
- dbarts.h stays untouched, so there is no ABI event. Its comment waits for the host-neutral header change
  (dec-B85).
- These stay as they are (see Calls made):
  - `parsePriors`' `node.prior` formal and its returned list names;
  - the `dbartsModel` slots `node.prior`, `node.hyperprior` and `node.scale`;
  - the classes `dbartsNodePrior` and `dbartsNodeHyperprior`;
  - tests that read `sampler$model@node.scale` and `@node.prior`;
  - engine identifiers: [`Chain::sampleNodeParametersFromPrior`](../../src/bartcore/chain.hpp),
    [`ForestCalibration`](../../src/bartcore/chain.hpp), the facade virtual and `nodeScaleFactors_`;
  - test file names, since test-calibration-*.R are cited by history links.
- These keep node, per the ruling:
  - `$getTrees`' rows;
  - `plotTree`'s `nodeWidth`, `nodeHeight` and `nodeGap`;
  - the tree prior's per-node split probability, and interactions.Rd's "per-node admissibility rule";
  - `gp()`'s `max.leaf.size`;
  - the reader's `leaf.model` attribute;
  - docs/architecture.md's "node values ride a RAWSXP", which is about saved tree nodes of either kind.
- Frozen records stay as written:
  - docs/decisions.md, docs/plans/archive/ and docs/plans/review-2026-08-24/;
  - landing notes, and designs and plans marked LANDED or COMPLETE;
  - NEWS sections before 1.0-0.

  Only a cite the rename breaks gets repaired.
- A concurrent predict/extract slice edits R/generics.R. Rebase onto it before gating.
- benchmarks/R/classic-compare.R keeps `node.prior` and `sigma`: it runs the same script under 0.9-34 and 1.0-0.

## Site inventory (dbarts at origin/bartcore)

R/:
- `dbarts` and [`dbartsSpec`](../../R/spec.R): `node.prior = normal` becomes `leaf.prior = normal` in place. `dbarts`
  gets a tombstone formal `node.prior = NULL` beside `sigma`; `dbartsSpec` gets none and refuses the name.
- `bart`: formal `node.prior = NULL` becomes `leaf.prior = NULL` in place, with no tombstone formal. `bart2` copies
  bart's formals. In `buildSamplerPriors`:
  - `matchedCall[["node.prior"]]` becomes `leaf.prior`, and so does the returned element;
  - the object name given to [`refuseColliding`](../../R/model.R) becomes "leaf.prior", so the k collision message
    names `leaf.prior`. A `node.prior` plus `k` call on `bart` never gets there, because
    `refuseForeignFrontDoorArgs` refuses it first.

  Also `buildHostSamplerCall`'s `samplerCall$node.prior`, `bartBT`'s `priors$node.prior` and `args$node.prior`, and
  the comments in `bart2Hurdle` and `detectAutoOrdinal`.
- `xbart`: formal `node.prior = NULL` becomes `leaf.prior = NULL`, with no tombstone. Also:
  - `matchedCall[["node.prior"]]` becomes `leaf.prior`;
  - `resolvedAs("node.prior", ..., "node prior specification")` gets leaf in both labels;
  - the locals `node.prior` and `node.spec` become `leafPrior` and `leafSpec`, not `leaf.prior`, which would shadow
    the formal;
  - the comments in [`xbartRunUnits`](../../R/xbart.R).
- [`refuseForeignFrontDoorArgs`](../../R/tombstones.R) gets a hint map, `list(node.prior = "leaf.prior")`. For a
  refused name listed there, the message adds "; the leaf prior is 'leaf.prior'". This is a refusal, not a
  tombstone: no registry entry, no warning and no NEWS text (dec-B128).
- `resolveSamplerSpec`: rename in a copy of `matchedCall` before redirecting, because `redirectCall` drops
  `leaf.prior`: `names(priorCall)[names(priorCall) == "leaf.prior"] <- "node.prior"`. Give
  `setDefaultsFromFormals` the formals with `priorFormals["node.prior"] <- callFormals["leaf.prior"]`, so the door's
  `leaf.prior` default fills an absent prior and the tombstone formal's NULL default does not.

  The refusal labels "a linear node prior" and "a Gaussian-process node prior" (two blocks) say leaf.
  `priors$node.prior` reads stay.
- `parsePriors`: the formal stays `node.prior`. `resolveSpec` gets the labels "leaf.prior" and "leaf prior", so
  every error names the user's argument. Strings in R/model.R that say "node prior" become "leaf prior":
  - the `dbartsModel` initialize columns refusal and the monotone refusal;
  - [`resolveLeafCovariates`](../../R/model.R), `linear` and `gp`;
  - "give at most one of 'sd' and 'scale' to a node prior";
  - the dbartsPriors comment.

  "no node scale is defined for family" stays with the slot.
- `dbartsSampler` (R/dbarts.R):
  - `sampleNodeParametersFromPrior` becomes `sampleLeafParametersFromPrior`, with the docstring "Draws leaf values
    from their prior; does not change tree structure.", plus the tombstone method.
  - `getCalibration` becomes `getLeafPrior` and `setCalibration` becomes `setLeafPrior`. That covers their
    docstrings, the `refuseCountsMutation` and `refuseAmplitudeMutation` labels, and the prior.mean refusal's hint.
  - The docstring cross-references in `setForestWeights`, `setForestBasis`, `getForestFits`,
    `getForestAmplitudes` and `getForestVariableCounts`.
  - growFromRoot's "linear and gp node priors", and [`samplePriorPredictive`](../../R/dbarts.R).
- Readers of the renamed method: [`packageBartResults`](../../R/bart.R), [`predictForest`](../../R/generics.R) and
  [`predictBlend`](../../R/generics.R). Comments in R/bartcore.R and inst/common/bartcoreHandle.R.
- R/tombstones.R gets two registry entries:
  - `node.prior`: kind "argument", owner dbarts, successor "leaf.prior";
  - `sampleNodeParametersFromPrior`: kind "rcMethod", owner dbartsSampler, successor
    "sampleLeafParametersFromPrior".

src/ (bridge only):
- Registrations in R_interface.cpp and definitions in R_interface_bartcore.hpp/.cpp: `bartcore_getCalibration`,
  `bartcore_setCalibration` and `bartcore_sampleNodeParametersFromPrior` become `..._getLeafPrior`,
  `..._setLeafPrior` and `..._sampleLeafParametersFromPrior`.
- The reader's column names `"node.scale.factor"` and `"node.scale.divisor"`.
- Strings saying "node prior": "a linear or Gaussian-process node prior", "scale of node prior", and the
  leaf-covariate and gp checks.
- Slot reads and "a non-default node scale" stay.

man/, vignettes/, README.md (manual prose follows):
- Usage and arguments: bart.Rd, dbarts.Rd, dbartsSpec.Rd and xbart.Rd get `leaf.prior`. dbarts.Rd
  also gets a trailing `node.prior = NULL`, written like dbarts.Rd's `sigma` item.
- dbartsSampler-class.Rd:
  - aliases, `\S4method` usage and every mention of the three methods and the two columns;
  - an alias and usage line for `sampleNodeParametersFromPrior` beside `startThreads`;
  - "a linear or gp node prior" in the growFromRoot and getTrees paragraphs.
- dbarts-deprecated.Rd: two new items.
- dbarts.Rd: the node prior mentions in the forests, family and prior-scale paragraphs.
- bart.Rd, in the multinomial paragraph:
  - "the usual node prior";
  - "the node prior's node.scale", which becomes "the leaf prior's scale";
  - "tree and node prior objects" and "tree and node priors" in the tree.prior/leaf.prior item;
  - the "End-node prior parameter k" and "end-node k" pointers.
- bartBT.Rd: the "End-node prior parameter k" subsection, "node parameters", "node prior standard deviation" and
  "end-node sensitivity".
- Also: dbartsPriors.Rd ("Normal prior on the node means"), dbarts-package.Rd, the `$setCalibration` row in
  dbarts-embedding.Rd, forest.Rd, samplePriorPredictive.Rd and xbart.Rd.
- Vignettes:
  - dbarts-as-a-component.Rmd;
  - working_with_saved_trees.Rmd;
  - gibbs_sampler_mixture_model.Rmd: the `node.prior` call and `$getCalibration`/`$setCalibration`. The logistic
    sentence becomes "so the leaf-value prior widens to match (its scale becomes pi * sqrt(3) in place of
    probit's 3)", without naming the slot.
- README.md's feature list.
- Where "calibration" names the prior in force, say leaf prior. "Calibration map" stays.
- inst/NEWS.Rd, 1.0-0 section only: every old spelling of the four renames (about 23 lines), plus the NEWS items
  below.

Tests and harnesses:
- 43 files in inst/tinytest:
  - test-active-rows-pins, test-argument-surface (formal lists), test-augmentation, test-bartcore;
  - test-bcf-creation (pins "linear node prior" and "Gaussian-process node prior"), test-bcf-family,
    test-bcf-loglik, test-bcf-r5-surface, test-boundary-inputs;
  - test-calibration-creation, test-calibration-midchain, test-calibration-prior-draws;
  - test-composition-sequences, test-data-handle, test-data-mixed-mutation, test-data-mixed, test-data-sparse;
  - test-embedding-recipes, test-engine-constants, test-family-objects, test-fits-without-offset,
    test-forest-basis-r5, test-gp-leaves, test-grow-from-root;
  - test-heteroscedastic-mutation, test-heteroscedastic, test-host-shell-pins, test-level-fibre,
    test-linear-leaves;
  - test-model-errors (pins "'node.prior' must be a node prior specification"), test-model-priors, test-monotone;
  - test-multinomial-r5-surface, test-mutate-then-serialize, test-predict-blend, test-predict-forest,
    test-predict-sparse;
  - test-sampler-model, test-sampler-prior, test-sampler-updateState, test-sparse-factor, test-spec,
    test-zero-weights.

  Also test-xbart-error's comment and inst/common/leafPriorChecks.R. A `node.prior` passed to `xbart` in
  test-gp-leaves and test-linear-leaves becomes `leaf.prior`. Slot reads stay.
- benchmarks/R: every file except classic-compare.R.
  - Exact gates: aft-exact, backfit-exact, bcf-exact, bcf-exact-weak, bcf-exact-restricted, bd-balance,
    categorical-exact, change-balance, linear-exact, logistic-reference, monotone-reference, negbin-exact,
    ordinal-exact, perturb-balance, rule-gibbs-balance, swap-balance, t-exact.
  - Scripts outside any gate: sbc.R and geweke-mc.R (both call `$getCalibration`), binary-hyperprior,
    composition-matrix, constant-gp-max-leaf-size, constant-linear-leaf-covariates, memory-footprint.
  - equivalence.R.
- tests/cpp, tools/, .github/ and _pkgdown.yml: nothing.

docs/:
- Update these standing references:
  - docs/architecture.md: "the tree and node prior objects";
  - design/public-surface.md: its prior-DSL and sampler-method prose;
  - design/INDEX.md: the linear, gp and nameable-calibration summaries;
  - design/prior-defaults.md, design/feature-matrix.md and design/bart-as-a-component.md;
  - plans/bartcore-landing/changes.md: the door's argument list;
  - plans/bartcore-review-tour.md.
- Leave the other live-tree docs with hits; they are records.
- Broken cites: design/multinomial-mutation-arc.md cites `getCalibration` in R/dbarts.R and the two bridge
  functions. Convert each to the history form, pinned at 38a9c7877ad60307bf98aa5b8ca22d0c739bd2f2.

## Tombstone mechanics

- `resolveRenamedLeafPrior(matchedCall, nodePriorSupplied, leafPriorSupplied, caller)` in R/tombstones.R returns
  the call.
  - `node.prior` absent: the call is returned unchanged.
  - Both present: `stop("'node.prior' and 'leaf.prior' name the same prior on '<caller>'; supply one", call. =
    FALSE)`, even when the two agree.
  - Otherwise: `warnOnce(paste0("tombstone.node.prior.", caller), "'node.prior' is now 'leaf.prior' on
    '<caller>'; the value was used. The old name is removed in dbarts 1.1-0.")`, then
    `matchedCall$leaf.prior <- matchedCall$node.prior; matchedCall$node.prior <- NULL`.

  The helper never forces the argument, because the prior vocabulary is NSE. Placement: in `dbarts`,
  right after `resolveConsolidatedArgs`, beside the `sigma` resolution. Read `missing()` of both
  formals before either is assigned. The returned call is used from then on, so the stored call carries
  `leaf.prior`.
- Method `sampleNodeParametersFromPrior = function(updateState = NA)`: `warnOnce(
  "tombstone.sampleNodeParametersFromPrior", "'$sampleNodeParametersFromPrior' is now
  '$sampleLeafParametersFromPrior'; this call was forwarded. The old name is removed in dbarts 1.1-0.")`, then
  `sampleLeafParametersFromPrior(updateState)`.
- There is no tombstone for `getCalibration`/`setCalibration` or for the columns. An old call fails with R's "not a
  valid field or method name".

## Steps

1. dbarts commit 1, "Say leaf where a name means a leaf": the R/, src/, man/, vignettes/, README.md, NEWS,
   tinytest, inst/common and benchmarks edits above. test-tombstones.R gains these checks:
   - `node.prior` on `dbarts` warns exactly once. Reset the key first and count with `withCallingHandlers`. The
     fit is bitwise the one `leaf.prior` gives under one seed.
   - `node.prior` on `dbartsSpec` is refused, the message naming `leaf.prior`.
   - Supplying both spellings is the error.
   - `bart(..., node.prior = )` and `xbart(..., node.prior = )` fail with "unused argument 'node.prior'" and name
     `leaf.prior`.
   - `$sampleNodeParametersFromPrior` warns once and leaves the state `$sampleLeafParametersFromPrior` leaves.
   - No leaks: default `bart`, `bartBT`, `xbart`, `dbarts`, `dbartsSpec` and a multinomial `bart` emit zero
     warnings.

   Install with `--preclean`, since the bridge registrations change.
2. dbarts commit 2 (docs only): the standing-reference updates and the multinomial-mutation-arc.md cite repair.
3. Residual prose sweep, before step 1 is reviewed:
   `git grep -n -i -E 'node prior|end-node|end node|node mean|node param' -- man vignettes README.md R src/*.cpp src/*.hpp`.
   The allowlist:
   - R comments and bridge strings about the kept slots ("node scale", "no node scale is defined",
     "a non-default node scale");
   - parsePriors' `node.prior` plumbing;
   - the tombstone and deprecation text.

   Every other hit becomes leaf.
4. stan4bart (branch bartcore). `dbartsSpec` refuses `node.prior`, so a `bart_args$node.prior` fails until stan4bart
   spells it `leaf.prior`. mvbart's `k` path calls `dbarts()`, which keeps the tombstone, so it warns
   too. The changes:
   - R/stan4bart_fit.R: the `k` shorthand writes `spec_call[["leaf.prior"]]`, and its collision check tests
     `leaf.prior` and `node.prior` with the message "bart_args cannot set both 'k' and 'leaf.prior'". Also the
     "end node priors" comment near the chi hyperprior binding.
   - R/mvbart.R: `dbCall$node.prior` becomes `dbCall$leaf.prior`.
   - man/stan4bart.Rd (three mentions) and man/mvbart.Rd ("the end-node prior scale").
   - NEWS.md: the unreleased 0.0-14 bart_args item.
   - inst/tinytest/test-09-bartArgs.R: three calls and the pinned message.
   - docs/plans/dbarts-spec-adoption.md is a record; leave it.
5. bartCause (branch dbarts-1.0): in R/bcf.R, `sampler$getCalibration(1L)` becomes `sampler$getLeafPrior(1L)`.
   Without it bartCause breaks, since there is no tombstone. bcf passes either prior spelling through `dbarts()`'s
   formals.
6. treatSens and bairrtt: no commit. treatSens's named `parsePriors(node.prior = )`, its `priors$node.prior` and
   `new("dbartsModel", node.scale = )` keep working because those internals keep their names. bairrtt has no site.
7. Release contact: add countSTAR to the maintainer-contact list in TODO's release block. Its CRAN release calls
   `dbarts(..., node.prior = normal(k))`, which the critic reported and which was not re-checked here. The
   tombstone keeps it working until 1.1-0 with a once-per-session warning; ask for `leaf.prior`, version-guarded
   while it supports 0.9-x.
8. Push order: dbarts and bartCause together, since bartCause breaks against either side alone, then stan4bart.
   The records commit carries this plan's Landing note and removes the TODO entry.

## Verification

Against the slice's own library (`R_LIBS=$LIB` on every R call), each gate on its own exit status:

```sh
R CMD INSTALL --preclean -l $LIB .
(cd tests/cpp && make && ./test_bartcore)
Rscript -e 'lintr::lint_package()'
air format --check .
Rscript tools/check-rc-codoc.R .
Rscript tools/check-win-drift.R .
Rscript tools/check-doc-freshness.R .
Rscript -e 'db <- tools:::.build_news_db_from_package_NEWS_Rd("inst/NEWS.Rd"); stopifnot(!is.null(db)); print(nrow(db))'
Rscript -e 'tinytest::test_package("dbarts")'
for f in sbc geweke-mc binary-hyperprior composition-matrix constant-gp-max-leaf-size \
  constant-linear-leaf-covariates memory-footprint; do Rscript -e "invisible(parse('benchmarks/R/$f.R'))" || exit 1; done
```

- NEWS entry count is the pre-change count + 1.
- `R CMD check --as-cran` on a tarball built from a clean copy staged outside the tree: no errors or warnings, and
  notes unchanged.
- Smoke run of the `$getLeafPrior` paths: sbc.R at its smallest documented setting, and geweke-mc.R's calibration
  check (`chainSampler$getLeafPrior()` against the generator's). Each runs without error, and geweke's identity
  holds.
- Bitwise equivalence. Build a reference install with `R CMD INSTALL --preclean
  --configure-args=--enable-reference-build -l $REFLIB .`, then under `R_LIBS=$REFLIB` run:
  - `benchmarks/R/equivalence.R compare benchmarks/baselines/equivalence-3900c989.rds --bitwise --strict-coverage`
  - `bcf-equivalence.R compare benchmarks/baselines/bcf-equivalence-d49e2103.rds --bitwise`
  - `multinomial-equivalence.R compare benchmarks/baselines/multinomial-equivalence-80b1c8d4.rds --bitwise`

  Count the per-scenario "identical draws (same RNG stream)" lines against the scenario count. There must be no
  "max |z|" line.
- The exact-gates.yaml gate loop in `quick` mode, because its scripts are edited.
- Name sweep:
  `git grep -n -E 'node\.prior|sampleNodeParametersFromPrior|getCalibration|setCalibration|node\.scale\.(factor|divisor)' -- R src man vignettes inst README.md benchmarks`
  may list only these:
  - the two tombstone formals and the helper;
  - the tombstone method and the registry;
  - the hint map;
  - parsePriors' formal, its list names and the `@node.prior` reads;
  - the bridge's slot strings;
  - dbarts.Rd's tombstone entry, and dbarts-deprecated.Rd;
  - the dbartsSampler-class.Rd alias;
  - NEWS's tombstone list and 0.9-x history;
  - test-tombstones.R and classic-compare.R;
  - src/bartcore identifiers.
- Sister packages, each installed into `$LIB` beside the new dbarts:
  - stan4bart: `tinytest::test_package("stan4bart")`.
  - bartCause and treatSens: `testthat::test_dir("tests/testthat", package = "<pkg>", load_package = "installed")`
    from each checkout. treatSens is checked unchanged.
  - bairrtt: `tinytest::test_package("bairrtt")`, as a smoke test.

  Every suite passes, with no "is now 'leaf.prior'" warning and no `sampleNodeParametersFromPrior` warning in its
  output.

## NEWS draft

UPGRADING, after the `sigest` rename item:

```
      \item \code{dbarts(node.prior = )}
            is renamed \code{leaf.prior}, the name \code{bart},
            \code{xbart} and \code{dbartsSpec} take it under too, and the sampler's
            \code{$sampleNodeParametersFromPrior} is
            \code{$sampleLeafParametersFromPrior}; both old names are
            still accepted for one release (tombstone, above). A name now
            says leaf where it means only the trees' terminal nodes and
            node where it means every node: \code{$getTrees}' rows,
            \code{plotTree}'s \code{nodeWidth}, \code{nodeHeight} and
            \code{nodeGap}, and the tree prior's per-node split
            probability keep node.
```

Tombstone list, after the `sigma` entry: `\code{node.prior} on \code{dbarts} (successor
\code{leaf.prior}); \code{dbartsSampler}'s \code{$sampleNodeParametersFromPrior} (successor
\code{$sampleLeafParametersFromPrior});`

## Calls made

1. Tombstone reach: `dbarts` only. `bart`, `xbart` and `dbartsSpec` refuse `node.prior` and name `leaf.prior`. Rejected: tombstones on all four doors, because dec-A02 and dec-B128
   cover only shipped names.
2. The internal plumbing and the model keep node: `parsePriors`' formal and list names, the `dbartsModel` slots and
   the `dbartsNodePrior` classes. `dbartsModel` is not exported and its slots are undocumented. 0.9-x fits are
   refused by `refuseLegacyState` anyway. treatSens and CRAN stan4bart 0.0-13 build `dbartsModel` with
   `node.scale =` and read `priors$node.prior`. Rejected: renaming them, which breaks those consumers and misreads
   saved objects.
3. Both spellings supplied is an error, following the `sigma` precedent. Rejected: letting `leaf.prior` win.
4. The bridge's three registrations are renamed to match the R methods; engine identifiers wait for the
   host-neutral work. Rejected: keeping the bridge names, which would break the 1:1 naming between R methods and
   bridge functions. A sampler saved from a pre-rename branch build, with the old method already cached, fails that
   call; nothing released is affected.
5. The manual says leaf prior where "calibration" named the prior in force; "calibration map" stays. Rejected:
   keeping calibration as the noun (dec-B129 notes no package uses it).
6. treatSens needs no change, because item 2 keeps its internals working. Rejected: a lockstep treatSens edit.
   Moving treatSens onto `dbartsSpec()` stays a separate treatSens item.
7. [`resolvePriorScale`](../../R/model.R)'s own `node.prior`/`node.hyperprior` formals stay: they are a private
   helper called positionally from [`resolveSamplerSpec`](../../R/spec.R) and [`xbart`](../../R/xbart.R) (the
   latter through its renamed `leafPrior` local), so the parameter spelling is invisible to every caller. Kept in
   the same "internal plumbing keeps node" bucket as item 2, rather than renamed for cosmetic symmetry. Rejected:
   renaming the two formals, an unforced diff with nothing observing it.

## Landing

dbarts: two commits, `Say leaf where a name means a leaf` (59d0dbbb) and
`Update standing docs for the leaf-vocabulary rename` (beb60ae1). stan4bart
(branch bartcore): `Follow dbarts's node.prior -> leaf.prior rename`
(42bfdaf). bartCause (branch dbarts-1.0): `Call the sampler's renamed
getLeafPrior` (1c497f4).

Verification: tinytest 9336/9336; `tests/cpp/test_bartcore` all green;
`lintr::lint_package()` and `air format --check .` clean; `check-rc-codoc.R`,
`check-win-drift.R` and `check-doc-freshness.R` all OK; `inst/NEWS.Rd` parses
(323 entries); the 25-gate `exact-gates.yaml` battery in `quick` mode all
PASS; `sbc.R discrete-selfcheck` and `geweke-mc.R quick` both OK, the
[`$getLeafPrior`](../../R/dbarts.R) paths exercised by each; `R CMD check
--as-cran` on a tarball built from a clean staged copy: no errors, warnings
or new notes. Bitwise equivalence against the MANIFEST baselines, reference
build: `equivalence-3900c989.rds` 53/53, `bcf-equivalence-d49e2103.rds`
15/15, `multinomial-equivalence-80b1c8d4.rds` 11/11, every scenario
"identical draws (same RNG stream)", no `max |z|` line - confirms the
neutral classification. stan4bart: 565/565 tinytest expectations.
bartCause: `testthat::test_dir` clean (no failures, the usual On-CRAN
skips only). Name sweep and the Step 3 prose sweep both clean against
their allowlists.

test-host-shell-pins.R's method census bumped from 50/45 to 51/46 for the
new `sampleNodeParametersFromPrior` tombstone method - not itself a plan
item, but the census is a drift detector and would otherwise fail on any
method count change.
