# front-door

agent: opus (S1 the rename, shim, legacy door and tombstones; S2 family
  objects, consolidation and na.action); sonnet (S3 xbart; S4 manual,
  NEWS, pkgdown, consumer ports). Serial: S1, S2, S3, S4. S2 and S3
  both edit R/xbart.R (the consolidation removes xbart's `dart` and
  family-only formals), so they do not run in parallel.
rng: neutral for S1, S2 and S4 (the same arguments reach `dbarts()` in the
  same order). S3 is shifting:
  the xbart worker redistribution and per-unit seeds move xbart's stream,
  so its snapshot file is replayed and the two xbart equivalence
  scenarios re-record with the fold oracle
  ([test-xbart-fold-oracle.R](../../inst/tinytest/test-xbart-fold-oracle.R))
  as the oracle the baseline manifest requires of any re-record that
  moves a draw (its rule P17). The other three slices expect the three
  bitwise equivalence baselines (gaussian, BCF, multinomial) IDENTICAL.
window: R-only; runs beside the header arc (docs/plans/pure-c-header.md).
  Names lock at the CRAN push (dec-B79), so nothing here needs a
  tombstone for a name chosen in this arc. Settle before lorax's
  maintainer is contacted (TODO, release block). The interfaces arc
  (docs/plans/interfaces-and-dependencies.md) rebases onto S2, which
  rewrites the model.frame call it also edits; the engine arc
  (docs/plans/engine-performance.md) owns any new `dbartsControl`
  formal (its within-chain threading argument), and this arc's
  consolidation does not touch `dbartsControl`'s list.
budget: S1 ~350 R (bart body moves, legacy door ~120, shim ~60,
  tombstones ~80) + ~400 tests; S2 ~300 R + ~350 tests; S3 ~150 R + ~120
  tests; S4 ~500 Rd/NEWS + consumer ports (bartCause ~20, stan4bart ~5).

Decisions, all in [docs/decisions.md](../decisions.md): dec-B75 (bart is
the modern door), B76 (tombstones), B77 and B102 (xbart), B79 (names
lock at the CRAN push), B80 (sigest), B81 (two family vocabularies),
B82 (hurdle tokens), B83 (the legacy door is strict), B98 and B101
(family objects, consolidation, residual law), B108 (na.action), B109
(saved 0.9-x fit). dec-B78 (factor default) is already the shipped
behaviour and only gains its manual paragraph here.

## Goal

`bart` is the package's front door with today's `bart2` interface and
defaults. The BayesTree-style function lives under a new name with
0.9-34's `bart` formals exactly, refusing a factor response of three or
more levels with a two-remedy message. For one release `bart2` is an
alias and a BayesTree-spelled call to `bart` forwards to the legacy door
with a once-per-session warning. Removed functions and arguments are
tombstones that name their successor and expire together. Family
settings ride family objects; the residual law is a family; the front
door has a standard `na.action` whose default drops missing-response
rows only. `xbart` cross-validates with k fixed or modelled, distributes
work over folds, and reproduces a seed at any thread count. A saved
0.9-x fit is refused by name.

## Context

- The two doors today: [`bart`](../../R/bart.R) with BayesTree's formals
  plus bartcore-era extras (`family` of three tokens, `resid.dist`,
  `subset`, `storage`, `prior.scale`) and [`bart2`](../../R/bart.R) with
  the modern formals, thirteen family tokens, and the feature objects.
  Neither forwards to the other; both build priors with
  [`buildSamplerPriors`](../../R/bart.R) and reach
  [`dbarts`](../../R/dbarts.R), bart2 through
  [`buildHostSamplerCall`](../../R/bart.R). 0.9-34's `bart` formals are
  the list in this plan's legacy-door step; they carry no `family`,
  `resid.dist`, `subset`, `storage` or `prior.scale`.
- Family-gated formals: [`familyGatingInventory`](../../R/utility.R)
  already names them (sigest, sigdf, sigquant, resid.prior for the
  gaussian-scale families; dispersion for nbinom; breaks and max.rows
  for the hazard trio). `twopart` folds to `hurdle.lognormal` at both
  doors; the composition is [`bart2Hurdle`](../../R/bart.R); `dbarts()`
  lists the hurdle tokens and refuses them.
- Residual law: [`student`](../../R/model.R), [`gaussian`](../../R/model.R)
  and the unexported [`dbartsResidDists`](../../R/model.R); the spec
  refuses Student-t outside the continuous gaussian response
  (["student residuals require a continuous gaussian response"](../../R/spec.R)),
  and the header already carries
  [`DBARTS_FAMILY_STUDENT`](../../inst/include/dbarts/dbarts.h).
- Family objects: none exist. The prior constructors
  ([`chi`](../../R/model.R), [`chisq`](../../R/model.R),
  [`dart`](../../R/model.R), [`gp`](../../R/model.R)) and the feature
  objects ([`interactions`](../../R/model.R), [`blocks`](../../R/model.R),
  [`varianceForest`](../../R/model.R)) are the house pattern to copy:
  a constructor returning a validated S4 object, listed in
  [`dbartsPriors`](../../R/model.R).
- NA handling: `stats::na.pass` is forced at every model.frame call
  ([`dbartsData`](../../R/data.R) and the multinomial and forest-term
  sites); the response is then refused
  (["response contains missing values"](../../R/data.R)); predictors go
  to the `missing` argument (incorporate or error). No `na.action`
  formal exists.
- sigma vs sigest: `dbarts(sigma = )` and `bart(sigest = )`, validated
  together in [`validateArgumentsInEnvironment`](../../R/dbarts.R); the
  sampler method [`dbartsSampler$setSigma`](../../man/dbartsSampler-class.Rd)
  is untouched.
- xbart: [`xbart`](../../R/xbart.R) chunks replications across workers
  ([`xbartRunChunk`](../../R/xbart.R) takes a replication range and one
  seed), so k-fold with one replication runs on one worker (2.1x slower
  wall than 0.9-34 at four threads, measured 2026-09-08); k is a numeric
  grid or, when absent, the fixed value 2 or the node prior's k; a
  hyperprior is held, not swept. `control=` already has no formal; a
  three-element n.burn errors by length.
- Tombstone infrastructure: none. [`warnOnce`](../../R/utility.R) and
  `onceWarnState` are the once-per-session primitive. rbart_vi,
  rngSeed, startThreads and stopThreads survive only in NEWS.
- Saved state: the `formatVersion` attribute is read only by the bridge
  ([`stateFormatVersion`](../../src/R_interface_bartcore.cpp)); the R
  restore path is [`dbartsSampler$getPointer`](../../man/dbartsSampler-class.Rd)
  and `setState`, and [`predict.bart`](../../R/generics.R) reaches it
  through `object$fit`. A 0.9-x fit's state carries no such attribute.
- Consumers: bartCause's dbarts-1.0 branch reads `formals(dbarts::bart2)`
  in three files and redirects calls to `dbarts::bart2`; stan4bart's
  bartcore branch filters arguments by dbarts formals; 17 CRAN reverse
  dependencies call `bart` in BayesTree spelling (surface-refusals plan,
  section 14).
- Prior records amended by this arc: public-surface design sections
  [3. Response families](../design/public-surface.md#3-response-families)
  and [3a. Prior specification](../design/public-surface.md#3a-prior-specification);
  robust-errors design
  [8. Resolution (VD, 2026-07-18)](../design/robust-errors.md#8-resolution-vd-2026-07-18);
  prerc-surface-freeze's D5 deprecation-shims ruling; hurdle design
  section 13; negative-binomial section 4; survival section 2.

## Decision

Open forks, each with a recommendation; they go to VD one at a time
before S1 starts.

1. The legacy door's name: `bartBT` (VD 2026-09-08, chosen over
   `bart_bt` and `bartCompat`).
2. A `bart(x, y)` call with no BayesTree-spelled argument (VD
   2026-09-08): modern defaults silently, a positional guard (during
   the transition release a fourth or later positional argument to
   `bart` stops with a message naming both doors, since a 0.9-x call
   would bind it to `subset` where 0.9-34 read `sigest`), and a
   package-load message in R's standard form (`packageStartupMessage`
   from `.onAttach`, so `suppressPackageStartupMessages` silences it)
   saying `bart` is the modern door and `bartBT` the BayesTree-style
   one, for one release. The shim fires on BayesTree spellings only.
3. The legacy door's extras (`family`, `resid.dist`, `subset`,
   `storage`, `prior.scale`, none of which shipped on CRAN): all five
   dropped, so `bartBT` is 0.9-34's argument list exactly (VD
   2026-09-08: "No one wrote any scripts against bartcore's branch.
   Drop all five.").
4. The Student-t spelling. Since Student-t errors exist only for the
   continuous gaussian response, the residual-law setting has one home.
   Options: (a) `family = student(df)` as its own family object and
   token, matching the engine's own family list; (b)
   `family = gaussian(errors = student(df))`, a setting on gaussian.
   Recommended: (a); `resid.dist` and `dbartsResidDists` retire.
5. The `na.action` default's name. It behaves like `na.exclude` on the
   response column only. Candidates: `na.excludeResponse` (exact,
   long), `na.dbarts` (short, says nothing). Recommended:
   `na.excludeResponse`.
6. xbart's "leave k modelled" spelling. Options: (a) `k` accepts a list
   whose entries are numbers or hyperprior objects, so
   `k = list(1, 2, chi())` sweeps two fixed values and the modelled
   default; (b) `NA` in the numeric grid means modelled. Recommended:
   (a); an NA token is obscure and cannot name the hyperprior.

Fixed by the decisions: the consolidation's scope (dec-B98 says any
remaining family-only or feature-only formal moves onto its object, so
`breaks`, `max.rows`, `dispersion` and `resid.dist` go to family
objects and the two feature duplicates `dart`, a logical that
duplicates `tree.prior = dart()`, and `levelGibbs`, a categorical-split
setting, go to the tree prior; `monotone`, `interactions`, `blocks` and
`variance` already take objects and stay as the slots those objects
fill; `factors`, `missing`, `warm.start` and `n.grow.sweeps` are data
and start settings and stay); the once-per-session warning; one-release
expiry for every tombstone; sigest everywhere; the factor refusal
message; no conversion of 0.9-x fits.

## Constraints

- Gates per slice as above. tinytest whole suite; `R CMD check --as-cran`
  from a clean tarball (R/ and man/ change); `lintr::lint_package()`
  (names move); pkgdown index check for every new topic; NEWS parses.
- `formals(dbarts::bart2)` must stay a real formal list during the alias
  release (bartCause evaluates its defaults), so the alias copies bart's
  formals and forwards a matched call rather than being `function(...)`.
- Every tombstone lives in one file (R/tombstones.R) with its expiry
  version, and one test asserts the registry matches the exports and
  the NEWS list, so the release after 1.0-0 removes them by deleting the
  file.
- The legacy door and the modern door share `dbarts()`; the legacy door
  keeps indicator expansion for factors (dec-B78) and `missing = "error"`.
- Out of scope: any engine or bridge change; the survival formula
  interface and sparse columns (docs/plans/interfaces-and-dependencies.md);
  the rbart_vi fallback of dec-B105 (a stan4bart condition, tracked in
  TODO); the factor split-mass weighting research (dec-B78).

## Steps

S1, rename, shim, legacy door, tombstones:

1. Move `bart2`'s body and formals to `bart`; `bart2` becomes the alias
   (same formals, forwards `match.call()` with the function replaced,
   once-per-session message naming `bart`, expiry 1.1-0).
2. Create `bartBT` with 0.9-34's formals:
   `x.train, y.train, x.test, sigest, sigdf, sigquant, k, power, base,
   splitprobs, binaryOffset, weights, ntree, ndpost, nskip, printevery,
   keepevery, keeptrainfits, usequants, numcut, printcutoffs, verbose,
   nchain, nthread, combinechains, keeptrees, keepcall, sampleronly,
   seed, proposalprobs, keepsampler`. Its body is today's `bart` body
   minus the extras; binary responses go to probit; a factor response of
   three or more levels stops with the message naming `bart(family =
   "multinomial")` and `as.integer(y) - 1L`. Its result object is class
   `bart` as today so `predict`, `extract` and friends work unchanged.
3. The shim: `bart`, `dbartsControl` and `xbart` gain `...` for the
   shim release, since none has it today and R refuses an unknown name
   before the body runs; a name that is neither a formal nor a registry
   entry is refused by the derived-refusal mechanism the surface
   refusals plan already uses ([`foreignArgsFor`](../../R/generics.R)).
   `bart` inspects the supplied argument names; any name from
   the legacy door's formals that is not a modern formal (`x.train`,
   `y.train`, `x.test`, `ntree`, `ndpost`, `nskip`, `keeptrees`,
   `nchain`, `nthread`, `numcut`, `usequants`, `keepevery`,
   `printevery`, `keeptrainfits`, `printcutoffs`, `combinechains`,
   `keepcall`, `sampleronly`, `keepsampler`, `splitprobs`,
   `proposalprobs`, `binaryOffset`) forwards the whole call to the
   legacy door after `warnOnce` names it; a fourth or later positional
   argument stops with a message naming both doors; `.onAttach` emits a
   `packageStartupMessage` naming both doors. All three expire at 1.1-0
   (registry entries).
4. Tombstones (R/tombstones.R): `rbart_vi` errors naming stan4bart and
   saying the group-spread prior differs there so results move
   (dec-B105); `predict`, `extract`, `fitted`, `residuals` methods for
   class `rbart` say the same; `dbartsSampler$startThreads` and
   `stopThreads` are no-op methods; `dbartsControl(rngSeed = )` and
   `bart(rngSeed = )` are accepted as `seed` with a once-per-session
   warning; `family = "twopart"` errors naming `hurdle.lognormal`. The
   registry lists each with its expiry; the test asserts the registry,
   NAMESPACE and NEWS agree.
5. sigest everywhere (dec-B80): `dbarts(sigma = )` and
   [`dbartsSpec`](../../R/spec.R)`(sigma = )` become `sigest`; `sigma`
   accepted with a once-per-session warning for one release (registry
   entry). `validateArgumentsInEnvironment` validates one name. The
   `dbartsData` class's `sigma` slot keeps its name: a slot is not an
   entry point, and stan4bart writes it directly on its bartcore
   branch.
6. Hurdle tokens (dec-B82): `hurdle.lognormal` and `twopart` leave
   `dbarts()`'s family list; `bart` intercepts `hurdle.lognormal` before
   forwarding and calls the composition; `dbarts()` refuses the token by
   name pointing at `bart`.
7. Saved 0.9-x fit (dec-B109): in `getPointer`'s re-creation branch and
   in `setState`, a non-NULL state with no `formatVersion` attribute
   stops with "this fit was saved by dbarts 0.9-x (no state format
   field); dbarts 1.0-0 cannot read it; refit with this version".
   `predict.bart` reaches it through `object$fit`; a 0.9-x fit saved
   without trees keeps today's keepTrees refusal. Test with a fixture
   built by hand (a list with the 0.9-x field names and no attribute).

S2, family objects, consolidation, na.action:

8. Family objects: a `dbartsFamily` S4 class with a token slot and a
   settings list; constructors `gaussian()`, `student(df = NULL)`,
   `probit()`, `logistic()`, `multinomial()`, `ordinal()`,
   `nbinom(dispersion = NA)`, `aft()`, `hazard(breaks = NULL, max.rows =
   1e7, link = c("probit", "logistic"))`, `hurdle.lognormal()`. `family`
   accepts a token string or an object; a string resolves to the
   object with defaults. The bridge-facing resolution stays in
   [`dbartsSpec`](../../R/spec.R). One exported list `dbartsFamilies`
   mirrors `dbartsPriors`. `gaussian()` and `student()` keep working as
   `resid.dist` values only through the tombstone.
9. Consolidation: `breaks`, `max.rows`, `dispersion`, `resid.dist`,
   `dart` and `levelGibbs` leave `bart`, `dbarts` and `xbart` with
   registry tombstones that name the object and argument. The
   `familyGatingInventory` warning shrinks to the four gaussian-scale
   formals. The ten-million row cap moves to `hazard()`'s `max.rows`
   default and is measured in the constants audit
   (docs/plans/engine-performance.md).
10. na.action (dec-B108): formal on `bart`, `dbarts` and `dbartsData`,
    default the package function of fork 5. Ingestion passes the
    function to `model.frame` at the main site in `dbartsData` and
    applies it to the (y, x) pair on the matrix path; the other four
    `na.pass` sites (the test frame in R/data.R, the multinomial frame
    in R/bart.R, the forest-term bases in R/formulaTerms.R and the
    basis site in R/model.R) keep `na.pass` and align to the rows the
    main frame kept; the package default drops rows
    with a missing response, keeps missing predictors for `missing`,
    and sets a `na.action` attribute of class `exclude` so `fitted` and
    the training-fit slots pad to the data's length through
    `stats::naresid`. `na.omit`, `na.exclude`, `na.fail`, `na.pass` keep
    their base meaning; under `na.pass` the response check errors as
    today. Test data are unchanged. Tests cover each function on the
    formula and matrix paths, the padding length, and insight's
    NA-response case.

S3, xbart (dec-B77, B102):

11. Distribute (replication, fold) units across workers instead of
    replication ranges; cells run in fixed order inside a unit so the
    k warm start within a fold is unchanged; warm starts across folds
    stay out.
12. Each unit draws its seed from the call's seed and its (replication,
    fold) index, so a seed reproduces at any `n.threads`; test at 1 and
    4 threads.
13. k axis per fork 6: the grid holds numbers and hyperprior objects;
    a hyperprior cell hands the object to the node prior and reports
    "modelled" in the result's dimnames. Tombstones for `control=` and a
    three-element `n.burn` name the flat fields.
14. Replay the xbart snapshot file; re-record the two xbart equivalence
    scenarios with the fold oracle as the P17 row.

S4, manual, records, consumers:

15. `bart.Rd` documents the modern door; the legacy door gets its own
    page carrying the BayesTree defaults, the indicator-factor note
    with why (dec-B78), the factor-response refusal and the "no new
    features" statement; `bart2` and every tombstone share one
    `dbarts-deprecated.Rd`. A `dbartsFamilies.Rd` page documents the
    family objects and holds the one mapping table (dec-B81): front-door
    token, family object, engine family, entry points that accept it.
    `na.excludeResponse` (or the chosen name) gets its page. pkgdown
    reference entries for each.
16. NEWS 1.0-0 UPGRADING: the rename, the shim, the legacy door, the
    tombstone list with expiry, sigest, family objects, na.action, the
    0.9-x fit message, xbart's k grid.
17. Consumer ports on their branches: bartCause dbarts-1.0 reads
    `formals(dbarts::bart)` and redirects to `dbarts::bart`; stan4bart
    filters on `dbartsControl` and `dbartsSpec` formals, so the rename
    costs it nothing and the `dbartsSpec` sigest rename is its one
    edit. The 17 BayesTree-spelled reverse
    dependencies are left to the shim (one warning each) and named in
    the release block's revdep sweep.
18. Records: the design sections named in Context get one paragraph
    each; feature-matrix rows that name `bart2` move to `bart`.

## Verification

```
R CMD INSTALL -l <lib> .
R_LIBS=<lib> Rscript -e 'tinytest::test_package("dbarts")'
R_LIBS=<lib> Rscript benchmarks/R/equivalence.R compare benchmarks/baselines/equivalence-fbff1989.rds
  # S1, S2, S4: 50 identical; S3: 48 identical, xbart and xbartmixed re-recorded
R_LIBS=<lib> Rscript -e 'lintr::lint_package()'
R CMD build <clean copy> && R CMD check --as-cran dbarts_*.tar.gz
R_LIBS=<lib> Rscript -e 'pkgdown::check_pkgdown(".")'
R_LIBS=<lib> Rscript -e 'tools:::.build_news_db_from_package_NEWS_Rd("inst/NEWS.Rd")'
```

Expected: `bart(y ~ ., df)` and `bart2(y ~ ., df)` give identical draws
at one seed; `bart(x.train = x, y.train = y)` gives the legacy door's
draws and one warning per session; the legacy door at 0.9-34's defaults
reproduces today's `bart` draws for gaussian and probit; a hand-built
0.9-x fit object fails `predict` with the version message.
