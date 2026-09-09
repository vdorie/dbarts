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

All six forks were put to VD on 2026-09-08 and are recorded here with
the choice; none remains open.

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
4. The Student-t spelling: `family = student(df)` as its own family
   object and token, matching the engine's family list; `resid.dist`
   and `dbartsResidDists` retire (VD 2026-09-08, "Use your
   recommendation", over a `gaussian(errors = student(df))` setting).
5. The `na.action` default's name: `na.keepPredictors` (VD
   2026-09-08, chosen for saying what differs from `na.omit`; rpart's
   `na.rpart` is the behavioural precedent and the help page says so).
   `missing` is retired outright, no tombstone (VD 2026-09-08): missing
   predictors are always modelled for the rows `na.action` keeps;
   `na.fail` covers the old `missing = "error"`.
6. xbart's "leave k modelled" spelling (VD 2026-09-08): `k` accepts a
   numeric vector as today or a list whose entries are numbers or
   hyperprior objects, so `k = list(1, 2, chi())` sweeps two fixed
   values and the modelled default; the modelled cell is labelled by
   its constructor call. An absent `k` runs one cell at the front
   door's default for the response type (fixed 2 continuous,
   `chi(1.5, 2)` binary), so a default `xbart` call scores the model a
   default `bart` call fits; this supersedes dec-B12's fixed-2 binary
   default, and a supplied `node.prior` with its own k still wins.

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
  keeps indicator expansion for factors (dec-B78) and 0.9-34's `na.omit` row rule.
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
    default `na.keepPredictors`; the `missing` formal leaves `bart`,
    `dbarts`, `dbartsData` and `xbart` with no tombstone, and the
    incorporate path is always on for kept rows (`bartBT` keeps
    0.9-34's `na.omit` behaviour). Ingestion passes the
    function to `model.frame` at the main site in `dbartsData` and
    applies it to the (y, x) pair on the matrix path; the other four
    `na.pass` sites (the test frame in R/data.R, the multinomial frame
    in R/bart.R, the forest-term bases in R/formulaTerms.R and the
    basis site in R/model.R) keep `na.pass` and align to the rows the
    main frame kept; the package default drops rows
    with a missing response, keeps missing predictors for the trees,
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
    the constructor call in the result's dimnames; an absent `k` is one
    cell at the response type's front-door default (the binary default
    cell therefore moves from fixed 2 to `chi(1.5, 2)`, a draw change
    covered by the xbart re-record). Tombstones for `control=` and a
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
    `na.keepPredictors` gets its page. pkgdown
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

## Landing note, S1 (2026-09-09)

LANDED at ecb319aa00a7d3d70a80131b11c24a081dcbcdc2, eight commits:

- 417a7a41f1b42a399f5d06df947872901fe4908b Rename the BayesTree-style door to bartBT and make bart the modern one
- 533cba8e8084dabb98338e834d16b7896e7c621c Add the tombstone registry, the bart2 alias and the startup message
- 282d18044b4f85539f7ec94e55052f00050ec05c Spell the creation-time estimate sigest, retire the hurdle aliases, refuse a 0.9-x state
- a01436cee459adac84c25a9c843b307a44860146 Update the test suite for the two doors and the tombstones
- 94588cc3e5e88f4d0a5e95a362f1f8fe85d5d778 Update the manual and the reference index for the rename
- 215925aaff9fab43a86ebfa372f128cc0c749ba9 Apply air formatting
- 6a393a74fd30bea3a0853cc2a1781bcc40180b05 Repoint the design cites the front-door rename moved
- ecb319aa00a7d3d70a80131b11c24a081dcbcdc2 Address review: hurdle component calls, the legacy door's formula-path factor refusal, the 0.9-x state shape test, and the remedy call forms

[`bart`](../../R/bart.R) carries the former `bart2` formals and body;
[`bart2`](../../R/tombstones.R) is an alias, real formals, forwarding a
matched call, once-per-session message, expiry 1.1-0. [`bartBT`](../../R/bart.R)
has 0.9-34's 31 formals byte-identical in name, order and default
(checked against the CRAN 0.9-34 tarball in review), none of the five
bartcore-era extras, result class `bart`, and refuses a 3+-level factor
response on both paths with the two-remedy message. The shim forwards
22 BayesTree spellings to `bartBT` after `warnOnce`; a fourth-or-later
positional `bart` argument stops naming both doors; `.onAttach` emits a
`packageStartupMessage`. `...` on `bart`, `dbartsControl` and `xbart`
still refuses an unknown name through [`foreignArgsFor`](../../R/generics.R).
[`dbartsTombstones`](../../R/tombstones.R) holds 16 entries at expiry
1.1-0 (`rbart_vi` and its four methods, `$startThreads`/`$stopThreads`
no-ops, `rngSeed` as `seed` on `bart` and `dbartsControl`, `twopart`
naming `hurdle.lognormal`, `sigma` as `sigest` on `dbarts` and
`dbartsSpec`, `bart2`); one test checks the registry against NAMESPACE,
its NEWS half deferred to S4 by an `exit_file` skip. `sigest` replaces
`sigma` on `dbarts()`/`dbartsSpec()`; [`dbartsData`](../../R/data.R)'s
`sigma` slot is unchanged. `hurdle.lognormal`/`twopart` leave
`dbarts()`'s family list, refused there by name;
[`bart2Hurdle`](../../R/bart.R) is reached from `bart` only.
[`getPointer`](../../R/dbarts.R) and `setState` refuse a 0.9-x state by
the plan's message, reached from [`predict.bart`](../../R/generics.R);
`pdbart`/`pd2bart` redirect to `bartBT`. man/bart.Rd and man/bartBT.Rd
swap and man/dbarts-deprecated.Rd is new (usage/alias/arguments correct
for check, prose rewrite left to S4); ten design-doc cites the rename
broke were repointed.

Real diff: 131 files including the header slice underneath; front-door
alone roughly R 665+/231-, tests 1070+/629- over 77 files, man
845+/808- - over budget because ~70 test files moved bart2 to bart and
two Rd files were renamed.

Gates, run independently by review and again on the merged tree with
the header slice: tinytest 8085/0 alone, 7877/0 on the stack;
equivalence 50/12/11 identical, 0 skipped, no "max |z|" line;
`R CMD check --as-cran` OK, 0 notes, clean tarball; `lint_package`,
`air format --check`, pkgdown check clean; NEWS parses;
check-doc-freshness.R, check-rc-codoc.R exit 0; `bartBT` at 0.9-34
defaults reproduces pre-slice `bart`'s gaussian/probit draws bitwise
against a throwaway install of 642f1a54; modern `bart` equals old
`bart2` bitwise; the shim equals `bartBT` bitwise, one warning per
session. Mutation probes: disabling the positional guard and the
shim's `warnOnce` each fail their tests, as do a no-op
`refuseLegacyState` and a dropped `x.train` from the shim's name list.

Review findings fixed before landing: [`bart2Hurdle`](../../R/bart.R)'s
two component calls forwarded through `bart2` and each raised its
warning (now call `bart` directly, zero warnings); the legacy door's
formula-path factor refusal was unreachable and the `dbarts()` message
pointed at `bart2` (fixed at the spec.R caller string, the data.R
suggestions, `bartBT`'s `tryCatch`); [`setState`](../../R/dbarts.R)
called any non-state object a 0.9-x fit (now gated on the 0.9-x shape);
remedy strings spelled `bart(x.train, y.train, ...)`, which the shim
would itself capture (now `bart(x, y, ...)`); one false sentence in
man/bartBT.Rd.

Remaining: S2 (family objects, consolidation, na.action), S3 (xbart),
S4 (manual, NEWS, records, consumer ports); the ~40 error strings and
~100 manual references still naming `bart2()` go with S4;
docs/design/multinomial-mutation-arc.md names man/bart2.Rd by line
number in backticked prose, not a cite - S4 decides what to do with it.
