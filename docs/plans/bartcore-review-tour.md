# bartcore: the merge review

Current at a92f8c9b (bartcore), 2026-09-23.

The case for merging bartcore into main, the step before the 1.0-0
release, is sections 1 to 6; Appendix A is a reading order for the code.
Differences are stated against dbarts 0.9-34, the release on main.

## 1. What the merge replaces

The sampling engine. 0.9-34's engine is deleted, with the C++ headers
other packages compiled against. In its place is a new C++20 engine and one
shipped header, `inst/include/dbarts/dbarts.h`, a C interface.

The R entry points. `bart` becomes the formula-first function that 0.9-34
called `bart2`, with its defaults. The BayesTree-style function with
0.9-34's argument names and defaults moves to `bartBT`. `bart2` stays for
one release as an alias of `bart`. `rbart_vi` is removed; grouped random
effects are stan4bart's.

What the new engine adds, none of which 0.9-34 could fit: Student-t,
logistic, ordinal, multinomial, negative binomial, log-normal accelerated
failure time, discrete-time hazard and hurdle responses; a variance forest
for heteroscedastic residuals; a multi-forest family of which Bayesian
causal forests are the two-forest case; monotone and interaction
constraints; linear and Gaussian-process leaves; the DART prior; sparse
predictors; missing predictors modelled rather than dropped; and per-draw
callbacks. The BCF fitting function itself, `bcf()`, lives in bartCause;
dbarts carries the engine it runs on.

The engine implements six likelihoods: gaussian, probit, logistic, aft,
ordinal and nbinom ([`ResponseFamily`](../../src/bartcore/model.hpp)). The
other models named above - Student-t, multinomial, hazard, hurdle, BCF,
heteroscedastic - are built from those six. `docs/design/feature-matrix.md`,
the table of what each model supports, has twelve rows, one per model a user
can pick.

## 2. Breaking changes for R users

The UPGRADING block at the head of `inst/NEWS.Rd`'s 1.0-0 section is the
complete list. These are the breaks most likely to bite a 0.9-34 script,
most likely first.

1. **Every fit gives different draws.** The engine draws a different random
   stream, the tree change move now satisfies detailed balance, the initial
   forest is drawn without empty leaves, the chi hyperprior on `k` samples
   the degrees of freedom asked for, and the default tree-move mixture drops
   swap from 0.1 to 0, giving its mass to birth/death (0.6, 0, 0.4 for
   birth/death, swap, change); `bartBT` keeps BayesTree's 0.5, 0.1, 0.4. A
   seeded script does not reproduce its 0.9-34 numbers. R's own random
   numbers change as well: an unseeded fit draws its chains' seeds from R's
   generator once, when it is created, and uses its own generators after
   that; a fit given `seed` leaves R's generator untouched. Where the two
   releases fit the same model with the same priors, the posteriors agree
   (section 4).
2. **The name `bart` now means a different function.** `bart` is the
   formula-first function 0.9-34 called `bart2`, so a 0.9-34 `bart` script
   runs under new defaults: 75 trees rather than 200, four chains of 500
   kept draws after 500 burn-in rather than one chain of 1000 after 100, the
   chains merged, and the factor, missing-value and binary-prior rules
   below. `yhat.train` comes back with 2000 rows where it had 1000. The
   arguments themselves still land where they were meant; what the user is
   told depends on how the call is written. A call that names a
   BayesTree-style argument (`x.train`, `ntree`, `ndpost`, `nskip`,
   `keeptrees` and the rest) is recognized and forwarded whole to `bartBT`,
   with a once-per-session warning; `bartBT` keeps 0.9-34's 31 arguments
   and their defaults, the tree-move mixture included.
   A call that names none, such as `bart(x.train, y.train, x.test)`, reads
   the same under both versions, so nothing can forward it: it runs under
   the new defaults, and the first such call in a session prints an
   informational message naming the change and `bartBT` (a call from
   another package's code does not). A fourth positional argument is
   refused, since 0.9-34 read it as `sigest` and `bart` would read it as
   `subset`. The forwarding, the message, the refusal and the package
   startup message are removed in 1.1-0.
3. **The binary leaf-scale prior moves.** For a binary response `k` is
   sampled under `chi(1.5, 2)`, a proper prior, where 0.9-34's `bart2` used
   the improper `chi(1.25, Inf)` and its `bart` fixed `k = 2`. `chi()`'s
   defaults move from `(1.25, Inf)` to `(1.5, 2)`, and with the corrected
   degrees of freedom even an explicit `chi(1.25, Inf)` is a different prior.
   `k = chi(1.5, Inf)` restores `bart2`'s old prior. `bartBT`, and a binary
   fit with monotone constraints or several forests, keep `k = 2`. Every
   other probit posterior moves.
4. **`bart2` becomes an alias, and its arguments move.** `bart2` forwards to
   `bart` with a once-per-session warning until 1.1-0. `combineChains` now
   defaults to `TRUE`, so `yhat.train` and the other draw arrays come back
   with chains and samples merged. In the merged `sigma` and `k` vectors all
   of chain 1's draws come first, then chain 2's; 0.9-34 alternated chains
   draw by draw. Six arguments move onto objects:
   `sigdf` and `sigquant` to `family = gaussian(sigma = chisq(df, quant))`,
   `power` and `base` to `tree.prior = cgm(power, base)`, `split.probs` to
   `tree.prior = cgm(split.probs = )`, and `proposal.probs` to
   `control = dbartsControl(proposal.probs = )`. The old names still work
   until 1.1-0, with a once-per-session warning. New arguments include
   `family`, `factors`, `na.action`, `tree.prior`, `node.prior` and
   `control`. `rngSeed` is spelled `seed`; the old name still works, but a
   package that passes on only the names `dbartsControl` itself accepts, as
   the released stan4bart, WeightIt and MatchIt do, drops it and runs
   unseeded.
5. **Factor predictors enter as one column.** An unordered factor is split
   on subsets of its levels and an ordered factor is split at thresholds
   between its levels, where 0.9-34 expanded both into indicator columns.
   The design, the `varcount` names and the draws all change.
   `factors = "indicators"` restores the expansion on `bart`, `dbarts`,
   `dbartsData` and `xbart`; `bartBT` expands as before.
6. **Rows with missing predictors are kept.** 0.9-34 silently dropped them.
   The default `na.action = na.keepPredictors` drops only rows with a
   missing response and models missing predictors, so the number of rows
   and the length of fitted values grow. `na.action = na.omit` restores the
   old rule, which `bartBT` keeps. A missing value in test data, in a column
   that was complete in training, is now an error.
7. **A factor response of three or more levels picks its own family.** Under
   the default `family = "auto"`, `bart` fits an unordered one as
   multinomial and an ordered one as ordinal, and says so in one line.
   `dbarts` also fits the ordered case; otherwise `dbarts`, `xbart` and
   `bartBT` refuse, naming the `bart` family to use. On the matrix
   interface 0.9-34 fit a continuous model to the integer level codes.
8. **Some misplaced arguments are refused.** Arguments that belong to a
   neighbouring method - `sample` passed to `predict`; `newdata`, `offset`,
   `weights` or `n.threads` passed to `extract`, `fitted` or `residuals` -
   are refused by name, where 0.9-34 dropped them silently. Other unknown
   names still pass silently, except on `predict`, which warns. A fractional
   count (`n.trees = 2.5`) is refused rather than truncated. A `weights` vector of the wrong length is an error rather than
   recycled. `fitted`'s third positional argument is now `ci.level`.
9. **Weights.** A probit fit refuses a weighted likelihood, which 0.9-34 fit
   incorrectly; a 0/1 vector marks rows in or out, and integer counts belong
   on `family = "logistic"`. On a gaussian fit, rows at weight zero no
   longer count toward the residual variance's degrees of freedom; on the
   comparison design in section 4 the posterior mean of `sigma` moves from
   0.29 to 0.72.
10. **Saved objects.** A fit or sampler state saved by 0.9-x is refused by
    name at restore and at `predict`; refit. A saved `dbartsData` keeps its
    old design and fits differently from a fresh one; rebuild it.

Less common: the sampler's `setResponse(y, TRUE)` now sets `updateScale`;
its mutators refresh `$state` only when passed `updateState = TRUE`; `run`
takes no per-run thread count (`numThreads` is an error, `n.threads` is
ignored). `xbart` renames `sigma` to `sigest`, takes a two-element `n.burn`,
and no longer carries a chain across folds, so reported losses rise.
`rbart_vi` stops with an error naming stan4bart. R 4.2.0 and a C++20
compiler are required.

## 3. Breaking changes for linked packages

A package that compiled against 0.9-34 used `R_C_interface.hpp`: C++ types
(`dbarts::BARTFit` and the rest), creation from `SEXP`s, and `void`
setters. All of it is gone. Instead:

- **One C header, no R types.** `dbarts.h` compiles as plain C and includes
  no R header; a consumer includes `Rinternals.h` (and anything else the old
  header pulled in, such as `Rversion.h`) itself.
- **No creation in C.** The sampler is built in R - `dbarts()`, or
  `methods::new("dbartsSampler", control, model, data)` from a
  `dbartsSpec()` triple - and the handle is `R_ExternalPtrAddr` of the
  object's `$getPointer()`. The handle dies with the R object and changes
  whenever the object re-creates its engine from a stored state, so keep the
  object reachable and re-read the handle after any R-side restore.
  Predictor and test-data updates, weights, active rows, per-forest
  settings, state save and restore, and tree extraction are R methods on
  that object. The C surface is 22 sampler entries plus three version
  queries ([`DBARTS_C_API_LIST`](../../inst/include/dbarts/dbarts.h)).
- **Setters copy.** `dbarts_sampler_setResponse` and `_setOffset` copy into
  buffers the sampler owns. Writing through the caller's array afterwards
  has no effect; call the setter again.
- **Five entries say whether the sampler supports the call.**
  `setResponse`, `setOffset`, `setSigma`, `getLatents` and `predict` return
  `int`: 1 means the call did its work, 0 means this kind of sampler cannot
  do it at all and nothing changed - `setSigma` on a probit sampler, say. A
  bad argument still raises. The answer is fixed per sampler, so probe once
  at setup; a caller that ignores a 0 runs on, conditioned on the old
  value.
- **Changed signatures.** `dbarts_sampler_predict` takes a
  `dbarts_predictor_source` struct (`dbarts_dense_predictor_source()` builds
  the dense case) and a per-call thread count. `numTrees` and `printTrees`
  take a forest index, 0 for a single forest. `dbarts_sampler_run` fills a
  caller-owned `dbarts_results` rather than returning a `Results*`.
- **Random numbers and threads.** Each chain has its own generator, seeded
  from R's stream when the sampler is created; a run never touches R's
  stream, so a consumer needs no `GetRNGstate` bracket, and
  `dbarts_setRNGState` is gone. `setNumThreads` replaces the thread start
  and stop entries.
- **Errors.** A refusal raises an R error, which long-jumps through the
  caller's frames; a consumer holding C++ objects across a call wraps it in
  `R_UnwindProtect`. Failures inside the engine unwind as C++ exceptions
  first, so nothing dbarts owns leaks. The new per-draw callback
  ([`dbarts_draw_callback`](../../inst/include/dbarts/dbarts.h)) runs on the
  thread running its chain - a worker thread on a multithreaded run, the
  caller's thread otherwise - and must not allocate R memory. It may raise
  an R error or throw a C++ exception only on the caller's thread; on a
  worker it returns nonzero to stop.
- **Loading and versions.** The consumer's `NAMESPACE` must import from
  dbarts; `Imports:` alone compiles and then fails at load. With
  `DBARTS_USE_STUBS` defined, the first call checks the major and minor
  version. Those read 1.0 and become the contract at the release; after it,
  CI fails an ABI change without a minor bump.

The four packages we maintain that use dbarts. Their full test suites
passed on 2026-09-23 against dbarts as of 2026-09-14; the C header has not
changed since, but R code has (argument forwarding through `...`, xbart
warnings).

| consumer, branch | how it uses dbarts | state |
|---|---|---|
| stan4bart, `bartcore` | flat C API through the stubs, sampler built in R | ported, 559 of 559 tests pass; its CRAN release links the deleted headers, so 0.0-14 ships with dbarts 1.0-0 |
| treatSens, `dbarts-1.0` | flat C API through the stubs; calls `bartBT`; reads two dbarts internals, `parsePriors` and `estimateSigmaFromLinearModel`, through `asNamespace` | ported, 186 of 186; its main branch links the deleted headers; not on CRAN |
| bartCause, `dbarts-1.0` | R functions only | 790 of 790; releases from this branch |
| bairrtt, `main` | R functions only | 206 of 206 unchanged; its posteriors move with item 1 of section 2 |

`TODO`'s `release` item re-runs this against the final header.

## 4. What is checked

| workflow | runs | what a pass shows |
|---|---|---|
| `check-standard` | every push | `R CMD check` clean on macOS, Windows and Linux (R devel, release, previous release); on Windows arm64 the NEON kernels agree bitwise with scalar; once a 1.x release tag exists, an ABI change without a minor bump fails |
| `cpp-tests` | every push | the C++ component tests pass on macOS arm64, and the build fails if code that handles response families leaves one out. Then, on the reference build - a configure option that runs every sum in a fixed order - the four seeded snapshot test files and the three equivalence baselines (53 main-corpus, 15 BCF and 11 multinomial scenarios) reproduce bitwise |
| `sanitizers` | every push | the whole test suite, at least 5200 results, passes under clang and gcc address and undefined-behaviour sanitizers with no finding |
| `exact-gates` | every push | 25 scripts in quick mode: 5 detailed-balance checks (birth/death, change, swap, perturb, rule_gibbs); 16 checks against an exact posterior computed by closed form or enumeration on a small design, the logistic one also compared with the independent BART package; a check that each tree's leaf draw is exact given the rest of a full-size forest; a recovery and calibration check for aft with a variance forest; and 2 draw-for-draw reductions (hazard to binary, hurdle to its parts). Then the BCF and multinomial baselines are compared across hosts |
| `lint`, `doc-freshness`, `pkgdown` | every push | style, citations in the docs resolve, the site builds; nothing about results |
| `equivalence` | weekly, from main | the 53 main-corpus scenarios agree statistically on Linux with the shipped build |
| `sbc` | weekly, from main | simulation-based calibration over eight arms - gaussian, ordinal, nbinom, Student-t, multinomial, aft, heteroscedastic gaussian, heteroscedastic aft - and 57 functionals, Bonferroni-corrected; nbinom's two dispersion functionals are waived by name as an identifiability ridge |
| `rchk` | weekly, from main | PROTECT balance in all compiled code |
| `valgrind` | nightly, from main | the whole suite under memcheck, with no leak or invalid access |
| `revdep-smoke` | monthly, from main | stan4bart, bartCause and treatSens pass `R CMD check` from their compat branches |

Reproducibility. On one host a seed gives bitwise-identical draws at every
SIMD dispatch level and every thread count
([Reproducibility contract](../architecture.md#reproducibility-contract)).
Across hosts, the BCF and multinomial baselines must match the draws to a
relative deviation of `1e-8`, and a scenario outside that bound still passes
if a weak statistical comparison cannot tell the runs apart; the main
corpus is compared across hosts statistically only.

Agreement with 0.9-34. One script in 0.9-34's vocabulary ran under both
releases, installed side by side, on 26 scenarios spanning continuous and
probit responses, weights, offsets, factors, one to 200 trees, 5000 rows,
four chains, a Gibbs loop, a predictor swap and crossvalidation, with every
moved default pinned and 20 seeds a side. 22 agree within Monte Carlo
error; the four that differ trace to decided changes - zero-weight rows in
the variance's degrees of freedom, the change move under unequal cut
counts, and crossvalidation no longer carrying a chain across folds
([What differs, and why](classic-compare.md#what-differs-and-why)). The
comparison resolves about a third of a posterior standard deviation.

## 5. What is not checked

The five scheduled workflows have never run on schedule: GitHub runs
schedules only from the default branch, so on bartcore each ran only when
forced. `equivalence`, `sbc` and `revdep-smoke` last ran green. The forced
`rchk` and `valgrind` runs failed - rchk on protection errors in the
model-matrix code, since fixed, valgrind on test assertions - and later
hand runs are clean: rchk except for the bridge's state-restore entry, too
large for it to analyse (as it will be for CRAN's run), and valgrind over
the whole suite on x86 (`docs/plans/valgrind-xbart.md`).

Things that could be wrong and would not be caught:

- **Calibration has gaps.** No SBC arm covers BCF, which scales each forest
  by a multiplier drawn with the trees. With a gaussian response, `sigma`
  and the prognostic forest's multiplier trade off and mix too slowly for
  SBC to judge; with a probit or logistic response, the multipliers stay
  correlated past lag 200 when the prognostic multiplier is large
  (`docs/plans/bcf-latent-evidence.md`). BCF has exact checks and bitwise
  baselines instead. Monotone constraints and ordered-factor
  predictors have no SBC arm. Logistic, probit, DART, linear and GP leaves
  have SBC records made by hand in August (`docs/plans/sbc-calibration.md`)
  that no workflow re-runs. Hazard and hurdle are covered by exact checks
  and by reducing draw for draw to the binary and gaussian fits they expand
  into.
- **Mixing at scale is measured, not guaranteed.** On the He and Hahn design
  at n = 10000, a single chain's minimum pointwise effective sample size is
  2 of 2500 draws and its 95 percent interval coverage is 0.82. The shipped
  four pooled chains reach 0.96 and 0.90 on the two mean functions, because
  the chains disagree and pooling widens the interval
  ([10.4 C1, the He and Hahn factorial](../design/benchmark-surfaces.md#104-c1-the-he-and-hahn-factorial)).
  Nothing in CI scores mixing.
- **The 0.9-34 comparison is run by hand**, against a hand-installed 0.9-34,
  so a later engine change could separate a scenario unnoticed. It does not
  reach the sampler's accessors, saved state, or prediction from a saved
  sampler.
- **Warm starts at two or more forests.** Growing the initial forest from
  the root works there and is tested, but no equivalence scenario or
  calibration arm reaches it; a donor warm start is refused there
  ([`refuseMultiForestWarmStart`](../../R/bartcore.R)).
- **Two settings live outside saved state.** A per-forest weight and an
  active-row mask are re-applied by the R sampler object, but a state
  installed from a donor into a fresh sampler silently starts without them,
  and comparing the saved states will not reveal it
  ([3. What engine state does not carry, and who reinstalls it](../design/bart-as-a-component.md#3-what-engine-state-does-not-carry-and-who-reinstalls-it)).
- **Mutation records are dated.** The C++ tests caught 63 of 80 planted
  engine mutations on 2026-08-24
  (`docs/plans/review-2026-08-24/mutation-B-findings.md`), not re-run since;
  the package-level mutation battery documents three malformed-state
  refusals no test catches. `benchmarks/R/composition-matrix.R`, which
  checks the feature matrix's cells, runs in no workflow.

## 6. Decided, open, and more expensive after the merge

Open before the merge:

- **The release-candidate declaration** (`TODO`'s `rc-gate`), after the
  maintainer's read of this document.
- **31 agent-made decisions carry no ruling on whether they stand**: the
  entries in section A of `docs/decisions.md` that say "Not yet ruled on",
  superseded ones aside. For 28 of them the maintainer has recorded that the
  choice was an agent's, but not yet whether it stands.
  Those fixing user-visible surface cost a deprecation cycle to change
  after release: the single `seed` and lost generator options (dec-A04),
  ordered-factor cuts at level midpoints (dec-A09), mutators that store
  state only when told (dec-A14), `fitted`'s `ci.level` (dec-A15),
  automatic response-family detection (dec-A16), fit objects whose
  component names vary (dec-A17), four common nouns exported (dec-A67),
  documented arguments that do nothing (dec-A68), and 1-based forest
  indices in R (dec-A69).
- **One known defect.** A fit made through a wrapper's `...` stores
  `..1` in its call, so `update()` on it fails (`TODO`'s
  `forwarded-call-storage`).
- **The release items the maintainer holds**: contacting lorax's
  maintainer (its example fits a three-level factor response, which 0.9-34
  coded as 0, 1, 2 and 1.0-0 refuses); contacting WeightIt's and
  MatchIt's maintainer (both call `bart2`, which is removed in 1.1-0, and
  both drop a user's `rngSeed`, item 4); closing GitHub issue #80; and
  submitting dbarts with stan4bart 0.0-14.

Decided, and scheduled after 1.0-0: real-valued nbinom dispersion and
weighted binary responses, which share one open question about approximate
Polya-Gamma draws; formal interaction heredity, with no fixed position; an
absolute scale for the residual prior, `chisq(df, scale = )`; a revisit of
the binary `k` prior once the sampler's mixing improves; C entries that
move a constant between the forest and a host's intercepts, which
stan4bart needs, as a minor header addition; and a size threshold for the
fused residual pass, which loses up to 8 percent on small fits. The mixing
research may also give the rule-Gibbs tree move a nonzero default weight,
before or after 1.0-0.

Decided for after the merge and before 1.0-0: the engine stops calling R
for its density functions, printing and error reporting, taking them
through hooks the host installs with draws unchanged under R, and the C
interface gains an entry that creates a sampler from a plain-C
specification, with an error contract that does not assume R. Today a compiled consumer creates
its sampler through R and a host without R cannot create one at all. It is
an interface change the sister packages build against, so it lands before
the C interface becomes the 1.0 contract, and they are re-verified after
it.

Decided for 1.0-0: a scale update on a response swap is refused on BCF and
other models with two or more mean forests. A heteroscedastic model has one
mean forest plus a variance forest, and there the update recalibrates the
variance forest.

More expensive after the merge, because the release fixes them:

- **R names** lock at the CRAN submission; after it, a rename costs a
  deprecation cycle. The names deprecated now expire in 1.1-0.
- **The C interface.** Version 1.0 becomes the contract; after the
  release, entries and struct fields can only be appended.
- **Saved state.** Fits saved under 1.0-0 carry a format version, and a
  later format change must keep reading them.
- **Defaults.** Changing a default prior or the tree-move mixture after
  release, the rule-Gibbs weight included, moves users' posteriors a second
  time.

## Appendix A. Reading order

For a reviewer who opens the code after sections 1 to 6. The stops follow
what a linked package can be broken by. In the four design documents only
the sections named describe the current design - about 4,900 of their
17,100 words; the rest can be skipped.

### A.1 Orientation

Read `docs/architecture.md` (about 4,560 words) whole before any code; it
states the current design, and outranks any paraphrase.

### A.2 The C interface

Open `inst/include/dbarts/dbarts.h`, the head comment's contract list and
then the entry table
[`DBARTS_C_API_LIST`](../../inst/include/dbarts/dbarts.h); then
`src/C_interface.cpp`, where
[`dbarts_sampler_run`](../../src/C_interface.cpp) shows the error path.
Then `docs/plans/pure-c-header.md`, its Goal and Decision.

Judge: whether every non-void entry says whether it returns a value or a
capability status, and whether a discarded capability 0 is an acceptable
failure mode - it leaves the sampler unchanged and the run conditioned on
what it held before.

### A.3 The engine

Open `src/bartcore/facade.hpp`, `sampler.hpp` and `chain.hpp`:
[`SamplerBase`](../../src/bartcore/facade.hpp) and its pure virtuals,
[`SamplerFacade`](../../src/bartcore/facade.hpp),
[`createSampler`](../../src/bartcore/facade.hpp) and its siblings;
[`Sampler`](../../src/bartcore/sampler.hpp),
[`predictColumns`](../../src/bartcore/sampler.hpp) fanning out through
[`fanOutPredictSlabs`](../../src/bartcore/sampler.hpp);
[`Chain`](../../src/bartcore/chain.hpp),
[`setActiveRows`](../../src/bartcore/chain.hpp),
[`columnMaskStateFeasible`](../../src/bartcore/chain.hpp). On random
numbers and threads, prefer `docs/architecture.md`.

The leaf model - the prior on a terminal node's value and its draw - is
fixed when the code is compiled, so the per-node sums run at full speed;
`SamplerBase` gives the bridge one interface over every compiled variant.
The response family is chosen when a chain is built
([`ResponseModel`](../../src/bartcore/model.hpp)).

Judge: the `ResponseFamily` switches, which carry no `default:` arm, and
that restoring a state is semantic, not bitwise.

### A.4 Multiple forests

Read the legality table first, then the code that enforces it, then the
weight it does not save.

- `docs/design/bart-as-a-component.md`,
  [2. Which mutations are legal between sweeps](../design/bart-as-a-component.md#2-which-mutations-are-legal-between-sweeps),
  [The mutation-legality table](../design/bart-as-a-component.md#the-mutation-legality-table)
  and
  [3. What engine state does not carry, and who reinstalls it](../design/bart-as-a-component.md#3-what-engine-state-does-not-carry-and-who-reinstalls-it),
  about 1,500 words: which mutations a multi-forest sampler admits, and
  the two settings saved state does not carry.
- `docs/design/multiplier-combiner.md`, the preamble's first paragraph, then
  [The model](../design/multiplier-combiner.md#the-model),
  [The amplitude layout](../design/multiplier-combiner.md#the-amplitude-layout),
  [The reparameterization](../design/multiplier-combiner.md#the-reparameterization),
  [The amplitude conditional](../design/multiplier-combiner.md#the-amplitude-conditional),
  [bcf as the K = 2 instance](../design/multiplier-combiner.md#bcf-as-the-k--2-instance),
  [Surfaces](../design/multiplier-combiner.md#surfaces) and
  [What this family does not do](../design/multiplier-combiner.md#what-this-family-does-not-do),
  about 1,560 words: the basis-and-amplitude family and where BCF sits in
  it.
- `src/bartcore/combiner.hpp`:
  [`ForestCombiner`](../../src/bartcore/combiner.hpp),
  [`AmplitudeForestCombiner`](../../src/bartcore/combiner.hpp), which saves
  the per-forest multipliers under the key `"glue"`, and
  [`MultinomialForestCombiner`](../../src/bartcore/combiner.hpp).
- `docs/design/bcf.md`, the preamble's model equation and
  [The multiplier snap and the per-forest weight (2026-08-10)](../design/bcf.md#the-multiplier-snap-and-the-per-forest-weight-2026-08-10),
  about 340 words: why a row can carry an exact-zero weight in one forest,
  and why that weight is not saved state.

Judge: which mutations the combiner refuses, and why.

### A.5 The R bridge

Open `src/R_interface_bartcore.cpp`:
[`bartcore_create`](../../src/R_interface_bartcore.cpp),
[`bartcore_run`](../../src/R_interface_bartcore.cpp), the setters,
[`bartcore_storeState`](../../src/R_interface_bartcore.cpp),
[`bartcore_setState`](../../src/R_interface_bartcore.cpp),
[`bartcore_installForests`](../../src/R_interface_bartcore.cpp),
[`bartcore_predict`](../../src/R_interface_bartcore.cpp),
[`bartcore_predictPerForest`](../../src/R_interface_bartcore.cpp),
[`bartcore_getTrees`](../../src/R_interface_bartcore.cpp); then the shared
guards [`refusedAmplitudeFamilyReason`](../../src/R_interface_bartcore.cpp),
[`refuseMultiForestMutation`](../../src/R_interface_bartcore.cpp),
[`refuseUndefinedTestFits`](../../src/R_interface_bartcore.cpp),
[`refusePinnedSigmaChange`](../../src/R_interface_bartcore.cpp) and
[`refuseNonBinaryMask`](../../src/R_interface_bartcore.cpp).
`tests/cpp/test_facade.cpp` is the facade's conformance test, one check per
`SamplerBase` virtual driven through the base.

Judge: the comment on `refusePinnedSigmaChange`, the clearest statement in
the source of why a guard keys on the family rather than an internal flag.

### A.6 Tree moves and data

Open `docs/design/empty-leaf-veto.md`,
[Where the constant is read](../design/empty-leaf-veto.md#where-the-constant-is-read),
[Is vetoed-vs-vetoed reachable? Yes; the veto is a RANK (2026-08-18)](../design/empty-leaf-veto.md#is-vetoed-vs-vetoed-reachable-yes-the-veto-is-a-rank-2026-08-18),
[What counts as empty: the weight law (2026-08-12)](../design/empty-leaf-veto.md#what-counts-as-empty-the-weight-law-2026-08-12)
and
[Which weights the predicate sees](../design/empty-leaf-veto.md#which-weights-the-predicate-sees),
about 1,470 words: why a leaf with no members vetoes a move outright while
a leaf with members but no weight is only penalized. Then the code:
[`metropolisJumpForTree`](../../src/bartcore/moves.hpp) and
[`resolveVetoRank`](../../src/bartcore/moves.hpp);
[`Tree`](../../src/bartcore/tree.hpp),
[`Tree::leafVetoRank`](../../src/bartcore/tree.hpp),
[`columnMaskSubtreeIsValid`](../../src/bartcore/tree.hpp);
[`scanOrdinalCuts`](../../src/bartcore/scan.hpp);
[`growTreeFromRoot`](../../src/bartcore/grow.hpp);
[`ColumnStore`](../../src/bartcore/data.hpp),
[`ScopedCutGrid`](../../src/bartcore/data.hpp),
[`ColumnKind`](../../src/bartcore/data.hpp) and the derived
[`kindSplitsBySubset`](../../src/bartcore/data.hpp).

Judge: detailed balance in the change move; the veto ranking; whether a
column's type (unordered, ordered, numeric) is kept apart from the rule for
how it splits, so that only the code that builds cut grids, checks input
and reports results looks at the type; and the doubled entry layout
`scanOrdinalCuts` uses for a node holding missing values.

### A.7 The capability grid

`docs/design/feature-matrix.md`, the one deep read: what each model can and
cannot do, and a Gaps section listing every missing cell. Judge the cell
values; only the citations are machine-checked.

### A.8 Build support

`configure`, `tools/`, `src/misc/` and `src/external/` are skim-only, except
`src/misc/simd.c`'s `cpuid`, which asks for subleaf 0 explicitly so AVX2 is
never read from a stale subleaf, as 0.9-34 allowed.

### A.9 Reference, as questions arise

`docs/design/INDEX.md` and `docs/plans/INDEX.md` list every design and plan
document; `docs/decisions.md` every decision and who made it;
`docs/plans/classic-compare.md` the 0.9-34 comparison; `TODO` the open
backlog, whose `release` item is the one ordered procedure.
