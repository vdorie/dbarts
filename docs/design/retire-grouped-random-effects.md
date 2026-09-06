# Retire grouped random intercepts

Status: ACCEPTED, 2026-09-06 (gate fired on two of three seeds; proceeding;
stan4bart tau-mixing bar and bartCause group.by route are release
prerequisites)

Remove grouped random intercepts from dbarts entirely - the engine-side
[[MOD#GroupedResponse]] decorator and its tau block, the bridge that reaches
them, `rbart_vi()` and its R Gibbs loop, `R/sliceSample.R`, the nine `.rbart`
S3 methods, and the two `dbarts_results` fields - and let stan4bart be the
home for multilevel structure. Grounded in the 2026-09-06 census of the whole
path, taken at 088098ff; every count below is that census's.

## The problem

The path costs 428 engine lines, ~247 bridge lines and ~2293 R lines
(`R/rbart.R` 1334, `R/sliceSample.R` 300, the nine methods ~656), 17
wholly-grouped tinytest files at 2297 lines with a grouped block in 36 more,
5 C++ component-test functions at 504 lines, 2 of 51 equivalence scenarios,
`man/rbart.Rd` at 216 lines with grouped content in 11 further man pages, two
wholly-grouped design docs (237 + 407 lines) and two wholly-grouped plans, and
59 of the 171 `\item` entries in `inst/NEWS.Rd`'s 1.0-0 section (6 of them
wholly grouped, 52 naming a surviving entry point too, 1 a false positive on
"grouped-GAMI"). 164 of the
2099 commits on this branch touch `R/rbart.R` or `src/bartcore/model.hpp`.

Against that:

- **One consumer.** bartCause, which already ships a `stan4bart::stan4bart`
  route in the same dispatch and whose tests fit the same grouped model both
  ways - though that route refuses `group.by` today, so the migration is real
  bartCause work rather than a reroute (see Consumers). treatSens and bairrtt
  do not use it; stan4bart uses it only as a regression baseline, in one
  `at_home()` test.
- **No decisive advantage, now measured.** The 2026-09-06 comparison (three
  simulated datasets, one chain, 200 trees, 1000 burn + 1000 kept) puts the
  two within 12-27 percent on wall time and identical to two decimals on f
  RMSE after centering; the two posterior-mean f vectors are statistically
  indistinguishable at this Monte Carlo error (correlation above 0.997, but
  the centered RMSE between them is ~0.34 against each fit's own ~0.49, so
  this is agreement at one chain's precision, not near-exactness).
  Group-effect mixing is a wash once the unidentified shared level is
  separated (centered median ESS 101/25/262 for dbarts against 142/73/185;
  the level itself mixes badly in both). stan4bart covers strictly more
  structure - intercepts and slopes, crossed and nested - on gaussian and
  binomial(probit), the two families that carry essentially all of the
  grouped traffic.
- **The tau chains do not target the same posterior.** dbarts's default is
  half-Cauchy(0, 2.5 sd(y)) - scale 12.8 on the gaussian designs; stan4bart's
  `decov()` default puts Exp(1) on tau and forms the block sd as
  `tau * re_scale * dispersion`, so its random-effect sd prior is light-tailed
  and coupled to sigma. E[tau] agrees (within 0.02 gaussian, 0.07 probit) and
  the second moment does not - on the K = 5 case sd(tau) differs by 1.98x and
  the intervals are [0.814, 3.928] against [0.809, 2.364], the signature of a
  heavy tail against a light one. Any tau ESS ratio between the two is
  therefore a comparison across different targets, and is treated below as
  evidence about the two samplers' autocorrelation, not about their
  efficiency on a common problem.
- **Unrun gates.** Zero of the 12 workflows names grouped: the four SBC
  configs in `benchmarks/R/sbc.R` are absent from `.github/workflows/sbc.yaml`'s
  matrix, there is no `grouped-exact.R` among the exact gates, and
  `benchmarks/R/grouped-mixing.R` runs in no workflow and now disagrees with
  its own header's ratios. The default half-Cauchy tau prior is
  SBC-intractable through the engine's slice sampler, so the recorded SBC
  arms are gamma-prior runs of a prior users do not get by default.
- **A NO-GO defect.** [[docs/design/forest-ranef-interweaving.md#0. The load-bearing finding: there is NO cheap interweave; the fix is a collapse]]
  measured forest-ranef confounding at 1.9-3.3x tau IACT and closed the door:
  a fix costs ~800-1300 engine lines and moves grouping off the decorator seam
  onto the constant-leaf hot path.
- **Mutual refusal with the variance forest**
  ([[spec.R#"does not support grouped random effects"]]), so the feature is
  already fenced off from the newest response decoration.
- **A seam with one instantiation.** `GroupedResponse` is the only
  `ResponseModel` that holds a `std::unique_ptr<ResponseModel>`; nine other
  families are siblings, not decorations.

## The decision

Retire it. Multilevel structure becomes stan4bart's, and dbarts keeps the
conduit that lets an outer sampler - stan4bart included - drive it.

**The gate fired on two of three seeds, and the decision stands anyway.** This
design set a threshold: stan4bart slower than `rbart_vi` by more than 5x in
tau ESS per second on the K = 20 gaussian case reopens the question. Paired by
seed on that case the ratios are 1.5x, 8.9x and 303x (tau ESS 20.4/13.2,
154.6/17.3, 636.7/2.1) - so the gate fires twice out of three, and its
magnitude is seed-sensitive to two orders of magnitude at one chain per cell
with no MCSE on the estimates. The earlier "10-50x" was a range over pooled
marginals, not over paired ratios; it is withdrawn.

What survives the seed noise, and does not depend on the target mismatch above,
is the autocorrelation itself: lag-1 of tau is 0.14-0.47 under dbarts's Gibbs
draw and 0.97-0.99 under WALNUTS, on every dataset and every seed. stan4bart's
tau chain is materially worse. That is the finding; the ratio is not.

It stands because the deficit is confined to that one chain. Wall time, the
group effects, f and the estimand-bearing posterior are equivalent within this
comparison's precision, and the gap is a defect in the receiving package's
sampler rather than a reason the capability belongs here.

**The cause is not yet diagnosed, and diagnosing it is stan4bart-side work.**
The obvious remedies are not available as stated. stan4bart's random-effect
block is *already* non-centered - it samples `z_b` and forms
`b_level = T_i * z_b_level` in `src/parametric_model.hpp`, the rstanarm decov
parameterization verbatim - so 0.98 is what non-centering produces here, not
something it would fix. And a conjugate Gibbs step for tau is not conjugate
under a Gamma prior on the sd: adopting one means changing stan4bart's prior
(to an inverse-gamma, or to the half-Cauchy scale mixture dbarts uses), which
is a modelling decision, not a sampler switch. What the fix is remains open.

**Release prerequisite 1, stan4bart-side: a tau-mixing bar.** Stated
absolutely, so it neither cites `rbart_vi` nor depends on scripts that stop
running once this lands. On the K = 20 gaussian design of the comparison
(n = 2000, Friedman f, tau = 1, one chain, 1000 kept draws after 1000 warmup),
stan4bart must reach **lag-1 autocorrelation of tau below 0.8 and tau ESS of
at least 100 per 1000 draws**. The reference the bar was set from, recorded
here so it survives the deletion: on the same design `rbart_vi` measured lag-1
0.141-0.150 and tau ESS 20.4 / 154.6 / 636.7 across three seeds, against
stan4bart's lag-1 0.967 and tau ESS 2.1 / 13.2 / 17.3. The comparison scripts
and the three simulated designs are to be committed to stan4bart's
`benchmarks/` when the item is filed there, together with those numbers.

**Release prerequisite 2, bartCause-side: a `group.by` route that does not go
through `rbart_vi`.** See Consumers - bartCause's stan4bart branch refuses
`group.by` today, so the one consumer has no working path until it changes.

Both are sister-repo work, filed from those repositories. Neither gates this
deletion; both gate the lockstep release, and dbarts ships first in that pair,
so the ordering is the thing to watch: once dbarts 1.0-0 is submitted, the
prerequisites can no longer be satisfied by delaying it.

The one capability moved rather than kept is a random intercept on an AFT
response: stan4bart is gaussian and binomial(probit) today. It is not
unconsidered on that side - stan4bart's TODO carries `aft-frailty` for
multilevel AFT survival, with `dbartsSpec(family = "aft", survival = )`
already building the specification and the open question being whether the
imputed censored response is reachable through `getLatents`. So the capability
is relocated and unbuilt, not abandoned. Accepted.

## What goes

**Engine** (428 lines).

- `src/bartcore/model.hpp`: the [[MOD#GroupedResponse]] class (227 lines with
  its doc block); the free-function block
  [[MOD#TauPriorKind, logTauPrior, logTauPosterior, sliceSampleOnce, drawTauCauchyExactIG, drawGroupEffects]]
  (129 lines - `sliceSampleOnce`'s only shipped caller is the gamma tau draw);
  the four base-class hooks
  [[MOD#ResponseModel::groupEffects, ResponseModel::numGroupEffects, ResponseModel::groupTau, ResponseModel::restoreGroupEffects]].
- `src/bartcore/chain.hpp` (~50): the five
  [[CH#SamplerOptions::groupIndices, SamplerOptions::numGroups, SamplerOptions::tauPriorKind, SamplerOptions::tauPriorScale, SamplerOptions::tauSliceSteps]]
  fields, [[CH#Results]]'s `tau` and `groupEffects` members, the decorator
  construction, the state save / shape check / restore arms, and the per-draw
  record and de-scale in [[CH#Chain::storeSample]].
- [[SAM#Sampler::numGroups]] and its two call sites (6);
  [[FAC#SamplerShape::numGroups]] (3);
  [[COM#ChainStateData::groupEffects, ChainStateData::groupTau]] (3).

**Bridge** (~247 lines).

- `src/R_interface_rbart.cpp` and `src/R_interface_rbart.hpp`, whole files -
  `rbart_getFitted`, [[src/R_interface_rbart.cpp:21@916271f3]] - plus its
  include and its [[src/R_interface.cpp#rbart_fitted]] registration. This and
  the two below are history cites deliberately: a symbol cite into a file the
  change deletes cannot resolve afterwards, and `retired:` does not rescue it
  (the guard resolves the path first).
- `src/R_interface_bartcore.cpp`: [[RIB#applyGroupAttribute]] (56);
  [[RIB#refuseGroupedScaleUpdate]] and its declaration in
  `src/R_interface_bartcore_common.hpp` (~39); the ordinal, nbinom and
  variance-forest composition refusals and the
  [[RIB#"grouped random effects"]] offender string (~23); the `"tau"` and
  `"ranef"` result channels (6); the two mutation guards and the
  [[RIB#"grouped random effects fix the data at creation"]] `setData` refusal
  (5); the `"ranef"`/`"tau"` state slots and their restore arm (~14).
- `src/C_interface.cpp`: the `bartcore.groups` note, the two
  [[C_interface.cpp#refuseGroupedScaleUpdate]] calls in
  `dbarts_sampler_setResponse` and `dbarts_sampler_setOffset`, and the
  layout asserts below (~9). Two `tau` / `groupEffects` entries also leave
  [[C_interface.cpp#DBARTS_RESULTS_FIELDS]], the list macro that drives both
  the alignment asserts and the hash fold - a compile error if missed.

**C API.** [[CAPI#dbarts_results]] loses `tau` and `groupEffects` and the
sentence in its doc block that conditions them on a grouped sampler; the
`logLikelihood` paragraph loses its grouped clause;
[[CAPI#DBARTS_RESULTS_INIT]] drops two of its twelve positional initializers.
The header's own field-discipline prose says fields "append monotonically
below the marked boundary and never reorder" - it does not authorise a
removal at all, so that paragraph is amended in the same edit to say what a
pre-1.0-0 removal does.

**R** (~2293 lines).

- `R/rbart.R`, whole: `rbart_vi`, `rbart_vi_run`, `rbart_vi_fit_bartcore`,
  `rbart_vi_fit`, `packageRbartResults` and `rbart.priors`,
  [[R/rbart.R:1-1334@916271f3]].
- `R/sliceSample.R`, whole: `sliceSample` and `rejectionSample`,
  [[R/sliceSample.R:1-300@916271f3]].
- The nine methods: [[generics.R#predict.rbart, extract.rbart, fitted.rbart, residuals.rbart, plotTree.rbart, print.rbart]],
  [[R/plot.R#plot.rbart]], [[bart.R#survivalProbabilities.rbart]], and
  [[R/diagnostics.R#summary.rbart, as_draws_array.rbart, as_draws_df.rbart]].
- Beyond the methods, `R/generics.R` also carries
  [[generics.R#rbartUnusedArgs]] and its four dispatch references, and the
  live user-facing string [[generics.R#"is the grouped (rbart_vi) fit's own predict argument"]],
  pinned by [[test-generics-errors.R#"group.by"]].
- Scattered awareness: the `"rbart"` entry in [[R/hooks.R#as_draws_array]]'s
  class list, [[spec.R#"dbarts()/rbart_vi()/bart()/xbart"]] and the
  variance-forest refusal comment, and comments or dispatch lists in
  `R/partialDependence.R`, `R/data.R`, `R/utility.R`, `R/xbart.R` and
  `R/dbarts.R`.

**NAMESPACE.** Ten lines: `export(rbart_vi)` and the nine
`S3method(*, rbart)` registrations (plot, fitted, extract,
survivalProbabilities, predict, residuals, print, summary, plotTree).

**Documentation.** `man/rbart.Rd`, whole (216 lines, six aliases). Grouped
content out of `summary.bart.Rd` (three `.rbart` aliases, `tau` in three
`vars =` usage lines), `survivalProbabilities.Rd` (the `.rbart` alias, the
`group.by` argument entry, the grouped-AFT section), `sampler.Rd` (the
`updateScale`-on-grouped refusal paragraph, the grouped response-mutation
paragraph, `ranef`/`tau` in `run()`'s value docs), `bart2.Rd`, `plotTree.Rd`
(its `.rbart` alias), `bart.Rd`, `dbarts.Rd`, `xbart.Rd`,
`dbarts-embedding.Rd`, `dbartsControl.Rd`, `dbarts-package.Rd`. The `- rbart`
reference entry in `_pkgdown.yml`. The user-facing feature bullet
[[./README.md#"Grouped random effects"]]. `docs/architecture.md`, which
describes `GroupedResponse` as engine architecture and names `rbart_vi` in the
`R/` layer map - outside `docs/design/`, so the thirty-doc sweep does not
reach it. [[docs/plans/bartcore-review-tour.md#2. Breaking changes for R users]] in
three places, including the note that `grouped-mixing.R` disagrees with its
own header.

The two `rbart_vi` paragraphs in `vignettes/gibbs_sampler_mixture_model.Rmd`
point readers at `R/rbart.R` as the package's worked embedded-Gibbs example.
The replacement is concrete rather than owed: inline the deleted loop's core -
a ~60-line random-intercept Gibbs driven through `$setOffset` - as the
vignette's own example. It is the exact pattern the retained conduit exists
for, and shipping it costs nothing.

**Tests.** The 17 wholly-grouped tinytest files (2297 lines) - the 14
`test-rbart-*.R` plus `test-grouped-swap.R`, `test-slice-sample.R` and
`test-reproducibility-rbart.R` - plus the shared fixture
`inst/common/rbartGroupData.R`; the grouped block in 36 further tinytest
files; the five [[tests/cpp/test_model.cpp#testGroupedMath, testGroupedEndToEnd, testGroupedBinary, testGroupedResponseSwap, testGroupedStateRoundTrip]]
functions and their registrations in `runModelTests`; the "grouped intercepts
delegate" sub-block of [[tests/cpp/test_sampler.cpp#testActiveRows]]; the
`numGroups` assertion in [[tests/cpp/test_shape.cpp#testConstantGaussian]];
and the `ChainStateData` equality helper's two grouped comparisons in
[[tests/cpp/common.cpp#groupEffects, groupTau]].

The ABI test is its own build-breaking site: [[inst/tinytest/capi/consumer.c#capi_run_grouped]]
is a whole grouped-run entry point setting `results.tau` / `results.groupEffects`
and naming `"tau"` / `"ranef"` back to R, and the same two fields are poisoned
in its structSize probes; it is driven from
[[test-capi.R#"capi_run_grouped"]], including a grouped-plus-variance-forest
refusal. Being a C file compiled against the header, it belongs with the ABI
edit, not with the R tests.

**Benchmarks.** Two whole files: `benchmarks/R/grouped-mixing.R` (209 lines)
and [[benchmarks/R/forest-ranef-collapse-proto.R#"collapse"]] (209 lines, the
isolation prototype behind the NO-GO doc being superseded). The `grouped` and
`grouped_aft` scenarios and [[benchmarks/R/equivalence.R#fitViaRbart]];
`probeRbart` and the `"groupedRanef"` entry in
[[benchmarks/R/composition-matrix.R#probeRbart]]; the grouped poison test in
[[benchmarks/R/mutation-battery.R#"grouped-intercept precision counts members"]]
(that file also pins the current equivalence baseline - see Gates).

`benchmarks/R/sbc.R`'s grouped surface is wider than the two functions its
plan step names: beyond the four configs,
[[benchmarks/R/sbc.R#sbcMakeGroupedFit, runSbcGrouped]] and the config
extender [[benchmarks/R/sbc.R#sbcAddGrouping]], it carries the
[[benchmarks/R/sbc.R#isGrouped]] predicate and its four dispatch branches, the
usage line in the module header, and four comment blocks. The Bonferroni
denominator is unaffected: [[benchmarks/R/sbc.R#sbcMatrixFunctionals]] is a
literal over the five matrix arms, none of them grouped, so no recorded SBC
verdict moves.

**Workflows.** No workflow names grouped. The one line that moves is the
baseline pin in [[.github/workflows/equivalence.yaml#"--strict-coverage"]].

**Tooling.** `tools/regenerate-snapshots.R` hardcodes
`"test-reproducibility-rbart.R"` in a file list; the pin-site note at the head
of `benchmarks/baselines/MANIFEST` names `test-rbart-loop-callback.R`; `.lintr`
carries two rationale comments naming `rbart_vi` (cosmetic, but they cite a
symbol that will not exist).

**Docs and backlog.** [[docs/design/grouped-random-effects.md#In-core grouped random effects]]
and [[docs/design/forest-ranef-interweaving.md#6. Recommendation: go/no-go]]
flip to RETIRED / SUPERSEDED with a one-line pointer here;
[[docs/plans/group-by-exposure.md#Goal]] and
[[docs/plans/tau-slice-review.md#VERDICT (summary; detail below)]] likewise. In
`docs/design/feature-matrix.md`: the `grouped` row in each of the five tables,
the flat-C grouped cell, footnotes [f8], [f13], [f14], [f31], [f32], [f37] and
[f44] entirely, the grouped clauses in [f1], [f3], [f27], [f30], [f40] and
[f50], and the per-family `grouped` prose block. Thirty further design docs
mention it in passing. In `TODO`: the `group-by-exposure` entry is dropped;
`sparse-extensions` loses its rbart_vi-on-sparse half; `correlated-outcomes`
restates "the rbart_vi pattern" as the `setOffset` conduit;
`negative-binomial`'s "grouped NB are companion doors", the multinomial
entry's `rbart_vi` family-token clause, and the ingestion-mode entry's
"grouped frailty" all go.

**NEWS.** 59 `\item` entries in the 1.0-0 section match on grouped text, and
they are not one class. **6** are wholly grouped and go
(`rbart_vi`'s `$fit` shape, `predict.rbart`'s removed `value` argument, the
in-core Gibbs item, the grouped-sampler mutability item, the parallel seed
restore, the custom-prior-by-name fix). **52** name a surviving entry point
too - `bart`, `bart2`, `xbart`, `dbarts`, `predict.bart`, `plot.bart` - and
are EDITED to drop the `rbart_vi` clause; deleting them would destroy the
release record for changes that survive this retirement. **1** is a false
positive on "grouped-GAMI decomposition" and is left alone. Plus one new
UPGRADING line (below).

Twenty-five files are deleted outright: the two bridge files, the two R
files, `man/rbart.Rd`, the shared fixture, the two whole benchmark scripts,
and the 17 tinytest files. Everything else on this list is an edit.

## What stays, and why

- **The per-sweep conduit.** `bartcore_setOffset`, `bartcore_setResponse`,
  `bartcore_setPredictor` and the whole R5 mutation surface are untouched.
  Grouped only ever *added a guard* to them - the three
  `refuseGroupedScaleUpdate` call sites go with the function, and the conduit
  itself is what every external Gibbs consumer, stan4bart included, drives.
- **The `ResponseModel` seam** and all eight concrete families. The decorator
  pattern is cited as load-bearing architecture past grouped
  ([[docs/design/forest-combiner.md#A standalone hierarchy, not a ResponseModel subclass]]):
  a combiner reads `response_->workingResponse()` whatever decorates it, and
  that property is what keeps a future decoration from needing a class per
  pairing. What retires is the one instantiation, not the seam.
- **`refreshLatents`'s sigma parameter**, added for grouped. The ungrouped
  families ignore it, so reverting it is an optional cleanup, not part of
  this change.
- **The `bartcore.*` control-attribute mechanism**, which also carries
  `bartcore.survival`, `bartcore.variance`, `bartcore.forests` and
  `bartcore.dispersion`. Only the `bartcore.groups` key retires.
- **The state-restore machinery**, minus two `ChainStateData` members and two
  R slots.

## ABI and saved state

`dbarts_results` loses two pointer members from the middle of the struct, so
the three fields below them shift by two pointer widths and `sizeof` shrinks
by the same: the per-field `static_assert` offsets in `src/C_interface.cpp`
re-index, and the size assert goes from 11 pointer members to 9.

The structSize contract cannot cover this in principle. It is safe only under
the header's stated rule - fields append below the boundary and never reorder,
which does not authorise a removal at all - so a stale binary would pass a
*larger* structSize and the library would fill what it believes are
`logLikelihood`, `dispersion` and `residualDf` two pointers off.

In practice both that hazard and its named mitigation are inert against the
only consumer that exists. stan4bart's `setCurrentPointers` fills `sigma`,
`train`, `test`, `varcount` and `k` and leaves `varprobs`, `logLikelihood`,
`dispersion` and `residualDf` NULL - every field at or below the removed pair
is NULL on both sides of the shift, so nothing is mis-written. Symmetrically,
stan4bart does not define `DBARTS_REQUIRE_EXACT_ABI`, so the load-time check
the header offers is not armed there either. What actually protects it is the
lockstep rebuild: the two release together, dbarts first, and stan4bart is
recompiled against the new header before it ships.

[[CAPI#DBARTS_C_API_HASH]] still moves, and should - it folds each struct's
size and each field's name paired with its offset, so a removed field re-bakes
the literal (currently `0xca7b56a64c812b8dULL`) while
`dbarts_apiSignatureToken` does not, since no entry point's signature changes.
It is the signal for the next consumer, not for this one. Neither
`DBARTS_C_API_MAJOR` nor `DBARTS_C_API_MINOR` moves: pre-1.0-0 the initial
field set is still being fixed, so no version constant is owed - and the
header prose is amended in the same edit to say that, since today it speaks
only of appends.

Saved states lose the `"ranef"` and `"tau"` per-chain slots. `stateFormatVersion`
in `src/R_interface_bartcore.cpp` is internal encoding and does **not**
increment pre-release; a state written by a grouped build is not something
1.0-0 owes compatibility with. What such a state does on restore is
unspecified today - the restore path is shape-checked, so it either ignores
two extra list elements or errors - and the implementation must determine
which and say so in one sentence, because users on this branch do hold such
states.

## Gates to re-record

The equivalence baseline is `benchmarks/baselines/equivalence-d4bca4ce.rds`,
51 scenarios. It loses `grouped` and `grouped_aft` and re-records at 49. This
is a subtraction, so the P17 oracle owed in the MANIFEST row is the neutrality
proof: the 49 survivors must reproduce `d4bca4ce` **bitwise** in a non-strict
compare against it, since nothing in the deletion touches a draw any of them
makes.

`docs/plans/README.md` states no place obligation - it covers RNG classes and
their gates only. The practice comes from the last three re-record commits,
and it is not the four places the earlier drafts of this note named. All three
touched exactly: the pin in `.github/workflows/equivalence.yaml`, the MANIFEST
row with its oracle line, the scenario count in feature-matrix [f39], and
`equivBaseline` in `benchmarks/R/mutation-battery.R`. The TODO ledger entry is
NOT one of them - only one of the three touched TODO, and there is no baseline
pin in TODO today, so step 4 must not go looking for a line to update. Four
places, with `mutation-battery.R` in place of TODO.

AFT loses its only equivalence coverage: `grouped_aft` is it ([f44]), and
`benchmarks/R/aft-exact.R` is not a MANIFEST entry. Recommend adding a
standalone `aft` scenario in the same re-record rather than shipping AFT with
exact-gate-only coverage; an addition stays bitwise-neutral for the rest.

SBC needs no matrix change - the four grouped configs were never in it. The
0.9-34 cross-implementation anchor rows (`rbart`, `rbart_sym`, and the E3
adjudication of the warmup re-anchoring divergence) are a frozen historical
record and stay exactly as written.

## Consumers

**bartCause, branch dbarts-1.0. The migration is real work, not a reroute.**
`use.ranef = TRUE` is the default on `bartc()`, and it dispatches to
`dbarts::rbart_vi` in `R/responseFit.R` (twice - the `fn` assignment and a
`redirectCall`) and `R/treatmentFit.R`. The stan4bart branch sits immediately
above each of those, but it is gated on `parametric` and its FIRST act is to
refuse `group.by` - "`group.by` must be missing or NULL if `parametric` is
supplied; for varying intercepts, add (1 | group) to parametric equation". The
two are mutually exclusive today, so there is no branch to reroute to. What
the change actually costs bartCause:

- synthesise a `parametric` formula carrying `(1 | g)` from `group.by`, and
  delete that refusal;
- keep `object$group.by` and `group.effects` populated, since six sites in
  `R/generics.R` and the summary/plot surface read them for the group-level
  estimands;
- reconcile the predict surface: the `use.ranef` branches call
  `predict(fit, x.new, group.by = , combineChains = FALSE)`, a signature no
  `stan4bartFit` method has - the stan4bart branch uses `combine_chains` and
  an `aperm` to reorder the chain margin;
- decide what `crossvalidate` does, since it refuses any `bartMethod` other
  than `"bart"`.

Two direct `rbart_vi` calls in `tests/testthat/` and the `rbart_vi` link in
`man/bartc.Rd` go with it. That work is bartCause's, on VD's side of the
release - but it is release prerequisite 2, and until it lands
`revdep-smoke.yaml`'s bartCause leg is expected RED against this branch. The
plan's Verification carries that as an accepted, time-bounded red.

**stan4bart, branch bartcore.** Its `docs/design/dbarts-capability-adoption.md`
records "Do not absorb the random-intercept-only case" as the package's
self-declared niche boundary, on the ground that dbarts has the case in-engine.
That premise is what this change removes, so the note reverses: the
random-intercept-only case becomes stan4bart's. It also has a live code
dependency, not only a documentary one: `inst/tinytest/test-02-binary.R` calls
`rbart_vi` inside an `at_home()` block as a deviance baseline, and dbarts
ships first, so that test breaks against the released dbarts before
stan4bart's own release. One enumerated edit, sister-repo side.

**treatSens and bairrtt.** No use; nothing owed. Neither grep matches
`rbart`, and stan4bart never names `groupEffects`.

## The break versus 0.9-34

`rbart_vi` and its methods are removed. The exact surface a 0.9-34 script can
lose is `rbart_vi(formula, data, ..., group.by, group.by.test, prior)` and the
returned object's `$ranef`, `$tau` and `$group.by`, reached through those
methods. Nothing on `bart()`, `bart2()` or `dbarts()` changes: grouping was
never exposed there on 0.9-34 either.

The UPGRADING `\item` to add, in `inst/NEWS.Rd`'s 1.0-0 section:

> `rbart_vi` and its methods are removed. Multilevel structure - random
> intercepts and slopes, crossed and nested grouping factors - is now
> stan4bart's: `stan4bart::stan4bart` fits the same models with an
> lmer-style formula term, and more of them.

## Considered and rejected

**Keep the R loop, drop the engine path.** It works: built-in priors resolve
to plain R functions before the in-core branch chooses, so the fallback loop
does not depend on the engine at all. Rejected because it keeps ~2293 R lines,
the whole user-facing surface and its 17 test files - most of the cost the
decision is about - and re-introduces the R loop's unweighted group mean and
its warmup re-anchoring, the adjudicated E3 divergence the in-core path fixed.

**Keep the engine path, drop the R surface.** The obvious rejection - "the
decorator would have no caller" - does not hold, and should not be the reason.
With the R surface gone the path is still reachable from the flat C API, and a
flat-C consumer is exactly what stan4bart is: `bartcore.groups` plus the
`tau` / `groupEffects` channels would hand it the very conjugate tau draw
whose lag-1 of 0.15 is the entire measured gap, without its writing one. The
real reasons to reject it are three, and none is caller absence:

- `GroupedResponse` is single-factor intercepts only. It does not serve
  stan4bart's slopes, crossed or nested cases - the structures that are the
  package's whole point - so adopting it would leave two grouping engines side
  by side rather than one.
- The channel comes with dbarts's tau prior. stan4bart's Gamma-on-sd has no
  conjugate draw, so taking the engine's draw means taking half-Cauchy too - a
  modelling change to the receiving package, decided here by accident.
- It keeps the 428 engine lines, the decorator seam's sole instantiation and
  both ABI fields alive indefinitely, for a consumer that has not asked. Most
  of the composition refusals would go with the R surface anyway, so the
  saving it offers is smaller than it looks.

**Delete now, but keep the R Gibbs loop as a vignette example.** Not rejected -
adopted, as part of the plan. It is listed here because it is the cheapest
answer to a gap the deletion otherwise opens: the embedded-Gibbs vignette
currently points at `R/rbart.R` as its worked example, and the deleted loop is
a ~60-line random-intercept Gibbs over `$setOffset`, the exact conduit this
change keeps. See the Documentation entry above.

## What would reverse this

**The gate has run, and it fired on two of the three seeds.** Paired ratios of
1.5x, 8.9x and 303x against a 5x threshold: fired, seed-sensitively, at a
precision that cannot rank the magnitude. The question was reopened and
re-argued on the measurement rather than waved through, and the re-argument is
recorded under "The decision" - the durable finding is stan4bart's tau lag-1
of 0.97-0.99 against dbarts's 0.14-0.47, the deficit is confined to that one
chain, its cause is undiagnosed, and the two tau chains are not even targeting
the same posterior. What the firing bought is release prerequisite 1, an
absolute bar on stan4bart's tau chain, which is now a live obligation rather
than a hypothetical. The bar is the thing to watch: dbarts ships first, so
once it is submitted the prerequisite cannot be satisfied by delay.

Three other things would reverse the decision itself, none of them true today:

- A **second consumer** with a random-intercept need stan4bart cannot serve.
- A named need for a **random intercept on an AFT response** - the one
  capability this loses outright, since stan4bart is gaussian and
  binomial(probit) only. The riAFTBART shape the TODO's `group-by-exposure`
  entry names is exactly that shape; the entry is decision-gated and has never
  found a consumer. Weaker than it was, since stan4bart's `aft-frailty` TODO
  item already covers it - but a real consumer arriving before that item does
  would still be dispositive.
- **The tau-mixing bar proving unreachable.** Neither named remedy is
  available as stated - the block is already non-centered, and a conjugate
  Gibbs step would require changing stan4bart's prior - so the fix is
  genuinely undiagnosed rather than merely unbuilt. If stan4bart cannot get
  tau's lag-1 under 0.8, a package whose only home for a random intercept
  cannot draw that parameter is not a replacement, and the decision would have
  to be taken again on that fact.
