# Landing memo: bartcore

## The decision and what I recommend

I am asking you for three decisions: declare dbarts 1.0-0 a release candidate,
merge bartcore, the branch this release was written on, onto main, and put
1.0-0 on CRAN beside stan4bart and bartCause. Four packages call dbarts and
all four merge in one sitting: those two, treatSens, which CRAN archived, and
bairrtt, your causal-inference package, never submitted to CRAN. Merge today
and leave the candidate for later. The merge is a fast-forward of 2024 commits
with nothing on main to reconcile, and every check that fires on the commit
you would merge is green, calibration included. Merging also starts five
schedules GitHub fires from the default branch alone, so five checks that have
never run on one begin to: the memory check (the suite under valgrind),
statistical equivalence (the comparison against 0.9-34 below), calibration
(whether a posterior fit to simulated data covers the truth at its stated
rate), the unprotected-pointer analysis (R objects left unprotected across an
allocation) and the reverse-dependency smoke test (the packages that call
dbarts installed and fitted). Declaring the candidate today would reverse your
ruling that the binary node hyperprior study comes first.

What the release gains. 0.9-34 fits a continuous or probit response with a
constant in each leaf. 1.0-0 adds eight more response families, listed in
Appendix B; a second forest for a variance that moves with the predictors, and
a causal forest whose second forest models the treatment effect; linear-ridge
and Gaussian-process leaves; monotone constraints on the constant leaf, and
limits on which predictors may interact; factors as one categorical column
instead of indicators; sparse predictors; and rows with a missing predictor,
which 1.0-0 models rather than drops. It repairs three defects: the change
move omitted its proposal-density ratio and accepted the wrong proposals, as
every released dbarts does; chi(df) has drawn twice the shape it names since
0.9-10; and since 0.9-32 the engine divided the sum of squared residuals by
the response range, not its square.

Sampler-in-a-loop use, the path you built the package for, runs faster. A
predictor swap runs 19 percent faster on x86 where it is accepted and 30 where
it is rejected, 25 and 37 on arm64. Three further figures come from one host
each: a weighted sweep runs 26 to 29 percent faster at n = 1e5, a four-chain
call of one sweep drops from 105 ms to 0.12 ms, and at that size peak memory
for a four-chain fit falls 1.6 GB. Every figure compares two builds of this
branch, not 1.0-0 against 0.9-34. Compiled callers now use one documented C
header of 25 entries with a version check, not 70 unprefixed C++ symbols.

## What could change my recommendation

**The whole comparison against 0.9-34 is nine scenarios, run once, and nothing
can widen it.** The record is statistical, not bitwise. Over 329 summaries the
largest absolute Welch z between the two engines was 3.83, and in no summary
did one engine's per-seed values all fall outside the other's range. The
harness fails a summary above 4 and warns above 3, so that reading passed with
a warning. Its z treats autocorrelated draws as independent, so it is not a
calibrated p value. The nine scenarios, in Appendix B, are continuous and
probit fits at 0.9-x defaults; none exercises a new response, a factor
predictor or a missing value. The BayesTree-descended engine 0.9-34 ships is
deleted here, so nothing wider can be recorded and no check can fail on a
drift. If nine are too thin for a rewritten sampler, no work after the merge
repairs that, and you would say what else to measure.

**The binary node hyperprior study has not started.** The binary default moved
from chi(1.25, Inf), 1.25 degrees of freedom and the old improper infinite
scale, to chi(1.5, 2.0). A probit fit samples k, the prior scale on leaf
values, by default, so its k posterior and its fitted probabilities differ
from 0.9-34's. The study behind it fit 7680 models: 48 cells, four simulated
data-generating processes crossed with three sample sizes, two predictor
counts and two base rates, at eight prior settings and 20 repetitions each.
Averaged over cells the two scales differ little: the infinite scale covers
the true probability 0.936 of the time against the scale of 2's 0.977. But in
weak-signal balanced cells at n = 100 the infinite scale's k rose to a median
of 997, at coverage 0.79. It rejects the infinite scale, prefers scale 2 over
5 narrowly, and says nothing about the degrees of freedom. You asked for three
things it lacks - a harness in the repository, the degrees of freedom varied,
and wider cases, real datasets included - and ruled it a pre-candidate
requirement, not the time-permitting item it was. No harness is here, and
nothing records what the wider study costs. It blocks the candidate and
through it the CRAN submission, not the merge.

**stan4bart has not met the mixing standard you set as the price of removing
rbart_vi.** You let rbart_vi go on the condition that stan4bart get lag-one
autocorrelation below 0.8 and an effective sample size of 100 per 1000 draws
on the group-spread parameter; if it cannot, the fallback is an R-only grouped
intercept through setOffset. Neither repository measures where it stands. The
work is stan4bart's and blocks the CRAN submission, not the merge.

**Three sites in the C bridge leak memory, and no memory-check run has
passed.** An error leaves the flat C API as an R longjmp, which skips C++
destructors, so their allocations are never freed. The sweep callback in the
run entry leaks when an error escapes the callback you pass, the chain RNG
vector when an allocation fails, and the predict fan-out and the per-chain
loop when the engine throws. No argument refusal triggers them, so no valgrind
run exercises them and a green memory check would not prove them fixed. The
check has run twice: the first lost 48 bytes in the predictor-source columns,
since fixed. In the second valgrind reported nothing lost and the job failed 7
of 8164 assertions: one convergence-diagnostic snapshot off in the fourth
decimal, and six assertions that a seeded xbart returns at 2, 3 and 4 threads
what it returns at one, which nothing explains. I fix the sites and re-run the
check before the CRAN submission, step 5 below.

**You can shorten two windows of breakage, not close them.** From your merge
until each dependent's own merge, a source install of stan4bart or treatSens
from its main fails to compile. From CRAN accepting dbarts 1.0-0 until it
accepts stan4bart 0.0-14 and bartCause, a user who updates dbarts alone breaks
both: stan4bart 0.0-13 at load and bartCause 1.0-10 at its first fit. CRAN
reviews one submission at a time, so the two cannot land together.

**lorax's CRAN checks fail once 1.0-0 is published.** lorax, a CRAN package
that reads rules out of fitted tree ensembles, suggests dbarts and fits bart
in two examples with a three-level factor response, which 1.0-0 refuses. At
the branch tip both examples error and 2 of its 724 tests fail; no lorax test
passes a factor response, so what those two are is unrecorded. Your ruling was
to let it break and ask the author to rewrite.

## Merging, and in what order

This branch replaces stan4bart's sampler back end, so the merge replaces that
package. bartCause calls dbarts from R only. treatSens compiles against the C
header and also calls two unexported dbarts functions by name, one resolving
its priors and one estimating the creation-time sigma. A rename of either
breaks it with no signature for a check to catch. You either export a public
route before the release or leave it calling internals; a decision row defers
that choice to after the release. The ids below name rows in three registers:
decisions (dec-), user-visible changes (chg-) and completion (cmp-). bairrtt
needed no branch: I renamed dbartsSampler$run's fourth argument and added
formals to dbarts(), but bairrtt calls both positionally.

The order is forced. Two ported branches include the new C header, which
dbarts main does not ship, and all four declare a dbarts floor of 1.0-0 that
no released dbarts meets. dbarts merges first, and you choose only how long
the others lag. Steps 1, 4 and 6 are mine; 2, 3 and 5 are yours.

1. Mine, before any merge: re-compute the header's checksum, reinstall
   stan4bart and treatSens from scratch, run all four suites, push the ported
   branches only. You settle treatSens's two internal calls here.
2. Yours: merge bartcore into dbarts main and push. Breakage starts.
3. Yours, that same sitting: merge each ported branch into its own main.
   Source-install breakage ends.
4. Mine: repoint dbarts's reverse-dependency job at the three mains it names,
   bairrtt not among them; unpin the consumer workflows; delete the branch
   note from the stan4bart and bairrtt readmes. The install-from-GitHub steps
   stay while CRAN serves 0.9-34.
5. Yours: submit dbarts 1.0-0 to CRAN; on acceptance submit stan4bart 0.0-14
   and bartCause the same day.
6. Mine: when CRAN serves all three, delete the install-from-GitHub steps from
   stan4bart's workflows and bairrtt's.

Nothing records what to do if you abandon the candidate after step 2. CRAN
serves 0.9-34 until step 5, so CRAN users are untouched. To back out you reset
dbarts main and each dependent's main, and whoever installed from GitHub keeps
1.0-0.

## What breaks, for whom, and whether you hear about it

### After you merge dbarts, before the dependents follow

A 0.9-x script errors outright on six calls, listed in Appendix B. Most
removed and renamed names do not error: they stay reachable for one release as
stubs on one deprecation page, each either warning once a session and
forwarding or refusing and naming what to write instead. The one-release rule
is mine; no row names who deletes the stubs later. rbart_vi refuses and names
stan4bart, whose group-spread prior differs, so results move rather than
reproduce. plot.rbart and print.rbart have no stub, and neither do the deleted
C++ headers.

Thirty-eight user-visible changes alter behaviour: 34 return a different
answer to the same call, four a message or a side effect. Appendix B lists
them, and the seven a 0.9-x script meets most often beside what each returns
now.

In this window stan4bart and treatSens fail to compile from their mains,
loudly, on deleted headers. bartCause 1.0-10 installs and then fails on every
response route, loudly too: its response fit assigns into the predictor matrix
of the data object it builds, which 1.0-0 refuses, so R reports "incorrect
number of subscripts on matrix"; its grouped route errors earlier still, on
rbart_vi. bairrtt installs and no call site it uses is affected, so the same
call returns a different posterior with nothing said. No step below touches
that; bairrtt is yours.

One break outlives this window for CRAN users. stan4bart matches control
arguments against dbartsControl's formals and drops those that do not match. I
renamed rngSeed to seed, so stan4bart 0.0-13 drops the engine seed you pass
and keeps its own: the fit looks reproducible while ignoring your seed. Its
ported branch remaps the old name, so this ends at step 3 for a GitHub
install, step 5 for a CRAN one.

The compiled handshake has one hazard: a stale binary reads the wrong numbers
rather than failing. A package can record a checksum of the header it compiled
against and refuse a dbarts whose checksum differs. You ruled that check off
in both packages that could use it, so each compares only the major and minor
version, 1.0 throughout. A binary built against an early 1.0-0 header loads
against the shipped one. Appended fields are safe, since each struct begins
with its own size. A removed field is not: the binary reads one quantity from
the slot another now occupies, with no error. Step 1 rebuilds both packages,
retiring every stale binary; the CI job that would catch a checksum moving
without a version bump stays dormant until the 1.0-0 tag exists.

### After CRAN accepts dbarts 1.0-0

stan4bart 0.0-13 fails at load, R reporting a function not provided: its init
resolves 20 of the old unprefixed symbols. bartCause 1.0-10 fails as above.
treatSens has no CRAN copy to break; a local 3.0 install resolves its dbarts
symbols inside its analysis driver, not at load, so it would fail at its first
analysis call, which no build here confirms. lorax's examples and tests error.

## Whether 1.0-0 computes what 0.9-34 computed

The nine scenarios above are the whole answer. Everything else green compares
1.0-0 against itself or an independent target.

The seven per-push workflows, all green at the tip, are the package check on
five platforms plus a Windows ARM64 vector-instruction leg, the C++ component
tests, 25 scripts checking draws against an exact posterior or a move against
detailed balance, three bitwise draw comparisons against recorded baselines,
one per sampler, over 52, 12 and 11 scenarios, two sanitizer builds, lint,
documentation freshness, and the pkgdown site.

Of the 25 exact scripts, 23 derive their target analytically or by
enumeration, one of those checking part of its work against the independent
BART package; the other two assert that an expanded hazard or hurdle design
fits identically to a probit or logistic call, dbarts against dbarts. The
three bitwise comparisons reproduce every recorded output series from a fresh
install, but their baselines come from one arm64 macOS host; draws are bitwise
only within one host architecture, so those comparisons run statistically
elsewhere in CI. A compiled test caller drives all 25 C entries. The
unprotected-pointer analysis would also test the tip and has not run there.

## What I decided without you, costliest first

I made 70 decisions that cost users something and you ruled on each: you
replaced 28, claimed 9 as yours and folded 1 into a deferral, leaving 32 of
mine standing. Separately, of the rows I filed as resting on your evidence you
marked 27 as not yours. You left 22 of those 27 with no replacement decision.
Appendix B lists both sets. These eight cost the most.

- I dropped rngKind and rngNormalKind and renamed rngSeed to seed, so no fit
  can match another package's stream and stan4bart 0.0-13 drops your engine
  seed (dec-A04).
- Ordered-factor predictors get K-1 cuts at level-code midpoints, so n.cuts no
  longer applies there, draws move for any design with one, and a saved 0.9-x
  data object has no upgrade path (dec-A09).
- Mutators store state only on updateState = TRUE, where NA deferred to the
  control slot, which run and three other methods still read, so a save after
  a mutation writes stale state (dec-A14).
- family = "auto" reads probit, ordinal or multinomial from the response's
  class where 0.9-x fit gaussian on the level codes, so the model comes from
  the data (dec-A16).
- The fit object drops NULL components, so names(fit) varies with family,
  forests and options, and every reader, the package's generics included,
  tests before it reads (dec-A17).
- Errors leave the C API as longjmps with no return codes, which is what
  leaves the three leaking sites above, and validation is partial, so a bad
  pointer crashes R (dec-A31, dec-A34).
- The bitwise gate compares only against this branch's own baselines, so with
  0.9-34's engine deleted no check can notice a drift (dec-A61).
- The multinomial defaults are mine and you ruled on none of them: a counts
  row with no trials errors rather than being dropped (dec-A65).

## What is not finished, and who finishes it

**Before you declare the release candidate.** Four items: I build and run the
hyperprior harness and you read its verdict; I bump the DESCRIPTION date; I
re-record at the candidate tip the draws the causal-forest comparison checks
against; and there I run the script that checks every documented combination
of family, leaf model and option against the table of which should fit and
which should be refused. That run is left manual.

**Before the CRAN submission.** stan4bart clears its mixing standard. I fix
the three leaking C sites and get a green memory-check run. I run the
submission battery, never run at a candidate: a package check on a fresh
tarball, the suite, and those three draw comparisons at full size. I write to
lorax's author. You close the open GitHub issue.

**At and after the merge.** Merging starts the schedules of the five workflows
above. Each also runs on a manual start and on a push touching its own file,
and all five have run here that way: calibration and the smoke test green at
the tip, the other three at earlier commits. You watch for their first
scheduled runs, or hand that to me. I repoint the reverse-dependency job and
drop the branch pins at step 4.

**Proven only in part at the tip, and due again before step 5.** Five checks,
in Appendix B with what each proved and what stays unproven. Three have not
run at the tip. Calibration and the smoke test are green there, over seven
cases and against the ported branches rather than the CRAN copies. Each is
mine to re-run except the CRAN reverse-dependency sweep, a package check over
the 27 CRAN packages in Appendix B, which needs a script nobody has written.

**Abandoned.** Four items. Three are in no backlog entry and so lost unless
you want them recorded: a non-conjugate move strategy for Gaussian-process
leaves under non-Gaussian likelihoods; the causal forest's joint rescale of
the treatment scale, derived and never built; and a set of additional C
entries planned for after the release. Only the fourth, a cross-repository
sanitizer job added and removed the same day, has a decision behind it.

**Two gaps in automation, both yours and in no step below.** bartCause runs no
CI on a push, so nothing catches a dbarts change there before a user does.
treatSens runs none either, and its 3.0-1 needs a fresh CRAN review.

## Appendix A. Claims and the register rows behind them

Prefixes: dec- is docs/decisions.md, chg- and cmp- are the changes and
completion registers beside this memo.

| claim | rows and anchors |
|---|---|
| fast-forward merge, nothing to reconcile | git rev-list --count bartcore..main = 0 |
| the branch is 2024 commits | git diff --shortstat main...bartcore: 2024 commits, 903 files, +305400/-36878; src/, R/ and inst/include alone 147 files, +57430/-28489 |
| every check that fires at the tip is green, calibration included | cmp-K06, cmp-S05, cmp-S03, cmp-S06, cmp-K04, cmp-K05, cmp-K02, cmp-S07, cmp-S08; the seven workflow files' on blocks; at 3ac244a7 runs 34664514341 (R-CMD-check), 34664514381 (cpp-tests), 34664514370 (exact-gates), 34664514351 (sanitizers), 34664514328 (lint), 34664514281 (pkgdown), 34664514353 (doc-freshness), 34664514352 (sbc, all seven arms) and 34665457872 (revdep-smoke, dispatched), every one success |
| what the release adds | chg-U52, chg-U53, chg-U54, chg-U55, chg-U79, chg-U80, chg-U89, chg-C03, chg-C23; bart's family formal in R/bart.R, eight tokens beyond gaussian and probit; DESCRIPTION Description field |
| monotonicity constrains the constant leaf | R/model.R's leaf vocabulary and its refusal of a monotone constraint with a linear or gp leaf |
| the two-part lognormal fits through bart on a matrix only | R/dbarts.R's refusal pointing at bart; R/bart.R's formula, dbartsData and dgCMatrix refusals |
| the three repaired defects, and when each entered | dec-B03, dec-B04, dec-B28, chg-U58, chg-U61, chg-U34; git log over changeRule.cpp, parameterPrior.cpp and the de-scaling site |
| the four speed and memory figures, and what each compares | dec-B117, chg-U97, chg-U94, chg-U93, chg-U92; TODO's setPredictor landing note, 19 percent on x86 accepted and 30 rejected against 25 and 37 on arm64; docs/plans/engine-performance.md's weighted table and its 105 ms / 0.12 ms latency measurement; docs/design/memory-footprint.md reference case 1, 1612.2 MB; benchmarks/baselines/MANIFEST |
| 70 registered C++ callables on main against 25 flat C entries | chg-C01, chg-C03, chg-C27, chg-C30; git show main:src/R_interface.cpp's C_callMethods table |
| max absolute z 3.83 over 329 summaries, and the harness's own thresholds | dec-A61, benchmarks/baselines/MANIFEST, benchmarks/R/equivalence.R's warn and fail cutoffs |
| the nine scenario names and their settings | the results names and meta block in benchmarks/baselines/equivalence-5430fdb.rds; Appendix B |
| the hyperprior study gates the candidate, and what it found | dec-B118, dec-B106, cmp-V09; docs/plans/archive/chi-default-research.md; man/dbartsPriors.Rd's chi entry, scale = Inf the old improper prior |
| stan4bart's mixing standard and the fallback | cmp-L05, dec-B105, cmp-L06 |
| three leaking C sites, no unwind-protect | cmp-U02, dec-A31; no R_UnwindProtect in src/C_interface.cpp |
| the two memory-check runs and the seven failing assertions | cmp-K01; the log of run 34349305160; the cross-thread expect_identical calls in inst/tinytest/test-xbart-reproducibility.R |
| the two windows of breakage, and CRAN's serial review | changes register section 6, steps 2, 3 and 5; cmp-L07 |
| CRAN is untouched until step 5, and nothing records a back-out | changes register section 6, whose six steps run forward only and name no reversal; no completion row names one either |
| a local treatSens 3.0 fails at its first analysis call, not at load | changes register section 5, treatSens row; cmp-U42, cmp-U41; treatSens master's 14 callables resolved lazily inside its analysis entry |
| lorax breaks | cmp-L08, dec-B83, chg-U70 |
| 38 behavioural changes, 34 of them answer-changing | changes register section 1, rows whose breaking column reads behavioural; Appendix B |
| 27 CRAN packages declare dbarts | cmp-S09; tools::package_dependencies("dbarts", reverse = TRUE, which = c("Depends", "Imports", "LinkingTo", "Suggests")); Appendix B |
| what each dependent calls, and bairrtt's positional calls | changes register section 5; bairrtt R/irt_causal_bart.R's run and dbarts calls |
| the order is forced; the six steps and their owners | changes register section 6 |
| the reverse-dependency job names three repositories | .github/workflows/revdep-smoke.yaml's matrix |
| treatSens's two internal calls, dbarts:::parsePriors and dbarts:::estimateSigmaFromLinearModel, and the open choice | cmp-U50; treatSens master's asNamespace calls |
| the loud refusals a 0.9-x script meets | chg-U01, chg-U24, chg-U70, chg-U44, chg-U21, chg-U75 |
| removed and renamed names stay reachable for one release, the rule mine and the deletion unowned | chg-U71, dec-B76, man/dbarts-deprecated.Rd, R/tombstones.R |
| rbart_vi refuses and names stan4bart | dec-B105, chg-U01, chg-U71, R/tombstones.R |
| plot.rbart and print.rbart have no stub | git show main:NAMESPACE against R/tombstones.R |
| the seven changes met most often, and the attach message | chg-U69, chg-U04, chg-U58, chg-U05, chg-U15, chg-U16, chg-U08, chg-U59; R/hooks.R's startup message |
| bartCause 1.0-10 errors on every response route | cmp-U45; bartCause R/responseFit.R's subassignment into the predictor matrix; dec-B35, chg-U36 |
| stan4bart drops the renamed seed, and its branch stops doing so | dec-A04, changes register section 5; stan4bart bartcore R/stan4bart_fit.R's remap |
| the checksum check is off, leaving the version pair, and its CI guard is dormant | dec-B111, cmp-U16, cmp-U17, chg-C23, chg-C32, inst/include/dbarts/dbarts.h, tools/check-api-hash.sh's tag test |
| struct appends protected by the size field, removals not | dec-B59, chg-C06 |
| stan4bart 0.0-13 fails at load on 20 old symbols | changes register section 5; stan4bart src/init.cpp |
| three bitwise draw comparisons of 52, 12 and 11 scenarios, one per sampler, and where they run | cmp-S04, cmp-S02, chg-I21, benchmarks/R/equivalence.R, benchmarks/R/bcf-equivalence.R, benchmarks/R/multinomial-equivalence.R, .github/workflows/cpp-tests.yaml, .github/workflows/exact-gates.yaml |
| 25 exact scripts, 23 of them independent of dbarts | cmp-S06, .github/workflows/exact-gates.yaml, benchmarks/R/hazard-reduction.R, benchmarks/R/hurdle-reduction.R |
| all 25 C entries driven by a compiled caller | cmp-K07, cmp-K08, inst/tinytest/capi/consumer.c |
| the five schedule-only workflows and their three triggers | cmp-L02, cmp-L03, cmp-S02, cmp-K09; the five workflow files' on blocks, each carrying schedule, workflow_dispatch and a push on its own path; run 34665457872, dispatched from this branch |
| 70 costed decisions of mine, 28 replaced, 9 yours, 1 deferred, 32 standing | decisions register section A, VD column; Appendix B |
| 27 rows marked not yours, 22 of them unsuperseded | decisions register section B, VD column; Appendix B |
| the eight costliest standing decisions, the seventh pairing two rows | dec-A04, dec-A09, dec-A14, dec-A16, dec-A17, dec-A31, dec-A34, dec-A61, dec-A65; dbartsControl's formals in R/dbarts.R against git show main:R/dbarts.R's rngKind, rngNormalKind and rngSeed; fillCutsAtLevelMidpoints in src/bartcore/data.hpp; the updateState formal on the sampler class in R/A_class.R and the four methods in R/dbarts.R that still read NA from the control slot; the family = "auto" branch in R/spec.R; R/bart.R's drop of NULL elements from the fit; the Rf_error exits in src/C_interface.cpp and the contract block in inst/include/dbarts/dbarts.h; benchmarks/baselines/MANIFEST; the "every 'counts' row must have at least one trial" refusal in R/data.R and R/A_class.R |
| the four pre-candidate items | cmp-V09, cmp-U14, cmp-U13, cmp-D07; benchmarks/R/composition-matrix.R against docs/design/feature-matrix.md |
| the submission items | cmp-U19, cmp-L05, cmp-U02, cmp-K01, cmp-L08, cmp-U18 |
| the five checks proven only in part at the tip | cmp-S01, cmp-S04, cmp-S07, cmp-S08, cmp-S09 |
| the four abandoned items | cmp-X01, cmp-X02, cmp-X03, cmp-X04, dec-A36 |
| bartCause and treatSens automation | cmp-U44, cmp-U45, cmp-U41, cmp-U42 |

## Appendix B. The lists the body counts

### The eight response families 1.0-0 adds

Student-t, logistic, ordinal, multinomial, negative binomial, accelerated
failure time, discrete-time hazard, and a two-part lognormal, that last one
through bart on a matrix only.

### The six calls a 0.9-x script now errors on

rbart_vi; a three-element n.burn in xbart; a factor response of three or more
levels through a BayesTree-spelled call; a fractional value where a count
belongs; a weighted probit fit, unless every weight is 1; and a sampler
reloaded from a 0.9-x session.

### The five checks whose evidence is incomplete at the tip

| check | what it did prove | what stays unproven |
|---|---|---|
| unprotected-pointer analysis | a local re-run at zero findings over 15534 functions, after eight in the model-matrix builder | no CI run since the fix |
| bitwise agreement with the baselines | they reproduce from a second install built from scratch, run by hand | no CI artifact for that recording, though CI runs the comparisons themselves |
| calibration | all seven cases green at the tip: gaussian, Student-t, ordinal, negative binomial, multinomial, accelerated failure time, and a discrete self-check | no verdict for the families outside the matrix, the causal forest, hazard, hurdle, heteroscedastic and monotone among them |
| reverse-dependency smoke test | the three ported branches install and fit at the tip, on a run I started by hand | it names those branches rather than the CRAN copies users will pair with 1.0-0 |
| CRAN reverse-dependency sweep | 22 of 24 packages OK against a much earlier engine | CRAN lists 27 today; the sweep kept no package list, so which three arrived since cannot be named |

### The nine scenarios in the comparison against 0.9-34

Continuous Friedman data; a probit response; weighted continuous; non-default
split probabilities; a sampled k under the chi hyperprior; four chains; a
sampler whose data is replaced between sweeps; weighted with an offset; and
quantile cut points. All at 0.9-x defaults: 20 seeds, 1000 draws after 500
discarded, 200 trees.

### The seven changes a 0.9-x script meets most often

| change | what the call returns instead |
|---|---|
| bart(x, y) with two or three positional arguments fits at 1.0-0's defaults instead of forwarding to bartBT, the legacy function that keeps 0.9-x's names and defaults | four chains of 500 draws from 75 trees where 0.9-34 gave one chain of 1000 from 200, and a posterior differing by more than the draw count |
| the default move mixture gives swap's share to birth and death, and the change move includes its proposal-density ratio | wherever some split is possible, a different tree-structure posterior from every dbarts and BayesTree fit ever run, and a single-tree fit that mixes worse until you set swap's probability back |
| the binary node hyperprior default moves | a probit fit that samples k returns a different k posterior and different fitted probabilities |
| a factor predictor becomes one categorical column split by level subset | a different predictor count and split prior, and a variable-count table with one row per factor where it had one per indicator |
| a row with a missing predictor is kept and modelled | n and the length of the fitted vector exceed 0.9-x's on the same data |
| the sigma and k draws flatten chain-major, as the matrix channels already did | a script that reshapes combined sigma draws into chains pairs draws with the wrong chain, and a per-chain summary of them mixes chains |
| the initial forest is rejection-sampled until no leaf is empty | a seeded run does not reproduce 0.9-x's numbers even where nothing else changed |

The package announces the first of these seven in the message it prints at
attach, and says nothing about the other six.

### The 38 user-visible changes that alter behaviour

The four marked "message only" change a message or a side effect; the other 34
return a different answer to the same call.

| row | change |
|---|---|
| chg-U04 | tree-move mixture: swap's mass moves to birth and death |
| chg-U05 | binary node hyperprior default chi(1.25, Inf) to chi(1.5, 2.0) |
| chg-U08 | sigma and k draws flatten chain-major |
| chg-U15 | factor predictors as single categorical columns |
| chg-U16 | rows with missing predictors kept and modelled |
| chg-U17 | bart keeps incomplete predictors where 0.9-x dropped the row |
| chg-U19 | named x with an unnamed test warns where it was silent (message only) |
| chg-U24 | xbart's n.burn loses its third element; no chains across folds |
| chg-U26 | xbart's default k follows the response type |
| chg-U27 | xbart's fold loop becomes R workers, so its draw stream differs |
| chg-U28 | startThreads and stopThreads become warning no-ops (message only) |
| chg-U30 | setResponse's second unnamed argument means updateScale |
| chg-U31 | fourteen mutators store state only on updateState = TRUE |
| chg-U34 | sum of squared residuals de-scaled by the range squared |
| chg-U39 | fitted.bart's third positional argument is ci.level |
| chg-U45 | the fit object drops NULL components |
| chg-U49 | plot.bart restores graphics parameters and errors without fits (message only) |
| chg-U50 | print.bart prints a synopsis after the call (message only) |
| chg-U55 | sparse columns, splits on subsets of factor levels, a direction at each rule for a missing value, and ordered-factor cuts at level midpoints |
| chg-U56 | categorical rules widen to 64 bits; level cap 65535 |
| chg-U58 | the change move includes its proposal-density ratio |
| chg-U59 | the initial forest is rejection-sampled to non-empty leaves |
| chg-U60 | sigma's degrees of freedom count positive-weight rows |
| chg-U61 | the chi hyperprior's shape changes and the draw is capped |
| chg-U62 | chains on standard threads; results no longer depend on count |
| chg-U63 | the model-matrix builder emits NA cells for NA factor codes |
| chg-U64 | the count of usable cores comes back as a missing value where it came back as -1 |
| chg-U68 | a BayesTree-spelled bart call forwards to bartBT |
| chg-U69 | two or three positional bart arguments fit at 1.0-0's defaults |
| chg-U74 | na.action defaults to dropping missing-response rows |
| chg-U76 | xbart seeds each unit, so a seeded call reproduces at any thread count |
| chg-U77 | xbart's k grid accepts hyperprior objects; an absent k moves |
| chg-U81 | an indicator fit stores past 100 levels sparse; leaf sums shift |
| chg-U83 | n.threads defaults to min(guessNumCores(), n.chains) and warns |
| chg-U86 | a sparse-stored column reorders leaf members, shifting last bits |
| chg-U89 | keepFits turns FALSE automatically when a callback is supplied |
| chg-U94 | weighted families get the single pass that updates residuals and leaf sums together, shifting last bits |
| chg-U97 | setPredictor partitions in place, shifting leaf sums' last bits |

### The 27 CRAN packages that declare dbarts

adrftools, bartCause, bartMan, bartXViz, bundle, butcher, CausalState,
countSTAR, EBcoBART, funcml, glossa, insight, lorax, marginaleffects, MatchIt,
mcmcsae, nlfh, orbital, riAFTBART, stan4bart, tidyAML, tidypredict,
tidytreatment, tmle, twoStageDesignTMLE, voi, WeightIt.

### The 32 standing decisions of mine that cost something

dec-A04, the RNG choice cut to two generators. dec-A09, ordered-factor cuts at
level-code midpoints. dec-A10, the initial forest rejection-sampled until no
leaf is empty. dec-A11, sigma's degrees of freedom counting positive-weight
rows. dec-A13, the sampled k capped at 1e6 with no warning. dec-A14, state
stored only when you ask for it. dec-A15, ci.level third in fitted. dec-A16,
family = "auto" reading the response's class. dec-A17, NULL components dropped
from the fit. dec-A18, two behaviour changes riding in the C bridge: NA cells
for NA factor codes, and a missing value from the core-count probe. dec-A31,
errors as longjmps with no return codes. dec-A32, the documented kinds of
non-void return. dec-A34, validation partial by intent. dec-A36, the
cross-repository sanitizer job dropped the day it landed. dec-A39, two
threading mechanisms side by side. dec-A44, the causal forest's optional ridge
on the treatment forest shipped switched off, with nothing in R able to switch
it on. dec-A51, about thirty test-only accessors compiled into the shipped
engine. dec-A52, a second copy of the R-to-engine handle code shipped to every
user so the tests can call unexported entries. dec-A53, the prior-predictive
sampler re-deriving the sigma calibration in R and building a fresh sampler
per call. dec-A54, three helpers from the R-versus-C++ work exported into the
user API. dec-A55, eight registered methods whose whole body is a refusal.
dec-A56, blocks() shipped beside the older interactions(groups=) idiom.
dec-A57, Windows ARM64 support shipped before any native probe, its
architecture string hedged three ways. dec-A60, seed-locked expected values in
four files, labelled a drift tripwire and regenerated wholesale. dec-A61,
equivalence measured only against this branch's baselines. dec-A62, gate
policy written as prose rather than code. dec-A63, the test-count floor of
5200. dec-A64, 24 superseded baselines kept in the tree, with the
      near-duplicate test file names beside them. dec-A65, user-facing choices
      settled at my discretion, the multinomial defaults among them. dec-A67,
      four constructors exported as bare top-level names while the priors are
      bundled. dec-A68, three documented arguments shipped inert. dec-A69,
      forests indexed from 1 in R and from 0 in C.

### The 22 rows you marked as not yours, with no replacement decision

dec-B11, xbart refusing an over-long n.burn by name. dec-B15, Student-t
degrees of freedom on a capped grid and integer-only negative-binomial
dispersion. dec-B16, monotone leaves drawn from the exact target, whose
quadrature costs more than the rest of a monotone fit. dec-B17, the variance
forest reusing the observation-weight channel, so no family that already uses
weights can have one. dec-B18, a half-Cauchy amplitude on the causal forest's
prognostic part. dec-B20, ordinal identification with one-at-a-time cutpoint
updates, dispatched on an ordered response. dec-B25, no 8-bit predictor
storage layer. dec-B26, the support library's two thread managers archived and
cut. dec-B28, the sum of squared residuals de-scaled by the range squared.
dec-B29, getTrees reporting a missing value and directions for a categorical
rule. dec-B32, state restore semantic rather than bitwise, fits rebuilt by
resumming trees. dec-B34, a new NA at predict on a column that had none
refused by name. dec-B38, probabilities that must sum to one snapped within
the square root of machine epsilon. dec-B39, the per-update cost of the R-side
mutation path accepted rather than an opt-out built. dec-B42, the reported log
likelihood under Student-t errors being the t density rather than the gaussian
one conditional on the augmentation. dec-B43, two sampler getters refusing a
result argument by name. dec-B50, the withdrawn claim that single-chain fits
at n above 1e5 are common. dec-B56, the ABI mechanism: a major and minor pair,
one generated entry list, and an in-header checksum. dec-B57, per-symbol
lookup kept instead of a dispatch table. dec-B59, ABI structs leading with a
size field and growing by appending. dec-B62, the per-sweep callback with a
chain index and a documented order, refused where chains would run on worker
threads. dec-B66, the pass bands for how often a multi-forest move is vetoed
ratified as judged rather than measured.