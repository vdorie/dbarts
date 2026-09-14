# Landing memo check sheet

Rubric item 16: items 3 to 12 checked against the code at the tip and the
decision register. Items 13 to 15 are the maintainer's read and are not
graded here. Every failing sentence below was fixed in the memo; the fixes
are listed at the end. Loud cases were run against the installed 1.0-0 and,
where the memo says 0.9-34 behaved differently, against a private 0.9-34
library.

3. FAIL, one sentence: "The cross-validation function is now an R driver
   over the sampler, so it reaches every family and reproduces at any
   thread count." xbart admits only auto, gaussian, probit and logistic
   (R/xbart.R's family formal; the first row of feature-matrix.md's Gaps),
   and the memo's own missing-items list says so. Everything else in the
   paragraph holds: the eight families, the variance forest under gaussian
   and aft, the multi-forest family, the two constraint kinds, the two
   non-constant leaf models, DART, sparse designs, missingness in model,
   callbacks, the active-row mask and the flat C header all exist at the
   tip; the four speed ratios are the timing run's 1.735, 2.555, 1.715 and
   1.607, and the 1.6 GB is NEWS's own figure for that fit.

4. FAIL, three sentences. (a) "Loud for bart, which refuses incomplete
   predictors." bart keeps the row and models the missingness - its
   na.action defaults to na.keepPredictors - and returns 60 fitted values
   for a 60-row design with one NA, where bartBT returns 59. Nothing
   refuses. (b) "...and binary fits use a fixed k = 2." An absent k
   resolves to the front door's own default for the response type, which
   for a binary one is chi(1.5, 2) (resolveKGrid into
   resolveNodeHyperprior); 0.9-34's own binary default was chi(1.25, Inf),
   also a hyperprior, so the clause is wrong on both releases. NEWS carries
   the same error at its k entry and contradicts itself at its xbart entry.
   (c) "rbart_vi and its methods point to stan4bart". rbart_vi is a
   tombstone that names stan4bart; plot.rbart and print.rbart are simply
   gone, with no stub. Every other loud case reproduces at the tip with the
   message the memo describes, every silent one is in the code or the
   register, and the four sister paragraphs check out against the sister
   repositories.

5. FAIL, two bullets. (a) "The posterior package is still a dependency."
   It is not. DESCRIPTION lists it under neither Imports nor Suggests,
   R/diagnostics.R computes split R-hat and bulk and tail effective sample
   size in package, draws() has replaced the as_draws methods in NAMESPACE,
   and no tinytest file or manual page calls it. All of dec-B99 has landed;
   only the TODO line saying S2 to S4 are open is stale. (b) "The survival
   exact-check script is not registered in the baseline manifest, so
   nothing in CI runs it." exact-gates.yaml runs aft-exact.R on every push
   and pull request to bartcore and main. What is missing is the MANIFEST
   row, which feature-matrix.md records as an unscheduled gap, not as
   pre-merge work. Everything else is in TODO, a design doc or the
   register with matching placement, and nothing already landed is listed
   as missing: the heteroscedastic SBC arms, the variance-forest state
   accessor, stan4bart's slice move and the hyperprior study are all
   described as done.

6. FAIL, two sentences. (a) "Twenty-three scenarios agree at the rate the
   null predicts. Three differ". classic-compare.md sets four rows aside -
   zero weights, both crossvalidation rows and the unequal-cut-point probe
   - and reports twenty-two agreeing. (b) "...at 200 replications, with
   poison runs that must fail." The poisons are opt-in through SBC_POISON,
   exist only for the heteroscedastic and BCF-latent arms, and sbc.yaml
   sets that variable on none of its matrix rows. The rest matches: the
   twenty-six scenarios, twenty seeds a side, the three defects with their
   controls, the third-of-a-posterior-standard-deviation resolution, and
   the per-family stand-ins against sbc.yaml's seven 200-replication arms,
   exact-gates.yaml's twenty-five gate scripts and the scripts in
   benchmarks/R.

7. PASS. No merge or CRAN procedure, no run ids, no commit hashes, no
   register-row counts and no tallies of who decided what appear in the
   body. The one ownership phrase, "the maintainer's own items", names two
   owed items rather than assigning steps in a procedure.

8. FAIL, one sentence, in the register rather than the memo: dec-B118's
   "The release candidate therefore waits on a study that has not started,
   and on the breadth that study is meant to have". The study ran and
   reported on 2026-09-14, recommends keeping chi(1.5, 2), and only the
   maintainer's confirmation is outstanding; dec-A07 reads the same way
   ("the study is extended"). Otherwise the register passes: every entry
   sampled states what was decided in plain words, the alternative, what a
   user notices and the ruling where there is one, and no codename, slice
   label or finding code appears inside a sentence - register ids sit at
   entry ends only.

9. PASS. Every ruled entry sampled leads with the ruling as the current
   state and gives the earlier agent position as the alternative not taken;
   no entry marked as the maintainer's says it is unruled. The posterior
   entries are the clearest case - both dec-A24 and dec-B99 state the
   removal as what the package does now.

10. FAIL: the same five package-behaviour claims named under items 3, 4
    and 5 - xbart's family reach, bart's handling of a missing predictor,
    xbart's binary k default, the posterior dependency, and rbart_vi's
    methods. Every other behaviour claim in the memo was either run against
    the installed tip or read off the code at this commit, including the
    forwarding warning, the positional bart defaults, the chain-major
    combined layout, the probit weight refusal, the three-element burn-in
    refusal, the renamed run argument, the setResponse positional warning,
    the rngSeed warning, the bartBT factor-response refusal, the dropped
    RNG-kind arguments, the thread-method no-ops, the dropped NULL
    components and fitted's third argument.

11. PASS. Each "cannot" is checked and carries its reason in the text: the
    saved sampler state cannot be restored because the stored format is
    opaque and has no conversion path; the scheduled workflows cannot run
    before the merge because GitHub binds both schedule and dispatch to the
    default branch; the multi-forest family has no calibration arm because
    its amplitude chains decorrelate too slowly; the sampled k does not
    converge at any affordable length. The sister-package claims - load
    failure, compile failure, the refused predictor assignment - were each
    reproduced or read in the sister repositories.

12. PASS. Every number in the body changes what the reader would do and
    carries its meaning: the four speed ratios with their design, the
    memory figure with the fit it describes, the moved bart defaults, the
    0.29 against about 0.7 residual scale, the 2.1 times slower
    crossvalidation, the 8 percent fused-pass loss, the 2^29 enumeration
    behind the cap of 30, and the comparison's own counts and resolution.
    No list is carried that only justifies a count.

## Changes made to the memo

1. "so it reaches every family and reproduces at any thread count" ->
   "so it reaches the sampler's own priors and settings and reproduces at
   any thread count".
2. "Silent for dbarts and the sampler" -> "Silent for bart, dbarts and the
   sampler", and "Loud for bart, which refuses incomplete predictors." ->
   "bartBT drops the row, as 0.9-34 did."
3. Dropped ", and binary fits use a fixed k = 2" from the crossvalidation
   bullet; the binary k move is already the Binary responses bullet.
4. "rbart_vi and its methods point to stan4bart" -> "rbart_vi points to
   stan4bart and its plot and print methods are gone".
5. Deleted the bullet beginning "The posterior package is still a
   dependency."
6. Deleted the bullet beginning "The NEWS speed sentence still quotes a
   July estimate", the NEWS sentence having been replaced by the measured
   statement.
7. Moved the survival exact-check item out of the pre-merge list and into
   the gaps with no decision yet, as "the survival family's exact-posterior
   check runs in continuous integration but carries no row in the baseline
   manifest".
8. "Twenty-three scenarios agree... Three differ, and each difference was
   traced to a defect" -> "Twenty-two scenarios agree... The other four
   differ in three ways, and each difference was traced to a defect".
9. Dropped ", with poison runs that must fail".
10. "both report 0.7" -> "both report about 0.7" (the measured values are
    0.72 and 0.73).
11. "a C entry to shift a constant... lands" -> "C entries to shift a
    constant... land"; the record carries two entries.

## Left for the maintainer

Two errors in inst/NEWS.Rd found while checking and not fixed here, being
outside this pass: its k entry says xbart "defaults to the fixed value
k = 2 on every family, where a binary response previously inherited the
hyperprior", which its own xbart entry and the code both contradict, and
that xbart entry's "the old fixed value of 2" misstates 0.9-34's binary
default, which was chi(1.25, Inf). Separately, the missing-predictor entry
says bart "refuses incomplete predictors outright" where the code keeps the
row, and two entries write bart where bartBT is meant.
