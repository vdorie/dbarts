# small-rulings-1108-engine: two draw-moving rulings of 2026-11-08

Status: LANDED 2026-11-10 (669f62cbd..462e3b24a; dec-B714 and dec-B690 for the rulings, dec-A322 for the calls made building them)

agent: opus implementer, one; blind critique taken 2026-11-09; one opus reviewer who runs the mutants below.
rng: POSTERIOR-CHANGING at the call only: a `$setForest` call that holds a released forest with no basis under a
drawn k (probit, logistic, or gaussian given a prior on k) and hands back its sd now leaves k drawn (item 1); a
weighted binomial fit's test-row predictive draws are Binomial(k + f, p) where they were one trial (item 2). NEUTRAL,
bit for bit, for every sampler draw and every other call.
window: before 1.0-0; engine slices stay serial.
budget: ~240 lines (R ~35 changed, -15 net: item 1 -25, item 2 +10; tests ~160, help and docs ~70). No src/ change.

## Goal

The two rulings are built: a setting given back to `$setForest` is compared with the value in force before the call,
whatever else the call changes (dec-B714, extending dec-B538), and a weighted binomial fit's test row draws at its own
k + f trials, a test-frame weight of 0 holding it out (dec-B690, revising the call that drew one trial). Both are R
code, but each moves draws, so the tier is "Changes draws" ([Process by risk](README.md#process-by-risk)).

## Items

### 1. setforest-roundtrip-sd (dec-B714)

Ruling: a setting given back to `$setForest` is compared with the value in force before the call, so a handed-back sd
is no change on every model, whatever else the call changes; to pin an sd, state a different one or the same in a
second call.

Premise, run on a library whose door functions deparse identical to main d3532bd4f (a released forest with no basis,
4 sweeps, `g <- s$getForest(1); g$amplitude <- fixed()` handed back): probit, sd 1.574 becomes a statement, k 1;
logistic the same (sd 2.85, k 1); gaussian under `leaf.prior = normal(k = chi())` the same (sd 0.80, k 1);
`fixed()` alone keeps k drawn (1.906) and the reader's sd NULL; gaussian at a fixed k, the sd is no change and k is 2.
A chi prior of infinite scale, as `chi(1.25, Inf)`, gives no size, so there the release is refused and the case is
not reached.

Change: in [`changedForestSettings`](../../R/forestDoor.R) the sd is compared with `inForce$sd`, the value before
the call, like every other name; the `sdUnder` argument and its branch go, and so does
`sdUnderCoefficient`, whose one caller is [`writeForestDoor`](../../R/forestDoor.R). Only
its rejoin branch under a `dbartsChiHyperprior` differed from `inForce$sd`, so only that case moves.

Draws: a `$setForest` call naming `amplitude = fixed()` and the sd in force on a released forest with no basis that
states no sd, on probit, logistic and a gaussian fit given a prior on k (the hazard, aft, ordinal and negative
binomial families refuse a forest() term with a basis, so none reaches it): from the call on, the forest follows the leaf prior with k drawn, as after `fixed()` alone, where it was
pinned at the sd with k 1. No equivalence scenario, exact gate or snapshot passes an sd beside an amplitude change
(bcf-equivalence.R's scenario e writes an sd alone, forest-prior-exact.R and probit-k-scale-exact.R an amplitude
alone).

Tests, in [test-coefficient-prior-live.R](../../inst/tinytest/test-coefficient-prior-live.R):
- "an sd beside a hold under probit is the forest's own" is rewritten, as
  ["an sd in force beside a hold is no change"](../../inst/tinytest/test-coefficient-prior-live.R), on probit, on logistic and on gaussian with `leaf.prior = normal(k = chi())`: the sd in force beside `fixed()` leaves the record's sd `NA`, `getK` drawn
  (not 1) and the reader's sd NULL; records and the draws after equal those of a twin given `fixed()` alone by the
  same route. The ruling's own call, the whole `getForest` list with `amplitude` set to `fixed()`, gives the same.
- Then, on the same three, a second call `$setForest(1, sd = sd)` states it: k pinned at 1, records and draws equal
  to the creation twin `forest(x1, sd = sd, amplitude = fixed())` the old pin used.
- A different sd beside `fixed()` under probit is stated in the one call, k 1.
- Which pins discriminate: the rewritten block fails before and passes after; the second-call and different-sd pins
  pass before too and guard the mutants named for them.
- Kept unchanged: the gaussian held-k cases ("dropped", the size kept), "an sd in force beside another name is
  dropped", the inert round trips and test-forest-door.R's.

Help: the setForest docstring in [dbarts.R](../../R/dbarts.R) (rc-codoc) and its paragraph in
[dbartsSampler-class.Rd](../../man/dbartsSampler-class.Rd), after "so does a value equal to the one in force": each
name is compared with the value in force before the call, whatever else the call changes, so an sd handed back
beside `amplitude = fixed()` is no change; to pin a forest at its sd, state a different one or the same in a second
call. Design note: [forest-sd-unit.md](../design/forest-sd-unit.md) item 6, "the number in force before the call".
[coefficient-prior-live.md](coefficient-prior-live.md)'s amendment 4 gains "Revised by dec-B714."

Size: ~110 lines (R -25, tests ~100, help ~10).

### 2. binomial-test-trials (dec-B690)

Ruling: on a weighted binomial fit a test row draws at its own k + f trials as its training twin would; a test-frame
weight of 0 holds the row out (NA) and other weights do not change the trials; `predict(weights = )` still gives new
trials (dec-B594).

Premise, run on a library whose `binomialFitTrials`, `fitPpdWeights`, `ppdRowsOut` and `sampleFromPPD` deparse
identical to main: 60 rows of 5 trials, `weights = m`, a test frame of k = 30, f = 70 with weights 0, 1 (ten rows)
and 3: `fit$trials.test` NULL and `extract(fit, "ppd", sample = "test")` draws 0/1 (means 0.42 to 0.59, the weight-3
row 0.48, the 0 row NA); `predict(fit, test, type = "ppd")` draws at 100 trials (means 43 to 54, the 0 row NA);
with `weights = 2` given to predict, 0 to 2; without weights the test sample draws at 100.

Change:
- [`binomialFitTrials`](../../R/bart.R) loses the condition `length(sampler$data@weights.test) == 0L`, so a test
  frame carrying both response columns gives `trials.test` whatever its weights.
- [`fitPpdWeights`](../../R/generics.R), at the test rows with `trials.test`, returns
  `out = weights.test == 0` beside the trials, where the fit read test weights. [`sampleFromPPD`](../../R/generics.R)
  already draws a held-out row at its trials and then sets it NA, so every other draw and the stream are as without
  the weights. [`ppdRowsOut`](../../R/generics.R) needs no change.
- The fit-time refusal of a test row missing a count now reaches weighted fits too (dec-B641), and its remedy changes
  from "state the test rows' trials as weights, or drop those rows from 'test'" to the remedies of a test frame's
  missing weight (dec-B698): "fill them, drop those rows from 'test', or fit without them and give
  predict(fit, newdata, na.action = na.pass) the rows to leave at NA draws" (call 3).

Draws: `extract(fit, "ppd", sample = "test")` and `fitted(fit, type = "ppd", sample = "test")` on a weighted
binomial fit whose test frame carries both response columns: Binomial(k + f, p)
where they were Bernoulli(p). Unweighted binomial fits, test frames without the columns (one trial, 0 NA, dec-B638),
predict, and every other family are unchanged; 0.9-34 fit no binomial response. No equivalence scenario fits one.

Tests, in [test-binomial.R](../../inst/tinytest/test-binomial.R):
- The block "stated test weights keep dec-B596's meaning, one draw a row" is rewritten, with test weights c(0, 1, 3)
  (its training weights stay `w = 1`, still a weighted fit): `trials.test` is `k + f`; the test draws are NA at row 1
  and match, at rows 2 and 3, the same fit's draws with `weights.test` set to NULL at the same seed (so a weight of 3
  changes nothing); their mean is within 3 of 75, as the dec-B638 pin's is.
- Which pins discriminate: the `trials.test` and mean pins fail before; the match at rows 2 and 3 passes before too
  (one-trial draws on one stream) and guards the mutants named for it.
- The `gapRows` refusal pin takes the new words, and a weighted fit with `gapRows` is refused in them.
- New: a weighted fit whose test frame lacks the count columns draws 0/1, its 0 row NA.
- Kept: `predict(..., weights = 5)` drawing at most 5, and test-predict-newdata-weights.R's binomial block.

Help: [bart.Rd](../../man/bart.Rd), the family item's sentence "A row of the test sample draws at the trials its
columns give; but where the training weights name a column of test ..." becomes: a row of the test sample draws at
the trials its columns give, as predict does, a test weight of 0 holding it out and other weights not changing the
trials; the weights item's dec-B596 sentence says the binomial exception. NEWS's binomial entry: a test frame's
weights only hold rows out. Design note: [frequency-weights.md](../design/frequency-weights.md) item 10 gains
"Amended by dec-B690". The predict refusal's "or state the trials in 'weights'" stays: predict's weights are trials.

Size: ~130 lines (R ~10, tests ~60, help, NEWS and docs ~60).

## Order of work

One commit per item, each with its test file green: 1, then 2. A plain install suffices (no src/ change).

## Gates

On the slice tip, against its own library, run independently of the implementer
([RNG classes and their gates](README.md#rng-classes-and-their-gates), posterior-changing):
- The full tinytest suite in one process (`at_home = TRUE`).
- Reference build (`--enable-reference-build`, `--preclean`): equivalence.R `compare --bitwise --strict-coverage`,
  bcf-equivalence.R and multinomial-equivalence.R `compare --bitwise`, against the MANIFEST's current files; the
  four arm64 reproducibility snapshot files.
- Every gate of exact-gates.yaml in quick mode on the shipped build (none reaches these calls; run for the class).
- `R CMD check --as-cran` from a clean tarball; `lintr::lint_package()`, air, rc-codoc, win-drift, doc-freshness;
  NEWS parses.
- Not owed: tests/cpp, sanitizers and `verify-anchors` (no src/ change, no new numerics); bench-sampler compare (no
  hot path). No timing measurement and no x86 leg can run until the maintainer frees the box; the x86
  leak-sanitizer leg before landing waits for it.

Reviewer's mutants, each of which must fail a test:
- the `sdUnder` branch restored (item 1's probit and logistic pins);
- the sd compared with NULL whenever amplitude changes (the gaussian held-k pins);
- an sd beside an amplitude change dropped whatever its value (the different-sd pin);
- the `sdUnder` branch kept on gaussian only, or on the binary families only (the gaussian chi() or the probit and
  logistic cases of the rewritten block);
- `binomialFitTrials` keeping the weights condition (`trials.test` and the draws);
- `out` left off in `fitPpdWeights` (row 1 not NA);
- the trials multiplied by the test weight (the weight-3 match);
- the 0 row dropped before drawing rather than set NA after (the match at rows 2 and 3).

## Expected verdicts

Against the MANIFEST: equivalence-38da480c 55 of 55, bcf-equivalence-38da480c 19 of 19 and
multinomial-equivalence-80b1c8d4 11 of 11 "identical draws (same RNG stream)", with no |z| line, so no re-record and
no MANIFEST change; the four snapshot files pass unchanged; every exact gate passes. Only the two rewritten pin blocks
see a draw move.

## Stop conditions

Stop and report, without working around, when:
- the diff passes ~360 lines, or one item passes 1.5x its size;
- any equivalence scenario or snapshot moves;
- a test outside the two named blocks changes its expectation;
- a route to item 1's case behaves otherwise than the ruling states.

## Interactions

- This sitting's R-surface slice may touch man/bart.Rd's weights item and R/dbarts.R docstrings; conflicts are
  textual.
- No consumer reaches either path: stan4bart, bartCause, treatSens and bairrtt call no `$setForest` and fit no
  `cbind` binomial response; their suites are not owed.

## Landing records

The ledger entry for the calls made; this plan's Status and Landing note; the TODO items setforest-roundtrip-sd and
binomial-test-trials removed; the INDEX row.

## Calls made

1. Item 2 is classed posterior-changing, its predictive draws taking another distribution though no sampler draw
   moves. Alternative: shifting. The gates differ only by the exact gates and the design note, both run.
2. `sdUnderCoefficient` and the `sdUnder` argument are removed, not left unused. Alternative: keep them for a later
   rule; dead code in a door's comparison invites a wrong revival.
3. The missing-count refusal's remedy loses "state the test rows' trials as weights": a test weight never gave a
   binomial test row trials (it drew one 0/1 row, and after this slice it only holds rows out), so the remedy never
   worked. It takes the missing test weight's remedies, dec-B641 treating a missing count as a missing weight is.
   Alternatives: keep the words; "complete them, or drop those rows from 'test'". Not marked: a remedy that does not
   work is a defect, and the new words are a sibling refusal's.
4. The refusal reaches a weighted fit's test frame whose count columns have some values missing, which today fits and
   draws one trial; dec-B641 refuses a missing count at a test frame, and dec-B698 a missing test weight on any row,
   so a weight-0 row missing a count is refused too. Alternative: skip the check at weight-0 rows, whose draws are NA
   anyway. Not marked: the rulings state it.
5. Item 1's pins cover probit, logistic and gaussian given `normal(k = chi())`, each after a release under
   `student(3)` and under `normal()`. No other family reaches the case: hazard, aft, ordinal and negative binomial
   fits refuse a forest() term with a basis. Corrected in the fix round, which found the earlier "the hazard families,
   the same branch" untrue. Alternative: none left to pin.

Nothing here changes what a user sees beyond what the two rulings state, so no call is marked for the maintainer.

## Landing note

Landed 2026-11-10 as 669f62cbd..462e3b24a on main, by fast-forward: the plan (05e68c0a2) and its amendment after a blind critique (708ddc9ca), the two items (de50f1cc2 for dec-B714, 3016427df for dec-B690), the plan's cites of the removed helper (f44a52bd4) and the review's two fix commits (7f4635cf2, 462e3b24a). The calls made building them are dec-A322.

Review: one opus review that ran the plan's mutants and two of its own, LAND after one fix round, checked. The fix round pinned a release under normal() as under student(3) (an old rule kept for a normal() release had survived), limited the help's "the same in a second call" to where k is drawn, and qualified the weights item of bart.Rd for a binomial fit.

Gates on 462e3b24a: the suite in one process 43706 results, 0 failed; on the reference build equivalence 55 of 55, bcf-equivalence 19 of 19 and multinomial-equivalence 11 of 11, all identical draws on the same stream, no re-record; the four arm64 snapshot files 27 results, 0 failed; the exact gates at quick 45 of 45; CI 10 of 10. No x86 leg and no tests/cpp, no C/C++ change.

Draws moved, as planned: a handed-back sd beside `fixed()` under a drawn k (probit, logistic, gaussian given `normal(k = chi())`) leaves k drawn; a weighted binomial fit's test-row draws (extract, fitted) are Binomial(k + f, p). Nothing else moved.
