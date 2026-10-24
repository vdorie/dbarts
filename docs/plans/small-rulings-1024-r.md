# small-rulings-1024-r: five small R-surface rulings of 2026-10-23 and 2026-10-24

Status: PLANNED 2026-10-24 (written to bartcore bbd4f0fd; dec-B478, dec-B481 to dec-B484 for the rulings)

agent: sonnet implementer, one, one commit per item; then one sonnet review.
rng: NEUTRAL. No engine draw moves. Item 2 changes which calls the bridge refuses, never what an accepted call draws. Item 3 computes a row it then sets to NA or drops from the result, so a posterior predictive draw at a training row under na.pass or na.omit is made over every training row and the kept rows' draws are not those of the call at newdata with the dropped rows removed; no sampler draw moves. The four seeded-drift snapshot files are checked on the reference build.
window: pre-release.
budget: ~900 lines with tests, upper figure 1500.

## Goal

Five maintainer rulings that change the R surface (and, for item 2, one refusal in the bridge) are built, each with help that states the rule, tests that pin it and the refusals that remain, a mutant per draw-free decision that the tests kill, and a NEWS line only where 1.0-0 differs from 0.9-34 for a 0.9-34 user.

## Context

Each item is a root TODO item of the same name; the ledger entry (docs/decisions.md) quotes the maintainer and states the rule. Item 1 revises dec-A213 and extends dec-A215; item 2 revises that call of dec-A224 and reads dec-B330; item 3 revises that call of dec-A225 and extends dec-B470 and dec-B472; item 4 confirms that call of dec-A225 with the message added; item 5 revises that call of dec-A225 and extends dec-B475 to the fit's tree table. The model for the plan and the build is small-rulings-1023-r.md. They were held for response-scale-rows slice D, which landed at bbd4f0fd, since items 1 and 2 touch setData's path.

## Constraints

R-surface tier of "Process by risk": the touched test files, the snapshots, lint, air, doc-freshness and the codoc checks. Item 2 touches src/R_interface_bartcore.cpp and the comment of inst/include/dbarts/dbarts.h, so tests/cpp is run plain and under ASAN/UBSAN and the header's hash is checked; no signature in dbarts.h changes. docs/decisions.md and TODO are the orchestrator's.

## Steps

1. setdata-names-refused (dec-B478). setData and setPredictor refuse names that are not the sampler's: a frame, or a matrix that carries column names, whose names are not the sampler's is refused, naming the sampler's names that are missing; the sampler's own names in another order are matched by name as setData did; a matrix with no column names is read by position with no message, as is any data handed to a sampler that has no names. Test: test-setdata-factor-levels.R and test-sampler-setData.R restated (y ~ x.new on a y ~ x sampler), and the new pins in test-setdata-records.R or a file of their own: setData and setPredictor over frames, named and unnamed matrices, a permutation, a subset, a superset, and the sampler left as it was. Consumers checked: stan4bart, bartCause, treatSens. NEWS restates the setData sentence.
2. updatescale-no-scale-families (dec-B481). The refusal of updateScale = TRUE on a gp sampler holding saved draws (src/R_interface_bartcore.cpp refuseSavedGPDrawReanchor, holdsSavedGPDraws) skips probit, logistic and ordinal, which have no response scale, and stays on gaussian, Student-t, aft and nbinom; the flat C API's setResponse and setOffset entries follow, with no signature changed. Test: probit, logistic and ordinal for setResponse, setOffset, setWeights, setActiveRows and setData taking the call, and the refusal kept on gaussian, in test-setdata-records.R and its neighbors; tests/cpp plain and under ASAN/UBSAN.
3. predict-training-na-action (dec-B482). predict with no newdata follows na.action for a missing weight or offset exactly as newdata does, by the one resolver (R/data.R resolvePredictRows) the newdata path runs: the default refuses with the same hint, na.pass gives those rows NA draws, na.omit drops them, na.exclude and na.fail as at newdata. The type-by-type reading of weights is as it was. bartBT's help drops the stated difference. Test: test-predict-training-rows.R, every type and family, against the call at newdata = the training data.
4. training-bases-null-hint (dec-B483). The refusal of a NULL entry for a forest that carries a basis, with no newdata, names fit$bases and the idiom b <- fit$bases; b[[2]] <- z1, with the forest's index. Test: test-predict-blend.R.
5. trees-forest-labels (dec-B484). The forest column of extract(fit, "trees") on a fit is a factor whose levels are the forests' labels in forest order, forest<i> where a forest has none; the sampler's $getTrees stays numbered. Every route that returns the fit's tree table is checked, and the tests that compare the column with integers are restated. Test: test-forest-labels.R and the restated files.
