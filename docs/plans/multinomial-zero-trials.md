# multinomial-zero-trials: accept a count row with no trials

Status: PLANNED 2026-09-28

agent: opus (one engine edit in the softmax coupling; the R, bridge, test and doc edits ride the same slice)
rng: neutral (every input accepted today draws bitwise as before; a zero-trial row was refused, so its draws are new)
window: before the RC (TODO multinomial-zero-trials)
budget: ~520 lines, one slice. Engine ~60 (combiner ~30 code, ~30 comment), bridge ~15 (mostly comments), R ~45,
man ~15, NEWS ~3, docs ~30, tinytest ~230 (a new file plus three rewritten pins), tests/cpp ~120.

## Goal

A multinomial count row whose trials n_i are 0 is accepted at creation and through `$setCounts`, contributes
nothing to the likelihood, keeps its leaf occupancy and its fitted probabilities, and the first such row in a
session warns once (dec-B133). A fit with added zero-trial rows draws the same posterior as the fit without
them, exactly when those rows add no cut points.

## Context

- Ruling: dec-B133 (maintainer, 2026-09-28), replacing the zero-trial half of dec-A65. Public surfaces follow
  base R (dec-B132): glm keeps a (0, 0) binomial row with a fitted value, zero prior weight and no likelihood
  contribution. NEWS covers only changes against 0.9-34 (dec-B128), and the multinomial family is new in 1.0-0.
- The softmax coupling: [`MultinomialForestCombiner`](../../src/bartcore/combiner.hpp). The mechanism a
  zero-trial row reuses: [The active-row mask](../design/active-rows-mask.md#the-active-row-mask), whose
  multinomial arm skips an inactive row's K Polya-Gamma draws and zeroes its composed precision.
- Refusals today, all removed or relaxed here:
  - R: the `dbartsData` validity method,
    [R/A_class.R:781-785](https://github.com/vdorie/dbarts/blob/1d873a779c0f9392f404090cab62d072f8ffc930/R/A_class.R#L781-L785);
    the `bart()` count-matrix route,
    [R/bart.R:1252-1257](https://github.com/vdorie/dbarts/blob/1d873a779c0f9392f404090cab62d072f8ffc930/R/bart.R#L1252-L1257);
    [`validateMultinomialCounts`](../../R/data.R), which serves both `dbartsData` and `$setCounts`,
    [R/data.R:1630-1632](https://github.com/vdorie/dbarts/blob/1d873a779c0f9392f404090cab62d072f8ffc930/R/data.R#L1630-L1632).
  - Bridge: [`createMultinomialCountsHolder`](../../src/R_interface_bartcore.cpp),
    [src/R_interface_bartcore.cpp:3645-3647](https://github.com/vdorie/dbarts/blob/1d873a779c0f9392f404090cab62d072f8ffc930/src/R_interface_bartcore.cpp#L3645-L3647),
    and [`bartcore_setCounts`](../../src/R_interface_bartcore.cpp),
    [src/R_interface_bartcore.cpp:4061-4063](https://github.com/vdorie/dbarts/blob/1d873a779c0f9392f404090cab62d072f8ffc930/src/R_interface_bartcore.cpp#L4061-L4063).
    The comment above the first,
    [src/R_interface_bartcore.cpp:3589](https://github.com/vdorie/dbarts/blob/1d873a779c0f9392f404090cab62d072f8ffc930/src/R_interface_bartcore.cpp#L3589),
    gives the reason this plan retires.
  - Engine: nothing refuses. [`MultinomialSpec`](../../src/bartcore/combiner.hpp) and the `trials_` member,
    [src/bartcore/combiner.hpp:1771](https://github.com/vdorie/dbarts/blob/1d873a779c0f9392f404090cab62d072f8ffc930/src/bartcore/combiner.hpp#L1771),
    document n_i >= 1 as an assumption.
  - C API: none. [dbarts.h](../../inst/include/dbarts/dbarts.h) has no multinomial creation or count entry
    (TODO multinomial-doors), so the header, its docs and `DBARTS_C_API_HASH` do not change.
  - Tests pinning the refusal:
    [inst/tinytest/test-multinomial-counts-mutation.R:305-310](https://github.com/vdorie/dbarts/blob/1d873a779c0f9392f404090cab62d072f8ffc930/inst/tinytest/test-multinomial-counts-mutation.R#L305-L310),
    [inst/tinytest/test-multinomial-r5-surface.R:141](https://github.com/vdorie/dbarts/blob/1d873a779c0f9392f404090cab62d072f8ffc930/inst/tinytest/test-multinomial-r5-surface.R#L141),
    [inst/tinytest/test-multinomial-surface.R:469-474](https://github.com/vdorie/dbarts/blob/1d873a779c0f9392f404090cab62d072f8ffc930/inst/tinytest/test-multinomial-surface.R#L469-L474).
    The fuzzer's count swap lifts every empty row to one trial (`OP_SET_COUNTS` in
    [test_fuzz.cpp](../../tests/cpp/test_fuzz.cpp)).
- Docs stating the refusal: the count bullet of [The surface](../design/multinomial.md#the-surface),
  [docs/design/multinomial.md:235-237](https://github.com/vdorie/dbarts/blob/1d873a779c0f9392f404090cab62d072f8ffc930/docs/design/multinomial.md#L235-L237);
  [docs/design/empty-leaf-veto.md:278](https://github.com/vdorie/dbarts/blob/1d873a779c0f9392f404090cab62d072f8ffc930/docs/design/empty-leaf-veto.md#L278);
  the `counts` item of [`dbartsData`](../../man/dbartsData.Rd), the multinomial item of
  [`dbarts`](../../man/dbarts.Rd), the count-matrix text of [`bart`](../../man/bart.Rd), the `setCounts` item of
  [`dbartsSampler`](../../man/dbartsSampler-class.Rd) and its docstring in
  [`dbartsSampler$setCounts`](../../R/dbarts.R); the 1.0-0 NEWS item on `$setCounts` ("at least 1").

## What a zero-trial row does

Category k's one-vs-rest conditional for row i is binomial(n_i, sigmoid(eta_ik)) with eta_ik = f_ik - C_ik. The
Polya-Gamma identity writes it as e^{kappa eta} times the integral of e^{-omega eta^2 / 2} against PG(n_i, 0),
kappa = y_ik - n_i/2. At n_i = 0, y_ik = 0: the likelihood factor is the constant 1, PG(0, .) is the point mass
at 0 and kappa = 0. So "inert" means, per quantity:

- latent: omega_ik = 0 exactly, drawn with no variate, so the generator stream is the compacted data set's;
- working precision 0 in every category: no leaf sufficient statistic, branch score or leaf draw of any forest
  sees the row, and the empty-leaf veto counts it absent;
- fits: the row keeps its leaf occupancy, so every sweep updates its fit and it reports K probabilities, the
  posterior of the regression surface at x_i, as glm's fitted value at a zero-weight row;
- level centering: likelihood-invariant, so nothing to change; it counts leaves by membership, as it does for
  the mask, and a leaf holding only zero-trial rows is weight-empty and vetoed;
- log-likelihood: the engine channel is NaN-flagged for this family (`logLikelihoodIsDefined` is false) and the R
  reader computes the multinomial log-pmf, which is exactly 0 at an empty row (log `dmultinom` of an empty row);
- ppd: an R-side draw of one category from the row's probabilities, as for every row; nothing trial-dependent.

Today's engine does something else. [`MultinomialForestCombiner::drawForestGlue`](../../src/bartcore/combiner.hpp)
takes one PG(1, psi) variate before its trials loop,
[src/bartcore/combiner.hpp:1465-1467](https://github.com/vdorie/dbarts/blob/1d873a779c0f9392f404090cab62d072f8ffc930/src/bartcore/combiner.hpp#L1465-L1467),
so an empty row gets omega ~ PG(1, psi) and working response C_ik at precision omega. Nothing divides by zero
and nothing goes NaN, but integrating omega out leaves a factor cosh(eta/2)^-1, proportional to sqrt(p(1 - p)):
half a success and half a failure per category, pulling every one-vs-rest log-odds toward 0. The R comment
that PG(0, .) "would break" the working response is right about the arithmetic (0/0) and the reason the fix
must not write 0 into omega; it is not what the code would reach.

Evidence, a private build of 1d873a77 with only the R and bridge refusals relaxed (library A), then with the
engine edit of step 1 (library B); n = 200 data rows, K = 3, 50 trees, one chain, seed 7:

- A, 800 zero-trial rows added: mean |posterior-mean p - truth| at the data rows 0.199, against 0.046 with the
  same rows masked through `$setActiveRows`. Entropy 0.857 against 0.763 (truth 0.732). Not inert.
- B, 40 zero-trial rows: identical at the data rows, and within 2e-14 at the empty ones, to a sampler whose
  empty rows carry other counts and are masked. `$setCounts` zeroing 20 rows mid-run matches `$setActiveRows`
  on the same rows the same way. (Empty rows against the same sampler with those rows also masked are bitwise
  too, but that is nearly true by construction and proves little.)
- B, trees drawn by `$sampleTreesFromPrior()` (the only caller of the veto vector): empty rows against the same
  rows masked with other counts differ by 0 at the data rows; a copy of the fix whose veto still reads the
  caller's mask differs by 0.59 there, and passes every other arm (independent critique).
- B against the fit without the rows: 800 empty rows at fresh x drawn within the data range, which leaves the
  default uniform cut grid unchanged, agree with the fit without them to 5e-14 at every data row; without the
  engine edit they differ by up to 0.99 (critique). Fresh x outside the range moves the draws (0.70 in the first
  spike) because the grid is built from every row's x, as under zero gaussian weights and the mask: a change of
  tree prior, not of likelihood. Duplicated x also agrees (6e-14), but a twin always shares its data row's leaf,
  so no leaf ever holds only empty rows and the arm cannot see the veto.
- Re-creation (`$copy()`, save and load) is not bitwise even without empty rows (8e-15, critique), so a
  re-created sampler is compared with a re-created one.
- B, no zero-trial row: draws bitwise library A's. All-zero counts run at creation and through `$setCounts`,
  every reported row a finite simplex. The C++ component suite passes on B; the multinomial tinytest files
  fail only the three refusal pins.

## Agent-made calls

The ruling settles accept, inert, and warn once at creation and through `$setCounts`. These are not settled by
it; each is a recommendation VD can overturn before implementation. Calls 2, 6 and 7 were put to the maintainer and are ruled:
put to VD before the slice starts.

1. Engine form: a zero-trial row is composed into the coupling's existing global mask, `effective = mask AND
   (n_i > 0)`, rather than a per-row trials test in the sweep loops. The hot loops keep their shape, a data set
   without empty rows serves the caller's mask (or none) exactly as today, and the row becomes, in every draw,
   the inactive row the mask tests already pin.
2. Ruled 2026-09-28 (maintainer: "Sure.", after asking what all-zero weights do - accepted at creation and mid-run, as an all-off row mask is). All-zero counts are accepted, at creation and through `$setCounts`, under the same
   warning. The all-zeros mask is accepted and runs for the same reason (a stratum that empties needs no special
   case). Base R: lm accepts all-zero weights; glm on all (0, 0) rows warns "no observations informative" and
   then errors; nnet's multinom refuses any empty row. Alternative: refuse at creation only.
3. Two warning sites, each after the point where the call can no longer be refused, so a refused call never
   uses up the once-per-session key: sampler creation, right after the creation `.Call` returns in
   [`dbartsSampler$initialize`](../../R/dbarts.R), which `dbarts`, `bart` and `dbartsSpec`'s fits reach; and
   `$setCounts`, right after the `.Call` in [`bartcoreSamplerSetCounts`](../../R/bartcore.R). A refusal after
   creation (a category offset `bart` installs next, say) still leaves the warning spent, correctly: the rows
   were accepted. [`validateMultinomialCounts`](../../R/data.R) itself stays silent, so `dbartsData(counts = )`
   alone does not warn; the bridge stays silent (a direct `.Call` is internal); and re-creating a sampler from
   its data object (`getPointer` after load, `$copy`, `$setState`) does not warn, since those paths do not run
   `initialize`'s creation call and no user action introduced the rows.
4. One warnOnce key, `multinomialZeroTrials`, shared by creation and `$setCounts`: one warning per session in
   all, reading "the first such row in a session" literally.
5. Class and wording (error-style R15, R2, R3, R6): `dbartsZeroTrialsWarning` under `dbartsWarning`, message
   "multinomial count rows with zero trials (%d of %d) contribute nothing to the likelihood and still receive
   fitted probabilities". No argument is named, since the matrix arrives as `y.train`, `data` or `counts`.
6. Ruled 2026-09-28 (maintainer: "Follow glm."): `residuals` follows glm at a zero-weight row - 0 for the default and pearson types, -p (observed 0 minus the fitted probabilities) for response residuals. Superseded proposal: `residuals` returns an NA row: there is no observed proportion. Today's arithmetic
   would give NaN. Base R: glm at a zero-weight row gives a response residual of -mu (it sets y = 0 at n = 0)
   and deviance and pearson residuals of 0, deviance being the default; `rstandard` drops the row and `nobs`
   excludes it. Alternatives: 0, following glm's default residual; -p, following its response residual.
7. Ruled 2026-09-28 (maintainer: "Yep, 0."). `extract(type = "loglik")` reports 0 at the row in every draw: the multinomial log-pmf,
   and what glm's summed logLik gets from a zero-weight row, so a sum over rows is the fit's log-likelihood.
   Multinomial fits never report NaN for a masked row's log-likelihood (the mask's NaN flag is the engine
   channel of other families), so 0 departs from no multinomial convention. Cost: loo reports Pareto k = Inf
   for a constant-zero column (checked), so the manual says to drop those columns; NaN would make loo refuse
   outright.
8. `fitted`, `fitted(type = "class")`, `predict`, the ppd, `summary`'s pooled probability, `plot`'s trace panel
   and variable counts are unchanged: the row has fitted probabilities and no response-dependent reader.
   `plot`'s second panel drops empty rows; today's single-trial branch would plot p(category 1) as the
   "observed" category of an empty row
   ([R/plot.R:207](https://github.com/vdorie/dbarts/blob/1d873a779c0f9392f404090cab62d072f8ffc930/R/plot.R#L207)).
9. No statistical comparison and no exact-gate arm: the in-range fresh-x test gives the fit without the rows
   to rounding, a stronger check than a quadrature arm. No design note (neutral class); the facts land in
   multinomial.md and active-rows-mask.md.
10. NEWS: the 1.0-0 `$setCounts` item's "at least 1" becomes "at least 0, an empty row entering no likelihood";
    no new item.
11. The fuzzer's count swap may emit empty rows, reaching the new composition under the fuzz invariants.
12. xbart and `survivalProbabilities` need nothing: xbart refuses count-carrying data before any loss, and the
    other is refused for this family.

## Constraints

- RNG neutral: no draw of any input accepted today moves. A mismatch in the gates below is a defect, never a
  re-record.
- Engine stays R-agnostic; no new virtual, no facade change, no `dbarts.h` change.
- Out of scope: the three multinomial doors (TODO multinomial-doors), weights, per-forest masking.

## Steps

1. Engine, [`MultinomialForestCombiner`](../../src/bartcore/combiner.hpp): add `effectiveRows_` and a private
   `composeEffectiveRows()` (copy of `activeRows_`; if any trials are 0, materialize ones where empty and write 0
   at those rows). Call it at the end of the constructor, `setCounts` and `setActiveRows`. `drawForestGlue`,
   `formForestResponse` and `formForestVetoWeights` read `effectiveRows_` where they read `activeRows_`. Never
   write 0 into `omega_`. Update the comments: the class header, `drawForestGlue` (n_i = 0 is PG(0, .) = 0, so the
   skip is the exact law), `setCounts` (no longer a pointer swap: an O(n) recompose), `setActiveRows` (clearing
   the mask leaves empty rows out), the `trials_` member and `MultinomialSpec` (n_i >= 0), and
   [`Chain::setCounts`](../../src/bartcore/chain.hpp).
2. Bridge: drop the two trials loops in [`createMultinomialCountsHolder`](../../src/R_interface_bartcore.cpp) and
   [`bartcore_setCounts`](../../src/R_interface_bartcore.cpp) (a negative cell is still refused per cell) and
   rewrite the n_i >= 1 comments.
3. R: remove the refusals in the `dbartsData` validity method, in `bart()` and in
   [`validateMultinomialCounts`](../../R/data.R). Add one helper that, given a count matrix with an empty row,
   calls `warnOnce("multinomialZeroTrials", warningCondition(<call 5>, class = c("dbartsZeroTrialsWarning",
   "dbartsWarning")))`, and call it at the two sites of call 3, after the `.Call` returns. [`residuals.bartMultinomial`](../../R/generics.R):
   -p at rows with zero trials (call 6). [`plot.bartMultinomial`](../../R/plot.R): the second panel over rows with
   trials only, in both branches. Rewrite the `$setCounts` docstring.
4. Tests, tinytest: a new `inst/tinytest/test-multinomial-zero-trials.R`, resetting the key first
   (`env <- dbarts:::onceWarnState; env[["multinomialZeroTrials"]] <- NULL`) and counting warnings with
   `withCallingHandlers`:
   - creation through `dbarts` and `bart` warns exactly once in the file, with the class; `dbartsData` alone
     does not; a second creation is silent; a creation refused after validation (an unsupported argument) does
     not spend the key; `$setCounts` after a creation warning is silent, and after a key reset warns once;
   - bitwise, leading: empty rows against a sampler whose empty rows carry other counts and are masked
     (identical at data rows, within 1e-12 at the empty ones); then against the same sampler with those rows
     also masked (`expect_identical`, all rows);
   - the veto: the same comparison after `$sampleTreesFromPrior()`, identical at the data rows;
   - against the fit without the rows: empty rows appended at fresh x drawn strictly within the data range (the
     default uniform cut grid unchanged), enough of them that leaves holding only empty rows arise, agree with
     the fit without them to 1e-12 at every data row;
   - mid-run: `$setCounts` emptying rows equals `$setActiveRows` on those rows over the next run; restoring the
     counts equals clearing the mask, to 1e-10 (drop to the phase that is bitwise if the restored phase is not);
     a mask installed over empty rows and then cleared leaves them out;
   - readers on a `bart()` fit with empty rows: `fitted` a finite simplex there, `residuals` -p rows,
     `extract(type = "loglik")` 0 there and finite elsewhere, ppd codes in 1..K, `summary` and `plot` run
     without warning;
   - all-zero counts: warns, runs, finite simplex, at creation and through `$setCounts`;
   - `$copy()` and a save/load re-creation keep the rows inert and do not warn: a re-created sampler with empty
     rows is bitwise a re-created sampler whose rows carry other counts and are masked, at the data rows.
     Re-created against original is not bitwise even without empty rows (8e-15), so never compare those.

   Rewrite the three refusal pins as acceptance checks.
5. Tests, tests/cpp: `testZeroTrialsMultinomialKernel` beside
   [`testActiveRowsMultinomialKernel`](../../tests/cpp/test_sampler.cpp), the same structure: empty rows against
   the compacted kernel bitwise, zero precision, finite response, the next uniform equal; composed with a mask
   over other rows; a `setCounts` that fills an empty row re-admits it. In `OP_SET_COUNTS`, drop the lift of an
   empty row to one trial.
6. Docs: [The surface](../design/multinomial.md#the-surface) count bullet; a sentence in
   [Per family](../design/active-rows-mask.md#per-family) (a zero-trial row is an inactive row, composed in the
   coupling) and its test list; [Which weights the predicate sees](../design/empty-leaf-veto.md#which-weights-the-predicate-sees);
   the class row in the R15 inventory of [error-style.md](../design/error-style.md).
   Rd: `dbartsData`'s `counts`, `dbarts`'s multinomial item, `bart`'s count-matrix text (n_i >= 0, the class
   sentence, residuals -p, loglik 0 and dropping those columns before loo), `dbartsSampler`'s `setCounts`. NEWS
   per call 10.
7. Landing: this plan's Landing note, its INDEX row's status, the TODO item removed. (The INDEX row itself must
   ride the commit that adds this plan, or doc-freshness fails.)

## Verification

- `R CMD INSTALL --preclean -l <lib> .`, then `R_LIBS=<lib> Rscript -e 'tinytest::test_package("dbarts")'`: zero
  failures, the new file's count reported.
- `cd tests/cpp && make && ./test_bartcore`: "all tests passed"; then the sanitizer build and the R-loaded ASAN
  run of the new tinytest file, per [Gate hygiene](README.md#gate-hygiene).
- Neutrality: `R CMD INSTALL --preclean --configure-args=--enable-reference-build -l <ref> .`, then
  `R_LIBS=<ref> Rscript benchmarks/R/multinomial-equivalence.R compare
  benchmarks/baselines/multinomial-equivalence-80b1c8d4.rds --bitwise`: every scenario "identical draws", no
  "max |z|".
- Discrimination, each mutation run and reverted: restoring the unconditional PG draw for empty rows fails the
  bitwise tinytest and cpp arms and the fresh-x arm; reading `activeRows_` in `formForestResponse` fails them
  too; reading `activeRows_` in `formForestVetoWeights` fails the `$sampleTreesFromPrior()` arm (0.59 in the
  critique) and passes every other, which is why that arm exists. `touch` each file after reverting.
- `air format --check .`; `lintr::lint()` on each touched R file; `Rscript tools/check-rc-codoc.R .`,
  `Rscript tools/check-win-drift.R .`, `Rscript tools/check-doc-freshness.R .`, each on its own exit status; the
  NEWS parse check; `R CMD check --as-cran` on a tarball built from a clean copy.
