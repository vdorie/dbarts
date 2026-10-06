# written-surface: a model of several forests is written one way, forest by forest

Status: PLANNED (dec-B266 to dec-B275; dec-B246 as revised).

agent: three pushes. Push 1 (two arguments): sonnet implementer, opus reviewer. Pushes 2 and 3 (the
grammar): opus implementer for the R code, sonnet for the respelled tests and the help once the code is
fixed, opus reviewer. The reason for opus on the grammar: every hard part is code read as code
(`substitute()` through `do.call`, the package's own evaluator that resolves `forest` by bare name, formula
environments, `terms()` with `.` and `-`, `model.frame` with `subset`, an na.action and stored `terms` at
new rows), and a slip there is a model fitted in silence, not a failure. The design's prototype, written
with care, still accepted `forest(x1, x2)` as a multiplier and `forest(x1, dose, 30)` as a size until a
critic probed it.
rng: stated per text, as the process contract defines the classes.
- NEUTRAL for every call in inst/tinytest, benchmarks/R, the baselines, the exact gates, the four snapshot
  files and the four consumer packages: each is either left as written or respelled into the same model,
  seed for seed. Proved by the pair script of Verification (34 of 35 pairs identical when the design's
  prototype stood in for the slice; the 35th is the one text below that is meant to move).
- POSTERIOR-CHANGING for six texts, each new in 1.0-0, each accepted before and after and each absent
  from every test, baseline, gate and consumer (searched at the tip). Push 3 carries all six.

  | text | before, measured at the tip | after |
  |---|---|---|
  | `basis = ~ dose + age` | one column, the sum | two columns; the model the tip fits for `~ cbind(dose, age)` |
  | `basis = ~ dose - 1` | one column, dose minus one | one column, dose; the tip's `~ dose` |
  | `basis = ~ 1 + dose` | one column, dose plus one | a constant column and dose; the tip's `~ cbind(1, dose)` |
  | `basis = dose` in a formula term, with another `dose` where the formula is written | the caller's `dose` (draws 3.04 apart from the data's) | the data's; the tip's `~ dose` |
  | `basis = dose` in a `forests` list, with another `dose` where the list is written | the caller's `dose` | the data's; the tip's `~ dose` |
  | `scale()` or `poly()` in a formula term under `subset` | centre and scale of the rows kept | of every row, as `lm`; equal to rounding to the constants written out |

  Each becomes a model another text fits at the tip, so its oracle is that text's draws and no exact gate
  is re-derived. The design listed seven; two are not changes of meaning (Calls made in planning).
- Refusals and acceptances, no draws involved: texts the tip accepts and the slice refuses, and the
  reverse, are listed under "Refused forms" and "Context".
window: pre-release, first of the forest-prior-args slices: the sd unit (N), the kind by class (K), the
multiplier law (L) and everything after them write their tests and help in this spelling. Serial with any
other work in [R/formulaTerms.R](../../R/formulaTerms.R), [`forest`](../../R/model.R),
[`resolveForests`](../../R/model.R), [`expandForestBasis`](../../R/model.R),
[`replayForestBasis`](../../R/model.R), [`setForestBasis`](../../R/dbarts.R) or the multi-forest block of
[`resolveSamplerSpec`](../../R/spec.R). No engine, bridge or header file changes, so it shares only two
files with [cut-points-undo.md](cut-points-undo.md) and [leaf-conversions.md](leaf-conversions.md)
(R/dbarts.R and man/dbartsSampler-class.Rd, other functions and items) and may land before, between or
after them. bartCause's edit lands the day push 1 does.
budget: ~4800 lines changed (R ~1750; tinytest ~2330, of which 587 removed with the two rewritten files,
~680 written in their place, ~900 new and ~160 changed elsewhere; help, vignette and NEWS ~400;
benchmarks ~20; design note, architecture and TODO ~300), upper figure 7000. By push: 1 ~520 (upper 800),
2 ~1800 (upper 2700), 3 ~2480 (upper 3600). The design estimated 3050 and planned for 5500; plans have run
1.5-2x low, and its prototype is 2020 changed lines of R alone with no test or help.

## Assumed, pending the maintainer

Three of the four points this plan was asked to assume were ruled while it was written or just after.

- Ruled, dec-B273: a star at the top of a basis is refused by name. Planned so; no alternative carried.
- Ruled, dec-B274, against the assumption given: a formula whose forests all have a basis is ACCEPTED, and
  a forest's default tree count and tree prior are to go by its kind. This slice does not refuse such a
  formula and does not build the defaults by kind: the forests of such a formula are taken in written
  order, exactly as a `forests` list of them is taken today, until the slice that builds dec-B274 lands
  (Out of scope). One test pins that (step 2.6); that slice turns it around.
- Ruled, dec-B275, as assumed: `extract(type = "leaf.prior.sd")` on several forests returns a list, one
  element per forest named by label, each a vector named by column. This slice builds no part of it; it
  fixes the two sets of names.
- Assumed (d): `setForestBasis(updateBasisScale = )`, `leaf.prior = normal(sd = )` on `forest()` and
  selecting a forest by label are not in this slice. If they are to follow the release instead of
  landing in their later slices: nothing here changes, since this slice already leaves `leaf.prior =`
  and `label =` as R's unused-argument error and a character `forest` refused, which is what makes each
  an addition. Only "Out of scope" changes where it points.

## Goal

A model of several forests is written as a sum of `forest()` terms, each one function of the predictors
named inside it: `y ~ forest(x1 + x2) + forest(x1 + x2, basis = dose, sd = 2)`. A forest's multiplier is
its `basis` argument, written as the right-hand side of a model formula with no tilde, and it means what
it means in `lm`: `dose + age` is two columns, `factor(z)` one per level, `scale(age)` and `poly(dose, 2)`
are rebuilt at new rows from the fitted rows. The same `forest()` is written in a `forests` list. The two
spellings that go, `z:forest(x)` and `update.amplitude`, are refused or unknown, and every spelling that
could later state a size per column is refused by name now, so that adding one changes no call.

## Context

Measured at the tip on the shipped build (R 4.6.1): 150 rows, 15 trees, 25 + 25 sweeps, one chain unless
said.

- The suite. 224 test files, 14274 assertions, none failing.
- Two spellings of a multiplier. `y ~ x1 + x2 + z:forest(x1 + x2)` and
  `y ~ x1 + x2 + forest(x1 + x2, basis = ~ z)` are one fit, draw for draw; so are `(dose + age):forest(x1)`
  and `basis = ~ cbind(dose, age)`. The colon takes a name, `factor(name)` or a parenthesised sum of
  names; `scale(age):forest(x1)` is refused ("is not a supported forest() modulator") and `z * forest(x1)`
  is refused naming the colon. [`walkFormulaTerms`](../../R/formulaTerms.R),
  [`desugarBasisOperand`](../../R/formulaTerms.R).
- The first forest cannot be written. `y ~ forest(x1 + x2) + forest(x1, basis = ~ dose)`,
  `y ~ forest(x1 + x2)` and `y ~ forest(x1, basis = ~ dose) + forest(x2, basis = ~ age)` all stop with "a
  formula whose only right-hand-side content is a forest() term has no predictors left for the base
  forest after rewriting". Plain terms beside `forest(x1 + x2)` stop with "forest 2 needs a 'basis'".
- A term's predictors. The unnamed argument is names joined by `+`, each among the plain terms:
  `forest(log(x1), ...)` and `forest(x3, ...)` beside `x1 + x2` are refused; a second unnamed argument is
  refused ("gives more than one unnamed argument"). `vars = c("x1")` and `vars = 1` by name are accepted
  in a term. [`processHit`](../../R/formulaTerms.R), [`finalizeTermForests`](../../R/formulaTerms.R).
- Around the terms. `offset(off)` and `0 +` among the plain terms are accepted; `- 1` after a term is
  refused ("must appear as a top-level additive term"); `.` with removals beside a term fits the columns
  left. `test` and `forests` beside a term are refused by name. Gaussian, probit and logistic take a term;
  aft, ordinal, nbinom, hazard and multinomial refuse it by name.
- Who takes a term. `bart` (and its alias `bart2`) and `dbarts`. `bartBT` stops with the refusal about
  `test`, which the caller did not give; `xbart`, `rbart_vi` and `dbartsData` stop inside R's model frame
  ("invalid type (list) for variable 'forest(x1, basis = ~z)'"); `pdbart` stops at its predict.
- A basis is R code. Columns of `data@bases[[2]]` for `basis = ~ <text>` in a list: `dose + age` one (the
  sum, 39.41 in row 1); `dose * age` one (27.79, the product); `dose - 1` one (-0.28); `1 + dose` one
  (1.72); `dose / 30`, `30 * dose`, `dose - age`, `dose^2`, `z + zl` and `offset(dose) + age` one each;
  `cbind(dose, age)` two. `dose:age` is refused ("'basis' must have the same length as 'y'", with R's
  warning that only the first element was used); `dose + factor(z)` is refused ("a 'basis' cannot be NA");
  `normal(dose, sd = 30)` stops with "could not find function". `lm(y ~ 0 + dose / 30)` is R's "invalid
  model formula in ExtractVars". [`evaluateForestBasis`](../../R/model.R).
- Column names. `~ cbind(dose, age)` gives `dose`, `age`; `~ poly(dose, 2)` gives `1`, `2`; `~ dose`,
  `~ scale(age)`, `~ factor(z)`, a character and a logical column give none. A value keeps what it has,
  partial names included: `cbind(1 - z, z)` is `""`, `"z"`, and `cbind(a = , a = )` is accepted.
- Where names are found. A formula term's `basis = ~ dose` is built in the fit's model frame, the data
  first. Every other form is evaluated where it was written and handed over as a value: `basis = dose`
  with `dose` only in `data` is "object 'dose' not found" at both doors, and with a second `dose` beside
  the call it is that one. `basis = "dose"` stops with "a 'basis' factor must have at least two levels".
  A constant beside the formula, `~ I(dose / k)`, is refused ("variable lengths differ (found for 'k')").
- `forest()`'s formals are `basis, vars, n.trees, base, power, sd, interactions, blocks,
  amplitude.prior.variance, update.amplitude`. In a list `forest(~ factor(z))` is a positional basis and
  `forest(~ factor(z), c("x1"))` a basis and a selection. No test, benchmark, help page or consumer writes
  a positional argument in a list (searched). A selection by value collapses a repeat: `c("x1", "x1")` is
  column 1.
- `sd`. At creation and through `$setLeafPrior(forests = )` alike: `"2"` is read as 2, `TRUE` and
  `factor("2")` as 1, `list(2)` and `matrix(2)` as 2, `as.Date("2026-01-01")` as 20454, `c(a = 1)` as 1
  with the name dropped; `c(1, 2)`, `NA`, `NaN`, `Inf`, 0 and -1 are refused with one text ("must be a
  single positive finite number"). The front door's `normal(sd = )` refuses a string and `c(1, 2)` and
  accepts `TRUE`, a Date, a list, a matrix and a named number.
  [`validateForestSd`](../../R/model.R), [`validateLeafSd`](../../R/model.R).
- The coefficient. `update.amplitude = FALSE` on both forests of a two-level factor model holds the
  amplitudes at 1, 0, 1. `amplitude = fixed()` stops with "could not find function \"fixed\"";
  `leaf.prior =`, `tree.prior =` and `label =` on `forest()` are R's unused-argument error.
- One written forest. `forests = list(forest(sd = 2))` and `list(forest(update.amplitude = FALSE))` are
  refused ("a single-forest 'forests' has none"); `list(forest(n.trees = 7))` gives 7 trees. In
  `bart(y ~ x1 + x2 + forest(x1, basis = ~ z), n.trees = 15)` the forests have 15 and 50 trees.
- All forests multiplied. `forests = list(forest(basis = ~ dose), forest(basis = ~ age))` is accepted;
  the first has the control's 15 trees and the second 50; the two swapped give draws 1.62 apart.
- Labels. A formula's forests and an unnamed list have none (`attr(fit, "forest.labels")` is `NULL`);
  `list(forest(), dose = forest(...))` gives `""`, `"dose"`; two list names alike and a second forest
  named `forest1` are accepted. `getLeafPrior("dose")` is refused: a forest is selected by position.
  `extract(type = "leaf.prior.sd")` on three forests is a vector named `forest1`, `forest2`, `forest3`.
- Rows. Under `subset`, a formula term's `~ scale(age)` is centred on the rows kept and a list's on every
  row. A factor column with a level no row has is refused at both doors, as is one whose level `subset`
  empties; `~ factor(g)` under the same `subset` is accepted, `factor()` being called on the rows kept. A
  missing value in a basis is refused. `predict` rebuilds a formula term's `scale()` at new rows; for a
  basis written in a list it stops ("no off-sample basis").
- `$setForestBasis(2, cbind(age = , dose = ))` on a forest created with `dose`, `age` is accepted and the
  stored names become `age`, `dose`: columns are taken by position and names follow the value.
- What is written where, counted in the tree at the tip:

  | spelling | where |
  |---|---|
  | a colon or star against `forest()` | 37 lines in 9 test files (28 in test-formula-terms.R), 18 assertions whose own call spells it; 3 lines in 3 benchmark scripts, one of them the equivalence scenario `bart2twoforest`; 8 lines of 2 help pages; one NEWS item |
  | `update.amplitude` | 17 lines in 6 test files (one a comment), 2 assertions whose own call spells it; 12 lines in 6 benchmark scripts, the four bcf exact gates and the bcf equivalence script among them; 2 lines of man/forest.Rd |
  | the two files that test the old grammar | test-formula-terms.R (263 lines, 78 assertions), test-forest-basis-terms.R (324 lines, 75 assertions) |
  | the 15 test files with a colon or `update.amplitude` | 1484 of the suite's 14274 assertions |
  | `basis = ~` | 183 times on 178 lines of 43 test files: `factor()` of a name 143, a bare name 32, `scale()` 3, one each of `cbind(a, b)`, a subscripted matrix and `rep(0, n)`, and 2 inside pinned message texts. The first three kinds, 178 of them, mean the same afterwards |
  | `basis = ` with no tilde, in a list | 20 lines in 7 test files, each a value today and code afterwards; none is written beside a data column of the same name, so each is the same model |
  | help | man/forest.Rd (all of it, 3 example lines), man/bart.Rd (three argument items, the section "Formula Terms", one example), man/dbarts.Rd, man/dbartsSampler-class.Rd (two items), man/dbartsForests.Rd (one sentence); one vignette line |
  | bartCause (dbarts-1.0, dc4e397) | R/bcf.R: `update.amplitude = update.a` and `= update.b`, 2 lines; tests: `update.amplitude = TRUE` on 2 lines of 2 files |
  | stan4bart (a9d081b), treatSens (dbarts-1.0, aecec71), bairrtt (3f57f61) | `dbarts::dbartsForests$forest(n.trees = )` and nothing else of this surface |

- Checked cites that name what this slice rewrites: ten symbols of R/formulaTerms.R, R/model.R and
  R/dbarts.R cited from ten documents under docs/, and the fragments
  ["expectSameForest"](../../inst/tinytest/test-formula-terms.R) and
  ["needs at least two forests"](../../inst/tinytest/test-bcf-creation.R).
- The design this plan follows is recorded in docs/design/written-surface.md, which push 3 writes; until
  then the rulings are its record. Its prototype is evidence that the grammar can be built, not code to
  copy: it carries option switches and stand-ins for two later slices.

## The grammar

    forest(vars = NULL, basis = NULL, sd = NULL, n.trees = NULL, base = NULL, power = NULL,
           amplitude = NULL, interactions = NULL, blocks = NULL, amplitude.prior.variance = NULL)

Those are the formals after push 3. `base` and `power` become `tree.prior` and `leaf.prior` arrives in
later slices; `amplitude.prior.variance` goes with the multiplier law. Only `vars` may be given unnamed.
"Door" below is a formula or a `forests` list; "plain terms" are a formula's terms outside any `forest()`.

| written | means |
|---|---|
| `y ~ x1 + x2`, `y ~ forest(x1 + x2)` | one forest; the same single-forest fit |
| `y ~ forest(x1 + x2) + forest(x1, basis = z)`, `y ~ x1 + x2 + forest(x1, basis = z)` | two forests; one model |
| `y ~ x1 + x2 + forest(x3, basis = z)` | forest 1 splits on x1, x2 and forest 2 on x3; the fit's predictors are x1, x2, x3 |
| `y ~ x1 + x2 + forest(basis = z)` | forest 2 splits on every predictor of the fit |
| `y ~ forest(x1, basis = z) + forest(x1 + x2)` | the forest with no basis is forest 1 wherever written |
| `y ~ forest(x1, basis = a) + forest(x2, basis = b)` | accepted (dec-B274); forests in written order, the first taking the fitting function's tree count and tree prior, as a list of the two does today |
| `y ~ offset(off) + forest(x1 + x2) + ...`, `y ~ 0 + ...`, `... - 1` | the offset and the intercept token are the fit's, as in `y ~ offset(off) + x1 + x2` |
| `y ~ . - z + forest(x1, basis = z)`, `forest(. - z)`, `forest(log(x1) + factor(g))` | a forest's predictors are the right-hand side of a model formula of its own: `.` is every column but the response, `-` removes a term |
| `forests = list(forest(), dose = forest(x1, basis = dose))` | the list door; there the first argument selects among the fit's predictors, by terms or by a value of names or positions |
| `basis = dose + age` | two columns, one coefficient each, one forest |
| `basis = I(dose + age)`, `log(dose)`, `scale(age)`, `poly(dose, 2)` | one term each, as in `lm`, at new rows too |
| `basis = factor(z)`, a character or logical column | one coefficient per level, every level that a kept row has |
| `basis = dose:age` | one column, the product, named `dose:age` |
| `basis = 1 + dose`; `0 + dose`, `dose - 1` | a constant column and dose; dose |
| `basis = ~ dose + age`, `basis = b` where `b <- ~ dose + age` | the same as without the tilde |
| `sd = 2` | one unnamed number, the size of the forest for every column of its basis; its unit is unchanged by this slice |
| `amplitude = fixed()`, `fixed(1)` | the coefficient is held; what `update.amplitude = FALSE` is today |

One written forest and the fitting function's own arguments:

| on the one forest with no basis | this slice | later |
|---|---|---|
| its terms | its predictors | |
| `n.trees` | served; refused by name when `bart()`'s own `n.trees` is given too | `dbarts()`'s own count: the control-migration arc (dec-B241) |
| `sd` | refused, naming `leaf.prior = normal(sd = )` on the fitting function | served once the two share a unit (slice N) |
| `amplitude`; a `basis` with no second forest | refused by name | |
| `interactions`, `blocks` | served; both given refused, as at the tip | |

## Refused forms, with their texts

Base R's style: no capital, no closing stop, the argument named in single quotes, the caller's own text
quoted, the form to write given. `<f>` is the forest's text; the push that adds each is in brackets.

    [2] 'dose:forest(x1 + x2)': a forest() is not crossed with another term; a forest's multiplier is its 'basis' argument: write forest(x1 + x2, basis = dose)
    [2] forest() takes one unnamed argument, the predictors the forest splits on, joined by '+' as forest(x1 + x2); every other argument is given by name: a multiplier is 'basis =' and a size is 'sd ='
    [2] forest()'s first argument is the predictors the forest splits on, written without '~', as forest(x1 + x2); a multiplier is 'basis ='
    [2] forest()'s first argument, 'vv', holds a formula; write its terms in place, as forest(x1 + x2), or give the predictors' names as a character vector
    [2] forest()'s first argument selects predictors of the fit and names 'x1' more than once; name each once. A multiplier is given as 'basis ='
    [2] 'nosuch' is not a predictor of this fit (x1, x2, x3); here a forest()'s first argument selects among them
    [2] '<f>': in a formula a forest's predictors are written as terms, as forest(x1 + x2), or named, as forest(c("x1", "x2")); a position is for a 'forests' list
    [2] the formula has plain terms (x1 + x2) and a forest() with no basis ('<f>'): each is the forest with no multiplier, and a model has one. Write the plain terms inside that forest(), or give it a 'basis'
    [2] the formula has 2 forest() terms with no basis; a model has one forest with no multiplier, and every other forest states a 'basis'
    [2] the formula names no predictors: write them as plain terms or inside a forest(), as forest(x1 + x2)
    [2] '<f>': an offset() is the fit's and no forest's; write it beside the forests, as y ~ offset(off) + forest(...)
    [2] '<f>': an intercept term (1, 0 or - 1) is the fit's and no forest's; write it beside the forests
    [2] a forest() term must appear as a top-level additive term, not inside '<text>'
    [2] 'n.trees' is given to the fitting function and to the forest with no basis ('<f>'), which are the same count; give one
    [2] a multi-forest model needs at least two forests, and this call's 'basis' declarations resolve to 1: a forest with a 'basis' stands beside another forest. Write the forest with no multiplier too, as y ~ forest(x1 + x2) + forest(x1 + x2, basis = z1) or forests = list(forest(), forest(basis = z1)), or use a single forest with linear() leaves; otherwise drop the basis
    [2] xbart() does not take a forest() term ('<f>'): the forests of a model are written in the formula of bart() or dbarts(), or in a 'forests' list
    [1] forest 'sd' must be a single number, not 2: a forest states one sd, for every column of its basis; to size the columns differently, rescale them in 'basis', as I(dose / 30)
    [1] forest 'sd' must not be named ("dose"): it is one number, for every column of a basis
    [1] forest 'sd' must be a number, not a string            (a logical, a factor, a Date, a list, a matrix)
    [1] forest 'sd' must not be NA; leave it out for the default            forest 'sd' must be positive and finite
    [1] forest 'sd' is stated for a forest of a model of several; a model of one forest states its size as the fitting function's leaf.prior = normal(sd = )
    [1] 'amplitude = fixed(2)': a held coefficient takes the value its forest's shape gives it, and fixed() takes no other here; write fixed(), and state the forest's size with 'sd'
    [1] a forest's 'amplitude' must be fixed(), which holds its coefficient, or left out, which draws it
    [1] 'amplitude' is the law of the coefficient that a model of several forests gives each of them; a model of one forest has none
    [3] 'basis' does not take '*' between its terms ('dose * age'): in a model formula it is both columns and their product. Write I(dose * age) for the product alone, or dose + age + I(dose * age) for all three
    [3] 'basis' term 'dose/30' divides a column by a number, which a model formula does not take; write I(dose/30) for the rescaled column
    [3] 'basis' term '30 * dose' multiplies a column by a number, which a model formula does not take; write I(30 * dose) for the rescaled column
    [3] 'basis' does not take '-' between its terms ('dose - age'): in a model formula it removes a term and subtracts nothing; write the arithmetic inside I(), as I(dose - age)
    [3] 'basis' does not take '^' between its terms ('dose^2'): in a model formula it crosses terms and raises nothing to a power; write the arithmetic inside I(), as I(dose^2)
    [3] 'basis' does not take '/' between its terms ('dose/age'): in a model formula it nests one term in another and divides nothing; write the arithmetic inside I(), as I(dose/age)
    [3] 'basis' does not take '%in%' ('dose %in% age'): in a model formula it nests one term in another
    [3] 'basis' does not take '|' ('dose | z')            'basis' does not take '.': name the columns the forest is multiplied by            'basis' does not take an offset() term
    [3] 'basis' has the number 2 as a term; only 1, a constant column, and 0, none, are terms. Write arithmetic on a column inside I()
    [3] 'basis' has the constant TRUE as a term; a term is a column, and only 1, a constant column, and 0, none, are written as numbers
    [3] 'basis' (0 + 1) is a constant column and nothing else, which multiplies the forest by a constant: leave 'basis' out for the forest with no multiplier
    [3] 'basis' term 'cbind(dose, age)': the columns of a basis are separated by '+'; write dose + age
    [3] 'basis' term 'normal(dose, sd = 30)' calls normal(), which is not a column: a prior is not stated on a basis term. Give the forest's size as 'sd' and its coefficient's law as 'amplitude'
    [3] 'basis' (dose + factor(z)) mixes a factor with other terms: a basis is one factor, a character or a logical vector, with one coefficient per level, or numeric columns, with one each. Give the other terms a forest() of their own, or write the factor as numbers
    [3] 'basis' (<text>) has a term that is a logical matrix: a basis is a factor, a character or logical vector, or numeric
    [3] 'basis' (dose > 0) is TRUE on every row the fit keeps, so its other level has no observations and the forest would be multiplied by a constant
    [3] 'basis' is 'v', which holds the string "dose": a basis is a column, not its name. Write basis = dose, or from a program do.call(forest, list(basis = as.name("dose"))), or hand over the column itself
    [3] 'basis' is the single value 2: a basis has a value for every observation
    [3] 'basis' (zBasis[idx, ]) must have the same length as the data: it has 20 rows and the data 40; a basis covers every row of the data and is cut by 'subset' and the na.action with it
    [3] 'basis' (nosuch): object 'nosuch' not found in 'data' or where forest() was called
    [3] 'basis' (<text>) has two columns named "a"; the columns of a basis are told apart by name
    [3] a 'basis' formula must be one-sided, as ~ dose
    [3] 'basis' has the columns of forest 2's basis in another order (age, dose; the forest has dose, age): columns are taken by position, so give them in the forest's order
    [3] 'forests' names two forests "a"; a label names one forest
    [3] 'forests' names forest 2 "forest1", which is the label an unnamed forest 1 has; choose another name

An argument that does not exist gets R's own error: `unused argument (update.amplitude = FALSE)` from push
1, and `by =`, `label =`, `leaf.prior =` and `tree.prior =` as today. Kept as they are: the refusals of
`test` and of `forests` beside a term, the families that take no term, a term on the left-hand side, and
"a 'basis' cannot be NA". Between pushes 2 and 3 the first text ends `basis = ~ dose`, the tip's spelling
of the same operand, and the at-least-two text writes `basis = ~ z1`; push 3 drops the tildes.

## The label rule

1. Every forest has a label, fixed when the sampler is created. `$setForestBasis` never changes it.
2. It is the list name where one is written.
3. Otherwise it is the text of the forest's basis: R's deparse of the code, after a name that holds a
   formula is replaced by that formula and the tilde is taken off. The same text gives the same label at
   either door: `dose`, `scale(age)`, `I(dose/30)` however it was spaced.
4. Otherwise (no basis, a basis handed over as a value, a basis from the data object) it is `forest<i>`,
   i the forest's position.
5. Labels are unique. Two list names alike, and a list name of the form `forest<digits>` that is not its
   own position's, are refused. A basis text that repeats an earlier label takes `make.unique`'s suffix
   (`dose`, `dose.1`), as the names of a data frame do; a formula has no other way to name a forest.

In this slice the labels are recorded (the `labels` entry of the forests' configuration,
`attr(fit, "forest.labels")`), `$setLeafPrior(forests = )` checks a list's names against them, and they
select nothing: a forest is still chosen by position, and the per-forest margins keep `forest1`,
`forest2`. Selecting by label is a later slice's; the rule is fixed now because the reader's and
`extract`'s lists will be named by it for good (dec-B272).

## The capture rule, and the help's warning

- `forest()` reads `vars` and `basis` unevaluated and keeps the code with the place it was called: the
  formula's environment for a term, the caller's frame for a list (not the package's own evaluation
  frame, which binds `forest`, `fixed`, `interactions` and `blocks`).
- A name in them is looked up in `data` first and then there, as `lm` does. With no data frame (the
  matrix interface, `dbartsSpec`) every name is looked up where the call was made.
- What is handed over as an object is a value and is used as it is, its class deciding its kind: through
  `do.call`, or in a call built with the value in it, as bartCause builds its forests. A name or a call
  handed over (`as.name("dose")`, `quote(scale(age))`) and a one-sided formula are code.
- A name that is not a column of `data` and holds a one-sided formula stands for that formula. A tilde
  written in place is accepted for good. A two-sided formula is refused.
- Code that gives `NULL`, and a name holding `NULL`, state no basis.
- Three slips of a caller who builds forests in a loop have their own texts: a string as a basis, a
  single number as a basis, a formula held in `vars`.
- Every other argument is evaluated as today, with `fixed`, `interactions` and `blocks` resolved by bare
  name inside the argument.

The help for `forest` carries, under its `basis` item, what `?subset` carries for its own arguments:

    'vars' and 'basis' are read as code, not as values. A name in them is looked up in 'data' first and
    then where forest() was called, so a column of 'data' hides a variable of yours with the same name,
    whether that variable holds a column, a set of names or a one-sided formula. In a function that
    builds forests, hand values over instead: do.call(forest, list(basis = w)) for a column you hold,
    list(basis = as.name(nm)) for a column of 'data' by its name, list(basis = f) for a one-sided
    formula, and list(vars = nms) for predictors by name. A number found where forest() was called, the
    k of I(dose / k), is looked up again by predict, as in lm: change it between the fit and the
    prediction and the prediction changes, with no message.

## What dec-B272 requires of this slice

Each is a requirement with the test that holds it; "what fails it" is under the step named.

| requirement | test | step |
|---|---|---|
| Only the first argument of `forest()` may be unnamed, so a later formal is an addition anywhere | `forest(x1, x2)`, `forest(x1, dose, 30)`, `forest(x1, ~ z)` in a formula, `forest(x1, dose)` in a list and `do.call(forest, list("x1", w))` are each refused with the one-unnamed-argument text | 2.2 |
| An sd longer than one, or named, is refused wherever a number can be stated | `c(1, 2)`, `c(dose = 1, age = 2)` and `c(a = 1)` at creation, through `$setLeafPrior(forests = )` and in the front door's `normal(sd = )` | 1.2 |
| A basis term that calls a reserved name is refused: `normal`, `fixed`, `student`, `cauchy`, `linear`, `gp`, `cgm`, `dart`, `chisq`, `chi`, `invchi`, `forest`, `varianceForest` | each of the thirteen as a call at the top of a basis, a caller's own `normal()` and `student()` included; and a caller's `weights2(dose)` is a column, pinned so that a wrapper under a new name fails this test | 3.2 |
| A term divided or multiplied by a number is refused; a star between terms is refused (dec-B273) | `dose / 30`, `30 * dose`, `dose * 30`, `dose * age`, `poly(dose, 2) * age`; `I(dose / 30)` is accepted with the column named `I(dose/30)` | 3.2 |
| Columns are named as `coef(lm(y ~ 0 + <basis>))` names them, by a complete rule | 18 shapes against `lm`; a value's names; two columns alike refused | 3.3 |
| Labels follow a complete rule | the literals of "The label rule" at both doors | 3.6 |
| One sd on several columns means the same size for each, now and later | the help's `sd` item says so; a two-column basis with `sd = s` is, draw for draw, the tip's `~ cbind()` of the two with `sd = s` (pairs 03 and 14) | 3.9 |
| A basis written as code is cut by `subset` before its levels are looked at, at both doors | the same text at the formula door and in a list, under `subset` and under a dropped missing response: identical bases and draws, four pairs | 3.4 |
| A star is settled before the release | dec-B273: refused; the test above | 3.2 |

The critic's sixth condition, that the engine take one stated sd per column, is the multiplier-law
slice's, in tests/cpp; nothing in R states one and this slice adds no bridge read.

## Constraints

- Every fit that the suite, the benchmarks, the baselines, the gates and the consumers make draws what it
  draws at the tip. No baseline is re-recorded and no snapshot file regenerated.
- No change to src/, to the flat C header, to the stored state or to the eight numbers per forest that
  the bridge reads; `fixed()` is the eighth number's zero, as `update.amplitude = FALSE` is.
- The unit of `sd` does not change, and the help states it as it is today, with figures in that unit.
- The kind of a multiplier is not re-decided: what expands to level columns today does afterwards. The
  data door ([`validateForestBases`](../../R/data.R)), `data@basis.levels`, the swap table of
  `$setForestBasis` and a value's emptied level are slice K's.
- The reader and `extract` keep today's shapes. None of the six returns the design pins as literals is
  built here: the multiplier-law slice builds the per-column `sd` entries and `extract`'s list, after the
  unit slice so that they are written once in the response's units, and the forest-leaf-prior slice makes
  `leaf.prior` a `normal()` object. `$setLeafPrior(forests = list(forest(sd = )))` stays the writer.
- An all-multiplied formula is accepted and no top-level argument is refused for it (dec-B274's refusal
  and defaults are a later slice's).
- Nothing of `tree.prior`, of `leaf.prior` on `forest()` or of `updateBasisScale` rides along.
- `dbartsForests$forest(n.trees = )` called outside any argument keeps working: three consumers call it.
- The two mutation-battery anchors in [`writeForestSpreads`](../../R/dbarts.R) are not edited.
- Each push leaves the help saying what the code does.
- Base R calls stay within DESCRIPTION's R floor.

## Pushes

Three, each gated on its own (Verification) and each a coherent tip.

1. Two arguments (dec-B271, dec-B269). `amplitude = fixed()` in place of `update.amplitude`; one check of
   `sd`. bartCause's edit lands the same day. NEUTRAL.
2. The forests of a formula (dec-B267, dec-B268). The colon goes, the first forest may be written, a
   forest's predictors are its own terms, one unnamed argument. A basis is still written and read as at
   the tip, tilde and all, so the help's `basis` item is untouched and true. NEUTRAL.
3. The basis (dec-B270, dec-B272, dec-B273). No tilde, a model formula, names, labels, both doors alike,
   `predict` at both doors. POSTERIOR-CHANGING for the six texts; the design note lands here.

Pushes 1 and 2 may be reviewed together if the coordinator prefers two reviews; push 3 stays apart.

## Steps

"Fails today" is what the tip does where the test expects otherwise. New functions are named for the
reader's sake; the implementer may name them otherwise.

### Push 1: two arguments

1.1 [`forest`](../../R/model.R) gains `amplitude` and loses `update.amplitude`.
    [`validateForestKnobs`](../../R/model.R) takes `fixed()`, `fixed(1)` and a bare `fixed`, refuses any
    other `fixed(v)`, a `normal()`, a string and a logical with the two texts, and stores a flag;
    [`forestParams`](../../R/model.R) reads the flag for its eighth number;
    [`resolveForests`](../../R/model.R) refuses `amplitude` and `sd` on a single forest with the two
    one-forest texts. [`FOREST_ARGUMENT_VOCABULARIES`](../../R/model.R) lets `fixed` resolve by bare name
    in `forests`, in a term ([`processHit`](../../R/formulaTerms.R)) and in the writer's `forests`.
    [`resolveForestSpreads`](../../R/dbarts.R) names `amplitude` in its fixed-at-creation text.
    Tests, a new file test-forest-arguments.R: `fixed()` on both forests leaves the amplitudes at their
    starting values over 50 sweeps and gives other draws than the drawn model (fails today: could not find
    function); `fixed`, `fixed(1)`, `fixed(1L)` and `NULL` accepted, `fixed(2)`, `fixed(c(1, 1))`,
    `fixed(TRUE)`, `normal()`, `"fixed"`, `FALSE` refused; `update.amplitude` is R's unused-argument error
    at both doors (fails today: accepted); `fixed()` resolves in a call built where nothing of dbarts is
    bound, as bartCause builds it; `fixed()` is the constructor whatever the caller has bound to the name.
1.2 One check of a stated sd: a bare numeric or integer of length one, unnamed, finite, positive,
    nothing coerced. [`validateForestSd`](../../R/model.R) is that check, with the texts marked [1].
    [`validateLeafSd`](../../R/model.R) refuses, ahead of the checks it has, what the forest's check
    refuses and it accepts today (a logical, a factor, a Date, a list, a matrix, a named number); its
    `NULL`, its `invchi()`, its other texts and the check it shares with `k` stay. Tests, same file: the
    15 values of Context at creation, through `$setLeafPrior(forests = )` and in `normal(sd = )`, one
    verdict at the three (fails today: seven accepted on a forest, and five of them at the front door);
    the text for each; `2L` accepted.
1.3 Respell: the 16 code lines in 6 test files and the 12 lines in 6 benchmark scripts become
    `amplitude = fixed()` or `amplitude = if (...) NULL else fixed()`; the two pins of
    ["single-forest 'forests' has none"](../../inst/tinytest/test-bcf-creation.R) and the three lines
    (seven assertions as run) that pin
    ["single positive finite number"](../../inst/tinytest/test-multiforest-leaf-prior-writer.R) take the
    new texts.
1.4 Help: man/forest.Rd's usage, its `update.amplitude` item rewritten as `amplitude`, one sentence in
    `sd` (one unnamed number; a longer or named one is refused);
    [`dbartsSampler$setLeafPrior`](../../man/dbartsSampler-class.Rd) where it lists what the writer
    refuses. docs/design/public-surface.md's one line.
1.5 bartCause, same day (its own commit on dbarts-1.0): in R/bcf.R each `update.amplitude = update.x`
    becomes `amplitude = if (update.x) NULL else quote(fixed())`, behind a check that `update.a` and
    `update.b` are each `TRUE` or `FALSE` (dbarts made that check until now; bartCause has none); in the
    two test files the line `update.amplitude = TRUE` goes; one new expectation for the check.

### Push 2: the forests of a formula

2.1 [`walkFormulaTerms`](../../R/formulaTerms.R): a hit is a `forest()` call at the top of the `+` chain
    and nothing else. A colon or star with a forest on either side, at any depth, is refused with the
    first text, the rewrite built from the operand as [`desugarBasisOperand`](../../R/formulaTerms.R)
    builds a basis today; that function and [`flattenPlusSymbols`](../../R/formulaTerms.R) go. A trailing
    `- 1` is the fit's. Plain terms that are an `offset()` or an intercept token stay in the fit's formula
    and count as no forest. Tests (test-formula-terms.R, rewritten): eight colon and star shapes refused,
    each with its rewrite (fail today: fitted); `I(forest(x1))`, a forest in a removal and on the
    left-hand side refused; six placements of `offset(off)`, `0 +` and `- 1` give identical draws (fails
    today: `- 1` refused).
2.2 [`forest`](../../R/model.R): formals in the order of "The grammar" less what push 3 adds, `vars`
    first; `vars` is captured as code with the place of the call. A second unnamed argument is refused at
    both doors and through `do.call`. Tests: dec-B272's first row (fails today: in a list the second is
    taken as `vars`, and a first as `basis`); `forest(basis = ~ z, x1 + x2)` accepted, the one unnamed
    argument written second; `dbartsForests$forest(n.trees = 5L)` outside any argument builds what it
    builds today.
2.3 A forest's predictors in a formula. A new reader of a right-hand side (`readRhsTerms`) runs `terms()`
    on each forest's first argument against the data and the response: names, `factor(g)`, `log(x)`, `.`,
    `-`; an interaction is refused by the fit's own message; `offset()` and an intercept token inside a
    forest are refused. The fit's formula is rebuilt once from every forest's terms, the unmultiplied
    forest's first and then the others' as written, each term once, and the predictor matrix is built by
    the code that builds it today. [`finalizeTermForests`](../../R/formulaTerms.R) gives each forest the
    columns of its own terms ([`resolveTermColumns`](../../R/model.R)); a forest naming every column is
    unrestricted. A value in a term may be names (`forest(c("x1", "x2"))`), not positions. Tests:
    `y ~ x1 + x2 + forest(x3, basis = ~ z)` has predictors x1, x2, x3 and forest 2 splits on x3 alone over
    300 sweeps (fails today: refused); `forest(log(x1) + factor(g), ...)` splits on those columns; `. - z`
    inside and outside a forest; `forest(basis = ~ z)` splits on all; no predictors anywhere, a tilde on
    the predictors, a held formula and a position are refused with their texts.
2.4 The forests of a formula. A `forest()` with no basis is the forest with no multiplier: alone it is
    the single-forest fit; with others it is forest 1 wherever written; beside plain predictor terms, or
    twice, it is refused. Tests: `y ~ forest(x1 + x2)` against `y ~ x1 + x2`, and
    `y ~ forest(x1 + x2) + forest(x1, basis = ~ z)` against the plain-terms form and against the two
    terms swapped: identical draws (all fail today: refused); the two refusals with their texts, the
    plain terms named and never an offset; `forest(x1 + x2, basis = ~ z)` alone refused with the
    at-least-two text; `interactions` on the written first forest served, and refused when given at the
    top too.
2.5 The tree count given twice. [`bart`](../../R/bart.R) refuses its own `n.trees` beside one on the
    written forest with no basis. Tests: the refusal; 7 trees with no top-level count; the count of
    another forest beside `bart()`'s is accepted (15 and 7).
2.6 The all-multiplied formula keeps written order. Test, the one that pins dec-B274's interim:
    `dbarts(y ~ forest(x1, basis = ~ dose) + forest(x2, basis = ~ age), d, control = <15 trees>)` is
    accepted (fails today: refused), has predictors x1, x2 and 15 and 50 trees, is draw for draw
    `dbarts(y ~ x1 + x2, d, forests = list(forest(x1, basis = ~ dose), forest(x2, basis = ~ age)))` under
    the same seed, and the two terms swapped give other draws. The help says in one sentence that the
    forest written first takes the fitting function's tree count and tree prior there, and tells the
    reader to state `n.trees` on each such forest.
2.7 The list door's first argument. [`resolveForests`](../../R/model.R) and the multi-forest block of
    [`resolveSamplerSpec`](../../R/spec.R) resolve `vars` as terms over the fit's predictors when its
    names are predictors, and otherwise evaluate it where the call was made and select by names or
    positions ([`resolveModerators`](../../R/model.R)); a repeat in a selection is refused. Tests:
    `forest(x1 + x3, basis = ~ z)`, `forest(. - x3, ...)`, `vars = c("x1", "x3")`, a variable holding
    that vector, positions, and a call built with the names in it give identical draws; an unknown name;
    `list(forest(z3), forest(basis = ~ dose))`, with z3 a vector of 1s and 2s per row, is refused with
    the repeat's text (today it is a basis; without the refusal it would be read as columns 1 and 2 in
    silence).
2.8 Doors that take no term. `xbart`, `rbart_vi`, `bartBT`, `pdbart`, `pd2bart` and a direct
    [`dbartsData`](../../R/data.R) refuse a formula with a `forest()` term by name
    ([`formulaHasForestTerm`](../../R/formulaTerms.R)). Tests: the text at each (fails today: four other
    messages).
2.9 Respell and records. The 9 colon lines in 8 other test files, the three benchmark lines (the
    equivalence scenario among them) and their two comments become `forest(<x>, basis = ~ <w>)`. Help:
    man/bart.Rd's section "Formula Terms" and its `formula`, `test` and `subset` items, its example;
    man/forest.Rd's usage, `vars` item (with the first half of the warning), its Details paragraph on
    terms and its See also; man/dbartsForests.Rd's sentence. inst/NEWS.Rd: the multi-forest item says
    `forest(x1 + x2, basis = ~ z)` in place of the colon. Cites: `retired:` before each checked cite of a
    function this push removes, with the prose saying it is gone; the fragment
    ["expectSameForest"](../../inst/tinytest/test-formula-terms.R) is kept in the rewritten file or its
    cite retired.

### Push 3: the basis

3.1 Capture. [`forest`](../../R/model.R) reads `basis` unevaluated, as "The capture rule" says;
    [`forestBasisDeclarations`](../../R/model.R) carries code and values apart. Tests (a new file
    test-forest-capture.R): `basis = dose` with `dose` only in `data` at both doors (fails today: not
    found); with a second `dose` beside the call the data's is used, at both doors (fails today: the
    caller's); the matrix interface finds it where called; a held formula, a tilde in place, a value
    through `do.call`, `as.name()` through `do.call` in an `lapply`, a caller's variable named `forest`;
    `NULL` held and computed; a string, a single number and a two-sided formula refused with their texts;
    a column `b` of the data beside `b <- ~ dose` uses the column, and `do.call` with the formula uses
    dose.
3.2 The grammar walk (`parseBasisGrammar`), before anything is evaluated: `+` separates, parentheses
    group, `:` between terms is kept, `1`, `0 +` and `- 1` are read as in `lm`; everything in "Refused
    forms" marked [3] from the star to the reserved names is refused here. Tests (test-forest-basis-terms.R,
    rewritten): every refusal with its text; the 57 operator shapes the design's critic tried at the top
    of a basis are each refused or are the columns `lm` gives; dec-B272's rows three and four;
    `base::cbind(dose, age)` is a matrix term as in `lm`; a column named `normal`, `fixed` or `forest`
    is a column.
3.3 Building. A code basis is `stats::model.frame` on `~ 0 + <basis>` with the fit's `subset` and
    na.action, then `stats::model.matrix`; one factor, character or logical term is expanded by
    [`expandForestBasis`](../../R/model.R), not by a contrast, and a mix is refused. The stored columns
    are named as `model.matrix` names them; a value keeps its names when every column has one and they
    differ, and is otherwise unnamed; two columns alike are refused. Tests: 18 shapes against
    `names(coef(lm(y ~ 0 + <basis>)))`; the six changed texts, each against the tip's spelling of its new
    meaning written as it must be written now (`I(dose + age)` and `dose + age` differ; `dose - 1` is
    `dose`; `1 + dose` has columns `(Intercept)`, `dose`); a constant column alone refused.
3.4 Rows, at both doors. A level that the kept rows leave empty is dropped for code, whatever emptied it;
    a logical left with one value is refused; a missing value is refused, as today. The list door hands
    the fit's `subset` to the basis's model frame ([`dbarts`](../../R/dbarts.R),
    [`dbartsSpec`](../../R/spec.R)); [`alignForestBasisToSubset`](../../R/model.R) keeps the value's
    rule. Tests: dec-B272's row eight (fails today: one door refuses, and `scale()` differs between the
    doors); `scale(age)` under `subset` equals `(age - mean) / sd` of every row at both doors (fails
    today at the formula door); a subscripted value refused with the same-length text.
3.5 New rows. The stored record is R's `terms` object with its environment, at both doors, and
    [`replayForestBasis`](../../R/model.R) rebuilds through `model.frame(terms, newdata, xlev = )`;
    [`forestBasisPredictCall`](../../R/formulaTerms.R) goes. Tests: `scale()`, `poly()`, both, and
    `I(age - mean(age))` at new rows, each equal to the constants written out; a basis written in a list
    is rebuilt (fails today: stops); a constant beside the call is used and looked up again; a per-row
    vector found beside the call is refused at predict, pointing to `bases =`; a level the fit never saw
    is refused.
3.6 Labels, by "The label rule": the multi-forest block of [`resolveSamplerSpec`](../../R/spec.R) records
    them, [`packageBartResults`](../../R/bart.R) attaches them always, and
    [`resolveForestSpreads`](../../R/dbarts.R) checks names against them. Tests: every literal of the
    rule at both doors (fails today: `NULL`), the two refusals, the suffix, a held formula's label, the
    labels unchanged by a swap.
3.7 [`setForestBasis`](../../R/dbarts.R): a value, or a one-sided formula read as a basis is; columns by
    position; the recorded names stay; the recorded names in another order are refused; other names are
    ignored. Tests: the four cases (fails today: the names follow the value).
3.8 Respell and repair. test-forest-basis-terms.R is rewritten around steps 3.2 to 3.5; the one
    `~ cbind(a, b)` elsewhere becomes `a + b`. Expected to need a new expectation, from the prototype's
    run of the suite: about 11 pins of `data@bases` against an unnamed matrix, 2 of labels, 2 comparisons
    of a stored description that now holds an environment, 5 texts about a basis's rows under `subset`,
    one emptied level and one predict that now rebuilds. Run the suite first and repair what it shows;
    do not respell the 178 `basis = ~` lines.
3.9 Help and records. man/forest.Rd rewritten around "The grammar": the `basis` item with the warning,
    `sd` (one number for every column, in today's unit), the examples without tildes and with the first
    forest written out; man/bart.Rd's section; man/dbarts.Rd's `forests`;
    [`dbartsSampler$setForestBasis`](../../man/dbartsSampler-class.Rd) and the method's docstring;
    the vignette's line. docs/design/written-surface.md with its index row: the grammar, the capture
    rule, the label and naming rules, the six changed texts with their oracles, and what dec-B272 fixes;
    docs/design/public-surface.md and docs/architecture.md where they describe the old doors. TODO:
    `basis-formula-leftovers` closes; the entries of "Out of scope" are added.
3.10 Mutations (Verification): apply each, install, run the named test, record the failing count, revert,
    `touch` the file.

## Verification

Every push, against the slice's own library (`R CMD INSTALL -l <lib> .`, `R_LIBS=<lib>` on every call;
check `dbarts:::buildInfo()$mode` and that the install postdates the source):

- `cd tests/cpp && make && ./test_bartcore`: unchanged and passing (nothing under src/ moves).
- The full tinytest suite on the shipped build, counted per file: no failure, no file stopping, and at
  least the tip's 14274 assertions less those removed with the two rewritten files plus the new ones;
  the landing note gives the three figures.
- On a reference build (`--preclean --configure-args=--enable-reference-build`): the four
  `test-reproducibility-*.R` files pass unchanged, and the three compares are bitwise, every scenario
  reporting identical draws, counted per scenario with no `max |z|` line: 55 against
  `equivalence-1b7d730c.rds`, 15 against `bcf-equivalence-1b7d730c.rds`, 11 against
  `multinomial-equivalence-80b1c8d4.rds`. Nothing is re-recorded.
- The pair script (below), run on the tip's library and on the slice's.
- The consumers, each suite whole against a private install of the slice, none failing: bartCause on
  dbarts-1.0 with its edit (1406 expectations at its last landing), stan4bart on bartcore (582),
  treatSens on dbarts-1.0 (306), bairrtt on main (207), the last three unedited. Before push 1's edit,
  also run bartCause unedited once and record that it fails only where `update.amplitude` is passed.
- `Rscript -e 'lintr::lint_package()'`, `air format --check .`, `Rscript tools/check-rc-codoc.R .`,
  `Rscript tools/check-win-drift.R .`, `Rscript tools/check-doc-freshness.R .`,
  `Rscript benchmarks/R/mutation-battery.R verify-anchors`, each on its own exit status;
  `pkgdown::check_pkgdown(".")`; `R CMD check --as-cran` on a tarball from a clean copy.

Per push, in addition:

- Push 1: `bcf-exact.R`, `bcf-exact-weak.R`, `bcf-exact-restricted.R` and `bcf-latent-exact.R` in `quick`
  (their calls are respelled), unchanged; `sbc.R`'s respelled arm constructs and runs its smallest
  setting. Pairs 18 and 30, and the old side of every other pair on both builds: identical.
- Push 2: `tools:::.build_news_db_from_package_NEWS_Rd("inst/NEWS.Rd")` is not `NULL` and its entry count
  is the tip's. `composition-matrix.R`'s respelled row and the respelled surface script each run their
  one fit. Every pair in the spelling push 2 accepts: the new term or list, with the basis as the tip
  writes it. The pairs that differ only in the basis (03, 05 to 09, 12, 19 to 21, 25 to 27, 32 to 35) are
  then one text on both builds and are still run.
- Push 3: every gate `.github/workflows/exact-gates.yaml` lists, in `quick`, unchanged (a fit carries
  another record and the class is posterior-changing); every pair in its final spelling; a test per
  changed text (step 3.3, 3.1, 3.4); the design note.

The pair script. Each row is one model, fitted in the old spelling on the tip's library and in the new
on the slice's with the same seed, on one fixture (150 rows; x1, x2, x3, dose, age, a 0/1 z, a logical
zl, a character g, an offset column; a second data frame of 40 new rows): `bart` fits compare
`yhat.train`, `sigma` and `predict` at the new rows; sampler fits compare the train draws, sigma and
`getForestAmplitudes()`. `identical()` on every row but the last. Run for this plan with the design's
prototype in the slice's place: 34 identical, the last 4.25e-13 apart.

| | old, on the tip | new, on the slice |
|---|---|---|
| 01 | `y ~ x1 + x2 + z:forest(x1)` | `... + forest(x1, basis = z)` |
| 02 | `factor(z):forest(x1)` | `forest(x1, basis = factor(z))` |
| 03 | `forest(x1 + x2, basis = ~ cbind(dose, age))` | `basis = dose + age` |
| 04 | `y ~ x1 + x2 + forest(x1 + x2, basis = ~ dose)` | `y ~ forest(x1 + x2) + forest(x1 + x2, basis = dose)` |
| 05, 06 | `basis = ~ scale(age)`; `~ poly(dose, 2)`, with predict | the same, no tilde |
| 07 | `~ cbind(scale(age), poly(dose, 2))`, with predict | `scale(age) + poly(dose, 2)` |
| 08 | `~ dose + age`, the sum | `I(dose + age)` |
| 09 | `~ I(age - mean(age))`, with predict | the same, no tilde |
| 10 | `(dose + age):forest(x1)` | `forest(x1, basis = dose + age)` |
| 11, 12 | `g:forest(x1)`, a character column; `basis = ~ zl`, a logical | `basis = g`; `basis = zl` |
| 13 | `factor(z):forest(x1, sd = 2, power = 2.5, base = 0.4, n.trees = 12L)` | the first forest written out, `forest(x1, basis = factor(z), <the same>)` |
| 14 | three forests, a colon term and a `~ scale(age)` term, each with `sd` | three `forest()` terms |
| 15 | `dose:forest(x1) + factor(z):forest(x2)` under `subset` | the two as `basis =` |
| 16 | list: `forest(basis = ~ dose, sd = 1.5)` | `dose = forest(basis = dose, sd = 1.5)` |
| 17 | matrix interface: `forest(basis = ~ factor(z), vars = c("x1", "x3"), n.trees = 9L)` | `forest(x1 + x3, basis = factor(z), n.trees = 9L)` |
| 18 | `update.amplitude = FALSE` on both forests | `amplitude = fixed()` on both |
| 19, 20 | `basis = ~ dose` | a held formula; a numeric vector through `do.call` |
| 21 | a term with the tilde written | the same text |
| 22 | `forest(basis = ~ factor(z), vars = c("x1", "x3"))` | `forest(vars = c("x1", "x3"), basis = factor(z))` |
| 23, 24 | probit and logistic, a colon term | `forest()` terms, the first written out |
| 25, 26, 27 | `~ dose * age`; `~ dose - 1`; `~ 1 + dose`, each one column of arithmetic | `I(dose * age)`; `I(dose - 1)`; `I(1 + dose)` |
| 28 | weights and `offset(off)` beside a colon term | beside `forest()` terms |
| 29 | `list(forest(n.trees = 11L), forest(basis = ~ factor(z), sd = 1.2))` | the same text |
| 30 | a data object with bases and `call("forest", vars = <names>, ..., update.amplitude = FALSE)` | `amplitude = quote(fixed())` |
| 31 | `.` with removals beside a colon term | beside `forest()` |
| 32 | `~ cbind(1, dose)` | `1 + dose` |
| 33 | `~ dose`, the data's | `basis = dose` with another `dose` beside the formula |
| 34 | `~ cbind(dose, age)` in a list, then `$setForestBasis` and 25 more sweeps | `dose + age`, the same swap |
| 35 | under `subset`, `~ I((age - m) / s)` with every row's m and s written as numbers | `scale(age)`: equal to 1e-10, not identical |

Mutations, each expected to fail the named test and no gate before it:

- push 1, `fixed()` leaves the eighth number at drawn: step 1.1's "amplitudes stay at their starting
  values";
- push 1, the check coerces with `as.double` again: step 1.2's `"2"` and `TRUE`;
- push 1, the check drops its test of names: step 1.2's `c(a = 1)` at the three places;
- push 2, a second unnamed argument is matched to `basis`: step 2.2's `forest(x1, x2)`;
- push 2, the forest with no basis is left where it was written: step 2.4's swapped terms;
- push 2, a forest's own terms are added to the predictors but its restriction is dropped: step 2.3's
  "splits on x3 alone";
- push 2, an `offset()` among the plain terms is counted as a predictor term: step 2.1's six placements;
- push 3, the list door builds a code basis before `subset`: step 3.4's four pairs;
- push 3, a name is looked up where the call was made before `data`: step 3.1's second `dose`;
- push 3, the stored terms lose their `predvars`: step 3.5's `scale()` at new rows;
- push 3, `+` at the top of a basis is evaluated as arithmetic: step 3.3's `dose + age`;
- push 3, a star is read as `lm`'s crossing: step 3.2's `dose * age`;
- push 3, `student` leaves the reserved names: step 3.2's thirteen;
- push 3, a label is the name of the variable that holds a formula: step 3.6's held formula;
- push 3, `$setForestBasis` takes the value's names: step 3.7.

Not a hot-path change: nothing a sweep runs is touched.

## NEWS

No new item: `forest()` and the multi-forest formula are new in 1.0-0 and nothing released changes. One
existing item, "Multi-forest models", shows the colon form and is respelled in push 2.

## Out of scope, and where it goes

- Defaults by kind (dec-B274): a multiplied forest's default tree count and tree prior wherever it is
  written, and the refusal of a top-level count or tree prior when no forest is without a basis. Its own
  slice, after this one; it turns step 2.6's test around and rewrites one help sentence.
- The unit of `sd` (slice N, directly after this one); the kind by class at the data door and in the swap
  (K); the default per column, the reader's per-column entries, `extract`'s list, the printed block and
  `amplitude.prior.variance` (L); `updateBasisScale` (U); `tree.prior` on `forest()` and on the variance
  forest (P, V); `leaf.prior = normal(sd = )` on `forest()`, `$setLeafPrior`'s single form and selecting
  a forest by label (S). TODO `forest-prior-args` lists them.
- A size per column, in any of its three spellings: after the release (dec-B272).
- `dbarts()`'s own `n.trees`, a bare `forest()` for `forests`, and the both-given refusal at that door:
  the control-migration arc (dec-B241).
- To TODO as new entries: a "restarting interrupted promise evaluation" warning beside the error when a
  wrapper passes an unevaluable basis through a promise; a hazard fit of several forests skips the check
  that a basis covers the data under `subset`, its rows being expanded; a missing value in a basis is
  refused where `lm` drops the row, a departure from "as in `lm`" that the help names.

## Calls made in planning

- Three pushes, cut between the arguments, the forests of a formula and the basis. The alternative, one
  push, is about 4800 lines in one review. Cutting the grammar in two leaves a tip on which the first
  forest is written `forest(x1 + x2)` and a basis still `~ dose`; that tip is coherent (the tilde stays
  valid for good, and its help is the tip's), and about 15 help lines are touched twice.
- Six changed texts, not the design's seven. Measured at the tip, `basis = ~ dose:age` is refused today
  (a length error), so it becomes accepted and changes no model; and `~ dose * age` is refused afterwards
  under dec-B273. One is added that the design's list did not carry: `basis = dose` in a `forests` list
  beside a data column of that name, where today the caller's is used.
- A repeat in a selection of predictors is refused. With `vars` first, an old positional basis in a list,
  `forest(z3)`, is read as a selection; when its values are whole numbers no larger than the number of
  predictors it would be accepted as columns, a forest with no multiplier, in silence (run on the
  prototype). The alternative, refusing a selection longer than the predictors, says less clearly what is
  wrong. Cost: `vars = c("x1", "x1")`, column 1 today, is refused.
- The star's text names the ruling's two forms, `I(dose * age)` and `dose + age + I(dose * age)`. A colon
  between two basis terms stays accepted as the product column: it has no arithmetic reading, and the
  critique kept it. The alternative is refusing it too; it would then be an addition later.
- A value's column names: kept only when every column has one and they differ, otherwise none. The
  design said "its colnames or none"; a value such as `cbind(1 - z, z)`, which bartCause passes, has one
  name of two, and the names are fixed for good. The alternative, refusing a partly named value, breaks
  that call.
- The sd refusals carry no unit ("for every column of its basis") where the design's text says "that many
  units of the response": the unit is the latent scale until slice N, and a message must not promise
  another. Slice N may add it.
- The front door's `normal(sd = )` refuses from push 1 what a forest's `sd` refuses, as early checks
  in front of the ones it has, so that its 12 pinned refusals in 4 test files keep their texts and the
  check it shares with `k` is not touched. The design had one validator with one set of texts at every
  place; that rewrites the shared check and moves those pins. The other alternative, `forest(sd = )`
  alone until the long form arrives, leaves `normal(sd = TRUE)` accepted and two spellings of one
  statement judged differently.
- Doors that take no `forest()` term refuse it by name (step 2.8). It is not in the design. With the
  first forest written as `forest(x1 + x2)` in the help's lead example, a user will write it in `xbart`,
  and today four doors answer with four unrelated messages. The alternative is a TODO entry.
- The all-multiplied formula: accepted in written order and pinned, per dec-B274 and the coordinator's
  instruction, where the design's critic and this plan's brief had it refused.
- The pair script is run at landing and not tracked: its old side needs the tip's build. What stays in
  the suite is the same identity between spellings the slice still accepts (steps 2.4, 2.7, 3.1).
- Help between pushes: push 2 writes the multiplier as `basis = ~ z` everywhere, push 3 drops the tilde.
  The alternative, landing the help once with push 3, leaves a tip whose help shows the colon.
- Opus for the grammar's R code. The design said sonnet with an opus review; the reasons are on the
  `agent:` line.
