# leaf-prior-k-or-sd: name the leaf prior by k or by sd

Status: LANDED 2026-09-29 (b53a637f, a8b36094, 7a52285f, 645c9356, e00f7acb); see Landing. The four maintainer questions are ruled (see Open maintainer questions), and the
writer was restated by a later ruling (see Landing).

agent: sonnet (R surface, bridge strings, manual; no engine change)
rng: neutral (every spelling that survives draws as today; each removed spelling's replacement reaches
the engine with identical inputs at the defaults, and the equivalence baselines use none of them)
budget: ~1000 changed lines in dbarts, ~490 of them test (R ~280, bridge ~10, man/NEWS/vignettes ~175,
standing docs ~55, benchmarks ~20, tests ~490), plus the answer-dependent increments named in each open
question (at most ~+170). The Step 0 script (~45) stays out of the tree. stan4bart 0, bartCause 0.

## Goal

A user names the leaf prior one of two ways, never both:

- `k`: shrinkage relative to the anchor the data fixes. A number, a string such as `"chi(1.5)"` (kept
  for 0.9-x), or a hyperprior on k, `chi(degreesOfFreedom, scale)`, as today.
- `sd`: the prior standard deviation of the forest total, on the scale the family's forest fits. A
  number, or a hyperprior on the sd itself, `invchi(nu, c)`. It takes no string form.

The anchor k is relative to, and the units sd is stated in, per family:

| family | sd is stated in | anchor at k = 1 |
|---|---|---|
| gaussian, Student-t | response units | half the training response range |
| aft | log survival time | half the observed log-time range |
| probit, weighted binary, ordinal | probit latent (ordinal: relative to the pinned first cutpoint) | 3 |
| logistic, nbinom | log-odds latent | pi * sqrt(3) |
| hazard (probit or logistic link) | that link's latent scale | 3 or pi * sqrt(3) |
| multinomial, hurdle, multi-forest | refused (see Calls made 10 and 11) | - |

`scale` leaves the user-facing surface: the `scale` argument of `normal()`, `linear()` and `gp()`,
`$setLeafPrior(prior.scale = )`, and, unless open question (b) keeps an anchor readout, the
`prior.scale` column of `$getLeafPrior()`. Every prior expressible today stays expressible, with one
exception stated in the mapping (hurdle). Maintainer ruling on dec-A105. `invchi` is a placeholder
for the name in open question (d); the implementer substitutes it everywhere.

## Context

- The design being replaced: [nameable-calibration.md](nameable-calibration.md), whose `prior.scale` is
  the forest total's prior sd at k = 1. The rename to leaf vocabulary: [leaf-vocabulary.md](leaf-vocabulary.md).
- Only the ratio of anchor to k enters a draw law. Under k ~ s * chi_nu (the engine's
  [`ChiKHyperprior`](../../src/bartcore/model.hpp) draws k^2 from a gamma with shape (M + nu) / 2 and a
  rate term 1 / (2 s^2)), the spread is sd = (anchor / s) / chi_nu: a scaled inverse chi on the sd
  scale, or equivalently sd^2 ~ scaled inverse chi-square with nu degrees of freedom. A named anchor and
  the hyperprior's own scale are not separately identified; only anchor / s is. That is why `scale`
  plus a k hyperprior collapses to one sd hyperprior. At s = Inf the kernel is the improper
  sd^-(nu + 1), which is the sd hyperprior at c = 0.
- R today: [`normal`](../../R/model.R), [`linear`](../../R/model.R) and [`gp`](../../R/model.R) take
  `k`, `sd` and `scale`; retired: [`resolveNamedScaleArgs`](../../R/model.R) kept at most one of `sd`
  and `scale`, and retired: [`resolvePriorScale`](../../R/model.R) turned `sd` into `sd * k` and refused
  it under a hyperprior (both removed by this plan); [`resolveLeafHyperprior`](../../R/model.R) turns `k` into a
  [`dbartsFixedHyperprior`](../../R/A_class.R) or passes a
  [`dbartsChiHyperprior`](../../R/A_class.R) through. [`resolveSamplerSpec`](../../R/spec.R) and
  [`xbart`](../../R/xbart.R) write the result into the [`dbartsModel`](../../R/A_class.R) slot
  `prior.scale` (the anchor in response units, NA when unnamed).
- [`bart2Hurdle`](../../R/bart.R) passes the one `leaf.prior` to both halves: a probit zero part, whose
  default k is chi(1.5, 2), and a gaussian positive part on log y. The [`hazard`](../../R/family.R)
  family is one binary model on the person-period expansion, so a named sd there has one meaning.
- Bridge: [`parseModel`](../../src/R_interface_bartcore.cpp) reads `prior.scale` into
  [`ParsedModel`](../../src/R_interface_bartcore.cpp) and, for a chi hyperprior, reads df and scale but
  never k, so a drawn k always starts at `ParsedModel`'s default of 2. The engine converts the anchor at
  [`Chain::resolvedNodeScale`](../../src/bartcore/chain.hpp). The engine and the bridge need no change:
  every form below reaches them as an anchor (or NA) plus a fixed k or a chi hyperprior.
- Sampler: [`getLeafPrior`](../../R/dbarts.R) returns the bridge's twelve columns (the first is
  `prior.scale`); [`setLeafPrior`](../../R/dbarts.R) takes `prior.scale` or `prior.sd` and writes an
  anchor through the bridge. The hyperprior is not saved state: it comes from the model at creation and
  at `setModel`, so the R5 `model` field always holds the hyperprior in force. The prior constructors
  are reachable only by vocabulary evaluation or through `dbartsPriors`, so a sampler method taking an
  argument evaluated normally cannot take an `invchi(...)` call; the writer takes numbers.
- The flat C header carries no leaf-prior entry (the dbarts.h freeze removed it), so this is not an ABI
  event.
- Released surface: CRAN and main are 0.9-34, whose `normal(k = 2.0)` has no `scale` and no `sd`, and
  no sampler reader or writer. Everything this plan removes is branch-only, so it is removed, not
  tombstoned (dec-B128). No released user was found.

## The mapping

Translation, in R, from the new forms to what the bridge receives (`model@prior.scale`, the hyperprior
slot). The reference k is 2 in both named sd rows.

| form | `prior.scale` | leaf hyperprior |
|---|---|---|
| `k` absent | NA | family default, as today |
| `k = k0` | NA | fixed k0 |
| `k = chi(nu, s)` | NA | chi(nu, s) |
| `sd = x` | 2x | fixed 2 |
| `sd = invchi(nu, c)`, c > 0 | 2c | chi(nu, 2) |
| `sd = invchi(nu, 0)` | NA | chi(nu, Inf) |

Why 2 and not the unit scale the ruling sketched: the bridge starts every drawn k at 2. With chi(nu, 2)
the chain starts at spread c, the named value, and the binary default and every old spelling at the
defaults reach the engine with identical inputs. With chi(nu, 1) the prior is the same but the chain
starts at c / 2 and matches no old spelling bitwise. Doubling is exact in floating point, so nothing
rounds. Measured on the tip (55808b01, reinstalled 2026-09-29): probit default against the scale-2
form, draws identical; against the unit-scale form, and against chi(1.5, 4) with the same ratio,
different.

Every combination the tip accepts, mapped. "Identical inputs" means the same `prior.scale` and
hyperprior bits; since the engine does not change, identical inputs give identical draws on any build.

| old spelling (tip) | new spelling | old inputs -> new inputs | draws |
|---|---|---|---|
| `scale = P`, k fixed k0 | `sd = P / k0` | (P, k0) -> (2P/k0, 2) | identical inputs at k0 = 2 (the continuous default); bitwise for a power-of-two k0 (measured, k0 = 4 and k0 = 1); otherwise the same prior, draws differ by one rounding of P / k0 (measured 4.7e-15 at k0 = 3) |
| `scale = P`, binary default k | `sd = invchi(1.5, P / 2)` | (P, chi(1.5, 2)) -> same | identical inputs |
| `scale = P`, `k = chi(nu, s)`, s finite | `sd = invchi(nu, P / s)` | (P, chi(nu, s)) -> (2P/s, chi(nu, 2)) | identical inputs at s = 2; otherwise the same prior, but the chain starts at spread P / s rather than P / 2 |
| `scale = P`, `k = chi(nu, Inf)` | `sd = invchi(nu, 0)` or `k = chi(nu, Inf)` | (P, chi(nu, Inf)) -> (NA, chi(nu, Inf)) | the same prior: the spread's prior is sd^-(nu + 1) whatever P is; P set only the starting spread (measured: two P, same seed, draws differ) |
| `sd = x`, k fixed k0 | `sd = x` | (x k0, k0) -> (2x, 2) | identical inputs at k0 = 2; otherwise as the first row |
| unnamed default, any family | unchanged, or `sd = invchi(1.5, anchor / 2)` on binary and `sd = range / 4` on gaussian | same | identical inputs (measured for probit, logistic and gaussian) |
| hurdle with `scale = P` (both halves) | refused | - | LOST: one sd cannot mean probit units on one half and log-y units on the other. The same posterior is two ordinary fits, a probit on 1{y > 0} and a gaussian on log y over the positive rows, each with its own sd, since the hurdle's halves are independent models |
| `$setLeafPrior(prior.scale = P)`, k fixed k0 | `$setLeafPrior(prior.sd = P / k0)` | writes P -> writes (P / k0) k0 | identical write for a power-of-two k0 |
| `$setLeafPrior(prior.scale = P)`, k ~ chi(nu, s), s finite | `$setLeafPrior(prior.sd.scale = P / s)` | writes P -> writes (P / s) s | identical write for a power-of-two s, the default 2 included |
| `$setLeafPrior(prior.scale = P)`, k ~ chi(nu, Inf) | `$setLeafPrior(prior.sd.scale = 0)`, a no-op | - | nothing lost: the old write changed only the next sweep's leaf draws; the spread's law ignores P |
| `xbart(k = grid, leaf.prior = normal(scale = P))` | per open question (c) | per cell as the first row | each cell's prior is kept; whether one call can still sweep them is question (c) |

Two things that are not priors change: the starting value of a drawn k when an old spelling paired a
named scale with a k-hyperprior scale other than 2, and (under answer c1) the one-call absolute-sd
sweep in `xbart`.

## Open maintainer questions

### (a) Does a named sd stay absolute when the response is re-anchored?

Today an absolute sd does not stay absolute. Measured on the tip: `dbarts(x, y, leaf.prior =
normal(k = 2, scale = 2))`, the new `normal(sd = 1)`, reads `prior.sd` 1; after `setResponse(10 * y,
updateScale = TRUE)` it reads 10, while the model's recorded intent stays 2. The engine holds the leaf
scale in internal units, and the re-anchoring channels (`setResponse` and `setOffset` at `updateScale =
TRUE`, and `setData`) move the transform under it.

- a1, document it (0 lines of code, ~10 of manual): a named sd is stated against the transform in
  force at creation and at `setModel`; a re-anchoring channel scales it; `setLeafPrior` restates it.
  Cost: the absolute spelling is absolute only between re-anchors, and a reader must know which
  channels re-anchor.
- a2, re-issue it (~+40 R, ~+40 test): after each re-anchoring channel, when the model names an sd
  (`model@prior.scale` finite), the R5 method writes the named anchor back through the existing bridge
  writer on every forest it applies to. `setLeafPrior` then also records its write in the R5 model, so
  the re-issue restores the latest intent rather than the creation one. Cost: a re-anchor on a
  named-sd sampler draws differently than today (the named path only; no default moves); reverses
  nameable-calibration's authority rule that the engine alone holds what is in force.

RULED a2 (maintainer, 2026-09-29, after asking whether it is right by design and what users expect, and shown that a fixed residual sd and the variance forest's prior already keep their meaning in response units across a re-anchor: "OK, proceed using it."). Step 4, item 3 implements a2.

### (b) What does k mean when sd is named?

The engine's k is relative to whatever anchor it holds. Under a named numeric sd that is 2x, so
`getLeafPrior`'s `k` reads 2; under an sd hyperprior a fit's `k` draws are 2c / sd; after a
`setLeafPrior` write, k is relative to the written anchor. Under the k spelling, k is relative to the
data's anchor. So one column means three things.

- b1, report the engine's k (0 lines): the manual states the three readings. Cost: `k` and `fit$k`
  are hard to read on any named-sd fit, and nothing on the reader shows the in-force anchor once
  `prior.scale` goes.
- b2, report a data-relative k (~+50 R, ~+40 test, ~+10 manual): R derives k = defaultLeafScale(family)
  * response.scale / prior.sd on the reader and divides `fit$k`'s draws the same way (in
  [`packageBartResults`](../../R/bart.R) and the extract path) whenever the model names an sd; an
  unnamed model reports the engine's k unchanged. The reader keeps one in-force anchor column, the data
  anchor, named `anchor` rather than `prior.scale` since it is no longer the named quantity. Cost: on
  named fits `k` is a derived value, not bitwise the engine's; one more column (14); multi-forest and
  multinomial are untouched (they refuse a named sd).

RULED (maintainer, 2026-09-29): "Yes, report in a way that matches what a user specified / will want to be able to use." Read as b2 extended: a k-named fit reports k (draws included) as today; an sd-named fit reports the spread - prior.sd on the reader and, under an sd hyperprior, the fit's draws of the spread instead of k; the reader keeps both, with k shown data-relative and an `anchor` column; every reported value writes back through setLeafPrior in the same terms. Step 4, item 1 implements this; the implementer adds the sd-draws channel to packageBartResults and extract.

### (c) xbart with a named sd and a k grid

Today `xbart(k = grid, leaf.prior = normal(scale = P))` sweeps the absolute sd P / k_i, held fixed
across folds. For the binary families the anchor is a constant of the latent scale, so an sd grid IS a
k grid (sd_i = anchor / k_i, identical inputs at power-of-two ratios, and invchi(nu, c_i) is chi(nu,
anchor / c_i)); nothing is lost there. For gaussian and aft the anchor is each fold's training range, so
a k grid is relative per fold and an absolute sd grid is a different sweep.

- c1, refuse `k` beside a named sd (~40 R): a named sd rides every cell as one k-axis cell; the manual
  gives the binary k-grid equivalent and, for gaussian, one `xbart` call per sd with a shared seed.
  Cost: the one-call absolute sweep with warm starts across sd values, which only the branch ever had.
- c2, accept a list of leaf priors as the grid axis (~+60 R over c1, ~+40 test): `leaf.prior =
  list(normal(sd = 0.5), normal(sd = 1))` sweeps them, labelled by each prior's sd form, and is refused
  beside a `k` grid. Cost: a second grid-axis vocabulary on `xbart`.

RULED c3 (maintainer, 2026-09-29: "(c) seems right."): xbart gains an `sd` grid argument beside `k`, exclusive with it, sweeping absolute spreads held fixed across folds with the same warm starts and unit parallelism the k grid has; a named sd in leaf.prior beside a k or sd grid is refused, naming the grid argument; the result's grid dimension is labelled sd. Step 3 implements c3 (about the cost of c2).

### (d) The sd hyperprior's name

The constructor states sd = c / chi_nu, c on the sd's scale, nu degrees of freedom. The name decides
which scale the user thinks on (sd or variance) and whether the parameters carry a sqrt(nu) or a square.

Why only this family: it is the one the engine can honor. k^2 is conjugate to the normal leaves, so
the engine draws it exactly from a gamma; any law on the sd that is a chi on k (scaled inverse chi on
the sd) comes free, and a half-t or half-Cauchy on the sd would need a new, non-conjugate engine
update. The manual says so.

What others call this family (per each package's documentation; re-check at implementation):

- Base R and the recommended packages: nothing. No inverse chi-square or inverse gamma density.
- Stan: `scaled_inv_chi_square(nu, s)` on the variance, s on the sd scale, variance = nu s^2 / chi^2_nu;
  also `inv_chi_square(nu)` and `inv_gamma(alpha, beta)`. brms uses Stan's names.
- geoR and LaplacesDemon: `dinvchisq(x, df, scale = 1/df)`, scale on the variance scale (s^2).
  extraDistr: `dinvchisq(x, nu, tau)`, tau = s^2. extraDistr, LaplacesDemon, MCMCpack and the invgamma
  package: `dinvgamma` with shape and scale (or rate).
- No R package found ships a density for the inverse chi on the sd scale.
- dbarts's own vocabulary: `chisq(df, quant)` is already a scaled inverse chi-square on sigma^2, named
  after the chi-square it inverts. And `chi`'s first argument is `degreesOfFreedom`, not `df`, and
  `chi(df = )` does not partial-match, so a new constructor can mirror `chi` or `chisq` in its first
  argument's name, not both.

Options, ranked:

1. `invchi(degreesOfFreedom, scale)` or `invchi(df, scale)` on the sd scale. The exact mirror of `chi`:
   k ~ chi(nu, s) is the same prior as sd = invchi(nu, anchor / s). Argument and variate share units,
   the parameters carry no sqrt(nu), and the binary default is exactly `invchi(1.5, 1.5)` on probit.
   Cost: an unfamiliar name, easily read as inverse chi-square; the word scale stays as a distribution
   parameter (as in `chi`, `dcauchy`, `dgamma`); the first argument's name must be chosen against
   `chi` or `chisq`.
2. `invchisq(df, s)` on the variance, Stan's convention (s on the sd scale). The best-known family
   name. Cost: a variance law inside an `sd =` argument; c = s * sqrt(nu), so the binary default is not
   exactly representable (sqrt(1.5)); geoR, LaplacesDemon and extraDistr use the same name with s^2;
   and it gives one family two names in one vocabulary, since `chisq` already is that law on sigma^2.
3. `invgamma(shape, scale)` on the variance. Widely known in R; shape = nu / 2, scale = c^2 / 2. Cost:
   halves and squares, a variance law in an `sd =` argument, and scale-versus-rate conventions differ
   across packages.
4. Reuse `chisq(df, quant)`. `quant` calibrates against `sigest`, a data-fixed value. Carried over,
   P(sd < data anchor) = quant is a RELATIVE statement, a respelling of the k form (the binary default
   `chi(1.5, 2)` is quant = 0.78), not the absolute sd form the ruling asks for. An absolute version
   needs a third argument or a second meaning for `chisq`.

Recommendation: option 1. Evidence that would change it: the maintainer preferring a recognized family
name over unit consistency (then option 2, Stan convention).

The ruled name also fixes the parameter names (placeholders `nu`, default 1.5 as `chi`'s, and `c`, no
default, since an absolute scale has no data-free default) and the vocabulary entry. The reader and
writer columns `prior.sd.df` and `prior.sd.scale` keep those names whatever the constructor's
parameters are called. The S4 class is `dbartsSdHyperprior` whatever the name.

RULED (maintainer, 2026-09-29): `invchi(df, scale)` on the sd scale ("yes, I think `invchi` is right"), with `df` as its first argument, as base R's dchisq and dt spell it. And, asked "Should we rename `degreesOfFreedom` to `df` in `chi`?", shown that the name shipped in 0.9-34: "Use (a)." So `chi(df, scale)` too, with `degreesOfFreedom` a tombstone (warns once per session, uses the value, removed in 1.1-0; registry entry, man, NEWS, test); positional and string calls are unchanged, and the internal class's slot keeps its name. The placeholder invchi is `invchi` throughout.


### dec-A106, ruled with the plan

One `sd` for every leaf model (maintainer, 2026-09-29, after an advisory panel of four independent reviewers - an applied user, a base-R package developer, a Bayesian methodologist and a newer user - all chose it: "Sure, sounds good."). The panel's refinements are part of the ruling and of Steps 4 and 6:
1. `sd` is defined as the standard deviation of the normal prior on the leaf model's own parameter (the leaf value, the coefficients, the GP amplitude); the bound it implies on the prior spread of f(x) is documented as a consequence, never as the definition.
2. Every constructor's entry says `sd` (and `k`) is for the whole forest, the sum of trees, not one tree.
3. For linear leaves the first line says "per standardized covariate".
4. Monotone leaves get an exact statement of which parameter `sd` scales, not only a bound direction (the implementer derives it from the engine and states it).
5. `getLeafPrior` gains a column naming what `prior.sd` refers to ("leaf value", "coefficient", "amplitude").
A prior-predictive helper returning the prior spread of f at given x is filed in TODO, not built here.

## Constraints

- RNG neutral. The equivalence trio and the tinytest snapshots use no removed spelling; the three exact
  gates that do (aft-exact, t-exact, logistic-reference) all run k = 2 and move to identical inputs.
  Answer a2 changes draws on the new named path only, after a re-anchor, where no baseline reaches.
- Engine, bridge logic and inst/include/dbarts/dbarts.h are unchanged. The bridge changes three refusal
  strings only.
- Internal names keep scale (Calls made 7): the `dbartsModel` slot `prior.scale`,
  `ParsedModel::priorScale`, the bridge's `bartcore_setLeafPrior` argument, and the engine's
  [`Chain::setForestPriorScale`](../../src/bartcore/chain.hpp) and
  [`ForestCalibration`](../../src/bartcore/chain.hpp).
- Out of scope: `forest(sd = )` and the reader's `amplitude.prior.scale` (the multi-forest amplitude
  channel, a different quantity), `chi()`'s own `scale` parameter.
- Frozen records stay as written: docs/decisions.md (the orchestrator's), landing notes, LANDED plans
  including nameable-calibration.md, and NEWS sections before 1.0-0. Landed design docs change only
  where they describe the present surface (Step 8).

## Steps

0. Before any edit, `R CMD INSTALL --preclean` the tip into a private library, then, in a scratch
   directory outside the tree (~45 lines, not committed), run a script that for every row of the
   mapping table builds the old spelling under one seed, two chains, gaussian and probit data, and saves
   (a) `sampler$model@prior.scale` and the hyperprior slots as `sprintf("%a")` strings and (b) the
   train and k draws. Since the engine does not change, (a) alone proves the draws; (b) guards the R
   path (a slot the bridge reads that the translation forgot). Run it on the reference build
   (`--enable-reference-build`) so (b) is comparable after the change. The (a) strings become the
   literals the Step 5 test pins.
1. The constructor and classes (R/A_class.R, R/model.R; ~70).
   - `dbartsSdHyperprior`: slots nu and c, validity nu > 0 and finite, c >= 0 and finite. Not a
     subclass of `dbartsLeafHyperprior`, which means a law on k, so `k = invchi()` is refused.
   - `invchi(nu = 1.5, c)`: `c` missing is an error. Added to [`dbartsPriors`](../../R/model.R).
   - `dbartsNormalPrior`, `dbartsLinearPrior`, `dbartsGPPrior`: drop the `prior.scale` slot; `prior.sd`
     becomes class "ANY", NULL for unnamed, a positive number, or a `dbartsSdHyperprior` (mirroring the
     `k` slot).
2. The translation and its refusals (R/model.R, R/spec.R, R/bart.R, bridge strings; ~90).
   - `normal(k = NULL, sd = NULL)`, `linear(columns, k = NULL, sd = NULL)`, `gp(columns, k = NULL,
     lengthscale = NULL, max.leaf.size = 256L, sd = NULL)`. Both `k` and `sd` given: "give either 'k'
     (relative to the data's anchor) or 'sd' (on the family's scale) to a leaf prior, not both". A chi
     object as `sd`: refused, naming `k = chi()` and `sd = invchi()`. A string as `sd`: refused ("'sd'
     must be a number or invchi(); unlike 'k' it takes no string form"). `scale =` is an unused
     argument, R's own error.
   - Replace `resolveNamedScaleArgs` and `resolvePriorScale` with one resolver returning the anchor and
     the hyperprior per the translation table; [`resolveLeafHyperprior`](../../R/model.R) keeps the k
     branch. Under `monotone`, an sd hyperprior is refused as a k hyperprior is.
   - [`parsePriors`](../../R/model.R) and [`resolveSamplerSpec`](../../R/spec.R) call the resolver.
     The multinomial refusal list labels the offender "a named leaf-prior 'sd'". A multi-forest fit
     (`forests =`, bases, a causal forest) refuses any named sd before translation, so the message
     names sd rather than "a 'k' hyperprior", and says: "the leaf prior's 'sd' is not a forest's
     'sd': each forest's scale is set by forest(sd = ), which states that forest's share of the
     combined location's prior (see ?forest)".
   - [`bart2Hurdle`](../../R/bart.R) evaluates `leaf.prior` in the prior vocabulary before the split
     and refuses a named sd, number or hyperprior: "a hurdle fit passes one leaf prior to a probit zero
     part and a log-scale positive part, so one 'sd' cannot state both; name 'k', or fit the two parts
     separately". A named sd reaching the zero part would also have turned its chi(1.5, 2) default
     into a fixed k. `hazard` needs no refusal (one binary model).
   - Bridge: the two "a named 'prior.scale'" offender strings and `bartcore_setLeafPrior`'s
     "'prior.scale' must be ..." become "a named leaf-prior sd" and "the leaf-prior anchor must be
     ..."; they backstop callers that skip R.
3. [`xbart`](../../R/xbart.R) (~40 under c1; ~+60 under c2).
   - c1: a leaf prior naming `sd` together with the `k` argument is refused: "'k' and the leaf prior's
     'sd' state one spread two ways; drop 'k'. On a binary family the sd grid is the k grid anchor /
     sd; otherwise run one xbart per sd value".
   - Either answer: a named sd with no `k` is a one-cell k axis holding the translated hyperprior and
     rides every n.trees, power and base cell. [`kGridLabel`](../../R/xbart.R) labels that cell by the
     sd form (the number, or `invchi(nu, c)`) rather than by the translated chi. The vocabulary
     `xbart` evaluates `leaf.prior` in gains `invchi`.
   - c2 adds: `leaf.prior` may be a list of leaf priors of one leaf model, each a grid cell on the
     k axis, refused beside a `k` grid.
4. The sampler (R/dbarts.R, R/bart.R; ~80 under b1 and a1; ~+50 under b2; ~+40 under a2).
   1. `getLeafPrior`: the bridge's first column is dropped in R and two columns follow `prior.sd`:
      `prior.sd.df` and `prior.sd.scale`, the sd law in force when this forest's k is drawn (nu, and
      the anchor divided by the in-force chi scale; 0 under an infinite chi scale), NaN when k is fixed.
      The hyperprior is read from the R5 `model` field and gated by the engine's per-forest
      `k.has.hyperprior`, so a BCF or multinomial forest (k pinned) reads NaN. Under b1, 13 columns
      and the stacked read is n.chains x 13 x n.forests. Under b2, an `anchor` column follows
      `k.has.hyperprior` (14 columns), and on a named-sd model the `k` column and `fit$k` are the
      data-relative k of open question (b).
   2. `setLeafPrior(prior.sd, prior.sd.scale, prior.mean, forest = 1L, updateState = NULL)`, exactly
      one of the first two, both numbers, matching the reader's columns. `prior.sd` under a fixed k
      writes `prior.sd * k` as today (the diverged-k refusal stays); under a drawn k it is refused,
      naming `prior.sd.scale` and `setModel`. `prior.sd.scale` under a drawn k chi(nu, s) writes
      `prior.sd.scale * s` when s is finite and the value positive; 0 under s = Inf is accepted and
      writes nothing (the round trip of a read); 0 under a finite s and a positive value under s = Inf
      are refused, naming `setModel` (they change the hyperprior itself); under a fixed k it is
      refused, naming `prior.sd`. The df cannot be written mid-chain; `setModel` changes it.
   3. Under a2 only: `setResponse` and `setOffset` at `updateScale = TRUE`, and `setData`, re-issue
      the named anchor per open question (a); `setLeafPrior` records its write in the R5 model.
   4. Docstrings of both, and the cross-references in the other per-forest methods where they describe
      the columns.
5. Tests (~490; see Tests).
6. Manual, NEWS, vignettes (~175).
   - man/dbartsPriors.Rd: `normal`, `linear`, `gp` usage and text. `k` or `sd`, never both; `sd` a
     number or `invchi`, no string form; the family table of units and anchors. The `sd` argument's own
     entry states what sd bounds per leaf model (dec-A106): constant, exactly the prior sd of f(x);
     linear, sd is the forest total's prior sd of each coefficient on the standardized leaf covariates,
     so it is a LOWER bound on sd(f(x)), attained at the covariate origin; gp, an UPPER bound, attained
     at leaf members and decaying away from the leaf's data; monotone, a LOWER bound in the interior.
     A new `invchi` item: the k equivalence (`chi(nu, s)` is `invchi(nu, anchor / s)`, and
     `chi(nu, Inf)` is `invchi(nu, 0)`), the binary defaults in sd form (probit `invchi(1.5, 1.5)`,
     logistic `invchi(1.5, pi * sqrt(3) / 2)`), the Stan and inverse gamma equivalents, and why no
     other family is offered (conjugacy). `chi`'s item gains the sd equivalent.
   - man/dbartsSampler-class.Rd: usage, the `prior.scale` and `prior.sd` argument items replaced by
     `prior.sd` and `prior.sd.scale`, the `getLeafPrior` value paragraph (column count, the `k`
     column's meaning per open question (b)), and the decomposition identity restated on `prior.sd`
     (the calibration map pins k at 1, so it is exact).
   - man/dbarts.Rd "Naming the leaf calibration" (and the re-anchoring sentence per open question (a))
     and the multinomial refusal list; man/xbart.Rd's named-calibration paragraph per open question
     (c); man/bart.Rd's hurdle paragraph (the refusal and the two-fit recipe);
     man/dbarts-embedding.Rd's table row and "Whose prior is in force"; man/forest.Rd's two
     `prior.scale` reads and a sentence that `forest(sd = )` is not the leaf prior's `sd`.
   - inst/NEWS.Rd 1.0-0: rewrite the calibration item and its reader and writer sentences in the new
     spelling. No item says `scale` was removed, since it never reached main.
   - vignettes/dbarts-as-a-component.Rmd and gibbs_sampler_mixture_model.Rmd: `normal(k = 2, scale =
     2)` becomes `normal(sd = 1)`, `setLeafPrior(prior.scale = 2)` becomes `prior.sd = 1`, and the
     budget split reads and writes `prior.sd`. Each is the identical write at the k those samplers run.
7. Benchmarks (~20): aft-exact, t-exact and logistic-reference `normal(k, scale = X)` become
   `normal(sd = X / k)`; backfit-exact reads `prior.sd`; composition-matrix's `setLeafPrior` call is
   restated at its sampler's k; sbc.R's comment.
8. Standing docs (~55): docs/design/prior-defaults.md's "prior.scale (naming the calibration)"
   section becomes "Naming the leaf prior: k or sd" (the identification argument, the translation, the
   start-value rule, and the answers to (a) and (b)); design/feature-matrix.md's two `prior.scale`
   cells; design/INDEX.md's nameable-calibration row; design/nameable-calibration.md's Status line
   gains "vocabulary superseded by docs/plans/leaf-prior-k-or-sd.md", its body left as the record;
   design/multiplier-combiner.md's reader passages (the induced-index formula and the decomposition
   identity read `$getLeafPrior(f)[, "prior.sd"]`, the method names follow the leaf-vocabulary rename);
   design/multinomial-mutation-arc.md's refusal list ("named prior.scale" becomes "a named leaf-prior
   sd"). docs/plans/INDEX.md's row.
9. After the change, rerun the Step 0 script with the new spellings on the same reference build:
   every identical-inputs row matches bitwise, train and k draws; the other rows differ as the table
   says and nowhere else. The counts go in the Landing note.
10. Sister packages: no edit. stan4bart uses only the k form (`chi(1.25, Inf)` through `bart_args`),
    bartCause reads `response.scale` and `response.shift` by name and passes `k = "chi(1, Inf)"`.
    Run their suites against the new dbarts.
11. Records: Landing note here, the design Status line, the INDEX row.

## Tests

- New inst/tinytest/test-leaf-prior-k-or-sd.R:
  - Translation pins: each new spelling's `model@prior.scale` and hyperprior are `identical()` to the
    Step 0 literals of its old spelling (host-independent; pure R arithmetic). `invchi(nu, 0)` gives
    NA and chi(nu, Inf).
  - Live bitwise oracles, two chains, every chain compared: probit and logistic defaults against
    `sd = invchi(1.5, anchor / 2)`; gaussian `normal()` against `sd = diff(range(y)) / 4`; gaussian
    `normal(k = 4)` against `sd = diff(range(y)) / 8`; `k = chi(1.25, Inf)` against `sd =
    invchi(1.25, 0)`; one linear, one gp, one monotone and one hazard row with a numeric sd against its k-form equivalent. Each
    compares train and k draws with `identical()`.
  - Refusals: k with sd on each constructor; `scale =` unused on each constructor and on
    `setLeafPrior`; chi as sd, `invchi` as k, a string as sd; `invchi` validity (nu non-positive,
    c negative, either non-finite or NA, c missing); sd hyperprior under monotone; named sd, number and
    hyperprior, on multinomial, `forests =` (the message names `forest(sd = )`) and a causal forest;
    on a hurdle `bart`, number and hyperprior, checked before either half is built; `xbart` with k and
    a named sd (c1) or a list grid beside k (c2).
- inst/tinytest/test-calibration-creation.R, -midchain.R, -prior-draws.R: every `scale =` and
  `prior.scale` rewritten in the new spelling; their oracles keep their tolerances. The pinned reader
  dims in test-calibration-midchain.R (`c(2L, 12L)` twice, `c(2L, 12L, 2L)`, `c(1L, 12L)`) move to 13
  (14 under b2). Mid-chain adds: get then set is bitwise inert on a default binary sampler
  (`setLeafPrior(prior.sd.scale = )` of the read value skips the write) and on a `chi(1.25, Inf)`
  sampler (0 writes nothing); every `setLeafPrior` refusal of Step 4; the reader's column set and
  order, NaN under fixed k, `prior.sd.scale` equal to c after `normal(sd = invchi(nu, c))` and to 3 / s
  under probit `k = chi(nu, s)`. Under a2: `setResponse(10 * y, updateScale = TRUE)` on `normal(sd =
  1)` leaves `prior.sd` at 1, and after `setLeafPrior(prior.sd = 2)` at 2. Under b2: `k` on a named-sd
  fit equals defaultLeafScale * response.scale / prior.sd, and an unnamed fit's `k` is bitwise the
  engine's.
- test-embedding-recipes.R, test-augmentation.R, test-bcf-family.R, test-forest-basis-r5.R,
  test-multinomial-r5-surface.R, test-bartcore.R, test-na-as-none.R, test-argument-surface.R,
  test-tombstones.R, test-model-errors.R: the spellings and pinned messages that name `prior.scale`
  or `scale =`.
- tests/cpp: nothing.

## Verification

Against the slice's own library (`R_LIBS=$LIB` on every R call), each gate on its own exit status:

```sh
R CMD INSTALL --preclean -l $LIB .
(cd tests/cpp && make && ./test_bartcore)
Rscript -e 'lintr::lint_package()'
air format --check .
Rscript tools/check-rc-codoc.R .
Rscript tools/check-win-drift.R .
Rscript tools/check-doc-freshness.R .
Rscript -e 'db <- tools:::.build_news_db_from_package_NEWS_Rd("inst/NEWS.Rd"); stopifnot(!is.null(db)); print(nrow(db))'
Rscript -e 'tinytest::test_package("dbarts")'
```

- NEWS entry count unchanged (a rewrite, not an addition).
- `R CMD check --as-cran` on a tarball built from a clean copy staged outside the tree: no errors or
  warnings, notes unchanged; the vignettes build.
- Bitwise equivalence on a reference install: the three baselines benchmarks/baselines/MANIFEST names
  as current, `--bitwise` (`--strict-coverage` on the gaussian one); count the per-scenario "identical
  draws (same RNG stream)" lines against the scenario count, no "max |z|" line.
- The exact-gates.yaml gate loop in `quick` mode, since three gate scripts change.
- Step 9's replay, counts recorded.
- Under c1, the `xbart` fold check: two `xbart` calls with the same seed and different named sd use the
  same fold assignment, so the per-sd recipe in the manual holds. If they do not, the recipe says how to
  fix the folds.
- Name sweep: `git grep -n -E 'prior\.scale|(normal|linear|gp)\([^)]*scale *=' -- R man vignettes inst/NEWS.Rd`
  lists only the internal slot and its bridge read, `amplitude.prior.scale`, and NEWS sections before
  1.0-0.
- stan4bart `tinytest::test_package("stan4bart")` and bartCause `testthat::test_dir("tests/testthat",
  package = "bartCause", load_package = "installed")` against `$LIB`: pass.

## Calls made

For the ledger (the orchestrator rewrites dec-A105's entry and adds these). The four open questions
are not calls; they are listed above for the maintainer.

1. The reference k is 2, not 1: a named sd reaches the engine as anchor 2x with k fixed at 2, and an sd
   hyperprior as anchor 2c with chi(nu, 2). The bridge starts a drawn k at 2, so this starts the chain
   at the named spread, and it gives the binary default and every old spelling at the defaults
   identical engine inputs. Rejected: the unit scale (same prior, chain starts at half the named
   spread, bitwise-identical to nothing); a start-k channel in the bridge (a bridge change for an
   initialization detail).
2. The reader drops `prior.scale` and adds the sd law in force (`prior.sd.df`, `prior.sd.scale`, NaN
   under a fixed k), computed in R from the model's hyperprior and the engine's anchor. Rejected:
   adding the hyperprior to the engine's `ForestCalibration` (an engine change for a value the R model
   already holds truthfully). Whether an anchor column stays is open question (b).
3. `invchi` accepts c = 0, translated to an unnamed anchor with chi(nu, Inf), so the sd form covers
   everything the k form does and the reader's 0 can be written back. Rejected: spelling the improper
   prior only as `k = chi(nu, Inf)` (the reader would report a value the writer refuses).
4. `setLeafPrior` takes numbers only: `prior.sd` under a fixed k, `prior.sd.scale` under a drawn k,
   matching the reader's columns; the df and the hyperprior's scale move only through `setModel`.
   Rejected: an `invchi` object argument (the constructor is not in scope where the method's
   arguments are evaluated); a k writer (not asked for).
5. Under answer c1, `xbart` refuses `k` beside a named sd (question (c) holds the alternative).
6. No tombstones: `scale =` and `prior.scale` never reached main (dec-B128). A call using them meets
   R's unused-argument error.
7. Internal names keep scale (the `dbartsModel` slot, the bridge and engine identifiers): they are the
   engine's anchor, invisible to users, and read by no sister package. Rejected: renaming them, an
   unforced bridge diff.
8. `invchi`'s c has no default and nu defaults to 1.5, as `chi`'s does; `sd` takes no string form,
   unlike `k`, whose strings exist only for 0.9-x.
9. Under `monotone`, an sd hyperprior is refused, as a k hyperprior is; a numeric sd is accepted.
10. Multinomial, multi-forest and causal-forest models refuse a named sd, number or hyperprior, as they
    refused `prior.scale` (dec-A107), with the offender named in the new vocabulary; the multi-forest
    message points to `forest(sd = )` and says the two differ.
11. A hurdle fit refuses a named sd, number or hyperprior: one value cannot state a probit-latent and
    a log-y spread, and it would silently fix the zero part's drawn k. Rejected: applying it to both
    halves, as `scale =` did. Cost: the old hurdle-plus-`scale` prior has no one-call spelling; the
    same posterior is two ordinary fits. `hazard` accepts a named sd in its link's latent units.

## Records at landing

This plan's Landing note, design/nameable-calibration.md's Status line, design/prior-defaults.md's
section, and the INDEX rows. The TODO entry, if the orchestrator opens one, closes with the records
commit.

## Landing

dbarts, five commits: the package change (b53a637f), the benchmarks
(a8b36094), the standing docs and the TODO entry (7a52285f), the leaf-scale
writer's arithmetic (645c9356), and the review follow-ups (e00f7acb). stan4bart
and bartCause: no change.

Verification: tinytest 10541/10541 on the shipped build; `tests/cpp` all
green; the lint chain clean; `inst/NEWS.Rd` parses with 332 entries,
unchanged; the 25-gate `exact-gates.yaml` battery in `quick` mode all PASS.
On the reference build: the four seeded-drift snapshot files pass unchanged,
and the three bitwise equivalence compares pass (53 of 53 scenarios
identical, no |z| line; bcf and multinomial every channel identical). Step
9's replay, 26 rows: 23 draw bitwise-identical train and k values, and the
three the table predicts differ (a chi scale other than 2 and chi(nu, Inf)
beside a named scale start the chain at another spread; the old write under
chi(nu, Inf) moved draws). The writer fix moved no draw on any of these.
stan4bart 491/491 and bartCause 1055/1055 against the new library.
`R CMD check --as-cran`: one ERROR, not this slice's (an `xbart` call
without `n.threads` in test-auto-family.R trips CRAN's core limit), and the
Date-field NOTE.

Departures from the plan:

- The writer follows the later ruling: `$setLeafPrior(leaf.prior)` takes the
  creation vocabulary, has no forest argument, and replaces
  `prior.sd`/`prior.sd.scale`/`prior.mean`. A new anchor under the law in
  force goes through the engine's leaf-scale writer; a new law goes through
  the model install, with a fixed sigma put back. The engine's writer now
  derives its internal scale with creation's arithmetic, so writing a named
  sd back is bitwise inert (an approved correctness fix; draw-neutral).
- The reader's label of what `prior.sd` refers to is the attribute
  `prior.sd.of`, not a column, while the reader's shape is with the
  maintainer. The reader has 14 columns.
- A fit named by an sd hyperprior carries `sd`/`first.sd` in place of
  `k`/`first.k`, read by `extract(type = "sd")`.
- `xbart` refuses a `k` inside `leaf.prior` beside an `sd` grid, as it
  refuses an `sd` there beside either grid.
- `xbart`'s `sd` grid takes `invchi` entries as its `k` grid takes `chi`
  ones, and a named sd with no grid is a one-cell `sd` axis.
- After `setState`, the next re-anchor re-applies the model's latest written
  sd, replacing the spread the installed state carried.
- The unnamed-default rows of the mapping reach the engine with different
  inputs (no anchor against a named one), with identical draws; the k = 3
  rows are bitwise identical rather than a rounding apart.

