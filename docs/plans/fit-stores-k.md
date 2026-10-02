# fit-stores-k: a fit stores k, extract converts, and a value held fixed comes back as one number

Status: PLANNED 2026-10-02 under dec-B192, dec-B193, dec-B194, dec-B198 and dec-B199 in
[decisions.md](../decisions.md). Lands before [state-not-model.md](state-not-model.md).

agent: sonnet implementer, one; opus reviewer.
rng: NEUTRAL. R only. No sampler call, draw or default moves: the packagers read the sampler's readers once
more, and the ordinal loop keeps a channel `run` already returns.
window: pre-release. The released package has no sigma, k, sd or shape type in extract, no `summary.bart` and no
`sd` on a fit, so nothing is tombstoned; the fit's `k`, `first.k`, `sigma` and `first.sigma` are released and
keep their meaning and layout.
budget: ~950 lines (R ~330, manual, vignette and NEWS ~130, tinytest ~400, benchmarks ~30, records ~60). Plans
have run 1.5-2x low.

## Goal

Every fit carries the k its sampler recorded, under every naming of the leaf prior, and the leaf prior it ran
under. `extract` answers `"k"` and `"leaf.prior.sd"` on every fit, returns one number for sigma, shape, k or
leaf.prior.sd when the fit held it fixed and draws when it sampled it, and returns 1 for sigma on a family with
no free residual scale. `summary` tabulates what was sampled and names what was fixed on a line under the table.
The fit's `sd` and `first.sd` are gone.

## Context

- The conversion is the anchor over k, the anchor being the reader's ([`reportLeafPrior`](../../R/dbarts.R)):
  the forest total's prior sd at k = 1, in the units the forest fits. No tree count enters. Measured from
  prior-only draws at 10 and 50 trees, the spread of the forest total over that number is 0.995 to 1.004 with
  a constant leaf on every family; with a linear leaf the number is each coefficient's sd, with a gp leaf the
  marginal sd, and under a monotone constraint the total's spread runs 3 to 5 percent above it.
- Today: a fit has `k`, `first.k` only when k is drawn under a k-named prior and `sd`, `first.sd` only when
  drawn under an sd-named one; extract refuses k and sd when the leaf scale is fixed or under the other naming;
  a fixed sigma or shape is a channel repeating one number, and extract returns the repeat; an ordinal fit drops
  a drawn k; a `forest()` fit has k pinned at 1 on each forest and refuses k as fixed; no fit carries an anchor.
- R's `$` matches a list component by prefix. A new component whose name starts with `k` would answer `fit$k`
  on every fit whose k is fixed, where the equivalence harness, two benchmark scripts and
  ["expect_null(fitFixedK$k)"](../../inst/tinytest/test-nbinom.R) take NULL to mean fixed.
- The sampler and the model can disagree about a fixed value: `fixed(0.3)` names a variance and the sampler
  holds its square root, and until state-not-model lands a warm start leaves the sampler at the donor's values.
  The sampler's readers say what the draws were made under.
- Per-draw readers of the channels, none of which may see a changed layout: predict, fitted, residuals, the
  log-likelihood, [`survivalProbabilities.bart`](../../R/generics.R), the posterior-predictive draw,
  `pdbart`'s result builder, rbart's binary test and bartCause's, which is `is.null` of the fit's `sigma`.
- A fit reports in the terms the user named (dec-A105), so which leaf-scale quantity `summary` shows by default
  follows the naming.

## The mapping

| Today | After |
|---|---|
| the fit's `sd`, `first.sd`; an ordinal fit's drawn k dropped | `k`, `first.k` as the sampler recorded them, on every naming; an ordinal fit keeps `k` (it and the negative-binomial fit have no burn-in channel) |
| no record of the prior | `leaf.prior`: the reader's list as the run ended - the prior as named, what its sd is the sd of, the anchor - one list on a single forest, a list of them on several |
| a fixed value seen only as a repeated channel, or not at all | `fixed`: the scalars the sampler held fixed, read from the sampler, among `sigma`, `shape`, `k`, `resid.df`; empty when none. The channels keep the repeat |
| `extract(type = "sd")` | `"leaf.prior.sd"`: the anchor over k, on every class but rbart |
| `extract(type = "k")` refused when fixed or sd-named | the sampler's k, on every class but rbart |
| sigma, shape, k fixed: the repeat, or a refusal | one number, with no chain margin under either `combineChains` |
| sigma refused or not a type on probit, logistic, hazard, ordinal, negative-binomial and multinomial fits | 1 |
| k and the sd on a fit with several forests: refused | one named number per forest, or the forest `forest =` selects; a multinomial fit, whose forests share one prior, one number |
| a hurdle fit's k: a list of the drawn parts | a list of both parts, each draws or one number; leaf.prior.sd the same, each in its part's units |
| `summary`: a fixed sigma or shape a constant row, a fixed k or df no row, the pinned first threshold a constant row | sampled parameters only; a line under the table names each fixed one and its value |
| `summary`'s default leaf-scale row: k or sd on bart and hurdle fits, k only on the others | by naming on every class: k on a k-named fit, leaf.prior.sd on an sd-named one; either on request |
| `plot` traces a constant sigma or shape; `print.bartNegbin` gives a fixed shape as a posterior mean | the panel is not drawn; print says the shape is fixed |

## Constraints

- `leaf.prior` and `fixed` are descriptors in dec-B126's sense: on every bart, negative-binomial, ordinal and
  multinomial fit and each hurdle part, whatever the run kept. `k` and `first.k` stay channels, present only
  when drawn. No new component's name starts with `k`, `sigma`, `shape`, `first`, `fit` or `y`.
- Fixed values and the anchor come from the sampler's readers, not the model, and a fit with a variance forest
  has no fixed sigma whatever its residual prior says: extract returns its per-observation surface as today.
- The naming is the reader's: an sd-named prior is one the reader states with an sd, `invchi(df, 0)` included.
- One helper builds the two descriptors for every packager; one serves the four types to every extract method.
  [`packageMultinomialResults`](../../R/bart.R) takes no sampler today and its callers hold one.
- A fit without the descriptors, one saved by the released package, is answered from its channels where they
  suffice and refused by name otherwise.
- The manual says once that k is the sampler's, on an sd-named fit relative to the anchor and not to the data;
  that leaf.prior.sd is in the units the forest fits - the response's on a gaussian fit, log time, the log mean
  or the link's latent scale on the others; what it is the sd of under each leaf model, and that under a
  monotone constraint the total's prior spread is a few percent wider; and that on several forests it is each
  forest's own total before its amplitude multiplies it.
- Out of scope: the sampler, `run` and the saved state; rbart, whose chains end on different anchors and which
  is removed in 1.1-0; pdbart's objects; a burn-in argument to extract (`first.k` stays on the fit and
  `summary` takes it by name; `first.sd` has no successor); the consumers' own extract methods; anything under
  `src/`.

## Steps

1. Packaging: the descriptor helper, [`packageBartResults`](../../R/bart.R) and the negative-binomial, ordinal
   and multinomial packagers beside it, with the ordinal k channel from [`bart2Ordinal`](../../R/bart.R).
2. extract: one shared reader in place of [`extractLeafSpread`](../../R/generics.R) and
   [`fitAllowsKHyperprior`](../../R/generics.R); the methods' type lists; sigma's 1; the per-forest answer and
   [`refuseForestSelectionOutsideForestArm`](../../R/generics.R), which calls k a model parameter today.
3. summary, print and plot: [`scalarFields`](../../R/diagnostics.R), [`drawsField`](../../R/diagnostics.R),
   [`plotSigmaTrace`](../../R/plot.R), the fixed line and its component on the summary object.
4. `man/bart.Rd`, `man/bartBT.Rd`, `man/summary.bart.Rd`, the vignette's sentence that summary reports sigma
   and k unconditionally, and `inst/NEWS.Rd` with final spellings only.
5. Tests: one new file over every fit class, drawn and fixed, chains combined and separate, with and without
   the kept sampler, each extract answer checked against the sampler's readers; the descriptor presence rows in
   inst/tinytest/test-fit-descriptors.R. Rewrite the blocks that pin today's behaviour:
   ["a fit reports the leaf prior in the terms"](../../inst/tinytest/test-leaf-prior-k-or-sd.R),
   ["leaf-prior sd was not sampled"](../../inst/tinytest/test-nbinom.R),
   ["k was fixed, not sampled"](../../inst/tinytest/test-convergence-diagnostics.R) and the fixed-df summary
   beside it, ["with a fixed component left out"](../../inst/tinytest/test-hurdle.R), and the hurdle k loop of
   inst/tinytest/test-one-chain-dimension.R.
6. Benchmarks: the three readers of a fit's k read it by exact name, with a comment saying why.
7. Records: the Landing note, the index row, the TODO item, and a ledger entry for the calls made here: the
   two descriptor names, the per-forest answer, the pinned threshold and a fixed df on the fixed line, the
   dropped plot panel, rbart left alone and the burn-in spread losing its name.

## Verification

Against a private library:

- `cd tests/cpp && make && ./test_bartcore` and `tinytest::test_package("dbarts")` pass with no new warning.
  Mutation, run once and reported: with the fixed branch of the shared reader removed the new tests fail.
- What a fit carries changes, so every gate `.github/workflows/exact-gates.yaml` lists passes in `quick` mode.
- On a reference build the three equivalence compares report identical draws for every scenario, counted per
  scenario with no skipped one, and the four seeded-drift snapshot files pass: a key gained or lost by a fit
  shows here.
- bartCause's suite passes against a private-library chain built on this tip.
- `lintr::lint_package()`, `air format --check .`, `tools/check-rc-codoc.R`, `tools/check-win-drift.R` and
  `tools/check-doc-freshness.R` pass, each on its own exit status; `inst/NEWS.Rd` parses; `R CMD check
  --as-cran` on a tarball from a clean export.
- `git grep -nE 'first\.sd|type = "sd"' -- R man inst vignettes benchmarks` lists nothing; the report gives
  the output.
