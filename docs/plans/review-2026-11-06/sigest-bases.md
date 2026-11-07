# Review 4, slice S-G: the starting sigma reads the forests' bases

Status: built, in review (branch review4-sigest-bases).

agent: implementer (opus), one opus review with mutants. rng: posterior-changing for gaussian fits with a basis
outside the span of the intercept and the predictors and no stated sigest; bitwise for every other fit. Tier: changes
draws (a default; R code only).

Finding multiforest-03 of the fourth whole-repo review, ruled by dec-B603 ("I think (a) sounds right"): with no
`sigest`, the regression whose residual sd centers the noise prior takes every forest's basis beside the predictors.
The design note is [The forests' bases](../../design/starting-sigma-sensitivity.md#the-forests-bases).

## Fixes as built

**The columns.** [`sigmaBasisColumns`](../../../R/utility.R) turns a data object's bases into the regression's added
columns, one per forest with a basis: a forest of numbers its column, a forest of levels the indicators of every
level but its first, which for the two levels a block holds is the second's, the one column a two-level factor
predictor enters as ([`emittedLevels`](../../../R/utility.R)). [`startingSigmaDesign`](../../../R/utility.R) takes
them as `bases = NULL` and puts them after the predictors' columns: on the dense route by `cbind` after
[`denseSigmaDesign`](../../../R/utility.R), on the sparse route as columns of `X` with `dense` TRUE and `term` NA, so
the sparse front end centers and pairs them as any column a caller gave dense. Collinearity with the intercept, a
predictor or another forest's basis is left to lm.fit's pivot and to the sparse routine's pivoted Cholesky, whose
rank sets the residual degrees of freedom on both routes; no column is dropped by hand.
[`estimateStartingSigma`](../../../R/spec.R) and the default `design` of
[`estimateSigmaFromLinearModel`](../../../R/utility.R) pass them; the verbose line's column count reads the design
and so counts them. xbart's positional `startingSigmaDesign(data@x, route)` calls are untouched (xbart refuses
bases).

**Every route.** [`resolveSamplerSpec`](../../../R/spec.R) installs a declaration's blocks (`bases`, `basisLevels`,
and the dropped `basis.scale` entry of a replaced forest) above `enforceWeightPolicy`, ahead of the estimate, where
it installed them after; so `dbartsSpec()`'s forests list, and its replacement of a data object's own bases, are
read by the estimate as `dbarts()`'s already were. Between the two positions only `declaredBases` read the bases,
and it reads the same blocks.

**The fallback line.** Where the design leaves no residual degrees of freedom the verbose line reads "the linear
model on N columns (factors as indicators, each forest's basis among them)" when a forest has a basis, and as before
otherwise.

**Help and comments.** dbarts.Rd's and bart.Rd's `sigest`: "the linear regression of the response on the predictors
and on each forest's basis (a factor's as the indicator of its second level)". bartBT has no basis door and xbart
refuses bases, so their pages keep the sentence. The comments over `estimateStartingSigma` and
`estimateSigmaFromLinearModel` say the same. No NEWS line: 0.9-34 had no bases.

## RNG class and reach

Posterior-changing for gaussian fits with a basis outside the span of the intercept and the predictors and no stated
sigest or fixed sigma; gaussian is the only family that takes a basis and estimates a starting sigma. Bitwise for
every other fit: no basis, a latent family, a stated or fixed sigma, a data object carrying its own sigma into
`dbartsSpec()`, and a basis already among the predictors (bartCause's bcf puts z among x).

## Tests

["the forests' bases enter the regression"](../../../inst/tinytest/test-starting-sigma.R), 30 results in the file
(22 new):

- dense route, `expect_identical` against `summary(lm(...))$sigma`: a factor basis (`y ~ x1 + x2 + factor(z)`) bare,
  with weights and with an offset; a numeric basis w likewise; three forests, one of each kind; the data object's
  own bases (`dbartsData(bases = )`); a `bart()` fit's `sigest`.
- a basis among the predictors, z with `factor(z)` and w with w: identical to the fit with no basis.
- sparse route by `sigest = "sparse"` (a fit with a basis refuses a sparse x), `expect_equal` at 1e-12: the factor
  basis bare and with weights, the numeric basis with an offset, and z among the predictors with `factor(z)`; and
  `dbartsSpec()` on a dgCMatrix x with a factor basis.
- `dbartsSpec()`: a forests list over a data object with no bases, and one replacing the object's `factor(z)` by w,
  identical to lm on the declared basis.
- a stated `sigest` taken as given; the verbose line's count (100 rows, 3 columns); the fallback line's text on 5
  rows and 4 columns.

Three fixtures elsewhere state a starting sigma, so that what they pin is not the default estimate's (Calls made, 9).

Mutants, patched into the namespace, each run against the file: the bases dropped (18 of 30 fail); the install
moved back after the estimate (2 fail, the two `dbartsSpec()` declarations); the residual degrees of freedom from
the column count on the dense route (6 fail) and on the sparse route (1 fails, the aliased sparse case).

## Gates

On the arm64 macOS development machine, at most two cores.

- The tinytest suite in one process at its default ncpu, shipped build: 41417 results, 0 failed. Its first run
  failed 3, the fixtures of Calls made, 9.
- bcf-equivalence.R, reference build, `--bitwise` against a6e83f38: exactly the 16 gaussian scenarios move, the 3
  latent ones bitwise. Partition (main's estimate stated as `sigest` at the tip): 19 of 19 bitwise. Oracle: 16 of 16
  estimates `identical()` to lm on the predictors and the basis column, with weights and offset. Re-recorded whole;
  the new file reproduces 19 of 19 bitwise, and `--cross-host` against it passes.
- equivalence.R, reference build, `--bitwise --strict-coverage` against a6e83f38: 54 identical, bart2twoforest
  moved (209 summaries, largest |z| 3.32, one above 3), 55 compared, 0 skipped. Partition: bart2twoforest at main's
  estimate 1.9650 reproduces a6e83f38 bitwise. Oracle: its estimate 1.9107 is `identical()` to lm. Merged
  re-record; the new file reproduces all 55 in one run.
- multinomial-equivalence.R against 80b1c8d4, reference build: 11 of 11 bitwise; `--cross-host` passes.
- The four seeded-drift snapshot files, reference build: 27 results, 0 failed.
- exact-gates.yaml's commands with `quick`, shipped build: 46 of 46 pass, the BCF gates among them (bcf-exact,
  bcf-exact-weak, bcf-exact-restricted, bcf-latent-exact).
- `R CMD check --as-cran --no-manual` on a tarball built from a clean export: one NOTE, the incoming-feasibility
  one (the Date field, and an expired certificate on a URL in bartBT.Rd and pdbart.Rd, pages the slice does not
  touch); the tests 215 s.
- `lintr::lint_package()` against the shipped library: no lints. `air format --check .`, tools/check-rc-codoc.R (53
  methods), tools/check-win-drift.R and tools/check-doc-freshness.R: each exit 0.
- bartCause's test-14-bcf.R against the slice's shipped library (z among the predictors): 20 tests, 186
  expectations, 0 failed.

## Calls made

1. A forest of levels adds the indicators of every level but its first, `block[, -1L]`, in place of the literal
   second column the plan names. The two are one column for the only size a basis takes (two levels); for a wider
   block, should one ever be allowed, this is lm's treatment-contrast span and the literal one would drop levels.
   Alternative: the second column; or `emittedLevels` over the block's column counts, as a factor predictor's
   columns are chosen.
2. The added columns are made once, as a matrix, by a helper of their own (`sigmaBasisColumns(bases, levels)`), and
   `startingSigmaDesign` takes that matrix. Alternative: `startingSigmaDesign` takes the data object or its bases
   list; that changes the predictor-source argument xbart passes positionally and builds the columns once per route.
3. The added columns carry no names on either route. Nothing reads a design column's name; the verbose line counts
   them. Alternative: name them by the forest's label.
4. The install's new place is directly above `enforceWeightPolicy`, after the family is resolved, as amendment 1
   says. Alternative: at the top of `resolveSamplerSpec`, beside `checkBasisRecord`; nothing between reads the
   bases either way.
5. The fallback line names the bases only when a forest has one. Alternative: the longer text on every fit,
   naming bases a fit without them does not have.
6. The partition injects main's estimate as a stated `sigest` into a scratch copy of each harness, the values
   recorded by tracing `estimateStartingSigma` while the same scratch copy ran on a reference install of main.
   Alternative: compute main's value at the tip from the no-basis path, which would lean on the slice's own code.
7. The oracle traces `estimateStartingSigma` at the tip and fits `lm()` on the data object's predictors and each
   forest's basis column (the second level's indicator for a forest of levels), its weights and offset.
   Alternative: rebuild each scenario's data by hand, which duplicates the harness.
8. The equivalence baseline takes bart2twoforest fresh, merged into a copy of its predecessor in scenario order with
   the other 54 carried as identical objects, as earlier single-scenario re-records did; the bcf baseline, whose
   harness runs all 19 at once, is recorded whole, and its three latent scenarios are identical objects to the
   predecessor's. Alternative: re-record all 55, a quarter hour at two cores for 54 identical objects.
9. Three fixtures whose pins the new estimate moved state a starting sigma: test-forest-door.R's kind `three`
   (`sigest = 0.4`), where at the estimate 0.321 a power written on the dose forest draws its ten sweeps as the
   default's (at three of seeds 1 to 8, the file's 5 among them), and moves them at main's 0.414; test-control-migration.R's three
   forests (`sigest = 0.5`), where at 0.354 the effect forest's moves draw six sweeps as no moves', and move them at
   every value tried from 0.3 to 1 and at seeds 1 to 3; and test-forest-scale-record.R's held swap pair
   (`sigest = 1`), whose swapped sampler was created with `dose` as its basis and so estimates on it, where the
   sampler it is compared with has none. The estimate is creation's, as a predictor swapped in later keeps it.
   Alternative: another seed, or longer runs, for the two moves-a-draw pins.

## Deviations from the plan

- A forest of levels enters as `block[, -1L]`, every level's indicator but the first, where the amendment names the
  block's second column; for the two levels a block holds they are the same column (Calls made, 1).
- Three fixtures in other test files state a starting sigma (Calls made, 9), where the plan's critique found no
  tinytest pinning a basis fit's default sigma: two pinned that a setting moves a draw at a seed the new estimate
  happens to leave unmoved, and one compared samplers created with different bases.
- The fallback line's text is pinned by a test the plan did not name, and the design index row for the sensitivity
  note gains a clause naming the new section.
- The plan's "a fit with no basis unchanged bitwise" is held by the file's existing identities against lm and by the
  bitwise compares (the 54 equivalence scenarios with no basis, the snapshot files); no new no-basis case was
  added.

## Not verified

- Off arm64 macOS: no x86 or Windows run; the slice is R code only, and the dense identities are lm's own routine.
- Frequency-weights S5's weighted read-offs, which build on this slice (its residual degrees of freedom being the
  sum of the weights less the rank the routine returns), are not built.
- No posterior study beyond the review's 20 seeds, and no multi-forest package other than bcf.
- The sparse route's cost of a dense basis column at large n (n doubles a forest) was not measured.
