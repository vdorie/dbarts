# Kept draws and seeded forests keep their function when the units under them change

Status: LANDED 2026-10-07 (the conversions as a8a72c6b to f668cc86, the Gaussian-process refusals in 21bbbfd7 to 11377982).
Plan: [leaf-conversions.md](../plans/leaf-conversions.md). Rulings: dec-B200, dec-B231, dec-B233 and dec-B237
in docs/decisions.md.

A stored leaf is a set of numbers read through two sets of units. One is the response transform, a multiplier
and a shift every leaf value is stored relative to. The other, on a linear or Gaussian-process (gp) leaf, is the
covariate standardization: a centre and a scale per leaf covariate, and on a gp leaf a lengthscale. A replay
of a kept draw reads whatever units the sampler holds at the time
([`Chain::addFlatPredictions`](../../src/bartcore/chain.hpp)). So a call that re-derived either set changed
what every kept draw returned, and a forest handed from one sampler to another was read in the second
sampler's standardization. The rule now is that a draw the sampler has kept is the function that was drawn,
and a seeded forest starts at its donor's function; the numbers are rewritten when the units move.

## What was measured

Before the change, on fixtures with a response sd of about 2:

| call | what moved | by how much | now |
|---|---|---|---|
| `setData`, leaf covariate 3 x + 5 | kept draws at the same raw points, linear leaf | up to 7.9 (8.1 on two chains) | 7e-15 |
| `setData`, 40 rows appended inside every column's range | live fit on the old rows, linear leaf | up to 1.10 | 6e-15 |
| `setResponse(3 y + 10, updateScale = TRUE)`, and `setData` with that response | every kept draw, any leaf | up to 20.7; equal to 3 f + 10 | 7e-15 |
| `setOffset(-5, updateScale = TRUE)` | every kept draw | 5 | 4e-15 |
| the same re-anchor under a variance forest | kept variance draws | times 9 | 1e-15 |
| a warm start, donor centre and scale -0.019, 0.973, recipient 4.94, 2.92 | seeded fit against the donor's | up to 6.30 | 8e-15 |
| a warm start from another cut grid | the recipient's centre and scale (and gp lengthscale) | replaced by ones re-derived from its rows | kept |

A constant leaf's kept draws never moved across `setData`. With `updateScale = FALSE` nothing moved and
nothing moves.

## The rules

- One leaf's coefficients go from centre m and scale s to m' and s' by
  [`convertLinearCoefficients`](../../src/bartcore/model.hpp): each slope is multiplied by s' / s and the
  intercept gains slope (m' - m) / s. The leaf is then the same function wherever the covariate is observed. A
  row missing the covariate is read at the centre in force, so its fit moves by that intercept term; the
  manual says so. A column whose centre and scale agree is not touched, and equal standardizations run no
  arithmetic, so such a call or install is bit for bit what it was.
- A covariate that holds a single value has no spread: it is centred at that value exactly
  ([`standardizationMomentsForColumn`](../../src/bartcore/data.hpp)) and divided by a placeholder scale of 1,
  so the leaf reads zero for it on every training row, its slope adds nothing to the fit and no observation
  informs it. The state stores the scale as `NA` beside a finite centre; a state that holds 1 reads as 1.
  Whether a column is one is read from what the leaf holds when it is asked - the scale is 1 and every
  gathered row is zero ([`standardizedColumnHasSpread`](../../src/bartcore/model.hpp)) - and is not
  remembered, so a placeholder column that `setPredictor` gave values, a constant moved to another constant
  and a state installed over other rows are each what they then are. When coefficients a chain will go on
  drawing from are converted, with spread on neither side nothing moves, the fit being the intercept before
  and after. With spread on one side only the formula runs with 1 for that side's scale (dec-B329), so the
  leaf stays the function it was: onto a constant it takes the function's value there, and off one a slope
  no observation informed is carried onto the real scale (0.13 against a constant column at 1000 becomes 40
  on the scale of 307 that varying data gave), the fit moving until the sampler runs. Until 2026-10-08 the
  one-sided cases set the slope to zero instead. A kept draw is only replayed, so it is converted by the
  formula with 1 in every case, both sides without spread included.
- `setData` on a linear leaf converts the live coefficients, keeping a column without spread on both sides,
  and every kept draw
  ([`Chain::applyNewData`](../../src/bartcore/chain.hpp)).
- A re-anchor - `setResponse` or `setOffset` with `updateScale = TRUE`, and every `setData` - rewrites the
  kept mean draws and the kept variance factors into the new transform
  ([`Chain::restateSavedDraws`](../../src/bartcore/chain.hpp)), by the arithmetic a state in other units is
  converted with ([`Chain::convertStateUnits`](../../src/bartcore/chain.hpp)). Only the slots the runs have
  filled are written ([`SavedDrawSlots`](../../src/bartcore/chain.hpp)), so a sampler that keeps nothing pays
  one comparison and a store nothing was recorded into is left as it was. The live chain keeps its internal
  values, so its fit follows the new range: that is what re-anchoring a live chain means and it did not
  change. An install's own move of the transform is not such a site; the draws it brings are converted
  already.
- A warm start converts the donor's coefficients, live or from a kept draw, from the standardization the
  donor's state records into the recipient's, as live coefficients are
  ([`Chain::convertDonorStandardization`](../../src/bartcore/chain.hpp), the blocks read by
  [`readWarmStartState`](../../src/R_interface_bartcore.cpp)). The recipient keeps its centre, scale and
  lengthscale on either grid: [`Chain::rebuildLiveForestRemapped`](../../src/bartcore/chain.hpp) no longer
  re-derives them (dec-B231, derived data changes only by a call that names it). A donor that records no
  standardization is read as before. A record that does not fit the leaf is refused with nothing touched.
- `setState` is unchanged: it installs the standardization a state holds.

## One sequence whose posterior changes

A warm start from a donor on another cut grid, into a linear or gp sampler whose leaf covariate was replaced
after creation by `setPredictor`. That install used to re-derive the recipient's standardization, and default
lengthscale, from its current predictors; the recipient now keeps the one it had, which shapes its slope prior
or kernel. For a recipient unchanged since creation re-deriving gave the same numbers, so nothing moves.
Everything else is either unchanged draw for draw or, on a linear leaf after `setData` or a warm start
between unequal standardizations, the same posterior reached from converted coefficients.

One kind of fit changes too. A leaf covariate constant at a value its mean missed by rounding (0.1 over 150
rows) was given that rounding as its scale, 2.5e-16, and read about 1 on every row; a test row at 0.11 was
predicted at 1.8e14. It is now centred at its value and divided by 1 like a constant at 1000, and such a fit
draws what the same fit on a constant at 1000 draws.

## Gaussian-process leaves

A kept gp draw is a set of kernel weights. It replays only under the centre, scale and lengthscale it was
drawn with, and it has no mean term that could carry a shift of the response. So it cannot be rewritten at
`setData`, nor at a re-anchor that moves the shift; where the shift holds, the ratio is carried like any other
leaf's. A gp sampler holding kept draws refuses `setData` and every `updateScale = TRUE` by name instead of
silently changing them ([`refuseSavedGPDrawReanchor`](../../src/R_interface_bartcore.cpp), on the R methods and
on both flat C entries). Two decisions stand behind that. dec-B237 rules the change of standardization, which
is `setData`'s refusal. That a re-derived response range is refused as well, and outright, is the
coordinator's decision of 2026-10-06, recorded in the plan
([Calls made in planning](../plans/leaf-conversions.md#calls-made-in-planning)) and not yet marked by the
maintainer. The refusal reads neither the new values nor the family. One that passed the re-anchors leaving the
midpoint alone would let a line of calling code through on one sweep and stop it on the next, so a response
already in force, one of twice the spread about the same midpoint and a probit response, each of which was
accepted and the first and last of which changed nothing, are refused with the rest. It is raised before a
value is read: the state, the data object and the kept draws are what they were. With `updateScale = FALSE`
both calls are served as before, and a gp sampler that holds no kept draw (trees not kept, nothing run yet, or
the store emptied by a warm start) takes every call as before. The cost falls on one idiom: a sampler that
keeps trees and re-anchors every sweep while warming up was served when each sweep was run as a kept sample,
and now stops at the second, a draw being kept by then. Run as burn-in the sweeps keep nothing and the loop is
served; the help says so. A gp sampler's per-row fits hold no kernel, so
a warm start copies them on the donor's grid and starts them at zero on another, as before.

## Where it is tested

[`testLinearCoefficientConversion`](../../tests/cpp/test_model.cpp),
[`testNoSpreadStandardization`](../../tests/cpp/test_model.cpp),
[`testReanchorKeepsSavedDraws`](../../tests/cpp/test_model.cpp),
[`testLinearLeafSetDataConversion`](../../tests/cpp/test_moves.cpp),
[`testWarmStartStandardization`](../../tests/cpp/test_state.cpp),
[`testNoSpreadLiveConversion`](../../tests/cpp/test_state.cpp), and from R
["kept draws across a re-anchor"](../../inst/tinytest/test-leaf-conversions.R),
["live coefficients across a covariate without spread"](../../inst/tinytest/test-leaf-conversions.R),
["gp leaves that hold kept draws"](../../inst/tinytest/test-leaf-conversions.R), and through the flat entries
["a sampler with gp leaves that holds saved draws refuses a re-derived range"](../../inst/tinytest/test-capi.R).

## Not here

The centre, scale, lengthscale and response range as records on the data object, `setData`'s arguments for
holding them, and a state install that compares a stored standardization and converts instead of installing
it: the TODO item state-frame-prior.
