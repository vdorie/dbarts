# Which rows a prior's scale is read from

Status: LANDED 2026-10-11 (d5b637ee). The plan is [response-scale-rows.md](../plans/response-scale-rows.md), which carries the
steps, the gates and the calls made in planning; this note is the standing statement of the rule.

What a fit takes from the response to place or size a prior - the gaussian and Student-t response range, a
count fit's center, the unit of a forest's default sd in a model of several forests, sigest's stand-in and its
floor, and a hazard fit's periods - is taken from the rows in the likelihood when the sampler is created. A
survival (aft) fit takes no weights, so it has no row out at creation and its range is every row's then; a mask
installed afterwards is where the rule reaches it, at a re-derivation. A row is IN when its weight is positive
and no mask switches it off, OUT otherwise. A row of weight 0 stays in the design: it occupies a leaf and
receives a fitted value, so a fit with zero weights is close to the fit on the remaining rows without being it.
What it no longer does is move the numbers the prior is stated against.

## 1. The rule

1. Rows. The range, the count center and the unit of a forest's sd are computed over the rows in at creation,
   unweighted: a row of positive weight counts in full.
2. Held. No later call moves them unless it carries `updateScale = TRUE`; `setWeights` and `setActiveRows`
   without it, `setState` in either form, a copy, a reload and a warm start do not. So the three ways to leave a
   row out differ by when it left: at creation (the numbers follow the rows in), by `setWeights` or
   `setActiveRows` afterwards (the numbers stay creation's, unless the call asks), by `subset` (the row is gone,
   cut points included).
3. Re-derived. `updateScale = TRUE` on `setResponse`, `setOffset`, `setWeights` and `setActiveRows` reads the rows
   in at that call; on the last two it comes after `updateState`, installs first and reads after, whether or not
   a row moved, and leaves the test offset alone. `setData` replaces the data and derives as creation does.
   Several forests, and a multinomial sampler on `setActiveRows`, refuse `updateScale = TRUE` before anything is
   installed, and so does a gp sampler holding saved draws.
4. Too little to read. Where a row is out and the rows in hold fewer than two distinct values (of the response
   less its offset for a range, of the count for a center):
   - at creation and at `setData`, each number is taken over every row, sigest included, under one warning of
     class `dbartsScaleFallbackWarning`; the starting sigma records no stand-in for that cause;
   - at a re-derivation through `setResponse`, `setOffset`, `setWeights` or `setActiveRows`, the transform and the
     center in force are kept: the call completes as with `updateScale = FALSE`, bit for bit, the weights or the
     mask installed, and the R methods warn once per sampler with the same class.
   A constant response with every row in is not this case: it takes the window of width 1 about the constant,
   at creation and at a re-derivation, with no warning.
5. Nothing out, nothing changed. With every row in every number has the bits it had before.
6. sigest stays the weighted regression over the rows in and is re-estimated by no call. Where that regression
   leaves no residual degrees of freedom, the sd that stands in for it and the floor under it read the rows in,
   unweighted.
7. `xbart`: each fold's range and sigma follow 1, 4 and 6 over the fold's training rows; the call warns once and
   a fold falls back in silence.
8. Hazard periods. The default periods, and under `breaks = K` its quantiles, are read from the times of the
   subjects in the likelihood; stated boundaries are checked against every subject. A subject at weight 0 is
   expanded on that grid, its time in the interval it falls in; one past the last period is at risk in every
   period with no event. With no subject in, the periods are every subject's under one warning.

## 2. What read which rows

Before the slice every one of these read every row (measured on a 200-row fixture, rows 101 to 200 at weight 0
with their response 50 higher):

| quantity | read by | before | rows in |
| --- | --- | --- | --- |
| response range | [`GaussianResponse::rescale`](../../src/bartcore/model.hpp), inherited by Student-t, and by aft, which meets a row out only at a re-derivation under a mask | [-2.165, 56.371] | [-2.165, 6.029] |
| count center | [`NBResponse::computeShift`](../../src/bartcore/model.hpp) | 0.2964 | 0.9895 (half the rows structural zeros) |
| unit of a forest's sd | [`scaledResponseSd`](../../src/bartcore/chain.hpp) | 25.25 in y | 1.466 |
| sigest's stand-in | [`floorSigmaEstimate`](../../R/utility.R) | 25.245 (three rows in, two predictors) | 0.890 |
| hazard periods | [`resolveHazardGrid`](../../R/dbarts.R) | 7 periods, 167 person-period rows | 5 and 134 |

Already over the rows in, and unchanged: sigest's regression itself ([`residualStandardError`](../../R/utility.R),
[`sparseResidualStandardError`](../../R/utility.R)), and every probit, logistic, ordinal and multinomial fit,
which take nothing from the response.

The named-sd case shows what is at stake. A stated `sd` is absolute, held as `k = k.scale / sd`. With the range
read over every row `k.scale` was 29.27 where the rows in give 4.10, so a leaf sd named at 0.65 became a prior
centered at 52 where the rows in sit near 2.

## 3. Why membership and not weight

The rows read are those of positive weight, each counting in full, and not the weights themselves. A weight is
a precision; the range and the center are properties of which responses the likelihood can see. Weighting the
range would make it vary with the scale of the weights, and a weight of 0.001 would nearly drop a row where a
weight of 0 drops it entirely, so the step is at 0. Student-t serves its Gaussian part `w * lambda * a`; it
counts rows in by `w * a`, since a product that underflowed would make served precision and membership disagree.

## 4. The keep and its status

The engine decides the keep once, in the sampler, from chain 0's response model, by the read the model's own
`rescale` uses ([`GaussianResponse::readRowsIn`](../../src/bartcore/model.hpp)). On a keep it calls every chain
with `updateScale` false and skips the anchor record, so the keep is that call, literally; no chain or response
model grows a branch. The decision travels as a returned status, [`Sampler::setResponse`](../../src/bartcore/sampler.hpp)
and `Sampler::setOffset` answering `ResponseUpdateResult::scaleKept`, and `Sampler::rederiveScale`, which
the bridge calls for `setWeights` and `setActiveRows` after their install, with the offset in force; the engine never warns.

- R: the bridge returns a length-one logical, "warn now", true on a sampler's first keep, and the R methods raise
  the warning as their last act ([`warnScaleKept`](../../R/bartcore.R)), after the mirror and the stored state, so
  a warning turned into an error leaves the object agreeing with the engine. The flag lives with the engine
  handle: a copy or a reload is a new sampler and warns again, `setData` is not.
- Flat C entries: [`dbarts_sampler_setResponse`](../../inst/include/dbarts/dbarts.h) and `dbarts_sampler_setOffset`
  return `DBARTS_SWAP_SCALE_KEPT` (2) on a keep, `DBARTS_SWAP_DONE` (1) on any other swap and `DBARTS_SWAP_NONE`
  (0) where the coupling admits none. The value is named by a fourth enumerator list, `dbarts_swap_status`,
  folded into the ABI token, so the hash moves and the version pair stays 1/0 before the first release. A caller
  tests for nonzero, not for 1; a keep warns nothing.

At creation and `setData` R also tests rule 4, for the warning and for sigest: one function,
[`scaleRowsFallBack`](../../R/spec.R), called by `resolveSamplerSpec`, `bartcoreSamplerSetData` and `xbart`; the
engine reads every row there by the same rule.

## 5. Creation under a mask

Of the families whose 0/1 weights become the mask only the count family takes anything from the response. Its
creation is the install of the mask with the engine's re-derivation against the offset it already holds, in one
bridge call ([`installWeightsMask`](../../R/bartcore.R)): the mask entry called with `updateScale` true, the same
call `$setActiveRows(mask, updateScale = TRUE)` makes, which
leaves the sampler where the bridge entry `$setOffset` uses, called with that offset and `updateScale` true,
leaves it; then the record. It draws nothing and touches nothing else: the data object, the test offset
in R and in the engine, and the stored state are as creation or `setData` left them. In rule 4's case the
mask is installed without the re-derivation and the every-row center stands. A sampler built by hand from
[`dbartsSpec`](../../R/spec.R) and `new("dbartsSampler", ...)` is created with every row in; calling
`$setActiveRows(spec$active, updateScale = TRUE)` on it makes it `dbarts()`'s, `identical()` in every reader and
draw, and the help documents that call and the loop form. It does not keep the every-row center. Where the mask
keeps fewer than two distinct counts the by-hand route has warned twice, `dbartsSpec` and the keep, and has spent
the sampler's one warning, so a later keep is silent on it and audible on `dbarts()`'s.

## 6. The hazard grid

[`hazardGridTimes`](../../R/dbarts.R) hands the grid the times of the subjects in; both callers of
[`expandDiscreteTimeHazard`](../../R/dbarts.R) pass it as `gridTime`. A subject past the last period gets the
response 0 on every row. A hazard fit with every subject in, with stated boundaries, or whose subjects at weight 0
hold only times some subject in holds, is unchanged.

## 7. What moves and what does not

Every sampler created with every row in is bit for bit what it was. The posterior changes for a sampler created
with rows out: a gaussian or Student-t fit with an extreme on a row out, a model of several forests with any row
out, a count fit under 0/1 weights, `xbart` with a weight of 0, any fit where the rows in leave sigest no residual
degrees of freedom, a hazard fit with a subject at weight 0 off the grid, and `updateScale = TRUE` on the two
response setters while rows are out. Of the recorded baseline scenarios one moves, `zeroweights`.

## 8. The four slices

A: the survival re-anchor reads observed times (landed). B: this note, creation. C: `updateScale` on
`setWeights` and `setActiveRows`, which extends the keep to them and documents the by-hand count sampler's call
(planned in [response-scale-rows-c.md](../plans/response-scale-rows-c.md)).
D: a linear or gp leaf's covariate center and scale, and the default lengthscale, over the rows in, which waits on
the record slices of the state.
