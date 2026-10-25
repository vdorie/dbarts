# Which rows a prior's scale is read from

Status: LANDED 2026-10-11 (d5b637ee); slice D LANDED 2026-10-24 (a3e65674..da962901 and 8338e66e..60746517). The plan is [response-scale-rows.md](../plans/response-scale-rows.md), which carries the
steps, the gates and the calls made in planning; this note is the standing statement of the rule.

What a fit takes from the response to place or size a prior - the gaussian and Student-t response range, a
count fit's center, the unit of a forest's sd in a model of several forests (the number a stated sd is divided by, and a default is a multiple of; see [forest-sd-unit.md](forest-sd-unit.md)), sigest's stand-in and its
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
   installed, and so does a gp sampler holding saved draws whose family has a response scale (gaussian,
   Student-t, aft, nbinom); a probit, logistic or ordinal one has none to derive and takes it.
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
9. Leaf covariates. A linear or gp leaf's covariate center and scale, and the lengthscale a gp leaf derives with
   them, are read over the rows in wherever the sampler derives them: at creation, under a mask family the rows
   the creation mask keeps; at a `setData` that derives them, over the new data's rows in; in an `xbart` fold,
   over the fold's training rows in. Unweighted, as in rule 1, and over every row where fewer than two are in,
   under every family (dec-B485), a hazard fit counting its subjects and not their period rows, so that one
   subject in reads every row however many periods it spans (dec-B489), or where the response's own read falls
   back to every row under rule 4 (fewer than two distinct values of it in, at creation, a `setData` that derives
   the range, or an `xbart` fold over its training rows), under the creation warning, which counts the rows in, or
   a hazard fit's subjects: one text with the response's where the response falls back too, as one row in always
   makes it under the families that read it, and with a hazard fit's periods where no subject is in (rule 8). The
   subjects are R's to count, from the period column's layout, each subject's rows together at periods 1, 2, 3 and
   on as the expansion writes them, and the engine is told that too few are in, a fit of one subject only where
   some of its rows are out; data that `setData` or `new()` is handed in another layout, and a period column
   `setPredictor` or `updatePredictorPerObservationJointly` could leave it in, whole or row by row, are refused,
   their subjects not being countable. Where the rows read hold fewer than two distinct observed values of a
   covariate, one value or none, two or more rows in at one value or every row at one value, a derivation is
   refused, there being no number in force to keep, unless a gp leaf states its lengthscale, which takes the
   placeholder: the one value, divided by 1. `setPredictor(updateStandardization = TRUE)` re-derives the replaced
   covariates over the rows in at that call, under the weights and the mask in force, and keeps one with too
   little, every one on a hazard fit where fewer than two subjects are in (dec-B491), warning once per sampler as
   rule 4 does, on the flag the response's keeps share. `updateStandardization = TRUE` on `setWeights` and
   `setActiveRows` re-derives every leaf covariate over the rows in after the install, keeping by the same rule,
   one row in included; a hazard fit counts its subjects there too, and fewer than two in after the call keep
   every covariate, the period column among them, under one warning once per sampler that counts subjects, R
   telling the re-derivation as it tells creation (dec-B490). Where `updateScale = TRUE` in the same call keeps
   the response's scale (rule 4), the prior in force is kept whole: every covariate keeps too, whatever its spread
   over the rows in, under one warning naming the response and the covariates, as the same rows at creation read
   both from every row; with `updateScale = FALSE` nothing reads the response and each covariate keeps or
   re-derives by its own values. A loop whose rows in change passes it while burning in and then stops, so every
   kept draw comes from one prior (dec-B465). Held otherwise, through the data object's record: a copy, a reload
   and a re-creation take it and read no row. Cut points still span every row.

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

Slice D adds the leaf covariates: on a fixture with x1 uniform on rows 1 to 100 and on [5, 15] on rows 101 to
200, which carry weight 0, a linear leaf's center and scale of x1 were 5.346 and 5.210, every row's, and are
0.518 and 0.268 over the rows in, read by [`standardizationMomentsForColumn`](../../src/bartcore/data.hpp) and,
for a gp leaf's default lengthscale, [`GPGaussianLeaf::medianPairwiseDistance`](../../src/bartcore/model.hpp)
over the same rows.

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

A leaf re-derivation decides its keep the same way, once, from chain 0's leaf, per covariate, over the rows the
weights and the sampler's copy of the mask in force leave in ([`Sampler::rederiveStandardization`](../../src/bartcore/sampler.hpp)),
and answers the covariates it kept; the bridge carries them on the predictor setters' answer, and on the row
setters' with whether the transform was kept beside them, with "warn now" from the holder's one flag, and the
R method names them in its last act. A derivation that read every row for
want of a row in marks those covariates in the engine's records (`DerivedRecords::leafFromAllRows`), which R
reads for the creation warning; a derivation with too little is the engine's typed error, worded by the bridge
with the column's name and the counts, and `setData` raises it before anything moves.

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
(landed, [response-scale-rows-c.md](../plans/response-scale-rows-c.md)).
D: a linear or gp leaf's covariate center and scale, and the default lengthscale, over the rows in, rule 9
([response-scale-rows-d.md](../plans/response-scale-rows-d.md); its first push landed 2026-10-24, dec-A226), and
`updateStandardization` on `setWeights` and `setActiveRows` (its second push, landed 2026-10-24, dec-A227).
