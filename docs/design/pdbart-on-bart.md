# pdbart and pd2bart fit through bart

Status: PLANNED 2026-10-04. Plan: [pdbart-on-bart.md](../plans/pdbart-on-bart.md).
Grown slice by slice; this text covers slice 1, the fit route.

## The fit route

A call that passes data is rewritten into a `bart` call and evaluated in the
caller's frame, so an argument written unevaluated, such as `k = chi(df,
scale)`, resolves in `bart`'s prior vocabulary rather than the caller's
([`pdbart.fitData`](../../R/partialDependence.R)). The rewrite removes
pdbart's own arguments, refuses the ones pdbart sets itself (`samplerOnly`,
`test`, `offset.test`, `keepTrees = FALSE`, under either spelling), and adds
`keepTrees = TRUE, keepSampler = TRUE`. The fit is then read exactly as a fit
passed in is: each grid value is predicted from the saved trees. So a data
call, `pdbart(bart(..., keepTrees = TRUE))` and the refit of a fit kept
without its trees are one model at one seed.

A fit kept without its trees or its sampler is refit from its stored call
through the door its argument names show it came from, `bart` or `bartBT`,
with both kept ([`pdbart.refit`](../../R/partialDependence.R)). A sampler
passed with saved trees predicts from them; one without is run over every
setting's rows stacked, the one place the old `samplerOnly` route survives
([`pdbart.drawsAt`](../../R/partialDependence.R)).

Families are refused by name before anything is fit
([`pdbart.refuseFamily`](../../R/partialDependence.R)): from an explicit
`family`, from the response under `"auto"` (a `Surv` response is aft, a
count matrix or a factor of three or more levels multinomial, an ordered one
ordinal), and from a fit or sampler passed in, a hazard sampler known by its
period grid. Multinomial and ordinal stay refused; negative binomial,
hurdle, aft and hazard are refused until their slices.

## The averaged rows

The rows averaged over are the fit's own, less those it gives a weight of 0
or masks out through its active rows, each carrying its stored offset into
its prediction ([`pdbart.averagedRows`](../../R/partialDependence.R)). A
value is therefore the row mean of `predict(fit, newdata, type = "bart")`
with the variable set. pd2bart predicts each grid point as one row only when
every averaged row is that row: two predictor columns, both varied, and one
offset for every row.

Draws come back merged, each chain's in turn. When the fit keeps its chains
apart, as its `varcount` shows, `fd` gains a leading chain margin, so the
argument means the same on `fd` as on the fit's own components; the plot
methods merge it before taking quantiles.

## The translation

Every name `bartBT` takes and `bart` does not is translated, less `x.test`
and `sampleronly`, which are refused
([`pdbartBayesTreeNames`](../../R/tombstones.R),
[`translatePdbartCall`](../../R/tombstones.R)). Each warns once per session
per name and function, package callers included. `power`, `base` and
`splitprobs` become one `tree.prior = cgm()`; `proposalprobs` is set on a
copy of the caller's control or on a fresh one; `sigdf` and `sigquant` reach
`bart` under their own names, since the residual prior rides a family that
cannot be built before the response is known. A setting under both
spellings is refused. `bart`'s defaults message and its `sigdf` and `sigquant` warnings are held back
inside pdbart, their keys restored as they were found
([`holdingBartNotices`](../../R/tombstones.R)); pdbart shows its own
defaults message instead ([`notePdbartDefaults`](../../R/tombstones.R)).
