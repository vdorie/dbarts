# Review 4, slice S-E: documentation

Status: LANDED on branch review4-docs (documentation only; no code and no behavior changed).

Each item is a finding of the fourth whole-repo review; the pages changed follow.

## Items done

- state-02, state-03: the `updateState` item and the seeding contract state what the first read of `state` does under each setting and that the mutators store after each write by default (dbartsSampler-class.Rd, dbarts-embedding.Rd, dbartsControl.Rd).
- state-04: a probit, ordinal or aft latent is on the offset-inclusive scale; an outer block that owns the offset regresses the latent less `getFitsWithoutOffset()` (dbarts-embedding.Rd, dbartsSampler-class.Rd twice, the recipe 6 prose in dbarts-as-a-component.Rmd).
- state-05: recipe 1 builds on `2 * y - 1` and says why; recipe 3 says the composed prior's center is K times the build response's midrange and how to set it; dbarts-embedding.Rd says the center is inherited (the vignette, the page, test-embedding-recipes.R pinning `prior.mean`).
- state-06: the aft `status` is written to the sampler's data object (dbartsSampler-class.Rd).
- state-07, state-V-01: creation handed `active` (dbartsSampler-class.Rd).
- families-03, families-V-01: the `offset` default at predict (bart.Rd, bartBT.Rd).
- docs-03, docs-05, docs-08 and docs-V-02: NEWS UPGRADING and padding sentences (NEWS.Rd, bartBT.Rd, na.keepPredictors.Rd).
- docs-06, docs-13: xbart.Rd `power` and the loss's `weights` argument; the seed comment in working_with_saved_trees.Rmd; the package page's See Also; the Matrix guard on the dbartsMixedMatrix.Rd example; the `binaryOffset` component text.
- docs-07: the forest column of `getTrees` is present only with more than one forest (four texts and the comment above the method).
- docs-09: the follow-up, limitation, unadjudicated and in-this-version clauses are cut to the refusal or limit they qualify (bart.Rd, dbarts.Rd).
- docs-10 (man, vignette and `$getShape` docstring share): "host" becomes the outer sampler; bart.Rd says "the temporary sampler it builds". dbarts.h and the R strings stay with their own slices.
- docs-11: every one-argument `\eqn{}` holding a TeX macro gets a text form (85 in ten pages).
- docs-12: dbarts-deprecated.Rd lists `sigest` beside a fixed residual scale and the three control slots.
- multiforest-02, ingest-03, rgen-07, aux-04: forest.Rd, dbartsData.Rd and the `resolveForests` header comment; extract.dbartsSampler.Rd; bartBT.Rd `type` (loglik for aft); makeind.Rd `drop`.

## Left for a maintainer call

docs-01, docs-02, rgen-05 and multiforest-03.

## Deviations from the plan's texts

- state-05: the recipe 1 sampler is renamed `outer` (docs-10) as the plan says; recipe 3 gets the text only, its chunk unchanged.
- families-V-01: a `predict` offset on a `bartBT` fit is added to the `binaryOffset`, not substituted for it, so the page says that.
- docs-09: rbart.Rd holds no such clause; nothing changed there.
- docs-13: the first tree printed under seed 56 has 7 nodes on this build, not 5, so the comment says 7.
- docs-12: the three slot rows are written as one item.
