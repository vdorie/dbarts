# interface-review

agent: sonnet readers, one dimension each; fable adjudicates
rng: neutral (findings only; any fix becomes its own gated item)
budget: findings recorded here; deliverable is a fix-before-release
        list, a 2.0-wishlist, and taste calls batched to VD

Review 2 of the retrospective program (retrospective-reviews.md).
Runs now because 1.0-0 freezes the R surface at submission: renames
and default changes cost nothing today and are breaking changes the
day after.

## Goal

Every exported entry point reads as one library rather than fifteen
years of accretion: shared concepts share names and defaults,
generics behave uniformly across fit types, failures speak the
user's language, and a newcomer can get from install to prediction
on the shipped docs alone.

## Scope

The exported surface per NAMESPACE: dbarts, bart, bart2,
dbartsControl, dbartsData, dbartsPriors, the dbartsSampler
reference class, rbart_vi, xbart, pdbart/pd2bart, predict/fitted/
extract/residuals/plot/print/summary methods, plotTree,
updatePredictorPerObservationJointly, makeind/
makeModelMatrixFromDataFrame/makeTestModelMatrix, guessNumCores.
dbarts.h is out of scope (its compat contract is already frozen);
internal BCF plumbing is out of scope (unexported, bartCause-gated).

## Method

Four read-only readers, one dimension each, findings only:

A. Argument and default coherence: the full argument matrix across
   entry points - names, casing, defaults, ordering; the same
   concept under different names (or worse, different concepts
   under one name); defaults that disagree between bart, bart2, and
   dbartsControl for the same knob; documented vs actual defaults.
B. Generics and methods: predict/fitted/extract/residuals/plot/
   print/summary across bart, rbart, and the sampler - argument
   conventions (type/value vocabularies), return shapes and
   dimension orders, combineChains coherence, behavior when the fit
   cannot support the request (keeptrees off, no test data), and
   whether the sampler's mutable surface matches its documentation.
C. Error-message and validation quality: the user-reachable stop()/
   warning()/Rf_error sites - do they fire early, name the
   argument, say what was received and what was expected; silent
   coercions and silently-ignored arguments; NA/type/length edge
   validation at the R boundary.
D. New-user walkthrough on shipped docs only (installed man pages,
   README, examples - no source reading, no docs/): regression,
   binary outcome, weights and offset, prediction on new data,
   the sampler as a conditional model, rbart_vi grouping, xbart
   cross-validation. Every stumble recorded verbatim at the moment
   it happens.

Fable adjudicates into: fix-before-release (cheap, surface-stable,
uncontroversial), 2.0-wishlist (breaking or opinionated), and taste
calls for VD (batched in chat with recommendations). Fixes become
their own items under the standard gates.

## Constraints

- Findings only; the tree does not change under this review.
- Readers are read-only, FOREGROUND only, and spawn nothing.
- Reader D judges only what ships in the installed package.

## Verification

- The three lists recorded in this file's Status; VD's decisions on
  the taste calls recorded beside them.

## Status

- 2026-07-08: started; readers A-D dispatched.
