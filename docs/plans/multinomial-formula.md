# multinomial-formula

agent: sonnet (a pure R-SURFACE arc: the engine is untouched, so no draw can
  move. The work is formula ingestion for family = "multinomial" - a factor LHS
  (and, if multinomial-counts has landed, a cbind(...) count-matrix LHS) - plus
  retaining the terms so predict codes a data.frame newdata. The care is
  matching the established bart2 formula idiom and threading levels/terms
  exactly, not any numerics).
rng: RNG-NEUTRAL, throughout. The engine, the bartcore .Call surface, and
  bart2Multinomial's run sequence are UNTOUCHED; only the response/predictor
  INGESTION in front of them changes. Every commit is trivially bitwise on all
  three standing anchors (single-forest equivalence 22/22 vs equivalence-ac6ec2c,
  BCF 5x6 vs bcf-equivalence-99205ee, the multinomial fixture vs its current
  baseline) and on the C7 reproduction gate - a formula fit and the equivalent
  matrix fit produce IDENTICAL draws at the same seed (same labels, same coded
  design). Run the anchors per discipline; a move means ingestion perturbed the
  design matrix or the seed path.
window: NONE. family = "multinomial" already ships (bb29d00); this arc lifts the
  formula refusal on it and adds the ingestion behind it. No new family slot, no
  dbarts.h change.
budget: ~200-260 lines. C1 ~200 (formula parsing for the factor LHS + the
  count-matrix LHS if counts landed, terms/xlevels retention, predict data.frame
  coding, tinytest). C2 ~40 (docs + Rd examples + pkgdown check). Chiefly
  R/bart.R, R/generics.R, man/bart2.Rd, inst/tinytest/test-multinomial-surface.R.
  No src/ change; no --preclean, no tests/cpp.

## Goal

bart2(family = "multinomial") accepts the FORMULA interface: a factor response
LHS (bart2(y ~ x1 + x2, data = df)) and, when multinomial-counts has landed, a
count-matrix LHS (bart2(cbind(c1, c2, c3) ~ x, data = df), the glm binomial
idiom). The predictor matrix is built from the RHS via the standard model-matrix
machinery; the fit retains its terms so predict codes a data.frame newdata. This
matches every other bart2 family's formula surface and closes the formula
follow-up the C7 landing filed (docs/plans/multinomial.md C7 scope-narrowings;
docs/design/multi-forest-models.md; TODO).

## Ordering

Lands AFTER multinomial-counts (that arc's Q1): the count-matrix LHS branch here
is gated on the count path existing. If VD routes formula first, C1 ships the
factor LHS only and the count LHS is added when counts land (an additive elif on
the LHS-type branch, not a rewrite). This plan is written for the counts-landed
order; the factor-only fallback is noted at C1.

## Context (seams, read in code)

- The current refusal: the multinomial branch rejects a formula or dbartsData
  outright (bart.R:466-470, "supports the matrix interface only this arc") and
  reads the response from the `data` positional as a factor (472-499), then calls
  bart2Multinomial (501-513). bart2Multinomial (654-718) takes the validated
  factor y, derives labels = as.integer(y) - 1L and K = nlevels(y) (667-668),
  builds the host sampler with the LABEL VECTOR as its response (686-696), and
  runs bartcoreMultinomialSampler + bartcoreRun. Nothing below the ingestion
  needs to change - the labels and the coded design x are the only inputs.
- The standard formula parse (data.R:454-522): model.frame -> `model.response(
  modelFrame, "numeric")` (464) COERCES the response to numeric, discarding
  factor levels - so the multinomial path CANNOT reuse dbartsData for its
  response; it must run its own model.frame/terms extraction that KEEPS the
  factor. The RHS -> makeModelMatrix(modelFrame[termLabels]) (498-509) is exactly
  the predictor matrix to reuse; the `test` term is pulled from the formula data
  the same way (511-522).
- Family is NEVER auto-detected for multinomial: dbarts() resolves "auto" ->
  probit for a 0/1 response (dbarts.R:285-294) and multinomial is not in its
  family list; bart2 dispatches multinomial only on the EXPLICIT family =
  "multinomial" (bart.R:429). A multi-level factor LHS under "auto" stays
  gaussian/probit as today - Q1.
- predict.bartMultinomial (generics.R:467-492) codes newdata through
  validateXTest(newdata, object$fit$data@x) (482) - correct for a MATRIX newdata,
  but a data.frame newdata from a formula fit needs the terms to build the model
  matrix (expand factors against the fit's xlevels) BEFORE coding. The fit must
  retain terms + xlevels; predict builds the model matrix from a data.frame
  newdata first. This mirrors predict.bart's formula handling.
- packageMultinomialResults (bart.R:770-814) and shapeMultinomialChannel
  (728-...) thread levels(y) onto the K margin and predictor names onto the p
  varcount margin; the formula path supplies the same levels (factor levels or
  count-matrix colnames) and colnames(model.matrix) as predictor names.

## Design

Extend the multinomial branch (bart.R:429) to accept a formula:

1. Remove the 466-470 refusal. When `formula` is a formula (or a one-sided
   response is given), run a self-contained model.frame/terms extraction (NOT
   dbartsData, which numeric-coerces the response): build the model frame from
   formula + data + subset + na.action, pull the response with model.response(mf)
   (NO type coercion), and the RHS predictor matrix with makeModelMatrix on the
   term labels (the data.R:498-509 recipe).
2. Branch on the response type:
   - a FACTOR (or character coerced to factor) LHS -> the existing factor path
     (labels, K, levels), unchanged below the extraction.
   - a MATRIX LHS (cbind(...), only if multinomial-counts landed) -> the count
     path (counts, K = ncol, levels from colnames), routing to the count analog
     of bart2Multinomial. Refuse the matrix LHS with a clear "requires the
     count-matrix likelihood" message if counts have not landed (factor-only
     fallback).
3. Retain the terms + xlevels on the fit (a $terms / $xlevels field, or reuse the
   host sampler's stored terms) so predict can code a data.frame newdata.
4. predict.bartMultinomial: if newdata is a data.frame and the fit carries terms,
   build the model matrix from newdata via the terms (delete.response, expand
   against xlevels) before validateXTest; a matrix newdata keeps the current
   path. Fitted/extract/print are unchanged (they read the stored probabilities).

The matrix interface (bart2(x.train, y.train, family = "multinomial")) is
untouched - both entries converge on the same labels/counts + coded design, so
the C7 reproduction gate (public == internal, bitwise) extends to formula ==
matrix at the same seed.

## Commits

C1. Formula ingestion, one gated commit (RNG-neutral). Sub-parts: lift the
   refusal; the self-contained factor-preserving model.frame/terms extraction;
   the LHS-type branch (factor now; count-matrix if multinomial-counts landed,
   else refused by name); terms/xlevels retention; predict data.frame coding;
   tinytest - a bart2(y ~ x1 + x2, data = df, family = "multinomial") fit
   produces the same probabilities as the equivalent matrix fit at a fixed seed
   (formula == matrix reproduction, the C7 pattern), predict on a data.frame
   newdata reproduces the fit-time test channel, a factor with unused levels is
   handled, and (if counts) a cbind(...) LHS reproduces the count matrix fit.
   Files: R/bart.R, R/generics.R, inst/tinytest/test-multinomial-surface.R. Gate:
   full tinytest (grows; no regen - RNG-neutral) + R CMD check man + all three
   standing anchors identical + the multinomial fixture identical + air. Size: L.
   Abort: a formula fit diverges from the equivalent matrix fit on the same seed
   (the ingestion changed the coded design or the seed path).

C2. Docs + Rd. Update man/bart2.Rd (the formula interface now supports family =
   "multinomial"; a factor-LHS example, and a count-matrix cbind example if
   counts landed); docs/design/multinomial.md "The surface" (formula ingestion
   defined); mark landed in docs/design/multi-forest-models.md and the TODO; this
   plan's landing notes. If a NEW exported topic is added (none expected - bart2
   is an existing topic), add the _pkgdown.yml entry + check_pkgdown; otherwise
   note pkgdown untouched. Files: man/bart2.Rd, docs/design/*, docs/plans/*, TODO,
   _pkgdown.yml (only if a new topic). Gate: full tinytest + R CMD check man
   (codoc/examples clean) + check_pkgdown if touched. Size: S.

## Verification

- Every commit: all three standing anchors identical (the engine is untouched, so
  this is a formality - but run them, a move means ingestion perturbed the design
  or seed).
- C1: the formula == matrix reproduction bit-for-bit at a fixed seed (the
  strongest gate - it proves the ingestion produces the identical labels/counts +
  coded design), predict on a data.frame newdata reproduces the fit's test
  channel, and (if counts) the cbind LHS == count matrix fit.
- dbarts.h unchanged -> no stan4bart lockstep, no rchk (no .Call change; ingestion
  is R-only).

## Open questions for VD

ADOPTED AS WORKING DEFAULTS (orchestrator, 2026-07-17, VD afk): both
recommendations below - explicit-only dispatch and the cbind LHS
idiom. Each is additive-relaxable later (auto-detect could be added,
a bare-matrix LHS could be accepted alongside cbind); VD may overrule
on return, before release.

- Q1 (auto-detect a multi-level factor LHS as multinomial?). RECOMMEND NO - keep
  family = "multinomial" EXPLICIT, never auto. A 2-level factor already resolves
  to probit under the binary path (dbarts.R:285-294), and silently routing a
  multi-level factor LHS to multinomial would collide with that idiom and
  surprise a user who expected an error or a different coding. AGAINST: some
  ecosystems infer multinomial from a factor response; but dbarts' "auto" is
  binary-vs-gaussian only, and multinomial is a deliberate opt-in. Gates C1's
  dispatch. RECOMMEND explicit-only.
- Q2 (count-matrix formula LHS syntax). RECOMMEND the cbind(c1, ..., cK) ~ x
  idiom, matching glm's binomial cbind(success, failure) form a user already
  knows; levels from the cbind column names. AGAINST: none material; a bare
  matrix symbol on the LHS (Y ~ x with Y a matrix in data) could also be
  accepted, but cbind is the conventional, self-documenting form. Only relevant
  once multinomial-counts has landed; gates C1's count branch (or a later
  additive commit under formula-first). RECOMMEND cbind.

## Landings

### C1 (2026-07-17)

As specced, both LHS branches together (multinomial-counts had already
landed). The refusal on a formula/dbartsData 'formula' argument is
split: a dbartsData object is still refused by name (out of scope,
unchanged message); a formula runs its own model.frame/terms extraction
(extractMultinomialFormulaData, R/bart.R) that pulls the response via
model.response(modelFrame) with no type argument, so a factor keeps its
levels and a cbind(...) matrix keeps its column names - the two then
fall through the SAME factor-vs-count-matrix dispatch and validation
the matrix interface already had, unduplicated.

One implementer finding worth keeping: the obvious design - code the
right-hand side here (makeModelMatrix on the term-labeled data frame,
matching dbartsData's own formula recipe) and install the resulting
matrix as the host sampler's matched-call 'formula' - does not work.
'factors = "categorical"', the default, produces a dbartsMixedMatrix
container even for all-numeric predictors (no factor columns), and the
matrix/vector top-level dispatch in dbartsData does not recognize that
container as a valid 'formula' argument (is.numeric/is.data.frame/
is.factor are all false for it) - only dbartsData's OWN
data-frame-as-x.train branch knows how to route it forward. The fix:
install the term-labeled predictor DATA FRAME (uncoded) as the matched
call's 'formula' instead, and let that existing branch call
makeModelMatrix itself. This means no coded matrix, and no terms/
xlevels object, is built or retained by this arc at all - the coding
and the newdata retention both ride attributes (term.labels,
factor.levels, drop) that makeModelMatrix already attaches, and
predict.bartMultinomial's existing validateXTest(newdata,
object$fit$data@x) call already reads them, unchanged.

Gates: full tinytest 2969/0 (+10, no regen - RNG-neutral); all three
standing anchors identical; the multinomial fixture identical vs
multinomial-equivalence-2bd34db.rds; air clean. tinytest additions
cover formula == matrix bit for bit for both the factor and cbind
count-matrix response forms, predict on a data.frame newdata against a
keepTrees formula fit, cbind column names as levels, an unused factor
level kept, character-LHS coercion, and the explicit-only default-family
refusal (unchanged from before this arc).

### C2 (2026-07-17)

docs/design/multinomial.md's surface section gained the formula-ingestion
bullet (the dispatch rule, the model.response no-coercion detail, the
data-frame-as-x.train routing, the reproduction gate); TODO's multi-forest
entry marks every multinomial follow-up landed. No man/bart2.Rd exists -
bart2 is documented under man/bart.Rd's family = "multinomial" paragraph,
updated there to describe the formula forms and drop the
matrix-interface-only sentence. No new exported topic, so no pkgdown
change. The Q1/Q2 defaults were adopted during VD's absence as annotated
above; neither was contradicted by implementation.
