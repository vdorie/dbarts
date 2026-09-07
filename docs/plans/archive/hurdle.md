# hurdle

agent: sonnet (R surface: ingestion, packaging, gates, docs); C2's combine +
  retransformation numerics to opus (the only report-time code a mistake biases).
rng: neutral - R-only composition, zero engine/bridge/dbarts.h/state code. Every
  existing equivalence anchor stays bitwise identical with NO re-record.
budget: ~450-700 lines R + gate scripts across 3 commits (design section 10):
  token+ingestion+wrapper; class+combine generics; gates+docs.

## Goal

family = "hurdle.lognormal" (alias "twopart") on dbarts()/bart2(): a semicontinuous
two-part fit composed R-side from a probit occupancy fit of 1{y > 0} over all n plus
a gaussian fit of log(y) over the y > 0 subset, glued at report time. A new
bartHurdle class + generics default to the natural (response) scale,
E[y | x] = Pr(y > 0 | x) * exp(f(x) + sigma^2 / 2). No engine, bridge, state-format,
or dbarts.h change; xbart/rbart_vi refuse the token. Spec: docs/design/hurdle.md,
section 13 authoritative (supersedes the in-body defaults of sections 3 and 6).

## Context - anchors, re-verified 2026-07-20

- Family token vectors + match.arg: dbarts() [R/dbarts.R:349-360](https://github.com/vdorie/dbarts/blob/09d617739e7f6d803c2f9248280c6230a2fa44cc/R/dbarts.R#L349-L360), [R/dbarts.R:370](https://github.com/vdorie/dbarts/blob/09d617739e7f6d803c2f9248280c6230a2fa44cc/R/dbarts.R#L370); bart2()
  [R/bart.R:412-424](https://github.com/vdorie/dbarts/blob/09d617739e7f6d803c2f9248280c6230a2fa44cc/R/bart.R#L412-L424), [R/bart.R:434](https://github.com/vdorie/dbarts/blob/09d617739e7f6d803c2f9248280c6230a2fa44cc/R/bart.R#L434). Refusal by omission: xbart [R/xbart.R:27](https://github.com/vdorie/dbarts/blob/09d617739e7f6d803c2f9248280c6230a2fa44cc/R/xbart.R#L27), [R/xbart.R:74](https://github.com/vdorie/dbarts/blob/09d617739e7f6d803c2f9248280c6230a2fa44cc/R/xbart.R#L74); rbart_vi
  [R/rbart.R:48](https://github.com/vdorie/dbarts/blob/09d617739e7f6d803c2f9248280c6230a2fa44cc/R/rbart.R#L48), [R/rbart.R:54](https://github.com/vdorie/dbarts/blob/09d617739e7f6d803c2f9248280c6230a2fa44cc/R/rbart.R#L54) (the nbinom precedent - their match.arg vectors ARE the refusal).
- R-composition precedent to mirror (discrete-time hazard): ingestion splits,
  remaps the token, refuses subset/test [R/dbarts.R:399-459](https://github.com/vdorie/dbarts/blob/09d617739e7f6d803c2f9248280c6230a2fa44cc/R/dbarts.R#L399-L459); parks its marker on a
  control attr [R/dbarts.R:794-795](https://github.com/vdorie/dbarts/blob/09d617739e7f6d803c2f9248280c6230a2fa44cc/R/dbarts.R#L794-L795), which packaging reads into $periods
  [R/bart.R:274-280](https://github.com/vdorie/dbarts/blob/09d617739e7f6d803c2f9248280c6230a2fa44cc/R/bart.R#L274-L280). bart2 family branch shape [R/bart.R:728-755](https://github.com/vdorie/dbarts/blob/09d617739e7f6d803c2f9248280c6230a2fa44cc/R/bart.R#L728-L755) (nbinom) ->
  bart2Negbin [R/bart.R:1553](https://github.com/vdorie/dbarts/blob/09d617739e7f6d803c2f9248280c6230a2fa44cc/R/bart.R#L1553); packageNegbinResults + `class(result) <- "bartNegbin"`
  [R/bart.R:1681](https://github.com/vdorie/dbarts/blob/09d617739e7f6d803c2f9248280c6230a2fa44cc/R/bart.R#L1681), [R/bart.R:1746](https://github.com/vdorie/dbarts/blob/09d617739e7f6d803c2f9248280c6230a2fa44cc/R/bart.R#L1746). Standard fit: packageBartResults -> class "bart"
  [R/bart.R:134](https://github.com/vdorie/dbarts/blob/09d617739e7f6d803c2f9248280c6230a2fa44cc/R/bart.R#L134), [R/bart.R:300](https://github.com/vdorie/dbarts/blob/09d617739e7f6d803c2f9248280c6230a2fa44cc/R/bart.R#L300); bart2 std path calls it [R/bart.R:858](https://github.com/vdorie/dbarts/blob/09d617739e7f6d803c2f9248280c6230a2fa44cc/R/bart.R#L858).
- Fit-class generics idiom: bartNegbin extract/fitted/residuals/predict/print
  [R/generics.R:760](https://github.com/vdorie/dbarts/blob/09d617739e7f6d803c2f9248280c6230a2fa44cc/R/generics.R#L760), [R/generics.R:788](https://github.com/vdorie/dbarts/blob/09d617739e7f6d803c2f9248280c6230a2fa44cc/R/generics.R#L788), [R/generics.R:796](https://github.com/vdorie/dbarts/blob/09d617739e7f6d803c2f9248280c6230a2fa44cc/R/generics.R#L796), [R/generics.R:806](https://github.com/vdorie/dbarts/blob/09d617739e7f6d803c2f9248280c6230a2fa44cc/R/generics.R#L806), [R/generics.R:856](https://github.com/vdorie/dbarts/blob/09d617739e7f6d803c2f9248280c6230a2fa44cc/R/generics.R#L856); S3 registrations NAMESPACE:52-56. The combine is
  NEW code keyed on the bartHurdle class; each component keeps its own $family
  ("probit"/"gaussian") so $family-dispatched helpers (probabilityFromLatents
  [R/generics.R:13-19](https://github.com/vdorie/dbarts/blob/09d617739e7f6d803c2f9248280c6230a2fa44cc/R/generics.R#L13-L19), pointwiseLogLikelihood [R/generics.R:35-87](https://github.com/vdorie/dbarts/blob/09d617739e7f6d803c2f9248280c6230a2fa44cc/R/generics.R#L35-L87)) stay correct per component.
- predict.bart type aliases "response"->"ev", "link"->"bart" [R/generics.R:231-235](https://github.com/vdorie/dbarts/blob/09d617739e7f6d803c2f9248280c6230a2fa44cc/R/generics.R#L231-L235);
  ppd via sampleFromPPD [R/generics.R:203-204](https://github.com/vdorie/dbarts/blob/09d617739e7f6d803c2f9248280c6230a2fa44cc/R/generics.R#L203-L204); heteroscedastic per-obs s(x) rides back as a yhat
  attribute [R/generics.R:183-193](https://github.com/vdorie/dbarts/blob/09d617739e7f6d803c2f9248280c6230a2fa44cc/R/generics.R#L183-L193), [R/generics.R:217-219](https://github.com/vdorie/dbarts/blob/09d617739e7f6d803c2f9248280c6230a2fa44cc/R/generics.R#L217-L219) (the per-observation sigma the
  retransformation must consume when variance = ~x is set on the positive part).
- y >= 0 validation precedent (nbinom count check) [R/dbarts.R:552-561](https://github.com/vdorie/dbarts/blob/09d617739e7f6d803c2f9248280c6230a2fa44cc/R/dbarts.R#L552-L561). seed ->
  rngSeed [R/bart.R:401](https://github.com/vdorie/dbarts/blob/09d617739e7f6d803c2f9248280c6230a2fa44cc/R/bart.R#L401), [R/bart.R:498-499](https://github.com/vdorie/dbarts/blob/09d617739e7f6d803c2f9248280c6230a2fa44cc/R/bart.R#L498-L499) (derive two independent seeds here).
- Reduction-gate SHAPE (not its tautology, hardening a): benchmarks/R/
  hazard-reduction.R (compareLink; markerOnly [R/bart.R:85](https://github.com/vdorie/dbarts/blob/09d617739e7f6d803c2f9248280c6230a2fa44cc/R/bart.R#L85)). Equivalence harness: scenario
  list [benchmarks/R/equivalence.R:60-561](https://github.com/vdorie/dbarts/blob/09d617739e7f6d803c2f9248280c6230a2fa44cc/benchmarks/R/equivalence.R#L60-L561), hazard scenario [benchmarks/R/equivalence.R:524-559](https://github.com/vdorie/dbarts/blob/09d617739e7f6d803c2f9248280c6230a2fa44cc/benchmarks/R/equivalence.R#L524-L559), fitViaHazard
  [benchmarks/R/equivalence.R:708-732](https://github.com/vdorie/dbarts/blob/09d617739e7f6d803c2f9248280c6230a2fa44cc/benchmarks/R/equivalence.R#L708-L732), fitSummaries dispatch [benchmarks/R/equivalence.R:801-814](https://github.com/vdorie/dbarts/blob/09d617739e7f6d803c2f9248280c6230a2fa44cc/benchmarks/R/equivalence.R#L801-L814), new-scenario "skipped/uncovered" policy
  [benchmarks/R/equivalence.R:428](https://github.com/vdorie/dbarts/blob/09d617739e7f6d803c2f9248280c6230a2fa44cc/benchmarks/R/equivalence.R#L428), [benchmarks/R/equivalence.R:1127-1152](https://github.com/vdorie/dbarts/blob/09d617739e7f6d803c2f9248280c6230a2fa44cc/benchmarks/R/equivalence.R#L1127-L1152). pkgdown reference sections [_pkgdown.yml:9-43](https://github.com/vdorie/dbarts/blob/09d617739e7f6d803c2f9248280c6230a2fa44cc/_pkgdown.yml#L9-L43) (no per-family Rd
  exists today - man/ has none for bartNegbin/bartOrdinal). air.toml present.

## Constraints

- Engine byte-neutral: no src/, inst/include/dbarts/dbarts.h, *.in, or state-format
  touch; no tests/cpp impact. Every existing equivalence anchor bitwise identical.
- Canonical token "hurdle.lognormal"; "twopart" is an accepted alias that resolves
  and PRINTS as hurdle.lognormal (design section 13 NAMING).
- INDEPENDENT deterministic per-component seeds derived from the user seed
  (hardening b): a shared seed correlates the two chains and biases the combined
  credible interval (not the mean).
- The positive fit receives the full-n x as x.test (hardening c) so in-sample
  fitted()/extract() carries E[y | y > 0, x] at the zero rows it never trained on.
  Both component fits keep trees when the hurdle fit does (predict replays both).
- Reporting default NATURAL scale via posterior-predictive Monte Carlo (section 13):
  E[y | y > 0, x] = exp(f + sigma^2 / 2) consuming PER-OBSERVATION sigma from the
  positive fit (single sigma^2 by default, s(x)^2 when variance = ~x), then
  E[y | x] = pi * E[y | y > 0]. Opt-ins type = "link"/"log"; type = "ppd" is bimodal
  (draw Bernoulli(pi), then the lognormal, hardening d). Duan smearing is a door.
- The wrapper is family-agnostic on the positive part internally (count door,
  section 9); v1 wires only gaussian-on-log-y.
- Out of scope: engine hurdle, zero-inflation/Heckman, count/gamma positive parts,
  and a SHARED variable-selection prior across parts (document it foreclosed by the
  R route, hardening e).

## Steps

1. C1 - token + ingestion + wrapper, sonnet. Add "hurdle.lognormal"/"twopart" to the
   dbarts and bart2 family vectors ([R/dbarts.R:349](https://github.com/vdorie/dbarts/blob/09d617739e7f6d803c2f9248280c6230a2fa44cc/R/dbarts.R#L349), [R/bart.R:412](https://github.com/vdorie/dbarts/blob/09d617739e7f6d803c2f9248280c6230a2fa44cc/R/bart.R#L412)), resolve the alias,
   validate y >= 0 ([R/dbarts.R:552](https://github.com/vdorie/dbarts/blob/09d617739e7f6d803c2f9248280c6230a2fa44cc/R/dbarts.R#L552) precedent), refuse weights/subset/test as hazard
   does. A bart2Hurdle branch ([R/bart.R:728](https://github.com/vdorie/dbarts/blob/09d617739e7f6d803c2f9248280c6230a2fa44cc/R/bart.R#L728) shape) fits the occupancy probit over all
   n and the gaussian over log(y[S]) with x.test = full-n x, at two derived seeds
   ([R/bart.R:498](https://github.com/vdorie/dbarts/blob/09d617739e7f6d803c2f9248280c6230a2fa44cc/R/bart.R#L498)); package a bartHurdle holding both component fits + a variant marker
   + minimal print. Gate: R CMD INSTALL; benchmarks/R/hurdle-reduction.R - each
   internal fit equals a standalone probit/gaussian fit at the SAME derived seed,
   bitwise (markerOnly), the hazard-reduction shape as a sanity FLOOR, not the
   correctness argument. ~150-250 lines.
2. C2 - class + combine/retransform generics, opus (the risk-bearing commit).
   extract/fitted/predict/residuals/print.bartHurdle ([R/generics.R:760](https://github.com/vdorie/dbarts/blob/09d617739e7f6d803c2f9248280c6230a2fa44cc/R/generics.R#L760) idiom;
   NAMESPACE:52 registrations). The combine layer: type = "ev"/"response" natural-
   scale E[y | x] = pi(x) * exp(f(x) + sigma^2 / 2) per draw, reading per-observation
   sigma via the positive fit's attribute path ([R/generics.R:183-193](https://github.com/vdorie/dbarts/blob/09d617739e7f6d803c2f9248280c6230a2fa44cc/R/generics.R#L183-L193)); "link"/"log"
   the positive linear predictor; "prob" pi(x) through the correct link; "ppd"
   bimodal. predict replays both saved forests at newdata and combines. Gate:
   R CMD INSTALL; tinytest - the ANALYTIC combine/retransform oracle (hand-set
   pi/mu/sigma checked in closed form, the only genuinely new code, hardening a) +
   predict-on-newdata + save/load; a recovery smoke recovering pi(x), E[y | y > 0, x],
   and combined E[y]. ~150-250 lines.
3. C3 - gates + docs, sonnet. New "hurdle" scenario in benchmarks/R/equivalence.R
   ([R/generics.R:524](https://github.com/vdorie/dbarts/blob/09d617739e7f6d803c2f9248280c6230a2fa44cc/R/generics.R#L524) shape) recording the occupancy channel + the positive channel + the combined
   predict; fitViaHurdle + fitSummaries dispatch ([R/generics.R:801](https://github.com/vdorie/dbarts/blob/09d617739e7f6d803c2f9248280c6230a2fa44cc/R/generics.R#L801)). tinytest: family routing, the
   "twopart" alias printing as hurdle.lognormal, y >= 0 validation, xbart/rbart_vi
   refusal. Document the foreclosed shared variable-selection prior + the smearing
   door. Gate: equivalence compare vs the current baseline - ALL existing scenarios
   "identical draws (same RNG stream)", hurdle uncovered/skipped (no re-record; the
   anchor re-records at landing, the nbinom/hazard trail); air format --check .;
   check_pkgdown only if a Rd topic lands (none expected); full tinytest. ~150-250.

## Verification

Per commit: R CMD INSTALL . (no --preclean - no header/Makevars/config change);
tinytest::test_package("dbarts") green. C1: hurdle-reduction.R prints BITWISE
IDENTICAL for both components. C2: the analytic oracle, predict-on-newdata,
save/load, and recovery tinytest pass. C3: equivalence compare reports every
pre-existing scenario "identical draws (same RNG stream)" and hurdle as uncovered;
air format --check . clean. No tests/cpp run needed (engine untouched); confirm the
git diff touches no src/, no inst/include/dbarts/dbarts.h, no *.in, and no state
format.

## Landing

Landed 2026-07-20 as three R-only commits; the engine binary is untouched and
every pre-existing equivalence scenario draws identically.

1. C1 901581e - the family token (hurdle.lognormal, "twopart" resolving to it)
   and the bart2Hurdle wrapper: an occupancy probit on 1{y>0} over all n plus a
   gaussian on log(y) over the y>0 subset, at two independently derived seeds,
   the positive fit given the full-n x as its x.test. dbarts() refuses the token
   (one sampler cannot compose two). bartHurdle class + minimal print. Gate:
   benchmarks/R/hurdle-reduction.R (each component reduces bitwise to a
   standalone fit - the sanity floor).
   Follow-up a70fec6: the token was added to the formals but not the Rd usage,
   so R CMD check flagged a codoc mismatch (CI runs error_on=warning); fixed in
   bart.Rd/dbarts.Rd. LESSON: R-touching commits need R CMD check locally, not
   just INSTALL + tinytest.
2. C2 6b11487 - extract/fitted/predict/residuals.bartHurdle: the per-draw,
   natural-scale, heteroscedasticity-aware retransformation
   E[y|x] = pi*exp(f + sigma^2/2) by posterior-predictive Monte Carlo (not a
   plug-in), draw-aligned across the two independent fits, reading per-obs sigma
   so a variance=~x positive part works unchanged; prob/link/log types and the
   bimodal ppd; predict replays both forests. Gate: an analytic combine/
   retransform oracle (hand-set pi/f/sigma in closed form) + predict-on-newdata
   + save/load + a recovery smoke.
3. C3 7903855 - the equivalence hurdle scenario, surface tinytest (routing, the
   twopart alias, the non-negative / require-a-zero / require-a-positive
   validation, the dbarts and xbart/rbart_vi refusals), the hurdle.lognormal
   family Rd paragraph in bart.Rd/dbarts.Rd, and a NEWS 1.0-0 bullet.

Gates: R CMD check 0 errors / 0 warnings (only the two pre-existing show/rnbinom
NOTES); the equivalence trio bitwise (all 26 pre-existing scenarios identical,
hurdle new); suite 3359/0; air clean. The equivalence baseline re-recorded to
equivalence-7903855.rds (27 scenarios; self-reproduces 27/27 under
--strict-coverage; f494156 demoted), MANIFEST and equivalence.yaml re-pinned.

Design-vs-code reconciliation (section 11): forest-combiner.md's hurdle bullet
framed hurdle as the engine model that breaks Chain's single response_ - correct
about what an ENGINE hurdle would need, but hurdle landed R-side so that
invariant break stays unbuilt (noted there). The "second two-leaf-type consumer"
heteroscedastic deferred to hurdle does not arrive via hurdle; the engine
two-response generalization waits for a genuinely coupled model (zero-inflation
or Heckman selection, section 9). Doors: count hurdle (needs a zero-truncated
count family first), gamma positive part, logistic occupancy, grouped hurdle,
Duan smearing, and the coupled cousins.
