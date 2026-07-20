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

- Family token vectors + match.arg: dbarts() R/dbarts.R:349-360,370; bart2()
  R/bart.R:412-424,434. Refusal by omission: xbart R/xbart.R:27,74; rbart_vi
  R/rbart.R:48,54 (the nbinom precedent - their match.arg vectors ARE the refusal).
- R-composition precedent to mirror (discrete-time hazard): ingestion splits,
  remaps the token, refuses subset/test R/dbarts.R:399-459; parks its marker on a
  control attr R/dbarts.R:794-795, which packaging reads into $periods
  R/bart.R:274-280. bart2 family branch shape R/bart.R:728-755 (nbinom) ->
  bart2Negbin R/bart.R:1553; packageNegbinResults + `class(result) <- "bartNegbin"`
  R/bart.R:1681,1746. Standard fit: packageBartResults -> class "bart"
  R/bart.R:134,300; bart2 std path calls it R/bart.R:858.
- Fit-class generics idiom: bartNegbin extract/fitted/residuals/predict/print
  R/generics.R:760,788,796,806,856; S3 registrations NAMESPACE:52-56. The combine is
  NEW code keyed on the bartHurdle class; each component keeps its own $family
  ("probit"/"gaussian") so $family-dispatched helpers (probabilityFromLatents
  R/generics.R:13-19, pointwiseLogLikelihood :35-87) stay correct per component.
- predict.bart type aliases "response"->"ev", "link"->"bart" R/generics.R:231-235;
  ppd via sampleFromPPD :203-204; heteroscedastic per-obs s(x) rides back as a yhat
  attribute R/generics.R:183-193,217-219 (the per-observation sigma the
  retransformation must consume when variance = ~x is set on the positive part).
- y >= 0 validation precedent (nbinom count check) R/dbarts.R:552-561. seed ->
  rngSeed R/bart.R:401,498-499 (derive two independent seeds here).
- Reduction-gate SHAPE (not its tautology, hardening a): benchmarks/R/
  hazard-reduction.R (compareLink; markerOnly :85). Equivalence harness: scenario
  list benchmarks/R/equivalence.R:60-561, hazard scenario :524-559, fitViaHazard
  :708-732, fitSummaries dispatch :801-814, new-scenario "skipped/uncovered" policy
  :428,1127-1152. pkgdown reference sections _pkgdown.yml:9-44 (no per-family Rd
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
   dbarts and bart2 family vectors (R/dbarts.R:349, R/bart.R:412), resolve the alias,
   validate y >= 0 (R/dbarts.R:552 precedent), refuse weights/subset/test as hazard
   does. A bart2Hurdle branch (R/bart.R:728 shape) fits the occupancy probit over all
   n and the gaussian over log(y[S]) with x.test = full-n x, at two derived seeds
   (R/bart.R:498); package a bartHurdle holding both component fits + a variant marker
   + minimal print. Gate: R CMD INSTALL; benchmarks/R/hurdle-reduction.R - each
   internal fit equals a standalone probit/gaussian fit at the SAME derived seed,
   bitwise (markerOnly), the hazard-reduction shape as a sanity FLOOR, not the
   correctness argument. ~150-250 lines.
2. C2 - class + combine/retransform generics, opus (the risk-bearing commit).
   extract/fitted/predict/residuals/print.bartHurdle (R/generics.R:760 idiom;
   NAMESPACE:52 registrations). The combine layer: type = "ev"/"response" natural-
   scale E[y | x] = pi(x) * exp(f(x) + sigma^2 / 2) per draw, reading per-observation
   sigma via the positive fit's attribute path (R/generics.R:183-193); "link"/"log"
   the positive linear predictor; "prob" pi(x) through the correct link; "ppd"
   bimodal. predict replays both saved forests at newdata and combines. Gate:
   R CMD INSTALL; tinytest - the ANALYTIC combine/retransform oracle (hand-set
   pi/mu/sigma checked in closed form, the only genuinely new code, hardening a) +
   predict-on-newdata + save/load; a recovery smoke recovering pi(x), E[y | y > 0, x],
   and combined E[y]. ~150-250 lines.
3. C3 - gates + docs, sonnet. New "hurdle" scenario in benchmarks/R/equivalence.R
   (:524 shape) recording the occupancy channel + the positive channel + the combined
   predict; fitViaHurdle + fitSummaries dispatch (:801). tinytest: family routing, the
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

(stub - fill at landing: commit SHAs, gate transcripts, any design-vs-code
deviations, and the section-11 reconciliation edits the orchestrator owns.)
