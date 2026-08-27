# ordinal-outcomes

agent: opus
rng: posterior-changing (new family)
budget: design note + ~400 lines

## Goal

Ordered categorical responses fit via cumulative probit: the existing
truncated-normal latent machinery plus a cutpoint block, surfaced as
family = "ordinal" (or ordered-factor response auto-dispatch - the
note decides).

## Context

- Probit latents already exist (ProbitResponse,
  src/bartcore/model.hpp); ordinal generalizes the truncation from
  {(-inf,0), (0,inf)} to per-category intervals (gamma_{k-1},
  gamma_k].
- The cutpoint update is the known hard part: Albert-Chib's Gibbs
  cutpoints mix badly; the note should evaluate the standard fixes
  (Cowles' MH block update, or fixing K-2 cutpoints via a latent-scale
  identification and sampling a scale). Identification also interacts
  with the response scaling and node.scale calibration.
- Ingestion: ordered factors currently fit as ordinal predictors;
  as a response they error - dbartsData response typing gains a path.

## Constraints

- Exact-posterior gate (single tree, few categories: enumeration +
  quadrature over latents and cutpoints on a grid) before shipping.
- rbart_vi/grouped decorator composition is desirable (ordinal
  mixed models) but not v1; refuse cleanly.
- Out of scope: unordered multinomial (multi-forest-models);
  ordinal-scale loss functions in xbart (follow-up).

## Steps

1. Design note (docs/design/ordinal.md): identification scheme,
   cutpoint sampler choice with a small mixing study, surface and
   prediction semantics (category probabilities vs latent scale -
   probabilityFromLatents generalizes).
2. OrdinalResponse + cutpoint block; ingestion typing; R surface,
   predict/fitted transforms, Rd.
3. Gates: exact-posterior, component tests (truncation intervals,
   cutpoint update), recovery, mutation smoke (setResponse re-draw
   semantics defined).

## Verification

- Exact-posterior gate to MC error; component tests; full tinytest;
  bench on existing families unchanged.

## Landings

Design 8777bfe (2026-07-18). docs/design/ordinal.md, written Opus-first
then independently critiqued (ACCEPT WITH CHANGES, seven findings
applied - the surface fork was reframed after the critique showed
bart2's auto path silently fits unordered multinomial on ordered
factors, and the designer's recommendation flipped to auto-dispatch).
All decisions resolved by VD: scheme A identification (sigma = 1,
gamma_1 = 0, K-2 free cutpoints; a pure location anchor, the
Albert-Chib/MCMCoprobit convention), Cowles-style marginal MH (one
cutpoint at a time, plain RW proposal with out-of-bounds rejection,
fixed count-derived scale - prototype mixing study in
benchmarks/R/ordinal-cutpoint-mixing.R: 15-54x ESS over Albert-Chib
Gibbs, n-stable), auto-dispatch on is.ordered() with announcement
(class-based detection, never level names). Shipped constants recorded
in the note: log-gap prior N(0, 1.5^2), proposal constant 2.5.

C1 1afe202 (2026-07-18). OrdinalResponse in model.hpp (marginal
cutpoint MH then doubly-truncated latent redraw; K = 2 bitwise probit
via one-sided-primitive routing and matching NaN fallback, locked by a
stream-identity component test); ext_rng_simulateTruncatedNormalScale1
(inverse-CDF, Robert rejection fallback) in external/random; three
RNG-insulated component tests. All anchors bitwise.

C2 227f46a (2026-07-18). ResponseFamily::ordinal + the
carriesCutpoints virtual quartet (the residualDf trio's vector
analog); SamplerOptions.numCategories; construction in chain.hpp;
by-name "cutpoints" state slot (conditional write, absence-tolerated
decode, stateIsValid length check; old states load unchanged); bridge
third response-shape channel via the bartcore.n.categories control
attribute, with weight, grouped, and off-shape family refusals. No
dbarts.h diff. State round-trip component test. All anchors bitwise.

C3 b929784 (2026-07-18). R surface: family = "ordinal" on
dbarts/bart2, auto-dispatch on 3+-level ordered factors with
announcement (2-level ordered stays probit); ordered factors split
out of detectAutoMultinomial - bart2's silent auto behavior (fit
unordered multinomial on an ordered response, order discarded) is
replaced, the one deliberate behavior change; bartOrdinal fit class
(levels, cutpoints, n x K probability draws, latent draws, type =
"link"); rbart_vi/xbart/bart refuse informatively. Gates:
benchmarks/R/ordinal-exact.R z-free quadrature reference (max
category-prob gap 0.0006 vs tol 0.012, gamma_2 gap 0.0011 vs 0.040);
suite grew 3098 -> 3168; new ordinal equivalence scenario, baseline
re-recorded as equivalence-227f46a.rds (named for the recording
engine HEAD; C3 is R-only) with the neutrality trail: all 23
pre-existing scenarios bitwise vs 31dc05a; equivalence.yaml
re-pinned. Deviation accepted: no engine per-sample cutpoint results
channel in v1, so the bart2 driver keeps one sample per run() call
and reads cutpoints from state - draw-neutral (locked by the
equivalence scenario), retired by the ordinal-results-channel TODO
entry. One implementer API failure mid-run, resumed cleanly from
transcript. ARC CLOSED; recorded doors: grouped ordinal, robit-style
links, xbart ordinal loss, dbarts.h exposure.

Follow-up 4727d1e (2026-07-18). The results channel: "cutpoints"
appended as a ninth run-result slot only when the family carries
cutpoints (Results.cutpoints + numCutpoints mirror the
numVariableCountForests pattern; per-chain stride in the sampler;
facade virtual; non-ordinal result lists byte-identical at 8 slots);
the bart2 driver collapsed to a single run(n.burn, n.samples) reading
the channel - the per-sample state reads are gone. Draw-neutrality
VERIFIED: the ordinal equivalence scenario, recorded under the old
one-sample-per-run driving, reproduces bitwise under the single run
(24/24 identical draws), proving the drivings equivalent. All seven
gates green; exact gate unchanged (0.0006 vs 0.012); dbarts.h
untouched. The C3 deviation is fully retired.
