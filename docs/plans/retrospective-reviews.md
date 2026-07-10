# retrospective-reviews

agent: fable orchestrates; per-review plans authored as each starts
rng: neutral unless a review surfaces a fix (fixes ride their own
     items under the standard gates)
budget: program of six reviews; each gets its own plan file and lands
        its findings as a design note or plan resolution

Greenlit by VD 2026-07-08 ("All 6, in any order you please"), with
binding execution constraints: reviews run STRICTLY SEQUENTIALLY - one
review at a time, modest fan-out inside a review, sub-agent spawning
by readers forbidden (the 2026-07-07 session-limit exhaustion came
from a reader fan-out of fan-outs) - and all delegated work runs on
Opus (engine/numerics/derivations) or Sonnet (mechanical/R/docs),
with meta review by the orchestrating Fable session (VD, 2026-07-08).

## The program

The bartcore branch closes ~15 years of development; this is the
pause-and-reflect pass. Six reviews, in execution order:

1. correctness-audit: adversarial re-derivation of every acceptance
   ratio and conjugate update (moves.hpp, model.hpp, chain.hpp -
   birth/death incl. node selection, change/swap cancellations,
   Polya-Gamma, chi-k, DART Dirichlet, grouped intercepts, BCF glue,
   weights per family); multiple independent derivers per formula
   block prompted to find errors, cross-checked. Highest value: the
   exact-posterior gates verify few configurations; equivalence
   verifies only sameness.
2. interface-review: naming/default coherence across bart/bart2/
   dbarts/rbart_vi/xbart, generics consistency, error-message
   quality, new-user walkthrough on shipped docs only. Deliverable:
   fix-before-release vs 2.0-wishlist, taste calls to VD. Runs early
   because 1.0-0 freezes the R surface at submission.
3. gate-blindspot-audit: coverage measurement, poison-check sweep
   (break each kernel deliberately, confirm a gate fails), and a
   feature-combination x gate matrix. Feeds item 4.
4. sbc-calibration: simulation-based calibration (Talts et al.) across
   the feature matrix (families x weights x DART x grouped x BCF x
   missingness), prioritized by item 3's uncovered combos; candidate
   new standing gate class.
5. roadmap-survey: 15 years of BART literature + competing packages
   (BART, bartMachine, SoftBart, stochtree/bartz) vs the shipped
   feature set; ranks the existing consumer-gated TODO items by
   demand evidence into a roadmap note. Decisions stay VD's.
6. architecture-critique + numerical-robustness (combined item if
   both stay small): fresh-eyes stress of the facade/concept
   decomposition against the five likely-next features; extreme-scale
   and log-space robustness, GP nugget conditioning, long-run drift.

Queued behind the two greenlit experiments (parallel-falsifiers, then
grow-from-root). Each review's plan file is authored against the code
when it starts; findings land as design notes or plan resolutions,
and any code fix becomes its own gated item.

## Status

- 2026-07-08: program recorded; parallel-falsifiers in flight.
- 2026-07-08: review 1 (correctness-audit) COMPLETE - all seven
  blocks adjudicated; three defects found and fixed (change-move
  balance, chi hyperprior df, zero-weight sigma df), two latent
  items filed and one already landed (chi-k-empty-leaf-count) with
  bcf-testfits-guard remaining; full records in
  correctness-audit.md.
- 2026-07-08: review 2 (interface-review) started; plan authored
  (interface-review.md).
- 2026-07-08: review 2 COMPLETE - four readers, findings verified
  by orchestrator (one refuted), F1-F11/D1-D6 landed 8744e77, all
  eleven VD taste decisions landed c1ea79a (both draw-neutral,
  equivalence identical). Full record in interface-review.md.
- 2026-07-09: review 3 (gate-blindspot-audit) COMPLETE - coverage
  matrix + 16-poison sweep. One TRUE survivor (BCF a-glue prior
  precision, zero gates) plus 5 single-point-of-failure gaps and a
  --preclean stale-.so footgun; deliverables filed as
  gate-hardening-1.0 and the SBC target list. Full record in
  gate-blindspot-audit.md. (Orchestrator downgraded Fable->Opus
  mid-review at VD's Fable-limit; no work lost.)
- 2026-07-09: gate-hardening-1.0 LANDED (0e9ccca; extended
  21-scenario anchor equivalence-0e9ccca.rds recorded 3ef0308).
- 2026-07-09: review 4 (sbc-calibration) COMPLETE and landed
  (25d3dc1): 18 configs across three tiers, ONE calibration defect
  (BCF glue-on sigma biased high, localized to the glue-updating
  path by a fixed-glue control; filed as
  bcf-glue-sigma-calibration with a mixing-vs-defect experiment
  plan). All review-3 priority-one NA/imputation paths calibrate;
  DART 1e-300 floor does not bias selection; the poison-13 a-glue
  channel calibrates (gate gap, not a latent bug - now gated by
  bcf-exact-weak.R). Residual SBC gaps recorded: sampled-k chi
  hyperprior has no SBC path for any leaf type; DART alpha held
  fixed; half-Cauchy grouped tau SBC-intractable through the slice
  sampler. Routed to review 6 (numerical robustness): the GP
  predict-vs-stored-fits jitter (~2e-3 at training rows), the tau
  slice-sampler hang under extreme ranef weights, the half-Cauchy
  tail behavior. Soft call for VD: the gaussian SBC baseline
  (~52s at R=200) is a candidate standing CI calibration gate.
- 2026-07-09: review 6 (architecture-critique + numerical-
  robustness) started BEFORE review 5 (roadmap survey needs
  sustained web access; the machine is on plane WiFi - VD
  authorized any order). Plan: architecture-numerical-review.md.
