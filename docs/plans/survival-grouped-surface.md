# survival-grouped-surface

agent: opus
rng: neutral - existing gaussian/binary rbart draws bit-identical (the
  aft branch fires only for a Surv LHS or family = "aft", and
  extractSurvivalResponse is pure); grouped aft is newly reachable, no
  existing snapshot moves. A grouped-aft equivalence scenario is ADDED,
  so the anchor re-records at landing (all existing scenarios must stay
  identical against the old anchor first).
budget: ~200 lines R + Rd + tinytest; ZERO C/C++ expected.

## Goal

Grouped (random-intercept) AFT survival fits reachable from R so
riAFTBART can migrate: rbart_vi grows a family = "aft" path and
survivalProbabilities.rbart replaces its error with real
S(t | x_i, group_i) draws that include the drawn intercepts. The
engine composition GroupedResponse(AFTResponse) is already built and
gated; this is R surface only.

## Context (verified against the code)

- Engine/bridge DONE: createHolder applies the group and survival
  attributes unconditionally (src/R_interface_bartcore.cpp ~1209 and
  ~1212); grouped-aft is gated at inst/tinytest/test-aft.R (the
  grouped block returning ranef); docs/design/survival.md's
  composition-check section. Engine delta expected ZERO; verify only.
- rbart_vi is formula-only, has no family argument. Its in-core fast
  path (built-in tau prior, no callback) pre-sets
  attr(control, "bartcore.groups") and calls dbarts(); the
  custom-prior/callback path runs an R Gibbs loop with a binary
  latent refresh via getLatents.
- Guards landed at f0efc03: dbarts() refuses family "aft" on the
  indirect path (R/dbarts.R ~224); dbartsData's formula path refuses
  a Surv response (R/data.R ~418). Both stay intact for dbarts/bart2.
- survivalProbabilities.bart computes draws from extract(type="bart")
  plus sigma; the rbart ev channel (f + b) is on the log scale and is
  NOT probability-transformed when sigma is present; predict.rbart
  draws N(0, tau) intercepts for unseen groups.

## Decisions (VD 2026-07-10, all approved)

- Ingestion shape A: Surv / two-column (time, status) response through
  the rbart_vi FORMULA LHS, lifting the Surv refusal ONLY on the
  rbart_vi path (rbart_vi has no matrix interface for the refusal to
  redirect to). A Surv LHS auto-selects aft from family "auto"; a bare
  two-column response needs explicit family = "aft".
- Family set minimal: family = c("auto", "gaussian", "aft"). The
  logistic asymmetry with bart2 stays a separate item.
- WIRE the custom-prior/callback aft path (mirror the binary
  getLatents latent refresh in the R loop, a few lines), do not
  refuse it.
- Add a grouped-aft equivalence scenario; the anchor re-records at
  landing (orchestrator step).

## Steps

1. rbart_vi: add the family argument. On a Surv LHS or
   family = "aft", extract (log.time, status) BEFORE dbartsData sees
   the Surv (parallel to the group.by pre-processing), rewrite the
   response to log-time so the formula path passes, thread status via
   attr(control, "bartcore.survival"), refuse weights and subset,
   pass family = "aft" into the dbarts() call.
2. dbarts() guard: permit family = "aft" on the indirect path when
   attr(control, "bartcore.survival") is already set (the internal
   rbart channel, ~3 lines); the public refusal stays otherwise.
3. Custom-prior/callback path: refresh the working response for aft
   as the binary branch does (getLatents mirror); thread family in.
4. survivalProbabilities.rbart: mirror the bart method, sourcing the
   linear predictor from extract(type = "ev", sample = "train") /
   predict(type = "ev", group.by = ) so the intercepts are included;
   sigma from object$sigma in the uncombined convention; args
   times / newdata / group.by / combineChains; gate on
   object$family == "aft"; unseen-group behavior inherited from
   predict.rbart.
5. Docs: rbart.Rd (family, Surv/two-column ingestion, log-time scale
   language mirroring bart.Rd's) and survivalProbabilities.Rd (rbart
   method arguments, drop the not-available note). NAMESPACE already
   exports both methods.
6. Equivalence: add a grouped-aft scenario to benchmarks/R/
   equivalence.R (seeded, censored, in-core fast path); confirm every
   EXISTING scenario stays "identical draws" against
   equivalence-de67cbb.rds (the new scenario reports as skipped).

## Verification

- Reduction gate (primary, free): rbart_vi family = "aft" on
  all-status-1 data == rbart_vi gaussian on log T, same seed, BITWISE
  (ranef, tau, sigma, train fits) - the grouped analogue of the
  ungrouped aft reduction gate.
- tinytest: survivalProbabilities.rbart dims (draws x times x obs;
  chain margin under combineChains = FALSE), values in [0, 1] and
  monotone in t, a high-intercept group shows uniformly higher S(t)
  (proves intercepts are included), newdata + group.by with an unseen
  group, non-aft refusal; grouped recovery-under-censoring smoke;
  custom-prior aft smoke (the wired R-loop path runs and returns
  finite draws); weights/subset refusal.
- Full tinytest; equivalence existing scenarios identical; air +
  lintr; pkgdown clean; tests/cpp untouched.
