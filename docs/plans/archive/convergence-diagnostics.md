# convergence-diagnostics

agent: sonnet
rng: neutral (accessors over existing draws)
budget: ~250 lines R + Rd

## Goal

Multi-chain fits answer "did my chains converge" without the user
hand-rolling array surgery: draws come out in posterior-package
shape, and summary() reports R-hat/ESS for the scalar parameters.

## Context

- Chain-dimensioned draws already exist on fits (sigma, k, varcount,
  varprobs, yhat.*) but as raw arrays with dbarts' native ordering;
  no as_draws support, no diagnostics anywhere.
- The posterior package is the ecosystem standard (rstan, brms, cmdstanr
  interoperate through it); it must be Suggests, not Depends.

## Constraints

- No hard dependency: as_draws methods registered conditionally;
  summary() computes rank-normalized split-R-hat/ESS via posterior when
  available, and degrades to a plain quantile summary with a note when
  not.
- Scalars first: sigma, k, tau (grouped), leaf-count/tree-depth
  summaries if cheaply derivable from stored trees; f(x) draws are
  large - provide the accessor, not automatic diagnostics over n.
- Out of scope: within-fit warnings ("your chains have not converged"
  hard errors); plotting.

## Steps

1. as_draws_array/as_draws_df methods (or a draws() extractor) for
   bart/bart2/rbart fits mapping the native (par x sample x chain)
   arrays to posterior's (iteration, chain, variable) convention -
   note this is the same transpose yhat-shape already pays; do it
   lazily.
2. summary() method gains an R-hat/ESS block for the scalar
   parameters; document the thresholds it flags (rhat > 1.01 noted,
   not enforced).
3. Rd + a vignette-refresh cross-reference (one paragraph there).
4. tinytest: shape round-trips; diagnostics agree with posterior's
   reference values on a fixed fit; graceful degrade without posterior
   (simulate via with_mock or an environment check).

## Verification

- Full tinytest with and without posterior installed.
- R CMD check clean with posterior in Suggests.

## Status (2026-07-07)

Landed (b17a78f). summary.bart/.rbart report rank-normalized
split-R-hat and bulk ESS for the scalar parameters via posterior
when available (verified equal to posterior's own values on the
extracted draws) and degrade to a plain quantile table with an
install note otherwise (degrade path tested via a namespace stub);
as_draws_array/as_draws_df map any chain-dimensioned field to
posterior's conventions through a `vars` argument (scalars by name,
array fields as "field[column]"), reconstructing chains from
combineChains-flattened fits. posterior is Suggests-only: the
methods register in .onLoad through setHook(packageEvent(...)) plus
an isNamespaceLoaded check - a review fix replacing the implementer's
requireNamespace form, which force-loaded posterior's dependency
chain at dbarts load and missed mid-session installs; the hook form
was verified end to end (dbarts load leaves posterior unloaded; a
posterior loaded afterward still dispatches). The hook lives in
R/hooks.R beside .onUnload (zzz.R folded in on review). Deviation:
leaf-count/tree-depth summaries skipped - no cheap tree-walk helper
exists and the plan made them conditional. Gates: tinytest 2458/0
(33 new), both diagnostics paths exercised, air/lintr clean,
checkRd clean, plain R CMD check at the pre-existing compiled-code
NOTE only.
