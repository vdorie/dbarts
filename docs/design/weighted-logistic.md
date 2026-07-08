# Weighted logistic responses

Status: LANDED 2026-07-05.

## Motivation

dbarts accepts observation `weights` for gaussian responses (the residual
precision scales as `sigma^2 / w_i`). For binary responses the engine has
always refused weights: the classic engine fit a weighted probit
incorrectly, and the bartcore bridge kept a blanket refusal over both
binary families. That over-refuses: logistic weights are tractable.

## The two interpretations

Weights admit two readings that coincide for gaussian but diverge for
logistic:

- Frequency weights ("`w_i` copies of observation `i`"). Gaussian: `w`
  copies of `N(f, sigma^2)` give joint likelihood proportional to
  `exp(-w (y - f)^2 / 2 sigma^2)`, identical to scaling the precision by
  `w`. So for gaussian, copies and precision weights are the same thing,
  and that is what dbarts already does. Logistic: `w` copies of a
  Bernoulli give `likelihood^w`, whose exact Polya-Gamma augmentation is
  `omega_i ~ PG(w_i, psi_i)` with working precision `omega_i` and working
  response `w_i (y_i - 1/2) / omega_i`.

- Analysis / precision weights (`w_i * omega_i`, `omega ~ PG(1, psi)`).
  Same posterior mean (the node numerator is `w_i (y_i - 1/2)` either
  way), but the latent precision has the wrong variance: `PG(w, psi)` is
  the sum of `w` independent `PG(1, psi)` draws (variance proportional to
  `w`), while `w * omega` scales a single draw (variance proportional to
  `w^2`). It is therefore not the posterior of `w` copies; only for
  gaussian, where the precision is deterministic, do the two coincide.

Decision (Vincent, 2026-07-05): support the common "multiple copies"
interpretation, so weights mean the same thing (copies) in both families.
That is `omega_i ~ PG(w_i, psi_i)`. For now this is restricted to positive
integer counts (see "Why integer only, for now" below); arbitrary real
weights are deferred.

## Sampling PG(w, psi)

The engine's sampler `ext_rng_simulatePolyaGamma(rng, psi)` draws
`PG(1, psi)` (Devroye's alternating-series method). For a positive integer
`w`, `PG(w, psi)` is exactly the sum of `w` independent `PG(1, psi)`
draws, so weighted logistic reuses the existing sampler with no new
numerical code. This is also the exact path taken by the established
integer-shape implementations (BayesLogit `rpg.devroye`, Makalic and
Schmidt's `pgdraw`).

The cost is O(w) draws per observation per sweep. For the counts use case
(small `w`) this is negligible; very large counts would motivate a normal
approximation to `PG(w, psi)` by the CLT, deferred until a consumer needs
it.

## v1 boundaries

- Weights for a logistic model must be positive integers (observation
  counts). Non-integer weights are not copies and are refused; a zero
  count is a dropped row, also refused. Continuous weighting of a binary
  outcome (e.g. inverse-probability weights) is a different, analysis-
  weight semantics and is out of scope.
- Probit weights stay refused: a weighted probit has no tractable latent-
  variable form.
- Logistic weights are fixed at sampler creation. The between-samples
  weight-mutation surface (`setWeights`, `setData` with new weights) stays
  refused for both binary families in v1, like other creation-fixed
  structure (sparse, grouped, linear designation).

## Why integer only, for now

Non-integer `w` invokes the power likelihood `likelihood^w`, whose exact
augmentation is `PG(w, psi)` for real shape. That is legitimate and is what
gaussian weights already are, but it costs a general real-shape PG sampler
(the engine has only the exact `PG(1, psi)` Devroye sampler), which is a
real numerical component with an approximation-or-large-port tradeoff and a
weaker validation story (you cannot replicate a row 2.5 times, so the
end-to-end "weighted == replicated" check does not exist for real `w`).

The gaussian consistency argument is also weaker than it looks: for
gaussian, "duplicate the observation", "multiply the log-likelihood
contribution", and "scale the precision" all coincide, so gaussian's real
weights are just the copies interpretation extended continuously, not a
separate case-weight design choice.

The use case that genuinely wants real `w` is survey / inverse-probability
weighting (weights `1/pi_i`, non-integer), including weighted-outcome-
regression doubly robust estimators in causal inference (BART's home turf).
But treating survey weights as likelihood weights is a trap: it gives a
consistent point estimate with badly overconfident uncertainty (the model
behaves as if it saw `sum w_i` units, not `n`), because correct survey /
DR inference draws its variance from the design or the efficient influence
function, not the outcome likelihood -- and IPW turns the likelihood into a
data-dependent pseudo-likelihood, so the Bayesian posterior is no longer
coherent or calibrated, which is the property BART causal inference exists
to provide. The mainstream DR recipes (AIPW, TMLE) keep the weights out of
the BART fit for exactly this reason. So refusing non-integer weights loses
little of real value and doubles as a "these are counts, not design
weights" signal; supporting real weights is left as a deferred option to be
paired with loud documentation and weight stabilization if a consumer needs
in-model IPW.

## Refusal matrix

| family   | weights at creation           | weight mutation |
| -------- | ----------------------------- | --------------- |
| gaussian | any positive real (unchanged) | allowed         |
| logistic | positive integers only        | refused (v1)    |
| probit   | refused (intractable)         | refused         |

## Landing notes

- Engine (model.hpp): `LogisticResponse` gained a borrowed `weights_`
  pointer. `refreshLatents` draws `omega_i` as the sum of `lround(w_i)`
  `PG(1, psi_i)` variates and sets the working response to
  `w_i (y_i - 1/2) / omega_i - offset`; `coldStart` seeds `omega` at
  `w_i / 4` (working response `4(y - 1/2) - offset`, weight-independent);
  `restoreLatents` and `setData` carry the weight. Null weights take the
  single-draw path, so unweighted logistic is bit-for-bit unchanged
  (equivalence "logistic" scenario stays identical). chain.hpp passes the
  case weights into the constructor, as gaussian already did.
- Bridge (R_interface_bartcore.cpp): `enforceBinaryWeightPolicy` at the two
  creation sites refuses probit weights and validates logistic weights as
  positive integers; the two mutation sites (setData, setWeights) keep
  refusing binary weight changes with a generalized message (logistic
  weights are creation-fixed in v1).
- R (dbarts.R): the same policy up front -- probit refused, logistic
  weights required positive integers, gaussian unrestricted. data.R
  coerces weights to double so integer count vectors are accepted.
- Tests: component test (deterministic under a fixed seed, recovers the
  monotone signal, positive omega) and test-weighted-logistic.R (fits,
  refusals, and a weight-w fit matching the physically replicated-rows fit
  up to Monte Carlo error); new equivalence scenario "wtlogistic". The
  unweighted path (null weights) is bit-identical to the pre-change code
  (1.0 * x == x, coldStart omega = 0.25 * 1), so the existing "logistic"
  equivalence scenario stays identical; all-ones weights are likewise
  bit-identical to the unweighted fit through the engine (verified via the
  R API across configs and seeds). The component test pins determinism
  (same weights + seed) and signal recovery rather than a weight-1 ==
  unweighted equality: an early two-sampler C++ harness comparison diverged,
  but that traced to the codebase's known heap-layout / SIMD-reduction-split
  sensitivity between two separately-constructed samplers, not real engine
  behavior.
