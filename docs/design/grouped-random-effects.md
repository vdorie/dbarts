# In-core grouped random effects

Wave-2 item, the last one (core-generalization.md section "Wave 2
models"): move rbart_vi's random-intercept Gibbs loop from R into the
engine, composing with any response family. Proposal only; nothing is
implemented.

## Model

What rbart_vi fits today, unchanged:

    y_i = f(x_i) + b_{g(i)} + offset_i + eps_i
    b_j ~ N(0, tau^2),  j = 1..J          (one intercept per group)
    tau ~ built-in prior (see below)

with f the BART forest and eps the response family's error (gaussian
sigma, probit/logistic unit latent). The Gibbs blocks per iteration, in
rbart_vi's order:

1. b | rest: conjugate normal per group. With working response z,
   working weights w, and current total fits F,
       prec_j = (sum_{i in j} w_i) / sigma^2 + 1 / tau^2
       mean_j = (sum_{i in j} w_i (z_i - F_i)) / (sigma^2 prec_j)
   The R loop uses the UNWEIGHTED group mean; in-core uses the working
   weights, which both fixes the (mild) inconsistency for weighted
   gaussian fits and makes logistic compose for free (Polya-Gamma omegas
   are the weights). For unweighted data the two coincide.
2. f | rest: one tree sweep with b_{g(i)} entering as part of the offset
   the trees condition on.
3. sigma | rest (gaussian only), conditioning on F + b.
4. tau | b: slice sampling (stepping-out + shrinkage, the R loop's
   sliceSample) on (0, Inf) of
       -J log(tau) - 0.5 sum(b^2) / tau^2 + log p(tau)
   The R loop takes n.thin slice steps and keeps the last; in-core
   reproduces that count.

Built-in tau priors, on the original response scale with
rel.scale = sd(y) (continuous) or 0.5 (binary):

    cauchy: dcauchy(tau; 0, 2.5 rel.scale), log scale
    gamma:  dgamma(tau; shape 2.5, scale 2.5 rel.scale), log scale

Custom R-function priors cannot cross into C; they keep the R loop (see
Surface).

## Engine

A ResponseModel decorator, per the phase-2 decomposition: grouped
random effects wrap any concrete response model and the backfitting
engine stays written against (z, w). No new template instantiation -
ResponseModel is runtime-virtual, so the decorator costs one virtual
hop on calls that are O(n) inside.

    class GroupedResponse final : public ResponseModel {
      std::unique_ptr<ResponseModel> base_;
      std::vector<std::uint32_t> groupIndex_;   // per observation, 0..J-1
      std::vector<double> groupEffects_;        // b, internal scale
      double tau_;                              // internal scale
      TauPrior prior_;                          // kind + scale, converted once
      std::vector<double> workingResponse_;     // base z minus b_{g(i)}
      ...
    };

- The decorator presents workingResponse() = base z - b_{g(i)} and
  passes workingWeights() through, so trees see the ranef exactly as
  rbart_vi's setOffset made them see it, without touching the user
  offset or the gaussian scale anchoring. offset() stays the base
  model's (user) offset: recorded fits keep their meaning, f(x) only,
  matching rbart_vi's yhat.train = train - ranef convention (the R loop
  subtracts because its ranef rides the engine offset; in-core it never
  enters the recorded fits in the first place).
- refreshLatents(rng, totalFits): shift fits by +b (scratch), let the
  base refresh its latents/PG weights against f + b, then draw b given
  the fresh (z, w) and fits, then tau by slice, then rebuild the exposed
  working response. This runs the blocks in rbart_vi's order relative
  to the tree sweep (b conditions on the previous sweep's fits).
- drawSigma conditions on totalFits + b via the same shifted scratch.
- Scaling: b and tau live on the gaussian internal scale; the prior's
  rel.scale converts once at construction (divide by sigmaScale()) and
  reported draws de-scale like sigma does. Probit/logistic have unit
  scale, so everything is literal.
- The slice sampler is a small free function (stepping-out with width
  ~1 prior sd, then shrinkage), component-tested against R's
  sliceSample on fixed inputs.
- Group structure is fixed at creation (indices 0..J-1 per training
  observation). setData/setResponse on a grouped sampler are refused at
  the bridge (rbart_vi never mutates mid-run; lift later if an embedded
  consumer wants it). The predictor-mutation surface is independent of
  the response model and stays live.

## Surface

rbart_vi keeps its signature and output shape; it gains the in-core
fast path and keeps the R loop as the general fallback:

- In-core when prior is one of the built-in symbols/names (cauchy,
  gamma) AND no callback is supplied. The whole chain then runs inside
  Chain::run - warmup and sampling in one call, chains on worker
  threads instead of PSOCK/FORK clusters, no per-iteration R.
- R loop when prior is a user function or callback is given: identical
  behavior to today, driven per iteration through setOffset as now.
- Results shift under the in-core path (different rng stream), like
  every engine migration this cycle; rbart's hardcoded regression tests
  regenerate. Statistical equivalence gates the switch.

Engine plumbing: SamplerOptions gains groupIndices/numGroups/
tauPriorKind/tauPriorScale/tauSliceSteps (borrowed pointer consumed at
construction like leafCovariateColumns); the factory wraps the family's
response model in GroupedResponse when groups are present. The C entry
points and dbarts()/bart2 surfaces do NOT expose grouping in this pass
- rbart_vi stays the only consumer, constructing through the internal
bridge. Public exposure (a group.by on dbartsData) waits for demand.

## State and reporting

- Results gain per-iteration channels for tau and b (J doubles): the
  bridge preallocates like sigma/varcount, rbart_vi collects
  samples$tau/samples$ranef from them unchanged.
- ChainStateData += ranef (J doubles) and tau; R state slots
  "ranef"/"tau". Restores are exact (values pass through untransformed),
  keeping the bitwise-continuation guarantee.
- predict.rbart is untouched: it already composes f-only predictions
  with per-group ranef draws R-side, drawing unseen groups from
  N(0, tau).

## Costs and risks

- Per iteration the decorator adds O(n) (group scatter/gather, working
  response rebuild) plus O(J) draws plus the slice steps - the same
  arithmetic the R loop does, minus the R overhead and the full
  setOffset scale-refresh churn. Expect the in-core path to dominate
  the R loop wall-clock, especially multi-chain.
- The decorator must not perturb ungrouped paths: it only exists when
  groups are present, and ungrouped construction compiles the same
  code as today (gates: equivalence identical draws, speed baseline).
- Warmup scale evolution: rbart_vi re-anchors the gaussian scale during
  warmup (setOffset updateScale = TRUE) because ranef rides the offset;
  in-core b never touches the response transform, so the anchoring is
  fixed at creation from y alone. This is a (deliberate) difference -
  the transform no longer depends on the warmup trajectory - and is
  part of why draws shift. Flagged for the statistical gate.

## Gates

- Ungrouped paths bitwise unchanged: full suite, equivalence compare
  (identical draws), speed baselines.
- Component tests: conjugate b update against hand-coded R reference
  (weighted and unweighted, probit latents); slice sampler vs
  R sliceSample on fixed rng streams; end-to-end recovery of b and tau
  on simulated grouped data (gaussian + probit + logistic); state
  round trip with ranef/tau; refusals (setData/setResponse).
- rbart_vi in-core vs R loop: statistical equivalence on posterior
  means/intervals of tau, ranef, and fits (z-tests as the equivalence
  harness does across engines).
- R CMD check --as-cran Status OK, zero NOTEs.

## Open decisions (Vincent)

1. Default engine for built-in priors: this proposal flips rbart_vi to
   in-core whenever the prior is built-in and no callback is given
   (results shift, snapshots regenerate). Alternative: opt-in argument
   first.
2. The weighted group update (in-core) vs the R loop's unweighted mean:
   proposal treats weighted as the correct behavior and documents the
   change. Keep the R loop's unweighted update untouched either way.
3. tau slice steps: reuse control@n.thin as the step count (the R
   loop's coupling) or a dedicated option with default equal to it.
4. Exposure beyond rbart_vi (dbartsData group.by, dbarts.h entry
   points): deferred here; any preference on naming reserves it.
