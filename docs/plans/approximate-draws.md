# approximate-draws: Polya-Gamma latents whose exact cost grows with their shape

Status: DESIGN, re-planned 2026-11-03 under dec-B583, dec-B584 and dec-B585, after a blind critique and an
independent check of the exactness derivation. No sizing run has started. What is built waits on the sizing
below; the defaults wait on a later study.

agent: opus for the scratch engine build, the sizing's statistics and each engine slice; one opus review that
runs mutants per engine slice; sonnet for the R-surface slice and its review.
rng: NEUTRAL for the sizing (a scratch build that never lands). Slices: the density, neutral; the exact
accept-or-keep route, shifting where it is selected; real-shape weights, a new capability; approximate
settings, posterior-changing where selected. Nothing moves by default until the maintainer rules the default
(Open call 1), and any default keeps every recorded draw at shapes of 64 or less bitwise.
window: the exact pieces could make 1.0-0 (Open call 2); the approximate settings are after it unless ruled
otherwise.
budget: sizing about 900 lines in a scratch build and about 5 core-hours; the build, if the sizing keeps every
survivor, about 1,900 lines over five slices and about 35 agent-hours ([What is built](#what-is-built)).

## Summary for the maintainer

The problem. A row's hidden Polya-Gamma value is drawn exactly today as one small draw per trial, about 160
nanoseconds each. A logistic fit of 200 rows spends about 2 microseconds a row per sweep on everything else, so
at 100 trials a row the sweep is 9 times slower and at 10,000 trials 800 times. Binomial responses, frequency
weights and negative binomial counts all carry such shapes, and real-valued weights and shapes have no exact
draw at all today.

What you ruled. An accept-or-keep step whose long-run distribution is correct is exact. A sampler beaten on both
accuracy and cost by another is dropped; the rest are sized by effective draws per second of whole fits, and
only those that change what a user gets per second are built, each chosen by a setting, until a study sets the
defaults.

The screen. Of the samplers in the literature and the ones this design constructed, four survive: today's exact sum;
the exact accept-or-keep route, whose cost does not grow with the trials; an approximation, the first few terms
of the hidden value's exact series plus one gamma for the rest; and a plain normal. Dropped: the saddle-point
sampler (the approximation with 10 terms is more accurate at every shape measured and about a third of its
cost), BayesLogit's sum of 1,000 gammas (1.2 to 7,600 times less accurate than the 1-term approximation at 500
times its cost), and the series without its correcting gamma.

The sizing. In a scratch copy of the engine that never ships, each survivor is timed on fits shaped like users'
(rows, trees, trial counts from 16 to 100,000, negative binomial counts with fixed and estimated shape, weights
below and above 1) and its effective draws per second measured on the fitted values, the probabilities and a
tail quantile. A survivor is built if, somewhere users fit, it gives at least 1.25 times the effective draws per
second of every survivor at least as accurate. About 5 core-hours. The expected result, from per-draw costs: the
exact route beats the sum from about 25 trials up, by about 250 times per sweep at 10,000; the approximation beats the exact
route by 2.2 to 2.9 times at 75 trees, less with more trees; the normal beats nothing.

What it may build. The exact route with its exact density (logistic, binomial, logistic hazard and negative
binomial), real weights wherever the exact route reaches them (given when the fit is created; changing a
fractional weight on a running sampler is refused until an exact fresh draw at such shapes exists, since none is
published below 1), and, if the sizing keeps them, the approximate settings. Multinomial fits keep the exact sum under every exact setting: their hidden values are temporary and
the accept-or-keep step does not hold there.

Your calls: the default until the study (recommended: the exact route above 64 trials, the sum below); what
ships in 1.0-0 (recommended: the exact pieces if their gates pass in time, the approximate settings after); and
the settings' name and values.

## Goal

Every Polya-Gamma family draws exactly at a cost that does not grow with the shape wherever the exact route
holds; real weights reach those families where the draw is exact (dec-B568); and the approximate samplers that
change what a user gets per second exist as settings, for a later study to choose defaults from.

## Context

Rulings, quoted, not reopened:

- dec-B583: an accept-or-keep latent is exact; the maintainer, 2026-11-03: "As far as a ruling goes, exact is
  exact, right? It either has the correct distribution or it doesn't." Where the step keeps the distribution
  "is the design's to show".
- dec-B584, revised by dec-B585: "Oh, I don't know if we need to build them all. Maybe we do? Are some strictly
  dominated? Maybe size that to see what we actually need for ESS / sec". A dominated sampler is dropped, the
  rest sized by effective draws per second of whole fits, only those that change what a user gets per second
  built, each chosen by a setting, until a study sets the defaults. Not ruled: the default until the study,
  whether the approximate settings ship in 1.0-0 against dec-B14, and the settings' shape.
- dec-B14: "Every draw in 1.0 stays exact: an approximate Polya-Gamma draw ... is declined."
- dec-B568: real frequency weights where the draw is exact; the Polya-Gamma families whole numbers "until a
  real-shape draw exists". The maintainer has noted that this work also unlocks some non-integer weights.
- dec-B567 (the binomial response, shape n w per row), dec-B563 (frequency weights as shape times w on the
  Polya-Gamma families), dec-B143, dec-B15 (integer dispersion), dec-B202 (the one-million cap, a consequence to
  revisit, not in scope).

Engine: the seam [`simulatePolyaGammaShape`](../../src/bartcore/model.hpp) sums round(b) draws of
[`ext_rng_simulatePolyaGamma`](../../src/external/random.c). [`LogisticResponse::refreshLatents`](../../src/bartcore/model.hpp)
(logistic, binomial, hazard.logistic) and [`MultinomialForestCombiner::drawForestGlue`](../../src/bartcore/combiner.hpp)
sum in their own loops. [`NBResponse::refreshLatents`](../../src/bartcore/model.hpp) draws r collapsed over omega
([`NBShapePrior`](../../src/bartcore/model.hpp)), then [`NBResponse::drawOmega`](../../src/bartcore/model.hpp) at the
new r. The same latent paths serve the setters (setWeights, setActiveRows, setResponse, reapplying weights) and
the bridge's [`drawAugmentation`](../../src/R_interface_bartcore.cpp).

Background: [nbinom-large-shape-pg.md](nbinom-large-shape-pg.md), the exact construction designed in September
(saddle-point envelope, certified density, WPS Conjecture 15), not built. Probes: scratch/apx-design/probes
(this design), scratch/apx-crit (the critique), scratch/apx-check (the independent check).

### Where large or real shape arises, and today's cost

Measured 2026-11-03 on the laptop (arm64, load 5 to 9, per-round minima), shipped build of main.

| source | shape per row | example users meet | cost today |
|---|---|---|---|
| binomial trials (dec-B567) | n w | NC SIDS: births per county 248 to 21,600; schools: 101 to 4,120 pupils | 160 ns per trial |
| logistic, hazard.logistic frequency weights | w | collapsed identical rows: Titanic's cells to 670 | the same |
| nbinom | (y + r) w | Migration (carData): counts to 99,400 | the same |
| multinomial trials | n_i per category | votes over parties | 480 ns per trial over 3 categories |
| real shape | real w: power priors, fractional counts; real r | five of the nbinom study's twelve datasets have a shape below 1 | refused |

Logistic fit, 200 rows, 75 trees, one chain: 2 us a row per sweep at w = 1, 18 at 100, 166 at 1,000, 1,630 at
10,000; nbinom fixed r 160 at 1,000, 1,590 at 10,000; multinomial, 3 categories, 480 at 1,000 trials.

## The dominance screen

TV is total variation distance to PG(b, psi), computed against the exact density by Fourier inversion
(probes laterr.R and screen.R, accurate to about 1e-13), at psi 0 and 8 (2 where shown). Per-draw costs on the
laptop: a Devroye PG(1) 146 to 170 ns, a gamma about 41 ns, a normal with its moments 35 ns. The series is
PG(b, psi) = sum_k g_k / d_k, g_k ~ Gamma(b, 1), d_k = 2 pi^2 (k - 1/2)^2 + psi^2 / 2 (Polson, Scott and Windle
2013); "series K" draws the first K terms exactly and the rest as one gamma with the rest's exact mean and
variance (the PG moments less the first K terms'). Series 0 is the moment-matched gamma.

| candidate | TV, b = 16 | 64 | 1e3 | 1e4 | cost per draw | verdict |
|---|---|---|---|---|---|---|
| exact sum | 0 | 0 | 0 | 0 | 160 ns times b | survives (exact; cheapest below about 25) |
| exact accept-or-keep route (below) | 0 | 0 | 0 | 0 | estimated 3 to 4.5 us, flat in b | survives (exact) |
| series 10 | 2.5e-7 | 1.2e-7 | 2.9e-8 | 9.1e-9 | about 450 ns | survives as the approximation's high setting |
| series 5 | 7.1e-6 (psi 8, b = 16) | 3.2e-6 | 8.0e-7 | 2.5e-7 | about 250 ns | survives |
| series 1 | 2.2e-3 | 1.0e-3 | 2.6e-4 | 8.2e-5 | about 80 ns | survives |
| series 0, moment-matched gamma | 1.5e-2 | 7.7e-3 | 1.9e-3 | 6.1e-4 | about 41 ns | survives with the normal |
| moment-matched normal (WPS 2014 sec 6; BayesLogit above 170; Python polyagamma above 50) | 6.3e-2 | 3.1e-2 | 7.8e-3 | 2.5e-3 | 35 ns | survives on cost by 6 ns; expected to lose on effective draws per second |
| saddle-point sampler (WPS 2014 sec 5) | - | 2.2e-5 | 3.5e-7 | 1.1e-8 | about 1.5 us (WPS Table 1: about 9 Devroye draws at n = 100; not run here) | dropped: series 10 is more accurate at every b measured, at under a third of the cost |
| truncated sum, 1,000 gammas, no remainder (BayesLogit's default from 3 to 13, per the critique's reading of its source) | 2.6e-3 | 5.2e-3 | 2.1e-2 | 6.5e-2 | about 41 us | dropped: series 1 is more accurate everywhere measured at 1/500 of the cost |
| truncated sum, 200 gammas, no remainder | - | 2.6e-2 | 0.10 | 0.32 | 8 to 13 us | dropped, as above |
| WPS alternate sampler (sec 4) | 0 | 0 | 0 | 0 | linear in b, faster than the sum by 1.2 to 2.7 for n 3 to 16 (WPS Table 1) | dropped for whole b: the exact route is flat in b; exact only if WPS Conjecture 8 holds (verified on h in [1, 4]) |
| the September exact construction | 0 | 0 | 0 | 0 | estimated 2 to 4 us | dropped: rests on WPS Conjecture 15, and gives independent draws only where the accept-or-keep route already serves |

At small real shapes the series needs more terms: at b = 0.5, series 10 is 1e-5 to 4e-5 and series 20 about
1.4e-6. The normal survives only because a normal draw is 6 ns cheaper than a gamma; at 2 us a row of other work
that is 0.3 percent of a sweep, against a TV four times the gamma's.

## The exact accept-or-keep route

Each sweep, after the trees, each row proposes omega' from series 3 re-matched at its current psi (K grows with
|psi|, below), and keeps omega' with probability min(1, f(omega') q(omega) / (f(omega) q(omega'))), otherwise
keeps omega.

- f is PG(b, psi), an exact exponential tilt of PG(b, 0): f's value at the current omega is cached per row as
  log PG(b, 0)(omega) and tilted to the current psi analytically. q is not a tilt: re-matching the remainder at
  psi changes its shape, and matching it once at psi = 0 and tilting is poor (probe tiltq.R: at b = 64,
  psi = 20, acceptance 0.66 tilted against 0.9993 re-matched; at b = 1e3, psi = 20, 0.08 against 0.9998). So q is
  evaluated at both omega' and omega each sweep: three density evaluations a row, f and q at omega' in one
  pass over shared nodes, q at omega in a second.
- Density: the trapezoid rule on the vertical line through the saddle point, 33 nodes, measured at 1.2 to 1.5
  us per evaluation, flat in b from 16 to 1e5 (scratch/apx-crit/pgdens.c, at load 19 to 25). Its accuracy is a
  bias floor, not a formality: an error d in the ratio moves the chain's distribution by about d, and with 33
  nodes d is 1e-11 at b >= 64 but 1e-7 at b = 16. The decision is made certified: the density comes with an a
  priori error bound (analyticity in a strip around the line, bounded truncated tails, as nbinom-large-shape-pg.md
  derives), the step accepts or rejects only when the uniform clears the ratio by more than the bound, and
  otherwise refines the nodes and decides again. The result is exact down to floating-point rounding, the
  standard Devroye's sampler meets. Node counts rise as b falls (the transform decays like exp(-b sqrt(2|t|)));
  the sizing measures them at b in {0.3, 0.7, 1.5, 4, 16}.
- The proposal's weight f/q is bounded for every b > 0, psi and K >= 1 (its first term carries PG's exact
  right-tail rate pi^2/2 + psi^2/2; the left tail of f falls faster than q's), but not uniformly: its right-tail
  limit is exp(b D(psi, K)), D > 0, which at K = 3 is 1.06 at b = 100, psi = 8 and above 290 at b = 1e4, psi = 8
  (scratch/apx-check tail-ratio.R). A chain landing far in that tail would stay a while; the bulk acceptance is
  unaffected (at least 0.9998 at b 16 to 1e3, psi to 20, K = 5, and 0.997 at b = 1; probe tiltq.R). K is therefore chosen per row so
  that b D(psi, K) <= 0.01: K = 10 holds that to about b = 5,000 at psi 8, and D falls roughly like K^-5. This
  keeps each step uniformly ergodic with a constant the design controls; the cost is K + 1 gammas, under 1.3 us
  at K = 30.
- Rejection from the same proposal with bound M would give independent exact draws, which would also serve
  multinomial and the setters below. It needs a certified M. The tail limit is not one: on a grid, sup f/q
  exceeds exp(b D) at b <= 1 (1.10 against 1.0002 at b = 0.3, psi = 8, about one sd below the mean), though not
  at b >= 16 where checked (probe supcheck.R). Recorded as research, not built.

Where the route keeps the posterior (scratch/apx-check report, with exact finite-state probes):

- Within a sweep it holds: rows are independent given the trees, and nothing between the trees and the latent
  step reads omega.
- Estimated-r nbinom holds under this design's rule: on a sweep where r moved, every row draws fresh and
  exactly; where r stayed, the step applies. r's draw never reads omega, so on the event that r stayed the old
  omega is a draw from its conditional at that r. Checked exactly: the rule's error about 1e-16, always applying
  the step 0.22 in TV (pcg-exact.R); the September critique's contrary reading was unconditional. The same
  holds for any future r step that leaves r's conditional invariant without reading omega.
- Multinomial does not hold: [`MultinomialForestCombiner::drawForestGlue`](../../src/bartcore/combiner.hpp) draws
  each category's omega as a temporary latent against the current margins just before that category's forest,
  and there is no persistent joint omega to keep. The route there gives TV 0.023 in the check's exact toy.
  Multinomial keeps the exact sum under every exact setting; an approximate setting, being a fresh draw, can
  serve it.
- Setters: a row whose conditional changes between sweeps is drawn fresh and exactly AT the setter, before the
  trees move (a draw at the next refresh is too late: the trees would consume the stale omega). Today's
  setWeights, setActiveRows, nbinom setResponse and reapplied weights already do; they must stay on an exact path
  separate from the route, since refreshLatents, NBResponse::drawOmega and the seam are shared with them and with
  the augmentation helpers. Cold and warm starts, installShape, setData and state restores are starting points
  and need no fresh draw.
- The per-row density cache is scratch, recomputed from omega; it is invalidated on restoreLatents, nbinom's
  cancel-restore, installShape, restoreShape, setData, any change of b and any cold start.
- Recorded, not planned here: LogisticResponse::setOffset, NBResponse::setOffset and Sampler::setPredictor keep
  omega today, which is exact only if the outer step that moved them read omega; that is a question about the
  outer sampler's contract, the same with exact draws, taken separately.

### Real shapes

The route works at any real b > 0 given the density, which the saddle-line inversion supplies at every x > 0
with more nodes as b falls (the alternating series cancels catastrophically for large x, and no gamma-kernel
series exists at non-integer b). That is the exact draw dec-B568 waits on, within sweeps and at creation (a cold
start). A fresh exact draw at a fractional shape is another matter: none is published below 1, and WPS's
alternate sampler above 1 rests on Conjecture 8. So, under every exact setting:

- Real weights on logistic, binomial and hazard.logistic are taken at creation; a setter that would need a fresh
  draw at a fractional shape (setWeights changing such a row, setActiveRows reactivating one) is refused by name,
  as dec-B568 takes real weights only where the draw is exact. Rejection from the proposal, once its bound is
  certified, would lift that.
- Real negative binomial dispersion (TODO negbin-real-dispersion) needs a fresh draw whenever r moves, so it
  needs either that rejection draw or an r step that reads omega (a Metropolis step on r given omega, whose
  ratio uses the same density); its own plan decides. Not in this build.

## The sizing

On a scratch build carrying every survivor behind a compile-time or environment switch, never landed. Fits, at
most 2 cores, run on an x86-64 Linux host for timing (quiet, alternating arms, per-round minima, load stamped):

| family | rows | trees | shape per row | arms |
|---|---|---|---|---|
| logistic with binomial rows | 200, 2,000 | 75, 200 | 16, 64, 256, 1e3, 1e4, 1e5 (the sum to 1e4 at 2,000 rows), and SIDS-like spread (lognormal, median 2,000) | sum, route, series 1, 5, 10, gamma, normal |
| logistic with real weights | 2,000 | 75 | uniform on (0.2, 0.8) and on (1.2, 3) | route, series 10 and 20, gamma |
| nbinom, fixed r = 1e3 | 200, 2,000 | 75 | counts mean 10 and 1e3 | sum, route, series 5, 10 |
| nbinom, estimated r | 200, 2,000 | 75 | counts mean 20 and 1e3 | sum, route (fresh on r moves), series 5 |
| multinomial, 3 categories | 200 | 75 | 1e3 trials | sum, series 5 |

Measured: seconds per sweep, all cells; effective draws per second (rank-normalized bulk ESS of f(x_i) and of
the probability, tail ESS of a 5 percent quantile; Vehtari et al. 2021) on 8 cells spanning the table, 2 chains
of 2,500 kept sweeps each; route acceptance and the share of sweeps where estimated r moves; the density's node
counts at the real shapes. Effective draws per sweep should not differ between arms beyond Monte Carlo error,
since the latent's distribution does not touch the trees' mixing beyond its own small error; if they do in the
8 cells, every cell is measured for ESS.

Stated before any run: a survivor is built if, in at least one cell, its effective draws per second (the
smallest of the three quantities) are at least 1.25 times those of every other survivor at least as accurate
(accuracy ordered exact, series 10, 5, 1, gamma, normal). 1.25 is about five times the per-sweep timing grain on
the host (about 3 percent single-run, less over 5 rounds) and well above the ESS estimates' error at these
lengths (about 5 percent), so a survivor passing it is faster than its betters beyond doubt; it is this plan's
bar, and the result tables show every ratio so the maintainer can read another. The exact sum is kept
regardless (multinomial, setters, the augmentation door). The route is built if it passes or if it is the only
exact draw for a cell (real weights).

Cost: about 900 lines in the scratch build (density with its certified bound about 300, the route with cache,
K rule and the moved-r rule about 250, the approximations about 100, switches about 50, tests/cpp about 200);
compute about 5 core-hours (timing 2, ESS cells 2, real shapes, nbinom and multinomial 1).

## What is built

Conditional on the sizing; each sampler chosen by one setting (Open call 3). Tiers by
[Process by risk](README.md#process-by-risk).

| slice | content | tier | lines | agent-hours |
|---|---|---|---|---|
| S1 | the certified PG density and its tests/cpp pins (agrees with inversion and the alternating series where that is accurate; integrates to 1; moments) | engine, neutral | 400 | 6 |
| S2 | the route at the seam for logistic, binomial, hazard.logistic, fixed- and estimated-r nbinom; the separate exact path for setters and augmentation; cache and invalidation; the K rule; the setting's engine side | engine, shifting where selected | 750 | 12 |
| S3 | real weights under the route on logistic, binomial, hazard.logistic; the fractional-setter refusals | engine plus R surface | 250 | 6 |
| S4 | the approximate settings kept by the sizing (series at its K, and the normal if kept), all PG families including multinomial | engine, posterior-changing where selected | 200 | 5 |
| S5 | the setting at the doors users reach, help, NEWS | R surface | 300 | 5 |

Gates for S2 and S3 (exactness), each with the exact sum as control, passing first:
- joint-distribution check (Geweke 2004) as [geweke-mc.R](../../benchmarks/R/geweke-mc.R) runs it, on logistic
  with binomial rows, fixed- and estimated-r nbinom and real weights; broken builds that must fail: the step
  applied after r moved, a stale cache after psi changes, q tilted from psi = 0, a setter's redraw routed
  through the step;
- two-sample KS of the route's omega at fixed (b, psi) against the exact sum, b <= 1e3, and moments to the
  third cumulant above;
- the exact single-tree posterior of [logistic-reference.R](../../benchmarks/R/logistic-reference.R) part 1
  with rows of b trials;
- simulation-based calibration through [sbc.R](../../benchmarks/R/sbc.R)'s logistic machinery with binomial
  rows, its exact arm passing first;
- the engine tier's battery: tests/cpp, full tinytest, sanitizers on tests/cpp and the R-loaded path, R CMD
  check, the arm64 snapshots, the equivalence baseline bitwise at shapes of 64 or less;
- mutants: the review runs the mutation battery's anchors plus mutants of the accept comparison, the K rule,
  the cache invalidation, the moved-r branch and the setter path.

S4's gates: its own latent TV pinned in tests/cpp against inversion; the exact gates stay on exact settings.

If the sizing drops the approximations, S4 is not built and S5 offers the exact settings only (about 1,700
lines, 30 agent-hours). If it drops the route at some shapes (it will below about 25), the default rule there
is the sum.

## The defaults study, after the build (outline)

It measures, on fits shaped like users' and on real covariates with real trials, each built setting's
effective draws per second and, for the approximate ones, their effect on what users read: coverage and width of
90 and 95 percent intervals for probabilities, rates and contrasts, RMSE and bias, against exact posteriors
(the one-parameter toy by quadrature, the single-tree enumeration) and semi-synthetic truths, at equal running
time. Choosing among exact settings changes no posterior and is a speed finding. Making an approximate setting a
default is a change users' results carry, and comes as a short research report under the studies standard
(dec-B559): at most 5 pages, readable alone, two referees first.

## Decision

### Open call 1: the default until the defaults study

Background. Once built, a fit's hidden values can be drawn three ways: the exact sum (today, cost growing with
the trials), the exact accept-or-keep route (cost flat, about 3 to 4.5 microseconds a row whatever the trials),
or an approximation. Both exact ways give the same posterior; they differ in speed and in the exact draws a seed
produces. The sum is cheaper below about 25 trials a row, the route above. Logistic weights, binomial responses
and negative binomial fits are not in 0.9-34, so no earlier release's draws are owed.

- (a) The route above 64 trials (or count plus shape) a row, the sum at or below. Every recorded draw stays as
  it is; a fit of 1,000 rows at 10,000 trials goes from about 1.8 core-hours to under a minute. Costs the
  gates of S2 before it can be the default.
- (b) The sum everywhere, the route by setting only. Nothing moves; large-trial fits stay slow unless the user
  finds the setting.
- (c) An approximation above a shape. Faster again by an estimated 2.2 to 2.9 times, but a default that is not exact, before
  the study that is to decide it.

Recommendation: (a). It is exact, faster wherever it applies, and needs nothing from users.

### Open call 2: what ships in 1.0-0

Background. Your ruling for 1.0 reads "Every draw in 1.0 stays exact"; the accept-or-keep route and real
weights under it are exact under your later ruling, so they fit 1.0-0 if built and gated in time. The
approximate settings do not.

- (a) S1 to S3 and S5's exact settings in 1.0-0 if their gates pass in time (about 30 agent-hours); the
  approximate settings after 1.0-0.
- (b) Everything in 1.0-0, the approximate settings as opt-in only, revising the 1.0 clause of the earlier
  ruling (about 35 agent-hours).
- (c) Nothing from this arc in 1.0-0; the cost stays documented as now.

Recommendation: (a). It keeps the 1.0 promise as written and gives binomial and weighted fits their speed and
real weights at release.

### Open call 3: the setting's name and values

Background. Elsewhere: BayesLogit exposes separate functions per method (its default switches between them
silently); the Python polyagamma package takes method = with "devroye", "alternate", "saddle" and "gamma",
None meaning its hybrid; R's glm and the other model fitters users know expose nothing of the kind. In dbarts,
how the sampler draws lives in the control settings, beside the move probabilities. Under the package's rule
for absent values, NULL is the default rule.

- (a) A control setting, latent.sampler = NULL, NULL for the default rule, otherwise "sum", "metropolis" or
  "series" ("normal" too if the sizing keeps it), reachable from bart, dbarts and xbart as the other control
  settings are. The number of series terms is chosen inside from an error target, not exposed.
- (b) Two settings: one choosing exact or approximate, one choosing the method within it. Harder to read; a user
  wanting speed sets two.
- (c) A family argument, such as logistic(latent = ...). Puts a sampler choice on the model, and the families
  without a Polya-Gamma latent would have nowhere to write it.

Recommendation: (a), named for what it chooses, as polyagamma's method does, with values naming the method.

## Constraints

- No package code from the sizing lands; slices follow their tier.
- Any default keeps every recorded draw at shapes of 64 or less bitwise.
- Draws come only from the chain's RNG; the cache is scratch, never state.
- Nothing from BayesLogit's source (GPL 3; dbarts is GPL 2 or later).
- Out of scope: the probit binomial, real negative binomial dispersion, the dec-B202 cap, the setOffset and
  setPredictor contract question.

## Steps

1. Build the scratch sizing build; review its density and route code before any run.
2. Run the sizing; write its tables and the build-or-not verdict per survivor under the stated rule.
3. Open calls to the maintainer with the measured ratios.
4. S1, then S2 and S3, then S4 if kept, then S5, each through its tier.
5. The defaults study, planned once the build lands.

## Verification

- The screen's numbers reproduce from scratch/apx-design/probes (laterr.R, screen.R, tiltq.R, supcheck.R) and
  scratch/apx-check/probes (tail-ratio.R, pcg-exact.R).
- The sizing's verdict follows from its tables by the rule above, mechanically.
- S2's broken builds fail the joint-distribution check; a gate passing one has no power and stops the slice.

## Taken and not taken

- From the critique (September design's revision): the exact-moment remainder and the re-ranking; the
  unbounded weight of the gamma proposal; certified meaning; exact draws at swaps.
- From the independent check: estimated-r nbinom holds under the moved-then-fresh rule (as this design argued,
  against the critique); multinomial excluded; fresh exact draws at the setter on a separate exact path; the
  cache's invalidation list; f and q both at the current psi (checked: tilting q is poor, so q is recomputed);
  the weight bounded but not uniformly, so K grows with b and |psi|; the density's error as a bias floor and
  the certified decision; real shapes (none exact below 1, the alternate's conjecture, the saddle-line
  inversion); the setOffset and setPredictor question recorded, not planned.
- Not taken: nothing in the check is disputed. Added here: rejection from the proposal as an independent exact
  draw, whose tail limit fails as a bound at b <= 1 (so it is research, not built), and the fractional-setter
  refusal that follows from it.

## Not verified

- The route's whole cost: three density evaluations a row is an estimate from one measured under heavy load.
- The density's node counts and error bound at real b below 16.
- How often estimated r moves per sweep, which sets estimated-r nbinom's cost under the route.
- The saddle-point sampler's cost (WPS Table 1) and BayesLogit's switching points (the critique's reading).
- That effective draws per sweep match across arms, which the sizing checks.
- The far-tail behavior of f/q at b <= 1 on finer grids: supcheck.R's largest ratios at b = 0.3, psi = 0 sit
  where the inversion loses accuracy and are not read.
