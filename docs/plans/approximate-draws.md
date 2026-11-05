# approximate-draws: Polya-Gamma latents whose exact cost grows with their shape

Status: DESIGN, re-planned 2026-11-03 under dec-B583, dec-B584 and dec-B585; amended 2026-11-05 from the sizing
study and its independent review, and for dec-B586 (augmentation variables drawn at the start of each sweep).
The sizing has run; what is built and when wait on the three open calls; the defaults wait on a later study.

agent: opus for each engine slice; one opus review that runs mutants per engine slice; sonnet for the R-surface
slice and its review.
rng: the density, neutral; the exact accept-or-keep route, shifting where it is selected; real-shape weights, a
new capability; approximate settings, posterior-changing where selected. Nothing moves by default until the
maintainer rules the default (Open call 1), and any default keeps every draw of a fit whose shapes are all whole
numbers of 64 or less as the build before it draws them. S2 lands after the latents-first plan's engine slices.
window: the exact pieces could make 1.0-0; whether the approximate settings do is Open call 2.
budget: the sizing ran at 970 lines in a scratch build and 1.6 core-hours. The build, with every survivor and the
recommended setting, about 2,350 lines over five slices and about 40 agent-hours; the exact pieces alone about
2,000 lines and 33 agent-hours ([What is built](#what-is-built)).

## Summary for the maintainer

The problem. Each row of a logistic, binomial, logistic-hazard or negative binomial fit carries a hidden value
that is redrawn every sweep. Today it is drawn exactly as one small draw per trial, 130 to 160 nanoseconds each.
Everything else a sweep does costs 0.13 to 1 microsecond a row at 75 to 200 trees, so from a few dozen trials a
row the hidden values are most of the sweep, and at 10,000 trials a fit is hundreds of times slower than it need
be. Real-valued weights and shapes have no exact draw today and are refused.

What you ruled. An accept-or-keep step whose long-run distribution is correct is exact. A sampler beaten on both
accuracy and cost is dropped; the rest are sized by effective draws per second of whole fits, and only those that
change what a user gets per second are built, each chosen by a setting, until a study sets the defaults. Since
then: the hidden values are drawn at the start of each sweep, before the trees, so a row whose distribution a
setter changed gets its exact draw there.

The candidates. Today's exact sum; the exact accept-or-keep route, which proposes from a close approximation and
keeps or refuses against the exact density, at a cost that does not grow with the trials; an approximation that
draws the first few terms of the hidden value's exact series and one gamma for the rest; and a plain normal. The
saddle-point sampler and BayesLogit's truncated sums were dropped as beaten on both counts.

What the sizing measured. On a scratch copy of the engine, each candidate was timed on 39 fits shaped like users'
(200 and 2,000 rows, 75 and 200 trees, 16 to 100,000 trials, negative binomial counts with fixed and estimated
shape, fractional weights, a 3-category multinomial) on one core of an x86 Linux machine, and its effective draws
per sweep measured on 8 of them. Effective draws per sweep did not differ between candidates beyond noise, so
speed per sweep decides. Found:

- The route costs 3 to 9 microseconds a row (more where an estimated dispersion moves often). The sum is 3 to
  4 times faster at 16 trials, the route twice as fast at 64, 31 to 38 times at 1,000 and about 300 times at
  10,000. Between 16 and 64 was not timed.
- The approximation with 10 exact terms is 3 to 7 times faster than the route at 64 trials and up, 58 times at
  estimated-shape negative binomial with 200 rows, and 173 times faster than the sum in the multinomial fit,
  which has no exact route.
- Built: the sum, the route, and the approximation (at 20, 10, 5 and 1 terms; the 5- and 1-term versions only
  just, and only at 75 trees). Dropped: the moment-matched gamma alone and the normal, which never gained a
  quarter over a more accurate candidate. Those verdicts hold for 75 trees and up; at fewer trees, which the
  grid did not measure, arithmetic suggests the gamma and the normal could pass.

What the review found. An independent review recomputed the verdict and confirmed it, and found the scratch
route wrong in a band of shapes no timed fit reached: shapes from about 2 to 12 with log-odds beyond about 80 in
size, where its density for the proposal was off while claiming to be accurate. The density the package ships
must therefore carry a true error bound for both densities at every shape and log-odds a fit can reach, and the
step must never decide on a value it could not compute; the gates now test that band. The approximation is also
far less accurate at small shapes and large log-odds than at the shapes first screened: with 10 terms its
distance from the exact distribution ranges over more than 10,000 times across rows.

What it builds. The certified density; the exact route for logistic, binomial, logistic hazard and negative
binomial; fractional weights under the route; the approximation as a setting, for every Polya-Gamma family
including multinomial; and the setting at the doors users reach. Under an exact setting, a change made between
sweeps that would leave a fractional-weight row needing a fresh exact draw (a new offset, say) is refused, since
none exists below a shape of 1.

Your calls: the default until the study (recommended: the route for rows above 64 trials or with fractional
shapes, the sum for the rest); what ships in 1.0-0 (recommended: everything, the approximation opt-in only); and
the setting's shape (recommended: one setting naming the method, or giving the approximation's accuracy as a
number, so that its meaning is the same in every row).

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
- dec-B586: augmentation variables, the Polya-Gamma precisions among them, are drawn at the start of each sweep,
  after the per-sweep hook and before the trees; it notes that an accept-or-keep latent "still needs one exact
  draw at the next sweep for each row a setter moved". The maintainer, 2026-11-03: "Ok, do it and make sure it is
  well documented". Its plan (latents-first, on its branch) states the constraint this plan meets
  ([Setters and marked rows](#setters-and-marked-rows)).
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
the bridge's [`drawAugmentation`](../../src/R_interface_bartcore.cpp); after latents-first the setters draw nothing.

Background: [nbinom-large-shape-pg.md](nbinom-large-shape-pg.md), the exact construction designed in September
(saddle-point envelope, certified density, WPS Conjecture 15), not built. Probes: scratch/apx-design/probes
(this design), scratch/apx-crit (the critique), scratch/apx-check (the independent check), scratch/apx-size (the
sizing: report, results, probes; scratch build on branch apx-size, never lands), scratch/apx-rev/probes (the
sizing's review), scratch/apx-amend/probes (this amendment).

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
10,000; nbinom fixed r 160 at 1,000, 1,590 at 10,000; multinomial, 3 categories, 480 at 1,000 trials. On the
x86 box under the scratch build the sum costs about 0.13 us per trial and everything but the latent 0.13 to
1.0 us a row (below).

## The dominance screen

TV is total variation distance to PG(b, psi), computed against the exact density by Fourier inversion
(probes laterr.R and screen.R, accurate to about 1e-13), the larger of psi 0 and 8 (2 where shown). psi is the
row's linear predictor (the log-odds for logistic). Per-draw costs on the laptop: a Devroye PG(1) 146 to 170 ns,
a gamma about 41 ns, a normal with its moments 35 ns. The series is PG(b, psi) = sum_k g_k / d_k,
g_k ~ Gamma(b, 1), d_k = 2 pi^2 (k - 1/2)^2 + psi^2 / 2 (Polson, Scott and Windle 2013); "series K" draws the
first K terms exactly and the rest as one gamma with the rest's exact mean and variance (the PG moments less the
first K terms'). Series 0 is the moment-matched gamma.

| candidate | TV, b = 16 | 64 | 1e3 | 1e4 | cost per draw | verdict |
|---|---|---|---|---|---|---|
| exact sum | 0 | 0 | 0 | 0 | 160 ns times b | survives (exact) |
| exact accept-or-keep route (below) | 0 | 0 | 0 | 0 | measured 3 to 9 us a row (the sizing) | survives (exact) |
| series 10 | 2.5e-7 | 1.2e-7 | 2.9e-8 | 9.1e-9 | about 450 ns | survives |
| series 5 | 7.1e-6 (psi 8, b = 16) | 3.2e-6 | 8.0e-7 | 2.5e-7 | about 250 ns | survives |
| series 1 | 2.2e-3 | 1.0e-3 | 2.6e-4 | 8.2e-5 | about 80 ns | survives |
| series 0, moment-matched gamma | 1.5e-2 | 7.7e-3 | 1.9e-3 | 6.1e-4 | about 41 ns | survives with the normal |
| moment-matched normal (WPS 2014 sec 6; BayesLogit above 170; Python polyagamma above 50) | 6.3e-2 | 3.1e-2 | 7.8e-3 | 2.5e-3 | 35 ns | survives on cost by 6 ns |
| saddle-point sampler (WPS 2014 sec 5) | - | 2.2e-5 | 3.5e-7 | 1.1e-8 | about 1.5 us (WPS Table 1: about 9 Devroye draws at n = 100; not run here) | dropped: series 10 is more accurate at every b measured, at under a third of the cost |
| truncated sum, 1,000 gammas, no remainder (BayesLogit's default from 3 to 13, per the critique's reading of its source) | 2.6e-3 | 5.2e-3 | 2.1e-2 | 6.5e-2 | about 41 us | dropped: series 1 is more accurate everywhere measured at 1/500 of the cost |
| truncated sum, 200 gammas, no remainder | - | 2.6e-2 | 0.10 | 0.32 | 8 to 13 us | dropped, as above |
| WPS alternate sampler (sec 4) | 0 | 0 | 0 | 0 | linear in b, faster than the sum by 1.2 to 2.7 for n 3 to 16 (WPS Table 1) | dropped for whole b: the exact route is flat in b; exact only if WPS Conjecture 8 holds (verified on h in [1, 4]) |
| the September exact construction | 0 | 0 | 0 | 0 | estimated 2 to 4 us | dropped: rests on WPS Conjecture 15, and gives independent draws only where the accept-or-keep route already serves |

The series' accuracy falls at small shapes and large |psi|, by far more than the columns above show
(scratch/apx-rev/probes tv-small.R, re-run 2026-11-05; the b = 16, psi 8 row reproduces the table above):

| b | psi | series 1 | series 5 | series 10 | series 20 |
|---|---|---|---|---|---|
| 0.3 | 0 | 5.0e-2 | 1.5e-3 | 9.5e-5 | 3.5e-6 |
| 0.3 | 8 | 7.0e-2 | 3.4e-3 | 2.3e-4 | 8.8e-6 |
| 0.3 | 50 | 4.3e-2 | 1.9e-2 | 4.9e-3 | 4.0e-4 |
| 2 | 0 | 8.3e-4 | 8.6e-7 | 2.8e-8 | 8.9e-10 |
| 2 | 8 | 1.1e-2 | 4.6e-5 | 1.7e-6 | 5.4e-8 |
| 2 | 50 | 1.5e-2 | 5.1e-3 | 9.7e-4 | 6.7e-5 |
| 16 | 0 | 3.9e-5 | 3.1e-8 | 1.0e-9 | 3.2e-11 |
| 16 | 8 | 2.2e-3 | 7.1e-6 | 2.6e-7 | 8.3e-9 |
| 16 | 50 | 5.2e-3 | 1.7e-3 | 3.0e-4 | 2.1e-5 |

So a fixed number of terms means accuracies more than 10,000 times apart across rows (series 10: 1e-9 at b 16,
psi 0; 5e-3 at b 0.3, psi 50), and any approximate default or fixed K rests on these rows, not the screen's
columns (Open call 3). The normal is worse than its TV column suggests at small shapes: its mean is off by 150
standard errors at b 0.2 in 1e5 draws (series-check.R), and it was dropped (below).

## The exact accept-or-keep route

At the start of each sweep, before the trees (dec-B586), each row proposes omega' from series K re-matched at
its current psi (K by the rule below), and keeps omega' with probability
min(1, f(omega') q(omega) / (f(omega) q(omega'))), otherwise keeps omega.

- f is PG(b, psi), an exact exponential tilt of PG(b, 0): f's value at the current omega is cached per row as
  log PG(b, 0)(omega) and tilted to the current psi analytically. q is not a tilt: re-matching the remainder at
  psi changes its shape, and matching it once at psi = 0 and tilting is poor (probe tiltq.R: at b = 64,
  psi = 20, acceptance 0.66 tilted against 0.9993 re-matched; at b = 1e3, psi = 20, 0.08 against 0.9998). So q is
  evaluated at both omega' and omega each sweep.
- Density (S1). An error d in the ratio moves the chain's distribution by about d, so the decision is made
  certified: f and q each come with a true error bound, the step accepts or rejects only when the uniform clears
  the ratio by more than the bound, and otherwise refines and decides again. Both f and q are certified at every
  (b, psi) the route can reach, b from 0.2 up and |psi| to 150 at least, with the bound shown to hold against
  independent references at every pin. What the evidence supports, by region:
  - b >= 12, f and q: the trapezoid rule on the vertical line through the saddle point (17 to 21 nodes at
    b >= 64, 48 at b 16; correct to |psi| 150 in the review's chains). Its a priori bound follows from
    analyticity in a strip around the line and bounded truncated tails (nbinom-large-shape-pg.md).
  - b < 12, f_0 at small x: the Polson-Scott-Windle left series, where it converges; an alternating series, its
    error is bounded by the first omitted term.
  - b < 12, f elsewhere and q everywhere: a contour through the saddle point with an a priori bound. The scratch
    build's Talbot contour shifted to the first singularity gave log f within 2.4e-10 of three independent
    references (line, left series, parabolic contour, agreeing to 1e-13) at b 0.2, but its error estimate is not
    a bound: for q at b 8, |psi| 100 to 150 it was wrong by 100 to 200 in log while estimating 0.04 to 1, and at
    b 3.5 to 8 it returned -inf at about half the points with |psi| >= 50. A parabolic contour with the
    trapezoid rule (geometric convergence in a strip, with a known bound) matched the references to 1e-13 in the
    review. S1 builds the contour whose bound it can derive and pins it, not an after-the-fact estimate.
  - Refinement adds nodes, never raises the rounding floor: the first scratch build refined by doubling the
    contour's nodes, and at 128 and 256 nodes log f was off by 4 to 67 (at its starting 32 nodes, by at most
    1.5e-7 over the bulk), which stranded the chain at shapes below 1. Adding 8 nodes a step fixed it. S1's
    refinement is pinned not to grow the error at any pin.
  - A non-finite value is never decided on: it escalates (more nodes, then the other rule for that region),
    and S1's pins require a finite value with a finite bound at every reachable point.
- The decision (S2). The step decides only once the bound clears the gap between the uniform and the ratio. It
  never decides on a non-finite or unresolved ratio: it refines until the decision resolves. The one exception is
  the rounding floor, where the gap is below about 1e-12 relative, an event of about that probability, the limit
  every floating-point sampler meets; it is tallied and must be zero in every gate cell. A fresh exact draw in
  place of an unresolved decision is not taken: it does not keep the distribution either (whether a decision
  resolves depends on the current value), and no fresh exact draw exists at fractional shapes. In the scratch
  build a NaN ratio (q -inf at both points) rejected silently; at b 4 that happens to be symmetric and kept the
  distribution, at b 2.5 and psi 150 (acceptance 0.011, the mean off by 45 standard errors) taking unresolved
  decisions did not.
- The K rule (S2). The proposal's weight f/q is bounded for every b > 0, psi and K >= 1, but not uniformly: its
  right-tail limit is exp(b D(psi, K)), D > 0 (scratch/apx-check tail-ratio.R). K is chosen per row so that
  b D(psi, K) <= 0.01, which keeps each step uniformly ergodic with a constant the design controls. The scratch
  build's table stopped at |psi| 64 and capped K at 48, where the bound fails above b of about 140 at |psi| 64,
  7 at 100 and 0.6 at 150 (scratch/apx-rev/probes kbound.R). S2 computes K from D's closed form at every |psi|:
  for b of 100 and more K is about 50 to 60 at |psi| 64, 80 to 90 at 100 and 110 to 120 at 150
  (scratch/apx-amend/probes kneed.R), K + 1 gammas at about 41 ns each, so 5 us at the extreme. Such |psi| arise
  in fits near separation, not ordinary ones.
- Rejection from the same proposal with bound M would give independent exact draws, which would also serve
  multinomial and fractional-shape rows a setter changed. It needs a certified M. The tail limit is not one: on a
  grid, sup f/q exceeds exp(b D) at b <= 1 (1.10 against 1.0002 at b = 0.3, psi = 8, about one sd below the
  mean), though not at b >= 16 where checked (probe supcheck.R). Recorded as research, not built.

Where the route keeps the posterior (scratch/apx-check report, with exact finite-state probes):

- Within a sweep it holds: rows are independent given the trees, and a row's step reads only that row's omega
  and psi.
- Estimated-r nbinom holds under this design's rule: on a sweep where r moved, every row draws fresh and
  exactly; where r stayed, the step applies. r's draw never reads omega, so on the event that r stayed the old
  omega is a draw from its conditional at that r. Checked exactly: the rule's error about 1e-16, always applying
  the step 0.22 in TV (pcg-exact.R). The same holds for any future r step that leaves r's conditional invariant
  without reading omega. Its cost: r moved on 38 to 49 percent of sweeps at 200 rows and 0.3 to 0.6 percent at
  2,000 in the sizing, and each move draws every row by the sum, so at 200 rows and counts near 1,000 the route
  costs 48 us a row and is 3 times the sum's speed, not 30. The remedy is an r step that reads omega (a
  Metropolis step on r given omega, its ratio from the same density), which never forces fresh draws; it is also
  what real dispersion needs (below), and is recorded for that plan, not built here.
- Multinomial does not hold: [`MultinomialForestCombiner::drawForestGlue`](../../src/bartcore/combiner.hpp) draws
  each category's omega as a temporary latent against the current margins just before that category's forest,
  and there is no persistent joint omega to keep. The route there gives TV 0.023 in the check's exact toy.
  Multinomial keeps the exact sum under every exact setting; an approximate setting, being a fresh draw, can
  serve it.
- The per-row density cache is scratch, recomputed from omega; it is invalidated on restoreLatents, nbinom's
  cancel-restore, installShape, restoreShape, setData, any change of b and any cold start.

### Setters and marked rows

The latents-first plan's constraint: an accept-or-keep latent is exact for an outer step that integrates it out
only where it is a fresh draw from its current conditional, so every row whose conditional a setter or a restore
moved since the last sweep is drawn exactly at the start of the next sweep, and its mark is cleared only after
that draw. The rows marked: every row for setOffset, setForestBasis, setResponse, setData, setModel's shape and
setState; the moved rows for setPredictor; the changed rows for setWeights and setActiveRows; rows changed inside
the per-sweep hook included. The setters themselves draw nothing.

So the separate exact path this plan had for setters becomes one branch at the sweep's start (S2): a marked row
at a whole-number shape is drawn by the sum, an unmarked one takes the step. The count of exact draws is the
same as drawing at the setter; only the point in the sweep, and so which draws a seed produces, differs. The
order changes no sizing number: one latent update per row per sweep either way, so the cost per sweep, the
acceptance and the effective draws per sweep stand. One consequence the sizing did not price: a larger sampler
that changes the offset every sweep marks every row every sweep, so its fit draws by the sum throughout and the
route's speed does not reach it; an approximate setting does.

At a fractional shape a marked row has no exact fresh draw ([Real shapes](#real-shapes)).

### Real shapes

The route works at any real b > 0 given the density, which S1 supplies at every x > 0. That is the exact draw
dec-B568 waits on, within sweeps and from a starting value. A fresh exact draw at a fractional shape is another
matter: none is published below 1, and WPS's alternate sampler above 1 rests on Conjecture 8. So, under every
exact setting:

- Real weights on logistic, binomial and hazard.logistic are taken at creation. A change between sweeps that
  would mark a fractional-shape row (setOffset, setPredictor, setForestBasis, setResponse, setWeights,
  setActiveRows, setModel's shape) is refused by name, as dec-B568 takes real weights only where the draw is
  exact; the message names the approximate setting as the way to continue. A starting point (creation, setData,
  setState) is served by the route from the installed value, as a cold start is: a start needs no exact draw,
  only a chain that converges from it. Rejection from the proposal, once its bound is certified, would lift the
  refusals.
- Integer plus fraction (PG(n + f) is PG(n) plus an independent PG(f)) would send only shapes below 1 to the
  route, where the review's checks passed. Not taken: it is exact only with the fraction's latent kept as its own
  chain state (saved states, restores, copies); it buys no speed (the route at fractional shapes costs about
  7 us a row with or without it); and the density must be certified below 1 regardless, which also covers
  1 to 12.
- Real negative binomial dispersion (TODO negbin-real-dispersion) needs a fresh draw whenever r moves, so it
  needs either that rejection draw or an r step that reads omega; its own plan decides. Not in this build.

## The sizing, as run

Run 2026-11-03 on a scratch build carrying every survivor behind an environment switch (970 changed lines,
never lands; its density as described above, with the defects found). Timing: an x86-64 Linux box, one core,
seconds per sweep, minimum over 5 alternating rounds; median round-to-round spread 0.4 percent, 95th percentile
2.2; load at or below 1.32. Effective draws: the laptop, 2 chains of 2,500 kept sweeps, rank-normalized bulk ESS
of f(x_i) and of the probability and tail ESS of a 5 percent quantile (Vehtari et al. 2021), on 8 cells; ten
seeds of one arm in three cells for the spread. Compute 1.6 core-hours.

| family | rows | trees | shape per row | arms |
|---|---|---|---|---|
| logistic with binomial rows | 200, 2,000 | 75, 200 | 16, 64, 256, 1e3, 1e4, 1e5 (the sum to 1e4 at 2,000 rows), and SIDS-like spread (lognormal, median 2,000) | sum, route, series 1, 5, 10, gamma, normal |
| logistic with real weights | 2,000 | 75 | uniform on (0.2, 0.8) and on (1.2, 3) | route, series 1, 5, 10, 20, gamma, normal |
| nbinom, fixed r = 1e3 | 200, 2,000 | 75 | counts mean 10 and 1e3 | sum, route, series 1, 5, 10, gamma, normal |
| nbinom, estimated r | 200, 2,000 | 75 | counts mean 20 and 1e3 | sum, route (fresh on r moves), series 1, 5, 10, gamma, normal |
| multinomial, 3 categories | 200 | 75 | 1e3 trials | sum, series 1, 5, 10, gamma, normal |

Costs measured:

- The sum: about 0.13 us per trial per row.
- The route: 3.1 to 3.9 us a row at shapes of 64 and up (about two density passes a step at 17 to 21 nodes,
  the Newton solve and K + 1 gammas); 8.8 at 16 (47.5 nodes on the line); about 7.0 at fractional weights (34
  nodes: the left series and a 32-node contour); 6.2 to 7.7 at estimated-r nbinom with counts near 20; 48 at
  200 rows with counts near 1,000, where r moved on 38 percent of sweeps. Acceptance 0.9998 to 1.0000 on every
  logistic and nbinom cell, 0.9936 and 0.9977 at the two real-weight cells.
- Everything else a sweep does: 0.13 to 0.40 us a row at 75 trees, 0.5 to 1.0 at 200 trees and in multinomial,
  not the roughly 2 us the plan had assumed. So under the route the latent is most of the sweep, and the
  approximations beat it by 3 to 7 times, not the 2.2 to 2.9 first estimated.
- Route against sum: the sum 3.2 to 4.0 times faster at 16 trials; the route 1.9 to 2.1 times faster at 64, 31
  to 38 at 1,000, 300 to 350 at 10,000. The crossover lies between 16 and 64 and was not timed in between.

Effective draws per sweep, smallest of the three quantities (always the fitted value's bulk ESS), min to max over
arms: logit 200/75/16 0.30 to 0.34; 200/75/1000 0.010 to 0.016; 200/200/1e4 0.039 to 0.049; 2000/75/256 0.0023
to 0.0047; SIDS-like at 2,000 rows 0.0010 to 0.0013; multinomial 0.013 to 0.018; estimated-r nbinom 200/75/20
0.26 to 0.28; real weights low 0.55 to 0.57. They do not differ between arms beyond Monte Carlo error, so the
rule's ratios are ratios of time per sweep. The limit of that finding: a single run separates two arms only when
they differ by 1.57 to 1.73 times on the smallest quantity (the ten-seed spread), far above the 1.25 bar; the
support is theory (the latent touches the trees' mixing only through its own error) and the two high-ESS cells,
where the arms sit within 14 and 5 percent in no consistent order by accuracy
(scratch/apx-rev/probes ess-power.R).

The rule, stated before the run: a survivor is built if, in at least one cell, its effective draws per second
are at least 1.25 times those of every survivor at least as accurate (exact, series 20, 10, 5, 1, gamma,
normal). Applied per survivor (the review recomputed it; median-of-rounds times change no verdict):

| survivor | largest ratio against every better survivor | cells passing | built |
|---|---|---|---|
| sum | 4.00 against the route (2,000 rows, 75 trees, 16 trials); 1.80 to 2.24 at nbinom counts near 20 | 6 | yes (kept regardless: multinomial, marked rows, the augmentation door) |
| route | 2.0 at 64 trials; 31 to 38 at 1e3; 3,170 at 1e5; the only exact draw at real weights | 28 | yes |
| series 20 | 5.84 to 5.86 against the route (real weights) | 2 of 2 | yes, for real weights |
| series 10 | 2.2 to 6.9 against the exact draw on logistic; 58 at nbinom estimated r, 200 rows, 1e3; 173 against the sum in multinomial | 39 of 39 | yes |
| series 5 | 1.48 (real weights); 1.33 to 1.40 at 2,000 rows, 75 trees; 1.13 to 1.22 at 200 trees | 23 | yes |
| series 1 | 1.64 (real weights high); 1.46 to 1.47 at 2,000 rows, 75 trees; 1.12 to 1.22 at 200 trees | 25 | yes |
| gamma | 1.17 | 0 | no |
| normal | 1.18 | 0 | no |

Limits of the verdict:

- Trees: the grid has 75 and 200. With fewer trees the latent is a larger share of the sweep. Scaling the other
  work by the tree count (arithmetic, not measured), at 2,000 rows the gamma would pass against series 1 below
  about 30 to 50 trees and the normal against the gamma below about 17 to 48 (drop-threshold.R). Dropping them
  holds for 75 trees and up; a timing check at 2,000 rows with 25 and 50 trees precedes S4 (step 6).
- Series 5 and 1 clear the bar by 6 to 18 percent, robust to timing noise, but only at 75 trees (they fail at
  200) and at the cost of R's gamma generator: halving its cost would bring them to about 1.26 (arithmetic).
- ESS at 2,000 rows with large shapes: 5 to 23 effective draws per 5,000 sweeps, BART itself mixing slowly
  there, so those cells rest on the timing and the equality above.
- The real-weights cell with weights 1.2 to 3 on 0/1 responses drifted outward over thousands of sweeps; the
  exact sum at weight 2 and the sum on every row entered twice drift alike, so it is the model (weight 2 is the
  duplicated-rows posterior), not the sampler (drift-rev.R). The first round, before the drift, was used.
- No Geweke check at real weights (0/1 responses have no generative model under fractional weights); moments,
  KS against the integrated density and the posterior-mean comparison (0.57027 against 0.57021 to 0.57036 for
  the approximations) stand in.
- Not run: arm64 timing, a faster gamma generator, real-r nbinom.
- The scratch route is biased in a band no sizing cell entered (b about 2 to 12 with |psi| of 80 or more: the
  mean 0.6 to 0.9 percent high at b 8 and 11.5, z 11.5 to 22; route-onset.R, route-b2.R, re-run 2026-11-05),
  from its uncertified q. Across 63 cells outside the band (b 0.2 to 1,000, psi 0 to 150) every |z| was 2.62 or
  less and the smallest KS p 0.015 (route-check.R). Its costs stand; its density is replaced by S1's.

## What is built

Tiers by [Process by risk](README.md#process-by-risk). S4 and S5's approximate values as Open call 2 rules; the
setting's shape as Open call 3 rules (the table assumes the recommendations).

| slice | content | tier | lines | agent-hours |
|---|---|---|---|---|
| S1 | the certified PG density, f and q, at every reachable (b, psi): the line from 12, the left series and a bounded contour below; true a priori bounds; refinement that never raises the floor; non-finite values escalated; tests/cpp pins on a grid of b 0.2 to 1e5 and \|psi\| 0 to 150 against independent references, the bound covering the error at every pin | engine, neutral | 650 | 9 |
| S2 | the route at the sweep's start for logistic, binomial, hazard.logistic, fixed- and estimated-r nbinom; marked rows drawn by the sum; the moved-r rule; cache and invalidation; the K rule from D's closed form at every \|psi\|; decisions only on resolved ratios, with tallies; the setting's engine side | engine, shifting where selected | 800 | 13 |
| S3 | real weights under the route on logistic, binomial, hazard.logistic; refusals of between-sweep setters that would mark a fractional row; starting points served by the route | engine plus R surface | 300 | 7 |
| S4 | the series as an approximate setting for all PG families including multinomial, its terms chosen per row from a TV target by a table computed once by inversion (0 terms, the gamma, at the loose end); the normal only if step 6's check passes it | engine, posterior-changing where selected | 300 | 6 |
| S5 | the setting at the doors users reach, help, NEWS | R surface | 300 | 5 |

The exact pieces alone (S1 to S3, and S5 without the approximate values): about 2,000 lines, 33 agent-hours.

Gates for S2 and S3 (exactness), each with the exact sum as control, passing first:
- joint-distribution check (Geweke 2004) as [geweke-mc.R](../../benchmarks/R/geweke-mc.R) runs it, on logistic
  with binomial rows, fixed- and estimated-r nbinom and real weights; broken builds that must fail: the step
  applied after r moved, a stale cache after psi changes, q tilted from psi = 0, a marked row taking the step in
  place of an exact draw;
- two-sample KS of the route's omega at fixed (b, psi) against the exact sum (an independent reference at
  fractional b), and moments to the third cumulant, on b in {0.2, 0.5, 1, 2, 2.5, 4, 8, 11.5, 12.5, 16, 64, 1e3}
  by |psi| in {0, 8, 50, 80, 100, 120, 150}, with zero decisions at the rounding floor in every cell; a broken
  build with q on the scratch build's uncertified contour must fail at b 8, |psi| 100 and b 2.5, |psi| 150;
- the exact single-tree posterior of [logistic-reference.R](../../benchmarks/R/logistic-reference.R) part 1
  with rows of b trials;
- simulation-based calibration through [sbc.R](../../benchmarks/R/sbc.R)'s logistic machinery with binomial
  rows, its exact arm passing first;
- the engine tier's battery: tests/cpp, full tinytest, sanitizers on tests/cpp and the R-loaded path, R CMD
  check, the arm64 snapshots, the equivalence baseline bitwise for fits whose shapes are whole numbers of 64 or
  less;
- mutants: the review runs the mutation battery's anchors plus mutants of the accept comparison, the K rule (a
  cap at 48), a non-finite ratio decided on, the cache invalidation, the moved-r branch and the marked-row
  branch.

S4's gates: its TV per row against inversion in tests/cpp, the chosen K meeting the target at every grid cell of
the table above and between; the exact gates stay on exact settings.

## The defaults study, after the build (outline)

It measures, on fits shaped like users' and on real covariates with real trials, each built setting's effective
draws per second and, for the approximate ones, their effect on what users read: coverage and width of 90 and
95 percent intervals for probabilities, rates and contrasts, RMSE and bias, against exact posteriors (the
one-parameter toy by quadrature, the single-tree enumeration) and semi-synthetic truths, at equal running time.
Its fractional-weight cells use fractional counts with proportions, or weights below 1, never 0/1 responses at
weights above 1: weight 2 on a 0/1 row is the posterior of that row entered twice, which drifts toward
separation under any sampler. It covers small shapes and |psi| to 50 at least, where the approximations are
least accurate, and tree counts below 75. Choosing among exact settings changes no posterior and is a speed
finding. Making an approximate setting a default is a change users' results carry, and comes as a short research
report under the studies standard (dec-B559): at most 5 pages, readable alone, two referees first.

## Decision

### Open call 1: the default until the defaults study

Until the study, which exact way should a fit draw its hidden values when the user sets nothing?

Background. Each row of a logistic, binomial, logistic-hazard or negative binomial fit carries a hidden value
redrawn every sweep. It can be drawn exactly two ways: the sum, one small draw per trial, whose cost grows with
the trials; or the accept-or-keep route, which proposes from a close approximation and keeps or refuses against
the exact density, at a cost nearly flat in the trials. Both give the same posterior; they differ in speed and in
which draws a seed produces. 0.9-34 had none of these families, so no release's draws are owed.

What was measured, and how. Whole fits on a scratch copy of the engine, timed on one core of an x86 Linux machine
(the fastest of 5 rounds; rounds differed by 0.4 percent in the median), at 200 and 2,000 rows, 75 and 200 trees
and 16 to 100,000 trials a row. Effective draws per sweep were the same for both ways within noise, so speed per
sweep is the comparison. The sum costs about 0.13 microseconds per trial; the route 3 to 4 microseconds a row
from 64 trials up, about 9 at 16 trials and about 7 at fractional weights. The sum is 3 to 4 times faster at 16
trials; the route twice as fast at 64, 31 to 38 times at 1,000 and about 300 times at 10,000. Limits: nothing
between 16 and 64 trials was timed, nothing on arm64; the route's density in the scratch copy was wrong at some
extreme shapes, which the build's density replaces, and its cost was measured on the scratch one.

- (a) Per row: the route when the row's shape (its trials, or count plus dispersion for the negative binomial)
  is above 64 or not a whole number, the sum otherwise. A fit of 1,000 rows at 10,000 trials over 5,000 sweeps
  goes from about 1.8 core-hours to under a minute, and a fit whose shapes are all whole numbers of 64 or less
  draws as it does today. Costs: S2's gates before it can be the default; between about 16 and 64 trials the sum
  is kept where the route may be faster, by at most twice; a negative binomial fit with estimated dispersion and
  few rows gains less (at 200 rows the dispersion moves on 38 to 49 percent of sweeps, each move redraws every
  row by the sum, and the route is 3 times the sum's speed at counts near 1,000 rather than 30); and a larger
  sampler that changes the offset every sweep gets the sum on every row.
- (b) The sum wherever it can draw, the route only by setting and for fractional shapes. Nothing moves; fits
  with hundreds of trials or more stay 30 to 350 times slower unless the user finds the setting.
- (c) The approximation above a shape. The 10-term series measured 3.0 to 6.9 times faster than the route at 64
  trials and up (most at 2,000 rows and 75 trees, least at 200 rows and 200 trees), but it is not exact, and it
  would be a default before the study that is to decide one.

Recommendation: (a). It is exact, faster wherever it applies, and asks nothing of users; the threshold sits where
the route was already twice as fast on the measured machine, erring toward the sum, which leaves those fits'
draws as they are.

### Open call 2: what ships in 1.0-0

Which of the samplers built here should ship in 1.0-0: the exact ones only, or also the approximate ones as
settings a user must ask for?

Background. Your ruling for 1.0 reads "Every draw in 1.0 stays exact: an approximate Polya-Gamma draw ... is
declined"; the accept-or-keep route is exact under your later ruling. The approximate samplers would be settings
only: no fit draws approximately unless its user asks, until the defaults study. Adding them after 1.0-0 is
possible as new values of the same setting.

What was measured, and how. The same timed fits. The 10-term approximation against the fastest exact way: 3.0 to
6.9 times faster than the route at 64 trials and up; 58 times faster for a negative binomial with estimated
dispersion at 200 rows and counts near 1,000, where the route is held back by the dispersion moving; 173 times
faster than the sum for a 3-category multinomial at 1,000 trials, which has no exact route (its hidden values are
drawn afresh before each category's trees, so there is nothing to keep or refuse). Its accuracy, as the
distance of each draw's distribution from the exact one (total variation): 1e-9 to 3e-7 at shapes from 16 up and
log-odds up to 8, but up to 3e-4 at log-odds 50 and up to 5e-3 at a shape of 0.3 and log-odds 50. Not measured:
what it does to intervals and estimates, which is the defaults study's job.

- (a) The exact pieces in 1.0-0 if their gates pass in time, the approximate settings after. About 33
  agent-hours. Keeps the 1.0 clause as written. Multinomial fits with many trials, estimated-dispersion negative
  binomial fits with few rows and embedded fits that change the offset every sweep have no faster draw at
  release.
- (b) Everything in 1.0-0, the approximate settings opt-in and documented as approximate, revising the 1.0
  clause. About 40 agent-hours, 7 more. Every default draw stays exact; the cases above get their only fast path
  at release. The settings reach users before the study has measured their effect on intervals, so their help
  can state only their accuracy per draw.
- (c) Nothing from this arc in 1.0-0. Large-trial fits stay 30 to 350 times slower than they need be, and
  fractional weights stay refused on these families.

Recommendation: (b). No default draw changes, the cost over (a) is about 7 agent-hours, and the fits that the
exact route cannot speed up, among them every multinomial fit, are where the measured gain is largest. The plan's
earlier recommendation was (a); the measured multinomial and estimated-dispersion ratios are what changed it.

### Open call 3: the setting's shape

How does a user choose the sampler and, for the approximation, its accuracy?

Background. Elsewhere: BayesLogit exposes separate functions per method (its default switches between them
silently); the Python polyagamma package takes method = with "devroye", "alternate", "saddle" and "gamma", None
meaning its hybrid; R's glm and the other model fitters users know expose nothing of the kind; base R's numerical
routines take a tolerance (integrate's rel.tol, optim's reltol). In dbarts, how the sampler draws lives in the
control settings, beside the move probabilities, reachable from bart, dbarts and xbart, and sigest already takes
either a rule's name or a number. Under the package's rule for absent values, NULL is the default rule.

What was measured, and how. The approximation at 20, 10, 5 and 1 exact terms each passed the build rule, each
step down about 1.2 to 1.7 times faster per sweep at 75 trees, less at 200; from 20 terms to the gamma alone,
about 4 times at fractional weights. Its accuracy at a fixed number of terms varies more than 10,000-fold across
rows, measured by inversion of the exact density: 10 terms give 1e-9 at shape 16 with log-odds 0 and 5e-3 at
shape 0.3 with log-odds 50. Limits: accuracy is per draw; nothing yet says what accuracy intervals need.

- (a) One setting, latent.sampler = NULL (the default rule), "sum", "metropolis" or "series", the series' terms
  chosen inside from one fixed accuracy target. The simplest to read; it gives up the speed range, so neither a
  user nor the defaults study can trade accuracy for speed.
- (b) As (a), plus the number of terms, latent.terms = 10. Keeps the range, but the same number means accuracies
  10,000 times apart in different rows and fits, so a finding at one shape, the study's included, does not carry
  to another.
- (c) One setting taking a method's name or a number: "sum", "metropolis", "series" (the series at a default
  accuracy, 1e-6 total variation per draw, between the 5- and 10-term series at shape 16), or a number, the series at
  that total variation per draw, its terms chosen per row; 0 terms, the gamma alone, is its loose end. As sigest
  takes a name or a number. Keeps the range with one meaning in every row; costs about 100 lines and an
  agent-hour in S4 for the table of terms and its pins, and rows at small shapes with large log-odds take more
  terms (20 terms reach 4e-4 at shape 0.3, log-odds 50), so their draws cost more there.

Recommendation: (c). One setting, as in (a), named for what it chooses as polyagamma's method is, with values
naming the method; and a number for the accuracy means the same in every row, which is what lets the defaults study's
finding become a default.

## Constraints

- No package code from the sizing lands; slices follow their tier.
- Any default keeps every draw of a fit whose shapes are all whole numbers of 64 or less as the build before it.
- Draws come only from the chain's RNG; the cache is scratch, never state.
- S2 runs after the latents-first plan's engine slices and meets its marked-row constraint.
- Nothing from BayesLogit's source (GPL 3; dbarts is GPL 2 or later).
- Out of scope: the probit binomial, real negative binomial dispersion and the r step that reads omega (recorded
  for that plan), the dec-B202 cap.

## Steps

1. Done: the scratch sizing build (970 lines), reviewed before its run.
2. Done: the sizing and its verdict per survivor under the stated rule; the independent review (sound with
   corrections) and this amendment.
3. Open calls to the maintainer with the measured ratios.
4. The latents-first engine slices land.
5. S1, then S2 and S3, each through its tier.
6. Before S4: time 2,000 rows at 25 and 50 trees on the scratch build (about 0.2 core-hours); the normal is
   built in S4 if it passes the rule there (about 30 lines).
7. S4 if Open call 2 includes it in this window; S5.
8. The defaults study, planned once the build lands.

## Verification

- The screen's numbers reproduce from scratch/apx-design/probes (laterr.R, screen.R, tiltq.R, supcheck.R),
  scratch/apx-check/probes (tail-ratio.R, pcg-exact.R) and scratch/apx-rev/probes (tv-small.R, series-check.R).
- The sizing's verdict follows from scratch/apx-size/results (tables.txt, analysis.rds) by the rule above,
  mechanically (scratch/apx-rev/probes rule-recompute.R).
- The scratch route's band and the first build's defect reproduce on the scratch library from
  scratch/apx-rev/probes (route-onset.R, route-b2.R, route-b16.R, q-b4.R, q-b8.R, old-refine.R); re-run
  2026-11-05, matching the review.
- S1's pins hold the bound above the error at every pin; S2's broken builds fail their gates; a gate passing one
  has no power and stops the slice.

## Taken and not taken

- From the critique (September design's revision): the exact-moment remainder and the re-ranking; the
  unbounded weight of the gamma proposal; certified meaning; exact draws at swaps.
- From the independent check: estimated-r nbinom holds under the moved-then-fresh rule; multinomial excluded;
  exact draws for rows a setter changed; the cache's invalidation list; f and q both at the current psi; the
  weight bounded but not uniformly, so K grows with b and |psi|; the density's error as a bias floor and the
  certified decision; real shapes (none exact below 1, the alternate's conjecture); the setOffset and
  setPredictor question, now settled by dec-B586 and the marked rows.
- Added by this design: rejection from the proposal as an independent exact draw, research; the
  fractional-setter refusal that follows.
- From the sizing: the verdict per survivor; the costs; the integer-plus-fraction split, considered and not
  taken ([Real shapes](#real-shapes)); its diagnosis of the first build's failure as the contour itself and its
  reading of the acceptance drop as a weak proposal, both replaced by the review's, which ran the decisive probes.

## Review amendments

The sizing's independent review (2026-11-05, sound with corrections) proposed 15 amendments:

1. S1 certifies q as well as f at every reachable (b, psi): applied, [The exact accept-or-keep route](#the-exact-accept-or-keep-route)
   (Density) and S1.
2. A true error bound, refinement never raising the rounding floor: applied, the same.
3. Never decide on a non-finite or unresolved ratio: applied, the route (The decision) and S2; the fallback to a
   fresh exact draw declined (not invariant either; none exists at fractional shapes).
4. The K rule past |psi| 64: applied, the route (The K rule) and S2, with K at |psi| 64 to 150 from a new probe.
5. Gate cells at b 1 to 12 and |psi| to 150: applied, S2 and S3 gates, with a broken build on the uncertified q.
6. The screen's TV table gains small-b and |psi| 50 rows: applied, [The dominance screen](#the-dominance-screen); it
   informs Open call 3 (c).
7. Open call 1 (c)'s ratios and the other work per row: applied, [The sizing, as run](#the-sizing-as-run) and Open
   calls 1 and 2.
8. The crossover between 16 and 64, the threshold of 64 standing: applied, the sizing and Open call 1.
9. Estimated-r nbinom's cap and the r step that reads omega as its remedy: applied, the route (where it keeps
   the posterior), Open call 1 and Constraints; the remedy recorded, not built.
10. Gamma and normal dropped for 75 trees and up only: applied, the sizing's limits, with the timing check of
    step 6; under Open call 3 (c) the gamma is the tolerance's loose end at no extra cost.
11. Series 5 and 1 pass only at 75 trees and R's gamma cost: applied, the sizing's limits.
12. Equal ESS per sweep rests on theory and two high-ESS cells: applied, the sizing.
13. The integer-plus-fraction split keeps the fraction's latent as chain state: considered, split not taken
    ([Real shapes](#real-shapes)); fractional setters stay refused.
14. No 0/1 responses at weights above 1 in the defaults study: applied, the defaults study.
15. The first build's defect is the refinement doubling the nodes; the plan prescribes the saddle line: applied,
    the route (Density), the sizing and Taken and not taken.

## Not verified

- The route's cost under S1's certified density: the scratch density chose its nodes from after-the-fact
  estimates, and a true bound may take more.
- The contour's a priori bound below shape 12: not yet derived (S1's work); the parabolic contour's accuracy was
  measured against references, not bounded.
- The crossover between 16 and 64 trials, and every timing on arm64.
- K past |psi| 64 is from the closed form of D (kneed.R), not from chains at those K.
- Gamma and normal below 75 trees: arithmetic only, until step 6.
- That effective draws per sweep match across arms beyond the two high-ESS cells and theory.
- The cost of the table of terms for Open call 3 (c), estimated, not built.
- The refusal of between-sweep setters at fractional rows and the treatment of starting points are this
  amendment's reading of the latents-first constraint, not checked by a probe.
- The saddle-point sampler's cost (WPS Table 1) and BayesLogit's switching points (the critique's reading).
- The far-tail behavior of f/q at b <= 1 on finer grids: supcheck.R's largest ratios at b = 0.3, psi = 0 sit
  where the inversion loses accuracy and are not read.
