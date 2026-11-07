# approximate-draws: Polya-Gamma latents whose exact cost grows with their shape

Status: DESIGN, re-planned 2026-11-03 under dec-B583, dec-B584 and dec-B585; amended 2026-11-05 from the sizing
study and its independent review, and for dec-B586 (augmentation variables drawn at the start of each sweep);
amended 2026-11-06 for the proof dec-B591 asked for and for the rulings dec-B590 to dec-B593, dec-B606 and
dec-B607; amended again 2026-11-06 from a blind critique of that amendment. The three open calls are ruled. The
build starts with S1, whose two gates (the close-call floor and the timing) decide whether S2 proceeds as planned;
the per-row threshold's dependence on |psi| comes back to the maintainer with S1's timings. The defaults come from
a study run before 1.0-0, right after the build (dec-B592).

agent: opus for each engine slice; one opus review that runs mutants per engine slice; sonnet for the R-surface
slice and its review.
rng: the density, neutral; the exact accept-or-keep route, shifting where it is selected; real-shape weights, a
new capability; approximate settings, posterior-changing where selected. Until the defaults study a fit draws by
dec-B590's per-row rule, its threshold set from S1's timings once ruled and never below 64
([The per-row rule](#the-per-row-rule)); a fit whose shapes are all whole numbers of 64 or less draws as the build
before it. The latents-first engine slices S2 waits on landed on main on 2026-11-05; the branch rebases onto main
before S1. S3 lands after frequency-weights.
window: everything built here ships in 1.0-0, the approximate settings opt-in (dec-B592); the defaults study runs
before 1.0-0, right after the build.
budget: the sizing ran at 970 lines in a scratch build and 1.6 core-hours. The build, with every survivor, about
4,950 written lines over five slices (the proof's source and scripts, moved into the tree, not counted) and about
75 agent-hours, with about 13 core-hours of local gates and about 20 of scheduled CI ([What is built](#what-is-built));
before the proof, 2,350 lines and 40 agent-hours.

## Summary for the maintainer

The problem. Each row of a logistic, binomial, logistic-hazard or negative binomial fit carries a hidden value
that is redrawn every sweep. Today it is drawn exactly as one small draw per trial, 130 to 160 nanoseconds each.
Everything else a sweep does costs 0.13 to 1 microsecond a row at 75 to 200 trees, so from a few dozen trials a
row the hidden values are most of the sweep, and at 10,000 trials a fit is hundreds of times slower than it need
be. Real-valued weights and shapes have no exact draw today and are refused.

What you ruled. An accept-or-keep step whose long-run distribution is correct is exact, and, being our own
construction, it is proven mathematically by independent agents before it is built. Samplers beaten on both
accuracy and cost are dropped and the rest sized by effective draws per second. Until a defaults study, a row
draws by the accept-or-keep route above a shape threshold or at a fractional shape, and by the exact sum
otherwise. Everything built ships in 1.0-0, the approximate sampler only when a user asks, and the study runs
before 1.0-0. A user names the method, and the approximate one's accuracy is a separate tolerance. The step's
close calls are kept, counted, and warned about, never silent. At fractional weights every change between sweeps
is taken, and the help says which larger samplers stay exact there, with an example.

The candidates. Today's exact sum; the exact accept-or-keep route, which proposes from a close approximation and
keeps or refuses against the exact density, at a cost that does not grow with the trials; an approximation that
draws the first few terms of the hidden value's exact series and one gamma for the rest; and a plain normal. The
saddle-point sampler and BayesLogit's truncated sums were dropped as beaten on both counts.

What the sizing measured. On a scratch copy of the engine, each candidate was timed on 39 fits shaped like users'
(200 and 2,000 rows, 75 and 200 trees, 16 to 100,000 trials, negative binomial counts with fixed and estimated
shape, fractional weights, a 3-category multinomial) on one core of an x86 Linux machine, and its effective draws
per sweep measured on 8 of them. Effective draws per sweep did not differ between candidates beyond noise, so
speed per sweep decides. The approximation with 10 exact terms was 3 to 7 times faster than the route at 64 trials
and up, and 173 times faster than the sum in the multinomial fit, which has no exact route. Built: the sum, the
route and the approximation; dropped: the gamma alone and the normal, for 75 trees and up. The route's costs there
(3 to 9 microseconds a row) were measured with a density whose error was estimated, not bounded.

What the proof found. The route is exact under the standard every sampler in the package meets. But a density
with a proven error bound costs far more than the sizing's. With the bound as first specified, a step's density
work alone costs about 30 to 110 microseconds a row at log-odds up to 8 and shapes to 10,000, 90 to 350 at log-odds 50 for shapes to
1,000, and 320 and up at log-odds 150. So the route beats the sum only from a few hundred trials at log-odds up to
8, from several thousand at log-odds 50, and never at log-odds 150. At fractional weights, where there is no sum,
a step costs 40 to 90 microseconds, and about 320 at log-odds 150, against the sizing's 7. A cheaper bound, proven
but not yet built, might bring the route under the sum from about 1,000 trials at every log-odds tested. So the
first slice builds that bound and times the route. The threshold between the sum and the route is then set from
those timings, by log-odds as well as by shape, never below your 64, and its dependence on log-odds comes back to
you with them. A close call, too close to decide
within the bound, came at most about once in 160 billion steps where tested.

What it builds. The certified density; the exact route for logistic, binomial, logistic hazard and negative
binomial; fractional weights under the route; the approximation as a setting, for every Polya-Gamma family
including multinomial; and the setting at the doors users reach.

What comes back to you: the threshold's dependence on log-odds, with the first slice's timings; that slice itself
if the route does not beat the sum by 1,000 trials at small log-odds, or if its close calls are more frequent than
its gate allows and tighter rounding does not fix that; whether the C header packages build against gains a reader
of the close-call count; and, with the setting's names, which your ruling left to the build, four smaller calls on
what a user is told or refused.

## Goal

Every Polya-Gamma family draws exactly at a cost that does not grow with the shape wherever the exact route
holds and pays; real weights reach those families where the draw is exact (dec-B568), with every between-sweep
setter taken (dec-B607); and the approximate samplers that change what a user gets per second exist as settings,
for the defaults study to choose from before 1.0-0.

## Context

Rulings, quoted, not reopened:

- dec-B583, revised by dec-B591: an accept-or-keep latent is exact if the chain's long-run distribution is the
  posterior; the maintainer, 2026-11-03: "As far as a ruling goes, exact is exact, right? It either has the
  correct distribution or it doesn't." dec-B591, 2026-11-05: "I want it validated heavily, and not just
  empirical. I want it proven mathematically. Use independent agents." It also asks for a "teensy" vignette with
  the proof's math and a plain-R comparison of the samplers "if you think your contribution to Bayesian
  statistics is significant enough".
- dec-B584, revised by dec-B585: "Maybe size that to see what we actually need for ESS / sec". A dominated
  sampler is dropped, the rest sized by effective draws per second of whole fits, only those that change what a
  user gets per second built, each chosen by a setting, until a study sets the defaults.
- dec-B586: augmentation variables, the Polya-Gamma precisions among them, are drawn at the start of each sweep,
  after the per-sweep hook and before the trees; an accept-or-keep latent "still needs one exact draw at the next
  sweep for each row a setter moved". Built by the latents-first plan, landed 2026-11-05 (dec-A263).
- dec-B589: setData keeps the ordinal cutpoints, Student-t degrees of freedom and negative binomial shape; the next
  sweep draws the latents from the kept values.
- dec-B590: until the study, the route where a row's shape is above 64 or not a whole number, the exact sum
  otherwise; "And yes, you can use your recommendation." S1's timed rule may raise the threshold by |psi|, never
  below 64 ([The per-row rule](#the-per-row-rule)); that dependence comes back to the maintainer.
- dec-B592: every sampler built ships in 1.0-0, the approximate ones opt-in, revising dec-B14's 1.0 clause and
  dec-B143's "for 1.0-0"; the defaults study runs before 1.0-0, right after the build: "Sure, we can do that."
- dec-B593: the method by name, NULL the automatic rule, and a separate tolerance, the target total variation
  per draw, read only by an approximate method; "OK, use (d)." The names are settled at the build.
- dec-B606: the route is built as the proof's option A (densities in double precision with rigorous error
  bounds, the lead term and the cached density in double-double, no new dependency); a close call keeps the row's
  old value and is counted; the count is reported on the fit and a nonzero count warns; the documentation quotes
  the floor the engine measures in S1; if S1's gate of 2.3e-11 a step fails, the build tightens the rounding or
  comes back to the maintainer before taking B'. The maintainer, 2026-11-06: "I don't want error silently
  accumulating and biasing the sampler. So I think that argues for A, within reason."
- dec-B607: at fractional shapes every between-sweep setter, and setState, is taken as at whole-number shapes,
  the rows it marks served by the step from their installed values; the help of the setters and of getLatents
  says when that is exact; the documentation carries a short example. The maintainer, 2026-11-06: "Yes, as long
  as there is a way of doing it correctly with the tools we have (or are cheap to build) and we can document it,
  that seems enough for me."
- dec-B14: "Every draw in 1.0 stays exact: an approximate Polya-Gamma draw ... is declined", revised by dec-B592
  for opt-in settings.
- dec-B568: real frequency weights where the draw is exact; the Polya-Gamma families whole numbers "until a
  real-shape draw exists". dec-B607 reads "where the draw is exact" at the setters as exact for a larger sampler
  that reads the latents.
- dec-B567 (the binomial response, shape n w per row), dec-B563 (frequency weights as shape times w on the
  Polya-Gamma families), dec-B143, dec-B15 (integer dispersion), dec-B202 (the one-million cap, a consequence to
  revisit, not in scope).

Engine: the seam [`simulatePolyaGammaShape`](../../src/bartcore/model.hpp) sums round(b) draws of
[`ext_rng_simulatePolyaGamma`](../../src/external/random.c). [`LogisticResponse::refreshLatents`](../../src/bartcore/model.hpp)
(logistic, binomial, hazard.logistic) and [`MultinomialForestCombiner::drawForestGlue`](../../src/bartcore/combiner.hpp)
sum in their own loops. [`NBResponse::refreshLatents`](../../src/bartcore/model.hpp) draws r collapsed over omega
([`NBShapePrior`](../../src/bartcore/model.hpp)), then [`NBResponse::drawOmega`](../../src/bartcore/model.hpp) at the
new r. Since latents-first the setters draw nothing and every active row's latent is drawn at the sweep's start;
no per-row mark exists (the latents-first plan's section on accept-or-keep latents: "Not built here"), so S2 builds
the marks. [`dbartsDrawLatents`](../../R/augmentation.R) is a fresh draw and keeps refusing fractional weights
(augVector's counts rule): none exists.

The proof (dec-B591): "The accept-or-keep Polya-Gamma step: where it keeps the posterior exact, and where it does
not", 2026-11-06, checked by two independent referees, the second confirming every point of its last round. Its
sections are cited below by name; its section "Changes the design needs" is taken in item by item ([The proof's
changes](#the-proofs-changes)). Its TeX source and the scripts its appendix "Numerical verification" lists land in
docs/design/approximate-draws-proof/ with S1, with a short README saying how to build the PDF (no PDF is checked
in; docs/ holds none). A design note, docs/design/approximate-draws.md, lands beside it, as
[RNG classes and their gates](README.md#rng-classes-and-their-gates) asks of the posterior-changing class (S4);
Verification cites both there. dec-B591's teensy vignette stays S5's.

Order with frequency-weights (on its branch; its interaction section with this plan): S3 lands after
frequency-weights, since it reads weights.type, the binomial row and that plan's Refusal 5; S1 and S2 may land on
either side. S3 edits that plan's interaction paragraph for dec-B607 (no refusal of setters at fractional rows)
and for the certified costs (it still quotes 3 to 9 us a row and a threshold of plain 64).

Background: [nbinom-large-shape-pg.md](nbinom-large-shape-pg.md), the exact construction designed in September
(saddle-point envelope, certified density, WPS Conjecture 15), not built. Probes: scratch/apx-design/probes
(this design), scratch/apx-crit (the critique), scratch/apx-check (the independent check), scratch/apx-size (the
sizing: report, results, probes; scratch build on branch apx-size, never lands), scratch/apx-rev/probes (the
sizing's review), scratch/apx-amend/probes (the first amendment).

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
| exact accept-or-keep route (below) | 0 | 0 | 0 | 0 | measured 3 to 9 us a row uncertified (the sizing); 30 us and up certified (the proof) | survives (exact) |
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
psi 0; 5e-3 at b 0.3, psi 50), which is why the approximate method takes a tolerance and chooses its terms per row
(dec-B593). The normal is worse than its TV column suggests at small shapes: its mean is off by 150 standard
errors at b 0.2 in 1e5 draws (series-check.R), and it was dropped (below).

## The exact accept-or-keep route

At the start of each sweep, before the trees (dec-B586), each row proposes omega' from series K re-matched at its
current psi, and keeps omega' with probability min(1, r), r = f(omega') q(omega) / (f(omega) q(omega')),
otherwise keeps omega. f is PG(b, psi); q is the density of what the sampler actually draws: its rates the doubles
fl(c_k + fl(psi^2 / 2)), the first rate rounded down so that r_1 <= d_1, and the remainder's (alpha, beta) the
doubles formed from the tail sums. Rounding r_1 down buys uniform ergodicity, not invariance: rounded up, the
weight f/q is unbounded (the proof, "The row kernel"). The proof shows the step leaves PG(b, psi) invariant for
every b > 0, psi and K >= 1, and the sweep the joint posterior ("The sweep").

The ratio (S1 forms it, S2 uses it). It is evaluated in the tilted, cancellation-free form: log f and log q
directly, each node's exponent formed from the node's offset Delta from the vertex, and the lead x t_0 + log L(t_0)
in double-double ("The tilted, cancellation-free form and its rounding floor"). The untilted form, from PG(b, 0)
and the untilted proposal, put the close-call floor at up to 1e-7 a step at the largest shapes; in plain double
the lead alone reaches 1e-10 far from the mean. q is evaluated at both points each sweep: re-matching the
remainder at psi changes its shape, and matching it once at psi = 0 and tilting is poor (probe tiltq.R: at
b = 64, psi = 20, acceptance 0.66 tilted against 0.9993 re-matched). f is evaluated fresh at the proposal only; the
current point's log f is the cached value moved to the new psi by the exact update

    log f_{b,psi'}(x) = log f_{b,psi}(x) + b [lambda(|psi'|/2) - lambda(|psi|/2)] - (x/2)(psi' - psi)(psi' + psi),

lambda = log cosh, carried in double-double (two double-double exponentials and log1p's, about ten double-double
operations). So a step makes three density evaluations and one update when the row's previous step left a valid
cache entry, and four otherwise; without the cache, always four, which the proof allows ("equivalent").

The cache (S2, conditional). It saves at most one of four evaluations, so it moves b* by at most a factor of 4/3,
and it costs the double-double update, about 56 bytes a row per chain, and a broken build and a mutant of its own.
It is built only if S1's three- and four-evaluation timings move b* materially; S1 states, before timing, the
change in b* that counts as material. As specified, if built: per row: log f at a point, in double-double, with its enclosure, its refinement level, the point
as an exact double and the (b, psi) it was computed at. It is used only when the point equals the row's omega_i
bit for bit and b is unchanged, and is then moved to the current psi; otherwise it is ignored and the step makes
four evaluations. Keying on (b, psi) alone would be wrong: setData with only the response changed, setState, a
fresh draw by the sum and a larger sampler writing the latents each replace omega_i with (b, psi) unchanged, and a
stale entry would give a wrong ratio silently. The key replaces the invalidation list this plan had (restores,
installShape, setData, cold starts): every read is a valid enclosure whatever the history. It does not hold
log q(omega). After a fresh draw (a marked row, the first sweep, a row the per-row rule moved from the sum, every
negative binomial row after r moved: about 38 percent of steps in the estimated-r cell) a step makes four. The
cache is scratch, never state.

The density (S1). One rule serves f and q at every shape and every x > 0: the trapezoid rule on the parabolic
contour through the saddle point, with an a priori bound proven in the proof ("Density bounds"): the
discretization bound from two edge lines (a few a_+ in (0, rho) and a_- > 0 tried, the line integrals bounded by
cells), truncation by each dropped node's majorant plus a monotone tail, and, for q, the tail from the combined
majorant, since at large K the per-term minima do not decay within reach. rho^2 = max(T_1 - pi^2/4, T_1/2); the
step ladder h in {2^-j, (3/4) 2^-j}, so nodes are exact; evaluation stops when the mathematical bound is below the
rounding radius or stops improving. The split at b = 12, the Talbot contour and its error estimate go. (The
scratch build's Talbot estimate was not a bound: for q at b 8, |psi| 100 to 150 it was wrong by 100 to 200 in log
while estimating 0.04 to 1.) The left series is optional at small x under its start condition; the vertical line
through the saddle is valid at every b but its cost grows like b^-2 at small b, and it is not used.

Rounding is carried rigorously in double precision (absolute plus relative radii, every operation charged, no
first-order step), as the proof's appendix "Rounding" states the rules; the lead and the cached value are in
double-double, with five double-double functions carrying stated error bounds (exp, log, log1p, sine, square
root): dec-B606's option A, no new dependency. Enclosures are intersected across refinement levels, so refining
never widens a value, which subsumes this plan's earlier "refinement never raises the floor". A non-finite value
is undecided, never decided on.

The cheaper certificate (S1). The certificate as first specified (cells with per-term majorants, dropped nodes one
by one, h found by walking the ladder) costs far more than the main sum. The proof ("A cheaper certificate")
proves the pieces of a cheaper one: an aggregated majorant on the edge lines with O(1) work per cell, in closed
form for f; a monotone tail in place of node-by-node truncation; and the step chosen from the bound instead of by
walking the ladder. Not built or measured. S1 builds it with its cell settings fixed before any evaluation,
proves the pieces the proof leaves open (how many cells, how large a J_0) before timing, checks it at the pins,
and then times it.

The cost of a certified step. Term evaluations per row step counted under the certificate as first specified
(the a-priori cascade, which stopped at its first setting everywhere), priced at kernel times in plain double
without the radius bookkeeping (116 ns an f node, 12.6 ns a q node term), so lower estimates; the gammas and saddle
solves are not included. Today's sum timed from the package's own sampler on the same machine. The density part
of a step in microseconds at the mean, in parentheses the main sum alone at the coarsest certified step, which no
certificate can go below:

| b | sum | psi 0 | \|psi\| 8 | \|psi\| 50 | \|psi\| 150 |
|---|---|---|---|---|---|
| 0.2 | none | 62 (35) | 53 (18) | 91 (37) | 323 (189) |
| 1 | 0.14 to 0.16 | 42 (16) | 43 (11) | 88 (23) | 323 (117) |
| 16 | 2.2 to 2.5 | 33 (5.7) | 35 (4.1) | 112 (13) | 358 (51) |
| 64 | 8.7 to 10 | 31 (3.6) | 38 (4.6) | 146 (22) | 535 (67) |
| 1e3 | 136 to 157 | 41 (3.0) | 57 (6.5) | 347 (27) | 1,660 (115) |
| 1e4 | 1,360 to 1,570 | 109 (3.7) | 100 (7.1) | 951 (72) | 5,120 (305) |
| 1e5 | 13,600 to 15,700 | 874 (9.4) | 299 (24) | 3,910 (276) | 24,000 (1,760) |

So the route as first specified crosses the sum near b 250 at psi 0 and near 400 at |psi| 8, loses to 1e3 and
wins from 1e4 at |psi| 50, and loses at every shape tested at |psi| 150, where the K + 1 terms of each q node
decide it. The main sum alone would beat the sum at every whole shape from 1e3 up at every |psi| tested, about even
at 1e3 and |psi| 150: the certificate, not the method, decides where the route pays. The sizing's 3.1 to 3.9 us
at 64 and up, 8.8 at 16 and 7.0 at fractional weights were measured with the uncertified density.

The decision (S2). The enclosure of log r is formed from all four density enclosures plus the rounding of their
sum. The uniform has at least 53 bits, held as the interval [U_-, U_+) its drawn bits determine, with more bits
drawn when the interval meets the enclosure; a 32-bit uniform alone puts a grid of 2.3e-10 under every decision,
above the floor. The bits are formed from the Mersenne-Twister words of the chain's generator, the algorithm every
chain is created with. The step accepts when l_- >= min(0, log U_+), keeps when l_+ <= log U_-, and otherwise refines,
re-evaluating f(omega) only when its cached enclosure is too wide. At the top level a close call keeps the row's
old value and is counted (dec-B606): the count is kept per chain and summed, exposed by a sampler accessor, warned
about by run() when nonzero, and carried on bart2, bart and rbart fits, whose creators warn as well. No gate can
provoke a close call at the real floor, so a test-only setting inflates every enclosure until close calls occur;
the count, report and warning tests and the count's mutant use it. Whether the shipped C header gains a reader of
the count, for packages that build against it, is an ABI event, as dec-B601's accessor was, and an open call to the
maintainer at S2 ([Calls left open](#calls-left-open)). A fresh exact draw in place of an undecided decision stays declined: it does not keep PG(b, psi)
(the proof's exact counterexample, a stationary distribution 0.0108 from the target), whether a decision resolves
depending on the current value, and none exists at fractional shapes.

The floor. At the proof's pins (b 0.2 to 1e5, |psi| to 150, x from 4 standard deviations below the mean to 12
above), in ball arithmetic with the lead in double-double, each log-density's enclosure has radius at most 1.4e-12
and a step is a close call with probability at most 6.2e-12 + 2^(1 - beta), beta the uniform's bits: about one fit
in 1,600 of 10,000 rows by 10,000 sweeps would see one. That is a maximum over the pins, not a bound over every
reachable point, and it is Arb's radius, which an engine following the rounding rules exceeds by a factor not
measured. So S1 measures the engine's own floor at the pins, next to Arb's; the documentation quotes the engine's
measured maximum and nothing before it, and, as recommended ([Calls left open](#calls-left-open)), says the floor
is unmeasured above |psi| 150, outside the pins, where the count still runs (a fractional row can reach there; a
whole-number one takes the sum); S1 does not pass if the floor gives a close-call probability above 2.3e-11
a step, a tenth of a 32-bit grid step (about one fit in 430 of that size). If the gate fails, the build tightens
the evaluator's rounding (compensated node sums, the slack in the rounding bounds) or comes back to the maintainer
before taking B', a second precision level for close calls (dec-B606). The standard is the package's: exact under
perfect arithmetic, as today's sum is; each of the sum's comparisons against a 32-bit uniform can flip with
probability up to about 2^-32, about one draw in 20 to 30 million at a shape of 64, uncounted ("Floating point: the
standard every sampler in the package meets").

The K rule (S2). A cost choice: with r_1 rounded down every K >= 1 is exact and uniformly ergodic. The weight f/q
has right-tail limit exp(b D(psi, K)), D > 0 in closed form; K is chosen per row as the smallest with
b D <= 0.01. That bounds the tail limit only, not sup f/q, which the rule does not control: certified lower bounds
on sup f/q exceed exp(b D) at small b (1.263 against 1.010 at b 0.3, psi 8, K 1; 35.6 against 28.7 at b 0.2,
|psi| 150). The smallest K, from D in 40-digit arithmetic:

| \|psi\| | b = 0.2 | 1 | 16 | 100 | 1e3 | 1e4 | 1e5 |
|---|---|---|---|---|---|---|---|
| 0 | 1 | 1 | 1 | 2 | 3 | 4 | 6 |
| 8 | 1 | 2 | 3 | 5 | 8 | 12 | 19 |
| 50 | 10 | 14 | 26 | 38 | 60 | 95 | 150 |
| 64 | 13 | 19 | 35 | 51 | 80 | 128 | 202 |
| 100 | 23 | 33 | 59 | 86 | 137 | 217 | 344 |
| 150 | 37 | 54 | 96 | 140 | 223 | 353 | 560 |

The first amendment's 110 to 120 at |psi| 150 came from a truncated sum that underestimated D. At 41 ns a gamma the
proposal's draws alone cost about 23 us a row at the extreme, and q at two points with K + 1 terms per node and
cell adds more (the cost table). The scratch build's cap at 48 broke the rule above b 81 at |psi| 64, 6.1 at 100
and 0.63 at 150. Such |psi| arise in fits near separation, not ordinary ones.

The cancel (S2). Polled between rows only, never inside a row's refinement: a refinement's length depends on
omega_i, so a cancel landing inside it would leave the row unchanged by a rule read from omega_i.

Rejection from the same proposal with a certified bound M on f/q would give independent exact draws, which would
also serve multinomial and fractional-shape rows a setter changed. The proof's explicit bound is finite but useless
(about 1e18 at b 1, psi 0, K 1, against a certified lower bound of 1.07), and no certified sharp bound exists.
Recorded as research, not built.

Where the route keeps the posterior (the proof's sweep, negative binomial, per-row rule, starting values and
multinomial sections, with exact finite-state examples):

- Within a sweep it holds: rows are independent given the trees, and a row's step reads only that row's omega and
  psi.
- Estimated-r nbinom holds under this design's rule: on a sweep where r moved, every row draws fresh and exactly;
  where r stayed, the step applies. psi_i = eta_i - log r, so a moved r changes both the shape and the tilt. r's
  draw never reads omega, so on the event that r stayed the old omega is a draw from its conditional at that r.
  Checked exactly: the rule exact, always applying the step 0.039 from the target in total variation. Its cost: r
  moved on 38 to 49 percent of sweeps at 200 rows and 0.3 to 0.6 percent at 2,000 in the sizing, and each move draws
  every row by the sum. The remedy is an r step that reads omega (a Metropolis step on r given omega), which never
  forces fresh draws; it is also what real dispersion needs, and is recorded for that plan, not built here.
- The per-row choice between the sum and the step keeps each row's conditional when it depends on anything but the
  row's own omega_i; a choice read from omega_i can break it (an exact four-state example: 0.258 in total
  variation, each kernel alone exact). The rule below reads b_i and psi_i, not omega_i.
- Any starting value is valid: a start needs no exact draw, only a chain that converges from it, with convergence
  from a fixed start assumed as for today's sampler.
- Multinomial does not hold: [`MultinomialForestCombiner::drawForestGlue`](../../src/bartcore/combiner.hpp) draws
  each category's omega as a temporary latent against the current margins just before that category's forest,
  and with persistent latents the forest steps are not invariant for the target the step needs (an exact example:
  0.0061 in total variation). Multinomial keeps the exact sum under every exact setting; an approximate setting,
  being a fresh draw, can serve it. A forced route on multinomial is refused naming the method, as recommended
  ([Calls left open](#calls-left-open)).

### The per-row rule

dec-B590 ruled the route at a row whose shape is above 64 or not a whole number, the sum otherwise. The 64 was
where the sizing's uncertified route was twice the sum's speed; under the certified cost the route loses to the
sum at 64 at every |psi| measured, while a cheap certificate could bring it below the sum there (the main sum alone
costs 3.6 us at b 64, psi 0, against the sum's 8.7 to 10). S1's timed rule refines the threshold (the proof,
"Changes the design needs", item 8, of the note dec-B606 rules on): at a whole-number row, the step if and only if
b_i > b*(|psi_i|), the sum otherwise; a fractional row always takes the step. b*(|psi|) is the timed value or 64,
whichever is larger, so dec-B590's guarantee, that a fit whose shapes are all whole numbers of 64 or less draws as
the build before it, holds whatever S1 measures; the proof's counts put the timed value at 250 and up. The
dependence on |psi| revises dec-B590's rule, so it comes back to the maintainer with S1's timings before S2 makes
it the default.

S1's timing gate. S1 times, on the x86 bench box, quiet, one core (arm64 reported, not used: the constants are
fixed in the code for every user, and dec-B590's 64 was set on that box), the route's per-row cost c_R(b, psi) and
the sum's c_S(b) at |psi| in {0, 8, 50, 150} and b in {0.2, 1, 2.5, 16, 64, 1e3, 1e4, 1e5} (0.2 and 2.5 for the
fractional-shape budget). c_R is the four-evaluation cost of a step (no valid cache entry): the density part with
its bookkeeping, the saddle solves and the K + 1 gammas, the maximum over z in {-3, 0, 3}; the three-evaluation
figure, with the update, is reported beside it and decides only whether the cache is built. So b* is decided from
S1's timings alone, with no share of steps or of rows the rule moves to estimate from fits, and errs toward the
sum, as dec-B590 did. The timed b*(|psi|) is the smallest whole b with c_R <= c_S, log-interpolated between grid
points and rounded up; at each |psi_i| the rule takes the larger b* of the bracketing grid points; b* is infinite
where the route never wins and at every |psi_i| above 150, so the sum is used there.
The gate passes if the timed b* <= 1e3 at |psi| <= 8, where the sizing's binomial and negative binomial cells lie (the proof's
counts give about 250 to 400). If it fails, the step has no case at whole-number shapes and S1 goes back to the
maintainer before S2. Fractional shapes have no sum to fall back on and keep the route whatever it costs; no
inexact fallback is introduced. If c_R there exceeds 70 us a row step, ten times the sizing's figure, S1 reports
the cells that do (in the proof's counts, |psi| >= 50) for the help of real weights.

psi_i moves every sweep, so a row near b* moves between the sum and the step, starting each stepped sweep after a
sum draw with no valid cache entry, which the four-evaluation c_R already prices; the rule's table is fixed
constants in the code, so which draws a seed produces does not depend on the machine.

### Setters and marked rows

The latents-first constraint: an accept-or-keep latent is exact for a larger sampler that integrates it out only
where it is a fresh draw from its current conditional, so every row whose conditional a setter or a restore moved
since the last sweep is marked and drawn at the start of the next sweep, its mark cleared only after that draw. The
rows marked: every row for setOffset, setForestBasis, setResponse, setData, setModel's shape and setState; the
moved rows for setPredictor; the changed rows for setWeights and setActiveRows; rows changed inside the per-sweep
hook included. The setters themselves draw nothing. setSigma, in the latents-first list, is left out: sigma is
fixed at 1 on these families, so it moves no row's conditional.

The marks do not exist on main, where every active row is drawn every sweep; S2 builds them. Per row, per chain:
set by setOffset, setForestBasis, setResponse, setData, setModel's shape and setState (every row), setPredictor
(the moved rows), setWeights and setActiveRows (the changed rows), inside the per-sweep hook included; copied with
the sampler; never state; each cleared only after the row's draw, so a cancel leaves it set.

So the separate exact path this plan had for setters is one branch at the sweep's start (S2): a marked row at a
whole-number shape is drawn by the sum, and a marked row at a fractional shape takes the step from its installed
value (dec-B607). The order changes no sizing number. A larger sampler that changes the offset every sweep marks
every row every sweep, so its fit draws by the sum throughout at whole shapes and the route's speed does not reach
it; an approximate setting does.

The proof ("Setters, marked rows and the two contracts") separates two kinds of larger sampler. One updates its
own unknowns from their conditional given the latents, reading getLatents: every setter is then exact with either
a fresh draw or the step at the rows it changed. The other integrates the latents out and does not read them: it
is exact with fresh draws at changed rows, and the step at a changed fractional row leaves a small bias (in the
proof's exact example, a stationary distribution 0.0154 from the target in total variation). Only a fresh draw
would correct that, and none exists at fractional shapes. An update that reads the latents while targeting the
conditional without them is exact in neither case, fresh draws or not.

dec-B607 takes every setter at fractional shapes and documents this (S3). The help of setOffset, setPredictor,
setForestBasis, setResponse, setWeights, setActiveRows, setModel, setData and setState, and of getLatents, says:
changing data mid-chain is exact when the larger sampler updates its own unknowns from their conditional given the
latents; when its update ignores them, a changed row with a fractional weight is left with a small bias, so such a
sampler keeps whole-number weights on the rows it changes; and an update that reads the latents while targeting
the conditional without them is exact in neither case. The bias's size it quotes is the gate's measured one (below).

Rows that enter the fit. The proof's result for the latent-reading sampler holds for a fixed set of rows: when the
active set is among the larger sampler's unknowns, the target defines omega only on the active rows, so a move that
brings row j in cannot condition on omega_j. An inactive row keeps its last omega, stale while the trees move (the
engine skips it). At whole-number weights the returning row is marked and drawn fresh, exact either way. At a
fractional weight, the proof's referee confirmed (reading the proof, the ruling and the engine; no runs):

- A row setData adds, and a weight moving from 0 to a value the larger sampler draws, start from a value the step
  forgets only over sweeps: exact under neither kind of larger sampler until then. The help says so, and such a
  sampler keeps whole-number weights on those rows.
- setActiveRows reactivating a row, and a weight toggled between 0 and a fixed w_j, are exact for a latent-reading
  sampler on a pseudo-prior extension of the target (Carlin and Chib 1995): each inactive row carries
  omega_j ~ PG(b_j, 0), b_j its shape when active, which does not depend on the trees, so the tree step is unchanged,
  the frozen omega is invariant for it, and the target is the extension's marginal. Such a sampler updates each
  row's activity from its conditional given every row's omega, the stale ones included, with the active factor over
  the inactive one, 2^(-b_j) exp(kappa_j psi_j - omega_j psi_j^2 / 2), read from getLatents (which returns the
  stale value at an inactive row) and getFitsWithoutOffset (fits are kept at inactive rows). No engine cost while
  inactive omega are never touched, as today; setState carries them. It does not cover a weight drawn from 0 to a
  new value (the ratio then needs a density evaluation) or a row setData adds.

The pseudo-prior route is our construction, so a second independent referee checks it before S3 (dec-B591); the
help carries it, with a two-line example of the membership step, only once that check confirms it. The proof, when
it moves into the tree in S1, gains the scope note on its result for latent-reading samplers: "for a fixed set of
rows; an activity change at a fractional shape is covered only on the pseudo-prior extension".

The example (S3), in the help page on running the sampler inside a larger one (man/dbarts-embedding.Rd), pointed
to from getLatents and the setters: a logistic random-intercept sampler with groups g and intercepts u_g ~ N(0, tau^2), tau fixed. Each
iteration runs one sweep, reads omega = getLatents() and the trees' fit f = getFitsWithoutOffset(), draws each u_g
from its Gaussian conditional given the latents and the fit,

    v_g = 1 / (1 / tau^2 + sum_{i in g} omega_i),   u_g ~ N(v_g sum_{i in g} (kappa_i - omega_i f_i), v_g),

kappa_i = w_i (y_i - 1/2) (w_i (y_i - n_i / 2) for binomial rows), and installs u through setOffset. It runs at
fractional weights. The tools exist on main.

The example's gate (S3). The check that fits is an exact posterior, not a joint-distribution check as
[geweke-mc.R](../../benchmarks/R/geweke-mc.R) runs it: that check draws y from its model given the parameters, and
0/1 responses under fractional weights have no such model (a fractional weight is a power of the likelihood). Nor
is the example's own posterior computable: its intercepts are continuous and cross leaves, so quadrature costs
(grid size)^G times the gate's per-leaf integrals, infeasible past one group. The exact gate for an outer step
through setOffset (benchmarks/R/outer-step-exact.R on main, from the latents-first plan: one tree on an enumerated
tree space, an offset gamma v with gamma on a grid, leaves and k by quadrature, long-run quantities by batch-means
z) already has the example's shape in its logistic-conditional arm, a latent-reading outer step through setOffset,
differing only in the conditional (a discrete grid). Its logistic-conditional and logistic-integrated arms gain
fractional-weight variants, w fixed and uniform on (0.2, 0.8): about 20 changed lines (w times the log-likelihood
in the exact posterior, kappa = w (y - 1/2) in the conditional). The fractional logistic-conditional arm must pass.
Broken arms: a mutant build whose marked fractional rows keep the old omega with probability 1/2 and otherwise take
the step must fail under the outer step that integrates the latents out and must pass under the one that reads
them; and the conditional arm reading latents saved before the sweep must fail. The real step under the
integrating outer step cannot fail at a feasible size (a changed row stays stale only with the step's rejection
probability, 0.002 to 0.006 at the sizing's real-weight cells, against a fully stale latent's |z| of 17.8 in quick
mode), so it is run and its largest |z| reported, not gated, and the help's "small bias" quotes it. The example
itself is run by R CMD check and by a tinytest comparing its intercept means at whole-number weights, with the
method forced to the route, against a run with the method forced to the sum. Compute: 32 rows at 40 to 90 us a row
step is 1.3 to 2.9 ms a sweep, 20 to 50 minutes an arm in quick mode against about 35 seconds today; five arm runs,
about 2.5 to 4 core-hours.

### Real shapes

The route works at any real b > 0 given the density, which S1 supplies at every x > 0. That is the exact draw
dec-B568 waits on, within sweeps and from a starting value. A fresh exact draw at a fractional shape is another
matter: none is published below 1, and WPS's alternate sampler above 1 rests on Conjecture 8. So, under every
exact setting:

- Real weights on logistic, binomial and hazard.logistic are taken at creation, and every between-sweep setter and
  setState is taken (dec-B607), the rows they mark at fractional shapes served by the step from their installed
  values: exact for a larger sampler that reads the latents ([Setters and marked rows](#setters-and-marked-rows)).
  A starting point is served by the step as a cold start is. A forced sum at fractional weights has no draw and is
  refused naming the method, as recommended ([Calls left open](#calls-left-open)). dbartsDrawLatents, a fresh
  draw, keeps refusing fractional weights. A fresh draw at fractional shapes, by rejection from
  the proposal once its bound is certified or by a transport of the distribution functions, would make the
  integrated-out case exact too; research, not built.
- At fractional shapes the step is the only exact draw, and costly: under the certificate as first specified, 40 to
  90 us a row step at |psi| <= 50 and about 320 at 150, against the sizing's 7. S1 measures it; its cells above
  70 us go into the help of real weights.
- Integer plus fraction (PG(n + f) is PG(n) plus an independent PG(f)) would send only shapes below 1 to the
  route. Not taken: it is exact only with the fraction's latent kept as its own chain state (saved states,
  restores, copies), and the density must be certified below 1 regardless, which also covers 1 to 12.
- Real negative binomial dispersion (TODO negbin-real-dispersion) needs a fresh draw whenever r moves, so it needs
  either that rejection draw or an r step that reads omega; its own plan decides. Not in this build.

## The sizing, as run

Run 2026-11-03 on a scratch build carrying every survivor behind an environment switch (970 changed lines,
never lands; its density uncertified, with the defects found). Its route costs are those of that density; the
certified costs are above, and S1 retimes them. Timing: an x86-64 Linux box, one core, seconds per sweep, minimum
over 5 alternating rounds; median round-to-round spread 0.4 percent, 95th percentile 2.2; load at or below 1.32.
Effective draws: the laptop, 2 chains of 2,500 kept sweeps, rank-normalized bulk ESS of f(x_i) and of the
probability and tail ESS of a 5 percent quantile (Vehtari et al. 2021), on 8 cells; ten seeds of one arm in three
cells for the spread. Compute 1.6 core-hours.

| family | rows | trees | shape per row | arms |
|---|---|---|---|---|
| logistic with binomial rows | 200, 2,000 | 75, 200 | 16, 64, 256, 1e3, 1e4, 1e5 (the sum to 1e4 at 2,000 rows), and SIDS-like spread (lognormal, median 2,000) | sum, route, series 1, 5, 10, gamma, normal |
| logistic with real weights | 2,000 | 75 | uniform on (0.2, 0.8) and on (1.2, 3) | route, series 1, 5, 10, 20, gamma, normal |
| nbinom, fixed r = 1e3 | 200, 2,000 | 75 | counts mean 10 and 1e3 | sum, route, series 1, 5, 10, gamma, normal |
| nbinom, estimated r | 200, 2,000 | 75 | counts mean 20 and 1e3 | sum, route (fresh on r moves), series 1, 5, 10, gamma, normal |
| multinomial, 3 categories | 200 | 75 | 1e3 trials | sum, series 1, 5, 10, gamma, normal |

Costs measured:

- The sum: about 0.13 us per trial per row.
- The scratch route: 3.1 to 3.9 us a row at shapes of 64 and up (about two density passes a step at 17 to 21
  nodes, the Newton solve and K + 1 gammas); 8.8 at 16; about 7.0 at fractional weights; 6.2 to 7.7 at estimated-r
  nbinom with counts near 20; 48 at 200 rows with counts near 1,000, where r moved on 38 percent of sweeps.
  Acceptance 0.9998 to 1.0000 on every logistic and nbinom cell, 0.9936 and 0.9977 at the two real-weight cells.
- Everything else a sweep does: 0.13 to 0.40 us a row at 75 trees, 0.5 to 1.0 at 200 trees and in multinomial.
- Scratch route against sum: the sum 3.2 to 4.0 times faster at 16 trials; the route 1.9 to 2.1 times faster at
  64, 31 to 38 at 1,000, 300 to 350 at 10,000.

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

Under the certified route the route's own row changes (it wins in fewer cells, and only S1's timing says which),
and the approximations' ratios against it grow; no survivor's verdict turns on it, since every approximate
survivor already passed against the uncertified route, and the route stays the only exact draw at real weights.

Limits of the verdict:

- Trees: the grid has 75 and 200. With fewer trees the latent is a larger share of the sweep. Scaling the other
  work by the tree count (arithmetic, not measured), at 2,000 rows the gamma would pass against series 1 below
  about 30 to 50 trees and the normal against the gamma below about 17 to 48 (drop-threshold.R). Dropping them
  holds for 75 trees and up; a timing check at 2,000 rows with 25 and 50 trees precedes S4 (step 10).
- Series 5 and 1 clear the bar by 6 to 18 percent, robust to timing noise, but only at 75 trees (they fail at
  200) and at the cost of R's gamma generator: halving its cost would bring them to about 1.26 (arithmetic).
- ESS at 2,000 rows with large shapes: 5 to 23 effective draws per 5,000 sweeps, BART itself mixing slowly
  there, so those cells rest on the timing and the equality above.
- The real-weights cell with weights 1.2 to 3 on 0/1 responses drifted outward over thousands of sweeps; the
  exact sum at weight 2 and the sum on every row entered twice drift alike, so it is the model (weight 2 is the
  duplicated-rows posterior), not the sampler (drift-rev.R). The first round, before the drift, was used.
- No joint-distribution check at real weights (0/1 responses have no generative model under fractional weights);
  moments, KS against the integrated density and the posterior-mean comparison (0.57027 against 0.57021 to 0.57036
  for the approximations) stand in.
- Not run: arm64 timing, a faster gamma generator, real-r nbinom.
- The scratch route is biased in a band no sizing cell entered (b about 2 to 12 with |psi| of 80 or more: the
  mean 0.6 to 0.9 percent high at b 8 and 11.5, z 11.5 to 22; route-onset.R, route-b2.R, re-run 2026-11-05),
  from its uncertified q. Across 63 cells outside the band (b 0.2 to 1,000, psi 0 to 150) every |z| was 2.62 or
  less and the smallest KS p 0.015 (route-check.R). Its density is replaced by S1's.

## What is built

Tiers by [Process by risk](README.md#process-by-risk). Everything ships in 1.0-0 (dec-B592); the setting's shape is
dec-B593's.

| slice | content | tier | lines | agent-hours |
|---|---|---|---|---|
| S1 | the certified PG density, f and q, at every x > 0 for every reachable (b, psi): the parabolic contour through the saddle with its a priori bound; the tilted, cancellation-free form; the lead and the cached value in double-double with five double-double functions with bounds; rounding carried rigorously; q from the sampler's own rates; the cheaper certificate, its settings fixed before any evaluation and its open pieces proven before timing and checked by an independent referee before the timing gate is read; pins against reference tables generated in Arb and checked in (tests/cpp cannot call Arb); the double-double functions' gate; the floor gate; the timing gate on the x86 box and b*; the proof's source and scripts moved into docs/design/approximate-draws-proof/ with a README, and the design note | engine, neutral | 1,950 | 27 |
| S2 | the per-row marks; the route at the sweep's start for logistic, binomial, hazard.logistic, fixed- and estimated-r nbinom (psi_i = eta_i - log r); marked rows at whole shapes by the sum; the per-row rule by b*(\|psi_i\|) once ruled; the moved-r rule; the cache keyed on the point, if S1's timings call for it; the K rule from D's closed form with r_1 rounded down; the decision on enclosures with an interval uniform of at least 53 bits from the Mersenne-Twister words; close calls kept, counted per chain, read by a sampler accessor, warned by run() and carried on bart2, bart and rbart fits, with a test-only setting that provokes them; the cancel between rows; the setting's engine side and an internal R switch forcing the method, ahead of S5's door; the joint-distribution harnesses for the Polya-Gamma families; an equivalence scenario of binomial rows above b* | engine, shifting where selected, with a small R piece for the count and the switch | 1,750 | 24 |
| S3 | real weights under the route on logistic, binomial, hazard.logistic; every setter and setState taken at fractional shapes, marked rows served by the step; the help passage of the setters and getLatents, with the clause on rows that enter the fit, and the pseudo-prior route with its two-line example once the second referee confirms it; the example, its tinytest and the outer-step gate's fractional arms; an equivalence scenario at fractional weights; the frequency-weights plan's interaction paragraph brought up to date | engine plus R surface | 450 | 10 |
| S4 | the series as an approximate method for all PG families including multinomial, its terms chosen per row from the tolerance by a table computed once by inversion (0 terms, the gamma, at the loose end); the normal only if step 10's check passes it; the design note's section for it | engine, posterior-changing where selected | 300 | 6 |
| S5 | the setting at the doors users reach: the method by name, NULL the automatic rule, and a separate tolerance read only by the approximate method (dec-B593), the names put to the maintainer with the four calls below; a forced route on multinomial and a forced sum at fractional weights refused by name; dec-B202's cap message corrected; help quoting the engine's measured floor, unmeasured above \|psi\| 150, and the close-call count; NEWS; dec-B591's teensy vignette if S1 leaves the route a case worth motivating | R surface | 500 | 8 |

The sizing, re-audited after the critique: about 4,950 written lines and 75 agent-hours, against the critic's
4,450 and 69 and the previous 3,700 and 56. Per slice:

- S1, 1,950 lines and 27 hours (critic 1,800 and 26; previously 1,500 and 20). Added: the double-double functions'
  gate against 200-digit references (about 100 lines, 1 hour); the reference tables generated in Arb and checked in
  (about 200 lines); the independent referee of S1's own proofs (3 hours); the design note (about 150 lines) and
  the proof's move with its README (2.5 hours with the note; the TeX and scripts, about 2,000 lines, copied rather
  than written, are not counted). The timing gate is simpler, its cost decided from S1's measurements alone. The
  critic's figure omits the design note and the move, which accounts for the difference.
- S2, 1,750 lines and 24 hours (critic 1,450 and 22; previously 1,000 and 15). Added: the marks (about 250 lines, 4
  hours); three joint-distribution harnesses extended from the gaussian one (about 350 lines, 3 hours, counted
  higher than the critic's); the test-only setting, the accessor and the fit fields beyond the earlier count (about
  100 lines, 1.5 hours); the internal switch and the new equivalence scenario under its oracle (about 60 lines, 1
  hour). Without the cache, about 100 lines and 1.5 hours less.
- S3, 450 lines and 10 hours (previously 450 and 8; the critic gave no figure). The quadrature arms go (about 130
  lines less); the 20-line fractional arms, the tinytest, the mutant arms, the fractional equivalence scenario and
  the frequency-weights edit come; the second referee's check of the pseudo-prior route (1 hour), its help passage
  and example, and the gate's longer runs add 2 hours net.
- S4, 300 lines and 6 hours, unchanged; a tolerance separate from the name costs as the single overloaded setting
  did.
- S5, 500 lines and 8 hours (previously 450 and 7): the two refusals by name, the cap message and the floor's
  caveat.

The defaults study is not sized here.

Compute per gate, one core unless named; about 2.3 core-hours in S1, 5 in S2, 5 in S3 and 0.7 before and in S4,
about 13 locally, and about 20 in CI:

| slice | gate | compute |
|---|---|---|
| S1 | reference tables for the pins and the double-double functions (Arb, offline; q at K up to 560 is the slow part) | about 2 core-hours |
| S1 | pins, floor gate, double-double gate in tests/cpp | minutes |
| S1 | timing gate on the x86 box: 96 points of (b, psi, z), both methods | about 0.3 core-hours |
| S2 | joint-distribution checks, three families at quick size (300 chains of 100 refreshes, 3 sweeps between, about 100 rows at 40 us a row step under the forced route): about 10 minutes a family; the five broken builds again, each a build | about 1.5 core-hours |
| S2 | KS cells, 84 cells of 1e5 steps each, dominated by the 48 cells at \|psi\| 80 and up at 0.3 to 1.7 ms a step | about 1 core-hour |
| S2 | the single-tree exact posterior with rows of b trials under the forced route | about 0.5 core-hours |
| S2 | simulation-based calibration with binomial rows under the forced route: at its default 150 rows, about 6 ms a sweep, some 20 core-hours at the logistic arm's size; CI's scheduled run, not a local gate, its size set in S2 against that | about 20 core-hours, in CI |
| S2, S3 | the engine tier's battery and the exact gates, quick | about 2 core-hours each |
| S3 | outer-step gate, five arm runs at fractional weights in quick mode, 20 to 50 minutes each | about 2.5 to 4 core-hours |
| S4 | the table of terms by inversion, and TV per row in tests/cpp | about 0.5 core-hours |
| before S4 | step 10's timing at 25 and 50 trees | about 0.2 core-hours |

Gates for S1:
- pins in tests/cpp: at every pin the enclosure contains an independent reference, over b in [0.2, 1e5],
  |psi| <= 150 and x from 4 standard deviations below the mean to 12 above, the far tails the proposal visits, with
  the proposal at the rule's K from the sampler's own rates; for q one reference independent of the parabolic
  contour (the vertical line at large total shape, Moschopoulos's series where the largest rate times x is small);
  the cheaper certificate checked at the same pins; the references are tables generated in Arb and checked in;
- the five double-double functions against checked-in references at 200 digits over their reachable ranges, the
  largest error reported against the stated bound (the proof assumes at most 2^-100, unproven for any library); a
  breach fails S1;
- S1's own proofs (the cheaper certificate's cells and J_0) checked by an independent referee before the timing gate
  is read (dec-B591);
- the floor gate: the engine's radius at every pin reported next to Arb's; the step floor it implies at most
  2.3e-11; on failure the rounding tightened, or back to the maintainer before B';
- the timing gate of [The per-row rule](#the-per-row-rule), on the x86 box: the timed b* <= 1e3 at |psi| <= 8, or
  back to the maintainer before S2; the fractional-shape cells above 70 us reported; the three-evaluation figures
  against the materiality threshold S1 stated, deciding the cache;
- the engine tier's battery, sanitizers on tests/cpp included.

Gates for S2 and S3 (exactness), each with the exact sum as control, passing first. Every gate forces the method
to the route through S2's internal switch: at b 16, say, the rule picks the sum.
- joint-distribution check (Geweke 2004), extended from [geweke-mc.R](../../benchmarks/R/geweke-mc.R) (gaussian
  only) to logistic with binomial rows and fixed- and estimated-r nbinom: the successive-conditional arm runs m >= 3
  sweeps between refreshes of y, with the method forced to the route, so sweeps 2 to m take the step from a
  sum-drawn or stepped omega (a refresh through setResponse marks every row, so with one sweep between refreshes the
  check would never run the step); broken builds, each shown to fail: the step applied after r moved, the point
  match removed (b alone; if the cache is built), q tilted from psi = 0, a marked whole-number row taking the step
  in place of a fresh draw, the per-row choice read from omega_i;
- two-sample KS of the route's omega at fixed (b, psi) against the exact sum (an independent reference at
  fractional b), and moments to the third cumulant, on b in {0.2, 0.5, 1, 2, 2.5, 4, 8, 11.5, 12.5, 16, 64, 1e3}
  by |psi| in {0, 8, 50, 80, 100, 120, 150}, with the close-call count zero in every cell; a broken build with q on
  the scratch build's uncertified contour must fail at b 8, |psi| 100 and b 2.5, |psi| 150;
- the exact single-tree posterior of [logistic-reference.R](../../benchmarks/R/logistic-reference.R) part 1
  with rows of b trials;
- simulation-based calibration through [sbc.R](../../benchmarks/R/sbc.R)'s logistic machinery with binomial
  rows, its exact arm passing first; run in CI on its schedule, sized in S2 against its compute (below);
- S3: the outer-step gate's fractional arms, the mutant arms and the example's tinytest
  ([Setters and marked rows](#setters-and-marked-rows)), the real step's largest |z| under the integrating outer
  step reported;
- close calls: under the test-only setting that inflates every enclosure, the count, its report on the sampler and
  the fits, and run()'s warning, each tested; the warning counted with withCallingHandlers;
- RNG and baselines: S2 is bitwise on every current baseline (checked: no scenario's shapes reach 64; logistic
  weights are 1 to 3, the nbinom counts small; the four snapshot files are probit, gaussian and xbart); S2 adds an
  equivalence scenario of binomial rows above b* drawn by the route and S3 one at fractional weights, each recorded
  with its oracle under the baselines MANIFEST's rule P17;
- the engine tier's battery: tests/cpp, full tinytest, sanitizers on tests/cpp and the R-loaded path, R CMD
  check, the arm64 snapshots, bench-sampler.R compare on a quiet machine (maintainer-run) for the hot path; and,
  since the fit carries the close-call count, the exact gates the workflow lists, quick;
- mutants: the review runs the mutation battery's anchors plus mutants of the accept comparison, the K rule (a
  cap at 48), a non-finite ratio decided on, the cache's point match (if built), the moved-r branch, the
  marked-row branch, a mark cleared before the row's draw, the per-row rule and the close-call count (with power
  under the test-only setting).

S4's gates: its TV per row against inversion in tests/cpp, the chosen K meeting the tolerance at every grid cell of
the screen's table and between; the exact gates stay on exact settings.

## The defaults study, after the build (outline)

It runs before 1.0-0, right after the build (dec-B592), and 1.0-0 ships the default it derives; dec-B590's rule
holds only for the development builds between. It measures, on fits shaped like users' and on real covariates
with real trials, each built setting's effective draws per second and, for the approximate ones, their effect on
what users read: coverage and width of 90 and 95 percent intervals for probabilities, rates and contrasts, RMSE and
bias, against exact posteriors (the one-parameter toy by quadrature, the single-tree enumeration) and
semi-synthetic truths, at equal running time. Its fractional-weight cells use fractional counts with proportions,
or weights below 1, never 0/1 responses at weights above 1: weight 2 on a 0/1 row is the posterior of that row
entered twice, which drifts toward separation under any sampler. It covers small shapes and |psi| to 50 at least,
where the approximations are least accurate, and tree counts below 75. Choosing among exact settings changes no
posterior and is a speed finding. Making an approximate setting a default is a change users' results carry, and
comes as a short research report under the studies standard (dec-B559): at most 5 pages, readable alone, two
referees first.

## Decision

The three open calls this plan put are ruled. Their records stay here, in short; the ledger entries carry the
options and the maintainer's words.

### Open call 1: the default until the defaults study

Ruled, dec-B590 (2026-11-05): per row, the route where the shape is above 64 or not a whole number, the sum
otherwise; the route takes this place only once proven (dec-B591) and its gates pass. The options as put: (a) that
rule, recommended; (b) the sum wherever it can draw, the route by setting and for fractional shapes; (c) the
approximation above a shape. S1's timed rule may raise the threshold by |psi|, never below 64, and that dependence
comes back to the maintainer ([The per-row rule](#the-per-row-rule)); the per-row form, the guarantee for shapes of
64 or less, the fractional rows and the forcing setting stand.

### Open call 2: what ships in 1.0-0

Ruled, dec-B592 (2026-11-05): everything, the approximate settings opt-in and documented as approximate, revising
dec-B14's 1.0 clause; the defaults study runs before 1.0-0, right after the build, and 1.0-0 ships its default.
The options as put: (a) the exact pieces in 1.0-0, the approximate settings after; (b) everything, recommended;
(c) nothing from this arc.

### Open call 3: the setting's shape

Ruled, dec-B593 (2026-11-05): the method by name (the exact sum, the exact accept-or-keep route, the approximate
series), NULL the automatic rule, and a separate tolerance, the target total variation per draw, read only by an
approximate method and ignored under the exact ones, as base R ignores approx's f under "linear". The names are
settled at the build. The options as put: (a) names only at one fixed accuracy; (b) names plus the number of terms;
(c) one setting taking a name or a number, first recommended; (d) a name and a tolerance, the revised
recommendation. S5 proposes the names and the tolerance's default (the plan's figure: 1e-6 total variation per
draw, between the 5- and 10-term series at shape 16).

### The proof's two questions

Ruled, dec-B606 (2026-11-06): what exact means in floating point, as the proof's option A
([The exact accept-or-keep route](#the-exact-accept-or-keep-route), the decision and the floor). Ruled, dec-B607
(2026-11-06): what the setters do at fractional shapes, every setter taken with the help passage and the example
([Setters and marked rows](#setters-and-marked-rows)); the planned refusals are dropped.

### Calls left open

To the maintainer, each when its evidence exists:

- With S1's timings, before S2 makes it the default: the threshold's dependence on |psi|, b*(|psi|) never below
  dec-B590's 64. Recommendation then: adopt it as timed.
- From S1, if its timing gate fails (the timed b* above 1e3 at |psi| <= 8): the step has no case at whole-number
  shapes; before S2.
- From S1, if its floor gate fails after the rounding is tightened: before taking B' (dec-B606).
- At S2: whether the shipped C header, dbarts.h, gains a reader of the close-call count, for packages that build
  against it. An ABI event, as dec-B601's accessor was. Recommendation: yes, named as the header's other readers
  are, since a consumer running chains through the C API otherwise sees no close call at all.
- At S5, with the setting's names and the tolerance's default (dec-B593), four smaller calls, each with its
  recommendation:
  - a forced route on multinomial, which has no exact route: refuse, naming the method (recommended), against
    quietly drawing by the sum;
  - a forced sum at fractional weights, which has no draw: refuse, naming the method (recommended), against
    quietly taking the route;
  - dec-B202's cap message ("each unit of a count is a Polya-Gamma draw per row per sweep"), false once the route
    serves large counts: correct the message, the cap itself unchanged (recommended);
  - the help saying the floor is unmeasured above |psi| 150, outside the pins, where the count still runs: say so
    (recommended).

Left to the build, with the evidence named:

- S1 states, before timing, the change in b* that counts as material; the cache is built only if the three- and
  four-evaluation timings move b* by that much.
- dec-B591's vignette is the agents' judgment by its own words; decided once S1's timings show where the route
  pays, and written in S5 if it does.

Checks by independent referees (dec-B591), before the step they gate: S1's own proofs, before its timing gate is
read; the pseudo-prior route for rows that re-enter at fractional weights, before S3, the help carrying it only
once confirmed.

The 53-bit uniform formed from the Mersenne-Twister words is a build detail, not a call.

## Constraints

- No package code from the sizing lands; slices follow their tier.
- Until the defaults study, a fit whose shapes are all whole numbers of 64 or less draws as the build before it:
  b* is never below 64.
- Draws come only from the chain's RNG, the uniform's extra bits included; the cache, if built, is scratch, never
  state, and is read only on a bit-for-bit match of the point.
- No new system dependency (dec-B606's option A).
- The cancel is polled between rows, never inside a row's step.
- A close call is never silent: counted per chain, read by an accessor, warned by run() and on the fits.
- S2 builds the marks and meets the latents-first marked-row constraint; the marks are never state and a cancel
  leaves them set.
- Order: S3 lands after frequency-weights; S1 and S2 on either side. Whichever of S2 and frequency-weights lands
  second rebases over the other's LogisticResponse and NBResponse hunks.
- Nothing from BayesLogit's source (GPL 3; dbarts is GPL 2 or later).
- Out of scope: the probit binomial, real negative binomial dispersion and the r step that reads omega (recorded
  for that plan), the dec-B202 cap, B' and a precision ladder.

## Steps

1. Done: the scratch sizing build (970 lines), reviewed before its run.
2. Done: the sizing and its verdict per survivor under the stated rule; the independent review (sound with
   corrections) and the first amendment.
3. Done: the open calls ruled (dec-B590, dec-B592, dec-B593).
4. Done: the latents-first engine slices, landed on main 2026-11-05.
5. Done: the proof of dec-B591, checked by two independent referees; its two questions ruled (dec-B606,
   dec-B607); the amendment for them; the blind critique of that amendment and the second amendment; the proof's
   referee on rows that re-enter the fit.
6. Rebase the branch onto main.
7. S1 through its tier: the proof and the design note into docs/design/; its own proofs checked by an independent
   referee; its double-double, floor and timing gates; then to the maintainer with the timings for the threshold's
   dependence on |psi|, and for S1 itself if either of its two decision gates fails.
8. S2 through its tier, with the C-header call put at its start.
9. The second referee's check of the pseudo-prior route; frequency-weights landed; then S3 through its tier.
10. Before S4: time 2,000 rows at 25 and 50 trees on the scratch build (about 0.2 core-hours); the normal is
    built in S4 if it passes the rule there (about 30 lines).
11. S4, then S5, with S5's calls put to the maintainer.
12. The defaults study, planned once the build lands and run before 1.0-0.

## Verification

- The screen's numbers reproduce from scratch/apx-design/probes (laterr.R, screen.R, tiltq.R, supcheck.R),
  scratch/apx-check/probes (tail-ratio.R, pcg-exact.R) and scratch/apx-rev/probes (tv-small.R, series-check.R).
- The sizing's verdict follows from scratch/apx-size/results (tables.txt, analysis.rds) by the rule above,
  mechanically (scratch/apx-rev/probes rule-recompute.R).
- The scratch route's band and the first build's defect reproduce on the scratch library from
  scratch/apx-rev/probes (route-onset.R, route-b2.R, route-b16.R, q-b4.R, q-b8.R, old-refine.R).
- The proof's figures (the cost table, the K table, the floor, the exact examples) reproduce from the scripts its
  appendix "Numerical verification" lists, in docs/design/approximate-draws-proof/ once S1 moves them there; its
  README says how to build the PDF.
- S1's pins hold every reference inside the enclosure; its floor and timing gates report their numbers; S2's and
  S3's broken builds fail their gates; a gate passing one has no power and stops the slice.

## Taken and not taken

- From the critique (September design's revision): the exact-moment remainder and the re-ranking; the
  unbounded weight of the gamma proposal; certified meaning; exact draws at swaps.
- From the independent check: estimated-r nbinom holds under the moved-then-fresh rule; multinomial excluded;
  exact draws for rows a setter changed; f and q both at the current psi; the weight bounded but not uniformly, so
  K grows with b and |psi|; the density's error as a bias floor and the certified decision; real shapes (none exact
  below 1, the alternate's conjecture).
- Added by this design: rejection from the proposal as an independent exact draw, research.
- From the sizing: the verdict per survivor; the costs (the route's replaced by the proof's); the
  integer-plus-fraction split, considered and not taken ([Real shapes](#real-shapes)).
- From the proof: every item of its "Changes the design needs" ([The proof's changes](#the-proofs-changes)),
  with two replaced by rulings: its recommendation to refuse setData and setState at fractional shapes, and the
  design's refusals at the other setters, both replaced by dec-B607's taking every setter; and its fallback to B'
  on a failed floor gate, which dec-B606 makes a return to the maintainer.
- From the blind critique of the amended plan: all 18 amendments ([The second critique](#the-second-critique)).
- From the proof's referee on rows that re-enter the fit: the critic's claim confirmed for rows setData adds and
  weights drawn from 0, refuted for reactivation and fixed-weight toggles, which the pseudo-prior route covers.

## Review amendments

The sizing's independent review (2026-11-05, sound with corrections) proposed 15 amendments, applied in the first
amendment: q certified as well as f; a true error bound; never deciding on a non-finite or unresolved ratio, the
fallback to a fresh draw declined; the K rule past |psi| 64; gate cells at b 1 to 12 and |psi| to 150, with a
broken build on the uncertified q; small-b and |psi| 50 rows in the screen; Open call 1 (c)'s ratios and the other
work per row; the crossover between 16 and 64; estimated-r nbinom's cap and the r step that reads omega; gamma and
normal dropped for 75 trees and up only; series 5 and 1 only at 75 trees and R's gamma cost; equal ESS per sweep
resting on theory and two high-ESS cells; the integer-plus-fraction split, not taken; no 0/1 responses at weights
above 1 in the defaults study; the first build's defect as the refinement doubling the nodes. The proof revises
three of them: the K figures past |psi| 64, the density by region (now one rule), and the refusals at fractional
setters (dropped by dec-B607).

## The second critique

A blind critique of the plan as amended for the proof (2026-11-06) proposed 18 amendments, all taken, some as the
orchestrator amended them:

1. The joint-distribution check extended to the Polya-Gamma families, m >= 3 sweeps between refreshes, the route
   forced, so the step runs: S2's gates.
2. The example gate's broken arm a mutant keeping the old omega at marked fractional rows half the time; the real
   step under the integrating outer step reported, not gated: [Setters and marked rows](#setters-and-marked-rows).
3. The gate simplified to fractional-weight variants of the outer-step gate's logistic arms, the example run by
   R CMD check and a tinytest: the same.
4. The per-row marks built by S2, setSigma left out: the same, and S2.
5. The proof's source and scripts into docs/design/approximate-draws-proof/ with a README, no PDF, and a design
   note beside it, with S1: Context.
6. S1's own proofs checked by an independent referee before the timing gate is read: S1's gates.
7. c_R the four-evaluation cost, the maximum over z; b* from S1's timings alone: [The per-row rule](#the-per-row-rule).
8. b* the timed value or 64, whichever is larger: the same, Constraints.
9. Rows that enter the fit at fractional weights: the claim put to the proof's referee, who confirmed it for added
   rows and drawn weights and found the pseudo-prior route for reactivation and fixed toggles, itself to be checked
   by a second referee before S3: [Setters and marked rows](#setters-and-marked-rows).
10. The close-call count testable through a test-only setting, read by an accessor, warned by run(), on the fits;
    the C header's reader an open call at S2: the route (the decision), [Calls left open](#calls-left-open).
11. Order with frequency-weights; that plan's paragraph brought up to date in S3; dbartsDrawLatents keeps refusing:
    Context, Constraints.
12. Baselines checked bitwise, new equivalence scenarios under P17, an internal switch forcing the method: S2's and
    S3's gates.
13. The sizing re-audited with compute per gate: [What is built](#what-is-built).
14. The double-double functions' own gate: S1's gates.
15. Timing on the x86 bench box: [The per-row rule](#the-per-row-rule).
16. The cache built only if S1's timings move b* materially, the threshold stated before timing: the route (the
    cache).
17. The broken build renamed: the point match removed (b alone): S2's gates.
18. Four user-visible calls listed with recommendations, for S5; the uniform from the Mersenne-Twister words, a
    build detail: [Calls left open](#calls-left-open).

## The proof's changes

Each item of the proof's "Changes the design needs", and where it is built:

1. One density rule, the parabolic contour for f and q at every shape; the split at 12 and the Talbot contour go:
   S1 ([The exact accept-or-keep route](#the-exact-accept-or-keep-route), the density).
2. The ratio in tilted, cancellation-free form; q at both points; f fresh at the proposal, the cached value moved
   by the exact update; three evaluations on a valid entry, four otherwise (always four without the cache): S1
   forms, S2 uses (the ratio).
3. The floor figure: 1.4e-12 a density and 6.2e-12 + 2^(1 - beta) a step at the pins, the engine's own measured
   before it is quoted, exact in the package's standard (dec-B606): S1 and S5 (the floor).
4. The proposal's own numbers, the first rate rounded down: S1 (q) and S2 (the proposal).
5. The decision counts all four densities and the ratio's arithmetic; enclosures intersected; non-finite values
   undecided: S2 (the decision).
6. The uniform, at least 53 bits held as an interval: S2 (the decision).
7. The cache keyed bit for bit on omega_i and b, moved in double-double: S2 (the cache), built only if S1's timings
   call for it.
8. The cost figures, the cheaper certificate built and timed in S1, the timing gate and its pass rule, b* infinite
   above |psi| 150: S1 (the cost, [The per-row rule](#the-per-row-rule)), with c_R the four-evaluation cost at the
   worst z and b* never below 64; the dependence on |psi| to the maintainer.
9. The K rule as a cost choice; the table at |psi| 150 of 140 to 560; about 23 us of gammas at the extreme: S2 (the
   K rule).
10. A fresh draw on an undecided decision stays declined: S2 (the decision).
11. Negative binomial with psi_i = eta_i - log r; the rule stands; real dispersion needs an r step that reads omega:
    S2 (where the route keeps the posterior); out of scope.
12. The setters and the two kinds of larger sampler, the help saying the guarantee covers those two only: S3,
    under dec-B607 ([Setters and marked rows](#setters-and-marked-rows)), with the scope note on a fixed set of
    rows and the pseudo-prior route for rows that re-enter.
13. The cancel polled between rows: S2 (the cancel).
14. Starts: any starting value valid: S2 and S3 (where the route keeps the posterior).
15. Multinomial keeps the sum: S2 and S4.
16. Pins with the sampler's own rates over the far tails, q against a reference independent of the contour: S1
    (gates).
17. Two S1 gates, the floor and the timing, with their failure paths: S1 (gates), the floor's under dec-B606.

## Not verified

- The route's cost in an engine: the proof counted terms and priced them at kernel times without the radius
  bookkeeping, so its figures are lower estimates; the cheaper certificate is proven but neither built nor
  measured.
- The engine's close-call floor: the proof's is Arb's radius at the pins, and an engine following its rounding
  rules exceeds it by a factor not measured; the double-double functions' error bounds and the math library's
  accuracy are assumed, not proven for a particular library.
- That 53 bits make every enclosure finite and narrow at every reachable point: checked at the pins, not proven.
- A sharp certified bound on f/q, so the rate of convergence has no usable number.
- The crossover between the sum and the route, and every timing on arm64.
- K past |psi| 64 is from the closed form of D, not from chains at those K.
- Gamma and normal below 75 trees: arithmetic only, until step 10.
- That effective draws per sweep match across arms beyond the two high-ESS cells and theory.
- The cost of the table of terms behind the tolerance, estimated, not built.
- How large the real step's bias under the integrating outer step is in the gate's design: S3 measures and
  reports it, and the help quotes it.
- The pseudo-prior route for rows that re-enter at fractional weights: one referee's reading, a second referee's
  check pending before S3.
- The compute figures per gate: arithmetic from the proof's per-step costs and the harnesses' sizes, not runs; the
  simulation-based calibration figure assumes the logistic arm's sizes.
- The saddle-point sampler's cost (WPS Table 1) and BayesLogit's switching points (the critique's reading).
