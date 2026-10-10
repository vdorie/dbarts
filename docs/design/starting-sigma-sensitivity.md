# Starting sigma: how much the fit depends on it, and what to do for a large sparse design

Status: FINDING, 2026-10-10 (dec-B403, dec-B421, dec-B422). A sensitivity study on dense data and an
LSQR prototype for sparse designs too large for an exact fit. Measured on arm64 macOS (M1 Max, R 4.6.1,
R's reference BLAS), two jobs at a time. Scripts and outputs are in scratch/ssens, untracked.

## Summary

The starting sigma (sigest) is the residual sd of a linear regression. It does two jobs: it calibrates the
default residual-variance prior, and it is the value the chain starts from.

**The fit depends on sigest at n = 200, hardly at n = 1000 and not measurably at n = 5000.** At n = 200 a
sigest wrong by a factor of 2 to 4 moves the posterior of sigma by as much as 19 to 57 points of the true
value and the coverage of a 90 percent predictive interval by as much as 3 to 6 points; a sigest that is
too small also costs up to 9 percent in RMSE where p is near n. At n = 1000 the same errors move coverage
by under half a point and sigma by under 2 points, a fourfold overestimate excepted (1 point and 6 points).
At n = 5000 no arm differs from the default by more than a refit of the default with another seed does,
give or take 0.2 points of coverage.

**All of it is the prior.** Starting the chain elsewhere under the same prior changes nothing.

**An iterative least-squares solver (LSQR) gives the sparse estimate at any size in seconds.** With a
design's known structure it agrees with the exact routine to 4e-5 or better; it takes 30 seconds where the
exact routine cannot allocate (m = 1e5, m the smaller of n and p + 1), and 1 to 2 seconds for cross
validation's 200 per-fold estimates at n = 1e4, p = 1000, where the exact routine takes 30.

**Recommendation.** Above the cap, use LSQR, silently. And lower the cap from m = 10,000 to about 2,000:
exactness matters only where n is small, and there m is small and the exact routine is fast. See
[Recommendation](#recommendation).

## Part 1: sensitivity

### Setup

- Data: Friedman's five-variable test function on uniform predictors, the other p - 5 columns noise.
  n = 200, 1000, 5000; p = n / 20 or 0.9 n; signal-to-noise (variance of the function over the noise
  variance) 0.5, 2 and 8; 75 trees (the default) and 200. Everything else is bart()'s default: four chains,
  500 burn-in sweeps and 500 kept each.
- Arms: sigest set to 0.25, 0.5, 1, 2 and 4 times the linear-model estimate, and to sd(y). The same data
  and the same chain seed in every arm.
- Yardstick: the 1x arm refit with another chain seed. A difference between arms no larger than the
  difference between those two fits is seed noise.
- 12 data sets at n = 200, 10 at n = 1000, 4 at n = 5000; 1000 test points each.
- Under this function the linear estimate is itself 1.0 to 1.1 times the true sigma at low signal, 1.2 at
  medium and 1.7 at high (the regression leaves the nonlinear part in its residual), and sd(y) is 1.2, 1.7
  and 3.0 times. So the arms span 0.25 to 7 times the truth.

### Results

Largest change against the 1x arm in any of the 12 cells at that n (mean over data sets; "seed" is the
yardstick's largest):

| measure | n | 0.25x | 0.5x | 2x | 4x | sd(y) | seed |
|---|---|---|---|---|---|---|---|
| RMSE against the true function, percent | 200 | 9.3 | 5.8 | 4.8 | 6.8 | 4.1 | 2.2 |
| | 1000 | 2.2 | 1.8 | 1.9 | 3.6 | 2.1 | 2.9 |
| | 5000 | 4.8 | 7.5 | 6.4 | 8.2 | 6.7 | 9.4 |
| coverage of a new observation, 90 percent interval, points | 200 | 6.2 | 3.0 | 2.6 | 5.1 | 1.7 | 0.4 |
| | 1000 | 0.3 | 0.2 | 0.3 | 1.0 | 0.2 | 0.3 |
| | 5000 | 0.2 | 0.4 | 0.4 | 0.7 | 0.5 | 0.5 |
| posterior mean of sigma, points of the true value | 200 | 33 | 19 | 23 | 57 | 18 | 1.2 |
| | 1000 | 0.8 | 1.2 | 1.8 | 5.7 | 1.6 | 0.8 |
| | 5000 | 3.7 | 2.6 | 2.2 | 3.2 | 2.2 | 4.5 |

The direction, in the most sensitive setting (p = 0.9 n, high signal, both tree counts pooled), as the
change against 1x:

| n | measure | 0.25x | 0.5x | 2x | 4x | sd(y) | seed |
|---|---|---|---|---|---|---|---|
| 200 | sigma, points | -21.6 | -13.0 | +19.6 | +49.3 | +15.0 | -1.0 |
| | coverage of a new observation, points | -3.9 | -1.7 | +2.1 | +4.3 | +1.4 | 0.0 |
| | coverage of the function, points | -2.6 | -0.8 | +1.1 | +1.3 | +0.5 | +0.3 |
| | RMSE, percent | +6.0 | +2.5 | -2.2 | -1.9 | -1.5 | -0.4 |
| 1000 | sigma, points | +0.2 | -0.2 | +1.3 | +4.8 | +1.4 | +0.5 |
| | coverage of a new observation, points | -0.2 | -0.1 | +0.1 | +0.4 | +0.1 | -0.2 |
| | RMSE, percent | +1.6 | +0.8 | +0.4 | +0.5 | +1.2 | +2.5 |
| 5000 | sigma, points | -1.5 | -0.5 | -0.4 | -0.7 | -0.5 | -1.9 |
| | coverage of a new observation, points | -0.1 | +0.1 | 0.0 | 0.0 | -0.1 | -0.2 |
| | RMSE, percent | -0.8 | -1.0 | 0.0 | -0.3 | +0.6 | -0.8 |

In plain words:

- **n = 200.** sigest matters, most at p near n, high signal and 200 trees. There the posterior mean of
  sigma is 0.24, 0.37, 0.56, 0.79 and 1.13 of the truth at 0.25x to 4x (0.75 at sd(y)), and the predictive
  interval covers 0.86, 0.89, 0.92, 0.94 and 0.96 (0.93 at sd(y)). A sigest that is too small is the
  harmful direction: the forest fits noise, sigma collapses, intervals undercover and RMSE rises 4 to 9
  percent. A sigest that is too large widens intervals by a few points and, at p near n, lowers RMSE by 1 to
  3 percent. With p well below n a sigest that is too small does little (sigma down 4 points, coverage
  under 1) while one that is too large still widens intervals (at 4x and high signal, sigma up 26 points
  and coverage up 4); RMSE does not move (under 2 percent).
- **n = 1000.** Only the 4x arm is visible: sigma up 1 to 6 points, coverage up 0.1 to 1.0. Halving,
  quartering or doubling sigest, or using sd(y), stays within seed noise on every measure.
- **n = 5000.** Nothing is visible: the largest candidate is the 4x arm's coverage, up 0.7 points in one
  cell where the yardstick's largest is 0.5. Seed noise itself is large here (RMSE differs by 3 to 12 percent
  between two seeds of the same fit), because 500 burn-in sweeps are not enough:
  in 19 percent of fits at least one chain's sigma had not reached its kept range by the end of burn-in,
  the same share in every arm (15 to 25 percent).
- **sd(y) in place of the linear estimate**, the fallback in question: at n = 200, sigma up 6 to 15 points
  and coverage up 1.1 to 1.4 points at high signal and at medium signal with p near n (under 1.5 and 0.4
  elsewhere), RMSE 1.5 percent lower at p near n; at n = 1000 and 5000, within seed noise.
- **Sweeps for sigma to reach its stationary range** (the first burn-in sweep inside the central 90 percent
  of the kept draws, median over chains and data sets): at n = 200, 1 to 11 with low signal and 6 to 114
  with high; at n = 1000, 1 to 264; at n = 5000, 29 to more than 500. It is set by n and the signal, not by
  sigest. The one visible arm effect is at n = 200, p near n, high signal, 200 trees: 114 sweeps at 0.25x,
  35 at 1x, 16 at 4x, because under a prior pulled low the forest has further to grow.
- **Prior or starting value?** At n = 200, p = 180, 200 trees, high signal (8 data sets), with the prior
  calibrated at one multiple and the chain started at another: posterior mean sigma over truth is 0.55 with
  both at 1x, 0.58 and 0.55 with only the start moved to 0.25x and 4x, and 0.24 and 1.12 with only the
  prior moved. The start is forgotten; the prior is the whole effect.

### What the dense path does at p >= n

A dense design with p + 1 >= n leaves the regression no residual degrees of freedom; the estimate is then
sd(y), with a warning. Just below that the estimate is used as it comes, with no warning, however few
degrees of freedom it rests on: on k of them its relative sd is about 1 / sqrt(2k), 16 percent at the 19
the n = 200, p = 180 cells had and 32 percent at 5.

## Part 2: LSQR

Sized by Part 1 as a feasibility check: above the cap n is over 10,000, where Part 1 says the fit cannot
tell an exact estimate from a rough one.

### The routine

LSQR (Paige and Saunders) on the regression of the response on an intercept and the design, in R, with
Matrix's sparse products only; Matrix is already imported, so there is no new dependency. Conjugate
gradients on the normal equations is the same iteration in exact arithmetic and less stable in floating
point, so it was not run separately. About 100 lines.

- Every column is centered and scaled to unit norm inside the operator, so the design is never copied or
  densified and the intercept is orthogonal to every column. Constant columns are dropped.
- Residual degrees of freedom: n less the intercept and the columns, plus one for each full indicator
  block (its columns sum to the intercept). Dependencies not known from structure are not found.
- Below 10 percent of n residual degrees of freedom the routine returns sd(y) without iterating.
- It stops when the normal-equations residual is below 1e-6 of the product of the operator and residual
  norms, or at an iteration cap. The tolerance is floored at 1e-12: run past rounding-level convergence
  the recurrence drifts (a design converged at 36 iterations was 138 percent off at 100).
- The reported sigma comes from the residual of the returned coefficients, computed directly.

**It can only overestimate.** Any coefficient vector's residual sum of squares is at least the minimum, so
stopping early errs upward, and so does a degrees-of-freedom count that misses a dependency. Stopping at a
fixed iteration, on the slowest designs tried (15 percent residual degrees of freedom): +27 to +33
percent at 5 iterations, +9 to +10 at 10, +0.5 to +0.6 at 25, +0.01 at 50, under 1e-7 at 100. On a
well-conditioned design: 2e-6 at 5. A missed dependency count d overestimates by the factor
sqrt((df + d) / df).

### Accuracy against the exact routine

Eleven designs under the cap, three draws each (one with weights, one with an offset); n = 5000 unless
noted. Largest relative difference in sigma and largest iteration count at the default tolerance:

| design | tolerance 1e-4 | 1e-6 | 1e-8 | iterations at 1e-6 |
|---|---|---|---|---|
| sparse numeric, p 500 | 3e-8 | 8e-12 | 2e-15 | 15 |
| one-hot, 20 factors x 50 levels | 9e-8 | 1e-11 | 3e-15 | 22 |
| one-hot, skewed levels, a third of them empty | 4e-8 | 1e-11 | 3e-15 | 23 |
| one-hot with two reference levels dropped | 3e-7 | 1e-11 | 4e-15 | 37 |
| one-hot, numeric, a timestamp column and a large-mean column | 2e-8 | 5e-10 | 1e-10 | 15 |
| 50 near-copies, differing at 1e-3 | 1.1e-2 | 4e-5 | 6e-9 | 69 |
| n 2000, p 1700 numeric (15 percent residual df) | 2e-5 | 4e-9 | 5e-13 | 254 |
| n 2000, 34 factors x 50 levels (17 percent) | 1e-5 | 3e-9 | 4e-13 | 231 |
| 50 exact linear combinations, not structural | 5.4e-3 | 5.4e-3 | 5.4e-3 | 18 |
| 50 near-copies, differing at 1e-6 | 5.5e-3 | 5.5e-3 | 5.3e-3 | 15 |
| n 2000, p 1840 (8 percent residual df) | sd(y) | sd(y) | sd(y) | 0 |

- With the structure known, the default tolerance agrees to 4e-5 or better.
- The 0.5 percent rows are the degrees of freedom, not the solver: 50 dependencies the structure does not
  show, counted as 50 fitted columns out of 4550 residual degrees of freedom. The second of them is the
  accepted rank band (dec-B421): the exact routine drops the near-copies and LSQR counts them.
- The last row is the 10 percent rule: sd(y) is 15 to 164 percent above the exact estimate there. The exact
  routine has no such rule; it gives an estimate down to one residual degree of freedom.

### Cost and memory

n = 1.5 m and p = m, so a third of n is left as residual degrees of freedom; numeric columns at the stated
density, or full one-hot blocks with 1 / density levels. Seconds and iterations at tolerance 1e-6:

| m | stored entries at 0.1 / 1 / 5 percent | numeric, seconds | one-hot, seconds | iterations |
|---|---|---|---|---|
| 1e4 | 1.5e5 / 1.5e6 / 7.5e6 | 0.04 / 0.15 / 0.97 | 0.04 / 0.19 / 0.92 | 46 to 56 |
| 5e4 | 3.75e6 / 3.75e7 / 1.9e8 | 0.5 / 5.4 / 43 | 0.7 / 6.3 / not run | 51 to 54 |
| 1e5 | 1.5e7 / 1.5e8 / 7.5e8 | 2.2 / 29 / not run | 3.0 / 29 / not run | 51 to 55 |

- Cost is linear in the stored entries and the iteration count: 3e-9 to 4e-9 seconds per entry per
  iteration. The points not run hold 2 to 9 GB of design; m = 1e5 at 5 percent would take about 2.5
  minutes at that rate.
- The iteration count rises as the residual degrees of freedom shrink: at m = 5e4 with 11 percent left
  (n = 56,000), 173 iterations and 11 s for numeric columns at 1 percent, 97 iterations and 1.4 s for
  one-hot blocks of 500 levels.
- Memory: the peak R heap, design included, was 1.5 to 2.8 times the design's own size (2.5 GB numeric and
  4.7 GB one-hot at m = 1e5 and 1 percent, where the design is 1.7 GB).
- The exact routine, for comparison: 14 s at m = 5000, 2 minutes and 1.8 GB at 1e4, about an hour at 2e4,
  7 hours or more and 40 GB at 5e4, an allocation failure at 1e5.

### Cross validation's per-fold estimate

xbart re-estimates on each fold's training rows: 40 replications of 5 folds, 200 estimates, by default.
Held-out rows are passed as zero weights, so no fold copies the design.

| design | all rows | one fold | 200 folds, LSQR | 200 folds, exact routine |
|---|---|---|---|---|
| n 1e4, p 1000, numeric at 1 percent | 0.06 s | 0.005 s | 1 s | 28 s |
| n 1e4, one-hot 20 x 50 | 0.06 s | 0.01 s | 2 s | 32 s |
| n 7.5e4, p 5e4, numeric at 0.1 percent | 0.6 s | 1.0 s | 3.4 min | not possible |
| n 7.5e4, one-hot 500 x 100 (3.75e7 entries) | 6 s | 9 s | 30 min | not possible |

A fold costs more than all rows at the large sizes because its training rows leave fewer residual degrees
of freedom (17 percent against 33), which raises the iterations 1.6 to 2 times. Starting each fold from the
all-rows coefficients saves 2 to 8 percent, not worth carrying. At the present cap the exact routine's 200
folds would take 7 hours at m = 1e4 and 47 minutes at m = 5000.

### Simpler alternatives

On n = 8000, p = 5000 (12 designs: two densities, sparse or dense signal, three draws), with screening to
the 2000 columns most correlated with the response and then the exact routine:

| method | difference from the exact estimate | cost |
|---|---|---|
| LSQR | 1e-10 | 0.02 to 0.08 s |
| screening, then exact | -11 to +7 percent | 0.9 s here; the exact routine's 2 minutes at 10,000 columns |
| sd(y) | +8 to +61 percent | none |
| exact on every column | 0 | 14 s |

Screening errs in both directions: low when the signal is sparse (-7 to -11 percent; columns chosen for
their correlation with the response also fit its noise) and either way when every column carries some
(-9 to +7). It costs as much as the cap allows on every call and is the least predictable of the three.
sd(y) is free and always an overestimate, by the whole linear signal: 1.15 to 1.75 times the linear
estimate in Part 1's settings.

## Recommendation

1. **Above the cap, use LSQR, with no warning**: structural degrees of freedom, tolerance 1e-6, a cap of
   1000 iterations, sd(y) below 10 percent residual degrees of freedom. It keeps dec-B403's promise (a
   sparse design gets its dense equivalent's estimate) to a percent or better on every design measured, it
   errs only upward, it runs in seconds to a few minutes where the exact routine takes hours or cannot
   run, and its remaining error is far below anything Part 1 says a fit with n over 10,000 can show.
   Screening is dominated. sd(y) alone would also pass Part 1 at that n, but it gives up dec-B403 to save
   seconds.
2. **Lower the cap to about m = 2000.** Part 1 says the exact routine's exactness is worth something only
   at small n, and m is never more than n, so a design with m over 2000 has n over 2000. At n = 1000, the
   nearest size measured below that, a sigest off by a factor of 2 is already within seed noise, and
   LSQR's error is a percent or less on the designs measured. Between m = 2000 and 10,000 the exact
   routine costs 1 s to 2 minutes a call and 350 MB to 1.8 GB, and up to 7 hours across xbart's 200 folds;
   LSQR costs under a second, with memory in proportion to the design's stored entries. At m = 2000 the
   exact routine's 200 folds take about 3 minutes.
3. **Keep the exact routine below the cap.** There the rank it finds matters (a wide design at small n can
   hide many dependencies, and LSQR's count would overestimate by sqrt((df + d) / df)), and there Part 1
   says sigest matters.

Two things this leaves for the maintainer. The 10 percent rule is new behavior above the cap only: the
exact routine and the dense path estimate down to one residual degree of freedom, and whether they should
is a separate question (Part 1's noisiest linear estimates, 19 degrees of freedom at n = 200, are where a
low draw hurts). And the dense path has no cap at all, though its regression costs more than the exact
sparse routine at the same size; the same solver would serve it, and that was not measured here.

Limits of the evidence: one test function, gaussian responses, unweighted fits, n no larger than 5000 in
Part 1 (larger n was not run; the trend from 200 to 5000 is monotone), and random sparse designs in the
cost tables, which converge in about 50 iterations; a badly conditioned real design would take several
times as many, at the same cost per iteration.
