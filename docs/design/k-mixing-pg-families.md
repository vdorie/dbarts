# k under logistic and negative binomial: how it mixes

Status: FINDING, 2026-10-09. The choice of a fix for logistic is the
maintainer's (TODO k-mixing-pg-families); this note is the evidence.

## The question

Under probit, k (the leaf prior's scale, drawn by default for binary
fits) mixed slowly where the true k is small or the response nearly
separated (probit-k-calibration.md). The probit fix rescales the latent
variables, the leaves and k together (probit-k-scale-move.md). Logistic
and negative binomial fits also draw k by default but use Polya-Gamma
latents, which are precisions rather than shifted responses, so that
step does not apply. Do they need a fix?

## How it was measured

The same generator for all three families: 150 rows, 3 predictors, 50
trees, k drawn from the default chi(1.5, 2) prior, trees from their
prior, the response from the model. Logistic dataset r shares probit
dataset r's true k. Each chain ran 30,000 burn-in sweeps and 200,000
recorded at thin 10 (negative binomial 5,000 and 50,000, its sweeps
costing about ten times more; its shape drawn from the grid prior and
estimated). The autocorrelation time tau is in sweeps. Probit numbers
are the sampler before the rescaling step.

A second check started each chain at the true trees, leaves and k (an
exact posterior draw once the latents are drawn given the response) and
asked whether k drifts away from its prior over up to 10,000 sweeps.

## Results

tau of k, in sweeps:

| family | datasets | median | 90th | 99th | worst | share above 5000 |
|---|---|---|---|---|---|---|
| probit, before the step | 100 | 524 | 6108 | 23278 | 52751 | 0.11 |
| logistic | 100 | 272 | 2928 | 27505 | 50605 | 0.07 |
| logistic, near-separable (true k < 0.4) | 20 | 5242 | 27521 | 32966 | 33610 | 0.50 |
| negative binomial | 84 of 100 | 191 | 611 | 3253 | 3947 | 0 |
| negative binomial, low counts | 27 | 286 | 712 | 2726 | 3063 | 0 |

Median tau of k by the true k:

| true k | probit | logistic | negative binomial |
|---|---|---|---|
| up to 0.3 | 12797 | 26837 | 3947 |
| 0.3 to 0.6 | 4744 | 3553 | 1833 |
| 0.6 to 1 | 2711 | 1370 | 481 |
| 1 to 1.5 | 1749 | 625 | 356 |
| above 1.5 | 293 | 169 | 159 |

- Logistic is slow on the same datasets as probit (log-tau correlation
  0.67; logistic's tau is a median 0.58 of probit's), and slows with
  separation as probit does: tau 121 where the true fit's sign matches
  the response on up to 70% of rows, 3798 above 95%.
- Where k is stuck, so is the average fit (logistic near-separable:
  median tau of the average fit 4330).
- Negative binomial's slow cases are mostly-zero data (more than 60%
  zeros: tau about 2000) and very small true k. Its cost per sweep grows
  with the counts, so its effective draws of k per second (median 9)
  trail probit's (33) despite the lower tau.
- No drift from a start at the truth: logistic -0.010 in log k (se
  0.008) over 2,800 replications, its draws matching the prior at every
  lag; negative binomial +0.005 (se 0.009). The slow chains are slow
  mixing, not a wrong sampler.

## Gap

15 of the 100 negative binomial prior datasets were not run: their
counts sum past 10,000 (some past 10^9), and the latent draw's cost is
linear in the counts. They are mostly the small-k datasets, where k is
slow in every family, so negative binomial's tail is undermeasured; the
measured trend puts it at roughly 1,000 to 4,000 sweeps. That about 15%
of datasets drawn from the default prior have such counts bears on
TODO approximate-distributions as much as on k.

## What a fix would target

Logistic needs a move on the fit's overall scale in the small-k,
near-separated region, as probit did. The candidate is a collapsed
move: rescale every leaf and 1/k together and accept on the logistic
likelihood with the latents integrated out, one likelihood pass a
sweep; the same move covers negative binomial. Interweaving with the
latents as weights is weak in exactly this tail, as it was under probit.
For negative binomial alone, k is not the priority: the latent draw's
cost and the counts' heavy tail under the default prior are.
