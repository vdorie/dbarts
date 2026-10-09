# k under probit: the calibration finding and what the fix must do

Status: FINDING, 2026-10-08 (dec-B371). The fix is the probit rescaling
step, [probit-k-scale-move.md](probit-k-scale-move.md); this note is the
evidence it starts from, and its section on what a correct sampler gives
describes the sampler before that step.

## The symptom

The probit-k arm of the calibration suite (benchmarks/R/sbc.R; k drawn
under the binary default chi(1.5, 2), R = 600 replications, 99 kept draws
at thin 1000, ranks in 20 bins) passes its band but puts too many of k's
ranks in the end bins on every full run: 48 in the lowest bin on the
laptop, 52 in the highest on CI, against 30 expected (chi-square p 0.012
and 0.013).

## What it is

The sampler is correct; k mixes too slowly for the thinning.

- **The k step draws its exact conditional.** Each recorded k, scored
  against the gamma conditional of k^2 given that sweep's leaves, is
  uniform over 24000 draws on six datasets (KS p 0.18), including datasets
  where mixing is slowest.
- **The arm draws from the sampler's model.** With every row inactive the
  sampler targets the prior: its trees average 73.07 splits against the
  generator's 73.06 (sd 5.90 and 5.89), and its k follows chi(1.5, 2)
  (KS p 0.53). The leaf scale and the cut grid agree.
- **Started at the truth, the chain does not drift.** A replication's true
  state is an exact posterior draw of its own data, so a correct sampler
  started there stays at the posterior at every lag, however slowly it
  mixes. Over 9080 replications, k at a fixed lag minus the true k is
  +0.007 in log k (se 0.004). The planted error the arm exists to catch, a
  half-unit error in the shape of k's conditional, shows +0.135 (z 5.6 at
  300 replications). The same holds within groups chosen by the data alone
  (slow chains, near-separated responses).
- **The slow part.** k's integrated autocorrelation time is a median of 500
  sweeps, 7100 at the 90th percentile and 25000 at the 99th, rising as the
  true k falls; a third of replications keep fewer than 90 effective draws
  of their 99. Those replications carry the excess; the rest rank
  uniformly (p 0.91).
- **Separation.** Below k of about 0.135 the forest's fit is large enough to
  separate y. That region is real posterior mass: on a one-tree problem
  with a known answer, P(k < 0.135 | y) = 0.103, the sampler gives 0.102 to
  0.107 from exact starts, and free chains of four million sweeps reach
  only 0.046 to 0.079. A chain that starts there stays for 100k sweeps; one
  that starts outside rarely enters (2 of 600 arm chains walked in from
  ordinary starts). This is the known slow mixing of probit's latent
  augmentation under separation.

Thinning does not cure it: the end bins fall from thin 100 to thin 300,
then stop falling; the slowest chains would need a thin near 25000.

## What a correct sampler gave before the rescaling step

Pooled over 2180 replications of the arm's own set-up, scaled to R = 600:
end bins averaging 47 and 47 against 30; the larger end bin a median of
50, 58 at the 90th percentile, 64 at the 99th (above 70 about one run in
a thousand); chi-square p with quartiles 3e-4, 0.004 and 0.03. The arm's
flag (p below 6e-4) therefore fires on about one correct run in three.
Until the fix lands the arm is read by its end bins (benchmarks/README.md):
defective only where an end bin passes 70 or the two ends differ by more
than about 25 (the planted error put 86 in the lowest bin).

## What the fix must do

Every binary fit draws k by default, as 0.9-34's bart2 did, so this is a
user-facing mixing problem, not a property of the check. The fix is one
cheap move a sweep that changes the fit's overall scale in one step:
interweaving on k (a non-centred draw of k given standardized leaves) or
a probit parameter expansion (a working scale on the latents and leaves).
The candidates are compared on k's autocorrelation time over the arm's
datasets and on how often free chains reach the small-k region of the
one-tree problem against its known 0.103; the winner carries an exact
gate. Once it lands the arm should rank k uniformly at a much smaller
thin, its flag back in force.
