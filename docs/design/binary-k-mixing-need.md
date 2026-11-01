# Does k have to move fast? The finding

## The answer

**Under logistic, yes, where the signal is strong: the shipped sampler's k does not get down to the small
values the data want within a default run, and intervals come out too narrow. On calibration, the collapsed
move removes that harm in every cell where it is shown. Under probit no harm from k is shown, though much is
unresolved.**

In the bulk of use (tier 1) no reading meets the rule for a problem, and 13 of 38 readings resolve under it.
That does not mean common use is clear. Reading every long run, which the rule does not do, sorts the 24 tier
1 readings left unresolved by failed long runs: 12 clear, 8 stay open and 4 sit at or past the 2-point bar. One
of the four looks like the same harm: the interior probability in the sampler loop under logistic, where the
shipped sampler is short by 3.5 points (se 1.3) over all its long runs and the move by 0.5, and where
calibration's difference from the fit told the true k is 1.9 points (0.8) per fit, 1.0 (0.5) pooled over rows.
The other three, stan4bart's group intercepts and group standard deviation with census regions and the
propensity, are shared by both samplers, so they are not k's. That reading was made after the results and
replaces no rule; the verdicts below do not use it. Under the rule itself, effects at moderate prevalence clear
in every tier 1 class, and stan4bart, which is probit, clears with its twins data.

Under probit, 45 of 137 readings clear (15 bounded, 3 measured, 27 agreeing with their one long run) and 92
are unresolved. Nothing there is shown to be k's, but not all of it is sound. Interior probabilities read 0.891 pooled over
rows and 0.879 (0.004) per fit, with the told fit at 0.884 and the move at 0.879, so about 2 points short per
fit in all three, and what a user sees on the default link. Ten times the length takes the told fit to 0.895 (0.005)
per fit, which points to the trees' cost at the default length and not the model's. Near 0 or 1 in near
separation, rare outcomes and many predictors the sampler under-covers by 6 to 10 points; the told fit does no
better, so the study cannot separate k from the trees there, and the move gains 3.4 points (se 1.4) near
separation and at most 1.3 elsewhere. The told fit climbs there with length, which also points to the trees. Two real datasets with very small k, mushroom and banknote, are the
exception: under probit too their long runs disagree in substance.

The largest miss in the study is not k's. Both samplers cover multiple treatments' contrasts at 0.57 to 0.72
against a nominal 0.90, on 10 replications of each of three designs. The move adds about 2 points under logistic
and nothing under probit, so this points to the trees and the model. Your colleagues' case did not run; there is
room for an amendment when it arrives.

The study does not say "k must mix fast". It says "logistic needs to reach small k". The simulated regimes
where that is shown are specific, but strong signal is not rare in real prediction data: 8 of 22 real
prediction datasets disagree in substance under logistic. It opens the pricing of the move, as a sequence with
each step's rule written before any fit: first that the move's small-k tail is right, or pricing stops; then that
it wins at equal time on coverage, on these regimes and on fits with no strong signal, and if so it is on for
every logistic fit, with no trigger; only if that costs too much, a check on the data; an opt-in if both fail.

The decisions this raised are ruled and listed at the end. Run ten times as long, the told fit reads on nominal under
every problem shown, and off nominal in two places under probit: inverse weights, which are clipped, so an exact interval
on them holds more than 0.90, and near 0 or 1 near separation, which comes within noise of nominal at forty times the
length. The plan's rule read off nominal as a fault that stops a regime. Neither reads as a fault here, so you ruled the
readings kept with a note and the rule reworded. Applied as first written to all of near separation, it would have
withdrawn three of the seven problems and left this answer as it stands.

## What was asked and how it was read

You ruled that k must move in binary fits, and asked whether any movement is enough. For each regime and
estimand, under each link, the study asked whether the shipped sampler at the default length gets interval
coverage wrong by at least 2 points, and whether the error is k's. Two measures decide. In *calibration*, the
truth is drawn from the model's own prior, so a sampler that reaches the posterior covers at nominal; a
companion fit that is told the true k runs the same tree moves for the same length, and a shortfall counts as k's only
where the shipped sampler falls at least 2 points behind it.
In the *long-run check*, a default-length interval is compared with a very long run that both samplers reach
and agree on: each long run has a shipped half and a move half, and passes acceptance when the two agree and
their precision checks hold. A cell's long-run reading counts only if at least 80 percent of its long runs pass
(the 80 percent rule) and it survives scoring each failed run against its move half. Each reading is one of:
problem shown; no problem, bounded (even a k frozen anywhere in its posterior could not cost 2 points); no
problem, measured (every available reading clears 2 points); agrees with its one long run (real data, below); or
unresolved.

**What "k's" means here.** Coverage alone names no cause, so a shortfall the told fit shares goes to the trees or the
model and is not counted against k, however large. All seven problems are on probabilities and bands on them; no effect
is among them. "k's" means that getting k right closes the gap at the default length with today's tree moves. The move
shows that a second way: it rescales the leaves and k with the trees held, and closes the same gaps. It does not clear
the trees. A wrong k has the trees searching under the wrong leaf scale, and k's draws depend on the leaves the trees give
it, so the gap holds both, and whether better tree moves would close it too was not run. Calibration's truths come from
the model's own prior, probably the trees' easiest target; multiple treatments, with designed truths, have no told fit.
Near 0 or 1 under probit the split fails, since fixing k also turns off the rescaling step, so the told fit there is not
the shipped sampler less k. The long-run check reads total error at the default length, trees included: it can clear a
reading, and charges k only through the move's test or a chain held at a fixed k.

| | bounded | measured | agrees with its one long run | problem shown | unresolved |
|---|---|---|---|---|---|
| tier 1, the bulk of use | 9 | 4 | 0 | 0 | 25 |
| tier 2, common in health data | 11 | 0 | 0 | 1 | 25 |
| tier 3, niche | 0 | 1 | 0 | 4 | 13 |
| real data | 0 | 0 | 46 | 0 | 18 |
| simultaneous bands | 0 | 0 | 1 | 2 | 48 |
| designs and datasets whose first long run failed |  |  |  |  | 30 |
| all | 20 | 5 | 47 | 7 | 159 |

A real dataset has one long run and no known truth, so at best its default interval agrees with that run. By the
written rule those 47 read "no problem, measured", for 52 in all. They are counted apart because they carry no
truth, which changes one count: measured falls from 52 to 5. 46 are real data's and one, lindner's band on
individual effects under logistic, sits in the bands row. Bands are apart because only calibration can read them,
and at these counts it rarely can.

One-value estimands read on long runs alone, as the plan says and the rule as built did not: one interval a fit
cannot resolve under calibration at 400 replications. That moves two readings to measured, the average effect at
moderate prevalence under probit (shortfall -0.001, upper bound 0.005) and stan4bart's group standard deviation
with twins (0.002, 0.019), both on 6 of 6 long runs. Calibration's figure stays beside them, and a large problem
there would still count against them.

Why the 159 are unresolved:

| | not run | a band | too few counts | failed, clear on all runs | failed and suspected | failed, still open | all |
|---|---|---|---|---|---|---|---|
| tier 1, the bulk of use |  |  | 1 | 12 | 4 | 8 | 25 |
| tier 2, common in health data |  |  | 1 | 10 | 4 | 10 | 25 |
| tier 3, niche |  |  | 4 |  | 8 | 1 | 13 |
| real data |  |  | 18 |  |  |  | 18 |
| simultaneous bands |  | 48 |  |  |  |  | 48 |
| designs and datasets whose first long run failed | 30 |  |  |  |  |  | 30 |
| all | 30 | 48 | 24 | 22 | 16 | 19 | 159 |

- **Not run**: the first long run failed, which stops that design's or dataset's other fits, so nothing can be read.
- **A band**: the long-run check has no reading for a band, and calibration's counts are too few.
- **Too few counts**: the long runs passed, and more runs or replications would settle it. 18 of the 24 are real
  datasets at one run each.
- **Failed**: too many long runs failed acceptance for the rule to read the cell, sorted by the reading of every
  long run below. *Clear* means the worse half's upper bound is under the bar; *suspected* means the worse
  half's estimate is at or past the bar; *still open* is neither. Near 0 or 1, all 8 failed readings are at or past
  the bar on that arithmetic and count as suspected only, since tails stay under the rule. The largest inverse
  weight's two readings count as open, being pinned at their bound.

The study took about 69 core-hours against a stop of 211, about 4.5 of them the long told fit below.

Two checks qualify these readings.

**Is the told fit a sound comparator?** For every problem shown, yes. The plan checks it by running the told fit ten
times as long: on nominal then means the default's gap is the trees' cost; off nominal, if a longer run leaves it off,
means a fault in the study's setup, which stops that regime's readings. Every regime owed ran, about 4.5 core-hours for the last of them, on 25
replications of each design and link (20 for rare outcomes at 10,000 rows under probit, which has no more), and 100 near
separation under probit. The sampler loop did not run, being owed only by a band. Off nominal means past 3.44 standard
errors, the noise bound for the 86 readings; without the 17 bands the bound is 3.38 and the same two readings pass it.

The problems' comparators hold. Under logistic the told fit at ten times the length reads 0.916 and 0.906 near
separation, 0.889 and 0.916 in small samples and 0.885 and 0.906 for rare outcomes at 500 rows, interior and near 0 or
1, all within noise of nominal. Since near separation under probit moved on going from 25 replications to 100, the same
extension was run on these three under logistic after the results were in, and none is past the bound: over 100
replications (fewer in the interior, which not every fit has), interior and near 0 or 1, 0.908 (0.008) and 0.900 (0.007)
near separation, 0.896 (0.007) and 0.907 (0.007) in small samples and 0.883 (0.009) and 0.894 (0.007) for rare outcomes
at 500 rows, with the two problem bands at 0.880 (0.037) and 0.883 (0.041). These and the two further runs below ran on
a second machine from the same seeds, which matched the first exactly on the 27 replications compared.

Of the 69 readings without a band, 67 are within the bound. 32 of those are one-value estimands, whose yes or no at 25
replications says little (the largest is 1.7 standard errors). Two are past it, and neither reads as a fault.

- **Inverse weights, propensity and weights, 0.929 (0.007), +3.9 standard errors.** The weights are clipped at a
  propensity of 0.025 and 0.975, and a central interval on a clipped quantity holds more than 0.90 even when it is
  exact, since a value sitting on the bound is inside the interval whenever the interval reaches it. On these
  replications 44 percent of the true weights sit on a bound, where the interval covers 0.991 of them, against 0.880
  elsewhere. The fit's own draws put 0.943 (0.005) of their mass inside their interval, which is what an exact sampler's
  coverage averages to, so 0.929 is 1.4 points under its own target and within noise. It does not move with length
  (0.929 on these 25 replications at both lengths, 0.934 over all 100 at the default), and the shipped sampler and the
  move read the same, 0.935. The propensity itself, on the same rows with no clipping, reads 0.881 (0.009), within
  noise, climbing from 0.850.
- **Near 0 or 1, near separation, under probit, 0.854 (0.013), -3.7 standard errors over 100 replications.** The told
  fit climbs with length and has not finished at ten times: 0.741 at the default, 0.797 with four times the draws and
  0.854 at ten. The shortfall sits where the true k is small, and a told fit should cover at nominal at every true k, so
  this grouping is fair to it. With a true k of 0.5 or more, 57 replications, it reads 0.886 at the default and 0.899
  (0.009) at ten times. Under 0.5, 43 replications, it reads 0.573 and 0.803 (0.023). A small k makes latent values
  large, and the told fit has no rescaling step to reach them, since fixing k turns that step off, the likely reason the
  shipped sampler beats it at the default (0.797 on the same replications). To separate slow trees from a fault, the
  same fits were run on the same seeds at forty times the length, after the results were in and under no rule of the
  plan: 0.885 (0.008), -1.8 standard errors, inside the bound, with 0.899 (0.009) at a true k of 0.5 or more and 0.869
  (0.013), still climbing, under it. Scaled squared errors, which should average 1, fall from 12.3 at the default to
  1.56 and then 1.10. A mismatch between the prior draw and the fit would stay at any length and would not spare the
  larger true k. So this reads as the trees' cost at more than ten times the default, a case the plan's two did not
  foresee. Interior probabilities in that regime read 0.888 (0.006) at ten times and 0.889 (0.006) at forty, within the
  bound. Whether the reading under a true k of 0.5 reaches nominal was not run, so a small fault at a small k is not
  excluded; no reading rests on it, since nothing in this regime under probit resolves.

Read as the plan first wrote it, either reading stops its regime. You ruled that not applied (under "What was decided"),
and the plan now asks that a reading also fail at a longer run. Neither regime's readings under probit change a count, since all of them are unresolved
already. Stopping near separation under both links would withdraw three of the seven problems, giving bounded 20,
measured 52, problem 4, unresolved 162. The answer would rest on small samples and rare outcomes at 500 rows and would
not change. Logistic's told fit there passes the check on its own.

Taken together the 37 readings in the table are not clean: their squared z-scores sum to 79.6 on 37 (p about 0.00006,
treating them as independent), with a weighted mean of 0.897 (0.002). Without the two above, 50.8 on 35 (p about 0.04,
or 0.11 giving each z the heavier tails of its few replications), with a weighted mean of 0.896 (0.002). So beyond those
two the excess is slight, and it is mostly small gaps in interior probabilities, such as 0.881 with many predictors
under logistic and in stan4bart with twins. The noise bound was set after the first results, and a regime read at 25
replications can miss what 100 find, as near separation under probit did (0.829 at 25). Bands are aside by the rule: 14
of the 17 read under nominal at ten times the length, with a mean of 0.824, and 0.693 (0.053) near separation under
probit.

The same shape runs through the other regimes near 0 or 1. At ten times the length, over every regime with a long told
fit, a true k of 0.5 or more reads 0.893 (0.005) under probit and 0.894 (0.010) under logistic. Under 0.5 it reads 0.818
(0.019) over 55 fits under probit and 0.895 (0.012) over 23 under logistic, leaving out rare outcomes at 10,000 rows,
where three fits at a small true k cover 0.17 under probit and 0.61 under logistic and carry that regime's low figures
(0.849 and 0.863; without them 0.891 and 0.892). So under probit a told fit at a small true k is a weaker comparator
than its guarantee suggests, in the direction that understates any gap it is compared with. Under logistic, where the
problems are, it holds.

| regime | link | estimand | told, default length, all its replications | the same, on the long run's replications | told, ten times as long | replications | standard errors from nominal |
|---|---|---|---|---|---|---|---|
| near separation | logistic | interior | 0.897 (0.008) | 0.904 | 0.916 (0.011) | 18 | +1.4 |
| near separation | logistic | near 0 or 1 | 0.883 (0.012) | 0.871 | 0.906 (0.015) | 25 | +0.4 |
| near separation | probit | interior | 0.877 (0.007) | 0.877 | 0.888 (0.006) | 76 | -1.9 |
| near separation | probit | near 0 or 1 | 0.741 (0.027) | 0.741 | 0.854 (0.013) | 100 | -3.7 |
| many predictors | logistic | interior | 0.889 (0.006) | 0.875 | 0.881 (0.007) | 24 | -2.6 |
| many predictors | logistic | near 0 or 1 | 0.874 (0.010) | 0.870 | 0.901 (0.007) | 19 | +0.1 |
| many predictors | probit | interior | 0.898 (0.005) | 0.881 | 0.889 (0.007) | 21 | -1.6 |
| many predictors | probit | near 0 or 1 | 0.815 (0.030) | 0.870 | 0.891 (0.008) | 20 | -1.1 |
| moderate control | logistic | interior | 0.905 (0.005) | 0.901 | 0.903 (0.006) | 25 | +0.5 |
| moderate control | logistic | near 0 or 1 | 0.897 (0.010) | 0.888 | 0.892 (0.017) | 14 | -0.5 |
| moderate control | probit | interior | 0.898 (0.004) | 0.895 | 0.901 (0.005) | 25 | +0.2 |
| moderate control | probit | near 0 or 1 | 0.886 (0.007) | 0.873 | 0.888 (0.016) | 20 | -0.8 |
| effects, moderate | logistic | individual effects | 0.885 (0.012) | 0.863 | 0.866 (0.025) | 25 | -1.4 |
| effects, moderate | probit | individual effects | 0.905 (0.013) | 0.931 | 0.936 (0.017) | 25 | +2.0 |
| stan4bart, twins | probit | interior | 0.896 (0.004) | 0.880 | 0.881 (0.007) | 25 | -2.6 |
| stan4bart, twins | probit | near 0 or 1 | 0.888 (0.007) | 0.862 | 0.894 (0.009) | 24 | -0.6 |
| stan4bart, twins | probit | group intercepts | 0.910 (0.008) | 0.885 | 0.889 (0.024) | 25 | -0.5 |
| rare, 500 rows | logistic | interior | 0.884 (0.009) | 0.893 | 0.885 (0.021) | 22 | -0.7 |
| rare, 500 rows | logistic | near 0 or 1 | 0.891 (0.006) | 0.902 | 0.906 (0.009) | 25 | +0.6 |
| rare, 500 rows | probit | interior | 0.895 (0.009) | 0.911 | 0.909 (0.013) | 22 | +0.7 |
| rare, 500 rows | probit | near 0 or 1 | 0.819 (0.021) | 0.892 | 0.900 (0.011) | 25 | +0.0 |
| rare, 10,000 rows | logistic | interior | 0.906 (0.007) | 0.905 | 0.915 (0.006) | 25 | +2.4 |
| rare, 10,000 rows | logistic | near 0 or 1 | 0.840 (0.030) | 0.825 | 0.863 (0.027) | 25 | -1.4 |
| rare, 10,000 rows | probit | interior | 0.878 (0.011) | 0.878 | 0.899 (0.008) | 20 | -0.1 |
| rare, 10,000 rows | probit | near 0 or 1 | 0.794 (0.048) | 0.794 | 0.849 (0.042) | 20 | -1.2 |
| small n, strong signal | logistic | interior | 0.894 (0.004) | 0.888 | 0.889 (0.009) | 23 | -1.3 |
| small n, strong signal | logistic | near 0 or 1 | 0.890 (0.005) | 0.911 | 0.916 (0.013) | 25 | +1.3 |
| small n, strong signal | probit | interior | 0.889 (0.008) | 0.894 | 0.896 (0.015) | 20 | -0.3 |
| small n, strong signal | probit | near 0 or 1 | 0.876 (0.009) | 0.864 | 0.893 (0.014) | 25 | -0.5 |
| effects, rare | logistic | individual effects | 0.890 (0.019) | 0.908 | 0.894 (0.029) | 25 | -0.2 |
| effects, rare | probit | individual effects | 0.916 (0.012) | 0.884 | 0.888 (0.034) | 25 | -0.4 |
| propensity and weights | probit | individual effects | 0.897 (0.015) | 0.895 | 0.900 (0.022) | 25 | +0.0 |
| propensity and weights | probit | propensity | 0.864 (0.012) | 0.850 | 0.881 (0.009) | 25 | -2.2 |
| propensity and weights | probit | inverse weights | 0.934 (0.004) | 0.929 | 0.929 (0.007) | 25 | +3.9 |
| stan4bart, regions | probit | interior | 0.895 (0.006) | 0.901 | 0.901 (0.011) | 24 | +0.1 |
| stan4bart, regions | probit | near 0 or 1 | 0.848 (0.016) | 0.843 | 0.886 (0.014) | 25 | -1.0 |
| stan4bart, regions | probit | group intercepts | 0.870 (0.017) | 0.876 | 0.907 (0.025) | 25 | +0.3 |

Standard errors in parentheses. The default's column over all its replications and the long run's are on different
counts, and a difference between them can be which replications ran: rare outcomes at 500 rows under probit near 0 or 1
reads 0.819 over its 100 but 0.892 over the 25 with a long fit, against 0.900 at ten times. The climb is the column
beside it.

**Are accepted long runs right as a class?** Under calibration, long runs that pass should themselves cover
the truth at nominal. Checked regime by regime, 2 of 67 readings stand past the noise bound: the sampler loop's
interior under probit (0.958 over 6 runs) and individual effects at rare prevalence under probit (0.987 over 5).
Both over-cover, and on the same replications the shipped and told fits do too (0.958 and 0.939; 0.980 and
0.988), so it looks like those few replications and not acceptance passing bad runs. The plan's written rule
withdraws such a regime's long-run readings. This departs from it: the readings are kept with this note, as you
ruled. Kept, the counts are bounded 20, measured 52 (5 and the 47 above), problem 7, unresolved 159; withdrawn,
they would be 13, 55, 7, 163, with no problem appearing or going. The regime-by-regime check is weak: with so few
runs the noise bound is 7.3 standard errors at 6 runs, 9.3 at 5 and 14.3 at 4, so only an extreme can stand past
it, and that the other 65 do not says little. On a tighter bound, the normal one at 3.4 standard errors, one more
stands out: rare outcomes at 500 rows under logistic near 0 or 1 (0.941 over 4 runs), on which no reading rests.
The stronger check pools over regimes: accepted long runs cover the truth at 0.894 (0.017) under logistic and
0.915 (0.010) under probit, over row quantities, on 41 and 75 runs. That supports keeping them, and supports
acceptance as a class.

An independent critique went through these results from the stored records. It verified the seven problem
readings, how calibration was built and that selection never saw the truth, all 20 bounded readings, the 8 of 22
real datasets, the small-k figures of the long runs, and the 5 to 15 percent time figure, noisily. It inferred,
without verifying: that labeling rows from the shipped chain did not create the harm, from the move sitting within
noise of the told fit or above it in every problem cell; that a projection of running on to 40,000 sweeps is right;
that the group-level gap in stan4bart is the trees' or the group scale's rather than a fault in the setup; and that
moving k during burn-in alone could suffice. The long told fit, which came after the critique and which it did not see,
bears on the third of these and on nothing else in the list: census regions read on nominal at ten times the length.

## Where the problem is

All seven are logistic and all are shown by calibration. Coverage of a nominal 90 percent interval at the
default length, not at equal time, pooled over rows as the rule reads it, with the difference per fit beside it:

| regime | estimand | tier | shipped | told the true k | the move | shipped minus told (se) | the same, per fit (se) | move minus told (se) | replications |
|---|---|---|---|---|---|---|---|---|---|
| small samples, strong signal | interior probability | 2 | 0.872 | 0.894 | 0.893 | -0.022 (0.005) | -0.037 (0.007) | -0.001 (0.004) | 361 |
| small samples, strong signal | probability near 0 or 1 | 3 | 0.764 | 0.890 | 0.893 | -0.125 (0.014) | -0.096 (0.011) | +0.003 (0.008) | 400 |
| small samples, strong signal | band on probabilities | 3 | 0.696 | 0.833 | 0.867 | -0.138 (0.030) | | +0.033 (0.023) | 240 |
| near separation | interior probability | 3 | 0.838 | 0.897 | 0.896 | -0.059 (0.012) | -0.116 (0.024) | -0.001 (0.006) | 77 |
| near separation | probability near 0 or 1 | 3 | 0.691 | 0.883 | 0.899 | -0.192 (0.028) | -0.183 (0.027) | +0.015 (0.011) | 100 |
| near separation | band on probabilities | 3 | 0.533 | 0.733 | 0.760 | -0.200 (0.054) | | +0.027 (0.062) | 75 |
| rare outcomes, 500 rows | probability near 0 or 1 | 3 | 0.811 | 0.891 | 0.915 | -0.079 (0.021) | -0.079 (0.021) | +0.024 (0.011) | 100 |

The first row is the thinnest: about 2 points, possibly just under. It meets the rule at 400 replications only
with the look that triggered the extension: the first 100 replications read -0.034 (0.009), the fresh ones after
them -0.019 (0.006). Per fit it is larger and steadier, -0.040 (0.015) and then -0.036 (0.007), because fits with
few interior rows are the most harmed and count for more that way. Any weighting by the data keeps an exact
sampler on nominal, so the per-fit figure is as fair as the pooled one, though the rule reads the pooled. Pricing's
sequence does not weigh tiers, so this cell no longer sets the path. Bands still under-cover with the move (0.760
near separation) because the told fit does too.

Pooling over rows also hides harm elsewhere. In rare outcomes at 500 rows the interior under logistic reads +0.010
pooled and -0.046 (0.022) per fit, with the move at +0.030 per fit. Under logistic as a whole the interior reads
0.900 shipped, 0.901 told and 0.908 the move pooled, and 0.854 (0.006), 0.891 (0.004) and 0.896 (0.004) per fit.

Beside these, the same shape appears without resolving: probabilities near 0 or 1 under logistic in the
sampler loop (0.816 shipped, 0.896 told, 0.903 the move), moderate control (0.862, 0.897, 0.906) and many
predictors (0.829, 0.874, 0.866).

The closest a common-use cell comes is the interior probability in the sampler loop under logistic, tier 1:
0.892 shipped, 0.903 told, 0.902 the move, a difference of -0.010 (0.005) pooled and -0.019 (0.008) per fit.
Its long-run shortfall is 0.013 (0.005) over the 8 of 12 long runs that passed, and it reads 0.035 (0.013) over all
12, next section. Under the rule it stays unresolved.

## Multiple treatments: the largest miss, and not k's

Scored against the designed truth, 90 percent intervals on the treatment contrasts cover far below nominal under
both samplers. Ten earlier replications of each design and link, with standard errors:

| design | link | shipped | the move | move minus shipped |
|---|---|---|---|---|
| first arrangement, moderate signal | probit | 0.627 (0.035) | 0.635 (0.029) | +0.008 (0.012) |
| first arrangement, near separation | probit | 0.721 (0.017) | 0.718 (0.017) | -0.003 (0.005) |
| second arrangement, moderate signal | probit | 0.609 (0.027) | 0.614 (0.023) | +0.006 (0.017) |
| first arrangement, moderate signal | logistic | 0.587 (0.032) | 0.611 (0.034) | +0.023 (0.007) |
| first arrangement, near separation | logistic | 0.697 (0.018) | 0.713 (0.018) | +0.016 (0.009) |
| second arrangement, moderate signal | logistic | 0.566 (0.040) | 0.582 (0.038) | +0.016 (0.008) |

This is the study's largest number and rests on the fewest replications. A move that mixes k far faster barely
changes it, so it points to the trees and the model, not k, and it goes with the tree-mixing work. No told fit exists
here to say which, since these designs have no calibration. The near-separation rows carry single-precision
probabilities. Each of these designs failed its first long run, so the rule has no reading for them, and the null
design has no such figure.

## Every long run read

This reading was made after the results were known. It decides no verdict and replaces no rule, and the counts
above do not use it. It exists because the rule reads only the long runs that pass, and those are the easy ones. Under
logistic, the interior shortfall over passing runs against the shortfall over runs that failed with little
disagreement, scored against each half, reads 0.013 against 0.079 and 0.065 in the sampler loop, 0.018 against
0.061 and 0.064 in small samples, 0.043 against 0.076 and 0.071 near separation, and -0.003 against 0.062 and 0.062
with many predictors.

Here every long run with a default fit to score counts, passing or not: 250 of 292, the rest being first long
runs where nothing else ran and the weighted effect's second pipeline, which enters through the move's column. Each
is scored against the shipped half and the move half separately. The two-half gap says how far those differ, and the
worse half is the one read, with the rule's own upper bound. The halves are close enough for that: over row quantities
the median gap is 0.0018 on passing runs, 0.0045 on runs that failed with little disagreement and 0.10 on those
that disagree in substance. Where a truth exists they also track it: in failed runs under logistic near 0 or 1 the
shipped fit covers the truth at 0.809 and holds 0.804 and 0.802 of the two halves, and in stan4bart's group intercepts
it covers at 0.852 and holds 0.857 and 0.865.

Tails stay under the existing rule. Failed runs cover tail truth at only 0.880 under logistic and 0.863 under probit
pooled over rows, against 0.913 and 0.900 for passing runs, though 0.906 and 0.912 when each run counts once, so they
are a weaker reference there. Nothing below reads them.

Tier 1, the 24 readings unresolved for failed runs (shortfalls; negative means the interval holds more than 0.90):

| regime | link | estimand | long runs passed | passing runs only | all runs, against the move half | against the shipped half | worse half (se) | upper bound | two-half gap | the move's own, all runs | reads |
|---|---|---|---|---|---|---|---|---|---|---|---|
| sampler loop | logistic | interior | 8 of 12 | 0.013 | 0.035 | 0.030 | 0.035 (0.013) | 0.073 | 0.006 | 0.005 | at or past the bar |
| stan4bart, regions | probit | interior | 9 of 22 | 0.008 | 0.010 | 0.009 | 0.010 (0.003) | 0.018 | 0.003 | 0.007 | clears |
| stan4bart, regions | probit | group intercepts | 9 of 24 | 0.010 | 0.031 | 0.025 | 0.031 (0.009) | 0.054 | 0.011 | 0.022 | at or past the bar |
| stan4bart, regions | probit | group sd | 9 of 24 | 0.004 | 0.029 | 0.024 | 0.029 (0.008) | 0.051 | 0.011 | 0.024 | at or past the bar |
| propensity and weights | probit | average | 4 of 6 | 0.002 | -0.002 | -0.001 | -0.001 (0.003) | 0.008 | 0.004 | -0.001 | clears |
| propensity and weights | probit | treated | 4 of 6 | 0.002 | -0.003 | -0.002 | -0.002 (0.003) | 0.009 | 0.003 | 0.000 | clears |
| propensity and weights | probit | weights effective size | 4 of 6 | 0.003 | 0.011 | 0.011 | 0.011 (0.010) | 0.046 | 0.005 | 0.012 | open |
| propensity and weights | probit | inverse-weighted effect | 4 of 6 | 0.005 | -0.002 | -0.004 | -0.002 (0.010) | 0.034 | 0.004 | -0.007 | open |
| propensity and weights | probit | inverse weights | 4 of 6 | -0.031 | -0.040 | -0.039 | -0.039 (0.008) | -0.009 | 0.000 | -0.041 | clears |
| propensity and weights | probit | propensity | 4 of 6 | 0.003 | 0.025 | 0.029 | 0.029 (0.026) | 0.119 | 0.004 | 0.015 | at or past the bar |
| propensity and weights | probit | subgroups | 4 of 6 | 0.000 | -0.004 | -0.002 | -0.002 (0.002) | 0.006 | 0.003 | -0.003 | clears |
| propensity and weights | probit | weighted effect | 4 of 6 | 0.004 | -0.002 | -0.001 | -0.001 (0.003) | 0.010 | 0.004 | 0.001 | clears |
| weighted effect, base rate 0.30 | probit | average | 1 of 6 | -0.004 | -0.001 | -0.004 | -0.001 (0.003) | 0.009 | 0.003 | -0.003 | clears |
| weighted effect, base rate 0.30 | probit | treated | 1 of 6 | -0.005 | 0.003 | 0.001 | 0.003 (0.002) | 0.011 | 0.003 | 0.001 | clears |
| weighted effect, base rate 0.30 | probit | largest weight | 1 of 6 | -0.100 | -0.100 | -0.100 | -0.100 (0.000) | -0.100 | 0.000 | -0.100 | pinned at its bound; not a reading |
| weighted effect, base rate 0.30 | probit | subgroups | 1 of 6 | -0.004 | 0.002 | -0.001 | 0.002 (0.003) | 0.013 | 0.004 | 0.000 | clears |
| weighted effect, base rate 0.30 | probit | weighted effect | 1 of 6 | -0.004 | 0.003 | 0.000 | 0.003 (0.003) | 0.012 | 0.003 | -0.003 | clears |
| weighted effect, base rate 0.05 | probit | average | 0 of 6 |  | -0.006 | -0.002 | -0.002 (0.008) | 0.025 | 0.005 | 0.006 | open |
| weighted effect, base rate 0.05 | probit | treated | 0 of 6 |  | -0.008 | -0.006 | -0.006 (0.004) | 0.010 | 0.002 | -0.001 | clears |
| weighted effect, base rate 0.05 | probit | largest weight | 0 of 6 |  | -0.100 | -0.100 | -0.100 (0.000) | -0.100 | 0.000 | -0.100 | pinned at its bound; not a reading |
| weighted effect, base rate 0.05 | probit | subgroups | 0 of 6 |  | 0.003 | 0.005 | 0.005 (0.008) | 0.031 | 0.003 | 0.012 | open |
| weighted effect, base rate 0.05 | probit | weighted effect | 0 of 6 |  | -0.008 | -0.005 | -0.005 (0.005) | 0.014 | 0.003 | 0.003 | clears |
| multiple treatments, null | probit | unit contrasts | 3 of 4 | -0.005 | -0.004 | 0.002 | 0.002 (0.008) | 0.039 | 0.007 | -0.003 | open |
| multiple treatments, null | probit | average contrasts | 3 of 4 | -0.004 | -0.004 | 0.001 | 0.001 (0.006) | 0.032 | 0.006 | -0.002 | open |

- **Clear (12).** Stan4bart's interior with census regions; propensity and weights' average, treated, subgroup and
  weighted effects, and its inverse weights; the weighted effect's average, treated, subgroup and weighted figures
  at base rate 0.30, and treated and weighted at 0.05. The inverse weights hold 3.9 points more than 0.90 under the
  shipped sampler and 4.1 under the move. That is not an error: they are clipped at a bound, and a central interval
  on a clipped quantity holds more than 0.90 even when it is exact, about 0.943 here (the check on the told fit).
- **At or past the bar (4).** The sampler loop's interior under logistic is the one that looks like k's: the move
  sits on nominal (0.005) where the shipped sampler misses, and two of its runs fail only on row Rhat with the shipped
  sampler short by 8 to 14 points and the move within 2.3. Stan4bart's gap on group intercepts and group standard
  deviation with census regions is shared by both samplers, 0.022 and 0.024 for the move over all runs and 0.034 on
  the failed runs, against 0.043 and 0.044 for the shipped sampler on the same runs, so it is not k's; it is a real
  question for stan4bart's group scale. The propensity has the same shape (0.041 for the move on its failed runs)
  and a bound of 0.119, so it says little.
- **Still unresolved (8).** The unit and average contrasts of the null multiple-treatments design; the weighted
  effect's average and subgroups at base rate 0.05; propensity and weights' inverse-weighted effect and weights'
  effective size; and the largest weight at both base rates, which reads -0.100 because it is pinned at its bound and
  is not a reading until that is examined.

Tiers 2 and 3:

| regime | link | estimand | long runs passed | passing runs only | all runs, against the move half | against the shipped half | worse half (se) | upper bound | two-half gap | the move's own, all runs | reads |
|---|---|---|---|---|---|---|---|---|---|---|---|
| propensity and weights | probit | individual | 4 of 6 | -0.006 | 0.001 | -0.003 | 0.001 (0.007) | 0.027 | 0.005 | 0.029 | open |
| weighted effect, base rate 0.30 | probit | individual | 1 of 6 | 0.005 | 0.003 | 0.005 | 0.005 (0.012) | 0.047 | 0.003 | 0.027 | open |
| weighted effect, base rate 0.30 | probit | log odds ratio | 1 of 6 | 0.000 | -0.001 | -0.004 | -0.001 (0.003) | 0.009 | 0.003 | -0.003 | clears |
| weighted effect, base rate 0.30 | probit | log risk ratio | 1 of 6 | 0.000 | -0.002 | -0.005 | -0.002 (0.003) | 0.010 | 0.003 | -0.005 | clears |
| weighted effect, base rate 0.05 | probit | individual | 0 of 6 |  | 0.016 | 0.017 | 0.017 (0.007) | 0.041 | 0.003 | 0.023 | open |
| weighted effect, base rate 0.05 | probit | log odds ratio | 0 of 6 |  | -0.011 | -0.012 | -0.011 (0.005) | 0.008 | 0.003 | 0.003 | clears |
| weighted effect, base rate 0.05 | probit | log risk ratio | 0 of 6 |  | -0.012 | -0.013 | -0.012 (0.006) | 0.008 | 0.002 | 0.002 | clears |
| rare, 500 rows | logistic | interior | 4 of 8 | 0.007 | 0.008 | 0.013 | 0.013 (0.007) | 0.035 | 0.008 | 0.005 | open |
| rare, 500 rows | logistic | log risk ratio | 4 of 12 | 0.000 | -0.001 | 0.005 | 0.005 (0.007) | 0.023 | 0.008 | 0.000 | open |
| rare, 500 rows | probit | interior | 7 of 11 | 0.009 | 0.014 | 0.014 | 0.014 (0.007) | 0.034 | 0.003 | 0.013 | open |
| rare, 500 rows | probit | log risk ratio | 7 of 12 | 0.012 | 0.004 | 0.004 | 0.004 (0.009) | 0.029 | 0.005 | 0.003 | open |
| rare, 10,000 rows | logistic | interior | 0 of 5 |  | -0.001 | 0.003 | 0.003 (0.003) | 0.016 | 0.004 | 0.008 | clears |
| rare, 10,000 rows | logistic | log risk ratio | 0 of 6 |  | -0.017 | -0.009 | -0.009 (0.009) | 0.025 | 0.011 | -0.007 | open |
| rare, 10,000 rows | probit | interior | 0 of 6 |  | 0.057 | 0.060 | 0.060 (0.021) | 0.135 | 0.004 | 0.051 | at or past the bar |
| rare, 10,000 rows | probit | log risk ratio | 0 of 6 |  | 0.030 | 0.034 | 0.034 (0.024) | 0.118 | 0.009 | 0.034 | at or past the bar |
| many predictors | logistic | interior | 3 of 5 | -0.003 | 0.023 | 0.022 | 0.023 (0.017) | 0.089 | 0.004 | 0.015 | at or past the bar |
| many predictors | probit | interior | 3 of 5 | -0.004 | 0.008 | 0.009 | 0.009 (0.010) | 0.050 | 0.004 | 0.005 | open |
| effects, moderate | logistic | log risk ratio | 5 of 6 | 0.002 | 0.002 | 0.002 | 0.002 (0.004) | 0.015 | 0.005 | -0.003 | clears |
| effects, rare | logistic | average | 4 of 6 | 0.000 | -0.005 | -0.004 | -0.004 (0.006) | 0.019 | 0.005 | 0.002 | clears |
| effects, rare | logistic | treated | 4 of 6 | -0.001 | -0.001 | -0.002 | -0.001 (0.004) | 0.013 | 0.007 | 0.005 | clears |
| effects, rare | logistic | individual | 4 of 6 | 0.001 | 0.075 | 0.036 | 0.075 (0.082) | 0.365 | 0.052 | 0.001 | at or past the bar |
| effects, rare | logistic | log odds ratio | 4 of 6 | -0.002 | -0.006 | -0.004 | -0.004 (0.005) | 0.013 | 0.005 | -0.003 | clears |
| effects, rare | logistic | log risk ratio | 4 of 6 | -0.001 | -0.003 | -0.003 | -0.003 (0.005) | 0.015 | 0.005 | -0.001 | clears |
| effects, rare | logistic | subgroups | 4 of 6 | -0.001 | 0.006 | 0.001 | 0.006 (0.007) | 0.032 | 0.009 | -0.002 | open |

| regime | link | estimand | long runs passed | passing runs only | all runs, against the move half | against the shipped half | worse half (se) | upper bound | two-half gap | the move's own, all runs | reads |
|---|---|---|---|---|---|---|---|---|---|---|---|
| near separation | probit | interior | 3 of 5 | 0.002 | 0.006 | 0.002 | 0.006 (0.005) | 0.025 | 0.005 | 0.029 | open |

Of tier 2's 24, 10 clear, 10 stay open and 4 are at or past the bar. Among those four, many predictors' interior
under logistic has the shape of k's (0.023 shipped, 0.015 the move), as does individual effects at rare prevalence
under logistic (0.075 shipped, 0.001 the move), though one run carries that and its standard error is 0.082. Rare
outcomes at 10,000 rows under probit is shared by both samplers. Readings that already clear hold up: of the 11 with a
failed run among them, none reaches the bar over all runs, and only one bound passes it (individual effects at rare
prevalence under probit, 0.024).

What would make this reading wrong:

- **Both halves wrong together.** A failed run is not a verified posterior. Two halves can agree and both be off,
  most of all where the same trees hold both. Failed references that are too wide overstate a shortfall and too
  narrow understate it. On the interior they cover the truth at 0.924 and 0.920, against 0.910 and 0.905 for passing
  runs, within noise but on the wide side.
- **Few runs, and one or two carrying them.** The sampler loop's 0.035 rests on 12 runs and largely on two. If those
  two have a poor reference, it falls toward the 0.013 of the passing runs.
- **"Not k's" relies on the move half.** A gap shared by both samplers is the trees' or the group scale's only if
  the setup is sound. For census regions the long told fit supports that, on 25 replications: group intercepts climb
  from 0.876 to 0.907 (0.025) and the group standard deviation reads 0.880 (0.067).
- **Chosen after seeing the results.** The worse-half form, the bound and the choice to leave tails out were all
  made with the numbers in view, so this is a description and not a test.

## Under probit

No harm from k shows in the interior or in effects: the told fit is near nominal pooled over rows and the move
differs from shipped by at most about 0.009. Effects stay on nominal per fit too (0.905 shipped, 0.906 told, 0.902 the
move). Interior probabilities do not. Pooled over rows they read 0.891 shipped, 0.894 told and 0.890 the move; per
fit 0.879 (0.004), 0.884 (0.004) and 0.879 (0.005). Pooled weighting gives little weight to fits with few interior
rows, and any weighting by the data keeps an exact sampler on nominal, so the 2 points are a real cost. All three
arms share it, so it is not k's. On the 233 fits with a long told fit, ten times the length moves the told fit per fit
from 0.876 to 0.895 (0.005), a climb of 1.8 points (0.5), which points to the trees' cost at the default length and
not the model's. It is largest in rare outcomes, near separation and census regions (0.867, 0.850 and
0.873 per fit shipped) and within a point of nominal in moderate control and the sampler loop, 0.895 and 0.891.

Near 0 or 1 the picture is less clear. Fixing k also turns off probit's rescaling step, so under
probit "told" is not "the trees alone", and its near-0-or-1 coverage sits far below logistic's: 0.741 against
0.883 near separation, 0.819 against 0.891 for rare outcomes at 500 rows, 0.794 against 0.840 at 10,000, 0.815
against 0.874 with many predictors. The difference from the told fit therefore cannot show k's harm there. The
move is the fairer contrast: +0.034 (0.014) near separation, above the bar though not past the rule, and at
most +0.013 elsewhere. All three samplers climb with length, which favors the trees' cost, and the long told fit now shows it
regime by regime on the same replications: 0.870 to 0.891 with many predictors, 0.892 to 0.900 for rare outcomes at 500
rows, 0.864 to 0.893 in small samples and 0.794 to 0.849 (0.042) at 10,000 rows. Those replications are not the worst ones: over all
replications the first two read 0.815 and 0.819 at the default. Near separation goes from 0.741 to 0.854, and to 0.885 at
forty times the length, as the check on the told fit reads.

## Why: it is about reaching small k

This section groups results by the true k, which the rule does not score under calibration, and the threshold
of 1 was chosen after seeing the data. So it explains the problems and decides nothing. Small k means a strong
signal. Replications below a true k of 1: logistic interior 309 of 765, near 0 or 1 349 of 781; probit 230 of
750, 279 of 746. Pooled over rows:

| link | true k | interior: shipped minus told | interior: move minus told | near 0 or 1: shipped minus told | near 0 or 1: move minus told |
|---|---|---|---|---|---|
| logistic | below 1 | -0.063 (0.008) | -0.019 (0.004) | -0.124 (0.022) | +0.015 (0.007) |
| logistic | above 1 | +0.007 (0.002) | +0.011 (0.002) | +0.009 (0.005) | +0.014 (0.003) |
| probit | below 1 | -0.021 (0.005) | -0.023 (0.006) | +0.031 (0.011) | +0.061 (0.016) |
| probit | above 1 | -0.002 (0.002) | -0.000 (0.002) | +0.013 (0.004) | +0.019 (0.005) |

The positive differences near 0 or 1 under probit below 1 are the told fit's own shortfall at a small true k, read under
the first check above. Some loss below 1 is expected of any sampler that does not know k: calibration holds on average over k, not at
each k. That is the -0.02 in the interior that the move and probit share. The excess under shipped logistic,
about 4 points in the interior and 14 near 0 or 1, is the slow k. It survives every way of cutting it.
Weighting replications equally makes it larger (interior -0.110 (0.010), move -0.011 (0.007)). Within regime
it is -0.136 near separation, -0.167 for rare outcomes at 500 rows and -0.094 in small samples below 1, and
about +0.01 to +0.04 above 1, where shipped is no better than the move, so it is not an artifact of grouping.

It also reaches beyond the named regimes, rarely. By replication below a true k of 1, interior: moderate
control -0.043 (0.004, 5 replications), sampler loop -0.127 (0.046, 12), and -0.089 (0.023, 27) with many
predictors added. Individual effects do not show it within noise: -0.018 (0.023, 21) at moderate prevalence
and -0.013 (0.032, 32) at rare prevalence.

The long runs show the same thing directly. Of the 36 logistic long runs where the move's k goes below 0.05,
the shipped sampler's gets there in 2; of the 36 such probit long runs, in 30. Across those logistic runs the
median lowest k is 0.138 shipped and 0.0015 for the move. Probit has a rescaling step that moves k and
logistic has none, which is the likely reason.

A longer run helps, at a poor price. Four times the draws moves small samples' interior from 0.872 to 0.880
(the told fit 0.895), which brings that row inside the bar, and tail coverage gains 4 to 5 points (0.764 to
0.806 in small samples, 0.691 to 0.741 near separation). The move does better for 5 to 15 percent more time a
fit.

A k frozen in place costs little only in easy fits. Over the long runs that passed, chains each held at a k
drawn from the posterior lose under 1 point in the interior in these regimes (0.7 near separation, 0.9 in
small samples), but near 0 or 1 they lose more: 2.4 points (0.4) in small samples under logistic and 5.2 (0.9)
under probit. On runs that failed acceptance, scored against the move half, interior: 2.5 (0.6, 6 runs) near separation
under logistic, 3.3 (0.8, 2) under probit, 2.6 (0.4, 4) for rare outcomes under probit. Near 0 or 1: 7.8 (1.5,
10) and 16.9 (4.6, 3) near separation under logistic and probit, 4.3 (0.6, 15) and 16.3 (3.0, 4) in small
samples. So the prediction written down before the study, that a frozen k stays inside the bar, held for
interior probabilities and effects in the regimes it could see, with small samples under logistic at half the
bar (0.009), and only over the long runs that passed.

How much of the harm is the start, so that moving k only during burn-in would fix it, is not settled. Over the
long runs that passed, near separation's shortfall is mostly a start-up cost the move removes (0.033 of
0.043), while small samples' interior shortfall persists in a shipped chain begun at stationarity (0.014,
against 0.002 for the move). But only 2 of 12 and 9 of 24 long runs passed there, and those are the easy ones.
Pricing should test burn-in-only directly.

## Long runs that failed acceptance

141 of 292 failed acceptance at 20,000 sweeps, and this is most of why cells stay unresolved: 48 of the 50 unresolved
tier 1 and 2 readings outside the real data. The rule stops at that length and they stay failed.

Failing is not the same as not converging. The limit is Rhat at most 1.01 on every one-value quantity and 1.02 at the
99th percentile of rows, over all 16 or 20 chains, and with that many chains it is mostly a demand for precision: a
perfectly stationary chain passes 1.01 about a quarter of the time at 1,200 effective draws and two-thirds at 1,600
over 16 chains, a third at 1,600 over 20, and takes about 2,000 to 2,400 to pass reliably, with the row limit asking
about 1,200 over 16 chains, where the check on effective draws itself asks 400. Relaxing it would not rescue most of
these runs: with both limits at 1.04 and the other checks kept, 29 of the 141 pass, 13 of which already passed every
check at their last block, and with Rhat dropped altogether, 43. Most fail on more than one count: row Rhat in 117,
row effective draws in 70, agreement between the halves in 69, one-value Rhat in 35.

They fall into three groups.

- **23 disagree in substance**: k disagrees, or 5 percent or more of quantities do. The shipped sampler does
  not reach what the move reaches, even at 20,000 sweeps. 21 are logistic and 2 probit. Among simulated truths
  they are almost all logistic with a true k below 1: of the 16 simulated runs past the cut used here (next
  item), 15 are logistic and 13 of those have a true k below 1, and a logistic long run with a true k below
  0.5 passed 1 time in 16. Of the 22 real prediction datasets, 8 disagree on 20 percent or more of quantities
  under logistic and 2 under probit (mushroom 93 percent, banknote 20 percent, both with k about 0.05 to
  0.07). No reading exists for these fits, and they look like the same harm the calibration shows.
- **10 are just past the cut**: 1 to 3 percent of quantities disagree and k agrees, 6 of them under probit.
- **108 failed with little disagreement**: k agrees, and fewer than 1 percent of quantities disagree badly. "Near
  miss" would oversell them. Agreement failures there run well above chance, with 4.4 percent of quantities past 3
  standard errors against 0.7 in accepted runs, and more length is not known to pass them. 13 of the 108 passed every
  check at their last block, and 6 fail only an agreement test that a sound run fails by chance 2 to 22 percent of
  the time. The other 89 fail Rhat or effective draws. They are most of the failures in the common-use cells (15 of 24
  in stan4bart with census regions, 21 of 24 in the weighted effect).

What holds those 108 back differs by link. Medians, with the runs that passed beside them:

| long runs | number | row effective draws, shipped half | move half | k's effective draws, shipped half | move half | k's Rhat, shipped half | move half |
|---|---|---|---|---|---|---|---|
| probit, failed with little disagreement | 69 | 1,287 | 1,105 | 203 | 437 | 1.04 | 1.02 |
| logistic, failed with little disagreement | 39 | 249 | 1,171 | 38 | 1,755 | 1.14 | 1.00 |
| probit, passed | 93 | 2,986 | 3,119 | 543 | 1,019 | 1.01 | 1.01 |
| logistic, passed | 58 | 3,009 | 3,716 | 162 | 1,659 | 1.05 | 1.00 |

Under probit they are near convergence in the trees and not k: the halves have similar row draws. Under logistic it
is the shipped half's k: its rows have a fifth of the move half's draws, and its k 38 effective draws against 1,755.
Across all failures the pattern is the same, k's effective draws 30 shipped and 1,446 the move under logistic (64 runs),
177 and 433 under probit (77). Row figures here are tail effective draws, a move half's being the joint figure less
the shipped half's, and k's are bulk.

So "k agrees" is a weak test where the shipped half has few effective draws, since those widen the test. It passes k's
quantiles 10 to 35 percentile points apart in 15 logistic runs and 6 probit runs, for example wdbc under logistic, with
median k 0.32 shipped and 0.19 the move at 3.46 standard errors. Accepted runs never differ by more than about 8. Among
the 108 it is mostly true: nine in ten are within 10 points, and the worst is 18.

Running these runs further was considered and dropped before any run, because no decision now hangs on them: pricing
no longer weighs tiers, and 42 of the 48 in tier 1 cells are probit, which the move as priced does not touch.

## Real data

Each dataset reads alone. It has one long run and no known truth, so the most it can say is that its default interval
agrees with that run: 47 readings do. Its own error cannot be estimated from one run, and it measures consistency, not
calibration. Of the 22 prediction datasets 12 failed acceptance on their long run and have no reading; the three effect datasets all
passed.

Pooled figures sit beside those verdicts and cover only the datasets that passed, so they leave out the fits in question,
and those 10 of 22 mostly lack strong signal. Counting each dataset as one replication:

| datasets that passed | link | interior shortfall (se) | near 0 or 1 (se) | what it would read over those alone |
|---|---|---|---|---|
| 10 of 22 prediction | probit | 0.003 (0.001) | 0.007 (0.001) | bounded |
| 10 of 22 prediction | logistic | 0.007 (0.002) | 0.012 (0.003) | bounded |

10 of 22 is under the 80 percent rule, so neither row is a verdict. Over the three effect datasets, individual effects
under probit (-0.001, se 0.001) would read bounded, subgroups under probit and the average and treated effects under
logistic measured, and the rest unresolved.

## Could a rule find the fits that need it?

This part is exploratory. Across every calibration replication, the difference from the told fit was modeled
additively on rows, predictors, outcome rate, events, the separation score (a cross-validated tree's AUC, from
the data alone) and the shipped run's own k split-Rhat.

- **Under probit, nothing to find**: 2 percent of the interior's variation explained.
- **Under logistic the separation score is the one steady signal.** Rows, predictors and events add little,
  and Rhat adds nothing once the score is in (interior 11.6 percent explained with it and without it). The
  long-run shortfalls that passed rise with the score under logistic, not under probit.
- **Rhat will not serve.** 18 of 2,160 fits are under 1.05 (medians 1.15 to 1.65 by regime), under both links,
  harmed or not. It also needs the run to have happened.
- **A threshold of 0.8 on the score separates these regimes**:

| logistic fits | interior: fits | interior: shipped minus told | near 0 or 1: fits | near 0 or 1: shipped minus told | the move minus told, interior and near 0 or 1 |
|---|---|---|---|---|---|
| score at or above 0.8 | 443 of 765 | -0.031 (0.004) | 507 of 781 | -0.145 (0.015) | -0.001, +0.008 |
| score below 0.8 | 322 of 765 | +0.004 (0.002) | 274 of 781 | -0.006 (0.006) | +0.009, +0.016 |

Those shares reflect the regimes chosen, not use. Three limits. The contrast is nearly one of regimes: near
separation and small samples hold 438 of the 443 interior fits at or above 0.8, and outside them 5 of 327
reach it, so the model does not predict a regime held out of it (correlation 0.03). The score misses rare
outcomes, whose tail harm (-0.081) sits below 0.8. And it is an ordinary range for real prediction problems:
12 of the 22 real prediction datasets score at or above 0.8, and 10 of those 12 failed acceptance on their long run, against
2 of the 10 below. The three real effect datasets score about 0.55. So a trigger on the data is a plausible
candidate for pricing, tested on its own terms. The 0.8 was picked after seeing these data and should not be
carried over.

## What this leads to

A problem shown opens the pricing of the move, a separate plan. It needs no estimate of how common strong signal is,
since the move touches only logistic fits at a modest cost. Pricing is a sequence, with each step's rule fixed before
any fit:

1. **The move's small-k tail is right**, tested on the simulated small-sample, near-separation and rare-outcome regimes
   with k drawn from the prior. If it fails, pricing stops.
2. **It wins at equal time on coverage**, calibration being the target and over-coverage counting as error, on those
   regimes and on fits with no strong signal, which also gives the real benchmark of its cost. If so, the move is on
   for every logistic fit, with no trigger.
3. **Only if that cost is too high**, a check on the data, its form fixed on the simulated regimes alone and tested on
   real logistic datasets not yet looked at, the 22 having been seen.
4. **An opt-in** if both fail.

It has not been sized: a rough guess of tens of core-hours, mostly the second step, with the move's code priced
separately at 1,400 to 2,900 lines and a slice that changes the posterior. What must be written down before it starts
is a backlog item; the independent critique's list there is a proposal and not a ruling.

Stored fits already point the way on the second step, read after the results and deciding nothing. Per fit, the
shipped sampler at equal time against the move at default covers the interior at 0.779 and 0.904 near separation, 0.839
and 0.914 for rare outcomes at 500 rows, 0.854 and 0.887 in small samples and 0.883 and 0.903 in the sampler loop, and
near 0 or 1 at 0.712 and 0.904 near separation and 0.803 and 0.898 in small samples.

For the cost, at the default length the move gives k about 7 to 8 effective draws becoming 50 in small samples
and 7 becoming 18 near separation and for rare outcomes, for typically 5 to 15 percent more time a fit, seen
in passing on a busy machine and not a benchmark. Effective draws per second of row estimands under both
kernels are not here and cannot be recovered from what was stored; pricing must measure them afresh, with a
few fits in small samples, near separation and rare outcomes under logistic, about a core-hour or a few.
Scalars and k survive.

The first step matters because the tail looks genuine and is barely checked. Under calibration
the mean posterior share under a cutoff should match the share of true k under it, and for the logistic move
over 90 long runs it is 0.84 percent under 0.05 (se 0.35) against 1 of 90 true values (5 of 1,040 over all
replications), where the shipped sampler has 0.005 percent. That is consistent, with little power on so few
true values, and 34 of the 36 logistic long runs with that tail failed acceptance, so no long-run agreement
checks it. The tail matters: among those 36 runs the move spends on average 7.7 percent of its draws under
0.05, and on some real data the bulk of them (tictactoe 82 percent, mushroom 54, banknote 42). Latent values
reach the hundreds there (from -1,406 to 435 in one fit with a 1.2 percent outcome rate). That fit's draws and
every stored interval are finite, but a log risk or odds ratio must be computed on the log scale, as a plain
mean of probabilities underflows. Any priced move inherits that tail for fits that did not reach it before.

What the study found about the trees, the gaps both samplers share, the multiple-treatments miss and the long runs that
fail for reasons other than k, goes to the tree-mixing work whatever pricing finds.

## Tier 1 cells unresolved

Shortfall is the long-run check over the passing runs, with its upper bound; "told" is calibration's
difference from the fit told the true k; "all runs" is the reading above, the worse half and its upper bound. Multiple
treatments and the weighted effect have no calibration: their truth is designed.

| regime | link | estimand | shortfall (bound) | told (se) | all runs (bound) | long runs passed | why |
|---|---|---|---|---|---|---|---|
| moderate control | logistic | interior | 0.015 (0.027) | +0.002 (0.004) | | 6 of 6 | too few: at its estimate 13 runs in all for bounded (7 more, about 20 core-minutes); 38 for measured (32 more, about 1.5 core-hours) |
| sampler loop | logistic | interior | 0.013 (0.027) | -0.010 (0.005) | 0.035 (0.073) | 8 of 12 | failed and suspected |
| stan4bart, regions | probit | interior | 0.008 (0.016) | -0.010 (0.004) | 0.010 (0.018) | 9 of 24 | failed, clear |
| stan4bart, regions | probit | group intercepts | 0.010 (0.025) | +0.009 (0.014) | 0.031 (0.054) | 9 of 24 | failed and suspected; shared by both samplers |
| stan4bart, regions | probit | group sd | 0.004 (0.014) | +0.030 (0.030) | 0.029 (0.051) | 9 of 24 | failed and suspected; shared by both samplers |
| propensity and weights | probit | average | 0.002 (0.012) | -0.020 (0.030) | -0.001 (0.008) | 4 of 6 | failed, clear |
| propensity and weights | probit | treated | 0.002 (0.012) | -0.030 (0.030) | -0.002 (0.009) | 4 of 6 | failed, clear |
| propensity and weights | probit | subgroups | 0.000 (0.005) | -0.014 (0.011) | -0.002 (0.006) | 4 of 6 | failed, clear |
| propensity and weights | probit | inverse-weighted effect | 0.005 (0.028) | -0.015 (0.021) | -0.002 (0.034) | 4 of 6 | failed, still open |
| propensity and weights | probit | weighted effect | 0.004 (0.013) | -0.010 (0.030) | -0.001 (0.010) | 4 of 6 | failed, clear |
| propensity and weights | probit | propensity | 0.003 (0.008) | +0.014 (0.011) | 0.029 (0.119) | 4 of 6 | failed and suspected; shared by both samplers |
| propensity and weights | probit | inverse weights | -0.031 (-0.006) | +0.001 (0.001) | -0.039 (-0.009) | 4 of 6 | failed, clear; holds more than 0.90 because the weights are clipped |
| propensity and weights | probit | weights effective size | 0.003 (0.013) | +0.030 (0.030) | 0.011 (0.046) | 4 of 6 | failed, still open |
| multiple treatments, null | probit | unit and average contrasts | -0.005, -0.004 (0.037) | | 0.002, 0.001 (0.039, 0.032) | 3 of 4 | failed, still open |
| multiple treatments, null | both | unit and overall bands | | | | | a band |
| multiple treatments, other three designs | both | all | | | | first run failed in each | not run |
| weighted effect, base rate 0.30 | probit | average, treated, subgroups, weighted | about -0.005 on the one that passed | | at most 0.003 (0.013) | 1 of 6 | failed, clear |
| weighted effect, base rate 0.05 | probit | treated, weighted | none passed | | -0.006, -0.005 (0.010, 0.014) | 0 of 6 | failed, clear |
| weighted effect, base rate 0.05 | probit | average, subgroups | none passed | | -0.002, 0.005 (0.025, 0.031) | 0 of 6 | failed, still open |
| weighted effect, both base rates | probit | largest weight | -0.100 on the one that passed | | -0.100 | 1 of 6 and 0 of 6 | failed; pinned at its bound and to be examined before it is ever read |
| real data: pima, kyphosis, transfusion | logistic | interior | 0.010, 0.015, 0.008 | | | one run each, all passed | too few: 14, 5 and 4 seeds at these spreads, minutes each, under the single-dataset reading |
| real data: 12 prediction datasets | both | all | | | | first run failed in each | not run |
| your colleagues' multiple-treatments case | | | | | | | waiting on the case; room for an amendment, about 5 core-hours |

## The map

Each entry is the reading and its tier: **P** problem shown, B bounded, M measured, A agrees with its one long run, u
unresolved. Blank means the estimand does not apply. A failed reading that clears on all runs is still u here.


### Logistic

Probabilities and groups:

| |interior | near 0 or 1 | band | log risk ratio | group intercepts | group sd |
|---|---|---|---|---|---|---|
| near separation | **P**3 | **P**3 | **P**3 |  |  |  |
| rare, 500 rows | u2 | **P**3 | u3 | u2 |  |  |
| rare, 10,000 rows | u2 | u3 | u3 | u2 |  |  |
| small n, strong signal | **P**2 | **P**3 | **P**3 |  |  |  |
| many predictors | u2 | u3 | u3 |  |  |  |
| moderate control | u1 | u3 | u3 |  |  |  |
| sampler loop | u1 | u3 | u3 |  |  |  |

Effects:

| |individual | band | average | treated | subgroups | log risk ratio | log odds ratio |
|---|---|---|---|---|---|---|---|
| effects, moderate | B2 | u3 | B1 | B1 | B1 | u2 | B2 |
| effects, rare | u2 | u3 | u2 | u2 | u2 | u2 | u2 |
| real: lindner | A2 | A3 | A2 | A2 | A2 |  |  |
| real: actg175 | A2 | u3 | A1 | A1 | A1 |  |  |
| real: lalonde | A2 | u3 | A1 | A1 | A1 |  |  |

The propensity and weighting estimands did not run under logistic: those regimes are probit only.

### Probit

Probabilities and groups:

| |interior | near 0 or 1 | band | log risk ratio | group intercepts | group sd |
|---|---|---|---|---|---|---|
| near separation | u3 | u3 | u3 |  |  |  |
| rare, 500 rows | u2 | u3 | u3 | u2 |  |  |
| rare, 10,000 rows | u2 | u3 | u3 | u2 |  |  |
| small n, strong signal | u2 | u3 | u3 |  |  |  |
| many predictors | u2 | u3 | u3 |  |  |  |
| moderate control | B1 | u3 | u3 |  |  |  |
| stan4bart, twins | B1 | M3 |  |  | B1 | M1 |
| stan4bart, regions | u1 | u3 |  |  | u1 | u1 |
| sampler loop | B1 | u3 | u3 |  |  |  |

Effects:

| |individual | band | average | treated | subgroups | log risk ratio | log odds ratio |
|---|---|---|---|---|---|---|---|
| effects, moderate | B2 | u3 | M1 | B1 | B1 | B2 | B2 |
| effects, rare | B2 | u3 | B2 | B2 | B2 | B2 | B2 |
| propensity and weights | u2 | u3 | u1 | u1 | u1 |  |  |
| weighted effect, base rate 0.3 | u2 | u3 | u1 | u1 | u1 | u2 | u2 |
| weighted effect, base rate 0.05 | u2 | u3 | u1 | u1 | u1 | u2 | u2 |
| real: lindner | A2 | u3 | A2 | A2 | A2 |  |  |
| real: actg175 | A2 | u3 | A1 | A1 | A1 |  |  |
| real: lalonde | A2 | u3 | A1 | A1 | A1 |  |  |

The propensity and weighting estimands:

| |weighted effect | inverse-weighted effect | propensity | inverse weights | weights effective size | largest weight |
|---|---|---|---|---|---|---|
| propensity and weights | u1 | u1 | u1 | u1 | u1 |  |
| weighted effect, base rate 0.3 | u1 |  |  |  |  | u1 |
| weighted effect, base rate 0.05 | u1 |  |  |  |  | u1 |

### Multiple treatments

The null design only; each of the other three designs failed its first long run under both links, so nothing else ran for them.

Logistic:

| |unit contrasts | average contrasts | unit band | overall band |
|---|---|---|---|---|
| null | M1 | M1 | u1 | u1 |

Probit:

| |unit contrasts | average contrasts | unit band | overall band |
|---|---|---|---|---|
| null | u1 | u1 | u1 | u1 |

### Real prediction datasets

| dataset | logistic: interior, near 0 or 1, band | probit: interior, near 0 or 1, band |
|---|---|---|
| adult | no reading | no reading |
| bank | no reading | no reading |
| banknote | no reading | no reading |
| biopsy | no reading | no reading |
| birthwt | A1, u3, u3 | A1, u3, u3 |
| cleveland | A2, u3, u3 | A2, A3, u3 |
| climate | no reading | no reading |
| creditapproval | u3, u3, u3 | u3, u3, u3 |
| german | A1, A3, u3 | A1, A3, u3 |
| haberman | A1, u3, u3 | A1, u3, u3 |
| infert | A1, u3, u3 | A1, u3, u3 |
| ionosphere | no reading | no reading |
| kyphosis | u1, u3, u3 | A1, A3, u3 |
| magic | no reading | no reading |
| mushroom | no reading | no reading |
| pbc | A1, u3, u3 | A1, A3, u3 |
| pima | u1, u3, u3 | A1, A3, u3 |
| sonar | no reading | no reading |
| spambase | no reading | no reading |
| tictactoe | no reading | no reading |
| transfusion | u1, u3, u3 | A1, A3, u3 |
| wdbc | no reading | no reading |

## Unresolved in tiers 2 and 3, with the counts each would need

"Long runs needed" assumes each of them passes and a true shortfall of 0, and "replications needed" is at 90
percent power; where the reason is failed long runs, more runs do not help unless they pass. A band has no
long-run reading. The reasons are those of the split above, with the 8 failed readings near 0 or 1 suspected only on
arithmetic.

| tier | regime | link | estimand | reason | long runs needed | replications needed |
|---|---|---|---|---|---|---|
| 2 | rare, 500 rows | logistic | interior | failed, still open | 10 | 747 |
| 2 | rare, 500 rows | logistic | log risk ratio | failed, still open | 6 | 3040 |
| 2 | rare, 500 rows | probit | interior | failed, still open | 11 | 202 |
| 2 | rare, 500 rows | probit | log risk ratio | failed, still open | 9 | 3040 |
| 2 | rare, 10,000 rows | logistic | interior | failed, clear on all runs |  | 14 |
| 2 | rare, 10,000 rows | logistic | log risk ratio | failed, still open |  | 3040 |
| 2 | rare, 10,000 rows | probit | interior | failed and suspected |  |  |
| 2 | rare, 10,000 rows | probit | log risk ratio | failed and suspected |  | 897 |
| 2 | small n, strong signal | probit | interior | too few | 37 |  |
| 2 | many predictors | logistic | interior | failed and suspected | 5 |  |
| 2 | many predictors | probit | interior | failed, still open | 8 |  |
| 2 | effects, moderate | logistic | log risk ratio | failed, clear on all runs | 3 | 3040 |
| 2 | effects, rare | logistic | average | failed, clear on all runs | 3 | 3040 |
| 2 | effects, rare | logistic | treated | failed, clear on all runs | 4 | 3040 |
| 2 | effects, rare | logistic | individual | failed and suspected | 1 | 381 |
| 2 | effects, rare | logistic | log odds ratio | failed, clear on all runs | 6 | 3040 |
| 2 | effects, rare | logistic | log risk ratio | failed, clear on all runs | 3 | 3040 |
| 2 | effects, rare | logistic | subgroups | failed, still open | 4 | 460 |
| 2 | propensity and weights | probit | individual | failed, still open | 5 | 113 |
| 2 | weighted effect, base rate 0.30 | probit | individual | failed, still open |  |  |
| 2 | weighted effect, base rate 0.30 | probit | log odds ratio | failed, clear on all runs |  |  |
| 2 | weighted effect, base rate 0.30 | probit | log risk ratio | failed, clear on all runs |  |  |
| 2 | weighted effect, base rate 0.05 | probit | individual | failed, still open |  |  |
| 2 | weighted effect, base rate 0.05 | probit | log odds ratio | failed, clear on all runs |  |  |
| 2 | weighted effect, base rate 0.05 | probit | log risk ratio | failed, clear on all runs |  |  |
| 3 | near separation | probit | band | band | 7 | 7316 |
| 3 | near separation | probit | interior | failed, still open | 3 |  |
| 3 | near separation | probit | near 0 or 1 | failed and suspected | 6 | 1470 |
| 3 | rare, 500 rows | logistic | band | band | 17 | 592 |
| 3 | rare, 500 rows | probit | band | band | 10 | 15510 |
| 3 | rare, 500 rows | probit | near 0 or 1 | failed and suspected | 19 | 892 |
| 3 | rare, 10,000 rows | logistic | band | band |  | 4988 |
| 3 | rare, 10,000 rows | logistic | near 0 or 1 | failed and suspected |  | 38 |
| 3 | rare, 10,000 rows | probit | band | band |  | 367 |
| 3 | rare, 10,000 rows | probit | near 0 or 1 | failed and suspected |  | 32 |
| 3 | small n, strong signal | probit | band | band | 11 | 903 |
| 3 | small n, strong signal | probit | near 0 or 1 | too few | 59 | 817 |
| 3 | many predictors | logistic | band | band | 12 | 4432 |
| 3 | many predictors | logistic | near 0 or 1 | failed and suspected | 7 | 399 |
| 3 | many predictors | probit | band | band | 1 | 6097 |
| 3 | many predictors | probit | near 0 or 1 | failed and suspected | 2 | 664 |
| 3 | moderate control | logistic | band | band | 11 | 2384 |
| 3 | moderate control | logistic | near 0 or 1 | too few | 3 | 162 |
| 3 | moderate control | probit | band | band | 29 | 3720 |
| 3 | moderate control | probit | near 0 or 1 | too few | 7 | 326 |
| 3 | effects, moderate | logistic | band | band | 26 | 6122 |
| 3 | effects, moderate | probit | band | band | 9 | 5604 |
| 3 | effects, rare | logistic | band | band | 36 | 6433 |
| 3 | effects, rare | probit | band | band | 23 | 547 |
| 3 | propensity and weights | probit | band | band | 45 | 377 |
| 3 | stan4bart, regions | probit | near 0 or 1 | failed and suspected | 5 | 395 |
| 3 | sampler loop | logistic | band | band | 7 | 773 |
| 3 | sampler loop | logistic | near 0 or 1 | failed and suspected | 13 | 294 |
| 3 | sampler loop | probit | band | band | 19 | 8369 |
| 3 | sampler loop | probit | near 0 or 1 | too few | 6 |  |
| 3 | weighted effect, base rate 0.30 | probit | band | band |  |  |
| 3 | weighted effect, base rate 0.05 | probit | band | band |  |  |
| 3 | real: pima | logistic | band | band | 16 |  |
| 3 | real: pima | logistic | near 0 or 1 | too few | 33 |  |
| 3 | real: pima | probit | band | band | 1 |  |
| 3 | real: infert | logistic | band | band | 24 |  |
| 3 | real: infert | logistic | near 0 or 1 | too few |  |  |
| 3 | real: infert | probit | band | band | 5 |  |
| 3 | real: infert | probit | near 0 or 1 | too few |  |  |
| 3 | real: kyphosis | logistic | band | band | 23 |  |
| 3 | real: kyphosis | logistic | near 0 or 1 | too few | 22 |  |
| 3 | real: kyphosis | probit | band | band | 2 |  |
| 3 | real: pbc | logistic | band | band | 17 |  |
| 3 | real: pbc | logistic | near 0 or 1 | too few | 4 |  |
| 3 | real: pbc | probit | band | band | 5 |  |
| 3 | real: birthwt | logistic | band | band | 87 |  |
| 3 | real: birthwt | logistic | near 0 or 1 | too few |  |  |
| 3 | real: birthwt | probit | band | band | 31 |  |
| 3 | real: birthwt | probit | near 0 or 1 | too few | 4 |  |
| 3 | real: german | logistic | band | band | 2 |  |
| 3 | real: german | probit | band | band | 1 |  |
| 3 | real: transfusion | logistic | band | band | 8 |  |
| 3 | real: transfusion | logistic | near 0 or 1 | too few | 6 |  |
| 3 | real: transfusion | probit | band | band | 6 |  |
| 3 | real: creditapproval | logistic | band | band | 4 |  |
| 3 | real: creditapproval | logistic | interior | too few | 19 |  |
| 3 | real: creditapproval | logistic | near 0 or 1 | too few | 12 |  |
| 3 | real: creditapproval | probit | band | band | 2 |  |
| 3 | real: creditapproval | probit | interior | too few | 4 |  |
| 3 | real: creditapproval | probit | near 0 or 1 | too few | 6 |  |
| 3 | real: haberman | logistic | band | band | 6 |  |
| 3 | real: haberman | logistic | near 0 or 1 | too few |  |  |
| 3 | real: haberman | probit | band | band | 6 |  |
| 3 | real: haberman | probit | near 0 or 1 | too few | 3 |  |
| 3 | real: cleveland | logistic | band | band | 12 |  |
| 3 | real: cleveland | logistic | near 0 or 1 | too few | 4 |  |
| 3 | real: cleveland | probit | band | band | 9 |  |
| 3 | real: lindner | probit | band | band | 6 |  |
| 3 | real: actg175 | logistic | band | band | 3 |  |
| 3 | real: actg175 | probit | band | band | 6 |  |
| 3 | real: lalonde | logistic | band | band | 13 |  |
| 3 | real: lalonde | probit | band | band | 29 |  |

## What was not run

- Your colleagues' multiple-treatments case, and every other fit of the 30 readings whose design or dataset failed its first
  long run.
- Any long run past 20,000 sweeps (you dropped the near-miss runs before any ran), so failed runs stay failed.
- The long told fit on the sampler loop. Multiple treatments, the weighted effect and real data have no told fit at all,
  having no drawn truth.
- The long told fit past 25 replications, except near separation under both links and small samples and rare outcomes at
  500 rows under logistic, and past ten times the length, except near separation under probit, whose forty times was not
  taken further.
- A rule's reading at equal time, effective draws per second of row quantities, moving k during burn-in alone and held-out
  scores: pricing's, or proposals for it.

## What was decided

Each decision this study raised has its ruling.

**Does an off-nominal long told fit stop its regime?** Only if it also fails at a longer run: "We can use your
recommendation." The plan said it did, reading off nominal as a fault in the setup; the evidence under "Is the told fit a
sound comparator?" reads neither case as one. The readings stand with a note, as you ruled for accepted long runs that
over-cover, and the rule is reworded. No compute, and pricing's path did not depend on it. Counts as kept: bounded 20,
measured 52, problem 7, unresolved 159. Stopped under probit only, the same. Stopped under both links near separation: 20,
52, 4, 162.

**How common is strong signal in logistic fits?** Pricing no longer needs that estimate. Shown a default scoped by a
trigger against an opt-in, weighed by that share, you asked: "How are you planning on pricing? That seems a bit vague -
on the existing simulations and real-world datasets?" Shown the sequence under "What this leads to": "Sure, let's try
it." (dec-B530)

**Run the long runs that failed narrowly further?** No. You first ruled a quarter first: "OK, then we can use your
recommendation: run a quarter first, with a stopping rule written now." An independent critique then found no decision
hanging on them, and on that recommendation: "I'm not really worried about compute time, but with that in mind we can
use whatever your recommendation would be." Nothing ran, so the stopping rule was withdrawn and not bent. The runs stay
failed, the reading of every long run sits beside them, and, compute being no concern, the long told run covers every
regime owed, rare outcomes at 10,000 rows included. (dec-B531, dec-B537)

**How to read real data?** One dataset at a time, with the pooled figures beside them labeled as covering only those
that passed: "Use your recommendation for this item." (dec-B532)

**May one-value estimands read on long runs alone?** Yes: "Read them on long runs alone." Two readings moved to
measured; calibration's figure stays beside them. (dec-B533)

**Keep the readings where accepted long runs over-cover?** Yes: "Keep the readings with a note." The note, the
departure from the written rule and the counts either way are in the second check at the top. (dec-B534)

**Your colleagues' multiple-treatments case.** It did not run: "I asked, it could be a while before we get an answer.
Let's leave room for an amendment and proceed without it." The slot stays, about 5 core-hours. The critique's proposal,
not a ruling: every multiple-treatments design with a real effect failed its first long run, so a long-run cell there
would probably read nothing, and an amendment should score it against the designed truth across arms. (dec-B529)
