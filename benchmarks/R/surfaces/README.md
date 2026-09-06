# Published surfaces

Test problems taken from the literature, each run against the SHIPPED
default sampler and reported beside the number its source paper published
for it. Every cell has an exact truth, a pre-registered primary statistic
chosen before the run, and a generating process transcribed from the primary
source with the citation in a comment above it.

Nothing here ships with the package. Everything runs against the *installed*
dbarts, so build first:

    R CMD INSTALL .

Each script takes an output directory as its first argument and writes one
`.rds` there; with no argument it writes to a scratch directory, never into
the working tree. Append `quick` for a smoke run at reduced replicate counts
and chain lengths.

    Rscript benchmarks/R/surfaces/P2-confounded-step.R /path/to/output
    Rscript benchmarks/R/surfaces/P6-diagonal-shelf.R /path/to/output
    Rscript benchmarks/R/surfaces/P5-checkerboard.R  /path/to/output
    Rscript benchmarks/R/surfaces/C1-he-hahn.R       /path/to/output quick

`surfaces-common.R` holds the generating processes, the matched-seed idiom
and the readouts; the cell scripts source it and are the only things meant to
be run.

## The cells

### P2, the confounded step function

Pratola (2016) section 2.3, taking the problem from Wu, Tjelmeland and West
(2007). Three hundred rows, three columns, one tree. The covariates are
drawn in blocks that confound x1 with x3, so root-on-x1 and root-on-x3 are
two exactly equiprobable representations of the same fitted function.

- Primary statistic: between-chain standard deviation of the fraction of
  draws whose root splits on x1, against the exact 0.5 symmetry the design
  supplies.
- Published reference: "the acceptance rate of tree moves (after the initial
  few steps of the sampler) was 0" (arXiv 1312.1895 section 2.3).
- Two arms on matched seeds, both reachable through `proposal.probs`: the
  shipped mixture, and birth/death only.
- Runs with its own null control, a design with two exactly duplicated
  predictor columns on a four-value grid, where the pooled fraction on the
  first of the pair is 1/2 by construction and switching must occur.

### P6, the diagonal shelf with targeted selection

Hahn, Murray and Carvalho (2020) Example 1. Two uniform covariates, 250
rows, a homogeneous treatment effect of -1, and a propensity that tracks a
prognostic shelf running along x1 = x2. The estimand is the outer one.

- Primary statistic: ATE bias and 95% interval coverage over 200
  replications.
- Published reference: standard BART bias 0.27, coverage 65%, RMSE 0.31;
  the propensity-augmented BCF prior 0.14, 95% and 0.21 (arXiv 1706.09523v4
  Table 1).
- The paper does not print its prognostic function. Three reconstructions
  are run: one calibrated to the paper's own Figure 4, and two widths of the
  symmetric near-step family its Figure 3 caption describes. No arm of this
  cell can carry a reproduction verdict on its own.

### P5, the checkerboard on an autocorrelated design

Zhu, Zeng and Kosorok (2015) scenario 3, at n = 1600 and p = 40. A pure
two-way interaction, f = 2 x5 x10 + 2 x15 x20, on predictors with covariance
0.9^|j-k|, so every true column sits between two decoys correlated 0.9 with
it.

- Primary statistic: between-chain standard deviation of time-averaged
  inclusion on {x5, x10, x15, x20} and on their immediate neighbours, read
  against a mixing null computed from the same draws.
- Published reference: none for dbarts. The oracle is the exact inclusion
  truth and its near-decoys.

### C1, the He and Hahn factorial

He and Hahn (2023) section 4.1, correlated-factor predictor arm, at
n = 10000, p = 30 and kappa = 1, on two of the four mean functions:
Trig+poly for the interaction and Single index for the rotated ridge.

- Primary statistic: 95% pointwise coverage of the true mean function.
- Published reference (arXiv 2002.03375v4 Table 4, BART column, kappa = 1):
  Trig+Poly coverage 0.74, interval length 2.89, RMSE 1.27; Single Index
  0.73, 4.62 and 2.08.
- The paper's section 5 does not say which predictor arm or how many trees
  its Table 4 used, so both are run as diagnostic arms beside the
  pre-registered correlated-design, default-tree-count one.

## Conventions

Matched seeds, the idiom the move-set grid already uses: the data seed is
indexed by cell, design and replicate; the sampler seed is the replicate
alone. Both are shared across arms, so every contrast within a replicate is
paired.

Structural readouts - inclusion proportions, root split variables - are not
label invariant and are read BETWEEN chains, never within one.

Effective sample size is `posterior::ess_basic`.
