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
`.rds` there (P1 one per rung); with no argument it writes to a scratch directory, never into
the working tree. Append `quick` for a smoke run at reduced replicate counts
and chain lengths.

    Rscript benchmarks/R/surfaces/P1-friedman.R        /path/to/output house
    Rscript benchmarks/R/surfaces/P2-confounded-step.R /path/to/output
    Rscript benchmarks/R/surfaces/P2-null-at-scale.R  /path/to/output
    Rscript benchmarks/R/surfaces/P6-diagonal-shelf.R /path/to/output
    Rscript benchmarks/R/surfaces/P5-checkerboard.R  /path/to/output
    Rscript benchmarks/R/surfaces/C1-he-hahn.R       /path/to/output quick

`P1-friedman-census.R` is the one script here that does not run against an
ordinary build: it reads per-move acceptance out of the instrumented library
`benchmarks/R/move-census.R` documents.

`surfaces-common.R` holds the generating processes, the matched-seed idiom
and the readouts; the cell scripts source it and are the only things meant to
be run.

## The cells

### P1, the low-noise Friedman emulator

Pratola (2016) section 2.2, the battery's known-positive control and its
absolute gate. Two rungs. Pratola's own: Friedman's function observed with
noise on ten uniform columns, n = 5000, 200 trees, 5000 burn-in and 5000
kept, at sigma^2 = 1 and sigma^2 = 0.1. The house rung: the cheaper cell the
move-set A/B measured, n = 2000 and sigma = 0.25, 1000 burn-in and 2000 kept.

- Primary statistic: 90% pointwise coverage of the true f, at the training
  points and on 1000 held-out rows drawn with each replicate.
- Published reference (arXiv 1312.1895 section 2.2, birth/death only):
  acceptance around 18% and coverage 81% at sigma^2 = 1; acceptance around
  4% and coverage 53.8% at sigma^2 = 0.1. Both are read as in-sample - the
  section reports them without qualification and its figure plots the
  intervals against the fitted settings, where the paper's section 6 names
  its out-of-sample coverage as such. The response is eta plus noise, not
  the deterministic simulator output.
- The house rung is the gate the rest of the battery depends on: its default
  arm must sit near the 0.71 that all three mixtures of the move-set A/B
  read, or no verdict from any other cell is valid.
- Three arms on matched seeds, all through `proposal.probs`: the shipped
  mixture, birth/death only (Pratola's own arm), and the historical mixture
  carrying swap at 0.1.
- Pratola prints the mean function with `10 sin(2 pi x1 x2)` where Friedman
  (1991) has `10 sin(pi x1 x2)`, and does not say which his figures were
  taken on. His rung is run on the printed one, with the sigma^2 = 0.1 cell
  run a second time on Friedman's own frequency under the birth/death arm so
  that the ambiguity is measured rather than assumed.
- `P1-friedman-census.R` carries the structural readout, one seed of each
  rung and arm under the instrumented build, reporting acceptance per move.

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
- Three arms on matched seeds, all reachable through `proposal.probs`: the
  shipped mixture (birth/death and change, swap at zero), birth/death only,
  and the same mixture with swap set to 0.1.
- Runs with its own null control, a design with two exactly duplicated
  predictor columns on a four-value grid, where the pooled fraction on the
  first of the pair is 1/2 by construction and switching must occur.
- `P2-null-at-scale.R` re-runs that null control at production tree counts
  (50 and 200, the shipped and swap-carrying arms), reading a tree-level root
  share instead of a chain-level one: does the single-tree stuck-on-x3
  representation survive an ensemble, or does the ensemble wash it out.

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
- A fifth arm, `correlated75grow`, adds the opt-in XBART grow-from-root
  warm start (`n.grow.sweeps = 5`) to the pre-registered arm's settings,
  paired on the same seeds. A sixth, `independent75grow`, pairs the same
  warm start with the independent design instead, against `independent75`.
- Three chain-configuration arms on the independent design at 75 trees,
  the same seeds: `independent75pool4` (four chains of 500 + 500, pooled -
  `bart2`'s shipped default), `independent75pool4long` (four chains of
  1000 + 2500, pooled) and `independent75long` (one chain of 1000 + 25000,
  with coverage also read on its first 2500 kept draws). Pooled arms report
  per-chain minimum ESS summed and a between-chain ratio: the median over
  the ESS points of the between-chain sd of each chain's posterior mean
  over the pooled posterior sd.
- Four move-set arms on that same four-chain configuration and the same
  seeds, varying only `proposal.probs`: `independent75pool4bd` (birth/death
  only), `independent75pool4swap` (the former default, carrying swap at 0.1),
  `independent75pool4perturbB` (perturb at 0.16 taken out of change alone)
  and `independent75pool4perturbMixed` (the same 0.16 taken out of change and
  birth/death both, which is a different mixture and not the arm above).
  `independent75pool4` is their paired control and is re-run beside them, so
  the earlier move-set grid's one-chain reading is retaken where the shipped
  chain default reads it. They report, beside the pooled readouts, the
  per-seed paired difference against the control in summed minimum ESS,
  coverage and RMSE, and for the control arm alone the per-chain 95%
  intervals at the 25 ESS points and how far they overlap. Both perturb
  levels are pilots, not a confirmatory run of the perturb design.

## Conventions

Matched seeds, the idiom the move-set grid already uses: the data seed is
indexed by cell, design and replicate; the sampler seed is the replicate
alone. Both are shared across arms, so every contrast within a replicate is
paired.

Structural readouts - inclusion proportions, root split variables - are not
label invariant and are read BETWEEN chains, never within one.

Effective sample size is `posterior::ess_basic`.
