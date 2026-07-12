# GPUs for BART: a design survey

Status: survey memo, 2026-07-07. No commitment to any mechanism. This document
weighs the design space and ends in a recommendation; any mechanism it endorses
gets its own plan and its own gate class. See docs/plans/gpu-bart.md for the
originating item, and docs/design/parallel-bart-frontier.md for the companion
note on not-yet-existing constructions and the CRAN conditional-backend axis.

Update (2026-07-10): the CPU cut-scan prototype that sections 2d and 4 gate the
GPU cut-scan experiment behind has landed, role (a) warm-start producer only
(src/bartcore/scan.hpp, src/bartcore/grow.hpp; docs/design/grow-from-root.md).
That prerequisite is satisfied. The GPU falsifying experiment itself remains
open work, sequenced behind frontier items 4-5 (the informed-kernel and
CG-leaf prototypes; parallel-bart-frontier.md section 5).

## 1. Why classic BART is GPU-hostile, and where that framing breaks

The default dbarts sampler is sequential backfitting (Chipman, George and
McCulloch 2010). One sweep visits the m trees in order; tree t is drawn against
the residual net of every other tree's current fit, so tree t+1 cannot start
until tree t's fit is retired (src/bartcore/chain.hpp, the treeY running-residual
loop). Three properties make this look wrong for a GPU:

- Serial in the tree dimension. The m trees do not update concurrently; the
  conditional dependence is the whole point of the Gibbs construction. At
  defaults that is ~200 strictly ordered steps per sweep.
- Per-tree work is O(n) but shallow-parallel. Each tree's cost is a partition
  of the node's rows (misc_partitionIndices, a two-pointer swap) and a per-leaf
  reduction for the sufficient statistic (misc/moments.c: sum w, sum wz,
  sum wz^2). The reductions are over leaf-sized index sets that shrink fast as
  the tree descends, so the exploitable width collapses near the leaves.
- Branchy structure logic. The Metropolis-Hastings moves (moves.hpp:
  birth/death/change/swap) walk a pointer tree, draw rules, and score integrated
  likelihoods over branch nodes. This is control-flow-heavy, divergent code, the
  opposite of a dense SIMT kernel.

The honest correction: this framing overstates the obstacle, and the literature
already refuted the strong version of it. Petrillo (2024, the JAX package bartz)
runs the exact single-chain BART MCMC on a GPU and reports up to 200x over a
single CPU core, GPU-competitive from n ~ 20,000 and dominant at n >= 1e6. The
recipe: represent every tree as a fixed-depth dense heap (max depth 6, above which
BART trees are rare at default priors), turning branchy routing into branchless
arithmetic; vectorize traversal and leaf sampling across all n observations and
all m trees at once; keep the one irreducibly sequential-in-trees piece, the
per-leaf residual reduction, on device. The serial tree dimension is real, but the
per-tree work over the observation axis is where the parallelism lives, and at
large n that axis is wide enough to fill a device without changing the sampler.
The obstacle is not the algorithm; it is that dbarts's engine uses the opposite
representation (sparse, arena-allocated pointer trees, unbounded depth, divergent
moves), so the bartz result is an existence proof, not a portable one.

## 2. The design space

Each direction below is rated by its GPU ceiling (how much of a device it can
actually keep busy), its fit to the current dbarts engine, and its cost.

### 2a. Prediction and test-fit offload

Mechanism: route the numTest x numTrees traversal, or saved-tree replay in
predict(), on the device. Rows are independent and writes are disjoint.

Ceiling: modest, and shrinking. The CPU version already landed
(test-fit-parallel, 2026-07-07: 1.48x wall at numTest = n = 1e5). A GPU wins only
when numTest x numTrees is very large and the host already holds the flattened
forest. Fit: clean, no sampler change. Cost: low, but so is the payoff, and it
does nothing for the sampling bottleneck that motivates this work.

### 2b. Many-chain parallelism (device-per-chain, chains-as-batch)

Mechanism: chains are already the coarse parallel axis (chain.hpp fans chains
across std::thread workers, results thread-count invariant). Put one chain per
device, or batch several chains as a leading tensor dimension SIMT-style. The
distributed-data variants of this idea are worked out in the literature: Pratola
et al. (2014) split observations across cores under an SPMD master-slave scheme,
and Entezari, Craiu and Rosenthal (2018) and Luo and Pratola (2023) shard the data
and combine sub-posteriors. All target cluster CPUs, not one device, and address
throughput rather than the single-chain latency this work is about.

Ceiling: bounded by the number of chains you actually want, a handful for one
dataset. Chains-as-batch also fights SIMT divergence: independent chains grow
differently shaped trees, so a batched kernel either pads to a common shape (the
bartz dense-heap trick, which then is really direction 2g) or diverges. Critically,
this axis does nothing for the confirmed-common flagship case, a single chain at
n >= 1e5 (within-chain-threading.md). Fit: good. Cost: moderate. Value for the
motivating workload: low.

### 2c. Blocked-jacobi noise splitting (b = m)

Mechanism: the construction in blocked-jacobi-trees.md. Augment with per-tree
pseudo-responses y_k ~ N(g_k, sigma^2 / b) constrained to the batch residual;
marginally the model is unchanged, so a Gibbs alternation of a cheap Gaussian
bridge draw and b independent single-tree updates is exact. At b = m every tree
updates concurrently, which is the one classic-BART formulation that manufactures
m-wide concurrent structure work per step.

Ceiling on paper: high concurrency. Ceiling in practice: undercut from two sides.
First, the price the plan already names: a structure move within a batch sees
precision b / sigma^2, so a birth changing fit by delta pays exp(-b delta^2 /
2 sigma^2) and grows conservative as b rises. ESS per second is the only honest
metric, and the conservative-birth tax debits it. A cautionary precedent:
Entezari, Craiu and Rosenthal's (2018) Likelihood Inflating Sampling Algorithm
raises each shard's likelihood to a power (the same precision-scaling move in
disguise) and documents that inflated-precision tree updates distort tree topology
enough to need a dedicated correction, the same failure mode b > 1 courts. Second,
and decisive for the GPU question: bartz shows the observation axis alone already
saturates a device at b = 1, so blocked-jacobi's entire reason to exist, to
synthesize concurrent per-step work, is obviated on a GPU by the parallelism
already present over n. It may still earn its keep as a CPU within-chain-threading
candidate (the m-wide axis maps onto a thread pool, m/b barriers per sweep instead
of ~3m), but as a GPU strategy it solves a scarcity of parallel work a GPU does not
have. Fit: a new exact kernel gated posterior-changing. Cost: high. GPU value: low
once bartz's lesson is admitted.

### 2d. Grow-from-root / XBART-style root-down sampling

Mechanism: He, Yalov and Hahn (2019); He and Hahn (2023). Regrow each tree from
the root each iteration by a stochastic recursive split, choosing among all cuts
of all variables at a node from one histogram pass (the planned misc_scanColumn
cut-scan kernel, kernel-vocabulary.md addition 5). This is exactly the histogram
primitive GPU gradient boosting is built on.

Ceiling: the highest among directions that fit the dbarts architecture. A cut-scan
evaluates every candidate split of every variable at a node in one sweep over the
node's rows, accumulating (count, sum, ssq) per bin. That is embarrassingly
parallel over rows and over bins, sparsity-aware, and reducible by histogram
subtraction, precisely the pattern GPU GBDT implementations exploit (Wen et al.
2019; Zhang, Si and Hsieh 2017; XGBoost's gpu_hist and LightGBM's histograms).
Fit: the engine was provisioned for it. MoveStrategy is the designed seam
(core-generalization.md), the cut-scan kernel is on the vocabulary's planned list,
and tree import/export exists. This is the realistic doorway per the backlog
(grow-from-root.md). Cost: it is a posterior-changing sampler, not MH-exact, so it
inherits the approximate-sampler validity story and its own gate; and the cut-scan
kernel must first prove out on CPU per that item before any GPU port is worth
discussing. GPU value: high, but bundled with a sampler change and sequenced behind
a CPU prerequisite.

### 2e. Particle and SMC tree samplers

Mechanism: Lakshminarayanan, Roy and Teh's top-down particle filter for decision
trees (2013) and Particle Gibbs for BART (2015) propose whole trees to fit the
residual rather than making local MH edits, running many particles that are
naturally parallel and mix better on deep trees.

Ceiling: moderate-to-high; particles are independent parallel work, and top-down
particle growth is histogram-shaped like grow-from-root. Fit: poor today, no PG
MoveStrategy exists and it is a larger sampler change than grow-from-root. Cost:
high. It shares grow-from-root's GPU-shaped inner kernel (the cut-scan) without
grow-from-root's provisioning, so as a GPU on-ramp it is dominated by 2d.

### 2f. Mixed CPU/GPU split (device owns O(n) kernels, host owns tree logic)

Mechanism: keep the branchy move logic on the host and offload the O(n) kernels,
partition and suffstat accumulation and the residual reduction, to the device per
node operation.

Ceiling: architecturally this matches the kernel-vocabulary boundary exactly (the
generic layer only picks which kernel on which index set). But the granularity is
wrong. A node's suffstat is O(n_leaf), which halves at every level, so a per-node
offload pays a host-device transfer of the index set and residual for work that
shrinks toward nothing near the leaves. bartz's central lesson refutes it: never
round-trip, keep the entire state resident for the whole chain. A design that
PCIe-hops per MH proposal loses to CPU on latency alone. Fit: tempting on paper.
Cost: moderate. GPU value: low; the transfer granularity defeats it.

### 2g. Dense-array single-chain rewrite (the bartz recipe)

Mechanism: adopt Petrillo's design wholesale, a from-scratch engine in an array
framework (JAX or CUDA) with fixed-depth dense-heap trees and everything resident
on device, running the exact single-chain MCMC.

Ceiling: the highest of any option and already proven at the exact regime that
motivates this (single chain, n >= 1e5, up to 200x, memory ~ n(p + m) bytes).
Cost: also the highest, and structurally incompatible with dbarts. It is a
parallel reimplementation, not a backend to the existing sparse pointer-tree
engine; it caps tree depth (a mild prior approximation); and it does not obviously
carry dbarts's distinguishing contract, data mutable between iterations for
embedding in an outer Gibbs sampler (the whole reason dbartsSampler exists). Do not
rebuild it: bartz already ships it as a Python package (stochtree ships XBART).
A user needing GPU BART today is better served by being pointed there than by a
CUDA backend grafted into this engine.

## 3. Constraints every proposal inherits

From the plan, and non-negotiable:

- CRAN ships no CUDA toolchain. Any GPU path is an optional backend probed at
  configure time and absent from the default build, or a separate package. It can
  never be on the critical path of a plain `R CMD INSTALL`.
- Thread-count-invariant reductions. The engine guarantees bitwise-identical
  results across thread counts by fixed-block reduction order (within-chain-
  threading.md). A GPU reduction has its own nondeterministic accumulation order;
  a GPU backend cannot claim the CPU path's invariance and must define its own
  determinism contract.
- Exact-equivalence testing culture. dbarts gates changes bitwise against saved
  baselines (benchmarks/R/equivalence.R) and against exact posteriors on
  enumerable problems (logistic-reference.R, categorical-exact.R). A GPU backend
  needs its own gate class: MH-exact offloads (2a, 2f) gate to tolerance against
  the CPU draws; posterior-changing samplers (2c at b > 1, 2d, 2e) gate against the
  exact posterior on small problems, never against the CPU sampler's draws.
- The motivating workload is real. Single-chain n >= 1e5 is the confirmed-common
  case (VD 2026-07-06). Any direction that does not serve it (2a, 2b) is off-target
  regardless of its GPU ceiling.

## 4. Recommendation

No direction earns a dbarts GPU prototype yet. The one direction with both a high
GPU ceiling and a real seam in this engine is grow-from-root (2d), but its own item
already sequences a CPU cut-scan kernel prototype ahead of any sampler work, so a
GPU cut-scan cannot be evaluated until that CPU kernel exists and is confirmed the
bottleneck: the GPU work here is blocked on a prerequisite that is not GPU work.
The highest-ceiling option overall (2g) is proven but is a from-scratch rewrite
incompatible with the engine and the embeddable-mutation contract, and its value is
already captured by bartz, so it does not earn a backend. Blocked-jacobi (2c) is
obviated as a GPU strategy by the observation-axis parallelism bartz demonstrates
at b = 1 and survives only as a CPU threading candidate. The rest (2a, 2b, 2e, 2f)
are off-target for the flagship workload or dominated on ceiling, fit or cost.

Ranked, if a slot ever opens: 2d (grow-from-root cut-scan) first, behind its CPU
prerequisite; 2g (adopt or point users to bartz) as the pragmatic answer for GPU
BART today; everything else below the line.

Cheapest falsifying experiment. The single measurement that would flip this "no"
to "prototype worth it": port the one cut-scan kernel (misc_scanColumn, all cuts of
all variables over a node's index set into per-bin count/sum/ssq) to CUDA or JAX,
and benchmark one node-expansion GPU vs the CPU scalar/SIMD kernel at n in
{1e5, 1e6}, INCLUDING the host-device transfer of the residual and index vectors.
If the GPU kernel does not beat CPU by a wide margin with transfer counted, then no
per-node-offload direction (2d as an offload, 2f) clears the bar, and the only
viable GPU path is a fully resident rewrite (2g), which is out of scope as a
backend, closing the question. This experiment should run only after
grow-from-root's CPU cut-scan prototype lands, since the kernel's on-CPU cost is the
prerequisite datum and the kernel would be built there first anyway.

What would change the ranking:

- A confirmed dbarts workload at n >= 1e6 where a single chain is DRAM-bound and
  profiling puts the cut-scan or suffstat reductions at a dominant share of
  runtime. That raises 2d and revives 2f.
- grow-from-root landing on CPU with its histogram kernel measuring as a large,
  device-shaped fraction of the sweep. That makes the GPU cut-scan experiment the
  natural next step rather than a speculative one.
- bartz or stochtree exposing a mutable-data C boundary that supports the
  embedded-Gibbs contract. Adoption (2g via a bridge) would then strictly beat any
  from-scratch backend, and dbarts's role becomes the R-side driver.
- A published exact GPU BART sampler that keeps unbounded-depth sparse trees and
  beats bartz's dense-heap design. None exists today; the field's one GPU BART
  result is the dense-heap one, the strongest evidence that dbarts's current
  representation is the wrong starting point for a device.

## References

- Chipman, George, McCulloch (2010). BART: Bayesian Additive Regression Trees.
  Annals of Applied Statistics 4(1):266-298.
- He, Yalov, Hahn (2019). XBART: Accelerated Bayesian Additive Regression Trees.
  AISTATS 2019 (PMLR); arXiv:1810.02215.
- He, Hahn (2023). Stochastic Tree Ensembles for Regularized Nonlinear Regression.
  Journal of the American Statistical Association 118(541):551-570.
- Lakshminarayanan, Roy, Teh (2013). Top-down Particle Filtering for Bayesian
  Decision Trees. ICML 2013 (PMLR v28); arXiv:1303.0561.
- Lakshminarayanan, Roy, Teh (2015). Particle Gibbs for Bayesian Additive
  Regression Trees. AISTATS 2015 (PMLR v38); arXiv:1502.04622.
- Petrillo (2024). Very Fast Bayesian Additive Regression Trees on GPU.
  arXiv:2410.23244. Python package bartz.
- Pratola et al. (2014). Parallel Bayesian Additive Regression Trees. Journal of
  Computational and Graphical Statistics 23(3):830-852; arXiv:1309.1906.
- Entezari, Craiu, Rosenthal (2018). Likelihood Inflating Sampling Algorithm.
  Canadian Journal of Statistics 46(1):147-175; arXiv:1605.02113.
- Luo, Pratola (2023). Sharded Bayesian Additive Regression Trees.
  arXiv:2306.00361.
- Wen, Shi, He, Chen, Ramamohanarao (2019). Exploiting GPUs for Efficient Gradient
  Boosting Decision Tree Training. IEEE Transactions on Parallel and Distributed
  Systems 30(12).
- Zhang, Si, Hsieh (2017). GPU-acceleration for Large-scale Tree Boosting.
  arXiv:1706.08359 (basis of ThunderGBM).
