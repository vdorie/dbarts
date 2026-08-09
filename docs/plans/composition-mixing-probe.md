# composition-mixing-probe

status: RE-REGISTERED (v2), NOT RUN. v1 was refuted as pre-registered by
  blind adversarial critique - eleven blocking findings, nine of them
  established by measurement on v1's own prototypes. All eleven adopted.
  This file replaces v1 rather than amending it; the gate architecture, the
  ridge kill, the power arithmetic and five comparability details are new.
  Adjudication record: `.claude/composition-probe/synthesis.md`.
  Critique: `.claude/composition-probe/critique.md`.
agent: opus (harness + analysis; no engine change)
rng: neutral - measurement only. No verdict here authorizes an engine
  change of any kind (see "What a GREEN licenses").
budget: harness only, under `.claude/composition-probe/harness/`; nothing
  under `R/`, `src/`, `inst/`, `benchmarks/`. ~1100-1300 lines across seven
  files. That is **M** on the record's own size scale
  (`tree-mixing-proposals.md:379-381`), not the "XS-to-S" that record's
  sec 4.1 assumed; the correction is noted rather than argued with.
window: none - no timing metric exists anywhere in this design, so the
  measurement is load-insensitive. Runs on the laptop. It must NOT share a
  machine with anything timing-gated (`benchmarks/R/bench-sampler.R`, the
  cheap-uniformity S4 falsifier), because this probe would perturb THEM.
tip: `b4b8614` (branch `bartcore`), engine `9929ede`. Every code anchor
  below re-verified at that tip.
sequencing: runs after the `multiforest-veto-rate-falsifier` measurement.
  No hurry (VD 2026-08-09).

## Goal

Answer, with the first number anyone has taken: **does absorbing a smooth
share of the signal into a parametric block improve the forest's own
sampling behaviour, and how big is the cross-block ridge?**

`docs/design/tree-mixing-proposals.md` sec 4.1 ranks parametric absorption
by composition first among all candidates, on the claim that

> signal that is multimodal in tree space is unimodal in coefficient space

- a smooth separable surface has combinatorially many renderings as a deep
tree ensemble and exactly one as a coefficient vector, so moving that share
into a parametric block removes the renderings rather than asking the chain
to navigate them. The record also states the hazard head-on: additive
competition puts a ridge between the two blocks and alternating updates
crawl along ridges, measured in this house at ~6x on the grouped surrogate
(`forest-ranef-interweaving.md`). The hazard is established in print three
times over (Hahn et al.'s regularization-induced confounding; BCF's mu/tau
aliasing; CSP-BART's shared-covariate non-identifiability). **The benefit
is measured nowhere**: a full-text audit of ten papers found no effective
sample size, autocorrelation, R-hat or Gelman-Rubin for any composed
parametric-plus-BART sampler, and stan4bart's own sec 4.5 poses the
question and answers it "More research will need to be performed to confirm
this" (2022, still open).

This probe produces that measurement. Its whole value is the number; no
engine change follows from it either way.

**Scope sentence, registered up front.** Even under the design below, the
claim on offer is *"adding an absorbable share to the DGP degrades this
arm's forest less than it degrades BART-alone's."* That is a within-arm
sensitivity, not "the composed model is sampled better". A0 and A2 target
different models; no design can turn a comparison of their sampling into a
comparison of one model sampled two ways. The within-arm dial is what makes
the question answerable at all.

## Context

- Source of truth: `docs/design/tree-mixing-proposals.md`, sec 4.1 (the
  candidate and its falsifier sketch), sec 3.1-3.5 (the five stickiness
  modes), sec 7 (why this runs), sec 9 (what the survey could not settle).
- `docs/design/forest-ranef-interweaving.md` sec 0, 2, 5, 6, 9 - the only
  in-house measurement of the hazard (56.1 vs 9.3 IACT alternating vs
  collapsed at 3 groups; ~3-9x on the real engine; ASIS made it *worse*).
- `docs/design/grow-from-root-default.md` sec 3, 4.4, 4.8, 4.9 and
  `docs/plans/grow-from-root-default-study.md` - the house design law this
  file is bound by, and the source of the one reusable statistic (C2).
- Form: `docs/plans/grow-from-root-default-study.md`,
  `docs/plans/multiforest-veto-rate-falsifier.md`.
- Feasibility log with every executed command:
  `.claude/composition-probe/feasibility.md` (gitignored). Prototypes under
  `.claude/composition-probe/proto/`.

Three house laws this design is bound by:

1. **At m = 75 no structural statistic can detect structural mode
   collapse** - the ensemble self-averages a label like "the root splits on
   x1" at rate `1/sqrt(#splits)` with no mixing required. Structure-switch
   counts are therefore reported-only at m = 75. **Depth and leaf count are
   exempt, and the exemption is this file's, argued not cited** (the record
   never mentions either statistic): the self-averaging objection is about
   a *discrete label* whose distribution over trees is what the average
   estimates, so the average converges without the chain moving. Mean
   realized depth is a per-tree *size*, not a label; the ensemble average of
   75 sizes does not converge to the posterior mean of one tree's size
   unless the sizes themselves move. A frozen forest has a frozen depth,
   and a forest that self-averages a root label can still have a depth that
   changes every sweep. So the objection does not transfer, and depth and
   leaf count are gated at m = 75.
2. **R-hat is not gateable** on this package's slow functional at feasible
   replicate counts. This design does not use R-hat. It does use a
   between-chain inefficiency factor built from the same variance
   decomposition; that statistic is gated as a *dial difference*, not as a
   level, and it carries a registered exit at FREEZE.
3. **Per-cell kill checks, never a pooled aggregate**; thresholds frozen
   against a pilot before any confirmatory contrast is looked at; a
   mandatory fresh-seed re-run of any single flagged cell.

## What changed from v1, in one place

| v1 | v2 | why |
|---|---|---|
| gate on raw cross-arm contrasts at D50 | gate on **within-arm dial slopes**, cross-arm only on the slopes | A1 is 20% shallower and 45% worse on T2 than A0 at s = 0, with nothing to absorb; the raw D50 contrast is 2.5x the transfer effect it claims to measure (critique B3, cross-target ruling) |
| dial pins total variance, `c_U` moves | dial fixes **`c_U = 1`**, `c_A` moves; a **placebo cell DP** pays for the SNR change | the DiD needs the un-absorbable target identical across the pair |
| D0 is a control whose failure voids the study | D0 is the **calibration anchor** of the DiD | measured, D0 voids the study as registered |
| D100 is a scientific control | D100 is a **harness validity check**, A3 excluded | unreachable on the `c_U = 1` dial, and A3's leaf prior is ~5x tighter there |
| T2 raw between-chain SD, gated | **E2** dimensionless inefficiency factor, gated as a log DiD | raw dispersion confounds sampling inefficiency with posterior width |
| T4 (XOR root dispersion), gated | reported, with the conditioning subsample size published | mechanically biased toward the tested conclusion |
| **KILL-2** on `R1 = IACT(d)/IACT(s)` against A1's | **no ridge kill.** Absolute IACTs reported with a pre-registered ladder | R1(A1) measured 4.9-33.5, not ~1; no cell calibrates a no-ridge null; the coordinate is degenerate at the dial anchor |
| `4 x SE <= margin`, margins set at FREEZE | alpha 0.05 one-sided, **power 0.90 at margins written now** | `4 x SE` is 50% power at the margin and the error direction is a spurious KILL |
| C = 4 chains | **C = 8** | the record's SMALL stratum ran 24 x 8; v1 contradicted itself |
| `R CMD INSTALL .`, no `--preclean` | **`--preclean`**, engine hash recorded | no header dependency tracking in `src/Makevars` |
| A3 as one C-chain fit | A3 as **8 single-chain fits**, `skip = 1L`, `cores = 1L` | reproducibility, memory, and V3's exposure identity |
| A3 in G2 | **A3 out of G2** | A3 has no sigma draw (`resid.prior = fixed(1)`) |

## The absorbable-share dial

Friedman's function decomposes into terms in disjoint coordinates, so its
variance components add exactly. Computed (MC at N = 4e6 for the
interaction term, closed form for the rest; `feasibility.md` sec 7):

| term | Var under U(0,1) |
|---|---|
| `10 sin(pi x1 x2)` - the interaction | 11.19601 |
| `20 (x3 - 0.5)^2` - smooth, univariate, nonlinear | 2.22222 (= 400/180) |
| `10 x4 + 5 x5` - exactly linear | 10.41667 (= 125/12) |
| **total** | **23.83490** (sd 4.88210) |

`V_U = 11.19601 + 2.22222 = 13.41823` (un-absorbable), `V_A = 10.41667`
(absorbable). The DGP is

```
f(x) = c_U * [10 sin(pi x1 x2) + 20 (x3 - 0.5)^2] + c_A * [10 x4 + 5 x5]
y    = f(x) + N(0, 1),   x ~ U(0,1)^p,  p = 10
```

**`c_U` is fixed at 1 at every gated cell**, and the absorbable share `s`
is dialed by `c_A` alone:

```
c_A(s) = sqrt( s/(1-s) * V_U / V_A )
```

| cell | c_U | c_A | Var(f) | share s | role |
|---|---|---|---|---|---|
| **D0** | 1.00000 | 0.00000 | 13.41823 | 0.000 | DiD anchor; calibration cell |
| **D50** | 1.00000 | 1.13497 | 26.83646 | 0.500 | the gated cell |
| **DP** | 1.41421 | 0.00000 | 26.83646 | 0.000 | placebo: SNR matched to D50, nothing absorbable |
| **D100** | 0.00000 | 1.13497 | 13.41823 | 1.000 | harness validity check only |

**Why `c_U = 1` and not pinned total variance.** The DiD's whole purpose is
that everything constitutive cancels between the two cells of a pair. That
requires **the un-absorbable target the forest must find to be
byte-identical** across the pair. v1's pinning moved `c_U` from 1.333 to
0.942, so the forest's own job changed with the dial and the DiD was only
first-order clean. v1's stated reason for pinning was that stickiness mode
3.2 (structure freezes as sigma falls) is an SNR mechanism; measured, that
mechanism does not appear to bind here at all - realized sigma sat at
0.999-1.032 across every arm and both dial levels - so the pinning was the
right instinct for the wrong reason. The SNR change is paid for by DP
instead, which measures it directly.

**The trilemma, stated.** Three desiderata - `c_U` identical, SNR
identical, share different - and two free coefficients. All three cannot
hold. So both achievable pairs run:

- **P1 = (D0, D50)**: `c_U` identical, SNR moves 13.42 -> 26.84 at
  sigma = 1. **PRIMARY.**
- **P2 = (DP, D50)**: SNR identical, un-absorbable scale moves
  sqrt(2) -> 1. **CORROBORATION.**

A transfer claim must survive both in sign. Disagreement is YELLOW with the
disagreement as the finding.

**D100 is not on this dial** - `c_A -> infinity` as `s -> 1` under
`c_U = 1` - so it is not a scientific cell. It runs as a harness check with
`c_A` byte-identical to D50's and `c_U = 0`, at reduced replication, on
A0/A1/A2 only (A3 excluded: its warmup rescale anchors a response transform
~5x tighter there, see "A3's exact discrepancies").

Note, worth stating: **the standard Friedman function already sits at a
linear-absorbable share of 0.437**, so D50 is very nearly the canonical
BART benchmark.

**The absorbable span is `Z = [x4, x5]` - exactly linear, no intercept, no
basis choice.** This is the mechanism's *weak* case: an exactly-linear
monotone term is close to the fewest-renderings case, while the record's
argument is about surfaces with many renderings. v1 chose it for tidiness
(no basis-specification freedom) and left the strong case ungated. v2 keeps
it primary - it is the only span with no basis freedom, and A1 can only be
designated on actual columns - and **promotes the strong case to
gated-secondary** as cell pair Q1.

**Covariate overlap is ON in the primary specification.** The forest sees
all p = 10 columns including x4 and x5 in every arm. Forced for A1
(designating leaf covariates does not remove them as split candidates), so
any disjoint specification would make A1 incomparable. It also happens to
be exactly the configuration stan4bart sec 4.5 defends as beneficial
parameter expansion and CSP-BART attacks as harmful non-identifiability.
The disjoint variant runs as a ridge-mechanism cell (DJ), **A2 only**,
because at DJ arm A0 is misspecified by construction and is not a valid
comparator.

**Z carries no intercept column, in any arm.** The BART component owns the
constant. Verified necessary: a prototype with a leading ones column
produced a lag-0 cross-block correlation of -0.95 that was pure
intercept-versus-forest-mean aliasing, exactly the trap CSP-BART names.
Every ridge statistic is computed on mean-centred projections.

## Arms

Four scientific arms plus three calibration/mechanism arms. Everything
already ships; nothing is built in the package.

| arm | what it is | what runs it | role |
|---|---|---|---|
| **A0** | BART alone | `dbarts()` sampler, constant leaves | reference slope |
| **A1** | inner composition | `dbarts(..., node.prior = linear(c("x4","x5")))` | scientific |
| **A2** | outer composition, conjugate normal block, R-level Gibbs | `dbarts()` + `$setOffset()` between sweeps | scientific |
| **A3** | outer composition, WALNUTS parametric block | `stan4bart(y ~ bart(x1+...+x10) + x4 + x5, ...)` | scientific |
| **A2r** | A2 with stan4bart's warmup rescale schedule | `$setOffset(offset, updateScale = TRUE)` on the decaying schedule during burn | calibration: isolates A3's rescale |
| **A2c** | A2 with componentwise beta draws | two univariate conditionals instead of the joint draw | mechanism: isolates the alternation the ridge is about |
| **A2f** | A2 at `resid.prior = fixed(1)` | reported cell DF | mechanism: removes the freeze confound by construction |

**A0.** `dbarts(y ~ ., df, control = dbartsControl(n.trees = 75L,
n.chains = 1L, n.samples = 1L, n.burn = 0L, keepTrees = FALSE,
updateState = FALSE, n.threads = 1L), seed = s)`, driven `run(0L, 1L)` for
B + K sweeps, one sampler object per chain.

**A1.** Identical, plus `node.prior = linear(c("x4", "x5"))`. Verified
surface: `linear` (`R/model.R:886`) is bound into the prior evaluation
environment by `parsePriors` (`R/model.R:90-127`) from `dbartsPriors`
(`R/model.R:1076-1085`) and resolves in call position on `dbarts()`. It is
NOT reachable from `bart()` or `bart2()`, which construct
`node.prior = normal(k)` internally (`R/bart.R:370-373`) - **which is why
no arm uses `bart2()` and the grow-from-root `fitArm()` wrapper cannot be
reused.** Designation limits, verified in source: at most 8 columns
(`model.hpp:913`, enforced `facade.hpp:531`), continuous only, no sparse
matrix, no variance forest (`facade.hpp:573-576`), no monotone constraint,
no `growFromRoot`.

**A2.** One `dbartsSampler` per chain. Per sweep:

```
r     <- s$run(0L, 1L)                       # one BART sweep
fhat  <- as.vector(r$train) - Z %*% beta     # strip the installed offset
sig   <- as.numeric(r$sigma)[1L]
Vn    <- chol2inv(chol(crossprod(Z)/sig^2 + V0inv))
beta  <- Vn %*% (crossprod(Z, y - fhat)/sig^2) + t(chol(Vn)) %*% rnorm(ncol(Z))
s$setOffset(as.vector(Z %*% beta))
```

`setOffset` is `R/dbarts.R:1004`, shipped and public. `r$train` **includes**
the installed offset, so the strip is correct; verified by a five-line
check (install `offset = 100`, run to stationarity, confirm `mean(train)`
returns to `mean(y)`; it does, 7.447 against 7.417). That check is V9.

`V0inv = diag(1e-2, 2)`. At n = 5000 each slope's posterior sd is
~`sigma / sqrt((Z'Z)_jj)` ~ 0.025 against a prior sd of 10, so the prior is
negligible and v1's registered sensitivity sweep over `diag(1)` /
`diag(1e-4)` would be a null check. It is downgraded to a one-line
assertion that the posterior sd is < 1% of the prior sd. **The sensitivity
that matters is joint versus componentwise drawing**, which is the
alternation the ridge is about; that is arm A2c.

**A3.** Eight independent single-chain fits per replicate:

```
stan4bart(y ~ bart(x1 + x2 + ... + x10) + x4 + x5, data = df,
          chains = 1L, cores = 1L, skip = 1L,
          iter = B + K, warmup = B, seed = s * 1000L + chainIndex,
          bart_args = list(n.trees = 75L, keepTrees = TRUE), verbose = -1L)
```

Eight single-chain fits rather than one eight-chain fit, for three reasons:
it matches A0-A2's one-chain-per-object discipline; it removes
`makeCluster`'s silent fallback (stan4bart's own test asserts `cores = 1`
and `cores = 2` differ at the same seed,
`stan4bart/inst/tinytest/test-05-rng.R:60`); and it caps
`extract(fit, "indiv.bart")` at 5000 x 2000 doubles (80 MB) instead of
640 MB.

`skip = 1L` is **pinned and asserted**: it wires into dbarts' `n.thin`
(`stan4bart/R/stan4bart_fit.R:474`, `:515`) and dbarts runs
`(numBurnIn + numSamples) * numThin` sweeps (`chain.hpp:917`), so V3's
exposure identity is silently false at any other value.

**Matched exposure, and its exact limits.** The exposure unit is the BART
sweep, and stan4bart runs exactly one per HMC iteration
(`stan4bart/src/init.cpp:642`, inside the loop at `:570`), with `iter`
counting warmup (verified: iter 60 / warmup 30 returns 30 kept draws). No
adaptation or thinning branch. So `iter = B + K, warmup = B` delivers
exactly B burn and K kept BART sweeps. Three exact discrepancies:

1. **A3 gets one extra BART sweep** at `init.cpp:245`, before the loop.
   Total exposure B + K + 1. Ruled **noise**: one sweep in 4001, preceded
   by `sampleTreesFromPrior` with zero-valued leaves, entirely inside burn.
2. **A3's sweep contains no sigma draw.** `resid.prior = fixed(1)`
   (`stan4bart/R/stan4bart_fit.R:553`) sets `sigmaIsFixed`; `chain.hpp:1085`
   guards the draw; sigma arrives from WALNUTS via
   `dbarts_sampler_setSigma` (`init.cpp:617`). So matched *sweeps* are not
   matched *Gibbs scans*: A0-A2 execute (structure, leaves, sigma), A3
   executes (structure, leaves) with sigma from another block under another
   conditional. **Registered consequences: A3 exits G2 entirely** (a
   chi-squared draw and a WALNUTS draw under different priors are not a
   harm comparison), and **A3 is registered as the shipped stan4bart
   sampler measured as a whole, not a controlled variant of A2.** Its DiD
   is interpretable within-arm; its cross-arm use is limited to G1.
3. **A3 rescales the BART response during warmup.** `init.cpp:634-635`,
   `update_scale_mod = 1 << (8 * iter / numIter)` with `numIter` the
   *phase* count (`:509`): at iter 4000 / warmup 2000 it fires **498 times,
   last at warmup iteration 1920, zero times after**, plus one at setup
   (`:233`). "Warmup-only" is true on firing and **false on consequence**:
   `updateScale` re-anchors the *persistent* internal response transform
   (`dbarts.h:363-367`), which sets the leaf prior scale via
   `nodeScale / sqrt(numTrees)` (`model.hpp:915`), which governs how
   readily the forest splits - i.e. depth and leaf count, the primary limb.

   **Ruled bias, and handled two ways.** (i) Under the `c_U = 1` dial the
   rescale's target is dial-invariant to first order: the offset-adjusted
   residual sd is `sqrt(V_U + 1) = 3.797` at *both* D0 and D50, because the
   absorbable block is exactly what the offset removes and `c_U` does not
   move. So the DiD differences it out to first order. (ii) The second-order
   residual - the anchor is set at warmup iteration 1920, when the
   parametric block has not fully converged - is **measured, not assumed**,
   by arm **A2r**, which drives A2 with the identical decaying schedule
   during burn (`updateScale = isWarmup && (iter %% (2^(8*iter/B)) == 0)`).
   The critic verified `setOffset(offset, updateScale = TRUE)` succeeds on
   a plain dbarts sampler at the live tip. **If |DiD(A2r, T1) -
   DiD(A2, T1)| exceeds the T1 margin, A3 exits the T1/T1b gate and the
   results say so in the title line.**

   D100 is where the first-order argument fails (`c_U = 0` makes the target
   ~5x tighter), so **D100 excludes A3**.

## Scenarios and cells

`n = 5000`, `p = 10`, `sigma = 1`, `m = 75` unless stated. **`C = 8`
chains** (the record's SMALL stratum; v1's C = 4 halved the between-chain
dof from 7 to 3 on statistics that are between-chain dispersions).
`B = 2000` burn sweeps, `K = 2000` kept sweeps, `n.thin = 1`. 200 held-out
points from the same DGP for every prediction readout. `n.burn = 0` in the
control; burn is executed by the driver loop so every arm executes an
identical number of `run(0L, 1L)` calls with identical per-sweep
bookkeeping.

| id | content | m | arms | R | role |
|---|---|---|---|---|---|
| **D0** | dial anchor, `c_A = 0` | 75 | A0 A1 A2 A3 | 24 | **primary pair P1** + calibration |
| **D50** | `c_A = 1.13497` | 75 | A0 A1 A2 A3 | 24 | **primary pair P1 and P2** |
| **DP** | placebo, `c_U = sqrt(2)`, `c_A = 0` | 75 | A0 A1 A2 A3 | 24 | **corroboration pair P2** |
| **DQ0 / DQ50** | Q1: absorbable span widened to `Z = [x3, x3sq, x4, x5]`, un-absorbable block = the interaction alone, `c_A' = 0 / 0.94119` | 75 | A0 A2 | 24 | **GATED-SECONDARY** |
| **D100** | `c_U = 0`, `c_A = 1.13497`, B = K = 500 | 75 | A0 A1 A2 | 4 | harness validity check |
| **DR0 / DR50** | A2r rescale calibration | 75 | A2 A2r | 8 | calibration for A3 |
| **DC0 / DC50** | A2c componentwise-beta mechanism | 75 | A2 A2c | 8 | reported |
| **DF0 / DF50** | fixed sigma, `resid.prior = fixed(1)` | 75 | A0 A2f | 8 | reported |
| **DJ0 / DJ50** | disjoint: forest sees `x1,x2,x3,x6..x10` only | 75 | A2 | 8 | reported ridge mechanism |
| **DG0 / DG50** | D0/D50 plus a 10-group random intercept, tau = 1 | 75 | A0 A3 | 8 | reported ranef bridge |
| **X0 / X50** | `4*XOR(x1>.5, x2>.5) + c_A''*(10 x4 + 5 x5) + eps`, `c_A'' = 0 / 0.61968` | **1** | A0 A1 A2 A3 | 8 | REPORTED |
| **N0 / N50** | `x1 == x2` exactly, `5 x1 + c_A_N*(10 x4 + 5 x5) + eps`, `c_A_N = 0 / 0.44721` | **1** | A0 A1 A2 A3 | 8 | reported estimator certificate for X |

Q1's coefficients: `V_U' = 11.19601` (interaction alone),
`V_A' = 2.22222 + 10.41667 = 12.63889`, so `c_A'(0.5) = 0.94119`; the
interaction carries coefficient 1 at both levels. **A1 is excluded from
Q1**, and the exclusion is registered with its reason: designating `x3sq`
as a leaf covariate requires adding it as a design-matrix column, which
changes the forest's split candidate set and so changes the DGP's
comparability. Excluding A1 is also what makes Q1 affordable at full R.

**X and N are REPORTED, not gated.** Four independent reasons, any one
survivable, together not: (i) adding absorbable terms voids the record's
0.3619 baseline, which this design says itself; (ii) the record's own
version of the statistic **left the gate** on floor-above-ceiling at 24 x 8
(`grow-from-root-default.md:408-410`), and the right mixing null for that
configuration is 0.0368, not 0.05 (`:406`); (iii) the conditional statistic
carries a mechanical bias toward the tested conclusion (see T4 below); (iv)
the null control N is near-certain to fire - the record's cold arm failed
"non-zero switches in every chain" at 4/64 and 10/192 - and nothing gated
should ride on that. Demoting X removes the need for a V6 contingency
entirely: **N is a reported caveat on X, and voids nothing.**

## Metrics

Indices: chain `c = 1..C`, kept sweep `t = 1..K`, tree `j = 1..m`,
column `q = 1..p`.

**Every gated quantity is a within-arm difference across a dial pair.** For
arm X, functional F and pair P = (cell_lo, cell_hi):

```
Delta_X(F, P) = F(X, cell_hi) - F(X, cell_lo)
DiD_X(F, P)   = Delta_X(F, P) - Delta_A0(F, P)
```

Gates read `DiD_X`. **No gate compares a level across arms** (one exemption,
G1, argued where it appears). Every constitutive arm difference - the leaf
model, the extra block, A3's setup sweep, A3's rescale to first order - is
present at both dial levels and cancels.

### Tree-space readouts

| id | quantity | role |
|---|---|---|
| **T1** | mean realized tree depth: max leaf depth per tree, averaged over trees, chains, kept sweeps | **PRIMARY, limb (a)** |
| **T1b** | mean leaf count per tree, same averaging | **limb (a) FALLBACK**; also reported as the full distribution |
| **E2** | inclusion-dispersion inefficiency factor, between-chain route | **PRIMARY, limb (b)** |
| **E3** | inclusion autocorrelation time, within-chain route | **limb (b) alternate** |
| **E4** | IACT of per-sweep mean tree depth | **limb (b) second alternate**; reported |
| **T3** | per-tree root-variable switch count | REPORTED at m = 75 (house law 1); reported at m = 1 |
| **T4** | at X: root-on-x1 statistics | REPORTED, see below |
| **T5** | IACT of held-out error, and of sigma (A0/A1/A2 only) | REPORTED - the record's secondary readouts, restored |

T1 and T1b come from a pre-order stack walk of `getTrees(current = TRUE)`
per sweep. Verified mechanics: rows are depth-first over
`(chain, sample, tree, n, var, value)`, `var == -1` marks a leaf, internal
`var` is 1-based; `current = TRUE` **drops the `sample` column always and
the `chain` column when `n.chains == 1`**, so the harness restores them
before any grouped operation.

**E2, definition.** With `v_{c,t,q}` the per-draw inclusion proportion
(`varcount` row normalized to sum 1),

```
vbar_{c,q} = mean_t v_{c,t,q}
B_q        = var_c(vbar_{c,q})          # C-1 denominator, C = 8
W_q        = mean_c var_t(v_{c,t,q})
E2_q       = sqrt( K * B_q / W_q )
E2         = mean_q E2_q
```

`E2` is dimensionless and equals ~1 at perfect mixing: `W_q / K` is the
variance the chain mean would have under iid draws, so the ratio is a
variance inflation and `E2_q^2` estimates the integrated autocorrelation
time by the between/within route. **Posterior width cancels**, which is
exactly the defect that removed v1's raw T2: a smaller raw between-chain SD
can simply mean the arm's inclusion proportions have a tighter posterior,
which is what absorbing x4 and x5 should produce. Gated on `log E2`, as a
DiD.

**E3, definition.** `E3 = mean_q ( mean_c IACT(v_{.,.,q}) )` with `IACT`
from `coda::effectiveSize` on the chains x draws matrix. Same functional as
E2, autocorrelation route. Gated on `log E3`, as a DiD.

**E4.** `IACT` of the per-sweep mean tree depth series, pooled over chains.
A second tree-space autocorrelation on an independent functional. Reported;
promoted into limb (b) only if E2 and E3 both exit at FREEZE.

**T4, and its registered bias.** At X, T4 is a *conditional* statistic -
the fraction of kept draws with root on x1 among draws whose root is in
{x1, x2}. This design's own rationale for X is that at m = 1 arm A0's
single tree must spend itself on the linear terms before it reaches the
interaction, which is precisely a statement that **A0's conditioning
subsample is smaller**, so its per-chain fraction is estimated from fewer
draws and its between-chain SD is mechanically larger before any mixing
difference exists. Registered treatment, since X is reported not gated:
publish `n_cond` per arm per chain beside every T4 number, and report the
**per-conditioned-draw root-switch rate** and the **root-indicator IACT**
as the primary X statistics, with the raw between-chain SD reported only
alongside `n_cond`.

**Structural acceptance rate by move type: DELETED, with reason.** The
record lists it as a PRIMARY tree-space readout (`:665`). It is **not
reachable**: no per-move acceptance tally is exposed on the R surface, in
`inst/include/dbarts/dbarts.h`, or by `run()`. Adding one is an engine
change, which GREEN's scope forbids. Recorded rather than silently dropped.

### The ridge readouts - MEASURED AND REPORTED, NOT GATED

**There is no ridge kill in this design.** Three independent reasons:

1. **No cell calibrates a no-ridge reference.** v1 nominated A1 as the
   empirical null on the theory that `model.hpp:1023-1040` integrates the
   leaf coefficients out of the structural score. Measured across four
   seeds with the reconstruction validated at `cor(recon, r$train) = 1` on
   12,000 sweeps: `R1(A1)` = 4.9 / 7.6 / 20.3 / 33.5, and `R2(A1)` = -0.52
   to -0.71 (i.e. the nominated null clears v1's own "confirms the ridge
   mechanism" bar in 4/4 seeds). The disjoint-covariate cell is not a null
   either (R1 1.6-5.2): **the two blocks compete through the shared
   residual, not through the design matrix.** Removing shared covariates
   does not remove the ridge. That is itself a finding, and it speaks
   directly to the CSP-BART / stan4bart sec 4.5 disagreement.
2. **The record's kill line is a cross-arm level comparison** of a
   coordinate that is not the same object across arms: in A1 `d_t` is a
   within-tree slope-versus-constant coordinate, in A2/A3 a between-block
   coordinate. It falls to this design's own no-cross-arm-levels rule.
3. **No within-arm DiD exists for it**, because `a_t` is identically zero
   at the dial anchor.

A fourth cell pair (`s = 0.25` vs `s = 0.5`) would restore a within-arm
ridge DiD. **Considered and declined**: it answers "does doubling the
absorbable share worsen the ridge?", not "does the ridge eat the transfer?",
and it costs a main cell.

**Consequence for the verdict, stated plainly: the probe cannot execute the
record's second kill line, and the inner-versus-outer question is NOT
adjudicated here.** A GREEN re-ranks the composition candidate as a whole;
no guidance following a GREEN may prefer linear leaves to `setOffset` or
the reverse on this evidence.

**Operational definition.** Let `U` be an orthonormal basis of the
mean-centred span of `[x4, x5]`, restricted to a fixed random subsample of
`n_r = 1000` rows (registered once per replicate, shared by every arm and
every sweep). Defining `U` from the **span**, not from the scaled mean
vector `c_A * (10 x4 + 5 x5)`, is required: the latter is 0/0 at `s = 0`.
For each basis direction `u_k`:

```
a_t^k = < parametric-block fitted vector at sweep t, u_k >
b_t^k = < forest fitted vector at sweep t (offset stripped, centred), u_k >
s_t^k = a_t^k + b_t^k     the total absorbable fit  - identified
d_t^k = a_t^k - b_t^k     the ridge coordinate
```

Projecting a two-dimensional block onto a single convenience direction
picks a direction for tidiness, not for maximal slowness. So the headline
takes the **worst case over three candidate directions**: `u_1`, `u_2`, and
the leading eigenvector of the lag-1 cross-covariance of `(a_t, b_t)`.

| id | quantity |
|---|---|
| **W1** | `IACT(d*)` at the worst-case direction, absolute, per chain, pooled |
| **W2** | `IACT(s*)` at the same direction |
| **W3** | `IACT(d*) / IACT(s*)` - reported, never a headline: it is a ratio and it deletes the effect. Measured absolute IACTs: A1 `s_t` 0.94-1.00 / `d_t` 4.9-31.6; A2 `s_t` 8.87-11.32 / `d_t` 19.1-172.7 - a large, consistent signal that the ratio reverses on one of four seeds |
| **W4** | `cor(a_t, b_t)` per direction, per chain |
| **W5** | `IACT(a_t)` - the parametric block's own mixing |
| **W6** | `cor(a_t, b_{t+k})`, k = 0..20, as a curve |

**Pre-registered directional interpretation**, on an absolute, model-free
scale (effective sample size in the kept window), per arm, at D50, over
R_ridge = 8 replicates:

- **MATERIAL**: `IACT(d*) >= 50` (ESS below 40 per 2000-sweep chain, below
  320 pooled over 8) in at least 2/3 of replicates, AND `W3 >= 3`.
- **IMMATERIAL**: `IACT(d*) < 20` in at least 2/3 of replicates.
- **INDETERMINATE** otherwise, which is a legitimate reading given the
  measured seed-to-seed spread.

Reported per arm in {A1, A2, A2c, A3} plus cell DJ, with the standing
caveat that `d_t` is not the same object across arms, so cross-arm ordering
is descriptive only. Also reported: the seed-to-seed spread of every IACT,
since `coda::effectiveSize` is an AR-spectral estimator and
`forest-ranef-interweaving.md` records single-seed excursions to 310
against a control maximum of 42 on exactly this class of coordinate.

**The A1 reconstruction.** `a_t` for A1 is the leaf-slope contribution,
`b_t` the leaf-constant contribution, both reconstructed from `getTrees`'
`beta.<col>` columns with per-observation leaf assignment. Two registered
facts: `beta.<col>` is on the **standardized-covariate and
internal-response** scale, not the raw scale (measured: for
`y = 10 x4 + 5 x5 + eps` at m = 1, `beta.x4 = 0.1659` against
`10 sd(x4) / diff(range(y)) = 0.1706`, and against a naive raw prediction
of 10), with standardization using the engine's own moments (sample sd,
n - 1; `data.hpp:62-80`, matching R's `scale()`); and the reconstruction is
gated by **V8**, `cor(reconstructed total fit, r$train) == 1` to machine
precision, which validates leaf assignment, the `x <= value` split
convention and the standardization in one number.

### Guardrails - composition must not buy mixing at accuracy's expense

| id | quantity | role |
|---|---|---|
| **G1** | RMSE of the pooled posterior mean against the true `f` on 200 held-out points, kept window | **gated harm check**, cross-arm level (exemption argued below) |
| **G2a** | mean sigma over the kept window, arm X above A0 | **gated harm check**, A0/A1/A2 only |
| **G2b** | the sigma trajectory, split at the last A3 rescale firing (warmup iteration 1920 at B = 2000), and the "fit improved, forest froze" conjunction | REPORTED |
| **G3** | recovery of the un-absorbable interaction: RMSE against `10 sin(pi x1 x2)` on a held-out grid where only x1, x2 vary and x3..x10 sit at medians, after removing the fitted grid mean | REPORTED |
| **G4** | 90% pointwise interval coverage of the true mean | REPORTED |

**G1's cross-arm exemption, argued.** The no-cross-arm-levels rule exists
because a dispersion under one posterior is not comparable to a dispersion
under another. G1 is not a posterior property: it is the error of an
estimator against a **known truth `f` that is identical across arms within
a cell**. Comparing two estimators of the same object is well posed
regardless of their posteriors. G1's DiD is also reported, for symmetry.

**G2's two directions, separated.** As a harm test, higher sigma is harm.
As the record's freeze diagnostic - the reason G2 exists - the signal is
*lower* sigma. One one-sided test cannot serve both, so G2a is the harm
test (one-sided, higher is harm) and G2b is the reported diagnostic.
Measured expectation: sigma sat at 0.999-1.032 across all arms and both
dial levels, so the freeze mechanism, which is keyed on *realized* sigma
(`tree-mixing-proposals.md:201-209`), does not appear to bind at sigma = 1.
Cell DF (`resid.prior = fixed(1)`, verified to work including alongside
`node.prior = linear(...)`) turns that expectation into a control.

**A3 exits G2 entirely** (no sigma draw; see "A3's exact discrepancies").

## Validity gates (Stage 0; failure blocks the metrics named)

- **V0 - depth walk.** `#nodes == 2 * #leaves - 1` for every tree of every
  arm in one small cell; the pre-order stack walk agrees **exactly** with
  `dbarts:::getTreeDepthAndSize` (`R/plotTree.R:39-55`, the recursive
  depth-and-size walk over the same frame, called from `R/dbarts.R:1457`)
  on 200 randomly chosen trees. v1 claimed no such function existed and
  said "verified"; it does, and it is the independent implementation.
  Failure blocks T1, T1b, T3, E4.
- **V1 - frame shapes, per arm.** `getTrees(current = TRUE)` column set
  matches the documented one for the chain count in use;
  `extract(s4bFit, "trees")` returns the same six columns; `varcount` is
  `chains x draws x p` from `run()` and `p x (chains*draws)` from
  stan4bart's `extract`, with the chain split asserted, not assumed. Any
  shape difference stops the run.
- **V2 - no intercept anywhere.** `rownames(extract(fit, "fixef"))`
  contains no intercept term (confirmed as a general rule:
  `stan4bart/R/lme4_functions.R:185-188` drops `(Intercept)`
  unconditionally, the BART term carries it), asserted per fit anyway
  because it is free. A2's `Z` has no ones column.
- **V3 - exposure identity.** Each arm executed exactly `B + K` BART
  sweeps, except A3 which executed `B + K + 1`; asserted from the driver's
  counters and from `dim(extract(fit, "trees"))`. `skip == 1L` and
  `cores == 1L` asserted on every A3 fit.
- **V4 - stan4bart still agrees with itself.** stan4bart's own tinytest
  suite passes against the freshly `--preclean`-installed dbarts tip.
  **Never run; it is the real check on "no rebuild needed".** 21 files,
  ~1-3 min at CRAN settings, 10-30 min with `NOT_CRAN=true`; run the full
  suite. `dbarts.h` has been byte-stable since `876a339` (verified at the
  live tip), but the ABI is not the whole contract. Failure triggers the
  three-arm degradation below.
- **V5 - the dial is calibrated.** On 20 fresh datasets per cell, realized
  `var(f)` within 2% of nominal (13.41823 / 26.83646 / 26.83646 / 13.41823
  for D0 / D50 / DP / D100) and realized absorbable share within 0.02 of
  nominal. Fails only on an arithmetic error.
- **V6 - reproducibility.** Two invocations at the same `BASE_SEED` produce
  identical tables for A0/A1/A2/A2r/A2c/A2f. A3's is stan4bart's to
  guarantee and is asserted: two fits at the same seed agree bitwise on
  `extract(fit, "sigma")`.
- **V7 - A3 at m = 1 runs.** Downgraded from a gate to an assertion:
  `n.trees` is forwarded untouched with only a `> 0` check
  (`R/A_class.R:326-328`) and stan4bart's own suite already fits at
  `n.trees = 1L` (`inst/tinytest/test-06-no_ranef.R:19`).
- **V8 - the A1 reconstruction identity.** `cor(reconstructed total fit,
  r$train) == 1` to machine precision on every sweep of one small cell.
  Failure voids the A1 ridge readout only.
- **V9 - offset semantics.** Install `offset = 100`, run to stationarity,
  confirm `mean(r$train)` returns to `mean(y)` rather than
  `mean(y) - 100`. Establishes that A2's `fhat <- r$train - Z %*% beta` is
  correct.
- **V10 - the D100 collapse.** A0/A1/A2 at D100: every composition arm's
  `T1b` falls below a threshold frozen at Stage 0 from A0's D100 value and
  the theoretical floor of 1. If the forest does not collapse when the
  signal lies exactly in Z, **the harness is wrong** and the run stops.
  This is the only survivor of v1's "controls" section, and it was always
  a harness check rather than a scientific control.

**Nothing voids the estimator family on a scientific reading.** v1
registered D0 and N as voiding controls; both were measured or predicted to
fire (A1 differs from A0 at s = 0 by 20% on depth and 45% on T2; the
record's cold arm failed N's switch clause at 4/64 and 10/192). D0 is now
the DiD anchor, whose job is precisely to measure those offsets, and N is a
reported caveat on a reported cell.

### Degradation contingency, pre-registered

**If V4 fails, or A3 cannot be built for a cell, the probe runs with three
arms in that cell and records the loss in the results title line**, not a
footnote. A three-arm probe still answers KILL-1; what it loses is the
HMC-block arm, i.e. the specific adjudication of stan4bart sec 4.5.

**If A2r shows the rescale moves T1's DiD by more than the T1 margin**, A3
exits the T1/T1b gate, also in the title line.

## Pre-registered decision lines

**Sidedness: every gated test is ONE-SIDED at alpha 0.05**, with the
direction stated per metric. A limb clears only if the one-sided test
rejects **and** the point estimate lies beyond the margin - the record's
actual practice (`grow-from-root-default.md:696`, `:708`;
`-study.md:157`). v1 used `4 x SE` as the decision threshold itself, which
is a ~4-sigma critical value with no precedent in the record and which
inverts the direction of conservatism.

**Every margin is written HERE, NOW, before Stage 0 measures anything.**
The two-bar rule (a SCIENTIFIC ceiling and a STATISTICAL floor) degenerates
to one self-referential bar if the ceiling is chosen after the standard
error is known.

| metric | scale | direction (benefit) | margin | scientific reasoning |
|---|---|---|---|---|
| T1 mean depth | levels | DiD negative | **0.10** | ~7% of the ~1.45 realized depth; below the span ordinary `k`/`base` jitter covers, so a smaller effect could not change advice about which surface to use. It sits between the critic's measured DiDs for A1 (~0.21) and A2 (~0.03), so it discriminates |
| T1b mean leaves | leaves | DiD negative | **0.15** | ~6% of the ~2.5 realized count; monotone companion to T1 |
| log E2 | log inefficiency | DiD negative | **0.223** = log 1.25 | 25% in an inefficiency factor is ~1.56x in ESS - the smallest change that would alter advice about run length. Below it, run 25% longer |
| log E3 | log IACT | DiD negative | **0.223** | same |
| log E4 | log IACT | DiD negative | **0.223** | same |
| G1 held-out RMSE | relative | X above A0 is harm | **0.060** | the record's frozen C2 margin |
| G2a mean sigma | absolute | X above A0 is harm | **0.05** | 5% at sigma = 1, just above the 0.999-1.032 arm spread measured |

### KILL-1 - the transfer does not happen

**KILL** if, on the primary pair **P1 = (D0, D50)**, NO composition arm in
{A1, A2, A3} satisfies BOTH

  (a) `DiD_X(T1) <= -0.10` with the one-sided test rejecting at alpha 0.05
      (fallback: `DiD_X(T1b) <= -0.15` if T1 exits at FREEZE), and

  (b) at least one of `DiD_X(log E2) <= -0.223`, `DiD_X(log E3) <= -0.223`
      (or `log E4` if both exit), same test,

at the frozen replication, **and** the confirmation re-run below also
fails.

**Confirmation re-run, mandatory before KILL-1 fires.** KILL-1 fires on the
*absence* of a flag, so a fresh-seed re-run of a flagged cell offers it no
protection. Instead: the arm with the most negative `DiD_X(T1)` is re-run
on a fresh seed block at R = 12, and KILL-1 fires only if that arm's
re-run point estimates also fall short of both limbs' margins.

This instantiates the record's kill line - "the transfer either does not
shrink the forest's job or does not translate into better tree-space
behaviour when it does" - with the record's raw cross-arm contrast replaced
by the dial DiD, because the raw contrast was measured to be 2.5x the
transfer effect it claims to measure.

**Corroboration on P2.** For any arm carrying a GREEN, the sign of
`DiD_X(F, P2)` must agree with `DiD_X(F, P1)` on both limbs. Disagreement
does not kill; it returns YELLOW with the SNR-versus-scale confound as the
finding.

**Gated-secondary: Q1.** The same two limbs, on pair (DQ0, DQ50), for A2
against A0 only. Q1 tests the mechanism's **strong** case - the smooth
nonlinear share, the one with the most tree renderings - where the primary
span tests its weak case. Q1 does not enter KILL-1: a Q1 flag with a
primary GREEN scopes the GREEN to the smooth share and says so; a Q1 flag
with a primary KILL is reported as the finding that the mechanism is
share-shaped.

### KILL-2 - REMOVED

**There is no ridge kill.** See "The ridge readouts". The record's second
kill line is not executable in this design, and the inner-versus-outer
question is not adjudicated by this probe. Registered as a scope loss, not
hidden.

### Harm - mandatory, not confirmatory

**KILL the composition recommendation entirely**, whatever KILL-1 returns,
on a per-cell harm rejection on **G1** or **G2a**: a one-sided harm test
(H0: contrast <= 0) rejecting under Holm across all G1/G2a x arm x cell
contrasts at family-wise 0.05, AND a point estimate beyond its margin, in
at least two cells. **One cell rejecting triggers a fresh-seed re-run of
it; confirmation counts as the second** (`-study.md:160-161`).

v1 dropped that escape hatch while citing, as validation, an incident in
which **exactly one cell flagged** (S4-5000, C2 = +11.10%) and the re-run
supplied the second (`grow-from-root-default.md:730`). Under v1's clause
that KILL would never have fired. Restored verbatim.

### GREEN, and exactly what it licenses

**GREEN** iff, on P1: KILL-1 does not fire; no harm rejection on G1 or G2a;
every validity gate holds; the P2 corroboration agrees in sign for the
arm(s) carrying the GREEN; and at least one composition arm shows both
limbs. **YELLOW** otherwise: tables to VD.

A GREEN **licenses**:

1. **The measurement itself, published in-repo.** It is the first mixing
   diagnostic for a composed parametric-plus-BART sampler anywhere, and it
   is the first ridge measurement on one. It informs stan4bart sec 4.5's
   open question and CSP-BART's opposite position without settling them.
2. **Documentation.** Composition guidance in `docs/design/linear-leaves.md`
   and on the `?dbartsPriors` / vignette surface, describing a surface that
   already ships. No code.
3. **A re-prioritization input** for the tree-space program, per the
   record's sec 7: "if the composition probe says depth transfer works and
   survives the ridge, then every tree-space candidate is worth less than
   it looks." It re-ranks the perturb study (sec 4.2 / 6) and the move
   census against each other. It authorizes neither.

A GREEN **does NOT license**:

- **Any engine change.** Not a default flip, not a new kernel, not the
  collapse (the 800-1300 line hot-path door already priced and declined in
  `forest-ranef-interweaving.md`).
- **Any per-second or cost claim.** No timing readout exists here, by
  construction. A1 costs ~8x A0 per sweep; any guidance following a GREEN
  must be paired with a separate cost measurement on a quiet machine.
- **Any preference between the inner and outer variants.** The ridge is
  reported, not gated, and no cell calibrates a reference against which to
  rank them.
- **Any generalization beyond a normal linear block on two designated
  continuous columns under this DGP family**, at n = 5000, p = 10, m = 75,
  sigma = 1, with overlapping covariate sets.

## Power and replication

Design condition, replacing v1's `4 x SE <= margin` (which delivers **50%**
power at a true effect equal to the margin, not the 0.99 v1 claimed - power
is not defined at the null, and the specificity that rule buys is
`1 - P(Z > 4) = 0.99997`):

```
(z_0.95 + z_0.90) x SE <= margin,   i.e.   2.926 x SE <= margin
R >= 9 x sd_rep^2 / margin^2        alpha 0.05 one-sided, power 0.90
```

**Planned R = 24 x C = 8 chains** per (arm, cell), matching the record's
SMALL stratum exactly (24 x 8; `-study.md:211-212`,
`grow-from-root-default.md:463`). v1's C = 4 contradicted its own two
citations of 24 x 8 and would have inflated every between-chain sd_rep by
~1.53x.

| metric | prior on sd_rep(DiD) | source | margin | status at registration |
|---|---|---|---|---|
| T1 | **none** | - | 0.10 | Stage 0. Expected comfortable: the statistic averages ~1.2M tree depths per (arm, cell, replicate) and the DiD pairs on a shared X and a shared noise draw |
| T1b | **none** | - | 0.15 | Stage 0 |
| log E2 | **none** | - | 0.223 | Stage 0. **The binding risk in the design** |
| log E3 | **none** | - | 0.223 | Stage 0 |
| log E4 | **none** | - | 0.223 | Stage 0 |
| G1 | 0.03998 relative | record SMALL C2, the study's **frozen Stage-0 table** (`grow-from-root-default.md:399`) | 0.060 | `R = 9 x 0.040^2 / 0.060^2 = 4.0`. **Registered as passing at R = 24** |
| G2a | **none** | - | 0.05 | Stage 0 |

**No sd_rep prior is imported except G1's**, and its provenance is the
record's re-measured frozen value. v1 imported 0.0011-0.0015 for T2 from a
scratch probe the record explicitly superseded (`-study.md:28-30`), read
two strata as a range, and paired a LARGE-n sd with a SMALL-n R; at the
record's own re-measured LARGE value with the four-chain inflation, that
metric failed floor-above-ceiling before Stage 0 ran.

**Seed and pairing discipline.** Within replicate `r`: draw `X` **once**
and the noise vector `eps` **once**, from `set.seed(BASE_SEED + r)`;
construct `y` at every dial level from that same `X` and `eps`. Sampler
seed `= r` for every arm and every cell in the replicate (A3:
`r * 1000L + chainIndex`). This matches the record's both-seeds pairing
(`-study.md:81-82`). Sampler streams are not matched *across arms* - the
models differ, so they cannot be - but they are matched **across dial
levels within an arm**, which is exactly what the DiD needs and what
shrinks `sd_rep(DiD)`. Fresh-seed blocks for re-runs: `BASE_SEED + 1000 + r`.
Replicates are independent; parallelism is over replicates; two invocations
at the same `BASE_SEED` produce identical tables.

**Escalation, with a cost ceiling.** For any metric whose Stage-0 sd_rep
exceeds `margin x sqrt(R) / 2.926`, either R rises to
`ceil(9 sd_rep^2 / margin^2)` or the metric leaves the gate under
floor-above-ceiling. Which, per metric, is decided at FREEZE and written
here. **Ceiling: `R_max = 48`, and a total escalation budget of +10 h
single core, whichever binds first.** Beyond either, the metric exits.

**Registered fallback ladder, so no gate can become unevaluable by
surprise:**

- limb (a): T1, then T1b. If both exit, KILL-1 runs single-limb on (b) and
  the verdict is scoped to mixing transfer without the forest-shrinkage
  leg, declared at FREEZE.
- limb (b): E2, then E3, then E4.
- **If limb (b) has no surviving metric, KILL-1 is VOID and the probe is
  declared an ESTIMATION STUDY at FREEZE, in writing, before any Stage 1
  contrast is read.** Under that declaration it reports DiDs with
  confidence intervals, licenses documentation and re-ranking inputs, and
  issues no kill and no GREEN. This is a registered fork, not a deviation.

**Spurious-kill probability of the composite.** KILL-1 errs by firing when
a real transfer exists, and that error re-ranks the whole tree-mixing queue
against composition on a false negative. At a true transfer exactly equal
to the margins in all three arms, treating the two limbs as independent
(the conservative direction, since positive correlation raises joint
power):

- per arm, both limbs: 0.90^2 = 0.81; arm fails: 0.19
- all three arms fail: 0.19^3 = **0.0069**
- times the confirmation re-run's failure probability (~0.5 on a
  point-estimate criterion at R = 12): **P(spurious KILL) ~ 0.004**
- at 1.5x the margins: **< 1e-5**

Limb (b) is an OR over three metrics, which only raises power, so 0.004 is
an upper bound.

**Specificity.** Under a true DiD of zero everywhere, a limb clears only if
it rejects at alpha 0.05 *and* the estimate exceeds the margin; the latter
has probability `P(Z > 2.926) = 0.0017`. Per arm both limbs ~3e-6, so
`P(GREEN | no transfer anywhere) ~ 2e-5`: **KILL-1 fires with probability
~0.99998 when the transfer is truly absent.** That is what makes this a
falsifier rather than an estimation study, and it is conditional only on
limb (b) surviving FREEZE.

**Multiplicity: none on the benefit side, deliberately.** The false-GREEN
probability across three arms x two limbs at these thresholds is ~2e-5. The
risk in this design is a false KILL, so the correct response is power, not
correction. The harm clause is Holm-corrected because its error direction
is the opposite.

## Harness, cost, machine

```
.claude/composition-probe/harness/
  common.R    dial constructors, cells, seeds, depthsPreorder, treeSummary,
              E2/E3/E4, ridge projections, iact, checkpoints
  arms.R      runA0 runA1 runA2 runA2r runA2c runA2f runA3 - one signature,
              one exposure
  stage0.R    V0-V10, sd_rep(DiD) for every readout at R = 8, the A2r
              rescale calibration, the IACT stability report
  freeze.R    writes every frozen margin, escalation and gate-exit into
              THIS file
  stage1.R    D0 / DP / D50 at full R; then DQ, D100, DR, DC, DF, DJ, DG,
              X, N
  analyze.R   per-cell tables, DiD beside its margin, verdicts
  results.md  tables + the verdict against this file
```

Cost basis, measured at n = 5000, p = 10, m = 75, 1 thread
(`feasibility.md` sec 8) - **planning input only; it appears in no
criterion and in no results table**: A0 1.05 ms/sweep, A0 + `getTrees`
1.25, A1 8.23, A2 1.08, A3 2.22 ms/iter, depth extraction 0.062 ms/sweep.
Effective per-sweep with `getTrees` and the depth walk on kept sweeps only:
**A0 1.18, A1 8.36, A2 / A2r / A2c / A2f 1.21, A3 3.50** (A3 raised from
2.22 to price `extract(fit, "trees")`, which materializes ~706k data-frame
rows per single-chain fit). At C = 8 and B = K = 2000 a replicate is 32,000
sweeps per (arm, cell).

Two items v1 left unpriced, now priced:

- **the A1 ridge reconstruction**: 2-4 ms/sweep at n = 1000, so 10-20 at
  n = 5000. Mitigated by the `n_r = 1000` row subsample and by running the
  ridge readouts only at D50 on a pre-registered `R_ridge = 8` subset;
  48 s/replicate.
- **A3's tree and `indiv.bart` extraction**: mitigated by the
  eight-single-chain-fits decision, which caps `indiv.bart` at 80 MB
  instead of 640 MB per replicate. `indiv.bart` is extracted **only** in
  the ridge subset.

| stage | content | est. wall (1 core) |
|---|---|---|
| Stage 0 | V0-V10; sd_rep(DiD) at R = 8 on D0 + D50, four arms, full scale; A2r; the IACT stability report | 2.5 h |
| FREEZE | write every margin, escalation and gate-exit into this file | - |
| Stage 1 main | D0 / DP / D50, 4 arms x 8 chains x 4000 sweeps x 24 reps | 9.1 h |
| ridge subset | D50, R = 8, all arms + A2c | 0.35 h |
| A2r calibration | DR0 / DR50, R = 8 | 0.17 h |
| Q1 gated-secondary | DQ0 / DQ50, A0 + A2, R = 24 | 1.0 h |
| reported cells | DF, DJ, DG | 1.2 h |
| structure cells | X, N at m = 1, R = 8 | 0.34 h |
| D100 | validity check, R = 4, reduced window | 0.10 h |
| KILL-1 confirmation | best arm, fresh seeds, R = 12 | 2.0 h budget |
| harm re-runs | fresh-seed confirmations | 1.5 h budget |
| report | tables, verdicts | 0.5 h |
| **total** | | **~19 h single core** |

Embarrassingly parallel over replicates: **~3 h on 8 cores** at the ~6.3x
efficiency the cost audit verified. With full escalation the hard cap is
~29 h single core, ~4.6 h on 8 cores. The measurement is load-insensitive
by construction - every gated quantity is a function of the draw sequence
and no readout is wall-clock - so it runs on the laptop under any load. It
must not share the machine with anything that IS timing-gated.
`dbarts-bench` is optional and buys nothing here.

**Memory.** ~0.5 GB per worker outside the ridge subset, up to ~1.2 GB
inside it. Registered: cap the ridge subset at 4-way parallelism on a
machine with under 16 GB.

**Cut order if the budget binds**, registered in advance: DG, then DJ, then
DF, then X/N. The gated cells (D0, DP, D50, DQ) and the calibration cells
(DR, D100) are not cuttable.

## Instrumentation

### Ships today; sufficient for every arm; no engine change

| need | shipped hook |
|---|---|
| A0 / A1 / A2* sampler | `dbarts(formula, data, control, seed, node.prior, resid.prior)`, `R/dbarts.R:328` |
| linear leaves | `node.prior = linear(c("x4","x5"))`, resolved by `parsePriors` (`R/model.R:90-127`, `:1076-1085`); **call position on `dbarts()` only** |
| fixed sigma (cell DF) | `resid.prior = fixed(1)` (`R/model.R:959`, in `dbartsPriors` at `:1083`), verified to work alongside `node.prior = linear(...)` |
| one sweep | `sampler$run(0L, 1L)` -> `sigma, train, test, varcount, k, varprobs, tau, ranef`, `R/dbarts.R:755` |
| offset exchange | `sampler$setOffset(offset, updateScale, updateState)`, `R/dbarts.R:1004` |
| the A2r schedule | the same call with `updateScale = TRUE`, verified to succeed on a plain sampler at the live tip |
| live tree structure | `sampler$getTrees(current = TRUE)` -> pre-order frame; `beta.<col>` columns under a linear leaf |
| independent depth walk (V0) | `dbarts:::getTreeDepthAndSize`, `R/plotTree.R:39-55` |
| A3 | `stan4bart(...)`, 0.0-14 on branch `bartcore` at `6ce0440` |
| A3 trees | `extract(fit, "trees")` - the identical six-column frame |
| A3 fits, sigma, varcount, fixef, ranef | `extract(fit, type = )` in `{indiv.bart, sigma, varcount, fixef, ranef}` (there is no `"bart"` type) |

**Not usable, registered so nobody starts:** stan4bart's WALNUTS
diagnostics. `divergent__`, `treedepth__` and `energy__` are constant-zero
placeholders and `check_sampler_diagnostics` is a documented no-op
(`stan4bart/R/stan4bart.R:280-295`).

### Reused from the grow-from-root harness

Present at `.claude/grow-from-root-default-research/harness/` - **gitignored,
so recoverable on this machine but not from git**;
`grow-from-root-default.md` sec 8 is correct that the supported path is
reconstruction from the pre-registration. Read and confirmed portable:
`friedman()` (`common.R:56`), the S5 / S6 constructors (`cells.R`),
`inclusionByChain()` (`:209`), `iact()` (`:244`), the root-extraction idiom
inside `rootStats()` (`:307`), the test helpers `harmTest` / `benefitTest`
(`:360-385`), the checkpoint helpers `ckptPath` / `withCkpt` (`:386-390`).

**Not portable:** `fitArm()` (`:150`), `armStats()` (`:257`) and
`rootStats()` (`:307`) are written against `bart2()` fit objects; this
probe drives `dbarts()` sampler objects sweep-by-sweep because A1 is
unreachable from `bart2()`. Their bodies port; their interfaces do not. The
old harness sources `common.R` by a hardcoded absolute path; the new one
must not.

### Built new (R only)

1. the seven-arm sweep driver with per-arm hooks (~180 lines);
2. `depthsPreorder()` and the per-sweep tree summary (~50 lines);
3. the conjugate normal block, its componentwise variant, and the A2r
   schedule (~40 lines; the joint version is prototyped and working);
4. E2 / E3 / E4 (~60 lines);
5. the ridge decomposition, including the A1 linear-leaf reconstruction
   with per-observation leaf assignment and the V8 identity (~140 lines -
   more than v1's budgeted 80, because the reconstruction needs leaf
   assignment through each tree on the standardized-covariate and
   internal-response scale);
6. the stan4bart adapter, with the `varcount` reshape and the
   `skip`/`cores` assertions (~70 lines);
7. the dial, Q1, DJ, DG, X and N constructors (~100 lines);
8. V0-V10 (~180 lines) and the DiD analysis / verdict layer (~300 lines).

## Steps

1. **`R CMD INSTALL . --preclean`**, and record `git rev-parse HEAD` in
   `results.md`. `src/Makevars` carries no header dependency tracking (no
   `-MMD`, no `.d` files) and `CLAUDE.local.md`'s own gotcha is that stale
   objects against changed `facade.hpp` virtuals bus-error; the probe is
   sequenced after another measurement, so more header commits will land
   before it runs. v1's justification for skipping ("nothing in `src/`
   changes") was true of the probe and irrelevant to the build state.
   **The whole probe runs against one installed build**; if the tree moves
   mid-run, restart or record a deviation.
2. Run stan4bart's full tinytest suite against that build (V4).
3. Write `common.R` and `arms.R`; run V0, V1, V2, V3, V5, V7, V8, V9 on one
   small cell each. Any categorical failure stops the run.
4. Run V10 (D100 collapse).
5. **Stage 0**: D0 + D50 at R = 8, full scale, four arms - `sd_rep(DiD)`
   for T1, T1b, log E2, log E3, log E4, G1, G2a; the A2r rescale
   calibration on DR0 / DR50; the ridge readouts and their seed-to-seed
   spread; V6 reproducibility; the `V0inv` negligibility assertion; the X
   and N reported statistics with `n_cond`.
6. **FREEZE.** Write into this file, before reading any Stage 1 contrast:
   every `sd_rep(DiD)`; which metrics are gated and which exited under
   floor-above-ceiling; any escalated R and the budget it consumed; the
   D100 threshold; the A2r verdict on A3's T1/T1b eligibility; and, if
   limb (b) has no surviving metric, **the estimation-study declaration**.
7. **Stage 1**: D0, DP, D50 at the frozen R; then DQ0/DQ50; then DC, DF,
   DJ, DG, X, N.
8. Apply, in order: the validity gates, KILL-1 with its confirmation
   re-run, the P2 corroboration, Q1, and the harm clause with its
   fresh-seed re-run.
9. **Report**: per-cell tables printing each DiD beside its margin and
   Holm-adjusted p; the T1 / T1b / E2 / E3 / E4 panel per arm per cell; the
   **calibration table** (each arm's D0 offset, and what fraction of the
   naive raw D50 contrast it accounts for); the ridge panel W1-W6 with the
   MATERIAL / IMMATERIAL / INDETERMINATE reading per arm; G1-G4 as the harm
   panel; X/N with `n_cond`; the reported cells; and the verdict, with the
   GREEN-licenses and GREEN-does-not-license text quoted verbatim.

## Verification

- Two invocations at the same `BASE_SEED` produce identical tables for
  every dbarts-side arm; A3 agrees bitwise on `extract(fit, "sigma")`.
- `#nodes == 2 * #leaves - 1` for every tree of every arm, and the depth
  walk agrees exactly with `getTreeDepthAndSize`.
- Realized sweep counts match `B + K` exactly (`B + K + 1` for A3);
  `skip == 1L` and `cores == 1L` on every A3 fit.
- V8's reconstruction identity is exactly 1 on every checked sweep.
- V10 holds; otherwise the harness is wrong and the run stops with a
  finding rather than a verdict.
- The harness exits nonzero only on a validity-gate failure, never on a
  study finding.

## What this measurement cannot tell us

1. **Whether a fixed posterior is sampled faster.** A0 and A2 target
   **different models**. The claim on offer is the within-arm sensitivity
   stated in the Goal, not "BART mixes better". The guardrails exist
   because a slope difference across different targets is only interesting
   if the fit is not worse.
2. **Which composition variant is better.** The ridge is reported, not
   gated, and no cell calibrates a no-ridge reference. Inner versus outer
   is not adjudicated here.
3. **Anything per second.** No timing readout exists, deliberately. A1's
   measured ~8x per-sweep cost means a per-sweep GREEN is not a
   recommendation.
4. **Whether HMC helps.** With `Z = [x4, x5]` the conjugate draw in A2 is
   exact and joint, so A3 has no within-block ridge to remove. A3 measures
   the shipped ecosystem pattern. DG is the ungated bridge to the regime
   where HMC would pay, and it compares two different model classes (A0
   takes `g` as a 10-level factor predictor; A3 takes it as a random
   intercept), so it is reported for that reason too.
5. **Whether the forest-versus-parametric ridge is as bad as the
   forest-versus-group-intercept ridge.** The house 6x is measured on a
   per-group intercept. `tree-mixing-proposals.md` sec 9 declines to bound
   this in either direction and so does this probe. W1 is not on a common
   scale with the 56.1 / 9.3 surrogate numbers.
6. **Anything about the collapse.** No arm here is a collapsed sampler.
   That door stays priced and declined.
7. **Structural mode collapse at m = 75.** House law 1. T3 at m = 75 is
   reported, never gated; the m = 1 cells are reported for the four reasons
   listed under Scenarios.
8. **Structural acceptance behaviour.** Not exposed by any shipped surface;
   the record's primary readout is deleted with that reason.
9. **Generalization past this DGP family.** One interaction plus one smooth
   univariate plus two exactly linear terms, at one `(n, p, m, sigma)`.
   Nothing here speaks to high-dimensional Z, categorical predictors,
   non-Gaussian families, or a misspecified parametric block.

## Thresholds and claims that are judgment calls

Flagged rather than buried. Mechanical items (the dial coefficients, the
variance components, matched exposure, the power arithmetic) follow from
arithmetic above.

1. **JUDGED - `Z = [x4, x5]` as the primary absorbable span.** It is the
   only span with no basis-specification freedom and the only one A1 can be
   designated on without changing the design matrix, but it is the
   mechanism's *weak* case. Q1 is gated-secondary precisely for that
   reason. An equally defensible design makes Q1 primary.
2. **JUDGED - the P1/P2 split.** The trilemma is real and the choice of
   which pair is primary is a choice. P1 is primary because the DiD's
   cancellation argument depends on the un-absorbable target being
   identical, which is what P1 buys.
3. **JUDGED - the T1 margin of 0.10 levels, and 0.223 on every log
   efficiency.** Both are substantive judgments about what would change
   advice, not measurements. Both are written before Stage 0, which is the
   point.
4. **JUDGED - the ridge interpretation ladder (50 / 20 sweeps, W3 >= 3).**
   Keyed on effective sample size in the kept window. No measurement
   anchors the exact numbers; nothing gated depends on them.
5. **JUDGED - `V0inv = diag(1e-2, 2)`.** Negligible at n = 5000; asserted
   rather than swept.
6. **JUDGED - overlap ON in the primary specification.** Forced by A1. The
   choice to gate the overlapping configuration is still a choice, and it
   is the configuration the literature disagrees about. Note the measured
   finding that removing shared covariates does **not** remove the ridge,
   which means overlap-versus-disjoint measures less than v1 claimed.
7. **JUDGED - `B = 2000` burn.** `bcf-sigma-residual.md`'s burn curve shows
   the structural bottleneck persisting to 36k-72k sweeps in the
   strong-scale regime. This probe is not in that regime, and measured
   sigma sits at its plateau in all arms, but 2000 is a convention. G2b
   reports the sigma trajectory split at warmup iteration 1920, where A3's
   last rescale fires, so a scale artifact is visible.
8. **JUDGED - X and N reported rather than gated, and DG comparing two
   model classes.** Both defensible the other way at higher cost.
9. **UNVERIFIED - V4 has still never been run.** "No stan4bart rebuild
   needed" rests on `dbarts.h` being byte-stable since `876a339` plus one
   smoke fit. `9929ede` changed engine behaviour without changing the ABI.
   V4 after a `--preclean` install is the real check.
10. **UNVERIFIED - the second-order part of A3's rescale.** The first-order
    dial-invariance argument is arithmetic; the warmup-time realized anchor
    is not, and A2r is the measurement that decides it.
11. **UNMEASURED - every effect size in this file.** The critic's numbers
    (A1 depth DiD ~ -0.21, A2 ~ -0.03; E-family none) are at n = 1000, 2-4
    seeds, and are sufficient to refute point claims and to show sign and
    rough size. They are **not** effect sizes, **not** sd_rep estimates, and
    nothing here may be frozen from them.

## Deviations from this pre-registration

(none yet - append here, with date and reason, before reading any gated
metric)
