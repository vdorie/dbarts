# probit-k-scale-move: a parameter-expansion step for k under probit

Status: PLANNED 2026-10-08 (dec-B371; TODO probit-k-mixing). Queued behind the engine small-rulings slice
(Interactions).

agent: opus implementer, one; blind critique of this plan first; one opus reviewer who runs the mutants below.
rng: POSTERIOR-CHANGING by the README's rule (a default changes: every single-forest probit fit with a drawn k
takes a new kernel); the stationary distribution is unchanged by construction and the exact gate is what shows
it. NEUTRAL, bit for bit, everywhere else (Scope).
window: before 1.0-0; engine slices stay serial.
budget: ~1750 lines (code ~300, tests ~330, gate and measurement ~420, docs ~380, help and records ~120; the
re-recorded baseline is data).

## Goal

Once a sweep, before the trees, a single-forest probit fit with a drawn k multiplies every active latent and
every occupied leaf by one factor alpha and divides k by it, alpha drawn from its exact conditional. k's
autocorrelation time on the probit-k arm's datasets falls from a median of 524 sweeps to about 71 and its 99th
percentile from 23,000 to about 800, the separated toy's small-k mass is reached at its known value, and the
probit-k calibration arm reads as an ordinary SBC again. `dbartsControl(treeScale = "never")` turns it off. The
tier is "Changes draws" ([Process by risk](README.md#process-by-risk)).

## Context

- dec-B371 rules a scale move for k under probit before 1.0-0, chosen by k's autocorrelation time over the arm's
  datasets and small-k visitation on the separated toy. The design compared interweaving (a), the parameter
  expansion carrying k (b) and the collapsed move (c) on a prototype; (b) wins, (c) is a door, (a) is dropped.
  The design and its evidence are in scratch/k-move/design.md (runs in scratch/k-move/out; the prototype patch
  0001-*.patch is on bartcore 0386245f-era code, switched by an environment variable, not for landing). Per
  dec-B389 the design and its evidence move into docs/design/probit-k-scale-move.md with the code, and the
  measurement scripts the evidence rests on move into benchmarks/R.
- The finding: [probit-k-calibration.md](../design/probit-k-calibration.md).
- The sweep: [One sweep](../architecture.md#one-sweep). The level step,
  [`drawLevelShift`](../../src/bartcore/chain.hpp), is the model for placement, eligibility and "a forest that
  skips consumes no generator draw". The cache rule ([Forests and combiners](../architecture.md#forests-and-combiners)):
  a multiplicative transform re-derives `totalFits` from the leaves before the sweep ends, which
  [`assertForestCachesMatchLeaves`](../../src/bartcore/chain.hpp) checks in debug builds.
- The k draw: [`ChiKHyperprior`](../../src/bartcore/model.hpp) (the engine's one k hyperprior; `invchi()` reaches
  it too). Defaults: [`drawsLeafKByDefault`](../../R/spec.R) draws k for probit, logistic and nbinom.
- The switch's idiom: `treeShift`, through [`resolveTreeShift`](../../R/dbarts.R),
  [`controlArgumentFromSlot`](../../R/dbarts.R), the `levelGibbs` slot in [A_class.R](../../R/A_class.R),
  [`parseControl`](../../src/R_interface_bartcore.cpp), [`printInitialSummary`](../../src/R_interface_bartcore.cpp),
  [`optionsFromParsed`](../../src/R_interface_bartcore.cpp) and [`SamplerOptions`](../../src/bartcore/chain.hpp).

## The move

Model, one probit forest: k ~ s chi(nu); each occupied leaf mu ~ N(0, (c / k)^2), c the leaf scale; latent
z_i ~ N(o_i + f_i, 1) on the active rows, y_i = 1{z_i > 0}. The group alpha > 0 acts as
(z, mu, k) -> (alpha z, alpha mu, k / alpha) on the active latents, the occupied leaves and k. With n active rows,
R = sum_i (z_i - f_i)^2, Q = sum_i o_i (z_i - f_i) over active rows, and C = k^2 / (2 s^2) (0 under an infinite
scale), the Jacobian alpha^(n + M - 1), the leaf prior's alpha^-M and the k prior's alpha^-(nu - 1) give, in
v = log alpha,

    log p(v) = (n - nu) v - (R / 2) e^(2v) + Q e^v - C e^(-2v)

The truncation is invariant for alpha > 0; empty leaves sit at zero and stay there; an inactive row is not in the
model, so its latent is neither counted nor scaled. Without an offset alpha^2 is generalized inverse Gaussian;
with one it is not, so v is drawn by one slice step with stepping out (Neal 2003) from v = 0, a fixed width, each
evaluation O(1).

Exactness, in brief: the step is a Gibbs draw of the orbit coordinate given the orbit-invariant coordinates
(k z, k mu, trees), the generalized Gibbs step of Liu and Sabatti (2000) on a group with Haar measure
d alpha / alpha; a slice step is a reversible kernel for that one-dimensional conditional, so the composite leaves
the posterior invariant. Every decline below depends only on orbit-invariant quantities, so mixing the step with
the identity keeps it invariant. Placed ahead of the tree loop, beside the level step, every recorded channel
(training and test fits, saved trees, k) is written after it. The design note carries the full derivation,
why (a) cannot reach the slow cases, why carrying k matters (b0), and (c) as a door.

## Scope

Taken every sweep, per chain, when all hold (read before any generator call; otherwise nothing is drawn):
`treeScale` is "auto"; one forest, no combiner, no variance forest; family probit; the plain constant leaf (not
monotone, linear or gp); the forest draws k and k is finite and positive; no tree's map is stale
([`leafOfStale`](../../src/bartcore/combiner.hpp), the sweep after `sampleTreesFromPrior` or a wholesale reset);
R > 0; and, under an infinite prior scale (C = 0), n > nu, where the conditional is otherwise improper.

Reached through: `bart`, `dbarts` and `xbart` on a 0/1 response, `hazard` (one forest), the zero part of
`hurdle.lognormal`, `rbart_vi`'s binary fit (its intercepts arrive as the offset, the Q term) - each where k is
drawn, the default. Measured on the current tip: modern `bart` and `dbarts` draw k on a binary response, as do
the hazard fit and the hurdle zero part; `bartBT` (and `bart` called with BayesTree names, forwarded to it)
fixes k.

Bitwise unchanged: every other family (logistic, nbinom, gaussian, Student-t, aft, ordinal, multinomial); probit
with a fixed k (bartBT's default included); BCF, amplitude-coupled and multi-forest hazard fits; linear, gp and
monotone leaves; `treeScale = "never"`; grow-from-root's sweeps (an initializer; the move is in the sampling
sweep only). Declines inside an eligible fit (stale maps, every row inactive, k infinite) consume no draw at that
sweep.

## Changes

Engine (chain.hpp, model.hpp; names are proposals, cited by symbol once landed):
- `SamplerOptions` gains `bool scaleExpansion = true`, commented beside `levelGibbs`.
- A free function template `sliceFromZero(ext_rng*, logDensity, width)` with a named width constant (0.25, the
  prototype's); the shrink loop terminates because v = 0 is always in the slice.
- `Chain::scaleExpansionApplies(forest)` (the predicate above, no generator) and `Chain::drawScaleExpansion(forest)`,
  called in `runSweeps` after the level step. One pass over the rows for n, R, Q (reading `latents()`, `offset()`,
  and `workingWeights()`, which is the mask under probit); one slice draw; then each tree's occupied bottom nodes
  multiplied by alpha (walked as `drawLevelShift` walks them), k divided by alpha, the latents scaled, and
  `totalFits` re-derived from the leaves ([`rebuildTotalFitsFromTrees`](../../src/bartcore/chain.hpp)), never
  multiplied in place.
- [`ResponseModel`](../../src/bartcore/model.hpp) gains `scaleLatents(double)`; the base refuses (unreachable past
  the predicate); [`ProbitResponse`](../../src/bartcore/model.hpp) multiplies its active latents and rebuilds the
  working response (working = alpha z - o, not alpha (z - o)), in place, so the sweep's `y` pointer stays valid.
- No new chain state, no state-format change, no facade virtual, nothing in dbarts.h (no ABI event).

Bridge: `parseControl` reads the slot into `ParsedControl`, `optionsFromParsed` into the option,
`printInitialSummary` prints "scale expansion step: auto|never" after the level-step line.

R: a `scaleExpansion` slot, logical TRUE ("auto") or FALSE ("never"), plain TRUE/FALSE validity; `dbartsControl`
gains `treeScale = c("auto", "never")` after `treeShift`, resolved by a `resolveTreeScale` beside
`resolveTreeShift`; `controlArgumentFromSlot` maps it back; `$setControl` adds it to the creation-fixed list
(restating the value is accepted, as for treeShift).

Help: [dbartsControl.Rd](../../man/dbartsControl.Rd) usage and a `treeScale` item (what the step does, that the
posterior is the same, where it applies, fixed at creation); the control-settings lists in
[bart.Rd](../../man/bart.Rd) (both places), [xbart.Rd](../../man/xbart.Rd) and the setControl paragraph in
[dbartsSampler-class.Rd](../../man/dbartsSampler-class.Rd).

## Tests

tests/cpp, in test_ensemble.cpp (beside the level-step conditional test, reached from
[`runEnsembleTests`](../../tests/cpp/test_ensemble.cpp)), through [`TestPeer`](../../tests/cpp/test_peer.hpp):
1. Conditional: a burned-in probit drawn-k sampler with an offset and a mask, state frozen and restored between
   2e5 draws; the drawn log alpha against the target CDF of the four constants by trapezoid quadrature, KS bound
   at p = 1e-3 (D < 1.95 / sqrt(N)).
2. Mapping, one move: k' alpha = k; mu' = alpha mu on occupied leaves, empty leaves 0; z' = alpha z on active
   rows, inactive latents unchanged; working = z' - o; `totalFits` equal to the tree-order sum bitwise; ratios to
   1e-14. This, not (1), catches k or the latents left unscaled after a correct draw.
3. Inertness: option off, fixed k, logistic, a stale map, every row inactive, k infinite, and C = 0 with n <= nu
   consume no generator draw (the next uniform is unchanged).
4. The slice helper on a one-dimensional density with a closed-form CDF, KS.

tinytest, a new file inst/tinytest/test-probit-scale-move.R: vocabulary
(the two words, a partial match, NA, length 2 and other words refused naming `treeScale`), the default, the
setControl refusal and restating, a control rebuilt by a front door keeps it; "auto" and "never" draw bitwise the
same on gaussian, logistic, nbinom, probit at fixed k, a linear-leaf probit with a drawn k, multinomial and a BCF
probit fit; they differ on a drawn-k probit fit through bart, dbarts, hazard and the hurdle zero part; an offset,
a mask installed mid-run and chi(1.5, Inf) run finite; storeState and setState, copy and a saveRDS reload
continue. Entries in [test-argument-surface.R](../../inst/tinytest/test-argument-surface.R),
[test-control-valuesAreUsed.R](../../inst/tinytest/test-control-valuesAreUsed.R) and
[test-na-as-none.R](../../inst/tinytest/test-na-as-none.R) beside treeShift's.

## The exact gate

benchmarks/R/probit-k-scale-exact.R, from scratch/k-move/gate.R and gate-fix.R. One tree, structure frozen,
probit, k ~ chi(1.5, 2), n = 150, the move on. Given the tree, p(k | y) is exact by one-dimensional quadrature a
leaf; the posterior of k is integrated over log k by the trapezoid or midpoint rule (counting whole grid cells
biased every arm by 0.001 to 0.003 and put correct moves at z -3 to -5). Arms: pure (every leaf separated, the
small-k mass), mixed (none separated), offset (o = 0.4, the Q term), mask (a fifth of the rows inactive).
Statistics: P(k < q_j) at the exact deciles 0.1, 0.25, 0.5, 0.75, 0.9 and E[mu_l] on unseparated leaves, each a
z against a 50-batch batch-means error, failing above 4.5. Mixing bound: the pure arm's batch-means error of each
P(k < q_j) under 0.01 in full mode (prototype: 0.0039 with the move, 0.0194 without); quick mode's bound is set
between the two from quick runs at landing and recorded in the header, so a move silently off fails per push.
Quick 1e6 sweeps an arm (~30 s), full 4e6 (~2 min); a `never` argument runs the move off, for the reviewer, and
must fail. Added to [exact-gates.yaml](../../.github/workflows/exact-gates.yaml)'s list and to the gate list and
a paragraph in benchmarks/README.md.

Prototype verdict to reproduce: with the move |z| <= 1.9 on all 20 decile statistics and every E[mu]; without
it the pure arm reaches z -5.2.

## Measurement script (dec-B389)

benchmarks/R/probit-k-mixing.R (measurement, not a gate), from mix.R, iat.R and stat.R, driven by the landed
control in place of the environment variable: `census` (k's and avg.f's integrated autocorrelation time over the
arm's prior-drawn datasets, the arm's start, 30,000 burn and 200,000 recorded sweeps) and `truth` (the
start-from-truth window statistic at fixed lags). A section in benchmarks/README.md. Its runs are this slice's
kill-criteria evidence and the design note's tables.

## Verdicts

Against the [MANIFEST](../../benchmarks/baselines/MANIFEST), on the reference build:
- equivalence-1b7d730c: 50 of 55 "identical draws (same RNG stream)". The movers are exactly the drawn-k probit
  scenarios - chik, maskprobit, hazard, hurdle and bart2probit (read from the scenarios and the baseline's k
  channels; probit runs through bartBT at a fixed k). Re-record as equivalence-<sha>, the five recorded fresh and
  merged into a copy, as the 1b7d730c row did; the z sizes the move and is not the oracle. ORACLE (P17): this
  gate's quick and full passes with the never arm's failure as its poison, and the start-from-truth statistic.
- bcf-equivalence-1b7d730c 15 of 15 and multinomial-equivalence-80b1c8d4 11 of 11 bitwise; no re-record.
- The four snapshot files pass unchanged: test-reproducibility-binaryResponse.R fits through bartBT at k = 4.5
  (checked: no k is drawn), so the design's expected regeneration does not arise. A snapshot that moves is a stop.
- Every existing exact gate passes; the probit ones (hazard-, hurdle-, mask-redraw-, categorical-exact) fix k at 2
  and so draw bitwise; hazard-reduction and hurdle-reduction move both sides together and must stay bitwise equal.
- Hot path: bench-sampler.R compare (maintainer-run, quiet machine). run-binary-n1000-p10-t75 takes the move and
  sits at a different k, so it may pass 1.05 from different trees; the move's own cost is criterion 3 below.

## The probit-k SBC arm

No harness change: the arm's dbarts sampler takes the default "auto". The arm keeps R = 600, 99 draws at thin
1000, 30,000 burn (Open call 2). On landing: benchmarks/README.md's "Until the k mixing move lands" paragraph and
the end-bin reading for k go; the probit-k comments in sbc.R ([`sbcConfigProbitK`](../../benchmarks/R/sbc.R) and
the burn default) and in sbc.yaml are rewritten. TODO workflow-text-edits adds `SBC_EXPECTED_FLAGS: k` to the
probit-k job until this lands: if it has landed, this slice deletes that line; if not, its probit-k clause is
struck from that TODO item. Either way this slice edits sbc.yaml, so its push starts the full matrix, which is the
confirmation run (criterion 5).

## Docs and records

- docs/design/probit-k-scale-move.md, from the scratch design: the derivation, the exactness argument, the
  tables, the other families and doors ((c) collapsed, ordinal by scaling the free cutpoints, linear and gp leaves,
  (a) for gaussian-family drawn k), the gate, and the kill criteria with the landed numbers. Its row in
  docs/design/INDEX.md. probit-k-calibration.md's Status points to it; its "today" section is marked as before the
  move.
- docs/architecture.md, One sweep: step 2 gains the move.
- inst/NEWS.Rd, NEW FEATURES, the tree-moves item (scope: against 0.9-34, users' view): "Binary fits under the
  probit link that draw \code{k} - the default for \code{bart}, \code{dbarts} and \code{xbart} - take a further
  step each iteration that rescales the latent variables, the leaf values and \code{k} together, drawn from its
  exact conditional distribution. \code{k} and the fit's overall scale mix far faster, most where the response is
  nearly separated, where a chain could otherwise hold \code{k} for tens of thousands of iterations; the posterior
  is unchanged. \code{dbartsControl(treeScale = )} turns it off." Parse-gated as the README's checklist says.
- TODO: probit-k-mixing removed; dec-B390's entry reads "now that probit-k-scale-move has landed"; a new item:
  "k-mixing-pg-families: measure k's mixing under logistic and nbinom (Polya-Gamma, k drawn by default) by
  benchmarks/R/probit-k-mixing.R's census on a logistic arm and an nbinom arm. (b) does not apply there; if k is
  slow, the collapsed move (c) of docs/design/probit-k-scale-move.md generalizes. Window: Open call 3."
- Ledger entry for the calls made; this plan's Status and Landing; INDEX rows.

## Order of work

1. Branch off origin/bartcore at a green engine tip carrying the small-rulings engine slice.
2. Engine and tests/cpp (plain and under ASAN/UBSan); commit.
3. Bridge, R, help, tinytest; `--preclean`; full tinytest; commit.
4. The gate, exact-gates.yaml, README; quick and full runs, the never arm, the poisons; commit.
5. The measurement script; criteria 2 to 4; commit.
6. Design note, docs, NEWS, TODO, sbc text; commit.
7. Reference build, `--preclean`: the three equivalence compares, the partition, the re-record and its MANIFEST
   row; the four snapshot files; commit.
8. Consumers: stan4bart's, bartCause's, treatSens's and bairrtt's suites against the build; a seeded expectation
   that moves on a drawn-k probit fit is re-recorded on that package's lockstep branch, the orchestrator pushing.

## Gates

On the slice tip, against its own library, independently of the implementer
([RNG classes and their gates](README.md#rng-classes-and-their-gates), posterior-changing): tests/cpp plain and
sanitized; the R-loaded ASAN path over the new and touched test files; the full tinytest suite; the reference-build
compares and snapshots above; every exact-gates.yaml gate in quick mode plus the new one in full; `R CMD check
--as-cran` from a clean tarball; lintr, air, rc-codoc, win-drift, doc-freshness; bench-sampler compare
(maintainer).

Mutants, each must fail a test or the gate: the Jacobian exponent off by two (p = n - nu + 2); the latents not
scaled; k not scaled; the Q term dropped; inactive rows counted in n or their latents scaled; working rebuilt as
alpha (z - o); `totalFits` multiplied in place instead of re-derived (the mapping test's bitwise sum); the move
declining always (the gate's mixing bound, the tinytest "differ" checks); a generator draw taken before an
ineligible forest declines (inertness, bitwise tinytests); the C = 0, n <= nu guard removed (tests/cpp 3).

## Stop conditions

Stop and report, without working around, when (criteria 1 to 4 are the design's, pre-registered):
1. The gate fails with the move on (any |z| > 4.5, or the pure arm's error above 0.01 in full) after one fix.
2. The census on the arm's 100 datasets puts k's tau at a 99th percentile above 2000 or a median above 150
   (prototype 784 and 71).
3. The move's own cost, timed alone on one frozen state at n = 5000 and 200 trees, passes 3 percent of a sweep.
4. The start-from-truth window statistic passes |z| 3 at R = 2000.
5. An equivalence scenario other than the five moves, or any snapshot moves.
6. The diff passes 2600 lines, or one part passes 1.5x its budget line.

After landing, criterion 5 of the design: the first arm run (R = 600, thin 1000) with an end bin above 45 or k's
chi-square p below 0.01 means the arm is not yet ordinary SBC; it is investigated, the reading note restored
meanwhile.

## Interactions

- Engine small-rulings (docs/plans/small-rulings-1008-engine.md, on its own branch): it adds `Chain::scaleDrawnK`
  writing a drawn k at a host call and edits chain.hpp (`SamplerOptions` beside `perturbProbability`, beside
  `setForestFixedK`, the conversion sites) and facade.hpp. Engine slices land serially; it lands FIRST. Its
  verdicts are pinned to the current MANIFEST with no re-record and unchanged snapshots, and they hold only if it
  lands on today's baselines; this slice re-records once, on a tip already carrying it, and its partition is then
  read against one baseline. It is smaller and its plan is written; this one still needs its critique. The hunks
  are textually disjoint and share no code: `scaleDrawnK` keeps the spread at a host call, this move scales k with
  the leaves and latents inside a sweep, and neither calls the other.
- workflow-text-edits: the SBC arm section.
- dec-B390: the binary hyperprior study is rerun after the main merge on chains carrying this move.

## Open calls

1. The switch's name and vocabulary. (a) `treeScale = c("auto", "never")`, word-valued, fixed at creation,
   guarded by setControl - treeShift's idiom, named for what the step does to the trees as treeShift is
   (recommended). (b) `scaleExpansion = TRUE/FALSE`, a logical like useQuantiles, naming the method. (c) A name
   on k: avoid, since `k.scale` already names the leaf prior's reference value (dec-B201). (d) No switch: the move
   is exact and costs nothing, but then no fit can reproduce the draws without it or isolate it in a diagnosis;
   not recommended. Whatever the name, there is no "always": an eligible forest has nothing to skip, an
   ineligible one has no move to force.
2. The arm's thin after landing. Keep 1000 (recommended: with the move every dataset keeps 50 or more effective
   draws of 99, median 98; at thin 300, 8 percent fall under 50, the end-bin excess returning) or 300 (a third of
   the 67-minute job). The TODO's "probably at a smaller thin" was written before the measurement.
3. The window for the logistic and nbinom measurement. Before 1.0-0 (recommended: both draw k by default, the
   same user-facing risk; a census of compute only, whose answer says whether a door is owed before release) or
   after.
