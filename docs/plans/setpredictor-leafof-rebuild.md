# setpredictor-leafof-rebuild

status: CLOSED 2026-08-05 (memo executed under VD's proceed-at-
        discretion grant: all three leafOf-scoped mechanisms declined
        by measurement, ceiling +9.7%; the mu[leafOf]-gather SIMD door
        in the Memo's RECOMMENDATION was TAKEN in a roll-only form,
        LANDED 9141274 2026-08-06 - see the landing note, which also corrects
        the Memo)
        (found 2026-08-04 by the constant-leaf-fits bench discharge,
        docs/plans/constant-leaf-fits.md compare record)
agent: opus
rng: neutral expected (leafOf is a derived cache; draws must not change)
budget: memo first, then ~100-200 engine lines

## Goal

setPredictor-accept returns to within the 5% gate of the pre-compaction
baseline. The constant-leaf-fits refactor made accepted predictor
mutations pay a full O(n * trees) leafOf rebuild (its own landing note
flags this), measured +22-28% on setPredictor-accept-n1000-t75 (x86,
reproducible, isolated to the refactor by an attribution run). This is
the embedded-Gibbs hot path - the package's distinguishing feature - so
the regression matters to per-sweep-mutation consumers (IRT-style
samplers, stan4bart).

## Context

- docs/plans/constant-leaf-fits.md: the leafOf/muByTree design, the
  per-tree leafOfStale lazy-rebuild mechanism (built for
  sampleTreesFromPrior), and the landing note conceding mutation keeps
  the full rebuild.
- The mutation transaction (docs/design/data-store.md): accepted
  setPredictor repartitions changed trees; a tree whose partitions did
  not change needs no leafOf rebuild, and one that changed locally may
  need only a subtree's rows remapped.

## Decision (memo first)

Measure where the +25% goes before building: (a) leafOf rebuild proper,
vs (b) rebuild running eagerly where the old code amortized, vs (c)
double work (rebuild plus a later stale-flagged rebuild). Then pick:
mark-stale-and-rebuild-lazily (reuse leafOfStale; cheapest, but the
next run() pays it), rebuild only trees whose partitions changed, or
per-subtree incremental remap. Evidence that would change the ranking:
if (a) dominates and most trees change on a typical accept, only the
incremental remap helps and the item should be weighed against its
complexity.

## Constraints

- Draws bitwise-unchanged (equivalence trio identical); leafOf equality
  against a freshly derived map stays asserted (the existing 25-sweep
  test).
- No public-surface change.
- Out of scope: the reject path (measured fine), sparse-column
  mutation internals.

## Steps

1. Memo: instrument the accept path on the bench box, attribute the
   +25%, pick the mechanism, record here.
2. Implement per the memo; extend the leafOf-consistency test with an
   accepted-mutation round.
3. Same-machine A/B on dbarts-bench vs 9047e05 and vs the pre-fix tip:
   setPredictor-accept within the gate, run metrics unchanged.

## Verification

- Gate battery per CLAUDE.local.md (install --preclean, tests/cpp from
  clean, full tinytest, equivalence trio bitwise).
- bench-sampler same-machine A/B as step 3; report the three-run table.

## Memo (2026-08-05, step 1 executed on dbarts-bench)

Method: isolated driver (the setPredictor-accept-n1000-t75 scenario
lifted verbatim from bench-sampler.R), pinned taskset -c 2, box idle
(loadavg 0.62-1.31 checked per round), min of 7 reps x 3 rounds;
rdtsc-instrumented builds for composition; md5-over-draws oracle
(20 accepts then run(0,50), train fits + sigma). base = 9047e05,
tip = a34fdff. The isolated driver reproduces +12.2%, not the suite's
+22-28%; the COMPOSITION is measured directly and is what the decision
rests on (if the true suite regression is ~28%, the shares hold).

Attribution of the +0.039 ms/update:

    (a)  leafOf rebuild proper (all 75 trees, every accept)   21%
    (a') mu[leafOf] gather replacing contiguous SIMD passes   79%
    (b)  eager-vs-amortized                                    0
    (c)  double work                                           0

(b): base ran the same two revalidate phases equally eagerly - nothing
was amortized before. (c): installLeafOfAndAddToTotal clears
leafOfStale (chain.hpp ~2677); only sampleTreesFromPrior sets it.
Census: 19/75 trees carry a rule on the mutated column. Cycle-level,
the whole regression is the per-tree fits-rebuild trio: base
subtract 37k + scatter 110k + SIMD add 26k = 173k cycles vs tip
subtract 140k + mu copy 3k + fused install 247k = 390k; of tip's 390k,
~128k is the leafOf write and ~262k is reading tree fits through
mu[leafOf] twice instead of streaming a contiguous row.

Bitwise findings (md5): base == tip (compaction draw-neutral, driver a
valid oracle). tip-noleaf (every leafOf rebuild removed, unfused adds)
== tip - the core move of (ii)/(iii) is bitwise-safe, and its measured
CEILING is +9.7% vs base, still outside the 5% gate. tip-aggr
(subtract/add pair elided for untouched trees) != tip - sigma differs
in the last bits, so (x - v) + v is not identity here and no mechanism
may elide the pair; this kills the only variant that beat base.
Structural constraint: repartitionSubtree(data_, 0) cannot be skipped
even for unaffected trees - misc_partitionRange re-canonicalizes the
index order that leaf-stat sums run over, so skipping changes draws;
only the leafOf write and the fits passes are skippable.

Ranking: (ii) rebuild-only-repartitioned is the best and still
declined - realizable ~+10.3% vs the +9.7% ceiling, recovering ~2 of
the 12 points for ~90-130 lines. (iii) is bounded by (ii)'s ceiling
and applies only to the 19/75 trees (ii) cannot skip; its premise
("most trees change") is false here. (i) mark-stale is invalid before
the economics: deferring forces totalFits to stay on the old map
through the roll - a different fp accumulation of exactly the pair
tip-aggr shows changes draws - and for the mutate-then-run consumer
the deferred rebuild is paid in the very next run() anyway.

RECOMMENDATION: close the item as not worth its complexity. What
would actually close the gap, if funded separately: vectorize the
mu[leafOf] gather in subtractTreeFitsFromTotal / addTreeFitsToTotal
(chain.hpp ~2921-2947) - elementwise and order-preserving, so
bitwise-safe by construction; it attacks the 79%, and the same kernel
serves rollTreeResidual (~2863-2884), the hot sweep path. A misc/SIMD
item with the cross-ISA tests/cpp gate, not this plan's engine lines.
Standing trade either way: at n=1000, t=75 an accepted single-column
mutation costs 0.358 ms against a 0.17 ms sweep; the compaction bought
14-18% on run at n=10000 and gave back 12% here.

## Landing note: the gather door, roll-only (LANDED 9141274 2026-08-06)

The RECOMMENDATION above was designed out and independently critiqued
before building. Both passes agreed on the shape, and both refuted the
premise the RECOMMENDATION was written for. What shipped is a manual
unroll-by-4 (n % 4 prologue first) of three constant-leaf gather loops,
in src/bartcore/chain.hpp and nowhere else:

    G1  rollTreeResidual, t == 0
    G2  rollTreeResidual, t > 0
    G3  finalizeTotalFits

No misc kernel, no dispatch slot, no new ISA translation unit, no
intrinsics. At -O2 the rolled loop compiles scalar on every ISA and no
compiler emits a hardware gather at any flag level (gcc's own cost model
refuses it on znver2, and docs/plans/simd-flag-multiversioning.md already
measured hardware AVX2 gather 1.15-2.2x SLOWER than scalar loads); the
unroll is what lets gcc assemble a 128-bit software gather plus packed
add and lets clang pair the index loads.

Measured, in situ (isolated driver, bench-sampler scenarios verbatim,
pinned core, min of 7-9 reps x 3 rounds, loadavg recorded):

- x86-64 Zen2, gcc 13.3: run-n1000 -3.0% / -3.7% and run-n10000 -3.7% /
  -3.9% in two independent measurements on two different cores, ranges
  disjoint (9 v 9 observations). setPredictor-accept -0.3%, inside the
  +-1% band it was pre-registered to stay in.
- LAYOUT LUCK CLOSED. The gap survives a pure layout perturbation
  (-falign-loops=32 on the engine TU moved the baselines 0.2-0.7% and
  left the gap at -3.5% / -4.3%), and separately a whole-package -mavx
  rebuild. The 3-4% is the unroll.
- Whole-package -mavx buys ZERO end-to-end on either base or roll
  (differences inside noise, both sizes). The dispatched-misc-kernel
  design's only differentiator is worth nothing; it was rejected.
- arm64 (Apple Silicon, clang 17): run-n1000 -5.9%, run-n10000 -12.7
  to -13.9% (quiet-box A/B, 3 rounds, base and roll back-to-back per
  scenario, min of 9, 1-min loadavg 1.75-1.88 throughout; base minima
  identical across all rounds). setPredictor-accept -0.8 to -1.6%,
  crossing the pre-registered +-1% band in the FAST direction;
  reattributed per protocol by call-site census: rollTreeResidual and
  finalizeTotalFits are called only from the sweep loop and
  growForestFromRoot, never from the setPredictor path, so the
  movement is code layout - the class the x86 -falign-loops probe
  bounded at 0.2-0.7% - not execution of the new code.
- Draws bitwise identical at every variant and width measured
  (md5 over train fits and sigma, n mod 4 in {0,1,2,3}, gaussian and
  probit and storage = "single"), and the landing gates below reconfirm.

fp32 IS GUARDED OUT. Each unrolled body sits behind
`if constexpr (std::is_same_v<ResidT, double>)` with the rolled loop kept
verbatim in the else. clang already vectorizes the ResidT = float form of
G1/G2 at tip with NEON lane-insert gathers (`ld1.d` + `fadd.2d`) - the
very mechanism the design memo had listed as REJECTED-but-plausible - so
a plain unroll would bolt a scalar prologue onto codegen that is already
good and inflate G3<float> for nothing. Verified at the object level: the
float instantiation's emitted code is instruction-identical pre and post
patch (`Chain<ConstantGaussianLeaf, float>::run`, 800 instructions, and
all 161 float-instantiated symbols in the engine TU).

Corrections to the Memo above, all measured:

1. "reading tree fits through mu[leafOf] twice" is wrong. The
   constant-leaf accept path reads mu[leafOf] ONCE - subtractTreeFitsFromTotal
   (G5, 140k cycles) - plus installLeafOfAndAddToTotal's SCATTER add
   (~119k). Only the 140k is addressable by a gather; 140/262 = 53%. The
   "79% of the regression is the gather" share double-counts by ~2x.
2. addTreeFitsToTotal (G4) is NOT on the constant-leaf accept path at
   all: chain.hpp puts installLeafOfAndAddToTotal in the
   `if constexpr (leafIsConstant)` arm and G4 in the `else`. The Memo
   names it as a primary target; for that path it is dead code.
3. Unrolling G4/G5 - the two sites the Memo names - measured 1.7-2.4%
   WORSE on setPredictor-accept, reproducible, ranges disjoint, in two
   independent runs. DECLINED. Do not re-propose. installLeafOfAndAddToTotal
   is a scatter with no portable SIMD form below AVX512, and the Memo's
   own +9.7% ceiling already keeps that item closed.
4. G2 - the hottest site, 74 of 75 trees per sweep - only vectorizes its
   ADD on x86. All four body subtractions stay scalar `subsd`. G1 and G3
   vectorize fully. The mechanism is half operative where most of the
   traffic is; the win is still real.
5. "A misc/SIMD item with the cross-ISA tests/cpp gate" - there is NO
   cross-ISA gate in tests/cpp. main.cpp calls misc_simd_init() once and
   never re-dispatches. The dispatch gate is inst/tinytest/test-simd.R,
   at the R level - and it CANNOT fail on a header change like this one,
   because it compares dispatch levels of a single binary and a header
   change moves every level identically. The real regression nets here
   are the value-hardcoded tinytests and the equivalence trio.
6. What the Memo got right and carried the item: the gather is
   elementwise and order-preserving, which makes the change bitwise-safe
   by construction, and the same code serves rollTreeResidual. That
   second half turned out to be the entire value.

Line references at the time of landing: rollTreeResidual and
finalizeTotalFits carry the WHY comments; G4/G5/S1 are untouched.
New coverage: tests/cpp testGatherTailShapes drives n = 1001 (a
one-element prologue) for both ResidT instantiations and pins
totalFits against a freshly gathered per-tree sum and the rolled
residual against a scalar one - the n % 4 prologue is the only new
failure mode, and a deliberately broken prologue was confirmed to
fail both assertions.
