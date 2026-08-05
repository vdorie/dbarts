# setpredictor-leafof-rebuild

status: MEMO RECORDED 2026-08-05 (see Memo below) - the plan's three
        mechanisms are all declined by measurement; decision pending
        VD: close, or scope the mu[leafOf] gather SIMD item instead
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
