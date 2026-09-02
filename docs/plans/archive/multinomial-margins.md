# multinomial-margins

agent: Opus
rng: shifting, multinomial channels ONLY (rounded psi changes the PG
  stream); every non-multinomial anchor stays bitwise
budget: ~200 changed lines (stop past ~300)

## Goal

Reduce the per-sweep multinomial margin cost from O(nK^2) to O(nK):
otherCategoryMargin recomputes a full K-term log-sum-exp per
observation per category; replace it with a per-sweep suffix table and
a running prefix merge, two-term LSE per (i, f).

## Context

- drawForestGlue(f) ([[combiner.hpp:645-659@23fff839]]) forms C_if =
  log sum_{j != f} exp(f_ij) via otherCategoryMargin
  ([[combiner.hpp:775-784@23fff839]]), O(K) exps per call -> O(nK^2) per sweep.
- The interleaving contract (comment at [[combiner.hpp:638-644@23fff839]]): when
  category f draws, forests[j].totalFits for j < f hold THIS sweep's
  fits, j > f hold LAST sweep's. Any restructure must reproduce
  exactly that mix - a whole-sweep precompute over one snapshot is
  wrong.
- Shape: at sweep start (f == 0) build suffix_[f][i] = LSE over j > f
  of the OLD fits (one O(nK) backward pass, n x K storage, or the
  equivalent rolling form); maintain prefix_[i] = LSE over j < f of
  the NEW fits by one two-way merge per category after forest f's
  update (O(n) each). Then C_if = LSE2(prefix_[i], suffix_[f][i]).
  The combiner has no post-forest-update hook per category except the
  next drawForestGlue call - folding the prefix merge of category
  f - 1's new fits into the top of drawForestGlue(f) is acceptable.
- combinedFits/combinedTestFits and the exact-gate blend are NOT in
  scope; they are already O(nK) through softmaxLocationMajor.

## Constraints

- Two-way LSE merges only: LSE2(a, b) = max + log1p(exp(-|a - b|)) or
  the max + log(sum exp) form - implementer's choice, but NEVER
  subtract-exp from a full-K LSE (catastrophic cancellation when the
  excluded category is the argmax; this is the recorded reason the
  naive fix is wrong).
- Empty-set identities must be exact: prefix at f == 0 and suffix at
  f == K - 1 are LSE over nothing (-inf); LSE2(-inf, x) must return
  exactly x (branch, do not compute through exp(-inf)). K == 2 then
  reduces C_if to the other category's raw fit - exactly what the old
  code produced.
- Values may differ from the old path only by rounding; the margins
  are the same mathematical quantity. Expect gate divergence CONFINED
  to multinomial: equivalence-ac6ec2c and bcf-equivalence MUST stay
  bitwise (any divergence there = STOP, report, revert nothing
  yourself). The multinomial-equivalence compare is EXPECTED to fail;
  report which channels moved (all should).
- benchmarks/R/multinomial-exact.R must PASS (statistical gate; the
  new margins are the same posterior). Run it. If tinytest has
  RNG-locked multinomial snapshots that now fail, do NOT regenerate
  them - report the failing files; snapshot replay is the
  orchestrator's landing step.
- No dbarts.h or inst/include change; no R change.
- Do not commit; leave the diff in the tree.

## Steps

1. Storage: suffix table (or rolling equivalent) + prefix vector on
   the combiner; sized once, reused per sweep.
2. Rewrite drawForestGlue's margin source; delete or demote
   otherCategoryMargin (keep it only if a cold path still needs it).
3. Verify the interleaving mix (new below f, old above) is preserved;
   the f-specific handoff to formForestResponse is unchanged.

## Verification

R CMD INSTALL --preclean .; tests/cpp make clean && make &&
./test_bartcore; tinytest (multinomial snapshot failures reported,
not regenerated); equivalence.R vs equivalence-ac6ec2c.rds 22/22
BITWISE (hard requirement); bcf-equivalence.R vs
bcf-equivalence-99205ee.rds BITWISE (hard requirement);
multinomial-equivalence.R vs multinomial-equivalence-2bd34db.rds
EXPECTED to diverge - report channels; multinomial-exact.R PASS.
Re-record of the multinomial fixture, MANIFEST update, and any
snapshot replay happen at landing, orchestrator-run.

## Landings

8c2b5fc (2026-07-17). Landed as planned: per-sweep suffix snapshot
(backward two-way merges from an empty top row) + running prefix
folded at the top of each in-order drawForestGlue, margin =
LSE2(prefix, suffix); logSumExp2 branches -HUGE_VAL exactly and
otherCategoryMargin is deleted. One addition beyond the sketch: a
lastF_ cursor detects fresh-sweep entry (or a direct out-of-order
call) and rebuilds; on the sampler path only f == 0 rebuilds, so the
call pattern is unchanged. The tests/cpp single-trial reduction
oracle re-inlined the OLD full-K arithmetic; its margin was updated
(by the reviewer, at landing) to the engine's pairwise form -
symmetric LSE2 makes the K = 3 static-fits margin one merge of the
two other fits - keeping the reduction assertion (one PG draw, same
psi, indicator response) bitwise-meaningful. Gates: 22/22 and bcf 5x6
BITWISE (shift confined as required); multinomial-exact all arms OK
(gaps 0.0000/0.0008/0.0012 vs tol 0.008/0.008/0.015); suite 3050/0
(no R-level multinomial snapshots exist); component tests pass after
the oracle update. Fixture re-recorded as
multinomial-equivalence-8c2b5fc.rds with the neutrality trail in the
MANIFEST: k2 bitwise-identical to 2bd34db on all channels, k3/
k3counts moved in train/test/forestFits only with both varcount
channels unmoved. Diff 49/12 + the oracle update.
