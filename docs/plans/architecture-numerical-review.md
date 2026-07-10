# architecture-numerical-review

agent: two opus readers (architecture critic; numerical auditor),
       read-only; fable adjudicates
rng: neutral (findings only; any fix becomes its own gated item)
budget: two reports recorded here; findings become fix items,
        design-note updates, or explicit accepted-risk notes

Review 6 of the retrospective program (retrospective-reviews.md),
run before review 5 (the roadmap survey needs sustained web access;
this one is purely local).

## Goal

Fresh-eyes stress of the engine now that it is functionally frozen:
(a) does the facade/concept decomposition actually carry the next
five features without deforming, and (b) where does the arithmetic
break first under extreme but reachable inputs - scale, sample
count, weight ratios, and long single chains.

## Method

Two independent read-only Opus readers, one dimension each.

A. ARCHITECTURE CRITIC: stress the decomposition (facade.hpp
   concept selection; Chain/Forest split; ResponseModel decoration;
   LeafModel concepts constant/vector/function; the mutation
   surface; the C bridge) against the five likely-next features,
   each already specified in TODO/design notes:
   1. grow-from-root-warm-start (a second tree-construction kernel
      beside metropolisJumpForTree, sharing a cut-scan primitive);
   2. multi-forest-models (multinomial/heteroscedastic/hurdle on
      the BCF two-forest template);
   3. negative-binomial + ordinal-outcomes (new response families:
      PG-with-r-update, cumulative probit with a cutpoint block);
   4. data-ownership (owned predictor storage replacing borrowed
      REAL(x), per-column widths);
   5. the parallel seam (blocked-jacobi trees / gpu-bart backend).
   For each: what bends, what breaks, what is already provisioned.
   Verdict per feature: FITS (seams named) / STRAINS (what deforms
   and the smallest re-cut) / BLOCKED (what must change first).
   Known couplings to weigh: BCF is constant-leaf-only by
   static_assert; grouped and BCF are mutually exclusive by
   construction; the empty-leaf veto is the single surrogate target
   all moves share; the proposal-equals-prior identity underpins
   birth/death cancellation.

B. NUMERICAL AUDITOR: where does the arithmetic degrade first.
   Seed list (verify, quantify the onset threshold, and judge
   severity; add what you find beyond it):
   - the constant-leaf centered sum-of-squares catastrophic
     cancellation under a large shared residual offset (block-5
     fragility note; score side only);
   - totalFits telescoped rounding drift in a pure run (~1e-12
     relative at 1e5 sweeps measured; single-chain n >= 1e5 is a
     standing workload - where does it matter first);
   - the GP path: 1e-6 conditioning jitter vs the nugget; the
     predict-vs-stored-fits discrepancy (~2e-3 at training rows -
     re-kriging through the cached factorization is not the
     recorded fit; review-4 finding) - defect, doc note, or fix;
   - the grouped tau slice sampler: hangs under extreme ranef
     scale (review-3 poison side effect; review 4 found the
     half-Cauchy tail makes prior draws stall it) - bound the
     stall region, propose the cheapest guard (step-out cap,
     bracket limit);
   - the 1e-9 BCF multiplier floor (block-6 note: magic constant
     unrelated to data scale);
   - log-space robustness of the acceptance ratios (change move's
     counted-ratio logs, chi-k gamma draws at extreme shapes -
     review 3 saw k run away to 1e153 under a poison; is the
     UNPOISONED k feedback loop bounded for extreme nu/M);
   - range scaling at extreme response scales (y ~ 1e15 or 1e-15:
     where do fitScale/fitShift/sigma round-trips lose precision);
   - setState occupancy validation accepting an n-mismatched
     gaussian state (block-7 note) - quantify what a mismatch
     corrupts.
   For each: REAL (onset threshold + severity + smallest fix),
   ACCEPTED RISK (why it cannot bite at reachable inputs), or
   DOC NOTE. Small numerical experiments in R against the
   installed package are encouraged; no engine edits.

## Constraints

- Findings only; the tree does not change under this review.
- Readers are read-only, FOREGROUND only, and spawn nothing.
- docs/design/*.md may be read for INTENT here (unlike the
  correctness audit, this review judges design against goals, so
  the design notes are evidence, not contamination).

## Verification

- Both reports recorded here; each REAL finding becomes a fix item
  or an accepted-risk note with VD visibility.

## Status

- 2026-07-09: plan authored; both readers dispatched.
