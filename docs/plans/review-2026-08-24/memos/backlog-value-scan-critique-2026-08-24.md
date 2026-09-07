# Adversarial critique: backlog-value-scan-2026-08-24 (branch bartcore, ee70e6f7)

Posture REFUTE. `git archive HEAD` staged under a private scratch prefix, `R CMD INSTALL
--preclean -l ../lib .` (EXIT 0); every probe `R_LIBS=lib Rscript`. Main checkout untouched.

## 1. Asserted defects
### (a) Heteroscedastic loglik and PPD - CONFIRMED, memo UNDERSTATES scope

`bart2(x, y, variance = 2L, ...)`, n=60, 20 draws (p1-hetero.R, p1b-ppd.R):

    family: gaussian  len(sigma): 20  unique(sigma): 1 = 7.806093 == diff(range(y))
    s.train range: 0.355949 3.94287
    max |loglik - dnorm(y, ev, scalar sigma)|: 0     |  - dnorm(y, ev, s.train)|: 2.894689
    sum loglik as reported: -3592.142    sum loglik correct: -2031.674
    sd(ppd - ev): 8.043541   scalar sigma: 7.806093   mean s.train: 1.402687
    cor(colwise sd(ppd-ev), colwise mean s.train): 0.0206

sigma is the fixed unit residual of parameterisation A ([docs/design/heteroscedastic.md:279-281](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/docs/design/heteroscedastic.md#L279-L281), `sigmaIsFixed_`)
rescaled by the y range - no posterior content at all. TWO CHANNELS THE MEMO MISSES, both live:

- `predict(type = "ppd")` ([R/generics.R:270-296](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/R/generics.R#L270-L296)) computes `s` into a local at [R/generics.R:273-279](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/R/generics.R#L273-L279), attaches it as
  attr "s" to the very object it returns, and still calls `sampleFromPPD` with the pinned scalar.
  p6-predppd.R:

      s at low-x2 rows: 1.277   s at high-x2 rows: 5.406
      sd(ppd - ev) low-x2: 21.38  high-x2: 23.16  | scalar sigma: 22.18

- `summary()` on a heteroscedastic fit (p8-scan2.R) prints exactly one row - the placeholder
  presented as the residual scale, nothing about s(x):

      variable  mean median  sd  mad   q5   q95  rhat
      sigma     22.2   22.2   0    0  22.2  22.2    NA

Hetero + `weights` is also ACCEPTED, loglik scoring `dnorm(y, ev, scalarSigma/sqrt(w))` exactly
(p11.R) with `mean(s.train) = 4.16` against `sigma = 19.6`. VERDICT CONFIRMED, scope 4 exported
channels not 2; `fitted`/`residuals`/`print` are clean. Not new to the repo - it is feature-matrix
[f28] ([R/generics.R:802-805](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/R/generics.R#L802-L805)), which the memo's sec 3 concedes.

### (b) Variance-forest setState column-mask gap - CONFIRMED, with a reproducer
p2c.R / p2d.R: donor `variance = TRUE` (both columns), recipient `variance = "v2"`, same data and
shape.

    donor variance split-variable table:     -1:100   1:32   2:28
    recipient own table:                     -1:86           2:46
    rec$setState(donor state)  ->  ACCEPTED
    post-setState recipient table:           -1:100   1:32   2:28
    after 5 further sweeps:                  -1:97    1:12   2:45
    rec$installTrees(donor)    ->  REFUSED: "warm-start donor holds a tree that splits on a
      variable outside this forest's allowed column set..."

32 forbidden splits install and 12 are still live five sweeps later;
`splitVariableLogProbability` ([src/bartcore/model.hpp:2135-2148](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/src/bartcore/model.hpp#L2135-L2148)) scores them against `collectAvailableVariables`,
which excludes the masked column, so s(x) depends on a column the declaration forbids. Cite drift
only: the false comment is [src/bartcore/chain.hpp:3345-3346](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/src/bartcore/chain.hpp#L3345-L3346), the variance pass [src/bartcore/chain.hpp:3384-3403](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/src/bartcore/chain.hpp#L3384-L3403).

### (c) dbartsData(bases =) - CONFIRMED, and worse than "ambiguous"
p3-bases.R. The two exported entry points accept precisely what the other refuses.

    Case A, N=10, subset = 1:6
      dbartsData(formula, subset, bases = 10-row)      REFUSED ("must have the same length as 'y'")
      dbartsData(formula, subset, bases = 6-row)       ACCEPTED
      dbarts(formula, subset, forest(basis = 10-row))  ACCEPTED
      dbarts(formula, subset, forest(basis = 6-row))   REFUSED (matches 'subset', not the full data)
    Case B, N=10, subset = rep(1:5, 2)  (equal counts)
      dbartsData basis rows: 1,2,3,4,5,6,7,8,9,10
      dbarts()  basis rows: 1,2,3,4,5,1,2,3,4,5

Case B builds different data from the same arguments, no diagnostic. Root cause [R/data.R:1067](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/R/data.R#L1067)
(post-subset count) vs [R/data.R:1165](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/R/data.R#L1165) (full length + subset), with [R/model.R:962-964](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/R/model.R#L962-L964) taking the
`n == full` branch first when full == length(index).

### (d) Multinomial twin-create - LIVE, but the SURFACE claim is REFUTED

Live at HEAD: [R/bart.R:1790](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/R/bart.R#L1790) + 1802 (ordinal), same pair at [R/bart.R:2056](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/R/bart.R#L2056) (nbinom), `adoptPointer` both. The
memo's "SURFACE: draw-moving at a fixed seed for two shipped families" is false. `createChainRngs`
([src/R_interface_bartcore.cpp:1864-1891](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/src/R_interface_bartcore.cpp#L1864-L1891)): with a seed the chain seeds come from a dedicated MT; without
one they are `unif_rand()` draws off R's stream. p4-twin.R:

    fixed seed=77, different R streams identical:     TRUE
    seed=NA, same R stream identical:                 TRUE
    seed=NA, stream advanced by one draw identical:   FALSE

[docs/design/multinomial-mutation-arc.md:389-398](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/docs/design/multinomial-mutation-arc.md#L389-L398) already says exactly this and already ruled: "MOVES DRAWS - but
only for `seed = NA` ... the 2-of-43 re-record is a property of `fitSummaries`' `set.seed()`
fixture, not of the change ... Do not spend the program's first re-record to reclaim 3ms. Defer."
The memo answers its own "which two, why" wrongly (the harness fixture is unseeded, no seeded user
moves) and relitigates a settled deferral with "now or never" and no new evidence. Mechanism
CONFIRMED, surface claim REFUTED, strike the recommendation.

## 2. PRE-RC verdicts 6-9
Threading cost across four landed formals (R/ hits / man files / test files): `monotone` 47/8/12,
`blocks` 57/10/13, `interactions` 52/7/5, `variance` 82/13/26. Nearest landed "open a family's
public surface" slice, `git show --stat 5a3bc276`: R/ +723, one test file +537, five man files.

6. group-by-exposure. Value clause SURVIVES: "a multi-site binary or count outcome cannot get a
   random intercept at all, rbart_vi being gaussian/aft only". Engine risk LOWER than the memo
   implies - `GroupedResponse` is a decorator applied unconditionally at [src/bartcore/chain.hpp:625-629](https://github.com/vdorie/dbarts/blob/5a3bc276693514ca8e59956e5d174df58d11aa49/src/bartcore/chain.hpp#L625-L629), after
   every family branch, so all six families are wrapped-capable. A NEW defaulted formal =
   ADDITIVE-LATER. ~200 R is optimistic 2-3x against 5a3bc276 and the threading counts, and the tau
   `rel.scale` for latent families is unadjudicated. KEEP as face, not lock; reprice.
7. survival `entry=`. Value clause SURVIVES ("staggered-entry registry/EHR cohorts are hand-expanded
   today"). New defaulted formal = ADDITIVE-LATER. ~60 R optimistic against 47-82 sites plus a
   family gate and `warnFamilyGatedArgs`; reprice ~150.
8. sparse-extensions halves. Refusals verified live ([R/rbart.R:302-307](https://github.com/vdorie/dbarts/blob/5a3bc276693514ca8e59956e5d174df58d11aa49/R/rbart.R#L302-L307), [R/model.R:215-220](https://github.com/vdorie/dbarts/blob/5a3bc276693514ca8e59956e5d174df58d11aa49/R/model.R#L215-L220) and [R/model.R:255-265](https://github.com/vdorie/dbarts/blob/5a3bc276693514ca8e59956e5d174df58d11aa49/R/model.R#L255-L265),
   engine reason [src/bartcore/data.hpp:287-292](https://github.com/vdorie/dbarts/blob/5a3bc276693514ca8e59956e5d174df58d11aa49/src/bartcore/data.hpp#L287-L292)). Removing a refusal is ADDITIVE-LATER by construction - the memo
   concedes this in its own last clause. NOT pre-RC on a window argument.
9. rbart_vi logistic token. "Grouped binary silently gets probit" is QUALIFIED - p10.R:
   `rbart_vi(family = "logistic")` is REFUSED by name ("'arg' should be one of auto, gaussian,
   aft"), and `rbart_vi(binary y, "auto") -> family: probit` is the documented default, not a
   silent downgrade. The real value clause is narrower and still good: logistic is the only binary
   family accepting weights ([R/spec.R:52-67](https://github.com/vdorie/dbarts/blob/5a3bc276693514ca8e59956e5d174df58d11aa49/R/spec.R#L52-L67)), so grouped WEIGHTED binary is unreachable. Widened
   formal value = ADDITIVE-LATER. Cheapest of the four.

Net: none of 6-9 is window-locked in the sense "adding it after 1.0 breaks callers" - all four are
clean additive formals or widenings, pre-RC only on the weaker "the 1.0 face should be symmetric"
argument. The memo should say so rather than let them outrank a wrong number.

## 3. Doors
- D1 per-forest row subsetting: MISSTATED. [TODO:221-223](https://github.com/vdorie/dbarts/blob/5a3bc276693514ca8e59956e5d174df58d11aa49/TODO#L221-L223) says the classes may want different SHAPES,
  not grids. The source fork is mask / ZERO WEIGHT / compaction; zero weight already SHIPPED
  ([docs/design/model-space-survey.md:307-312](https://github.com/vdorie/dbarts/blob/5a3bc276693514ca8e59956e5d174df58d11aa49/docs/design/model-space-survey.md#L307-L312)) and the per-forest MASK was permanently REFUSED as redundant with
  `setForestWeights` ([docs/plans/latent-subset-mask.md:323-326](https://github.com/vdorie/dbarts/blob/5a3bc276693514ca8e59956e5d174df58d11aa49/docs/plans/latent-subset-mask.md#L323-L326), [docs/plans/latent-subset-mask.md:488-492](https://github.com/vdorie/dbarts/blob/5a3bc276693514ca8e59956e5d174df58d11aa49/docs/plans/latent-subset-mask.md#L488-L492), arc COMPLETE). The memo drops the
  shipped arm and conflates mask with zero weight. Honest residue: PHYSICAL COMPACTION only, already
  deferred by VD with a named trigger ([TODO:223-224](https://github.com/vdorie/dbarts/blob/5a3bc276693514ca8e59956e5d174df58d11aa49/TODO#L223-L224)). As posed it reopens two settled arms.
- D2 approximate-Polya-Gamma: correct and genuinely open ([docs/design/negative-binomial.md:105-108](https://github.com/vdorie/dbarts/blob/5a3bc276693514ca8e59956e5d174df58d11aa49/docs/design/negative-binomial.md#L105-L108) resolved fork
  A only; [docs/design/negative-binomial.md:646-654](https://github.com/vdorie/dbarts/blob/5a3bc276693514ca8e59956e5d174df58d11aa49/docs/design/negative-binomial.md#L646-L654) holds fork B). Two fixes: "explicit opt-in" is the TODO's gloss, never
  adjudicated - the record says a "project-level identity decision ... its own arc" ([docs/design/negative-binomial.md:309-313](https://github.com/vdorie/dbarts/blob/5a3bc276693514ca8e59956e5d174df58d11aa49/docs/design/negative-binomial.md#L309-L313),
  [docs/design/negative-binomial.md:648-651](https://github.com/vdorie/dbarts/blob/5a3bc276693514ca8e59956e5d174df58d11aa49/docs/design/negative-binomial.md#L648-L651)); and [TODO:498](https://github.com/vdorie/dbarts/blob/5a3bc276693514ca8e59956e5d174df58d11aa49/TODO#L498)'s FIRST half (integer-weight probit by replicated latents) is exact and
  cost-gated, not part of this door ([docs/plans/weighted-binary.md:18-32](https://github.com/vdorie/dbarts/blob/5a3bc276693514ca8e59956e5d174df58d11aa49/docs/plans/weighted-binary.md#L18-L32)).
- D3 hurdle samplerOnly: correct and open ([docs/design/multinomial-mutation-arc.md:956-961](https://github.com/vdorie/dbarts/blob/5a3bc276693514ca8e59956e5d174df58d11aa49/docs/design/multinomial-mutation-arc.md#L956-L961) recommends defer;
  hurdle.md carries no ruling). Live-code fix: THREE families opt back in via `allow.samplerOnly` -
  multinomial ([R/bart.R:919](https://github.com/vdorie/dbarts/blob/5a3bc276693514ca8e59956e5d174df58d11aa49/R/bart.R#L919)), ordinal ([R/bart.R:1072](https://github.com/vdorie/dbarts/blob/5a3bc276693514ca8e59956e5d174df58d11aa49/R/bart.R#L1072)), nbinom ([R/bart.R:1105](https://github.com/vdorie/dbarts/blob/5a3bc276693514ca8e59956e5d174df58d11aa49/R/bart.R#L1105)); the comment at [R/bart.R:603](https://github.com/vdorie/dbarts/blob/5a3bc276693514ca8e59956e5d174df58d11aa49/R/bart.R#L603) saying
  "ordinal, nbinom" went stale at S2+S3, a small independent defect.
- BCF setData smuggling: NO. D1 is fixed-n membership, Door 1 is n-free joint x/y/z, and the survey
  lists the mixture-cure AFT BCF among the fixed-n classes. But the memo's PROOF is anchored wrong:
  [src/C_interface.cpp:665-673](https://github.com/vdorie/dbarts/blob/5a3bc276693514ca8e59956e5d174df58d11aa49/src/C_interface.cpp#L665-L673) and [src/C_interface.cpp:694-698](https://github.com/vdorie/dbarts/blob/5a3bc276693514ca8e59956e5d174df58d11aa49/src/C_interface.cpp#L694-L698) are the setPredictor/updatePredictor guard-retirement
  comments and `grep -c setData src/C_interface.cpp` = 0 - there is no setData in the flat C API at
  all. Re-anchor to [src/R_interface_bartcore.cpp:4631](https://github.com/vdorie/dbarts/blob/5a3bc276693514ca8e59956e5d174df58d11aa49/src/R_interface_bartcore.cpp#L4631).

## 4. Stale-premise findings
- dimnames ([TODO:265-267](https://github.com/vdorie/dbarts/blob/5a3bc276693514ca8e59956e5d174df58d11aa49/TODO#L265-L267)): substance CONFIRMED (e8fcf8be; 211 commits before HEAD;
  [R/bartcore.R:293-307](https://github.com/vdorie/dbarts/blob/e8fcf8be106fe14e33944e79342e60285c3baa8e/R/bartcore.R#L293-L307), [man/dbartsSampler-class.Rd:373](https://github.com/vdorie/dbarts/blob/e8fcf8be106fe14e33944e79342e60285c3baa8e/man/dbartsSampler-class.Rd#L373), [inst/tinytest/test-sampler-predictors.R:75-105](https://github.com/vdorie/dbarts/blob/e8fcf8be106fe14e33944e79342e60285c3baa8e/inst/tinytest/test-sampler-predictors.R#L75-L105)). CHRONOLOGY
  REFUTED: `git blame` dates the "fix is UNTICKETED" clause to 96b3f0f18, 2026-08-19, THREE DAYS
  AFTER the 2026-08-16 fix - not "the pre-e8fcf8be world" but a post-fix sweep that saw the missing
  ticket and not the fix. Sparse/container residual CONFIRMED ([R/bartcore.R:331](https://github.com/vdorie/dbarts/blob/e8fcf8be106fe14e33944e79342e60285c3baa8e/R/bartcore.R#L331) on; mixedMatrix.R
  never reads the argument's names).
- feature-matrix [f19] / [R/bartcore.R:1039](https://github.com/vdorie/dbarts/blob/e8fcf8be106fe14e33944e79342e60285c3baa8e/R/bartcore.R#L1039): CONFIRMED (4975c20b, c3af16a1 real and as described;
  [R/generics.R:66-93](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/R/generics.R#L66-L93) guards and scores the t marginal). Nit: [f19]'s "per-row lambda_i is not stored"
  clause is still true, so it is stale in its verdict only.
- INDEX.md nameable-calibration row: CONFIRMED ([docs/design/INDEX.md:82](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/docs/design/INDEX.md#L82) vs [docs/design/feature-matrix.md:1108-1112](https://github.com/vdorie/dbarts/blob/7a42305e17bdb90d6e82a7fd6f5e3e2b7622ef30/docs/design/feature-matrix.md#L1108-L1112); ab3aa2fa
  real). MISSED: [docs/design/INDEX.md:73](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/docs/design/INDEX.md#L73) still calls it one of "the two designed-but-unbuilt arcs".
- bcf-bartcause-relocation anchors: MIXED. Path is docs/PLANS/, not docs/design/. Stale [docs/plans/bcf-bartcause-relocation.md:80-143](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/docs/plans/bcf-bartcause-relocation.md#L80-L143)
  and [docs/plans/bcf-bartcause-relocation.md:497-500](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/docs/plans/bcf-bartcause-relocation.md#L497-L500) cites CONFIRMED; vanished bart2-argument-consolidation entry CONFIRMED (a dangling
  ordering reference survives at [docs/plans/bcf-bartcause-relocation.md:405](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/docs/plans/bcf-bartcause-relocation.md#L405)); S13 172523e6 door stale CONFIRMED. "Two of its doors are
  answered too" REFUTED - the per-forest replay door at [docs/plans/bcf-bartcause-relocation.md:1271-1274](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/docs/plans/bcf-bartcause-relocation.md#L1271-L1274) ALREADY records "LANDED
  2026-08-20 at dbarts 63df524e", landing note at [docs/plans/bcf-bartcause-relocation.md:1613](https://github.com/vdorie/dbarts/blob/0047460641d0ff14bc0dd336ef3788022334c07b/docs/plans/bcf-bartcause-relocation.md#L1613). One stale door, not two. [R/generics.R:505](https://github.com/vdorie/dbarts/blob/63df524e8b54948c6bdca1820a949f2c5eb17385/R/generics.R#L505) is
  the comment above `predictForest` ([R/generics.R:516](https://github.com/vdorie/dbarts/blob/63df524e8b54948c6bdca1820a949f2c5eb17385/R/generics.R#L516)).

## 5. New findings
- DEFECT x2: `predict(type = "ppd")` and `summary()` on a heteroscedastic fit - sec 1a.
- SURFACE: `summary()` has methods for bart, bartMultinomial and rbart only. On ordinal, nbinom or
  hurdle it falls through to `summary.default` and prints a str-style Length/Class/Mode table
  (p9-fams.R). Additive-later, but an ugly 1.0 face. Not in TODO, not in the memo.
- MINOR: `xbart(..., variance = )` dies with R's raw "unused argument (variance = 2)", the only gate
  found that does not follow error-style.md.
- MINOR: the memo's header says "1 PROCESS" while sec 1 marks two ([R/generics.R:372](https://github.com/vdorie/dbarts/blob/63df524e8b54948c6bdca1820a949f2c5eb17385/R/generics.R#L372) and the release block
  [R/generics.R:510](https://github.com/vdorie/dbarts/blob/63df524e8b54948c6bdca1820a949f2c5eb17385/R/generics.R#L510)). TODO holds 41 top-level entries, all covered, so nothing is missed - only the count line.
- Checked and REFUTED as candidates: aft + weights refused outright; `extract(type = "loglik")` on
  nbinom/ordinal/hurdle refused by `match.arg`, not silently scored; `fitted`/`residuals`/`print` on
  a heteroscedastic fit correct.

## 6. Corrected PRE-RC list
1. Heteroscedastic scale channels. Score and draw at s(x)/sqrt(w) in `pointwiseLogLikelihood`,
   `sampleFromPPD` (extract AND predict call sites) and `summary`; refusal as fallback. SURFACE four
   exported channels, NEWS-visible. VALUE loo/waic, PPD checks and the first summary a user prints
   stop being wrong by 4-6x. COST ~120 R + tests. TIER Sonnet, Opus read on the draws axis.
2. Variance-forest setState column-mask gap. A variance pass in `stateIsValid` plus the
   `rebuildVarianceForest` backstop the mean side already has. SURFACE refusal on a shipped method.
   VALUE correct s(x) after any restore; unblocks the SBC hetero setState lift. ~40 engine, Opus.
3. dbartsData(bases =) alignment. Converge dbartsData's formula path on the dbarts() full-data rule;
   refuse the ambiguous shape by name at BOTH entries. SURFACE an exported constructor's contract,
   today disjoint between entries. COST ~25 R + tests. Sonnet.
4. ordinal-scan-missing-rows ([TODO:305](https://github.com/vdorie/dbarts/blob/63df524e8b54948c6bdca1820a949f2c5eb17385/TODO#L305)). Independently CONFIRMED: [src/bartcore/scan.hpp:100](https://github.com/vdorie/dbarts/blob/63df524e8b54948c6bdca1820a949f2c5eb17385/src/bartcore/scan.hpp#L100) skips `naCode` and
   totals surviving bins ([src/bartcore/scan.hpp:105-108](https://github.com/vdorie/dbarts/blob/63df524e8b54948c6bdca1820a949f2c5eb17385/src/bartcore/scan.hpp#L105-L108)) while [src/bartcore/grow.hpp:181-186](https://github.com/vdorie/dbarts/blob/63df524e8b54948c6bdca1820a949f2c5eb17385/src/bartcore/grow.hpp#L181-L186) builds the no-split term from the FULL
   node stats. Sole caller of `scanOrdinalCuts` is [src/bartcore/grow.hpp:268](https://github.com/vdorie/dbarts/blob/63df524e8b54948c6bdca1820a949f2c5eb17385/src/bartcore/grow.hpp#L268), so reach is `n.grow.sweeps > 0`
   (default 0L, [R/bart.R:677](https://github.com/vdorie/dbarts/blob/63df524e8b54948c6bdca1820a949f2c5eb17385/R/bart.R#L677)) - opt-in, but it hits EVERY family, this being an ordered COLUMN scan,
   not family = "ordinal". SURFACE draw-law change on that path. COST ~60 engine + tests. TIER Opus.
5. rbart_vi logistic token. Widen the family list, leave "auto" alone. SURFACE one widened formal
   value, additive-later. VALUE grouped WEIGHTED binary is unreachable today. ~30 R, Sonnet.
6. survival `entry=` ([TODO:427](https://github.com/vdorie/dbarts/blob/63df524e8b54948c6bdca1820a949f2c5eb17385/TODO#L427)). Additive-later, face only. ~150 R + Rd + tests, Sonnet.
7. group-by-exposure ([TODO:137](https://github.com/vdorie/dbarts/blob/63df524e8b54948c6bdca1820a949f2c5eb17385/TODO#L137)). Additive-later, largest missing face piece, engine decorator
   already family-agnostic. COST ~500 R + Rd + tests plus a tau `rel.scale` decision. TIER Sonnet +
   Opus design read.
8. sparse-extensions, rbart_vi half ([TODO:419](https://github.com/vdorie/dbarts/blob/63df524e8b54948c6bdca1820a949f2c5eb17385/TODO#L419)). Additive-later widening. TIER Sonnet.
9. summary methods for ordinal / nbinom / hurdle (new). Additive-later, cheap face repair. Sonnet.

STRUCK from the memo's pre-RC list: multinomial twin-create deletion (its rank 5) - draws move only
at `seed = NA`, the win is 3ms, and [docs/design/multinomial-mutation-arc.md:389-398](https://github.com/vdorie/dbarts/blob/63df524e8b54948c6bdca1820a949f2c5eb17385/docs/design/multinomial-mutation-arc.md#L389-L398) already ruled "Defer".
DEFERRED as the memo has them: interaction-constraints heredity and gp-followups sampled
lengthscales, both additive on exported helpers and both large.
