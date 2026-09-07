# Lens 2 - backlog and doors value scan (bartcore 7a8c7286, verified; tree clean)

Read in full: TODO; feature-matrix.md's gap/door/refused/open cells + Gaps; every docs/design
doors/deferred/declined/post-1.0/open passage; consolidated-report.md + decision-brief.md; both backlog-value-scan memos; the
newest three landing notes' "Doors left"; NEWS 1.0-0. Every status claim re-verified against live code at this tip.

## What changed since backlog-value-scan-2026-08-24 (ee70e6f7 -> 7a8c7286, 82 commits)

That scan is FULLY ADJUDICATED - nothing in it is open and undecided. Its corrected pre-RC list 1-4 and 9 LANDED (221ec7af
hetero loglik/ppd/summary at s(x); c95a5e83 variance setState column mask; 47cdb96a dbartsData(bases=); b4b9119d+044a9098
ordinal grow-from-root missing rows; 1583140b summary for ordinal/nbinom/hurdle). Its 5-8 were DECIDED by VD 2026-08-24
([docs/plans/release-candidate-review.md:1036-1042](https://github.com/vdorie/dbarts/blob/1583140b4d0a48a76cf2eb03e487a295818629cd/docs/plans/release-candidate-review.md#L1036-L1042)): group.by, survival entry=, sparse-extensions' two halves, rbart_vi's logistic token
deferred post-1.0 "window not lock-in"; fit-time test basis deferred; approximate Polya-Gamma DECLINED for 1.0; hurdle
samplerOnly stays refused; twin-create STRUCK. The second review's 8 blockers and 42/44 majors landed in waves 1-3, so all of
section 1 is new.

## 1. PRE-RC by rule (b) - surface-bearing with nameable value

P1 Formula `*` / `:` terms die with R's own "undefined columns selected". [R/data.R:1241](https://github.com/vdorie/dbarts/blob/1583140b4d0a48a76cf2eb03e487a295818629cd/R/data.R#L1241)
   `makeModelMatrix(modelFrame[termLabels])`: termLabels carries "x1:x2", the model frame does not; poly(), ns(), log(), ".",
   offset() all work. ENABLES honouring product terms, or refusing by name - the unnamed internal error is the class NEWS
   1.0-0's BUG FIXES swept everywhere else. SURFACE yes (the formula input contract of dbarts/bart2/dbartsData). ~30 R +
   tests, Sonnet. NOT RECORDED. TAKE.
P2 predict() refuses on a default fit: keepTrees = FALSE ([R/bart.R:679](https://github.com/vdorie/dbarts/blob/1583140b4d0a48a76cf2eb03e487a295818629cd/R/bart.R#L679), [R/bart.R:2685](https://github.com/vdorie/dbarts/blob/1583140b4d0a48a76cf2eb03e487a295818629cd/R/bart.R#L2685)), refusal [R/generics.R:271](https://github.com/vdorie/dbarts/blob/1583140b4d0a48a76cf2eb03e487a295818629cd/R/generics.R#L271). ENABLES the first
   workflow a new user runs. SURFACE a DEFAULT - what the release locks hardest - with a real cost on the other side (every
   fit carries its trees). FORK for VD: (a) flip it, (b) keep it and make the refusal name the one-argument cure, (c) keep as
   is. ~10 R either way. NOT RECORDED. Put (a)/(b) to VD now; after 1.0 only (b) stays reachable.
P3 Arguments on public generics that silently do nothing: predict(bases=), predict(sample=), fitted(combineChains=) on the
   own-class fits; as_draws_array.bartMultinomial(vars=) ignores a non-meanProb value ([docs/plans/release-candidate-review.md:846-850](https://github.com/vdorie/dbarts/blob/1583140b4d0a48a76cf2eb03e487a295818629cd/docs/plans/release-candidate-review.md#L846-L850),
   explicitly outside the judgement table). ENABLES a caller who passes them stops getting a silently wrong answer. SURFACE
   yes - refusing later breaks whoever passes them today. ~40 R + tests, Sonnet. TAKE.
P4 Fractional n.threads truncates silently on the six predict generics. [R/generics.R:234-247](https://github.com/vdorie/dbarts/blob/1583140b4d0a48a76cf2eb03e487a295818629cd/R/generics.R#L234-L247) refuses non-numeric / length!=1 /
   NA / <1 by name, then as.integer() takes 2.7 as 2 while the message says "a single positive integer". The formal is NEW IN
   1.0. SURFACE yes - a validation contract on a brand-new formal. ~5 R + 1 test, Sonnet. Recorded only as "reported, not
   fixed". TAKE.
P5 Two compositions CONSTRUCT unadjudicated: grouped + variance=, and hetero + group.by ([docs/design/feature-matrix.md:727-731](https://github.com/vdorie/dbarts/blob/1583140b4d0a48a76cf2eb03e487a295818629cd/docs/design/feature-matrix.md#L727-L731); [src/bartcore/chain.hpp:641](https://github.com/vdorie/dbarts/blob/1583140b4d0a48a76cf2eb03e487a295818629cd/src/bartcore/chain.hpp#L641)
   decorates before [src/bartcore/chain.hpp:742](https://github.com/vdorie/dbarts/blob/1583140b4d0a48a76cf2eb03e487a295818629cd/src/bartcore/chain.hpp#L742) builds the variance forest). ENABLES 1.0 not shipping acceptance of a composition nobody has
   checked. SURFACE yes - acceptance IS the contract. VD's own FX2 posture is the ready-made answer: refuse as a validation
   error, formal stays, no interface friction, door memo. ~20 R + memo, or an Opus adjudication. APPLY FX2's POSTURE.
P6 pdbart/pd2bart reach only bart/bart2 gaussian-binary fits and the sampler ([R/partialDependence.R:58](https://github.com/vdorie/dbarts/blob/1583140b4d0a48a76cf2eb03e487a295818629cd/R/partialDependence.R#L58), [R/partialDependence.R:81-84](https://github.com/vdorie/dbarts/blob/1583140b4d0a48a76cf2eb03e487a295818629cd/R/partialDependence.R#L81-L84)). rbart is
   class "rbart" ([R/rbart.R:1307](https://github.com/vdorie/dbarts/blob/1583140b4d0a48a76cf2eb03e487a295818629cd/R/rbart.R#L1307)) and the five own classes do not inherit "bart" ([R/bart.R:1729](https://github.com/vdorie/dbarts/blob/1583140b4d0a48a76cf2eb03e487a295818629cd/R/bart.R#L1729)/[R/bart.R:1975](https://github.com/vdorie/dbarts/blob/1583140b4d0a48a76cf2eb03e487a295818629cd/R/bart.R#L1975)/[R/bart.R:2205](https://github.com/vdorie/dbarts/blob/1583140b4d0a48a76cf2eb03e487a295818629cd/R/bart.R#L2205)/[R/bart.R:2330](https://github.com/vdorie/dbarts/blob/1583140b4d0a48a76cf2eb03e487a295818629cd/R/bart.R#L2330)), so all six
   get a generic "must be a matrix, data.frame, formula, fitted bart model, or dbartsSampler". plotTree and
   survivalProbabilities grew by-name refusals in the review; pdbart grew neither. ENABLES BART's primary interpretation
   surface reaching the families 1.0 adds. SURFACE yes. Split: by-name refusal ~30 R pre-RC (Sonnet), real support ~200 R
   post-1.0. NOT RECORDED. TAKE THE REFUSAL.
P7 bcf equivalence baseline predates the statistical mode, and the cross-host channel decision is untaken.
   baselines/bcf-equivalence-6e3b9fb8.rds (Aug 16) carries no summaries; [benchmarks/R/bcf-equivalence.R:500-509](https://github.com/vdorie/dbarts/blob/ecdfb9454385326580678b2d72c627bfb108990b/benchmarks/R/bcf-equivalence.R#L500-L509) degrades loudly;
   [.github/workflows/equivalence.yaml:89](https://github.com/vdorie/dbarts/blob/ecdfb9454385326580678b2d72c627bfb108990b/.github/workflows/equivalence.yaml#L89) pins it and binds to the default branch, so it goes live at the merge. [TODO:68-86](https://github.com/vdorie/dbarts/blob/ecdfb9454385326580678b2d72c627bfb108990b/TODO#L68-L86) adds that a re-record
   alone is insufficient - the snapshot channels (mu, tau, glue, varcount, forestFits, accepted, installed) mismatch
   cross-host by construction. ENABLES the BCF reproducibility gate actually running on x86 CI. SURFACE the baseline FILE
   FORMAT (exempt snapshot channels under a cross-host flag, vs convert to draws-axis recordings). Re-record valid only from
   the recording host. DECIDE FORMAT NOW, re-record at the RC tip.
P8 plot.bart and plot.rbart leak par(mfrow): plotSigmaTrace ([R/plot.R:9-12](https://github.com/vdorie/dbarts/blob/ecdfb9454385326580678b2d72c627bfb108990b/R/plot.R#L9-L12)) sets it with no save; callers at [R/plot.R:53](https://github.com/vdorie/dbarts/blob/ecdfb9454385326580678b2d72c627bfb108990b/R/plot.R#L53) and [R/plot.R:119](https://github.com/vdorie/dbarts/blob/ecdfb9454385326580678b2d72c627bfb108990b/R/plot.R#L119) do
   not restore, while plot.bartHurdle ([R/plot.R:415-416](https://github.com/vdorie/dbarts/blob/ecdfb9454385326580678b2d72c627bfb108990b/R/plot.R#L415-L416)) and five other sites do. CRAN policy; recorded door, location not placeable at this sha - R/plot.R is only 549 lines there (unresolved: [R/plot.R:607](https://github.com/vdorie/dbarts/blob/ecdfb9454385326580678b2d72c627bfb108990b/R/plot.R#L607)). ~8 R + sentinel,
   Sonnet. TAKE - the last residue of a pass that fixed six sites and missed the two most-used ones.

## 2. POST-1.0 additive (nameable value, additive-later by construction)

Settled deferrals (VD 2026-08-24): group.by on bart2/dbarts; survival entry=; sparse-extensions' rbart_vi and linear-leaf
halves; rbart_vi logistic token; fit-time test basis. Families/links: Poisson ([docs/design/negative-binomial.md:635](https://github.com/vdorie/dbarts/blob/ecdfb9454385326580678b2d72c627bfb108990b/docs/design/negative-binomial.md#L635)); grouped
ordinal/nbinom/hazard ([f31], [docs/design/ordinal.md:490](https://github.com/vdorie/dbarts/blob/ecdfb9454385326580678b2d72c627bfb108990b/docs/design/ordinal.md#L490), [docs/design/negative-binomial.md:641](https://github.com/vdorie/dbarts/blob/ecdfb9454385326580678b2d72c627bfb108990b/docs/design/negative-binomial.md#L641), [docs/design/survival.md:706](https://github.com/vdorie/dbarts/blob/ecdfb9454385326580678b2d72c627bfb108990b/docs/design/survival.md#L706)); NB integer frequency weights ([docs/design/survival.md:455](https://github.com/vdorie/dbarts/blob/ecdfb9454385326580678b2d72c627bfb108990b/docs/design/survival.md#L455));
hurdle.nbinom, gamma part, logistic occupancy, grouped hurdle, Duan smearing ([docs/design/hurdle.md:320-329](https://github.com/vdorie/dbarts/blob/ecdfb9454385326580678b2d72c627bfb108990b/docs/design/hurdle.md#L320-L329)); competing risks
([docs/design/survival.md:692](https://github.com/vdorie/dbarts/blob/ecdfb9454385326580678b2d72c627bfb108990b/docs/design/survival.md#L692), moves with multinomial D2); time-varying covariates / long-format ingestion ([docs/design/survival.md:688](https://github.com/vdorie/dbarts/blob/ecdfb9454385326580678b2d72c627bfb108990b/docs/design/survival.md#L688)); probit multinomial path
([docs/design/multinomial.md:383](https://github.com/vdorie/dbarts/blob/ecdfb9454385326580678b2d72c627bfb108990b/docs/design/multinomial.md#L383)); xbart ordinal ranked-probability loss ([docs/design/ordinal.md:500](https://github.com/vdorie/dbarts/blob/ecdfb9454385326580678b2d72c627bfb108990b/docs/design/ordinal.md#L500)). Amplitude coupling
([docs/design/multiplier-combiner.md:572-589](https://github.com/vdorie/dbarts/blob/ecdfb9454385326580678b2d72c627bfb108990b/docs/design/multiplier-combiner.md#L572-L589)): aft/ordinal/nbinom under the coupling; variance forest under it; the combining-sampler test
surface (bairrtt named); per-forest nameable calibration on combiners. C API additive, all bump MINOR and re-bake by the
header's own rule ([inst/include/dbarts/dbarts.h:133-137](https://github.com/vdorie/dbarts/blob/ecdfb9454385326580678b2d72c627bfb108990b/inst/include/dbarts/dbarts.h#L133-L137)): multinomial creation + K-aware predict ([docs/plans/archive/c-api-growth.md:695-722](https://github.com/vdorie/dbarts/blob/758bccdda099af6f8184eca3d9fcac4aecb65094/docs/plans/archive/c-api-growth.md#L695-L722)); per-observation
predictor updates, setCutPoints, setData, observer callback ([docs/design/public-surface.md:408-412](https://github.com/vdorie/dbarts/blob/ecdfb9454385326580678b2d72c627bfb108990b/docs/design/public-surface.md#L408-L412)); forest-indexed predict
([docs/design/bart-as-a-component.md:167-174](https://github.com/vdorie/dbarts/blob/ecdfb9454385326580678b2d72c627bfb108990b/docs/design/bart-as-a-component.md#L167-L174)). Workflow layer, NOT ON ANY LIST: ICE / centered PDP / >2-variable PDP; update() to extend a
chain; an xbart result class with a best-cell extractor ([man/xbart.Rd:134](https://github.com/vdorie/dbarts/blob/ecdfb9454385326580678b2d72c627bfb108990b/man/xbart.Rd#L134) states the bare array deliberately); interaction
DETECTION over extract(type="trees"); variable-selection inference; as.matrix/as.data.frame on a fit; na.action as a formal
([R/data.R:1126](https://github.com/vdorie/dbarts/blob/ecdfb9454385326580678b2d72c627bfb108990b/R/data.R#L1126) hard-sets na.pass); an unseen-level escape through the shipped MIA code; logLik/loo wrappers plus the hazard
per-person-period LOO trap ([docs/plans/bartcore-review-tour.md:290](https://github.com/vdorie/dbarts/blob/ecdfb9454385326580678b2d72c627bfb108990b/docs/plans/bartcore-review-tour.md#L290), "stated nowhere"). Diagnostics/perf: per-sweep vetoed-leaf counter
(none exists in src/); monotone leaf branch-fill; negbin rootogram and burn-in dispersion channel; ordinal log_diff_exp tail
precision; nbinom per-sweep loop collapse (owes bench-sampler); run()'s per-call thread override D1 and a predict-shaped
thread default D2 ([docs/design/threaded-predict.md:245](https://github.com/vdorie/dbarts/blob/ecdfb9454385326580678b2d72c627bfb108990b/docs/design/threaded-predict.md#L245)); a per-sweep run callback ([docs/design/correlated-outcomes.md:110](https://github.com/vdorie/dbarts/blob/ecdfb9454385326580678b2d72c627bfb108990b/docs/design/correlated-outcomes.md#L110)). Data handle: serialization
and public exposure of the standalone container stay open as decided ([docs/design/public-surface.md:492-501](https://github.com/vdorie/dbarts/blob/ecdfb9454385326580678b2d72c627bfb108990b/docs/design/public-surface.md#L492-L501)) - surface on two axes, a new
exported class and a serialized format; and shared MUTABLE codes across attached samplers ([docs/design/data-ownership.md:205-207](https://github.com/vdorie/dbarts/blob/ecdfb9454385326580678b2d72c627bfb108990b/docs/design/data-ownership.md#L205-L207)), which
would collapse bairrtt's two-copy workaround. Research: gp sampled lengthscales; interaction heredity; per-forest row
subsetting - COMPACTION arm only (mask arm permanently refused, zero-weight arm shipped); exact AR-1; tree-mixing perturb /
rotation / heated companion chains; GPU cut-scan; python bindings; informed grow-from-root. Evidence/calibration: hetero SBC
setState lift; the aft SBC arm is BUILT (0045507c) with no run recorded; hurdle's combine/retransform analytic oracle;
runSbcBCF repair; BCF probit/logistic equivalence and SBC arms; exact gates are single-tree in 16 of 16 (both "exceptions"
hold the predictor constant). Undecided singles from the second review: M20, D7, A15, A16, M8-gen (setForestBasis(k, ~var)
evaluates in environment(basis)), C8.

## 3. DECLINE / CLOSED, with the recorded reason

Approximate Polya-Gamma - real dispersion AND real binary weights, one decision: DECLINED for 1.0, location not placeable at this sha - data-ownership.md is only 340 lines there (unresolved: [docs/design/data-ownership.md:1040](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/docs/design/data-ownership.md#L1040)). hurdle samplerOnly:
stays refused. Multinomial twin-create deletion: STRUCK as a relitigation (draws move only at seed = NA; 3ms). $getLatents
build: written decline, VD 2026-08-20. forest-ranef-interweaving: NO-GO with a benchmark gate and a VD-sign-off reopen clause.
BCF whole-data setData (door 1) and setData on CSC/mixed (door 3): KEEP UNDESIGNED, gated on a considered failure to find an
enabled model - the licensed form. Per-observation CSC mutation: refused with a model-space reason. bcf() verb: relocated to
bartCause. Flat forests= on bart2(): declined-but-addable (XOR against the formula route). Per-draw amplitude channel in flat
C: DECLINED. Multinomial per-forest off-sample replay: standing refusal with a named re-open trigger. Also closed: annealed
burn-in; continuous-time birth-death; fp32 scratch bundle; leafOf uint16; blocked-jacobi; grow-from-root-default; x86-simd as
a perf lever; within-chain threading (0.91x x86 / 1.10x arm); soft trees; zero-inflation; cross-version state migration shims;
cloglog and continuous-time survival.

## 4. NEW - not on any list before this scan

Live code: P1, P2, P4, P6 above.
N1 NA in newdata on a column COMPLETE in training routes left, silently. [src/bartcore/data.hpp:702](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/src/bartcore/data.hpp#L702) codes any NA to naCode at predict;
   [src/bartcore/tree.hpp:160](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/src/bartcore/tree.hpp#L160) sends naCode to missingGoesRight(); [src/bartcore/tree.hpp:1559](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/src/bartcore/tree.hpp#L1559) refuses to restore a set missing flag on an NA-free column,
   so the bit is always 0. [man/dbarts.Rd:117](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/man/dbarts.Rd#L117) claims the opposite. SURFACE a new refusal or new draws; ~20 R + ~10 C++ to
   refuse. Only the mid-run mutation case is recorded ([docs/design/mia-missingness.md:123](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/docs/design/mia-missingness.md#L123)). Refuse or document, pre-RC.
N2 misc D1: the alpha == 0.0 short-circuit exists only in the scalar fallback ([src/misc/linearAlgebra.c:67](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/src/misc/linearAlgebra.c#L67), [src/misc/linearAlgebra.c:85](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/src/misc/linearAlgebra.c#L85)), in none of the
   neon/avx/sse2 variants, so 0*x diverges for non-finite x and signed zero - a hole in the stated within-host
   bitwise-across-SIMD invariant. Reachability unassessed; cheap.

Stale records a human reviewer will read as open work.
N3 [docs/design/model-space-survey.md:581-613](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/docs/design/model-space-survey.md#L581-L613) still reads as an open RELEASE BLOCKER (variance-forest mutation routing, an out-of-bounds
   write from dbarts(variance=~1) then setData). FIXED AT HEAD: [src/bartcore/chain.hpp:2540-2551](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/src/bartcore/chain.hpp#L2540-L2551) resizes and re-anchors, [src/bartcore/chain.hpp:2612](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/src/bartcore/chain.hpp#L2612) covers
   forceRefreshTrees, [src/bartcore/chain.hpp:3699](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/src/bartcore/chain.hpp#L3699) the donor path, resizeVarianceStorage:4402 names the seven allocations. Most alarming stale
   passage in docs/. Fix first.
N4 feature-matrix [f19]:555 and [f28]:694 claim no test pins the student / hetero loglik channel, so those cells stay "?".
   Both are pinned BY VALUE ([inst/tinytest/test-pointwise-loglik.R:386-398](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/inst/tinytest/test-pointwise-loglik.R#L386-L398); [inst/tinytest/test-heteroscedastic-channels.R:37-53](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/inst/tinytest/test-heteroscedastic-channels.R#L37-L53), tol 1e-12). Only 2 of the
   4 "?" cells are real - P5.
N5 Other verified staleness. r-c-division.md carries two: [docs/design/r-c-division.md:329-336](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/docs/design/r-c-division.md#L329-L336) marks the getLatents docs slice "STILL OPEN" though
   $getFitsWithoutOffset() ships and its own text names the trap ([R/dbarts.R:1696-1705](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/dbarts.R#L1696-L1705), [man/dbartsAugmentation.Rd:31](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/man/dbartsAugmentation.Rd#L31)) and
   [inst/include/dbarts/dbarts.h:770-778](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/inst/include/dbarts/dbarts.h#L770-L778) documents the per-family semantics; [inst/include/dbarts/dbarts.h:292-296](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/inst/include/dbarts/dbarts.h#L292-L296)'s adopted flat rename is DONE ([inst/include/dbarts/dbarts.h:497-510](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/inst/include/dbarts/dbarts.h#L497-L510) carries
   setForestBasis/numForestAmplitudes/forestAmplitudes with setForestWeights kept). Also: [docs/design/feature-matrix.md:503-506](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/docs/design/feature-matrix.md#L503-L506) reads as a
   live apiHash defect (impossible - [src/C_interface.cpp:465](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/src/C_interface.cpp#L465) static_asserts it); [docs/plans/archive/c-api-growth.md:709-712](https://github.com/vdorie/dbarts/blob/758bccdda099af6f8184eca3d9fcac4aecb65094/docs/plans/archive/c-api-growth.md#L709-L712) prices the multinomial C
   door against two superseded ABI literals (live 0x66d33f1613892406, 0xcb83367ee0c4175b); [docs/design/public-surface.md:220-221](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/docs/design/public-surface.md#L220-L221) still
   asks as "Open:" whether family belongs on xbart, which [R/xbart.R:26](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/xbart.R#L26) has carried since (auto/gaussian/probit/logistic,
   probit the binary default); [TODO:233-234](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/TODO#L233-L234) still lists the review-tour refresh as remaining though it landed at 0b89ab8b;
   [TODO:153-154](https://github.com/vdorie/dbarts/blob/0b89ab8be60ec469ae95971b2919c8242c543759/TODO#L153-L154) calls the statesAgree above-chain gap "worked around" while [tests/cpp/test_fuzz.cpp:164-168](https://github.com/vdorie/dbarts/blob/0b89ab8be60ec469ae95971b2919c8242c543759/tests/cpp/test_fuzz.cpp#L164-L168) documents it as deliberate.
   Verified clean: 70 freshness advisories (tool exits OK), zero TODO/FIXME/XXX in shipped dirs.
N6 nameable-calibration.md's Status line says PARTIAL, flat-C half pending (its location
   could not be placed at 0b89ab8b: unresolved: [docs/design/nameable-calibration.md:3](https://github.com/vdorie/dbarts/blob/0b89ab8be60ec469ae95971b2919c8242c543759/docs/design/nameable-calibration.md#L3)); [docs/design/feature-matrix.md:1024](https://github.com/vdorie/dbarts/blob/0b89ab8be60ec469ae95971b2919c8242c543759/docs/design/feature-matrix.md#L1024) says ARC COMPLETE, four slices landed;
   docs/plans/INDEX.md (its location could not be placed at 0b89ab8b: unresolved: [docs/plans/INDEX.md:73](https://github.com/vdorie/dbarts/blob/0b89ab8be60ec469ae95971b2919c8242c543759/docs/plans/INDEX.md#L73))
   still calls it and latent-subset-mask "the two designed-but-unbuilt arcs" (both ARC COMPLETE).
N7 docs/plans/INDEX.md: 3 stale rows of 151 - multinomial-level-centering (OPEN vs own LANDED ec2a3d0), grouped-equivalence
   (RESEARCH-OPEN vs own CLOSED), setpredictor-leafof-rebuild (OPEN vs own CLOSED). check-doc-freshness.R checks docs/design
   labels only (37 checked).
N8 Variable-selection inference and random-effects breadth (slopes, crossed, nested) are recorded ONLY in
   docs/plans/archive/roadmap-survey.md, which the TODO does not reference - a TODO-driven census misses them. Both large and
   post-1.0.

## 5. VD-HELD (recommendation attached, not an action)

binary-kforest-k1-reachability - gated on acceptance evidence for a new shipped configuration; the engine takes K=1, the
refusal is [R/spec.R:583-591](https://github.com/vdorie/dbarts/blob/ecdfb9454385326580678b2d72c627bfb108990b/R/spec.R#L583-L591). The gate as written is the licensed form; leave it, noting a K=1 configuration is surface, so if
ever taken it wants the window. BCF whole-data setData (door 1) - keep shut; the gate is licensed. c-api-growth's multinomial
C door - the gate reads "open it when stan4bart asks", which is consumer-absence and not licensed by rule (a). Outcome
unchanged (additive post-1.0), but restate the gate as cost-and-schedule. Also VD-held: the RC declaration; the coordinated
merge that first registers equivalence/rchk/revdep-smoke/sbc/valgrind (none has ever run anywhere); the DESCRIPTION Date; this
scan itself.
