# Lens 2 - backlog and doors value scan (bartcore 7a8c7286, verified; tree clean)

Read in full: TODO; feature-matrix.md's gap/door/refused/open cells + Gaps; every docs/design
doors/deferred/declined/post-1.0/open passage; consolidated-report.md + decision-brief.md; both backlog-value-scan memos; the
newest three landing notes' "Doors left"; NEWS 1.0-0. Every status claim re-verified against live code at this tip.

## What changed since backlog-value-scan-2026-08-24 (ee70e6f7 -> 7a8c7286, 82 commits)

That scan is FULLY ADJUDICATED - nothing in it is open and undecided. Its corrected pre-RC list 1-4 and 9 LANDED (221ec7af
hetero loglik/ppd/summary at s(x); c95a5e83 variance setState column mask; 47cdb96a dbartsData(bases=); b4b9119d+044a9098
ordinal grow-from-root missing rows; 1583140b summary for ordinal/nbinom/hurdle). Its 5-8 were DECIDED by VD 2026-08-24
(release-candidate-review.md:1036-1042): group.by, survival entry=, sparse-extensions' two halves, rbart_vi's logistic token
deferred post-1.0 "window not lock-in"; fit-time test basis deferred; approximate Polya-Gamma DECLINED for 1.0; hurdle
samplerOnly stays refused; twin-create STRUCK. The second review's 8 blockers and 42/44 majors landed in waves 1-3, so all of
section 1 is new.

## 1. PRE-RC by rule (b) - surface-bearing with nameable value

P1 Formula `*` / `:` terms die with R's own "undefined columns selected". R/data.R:1241
   `makeModelMatrix(modelFrame[termLabels])`: termLabels carries "x1:x2", the model frame does not; poly(), ns(), log(), ".",
   offset() all work. ENABLES honouring product terms, or refusing by name - the unnamed internal error is the class NEWS
   1.0-0's BUG FIXES swept everywhere else. SURFACE yes (the formula input contract of dbarts/bart2/dbartsData). ~30 R +
   tests, Sonnet. NOT RECORDED. TAKE.
P2 predict() refuses on a default fit: keepTrees = FALSE (R/bart.R:679, :2685), refusal R/generics.R:271. ENABLES the first
   workflow a new user runs. SURFACE a DEFAULT - what the release locks hardest - with a real cost on the other side (every
   fit carries its trees). FORK for VD: (a) flip it, (b) keep it and make the refusal name the one-argument cure, (c) keep as
   is. ~10 R either way. NOT RECORDED. Put (a)/(b) to VD now; after 1.0 only (b) stays reachable.
P3 Arguments on public generics that silently do nothing: predict(bases=), predict(sample=), fitted(combineChains=) on the
   own-class fits; as_draws_array.bartMultinomial(vars=) ignores a non-meanProb value (release-candidate-review.md:846-850,
   explicitly outside the judgement table). ENABLES a caller who passes them stops getting a silently wrong answer. SURFACE
   yes - refusing later breaks whoever passes them today. ~40 R + tests, Sonnet. TAKE.
P4 Fractional n.threads truncates silently on the six predict generics. R/generics.R:234-247 refuses non-numeric / length!=1 /
   NA / <1 by name, then as.integer() takes 2.7 as 2 while the message says "a single positive integer". The formal is NEW IN
   1.0. SURFACE yes - a validation contract on a brand-new formal. ~5 R + 1 test, Sonnet. Recorded only as "reported, not
   fixed". TAKE.
P5 Two compositions CONSTRUCT unadjudicated: grouped + variance=, and hetero + group.by (feature-matrix.md:727-731; CH:641
   decorates before CH:742 builds the variance forest). ENABLES 1.0 not shipping acceptance of a composition nobody has
   checked. SURFACE yes - acceptance IS the contract. VD's own FX2 posture is the ready-made answer: refuse as a validation
   error, formal stays, no interface friction, door memo. ~20 R + memo, or an Opus adjudication. APPLY FX2's POSTURE.
P6 pdbart/pd2bart reach only bart/bart2 gaussian-binary fits and the sampler (R/partialDependence.R:58, :81-84). rbart is
   class "rbart" (R/rbart.R:1307) and the five own classes do not inherit "bart" (R/bart.R:1729/:1975/:2205/:2330), so all six
   get a generic "must be a matrix, data.frame, formula, fitted bart model, or dbartsSampler". plotTree and
   survivalProbabilities grew by-name refusals in the review; pdbart grew neither. ENABLES BART's primary interpretation
   surface reaching the families 1.0 adds. SURFACE yes. Split: by-name refusal ~30 R pre-RC (Sonnet), real support ~200 R
   post-1.0. NOT RECORDED. TAKE THE REFUSAL.
P7 bcf equivalence baseline predates the statistical mode, and the cross-host channel decision is untaken.
   baselines/bcf-equivalence-6e3b9fb8.rds (Aug 16) carries no summaries; bcf-equivalence.R:500-509 degrades loudly;
   equivalence.yaml:89 pins it and binds to the default branch, so it goes live at the merge. TODO:68-86 adds that a re-record
   alone is insufficient - the snapshot channels (mu, tau, glue, varcount, forestFits, accepted, installed) mismatch
   cross-host by construction. ENABLES the BCF reproducibility gate actually running on x86 CI. SURFACE the baseline FILE
   FORMAT (exempt snapshot channels under a cross-host flag, vs convert to draws-axis recordings). Re-record valid only from
   the recording host. DECIDE FORMAT NOW, re-record at the RC tip.
P8 plot.bart and plot.rbart leak par(mfrow): plotSigmaTrace (R/plot.R:9-12) sets it with no save; callers at :53 and :119 do
   not restore, while plot.bartHurdle (:415-416) and five other sites do. CRAN policy; recorded door at :607. ~8 R + sentinel,
   Sonnet. TAKE - the last residue of a pass that fixed six sites and missed the two most-used ones.

## 2. POST-1.0 additive (nameable value, additive-later by construction)

Settled deferrals (VD 2026-08-24): group.by on bart2/dbarts; survival entry=; sparse-extensions' rbart_vi and linear-leaf
halves; rbart_vi logistic token; fit-time test basis. Families/links: Poisson (negative-binomial.md:635); grouped
ordinal/nbinom/hazard ([f31], ordinal.md:490, negative-binomial.md:641, survival.md:706); NB integer frequency weights (:455);
hurdle.nbinom, gamma part, logistic occupancy, grouped hurdle, Duan smearing (hurdle.md:320-329); competing risks
(survival.md:692, moves with multinomial D2); time-varying covariates / long-format ingestion (:688); probit multinomial path
(multinomial.md:383); xbart ordinal ranked-probability loss (ordinal.md:500). Amplitude coupling
(multiplier-combiner.md:572-589): aft/ordinal/nbinom under the coupling; variance forest under it; the combining-sampler test
surface (bairrtt named); per-forest nameable calibration on combiners. C API additive, all bump MINOR and re-bake by the
header's own rule (dbarts.h:133-137): multinomial creation + K-aware predict (c-api-growth.md:695-722); per-observation
predictor updates, setCutPoints, setData, observer callback (public-surface.md:408-412); forest-indexed predict
(bart-as-a-component.md:167-174). Workflow layer, NOT ON ANY LIST: ICE / centered PDP / >2-variable PDP; update() to extend a
chain; an xbart result class with a best-cell extractor (man/xbart.Rd:134 states the bare array deliberately); interaction
DETECTION over extract(type="trees"); variable-selection inference; as.matrix/as.data.frame on a fit; na.action as a formal
(R/data.R:1126 hard-sets na.pass); an unseen-level escape through the shipped MIA code; logLik/loo wrappers plus the hazard
per-person-period LOO trap (bartcore-review-tour.md:290, "stated nowhere"). Diagnostics/perf: per-sweep vetoed-leaf counter
(none exists in src/); monotone leaf branch-fill; negbin rootogram and burn-in dispersion channel; ordinal log_diff_exp tail
precision; nbinom per-sweep loop collapse (owes bench-sampler); run()'s per-call thread override D1 and a predict-shaped
thread default D2 (threaded-predict.md:245); a per-sweep run callback (correlated-outcomes.md:110). Data handle: serialization
and public exposure of the standalone container stay open as decided (public-surface.md:492-501) - surface on two axes, a new
exported class and a serialized format; and shared MUTABLE codes across attached samplers (data-ownership.md:205-207), which
would collapse bairrtt's two-copy workaround. Research: gp sampled lengthscales; interaction heredity; per-forest row
subsetting - COMPACTION arm only (mask arm permanently refused, zero-weight arm shipped); exact AR-1; tree-mixing perturb /
rotation / heated companion chains; GPU cut-scan; python bindings; informed grow-from-root. Evidence/calibration: hetero SBC
setState lift; the aft SBC arm is BUILT (0045507c) with no run recorded; hurdle's combine/retransform analytic oracle;
runSbcBCF repair; BCF probit/logistic equivalence and SBC arms; exact gates are single-tree in 16 of 16 (both "exceptions"
hold the predictor constant). Undecided singles from the second review: M20, D7, A15, A16, M8-gen (setForestBasis(k, ~var)
evaluates in environment(basis)), C8.

## 3. DECLINE / CLOSED, with the recorded reason

Approximate Polya-Gamma - real dispersion AND real binary weights, one decision: DECLINED for 1.0 (:1040). hurdle samplerOnly:
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
N1 NA in newdata on a column COMPLETE in training routes left, silently. data.hpp:702 codes any NA to naCode at predict;
   tree.hpp:160 sends naCode to missingGoesRight(); tree.hpp:1559 refuses to restore a set missing flag on an NA-free column,
   so the bit is always 0. man/dbarts.Rd:117 claims the opposite. SURFACE a new refusal or new draws; ~20 R + ~10 C++ to
   refuse. Only the mid-run mutation case is recorded (mia-missingness.md:123). Refuse or document, pre-RC.
N2 misc D1: the alpha == 0.0 short-circuit exists only in the scalar fallback (linearAlgebra.c:67, :85), in none of the
   neon/avx/sse2 variants, so 0*x diverges for non-finite x and signed zero - a hole in the stated within-host
   bitwise-across-SIMD invariant. Reachability unassessed; cheap.

Stale records a human reviewer will read as open work.
N3 model-space-survey.md:581-613 still reads as an open RELEASE BLOCKER (variance-forest mutation routing, an out-of-bounds
   write from dbarts(variance=~1) then setData). FIXED AT HEAD: chain.hpp:2540-2551 resizes and re-anchors, :2612 covers
   forceRefreshTrees, :3699 the donor path, resizeVarianceStorage:4402 names the seven allocations. Most alarming stale
   passage in docs/. Fix first.
N4 feature-matrix [f19]:555 and [f28]:694 claim no test pins the student / hetero loglik channel, so those cells stay "?".
   Both are pinned BY VALUE (test-pointwise-loglik.R:386-398; test-heteroscedastic-channels.R:37-53, tol 1e-12). Only 2 of the
   4 "?" cells are real - P5.
N5 Other verified staleness. r-c-division.md carries two: :344-353 marks the getLatents docs slice "STILL OPEN" though
   $getFitsWithoutOffset() ships and its own text names the trap (R/dbarts.R:1696-1705, man/dbartsAugmentation.Rd:31) and
   dbarts.h:770-778 documents the per-family semantics; :292-296's adopted flat rename is DONE (dbarts.h:497-510 carries
   setForestBasis/numForestAmplitudes/forestAmplitudes with setForestWeights kept). Also: feature-matrix.md:503-506 reads as a
   live apiHash defect (impossible - C_interface.cpp:465 static_asserts it); c-api-growth.md:709-712 prices the multinomial C
   door against two superseded ABI literals (live 0x66d33f1613892406, 0xcb83367ee0c4175b); public-surface.md:220-221 still
   asks as "Open:" whether family belongs on xbart, which R/xbart.R:26 has carried since (auto/gaussian/probit/logistic,
   probit the binary default); TODO:233-234 still lists the review-tour refresh as remaining though it landed at 0b89ab8b;
   TODO:153-154 calls the statesAgree above-chain gap "worked around" while test_fuzz.cpp:164-168 documents it as deliberate.
   Verified clean: 70 freshness advisories (tool exits OK), zero TODO/FIXME/XXX in shipped dirs.
N6 nameable-calibration.md:3 says PARTIAL, flat-C half pending; feature-matrix.md:1024 says ARC COMPLETE, four slices landed;
   INDEX.md:73 still calls it and latent-subset-mask "the two designed-but-unbuilt arcs" (both ARC COMPLETE).
N7 docs/plans/INDEX.md: 3 stale rows of 151 - multinomial-level-centering (OPEN vs own LANDED ec2a3d0), grouped-equivalence
   (RESEARCH-OPEN vs own CLOSED), setpredictor-leafof-rebuild (OPEN vs own CLOSED). check-doc-freshness.R checks docs/design
   labels only (37 checked).
N8 Variable-selection inference and random-effects breadth (slopes, crossed, nested) are recorded ONLY in
   docs/plans/roadmap-survey.md, which the TODO does not reference - a TODO-driven census misses them. Both large and
   post-1.0.

## 5. VD-HELD (recommendation attached, not an action)

binary-kforest-k1-reachability - gated on acceptance evidence for a new shipped configuration; the engine takes K=1, the
refusal is spec.R:583-591. The gate as written is the licensed form; leave it, noting a K=1 configuration is surface, so if
ever taken it wants the window. BCF whole-data setData (door 1) - keep shut; the gate is licensed. c-api-growth's multinomial
C door - the gate reads "open it when stan4bart asks", which is consumer-absence and not licensed by rule (a). Outcome
unchanged (additive post-1.0), but restate the gate as cost-and-schedule. Also VD-held: the RC declaration; the coordinated
merge that first registers equivalence/rchk/revdep-smoke/sbc/valgrind (none has ever run anywhere); the DESCRIPTION Date; this
scan itself.
