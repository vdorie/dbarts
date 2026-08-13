# multiforest-extension-surface

agent: M0 sonnet (docs + vignette + pins, no src). M1 sonnet (one R5 method and
  its mirror). M2 opus (creation-surface replacement, one resolver, a wide
  refusal list). M3 opus, folded INTO dbarts-h-reshape S1 (header re-sign, one
  re-bake). M4 opus throughout (engine), RESOLVED pre-release arc (VD
  2026-08-11, fork 1), own slice list, scheduled after the reshape.
  Serialized: one implementer, each slice lands before the next starts.
rng: M0, M1 NEUTRAL (no src). M2 NEUTRAL and BITWISE-GATED: the new creation
  route resolves to the SAME spec the treatment= route resolves to, so BCF must
  come out bitwise identical on all six bcf-equivalence channels; a divergence
  is a LEAK, never a re-record. M3 NEUTRAL (header and C_interface only; the
  trio is a formality, labelled one, per that arc's rng note). M4 owns the only
  re-record in this arc and only if the amplitude conditional cannot be made
  bitwise at K = 2; every other M4 slice gates the trio bitwise.
window: M0, M1, M2 land INSIDE the pre-release breaking window. M3 rides
  dbarts-h-reshape S1, which is the SECOND AND LAST re-bake of that window, so
  the mean channel's flat spelling is decided there or never (pre-release).
  M4 is a PRE-RELEASE arc, scheduled after M3/the reshape (RESOLVED
  2026-08-11), and needs no header movement if M3 lands. Everything is
  breakable and the four sister packages migrate in lockstep at the freeze,
  once, against the final header (VD 2026-08-10).
budget: M0 ~220 design doc + ~180 vignette + ~200 test. M1 ~75-85 R + ~50-60
  man + ~140 test, ~270-290 total (RESTATED 2026-08-13 pre-M1, superseding the
  ~60 + ~40 + ~140 committed 2026-08-11; the carried items are in M1 item 5).
  M2 ~180 R + ~60 man + ~260 test (test-bcf-creation.R rewritten).
  M3 ~40 header + ~80 C_interface + ~70 consumer.c + ~60 test-capi.R, inside
  dbarts-h-reshape S1's own budget. Total in-repo before the freeze ~1630
  (~1590 before M1's 2026-08-13 restatement).
  M4 priced in its own section (~600-900), not committed here.
tip: every anchor below re-read at db81bfe in this worktree, with src/ read
  through `git show HEAD:` because a code implementer is editing
  src/{R_interface_bartcore.cpp,bartcore/{chain,combiner,facade,sampler}.hpp}
  concurrently (bcf-public-surface S4, in flight). Design artifacts (memo,
  blind refuting critique) sit at `.claude/multiforest-extension-surface-design/`,
  which is GIT-IGNORED; this file is the only tracked record and carries the
  load-bearing facts rather than pointing at them.
landed: 2026-08-11 at e2cc1de by the orchestrator, verbatim from the synthesis
  except this note. bcf-public-surface S4, in flight when the tip note was
  written, landed 1df9c0c with CI all six green. The TODO door sentence cited
  below at :625 sat at :629 at the landing tip, and is now the consolidated
  multiforest-extension-surface TODO entry at :154-178, rewritten again by
  the records commit below; at 934a02d5 that entry sits at TODO:190-214. All
  eight cross-plan amendments applied at this commit.
records: anchor-refresh pass applied 2026-08-11, in the records commit that
  also lands the four VD forks (see "Open decisions" below) and the R/C-
  division principle. 31 line-number corrections applied against the live
  tree at 33f6fdc (table authored in
  scratchpad/extension-surface-anchor-refresh.md): every chain.hpp/
  combiner.hpp anchor drifted by the in-flight bcf-public-surface S4 landing
  is corrected below to its 33f6fdc position; dbarts.h,
  R_interface_bartcore.cpp, R_interface_bartcore_common.hpp and
  C_interface.cpp anchors were RE-VERIFIED against 33f6fdc (the stop-loss
  slice that landed there touched all four files) and none had drifted.
amended: 2026-08-13, PRE-M1, by the orchestrator. Every anchor in this file
  re-verified BY SYMBOL against 934a02d5 and corrected where it had drifted:
  the latent-subset-mask arc (S0-S4, through 414c50a6) and reshape S0
  (a262cd26) moved chain.hpp, sampler.hpp, facade.hpp and
  R_interface_bartcore.cpp under the 33f6fdc numbers. `chain.hpp:870` was a
  FALSE-POSITIVE anchor: at HEAD it resolves to `supportsResponseMutation`,
  and `Chain::setForestWeights` is `:966`. The M1 section additionally gains
  the SETTLED forest-index base, the mirroring mechanism, the three mask-arc
  interactions, the multinomial refusal siting and an explicit restated
  budget; nothing else about the arc's content or ordering moved.

## The question, answered

VD asked: *"Is there a way to ship a multiforest ABI or API that model writers
would want to use?"*

**Yes, and it is already almost shipped.** The surface model writers use is the
one where the CALLER owns the outer loop and dbarts owns one whole coupled
sweep. Four of four sister packages drive dbarts that way; the reference
ecosystem's only two extension surfaces have that shape (stochtree's
`ForestModel$sample_one_iteration`, SoftBart's `forest$do_gibbs`); zero BART
packages ship extension headers; and bcf-public-surface S1 to S3 already made a
multi-forest sampler creatable and drivable from both R and C. What is missing
is not mechanism but a CONTRACT (what may be mutated between sweeps, what a
rebuild drops, what the driver loop costs) and two channels whose reach is
inconsistent between R and C.

The second surface, the one that carries the published models, is a
DECLARATIVE per-forest basis: each forest's contribution to the mean is scaled
by caller-supplied data columns with their own coefficients. Its Gaussian half
is now measured to be ergonomics over a pure-R composition that works today
(1.39-1.43x at K = 2, and a slightly different prior). Its irreducible engine
content is the NON-GAUSSIAN case, where the latents are drawn once against the
COMBINED fit, and that is the content worth building: with an engine-provided
family selected by `family =`, an author declares a non-Gaussian multiforest
model in R by supplying basis COLUMNS as data instead of writing C++.

Both surfaces sit on the additive side of VD's boundary (binding decision 1):
R-first covers additive mean structure, while response shape and priors are
engine axes. Shipped C++ extension headers are refuted on the build and are not
the graduation step for either: `dbarts.h` is the graduation step for
performance, and contributing a family or a leaf model to dbarts is the path for
authoring one.

## Binding decisions inherited (do not reopen)

1. **The prototyping principle (VD 2026-08-10), near verbatim:** *"to the degree
   possible I'd like R to be a way that authors might rapidly prototype and
   develop, with C/C++ as a future step for performance gains. That might shape
   what and how things are supported."* Applications, binding on every slice
   below: (a) the R surface is the on-ramp and must never be LESS capable than
   the flat C surface; (b) a documented R composition cost is a prototyping
   price, not a defect; (c) the graduation target for performance is
   `dbarts.h`, not C++ headers; (d) anything reachable only from C is a gap to
   record.
1b. **Its SCOPE LIMIT (VD 2026-08-10), near verbatim:** *"I can see where R-first
   doesn't work well if there are model details that are not just additive, like
   the shape of the response or the priors."* This aligns the principle with
   binding decision 2 and draws the design boundary this plan is built to, in
   three categories:
   - **Additive composition is R's side.** Offsets, observation and per-forest
     weights, basis or multiplier columns, residual and latent-covariate
     handoffs, and the outer loop that sequences them. This is exactly where the
     measured pure-R path works, and where the on-ramp is documented (M0).
   - **CONFIGURING an engine-provided family or prior is R's side too**, because
     it is argument passing: `family =`, `node.prior = normal(k)`,
     `tree.prior = cgm(base, power)`, `n.trees`, and the per-forest versions of
     those inside `forest()`.
   - **AUTHORING a response shape or a prior is the engine's side**, with C++ as
     its path: a new family, link or latent structure; a new leaf model; a new
     hyperprior law; a new coupling. The path is to add it to dbarts (the
     internal seams are the `model.hpp` concepts and the `facade.hpp` factories),
     not to reach it from a consumer package.
   Consequences applied below, so the design is not contorted against this line:
   the per-sweep R hook is NOT priced as a route to prototyping families or
   priors (its honest scope is additive and glue-adjacent per-sweep
   intervention); the declarative K-forest surface is NOT specified to express
   caller-authored priors; and supplying basis columns against an
   engine-provided family stays on R's side because it is DATA, while authoring
   the family is not.
2. **dbarts stays an engine.** Response models switched by `family`, leaf and
   parameter models switched by prior arguments, mean structure declared as
   data. No estimands, no causal vocabulary, no domain fit functions.
3. **`treatment =` / `moderators =` / `treatmentForest =` on `dbarts()` and
   `dbartsSpec()` are PROVISIONAL and scheduled for removal** once BCF finds its
   home, sequenced AFTER a replacement creation route exists
   (bcf-public-surface.md, Open decision, AMENDED block, VD 2026-08-10). The
   S1/S2 mechanism survives re-skinning; only public names move.
4. **No version-constant increments pre-release** (VD 2026-08-10):
   `DBARTS_C_API_MAJOR` stays 1, `DBARTS_C_API_MINOR` stays 0, and
   `dbarts_apiHash()` is the only runtime lockstep signal.
5. **Everything pre-release is breakable; the sister packages migrate in
   lockstep; consumer compatibility is a cost to ENUMERATE, never a design
   constraint** (VD 2026-08-10). Consequence used below: `treatment =` needs no
   deprecation cycle, because nothing is released and bartCause has not adopted
   it (its `bcf` slot is still one commented line, `R/bartc.R`).
6. **Ship what makes the package useful** (VD), read together with the
   enabling-value gate: gate an item only if nothing valuable can be named that
   it enables. Absence of a consumer today is never the gating fact.
7. **bcf-public-surface S3 LANDED (1622eb9).** `inst/include/dbarts/dbarts.h`
   carries `dbarts_sampler_numForests` (:264), `setTreatment` (:266),
   `forestFits` (:268), `bcfGlue` (:271), the widened
   `setResponse(sampler, y, updateScale)` and `DBARTS_C_API_HASH
   0x1a911c00bb26dcd7ULL` (:83). S4 (per-draw reporting) is in flight; S5
   (`bcf()` and the fit class) is PAUSED on this arc.
8. **The committed queue is not displaced.** bcf-public-surface S4, then
   multiforest-predictor-mutation S0-S4, then dbarts-h-reshape S0-S2, then the
   1.0-0 freeze. This arc's pre-freeze slices slot into that order (see
   "Sequencing"), and M3 is an item INSIDE the reshape's S1.

## Adjudication: the critique's six blocking findings

Every anchor re-verified here at db81bfe, and again at 934a02d5 (2026-08-13);
where a document's own line numbers had drifted the current ones are given.

**B1 (route 3's held door for a state-dependent multiplier callback): ADOPTED
in full, and the door is DELETED.** Verified independently: `Chain::run`
evaluates `onSweep` at the TOP of each sweep, before the forest loop
(`chain.hpp:1110-1111`, `run` at `:1095`, forest loop at `:1145`), and the
previous sweep ended with `combinedFits` (`:1261`), `drawGlue` and
`afterCombine` (`:1275-1276`, line numbers re-verified at 934a02d5). So the
memo's stated hook point ("after `combinedFits` and before the next
`formForestResponse`") IS the sweep boundary, which the shipped
`dbarts_sampler_callback` already occupies. The ABI does not refuse it; it
charges a documented price (`sampler.hpp:260-263` inline-only;
`dbarts.h:372-378` the same refusal predicate). Nothing needs building inside
the sweep, so the memo's per-observation-dispatch objection describes a hook
nobody needs. The one class the door was held for (LongBet-style panel with a
sampled time factor) needs a combiner-internal draw, not a host hook. RE-SCOPED
as an M4 spec question: "a per-forest amplitude with a caller-supplied
covariance kernel, drawn in the combiner", structurally what BCF's `a` and its
half-Cauchy auxiliary already are (`combiner.hpp:589-607`).

**B2 (open decision 5's vector-leaf pricing): ADOPTED, and the door is
re-priced and narrowed.** Verified: `VectorLeafModel` requires
`{ leaf.fitForObservation(v, i) } -> std::same_as<double>` (`model.hpp:60-67`),
so dbarts' "vector leaf" is a vector of PARAMETERS producing a SCALAR fit;
`ResponseFamily` has no multivariate member (`model.hpp:2523`). McJames et al.'s
vector axis is the OUTCOME dimension, so lifting a leaf assert cannot buy
multivariate-outcome BCF: a d-dimensional outcome needs a d-dimensional
response, d-dimensional fits, a d x d residual covariance and a Results layout.
**That half of the door is struck as a category error.** There are TWO asserts
(`combiner.hpp:502` BCF, `:825` multinomial), and the naive lift is
undefined behaviour, not a configuration change: `muByTree` is documented
"Empty for vector and function leaves, which keep the dense `treeFits` slab"
(`combiner.hpp:174-181`) while `BCFForestCombiner::afterCombine` indexes
`forest.muByTree[t]` unconditionally (`:651`) and rescales it (`:682-683`), as
does the multinomial combiner (`:1166`, `:1189`), and the saved-slot rescale
walks `FlatNode`s carrying one scalar `value`. ADDED here, and it matters for
VD's fork 1: multi-arm treatment on a vector leaf and the per-forest basis
channel are DIFFERENT arcs serving different classes (one partition shared
across arms versus one partition per ensemble), so neither subsumes the other
and the vector-leaf arc is not an alternative to the family.

**B3 (the 1.04x attribution): ADOPTED in full.** The receipt
`bcf-public-surface.md:85-93` measures `run(0,1)` in a loop against `run(0,n)`
on ONE BCF sampler, not K single-forest samplers. The plan's numbers are now
the critique's measured ones: engine per-sweep drive **0.89-0.95x** of a
batched run (free, and inside noise), pure-R K-sampler composition
**1.39-1.43x** at K = 2, growing with K (K stores, K `.Call` round trips, O(nK)
R-level residual arithmetic per sweep). Every occurrence of the memo's 1.04x
attributed to an R composition is void.

**B4 (the amendments to the committed plans): ADOPTED, and the plan conflict is
ADJUDICATED.** (a) The memo's "S3 has not started" amendment and its Q7 are
MOOT: S3 landed. (b) The memo cited `dbarts-h-reshape.md` binding decision 1
(reshape-by-replacement) as authority for re-signing entries that decision 3
carves out by name; verified verbatim, decision 3: *"bcf-public-surface S3's
entries are adopted VERBATIM ... This arc re-signs none of them and re-relaxes
none of their guards."* (c) The conflict is real:
`bcf-public-surface.md`'s AMENDED block says *"The flat C names S3 shipped
(`setTreatment`, `bcfGlue`) are likewise renameable at the queued dbarts.h
reshape re-bake if the settled surface uses engine vocabulary."*

**ADJUDICATION: `dbarts-h-reshape.md` decision 3 YIELDS, narrowly and
conditionally.** Three grounds, in order of weight: (i) the bcf-public-surface
amendment records a VD decision dated 2026-08-10, while reshape decision 3 is
an orchestrator sequencing convenience written before that amendment existed, so
the later VD record governs; (ii) the reshape's own decision 1 makes re-signing
in place the sanctioned pre-release mode, and its decision 2 plus its window
make consumer compatibility a cost rather than a constraint; (iii) the reshape's
S1 is the LAST re-bake, so "never re-signed" means "post-release MINOR bump plus
a lockstep bump of stan4bart's floor" - the exact cost that plan's own resolved
question 2 refused to pay in order to build `setForestWeights` in-window.
Symmetry forces the same answer for the mean channel. The yield is NARROW: only
the two BCF-SPECIFIC names and the creation contract may be re-signed;
`numForests` and `forestFits` are already engine vocabulary, and the widened
`setResponse` and the relaxed guards are adopted verbatim exactly as decision 3
says. It is CONDITIONAL: if VD's answer to fork 3 has not landed when the
reshape's S1 starts, the shipped names stay and the rename is dropped, with the
post-release append recorded as the accepted cost. Amendment text for both
documents is in "Cross-plan amendments".

**The rename surface is THREE things, not two.** Beyond `setTreatment` and
`bcfGlue`, `dbarts.h:348-357` documents BCF creation as *"the R specification
dbartsSpec(data, control, treatment = z)"*, so the causal vocabulary is in the
shipped header's CREATION contract; `dbarts.h:43` also documents "setTreatment
copies the 0/1 values it is handed". bcf-public-surface's amendment misses this
and is amended below.

**B5 (the naming collision is a schedule blocker): the FINDING is adopted, its
CONSEQUENCE is narrowed.** Verified: the precision channel is entirely
`dbarts:::`-only today - `Chain::setForestWeights` (`chain.hpp:966`; the
`:870` this plan cited before the 2026-08-13 refresh is a FALSE-POSITIVE
anchor at HEAD, where it resolves to `supportsResponseMutation`), bridge
`bartcore_setForestWeights` (`R_interface_bartcore.cpp:3688-3728`), R
`bartcoreSetForestWeights` (`R/bartcore.R:999`,
absent from `NAMESPACE`, with no R5 method) - and the flat entry does not exist:
`zero-weight-exactness.md` only reserves it and `dbarts-h-reshape.md` S1 item 5
BUILDS it, in the window's last re-bake. So the precision channel's public
spelling does freeze there. NARROWED: the collision the memo feared
(`setForestWeights` versus `setForestMultipliers`, one word apart) is dissolved
by naming the mean channel **basis** rather than multiplier or weights.
"Weights" in dbarts already means precisions - `setWeights` installs
observation precisions, and `setForestWeights` is its per-forest sibling,
composing as `(w_i a_i) m_f^2 s_i` (`chain.hpp:948-950`,
`composeForestWeights` at `:3431-3447`, fn at `:3439`; the `a_i` term is the
active-row mask, which landed after this plan was written - see M1 item 3a) -
so renaming it would cost the parallel
that explains it. Therefore fork 3's mean-channel spelling does NOT block
reshape S1, provided the mean channel never asks for the word "weights", and
under this plan's spec it must not. What DOES remain schedule-bound, and is put
to VD as part of fork 3: if the precision channel is to be renamed at all
(`setForestRowPrecision`, `setForestCaseWeights`), that decision must land
before reshape S1. Recommendation: keep `setForestWeights`.
**Second finding, from the same anchors:** at reshape S1 the per-forest weight
channel becomes PUBLIC in C while remaining `dbarts:::`-only in R. Under binding
decision 1(a) that is backwards, and M1 closes it.

**B6 (the Linero quotation): ADOPTED.** The sentence *"we only require sharing
the tree topologies across model components, while the BCF model shares the
entire function h(x)"* is the source CONTRASTING shared forests with BCF, whose
`h(x)` is shared across potential outcomes; it is not the literature conceding
that shared-plus-group-specific decompositions ride a vector leaf. The correct
receipt for Linero's own model is its `psi` sentence
(`psi_{t,l}^(m)(x) == psi_{t,l}(x)`). This plan cites the class without the
contrast sentence, and the phrase "the literature says so in one sentence" is
dropped.

## Adjudication: corrected facts this plan is built on

Confirmed by the critique's probes and re-verified structurally here:

1. **`setForestWeights` is a PRECISION channel that `combinedFits` never
   reads.** `BCFForestCombiner::combinedFits` (`combiner.hpp:567-576`) reads
   `glue_.a`, `glue_.b0/b1` and `z`, never `forestWeights_`; the only consumer
   of the channel is `composeForestWeights`. Numerically: with `s = 0` on every
   row of forest 1, the reported location still contains that forest's full
   contribution to 1.8e-15. Any premise that the per-forest weight enters the
   combined fit is retired.
2. **dbarts already ships an R per-sweep callback, single-chain.**
   `bartcore_runWithCallback` refuses multi-chain outright
   ("bartcore_runWithCallback requires a single chain"), carries no
   `GetRNGstate`/`PutRNGstate` bracket by design (R owns `.Random.seed` because
   the closure draws from R's stream while the chain's generator does not), and
   evaluates the closure under `R_tryEval` so an error becomes a cooperative
   stop instead of a longjmp through C++ frames. "R callbacks are impossible" is
   false and must never be written.
3. **Rebuild plus `setState` silently drops a per-forest weight while the two
   stored states compare equal.** Real (reproduced: 0.38 max difference in tau
   fits, identical states) and already a RECORDED DECISION, not a discovery:
   `chain.hpp:962-965` documents it, `zero-weight-exactness.md:121-126` decides
   it, and `inst/tinytest/test-forest-weights.R` pins the same-holder round trip
   and records the cross-sampler case as documented rather than asserted. It is
   a contract item for M0's doc, not a defect.
4. **The non-bartcore header closure is 10, not 11, and `misc/simd.h` is not on
   it.** Recomputed here from the 11 bartcore headers: `external/io.h`,
   `external/random.h`, `external/stats.h`, `external/R.h`,
   `misc/linearAlgebra.h`, `misc/partition.h`, `misc/stats.h`, `misc/thread.h`,
   `misc/stddef.h`, `misc/types.h` (configure-generated from `types.h.in`).
   `misc/simd.h` appears only in `src/R_interface.cpp`, `src/misc/simd.c`,
   `src/misc/moments.c` and `src/misc/Makefile`.
5. **A pure-R Gaussian K-forest composition is a DIFFERENT PRIOR, and
   statistically close.** Measured (n = 400, p = 4, sigma fixed, BCF with its
   glue pinned so its model is exactly `y = mu + z tau`): engine RMSE(CATE)
   0.1066, R composition 0.0781 and 0.0773, pointwise max difference 0.187 on a
   CATE of about 2.0. The mechanism is the pinned response transform: each
   single-forest sampler's reported function carries the transform's `shift`, so
   the composition attributes `shift` K times (multiplier-weighted) where the
   engine attributes it once, and forest f's internal response leaves the
   [-0.5, 0.5] band the leaf prior is calibrated for by O(1/m_f). So the
   BLOCKING is identical; the leaf and amplitude priors are not. That residual
   sits exactly on the engine side of binding decision 1b's line: the additive
   decomposition ports to R perfectly, while the pinned response transform and
   the leaf calibration are prior shape, which is the part R-first does not
   carry. It is the sharpest available illustration of VD's qualification, and
   M0 documents it as such rather than as an anomaly.
6. **The honest cost numbers** are B3's: 0.89-0.95x engine per-sweep drive,
   1.39-1.43x R K-sampler composition at K = 2.
7. **Citations, corrected.** Ting and Linero: JASA 120(551):1400-1413, 2025,
   doi:10.1080/01621459.2025.2491155. Kim and Zigler: Biometrics 81(1), article
   ujaf024, 2025, doi:10.1093/biomtc/ujaf024. VCBART (Deshpande, Bai, Balocchi,
   Starling, Weiss): Bayesian Analysis 21(1):281-308, 2026 (the 2024-10-04 date
   is the advance publication). Hu et al., SMMR 29(11):3218-3234, 2020, is
   resolved and is a THIRD structure: one scalar-leaf ensemble with the arm as a
   covariate, neither K forests nor a vector leaf. Woody, Carvalho, Hahn, Murray
   (arXiv:2007.09845) is not merely venue-unverified, it is UNPUBLISHED (v1
   only, no journal reference) - load-bearing, because it is the only source for
   the continuous-multiplier class that fork 1 option (b) rests on. VCBART has
   NO amplitude: its leaf jumps are `N(0, tau_j^2)` with `tau_j` a fixed
   hyperparameter, never sampled.
8. **Two quote defects, not repeated here.** The `BART` package's
   `inst/cxx-ex/README.txt` says "This directory contains the BART source code
   and an example of how to create an executable for BART that does not require
   R"; the memo's "creating an executable independent of R" was a paraphrase
   inside quotation marks, and "independent of R" is itself overstated (the
   Makefile shells out to `R CMD config` and links `-lRmath`). LongBet's PAPER
   states no prior and no MCMC scheme for its time factor; the sampling facts
   belong to the reference implementation `google/longbet`
   (`src/mcmc_loop.cpp`, a conjugate draw once per sweep after both forest
   loops). This plan cites neither claim to the wrong artifact.

New findings made in synthesis, each a departure from the memo (receipts in
"Departures"):

9. **BCF is NOT an instance of the memo's declarative form, so the memo's
   spec and its falsifier F2 are both wrong.** `combinedFits` uses `b_{z_i}`
   and `drawGlue` draws `b0` and `b1` as two SEPARATE conjugate coefficients
   accumulated over the `z = 0` and `z = 1` subsets (`combiner.hpp:567-576`,
   `:609-622`). Since `b_{z_i} = b0 (1 - z_i) + b1 z_i`, forest 1 carries a
   TWO-COLUMN indicator basis with two coefficients, not one column times one
   amplitude. "Amplitude times one multiplier column" therefore does NOT
   contain BCF (it is BCF with `b0` pinned at 0), and the memo's F2 ("a K = 2
   declarative sampler with multiplier (1, z) is bitwise identical to a BCF
   sampler") is false as written. The channel must be a per-forest basis
   MATRIX. This is what "spec the actual channel it rides" resolves to, below.
10. **The general amplitude conditional is a small dense solve and the
    machinery ships.** `choleskyDecompose`, `solveLowerTriangular` and
    `solveLowerTriangularTransposed` (`model.hpp:869-895`) already serve
    `LinearGaussianLeaf`'s normal equations and `GPGaussianLeaf`'s kernel
    solves. BCF's indicator basis is ORTHOGONAL, which makes its Gram matrix
    diagonal, which is exactly why the shipped code draws two scalars.
11. **The memo's D2 is half moot.** The per-forest READERS were promoted at
    bcf-public-surface S2: `$getForestFits`, `$getBCFGlue` and
    `$getForestVariableCounts` are public R5 methods (`R/dbarts.R:1344`,
    `:1349`, `:1354`). Only the per-forest weight SETTER remains internal.
12. **Non-Gaussian multiforest is a refusal plus a calibration map, not missing
    machinery.** The latent draw already runs against the combined fit
    (`chain.hpp:1261-1262`) and the multinomial K-forest chain already runs a
    non-Gaussian response (its constructor sets `family_ = logistic` with
    Polya-Gamma augmentation). What forces Gaussian on the BCF path is the
    constructor itself: it builds a `GaussianResponse` and sets
    `family_ = ResponseFamily::gaussian` unconditionally (`chain.hpp:702-705`),
    and BCF's calibration map is stated through `scaledResponseSd()`, a
    range-scaled Gaussian notion. So the memo's "the non-Gaussian case is the
    family's irreducible engine content" STANDS, and its price is lower than
    "engine-only" implies: wire the family enum through the K-forest
    constructor and define the calibration map against each family's latent
    scale.
13. **The on-ramp's user-facing document already exists.**
    `vignettes/gibbs_sampler_mixture_model.Rmd` is titled "Building a Gibbs
    Sampler with dbarts" and opens *"For users who are comfortable with writing
    their own posterior samplers, `dbarts` makes it easy to incorporate a
    linear BART component in a Gibbs sampler"*; `man/dbartsSampler-class.Rd`
    carries an outer-Gibbs response-swap example. The memo's "no document
    states dbarts as a component as a supported contract" is half wrong: what
    is missing is the CONTRACT (mutation legality, state carriage, cost) and
    the K-forest recipe, not the introduction.

## The channel, specified

This is what the declarative recommendation rides. It is one channel, distinct
from the precision channel, and it contains BCF exactly.

**Model.** Each forest f carries a per-forest basis matrix `B_f` (n x q_f,
column-major) and a coefficient vector `a_f` of length q_f. The mean is

    E[y_i] = sum_f ( sum_k a_{f,k} B_f[i,k] ) f_f(x_i)

so the per-row effective multiplier `m_{f,i} = dot(a_f, B_f[i,])` stays a
SCALAR. `formForestResponse`'s reparameterization (response `r_i / m_{f,i}`,
weight `w_i m_{f,i}^2`) and the `0x1p-26` zero-multiplier snap
(`combiner.hpp:525-565`) are unchanged in shape, which is why this is a
generalization rather than a second mechanism.

**Instances.** BCF: forest 0 with an intercept basis and free amplitude `a`;
forest 1 with the two-column indicator basis `(1 - z, z)` and amplitudes
`(b0, b1)`. Continuous or dose exposure: forest 1 with the single column `z`.
VCBART: forest j with the single column `X_j`. Mediation: forests with columns
`1`, `a` and the observed mediator `m`. Principal stratification: forests with
column `a_i`.

**R spelling.** A numeric basis column enters linearly; a FACTOR basis expands
to level indicators with one coefficient per level. That is R's own convention
for a model matrix, and it is what makes BCF's `b0`/`b1` explicable as "the
coefficient vector of a factor basis":

    dbarts(y ~ x1 + x2, data = d,
           forests = list(forest(),
                          forest(basis = ~ factor(z), vars = c("x1", "x2"),
                                 n.trees = 50L, base = 0.25, power = 3)))

Every word is structural (`forests`, `forest`, `basis`, `vars`, `n.trees`,
`base`, `power`, `interactions`, `blocks`) and `forests = NULL` is exactly
today's single-forest path, so the argument is byte-neutral when absent (the
`variance =` precedent).

**Flat C spelling**, replacing the two BCF-specific names at the reshape
re-bake:

    X(int, dbarts_sampler_setForestBasis,
      (dbarts_sampler*, size_t forest, const double* basis, size_t numColumns),
      (sampler, forest, basis, numColumns))
    X(size_t, dbarts_sampler_numForestAmplitudes,
      (const dbarts_sampler*, size_t forest), (sampler, forest))
    X(int, dbarts_sampler_forestAmplitudes,
      (const dbarts_sampler*, size_t forest, double* out), (sampler, forest, out))

`basis` is column-major n x numColumns; `forestAmplitudes` fills
numForestAmplitudes(forest) x numChains; `1 = accepted, 0 = refused` (the house
convention). A raw pointer plus a count, never a struct, because the FNV-1a
token hashes signatures and is blind to struct layout
(`dbarts-h-reshape.md:198-207`, measured). **Ownership is COPY**, matching the
shipped `setTreatment` ("setTreatment copies the 0/1 values it is handed",
`dbarts.h:43`) and the bridge's `ownedForestWeights` pattern, and it must be
stated in the Doxygen: a continuous basis cannot be coerced-and-copied
incidentally the way a 0/1 `z` can, so the ownership sentence is load-bearing
rather than boilerplate.

**R5 spelling.** `$setForestBasis(forest, basis)` writing the engine AND the
`dbartsData` slot, on the mechanical mirror bcf-public-surface's state-carriage
decision installed for z; `$getForestAmplitudes(forest)`.

**Staged guards, which is what makes the signature shippable before the
family.** The body accepts, TODAY, only what today's engine honours: forest 1,
a two-column complementary 0/1 indicator basis, Gaussian family. Everything
else refuses at the entry, naming the capability. The family's arrival relaxes
guard bodies and moves NO header. This is the same posture the repo already
ships for the multi-forest guards ("unreachable today, guarded defensively") and
the same refuse-rather-than-drop discipline bcf-public-surface S1's F4 pinned.

**Two channels, forever.** The precision channel (`setForestWeights`,
`s_{f,i}`, forest f's own leaf conditionals) and the mean channel
(`setForestBasis`) stay separate entries. "Widen `setForestWeights` with a
basis" is the obvious wrong idea and is refused on the rule
zero-weight-exactness already applied to `setWeights` ("must NOT be widened
with a forest index, since the two channels differ on the sigma df").

**What this channel deliberately does NOT express** (binding decision 1b). It
carries DATA (basis columns) and the SELECTION of engine-provided priors per
forest: tree count, `base`, `power`, the leaf scale, the moderator column mask,
per-forest interactions and blocks, all of which `BCFForestSpec` already carries.
It carries no caller-authored prior and no caller-authored response shape: the
amplitude's prior law, the leaf model, the family and the link are the engine's,
selected by argument. So a spec asking for a prior form the engine does not
implement REFUSES at creation rather than being approximated - the F4 discipline
- and requests of that kind are engine arcs (a new leaf model, a new family, a
new hyperprior), never extensions of this channel. Concretely: BCF's asymmetric
amplitude prior (half-Cauchy on forest 0 via a scale-mixture auxiliary, normal
on the rest) is an engine choice with a default and a knob, not a slot a caller
fills with a law of their own; pinned rather than sampled amplitudes are a
DOOR, not a caller-supplied prior.

## Layering: three homes, and what replaces `treatment =`

**dbarts, the engine.** Owns: every response family (already a runtime enum);
every leaf model (five instantiations in `facade.hpp`); mean-structure
declarations (`blocks =`, `interactions =`, `monotone =`, `variance =`, and
`forests =`); the couplings a caller cannot compose (non-Gaussian latents
against the combined fit; the per-forest ASIS rescale, which writes leaf values
and for which no setter exists or should; exactness at near-zero multipliers;
one store, one cut grid, one state); the driver contract and its callback. Does
NOT own `bcf()`, `treatment =`, `moderators =`, or any estimand. Cost: M0 to M3
below, about 1630 lines pre-freeze, plus M4 if VD wants the family.

**stan4bart, the WALNUTS semiparametric collection (VD's question 2).** Its
boundary rule: a model belongs here when its non-BART component is PARAMETRIC
(a coefficient vector, random effects, a spline basis, a latent vector with a
gradient-friendly target) or when the coupling is through a DESIGN COLUMN
rather than a shared additive location. Two channels, both verified in code:
the additive offset (`setOffset` plus `setSigma` per iteration, stan4bart
`src/init.cpp`) and latent-covariate predictor mutation
(`updatePredictorPerObservationJointly` plus `$setPredictor(forceUpdate =
TRUE)`, bairrtt `R/irt_causal_bart.R`). What else fits, answering VD directly:
multilevel and random-effects BART (offset; `rbart_vi` is the in-engine
degenerate case); partially linear and semiparametric BART (offset; Zeldow, Lo
Re III, Roy, AOAS 13(3):1989-2010, 2019); CSP-BART and orthogonality-constrained
semiparametrics (offset; both cited through
`docs/design/tree-mixing-proposals.md` rather than re-fetched here); splines or
GAM plus BART (offset); AR-1 errors (offset, deferred in
`correlated-outcomes.md`); multivariate and SUR outcomes (ALREADY SHIPPED as
`mvbart()` with zero dbarts engine change, `docs/design/INDEX.md:35`, commit
e27a7c3); IRT plus BART (latent covariate; bairrtt, live); latent-confounder and
measurement-error sensitivity (latent covariate; the treatSens class);
sequential covariate imputation inside a fit (latent covariate); instrumental
variables (the first-stage fit feeds the second). Hurdle and two-part models
belong to NEITHER channel: the split is observed, so they are two conditionally
independent ordinary fits, shipped R-side with zero engine code
(`INDEX.md:31`). Cost: zero dbarts change for the offset channel; the
latent-covariate channel needs multiforest-predictor-mutation only if the BART
side becomes multi-forest. One correction that must travel with this home:
WALNUTS does NOT fix the cross-block ridge - it removes ridges INSIDE the
parametric block, and the parametric-versus-forest coupling remains an
alternating two-block Gibbs (`tree-mixing-proposals.md:558-572`, verified in
stan4bart's own loop).

**bartCause, the causal fit layer.** Owns `bcf()`, the estimands, common
support, the two counterfactual surfaces and the option-A names
(`treatment`, `moderators`, `mu`, `tau`, `glue`, `mu.hat.obs`, `mu.hat.cf`),
which stay correct for a fit function even though they leave the engine. Its
estimator layer is already written against two counterfactual SURFACES, and
under BCF both follow from mu, tau and glue with no test matrix at all
(`bcf-public-surface.md:503-510`; the `:479-484` this plan cited before the
2026-08-13 refresh now lands on S4's rng note and the S5 heading). Cost:
uncomment the reserved
`#bcf = redirectCall(...)` line in its response-method switch, write a `bcf`
response fitter against the S5 output contract, and add a moderator EXCLUSION so
its propensity-score-as-covariate option enters mu only. It has no `src/` and no
`LinkingTo`, so it can host no C++ - which is why the replacement creation route
must be a PUBLIC R route in dbarts.

**A third package: NO.** Sketch it as `y = X beta + sum_f a_f m_f f_f(x) + eps`
and it consumes only public dbarts primitives for the Gaussian case, degrades to
the R composition at 1.39-1.43x and a different prior (fact 5, not the memo's
1.04x), and needs three things that live only in `chain.hpp`/`combiner.hpp`: the
latent draw against the combined fit, the ASIS rescale that writes leaf values,
and the near-zero snap. Its Gaussian half collapses into being a stan4bart model
class (the `mvbart()` precedent adds a model class there with zero dbarts
change); its causal reporting layer is bartCause; its only engine content argues
for dbarts. Declined.

**What replaces `treatment =`: `forests =` / `forest()`, and this is DERIVED,
not chosen.** bartCause has no `src/` and no `LinkingTo`, so `bcf()` in its home
needs a public R creation route in dbarts; under binding decision 2 that route
must be engine vocabulary; the engine-vocabulary spelling of "the mean is a
weighted sum of K ensembles whose weights are data" is the basis declaration
above. Sequencing, and it satisfies VD's "removal only after a replacement
exists":

1. M2 lands `forests =` on `dbarts()`/`dbartsSpec()`, resolving to the SAME
   spec machinery `treatment =` resolves to, and REMOVES
   `treatment =`/`moderators =`/`treatmentForest =` in the same slice. No
   deprecation cycle: binding decision 5 (nothing released, lockstep
   migration, and bartCause has not adopted the argument).
2. M3 re-signs the flat names and the header's creation contract to match,
   inside the reshape's re-bake.
3. `bcf()` ships in bartCause over `forests =` (fork 4), and
   bcf-public-surface S5's `bartBCF` class relocates with it.
4. The R5 `$setTreatment` becomes `$setForestBasis`; `$getBCFGlue` becomes
   `$getForestAmplitudes`.

## Sequencing against the committed queue

    bcf-public-surface S4  (in flight, untouched)
    M0  on-ramp contract          docs + vignette + tests, no src
    M1  R/C parity on the weight channel   R5 + tests, no src
    multiforest-predictor-mutation S0-S4   (committed, untouched)
    M2  forests= replaces treatment=       R + tests, no src, bitwise-gated
    dbarts-h-reshape S0-S2, with M3 as an item inside S1   (one re-bake)
    M4  the general basis family           RESOLVED pre-release (fork 1,
                                            VD 2026-08-11); FA1/FA2 run
                                            first, design-informing
    1.0-0 freeze

M0 and M1 must land after bcf-public-surface S4 rather than beside it: S4 is
editing `inst/NEWS.Rd` and `man/dbartsSampler-class.Rd`, which both slices
touch. Nothing in M0 to M2 moves the API hash, so the window still contains
exactly two re-bakes (bcf S3, already landed, and the reshape's S1).

## M0. The on-ramp contract. Docs, vignette, pins. No src.

Under binding decision 1 this is the arc's headline deliverable, not an
ergonomics slice: it is the documented, tested stability story for the surface
model authors already use.

1. **`docs/design/bart-as-a-component.md`** (new), the internal contract. It
   COMMITS exactly three things, each already enforced in code:
   (i) a per-sweep driver loop is BITWISE identical to a batched run, at
   0.89-0.95x the batched cost;
   (ii) which mutations are legal between sweeps and under what predicate -
   `refuseMultiForestResponseMutation` (bridge and flat C share the one
   helper), `Chain::supportsResponseMutation` (a combiner that admits it AND
   `family_ == gaussian`) together with `updateScale == FALSE`,
   `refuseBCFTestSurface`, `refuseMultiForestMutation`,
   `refuseMultiForestTransactionalUpdate`;
   (iii) what the engine state does NOT carry - raw conditioning vectors (y,
   weights, offset, x, z, per-forest weights) as against derived quantities it
   cannot recompute (the quantized cut grid, the pinned response transform, tree
   structure, rng position) - and WHO reinstalls them: the R5 layer mirrors z
   through `data@treatment` and observation weights through `data@weights`, the
   per-forest weight rides an R5 field (M1), and a flat C consumer owns all of
   them.
   It COMMITS NOTHING about: per-forest decomposition of test fits
   (`testFitsAreDefined() == false` under BCF); mid-sweep hooks (refused on
   invariant grounds as well as threading grounds - during the tree loop
   "totalFits is stale until rebuilt after the loop", `chain.hpp:1183`);
   per-forest saved-tree replay; cross-host bitwise reproducibility
   (within-host across any SIMD dispatch only).
   It RECORDS the sweep-boundary hooks that exist and their constraints: the
   flat `dbarts_sampler_setCallback` (refused when `numThreads > 1 &&
   numChains > 1`; inline multi-chain runs chains sequentially) and the internal
   single-chain R hook (`bartcore_runWithCallback`, R owns `.Random.seed`,
   `R_tryEval` makes an error a cooperative stop).
2. **Extend `vignettes/gibbs_sampler_mixture_model.Rmd`** with a K-forest
   section: the R composition recipe (forest f gets
   `$setResponse(r_f / m_f, updateScale = FALSE)`, `$setWeights(w * m_f^2)`,
   `$setSigma`, `$run(0L, 1L)`; sigma drawn in R against the combined
   residual), its measured price (1.39-1.43x at K = 2, growing with K), and the
   five honest losses stated as prices rather than caveats: no ASIS ridge
   rescale (the in-house surrogate is a per-group INTERCEPT with a closed-form
   collapse, where alternating conditional blocks give IACT 56.1 against 9.3
   for a joint move at 3 groups with strong confounding, and the two are level
   at 10 groups with mild confounding, 13.1 against 11.7 - so it is an order of
   magnitude, not a forest-versus-forest measurement); the caller divides by
   `m_f` and so owns the
   `0x1p-26` condition-number cap itself; K quantized stores and K states; K
   leaf calibrations; a serial outer loop.
3. **The graduation section**, which is the principle's semantic-continuity
   answer AND the place binding decision 1b's boundary is stated to users: the
   additive decomposition ports exactly, and what does not port by itself is
   prior shape. The R recipe and the engine target the same posterior only if the
   caller (a) centres y at its midrange or absorbs the transform's `shift` with
   `setOffset`, because the composition attributes `shift` once per forest and
   the engine attributes it once; (b) pins the leaf prior identically across
   samplers (`node.prior = normal(k)`, and construct every sampler on a common
   y with `updateScale = FALSE`); (c) accepts that with a multiplier far from 1
   the internal response leaves the band the leaf prior is calibrated for. State
   plainly that without (a) and (b) the two routes are different priors that
   agree closely, and that with them the difference is the leaf calibration
   only - and that when `forests =` lands, the per-forest spec carries the same
   scale knobs, so a prototype ports without changing the model it targets.
   State equally plainly what the recipe does NOT reach and never will: a
   response family, link or latent structure the engine does not provide, and a
   leaf or hyperprior LAW the engine does not provide. Those are selected by
   argument, not composed in R, and adding one is an engine change.
4. **Also record the graduation-path debt** the flat surface still has relative
   to R5: no flat per-observation predictor update or joint session, no flat
   `setCutPoints`/`setData`, no flat `predictVariance`, no forest-indexed
   `predict`. All four are doors already recorded in `dbarts-h-reshape.md`;
   naming them here is what makes them on-ramp debt rather than trivia.

Falsifiers: **FD1** the bitwise identity, both halves - extend
bcf-public-surface S0's BCF pin to a SINGLE-forest sampler, and a deliberately
mismatched thinning must go red. **FD2** a refusal matrix asserting (ii) entry
by entry, R5 and flat C. **FD3** the rebuild divergence of fact 3 pinned as
DOCUMENTED behaviour with its receipt in the test comment. **FD4, the
continuity falsifier, runnable today with no engine work**: a BCF sampler with
`update.a = update.b = FALSE` (glue pinned, model exactly `y = mu + z tau`) and
the R recipe with the section-3 remedies agree on CATE to Monte Carlo error;
NEGATIVE HALF, mandatory: without the remedies the difference is OBSERVABLE
(about 9% pointwise at the critique's configuration), which is what makes the
documented price honest rather than reassuring.
rng: NEUTRAL (no src). Gates: full tinytest; `R CMD check` so the vignette
builds; `air format --check .`; `lintr` on touched R.

## M1. R/C parity on the per-forest weight channel. No src.

`$setForestWeights(forest, weights)` on `dbartsSampler`, with the R5 field
re-applied after re-creation (the mirroring rule bcf-public-surface recorded for
this channel), Rd naming BOTH channels and cross-referencing them, and the
rebuild hazard documented at the method. Reason it is not optional: at reshape
S1 the flat entry becomes public, and a channel reachable from C but not from R
inverts binding decision 1(a). It also gives the on-ramp its
zero-information-row affordance (annealed or partial membership) without a
`dbarts:::` reach-in.

**Anchors, re-verified BY SYMBOL at 934a02d5 (added 2026-08-13, pre-M1).**
Engine: `Chain::setForestWeights` `chain.hpp:966` with its semantics doc
`:944-965` (the composition sentence `:948-950`, the rebuild-drop sentence
`:962-965`); `Chain::supportsForestWeights` `:878-880`;
`Chain::composeForestWeights` `:3431-3447`, fn `:3439`, called from the forest
loop at `:1160`; the capability is true only on `BCFForestCombiner`
(`combiner.hpp:738`). Fan-out `sampler.hpp:1170`; facade virtual
`facade.hpp:279`, shape field `:80`, populated `:358`. Bridge
`bartcore_setForestWeights` `R_interface_bartcore.cpp:3688-3728`. Internal R
wrapper `bartcoreSetForestWeights` `R/bartcore.R:999` (its semantics comment
`:984-998`), absent from `NAMESPACE`, with no R5 method. The R5 method's site
is between `$setActiveRows` (`R/dbarts.R:1109`) and `$setTreatment` (`:1132`).
Rd: aliases `man/dbartsSampler-class.Rd:18-19`, usage `:62-63`,
`\item{weights}` `:144-146`, `\item{active}` `:147-161`, `\item{forest}`
`:165-167`. The ENGINE channel is already pinned by
`inst/tinytest/test-forest-weights.R` (293 lines), including the all-ones
identity (`:73-79`), `storeState`/`setState` and `installForests` round trips
(`:276`, `:292`).

**1. Forest indexing: 1-BASED, converted at the boundary. SETTLED (VD
2026-08-13; added 2026-08-13, pre-M1.)** R code holds 1-based throughout, with
ONE conversion point at the bridge call site - the `resolveForestIndex` pattern
(`R/bartcore.R:1051`, which rejects anything below 1 and returns `forest - 1L`).
The C API, the bridge and the engine stay 0-based. So `$setForestWeights(forest,
weights)` routes through `resolveForestIndex`, NOT through the raw 0-based
internal wrapper `bartcoreSetForestWeights` (which is 0-based, takes a bartcore
handle rather than an R5 pointer, and is therefore not reusable here anyway):
the method calls `C_dbarts_bartcore_setForestWeights` directly with the
converted index, as `$getCalibration`/`$setCalibration` already do
(`R/dbarts.R:1362`, `:1393`). A BCF sampler's treatment forest is therefore
`2L` in every M1 example, Rd sentence and test. Relation to the open
`r5-forest-indexing` ticket (`TODO:413-424`): M1 MUST NOT extend the 0-based
split onto a new axis - it lands on the 1-based side, leaving the ticket with
exactly the two 0-based getters it already names (`$getForestFits`,
`$getForestVariableCounts`), which the dbarts-h-reshape re-bake still
normalizes. The Rd's `\item{forest}` paragraph (`:165-167`) documents both
conventions today and gains a third clause naming `setForestWeights` on the
1-based side.

**2. The mirroring mechanism: a NEW R5 field, re-applied at BOTH re-creation
sites (added 2026-08-13, pre-M1).** `$setTreatment` is not the template it
looks like. It mirrors through `data@treatment` (`R/dbarts.R:1144-1146`), a
`dbartsData` slot (`R/A_class.R:518`) that `bartcore_create` reads back for
free on every re-creation. `dbartsData` has NO per-forest weight slot
(`R/A_class.R:478-541`, verified slot by slot), and the R5 fields are exactly
six - `pointer`, `control`, `model`, `data`, `state`, `hostFor`
(`R/dbarts.R:711-722`) - none of which can carry one. So M1 must ADD an R5
field holding the installed per-forest weights (a list or a matrix indexed by
1-based forest, `NULL` where none is installed) and RE-APPLY it after every
re-creation, at BOTH live sites:
   - `getPointer()` (`R/dbarts.R:1428-1456`), after its
     `C_dbarts_bartcore_create` + `C_dbarts_bartcore_setState` pair;
   - `setState()` (`R/dbarts.R:1458-1482`), after its own
     `C_dbarts_bartcore_setState`, which re-creates on an invalid pointer by
     the same two calls.
   Re-application must not recurse through `getPointer()`. FW1 as the plan
   originally wrote it exercises only the first site; **FW1 is extended to
   exercise BOTH**, each with the field-removed negative half. `$setActiveRows`
   is NOT a counter-precedent to mirroring: it installs no mirror at all and
   the Rd says so (`man/dbartsSampler-class.Rd:158`, the mask does not ride a
   saved state), which is precisely the asymmetry M1 is closing for this
   channel under bcf-public-surface's recorded mirroring rule.

**3. Interactions with the latent-subset-mask arc, which this plan predates
(added 2026-08-13, pre-M1).**
   (a) **Composition is now `(w_i a_i) m_f^2 s_i`, not `w_i m_f^2 s_i`.** The
   active-row mask lands inside the response family: `workingWeights()` serves
   `w * a` while a mask is installed and `w` by identity when it is not
   (`docs/plans/latent-subset-mask.md`, "Composition with setWeights",
   `:436-448`), the chain reads it at `chain.hpp:1101`, the combiner forms
   `(w a) m_f^2` in `formForestResponse` (`combiner.hpp:562`), and
   `composeForestWeights` multiplies `s` on top afterwards (`chain.hpp:1160`,
   `:3439`). The two channels are independent and order-free (the Rd already
   states `w_i * a_i` for the mask, `man/dbartsSampler-class.Rd:148`). M1's Rd
   states the FOUR-factor form. The internal comments still say the three-factor
   `w_i m_f^2 s_i` - `chain.hpp:948-950` and `R/bartcore.R:984-998`; refreshing
   the R one is free for M1, the `chain.hpp` one is comment-only src and rides a
   later slice rather than breaking M1's no-src rule.
   (b) **The Rd's degenerate-value contrast becomes THREE-way.** Verified:
   `setWeights` all-ones INSTALLS and is measurably distinct from carrying no
   weights at all; `setActiveRows` all-ones reports success, installs NOTHING
   and CLEARS any installed mask (both stated at
   `man/dbartsSampler-class.Rd:148`); `setForestWeights` all-ones INSTALLS -
   it takes the composed path (scratch buffer, an O(n) multiply) where not
   installing passes the combiner's own pointer through - and is BITWISE
   identity, pinned as the identity arm at
   `inst/tinytest/test-forest-weights.R:73-79`, while staying in force across a
   state round trip (`:276`, `:292`). M1's Rd must state all three, and must
   not flatten the third into "no-op": it installs, and a later round trip can
   tell.
   (c) **`setData` vs an installed per-forest weight: the interaction is
   UNREACHABLE, and that is the answer M1's Rd documents.** `setData` clears
   the active-row mask because the mask is length-n and `setData` may change n.
   The per-forest weight needs no such rule: `bartcore_setData` refuses on ANY
   multi-forest sampler before it does anything (`R_interface_bartcore.cpp:4337`
   calling `refuseMultiForestMutation`, whose body is at `:2503-2508` - "a
   multi-forest sampler fixes its data at creation; make a new sampler
   instead"), and `supportsForestWeights` is true only on `BCFForestCombiner`
   (`combiner.hpp:738`), i.e. only on a sampler with two forests. So a sampler
   that can HOLD a per-forest weight can never TAKE a `setData`, and nothing
   clears anything. The Rd says the channel is BCF-only and that `setData` is
   refused on such a sampler, rather than repeating the mask's "setData clears"
   sentence. DOOR: if a later slice makes `setData` reachable on a multi-forest
   sampler, that slice owns the clearing rule for both the engine pointer and
   this R5 field.

**4. The multinomial refusal sits INSIDE the entry M1 wraps; the R5 posture is
PASS-THROUGH (added 2026-08-13, pre-M1).** Mask S3 put the softmax-margin
refusal in `bartcore_setForestWeights` itself
(`R_interface_bartcore.cpp:3694-3704`: the capability probe comes first, and
under `supportsCountsMutation` it raises "a multinomial mask applies to every
category: the softmax margin reads all K forests, so a row cannot leave one
category's likelihood alone" before falling through to
"bartcore_setForestWeights requires a BCF sampler"). POSTURE: the R5 method
adds NO pre-refusal of its own and lets the bridge raise. This matches the
house posture of its neighbours - `$setActiveRows` (`R/dbarts.R:1109-1131`)
validates shape and length in R and leaves every family/capability refusal to
the bridge, and `$setCalibration`'s only R-side pre-refusal
(`refuseBCFMutation`, `:1373-1378`) is a BCF-calibration refusal, not a family
one. Like every R5 mutator, the method opens with
`refuseHostMutation("$setForestWeights")` (`:770-783`) and validates length,
`NA` and non-negativity in R before the call, on the `$setWeights` pattern.
   **One correction, verified rather than assumed: the multinomial message is
   NOT reachable through the public R5 route at HEAD, so M1 must not budget a
   test that pins it there.** A multinomial engine only
   ever exists as a bare-environment bartcore handle built by
   `C_dbarts_bartcore_createMultinomial` (`R/bartcore.R:797-806`); the
   `dbartsSampler` a multinomial fit carries is a HOST SHELL over an ordinary
   single-forest engine (`R/bart.R:1258`, `:1336`), every mutator on which is
   refused first by `refuseHostMutation`. The R5 route's reachable refusals are
   therefore the host-shell one and, on an ordinary sampler,
   "bartcore_setForestWeights requires a BCF sampler". M1 pins those two
   through the PUBLIC route; the multinomial message keeps its existing
   internal-route pin at `inst/tinytest/test-active-rows-pins.R:501-504`
   (the call at `:502`), which STAYS. If a later slice gives the R5 layer a
   multinomial-shaped pointer, the public pin lands there.

**5. Budget, RESTATED with the carried items included (added 2026-08-13,
pre-M1; supersedes the ~60 R + ~40 man + ~140 test committed 2026-08-11, and
carried into this plan's header).**
~270-290 dense-equivalent total: **R ~75-85** (the method, plus the new field
and its re-application at both re-creation sites), **man ~50-60** (the method's
own entry, the three-way degenerate-value contrast, the `setData` clause, the
1-based `\item{forest}` clause, the rebuild hazard, the two-channel
cross-reference), **test ~140 sited on the R5 LAYER** - the mirror at BOTH
re-creation sites, the R5-route refusals, FW2, and the 1-based boundary
conversion (that `2L` reaches forest 1 and that `0L` is refused) - and NOT on
re-pinning the engine channel, which `test-forest-weights.R` already covers end
to end. The calibration-S2 lesson applies: this budget moves with the
obligations, and a stale one beside added obligations is the failure mode.

Falsifiers: **FW1** the mirror, both halves, at BOTH re-creation sites - a
`storeState`/`getPointer()` re-creation and a `setState()` re-creation each
keep the weight in force, and with the field removed the difference is
observable in each. **FW2** the channel does not move the mean: with
`s = 0` on every row of one forest, the reported location still contains that
forest's full contribution (the critique's 1.8e-15 probe promoted to a shipped
test), so a reader cannot mistake it for an exclusion from the fit.
rng: NEUTRAL (no src). Gates: `R CMD INSTALL` into a private library; full
tinytest with NO snapshot regenerated; `saveRDS`/`readRDS` round trip;
`R CMD check`; `air format --check .`; `lintr` on touched R; pkgdown check,
with NO new Rd topic and NO `_pkgdown.yml` edit - the alias rides
`dbartsSampler-class`, already listed at `_pkgdown.yml:26`. The trio is
expected IDENTICAL (no src), against baselines `equivalence-8b047f8b` (37
scenarios), `bcf-equivalence-8b047f8b` (12) and
`multinomial-equivalence-1027be5` (10).

## M2. `forests =` replaces `treatment =`. R only, bitwise-gated.

1. `forests = list(forest(...), ...)` on `dbarts()` and `dbartsSpec()`, with
   `forest(basis =, vars =, n.trees =, base =, power =, interactions =,
   blocks =, sd =)` as the single constructor (the
   `interactions()`/`blocks()`/`treatmentForest()` precedent: fourteen knobs
   ride one constructor, so `dbarts()` grows exactly ONE argument, down from
   three).
2. `resolveSamplerSpec` resolves `forests =` to the SAME
   `attr(control, "bartcore.bcf")` payload plus `data@treatment` that
   `treatment =` resolves to today, expanding a factor basis to its level
   indicators in R so the bridge sees exactly what it sees now.
3. `treatment =`, `moderators =`, `treatmentForest =` are REMOVED, and
   `$setTreatment`/`$getBCFGlue` are renamed to
   `$setForestBasis`/`$getForestAmplitudes` with the old spellings gone.
   `inst/tinytest/test-bcf-creation.R` is rewritten; the other ten files
   touching `treatment` have their creation calls updated.
4. Every declaration today's engine cannot honour REFUSES at creation, one
   assertion each: K > 2 forests; a numeric (non-factor) basis; a factor with
   more than two levels; a basis on forest 0; a non-Gaussian family; a variance
   forest alongside; plus the eleven S1 refusals inherited unchanged.
5. Docs: `docs/design/bcf.md` and `docs/design/public-surface.md` gain the
   engine-vocabulary surface; `inst/NEWS.Rd`; `_pkgdown.yml` if a new topic
   appears (the recorded new-exported-Rd-topic lesson).

Falsifiers: **FS1, load-bearing** - with `control@rngSeed` SET, a `forests =`
BCF is BITWISE IDENTICAL to the removed `treatment =` route (reconstructed from
the pinned tests) on all six bcf-equivalence channels; NEGATIVE HALF: perturb
the factor expansion (swap the two indicator columns) and it must go red.
**FS2** the refusal list of item 4, one assertion each, each raising an R
condition with no handle escaping. **FS3** `forests = NULL` is byte-neutral: a
single-forest fit is bitwise unchanged.
rng: NEUTRAL on every existing path. Gates: `R CMD INSTALL` into a private
library; full tinytest with NO snapshot regenerated; the trio BITWISE; air;
lintr; pkgdown check. ABORT: any trio divergence.

## M3. The flat mean channel. An item INSIDE dbarts-h-reshape S1.

Not a separate slice and not a separate re-bake. Conditional on fork 3 landing
before that slice starts; if it has not, this item is dropped and the shipped
BCF names go to release.

1. Re-sign `dbarts_sampler_setTreatment` as `dbarts_sampler_setForestBasis`
   with the signature and staged guards specified above; the guard body accepts
   forest 1, a two-column complementary 0/1 basis, Gaussian family, and refuses
   the rest naming the capability.
2. Re-sign `dbarts_sampler_bcfGlue` as `dbarts_sampler_numForestAmplitudes` plus
   `dbarts_sampler_forestAmplitudes`, so a reader can size its own buffer for a
   ragged amplitude vector (BCF: one for forest 0, two for forest 1).
3. Re-word the creation Doxygen at `dbarts.h:348-357` and the ownership
   sentence at `:43` in engine vocabulary, and state the COPY contract for a
   basis explicitly.
4. `consumer.c` and `test-capi.R` gain the accept and refusal legs; the hash is
   re-baked once, together with the reshape's own edits, and both version
   constants stay 1 and 0.

Falsifiers: **FC1** a flat-created sampler driven through `setForestBasis`
reproduces the R5 path bitwise at a shared seed. **FC2** the refusal matrix run
from `consumer.c`, outcome AND message text checked in C (the S3 harness
pattern). **FC3** the hash moved and neither version constant did.

## M4. The general basis family. RESOLVED pre-release (VD 2026-08-11, fork
1), scheduled after the dbarts.h reshape.

VD adopted fork 1: BUILD, pre-release, not post-freeze as recommended. FA1
and FA2 still run BEFORE any engine work, now as design-informing probes
rather than go/no-go gates - the orchestrator's reading of the decision,
stated to VD and uncorrected - because they can still shape the amplitude
and ASIS design even though a null result no longer collapses this slice
into option (c). Falsifier FA5 stands as a required respec at scheduling
time: its committed form below tests a strawman (K independent probit
samplers); the decisive arm is K GAUSSIAN samplers with host-drawn latents
against the combined fit, which measurement predicts will AGREE.

- **M4.0 (tests only).** Pin today's combiner at the seam that generalizes:
  `forestMultiplier`'s two values, `combinedFits`' blend, `drawGlue`'s
  conditional order (a, aVariance, b0, b1, then the ridge draw), and
  `afterCombine`'s applied scale.
- **M4.1 (engine, byte-identical).** `forestMultiplier` becomes
  `dot(a_f, B_f[i,])`; `combinedFits` becomes a K loop. Ship with BCF's basis
  synthesized internally from z so BCF is bitwise unchanged. This is the whole
  architectural risk, isolated. ABORT on any bcf-equivalence divergence.
- **M4.2 (engine).** The q-variate amplitude conditional, using the shipped
  Cholesky helpers; `afterCombine` becomes a per-forest ASIS rescale of the
  whole amplitude vector along the likelihood-invariant orbit. BCF's K = 2
  orthogonal case must come out bitwise, or its specialized two-scalar path is
  kept explicitly and the general path is gated statistically. This is the one
  place in the arc that could force a bcf-equivalence re-record; say so in the
  slice's rng note rather than discovering it.
- **M4.3 (spec, factory, refusal relaxations).** A K-length spec vector of
  `BCFForestSpec`-shaped entries plus per-forest basis and amplitude priors;
  `BCFSpec` becomes a thin adapter so `bartcoreBCFSampler` and the fixture keep
  working; the M2 and M3 guards relax in lockstep, with no header movement.
- **M4.4 (non-Gaussian, the family's justification).** Wire the family enum
  through the K-forest constructor and define the calibration map against each
  family's latent scale; probit and logistic in v1, the rest doors.
- **M4.5 (docs).** `docs/design/multiplier-combiner.md`;
  `forest-combiner.md`'s "What still re-carves" updated;
  `model-space-survey.md` gains the four verified classes it lacks.

Falsifiers: **FA0** a K = 2 factor-basis sampler is bitwise identical to
today's BCF at the same seed (this REPLACES the memo's F2, which fact 9
falsifies). **FA1, the amplitude falsifier** (the critique's A3, adopted): a
K-forest family WITHOUT per-forest amplitudes, in VCBART's own
parameterization, against one WITH amplitudes plus the ASIS move, same DGP,
IACT and RMSE both reported; if the no-amplitude arm mixes comparably, the
amplitude and its ASIS remedy are decoration and the Gaussian half has no
engine justification. **FA2** with the per-forest ASIS rescale removed, the
amplitude's IACT degrades measurably; if it cannot go red, the remedy is
decoration. **FA3** a continuous single-column basis recovers a known
closed-form posterior to Monte Carlo error, in `bcf-exact.R`'s idiom. **FA4**
a K = p+1 VCBART-shaped sampler recovers known coefficient surfaces better than
one `LinearGaussianLeaf` fit on a DGP where the coefficients break at
DIFFERENT modifier thresholds; a null means this should have been a leaf model.
**FA5** the non-Gaussian claim, MEASURED not asserted: a probit K-forest
sampler draws its latents against the combined fit, and K independent probit
samplers composed in R are a measurably DIFFERENT model on the same DGP.
**FA6** every refused creation raises an R condition and no external pointer
escapes.

Priced honestly, and higher than the memo's estimate because of fact 9: ~450 to
650 engine lines, ~200 bridge, ~250 R, plus the M4.0 pins and one possible
bcf-equivalence re-record. Sequenced after multiforest-predictor-mutation S0-S4
so the window holds one re-record, and after the dbarts.h reshape (RESOLVED
pre-release, VD 2026-08-11 - not after the freeze as this plan recommended).

## Cross-plan amendments (apply verbatim)

**1. `docs/plans/dbarts-h-reshape.md`, "Binding decisions inherited (do not
reopen)", item 3.** Replace the item with:

> 3. **bcf-public-surface S3's entries are adopted VERBATIM, with ONE carve-out
>    (amended 2026-08-10 by multiforest-extension-surface, which adjudicated a
>    conflict between this item and bcf-public-surface's AMENDED Open decision
>    block).** The widened `dbarts_sampler_setResponse(sampler, y,
>    updateScale)`, `numForests` and `forestFits` are adopted verbatim and
>    re-signed never; no guard S3 relaxed is re-relaxed. CARVE-OUT: the two
>    BCF-SPECIFIC names `setTreatment` and `bcfGlue`, and the creation contract
>    documented at `dbarts.h:348-357` plus the ownership sentence at `:43`, MAY
>    be re-signed to engine vocabulary in this arc's S1, because
>    bcf-public-surface's later VD amendment authorizes exactly that ("The flat
>    C names S3 shipped (`setTreatment`, `bcfGlue`) are likewise renameable at
>    the queued dbarts.h reshape re-bake if the settled surface uses engine
>    vocabulary") and because this arc's S1 is the window's LAST re-bake, so
>    declining costs a post-release MINOR bump plus a lockstep bump of
>    stan4bart's floor - the exact cost resolved question 2 refused to pay for
>    `setForestWeights`. The re-sign is CONDITIONAL: it happens only if
>    multiforest-extension-surface's naming fork is answered before this slice
>    starts, and its content is that plan's M3. Otherwise the shipped names go
>    to release unchanged.

**2. `docs/plans/dbarts-h-reshape.md`, S1, as a new item 5b after item 5.**
Insert:

> 5b. **The mean channel, CONDITIONAL on multiforest-extension-surface's fork 3
>    (see binding decision 3's carve-out).** Re-sign
>    `dbarts_sampler_setTreatment` as
>    `X(int, dbarts_sampler_setForestBasis, (dbarts_sampler*, size_t forest,
>    const double* basis, size_t numColumns), (sampler, forest, basis,
>    numColumns))`, `basis` column-major n x numColumns; replace
>    `dbarts_sampler_bcfGlue` with
>    `X(size_t, dbarts_sampler_numForestAmplitudes, (const dbarts_sampler*,
>    size_t forest), (sampler, forest))` and
>    `X(int, dbarts_sampler_forestAmplitudes, (const dbarts_sampler*, size_t
>    forest, double* out), (sampler, forest, out))` filling
>    numForestAmplitudes(forest) x numChains. 1 = accepted, 0 = refused. The
>    body accepts only what today's engine honours - forest 1, a two-column
>    complementary 0/1 basis, Gaussian family - and refuses the rest naming the
>    capability, so the family relaxes guard bodies later and moves no header.
>    Ownership: the entry COPIES, matching `setTreatment` (`dbarts.h:43`); state
>    it in the Doxygen, because a continuous basis cannot be coerced-and-copied
>    incidentally the way a 0/1 z can. Re-word the creation Doxygen
>    (`dbarts.h:348-357`) in engine vocabulary. Budget ~40 header + ~80
>    C_interface + ~70 consumer.c + ~60 test-capi.R.

**3. `docs/plans/dbarts-h-reshape.md`, S1 item 5 (`setForestWeights`).** Append
one sentence:

> The precision channel and the mean channel stay TWO entries forever: `s_{f,i}`
> scales forest f's own leaf conditionals and never enters `combinedFits`
> (measured, 1.8e-15), while a basis scales forest f's contribution to the mean.
> Widening this entry with a basis is refused on the rule
> zero-weight-exactness applied to `dbarts_sampler_setWeights`.

**4. `docs/plans/bcf-public-surface.md`, Open decision, at the end of the
AMENDED block.** Append:

> **RESOLVED, in part, by multiforest-extension-surface (2026-08-10).** The
> replacement creation route is `forests = list(forest(basis = ...), ...)` on
> `dbarts()`/`dbartsSpec()`, resolving to THIS arc's spec machinery unchanged;
> `treatment =`, `moderators =` and `treatmentForest =` are REMOVED in the same
> slice rather than deprecated, since nothing is released and bartCause has not
> adopted them. The rename surface is THREE things, not two: this block names
> `setTreatment` and `bcfGlue`, and misses the creation contract documented at
> `dbarts.h:348-357` ("dbartsSpec(data, control, treatment = z)") plus the
> ownership sentence at `:43`. The conflict between this block and
> dbarts-h-reshape's binding decision 3 is adjudicated in favour of this block:
> that decision now carries a conditional carve-out for the two BCF-specific
> names and the creation contract.

**5. `docs/plans/bcf-public-surface.md`, S5.** Prepend to the slice:

> **PAUSED, and RELOCATED under multiforest-extension-surface's home decision
> (conditional on VD's fork 4).** `bcf()` and the `bartBCF` class ship in
> bartCause, not in dbarts: the engine exports no fit function with domain
> semantics, and bartCause already owns the estimands, common support and the
> two-counterfactual-surface output contract this slice was written against.
> What survives verbatim is the CONTRACT - per-draw mu, tau, glue, sigma,
> per-forest varcount and both counterfactual surfaces, forest-INDEXED rather
> than positional, with the option-A element names, which stay correct for a
> fit function. If VD keeps `bcf()` in dbarts, this slice resumes as written.

**6. `docs/plans/bcf-public-surface.md`, "Doors held open", the "Public
`setForestWeights`" line.** Replace with:

> - **Public `setForestWeights`** - SCHEDULED as multiforest-extension-surface
>   M1, which must land before dbarts-h-reshape S1 makes the flat entry public;
>   a channel reachable from C but not from R inverts VD's prototyping
>   principle. It inherits the mirroring rule on exposure and carries the
>   documented rebuild hazard at the method.

**7. `docs/plans/multiforest-predictor-mutation.md`, "Binding decisions
inherited (do not reopen)", after item 5.** Insert:

> 5b. **The large-K price is UNPRICED and is recorded, not assumed away**
>    (multiforest-extension-surface, 2026-08-10). The veto rate scales with the
>    j-splitting trees summed over ENSEMBLES, and the per-forest column mask
>    does not save a VCBART-shaped model, where all p+1 ensembles are functions
>    of the SAME modifiers, so the mask is vacuous there rather than an opt-out.
>    No measurement exists above K = 4. If a general K-forest family lands, that
>    arc owns a fresh veto-rate measurement at its own K.

**8. `TODO`, the `multiforest-extension-surface` door sentences.** APPLIED,
records commit (2026-08-11). The sentences originally cited at :86 and :625
were already consolidated to one block at TODO:124-144 by e2cc1de (past this
item's own literal replacement text - it named the completed R/C-division
review); the records commit rewrites that block again, then at TODO:154-178
and at 934a02d5 at TODO:190-214,
recording: plan path `docs/plans/multiforest-extension-surface.md`,
artifacts under `.claude/multiforest-extension-surface-design/`, the answer
(the driver surface is the shape and is nearly shipped; the declarative
per-forest basis is the model-carrying surface and its justification is the
non-Gaussian case), the slices M0 to M3 pre-freeze plus M4 (the general basis
family) ALSO pre-release, scheduled after the dbarts.h reshape (not
post-freeze as this item originally proposed), the three-home layering
(dbarts engine, stan4bart WALNUTS collection, bartCause causal fit layer, no
third package), VD's R-first principle with its scope limit (R-first for
additive composition and for configuring engine-provided families and
priors; authoring a response shape or a prior law is an engine change with
C++ as its path), and the four VD forks, ALL RESOLVED 2026-08-11 (see "Open
decisions" below).

## Open decisions (VD)

**Fork 1. Build the general per-forest basis family, or stop at the driver
surface?**

The family is the declarative surface specified above: each forest's
contribution to the mean is scaled by caller-supplied data columns with their
own coefficients, so BCF, dose-response BCF, VCBART's varying coefficients,
heterogeneous mediation and principal stratification are all one engine with
different basis columns.

- **(a) Build it, post-freeze, entry-gated (RECOMMENDED).** Cost ~900 to 1100
  lines across engine, bridge, R and tests, plus M4.0's pins and possibly one
  bcf-equivalence re-record. What it buys: five published classes (BCF,
  Bayesian Analysis 2020; continuous exposure, unpublished preprint; VCBART,
  Bayesian Analysis 21(1), 2026; mediation, JASA 120(551), 2025; principal
  stratification, Biometrics 81(1), 2025); the demand datum that a model writer
  with a two-forest need REIMPLEMENTED a BART engine rather than reuse one
  (Sun and Song's AFT mixture-cure model, whose README says its tree code was
  "constructed following the framework of R package bcf"); and, decisively under
  the prototyping principle as qualified by binding decision 1b, it keeps
  NON-GAUSSIAN multiforest models declarable in R without crossing the boundary:
  the author supplies basis COLUMNS as data against a family the engine already
  provides, which is composition, not response-shape authorship. The gate: run
  FA1 and FA2 first. They are cheap, and a null on both plus a Gaussian-only
  scope makes (c) the right answer on evidence.
- **(b) Continuous single-column basis only (dose-response BCF).** ~150 to 200
  lines on top of M4.1; the numerics already ship (the `0x1p-26` snap and its
  condition-number argument). Cost of choosing it: the only source for the class
  is an UNPUBLISHED preprint, and the general family then arrives as a second
  refactor over the same code.
- **(c) Decline the family. Ship M0 to M3 and stop.** Coherent, cheapest,
  serves every consumer that exists today, and honest: the Gaussian half is
  measured ergonomics (1.39-1.43x, and a documented prior difference). Cost: the
  five published classes, and non-Gaussian multiforest models stay unreachable
  from R at all - a K-forest probit model would then require authoring engine
  code for what is, on the boundary of binding decision 1b, a composition
  question.

**Recommendation: (a), entry-gated on FA1 and FA2, with its justification
stated as the non-Gaussian case rather than Gaussian ergonomics.** What would
change it: FA1 coming back null AND a decision that non-Gaussian multiforest is
out of scope, which collapses (a) into (c).

**RESOLVED (VD 2026-08-11): (a), BUILD the family - PRE-RELEASE, not
post-freeze as recommended.** Scheduling within the backlog at orchestrator
discretion -> placed after the dbarts.h reshape (M3 carries its header
entries into reshape S1). FA1 and FA2 still run first, but as
DESIGN-INFORMING probes rather than go/no-go gates (orchestrator
interpretation, stated to VD, uncorrected) - a null result now shapes the
amplitude/ASIS design rather than collapsing this slice into option (c).

**Fork 2. Does the family's v1 cover non-Gaussian responses?**

This is the axis that decides whether the family has engine content at all. For
Gaussian there are no latents to refresh, so the composition is expressible in R
today; for probit, logistic, Student-t, ordinal, negative-binomial and the
softmax the latents are drawn ONCE against the COMBINED fit
(`chain.hpp:1261-1262`), so K independent R samplers are a DIFFERENT model, not
a slower one.

- **(a) Probit and logistic in v1, the rest doors (RECOMMENDED).** Cost is lower
  than "engine-only" implies (fact 12): the latent machinery already runs
  against the combined fit and the multinomial K-forest chain already runs a
  non-Gaussian response; the work is wiring the family enum through the K-forest
  constructor, which today hardcodes `GaussianResponse` and
  `family_ = gaussian`, and defining the calibration map against the latent
  scale in place of `scaledResponseSd()`. Probit and logistic carry the biggest
  multiforest literature (principal stratification's intermediate model,
  mediation with a binary mediator).
- **(b) Gaussian-only v1, documented as such.** Cheaper by the calibration map
  and its gates. Cost: `family = "probit"` with `forests =` must REFUSE, and a
  user will expect it to work; and the family then ships with no irreducible
  engine content, which is exactly the line fork 1's option (c) is drawn on.

**Recommendation: (a).** What would change it: a decision that the calibration
map for latent families needs its own research arc, in which case (b) ships
first and (a) becomes the second slice.

**RESOLVED (VD 2026-08-11): (a), probit and logistic in v1, other families as
doors** - ADOPTED as recommended, with the escape hatch kept live: if the
latent-scale calibration map needs its own research arc, Gaussian-only ships
first and non-Gaussian is the immediate second slice, both still
pre-release.

**Fork 3. The flat naming, and which committed plan yields.**

Two committed plans disagree about who may rename the BCF-specific flat entries:
`dbarts-h-reshape.md` binding decision 3 says they are adopted verbatim and
re-signed never; `bcf-public-surface.md`'s VD amendment says they are renameable
at the reshape re-bake. This plan adjudicates for the amendment and emits both
edits; VD's decision is what the rename CONTENT is, and it is schedule-bound to
the reshape's S1.

- **(a) Re-sign at the reshape re-bake to the general staged signatures
  (RECOMMENDED).** `setTreatment` becomes `setForestBasis(sampler, forest,
  basis, numColumns)`; `bcfGlue` becomes `numForestAmplitudes` plus
  `forestAmplitudes`; the creation Doxygen goes to engine vocabulary. Cost ~250
  lines inside a re-bake that is happening anyway. What it buys: the family
  later needs NO header change ever, and the flat surface stops naming a causal
  model at the moment the R surface stops.
- **(b) Ship the BCF names; the mean channel's flat entry becomes a
  post-release APPEND.** Cost: a MINOR bump plus a lockstep bump of stan4bart's
  Depends/LinkingTo floor in the same release - the recorded, non-trivial cost
  that dbarts-h-reshape's resolved question 2 refused for `setForestWeights`.
- **(c) Ship the BCF names and never generalize the flat surface** (the family
  stays R-only). Coherent only if fork 1 is (c), and it makes the flat API
  permanently less capable than the R surface, against binding decision 1(a).

**Sub-question, also schedule-bound to reshape S1: keep or rename the PRECISION
channel?** `setForestWeights` is `dbarts:::`-only today and its flat entry is
built at that slice, after which the spelling is frozen. RECOMMENDATION: KEEP
it. "Weights" already means precisions in dbarts, and `setForestWeights` is the
per-forest sibling of `setWeights`, composing as `w_i m_f^2 s_i`; renaming costs
the parallel that explains it, and naming the mean channel `basis` dissolves the
collision the memo feared without touching the old channel. Alternatives if VD
disagrees: `setForestRowPrecision` (accurate, loses the parallel),
`setForestCaseWeights` (keeps the parallel, adds a word).

**Recommendation: (a) for the rename, KEEP for the precision channel.** What
would change (a): a decision that fork 1 is (c) and the family will never land,
which makes the causal names accurate and (c) the honest choice.

**RESOLVED (VD 2026-08-11): (a) for the rename, KEEP for the precision
channel** - both ADOPTED as recommended. `setTreatment` -> `setForestBasis
(sampler, forest, basis, numColumns)`; `bcfGlue` -> `numForestAmplitudes`
plus `forestAmplitudes`; the creation Doxygen moves to engine vocabulary -
all at the dbarts.h reshape re-bake. `setForestWeights` keeps its name. The
schedule-bound carve-out (dbarts-h-reshape.md binding decision 3) is
satisfied.

**Fork 4. Where does `bcf()` live?**

- **(a) bartCause (RECOMMENDED).** It already owns the causal vocabulary, the
  estimands, common support and the two-counterfactual-surface output contract
  that bcf-public-surface S5 was written against, and its `bcf` slot is one
  commented line from consuming a fitter. It has no `src/`, so it costs nothing
  to compile and can host no C++, which is fine because the sampler stays in
  dbarts. Cost: uncomment the slot, write the response fitter, add a moderator
  EXCLUSION for the propensity-score column, and relocate S5's class (amendment
  5).
- **(b) dbarts.** Coherent, and it is what bcf-public-surface S5 currently
  specifies. Cost: it reintroduces causal vocabulary at the engine's front door,
  which is what VD's removal decision was aimed at, and it puts an estimand
  layer in a package whose stated rule is that it exports none.
- **(c) stan4bart.** Correct only when a parametric block is present (multilevel
  BCF); as the general home it puts a fit function in the WALNUTS collection
  whose boundary rule is "the non-BART component is parametric".
- **(d) A third package.** Declined on the record above: every capability it
  wanted is either already stan4bart's shape or has to be added to dbarts
  anyway, and its one honest form is bartCause.

**Recommendation: (a).** What would change it: a judgment that dbarts should
ship one turnkey causal fit function for discoverability, which is coherent but
inverts binding decision 2.

**RESOLVED (VD 2026-08-11): (a), bartCause** - ADOPTED as recommended, on its
dbarts-1.0 branch. bcf-public-surface S5 relocates there; the moderator
EXCLUSION for the propensity-score-as-covariate column is part of that cost,
not a separate decision.

## What NOT to build, recorded so it is not re-proposed

- **Curated C++ extension headers via `LinkingTo`.** Six facts, all
  re-verified: there is no portable mechanism for one R package to link
  another's static archive, and dbarts links `misc.a`, `external.a`, `rc.a`;
  `misc.a` dispatches SIMD through function-pointer globals only it defines;
  TEN non-bartcore headers would ship, one of them configure-generated, so the
  shipped set would encode one machine's type and ISA decisions; the headers are
  not R-free (`ext_printf` is `Rprintf`, `external/stats.h` pulls `Rmath.h`); a
  user combiner must be injected at `Chain` construction, producing a second,
  disjoint engine that dbarts' own handle, R5 surface and `dbarts.h` cannot
  address, with about 17,000 lines re-parsed at CXX20 per consumer translation
  unit; and the stability commitment is the whole engine's internal shape,
  against the standing ruling that the 1.0-0 contract is the R surface plus
  `dbarts.h`. Under the prototyping principle as qualified by binding decision
  1b this is not any step: `dbarts.h` is the graduation step for PERFORMANCE,
  which is what stan4bart already uses, and the path for AUTHORING a response
  shape, a leaf model, a prior law or a combiner is to add it to dbarts itself -
  the same conclusion stochtree's own documentation reaches for itself ("we add
  our model to the `ModelType` enum"). Extension headers would be a third
  surface serving neither, with the whole engine's internals as its contract.
- **A per-observation multiplier callback in any form.** It would sit inside
  `formForestResponse`'s `i` loop and its inner `g` loop, violating the
  combiner design's own rule that combiner work is per-SWEEP, never inside the
  per-observation kernels.
- **A new state-dependent-multiplier hook (the memo's held door).** Deleted per
  B1: the sweep-boundary hook it describes already ships.
- **An R combiner callback fired inside the forest loop.** Refused on
  invariant grounds (`totalFits` is stale during the tree loop) as well as
  threading grounds. Never write "R callbacks are impossible": a single-chain
  sweep-boundary R hook ships.
- **Exposing a per-sweep R callback as an on-ramp affordance.** Reconsidered
  under the prototyping principle and DECLINED with a receipt: the loop it
  would replace is measured FREE (`run(0,1)` per sweep is 0.89-0.95x of a
  batched run and bitwise identical), so the hook buys nothing an R prototyper
  needs. The internal single-chain hook stays what it is - rbart_vi's mechanism
  for interleaving its own parametric draw with caller-owned per-sweep buffers -
  and M0 records it rather than promoting it. Its honest scope is ADDITIVE and
  glue-adjacent intervention, per binding decision 1b, and rbart_vi is the
  boundary case that shows the line: its random intercepts reach the forest as an
  OFFSET and the quantity it draws in R is its own block's variance, so nothing
  about the response shape or the forest's leaf prior is authored in R. It is
  therefore NOT a route to prototyping a new family, link or prior law, and must
  not be priced as one. Widening it to inline multi-chain is a DOOR, priced: its
  results buffers are single-chain by construction and its lambda drops
  `chainIndex`, so it is roughly 60 bridge lines plus an R wrapper and tests,
  not a predicate change.
- **A multiplier folded into `setForestWeights`.**
- **A third package.**
- **Multi-arm treatment as K forests.** It is a vector leaf sharing one
  partition. The arc that buys it is "lift the two constant-leaf asserts so a
  combiner may carry a vector leaf", which per B2 also requires rewriting BOTH
  `afterCombine` bodies (or documenting a refusal of the ASIS move under a
  vector leaf) and deciding saved-tree replay, and which does NOT buy
  multivariate-outcome BCF. A door with a corrected price, not a cheaper
  alternative to fork 1.

## Doors held open (recorded, not scheduled)

- **Multi-arm treatment on a vector-leaf combiner**, priced as above.
- **Multivariate-outcome BCF**, which needs a d-dimensional response, fits,
  residual covariance and Results layout: a response-model arc, not a leaf arc.
- **A per-forest amplitude with a caller-supplied covariance kernel drawn in
  the combiner** (LongBet-shaped panel), the re-scoped remains of B1's door: a
  spec question for M4, not an ABI question.
- **Per-forest saved-tree replay**, which is what a component `predict` on new
  rows and bairrtt's filter need, and what would retire
  `testFitsAreDefined() == false`.
- **Pinned (non-sampled) amplitudes**, VCBART's own parameterization, if FA1
  says the amplitude is decoration for non-BCF instances.
- **A test basis vector**, which would make a K-forest sampler's blended test
  fits defined.
- **Grouped x K-forest** (stan4bart's multilevel BCF): the combiner composes
  with `GroupedResponse` by design, the constructor does not build the
  decorator.
- **Inline multi-chain for the internal R sweep hook**, priced above.
- **The flat surface's graduation debt**: per-observation predictor update and
  the joint session, `setCutPoints`, `setData`, `predictVariance`,
  forest-indexed `predict` - all recorded in `dbarts-h-reshape.md`, listed in
  M0 as on-ramp debt.

## Departures from the memo and the critique (record)

1. **The memo's declarative spec is OVERTURNED and replaced** (fact 9). Its
   "amplitude times an n x K multiplier matrix" does not contain BCF, because
   `combinedFits` uses `b_{z_i}` and `drawGlue` draws `b0` and `b1` as two
   conjugate coefficients over the two z subsets (`combiner.hpp:567-576`,
   `:609-622`), i.e. forest 1 carries a two-column indicator basis. The channel
   is a per-forest basis MATRIX with a coefficient VECTOR. Its falsifier F2 is
   struck and replaced by FA0.
2. **The memo's cost estimate for the family is raised**, from 350-500 engine
   lines to 450-650, on departure 1 plus the q-variate amplitude conditional -
   which is cheaper than it sounds because `choleskyDecompose` and the
   triangular solves already ship (`model.hpp:869-895`).
3. **The memo's 1.04x is void wherever it is attributed to the R composition**
   (B3, adopted): 1.39-1.43x measured at K = 2, and 0.89-0.95x for the engine's
   own per-sweep drive.
4. **The memo's "That is not an approximation" is corrected**: the BLOCKING is
   identical, the leaf and amplitude priors are not, and the difference has a
   named mechanism (the pinned transform's `shift` attributed once per forest
   rather than once) and named remedies, which M0 documents and FD4 tests.
5. **The memo's held route-3 door is DELETED** (B1) and its LongBet sampling
   claim is not repeated against the paper.
6. **The memo's open decision 5 is re-priced and half struck** (B2).
7. **The memo's amendments to the two committed plans are replaced**: its
   amendment 1 and its Q7 are moot (S3 landed); its 7.3 amendment 2 cited the
   wrong binding decision; the real item is the plan conflict, adjudicated here
   with both edits emitted.
8. **The memo's D2 is half moot** (fact 11): only the SETTER is internal.
9. **The memo's "no document states dbarts as a component as a supported
   contract" is half wrong** (fact 13): the vignette and an Rd example exist;
   what is missing is the contract and the K-forest recipe. M0's scope shrank
   accordingly, from a new document to a new design doc plus a vignette
   section.
10. **The memo's 2.1 falsifier is re-phrased** per the critique's A1: it must
    read "no line reads a per-observation multiplier VALUE from caller data",
    since `forestMultiplier` does read `glue_.z[i]`, caller-supplied borrowed
    data, to select `b1` over `b0`.
11. **The memo's "borrowed pointer" wording is corrected** per the critique's
    A4: flat `setTreatment` COPIES (`dbarts.h:43`) and the bridge copies into
    holder-owned storage. The recommendation (a raw pointer plus a count, not a
    struct) survives; the ownership question is answered explicitly (COPY) in
    the channel spec rather than left implicit.
12. **The memo's non-Gaussian pricing is corrected DOWNWARD** (fact 12): the
    combined-fit latent draw and a non-Gaussian K-forest chain both already
    exist; the BCF-path constructor's hardcoded `GaussianResponse` and the
    Gaussian calibration map are the actual work.
13. **The critique's B5 consequence is NARROWED, not adopted.** Its finding
    (the precision channel's spelling freezes at reshape S1) is verified and
    adopted; its conclusion (the new channel's name must be answered before that
    slice) does not follow once the mean channel is spelled `basis`, because
    nothing about a mean channel wants the word "weights", which in dbarts means
    precisions. What survives as schedule-bound is only the question of
    renaming the OLD channel, put to VD in fork 3 with a recommendation to keep
    it.
14. **The critique's implied ordering of the driver surface as a thin
    docs slice is superseded by VD's prototyping principle**, which arrived
    after both documents: M0 is the arc's headline deliverable, its measured
    composition cost is a prototyping price rather than a defect, and it gains a
    semantic-continuity section and FD4, neither of which either document
    contemplated. The recommendation's ORDERING did not change (driver surface
    first, family second) because both documents already reached it; what
    changed is the family's JUSTIFICATION, which is now the non-Gaussian case
    keeping multiforest prototyping in R rather than Gaussian ergonomics, and
    the scope of M0, which grew a graduation contract.
15. **VD's scope limit on the R-first principle (binding decision 1b) is adopted
    with NO disagreement and ONE refinement.** The cut it draws - additive
    composition on R's side, response shape and priors on the engine's - is the
    cut the evidence already supported: the measured pure-R path reproduces the
    additive BLOCKING exactly and misses only the pinned response transform and
    the leaf calibration (fact 5), which are precisely "the shape of the response
    or the priors". The refinement is a middle category the two-way phrasing
    leaves implicit: CONFIGURING an engine-provided family or prior by argument
    (`family =`, `node.prior = normal(k)`, the per-forest knobs inside
    `forest()`) is on R's side, while AUTHORING a new family, link, latent
    structure, leaf model or hyperprior law is on the engine's. Without that
    middle category the boundary would read as excluding prior arguments from the
    on-ramp, which no evidence supports and which would make `forest(base =,
    power =, sd =)` look boundary-crossing. Three concrete effects on this plan,
    all recorded above: the per-sweep R hook is explicitly not priced as a route
    to prototyping families or priors (rbart_vi is the boundary case and sits on
    the additive side); the basis channel is specified to carry data plus prior
    SELECTION and to refuse a prior form the engine does not implement; and the
    family's non-Gaussian justification is restated as supplying columns against
    an engine-provided family, which is composition, not response-shape
    authorship. Nothing in the recommendation, the slice list, the falsifiers or
    the four forks changed as a result; the scope limit tightened the wording and
    added the boundary statement to M0's graduation section.
16. **Small corrections not worth a finding**: the non-bartcore header closure
    is 10 and excludes `misc/simd.h`; the memo's per-header line counts are one
    high on five headers (they came from out-of-range reads), and bartcore is
    about 17,000 lines across 11 headers; four of the memo's UNVERIFIED venues
    are published and one supposedly-unverified venue does not exist at all
    (Woody et al. is unpublished, which weakens fork 1 option (b), not the
    census).

## Provenance and honest limits

    repo      /Users/vdorie/Repositories/dbarts, worktree .claude/worktrees/
              zero-weight, branch wt/zero-weight, tip db81bfe. src/ read
              through `git show HEAD:` because a code implementer is editing
              it concurrently (bcf-public-surface S4). Nothing was built here.
    verified  every B1-B6 anchor; the reshape's binding decision 3 and its S1
              item 5; bcf-public-surface's AMENDED block; dbarts.h's four BCF
              entries, hash, creation Doxygen and ownership sentence; the
              chain sweep loop and its onSweep placement; both constant-leaf
              asserts and the muByTree/afterCombine coupling; VectorLeafModel
              and the ResponseFamily enum; setForestWeights at all three
              layers and its absence from NAMESPACE; the S2 R5 readers; the
              BCF and multinomial chain constructors; the Cholesky helpers;
              the header include closure; the vignette and the Rd example.
    taken as  the critique's four numerical probes (the 1.8e-15 combined-fit
    fact      probe, the rebuild divergence, the timing trio, the F6 run) and
              its citation audit. Not re-run here.
    not read  docs/plans/capi-callbacks.md past line 70; docs/plans/
              facade-shape.md. Both bear on the capability-probe pattern any
              multiforest API should use instead of a forest count.

**Four things this plan cannot tell you.** (1) How bad a forest-versus-forest
ridge is: nobody has measured one. Do not cite `composition-mixing-probe.md`'s
IACT 850-880 panel as a forest-versus-forest result - every arm there is a
parametric block against a forest, the panel is reported rather than gated, and
its coordinate is not comparable across arms; and its own open item
(`:1087-1091`) declines to bound even forest-versus-parametric against the
group-intercept surrogate in either direction.
(2) Whether the amplitude and its ASIS remedy earn their place for any class but
BCF: FA1 is the pre-registration and has not been run. (3) What the family's
K > 2 instances cost under multiforest-predictor-mutation's veto: no
measurement exists above K = 4, and the per-forest column mask is vacuous for
VCBART-shaped models. (4) Whether the q-variate amplitude conditional can be
made bitwise identical to BCF's two scalar draws; if it cannot, M4.2 either
keeps BCF's specialized path or owns a bcf-equivalence re-record.
