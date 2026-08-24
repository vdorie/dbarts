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
window: M0, M1, M2 land INSIDE the pre-release breaking window. M3 LANDED
  riding dbarts-h-reshape S1 (ab3aa2fa), the SECOND AND LAST re-bake of that
  window - the mean channel's flat spelling was decided there.
  M4 is a PRE-RELEASE arc, scheduled after M3/the reshape (RESOLVED
  2026-08-11), and needs no header movement, M3 having landed. Everything is
  breakable and the four sister packages migrate in lockstep at the freeze,
  once, against the final header (VD 2026-08-10).
  "No header movement" is VERIFIED, not forecast (amended 2026-08-13, pre-M4):
  creation crosses as SEXP (`dbarts_sampler_create(SEXP, SEXP, SEXP, const
  char* family)`, `dbarts.h:506-507`), `setForestBasis(..., forest, basis,
  numColumns)` (`:700-701`) is already general, `numForestAmplitudes(...,
  forest)` (`:715-716`) plus `forestAmplitudes(..., forest, out)` (`:724-725`)
  are already ragged, and the family selector is a `const char*`. The limiting
  Doxygen sits OUTSIDE `DBARTS_C_API_LIST` (`dbarts.h:297-420`), so relaxing it
  moves no hash - reshape S1 proved that empirically for three such edits
  (`dbarts-h-reshape.md:1747-1749`). SIX paragraphs, not two, change meaning
  under K > 2 and must be swept by M4.5: `:505` ("Gaussian responses only."),
  `:696-699` ("Today's engine honours exactly one basis..."), `:133-134`
  (`logLikelihood`'s "a BCF two-forest fit"), `:618-620` (the saved-tree
  `forest` argument), `:681-683` (`numForests`, "2 for a Bayesian causal
  forest"), `:779-780` (`setForestPriorScale`'s "a two-forest or multinomial
  sampler"). None is a signature; none is inside the hashed macro.
  The six SPLIT, re-scoped 2026-08-13 on the probe verdicts
  (`probes-2026-08-13.md`): five ride M4.5, but `:505` ("Gaussian responses
  only.") stays ACCURATE until M4.4 and moves with M4.4, because M4.5 now lands
  Gaussian-complete AHEAD of it (fork 2's escape hatch, taken). Two
  near-misses are explicit NON-GOALS: a per-draw amplitude channel in flat C
  (declined, see the M4 section) and any re-sign of
  `dbarts_sampler_setForestBasis`'s error channel (a guard-body change only -
  its two-channel Doxygen at `:693-699` already anticipates the widening
  verbatim).
budget: M0 ~220 design doc + ~180 vignette + ~200 test. M1 ~75-85 R + ~50-60
  man + ~140 test, ~270-290 total (RESTATED 2026-08-13 pre-M1, superseding the
  ~60 + ~40 + ~140 committed 2026-08-11; the carried items are in M1 item 5).
  M2 ~200-240 R + ~120-160 man + ~435-475 test, ~775-895 total (RESTATED
  2026-08-13 pre-M2, superseding the ~180 + ~60 + ~260 committed 2026-08-11;
  the carried items are in M2 item 12).
  M3 ~40 header + ~80 C_interface + ~70 consumer.c + ~60 test-capi.R, inside
  dbarts-h-reshape S1's own budget. Total in-repo before the freeze ~2145-2265
  (~1630 before M2's 2026-08-13 restatement, ~1590 before M1's).
  M4 ~1850-2400, RE-PRICED 2026-08-13 pre-M4 and superseding both the ~600-900
  written here and the ~900-1100 the M4 section and fork 1 option (a) carried;
  the split and what each addition is are in the M4 section's budget paragraph.
  The total above EXCLUDES M4, which is priced in its own section.
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
amended: 2026-08-13, PRE-M2, by the orchestrator, against 4c7aa200 (M1 landed
  05ac3b4b plus its records commit). The M2 section's anchors are re-verified
  BY SYMBOL and given with their stale values - M1 added ~67 lines to
  `R/dbarts.R`, moving every R5 anchor the section carried - and the section
  gains seven items: the surface bookkeeping M2 owes `NAMESPACE`,
  `_pkgdown.yml` and `man/` in both directions (item 6); the settled 1-based
  rule applied, with `$getBCFGlue`'s argument-FREE signature recorded as the
  fact that keeps `$getForestAmplitudes` argument-free too (item 7); the three
  M1 interactions, all clean (item 8); the verified migration list, with the
  finding that `benchmarks/R/bcf-equivalence.R` needs NO edit because it drives
  the internal route (item 9); the refusal count corrected from eleven to
  seventeen inherited plus one retired, and the knob map found INCOMPLETE (item
  10); `dbartsData(treatment = )` FLAGGED as an unsettled fourth site (item
  11); and a restated budget (item 12). Nothing about the arc's content,
  ordering or decisions moved.
  FOLLOW-UP the same day, still pre-M2: two of those three flags are now
  DECIDED by the orchestrator under the standing grant, VD may veto. (1) The
  knob map is settled and total - the three unhoused `treatmentForest()` knobs
  become `forest(amplitude.prior.variance =)` and `forest(update.amplitude =)`,
  with `amplitude.prior.variance` legal only on a forest given `basis =` and
  `update.amplitude` legal on every forest of a `forests =` spec (the literal
  basis-gate had to be split there or `update.a` would be unreachable and M0's
  FD4 unexpressible; item 10 records the split and the reasoning). Item 1
  carries the complete constructor, item 4 goes from six refusals to nine, and
  item 12 gains ~15 test lines for the FD4-expressibility pin - the three new
  refusals themselves ride the refusal matrix item 4 and FS2 already own.
  (2) `dbartsData(treatment = )` and the `data@treatment` slot STAY AS-IS at
  M2 as recorded data-side debt, with `forests =` populating the slot
  internally exactly as `treatment =` did; item 11 is now a settled
  disposition rather than two options, the escalation sentence is struck (the
  split surface is shippable), and the M2 RECORDS commit opens the TODO ticket
  for the debt. The third flag, the refusal count, was already corrected in
  place at 3fb77060 and stands: 17 inherited, 1 retired.
amended: 2026-08-13, PRE-M4, by the orchestrator, against 25a21d3b (M1
  05ac3b4b, M2 64b13b98, reshape S0 a262cd26, item 9 a14040de, S1 with M3
  folded inside it ab3aa2fa, S2 1bf2e69c, arc close 25a21d3b - all landed
  since the M4 section was authored 2026-08-11 against a tree where
  `treatment =` was the creation route, `setTreatment`/`bcfGlue` were the flat
  C names, and the API hash was `0x1a911c00bb26dcd7`). Artifacts:
  `.claude/m4-basis-design/pre-m4-scoping-2026-08-13.md` (scoping),
  `critique-2026-08-13.md` (adversarial critique), `synthesis-2026-08-13.md`
  (the adjudications, one line each) - all GIT-IGNORED, so the load-bearing
  facts are carried here rather than pointed at. 31 anchors/claims re-verified
  BY SYMBOL: 22 verified, 8 line-number drifts, 1 name GONE. No M4 claim is
  falsified by the tree and no slice loses its target. Where the scoping pass
  and the critique conflict the critique wins, EXCEPT on three anchors where
  the critique's own correction is refuted by the tree and the verified value
  is used instead: `R/spec.R`'s composition refusals are the vector `:429-456`
  emitted at `:457-463` (not `:428-463`/`:464-469`), and
  `test-bcf-reporting.R`'s three-chain glue pin is `:114` (not `:113`). What
  moved: the window paragraph and the M4 budget above; the whole M4 section
  (preamble, sequence, all six slices, all seven falsifiers, the budget); the
  two stale `chain.hpp:1261-1262` anchors, both now `:1266-1267`. Nothing about
  the arc's decisions or its ordering of M4.0-M4.5 moved.
  DRIFT SWEEP, `dbarts.h`, applied once here rather than by editing landed
  records: reshape S1 moved every citation this file carries. `:43` (the
  ownership sentence) is now `:45-49`, and it names `setForestBasis` explicitly
  ("setForestBasis copies the basis columns it is handed", `:48-49`);
  `:348-357` (the creation contract) is now `:484-505`, re-worded end to end in
  `forests =`/`forest(basis =)` vocabulary; `:372-378` (the callback refusal
  predicate) is now `:520-529`; `:264-271` (the four BCF entries) is now
  `:389-402`, and TWO of the four names are GONE - the live set is
  `numForests`, `setForestBasis`, `forestFits`, `numForestAmplitudes`,
  `forestAmplitudes`. The claim cited via `dbarts-h-reshape.md` that the hash
  is blind to struct layout re-verifies at `dbarts.h:95-99` plus
  `C_interface.cpp:280,334-344`: it is FNV-1a over `DBARTS_C_API_DECLS`,
  signatures only. The M3 section and the cross-plan amendment block are
  LANDED RECORDS and keep their pre-re-bake numbers, marked as such in place.
re-scoped: 2026-08-13, POST-PROBES, by the orchestrator, on the verdicts in
  `.claude/m4-basis-design/probes-2026-08-13.md`. That file is GIT-IGNORED, so
  its numbers are carried into this one verbatim, each with the protocol that
  produced it. **Provenance of the run:** tree sha 47c1fbe1, equal to
  `origin/bartcore`, worktree clean, no tracked file edited; the package is
  the pre-existing verified install `.claude/privlib-s2-docs`, confirmed by
  `DBARTS_C_API_HASH == 0xcd88efcd67de55d7`; base seed 20260813 for FA1a and
  FA1b (per replicate s in 1..8 the DGP seed and `rngSeed` are both 20260813 +
  s, so the arms are paired), dataset seed 20260813 for FA5 with chain seeds
  5101-5108 (arms), 5901-5908 (leg-P reference) and 5801-5808 (leg-G reference
  validation); harnesses, logs and `.rds` results gitignored under
  `.claude/m4-basis-design/harness/`. **Three things moved and nothing else.**
  (1) FA1's pre-registered re-scope is NOT licensed - both conjuncts fail - so
  M4.2 keeps the q-variate amplitude conditional and the per-forest ASIS
  rescale in FULL scope, and FA1 is not promoted to `benchmarks/R/`.
  (2) FA5 FALSIFIED M4.4's headline ground (arm B AGREES on all 12 functionals
  at max |z| = 2.54), so fork 2's Gaussian-first escape hatch is TAKEN: the arc
  order is M4.0 -> M4.1 -> M4.2 -> M4.3 -> M4.5, with M4.4 the IMMEDIATE
  follow-on slice, still pre-release, re-justified on the surviving grounds
  only. (3) The 1.39-1.43x composition tax this file asserted in eight places
  is STRUCK as unsourced and replaced by measured, protocol-stated figures -
  the erratum is "Departures" item 3. Nothing about the arc's decisions moved
  beyond the M4.4/M4.5 order.

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
(re-scoped 2026-08-13 on the probe verdicts: **5.14x at K = 2 and 5.43x at
K = 8**, per-sweep wall time of the K-sampler R composition over a SINGLE
batched engine sampler at matched total trees, same n, back to back -
`probes-2026-08-13.md` FA1b; the 1.39-1.43x this sentence carried is struck,
"Departures" item 3 - and a slightly different prior). Its irreducible engine
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
1.39-1.43x at K = 2, growing with K (K stores, K `.Call` round trips, O(nK)
R-level residual arithmetic per sweep). Every occurrence of the memo's 1.04x
attributed to an R composition is void.
**ERRATUM, 2026-08-13 (post-probes).** The composition half of this adoption is
STRUCK: **1.39-1.43x has no receipt anywhere in the tree** - the grep finds the
assertion, never a measurement - and the tool-verified-claims discipline applies
to this plan's own numbers, not only to the memo's. The engine half keeps its
receipt (`bcf-public-surface.md:85-93`) and stands. The replacement is FA1b's
measured figure with its protocol stated: **5.14x at K = 2, 5.36x at K = 4,
5.43x at K = 8**, per-sweep wall time of the K-sampler R composition over a
single BATCHED engine sampler carrying the same total tree budget (K * 50) on
the same n, measured back to back on the same box - and the tax is FLAT in K
(1.06x growth from K = 2 to K = 8), so "growing with K" is struck too. Full
erratum in "Departures" item 3.

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
6. **The honest cost numbers** are B3's engine figure, 0.89-0.95x per-sweep
   drive, plus - re-scoped 2026-08-13 on the probe verdicts, B3's composition
   figure being unsourced - FA1b's MEASURED R K-sampler composition tax:
   **5.14x at K = 2, 5.43x at K = 8**, against a single batched engine sampler
   at matched total trees, flat in K (1.06x growth).
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
    (`chain.hpp:1266-1267`) and the multinomial K-forest chain already runs a
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
ROW-major - corrected 2026-08-14 at M4.5; this line said column-major, and the
shipped contract is row-major at every layer, `combiner.hpp:380-384` and
`dbarts.h:707-712`, M4.3 item 5 having resolved the transpose in the doc
direction) and a coefficient vector `a_f` of length q_f. The mean is

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
the R composition at 5.14x (K = 2) to 5.43x (K = 8) against a batched engine at
matched total trees and a different prior (fact 5 and fact 6 as re-scoped
2026-08-13, not the memo's 1.04x), and needs three things that live only in
`chain.hpp`/`combiner.hpp`: the latent draw against the combined fit, the ASIS
rescale that writes leaf values, and the near-zero snap. Its Gaussian half
collapses into being a stan4bart model class (the `mvbart()` precedent adds a
model class there with zero dbarts change); its causal reporting layer is
bartCause; its only engine content argues
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
                                            VD 2026-08-11); probes RAN
                                            2026-08-13; M4.0 LANDED 562ee684,
                                            M4.1 LANDED 1458328c, M4.2 LANDED
                                            1a2aaedc, M4.3 LANDED 9c63e9d8,
                                            next step M4.5, then M4.4 (hatch
                                            order)
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
   `refuseUndefinedTestFits`, `refuseMultiForestMutation`,
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
   residual), its measured price (re-scoped 2026-08-13 on the probe verdicts:
   **5.14x at K = 2, 5.43x at K = 8**, against a single batched engine sampler
   at matched total trees - flat in K, NOT growing with it; the vignette states
   that denominator, since a price quoted without one is what the struck
   1.39-1.43x was), and the
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
   `forest(basis =, vars =, n.trees =, base =, power =, sd =, interactions =,
   blocks =, amplitude.prior.variance =, update.amplitude =)` as the single
   constructor (the `interactions()`/`blocks()`/`treatmentForest()` precedent:
   fourteen knobs ride one constructor, so `dbarts()` grows exactly ONE
   argument, down from three). The last two are the amplitude-prior knobs
   DECIDED 2026-08-13 pre-M2 (orchestrator discretion under the standing
   grant, VD may veto); the map is complete and total, and item 10 walks it
   knob by knob with the two refusals it adds.
2. `resolveSamplerSpec` resolves `forests =` to the SAME
   `attr(control, "bartcore.bcf")` payload plus `data@treatment` that
   `treatment =` resolves to today, expanding a factor basis to its level
   indicators in R so the bridge sees exactly what it sees now.
3. `treatment =`, `moderators =`, `treatmentForest =` are REMOVED, and
   `$setTreatment`/`$getBCFGlue` are renamed to
   `$setForestBasis`/`$getForestAmplitudes` with the old spellings gone.
   `inst/tinytest/test-bcf-creation.R` is rewritten; the "other ten files
   touching `treatment`" this item named before the 2026-08-13 pre-M2 refresh
   was never enumerated and is wrong in both directions - item 9 replaces it
   with the verified list.
4. Every declaration today's engine cannot honour REFUSES at creation, one
   assertion each. NINE, up from the six this item carried before the
   2026-08-13 pre-M2 refresh: K > 2 forests; a numeric (non-factor) basis; a
   factor with more than two levels; a basis on forest 0; a non-Gaussian
   family; a variance forest alongside; plus the three the settled knob map
   adds (item 10) - `amplitude.prior.variance` on a forest with no `basis =`,
   any amplitude knob on a K = 1 `forests =` spec, and the top-level
   `interactions =`/`blocks =` supplied together with the same knob on forest
   0. On top of those, the shipped refusal set is inherited unchanged (item 10
   corrects "the eleven S1 refusals" this item carried, and names the one
   refusal M2 RETIRES rather than inherits).
5. Docs: `docs/design/bcf.md` and `docs/design/public-surface.md` gain the
   engine-vocabulary surface; `inst/NEWS.Rd`; `_pkgdown.yml` if a new topic
   appears (the recorded new-exported-Rd-topic lesson). Item 6 settles that a
   new topic DOES appear, and that an old one goes.

**Anchors, re-verified BY SYMBOL at 4c7aa200 (added 2026-08-13, pre-M2).**
Every number this section carried predates M1 (05ac3b4b, ~67 lines added to
`R/dbarts.R`) and the mask arc; these are the live ones, and where the section
cited a number it is given with its stale value.
Creation formals: `dbarts()`'s three at `R/dbarts.R:350-352`, reaching
`dbartsData()` through `dataCall <- redirectCall(matchedCall, ...)` (`:539`)
and `resolveSamplerSpec` at `:562-566`, whose comment (`:562-563`) records why
`treatment = NULL` is passed there - the vector already rides the data object,
because `dbartsData()` is the one place that knows which rows `subset` kept.
That redirect is load-bearing for M2: dropping `dbarts()`'s formal alone does
not remove the argument if `dbartsData()` keeps one of the same name (item 11).
`dbartsSpec()` `R/spec.R:471`, its formals `:486-488`, forwarding `:567-569`;
`resolveSamplerSpec` `R/spec.R:15`, formals `:32-34`. The resolution path is
`R/spec.R:363-460`: the design comment `:363-371`, `data@treatment <-
validateTreatment(...)` `:373`, the gaussian gate `:375-384`, the `unsupported`
vector `:391-420`, its `stop` `:421-426`, `resolveTreatmentForest` `:427`,
`resolveModerators` `:428`, the `attr(control, "bartcore.bcf")` payload
`:429-454`, the orphan-argument gate `:455-460`. Constructor and resolvers:
`treatmentForest()` `R/model.R:1238-1264` with its doc `:1228-1237`,
`resolveTreatmentForest` `:651-681`, `resolveModerators` `:620-646`. Data half:
`dbartsData()`'s own `treatment` formal `R/data.R:697`, `validateTreatment`
`:631-654` called at `:880`, `:962`, `:1028`; the slot `R/A_class.R:518`
(comment `:513`, prototype `:537`, validity `:582-591`). R5: `$setTreatment`
`R/dbarts.R:1185-1212` (was cited `:1132`), its docstring `:1186`, the
`data@treatment` mirror `:1196-1197`, the `.Call` `:1199`; `$getBCFGlue`
`:1402-1406` (was cited `:1349`), between `$getForestFits` `:1397` (was
`:1344`) and `$getForestVariableCounts` `:1407` (was `:1354`);
`$getCalibration` `:1412` (was `:1362`) and `$setCalibration` `:1417` (was
`:1393`), the two `resolveForestIndex` precedents, at `:1415` and `:1446`.
Rd: aliases `man/dbartsSampler-class.Rd:21` and `:34` - the block gained
`setForestWeights` at `:20` in M1, so both stale by one - usage `:66` and
`:81`, the mutator list naming `setTreatment` `:104`, the `updateScale`
refusal sentence `:131`, `\item{z}` `:166-168`, `\item{forest}` `:169-171`,
`\value` `:355`; `man/dbarts.Rd` usage `:20` and the three `\item` blocks
`:83-85`, `:86-88`, `:89-91`, plus the calibration sentence `:153`;
`man/dbartsSpec.Rd:17` and `:28`; `man/dbartsData.Rd:13` and `:16`;
`man/treatmentForest.Rd`, 92 lines; `inst/NEWS.Rd:85-95`, the R5 bullet naming
both renamed methods, and `:693`, the FLAT entry, which is M3's and not M2's.
One anchor outside R and decisive for item 11: the bridge reads the
`"treatment"` slot BY NAME at `src/R_interface_bartcore.cpp:1112` (length check
`:1117`), so `data@treatment` cannot move inside an R-only slice.

**6. Surface bookkeeping, made explicit (added 2026-08-13, pre-M2).** The
section hedged this ("`_pkgdown.yml` if a new topic appears"); scoping settles
it in both directions, and each half has a gate consequence.
   (a) `treatmentForest` is EXPORTED (`NAMESPACE:18`, verified) and is a listed
   pkgdown reference topic (`_pkgdown.yml:33`, verified, inside "The mutable
   sampler"). M2 removes it from BOTH, and deletes or rewrites
   `man/treatmentForest.Rd`. An export removed without its `_pkgdown.yml` line
   leaves `check_pkgdown` pointing at a topic that no longer exists.
   (b) `forest()` is a NEW exported constructor with a NEW Rd topic, so M2 adds
   the `NAMESPACE` export, writes `man/forest.Rd`, and adds the
   `_pkgdown.yml` entry beside `interactions`/`blocks` where `treatmentForest`
   sat. This is the new-exported-Rd-topic lesson's fifth occurrence; its fourth
   was caught at review rather than by a gate (`bcf-public-surface.md:886-888`),
   which is why `pkgdown::check_pkgdown` is a NAMED gate below rather than
   "pkgdown check", and why `R CMD check` (undocumented-export and codoc) joins
   M2's gate line - M1 needed neither, having no new topic.
   (c) No name collision: `forest` is unused as a symbol anywhere in `R/`
   (verified), so the export is free. It does collide with the loop variable in
   `reapplyForestWeights` (`R/dbarts.R:1483`) only lexically, which is harmless
   and worth neither a rename nor a comment.

**7. Forest indexing: the SETTLED 1-based-at-the-boundary rule applies to both
renamed methods (added 2026-08-13, pre-M2).** Cited, not restated: M1 item 1
(VD 2026-08-13). `$setForestBasis(forest, basis)` takes a forest index and
routes it through `resolveForestIndex` (`R/bartcore.R:1051`), as
`$setForestWeights` and `$getCalibration`/`$setCalibration` do; the treatment
forest is `2L` in every M2 example, Rd sentence and test.
`$getForestAmplitudes` is the exception, and it is a FACT not a choice:
`$getBCFGlue` at HEAD takes NO argument (`getBCFGlue = function()`,
`R/dbarts.R:1402`) and returns the whole 3 x n.chains `(a, b0, b1)` matrix, so
the plain rename is argument-free too. Giving it a `forest` argument is a
SURFACE CHANGE, not a rename: it would have to return forest 0's one amplitude
or forest 1's two, i.e. the ragged shape M3's `numForestAmplitudes` plus
`forestAmplitudes` pair exists to size. M2 keeps the argument-free
`$getForestAmplitudes()` and the 3-row matrix; the forest-indexed reader
arrives with M3's flat pair or with M4's K-length amplitude vector, whichever
lands first. Recorded so a later slice does not read the rename as having
settled the reader's shape. `$getForestFits`/`$getForestVariableCounts` stay
0-based, untouched, still the `r5-forest-indexing` ticket's two named cases.

**8. Interactions with M1, all three checked (added 2026-08-13, pre-M2).**
   (a) **The `data@treatment` mirror SURVIVES the rename, unchanged in
   mechanism.** `$setTreatment` writes the engine and the slot under a
   `tryCatch` that rolls the slot back on a bridge error
   (`R/dbarts.R:1196-1205`), which is what makes `getPointer()`'s transparent
   re-creation carry the current assignment. `$setForestBasis` keeps exactly
   that body: the slot is read back by name from C at every create
   (`src/R_interface_bartcore.cpp:1112`), so no R-only slice can replace the
   mechanism, and M0 item 1(iii) commits it as contract. What moves is the
   ARGUMENT's spelling - `z` becomes the two-column indicator `basis`, and the
   method reduces it to the 0/1 column the slot and the bridge still take, on
   the same expansion rule item 2 uses at creation. The Rd's `\item{z}`
   (`man/dbartsSampler-class.Rd:166-168`) becomes `\item{basis}` and keeps its
   mirroring sentence verbatim.
   (b) **M1's `forestWeights` field and its three re-application sites are
   name-independent of both renames. Verified, no collision.** The field
   (`R/dbarts.R:725-728`, initialized `:745`), its install
   (`:1167-1172`), `reapplyForestWeights` (`:1481-1489`) and its three call
   sites - `$copy` `:906-907`, `getPointer()` `:1518`, `setState()` `:1544` -
   name neither `setTreatment` nor `getBCFGlue` and touch neither
   `data@treatment` nor the BCF control attribute. Two Rd adjacencies, both
   benign: M1's `setForestWeights` paragraph (`:149`) says "built with
   `treatment = `" and the method's own docstring (`:1150`) says the same, so
   both need the one-phrase re-word to `forests = ` that every other
   `treatment = ` reference needs; and `\item{forest}` (`:169-171`) already
   carries three clauses (the 0-based readers, the 1-based
   calibration pair, M1's `setForestWeights`) and gains a fourth for
   `setForestBasis` on the 1-based side. Paragraph ORDER does not move.
   (c) **`inst/tinytest/test-forest-weights-r5.R` builds its BCF samplers with
   `treatment = ` (`:34`) and reads `$getBCFGlue` (`:166`), so M1's own test
   file migrates inside M2.** It is one of the seven in item 9; the budget
   covers it there.

**9. The migration, enumerated (added 2026-08-13, pre-M2).** Two routes create
a BCF sampler and only ONE of them moves. The PUBLIC formals go; the internal
`dbarts:::bartcoreBCFSampler` route (`R/bartcore.R:622-680`, with its own
`n.trees.treatment =` and `moderators =` formals) is untouched by M2 - it is
`dbarts:::`-only, absent from `NAMESPACE`, and its flat siblings are M3's.
   PUBLIC-formal creation calls that M2 must migrate, all under
   `inst/tinytest/`, verified by grep at 4c7aa200: `test-bcf-creation.R`
   (about 30 hits, REWRITTEN wholesale, 385 lines today, 26 `expect_error`);
   `test-bcf-r5-surface.R` (4 creation calls `:33`, `:73`, `:131`, `:179`, plus
   2 `$setTreatment` `:111`, `:134` and 2 `$getBCFGlue` `:38`, `:147`, plus
   4 header/section comments); `test-bcf-reporting.R` (3 calls `:35`, `:82`,
   `:99`, plus `$getBCFGlue` `:91`); `test-calibration-creation.R` (2, `:198`,
   `:208`); `test-active-rows-pins.R` (2, `:94`, `:103`);
   `test-calibration-midchain.R` (1, `:333`); `test-capi.R` (1 call `:522-523`,
   the only migrating site that also passes `treatmentForest()`);
   `test-forest-weights-r5.R` (1 call `:34`, plus `$getBCFGlue` `:166`).
   SEVEN files plus the rewritten one; 14 creation calls and 5 R5-method call
   sites outside the rewrite.
   NOT migrating, because they drive the internal route: `test-bcf.R` (19
   hits), `test-forest-weights.R` (7), `test-interactions.R` (4),
   `test-blocks.R` (4), `test-multi-forest-seam.R` (2),
   `test-bcf-zero-multiplier.R` (2), `test-bcf-mutation-pins.R` (2),
   `test-multinomial-test-offset.R` (1), `test-multinomial-counts-mutation.R`
   (1), `test-multinomial-category-offset.R` (1) - the three multinomial ones
   build a BCF sampler only as a refusal fixture. `vignettes/` and
   `inst/common/` have ZERO hits (verified); no benchmark under
   `benchmarks/R/` uses a public formal.
   **The bcf-equivalence harness needs NO edit, and that is why the neutrality
   compare stays runnable.** `benchmarks/R/bcf-equivalence.R` drives
   `dbarts:::bartcoreBCFSampler` throughout (its header says so at `:11`); its
   `dbarts()` calls build a plain single-forest HOST and pass no BCF argument;
   `n.trees.treatment =` and `moderators =` there are that internal
   constructor's formals, not `dbarts()`'s. `settingsList()` (`:386-395`)
   carries only `quick`, `n.threads`, `n.burn`, `n.samples`, `n.trees.mu`,
   `n.trees.tau` and the `seeds` vector - no creation vocabulary at all - and
   is re-checked against the baseline meta at compare (`:419-420`), so it stays
   byte-identical with zero work. Scenario names and seeds are untouched.
   Consequence for FS1: the harness is NOT the vehicle. FS1 runs in
   `test-bcf-creation.R`, where the ORACLE is that same internal route, exactly
   as the shipped F1 pin does today (`:33-70` positive, `:76-93` negative) -
   so FS1 needs nothing "reconstructed from the pinned tests": it re-points the
   existing public-versus-internal bitwise pin at `forests =`, keeps
   `control@rngSeed = 17L`, and keeps the unseeded negative half (the internal
   route builds and discards a host engine, drawing `n.chains` `unif_rand()`s
   off R's stream, so the two must differ). The six channels are `train`,
   `sigma`, `varcount`, both forests' `bartcoreForestFits`, and
   `bartcoreBCFGlue`.

**10. The refusal set, counted, and the one refusal M2 RETIRES (added
2026-08-13, pre-M2).** "The eleven S1 refusals" in item 4 is stale: eleven was
the PRE-REGISTERED count, and S1 landed 17 firing assertions
(`bcf-public-surface.md:892-895`). At HEAD the shipped set is the 16-entry
`unsupported` vector (`R/spec.R:391-420`: DART prior, `split.probs`,
`monotone`, linear node prior, GP node prior, `k` hyperprior, non-default `k`,
non-default `node.scale`, named `prior.scale`, non-default `proposal.probs`,
Student-t residuals, grouped random effects, `variance`,
`storage = "single"`, per-column `n.cuts`, test predictors) plus the
gaussian-family gate
(`:375-384`) - 17 inherited unchanged. The eighteenth, the orphan-argument gate
(`:455-460`, "'moderators' and 'treatmentForest' configure the treatment forest
a 'treatment' vector selects"), is RETIRED by M2 along with the arguments it
guards; its replacement is whatever `forests =` needs in its place, which under
item 1's one-argument shape is nothing - there is no second argument left to
orphan. Say so in the slice rather than porting a dead message.
   **The knob map, walked and SETTLED (decided 2026-08-13 pre-M2, orchestrator
   discretion under the standing grant, VD may veto; item 1's earlier
   eight-name sketch was not total).** `treatmentForest()` carries TEN knobs
   (`R/model.R:1238-1249`). Each one's landing place, in order:
   `n.trees` -> `forest(n.trees =)`, per forest; `base` -> `forest(base =)`,
   per forest; `power` -> `forest(power =)`, per forest; `sd.control` ->
   `forest(sd =)` on forest 0 and `sd.moderate` -> `forest(sd =)` on forest 1
   (they are per-forest leaf scales, so the pair collapses into one per-forest
   name); `b.prior.variance` -> `forest(amplitude.prior.variance =)` on the
   BASIS forest; `update.a` -> `forest(update.amplitude =)` on forest 0 and
   `update.b` -> `forest(update.amplitude =)` on forest 1 (again one
   per-forest name for a pair that was only ever per-forest);
   `interactions` -> `forest(interactions =)` and `blocks` ->
   `forest(blocks =)`, per forest. The removed `moderators =` lands as
   `forest(vars =)` on the basis forest, completing the three removed formals.
   The two new names match the fork-3 `forestAmplitudes` vocabulary, and the
   whole map is a mapping rather than a new mechanism: the payload at
   `R/spec.R:429-454` takes the same eight doubles in the same order.
   **The legality rule, and where the "only with `basis =`" form had to be
   split.** `amplitude.prior.variance` is the `N(0, .)` variance of the
   `b0`/`b1` glue alone, so it is LEGAL ONLY on a forest given `basis =` and
   REFUSED at creation on a forest without one - forest 0's amplitude `a`
   carries the half-Cauchy scale-mixture prior, which the channel spec above
   fixes as an engine choice with no caller slot. `update.amplitude` cannot
   take that rule literally without becoming unreachable: `update.a` is forest
   0's flag, forest 0 takes NO explicit `basis =` (item 4 refuses one there,
   its basis being the implicit intercept), so a strict basis-gate would make
   `update.a` unexpressible and break M0's FD4 - the very failure this decision
   exists to prevent. It is therefore legal on EVERY forest of a `forests =`
   spec, which is the same rule read against the thing that actually gates it,
   an amplitude rather than an explicit basis, and refused where no amplitude
   channel exists at all: a K = 1 `forests =` spec. Both refusals are
   reachable and both are new (item 4). FD4's public fit is expressible and
   must be PINNED positively:
   `forests = list(forest(update.amplitude = FALSE), forest(basis = ~
   factor(z), update.amplitude = FALSE))` reproduces today's
   `treatmentForest(update.a = FALSE, update.b = FALSE)`.
   **One disposition falls out of the walk rather than being invented here.**
   `interactions`/`blocks` are per-forest above, but `dbarts()` ALSO has
   top-level arguments of those names, and today they are forest 0's
   constraints (`mu.interactions`/`mu.blocks`, `R/spec.R:441-452`) while the
   constructor's are forest 1's. Under `forests =` the two spellings address
   the same slot for forest 0, so supplying both is refused as ambiguous - the
   third new assertion in item 4. The top-level pair keeps its meaning
   untouched on the single-forest path, so FS3's byte-neutrality is unaffected.

**11. `dbartsData(treatment = )` is a FOURTH site, and it STAYS AS-IS at M2 as
recorded data-side debt. SETTLED 2026-08-13 pre-M2 (orchestrator discretion
under the standing grant, VD may veto).** Item 3 names `dbarts()` and
`dbartsSpec()`. `dbartsData()` has its own `treatment` formal (`R/data.R:697`,
documented `man/dbartsData.Rd:13`, `:16`, its own pkgdown topic
`_pkgdown.yml:28`), it is where `subset` alignment happens
(`validateTreatment(treatment, initialNumObservations, subset)`,
`R/data.R:962`, `:1028`), and it is REACHED BY `dbarts()`'s own redirect
(`R/dbarts.R:539`), so the two could not be decided separately. Two verified
constraints forced the disposition: the SLOT `data@treatment` cannot move at M2
(the bridge reads it by name, `src/R_interface_bartcore.cpp:1112`, and M2 is
no-src), and `dbartsSpec(dbartsData(x, y, treatment = z), control)` is a
shipped, tested creation route (`test-bcf-creation.R:272`) resolving off the
slot rather than off a `dbartsSpec()` argument (`R/spec.R:375`).
   The settled path, in four sentences the implementer can build to.
   `dbarts(forests = ...)`'s resolution populates `data@treatment` INTERNALLY,
   with exactly the values and at exactly the point `treatment =` populated it
   (`R/spec.R:373`), so the bridge, the payload cross-check and the R5 mirror
   of item 8(a) all see what they see today - which is also what makes FS1
   bitwise rather than merely close. The public `dbartsData(treatment = )`
   formal, its Rd entry and its pkgdown topic REMAIN, untouched and
   undeprecated; item 6's bookkeeping covers `treatmentForest` only.
   `dbarts()` and `dbartsSpec()` still lose their three formals, so the causal
   vocabulary leaves the two FITTING functions in this slice, which is what
   binding decision 3 sequenced and what bartCause's `bcf()` needs. The debt is
   the remaining split - a `forests =` model whose column still arrives under a
   `treatment =` name - and it retires at a later slice, the candidates being
   M3, M4 (whose n x q_f basis is the honest carrier and would absorb it whole)
   or the reshape re-bake window; the M2 RECORDS commit opens the TODO ticket
   for it, and this plan does not edit `TODO` itself.
   No escalation: the split surface is shippable. It is one exported
   data-constructor argument naming a column, against a fitting surface that no
   longer names one, inside a pre-release window where nothing is frozen
   (binding decision 5) - and building the basis-shaped carrier here instead
   would buy a second refactor over the same code when M4 lands the real one,
   which is fork 1 option (b)'s recorded failure mode.

**12. Budget, RESTATED with the migration and the bookkeeping included (added
2026-08-13, pre-M2; supersedes the ~180 R + ~60 man + ~260 test committed
2026-08-11).** ~775-895 dense-equivalent total.
   **R ~200-240**: `forest()` and its validator (`treatmentForest()` +
   `resolveTreatmentForest` are 27 + 31 lines today and are re-skinned rather
   than rewritten, but gain `basis`, `vars` and the three amplitude knobs of
   item 10); the `forests =` list resolver (arity, per-forest dispatch, the
   K > 2 and forest-0-basis refusals); the basis EXPANSION - evaluating a
   formula in the data environment, expanding a factor to level indicators, and
   aligning it with `subset` - which is the one genuinely new mechanism in the
   slice and has no ancestor in the `treatment =` path; the two R5 renames
   (`$setTreatment`'s 28-line body plus `$getBCFGlue`'s 5); and the removals
   across `R/dbarts.R`, `R/spec.R`, `R/model.R`.
   **man ~120-160 raw**: a new `man/forest.Rd` in the shape of
   `man/treatmentForest.Rd` (92 lines) with `basis`/`vars` added and the
   amplitude knobs documented; the three `\item` blocks and the usage line of
   `man/dbarts.Rd` collapsing to one; `man/dbartsSpec.Rd`; five sites in
   `man/dbartsSampler-class.Rd` (two aliases, two usage lines, the mutator
   list, `\item{z}` -> `\item{basis}`, `\item{forest}`'s fourth clause,
   `\value`) plus M1's two `treatment = ` phrases (item 8(b)); `inst/NEWS.Rd`.
   The dense single-paragraph convention applies, as it did to M1's 4 raw
   lines - these are raw line counts, not prose volume.
   **test ~435-475**: `test-bcf-creation.R` rewritten (385 lines today, and it
   GROWS - item 4 adds NINE refusals to the 26 assertions already there); FS1's
   two halves re-pointed at `forests =` with the internal-route oracle
   unchanged; FS3's byte-neutrality arm; and the seven-file migration of item 9
   (14 creation calls, 5 R5-method call sites, ~40 lines of edits carrying no
   new assertions). The settled knob map moves this by about +15 over the
   ~420-460 first written on 2026-08-13: its three new refusals ride the
   rewritten file's refusal matrix at no extra cost, since item 4 and FS2
   already own that matrix, but item 10's FD4-expressibility arm is NEW - a
   positive pin that `forest(update.amplitude = FALSE)` on both forests
   reproduces `treatmentForest(update.a = FALSE, update.b = FALSE)`, without
   which a knob can be accepted and dropped in silence.
   The calibration-S2 lesson applies, as it did at M1: this budget moves with
   the obligations, and a stale one beside added obligations is the failure
   mode.

Falsifiers: **FS1, load-bearing** - with `control@rngSeed` SET, a `forests =`
BCF is BITWISE IDENTICAL on all six bcf-equivalence channels to the
`dbarts:::bartcoreBCFSampler` oracle the shipped F1 pin already compares
against (item 9 - the removed `treatment =` route needs no reconstruction, and
the phrase "reconstructed from the pinned tests" is struck); NEGATIVE HALF, two
arms: perturb the factor expansion (swap the two indicator columns) and it must
go red, and the shipped unseeded arm must keep differing.
**FS2** the refusal list of item 4 as corrected by item 10, one assertion each,
each raising an R condition with no handle escaping. **FS3** `forests = NULL`
is byte-neutral: a single-forest fit is bitwise unchanged.
rng: NEUTRAL on every existing path. Gates: `R CMD INSTALL --preclean` into a
private library; full tinytest with NO snapshot regenerated; `R CMD check` on a
clean-copy tarball (NEW relative to M1: M2 moves `NAMESPACE`, so
undocumented-export and codoc are live); `air format --check .`; `lintr` on
touched R; `pkgdown::check_pkgdown` (NAMED, not "pkgdown check" - it is the gate
that catches item 6's missing `_pkgdown.yml` entry); the trio expected BITWISE
against baselines `equivalence-8b047f8b` (37 scenarios),
`bcf-equivalence-8b047f8b` (12) and `multinomial-equivalence-1027be5` (10), all
three verified current in `benchmarks/baselines/MANIFEST` at 4c7aa200. ABORT on
any trio divergence: M2 is R-only against an untouched engine and an untouched
internal creation route, so a divergence is a LEAK in the new resolver, never a
re-record.

## M3. The flat mean channel. An item INSIDE dbarts-h-reshape S1.

**LANDED at ab3aa2fa** (reshape S1, the window's second and last re-bake); the
four items below are a record of what shipped, not forward work. Every
`dbarts.h` line number in this section and in the cross-plan block below is
PRE-RE-BAKE and stale; the drift is swept once in the `amended:` note at the
top of this file rather than by editing text that was applied verbatim.

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

**Everything below this line is AMENDED 2026-08-13, pre-M4** (artifacts
`.claude/m4-basis-design/pre-m4-scoping-2026-08-13.md` and
`critique-2026-08-13.md`; adjudications in `synthesis-2026-08-13.md`). The
section was the last place in this plan still speaking the retired pre-fork-3
vocabulary; it now uses the shipped names - `basis` (`forest(basis =)`,
`$setForestBasis`, `dbarts_sampler_setForestBasis`) and `amplitude`
(`forest(amplitude.prior.variance =)`, `forest(update.amplitude =)`,
`$getForestAmplitudes()`, `numForestAmplitudes`/`forestAmplitudes`).
`setTreatment` and `bcfGlue` are GONE as flat C names and no bridge message
spells `treatment` (reshape S1's landing note). M3 no longer exists as a
slice: it landed INSIDE reshape S1, so "the M2 and M3 guards relax in
lockstep" below means the guards now sited in M4.3's own list.

**Sequence, amended.** The committed order put a VD fork before the probes on
a cost argument that does not survive: `benchmarks/baselines/
bench-sampler-235bebc.csv:3` records 0.47 ms/sweep at
`run-n1000-p10-t200`, exactly FA1a's shape, so FA1a's 32 fits x 8000 sweeps is
~3 MINUTES of engine time, not the hours the scoping pass estimated;
`bcf-ridge-interweaving.md:481-483` independently records 0.36 ms/sweep for
BCF at n = 500. The probes are minutes, so they run FIRST and nothing forks
ahead of them. **Amend the plan -> run FA1a, FA1b and the revised FA5 ->
re-scope M4.2 and M4.4 on the verdicts -> M4.0, M4.1, M4.2, M4.3, M4.4,
M4.5**, serialized, one implementer, each slice landing before the next
starts. **STATUS, re-scoped 2026-08-13 on the probe verdicts
(`probes-2026-08-13.md`): steps 1 to 3 are COMPLETE** - the plan was amended
(the PRE-M4 block), all three probes ran in the foreground at 47c1fbe1 (no arm
crashed, none was UNOBTAINABLE), and this commit is the re-scope. **M4.0
LANDED 562ee684**, whose pin scope - BOTH `afterCombine` overrides - the
verdicts left unchanged, per the Landing notes below. **M4.1 LANDED
1458328c**, follow-up e48fc5de. **M4.2 LANDED 1a2aaedc**, amended once after
independent review. **M4.3 LANDED 9c63e9d8**, amended once after independent
review (LAND-AFTER-CHANGES, four dispositions applied). **The next step is
M4.5, then M4.4 (hatch order).**
**The slice order is now M4.0 -> M4.1 -> M4.2 -> M4.3 -> M4.5,
Gaussian-complete, with M4.4 as the IMMEDIATE follow-on slice**, still
pre-release (fork 2 holds both slices inside the window; the escape hatch
REORDERS, it does not cancel). FA5 licensed that hatch and it is TAKEN; see
M4.4 and FA5's teeth. The probes are STATISTICAL (IACT and posterior agreement,
both load-independent), so unlike `bench-sampler.R` they needed no quiet
machine and could share the box. Their harnesses live GITIGNORED under
`.claude/m4-basis-design/harness/` per the house convention
(`composition-mixing-probe.md:13-16`, `grow-from-root-default.md:6,55`);
the binding corollary is that their NUMBERS are carried into this section
verbatim, never by reference to a directory that does not arrive with a clone
(`multiforest-predictor-mutation.md:71`). PROMOTE FA1 to
`benchmarks/R/basis-amplitude-mixing.R` if and only if its verdict re-scopes
M4.2 - that is exactly the condition `benchmarks/R/grouped-mixing.R` was
promoted under. FA5's numbers are a one-time adjudication and stay gitignored.
**NO PROMOTION (2026-08-13):** FA1's verdict does NOT re-scope M4.2, so the
condition is not met and `fa1a-basis-amplitude-mixing.R` /
`fa1b-composition-scaling.R` stay gitignored with FA5's; the numbers recorded
in FA1 and FA5 below ARE the record.

**Four candidate VD items, ALL resolved WITHOUT a VD question (orchestrator,
2026-08-13, under the standing grant; VD may veto any of them).**

1. **The b-ridge fork is DISSOLVED - there is no fork.** It was to be put on
   the b-move's "shot at the open `bcf-sigma-residual` flag". That flag is
   RESOLVED (`docs/plans/bcf-sigma-residual.md:2-6`, `INDEX.md:20`; no `TODO`
   entry survives) and the b-move was EXONERATED as its carrier by the memo's
   own controls. Both halves of the trade are gone, so the item is a design
   consequence to record rather than a decision to take. See M4.2. **No
   `bcf-sigma-residual` ticket is opened or re-opened by M4.**
2. **The amplitude-free default was PROBE-CONTINGENT on FA1. CLOSED
   2026-08-13: NO amplitude-free default.** The pre-registered antecedent was
   that FA1a's WITHOUT arm match from its best scale cell AND FA1b show no
   K-scaling degradation; **both conjuncts FAIL** (FA1a: WITHOUT is 1.08x-1.72x
   worse on IACT on 7 of 8 functionals from its own best cell and 12-17% worse
   on muRMSE in every cell, winning 0 of 8 seeds; FA1b: per-forest IACT grows
   2.8x and sigma IACT 3.0x from K = 2 to K = 8 - numbers and protocols in FA1
   below). So the general path does NOT ship amplitude-free, amplitudes are NOT
   demoted to opt-in, and **M4.2 keeps the q-variate amplitude conditional and
   the per-forest ASIS rescale in FULL scope**. The saving "worth roughly
   M4.2's whole engine budget" is NOT available; the re-priced budget below
   never took it, so it stands unchanged. Reported to VD as a fact, no
   question.
3. **M4.4's justification was PROBE-CONTINGENT on the revised FA5. CLOSED
   2026-08-13: the headline ground FALLS and the ESCAPE HATCH IS TAKEN.** Arm B
   - two single-forest GAUSSIAN samplers with the latents drawn in R against the
   COMBINED fit - AGREES with the reference on all 12 functionals (max
   |z| = 2.54 against a threshold of 3.0), with the power precondition MET (arm
   A, the strawman, differs at max |z| = 772, 12 of 12 over 4). So "the caller
   cannot compose non-Gaussian latents against the combined fit" is FALSE, and
   M4.4 is re-justified on the surviving grounds only (see M4.4). Fork 2 had
   already licensed both branches, so this selects among options VD holds: the
   arc ships Gaussian-complete through M4.5 and M4.4 follows IMMEDIATELY, both
   pre-release. Reported, not asked.
4. **A per-draw amplitude channel in flat C is DECLINED at orchestrator
   discretion.** `dbarts_results` (`dbarts.h:138-152`) has no amplitude member
   and the header contains zero "glue" tokens, while the R bridge does report
   per-draw amplitudes - a real R/C asymmetry, and the exact
   binding-decision-1(a) inversion M1 existed to close. But it PREDATES M4,
   M4 does not need it, and closing it would append to `dbarts_results`: not a
   hash move (struct layout is unhashed, `dbarts.h:95-99`) but a
   `DBARTS_C_API_MINOR` 0 -> 1 bump by the header's own append rule
   (`:122-123`), plus a break of the `sizeof(dbarts_results)` assert at
   `C_interface.cpp:274`. Binding decision 8 pins both constants at 1/0 through
   the window and forbids it. Recorded as a DOOR on the TODO, not a fork.

- **M4.0 (tests only).** Pin today's combiner at the seam that generalizes:
  `forestMultiplier`'s two values (`combiner.hpp:762`, private, call sites
  `:552` and `:560` - so the pins live in `tests/cpp`, not tinytest);
  `combinedFits`' blend (virtual `:339`, BCF override `:567-576` reading
  `forests[0]`/`forests[1]` by literal index); and the draw order, which the
  committed text folds into one method and which is TWO. `drawGlue`
  (`:582-623`) draws `a` (`:598-599`, standard normal), `aVariance` (`:606`,
  gamma), `b0` (`:620`), `b1` (`:621`) - and nothing else. The ridge draw is
  the GIG at `:672` inside `afterCombine` (`:638-709`), which `Chain` calls
  immediately after `drawGlue` (`chain.hpp:1280-1281`). Pin BOTH, in that
  order, through the hook that already exists for exactly this -
  `Chain::interweaveGlueRidge(record, sampleNum)`, `chain.hpp:1063-1067`,
  documented "for the component tests". **The pin scope covers BOTH live
  overrides of `afterCombine`, not just BCF's** (amended 2026-08-13):
  `MultinomialForestCombiner::afterCombine` (`combiner.hpp:1152-1197`) is a
  second implementation of the same virtual in the same sweep slot, drawing an
  ADDITIVE level-centering shift that consumes RNG at `:1180-1181`
  (`ext_rng_simulateStandardNormal`), applying it to every forest's
  `totalFits` and `muByTree` at `:1183-1195`, and returning 1.0
  unconditionally at `:1196`. M4.2 redefines that virtual; its multinomial
  implementation must be pinned on RNG consumption and return value before
  that happens, or the arc rewrites a virtual with half its implementations
  unpinned. Also pin BCF's `1.0`-on-every-skip returns (`:640`, `:662`,
  `:674`, `:676`) and its returned `c` (`:708`).
- **M4.1 (engine, byte-identical).** `forestMultiplier`
  (`combiner.hpp:761-765` - the whole K = 2 basis map is five lines) becomes
  `dot(a_f, B_f[i,])`; `combinedFits` (`:567-576`) becomes a K loop.
  `formForestResponse` (`:546-565`) needs NO change: it already loops every
  forest through `forestMultiplier`, and it is generic in q as well as K - the
  amplitude vector contracts with the basis row to a per-row SCALAR, so
  `resid / m`, `w m^2` and the `0x1p-26` snap (`:525`, `:553`) all survive a
  wide basis with the condition-number cap argument at `:538-541` intact. That
  is this plan's own "generalization rather than a second mechanism" claim, and
  it holds by symbol. Ship with BCF's basis synthesized internally from z so
  BCF is bitwise unchanged - which only stays true if the R route keeps
  producing exactly today's column, so `expandForestBasis`
  (`R/model.R:686-708`) must not move here. This is the whole architectural
  risk, isolated. ABORT on any divergence in the gating trio, whose CURRENT
  baselines are `bcf-equivalence-6e3b9fb8.rds` (12 scenarios,
  `benchmarks/baselines/MANIFEST:42`), `equivalence-8b047f8b.rds` (37
  scenarios, `MANIFEST:16`) and `multinomial-equivalence-1027be5.rds` (10
  scenarios, `MANIFEST:49`) - the last of these is a GATING baseline for M4.1
  and M4.2 alongside `bcf-equivalence`, because the multinomial combiner is a
  live K-forest model overriding the same virtuals.
- **M4.2 (engine). Scope CONFIRMED IN FULL, re-scoped 2026-08-13 on the probe
  verdicts** (`probes-2026-08-13.md`): FA1's conjunction fails on both
  conjuncts, so neither the amplitude conditional nor the ASIS rescale is
  demoted to opt-in and no BCF-specific carve-out is taken on that ground. The
  q-variate amplitude conditional, using the shipped
  Cholesky helpers (`model.hpp:869` `choleskyDecompose`, `:886`
  `solveLowerTriangular`, `:895` `solveLowerTriangularTransposed`; in-place,
  row-major p x p, NO failure path because callers guarantee definiteness - so
  the conditional must supply its own ridge/prior term to keep the crossproduct
  PD). `afterCombine` becomes a per-forest ASIS rescale of the whole amplitude
  vector along the likelihood-invariant orbit, and it must RE-STATE the base
  virtual's contract rather than inherit it: the Doxygen at
  `combiner.hpp:425-427` ("returns the scale its move applied (1.0 when it
  makes none)") is already false for the multinomial implementation, which
  makes an additive move and returns 1.0 by convention (`:1261`, `:1306`).
  BCF's K = 2 orthogonal case must come out bitwise, or its specialized
  two-scalar path is kept explicitly and the general path is gated
  statistically - that in-slice decision mechanism is unchanged. This is the
  one place in the arc that could force a `bcf-equivalence` re-record; say so
  in the slice's rng note rather than discovering it. **Standing rule for any
  re-record:** it bumps `.github/workflows/equivalence.yaml`'s baseline lines
  in the SAME commit (`:61` equivalence, `:87` bcf-equivalence, `:113`
  multinomial-equivalence). A re-record that leaves the workflow pointing at
  the old file passes locally and fails CI - or worse, passes both against a
  superseded baseline.

  **The containment, recorded (amended 2026-08-13).** The general per-forest
  ASIS applied at the TREATMENT forest IS the b-move that
  `docs/plans/bcf-b-ridge.md` derives - not a second mechanism M4.2 adds, but a
  design consequence of the generalization. The memo's own general rule makes
  the identity exact: "rescaling k leaf params by c against d glue scalars by
  1/c gives p = (k-d)/2 (a: k=L, d=1; b: k=Lt, d=2)" (docs/design/multiplier-
  combiner.md, "The exponent rule" - moved there from this memo's `:192-194`
  since), which IS M4.2's q-variate exponent, with `B = ||a_f||^2 / priorVar_f`
  for a fixed-variance amplitude prior and `a^2/aVariance` for the scale
  mixture (same section, moved from `:183-196`). Instantiated at BCF's forest
  1 (q = 2, `bPriorVariance` fixed at `combiner.hpp:300`) the general move is
  the b-move exactly. Three consequences:
  - **The sigma-motivated build stays REJECTED.** `bcf-b-ridge.md`'s own
    controls settle it: `:372-381` (control (a), fixed-b) "EXONERATES the
    b-ridge as the carrier of the residual sigma bias: the a-only memo's
    'prime suspect' is WRONG"; `:396-406` (control (b), burn-2x) "The residual
    sigma bias is a BURN TRANSIENT, confirmed"; `:420-421` "ROUTE: FIX
    INITIALIZATION / BURN, do NOT implement the b-move for the sigma flag";
    `:452-454` "Do not gate the b-move on sigma, and do not expect the b-move
    to substitute for the burn/init fix." `bcf-sigma-residual` is RESOLVED
    (that plan's `:2-6`, `INDEX.md:20`); the only survivor is
    `TODO`'s `bcf-sigma-tail-mixing` OPTIONAL research door, which routes to
    tree-structure mixing at high SNR, not to a glue move. **M4.2 opens no
    sigma ticket and claims no sigma payoff.**
  - **`bcf-b-ridge.md` is an M4.2 INPUT, and prices the slice DOWN.** It
    already delivers, for free, three of M4.2's most expensive deliverables:
    the general `p = (k-d)/2` exponent (docs/design/multiplier-combiner.md,
    "The exponent rule") with the naive move-map shortcut identified as off by
    one (same section); a PASSED adversarial R prototype across three `Lt` and
    two parameterizations, whose DISCRIMINATION arm rejects the off-by-one
    exponent at KS = 1.6e-21 (same section); and the complete
    rescale-consistency checklist - `treeFits`, `totalFits`,
    `totalTestFits`/`currTestFits`, and the keepTrees flattened-slot sharp
    edge (`bcf-b-ridge.md:246-289`) - each of which the shipped `afterCombine`
    implements at `combiner.hpp:680-707` and each of which the q-variate
    generalization must re-derive. Cite docs/design/multiplier-combiner.md's
    "The exponent rule" and bcf-b-ridge.md section 4 as the derivation and
    prototype, not merely as evidence that a move is unbuilt.
    The GIG generator (`ext_rng_simulateGeneralizedInverseGaussian`) already
    ships.
  - **The move's acceptance is its OWN, per that memo's `:438-449`:** judge it
    on its own IACT payoff (`|b1-b0|` or tau-amplitude IACT on a
    strong-treatment-signal DGP, mirroring the a-move's `|a|` IACT 2.5x), with
    correctness acceptance `bcf-exact` mode-2b stays exact AND a keepTrees BCF
    round-trip tracks. This is FA2's q-variate half, run as an M4.2-internal
    gate (see the falsifiers).

  **Handoff from M4.1 (LANDED 1458328c/e48fc5de), recorded at the M4.1
  landing.** Four items MANDATORY before M4.2 starts, each verified against
  the live file:
  (a) Where this section says "BCF's K = 2 orthogonal case must come out
  bitwise": the summation-association constraint is now EXPLICIT, not
  implicit - a re-association in `combinedFits` breaks
  `bcf-equivalence-8b047f8b` on all 12 scenarios SILENTLY past the seam pins
  (`combiner.hpp:617-632`, the load-bearing reverse accumulation at
  `:624-629`); `testCombinedFitsAssociation`
  (`tests/cpp/test_sampler.cpp:3206`) is the ONLY in-process guard.
  (b) `BCFState`'s `a()`/`b0()`/`b1()` accessors hard-code amplitude indices
  0/1/2 (`combiner.hpp:330-335`) - correct today, but they alias the wrong
  forest's block the moment `q0 != 1`. M4.2 must move to per-forest offset
  indexing (`glue_.amplitudeOffset`, already carried for exactly this).
  (c) THREE accumulation conventions now coexist: `combinedFits`' reverse
  accumulation (`:624-629`), `forestMultiplier`'s dot product forward
  (`:833-835`), `formForestResponse`'s residual forward (`:592-593`). Only
  the first is documented load-bearing; the other two become observable at
  q > 1 / K > 2. M4.2 must decide and document each.
  (d) `installTreatment` pins `basis.resize(2)` (`:855`) while `combinedFits`
  is general in K - a K = 3 BCF chain would read `basis[2]` OOB.
  Unconstructible today; M4.2 generalizes the sizing with the K-length spec.
- **M4.3 (spec, factory, refusal relaxations). RE-SCOPED WHOLESALE
  2026-08-13** - most of what it was written to design has SHIPPED, and most of
  what remains it priced at zero. **AMENDED 2026-08-14**: the amendment block
  at the end of this slice discharges both entry preconditions (the
  interaction spec, the test re-price), supersedes the tests paragraph and the
  M4.2 handoff item below, and its anchor list corrects every stale
  `combiner.hpp` citation in this section. **M4.3 LANDED 9c63e9d8, and
  everything in this bullet is the PRE-SLICE record** - written against
  `2e538430`/`e7708b7c`, true then, not retro-edited since. The landing note at
  the end of this file is what is current, and where M4.4's section 4 (every
  anchor re-derived by symbol at `efec6ba2`) disagrees with anything here,
  section 4 governs.

  **Already shipped, strike it from the slice:** the public creation route
  `forests = list(forest(...))` with the complete ten-knob map (M2), so M4.3
  relaxes a `length(forests) > 2L` test rather than designing a surface; the
  per-forest amplitude priors `forest(amplitude.prior.variance =)` and
  `forest(update.amplitude =)` with their settled legality split; the whole
  flat C mean channel (`setForestBasis`, `numForestAmplitudes`,
  `forestAmplitudes`) with general, staged-guard signatures; the R5 mirror
  `$setForestBasis(forest, basis)` and `$getForestAmplitudes()`, 1-based via
  `resolveForestIndex` (`R/bartcore.R:1051`) since a14040de; the refusal
  discipline (M2's FS2 landed 42 message-asserting refusals; FC2 landed the
  flat C refusal matrix in `consumer.c` with an `18L` `LEG_COUNT` pin).

  **What remains, by site.** R creation, `R/model.R`: `:772-776` the `> 2L`
  arity refusal; `:650-657` `forestBasisDeclaration` (hardcodes
  `length(forests) != 2L`, must return a per-forest basis LIST); `:686-708`
  `expandForestBasis` (the whole function is the K = 2 expander - the
  factor-only guard at `:693` with its message at `:695-696`, "exactly two
  levels" at `:699-702`, the single 0/1 column collapse at `:708`, plus `:670`
  "a 'basis' formula must be one-sided"); `:784-818` `resolveForests`'s six
  positional refusals, each written against "first forest"/"second forest";
  `:840-854` `forestParams`, whose positional length-8 write (`:844-853`,
  `declared(tau$sd, 1)` at `:849`) becomes a ragged per-forest transport -
  **this is the single most K = 2-shaped thing in the arc**; `:622-645`
  `resolveModerators`, one column-index vector rather than per-forest. R
  creation, `R/spec.R`: `:412-419` the gaussian-only refusal (M4.4's R-side
  gate); `:429-456`'s SIXTEEN named composition refusals emitted at `:457-463`
  - M4 relaxes the family one, the other fifteen and their pins STAY;
  `:466-486` the `mu.`/`tau.`-prefixed control attribute. R data layer:
  `R/data.R:643-655` `validateTreatment`'s 0/1 coding check, the K = 2 pin
  every route reaches. R5: `R/dbarts.R:1218-1223` `$setForestBasis`'s
  `index != 1L` refusal and `:1237` its single-slot `data@treatment` mirror.
  Bridge (`src/R_interface_bartcore.cpp`): the length-8 params check at TWO
  sites (`:2134-2135`, `:2913-2914`) and the positional unpack `:2138-2148`;
  `applyBCFSpec`'s two-forest shape `:2120-2126` and `applyBCFAttributes`
  `:2200-2205`; the `data@treatment` read by name `:1112` and the 0/1
  coercion scans `:2843-2845`, `:2933-2934`; the two-way creation cross-check
  `:2831-2836`; `double amplitudes[3]` at `:3673` and the stride-3
  `Rf_allocMatrix` at `:3769-3772`. Flat C (`src/C_interface.cpp`): guard
  bodies only as designed (`:768`, `:773`, `:775-778`, `:780-787`, `:790-791`,
  `:814-815`, `:828-832`), plus the THREE `double amplitudes[3]` stack buffers
  at `:767`, `:812`, `:824`, which a K-amplitude model overflows - a
  correctness item, not a signature item. Engine: `BCFForestSpec`
  (`combiner.hpp:211-241`, WIDER than this plan describes - it now carries
  `columns`/`numColumns`, the three interaction fields and the three block
  fields, so a K-length vector of these is a bigger object than priced) ->
  the K-length entry type; `BCFSpec` (`:244-251`) -> the thin adapter;
  `BCFState`'s `a, b0, b1, aVariance` (`:290-302`) and the same four in
  `ChainStateData` (`:92-94`) -> per-forest amplitude vectors; the BCF
  constructor `chain.hpp:692-723` with its two asymmetric `buildBCFForest`
  calls (`:716-719`) and `buildBCFForest` itself (`:4454`); `facade.hpp:272-274`
  still spelling the retired `setTreatment` in the engine vtable, and
  `createBCFSampler` (`:750-762`). The constant-leaf `static_assert`
  (`combiner.hpp:502-503`) STAYS.

  **Six collisions this slice collides with and priced at zero.**
  (a) **The `bcfGlue` out-parameter is a FOUR-LAYER vtable, not a bridge
  concern:** base virtual `combiner.hpp:386-389` `bcfGlue(double&, double&,
  double&)`, `Chain::bcfGlue` `chain.hpp:1018-1020`, `Sampler::bcfGlue`
  `sampler.hpp:1242-1244` unpacking `out[0..2]`, `facade.hpp:315` carrying no
  length at all. Re-signing it is ENGINE work.
  (b) **The state format's `bcf` block.** Four fixed scalars serialized at
  `R_interface_bartcore.cpp:5812-5818`, read with an `Rf_xlength(bcfExpr) != 4`
  check and a "malformed bcf glue in bartcore state" message at `:6163-6167`.
  A ragged amplitude vector is non-additive to that encoding. **RE-ENCODE
  "bcf" IN PLACE and bump nothing.** The registry's append-only rule
  (`:5645-5655`) is real but the paragraph after it waives this case:
  `:5658-5665` records that version 1 is the FIRST shipped format and "the
  pre-release development increments (the forests-list/BCF restructure
  included) are collapsed into it - no released reader or writer ever saw
  them. Pre-1.0 states are not a compat target." The "bcf" block was itself
  ADDED by exactly such a pre-release restructure. Do NOT add a parallel
  optional block: that would carry a dead K = 2 block into the released 1.0-0
  format permanently, plus read-side precedence rules for four states.
  (c) **The run-result `glue` slot is 3-row and pinned.**
  `R_interface_bartcore.cpp:4100-4102` allocates `3 x numSamples` /
  `3 x numSamples x numChains`; pinned at `test-bcf-reporting.R:41` and `:114`
  with the reconstruction identity at `:57-71`.
  (d) **`$getForestAmplitudes()` is the one shipped R5 signature M4 must
  BREAK.** Argument-free and documented "(a, b0, b1) ... as a 3 x n.chains
  matrix" (`R/dbarts.R:1442-1446`), pinned `test-bcf-r5-surface.R:44,50`.
  Under a ragged K-length vector it takes a `forest` argument (matching the
  flat pair) or returns a list.
  (e) **The `dbartsData(treatment = )` debt retires HERE and is entangled
  four ways.** `TODO:90-105` routes it to M4 ("whose n x q_f basis is the
  honest carrier and would absorb it whole") and it has already MANIFESTED -
  bartCause `dbarts-1.0@695c603`'s suite is red on it. Beyond data carriage
  the slot is load-bearing for: the bridge's read by name
  (`R_interface_bartcore.cpp:1112`), the R5 `$setForestBasis` mirror
  (`R/dbarts.R:1237`), the two-way creation cross-check (`:2831-2836`), and -
  decisively - `isBCFSampler <- !is.null(sampler$data@treatment)`
  (`R/bartcore.R:20-22`), on which the ENTIRE R5 refusal system depends
  (`refuseBCFMutation` `:30-34` and its five call sites: `setResponse` `:285`,
  `setOffset` `:306`, `setData` `:383`, `setModel` `R/dbarts.R:1044`,
  `setCalibration` `:1484`). **The successor is the CAPABILITY probe the
  shipped code already names**, generalized from "carries BCF glue" to
  "carries amplitudes": `R/bartcore.R:18-19` documents the slot test as "the
  R5-layer capability probe, cheaper than a round trip through the bridge's
  own (`Chain::bcfGlue`-based) one". **A `numForests() > 1` probe is REFUSED
  by the shipped code twice** - `R_interface_bartcore.cpp:3668-3670` ("a
  K-forest multinomial (K >= 2) defeats a numForests test") and
  `C_interface.cpp:764-766` independently. Under M4 that matters more, not
  less: the general family is a THIRD multi-forest combiner, so a forest-count
  probe misfires on two of three.
  (f) **`bcf-equivalence.R` drives the INTERNAL route** and needs no edit -
  `:11`, `:115`, `:130`, `:147` all call `dbarts:::bartcoreBCFSampler`, and
  `:95-101` read through `bartcoreForestFits`/`bartcoreBCFGlue`. That holds
  for M4 only if the "thin adapter" keeps that function's signature, which is
  why the adapter clause is LOAD-BEARING rather than courtesy: the bitwise
  gate runs THROUGH it. `bartcoreBCFSampler` (`R/bartcore.R:622`, internal) is
  also still the oracle for six test files.

  **Tests the relaxation rewrites** (existing pins, not new tests). RE-COUNTED
  2026-08-14 against the live files; this paragraph supersedes the count and
  the citations the slice was re-scoped with, and every line number here is the
  `expect_error(` ASSERTION line, where the superseded list cited MESSAGE
  lines (`:550`, `:563`, `:567`, `:579`, ...). `test-bcf-creation.R` carries 45
  `expect_error` (84 `expect_*`), not M2's 42: `git blame` attributes the three
  `"length of 'basis'"` pins `:321`, `:329`, `:333` to reshape S1 (M3), which
  parameterized `validateTreatment`'s message. All 45, by disposition:
  - **RELAX outright, 6** (the refusal becomes a positive route): `:556` "at
    most two forests", `:565` "not a numeric basis column", `:569` "exactly two
    levels", `:581` "first forest takes no 'basis'", `:590` "'vars' restricts a
    basis forest", `:599` "the first forest has no 'basis'".
  - **RESTATED for K forests, 1** (the refusal survives): `:656` "second forest
    needs a 'basis'".
  - **REWRITTEN by the treatment-slot retirement, 4**: `:316` "coded 0", `:317`
    "length of 'treatment'", `:808` "no treatment forest was configured",
    `:819` "carry no treatment vector".
  - **STAND, 34**: the FIFTEEN composition refusals `:426-548` (the gaussian
    one is M4.4's, not M4.3's); the 6 knob validators `:694-753`; `:754`
    one-sided; the 3 shape refusals `:666`, `:670`, `:679`; the 2 top-level
    ambiguity refusals `:624`, `:634`; the 2 K = 1 amplitude refusals `:611`,
    `:615`; the 2 bridge backstops `:777`, `:793`; the 3 `"length of 'basis'"`
    `:321`, `:329`, `:333`.

  The composition refusals are FIFTEEN as pinned, not sixteen: `R/spec.R:429-
  456` emits sixteen, but only 15 `expect_error` land in [426, 548], and the
  STAND row balances at 34 only at 15. So 11 of 45 pins move and 34 stand, plus
  4 positional-`params` assertions the ragged transport rewrites though they
  are not refusals: `:343-346`, `:367` (`params[4L]`), `:374` (`params[7L]`),
  `:383` (`params[1L]`).

  Other files. `test-bcf-r5-surface.R` (10 `expect_error`, 40 `expect_*`): the
  3-row `dim(glue)` pin `:50`, the `"only on the second"` refusal `:135`, the
  whole F3 mirror block `:153-175`. `test-bcf-reporting.R`: the dims `:41` and
  `:114`, and the `reconstructionError` helper `:57-71`, which reads
  `glue[1L,]/[2L,]/[3L,]` by position. `inst/tinytest/capi/consumer.c`: FOUR
  legs, not three - `:897` `setForestBasis` itself and `:898-900`
  `.forest0`/`.width`/`.values` - plus their bodies `:1003-1032` and the `18L`
  `LEG_COUNT`; TWO of the four INVERT, and an inverted refusal leg is exactly
  where a capability answer and a malformed answer can quietly swap channels:
  `LEG_BASIS_FOREST_0` (`:1006-1018`) asserts `setForestBasis(sampler, 0, ...)
  == 0` and becomes an ACCEPTANCE (forest 0 takes a basis; its forest 2 probe
  stays a capability answer on a K = 2 sampler), and `LEG_BASIS_VALUES`
  (`:1022-1032`) asserts a raise on a 0.25/0.75 pair that becomes LEGAL (see
  the amendment's draw-path item). Both leg bodies also re-lay their bases
  ROW-major under the transpose item. `test-capi.R:1165-1179` (the ragged 1/2
  amplitude counts); `test-bcf-mutation-pins.R` (6), `test-bcf.R:57` and
  `test-multi-forest-seam.R:285` follow the route rename.
  `test-bcf-zero-multiplier.R` never mutates a basis and is STRUCK from this
  list. In `tests/cpp` the pins whose CONTENT moves are THREE:
  `testForestBasisSynthesis` (`test_sampler.cpp:3296`, superseded by the
  amendment's arms), `testAmplitudeOffsetIndexing` (`:3670`, which calls
  `installForestBasis` directly and gains the surface route), and
  `testBCFTwoForest`'s `sampler.setTreatment` call (`:2198`) with its "BCF
  setTreatment refresh runs" message (`:2207`). The other 22 of the 25
  `BCFSpec` fixtures (`test_sampler.cpp` x 21, `test_fuzz.cpp` x 3,
  `test_model.cpp` x 1) stay untouched IF collision (f)'s thin-adapter clause
  is honored - which is why that clause is load-bearing for the C++ suite
  exactly as it is for `bcf-equivalence.R`.

  **Handoff from M4.2 (LANDED 1a2aaedc), recorded at the M4.2 landing.** ONE
  item MANDATORY before M4.3 starts: define the `setTreatment` x
  `installForestBasis` interaction. `setTreatment` silently discards a
  widened basis today - M4.2's review pass only documented the fact on
  `installTreatment`'s Doxygen, it did not resolve it. M4.3 wires the
  per-forest basis surface `installForestBasis` is the engine half of, and
  must specify which route wins, with a test. **DISCHARGED 2026-08-14 by item
  1 of the amendment below, with the test shape at item 2.**

  **AMENDED 2026-08-14**, on `.claude/m4-basis-design/pre-m43-2026-08-14.md`
  (the pre-slice scoping pass, read at e7708b7c) and
  `.claude/m4-basis-design/critique-m43-2026-08-14.md` (its adversarial
  critique, read at the same commit). Where the two conflict the CRITIQUE's
  evidence wins: it does at items 1, 3, 4, 6 and 7 below and in the tests
  paragraph above; items 2, 5, 8, 9 and 10 are the scoping pass as written,
  certified. Both of M4.3's entry preconditions are hereby discharged: the
  interaction spec (item 1) and the fresh test re-price (item 6). M4.3 may
  start.

  **1. The interaction spec: SUBSUME. `setTreatment` RETIRES as a mutator.**
  Basis synthesis becomes construction-only and `installForestBasis`
  (`combiner.hpp:607-613`) becomes the sole mutation route, at any forest, at
  any width. It wins by being the only operation there is.
  - `installTreatment` (`:1119-1142`) becomes the CONSTRUCTOR's private
    synthesis: `spec.z` -> `basis[0]` all-ones (`:1124-1125`), `basis[1]` the
    (1 - z, z) pair (`:1127-1134`), forests 2.. all-ones (`:1136-1139`).
    Nothing outside the constructor calls it. Today it rebuilds BOTH
    synthesized bases from scratch, silently discarding any widening - the
    defect this spec closes, and it closes it by DELETING the operation rather
    than documenting a RESET.
  - The `setTreatment` virtual is re-signed to `setForestBasis(f, values, q)`
    through all four layers (`combiner.hpp:599`, `chain.hpp:863-865`,
    `sampler.hpp:1158-1160`, `facade.hpp:274`/`:521`) - the same four-layer
    edit shape over the same four files that collision (a) already prices for
    `bcfGlue`. Do both in one pass.
  - `dbarts_sampler_setForestBasis` (`C_interface.cpp:759-794`) and
    `$setForestBasis` (`R/dbarts.R:1214-1253`) route to `installForestBasis`
    directly; the 2-column collapse-to-z at `C_interface.cpp:790-791` is
    DELETED, not generalized (see the `glue_.z` clause), and
    `$setForestBasis`'s `index != 1L` refusal (`R/dbarts.R:1218-1223`)
    relaxes. 1-based `forest` at the R5 layer is untouched
    (`resolveForestIndex`, `R/bartcore.R:1051`).
  - **Ordering is LAST INSTALL WINS, per forest**, and both orderings collapse
    to that one rule because there is only one operation. widen-then-swap: the
    "swap" IS an install on the forest swapped, so forest 0's widening
    survives it. swap-then-widen: the widening replaces exactly that forest's
    basis, and the swapped forest's block is carried to its new offset. The
    two commute exactly - `rebuildAmplitudeLayout` (`:1151-1171`) derives
    offsets as a pure prefix sum of the width vector and carries by position at
    `min(newWidth, oldWidth)`.
  - **Amplitudes PRESERVE-and-remap**, which is what `rebuildAmplitudeLayout`
    already does: a width-preserving install early-returns at `:1159` and is
    the BITWISE IDENTITY on every amplitude (that is bcf's mid-life z swap, and
    it is baseline-gating); a width change carries every other forest's block
    bitwise to its new offset and enters new coordinates at the neutral 1.0.
  - **`glue_.z` RETIRES with `setTreatment`.** It has exactly ONE writer,
    `combiner.hpp:1121` inside `installTreatment`, and two readers, `:836` and
    `:869` inside `drawShippedGlue`; the pointer is the borrowed
    `BartcoreHolder::ownedTreatment` (`R_interface_bartcore_common.hpp:33-34`)
    whose contents only `R_interface_bartcore.cpp:3679-3684` and
    `C_interface.cpp:790-792` refresh. Delete the collapse-to-z while keeping
    `glue_.z` and it freezes at its construction value forever: a
    width-preserving z swap would install the new pair into `basis[1]`, leave
    the layout unmoved, and then draw `b0`/`b1` partitioned by the OLD z at
    `:869` while `forestMultiplier` (`:1081-1090`) contracts the NEW basis - a
    silently different model on the one route this spec calls bitwise, and one
    the amplitude assertions cannot see. So `drawShippedGlue` must partition
    from `basis[1]` column 1 (`values[2 * i + 1]`, row-major) and `glue_.z`
    goes. The stored values are the same values `ownedTreatment` holds, so
    this is expected to be BITWISE on the baselines - a GATED claim, not an
    assumption: the trio (`bcf-equivalence`, the multinomial baseline,
    `tests/cpp`) proves it before the slice lands. The bridge note at
    `R_interface_bartcore.cpp:2489-2493` (an M4.1 review correction) records
    the two objects as SEPARATE and must be rewritten with them merged. The
    same clause answers the general route: with no `glue_.z` there is no
    `nullptr` for a K-forest spec that carries per-forest bases and no z.
  - **Draw-path selection is a per-forest IS-CANONICAL predicate, not a width
    test.** `drawGlue` (`:732-738`) today selects `drawShippedGlue` on
    `forests.size() == 2 && !generalAmplitudeDraw_`, i.e. on the forest count.
    Selecting instead on widths (K == 2, q0 == 1, q1 == 2) is WRONG, because
    M4.3 opens legal CONTINUOUS two-column bases into exactly that shape -
    `consumer.c`'s `LEG_BASIS_VALUES` (`:1022-1032`) installs a 0.25/0.75 pair
    that today is pinned to raise - and `drawShippedGlue` is not a general
    q = 2 conditional: it never reads `basis[1]` at all, it forms two disjoint
    group-precision accumulators keyed on the indicator (`:865-872`). The
    predicate is therefore a VALUE property: forest f is CANONICAL when its
    basis is exactly the constructor's shape (all-ones; or, at the treatment
    forest, a complementary 0/1 pair - each entry in {0, 1}, each row summing
    to 1). A non-canonical basis at ANY forest forces the general path
    (`drawAmplitudes`) for that draw.
  - **The predicate is RECOMPUTED, never serialized**: one `bool` per forest, a
    pure function of the basis values, recomputed at install and at state
    restore - O(n) once per install, not once per sweep. **This corrects the
    scoping report's third decisive argument**, which claimed SUBSUME needs
    "the LEAST new state - none". It does not. It needs exactly this one
    recomputable per-forest bool: the minimal state of the three options, and
    strictly less than the persistent installed/synthesized flag REFUSE needs,
    but not zero. The argument stands as re-stated here, not as written there.
  - **`installForestBasis` becomes the terminus of both public entries and must
    GUARD.** Today it checks nothing - `glue_.basis[f]` indexes a vector sized
    only by the constructor, no `numColumns >= 1` check, no finiteness check -
    which is fine while two `tests/cpp` fixtures are its only callers. As the
    sole mutator it owns the dimension guard (`f` in range, `numColumns >= 1`),
    the finiteness check, and the canonical-predicate recompute. Same class of
    item as the `amplitudes[3]` overflow: correctness, not signature.
  - **Persistence is MANDATORY, and its ORDER is part of the contract.**
    Nothing about a basis survives serialization today (`serializeGlue`/
    `restoreGlue` `:803-816` carry `a, aVariance, b0, b1` only, through
    `ChainStateData` `:92-94`); the basis is re-synthesized at re-creation from
    `data@treatment`. The successor is M1's pattern - an R5 field mirroring
    every installed basis, replayed at the same THREE sites
    `reapplyForestWeights` is (`R/dbarts.R:1564` `getPointer()`, `:1590`
    `setState()`, `:936` `$copy()`) - but NOT with the same ORDERING. Those
    three sites all reapply AFTER the state install, which is right for
    WEIGHTS precisely because weights are deliberately not in the state
    (`R/bartcore.R:997-999`). Amplitudes ARE state. Reapply a widening after
    `setState` and the restore either refuses on length or writes a ragged
    vector into the narrow creation layout, after which the widening enters the
    added coordinates at the neutral 1.0 and every restored amplitude the
    widening covered is lost. **The contract is RESTORE-THEN-WIDEN**: widening
    after the restore must preserve-and-remap the RESTORED amplitudes. Either
    the bases ride creation (the way `data@treatment` does today) or the basis
    reapply runs BEFORE `setState` at all three sites - which makes it a
    different helper from `reapplyForestWeights`, not the same hook. Whichever
    the implementer picks, the arm asserts AMPLITUDE VALUES across the round
    trip, not just the surviving basis width.
  - **The persistence GUARANTEE survives the slot retirement.** Collision (e)
    retires `data@treatment` as `$setForestBasis`'s mirror;
    `man/dbartsSampler-class.Rd:167` documents, and
    `test-bcf-r5-surface.R:167` pins, that a re-created sampler carries the
    current assignment. The SLOT may go; the GUARANTEE may not - the successor
    mirror carries it. A tightening of collision (e), not a new public
    semantic.
  - **The not-taken alternative, recorded: REFUSE, never RESET.** If an
    implementer keeps `setTreatment` as a live mutator, the only honest rule is
    that it REFUSES on a combiner any of whose bases was installed rather than
    synthesized. That is a cheaper slice by the four-layer re-sign, but it
    costs a persistent per-forest installed flag and leaves two mutators whose
    orderings must BOTH be pinned. NOT TAKEN, on two grounds: (i) both shipped
    public basis entries route to `setTreatment` (`C_interface.cpp:792`,
    `R/dbarts.R:1239`), so under RESET the documented behavior of one public
    entry is "silently undoes the other" - a surface whose two halves fight;
    (ii) a width change now SELECTS THE CONDITIONAL (the canonical predicate
    above), so it is a MODEL fact, and a model fact cannot be revertible by a
    data swap - which is what `installForestBasis`' own Doxygen already asserts
    (`:604-606`). A documented RESET is the one option that must not ship.

  **2. The test shape that discharges the mandate** (both orderings,
  red/green-able). Engine, `tests/cpp/test_sampler.cpp`, replacing
  `testForestBasisSynthesis` (`:3296-3352`, today the ONLY component guard on
  either half) as `testForestBasisOrdering`:
  - **Arm 1, widen-then-swap.** Install a q0 = 2 prognostic basis with a
    DISCRIMINATING second column (nonzero, row-varying - unit values vacate
    pins, the M4.0 lesson), restore distinct amplitudes, then swap forest 1.
    Assert `basis[0]` is still two columns, the prognostic multiplier still
    reads its second term on every row, and forest 1's block is bitwise its
    swapped value. RED under today's reset.
  - **Arm 2, swap-then-widen.** Swap forest 1, then widen forest 0: the swapped
    values survive the layout move and forest 1's two coordinates are carried
    bitwise to their new offsets. RED if the install re-synthesizes.
  - **Arm 3, width-preserving canonical reinstall is the BITWISE IDENTITY** on
    every amplitude (the `:1159` early return). This is the pin that stands
    between M4.3 and a `bcf-equivalence` re-record.
  - **Arm 4, the z divergence** - the leg Arm 3 structurally cannot see, since
    Arm 3 asserts amplitudes and the defect is in the NEXT draw. A
    width-preserving swap to a DIFFERENT complementary pair must change that
    draw's partition: draw, and assert `b0`/`b1` reflect the new grouping. RED
    against a retained `glue_.z` (old-z-new-basis), GREEN once
    `drawShippedGlue` partitions from `basis[1]`.
  - **Arm 5, draw-path selection.** (i) q0 = 2 does not take the two-scalar
    path; (ii) a NON-CANONICAL width-2 basis at forest 1 (the 0.25/0.75 pair)
    does not either - by rng consumption count, or by asserting the drawn
    a-block's second coordinate moves. 5(i) is RED today; 5(ii) is the arm a
    width-only predicate passes wrongly.

  R5, `inst/tinytest/test-bcf-r5-surface.R` (or a new
  `test-forest-basis-r5.R`):
  - **Arm 6, both orderings at the surface.** `$setForestBasis(1, B)` then
    `$setForestBasis(2, C)` and the reverse: the resulting per-forest fits and
    amplitudes agree, and neither install disturbs the other forest.
  - **Arm 7, persistence.** Widen, `storeState`, `saveRDS`/`readRDS`, force
    `getPointer()`'s re-creation; assert the widened basis AND the restored
    amplitude VALUES survive. Run M1's per-site falsifier discipline: comment
    out each of the three reapply sites in the installed copy and confirm
    EXACTLY its own arm goes red.

  **3. The state `bcf` block must carry PER-FOREST WIDTHS.** Collision (b)
  re-encodes the length-4 block in place; a total length is not a layout -
  `q = (1, 3)` and `q = (2, 2)` both serialize four amplitudes, and both
  readers today check `Rf_xlength(bcfExpr) != 4` (`:6167`, `:6366`), whose
  natural generalization to `sum q_f` cannot tell those apart. `restoreGlue`
  (`combiner.hpp:810-816`) writes THROUGH `amplitudeOffset` (`:369-376`), so a
  mis-laid vector is silently mis-assigned rather than refused, and combined
  with the ordering hazard in item 1 a save/load across differently-shaped
  bases would silently permute forest blocks. NAME THE ENCODING: a K-length
  width vector followed by `sum q_f` amplitudes (or a reader that validates the
  incoming layout against the live `amplitudeOffset`). And there are THREE
  state sites, not the two collision (b) names: write
  `R_interface_bartcore.cpp:5814-5820`, read `:6165-6175`, and the WARM-START
  DONOR read `:6364-6374`, which carries its own length check and its own
  "malformed bcf glue in warm-start donor" message. The round-trip arm pins a
  REDISTRIBUTION at constant total, not only a malformed length, at both
  readers.

  **4. The route enumeration is THREE R routes plus the flat one.** Beyond
  `$setForestBasis` (`R/dbarts.R:1214-1253`, `.Call`ing
  `C_dbarts_bartcore_setTreatment` at `:1239`) and the flat
  `dbarts_sampler_setForestBasis` (`C_interface.cpp:759-794`), there is
  `dbarts:::bartcoreSetTreatment` (`R/bartcore.R:979-981`), which calls the
  bridge entry directly and is ABSENT from the R site list above - add it. Its
  callers: `benchmarks/R/bcf-equivalence.R:190`, `test-bcf.R:58,78,138`,
  `test-bcf-mutation-pins.R:74`, `test-multi-forest-seam.R:177,285`. The first
  is the mutation call inside equivalence scenario (e), which runs, swaps and
  RUNS AGAIN, so **the harness must keep working or be migrated IN THE SAME
  COMMIT, with the baseline UNTOUCHED** - this slice does not license a
  re-record. The same holds for the `bartcoreBCFGlue` adapter
  (`R/bartcore.R:1013`, read at `bcf-equivalence.R:97` inside
  `recordChannels`), whose 3 x n.chains SHAPE the bridge's `getBCFGlue`
  per-forest rewrite moves. So collision (f)'s thin-adapter clause extends BY
  NAME to `bartcoreBCFGlue`'s 3-row reading and `bartcoreSetTreatment`'s
  z-argument signature; either they are preserved as adapters, or they move in
  the landing commit with the recorded channels bitwise either way.

  **5. The column-major/row-major transpose is a named BRIDGE item, and a
  silent-wrong-answer one.** `inst/include/dbarts/dbarts.h:686-687` documents
  the flat basis as COLUMN-major numObservations x numColumns and
  `C_interface.cpp:781` reads it that way; `ForestBasis` is ROW-major by
  contract (`combiner.hpp:303-317`, "Stored ROW-major (row i at
  i * numColumns) because that contraction is the only read") and
  `installForestBasis` assigns straight through (`:611`). Nothing is red today
  only because the collapse-to-z at `:790-791` discards the layout entirely.
  **RESOLVE IN THE DOC DIRECTION**: the header text changes to say ROW-major,
  matching the engine; no signature moves and no hash re-bake. The flat entry
  then copies through, `consumer.c`'s basis-building legs re-lay ROW-major, and
  the R bridge transposes, since R matrices are column-major. M2's R5 route has
  no live mismatch - a two-level factor collapses to a single column
  (`R/model.R:707`) and one column is layout-invariant - but it INHERITS the
  fix the moment q > 1.

  **6. The re-price. TESTS ~645 dense-equivalent (band ~590-660)**, stated as
  M4.3's OWN budget: "the remainder is M4.3's rewrite share" is STRUCK, since
  the ~350-450 band is spent at ~797 (item 8). Density convention, stated
  because the first pass had none and was optimistic without it: a refusal
  assertion in `test-bcf-creation.R` measures 9.2 lines (multi-line `dbarts()`
  calls), a creation-shaped POSITIVE arm costs more than the refusal it
  replaces (a spec read, a run, then assertions) at ~14 lines/assertion, a
  simple `expect_*` on an existing object ~5, and a `tests/cpp` arm ~10.
  - T1 ~130, the creation-pin relaxations (6 refusals -> positive routes, 1
    restated, 4 slot-retirement rewrites, 4 params assertions). Priced ~75 at
    the first pass, which is BELOW the ~95 lines of live text it REPLACES:
    `:556-610` = 55, `:656-664` = 9, `:316`/`:317` + `:808-818` + `:819-829` =
    24, the 4 params assertions = 7. Re-priced at the convention above.
  - T2 ~95, the K > 2 / wide-basis positive surface: creation at K = 3-4, a
    continuous basis, a >2-level factor, run sanity, and the FA0 K = 2 bitwise
    reduction through the new spec route (re-checked at the same convention).
  - T3 ~75, the interaction spec, arms 1-7 (cpp ~45, tinytest ~30).
  - T4 ~55, `forestParams` ragged transport: per-forest resolution, both bridge
    routes agreeing, the retired length-8 message.
  - T5 ~55, state block round trip at q > 1 / K > 2: both readers (`:6165`,
    `:6364`), the redistribution-at-constant-total case, malformed length,
    keepTrees continuation.
  - T6 ~40, ragged run-result `glue` slot: dims plus the reconstruction
    identity at K = 3, three sites.
  - T7 ~40, `$getForestAmplitudes(forest)` signature break: per-forest read,
    ragged widths, out-of-range refusal, the `test-bcf-reporting.R` helper.
  - T8 ~55, treatment-slot retirement + capability-probe successor: the five
    `refuseBCFMutation` sites stay red (`R/bartcore.R:285`, `:306`, `:383`,
    `R/dbarts.R:1044`, `:1484`), the positive probe, the multinomial
    non-misfire.
  - T9 ~65, flat C: the four `consumer.c` legs (two inverting), a wide basis
    accepted, a ragged amplitude read, the `LEG_COUNT` bump, `test-capi.R`.
  - T10 ~35, `bcfGlue` four-layer re-sign coverage plus the `amplitudes[3]`
    overflow (`C_interface.cpp:767`, `:812`, `:824`).

  **7. The slice's non-test price: ~970-1110 dense-equivalent.**
  - **Engine ~250-300.** K-length forest spec vector + `BCFSpec` thin adapter
    (~60); the constructor's K `buildBCFForest` loop and calibration map
    (`chain.hpp:692-724`, `:4457`) (~40); construction-only synthesis plus
    `installForestBasis` as sole mutator WITH its guards, and the four-layer
    `setTreatment` -> `setForestBasis` re-sign (~55); the `glue_.z` retirement
    and the canonical predicate (~15); `bcfGlue` four-layer re-sign to a
    length-carrying read (~45); `serializeGlue`/`restoreGlue` +
    `ChainStateData` ragged with widths (~45); reporting-channel counts (~20).
  - **Bridge ~280-320.** Ragged params transport at both sites (~70); the state
    block re-encoded in place at all three sites (~55); ragged `glue` result
    slot (~35); `getBCFGlue` -> per-forest (~30); the basis entry -> per-forest
    install INCLUDING the transpose (~45); treatment-slot retirement + the
    two-way cross-check rewrite `:2833-2838` (~40); `amplitudes[3]` buffers
    (~10). Band ~350-450, nothing consumed.
  - **R ~360-410 + man ~35.** As scoped (`forestParams` ragged ~55;
    `forestBasisDeclaration` + `expandForestBasis` per-forest ~80;
    `resolveForests`' six positional refusals -> per-forest loop ~70;
    `resolveModerators` ~25; basis validation generalized off
    `validateTreatment` ~30; `$setForestBasis` + `$getForestAmplitudes(forest)`
    ~45; `isBCFSampler` successor + its five call sites ~30; the basis mirror +
    its reapply at three sites ~30; Rd deltas ~35), plus ~10 for the internal
    `bartcoreSetTreatment` route (item 4). Band ~400-500, nothing consumed.
  - **Flat C ~45.** Guard bodies at `C_interface.cpp:768-792`, `:815-816`,
    `:828-832`, and the three `amplitudes[3]` stack buffers -> K-sized.
  - `dbarts.h` MOVES NO SIGNATURE: `setForestBasis` already takes `numColumns`,
    `numForestAmplitudes`/`forestAmplitudes` already take `forest`, and
    `:698-699` says so in the header ("widening it relaxes this guard and moves
    no signature"). **No hash re-bake, no version-constant move.** Doxygen
    sweeps ride M4.5; the transpose sentence (item 5) is the one exception and
    is a doc-only edit.
  - **Engine-band arithmetic, recorded AS ARITHMETIC and not as a verdict.**
    The engine band is ~500-700 with ~228 consumed (M4.1 +48, M4.2 ~180);
    M4.3 at ~250-300 puts it at ~478-528, leaving ~172-222. **M4.4's engine
    work was UNPRICED here; it is PRICED in M4.4's own bullet** (49-82 dense,
    stop 123, against ~393 actually consumed through M4.5), and the band
    holds at the stop. There is no breach to record.

  **8. Erratum: the consumed tests band is ~797, not ~732** - M4.0 ~348 +
  M4.1 ~98 + M4.2 ~351 (M4.2's own landing note reads "tests ~351", i.e.
  286 + 65) = 797. Corrected in place in the STATUS paragraph and the Budget.

  **9. Collision status, re-verified at e7708b7c.** (a) OPEN in full, all four
  layers unchanged (`combiner.hpp:463`, `chain.hpp:1021-1023`,
  `sampler.hpp:1242-1244`, `facade.hpp:315`/`:546`) - do it in the same pass as
  item 1's re-sign. (b) OPEN, and it is THREE sites not two (item 3); its
  `amplitudeOffset` prerequisite is DONE (M4.2 handoff (b)). (c) OPEN, and it
  is THREE glue-stride sites not one: add `chain.hpp:4661-4662`
  (`results.glue + sampleNum * 3`) and `sampler.hpp:321` (`results.glue +
  c * numSamples * 3`) beside the bridge alloc `:4098-4102`. (d) OPEN; its
  twin, `installTreatment`'s `basis.resize(2)` K-sizing, is DONE - every
  per-forest array is now sized by the chain's forest count
  (`combiner.hpp:1122`, `:1152-1153`; M4.2 handoff (d), with its ASAN
  red/green). (e) OPEN, tightened by item 1's persistence clause. (f) OPEN but
  zero-cost IF the adapter clause is honored, and it is load-bearing on the C++
  suite's ground as well as `bcf-equivalence.R`'s: 25 `tests/cpp` fixtures
  construct `BCFSpec`, so the clause is not a courtesy to `bartcoreBCFSampler`
  alone. Item 4 extends it by name to two more functions.

  **10. Stale `combiner.hpp` anchors** (M4.2 moved the file again). Verified BY
  SYMBOL at e7708b7c; these correct the M4 text above wherever it cites them.
  `BCFSpec` `:244-251` -> `:244-263` (M4.2 added `ridgeA`/`ridgeB`/
  `generalAmplitudeDraw`); `BCFState`'s `a, b0, b1, aVariance` `:290-302` ->
  struct `:342-377` with accessors `:369-376`; the constant-leaf
  `static_assert` `:502-503` -> `:583-584`; the base `bcfGlue` virtual
  `:386-389` -> `:463`; `forestMultiplier` `:761-765` -> `:1081-1090`;
  `combinedFits` `:567-576` -> `:691-706`; `installTreatment` `:855` ->
  `:1119-1142` (the resize at `:1122`); `drawGlue` `:582-623` -> `:732-738`
  (plus `drawShippedGlue` `:824-877`); `afterCombine` `:638-709` -> `:758-767`;
  `MultinomialForestCombiner::afterCombine` `:1152-1197` -> `:1560-`;
  `buildBCFForest` `chain.hpp:4454` -> `:4457`; `Chain::bcfGlue`
  `chain.hpp:1018-1020` -> `:1021-1023`; the bridge's length-8 check #2
  `:2913-2914` -> `:2915-2916`; its cross-check `:2831-2836` -> `:2833-2838`;
  `amplitudes[3]` `:3673` -> `:3675` and the stride-3 alloc `:3769-3772` ->
  `:3771-3774`; the state write `:5812-5818` -> `:5814-5820` and read
  `:6163-6167` -> `:6165-6169`; `R/spec.R`'s control attribute `:466-486` ->
  `:467-483`. UNCHANGED **AT e7708b7c, WHICH IS ALL THIS LIST EVER CLAIMED** -
  M4.3 and M4.5 have landed since and moved several R-side entries (`R/data.R`'s
  `validateTreatment` is GONE outright), so read it as a stamp and not as
  current: `BCFForestSpec` `:211-241`, every `R/model.R` site, every other
  `R/spec.R` site, `R/data.R:643-655`, `R/bartcore.R:18-22`, `:30-34`, `:285`,
  `:306`, `:383`, `:622`, `:1051`, `R/dbarts.R:1044`, `:1484`, `:1442-1446`,
  all three `C_interface.cpp` `amplitudes[3]` buffers and its guard bodies,
  `facade.hpp:274`/`:315`/`:750`, `sampler.hpp:1242-1244`,
  `R_interface_bartcore.cpp:1112`, `:2134-2135`, `:2138-2148`, `:4098-4102`.
  The two entries M4.4 reaches are corrected here rather than left to it:
  `$setForestBasis` `R/dbarts.R:1218-1223` is now DEFINED at `:1223` with its
  body to `:1261`, and its `index != 1L` refusal is GONE - M4.3 removed it, and
  the docstring at `:1224` reads "at any forest and any width"; its `:1237`
  single-slot `data@treatment` mirror is now the `data@bases` mirror at `:1245`.
- **M4.4 (non-Gaussian, the family's justification).** Wire the family enum
  through the K-forest constructor and define the calibration map against each
  family's latent scale; probit and logistic in v1, the rest doors.

  **The family gate is SIX sites, not four** (the count was amended to four
  2026-08-13; ANCHORS RE-DERIVED BY SYMBOL 2026-08-14, when ten of the eleven
  this bullet used to cite came back stale). Every anchor below and throughout
  this bullet was located by opening the file at `efec6ba2`; none is an offset
  of an earlier number. The struck citations are kept in parentheses so a
  reader of the old text can find where it went.
  - `R/spec.R:423-431` - the `family != "gaussian"` refusal off `data@bases`:
    guard `:423`, predicate `:424`, `stop()` `:425-430` (its message literals
    `:426-429`). (Was `:412-419`, now a comment inside the BCF prose block.)
  - bridge `R_interface_bartcore.cpp:2299-2300`, the first test inside
    `refuseUnsupportedBCFComposition` (`:2295-2337`), called from exactly one
    place, `createHolder:2916`. (Was `:2222-2223`.)
  - bridge `R_interface_bartcore.cpp:2993-2994`, inside `createBCFHolder`
    (`:2981-3038`), reached from `bartcore_createBCF` (`:3526-3534`) and the
    internal R helper `bartcoreBCFSampler` (`R/bartcore.R:626-708`).
    (Was `:2925-2926`.)
  - engine `chain.hpp:705` (`family_ = ResponseFamily::gaussian`), with
    `GaussianResponse` hardcoded at `:702-704` and the EXECUTABLE calibration
    map at `:722-727` under its Doxygen `:709-721`. (Was `:704` and
    `:716-719`; the latter lands mid-comment.)
  - engine `facade.hpp:777-789`, `createBCFSampler` (Doxygen `:774-776`).
    **This bullet MISCHARACTERIZED the site.** Its single
    `SamplerFacade<ConstantGaussianLeaf>` instantiation is the LEAF
    instantiation, not a family one: the family is a RUNTIME enum dispatching
    to a virtual `ResponseModel`, and the same specialization already serves
    all six families on the single-forest path (`facade.hpp:640-652`). M4.4
    adds ZERO template instantiations, ZERO compile surface, no object-size
    growth and no virtual-table move - so the stale-object bus-error gotcha is
    not in play beyond the standing header rule. What is actually baked in
    here is that the function takes no family argument at all. (Was
    `:750-762`, now `createSamplerOverStore`.)
  - engine `sampler.hpp:157`, `family_(ResponseFamily::gaussian)` in the
    K-forest `Sampler` constructor. The SIXTH site, absent from every earlier
    list including the calibration respec's own. `Sampler::family()`
    (`sampler.hpp:1123`) is what every family-keyed bridge predicate reads
    through `shape.family`; left pinned, a probit K-forest sampler LIES to
    `refusePinnedSigmaChange` (`R_interface_bartcore.cpp:2746-2757`, family
    test `:2752-2753`), `refuseBinaryWeightChange` (`:1626-1639`) and the
    `drawsSigma` branch in `bartcore_setModel` (`:4736-4737`).

  `model.hpp:2577`'s `enum class ResponseFamily` is NOT a gate at all: `probit`
  and `logistic` are already members. (Was `model.hpp:2523`.)

  **RE-SCOPED 2026-08-13 on the probe verdicts (`probes-2026-08-13.md`): this
  slice moves to the END of the arc, as the IMMEDIATE follow-on to M4.5, and
  its justification is rewritten.** FA5's arm B AGREES with the reference on
  all 12 functionals (max |z| = 2.54, threshold 3.0) with the power
  precondition met, and it did so from SHIPPED knobs only (`resid.prior =
  fixed(1)`, `sigma =`, `node.prior = normal(k, scale)`, creation-time
  `weights =`, `$setResponse(..., updateScale = FALSE)`, `$run(0L, 1L)`,
  `$getCalibration()`). **The headline ground - "the couplings a caller cannot
  compose" - is FALSIFIED and is STRUCK from this slice.** What remains, and
  all this slice may now be argued on:
  - **The measured composition tax.** **5.43x at K = 8** (5.14x at K = 2, 5.36x
    at K = 4; the K = 2 figure RE-MEASURED at **5.11x** 2026-08-14, below),
    where the denominator is stated explicitly because the number it
    replaces was not: per-sweep wall time of the K-sampler R composition over a
    SINGLE batched engine sampler carrying the same total tree budget (K * 50)
    on the same n, measured back to back on the same box (FA1b). The tax is
    FLAT in K (1.06x growth K = 2 -> K = 8), so it is a level, not a scaling
    argument. Plus an O(nK) R-level latent draw per sweep on top, for the
    non-Gaussian case specifically.
  - **Correctness unreachable from the R composition:** the per-forest ASIS
    rescale writes LEAF values, and no setter exists or should; and the
    `0x1p-26` snap's exactness, which the R host can only mirror by hand
    (arm B did mirror it - the constant is `combiner.hpp:787`, its snap branch
    `:823-827`, split across two places and NOT at the `:534-537` this bullet
    used to cite - and that is the point: the host must know the constant).
  - **Composition ergonomics:** K stores, K states, K `.Call` round trips, K
    leaf calibrations, a serial outer loop, and a host that must own the
    latent draw against the combined fit correctly to get arm B's answer.

  That is a WEAKER justification than this plan asserted, and it is written
  down as such rather than argued back up. **Entry criteria, added
  2026-08-13:** (i) RE-MEASURE the composition tax at K = 2 under the SAME
  protocol against the then-current engine at respec time - the harness exists
  (`.claude/m4-basis-design/harness/fa1b-composition-scaling.R`), so this is
  minutes, and a tax quoted without its denominator is exactly the defect the
  erratum records; (ii) respec the calibration map against each family's latent
  scale, as already written above.

  ---

  **AMENDED 2026-08-14 on the pre-slice pricing pass.** Five findings under
  `.claude/m4-basis-design/` (`pre-m44-fa1b-remeasure.md`,
  `-engine-pricing.md`, `-calibration-respec.md`, `-anchor-simulation.md`,
  `-open-questions.md`, indexed by `pre-m44-2026-08-14.md`). That directory is
  GIT-IGNORED, so this bullet carries every conclusion it depends on rather
  than pointing at it. **Both entry criteria are DISCHARGED, the engine work is
  PRICED and FITS, and the slice is CLEARED TO PROCEED.**

  **1. Entry criteria, discharged.**
  - (i) The composition tax at K = 2 re-measured **5.11x** against `efec6ba2`,
    versus 5.14x recorded 2026-08-13 against `47c1fbe1` - a 0.62% drift, an
    order of magnitude inside the +/- 15% materiality band and smaller than
    either run's seed-to-seed spread. Numerator and denominator both moved
    2-3% TOGETHER (common-mode box state), which is why the ratio is what is
    quoted. All 80 load-independent K = 2 statistics came back BITWISE
    identical to the recorded run, confirming M4.1/M4.2/M4.3/M4.5 left the
    shared single-forest draw path untouched. One honest caveat that travels
    with the number: neither arm of that harness exercises the K-forest engine
    surface, so what criterion (i) establishes is that the BASELINE both arms
    ride is unchanged, not that the native K-forest family costs what the tax
    predicts. The native arm still owes its own measurement, and arm E is
    where it lands.
  - (ii) The calibration map is RESPEC'D against each family's latent scale.
    Its ALGEBRA survives verbatim as a construction-time constant under both
    families; only its ANCHOR moves. See item 2.

  **2. The calibration anchor: SETTLED, Option L.** `s` becomes a family-keyed
  constant:

      gaussian : s = scaledResponseSd()          // UNCHANGED, bitwise
      probit   : s = 1.0
      logistic : s = M_PI / std::sqrt(3.0)       // 1.8137993642342178

  the link's own error sd in each family. The rejected alternative, Option C,
  is CGM's forest-total plausible-index budget (probit `3.0`, logistic
  `pi*sqrt(3)` - the shipped single-forest `node.scale` at `R/spec.R:301` and
  `:308`), exactly `3 * L` in both families, so the fork was a single factor of
  three and one decision covered both. A third option, keeping the empirical
  `scaledResponseSd()`, is REJECTED outright (item 2c).

  **2a. The decisive argument is ANALYTIC and FAMILY-FREE.** Take the model
  shape in which nothing adapts - every forest carries a basis, so every
  amplitude block is fixed-variance. It is the only honest test of an anchor
  because nothing in the sampler can compensate for it, and it is REACHABLE:
  `dbartsData(x, y, bases = list(matrix(1, n, 1), cbind(1 - z, z)))` resolves
  BOTH forests to `(50, 0.25, 3, 1, 0.674, 0.5, 0, 1)`, MEASURED against a
  private-lib build. **It has NO FIXTURE ANYWHERE.** Repo-wide,
  `grep "bases = list(" inst benchmarks | grep -v "list(NULL"` returns
  nothing, and `inst/tinytest/test-bcf-creation.R:392-425` - which this bullet
  used to cite as the exercise - is the `bases = list(NULL, zBasis)`
  hasBasis-escape block, i.e. the HALF-CAUCHY shape, the opposite arm.
  Building the fixture is MANDATED work (checklist item 23).

  Two basis forests at the shipped defaults give
  `sd(eta) = sqrt(2 * 0.5) * (s / 0.674) = 1.4837 * s`. The comparator is
  dbarts's own shipped SINGLE-FOREST binary prior on the same index, and it
  must be quoted TWICE, because the shipped binary default is not the `k = 2`
  an earlier version of this bullet assumed. `resolveNodeHyperprior`
  (`R/model.R:392-400`) resolves an unsupplied `k` to `chi(1.5, 2.0)` whenever
  `binary` is true; only `monotone || !binary` takes the fixed 2.0. At the
  prior `ChiKHyperprior` (`model.hpp:2500-2524`) draws
  `k = sqrt(Gamma(nu/2, rate 0.5/scale^2)) = sqrt(Gamma(0.75, rate 0.125))`,
  and `E[k^-2]` does not exist for a shape below 1, so **the shipped
  single-forest binary index prior has NO SECOND MOMENT** - the same
  Cauchy-tail property 2b uses to refute the "(-3, 3)" framing for the
  K-forest. Robust sd 1.587 at the shipped default, 1.500 at a fixed k = 2;
  both 2e6-draw simulations (seed 20260814), re-run for this revision.

  Probit below; logistic is every column times `pi/sqrt(3)`, and every ratio
  is identical because `node.scale` is `3 * (latent sd)` in BOTH families
  (`R/spec.R:301`, `:308`).

  | arm | L | C | L vs shipped k | C vs shipped k | L vs k = 2 | C vs k = 2 |
  |---|---|---|---|---|---|---|
  | K = 2 all-basis, `sd(eta)` | 1.4837 | 4.4510 | **0.934x** | **2.80x** | 0.989x | 2.97x |
  | K = 2 shipped shape, robust sd | 2.254 | 6.763 | 1.42x | 4.26x | 1.50x | 4.51x |

  **"Structural identity" and "within 1.1%" are STRUCK, and must not be
  restored.** `(1/0.674)/(3/2) = 0.9891` is an identity against a k-FIXED
  configuration the binary surface does not default to; against the shipped
  default there is no identity at all, only a simulated 0.934x. **What carries
  the verdict is the RATIO L : C, which is exactly 3 under every comparator
  and in both families**, because Option C's `s` is 3x Option L's by
  construction. Read off the corrected numbers: L lands within 7% of a prior
  dbarts ALREADY SHIPS for binary BART, and C lands 2.8x wide of it. A
  2e6-draw prior-predictive simulation (seed 20260814) against the shipped
  K = 2 shape, the all-basis shape and K in {2, 4, 8} points the same way on
  every arm and the other way on none. None of the struck arithmetic may be
  copied into `latentScaleAnchor`'s Doxygen, `multiplier-combiner.md`'s map
  section or `nameable-calibration.md`: it would ship a false statement about
  dbarts's own default in the one comment a future reader will trust.

  **2b. The "lands inside CGM's (-3, 3)" framing is REFUTED. Do not restore
  it.** The respec argued that Option L "lands on CGM's own (-3, 3) probit
  target once the amplitude is counted." It does not. On the shipped K = 2
  shape the index has CAUCHY tails (a half-Cauchy amplitude times a normal
  forest total), so `P(|eta| > 3) = 0.307` under L and 38% of the induced prior
  on p sits outside (0.01, 0.99); eta's 2.5/97.5 quantiles are -20.4 and 20.2,
  not -4.4 and 4.4. The plug-in numbers the respec computed are RIGHT (they
  evaluate the amplitude at its prior median) and the conclusion drawn from
  them is WRONG. **The case for L is COMPARATIVE, not absolute**, and it is
  the two-comparator table in 2a. Written that way, the case is stronger than
  the refuted version, not weaker.

  **2c. What a naive M4.4 would ship, quantified.** If the family enum is wired
  and `chain.hpp:722` is left alone, `s` is the sd of the COLD-START working
  response - `2y - 1 - offset` for probit (`model.hpp:3057-3061`),
  `4(y - 0.5) - offset` for logistic (`:3657-3664`) - i.e.
  `2 sqrt(p(1-p) n/(n-1))` for probit and `4 sqrt(p(1-p) n/(n-1))` for
  logistic (2.001 at p = 0.5, against Option L's 1.8138; the two families do
  NOT share one formula), a base-rate artifact. The prior narrows 5.03x
  between p = 0.5 and p = 0.01 (ratio `sqrt(0.25/0.0099)`, independent of n),
  in the direction that WITHDRAWS prior support from the region the data
  demands: at p = 0.01 the required index is `qnorm(0.01) = -2.326` and the
  prior's IQR is [-0.303, 0.302], with prior mass below the required index
  falling from 0.189 (p = 0.5) to 0.043 (p = 0.01). The failure signature is
  the worst possible: probit's naive `s` equals Option L's `1.0` at p = 0.5 up
  to the Bessel factor, so the defect is CORRECT on balanced data and will pass
  every balanced fixture, with no error and no warning.

  **2d. Two facts about the shipped prior that must travel with the anchor,
  and are new.** (a) **Adaptivity is capped at one forest, for any K.**
  `resolveForests` (`R/model.R:772-840`, refusal `:799-806`) requires every
  forest past the first to carry a basis, and `forestParams` writes the LITERAL
  `0` for `amplitudePriorScale` whenever a basis is present (`R/model.R:874`),
  from which the bridge derives `forest.ridge = false`
  (`R_interface_bartcore.cpp:2207`). So forests 2..K are ALWAYS fixed-variance,
  a half-Cauchy on a basis forest is unreachable through any shipped surface,
  and every forest added past K = 2 lands entirely on the non-adaptive channel.
  This makes the anchor MORE load-bearing at larger K, not less. (b) There is
  **no per-K renormalization anywhere in the map**: `sd.control`/`sd.moderate`
  are counts of one anchor unit PER FOREST and the K-forest total is their
  root-sum-square, so the index disperses as `sqrt(K)` by construction, with
  the exponent EXACTLY 1/2 (`ForestSpec` carries the factor and divisor per
  forest and `expandForestSpecs` never reads K; simulated log-log slope
  0.50025 over K in {1,2,4,8,16}). On the all-basis shape at unit row norms
  `sd(eta) = 1.04912 * sqrt(K) * s`, i.e. `1.4837 s` at K = 2 and `2.9674 s`
  at K = 8 - the latter 1.978x the fixed-k comparator, two thirds of the
  ratio that condemns Option C at K = 2. On the shipped K = 2 shape robust sd
  is 2.26 at K = 2, 2.96 at K = 4, 3.95 at K = 8. **`2.9674` is Option L's
  K = 8 index SCALE and, separately, Option C's K = 2 RATIO; they are not the
  same object and neither is "exactly Option C's K = 2 value" (which is
  `4.4510 s`).** Both facts are design statements owed a sentence beside the
  anchor table; the `sqrt(K)` law itself is DEFERRED, see 2f.

  **2e. The anchor is the slice's single biggest risk, and checklist item 25 is
  the gate on it.** The lines are cheap; the NUMBER in them is argued, not
  derived from the repo. Getting it wrong produces a sampler that RUNS and MIXES
  and is silently mis-calibrated, and no STANDING gate sees that:
  `bcf-equivalence` compares gaussian only, the C++ component tests pin
  mechanism not calibration, and arm E is a mixing-and-coverage comparison whose
  power against a mis-scaled prior is unmeasured. So the anchor gets its OWN
  gate rather than a concession - item 25 asserts the induced index prior
  against its analytic value and is RED on the naive cold-start anchor, on
  Option C and on a dropped basis row-norm divisor. The anchor derivation still
  gets a reviewer BEFORE the switch is written.

  **2f. The anchor's UNIT was wrong, and M4.4 SHIPS the fix (options A + E).**
  The anchor is stated PER FOREST and consumed PER UNIT OF BASIS ROW NORM, and
  nothing made the two agree. From the code (`forestMultiplier`
  `combiner.hpp:1304-1313`, `buildBCFForest` `chain.hpp:4495-4511`, the map
  `chain.hpp:722-727`), on the all-fixed-variance shape:

      sd(eta_i) = s * sqrt( sum_f (factor_f/divisor_f)^2 v_f ||B_f(i,.)||^2 )

  A basis-free forest's term is Cauchy and has no sd, so the formula prices
  the fixed-variance forests and the rest is Cauchy. `||B_f(i,.)||` is
  UNCONSTRAINED on every route: `validateForestBases` (`R/data.R:635-666`)
  checks list-ness, numeric-or-logical, row count, subset alignment and
  finiteness; `expandForestBasis` (`R/model.R:687-720`) passes a numeric
  matrix through; `dbarts_sampler_setForestBasis` (`C_interface.cpp:761-781`)
  scans finiteness only. MEASURED against a private-lib build:
  `bases = list(NULL, matrix(rnorm(n, 50, 10)))` resolves to
  `(50, 0.25, 3, 1, 0.674, 0.5, 0, 1)` and RUNS to completion; under Option L
  probit that is a prior sd on `eta` of 52.46 latent sd at `||b|| = 50`, with
  `P(0.01 < p < 0.99) = 0.119` for that term alone and `0.099` for the full
  K = 2 shape, against 0.879 for the shipped single-forest comparator. (The
  product-normal values. `eta` is a product of two normals, not a normal, so
  the normal approximation's 0.035 understates central mass by 3x, in the
  direction that makes the trap look worse than it is. Use the simulated
  numbers.) Under gaussian the posterior absorbs it - the same construction's
  basis-forest amplitude falls by ~41x against its unit-norm twin - which is
  exactly why the defect is invisible today and goes live when sigma is
  pinned.

  **VD decision 2026-08-14: SHIP THE FIX.** The map becomes

      nodeScale_f = factor_f * s / (divisor_f * nbar_f)

  with `nbar_f` the MEDIAN OF THE NONZERO row norms of forest f's basis, and
  `1.0` by convention when `spec.forests[f].basis == nullptr`. **The median,
  not an RMS, and the choice is named here rather than left to the
  implementer:** the map's own convention is already a median (`divisor =
  0.674` is `qnorm(0.75)`, `expandForestSpecs` `combiner.hpp:334`), the median
  is sparsity-invariant where an all-rows RMS inflates a bare 0/1 indicator
  column's scale by `1/sqrt(p)` - a legitimate treatment-only forest - and it
  is outlier-robust where a max is not. Statistically this is identical to
  normalizing the basis at ingestion, but the constant lives in the node scale
  rather than in the stored basis, so `$getForestAmplitudes` keeps reporting
  in the USER'S OWN basis units, no getter contract moves, no new state slot
  appears and `structSize` does not move.

  **`nbar_f == 1.0` EXACTLY on every shipped route**, which is why this should
  be free: the internal `bartcoreBCFSampler` route leaves the pointer null
  (`applyForestBases`, `R_interface_bartcore.cpp:2246-2266`, installs nothing
  when `data.bases[f]` is NULL, and the combiner synthesizes), and on the
  public `forests =` / `bases =` routes the basis is the `(1-z, z)` pair, a
  one-hot factor expansion or an all-ones column - every row norm exactly
  `1.0`, whose median is exactly `1.0`, and `divisor * 1.0 == divisor` in
  IEEE-754.

  **And the fix is CLOSED under `$setForestBasis`.** `nbar_f` is otherwise a
  construction-time constant and no mutation re-derives a K-forest leaf scale
  (`setModel`, `chain.hpp:1440-1507`, touches `forests_[0]` only), so a
  mid-sample basis swap would leave the divisor stale. That swap is a SHIPPED
  knob on three routes (`$setForestBasis` `R/dbarts.R:1223`,
  `Chain::setForestBasis` `chain.hpp:875-878`, flat C
  `dbarts_sampler_setForestBasis` `C_interface.cpp:761`), so a divisor that
  goes stale there would introduce a defect while fixing one. Checklist items
  21 and 22.

  **A DECISION POINT the implementer hits, not an open question.** E is GATED
  on the equivalence trio coming back BITWISE under `--preclean` on a build
  carrying E and nothing else. If any leg moves, **FALL BACK to option A
  alone**: ship L unchanged, write the identical contract sentences, record
  `forest(sd = )`'s per-unit-of-row-norm reading as a KNOWN GAP between doc
  and code in the landing note, and hand E to `binary-kforest-prior-default`.
  The fallback is still a strict improvement, because `man/forest.Rd:52`
  promises the unqualified version today.

  **The contract sentence, written in the form the fix makes TRUE.** Four
  sentences, all load-bearing, identical under A and under E, which is what
  makes the code fix separable from the contract. They go into
  `latentScaleAnchor`'s Doxygen, `combiner.hpp:258-273`,
  `multiplier-combiner.md`'s map section, `man/forest.Rd:52` and `:71`, and
  `nameable-calibration.md` alike.
  1. `forest(sd = )` and `amplitude.prior.variance` are stated PER UNIT OF
     BASIS ROW NORM: a forest whose basis rows have median nonzero norm `c`
     contributes the scale the argument names, and the map divides `c` out.
  2. The induced prior sd on the index at row `i` is
     `sqrt( sum_f prior.scale_f^2 v_f ||B_f(i,.)||^2 )`, with `prior.scale_f`
     read from `$getCalibration(f)[, "prior.scale"]`, `v_f` the forest's
     `amplitude.prior.variance` (default `0.5`) and `B_f` from
     `data@bases[[f]]`; valid on the fixed-variance forests, and a basis-free
     forest's contribution is Cauchy with no sd.
  3. Under probit and logistic that index is in LATENT sd units and sigma is
     PINNED, so nothing in the sampler absorbs a mis-scaled basis; under
     gaussian it is in `sd(y)` units and a drawn sigma partly does.
  4. The index disperses as `sqrt(K)`: `1.04912 sqrt(K) s`, so `1.4837 s` at
     K = 2 and `2.9674 s` at K = 8, at unit row norms and shipped defaults.

  **Two things the fix does NOT touch. Both are DEFERRED BY NAME to
  `binary-kforest-prior-default`, and M4.4 must NOT do either.**
  (i) **`sqrt(K)`.** No basis-side option changes the exponent - A, C and E
  all leave `1.04912 sqrt(K) s` exactly as it is - and the only option that
  does (normalizing the INDEX rather than the forest) is rejected on its own
  grounds: it shrinks every declared moderator's prior by `1/sqrt(K/2)`,
  which is the wrong direction for a varying-coefficient model, makes
  `forest(sd = 1)` non-composable, has no finite variance to sum on the
  half-Cauchy channel, and forks the gaussian and latent maps. Any K fix is
  therefore a change of DEFAULTS, which no slice's surface pins. **M4.4 owes
  the exponent and the two endpoints in PROSE and owes no code.**
  (ii) **Observability.** `v_f` is readable from NO getter
  (`$getForestAmplitudes`, `R/dbarts.R:1450`, returns draws;
  `$getCalibration`, `:1468`, returns `prior.scale` and no amplitude prior),
  so a user who took the default cannot reconstruct their induced prior from
  the public surface at all. Closing that (the costing study's option F: a
  reader reporting `prior.scale_f`, `v_f`, a row-norm summary and the implied
  `sd(eta_i)` quantiles, R 8-14 dense) is the cheapest thing that makes any
  anchor checkable by the person who gets it wrong - and it is NOT M4.4's.

  **One arithmetic correction, recorded so it is not re-derived wrong.**
  `forest(amplitude.prior.variance = 2)` at K = 2 all-basis gives
  `sqrt(2*2)/0.674 = 2.9674 s`, which is the K = 8 index scale, NOT Option
  C's. The value that reproduces Option C from Option L's anchor is
  **`v = 4.5`**, since `sqrt(2v)/0.674 = 3/0.674` needs `2v = 9`. `v` is a
  documented knob (`man/forest.Rd:71`) and a free multiplier on the index
  variance, so any anchor is reachable from any other through it. That is not
  a defect and not a decision - one sentence in the docs.

  **3. Price, by layer. All figures DENSE-EQUIVALENT non-comment lines** (the
  repo is air/clang-format wrapped at 80 columns, so raw diff counts inflate
  roughly 2x; comment lines are quoted separately and never folded in). Every
  site was opened and read.

  | layer | dense-equivalent | band | consumed before M4.4 |
  |---|---|---|---|
  | **Engine** (`src/bartcore/*.hpp`) | **49-82** | ~500-700 | ~393 (M4.1 ~48 + M4.2 ~180 + M4.3 ~165 + M4.5 0) |
  | **Bridge** (`src/R_interface_bartcore*.{cpp,hpp}`) | **14-28** | ~350-450 | 0 |
  | **R** (`R/*.R`) | **15-32** | ~400-500 | 0 |
  | **Flat C** (`src/C_interface.cpp`) | **0** | ~45 booked at M4.3 | - |

  Engine, site by site. FAMILY WORK, 29-55: `chain.hpp:705-706` family + sigma
  pin 3; `chain.hpp:702-704` the 3-arm response switch 10-14;
  `chain.hpp:722` the family-keyed anchor helper 8-12;
  `combiner.hpp:294-317` `BCFSpec` gains the family member 1-2;
  `sampler.hpp:157` 1; `facade.hpp:777-789` door refusal 2-8;
  `chain.hpp:1408` `combinedFits()` 2-3; `chain.hpp:887` the conjunct 2-4;
  `chain.hpp:1514` the residuals fix 2-4 (a guarded substitution, item 6, not
  a bare one); `model.hpp:2577` 0. SCALE FIX (2f), 20-27:
  `chain.hpp:4443` the `basisRowNorm` helper beside `scaledResponseSd` 14-17
  (median form; an all-rows RMS is 11-13 and is REJECTED in 2f);
  `chain.hpp:724-727` the changed map expression 0 (same line);
  `chain.hpp:875-878` the `$setForestBasis` recomputation plus the retained
  per-forest factor/divisor/`s` members 6-10. Adjacent and NOT folded in:
  engine ~50-70 COMMENT lines, bridge ~10-15, `dbarts.h` ~3-8.

  Flat C is **0 non-comment**, MEASURED four ways: `dbarts_sampler_create`
  (`inst/include/dbarts/dbarts.h:514-515`) ALREADY takes the family string and
  its Doxygen at `:495-513` already says the K-forest family is created through
  this same entry point; the four K-forest guard entries carry no family
  parameter and no family guard (`dbarts.h:724`, `:731`, `:739`, `:750`; TWO
  of them gate on `totalAmplitudes()` - `setForestBasis`
  `C_interface.cpp:768` and `forestAmplitudes` `:812-813` - while `forestFits`
  `:788` and `numForestAmplitudes` `:798-799` gate only on
  `forest >= shape().numForests`, which does not move the conclusion); no
  flat-C struct gains a member, so `structSize` does NOT move
  (`C_interface.cpp:266`, `:284`, `:308`); and `DBARTS_C_API_STRINGIZE`
  stringizes return type, name and parameter list only, with the preprocessor
  stripping comments before `#` applies. **No `dbarts.h` signature moves, no
  hash re-bake** - `DBARTS_C_API_HASH` (`dbarts.h:101`, `0xcd88efcd67de55d7`),
  `DBARTS_C_API_MAJOR` (`:80`) and `DBARTS_C_API_MINOR` (`:81`) all stay put,
  and no stan4bart floor bumps. The M4.3 `structSize` tripwire
  `static_assert(sizeof(ChainStateData) == 344)` at `tests/cpp/common.cpp:59`
  does NOT move: latents are already a generic chain-level state slot
  (`chain.hpp:2647-2651`, restored `:3240-3241`, validity-checked
  `:2813-2814`, written
  `R_interface_bartcore.cpp:5879-5886` with no family test and no forest-count
  test) and the "bcf" state block is purely geometric.

  **Why it is this cheap, stated so the number can be attacked.** The engine is
  already family-polymorphic (a virtual `ResponseModel` selected by a runtime
  enum, SIX concrete families - `enum class ResponseFamily` has six members
  and `facade.hpp:640-651` six arms); the sweep already hands the response the
  COMBINED location, not forest 0's (`chain.hpp:1303-1304`, `:1313-1314`); the
  combiner's reparameterization `(r/m, w m^2)` is scale-free and null-weight
  safe at all four of its `w` reads (`combiner.hpp:832`, `:1057`, `:1087`,
  `:1171`); the ASIS ridge is likelihood-invariant and reads no sigma and no
  `w`; and `MultinomialForestCombiner` already SHIPS a K-forest chain running
  under a fixed-scale non-Gaussian likelihood (`chain.hpp:742-763`). What M4.4
  actually buys is a 3-arm switch, a family-keyed anchor, four one-line pins,
  and two refusal relaxations.

  **The multinomial precedent covers LESS than that last clause implies, and
  the gap is exactly what arm E must have power against.**
  `ForestCombiner::totalAmplitudes()` returns 0 in the base
  (`combiner.hpp:558`) and is overridden ONLY by `BCFForestCombiner`
  (`:772-774`), so the multinomial path has no amplitude, no basis, no
  `forestMultiplier`, no `(r/m, w m^2)` reparameterization and no `2^-26`
  snap. What it establishes is that a K-forest chain runs under a pinned
  sigma; what it does NOT establish is the interaction M4.4 introduces - a
  per-observation multiplier `m_f(i)` dividing a working response whose
  precision is a per-sweep PG draw, with a snap - which has NO shipped
  precedent anywhere in the engine. The price conclusion still holds; the
  reason is narrower than "already ships".

  **4. The work-item checklist.** Every anchor re-derived by symbol at
  `efec6ba2`. **F<n>** marks a fork; all seven are RESOLVED in item 5.

  ENGINE
  1. `chain.hpp:692-695` - the K-forest `Chain` constructor gains the family.
     It rides `BCFSpec` as a defaulted member, NOT a positional parameter
     (checklist item 9).
  2. `chain.hpp:702-704` - 3-arm response switch (gaussian / probit /
     logistic), mirroring the single-forest 6-arm switch at `:537-577`;
     `chain.hpp:705` `family_ = spec.family`. **F1: inline vs factored.**
  3. `chain.hpp:706` - adopt the single-forest `sigmaIsFixed_` rule, which is
     at `:606-608`: `(family != gaussian && family != aft) ||
     options.sigmaIsFixed`. Left alone, a probit K-forest enters `run` with
     `sigmaIsFixed_ == false` and calls `drawSigma` every sweep; the VALUE is
     harmless (`ProbitResponse::drawSigma` returns its argument,
     `model.hpp:3099-3101`) but a semantic gate set wrong is a defect whether
     or not it currently bites.
  4. `chain.hpp:722` - `double s = latentScaleAnchor(family);`, a new private
     helper, three arms. **Two properties the reviewer must check.** The
     gaussian arm calls `scaledResponseSd()` (`chain.hpp:4443-4455`) through
     the IDENTICAL expression - a refactor that computes `s` once and
     multiplies by `1.0` on the gaussian arm is NOT acceptable; `x * 1.0` is
     bitwise but the reassociation risk in `nodeScaleFactor * s /
     nodeScaleDivisor` is not worth the elegance, and `chain.hpp:718-720`
     already records that these expressions being "written exactly as bcf's
     were" is what keeps K = 2 bitwise. And the non-gaussian arms must NOT
     call `scaledResponseSd()` at all, not even to discard the result: it is
     O(n), and its presence is what a later reader would mistake for the
     anchor. Its Doxygen states the ANCHOR, and the logistic arm says "the
     logistic law's standard deviation" - NOT "the link's own error sd",
     which reads as a property of the augmentation and is false of one: the
     engine implements Polya-Gamma (`model.hpp:3501-3508`), whose working
     residual sd at `psi = 0` is 2 (`omega = 1/4`, `coldStart` `:3657-3663`),
     not 1.8138. `pi/sqrt(3)` is still the right anchor - it is the sd of the
     logistic law the log-odds index is stated against, and it is what makes
     `R/spec.R:308`'s `3 * pi/sqrt(3)` equal `3 * (latent sd)`.
  5. `chain.hpp:1408` and `:1514` - `combinedFits()` in place of
     `forests_[0].totalFits.data()`. **F2 (m-i / m-ii).** `:1514` is the
     FOLDED-IN residuals defect and is NOT a bare substitution: see SECTION 6
     for what it does per response family.
  6. `chain.hpp:887` - drop the `family_ == ResponseFamily::gaussian` conjunct
     from `supportsResponseMutation` (`:885-888`), and rewrite its Doxygen at
     `:879-884` (the block opens at `:879`, "Whether this chain's forest
     coupling permits"), which names item 5's site as the reason the conjunct
     exists. Paired with 5's (m-ii). **The governing rationale is NOT there,
     it is at `combiner.hpp:953-955`** - the two conditions
     `BCFForestCombiner::supportsResponseMutation() == true` (`:963`) rests
     on, the first of which is "the gaussian response re-maps y through the
     pinned (min_, range_) and touches no forest". Rewriting the chain-side
     Doxygen and leaving the combiner-side conditions saying "gaussian"
     leaves the argument false where it actually lives. It carries, but it
     must be RE-ARGUED there: the second condition ("this combiner caches
     nothing per-forest across sweeps - `formForestResponse` re-derives every
     per-forest residual and precision from y and w each sweep") needs
     restating for a family whose y and w are the latent working values
     refreshed at `chain.hpp:1304`.
  7. `facade.hpp:777-789` - `createBCFSampler` reads `spec.family`; add a door
     refusal (`return nullptr` for aft / ordinal / nbinom) beside the existing
     `numVarianceTrees` one at `:785`. Signature UNTOUCHED.
  8. `sampler.hpp:157` - `family_(ResponseFamily::gaussian)` becomes
     `family_(spec.family)`.
  9. `combiner.hpp:294-317` - `BCFSpec` gains
     `ResponseFamily family = ResponseFamily::gaussian;` (`model.hpp` is
     already included at `combiner.hpp:17`). Defaulted, so it touches ZERO of
     the ~37 existing `tests/cpp` `BCFSpec` fixture constructions - the M4.3
     adapter clause honored.

  BRIDGE
  10. `R_interface_bartcore.cpp:2299-2300` - the first family gate becomes an
      allowlist over {gaussian, probit, logistic}, refusing the doors by name
      with a per-family reason.
  11. `R_interface_bartcore.cpp:2993-2994` - the second family gate, same
      allowlist. **F3: in place vs delegate.**
  12. `R_interface_bartcore.cpp:2990-2991` - `createBCFHolder` hardcodes the
      family name to `""`, so `resolveFamily` (`:1558-1596`, the binary arm at
      `:1563-1564`) maps it to `probit` on a binary response and the route can
      NEVER reach logistic. **F4.**
  13. `R_interface_bartcore.cpp:2311` - key the node-scale gate to the family
      default rather than the literal `0.5`. **F5**, shared with item 16.
      `refuseUnsupportedBCFComposition` already takes `family` as its first
      parameter (`:2295`), so the family is in hand.
  14. **NON-ITEM, do not add it.** `enforceBinaryWeightPolicy`
      (`R_interface_bartcore.cpp:1604-1618`) on the two K-forest creation
      paths was listed as required work by the respec and is NOT: it is
      already called from `parseSamplerSpecification` at `:2518`, and BOTH
      `createHolder` (`:2851-2853`) and `createBCFHolder` (`:2990-2991`) call
      `parseSamplerSpecification`. Probit's outright weight refusal and
      logistic's positive-integer-count check therefore fire the moment
      `family` becomes non-gaussian on either route, with zero new code, and
      `refuseBinaryWeightChange` (`:1626-1639`) already covers the
      post-creation side at `bartcore_setData:4467` and
      `bartcore_setWeights:4659`. **A K-forest probit/logistic sampler
      inherits the whole weight policy for free.**

  R
  15. `R/spec.R:423-431` - the refusal becomes an allowlist with per-family
      reason text. Both of today's clauses are wrong under M4.4: the admitted
      set widens, and "has its own fixed error scale" is now the REASON probit
      and logistic work rather than a reason they cannot. `%not_in%` already
      exists (`R/dbarts.R:416`).
  16. `R/spec.R:452` - `"a non-default 'node.scale'" = model@node.scale != 0.5`
      fires on the plainest possible binary K-forest call, because the family
      switch at `R/spec.R:297-309` gives probit `3.0` (`:301`) and logistic
      `pi * sqrt(3.0)` (`:308`). **F5.**
  17. `R/model.R:394` - `k <- if (monotone || !binary) 2.0 else chi(1.5, 2.0)`
      gains a multi-forest disjunct, so that `R/spec.R:446`
      (`is(priors$node.hyperprior, "dbartsChiHyperprior")`) stops firing on a
      hyperprior the user never asked for. **NOT a fork - it has a clean
      shipped precedent in the SAME expression**, the monotone conjunct,
      documented at `R/model.R:386-391`: a structural feature makes the chi-k
      update meaningless, the fix is applied to the DEFAULT RESOLUTION rather
      than to the downstream guard (an unsupplied `k` resolves to the fixed
      2.0, so the guard sees a default and stays silent), an EXPLICIT
      non-default still refuses by name (`R/model.R:395-400`), and the forced
      value is not a modelling statement because it is discarded downstream.
      Point four is why the precedent transfers EXACTLY here: `buildBCFForest`
      pins `forest.k = 1.0` (`chain.hpp:4506`) and
      `forest.updateK = false` (`:4504`), so forcing a binary K-forest's
      default `k` to 2.0 changes no fitted model. The caller has what it needs
      at `R/spec.R:206-216` (`declaredBases`, `forestSpec`), resolved BEFORE
      `parsePriors` runs at `R/spec.R:257`, and
      `resolveNodeHyperprior` is called with `control@binary` at
      `R/model.R:175-179`. **It is 4-6 dense, not 1, and it needs a new
      formal:** the multi-forest fact is NOT in `parsePriors`' scope
      (`R/model.R:92-101`, formals `control, data, tree.prior, node.prior,
      resid.prior, resid.dist, monotone = NULL, parentEnv`), which
      `R/spec.R:257` reaches through `eval(parsePriorsCall)`. So: a new
      DEFAULTED formal on `parsePriors`, a `parsePriorsCall$<flag> <- ...`
      assignment in `R/spec.R` on the `monotone` pattern at `:244`, a new
      argument at `R/model.R:175-179`, and the disjunct.
      `resolveNodeHyperprior` has exactly two call sites; the second is
      `R/xbart.R:251` (`resolveNodeHyperprior(node.prior@k, binary = FALSE)`),
      which a DEFAULTED parameter leaves alone and a REQUIRED one breaks at
      load. Defaulted is silently correct - xbart fits no K-forest - and that
      is stated here rather than assumed.
  18. **ALL FIVE `refuseBCFMutation` messages** (`R/bartcore.R:34-38`), not
      just the one this bullet used to name. Four carry K = 2 "both forests"
      language on a K-length surface: `R/bartcore.R:291-292`
      (`setResponse(updateScale = TRUE)`) and `:312-313`
      (`setOffset(updateScale = TRUE)`), both "both forests keep leaf
      calibrations stated against the response transform fixed at creation";
      `:389-390` (`setData`), "both forests are calibrated against the data at
      creation"; `R/dbarts.R:1055-1056` (`setModel`), "both forests' node and
      tree priors are calibrated at creation"; and `R/dbarts.R:1496-1501`
      (`setCalibration`), "both forests' leaf scales come from the two-forest
      calibration map". The first two are the `setResponse` / `setOffset`
      conduit item 5 / F2 is opening, which is a STRONGER claim to M4.4's
      ownership than the calibration one this bullet argued from.

  ADDITIONAL
  19. Per-forest weight channel semantics under a latent family. **F6.**
  20. `combiner.hpp:800-802`'s "half the mantissa" claim. **F7.**

  ENGINE, the scale fix (2f). GATED: build these ALONE first and run the
  equivalence trio under `--preclean` before the family work goes on top. If
  any leg moves, DROP 21 and 22 and take the option-A fallback 2f names.
  21. `chain.hpp:4443` - `basisRowNorm(const ForestSpec&, std::size_t n)`, a
      new private helper beside `scaledResponseSd`, returning the MEDIAN of
      the nonzero row norms of `spec.basis` (row-major, `spec.numBasisColumns`
      wide) and `1.0` when `spec.basis == nullptr`. `ForestSpec` already
      carries `basis` and `numBasisColumns` (`combiner.hpp:285-286`), and the
      ctor builds the forests at `chain.hpp:724-727` BEFORE the combiner at
      `:729`, so the borrowed pointer is in hand at exactly the right moment.
      The map expression at `:724-727` becomes `forestSpec.nodeScaleFactor *
      s / (forestSpec.nodeScaleDivisor * basisRowNorm(forestSpec, n))`.
  22. `chain.hpp:875-878` - `Chain::setForestBasis` re-derives forest f's leaf
      scale from the NEW basis on a successful install, so the divisor cannot
      go stale on the one mutation route the multi-forest surface exists for.
      This needs the map's per-forest `nodeScaleFactor` / `nodeScaleDivisor`
      and the construction-time `s` RETAINED AS MEMBERS - they are ctor locals
      today. **Retain `s` rather than recomputing `scaledResponseSd()`:** the
      K-forest map is pinned at creation (`setResponse(updateScale = TRUE)` is
      refused, item 18), so recomputing would recalibrate onto a scale the
      other forests are not on, and reusing the stored double is what keeps a
      norm-preserving re-install bitwise. Write the recomputation as the SAME
      expression in the SAME order as item 21's, for the same reason
      `chain.hpp:718-720` gives.

  TESTS
  23. The all-basis K-forest fixture 2a says does not exist. Gaussian, at
      K = 2 and K = 3, through `dbartsData(x, y, bases = list(...))` with
      every entry non-NULL, asserting both forests resolve to the
      fixed-variance channel (`params[[f]][7] == 0`, `[[f]][5] == 0.674`) and
      that the sampler runs. **It must carry a NON-UNIT-NORM basis in at
      least one cell**, so 2f's sensitivity is PINNED rather than assumed -
      and under item 21 that cell also pins `nbar_f != 1` behavior, which is
      the only place in the suite that does.
  24. **Item 22's staleness closure has NO gate as the slice stands.** No
      equivalence leg calls `$setForestBasis` at all (`bcf-equivalence.R`
      mutates z only through creation), and `test-forest-basis-r5.R` asserts
      widths and amplitudes and never a calibration, so item 22's correctness
      AND its bitwise claim would ship on argument. Gate both, in that file, on
      a probit K-forest. (i) STALENESS: swap a forest to a basis whose median
      nonzero row norm is 4x the old one and assert
      `$getCalibration(f)[, "prior.scale"]` moved to
      `factor_f / (divisor_f * nbar_new)`. RED against a divisor left at its
      construction value, which is the defect item 22 exists to prevent.
      (ii) BITWISE: re-install the SAME basis values mid-run and assert
      `identical()` on `prior.scale`, and on `$getForestFits`,
      `$getForestAmplitudes` and `sigma` after a fixed-seed `$run`, against a
      twin that never swapped - that is the claim item 22's "retain `s` rather
      than recompute `scaledResponseSd()`" rests on, and nothing else in the
      suite or the trio can see it. 14-24 dense.
  25. **THE ANCHOR'S OWN GATE (2e).** The only assertion in the slice that can
      fail on a wrong `s`, and it is exact rather than statistical. Under a
      latent family `fitScale() == 1` and `fitShift() == 0`
      (`model.hpp:3137-3138` probit, `:3616-3617` logistic), and
      `priorScaleFactor` (`chain.hpp:3472-3474`) cancels the `sqrt(m)`
      `buildBCFForest` divides in at `:4507-4508`, so
      `$getCalibration(f)[, "prior.scale"]` (`chain.hpp:1005`) IS the map's
      `nodeScale_f = factor_f * s / (divisor_f * nbar_f)`, in latent units, with
      no conversion at all. Three assertions over item 23's all-basis fixture
      carried to probit and logistic:
      - **(a) BASE-RATE INVARIANCE, and it must LEAD.** Two samplers differing
        only in the response's base rate (p = 0.5 and p = 0.1) report
        `identical()` `prior.scale`. This is first because 2c's failure
        signature defeats the obvious fixture: probit's naive `s` equals Option
        L's `1.0` at p = 0.5, so a balanced fixture passes the defect.
      - **(b) THE ABSOLUTE ANCHOR.** `prior.scale_f` equals
        `factor_f * s / (divisor_f * nbar_f)` with `s` written as the LITERAL
        `1.0` / `pi/sqrt(3)` (never read back from the sampler) and `nbar_f`
        recomputed R-side as the median nonzero row norm, on a fixture carrying
        a non-unit-norm cell.
      - **(c) THE INDUCED INDEX**, contract sentence 2 made executable:
        `sqrt(sum_f prior.scale_f^2 v_f ||B_f(i,.)||^2)` equals
        `1.04912 sqrt(K) s` at unit row norms, at K = 2 and K = 3 - `1.4837 s`
        and `1.8171 s` - which is contract sentence 4's own arithmetic.
      RED on each of the three wrong anchors, by construction: the naive
      cold-start `s` (0.600x at p = 0.1, and (a) catches it with no tolerance at
      all), Option C (exactly 3x, on every arm and in both families), and a
      dropped `nbar` divisor (4x on the fixture's non-unit-norm cell, which (b)
      and (c) both see). GREEN on the shipped design. Tolerance `1e-12` relative
      on (a)/(b) and `1e-4` on (c), and it is NOT vacuous at that width: there
      is no Monte Carlo anywhere in the test, so no simulation error sets a
      floor; `1e-4` is there only because `1.04912` is a six-digit rounding of
      `sqrt(0.5)/0.674`; and the SMALLEST discriminating gap is 40%, four orders
      above it. 18-30 dense.

  **5. The seven forks, RESOLVED.** Options, choice, argument.

  **F1 - the response switch: INLINE, not factored.** The alternative is
  lifting a shared `makeResponse(family, ...)` out of the single-forest
  constructor's family switch (`chain.hpp:537-577`) and calling it from both,
  at ~30-35 dense
  against the inline 10-14. **Take inline**, on two grounds. First, the
  factored form puts the diff INSIDE the constructor `bcf-equivalence` covers
  bitwise, and the M4.1 reviewer's standing trap applies: a perturbed install
  has reported bitwise-identical WITHOUT `--preclean`. A mechanical move
  should be bitwise, but "should be" is what that trap eats. Second, the
  dedup is smaller than it looks: the single-forest switch is SIX arms
  carrying option plumbing the K-forest path deliberately does not have
  (`residualDf`, `survivalStatus`, `numCategories`, `dispersion`), so a shared
  helper would take four parameters the K-forest caller passes as dead
  defaults and three arms it can never reach. Inline is also what the price
  table quotes; the factored upper bound (engine ~50-75) is recorded as a
  contingency only.

  **F2 - `setResponse` against forest 0: (m-ii), fix both sites and drop the
  conjunct.** (m-i) keeps the refusal shut, which is cheapest and ships a
  latent K-forest that runs standalone and cannot be composed - a direct
  contradiction of the only justification this slice has left, which is
  composition ergonomics. **Take (m-ii).** They are ONE change, not two: with
  the conjunct dropped and `:1408` unfixed the slice SHIPS a probit K-forest
  whose latents are refreshed against the prognostic forest's bare fits after
  every `$setResponse`, silently, with no crash and no refusal. And with the
  residuals defect folded in (section 6 below), `:1514` is being edited
  anyway, so `:1408` is the same substitution in the same hunk. Note that the
  two sites are NOT symmetric today: `:1408`'s pointer is inert under gaussian,
  because
  `GaussianResponse::setResponse` (`model.hpp:2807-2830`) leaves its third
  parameter UNNAMED and neither branch of its body reads it - so `:1408` is
  genuinely M4.4's own, while `:1514` is a live defect the slice inherits.

  **F3 - the second bridge gate: SAME ALLOWLIST IN PLACE, do not delegate.**
  The pricing census called delegating to `refuseUnsupportedBCFComposition`
  "better", on the reading that it deletes a duplicate check. **That reading
  is wrong and the census did not verify it.** `createBCFHolder` does NOT call
  `refuseUnsupportedBCFComposition` at all - only `createHolder` does, at
  `:2916`. The two checks are not duplicates: `refuseUnsupportedBCFComposition`
  is a SUPERSET carrying FOURTEEN further refusals (DART, split probabilities,
  monotone, linear/GP leaves, k hyperprior, non-default k, node scale, named
  `prior.scale`, proposal probabilities, Student-t, grouped, variance forest,
  fp32, test predictors, per-column cut counts). Delegating would therefore
  ADD all of those to the internal route, which is not a cleanup but a
  behavior change - and that route is what the ENTIRE gate infrastructure runs
  on: `benchmarks/R/bcf-equivalence.R` (all 12 scenarios), `bcf-exact.R`,
  `bcf-exact-restricted.R`, `bcf-exact-weak.R`, `benchmarks/R/sbc.R:991`, and
  ELEVEN tinytest files (`test-bcf-creation.R`, `test-bcf-mutation-pins.R`,
  `test-bcf-zero-multiplier.R`, `test-bcf.R`, `test-blocks.R`,
  `test-forest-weights.R`, `test-interactions.R`, `test-multi-forest-seam.R`,
  and the three `test-multinomial-*.R`). Turning the slice's own bitwise gate
  red for reasons unrelated to the slice is exactly the failure this costs. It
  would
  also require reordering `optionsFromParsed` above the check, since
  `createBCFHolder` builds `options` after it. **Keep the two gates separate,
  pay the extra 1-2 dense, and record the divergence** (the internal route
  does not run the offender cascade) as a named item that is NOT M4.4's.

  **F4 - the `bcf()` family: READ IT OFF THE MODEL, taking neither stated
  option.** The fork was posed as "thread a family string through
  `createBCFHolder` / `bartcore_createBCF` / `bartcoreBCFSampler`" versus
  "document the route as gaussian-or-probit-only". **Both are worse than a
  third answer the sources did not see:** `dbartsModel` already carries a
  `family` slot (`R/A_class.R:400`, default `"auto"`, validated `:451-467`),
  and `bartcoreBCFSampler` already passes `sampler$model` as
  `bartcore_createBCF`'s second argument (`R/bartcore.R:690`). So
  `createBCFHolder` should DERIVE its family name from `modelExpr`'s `family`
  slot, applying the same `"auto" -> ""` mapping the R5 route applies at
  `R/dbarts.R:792`, instead of hardcoding `""` at `:2990-2991`. No `.Call`
  arity move (`DEF_FUNC("dbarts_bartcore_createBCF", bartcore_createBCF, 8)`,
  `src/R_interface.cpp:182`, stays at 8), no three-file edit, and nothing in
  `dbarts.h`.

  Three arguments. **(a) Value, named.** The internal route is not a shorthand
  whose reach is optional - it is the route every K-forest GATE is written
  against (F3's list). A family the gate harness cannot construct is a family
  whose bitwise, exact-posterior and SBC coverage cannot be written in the
  shipped harness shape; with this fix, adding a probit `bcf-equivalence`
  scenario is a one-line change to the host sampler. **(b) It removes a class
  of mismatch rather than adding one.** Threading a separate string lets a
  caller state the family twice and disagree with itself; reading it off the
  model makes agreement structural. Today the two ALREADY disagree - passing a
  probit host sampler makes `resolveFamily` see `""` on a binary response,
  resolve to probit, and get refused one line later - which is precisely the
  silent shape-decides-the-family behavior `resolveFamily`'s own design comment
  (`:1551-1557`, "Each shape refuses the others' family names by name") exists
  to prevent. **(c) It is a strict superset with an identical default.** The
  COMMON case is `"gaussian"`, not `""`: `dbartsSpec` already resolves the
  family away from `"auto"` (`R/spec.R:85-86`) and writes the resolved value
  into the model (`:290`), so every existing gate and tinytest derives
  `"gaussian"` - and `resolveFamily(control, "gaussian")` on a continuous
  response returns gaussian, so behavior is identical. The rare hand-built
  `"auto"` model maps to `""`, also identical.

  **The sibling convention, which the fork's sources did not raise and which
  F4 must not break by omission.** The other two internal creation helpers
  take an EXPLICIT family beside the model that already carries one:
  `bartcoreSampler(sampler, family = "")` (`R/bartcore.R:541`) and
  `bartcoreSamplerFromHandle(..., family = "", ...)` (`:584-591`), and
  `R/bart.R:1523` uses it (`family = "ordinal"`) against a spec that set
  `family = "ordinal"` at `:1511`. Under F4 alone `bartcoreBCFSampler` would
  become the ONLY internal creator whose family is invisible at the call site,
  which is precisely the payoff sentence in (a) read as a cost. **Take both,
  at ~2-3 R dense:** the C-side derivation from the model slot stays (it is
  what covers the flat and internal routes and makes agreement structural),
  AND `bartcoreBCFSampler` gains `family = NULL`, which when supplied is
  written into the model object passed as `.Call` argument 2 before the call
  (`model <- sampler$model; if (!is.null(family)) model@family <- family`).
  One source of truth C-side, visible at the call site R-side, no arity move,
  and an unknown name is still refused by `resolveFamily`. Argument (b)'s
  "state it twice and disagree" objection does not apply, because there is no
  second channel to disagree with.

  **F5 - the `node.scale` gate: KEY IT TO THE FAMILY DEFAULT (option a).** The
  alternative (b) is to record whether the user NAMED a node scale and gate on
  that. It is truer to the guard's stated intent, and it is still the wrong
  change here. **The monotone precedent at `R/model.R:394` does NOT transfer**
  and the implementer must not pretend it does: `node.scale` has no unsupplied
  arm to redirect - it is written UNCONDITIONALLY by the family switch at
  `R/spec.R:297-309`. The argument for (a) is independent and it is about what
  the shipped guard already MEANS. `!= 0.5` is not a test for "the user named
  something"; `0.5` is GAUSSIAN'S DEFAULT, and the guard has always read
  "differs from the family default" - the literal was inlined only because the
  family was always gaussian. (a) generalizes the rule the guard already
  encodes. (b) changes it, and changes it for gaussian too: a gaussian caller
  who explicitly writes `node.scale = 0.5` is silent today and would newly be
  refused, which is a behavior move on the shipped path for no gain to this
  slice. (a)'s one false negative - a caller who explicitly names exactly the
  family default gets silence, and the map ignores their number anyway - is
  the same false negative gaussian has shipped since the guard was written.

  Mechanically: the family-keyed switch is written out TWICE today
  (`R/spec.R:297-309` and `R/xbart.R:279-284`, the latter a three-family
  subset), so this wants a shared helper rather than a third copy. Introduce
  `defaultNodeScale(family)`, use it at both existing sites and at
  `R/spec.R:452`, and give the bridge twin at
  `R_interface_bartcore.cpp:2311` the matching C-side switch.
  `R/spec.R:456`'s `prior.scale` gate (`!is.na(model@prior.scale)`) is CORRECT
  as written for all three families and needs no change - it stays FALSE
  unless a calibration is NAMED.

  **F6 - per-forest weights under a latent family: ACCEPT, do not refuse; one
  Doxygen sentence.** The pricing census worried this composes "precision-weight
  semantics on a family whose weights are declared to mean COPIES". **Reading
  the channel's own declared contract dissolves the worry.** `setForestWeights`
  (`chain.hpp:981-986`, Doxygen `:959-980`) declares `s` as "a multiplicative
  PRECISION factor on forest f's own leaf conditionals, composing with the
  observation weight so that forest f's draws see `w_i m_f^2 s_i`", explicitly
  NOT a removal from occupancy, from the combination, or from the residual
  sigma degrees of freedom. It is consumed at exactly one place,
  `composeForestWeights` (`chain.hpp:3485-3493`), which multiplies it into the
  vector the combiner just formed - DOWNSTREAM of `refreshLatents`. It never
  reaches `LogisticResponse`, so it cannot move the PG shape and cannot be
  confused with the copy count: by the time `formForestResponse` reads
  `workingWeights()`, logistic's `w` is `omega_i`, a PRECISION
  (`model.hpp:3529-3531`), and the copy count has already been spent as the PG
  shape parameter. Every family serves a working precision; the channel is
  declared against a working precision; there is no semantic clash to gate on.
  Nor is it an end-run around `refuseBinaryWeightChange`, which guards the
  observation weight (the creation-fixed copy count), a different quantity.
  Cost 0 dense. **The obligation is a Doxygen sentence at `setForestWeights`
  that carries the COHERENCE point as well as the copies-versus-precision
  one.** `composeForestWeights` applies the factor to the LEAF conditionals
  only; it is NOT applied in `drawGlue` / `drawForestAmplitude`, which read
  the response's own weights (`combiner.hpp:1057`, `:1087`, `:1171`), and NOT
  in `rescaleAmplitudeRidge`, which reads no weights at all. So `m_f` and
  `f_f` - the two factors of one product - are drawn under three different
  precision vectors. That is shipped behavior under gaussian; under a latent
  family it is newly reachable, on the one channel whose Doxygen this slice is
  rewriting anyway. Plus one line in the docs delta.

  **F7 - the "half the mantissa" claim: NARROW THE COMMENT.** The alternative
  is re-deriving a small-`omega` quantile bound on `PG(1, psi)` for the
  logistic arm. **Narrow it.** Decompose what `combiner.hpp:800-802` asserts:
  the amplification cap ("the division amplifies by at most 2^26") is a
  property of the constant `combiner.hpp:787` and the snap branch `:823-827`,
  holds under any family, any weight vector and any K, and transfers verbatim.
  The absolute-precision conclusion ("the cancellation downstream stays inside
  half the mantissa") additionally needs a bounded NUMERATOR, which today comes
  from gaussian range-anchoring `workingResponse()` to `[-0.5, 0.5]`
  (`model.hpp:2791`) - probit's is a truncated normal on the latent scale and
  logistic's is `O(1/omega_i)` with `omega_i` unbounded below, so neither
  family supplies it by construction. Re-deriving the logistic bound is genuine
  analytic work with NO gate consequence, because the amplification cancels
  ANALYTICALLY in what the node kernels accumulate: `sumWeights = sum_i w_i
  m_i^2` and `sumWeightedResponse = sum_i (w_i m_i^2)(r_i/m_i)`, whose exact
  value is `sum_i w_i m_i r_i`; and `combinedFits` (`:854-869`) and
  `drawForestAmplitude` (`:1153-1198`) deliberately keep the EXACT multiplier
  (`:851-853` states this), so nothing downstream inherits the amplified value.
  Keep the family-free half, drop the absolute-precision conclusion, and state
  the cancellation instead. The small-`omega` bound is recorded as an open
  numerical-hygiene question, out of M4.4.

  **6. FOLDED IN: `$getSumsOfSquaredResiduals()` returns a wrong number today.**
  VD decision 2026-08-14: this is M4.4's, not a separate slice.
  `Chain::sumOfSquaredResiduals` (`chain.hpp:1511-1516`) reads
  `forests_[0].totalFits.data()` at `:1514`, and unlike its twin at `:1408` the
  pointer IS read - by `misc_computeSumOfSquaredResiduals`, unconditionally,
  for every family INCLUDING gaussian. The path carries no multi-forest guard
  anywhere: `$getSumsOfSquaredResiduals` (`R/dbarts.R:1440-1444`) has none,
  `bartcore_getSumsOfSquaredResiduals` (`R_interface_bartcore.cpp:4799-4808`,
  call at `:4805`) has none. Its documented contract is
  `sum((y - yhat)^2)` (`man/dbartsSampler-class.Rd:351`), and on a K-forest
  `yhat` is `sum_f m_f(i) f_f(x_i)`, not `f_0(x_i)`. MEASURED on four chains
  after `$run(20L, 5L)` on a K = 2 gaussian K-forest, with the combination
  reconstructed from public getters only: the reported value equals the
  forest-0 residual TO THE LAST PRINTED DIGIT on every chain - with the
  amplitude not even applied - and differs from the correct combined residual
  by 1% to 49%. **This is a shipped wrong answer on the public R5 surface,
  gaussian, today.** Nothing caught it because every test that exercises the
  getter (`test-sampler-residuals.R:17`, `test-data-mixed-mutation.R`,
  `test-bartcore.R:912`, `test-mutate-sparse-valued.R`) is single-forest, and
  grep finds no internal caller (`rbart`, `bart2`, `xbart` do not use it).
  Corroborating, and the shape of the fix: `storeSample` ALREADY routes the
  reported TRAIN channel through `combinedFits` (`chain.hpp:4636-4637`), so
  `$run()$train` is right on the same object where this getter is wrong.

  **The substitution is GUARDED, not blind, and here is what it does per
  response family.** `Chain::combinedFits()` (`chain.hpp:4589-4592`)
  dispatches virtually, and on a MULTINOMIAL chain (`chain.hpp:742-763` - a
  K-forest chain with `combiner_` set) `MultinomialForestCombiner::
  combinedFits` (`combiner.hpp:1725-1730`) returns the **K x n softmax
  PROBABILITY slab**, channel-major, not a location and not length n. A bare
  substitution would silently change the number this getter reports there
  (`MultinomialResponse::workingResponse()` is n zeros, `model.hpp:3691`,
  `:3697`, so it would go from `sum(totalFits_0[i]^2)` to `sum(p_0i^2)`), at
  a site whose SIBLING carries an explicit layout guard warning about exactly
  this K-versus-1 channel confusion (`model.hpp:3701-3704`).
  - **gaussian, probit, logistic, and every single-forest chain:** take
    `combinedFits()`. Off a combiner it is a LITERAL identity
    (`chain.hpp:4590-4591` returns `forests_[0].totalFits.data()`), and on a
    BCF coupling it is the fix.
  - **multinomial (any multi-location combiner):** the quantity is not
    defined - there is no single location and no working response - so gate
    the substitution on `numReportedLocations() == 1`
    (`combiner.hpp:612`, overridden `:1493`) and return `NaN` on the other
    branch, with a Doxygen sentence saying why. Today's value there is
    equally meaningless, so nothing honest is lost: the path is unreachable
    from R (`bartcoreMultinomialSampler`, `R/bartcore.R:792-828`, returns a
    bare environment, not a `dbartsSampler` R5, and
    `$getSumsOfSquaredResiduals` `R/dbarts.R:1440-1444` is the only R
    caller) and absent from `dbarts.h` and `C_interface.cpp` entirely, but it
    IS reachable from C++ through `facade.hpp:516-517`, and the two
    `tests/cpp` users
    (`test_grow.cpp:373` finiteness, `test_sampler.cpp:2021` self-recomputed)
    are both single-forest gaussian and unaffected. Do not leave this
    implicit.

  **Its own acceptance criterion, separate from the rest of the slice.** Take
  the CHEAP oracle, not a reconstruction: `storeSample` routes train through
  `combinedFits` and runs last in the sweep, so on a GAUSSIAN K-forest

      res <- s$run(nBurn, 1L)
      expect_equal(s$getSumsOfSquaredResiduals(),
                   colSums((y - res$train[, 1L, ])^2))

  is the shape `inst/tinytest/test-sampler-residuals.R:16-19` already ships
  (17 non-comment lines for the whole file), it needs no basis, and it is
  MEASURED to fail on `efec6ba2` (K = 2, four chains: reported 548.0 / 447.3 /
  240.4 / 256.5 against the oracle's 201.7 / 197.1 / 204.6 / 204.3; K = 3 with
  a two-column basis, two chains: 472.6 / 761.1 against 203.8 / 200.3). At
  K = 3 with a multi-column basis, run it. `res$train` is `n x nSamples x
  nChains`, so the `[, 1L, ]` drop is required at `nSamples = 1`; at one chain
  it is `n x 1`. An offset is FINE and needs no special handling - the
  recorded train channel includes it by the original-scale convention and the
  identity still holds, as `tests/cpp/test_sampler.cpp:2005-2022` pins in C++
  with an offset present. **The reconstruction oracle this bullet used to
  mandate is WRONG and must not be written:** `$getForestAmplitudes()`
  (`R/dbarts.R:1450-1458`) returns AMPLITUDES, the ragged stacked
  `sum(q) x n.chains` block, while the combination uses the per-observation
  multiplier `m_f(i) = dot(a_f, B_f(i,.))` (`combiner.hpp:854-869`), and NO
  getter exposes the basis. A literal `sum_f a_f * fits_f` is accidentally
  forgiven by the K = 2 two-column indicator and is simply wrong at the K = 3
  the criterion mandates - and "must FAIL before, pass after" does not
  discriminate, because a wrong oracle also fails before and can be "fixed"
  by matching the buggy getter. It must FAIL on `efec6ba2` and pass after.
  This assertion is gaussian-only and does not depend on any other part of
  M4.4, which is the point: it is the evidence that the fix landed independent
  of whether the family work does.

  **The getter's CONTRACT is still wrong on the surface M4.4 newly opens, and
  the docs delta must pick it up.** `man/dbartsSampler-class.Rd:351` promises
  `sum((y - yhat)^2)` on the original scale, "a binary-response sampler
  reports on the latent scale instead". That covers probit. It does NOT cover
  logistic: `workingResponse()` is `w_i(y_i - 0.5)/omega_i - offset_i`
  (`model.hpp:3609-3613`) with `omega_i` a fresh PG draw each sweep, so the
  returned number is a function of this sweep's auxiliary variables and is not
  a residual sum of squares in any scale. Say so in the Rd rather than letting
  M4.4 make an unstated case reachable.

  **7. Docs deltas M4.4 inherits, by name.**
  - `inst/include/dbarts/dbarts.h:513` - "Gaussian responses only.", the LAST
    sentence of `dbarts_sampler_create`'s K-forest paragraph (`:495-513`) and
    the ONLY sentence IN THE HEADER asserting a gaussian-only K-forest
    constraint (the other gaussian hits at `:129`, `:282`, `:488`, `:489`,
    `:542`, `:544`, `:555`, `:566`, `:570` are single-forest or
    family-agnostic). RE-VERIFIED at `:513` (this bullet used to cite `:505`;
    it moved once already, by the paragraph's own growth at M4.5, which is
    exactly why it must be re-derived and not offset). Also `dbarts.h:560`
    (`dbarts_sampler_setWeights`, "gaussian responses only"), which becomes
    ambiguous once logistic K-forest weights are accepted at creation, and the
    family paragraph at `:485-493`, which describes families only for the
    single-forest case. Comment-only; the hash does not move.
  - **NINE further gaussian-only K-forest assertions OUTSIDE the header, all
    falsified by M4.4.** The header's uniqueness claim above is true of the
    header and was wrongly read as an inventory of the repo.

    | site | what it says |
    |---|---|
    | `chain.hpp:688` | the K-forest `Chain`'s own Doxygen, "gaussian response as y = sum_f ..." |
    | `sampler.hpp:150` | the K-forest `Sampler`'s own Doxygen, "gaussian only" - and still K = 2 language |
    | `combiner.hpp:679` | `BCFForestCombiner`'s Doxygen, "combined on a gaussian response" |
    | `combiner.hpp:953-955` | the GOVERNING RATIONALE for `supportsResponseMutation() == true` (`:963`); see checklist item 6 - it must be re-argued, not swept |
    | `R/bartcore.R:612` | "internal and gaussian only" - F4's own site |
    | `R_interface_bartcore.cpp:3522` | "internal, gaussian only" - F4's own site |
    | `man/dbartsSampler-class.Rd:147` | the R-side twin of `dbarts.h:560` |
    | `man/dbartsSampler-class.Rd:152` | `setActiveRows`, "its two forests share one gaussian response" |
    | `docs/design/feature-matrix.md:496-497` (`[f26]`) | "A BCF sampler's response IS a `GaussianResponse`" |
  - `docs/design/multiplier-combiner.md` - the scope sentence at `:14-17`
    ("Scope: GAUSSIAN responses only", naming `R/spec.R:423-431` and
    `chain.hpp:702-705`) becomes gaussian, probit and logistic; the map
    section `:321-341` gains the anchor table, items 2d(a)/(b) and 2f's four
    contract sentences; the Status section `:451-455`, which names
    `dbarts.h:513` as M4.4's, is discharged. **Plus M4.4's own landing note**,
    in the file's per-slice format - every other slice on the arc has one and
    this bullet had not booked it.
  - **The "sd(y) units" phrasing is a REQUIRED edit, not an optional sweep,
    and it is NOT confined to the design docs.** Under a latent family `sd(y)`
    is UNDEFINED, so every one of these becomes WRONG rather than merely
    incomplete: `multiplier-combiner.md:325-335` and its third instance at
    `:25` ("the sd(y)-unit calibration"); the Doxygen twins
    `combiner.hpp:262-273` and `:297-299`; `facade.hpp:774-776`
    ("constant-leaf and gaussian only"); and, on the surface a USER actually
    reads and where the units are actually chosen -
    **`man/forest.Rd:52-59`** (`\item{sd}`, "This forest's prior scale, in
    standard deviations of the response"), **`R/model.R:853`** ("the node
    scale stays at the response sd"), **`:856`** ("sits at sd sd(y)"),
    **`:1428`** ("in sd(y) units") and **`R/bartcore.R:616`** ("in sd(y)
    units"). `man/forest.Rd:52-59` is the worst of them: `forest(sd = )` is
    the ONE shipped knob, a probit K-forest caller writing `sd = 1.5` reads
    1.5 response standard deviations and gets 1.5 units of the LINK's error sd
    with no response sd existing anywhere in the model, and the deferred
    `binary-kforest-prior-default` slice names that same knob as its
    correction lever. `R/model.R:851-856` is `forestParams`' Doxygen - the
    function that WRITES `nodeScaleFactor` and `amplitudePriorScale`, the two
    numbers `latentScaleAnchor(family)` multiplies - so leaving it saying
    "stays at the response sd" guarantees the next reader reconciles the
    payload against an anchor the engine no longer uses. All five take 2f's
    contract sentences 1 and 3.
  - `docs/design/nameable-calibration.md` - two new rows after the units
    table's last row at `:76` (`| BCF | response; k fixed at 1 by the map |
    range_ | range_*0.5 + min_ |`), for BCF-probit and BCF-logistic, both
    `1` / `0`. Plus the sentence that `prior.scale` is the prior sd of forest
    f's own `f_f`, NOT of its contribution `m_f f_f` to the index - already
    true of the gaussian K-forest, but under a latent family the index is the
    only scale in the model, so the gap is the difference between a number a
    reader can act on and one they cannot. Point at `$getForestAmplitudes()`
    as the other half. **Anchor correction while here:** "refused at three
    creation sites and again mid-chain" is at `:134-136`, inside the refusals
    section `:125-145` (`## 5. Refusals` at `:125`, item 7 ending at `:145`,
    `## 6. Exactness` at `:147`); the respec's `:125-136`, the open-questions
    file's `:132-134` / `:124-140`, and this bullet's own earlier `:125-140`
    are all WRONG, verified by opening the file.
  - `docs/design/feature-matrix.md` - **MANDATE, and M4.4 cannot discharge it
    the way M4.5 did.** M4.5's landing note records "feature-matrix.md moves NO
    CELL: the family is still gaussian-only, so the new K-length creation route
    reaches no new response family". **M4.4 is exactly what changes that**, so
    it MUST move cells and update the matrix in place. Cited by ROW and COLUMN
    rather than by line, since the file's own edits shift its lines. At
    minimum: the `bcf` row's **SEVEN** section-2 mutation cells (`:119`, not
    six - `setResponse`, `setOffset`, `updateScale = TRUE`,
    `setPredictor (+ per-obs)`, `setWeights`, `setSigma`, `test surface`)
    become family-dependent (`setWeights` and `setSigma` are REFUSED for
    probit/logistic through `refuseBinaryWeightChange` and
    `refusePinnedSigmaChange`, while `setResponse`/`setOffset` OPEN under item
    5's F2, and `updateScale = TRUE` takes the shipped latent convention
    `- [f9]`, ignored rather than refused, since latent families have
    `fitScale() == 1` / `fitShift() == 0`); the `bcf` x `getLatents` cell in
    section 3 moves off `-` for the latent families, and footnote `[f18]`
    ("No latent vector exists: gaussian, BCF and heteroscedastic all
    leave...") must be rewritten; `bcf` x `zero-weight row subset` (`:142`)
    becomes family-dependent and `[f17]`'s closing "The same fix covers BCF
    and the Student-t row" is gaussian-scoped; `bcf` x `active-rows mask`
    (`:142`) carries `[f26]`, which is the ninth gaussian-only assertion
    above; the `bcf` x `variance forest` cell in section 4 gains checklist
    item 7's door refusal beside it; the `bcf` x `DART` cell cites
    `spec.R:425`, which M4.4's own edit at `:423-431` moves (and which is the
    wrong anchor anyway - the DART refusal for a K-forest is the `unsupported`
    entry at `R/spec.R:441`, fired at `:468-474`); `[f23]` (`:467-475`) cites
    `spec.R:440` for the R-side BCF composition refusal, which the same edit
    shifts; `[f3]` (`:204-209`, "The header's `family` documentation
    (CAPI:338-357) names only probit, logistic, gaussian, aft and BCF") is
    exactly what M4.4 changes AND its anchor is already wrong -
    `dbarts.h:338-357` is the `DBARTS_C_API_LIST` X-macro body and the family
    documentation is `:485-493`; the `bcf` row in section 5 gains the new
    probit/logistic scenarios and test files; the Rows prose citing `MOD:2524`
    for `enum class ResponseFamily` is stale (`model.hpp:2577`) and M4.4
    changes that enum's REACH, so it re-anchors; and the section-1 prose on
    what `dbartsSpec()` resolves gains the family statement. Add an "Exception
    on record" entry in the file's own format naming what M4.4 re-verified BY
    SYMBOL and what it did not, and do NOT bump the Status line's "current at"
    stamp unless a whole-file pass actually happens - the outstanding full pass
    is the root TODO's `feature-matrix-anchor-refresh`.
  - `man/dbarts.Rd:84` (inside `\item{forests}`) and `man/forest.Rd:93` both
    carry "Gaussian responses only." Rd is not air-wrapped, so raw is about
    dense there.
  - **State plainly rather than inherit silently:** BCF leaves
    `testFitsAreDefined()` and `logLikelihoodIsDefined()` false
    (`combiner.hpp:945-946`), so a probit/logistic K-forest sampler still
    refuses the whole test surface and reports no log-likelihood. Unchanged by
    M4.4, and a reader of the new family will look for it.

  **8. Gates.**
  - **`bcf-equivalence` MUST stay bitwise on the gaussian arm, all 12
    scenarios, all channels.** Re-check under `--preclean` regardless of which
    variant of checklist item 2 is taken:
    `docs/design/multiplier-combiner.md:415-417` records a perturbed install
    reporting bitwise-identical WITHOUT it. The
    single-forest `equivalence` (37 scenarios, `--strict-coverage`) and
    `multinomial-equivalence` (10) legs must stay bitwise too.
  - **The residuals fix IS draw-neutral, and the STATIC argument is the
    proof.** The reason matters, because the reason this bullet used to give
    was wrong. "It is a getter, so it cannot shift draws" does NOT hold:
    `BCFForestCombiner::combinedFits` (`combiner.hpp:854-869`) is not side
    effect free - it resizes and WRITES `glue_.combined` and
    `glue_.fitsByForest` (declared `:435-438`), combiner-owned scratch the
    sweep also uses, and returns a pointer INTO the first, which the sweep
    holds live from `chain.hpp:1303` across `refreshLatents` (`:1304`),
    `drawSigma` (`:1314`) and `drawGlue`/`afterCombine` (`:1317-1318`). The
    property that actually holds is three-part and each part is checkable:
    (a) both buffers are FULLY overwritten on every call (`:857-861` sizes
    them, `:863-868` writes every element), so no stale value can be read;
    (b) there is NO reachable call site inside `run()` or
    `growForestFromRoot` - `Chain::sumOfSquaredResiduals` has exactly one
    caller chain and it is entirely outside the sweep, `sampler.hpp:882-884`
    -> `facade.hpp:516-517` (virtual, declared `:262`) ->
    `R_interface_bartcore.cpp:4805`, and `Chain::setResponse`
    (`chain.hpp:1407`) likewise only `facade.hpp:391-392` (both `facade.hpp`
    hops were missing from this bullet's prescribed grep); (c) neither branch
    of `chain.hpp:4590-4591` consumes rng. Record (a) and (b) explicitly, not
    "consumes no rng" alone: a future combiner that CACHED across calls, or a
    future in-sweep caller, breaks the weaker claim silently.
  - **The equivalence trio here is a REGRESSION CHECK, not the proof, and the
    bullet used to have that ranking backwards.** Run it - under `--preclean`,
    on a build carrying ONLY checklist items 5 and 6, before the family switch
    exists - but do not treat green as evidence, because it cannot come back
    anything else: neither fix site is reachable from `run()`; on the 37
    single-forest scenarios the substitution is a LITERAL identity
    (`chain.hpp:4590-4591`); no harness records the channel (zero hits for
    `SumsOfSquaredResiduals` anywhere under `benchmarks/`); and
    `bcf-equivalence` never calls `$setResponse` at all. What DOES have teeth
    is the belt: `$getSumsOfSquaredResiduals()` is on no equivalence channel,
    so assert that a fixed-seed gaussian K-forest's `$getForestFits`,
    `$getForestAmplitudes` and `sigma` after `$run` are unchanged against the
    pre-fix build. **The `combinedFits()` FMA-association contract is still
    the live hazard for anything that RECOMPUTES the blend inside the sweep**
    (`combiner.hpp:840-849`: accumulating forward instead of seeding with the
    last forest moves ~30% of rows by one ulp and turns all 12 scenarios red),
    which is why item 21's map expression must be written in the stated order.
  - **The 2f scale fix gets its OWN trio run, and it is a real gate rather
    than a formality**, because item 21 changes an expression the ctor
    evaluates on every route. Build items 21 and 22 alone, run
    `bcf-equivalence` (12), `equivalence` (37, `--strict-coverage`) and
    `multinomial-equivalence` (10) under `--preclean`, and expect bitwise on
    the `nbar == 1.0` argument in 2f. **If any leg moves, take 2f's
    pre-registered option-A fallback** - do not re-record to accommodate a
    fix that was chosen for being free. Also re-run
    `inst/tinytest/test-forest-basis-r5.R`, whose `:47` (`wide <- cbind(1,
    x[, 4L])`, row norms in `[1, sqrt(2)]`) and `:101` (`cbind(1 - z, z,
    x[, 1L])`) fixtures change MODEL under item 21; neither pins a value, so
    both should still pass, and a failure there is a finding, not a snapshot
    to regenerate. **The trio does not reach item 22 at all** - no equivalence
    leg calls `$setForestBasis` - so its staleness closure is gated by
    checklist item 24 and by nothing else.
  - **If it turns out NOT to be draw-neutral,** the slice does not promise a
    bitwise verdict and instead RE-RECORDS. Re-record `bcf-equivalence`
    (`benchmarks/baselines/bcf-equivalence-8b047f8b.rds`) and, only if the
    single-forest or multinomial legs also move, `equivalence-8b047f8b.rds`
    and `multinomial-equivalence-1027be5.rds`. **A re-record touches FOUR
    places, not one, and all of them land IN THE SAME COMMIT as the new
    `.rds` files** - a workflow pointing at a deleted baseline is a red gate
    that looks like a regression, and a ledger naming a deleted file lies
    about which baseline is current.
    (1) the CI workflow lines `.github/workflows/equivalence.yaml:61`, `:87`,
    `:113`; (2) **`benchmarks/baselines/MANIFEST`** (`:15`, `:16`, `:42`,
    `:48`) - the AUTHORITATIVE ledger, with a `current`/`historical` role
    column and a per-baseline narrative entry recording scenario counts, the
    neutrality partition and the superseding hash, in the format the existing
    entries set; (3) `TODO:258`, `:387`, `:389`; (4)
    `docs/design/feature-matrix.md:627-628` (`[f39]`, "Current baselines"),
    plus `:308`, `:441`, `:507`, `:755`. **The same MANIFEST obligation
    attaches to F4's own stated payoff** - adding a probit `bcf-equivalence`
    scenario is a baseline re-record, and if the slice takes that payoff it
    must budget it rather than treating it as a one-line change. Any
    tinytest with a hardcoded expected value that shifts is regenerated by
    REPLAYING THE WHOLE TEST FILE, since the values depend on the file's full
    execution history and not just the preceding seed. Whichever way it goes,
    the finding is recorded in the landing note explicitly - "draw-neutral,
    proven at commit X" or "NOT draw-neutral, baselines re-recorded" - because
    the next slice will inherit the assumption either way.
  - **M4.4's OWN acceptance gate is arm E**, specified in the FA5 block below:
    the K-forest probit sampler, compared against FA5's decisive arm B (K
    gaussian samplers with host-drawn latents against the combined fit), which
    AGREED with the reference on all 12 functionals at max |z| = 2.54 against
    a threshold of 3.0. The harness exists,
    `.claude/m4-basis-design/harness/fa5-latent-coupling.R`, and arm E is a
    new arm in it. It must be a STATISTICAL comparison and not a unit test,
    because the open question it answers - whether the pinned sigma costs
    mixing on the FIXED-VARIANCE basis forest, which under gaussian would be
    absorbed by sigma growing - can only be answered by running it. The
    pre-registered rule already stated for arms A and B applies unchanged.
    **Arm E is NOT the anchor's gate and the two are not substitutes:** it
    measures mixing and coverage under a pinned sigma, with no measured power
    against a mis-scaled prior. Checklist item 25 is the anchor's gate.
  - Positive builds, not just the absence of the old refusal. **An implementer
    who tests the family relaxation by deleting `expect_error` and never
    fitting will not see the adjacent-guard misfire**, which is invisible
    until the first fit: relaxing `R/spec.R:424` and bridge `:2299` alone
    produces a sampler that refuses on EVERY probit fit. Both guards were
    confirmed to fire, for both families, by reading the expressions AND by
    running real probit and logistic samplers (probit `node.scale` 3, logistic
    5.441398, both `node.hyperprior` a `dbartsChiHyperprior`, with the user
    having asked for nothing). **Do not grep for a quoted prefix:**
    `unsupported` is a named logical vector (`R/spec.R:440-467`) and
    `:469-473` pastes EVERY true name, so with both `:446` and `:452` firing
    the actual message is "a treatment forest does not support a 'k'
    hyperprior, a non-default 'node.scale'; drop it or fit a single-forest
    model".
  - **Name the silent behavior flip in the landing note.** After checklist
    item 15, `dbarts(x, y01, forests = list(forest(), forest(basis = ~
    factor(z))))` stops raising "a treatment forest requires a continuous
    (gaussian) response" and silently builds a PROBIT K-forest:
    `R/spec.R:85-86` resolves a numeric 0/1 response to `"probit"` with no
    message, and unlike the categorical branch (`R/data.R:461`, `:525`, both
    calling `announceAutoFamily`) the numeric branch announces nothing. This
    is longstanding shipped single-forest behavior, so it is a flip M4.4
    inherits rather than a new inconsistency - but a causal-forest user with
    a binary outcome who wanted a linear-probability two-forest fit now gets
    latent-scale amplitudes with nothing said.
  - Retarget the one existing pin, `inst/tinytest/test-bcf-creation.R:590-594`
    (`expect_error(..., "gaussian")`), at a DOOR family. It is the only test
    anywhere that pins this refusal; nothing matches "treatment forest
    requires" or "requires a continuous", and the two C-side doors at
    `:2300` / `:2994` are unpinned by any R test.
  - `R CMD INSTALL --preclean` is mandatory (headers move), and the
    `benchmarks/kernels` binaries must be deleted by hand - they carry no
    header dependency tracking.
  - Making `family_` honest at checklist item 8 flips what EVERY family-keyed
    bridge predicate answers for a K-forest sampler, all at once:
    `refusePinnedSigmaChange`, `refuseBinaryWeightChange`, the `drawsSigma`
    branch in `bartcore_setModel`, `validateResponseSupport` at `setResponse`
    (`:4424`, `:4469`, `C_interface.cpp:464`) and the aft branch at `:4449`.
    Each flip is believed CORRECT; the engine cost is one line and the test
    cost is not. Mildly de-risking, and worth knowing before writing the
    tests: they are not entirely unexercised today, because the MULTINOMIAL
    K-forest sampler already reports a non-gaussian family
    (`sampler.hpp:187`, `chain.hpp:753`), so those four predicates already see
    a non-gaussian K-forest. What is untested is the BCF-coupling x latent
    family combination specifically.

  **9. Budget, dense-equivalent, with STOP conditions at 1.5x.**

  | layer | budget | STOP at |
  |---|---|---|
  | Engine | 49-82 | 123 |
  | Bridge | 14-28 | 42 |
  | R | 18-38 | 57 |
  | Flat C | 0 | 2 |
  | Docs | ~100-175 | 265 |
  | Tests | ~155-260, plus arm E's harness arm | 390 |

  Arc arithmetic: `~393 + 82 = ~475` of the `~500-700` engine band at the top
  of the inline budget, `~393 + 123 = ~516` at the engine stop - so **even a
  1.5x engine overrun does not breach the band**, and the stop is a signal to
  re-scope rather than a budget breach. Say so plainly rather than letting a
  reviewer infer a crisis. On hitting a stop: HALT and report, do not push
  through. The named contingency for an engine stop is the probit-only v1
  (one fewer switch arm, one fewer anchor arm, bridge and R unchanged since
  the allowlist has two entries instead of three), which saves only ~5-10
  engine and is **NOT recommended** - logistic costs one switch arm, one
  anchor constant, and inherits the entire weight policy for free. The 2f
  scale fix (items 21-22, engine 20-27) has its own pre-registered fallback
  and is the other lever. A split is not required and none is recommended.

  **The test band is re-derived from the MANDATED ORACLE and from comparable
  LANDED files, not from the engine delta.** The earlier `~45-100 / stop 150`
  scheduled a HALT on this slice's NORMAL path and must not be restored.
  Repo calibration, non-comment lines, measured: one comparable new K-forest
  R5 test file in this repo lands at 107-234 (`test-bcf-reporting.R` 107,
  `test-forest-weights-r5.R` 138 for ONE new R5 method,
  `test-bcf-r5-surface.R` 167, `test-forest-basis-r5.R` 234, M4.3's one new
  whole file) - so the old band was less than one such file and its stop was
  below the largest. Arc test-to-engine ratios from this plan's own landing
  notes: M4.1 98/48 = 2.0x, M4.2 351/180 = 2.0x, M4.3 ~645/165 = 3.9x; the
  old band was 0.8x, and M4.4 is the FIRST slice on the arc to move engine AND
  bridge AND R AND the public family surface (M4.1 and M4.2 were engine-only,
  with zero tinytest). Bottom-up against the mandated work: residuals oracle
  at K = 2 and K = 3, every chain, in the cheap form section 6 specifies 15-25;
  the all-basis fixture (item 23) 12-20; retargeting the one existing pin at a
  door family 3-6; positive probit and logistic K-forest builds, which no
  longer carry the anchor assertion because item 25 does, 25-45; the five
  flipped
  family-keyed bridge predicates x 2 families 30-55; guard relaxations
  positive AND negative x 2 families (node.scale default silent / non-default
  still refused; default k silent / explicit `chi()` still refused;
  `prior.scale` unchanged) 20-35; door refusals aft/ordinal/nbinom on the R
  route and both bridge routes 10-20; F4's internal-route derivation and its
  `family =` formal 8-15; item 24's staleness and bitwise arms on
  `$setForestBasis` 14-24; item 25's anchor gate 18-30. Total ~155-275,
  midpoint ~215, which is where ~155-260 comes from. The last two lines are
  the two gates the earlier draft left implicit inside its 40-70 builds line
  and its "items 21-22 `nbar` behavior 10-18" line; both are now named
  checklist items with their own red and green cases, and their old carriers
  were reduced rather than left standing. Arm E's new arm in
  `.claude/m4-basis-design/harness/fa5-latent-coupling.R` (490 lines today) is
  plausibly 60-120 more and is excluded from the number but NOT from the work.
  **If the band is the binding constraint, the split to take is the residuals
  fix**, which has its own independent acceptance criterion at 15-25 dense and
  no dependence on the family work.

  **The docs band is COUNTED BOTTOM-UP from section 7's named deliverables**,
  which is what the two estimates before it were not. `~85-135 / stop 202` was
  sized on a smaller inventory; `~160-250 / stop 375` grew it by the nine
  gaussian-only sites, the five "sd(y) units" sites, the feature-matrix cells,
  M4.4's own landing note and 2f's four contract sentences, but did it as an
  estimate and in RAW lines, which is the currency error M4.3's landing note
  records ("most of the raw is Doxygen"). Every row of this table is
  DENSE-EQUIVALENT: wrapped Doxygen and markdown count at about half their raw
  lines, Rd at raw. By site: engine, bridge and header Doxygen **37-64** (the
  price table's own ~50-70 + ~10-15 + ~3-8 RAW comment lines, which is 32-47
  dense, plus `latentScaleAnchor`'s and `basisRowNorm`'s new blocks, item 6's
  `chain.hpp:879-884` rewrite, F6's `setForestWeights` sentence, F7's narrowing
  at `combiner.hpp:800-802`, and `combiner.hpp:953-955` RE-ARGUED rather than
  swept); Rd **14-25** (`man/forest.Rd:52-59`, `:71`, `:93`, `man/dbarts.Rd:84`,
  `man/dbartsSampler-class.Rd:147`, `:152`, `:351`); `multiplier-combiner.md`
  **26-41** (scope `:14-17`, the map section `:321-341` with the anchor table,
  2d(a)/(b) and the four contract sentences, the `:25` sd(y) phrasing, the
  Status discharge `:451-455`, and M4.4's own landing note - M4.3's is 20 raw
  lines in that file's per-slice format); `nameable-calibration.md` **3-7**;
  `feature-matrix.md` **13-31**. **Total ~93-168, midpoint ~130**, and the band
  is `~100-175 / stop 265`. EXCLUDED and conditional, on the same convention
  the tests row uses for arm E's harness: the MANIFEST and baselines-ledger
  obligation, 0-14 dense, which fires only if the trio moves.

  **10. Explicitly NOT M4.4's.** Each recorded so it is not rediscovered as a
  surprise.
  - The K-forest binary prior default (`aPriorScale = 2`). Under Option L, 38%
    of the K = 2 probit prior on p still lies outside (0.01, 0.99) and 16%
    within 1e-9 of a boundary, and the source is the PROGNOSTIC channel alone,
    not the anchor: `a * mu` has median |.| 1.175 but 95th percentile 20.17.
    That default was set for a gaussian response where a DRAWN sigma absorbs
    the excess; under probit and logistic sigma is pinned and it cannot.
    **M4.4 must NOT move it**, and must DOCUMENT that the shipped binary
    K-forest prior is more diffuse than the shipped single-forest one
    (`P(p < 0.01 or p > 0.99)` = 0.376 versus 0.238) and point at the deferred
    slice. VD decision 2026-08-14: its own slice, scheduled after M4.4's arm E
    and before 1.0-0; root TODO `binary-kforest-prior-default`. It does not
    disturb the L-versus-C verdict - the ratio stays exactly 3 regardless, and
    the decisive all-basis arm does not involve `aPriorScale` at all. **Two
    more items join that same slice, both from 2f, and M4.4 must not do
    either:** the `sqrt(K)` dispersion law (no basis-side option touches the
    exponent, and any K fix is a DEFAULTS change; M4.4 owes the exponent in
    prose only), and the OBSERVABILITY gap (option F - `v_f` is readable from
    no getter, so a defaulting user cannot reconstruct their induced prior).
    The TODO entry names all three.
  - `bart2()` has no `forests` argument (grep-verified), so a binary K-forest
    fit is reachable only through `dbarts()` / `dbartsSpec()`. Widening
    `bart2` is a separate, larger surface.
  - `dbarts_sampler_setWeights` (`src/C_interface.cpp:483-490`) does not call
    `refuseBinaryWeightChange`, so the flat layer permits a post-creation
    weight swap on a binary sampler that the R bridge refuses, unenforced
    against `dbarts.h:560`. Pre-existing; +1-2 flat-C dense if taken.
  - The divergence F3 records: `createBCFHolder` does not run the offender
    cascade `createHolder` runs.
  - The small-`omega` `PG(1, psi)` tail bound F7 scopes.
  - The wider doors - aft, ordinal, nbinom, grouped, multinomial - stay
    refused by program Fork 2, and the refusal messages should name what each
    is missing rather than restating the refusal. aft is the CHEAPEST door
    (its map is the gaussian one verbatim, since `AFTResponse` delegates
    `fitScale()` to an internal `GaussianResponse` and its working response is
    the range-scaled log time; it needs the survival-status channel threaded
    to the K-forest creation path and `sigmaIsFixed_` handled, since aft
    DRAWS sigma). ordinal reuses probit's map exactly; nbinom reuses
    logistic's; grouped would pick up the group draws through an empirical
    anchor - a second instance of item 2c's defect - and needs its tau block
    shown to interleave with the amplitude block; multinomial has its own
    combiner and its own anchor and is structurally a different design, not a
    door on this one. Say that last one rather than leaving it implied.

  **11. Discipline for the implementer.** This slice edits files AND cites them
  by line number, which is what cost the previous slice two review rounds.
  Re-derive EVERY anchor by OPENING THE TARGET and locating the SYMBOL; never
  apply an arithmetic offset, not even when everything around it moved by a
  known amount. Run the anchor pass LAST, after every content edit is final.
  **The failure this rule exists to prevent has already happened once on this
  bullet, so the signature is known:** an earlier pass moved all seven
  `docs/design/multiplier-combiner.md` -> plan citations by EXACTLY the plan's
  net delta (+731) and reported them as re-derived. All seven landed correct
  by coincidence, and FIVE more citations above the edit region were left
  untouched and became wrong - one of them (`:319`) in the SAME SENTENCE as a
  citation the same pass did update. **A per-citation delta table in which
  every delta is identical is the signature of an offset pass**, and it is
  what a reviewer should ask for. `multiplier-combiner.md` carries EIGHTEEN
  such citations on SEVENTEEN lines (`:35`, `:63`, `:132`, `:134`, `:203`,
  `:233`, `:262`, `:263`, `:285`, `:286`, `:303`, `:318`, `:319`, `:348` -
  which carries TWO - `:382`, `:420`, `:445`); M4.4 edits that file, so every
  one of them must be re-derived by CONTENT after the edits land. **The
  "thirteen" this bullet used to state omitted `:262`, `:285`, `:348` and
  `:382`**, which is the same defect one level up: an inventory of anchors is
  itself an anchor, and is re-derived rather than inherited.
  And keep plan and process references OUT of code comments - a `docs/design`
  citation in a Doxygen block is fine, a `docs/plans` slice reference is not;
  code comments carry self-contained rationale.
- **M4.5 (docs). Now the LAST slice of the Gaussian-complete arc** (re-scoped
  2026-08-13 on the probe verdicts), landing BEFORE M4.4 rather than after it,
  so it documents a K-forest Gaussian family and sweeps FIVE of the six
  `dbarts.h` Doxygen paragraphs in the `window:` note - `:505` ("Gaussian
  responses only.") is still true at this point and rides M4.4. **That
  sentence is now at `dbarts.h:513`**, shifted by its own paragraph's growth
  when M4.5 swept the other five; M4.4 must re-derive it rather than take
  either number on faith.
  Deliverable status verified 2026-08-13:
  `docs/design/multiplier-combiner.md` **does not exist - NEW FILE** (its only
  mention anywhere is this line); `docs/design/model-space-survey.md` exists,
  so "gains the four verified classes it lacks" is an EDIT;
  `docs/design/forest-combiner.md`'s "What still re-carves when a second
  combiner lands" is at `:146-189`, an EDIT. Plus the SIX meaning-moving
  `dbarts.h` Doxygen paragraphs listed in the `window:` note above (comment
  edits only - no re-bake, no version-constant move). `docs/design/bcf.md:
  325-329` was factually wrong (it named `setTreatment`/`bcfGlue` as public
  `dbarts.h` entries and cited `dbarts.h:264-271`, all retired at reshape S1);
  **FIXED 2026-08-13 in the pre-M4 amendment commit**, so M4.5 inherits a
  correct paragraph rather than a deliverable.

  M4.3 LANDED 9c63e9d8 without touching either open doc item, so M4.5 now
  inherits both by name: the `docs/design/bcf.md` a-move-as-prognostic-only
  correction, stale since M4.2 recorded it and deferred it here, and the SIX
  meaning-moving `dbarts.h` Doxygen paragraphs listed above, still open
  after M4.3. M4.5 is the docs sweep and closes both.

Falsifiers, all AMENDED 2026-08-13 except FA3, FA4 and FA6.

**FA0** a K = 2 factor-basis sampler is bitwise identical to today's BCF at the
same seed (this REPLACES the memo's F2, which fact 9 falsifies). Its CREATION
half is DISCHARGED: M2's FS1 pinned the `forests =` route bitwise on all six
`bcf-equivalence` channels, because that route resolves to the SAME payload
`treatment =` did (`resolveForests` `R/model.R:759`, `forestParams` `:840-854`).
FA0's live content is now exactly that the post-M4.1 GENERALIZED combiner still
reproduces it.

**FA1, the amplitude falsifier** (the critique's A3, adopted), RESPECIFIED as
a binary-basis proxy. **CORRECTED 2026-08-14: the stated reason was FALSE and
is struck.** This plan claimed that "a strict VCBART DGP wants a CONTINUOUS
basis column, and today's engine refuses one by name in two places
(`R/model.R:693` with its message at `:695-696`; `C_interface.cpp:776-789`)".
Neither site refuses one: `R/model.R:694` is the logical-to-factor coercion
and `:697-698` is `stop("a 'basis' cannot be NA")`, and
`C_interface.cpp:777-780` is the finiteness scan. A continuous basis column is
ACCEPTED on every route and RUNS to completion, MEASURED against a private-lib
build (see 2f, where the same fact is what makes the basis-scale trap
reachable). FA1 stays a binary-basis proxy because that is the shape the
shipped `bcf` fixtures and the `forest(basis = ~ factor(z))` route give for
free, not because the alternative is refused. Two arms, both of which run
TODAY on shipped surface.
- **FA1a, the K = 2 binary-basis proxy (engine).** Both arms
  `dbarts(y ~ ., forests = list(forest(...), forest(basis = ~ factor(z), ...)))`.
  WITH: defaults - `a` free (half-Cauchy via `aVariance`), `b0`/`b1` free, the
  ASIS ridge move active on forest 0. WITHOUT: `forest(update.amplitude =
  FALSE)` on BOTH forests, which pins `a = 1, b0 = 0, b1 = 1` at their
  constructed values (`combiner.hpp:319`), makes `drawGlue` consume NO rng
  (both blocks gated at `:589` and `:609`) and makes `afterCombine` an
  immediate no-op (`:640`). The model is then `y = mu(x) + z tau(x) + eps` with
  FIXED leaf-scale hyperparameters - VCBART's parameterization exactly, at
  K = 2, with no engine change. **The confound and its handling:** forest 0's
  node scale is `s = scaledResponseSd` unconditionally (`chain.hpp:718`) and
  there is no knob for it on a BCF sampler (`$setCalibration` refuses,
  `R/dbarts.R:1484-1489`), so the WITHOUT arm's prior on the mu total is
  `N(0, s^2)` while the WITH arm's is `a` times that with
  `a ~ C+(0, aPriorScale = 2)`. Forest 1's scale IS tunable via `forest(sd =)`
  -> `sdModerate` (`R/model.R:849`, `R_interface_bartcore.cpp:2145`), so run
  the WITHOUT arm over a small grid of `forest(sd =)` values and report the
  arm's BEST cell, with forest 0's fixed scale stated as a recorded
  limitation. A no-amplitude arm that wins from its best cell is a strong null;
  one that loses only at a badly chosen scale proves nothing.
  **DGP:** Friedman (1991) for mu (the house convention,
  `benchmarks/R/grouped-mixing.R:56-58`), a tau surface breaking at DIFFERENT
  modifier thresholds than mu (FA4's spirit, so the two forests are not
  near-collinear), binary z with ~40% treated, n = 1000, sigma = 1.
  **Functionals, all identified:** `sigma`; `mean_i(tauhat_i)` (the WITH arm's
  `(b1 - b0) * mean tau`, the WITHOUT arm's `mean tau`); `||mu_total||_2`;
  fitted values at 5 fixed held-out rows. **Reported:** IACT and out-of-sample
  RMSE against the true `mu` and `(b1-b0)tau` surfaces, averaged over seeds
  with the spread printed. **IACT estimator, named:** `coda::effectiveSize`
  (AR-spectral), IACT = N/ESS, per `grouped-mixing.R:68-77`; `coda` is NOT in
  `DESCRIPTION` Suggests, so the harness `requireNamespace`-guards it and stops
  with a message as `grouped-mixing.R:39-41` does (the only ESS estimator with
  a declared dependency is `posterior::ess_bulk`, `R/diagnostics.R:95,100,127`
  - a promoted harness is the moment to consider switching). Report lag-1
  autocorrelation beside IACT; single-run IACT is recorded as unstable in two
  places, which is why seed replication is not optional. **Shape:** K = 2;
  n = 1000 train + 1000 test; 1 chain per fit (single-chain ESS, the house
  convention); 2000 burn + 6000 kept, thin 1; 8 seeds per arm per grid cell;
  3 grid cells for the WITHOUT arm. 32 fits, ~3 minutes of engine time at the
  0.47 ms/sweep the current baseline records.
  **RESULT, 2026-08-13** (`probes-2026-08-13.md`; ran as specced, 32 fits in
  298 s, 9.2-10.0 s per fit; base seed 20260813, both arms seeing IDENTICAL
  data at each replicate; the response-scale decomposition identity of
  `test-bcf-reporting.R:57-71` holds to 7.1e-15 on all 32 fits). **The WITHOUT
  arm does NOT match.**
  - **IACT**, each functional allowed its OWN best WITHOUT cell, mean over 8
    seeds: `sigma` 65.5 -> 107.7 (**1.65x**), `mean tauhat` 26.5 -> 40.7
    (**1.54x**), the five held-out fitted values 58.4-66.4 -> 71.7-100.4
    (**1.08x-1.72x**). Seven of the eight functionals are WORSE, and paired per
    seed WITHOUT wins at most 4 of 8 on any of them. Lag-1 autocorrelation
    orders the arms the same way throughout, so no verdict rests on the
    spectral estimator alone. The one functional WITHOUT wins is
    `||mu_total||_2` (0.08x), and it is not like-for-like: in the WITH arm that
    IS the amplitude `a`'s own mixing (IACT 235.7 - the bottleneck the ASIS
    ridge move exists to relieve, `bcf-ridge-interweaving.md:494-499`), while
    in the WITHOUT arm `a == 1` identically and there is no amplitude to mix.
  - **RMSE** on the held-out rows: muRMSE **0.6027** (WITH) against 0.6836 /
    0.7038 / 0.6934 (sd 0.5 / 1 / 2), i.e. 12-17% worse in EVERY cell, with
    WITHOUT winning **0 of 8** seeds - including when the best cell is chosen
    per seed after the fact (0.6750, ratio 1.120, still 0/8). tauRMSE TIES:
    0.3687 (WITH) against 0.3685 at the best cell, and 0.969 winning 6 of 8
    under per-seed post-hoc selection.
  - **The two constructibility caveats, recorded beside the numbers they
    qualify.** (i) *The muRMSE confound the spec anticipated is REAL and is
    carried wherever muRMSE is cited here:* forest 0's node scale is
    `s = scaledResponseSd` unconditionally (`chain.hpp:718`) and `forest(sd =)`
    on the FIRST forest maps to `spec.aPriorScale`
    (`R_interface_bartcore.cpp:2145`), not to a node scale, so the WITHOUT arm
    cannot retune the mu prior at ALL. It does NOT reach the IACT gap - a
    looser prior mixes worse, not better - and the tau tie runs the same
    direction. (ii) *The out-of-sample mechanism:* a BCF sampler REFUSES test
    predictors (`R/spec.R:455`), so the held-out half was carried as 1000
    additional design rows with WEIGHT 0. That is exact, not a workaround -
    zero-weight rows drop from the weighted SSR, every leaf conditional and the
    posterior df (`model.hpp:2481-2483`) while still being assigned to leaves
    and receiving fitted values (`moves.hpp:69`) - so the likelihood saw
    exactly 1000 rows in both arms, and the run is an incidental validation of
    the zero-weight machinery.
- **FA1b, the amplitude-free family at K > 2 (pure R composition).** FA1a
  cannot say whether the amplitude matters AT SCALE, which is the VCBART
  regime. Build the amplitude-free K-forest sampler as K single-forest gaussian
  dbarts samplers composed in R (exactly correct for gaussian - this plan's own
  fact 5), leaf priors matched via `$setCalibration`, host owning the residual
  bookkeeping and the sigma draw. Measure IACT and RMSE at K in {2, 4, 8} on a
  VCBART-shaped DGP (`y = sum_j beta_j(u) X_j + eps`, continuous `X_j`, legal
  here because the HOST applies the basis). Report the composition tax at each
  K; that number is also an input to FA5's teeth.
  **RESULT, 2026-08-13** (`probes-2026-08-13.md`; 24 fits, 8 seeds, K in
  {2, 4, 8}, n = 1000 fitted + 1000 test, 50 trees per forest, the beta library
  and X held FIXED across K within a seed so beta_1 and beta_2 are the same
  estimand at every K). **Degradation is PRESENT, and it is small.**
  - **IACT grows:** per-forest mean 1.75 -> 2.58 -> 4.94 across K = 2, 4, 8
    (**2.8x**) and `sigma` 3.7 -> 5.4 -> 11.2 (**3.0x**), monotone and present
    in every seed. That IS K-scaling degradation, and it is why this conjunct
    fails as written.
  - **It grows off a low base and does not reach accuracy:** every per-forest
    IACT at K = 8 is under 6; predictive RMSE is flat (beta1 +3%, beta2 +6%
    from K = 2 to K = 8 despite the DGP carrying four times as many terms); and
    the combined-location IACT is flat (11.2 -> 10.0 on the four non-outlier
    test rows - `loc4` at K = 8 is a single-seed outlier, max 243.1 against a
    median near 10, and is reported separately rather than folded into a mean).
  - **Composition tax, with its denominator stated:** **5.14x at K = 2, 5.36x
    at K = 4, 5.43x at K = 8** - per-sweep wall time of the K-sampler R
    composition over a SINGLE batched engine sampler carrying the same total
    tree budget (K * 50) on the same n, back to back on the same box. Growth
    from K = 2 to K = 8 is **1.06x**, so the tax is FLAT in K at a fixed total
    tree budget, not "growing with K". This is the figure that replaces the
    struck 1.39-1.43x ("Departures" item 3), and it is FA5's teeth input.
- **Verdict rule, PRE-REGISTERED.** If FA1a's WITHOUT arm matches on IACT and
  RMSE from its best scale cell AND FA1b shows no K-scaling degradation, the
  amplitude and its ASIS remedy are BCF-SPECIFIC, not family content.
  Consequence: M4.2 keeps BCF's specialized two-scalar path and the general
  path ships amplitude-FREE by default, amplitudes opt-in. That re-scope is
  worth roughly the whole of M4.2's engine budget, which is why FA1 runs first.
- **VERDICT, 2026-08-13: the rule is a CONJUNCTION and BOTH conjuncts FAIL.**
  FA1a's WITHOUT arm does not match (1.08x-1.72x worse IACT on 7 of 8
  functionals from its own best cell; 12-17% worse muRMSE in every cell on 0 of
  8 seeds), and FA1b measures 2.8x-3.0x IACT growth from K = 2 to K = 8. **The
  amplitude and its ASIS remedy are NOT shown to be BCF-specific, the re-scope
  is NOT licensed, and M4.2 keeps full scope** - the general path does not ship
  amplitude-free, amplitudes are not demoted to opt-in, and the budget saving
  is unavailable. FA1 is NOT promoted to `benchmarks/R/`, its promotion being
  gated on exactly the re-scope that did not happen.

**FA2** with the per-forest ASIS rescale removed, the amplitude's IACT degrades
measurably; if it cannot go red, the remedy is decoration. **Its
single-amplitude half is ALREADY SATISFIED and the measurement is in the
repo:** `docs/plans/bcf-ridge-interweaving.md:494-499`, from the arc that BUILT
the a-move - initial-monotone IACT, mean over 4 seeds, strong prognostic
signal, n = 200 - `|a|` IACT **321 -> 130 (2.5x)** (base spread 198-544, new
100-156) and sigma IACT **64 -> 37 (1.7x)** (base 17-112, new 28-50); with the
sharper control at `:488-492`, where refreshing `aVariance` after the move
(breaking the ASIS conditioning) HURT `|a|` mixing 69 -> 196. So the remedy is
not decoration and the mechanism is specifically the ASIS conditioning. **The
pre-arc re-run is RETIRED at explicit orchestrator discretion** - this is a
departure from this section's own "FA1 and FA2 still run BEFORE any engine
work" reading and is recorded as one rather than carried silently; the cost
saved is small (one `--preclean` build plus minutes), and what it buys is a
number the repo already has. **What genuinely remains is the q-VARIATE half,
and it is an M4.2-INTERNAL gate, not a pre-arc probe:** run it with the
per-forest move in and out, on the same two-build method the a-move used. Its
acceptance question is its OWN IACT payoff (`|b1-b0|` or tau-amplitude IACT on
a strong-treatment-signal DGP) plus `bcf-exact` mode-2b staying exact and a
keepTrees BCF round-trip tracking, per `bcf-b-ridge.md:438-449`. **The
`bcf-sigma-residual` acceptance question is STRUCK** - that flag is RESOLVED
and the move was exonerated as its carrier (see M4.2).

**FA3** a continuous single-column basis recovers a known closed-form posterior
to Monte Carlo error, in `bcf-exact.R`'s idiom - whose header `:1-20` states
the enumeration plus quadrature construction and its three modes, and whose
mode 2b comment (`:35-38`) independently names "the bscale ridge".

**FA4** a K = p+1 VCBART-shaped sampler recovers known coefficient surfaces
better than one `LinearGaussianLeaf` fit on a DGP where the coefficients break
at DIFFERENT modifier thresholds; a null means this should have been a leaf
model. Bounded at p <= 8 on the leaf side: `LinearGaussianLeaf`
(`model.hpp:918`) has `maxNumCovariates = 8` (`:924`).

**FA5, the non-Gaussian claim, RESPECIFIED.** The committed form has no teeth:
its comparison arm is a composition no competent host would write (K probit
samplers each drawing latents against its OWN partial fit), so a positive
result licenses nothing, and its engine arm does not exist until M4.4, so the
falsifier cannot run before the work it is supposed to justify. The respec is
TWO legs, split on where an exact posterior is constructible at all.

- **Leg G (Gaussian, EXACT, runs pre-arc).** `bcf-exact.R`'s exactness rests
  entirely on Gaussian conjugacy: "conditional on the glue (a, b0, b1) and
  sigma the leaf parameters integrate in closed form (a Gaussian block
  marginal), leaving a low-dimensional quadrature over the glue and a 1-D
  quadrature over sigma" (`:7-10`). Under probit, observation i contributes
  `Phi(mu_{l(i)} + z_i tau_{m(i)})`, which COUPLES its mu leaf to its tau leaf;
  the leaf marginal becomes a joint integral of dimension L_mu + L_tau,
  blocked only over connected components of the leaf-overlap bipartite graph
  and NOT factorizable per leaf - up to 6 dimensions per enumerated tree pair
  at this DGP. So the exact arm is restricted to the GAUSSIAN legs, where it
  belongs. Leg G's job is to validate the R-composition DRIVER against the
  exact posterior on the same enumerable design (`bcf-exact.R:47-49`: K = 3
  cells, single-tree forests): two single-forest Gaussian samplers, host-owned
  residual bookkeeping, must recover the exact posterior. It also FIXES the
  tolerance unit used downstream, and it proves the loop is not itself
  perturbing the chain via `test-bcf-reporting.R:86-99`'s batch-vs-loop
  `expect_identical` template.
- **Leg P (probit, NO exact reference, two arms plus a control).** DGP: the
  same degenerate enumerable design on a probit likelihood,
  `y_i ~ Bernoulli(Phi(mu(x_i) + z_i tau(x_i)))`, amplitudes PINNED
  (`update.amplitude = FALSE` in every arm) so the probe isolates the LATENT
  coupling and nothing else - the amplitude question is FA1's.
  - **Arm B (the decisive R-composed arm).** Two single-forest GAUSSIAN dbarts
    samplers over the same design, each created with `resid.prior = fixed(1)`
    so sigma is pinned at 1 across sweeps (the value is a VARIANCE, not an sd;
    recorded footgun, `docs/design/correlated-outcomes.md:161-165`). Per sweep,
    in R: (1) read both forests' fits and form the combined location
    `F_i = mu_i + z_i tau_i`; (2) draw `w_i | y_i, F_i ~ TN(F_i, 1)` truncated
    to the sign of `y_i` - **against the COMBINED fit**, which is the whole
    point; (3) per forest, `$setResponse(w - other_f, updateScale = FALSE)`,
    `$setWeights(m_f^2)`, `$run(0L, 1L)`. Leaf priors matched exactly through
    `$setCalibration(forest, prior.scale =)` (`R/dbarts.R:1461-1527`; legal on
    a single-forest sampler, refused only on BCF). **The location lever must be
    matched too, and it is easy to miss:** `prior.mean` is the response
    transform's SHIFT and an offset of `-prior.mean` re-centers the modelled
    quantity (`test-calibration-creation.R:56-63`, applied as
    `setOffset(rep_len(-priorMean, rows), updateScale = FALSE)`). With two
    samplers each taking a per-forest RESIDUAL response the shift is
    per-sampler and no longer cancels; omitted, an arm-B disagreement is
    uninterpretable, because it is the very confound the DIFFERS row sends the
    implementer to chase. **Predicted: AGREES.** Arm B is cheap because the
    K = 1 version already ships as a test:
    `inst/tinytest/test-calibration-creation.R:53-80` is this loop at K = 1.
  - **Arm A (the committed strawman, kept as the POSITIVE CONTROL).** Two
    independent PROBIT dbarts samplers, each drawing its own latents against
    its own fit, composed additively in R. **Predicted: DIFFERS.**
  - **Arm E (the engine), POST-M4.4.** The K-forest probit sampler, run at
    M4.4 as its acceptance gate. The decisive comparison for M4.4's own
    acceptance is arm B vs arm E.
  - **Reference for leg P:** a long INDEPENDENT MCMC written directly against
    the enumerated probit model with no dbarts sampler in the loop. It is
    APPROXIMATE, so it carries its own convergence gate (>= 4 chains, R-hat <=
    1.01 on every reported functional) and its own MC-error budget; nothing in
    leg P is called "exact".
- **Measurement and the PRE-REGISTERED decision rule.** Functionals:
  `E[mu(x)]` and `E[tau(x)]` at each of the 3 cells and `P(y=1 | x, z)` at each
  of the 6 (x, z) cells - 12 per arm, 24 z-statistics across arms A and B.
  Per arm per functional: posterior mean and its Monte Carlo standard error
  over >= 8 independent seeds, the z-statistic against the reference
  (`z = (arm - ref) / sqrt(se_arm^2 + se_ref^2)`), plus IACT
  (`coda::effectiveSize`, IACT = N/ESS) and per-sweep wall time (reported, not
  gated). **Threshold:** a functional AGREES at `|z| <= 3`. **Multiplicity:**
  at 24 comparisons an isolated `|z|` in (3, 4] is expected (~0.065 spurious
  flags at 3, ~1.5e-3 at 4), so an arm is declared to DIFFER only if some
  `|z| > 4` OR two or more functionals exceed `|z| = 3`; a single flag in
  (3, 4] with no second is AMBIGUOUS and is RESOLVED BY RE-RUNNING at 4x the
  seeds, never adjudicated after the fact. **Power precondition, evaluated
  BEFORE the outcome table is read:** arm A, the positive control, must DIFFER
  from the reference at `|z| > 6` on at least one functional. If it does not,
  the DGP has too little cross-forest signal to detect anything and **the run
  is VOID** - re-run with stronger cross-forest signal. This is a precondition,
  not a verdict; it is what keeps a weak-DGP artifact from being read as a
  null.
- **Outcomes, given the power precondition holds.** Two, and only two.
  - **Arm B AGREES on all 12.** The composition IS buildable in R by a
    competent host; the naive one is not. **M4.4's "the caller cannot compose
    non-Gaussian latents against the combined fit" ground FALLS.** M4.4 is
    re-justified on the residual, which is nameable: the per-forest ASIS
    rescale (it writes leaf values, and no setter exists or should), the
    `0x1p-26` snap's exactness, the K-sampler composition tax (MEASURED by
    FA1b: 5.14x at K = 2 and 5.43x at K = 8 against a single batched engine
    sampler at matched total trees, flat in K at 1.06x growth - the 1.39-1.43x
    written here is struck, "Departures" item 3 - plus an O(nK) R-level latent
    draw per sweep), and ergonomics. That is a WEAKER justification than this
    plan asserts and must be written down as such. Predicted outcome.
  - **Arm B DIFFERS.** The composition CANNOT be built correctly from the
    shipped R surface. **M4.4 is JUSTIFIED, and the probe hands the
    implementer the reason.** Mandatory follow-up: identify WHICH
    shipped-surface gap causes it - sigma pinning across sweeps, the response
    transform's location lever, the ordering of the amplitude draw relative to
    the latent draw, or the near-zero multiplier snap. If the gap is a missing
    R knob costing ~20 lines, M4.4's justification shrinks to that knob and
    the non-Gaussian half re-prices accordingly.
- **RESULT, 2026-08-13** (`probes-2026-08-13.md`; dataset seed 20260813, ONE
  dataset with only the MCMC seed varying, so the across-replicate spread is
  pure Monte Carlo error; amplitudes pinned in every arm). **The first outcome
  obtained.**
  - **Leg G: PASS, max |z| = 1.35.** The R-composition driver recovers the
    EXACT posterior (5 x 5 enumerated tree pairs, leaves integrated in closed
    form, sigma fixed and amplitudes pinned, so no Monte Carlo enters the
    reference). G0's batch-vs-loop `expect_identical` template holds - the
    driver loop does not perturb the chain - and the leg-P reference validated
    against the same exact posterior on its Gaussian branch at max |z| = 1.16.
    The tolerance unit this fixes is an MC standard error of 7e-4 to 2e-3 at 8
    seeds x 6000 kept.
  - **Power precondition: MET.** Arm A, the strawman (each forest drawing its
    own latents against its OWN contribution - a strictly sharper control than
    two independent probit samplers, and recorded as that departure), DIFFERS
    from the reference at max |z| = **772.49**, with 12 of 12 functionals over
    |z| = 4. The DGP carries ample cross-forest signal; the run is not void and
    the outcome table may be read. The leg-P reference itself passes its
    convergence gate at max R-hat 1.0003 against 1.01.
  - **Arm B AGREES on all 12**, max |z| = **2.54** against the pre-registered
    threshold of 3.0, zero functionals over 3, zero over 4 - so no AMBIGUOUS
    flag and no 4x-seed re-run. Its IACT tracks the reference's within about
    1.5x on every functional. **No arm-B functional required a knob the package
    does not ship**, and the location lever was HANDLED rather than assumed
    (centered construction vector, `prior.mean == 0` and `prior.sd == leafSd`
    ASSERTED on every sampler built).
- **The teeth, stated plainly.** M4.4's headline justification - "the couplings
  a caller cannot compose" - is FALSIFIED if the power precondition holds and
  arm B AGREES on all 12 functionals. **Fork 2's already-live Gaussian-first
  escape hatch is LICENSED, and should be taken, if all three hold: arm B
  agrees; FA1b's measured composition tax at K = 8 is under 2x; and no arm-B
  functional required a knob the package does not ship.** Only a measured
  arm-B failure preserves the claim as this plan writes it. **Buildability is
  RESOLVED, which is what sharpens the falsifier:** the one thing that could
  have made arm B unbuildable - holding sigma across sweeps - is a PUBLIC knob
  and the supported idiom. `resid.prior = fixed(v)` sets `model.sigmaIsFixed`
  (`R_interface_bartcore.cpp:1381-1387`), gating the end-of-sweep redraw at
  `chain.hpp:1276-1277`; `$setSigma()` ALONE does not suffice, because it
  installs unconditionally (`chain.hpp:1387`) and an unfixed gaussian sampler
  overwrites it every sweep. The bridge's own comment (`:2655-2680`) names the
  combination as intended and `test-sampler-errors.R:159-168` asserts it is
  permitted. `fixed()` is reachable publicly as `dbartsPriors$fixed(v)` or by
  bare name inside `resid.prior =`. So an arm-B failure cannot be blamed on a
  missing knob.
  **BITTEN, 2026-08-13.** The precondition held and arm B agreed on all 12, so
  M4.4's headline justification is FALSIFIED and is struck from that slice. The
  hatch's three conditions: (1) arm B agrees - YES, max |z| = 2.54; (2) no
  arm-B functional needed an unshipped knob - YES, enumerated in the result
  above; (3) FA1b's measured tax at K = 8 under 2x - **the clause's 2x anchor
  was calibrated against the 1.39-1.43x figure that has NO receipt, so the
  anchor falls with the figure** (erratum, "Departures" item 3), and the
  measured 5.43x is read as a REASON to build M4.4 rather than a bar to
  reordering it. **The hatch is TAKEN, by the orchestrator under fork 2's
  pre-authorization:** M4.0 -> M4.1 -> M4.2 -> M4.3 -> M4.5, then M4.4
  IMMEDIATELY, both slices pre-release. The hatch REORDERS; it does not cancel.
  Arm E is unchanged and is still M4.4's own acceptance gate, arm B against
  arm E.

**FA6** every refused creation raises an R condition and no external pointer
escapes. Already DISCHARGED as discipline: M2's FS2 landed 42
message-asserting refusals.

**Budget, RE-PRICED 2026-08-13 pre-M4**, superseding both the ~450/200/250
below and the ~900-1100 fork 1 option (a) carries. **~1850-2400 total:** engine
~500-700, bridge ~350-450, R ~400-500, tests ~350-450, docs ~250-300. What
moved and why. ENGINE grows only ~50 over the committed ~450-650, for the
four-layer `bcfGlue` vtable re-sign and the base virtual's return contract
re-stated over BOTH `afterCombine` implementations; **M4.2's own share moves
DOWN, not up**, because `bcf-b-ridge.md` already delivers the `p = (k-d)/2`
exponent, a passed adversarial prototype with a sharp discrimination arm, and
the complete rescale-consistency checklist. BRIDGE grows from ~200 for the
ragged params transport at two sites, the 3-row `glue` result slot, the
length-4 `bcf` state block re-encoded in place, the `amplitudes[3]` buffers,
and the treatment-slot retirement with its capability-probe replacement. R
grows from ~250 for `forestParams`' ragged transport, `forestBasisDeclaration`
and `expandForestBasis` going per-forest, `resolveForests`' six positional
refusals, `validateTreatment`, `$getForestAmplitudes()`'s broken signature and
`isBCFSampler`'s replacement. TESTS were priced only as "plus the M4.0 pins"
and are now priced: M4.0's `tests/cpp` pins over both `afterCombine`
overrides, plus the existing pins the relaxation rewrites. **M4.0 LANDED
562ee684 and consumed ~348 of the ~350-450 tests band (reviewer-counted); the
remainder is M4.3's rewrite share, which MUST be re-priced at its own
pre-slice check before M4.3 starts, not assumed from this figure.** **M4.1
LANDED 1458328c/e48fc5de and added ~98 tests; the band now reads ~446 of the
~350-450 tests band consumed (M4.0 ~348 + M4.1 ~98) - at or past the top of
the range, so the M4.3 re-price called for above is now MANDATORY, not a
contingency.** **M4.2 LANDED 1a2aaedc and added ~286+65 tests; the band now
reads ~797 dense-equivalent consumed (M4.0 ~348 + M4.1 ~98 + M4.2 ~351)
against the original ~350-450 band - well past the top of the range, so the
M4.3 RE-PRICE IS NOW A PRECONDITION: M4.3 does not start until its own test
budget is stated fresh.** (**~797, ERRATUM 2026-08-14**: this paragraph read
~732, which does not sum - 348 + 98 + 351 = 797, and M4.2's own landing note
reads "tests ~351".) **The M4.3 re-price is DONE (amendment 2026-08-14, in the
slice): tests ~645 dense-equivalent, band ~590-660, as M4.3's OWN line item
rather than a share of the spent band; non-test ~970-1110, of which engine
~250-300.** That last figure puts the ENGINE band at ~478-528 of ~500-700
consumed after M4.3, leaving ~172-222; **M4.4's engine work is UNPRICED and
its pre-slice check must price it** - the arithmetic is recorded here as
arithmetic, not as a verdict on the band. DOCS (M4.5) were
unpriced; `multiplier-combiner.md` is a NEW file. Plus one possible
`bcf-equivalence` re-record, with its `equivalence.yaml` bump in the same
commit. **Both scheduling preconditions are MET:** multiforest-predictor-
mutation S0-S4 landed (so the window holds one re-record) and the dbarts.h
reshape landed (S1 ab3aa2fa, S2 1bf2e69c, arc closed 25a21d3b), pre-release as
VD resolved 2026-08-11 rather than after the freeze as this plan recommended.
**UNCHANGED by the probe verdicts (2026-08-13):** this pricing already carried
M4.2 in full, so FA1's failure to license the re-scope removes a saving that
was never booked; and FA5 reorders M4.4 after M4.5 without changing either
slice's line count, beyond M4.4 carrying the one `dbarts.h:505` Doxygen
paragraph M4.5 would otherwise have swept.
**M0 remains DEFERRED at orchestrator discretion** (`TODO:236-237`) and is NOT
an M4 precondition; after M4.5 is the cheapest place to land it, since M4 moves
the vocabulary M0 documents.

## Cross-plan amendments (apply verbatim)

**ALL EIGHT APPLIED at e2cc1de; this block is a LANDED RECORD, not forward
text** (marked 2026-08-13, pre-M4). The quoted text is left exactly as it was
applied, including its pre-re-bake `dbarts.h` line numbers; the drift is swept
once in the `amended:` note at the top of this file.

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
  measured ergonomics (5.14x at K = 2 and 5.43x at K = 8 against a batched
  engine at matched total trees, re-scoped 2026-08-13 from the struck
  1.39-1.43x, and a documented prior difference). Cost: the
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
(`chain.hpp:1266-1267`), so K independent R samplers are a DIFFERENT model, not
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
  says the amplitude is decoration for non-BCF instances. **The condition FAILED
  (2026-08-13): FA1 says the opposite** - pinning costs 1.08x-1.72x on IACT and
  12-17% on muRMSE at K = 2 - so this stays a door in its ALREADY-SHIPPED opt-in
  form (`forest(update.amplitude = FALSE)`), never a default.
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
   (B3, adopted): 0.89-0.95x is the engine's own per-sweep drive, and it keeps
   its receipt (`bcf-public-surface.md:85-93`).
   **ERRATUM, 2026-08-13 (post-probes), in this plan's own house style: the
   composition figure this item asserted - 1.39-1.43x at K = 2, "growing with
   K" - is STRUCK, in all EIGHT places it appeared** (this item, the answer
   section, the B3 block, corrected fact 6, the third-package paragraph, M0's
   vignette item, FA5's arm-B outcome bullet, and fork 1 option (c)). Grounds:
   **no receipt exists anywhere in the tree** - the grep finds the assertion at
   each of those sites and never a measurement, and the figure was adopted from
   the critique without one. The tool-verified-claims discipline this plan
   applies to the memo, the critique and the shipped code applies to ITS OWN
   numbers; a cost quoted with no denominator and no receipt is exactly what
   that discipline exists to catch, and it had already propagated into a
   scheduling clause (FA5's "under 2x" hatch condition), which is how an
   unsourced number does damage. **Replacement, with the protocol stated so the
   denominator travels with the figure:** per-sweep wall time of the K-sampler
   R composition over a SINGLE batched engine sampler carrying the same total
   tree budget (K * 50) on the same n, measured back to back on the same box -
   **5.14x at K = 2, 5.36x at K = 4, 5.43x at K = 8**, K-relative growth
   **1.06x**, so the tax is FLAT in K rather than growing with it
   (`probes-2026-08-13.md`, FA1b, 8 seeds, gitignored harness). A stricter
   denominator than whatever produced 1.39-1.43x, and stated rather than
   implied: any future re-quote must carry it or re-measure.
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

## Landing notes

M1 LANDED 05ac3b4b, 2026-08-13 (implemented as 78b080c5, amended during
independent review). Public R5 `$setForestWeights(forest, weights)` on
`dbartsSampler`, forest 1-BASED and converted at the boundary via
`resolveForestIndex` (the settled decision, item 1 above; a BCF sampler's
treatment forest is `2L`); a new `forestWeights` R5 field mirroring every
installed per-forest weight, re-applied at THREE re-creation sites -
`getPointer()`, `setState()`, and `$copy()` (the third added during review,
beyond the plan's two named sites); Rd on the class page (the three-way
all-ones contrast, the four-factor composition `(w_i a_i) m_f^2 s_i`, a
1-based `\item{forest}` clause, the rebuild hazard, the unreachable-`setData`
disposition); the existing NEWS bullet amended in place (the `dbarts:::`-only
caveat removed); new `inst/tinytest/test-forest-weights-r5.R` (12 results).
R-only, no `src/` touch, no baseline change. Clears reshape S1's HARD
precondition (`bcf-public-surface.md`'s M1-before-the-flat-entry rule).

Review: independent Opus reviewer, full battery from scratch, verdict
LAND-AFTER-CHANGES; after the fix pass, final verdict LAND, all five findings
resolved. The two blockers: (1) both FW1 mirror arms were VACUOUS - they
compared a weighted sampler against a DIFFERENT plain sampler, so the
stored-state difference satisfied the assertion even with the mirror
mechanism deleted (proven empirically); (2) the Rd's all-ones leg claimed a
`setForestWeights` all-ones vector is "bitwise distinct from carrying none" -
false; it INSTALLS (a round trip reports it) but is the BITWISE IDENTITY on
draws (plan item 3(b) above, `test-forest-weights.R:73`). Red-team finding,
fixed as a review-mandated scope extension beyond the plan's two named sites:
`$copy()` silently DROPPED installed per-forest weights (a fresh R5 object,
`setState` reapplies the dupe's empty list); fixed with an explicit field
carry plus reapply onto the dupe's pointer, alias-safe by copy-on-modify,
both shallow and deep. Also fixed: FW2's materiality guard was float noise
(`sd(tau0) = 1.5e-16` on a never-splitting zero-weight forest) ->
`abs(mean(tau0)) > 0.1` (measured 0.456).

Fix verification, record this discipline: per-site falsifiers - commenting
out each of the three `reapplyForestWeights` call sites in the installed
copy turned EXACTLY its own test arm red with the other eleven results
green (the fixer ran all three; the reviewer independently re-ran the
`setState` falsifier in its own scratch tree - its arm red, 11/12 green).
Two design facts the fixer discovered along the way, worth the record (they
shaped the FW1 arms as symmetric fresh-pointer pairs sharing one donor
state and differing only in the field): `setState` on an already-valid
pointer does NOT clear a prior real forest-weight install; and a
fresh-pointer `setState` restore never reproduces a live pointer's
RNG-stream continuation (deliberate, per `inst/common/stateContinuation.R`).

Budgets: R 48 + ~4 (the copy fix) non-comment, man 4 raw (the dense
single-paragraph convention), NEWS 9, test ~137 - comfortably inside the
amended ~270-290 dense-equivalent envelope.

Gates: implementer, fixer and reviewer batteries all green from scratch
(preclean private-lib installs; `tests/cpp` plain all-pass confirming no
`src` drift; full tinytest 4436/0, no snapshot regenerated; the trio
IDENTICAL - `equivalence-8b047f8b` 37/37 strict, `bcf-equivalence-8b047f8b`
12/12, `multinomial-equivalence-1027be5` 10/10, no max |z| anywhere; air
clean; lintr zero lints; `R CMD check` on a clean-copy tarball Status OK
zero E/W/N; pkgdown no problems; an independent `saveRDS`/`readRDS` round
trip). NO x86 leg: an R-only slice with no baseline change - the x86
standard applies to engine slices. CI six-green on the push.

M2 LANDED 64b13b98, 2026-08-13 (implemented as 639ceea9, amended during
independent review). `dbarts()`/`dbartsSpec()` gain `forests =
list(forest(...))`, replacing the removed `treatment =`/`moderators =`/
`treatmentForest =` (formals, `NAMESPACE` export, `_pkgdown.yml` entry, Rd);
`$setTreatment` becomes `$setForestBasis(forest, basis)`, 1-based via
`resolveForestIndex` on the same side as `setForestWeights` and the
calibration pair (a BCF sampler's treatment forest is `2L`); `$getBCFGlue`
becomes the argument-free `$getForestAmplitudes()`; a new exported `forest()`
constructor is its own Rd topic and `_pkgdown.yml` entry, carrying the
complete ten-knob map settled in the pre-M2 amendment - `basis`, `vars`,
`n.trees`, `base`, `power`, `sd`, `interactions`, `blocks`,
`amplitude.prior.variance` (legal only on a forest given `basis =`) and
`update.amplitude` (legal on every forest of a `forests =` spec, refused only
at `K = 1` where no amplitude channel exists at all - so `forest(sd = )`
beside `dbartsData(treatment = )` configures rather than refuses).
`resolveSamplerSpec` resolves `forests =` to exactly the `treatment =`
payload. `test-bcf-creation.R` is rewritten (80 results); seven more test
files migrate off the removed formals; benchmarks are untouched, since the
bcf harness drives the internal constructor route and `settingsList()`
carries no creation vocabulary. R-only, no `src/` touch, no baseline change.

Review: independent Opus reviewer, full battery from scratch, verdict
LAND-AFTER-CHANGES; after the fix pass, final verdict LAND. The blocker:
`forest(n.trees = )`/`base = `/`power = ` on the first forest silently
overrode an explicitly declared `control@n.trees` / tree-prior `base`/`power`
value, undocumented and unpinned; refusal was ruled not cleanly
implementable, since an explicitly passed default is indistinguishable from
the default. Resolved as DOCUMENTED PRECEDENCE - the `forest()` value
governs, as the more specific declaration, a sentence added to
`man/forest.Rd` - plus a compound disagreement pin proving the rule
knob by knob: flipping precedence on `n.trees`, `base`, or `power` alone each
breaks it.

Four implementer deviations ACCEPTED. (i) The refusal set grew from the
planned nine to a larger set by item 4's own by-name-refusal rule: variables
on forest 0, a second forest with no basis anywhere, a non-`forest()` list
element, an empty list, a two-sided basis formula, and per-knob range
renames each earned their own refusal and their own pin - spec-consistent,
not scope creep. (ii) `n.trees`/`base`/`power` on forest 0 are implemented as
RESTATING the fit's own `control`/tree-prior slots, since the engine carries
no separate payload slot for forest 0 the way it does for forest 1's;
ruled the right call over a silent drop or a tenth refusal. (iii) the K = 1
nuance in the amendment's `update.amplitude` legality rule. (iv)
`resolveModerators` gained an argument parameter so PUBLIC errors say
`'vars'` while the internal route keeps `'moderators'`.

Review round two added pins for facts the first pass verified but left
unpinned - the `hasBasis` escape (`sd` reaching `params[4]`,
`update.amplitude` reaching `params[7]`, a K = 2 no-basis spec still
honouring `n.trees`) and the positive `interactions`/`blocks` routes - and a
suspected `tau.blocks` pin, tested and empirically falsified, was left
standing as written.

Budget ACCEPTED: about 791 dense-equivalent lines, inside the amended
775-895 band - R landed at 246, six over its own sub-band, bought back by
the accepted deviations.

Gates: implementer, fixer and reviewer batteries all green from scratch
(preclean private-lib installs; `tests/cpp` plain all-pass confirming no
`src` drift; full tinytest 4475/0, no snapshot regenerated; the trio
IDENTICAL - `equivalence-8b047f8b` 37/37 strict, `bcf-equivalence-8b047f8b`
12/12, `multinomial-equivalence-1027be5` 10/10; air clean; lintr zero lints;
`R CMD check` on a clean-copy tarball Status OK zero E/W/N;
`pkgdown::check_pkgdown` no problems, `forest` added and `treatmentForest`
removed; FS1 positive bitwise on all six channels plus both negative arms;
FS2 42 message-asserting refusals; FS3 single-forest byte-neutrality). NO
x86 leg: R-only, no baseline change. CI six-green on the push.

Recorded for M3 (the reviewer's finding 4, src-owned, non-blocker): removed
vocabulary still leaks from bridge messages - "bartcore_setTreatment
requires a BCF sampler", "bartcore_getBCFGlue requires a BCF sampler", and a
length-mismatched basis under `forests =` still erroring "length of
'treatment' must equal length of 'y'". M3 (riding `dbarts-h-reshape` S1)
retires these alongside the flat-C renames.

M3 LANDED at dbarts-h-reshape S1, ab3aa2fa, 2026-08-13 (implemented as
3a977b6d, amended during independent review). `dbarts_sampler_setTreatment`
re-signed to `dbarts_sampler_setForestBasis(sampler, forest, basis,
numColumns)` (C_interface.cpp:759); `dbarts_sampler_bcfGlue` replaced by
`dbarts_sampler_numForestAmplitudes` (:806) and `dbarts_sampler_forestAmplitudes`
(:819); the creation Doxygen (dbarts.h:484-505) re-worded end to end in
engine vocabulary - `forests = list(forest(), forest(basis = ...))`,
`setForestBasis`, `forestFits`, `numForestAmplitudes`/`forestAmplitudes` -
with no `treatment`/`bcfGlue` token left in it. All three message leaks
recorded above are retired: `bartcore_setTreatment` now errors "a forest
basis requires a sampler whose forests carry amplitudes"
(R_interface_bartcore.cpp:3675-3676), `bartcore_getBCFGlue` now errors
"forest amplitudes require a sampler whose forests carry them" (:3774), and
`validateTreatment`'s length-mismatch message is parameterized by an
`argument` name that the `forests =` creation path now passes as `"basis"`
(R/spec.R:404-410, `data.R:648` reads `"length of '", argument, "' must
equal length of 'y'"`) - no `'treatment'` string reaches a `forests =` user
anywhere. Falsifiers FC1-FC3 all green: FC1 a flat-created sampler through
`setForestBasis` reproduces the R5 path bitwise at a shared seed; FC2 the
refusal matrix, outcome and message text both checked in C; FC3 the hash
moved and neither version constant did. Budget, review verdict and gates are
recorded at the dbarts-h-reshape S1 landing note, not repeated here - this
item priced ~40 header + ~80 C_interface + ~70 consumer.c + ~60 test-capi.R
of that slice's ~1310 re-priced total.

M4.0 LANDED 562ee684, 2026-08-13, CI six-green sanitizers included, amended
ONCE after independent review (LAND-AFTER-CHANGES -> applied). Tests-only:
tests/cpp (3 files), inst/tinytest (1), TODO trim. rng NEUTRAL confirmed by
two independent trio runs.

Pins: testBCFCombinerSeam (test_sampler.cpp:2818-3160) and
testMultinomialCombinerSeam covering the full amended mandate -
forestMultiplier (combiner.hpp:762, both call sites :552/:560), combinedFits
(:567) incl. the snap-band arm (snap pinned in the reparameterization, exact
multiplier in the blend), drawGlue (:582) order a/aVariance/b0/b1 + both
block skips, BCFForestCombiner::afterCombine (:638) GIG rescale c =
sqrt(GIG((L-1)/2, ...)) with the three reachable 1.0 skips inert on the rng
stream, MultinomialForestCombiner::afterCombine (:1152) additive shift +
single normal (:1180) + the returns-1.0-having-moved convention (its pin is
the SOLE guard on that convention). Every pin states the wrong wiring it
catches; the reviewer confirmed no over-constraint - each BCF pin is a
quantity M4.1/M4.2 must reproduce as the K = 2 / q = 1, d = 1 special case.

The review's REQUIRED fix, recorded as the slice's lesson: the multinomial
fixture originally set numTrees = 1 on all forests, collapsing every
per-tree factor m_k to 1 and leaving four sites (:1160, :1175, :1176, :1186)
undiscriminated - the reviewer MEASURED perTree = c staying green. Fixed
with unequal per-forest tree counts ({2, 1, 1}, forest 0's leaf sum split),
plus a fixture ran-at-all guard; perTree = c and return c both re-proven RED
under the new fixture. LESSON, stated as a standing rule: a pin fixture must
give every factor in the pinned expression a value that discriminates it -
unit values silently vacate pins.

Red/green demonstrations: implementer one per pinned method (b0/b1 swap,
cross-amplitude, draw reorder, gigP -> L/2, shift doubling); reviewer
independently reproduced three (gigP -> L/2; b0/b1 order swap - not
implementer-emphasized; return c, which failed EXACTLY one pin, proving
non-vacuity).

Advisories accepted, no action: combiner.hpp:674/:676 non-finite returns are
in the mandate's skip list but unpinned - unreachable from a fixture
(near-true: an overflowing M could reach them), recorded as a deviation; the
hook arm exercises interweaveGlueRidge only with updateA = false, so the
chain.hpp:1280-1281 sweep order rests on the pre-existing testBCFInterweave
(:2471) and the trio.

Gates: implementer AND reviewer batteries both clean from scratch (tinytest
4582/0 twice; trio 37-strict + 12 + 10 bitwise twice, zero max-|z| lines;
tests/cpp clean builds; reviewer ASAN/UBSAN leg zero reports; air exit 0;
lintr clean).

Carried items closed in this commit's slice: the three S2 test-hardening
items (intent verbatim per grow-from-root-categorical-scan.md:924-930 -
20-cell ran-at-all counter, belowCap tally asserted both ways,
all-zero-weights present-category fixture with its own zeroed > 0 guard)
and the TODO residue sentence trimmed; the mask-arc carried cosmetics (test
labels, dash rulers) done.

Budget: 348 dense-equivalent added (reviewer-counted, claim accurate).

Arc status: M1 and M2 landed; M3 LANDED (ab3aa2fa); M4.0 LANDED (562ee684).
Only M0 (the on-ramp docs/vignette slice) remains, deferred at orchestrator
discretion; M4 (the general basis family) is RESOLVED pre-release, its
sequence now M4.0 (LANDED) -> M4.1 -> M4.2 -> M4.3 -> M4.5 -> M4.4, scheduled
after the reshape, whose S1 - the arc's last header-touching slice - landed
at ab3aa2fa, and whose S2 (consumer rebuilds and docs, moving no header) has
now landed too - the reshape arc is complete
(docs/plans/dbarts-h-reshape.md).

M4.1 LANDED 1458328c, follow-up e48fc5de, both CI six-green sanitizers
included. Engine: forestMultiplier -> dot(a_f, B_f[i,]), combinedFits -> K
loop (general in K, no K=2 special case), BCF's two bases synthesized
internally from z at installTreatment; 15 glue_.a/b0/b1 accessor sites
converted; state serializeGlue/restoreGlue round-trip and the flat C
bcfGlue/setForestBasis unchanged. Engine +48 net dense-equivalent, tests +98.

The slice's defining event: the first implementer (killed mid-run by an API
error) had a 12/12 bcf-equivalence divergence; the successor kept its patch,
proved unmodified HEAD reproduces all 59 baseline scenarios on this box, and
root-caused the divergence to FMA CONTRACTION ASSOCIATION - HEAD's
single-expression a*mu[i] + b_z*tau[i] contracts forest 0's product into the
closing add, a forward K-loop accumulator pre-rounds it instead (1 ulp on
~30% of rows, chain contamination within ~40 sweeps). Fix: accumulate FROM
THE LAST FOREST DOWN, leaving forest 0's product as the closing bare
multiply - shape-guaranteed at any K. Verified with -ffp-contract=off and a
three-way std::fma probe. The combination's summation association is now a
LOAD-BEARING, baseline-gating property.

Two tests added, each the SOLE guard on its fact: testCombinedFitsAssociation
(20 of 64 rows discriminating; the M4.0 seam pin CANNOT see association -
its reference expression inherits the test compiler's own contraction; the
follow-up hardened the pin to classify fusionVisible rows and skip the
direction assertion on never-contracting targets, hard-failing only on
absorbing the last forest's product) and testForestBasisSynthesis (mid-life
setTreatment rebuilds the basis, amplitudes preserved).

Review verdict LAND, the association story independently reproduced end to
end (forward accumulation -> 12/12 bcf mismatch with both seam pins green;
HEAD's expression restored verbatim -> association pin passes, so the pin
fixes the property, not the spelling). Reviewer's trap note, recorded as
standing: a perturbed install WITHOUT --preclean reported bitwise identical
- the no-header-tracking trap; every re-check must --preclean. Pins
byte-identical at both revs (extracted and hashed). No stale-basis mutation
route exists (all z routes funnel through installTreatment; setData/
setWeights/setModel refused multi-forest; per-observation and transactional
setPredictor touch x only - the bcf-equivalence per_observation scenarios
are the evidence). Follow-up also corrected the setData refusal's mechanism
comment (post-M4.1 only drawGlue indexes borrowed z; combinedFits/
forestMultiplier read the install-time basis snapshot); no message pins were
affected (the user-facing string is untouched, verified against its two
pins).

BENCH (VD-granted quiet window, 2026-08-14): B/A per-sweep ratios - bcf
n=20000 1.0098 +/- 0.0051 (12/12 rounds slower), bcf n=2000 1.0105 +/- 0.0025
(12/12), supplementary bcf n=200000 1.0112 +/- 0.0053 (4/4), single-forest
control 1.0030 +/- 0.0123 (sign-split, neutral - methodology clean). The
ratio is FLAT across 100x in n: per-element compute/indirection in the
combiner, NOT bandwidth. VERDICT: the ~1% BCF-only cost is ACCEPTED as the
price of the general multiplier; M4.2, which owns this code next, MAY claw
it back but bitwise identity governs - no re-record for an optimization.
Protocol/provenance: .claude/m4-basis-design/bench-m41-2026-08-14.md (12
alternating rounds, loadavg 2.11-3.68, per-arm git-archive trees and
privlibs; note there is NO checked-in armab script - armab is a protocol
name, and bench-sampler.R carries no BCF scenarios, so the harness was
written in its idiom, gitignored).

M4.2 LANDED 1a2aaedc, CI six-green sanitizers included, amended once after
independent review (LAND-AFTER-CHANGES -> all four dispositions applied).
Engine: per-forest q-variate amplitude conditional (P = I/priorVar + sum w x
x'/sigma^2 on the basis-scaled fits, residual net of the other forests)
drawn through a NEW square-root-free unit-lower L D L' beside the shipped
Cholesky (the two-sqrt solve path gives x/sqrt(d)/sqrt(d) != x/d, breaking
the q=1 bitwise reduction); afterCombine is ONE general per-forest ASIS
rescale at p = (L-q)/2 (bcf-b-ridge's rule, its q=1 instance
bitwise-identical to BCF's a-move). Engine ~180 net dense-equivalent, tests
~351.

THE IN-SLICE RULE'S OUTCOME: bitwise on the ASIS half; SPECIALIZED-PATH-KEPT
on the conditional half - bitwise is IMPOSSIBLE for the general drawGlue
loop, a compiler fact: the a-block accumulates in one fused statement while
the b-block forms per-row products before a branch and fuses unevenly; under
-ffp-contract=off all four accumulators agree. The reviewer CONFIRMED
impossibility over 21 accumulation variants against the real engine (its own
first standalone probe false-alarmed - vectorization split the multiplies in
the replica, a recorded methodology lesson: contraction probes must run
against the real engine build, not replicas). The measured split, recorded
exactly at the branch-deletion trigger comment (combiner.hpp ~:714-726,
corrected at review): all four PRECISIONS reproduce; the divergence is in
the MOMENTS - unweighted n1 reproduces/n0 differs, weighted both differ. BCF
keeps its two-scalar draw explicitly; BCFSpec::generalAmplitudeDraw is the
one-line switch a future re-record flips. NO baseline moved.

The b-move ships as code, held OFF for BCF (BCFSpec::ridgeB): enabling costs
a GIG draw/sweep = a re-record, and its own acceptance gate
(bcf-b-ridge.md:438-449 - IACT payoff, bcf-exact mode-2b, keepTrees
round-trip) was not run here. DOOR, not a fork: flips only if a measured
mixing case is named; bcf-equivalence's 12/12 bitwise on the glue channel is
the standing pin that no GIG draw entered the BCF stream. Reviewer concurred
with OFF.

Handoff items (a)-(d) all closed: combinedFits untouched (body
byte-identical); a()/b0()/b1() index through amplitudeOffset (aliasing
red/green run by implementer AND reviewer - only the multiplier arm
discriminates, the round-trip arm's limitation is recorded in-test); all
three accumulation directions documented as contracts in place; every
per-forest array sized by the chain's forest count (pre-M4.2 resize(2)
reproduces an ASAN heap overflow at K=3, current clean - both ran it).

Review's required changes, both applied: the trigger-comment correction
above, and TEETH for the LDL' helper - testUnitLowerFactorization
(test_model.cpp, p=1 bitwise-equals-b/P with the Cholesky route asserted to
differ, plus a p=2 orthogonal arm) AND an engine-level assertion that the
general path's q=1 block is BITWISE the shipped a draw (test_sampler.cpp
~:3496-3503); the red/green showed a substituted Cholesky draw turns
EXACTLY that assertion red with the other 241 ok lines green - previously
invisible to every suite. Also from review: drawShippedGlue's
prior-variance refresh now carries the same halfCauchyScale > 0.0 guard as
drawAmplitudes (at the reachable aPriorScale = 0 edge the two paths were
different MODELS; guard alignment verified rng-neutral - no baseline sets
it); installTreatment's Doxygen now states that the setTreatment route
rebuilds basis and layout from scratch, silently resetting an
installForestBasis widening.

Deviation accepted: installForestBasis (engine-internal, unreachable from
any public surface - reviewer verified the shipped setForestBasis routes to
setTreatment) exists because the q0 != 1 fixture is otherwise
unconstructible; it is the engine half M4.3 wires up.

docs/design/bcf.md still describes the a-move as prognostic-only - left for
M4.5 by design.

Battery: implementer and reviewer both from clean - tinytest 4582/0 twice,
trio bitwise twice (37 strict/12/10), ASAN zero twice, seven protected pin
bodies byte-identical across revs.

M4.3 LANDED 9c63e9d8, CI six-green sanitizers included, amended once after
independent review (LAND-AFTER-CHANGES -> four dispositions applied). The
K-length spec surface: setTreatment RETIRED as a mutator at all four layers
(combiner/chain/sampler/facade) -> setForestBasis(f, values, q);
installTreatment became the constructor-only synthesizeIndicatorBasis;
installForestBasis is the SOLE mutator with index/width/finiteness guards;
glue_.z DELETED (drawShippedGlue partitions from basis[1] - the trio proved
the bitwise claim); draw-path selection by the per-forest basisIsCanonical
VALUE predicate, recomputed at install and restore, never serialized;
bcfGlue re-signed to totalAmplitudes/numForestAmplitudes/amplitudes with a
K=2 thin adapter; BCFSpec gained the K-length forests vector with
expandForestSpecs (all 25 BCFSpec fixtures untouched); data@treatment
retired for a data@bases LIST riding creation - which is how
restore-then-widen is met with no fourth reapply hook; the persistence
guarantee survives (the adapted pin asserts VALUES through the combination,
review-verified not weakened). The bcf-equivalence harness migrated
IN-COMMIT; all baselines untouched.

The seven interaction arms: 1-5 in testForestBasisOrdering (replacing
testForestBasisSynthesis, coverage preserved), 6-7 plus the state/glue/probe
legs in the new inst/tinytest/test-forest-basis-r5.R (48 results after the
review round).

Refusal ledger AS THE REVIEW RE-DERIVED IT (the implementer's 6/1/4/34 claim
corrected): 5 relax, 3 restate, 4 rewritten by the slot retirement, 33
stand, and ONE refusal dropped outright with licensed successors - the
treatment-coding refusal, replaced by length/finiteness refusals at
creation. consumer.c: 4 legs moved (2 inverting), LEG_COUNT 18 -> 19, plus
the review's REQUIRED fix - the legRefusals initializer carried 20 elements
for LEG_COUNT 19, a silent-vanish hazard under -Werror toolchains in the
LinkingTo canary (one line, deleted).

Protected-pin record, corrected: the plan's protected set is the 22
non-rewritten BCFSpec fixtures; 11 cpp pin bodies verified byte-identical by
BOTH implementer and reviewer; exactly the three licensed rewrites moved;
three additional moves were FORCED by the licensed ragged ChainStateData and
are claim-accuracy corrections, not lost guards - test_fuzz.cpp
(state.chains a -> amplitudes), common.cpp structSize 272 -> 344 (the
tripwire designed to force it), statesAgree.

Review's other required fix, the M4.0 lesson recurring: Arm 5(iii) was a
DEAD PIN (canonicalA/canonicalB from identical expressions, true under any
implementation) - rewritten to hold the treated column fixed and move the
control column off its determined 1-z, with the red/green shown. Coverage
additions from the review: the ridge derivation pinned BY NAME (basis-free
(1, sd) vs basis-carrying (apv, 0), proven read by flipping the transported
scale); K-length numTrees (a (13,31,7) spec reports 13); a positive
continuous dbartsData(bases=) arm asserted through fitted-value
reconstruction.

The ridgeB question RESOLVED at review: no silent enablement is possible -
on every creation route scale-mixture holds if and only if a forest is
basis-free (R/model.R:874 writes the amplitude prior scale as 0 whenever a
basis is declared; R_interface_bartcore.cpp:2205 reads it into
amplitudePriorScale and :2207 gates ridge on it), so every basis-carrying
forest gets ridge = false; the one reachable q>1 scale-mixture state
(post-creation widening of a basis-free forest) consumes the SAME single
GIG draw already taken at q=1 at the M4.2-validated exponent - no stream
moves. The M4.2 door stays shut; the trio is the standing evidence.

Deviations adjudicated: $getForestAmplitudes(forest = NULL) SATISFIES the
collision license (takes the argument; the universal 3 x n.chains Rd claim
is gone; the unchanged K=2 shape left the two shipped pins green - more
evidence, not less); the K-length params transport derives each forest's
ridge from scale-mixture-ness (the correct trigger per the resolution
above); the sampler numTrees fix reads the expanded spec (caught by a red
setControl pin, not the plan).

BUDGET UNITS CORRECTION: the implementer's 604/632/314 were RAW nets; the
plan's bands are DENSE-EQUIVALENT. Dense engine net is ~165 vs the ~250-300
band (85-135 UNDER - most of the raw is Doxygen); tests raw 603 in the
590-660 band under item 6's own raw-span convention. No overage existed;
record the units so the next slice prices in one currency.

Minor open items recorded (review finding 5 residue): the restore-recompute
clause of the canonical predicate is satisfied but VACUOUS (restoreGlue
never touches basis values; the load-bearing recompute is the install-time
one) - a doc-level note; the dbarts.h transpose fix landed doc-only, hash
verified three ways (static_assert, consumer.c legs, the test-capi.R
literal pin; 188/0).

Battery: implementer and reviewer both from clean - tinytest 4646/0 final
(4633 pre-amend), trio bitwise twice through the migrated harness, ASAN
zero twice, R CMD check --as-cran OK 0/0/0 on a git-archive tarball, pkgdown
clean, no new Rd topic.

M4.5 LANDED fcfd4e29, 2026-08-14, the last Gaussian-complete slice of the
arc, amended once after independent review (first-pass verdict LAND AFTER
FIXES on fourteen defects; fix-round verdict LAND, four residual NITs, none
touching src/). Docs sweep for the K-forest amplitude family M4.0-M4.3
built. New: docs/design/multiplier-combiner.md, the design doc for the
general family - the model equation, the ragged amplitude layout, the
reparameterization and its exact-zero snap, the q-variate amplitude
conditional, the per-forest ASIS rescale, the canonical-basis VALUE
predicate, install and persistence semantics, the accumulation contracts,
and bcf as the K = 2 instance. Swept: the two-forest language across
forest-combiner.md (its re-carves section and its "first instance"
section, now a pointer at multiplier-combiner.md), model-space-survey.md
(a new closed door D4, the four multiplier classes that are the family's
five published instances minus bcf), bcf.md (the a-move-as-prognostic-
only correction plus two further stale paragraphs), and feature-matrix.md
(the bcf row label plus a bounded anchor re-verification). Ten dbarts.h
Doxygen paragraphs changed meaning under K > 2; all comment-only, all
outside DBARTS_C_API_LIST.

feature-matrix.md moves NO CELL: the family is still gaussian-only, so the
new K-length creation route reaches no new response family -
R/spec.R:423-431 refuses data@bases off any non-gaussian family, and
chain.hpp:702-705 builds GaussianResponse and sets
family_ = ResponseFamily::gaussian unconditionally on every K-forest route.
That discharges the matrix's per-landing mandate for this slice; only the
bcf row LABEL and the spec-named COM/CH anchors moved, both non-cell edits.

The five orchestrator decisions. Q1: M4.5 touches src/ after all -
dbarts.h is already outside the docs-only no-workflow path, so the
two-line combiner.hpp:1206-1210 comment fix (a dangling docs/design/bcf.md
IACT citation, repointed to multiplier-combiner.md) cost nothing extra and
was taken. Q2: all of forest-combiner.md, not just the plan's named
re-carves section - leaving the file documenting one combiner instance
right and one wrong was not acceptable. Q3: a BOUNDED anchor pass -
feature-matrix.md carries 59 COM/CH occurrences (~31 distinct), over the
15-anchor threshold for a full pass, so only the spec-named four anchors
plus the bcf row label were re-verified BY SYMBOL; the alias-table
sentence now reads "EXCEPT where a cell or footnote says otherwise," the
Status line's "current at" stamp was deliberately NOT bumped, and TODO
gained feature-matrix-anchor-refresh for the outstanding full pass. Q4:
all four classes written, evidence labelled honestly - continuous/
dose-response BCF's only source (arXiv:2007.09845) is UNPUBLISHED by the
plan's own correction, so D4 states that plainly rather than inheriting
the originating phrase's "verified." Q5: the naming debt RECORDED, not
acted on - multiplier-combiner.md states the family is general and the
spelling is bcf's (BCFSpec, BCFForestCombiner, createBCFSampler,
ChainStateData::hasBCF, the "bcf" state block), prices a rename
(--preclean, a structSize move, a state-block key needing the same
in-place re-encode M4.3 used), and TODO gained bcf-naming-generalization;
no rename was started.

Scope EXPANDED. The plan's M4.5 bullet named five deliverables - the new
file, and edits to forest-combiner.md, model-space-survey.md, bcf.md, and
dbarts.h. The scoping pass found eight further defects, F1-F8, all taken.
The single most important: dbarts_sampler_create's Doxygen stated
numForests == 2 as an invariant against a shipped K = 3 route
(test-forest-basis-r5.R:145, :236-238 build one through dbarts(...,
forests = list(forest(), forest(basis = ), forest(basis = )))) - a false
invariant in the only shipped header, the LinkingTo entry point stan4bart
and its kin call. Three design docs were found contradicting THEMSELVES:
forest-combiner.md's heteroscedastic paragraph called a decision
still-deferred that its own "Anticipated" section already recorded
landed; bcf.md said two setTreatment spellings "did NOT move" when both
had, one of which the same file already recorded moving elsewhere;
model-space-survey.md proposed the 1e-9 multiplier floor as future work
twice while its own earlier update, in the same file, already recorded
the floor replaced by the exact-zero snap. All three corrected.

Budget, in dense-equivalent terms. The plan's docs band is ~250-300, and
the plan's own M4 budget arithmetic records that the band never priced
multiplier-combiner.md as a new file. The content spec priced the six
mandated deliverables (the five plan items plus the feature-matrix check)
at ~371-460 dense-equivalent; the orchestrator's brief folded in all eight
findings and the two TODO entries and set the working budget at ~450-565.
Both figures exceed the docs band on their own; the slice took the full
scope named in the brief, and neither implementer nor reviewer reported an
overage against it.

Gates, as observed locally: tinytest 4646/0, test-capi.R 188/0, tests/cpp
all tests passed. dbarts_apiHash() 0xcd88efcd67de55d7, verified three ways
(the static_assert, the consumer.c raw leg, and the test-capi.R:59 literal
pin). Equivalence trio bitwise: 37 compared / 0 skipped under
--strict-coverage, bcf 12/12 on all channels, multinomial 10/10 on all
channels, zero "max |z|" lines. Every src/ and inst/include/ hunk verified
comment-only (git diff -U0, filtered for non-"///" lines: zero hits) and
outside DBARTS_C_API_LIST, so no version constant moved. ASCII clean
across all ten changed/new files. "Gaussian responses only." survives
verbatim, now at dbarts.h:513 (shifted from :505 by the paragraph's own
growth), and moves with M4.4. CI SIX-GREEN on fcfd4e29, sanitizers
included - this slice touches inst/include/ and src/, which the path
filters (paths-ignore: docs/**, TODO, **.md, benchmarks/**) do not
exclude, so unlike a pure docs push it fired the full battery and the
evidence is runner-executed rather than inferred from byte-identity.

THE STANDING LESSON: a slice that edits files AND cites them by line
number invalidates its own citations mid-flight. Independent review caught
three clusters of this on the first pass - the dbarts.h header edits
shifted every anchor past the first hunk by up to +21; a three-line
correction to this plan (the column-major -> row-major fix at :557-560)
shifted every plan anchor above :558 quoted in the new prose by +3; and
model-space-survey.md's OWN earlier edits, made earlier in the same
commit, shifted its own new D4 section's citations by +8 and +9. The fix
pass re-derived 141 citations and corrected 60 - and was not itself
immune: a single Doxygen comment fix in combiner.hpp (:898-902, four
lines out, five in) shifted MultinomialForestCombiner::setActiveRows off
the exact line the review had just given it. The discipline this leaves:
re-derive every anchor by OPENING THE TARGET, never by applying an
arithmetic offset - two plan ranges in the fix round correctly did NOT
move by +3, the evidence re-derivation happened rather than arithmetic -
and run the anchor pass LAST, after every content edit is final. The
reviewer's own hand sample, checked the same way, was 62 citations across
12 files with 60 landing exactly.

M4.4 LANDED 4da3bd8a, 2026-08-14. The arc's last slice, and the one that
makes the general K-forest basis/amplitude family non-Gaussian: `probit`
and `logistic` now reach it, at creation and through every family-keyed
predicate, leaving `aft`, `ordinal` and `nbinom` refused by name at the
door. The calibration map gains `latentScaleAnchor(family)` - the settled
Option L, `s = 1` under probit and `s = pi/sqrt(3)` under logistic, and
`sd(y)` under gaussian - so the map's node scale is
`nodeScaleFactor * s / (nodeScaleDivisor * basisRowNorm)` at every family
rather than at one. Under a latent family the forests combine into the
INDEX rather than the mean, sigma is pinned, the response transform is the
identity, and there is no response standard deviation for a prior scale to
be stated in; the "sd(y) units" phrasing was therefore a required edit and
not a sweep, on the one shipped knob (`forest(sd = )`) above all. Folded in
by VD decision rather than split out: the live `$getSumsOfSquaredResiduals()`
defect, which read `forests_[0].totalFits` and so reported forest 0 bare
with the amplitude never applied, 1 to 49 percent off and contradicting its
own Rd. `dbarts.h` is comment-only, so no signature moved, no hash re-baked
and no `structSize` moved; `SamplerFacade<ConstantGaussianLeaf>` is a LEAF
instantiation and the family is runtime, so the slice adds no instantiation
at all.

ARM E, the slice's OWN statistical acceptance gate, PASSED and is written
up in `.claude/m4-basis-design/arm-e-2026-08-14.md` (gitignored; the numbers
below are the record). The native K-forest probit sampler AGREES with FA5's
independent reference on all 12 functionals at max |z| = 1.46, and with arm
B - the R-composed arm that FALSIFIED this slice's original headline
justification - at max |z| = 1.73. Zero functionals over |z| = 3 on either
comparison, so the pre-registered AMBIGUOUS branch is not entered and no
4x-seed re-run is owed. The power precondition is MET with margin: arm A
differs at |z| = 772.49, and at 588.64 re-standardized in arm E's OWN
standard errors, so the verdict is not inherited from an arm that mixes
differently. E0, the batch-vs-loop identity, holds. The open question arm E
existed to answer is ANSWERED: the pinned sigma costs nothing measurable on
the fixed-variance basis forest. Arm E's IACT is at or below arm B's on 8 of
12 functionals (mean ratio 0.98) and 1.11x the reference's on average; the
one lag, Etau3, is shared in equal measure by BOTH dbarts arms, so it is the
MH tree move against the reference's exact enumerated conditional and not
the engine path and not sigma. That is the mixing and coverage evidence
`binary-kforest-prior-default` was ticketed to wait for, so its TRIGGER IS
MET. What arm E does NOT establish, stated because a later reader will
over-read it: it pins `sd` at exactly the point where the map returns
N(0, 1) and then asserts the map returned it, so a map wrong everywhere
except that point would pass it untouched; it exercises probit's anchor
only, at basis row norm exactly 1, with `update.amplitude = FALSE` in every
arm so no amplitude conditional and no ASIS ridge ever ran. Arm E is NOT the
anchor's gate. The anchor's gate is checklist item 25, which ships in
`test-bcf-family.R`.

ITEM 25'S GATE STRENGTH IS SUBSTANTIATED, by mutation builds run
independently rather than asserted: GREEN as shipped at exactly 76/76, and
RED under all three mutations, each failing on the assertion designed to
catch it and none failing incidentally - the naive anchor (every family
taking the gaussian answer) moves 14 assertions, Option C moves 12, and a
dropped row-norm divisor moves 4, all inside the anchor blocks. The three
assertion types are non-redundant: Option C correctly does NOT trip base-rate
invariance, since it is still a family constant, and is caught instead by the
absolute anchor and the induced-index check. Two limits of the gate are
recorded because a later reader will lean on it. First, the row-norm coverage
hangs on ONE deliberately planted fixture cell, `scaledBases` carrying a
`4 * zBasis` whose row norm is 4; delete that cell and the dropped-divisor
mutation goes fully GREEN. That is the arc's own "unit-valued fixture factors
silently vacate pins" lesson, defended against on purpose - do not let a
later tidy-up remove it. Second, the naive anchor's LOGISTIC arm misses the
absolute anchor by only 0.865 percent (2.66782 against 2.69110), so it is
caught only by the 1e-12/1e-4 tolerances plus the 33.8 percent divergence in
base-rate invariance; a looser tolerance would let the naive anchor through
on that arm.

GATES, re-run independently on the merged tree into a FRESH private library
rather than the implementer's: `R CMD INSTALL --preclean` clean with no
compiler warnings; `tests/cpp` 241 ok / 0 fail from `make clean`, and again
241 ok under ASAN+UBSAN with zero sanitizer reports; tinytest 141 files,
4732 tests, 0 failures; full `R CMD check` from a clean copy staged outside
the tree, Status OK with 0 ERRORs, 0 WARNINGs and 0 NOTEs, examples,
`tinytest.R` and the vignette rebuild included; `air format --check` clean
tree-wide (with the checker itself proven live against a deliberately
mis-formatted scratch file); `NEWS.Rd` parses; `pkgdown::check_pkgdown` finds
no problems, and the slice adds no NEW exported Rd topic so no `_pkgdown.yml`
entry was owed. **THE EQUIVALENCE TRIO IS BITWISE ON ALL THREE SUITES with
NO "max |z|" line anywhere: `equivalence-8b047f8b` 37/37 identical draws at
37 compared / 0 skipped, `bcf-equivalence-8b047f8b` 12/12 identical on every
channel, `multinomial-equivalence-1027be5` 10/10 identical on every channel.
No baseline was re-recorded.** So the answer the gates bullet demanded is
recorded plainly for the next slice to inherit: M4.4 IS DRAW-NEUTRAL, proven
at this commit, the folded-in residuals fix included - it moves a REPORTING
path, not a draw path. `lintr` on `R/spec.R`, the only touched R file,
reports NO LINTS against a private library carrying the slice's own build,
identically on the working tree and on `git show HEAD:R/spec.R`; the nine
`object_usage_linter` warnings an earlier run showed were the documented
FALSE-POSITIVE class, reproduced only by linting against a stale installed
dbarts that lacks `resolveForests`, `forestParams` and five more of this
arc's own internals. Speed was deliberately NOT measured and no speed claim
is made: the machine was loaded, and `bench-sampler` needs genuine quiet.

BUDGET against item 9, measured by independent review in DENSE-EQUIVALENT
non-comment net terms: engine 67 against a 49-82 band and a 123 stop, PASS;
bridge about 33 dense against a 14-28 band and a 42 stop, OVER BAND and under
the stop - and the currency matters here, because 42 is the bridge's RAW net
and reporting it as dense would have put the slice AT its own halt condition
on a figure in the wrong unit; R 27 against 18-38, PASS; flat C 0, PASS; docs
about 100 (markdown at raw/2 = 74, Rd raw = 26) against 100-175, PASS; tests
348 against a 155-260 band and a 390 stop, OVER BAND and under the stop. The
two over-band layers are both the arc's known under-pricing of test and
bridge work against a mandated oracle, not scope creep. Checklist items 1
through 25 are all present and spec-conformant, item 14 correctly shipping no
code. Two DELIBERATE deviations, recorded so neither reads as a miss: item
21's helper is `basisRowNorm(const double*, numColumns, n)` rather than the
specified `(const ForestSpec&, n)`, which item 22's call site requires; and
item 7's row for `man/dbartsSampler-class.Rd:147` was NOT edited because the
Rd's claim is not in fact falsified - post-creation `setWeights` is still
gaussian-only and the page already defers to `dbarts` for creation-time rules
- so the SPEC's row was wrong there, not the code. A weak spot in the new
tests is recorded rather than papered over: no test exercises a non-default
`forest(sd = )` under a latent family, so `nodeScaleFactor` is 1 throughout
the anchor assertions and a factor-versus-anchor confusion would go unpinned.
That is spec-conformant - item 25 specified this fixture - and it is added to
`binary-kforest-prior-default`, which is the slice that will move `sd`.

THIS PLAN'S OWN TEXT WAS WRONG IN NINE PLACES, all found by opening targets
rather than by trusting the citation, and all corrected HERE rather than in
the bullet, because `multiplier-combiner.md` carries eighteen citations into
this file on seventeen lines and an in-place edit would falsify them. One is
a behavioral prediction, not an anchor: item 7 predicted that
`bcf` x `updateScale = TRUE` would take the shipped latent convention and be
IGNORED rather than refused, since latent families have `fitScale() == 1`
and `fitShift() == 0`. It is REFUSED, under every family. `refuseBCFMutation`
keys on the sampler carrying bases and `refuseMultiForestResponseMutation` on
`numForests >= 2`; neither has ever consulted the family, so the cell does
not move and the reason is that the refusal was never family-keyed. M4.4's
own new test already pins it. The other eight are anchors: item 7's
"sd(y)-unit calibration" third instance is at `multiplier-combiner.md:28`,
not `:25`; `man/dbartsSampler-class.Rd:147` is listed among the nine
assertions M4.4 falsifies but is already FAMILY-scoped rather than
forest-count-scoped, so it was correct as written and took no edit; 2d(a)'s
`R/model.R:772-840`, `:799-806` and `:874` and
`R_interface_bartcore.cpp:2207` are really `:810-878`, `:837-843`, `:915`
and `:2208`; `[f23]`'s `spec.R:440` named the non-default-`k` entry rather
than the `prior.scale` refusal; `[f36]`'s `CH:1589` is the variance-forest
pre-step and not the forest loop; the DART refusal is at `spec.R:451` fired
at `:481-487`, not `:441`/`:468-474`; and the header's family documentation
is `dbarts.h:485-494` with the K-forest sentence M4.4 replaced at `:514-518`.

A COMMENT-ONLY residue pass was folded in rather than ticketed, on the
finding that six in-code sites still asserted what M4.4 falsifies - all
outside the bullet's tabled nine, and all within about forty lines of edits
the slice did make, which is the near-miss pattern this arc keeps producing.
The one that mattered is `facade.hpp`'s `SamplerShape::supportsResponseMutation`
docstring, a CAPABILITY docstring a `LinkingTo` host reads, still saying the
flag is false for a non-gaussian response after the slice dropped exactly
that conjunct. The others: `BCFForestSpec`'s "derives them from the response
sd", the K-forest constructor's "at aPriorScale sd(y)" and "sits at
sdModerate sd(y)", `BCFState`'s "y = a mu + b_z tau + eps", and two
`R/spec.R` rationales stating leaf scales come from the response sd. Zero
executable tokens moved: `git diff -U0` yields no changed non-comment line
and each file's non-comment content is md5-identical to its pre-pass state.

THE ANCHOR SWEEP FOUND ROT OLDER THAN THIS SLICE. Of about 199 code
citations checked across the touched design docs, 105 changed, and roughly
seventy of those in `multiplier-combiner.md` were ALREADY WRONG AT HEAD by
up to +78 lines - the file was written against a pre-M4.4 `combiner.hpp` and
no pass since M4.1 had swept it. Independent review then opened 50 of the
moved citations across `combiner.hpp`, `chain.hpp`, `facade.hpp`, `R/spec.R`
and `dbarts.h` and found ZERO mismatches. The deltas are not identical -
fourteen distinct values in `multiplier-combiner.md` alone (+1, +2, +8, +15,
+16, +18, +19, +25, +26, +27, +31, +35, +78 and -4), twenty-two across both
files, ten into `combiner.hpp` by itself - and endpoint deltas differ from
start deltas inside untouched regions, so this was neither an offset pass nor
a hunk-line-map pass. One anchor changed FILE (`COM:963` to `CH:955`), which
no arithmetic could produce. **The sweep's own "31 citations left unmoved"
figure did NOT reproduce and is struck**; measured, it is 18 unmoved in
`multiplier-combiner.md` (4 of them into files this slice edited) and 423 in
`feature-matrix.md` (66 into edited files). The review's re-count is what
turned up the two anchors this slice itself BROKE (`COM:1581` and `COM:486`,
correct at HEAD~1) and one the by-symbol pass had missed
(`multiplier-combiner.md:403`), all three fixed here.
`nameable-calibration.md`, the `man/` files and `inst/NEWS.Rd` carry no
source line citations at all.

`docs/design/feature-matrix.md` MOVED CELLS, discharging the standing
mandate that M4.5 could only record as vacuous. Fourteen cells in the `bcf`
row, seven footnotes rewritten (`[f3]`, `[f17]`, `[f18]`, `[f23]`, `[f26]`,
`[f33]`, `[f36]`), a new `[f48]` carrying the whole family-reach story, and
an "Exception on record" entry naming what was re-derived BY SYMBOL and what
was not. `[f26]`'s premise - "a BCF sampler's response IS a
`GaussianResponse`" - is falsified; its conclusion survives, with its
measurements now flagged gaussian-only. The Status stamp was deliberately
NOT bumped: the outstanding whole-file pass is still
`feature-matrix-anchor-refresh`, whose remaining scope the sweep enumerated
by anchor rather than leaving it as a count.

TICKETED, NOT FOLDED IN: `R_interface_bartcore.cpp`'s `calibrationMapName`
returns the literal "two-forest calibration map" for every non-softmax
coupling, and it reaches a user through the `setForestPriorScale` refusal -
a K = 5 sampler is told its scale comes from the two-forest map "which owns
both halves of its calibration". Wrong for K != 2, but wrong since the
general K-forest family shipped and not falsified BY M4.4, so it joins
`bcf-naming-generalization` rather than widening this slice.

THE STANDING LESSON, restated because this slice re-learned it in a new
shape: M4.5's rule was that a slice editing files it also cites invalidates
its own citations mid-flight. M4.4 found the CONCURRENCY corollary. Two docs
agents were deriving anchors into `combiner.hpp`, `chain.hpp`, `facade.hpp`
and `R/spec.R` at the same time a comment-only pass was editing those four
files, so citations were correct when written and stale when read, and one
file's anchors drifted twice inside a single session. The scheduling error
was the orchestrator's. It cost one extra sweep and no correctness, because
the anchor pass runs LAST by rule and the rule held - but the cheap fix is
to serialize any pass that edits a cited file against any pass that cites
it, rather than to rely on the final sweep to catch it.
