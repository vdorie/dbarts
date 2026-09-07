# Backlog value scan, 2026-08-24 (branch bartcore, tip ee70e6f7)

Read-only. Every status claim checked against live code or git at ee70e6f7; all 43 shas cited
in TODO exist. Verdicts: PRE-RC (surface + named value), POST-1.0 (additive, no lock-in),
NO-SURFACE, KILL, DOOR (VD decides), PROCESS. Counts, TODO entries only: 17 closed, 7 PRE-RC,
3 DOOR, 4 POST-1.0, 8 NO-SURFACE, 1 PROCESS. Sec 3 adds 4 PRE-RC and 1 DOOR found outside the
TODO.

## 1. Entry by entry

### Closed, verified

- augmentation-helpers ([TODO:31](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/TODO#L31)) 890efd3d; valgrind residue rides the merge.
- bcf-naming-generalization ([TODO:39](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/TODO#L39)) 8215a53c + be1d06e4; two dead static_assert strings still
  read "BCF is a constant-leaf model" - sweep on next touch.
- blocked-jacobi-trees ([TODO:80](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/TODO#L80)) KILLED df13c04 + 77597cc.
- cheap-uniformity ([TODO:84](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/TODO#L84)) S3 06c0254 + S4 9929ede, ARC COMPLETE; sparse predict SHIPPED (the
  entry's "all open" means open for use).
- feature-matrix-prose-staleness ([TODO:118](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/TODO#L118)) 6d82ac86; new staleness found, sec 3.
- forest-ranef-interweaving ([TODO:123](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/TODO#L123)) NO-GO with a reopen clause.
- grow-from-root-default ([TODO:145](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/TODO#L145)) KILLED 2026-08-08.
- host-shell-read-guards ([TODO:150](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/TODO#L150)) OBVIATED; zero hostFor / $bc hits in R/.
- latent-family-weight-channel ([TODO:164](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/TODO#L164)) d0701a6a; its door is in flight, skipped.
- multiforest-extension-surface ([TODO:170](https://github.com/vdorie/dbarts/blob/d0701a6a829c289c65965ed8af14e30704420010/TODO#L170)) complete; FA5 respec is design-only.
- multiforest-predictor-mutation ([TODO:236](https://github.com/vdorie/dbarts/blob/d0701a6a829c289c65965ed8af14e30704420010/TODO#L236)) S4 031184ce; one residue STALE, sec 3.
- multinomial-counts-mutation ([TODO:279](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/TODO#L279)) complete; its doors in flight, skipped.
- r-c-division ([TODO:327](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/TODO#L327)) guide accepted, slate landed.
- rc-gate ([TODO:336](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/TODO#L336)) slate complete 2026-08-19; declaration VD-held.
- setpredictor-leafof-rebuild ([TODO:415](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/TODO#L415)) closed; test-bartcore-filtered-run ([TODO:436](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/TODO#L436)) e74380f7;
  x86-simd ([TODO:502](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/TODO#L502)) closed as a perf lever.

### Open

- NO-SURFACE, verified open and none of them adds or moves public surface:
  bcf-sigma-tail-mixing ([TODO:58](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/TODO#L58), [docs/plans/bcf-sigma-residual.md:164-166](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/docs/plans/bcf-sigma-residual.md#L164-L166) - a regime arising only from a
  stale build scale under setResponse(updateScale=FALSE)); correlated-outcomes ([TODO:88](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/TODO#L88) - AR-1 is
  an engine door, panel/ITS compose in the outer sampler today); gpu-bart ([TODO:134](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/TODO#L134), waits on
  frontier 4-5); grow-from-root ([TODO:143](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/TODO#L143), frontier 3.1 informed proposals); python-bindings
  ([TODO:320](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/TODO#L320), no python tree, two host contracts named); sbc-family-tiers ([TODO:412](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/TODO#L412), BUILT - open doors
  are the monotone arm, the heteroscedastic setState lift blocked by the same defect as
  shortlist 2, and full-R ladder re-runs); tree-mixing-proposals ([TODO:441](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/TODO#L441), research complete -
  standing note: a proposal move shipping on by default moves every user's draws at a fixed
  seed, so if one is ever taken it is a pre-1.0-shaped change, and the research is not ripe
  enough to schedule that now).
- POST-1.0 per standing calls, all verified: binary-kforest-k1-reachability ([TODO:67](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/TODO#L67), a designed
  refusal in spec.R; surface yes, a new shipped configuration);
  equivalence-harness-statistical-mode ([TODO:96](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/TODO#L96) - CONFIRMED including the awkward half:
  bcf-equivalence-6e3b9fb8.rds (Aug 16) predates the summaries [benchmarks/R/bcf-equivalence.R:495-505](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/benchmarks/R/bcf-equivalence.R#L495-L505)
  needs, and that path degrades loudly; gate infra); repo-modernization ([TODO:410](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/TODO#L410), not started -
  inst/NEWS.Rd, no NEWS.md, no @export); wave5-review-followups ([TODO:466](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/TODO#L466) - no vetoed-leaf counter
  exists in src/bartcore, so item 5's first half is open as a new but additive R-readable
  channel).
- gp-followups ([TODO:127](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/TODO#L127)). CONFIRMED: gp(columns, k, lengthscale, ...) takes a fixed positive
  numeric only ([R/model.R:1454-1461](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/R/model.R#L1454-L1461)). Surface yes (lengthscale accepts a prior + a per-draw
  channel). Value: GP-leaf BART as a sequential-design surrogate. Large Opus arc
  (non-conjugate MH on theta). PRE-RC.
- group-by-exposure ([TODO:137](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/TODO#L137)). CONFIRMED: group.by is a formal on exactly three rbart surfaces
  ([R/rbart.R:17-18](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/R/rbart.R#L17-L18), [R/generics.R:1495](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/R/generics.R#L1495), [R/bart.R:2535](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/R/bart.R#L2535)); the engine sees attr(control,
  "bartcore.groups") ([R/rbart.R:388](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/R/rbart.R#L388), read at [src/R_interface_bartcore.cpp:1910](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/R_interface_bartcore.cpp#L1910)). Surface yes - a
  NEW formal on bart2()/dbarts(). Value: clustered/multi-site causal survival and grouped
  discrete-time hazard, neither reachable through rbart_vi (gaussian/aft only). ~200 R + Rd +
  tests, Sonnet with an Opus read on composition. PRE-RC.
- integer-predictor-storage ([TODO:153](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/TODO#L153)). CONFIRMED absent: PredictorSource carries const double*
  denseValues only ([src/bartcore/data.hpp:188](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/data.hpp#L188)), no integer ColumnSourceKind ([src/bartcore/data.hpp:129-135](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/data.hpp#L129-L135)), bridge
  and R both coerce ([src/R_interface_bartcore.cpp:1046](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/R_interface_bartcore.cpp#L1046), [R/data.R:204](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/R/data.R#L204)). storage= is residual
  precision, not predictor storage. NO-SURFACE.
- interaction-constraints ([TODO:157](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/TODO#L157)). CONFIRMED: shipped formal is interactions(max.order, groups,
  forbid) ([R/model.R:1672](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/R/model.R#L1672)); zero "heredity" hits in R/ src/ man/ inst/. Surface yes - a FOURTH
  formal plus a prior change, not a value on an existing one. Value: functional-ANOVA BART
  (dbarts ships the hard-constraint half). Soft penalties: gate stands, no model named. Opus
  engine + R, medium-large. PRE-RC (heredity half only).
- multiforest-mutation-gaps ([TODO:208](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/TODO#L208)). All four claims CONFIRMED. Flat C
  setPredictor/updatePredictor present and unguarded BY DESIGN with the rationale in place
  ([src/C_interface.cpp:665-673](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/C_interface.cpp#L665-L673), [src/C_interface.cpp:694-698](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/C_interface.cpp#L694-L698)); monotone 1e-9 deliberate; BCF whole-data setData
  postponed, not raised. Per-forest row subsetting: value lifted, shape unresolved. DOOR (D1).
- negbin-real-dispersion ([TODO:296](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/TODO#L296)). CONFIRMED live: [R/data.R:650-655](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/R/data.R#L650-L655) refuses real dispersion, and
  [R/data.R:634-639](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/R/data.R#L634-L639) already records that admitting it is "a validation relaxation, not a signature
  change". DOOR (D2, with weighted-binary).
- ordinal-scan-missing-rows ([TODO:305](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/TODO#L305)). CONFIRMED live and exactly as described. Split side skips
  missing (`if (code == naCode) continue;`, [src/bartcore/scan.hpp:100](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/src/bartcore/scan.hpp#L100)) and totals surviving bins
  ([src/bartcore/scan.hpp:105-108](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/src/bartcore/scan.hpp#L105-L108)); the no-split term takes the FULL node stats ([src/bartcore/grow.hpp:183-193](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/grow.hpp#L183-L193), fed by
  setNodeAverages over the whole index buffer, [src/bartcore/chain.hpp:1951](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/chain.hpp#L1951)). Categorical exempt by design
  ([src/bartcore/scan.hpp:217-220](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/scan.hpp#L217-L220)). Reachable only via $growFromRoot() / bart2(n.grow.sweeps > 0) (default
  0L, [R/bart.R:677](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/R/bart.R#L677)) on an ordinal column with NAs. Surface: NEWS-visible and DRAW-LAW-CHANGING
  there. Value: MIA + ordinal + grow-from-root is a shipped combination now scoring splits
  against a mismatched no-split term. ~60 engine + tests, Opus. PRE-RC.
- release-candidate-review ([TODO:372](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/TODO#L372)). All six VD forks executed (slice F 621da478, FX1-channel
  c3af16a1, FX2, L, P12, merge deferred); waves 0-5b landed; the remaining legs are in flight
  elsewhere. PROCESS.
- sparse-extensions ([TODO:419](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/TODO#L419)). CONFIRMED, both refusals live: rbart_vi at [R/rbart.R:302-307](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/R/rbart.R#L302-L307),
  linear/gp leaf covariates at [R/model.R:215-220](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/R/model.R#L215-L220) and [R/model.R:257-263](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/R/model.R#L257-L263) (engine reason
  [src/bartcore/data.hpp:287-295](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/data.hpp#L287-L295)). Surface: no formal moves, a NEWS-visible widening of
  accepted inputs. Value: grouped random effects on high-dimensional sparse designs. Perf
  halves: gate stands, no numbers. PRE-RC (rbart_vi + linear-leaf halves only).
- survival-followups ([TODO:427](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/TODO#L427)). CONFIRMED: expandDiscreteTimeHazard ([R/dbarts.R:121-129](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/R/dbarts.R#L121-L129)) has no
  entry-time argument and starts every subject at period 1 ([R/dbarts.R:158-159](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/R/dbarts.R#L158-L159)). It is INTERNAL, so the
  entry's "one-argument expander change" understates it - the user needs a NEW public formal
  (entry=) beside breaks/max.rows. Value: registry/EHR staggered entry. Competing risks moves
  with multinomial; cloglog and continuous-time gates stand. ~60 R + Rd + tests, Sonnet.
  PRE-RC (left truncation only).
- variance-forest-mutation-routing ([TODO:456](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/TODO#L456)). Three doors, all CONFIRMED, not the same kind. (a)
  setState variance COLUMN-MASK gap, LIVE and silent: columnMaskSubtreeIsValid ([src/bartcore/tree.hpp:536](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/tree.hpp#L536))
  has five call sites, all mean-forest ([src/bartcore/chain.hpp:3191](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/src/bartcore/chain.hpp#L3191), [src/bartcore/chain.hpp:3379](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/src/bartcore/chain.hpp#L3379), [src/bartcore/chain.hpp:3400](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/src/bartcore/chain.hpp#L3400), [src/bartcore/chain.hpp:3442](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/src/bartcore/chain.hpp#L3442), [src/bartcore/chain.hpp:3530](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/src/bartcore/chain.hpp#L3530));
  stateIsValid's variance branch ([src/bartcore/chain.hpp:3231-3275](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/src/bartcore/chain.hpp#L3231-L3275)) checks well-formedness, positivity and occupancy
  but never the mask, and rebuildVarianceForest ([src/bartcore/chain.hpp:4285-4316](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/src/bartcore/chain.hpp#L4285-L4316)) does not either; the one variance
  mask check ([src/bartcore/chain.hpp:3387-3403](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/src/bartcore/chain.hpp#L3387-L3403)) is reached from the warm-start pre-flight ([src/bartcore/sampler.hpp:918](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/sampler.hpp#L918)), not
  from setState ([src/bartcore/sampler.hpp:744-758](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/sampler.hpp#L744-L758)). (b) scale-leaf staleness under updateScale and (c)
  heteroscedastic prior-predictive are both fenced by clean refusals
  ([src/R_interface_bartcore.cpp:2692-2704](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/R_interface_bartcore.cpp#L2692-L2704); [R/dbarts.R:813-824](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/R/dbarts.R#L813-L824)). PRE-RC for (a), POST-1.0 for (b)
  and (c).
- weighted-binary ([TODO:498](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/TODO#L498)). CONFIRMED live: probit refuses weights outright ([R/spec.R:52-58](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/R/spec.R#L52-L58)),
  logistic refuses non-integer ([R/spec.R:60-67](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/R/spec.R#L60-L67)). DOOR (D2).
- release block ([TODO:509-550](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/TODO#L509-L550)). PROCESS, VD-triggered; 11 workflow files exist, registration waits
  on the coordinated merge.

## 2. Ranked shortlist

Ranked by lock-in severity (a wrong number, or a contract we would have to break later), then
named value. Costs are budgets, not gates.

1. Heteroscedastic loglik and PPD are silently wrong. NOT IN TODO. A variance-forest fit
   records family = "gaussian", so extract(type = "loglik") scores dnorm(y, ev, sigma/sqrt(w))
   with the SCALAR object$sigma ([R/generics.R:76-95](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/R/generics.R#L76-L95)), ignoring the per-observation s.train the
   same fit stores ([R/bart.R:243-260](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/R/bart.R#L243-L260)) - and sigma is structurally pinned under a variance
   forest, so this is not near-right. sampleFromPPD ([R/generics.R:2052-2078](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/R/generics.R#L2052-L2078)) has the matching
   hole: it guards resid.dist but not the variance surface, so type = "ppd" adds noise at the
   wrong scale. SURFACE: two exported channels change their numbers or gain a named refusal,
   NEWS-visible. VALUE: loo/waic and PPD checks on heteroscedastic fits stop being wrong. COST
   ~60 R + tests, Sonnet with an Opus read on the draws-axis pairing. RECOMMEND: score it (sd
   = s(x)/sqrt(w), s.train being on the same axis), refusal as the fallback - the posture VD
   took at FX1 and FX2.
2. Variance-forest setState column-mask gap ([TODO:456](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/TODO#L456)a). Public route: $setState /
   dbarts_sampler_setState ([src/C_interface.cpp:867](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/C_interface.cpp#L867)) on a variance = <column subset> fit (mask
   route [src/R_interface_bartcore.cpp:2021-2033](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/R_interface_bartcore.cpp#L2021-L2033) -> [src/bartcore/chain.hpp:4126-4132](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/chain.hpp#L4126-L4132)). The install succeeds; a
   live variance tree then splits on a column its own mask forbids, and
   splitVariableLogProbability ([src/bartcore/model.hpp:2135-2148](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/model.hpp#L2135-L2148)) scores it against a variable set
   excluding that column - wrong MH ratios, wrong s(x), nothing reported (-inf if split.probs
   zeroes it); [src/bartcore/chain.hpp:3347-3349](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/chain.hpp#L3347-L3349) asserts a guarantee that does not hold. SURFACE: a refusal
   or repair on a shipped method; adding it after 1.0 breaks callers. VALUE: state restore and
   warm starts on a column-subset variance forest give correct draws, and the SBC
   heteroscedastic setState lift unblocks. COST ~40 engine + tests, Opus. RECOMMEND: take; the
   backstop is one call.
3. ordinal-scan-missing-rows ([TODO:305](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/TODO#L305), evidence in sec 1). SURFACE: draw-law change on the
   opt-in grow path. VALUE: correct grow-from-root warm starts for ordinal designs with
   missing predictors. COST ~60 engine + tests, Opus; re-record only if a scenario uses
   n.grow.sweeps. RECOMMEND: take - a draw-law fix is a re-record now and a broken seed
   promise after 1.0.
4. dbartsData(bases =) equal-count ambiguity. NOT IN TODO; scoped out of S13
   ([docs/plans/bart2-argument-consolidation.md:1931-1934](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/docs/plans/bart2-argument-consolidation.md#L1931-L1934)). dbarts()'s formula path now requires a
   FULL-DATA basis and subsets it by the model frame's index ([R/dbarts.R:650-659](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/R/dbarts.R#L650-L659));
   dbartsData()'s formula path still checks the basis against the POST-subset count
   ([R/data.R:1067](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/R/data.R#L1067)) while its own x/y path checks full length and subsets ([R/data.R:1165](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/R/data.R#L1165)). The same
   (formula, subset, bases) call means different things at two exported entry points, and at
   equal row counts nobody can tell which was meant. SURFACE: an exported constructor's input
   contract; tightening it after 1.0 breaks callers. VALUE: one alignment rule for
   multi-forest users. COST ~25 R + tests, Sonnet. RECOMMEND: take - converge on the dbarts()
   rule, refuse the ambiguous shape by name.
5. Multinomial arc Door 1, twin-create deletion. NOT IN TODO
   ([docs/design/multinomial-mutation-arc.md:577](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/docs/design/multinomial-mutation-arc.md#L577)). Ordinal and nbinom still build a host
   sampler then a second engine and adopt the second ([R/bart.R:1802-1807](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/R/bart.R#L1802-L1807) says so); deleting
   the redundant bartcore_create MOVES 2 of 43 equivalence scenarios. SURFACE: draw-moving at
   a fixed seed for two shipped families. VALUE: one engine per fit, and the move is a
   re-record today. COST ~60/~40 plus a re-record, Opus. RECOMMEND: now or never.
6. group-by-exposure ([TODO:137](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/TODO#L137)). SURFACE: a new group.by formal on bart2()/dbarts(), spelling
   already reserved by rbart_vi. VALUE: clustered / multi-site causal survival (riAFTBART's
   shape) and grouped discrete-time hazard, unreachable via rbart_vi's gaussian/aft-only
   family list. COST ~200 R + Rd + tests, Sonnet + an Opus read. RECOMMEND: take if the RC has
   room - the largest missing piece of the 1.0 public face.
7. survival left truncation ([TODO:427](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/TODO#L427)). SURFACE: a new entry= formal on dbarts()/bart2() for
   family = "hazard". VALUE: registry and EHR cohorts with staggered entry, hand-expanded
   today. COST ~60 R + Rd + tests, Sonnet. RECOMMEND: take - the cheapest surface item here.
8. sparse-extensions, rbart_vi and linear-leaf-on-sparse ([TODO:419](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/TODO#L419)). SURFACE: two refusals
   removed, NEWS-visible, no formal moves. VALUE: grouped random effects on high-dimensional
   sparse designs; treed linear leaves on sparse. COST rbart_vi half is R plumbing + tests
   (Sonnet), the leaf half needs engine raw reads against CSC (Opus). RECOMMEND: take the
   rbart_vi half, price the leaf half apart - a widening stays compatible, so it is face, not
   lock.
9. rbart_vi has no logistic token. NOT IN TODO (feature-matrix logistic gap). family =
   c("auto", "gaussian", "aft") ([R/rbart.R:49](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/R/rbart.R#L49)), so grouped logistic is engine-reachable and
   R-unreachable and grouped binary silently gets probit. SURFACE: one widened formal value +
   Rd + NEWS. VALUE: multi-site binary outcomes on the logit scale; logistic is the only
   binary family with weights. COST ~30 R + tests, Sonnet. RECOMMEND: take - cheapest item
   here, and it closes an asymmetry that reads as an oversight.
10. interaction-constraints formal heredity ([TODO:157](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/TODO#L157)). SURFACE: a fourth interactions()
    formal plus an ensemble lattice prior. VALUE: functional-ANOVA BART. COST Opus engine + R,
    medium-large. RECOMMEND: only if the RC slips - additive on an exported helper, low window
    cost.
11. gp-followups sampled lengthscales ([TODO:127](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/TODO#L127)). SURFACE: lengthscale accepts a prior, plus a
    per-draw channel. VALUE: sequential design of computer experiments. COST large Opus arc.
    RECOMMEND: post-RC in practice; the widening stays compatible and the arc does not fit the
    window.

DOORS - the decision in plain words:

- D1, per-forest row subsetting ([TODO:208](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/TODO#L208)). DECIDE: does a forest seeing fewer rows get a MASK
  (rows present, zero weight, shared cut grid), a physically COMPACTED store with its own cut
  grid, or nothing pre-1.0? Value is lifted (AFT mixture-cure BCF, whose authors reimplemented
  bcf's tree code to get it); shape is open, the two interested classes may want different
  grids. My read: settle the cut-grid question first; its own grid means not pre-RC.
- D2, the approximate-Polya-Gamma posture ([TODO:296](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/TODO#L296) + [TODO:498](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/TODO#L498), one decision). DECIDE: does dbarts
  ship its first documented-approximate family (truncated gamma-sum, bias ~1e-3 at K = 200)
  behind an explicit opt-in, or keep refusing real dispersion and real binary weights? Both
  refusal sites say the relaxation needs no signature change, but the TODO's "explicit opt-in"
  framing would ADD a formal - the part wanting the window. My read: answer the posture now
  even if the code is post-1.0; it decides whether that formal is in the 1.0 face.
- D3, hurdle samplerOnly ([docs/design/multinomial-mutation-arc.md:581](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/docs/design/multinomial-mutation-arc.md#L581), NOT IN TODO). DECIDE: may a
  two-sampler hurdle $fit be returned unrun? Refused today while ordinal and nbinom opt back
  in via allow.samplerOnly ([R/bart.R:597-616](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/R/bart.R#L597-L616)). My read: keep refusing and say so in the Rd - a
  pair of samplers has no single unrun object.

## 3. Observed while verifying, beyond the TODO

- STALE PREMISE, [TODO:265-267](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/TODO#L265-L267). The sampler$data@x dimnames-drop limitation is FIXED, not
  pending: [R/bartcore.R:292-307](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/R/bartcore.R#L292-L307) captures and reinstalls dimnames (e8fcf8be, 211 commits before
  HEAD), [man/dbartsSampler-class.Rd:373](https://github.com/vdorie/dbarts/blob/e8fcf8be106fe14e33944e79342e60285c3baa8e/man/dbartsSampler-class.Rd#L373) documents the carrying rule,
  [inst/tinytest/test-sampler-predictors.R:75-100](https://github.com/vdorie/dbarts/blob/e8fcf8be106fe14e33944e79342e60285c3baa8e/inst/tinytest/test-sampler-predictors.R#L75-L100) pins it. The "fix is UNTICKETED" clause
  describes the pre-e8fcf8be world. Separate minor residual: the sparse/container branch
  ([R/bartcore.R:332-352](https://github.com/vdorie/dbarts/blob/e8fcf8be106fe14e33944e79342e60285c3baa8e/R/bartcore.R#L332-L352)) keeps the container's names and ignores names on the incoming
  argument.
- STALE PROSE, feature-matrix.md [f19] ([docs/design/feature-matrix.md:671-675](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/docs/design/feature-matrix.md#L671-L675)) and its Gaps student line ([docs/design/feature-matrix.md:1039](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/docs/design/feature-matrix.md#L1039)) still say
  a Student-t fit scores the gaussian density "with no guard, no doc, no TODO" - false since
  FX1-guard 4975c20b and FX1-channel c3af16a1; [f28] (heteroscedastic loglik) is by contrast
  still true, shortlist 1.
- STALE ROW, docs/design/INDEX.md lists nameable-calibration.md as PARTIAL with the flat-C
  half "riding the dbarts.h reshape's S1"; that S1 landed (ab3aa2fa) and feature-matrix.md
  calls the arc COMPLETE; source is that doc's own Status line ([docs/design/nameable-calibration.md:3-4](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/docs/design/nameable-calibration.md#L3-L4)).
- STALE ANCHORS, [docs/plans/bcf-bartcause-relocation.md:1255-1295](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/docs/plans/bcf-bartcause-relocation.md#L1255-L1295) cites [TODO:80-143](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/TODO#L80-L143) and [TODO:497-500](https://github.com/vdorie/dbarts/blob/ee70e6f7b1926a73f888540c52362b13ff863df6/TODO#L497-L500) for
  entries no longer there (80 is blocked-jacobi-trees, 498 is weighted-binary) and names a
  bart2-argument-consolidation TODO entry that no longer exists (arc closed 2026-08-17). Two
  of its doors are answered too: formula-path bases subsetting landed (S13 172523e6) and
  per-forest saved-tree replay landed (63df524e), the R half shipped as predict(type =
  "forest") ([R/generics.R:505](https://github.com/vdorie/dbarts/blob/63df524e8b54948c6bdca1820a949f2c5eb17385/R/generics.R#L505)).
- The multinomial arc's five doors have no TODO entry. Doors 1 and 5 are shortlisted; Door 2
  (flat ABI) carries the standing not-pre-1.0 call; Door 3 (nbinom's per-kept-draw bartcoreRun
  loop, [R/bart.R:2074-2076](https://github.com/vdorie/dbarts/blob/63df524e8b54948c6bdca1820a949f2c5eb17385/R/bart.R#L2074-L2076)) is perf-only, POST-1.0, and owes bench-sampler; Door 4
  ($getLatents) is declined in writing.
- Forest 1's sd and update.amplitude have no formula-term spelling
  ([docs/plans/bart2-argument-consolidation.md:1935-1937](https://github.com/vdorie/dbarts/blob/63df524e8b54948c6bdca1820a949f2c5eb17385/docs/plans/bart2-argument-consolidation.md#L1935-L1937)), the recorded cost of the XOR against a forests
  = formal. Additive; POST-1.0 unless VD wants the 1.0 face symmetric across forests.
- NOT a gap, contrary to feature-matrix [f36]: grow-from-root IS multi-forest-aware -
  growForestFromRoot loops forests_ and routes each through the combiner
  ([src/bartcore/chain.hpp:1907-1960](https://github.com/vdorie/dbarts/blob/63df524e8b54948c6bdca1820a949f2c5eb17385/src/bartcore/chain.hpp#L1907-L1960)). What two forests lack is test coverage.
