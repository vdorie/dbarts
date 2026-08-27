# roadmap-survey

agent: one opus surveyor with web access; fable adjudicates
rng: neutral (findings only; a roadmap note, no code)
budget: one survey pass -> a ranked roadmap note; decisions stay
        VD's

Review 5 of the retrospective program (retrospective-reviews.md),
run last (needs sustained web access).

## Goal

Fifteen years of BART method development and the competing package
ecosystem, held against dbarts's shipped feature set and backlog:
which consumer-gated TODO items have real demand evidence, what is
missing from the backlog entirely, and what the field has moved on
from. The deliverable ranks; VD decides.

## Method

One Opus surveyor, three sweeps:
1. PACKAGES: the feature sets, defaults, and apparent user bases of
   BART (the CRAN survival-capable one), bartMachine, SoftBart,
   stochtree, bartz, pymc-bart, and dbarts's own reverse
   dependencies (bartCause, stan4bart, embarcadero, others found on
   CRAN). While there: record every package's BINARY-outcome
   leaf-prior scaling default (fixed k vs hyperprior, values, any
   published rationale) - this fulfills part 1 of the staged
   chi-default-research item.
2. LITERATURE: the major post-2010 BART method families and their
   traction (citations, implementations): sparsity/DART, soft
   trees, XBART/grow-from-root, causal (BCF and successors),
   monotone/shape-constrained, heteroscedastic, multinomial,
   SURVIVAL (notably absent from the dbarts backlog), longitudinal/
   random effects, targeted smoothing/GP leaves, model-trees/linear
   leaves, variable-selection inference.
3. DEMAND: dbarts GitHub issues (feature requests), reverse-dep
   usage patterns, what questions recur on forums; which backlog
   items (TODO: negative-binomial, ordinal-outcomes, multi-forest-
   models, monotone-bart, interaction-constraints, gp-followups,
   group-by-exposure) have external pull vs are speculative.

Deliverable: a roadmap note ranking the backlog by demand evidence,
flagging genuine gaps (features the ecosystem has that dbarts and
its backlog both lack), and noting where dbarts already leads
(mutable-data embedding, linear/GP leaves, in-core grouped).

## Constraints

- Findings only; web reads are fine, no code changes.
- One review at a time; the surveyor is the whole fan-out.
  FOREGROUND only, no sub-agents.
- Ranking is evidence-based (name the evidence per item); the
  note recommends, VD decides.

## Verification

- The roadmap note recorded here (or as docs/design/roadmap-2026.md
  if it outgrows a Status section); chi-default-research part 1
  results recorded on that TODO item.

## Status

- 2026-07-09: plan authored; surveyor dispatched.
- 2026-07-10: survey COMPLETE (citation counts via OpenAlex -
  Semantic Scholar was rate-limited; absolute numbers run low but
  ordering is reliable; a few package doc pages 403'd and were
  recovered from source). Review 5 done; the program's decisions
  stay VD's.

PACKAGES vs dbarts: BART (survival suite, multinomial, 3062 dl/mo),
bartMachine (variable-selection thresholds, 3220), SoftBart (soft
trees, 803), stochtree (XBART grow-from-root, heteroscedastic
variance forests, R+Python, serialization, inner-loop custom-
sampler API - the fast-growing one to watch, 768), bartz (GPU JAX,
n=100k demos), pymc-bart (PPL composition, most-starred at 153).
dbarts is the CRAN LEADER at 6181 dl/mo (~2x BART/bartMachine) and
the de facto engine layer: 26 CRAN reverse deps incl. stan4bart,
bartCause (558/mo), riAFTBART (581/mo), tidytreatment, voi
(656/mo), countSTAR (317/mo).

BINARY k DEFAULTS (chi-default-research part 1, fulfilled): dbarts
is the ONLY package whose binary default is a hyperprior on k
(chi(1.5, Inf), unpublished simulation origin). The field: fixed
k=2 (BART, bartMachine, bartz), SoftBart alone at fixed k=1
(unexplained, paired with soft trees + 20 trees), stochtree adapts
leaf VARIANCE via a sampled IG (not k, and for continuous too),
pymc-bart exposes no knob. No published rationale anywhere for a
non-2 binary default. Part 2's simulation breaks new ground; the
conservative match-the-field answer is fixed k=2.

LITERATURE (traction order): causal/BCF (dominant - SHIPS),
SURVIVAL (Sparapani 2016 out-cites DART; suites in BART, nftbart
309/mo, CDC BART-Survival with JOSS paper - GAP), sparsity/DART
(ships), variable-selection inference (Bleich 2014, 145 cites -
gap: dbarts has varcounts but no inference layer), soft trees,
longitudinal/random effects (dbarts ships intercepts; slopes/
multi-RE are gaps), heteroscedastic (stochtree ships it),
monotone (thin), XBART/grow-from-root (memo landed), multinomial
(modest), GP smoothing (SHIPS), model-tree leaves (SHIPS).

DEMAND-RANKED BACKLOG:
1. Random-effects breadth - the largest persistent issue cluster
   (#17 new-group prediction, #20/#23 slopes + multiple REs, #18,
   #63, #46); the asks go BEYOND group.by exposure - scope any
   memo to slopes/multi-RE/new-group prediction.
2. interaction-constraints - two independent named asks (#9
   ecology w/ xgboost parity, #67 wants setSplitProbabilities for
   Linero hyperpriors) on the designed SplitSelector seam; the
   consumer gate is arguably already satisfied.
3. prior-predictive - #31 is the most-discussed feature thread.
4. weighted-binary - #79 is a filed CORRECTNESS complaint (binary
   weights interpreted inconsistently), not a wish; resolve the
   semantics regardless of the full real-weight PG sampler.
5. pointwise-loglik - #25/#37 diagnostics demand; loo/WAIC is
   table stakes now; cheap Results field.
6. negative-binomial (indirect pull via countSTAR),
7. ordinal-outcomes (NO mainstream implementation exists anywhere
   - a lead-not-follow opportunity, zero current pull),
8. multi-forest-models (competitive vs stochtree, no user pull,
   gated on forest-combiner), 9. monotone/gp-followups (parked).

GAPS THE BACKLOG ENUMERATION MISSED:
- SURVIVAL: the strongest missing candidate; docs/plans/
  survival-models.md ALREADY EXISTS in the repo (AFT log-normal +
  discrete-time hazard, riAFTBART named as the landing consumer) -
  the gap was in this review's enumerated list, not VD's thinking.
  Recommendation: promote survival-models to first tier; AFT
  log-normal first (probit truncated-normal machinery reused;
  uncensored reduces exactly to gaussian = a free gate);
  riAFTBART (581/mo, currently building survival ON TOP of
  dbarts) as the consumer.
- Variable-selection inference layer (permutation thresholds /
  posterior inclusion; bartMachine and revdep embarcadero both
  provide it): small, R-level, high visibility.
- Convergence diagnostics (#25): every competitor bolts this on.

WHERE DBARTS LEADS: all three leaf models in one engine
(constant+linear+GP - unmatched); the INTEGRATION (BCF + MIA +
grouped RE + per-obs mutation in one engine) is the moat; CRAN
reach. CAVEAT: the embedded-conditional-model niche is now
CONTESTED (stochtree inner-loop API, SoftBart general_, pymc-bart
PPL) - dbarts's mutable-data path is the most mature and
battle-tested, but frame it as leading, not uncontested.
