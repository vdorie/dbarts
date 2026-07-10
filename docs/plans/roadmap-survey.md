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
