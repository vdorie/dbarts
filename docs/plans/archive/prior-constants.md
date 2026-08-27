# prior-constants

agent: sonnet
rng: neutral (documentation; one optional argument)
budget: ~150 lines, mostly Rd and one design note

## Goal

Every adopted default constant has recorded provenance and, where one
exists, a derivation - in the man pages users read and one short design
note. The tau-steps = n.thin coupling gets an explicit rationale or a
dedicated argument.

## Context

- Unargued defaults on the new surface (audit, 2026-07-06): k = 2,
  power = 2, base = 0.95 (CGM/BayesTree), sigdf = 3 / sigquant = 0.9
  (CGM 2010; src/bartcore/model.hpp:1613-1615 attributes to "the
  classic engine"), node.scale 3.0 probit (bare anchor; logistic
  pi*sqrt(3) is derived FROM it, R/dbarts.R:295-302), n.trees 200
  (BayesTree) / 75 (dbarts historical), dart update delay = half
  burn-in (BART package startdart convention, R/model.R:123-127).
- Grouped tau slice sampling takes exactly n.thin steps because the
  rbart_vi R loop did (docs/design/grouped-random-effects.md, "the R
  loop's coupling").
- man/bart.Rd:336-345 cites CGM 2006/2009 generically, tying no
  constant to a source.

## Constraints

- No default changes. This item documents; changing a default is a
  posterior-changing item needing its own plan.
- Derivations only where honest: k = 2 has the CGM 95%-interval
  argument; power/base are empirical CGM recommendations - say so
  rather than inventing a derivation.
- Out of scope: range-scaling (its own decision plan).

## Steps

1. docs/design/prior-defaults.md (short): one entry per constant -
   value, source (paper section or package), derivation or "empirical
   recommendation", what it interacts with (k with node.scale and the
   response scaling).
2. man pages: one provenance sentence per parameter where users choose
   values (bart.Rd, dbarts.Rd, dbartsPriors.Rd), citing the note.
3. tau steps: add a documented rationale line to
   grouped-random-effects.md, or (if VD prefers) an explicit steps
   argument on the internal grouping attribute defaulting to n.thin -
   ask via the review, default to documentation only.

## Verification

- R CMD check: Rd cross-references and codoc clean.
- No test or draw changes (git diff shows only .Rd/.md and comments).

## Landing note (2026-07-07, 6d4f7d2)

Landed: docs/design/prior-defaults.md records every default's source
and interactions (k with the CGM variance derivation; power/base;
sigdf/sigquant; node.scale incl. the logistic pi/sqrt(3) widening;
n.trees 75 vs 200; dart update.delay; tau slice steps documented as
n.thin-reused in grouped-random-effects.md - the plan's
documentation-only default, no new argument). Absorbed the signed-off
range-scaling documentation: lineage rationale, outlier caveat with
the bartMachine log/winsorize workaround, the chi(1.25, Inf)
hyperprior escape, and the updateScale locking semantics, mirrored as
a Response scaling subsection in man/dbarts.Rd with provenance
sentences in bart.Rd/dbartsControl.Rd/dbartsPriors.Rd. Docs/Rd only.
Gates: tinytest 2470 ok, codoc clean, lint zero. Range-scaling's
step-2 closure is complete with this landing.
