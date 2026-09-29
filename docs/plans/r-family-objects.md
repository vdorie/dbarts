# r-family-objects: accept base R's family objects

agent: sonnet
rng: neutral (a mapped family builds the same sampler as its dbarts spelling)
budget: ~200 lines (R ~60, tests ~100, man and NEWS ~40)

Status: LANDED 2026-09-28 (dd383e55, 7f04b6ba)

## Goal

Every entry point's `family` argument accepts base R's family spellings as
`glm` does, maps the ones dbarts has a model for, and refuses the rest by
name (dec-B134).

## Mapping

- `gaussian(link = "identity")` -> gaussian; any other gaussian link refused.
- `binomial(link = "probit")` -> probit; `binomial(link = "logit")` ->
  logistic; any other binomial link refused ("binomial(link = \"cloglog\")
  is not supported; dbarts fits the probit and logit links").
- The word `"binomial"` and the bare function `binomial` mean
  `binomial()`, logit, as in `glm`. The word `"gaussian"` keeps its current
  meaning.
- Every other `stats` family, quasi family or foreign `"family"` object
  (poisson, Gamma, quasibinomial, MASS::negative.binomial, ...) is refused
  naming the family and pointing to ?dbartsFamilies.
- The resolved object is the ordinary dbarts family, so `family(fit)`,
  `update()` and printing are unchanged.

## Context

`resolveFamily` in R/family.R evaluates the argument in the family
vocabulary (dbarts's constructors shadow `gaussian` there, so a bare
`gaussian()` is already dbarts's); `binomial`, `stats::gaussian()` and the
word `"binomial"` reach its refusal ("'family' must be a family name or a
family object"), which is misleading for an R family object. Check every
entry point that takes `family` (bart, dbarts, dbartsSpec, xbart, bartBT if
it has one, rbart_vi if it takes one) and route them through the mapping;
xbart has its own token list.

## Steps

1. A mapping helper beside `resolveFamily`: a value inheriting "family"
   (or a function returning one, i.e. `binomial`) maps by `$family` and
   `$link`; the word "binomial" maps before `match.arg`. Refusals name the
   family and link.
2. The entry point's own admissible list still applies after mapping (a
   mapped logistic on an entry point without logistic is refused as
   logistic is today).
3. man/dbartsFamilies.Rd and the `family` argument of each entry point:
   one paragraph on R's family objects, the logit default of `binomial()`
   against dbarts's probit default for a 0/1 response, and the refusal.
4. NEWS: fold into the existing family-objects item (the family argument
   is new in 1.0-0 on these entry points; dec-B128).
5. Tests (inst/tinytest/test-family-r-objects.R): each accepted spelling
   gives the same sampler settings as its dbarts spelling (compare the
   resolved family object and, on one case, seeded draws); each refused
   family and link names itself; the word and function forms of binomial.

## Verification

Full tinytest (unwrapped); lint gates per CLAUDE.local.md; R CMD check
--as-cran --no-manual; equivalence unaffected (neutral).

## Agent-made calls

`"binomial"` and `binomial` follow glm's logit default; foreign family
objects with a `$family` dbarts could fit under another name (e.g.
MASS::negative.binomial with integer theta) are refused rather than mapped.

## Landing note (2026-09-28)

LANDED-pending. `resolvedFamily` in R/family.R maps base R's family objects
inside `resolveFamily`, so bart, dbarts, dbartsSpec and xbart take them;
bartBT and rbart_vi have no `family` argument. gaussian(identity),
binomial(probit) and binomial(logit) map to gaussian, probit and logistic;
the word "binomial" and the bare function `binomial` are logit; other
links and every other family or foreign "family" object are refused by
name. Manual pages (dbartsFamilies, bart, dbarts, xbart) and the existing
NEWS family-objects item carry the paragraph.

Gates, macOS arm64, private library: full tinytest unwrapped 9565 pass, 0
fail; lintr zero lints; air format clean; check-rc-codoc, check-win-drift
and check-doc-freshness OK; R CMD check --as-cran --no-manual --no-tests
on a tarball built with --no-build-vignettes: no new WARNING or NOTE (the
vignette warnings follow from skipping the build; the Date note was there
before). Sister packages: no call passes a stats family object into
dbarts (bartCause, treatSens use binomial with glm/glmer only; stan4bart's
test uses it with glmer; bairrtt none).

Lines: R 60, tests 142, man and NEWS 21.
