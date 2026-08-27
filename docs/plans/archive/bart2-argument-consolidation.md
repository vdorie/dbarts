# bart2 argument consolidation

Status: COMMITTED SPEC (2026-08-16). **ALL EIGHT FORKS DECIDED** (VD
walkthrough 2026-08-16) plus the term token; section 8 records every
resolution. Implementation is slices S1-S14 (section 7). Session evidence
(three censuses, three blind critique rounds, the conventions survey, four
verification rounds, the walkthrough decision log) is untracked session
evidence; everything load-bearing
from it is restated here. Section 8 is
Resolutions. Section 5 is fork 6's committed spec, rewritten under scoped
critique round 3 (3 BLOCKER / 12 MAJOR / 10 MINOR / 2 COMPLIANCE), whose
findings are accepted in full - including both compliance findings.

Revision history: rev 1 draft; rev 2 under critique r1; rev 3 under critique r2
(which refuted rev 2's gating route and loudness level); rev 4 the walkthrough
fold-in; rev 5 under scoped critique r3 on the fork-6 spec, S10, and 3.c.5.

Measured against branch `bartcore`, tip **4729313c**, built into a private lib:

```
R CMD INSTALL --preclean -l <private-lib> .
R_LIBS=<private-lib> Rscript <probe>
```

**MEASURED(tip)** = executed against that build. **MEASURED(r3)** = the round-3
critic's independent build, re-executed here before adoption. **MEASURED(rev5)**
= new this revision. **NEEDS MEASUREMENT** = unsettleable before the build.

**Rev 7 (this revision) closes scoped verification round 4** (0 blockers, 0
compliance failures, 5 MAJOR, 7 MINOR; the 10.16 amendment passed its full
eight-receipt audit and compliance with the decisions file scored 12/12). Every
finding was additive spec text over the term delta, and all of it is folded in:
the AST walk now covers the LHS (M-1) and recurses into the non-forest operand
(M-4); a factor member inside a compound operand is refused rather than
silently coerced to integer codes (M-2); the `:`'s left operand gains a closed
grammar that refuses multi-way colon chains (M-3); 5.3 step 4's contradiction
is resolved in favour of "the symbolic slot names a subset of the design, and a
name outside it is a refusal" (M-5); 5.6 grows from nine refusals to twelve;
and the S12 block-C mapping, the `bartcore.bcf` comparand, the dot-separated
indicator names, the `y ~ 0` residue, the ancestor-chain check, `-`'s place in
item 5, and the `x.train` quotation are corrected. No decision is reopened and
no settled section is touched.

**Rev 6: the term token is DECIDED** (VD
2026-08-16, the last open design question) and folded in - head `forest()` with
colon sugar, an AST-walk extraction, the term-context symbolic-vars grammar,
the refusal matrix it implies, S12's three-block test matrix, and Resolution 7.
Section 5.9 is now EMPTY. One receipt-driven amendment inside the decision: the
desugar TARGET expression is `basis = ~ z` / `~ cbind(a, b)` rather than the
decision file's `~0+z` / `~0+a+b`, because the shipped `basis =` channel
evaluates an EXPRESSION rather than building a model matrix - taken literally,
`~0+zf` ERRORS and `~0+a+b` collapses to one summed column (5.2.1, MEASURED).
The decided SEMANTICS are unchanged. Settled sections are untouched.

Rev 5 changes, one line each: the `forests =` formal is **removed from every
site** - VD decided an XOR and the end state is back to **54 formals** (10.11);
5.8 is rewritten on a new measurement that answers the prior question the
shift-attribution debate rested on (5.8.4); the packaged forest margin moves to
the TRAILING position to match both in-package precedents (5.8.1); T-E's
tolerance becomes `< 1e-12` (never "bitwise" for a float sum); `terms()` gains
`data =`; three new refusals join 5.6; S10's `family` default becomes `"auto"`
and loses `"hazard"`; S13 is reclassified as behavior-changing; S12's A/B
matrix is enumerated; and 3.c.5's one consumer-grounded clause is regrounded.

---

## 0. Owner constraints (standing, cited)

**VD direction 2026-08-16.** "bart and bart2 should have the same
functionality, but bart should reproduce BayesTree's defaults while bart2
should be our modern recommendations."

**VD direction 2026-08-16, amendment.** Capability parity is the historical
posture, not a requirement: each gap gets a disposition. bart-vs-bart2 default
divergence is INTENTIONAL; bart2-vs-xbart-vs-rbart_vi divergence is ACCIDENTAL
and is the reconciliation target.

**VD walkthrough 2026-08-16: all eight forks decided** (section 8). Two
directions bind this document beyond the decisions:

- **Fork 3's justification may not cite consumer behavior.** VD owns bartCause
  and rejects consumer behavior as a design input or precedent. The loudness
  level is grounded in R standards (3.c.5); wrapper consequences live in the
  migration map only (6.4).
- **Fork 7 was delegated to orchestrator judgment** and decided on the read
  recorded in section 8.

**Fork 6 is an XOR.** The door it closes was recorded as one -
`docs/plans/archive/bcf-bartcause-relocation.md`:1252-1258, "whether `bart2` reaches a
multi-forest model through a flat `forests =` formal **or** through
formula-level term syntax" - and VD picked formula terms. Rev 4 recommended
both and, worse, implemented the undecided branch in the counts, the signature,
the invariants, the migration map and two slices. That is reverted (10.11); the
recommendation survives as one door sentence in 5.5.3 and nothing downstream
depends on it.

**VD mandate 2026-08-15** (`TODO`:80-145): bart2's names and grouping lock at
1.0-0.

**Standing rules.** Sister packages migrate in lockstep; compat is a migration
cost to enumerate, never a design input. Gate an item only on failure to name
value it might enable.

**Freeze scope.** SPELLINGS and DEFAULTS only; it says nothing about
diagnostics (`bart()` already errors where BayesTree did not).

**Monotonicity rule.** *Diagnostics only get louder.* 3.c.2 audits it.

---

## 1. Problem, goals, non-goals

### 1.1 Problem

`bart2` carries 52 formals (`R/bart.R`:541-607) reaching seven destinations by
five mechanisms, and grows by accretion: no single grouping model; names that
diverge from their destinations and across entry points; defaults that diverge,
several INERT; a three-name dots backdoor; four inconsistent loudness levels;
and capability holes in the front door - `linear()`/`gp()` leaves and non-chisq
residual priors unreachable, and a multi-forest fit returning an object with no
per-forest channel.

### 1.2 Goals

G1. One channel per concept: exactly one supported ARGUMENT per knob on a given
    entry point, and the same spelling across the modern entry points. An
    argument may accept several TYPES (`dart`, `k`, `test`, `family` already
    do); that is one channel with several shapes.
G2. The signature is true: no default the function does not implement, no
    argument accepted and silently ignored.
G3. Aspiration with its gap named: every engine capability reachable from
    `dbarts()` SHOULD be reachable from `bart2()`. Fork 6 closes the
    multi-forest half for forests 2..K; **forest 1's `sd` and
    `update.amplitude` remain `dbarts()`-only under terms-only** (5.7), as does
    per-forest test-basis support (5.4) and out-of-sample per-forest replay
    (5.8.6).
G4. Growth is absorbed by existing vocabularies where the concept belongs to
    one, not by new flat formals per knob.
G5. The bart-vs-bart2 contract is explicit, maintained and tested.
G6. Argument-surface work is draw-neutral, and each slice names the gate that
    actually proves it (section 7).

### 1.3 Non-goals

N1. Reducing the formal count for its own sake (the end state is 54).
N2. Changing what the engine samples. No C/C++ change except two string
    literals for the `rngSeed` rename - no `dbarts.h` touch, no ABI hash move.
    **Fork 6 adds none**: the per-forest run channels already exist (5.1).
N3. Renaming `bart2` itself.
N4. Re-litigating `docs/plans/archive/interface-review.md`'s T1-T11.
N5. `bart()`'s BayesTree names and defaults; fork 7 closes CAPABILITY gaps only,
    by appending formals whose defaults reproduce today's behavior exactly.
N6. Refusing arguments inert because of another argument's value. Only
    inertness by RESPONSE FAMILY is in scope.
N7. Changing any existing refusal's severity.
N8. **Out-of-sample per-forest prediction** and **per-forest test bases**: both
    need engine-side doors (5.4, 5.8.6).
N9. **A `forests =` formal on bart2.** VD chose formula terms; the alternative
    is recorded as a door (5.5.3) and nothing in this design depends on it.

---

## 2. The current surface, in one page

**bart2**: 52 formals, `R/bart.R`:541-607. Destinations: 15 on a
`dbartsControl` SLOT (13 by identical name, `seed`->`rngSeed`, plus `...`), 5
on a control ATTRIBUTE, 9 on `dbartsData`, 15 on a model/prior object, 4 on a
sampler method, 17 inline.

**Reconstruction** (`R/bart.R`:655-690): an unknown-argument check over
`argNames <- names(matchedCall)[-1L]` (:655-665), `redirectCall` to
`dbartsControl` (:667), a backfill for the 13 shared names (:668-680), the hand
rename `seed -> rngSeed` (:681-683), three `%/% n.thin` divisions (:687-689).
`redirectCall` (`R/utility.R`:146-176) keeps an argument iff its name matches a
formal of the target exactly and drops the rest silently.

**Dispatch order inside bart2** (load-bearing for 3.c and 5.3): `argNames` is
captured at :655, BEFORE the multinomial (:707), ordinal (:906), nbinom (:935)
and hurdle (:965) branches and before the standard path's build (:1005-1024).

**Siblings**: `bart` 33; `xbart` 28 with `control =` at `R/xbart.R`:23 and six
slots overwritten (:52-57); `rbart_vi` 42; `dbarts` 31; `dbartsControl` 16;
`dbartsData` 10; `dbartsSpec` 21. MEASURED(tip), re-derived by critique r3.

**Vocabularies**: `forest()` 10 knobs (`R/model.R`:1547-1557; "Every knob is
per forest, so the fitting functions grow exactly one argument however many
forests a model has", :1532-1534); `interactions()` 3; `blocks()` 2;
`dbartsPriors`; `dbartsResidDists`.

**Doc surface**: `man/bart.Rd` (584 lines) documents bart AND bart2 in one
`\usage` and one `\arguments` (58 items from :146). `_pkgdown.yml`:16-19 has no
bart2 entry.

**Consumers.** bartCause is the only bart2 caller; stan4bart forwards
`bart_args` against `formals(dbartsControl)`/`formals(dbartsSpec)`; treatSens
and bairrtt have zero bart2 contact. 22 CRAN revdeps, 5 at risk; 17 are SAFE
because they are `bart()`-only.

---

## 3. Design dimensions

### 3.a Grouping model - DECIDED A1 (fork 1)

Facts: 13 bart2 formals are already `dbartsControl` formal names; `n.samples`
means sweeps in bart2 and kept draws in the control (MEASURED(tip):
`bart2(n.samples = 100, n.thin = 5)` returns 20 draws); **a supplied control
carries state** - every family/structure extension attaches to the CONTROL
(`R/spec.R`:79-80, :353, :361, :380) and `setControl` copies every
`^bartcore\.` attribute forward (`R/dbarts.R`:1050-1056), so
`bart2(control = <an earlier fit's control>)` would silently inherit a variance
forest, dispersion or survival status; xbart needs three precedence mechanisms
for one `control =` (MEASURED(tip)); `forest()`'s knobs belong to the
basis/amplitude family, not the variance forest; constraints resolve against
post-expansion model-matrix columns (F-a6).

**DECIDED: A1** - flat sampling dimensions; the variance quartet collapses;
`tree.prior`/`node.prior`/`resid.prior` admitted with by-name collision
refusals; no `control =` on bart2. A2' is recorded as available-later; its true
against-ground is the stale-attribute hazard, not decidability.

### 3.b Naming scheme - DECIDED

**b1.** Rename `dbartsControl`'s `rngSeed` formal and slot to `seed` (six
user-facing surfaces already say `seed`); sites in 6.1.
**b2.** Fitting functions standardize on `sigest`; sampler constructors keep
`sigma`. So `xbart(sigma =)` -> `sigest`, nothing else moves.
**b3.** The variance forest's engine surface is four values
(`src/R_interface_bartcore.cpp`:2011-2049), there is no interaction or block
channel for it, and the two selectors resolve differently (MEASURED(tip), three
divergences). **DECIDED: `varianceForest(vars, n.trees, base, power)`**, with
the plain selector retained as the same argument's other accepted type.
**b4.** Codify the casing rule; rename nothing.
**b5.** Do not change `bart2` itself, the 17 names consumers hardcode through
`formals()`-derived filters, or `updateState`.

### 3.c Family gating and loudness - DECIDED (forks 3 and 4)

#### 3.c.1 The rule

*An argument whose only effect is on a RESPONSE FAMILY this fit is not, and
which the user supplied, is diagnosed by name.* N6 is the complement; N7 fixes
the severity floor.

#### 3.c.2 The inventory, audited against what is already loud

| argument | inert when | today | in scope? |
|---|---|---|---|
| `sigest` | any fixed-unit-scale family | SILENT (MEASURED(tip)) | YES |
| `sigdf`, `sigquant` | same | SILENT (`R/spec.R`:251-260) | YES |
| `resid.prior` (new) | same | would be SILENT (MEASURED(tip): `fixed(2)` -> value 1) | YES, with S7 |
| `dispersion` | non-nbinom | SILENT (`R/spec.R`:79-80) | YES |
| `breaks`, `max.rows` | non-hazard | SILENT (`R/dbarts.R`:460-468) | YES |
| `prior.scale` | multinomial | **ALREADY A LOUD ERROR** (MEASURED(tip)) | NO |
| `samplerOnly`, `warm.start`, `n.grow.sweeps`, `keepTrainingFits = FALSE` | alternate families | ALREADY LOUD (`R/bart.R`:511-539) | NO |
| `weights`, `subset`, `offset`, `offset.test`, `test` | multinomial/hurdle | ALREADY LOUD | NO |
| `monotone` x non-default `proposal.probs`; `student` off gaussian; `variance` off gaussian or with monotone | - | ALREADY LOUD | NO |

Residual audit item for S2: `prior.scale` under
probit/logistic/ordinal/nbinom is unmeasured.

#### 3.c.3 Route (ii) is REFUTED (recorded, not dropped)

`bart()` reaches `dbarts()` through a literal list naming `resid.prior` and
`sigma` UNCONDITIONALLY (`R/bart.R`:2417-2436), so a `resolveSamplerSpec`-side
check fires for every `bart()` fit and therefore for the 17 revdeps that are
SAFE because they are bart()-only; `buildHostSamplerCall` is called with
`sigest = sigest` unconditionally on the standard path (:1017-1022), so a
name-presence test over-fires on every bart2 probit fit and under-fires on
ordinal/nbinom; and `resolveSamplerSpec`'s `matchedCall` is `dbarts()`'s
(`R/spec.R`:15-34), so d1's forwarding would erase the suppliedness bit.

#### 3.c.4 Route (i), the DECIDED site (fork 4)

**A standalone helper called at the six family-resolution sites**, keyed on the
`matchedCall` snapshot at `R/bart.R`:655. The snapshot precedes every family
branch and the standard path; bart2's two `matchedCall` mutations (:724-738,
:2025/:2033) happen after it; and d1/S3's forwarding changes what bart2
INSTALLS, never `matchedCall`. So S2 and S3 are independent.

| path | site | family known as |
|---|---|---|
| multinomial | :707, beside `checkFamilyUnsupportedArgs` | literal `"multinomial"` |
| ordinal | :906 | `"ordinal"` |
| nbinom | :935 | `"nbinom"` |
| hurdle.lognormal | :965 (components re-enter bart2 under their own families) | `"hurdle.lognormal"` |
| standard | after `sampler <- eval(samplerCall, ...)` at :1024, BEFORE the `samplerOnly` return at :1025-1027 | `sampler$model@family` |
| rbart_vi | after its own family resolution (`R/rbart.R`:389) | its local `family` |

`bart()` never calls the helper, so its silence is preserved by construction.
Defect D1 falls out with no new forwarding.

#### 3.c.5 Loudness level - DECIDED L1 warn, grounded in R standards (fork 3)

Three receipts, MEASURED(tip) and independently re-verified by critique r3:

1. **The stats package warns and proceeds for exactly this class.**
   `lm.fit(x, y, bogus = 1)` warns `extra argument 'bogus' will be
   disregarded`; `optim(par, fn, control = list(zzz = 1))` warns `unknown names
   in control: zzz`. Both complete normally. An argument that cannot act is
   base R's warning case. (Contrast `glm(y ~ x, bogus = 1)`, which errors
   `unused argument` - the UNKNOWN-NAME class, which bart2 already errors on.)
2. **Severity policy belongs to the caller, and only for warnings.**
   `options(warn = 2)` promotes (MEASURED: "(converted from warning)");
   `suppressWarnings()` and calling handlers demote. Nothing demotes an error.
   Choosing a warning leaves both severities available to every caller;
   choosing an error removes one for everybody.
3. **The diagnostic is a CLASSED condition.** MEASURED: a condition built with
   `warningCondition(msg, class = c("dbartsFamilyGatedWarning",
   "dbartsWarning"))` is muffled by `withCallingHandlers` on either class,
   caught by `tryCatch` on either, and demoted by `suppressWarnings(classes =)`
   on either, while ordinary warnings pass. **And base R already classes
   exactly this diagnostic**: `lm.fit`'s comes from `chkDots()` and carries
   class `chkDotsWarning` (MEASURED(r3)), so receipt 3 is a second precedent
   rather than a package invention.

**Class: `dbartsFamilyGatedWarning`**, constructed as
`warningCondition(msg, class = c("dbartsFamilyGatedWarning", "dbartsWarning"))`
- `warning` and `condition` are supplied by `warningCondition` itself and must
NOT be passed (MEASURED(r3): passing them duplicates the entries). The
resulting chain is `dbartsFamilyGatedWarning, dbartsWarning, warning,
condition`. This is the package's first classed condition, so the arc sets the
convention: `dbarts` prefix + CamelCase, with a `dbartsWarning` parent so a
future second class is mutable as a family.

Dead options, with receipts: **error-always** - by receipt (1) the stats
package's own treatment of an argument it will disregard is a warning, not an
error, and by receipt (2) an error removes the caller's ability to choose the
severity while a warning preserves both; **L4 warn-if-inferred/error-if-named**
- its discriminator is corrupted at route (ii)'s chokepoint (bart2 installs
`family = "ordinal"` itself on the auto path, MEASURED), and though it is
recoverable under route (i) (10.2) it ships two severities for one mistake;
**L5 document-only** - leaves G2 unmet for the class that produced D1/D5/D7.

Emission discipline: **one warning per call**, naming the argument, the
resolved family and why the argument cannot act.

### 3.d Defaults reconciliation

**d1.** In scope: `factors`, `missing` (token vectors, resolved and forwarded)
and `proposal.probs` (bart2's `NULL` vs `dbarts()`'s named 4-vector, the
mechanism of D4). Out of scope with reasons: `resid.dist` (its default is the
bare symbol `gaussian`, resolved inside `parsePriors`'s shadow environment);
`dispersion`, `breaks`, `max.rows` (default texts already match, so the
inertness is invisible - and this disposes of the `call$breaks <- NULL`
deletion trap).
**d2.** `n.samples`: sweeps in bart2/rbart_vi, kept draws in the control,
BayesTree's in bart. Keep both, document the boundary.
**d3.** Leave the deliberate numeric divergences; change rbart_vi's `k` to
`NULL` (D5) and both `split.probs` to `NULL` (d4).
**d4.** `split.probs = NULL` is IDENTICAL by construction
(`R/model.R`:295-301, :384-385); confirmed bitwise twice.

### 3.e Dots semantics - DECIDED rejection-only (fork 5)

`...` is not a value channel. Every name in it is an error carrying a
nearest-formal suggestion and a retired-name map; `storage` and `updateState`
become bart2/rbart_vi formals; `seed` covers `rngSeed`; the accepted set
becomes exactly `formals(bart2)`. VD recorded a **stated sunset** to
deletion-equivalent, so nothing may be designed to depend on the channel. The
shared helper also fixes the unqualified `formals()` lookups
(`R/rbart.R`:62, :63; `R/bart.R`:668-676).

### 3.f xbart and rbart_vi harmonization

**f1.** Remove xbart's `control =`; add `n.cuts`, `useQuantiles`, `n.thin`,
`storage` (four, not six: `printEvery`/`printCutoffs` act only under verbose,
which xbart forces FALSE).
**f2.** Keep the grid/schedule names; add shape checks.
**f3.** `sigma` -> `sigest`.
**f4.** Add `tree.prior` to xbart under the "grid axis overrides the object"
rule `node.prior` follows (`R/xbart.R`:168-173, :190-209).
**f5.** rbart_vi gets naming/defaults/loudness only; its model-surface gap
becomes `rbart-model-surface-parity`.

### 3.g The multi-forest front door - DECIDED, specified in section 5

**Level 1 = BUILD the in-sample per-forest output channel in-arc; Level 2 =
FORMULA TERMS**, an XOR against the `forests =` formal. The token is
`forest()` with colon sugar - `y ~ x1 + x2 + z:forest(x1 + x2)` (5.2).
Out-of-sample per-forest prediction stays doored.

---

## 4. Proposed end state

### 4.1 bart2

Grouping is PRESENTATIONAL. **Position rule, every entry point: APPEND AT THE
END.** The implemented order keeps today's positions 1-51 and appends.

```r
bart2(
  formula, data, test, subset, weights, offset, offset.test = offset,

  sigest = NA_real_, sigdf = 3.0, sigquant = 0.90,   # [~] diagnosed on fixed-scale families
  k = NULL, prior.scale = NA_real_,
  power = 2.0, base = 0.95, split.probs = NULL,      # [~] was 1/num.vars, identical
  dart = FALSE,

  n.trees = 75L, n.samples = 500L, n.burn = 500L,
  n.chains = 4L, n.threads = min(dbarts::guessNumCores(), n.chains),
  combineChains = TRUE, n.cuts = 100L, useQuantiles = FALSE, n.thin = 1L,
  keepTrainingFits = TRUE, printEvery = 100L, printCutoffs = 0L, verbose = TRUE,
  keepTrees = FALSE, keepCall = TRUE, samplerOnly = FALSE, seed = NA_integer_,

  proposal.probs = c(birth_death = 0.5, swap = 0.1, change = 0.4, birth = 0.5),  # [~] was NULL
  monotone = NULL, interactions = NULL, blocks = NULL,
  variance = NULL,                              # [~] selector OR varianceForest(...)
                                                # [-] n.trees.variance, power.variance, base.variance
  keepSampler = keepTrees, warm.start = NULL, n.grow.sweeps = 0L,
  factors = c("categorical", "indicators"),     # [~] resolved and forwarded
  family = c("auto", "gaussian", "probit", "logistic", "aft", "multinomial",
             "ordinal", "nbinom", "hazard", "hazard.probit", "hazard.logistic",
             "hurdle.lognormal", "twopart"),
  missing = c("incorporate", "error"),
  resid.dist = gaussian,
  dispersion = NA_real_,                        # [~] diagnosed unless nbinom
  breaks = NULL, max.rows = 1e7,                # [~] diagnosed unless hazard

  tree.prior = NULL, node.prior = NULL, resid.prior = NULL,   # [+] appended
  storage = c("double", "single"),                            # [+] appended
  updateState = TRUE,                                         # [+] appended
  ...                                           # [~] rejection-only diagnostic channel
)
```

**54 formals** (53 named + dots): 52 - 3 (variance quartet) + 3 (prior objects)
+ 2 (`storage`, `updateState`). This is the number VD approved. Multi-forest
reach arrives through formula terms (section 5), which add no formal.

**Partial-matching regression** (MEASURED(rev5) over live `formals()`, and
identical to critique r3's independent derivation):

| entry point | abbreviation broken | was matching | broken by |
|---|---|---|---|
| bart2 | `t=` | `test` | `tree.prior` |
| bart2 | `u=` | `useQuantiles` | `updateState` |
| bart2 | `r=`, `re=`, `res=`, `resi=`, `resid=`, `resid.=` | `resid.dist` | `resid.prior` |
| xbart | `n.th=` | `n.threads` | `n.thin` |
| rbart_vi | `u=` | `useQuantiles` | `updateState` |
| bart | none | - | - (re-measured for `subset`, `storage`, `family`, and also with `breaks`/`max.rows` added: still none) |

**Ten** collisions, not rev 4's twelve: `fo=`/`for=` (which shadow `formula`,
the first formal) existed only because of the `forests =` formal rev 4 should
not have had.

### 4.2 Group objects

**`varianceForest(vars = NULL, n.trees = 40L, base = NULL, power = NULL)`** -
four knobs matching the engine surface. `vars` resolves through
`resolveVarianceColumns`; `vars = NULL` ON THE OBJECT means every column, while
`variance = NULL`/`FALSE` means no variance forest; `base`/`power` NULL still
inherits the mean forest's. A `print`/`format` method rides with it.

**Prior objects.** `tree.prior`/`node.prior`/`resid.prior` take the
`dbartsPriors` vocabulary resolved by bare name inside those arguments, so
`node.prior = linear(columns = 1:3)`, `gp(...)` and `resid.prior = fixed(1)`
become reachable.

**Collision rule:** a prior object plus a shorthand that would build it is an
error naming both - `tree.prior` with `power`/`base`/`split.probs`/`dart`;
`node.prior` with `k`/`prior.scale`; `resid.prior` with
`sigdf`/`sigquant`/`sigest`.

**Family diagnosis:** `resid.prior` under a fixed-unit-scale family joins the
3.c.2 inventory (S7); `tree.prior` is honored on every family.

### 4.3 xbart

`control =` removed; `n.cuts`, `useQuantiles`, `n.thin`, `storage`,
`tree.prior` appended; `sigma` -> `sigest`; shape checks. Net **32 formals**.

### 4.4 rbart_vi

Dots rejection-only; `storage`/`updateState` appended; `k` -> `NULL`;
`split.probs` -> `NULL`; live token defaults; shared helpers. Net **44**.

### 4.5 bart

Three appended formals (fork 7, revised under critique r3 F1/F13/F14):
`subset = NULL`, `storage = c("double", "single")`, and
`family = c("auto", "logistic", "aft")`. Net 33 + 3 = **36 formals**. Every
default reproduces today's behavior: `subset = NULL` is what :2421 hardcodes,
`storage`'s first value is the current default, and **`family`'s first token is
`"auto"`** - MEASURED(rev5), `dbarts()`'s own default is `"auto"`, `bart()`
passes no `family` today, and `dbarts(x, binary_y)` resolves to probit while
`dbarts(x, binary_y, family = "gaussian")` resolves to gaussian, so a
`"gaussian"`-first token would have converted every binary `bart()` fit. See
8.7 and slice S10.

### 4.6 The maintained contracts (tests, not prose)

`inst/tinytest/test-argument-surface.R`:

- **T-A.** Every `dbartsControl` formal is a bart2 formal, spelled identically,
  asserted over an explicit list of the 16 slots as of 1.0-0 - a freeze-time
  snapshot, not a ratchet.
- **T-B.** For every name shared by `bart2` and `dbarts`, the deparsed default
  expressions agree, except:

  | name | bart2 | dbarts | why |
  |---|---|---|---|
  | `verbose` | TRUE | FALSE | fitters announce, constructors do not |
  | `n.samples` | 500L | 800L | different roles; semantics differ (d2) |
  | `family` | 13 tokens | 12 tokens | `multinomial` is a bart2-only composition |
  | `tree.prior` | `NULL` | `cgm` | bart2's NULL means "build from the shorthands" |
  | `node.prior` | `NULL` | `normal` | same |
  | `resid.prior` | `NULL` | `chisq` | same |

- **T-C.** The bart -> bart2 concept map (8.7's table as a data.frame) is total.
- **T-D.** One diagnosis test per row of 3.c.2: the classed warning fires under
  a family that ignores the argument and not under one that uses it; plus
  monotonicity assertions that the nine already-loud refusals still error.
- **T-E (fork 6).** The per-forest reconstruction identity (5.8.5) holds to
  `< 1e-12` - the tolerance the package's own assertion of this identity uses
  (`inst/tinytest/test-bcf-reporting.R`:64-76) - on BOTH a binary basis and a
  3-level-factor basis. Not "bitwise": the engine associates the sum
  differently from any R-side re-derivation, and MEASURED residuals are
  ~1e-15 with one chain's `cor` at 0.99999999999999967.

### 4.7 Honest arithmetic

52 -> 54. What consolidates is the number of CHANNELS.

---

## 5. Multi-forest front door (fork 6, committed spec)

Rewritten under scoped critique r3, whose central verdict was that 5.1's
load-bearing measurement CONFIRMS but that the spec built on it needed work.
Every claim is MEASURED or flagged.

### 5.1 What already exists, and the gate's actual predicate

MEASURED(tip), confirmed independently by critique r3 including K > 2 and
ragged-glue extensions: `run()` already emits both per-draw channels.

```
names: sigma, train, test, varcount, k, varprobs, tau, ranef, forestFits, glue
  train      n x n.samples x n.chains
  forestFits n x n.forests x n.samples x n.chains
  glue       sum_f q_f x n.samples x n.chains
```

Three forests give `forestFits` n x 3 x draws and `glue` 4 x draws
(q = 1 + 2 + 1); a 3-level factor basis gives `glue` 4 x draws. The
sampler-level readers `$getForestFits(forest)` (`R/dbarts.R`:1488-1491) and
`$getForestAmplitudes()` (:1499-1508) already ship.

**The gate is the COUPLING, not the forest count.**
`forestReportingIsDefined()` is true only on `BCFForestCombiner`
(`src/bartcore/combiner.hpp`:992; base false at :631; read at
`src/R_interface_bartcore.cpp`:4243), so a K-forest MULTINOMIAL carries neither
channel (`inst/tinytest/test-bcf-reporting.R`:144-146). Every reachability
statement in this section is keyed on **amplitude-coupled**, never on
"multi-forest", and S11's refusal is written "on a fit without forest
reporting", not "on a single-forest fit".

What is missing is entirely R-side: `packageBartResults` (`R/bart.R`:158) drops
both channels, so a bart2 multi-forest fit carries only `n.forests` (:327) and
a forest-widened `varcount` (:213-215, via `shapeMultinomialChannel` - anchor
corrected from rev 4's ":204"). **Level 1 is packaging plus generics, with no
engine, bridge or `dbarts.h` change.**

### 5.2 The term spelling - DECIDED (VD 2026-08-16)

**Head = `forest()`. Colon sugar is the canonical spelling.** Decided after the
conventions survey (20
packages opened, every claim measured) and two probe rounds.

```r
bart2(y ~ x1 + x2 + z:forest(x1 + x2), data)             # canonical
bart2(y ~ x1 + x2 + factor(z):forest(x1 + x2), data)     # canonical
bart2(y ~ x1 + x2 + forest(x1 + x2, basis = ~ z), data)  # general named form
```

Why this shape, from the survey: the head names the LEARNER and never the
modulation (survey F1 - `s`, `gp`, `bbs`, `btree`, `pspline`, `vc`, `lf`/`af`
all do this and nothing in 20 packages puts the modulation in the head); the
multiplier is spelled in a NAMED slot by all five varying-coefficient packages
(`by =` in mgcv, mboost, brms, gamlss, vcrpart), which `basis =` fills with
dbarts' own keyword instead of importing `by =`'s factor baggage (mgcv/brms
replicate the FUNCTION per factor level, mboost drops the reference level, and
dbarts does neither - one forest, one amplitude per level, no level dropped);
and `forest()` is textually identical to the shipped
`forests = list(forest(), forest(basis = ~z))` spelling, so the term route and
the constructor route are one vocabulary and one validator.

Rejected candidates, one line each (survey section 4): `bases(~z)` - the
`bs`/`ns` class-`basis` prior makes it read "expand z into predictor columns",
an additive-vs-multiplicative inversion, and the plural belongs to the
list-across-forests channel; `basis(~z)` - fixes the number, worsens the
homograph (it is literally `bs()`'s return class); `bart(~z)` - exported, so
`dbarts::bart(~z)` would be CALLED for real, and stan4bart's shipped
`y ~ bart(x_1 + x_2)` means the opposite thing to this exact audience;
`amplitude(~z)` - mis-names its referent (the amplitudes are the sampled
coefficients, z is what they multiply); `moderated`/`moderator` - disqualified
in-repo, since dbarts already spells the split-variable subset `moderators`;
`vc(~z)` - vcrpart exports it with the moderators in the unnamed slot, exactly
backwards; `modulate(~z)` - best fresh word, but prose vocabulary rather than
machinery vocabulary.

#### 5.2.1 The desugar rules

`:` with a `forest()` operand desugars to a basis on that forest. One semantic
channel - `forest(basis =)` - with sugar over it, the fork-2b pattern.

| sugar | desugars to | basis columns |
|---|---|---|
| `z:forest(X)`, numeric z | `forest(X, basis = ~ z)` | 1 (z) |
| `zf:forest(X)`, factor zf | `forest(X, basis = ~ zf)` | one per LEVEL, no reference dropped |
| `factor(z):forest(X)` | `forest(X, basis = ~ factor(z))` | one per level |
| `(a + b):forest(X)`, every member numeric or logical | `forest(X, basis = ~ cbind(a, b))` | 2 - ONE forest with a two-column basis |
| `forest(X):z` | same as `z:forest(X)` | operand order is immaterial |

#### 5.2.1a The left-operand grammar, and what it refuses

The `:`'s NON-forest operand needs a grammar as explicit as the unnamed slot's
(5.2.3(a)), because the desugar hands it to an EXPRESSION evaluator and R's
operators do not all mean there what they mean in a formula. **Accepted in v1:
a bare symbol; `factor(<symbol>)`; or a `(`-wrapped `+` chain of bare symbols.**
Everything else is refused by name, quoting the operand and pointing at the
general channel `forest(X, basis = ~ <expr>)`, which exists precisely for
arbitrary bases. Three shapes are singled out because they are reachable and
degenerate rather than merely unsupported:

- **A compound operand containing any NON-NUMERIC, NON-LOGICAL member** -
  `(a + zf):forest(X)` with `zf` a factor, `(a + zc):forest(X)` with `zc`
  character - is REFUSED: *"a compound modulating operand takes numeric or
  logical variables only; give the factor or character variable its own term,
  or use forest(basis = ) with an explicit basis."* MEASURED(v4, reproducing
  rev 6's own channel receipts): `~ cbind(a, zf)` yields 60x2 whose second
  column is zf's INTEGER CODES, so the compound target would silently fit an
  unordered factor as a numeric score on one amplitude instead of one amplitude
  per level; a CHARACTER member is worse still - `~ cbind(a, zc)` explodes to
  120x63 under `cbind`'s character coercion, which is neither the decided
  semantics nor a recognizable failure. **Timing, parallel to item 2's:** the
  syntactic grammar above admits `(a + zf)` - it is a `(`-wrapped `+` chain of
  bare symbols - so the member TYPES are unknowable until the variables are
  evaluated, and this refusal necessarily fires after the model frame exists,
  where the term's variables are columns with types. It is the one term refusal
  that is not decidable from the AST alone. This is also the one shape on which
  the 10.16 amendment is a REGRESSION rather than a fix - the decision's
  literal `~0+a+zf` would have errored loudly - so the refusal restores the
  loudness the amendment would otherwise have spent. The alternative (desugar
  member-wise through `expandForestBasis` and column-bind) is rejected for v1
  because its target is not an expression a user could have written, which
  would break the one-channel property the fork-2b pattern rests on; recorded
  as a door.
- **A multi-way colon chain** - `z:w:forest(X)`, which parses as
  `(z:w):forest(X)` - is REFUSED: *"modulate by one operand; combine variables
  via forest(basis = )."* MEASURED(v4): desugared naively to `~ z:w` the
  shipped evaluator applies R's SEQUENCE operator, giving a length-1 value and
  a 1x1 basis (with `numerical expression has 60 elements: only the first
  used`), and the user then meets a row-count complaint naming neither the term
  nor the nesting. Door: a future version may fold a chain into a
  column-crossed basis.
- **A `forest()` call on BOTH sides** - `forest(A):forest(B)`, and the buried
  form `z:forest(A):forest(B)` - is REFUSED by name. MEASURED(v4): under a
  one-of-two-operands predicate the inner `forest(A)` is absorbed into the
  basis expression, never walked again, never stripped, and reaches
  `evaluateForestBasis`, which calls the exported `forest()` for real and fails
  on types without naming a term. The walk therefore recurses into the
  non-forest operand as well (5.3 step 1) so a buried second head is found.

Literals (`3:forest(X)`), transformations (`log(z):forest(X)`), `-` and `.` in
the operand are refused by the same closed grammar. `log(z):forest(X)` happens
to evaluate correctly today; it is refused anyway, because accepting it by
accident would make the grammar's boundary undiscoverable.

**Correction to the decision's written target expression, with receipts.** The
decision file writes the target as `basis = ~0+z` / `~0+a+b`, which is
model.matrix idiom. The shipped `basis =` channel is not a model.matrix
channel: `evaluateForestBasis` (`R/model.R`:751-763) EVALUATES the formula's
RHS as an expression and `expandForestBasis` (:772-...) then applies the
expansion rule to the RESULT (factor -> level indicators with no reference
dropped, character -> factor -> indicators, logical -> two columns, numeric
vector/matrix -> itself). MEASURED(rev6) through the shipped channel:

```
basis = ~z            -> 60x1                 basis = ~0+z      -> 60x1  (accident of arithmetic)
basis = ~zf           -> 60x3                 basis = ~0+zf     -> ERROR "a 'basis' cannot be NA"
basis = ~factor(z)    -> 60x2                                       (Ops.factor: '+' not meaningful)
basis = ~cbind(a, b)  -> 60x2                 basis = ~0+a+b    -> 60x1  (the ELEMENTWISE SUM)
```

So the literal `~0+` target would ERROR on `factor(z):forest(X)` - one of the
two canonical spellings VD wrote - and would silently collapse
`(a+b):forest(X)` to a single summed column. The SEMANTICS the decision states
are preserved exactly (a no-intercept basis; a compound left operand folding
into one multi-column basis on one forest); only the target EXPRESSION changes,
to the one the shipped evaluator implements. Note the semantics coincide: the
shipped rule already emits no intercept column for a numeric operand and all
K indicators for a factor, which is what "no-intercept basis" means in each
case. The one shape where the amendment would have traded a loud failure for a
silent wrong answer - a factor member inside a compound operand - is closed by
5.2.1a's refusal rather than left to the evaluator.

#### 5.2.2 `*` with a forest operand is REFUSED

`z*forest(X)` errors, naming both explicit forms: *"'z * forest(...)' is not
supported: write z:forest(...) to modulate the forest by z, or
z + z:forest(...) to also include z as a main effect."* Honoring R's expansion
(`z + forest(X) + z:forest(X)`) would require an unmodulated ADDITIONAL forest
alongside the modulated one - unverified engine territory, and not what the
user is likely to mean.

#### 5.2.3 Term-context grammar: the unnamed slot is SYMBOLIC

**In term context the unnamed slot of `forest()` is a symbolic predictor set** -
`forest(x1 + x2)` means the variable set {x1, x2}, a la stan4bart's
`bart(x1 + x2)` - interpreted from the AST and **NEVER evaluated** (evaluated,
`x1 + x2` would be an elementwise sum). This DIVERGES from the constructor,
whose first positional formal is `basis`.

(a) **The grammar, v1, minimal with named refusals.** Accepted: a bare symbol,
or bare symbols joined by `+`. MEASURED(rev6) that each of the rejected shapes
is detectable in the AST before any evaluation: `log(x1) + x2`, `x1 - x2` and
`.` are all distinguishable from the accepted grammar by inspection.
Refused by name in v1: transformations (`log(x1)`), `-` removals, `.`, `:`/`*`
inside the slot, and literals. The refusal names the offending expression and
points at `vars = ` for anything the grammar cannot say. Rationale for the
minimal grammar: the slot maps onto `forest(vars =)`, which resolves against
model-matrix COLUMN NAMES (`resolveModerators`, `R/model.R`:698-721) and has no
expression vocabulary at all - a richer term grammar would promise something
the destination cannot honor.

(b) **Named arguments inside a term-context `forest()` evaluate NORMALLY**:
`n.trees`, `base`, `power`, `sd`, `interactions`, `blocks`,
`amplitude.prior.variance`, `update.amplitude`, and `basis` are ordinary
arguments, evaluated as written. Only the unnamed slot is symbolic.

(c) **Positional `basis` matching does NOT apply in term context**, and the
docs spell `basis =` by name everywhere - including in every `forests =`
example - so no reader is taught the collision. A positional SECOND unnamed
argument in a term is refused by name rather than silently matched to `vars`.

(d) **Symbolic vars and factor expansion.** `vars` resolves against
POST-EXPANSION model-matrix column names (F-a6). MEASURED(rev6, column names corrected under v4 m-1): with
`factors = "categorical"` a factor `zf` is ONE column named `zf`, so
`forest(x1 + zf)` resolves; with `factors = "indicators"` the columns are
**`zf.u, zf.v, zf.w`** - DOT-separated, and the dot is load-bearing, because
`resolveTermColumns` expands via
`which(startsWith(columnNames, paste0(name, ".")))` (`R/model.R`:527) and an
implementer working from a dotless spelling would get `integer(0)` rather than
an error (its own doc comment warns that a recognized term with no indicator
columns yields `integer(0)`, not NULL). The bare name `zf` does NOT match -
`resolveModerators` errors `'vars' name not found in the design's column
names`. So a symbolic factor name maps to **all of that factor's expanded
columns**: the term's resolver expands a symbolic name to every model-matrix
column derived from it, using the same `resolveTermColumns` term-expansion the
variance selector already uses (`R/model.R`:667-671), rather than the bare
`match()` `resolveModerators` performs today. That is a term-route resolver
detail, not a change to `forest(vars =)`'s own contract.

### 5.3 Ingestion

The rule: **a term is evaluated INSIDE the model frame, so `subset`, NA
handling and `data` scoping apply to it exactly as they do to predictors.**

Why it matters (MEASURED(tip)): today
`dbarts(y ~ a, d, subset = 1:20, forests = list(forest(), forest(basis = ~
tr)))` is refused - `length of 'basis' must equal length of 'y'` - because the
basis is evaluated against the pre-subset `data` while the model frame owns row
selection on the formula branch (`R/data.R`:895 validates without a subset;
:977 and :1043 are the two x/y sites that do pass one). Terms sidestep the
ambiguity the recorded door names.

**Extraction is an AST WALK over the formula, not `terms(specials =)`** (VD:
"you can just parse the formula"). The survey records that `terms(specials =)`
is the minority mechanism - of 20 packages only mgcv, survival and refund use
it, while lme4 and stan4bart walk the AST, mboost evaluates the RHS with a
shadowed `+`, brms regexes deparsed labels, gamlss tags attributes. The walk
also dissolves two hazards rev 5 had to design around: the qualified-head
escape (`terms(specials =)` matches only a bare symbol, so `dbarts::forest(...)`
was invisible to it) and the `.`-formula failure (`terms()` needs `data =` to
resolve a dot). MEASURED(rev6): the walk sees `y ~ .` with no `data` argument
and reports zero hits, no error.

Steps, in bart2 and in `dbarts()` (shared ingestion):

1. **Walk BOTH SIDES of the formula's AST.** A node is a term hit when it is a
   `:` or `*` call one of whose two operands is a `forest()` call - bare or
   `dbarts::`-qualified - or when it is a bare `forest()` call in a top-level
   additive position. Two coverage rules the canonical set does not exercise:
   - **The LHS is walked too.** `forest(x1) ~ x2` is a legal formula
     (MEASURED(v4): `deparse(f[[2L]])` is `"forest(x1)"`), and a walk that
     starts at the RHS would leave 5.6 item 5's "a term on the LHS" refusal -
     and S12 cell (21) - with no mechanism, so the shape would surface only as
     `model.frame`'s `invalid type (list) for variable`, a message no part of
     this spec owns. Any `forest()` found on the LHS is an immediate item-5
     refusal.
   - **The walk recurses into the non-forest operand** of a `:`/`*` hit, so a
     second, buried head is found rather than absorbed into the basis
     expression (5.2.1a's both-sides refusal).

   MEASURED(rev6) on the canonical set, all as specced:

   ```
   y ~ x1 + x2 + z:forest(x1 + x2)          -> 1 hit, op ':', forest right, other z
   y ~ x1 + x2 + factor(z):forest(x1 + x2)  -> 1 hit, op ':', other factor(z)
   y ~ x1 + forest(x1 + x2):z               -> 1 hit, op ':', forest LEFT, other z
   y ~ x1 + z*forest(x1 + x2)               -> 1 hit, op '*'  (-> the 5.2.2 refusal)
   y ~ x1 + (a + b):forest(x1 + x2)         -> 1 hit, other '(a + b)'
   y ~ x1 + z:dbarts::forest(x1 + x2)       -> 1 hit, qualified head matched
   y ~ x1 + forest(basis = ~z)              -> 1 hit, bare general form
   y ~ .                                    -> 0 hits, no error
   ```

   **Both operand orders are preserved distinct in the AST** (MEASURED(rev6):
   `z:forest(x1)` and `forest(x1):z` are non-`identical()` language objects and
   `terms()` labels them differently), so the walk reads the forest side rather
   than assuming a position; the desugared result is the same forest either way.
   Compound left operands (`(a + b)`) arrive as a single `(`-wrapped call and
   are flattened over `+` to the basis column set.
2. **Context is tracked, not just presence, and over the whole ANCESTOR CHAIN
   rather than the immediate parent.** A `forest()` call is accepted only when
   every ancestor up to the `~` is one of `+`, `:`, `*` (and the `~` itself);
   `-`, `^`, `I()` and any other call anywhere in that chain is an item-5
   refusal. MEASURED(rev6): a context-free walk descends into
   `I(z:forest(x1))` and reports a hit. MEASURED(v4): a PARENT-ONLY check also
   accepts `y ~ (x1 + z:forest(X))^2`, whose forest's immediate parent is a
   legal `:` - hence the chain, not the parent. Note the useful negative shape
   is unaffected: in `y ~ . - z + z:forest(x1)` the `-` sits in a sibling
   subtree, so the forest's chain is `:` -> `+` -> `~` and it is accepted; only
   a forest INSIDE a removal (`y ~ x1 - z:forest(X)`) is refused.
3. **Rewrite the formula by AST surgery**: delete each hit's node from the
   additive chain and keep everything else byte-for-byte. Because the rewrite
   is structural rather than a `reformulate` over `term.labels`, it preserves
   the intercept term, any in-formula `offset()`, back-ticked names and the
   environment by construction - the three things rev 5's label-based rebuild
   had to special-case (MEASURED(rev5): `offset()` lives on
   `attr(tt, "offset")` and never appears in `term.labels`; `reformulate`
   defaults to `intercept = TRUE`).
4. Build `x` from the rewritten formula through the ordinary `model.frame` /
   `dbartsData` path, unchanged. **The rule, stated once: the symbolic slot
   NAMES A SUBSET OF THE DESIGN the RHS already builds; a symbolic name that is
   not among the rewritten RHS's terms is a NAMED REFUSAL**, quoting the name
   and pointing at the RHS. Nothing is added to the design matrix by a term.
   Rev 6 stated both halves of a contradiction here; this is the resolution.

   Why this branch: `vars` resolves against `colnames(data@x)` and
   `attr(data@x, "term.labels")` - the DESIGN matrix, not the model frame
   (`resolveModerators` `R/model.R`:706, `resolveTermColumns` :521-529), so a
   slot variable must be a design column for the slot to mean anything at all
   (F-a6). The alternative - silently adding the slot's variables to `data@x` -
   would make `y ~ x2 + z:forest(x1 + x2)` give forest 1 a predictor `x1` the
   user never put on the RHS, which is a side effect of declaring forest 2 and
   the opposite of what the term is for. Under the chosen rule that call is
   refused by name and the user writes `y ~ x1 + x2 + z:forest(x1 + x2)`, the
   canonical spelling, in which the slot is exactly a restriction of forest 2's
   split set. Steps 1, 2 and 5 are unaffected; `subset`/NA handling reaches the
   slot's variables because they are RHS predictors.
5. Evaluate each basis against the MODEL FRAME (post-subset, post-NA) through
   `evaluateForestBasis`/`expandForestBasis` (`R/model.R`:751-763, :772-...),
   whose rule - a factor expands to its level indicators with NO reference
   level dropped - is the one `setForestBasis` documents and the one 5.2.1's
   desugar targets are written against. Raw `model.matrix(~ z)` is NOT the rule
   (MEASURED(tip): it emits `(Intercept), zv` where the engine wants
   `zu, zv`).
6. Hand the expanded matrices to `dbartsData(bases = )`, which on the formula
   branch validates against the post-subset row count (`R/data.R`:895), and the
   per-forest knobs to `forests =` as `forest()` objects.

**The un-stripped error path.** If a `forest()` call ever reaches
`model.frame` - i.e. the walk missed it, or a consumer evaluates the stored
formula itself - the failure is loud, not silent: the exported `forest()` is
called for real and returns a `dbartsForest` list, and `model.frame` errors
`invalid type (list) for variable 'forest(...)'`, naming the term
(MEASURED(tip)). That is the message-quality cost of an exported head, and it
is strictly better than `bart(...)`'s `argument "x.train" is missing, with no
default` (the formal corrected under v4 m-7). It is not, however, a substitute
for a refusal: 5.6's promise is that every refusal names the term AND what it
wanted, which is why the LHS is walked (step 1) rather than left to this path.

**NEEDS MEASUREMENT (build-time):** that a term-built fit is bitwise identical
to the equivalent hand-written `forest(X, basis = ~ z)` call, and to the
`dbartsData(bases =)` + `forests = list(...)` route, at the same seed. This is
S12 block B, and byte-identity is the right bar there (unlike T-E's float sum)
because the sugar rewrites to the named call before anything numeric runs.

### 5.4 subset, weights, offset, test, missing

- **`subset`**: handled by construction (5.3 step 3).
- **`weights` / `offset`**: untouched as model-frame specials, provided step 3
  preserves the `offset` attribute.
- **`test`: REFUSED by name, decided by measurement, not open.** MEASURED(r3
  and rev5): `dbartsData` has no test-basis channel at all, and an
  amplitude-coupled fit refuses test predictors as a WHOLE-MODEL by-name error
  - `a treatment forest does not support test predictors; drop it or fit a
  single-forest model` - consistent with
  `BCFForestCombiner::testFitsAreDefined() -> false`
  (`src/bartcore/combiner.hpp`:986) and the bridge refusal at
  `src/R_interface_bartcore.cpp`:2833-2836. So a term-bearing formula with
  `test` is refused, restating that message in the term's spelling. Supporting
  test bases is the engine-side door
  `docs/plans/archive/bcf-bartcause-relocation.md`:1274-1276 records as "A test
  treatment vector / test basis [...] A modelling decision first" - not cheap,
  and not this arc's. S12 loses the test-frame work rev 4 budgeted.
- **`missing = "incorporate"`**: an NA in an expanded basis has no defined
  multiplier; `validateForestBases`'s finiteness check (`R/data.R`:662-664)
  already refuses it and the term route inherits that.

### 5.5 How terms reach the forests machinery

#### 5.5.1 The two halves

A multi-forest fit needs the bases (conditioning DATA, riding `data@bases`) and
the per-forest SPEC (tree counts, priors, constraints, riding
`forests = list(forest(), ...)` into `resolveForests`, `R/spec.R`:193-215).
Terms supply the first directly and expand to a `forest()` object internally
for the second, so the existing resolution path is unchanged.

#### 5.5.2 Knobs on a term

**A term's named arguments are exactly `forest()`'s** - `n.trees`, `base`,
`power`, `sd`, `interactions`, `blocks`, `amplitude.prior.variance`,
`update.amplitude`, plus `basis` in the general form - evaluated normally
(5.2.3(b)). One vocabulary, one validator, closed by construction: a named
argument outside `forest()`'s formals is refused by `forest()` itself.

The unnamed slot supplies **`vars`**, symbolically (5.2.3). So `vars` has two
spellings in term context - the symbolic slot and `vars = ` by name - and
giving both is a collision refusal, as is giving both the `:` sugar and
`basis = ` (each pair is two spellings of one knob on one call, the same shape
as 4.2's prior-object collisions).

#### 5.5.3 The `forests =` formal: NOT taken (door)

VD decided formula terms as an XOR against a flat `forests =` formal, so bart2
does not gain one and nothing in this design depends on it. **A `forests =`
formal remains addable later; the dual spelling was declined 2026-08-16.** The
strongest argument for revisiting it is not the programmatic one rev 4 gave
(a term is constructible programmatically in one line via
`do.call(bart2, c(list(formula = reformulate(...)), args))`, and
`dbartsData(bases =)` remains open) but 5.7's: **forest 1's `sd` and
`update.amplitude` have no term spelling**, because forest 1 has no basis.

#### 5.5.4 Coexistence with flat monotone/interactions/blocks - A1 stands

Constraints stay FLAT and do not move into terms. Coherent because: the flat
`interactions`/`blocks` arguments are documented as the FIRST forest's
(`R/model.R`:1541-1542; `R/spec.R`:221-226 installs them there) and forest 1
has no basis and so no term to attach to; constraints resolve against
POST-EXPANSION model-matrix columns, which a formula does not name after factor
expansion; and fork 1 chose A1, of which moving constraints into the formula
(A4) was the rejected branch. A per-forest constraint for forest f >= 2 is a
term knob (5.5.2).

### 5.6 Error behavior

Twelve refusals, each naming the term and what it wanted. S12 block C covers
every one of them except item 7, which is enforced by construction and stated
as such:

1. A basis expanding to zero columns, or to an all-zero column. Half of this is
   already enforced and half is new work: MEASURED(v4),
   `expandForestBasis(matrix(0, n, 0))` errors `a 'basis' must have at least
   one column`, but `expandForestBasis(rep(0, n))` returns an n x 1 matrix with
   NO complaint, so the all-zero column needs a refusal that does not exist
   today.
2. **A term on a family the multiplier model does not support - refused at
   family resolution, by name, for the full matrix.** MEASURED(rev5), today:
   gaussian/probit/logistic ACCEPT; aft/ordinal/nbinom are refused by name
   (`a treatment forest does not support family "aft"`); and **`hazard` is NOT
   refused by name** - it fails with `length of 'basis' must equal length of
   'y'`, because the person-period expansion changes the row count after the
   basis is built. Under the term route that message would name neither hazard,
   nor periods, nor the term. So S12 adds an explicit by-name refusal keyed on
   the resolved family covering `hazard`, `hazard.probit`, `hazard.logistic`,
   `multinomial`, `hurdle.lognormal`/`twopart`, and it fires where the family
   is resolved, not where row counts are compared. This is slice content, not a
   doc sentence.
3. A term plus a pre-built `dbartsData` as `formula` - the shape `dbarts()`
   already refuses at `R/dbarts.R`:554-565, restated for terms.
4. `test` with a term (5.4).
5. **A term appearing anywhere other than as a top-level additive term or a
   `:` operand** - refused by name, with the offending expression quoted.
   MEASURED(rev6): a context-free walk descends into `I(z:forest(x1))` and
   reports a hit, so this refusal is the walk's ancestor-chain check (5.3
   step 2), not an afterthought. Family: `I(...)`, `^`, **`-`** (a forest
   inside a removal, `y ~ x1 - z:forest(X)`), a term inside any other function
   call, and **a term on the LHS** (`forest(x1) ~ x2`), which step 1's
   both-sides walk is what makes reachable as a refusal at all.
5a. **`*` with a `forest()` operand** - refused with the two-explicit-forms
   message (5.2.2). MEASURED(rev6) that the walk distinguishes `*` from `:` at
   the node, so the refusal fires before any evaluation.
5b. **An unsupported symbolic-vars expression** in the unnamed slot -
   `log(x1)`, `x1 - x2`, `.`, `:`/`*`, a literal - refused by name, quoting the
   expression and pointing at `vars = ` (5.2.3(a)). MEASURED(rev6): each is
   distinguishable in the AST before evaluation.
5c. **A symbolic-slot name that is not a term of the rewritten RHS** - refused
   by name, quoting the name and pointing at the RHS (5.3 step 4).
5d. **An unsupported LEFT-operand expression** - anything outside the closed
   grammar of 5.2.1a: a non-numeric, non-logical member inside a compound
   operand (the one term refusal that fires after the model frame rather than
   from the AST - 5.2.1a), a multi-way
   colon chain, a literal, a transformation, `-`, `.`. Each refusal quotes the
   operand and points at `forest(X, basis = ~ <expr>)`.
5e. **A `forest()` on both sides of a `:`/`*`** - refused by name (5.2.1a),
   found by step 1's recursion into the non-forest operand.
6. **A formula whose only RHS content is a term** (`y ~ z:forest(x1)`): the
   rewrite leaves a predictor-free RHS - `y ~ 1`, or `y ~ 0` when the user
   suppressed the intercept (MEASURED(v4): the specced AST surgery turns
   `y ~ 0 + z:forest(X)` into `y ~ 0`, and preserves the intercept form on both
   no-intercept spellings, which independently confirms step 3). Refused by
   name in either shape - a forest with no predictors to split on is not a
   model the user meant. (An intercept-only forest 1 is expressible through
   `dbarts()`.)
7. A term knob outside `forest()`'s vocabulary - caught by construction
   (5.5.2), which is why rev 4's "a variance forest declared on a term" example
   is deleted: it is unspellable, so nothing could fire on it.

### 5.7 What terms do NOT close (corrected)

Rev 4 claimed terms close G3's gap for every `forest()` knob. That is wrong,
and the correction matters because it is the honest cost of the XOR: **a term
declares an ADDITIONAL forest by its basis, and forest 1 - the fit's own - has
no basis and therefore no term spelling.** Forest 1's knobs are live, not
decorative: MEASURED(rev5), a two-forest fit's draws differ at the same seed
under `forest1 sd = 0.05`, under `update.amplitude = FALSE`, and under
`base = 0.5`; `amplitude.prior.variance` on forest 1 is refused by name (it has
no basis).

Of those, `n.trees`, `base`, `power`, `interactions`, `blocks` and the
split-variable restriction are ALREADY reachable from bart2 as flat arguments
(`n.trees`, `base`, `power`, `interactions`, `blocks`) - they are the fit's own
knobs and always were. The genuinely unreachable ones under terms-only are
**forest 1's `sd` and `update.amplitude`**, plus per-forest test bases (5.4)
and out-of-sample replay (5.8.6). Recorded in G3 and in the 5.5.3 door.

### 5.8 The in-sample output channel

#### 5.8.1 Shape: the forest margin goes LAST

MEASURED(rev5), the two in-package precedents both put the widened margin
trailing: a bart2 multinomial fit's `yhat.train` is draws x n x K with
`dimnames[[3]] = p, q, r`, and the shipped forest-widened `varcount` on a
two-forest bart2 fit is draws x p x forest with
`dimnames[[3]] = forest1, forest2`. Rev 4 proposed the forest margin in the
MIDDLE, which would have put the forest axis at position 3 for `varcount` and
position 2 for `forestFits` on the same object.

| element | raw run shape | packaged (combineChains = TRUE) |
|---|---|---|
| `forestFits` | n x n.forests x n.samples x n.chains | draws x n x n.forests |
| `glue` | sum_f q_f x n.samples x n.chains | draws x sum_f q_f |

`shapeMultinomialChannel` - the helper the widened `varcount` already goes
through (`R/bart.R`:213-215) - produces the trailing form, so S11 is one
existing helper and one call rather than a new reshape.

#### 5.8.2 Margin naming

The shipped `varcount` names its forest margin `forest1..forestK`
unconditionally (MEASURED(rev5)). To avoid two vocabularies on one object,
`forestFits` uses **the same `forest1..forestK` names**, and the declaration's
own labels (a term's text, or `names(forests)` on the `dbarts()` route) ride as
a `forest.labels` attribute on the fit. The `glue` margin is RAGGED and carries
a `forest` attribute mapping each row to its forest, per
`$getForestAmplitudes()`'s documented forest-major order.

#### 5.8.3 What `forestFits` contains - the prior question, MEASURED

Critique r3 was right that rev 4's shift-attribution recommendation was
internally contradictory, and right that the prior question had never been
asked. It is now answered.

MEASURED(rev5), two-forest fit, three draws, `getCalibration()` giving
`response.scale = 8.003532`, `response.shift = 1.861575`:

```
max| forestFits[,1,s] - getForestFits(1) |            = 0          (bitwise equal)
max| forestFits[,2,s] - getForestFits(2) |            = 0          (bitwise equal)
max| forestFits[,2,s] - g_z * getForestFits(2) |      = 0.2776613  (NOT the glued form)
UNGLUED  max| shift + scale*(f1 + f2)        - train | = 2.95 / 3.89 / 3.35
GLUED    max| shift + scale*(g1*f1 + g_z*f2) - train | = 4.4e-16 / 8.9e-16 / 8.9e-16
```

**`forestFits[, k]` holds the RAW forest total `f_k`, glue-free** - bitwise
identical to `$getForestFits(k)`, the shipped R5 accessor. So the additive
identity requires the glue applied outside the element, and any recommendation
that promises `yhat = shift + sum_k element_k` for glue-free elements is false.

#### 5.8.4 The consequent design choice

The three rev-4 candidates restated against that truth, plus the honest trade:

- **(a) rev 4's "scale-multiplied contributions plus a shift attribute"** -
  REFUTED. It presupposed the element already carried its multiplier.
- **(b) attribute the shift to forest 1** - still available, but it inherits
  (a)'s defect unless the glue is folded in, and folding it in makes element 2
  carry a per-observation multiplier (range -0.518 to 1.010 in one draw), i.e.
  not a function of forest 2 alone.
- **(c) ship the raw internal scale untouched** - self-consistent but forces
  every consumer to fetch `response.scale`/`response.shift` themselves.

There is **no single channel that is both additively reconstructing and a pure
function of one forest** - the multiplier is genuinely per-observation. Stating
that plainly is the honest outcome. Options for resolving it: ship both channels
(rejected on memory: `forestFits` is n x K x draws x chains, and at n = 1e5,
K = 2, 1000 draws that is ~1.6 GB per copy); or ship the raw channel and make
the contribution derivable.

**RECOMMENDATION: ship the raw forest total on the RESPONSE scale**
(`response.scale * f_k`, keeping `$getForestFits`'s meaning up to that scalar
so one name means one thing across the fit object and the sampler), ship `glue`
as the multiplier channel it already is, carry the expanded bases on the fit
(they have no draw axis, so the cost is sum_k n*q_k doubles - negligible), and
document the identity of 5.8.5. `extract(type = "forest")` then returns the raw
total by default and the CONTRIBUTION under `contribution = TRUE`, computed on
demand from raw x basis x glue rather than stored. That gives both readings,
each named for what it is, at one array's memory.

#### 5.8.5 The identity, in its general form

Rev 4 asserted the 2-forest binary-basis instance. The general form, MEASURED
(rev5) at 8.9e-16 on a 3-level factor basis (`glue` 4 x draws, q = 1 + 3):

```
yhat = response.shift + sum_k (B_k %*% g_k) * (response.scale * f_k)
```

with `B_k` forest k's expanded basis (forest 1's implicit intercept is a column
of ones) and `g_k` its slice of the ragged glue. `B_k` is NOT in the run
result - it lives on `data@bases` - which is why 5.8.4 carries the bases onto
the fit: otherwise the documented identity is not evaluable from the fit
object. T-E asserts it at `< 1e-12` on both a binary and a 3-level-factor
basis.

#### 5.8.6 `extract()` and prediction

Add **`type = "forest"`** with a `forest =` selector (index or margin name;
`NULL` = every forest) and `contribution = FALSE` (5.8.4), returning the
packaged slice under the same `combineChains` conventions as the other types
and refusing by name on a fit without forest reporting (5.1's predicate).
**`sample = "test"` is refused by name**: an amplitude-coupled fit cannot have
test fits at all (5.4), so "the same `sample` conventions" would be a promise
the model cannot keep.

`predict()` is NOT extended. Out-of-sample per-forest values need the
engine-side per-forest saved-tree replay door
(`docs/plans/archive/bcf-bartcause-relocation.md`:1271-1273), which this arc does not
open (N8). The door's language is updated at the records commit to record that
the in-sample half landed here.

### 5.9 Underdetermined

**Empty.** The term token was the last open design question and VD decided it
2026-08-16 (5.2). Rev 4's items 2 (the `forests =` formal) and 4 (test bases)
are decided by VD's XOR and by measurement; item 3 (shift attribution) is a
recommendation backed by the 5.8.3 measurement; item 5 (term-vs-named-form
equivalence) is a build-time GATE on S12, not an open question. The one
receipt-driven amendment this revision makes - the desugar target's expression,
5.2.1 - preserves the decided semantics exactly and is recorded in 10.16 rather
than reopened.

---

## 6. Migration map

### 6.1 `dbartsControl` `rngSeed` -> `seed`

| target | sites | class |
|---|---|---|
| dbarts R | `R/dbarts.R`:229, :249, :530-535, :685, :1064; `R/bart.R`:681-683; `R/bart.R`:2386; `R/spec.R`:662; `R/A_class.R`:246, :266, :315-316 | in-arc |
| dbarts C | `src/R_interface_bartcore.cpp`:441, :443 - rebuild only | in-arc |
| dbarts docs | `man/dbartsControl.Rd`:14, :59, :77; `man/dbarts.Rd`:96; `man/samplePriorPredictive.Rd`:70, :73, :80; `man/dbartsSampler-class.Rd`:114 | in-arc |
| design record | `docs/design/public-surface.md`:76, :101 | amend in place |
| dbarts tests | 66 occurrences across 32 files (re-derived by both critiques) | largest mechanical edit |
| bartCause | `R/bcf.R`:205-209 | loud, 1 line |
| stan4bart | forwardable set; `inst/tinytest/test-02-binary.R`:57; after the rename `seed` clears BOTH forwarding loops and is applied twice (benign, same slot); `bart_args = list(rngSeed = 5)` is dropped silently | test edit + 2 recorded consequences |
| CRAN | consumers with `formals()`-derived filters drop their users' `rngSeed` silently | user-facing silent |

### 6.2 xbart

`R/xbart.R`:23-24, :36-47, :49-57, :71-72, :93-95, :138-141, :168-173,
:190-209; `man/xbart.Rd`:86-90. bartCause `R/bayesOpt.R`:30-47 becomes a DROP
filter (`keepTrainingFits`, `keepTrees`, `updateState`, `printEvery`,
`printCutoffs` have no xbart formal and xbart has no dots); :22, :27, :28
survive. Zero CRAN xbart consumers; `iiasa/ibis.iSDM` (non-CRAN) breaks loudly
on `control =`.

### 6.3 Variance quartet -> `varianceForest()`

`R/bart.R`:579-582; `R/dbarts.R`:346-349; `R/spec.R`:612-615, :368-393; the
constructor plus `print` in `R/model.R`; `NAMESPACE`; four Rd files plus a new
`man/varianceForest.Rd`; 12 tinytest files (15 hits). Gate: byte-identical
`attr(control, "bartcore.variance")`. **CRAN row: not "none"** - WeightIt and
MatchIt expose every bart2 formal to their users through `formals()`-derived
accept-sets, so `weightit(..., n.trees.variance = 10)` works today and is
silently dropped afterwards; NEWS.

### 6.4 Family gating, dots, promotions

| target | sites | class |
|---|---|---|
| dbarts | the gating helper plus six sites; the shared unknown-argument helper; two new formals on bart2 and rbart_vi | in-arc |
| **bart()** | **none - by construction** | none |
| any wrapper forwarding one argument set to two families | sees a `dbartsFamilyGatedWarning`, mutable by class. A wrapper that prefers silence pre-filters the gated names or handles the class; that is the wrapper's choice and no part of this design depends on it | user-facing warning |
| bartCause | `R/responseFit.R`:213-216 may drop its `dbartsControl` clause (2 lines). `storage`/`updateState` become reachable through its `formals(bart2)`-keyed treatment filter (`R/utility.R`:165), so `bartc(..., storage = "single")` starts applying reduced-precision storage to the propensity fit too - a numeric change; NEWS both sides | 2 lines + a recorded numeric change |

### 6.5 Prior objects and formula terms

Prior objects are additive. **Formula terms are additive for every censused
consumer**: a term appears only in a formula the USER writes. bartCause builds
its formulas itself (`R/argParse.R`:78, :212) and never constructs a term;
stan4bart, treatSens and bairrtt never call bart2; WeightIt and MatchIt pass
`formula` through from their own callers, so a term written by their user
reaches bart2 unchanged and works. No consumer's accept-set changes, because no
formal is added (N9), and no existing formula's meaning changes, because a
formula without a `forest()` operand is not rewritten at all (S12 block A).

One ecosystem note for the Rd and NEWS, not a migration cost: **stan4bart ships
`y ~ bart(x_1 + x_2)`**, where the token names the predictors the BART component
is built on and `bart` is a pure marker absent from its namespace (survey
section 1, MEASURED there). dbarts' term is a different object - it declares an
ADDITIONAL forest and modulates it - which is why the head is `forest()` and not
`bart()`. The one idiom the two share is the symbolic unnamed slot (5.2.3), and
that sharing is deliberate: a stan4bart user reading `forest(x1 + x2)` gets the
right reading of the slot ("these variables") and needs to learn only what the
head and the colon add.

### 6.6 The bart() shim and fork 7's closes

`bart()` builds the control directly (`R/bart.R`:2374-2387) and `do.call`s
`dbarts()` on a literal list (:2417-2436), sharing three helpers. Fork 7's
closes edit that list (`subset` and `storage` become forwarded formals rather
than hardcoded values; `family` is forwarded), `man/bart.Rd`'s
`\usage`/`\arguments`, and NEWS. Every appended formal's default reproduces
today's behavior, and MEASURED(rev5) the final set breaks no existing `bart()`
abbreviation.

### 6.7 Documentation surface

Split `man/bart2.Rd` out of `man/bart.Rd`; add `bart2` to `_pkgdown.yml`;
rewrite `man/rbart.Rd`:87-89's `k` item; add `man/varianceForest.Rd`; document
the formula terms in `man/bart2.Rd` with a `\seealso` from `man/forest.Rd`.

### 6.8 Partial-matching regression

**Ten** abbreviations stop working (4.1's table): eight on bart2, `n.th` on
xbart, `u` on rbart_vi; `bart()`'s closes break none. No consumer in any census
abbreviates. NEWS-worthy because the failure ("argument N matches multiple
formal arguments") does not name the abbreviation's former target.

---

## 7. Slice plan

Every slice is independently gateable and draw-neutrality-classified.
**The equivalence harness cannot prove neutrality for most of these slices**:
`benchmarks/R/equivalence.R` has zero `xbart` occurrences and reaches `bart2`
only from the four alternate-family scenarios (calls at :942, :971, :1001,
:1039 - anchor corrected from rev 4's scenario-definition lines).

**The in-slice A/B pattern**, used wherever harness coverage is absent: a probe
script fits a fixed matrix of configurations at pinned seeds on the PRE-slice
build, saves draws to a scratch RDS, and re-fits on the post-slice build
asserting `identical()` element by element.

**Recorded harness-coverage note (not this arc's work):** the baseline should
gain a bart2 gaussian scenario, a bart2 probit scenario and an
amplitude-coupled two-forest bart2 scenario; they ride the NEXT scheduled
re-record.

| # | slice | content | draw-neutrality gate |
|---|---|---|---|
| S1 | defect batch | D2, D4, D5, D3, both unqualified `formals()` lookups | harness (four alternate-family scenarios cover D2) + in-slice A/B for D3/D5 |
| S2 | family gating | the helper plus six sites, `samplerOnly` reordering, `dbartsFamilyGatedWarning`, the `prior.scale` residual audit | in-slice A/B; harness on the four scenarios |
| S3 | inert defaults | d1 narrowed to `factors`/`missing`/`proposal.probs`; T-B | in-slice A/B incl. a defaulted-vs-explicit `proposal.probs` pair |
| S4 | dots + promotions | rejection-only dots; `storage`/`updateState` appended; shared helper; NEWS | in-slice A/B |
| S5 | naming | `rngSeed` -> `seed`; `xbart(sigma =)` -> `sigest`; T-A | in-slice A/B plus the full harness; `--preclean` |
| S6 | variance collapse | `varianceForest()` + `print`, three formals removed on all three surfaces, docs, 12 test files | byte-identical control attribute; in-slice A/B for the bart2 path |
| S7 | prior objects | three formals appended; four collision refusals; the `resid.prior` family diagnosis | in-slice A/B |
| S8 | xbart | `control =` removed, four flat knobs, `tree.prior`, shape checks | `inst/tinytest/test-xbart-*.R` |
| S9 | split.probs default | bart2/rbart_vi -> `NULL` | in-slice A/B; identity proved by construction |
| **S10** | **bart() parity** | append `subset`, `storage`, `family = c("auto", "logistic", "aft")`; forward all three in the literal list; **by-name refusals pointing at bart2 for the four own-class tokens** (`match.arg`'s generic "'arg' should be one of" names neither bart2 nor why, and the existing ordinal treatment at :2362-2372/:2440-2445 is ordinal-specific while multinomial's refusal lives in `dbarts()`, not `bart()` - MEASURED(r3)); 2 Rd items; the rule paragraph; tests | in-slice A/B on `bart()` (harness covers it as the default route) **plus a regression-shaped pair: a BINARY `bart()` fit at HEAD and at the slice must be bitwise identical** - the F1 failure mode |
| **S11** | **fork 6a: per-forest output channel** | thread `samples$forestFits`/`samples$glue` through `packageBartResults` via `shapeMultinomialChannel` (trailing margin); `forest1..K` names plus a `forest.labels` attribute; carry the expanded bases; response-scale conversion; `extract(type = "forest", forest =, contribution =)`; refusal on fits without forest reporting; T-E | existing single-forest fits bitwise unchanged (in-slice A/B); T-E's identity `< 1e-12` on binary AND 3-level-factor bases |
| **S12** | **fork 6b: formula-term ingestion** | the both-sides AST walk with ancestor-chain context tracking and non-forest-operand recursion (5.3 steps 1-2), the colon desugar (5.2.1) and the left-operand grammar (5.2.1a), the symbolic-vars grammar and its subset-of-the-design rule (5.2.3, 5.3 step 4), AST-surgery rewrite, post-subset basis evaluation, the twelve refusals of 5.6 (including `*`, the ancestor-chain check, both grammars, the factor-member and colon-chain cases, and the family matrix at family-resolution time), docs | **draw-neutral for every non-term call**, gated by the A/B matrix below; term-using calls are NEW paths gated by tests, including the desugar-equivalence and refusal cells below |
| **S13** | **fork 6c: formula-path bases subsetting** | converge `dbarts(forest(basis =))` onto the term route's post-subset rule | **BEHAVIOR CHANGE, not a defect fix** - see below |
| S14 | docs | `man/bart2.Rd` split, `_pkgdown.yml`, the defaults statement, T-C, term documentation | none (no code) |

**S12's test matrix**, in three blocks.

*(A) Draw-neutrality A/B - non-term formulas, run with and without the slice,
asserting `identical()`:* (1) `y ~ a + b`; (2) `y ~ .`; (3) `y ~ . - z`;
(4) `y ~ a + z:b` (a `:` with NO forest operand - the walk must ignore it);
(5) `y ~ a - 1`; (6) `y ~ a + offset(o)`; (7) a back-tick-requiring name;
(8) `y ~ a` with `subset`; (9) `y ~ a` with `weights` and `offset`; (10) the
x/y matrix interface, which the walk never touches. Cells 2, 4, 5 and 6 are the
ones a label-based recipe would have broken; the AST walk handles 2 and 6 by
construction (MEASURED(rev6)) and 4 by the forest-operand predicate.

*(B) Desugar equivalence - the sugar must be BYTE-IDENTICAL in sampler
configuration to the hand-written named form at the same seed:*
(11) `z:forest(x1 + x2)` vs `forest(x1 + x2, basis = ~ z)`, numeric z;
(12) `zf:forest(x1 + x2)` vs `forest(x1 + x2, basis = ~ zf)`, factor zf -
asserting the basis is n x nlevels with no reference dropped;
(13) `factor(z):forest(...)` vs `forest(..., basis = ~ factor(z))` - n x 2;
(14) `(a + b):forest(...)` vs `forest(..., basis = ~ cbind(a, b))` - ONE forest,
n x 2, NOT the elementwise sum;
(15) `forest(x1 + x2):z` vs `z:forest(x1 + x2)` - operand order immaterial;
(16) `dbarts::forest(...)` qualified head reaches the same configuration;
(17) the symbolic slot vs `vars = ` by name;
(18) a symbolic factor name expanding to all its indicator columns under
`factors = "indicators"`.
Assertion shape, with every comparand named as concretely as S6's:
`identical(data@bases)`, **`identical(attr(control, "bartcore.bcf"))`** - the
resolved per-forest spec object, a list of `params` (one length-8 numeric per
forest), `vars`, `interactions` and `blocks`, set at `R/spec.R`:559 and the
same attribute Resolution 3 / S6 already names as fork 2b's byte-identity gate
(rev 6 wrote "resolved `forests` specs", which names no object: MEASURED(v4),
`control` has no forests slot) - and identical draws at a pinned seed.
Byte-identity is available here (unlike T-E's float sum) because the sugar
rewrites to the named call before anything numeric happens.

*(C) Refusal cells, one per 5.6 item except item 7 (enforced by construction),
plus the two 5.5.2 collisions:* (19) `z*forest(...)` [5a]; (20)
`I(z:forest(x1))` and `(x1 + z:forest(X))^2` [5, ancestor chain]; (21) a term
on the LHS, `forest(x1) ~ x2` [5]; (21a) `y ~ x1 - z:forest(X)` [5, `-`];
(22) `forest(log(x1))`; (23) `forest(x1 - x2)`; (24) `forest(.)` [22-24: 5b];
(24a) a symbolic name absent from the RHS, `y ~ x2 + z:forest(x1 + x2)` [5c];
(24b) `(a + zf):forest(X)` factor member AND `(a + zc):forest(X)` character
member; (24c) `z:w:forest(X)`, colon chain;
(24d) `forest(A):forest(B)` and `z:forest(A):forest(B)` [24b-24d: 5d, 5e];
(25) both `:` sugar and `basis = ` [5.5.2]; (26) both the symbolic slot and
`vars = ` [5.5.2]; (27) `y ~ z:forest(x1)` alone AND `y ~ 0 + z:forest(x1)`
[6, both residues]; (27a) an all-zero basis column [1 - new work, since
`expandForestBasis(rep(0, n))` is accepted today]; (28) a term with `test` [4];
(29) a term with a pre-built `dbartsData` [3]; (30) each family item 2's
refusal is specced to cover - `hazard`, `hazard.probit`, `hazard.logistic`,
`multinomial`, `hurdle.lognormal`, `twopart` - plus the three already refused
by the engine (`aft`, `ordinal`, `nbinom`), each naming the family and the
term [2].

**S13's reclassification.** MEASURED(rev5), BOTH of these are ACCEPTED today,
because `validateForestBases` on the formula branch checks only the COUNT
(`R/data.R`:895):

```
dbarts(y ~ a, d, subset = 1:20,  forests = list(forest(), forest(basis = <20-row>)))  -> n = 20
dbarts(y ~ a, d, subset = 21:40, forests = list(forest(), forest(basis = <20-row>)))  -> n = 20
```

The second is the ambiguous case the recorded door names ("it makes a shipped
argument's contract ambiguous at equal row counts",
`bcf-bartcause-relocation.md`:1283-1284): a 20-row basis supplied against a
20-row subset of a 40-row frame is silently row-aligned by position. Under the
post-subset rule it either refuses or selects different rows. **The decision:
with `subset` present the basis must be FULL-DATA length and is subset by the
same index; an equal-count-but-not-full-length basis is refused by name.** That
changes a working call, so S13 carries a NEWS entry and an A/B over
{subset present} x {basis rows == subset length} x {basis rows == data rows},
not "tests only".

Sequencing: S1 independent. S2 and S3 independent. S5 before S4's T-A
assertion. S6, S7, S9 mutually independent. **S11 before S12** (the channel
exists before a new way to create the fits that fill it) and **S12 before S13**
(S13 converges the older route onto S12's rule). S14 last. No slice re-records
a baseline.

Count: **14 slices**, of which 3 are fork 6's and 1 is fork 7's.

---

## 8. Resolutions (VD walkthrough, 2026-08-16)

1. **Fork 1 = A1, flat-plus-collapse.** No `control =` on bart2; A2' recorded
   as available-later, with stale `bartcore.*` attribute inheritance as the
   true against-ground. (3.a; slices S6, S7.)
2. **Fork 2 = Complete, 54 formals.** The three prior-object formals land,
   making `linear()`, `gp()` and `fixed()` reachable from bart2. (4.1, 4.2;
   slice S7.) *Rev 4 recorded this as "55 after fork 6's `forests`"; that was a
   drift and is reverted - see 10.11.*
3. **Fork 2b = `varianceForest(vars, n.trees, base, power)`**, the plain
   selector retained as the same argument's other accepted type. (3.b b3, 4.2;
   slice S6.)
4. **Fork 3 = L1 warn, grounded in R standards** - the stats-package
   disregard precedent, the caller-owned severity policy, and the classed
   condition `dbartsFamilyGatedWarning` (whose base-R analogue is
   `chkDotsWarning`). Monotonicity stands. (3.c.5; slice S2.)
5. **Fork 4 = route (i)**, a standalone helper at the six family-resolution
   sites, keyed on the `matchedCall` snapshot at `R/bart.R`:655. Route (ii)
   REFUTED and recorded. (3.c.4; slice S2.)
6. **Fork 5 = rejection-only dots**, with a stated sunset to
   deletion-equivalent. (3.e; slices S4, S5.)
7. **Fork 6 = Level 1 BUILD the in-sample per-forest output channel in-arc;
   Level 2 = FORMULA TERMS**, an XOR against a flat `forests =` formal.
   Out-of-sample per-forest prediction stays doored. (Section 5; slices S11,
   S12, S13.)
   **TERM TOKEN DECIDED (VD 2026-08-16, after the conventions survey and two
   probe rounds): head = `forest()`, with COLON SUGAR as the canonical
   spelling** - `y ~ x1 + x2 + z:forest(x1 + x2)` and
   `factor(z):forest(x1 + x2)`. `:` desugars to a no-intercept basis on that
   forest; a compound left operand folds into ONE forest with a multi-column
   basis; `*` with a forest operand is refused naming both explicit forms;
   `forest(basis = ~z)` remains the general named form and the desugar target,
   so there is one semantic channel with sugar over it (the fork-2b pattern);
   extraction is an AST walk, not `terms(specials =)`, matching bare and
   `dbarts::`-qualified heads. In TERM CONTEXT the unnamed slot is a SYMBOLIC
   predictor set (`forest(x1 + x2)` means vars {x1, x2}, a la stan4bart's
   `bart(x1 + x2)`), interpreted from the AST and never evaluated - a
   deliberate divergence from the constructor, whose first positional formal is
   `basis`, which the docs neutralize by always spelling `basis =` by name.
   (5.2, 5.2.1-5.2.3, 5.3.)
8. **Fork 7 = CLOSE `subset`, `storage`, and `family` restricted to the
   families that both package as ordinary `"bart"` objects AND are not already
   reachable from the response** - i.e. `"logistic"` and `"aft"`, with
   `"auto"` as the first token. Revised under critique r3 from rev 4's list:
   `"gaussian"`/`"probit"`/`"hazard.probit"` are dropped as tokens because
   `"auto"` plus the response already reaches them, and **`"hazard"` and
   `"hazard.logistic"` are dropped from the close** because their required
   companions `breaks`/`max.rows` are not `bart()` formals (MEASURED(rev5)) and
   the person-period expansion changes the row margin of a `yhat.train` that
   `bart()`'s Rd documents as n-column. Closing hazard would need five appended
   formals and a return-contract amendment; recorded as a door instead. The
   four own-class families are refused BY NAME pointing at bart2 (S10), not by
   `match.arg`. LEAVE `offset.test`, `factors`, `missing`, `updateState`, the
   own-class families, and the bartcore-era block. Position rule: APPEND AT
   END. (4.5, 8.7; slice S10.)
9. **Fork 8 = parity in-arc as S10.**

### 8.7 The fork-7 disposition table as decided

| bart2 capability | disposition | note |
|---|---|---|
| `subset` | **CLOSE** | append; `subset = NULL` reproduces :2421 |
| `storage` | **CLOSE** | append; reduced-precision STORAGE (`src/R_interface_bartcore.cpp`:394) - it changes the numbers a fit returns, which the Rd must say |
| `family = "logistic"` | **CLOSE** | packages as `"bart"` (:348); not reachable from the response |
| `family = "aft"` | **CLOSE** | packages as `"bart"`; the `Surv`/two-column response already ingests |
| `family = "gaussian"`, `"probit"`, `"hazard.probit"` | LEAVE as tokens | reachable today via `"auto"` + the response; adding tokens would be spelling, not capability |
| `family = "hazard"`, `"hazard.logistic"` | **LEAVE (revised)** | needs `breaks`/`max.rows`, which `bart()` lacks, and breaks the documented row margin |
| `family = "multinomial"`/`"ordinal"`/`"nbinom"`/`"hurdle.lognormal"` | LEAVE | own classes; refused BY NAME pointing at bart2 |
| `offset.test`, `factors`, `missing`, `updateState` | LEAVE | recorded |
| `dart`, `monotone`, `interactions`, `blocks`, variance forest, `warm.start`, `n.grow.sweeps`, `dispersion`, prior objects, multi-forest | LEAVE | bartcore-era / modern surface |

Rejected for the whole table: a `...` on `bart()` forwarding into `bart2`.

---

## 9. Census defect dispositions and open questions

### 9.1 Dispositions

**D1.** `sigest` unforwarded on ordinal/nbinom; diagnosed by route (i) from the
snapshot, no new forwarding. (S2.)
**D7.** `sigest`/`sigma` inert on EVERY fixed-unit-scale family; LIVE on aft and
gaussian. (S2.)
**D2.** `keepSampler` ignored on every alternate-family path; retain `$bc`
under `control@keepTrees` and `$fit` under `keepSampler`. (S1.)
**D3.** The burn-in branch tests the RAW `n.burn`; no observable consequence;
code-clarity fix. (S1.)
**D4.** `monotone` with an explicit `proposal.probs = NULL` errors where the
defaulted call succeeds; d1's alignment plus a one-line guard at
`R/spec.R`:278. (S1/S3.)
**D5.** rbart_vi's `k = 2.0` is INERT and `man/rbart.Rd`:87-89 documents the
opposite. (S1, draw-neutral.)
**D6.** Forwarded formal defaults are inert generally; only `proposal.probs`
differs textually. (S3 plus T-B.)
**N-new-1.** Formula + `subset` + a declared basis: MEASURED(tip),
`dbarts(y ~ a, d, subset = , forests = list(forest(), forest(basis = ~ tr)))`
errors "length of 'basis' must equal length of 'y'"; the same call without
`subset` works; `dbartsData(x, y, subset =, bases =)` works on the x/y branch
(`R/data.R`:977, :1043) and errors on the formula branch. S13 closes it, as a
behavior change (section 7).

### 9.2 Open questions

**Q1.** CLOSED by `bayestree-verification.md`.
**Q2.** `keepTrees`/`keepSampler` are two flags for one idea; a door.
**Q3.** `combineChains` on a fitting function: leave.
**Q4.** Is `n.trees = 75` the modern recommendation? A live statistical
question; no slice depends on it.
**Q5.** xbart/rbart_vi model-surface gaps: capability tickets.
**Q6.** Consumers with `formals()`-derived filters drop their users' unknown
options with no diagnostic; the lever is NEWS.
**Q7.** stan4bart's `bart_args` has no name validation - a stan4bart ticket.
**Q8.** A reused control's stale `bartcore.*` attribute is a family-gated
inertness rule L cannot see; bears only on A2'.
**Q9.** dbarts resolves the family from response VALUES at runtime, so a
source-level audit of a consumer's family is never sufficient.
**Q10.** The term token (5.9) - the one open design question.
**Q11 (new).** Whether `bart()` should ever gain `breaks`/`max.rows` and the
hazard families (8.7's revised row) - a door, not a gap.

---

## 10. Adjudication appendix

**10.1 r2 B1/B2/B3 (route ii) - ACCEPTED IN FULL**, re-verified at tip.
**10.2 r2 B4 (L4) - ACCEPTED for route (ii); NARROWED under route (i)**, where
the discriminator is recoverable from the snapshot. Not adopted: two severities
for one mistake, and the error branch protects nothing the warning does not.
**10.3 r2 M6 - ACCEPTED**, and it produced the monotonicity rule.
**10.4 r2 M3 (equivalence coverage) - ACCEPTED**; it changed the slice plan.
**10.5 r2 M4 (partial matching) - ACCEPTED**, re-run this revision for the
final formal sets (ten collisions, not twelve, once `forests` is gone).
**10.6 r2 M7 (variance duplication) - ACCEPTED**, adjudicated in G1.
**10.7 r2 M8 (T-A as a ratchet) - ACCEPTED**; T-A is scoped to 1.0-0's slots.
**10.8 r2 M2 - ACCEPTED**, and rev 4 corrected Level 1's pricing in the
opposite direction: the channels already exist.
**10.9 r1 9.2 (b3 selectors) - NARROWING RETAINED**, reproduced by r2.
**10.10 r2 m1-m8 - ALL ACCEPTED** and carried.

**10.11 r3 F3 (COMPLIANCE, `forests =`) - ACCEPTED IN FULL; rev 4's propagation
reverted.** VD decided fork 6 as an XOR and rev 4 implemented the branch he did
not choose across eight sites, including Resolution 2, whose stated purpose is
that nothing be re-litigated. Reverted at every site: the signature (4.1), the
arithmetic (54, not 55), the abbreviation table (eight bart2 collisions, not
ten - `fo`/`for`, which rev 4 itself called "the sharpest", existed only
because of it), T-B, the migration map (6.5), 6.8's count, S12's content, and
Resolution 2. What survives is one door sentence in 5.5.3, and rev 5 additionally
supplies the argument FOR revisiting it that rev 4 did not have: forest 1's `sd`
and `update.amplitude` have no term spelling (5.7). I did not find the
formula-terms route unimplementable without a `forests =` formal - it is
implementable, at the cost recorded in 5.7 and G3.

**10.12 r3 F19 (COMPLIANCE, consumer residue) - ACCEPTED.** The error-always
rejection's first clause ("at any wrapper that forwards one argument set to two
models") was a consumer-behavior argument in a section VD ruled must not use
one. Rewritten on receipts (1) and (2) alone; the wrapper consequence stays in
6.4, where it is explicitly disclaimed as non-load-bearing.

**10.13 r3 F1 (BLOCKER, `family` default) - ACCEPTED.** MEASURED(rev5)
reproduces the critic's receipts exactly. The first token is `"auto"`, the
token set is narrowed to what `"auto"` cannot reach, and S10 gains the
binary-`bart()`-bitwise regression pair. The criterion is now stated out loud
in 8.7 rather than implied.

**10.14 r3 F2 (BLOCKER, shift attribution) - ACCEPTED, and the prior question
is now measured.** `forestFits[, k]` holds the RAW forest total (bitwise equal
to `$getForestFits(k)`), not the glued contribution (5.8.3). Rev 4's
recommendation (a) is refuted on its own terms. The honest finding - that no
single channel is both additively reconstructing and a pure function of one
forest - is stated, and the recommendation is raw + glue + bases with the
contribution derived on demand in `extract()` (5.8.4).

**10.16 The desugar target's expression - AMENDED WITH RECEIPTS, semantics
preserved.** The decision file writes `z:forest(X) -> forest(X, basis = ~0+z)`
and `(a+b):forest(X) -> basis = ~0+a+b`, which is model.matrix idiom. The
shipped channel is not a model.matrix channel: `evaluateForestBasis`
(`R/model.R`:751-763) evaluates the RHS as an EXPRESSION and `expandForestBasis`
(:772-...) expands the RESULT. MEASURED(rev6) through the shipped path:
`~0+z` -> 60x1 (right, by arithmetic accident), `~0+zf` -> **ERROR "a 'basis'
cannot be NA"** with an `Ops.factor` warning, `~0+a+b` -> 60x1, the ELEMENTWISE
SUM. Taken literally the target would break `factor(z):forest(X)` - one of the
two canonical spellings the decision itself lists - and silently give one column
where two were meant. The targets are therefore written `~ z`, `~ zf`,
`~ factor(z)`, `~ cbind(a, b)` (MEASURED 60x1 / 60x3 / 60x2 / 60x2), which
deliver the decided semantics exactly: no intercept column for a numeric
operand, all K indicators with no reference dropped for a factor, and one forest
with a multi-column basis for a compound operand. Recorded here rather than
treated as a reopened fork, because nothing about the decision's meaning
changes - only the spelling of the expression the sugar rewrites to.

**10.15 r3 F4-F25 - ALL ACCEPTED** and folded in: trailing forest margin via
`shapeMultinomialChannel` (F4); `< 1e-12`, never "bitwise", for a float sum
(F5); `terms(data =)` (F6); the `term.labels` scan plus the `specials`-index
referent, and the term-in-expression refusal (F7); the family x term matrix
with hazard named explicitly at family-resolution time (F8); `test` decided by
measurement as a by-name refusal (F9); S13 reclassified as behavior-changing
with its A/B cells and a stated resolution of the equal-count ambiguity (F10);
5.7 corrected (F11); the general identity with `B_k`, and the bases carried on
the fit so it is evaluable (F12); by-name own-class refusals in S10 (F13);
hazard dropped from the close (F14); S12's A/B matrix enumerated and the
attribute-preserving rebuild specified (F15); the namespace-qualified-escape
note (F16); `warningCondition`'s class argument (F17); `chkDotsWarning` as a
second precedent (F18); 5.6 item 6's vacuous example replaced (F20);
`forest1..K` everywhere plus `forest.labels` (F21); `sample = "test"` refused
(F22); N-new-1's two textual slips (F23); the amplitude-coupled predicate
(F24). F25 found nothing to fix.

---

## Landing notes (append-only)

S1 (defect batch) LANDED ebc57af7, 2026-08-16. D2: keepSampler retains
$fit on multinomial/multinomial-counts/ordinal/nbinom independent of
keepTrees ($bc and raw channels still keyTrees-gated); hurdle verified
already-forwarding, untouched. D3: runWithBurnIn keys on
control@n.burn, shared fix reaching bart(). D4: !is.null guard at the
monotone proposal-probs comparison. D5: rbart_vi k default 2.0 -> NULL
+ man/rbart.Rd rewritten. Both unqualified formals() lookups
qualified. Gates: implementer battery + independent gate-runner, both
green (5337/0 tinytest; trio 37+12+10 bitwise; R CMD check OK; air +
lintr clean; NEWS 246-entry parse). Mutation-proven: D2 reverts fail
exactly the three new assertions with hurdle staying green; D4 revert
reproduces the error. Draw-neutrality: equivalence trio bitwise plus
four A/B arms (bart2 and bart() burn corners, binary rbart_vi, the
two proposal.probs spellings) all identical() - the bart()-side
corner arm was the gate-runner's own addition.

S2 (family gating) LANDED 2eb6a1a4, 2026-08-16. warnFamilyGatedArgs +
a data-driven inventory (R/utility.R) wired at the six 3.c.4 sites;
ONE classed warning per call (dbartsFamilyGatedWarning, a
dbartsWarning), suppressWarnings(classes=)-silenceable; bart() silent
by construction. Standard site and rbart_vi resolve "auto" locally
(sampler$model@family strips hazard's remap;
resolveClassificationFamily leaves numeric responses unresolved).
Orchestrator review caught a hurdle leak: redirectCall forwarded
gated names into both recursive component calls (false probit sigest
warning; dispersion warned 3x) - fixed by stripping gated names from
the component calls, trio kept on the live positive half; warning
counts 0/1 mutation-proven (revert -> 1/3). prior.scale residual
audit: MEASURED LIVE under probit/logistic/ordinal/nbinom - stays out
of the inventory. T-D lands in inst/tinytest/test-argument-surface.R
(the arc's maintained-contracts home; T-A/B/C/E ride their slices).
Gates: implementer + independent runner both green (5353/0; trio
37+12+10 bitwise; 11/11 behavior probes; mutation 8-of-16 fail
exactly; A/B four configs identical; check OK with zero gated-warning
hits in the log).

S3 (inert defaults) LANDED 6f1ba79e, 2026-08-16. factors/missing
match.arg-resolved in bart2's/rbart_vi's own frames and stamped onto
matchedCall unconditionally AFTER the family-gating argNames snapshot
(the S2 suppliedness contract holds - probe-verified both directions);
proposal.probs default NULL -> dbarts()'s named 4-vector verbatim,
always forwarded resolved. T-B lands in test-argument-surface.R with
the plan's 6-row exceptions table (three rows inert until S7);
mutation-proven (a reordered vector fails T-B naming proposal.probs).
man/bart.Rd usage updated. Gates: implementer + independent runner
green (5363/0; trio 37+12+10 bitwise; 9-config A/B incl. hurdle
components identical; check OK twice after one transient
environment flake disproven by reruns).

S4 (dots + promotions) LANDED ec59caf6, 2026-08-16. storage/
updateState appended as formals on bart2 (54 incl. dots) and rbart_vi
(44), flowing through the existing redirect/backfill machinery; both
hand-written union checks replaced by shared rejectUnknownDotsArgs
(R/utility.R): every dots name errors with a nearest-formal agrep
suggestion; retired-name table starts empty (the fork-5 sunset
mechanism); rngSeed keeps flowing via a passthrough S5 deletes.
setdiff(control formals, bart2 formals) == "rngSeed" exactly - T-A
ready, lands with S5. Gate anomaly resolved during the run: the
implementer's equivalence deviation traced to comparing against the
STALE equivalence-a825263.rds without --strict-coverage (pre
empty-leaf-veto re-record; zeroweights moved legitimately at
21fc29c3) - the canonical 8b047f8b run is 37/37 bitwise, 0 skipped,
on the same build. Gates: both runners green (5377/0; trio bitwise;
mutation exactly 3 rejection tests; A/B incl. dots-vs-formal
updateState identical; check OK). Note: storag= now partial-matches
the storage formal by R's own semantics - inherent, recorded.

S5 (naming) LANDED ba4f761d, 2026-08-16. rngSeed -> seed everywhere
the 6.1 map lists plus live-call sites the sweep found beyond it
(benchmarks/R harness scripts, vignette, workflow prose,
docs/architecture.md); the two C string literals flipped, the
internal ParsedControl::rngSeed field untouched (N2's scope); bart2's
hand-rename block deleted; S4's passthrough retired into the first
retiredDotsNames row ("'rngSeed' was renamed to 'seed'"). xbart sigma
-> sigest with man/xbart.Rd restating the b2 rule; the shared
validateArgumentsInEnvironment gained a parallel sigest formal so
dbarts() keeps its own sigma validation (an implementation finding
not in the map, recorded). T-A LANDS: control formals are a subset of
bart2's and rbart_vi's, spelled identically - mutation-proven both
runs. Gates: both runners green (5381/0; trio 37+12+10 bitwise
canonical files; sweep verified per-area - 57 rngSeed remains, all
policy areas: plans/NEWS/retired-row/struct-field/negative-test plus
the two design-record in-place amendments; 4-config A/B identical
across spellings; check OK). Lockstep pairs LANDED same
day: bartCause a3aa98d (bcf.R + two test dbartsControl calls; suite
green, 1 recorded environmental stan4bart-ABI error untouched by the
diff), stan4bart 54e157b (test-02 spelling + a NEWS breaking-changes
bullet on the rename and the bart_args silent-drop consequence; full
455-result suite green against a rebuilt S5 dbarts - a whole-ABI
validation in passing).

S6 (variance collapse) LANDED 5aa4d292, 2026-08-16. varianceForest(
vars, n.trees, base, power) exported beside forest() with print/
format; the three .variance formals removed from bart2/dbarts/
dbartsSpec; the resolver's object branch routes vars through the
SHARED resolveVarianceColumns (object vars=NULL -> every column,
distinct from variance=NULL's no-forest - the two NULLs are
opposites, translated before the shared resolver); knobs land where
the flat formals' values landed. 12 test files + 2 benchmark scripts
migrated - the harness's own het scenarios now spell varianceForest
and still compare BITWISE against the untouched baseline, the
strongest round-trip proof available. New Rd + _pkgdown.yml +
pkgdown check; C error strings reworded off the retired names.
Gates: both runners green (5400/0; trio 37+12+10 bitwise; byte-
identity grid incl. full-knob attr verbatim; mutation exactly 2;
check OK).

S7 (prior objects) LANDED de2212bc, 2026-08-16. tree.prior/node.prior/
resid.prior = NULL on bart2 in 4.1's literal ordering (after
max.rows); forwarded UNEVALUATED through buildSamplerPriors'
matchedCall reads (the k/nodeK idiom - bare linear()/gp()/fixed()
resolve in the caller's frame); presence-based collision refusals at
four sites (resid.prior's set includes sigest per m6); explicit NULL
behaves as absent, documented in place. resid.prior joined the
family-gating inventory as a data row (liveIn gaussian/aft/
hurdle.lognormal); hurdle strips it occupancy-side, tree/node priors
flow to BOTH components (liveness-verified). rbart_vi deliberately
unchanged (rbart-model-surface-parity ticket). G3 CLOSED for the
prior family: linear()/gp()/fixed() reachable from bart2, verified
draw-identical against the dbarts() route. Two more abbreviations
break (t=, resid.=), asserted + NEWS. Gates: both runners green
(5415/0; trio 37+12+10 bitwise; mutation isolates exactly the 4
collision tests; 5-config A/B identical incl. warn-count parity).

S8 (xbart) LANDED d54c2791, 2026-08-16. control= removed; n.cuts/
useQuantiles/n.thin/storage appended spelled and defaulted exactly as
dbartsControl's (four, not six, per f1) plus tree.prior = NULL; 32
formals, no dots, control= now a native unused-argument error. xbart
builds its control internally - the knobs land where the flat values
landed (the S6 precedent), the six forced fields set at construction,
the validate call carrying the constructed object explicitly.
tree.prior follows the grid-axis-overrides-the-object rule node.prior/
k established: power/base/k stay legal alongside it (grids override
the object per cell, proven byte-identical against an object naming
different power/base), while dart/split.probs collide by name - a
deliberate divergence from bart2's collision set, where power/base
are ordinary scalars (recorded in code comment + NEWS).
Implementation findings: the unnamed-n.trees-from-control fallback
deleted as unreachable once control= is gone (A/B-covered); storage =
"single" is unconditionally refused on xbart's shared per-fold
data-handle path regardless of family (pre-existing engine fact,
error identical on both A/B arms; man/xbart.Rd states this plainly
rather than echoing dbartsControl's gaussian-only wording). n.th=
abbreviation breaks (6.8), asserted + NEWS. Gates: both runners green
(5440/0; trio 37+12+10 bitwise canonical files; 9-config implementer
A/B + 7-probe gate-runner A/B all identical incl. old-vs-new
spellings per arm; mutation isolates exactly the n.thin plumbing,
n.cuts shape-check, and grid-override assertions; air/lintr clean;
NEWS 254-entry parse; check OK twice).

S9 (split.probs default) LANDED 8b263ae5, 2026-08-16. bart2's and
rbart_vi's split.probs default respelled 1 / num.vars -> NULL (d4:
identical by construction). T-B disposition: dbarts() has no
split.probs formal of its own (tree.prior = cgm carries it), so T-B
never compared the name - no exception row exists or is needed; slice
coverage asserts both NULL formals. Usage blocks updated; the shared
Rd prose already documented NULL. Gates: both runners green (5417/0
own-base; trio 37+12+10 bitwise; 4/5-cell A/B bitwise with the
process-unique sampler pointer excluded; mutation isolates exactly
the bart2 assertion, rbart_vi's stays green; check OK twice). Landed
rebased onto the S8 chain (own-base gated sha c96af760).

S10 (bart() parity) LANDED 726dab10, 2026-08-16. subset = NULL,
storage = c("double", "single"), family = c("auto", "logistic",
"aft") appended - 36 formals - and forwarded through the literal
dbarts() list; family's first token "auto" (F1), the binary bart()
regression pair bitwise in both runners. Implementation finding
beyond the migration map: the literal list forwarded subset = NULL
EXPLICITLY, and dbarts()'s survival path refuses on bare
!missing(subset) (R/dbarts.R:472/:530), which would have foreclosed
the family = "aft" close - subset is now omitted from the list unless
non-NULL so missing() propagates; aft completes, aft + subset still
refuses cleanly. The four own-class families refused BY NAME pointing
at bart2 through one shared helper; both ordinal-response sites
(matrix pre-check, formula post-check) route through it, so the
refusal vocabulary is uniform. T-C lands: 8.7's disposition table as
a data.frame with totality, CLOSE-formal-existence, and token-set
assertions. No abbreviation breaks (asserted live). Gates: both
runners green (5437/0 own-base; trio 37+12+10 bitwise; F1 pair +
subset==manual-presubset bitwise; mutation isolates family-order/
ordinal-refusal/subset-forwarding exactly; check OK twice). Landed
rebased (own-base gated sha 6f0eb9c5).

S11 (fork 6a, per-forest output channel) LANDED 00abf336,
2026-08-16. packageBartResults packages forestFits - the
RESPONSE-scale raw totals response.scale * f_k, no shift, no glue
folded in - through shapeMultinomialChannel with the forest margin
LAST and forest1..K names, and glue as the ragged multiplier channel
with a forest-major `forest` attribute; the expanded bases ride the
fit so the 5.8.5 identity is evaluable from the fit alone;
declaration labels (names(forests)) ride as a forest.labels
attribute, display-only, never a second selector vocabulary.
extract(type = "forest", forest =, contribution =) serves the raw
slice or the on-demand contribution (basis %*% glue) * raw under
both combineChains conventions via the new reshapeChainedChannel
(round-trip verified against shapeMultinomialChannel's own forms);
refusals by name on fits without forest reporting and on sample =
"test"; predict() untouched (N8). Orchestrator review caught the
packaging gating on the forest COUNT where 5.1's predicate is the
COUPLING (a K-forest multinomial run carries neither channel) -
amended to numForests > 1 && !is.null(samples$forestFits), the
count/coupling split proven by a defect-reproduction probe on a
nulled-channel two-forest sampler. T-E LANDS: the reconstruction
identity < 1e-12 on binary AND 3-level-factor bases (8.9e-16 both,
independently re-derived by the gate-runner from the identity alone).
Gates: both runners green (5440/0 own-base; trio 37+12+10 bitwise;
single-forest, multinomial, bart(), rbart_vi fits bitwise unchanged;
two-forest pre-existing channels bitwise; mutation isolates
scale/margin/bases exactly; check OK twice).

BATCH NOTE. S9+S10+S11 landed as a batch under the shared
merged-tree battery clause, each slice individually gated on its own
base 3a5d1dc0 by implementer plus independent runner first. Merged
tree at 00abf336: 5489/0 tinytest; trio 37+12+10 bitwise; cross-slice
A/B vs the landed S8 tip 0adbf59a (binary/continuous bart(), bart2,
rbart_vi, xbart, two-forest channels, subset==manual, T-E 4.44e-16)
clean; NEWS 257-entry parse; check OK. S8's OWN CI was six-green
before the batch pushed. Pre-existing anomaly recorded, NOT
batch-caused: an S6-era print/format assertion sits displaced at
test-argument-surface.R:788 (present in 0adbf59a verbatim) - fold
into the arc-closure comment sweep. bartCause lockstep for S8 landed
same day: 31519d22 on dbarts-1.0 (bayesOpt.R's move-into-control
block replaced by a formals()-driven drop filter; probe pair proves
old-errors/new-completes; suite 724 pass with only the recorded
environmental stan4bart-ABI failure, proven pre-existing by an
identical-baseline check).

S12 (fork 6b, formula-term ingestion) LANDED 4b179585, 2026-08-17.
forest() is now a formula term on the shared dbarts()/bart2 formula
path: R/formulaTerms.R carries the both-sides AST walk with
whole-ancestor-chain legality (only +/:/*/~ above a hit; LHS,
removals, I(), ^, any other call refuse by name), hit detection for
bare additive terms and :/* nodes with a forest operand including the
dbarts::-qualified head, recursion into non-forest operands (both-
sides heads found), AST surgery preserving intercept/offset()/
back-ticks/environment by construction, the 5.2.1a closed left-
operand grammar with the desugar targets ~ z / ~ factor(z) /
~ cbind(a, b) through the shipped evaluateForestBasis channel, the
symbolic unnamed slot (never evaluated, subset-of-the-rewritten-RHS,
factor names expanding to all indicator columns via
resolveTermColumns), post-subset basis evaluation via a dedicated
model frame, the family matrix refused at resolution time (hazard*,
multinomial, hurdle.lognormal via bart2-side guards, aft/ordinal/
nbinom uniformly reworded), test= and pre-built-data refusals, both
5.5.2 collisions, the y ~ 1 / y ~ 0 residues, and the all-zero
basis-column refusal (in expandForestBasis, so the forests= route
gains it too - a shared-enforcement finding). TWO REVIEW-DRIVEN
ADDITIONS: a term-bearing formula given TOGETHER with forests= is a
by-name collision refusal (the plan left the combination
unspecified; recorded here as the resolution), and nested colon
chains in BOTH associativity directions (forest(x1):z:a,
a:(b:forest(x1))) refuse by name rather than falling to
model.frame's generic error - the walk originally skipped legal-
context :/* nodes without a direct forest operand. bart2 gains NO
forests= formal (the decided XOR). Docs: a term-context paragraph in
man/forest.Rd; the full bart2 term documentation rides the docs
slice. Gates: implementer + independent runner green (5569/0 full
suite, 80/80 term file; block A 9-10 cells cross-build identical;
block B 8 cells byte-identical on data@bases + bartcore.bcf +
draws; block C 30+ named refusals incl. the nested-chain pair;
subset-semantics draws identical to pre-subsetted fits; trio
37+12+10 bitwise canonical; mutation isolates exactly the ancestor-
chain/desugar/RHS-subset cells; NEWS 258; check OK). The final PASS
came from a respawned runner re-running the WHOLE battery after a
usage-limit kill mid-mutation; its audit found and restored the
predecessor's half-applied mutation before rebuilding its libs from
scratch.

S13 (fork 6c, formula-path bases subsetting) LANDED 172523e6,
2026-08-17. The decided BEHAVIOR CHANGE: on dbarts()'s formula path
with subset present, a forests= basis (raw matrix or formula-
evaluated) must be FULL-DATA length and is subset by the same index
the model frame uses (resolveFormulaBasisSubset reads the count from
the formula's own variables and the index by model.frame's subscript
reading; alignForestBasisToSubset applies it); the equal-count-but-
not-full-length shape a count-only check once silently row-aligned
is refused naming the forest and both counts; any OTHER count still
falls through to the pre-existing length-of-basis error.
Implementation finding: basis FORMULAS were evaluated eagerly
against the full pre-subset data all along, so under subset the old
count check made every formula basis UNUSABLE (always errored) - the
slice makes exactly the full-length shape work, for matrices and
formulas alike, converging on the term route's post-subset rule.
Scope containment verified: validateForestBases untouched at all
three call sites; the x/y interface ALREADY subset full-length bases
correctly (the asymmetry was formula-path-only - regression-
asserted); direct dbartsData(bases =) keeps its count-only contract
BY SCOPING, regression-asserted - the ambiguous shape survives on
that lower-level constructor, recorded for the arc review; the S12
term route is bitwise unchanged. NEWS carries the behavior change.
Gates: both runners green (5581/0; trio 37+12+10 bitwise canonical;
the six-cell behavioral matrix verified base-vs-slice on matrices
AND formulas; edge probes: all-rows-subset identity, neither-count
fallthrough, negative-index alignment vs a manual fit; mutation
isolates exactly the 2 refusal assertions; NEWS 259; check OK
twice).

S14 (docs) LANDED 9031b348, 2026-08-17. man/bart2.Rd split out of
man/bart.Rd (363 lines, ~210 moved family/Value text + ~150 written:
the Formula Terms section, the reworked argument items, one runnable
example); bart.Rd 609 -> 434 keeping bart and the shared generics;
every base alias still resolves (124/124 audited) and ?bart2 is its
own topic. _pkgdown.yml gains bart2; pkgdown clean. The d2
n.samples boundary documented on both sides (dbartsControl's
per-run return count vs the fitters' thinned sweep budget).
forest.Rd seealso points at the full term grammar. Accuracy-pass
findings fixed while splitting: offset.test was documented as NULL
but has always tracked offset (NEWS bug-fix entry); the bart.Rd dots
item claimed dbartsControl forwarding (actually plot/getTrees); the
Value section's "forestFits not reachable from bart2" sentence was
stale since the term route landed. Verified already-done, not
duplicated: varianceForest.Rd, rbart.Rd's k item, T-C-is-a-test.
Gates: both runners green (scope man/+pkgdown+NEWS only; tinytest
5581/0 unchanged; usage-vs-formals audit exact on
bart2/rbart_vi/xbart/dbarts/dbartsControl/forest/varianceForest,
bart's 11 integer-literal cosmetic mismatches proven pre-existing
byte-identical to base; content spot-checks vs the landed code found
no contradiction; pkgdown clean; NEWS 260-entry parse; check
--as-cran Status OK with the implementer's vignette-warning claim
REFUTED by a byte-identical base-vs-slice log diff - its own
library scope lacked the vignette toolchain). ARC SLICES S1-S14
COMPLETE; closure records and the comment sweep follow.

## ARC CLOSURE, 2026-08-17

S1-S14 all LANDED, each independently gated by an implementer plus an
independent runner, every gate CI six-green: S1 ebc57af7, S2 2eb6a1a4,
S3 6f1ba79e, S4 ec59caf6, S5 ba4f761d, S6 5aa4d292, S7 de2212bc,
S8 d54c2791, S9 8b263ae5, S10 726dab10, S11 00abf336 (S9+S10+S11 landed
as a batch under the shared merged-tree battery clause, own-base gated
first, merged tree re-verified at S11's tip - batch note above), S12
4b179585, S13 172523e6, S14 9031b348. The comment sweep dcc8262e closed
the arc's scaffolding per the owner's comments-state-constraints rule
(records may cite docs/plans; shipped code may not) - it swept both this
arc's own plan/slice tags (S4-S11 tags in R/utility.R, R/bart.R,
R/rbart.R, R/xbart.R, T-A..T-E citations in the tinytest files) AND
prior arcs' residue left standing in the same files it touched:
multiforest-extension-surface's M4.1/M4.3/M4.4 comments, binary-kforest-
prior-default's leg-(c)/fork-3b comments, and bcf-b-ridge.md's derivation
citations (its exponent rule, its b-move rationale) in R/bartcore.R and
the bartcore headers.

Lockstep consumer state: bartCause dbarts-1.0 = 31519d22 (S8's xbart
control= removal - bayesOpt.R's move-into-control block replaced by a
formals()-driven drop filter, own probe pair, suite green), which is
also bartCause's current dbarts-1.0 tip. The S5-era pairs already
recorded above: bartCause a3aa98d (rngSeed -> seed at bcf.R + two test
dbartsControl calls), stan4bart 54e157b (rngSeed -> seed at test-02 plus
a NEWS breaking-changes bullet).

Items the arc RECORDED for later rather than resolving:

- `dbartsData(bases =)`'s count-only contract at the lower-level
  constructor stays ambiguous at equal-but-not-full row counts, by
  SCOPING rather than by fix - the formula path's own ambiguity closed
  at S13, this one did not (S13 note above).
- Forest 1's `sd` and `update.amplitude` have no term spelling (5.7):
  the cost of VD's XOR against a `forests =` formal, which remains
  addable later and is recorded as a door rather than taken (5.5.3).
- The K-forest front-door spelling door `docs/plans/bcf-bartcause-
  relocation.md`'s Doors held section left riding this arc's review is
  now ANSWERED: formula-level terms with colon sugar, no `forests =`
  formal (5.2, 10.11) - that plan's own 2026-08-16 post-arc addendum
  already records this against section 5 of this plan, so nothing
  further is owed there.
- Out-of-sample per-forest replay (N8) stays a door: S11 shipped only
  the in-sample half (`extract(type = "forest")`); `predict()` needs the
  engine-side per-forest saved-tree replay door recorded at
  `docs/plans/archive/bcf-bartcause-relocation.md`:1271-1273, which this arc did
  not open. 5.8.6 says that door's language is updated "at the records
  commit" - recorded here rather than by editing that file, which is
  outside this closure's scope.
- The harness coverage note: the equivalence harness has no bart2
  gaussian/probit scenario and no two-forest-via-formula scenario, so
  several slices carried their own in-slice bitwise A/B gates instead;
  a harness coverage note rides the next baseline re-record (1.3's own
  text, unchanged by this closure).

One recommendation surfaced by the sweep: a handful of door memos in
`docs/plans/` carry rationale that SHIPPED code still relies on
(`docs/plans/archive/bcf-b-ridge.md`'s exponent rule was cited from a live
comment until this sweep removed the citation, not the reliance) -
these read more like standing design records than closed proposals.
Promoting the load-bearing ones to `docs/design/` (or cross-referencing
them from the relevant design doc) would give them the durability
`docs/design/*.md` is for; route this to the release-candidate-review
program rather than deciding it here.
