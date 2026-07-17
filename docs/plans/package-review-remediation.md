# package-review-remediation

scope: remediation of the 2026-07-17 seven-reviewer package review (engine
  correctness, bridge/memory safety, performance, R surface + tests, docs,
  writing, fresh-eyes user), run at VD's request, addressed at orchestrator
  discretion (VD 2026-07-17). The panel's verdict: the engine and bridge are
  exceptionally strong (zero high-severity engine findings; the delicate MH
  arithmetic was hand-cleared - change-move proposal correction, swap
  involution, BCF interweaving, rollback machinery); the defects live at the
  R-surface boundary and the release-documentation surface.

decisions adopted under the discretion grant (VD may overrule):
  - family = "auto" DETECTS factor responses (2 levels -> probit, 3+ ->
    multinomial) and prints a one-line verdict, rather than refusing or
    (the current bug) silently fitting gaussian on level codes.
  - keepTrees fits store sampler state automatically before serialization
    (no more undocumented fit$fit$state ritual).

## Arcs

A - pre-release correctness (this doc's C1-C4). B - release-surface docs
(C5). C - post-release performance + latent lows (recorded as TODO items,
not implemented in this arc).

## Findings and disposition

### Tier 1 - verified correctness bugs (orchestrator-reproduced)

T1a. rbart_vi ranef shape + predict. fit$ranef is samples x OBSERVATIONS
  (verified: 10 groups, 200 obs -> ranef.mean length 200, dim(ranef)
  20 x 200) where docs promise per-group; predict's level-matching against
  dimnames therefore never succeeds and it silently re-draws intercepts
  from N(0, tau) even for TRAINED groups (verified cor(fitted, predict)
  0.93 on identical data; reviewer saw 0.72 with a "levels not present in
  training" warning on trained groups). Also: predict(group.by =) crashes
  ("subscript out of bounds") on numeric/character input that fitting
  accepts. FIX (C1): per-group ranef/ranef.mean as documented; predict
  looks up trained groups exactly (round-trip cor ~1 up to tree-recompute
  noise) and draws from the prior ONLY for genuinely new levels, naming
  them; group.by coerced as at fit time. Draw-neutral (output shape/
  assembly only) - but if the equivalence fixtures record a ranef channel
  in the buggy shape, re-record per the demote-with-neutrality-trail
  precedent (underlying values must reproduce, reshaped).

T1b. Formula + factor response hard-errors ("'range' not meaningful for
  factors", after a "type = numeric ignored" warning) in bart2/dbarts/
  rbart_vi, while help promises probit and the x/y interface delivers it.
  FIX (C2): the formula path extracts the response without numeric
  coercion and routes factors as the x/y path does.

T1c. family = "auto" with a 3+-level factor via x/y silently fits gaussian
  on integer level codes (verified: class bart, yhat over code range).
  FIX (C2): per the adopted decision, auto detects (2 -> probit, 3+ ->
  multinomial, printed verdict); logical and 2-level character responses
  coerce to probit (both currently rejected: "when 'formula' is numeric,
  'data' must be numeric as well").

### Tier 2 - surface/contract defects (reviewer-verified live)

T2a. Multinomial generics layer (C3): residuals() falls to default and
  silently returns NULL; plot() falls to plot.default ("'x' is a list...");
  summary() falls to summary.default (structure dump); predict has no type
  arg (no PPD categories for newdata; extract's ppd covers stored channels
  only); print.bartMultinomial reports the COMBINED draw total as "per
  chain" under combineChains = TRUE (the default) - fitSynopsis already
  divides by n.chains correctly, route through it or divide.

T2b. Serialization: keepTrees + saveRDS + predict errors with the
  storeState ritual message (docs say fit$state, message says
  storeState() - also inconsistent). FIX (C4): auto-store per the adopted
  decision.

T2c. Error quality (C4): Inf in response -> "missing value where
  TRUE/FALSE needed" (NaN gets a clean message); 0-row data -> "subscript
  out of bounds"; n.samples = 0 -> "dim(X) must have a positive length".
  Validate up front with named errors. Unnamed-matrix newdata matches by
  position with only a warning while returning badly wrong numbers when
  columns are swapped - escalate when the fit had named predictors
  (implementer judgment on error vs stronger warning; do not break the
  legitimately-unnamed-trained case). Missing newdata variable named like
  a base function -> "invalid type (builtin) for variable 'c'" - name the
  missing variable instead.

T2d. Bridge (C4, C-side): assignInPlace (R_interface.cpp:77-80,122)
  writes R-supplied index/memcpy with no bounds/type/length validation on
  the per-iteration rbart hot path - add checks (defense-in-depth; callers
  currently valid). rc_get*At off-by-one (rc/bounds.c:55,271,487: "i >"
  should be "i >="; unused API, latent). simd.c:198 AVX512 CPUID max-leaf
  constant 0x60 should be 0xD (dead path today - no AVX512 kernels
  dispatch). C_interface.cpp:158-192 GetRNGstate/PutRNGstate bracket not
  restored on Rf_error longjmp - unwind-protect or drop where engine uses
  only ext_rng. createHolder/bartcore_run/unwindProtect OOM-and-longjmp
  leak windows (R_interface_bartcore.cpp:1545-1572, 2340-2365, 224,
  1844-1850) - tighten per the reviewer's RAII suggestion where cheap.
  simd runtime re-dispatch documented as init-only/main-thread (comment).

### Tier 3 - release-surface docs (C5)

NEWS (inst/NEWS.Rd) + README omit family = "aft" AND family =
"multinomial" entirely (NEWS last touched 387 commits before they landed).
bart.Rd Usage omits "multinomial"; description never mentions bart2's
family surface; xbart's narrower family set unnoted. dbarts.Rd claims
sparseFactor construction is refused (x/y path accepts it - verified;
formula path refuses). xbart.Rd offset wrongly binary-only (code passes
through unconditionally). dbartsSampler-class.Rd missing getSigmas/
getLatents/getSumsOfSquaredResiduals (NEWS references the latter as
documented). docs/design/survival.md still "Status: proposed" though
shipped (f0efc03 etc.) - update to LANDED like multinomial.md. Typos:
"Regreesion" (dbarts.Rd:5), "continous" (bart.Rd:110), "reponses" +
"or the a" (bart.Rd:125), "interpretted" (bart.Rd:270). pdbart.Rd stale
author emails (sync to DESCRIPTION). Writing batch: the wrong duplicated
family comment (R_interface_bartcore.cpp:1037-1039, 1512-1515 - omits
aft; fix first site, point second at it); R/dbarts.R:59 cryptic defaults
comment (state the constructor-vs-prototype n.threads divergence);
R/A_class.R:1 joke -> state the load-order invariant; delete debug
leftovers (R/multipleAssignment.R:48-145 commented cats,
R/sliceSample.R:72-75 abandoned branch); sampler.hpp:884/938 duplicated
rollback comment (shrink second to a pointer); dbarts.h x_test/
offset_test -> xTest/offsetTest (parameter names are not ABI; keep the
X-macro hash re-bake + MINOR bump discipline in mind - a parameter RENAME
does not change signatures/hash input? verify: the hash covers stringized
signatures INCLUDING parameter names, so a rename re-bakes the hash and
bumps MINOR with the stan4bart floor moving in lockstep at release; if
that cost is judged disproportionate, defer the rename to the next
planned hash-moving change and note it).

### Tier 4 - engine notes (all low; C-arc/TODO unless trivial)

E1. Multinomial level-centering precision (combiner.hpp:735-745):
  observation-count vs exact leaf-count precision - confirm the exact
  gate exercises unequal per-forest shrinkage, else switch. Recorded for
  investigation, not blind change.
E2. Constant-leaf centered-SS cancellation (model.hpp:117-118) -
  documented tradeoff, no action; Welford path recorded as speculative.
E3. combinedFits re-implements softmaxLocationMajor inline
  (combiner.hpp:668-682) - route both paths through the shared helper
  (bitwise-neutral if arithmetic identical; verify).
E4. Scan-order permutation drawn from chain 0's stream (facade.hpp:405) -
  document as intended.
E5. GaussianResponse::rescale zero-row OOB (model.hpp:2074) - guard or
  note upstream reject.

### Tier 5 - performance opportunities (TODO items, post-release)

P1. Constant-leaf treeFits scatter elimination: per-leaf mu + per-obs
  leafOf byte, incremental maintenance in moves; kills the recorded
  ~18.5% setIndexedVectorToConstant share and shrinks the hot slab 8x
  (the threading no-go's revival precondition). Claimed rng-neutral;
  distinct from the killed block-fusion (which kept the scatter).
  Constant-leaf only. Est. 6-12% net, memory-wall-capped.
P2. Drop sumWeightedResponseSq from the leaf suffstat (tree.hpp:500-523,
  moments.c, model.hpp:117): cancels in every MH ratio (scan.hpp already
  proves and exploits this by passing 0); bitwise-identical mechanical
  change, 1-3% + node shrink.
P3. Multinomial margins O(nK^2) -> O(nK): running log-sum-exp maintained
  across the interleaved loop (prefix/suffix LSE, NOT subtract-exp -
  cancellation when excluded category is argmax). Posterior-defining
  (rounding shift) -> exact-gate re-record. Found independently by two
  reviewers. Plus K-contiguous fit locality (P3b) and a shared exp cache.
P4. Birth-move child-by-subtraction + partition/suffstat gather fusion
  (tree.hpp:770-783, 758-767): 1-3% each, rng-neutral, order-careful.

## Commits

C1. rbart_vi ranef/predict fix (Opus). Gates: full tinytest; equivalence
  anchors bitwise EXCEPT any rbart ranef channel recorded in the buggy
  shape - if so, re-record with the neutrality trail; a new tinytest
  asserting per-group shapes, the predict round-trip (cor ~1 on trained
  groups), the new-level prior-draw path naming levels, and group.by type
  coercion at predict.
C2. Factor-response ingestion (Opus): formula factor routing, auto
  detection + verdict, logical/character coercion. Gates: suite + all
  anchors bitwise (existing paths untouched); tinytest for each new
  routing and the verdict line.
C3. Multinomial generics (Sonnet): residuals/plot/summary/print/predict
  type. Gates: suite; fixture bitwise; tests per method incl. the
  print per-chain arithmetic.
C4. Error-quality + serialization + bridge lows (Opus, R + C). Gates:
  suite + anchors bitwise; tests for each named error; serialization
  round-trip test without the ritual; tests/cpp + preclean where C
  touched.
C5. Docs + writing batch (Sonnet). Gates: R CMD check Rd checks; suite
  untouched; the dbarts.h rename ONLY with its hash/MINOR implications
  resolved per Tier 3.

Each commit lands with its own landing note appended here.

## Landings

### C1 (2026-07-17, ca56faf)

The diagnosis inverted T1a's hypothesis: the ranef assembly and the C
channel were always per-group and correct. The root cause was symbol
resolution - data[[which.max(names(data) == sym)]] binds to column 1
when the symbol is absent (which.max of all-false is 1), so a
standalone group.by grabbed the response and manufactured per-value
groups; predict's missing coercion was an independent second bug.
Fixed by-name at all three sites (idiom now eradicated repo-wide,
grepped); predict coerces to factor keeping unused levels (deliberate:
out-of-sample groups on a supplied factor still get prior draws) and
names defaulting levels in the warning. Round-trip cor 0.93 -> 1.0
exactly. Equivalence: case 1 - all 22 scenarios identical (the
fixtures pass group.by via $, never taking the buggy branch), BCF and
multinomial bitwise, NO re-record. Suite 2979/0 with a new
test-rbart-groupby.R. One pre-existing test expectation updated: it
had asserted the buggy silent-grab's downstream error for an
undefined group.by; the fix errors "'group.by' not found" instead.
Not touched: the fit-time group.by.test warnings do not name levels
(scoped to predict; fixture message-matching risk noted).

### C2 (2026-07-17)

Response-type routing rooted at the data layer: dbartsData records the
response's original type and level count in two new prototype-defaulted
slots (by family-resolution time y was already coded doubles, so the
information had to be captured at ingestion), with shared classify/code
helpers whose numeric path is byte-identical as.double. Formula
extraction drops the "numeric" coercion that crashed on factors. Routing
per the adopted decision: 2-level factor/logical/2-level character ->
probit with a one-line verdict everywhere; 3+-level -> multinomial in
bart2 (auto == explicit, bitwise) and an informative
requires-bart2-multinomial error in dbarts/bart/rbart_vi/xbart;
explicit contradictions (gaussian + factor) error instead of fitting
codes; ordered factors route as plain with the verdict noting it.
Every reproduction bitwise (formula == x/y, logical == 0/1-numeric,
auto == explicit); all three anchors identical; suite 3004/0 with a
new test-factor-response.R (24 assertions). Two tests that had
asserted the bugs updated. Rd accuracy touched in four topics (C5 owns
prose). Noted, untouched: count matrices stay explicit-only; the
auto-multinomial + weights two-step error; the pre-built dbartsData
generic error path.

### C3 (2026-07-17)

The five multinomial generics: print divides by n.chains directly (a
recorded deviation from routing through fitSynopsis, whose
keepTrees-FALSE fallback would misread the K-widened varcount margin);
residuals returns the n x K observed-proportion-minus-fitted matrix
(indicator at n_i = 1, y/rowSums for counts) - the direct
generalization of residuals.bart, never silent NULL; plot traces
per-category mean predicted probability (the closest analog of
plot.bart's sigma trace for a family with no residual scale); summary
mirrors summary.bart's shape over the same channel; predict gains
type = c("ev", "ppd") with ppd factored into a shared helper so
predict and extract draw categories identically. Default paths
RNG-neutral (asserted on .Random.seed) and bitwise-unchanged; fixture
identical; suite 3034/0 with 30 new assertions.

### C4 (2026-07-17)

Serialization: keepTrees fits auto-store sampler state at fit-object
assembly across bart/bart2/rbart_vi (byte-identical to the manual
ritual's payload; storeState() is now the one documented name). SIZE
CAVEAT recorded: the state is ~2.8x yhat.train (scales with trees x
samples x chains), an eager in-session cost accepted for ritual-free
round-trips - VD may prefer lazy, reversible. Error quality: named
errors for non-finite response, zero-row data, non-positive
n.samples, and missing test variables (previously a builtin-collision
message); the positional-matching warning spells out the assumed
column mapping. Bridge: assignInPlace fully validated (spec's
assignInPlaceDiagonal does not exist - one function, both modes
hardened); the three rc bounds off-by-ones fixed; the AVX512 CPUID
constant corrected; the C API's four GetRNGstate brackets DROPPED
after verifying every chain rides a private MT seeded once at
creation (test-capi bitwise green); simd re-dispatch documented
init-only. OOM-leak windows deliberately deferred per plan. Gates
independently re-run: preclean install, tests/cpp clean-build all
pass, tinytest 3050/0 (+16), equivalence 22/22 + bcf + multinomial
bitwise. Also folded: air format across the tree (three C2-era files
had drifted, failing the lint workflow - VD request).
