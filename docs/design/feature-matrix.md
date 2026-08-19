# Response-model feature matrix

Status: living reference, anchors current at ad4a131b (2026-08-19). The
"current at" commit stamps the ANCHORS and nothing else: the 2026-08-14 pass
relocated every anchor by symbol and re-adjudicated no cell value, the
2026-08-15 pass did the same, scoped to the files the
`binary-kforest-prior-default` arc's three landing commits touched (every
other anchor carries over from 2026-08-14 unmoved, since its host file is
untouched by the arc), the 2026-08-15 `adoption-slate` pass did the same
again, scoped to the files that arc's eight landing commits touched, the
2026-08-18 pass did the same as a whole-file walk, and the 2026-08-19 pass
did another whole-file walk for the pre-RC resync. Carries
no landing date and is not a design proposal - the orchestrator updates it in
place at every landing that changes a cell, and VD uses it to schedule
feature completion.

What each shipped response model can and cannot do, one row per model and one
column per capability that bears on scheduling. Every SHIPPED and REFUSED cell
carries an anchor verified against the live tree; a wrong cell misdirects
scheduling, so a cell that cannot be verified is marked `?` rather than
guessed.

## Legend

| code | meaning |
|---|---|
| `S` | SHIPPED. Works today; the anchor is the site that makes it work. |
| `R` | REFUSED on model or identification grounds - the refusal is part of the model, not a hole. The anchor is the refusal site. |
| `P` | PLANNED. A named design arc covers it; the slice is named in the cell. |
| `M` | MISSING. Not built, no schedule. An anchor, when given, is a guard that errors *because the thing is unbuilt* - a recorded door, not a model refusal. |
| `-` | N/A. The concept does not apply to this row; the row footnote says why. |
| `?` | UNVERIFIED. Constructs today with no refusal site, but no test, doc or adjudication backs it. Do not schedule against the cell until it is settled; every `?` is listed under "Gaps". |

Path aliases used in anchors (line numbers are at ad4a131b):

    RIB   src/R_interface_bartcore.cpp      CAPI  inst/include/dbarts/dbarts.h
    MOD   src/bartcore/model.hpp            CH    src/bartcore/chain.hpp
    FAC   src/bartcore/facade.hpp           COM   src/bartcore/combiner.hpp
    MOV   src/bartcore/moves.hpp            SAM   src/bartcore/sampler.hpp
    bart.R, dbarts.R, spec.R, rbart.R, xbart.R, data.R, generics.R,
    A_class.R, bartcore.R      -> R/<name>
    sampler.Rd -> man/dbartsSampler-class.Rd     bart.Rd -> man/bart.Rd
    bart2.Rd -> man/bart2.Rd

Anchor pass on record, 2026-08-14 (root TODO `feature-matrix-anchor-refresh`,
closed by it). EVERY anchor in this file was re-located BY SYMBOL against
55cc1756: target file opened, cited symbol found, number moved only where the
symbol had moved. No anchor was derived by applying an offset to another. The
two whole-column faults M4.4 recorded both reproduced - the section-1 flat-C
cells indexed the `DBARTS_C_API_LIST` X-macro body rather than the family prose
on `dbarts_sampler_create`, and the MOD column had not followed `model.hpp`
through M4.1-M4.3 - and both are fixed here, as are the anchors M4.5 and M4.4
left standing outside the `bcf` row. Those two partial-pass exceptions are
therefore superseded and are no longer carried. Two numbers are HISTORICAL
rather than live and were left as written: [f36]'s "the earlier CH:1589" and
[f23]'s "the earlier `spec.R:440`", each naming a superseded anchor rather than
a site. Cell values were NOT re-adjudicated - this pass moved anchors only, so
a cell whose VALUE was wrong before it is wrong still.

Anchor pass on record, 2026-08-15, SCOPED to the `binary-kforest-prior-default`
arc's three landing commits (4dbf2dbc, 0faeb416, e623fbf3; branch tip
63524e5e), not a whole-file walk. Only anchors targeting the files those three
commits touched were re-located by symbol against 63524e5e: CH, CAPI, RIB,
`C_interface.cpp`, `dbarts.R`, `bartcore.R`, `model.R`, `spec.R`, `A_class.R`,
COM, and the test files the arc edited. MOD (`model.hpp`) was NOT touched by
the arc and was not re-walked. Every anchor outside that set carries over from
the 2026-08-14 whole-file pass unchanged - its host file is byte-identical
between 55cc1756 and 63524e5e, so 55cc1756's numbers still hold. Cell values
were again NOT re-adjudicated.

Anchor pass on record, 2026-08-15, SCOPED to the `adoption-slate` arc's eight
landing commits (5bedf923, da3c76f9, eeedc07c, 890efd3d, 37d9ec81, 55811082,
cedf4c34, 13393350; branch tip c05322a8), not a whole-file walk. Only anchors
targeting the files those eight commits touched were re-located by symbol
against c05322a8: RIB, CAPI, MOD, CH, FAC, SAM, `bart.R`, `dbarts.R`,
`bartcore.R`, `bart.Rd`, `sampler.Rd`, `C_interface.cpp`, TODO,
`tests/cpp/test_model.cpp`, `tests/cpp/test_sampler.cpp` and
`inst/tinytest/test-capi.R`. COM (`combiner.hpp`), MOV (`moves.hpp`),
`spec.R`, `rbart.R`, `xbart.R`, `data.R`, `generics.R` and `A_class.R` were
NOT touched by the arc and were not re-walked. Every anchor outside that set
carries over from the 2026-08-15 `binary-kforest-prior-default` pass
unchanged - its host file is byte-identical between 63524e5e and c05322a8, so
63524e5e's numbers still hold. Cell values were again NOT re-adjudicated.

Value audit on record, 2026-08-17, SCOPED to the `bart2-argument-consolidation`
arc's fourteen landing commits (ebc57af7..9031b348, S1-S14) plus the dcc8262e
comment sweep - unlike the two anchor-only passes above, this pass
RE-ADJUDICATED (not just re-located) the cells the arc could plausibly have
changed: `bart()`/`bart2()`/`xbart()`/`rbart_vi()` argument surface, prior
objects, variance forest spelling, the per-forest output channel, formula
terms, multi-forest declaration routes. Findings and fixes: `bart()` gained a
narrow `family` formal (S10), which changed section 1's `bart()` column at
logistic and aft (R/`?` to `S`) and at multinomial/nbinom/hurdle (`-` to `R`,
the concept now applies and is refused rather than absent), plus footnotes
f1/f5 and the logistic/aft Gaps entries; `bart2()`'s formula interface (S12)
now reaches the bcf row's general K-forest amplitude capability through a
`forest()` term, changing that row's `bart2()` cell, footnote f7, the bcf Gaps
paragraph, and adding a note after the section-1 table and to footnote f48
(the in-sample per-forest output channel, S11). Prior objects (S7), variance
forest spelling (S6, the `variance =` argument name itself is unchanged - only
its former sub-formals collapsed into `varianceForest()`) and S13's
formula-path bases-subsetting behavior change introduced no stale claims found
in this file. This pass did NOT re-walk every anchor into the arc's touched
files (`bart.R`, `dbarts.R`, `rbart.R`, `xbart.R`, `spec.R`, `model.R` shifted
substantially across fourteen commits touching them repeatedly) - unlike the
two passes above, it is not a line-accurate anchor refresh, and the `Status`
line above is deliberately NOT advanced to this arc's tip; a scoped anchor-only
pass over these files, on the `binary-kforest-prior-default`/`adoption-slate`
pattern, remains owed.

Anchor pass on record, 2026-08-18, a WHOLE-FILE walk against 8e1e674c over
the full alias namespace - the line-accurate refresh the value audit above
recorded as owed, extended to every anchor in the file. Every anchor was
re-located by symbol: target file opened, cited symbol found, number moved
only where the symbol had moved, with range endpoints and comma-list
elements resolved independently. No anchor was derived by applying an offset
to another or reading a delta off a diff. MOV (`moves.hpp`) carries over
unchanged, byte-identical between c05322a8 and 8e1e674c - the established
carry-over rule - as do the anchors into the byte-identical
bcf-public-surface.md, sbc-family-tiers.md, sbc-calibration.md,
variance-forest-mutation-routing.md and tests/cpp/test_model.cpp (where
[f15]'s `:5449` and `:5541`, stale at c05322a8 itself, are corrected to the
symbols' true lines). The `bart2.Rd -> man/bart2.Rd` alias was added to the
table above (man/bart2.Rd, split out of man/bart.Rd at
`bart2-argument-consolidation` S14, postdated the table); NO anchor changed
file, to it or otherwise - every claim's cited content still lives in its
cited file, bart.Rd:151's probit auto-detection included (now bart.Rd:76).
The two historical numbers - [f36]'s "the earlier CH:1589" and [f23]'s "the
earlier `spec.R:440`" - were left as written. Cell values were NOT
re-adjudicated - this pass moved anchors only, so a cell whose VALUE was
wrong before it is wrong still.

Anchor pass on record, 2026-08-19, a WHOLE-FILE walk against ad4a131b (the
pre-RC resync). Every anchor was re-located BY CONTENT rather than by
symbol-grep alone: the construct cited at 8e1e674c was read from that
commit's blob, then found in the live tree by its own content (a function
signature, an assert message, a branch condition, a doc comment, a field),
with range endpoints and comma-list elements resolved independently -
several turned out to move by a DIFFERENT amount than their neighbor on the
same line, and one (`generics.R`'s `hurdleSigmaVec` doc comment) grew from
eight lines to eleven, which a same-length assumption would have missed. No
anchor was derived by applying an offset to another or to a whole block.
`src/R_interface_bartcore.cpp`, `src/bartcore/model.hpp` and
`src/bartcore/chain.hpp` (stale since the hetero-state-carry landing),
`inst/include/dbarts/dbarts.h` (a struct-layout/ABI-enum/callback reshape),
`src/C_interface.cpp`, `src/bartcore/combiner.hpp`,
`src/bartcore/sampler.hpp`, `bart.R`, `generics.R`, `data.R`, `TODO` and
`inst/tinytest/test-capi.R` all moved at least one cited anchor, by deltas
from 0 up to +175 depending on where in each file the construct sits -
never a uniform per-file offset. `tests/cpp/test_sampler.cpp` moved two of
its four cited sites and left two unchanged. `tests/cpp/test_model.cpp`'s
three cited test functions and `man/bart.Rd`'s, `man/bart2.Rd`'s,
`man/dbarts.Rd`'s and `man/dbartsSampler-class.Rd`'s single cited lines
each reproduced UNCHANGED despite their host files changing elsewhere -
the last two took in-place text expansions AT other lines, and
`man/dbartsSampler-class.Rd:151` itself grew in place without moving.
Two continuation forms the earlier passes' extractor regex cannot see were
resolved by hand this time: a second cell number following a first that
already carries the alias token (e.g. `setPredictor`'s per-observation
variant, the second number of the old RIB 4938/5094 cite), and
a bare parenthetical with no adjacent token at all (the section-5
`test-capi.R` attribute cites `(:603)`, `(:220)`, `(:1266)`, and
`hurdleSigmaVec`'s "definition :N" form). `spec.R` and `rbart.R` are
CHANGED files (a reworded comment and a `$` to `[[` rewrite, respectively,
neither touching a cited line) whose every cited anchor nonetheless
reproduced byte-identical at its old line - reverified by content, not
assumed from the small diff stat. Unchanged and carried over byte-identical
since 8e1e674c: MOV (`moves.hpp`, which changed substantially this cycle
but carries zero live line anchors in this file to relocate), FAC,
`dbarts.R`, `xbart.R`, `A_class.R`, `bartcore.R`, `model.R`,
`formulaTerms.R`, `tests/cpp/test_shape.cpp`, `benchmarks/R/equivalence.R`,
`inst/tinytest/test-active-rows-pins.R`, `test-calibration-prior-draws.R`,
`test-calibration-midchain.R`, `test-bcf-family.R`, `test-formula-terms.R`
and the cited design/plan docs other than
`variance-forest-mutation-routing.md` and `empty-leaf-veto.md` (the latter
changed 143 lines but is cited only by section title, no line to relocate);
within CHANGED files,
`benchmarks/baselines/MANIFEST`'s three cited lines ([f39]) also
reproduced unmoved. No CONTENT mismatch was found - every cell's cited
construct still matches the claim it anchors. The two historical numbers -
[f36]'s "the earlier CH:1589" and [f23]'s "the earlier `spec.R:440`" - were
again left as written. Cell values were NOT re-adjudicated - this pass
moved anchors only, so a cell whose VALUE was wrong before is wrong still.

## Rows

Thirteen rows. The first ten are response models proper; the last three are
couplings or decorations over a base response that a user selects the same way
and schedules against the same way, so they earn rows.

| key | model |
|---|---|
| gaussian | Gaussian (`ResponseFamily::gaussian`, `GaussianResponse` MOD:2777) |
| student | Gaussian + Student-t residuals (`resid.dist = student()`, `TResponse` MOD:4015) |
| probit | Binary probit (`ProbitResponse` MOD:3052) |
| logistic | Binary logistic, weights = observation counts (`LogisticResponse` MOD:3511) |
| ordinal | Ordered categorical, cumulative probit (`OrdinalResponse` MOD:3191) |
| nbinom | Negative binomial, positive-integer dispersion (`NBResponse` MOD:4312) |
| multinom | Multinomial softmax, K forests (`MultinomialResponse` MOD:3690 + combiner) |
| aft | AFT survival, log-normal (`AFTResponse` MOD:3768) |
| hazard | Discrete-time hazard (person-period sugar, dbarts.R:487-536) |
| hurdle | Hurdle / two-part semicontinuous (R-side composition, bart.R:2206) |
| bcf | K-forest amplitude family, bcf's two forests being its K = 2 instance (`BCFForestCombiner` COM:743) |
| grouped | Grouped random intercepts (`GroupedResponse` MOD:4680) |
| hetero | Heteroscedastic variance forest (CH:724) |

The engine's `ResponseFamily` enum has only six tokens (MOD:2580: gaussian,
probit, logistic, aft, ordinal, nbinom); student, hazard, hurdle, bcf, grouped
and hetero are all reached some other way, which is exactly why they need rows
here rather than an enum read. That enum now REACHES the bcf row, which is why
so much of that row is family-dependent below: since M4.4 the K-forest chain
selects its response model off `BCFSpec::family` (COM:323, `switch (spec.family)`
CH:755) instead of building an unconditional `GaussianResponse`, and the
K-forest `Sampler` constructor takes `family_(spec.family)` (SAM:160) instead of
pinning gaussian. Leaf models (constant, monotone, linear, GP) are
an orthogonal axis and are not rows; where a leaf model gates a capability the
cell says so.

## 1. Construction surfaces

| model | `bart()` | `bart2()` | `dbarts()` + R5 | `rbart_vi()` | `xbart()` | flat C `dbarts.h` |
|---|---|---|---|---|---|---|
| gaussian | S bart.R:2721 | S bart.R:671 | S dbarts.R:378 | S rbart.R:49 | S xbart.R:26 | S CAPI:647 |
| student | S bart.R:2591 | S bart.R:685 | S dbarts.R:365 | M rbart.R:62 | M xbart.R:2-33 | S RIB:2558 [f2] |
| probit | S bart.Rd:76 | S bart.R:672 | S dbarts.R:379 | S data.R:555 | S xbart.R:26 | S CAPI:645 |
| logistic | S bart.R:2594 [f1] | S bart.R:673 | S dbarts.R:380 | R rbart.R:49 | S xbart.R:26 | S RIB:1579 |
| ordinal | R bart.R:2603, 2634 [f1] | S bart.R:676 | S dbarts.R:382 | R data.R:500 | R data.R:500 | S RIB:1587 [f3] |
| nbinom | R bart.R:2603 [f1] | S bart.R:677 | S dbarts.R:383 | M rbart.R:49 | M xbart.R:26 | S RIB:1594 [f3] |
| multinom | R bart.R:2603 [f1] | S bart.R:675 | R dbarts.R:376 | R data.R:500 | R data.R:500 | M [f4] |
| aft | S bart.R:2594 [f1] [f5] | S bart.R:674 | S dbarts.R:381 | S rbart.R:49 | M xbart.R:26 | S CAPI:648 |
| hazard | - [f1] | S bart.R:678 | S dbarts.R:384 | M rbart.R:49 | M xbart.R:26 | M [f6] |
| hurdle | R bart.R:2603 [f1] | S bart.R:681 | R dbarts.R:445 | M | M | M [f6] |
| bcf | - [f1] | S R/formulaTerms.R (ingestFormulaTerms) [f7] | S dbarts.R:371 | R rbart.R:62 | M | S CAPI:668-690 |
| grouped | - [f1] | M [f8] | M [f8] | S rbart.R:388 | M | S RIB:1926 [f3] |
| hetero | - [f1] | S bart.R:664 | S dbarts.R:370 | R rbart.R:62 | M | S RIB:2016 [f3] |

`dbartsSpec()` (spec.R:684-692) resolves the seven single-forest tokens - auto,
gaussian, probit, logistic, aft, ordinal, nbinom - plus BCF through its
`forests =` argument (spec.R:681, `forest(basis = ...)` replacing the removed
`treatment =`, multiforest-extension-surface M2) and a variance forest through
`variance =` (spec.R:680); it does not reach multinomial, hazard, hurdle or
grouped. A `forests =` fit resolves **gaussian, probit or logistic** since M4.4;
aft, ordinal and nbinom are refused there by name, each stating what it is
missing (spec.R:505-527), with the same three-family gate at the bridge
(`refusedBCFFamilyReason` RIB:2281, called from both creation routes at
RIB:2329 and RIB:3137) and at the factory (`createBCFSampler` FAC:822-825).

Since S12 (`bart2-argument-consolidation`), `bart2()`'s formula interface
reaches the same `forests =` machinery through a `forest()` term rather than a
formal of its own - VD decided the term route as an XOR against a flat
`forests =` on `bart2`, recorded as a declined-but-addable door (that plan's
5.5.3) - so this is a THIRD construction route into the bcf row, alongside
`dbarts()`/`dbartsSpec()`'s `forests =` and bartCause's `bcf()`; the family gate
is identical (gaussian, probit or logistic only) since it resolves through the
same `forests =` machinery. See [f7].

## 2. Mutation channels on the R5 `dbartsSampler`

The channels that make dbarts a conditional model inside an outer sampler.
`updateScale` is broken out because it is refused independently of the setter
it rides on.

| model | `setResponse` | `setOffset` | `updateScale = TRUE` | `setPredictor` (+ per-obs) | `setWeights` | `setSigma` | test surface |
|---|---|---|---|---|---|---|---|
| gaussian | S MOD:2810 | S MOD:2890 | S MOD:2890 | S RIB:4989, 5150 | S MOD:2859 | S RIB:2871 | S RIB:4666, 4719 |
| student | S MOD:4084 | S MOD:4091 | S MOD:4091 | S RIB:4989 | S MOD:4108 | S RIB:2871 | S RIB:4666 |
| probit | S MOD:3105 | S MOD:3111 | - [f9] | S RIB:4989 | R RIB:2763 | R RIB:2873 | S RIB:4666 |
| logistic | S MOD:3579 | S MOD:3585 | - [f9] | S RIB:4989 | R RIB:2756 [f10] | R RIB:2873 | S RIB:4666 |
| ordinal | S MOD:3242 | S MOD:3250 | - [f9] | S RIB:4989 | R RIB:2763 | R RIB:2873 | S RIB:4666 |
| nbinom | S MOD:4379 | S MOD:4386 | - [f9] | S RIB:4989 | R RIB:2763 | R RIB:2873 | S RIB:4666 |
| multinom | - dbarts.R:1223 [f11] | - [f11] | - [f11] | - [f11] | - [f11] | - [f11] | - [f11] |
| aft | S MOD:3843 | S MOD:3856 | S MOD:3856 | S RIB:4989 | R RIB:2763 | S RIB:2871 | S RIB:4666 |
| hazard | S MOD:3105 [f6] | S MOD:3111 | - [f9] | S RIB:4989 | R RIB:2763 | R RIB:2873 | S RIB:4666 |
| hurdle | - [f12] | - [f12] | - [f12] | - [f12] | - [f12] | - [f12] | - [f12] |
| bcf | S CH:1045 [f48] | S CH:1045 [f48] | R bartcore.R:305, 326 [f48] | S RIB:4989, 5150 | S RIB:4816 [f48] | S RIB:4591 [f48] | R RIB:2843 |
| grouped | S MOD:4762 [f13] | S MOD:4773 | R RIB:2727 [f13] | S RIB:4989 | S MOD:4788 [f14] | S RIB:2871 [f14] | S RIB:4666 |
| hetero | S RIB:2701 | S RIB:2701 | R RIB:2701 | S RIB:4989, 5150 | S RIB:2696 | R RIB:2869 | S RIB:4666 |

`setData` (whole-data replacement, n free) is single-forest and dense-store
only (RIB:4600-4601) and is refused for grouped (RIB:4602) and aft (RIB:4605);
BCF/multinomial whole-data `setData` stays undesigned by the model-space
survey's verdict (model-space-survey.md doors 1 and 3).

## 3. Row subsetting, latents, calibration

| model | zero-weight row subset | active-rows mask [f15] | `getLatents` | pointwise loglik | nameable calibration [f16] |
|---|---|---|---|---|---|
| gaussian | S sampler.Rd:151, MOD:2785 [f17] | S MOD:2870 | - RIB:5884 [f18] | S generics.R:57 | S dbarts.R:1594, 1599 [f16] |
| student | S MOD:4064-4071 [f17] | S MOD:4119 | S MOD:4131 | ? generics.R:57 [f19] | S dbarts.R:1594, 1599 [f16] |
| probit | R RIB:1626 | S MOD:3095 | S MOD:3131 | S generics.R:102 | S dbarts.R:1594, 1599 [f16] |
| logistic | R RIB:1630 [f20] | S MOD:3531 | S MOD:3605 | S generics.R:102 | S dbarts.R:1594, 1599 [f16] |
| ordinal | R RIB:2549 | S MOD:3227 | S MOD:3272 | M generics.R:129 | S dbarts.R:1594, 1599 [f16] |
| nbinom | R RIB:2555 | S MOD:4359 | S MOD:4412 | M generics.R:129 | S dbarts.R:1594, 1599 [f16] |
| multinom | R RIB:3204 | S COM:1678 [f21] | M MOD:3690 [f22] | M generics.R:129 | R [f23] |
| aft | R RIB:2545 | S MOD:3825 | S MOD:3886 | S generics.R:107 | S dbarts.R:1594, 1599 [f16] |
| hazard | R RIB:1626 [f6] | S MOD:3095 [f6] | S MOD:3131 | S generics.R:102 [f24] | S dbarts.R:1594, 1599 [f6] |
| hurdle | R bart.R:1123 | - [f12] | - [f12] | M generics.R:129 [f25] | - [f12] |
| bcf | S COM:844, 887-891 [f17] [f48] | S MOD:2870, CH:1410 [f26] | S CH:1650 [f18] | M generics.R:129 | R [f23] |
| grouped | S MOD:4642-4664 | S MOD:4799 [f27] | S MOD:4808 | S generics.R:1545 | S MOD:4830 [f27] |
| hetero | S CH:4085, MOD:306 | S CH:4080 [f27] | - [f18] | ? generics.R:57 [f28] | S test-calibration-prior-draws.R:253 [f29] |

## 4. Model composition

| model | variance forest | grouped ranef | DART | warm start | grow-from-root |
|---|---|---|---|---|---|
| gaussian | S FAC:737 | S CH:626 | S CH:570 | S bart.R:1206 | S dbarts.R:995 |
| student | R spec.R:423 [f30] | S CH:626 | S CH:570 | S bart.R:1206 | S dbarts.R:995 |
| probit | R spec.R:410 | S CH:626 | S CH:570 | S bart.R:1206 | S dbarts.R:995 |
| logistic | R spec.R:410 | S CH:626 | S CH:570 | S bart.R:1206 | S dbarts.R:995 |
| ordinal | R spec.R:410 | M RIB:2986 [f31] | S CH:570 | R bart.R:608 | R bart.R:608 |
| nbinom | R spec.R:410 | M RIB:2991 [f31] | S CH:570 | R bart.R:608 | R bart.R:608 |
| multinom | R bart.R:876 | M RIB:1923 [f32] | R bart.R:876 [f33] | R bart.R:608 | R bart.R:608 |
| aft | R spec.R:410 | S CH:626 | S CH:570 | S bart.R:1206 | S dbarts.R:995 |
| hazard | R spec.R:410 | M rbart.R:49 [f6] | S CH:570 | S bart.R:1206 | S dbarts.R:995 |
| hurdle | R spec.R:410 [f34] | M | S bart.R:2247, 2256 [f35] | R bart.R:608 | R bart.R:608 |
| bcf | R FAC:817 [f48] | R RIB:2352 | R spec.R:537, 567-573 | S SAM:733 [f36] | S CH:1926 [f36] |
| grouped | ? [f30] | - | S rbart.R:591 | M rbart.R:45-51 [f37] | M rbart.R:45-51 [f37] |
| hetero | - | ? [f30] | S CH:570 [f38] | S SAM:733 | S CH:1921 |

Grow-from-root is gated by the LEAF model, not the family: linear and GP leaves
are refused at dbarts.R:998-1007 and no-op at CH:1904, so every family above
reads "constant leaf" in that column.

## 5. Evidence

| model | equivalence baseline [f39] | SBC verdict [f40] | dedicated tinytest files |
|---|---|---|---|
| gaussian | 21 scenarios (friedman, weighted, splitprobs, chains, setdata, wtoffset, quants, categorical, missing, dart, linear, gp, zeroweights, sparse, wtgp, chik2, pred*) | PASS 7/7 | ~20 (test-sampler-*.R, test-bart-bart2.R, test-zero-weights.R) |
| student | `student` | PASS 4/4 | test-robust-errors.R only |
| probit | `probit`, `chik` | PASS | test-binaryResponse-hyperprior.R, test-family.R, test-weighted-binary-ppd.R |
| logistic | `logistic`, `wtlogistic` | PASS 6/6 | test-weighted-logistic.R, test-family.R |
| ordinal | `ordinal` | 9/10 [f41] | test-ordinal.R only |
| nbinom | `nbinom` | 1/3 [f42] | test-nbinom.R only |
| multinom | 10 scenarios, own harness | aggregate PASS, raw `f_ik` OPEN [f43] | 5 (test-multinomial-*.R) |
| aft | `grouped_aft` only [f44] | OUT [f45] | test-aft.R, test-rbart-aft.R |
| hazard | `hazard` | OUT [f45] | test-hazard.R only |
| hurdle | `hurdle` | OUT [f45] | test-hurdle.R, test-hurdle-surface.R |
| bcf | 12 scenarios, own harness, gaussian only [f48] | PASS, gaussian only [f46] | 7 (test-bcf*.R) |
| grouped | `grouped`, `grouped_aft` | PASS (tier A) | 14 (test-rbart-*.R) |
| hetero | `hetforce`, `hetswap`, `hetpartial` | OUT [f47] | test-heteroscedastic.R, -mutation.R |

Flat-C test coverage is thinner than family reach: inst/tinytest/test-capi.R
drives only the `""`/`"probit"` tokens plus grouped (:672), heteroscedastic
(:223) and BCF (:1335) by control attribute. logistic, aft, ordinal and nbinom
are reachable through dbarts.h and untested there.

## Footnotes

[f1] Since S10 (`bart2-argument-consolidation`, 726dab10) `bart()` carries a
narrow, appended `family` formal, `c("auto", "logistic", "aft")` (`bart.R`,
`bart`'s formals), forwarded to `dbarts()` verbatim; `"auto"` (the default)
still infers gaussian/probit from the response, and `resid.dist` remains the
separate Student-t lever, untouched by `family`. Four bart2-own-class tokens -
`"multinomial"`, `"ordinal"`, `"nbinom"`, `"hurdle.lognormal"` - are refused BY
NAME ahead of `match.arg`'s generic message, by the shared
`refuseBartOwnClassFamily`/`bartOwnClassFamilies` (`bart.R`); ordinal's own
auto-detection off an ordered-factor response routes through the same helper
now, no longer a hand-rolled backstop. The remaining six tokens - `"gaussian"`,
`"probit"`, `"hazard.probit"`, `"hazard"`, `"hazard.logistic"`, `"twopart"` -
are simply not in the formal's vocabulary at all and fall to `match.arg`'s
generic "should be one of" message if typed: the first three add no capability
`"auto"` does not already give, and `"hazard"`/`"hazard.logistic"` need
`breaks`/`max.rows`, which `bart()` does not have (both reasons stated in
man/bart.Rd's `family` item). `bcf`, `grouped` and `hetero` are not `family`
tokens at all and stay out of reach by signature, unaffected by S10.

[f2] Student-t is not a token anywhere. It is selected by a finite `resid.df`
attribute on the model SEXP (RIB:2558-2569, gaussian-only gate) and refused for
every non-gaussian family R-side at spec.R:325-331. The engine family stays
`gaussian`, which is why the whole gaussian row applies to it in tables 2-4.

[f3] Reachable through `dbarts_sampler_create` but not discoverable from the
shipped header: ordinal needs the `bartcore.n.categories` control attribute
(RIB:456), nbinom `bartcore.dispersion` (RIB:469), grouped `bartcore.groups`
(RIB:1926), heteroscedastic `bartcore.variance` (RIB:2016). The header's
`family` documentation (CAPI:643-652, on `dbarts_sampler_create`; NOT
CAPI:443-462, which is the `DBARTS_C_API_LIST` X-macro body and never carried
it) names only probit, logistic, gaussian and aft, and the K-forest paragraph
beside it (CAPI:668-690) now names gaussian, probit and logistic, M4.4 having
replaced its "Gaussian responses only" at CAPI:686-690.

[f4] Multinomial has no dbarts.h creation path at all; it is reached only by
the R-internal `.Call` entries `C_dbarts_bartcore_createMultinomial` /
`...Counts` (bartcore.R:846, 909).

[f5] SUPERSEDED by S10 (`bart2-argument-consolidation`, 726dab10) and S14
(9031b348): `family = "aft"` is now an explicit, appended token on `bart()`
(`c("auto", "logistic", "aft")`), documented in man/bart.Rd's `family` item,
rather than an undocumented `Surv()` auto-dispatch quirk - see [f1]. The
underlying `Surv()`/two-column-`y.train` detection this footnote used to flag
is unchanged.

[f6] `family = "hazard"` / `"hazard.probit"` / `"hazard.logistic"` is
person-period ingestion sugar: dbarts.R:487-536 expands the design and remaps
the token to `"probit"` or `"logistic"` at dbarts.R:533 before any model is
built. The resulting sampler *is* an ordinary binary one, so its whole row
equals the probit (or logistic) row, and the fit records `family = "probit"`.
No engine code, hence no C-API token and no SBC arm. Refusals inside the
expander: no formula interface (:488), no `subset` (:503), no `test` (:506).

[f7] `treatment` is still not a `bart2()` formal - it is refused as an unknown
dots argument by the shared `rejectUnknownDotsArgs` (bart.R, called on bart2's
own dots). But since S12 (`bart2-argument-consolidation`, 4b179585) the general
K-forest amplitude capability this row names IS reachable from `bart2()`'s
FORMULA interface: a `forest()` term - `z:forest(x1 + x2)`, or the general
`forest(x1 + x2, basis = ~z)` form - rewrites the formula and feeds the same
`forests = list(forest(), forest(basis = ))` channel `dbarts()`/`dbartsSpec()`
already used (`ingestFormulaTerms`, `R/formulaTerms.R`, called from `dbarts()`'s
formula path; `bart2` dispatches through it for every family except
multinomial, which refuses a `forest()` term by name since it has no
amplitude-coupled slot). A two-forest fit built this way carries the same
`bartcore.bcf` control attribute as one built through `forests =`
(inst/tinytest/test-formula-terms.R Block B, byte-identical at the same seed).
What remains absent is the NAMED causal verb, not the general capability:
`bcf()` and `bartBCF` still ship in **bartCause**, not dbarts
(bcf-public-surface.md:483-496, fork 4 RESOLVED VD 2026-08-11; echoed
bcf.md:384-409). `forests =` on `bart2` itself was considered and declined as a
second spelling (bart2-argument-consolidation 5.5.3) - the formula route is the
only one bart2 gets; forest 1's `sd`/`update.amplitude` have no term spelling
(same plan, 5.7), which is why `bart2()`'s own knobs still don't reach a term's
forest the way `forests =` on `dbarts()` would.

[f8] `bartcore.groups` is written at exactly one site, rbart.R:388, and no other
entry point carries a `group.by` formal, so grouped random effects are an
`rbart_vi()`-only surface.

[f9] `updateScale` re-derives the internal response transform. The latent
families have `fitScale() == 1` and `fitShift() == 0` by definition, so there
is no transform to re-anchor and the flag is ignored rather than refused.

[f10] Logistic weights are the observation counts its Polya-Gamma latents were
built from, so they cannot change after creation (RIB:2756). This is the one
weight refusal that is "unbuilt" rather than "incoherent": creation accepts
positive-integer weights at RIB:1630.

[f11] `bart2(family = "multinomial")` returns a `dbartsSampler` that is only a
HOST SHELL - it carries the design and priors but not the model, which lives on
the separate bartcore handle `$bc` (bart.R:1451, 1534). Every R5 mutator is
refused by `refuseHostMutation` (dbarts.R:914-925). At the handle level:
`setResponse` R RIB:2667 and `setOffset` R RIB:2664 (both redirect to the counts
channel), `setWeights` R RIB:2674, `setSigma` R RIB:2873, `setTestOffset` R
RIB:2419 (a flat offset leaves the simplex), `setPredictor` S RIB:4989,
`setTestPredictor` S RIB:4666. Its own channels - `setCounts` RIB:3741,
`setCategoryOffset` RIB:3819, `setCategoryTestOffset` RIB:3858 - are S at the
bridge but unexported (`dbarts:::`) and absent from the R5 object.

[f12] Hurdle has no sampler of its own: `dbarts()` refuses construction at
dbarts.R:445 and `bart2Hurdle` (bart.R:2206) composes two ordinary `bart2()`
fits - an occupancy probit (bart.R:2251) and a lognormal positive part
(bart.R:2261) - glued at report time. The channel questions resolve on the
probit and gaussian rows of the two components.

[f13] Since adoption-slate S3 (eeedc07c) a grouped sampler accepts a
same-length setResponse and setOffset at the pinned scale - faithful
delegation, MOD:4762 - and refuseGroupedScaleUpdate (RIB:2727) refuses
updateScale != FALSE only under a base family with a data-derived transform
(gaussian, which is Student-t's report, and aft): b and tau are held on the
base's internal scale and converted by nothing, so a re-anchoring swap would
silently restate both in response units. Grouped probit and logistic take
updateScale = TRUE as the no-op it always was. The flat C API guards through
the same call (C_interface.cpp:585, 607). setData stays refused.

[f14] Reads off the BASE family: grouped gaussian takes `setWeights` (MOD:4788)
and `setSigma`; grouped probit is refused on both (RIB:2763, RIB:2873); grouped
aft takes `setSigma` and refuses `setWeights`.

[f15] Arc `latent-subset-mask` (TODO:469), design FINAL, ARC COMPLETE (S0
through S4 LANDED); artifacts .claude/latent-subset-mask-design/. A first-class
0/1 `setActiveRows` channel each family composes into its own precision vector,
with the latent draw skipped for inactive rows. Slices: **S0** pins (no engine
change); **S1** the channel plus gaussian, Student-t, probit, ordinal; **S2**
logistic, nbinom, aft; **S3** multinomial (global only); **S4** surface,
records, baselines. S0 landed at dc11a805 (the pins, now
inst/tinytest/test-active-rows-pins.R). S1 landed at 6db22aee: the engine
channel - `Chain::setActiveRows` CH:1602, which owns the single validating and
normalizing scan, `Sampler` SAM:1242, the facade's pure virtual FAC:323 and its
shape probe FAC:102 - plus gaussian, Student-t, probit and ordinal, the R5
`$setActiveRows` (dbarts.R:1274) and the bridge entry (RIB:3988). S2 landed at
87d370ea: logistic (`workingWeights()` MOD:3531) and nbinom
(`workingWeights()` MOD:4332) serve a SEPARATE a_i omega_i composite rather
than writing the zero into omega_ itself, since the working response divides
by it and 0 * inf in the node kernels is a NaN; nbinom's `setActiveRows`
(MOD:4359) additionally restricts the collapsed statistic S the dispersion
grid draw reads and REBUILDS the count-histogram kernel behind L_k
(`NBDispersionPrior::computeKernel` MOD:4254) over the active rows at every
mask change - the channel's one per-install cost; aft's `setActiveRows`
(MOD:3825) composes into its contained Gaussian, inheriting the sigma
degrees-of-freedom recount, and skips the censored redraw at an inactive row
(MOD:3809). All three report NaN pointwise log-likelihood at an inactive row.
Oracles: per-family kernel comparisons against the compacted arm, bitwise in
value and in RNG stream (`testActiveRowsLogisticKernel`
tests/cpp/test_model.cpp:5443, `testActiveRowsNBKernels` :5529,
`testActiveRowsAFTCensored` :5621 - each latent being a rejection sampler
means a discard-rather-than-skip at an inactive row fails the arm outright),
plus a sampler-level conditional independence oracle under substituted
inactive responses (inst/tinytest/test-active-rows-pins.R's S2 block:
substituting arbitrary in-support values at the inactive rows leaves every
active row's recorded draw bitwise). The heteroscedastic and grouped pins
move from FINITENESS to bitwise `setWeights(w * a)`
(test-active-rows-pins.R's hetero block; tests/cpp/test_sampler.cpp's grouped
block in `testActiveRows`), and ordinal - already S1 - gains a sampler-level
independence arm of its own beside the kernel-level coverage. S3 landed at
8b047f8b: multinomial's mask is GLOBAL and lands on the softmax coupling
rather than the response, which holds no precisions of its own -
`MultinomialResponse::setActiveRows` (MOD:3726) is a pass-through that only
advertises the capability (MOD:3725), and `Chain::setActiveRows` forwards the
mask to `MultinomialForestCombiner::setActiveRows` (COM:1678) after the
response's own install (`ForestCombiner::setActiveRows` COM:721 is the inert
default every additive coupling relies on instead). An inactive row's K
interleaved Polya-Gamma draws are SKIPPED, not drawn and discarded, in
`drawForestGlue` (COM:1749), and its composed precision is zeroed in every
category in `formForestResponse` (COM:1808-1809); the row keeps its leaf
occupancy and its reported softmax probabilities, and omega is never zeroed
since the working response divides by it. PER-FOREST masking is refused
permanently on model grounds at the only reachable per-forest,
per-observation channel, `bartcore_setForestWeights` (RIB:3948-3952) - see
[f21] for the full statement. The bridge's active-row refusal (RIB:3997) no
longer names multinomial: the old per-family `activeRowsFamilyName` helper is
deleted, and the message is now family-generic, reached only by a future
family that does not override the base refusal. Oracles: the kernel-level
`testActiveRowsMultinomialKernel` (tests/cpp/test_sampler.cpp:5557) pins the
SKIP semantics and the zeroed composed precision bitwise against a compacted
combiner, in both value and Polya-Gamma stream; the shape probe flips in
`testMultinomial` (tests/cpp/test_shape.cpp:379); and
inst/tinytest/test-active-rows-pins.R's S3 block adds the same
substituted-response independence arm the S2 families got (moving successes
AND trial counts at the inactive rows, since PG(n_i, .) sums n_i variates),
plus an all-zeros mask run (every category forest at its prior, every row
still reporting a simplex) and the `setForestWeights` model-grounds refusal
(also pinned in test-forest-weights.R). S4 landed at 93afd635
(implemented as 76fd3ba6, amended during independent review): Rd
(man/dbartsSampler-class.Rd), NEWS (inst/NEWS.Rd), a named recipe
(man/bart.Rd), the dbarts.h reservation (docs/plans/c-api-growth.md),
two new equivalence.R scenarios (maskprobit, maskordinal) and one
bcf-equivalence.R scenario (masked, pinning BCF - see [f26]). The
flat-C entry, `dbarts_sampler_setActiveRows`, LANDED at dbarts-h-reshape
S1 (ab3aa2fa, 2026-08-13; body `C_interface.cpp:1053-1062`): the
capability probe on `shape.supportsActiveRows` runs first and never
switches on family, so a probit sampler is reachable from C too - an
amend over this section's own proposed probit refusal, since every
`ResponseModel` subclass now reports `supportsActiveRows`
(`inst/tinytest/test-capi.R:1265-1305` pins a genuine mask, an all-ones
no-op, a fractional refusal, a NULL clear, and a probit mask moving
draws). ARC FULLY COMPLETE, R and flat C alike.

[f16] Arc `nameable-calibration` (TODO:611), design AMENDED FINAL, ARC
COMPLETE; artifacts .claude/nameable-calibration-design/. Names the
per-forest prior ANCHOR (`prior.scale`, the forest-total prior scale at k = 1,
in response units) rather than an sd, with a `$getCalibration` /
`$setCalibration` pair. Slices: **S0** signature freeze, LANDED 4c866286;
**S1** creation half, LANDED c2a7e89b; **S2** mid-chain get/set, LANDED
d809b944 (+ a records correction, 7da36dc3); **S3** the flat-C half, LANDED
inside dbarts-h-reshape S1 (ab3aa2fa). The R surface was already COMPLETE:
`$getCalibration`/`$setCalibration` read and write every chain of any
single-forest sampler, with a 1-based `forest` arg (`resolveForestIndex`,
bartcore.R:1128) mapped onto the engine's 0-based one. S1 names the model's
`prior.scale` slot (A_class.R:398), resolved from `node.prior`'s `scale =` /
`sd =` spelling at `dbartsSpec()` (spec.R:316, `resolvePriorScale` in
R/model.R), and converted against the response transform by a private
`Chain::resolvedNodeScale` helper (CH:3855) shared by the single-forest
constructor and every `setModel` reinstall, so a round trip through the model
SEXP no longer reverts a named calibration. S2 adds the reader,
`Chain::forestCalibration` (CH:1165) - the AUTHORITATIVE report of what is in
force, independent of the model's recorded intent, so a `setResponse` /
`setOffset` at `updateScale = TRUE` or a `setData` shows up as a move rather
than staying silent - and the writer, `Chain::setForestPriorScale` (CH:1207),
sharing one `priorScaleFactor` conversion (CH:3865) with the reader so neither
direction can drift from the other; both are total over the four leaf models
and carry no family switch (facade FAC:309, 316; `Sampler` SAM:1220, 1229; R5
`dbartsSampler$getCalibration`/`$setCalibration` dbarts.R:1594, 1599; bridge
`bartcore_getCalibration`/`bartcore_setCalibration` RIB:4105, 4160). Refused
under a `k` hyperprior (the `sd` spelling only, since a sampled `k` has no
single value to divide by, or once the chains' `k` have diverged) and for
BCF/multinomial forests at creation and again mid-chain (see [f23]);
`prior.mean` is refused as not writable, naming the `setOffset` recipe. `NaN`
is refused as a malformed value rather than read as the unnamed spelling, both
at creation (R/model.R:1455) and mid-chain (`validateLiveScale`,
R/model.R:1464). Shipped tests: inst/tinytest/test-calibration-creation.R (two
composed probit arms at construction ranges 16x apart agree to 1e-12 under a
shared name, against 8.6 and 2.5 unnamed), inst/tinytest/test-calibration-prior-draws.R
(what the named quantity means per leaf model - exact for the constant leaf,
bounded inequalities for linear/gp/monotone - across ten family and decoration
paths, [f29] adding the tenth), inst/tinytest/test-calibration-midchain.R (a
read-then-write is BITWISE inert on the internal scale, not merely on the
report, since the setter skips a write reproducing either spelling of what is
in force; a write-then-read is pinned at ULP level rather than bitwise, since
`(P / f) * f != P` for about a tenth of positive pairs - 7da36dc3 corrected
which of the fixture's four `(m, P)` pairs actually rounds and added the
rounding cell, `m = 25` at `P = 0.25`, so the tolerance is exercised rather
than merely permitted; plus the refusal matrix, every mutation channel the
reported value must not surprise on, and the save/load round trip), and a
component test at the engine boundary (`testForestCalibration`,
tests/cpp/test_sampler.cpp:6069). **S3**, the flat C entries
(`dbarts_sampler_forestCalibration`/`setForestPriorScale`,
`inst/include/dbarts/dbarts.h`), LANDED at dbarts-h-reshape S1 (ab3aa2fa,
2026-08-13) - the R-user-facing capability above was already `S` before S3,
since the flat-C gap never gated a `dbartsSampler` caller; S3 closes the
flat-C gap itself, with the two carried items from S2 shipping alongside:
the engine bounds check on `Chain::forestCalibration` (`chain.hpp:1170`,
returning a default-constructed calibration rather than reading past the
last forest) and the `refuseBCFMutation` reorder in `$setCalibration`
(`dbarts.R:1607-1627`, argument validation before the BCF refusal).
EXTENDED by `binary-kforest-prior-default` S1: the reader reports five further
columns - `amplitude.prior.variance` and `amplitude.prior.scale` (exclusive
per forest), `node.scale.factor`, `node.scale.divisor` and `basis.row.norm` -
the multi-forest calibration map's own decomposition of `prior.scale`, NaN on
every forest whose scale that map does not own, with the matching fields
appended to `dbarts_forest_calibration` below its 1.0-0 boundary (`sizeof`
moves, and so, since the token folds every ABI struct's layout, does the
apiHash - it did not when this landed). They are TRUE after a state install
rather than
a spec echo: a donor leaf scale differing bitwise from the one in force sends
both `node.scale` columns to NaN until `$setForestBasis` re-imposes the map,
while the amplitude prior follows the installed state. Carried by
test-calibration-midchain.R (the exact 12-column shape, the off-map NaN on the
single-forest and multinomial routes, and the anchor recovered by the
identity), test-bcf-family.R (the factors themselves under both latent
families, where only the logistic arm sees a dropped anchor),
test-forest-basis-r5.R (the row norm following a swap, and the four
truthfulness arms), test-capi.R (the appended fields, and a PRE-APPEND caller
whose five buffers are left untouched), and `tests/cpp`
`testBCFCalibrationMap`.

[f17] Zero weights are accepted, not refused (A_class.R:578-582 errors only
below zero and warns that zeros are ignored; bridge RIB:4823). The conditionals
are exact - leaf suffstats multiply by `w` (MOD:314, 1178), and the sigma
posterior counts only positive-weight rows (`numPositiveWeights_` MOD:2785,
recounted on every install at MOD:2974, consumed MOD:2799-2803). The one named
inexactness against a true subset fit is CLOSED (`empty-leaf-veto-fix`,
2026-08-12): the empty-leaf veto counts
POSITIVE-WEIGHT members, so a leaf held alive only by zeroed rows is empty and
its branch is vetoed, on the conjugate path (MOV) and the constrained-leaf path
(MOD) alike. Occupancy elsewhere - the birth scan's `count`,
`collapseEmptyNodes`' trigger, `stateIsValid` - still counts members
deliberately, so this does NOT make zero-weight occupancy match a compacted fit;
see docs/design/empty-leaf-veto.md, "What counts as empty". The same fix covers
the Student-t row and a GAUSSIAN K-forest. It says nothing about a probit or
logistic one, where a zero weight cannot exist to begin with: probit refuses
weights entirely and logistic holds them to positive integer counts, both at
creation, so the cell is family-dependent ([f48]).

[f18] For gaussian and heteroscedastic no latent vector exists: both leave
`ResponseModel::latents()` at its nullptr default (MOD:2692), and the bridge
returns `R_NilValue` (RIB:5886). A K-forest sampler is no longer one of them.
`Chain::latents()` (CH:1650) is a bare delegation to `response_->latents()`
carrying no coupling gate and no family switch, and `bartcore_getLatents`
(RIB:5884) gates only on that pointer being null, so since M4.4 a probit
K-forest reports its truncated normals (`ProbitResponse::latents` MOD:3131) and
a logistic one its Polya-Gamma omegas (`LogisticResponse::latents` MOD:3605). A
GAUSSIAN K-forest still reports none, which is why this cell is
family-dependent rather than plainly `S`.

[f19] A Student-t fit records `family = "gaussian"`, so
`extract(type = "loglik")` takes the gaussian branch and evaluates `dnorm`
against the scalar `object$sigma` (generics.R:57-62). The t marginal is not
distinguished and the per-row `lambda_i` is not stored on the fit. No guard, no
doc, no TODO covers this - flagged, unadjudicated.

[f20] Weights on logistic are PG copy counts and a zero count is refused by
name at creation ("drop zero-count rows", RIB:1630-1635; R mirror
spec.R:61-69), so zero-weight subsetting is foreclosed for this family by the
weight semantics themselves - it is exactly the hole `latent-subset-mask` S2
filled (87d370ea), by the mid-chain `setActiveRows` channel rather than by any
change to the zero-count creation refusal, which stands.

[f21] The GLOBAL channel shipped at 8b047f8b, landing on the softmax coupling
rather than the response, which holds no precisions of its own to compose a
mask into: `MultinomialResponse::setActiveRows` (MOD:3726) is a pass-through
that only advertises the capability (`supportsActiveRows` MOD:3725), and
`Chain::setActiveRows` forwards the mask to
`MultinomialForestCombiner::setActiveRows` (COM:1678) after the response's own
install. An inactive row's K interleaved Polya-Gamma draws are SKIPPED rather
than drawn and discarded, in `drawForestGlue` (COM:1749), and its composed
precision is zeroed in every category in `formForestResponse`
(COM:1808-1809); the row keeps its leaf occupancy and its reported softmax
probabilities, and omega is never zeroed since the working response divides
by it. PER-FOREST masking stays REFUSED, permanently and on model grounds:
the softmax margin is a log-sum-exp over the other K-1 forests, so a row
absent from category k's forest is still in every other category's
likelihood, and "row i is out of category k only" restricts no likelihood at
all. The refusal lands at the only reachable per-forest, per-observation
channel, `bartcore_setForestWeights` (RIB:3948-3952), naming the model reason
rather than "unbuilt". BCF's per-forest weight acceptance at that same
channel stands unaffected - a different (additive) coupling where the
per-forest mask is redundant with, not incoherent under, the combined
likelihood (see [f26]).

[f22] Multinomial's omegas live in the combiner, not the response model, and
`MultinomialResponse` does not override `latents()` (MOD:3690-3767), so
`getLatents` returns NULL. No accessor exposes them.

[f23] A named `prior.scale` is refused for BCF and multinomial forests both AT
CREATION and MID-CHAIN, by design - their per-forest leaf scales come from a
calibration map that owns them (map at COM:275-300), so a named value has
nowhere to land. Three creation-time refusal sites shipped at c2a7e89b: R-side
`dbartsSpec()`'s BCF composition (spec.R:555, the `"a named 'prior.scale'"`
entry of the `unsupported` vector; the earlier `spec.R:440` named the
non-default-`k` entry instead), the engine's own
BCF-composition gate (`refuseUnsupportedBCFComposition`, RIB:2346), and the
multinomial forest builder (`buildMultinomialSampler`, RIB:3255). S2
(d809b944) adds the mid-chain refusals, at TWO independent sites rather than
one shared gate: `$setCalibration`'s R5 method refuses BCF through
`refuseBCFMutation` (dbarts.R:1622, MEASURED "two-forest calibration map",
test-calibration-midchain.R:399-402) before ever reaching the bridge, and
refuses a multinomial fit's host shell through `refuseHostMutation`
(dbarts.R:1607, MEASURED "host sampler of a bart2", line 460-463); underneath
both, the engine-level gate any DIRECT low-level call still hits -
`Chain::setForestPriorScale` returning false whenever `combiner_ != nullptr`
(CH:1207), surfaced as `Rf_error(...calibrationMapName...)` at the bridge
(RIB:4172) - is what the unexported `dbarts:::bartcoreSetForestPriorScale`
hits on a multinomial forest's low-level handle (MEASURED "softmax calibration
map", line 441-444); the R5 layer never routes a BCF sampler there since
`refuseBCFMutation` refuses first, so only the multinomial arm exercises the
bridge gate directly. These cells stay `R`.

[f24] Evaluated per PERSON-PERIOD row, not per subject, since the fit's
response is the expanded binary indicator.

[f25] The composed hurdle fit has `family = "hurdle.lognormal"` and errors at
generics.R:129-134, but each component fit (`$occupancy` probit, `$positive`
gaussian) supports `extract(type = "loglik")` on its own.

[f26] SHIPPED, pinned bitwise at mask S4 (93afd635). Nothing on the path gates
on the coupling: the shape probe (FAC:102) reports `supportsActiveRows`, the
mask composes into whatever precision the installed response owns, and then
into the per-forest weights at `composeForestWeights` (CH:1410). M4.4
falsified the reason this footnote used to give - "a BCF sampler's response IS
a `GaussianResponse`" - without touching the conclusion: the K-forest chain now
builds a `ProbitResponse` or a `LogisticResponse` too (CH:755), and each
overrides `setActiveRows` on its own terms (MOD:3095, MOD:3564) exactly as it
does off a coupling, so [f15]'s S1 and S2 arms carry the latent K-forest with
no edit of their own. Gaussian's composition into the case weights, which is
what inherits the sigma df, is `GaussianResponse::setActiveRows` MOD:2870. The
measurements below are all GAUSSIAN two-forest ones; no latent K-forest mask is
measured anywhere. MEASURED at 6db22aee on a 200-row two-forest
sampler: `$setActiveRows(a)` and the bridge `bartcoreSetActiveRows` are both
accepted; on a sampler carrying `w` the mask is BITWISE `setWeights(w * a)` in
`train` and in `sigma`; an all-zeros mask runs finite; a fractional element is
refused. PINNED at mask S4: `inst/tinytest/test-active-rows-pins.R:84-112`
(masked-bcf, bitwise vs `setWeights(w * a)` on train and sigma) and the
`bcf-equivalence.R` `masked` scenario, carried by the current
`bcf-equivalence-6e3b9fb8.rds`. A per-forest mask is refused as REDUNDANT rather
than unbuilt: `setForestWeights` (RIB:3937) already expresses it - though note
that channel is deliberately NOT row removal (CH:1120-1150: it does not remove
the row from occupancy, the combination or the sigma df; it DOES reach that
forest's empty-leaf veto, which counts positive composed weights). It is now a
PUBLIC R5 method, `dbartsSampler$setForestWeights` (dbarts.R:1297, landed
multiforest-extension-surface M1, 05ac3b4b), 1-based via `resolveForestIndex`
(a BCF basis forest is `2L`) and mirrored across re-creation on a
`forestWeights` field re-applied at `getPointer`, `setState` and `copy`
(dbarts.R:1665); the unexported `bartcoreSetForestWeights` (bartcore.R:1056)
stays the 0-based internal wrapper the R5 method does not call.

[f27] Delegating / decorating, and that is what 6db22aee landed for the
active-rows column: neither row took an edit of its own. `GroupedResponse`
forwards `setActiveRows` to its base (MOD:4799) exactly as it forwards
`setWeights` (MOD:4788), advertising the base's capability (MOD:4796), and
`drawGroupEffects` already weights its per-group sums by `workingWeights()`
(MOD:4724), so an inactive row leaves its group's mean and precision and an
all-inactive group falls back to its prior through the same formula. The
heteroscedastic `formMeanWeights` (CH:4080-4087) reads
`response_->workingWeights()` at CH:4083 - the COMPOSED `w * a` while a mask is
installed - and divides by `s^2(x_i)`, so a zero stays a zero. Both are
PINNED: grouped at tests/cpp/test_sampler.cpp:1703 (an entirely inactive group
draws its effect from the prior, finite), heteroscedastic at
inst/tinytest/test-active-rows-pins.R:244-269. Both pins were STRENGTHENED from
FINITENESS-only to bitwise masked-vs-`setWeights(w * a)` at mask S2
(87d370ea): grouped's effects, training fits and sigma all agree bitwise
against a composed-weight sampler (test_sampler.cpp:1705-1708); heteroscedastic
likewise for train and varcount (test-active-rows-pins.R:267-268). The same
delegation carries BOTH halves of the
nameable-calibration column for grouped: `GroupedResponse::fitScale()`/`fitShift()`
forward to their base (MOD:4830), so `Chain::resolvedNodeScale` (CH:3855) at
creation and `Chain::forestCalibration`/`setForestPriorScale` (CH:1165, 1207)
mid-chain all convert a named `prior.scale` exactly as they do for the
undecorated family, with no edit of grouped's own. Creation-time, grouped is
one of the nine family/decoration paths c2a7e89b's own test measures
(inst/tinytest/test-calibration-prior-draws.R:223, MEASURED 0.74210 vs a 0.75
target); the mid-chain half rides the same generic mechanism every
non-coupled family does, with no dedicated grouped arm of its own in
test-calibration-midchain.R.

[f28] A heteroscedastic fit also records `family = "gaussian"` and takes the
same gaussian loglik branch with the SCALAR `object$sigma`, ignoring the
per-observation `s.train` surface the fit stores separately (bart.R:243-260).
Same standing as [f19]: no guard, no doc, unadjudicated.

[f29] The variance forest is a separate leaf model entirely, outside
`forests_`, and is not addressable by the mid-chain `setCalibration`: a
SHIPPED door (nameable-calibration synthesis 2.6 item 7), not an open one -
`Chain::setForestPriorScale`'s bounds check `f >= forests_.size()` (CH:1207)
never sees the variance forest, so a heteroscedastic sampler's
`shape.numForests` is 1 and `forest = 2` is refused by the ORDINARY
out-of-range check every single-forest sampler hits (`bartcore_getCalibration`
RIB:4108, `bartcore_setCalibration` RIB:4164), not by a hetero-specific gate.
The MEAN forest's own calibration - both halves - is not gated by that door:
`Chain::resolvedNodeScale` (CH:3855) runs at forest.leaf.scale assignment
(CH:633-635) before the variance-forest branch (CH:724), and the mid-chain
reader/writer (`forestCalibration`/`setForestPriorScale`, CH:1165, 1207) read
the same `response_->fitScale()`/`fitShift()` - none of the three reads any
family flag. Creation-time it is now PINNED rather than merely constructing:
7da36dc3 added a heteroscedastic arm, the tenth of
`test-calibration-prior-draws.R`'s `anchorSamplers` sweep
(inst/tinytest/test-calibration-prior-draws.R:253-257), measuring a named
`scale = 1.5` against the shared 0.75 target inside the same 9% band
(`familyBand`, line 180) every other family in that loop is held to (line
259-265) - closing the earlier gap where the conversion ran unguarded but no
test exercised it. Mid-chain, the mean forest rides the same generic S2
mechanism as every other single-forest family, with no dedicated
heteroscedastic arm of its own beyond that shared mechanism and the
`forest = 2` refusal above.

[f30] Split verdicts, measured by benchmarks/R/composition-matrix.R.
`resid.dist = student()` + `variance =` now REFUSES at `spec.R:423` ("a
variance forest does not support Student-t residuals: the two are not yet
shown to compose") - a validation error only, the formal stays, and the door
memo records that adjudicating the composition of two scale mixtures on the
same precision channel reopens it. grouped + `variance =` and hetero +
grouped still CONSTRUCT (CH:626 decorates before CH:724 builds the variance
forest); whether those compositions are models anyone wants is not
adjudicated anywhere.

[f31] Recorded but UNBUILT doors, refused with that reason in the comment:
grouped ordinal because the cutpoint block and the group block are not yet shown
to interleave (RIB:2982-2987, ordinal.md section 8), grouped nbinom the same for
the dispersion block (RIB:2988-2992, negative-binomial.md section 7).

[f32] No surface at all: `applyGroupAttribute` (RIB:1923) is called from exactly
one site, RIB:2980, on the single-forest holder path, so `bartcore.groups` is
never read for a multinomial sampler.

[f33] Formerly a defect: `buildMultinomialForest` hard-sets
`forest.useDart = false` (CH:5076), and `buildMultinomialSampler`
(RIB:3238-3294) copies only power/base/proposal-probability fields, so a DART
tree prior built from either the `dart` argument or a `tree.prior` object
never reached the K-forest engine. `bart2` now refuses both routes by name
(bart.R:876), matching BCF's own named refusal (spec.R:537,
`buildBCFForest` CH:5027) before either reaches the host sampler.

[f34] `bart2Hurdle` builds both component calls with `redirectCall`
(bart.R:2247, 2256), so a user's `variance =` is forwarded to BOTH - including
the occupancy component, which then sets `family = "probit"` and hits the
non-gaussian variance refusal at spec.R:410 before either component fits.
That refusal is deliberate, not a bug: `hurdleSigmaVec`'s comment
(generics.R:1032-1042, definition :1043) states the positive fit is always
homoscedastic because the gate makes a heteroscedastic component
unreachable. The Rd side was the wrong one and is now corrected
(bart2.Rd:224, dbarts.Rd:108) to match.

[f35] `dart` is forwarded to both components (bart.R:2247, 2256), each of
which is an ordinary single-forest chain that takes it.

[f36] No family gate: `installForests` checks shape, grid, DART and
variance-forest presence only (SAM:701-751), matching donor forest counts at
SAM:733, and `growForestFromRoot` loops every forest (CH:1926, the loop; the
earlier CH:1589 named the variance-forest pre-step above it). Neither is
exercised by a BCF test, and BCF has no `bart2()` surface, so both are reached
only through the R5 `$installTrees` (dbarts.R:1738) / `$growFromRoot`
(dbarts.R:995).

[f37] `rbart_vi()` carries no `warm.start` or `n.grow.sweeps` formal
(rbart.R:45-51) and its unknown-argument check (rbart.R:59-62) rejects them. The
underlying R5 sampler carries no group gate on either path, so this is a surface
gap, not an engine one.

[f38] The MEAN forest keeps DART; the variance forest never takes it
(`buildVarianceForest` CH:4044 never sets `useDart`, default false at CH:125).

[f39] Current baselines: `equivalence-d15a2bfb.rds` (42 scenarios),
`bcf-equivalence-6e3b9fb8.rds` (12), `multinomial-equivalence-4d9a3337.rds` (11)
- benchmarks/baselines/MANIFEST:30, 61, 68. Scenario names are the keys in
`makeScenarios()`, benchmarks/R/equivalence.R:60.

[f40] docs/plans/sbc-family-tiers.md (status BUILT) plus
docs/plans/sbc-calibration.md (DONE). The A/B/C "tiers" in the latter are
FEATURE tiers (A baselines/DART/grouped/weighted/BCF, B linear leaf, C GP leaf),
not family tiers - there is no per-family tier ladder, only per-family recorded
verdicts, which is what this column carries.

[f41] gamma3 flagged in one stream and RESOLVED as the cutpoint-vs-mean-level
ridge mixing slowly (does not reproduce across streams; 0.31 of the band at 3x
the chain length). The cutpoint block, the latent eta and all K category
probabilities calibrate.

[f42] `avg.mu` - the identified mean - passes cleanly; `r` and `agg.psi` flag at
thin = 30 and cross into the band at 5x the spacing. Read as the r-vs-psi ridge
mixing slowly (H-MIX), on two ladder points rather than three; the recorded
full-R third point is still owed.

[f43] Aggregate `p_k(x*)` passes at both chain lengths; the three raw per-forest
`f_ik` cells carry a persistent U the ladder does not shrink. Pre-registered
suspect: `MultinomialForestCombiner::afterCombine`'s level-centering draw. OPEN.

[f44] AFT is exercised only in combination, through the `grouped_aft` scenario
(equivalence.R:476). There is no standalone AFT equivalence scenario; the
separate exact oracle benchmarks/R/aft-exact.R is not a MANIFEST entry.

[f45] Out of the SBC matrix by scope, each for its own recorded reason
(sbc-family-tiers.md:43-51): aft because its censoring status is fixed at
creation, so a prior-draw replication cannot vary it (the enabler is a status
setter); hazard and hurdle because their person-period / two-part designs depend
on `y0`, which breaks exchangeability, and because neither owns any sampling
code.

[f46] Tier A PASS with the sigma channel resolved as H-MIX on the (a, mu) ridge
(sbc-calibration.md:650-660). Explicitly out of the family-tiers matrix, and
`runSbcBCF` errors at that plan's HEAD (repair tracked separately in
docs/plans/runsbcbcf-repair.md).

[f47] OUT but DEFERRED rather than blocked: prior draws never reach
`varianceForest_` today, and the arm is liftable R-side through `setState`
(sbc-family-tiers.md:50-51).

[f48] The K-forest coupling's family reach, landed at
multiforest-extension-surface M4.4 (625794fd). gaussian, probit and logistic
build; aft, ordinal and nbinom are refused at all three creation routes, each
naming what it is missing (spec.R:505-527, `refusedBCFFamilyReason` RIB:2281,
`createBCFSampler` FAC:822-825 - the last sitting directly beside the
variance-forest door FAC:817, which is unchanged and family-independent). The
calibration map's anchor is now family-keyed, `latentScaleAnchor` (CH:4910):
sd(y) under gaussian, 1 under probit, pi/sqrt(3) under logistic, and stated per
unit of basis row norm (`basisRowNorm` CH:4939). Cell by cell in section 2:
`setResponse`/`setOffset` are OPEN under every family, because
`Chain::setResponse` (CH:1630) now hands the response `combinedFits()` rather
than forest 0's bare totals, which is what let the gaussian conjunct come off
`Chain::supportsResponseMutation` (CH:1045); the combiner's own opt-in (COM:1059)
is unchanged. `setWeights` and `setSigma` are REFUSED for probit and logistic
and open only for gaussian, through the ORDINARY single-forest guards now that
the sampler answers `shape.family` for itself (`refuseBinaryWeightChange`
RIB:2752, `refusePinnedSigmaChange` RIB:2865); at creation the shared
`enforceBinaryWeightPolicy` (RIB:1616) refuses a probit weight outright and
holds a logistic one to positive integer counts, which is what makes the
zero-weight-subset cell family-dependent too. `updateScale = TRUE` stays
REFUSED under EVERY family - NOT the latent convention `- [f9]` M4.4's own plan
bullet predicted: `refuseBCFMutation` (bartcore.R:34) keys on the sampler
carrying bases, never on the family, and the bridge's
`refuseMultiForestResponseMutation` (RIB:2646) keys on `numForests >= 2`, so a
probit K-forest is refused too, though its transform is the identity and the
re-anchoring the refusal guards against cannot occur. PINNED, with the open
conduit and both refusals above, at inst/tinytest/test-bcf-family.R:406-422.
That file is the whole of the latent K-forest's evidence: no equivalence
scenario and no SBC arm reaches one.

The IN-SAMPLE per-forest OUTPUT channel (as opposed to construction, above)
LANDED at bart2-argument-consolidation S11 (00abf336): `packageBartResults`
packages `forestFits`/`glue` - response-scale raw per-forest totals and the
ragged multiplier channel - onto any bcf-shaped fit regardless of how it was
built, with `forest1..K` names and a `forest.labels` attribute, and
`extract(type = "forest", forest =, contribution = )` serves it (raw slice or
the on-demand `basis %*% glue` contribution) under both `combineChains`
conventions; refused by name on a fit without forest reporting or with
`sample = "test"`. S12/S13 add `inst/tinytest/test-formula-terms.R`'s own
identity checks (Block B: a term-built two-forest fit is byte-identical to the
equivalent `forests =` one at the same seed) as further bitwise evidence for
the gaussian/probit/logistic construction route, though it is not an
equivalence-harness scenario and does not change the counts in table 5.
Out-of-sample per-forest replay (`predict()`) remains a door - S11 did not open
it (bart2-argument-consolidation N8; `docs/plans/bcf-bartcause-relocation.md`
:1604-1605).

## Gaps

Every MISSING (`M`) and UNVERIFIED (`?`) cell above, as candidate work items,
grouped by family. Scheduling is VD's and the orchestrator's; nothing here
carries a schedule, and REFUSED cells are deliberately absent - they are part of
the models.

**gaussian.** None; every column is S or an intentional `-`. The one standing
inexactness (zero-weight rows survived the empty-leaf veto) closed at
`empty-leaf-veto-fix` ([f17]).

**student (Gaussian + Student-t residuals).** No `rbart_vi()` surface
(rbart.R:62); no `xbart()` surface (xbart.R:2-33). Pointwise loglik evaluates
the Gaussian density, not the t marginal, with no guard ([f19]). Composition
with a variance forest is refused by name ([f30]). Only one
dedicated tinytest file.

**probit.** None.

**logistic.** No `rbart_vi()` token (rbart.R:49), so grouped logistic is
engine-reachable but not R-reachable. No flat-C test coverage.

**ordinal.** No `xbart()` or `rbart_vi()` reach (both refused at data.R:500 as
unsupported response shapes). Pointwise log-likelihood unbuilt
(generics.R:129-134). Grouped ordinal is a recorded unbuilt door (RIB:2986).
`warm.start` and `n.grow.sweeps` unbuilt for the arc (bart.R:608). Its selecting
control attribute is undocumented in the shipped header ([f3]). One dedicated
tinytest file. SBC gamma3 resolved but not re-run at full R.

**nbinom.** As ordinal, plus: no `bart()` token, no `xbart()` token, no
`rbart_vi()` token. Grouped nbinom a recorded unbuilt door (RIB:2991). Pointwise
loglik unbuilt. Header attribute undocumented. One dedicated tinytest file. SBC
`r`/`agg.psi` flag standing (H-MIX read, third ladder point owed). Real-valued
dispersion remains a recorded door (TODO `negbin-real-dispersion`).

**multinomial.** No flat-C creation path at all ([f4]). No mutable
`dbartsSampler` - only a host shell ([f11]) - and its three real channels
(`setCounts`, `setCategoryOffset`, `setCategoryTestOffset`) are unexported with
no R5 methods. No `getLatents` ([f22]). No pointwise loglik. No grouped surface
([f32]). No warm start / grow-from-root. DART, `split.probs`, `monotone`, and
`variance` are all refused by name rather than silently dropped ([f33]).
SBC raw `f_ik` OPEN.

**aft.** No `xbart()` token. Pointwise loglik ships, but `setWeights` is refused
and the censoring status is fixed at creation, which is also what keeps AFT out
of the SBC matrix - a status setter is the named enabler ([f45]). No standalone
equivalence scenario ([f44]).

**hazard.** No flat-C token and no engine code of its own ([f6]); no
`rbart_vi()`, `xbart()` or `dbartsSpec()` reach. Out of SBC by design. One
dedicated tinytest file.

**hurdle.** No `dbarts()` sampler by construction ([f12]); no `rbart_vi()`,
`xbart()` or flat-C reach. No top-level pointwise loglik ([f25]). No warm start
/ grow-from-root. A heteroscedastic positive part is REFUSED via the occupancy
component's own gate, deliberately, not a partial feature ([f34]).

**bcf.** No `bart2()` surface for the NAMED `bcf()` causal verb, by the
resolved fork that puts it in bartCause; the general K-forest amplitude
capability itself IS reachable from `bart2()`'s formula interface via a
`forest()` term since S12 ([f7]). No pointwise loglik. Warm start and
grow-from-root are unrefused and untested for two forests ([f36]). Whole-data
`setData` stays undesigned (door 1 of the model-space survey).
`setForestWeights` now has a public R5 method (multiforest-extension-surface
M1, 05ac3b4b; [f26]) - no longer a gap. The probit and logistic arms (M4.4) have
no equivalence scenario, no SBC arm and no measured active-rows mask ([f48],
[f26]); aft, ordinal and nbinom are recorded doors, not gaps.

**grouped.** `rbart_vi()`-only surface ([f8]): no `bart2()`, `dbarts()`,
`xbart()` or `dbartsSpec()` reach, and no `warm.start` / `n.grow.sweeps` formals
([f37]) though the engine paths carry no group gate. The `setResponse` gap
CLOSED at adoption-slate S3 ([f13]; the section-2 cell is the record).
Composition with a variance forest constructs unrefused and untested ([f30]).

**heteroscedastic.** No `xbart()` reach. Pointwise loglik ignores the
per-observation `s(x)` surface it stores ([f28]). Selecting attribute
undocumented in the header ([f3]). Out of the SBC matrix, deferred not blocked
([f47]). One recorded unbuilt door from its own arc: the `setState` variance
column-mask gap (variance-forest-mutation-routing.md:500-501).

**Cross-cutting.** `nameable-calibration` is ARC COMPLETE, all four slices
LANDED (S0 4c866286, S1 c2a7e89b, S2 d809b944 + 7da36dc3, S3 at
dbarts-h-reshape S1, ab3aa2fa) - the flat-C entries shipped and gate nothing
further; the R surface was already complete before S3, since the flat-C gap
never gated a `dbartsSampler` caller. `latent-subset-mask` is likewise ARC
FULLY COMPLETE (S0 dc11a805, S1 6db22aee, S2 87d370ea, S3 8b047f8b, S4
93afd635, flat-C entry at dbarts-h-reshape S1, ab3aa2fa) - the mask covers
every response family and every surface, R and flat C alike; table 3 now
carries no `P` cell.
