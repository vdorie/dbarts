# Checking sheet for the decision register

One line per rubric item 8 to 12, checked against docs/decisions.md as
installed, the old table at HEAD, and the code at d36910af. Items 13 to 15 are
the maintainer's to judge and are not graded here.

## Verdicts

**Item 8 - each entry readable on its own.** PASS. Every entry names what was
decided, the alternative, what a user of the package notices, and the ruling
where there is one. No slice label, finding code or codename survives; the two
entries that closed "not put to VD" now say "not put to the maintainer"
(dec-B119, dec-B120).

**Item 9 - a ruled entry shows the ruling as the current state.** FAIL as
received, in seven entries; fixed. dec-A61 kept the superseded reason as the
live one ("no gate could notice a shift against the released package", which
rested on the deleted engine); dec-A43 and dec-A46 still said the tuning
constants could not be moved after the ruling that exposed them; dec-A27 left
the residual-law bundle as a conditional ("if that pass leaves the argument in
place") after the pass; dec-B91 gave the ruling without saying it had landed;
dec-B121 stated an unbuilt accessor as current behaviour; dec-A31's Marked line
dropped the supersession by dec-B119.

**Item 10 - claims about the package checked against the code.** FAIL as
received, in nine sentences; fixed. dec-B75: "move to a new function whose name
is still open" - the function ships as [`bartBT`](../../../R/bart.R), exported
in NAMESPACE. dec-A43: "prediction below ten million cells stays serial" and "No
caller can move any of them" - the cutoff is
dbartsControl(predictParallelCutoff = 50000L) and all three are control slots. dec-A46: the categorical limit is
categoricalExhaustiveCap on the control. dec-A27: resid.dist and both bundles
are gone; student(df) is a family. dec-B121: the variance-surface accessor is
not among dbartsSampler's methods. dec-B86: "offset ... wait for a consumer" -
dbarts_sampler_setOffset ships in the header. dec-A63: "61 percent of the 8539
actually measured" - no such count is recorded anywhere in the tree and the
workflow's own comment reads differently; the number is gone and the gap is
stated in words. dec-A70: "the downstream package insight ... works again" -
not measured here, softened to what the sweep recorded. dec-B50: "the earlier
ruling that within-chain threading ships" - dec-A37 is the opposite ruling, that
threads run chains only.

**Item 11 - "cannot" claims checked against everything to hand.** FAIL as
received, in three sentences; fixed. dec-B15: "A user cannot fit degrees of
freedom below 3" - student(df = 2) is accepted and fits (R/family.R validates
only positivity; the floor of 3 is the estimation grid's, model.hpp). dec-A04:
"Nothing can now match another package's generator stream" named no route; the
engine keeps user-supplied uniform and normal arms in the support library and
the entry now says that nothing in R reaches them. dec-A61: "nothing wider can
be recorded" is disproved by the 26-scenario comparison against an installed
0.9-34 (dec-B120), and the entry now reads as superseded.

**Item 12 - numbers only where they change what the reader would do.** FAIL as
received, in seven entries; fixed. Removed as count-justifying: dec-A01's "nine
S3 methods ... six of those"; dec-A25's four token tallies; dec-A29's "ten
families"; dec-A33's "five setters"; dec-A20's second count of the same
disjunction; dec-A63's 8539 and 61 percent. dec-A42 carried a superseded
estimate (3 to 4 percent) beside the measurement that replaced it; the estimate
is now named as an estimate. Numbers kept are the ones a maintainer would act
on: cap and limit values, measured speed changes, line and byte counts that
price a removal, and the z of the cutover comparison.

## Check lines resolved

The six "Check:" lines the writers left are resolved against the code, git
history and the cited records, and deleted.

- dec-A36. Both halves explained in the entry: the record calls the drop the
  maintainer's own but quotes no words and puts no fork, which under the
  register's own attribution rule is what places the entry in section A; the
  commit the record names is the pre-rebase one, and on this branch the landing
  and the drop are 9fe39856 and 99b356d8.
- dec-A42. The 3-to-4-percent figure was the estimate that motivated the work;
  the measurement that closed it put the same kernels under one percent on
  arm64 and within noise on x86. The measurement stands and the estimate is
  labelled as one.
- dec-A44. Not a contradiction. The bridge sets each forest's ridge from
  whether that forest's amplitude prior carries a positive half-Cauchy scale
  ([`applyAmplitudeSpec`](../../../src/R_interface_bartcore.cpp)); the
  treatment forest's prior is a fixed variance, so the move never turns on. The pair of ridge flags in the engine
  struct belongs to a fixture initializer the shipped path does not reach.
- dec-A48. DESCRIPTION requires R 4.2.0 and nothing enforces 4.3;
  docs/design/core-generalization.md still says 4.3 for the C++20 toolchains,
  which the entry now records as a stale statement rather than a second floor.
- dec-A61. Superseded, as above.
- dec-B59. Not a contradiction: the two rules cover different additions, and
  the header states both. A function added after 1.0-0 arrives under a new name
  with a minor bump; a struct still grows by appending a field, and that append
  bumps the minor version and re-bakes the hash.

## What changed from the writers' text

Twenty-six entries were edited. The three largest:

1. dec-A61 was rewritten end to end. It had kept a false premise as its
   reasoning; it now says the premise was wrong, marks the entry superseded,
   and keeps the part that survives - that no gate fires on a shift against the
   released package, because the wider comparison is run by hand.
2. dec-A43 and dec-A46 were rewritten to the state the ruling left: the
   thresholds and the categorical enumeration cap are control settings, the
   prediction cutoff was recalibrated to fifty thousand cells, and what stays
   fixed (cut-code width, leaf-regression columns, perturb width) is named.
   Their headings changed with them.
3. dec-B15 and dec-A04 lost their unqualified "cannot". Each now names the
   route that does exist - a supplied Student-t degrees of freedom below 3, and
   the engine's user-supplied generator arms - and says why it does or does not
   help.

The rest are single-sentence corrections: dec-B75 names bartBT, dec-B121 says
the accessor is ruled and unbuilt, dec-B91 says which settings landed, dec-B86
says test offset, dec-B50 names the entry it sits against, dec-A27 says the
consolidation happened, dec-A70 drops an unmeasured downstream claim, dec-A31
carries its supersession on the Marked line, dec-A63 drops an unrecorded count,
and the register's opening keeps the old header's ordering note and names the
Marked line. The writers' register is otherwise untouched.
