Grade sheet: landing memo

This sheet grades memo.md, as committed alongside it, against rubric.md at the
same commit. Every item is pass or fail; each failure is given as the quoted
line. 17 pass, 1 run separately.

## 1. Opening: decision, recommendation, gains

PASS. The first paragraph asks for three named decisions, recommends merging
today and leaving the candidate, and "What the release gains." follows it.

## 2. The facts that could change the recommendation

PASS. Six facts stand directly after the opening, each whole, descending from
the one that bears on the merge itself to the two that block nothing. Each says
in its own last sentence what it blocks. The facts picked up again later
(the nine scenarios, the leaking sites, calibration) are referred back, not
re-argued.

## 3. The questions the decision turns on

PASS. All six groups under "What is not finished, and who finishes it" name an
owner, a step, or both: the automation gaps now read "both yours and in no step
below", and the five other groups name you or me and the step each is due
before.

## 4. Sets listed in full before use

PASS. The one use of a step before the list names the step in words, "before
the CRAN submission, step 5 below". The other sets are listed where they are
first used: the four consumer packages and the five schedule-only workflows in
the opening, the seven per-push workflows and the six steps at the head of
their own sections.

## 5. Counts and measurements

PASS. Every count in the body arrives with what it counts and with its list,
inline or in Appendix B: the eight families, the six loud refusals, the 38
behavioural changes split 34 and 4, the seven most-met changes, the nine
scenarios, the 27 CRAN packages, the 32 and 22 decision sets, the five checks
and the 7680 fits decomposed into their 48 cells. Each measurement carries what
it means for the decision: the z of 3.83 with the harness's own cutoffs and with
what its z does not mean, the coverage figures with the verdict they support,
the speed and memory figures with what each compares.

## 6. Claims against the tree

PASS. Re-checked at the tip: the nine CI runs at the merge commit, all success,
calibration green on all seven arms; the five schedule-only workflows' triggers
and their runs on this branch; the eight family tokens beyond gaussian and
probit; the two-part lognormal's matrix-only refusals; the monotone constraint's
refusal of a linear or Gaussian-process leaf; the 25 entries in the header's
entry list; the probit weight refusal for any weight other than 1; the auto
family's routing off the response class; the 25 exact scripts and their 23, 1
and 2 split; the 52, 12 and 11 bitwise scenarios; the harness's warn-above-3 and
fail-above-4 cutoffs; the five check platforms and the Windows ARM64 leg; the
change move's missing proposal-density ratio on main; the chi hyperprior's entry
at 0.9-10 and the de-scaling site's at 0.9-32; the 38 behavioural rows, which
match Appendix B exactly. The claims using every, only, none and nothing were
tested for the case that would break them and none exists.

## 7. Decisions in their standing form

PASS. The feature-matrix run now reads "That run is left manual." and is
attributed to no one, which is what the completion row supports. The rest of
the item is met as before: the eight bullets are in the register's cost order,
each in its row's standing text; every ruling given as yours is a row that says
so; and the superseded hyperprior row is described from the row that supersedes
it, with the reversal named.

## 8. Completeness claims, gates and anchors

PASS. All eight costliest decisions are anchored in code. The eighth bullet
now states the multinomial defaults themselves, and the appendix anchors it in
the "every 'counts' row must have at least one trial" refusal in R/data.R and
R/A_class.R. Every completeness claim states whether its gate ran at the tip,
gates are named by what they check, and the commit and run ids sit in the
appendix and nowhere in the body.

## 9. Body, tables and appendix agree

PASS. The five schedule-only workflows, the seven per-push workflows, calibration
at the tip, the smoke test's branches, the 34 and 4 split, the six steps and
their owners, and the CRAN set of three all read the same in the body, the
tables and the appendix.

## 10. The evidence appendix maps the document's claims

PASS. Every cited id exists at the tip and its status supports its claim, and
every row maps a claim the document makes. The row about the commits since the
last engine change is gone with the sentence it served, and the insight row
went with the insight sentence.

## 11. You and I, and who does each action

PASS. "the maintainer" does not appear. The steps name owners, each unfinished
item names you or me, and the two unowned actions, deleting the stubs and
writing the CRAN sweep's script, are stated as unowned.

## 12. Plain words, glosses, ids, actors

PASS.

Register ids inside sentences: none. All eight in the body close a decision
bullet, which the item allows; the rest are in the appendix. The one sentence
naming the prefixes introduces them as prefixes, not as rows.

Terms: the branch name, the three registers and their prefixes, each of the
five scheduled checks, the per-push exact and bitwise gates, the submission
battery, chi(1.25, Inf), lorax and the step numbered before its list are all
glossed where they first appear; the coined "tombstone" and "door" are gone
from the body, and bartBT is glossed where the appendix first uses it.

Actors: "I renamed dbartsSampler$run's fourth argument and added formals to
dbarts()" and "To back out you reset dbarts main" name the actor as the
subject; no other sentence separates them.

## 13. The state at the tip

PASS. No date appears in the body: treatSens now reads "which CRAN archived".
No commit hash, run id, finding code or agent name appears there either.

## 14. Facts and consequences only

PASS. "honour" is now "read", "takes a manual start" is "runs on a manual
start", "sits on" and "takes the modern door" are gone, and "reaches them" is
"triggers them". The two remaining uses of "cover" are the statistical sense,
coverage of the true probability. No judgement of quality appears, and no
paragraph was found whose deletion would leave the reader's judgement
unchanged.

## 15. Tables

PASS. No two rows of any table share a cell in any column, and no cell carries a
second complete sentence. Every silent-breakage case says what the wrong answer
looks like: bairrtt's moved posterior, the dropped engine seed, the stale
binary's misread slot, and each of the seven most-met changes.

## 16. Length and appendices

PASS. The body is 2998 words. Appendix A is the evidence mapping; Appendix B
holds nine lists, each one the body counts, and nothing else.

## 17. ASCII only

PASS. No byte outside 0x20 to 0x7e, no em-dash, no en-dash, no arrow. The only
double hyphens are command flags in Appendix A and table rules.

## 18. Exit test

Run separately; see exit-final.md.
