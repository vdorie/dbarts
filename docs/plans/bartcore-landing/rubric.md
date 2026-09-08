# Rubric for maintainer-facing documents

Written for the bartcore landing memo; the standard for any document the maintainer reads, including the docs/ triage rewrite. Item 14 is the exit criterion: a fresh reader with only the document answers its questions correctly.

The document under review is named at the top of each review. Its reader is the
package maintainer (VD): a statistician who knows BART, R, and dbarts
0.9-x, has not read anything under docs/, and will read the memo once,
in one sitting, to decide whether to merge the bartcore branch onto main.
Its sources of truth are docs/decisions.md and
docs/plans/bartcore-landing/{changes,completion}.md (the registers).

Every item is pass/fail. A reviewer reports each failure as the quoted
line plus the item number. "Improve clarity" is not a finding.

## Content

1. The memo answers exactly four questions, in this order, each under
   its own plain noun-phrase heading: (a) can bartcore merge onto main,
   and in what order with the consumer packages; (b) what breaks for
   GitHub and CRAN installers in each intermediate state; (c) what did
   agents decide that the maintainer did not, and what does each cost;
   (d) what is unfinished, abandoned, or proven only by a stale gate,
   and why.
2. Every completeness claim traces to a completion-register row whose
   gate evidence is not "none".
3. Every decision named in the memo traces to a decision-register row,
   states who made it per that row's evidence class, and names what it
   forecloses or costs. No decision is presented as the maintainer's
   unless its row is class quote or choice.
4. No claim rests only on a prior design or plan document. Each has a
   register row, and each register row has a code or git anchor.
5. The merge order and the installer breakage appear as one table each.
   Every intermediate state names who is broken and whether the failure
   is loud (compile or load error) or silent (wrong result).
6. Nothing in the memo is process history: no "we first tried", "now",
   "no longer", "retired", review provenance, run numbers, finding
   codes, hashes, agent names, or dates except a dated maintainer ruling
   that a decision row quotes. Exception: in the decisions section only,
   a decision-register id in a short parenthetical such as "(A06)" is
   allowed, because the maintainer adjudicates that register with the
   memo open. Completion-register ids appear only in the appendix.
7. Nothing in the memo is a judgement of quality: no robust, clean,
   comprehensive, thorough, careful, ensures, seamlessly, leverages,
   elegant, powerful, or their kin. Facts and consequences only.

## Form

8. At most 2500 words in the body, at most four numbered sections plus
   one lettered evidence appendix. No other appendix. The appendix maps
   memo claims to register ids and holds nothing else.
9. Every named thing (function, argument, package, gate, branch) is
   either something the maintainer would type in R or at a shell, or is
   glossed in plain words on first mention. No codename, slice label,
   gate name, or plan title appears unglossed. No path into src/, R/,
   tests, or benchmarks appears in the body.
10. No em-dashes, no arrows, no parentheticals longer than a few words,
    ASCII only. Sentences about 20 words, one idea each, with a verb.
    No semicolon-joined clauses. No bullet that runs longer than two
    sentences. Bold only the first few words of a bullet, never a
    sentence.
11. Numbers appear in tables or on their own line, and only where they
    change what the reader would do. No number in running prose except
    a version number the reader must type.
12. No triplets for rhythm, no "not X but Y" framing, no headings that
    argue, no closing summary or offer, no restating an earlier section.
13. Per-paragraph test: if the paragraph were deleted, would the reader
    judge anything differently? If not, it fails.

## Cold-reader test

14. A reader who has only the memo, no repository, can answer the four
    questions in item 1 and their answers agree with the registers.
    Any disagreement is a memo failure, not a reader failure.
