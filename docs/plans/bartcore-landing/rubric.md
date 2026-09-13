# Rubric for documents the maintainer reads

Applies to the landing memo, the decision ledger, and any brief written for
the maintainer to act on. The maintainer reviews this rubric before it is
used and can change it at any time.

## Who reads and why

The reader is the maintainer: a statistician who knows BART, R and dbarts
0.9-x, has not read anything under docs/, and reads once to decide what to
do next. The document exists to show what the maintainer needs to know and
what the agents assumed on their behalf, so that gaps can be seen.

## Before writing

1. The writer puts the questions the document will answer to the maintainer
   in chat, as a short list, and gets them confirmed. The confirmed list
   opens the document. Nothing else appears in the body.
2. Every version is written from a blank page against that list. A refresh
   is a rewrite, not an edit of the last version. A section that answers no
   question on the list is not written.

## What the landing memo answers

3. What 1.0-0 is, in one paragraph: what it adds, what it repairs, and how
   it is faster, with the one or two numbers a user would quote.
4. What breaks for someone using 0.9-34 as they use it today. For each
   break: whether it is loud (an error at install, load or the call) or
   silent (a different answer), and for a silent one what the different
   answer looks like. The most common cases come first. Downstream packages
   are covered the same way, by package.
5. What is missing or wrong that should be fixed before a release
   candidate, and before 1.0-0: features not built, known defects,
   measurements not made, and downstream packages not yet working. Each
   item says whether it goes before or after the merge to main and why in
   one clause.
6. How much is known about whether 1.0-0 agrees with 0.9-34, and what
   it would take to know more.
7. The work ledger: every open item, with its status. Statuses are done,
   open, and deferred with a stated reason. Nothing is abandoned; an item
   with no home gets one in TODO.
8. Not in the memo: how to merge, how to submit to CRAN, who does which
   step, run ids, commit hashes, counts of register rows, tallies of who
   decided what.

## What the decision ledger must do

9. Each entry can be read on its own by the reader above: what was decided
   in plain words, what the alternatives were, what a user of the package
   notices, and the maintainer's ruling if there is one. Codenames, slice
   labels and finding codes are unpacked or dropped.
10. An entry the maintainer has ruled on shows the ruling as the current
    state. The earlier agent position appears only as the alternative that
    was not taken.

## What must be true

11. A claim about what the package does, refuses, returns or defaults to is
    checked against the code at the commit the document describes.
12. A claim about what can or cannot be done is checked against everything
    the maintainer has to hand: git history, CRAN, R itself and the sister
    packages. "Impossible", "nothing can", and "no fix" are used only after
    the obvious route has been tried or named and ruled out in the text.
13. A number appears only where it changes what the reader would do, and
    then with what it means. Lists the reader must act on go in an
    appendix; lists that only justify a count are not included.

## How it must read

14. Plain prose a statistician follows in one read. Headings are phrases
    the reader could say aloud. Each sentence has a clear subject, usually
    the package, a function, a user, or the maintainer. Anything not met in
    a released dbarts, in R, or earlier in the document is explained where
    it first appears. Register ids may close a ledger entry or a bullet;
    they do not appear inside sentences.
15. The document describes the current state, not how it got there. No
    review history, agent names or dates other than a maintainer ruling.
16. The body of the landing memo is at most 2000 words. ASCII only.

## Checking

17. An agent checks items 3 to 13 against the code and the ledger and
    leaves a written sheet beside the document, one line per item, naming
    any failing sentence. An agent's pass on items 14 to 16 does not count;
    readability is judged by the maintainer's read, and what the maintainer
    flags is fixed by rewriting the section, not the sentence.
