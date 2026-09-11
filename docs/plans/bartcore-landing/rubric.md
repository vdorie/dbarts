# Rubric for maintainer-facing documents

The standard a document must meet before the maintainer reads it. Written
for the bartcore landing memo; it applies to any document the maintainer
will act on. Every item is pass/fail. A reviewer reports each failure as
the quoted line plus the item number, and the completed grade sheet is
committed beside the document. A refresh regrades every item, not only
the lines the refresh touched. "Improve clarity" is not a finding.

The reader is the package maintainer: a statistician who knows BART, R
and dbarts 0.9-x, who has read nothing under docs/, and who reads the
document once, in one sitting, to make a decision. They adjudicate the
decision register with the document open, so a decision bullet may end
with its row id. The document is graded on what that reader can do after
one reading, not on what it contains.

## What the document must do

1. The opening paragraph states the decision being asked for, the
   writer's recommendation, what the release gains, and what the reader
   would have to learn to change the recommendation.
2. Directly after the opening, the facts that could change the
   recommendation are listed, each stated whole, ordered by how much of
   the release each one blocks. If a fact is split across sections, it is
   stated whole where it bears on the decision and referred back to
   elsewhere.
3. The body answers the questions the decision turns on. For the landing
   memo: whether the branch can merge and in what order with the consumer
   packages; what breaks for whom at each intermediate state, and whether
   loudly (a compile or load error) or silently (a wrong result); how much
   is known about whether the release computes the same answers as
   0.9-x, and where that is not known; what the agents decided without
   you and what each costs; and what is unfinished, abandoned, or proven
   only by a gate not run at the tip, with who does it and before which
   step.
4. Every set the document reasons over (consumer packages,
   prerequisites, gates, intermediate states) is listed in full before it
   is used.
5. A count appears only with what it counts and where the list is.

## What must be true

6. The reviewer re-checks against the tree at the tip every claim about
   what the package exports, refuses, errors on, returns or defaults to,
   and what any gate or check does and when it runs. A register row is
   evidence, not truth: where they disagree, the document states the tree
   and the row is corrected in the same commit. For any claim using
   every, only, none or nothing, the reviewer names the case that would
   break it and confirms it does not exist.
7. Every decision named traces to a decision-register row and is
   described in its standing form. A row marked superseded is described
   from the superseding row, with the reversal named. A row you have
   marked as yours is not presented as agent-made, and a decision is
   presented as yours only when its row says so. Decision bullets are
   ordered by cost.
8. Every completeness claim states whether its gate ran at the tip and
   names the commit or run it did prove. Every other claim has a register
   row, every cited row has a code or git anchor, and no claim rests only
   on a design or plan document.
9. The body, the tables and the appendix agree on every fact they share.
10. The evidence appendix maps each claim to register ids that exist at
    the tip and whose status supports the claim.

## How it must read

11. The reader is addressed as you and the writer speaks as I. The phrase
    "the maintainer" does not appear. Every action names who does it.
12. Ordinary words, one idea per sentence, each with its verb. A term
    coined for this project is replaced by plain words or defined where
    it first appears. Anything the reader has not met in a released
    dbarts, in R, or earlier in the document is glossed; anything they
    have met is not. A register id never appears inside a sentence as its
    subject or object; it may close a decision bullet, and otherwise
    lives in the appendix.
13. The document describes the state at the tip, not how it got there:
    no review provenance, finding codes, agent names, or dates except a
    maintainer ruling quoted from the register.
14. Facts and consequences only. No judgement of quality (robust, clean,
    thorough, careful, ensures, and their kin). No general verb (carry,
    hold, cover, admit, discharge) where a specific one (needs, costs,
    breaks, fails) would say it. No paragraph whose deletion would leave
    the reader's judgement unchanged.
15. A table only where each row differs in every column and every cell is
    a phrase; otherwise short headed sections. Every silent-breakage case
    says what the wrong answer looks like.
16. The body is at most 3000 words. Appendices hold the evidence mapping
    and any list the body counts, and nothing else.
17. ASCII only, no em-dashes, no arrows.

## Exit test

18. Run once before the document is sent, by a reader matching the
    profile who has not read the registers. With only the document, they
    write down what they would do, the recommendation, the facts of
    item 2, and the answers to item 3, without re-reading any sentence,
    and list every term they had to guess. Each answer they cannot
    produce, produce wrongly, or produce only after re-reading is a
    numbered failure, and "I cannot decide" is a failure. Their answers
    are then checked against the tree and the registers, and any
    disagreement is a failure of the document, not of the reader.
