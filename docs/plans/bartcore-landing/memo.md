# Landing memo: bartcore

## 1. Merge feasibility and order

main carries no commit the branch lacks. The merge brings exactly the branch's
own changes, and nothing has to be reconciled.

One gate stands in front of the merge itself, and it is the maintainer's own:
reading the branch review and declaring a release candidate. The CRAN
submission carries further prerequisites. Some were recorded when the
grouped random-effects fit was removed. Another is a CRAN package, lorax,
whose example fits bart() with a factor response that 0.9-x silently coded
as numbers. 1.0-0 refuses it, so its author has to adjust the example.

The order is forced. This branch ships a C header that dbarts main does
not carry, and some compatibility branches already include it. Every
consumer version floor names a dbarts version that no branch of dbarts
satisfies until the merge. dbarts merges first, and the only free variable is
the consumer mains' lag.

stan4bart's compatibility branch replaces its sampler back end wholesale.
bartCause and bairrtt call dbarts only from R, and bairrtt needs no branch.

The exact-hash flag below is a compile-time switch that folds the header's
entry-point signatures, enumerations, callback and struct layouts into a
token baked into each consumer at build time. A mismatch against the
installed package's own token is refused loudly.

| step | action |
|---|---|
| 1 | Pre-merge, pushing to no main: re-bake the header hash, reinstall and test each consumer, push compatibility branches only |
| 2 | Merge bartcore into dbarts main and push. The GitHub window opens |
| 3 | In the same sitting, merge each consumer's compatibility branch into its own main. The GitHub window closes |
| 4 | Point the consumer-build check at each consumer's main, take the consumer workflows off their branch pins and the branch notes out of the readmes, keep the temporary install steps |
| 5 | Clear the three release prerequisites, then submit dbarts 1.0-0 to CRAN. On acceptance, submit stan4bart 0.0-14 and bartCause the same day. The CRAN window opens |
| 6 | When CRAN carries all of them, delete the temporary install and reinstall steps |

CRAN serves dbarts 0.9-34 until the submission is accepted, so the temporary
install steps must survive a consumer resolver installing that older copy.

The first step's obligation is not yet discharged: every consumer must be
re-verified against the final re-baked header. bartCause hands out a thread
budget the new default would warn about, and treatSens reaches dbarts
internals by name. Each has been run only against an intermediate header,
none against the final one.

The exact-hash flag is already out of every consumer file, so the version
pair is the only guard, and a consumer built against a different header
sharing that pair is admitted. A check on dbarts's own side replaces the
flag: a build fails when the baked token moves and the version pair does not.

## 2. Installer breakage in the intermediate states

A window is a span in which an installer gets a broken pairing. Loud means a
compile or load error. Silent means the call succeeds and returns a different
answer.

| state | who breaks | loud or silent | what happens |
|---|---|---|---|
| after the dbarts merge, GitHub | a 0.9-x user script | loud | rbart_vi, a three-element burn-in, a factor response of three or more levels and fractional counts error |
| after the dbarts merge, GitHub | a 0.9-x user script | silent | the proposal mixture, the binary k prior, categorical factors, missing-data incorporation and combined chains change the answer, and an unnamed-argument fitting call takes bart(), the main fitting function, at the new defaults |
| after the dbarts merge, GitHub | a 0.9-x user script | warns once per session | a second positional argument to setResponse means updateScale. A renamed argument, an inert thread method and a call in the old vocabulary each name their successor and then do the right thing |
| after the dbarts merge, GitHub | stan4bart main, treatSens main | loud | a source install fails to compile against a deleted header |
| after the dbarts merge, GitHub | bartCause main | mixed | loud on the grouped route and the cross-validation call. Otherwise it installs, records one chain, and returns draws moved by the default changes |
| after the dbarts merge, GitHub | bairrtt main | silent | it installs and every signature it calls is unchanged, and its draws move with the defaults |
| after the dbarts merge, a CRAN copy | stan4bart 0.0-13 | loud | it fails at load, R reporting a function not provided |
| after the dbarts merge, a CRAN copy | bartCause 1.0-10 | mixed | loud on the grouped route and cross-validation, silently one chain otherwise |
| after the dbarts merge, a local copy | treatSens 3.0 | loud | it fails at the first analysis call, its lookups being lazy |
| after the consumer merges, GitHub | unverified | unknown | no consumer has been built against the final header, and stan4bart's branch carries an open defect |
| after CRAN accepts dbarts 1.0-0 | stan4bart 0.0-13, bartCause 1.0-10 | as above | a package update alone breaks them until the consumer submissions are accepted |
| after CRAN accepts dbarts 1.0-0 | lorax 0.1.0 | loud | its example and tests error on a factor response until its author changes them |
| after CRAN accepts dbarts 1.0-0 | insight 1.5.4, off CRAN only | loud | its test expects NA-response rows dropped and gets an error. The test is skipped on CRAN |
| after CRAN accepts the consumers | nobody | neither | closed |

The open defect is stan4bart writing n.cuts, an integer slot here, through a
control attribute, so a caller passing a double breaks on the next
assignment.

A break survives both windows: stan4bart matches control arguments against
the control constructor's formals, dropping the renamed seed and losing the
caller's engine seed. The filter is unchanged on the compatibility branch.

User-visible changes that return a different answer: 36

Neither window can be eliminated, only shortened: the GitHub window is as
short as a sitting, and the CRAN window stays open because no published
consumer bounds dbarts above and CRAN processes submissions serially.

The removals and renames behind the loud rows are agent-made, and this
breakage is their cost, as are categorical factors, combined chains and the
binary k prior. Missing-data incorporation and the proposal mixture carry
maintainer quotes.

## 3. Agent-made decisions and their costs

The registers behind this memo are files of tables, one row per decision,
change or open item. A decision is the maintainer's only where a row quotes the
maintainer choosing, or lays out alternatives and names the pick. An approval
marker, a standing grant, or no record at all leaves it agent-made.

Agent-made decisions carrying an identified cost: 70

Decisions with maintainer evidence: 117

- **Grouped random effects** (A01): rbart_vi and its S3 methods are deleted,
  breaking bartCause main at load. The maintainer confirmed the removal on a
  condition (B105): stan4bart clears its mixing bar, or an R-only grouped
  intercept stands in.
- **Equivalence evidence** (A61): the bitwise check compares only against a
  baseline this branch recorded, and the deleted engine forecloses recording
  any wider one.
- **Fixed engine constants** (A46): the leaf regression column cap, perturb
  width and cut width stay reachable from no argument. Other engine limits
  became settings (B91), kept control-only under B116, reached from bart()
  and xbart() through control = dbartsControl().
- **Factor predictors** (A06): bart(), the modern fitting function, defaults
  factor predictors to single categorical or ordinal columns, while the
  BayesTree-style bartBT() keeps indicator expansion (B78). The two fit
  different models on the same data by default.
- **Binary node hyperprior** (A07): the default moves, and so does the scale
  default of chi(), so every probit fit's posterior moves. The maintainer
  keeps the new default and has scheduled the evidence it rests on (B106),
  which is not yet in the tree.
- **Combined chains** (A08): bart() defaults combineChains to TRUE and
  flattens results chain-major, so code indexing the chain margin breaks.
  The recorded decision was a default of FALSE family-wide.
- **RNG surface** (A04): rngKind and rngNormalKind are dropped and rngSeed
  becomes seed, warned through the next release. Only Mersenne-Twister is
  reachable, so no fit can match another package's stream.
- **No deprecation cycle** (A02): the maintainer's rule (B76) keeps a removed
  or renamed name reachable through the next release as a marker naming its
  successor, and most ship. rbart_vi still errors, and the old C++ interface
  vanishes with no marker.
- **Cross-validation redesign** (A05): xbart drops its control argument and
  replaces its threaded fold loop with a cluster, breaking a consumer that
  passes a prebuilt control. control later returns as a formal under B116,
  merging a prebuilt control in under that formal's rule instead of breaking,
  unless it carries fit-state attributes.
- **Family detection** (A16): family = "auto" reads a factor response and
  routes it to probit, ordinal or multinomial, where 0.9-x fit gaussian on
  the same level codes.
- **Silent state and a silent cap** (A14, A10, A13): mutators store state
  only on an explicit updateState = TRUE, so a later save can write stale
  state. The cap on the sampled k is silent, unlike the cap on initial-forest
  draws, which errors.
- **Chains-only parallelism** (A37): a single-chain run gets no sampling
  parallelism at any thread count. Closed on remeasurement (B115): n.threads
  is capped at the chain count by default, and a larger budget only warns.
- **Scalar draw path** (B90): the maintainer's ruling keeps a configure flag
  selecting a scalar, fixed-order reference build for development, gated
  statistically before release. The vector kernels that would make the
  shipped build differ showed no meaningful gain and did not ship (B113),
  so both builds stay identical code, with a fused pass over the weighted
  families landed instead that moves those draws in the last bits.
- **The published C boundary** (A29): the header exposes calls over an opaque
  handle, families selected by string-named attributes. The maintainer's
  ruling (B84) keeps it pure C: a consumer builds the sampler through R, and
  no non-R host can create one until a post-release plan (B85).
- **Install floor and platforms** (A48, A49): the R floor rises from 3.1-0 to
  4.2.0, and Solaris and big-endian support are removed, excluding every
  earlier R release that worked under 0.9-x.
- **Standing discretion grant** (A65): a standing permission let user-facing
  forks be settled without reaching the maintainer, so a class of them
  carries no maintainer attribution.
- **setPredictor partition** (B117): the maintainer chose to land a faster
  re-partition step after dbartsSampler$setPredictor swaps a predictor
  column, at the cost of a re-record of the equivalence baselines for every
  scenario that calls it, since the change moves those draws in the last bit.

The register marks the rows below as the maintainer's, whose evidence covers
less than the row claims. Confirm them first.

- **Gaussian-process leaves** (B22): the ruling keeps them, leaving open
  whether they belong in the release, so an engine component and a man page
  freeze in. B110 narrows this: a fit that falls back to ordinary leaves
  warns.
- **Release-review rulings** (B40, B41): recorded as rulings in a batch, with
  no quoted words. One refuses a mid-chain rescale of a heteroscedastic fit,
  the other superseded by the fallback warning above.

One row the register marks as the maintainer's does not hold up: the quote
behind large-n workloads (B50) is an agent paraphrase the maintainer
disowns, leaving no shipped path for that priority.

## 4. Unfinished, abandoned and stale-gate work

Unfinished, with no decision recorded: 49

Abandoned: 4

Blocked: 8

Deferred with no maintainer evidence: 25

Deferred by the maintainer: 10

Done, but proven only by a gate that predates the tip: 5

Sites in the C API hold C++ objects
that own heap memory across R's error longjmp and leak, unreachable by an
ordinary refusal. The memory-check gate runs the suite under a leak detector
but has only run once, so the merge is its first real test, and it has never
covered the exported data-augmentation helpers. bartCause has no per-push
check of its own and carries the silent failure mode. treatSens relies only
on the monthly reverse-dependency smoke test.

Most abandoned items carry no decision: a non-conjugate move strategy for
Gaussian-process leaves, a causal-forest treatment-scale rescale, and a
slate of post-release flat-C readers, none of them in the TODO file. Another
item, a cross-repository consumer build under the sanitizers, was added and
deleted the same day under an agent-made decision.

The maintainer's own read of the branch and the release-candidate declaration
block the merge, and also block registering the scheduled workflows for
memory checks, equivalence, calibration, protect balance and the reverse-
dependency smoke test, none of which can fire from a non-default branch. A
nonzero share for the leaf-pair rule draw was rejected: it clears the
sample-size bar but regresses coverage past the margin, so the maintainer
keeps it at zero through the merge, its weight to be set later by the mixing
research. A bartCause route for grouped data is blocked because the deleted
function was its only path, and stan4bart owes a mixing bar on the
random-effect scale that the same removal named as a release prerequisite.

The deferred rows without maintainer evidence include a Python binding whose
stated premise the code contradicts, and a size threshold below which the new
fused per-sweep pass should be declined, which cannot land until the bitwise
corpus carries a fit large enough to exercise it. The rows deferred by the
maintainer include real-valued negative-binomial dispersion and arbitrary
real binary weights, both waiting on an approximate data-augmentation draw
the maintainer declined.

The stale gates prove their items only at a commit before the branch tip. The
protect-balance check, which verifies that R objects stay protected from
garbage collection, failed, and its fix carries a local re-run alone. The
check that the sampler recovers parameters drawn from its own prior sits
earlier still. The equivalence comparisons, the exact-posterior gate
scripts, the test suite, the package check, the C++ tests and the sanitizer
builds have since rerun and passed at the branch tip, closing those gaps.

The CRAN reverse-dependency sweep is one of the stale gates. It ran at an
earlier engine, with no script in the tree, and found the breaks named
above. More packages have declared dbarts since and have never been checked.

A record compares this engine against the old one, statistical rather than
bitwise, taken at the engine cutover. Its manifest marks it evidence and not
a gate.

Scenarios it covers: 9

It has never been re-run against the current engine, so a drift inside the
statistical band would be caught nowhere.

Every entry point in the shipped C header is exercised by a compiled consumer
test, including the version and hash mismatch handshake.

## Appendix A. Claims and register rows

Prefixes: dec- is the decision register, chg- the changes register, cmp- the
completion register.

| claim | rows |
|---|---|
| main carries no commit the branch lacks | git rev-list --count bartcore..main |
| the merge gate and the three release prerequisites | cmp-L01, cmp-L05, cmp-L06, cmp-L07, cmp-L08, dec-B74 |
| the insight break | cmp-L09, dec-A70 |
| the reverse-dependency sweep | cmp-S09 |
| the order is forced and dbarts merges first | chg-C02, chg-C03, cmp-K09, changes register, merge order preamble |
| consumer branches, versions and what each calls | changes register, consumer readiness, the stan4bart, bartCause, treatSens and bairrtt rows |
| the exact-hash flag and its removal | chg-C23, cmp-U16 |
| obligations at the first step | cmp-U15, cmp-U16 |
| the six-step order | changes register, merge order, steps 1 through 6 |
| the loud and silent rows | changes register, merge order, steps 2 through 6 |
| setResponse warns once per session | chg-U30, dec-B06 |
| bairrtt installs, with draws moved by the defaults | changes register, consumer readiness, the bairrtt row, and merge order, step 2 |
| n.cuts as an integer slot, and stan4bart's attribute write | chg-U03, cmp-U40 |
| no consumer verified against the final header | cmp-U15 |
| the seed dropped by a formals-derived filter | dec-A04, changes register, consumer readiness, the stan4bart row |
| count of user-visible changes that return a different answer | changes register, section 1, the rows marked behavioural |
| loud breakage causes | dec-A01, dec-A02, dec-A03, dec-A04, dec-A05 |
| silent breakage causes | dec-A06, dec-A07, dec-A08, dec-B01, dec-B02 |
| the attribution rule and the two counts | dec-A01 through dec-A70, dec-B01 through dec-B117 |
| capability foreclosed | dec-A01, dec-A61, dec-A46 |
| results changed | dec-A06, dec-A07, dec-A08, dec-A04, dec-A02, dec-A05, dec-A16, dec-A14, dec-A10, dec-A13, dec-B04 |
| speed and build time | dec-B90, dec-B113 |
| maintenance surface | dec-A29, dec-A48, dec-A49, dec-A65 |
| maintainer rows to confirm first | dec-B22, dec-B40, dec-B41 |
| the misattributed large-n quote | dec-B50, cmp-D24 |
| the six status counts | cmp-U01 through cmp-U46 and cmp-K01, cmp-X01 through cmp-X04, cmp-L01 through cmp-L07, cmp-D01 through cmp-D24, cmp-V01 through cmp-V07, cmp-S01 through cmp-S08 |
| unfinished items touching the merge | cmp-K01, cmp-U01, cmp-U02, cmp-U15, cmp-U41, cmp-U44, cmp-U45 |
| abandoned items, and the decision behind the fourth | cmp-X01, cmp-X02, cmp-X03, cmp-X04, dec-A36 |
| what is blocked and on what | cmp-L01, cmp-L02, cmp-L03, cmp-L04, cmp-L05, cmp-L06 |
| deferred without maintainer evidence | cmp-D10, cmp-D24 |
| deferred by the maintainer | cmp-V02, cmp-V03 |
| the stale gates | cmp-S01, cmp-S04, cmp-S07, cmp-S09 |
| the gates reconfirmed at the branch tip | cmp-S02, cmp-S03, cmp-S05, cmp-S06 |
| the cross-engine statistical record | dec-A61 |
| claims that hold at full strength | cmp-K07, cmp-K08 |
| the two doors and the shim between them | chg-U66, chg-U67, chg-U68, chg-U69, dec-B75, dec-B83 |
| a removed name reaches a marker rather than an error | chg-U71, dec-B76 |
| the exact-hash flag is already out of all five consumer files | cmp-U16, chg-C32, dec-B111 |
| the two named items inside the sister-verification step | cmp-U49, cmp-U50 |
| the four engine settings, and their control-only resolution | chg-U84, chg-U85, chg-U96, cmp-U47, dec-B91, dec-B116 |
| the fused pass, and the size gate it owes | chg-U94, cmp-D26, dec-B113 |
| within-chain threading closed, and the thread default capped | chg-U83, cmp-K21, dec-B115 |
| the setPredictor partition, its cost and the re-recorded baselines | chg-U97, chg-I21, cmp-K36, cmp-D25, dec-B117 |
