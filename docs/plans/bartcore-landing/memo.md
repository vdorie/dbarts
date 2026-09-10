# Landing memo: bartcore

## 1. Merge feasibility and order

main carries no commit the branch lacks. The merge brings exactly the branch's
own changes, and nothing has to be reconciled.

One gate stands in front of the merge itself, and it is the maintainer's own:
reading the branch review and declaring a release candidate. The CRAN
submission carries three further prerequisites. Two were recorded when the
grouped random-effects fit was removed. The third is a CRAN package, lorax,
whose example fits bart() with a three-level factor response. 0.9-x silently
coded the levels as numbers and 1.0-0 refuses, so its author has to adjust the
example.

The order is forced. This branch ships a single C header that dbarts main does
not carry, and two of the three compatibility branches already include it. All
four consumer version floors name a dbarts version that no branch of dbarts
satisfies until the merge. dbarts merges first, and the only free variable is
the consumer mains' lag.

stan4bart has a compatibility branch at 0.0-14, with its sampler back end
replaced wholesale. bartCause has one at 1.0-10, the same version string its
main carries. treatSens has one at 3.0-1, and CRAN carries no copy of it.
bartCause and bairrtt call dbarts only from R, and bairrtt needs no branch.

The exact-hash flag below is a compile-time switch. It folds the entry-point
signatures, enumerations, callback and struct layouts the header declares into
one token, baked when the consumer is built. The installed package reports its
own token when the consumer first resolves a stub, and a mismatch is refused
loudly.

| step | action |
|---|---|
| 1 | Pre-merge, pushing to no main: re-bake the header hash, reinstall and test each consumer, push compatibility branches only |
| 2 | Merge bartcore into dbarts main and push. The GitHub window opens |
| 3 | In the same sitting, merge each consumer's compatibility branch into its own main. The GitHub window closes |
| 4 | Point the consumer-build check at each consumer's main, take the consumer workflows off their branch pins and the branch notes out of the readmes, keep the temporary install steps |
| 5 | Clear the three release prerequisites, then submit dbarts 1.0-0 to CRAN. On acceptance, submit stan4bart 0.0-14 and bartCause the same day. The CRAN window opens |
| 6 | When CRAN carries all of them, delete the temporary install and reinstall steps |

CRAN serves dbarts 0.9-34 until the submission is accepted. A consumer's own
resolver would install that older copy, so the temporary install steps have to
survive the pin flip.

The first step carries one obligation, and it is not discharged. Every
consumer has to be re-verified against the final re-baked header. The checklist
names two items inside that step: bartCause hands out a thread budget the new
default would warn about, and treatSens reaches two dbarts internals by name.
Each consumer has been run against an intermediate header, none against the
final one.

The exact-hash flag is already out of all five consumer files, so the version
pair is the only guard. A consumer built against a different header sharing
that pair is admitted. A check on dbarts's own side replaces the flag: a build
fails when the baked token moves and the version pair does not.

## 2. Installer breakage in the intermediate states

A window is a span in which an installer gets a broken pairing. Loud means a
compile or load error. Silent means the call succeeds and returns a different
answer. Four rows carry neither, and each says what it does instead.

| state | who breaks | loud or silent | what happens |
|---|---|---|---|
| after the dbarts merge, GitHub | a 0.9-x user script | loud | rbart_vi, xbart's control argument, a three-element burn-in, a factor response of three or more levels and fractional counts error |
| after the dbarts merge, GitHub | a 0.9-x user script | silent | the proposal mixture, the binary k prior, categorical factors, missing-data incorporation and combined chains change the answer |
| after the dbarts merge, GitHub | a 0.9-x user script | silent | a fitting call with two or three unnamed arguments takes the new front door, at the new defaults, rather than the old function |
| after the dbarts merge, GitHub | a 0.9-x user script | warns once per session | a second positional argument to setResponse means updateScale; a renamed argument, a retired thread method and a call in the old vocabulary each name their successor and then do the right thing |
| after the dbarts merge, GitHub | stan4bart main, treatSens main | loud | a source install fails to compile against a deleted header |
| after the dbarts merge, GitHub | bartCause main | loud | the grouped route and the cross-validation call error |
| after the dbarts merge, GitHub | bartCause main | silent | it installs, records one chain, and returns draws moved by the default changes |
| after the dbarts merge, GitHub | bairrtt main | silent | it installs and every signature it calls is unchanged, and its draws move with the defaults |
| after the dbarts merge, a CRAN copy | stan4bart 0.0-13 | loud | it fails at load, R reporting a function not provided |
| after the dbarts merge, a CRAN copy | bartCause 1.0-10 | mixed | loud on the grouped route and cross-validation, silently one chain otherwise |
| after the dbarts merge, a local copy | treatSens 3.0 | loud | it fails at the first analysis call, its lookups being lazy |
| after the consumer merges, GitHub | unverified | unknown | no consumer has been built against the final header, and stan4bart's branch carries an open defect |
| after CRAN accepts dbarts 1.0-0 | stan4bart 0.0-13, bartCause 1.0-10 | as above | a package update alone breaks them until the consumer submissions are accepted |
| after CRAN accepts dbarts 1.0-0 | lorax 0.1.0 | loud | its example and tests error on a factor response until its author changes them |
| after CRAN accepts dbarts 1.0-0 | insight 1.5.4, off CRAN only | loud | its test expects NA-response rows dropped and gets an error; the test is skipped on CRAN |
| after CRAN accepts the consumers | nobody | neither | closed |

The open defect is stan4bart writing n.cuts through an attribute on the control
object. n.cuts is an integer S4 slot here, so a caller who passes a double
lands one in that slot and the next assignment refuses it.

One break survives both windows. stan4bart picks its control arguments by
matching the control constructor's formals, so the renamed seed is dropped and
the engine seed the caller passed is lost. The filter is unchanged on the
compatibility branch.

User-visible changes that return a different answer: 36

The register lists every such change.

Neither window can be eliminated, only shortened. The GitHub window is as short
as one sitting makes it. The CRAN window cannot be closed, because neither
published consumer bounds dbarts above and CRAN accepts one submission at a
time.

The removals and renames behind the loud rows are agent-made, and this breakage
is their cost. So are categorical factors, combined chains and the binary k
prior. Missing-data incorporation and the proposal mixture carry maintainer
quotes.

## 3. Agent-made decisions and their costs

The registers behind this memo are three files of tables, one row per decision,
change or open item. A decision is the maintainer's only where a row quotes the
maintainer choosing, or lays out alternatives and names the pick. An approval
marker, a standing grant, or no record at all leaves it agent-made.

Agent-made decisions carrying an identified cost: 70

Decisions with maintainer evidence: 115

- **Grouped random effects** (A01): rbart_vi and its S3 methods are deleted,
  so bartCause main breaks at load. The removal's own gate compared the
  replacement's mixing against the deleted path and failed. The maintainer has
  since confirmed it on a condition (B105): stan4bart clears its mixing bar
  first, or an R-only grouped intercept stands in.
- **Equivalence evidence** (A61): the bitwise check compares only against a
  baseline this branch recorded. The old engine is deleted, so no further
  baseline can be recorded.
- **Fixed engine constants** (A46): the leaf regression column cap, the
  perturb width and the cut width are reachable from no argument, and above the
  categorical enumeration limit a split proposal cuts the sorted level list
  rather than choosing a subset. Four constants became settings on the
  maintainer's ruling (B91). Whether they belong on the fitting function instead
  is an open agent-made fork.
- **Empty-leaf veto law** (A12): a leaf holding no positively weighted row is
  vetoed by rank rather than penalised. A vetoed chain mixes
  at constant likelihood, and the draw law differs from 0.9-x.
- **Factor predictors** (A06): the modern front door, the data constructor and
  cross-validation default to single categorical or ordinal columns, while the
  legacy door keeps indicator expansion (B78). The two doors fit different
  models on the same data by default, and the manual says so.
- **Binary node hyperprior** (A07): the default moves, and so does the scale
  default of chi(), the node hyperprior constructor. Every probit fit's posterior
  moves. The maintainer keeps the new default for now and has scheduled the
  evidence it rests on (B106). No such harness is in the tree.
- **Combined chains** (A08): the modern front door defaults combineChains to
  TRUE and flattens chain-major, so code indexing the chain margin breaks. The
  recorded decision was a default of FALSE family-wide.
- **RNG surface** (A04): rngKind and rngNormalKind are dropped and rngSeed
  becomes seed, the old spelling accepted with a warning for one release. Only
  Mersenne-Twister is reachable, so no fit can match another package's stream.
- **No deprecation cycle** (A02): the rule is replaced. The maintainer's own
  removal rule (B76) keeps a removed or renamed name reachable for one release
  as a marker naming its successor, and 31 of them ship. rbart_vi still errors,
  and the old C++ interface vanishes without one.
- **Cross-validation redesign** (A05): xbart drops its control argument and
  replaces its threaded fold loop with a cluster. A consumer passing a prebuilt
  control breaks. The maintainer kept the redesign and required a deterministic
  stream per grid cell (B77), so results no longer depend on the thread count.
- **Family detection** (A16): family = "auto" reads a factor response and
  routes it to probit, ordinal or multinomial. 0.9-x fit gaussian on the same
  level codes.
- **Silent state and a silent cap** (A14, A10, A13): mutators store state only
  on an explicit updateState = TRUE, so a later save writes stale state. The cap
  on the sampled k is silent, while the cap on initial-forest draws raises an
  error naming the attempt count.
- **Weighted sigma** (A11): the residual-variance posterior degrees of freedom
  count rows with positive weight. A fit carrying zero-weight rows moves against
  0.9-x, and one with all-positive weights does not.
- **Chains-only parallelism** (A37): a single-chain run gets no sampling
  parallelism at any thread count. Within-chain threading was re-measured and
  closed (B115). The thread budget keeps its own meaning, its default is capped
  at the chain count, and a larger budget warns rather than being ignored.
- **Worker run latency** (A38): fixed on the maintainer's word (B88). The run
  loop returns when the last chain finishes rather than on a fixed tick, and a
  four-chain single-sweep call is a thousand times faster.
- **Scalar draw path** (A42): the toggle the maintainer asked for (B73) is
  built, as a development build flag carrying the bitwise gates (B90). The
  vector kernels were then measured at under one percent and did not ship
  (B113). A fused pass over the weighted families landed instead, and moves
  those families' draws in the last bits.
- **The published C boundary** (A29): the header exposes a flat set of calls
  over an opaque handle, with families selected by string-named attributes on
  the R objects. On the maintainer's ruling (B84) it is pure C and creates
  nothing: a consumer builds the sampler through R. A non-R host cannot create
  one at all, and the plan for that is post-release (B85).
- **Error reporting** (A31): there are no error return codes, only R's error
  longjmp. A caller must sit in a frame safe to unwind, and sites in the C API
  leak across it.
- **Engine dependence on R** (A50): the model header calls Rmath density
  functions, and the sampler reaches R's print and error entry points. The
  TODO file's Python-binding entry asserts the opposite.
- **Install floor and platforms** (A48, A49): the R floor rises from 3.1-0 to
  4.2.0, and Solaris and big-endian support are removed. Every R release below
  4.2.0 that worked under 0.9-x is excluded.
- **Kernels at weight zero** (A47): two move kernels compile into every build
  at a default weight of zero, so no default fit reaches them. One moves a
  node's cut point, and one is the rule draw at a node whose children are
  leaves.
- **Standing discretion grant** (A65): a standing permission let user-facing
  forks be settled without reaching the maintainer. A class of them therefore
  carries no maintainer attribution.

The register marks the rows below as the maintainer's, whose evidence covers
less than the row claims. Confirm them first.

- **Gaussian-process leaves** (B22): the ruling keeps them and flags a
  deliberate decision on whether they belong in the release. That fork is open,
  so an engine component and a man page freeze in. What has been settled is
  narrower (B110): a fit that falls back to ordinary leaves now warns.
- **Release-review rulings** (B40, B41): both are recorded as rulings in a
  batch, with no quoted words. One refuses a mid-chain rescale of a
  heteroscedastic fit; the other is superseded by the fallback warning above.
- **Enabling value as a gate** (B48): the quoted standing rule licenses public
  entry points that freeze at release. It covers the uncalled C entries and the
  broad R argument list.
- **The adoption slate** (B53): the quote settles that the R-versus-C++ slate lands and how it is framed. Its arcs add public entry points that freeze at the release, and the helpers' export shape is agent-made (A54).
- **Single-chain workloads** (B50): the quoted standing fact keeps large-n
  items live. The vectorized draw-path loops it pointed at were measured and
  declined. What serves it is the fused pass and the memory reductions.
- **Veto-rate thresholds** (B66): the acceptance bands were ratified on a
  recommendation, over a recorded caveat that they were judged and not measured.
  They ran once, and their harness is not in the tree.

The register holds the rest.

## 4. Unfinished, abandoned and stale-gate work

Unfinished, with no decision recorded: 50

Abandoned: 4

Blocked: 8

Deferred with no maintainer evidence: 26

Deferred by the maintainer: 9

Done, but proven only by a gate that predates the tip: 9

Some unfinished items touch the merge. Sites in the C API hold C++ objects that
own heap memory across R's error longjmp, and leak. None is reachable by an
ordinary refusal. The memory-check gate runs the suite under a leak detector.
It failed on its only run and has not run since, so the merge is its first real
run. It has never covered the two exported data-augmentation helpers, which are
new C code. bartCause has no per- push check of its own and carries the silent
failure mode. treatSens has no checks of its own, only the monthly reverse-
dependency smoke test run from here.

Three abandoned items carry no decision. A non-conjugate move strategy for
Gaussian-process leaves was designed and never built. A causal-forest
treatment-scale rescale was derived and never implemented. A slate of post-
release flat-C readers was written down once. The TODO file carries an entry
for none of them. The fourth item, a cross-repository consumer build under the
sanitizers, was added and deleted the same day under an agent-made decision.

The maintainer's own read of the branch and the release-candidate declaration
block the merge. They also block registering the five scheduled workflows for
memory checks, equivalence, calibration, protect balance and the reverse-
dependency smoke test, whose triggers cannot fire from a non-default branch. A
nonzero share for the leaf-pair rule draw is settled the other way: it clears
the sample-size bar and regresses coverage past the margin, and the maintainer
has ruled that it stays at zero through the merge, its weight to be set with a
rationale by the mixing research that resumes afterward. A bartCause
route for grouped data is blocked because the deleted function was its only
path. stan4bart owes a mixing bar on the random-effect scale, which the same
removal named as a release prerequisite.

The deferred rows without maintainer evidence include a Python binding whose
stated premise the code contradicts, and a size threshold below which the new
fused per-sweep pass should be declined, which cannot land until the bitwise
corpus carries a fit large enough to exercise it. The rows deferred by the
maintainer include real-valued negative-
binomial dispersion and arbitrary real binary weights. Both wait on an
approximate data-augmentation draw the maintainer declined.

The stale gates prove their items only at a commit before the branch tip. The
protect-balance check, which verifies that R objects stay protected from
garbage collection, failed, and its fix carries a local re-run alone. The
comparisons against recorded baselines have since re-run and were green, three
engine commits before the tip. The check that the sampler recovers parameters
drawn from its own prior sits earlier still. Every slice records its own run of
the test suite, the package check, the C++ tests and the sanitizer builds, so
no engine change is unexercised.

The CRAN reverse-dependency sweep is one of the stale gates. It ran in July at
an earlier engine, with no script in the tree, and found the two breaks above.
Three packages have declared dbarts since and have never been checked.

One record compares this engine against the old one. It is statistical rather
than bitwise, and was taken at the engine cutover. Its manifest marks it
evidence and not a gate.

Scenarios it covers: 9

It has never been re-run against the current engine. No standing gate compares
the two engines, so a drift inside the statistical band would be caught
nowhere.

Every entry point in the shipped C header is exercised by a compiled consumer
test. The version and hash mismatch handshake is itself tested.

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
| the attribution rule and the two counts | dec-A01 through dec-A69, dec-B01 through dec-B73 |
| capability foreclosed | dec-A01, dec-A61, dec-A46, dec-A12 |
| results changed | dec-A06, dec-A07, dec-A08, dec-A04, dec-A02, dec-A05, dec-A16, dec-A14, dec-A10, dec-A13, dec-A11, dec-B04 |
| speed and build time | dec-A37, dec-A38, dec-A42, dec-B73 |
| maintenance surface | dec-A29, dec-A31, dec-A50, dec-A48, dec-A49, dec-A47, dec-A65 |
| maintainer rows to confirm first | dec-B22, dec-B40, dec-B41, dec-B48, dec-B53, dec-B50, dec-B66, dec-A54 |
| the six status counts | cmp-U01 through cmp-U46 and cmp-K01, cmp-X01 through cmp-X04, cmp-L01 through cmp-L07, cmp-D01 through cmp-D24, cmp-V01 through cmp-V07, cmp-S01 through cmp-S08 |
| unfinished items touching the merge | cmp-K01, cmp-U01, cmp-U02, cmp-U15, cmp-U41, cmp-U44, cmp-U45 |
| abandoned items, and the decision behind the fourth | cmp-X01, cmp-X02, cmp-X03, cmp-X04, dec-A36 |
| what is blocked and on what | cmp-L01, cmp-L02, cmp-L03, cmp-L04, cmp-L05, cmp-L06 |
| deferred without maintainer evidence | cmp-D10, cmp-D24 |
| deferred by the maintainer | cmp-V02, cmp-V03 |
| the stale gates | cmp-S01, cmp-S02, cmp-S03, cmp-S04, cmp-S05, cmp-S07, cmp-S09 |
| the cross-engine statistical record | dec-A61 |
| claims that hold at full strength | cmp-K07, cmp-K08 |
| the two doors and the shim between them | chg-U66, chg-U67, chg-U68, chg-U69, dec-B75, dec-B83 |
| a removed name reaches a marker rather than an error | chg-U71, dec-B76 |
| the exact-hash flag is already out of all five consumer files | cmp-U16, chg-C32, dec-B111 |
| the two named items inside the sister-verification step | cmp-U49, cmp-U50 |
| the four engine settings, and the open fork above them | chg-U84, chg-U85, cmp-U47, dec-B91 |
| the memory reductions, and the copy that is left | chg-U92, chg-I15, cmp-D27 |
| the fused pass, and the size gate it owes | chg-U94, cmp-D26, dec-B113 |
| within-chain threading closed, and the thread default capped | chg-U83, cmp-K21, dec-B115 |
