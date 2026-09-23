# Decision register

This register lists every design decision the 1.0-0 work carries forward, so
the maintainer can see what was decided on their behalf. Section A holds
decisions the agents made that carry a cost and need the maintainer's
adjudication. Section B holds decisions with maintainer evidence behind them.
Section C holds agent-made decisions with no identified cost. Section A runs
roughly from the costliest entry down, with the ones a later sweep added at the
end. A decision counts as the maintainer's only where the record quotes them
choosing, or lays out alternatives and records their pick; everything else,
including approval after the fact and silence, is tacit and counts as
agent-made. The maintainer marks each entry mine, not mine, or revisit on its
Marked line. Later documents cite an entry's id rather than saying the
maintainer decided.

## A. Agent-made decisions with a cost

**Grouped random effects leave the package**
The package no longer fits BART with grouped random intercepts: rbart_vi and the S3 methods registered alongside it are gone, and that structure now lives in stan4bart. No alternative was weighed when it was removed; the later ruling considered an R-only grouped intercept driven through the sampler's offset each sweep and rejected it. A user of 0.9-x who called rbart_vi has no runnable path in the package, bartCause fails at load until it is rebuilt against stan4bart, and the check on the group standard deviation's mixing that was meant to gate the removal failed on two of three seeds while the removal was kept. The maintainer has ruled that the function stays removed but as a tombstone, an exported name whose body stops with an error pointing at stan4bart, with one line saying the group-spread prior differs there so results move; shipping is conditioned on stan4bart's own mixing work, which gates the submission to CRAN rather than the merge to main. See also: [dec-B105].
Record: docs/design/retire-grouped-random-effects.md, which quotes no maintainer words on the removal. Marked: blank. [dec-A01]

**How a removed function reaches users**
A function removed before 1.0-0 stays exported for one further release as a tombstone: the name still resolves, and calling it stops with an error naming its successor. The alternative, which the agents had adopted, was that nothing is released yet so anything may be removed outright with no warning cycle, and that rule is what took out rbart_vi, the sampler's two thread methods, the cross-validation function's control argument and the whole old compiled C++ interface that other packages linked against. A user upgrading a 0.9-x script gets an error that says where to go instead of an object-not-found error. The maintainer ruled the tombstone rule in, kept the thread methods as no-ops, and left a renamed argument to the implementers' judgement, with the old rngSeed name accepted for one release with a warning. See also: [dec-B76].
Record: docs/plans/archive/bcf-public-surface.md, docs/plans/archive/multiforest-extension-surface.md. Marked: blank. [dec-A02]

**The sampler's thread methods stay as no-ops**
The sampler object keeps startThreads and stopThreads, which do nothing. The agents' alternative, not taken, was to delete both. A user calling either sees nothing happen; in the released package they drove the lifetime of the threading machinery, and in this work they had already become no-ops, so deleting them would have removed two methods documented in the manual and announced in the 0.9-31 news with no news entry of their own. The maintainer ruled they stay. See also: [dec-B76].
Record: docs/plans/archive/r5-cleanup.md; docs/design/public-surface.md. Marked: blank. [dec-A03]

**The sampler's generator options shrink**
The sampler takes one seed argument, spelled seed, and the engine offers only R's Mersenne-Twister or a user-supplied uniform generator, with inversion or a user-supplied normal; the rngKind and rngNormalKind arguments are gone. No alternative was weighed. A user who passed rngSeed by name now passes seed, and a package that builds its dbarts call from the formals silently drops the seed, so the engine's stream is no longer the caller's while the caller's own R seed still applies; the released stan4bart 0.0-13 is affected until its next release, its development branch remapping the old name with a once-per-session warning and stopping on an unmatched argument name. No argument selects a generator any more, so a caller cannot make dbarts follow another package's stream: the engine keeps its user-supplied uniform and normal arms in the support library, but nothing in R reaches them and the bridge always builds Mersenne-Twister chains. The values of the C enumeration are renumbered. Not yet ruled on.
Record: docs/plans/archive/bart2-argument-consolidation.md, under a naming-scheme heading marked decided with no maintainer quote. Marked: not mine. [dec-A04]

**Cross-validation gets a flat signature**
The cross-validation entry point takes its settings as ordinary arguments rather than a prebuilt control object, runs its folds on a cluster from the parallel package instead of a C++ thread loop, and starts each fold and replication fresh rather than carrying a chain over. No alternative was weighed when this was adopted; the maintainer was later given three. A caller that passed a prebuilt control breaks, the design matrix is serialized to every worker, per-cell cost rises because each fold burns in again, and results depend on the thread count even with a seed fixed. The maintainer kept the redesign and ruled that each grid cell gets its own deterministic stream so a seed reproduces at any thread count, that the removed arguments get tombstones, and that warm starts across folds return only if timing against 0.9-34 is bad, and then opt-in with the leakage stated; a later ruling on the entry points' signatures gives the function a control argument again. See also: [dec-B77], [dec-B116].
Record: docs/design/public-surface.md, which decides to break freely now and patch at the end, with no maintainer quote. Marked: blank. [dec-A05]

**Factors become one column split by subset**
A factor predictor enters the design matrix as a single column that trees split by choosing a subset of its levels, on the modern entry points, while the BayesTree-compatible function keeps the old indicator expansion. The maintainer was given three options and chose subset splits as the modern default. Every fit with a factor predictor differs from 0.9-x: the number of columns, the prior over splits and the names in the variable counts all change, variable importance is counted per factor rather than per level, and the two front ends fit different models on the same data by default. The maintainer ruled this in, with the manual stating the difference from BayesTree and why, and left weighting a factor column's split mass by its level count as post-release research. See also: [dec-B78].
Record: docs/design/public-surface.md, which states the major-version license for the change with no maintainer quote, and docs/plans/archive/interface-review.md, which records the flip as a side effect of the engine cutover, before the review meant to decide it. Marked: blank. [dec-A06]

**The binary leaf-scale prior gets a new default**
For binary outcomes the leaf-scale parameter k takes a chi hyperprior with 1.5 degrees of freedom and scale 2, and the chi constructor's own scale default moves from infinite to 2. Every probit fit's posterior moves against 0.9-x, and an explicit call giving only the degrees of freedom moves too. A study in July rejected the infinite scale clearly and preferred scale 2 over 5 narrowly, saying nothing about the degrees of freedom, and the relabel of the constructor's first argument is the maintainer's. The alternative not taken was any other arm of a wider grid: an extended study, run on converged chains with the scale grid extended below 1 and weighed against real as well as simulated data, found the coverage criterion has an interior optimum sitting where chi(1.5, 2) already is, with the largest available improvement anywhere on the grid costing more than it gains in log score or worst-cell regret. The maintainer confirmed chi(1.5, 2) as the shipped default on 2026-09-14, on that evidence (docs/plans/binary-hyperprior.md). A revisit stays scheduled for after the mixing research, since the sampled k's own non-convergence leaves the optimum's location resolved only to about a factor of two in scale. See also: [dec-B106], [dec-B118].
Record: docs/plans/archive/chi-hyperprior-df.md for the relabel and docs/plans/archive/chi-default-research.md for the scale default, the second with no maintainer quote; docs/plans/binary-hyperprior.md for the confirming study. Marked: mine. [dec-A07]

**bart2 returns its chains already combined**
bart2 combines chains by default, and the combined draws are laid out chain by chain rather than sample by sample. The recorded alternative, and the earlier decision, was to keep chains separate across the whole family of entry points. Code that indexes the chain margin of a bart2 result breaks, the order in which the combined per-sample scalars line up changes for everyone, and bart2 now disagrees with the sampler's own run method about the default. The only maintainer statement on the row settles that such changes belong in this release, not which way this one goes. The maintainer has claimed the decision as their own.
Record: docs/plans/archive/interface-review.md for the timing and commit 5287837a for the layout order. Marked: mine. [dec-A08]

**Ordered factors get cuts between level codes**
An ordered factor predictor is split at midpoints between its level codes, giving one fewer cut than it has levels, and the per-column cut cap is raised to fit them rather than the grid being thinned. No alternative was weighed. The n.cuts argument no longer has any effect on an ordered factor, draws move for any design containing one, and a data object saved under 0.9-x has no upgrade path. Not yet ruled on.
Record: docs/design/public-surface.md states the grid, and docs/plans/column-kind-consolidation.md marks the consolidation's rulings as taken at the implementers' discretion under a standing grant. Marked: not mine. [dec-A09]

**The first forest is rejection sampled**
The sampler draws its starting forest by sampling whole trees over and over until no leaf is empty, up to ten thousand attempts, instead of drawing once from the prior and collapsing the empty nodes. No alternative was weighed. The law governing the state at sweep zero differs from every released dbarts, and the attempt cap is not silent: exhausting it raises an error naming the count, so the one silent cap in the engine is the clamp on the drawn leaf scale. Not yet ruled on. See also: [dec-A13].
Record: code only, in the chain header's prior-draw path. Marked: not mine. [dec-A10]

**Zero-weight rows leave the degrees of freedom**
The posterior for the residual variance counts only rows with a positive weight in its degrees of freedom, where the released package added the full observation count. No alternative was weighed. A weighted fit that carries zero-weight rows gets a different variance posterior than under 0.9-x; when every weight is positive the two counts agree and nothing moves. Not yet ruled on.
Record: docs/plans/archive/sigma-df-zero-weights.md, with no maintainer quote. Marked: not mine. [dec-A11]

**Empty leaves are vetoed, not penalized**
A proposed tree move that would leave a leaf with no positive-weight observations is refused outright: proposals are ranked first on how many such leaves they leave and only then on likelihood, rather than being charged a finite penalty. No alternative is recorded as weighed, and the maintainer's stamp on the row covers when the fix lands, not its shape. A chain in a vetoed state moves at constant likelihood, driven by the prior and the proposal alone, and the draw law differs from 0.9-x. The maintainer has claimed the decision as their own.
Record: docs/design/r-c-division.md, whose adoption slate carries a maintainer stamp, and docs/design/empty-leaf-veto.md. Marked: mine. [dec-A12]

**The drawn leaf scale is capped silently**
When the leaf-scale parameter k is given a hyperprior and drawn each sweep, the draw is clamped at one million, with no warning and no entry in the news file. No alternative was weighed, although the question was routed to the maintainer. A user sees nothing: the clamp changes behaviour only in the runaway regime it exists for, and it is undocumented. Not yet ruled on.
Record: docs/plans/archive/chi-k-runaway.md, which records a cap on the sampled k alone with no warning path and no maintainer quote. Marked: not mine. [dec-A13]

**Setters store state only when told**
The sampler's mutating methods write their state back to the R object only when the caller passes updateState = TRUE; leaving it missing no longer falls back to the setting on the control object, which the run method and the prior-sampling and grow-from-root methods still honour. No alternative was weighed. One object now carries two conventions, and saving a sampler after a mutation silently writes stale state. Not yet ruled on.
Record: docs/plans/archive/interface-review.md, which records that no ambiguity surfaced needing the maintainer. Marked: not mine. [dec-A14]

**fitted gains a confidence level argument**
The fitted method for a bart fit takes ci.level as its third positional argument. No alternative was weighed. A 0.9-x call that passed a train-or-test string in the third position errors rather than quietly meaning something else, because the interval helper matches that string against its own argument names and refuses it; the pre-release surface freeze rules on the predict methods and does not name this one. Not yet ruled on.
Record: docs/plans/prerc-surface-freeze.md and commit 7b3ac6bf. Marked: not mine. [dec-A15]

**The response type is detected automatically**
Asking for family "auto", the default, inspects the response: two levels go to probit, an ordered factor to the ordinal family, three or more unordered levels to multinomial, and the call prints one line saying what it chose. No alternative was weighed. The model is therefore chosen from the data by default, where 0.9-x would have fitted a gaussian model to the same level codes. Not yet ruled on.
Record: docs/plans/archive/package-review-remediation.md, whose decisions are adopted under a standing discretion grant the maintainer may overrule. Marked: not mine. [dec-A16]

**The fit object's component set varies**
A fit object drops any component whose value is empty, so which names it carries depends on the family, the number of forests, the residual law and the options given. No alternative was weighed; the rule is decided nowhere. A user cannot rely on the names of a fit, and downstream code and the package's own methods must guard every access. Not yet ruled on.
Record: docs/plans/archive/bart2-argument-consolidation.md tabulates the packaging; the rule itself is recorded nowhere. Marked: not mine. [dec-A17]

**Two output changes ride in the bridge**
Two behaviour changes came in with the C bridge and have nothing to do with the interface rewrite: the model-matrix builder writes NA cells for a missing factor code, where the released package wrote through an out-of-bounds index, and the core-count probe returns NA rather than -1 on its generic Unix fallback. No alternative was weighed. A user building a model matrix from data with missing factor codes gets NA instead of whatever the stray write produced, and a user on a platform without a native core count sees NA; the matrix case has a test, the core-count case does not. Not yet ruled on.
Record: code only, in the model-matrix builder and the core-count probe. Marked: not mine. [dec-A18]

**Hurdle models are composed in R**
A hurdle model is fitted as two independent samplers composed in R, with no engine support, and one token, hurdle.lognormal, names it on the modern entry points. The alternative, which shipped first, advertised both hurdle.lognormal and twopart on the family argument while both always errored, kept only so that forwarding from bart2 would not choke; the maintainer was given three options and took the recommendation. A user who asked for either token got an error from a documented value; now the entry point intercepts the hurdle case before forwarding, and a caller of twopart gets a tombstone for one release. A correlated hurdle, with shared trees or correlated leaves, is a post-release engine family. See also: [dec-B82].
Record: docs/design/hurdle.md recommends the R composition, and docs/plans/review-2026-08-24/matrix-review-entries.md leaves the erroring tokens open as a maintainer judgement. Marked: blank. [dec-A19]

**One message covers six refusals**
The C entry that built a sampler could fail only by returning a null pointer, so six distinct reasons for refusing reached the user as one message listing all of them, growing with every family added. No alternative was weighed. A user handed that message had to work out which of the six applied. The maintainer ruled the creation entry out of the shipped C header before the release, so the refusal channel goes with it and a compiled consumer builds a sampler through the R interface instead. See also: [dec-B84].
Record: docs/plans/review-2026-08-24/consolidated-report.md, which leaves the refusal shape open as a maintainer item. Marked: blank. [dec-A20]

**Survival ships without a formula interface**
The discrete-time hazard family shipped with no left-hand side on the formula interface, no subset argument and no test set, so the DESCRIPTION advertised survival while the package's primary interface refused it. The maintainer was given three options and took the recommendation: the survival interface is completed before the release, with a Surv left-hand side, subset honoured for both survival families, and a test path that expands held-out subjects, all in the R layer with no engine change. A user can then fit a survival model the way they fit any other. See also: [dec-B97].
Record: docs/design/survival.md records three forks put to the maintainer, none of which covers the missing formula and subset surface. Marked: blank. [dec-A21]

**Hazard-only arguments move onto the family**
Settings only the discrete-time hazard expander reads, the interval breaks and a cap on expanded rows fixed at ten million, ride on family objects such as hazard(breaks, max.rows), in base R's glm idiom. The alternative, which shipped first, left both as top-level arguments on the two main entry points, inert for every family but one, with the cap unexplained. A user passing either at the top level gets a tombstone for one release, and the consolidation of the entry points' arguments is reopened so that any remaining family-only or feature-only argument moves onto its family object before the push to CRAN. The maintainer ruled for family objects and sent the ten-million cap to the audit of fixed constants. See also: [dec-B98].
Record: docs/plans/archive/bart2-argument-consolidation.md inventories both as inert and in scope; the placement is decided nowhere. Marked: blank. [dec-A22]

**The level step defaults to automatic**
An extra Gibbs step that redraws a forest's overall level runs automatically, but only where that forest's tree-structure proposals are frozen. No alternative was weighed. No shipped default reaches the step, yet it costs a control slot, an argument on bart2, a slot and a validity clause on the S4 object, a parse in the bridge and a manual entry; a user who sets the structural move weights to zero silently gains a Gibbs step. The maintainer has claimed the decision as their own.
Record: docs/design/level-fibre.md, which records the step as kept and automatic under a frozen mixture, naming no ruler. Marked: mine. [dec-A23]

**Diagnostics stop depending on another package**
Split R-hat and effective sample size are computed inside the package, and draws come out through a base-R extractor returning an iterations by chains by variables array with dimnames, which the posterior package's own constructors accept if a user has it. The alternative, which shipped first, suggested the posterior package and let the summary and draws methods return a different column set depending on whether it was installed, so the shape of a returned object depended on the machine. Matrix and survival, both recommended packages that ship with R, stay suggested behind their features, so sparse and survival work is unavailable without them. The maintainer ruled the posterior dependency out entirely, at the cost of a short diagnostics implementation of our own. See also: [dec-B99].
Record: docs/plans/archive/convergence-diagnostics.md states the suggests-only rule; no record argues suggested against required. Marked: blank. [dec-A24]

**Two vocabularies for the family argument**
There are two lists of family tokens: the one users type at the entry points, resolved at the bridge into the families the engine supports, and the engine's own internal list, with a single mapping table in the manual and the BayesTree-compatible function taking more of the user-facing tokens within reason. The alternative, which shipped first, carried three lists, one on the model object, one on the modern entry points and a third on the compatible function, with hazard tokens remapped before creation, twopart aliasing a composition and the compatible function refusing most of the rest by name, so a user could not form one mental model of family. The mapping table costs a maintenance item per new family. The maintainer ruled for the two lists. See also: [dec-B81].
Record: docs/plans/dbarts-h-freeze.md enumerates the three, and docs/plans/review-2026-08-24/consolidated-report.md parks the question as a maintainer item. Marked: blank. [dec-A25]

**Sparse factors reach the formula interface**
A sparse matrix or sparse factor is assigned into the data frame as a column and named in the formula like any other, with no marker; ingestion detects such columns by their class and lifts them around the model frame, and wide factors are stored sparse by the engine without the user asking. The alternative, which shipped first, refused sparse factor columns in the formula interface, leaving an exported class with its own constructor, methods and manual page reachable only through the x and y interface. The maintainer was given three options, ruled for accepting them in the formula, and asked why a sparse matrix already in a data frame would need marking; it lands before the release, in the R layer only. See also: [dec-B100].
Record: code only, in the mixed-matrix formula path. Marked: blank. [dec-A26]

**Residual laws move onto family objects**
The residual law rides the family object: student(df) is a family, so the resid.dist argument and both residual-law vocabulary bundles are gone from the user-facing signature, and no residual-law bundle is exported. The alternative, which shipped first, defined a residual-law bundle exactly like the prior bundle but did not export it, so a user who found the bundled priors looked for the bundled residual laws and did not find them. The maintainer was given three options and took the recommendation; exporting the bundle was the fallback had the consolidation left the argument in place, and it did not. See also: [dec-B101].
Record: docs/plans/archive/bart2-argument-consolidation.md lists both without noting the asymmetry. Marked: blank. [dec-A27]

**The BayesTree drop-in claim is withdrawn**
The DESCRIPTION says the package provides a BayesTree-compatible interface rather than a drop-in replacement, and the three removed configure arguments stay for one release as stubs that stop with a message naming the removal. The alternative, which shipped first, kept the old sentence while the build flag that made the streams match was deleted and autoconf silently ignored the three arguments, so a caller still passing them lost the BayesTree-matching stream with no message and the claim was no longer checkable in the tree. The maintainer ruled both halves, observing that the value in a drop-in replacement was that BayesTree was slow and unmaintained at the time, which is no longer true. See also: [dec-B107].
Record: code only, in configure.ac and DESCRIPTION. Marked: blank. [dec-A28]

**The published C boundary stops speaking R**
The one header the package ships for compiled consumers is pure C: the four entries that took or returned R objects, for creation, storing and restoring state and reading trees, are removed, and a consumer builds a sampler through the R interface and takes the handle from the external pointer on the sampler object, with R owning its lifetime; state and trees go through the R methods, and an opaque C state blob is added only if a consumer needs to store and restore from C. The alternative, which shipped first, published the R representation itself, one flat header over an opaque handle with R objects at creation, state and trees, and its families and compositions selected only by string-named attributes on the control and model objects, which made the S4 slot layout and those attribute names part of the binary interface even though its hash does not cover them. A compiled consumer still cannot create a sampler without R, construct or inspect engine objects, subclass results or read struct fields, and there is no C++ header and no C++ binary compatibility promise. The maintainer ruled the R-object entries out before the release, at the cost of one more header change: the hash is rebaked, stan4bart edited at three sites and treatSens at two. See also: [dec-B84].
Record: docs/design/public-surface.md, whose decided bullets carry no maintainer quote anywhere in the file. Marked: blank. [dec-A29]

**The old compiled interface retires at once**
The old compiled interface goes away in the same release that brings the new engine, with no transition release carrying both. No alternative was weighed beyond the port of stan4bart. A consumer that has not ported has no dbarts release that works with either surface. The maintainer has claimed the decision as their own.
Record: docs/design/public-surface.md, which resolves that no transition release is needed beyond the stan4bart port, with no maintainer quote. Marked: mine. [dec-A30]

**Errors leave the C entry points cleanly**
A C++ exception raised in a callback or in the engine is caught where it happens, under R's unwind protection so the jump unwinds through the callback's own frame rather than across it, and is rethrown once that frame has returned; it becomes an R error only at the bridge entry point, never as a raw jump out of engine or callback code. The alternative, which shipped first, let a raw R error jump straight out of engine and callback code, which left three sites holding C++ heap owners across the jump and leaking; the same change fixes those three sites. An R error raised at a C entry point still unwinds the caller's frames, so a compiled consumer must call from a frame that is safe to unwind, which the header says. A compiled consumer still may not jump to a saved position of its own from inside a callback, though raising an R error or throwing is safe. Not yet ruled on. See also: [dec-B119].
Record: docs/design/public-surface.md, which documents in the header that errors longjmp. Marked: not mine; superseded by dec-B119. [dec-A31]

**A returned integer means one of two things**
A non-void return from a C entry point is either a value or a status saying whether the capability exists, documented entry by entry. No alternative was weighed. A compiled consumer must read each entry's documentation to know which, since an integer is not uniformly a status; the header's contract block now names two kinds and enumerates two, the third kind having left with the entries that returned it. Not yet ruled on.
Record: the shipped header's contract block and commit 9df0cb50. Marked: not mine. [dec-A32]

**Every C setter copies what it is given**
Each C setter copies the caller's data into buffers the sampler owns and allocates once at creation; the sampler never keeps a caller's pointer and never allocates after creation, and the owned copies are the response, offset, test offset and weight vectors, predictors being re-encoded into cut codes as before. The alternative, which shipped first, was per-entry ownership: some setters kept the caller's pointer and the rest did not, a value written through a kept pointer was not guaranteed to be seen, and the forest-weight setter borrowed where the R bridge copied into its own buffer, so a bug found through one front end might not reproduce through the other. A compiled consumer now reads one sentence of contract instead of a table. The maintainer ruled for copying while asking that the footprint at large sample sizes be audited rather than assumed, since fitting larger datasets matters and the package has had memory trouble before. See also: [dec-B87].
Record: docs/plans/latent-subset-mask.md states the asymmetry as inherited fact and adjudicates neither half. Marked: blank. [dec-A33]

**Validation is partial by design**
The C entry points check struct sizes, declared shapes, family support and capability, but dereference the sampler handle, the output buffers and the required input vectors exactly as handed over. No alternative was weighed. A consumer passing a null, destroyed or short pointer crashes rather than getting an error, and the registry of valid pointers the released package kept is gone from the C path. Not yet ruled on.
Record: the shipped header's contract block; no record in either evidence sweep. Marked: not mine. [dec-A34]

**Sparse predictors densify on mutation**
Replacing predictors on a live sampler expands a sparse source to dense before the engine sees it, while reading in a test set and predicting both consume sparse data as it is. No alternative was weighed. A user swapping in a sparse matrix between sweeps gains validation and a uniform argument but no memory saving, and the asymmetry is documented rather than removed. The maintainer has deferred per-observation sparse mutation to after the release.
Record: docs/plans/archive/data-ownership-5-sparse.md states the mutation surface without adjudicating it. Marked: blank. [dec-A35]

**The cross-package sanitizer job is dropped**
The continuous-integration job that built and ran a downstream consumer against the development package under sanitizers was removed the day it landed, leaving a monthly reverse-dependency smoke test in its place. No alternative was weighed. Nothing before a merge now checks that a downstream consumer still builds and runs clean against development dbarts. Not yet ruled on.
Record: docs/plans/archive/capi-dispatch-table.md, which calls the drop the maintainer's own decision but quotes no words and records no fork put to them, which is why the entry sits here rather than in section B; the commit it names for the landing is the pre-rebase one, and on this branch the landing and the drop are 9fe39856 and 99b356d8. Marked: not mine. [dec-A36]

**Threads run chains, not within a chain**
Sampling is parallel across chains only, with a separate pool for the test fit and a fan-out for prediction. Within-chain threading is the alternative and it does not ship: re-measured on the current engine it reached at best 1.03 times, at four workers and a million observations, and lost about 5 percent at two workers and 15 percent at eight at a hundred thousand, so the maintainer ruled it closed and archived, with the correctness half of the prototype, byte-identical draws across worker counts, banked but not revived. A single-chain run therefore gets no sampling parallelism at any thread count, and the thread count is made honest instead: it keeps its own meaning as a total budget distinct from the number of chains, the control default becomes the smaller of the core count and the chain count, and a budget above the chain count warns once per fit, naming both counts, that tree sampling uses at most one thread per chain and the excess reaches only the test-fit pool and prediction. This sits against a standing fact the agents had recorded, that single-chain runs at a hundred thousand observations or more are common. See also: [dec-B50], [dec-B115].
Record: docs/design/within-chain-threading.md, closed on measurement, with the maintainer framing the workload but not ruling the closure. Marked: blank. [dec-A37]

**The run loop no longer sleeps between checks**
A multi-chain run waits for its chains to finish and wakes as soon as the last one does, keeping a hundred-millisecond timeout only for the interrupt poll and the progress flush. The alternative, which shipped first, blocked the calling thread in hundred-millisecond sleeps until every chain finished, adding up to a hundred milliseconds of latency to every multi-chain run call, including the single-sweep call the embedding pattern makes inside an outer loop. The maintainer ruled the fixed sleep replaced before the release, at the cost of one small engine rebuild and its gates and no change to the draws. See also: [dec-B88].
Record: docs/architecture.md documents the inline path's throttle and not the worker path's sleep. Marked: blank. [dec-A38]

**Two threading mechanisms live side by side**
The package runs threads two ways: the C thread manager in the support library, used for the test-fit fan-out alone, and raw threads with signal masking for everything else. No alternative was weighed. Nothing is visible to a user; the thread manager, its queue and its shared header are maintained for one call site, and the signal-mask call sits in a header-only C++ engine behind a platform guard. Not yet ruled on.
Record: docs/architecture.md calls the pool the one sanctioned exception, with no attribution. Marked: not mine. [dec-A39]

**The engine compiles into two translation units**
The engine is header-only and compiled into two translation units, one of which is the eight-thousand-line bridge. No alternative was weighed, and no record of a translation-unit decision exists. Nothing is visible to a user: touching any header recompiles the whole engine serially, and the bridge object is 14 MB unstripped against a 1.9 MB linked library, with no compile-time or peak-memory measurement on record. The maintainer ruled the shape kept for the release but provisional, with an independent review of what the engine's generic axes should be written for later reading and a measured research item after that. See also: [dec-B92].
Record: code only; no record of a translation-unit decision in either evidence sweep. Marked: blank. [dec-A40]

**Five copies of the sampler stack ship**
Five full instantiations of the sampler, chain and move code are compiled in, with the leaf model chosen at compile time and the response family reached through virtual calls. No alternative was weighed; the maintainer's stamp covers only the reduced-precision arm. Nothing is visible to a user; each new leaf model multiplies the whole stack again in binary size and compile time. The maintainer ruled the shape kept for the release but explicitly provisional, these being abstractions chosen very early from the general shape of the problem, and asked for a deep dive on what the generics should be. See also: [dec-B92].
Record: docs/design/reduced-precision-storage.md carries a maintainer directive for the reduced-precision arm; the other four instantiations have no record. Marked: blank. [dec-A41]

**The draw path stays scalar and fixed-order**
Reductions on the draw path are scalar and summed in a fixed order, so draws are reproducible bit for bit within a host. The recorded alternative is its opposite: an earlier maintainer decision made bit-identity a build-time toggle, with a vectorized fast path by default and the scalar kernel forced as the bit-identical reference, and it was never built. A user gets identical draws from identical seeds. The maintainer ruled for vectorizing with a scalar reference build for development; the kernels were then measured at under one percent of a weighted fit on arm64 and within noise on x86, so they did not ship. The estimate that had motivated the work, 3 to 4 percent of total runtime, was an estimate and the measurement is what stands. The reference-build flag is therefore inert, both builds identical, and the speed went instead into a fused residual-and-sum pass for weighted families, worth 26 to 29 percent at a hundred thousand observations. See also: [dec-B73], [dec-B90], [dec-B113].
Record: docs/plans/x86-simd-plan.md records the toggle decision as the maintainer's and marks it unimplemented; the shipped rule is stated in docs/architecture.md. Marked: blank. [dec-A42]

**Three uncalibrated tuning thresholds**
Three cutoffs decide when work goes parallel or sparse: a test fit below 65536 rows stays serial, prediction below a cell count stays serial, and a column is stored sparse below a density of 0.2. All three shipped fixed in code, uncalibrated and movable by nothing but a unit test. The maintainer ruled that such constants are documented with their origin, measured where measurement is possible so that defaults and recommended values are grounded, and that those that matter become control settings; all three now are, as testFitParallelCutoff, predictParallelCutoff and sparseDensityThreshold on the control object, and the prediction cutoff was recalibrated from ten million cells down to fifty thousand. See also: [dec-B91].
Record: docs/design/threaded-predict.md states the predict cutoff is correct at any value and leaves it open, and docs/design/sparse-columns.md states the density threshold without attribution. Marked: blank. [dec-A43]

**The treatment forest's rescaling move ships off**
In the causal-forest family, which fits separate prognostic and treatment forests, the treatment forest's ridge move, a rescaling that leaves the likelihood unchanged, never runs. It is not switched off by a flag: the bridge turns a forest's ridge on exactly when that forest's amplitude prior carries a positive half-Cauchy scale, and the treatment forest's prior is a fixed variance, so the derivation never enables it. The pair of ridge switches in the engine's own struct belong to a fixture the shipped path does not reach, so nothing in R can turn the move on either. No alternative was weighed. Nothing is visible to a user beyond the move not running; an implemented move stays dark to avoid re-recording a baseline. Not yet ruled on.
Record: code only, in the combiner header and the run bridge; docs/plans/archive/bcf-b-ridge.md describes a different joint rescale that was never built. Marked: not mine. [dec-A44]

**The prediction thread argument does real work**
The thread count on predict drives a real fan-out over chains and draws rather than being removed as an inert argument. No alternatives were recorded. A user asking for threads on a large prediction gets them, at the cost of a new parallel code path and a change to the C interface to carry the count. The maintainer ruled it stays wired, with its serial cutoff calibrated in the audit of fixed constants rather than left at a fixed value. See also: [dec-B93].
Record: docs/design/threaded-predict.md records the maintainer ruling the argument wired to real threading, with no alternatives recorded. Marked: blank. [dec-A45]

**Engine limits fixed in the code**
Several engine constants are fixed in code: cut codes are 16 bits, so a column has at most 65533 cut points and 65535 levels; exact enumeration of categorical partitions runs only up to ten present levels, above which the scan-based proposals see prefix splits rather than subsets; a leaf regression may designate at most eight covariate columns, so a nine-column one is refused; and the perturb move shifts a cut by exactly one grid position, any other width needing a private build. No alternative was weighed. The maintainer ruled that all of them are documented with their origin, measured where possible, and that those that matter, the categorical enumeration limit, the two parallel cutoffs and any the measurements show binding, become control settings, having said they seem arbitrary and may be limiting what can be fitted. The categorical limit is now the control's categoricalExhaustiveCap, still ten by default; the cut-code width, the leaf-regression column limit and the perturb width remain fixed. The maximum leaf size for Gaussian-process leaves is not one of these, being an argument of the gp() constructor. See also: [dec-B91].
Record: code only, in the data, scan, model, tree and moves headers; docs/design/public-surface.md states the opposite for a cut-width build option that was removed. Marked: blank. [dec-A46]

**Two tree moves ship at weight zero**
Two tree-move kernels, the perturb move and the rule redraw at nodes with no grandchildren, are compiled into every build with a default weight of zero. No default fit reaches either, yet about 330 lines plus an enumerator ship, with a permanent documentation and test surface for code no user reaches. The alternative not taken is removing both now that each has cleared its own kill criterion at zero measured benefit for a nonzero default share. The maintainer ruled on 2026-09-14 that both stay in the package at weight zero, unchanged, until the mixing research now scheduled has run; their removal is not an open pre-release question.
Record: the TODO file's tree-mixing entry and docs/design/level-fibre.md. Marked: mine. [dec-A47]

**The R version floor rises to 4.2.0**
The package requires R 4.2.0 or newer, where the released package required 3.1-0. No alternative was weighed, and the number is argued only in the commit message that set it: docs/design/core-generalization.md still says the C++20 toolchains need R 4.3, which nothing in the build enforces. Every user on R 3.1-0 through 4.1.x who worked under 0.9-x can no longer install this version. The maintainer has claimed the decision as their own.
Record: commit 72fc8b3e, the only place the number is argued. Marked: mine. [dec-A48]

**Solaris and big-endian support are dropped**
The package no longer supports Solaris or big-endian hosts, and eight autoconf macro files and the hardware-capability mapfile logic are deleted with them. No alternative was weighed. A user on such a host cannot build the package at all. The maintainer has claimed the decision as their own.
Record: docs/plans/archive/autoconf-dead-code.md, with no maintainer quote. Marked: mine. [dec-A49]

**The engine still calls into R**
Below the bridge the engine is not free of R: the model header calls R's math library for densities, the sampler and chain headers reach R's print and error entry points, two struct invariants exist only because an R error jumps out of the bridge, and the struct describing the sampler's shape carries the sizes of the bridge's output arrays. No alternative was weighed. Nothing is visible to an R user, but a non-R host would have to supply R's math library and its print and error entry points, no shape or source field may ever own storage, and a new family means a new field on a struct every sampler carries. The maintainer ruled a staged route after the release: the bridge gains an internal plain-C specification that the R parse fills and the engine is built from, with family selection a field of it, the engine's five R touchpoints go through hooks the host installs so that draws stay bit-identical, and a plain-spec creation entry and a host-neutral error contract arrive in a 1.x release once the struct settles. See also: [dec-B85].
Record: code only; the TODO file's Python-binding entry asserts the opposite. Marked: blank. [dec-A50]

**Test accessors compile into the shipped engine**
About thirty accessors that exist only for tests are compiled into the shipped engine, two of them virtual on the response base class. No alternative was weighed, and no record discusses it. Nothing is visible to a user; the production vtables and object layout are shaped by the test harness, and a reader cannot tell the model from the scaffolding. Not yet ruled on.
Record: no record in either evidence sweep. Marked: not mine. [dec-A51]

**A second handle layer ships for tests**
A second R-level handle layer, with its own validation, ships in the installed package so that tests can reach entry points the namespace does not export. No alternative was weighed. Every user gets two copies of one validation path that can drift apart, and 62 of about 150 test files are coupled to internals. Not yet ruled on.
Record: docs/plans/review-2026-08-24/consolidated-report.md parks the handle layer's fate as a maintainer item. Marked: not mine. [dec-A52]

**The prior calibration exists twice**
Drawing from the prior predictive re-derives the engine's calibration of the variance prior in R and builds a fresh sampler for each call. No alternative was weighed; asking the engine for the calibration was never considered. Nothing is visible to a user unless the two drift apart, but one calibration now has two implementations in two languages that must agree. Not yet ruled on.
Record: docs/plans/archive/prior-predictive.md spells out the R derivation as the design. Marked: not mine. [dec-A53]

**Three internal helpers become public**
Three helpers written for the work that moved some computation from C++ into R are exported: a validator for model compositions and two primitives that augment data one observation at a time. No alternative was weighed; the maintainer's stamp covers that the work lands, not what it exports. A user gains a general harness for simulation-based calibration and two augmentation primitives as permanent documented surface, the augmentation helper staying scalar because the C entry takes a scalar noise scale. The maintainer evaluated the exports on 2026-09-08 and ruled that they stay, with the manual stating the single-sweep function's contract as the part that will not change.
Record: docs/design/r-c-division.md, whose adoption slate carries a maintainer stamp, with the export shape resolved under a delegated grant. Marked: not mine. [dec-A54]

**Eight methods exist only to refuse**
Eight registered S3 methods have bodies that do nothing but stop with an error, and the generic that dispatches them takes a type argument accepting exactly one token. No alternative was weighed. A user calling one gets a clear refusal; the package carries eight namespace lines and eight bodies to keep aligned with every new fit class, and an abstraction that looks extensible has one instantiation. Not yet ruled on.
Record: docs/plans/review-2026-08-24/consolidated-report.md parks both as maintainer items. Marked: not mine. [dec-A55]

**A second way to constrain interactions**
The package exports blocks() alongside the existing idiom of naming groups in the interactions constructor. No alternative was weighed and the stamp on the row is tacit. A user now chooses between two constraint constructors, each with a manual page and tests: blocks fixes each group's tree capacity, while naming groups in interactions only forbids splits across groups and lets the allocation float, so the two express different priors over allocation rather than one approximating the other. Not yet ruled on.
Record: docs/design/interaction-constraints.md, which says it was built because the maintainer wanted the fixed-capacity guarantee, not because the adaptive path was insufficient. Marked: not mine. [dec-A56]

**Windows on ARM ships unproven**
The package ships NEON support for Windows on ARM64, written before any native probe existed and hedging the architecture string across three spellings. No alternative was weighed. A user on that platform may or may not get a working build: only the Rtools gcc and clang path is proven, a full package check has never run there, and the MSVC intrinsic arm is unexercised. Not yet ruled on.
Record: code only, commits a1cf3c60 and c937f394. Marked: not mine. [dec-A57]

**The documentation site deploys from the branch**
The public documentation site is rebuilt from the development branch on every push that is not a pull request and does not touch the ignored paths, which are the documentation directory, the backlog file and the benchmarks. No alternative was weighed. Anyone reading the site sees unreleased 1.0-0 API descriptions and suggested dependencies before the coordinated release. The maintainer has claimed the decision as their own.
Record: code only, in the site workflow's branch list. Marked: mine. [dec-A58]

**Citation grammar is a hard gate**
The grammar for code citations in the documentation is checked on every push and fails the build unconditionally. No alternative was weighed. Nothing is visible to a user, but every future design or plan edit that adds or changes a citation must conform exactly, and whether to tighten the rule further is open and undecided. The maintainer has claimed the decision as their own.
Record: the TODO file, which records the tightening as a maintainer decision with no plan document until it is ruled. Marked: mine. [dec-A59]

**Seed-locked expected values are a tripwire**
Expected values that depend on the random stream live in exactly four test files, labelled a drift tripwire rather than a correctness test, and are regenerated wholesale by a script. No alternative was weighed. Nothing is visible to a user, but any change that shifts draws costs a regenerate-and-read pass, and a careless regeneration silently blesses a defect. Not yet ruled on.
Record: docs/plans/README.md and docs/plans/archive/snapshot-tests.md, neither naming the maintainer. Marked: not mine. [dec-A60]

**Equivalence is measured against ourselves**
Bitwise equivalence is required only against baselines recorded from the new engine, never against the old one, running on every push on a pinned arm64 reference build with cross-host comparisons and a statistical mode on a schedule. No alternative was weighed. The reason given at the time was that deleting the old engine left nothing wider recordable than the nine-scenario statistical comparison taken at the cutover, whose largest absolute z was 3.83 over 329 summaries and which the manifest marks as evidence and not as a gate. That reason is wrong, and the entry stands superseded: installing dbarts 0.9-34 from CRAN alongside the branch records a comparison as wide as one cares to make it, and a 26-scenario one now exists. What survives is the narrower point, that no gate fires on a shift against the released package, because the wider comparison is run by hand against a hand-installed old library rather than in continuous integration, so a future engine change could separate a scenario again and go unnoticed between runs. Not yet ruled on. See also: [dec-B120].
Record: the baselines manifest header, docs/plans/archive/equivalence-ci.md and docs/architecture.md. Marked: not mine; superseded by dec-B120. [dec-A61]

**The gate policy is prose, not script**
Gates are ordered by a class each change declares for itself, saying whether it moves the random stream, and a re-recorded baseline that moves any draw must name an oracle in its manifest row. No alternative was weighed. Nothing is visible to a user; neither rule is enforced by a script, so a change misclassified as stream-neutral skips the equivalence gate entirely. Not yet ruled on.
Record: docs/plans/README.md and the baselines manifest header. Marked: not mine. [dec-A62]

**The test-count floor sits far below the count**
Continuous integration fails if the test suite runs fewer than 5200 assertions, with six per-file floors underneath. No alternative was weighed. Nothing is visible to a user; the floor sits far enough below what the suite actually runs that a refactor could empty a large part of it and still come up green, which is the opposite of what the floor is for. Not yet ruled on.
Record: code only, in the sanitizer workflow. Marked: not mine. [dec-A63]

**Superseded baselines stay in the tree**
The repository keeps its superseded equivalence baselines, 24 of them, and keeps the near-duplicate test file names left behind by the suite consolidation. No alternative was weighed. Nothing is visible to a user; about 15 MB of superseded baselines ship in the repository, and a contributor faced with a pair of near-identical names cannot tell which file to edit. Not yet ruled on.
Record: the baselines manifest role column and docs/plans/review-2026-08-24/gate-ledger-read.md, which leaves retention open as a maintainer judgement. Marked: not mine. [dec-A64]

**A standing grant settled user-facing forks**
A standing permission to proceed at the implementers' discretion let user-facing forks be settled without being put to the maintainer, and the multinomial family's surface defaults were adopted that way while the maintainer was unavailable. No alternative was weighed. A whole class of user-facing decisions therefore carries no maintainer attribution, the plans recording only that the maintainer may veto. On 2026-09-12 the maintainer ruled that the multinomial defaults stand: a row of counts with zero trials is refused, category names come from the count matrix's column names or else the numbers 1 to K, and a row that must leave the likelihood mid-run does so through the active-row mask.
Record: docs/plans/archive/multiforest-predictor-mutation.md, docs/plans/archive/multiforest-extension-surface.md and docs/plans/archive/multinomial-counts.md, whose defaults were adopted as working defaults by the implementers while the maintainer was away. Marked: blank. [dec-A65]

**One spelling for the noise estimate**
The estimate of the residual standard deviation supplied when a sampler is created is spelled sigest on every entry point, dbarts included. The alternative, which shipped first, kept two spellings for one concept, sigma on dbarts and sigest on the other three entry points. A user passing the old name to dbarts gets a warning for one release; the sampler's setter for the parameter itself is untouched, which is the distinction the maintainer asked be preserved. See also: [dec-B80].
Record: docs/plans/archive/bart2-argument-consolidation.md, under its own decided heading with no maintainer quote. Marked: blank. [dec-A66]

**Four common nouns enter the search path**
The package exports interactions, blocks, forest and varianceForest as bare top-level names, while the prior constructors are bundled into one list to keep generic names out of a user's search path. No alternative was weighed, and only the blocks half carries even a tacit stamp. A user attaching the package takes four ordinary English nouns into the search path, under a rule the same design document states in the opposite direction about a hundred lines away. Not yet ruled on.
Record: docs/design/interaction-constraints.md for blocks, and docs/design/public-surface.md for the priors, which decides to evolve in place and keep the no-pollution property. Marked: not mine. [dec-A67]

**Some documented arguments do nothing**
Three documented arguments are inert rather than refused: the run method's thread count, the multinomial draws method's variable selector, and the entry for tau in the default variable list of the summary method and the four draws pairs. No alternative was weighed beyond refusing them. A user reading the manual finds arguments documented as doing nothing, the last of them in eight public signatures; those draws methods were named as_draws until the ruling that removed the posterior dependency replaced them with a plain draws extractor. Not yet ruled on. See also: [dec-B99].
Record: docs/plans/surface-refusals.md, which records that they stay inert and documented rather than refused, attributed to the maintainer; the tau default has no record and the family that carried it has been deleted. Marked: not mine. [dec-A68]

**Forest indices start at one in R**
The forest index on the reference-class sampler is 1-based while the C interface is 0-based, and the bridge converts between them. No alternative was weighed. A user indexes forests from 1 as they would anywhere else in R; the package carries two index origins and every bridge call converts. Not yet ruled on.
Record: docs/plans/archive/multiforest-extension-surface.md, which settles 1-based indexing converted at the boundary, attributed to the maintainer with no fork recorded. Marked: not mine. [dec-A69]

**A missing response no longer stops the fit**
The modern entry points take a standard na.action argument whose default drops rows with a missing response, keeps missing predictors for the trees, and records the dropped rows so that training-set fitted values pad back to the length of the data as na.exclude does; the base functions keep their usual meanings, so na.omit drops any row with any missing value, na.fail errors, and na.pass keeps everything and then the response check errors. The alternative, which shipped first, was to error on any missing response where 0.9-x dropped those rows, the archived sweep leaving the choice between letting it break and softening it to the maintainer with no ruling recorded. A 0.9-x call on data with missing responses fits the complete cases again instead of erroring; the reverse-dependency sweep found the package insight breaking on exactly this, and that break is what the default addresses. The maintainer ruled the na.action argument in, to land in the family-objects consolidation pass. See also: [dec-B108].
Record: docs/plans/archive/cran-readiness.md, the reverse-dependency sweep run with suggested packages. Marked: not mine. [dec-A70]

**Gaussian-process leaves get no general move**
The move strategy for non-conjugate Gaussian-process leaves, the sixth phase of the founding engine design, stays designed but not built and not scheduled. No alternative was weighed. A user fitting Gaussian-process leaves under a general likelihood has no move strategy for them until a consumer or the maintainer asks for one. Not yet ruled on; it was never put to the maintainer and is recorded as not scheduled.
Record: docs/design/gp-leaves.md. Marked: blank. [dec-A71]

**The treatment-scale rescale move is unbuilt**
The joint rescale move for the treatment forest's scale in the causal-forest family stays derived and checked in a prototype but not implemented and not scheduled. No alternative was weighed. A user fitting causal forests does not get for the treatment forest the mixing improvement that the equivalent landed move gives the prognostic forest. Not yet ruled on; it was never put to the maintainer and is recorded as not scheduled.
Record: docs/plans/archive/bcf-b-ridge.md. Marked: blank. [dec-A72]

**Three more C entry points stay unbuilt**
The remaining additive work on the C interface, a plain-C getter for the noise scale, saving and restoring a data handle, and variable-selection inference, stays recorded but not built and not scheduled. No alternative was weighed. A compiled consumer wanting any of the three gets no C entry point until the maintainer asks for one. Not yet ruled on; it was never put to the maintainer and is recorded as not scheduled.
Record: docs/plans/prerc-surface-freeze.md. Marked: blank. [dec-A73]

## B. Decisions with maintainer evidence

**Missing predictors are modelled, not refused**
A missing predictor value no longer stops a fit: the sampler sends the missing cases down one side of a split, chosen as part of the rule, so the missingness itself can carry signal. The alternative, refusing any design with missing values, survives as a strictness setting named "error". A user fitting a design with missing predictors gets a fit rather than an error, and both the number of rows used and the fitted values differ from 0.9-x. The maintainer set the name and the default to "incorporate" and kept "error" as the strictness escape. At the tip that escape is reachable only by writing the slot on the model object, since no front-door argument spells it; the modern front door instead takes an na.action argument that drops rows with a missing response, and the BayesTree-style door drops any row with a missing value. See also: [dec-B108], [dec-B75].
Record: docs/design/mia-missingness.md. Marked: mine; superseded in part by dec-B108: na.action replaced the front-door spelling, so the escape is slot-only and bart no longer rejects NAs. [dec-B01]

**The swap tree move ships switched off**
The tree kernel keeps the swap proposal, which exchanges split rules between a node and its parent, but its default probability is zero and its share of the proposal mass goes to the birth and death moves. The alternative, and the first ruling, was deleting the move from the kernel outright; that removal was then partly reversed. Every fit's draws differ from 0.9-x because the mixture of proposals changed, and a user fitting a single tree who wants swap back has to set it. The maintainer first ruled that "The swap move is removed from the MCMC kernel before 1.0", then that "The removal is PARTIALLY reversed", which leaves the move present at zero weight.
Record: docs/design/swap-removal.md. Marked: mine. [dec-B02]

**The change move's acceptance ratio is repaired**
The change move, which redraws the split rule at an internal node, was missing the proposal-density term in its acceptance ratio; the repair uses restricted proposals with explicitly counted ratios for ordinal variables and draws from the prior for categorical ones. The alternatives were other variants measured against each other, and the maintainer picked this hybrid. Every 1.0-0 fit therefore differs from every dbarts and BayesTree fit since the package began, with no switch back to the old behaviour. The hybrid "is now THE change move".
Record: docs/design/change-move-balance.md; docs/plans/archive/change-move-fix.md. Marked: mine. [dec-B03]

**chi() means the degrees of freedom it names**
The chi hyperprior on the leaf-scale parameter k drew from a distribution with twice the shape its argument named; the constructor now samples the stated degrees of freedom, and the binary default is relabelled so the prior it names is the prior it draws. The alternative was relabelling the argument to match the old behaviour, a documentation fix rather than a code fix; the maintainer picked the code fix. A user who wrote an explicit chi with 1.25 degrees of freedom in 0.9-x gets a different prior now, while the relabelled default draws the same values it always did.
Record: docs/plans/correctness-audit.md. Marked: mine. [dec-B04]

**bart and bart2 were to differ by intent**
The package was to keep two entry points that differ permanently: bart reproducing BayesTree's defaults and bart2 carrying the current recommendations, with the same functionality behind both. No alternative was weighed. The same data would give different answers depending on which name a user called. That is no longer the shape: a later ruling makes bart itself the modern front door and moves the BayesTree-style names and defaults to a new function. See also: [dec-B75].
Record: docs/plans/archive/bart2-argument-consolidation.md. Marked: superseded by dec-B75. [dec-B05]

**setResponse keeps its argument order**
The sampler's setResponse method keeps the argument order it shipped with in 0.9-34, and warns once per session when the second argument is supplied positionally. That order matches the offset setter's, which is the reason it was kept; reordering it was the alternative. An embedding loop written against 0.9-x that passes the second argument positionally silently changes meaning after that single warning. The maintainer ruled to "Keep setResponse's argument order ... and warn once per session when the second argument is supplied positionally."
Record: docs/plans/release-candidate-review.md. Marked: mine. [dec-B06]

**One argument order for every predict method**
Every predict method in the package takes its arguments in the same order: the fitted object, the new data, the kind of prediction, then anything further. The alternative, rejected, was leaving the orders as they were and only refusing unrecognized argument names. A 0.9-x positional call silently changes meaning wherever the order moved, so a user who relied on position has to check. The maintainer chose the single order.
Record: docs/plans/prerc-surface-freeze.md. Marked: mine. [dec-B07]

**Fractional counts are refused, not truncated**
Wherever a caller can write a count, a column index or a group selector by hand, a fractional value is refused by name rather than silently truncated. No alternative was weighed; the ruling was a directive to finish the rule everywhere. A 0.9-x script that computes a count and passes it without rounding now errors. A later ruling narrows the rule to counts and selectors only, never to weights. See also: [dec-B103].
Record: docs/plans/release-candidate-review.md. Marked: superseded by dec-B103. [dec-B08]

**Unknown arguments warn; one is refused**
Each predict method refuses an argument named value by name, and any other unrecognized name in the dots produces a warning rather than an error, through a check the package implements itself. The alternatives were erroring on any unknown name or ignoring them silently; the maintainer chose warning, "grounded in R standards", including when the name arrives through a subclass's own method dispatch. A user whose wrapper forwards its own dots into a dbarts call sees a warning, and code that relied on the 0.9-x deprecation warning for value now gets an error. The ruling named base R's chkDots for the job; the package's own utility source states why it does not use it.
Record: docs/plans/release-candidate-review.md; docs/plans/archive/bart2-argument-consolidation.md. Marked: mine. [dec-B09]

**The modern front door stops accepting dots**
The modern front door takes no dots channel: an argument it does not name is refused, so every control knob has to be a named argument of its own. The alternatives ranged from keeping dots as a pass-through to deleting them outright, and the maintainer picked refusal with a stated sunset toward deletion. Any control field not promoted to a named argument was therefore unreachable from that entry point. That rule no longer holds: the front door and the cross-validation entry both take a control object again, so control settings are reachable without one argument apiece. See also: [dec-B116].
Record: docs/plans/archive/bart2-argument-consolidation.md. Marked: mine; dec-B116 reverses the formals-only rule, giving bart and xbart control = dbartsControl(). [dec-B10]

**xbart refuses a three-element n.burn**
The cross-validation entry point refuses a burn-in vector longer than two elements by name, where 0.9-34 accepted three and ignored the third. No alternative was weighed. A user who copies the documented 0.9-x default gets an error.
Record: docs/plans/release-candidate-review.md. Marked: not mine. [dec-B11]

**The cross-validated k grid for binary fits**
Cross-validation over the leaf-scale parameter k for a binary outcome was to use a grid of fixed values, matching the gaussian arm, rather than the hyperprior the modern default fits. The alternative was leaving the binary axis on its own semantics. The cross-validated k then no longer matches what the front door fits by default, and 0.9-x results on that axis are not comparable. A later ruling reopens it: the grid now accepts fixed values and an entry that leaves k modelled. See also: [dec-B102].
Record: docs/plans/release-candidate-review.md. Marked: superseded by dec-B102. [dec-B12]

**Weighted binary fits mean counts, not scaled latents**
The weighted probit path is gone rather than carried into the new engine, because scaling the latent draws by one over the square root of the weight is not a coherent weighted likelihood. Weights on a logistic binary fit are supported under the common reading of a weight as a number of identical observations, restricted to positive integers. A 0.9-x weighted probit call now errors, and the only migration is an integer-count logistic fit. The maintainer's words on the removal: it "is not a coherent weighted-likelihood model".
Record: docs/design/core-generalization.md; docs/design/weighted-logistic.md. Marked: mine. [dec-B13]

**No approximate sampling in this release**
Every draw in 1.0 stays exact: an approximate Polya-Gamma draw, the data augmentation that would make a logistic likelihood conjugate under arbitrary real weights, is declined. The alternative was accepting that approximation to unlock the models behind it. A real-valued negative-binomial dispersion and arbitrary real weights on a binary fit therefore cannot be fit in 1.0-0, and both wait on that one unlock. The ruling also keeps the hurdle model's sampler-only mode refused, and struck a proposal to delete the twin-create path as relitigation of a settled point.
Record: docs/plans/release-candidate-review.md. Marked: mine. [dec-B14]

**Heavy tails and overdispersion on capped grids**
A Student-t residual law estimates its degrees of freedom on a capped grid with a floor of 3, and negative-binomial dispersion is drawn on a capped grid of integers. The alternative for dispersion was a continuous parameter, which needs an approximate draw; the resolution on record is "integer dispersion, fully exact". A user cannot have the degrees of freedom estimated below 3, which is where the residual variance stops existing, though a value below 3 supplied outright is accepted and fit; a non-integer or very large dispersion has no such escape.
Record: docs/design/robust-errors.md; docs/design/negative-binomial.md. Marked: not mine. [dec-B15]

**Monotone leaves target the exact posterior**
Leaves constrained to be monotone in a predictor sample the exact target, with a constrained joint marginal at the seam where a tree move changes which constraints apply. The alternatives traded that exactness for a simpler move at the seam. The cost is 600 to 1000 lines of engine code and an adaptive quadrature that dominates the run time of a monotone fit. The record has the maintainer resolving for the exact shape and for proceeding with the monotone work at that point.
Record: docs/design/monotone.md. Marked: not mine. [dec-B16]

**The variance forest rides the weight channel**
A second forest modelling the residual variance enters the sampler through the observation weight channel, as a separately typed forest that may be absent. The alternative was giving it a channel of its own. Because the weight channel is then taken, a variance forest is refused for every family and model composition that already owns weights. The record has the maintainer resolving for this route and for implementing it at that point.
Record: docs/design/heteroscedastic.md. Marked: not mine. [dec-B17]

**The causal model's prognostic scale is half-Cauchy**
In the two-forest causal parameterization, which fits a prognostic forest and a treatment-effect forest, the prognostic forest's amplitude gets a half-Cauchy prior. The alternative was a different scale prior; parity with the bcf package was the stated goal. The R-side default is family-aware while the engine's own fixture default is not, so the two disagree by design and a user reaching the engine directly sees the plainer one.
Record: docs/design/bcf.md. Marked: not mine. [dec-B18]

**Causal forests live in bartCause, not dbarts**
dbarts exports no causal-forest entry point: the causal fit class and the causal argument vocabulary move to bartCause, and what stays in dbarts is the engine vocabulary for models with more than one forest. The alternative, a public creation route inside dbarts, existed during the pre-release and was withdrawn before it shipped. A dbarts user cannot fit a causal forest by name. The maintainer: "I don't think bcf belongs in the dbarts function."
Record: docs/design/bcf.md; docs/plans/archive/bcf-public-surface.md. Marked: mine. [dec-B19]

**An ordered factor now fits an ordinal model**
An ordered-factor response is detected by its class and fit as an ordinal model, identified by fixing the scale, with the cutpoints updated one at a time by a Cowles-style Metropolis step. The alternatives were requiring the user to name the ordinal family and inferring it from the level names; the condition on record is detection by a concrete class, never by level names, with the fit announcing what it has done. A user whose 0.9-x code passed an ordered factor got a continuous fit on the integer level codes and now gets a different model class, announced at the call. The BayesTree-style door refuses such a response by name rather than changing the model under the user.
Record: docs/plans/archive/ordinal-outcomes.md; docs/design/ordinal.md. Marked: not mine. [dec-B20]

**Latent draws are not readable for multinomial fits**
The sampler's latent accessor declines in writing when the fit is multinomial, instead of returning the augmentation draws. The alternative was returning them. A host embedding a multinomial sampler inside a larger model cannot read those latents. The maintainer's criterion was usefulness: "latents are not persisted absent a compelling use case."
Record: docs/design/multinomial-mutation-arc.md. Marked: mine. [dec-B21]

**Gaussian-process leaves ship with no users**
The Gaussian-process leaf model, in which a leaf holds a smooth function rather than a constant, stays in 1.0-0 although nothing outside the package uses it. The alternative was cutting it before the release. Keeping it freezes 773 lines of engine code, a factorization whose cost is cubic in leaf size, and a manual page into the release. The maintainer ruled to keep it while flagging that whether it belongs in 1.0 is a conversation still owed.
Record: docs/plans/pre-review-cleanup.md. Marked: mine. [dec-B22]

**Single-precision residual storage is optional**
The sampler can hold residuals in single precision, asked for by a storage argument set to "single", with double precision still the default. The alternative was not building it. A user who asks for it gets different numbers back, which the manual page has to say, and the knob is valid for one combination of leaf model and family and refused for the rest. The maintainer directed: "Build it, keep it optional, think about other ways to optionally decrease storage sizes at the same time." The cost is one more full instantiation of the engine stack.
Record: docs/design/reduced-precision-storage.md; docs/plans/archive/bart2-argument-consolidation.md. Marked: mine. [dec-B23]

**Observation indices narrow to 32 bits**
The engine stores the indices it gathers observations through as 32-bit integers by default, and the narrowing preserves every draw bit for bit. No alternative was weighed; the maintainer directed it. The number of observations is capped just under four billion, and every support-library signature that takes an index changed with it.
Record: docs/design/reduced-precision-storage.md. Marked: mine. [dec-B24]

**No eight-bit predictor codes**
Quantized predictor values stay in 16-bit codes: an eight-bit layer for the hot data and per-column code widths are not pursued. The alternative was a standalone phase of work to build them. Nothing a user sees changes, and the memory-bound regime at large n keeps the wider codes. The go or no-go on record reads "NO standalone phase 2."
Record: docs/plans/archive/hot-layer-u8.md. Marked: not mine. [dec-B25]

**The support library's thread managers are cut**
The two thread managers in the support library and the threaded wrappers around the moment calculations are removed from the tree, after being archived on a branch. The alternative was keeping them against a future scheme that threads within a single chain. About 2470 lines of within-chain reduction machinery leave, so such a scheme would start over. Nothing a user sees changes.
Record: docs/plans/pre-review-cleanup.md. Marked: not mine. [dec-B26]

**The response is still scaled by its range**
A continuous response is still shifted and scaled to the interval from minus one half to one half before fitting, rather than standardized by its standard deviation. The alternative was standardization, which is less sensitive to a single extreme response value. Range scaling keeps that sensitivity and the internal-scale bookkeeping it forces, and was kept for compatibility with the package's lineage. The maintainer signed off to "keep, document", with a setter added on the sampler for the response.
Record: docs/plans/archive/range-scaling.md. Marked: mine. [dec-B27]

**The residual sum of squares is rescaled**
The sum of squared residuals the sampler reports is de-scaled by the square of the response range, which is the correct conversion back to the data's own units. No alternative was weighed; it was a units slip with no consumers. The value returned differs from every released dbarts. It was fixed in both engines so that the two agreed.
Record: docs/design/core-generalization.md. Marked: not mine. [dec-B28]

**Categorical rules report directions, not a mask**
The tree reader returns a missing value in the split-value column for every categorical rule and gives the per-level directions instead, dropping the raw bit-mask value that narrow categorical rules used to carry. The alternative was keeping two vocabularies, one for narrow masks and one for wide ones. A 0.9-x script that read the value column for a narrow categorical rule now gets NA there. The maintainer signed off to "unify".
Record: docs/plans/archive/flat-format-v2.md. Marked: not mine. [dec-B29]

**family moves from control to model**
The response family is a slot on the model object rather than on the control object. No alternative was weighed. Any consumer that reads the family off the control object breaks, and the layout freezes at 1.0-0. The maintainer signed off to "move".
Record: docs/plans/archive/family-on-model.md. Marked: mine. [dec-B30]

**Saved sampler state carries a format version**
A saved sampler state records the format it was written in and is refused when that does not match the running package, with no shim to migrate an old one. The alternative was such a shim. A state written by any pre-release build cannot be loaded. The maintainer scoped the policy down on the ground that nothing from this line of work has been released, so there is no cross-version guarantee apparatus and nothing for binary formats. A later ruling adds the one recognition that matters, a 0.9-x fit refused by name with a message. See also: [dec-B109].
Record: docs/plans/archive/state-format-policy.md. Marked: not mine; superseded by dec-B109. [dec-B31]

**Restoring a state is semantic, not bitwise**
Loading a saved sampler state restores the model's meaning rather than its bits: the residual standard deviation rides the original response scale and fitted values are rebuilt by resumming the trees. The alternative was a bitwise restore. A consumer cannot check reproducibility by requiring bitwise equality across a save and a load. The maintainer signed off to "drop to semantic restore".
Record: docs/plans/archive/state-continuation.md. Marked: not mine. [dec-B32]

**A fit does not store its sampler state**
Trees are not kept by default, a fit does not capture sampler state on its own, and predicting from a reloaded fit errors with a message naming what to call first. The alternatives were flipping the default to keep trees, rejected, and the eager state capture that was built and then reverted. A user who saves a fit and reloads it hits that error unless the state was stored before saving. The maintainer overruled the eager capture because the resulting bloat, roughly 2.8 times the size of the fitted values, is not acceptable as a default, and accepted the requirement to touch the state as a known cost. A later ruling restates this with the manual wording and a fallback. See also: [dec-B104].
Record: docs/plans/archive/package-review-remediation.md; docs/plans/prerc-surface-freeze.md. Marked: superseded by dec-B104. [dec-B33]

**A new missing value at predict is refused**
Predicting on data whose column carries a missing value where the training column had none is refused by name. The alternative, documenting that such a row silently goes down the left branch, was rejected. A user's test frame with a new missing value errors where 0.9-x dropped the row.
Record: docs/plans/prerc-surface-freeze.md. Marked: not mine. [dec-B34]

**The data object's x slot is the source**
The x slot on the data object accepts any type and means the predictors as they were handed in, not a maintained view of the engine's quantized state. The alternative was keeping it a matrix that mirrors what the engine holds. The S4 slot no longer states the contract, so a consumer that assumes a matrix is unguarded. The maintainer approved the reconciliation conditional on the model being a collection of mutations and explicitly "NOT a maintained public view of the engine's quantized state".
Record: docs/plans/archive/data-ownership-3-mutation.md. Marked: mine. [dec-B35]

**Predictor matrices are borrowed, never copied**
The engine borrows a predictor matrix read-only while the sampler is being constructed, and never afterwards; every later path that needs the raw values gets them handed back at call time. The alternative, copying the raw values into the sampler, was rejected. Mutation, cut-point selection and tree replay all then require the R layer to supply the values again, which is work the R layer has to do. The maintainer rejected the copy and set both the construction-only borrow and the call-time supply.
Record: docs/design/data-ownership.md. Marked: mine. [dec-B36]

**sparseFactor ships under a Matrix-style name**
The package exports a constructor for a sparse factor predictor, named to follow the Matrix package's own naming convention so that it reads familiarly to users of sparse data. No alternative spelling was weighed once the convention was chosen. The cost is one more generic-sounding name at the top level of the package.
Record: docs/design/data-ownership.md. Marked: mine. [dec-B37]

**Probability vectors snap to sum to one**
Split and proposal probability vectors that sum to one within the square root of machine epsilon are snapped to exactly one. The alternative was the tighter fixed tolerance the released package uses. This tolerance is looser, so nothing that validated under 0.9-x is refused now, and some input that would have been refused is accepted. The maintainer asked for "the default tolerances from almost-equal type functions" and accepted that the change moves which user input is admitted. A separate site snaps near-zero multipliers to exact zero and shares only the constant with this one.
Record: docs/plans/archive/zero-weight-exactness.md. Marked: not mine. [dec-B38]

**The R mutation path pays its cost**
Changing predictors, a response, an offset or weights between draws through the R layer collects the data on every update, and that cost is documented rather than removed. The alternative, a flag that opts out of the collection, was considered and not built. An embedding loop pays the collection cost per mutation, and the documented escape is to drive the sampler through the C interface instead. The decision on record is to accept and document.
Record: docs/plans/archive/data-ownership-3-mutation.md. Marked: not mine. [dec-B39]

**A heteroscedastic fit refuses a response rescale**
A sampler whose leaves carry a variance scale refuses an update of the response scale rather than rescaling itself. The alternative was performing the rescale, which is left to an additive entry point after the release. A host embedding a heteroscedastic fit cannot re-anchor its scale in the middle of a chain. The maintainer ruled for the refusal after the alternatives, the tradeoffs and a recommendation were put.
Record: docs/plans/release-candidate-review.md. Marked: mine; superseded by dec-B122. [dec-B40]

**Gaussian-process leaves want few trees**
The manual recommends 10 to 25 trees for a fit with Gaussian-process leaves, rather than the ensemble-scale tree count the package defaults to. The alternative was letting the general default stand for this leaf model too. The package's own default is then wrong for one shipped leaf model, and a user has to override it. A later ruling adds a self-report as well: the fit counts leaf evaluations that fell back to constant leaves and warns when that share is high. See also: [dec-B110].
Record: docs/plans/release-candidate-review.md. Marked: not mine; superseded by dec-B110. [dec-B41]

**Student-t log-likelihood is the t marginal**
Under Student-t residuals the reported log-likelihood is the observation-level t density, marginal over the augmentation variables, not the gaussian density conditional on them. The alternative was that conditional value. The marginal is what keeps the channel comparable across residual laws, being the observation-level density that the widely applicable information criterion and importance-sampling leave-one-out are defined on; the conditional would not be comparable with the gaussian channel. A user computing either criterion from the channel gets a number they can compare across families.
Record: docs/plans/release-candidate-review.md. Marked: not mine. [dec-B42]

**Two sampler accessors refuse a result argument**
The sampler's accessor for the residual standard deviations and its accessor for the sums of squared residuals refuse an argument named result by name. The alternative was continuing to ignore it. Code from 0.9-x that passes that argument now errors instead of having it silently dropped.
Record: docs/plans/surface-refusals.md. Marked: not mine. [dec-B43]

**makeind keeps a formal that does nothing**
The design-matrix helper keeps the signature BayesTree gave it, including an argument named all that has no effect. The alternatives were implementing the argument or dropping it, and neither was done this round. A documented no-op argument ships, so a user who sets it sees no change in the result.
Record: docs/plans/pre-review-cleanup.md; docs/plans/surface-refusals.md. Marked: mine. [dec-B44]

**The C header renames ordinal thresholds**
In the shipped C header the ordinal model's latent thresholds are named as thresholds, leaving the term cut points for the split grid alone. The alternative was leaving the older name in place. A consumer compiled against the old symbol breaks, and the header's identity hash has to be recomputed. The maintainer ruled yes, rename.
Record: docs/plans/pre-review-cleanup.md. Marked: mine. [dec-B45]

**DESCRIPTION drops credits for deleted code**
The package description no longer credits authors whose code is not in the tree. The alternative was keeping the credits as a courtesy. Four contributors of autoconf macros and the author of a radix-tree implementation lose their credit line. The maintainer ruled to "REMOVE the credits for code that is gone".
Record: docs/plans/pre-review-cleanup.md. Marked: mine. [dec-B46]

**The forest block was to ship whole**
The forest and amplitude block of the flat C interface was to ship complete, including two entry points no known consumer calls, on the reasoning that consumer absence is not the gating fact. The alternative was shipping only what a consumer calls. Twenty-three entry points that nothing known calls would freeze for the life of the interface. That is no longer the shape: a later ruling ships the entries a consumer calls plus the small queries a host needs, and holds the multi-forest block back until a consumer or a named host design asks for it. See also: [dec-B86].
Record: docs/plans/pre-review-cleanup.md; docs/plans/archive/cheap-uniformity.md. Marked: superseded by dec-B86. [dec-B47]

**An item ships unless nothing valuable follows**
The standing gate on what to build is whether anything valuable could follow from it: an existing consumer is sufficient evidence of value and never necessary. The alternative is requiring a named consumer before building. The rule licensed the uncalled C entry points and the broad R surface, all of which freeze at release and have to be supported afterwards. The maintainer's standing rule: "gate an item IF we cannot think of anything valuable it might enable". For the C header alone it narrows later, because an entry put there cannot be taken back. See also: [dec-B86].
Record: the root TODO file; docs/plans/archive/cheap-uniformity.md. Marked: superseded by dec-B86. [dec-B48]

**Sister-package behaviour does not argue a design**
A design fork is not argued from what a consuming package happens to do; breakage of the sister packages appears only in migration maps. No alternative was weighed. The maintainer owns bartCause and rejects consumer behaviour as a design input or precedent, so a fork's justification may not cite it.
Record: docs/plans/archive/bart2-argument-consolidation.md. Marked: mine. [dec-B49]

**Large datasets are common; chains stay the default**
The register recorded agents' claim that single-chain workloads at a hundred thousand observations or more are common; the maintainer says that was not said. The standing fact is that large datasets are common and that multiple chains remain the default. No shipped path serves a single-chain speed-up, and the vectorized kernels that were the live remnant of the argument were measured and declined, so nothing is held pending on it. The apparent conflict with the entry that gives a single chain no sampling parallelism at any thread count is a difference of priority against measurement, not two contradictory decisions. See also: [dec-B113], [dec-A37].
Record: the root TODO file. Marked: not mine. [dec-B50]

**New surface only where function is missing**
The package adds an accessor when a capability is otherwise unreachable, not when it would merely save the user a step. No alternative was weighed. The accessor that reports fits without the offset ships with no with-offset twin: the identity that recovers the other is documented rather than built, so a user adds the offset back themselves. The maintainer's principle as stated: "add surface where FUNCTIONALITY is missing, not where recovery is convenient".
Record: docs/plans/adoption-slate.md. Marked: mine. [dec-B51]

**R does the conditionals, C++ does the integrand**
As a standing aspiration, the R layer handles the conditioning and setup while the C++ engine handles the integration, so every engine capability is meant to have an R route to it. No alternative was weighed. The rule binds every later design and obliges that R route for each engine capability, which is work each time. The maintainer asked that the principle "cover what we aspire to and not what we have currently implemented".
Record: docs/design/r-c-division.md. Marked: mine. [dec-B52]

**The whole adoption slate lands before 1.0**
Five surface-adding arcs all land before the release, chosen on utility, with cost used to size a budget and never to rank or gate them. The alternative was deferring some of them past the release on cost. They become permanent 1.0 surface that a user can call and the package then has to keep. The maintainer's mandate put all items pre-release, framed by utility: "price sizes a budget, it never ranks or gates".
Record: docs/design/r-c-division.md. Marked: mine. [dec-B53]

**Multi-forest models are formula terms**
A model with more than one forest is written as terms in the formula rather than through an argument listing forests, and the front door's argument names were to lock at 1.0-0. The alternative was such an argument. Two per-forest settings have no term spelling, so a user cannot reach them that way. The lock has since moved: the names lock when 1.0-0 is pushed to CRAN, not before. See also: [dec-B79].
Record: docs/plans/archive/bart2-argument-consolidation.md. Marked: superseded by dec-B79. [dec-B54]

**forest() means two things by context**
Inside a formula the token forest() marks a forest, with a colon spelling as its canonical short form, and its unnamed slot there is a symbolic set of predictors rather than the data that the constructor of the same name takes. The alternatives were other head tokens, surveyed against R conventions and tried in two rounds of probes. One exported name therefore means two different things depending on whether it appears inside a formula, which a user has to learn. The maintainer chose the token and the colon spelling, and the divergence from the constructor is deliberate.
Record: docs/plans/archive/bart2-argument-consolidation.md. Marked: mine. [dec-B55]

**The C interface detects its own drift**
The shipped header carries a major and a minor version number, generates its list of entry points from one macro-driven source, and computes a hash of that list at compile time which it compares against a value baked into the header. The alternative shapes were a runtime table and a hand-maintained list. The header is dense with macros that every consumer's compiler expands, the hash moves even on a purely additive append so it can only gate consumers that ship in lockstep, and every edit to the interface needs the baked value recomputed by hand. All three parts were decided together.
Record: docs/plans/archive/capi-dispatch-table.md. Marked: not mine. [dec-B56]

**Consumers still look up each symbol**
A compiled consumer resolves each entry point by symbol, as before, rather than fetching one table of function pointers from a single query. The alternative was that table, and it was rejected. Nothing that a user or a consumer loses was identified.
Record: docs/plans/archive/capi-dispatch-table.md. Marked: not mine. [dec-B57]

**The strict header check is opt-in**
The exact hash check on the C interface is a documented opt-in for consumers that ship in lockstep with dbarts, not the default. The alternative, leaving lockstep as the default, was rejected. With the check off, a consumer built against any 1.x header sharing the major and minor number is admitted, and since those constants have never moved that window covers the whole pre-release history. A later ruling drops the flag from the sister packages and leaves it off, with the version pair as the guard. See also: [dec-B111].
Record: docs/plans/prerc-surface-freeze.md; docs/plans/dbarts-h-freeze.md. Marked: not mine; superseded by dec-B111. [dec-B58]

**Structs crossing the interface carry their size**
Every struct that crosses the C interface begins with a size field and, before 1.0-0, grows by appending fields at the bottom. The alternative was freezing the layouts at once and adding a new entry point for anything further. Removing a field is unprotected and possible only before 1.0-0, and every consumer must either set the size field or use the initializer macro the header provides. The resolution on record was to extend the struct then, with a freeze at release time, and later to fold the layouts together. The two halves of that freeze cover different kinds of addition, which the header states: a new function after 1.0-0 arrives under a new name and a minor bump, while a struct still grows by appending a field, and such an append bumps the minor version and re-bakes the header's hash together.
Record: docs/plans/archive/capi-callbacks.md; docs/plans/release-candidate-review.md; inst/include/dbarts/dbarts.h. Marked: not mine. [dec-B59]

**The interface version stays put before release**
The major and minor version constants on the C interface do not move during the pre-release. No alternative was weighed, since no version of this work has been released. The constants therefore carry no information about the pre-release history, and the baked hash is the only detector of drift. The maintainer: "No need to increment versions".
Record: docs/plans/archive/dbarts-h-reshape.md; docs/plans/archive/bcf-public-surface.md. Marked: mine. [dec-B60]

**Two header calls: scope and a parameter name**
Two open questions about the header's shape were settled: the first cleanup item covers seven entries rather than the five originally named, and the basis setter's parameter is renamed to say that the data is row-major instead of transposing the contract to match the old name. The alternatives were the narrower five-entry scope and the transposition. A caller who lays the data out the other way still gets no error, only a parameter name that says which way round it goes. A later ruling took both entries this shaped out of the header. See also: [dec-B86].
Record: docs/plans/capi-shape.md. Marked: not mine; superseded by dec-B86, which took both entries it shaped out of the header. [dec-B61]

**A callback can run between sweeps**
A host can register a callback that fires once per sweep and receives the chain index, with the chains run inline and in a documented sequential order. The alternative was a callback that also works when chains run on worker threads. A host that wants per-sweep conditioning gives up multi-threaded multi-chain runs and cannot see interleaved progress, and on the worker path a registered callback is dropped silently rather than refused.
Record: docs/plans/archive/capi-callbacks.md. Marked: not mine. [dec-B62]

**Error messages follow base-style R practice**
New error and warning messages follow what highly regarded base-style R packages measurably do, and the tidyverse style guides carry no authority. The alternative was adopting those guides. The cost is a repository-wide sweep of the existing messages and a permanent style rule that reviewers hold to. The maintainer resolved for "best practices from highly regarded R packages" and refined it the same day to say that tidyverse practices carry no authority weight.
Record: docs/design/error-style.md. Marked: mine. [dec-B63]

**Protection-balance findings are preempted**
The repository runs the static checker CRAN uses to find R protection-balance errors and fixes what it reports, rather than waiting for CRAN's incoming inspection to raise them; the gate carries no list of suppressed findings. The alternative was risking them at submission. Any protect-balance tag fails the job, so a false positive costs code churn to silence. The maintainer redirected mid-plan: CRAN's additional-issues report gets preempted, not risked.
Record: docs/plans/release-candidate-review.md. Marked: mine. [dec-B64]

**The slow checks start at the merge**
The slow and statistical continuous-integration checks run only on a schedule or on manual dispatch, apart from a push trigger on each check's own workflow file, and they are registered on the default branch at the coordinated merge rather than before it. The alternative was registering them earlier. Five gates cannot fire on their schedules until then, so the merge is their first scheduled run. The maintainer held the default branch as it is until reviewing the branch, and deferred registration to the merge.
Record: the root TODO file; docs/plans/release-candidate-review.md. Marked: mine. [dec-B65]

**A gate's thresholds are judgement, not measurement**
The acceptance thresholds for the multi-forest veto rate, the share of proposals such a model rejects outright, are ratified as judged rather than measured. No alternative was weighed. The bands gated one measurement before an arc and nothing standing: no benchmark script, test file or workflow job re-runs them, so no build can fail on them. The maintainer ratified them on the recommendation, with an explicit caveat about where they came from.
Record: docs/plans/multiforest-veto-rate-falsifier.md. Marked: not mine. [dec-B66]

**The package publishes no coverage badge**
Test coverage is measured locally and on demand, with no badge published. The alternative was a coverage service and its badge. There is no public coverage signal for the package. The maintainer resolved: no badge, coverage stays local and on demand.
Record: docs/plans/repo-modernization.md. Marked: mine. [dec-B67]

**The release procedure has no submission comments**
The release procedure carries no step that writes comments to accompany a CRAN submission. No alternative was weighed. Nothing conveys submission context to CRAN. The maintainer's reason: CRAN does not read it.
Record: commit 6a97236d. Marked: mine. [dec-B68]

**The BayesTree-matching build flag is gone**
The build flags that made the old engine reproduce BayesTree's random number stream are not carried over; they stay with the old engine and die when it is deleted. No alternative was weighed. A user who built the package in the documented BayesTree-matching mode has no equivalent in 1.0-0.
Record: docs/design/core-generalization.md. Marked: mine. [dec-B69]

**The package requires C++20**
Building dbarts requires a C++20 toolchain, and the minimum R version follows from that. The alternative, a downgrade to C++17, was weighed during the release review and declined, even though the only C++20 features used are concept declarations and one bit-counting intrinsic, both of which could be expressed in C++17. A user on an older toolchain cannot build the package from source. The maintainer declined the downgrade because "the concept layer's if-constexpr seams and exact-match static_asserts do real work".
Record: docs/design/core-generalization.md. Marked: mine. [dec-B70]

**Interaction heredity was to be the next arc**
Formal heredity, the rule that an interaction can enter a tree only when its main effects are present, was scheduled as the first arc after 1.0, with soft path-dependent split penalties gated out. The alternative was building it before the release candidate. It would promise a fourth interaction argument and a change of prior in the next release. The schedule has since loosened: heredity stays on the roadmap with no position in it, and the manual states a two-stage workaround meanwhile. See also: [dec-B95].
Record: the root TODO file; docs/plans/pre-review-cleanup.md. Marked: superseded by dec-B95. [dec-B71]

**A cut-only rule draw, adopted at 0.16**
A proposal that redraws only the cut point of a rule, exactly rather than by Metropolis, was adopted at a mixture weight of 0.16 taken out of the change move's share, to land after the first release. The alternative was leaving the kernel as it is. Adopting it would commit the next release to a change in the default draws and a full re-record of the recorded baselines. The maintainer later recorded that the 0.16 came from an agent's suggestion and was approved without a rationale, so the weight ships at zero and the mixing research that resumes after the merge sets it. See also: [dec-B96].
Record: docs/design/nog-gibbs.md. Marked: superseded by dec-B96; the 0.16 is agent-chosen. [dec-B72]

**Bit-identity was to be a toggle**
The draw path was to vectorize by default, with a switch forcing the scalar kernel as a bit-identical reference, the two agreeing only in distribution. The alternative it replaced was fixing the vector lane count so that one kernel served both purposes. It was never built, and the rule that shipped instead keeps a single bit-identical kernel, so no reference-versus-fast split exists and the measured win stays unclaimed. A later ruling restates the shape as a build flag selecting a scalar reference build, and a measurement after that declined the vectorized kernels outright. See also: [dec-A42], [dec-B90], [dec-B113].
Record: docs/plans/x86-simd-plan.md. Marked: not mine; superseded by dec-B90. [dec-B73]

**A multi-level factor response is refused**
A factor response with three or more levels is refused, with a message pointing at the multinomial family on the modern front door, where 0.9-x silently coded the levels as 0, 1, 2 and fit them as numbers. The alternative was continuing that silent numeric coding. One package on CRAN has an example and tests that break this way, and its author has to adjust them. The maintainer accepted that cost: "we'll have to contact a package author who imports dbarts and have them adjust their example or test". The refusal now lives on the BayesTree-style door and names both remedies. See also: [dec-B83].
Record: this register; docs/plans/archive/cran-readiness.md records the break. Marked: superseded by dec-B83. [dec-B74]

**bart becomes the modern front door**
The name bart now carries the modern interface and defaults, bart2 becomes an alias for it for one release, and the BayesTree-style names and defaults move to a new function, which ships as bartBT. For one release bart detects a call written in the legacy style and forwards it to that function with a once-per-session warning naming it. The alternative was keeping the two entry points permanently different. Every 0.9-x script that calls bart sees that warning, reverse dependencies that call only bart keep working, and the package carries one release of a bart that answers to two vocabularies. The maintainer: "make bart our ideal front-door, redirect bart2 to bart, and move the BayesTree-similar front-door to a new function", with a shim version logged. This reopens the earlier rulings on removals, on factor predictors, on when the signature locks, on the name of the supplied residual-scale estimate, on family vocabulary and on the multi-level factor response, each for the new door. See also: [dec-B05], [dec-A02], [dec-A06], [dec-B54], [dec-A66], [dec-A25], [dec-B74].
Record: this register. Marked: mine. [dec-B75]

**A removed function stays as a tombstone**
A function that is removed stays exported for one release with a body that errors and names its successor, so the removed variable-intercept fit points at stan4bart. The sampler's thread start and stop methods stay as no-ops, and a renamed argument is handled at discretion, with the old seed argument accepted for one release under a warning. The alternative was deleting the names outright. A user who calls a removed function gets a clear error naming where to go, at the cost of a dozen stub lines carried for one release. The maintainer: "tombstones for rbart_vi makes sense, the thread methods can be no-ops ... The rngSeed name error is obscure and I don't care how it is handled". The cross-validation entry's control argument is settled separately. See also: [dec-A02], [dec-A03], [dec-A05].
Record: this register. Marked: mine. [dec-B76]

**Cross-validation keeps its flat redesign**
The cross-validation entry keeps its redesign: flat arguments, an R worker cluster, and no chain carried from one fold to the next. Tombstones cover the removed control argument and the three-element burn-in vector, and each grid cell gets its own deterministic random stream, so a seed reproduces the same result at any number of threads. Restoring warm starts across folds is considered only if the timing against 0.9-34 turns out bad, and would then be opt-in with the leakage stated. A user pays a fresh burn-in per fold and gets results that do not depend on the thread count. Three options were put and the maintainer said "go ahead with your recommendation". The control tombstone has since been reversed, since the entry takes a control object again. See also: [dec-A05], [dec-B116].
Record: this register. Marked: mine; the control tombstone is reversed by dec-B116, which returns control = dbartsControl() to xbart. [dec-B77]

**Factors split on subsets of their levels**
On the modern front door a factor predictor stays a single column and splits on subsets of its levels, while an ordered factor splits at a threshold on the level codes; the BayesTree-style function keeps expanding a factor into indicator columns. The alternative was indicator expansion everywhere. Every fit with a factor predictor differs from 0.9-x, and variable-importance counts are reported per factor rather than per level. The manual states the difference from BayesTree and why. Weighting a factor column's split mass by its number of levels is a research item for after the release. The maintainer chose "Subset splits as the modern default."
Record: this register. Marked: mine. [dec-B78]

**Argument names lock at the CRAN push**
The front door's argument names lock when 1.0-0 is pushed to CRAN, not before; until then they may change within reason. The alternative was locking them at 1.0-0 on the branch. Nothing a user sees changes before the push, and after it a rename needs the tombstone treatment. The front door takes multi-forest models as formula terms and never as an argument listing forests. The maintainer: "They lock now-ish. We should feel free to play with them, but they'll lock for certain when we push 1.0 to CRAN".
Record: this register. Marked: mine. [dec-B79]

**The supplied residual-scale estimate is sigest**
The estimate of the residual standard deviation supplied when a fit or a sampler is created is spelled sigest on every entry point, the sampler constructor included. The alternative was leaving the constructor's older spelling in place. A user who wrote the old name gets a warning for one release. The setter that sets the parameter itself on a live sampler is untouched, which was the maintainer's condition: "we also can set sigma itself with the sampler so we want to make sure we're just renaming the estimate supplied at creation".
Record: this register. Marked: mine. [dec-B80]

**Two family vocabularies, one mapping table**
The front door has its own list of family tokens, which the bridge resolves into the families the engine actually supports, while the engine's list stays internal. The alternative was one shared list. A user reads one mapping table in the manual, and that table is a maintenance item for each new family. The maintainer: "it makes sense to have two lists, one for the front door that is massaged into just those the engine supports at the bridge". The clause allowing the BayesTree-style function to take more of the front-door tokens is narrowed to none by the later ruling that makes that door a strict compatibility mode. See also: [dec-B83].
Record: this register. Marked: mine. [dec-B81]

**One token for the hurdle model**
The hurdle family tokens come off the sampler constructor's list, and the front door intercepts a hurdle request before forwarding it, so the package carries one token per model, the log-normal hurdle. The alternatives were leaving the tokens on the constructor and keeping both spellings of the name. A caller of the dropped two-part alias gets a tombstone for one release. The shipped hurdle stays an exact composition of two independent samplers driven from R; a hurdle with shared trees or correlated leaves would be a new engine family after the release.
Record: this register. Marked: mine. [dec-B82]

**The legacy door refuses and names two remedies**
The BayesTree-style function refuses a factor response with three or more levels, and its message names both remedies: the modern front door for a multinomial fit, or an explicit integer coding to get the old numeric behaviour. The alternative was keeping the old silent coding on that door. One CRAN package's author changes one line. The maintainer: "Refuse, with the two-remedy message. And I guess if we do that, we're being strict about BayesTree compatibility mode and not adding too many new features", which makes that door a strict compatibility mode gaining no new features. What to do with its existing extras beyond BayesTree, the logistic and accelerated-failure-time families and the chain and thread arguments, is settled when the rename is planned. See also: [dec-B74], [dec-B81].
Record: this register. Marked: mine. [dec-B83]

**The C header becomes pure C**
Before the release the shipped header drops the four entries that take or return R objects: sampler creation, storing state, restoring state and reading trees. A compiled consumer creates the sampler through the package's R interface and takes the handle out of the sampler object's external pointer, with R owning the lifetime, while state and trees go through the R methods; an opaque C state blob is added only if a consumer needs to store and restore from C. The alternative was keeping the R-typed entries, which the maintainer judged "a step in the wrong direction" for the goal of a non-R host, ruling "Yes, before release." The cost is one more pre-release header change: the identity hash is recomputed, stan4bart changes at three sites and treatSens at two, and the sampler factory's null-pointer refusal path disappears with the entry that carried it. The header stays the only promised interface, with no promise about the C++ one. See also: [dec-A29], [dec-A20].
Record: this register. Marked: mine. [dec-B84]

**A plain-C path to a non-R host, later**
After the release the bridge gains an internal plain-C specification that the parsing of R objects fills and the engine is built from, with the family as one of its fields, and the engine's five remaining R touchpoints, three density functions, printing and error reporting, go through hooks the host installs, so that draws stay bit for bit the same. A creation entry taking that specification and a host-neutral error contract are added in a later 1.x release once the structure settles. The alternative was doing this before the release, or not at all. Nothing a user sees changes; the cost is engine edits after the release. The maintainer adopted the staged shape alongside the header decision. See also: [dec-A50].
Record: this register. Marked: mine. [dec-B85]

**The header ships what a consumer calls**
The shipped header carries the entries a consumer actually calls plus the small queries a host needs to size its buffers. The multi-forest block, the predictor mutation pair, active rows, case weights, test predictors and their offset, and thinning all wait for a consumer or a named host design, while the latent reader stays in and a per-draw callback later shipped on named consumers. The alternative was shipping the whole block under the general rule that an item ships unless nothing valuable could follow from it; for the header alone that rule narrows to shipping when a consumer or a named host design calls for it, because an entry there cannot be taken back. Additions are batched into planned 1.x releases, so a consumer feature that needs a new entry waits for a dbarts release, or ships in lockstep from GitHub as stan4bart does today. The maintainer: "you can use your recommendation, but ... having to upgrade dbarts to ship a new feature in another package just for ABI reasons is a bit of a pain". See also: [dec-B47], [dec-B48], [dec-B61], [dec-B114].
Record: this register. Marked: mine. [dec-B86]

**The sampler owns every buffer it reads**
Every setter on the C interface copies into buffers the sampler allocated once at creation; the sampler never retains a caller's pointer and never allocates afterwards. The alternative was borrowing the caller's memory. Predictors are still re-encoded into cut codes as before, so the owned copies are the response, the offset, the test offset and the weights. The maintainer approved it with a standing concern attached: "proceed but know that I want to be able to fit larger datasets where possible, and we've had memory issues in the past", so the footprint at large n is audited rather than assumed. See also: [dec-A33].
Record: this register. Marked: mine. [dec-B87]

**Multi-worker runs stop sleeping**
The multi-worker run loop waits on chain completion and wakes as soon as the last chain finishes, keeping its hundred-millisecond timeout only for polling interrupts and flushing progress. The alternative was the fixed sleep it replaces. A user of a multi-chain run gets the result without waiting out a final sleep interval, and no draw changes. The maintainer: "yes, we can replace the fixed sleep". It lands before the release.
Record: this register. Marked: mine. [dec-B88]

**Within-chain threading was to return**
Within-chain threading was to come back before the release as an explicit opt-in revived from the archived prototype, because its absence is a regression against 0.9-34 and regressions are not deferred past the release. The thread count was to be made honest at once: a warning when the thread budget exceeds the number of chains, and a manual note that within-chain threading had been measured. The alternative was a post-release opt-in, which the maintainer revised after noting "memory is at a bit of a premium right now while cores are getting cheaper", since threads cost no memory while chains do. Re-measurement then closed it: the opt-in is archived rather than revived, and only the honest thread count shipped. See also: [dec-B115].
Record: this register. Marked: mine; superseded by dec-B115. [dec-B89]

**A reference build for bit-identity**
The shipped build was to vectorize the draw path while a build flag selected the scalar fixed-order kernel as a bit-identical reference for development, with the bitwise, exact-posterior and seed-locked snapshot checks run on the reference build and the test suite, the sanitizers and statistical equivalence run on the shipped one, plus a deterministic unit test that the vector and scalar node sums agree within tolerance. The alternative was shipping the scalar kernel alone, which is a regression against 0.9-34's vectorized sums. The cost is two continuous-integration builds, snapshot tests that skip on the shipped build, and about 3 to 4 percent of run time regained. The maintainer: "It makes sense to me to have a bit-identical path using compilation flags that can be used during development, with statistical gates before releases and whenever a set of changes justify it", adding "I'm getting a little tired of all the things we are deferring to post-release, particularly when they are regressions". Measurement later took the vectorization out, so the reference flag is inert and the two builds are identical. See also: [dec-B113].
Record: this register. Marked: mine; dec-B113 measured the vectorization out, so at the tip the reference-build flag is inert and the two builds are identical. [dec-B90]

**Seven fixed engine limits get measured and exposed**
Seven constants hard-coded in the engine are documented as limits with their origin, measured wherever measurement is possible, and the ones that matter become control settings: the size above which categorical splits stop being exactly enumerated, the cutoffs at which test fits and prediction go parallel, and any others the measurements show to bind. The alternative was leaving them fixed and undocumented. A user can set the ones that bind and read what the others mean. The control object now carries four of them: the categorical enumeration cap, the two parallel cutoffs and the sparse density threshold. The maintainer: "Document all seven, measure where it is possible so we can have sensible defaults and recommended values, and expose those that matter ... They seem arbitrary and I'm extremely worried that they're artificially limiting".
Record: this register. Marked: mine. [dec-B91]

**The engine's generic axes are kept but doubted**
The engine keeps its current shape for the release, with the leaf model fixed at compile time, the response family dispatched through a virtual interface and the code split across two translation units, but the shape is treated as provisional. An independent review of what the engine's generic axes should be, argued from the shape of the model space rather than from the current code, is written and set aside for the maintainer to read later, making it a research item with measurements rather than work now. The alternative was restructuring before the release. Nothing a user sees changes, a post-release refactor stays possible, and the parallel build proposed as the compile-latency mitigation is not in the build configuration, whose sublibrary loop is serial. The maintainer: "These were abstractions chosen very early on by just the general shape of the problem ... I wouldn't mind having a deep dive done on what the generics should be. So keep it, but be skeptical of it".
Record: this register. Marked: mine. [dec-B92]

**predict's thread argument does something**
The thread count passed to predict drives a real fan-out over chains and draws, and the size below which it stays serial is calibrated in the constants audit rather than left fixed. The alternative was leaving the argument inert. A user predicting from a large fit with several threads gets the answer sooner. Three options were put and the maintainer said "Use your recommendation". The cost is one small independent code path to run under the sanitizers.
Record: this register. Marked: mine. [dec-B93]

**Chains still start from the prior**
A chain starts from trees drawn from the prior by default, and the alternative start, grown by sweeping before the chain proper, stays opt-in through an argument for the number of growing sweeps. That choice rests on measurements showing the grown start's advantage plateaus in noise-heavy and large-n designs. A user sees no change from 0.9-x, and the manual states which families refuse the option. Three options were put and the maintainer said "You can use your recommendation".
Record: this register. Marked: mine. [dec-B94]

**Heredity happens, but not on a schedule**
Formal heredity stays on the post-release roadmap with no position in it: it will happen, and when is open, while soft path-dependent split penalties stay gated out. The alternative was fixing it as the first post-release arc. Meanwhile the manual states a two-stage workaround a user can run themselves. The maintainer: "it doesn't have to be the first post-1.0 item. I actually don't care when it happens, just that it does". See also: [dec-B71].
Record: this register. Marked: mine. [dec-B95]

**The exact rule draw ships at weight zero**
The exact rule draw stays in the tree with a default weight of zero: the shipped kernel draws the variable and the cut jointly, and the cut-only variant is reachable only in a private build. The alternative was shipping it at the weight of 0.16 adopted earlier, which the maintainer records as an agent's suggestion taken without a rationale and not a maintainer choice. A user sees no change in the draws for the release candidate, and a measured mixing gain stays unavailable until the weight is set. The mixing research that resumes after the merge sets the default weight with a rationale, before 1.0-0 if time allows and otherwise after, and the manual documents the kernel and its measured cost. The maintainer: "I approved of 0.16 without much of a rationale ... And I don't care where it happens." See also: [dec-B72].
Record: this register. Marked: mine. [dec-B96]

**The survival interface is finished before release**
The formula interface accepts a survival response on its left-hand side, the subset argument is honoured for both survival families, and the discrete-time hazard family gains a test path that expands held-out subjects into their time intervals. All of it is in the R layer, with no engine change. The alternative was deferring some or all of it past the release. A user can fit and predict survival models through the ordinary formula interface. Three options were put and the maintainer said "Use your recommendation".
Record: this register. Marked: mine. [dec-B97]

**Family settings ride family objects**
The family argument accepts a token or a call, so a hazard family carries its time breaks and its row cap and a negative-binomial family carries its dispersion, following base R's idiom for generalized linear models. Those settings leave the front door's own argument list with tombstones for one release, and the consolidation is reopened so that any remaining family-only or feature-only argument moves onto its object before the CRAN push. The alternative was keeping them as front-door arguments. A user writes the setting inside the family call rather than alongside it, and the ten-million-row cap joins the constants audit. The maintainer: "I thought we consolidated the interface? If we haven't, we should. And yes, family objects."
Record: this register. Marked: mine. [dec-B98]

**The posterior package is removed, not suggested**
The package neither depends on nor suggests the posterior package: the fit summary computes split R-hat and effective sample size itself, so the shape of that output is fixed by dbarts, and the draws-conversion methods are replaced by a base-R extractor returning an iterations by chains by variables array with dimnames, which that package's own constructors accept if a user has it installed. The alternative was keeping it as a suggested dependency. Matrix and survival stay suggested behind their features, both being recommended packages that ship with R. The cost is a short diagnostics implementation of the package's own, plus tests and a manual page rewritten without the dependency. The maintainer: "I hate dependencies ... Absolutely not. I want it removed. Whatever we were doing with it, find another way."
Record: this register. Marked: mine. [dec-B99]

**Sparse columns need no formula marker**
A sparse matrix or a sparse factor is assigned into the data frame as a column and named in the formula like any other predictor, with no marker; ingestion detects such columns by class and lifts them around the model frame, and wide factors are stored sparse by the engine without the user asking. The alternative was a formula term that marks a column as sparse. A user writes the ordinary formula and gets the sparse storage. The maintainer left the interface open and asked "if there is a sparse matrix in a data frame, why do we need to identify it as sparse in the formula?" The work is in the R layer only, before the release.
Record: this register. Marked: mine. [dec-B100]

**Residual laws fold into family objects**
The residual law folds into the family object, either as a family variant or as a setting on the gaussian family, retiring the separate argument for it and both of its vocabulary bundles from the front door. The alternative was keeping that argument. If the consolidation leaves it in place, the constructor list for residual laws is exported so that it matches the exported list for priors and a user can see the choices. Three options were put and the maintainer said "Use your recommendation".
Record: this register. Marked: mine. [dec-B101]

**Cross-validation can leave k modelled**
The cross-validation grid for the leaf-scale parameter k accepts both fixed values and an entry meaning leave k modelled under the family's hyperprior, so a user can cross-validate the other parameters with k sampled, or fix it, in one call. The alternative was a grid of fixed values only. The fixed grid stays the ordinary case, and the cost is one non-numeric grid entry and a path that hands the hyperprior through. The maintainer: "originally I wanted to be able to leave k modeled and crossvalidate over the other parameters, as well as be able to set it fixed." See also: [dec-B12].
Record: this register. Marked: mine. [dec-B102]

**Counts must be whole; weights need not be**
A fractional value is refused by name where an argument is a count, a column index or a group selector, and never where it is a weight: weights may be fractional wherever the family supports real weights. The alternative was extending the refusal to weights, which the maintainer suspected was a conflation of the two. A user whose count is computed needs a rounding call, and the upgrading section of the manual says that counts must be whole numbers. The cases where a weight must be an integer come from the negative-binomial and binary-count rulings, not from this one. The maintainer: "I wasn't sure if you were conflating counts and weights. You can refuse, as ruled." See also: [dec-B08], [dec-B15].
Record: this register. Marked: mine. [dec-B103]

**A saved fit needs its state stored first**
Trees are still not kept by default and sampler state is stored only when a user asks for it, since the stored state runs about 2.8 times the size of the training predictions. The alternative was storing it automatically. The manual entry for the tree-keeping argument says in its first sentence that a saved fit needs the state stored, and the refusal at predict names the call that does it. A three-valued setting for that argument is the fallback if the support burden proves real. The maintainer: "I was originally about memory size but wanted to double check ... I think that was the right call." See also: [dec-B33].
Record: this register. Marked: mine. [dec-B104]

**The multilevel fit moves to stan4bart**
The variable-intercept fit stays out of dbarts, with a tombstone naming stan4bart and one line saying that the group-spread prior differs there, so results move. The alternative, an R-only grouped intercept that drives the plain sampler by setting the offset each sweep with no engine work, is rejected. The condition on the removal is now met: stan4bart clears the group-spread mixing bar at its package defaults, since its closed-form slice move along the scale ridge landed on 2026-09-13, taking lag-one autocorrelation on the reference design from 0.96 to 0.09 and effective draws per 1000 from about 14 to over 500, with posterior means unmoved and wall time within one percent. The refreshed WALNUTS and its local refresh patch are recorded in stan4bart's own vendoring note. What remains in stan4bart is the shared level of the group intercepts, a separate failure outside the bar; that needs a level-shift C API entry, recorded in TODO as post-1.0 work. A user who wants a multilevel BART fit installs stan4bart, and bartCause builds its grouped route through that package before the release.
Record: this register. Marked: mine. [dec-B105]

**The binary k default waits on a study**
The relabelled chi constructor stands, and a chi prior with 1.5 degrees of freedom and scale 2 stays the binary default for now. The study behind it rejects an infinite scale clearly, prefers scale 2 over 5 only narrowly, and says nothing about the degrees of freedom, so a further evaluation is scheduled: a checked-in harness, the degrees of freedom varied, and a much wider range of cases than four simulated data-generating processes, including real datasets and the regimes users actually bring. The alternative was treating the recommendation as final. A user of a binary fit may see the default prior change once that study reports. The maintainer: "Schedule some additional evaluation. I'm also skeptical when these numbers run of where the evaluation was done - we should make sure that it works in a wide variety of use cases, not just 4 data generating processes". The timing has since hardened, so the study now completes before a release candidate. See also: [dec-B118].
Record: this register. Marked: mine; the timing is superseded by dec-B118, the default and the study design stand. [dec-B106]

**The drop-in replacement claim is dropped**
The package description says that dbarts provides a BayesTree-compatible interface, dropping the claim that it is a drop-in replacement. The alternative was keeping the older claim. Three removed build-configuration options stay as stubs for one release that stop with a message naming the removal, so a user who passes one is told rather than ignored. The maintainer: "the value in having a drop-in replacement for BayesTree was that it was slow and not maintained at the time. That's no longer true."
Record: this register. Marked: mine. [dec-B107]

**The front door gains na.action**
The front door takes a standard na.action argument whose default is a package function that drops rows with a missing response, keeps missing predictors for the trees to handle, and records the dropped rows so that training-set fitted values pad back out to the length of the user's data. The base functions keep their usual meanings: one drops any row with a missing value anywhere, one errors on any missing value, and one keeps everything and then trips the response check. The alternative was a bespoke argument, or none at all. A user with missing responses gets a fit on the complete cases and fitted values aligned to their data. The maintainer, asked whether missing responses and missing predictors have to be distinguished, said "Yes, go ahead."; the work lands in the family-objects consolidation pass. See also: [dec-A70].
Record: this register. Marked: mine. [dec-B108]

**An old saved fit is recognized and refused**
A saved 0.9-x fit is recognized by the absence of the state format field and refused with a message naming the version and pointing at a refit; there is no conversion. The alternative was a converter. A user who loads an old fit and predicts gets a message that says why, rather than an obscure failure. The maintainer: "Won't it just be obvious if it doesn't have a format field? But yes, recognize and say so." The few lines involved are removable once such fits are unlikely. See also: [dec-B31].
Record: this register. Marked: mine. [dec-B109]

**A degenerate GP fit reports itself**
A fit with Gaussian-process leaves counts the leaf evaluations that fell back to constant leaves and warns when that share is high, so a degenerate fit says so whatever the cause. The alternatives were guidance in the manual alone, or a hard refusal. A user who leaves the tree count at the ensemble default sees a warning rather than a quietly constant fit; the 10 to 25 tree guidance stays in the manual and the maximum leaf size joins the constants audit. Three options were put and the maintainer said "Use your recommendation." See also: [dec-B41].
Record: this register. Marked: mine. [dec-B110]

**The version pair guards the interface**
The exact hash flag is dropped from the sister packages at the coordinated merge and stays off; the major and minor version pair is the guard, and the package's continuous integration asserts that a changed interface hash arrives with a bump of the minor version, so the discipline the pair depends on is enforced. The alternative was keeping the strict check on for consumers that ship in lockstep. Nothing a user sees changes. The assertion is dormant until the 1.0-0 tag exists, since the check skips unless a 1.x release tag is present and the repository's tags are 0.8-7 and a pre-CRAN rebase marker. Three options were put and the maintainer said "Use your recommendation." See also: [dec-B58].
Record: this register. Marked: mine. [dec-B111]

**The C interface keeps a destroy entry**
The header keeps its destroy entry even under the rule that the header is pure C: an R-hosted consumer lets the garbage collector free the sampler, while a consumer written in C keeps a way to invalidate its own object. Destroying releases the engine sampler behind the handle and marks the R object's pointer dead; the R object's dead-pointer path then re-creates from a stored state or refuses, and a second destroy does nothing. The alternative, which the plan recommended, was removing the entry. The maintainer: "I assume that this ruling only applies to R consumers - they should let the objects be garbage collected. A C consumer should have a path to invalidate a C object." See also: [dec-B84].
Record: docs/plans/pure-c-header.md. Marked: mine. [dec-B112]

**Vectorized kernels are declined; the fused pass ships**
The vectorized sufficient-statistic kernels do not ship: they measured under one percent of a weighted fit on arm64 and within noise on x86 with AVX2, because the hot kernel is limited by gathering scattered values rather than by arithmetic. The code is archived with its record and the reference-build flag stays inert, so the two builds are identical. Instead the fused pass that computes residuals and their sums in one sweep is extended to weighted families, which measured 26 to 29 percent faster at a hundred thousand observations. The maintainer set the order of measurement, "Why don't we see how it works on weights first before deciding" and "Try the AVX2.", and read the result as shipping the fused pass and archiving the kernels. See also: [dec-B90].
Record: docs/plans/engine-performance.md. Marked: mine. [dec-B113]

**A callback fires on every saved draw**
A host can register a callback that fires once per saved draw per chain, on that chain's worker thread, over const pointers into the engine's own buffers; the shipped header carries the draw structure, the callback type and the setter for it. A new logical controlling whether the per-observation channels are kept is set to false automatically when a callback is supplied, and a nonzero return from the callback aborts the run through a shared cancel flag. The alternatives, settled fork by fork, included a narrower per-channel control and a callback with no return value. A user who supplies a callback silently stops getting test predictions, which nothing drops today; a mistake inside a callback crashes the session with no condition to catch; and an interrupt cannot land while a call is running. The maintainer asked for the work, ruled "Ship the header entry and its two types now.", approved the keep-fits control over every per-observation channel, and picked the integer return with its cancel plumbing, leaving the remaining forks to the recommendations. The costs are one irreversible header entry, a recomputed interface hash, and a control frozen at 1.0. See also: [dec-B62], [dec-B86].
Record: docs/plans/per-draw-callbacks.md; docs/design/per-draw-callbacks.md. Marked: mine. [dec-B114]

**Within-chain threading closes; the thread count gets honest**
Within-chain threading does not ship: re-measured on the current engine it reached at best 1.03 times the speed anywhere, with four workers at a million observations, and lost about 5 percent at two workers and 15 percent at eight at a hundred thousand, so the prototype is archived rather than revived. The thread count keeps its own meaning as a total thread budget distinct from the number of chains, rather than collapsing onto it, and its default becomes the smaller of a guess at the core count and the number of chains. A budget above the chain count warns once per fit, naming both numbers, that tree sampling uses at most one thread per chain and that the excess reaches only the test-fit pool and prediction's fan-out. The alternative was reviving the opt-in the earlier ruling called for. No single-chain speed lever ships before the release, the prototype's correctness half, byte-identical draws across worker counts, is banked rather than revived, and a caller who asks for more threads than chains gets a warning instead of a faster sweep. The maintainer: "Close it (archive it?). I dislike treating n.threads as n.chains and would like to keep that distinction, but we can simple warn that excess threads won't be used and default to setting n.threads to n.chains." See also: [dec-B89].
Record: docs/plans/engine-performance.md; docs/design/within-chain-threading.md. Marked: mine. [dec-B115]

**Prior settings move onto prior objects**
The front door sheds seven scalar arguments that duplicated slots on the prior objects: the tree prior's power and base and the split probabilities move onto the tree-prior constructors, the leaf prior's scale onto the normal leaf prior, and the residual prior's degrees of freedom and quantile onto the chi-squared residual prior. The mixture over tree moves becomes a control slot that the bridge reads, while the leaf-scale parameter k stays as the one promoted scalar and the supplied residual-scale estimate stays on the data side. The alternatives were leaving these reachable only through the sampler constructor, and promoting four of them to arguments with the tree prior carrying the cap. The retired names ride the dots through the existing tombstone mechanism and expire with the rest at 1.1-0, so every example and consumer test written against them now warns; the front door and the cross-validation entry both gain a control argument under one precedence rule, a supplied flat argument beating the control slot, and a control carrying fit-state attributes is refused. That reverses the retired control argument on cross-validation before its removal ever ships, and moving the residual prior onto the gaussian family is a recorded follow-on. The maintainer: "Great, then let's do it." See also: [dec-B77], [dec-B10].
Record: docs/plans/front-door-formals.md. Marked: mine. [dec-B116]

**Updating a predictor repartitions in place**
After a predictor update the engine partitions a dense root in place instead of rewriting the identity order, because the span it receives is already partitioned under the live rules, so the scan writes nothing where the rewrite wrote the whole span. Predictor updates run 19 to 37 percent faster at no cost to sampling, whose root dispatch is untouched, and that is the package's distinguishing path, a sampler driven inside a larger model. The member order differs, so leaf sums reassociate in the last bit and every predictor-update scenario has to be re-recorded; the sparse layout stays non-bitwise, since making it bitwise would cost sparse-tier fits 3 to 7 percent for a property no user observes. The alternative was landing it after the release to avoid that re-record; both directions were measured at the maintainer's request and replicated independently on two machines. The maintainer: "Land it before release."
Record: docs/plans/setpredictor-partition.md; the root TODO's engine constants audit item. Marked: mine. [dec-B117]

**The binary prior study gates the release candidate**
The binary node hyperprior evaluation completes before a release candidate is declared, rather than before 1.0-0 if time allows. Only the timing moved from the earlier ruling: the study design and the shipped default stood throughout. That study ran and reported on 2026-09-14 (docs/plans/binary-hyperprior.md: 28 priors, three chain lengths, 162 simulated cells and 22 real datasets including sixteen from UCI), recommending that chi(1.5, 2) stay, and the maintainer confirmed that default the same day on the study's evidence. The maintainer, on 2026-09-11: "I'd like it to be considered pre-release, and even pre-release candidate". See also: [dec-A07], [dec-B106].
Record: this register. Marked: mine. [dec-B118]

**Errors unwind instead of jumping**
A C++ exception thrown inside a callback is caught at the call and rethrown only after the callback's own frame has returned, the jump being made under R's unwind protection so that it unwinds through that frame rather than across it; an exception the engine itself raises travels the same path, and either becomes an R error only at the bridge entry point, once the unwind has run. No raw error call and no long jump leaves engine or callback code. No alternative was weighed. A host written in C++ still may not jump to a landing pad of its own from inside a callback: only raising an R error or throwing is safe. The same change fixes three sites that leaked heap memory, since every path now unwinds through a frame instead of jumping past it. See also: [dec-A31].
Record: inst/include/dbarts/dbarts.h; src/R_interface_bartcore_common.hpp. Marked: not put to the maintainer; landed 2026-09-13. [dec-B119]

**A wider comparison against 0.9-34**
A 26-scenario statistical comparison against an installed dbarts 0.9-34 widens the cross-engine record past the nine-scenario snapshot taken at the cutover: 22 of the 26 agree at the rate the null predicts, and the other four are diagnosed as decided engine changes rather than regressions, namely the degrees of freedom and veto rank used for the residual scale under positive weights, the absence of a chain carried across cross-validation folds, and the repaired acceptance ratio in the change move. No alternative was weighed. It is a measurement run by hand, not a gate: it needs a hand-installed 0.9-34 library, while the per-push equivalence check still compares only against baselines recorded on this branch. A future engine change could therefore separate a scenario again and go unnoticed between manual runs. See also: [dec-A61].
Record: docs/plans/classic-compare.md. Marked: not put to the maintainer; landed 2026-09-13. [dec-B120]

**The drawn variance surface is readable**
A sampler fitting a heteroscedastic model is to expose the variance surface it has drawn through a current-state accessor on the sampler object, returning what a run's variance and test-variance channels report and NULL when the sampler is homoscedastic, and the refusal that blocks prior-predictive draws for such samplers is lifted with it. The ruling is from 2026-09-13, and the accessor and the prior-predictive draw landed the same day. The two alternatives not taken were rebuilding the surface in R from the reported trees, which duplicates engine arithmetic inside the check meant to test it, and deferring both arms of the simulation-based calibration check, which would ship with heteroscedastic calibration unchecked. A user can read the variance surface a heteroscedastic fit drew and can draw from its prior predictive. Three options were put and the maintainer said "Use option 1".
Record: docs/design/aft-status-setter.md landing note. Marked: mine. [dec-B121]

**The variance forest's leaf prior recalibrates on a rescaling swap**
A setResponse or setOffset call with updateScale = TRUE on a heteroscedastic sampler recalibrates the variance forest's leaf prior to the rescaled response, reusing the recalibration the engine already performs on a setModel swap. The leaf prior is calibrated once, at creation, from the residual-prior triple on the working scale, and nothing had rescaled it when a later swap moved that scale: the prior stayed pinned to the old scale while the drawn surface adapted to the new one over the following sweeps, a defect no default configuration hits, since updateScale defaults to FALSE. The alternatives were refusing the swap outright for a heteroscedastic sampler under updateScale = TRUE, and documenting the staleness without fixing it. A user swapping in a rescaled response under updateScale = TRUE would otherwise see the leaf prior and the drawn surface disagree by a growing amount until the surface's own adaptation caught up. The maintainer ruled the recalibration in on 2026-09-14, gated on an identity test: a sampler swapped to a rescaled response must match one built fresh on it, with refusal as the fallback had the test failed. The test held, bitwise in the leaf prior and to rounding in the draws, and the recalibration shipped for setResponse and setOffset. The maintainer then extended the ruling the same day to the whole-data swap, setData, which had left sigma pinned to the creation-time scale on a heteroscedastic sampler, a defect a user saw as a sigma channel scaled by the ratio of the old and new response ranges; setData now re-anchors the variance forest by the same route and matches a fresh sampler under the same test, and an equivalence scenario exercises it.
Record: docs/plans/archive/variance-forest-mutation-routing.md, "Doors held open". Marked: mine. [dec-B122]

**The residual prior has one home, and sigest belongs to the chisq one**
The residual scale's prior is reached only through the family object that draws that scale - gaussian(sigma = chisq(df, quant)) or gaussian(sigma = fixed(value)), and likewise on student, aft and hurdle.lognormal - on every entry point. The flat resid.prior argument that the sampler constructors, the specification builder and the crossvalidator had kept beside the family's own setting is retired into the same one-release tombstone the fitter's already carried: it warns once per session, the value is still used, and it is gone at 1.1-0. The alternative was the surface as it stood, two homes on the constructors with a silent precedence rule that let the flat argument beat the family's setting; it was rejected because a caller reading a call cannot tell which of the two spellings the fit used. A call that writes the prior both ways is now refused where the two disagree, naming both spellings and saying which to delete, and is accepted in silence where they say the same thing. Beside the prior, sigest stays a flat argument: it is the residual-scale estimate a chisq prior's quantile is calibrated against, so it stands beside chisq and is refused beside fixed, which is the residual scale itself and would have the estimate overwritten by its square root. A user who wrote resid.prior on a sampler constructor sees one warning and the same fit; a user who wrote the prior twice and meant two different things gets an error instead of a silent winner; a user who wrote sigest beside a fixed scale gets an error instead of an argument that did nothing. The maintainer ruled both on 2026-09-14.
Record: docs/plans/front-door-formals.md, "the residual prior's one home" landing note. Marked: mine. [dec-B123]

## C. Agent-made decisions with no identified cost

**C entry points register under full names**
The compiled entry points are registered under their full symbol names rather than the unprefixed method names the released package used. No alternative was weighed. Nothing is visible to an R user; the naming is what makes an old compiled consumer fail at lookup instead of resolving a same-named symbol with a different signature. Not yet ruled on.
Record: docs/design/public-surface.md, which calls handle-prefixed names defensible. Marked: not mine. [dec-C01]

**The compiled interface is checked at build**
The package asserts its own side of the compiled interface when it builds: a signature token, a recomputed hash, and the offsets, sizes and alignments of the three structs the interface exposes. No alternative was weighed. Nothing is visible to a user; drift inside the package is caught at compile time. Not yet ruled on.
Record: code only, in the flat C interface implementation. Marked: not mine. [dec-C02]

**Internal call registrations get a new prefix**
The registrations R uses to reach compiled code are renamed to an engine-specific prefix, 58 of 65, with seven utility entries keeping the plain dbarts prefix and the dynamic-library line unchanged. No alternative was weighed. Nothing is visible to a user, and no R-visible name changes. Not yet ruled on.
Record: code only, in the registration table. Marked: not mine. [dec-C03]

**The support library no longer needs R**
The support library is host-agnostic: printing goes through function pointers the host installs, so the library carries no standard input-output and no R symbol. No alternative was weighed. Nothing is visible to a user. The maintainer adopted the same mechanism for the engine's remaining calls into R, staged for after the release. See also: [dec-B85].
Record: code only, in the support library's input-output unit. Marked: blank. [dec-C04]

**The grown warm start stays opt-in**
A sampler starts from a forest drawn from the prior, and the warm start that grows trees from the root first is available through an argument. The alternative, defaulting that warm start on, was measured against a criterion fixed in advance and killed, the measurement showing a cost in noise-heavy and large-sample designs. A user who wants it asks for it by name, and the manual says which families refuse the option. The maintainer was given three options and took this one. See also: [dec-B94].
Record: docs/design/grow-from-root-default.md, whose status records the measured kill. Marked: blank. [dec-C05]

**Two research switches ship unbuilt**
Two compile-time switches meant for research measurements ship in the source and are defined by no shipped build. No alternative was weighed. Nothing is visible to a user; both say in the code that they are for private builds only. Not yet ruled on.
Record: code only, in the moves header and the C++ test makefile. Marked: not mine. [dec-C06]

**Two developer tools ship unwired**
The stale-install detector and the snapshot regeneration script ship as local developer tools, wired into no workflow. No alternative was weighed. Nothing is visible to a user; both say so in their own headers. Not yet ruled on.
Record: code only, in the tools directory. Marked: not mine. [dec-C07]
