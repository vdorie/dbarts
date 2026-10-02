Review 3, lens: docs (tree 01dee4b4, library r3-lib)

Covered
- inst/NEWS.Rd 1.0-0 section: every UPGRADING / USER-VISIBLE / NEW FEATURES claim with a cheap probe was probed against r3-lib, and against a 0.9-34 build from `git archive main` (scratch lib r3-docs-mainlib) where the claim is about 0.9-34.
- Every man/*.Rd example extracted (tools::Rd2ex, dontrun/donttest included) and run; the three vignettes purled and run (guessNumCores pinned to 2).
- tools::codoc / undoc / checkDocFiles on the installed package: clean.
- Read and probed: bart.Rd (arguments, Value headline items), dbarts.Rd (arguments), dbartsControl.Rd, xbart.Rd, monotone.Rd, interactions.Rd, blocks.Rd, dbartsFamilies.Rd, dbartsPriors.Rd (normal/chi/details), dbarts-deprecated.Rd, na.keepPredictors.Rd, summary.bart.Rd, dbarts-package.Rd, plotTree.Rd usage, survivalProbabilities.Rd value, dbartsSampler-class.Rd Value section (run/predict/getTrees shapes, getK, getSSR, getFitsWithoutOffset, getLatents).

Not covered
- bartBT.Rd beyond binaryOffset/sigest/Saving; dbartsSampler-class.Rd arguments and method sections (98 KB) beyond the Value section; forest.Rd, dbartsForests.Rd, dbartsSpec.Rd, dbartsData.Rd, dbarts-embedding.Rd, dbartsValidateComposition.Rd, sparseFactor.Rd, varianceForest.Rd, pdbart.Rd prose, rbart.Rd (deprecated) prose; the per-family paragraphs of bart.Rd/dbarts.Rd (multinomial, ordinal, nbinom, hazard, hurdle) beyond the class/auto-detection claims; vignette prose (code only); NEWS performance and memory numbers (not reproducible here).

Findings: 2 MAJOR, 12 MINOR, 0 BLOCKER.

docs-01  MAJOR
Location: inst/NEWS.Rd (NEW FEATURES, first item and the xbart item); man/bart.Rd \item{storage}.
Claim: NEWS says `storage = "single"` is "a dbartsControl, bart, bartBT and xbart argument" giving 10 to 30 percent more speed, and that xbart "takes ... storage"; bart.Rd says storage is "a formal on both bartBT and bart". bartBT has no storage formal at all, and xbart refuses "single" unconditionally (xbart.Rd itself says so).
Probe:
  bartBT(df$a, df$y, storage = "single", ndpost = 10, nskip = 5, verbose = FALSE)
  #  Error: unused argument (storage = "single")
  "storage" %in% names(formals(bartBT))       # FALSE
  xbart(y ~ a + b + c, df, storage = "single", n.samples = 10, n.burn = c(5,5), n.reps = 2, n.threads = 1, n.trees = 5)
  #  Error: storage = "single" is not supported for this sampler
Why gates missed: codoc compares usage blocks with formals, and bartBT.Rd's usage correctly omits storage; nothing checks prose or NEWS against formals (TODO doc-freshness-guard-tightening covers neither man/ nor NEWS).
Fix: NEWS: "a dbartsControl and bart argument" and drop storage from the list of things xbart takes (or say xbart accepts only "double"); bart.Rd: delete "a formal on both bartBT and bart" (bartBT users reach it only by building their own sampler).

docs-02  MAJOR
Location: inst/NEWS.Rd USER-VISIBLE CHANGES, last item; R/utility.R warnOnce, R/tombstones.R (every tombstone), man/*.Rd.
Claim: NEWS says "Every warning the package raises has a class inheriting from dbartsWarning, documented on the help page of the function that raises it, so it can be caught or muted without matching its text." The once-per-session tombstone warnings - the ones every 0.9-34 user meets first - are plain simpleWarning, and several classes are documented nowhere; dbartsWarning itself is named on no help page.
Probe:
  withCallingHandlers(bart(x.train = x, y.train = y, ntree = 10, ndpost = 10, nskip = 5, verbose = FALSE),
    warning = function(w) print(class(w)))
  #  "simpleWarning" "warning" "condition"
  # same for: bart(..., power = 1), sigdf, proposal.probs, resid.prior (bart/dbarts/xbart), node.prior,
  # sigma (dbarts/xbart), rngSeed, seed = NA, sigest = NA, bart2(), $startThreads, $run(.., n.threads),
  # $sampleNodeParametersFromPrior
  grep over man/: dbartsExcessThreadsWarning (raised by bart/dbarts/bartBT whenever n.threads > n.chains),
  dbartsDuplicateNameWarning, and dbartsGPFallbackWarning (dbartsPriors.Rd describes the warning but not its
  class) appear on no page; "dbartsWarning" appears only in NEWS.
Why gates missed: test-tombstones.R asserts message text and once-ness, not class; no gate asserts that warning classes appear in man/.
Fix: Give warnOnce a class argument (e.g. c("dbartsDeprecatedWarning", "dbartsWarning")) used by every tombstone, document dbartsWarning and the three missing subclasses, or narrow the NEWS sentence to what is true.

docs-03  MINOR
Location: man/bart.Rd \item{keepTrees}.
Claim: "A fit saved to disk and reloaded needs its sampler's state stored first, with fit$storeState()". On a bart fit that call fails with an unhelpful error; the working call is fit$fit$storeState() (as bartBT.Rd's Saving section and bart.Rd line 463 correctly say).
Probe:
  f <- bart(y ~ ., df, n.samples = 10, n.burn = 5, n.chains = 2, n.threads = 1, verbose = FALSE, keepTrees = TRUE)
  f$storeState()        # Error: attempt to apply non-function
  f$fit$storeState(); saveRDS(f, tf); dim(predict(readRDS(tf), df[1:3, ]))   # 20 3
Why gates missed: prose only.
Fix: write \code{fit$fit$storeState()}.

docs-04  MINOR
Location: man/dbarts.Rd \item{test}.
Claim: "when test's columns cannot all be matched by name to data's, they are matched by position instead, with a warning (class dbartsPositionalArgsWarning)". When both are named and the names differ it is an error (R/data.R validateXTest); the positional warning is only for an unnamed side.
Probe:
  x <- matrix(rnorm(200), 100, dimnames = list(NULL, c("a","b"))); xt <- x[1:5, ]; colnames(xt) <- c("p","q")
  dbarts(x, rnorm(100), test = xt)
  #  Error: column names of 'test' do not match those of 'x': 'a, b' present in 'x' but not in 'test' (whose columns are 'p, q')
Why gates missed: prose only; the behaviour itself is tested.
Fix: "if both are named, test's names must cover data's (an error otherwise); if only one is named, columns are matched by position with a warning".

docs-05  MINOR
Location: inst/NEWS.Rd (USER-VISIBLE CHANGES / UPGRADING, "New data given to predict"); R/data.R validateXTest.
Claim: A user-visible change is missing from NEWS: an unnamed test or newdata matrix against a fit whose predictors have names now warns on every call (0.9-34 was silent). This includes the most common matrix-in-formula pattern, where the names were synthesized by dbarts, and the message then wrongly says the user's 'x' "had named predictors".
Probe (same script on both builds):
  x <- matrix(runif(100), 50, 2); y <- x[, 1] + rnorm(50, 0, .1)
  f <- bart(y ~ x, n.samples = 5, n.burn = 5, n.chains = 1, verbose = FALSE, keepTrees = TRUE); predict(f, x[1:4, ])
  # 0.9-34 (bart2): no warning
  # 1.0-0: Warning: 'test' is unnamed but 'x' had named predictors, matched to 'x' by position
  #        (column 1 = 'x.1', column 2 = 'x.2'); supply 'test' with column names to match by name instead
  # (every call, not once; also fires in vignette working_with_saved_trees, extract(..., newdata = x[1:5, ]))
Why gates missed: NEWS has no coverage gate against behaviour diffs.
Fix: add a NEWS line under "New data given to predict"; in the message, say "the fit's predictors are named" (and consider not warning when the names were generated from a matrix term).

docs-06  MINOR
Location: inst/NEWS.Rd UPGRADING, last item.
Claim: "names(fit) no longer lists ... on a bartBT fit without an offset binaryOffset, when the fit has none". A default binary bartBT fit (no offset given) still lists binaryOffset, as a zero vector, as 0.9-34's bart did; what changed is the bart2-style door, where 0.9-34 listed binaryOffset = NULL.
Probe:
  f <- bartBT(x, yb, ndpost = 5, nskip = 5, verbose = FALSE); "binaryOffset" %in% names(f)   # TRUE
  str(f$binaryOffset)   # num [1:100] 0 0 0 ...
  "binaryOffset" %in% names(bart(x, yb, n.samples = 5, n.burn = 5, n.chains = 1, verbose = FALSE))   # FALSE
Why gates missed: prose only.
Fix: "or, on a binary bart fit without an offset, binaryOffset".

docs-07  MINOR
Location: inst/NEWS.Rd NEW FEATURES ("Tree moves" item).
Claim: "bart and xbart gain a control argument". xbart already had `control = dbarts::dbartsControl()` in 0.9-34 (git show main:R/xbart.R). Relatedly, man/dbarts-deprecated.Rd says "control on xbart is no longer a tombstone ... the earlier retirement (which refused the argument outright) is reversed": that retirement never shipped, so it is development history on a user page.
Probe: git show main:R/xbart.R | sed -n '/^xbart <- function/,/^)/p'   # ... resid.prior = chisq, control = dbarts::dbartsControl(), sigma = NA_real_, ...
Why gates missed: no NEWS-vs-main surface check.
Fix: "bart gains a control argument"; delete the xbart-control paragraph from dbarts-deprecated.Rd.

docs-08  MINOR
Location: man/dbartsFamilies.Rd \description, third paragraph.
Claim: "resid.prior is retired that way [accepted for one release with a once-per-session warning] on dbarts, dbartsSpec and xbart as well". dbartsSpec refuses it, as dbarts-deprecated.Rd correctly says; the two pages contradict each other.
Probe:
  dbartsSpec(dbartsData(y ~ a + b + c, df), resid.prior = chisq(5, .9))
  #  Error: unused argument 'resid.prior' passed to 'dbartsSpec'; the residual prior rides its family: write family = gaussian(sigma = )
Why gates missed: prose only.
Fix: "on dbarts and xbart as well, and dbartsSpec refuses it".

docs-09  MINOR
Location: man/pdbart.Rd examples (third block); vignettes/working_with_saved_trees.Rmd chunk fitModel.
Claim: Shipped examples call bart with BayesTree names (keepevery, ntree, keeptrees; ndpost, nskip, nchain, nthread), so they run through the tombstone forwarding to bartBT, print its deprecation warning, and stop working in 1.1-0. The vignette's own bullet list says "For bart: keepTrees = TRUE" two lines above the call that uses keeptrees.
Probe (Rd2ex + source, and knitr::purl + source):
  pdbart.R: Warning: 'keepevery' is dbarts 0.9-x's BayesTree-style 'bart' argument; that function is now 'bartBT' and this call was forwarded to it. ...
  working_with_saved_trees.R: Warning: 'ndpost' is dbarts 0.9-x's BayesTree-style 'bart' argument; ...
Why gates missed: R CMD check does not fail on warnings from examples or vignettes.
Fix: call bartBT in both (or rewrite with bart's names); give the vignette's newdata matrix column names (docs-05).

docs-10  MINOR
Location: man/dbartsPriors.Rd examples (gp block); R/bartcore.R warnOnGPFallback.
Claim: The shipped gp() example warns that its own fit is mostly not a Gaussian process, and the message says "most of this fit" when the threshold is a quarter (here 35.6 percent).
Probe: running the extracted example:
  fit.gp <- dbarts(y.gp ~ x1 + x2, df.gp, leaf.prior = gp("x1", max.leaf.size = 30), control = dbartsControl(n.trees = 10, n.chains = 1, n.threads = 1)); fit.gp$run(20, 20)
  #  Warning: 35.6% of Gaussian-process leaf evaluations fell back to a constant leaf ..., so most of this fit is not a Gaussian process; ...
Why gates missed: example warnings are not gated.
Fix: raise max.leaf.size (or the tree count) in the example so it demonstrates a working gp fit; reword to "so much of this fit" or only say "most" above one half.

docs-11  MINOR
Location: man/dbartsPriors.Rd \details; man/bart.Rd \item{k}; R/bart.R (k evaluation).
Claim: dbartsPriors says prior constructors resolve by bare name inside the fitting functions' prior arguments, including "for an argument a wrapper forwards through its \dots". bart's k (documented as taking chi(df, scale)) does not: forwarded through a wrapper's dots it fails, while family, tree.prior and leaf.prior forwarded the same way work.
Probe:
  W <- function(...) bart(..., n.samples = 20, n.burn = 5, n.chains = 1, n.threads = 1, verbose = FALSE)
  W(yb ~ ., df, k = chi(1.5, 2))
  #  Error: could not find function "chi"; outside the argument that takes it, write dbartsPriors$chi(...)
  W(y ~ ., df, leaf.prior = normal(sd = 1)); W(y ~ ., df, tree.prior = cgm(2, .5)); W(y ~ ., df, family = gaussian(sigma = chisq(3, .9)))   # all run
  bart(yb ~ ., df, k = chi(1.5, 2), ...)   # runs when called directly
Why gates missed: test-constructor-vocabulary.R (17) covers wrapper forwarding for interactions etc., not k.
Fix: evaluate k in the prior vocabulary over the written environment as the other prior arguments are, or document that k = chi() must be written in the direct call (or as dbartsPriors$chi).

docs-12  MINOR
Location: man/dbartsControl.Rd \item{proposal.probs}; R/model.R resolveProposalProbs.
Claim: A proposal.probs vector with no recognized names - fully unnamed (0.9-x code often wrote c(0.5, 0.1, 0.4, 0.5)) or misspelled - silently becomes the default mixture. Every other unknown name in the package is refused by name; here the user's mixture is dropped with no message. Same in 0.9-34, so not a regression, but undocumented.
Probe:
  dbartsControl(proposal.probs = c(0.5, 0.1, 0.4, 0.5))@proposal.probs          # 0.6 0 0.4 0 0 0.5
  dbartsControl(proposal.probs = c(birth.death = 0.9, chnage = 0.1))@proposal.probs   # 0.6 0 0.4 0 0 0.5
Why gates missed: tests cover the fill rules over valid names only.
Fix: refuse names outside the six (and an unnamed vector) by name.

docs-13  MINOR
Location: R/spec.R (variance-forest family check), user-facing message.
Claim: The refusal names the requested family token instead of the resolved one, so a 0/1 response under the default reads as if "auto" were the problem.
Probe:
  bart(yb ~ x1 + x2, df, variance = TRUE, n.samples = 20, n.burn = 20, n.chains = 2, n.threads = 1, verbose = FALSE)
  #  Error: a variance forest requires family = "gaussian" or "aft"; family "auto" routes precision through its own latent channel instead
Why gates missed: tests match the first clause of the message.
Fix: name `family` (the resolved token, here "probit") rather than requestedFamily.

docs-14  MINOR
Location: man/sparseFactor.Rd details ("supply dbarts's sigma (sigest in bart)"); man/bartBT.Rd \item{sigest} ("Same concept as sigma in dbarts").
Claim: Both point users at dbarts's sigma, which is now the retired spelling (it warns and is removed in 1.1-0); dbarts's argument is sigest.
Probe: dbarts(y ~ ., df, sigma = 1)   # Warning: 'sigma' is now 'sigest' on 'dbarts'; the value was used. The old name is removed in dbarts 1.1-0.
Why gates missed: prose only.
Fix: say sigest in both places.

Checked and found correct
- NEWS UPGRADING: forwarding of BayesTree-named bart calls (warns once); one-time front-door message on bart(x, y, x.test); fourth positional argument refused; combineChains = TRUE default; bart2 forwards; rngSeed / sigma / node.prior / degreesOfFreedom / power / sigdf / proposal.probs / resid.prior tombstones forward with a message; rngKind refused; $startThreads/$stopThreads no-ops; $run thread count ignored; $sampleNodeParametersFromPrior forwarded; seed/sigest/n.samples NA read as NULL where 0.9-34 defaulted NA, refused where new; probit weights refused except all-1 and 0/1; wrong-length x/y weights an error; negative weights refused; fractional counts refused; unseen factor level an error; NA in a column complete at training refused by default and NA under na.pass; mismatched newdata names an error; $setResponse(y, TRUE) warns about the new second argument.
- NEWS USER-VISIBLE: default proposal.probs and bartBT keeping 0.5/0.1/0.4; chi(1.5, 2) binary default (fit and xbart drop = FALSE label); factors = "categorical"/"indicators" on bart, dbarts, dbartsData, xbart; na.keepPredictors default, bartBT dropping rows, xbart keeping them; auto-family detection (two-level factor/logical/character to probit, unordered 3+ to multinomial on bart only, ordered to ordinal on bart and dbarts, count matrix to multinomial) and verbose = FALSE silencing; excess-threads warning; xbart three-element n.burn refused and loss "log" refused on a continuous response; one-chain combineChains = FALSE keeps the chain margin on extract and predict; extract(type = "trees") always has chain; family-gated sigest warning.
- Seeds: seeded fits identical at n.threads 1 and 2; a single chain reproduces chain 1 of a multi-chain fit; set.seed reproduces and sampling leaves R's stream alone; a supplied seed leaves .Random.seed untouched.
- dbartsControl: fill rules for c(birth_death = 0.7) and c(birth_death = 0.5, change = 0.4); zero-default-only names refused; fractional seed refused; n.samples = 0 accepted by dbarts.
- monotone: direction vocabulary, case sensitivity, unordered factor refused, ordered factor accepted; birth/death-only forcing and refusal of a non-default control mixture; k hyperprior refused; linear leaf and variance forest refused; all-zero spec unconstrained; monotone.prior on print and summary; predictions monotone along a grid, including with missing values in another predictor.
- interactions / blocks: factor names under both factor codings, groups and blocks hold on every saved tree, max.order = 1 additive, total-partition and trees.per.group sum checks.
- dbartsFamilies: base gaussian/binomial mapping (bare binomial is logistic), poisson and non-identity gaussian refused, entry-point table rows checked for xbart student, dbarts hurdle, dbartsSpec hazard.
- na.keepPredictors: fitted/residuals padded to nrow(data) with row names; predict under the default, na.pass, na.omit, na.exclude, na.fail, NULL and zero surviving rows.
- dbartsSampler-class Value: run() shapes at one and two chains, predict with and without keepTrees, getTrees chain column only at > 1 chain, getSumsOfSquaredResiduals is the raw SSR, getFitsWithoutOffset and getLatents shapes.
- Multi-forest formula term (probit and logistic accepted; test refused; type = "forest"), variance forest s.train/s.test for gaussian and aft, survivalProbabilities shape, plotTree chainNum rule, samplePriorPredictive.
- All other Rd examples and the other two vignettes run without error or unexpected warning; codoc/undoc clean.
- Not reported because base R behaves the same: bart/dbarts/xbart called through `function(...) f(...)` with weights or subset fail with "..3 used in an incorrect context", as lm does, and as 0.9-34 did.
