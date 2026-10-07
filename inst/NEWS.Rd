\name{NEWS}
\title{dbarts News}
\encoding{UTF-8}
\section{CHANGES IN VERSION 1.0-0}{
  \subsection{UPGRADING}{
    What can stop code written for 0.9-34 from running, and what to do
    about it. Changes that alter results are under USER-VISIBLE CHANGES.
    \itemize{
      \item R 4.2.0 or newer is required (0.9-34 needed 3.1.0), as is a
            compiler supporting C++20.
      \item \code{bart} is now the formula-first interface that 0.9-34
            called \code{bart2}, with its argument names and defaults (75
            trees, four chains, 500 draws after 500 burn-in). 0.9-34's
            \code{bart}, with BayesTree's argument names and defaults
            (200 trees, one chain, factors expanded to indicator columns,
            incomplete rows dropped), is now \code{\link{bartBT}}; it
            returns an ordinary \code{"bart"} fit. For this release a
            \code{bart} call that uses a name only \code{bartBT} has
            (\code{x.train}, \code{ntree}, \code{ndpost}, and so on) is
            forwarded to \code{bartBT} with a warning. A call that uses no
            such name, such as \code{bart(x, y, x.test)}, runs under the new
            \code{bart}'s defaults and prints a one-time message; rename it
            to \code{bartBT} to keep 0.9-34's defaults. A fourth or later
            positional argument to \code{bart} is an error, since 0.9-34
            read it as \code{sigest}.
      \item \code{bart} defaults to \code{combineChains = TRUE}, where
            \code{bart2} defaulted to \code{FALSE}, so the draws of a
            multi-chain fit come back with the chains merged. Pass
            \code{combineChains = FALSE} for the previous shape.
      \item \code{pdbart} and \code{pd2bart} fit through \code{bart}, with
            its names and defaults (75 trees, four chains merged, so
            \code{fd} has 2000 rows by default; a factor is one predictor).
            Their data arguments are \code{formula} and \code{data};
            positional calls are unchanged. A BayesTree name
            (\code{ntree}, \code{keepevery}, \code{x.train}, ...) is
            translated with a warning until 1.1-0, and an argument neither
            takes is an error where it was ignored. The 0.9-34 model is
            \code{pdbart(bartBT(x, y, keeptrees = TRUE))}, its results still
            subject to the changes below. In a formula fit \code{xind}
            names variables of the data, where it named model-matrix
            columns. Multinomial and ordinal fits are refused, pointing at
            \code{predict}.
      \item Every name and argument kept reachable for one release past its
            removal or rename (a tombstone) warns once per session, names
            its successor, and is removed in dbarts 1.1-0; the list is on
            \code{\link{dbarts-deprecated}}. \code{bart2} forwards to
            \code{bart}. \code{rbart_vi} and its methods
            (\code{predict.rbart}, \code{extract.rbart},
            \code{fitted.rbart}, \code{residuals.rbart}, \code{plot.rbart},
            \code{print.rbart}) are deprecated; see the next item. On
            \code{bart}, \code{power}, \code{base} and \code{split.probs}
            become \code{tree.prior = cgm(power, base, split.probs)}, and
            \code{sigdf} and \code{sigquant} become
            \code{family = gaussian(sigma = chisq(df, quant))}; on
            \code{bart}, \code{dbarts} and \code{xbart},
            \code{resid.prior} becomes \code{family = gaussian(sigma = )};
            on \code{bart} and \code{dbarts}, \code{proposal.probs} becomes
            \code{control = dbartsControl(proposal.probs = )}.
            \code{rngSeed} (\code{bart}, \code{dbartsControl}) is now
            \code{seed}, \code{sigma} (\code{dbarts}, \code{xbart}) is now
            \code{sigest}, \code{node.prior} (\code{dbarts}) is now
            \code{leaf.prior}, and \code{chi()}'s \code{degreesOfFreedom} is
            now \code{df}. On the sampler,
            \code{$sampleNodeParametersFromPrior} is now
            \code{$sampleLeafParametersFromPrior}, \code{$startThreads} and
            \code{$stopThreads} do nothing, and a thread count passed to
            \code{$run} is ignored (use \code{$setControl}). \code{NULL}
            rather than \code{NA} now means \dQuote{not given}; where
            0.9-34's default was \code{NA}, an explicit \code{NA} is still
            read as \code{NULL} with a warning.
      \item An \code{n.cuts} with more values than there are predictor
            columns is an error; 0.9-34 dropped the extra values silently.
            A shorter one is recycled, as before.
      \item A package that forwards arguments to dbarts by filtering them
            against \code{formals(dbartsControl)} will now silently drop
            \code{rngSeed} and fit unseeded; accept \code{seed} instead.
      \item \code{rbart_vi} still runs for this release, with the fixes
            listed under BUG FIXES, but multilevel models belong to
            \code{stan4bart::stan4bart}, which fits the same random
            intercepts and slopes, and more, through an lmer-style formula.
            Its prior on the group spread differs, so results move. A
            binary \code{rbart_vi} fit now refuses weights other than 0 and
            1, and takes 1.0-0's defaults for \code{k} and the tree moves.
      \item Removed without replacement: \code{dbartsControl}'s
            \code{rngKind} and \code{rngNormalKind} (every chain now runs
            its own generator), and the \code{configure} options
            \code{--enable-match-bayes-tree},
            \code{--enable-thread-safe-unload} and \code{--with-xint-size},
            which stop the build naming the removal.
      \item Sampler states saved by 0.9-34, including those inside saved
            fits, cannot be restored, so such a fit cannot \code{predict};
            refit it. A \code{dbartsData} object
            saved by 0.9-34 keeps its old factor coding and fits
            differently from one built fresh; rebuild it with
            \code{\link{dbartsData}}.
      \item Weights. A \code{probit} fit refuses case weights: 0.9-34 fit a
            weighted probit, which was incorrect. Weights that are all 1
            are ignored, and weights of 0 and 1 mark which rows enter the
            likelihood; count weights belong on
            \code{family = "logistic"}. On the \code{x}/\code{y} interface
            a \code{weights} vector of the wrong length is an error rather
            than recycled, and negative or non-finite weights are refused
            everywhere.
      \item New data given to \code{predict}. Test column names that do not
            cover the training columns are an error (0.9-34 warned and
            matched by position). A factor level that training never saw is
            an error (0.9-34 silently predicted it as another level); under
            \code{factors = "indicators"} and in \code{bartBT} that includes
            a level declared but with no training rows. There, a test factor
            is now matched to the training levels by label, as in
            \code{lm}, so one declaring only some of the levels predicts as
            the full-level factor would, where 0.9-34 failed. A
            missing value in a column that had none in training is an error
            by default, naming the column (0.9-34 sent it down an arbitrary
            branch); the new \code{na.action} argument of every
            \code{predict} method chooses otherwise, \code{na.pass} returning
            \code{NA} for such rows. 0.9-34 also dropped a data frame's
            incomplete rows, returning fewer rows than given; see
            \code{\link{na.keepPredictors}}.
      \item An unnamed \code{test} or \code{newdata} matrix given to a fit
            whose predictors are named warns on every call, saying which
            columns it was matched to by position; 0.9-34 matched silently.
            Give the matrix column names to match by name.
      \item Argument checking is stricter. A fractional value for a count
            (\code{n.trees = 100.5}, say) is an error rather than truncated.
            \code{predict}, \code{extract}, \code{fitted} and
            \code{residuals} warn about, or refuse, arguments they do not
            use instead of ignoring them. On \code{predict.bart},
            \code{type} now comes before \code{offset} and \code{weights},
            and \code{n.threads} is last; pass them by name.
            \code{dbartsControl}'s \code{proposal.probs} must name every
            entry from its moves; an unnamed vector or an unknown name is
            an error listing them, where 0.9-34 dropped such entries and ran
            the default mixture. New arguments make three shortened names
            ambiguous, so they are errors: \code{t =} and \code{u =} on
            \code{bart} (0.9-34's \code{bart2} read them as \code{test} and
            \code{useQuantiles}) and \code{n.th =} on \code{xbart} (read as
            \code{n.threads}); write the name in full.
      \item The sampler's \code{$setResponse} gains \code{updateScale} as its
            second argument, pushing \code{updateState} to third, so pass
            \code{updateState} by name. It defaults to \code{FALSE}: the
            response's scaling stays as set at creation when the response
            is swapped, as a conditional update inside a larger sampler
            needs. 0.9-34 re-derived it from each new continuous response;
            \code{updateScale = TRUE} restores that.
            \code{$getSumsOfSquaredResiduals} returns the sum of squares
            itself rather than that sum divided by the response range.
      \item \code{names(fit)} no longer lists \code{yhat.test}, or, on a
            binary \code{bart} fit without an offset, \code{binaryOffset},
            when the fit has none; test with \code{is.null} rather than \code{names}.
      \item A formula calling a function from no package, such as one
            defined in the session, is refused, naming it; compute the
            column in \code{data} instead. A fit stores its formula's terms
            and \code{offset} expression without their environment, with
            every name that is not a data column fixed at its value when
            fitting, so a fit predicts from \code{newdata} alone, the same
            after a reload, and carries no frame of the function it was made
            in.
      \item A fit made with \code{keepCall = FALSE} stores \code{call =
            NULL}, as base R fits without a call do, so \code{update}
            refuses it as it refuses such an \code{lm} fit; 0.9-34 stored
            the placeholder call \code{NULL()}.
      \item The sampler's \code{$setControl} refuses a control whose
            \code{n.trees} differs from the sampler's; 0.9-34 resized the
            sampler in place. Create a new sampler with the count wanted.
    }
  }
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item Results differ from 0.9-34, including seeded ones. The new
            engine draws a different random stream, and the corrections
            under BUG FIXES change the posterior itself. Every chain now
            runs its own generator: seeded results do not depend on the
            thread count, \code{set.seed} alone makes any run reproducible,
            and sampling no longer advances R's random stream, so R-level
            draws made after a fit also change.
      \item The default tree-move mixture drops the swap move:
            \code{proposal.probs} defaults to \code{c(birth_death = 0.6,
            swap = 0, change = 0.4, perturb = 0, rule_gibbs = 0, birth =
            0.5)} rather than \code{c(birth_death = 0.5, swap = 0.1, change
            = 0.4, birth = 0.5)}. At usual forest sizes swap almost never
            changes anything, and dropping it costs nothing measurable in
            accuracy, coverage or mixing. A single-tree fit
            (\code{n.trees = 1}) should still set it positive.
            \code{bartBT} keeps 0.9-34's mixture.
      \item The default prior on \code{k} for a binary response, on
            \code{bart}, \code{dbarts} and \code{xbart}, is the proper
            \code{chi(1.5, 2)}, where 0.9-34 used the improper
            \code{chi(1.25, Inf)}. It centers \code{k} near the usual fixed
            value of 2 while still adapting to the data. \code{chi()}'s own
            \code{scale} now defaults to 2 rather than \code{Inf}. Pass
            \code{k = chi(1.5, Inf)} for the prior 0.9-34 actually sampled
            (see the \code{chi} fix under BUG FIXES).
      \item Unordered factor predictors are split on subsets of their
            levels, and an ordered factor is one ordinal predictor split
            between its levels, where 0.9-34 expanded both into indicator
            columns. The design matrix, the \code{varcount} column names and
            the draws change for any fit with a factor predictor. Pass
            \code{factors = "indicators"} (\code{bart}, \code{dbarts},
            \code{dbartsData}, \code{xbart}) for the previous expansion;
            \code{bartBT} expands as before.
      \item Rows with missing predictors are kept and their missingness
            modelled, where 0.9-34 dropped them. The new \code{na.action}
            of \code{bart}, \code{dbarts} and \code{dbartsData} defaults to
            \code{\link{na.keepPredictors}}, which drops only rows missing
            the response; \code{na.action = na.omit} restores the old rule.
            \code{xbart} always keeps such rows, and \code{bartBT} drops
            them as before.
      \item \code{family = "auto"} recognizes more responses. A two-level
            factor, logical or character response fits as probit. A factor
            with three or more levels fits as \code{"multinomial"} on
            \code{bart} if unordered and as \code{"ordinal"} on \code{bart}
            and \code{dbarts} if ordered, and is an error naming the family
            elsewhere; a matrix of three or more count columns fits as
            multinomial. 0.9-34 fit a factor's level codes as a continuous
            response, or failed on the formula interface. A one-line
            message reports the family chosen; \code{bart} prints it by
            default, and \code{verbose = FALSE} silences it.
      \item \code{n.threads} keeps its meaning as a total thread budget, but
            each chain samples on one thread, and \code{dbartsControl}'s
            default is now \code{min(guessNumCores(), n.chains)}.
            \code{bart} and \code{dbarts} warn when \code{n.threads} exceeds
            \code{n.chains}; the extra threads serve test fits and
            \code{predict} only.
      \item \code{xbart}. With \code{k} not given it scores the default a
            \code{bart} fit would use: fixed 2 for a continuous response and
            \code{chi(1.5, 2)} for a binary one, where 0.9-34 used
            \code{chi(1.25, Inf)}; with \code{drop = FALSE} the result
            then has a labelled \code{k} dimension. \code{n.burn} takes two
            values; a third is an error. A chain is no longer carried from
            one fold or replication to the next (0.9-34 could score a
            setting on rows it had effectively been trained on), the
            default \code{sigest} is estimated from each fold's training
            rows rather than from all rows, held-out responses included,
            and every fold uses the full data's cut points, so results
            change. Seeded
            runs are reproducible at any \code{n.threads} and leave
            \code{.Random.seed} as they found it, and warnings from the
            workers are now reported.
      \item A predictor holding \code{Inf} or \code{-Inf} has its cut points
            placed over its finite values, the infinite ones falling past the
            end cut points; 0.9-34 stretched the cut points to infinity,
            leaving the column unused, or under \code{useQuantiles = TRUE}
            placed an infinite cut point.
            When the starting \code{sigest} cannot be estimated because of
            such a value, the error names the column.
      \item \code{pdbart} and \code{pd2bart} evaluate a factor predictor
            at every level by default, take and report its levels by name
            in \code{levs}, and plot one point per level.
      \item \code{pdbart}'s and \code{pd2bart}'s \code{fd} includes the
            fit's offset, as BayesTree's included \code{binaryOffset}; rows
            a fit gives a 0 weight are left out of the average. The result
            carries \code{n.chains} and its \code{type}, and
            \code{combineChains = FALSE} gives \code{fd} a chain margin.
      \item Output shapes. Per-observation results carry the data's row
            names. Under \code{combineChains = FALSE}, \code{extract} and
            \code{predict} keep the chain dimension on a one-chain fit, so a scalar such as \code{sigma}
            is a 1 x samples matrix rather than a vector.
            \code{extract(type = "trees")} always has a \code{chain}
            column.
      \item \code{predict} applies the fit's \code{offset} argument at
            \code{newdata}, as \code{predict.lm} does: the expression it was
            given as, such as \code{offset = log(exposure)}, is evaluated on
            \code{newdata}, and an \code{offset} given to \code{predict}
            adds to it. One given as a plain vector for the training rows
            applies only to a \code{newdata} of as many rows, as in
            \code{lm}, and is otherwise refused unless \code{predict} is
            given an \code{offset} for the new rows.
            0.9-34 predicted without the fit's offset. The default
            \code{offset.test} is formed the same way on \code{test}.
      \item An \code{offset()} term in a formula is used, as in \code{lm}:
            it is added to any \code{offset} argument at fit time and
            evaluated on \code{test} and on \code{predict}'s
            \code{newdata}, which must then be data frames carrying its
            variables. 0.9-34 ignored it unless an \code{offset} vector was
            also given. An explicit
            \code{offset.test} is one value or one per test row, where
            0.9-34 recycled or cut it, and a bare name given as
            \code{offset.test} is looked up in \code{test} before
            \code{data}.
      \item The warnings a caller has a reason to catch apart from the rest
            carry a class inheriting from \code{dbartsWarning}, so they can
            be caught or muted without matching their text: deprecations
            (\code{dbartsDeprecatedWarning}, also a base R
            \code{deprecatedWarning}), the slow-count warning of a monotone
            fit (\code{dbartsSlowCountWarning}), and warnings that a supplied
            input was replaced by another (\code{dbartsFallbackWarning},
            \code{dbartsSigmaFallbackWarning},
            \code{dbartsSparseSigmaFallbackWarning} and
            \code{dbartsGPFallbackWarning}). Every other warning is an
            ordinary one. \code{bart} now warns when a supplied
            \code{sigest} has no effect because the family fixes the
            residual scale; 0.9-34 ignored it silently.
      \item Draws a sampler has saved (\code{keepTrees}) keep their
            values when the response scale is re-derived: after
            \code{setResponse} or \code{setOffset} with
            \code{updateScale = TRUE}, or \code{setData}, \code{predict}
            returns for them what it returned before the call. 0.9-34
            replayed saved trees on the scale in force, so they came back
            in the new response's units.
    }
  }
  \subsection{NEW FEATURES}{
    \itemize{
      \item \code{pdbart} and \code{pd2bart} take \code{type}
            (\code{"auto"}: the link scale, 0.9-34's, except the mean
            response on hurdle fits and survival at the median survival
            time on accelerated failure time and hazard fits, with
            \code{"event"}, \code{"cumhaz"} and \code{times}),
            \code{newdata}, \code{n.average.rows}, \code{average.weights}
            and \code{n.max.predictions}; they serve negative binomial,
            hurdle, accelerated failure time and hazard fits, and their
            plots label the scale and take \code{plot.type}.
      \item Sampling is faster. Against 0.9-34 on a four-core x86-64
            machine, fitting through the BayesTree-style interface with 200
            trees and 1000 draws after 500, with defaults that changed
            pinned to 0.9-34's: 1.7 times faster at 1,000 observations, 2.6
            times at 10,000, 1.7 times for four chains on four threads, and
            1.6 times for a probit fit. The gap widens with the number of
            observations. \code{storage = "single"} (a
            \code{\link{dbartsControl}} and \code{bart} argument) keeps the running residual in single
            precision for a further 10 to 30 percent at a million
            observations; it is opt-in, for continuous responses with
            constant leaves.
      \item Large fits use less memory: packaging a fit no longer holds an
            extra copy of the training and test draws, which lowers peak
            memory by 1.6 GB on a 100,000-observation, four-chain, 500-draw
            fit. \code{\link{dbarts-package}} has a \sQuote{Memory} section
            giving a fit's footprint and the options that control it. A
            per-draw compiled \code{callback} with \code{keepFits = FALSE}
            (\code{bart}, \code{dbarts}, \code{$run}) lets a caller reduce
            draws as they are made without storing them.
      \item Response families. \code{bart}, \code{dbarts} and \code{xbart}
            take a \code{family}, as a name or as a family object with its
            own settings (not every family on every function);
            \code{\link{dbartsFamilies}} lists them, and base
            R's \code{gaussian()} and \code{binomial()} are accepted as
            \code{glm} takes them. New are \code{"logistic"} (Polya-Gamma
            augmentation, accepting integer count weights),
            \code{student()} (robust Student-t errors),
            \code{"multinomial"}, \code{"ordinal"} (cumulative probit),
            \code{"nbinom"} (negative-binomial counts on a log link, with
            an integer shape fixed or sampled on a grid), \code{"aft"}
            (log-normal accelerated failure time for a right-censored
            \code{Surv} response), \code{"hazard"} and
            \code{"hazard.logistic"} (discrete-time survival), and
            \code{"hurdle.lognormal"} (a non-negative response with exact
            zeros, \code{bart} only). \code{\link{survivalProbabilities}}
            gives survival-curve draws from a survival fit. Each new fit
            class has \code{fitted}, \code{predict}, \code{extract},
            \code{residuals}, \code{summary} and \code{plot} methods, and
            \code{family()} returns the family of any fit. See
            \code{\link{bart}}.
      \item Leaf models. \code{leaf.prior = linear(columns)} fits a linear
            function of chosen predictors in each leaf, and
            \code{leaf.prior = gp(columns)} a Gaussian-process function.
            \code{normal(sd = )} states the leaf prior's spread directly on
            the response (or latent) scale instead of through \code{k}, as a
            number or under the new \code{invchi()} hyperprior. Available on
            \code{bart}, \code{dbarts} and \code{xbart}; see
            \code{\link{dbartsPriors}}.
      \item Tree priors. \code{tree.prior = dart()} runs DART (Linero 2018),
            a sparsity-inducing prior on which predictors the trees split
            on, and \code{cgm(split.probs = )} sets fixed, unequal
            split probabilities. \code{bart} gains \code{tree.prior} and
            \code{leaf.prior}, taking the same prior objects as
            \code{dbarts}, which are listed in \code{\link{dbartsPriors}}.
      \item Multi-forest models. \code{dbarts(forests = )} takes a list of
            \code{\link{forest}} specifications, each with its own trees,
            priors and predictors, and fits the mean as their sum, a forest
            given a \code{basis} entering multiplied by it. The Bayesian
            causal forest is the two-forest case with a treatment
            indicator as the basis. In a \code{bart} or \code{dbarts}
            formula, each forest is a term:
            \code{y ~ forest(x1 + x2) + forest(x1 + x2, basis = z)}
            declares the same thing. Gaussian, probit and logistic responses;
            \code{predict} and \code{extract} report each forest separately
            (\code{type = "forest"}) or combined.
      \item Heteroscedastic BART (Pratola et al. 2020): \code{variance = ~ x1
            + x2}, \code{variance = TRUE} or a
            \code{\link{varianceForest}} in \code{bart} and \code{dbarts}
            adds a forest for the residual variance, reported as
            \code{s.train} and \code{s.test}, and returned by
            \code{extract} and \code{predict} as \code{type = "sigma"}
            in place of a scalar sigma. For gaussian and \code{"aft"}
            responses.
      \item Constraints. \code{monotone = c(x1 = "increasing", x2 =
            "decreasing")}, or \code{1}, \code{-1} and \code{0} by name
            or position, makes the fit monotone in those predictors
            (Chipman et al. 2022). \code{\link{monotone}(directions,
            prior = )} also chooses the prior: \code{"joint"}, the
            default, conditions the trees and their leaf values on
            monotonicity together; \code{"leaf"} keeps the ordinary tree
            prior and restricts only the leaf values, at a cost in run
            time with very few trees.
            \code{interactions = \link{interactions}(max.order = , forbid =
            , groups = )} limits which predictors may interact, and
            \code{blocks = \link{blocks}(...)} makes the fit an exact sum
            of functions of separate predictor groups. All three on
            \code{bart} and \code{dbarts}; \code{interactions} and
            \code{blocks} also hold on each category forest of a
            multinomial fit.
      \item Predictors. Factors split on level subsets, with up to 65535
            levels; missing predictor values are modelled in place, each
            split learning which way they go (Twala et al. 2008). Sparse
            \code{Matrix::dgCMatrix} predictors are used without being
            expanded, and a data frame may mix
            ordinary columns with sparse ones, including the new
            \code{\link{sparseFactor}} column. Any other
            \code{Matrix} class, dense or sparse, and a \code{sparseVector}
            are converted, and integer and logical
            matrices are taken as numbers, in fitting and in
            \code{predict}; 0.9-34 refused an integer matrix.
      \item Starting points. \code{bart}'s \code{warm.start} starts the
            chains from a previous fit's trees, and \code{n.grow.sweeps}
            builds the starting trees by XBART-style grow-from-root sweeps
            (He et al. 2019); the sampler offers both as
            \code{$installTrees} and \code{$growFromRoot}. Neither changes
            the posterior, only where the chains start.
      \item Tree moves. \code{proposal.probs} gains two optional moves,
            \code{"perturb"} (shift one split point) and
            \code{"rule_gibbs"} (redraw a split rule from its conditional
            distribution), both off by default, and accepts all-zero
            structural probabilities, holding the trees' shapes fixed while
            everything else is sampled. \code{dbartsControl(treeShift = )}
            adds a step that improves mixing of the leaf values when the
            shapes are fixed. \code{bart} gains a
            \code{control} argument for these and other
            \code{\link{dbartsControl}} settings.
      \item Prediction and summaries. \code{fitted} and \code{predict} take
            \code{ci.level} to return credible or prediction intervals.
            \code{predict}'s \code{n.threads} splits the work across draws,
            so single-chain fits parallelize too, with identical results.
            \code{\link{summary.bart}} tabulates the parameters a fit
            sampled, \code{sigma} and the leaf scale (\code{k}, or
            \code{leaf.prior.sd} on a fit whose leaf prior was named by
            \code{sd}), with split R-hat and effective sample sizes, and
            names every one it held fixed, a Student-t df included, on a line
            under the table.
            \code{extract(type = "loglik")} returns the pointwise
            log-likelihood, for use with WAIC or PSIS-LOO (for example
            \pkg{loo}), and \code{extract} also reads \code{"sigma"},
            \code{"k"}, \code{"leaf.prior.sd"} and \code{"varcount"}:
            draws for a parameter the fit sampled and one number for one it
            held fixed, \code{"sigma"} being 1 on a family with no residual
            scale, and \code{"k"} and \code{"leaf.prior.sd"} one number
            per forest on a fit with several. A fit stores the \code{k} its
            sampler recorded, the \code{leaf.prior} it ran under and the
            scalars it held \code{fixed}. \code{\link{plotTree}} plots a
            tree directly from a fit. \code{\link{samplePriorPredictive}}
            draws from a model's prior before fitting it. Runs can be
            interrupted with Ctrl-C.
      \item \code{xbart} runs over the same sampler as \code{bart}, so it
            takes \code{family}, \code{factors}, \code{tree.prior},
            \code{leaf.prior} and sparse predictors; it accepts only
            \code{storage = "double"}. Its
            \code{k} may mix fixed values and hyperpriors
            (\code{k = list(1, 2, chi())}), \code{sd} crossvalidates the
            leaf spread, and \code{parallel} and \code{cl} choose forked,
            socket or caller-supplied workers.
      \item Using dbarts inside a larger sampler. New sampler methods include
            \code{$setActiveRows} (drop rows from the likelihood between
            draws), \code{$getFitsWithoutOffset}, \code{$getLeafPrior},
            \code{$setLeafPrior}, \code{$getK}, \code{$getVariance}, and per-forest
            readers and setters; an \code{"aft"} sampler's censoring can be
            changed through \code{$setResponse(y, status = )}.
            \code{\link{dbartsDrawLatents}} and \code{dbartsWorkingResponse}
            expose a family's augmentation draws to R. \code{\link{dbartsSpec}}
            resolves a model specification without building a sampler, for
            packages supplying their own design.
            \code{\link{dbartsValidateComposition}} checks by
            simulation-based calibration that a composed sampler targets the
            posterior it claims. The vignette
            \code{vignette("dbarts-as-a-component")} and the help page
            \code{\link{dbarts-embedding}} show worked compositions.
    }
  }
  \subsection{C API}{
    \itemize{
      \item The C++ headers (\code{dbarts/R_C_interface.hpp} and those
            beside it) are removed with the engine they exposed. A package
            compiling against dbarts uses the flat C header
            \code{dbarts/dbarts.h}, whose entry points are reached through
            \code{R_GetCCallable} (or the stubs enabled by defining
            \code{DBARTS_USE_STUBS}) and versioned by
            \code{DBARTS_C_API_MAJOR} and \code{DBARTS_C_API_MINOR}, with
            \code{dbarts_apiHash()} as an opt-in exact-ABI check. The header
            declares no R type and no longer defines \code{USE_FC_LEN_T} or
            includes \code{Rversion.h}; a consumer that needs them includes
            them itself.
      \item The header creates no sampler: build it in R
            (\code{\link{dbarts}} or \code{\link{dbartsSpec}}) and read its
            handle out of the object's external pointer with
            \code{R_ExternalPtrAddr}. The C entries run the sampler into a
            caller-owned \code{dbarts_results} (initialize it with
            \code{DBARTS_RESULTS_INIT}), swap its response, offset and
            residual scale (copying them), read its latents, predict, store
            and print trees, and register a per-draw callback; saving state,
            reading trees and changing predictors or weights go through the
            R methods. See \code{\link{dbarts-embedding}}.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item The tree change move now satisfies detailed balance. It omitted
            a proposal-density ratio, which biased the trees toward
            predictors with fewer cut points whenever predictors differed in
            their number of cut points, so variable selection in mixed
            continuous and categorical designs, and under DART, is now
            correct.
      \item With \code{useQuantiles = TRUE} (\code{usequants} on
            \code{bartBT}), a predictor with more distinct values than its
            cut points can separate has them spread evenly over the whole
            column. They had been taken at a fixed stride from the bottom,
            which could stop well short of the top, where no tree could
            then split: at the default 100 cut points, a third of the
            distinct values lay above the highest one at 150 distinct
            values, half at 199 and a tenth at 999. Results change for
            nearly every such fit, \code{bartBT}'s included, and likewise
            for a
            \code{$setPredictor(updateCutPoints = TRUE)} onto more distinct
            values than the column's cut points can separate, which had
            kept only the lowest ones, with a warning. A predictor with no
            more than \code{n.cuts + 1} distinct values, and the default
            \code{useQuantiles = FALSE}, are unaffected.
      \item The \code{chi(df)} hyperprior on \code{k} samples the degrees of
            freedom asked for; it had sampled \code{chi(2 * df - 1)}.
      \item Observations of weight zero no longer count in the residual
            variance's degrees of freedom, which had biased \code{sigma}
            low. Such an observation is still in the design, as before: it
            occupies a leaf and receives a fitted value.
      \item The starting forest, and trees drawn from the prior by
            \code{$sampleTreesFromPrior}, now follow the sampler's tree
            prior, which excludes empty leaves; they had been drawn without
            that restriction and pruned afterward.
      \item A fixed residual scale, \code{fixed(value)}, given when a
            sampler is created is now held; it had been ignored, and one
            installed through \code{$setModel} had been read as a standard
            deviation.
      \item Prior arguments no longer see internal variables named
            \code{control} and \code{data}, which had shadowed the caller's
            own (as in \code{tree.prior = cgm(2, 0.95, data)}).
      \item An \code{offset.test} expression is evaluated in \code{test},
            then in \code{data}, each with the formula's environment, as
            \code{model.frame} evaluates its arguments, and failing that in
            the caller's frame; without a formula, in the caller's frame. A
            formula fit's \code{test} and \code{weights} are evaluated in
            \code{data} with the formula's environment, and otherwise taken
            as given. They had seen internal variables such as \code{x}, and
            then the global environment ahead of the caller's own objects, so
            \code{offset.test = d$o} inside a function could read a global
            \code{d}.
      \item \code{xbart}'s weighted RMSE loss is on the scale of the
            response, and \code{loss = "log"} on a continuous response is an
            error rather than silently computing RMSE.
      \item \code{bartBT} (0.9-34's \code{bart}) keeps its \code{weights}
            on the fit, so \code{extract(type = "ppd")} on a weighted fit is
            no longer computed unweighted.
      \item \code{update} works on a fit made through a wrapper that passes
            its \code{\dots} on to \code{bart}, \code{bartBT}, \code{dbarts}
            or \code{xbart}; it failed with \dQuote{..1 used in an incorrect
            context}.
      \item A data-dependent formula term, such as \code{poly()},
            \code{ns()} or \code{scale()}, is rebuilt on \code{test} and on
            \code{predict}'s \code{newdata} from the training data, as
            \code{predict.lm} rebuilds it; 0.9-34 recomputed it from the new
            rows, so a prediction depended on the other rows given with it.
      \item A formula with no response is refused by \code{bart} and
            \code{xbart}; 0.9-34 failed, or silently fit an all-zero
            response when \code{sigest} was given.
      \item Several inputs that failed with an unrelated low-level error now
            get a clear one: a non-finite response, zero rows, a multi-column
            response to \code{dbartsData}, a sample count that thins to
            zero, a \code{:} or \code{*} term in a formula, heavy thinning
            that drove \code{printEvery} to zero, and \code{predict} on a
            multi-chain fit, or \code{predict(type = "ppd")} on any fit,
            kept without \code{keepTrees}.
      \item The starting estimate of \code{sigma} no longer depends on the
            machine's linear algebra library: an exactly fit response could
            give an estimate of 0, which the sampler then refused.
      \item \code{plot.bart}, \code{plot.pd2bart} and tree plots restore the
            caller's graphical parameters.
      \item Sampler methods. \code{$setResponse} and \code{$setData} check
            that a probit response is 0 or 1. A whole-matrix
            \code{$setPredictor} keeps column names, so later by-name column
            updates work. A sampler run more than once with
            \code{keepTrees} replays its saved trees in the order they were
            drawn; they had been rotated against the other results.
      \item \code{pd2bart} on a model with exactly two predictors returns
            one column per grid point, where it averaged over the whole
            grid and \code{plot} then failed.
      \item \code{plot} on a \code{pdbart} result accepts \code{type},
            \code{xlab} and \code{ylab}, which collided with its own.
      \item Arguments given beside a \code{dbartsData} object, which it
            already fixes, such as \code{subset} and \code{weights}, are
            reported in a warning rather than ignored silently.
      \item \code{$setCutPoints} refuses a \code{NaN} cut point.
      \item \code{rbart_vi}'s \code{prior} may name the caller's own
            function, or be passed on by a wrapper, and an unknown name is
            refused by name; only the built-in names worked.
            \code{predict} takes a character or numeric \code{group.by},
            as the fit does; random effects came back empty and unseen
            groups failed.
      \item An \code{n.cuts} above 65533, the most a predictor's grid can
            hold, is an error; 0.9-34 accepted it and could hang.
      \item \code{predict} on a fit kept with \code{keepSampler = TRUE}
            and \code{keepTrees = FALSE} is an error when the fit has more
            than one chain, asking for \code{keepTrees = TRUE}; 0.9-34
            returned the first chain's values and arbitrary numbers for
            the others. A one-chain fit still predicts from its current
            trees.
      \item \code{rbart_vi}: the intercepts are drawn with precision weights
            and start from a sweep of the trees, where a start from prior
            trees could inflate the posterior of \code{tau} many times over
            on a skewed response; with an offset and a test set of the same
            size, a group's intercept is no longer counted twice in the test
            fits; a fit saved and reloaded predicts as before;
            \code{predict} with several chains handles groups the fit did
            not see; a one-chain \code{extract(type = "trees")} works;
            a \code{group.by} naming a variable that is not a column of
            \code{data} is found in the calling environment, as documented,
            where 0.9-34 grouped by the first column of \code{data} without
            a warning; a seeded call leaves
            the random stream as it found it; and \code{set.seed}
            reproduces a multithreaded run at the same \code{n.threads}.
    }
  }
}
\section{CHANGES IN VERSION 0.9-34}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item `setPredictor` accepts `forceUpdate = "partial"`, installing a column
            one observation at a time and rolling back only those that would empty
            a node.
      \item `updatePredictorPerObservationJointly` applies a `"partial"`-style update across several
            samplers at once, installing an observation only when it keeps every
            sampler's nodes non-empty.
      \item The `getTrees` method of a sampler accepts `current = TRUE` to return
            the live working trees (the current sampler position) even when
            `keepTrees` is `TRUE`, rather than the saved samples.
      \item `getTrees` and `extract(., "trees")` accept a `newdata` matrix, routing
            those observations through each tree so the `n` column counts them
            instead of the training data.
    }
  }
  \subsection{C API}{
    \itemize{
      \item Added `dbarts_getTreesCountingData`, which routes a supplied matrix
            through the trees for the `newdata` feature above.
      \item `dbarts_getTrees` gained a trailing `useLiveTrees` argument (for the
            `current = TRUE` feature); callers linking to dbarts must update.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item `copy` now works for samplers with a stored state.
      \item Fixed reading a saved tree back from binary IO, which recursed with
            the writing routine instead of the reading one.
      \item Fixed installation under R-devel, where \code{RNGkind} gained an
            additional \code{binom.kind} argument and the code that parsed
            its body via \code{deparse} produced invalid R once the
            deparsed header spanned two lines. Fix thanks to Brian King
            (PR #81).
      \item Replaced a \code{std::bit_cast} cast used when registering the
            \code{.Call} interface with a plain \code{DL_FUNC} cast; Apple
            clang 14.0.3 on Intel macOS accepts \code{-std=gnu++20} but its
            libc++ lacks \code{std::bit_cast}, so compilation failed. Fix
            thanks to Brian King (PR #81).
      \item Fixed multi-chain and multi-threaded runs under R-devel, which
            warned that \code{.Random.seed} was not a valid integer and fell
            back to seeding from the clock, so that \code{set.seed} no longer
            controlled their results. R's new \code{binom.kind} field pushed
            the encoded value past the bound used to validate it; the uniform
            and standard normal kinds that dbarts uses are now validated
            individually and any further fields R encodes are ignored.
      \item The fallback used to validate \code{rngKind} when the arguments of
            \code{RNGkind} cannot be read no longer fails in a session that has
            yet to use the random number generator, and now restores the
            generator's state and all of its kinds afterwards.
      \item Slot and attribute temporaries in the \code{.Call} interface are now
            protected across the calls that can allocate, as reported by
            \code{rchk}'s static analysis of the protection stack.
    }
  }
}
\section{CHANGES IN VERSION 0.9-31}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item \code{predict} now accepts a \code{n.threads} argument and will use native
            threads to parallelize the saved-tree replay, one task per chain.
            (The chain axis is the 0.9-x one; 1.0-0 partitions per chain and
            posterior draw instead, and the argument's spelling is unchanged.)
    }
  }
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item When running with multiple threads, threads no longer persist after
            \code{run}ing or \code{predict}ing and will be re-created on the next call.
            Prior to this, threads were initialized at object creation and
            would run in the background for the life-cycle of the sampler
            object, waiting for work. To prevent frequent stopping/starting of
            threads, new sampler-object functions \code{startThreads()} and
            \code{stopThreads()} have been added that can bookend sampling or
            prediction loops. 
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item \code{rbart_vi} no longer fails when there are multiple out-of-sample
            grouping factors in test data. Bug report thanks to github user
            drcthal.
      \item Should now correctly install on non-Mac OS arm64 systems. Bug
            report thanks to Bob Dröge.
      \item \code{rbart_vi} now works when \code{k} is a variable. Bug report thanks
            to Bruno Tancredi.
    }
  }
}
\section{CHANGES IN VERSION 0.9-29}{
  \subsection{BUG FIXES}{
    \itemize{
      \item \code{dbarts} now works to initialize a sampler for binary outcomes when
            \code{resid.prior} is left at its default. Bug report thanks to github
            user LennMass.
      \item Installs correctly on Alpine Linux. Bug report thanks to Sebastian
            Meyer.
    }
  }
}
\section{CHANGES IN VERSION 0.9-27}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Fixed an error in sampler when setting test offsets when there
            was only a single test observation.
    }
  }
}
\section{CHANGES IN VERSION 0.9-24}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item Issues a warning instead of failing when weights are present in
            training data but not present in test. Suggestion thanks to
            github user Pentaonia (Loubert).
    }
  }
}
\section{CHANGES IN VERSION 0.9-23}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item No longer depends on gfortran.
      \item Uses SIMD instructions on M1 Macs.
      \item Added experimental \code{callback} functionality to
            \code{rbart_vi}.
    }
  }
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item Custom loss functiosn for \code{xbart} now require an additional
            weights argument.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Fixed a multithreaded issue leading to inconsistent results
            with \code{xbart}.
      \item \code{rbart_vi} should now correctly use default arguments.
      \item \code{rbart_vi} now works with \code{keepTrainingFits} as false.
      \item Weighted binary responses sample latent variables from the
            correct distribution.
      \item Extracting the values from the posterior predictive distribution
            for models with weights now incorporates them into the variance.
      \item Weighted values are considered in loss functions for
            crossvalidation.
    }
  }
}
\section{CHANGES IN VERSION 0.9-21}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item \code{extract} now accepts as a type \code{"trees"}, which
            allows for easier inspection of models fit with
            \code{"keepTrees"} as \code{TRUE}.
      \item \code{print} generics now exist for \code{bart} and \code{rbart}
            fits; implementation thanks to Emil Hvitfeldt.
      \item \code{xbart} now accepts a \code{seed} argument to enhance
            reproducibility.
      \item \code{bart}/\code{bart2} (and \code{dbarts} through its
            \code{tree.prior} argument) accept \code{splitprobs}/
            \code{split.probs} which controls the prior probability that any
            variable is used when splitting observations.
    }
  }
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item \code{fitted} for \code{rbart_vi} models now uses a C++
            implementation for the expected value that uses less memory
            and is faster.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item \code{xbart} for binary outcomes with log loss no longer returns
            NaN when some subset of the response is perfectly predicted by the
            covariates. Bug report thanks to Marcela Veselkova.
    }
  }
}
\section{CHANGES IN VERSION 0.9-20}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item \code{dbarts} now exposes access to the underlying proposal
             rules and their probabilities through its \code{proposal.probs}
             argument. \code{bart2} response to the same argument, while
             \code{bart} uses \code{proposalprobs}.
      \item \code{bart}, \code{bart2}, and \code{rbart_vi} accept a
             \code{seed} argument that will yield reproducible results, even
             when running with multiple threads and multiple chains.
    }
  }
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
       \item The interface registered under \code{R_RegisterCCallable} has
             changed to reflect proper fixed hyperpriors for \code{k}.
       \item Samples of the end-node sensitivity parameter, \code{k},
             are returned by \code{rbart_vi} when it modeled.
       \item Burn-in samples of the end-node sensitivity parameter,
             \code{k}, are included in the results of \code{bart},
             \code{bart2}, and \code{rbart_vi}.
       \item \code{rbart_vi} will now look for \code{group.by} and
             \code{group.by.test} in the \code{data} and \code{test} arguments
             before looking in the \code{formula} or calling environments.
     }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Fix for \code{k} mixing across chains when running multithreaded
            and with \code{k} being modeled. Bug report thanks to Noah
            Greifer.
      \item Fix for \code{xbart} with \code{method = "k-fold"} when data
            not evenly divided by number of folds. Rug report thanks to Jesse
            (@ALEXLANGLANG on Github).
      \item Sampler method \code{getLatents} and corresponding C function now
            add user supplied offset to result.
      \item Saved, flattened trees now correctly partition observations on
            left and right.
    }
  }
}
\section{CHANGES IN VERSION 0.9-19}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Samplers now have method \code{sampleNodeParametersFromPrior}.
            When used in conjunction with \code{sampleTreesFromPrior} allow
            the model to fully make predictions from the prior distribution.
      \item \code{dbartsControl} (and now \code{bart}/\code{bart2} through
            \code{...}) now accept \code{rngSeed} argument. This can be used
            to generate reproducible results with multiple threads. It should
            only be used for testing, as the thread-specific pRNGs are seeded
            using sequential draws from a pRNG created with the user-supplied
            seed.
      \item C interface supports \code{dbarts_createStateExpression} and 
            \code{dbarts_initializeState} which can be used to re-create
            samplers that were allocated using forked multithreading.
      \item C interface also supports \code{dbarts_predict},
            \code{dbarts_setControl}, and \code{dbarts_printTrees}.
      \item Exports \code{makeTestModelMatrix} to allow package authors to
            create test data at a later point from training data.
    }
  }
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item \code{varcount} for \code{bart} fits now has dimnames set.
      \item \code{residuals} generic added to \code{bart} and \code{rbart_vi}.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Parallelization for \code{rbart} now creates the correct
            number of chains.
      \item Should now compile on non-x86 architectures. Report thanks to
            Lars Viklund.
      \item Fixed hang when \code{verbose = TRUE} for multiple threads and
            multiple chains. Report thanks to Noah Greifer.
      \item Fixed potential memory access errors when recreating sample from
            saved state.
      \item Correctly de-serializes saved tree structure.
    }
  }
}
\section{CHANGES IN VERSION 0.9-18}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Sampler now explicitly supports \code{setSigma} for use in 
            hierarchical models.
      \item Sampler function \code{setOffset} has an additional argument
            of \code{updateScale}. When the response is continuous and
            \code{updateScale} is \code{TRUE}, the implicit scaling,
            effecting the node parameters' variance, is adjusted to match
            the range of the new data. This optionally reverts the change
            of version 0.9-13 with the intention of being used only during
            warmup when using an offset that is itself being sampled.
    }
  }
  \subsection{BUG-FIXES}{
    \itemize{
      \item Extraneous print line from debugging 0.9-17.
      \item Eliminated two race conditions from multithreaded crossvalidation.
            Report thanks to Ignacio Martinez.
      \item Eliminated garbage read on construction of crossvalidation sampler,
            removing inconsistencies across multiple runs with the same starting
            seed.
      \item \code{makeModelMatrixFromDataFrame} now converts character vectors
            to factors instead of dropping them. Report thanks to Colin Carlson.
    }
  }
}
\section{CHANGES IN VERSION 0.9-17}{
  \subsection{BUG-FIXES}{
    \itemize{
      \item Memory leak for \code{predict} when \code{keepTrees} is \code{FALSE}.
    }
  }
}
\section{CHANGES IN VERSION 0.9-16}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Added \code{extract} and \code{fitted} generics for \code{bart}
            models. Respects \code{"train"} and \code{"test"} sets of
            observations while returning \code{"ev"} - samples from the
            posterior of the individual level expected value, \code{"bart"}
            - the sum of trees component; same as \code{"ev"} for linear
            models but on the probit scale for binary ones, and \code{"ppd"}
            - samples from the posterior predictive distribution. To synergize
            with \code{fitted.glm}, \code{"response"} can be used as a synonym
            for \code{"ev"} and \code{"link"} can be used as a synonym for
            \code{"bart"}.
    }
  }
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item \code{predict} for \code{bart} models with binary outcomes returns
            a result on the probability scale, not probit. The argument
            \code{value} is deprecated - use \code{type} instead.
      \item \code{predict} further conforms to the same system of arguments as
            \code{extract} and \code{fitted}.
    }
  }
  \subsection{BUG-FIXES}{
    \itemize{
      \item \code{xbart} with a k-hyperprior should no longer crash. Report thanks
             to Colin Carlson.
    }
  }
}
\section{CHANGES IN VERSION 0.9-14}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Fits from \code{rbart_vi} now work with generics \code{fitted},
            \code{extract}, and \code{predict}. \code{extract} retrieves
            samples from the posterior distribution for the training and test
            samples, \code{fitted} applies averages across those samples,
            while \code{predict} can be used to obtain values for completely
            new observations.
    }
  }
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item \code{predict} for \code{rbart_vi} takes value "ev" instead of "post-mean"
            to clarify what is being returned, i.e. samples from the posterior 
            distribution of the observation-level expected values.

    }
  }
  \subsection{BUG-FIXES}{
    \itemize{
      \item \code{save}/\code{load} should work correctly. Report thanks to Jeremy Coyle.
    }
  }
}
\section{CHANGES IN VERSION 0.9-13}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item \code{predict} now works when trees aren't saved, for use in testing
            Metropolis-Hasting proposals.
      \item The \code{offset} slot no longer changes the relative scaling of the
            response. This stabilizes predictions across iterations. For a semantic
            where the scaling does change, use \code{setResponse} instead.
    }
  }
}
\section{CHANGES IN VERSION 0.9-12}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Varying intercepts model for probit regression.
    }
  }
}
\section{CHANGES IN VERSION 0.9-10}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item A hyperpriors for \code{k} has now been implemented. Passing
            \code{k = chi(degreesOfFreedom, scale)} now penalizes small values of
            \code{k}, encouraging more shrinkage.
    }
  }
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item Hyperprior of \code{chi(1.25, Inf)} is now default for \code{bart2}
            with binary outcomes. The default accuracy should improve substantially.
    }
  }
  \subsection{BUG-FIXES}{
    \itemize{
      \item \code{xbart} divides data correctly with random subsampling.
    }
  }
}
\section{CHANGES IN VERSION 0.9-9}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item More control over cut points has been added. It is now possible to specify
            the cut points for a variable once and subsequently change that predictor
            without also modifying the cuts using \code{sampler$setCutPoints} and
            \code{sampler$setPredictor}.
      \item \code{sampler$getTrees} implemented to get a flattened, depth-first down
            left traversal of the trees.
    }
  }
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item For \code{sampler$setPredictor}, an argument specifies whether or not to
            rollback or force the change if the new data would result in a leaf
            having 0 observations.
      \item \code{pdbart} and \code{pd2bart} now work with formula/data specifications,
            as well as taking models or samplers that have previously stored trees.
    }
  }
  \subsection{OPTIMIZATIONS}{
    \itemize{
      \item Stores \code{x} as integer matrix of the max of which cut point an observation is
            to the left of, by default using 16 bit integers. Limited to 65535 cut points.
            That can be increased with some special compilation instructions.
      \item Uses CPU dispatch and SIMD instructions for some operations. This and the integer
            \code{x} make BART about 30\% faster on datasets of around 10k observations.
      \item Saved trees are stored using significantly less memory.
    }
  }
}
\section{CHANGES IN VERSION 0.9-8}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item \code{plot} now works for fits from \code{rbart_vi}.
    }
  }
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item \code{rbart_vi} new reports \code{varcount}.
      \item \code{bart2} now defaults to not storing trees due to the memory cost.
      \item \code{bart2} now defaults to using quantile rules to decide splits.
    }
  }
  \subsection{BUG-FIXES}{
    \itemize{
      \item \code{predict} for binary outcomes now correct.
      \item Fix for verbose multithreading on Linux, reported by @ignacio82 on github.
      \item General improvements to slice sampler in \code{rbart_vi} thanks to reports from Yutao Liu.
      \item \code{sampler$plotTree} now handles multiple chains correctly.
      \item Negative log loss for \code{xbart} with binary outcomes should now be computed correctly.
    }
  }
}
\section{CHANGES IN VERSION 0.9-2}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item \code{rbart_vi} fits a simple varying intercept, random effects model.
    }
  }
}
\section{CHANGES IN VERSION 0.9-0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Now natively supports multiple chains running in parallel.
      \item Objects fit by \code{bart} can be used with the predict generic
      when instructed to save the trees.
      \item New function \code{bart2} introduced, similar to \code{bart} but with
      more efficient default parameters.
    }
  }
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item \code{dbartsControl} has had two parameters renamed: \code{numSamples}
      is now \code{defaultNumSamples} and \code{numBurnIn} is now
      \code{defaultNumBurnIn}.
      \item \code{dbartsControl} supports parameters \code{runMode},
      \code{n.chains}, \code{rngKind} and \code{rngNormalKind}.
      \item In the C interface, a new function (\code{setRNGState}) has been
      added to specify the states of the random number generators, of which there
      is now one for every chain.
      \item State objects saved by the handles no longer contain the total
      fits, since they can be rebuild from the tree fits. States are also lists
      of objects now, with one corresponding to each chain. Tree fits and strings
      are matrices corresponding to the number of trees and saved samples.
    }
  }
}
\section{CHANGES IN VERSION 0.8-6}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item random subsampling crossvalidation (\code{xbart}) has been implemented
      in C++. Refits model using current set of trees for changes in
      hyperparameters \code{n.trees}, \code{k}, \code{power}, and \code{base}.
      Natively parallelized.
      \item Rudimentary tree plotting added to sampler (\code{sampler$plotTree}).
      \item Exported \code{dbartsData} as a way of constructing data objects
      and setting the data seen by the sampler all at once. Sampler now supports
      \code{sampler$setData()}.
    }
  }
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item \code{keepevery} argument to \code{bart} matches \code{BayesTree}.
      \item \code{bart} now has argument \code{keepcall} to suppress
      storing the call object.
      \item \code{bart} now accepts a \code{weights} argument.
      \item \code{MakeModelMatrixFromDataFrame} now implemented in C, supports
      an argument for tracking/keeping dropped values from factors.
    }
  }
  \subsection{BUG-FIXES}{
    \itemize{
      \item Usage of weights was causing incorrect updates to posterior for
      \eqn{\sigma^2}.
      \item Should now JIT byte compile correctly.
      \item Cuts derived from quantiles should now be valid.
    }
  }
}    
\section{CHANGES IN VERSION 0.8-4}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Uses a rejection sampler to simulated binary latent variables
      (CP Robert 2009, http://arxiv.org/pdf/0907.4010.pdf). Code thanks to
      Jared Murray.
      \item Now encapsulates its own random number generator, so that the C++
      objects can safely be used in parallel. Shouldn't affect pure-R users 
      unless their RNG has non-exported state (i.e. Box-Muller normal kind).
      \item Includes a \code{offset.test} vector that can be controlled 
      independently of the \code{offset} vector, but in general inherits
      behavior from it. Set at creation with \code{dbarts()} or after
      with \code{setTestOffset} or \code{setTestPredictorAndOffset}.
    }
  }
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item By default, no longer attempts to obtain identical results as
      BayesTree. To recover this behavior, compile from source with 
      \code{configure.args = "--enable-match-bayes-tree"}.
      \item Changing the entirety of the test matrix using \code{setTestPredictor}
      no longer allowed. Use \code{setTestPredictors} instead.
      \item Changing the predictor can now result in failure if the covariates
      would leave an end-node empty. \code{setPredictor} returns a logical
      as to success.
      \item Saved \code{dbarts} objects may not be compatible and should be
      re-created to be sure of valdity.
      \item Now requires R versions >= 3.1.0.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Corrected binary latent variable sampler and no longer multiply
      adds offset (reported by Jared Murray).
      \item Relatively embarassing bug related to loop-unrolling when \code{n mod 5 != 0}
      fixed.
      \item Correct aggregation of results for multithreaded variance calculations.
      \item More equitably distributed tasks across multiple threads.
      \item Makevars tweaked to allow compilation on Ubuntu.
    }
  }
}
\section{CHANGES IN VERSION 0.8-3}{
  \itemize{
  	\item Initial public release.
  }
}
