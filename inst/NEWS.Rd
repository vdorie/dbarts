\name{NEWS}
\title{dbarts News}
\encoding{UTF-8}
\section{CHANGES IN VERSION 1.0-0}{
  \subsection{UPGRADING}{
    A summary of the changes that can break existing code or alter
    results; each is detailed under the other headings.
    \itemize{
      \item dbarts now requires R 4.2.0 or newer; the sampler core is C++20.
      \item Results differ from 0.9-x. The rewritten engine draws a different
            stream, and several corrections move the posterior itself: the
            tree change move now satisfies detailed balance, the initial
            forest is drawn from the tree prior conditioned on carrying no
            empty leaf, the residual-variance degrees of freedom no longer
            count zero-weight observations, and the \code{chi} hyperprior on
            \code{k} samples the degrees of freedom asked for. Seeded
            multi-chain results no longer depend on the thread count, and
            unseeded ones follow \code{set.seed}. Sampling no longer advances
            R's stream (chain generators are seeded from it at creation), so
            R-level draws made after a fit also change.
      \item The binary default for \code{k} on \code{dbarts}, \code{bart2},
            and \code{rbart_vi} is the proper \code{chi(1.5, 2)} rather than
            the improper 0.9-x default, spelled \code{chi(1.25, Inf)} there
            and, under the corrected degrees-of-freedom count above, the
            same law as \code{chi(1.5, Inf)} now. \code{chi()}'s own
            \code{scale} argument defaults to \code{2} rather than
            \code{Inf}, so an explicit \code{chi(df)} moves as well. Pass
            \code{k = chi(1.5, Inf)} for the old prior. \code{xbart} instead
            defaults to the fixed value \code{k = 2} on every family, where a
            binary response previously inherited the hyperprior.
      \item Sampler states saved by earlier versions cannot be restored;
            refit anything you need to continue.
      \item A \code{dbartsData} object saved by an earlier version keeps the
            design and column types it was built with; loading one does not
            upgrade it, and it fits differently from an object freshly built
            over the same data frame. There is no upgrade path for an ordered
            factor in particular: a column's type rides the \code{varTypes}
            slot, and the code an ordered factor carried before this version
            is the code a numeric column carries, so nothing in the stored
            object tells the two apart. Rebuild any saved data object with
            \code{dbartsData} to fit under the current coding.
      \item Packages compiling against dbarts must port to the flat C API in
            \code{dbarts/dbarts.h}. The C++ headers -
            \code{dbarts/R_C_interface.hpp}, \code{dbarts/bartFit.hpp} and
            the rest of that tree - are removed with the engine they exposed.
      \item \code{bart2} and \code{rbart_vi} now default to
            \code{combineChains = TRUE}, so \code{yhat.train} and the other
            draw channels of a multi-chain fit come back with the chain and
            sample margins merged rather than separate. Pass
            \code{combineChains = FALSE} for the previous shape.
      \item Unordered factors are split on level subsets by default in
            \code{dbarts}, \code{dbartsData}, \code{bart2}, \code{xbart},
            and \code{rbart_vi}, and an ordered factor becomes a single
            ordinal column split at a threshold, where 0.9-x expanded a
            factor of either kind into indicator columns. The design
            matrix, the \code{varcount} column names, and the draws all
            move for a fit with any factor predictor. \code{factors} is a
            new argument;
            pass \code{factors = "indicators"} for the previous expansion.
            \code{bart} does not take it and expands as before.
      \item 0.9-x silently dropped every row with a missing predictor, on
            the formula and matrix interfaces alike. \code{missing}, a new
            argument on \code{dbarts}, \code{dbartsData}, \code{bart2},
            \code{xbart}, and \code{rbart_vi}, defaults to
            \code{"incorporate"}: the rows are kept and their missingness is
            modeled (Missingness Incorporated in Attributes; Twala et al.
            2008). \code{missing = "error"} rejects them instead; neither
            setting drops anything. \code{bart} takes no \code{missing}
            argument and refuses incomplete predictors outright. An
            \code{NA} in test predictors or \code{newdata} is an error
            wherever the corresponding training column was complete: a split
            rule learns a missing-direction route only where its own column
            had training \code{NA}s, so such a value had no route to take and
            was previously answered silently.
      \item Where the training and test matrices both carry column names and
            the design has more than one predictor, test column names that
            do not cover the training design's are an error naming the
            missing columns; previously the mismatch warned and the columns
            were matched by position.
      \item An unordered factor response with three or more levels is
            multinomial under \code{family = "auto"}: \code{bart2} fits
            \code{family = "multinomial"}, and \code{bart}, \code{dbarts},
            \code{rbart_vi}, and \code{xbart} refuse it, naming
            \code{bart2(family = "multinomial")}. An ordered one is ordinal
            instead: \code{dbarts} and \code{bart2} fit
            \code{family = "ordinal"}, and the other three refuse it, naming
            \code{bart2(family = "ordinal")}. Passed
            through the \code{x}/\code{y} interface either previously fit a
            continuous response on the integer level codes.
      \item \code{dbartsControl} loses \code{rngKind} and
            \code{rngNormalKind}, and renames \code{rngSeed} to \code{seed}.
      \item \code{xbart} loses its \code{control} argument and renames
            \code{sigma} to \code{sigest}. \code{n.cuts},
            \code{useQuantiles}, and \code{n.thin} are ordinary formals,
            spelled and defaulted exactly as \code{dbartsControl}'s own.
            \code{n.burn} takes two values rather than three: the
            between-replication burn-in is gone, and a length greater than
            two is refused by name instead of the third (and any further)
            element being silently dropped - 0.9-34's own default was three
            elements and honored the third as the per-replication burn-in.
      \item A \code{probit} fit refuses \code{weights}. 0.9-x accepted them
            and fit a weighted probit, which was incorrect and is not carried
            forward; a vector identically 1 is treated as absent, and integer
            count weights belong on \code{family = "logistic"}.
      \item On the \code{x.train}/\code{y.train} interface a \code{weights}
            vector whose length is neither the number of observations nor
            one is an error instead of being recycled; the formula interface
            always refused it. \code{resid.prior = fixed(value)} now holds
            the residual variance at \code{value}, at creation and through
            \code{setModel} alike; previously a fixed prior given at
            creation was ignored, sigma sitting at its starting estimate
            for the whole run, and one installed through \code{setModel}
            was read as a standard deviation.
            \code{getSumsOfSquaredResiduals} returns the sum of
            squares itself, where it previously returned that sum divided by
            the response range.
      \item The sampler's mutators - \code{setData}, \code{setResponse},
            \code{setOffset}, \code{setWeights}, \code{setSigma},
            \code{setPredictor}, \code{setCutPoints} - refresh the local
            \code{state} cache only on an explicit \code{updateState = TRUE}.
            \code{NA} no longer defers to the control's own
            \code{updateState}, so the sampler's \code{$state} - and with it
            a fit written by \code{saveRDS} after the usual
            \code{invisible(fit$fit$state)} read - does not reflect a
            mutation made after \code{run} unless that mutator was passed
            \code{updateState = TRUE}. \code{run},
            \code{sampleTreesFromPrior}, and
            \code{sampleNodeParametersFromPrior} keep the old convention.
      \item The sampler's \code{run} renames its \code{numThreads} argument
            to \code{n.threads}; the old spelling is an unused-argument
            error.
      \item The sampler's \code{setResponse} gains an \code{updateScale}
            argument defaulting to \code{FALSE}: the response transform -
            the location and range that map the response into the tree
            prior's scale - stays as it was set at creation when the
            response is swapped, as a conditional update inside a larger
            sampler requires. A continuous-response sampler previously
            re-derived that transform from every new response, while a
            binary one re-derived nothing and the residual-variance prior's
            scale was kept in the original units either way; pass
            \code{updateScale = TRUE} for the old behavior.
            \code{updateScale} is the second positional argument and
            displaces \code{updateState} to third, so
            \code{setResponse(y, TRUE)} now sets \code{updateScale} where it
            used to set \code{updateState} - pass \code{updateState} by name.
            Supplying the second argument positionally now also warns once
            per session, naming the reorder, classed
            \code{dbartsPositionalArgsWarning}.
      \item An \code{rbart_vi} fit built with a built-in \code{prior} and no
            \code{callback} - the default - carries a length-one \code{$fit}
            list holding one multi-chain sampler, rather than one
            single-chain sampler per chain; code indexing \code{$fit[[i]]} by
            chain must be reworked. A custom prior or callback keeps the
            per-chain list.
      \item \code{predict.bart} and \code{predict.rbart} now share their
            first five formals, \code{(object, newdata, type, offset,
            weights)}. On \code{predict.bart}, \code{type} moves ahead of
            \code{offset} and \code{weights}, which previously preceded it,
            and \code{n.threads} moves from seventh to tenth. On
            \code{predict.rbart}, \code{type} moves from fifth to third,
            ahead of \code{offset}; \code{weights} and \code{n.threads} are
            new; and \code{group.by} moves after \code{\dots}, so it is
            matched by name only - passed positionally or omitted, it is
            refused by name rather than silently binding \code{type}.
            \code{offset.test} on either method is refused by name rather
            than dropped into \code{\dots}.
      \item \code{fitted}'s positional third argument is \code{ci.level} on
            \code{bart} and \code{rbart} fits, where it was \code{sample};
            \code{sample} moves to the fourth position.
      \item \code{predict.rbart}'s \code{value} argument - a deprecated
            \code{\dots} alias for \code{type}, never a formal - and its
            \code{type = "post-mean"} alias are removed. \code{value} is
            refused by name, pointing at \code{type}; \code{"post-mean"}
            falls into \code{type}'s own \dQuote{must be in} message, which
            lists the types it does take. Use \code{type = "ev"}.
      \item These names, foreign to the method actually called, are refused
            by name instead of being silently discarded through
            \code{\dots}: \code{sample} and \code{group.by} on
            \code{predict.bart}; \code{newdata}, \code{offset},
            \code{weights}, \code{n.threads}, and \code{group.by} on
            \code{extract}; the same five plus \code{sample} and
            \code{combineChains} on \code{fitted}; and all of those on
            \code{residuals}. \code{extract(type = "trees")}
            refuses \code{sample} and \code{combineChains}, where a numeric
            \code{sample} previously partial-matched the sampler's
            \code{getTrees(sampleNums)} and turned the train/test selector
            into a draw selector. An extra argument passed positionally is
            refused alongside the named ones.
      \item A fractional double is refused, naming the argument, rather than
            silently truncated on these count arguments: \code{n.samples},
            \code{n.trees}, \code{n.burn}, \code{n.thin}, \code{n.cuts},
            \code{n.chains}, \code{n.threads}, \code{printEvery},
            \code{printCutoffs}, and \code{seed} on \code{dbartsControl},
            \code{bart2}, and \code{rbart_vi}; \code{n.trees},
            \code{n.burn}, \code{n.thin}, \code{n.cuts}, \code{n.threads},
            \code{n.reps}, and, under \code{method = "k-fold"},
            \code{n.test} on \code{xbart}; and \code{bart}'s own lowercase
            spellings \code{ntree}, \code{ndpost}, \code{nskip},
            \code{nchain}, \code{nthread}, \code{keepevery}, \code{numcut},
            \code{printevery}, and \code{printcutoffs}. The rule also
            reaches these further counts and column/group selectors a
            caller can write directly:
            \code{dbarts} and \code{xbart}'s own \code{n.samples};
            \code{dbarts}'s single-value \code{breaks} (a hazard period
            count); \code{varianceForest}'s \code{n.trees}; \code{bart2}'s
            \code{n.grow.sweeps}; \code{bart}'s own
            \code{seed} (previously pre-truncated before reaching
            \code{dbartsControl}'s own check); the \code{forest} argument
            shared by \code{$getCalibration}, \code{$setCalibration},
            \code{$getForestFits}, \code{$getForestAmplitudes},
            \code{$getForestVariableCounts}, \code{$setForestWeights},
            \code{$setForestBasis}, and \code{extract}/\code{predict}'s
            own \code{forest} selector; \code{plotTree.rbart}'s
            \code{chainNum}; the \code{column} argument shared by
            \code{$setPredictor}, \code{$setCutPoints},
            \code{$setTestPredictor}, and
            \code{updatePredictorPerObservationJointly}; \code{linear} and
            \code{gp}'s numeric \code{columns}; \code{varianceForest} and
            \code{forest}'s \code{vars} selector; \code{forest}'s
            \code{n.trees}; \code{interactions}'s \code{max.order},
            \code{forbid}, and \code{groups}; \code{blocks}'s
            \code{trees.per.group} and \code{groups};
            \code{extract.rbart}'s own
            \code{chainNums} selector (a distinct meaning from
            \code{$getTrees}'s); and, on the sampler,
            \code{$run}'s \code{numBurnIn}/\code{numSamples},
            \code{$growFromRoot}'s \code{n.sweeps}, \code{$installTrees}'s
            \code{samples}, and \code{$getTrees}'s/\code{$printTrees}'s
            \code{treeNums}/\code{chainNums}/\code{sampleNums}. The sampler's
            \code{$getSigmas} and \code{$getSumsOfSquaredResiduals} now
            refuse their vestigial \code{result} argument by name instead of
            silently accepting and ignoring it; both readers always allocate
            their own vector.
    }
  }
  \subsection{NEW FEATURES}{
    \itemize{
      \item The rewritten engine replaces the pre-1.0 engine: every
            fitting function runs on it. Fits are not draw-for-draw
            identical to the pre-1.0 engine's, and the corrections listed
            under BUG FIXES move the posterior itself, so seeded results
            change. \code{bart} keeps its
            historical behavior apart from the engine (indicator
            expansion of factors, probit for binary responses).
            \code{dbartsControl} loses the \code{rngKind} and
            \code{rngNormalKind} arguments; explicit generator kinds were a
            pre-1.0-engine feature. Saved sampler states from previous
            versions cannot be restored.
      \item Sampling is faster: the benchmark run scenarios measure 25 to
            30 percent faster per iteration than the pre-1.0 engine,
            the gap widening with the number of observations. The run
            loop keeps a running residual across each sweep, and change
            and swap proposals score the tree prior on the modified
            subtree alone. For large problems, \code{storage = "single"}
            in \code{dbartsControl}, \code{bart}, \code{bart2},
            \code{rbart_vi}, and \code{xbart} additionally stores
            the running residual in single precision (fp32) while
            accumulating every reduction in double, roughly halving the
            memory traffic of the dominant per-sweep sufficient statistic
            pass; it is opt-in and currently supported only for gaussian
            responses with constant leaves (an error is raised otherwise),
            and the leaf draws stay double precision, so posteriors are
            statistically indistinguishable from the default
            double-precision fit. The speedup is largest where memory
            bandwidth binds - large \eqn{n}{n} and multi-chain runs -
            roughly 10-30\% at \eqn{n = 10^6}{n = 1e6} in testing, growing
            with the number of chains. The default \code{storage = "double"}
            changes nothing: it is the path every other fit takes. Separately, the
            internal per-observation index storage was narrowed, lowering
            sampler memory use with no change to results. See
            \code{?dbartsControl}.
      \item The rewritten engine implements the full \code{dbartsSampler}
            interface: every method works, including \code{keepTrees},
            \code{getTrees}, \code{predict}, \code{plotTree}, \code{setWeights}, test offsets,
            \code{printTrees}, sampling from the prior, \code{setControl},
            \code{setModel}, verbose progress output, and saving/loading
            sampler state (a restored sampler continues the same chain,
            reproducing an uninterrupted run to within floating-point
            reassociation rather than bit for bit; the stored state is opaque
            and incompatible with states saved by previous versions). Not
            supported: weights with a \code{probit} response, and
            \code{setControl} changes to settings fixed at creation
            (chain/tree counts, generators, and \code{useQuantiles}; a changed
            \code{n.cuts} is ignored, the grid being fixed with the data).
            \code{\link{extract.dbartsSampler}} returns the
            predictor matrix as the sampler currently sees it - the same
            numeric coding \code{getTrees} replays saved trees against -
            reflecting any \code{setPredictor} mutation since creation.
      \item Every chain runs its own generator. \code{seed} makes results
            reproducible for any number of chains without touching R's
            stream, seeded results no longer depend on the thread
            count, and a single-chain run with a given seed reproduces
            chain 0 of a multi-chain run with the same seed. Without a
            seed, chain generators are seeded from R's stream at
            creation, so \code{set.seed} alone suffices for reproducibility
            and sampling itself never advances R's stream (the pre-1.0
            engine ran through R's generator whenever a run was single-chain
            or single-threaded, and otherwise, absent an explicit
            \code{rngSeed}, seeded chains from the clock).
      \item \code{\link{rbart_vi}} with a built-in \code{prior}
            (\code{cauchy}, \code{gamma}) and no \code{callback} - the
            default - runs its whole Gibbs sampler, the random intercepts and
            their standard deviation included, inside the engine: one sampler
            with \code{n.chains} chains on worker threads, in place of
            per-chain samplers driven iteration by iteration from R over a
            PSOCK or FORK cluster. A custom prior function or a
            \code{callback} keeps the R loop. Results are statistically
            equivalent to that loop's but not draw for draw identical: the
            response transform is fixed at creation from the data alone,
            where the loop re-anchors it during warm-up with the current
            random effects folded into the offset, and with case weights the
            conjugate intercept update weights each group's contribution,
            where the loop's takes unweighted group means. Fitted objects keep
            their shape, and under \code{keepSampler} the \code{fit} element
            is a length-one list holding that one multi-chain sampler rather
            than one sampler per chain.
      \item \code{family} is a new argument on \code{dbarts}, \code{bart},
            \code{bart2}, \code{rbart_vi}, and \code{xbart}, and the new
            \code{dbartsSpec} takes it too, each over the tokens it can fit.
            \code{"auto"} keeps
            the existing dispatch (gaussian for continuous responses, probit
            for 0/1); \code{"gaussian"} forces a continuous fit even for a 0/1
            response; and \code{"logistic"} fits a logistic latent model via
            Polya-Gamma augmentation. \code{bart}'s own set is \code{"auto"},
            \code{"logistic"}, and \code{"aft"}, and \code{rbart_vi}'s is
            \code{"auto"}, \code{"gaussian"}, and \code{"aft"}.
            \code{xbart}'s built-in binary losses transform test
            predictions through the resolved family's link.
            A logistic
            model accepts positive-integer count weights, a weight-w
            observation contributing its likelihood to the w-th power;
            a probit refuses weights (the weighted probit the pre-1.0
            engine fit was incorrect and is not carried forward), except
            that weights identically 1 are treated as absent, as
            SuperLearner-style callers pass unit weights unconditionally.
            A logistic sampler's counts can be replaced between runs with
            \code{$setWeights} or \code{$setData}, which redraw the
            Polya-Gamma latents against the new counts before returning,
            so an outer sampler can vary exposure from sweep to sweep.
            \code{probit}, \code{ordinal}, \code{nbinom}, and
            \code{multinomial} take no case weights, unit weights aside, and
            \code{aft} none at all; none of the five accepts a weight change
            after creation, on any conduit.
            Fits returned by \code{bart}, \code{bart2}, and \code{rbart_vi} record the
            resolved family in a \code{family} element, and \code{predict},
            \code{extract}, \code{fitted}, and \code{plot} transform latent draws
            through that family's link; fits saved by earlier versions
            lack the element and are treated as probit. \code{family = "auto"}
            extends past probit and gaussian to detect a factor, logical,
            or character response directly: two levels resolve to probit,
            three or more ordered levels to \code{"ordinal"} (\code{dbarts}
            and \code{bart2}; see below), and three or more unordered levels
            to multinomial (\code{bart2} only; elsewhere an error names the
            limitation), each reporting the choice in a one-line message.
            Previously the x/y interface coerced a
            factor response to its zero-based level codes and fit those: at
            two levels its own 0/1 binary detection then fit them as probit -
            the same choice made here - and at three or more it fit them as a
            continuous response. A logical or character response the x/y
            interface refused outright; the formula interface refused a
            factor response and coerced a logical one to 0/1, which its own
            binary detection then fit as probit.
      \item Ordinal (ordered categorical) responses:
            \code{family = "ordinal"} in \code{dbarts} and \code{bart2}
            fits a cumulative probit - probit's truncated-normal latents
            generalized to K ordered categories split by sampled
            thresholds (\eqn{\gamma_1 = 0}{gamma_1 = 0} pinned; the free
            thresholds move by a marginal Metropolis update with the
            latents integrated out, under a log-gap
            \eqn{N(0, 1.5^2)}{N(0, 1.5^2)} prior). At K = 2 the model is
            probit exactly. An ordered factor response with three or
            more levels is auto-detected under \code{family = "auto"},
            announced; the level order defines the category order.
            \code{bart2} returns a \code{"bartOrdinal"} fit carrying the
            levels, the n x K category-probability draws, the latent
            draws, and the threshold draws, with \code{fitted},
            \code{extract}, \code{predict}, and \code{residuals}
            methods. Weights are not supported; \code{rbart_vi} and
            \code{xbart} refuse ordinal responses (grouped ordinal
            models and ordinal cross-validation losses are not
            supported). See \code{?bart2}.
      \item Negative-binomial count responses:
            \code{family = "nbinom"} in \code{dbarts} and \code{bart2}
            fits non-negative integer counts by the Polya-Gamma
            negative-binomial augmentation (the forest fits a log-odds
            latent \eqn{\psi = f(x) + o}{psi = f(x) + o}; the count law
            is \eqn{y \sim NB(r, \mathrm{plogis}(\psi))}{y ~ NB(r,
            plogis(psi))} with mean \eqn{r e^{\psi}}{r exp(psi)}, so the
            offset is a log-exposure). The response must be a
            non-negative integer count, and \code{family = "nbinom"} is
            always explicit - a count carries no unambiguous class to
            auto-detect. The dispersion \code{r} is estimated by default
            on a capped positive-integer grid
            (\{1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 30, 50\}) under a
            renormalized
            \eqn{\mathrm{gamma}(2, 0.1)}{gamma(2, 0.1)} prior; supplying
            \code{dispersion} as a positive integer fixes it (v1 ships
            the exact integer envelope, so a real fixed \code{dispersion}
            is refused). The dispersion in force at each draw is a
            first-class channel: \code{run()} gains a \code{dispersion}
            element, shaped exactly as \code{sigma}, which is how
            \code{bart2(family = "nbinom")} itself reports it, and
            \code{$getDispersion()} reads the value mid-sweep - a pure
            read of state the sweep already settled on, consuming no
            random numbers. Both are \code{NULL} on every other
            family, so a caller tests \code{!is.null()} rather than
            comparing. \code{bart2} returns a \code{"bartNegbin"} fit
            carrying the mean-count draws, the latent draws, and the
            per-draw dispersion, with \code{fitted}, \code{extract},
            \code{predict}, and \code{residuals} methods. Weights are not
            supported (exposure belongs in the offset); \code{rbart_vi}
            and \code{xbart} do not take the \code{"nbinom"} token (grouped
            negative-binomial models, real dispersion, and a Poisson
            family are not supported). See \code{?bart2}.
      \item \code{dbarts}, \code{bart}, \code{bart2}, and \code{rbart_vi}
            accept \code{family = "aft"}: an
            accelerated failure time, log-normal survival model (the forest
            models the conditional mean of log survival time, with \code{sigma} the
            log-scale residual standard deviation). The response is a
            \code{survival::Surv} object (right-censoring only) or a two-column
            \code{(time, status)} matrix or data frame; under \code{family = "auto"} a
            \code{Surv} response selects \code{"aft"} on its own. \code{rbart_vi} adds a
            random intercept per group to the same model, through its formula
            interface; \code{dbarts}, \code{bart}, and \code{bart2} take the
            response through their x/y interface. None of them support
            \code{weights} or \code{subset} for a survival fit.
            \code{predict}, \code{extract}, and \code{fitted} return the
            log-time-scale linear predictor, never the time scale; the
            \code{\link{survivalProbabilities}} generic draws the survival probability
            over a user-supplied time grid instead. The same response also
            fits as a discrete-time hazard model:
            \code{family = "hazard"} (probit link) and
            \code{family = "hazard.logistic"} in \code{dbarts} and \code{bart2}
            expand each subject into its at-risk rows with an ordinal
            period column and fit the discrete hazard
            \eqn{h(k \mid x) = g(f(x, k) + o)}{h(k | x) = g(f(x, k) + o)}
            as an ordinary binary model on the expanded rows - so a hazard
            fit is byte-identical, draw for draw, to the binary fit on the
            hand-expanded design. The time grid defaults to the distinct
            observed times; \code{breaks} coarsens it (an integer quantile
            count or an explicit boundary vector) and \code{max.rows}
            guards the expansion size. The fit is an ordinary
            \code{"bart"} object whose \code{$family} records the binary
            link and whose \code{$periods} carries the grid;
            \code{survivalProbabilities} again produces the survival-curve
            draws (requiring \code{keepTrees}). The matrix interface is
            required for a hazard fit, and weights follow the chosen
            binary family's policy. \code{rbart_vi} refuses the hazard
            forms and \code{xbart} refuses both (grouped or frailty hazard, a
            cloglog link, and cross-validation losses for either are not
            supported).
      \item Semicontinuous two-part (hurdle) responses:
            \code{family = "hurdle.lognormal"} (alias \code{"twopart"}) in
            \code{bart2} fits a non-negative response with
            exact zeros by composing, R-side, an occupancy probit fit of
            \eqn{z = 1\{y > 0\}}{z = 1{y > 0}} over all n observations with a
            lognormal positive-part fit - an ordinary gaussian fit of
            \eqn{\log y}{log y} - over the \eqn{y > 0} subset, glued at
            report time. The two parts share no parameters (no engine code
            is added) and fit at independently derived seeds. Predictions
            default to the natural scale,
            \eqn{E[y \mid x] = P(y > 0 \mid x) e^{f(x) + \sigma^2/2}}{
            E[y | x] = Pr(y > 0 | x) * exp(f(x) + sigma^2/2)}, formed at each
            posterior draw from the positive part's single
            \eqn{\sigma}{sigma} for that draw; \code{type = "prob"} gives the
            occupancy probability, \code{"link"}/\code{"log"} the positive
            part's log-scale linear predictor, and \code{"ppd"} a proper
            bimodal predictive (a point mass at zero, else a lognormal
            draw). \code{bart2} returns a \code{"bartHurdle"} fit carrying
            both component fits, with \code{fitted}, \code{extract},
            \code{predict}, and \code{residuals} methods. Fits currently use
            the matrix interface only; \code{weights}, \code{subset},
            \code{offset}, and \code{test} are refused (the positive part is
            instead given the full training \code{x} as its own
            \code{x.test}). \code{dbarts} does not fit this family - it
            composes two samplers, which only \code{bart2} builds - and
            \code{rbart_vi} and \code{xbart} do not fit it either (a shared
            variable-selection prior across the two parts, Duan smearing,
            and a heteroscedastic positive part are not supported: the
            positive part is always homoscedastic).
            See \code{?bart2}.
      \item \code{bart2} fits \code{family = "multinomial"}, and
            \code{dbarts} builds the same sampler directly, as does
            \code{dbartsSpec} for the flat C entry point: a K-category softmax
            classifier, K forests coupled through an interleaved Polya-Gamma
            one-vs-rest augmentation. The response is a factor (or character)
            with three or more levels, an n x K matrix of nonnegative integer
            trial counts (one column per category), or the equivalent
            \code{cbind(c1, ..., cK) ~ x} formula; under \code{family = "auto"} a factor
            (or character) response with three or more levels is detected and
            routed here, reporting the choice in a one-line message. The
            fit's class is \code{"bartMultinomial"}, with its own \code{residuals},
            \code{plot}, \code{summary}, \code{print}, and \code{predict} (\code{type = "ev"}/\code{"ppd"})
            methods, separate from \code{bart}'s. The sampler surface a
            multinomial fit's \code{$fit} exposes is its own:
            \code{$setCounts} replaces the n x K count matrix (trials
            \eqn{n_i = \sum_k \mathrm{counts}_{ik}}{n_i = sum_k counts[i, k]}
            at least 1), \code{$setCategoryOffset} installs or clears an
            n x K offset entering the log-sum-exp margins before the
            softmax - the response-side counterpart of \code{setCounts},
            distinct from the flat, softmax-null \code{setOffset} - and
            \code{$setCategoryTestOffset} does the same for the recorded
            test channel. \code{n} and \code{K} are fixed at creation; only
            the values change, and the trees carry over, fitted to the
            previous counts, exactly as \code{setResponse} leaves a
            single-forest sampler's. The predictor family
            (\code{$setPredictor}, \code{$setCutPoints},
            \code{$setTestPredictor}) and the global \code{$setActiveRows}
            stay open; the rest of the mutation surface is
            refused by name, each refusal naming the channel that serves the
            caller instead. \code{predict} on a fitted model takes its
            \code{offset} as the per-category shift at the predicted rows:
            an \code{nrow(newdata)} x K matrix entering the raw fits before
            the softmax, in the same layout as the training \code{offset}. A
            fit trained with an \code{offset} requires it - the predicted rows
            are not the training rows, so no resident offset describes them -
            and passing the training offset back at the training rows
            reproduces \code{yhat.train}.
      \item Robust (Student-t) residuals: \code{resid.dist = student(df)}
            in \code{dbarts}, \code{bart}, and \code{bart2} fits a
            continuous gaussian response with outlier-robust Student-t
            errors by the classic Gaussian scale-mixture augmentation
            (per-observation latent precisions, drawn conjugately each
            sweep, that downweight gross outliers). \code{student(df =
            nu)} fixes the degrees of freedom (\code{nu = 4} a
            conventional robust default); \code{student()} estimates them
            on a capped grid under a proper prior. The default remains
            \code{gaussian()}. Only continuous gaussian responses carry
            it; the binary and survival families keep their fixed latent
            scale. Note that \eqn{\sigma}{sigma} is then the conditional
            scale (marginal residual variance
            \eqn{\sigma^2\nu/(\nu-2)}{sigma^2 nu / (nu - 2)}), and that a
            flexible tree-forest mean partially confounds tail inference,
            so an estimated \eqn{\nu}{nu} should be posterior-checked. The
            degrees of freedom in force at each draw are recorded on the fit
            as \code{$resid.df} (\code{run()$resid.df} on a sampler, and
            \code{dbarts_results::residualDf} in the C API), and its
            pointwise log-likelihood scores the marginal
            \eqn{t_\nu}{t_nu} density from it rather than the gaussian one.
      \item Multi-forest models: \code{dbarts} and \code{dbartsSpec} accept a
            \code{forests} list of \code{\link{forest}} specifications, which
            fits the mean as a weighted sum of the declared ensembles instead
            of a single one. Each forest carries its own tree count, structure
            prior, leaf scale, column restriction (\code{vars}) and
            \code{interactions}/\code{blocks} constraints, and a forest given
            a \code{basis} enters multiplied by the amplitudes on that basis's
            columns. Any number of forests may be declared, each with a
            basis of any width - a factor of any width expands to one
            amplitude per level, and a numeric vector or matrix is already
            those columns - with the two-forest, two-level-factor case being
            the Bayesian causal forest
            \eqn{y = a \mu(x) + b_z \tau(x) + \epsilon}{y = a mu(x) +
            b_z tau(x) + eps}. Every forest past the first needs a
            \code{basis}, the amplitudes multiplying it being what distinguish
            it from the first. The result is an ordinary
            \code{dbartsSampler}: it stores, saves, and re-creates like any
            other, and the columns a basis expands to ride the data object
            beside the weights they mirror. Gaussian, \code{"probit"}, and
            \code{"logistic"} responses; under the latter two the forests
            combine into the index rather than the mean, on the link's own
            fixed scale, so sigma is pinned and every forest's \code{sd} is
            stated in latent standard deviations rather than
            \eqn{sd(y)}{sd(y)}; and, under every family, each forest's
            \code{sd} is stated per unit of basis row norm, the median
            of that basis's nonzero row norms dividing it out, so an
            unnormalized basis does not multiply the prior. \code{"aft"},
            \code{"ordinal"}, \code{"nbinom"} and \code{"multinomial"} are
            refused at creation, naming what each is missing, as are the
            options a multi-forest model does not read - among them
            \code{monotone}, a variance forest, Student-t residuals, a DART
            tree prior, a linear or Gaussian-process node prior,
            \code{storage = "single"}, and test predictors - each named
            rather than dropped in silence.
            The prior defaults a forest takes
            when it says nothing are stated in terms a caller can defend: a
            forest with no \code{basis} takes a half-Cauchy median of 1
            rather than 2 under \code{"probit"} and \code{"logistic"}
            (unchanged at 2 under a gaussian response, where a drawn
            \code{sigma} absorbs the difference), putting the prior on the
            latent index where the single-forest binary default already
            puts it, and a forest with a \code{basis} takes
            \eqn{\sqrt{2/K}}{sqrt(2/K)} rather than 1 for \eqn{K}{K}
            forests, so the prior on the combined location does not widen
            without bound as more forests are declared: with every forest
            carrying a basis it is the same at every \eqn{K}{K}, and
            with a basis-free first forest the basis forests' own
            contribution rises toward that value instead of past it.
            \code{forest(sd = )}
            overrides both and keeps its per-forest reading at every
            \eqn{K}{K} (two forests is the fixed point of the basis default,
            so nothing at two forests moves). \code{$setForestBasis} is the
            sole route by which a basis changes after creation, at any
            forest and any width, mirroring the matrix into
            \code{data@bases} as \code{$setWeights} does for weights; the
            amplitudes are preserved and remapped, so a width-preserving
            install leaves every one of them unchanged and a width change
            carries each forest's block to its new offset.
            \code{$setForestWeights} installs a per-forest, per-observation
            weight - a multiplicative precision factor on the named forest's
            own leaf conditionals, composing with \code{weights} and
            \code{active} as \eqn{(w_i a_i) m_f^2 s_i}{(w_i * a_i) * m_f^2 *
            s_i} rather than widening either channel, so it excludes a row
            from that forest's own leaf conditionals - and from that forest's
            empty-leaf veto - leaving its occupancy, its place in the
            combination, and the residual degrees of freedom untouched.
            \code{$getForestFits},
            \code{$getForestAmplitudes}, and \code{$getForestVariableCounts}
            read the forests' internal-scale fitted values, amplitudes, and
            split counts separately, each taking an optional \code{forest};
            the amplitude vector is ragged across forests, which is why one
            can be named. A run reports every surface for
            every draw: \code{$run} gains a \code{forestFits} element
            holding each forest's fitted values (n.obs x n.forests x
            n.samples) and a \code{glue} element holding the amplitudes
            that draw combines them through, stacked forest-major (sum(q) x
            n.samples); the \code{varcount} element similarly widens with a
            forest axis, forest-major within a draw and the prognostic
            forest first. Neither \code{forestFits} nor \code{glue} is
            present for any other model; a single-forest sampler's arrays
            are unchanged, and so is the flat C API's, whose
            \code{dbarts_results} declares no forest count and receives the
            prognostic forest's counts.
            \code{dbartsData} gains a \code{bases} argument naming the same
            matrices directly (a single basis, from either route, is refused
            when the sampler is built - a lone forest has no second ensemble
            for its amplitudes to distinguish it from), and predicting at new
            rows works both per forest and combined:
            \code{predict(fit, newdata, type = "forest", forest = )} replays
            each forest's saved trees and reports its own
            response-scale total - the out-of-sample twin of
            \code{extract(type = "forest")}, with the sampler's
            \code{$predictForests} serving the same read on the forests'
            internal scale - with no offset, response shift, or amplitude
            glue folded in, those being the caller's to add, given a
            \code{basis} of their choosing (\code{forest = NULL}, the
            default, returns every forest). Meanwhile
            \code{predict(fit, newdata, type =
            "ev"/"ppd"/"bart")} recombines that per-forest replay with
            each draw's own amplitudes and the response shift,
            \code{shift + sum_k (basis_k \%*\% glue_k) * forest_k}, taking
            the bases at those rows through a new \code{bases} argument (a
            length-K list, or the bare value when exactly one forest
            carries a basis) or re-deriving them from \code{newdata} when
            the fit declared them as a \code{bart2} \code{forest()} term;
            an \code{offset} at the new rows enters before the family's
            link on the combined arm only. Both amplitude arms require
            \code{keeptrees}/\code{keepTrees} (without the tree store the
            replay reports the current trees for every draw, which cannot
            be paired with the amplitudes draw by draw). A multi-forest
            sampler takes no test data at creation and no test surface
            afterward, its combined test fits having no off-sample basis; the
            per-forest arm separately refuses an offset, which shifts the
            recombination rather than any one forest's own total. Finally,
            \code{$setControl} carries a sampler's model-configuration
            attributes forward, so a multi-forest structure or a variance
            forest survives a control swap.
      \item Interaction and block constraints restrict which predictors may
            jointly shape the fit. The new \code{interactions} argument of
            \code{dbarts} and \code{bart2}, built with \code{\link{interactions}},
            offers three combinable rules: \code{interactions(max.order = K)} caps
            the distinct predictors on any root-to-leaf path (\code{K = 1} is
            additive), \code{interactions(forbid = list(c("a", "b")))} bars named
            predictors from co-occurring, and \code{interactions(groups = )}
            confines interactions to declared groups. The new \code{blocks}
            argument, built with \code{\link{blocks}}, is stricter still:
            it confines each whole tree to one declared group of predictors,
            so the ensemble is exactly a sum of per-group functions
            \eqn{f = \sum_G f_G}{f = sum_G f_G} (a functional-ANOVA /
            grouped-GAMI decomposition) over a partition that must be total
            and disjoint; \code{trees.per.group} optionally fixes each
            group's tree count (an even split by default).
            \code{interactions(groups = )} over the same kind of partition
            gives the block-additive property adaptively instead, with each
            tree's group chosen by the sampler rather than fixed. Both
            constraints are hard and apply per forest, so a multi-forest model
            can hold one forest additive-or-low-order while another stays
            free; on a single-forest model they compose with \code{monotone},
            which a multi-forest model refuses. Fits made without either
            argument are unchanged. See
            \code{?interactions} and \code{?blocks}.
      \item Monotone constraints: \code{monotone = c(x1 = "+", x3 = "-")}
            in \code{dbarts} and \code{bart2} constrains the fit to be
            monotone increasing or decreasing in named predictors
            (monotone BART; Chipman, George, McCulloch, and Shively
            2022). Directions accept the sign glyphs, the words
            \code{"increasing"}/\code{"decreasing"}, or \eqn{+1}/\eqn{-1};
            an unnamed length-p vector assigns them positionally. Only
            numeric and ordered predictors are eligible (a direction on
            an unordered factor is an error). A constraint forces
            birth/death-only proposals and a fixed \code{k} - 2 unless one
            is given, an explicit \code{k} hyperprior being an error; linear
            and Gaussian-process leaves are not supported under it.
            Unconstrained fits are unchanged. See \code{?dbarts}.
      \item Heteroscedastic BART: \code{variance = ~ x1 + x2} (or
            \code{variance = TRUE}) in \code{dbarts} and \code{bart2}
            adds a second forest modeling the residual variance surface
            \eqn{s^2(x)} as a product of scaled-inverse-chi-squared
            leaves (HBART; Pratola, Chipman, George, and McCulloch 2020).
            The shorthand takes 40 trees under the mean forest's own
            tree-structure prior; a \code{\link{varianceForest}(vars = ,
            n.trees = , base = , power = )} object in its place
            overrides any of them, its \code{vars} resolved through the
            same column selector the shorthand uses. The fit gains
            \code{s.train}/\code{s.test} draws of \eqn{s(x)}, and
            \code{predict} attaches \eqn{s(x)} at new data as an
            \code{"s"} attribute, which requires \code{keepTrees}. Gaussian
            responses only.
      \item Two new node priors fit smooth, rather than constant, leaf
            functions, through the \code{node.prior} argument of
            \code{dbarts}, \code{bart2}, \code{xbart}, and
            \code{dbartsSpec}. \code{node.prior = gp(columns, k,
            lengthscale, max.leaf.size)} fits each bottom node with a
            smooth Gaussian-process function of the designated continuous
            predictor columns under a squared-exponential kernel, the
            prior scale tied to \code{k} exactly as the other leaf priors;
            \code{k} may be fixed or given the \code{chi} hyperprior (the
            binary default), sampled from the drawn functions'
            standardized magnitudes; lengthscales fix at creation, per
            column or by a median-distance heuristic; leaves larger than
            \code{max.leaf.size} fall back to constant fits, confining the
            cubic kernel cost to small leaves; \code{getTrees} reports
            \code{NA} leaf values and no coefficient columns for these
            function-valued fits, and \code{keepTrees} storage grows with
            leaf sizes. \code{node.prior = linear(columns, k)} instead
            fits each bottom node with an intercept plus a linear term in
            the designated continuous predictor columns (standardized
            internally, every coefficient sharing the \code{normal(k)}
            prior), so a forest models smoothly-varying coefficients -
            treed regression in the sense of Chipman, George, and
            McCulloch, or the varying-slope structure of Bayesian Causal
            Forests when the covariate is a treatment indicator;
            \code{getTrees} reports leaf intercepts in \code{value} plus
            one \code{beta.<column>} slope column per covariate, and
            \code{plotTree} labels linear leaves with their coefficients.
            Both compose with the full sampler surface - weights, probit
            and logistic responses, \code{keepTrees} with \code{predict}
            (saved replays reproduce recorded test fits exactly), state
            serialization, and the predictor-mutation methods - and
            \code{xbart} takes them through its own \code{node.prior},
            sweeping its \code{k} grid under them; for linear leaves,
            folds standardize the leaf
            covariates with full-data constants, matching the final fit's
            calibration the way folds bin on the full-data cut grid.
      \item \code{tree.prior = dart(...)} runs DART (Linero 2018), a Dirichlet
            prior over the split-variable probabilities inducing variable
            selection, with optional sampling of the concentration.
            \code{bart2}, \code{rbart_vi}, and \code{xbart} expose it through a
            \code{dart} argument taking \code{TRUE} or a full prior specification,
            along with \code{split.probs} for fixed, non-uniform variable
            weights. Runs return the per-sample probabilities alongside
            the variable counts, as \code{varprobs} in the sampler's results
            and on \code{bart2} and \code{rbart_vi} fits.
      \item Two ways to start a forest other than a fresh draw from the
            tree prior. A warm start (\code{bart2}'s \code{warm.start}, or a
            sampler's \code{installTrees}) accepts a donor fit on a
            different cut grid (differing in \code{n.cuts},
            \code{useQuantiles}, or its predictor values): the donor's
            splits are remapped onto the destination's grid, collapsing
            any the grid starves, the same way replacing a sampler's data
            remaps its existing splits; a donor with a different tree
            count or DART setting is refused. \code{n.grow.sweeps}
            (\code{bart2}) and the sampler's \code{growFromRoot} method
            instead build the initial forest by XBART-style grow-from-root
            (He, Yalov and Hahn 2019): each requested sweep rebuilds every
            tree from the root against the current residual, in place of a
            draw from the tree prior, before the ordinary
            Metropolis-Hastings sweeps take over. It is opt-in (the
            default \code{n.grow.sweeps = 0L} leaves fits unchanged),
            mutually exclusive with \code{warm.start}, and leaves the
            targeted posterior unchanged - only where the chain starts.
            The initializer places real categorical split rules rather
            than skipping categorical predictors: at ten or fewer categories
            present at a node the node's choice is exact over all category
            subsets, and above ten a sorted-prefix family carries the same
            total prior mass, so a variable's realized split mass stays
            continuous in its level count. Where a node holds missing rows of an
            ordinal predictor it scores the two missing directions of each
            cut separately, each with those rows in the child it names; where
            it holds none the direction moves no statistic and is drawn after
            the fact. Either way the pair carries the prior mass of both
            rules, so such a column enters the node's choice at its full
            weight relative to not splitting. With \code{dart = TRUE},
            categorical predictors accordingly receive real split counts
            during the grow sweeps, so factor
            signals inform the sampled split-variable probabilities from
            the first regular sweep.
      \item Factor columns in data frames now default to
            \code{factors = "categorical"} in \code{dbarts}, \code{dbartsData}, \code{bart2},
            \code{xbart}, and \code{rbart_vi} (a new argument on all five):
            each unordered factor stays a single predictor with subset
            splits instead of expanding into indicator columns. Pass
            \code{factors = "indicators"} for the previous behavior. Under
            \code{factors = "categorical"}, splits send a subset of a
            factor's levels down each branch (at most 65535, see below);
            ordered factors are coded as ordinal. Level tables are
            retained so test data are coded identically, with unseen
            levels an error. \code{getTrees} gains, when the sampler has
            categorical predictors, a \code{directions} column decoding
            each categorical rule into one "L"/"R" character per level;
            such a rule's \code{value} is \code{NA}, as it has no scalar
            cut point. \code{plotTree} labels categorical rules with the
            levels sent down the left branch.
      \item Missing values in the predictors are supported. Under
            \code{missing = "incorporate"}, the default on \code{dbarts},
            \code{dbartsData}, \code{bart2}, \code{xbart}, and \code{rbart_vi}, every split
            rule learns a direction for missing values of its variable,
            so an \code{NA} observation follows that rule's chosen branch
            (Missingness Incorporated in Attributes; Twala et al. 2008).
            The route is learned per column, and only where that column's
            training values were missing: an \code{NA} in test predictors or
            \code{newdata} is refused, naming the column, wherever that
            column was complete in training. \code{missing = "error"} rejects
            incomplete predictors. The response, weights, and offset must
            be complete, and predictor columns cannot be entirely
            missing; both interfaces previously dropped incomplete rows
            silently.
            \code{getTrees} gains, where a predictor carried training
            \code{NA}s, a \code{missing}
            column reporting each rule's missing route, and \code{plotTree}
            appends it to rule labels. \code{bart}
            refuses missing data outright, where it previously dropped the
            incomplete rows silently. On complete data the new default draws
            exactly as \code{missing = "error"} does.
      \item Sparse predictor matrices: the x/y interfaces of
            \code{dbarts}, \code{dbartsData}, \code{bart}, \code{bart2},
            and \code{xbart} accept a \code{Matrix::dgCMatrix}. Columns
            are ordinal only; per column, the engine stores codes in a
            rank-bitmap layout when at most 20 percent of rows are
            nonzero (roughly a seventh to a tenth of dense memory at real
            sparsity)
            and densifies otherwise, so the sparse object is never
            expanded. Missing values enter as stored \code{NA} entries.
            A \code{dgCMatrix} test set stays resident too, coded against the
            training cuts rather than densified, though it carries no factor
            levels and so cannot supply a categorical training column's
            values. A \code{dgCMatrix} design accepts
            between-sweep column mutation: \code{sampler$setPredictor(x,
            column = j)} (and the flat C API's \code{setPredictor} /
            \code{updatePredictor}) re-quantize a sparse column in place
            when its nonzero pattern is unchanged and rebuild its layout
            when the set of nonzeros changes, with the transactional
            roll-back a dense column has. Whole-matrix replacement
            (\code{sampler$setPredictor(x)}) is accepted as well, and
            leaves the design a \code{dgCMatrix}; a replaced column's
            storage densifies permanently, since every row then differs
            from the implicit zero. Per-observation replacement of a
            sparse-backed column and \code{setData} stay fixed
            at creation. State serialization and \code{predict} work as
            usual. \code{rbart_vi} still requires dense inputs, and
            \code{node.prior = linear()} does not designate a sparse-backed
            column. A data frame passed to the x/y
            interfaces may also hold \code{Matrix::sparseVector} and
            \code{dgCMatrix} columns (assign them into the frame; they
            do not survive \code{data.frame(...)} or \code{I()})
            alongside ordinary numeric and factor columns, mixing dense
            and sparse predictors: sparse columns pass through unexpanded
            and tier exactly as an all-sparse design does, while
            dense-backed columns keep the full feature set - factors
            still split categorically and \code{node.prior = linear()}
            may designate dense-backed columns. A mixed design accepts
            between-sweep column mutation over dense-backed and
            sparse-backed columns alike, including a sparse
            \code{sparseFactor} column, and one call may name both kinds,
            with per-observation updates (\code{forceUpdate = "partial"})
            open on dense-backed columns. Whole-matrix replacement leaves
            \code{data@x} a container; a replaced sparse column's storage
            densifies permanently, as for an all-sparse design. A test
            data frame's sparse columns ride unexpanded as well, its dense
            columns materializing as usual, and a predictor
            container built from other data is re-coded against the
            training factor levels before it is used, so its codes mean
            what the training design's codes mean. The sparse column
            constructor is exported as \code{\link{sparseFactor}}; the
            container a mixed frame assembles into is documented as
            \code{\link{dbartsMixedMatrix}}, with its own \code{dim},
            \code{dimnames}, \code{[}, and \code{as.matrix} methods.
      \item An ordered factor predictor splits between every adjacent pair of
            its levels. Such a column is coded ordinal on its level codes,
            and its split grid is the midpoints between consecutive
            \emph{declared} levels - one candidate threshold per level
            boundary - so every pair of neighboring levels can be told apart,
            each distinct threshold enters the split prior with the same
            weight as any other, and a level that appears in test data but
            not in training still has its own boundary. \code{n.cuts} does
            not apply to such a column, exactly as it does not to an
            unordered factor: a factor's grid follows its level table, and a
            column wanting a fixed number of cut points over a range is a
            numeric one. An ordered factor is limited to 65534 levels, one
            below the unordered limit below, and a wider one is an error at
            every entrance that can declare a level table: naming the column
            for a factor, and naming the kind and its limit for a
            \code{sparseFactor} or a \code{varTypes} declaration.
      \item Creating a sampler over a data frame no longer assembles its
            factor columns into a temporary block of doubles: their integer
            level codes cross to the engine as the integers they already
            are, so each factor column's temporary is 4 bytes a row rather
            than 8. This covers a frame carrying a sparse column and a test
            frame as well as one of ordinary columns; a design supplied as a
            numeric matrix is unaffected. Nothing about the fit changes; the
            same codes reach the same store.
      \item The engine no longer keeps a factor column's values as doubles
            at all. Such a column's cells are level codes, which the
            quantized codes the engine already holds carry, and its cut grid
            follows its level table rather than its values, so the doubles
            duplicated information and served nothing. Two designs retained
            them: a test set, and a design supplied as a data frame carrying
            a sparse column. Both now hold 8 bytes a cell less per factor
            column for the life of the sampler - measured as peak resident
            memory, 160 MB on a test set of a million rows with twenty
            factor columns, and 233 MB on a data frame of a million rows
            with twenty factor columns beside a sparse one, where the
            temporary above accounts for the rest. Results are identical bit
            for bit, and a factor column designated as a linear- or
            GP-leaf covariate keeps its own copy, so those fits are
            unaffected too.
      \item Predicting over a resident test set - one the engine routes in
            place rather than densifying, which is what a data frame
            carrying a \code{sparseFactor} column gives - reads its factor
            columns as the integer level codes they already are instead of
            widening them first. The temporary those columns occupy halves:
            a 4-byte level code in place of an 8-byte double for every
            factor cell of the test set. The predictions are identical bit
            for bit, and the saved-tree replay behind \code{getTrees} reads
            the same way.
      \item Categorical predictors support up to 65535 levels, where a
            0.9-x categorical rule held its directions in a single 32-bit
            word: a rule's mask lives inline up to 63 levels and is pooled
            per tree above that. In \code{getTrees} such a rule's
            \code{value} is \code{NA}, as for any categorical rule, and the
            \code{directions} column carries the decode as before. Everything
            else composes unchanged, including missing values, saved
            trees with \code{predict}, and state serialization.
      \item Prior specifications are ordinary objects built by the
            constructors in the new exported list \code{dbartsPriors} (\code{cgm},
            \code{dart}, \code{normal}, \code{linear}, \code{gp},
            \code{chi}, \code{chisq}, \code{fixed}), usable
            programmatically and passable to the
            \code{tree.prior}/\code{node.prior}/\code{resid.prior} arguments;
            \code{dart}, \code{linear}, and \code{gp} are themselves new.
            The bare-name forms (\code{node.prior = normal(chi(1.5))}) keep
            working and, as in 0.9-x, resolve against dbarts' own vocabulary
            inside those arguments, so attaching packages that export
            \code{normal} or similar cannot change their meaning.
      \item The binary default for \code{k} - used by \code{dbarts},
            \code{bart2}, and \code{rbart_vi} when \code{k} is left
            unspecified - is now
            the proper \code{chi(1.5, 2)} rather than the improper 0.9-x
            default, spelled \code{chi(1.25, Inf)} there and, under the
            corrected degrees-of-freedom count listed under BUG FIXES, the
            same law as \code{chi(1.5, Inf)} now. It centers the sampled
            \code{k} near the field-standard fixed value of 2 (prior median
            1.9) while still adapting to the data. \code{chi()}'s own
            defaults move with it, from \code{chi(1.25, Inf)} to
            \code{chi(1.5, 2)}, so a bare
            \code{chi(df)} call is proper unless \code{Inf} is passed
            explicitly. Pass \code{k = chi(1.5, Inf)} for the
            old prior; seeded default binary fits otherwise draw
            differently.
      \item A model's leaf-prior calibration can be named at creation, in
            response units, rather than inherited from the range of the
            vector the sampler happened to be constructed on:
            \code{prior.scale} (the prior standard deviation of the forest
            total at \code{k = 1}) as a \code{dbartsModel} slot, as
            \code{node.prior = normal(scale = )} on \code{dbarts} and
            \code{xbart} (also \code{sd = }, the same quantity at the
            resolved \code{k}, refused under a \code{k} hyperprior), and as
            a \code{prior.scale} argument of \code{bart}, \code{bart2}, and
            \code{rbart_vi}; it is re-derived at every \code{setModel}, so a
            round trip through the model object holds it, and unset,
            nothing changes. A sampler reports the calibration actually IN
            FORCE, per chain and per forest, with \code{$getCalibration}:
            \code{prior.scale}, \code{prior.sd}, \code{prior.mean},
            \code{k}, \code{k.has.hyperprior}, the response transform, and
            a \code{leaf.model} tag. It reads the engine rather than the
            model object, so it is the first surface to show a
            re-anchoring channel - a \code{setResponse} or
            \code{setOffset} at \code{updateScale = TRUE}, or a
            \code{setData} - moving what is in force while the named intent
            on the model stays where it was put. For the monotone, linear, and
            Gaussian-process leaf models the reported \code{prior.sd} is a
            stated bound rather than an equality. \code{$setCalibration}
            rewrites its scale half on every chain, mid-chain, touching
            nothing else - unlike \code{setModel}, which also re-pins
            \code{sigma} and refuses a DART sampler outright; its
            \code{prior.sd} spelling meets the same sampled-\code{k}
            refusal as at creation. On a multi-forest model
            those five columns carry the calibration map the combination
            runs under, which no other R surface
            shows: the forest's amplitude prior variance or half-Cauchy
            scale (\code{amplitude.prior.variance} and
            \code{amplitude.prior.scale}, exclusive per forest), its
            \code{node.scale.factor} and \code{node.scale.divisor}, and the
            \code{basis.row.norm} the map divides out - so the induced
            prior on a multi-forest model's combined location can be
            reconstructed from the sampler rather than from the call,
            including after a \code{$setForestBasis} swap re-derives the
            row norm (\code{NaN} on a forest whose scale the map does not
            own, as on every single-forest sampler).
            \code{dbarts_sampler_getForestCalibration} reports the same
            quantities to a C consumer through a
            \code{dbarts_forest_calibration} struct that carries its own
            \code{structSize} and a documented 1.0-0 field boundary, below
            which every later append goes, so a caller compiled against the
            1.0-0 header keeps reading as before.
      \item \code{\link{samplePriorPredictive}} draws repeatedly from a
            sampler's prior - tree structure, end-node parameters, and
            optionally response noise - and evaluates the resulting forest
            at \code{x.test}, for inspecting a prior before ever fitting to
            data (e.g. the implied prior on a treatment effect
            \eqn{f(x_1) - f(x_0)}{f(x1) - f(x0)}). It runs on a freshly
            constructed private sampler built from the caller's control,
            model, and data, so the sampler passed in is never modified and
            successive calls draw independently unless the sampler it is
            given was created with \code{dbartsControl(seed = )}.
            \code{type = "ev"} draws the forest sum (or the
            link-transformed probability); \code{type = "ppd"} adds the
            observation model on top, sigma freshly drawn from its own
            prior each sample for a gaussian sampler. A heteroscedastic
            sampler refuses \code{type = "ppd"} rather than returning a
            homoscedastic prior predictive.
      \item Several small fitting-surface additions. \code{fitted} and
            \code{predict} gain a \code{ci.level} argument: given a level
            in (0, 1) they report \code{est}, \code{ci.lower},
            and \code{ci.upper} per observation - a matrix, or a trailing
            margin of three on a per-category \code{type} - in place of
            \code{fitted}'s posterior mean and \code{predict}'s draws,
            the interval kind following \code{type} - \code{"ev"} a
            credible interval for the mean (a probability for binary
            responses), \code{"ppd"} a wider prediction interval carrying
            the residual noise, and \code{"bart"} a credible interval on
            the latent scale. \code{dbarts} gains a \code{seed} argument
            mirroring \code{dbartsControl(seed = )}; \code{residuals}
            gains a \code{type} argument that reaches \code{fitted} on
            \code{bart}, \code{rbart}, and \code{bartHurdle} fits, so
            link-scale residuals are available (the multinomial, ordinal, and
            count fits refuse it by name); and \code{plotTree} is now
            a generic that plots a tree from any kept fit - \code{bart} and
            \code{rbart}, and the multinomial, ordinal, count, and hurdle
            classes - without the caller reaching into the sampler
            stored on it. \code{bart2} and \code{rbart_vi} default to
            \code{combineChains = TRUE}: a multi-chain fit's draw channels
            arrive with the chain and sample margins already merged,
            matching what \code{summary}, \code{as_draws_df}, and most
            downstream consumers expect without a reshape;
            \code{combineChains = FALSE} keeps them separate, as the
            previous default did. Finally, MCMC runs are now
            interruptible: pressing \code{Ctrl-C} during a \code{bart},
            \code{bart2}, \code{rbart_vi}, \code{xbart}, or sampler
            \code{run} stops it between iterations - joining any worker
            threads first - and returns to the prompt with an error,
            instead of deferring the interrupt until the whole run
            finishes; the check is polled about ten times a second and
            does not affect the draws.
      \item \code{\link{xbart}} is reimplemented as an R-level driver over
            the sampler interface, so cross-validation reaches the sampler's
            own features rather than a separate engine path: \code{family},
            \code{factors}, \code{missing}, \code{dart},
            \code{node.prior}, \code{tree.prior}, \code{resid.prior}, and
            \code{storage} are all arguments of it, and it takes a sparse
            (\code{dgCMatrix}) design. Chains warm-start only
            across parameter settings that share a tree count over one data
            split, never across tree counts, folds, or replications: the
            retired implementation carried chains between
            splits whose held-out rows the previous training set contained,
            letting a slowly mixing setting - a small ensemble especially -
            score against data it had effectively seen. Parallelism is a
            \pkg{parallel} cluster over replications, and results are
            reproducible for a fixed \code{seed} and \code{n.threads}.
      \item \code{predict}'s \code{n.threads} parallelizes for real, on every
            predict method (\code{bart}, \code{bartMultinomial},
            \code{bartOrdinal}, \code{bartNegbin}, \code{bartHurdle},
            \code{rbart}) and on the sampler's own \code{predict} and
            \code{predictForests}. 0.9-x partitioned a replay across chains
            only; this partitions the saved-tree replay per (chain,
            posterior draw), so a single-chain fit parallelizes too. Each
            partition writes its own rows and nothing is reduced across
            workers, so the answer is identical bit for bit at every value;
            a replay too small to pay for the partition runs inline. The
            argument defaults to the fit's own thread count, does not
            persist, and is the LAST positional formal on every method. A value
            that is not a single positive integer is an error. \code{run}'s
            \code{n.threads} argument remains inert.
      \item Reporting and diagnostics grew two surfaces with no 0.9-x
            counterpart. Every family's pointwise log-likelihood is now
            available for model comparison: R's
            \code{\link{extract}(type = "loglik")}, and the flat C API's
            \code{logLikelihood} results field for every model whose
            location the response model can see (it is \code{NaN} wherever
            the combined location is not visible to the response model to
            score - any amplitude-combined sampler, at any forest count, and
            the multinomial softmax), score each training observation at
            each posterior draw under the family's own likelihood - the
            Bernoulli mass for a binary family, the marginal
            \eqn{t_\nu}{t_nu} density for
            \code{resid.dist = student()}, the multinomial density of the
            observed count row, and so on, each documented under its
            family's own fit method; when chains are combined the result is
            a samples-by-observations matrix directly consumable by
            WAIC/PSIS-LOO implementations such as those in the \pkg{loo}
            package (a per-chain array, \code{combineChains = FALSE}, is
            reordered to the draws-by-chains-by-observations shape
            \code{loo::relative_eff} expects with
            \code{aperm(x, c(2, 1, 3))}). Separately,
            \code{\link{summary.bart}} (and \code{summary.rbart}) reports a
            per-variable posterior summary of a fit's scalar parameters
            (\code{sigma}, \code{k}, \code{tau} for \code{\link{rbart_vi}}),
            with split-\eqn{\hat{R}}{Rhat} and effective sample size when the
            \pkg{posterior} package is installed, and a built-in fallback
            summary when it is not; \code{as_draws_array} and
            \code{as_draws_df} methods for \code{bart}/\code{rbart} fits
            convert a fit's chain-dimensioned draws to \pkg{posterior}'s
            array/data frame conventions. Those two are methods for
            \pkg{posterior}'s own generics and are registered only when it
            loads, \pkg{posterior} being a \code{Suggests} rather than a
            hard dependency; \code{summary} is always available.
      \item A flat C interface for packages compiling against dbarts,
            \code{dbarts/dbarts.h}: an opaque sampler handle over the rewritten
            engine, created from the R specification objects, with entry
            points for running into caller-owned buffers, the mutation
            surface (response, offset, weights, sigma, predictors, test
            data), latents, prediction, tree storage and introspection,
            and state save and restore. Entry points are reached through
            \code{R_GetCCallable} and versioned by \code{DBARTS_C_API_MAJOR}
            and \code{DBARTS_C_API_MINOR}, with \code{dbarts_apiHash()} as an
            opt-in exact-ABI check, which moves on additive releases as well.
      \item New \code{dbartsSpec} resolves the \code{control}/\code{model}/\code{data} triple and
            family token that describe a sampler, without constructing one,
            for packages that embed dbarts and supply their own design
            matrix. It shares its resolution with \code{dbarts}, so a sampler
            built from its output draws identically to the equivalent
            \code{dbarts} call given the same data object (the two can still
            resolve a family differently, \code{dbarts} ingesting a raw
            response this function requires already shaped), and its returned
            family token names the \code{dbarts_family} enumerator
            \code{dbarts_sampler_create} takes. Together with \code{dbartsData} it
            makes the whole single-forest feature set - response families,
            monotonicity, interaction and block constraints, a variance
            forest, robust residuals - reachable from a \code{LinkingTo: dbarts}
            package without reaching into an unexported function.
      \item A sampler accepts a per-observation 0/1 active-row mask
            (\code{$setActiveRows}), which removes a row from the likelihood
            and from the latent draws of any sampler - gaussian, Student-t,
            \code{probit}, \code{ordinal}, \code{logistic}, \code{nbinom},
            \code{aft}, and multinomial alike - between draws while keeping
            its fitted value and its
            leaf occupancy - the row stays in the forest, just not in the
            model this sweep. \code{NULL} clears; an all-ones mask is
            treated as no mask; fractional or \code{NA} values are refused.
            Independent of \code{setWeights} and composes with it.
      \item A sampler carrying grouped random effects - the one
            \code{\link{rbart_vi}} returns in \code{$fit} when its prior is
            built in - is mutable on
            the response side: \code{$setResponse} and \code{$setOffset}
            accept a same-length replacement at the pinned scale
            (\code{updateScale = FALSE}, the default), so a grouped BART fits
            inside a larger Gibbs sampler on the same terms as an ungrouped
            one. The group effects and their scale are carried across the swap
            as the Gibbs block they are, and the per-observation group indices
            are unchanged, so every row stays in its group.
            \code{updateScale = TRUE} is refused outright, rather than
            defaulted away from, on both \code{$setResponse} and
            \code{$setOffset} under a gaussian, Student-t, or \code{"aft"}
            response: the random intercepts and their scale are held against
            the transform fixed at creation and nothing converts them, so
            re-anchoring would silently restate both in response units. A
            grouped \code{"probit"} or \code{"logistic"} sampler accepts it,
            where it is a no-op.
            \code{$setData} is refused: it may change the number of rows,
            and the grouping is fixed at creation, so the new rows would have
            no group membership to read - re-create the sampler instead.
      \item A sampler reports its combined per-observation location on the
            response scale WITHOUT the installed offset
            (\code{$getFitsWithoutOffset()}), which is what an outer block
            conditioning on \eqn{f(x)}{f(x)} needs: \code{run()$train}
            carries the offset, so differencing it against
            \code{$getLatents()} is biased. On a Bayesian causal forest it is
            the only route from R to the combined fit at the sampler's
            current state - \code{$predict} is refused and
            \code{$getForestFits} reports the forests' internal scales, while
            \code{run()$train} reports the same quantity at recorded draws,
            with the offset folded in. Refused on a
            multinomial sampler, whose channels are per-category softmax
            probabilities and not one additive location; there is no
            with-offset twin, the caller adds back the offset it installed.
            \code{?dbartsSampler} gains a table of which method reports which
            scale, with or without the offset, and from the current state or
            from stored draws. Built on that same offset-free location, the
            augmentation draws a non-gaussian family runs inside each sweep
            are callable from R too:
            \code{\link{dbartsDrawLatents}} draws the latent at a supplied
            fit, and \code{\link[=dbartsDrawLatents]{dbartsWorkingResponse}}
            turns it into the quantity a host regresses on - which for
            \code{"logistic"}, \code{"nbinom"} and
            \code{resid.dist = student()} is not the latent at all, those
            three reporting precisions. The draws come from R's own stream, so
            \code{\link{set.seed}} reproduces them and they compose with any
            other R draw; a sampler's chain generator is untouched, and both
            take and return the location WITHOUT the offset -
            \code{dbartsDrawLatents}'s \code{fit} as
            \code{$getFitsWithoutOffset()} reports it, and the working
            response with the offset subtracted back out.
      \item Six named compositions ship as a vignette,
            \code{vignette("dbarts-as-a-component")}: data augmentation written
            in R, an outer block on the offset channel and its calibration
            check, a K-sampler decomposition sharing one prior budget, a latent
            predictor shared across samplers, an outer-owned sigma, and the
            incremental read. Each is a runnable outer sampler and each is also
            a test cell. \code{\link{dbarts-embedding}} is the reference page
            behind them: which channel an outer block writes, which fit it
            reads back, and whose prior is in force after a swap. Checking
            that a composition actually targets the posterior it claims is
            itself a new tool: \code{\link{dbartsValidateComposition}} runs
            simulation-based calibration over a one-sweep step the caller
            supplies, for a model the caller states as a prior draw and a
            simulation, and reports per scalar functional whether the
            composed kernel is well calibrated. It catches the defect class
            in which a posterior mean is reported where a draw belongs, and
            needs no sampler: the functionals may be derived, the ranks are
            tie-broken, and the uniformity band is corrected over the
            functionals ranked in the call.
    }
  }
  \subsection{C API}{
    \itemize{
      \item The C++ interface in \code{dbarts/R_C_interface.hpp} is removed with
            the pre-1.0 engine it exposed; packages compiling against
            dbarts use the flat C interface in \code{dbarts/dbarts.h} instead.
      \item \code{dbarts_sampler_run} now errors when a non-null \code{dbarts_results}
            has a zero \code{structSize} (the caller allocated the results struct
            but forgot to set its size) rather than silently filling no fields
            and returning an uninitialized buffer. A new \code{DBARTS_RESULTS_INIT}
            initializer in \code{dbarts/dbarts.h} sets \code{structSize} and zeroes the
            field pointers, so the size cannot be forgotten:
            \code{dbarts_results r = DBARTS_RESULTS_INIT;}.
      \item The C API reaches Bayesian causal forests:
            \code{dbarts_sampler_numForests} reports the model's forest count,
            \code{dbarts_sampler_setForestBasis} swaps a forest's basis, and
            \code{dbarts_sampler_getForestFits} and the amplitude pair
            (\code{dbarts_sampler_numForestAmplitudes},
            \code{dbarts_sampler_getForestAmplitudes}) read the per-forest fits
            and the amplitudes the two combine through.
            \code{dbarts_sampler_setResponse} takes an explicit
            \code{updateScale} argument matching
            \code{dbarts_sampler_setOffset}, in place of always re-deriving the
            response transform, and the mutation entries now carry the same
            per-model rules as their R counterparts: a two-forest sampler takes
            a scale-pinned response, offset or weight swap and refuses its whole
            test surface.
      \item The flat C API's predictor entries
            (\code{dbarts_sampler_setPredictor},
            \code{dbarts_sampler_updatePredictor},
            \code{dbarts_sampler_setTestPredictors} and
            \code{dbarts_sampler_predict}) take a single self-describing source
            struct, \code{dbarts_predictor_source}, so a C consumer hands the
            sampler compressed-column storage for prediction and test data
            without densifying it and without an R container class, and every
            predictor argument declares its own width and its own CSC column
            count.
      \item \code{dbarts_sampler_numTrees}, \code{dbarts_sampler_getTrees} and
            \code{dbarts_sampler_printTrees} take a forest index, so they are
            unambiguous on a two-forest sampler; \code{dbarts_apiHash()} gives
            a consumer an opt-in exact-ABI lockstep check alongside the
            major/minor handshake; and \code{dbarts_sampler_setForestWeights},
            \code{dbarts_sampler_setActiveRows},
            \code{dbarts_sampler_getForestCalibration} and
            \code{dbarts_sampler_setForestPriorScale} reach the per-forest
            weight, active-row mask and leaf-prior calibration channels from C.
      \item A C host reads the negative-binomial dispersion \eqn{r}{r} without
            an R round trip: \code{dbarts_results} gains a \code{dispersion}
            field shaped as \code{sigma}, and \code{dbarts_sampler_getDispersion}
            reads the value in force mid-sweep, returning 0 and touching
            nothing off a family carrying none. The field is appended, so a
            consumer compiled against the smaller struct is unaffected.
      \item \code{dbarts_drawLatents} and \code{dbarts_workingResponse} mirror
            \code{\link{dbartsDrawLatents}} and \code{dbartsWorkingResponse}
            over caller-supplied arrays: the same laws, the same family rules
            and response-support validation, drawn from R's own stream.
      \item \code{dbarts_sampler_predict} takes a \code{numThreads} argument
            before its output pointer: a per-call worker count that does not
            persist, where 0 means the sampler's own count and a resolved count
            below 1 - including a sampler whose own count was set to 0 - is
            treated as 1, so the output buffer is always written. The replay is
            identical bit for bit at every value.
      \item \code{DBARTS_C_API_HASH} is re-baked for the entry points above.
            All but \code{dbarts_sampler_predict} are append-only, so only that
            one's call sites change; every consumer checking
            \code{dbarts_apiHash()} against the header constant must be
            rebuilt - the lockstep handshake as designed.
      \item A \code{dbarts_family} enum replaces the stringly-typed
            \code{family} argument of \code{dbarts_sampler_create},
            \code{dbarts_drawLatents} and \code{dbarts_workingResponse}, and a
            new \code{dbarts_sampler_family} accessor reports the family a
            sampler was built with (or what \code{DBARTS_FAMILY_AUTO} resolved to).
            \code{dbarts_sampler_getForestFits},
            \code{dbarts_sampler_getForestAmplitudes},
            \code{dbarts_sampler_getDispersion} and
            \code{dbarts_sampler_getForestCalibration} are renamed from
            \code{dbarts_sampler_forestFits},
            \code{dbarts_sampler_forestAmplitudes},
            \code{dbarts_sampler_dispersion} and
            \code{dbarts_sampler_forestCalibration}, matching the \code{get}
            verb every other reader already carries.
            \code{dbarts_sampler_setVerbose}'s \code{printEvery} widens from
            \code{uint32_t} to \code{size_t}, and
            \code{dbarts_sampler_printTrees} gains a \code{useLiveTrees}
            argument matching \code{dbarts_sampler_getTrees}'s. The
            \code{DBARTS_USE_STUBS} handshake now enforces
            \code{dbarts_apiMajorVersion}/\code{dbarts_apiMinorVersion} by
            default on every stub's first resolution; the exact-ABI hash check
            is opt-in, behind \code{DBARTS_REQUIRE_EXACT_ABI} defined before
            including the header. \code{DBARTS_C_API_HASH} is re-baked again
            for this surface.
      \item Seven entry points that refused by raising now answer in a return
            value where the refusal is a fixed property of the sampler:
            \code{dbarts_sampler_setResponse}, \code{dbarts_sampler_setOffset},
            \code{dbarts_sampler_setWeights}, \code{dbarts_sampler_setSigma},
            \code{dbarts_sampler_setTestPredictors},
            \code{dbarts_sampler_setTestOffset} and
            \code{dbarts_sampler_predict} return \code{int} in place of
            \code{void}, giving 0 without touching the sampler where it cannot
            carry the channel at all - a probit or heteroscedastic sampler
            handed a sigma, a family with no case weights, a coupling whose
            test blend is undefined - so a host driving a sampler it did not
            build probes instead of unwinding through its own frames. What a
            different argument would have fixed still raises: a response
            outside its family's support, a logistic count that is not a
            positive integer, a scale update against a pinned transform, a
            malformed source. On the same rule
            \code{dbarts_sampler_setActiveRows} moves its non-binary mask
            refusal from 0 to an error, leaving its 0 to mean the family
            implements no mask. \code{dbarts_sampler_getTrees} and
            \code{dbarts_sampler_printTrees} take \code{forest} as the
            argument after the sampler, where every other forest-naming entry
            has it and where it precedes the tree indices it qualifies, and
            \code{dbarts_sampler_setForestBasis} renames its \code{basis}
            parameter to \code{basisRowMajor}, the header's one documented
            exception to column-major storage. The header no longer defines
            \code{USE_FC_LEN_T}: it declares no Fortran routine, and the define
            turned \code{FC_LEN_T} on for a consumer's whole translation unit
            through \code{Rconfig.h} with no way to restore it - a consumer
            calling BLAS or LAPACK defines the macro itself, as R's own
            documentation directs, and one that relied on this header to pull
            in \code{<Rversion.h>} now includes it itself.
            \code{DBARTS_C_API_HASH} is re-baked for these signatures.
      \item \code{dbarts_drawLatents}'s \code{cutpoints}/\code{numCutpoints}
            parameters are renamed \code{ordinalThresholds}/
            \code{numOrdinalThresholds}, disambiguating the ordinal family's
            latent cut points from the engine's own \code{cutPoints} (the
            tree split grid's per-column candidate list); the two were never
            the same quantity. \code{DBARTS_C_API_HASH} is re-baked for this
            signature.
      \item \code{dbarts_column_type} gains a third enumerator,
            \code{DBARTS_COLUMN_ORDERED_FACTOR = 2}, which a
            \code{dbarts_predictor_source} may declare for a column whose
            level codes are ordered: such a column splits at a threshold on
            those codes rather than by a subset mask, as an ordinal column
            does, and \code{dbarts_sampler_setPredictor} and its siblings now
            name it in the refusal they raise for an unknown column type. The
            two existing enumerators keep their values.
            \code{DBARTS_C_API_HASH} is re-baked for the appended enumerator;
            every consumer checking \code{dbarts_apiHash()} against the header
            constant must be rebuilt.
      \item \code{dbarts_predictor_source} gains a \code{denseCodes} block and
            the \code{numDenseCodeColumns} width that bounds it, so a caller
            whose factor columns are already integer level codes hands them
            over as \code{int32_t} instead of building the double block the
            entries used to require of it; \code{INT_MIN} (R's
            \code{NA_INTEGER}) marks a missing code. The mutation entries
            widen that block once internally, so there the field saves the
            caller the loop and the buffer rather than the conversion; the
            test-side entries read the codes where they lie and build no
            double block at all. Both channels are indexed by the same
            \code{columnSources[j]}, within the channel the sampler's kind
            for that column selects: the
            double index is bounded by \code{numColumns} and the code index
            by \code{numDenseCodeColumns}, and an entry refuses a source that
            breaks either bound or that leaves out a channel its columns
            read. A null \code{denseCodes} behaves exactly as every caller
            written before the field did. Both fields are appended below the
            struct's 1.0-0 boundary, so a consumer compiled against the
            smaller struct is unaffected, but \code{DBARTS_C_API_HASH} is
            re-baked for them; every consumer checking
            \code{dbarts_apiHash()} against the header constant must be
            rebuilt.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item The \code{\link{predict.glm}} synonyms documented for \code{type}
            work on the \code{bart2} extended-family generics: \code{extract},
            \code{fitted} and \code{predict} for \code{bartMultinomial},
            \code{bartOrdinal} and \code{bartNegbin} fits take
            \code{"response"} for \code{"ev"} and \code{"link"} for
            \code{"bart"}, as the \code{bart} methods always have. Previously
            each refused both, naming only the types it accepted.
      \item \code{type = "bart"} and \code{type = "forest"} on a multinomial
            fit error naming the reason - a category's forest is a latent
            whose level is reproducibly structured yet not identified, so the
            identified content is the log-ratio of the reported probabilities
            - on \code{predict} and \code{fitted} as well as \code{extract},
            which alone gave a reason before.
      \item A refused state install no longer leaves a sampler holding a
            live but unfitted engine. \code{getPointer}'s transparent
            re-creation after save and load, and \code{setState} on a
            sampler whose pointer has died, now bind the re-created engine
            only once the state has installed into it; previously a refusal
            left the empty engine in place, so the next \code{run} sampled
            stumps and \code{storeState} overwrote the stored state with
            them. A refused install leaves the object exactly as it was.
      \item Restoring a saved state into a \code{logistic}-family sampler
            carrying different case weights than the state was stored under no
            longer pairs the stored Polya-Gamma latents with the new counts;
            the latents are re-derived against the weights in force, as
            \code{setWeights} does. Gaussian, Student-t, and every
            weight-refusing family are unaffected.
      \item A leaf containing only observations of weight zero now counts as
            empty for the tree-move vetoes, so no tree can keep a leaf that
            no observation in the likelihood reaches. Previously such a leaf
            was treated as occupied and drew its value from the prior. Fits
            with no weights, or with strictly positive weights, draw exactly
            as before.
      \item Trees drawn from the prior - by \code{sampleTreesFromPrior},
            \code{samplePriorPredictive}, the C API's
            \code{dbarts_sampler_sampleTreesFromPrior}, and the initial
            forest every fit starts from - now realize the sampler's actual
            tree prior, the Chipman-George-McCulloch law conditioned on
            carrying no empty leaf, which is what the tree moves price.
            Previously they were an unrestricted draw with its empty leaves
            collapsed away afterward, a different distribution that also
            admitted leaves holding only zero-weight observations, out of
            which a birth scored \code{NaN} and was silently rejected.
            Seeded results move with the initial forest.
      \item That conditioning is now applied per forest, against the weights
            the forest's own tree moves veto with - the coupling's
            per-observation multipliers and any weight installed by
            \code{$setForestWeights} - rather than the sampler's shared
            weights. A Bayesian causal forest is the case that moves by
            default: its treatment forest enters with amplitudes
            \code{(b0, b1) = (0, 1)}, so a control row carries no weight
            there, and its initial trees held leaves that no treated row
            reached. \code{$growFromRoot} keeps the same law, its scan
            having admitted a split whose side carried members but no weight.
            Where that vector is entirely zero - an all-zero
            \code{$setActiveRows} mask, or an all-zero
            \code{$setForestWeights} - there is no distribution to draw
            from, and the forest now initializes to single-leaf trees with
            zero parameters, the one structure that stays legal whatever
            weights are restored afterward. Drawing from the prior in that
            state previously signalled an error, on a sampler that was
            otherwise accepted and could be run.
      \item Mutating a dense-backed column of a mixed dense/sparse predictor
            container now updates the values the sampler re-reads, not only
            its bins. Cut point installation (\code{setCutPoints}), a state
            replay that reinstalls cut points (\code{setState}), and
            \code{linear} and \code{gp} leaf covariates all see the mutated
            column; each previously read the values the container was built
            from, so the next requantization silently reverted the column's
            bins and a leaf kept regressing on the original covariate. A
            rejected transactional update leaves the values as they were.
            Dense designs and pure \code{dgCMatrix} designs were unaffected.
      \item Under \code{factors = "categorical"}, a categorical predictor's
            category count is now the declared level table rather than the
            largest code the training rows happen to carry. A level declared
            but unobserved in training keeps its own bin, so test or
            replacement data carrying it are accepted at every entrance -
            creation, \code{setTestPredictor}, \code{predict}, and
            \code{setPredictor} - instead of being rejected as an existing
            category code. Dense factors, \code{sparseFactor} columns, and
            dense factors inside a mixed dense/sparse frame now agree; only
            \code{sparseFactor} columns previously honored their declared
            table. Designs whose factors observe every level are unaffected
            and draw identically. Call \code{\link{droplevels}} before
            fitting to model only the observed levels.
      \item \code{setPredictor} on a mixed or \code{dgCMatrix} design now
            refuses a whole-matrix replacement argument whose shape
            disagrees with the design instead of reinterpreting it
            column-major; previously only \code{is.matrix} arguments (never
            a \code{Matrix} class) had their dimensions checked, so a
            transposed argument of the same total length was silently
            accepted.
      \item A \code{dbartsMixedMatrix}'s per-sparse-column reference
            metadata (\code{sparseReference}, \code{sparseCategoryCount})
            is now validated at creation to the same standard
            \code{setPredictor}/\code{updatePredictor} always enforced:
            a container whose fields are missing, non-integer, or not
            sized one entry per sparse column is refused by name
            (\dQuote{malformed mixed predictor container}) when the
            sampler is built, rather than only once a later mutation
            or a resolved categorical column reads it.
      \item \code{missing = "error"} on \code{setPredictor},
            \code{setTestPredictor}, and \code{setTestPredictorAndOffset}
            now inspects a mixed or sparse container's actual values rather
            than the container's underlying list, which read a
            single-sparse-column container's reference-level metadata as a
            missing value that was not there, and could miss a real missing
            value stored in a container with two or more sparse columns.
      \item \code{setPredictor} accepts \code{Matrix::dgCMatrix},
            \code{Matrix::sparseVector}, and mixed-container values against
            a sparse or mixed design without the caller densifying them
            first.
      \item A single column of a sparse or mixed \code{x.test} can now be
            replaced with \code{setTestPredictor(x, column = j)}; previously
            any column update on a resident sparse or mixed test container
            was refused outright.
      \item A sparse test column declaring a reference level for a training
            column that is not categorical is now refused, at both
            \code{predict} and \code{setTestPredictor}. The two previously
            read such a column differently - \code{predict} densified with
            the declared reference, \code{setTestPredictor} stayed resident
            with an implicit zero - and disagreed silently.
      \item \code{predict} and \code{getTrees(newdata =)} route a sparse or
            mixed test set through the trees off its own storage; neither
            materializes the full dense test matrix any more. A designated
            leaf covariate column must still arrive dense.
      \item Seeded parallel \code{rbart_vi} now restores R's stream after
            drawing its chain seeds; the restore assigned a local
            variable instead of the global state.
      \item The sampler's \code{getSumsOfSquaredResiduals} method now returns
            residual sums of squares on the original response scale as
            documented; it previously descaled by a single factor of the
            data range instead of its square. Binary responses report on
            the latent scale.
      \item A \code{weights} vector whose length matches neither the response
            nor one is now an error for matrix and sparse inputs; it was
            previously recycled silently. Scalars still broadcast.
      \item \code{rbart_vi} accepts a custom prior function passed by name;
            any symbol-valued \code{prior} used to enter the built-in lookup
            and fail.
      \item Prior arguments no longer evaluate against internal
            \code{control} and \code{data} bindings, which silently
            shadowed a caller's own variables of those names: in
            \code{tree.prior = cgm(2, 0.95, data)}, \code{data} now
            resolves to the caller's object rather than the sampler's
            data under construction. Only \code{num.vars} and
            \code{numvars} remain bound alongside the prior
            constructors.
      \item The \code{chi(degreesOfFreedom = nu)} hyperprior on \code{k} now samples
            a chi distribution with \code{nu} degrees of freedom as documented;
            its Gibbs update used a shape corresponding to \code{chi(2 * nu - 1)},
            a dropped Jacobian in the k-squared transform. The binary
            default is relabeled from \code{chi(1.25, Inf)} to \code{chi(1.5, Inf)},
            which draws identically under the corrected shape, so seeded
            default binary fits are unchanged; fits that set \code{k} to a
            non-default \code{chi()} now draw the degrees of freedom requested
            and their results change accordingly.
      \item The tree change move's acceptance now satisfies detailed
            balance. It previously omitted the proposal-density ratio, a
            defect inherited from the original CGM-lineage implementation,
            which biased the tree-structure posterior toward
            low-cardinality and descendant-constrained split variables
            whenever a design's effective cut counts were unequal. Variable
            selection is now correct for mixed continuous and categorical
            designs and under DART. All seeded results change, and
            variable-inclusion summaries in unequal-cut designs shift
            accordingly.
      \item The residual-variance (\code{sigma}) posterior no longer counts
            zero-weight observations in its degrees of freedom. Such rows
            are documented as ignored and contribute nothing to the
            weighted residual sum of squares or any leaf, but the
            posterior degrees of freedom still counted them, which
            over-concentrated \code{sigma} and biased it low - severely when
            zero-weight rows were numerous. Fits that pass \code{weights} of 0
            change accordingly; all-positive-weight fits, including the
            unweighted default, are unaffected.
      \item \code{rbart_vi}'s \code{group.by}/\code{group.by.test} are now looked up by name
            instead of a positional scan that silently grabbed the wrong
            column (typically the response) when the given symbol was absent;
            \code{predict} also coerces \code{group.by} to a factor exactly as at fit
            time, so it now recognizes trained groups instead of silently
            re-drawing their random intercepts from the prior.
      \item A factor (or character) response no longer hard-errors
            ("'range' not meaningful for factors") when passed through a
            formula to \code{dbarts}, \code{bart2}, or \code{rbart_vi}; the formula path now
            routes it the same way the x/y interface always has.
      \item Several inputs that previously failed with an unrelated low-level
            error now raise a named one instead: a non-finite response
            (previously "missing value where TRUE/FALSE needed"), zero-row
            data ("subscript out of bounds"), and a non-positive \code{n.samples}
            ("dim(X) must have a positive length"). A \code{newdata} variable
            missing from the fit's predictors and sharing its name with a
            builtin function is now named in the error instead of reported as
            "invalid type (builtin)". An unnamed-matrix \code{newdata} with
            swapped columns now warns naming the assumed column mapping.
      \item A warm start (\code{$installTrees}) or \code{$setState} into a
            column-restricted forest (a column-restricted variance forest or
            a \code{blocks}-constrained forest) is now
            refused if the donor holds a tree that splits outside the
            destination's allowed columns, instead of installing it
            silently; such a tree previously scored its splits against an
            availability menu that excluded the very column it used,
            corrupting the tree-prior density with no visible symptom. The
            refusal happens
            before any live state is touched, so the destination sampler is
            unaffected and usable afterward. The two entries run one
            predicate over every forest the state carries - each mean
            forest, a \code{blocks} row, and the variance forest - and
            report the same message, so neither admits a donor the other
            refuses; \code{$setState} used to hold only the mean forests to
            the rule, and installed forbidden variance splits silently.
      \item \code{extract(type = "loglik")} and the posterior predictive
            draws of \code{extract(type = "ppd")} and \code{predict(type =
            "ppd")} evaluate a heteroscedastic fit (\code{variance = ~ x})
            at its own per-observation \eqn{s(x)}, so
            \eqn{y_i \sim N(f(x_i), s^2(x_i) / w_i)} throughout. All three
            used the scalar \code{sigma}, which under the variance forest's
            parameterization is a fixed unit residual times the range of the
            response - a constant carrying no posterior content - so the
            reported log-likelihood was wrong by a factor that grew with the
            spread of the variance surface, and posterior predictive draws
            added homoscedastic noise at the wrong level. Log-likelihoods
            feeding \code{loo}/\code{waic} and posterior predictive checks on
            such a fit move accordingly; homoscedastic fits are unchanged,
            draw for draw. \code{predict(type = "ppd")} already carried
            \eqn{s(x)} on its \code{"s"} attribute and now draws at it.
            Posterior predictive sampling at test rows of a heteroscedastic
            fit that stores no \code{s.test} is refused rather than drawn at
            the constant.
      \item \code{summary} on a heteroscedastic fit reports \code{mean.s},
            the mean of \eqn{s(x)} over the training observations for each
            draw, in place of \code{sigma}. It previously printed one
            \code{sigma} row holding that same constant, with a posterior
            standard deviation of exactly zero and an undefined R-hat, and
            said nothing about the variance surface. \code{as_draws_array}
            and \code{as_draws_df} resolve the \code{"sigma"} token the same
            way. Homoscedastic fits still summarize \code{sigma}.
      \item A heteroscedastic sampler (\code{variance = TRUE}) accepts
            transactional and per-observation predictor mutation, which it
            refused rather than silently reporting a variance surface
            routed by the previous predictors: every variance tree is now
            re-routed with the mean forest's, and a change that would leave
            a variance leaf empty rolls the whole transaction back, or
            declines the row, instead of installing. \code{$setState}
            validates variance tree occupancy as it has always validated
            it for mean trees, so a state carrying an unoccupied variance
            leaf is refused where it previously installed.
      \item \code{$setPredictor} and \code{$updatePredictor} with
            \code{forceUpdate = TRUE}, and \code{$setCutPoints}, re-route a
            heteroscedastic sampler's variance forest, so the reported
            variance surface and \code{s.test} follow the new predictors
            instead of the ones the forest was built with. Splits the new
            cut grid can no longer address collapse into their parent,
            merging their scale factors by a weight-weighted geometric
            mean.
      \item \code{$setData} works on a heteroscedastic sampler, including a
            replacement data set with a different number of observations:
            the variance forest's per-observation storage resizes with it
            and its trees re-route onto the rebuilt cut grid. A changed
            observation count previously overran that storage and the
            following \code{$run} corrupted memory.
      \item \code{$setModel} no longer unpins the residual standard
            deviation of a heteroscedastic sampler, and a warm start no
            longer overwrites its pinned value; a changed
            \code{resid.prior} recalibrates the variance forest's scale
            leaf, where it addresses the residual scale, rather than being
            silently ignored.
      \item \code{predict} on a heteroscedastic fit restored from a saved
            state - \code{saveRDS}/\code{readRDS}, a new session, or
            \code{$copy()} - reports the fitted \code{s(x)} again. The
            variance forest's saved trees were not part of the serialized
            state, so a re-created sampler replayed empty slots and returned
            \code{s(x)} of exactly zero, with the mean fit unaffected and no
            error or warning. A state whose saved variance block is absent
            while the sampler keeps trees is now refused rather than restored
            from those empty slots.
      \item A variance forest over a predictor with more than 63 levels no
            longer crashes when its trees are flattened, which a
            \code{keepTrees} run did every recorded sample and
            \code{$storeState} did on any heteroscedastic sampler.
      \item \code{$installTrees} refuses a donor and destination that
            disagree on whether they carry a variance forest.
      \item A warm start carries the donor's variance trees instead of
            cold-starting the variance surface, remapping them onto the
            destination's cut grid when the two differ, and
            \code{n.grow.sweeps} initializes against the variance-weighted
            precisions. A donor variance tree that splits outside a
            restricted destination's \code{variance = ~ x} columns, or that
            leaves a leaf empty on the destination's data, is refused.
      \item A warm start from a donor's saved sample now installs that
            sample's own variance surface. \code{$installTrees(donor,
            samples = k)}, and \code{bart2}'s \code{warm.start} on a donor
            fit with \code{keepTrees = TRUE}, paired sample \code{k}'s mean
            forest with the donor's FINAL variance surface whatever
            \code{k} named, so a heteroscedastic warm start kept only half
            the promise the \code{samples} argument makes, and every chain
            of a multi-chain start began from one shared scale surface
            instead of the overdispersed spread. A donor whose saved
            variance block does not cover the named sample is refused
            rather than falling back to its live one.
      \item \code{samplePriorPredictive(type = "ppd")} refuses a
            heteroscedastic sampler rather than returning a homoscedastic
            prior predictive; \code{type = "ev"} is unaffected.
      \item \code{$setPredictor} and \code{$updatePredictor} accept
            transactional (rollback) updates on BCF and multinomial
            samplers, which previously refused them. A replacement is
            installed only if no leaf of any tree of any forest of any
            chain would empty; otherwise the whole change rolls back and
            the call reports \code{FALSE}, exactly as on a single-forest
            sampler.
      \item The per-observation update session, including
            \code{updatePredictorPerObservationJointly} and
            \code{$setPredictor(forceUpdate = "partial")}, accepts BCF and
            multinomial samplers. A row is installed only if it empties no
            leaf of any tree of any forest of any chain, and the returned
            mask reports which rows were.
      \item A BCF forest whose scale multiplier is zero now contributes
            exactly zero weight and exactly zero response for the affected
            rows instead of a near-zero floor, so the treatment surface
            reported at control rows under a fixed glue is exact rather
            than accurate to about seven digits.
      \item A BCF sampler accepts a per-forest, per-observation weight,
            so a caller can exclude rows from one forest's leaf
            conditionals without excluding them from the model, via
            the sampler's \code{$setForestWeights(forest, weights)}.
            \code{forest} indexes from \code{1}, as
            \code{$getCalibration}/\code{$setCalibration} already do
            (a Bayesian causal forest's treatment forest is
            \code{2}). The weight does not ride the sampler's saved
            state, so it is mirrored and reinstalled automatically
            whenever \code{$getPointer} or \code{$setState}
            re-creates the underlying engine.
      \item A multinomial (softmax) sampler's count response can be
            replaced between runs, at the observation and category counts
            it was built with, so the model can serve as a conditional
            inside a larger sampler as every other family already could.
            The trees carry over, fitted to the previous counts. The
            response-side refusals name this channel instead of reporting
            that the sampler fixes its response at creation.
            \code{dbarts:::}-only; no public surface in this release.
      \item A multinomial (softmax) sampler accepts an observation by
            category offset matrix, given at creation or installed
            between runs. It enters each category's linear predictor
            before the softmax, so it shifts the model rather than the
            reported probabilities. \code{dbarts:::}-only; no public
            surface in this release.
      \item The multinomial category offset extends to test data and to
            out-of-sample prediction. Test rows take their own observation
            by category matrix, given at creation or installed between
            runs, and prediction takes one per call for the rows it is
            given rather than reusing either resident offset. A sampler
            that holds a category offset refuses to predict without one,
            an all-zero matrix being how the offset-free surface is asked
            for. \code{dbarts:::}-only; no public surface in this release.
      \item \code{bart2(family = "multinomial")} accepts an \code{offset}:
            an observation by category numeric matrix, entering each
            category's linear predictor before the softmax as it does on
            the \code{dbarts:::}-only channel above. A vector \code{offset}
            is refused by name, since a common per-observation shift is the
            softmax's own null direction and is identically inert.
      \item Probability vectors for the tree prior are validated with the
            standard \code{all.equal} tolerance instead of a hand-picked
            1e-10.
      \item \code{$setResponse} and \code{$setData} check that the
            replacement response lies in the family's support, as sampler
            creation always has: 0/1 for \code{probit} and
            \code{logistic}, an integer category index for
            \code{ordinal}, and a finite non-negative integer count for
            \code{nbinom}. A count response with a negative element
            previously sized the dispersion kernel's histogram from it and
            crashed the R session uncatchably; a binary or ordinal
            response off its support previously drew latents against it
            and returned silent garbage. The flat C API's
            \code{dbarts_sampler_setResponse} takes the same rule, and so
            does sampler creation through \code{dbarts_sampler_create},
            which has no R layer ahead of it.
      \item A heteroscedastic sampler (\code{variance = TRUE}) refuses
            \code{$setResponse} and \code{$setOffset} with
            \code{updateScale = TRUE}. The variance forest's scale leaf is
            calibrated once, against the response transform fixed at
            creation, so re-anchoring that transform left every
            \eqn{s^2(x)}{s^2(x)} measured on a scale the model no longer
            used and the fit ran away while \code{getSigmas} read
            unchanged. \code{updateScale = FALSE}, which pins the
            transform, is the supported response and offset swap.
      \item The refusal to change case weights after creation names the
            family that refuses - \code{probit}, \code{ordinal},
            \code{aft} or \code{nbinom} support no weights at all, and
            \code{logistic} weights are the observation counts its latents
            were built from - instead of telling every one of them about
            "a binary response".
      \item \code{$getForestFits} and \code{$getForestVariableCounts} now
            index \code{forest} from \code{1}, as
            \code{$setForestWeights}, \code{$setForestBasis}, and
            \code{$getCalibration}/\code{$setCalibration} already did (a
            Bayesian causal forest's treatment forest is \code{2});
            \code{0} is refused rather than silently naming the
            prognostic forest. The R5 sampler's forest index is now
            1-based everywhere.
      \item A whole-matrix \code{setPredictor} call (\code{column}
            missing) now carries column names onto the installed
            \code{sampler$data@x}: the replacement's own \code{dimnames}
            when it has them, the sampler's existing column names
            otherwise. Previously the replacement was built as a bare
            matrix with no names at all, so a later by-name column
            reference on that sampler - \code{setPredictor}'s own
            \code{column} argument, or
            \code{\link{updatePredictorPerObservationJointly}}'s
            shared-column match - failed even on a sampler created with
            named columns.
      \item \code{keepSampler} is now honoured on the alternate-family
            \code{bart2} fits (\code{family = "multinomial"},
            \code{"ordinal"}, \code{"nbinom"}, and \code{"hurdle.lognormal"}):
            \code{keepSampler = TRUE} retains \code{$fit} even when
            \code{keepTrees} is left \code{FALSE}. Previously these
            families keyed \code{$fit} on \code{keepTrees} alone, so
            \code{keepSampler} had no effect. \code{keepTrees} still also
            retains the saved forests \code{predict} replays through, and
            ordinal's/nbinom's own \code{thresholds.raw}/\code{dispersion.raw}
            channels, as before.
      \item \code{monotone} combined with an explicit
            \code{proposal.probs = NULL} no longer errors; it is now
            treated the same as an omitted \code{proposal.probs} and
            forced to the birth/death-only proposal set the constraint
            requires.
      \item \code{rbart_vi}'s \code{k} formal now defaults to \code{NULL},
            matching \code{bart2}: a binary fit with no \code{k} takes the
            \code{chi(1.5, 2)} node hyperprior, not a fixed \code{k = 2}.
            The old \code{k = 2.0} default was inert (only an explicit
            \code{k} was ever read), so this is a documentation
            correction, not a behavior change. See \code{?rbart_vi}.
      \item \code{bart2} and \code{rbart_vi} now warn, by name, when a
            supplied argument has no effect under the fit's resolved
            family: the residual-scale trio \code{sigest}/\code{sigdf}/
            \code{sigquant} off every fixed-scale family (\code{"probit"},
            \code{"logistic"}, \code{"ordinal"}, \code{"nbinom"}, the
            discrete-time hazard families, and \code{"multinomial"}),
            \code{dispersion} off any family but \code{"nbinom"}, and
            \code{breaks}/\code{max.rows} off any family but the
            discrete-time hazard ones. Previously these arguments were
            silently dropped. The warning is classed
            \code{dbartsFamilyGatedWarning} (itself a
            \code{dbartsWarning}), the package's first classed condition,
            so a caller can catch, mute with
            \code{suppressWarnings(classes = )}, or promote it with
            \code{options(warn = 2)} without matching on message text.
            \code{bart} never emits it. Every warning the package raises
            now carries a class under \code{dbartsWarning} on this same
            pattern; each is documented on the help page of the function
            whose call raises it.
      \item \code{bart2}'s \code{factors}, \code{missing}, and
            \code{proposal.probs} formals are now genuinely live: each is
            resolved in \code{bart2}'s own frame and always forwarded to
            \code{dbarts}, rather than being dropped from the call (and
            \code{dbarts}'s own default silently applying instead) whenever
            the caller left it unsupplied. \code{proposal.probs}'s default
            is now stated directly in the signature - the named vector
            \code{c(birth_death = 0.5, swap = 0.1, change = 0.4, birth =
            0.5)} \code{dbarts} already used - rather than \code{NULL}.
            \code{rbart_vi} forwards \code{factors}/\code{missing} the same
            way; it has no \code{proposal.probs} formal. Every resolved
            value already matched what applied before, so no fit's draws
            change.
      \item \code{bart2} and \code{rbart_vi} gain \code{storage} and
            \code{updateState} formals, appended at the end of each
            signature, spelled and defaulted exactly as
            \code{\link{dbartsControl}}'s own; both already reached the
            control silently through \code{\dots} and now do so as ordinary,
            documented arguments. Previously an unknown name that happened
            to match a \code{dbartsControl} formal was silently accepted
            and forwarded; now only
            \code{storage} and \code{updateState} are. Because
            \code{updateState} is now a formal, the partial match
            \code{u = } no longer resolves to \code{useQuantiles} on either
            function.
      \item \code{dbartsControl}'s \code{rngSeed} argument and slot are
            renamed \code{seed}, matching \code{dbarts}, \code{bart2},
            \code{xbart}, and \code{rbart_vi}, which already spelled it that
            way; \code{dbarts}'s own convenience \code{seed} argument
            (which overrides the control's) is unchanged. A direct
            \code{dbartsControl(rngSeed = )} call now gets an
            unused-argument error. A consumer that forwards its own
            arguments through a \code{formals()}-derived accept filter -
            stan4bart's \code{bart_args}, and CRAN packages such as WeightIt and
            MatchIt that expose the dbarts surface the same way - drops a
            user's \code{rngSeed} silently instead: the filter no longer
            recognizes the name, so the seed is never forwarded and the fit
            falls back to an unseeded stream.
      \item \code{xbart}'s \code{sigma} argument is renamed \code{sigest},
            matching \code{bart2} and \code{rbart_vi}: the fitting
            functions all spell the starting residual-scale estimate
            \code{sigest}, while the sampler constructors (\code{dbarts},
            \code{dbartsSpec}) keep \code{sigma}. A caller still naming
            \code{sigma} now gets an unused-argument error from
            \code{xbart}.
      \item The sampler's \code{$run} spells its per-call thread count
            \code{n.threads}, matching every other entry point; the old
            \code{numThreads} is an unused-argument error.
      \item \code{bart2} gains \code{tree.prior}, \code{node.prior}, and
            \code{resid.prior} formals (\code{NULL} by default, appended
            after \code{breaks}/\code{max.rows}), taking the same full
            \code{\link{dbartsPriors}} objects \code{\link{dbarts}} already
            does: \code{node.prior = linear(columns = )} or \code{gp(...)}
            leaf models and \code{resid.prior = fixed(value)} are now
            reachable from \code{bart2}, not only from \code{dbarts}.
            Supplying one of these alongside a shorthand argument that would
            otherwise help build the same prior is now an error naming both:
            \code{tree.prior} with any of \code{power}/\code{base}/
            \code{split.probs}/\code{dart}; \code{node.prior} with \code{k}/
            \code{prior.scale}; \code{resid.prior} with \code{sigdf}/
            \code{sigquant}/\code{sigest}. Leaving all three \code{NULL}
            builds the priors from the shorthand arguments exactly as
            before this change, so no existing call's draws move.
            \code{resid.prior} joins the \code{sigest}/\code{sigdf}/
            \code{sigquant} family-gating row: supplied under a
            fixed-unit-scale family, it is silently overwritten with
            \code{fixed(1)} regardless (unchanged engine behavior) and now
            warns, classed \code{dbartsFamilyGatedWarning}, instead of
            saying nothing. \code{tree.prior}/\code{node.prior} are honored
            on every family, including both component fits of
            \code{family = "hurdle.lognormal"}. \code{rbart_vi} does not
            gain these three formals. Two more \code{bart2} abbreviations stop
            working: \code{t=} (ambiguous between \code{test} and
            \code{tree.prior}) and \code{resid.=} (ambiguous between
            \code{resid.dist} and \code{resid.prior}).
      \item \code{xbart} loses its \code{control} argument; \code{n.cuts},
            \code{useQuantiles}, \code{n.thin}, and \code{storage} become
            ordinary formals instead, spelled and defaulted exactly as
            \code{\link{dbartsControl}}'s own (\code{xbart(control = )} is now
            an unused-argument error). A new \code{tree.prior} formal takes
            the same \code{cgm}/DART objects \code{dbarts} does; unlike
            \code{bart2}'s, it does not collide with \code{power}/\code{base}
            - they are \code{xbart}'s own swept grid, overriding a supplied
            object's \code{power}/\code{base} every cell exactly as \code{k}
            already overrides a supplied \code{node.prior}'s \code{k} - but
            does collide with \code{dart}/\code{split.probs}, which would
            only duplicate it. Every default reproduces prior behavior;
            \code{n.th=} no longer partial-matches \code{n.threads}.
      \item \code{xbart}'s default \code{k} is the fixed value 2 for every
            response family. Binary responses previously inherited
            \code{bart2}'s Chi hyperprior, which is held rather than swept
            and is drawn every sweep, so the \code{k} axis collapsed to a
            single cell scored under a shrinkage that moved within the fit;
            with \code{drop = FALSE} the result now carries the length-one
            \code{k} dimension those fits omitted. Binary crossvalidations
            that name \code{k} are unaffected, as is \code{bart2}'s own
            binary default.
      \item \code{bart2} and \code{rbart_vi}'s \code{split.probs} formal now
            defaults to \code{NULL} rather than \code{1 / num.vars}; both
            spellings already resolved to the same uniform prior, so this
            is a default respelling with no change to any fit's draws.
      \item \code{bart} gains \code{subset}, \code{storage}, and
            \code{family} formals, appended at the end of its signature and
            forwarded to \code{\link{dbarts}}; \code{subset} and
            \code{storage} replace what \code{bart} already hardcoded
            (\code{NULL} and \code{"double"} respectively), so no fit's
            draws change by default. \code{family} is \code{c("auto",
            "logistic", "aft")}, first-token default \code{"auto"} exactly
            as before this formal existed; \code{bart} newly reaches
            \code{family = "logistic"} and \code{family = "aft"} directly,
            both packaging as an ordinary \code{"bart"} object. The four
            own-class families - \code{"multinomial"}, \code{"ordinal"},
            \code{"nbinom"}, and \code{"hurdle.lognormal"} - are refused by
            \code{bart} BY NAME, naming the token and pointing to
            \code{bart2}, rather than falling through to
            \code{\link{match.arg}}'s generic \dQuote{should be one of}
            message. No existing \code{bart} abbreviation is affected.
      \item A multi-forest fit built through \code{\link{dbarts}}'s
            \code{forests =} now packages per-forest channels onto the
            returned fit: \code{forestFits} (RESPONSE-scale raw totals,
            trailing \code{forest1}..\code{forestK} margin), the ragged
            \code{glue} amplitudes, and the expanded \code{bases}, which
            together reconstruct the training fit (see \code{\link{bart}}'s
            \sQuote{Value}). \code{extract(type = "forest", forest = ,
            contribution = )} reads them, raw or as each forest's
            contribution; refused by name on any other fit or with
            \code{sample = "test"}.
      \item A \code{\link{forest}} call written inside a
            \code{\link{bart2}}/\code{\link{dbarts}} formula declares a
            second, amplitude-coupled forest without a separate
            \code{forests =} list: \code{z:forest(x1 + x2)} and
            \code{factor(z):forest(x1 + x2)} modulate the forest by
            \code{z}, a \code{(a + b):forest(x1)} compound operand (numeric
            or logical members only) gives it a two-column basis, and
            \code{forest(x1 + x2, basis = ~ z)} is the general named form
            the sugar rewrites to; \code{z * forest(x1)} is refused, naming
            both explicit spellings. Inside a term \code{forest()}'s unnamed
            slot is the split-variable set, symbolically, rather than
            \code{basis}. Every other shape - a term outside a top-level
            additive position or a \code{:} operand, a family the
            multiplier model does not support, \code{test} on a term-bearing
            formula, both a sugar and an explicit \code{basis =}, and a
            \code{forests =} argument alongside a term - is refused by name.
            \code{bart2} gains no \code{forests =} formal; the formula term
            is its only route to a second forest.
      \item \code{bart2(family = "multinomial")} now refuses \code{dart}
            (or a DART \code{tree.prior}), \code{split.probs},
            \code{monotone}, and \code{variance} by name. The K-forest
            engine copies only power/base/proposal-probability fields out
            of the host sampler these resolve against, so a DART tree
            prior, fixed split probabilities, and a variance forest
            previously reached no slot in it at all; \code{monotone} was
            worse than dropped, since its proposal-probability rewrite did
            reach the forest while the direction constraints did not,
            silently fitting a birth/death-only model nobody asked for.
      \item A gaussian-family sampler now refuses \code{variance} together
            with \code{resid.dist = student()} residuals: the two
            previously constructed and ran with no adjudication of whether
            the variance forest's weight-channel routing composes with the
            scale-mixture reweighting \code{student()} installs. A
            \code{\link{rbart_vi}} (grouped) fit was already unreachable
            with \code{variance}, since it declares no such formal.
      \item \code{bart2}'s (and \code{rbart_vi}'s, \code{dbarts}'s)
            \code{offset.test} argument was documented as defaulting to
            \code{NULL}; it has always defaulted to tracking \code{offset}
            (a scalar, or a vector already matching the test row count, is
            reused directly, any other length refused by name), matching
            the \code{offset.test = offset} formal default. Documentation
            only; no behavior changes.
      \item \code{bart2}'s (and \code{dbarts}'s) \code{family =
            "hurdle.lognormal"} documentation described the
            posterior-predictive mean as consuming the positive part's
            per-observation \eqn{\sigma(x)} when a heteroscedastic
            \code{variance = ~x} surface was enabled on it; that state is
            not reachable (a variance forest requires \code{family =
            "gaussian"}, and the hurdle's occupancy component is always
            \code{"probit"}), and the positive part has always been
            homoscedastic, one \eqn{\sigma} per draw recycled across
            observations. Documentation only; no behavior changes.
      \item \code{$getForestVariableCounts} now names its rows from
            \code{data@x}'s column names when the training matrix carries
            them; previously the read was always unnamed regardless of the
            sampler's own predictor names. \code{run()$varcount} is
            unaffected.
      \item \code{bart} now attaches \code{$weights} (and
            \code{$weights.test}, when given) to its returned fit,
            matching \code{bart2} and \code{rbart_vi}. Previously a
            weighted \code{bart} fit carried no \code{$weights}, so
            the pointwise log-likelihood and \code{extract(type =
            "ppd")} silently computed an unweighted likelihood or
            predictive distribution. \code{bart} also now refuses by
            name a call whose sample count thins to zero (e.g.
            \code{ndpost} smaller than \code{keepevery}), rather than
            faulting deep in the result-array reshape.
      \item Every fit built through \code{bart}/\code{bart2} now records
            its residual law (\code{$resid.dist}, \code{"gaussian"} or
            \code{"student"}). The posterior-predictive draw
            (\code{extract(type = "ppd")}) refuses a \code{student()} fit
            by name instead of silently drawing gaussian noise regardless of
            \code{resid.dist}. A fit predating the field is read as
            gaussian, unaffected.
      \item \code{xbart} now enforces the binary weight policy at the R
            level, matching \code{dbarts}/\code{bart2}/\code{rbart_vi}
            instead of falling through to the C bridge's own backstop. A
            \code{family = "probit"} run with weights identically 1 is now
            correctly treated as unweighted; previously it was refused,
            since the bridge backstop that caught it has no such courtesy.
            Other weighted probit calls, and non-integer-count
            \code{family = "logistic"} weights, still refuse but now
            report the message the other fitting functions use rather
            than the bridge's own wording. A binary-family run whose
            starting sigma cannot be estimated reports that same shared
            message (\dQuote{unable to obtain a starting estimate of
            sigma; provide one instead}) instead of a bare \code{lm}
            error, and a \code{resid.prior} supplied under a binary family
            is now always overridden to \code{fixed(1)} rather than only
            when left unsupplied.
      \item \code{xbart}'s \code{n.burn} takes one or two non-negative
            integers instead of up to three; its third element was never
            read. A length-3 (or longer) call still runs, its extra
            elements silently dropped, matching how \code{n.test} has
            always been reduced to its first element regardless of length.
      \item \code{bart2}, \code{xbart}, and \code{rbart_vi} now restore the
            caller's \code{.Random.seed} when a call given a \code{seed}
            argument errors partway through. Previously the error left
            \code{.Random.seed} set from the fixed seed instead of the
            caller's own stream; a successful call was already unaffected.
      \item Error and warning messages across the package now follow one
            consistent style (argument quoting, capitalization,
            punctuation, and refusal phrasing); their information content
            is unchanged.
      \item \code{sampler$setWeights} and \code{sampler$setForestWeights}
            now refuse an \code{NA} or \code{NaN} weight at the R surface by
            the same guard (and, for \code{setForestWeights}, the same
            guard as an infinite weight) that already refused a negative
            one, matching the bridge's own validation. Previously an
            \code{NA}/\code{NaN} case weight was still refused, but by a
            separate check reporting \dQuote{'weights' cannot be NA}; an
            infinite forest weight was not refused until the bridge, deep
            inside an internal \code{tryCatch}.
      \item \code{xbart} now refuses a non-positive \code{n.samples} by
            name, as the other entry points do, instead of running a cell
            with no draws and failing in the loss on an empty sample
            matrix (\dQuote{'x' must be an array of at least two
            dimensions}).
      \item The starting \code{sigest} estimate now warns when it falls back
            to the marginal response sd because the predictor set has
            sparse-backed columns (a \code{Matrix::dgCMatrix}, or a data
            frame mixing dense columns with \code{Matrix::sparseVector}/
            \code{dgCMatrix} ones); the fallback estimate itself is
            unchanged.
      \item A \code{Matrix::sparseVector} or \code{dgCMatrix} column named
            in a \code{formula} now refuses cleanly - the same guarded
            message \code{\link{sparseFactor}} columns already get -
            instead of raw-erroring inside \code{\link{model.frame}}.
      \item Under \code{factors = "indicators"} - \code{bart}'s \code{x}/
            \code{y} interface, \code{bart2(factors = "indicators")}, and
            \code{predict} on such a fit - a test factor declaring more
            levels than the training one is now refused by name. The
            expansion replays the training column's per-level table
            positionally, and the extra levels indexed past the end of it:
            whether the refusal fired at all depended on the values that
            read back, so the same data could be rejected on one run and
            silently fit on the next, with the undeclared levels' rows
            carrying indicators no training row ever showed. Matrix columns
            wider than the training design's are refused the same way, and
            \code{makeModelMatrixFromDataFrame} refuses either directly when
            handed a mismatching \code{drop} pattern. A test factor with
            fewer levels is unaffected.
      \item The \code{$fit} of a \code{bart2(family = "ordinal")} or
            \code{bart2(family = "nbinom")} fit is no longer a host shell:
            it is now the sampler whose engine actually ran, so
            \code{$setResponse}, \code{$run} and the rest work as they do
            on any sampler, and \code{$getDispersion()} on the count fit
            answers with the fit's own \eqn{r}{r} rather than refusing.
            \code{$fit$storeState()} followed by \code{\link{save}}/
            \code{\link{load}} now restores a sampler \code{predict} can
            replay through; previously \code{predict} on a reloaded fit
            failed with a NULL external pointer error for all three
            host-shell families, with no workaround.
            \code{samplerOnly} is supported for these
            two families too, returning the same sampler unrun.
      \item \code{dbarts} accepts \code{family = "multinomial"}, returning an
            ordinary \code{\link{dbartsSampler}} over the K-category softmax
            model that only \code{bart2} could reach before, and only as a
            one-shot fit. The response is an n x K matrix of non-negative
            integer counts, or - the single-trial case - a factor, character
            or integer-code vector one-hot expanded to one. It rides the new
            \code{dbartsData(counts = )} argument, with the per-category
            shifts riding a matrix \code{offset}/\code{offset.test} paired
            with it, so a sampler re-created from its own
            (control, model, data) triple after \code{storeState} and
            \code{\link{save}}/\code{\link{load}} carries its response and
            offsets without further discipline. \code{data@y} holds the
            trials \eqn{n_i}{n_i} the counts imply. \code{dbartsSpec} reaches
            the same family, and a data object carrying counts resolves to it
            from \code{family = "auto"}.
      \item The sampler gains \code{$setCounts},
            \code{$setCategoryOffset} and \code{$setCategoryTestOffset}, the
            three multinomial response channels, previously reachable only as
            unexported internal functions against a classless
            environment. n and K are fixed at creation; each write mirrors
            into the data object, so a re-created sampler carries it.
            \code{$predict} takes the matching per-category offset, one row
            per predicted row. Every channel the softmax gives no meaning to -
            \code{$setResponse}, \code{$setOffset}, \code{$setWeights},
            \code{$setSigma}, \code{$setData}, \code{$setModel},
            \code{$setCalibration}, \code{$setForestWeights},
            \code{$setForestBasis} and \code{$getFitsWithoutOffset} - is
            refused by a message naming the capability, and the channel that
            serves the caller instead where one exists. The predictor family
            (\code{$setPredictor}, \code{$setCutPoints},
            \code{$setTestPredictor}) and the global \code{$setActiveRows}
            stay open; per-forest masking is permanently refused, since a
            category's margin is a log-sum-exp over the other K - 1.
      \item \code{bart2(family = "multinomial")} now constructs its sampler
            directly, through the same dispatch
            \code{dbarts(family = "multinomial")} reaches: one create instead
            of a discarded host sampler plus a second, real K-forest engine.
            \code{$fit} is no longer a host shell for any of the three
            alternate families - it is now the sampler whose engine actually
            ran, mutable and readable on the same terms as a direct
            \code{dbarts(family = "multinomial")} sampler, and
            \code{samplerOnly} is supported (returning that sampler unrun).
            \code{$bc} is removed from all three families' fit objects
            (\code{thresholds.raw} and \code{dispersion.raw} keep their
            \code{keepTrees} gate); \code{predict.bartMultinomial} now reads
            through \code{$fit}, so save/reload/predict round-trips for
            multinomial fits the same way it already does for ordinal and
            nbinom ones. This moves the seeded draw stream for
            \code{bart2(family = "multinomial")} fits only - one
            sampler creation consumes R's RNG stream differently
            than two - so an unseeded multinomial fit made with this version
            will not reproduce one made with an earlier version; every other
            family is unaffected.
      \item Saved trees are now read oldest draw first, counting back from the
            store's write position, by \code{$predict},
            \code{$predictForests}, \code{$getTrees}, \code{$printTrees}, and
            the \code{installTrees}/\code{warm.start} donor pool. Previously
            every reader walked the circular store in slot order, so a
            hand-driven sampler that recorded more than one \code{$run} - a
            burn-in \code{$run(0, n.burn)} followed by
            \code{$run(0, n.samples)} with \code{keepTrees} left on - replayed
            its draws rotated against the \code{train}, \code{sigma} and
            \code{glue} channels those runs returned. A store holding fewer
            draws than \code{n.samples} now reports the draws it holds instead
            of padding the sample axis with slots nothing wrote, and one
            holding none is refused rather than answered with a constant
            surface. \code{bart}, \code{bart2} and \code{rbart_vi} fits were
            never affected. Sampler states carry the store's draw count and so
            move to format version 3; states written by earlier development
            versions are refused.
      \item \code{summary} now has methods for \code{bart2}'s
            \code{"ordinal"}, \code{"nbinom"}, and \code{"hurdle.lognormal"}
            fits. Previously these fell through to \code{summary.default}
            and printed a raw \code{Length}/\code{Class}/\code{Mode} table.
            An ordinal fit's summary reports its thresholds, a
            negative-binomial fit's its dispersion draws, and a hurdle
            fit's both component fits' summaries under their own headers -
            alongside whatever mean-function scale (\code{sigma}, \code{k},
            \code{tau}) each carries, exactly as \code{summary.bart}
            already reports for a plain \code{bart}/\code{bart2} fit.
      \item \code{extract(fit, type = "trees", ...)} no longer corrupts its
            own call. Passing \code{sample}, \code{combineChains},
            \code{forest}, or \code{contribution} - none of which the
            sampler's \code{getTrees} accepts - used to forward into
            \code{getTrees} through R's partial argument matching: a
            numeric \code{sample} silently filtered the returned trees to
            the wrong posterior sample, \code{sample = "train"}/\code{"test"}
            (and the equivalent positional call) died on an internal
            \code{NA}, and the other two raised a raw
            \code{unused argument}. Each is now refused by name, naming the
            \code{chainNums}/\code{sampleNums}/\code{treeNums}/
            \code{current}/\code{newdata} vocabulary \code{getTrees} does
            accept. \code{extract.rbart} had the same defect for
            \code{sample} and \code{combineChains} and is fixed the same
            way.
      \item \code{residuals(fit, sample = "train")} on a \code{bart},
            \code{rbart_vi}, or \code{bart2(family = "hurdle.lognormal")}
            fit no longer raises the raw R error
            \code{formal argument "sample" matched by multiple actual
            arguments}. \code{residuals} is always computed against the
            training response, so a caller-supplied \code{sample} is now
            refused by name instead of colliding with the fixed value
            \code{residuals} itself passes to \code{fitted}.
      \item \code{survivalProbabilities()} on a discrete-time hazard fit
            (\code{bart2(family = "hazard")} and its link aliases) no
            longer fails whenever the training predictors carry column
            names - the ordinary matrix-interface spelling. The
            re-expanded design's appended period column is now named
            \code{"period"} under the same rule the training design used,
            on both the default (training) and matrix \code{newdata}
            branches; previously the column went unnamed and
            \code{predict}'s column-name match refused it, which also
            broke a data.frame training input once it was converted to a
            named matrix.
      \item \code{fitted} and \code{residuals} on an \code{rbart_vi} fit
            whose \code{ranef} array carries no group names - or names that
            do not cover every group - now refuse the call by name instead
            of crashing the R session. The group labels were matched against
            the \code{ranef} dimnames and the unmatched \code{NA} indices
            passed to compiled code, which read outside the array. That entry
            point bounds-checks the indices as well, for a direct
            \code{.Call} that skips the R-level check.
      \item \code{bart2} and \code{rbart_vi} no longer accept a rejection-
            only \code{\dots}: an unrecognized argument now raises R's own
            \dQuote{unused argument} error, the same wall \code{dbarts},
            \code{dbartsSpec} and \code{xbart} already stood behind. The
            retired \code{rngSeed} spelling's pointed message is gone with
            it; \code{seed} is the only spelling now.
      \item \code{bart}'s \code{family} argument refuses every one of
            \code{bart2}'s ten other-family tokens by name, including
            \code{"twopart"} (previously silent, falling through to
            \code{\link{match.arg}}'s generic message); every refusal
            echoes the token as typed, never a resolved alias such as
            \code{"hurdle.lognormal"} or the binary link a hazard token
            remaps to.
      \item \code{dbartsData}'s \code{offset.category}/
            \code{offset.category.test} arguments are gone: a matrix
            \code{offset}/\code{offset.test} now installs the multinomial
            per-category shift when paired with \code{counts}, and is
            refused (naming the shape supplied) on any other data object;
            a flat \code{offset} on a \code{counts}-carrying object is
            refused instead of silently cancelling. The internal slots and
            the \code{$setCategoryOffset}/\code{$setCategoryTestOffset}
            sampler methods are unchanged.
      \item \code{dbartsData} (and every entry point built on it) now
            refuses a survival (\code{Surv}) or multi-column response by
            name instead of raising a false \dQuote{'x' must have the same
            number of observations as 'y'} - the response was being
            silently flattened column-major before the row count was
            compared.
      \item A named \code{test} matrix whose column names do not cover a
            named training design's is an error naming the columns it
            lacks, where the mismatch previously warned and the columns
            were then matched by position - silently pairing the wrong
            predictors. An unnamed \code{test} against a named design still
            matches by position, now warning with the mapping that match
            assumes.
      \item \code{monotone}'s per-element spelling no longer accepts the
            undocumented abbreviations \code{"inc"} and \code{"dec"}
            (\code{"0"}, the character form of the documented unconstrained
            \code{0}, is still accepted); matching the documented tokens is
            case-insensitive, now stated in \code{?dbarts} and
            \code{?bart2}.
      \item \code{bart2}, \code{xbart} and \code{rbart_vi} share one
            message refusing an \code{n.samples} that thins to zero draws;
            \code{dbarts}/\code{dbartsControl} still accept it, for a
            sampler driven by a host loop.
      \item A \code{keepSampler} fit from \code{bart}, \code{bart2},
            \code{rbart_vi} or a \code{family = "hurdle.lognormal"} fit now
            always carries \code{$n.chains}, not only when the sampler
            itself was dropped.
      \item \code{xbart} names a bad \code{family} before ingesting the
            response rather than after, and refuses a non-scalar
            \code{n.threads} by name instead of raising R's coercion error.
      \item \code{names(fit)} on a \code{bart}/\code{bart2} fit no longer
            lists \code{yhat.test}/\code{s.train}/\code{s.test} when the fit
            carries none of them, nor \code{binaryOffset} on a binary fit
            given no \code{offset} (\code{bart}, whose \code{binaryOffset}
            defaults to 0, always carries one). The elements were listed
            with a \code{NULL} value before.
      \item \code{extract}, \code{fitted}, \code{predict}, \code{residuals},
            \code{plotTree}, and \code{survivalProbabilities} on a
            \code{bart2(family = "multinomial"/"ordinal"/"nbinom"/
            "hurdle.lognormal")} fit silently accepted several
            \code{\link{bart}}-family arguments they do not implement -
            \code{extract(forest = )}, \code{extract(combineChains = )},
            \code{fitted(sample = )}, \code{residuals(type = )},
            \code{summary(vars = )} on a multinomial fit, and
            \code{forest}/\code{contribution} on \code{fitted},
            \code{residuals} and \code{predict} - because each
            fell through to an unused \code{\dots}. Every one is now
            refused by name; \code{extract}'s \code{combineChains} and
            \code{fitted}/\code{predict}'s \code{ci.level} are honoured
            instead of ignored, matching \code{\link{bart}}'s own
            generics; \code{plotTree}/\code{survivalProbabilities} on
            these four classes now refuse by name (naming the sampler(s)
            that do carry the trees) rather than raising a raw
            \sQuote{no applicable method}. \code{predict}'s multinomial
            argument \code{offset.category.test} is renamed \code{offset},
            matching \code{predict.bart}'s own new-row shift.
      \item \code{extract(type = "loglik")} is implemented for the four
            \code{bart2} own-class families, and \code{plot} for the
            three that lacked one (\code{bartOrdinal}, \code{bartNegbin},
            \code{bartHurdle}; \code{bartMultinomial}'s gains a second
            panel). \code{as_draws_array}/\code{as_draws_df} are
            registered for all four when \pkg{posterior} is installed,
            previously failing inside \pkg{posterior} with an opaque
            \sQuote{All list elements must be lists themselves}. See
            \code{\link{bart2}}'s \sQuote{Value} section for each
            family's formula and shape.
      \item \code{extract(type = "trees")} on a fit with no saved trees
            but a kept sampler (\code{keepSampler = TRUE} without
            \code{keeptrees}/\code{keepTrees}) silently returned the
            sampler's CURRENT trees, with no \code{chain}/\code{sample}
            column, undocumented; this is \code{\link{plotTree}}'s own
            fallback and is kept, now disclosed under \code{\link{bart}}'s
            \sQuote{Extracting Trees}.
      \item \code{plotTree(fit, sample = )}/\code{plotTree(fit, chain = )}
            partial-matched the wrong formal (\code{sampleNum}/
            \code{chainNum}) instead of raising an error or reaching
            \code{\dots}, silently drawing a different tree than
            requested; both are now refused by name.
      \item The four own-class \code{extract} methods validated
            \code{sample} with a bare \code{match.arg}, giving
            \sQuote{'arg' should be one of ...} where \code{\link{bart}}
            and \code{rbart_vi} fits say \sQuote{sample must be in
            'train', 'test'}; all eight sites now share the one wording.
      \item The starting sigma estimate a linear fit gives \code{\link{dbarts}}
            no longer depends on the host. A response the fit reproduces
            exactly previously returned whatever rounding noise the host's
            BLAS kernel left behind - occasionally exactly 0 - which the
            sampler then refused outright; that estimate is now floored at
            a relative epsilon instead. A design with no residual degrees
            of freedom previously returned an undefined estimate for the
            same reason; it now falls back to the marginal response sd,
            with a warning classed \code{dbartsSigmaFallbackWarning} (a
            \code{dbartsWarning}), exactly as the existing sparse-predictor
            fallback already did - classed
            \code{dbartsSparseSigmaFallbackWarning}, itself a
            \code{dbartsSigmaFallbackWarning} - and the fallback estimate
            is floored the same way.
      \item \code{\link{xbart}} again sweeps a numeric \code{k} grid largest
            to smallest and reports it back in the order it was given, so
            the result no longer depends on the order \code{k} is listed
            in; the rewrite had dropped this sort.
      \item \code{plot.bart}, \code{plot.rbart}, \code{plot.bartMultinomial},
            \code{plot.bartOrdinal}, \code{plot.bartNegbin},
            \code{plot.bartHurdle}, \code{plot.pd2bart}, and
            \code{\link{plotTree}} changed \code{par()} settings (layout,
            margins) without restoring them; each now saves and restores
            the caller's graphical parameters before returning.
      \item \code{\link{bart}}, \code{bart2}, \code{\link{xbart}}, and
            \code{rbart_vi} refuse a formula with no left-hand side
            (\code{~ x1 + x2}) instead of silently fitting an all-zero
            response; \code{\link{dbartsData}} keeps defaulting the
            response to zero, since a composed sampler's response can be
            set after creation.
      \item \code{family = "ordinal"} refuses a response that is neither a
            factor nor integer-valued, instead of treating each distinct
            continuous value as its own category.
      \item The \code{n.threads} argument on every \code{predict} method
            refuses a fractional value instead of silently truncating it,
            matching its documented single-positive-integer contract.
      \item \code{\link{pdbart}} and \code{pd2bart} refuse an
            \code{rbart_vi} fit or one of the four own-class family fits
            (\code{bartMultinomial}, \code{bartOrdinal}, \code{bartNegbin},
            \code{bartHurdle}) by name, instead of the generic
            \sQuote{'x.train' must be a matrix, data.frame, formula, fitted
            bart model, or dbartsSampler} message.
      \item A formula \code{:} or \code{*} term (e.g. \code{y ~ x1 * x2})
            is refused by name at the point the formula is parsed, instead
            of dying inside \code{makeModelMatrix} with R's own
            \sQuote{undefined columns selected}; \code{poly()}, \code{ns()},
            \code{log()}, \code{.}, and \code{offset()} terms are
            unaffected.
      \item \code{predict.bart} and \code{predict.rbart} refuse an
            unrecognized \code{offset.test} argument by name - the
            out-of-sample offset formal on these two methods is
            \code{offset} - instead of silently dropping it into
            \code{\dots}.
      \item Every \code{predict} method now refuses \code{value} by name
            (\code{predict.rbart} already did, with its own more specific
            message naming \code{type}); any other name that is neither a
            formal nor a recognized foreign one now warns - once per call,
            naming every such argument and the method it reached, classed
            \code{dbartsUnusedArgsWarning} (a \code{dbartsWarning}) -
            instead of vanishing silently into \code{\dots}.
      \item \code{Depends} now states \code{R (>= 4.2.0)}, where 0.9-34's
            floor was \code{R (>= 3.1-0)}: the C++20 \code{concept}
            declarations in the sampler's core need a compiler newer than
            R's pre-4.2 reference Windows toolchain (Rtools40, gcc 8.3)
            provides, so a pre-4.2 Windows toolchain cannot build the core
            at all.
      \item \code{\link{bart}}, \code{bart2}, and \code{rbart_vi} no longer
            drive \code{control@printEvery} to 0 under heavy thinning
            (\code{n.thin} greater than \code{printEvery}), which the
            sampler then refused, naming \code{printEvery} - an argument the
            caller may never have set; \code{printEvery} now floors at 1
            instead. \code{dbartsControl}'s validity check refuses
            \code{printEvery = 0} directly (previously accepted as
            \sQuote{non-negative}, contradicting both the sampler and
            \code{\link{dbartsControl}}'s own documentation).
      \item \code{dbarts_sampler_setVerbose}'s \code{printEvery = 0} is now
            refused (\sQuote{printEvery must be at least 1}) instead of
            dividing by zero inside the sampler's print condition.
      \item Gaussian case weights that are negative or not finite are now
            refused wherever weights are installed, including
            \code{dbarts_sampler_setWeights} in the flat C API, which
            previously installed them unchecked; a negative or non-finite
            weight enters a leaf's sufficient statistics as a corrupt
            precision, silently.
      \item \code{DESCRIPTION}'s \code{Authors@R} no longer credits code
            that is not in the tree: the classic engine's adaptive radix
            tree contributor, and the \code{m4} macros
            \code{ax_cxx_namespace_std.m4} and
            \code{ax_func_posix_memalign.m4}, which the build does not use;
            the remaining \code{m4} contributors' credited file lists are
            trimmed to the files still present.
      \item \code{tools/*.R} - the doc-freshness, codoc, and Windows-config
            drift checkers, developer-only scripts with no role in building
            or running the package - are \code{.Rbuildignore}d and no
            longer ship in the source tarball.
      \item \code{resid.prior = fixed(value)} now follows its documentation:
            the residual standard deviation holds at \code{sqrt(value)}, at
            creation and through \code{setModel} alike. The pre-1.0 engine
            ignored the value at creation, holding sigma at its starting
            estimate instead, and installed it as a standard deviation
            (rather than a variance) in \code{setModel}.
      \item A \code{"nbinom"} response is capped at \eqn{10^6}{1e6} and a
            larger count refused, at creation and on every conduit that swaps
            the response, on the R and the flat C surfaces alike. The
            dispersion grid's count histogram is sized from the largest
            count, so a bigger one allocated without bound - \code{y = 1e9}
            asked for 8 GB and took the process down where no R error could
            be raised. Counts at or under the bound draw exactly as before.
      \item \code{dbarts_sampler_setWeights} now refuses a weight change on
            any family but \code{gaussian} and \code{logistic}, the rule the
            R sampler's \code{$setWeights} and the header's own contract
            already state. A \code{"probit"}, \code{"ordinal"}, \code{"aft"}
            or \code{"nbinom"} sampler previously installed a vector nothing
            reads, and a \code{"logistic"} sampler carried a negative weight
            into a division by zero in its working response.
      \item The pointwise log-likelihood channel - the flat C API's
            \code{logLikelihood} results field, and, for a stored zero-weight
            \code{gaussian} fit, R's
            \code{\link{extract}(type = "loglik")} - now reports \code{NaN}
            at any row whose composed weight is zero, rather than a finite
            value or \code{-Inf}: a row not in the model has no likelihood to
            report. Previously a zero-weight \code{gaussian} row reported
            \code{-Inf}, and a zero-weight \code{probit} or \code{ordinal}
            row (reachable only through the active-row mask) reported a
            finite value. A multinomial sampler and any
            amplitude-combined one, at any forest count,
            fill that engine-side field with \code{NaN} instead, with or
            without a mask: their per-observation location is a combination
            of forests the response model does not see.
            \code{extract(type = "loglik")} scores such a fit in R.
      \item \code{xbart}'s weighted RMSE loss was not on the scale of the
            response, and \code{loss = "log"} on a continuous response now
            errors instead of silently computing RMSE. Folds now share the
            full data's cut grid through an internal data handle instead
            of each deriving a grid from its own training rows, so every
            fold bins the predictors identically to a fit on all of the
            data and creation does not re-quantize per fold; numeric
            results shift accordingly.
      \item A heteroscedastic variance forest combined with Student-t
            residuals is now refused at sampler creation on every entrance,
            the flat C API included, rather than only by the modelling
            functions. The two route the residual scale through the same
            per-observation weight channel and are not yet shown to compose;
            because a Student-t sampler reports the \code{gaussian} family,
            the engine's own refusal had let the pair through and fit an
            uncomposed model.
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
