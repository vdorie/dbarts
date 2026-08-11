\name{NEWS}
\title{dbarts News}
\encoding{UTF-8}
\section{CHANGES IN VERSION 1.0-0}{
  \subsection{UPGRADING}{
    A summary of the changes that can break existing code or alter
    results; each is detailed under the other headings.
    \itemize{
      \item Seeded results differ from 0.9-x: the rewritten engine draws
            statistically equivalent but not draw-for-draw identical
            chains. Seeded multi-chain results no longer depend on the
            thread count, and unseeded ones follow \code{set.seed}.
            Sampling no longer advances R's stream (chain generators
            are seeded from it at creation), so R-level draws made
            after a fit also change.
      \item Sampler states saved by earlier versions cannot be restored;
            refit anything you need to continue.
      \item Packages compiling against dbarts must port from the removed
            C++ interface (\code{dbarts/R_C_interface.hpp}) to the flat
            C API in \code{dbarts/dbarts.h}.
      \item Unordered factors split on level subsets by default in
            \code{dbarts}, \code{dbartsData}, \code{bart2}, \code{xbart},
            and \code{rbart_vi}; pass \code{factors = "indicators"} for
            the previous dummy expansion. \code{bart} is unchanged.
      \item Missing predictor values are modeled by default
            (\code{missing = "incorporate"}) rather than dropped
            (formulas) or silently miscoded (matrices); use
            \code{missing = "error"} to reject them. Complete data draw
            as before.
      \item \code{dbartsControl} loses the \code{engine}, \code{rngKind},
            and \code{rngNormalKind} arguments; a \code{weights} vector
            whose length matches neither the response nor one is an error
            instead of recycling; \code{resid.prior = fixed(value)} now
            holds the residual variance at \code{value}; and
            \code{getSumsOfSquaredResiduals} reports on the documented
            scale.
      \item The sampler's \code{setResponse} gains an \code{updateScale}
            argument defaulting to \code{FALSE}: the prior scale set at
            creation stays locked when the response is swapped, as a
            conditional update inside a larger sampler requires.
            Earlier versions always re-derived the prior location and
            scale from each new response; pass
            \code{updateScale = TRUE} for that behavior.
    }
  }
  \subsection{NEW FEATURES}{
    \itemize{
      \item Bayesian causal forests: \code{dbarts} and \code{dbartsSpec}
            accept a \code{treatment} vector of 0/1 assignments, which fits
            the two-forest model
            \eqn{y = a \mu(x) + b_z \tau(x) + \epsilon}{y = a mu(x) +
            b_z tau(x) + eps} instead of a single ensemble. The companion
            \code{moderators} argument restricts the treatment forest to a
            subset of the design, and \code{\link{treatmentForest}} carries
            its tree count, structure prior, prior scales and per-forest
            constraints. The result is an ordinary \code{dbartsSampler}: it
            stores, saves, and re-creates like any other, and the treatment
            vector rides the data object beside the weights it mirrors.
            Gaussian responses only; every option the two-forest model does
            not read is refused at creation rather than ignored. A Bayesian
            causal forest's \code{$setTreatment} mutates the treatment
            vector in place, mirroring it into \code{data@treatment} as
            \code{$setWeights} does for weights; \code{$getForestFits},
            \code{$getBCFGlue}, and \code{$getForestVariableCounts} read the
            prognostic and treatment forests' fitted values, combining
            coefficients, and split counts separately. \code{$setControl} no
            longer discards a sampler's \code{bartcore.*} model
            configuration (the two-forest structure, a variance forest, and
            similar), which used to make the next re-creation refuse loudly.
            A run reports both surfaces for every draw rather than only the
            current one: \code{$run} gains a \code{forestFits} element holding
            each forest's fitted values (n.obs x 2 x n.samples) and a
            \code{glue} element holding the
            \eqn{(a, b_0, b_1)}{(a, b0, b1)} that draw combines them through
            (3 x n.samples), so a caller reads mu, tau and the treatment glue
            from one call instead of driving the sampler one draw at a time.
            Neither element is present for any other model. The
            \code{treatment}, \code{moderators}, and \code{treatmentForest}
            argument names are provisional and are expected to be replaced
            by an engine-vocabulary \code{forests} argument once Bayesian
            causal forests find a permanent home outside \code{dbarts}.
      \item Interaction constraints: the new \code{interactions} argument of
            \code{dbarts} and \code{bart2}, built with \code{\link{interactions}},
            restricts which predictors may jointly shape the fit.
            \code{interactions(max.order = K)} caps the distinct predictors on
            any root-to-leaf path (\code{K = 1} is additive),
            \code{interactions(forbid = list(c("a", "b")))} bars named
            predictors from co-occurring, and \code{interactions(groups = )}
            confines interactions to declared groups; the three combine. The
            constraint is hard and applied per forest, so a multi-forest model
            can hold one forest additive-or-low-order while another stays
            free. An unconstrained fit is bit-for-bit
            identical to previous builds. See \code{?interactions}.
      \item Block-additive constraints: the new \code{blocks} argument of
            \code{dbarts} and \code{bart2}, built with \code{\link{blocks}},
            confines each whole tree to one declared group of predictors, so
            the ensemble is exactly a sum of per-group functions
            \eqn{f = \sum_G f_G}{f = sum_G f_G} (a functional-ANOVA /
            grouped-GAMI decomposition). The partition must be total and
            disjoint; \code{trees.per.group} optionally fixes each group's
            tree count (an even split by default). Applied per forest, and
            composes with \code{interactions} and \code{monotone}. A related
            idiom is
            already available: \code{interactions(groups = )} over a total
            disjoint partition gives the same block-additive property
            adaptively, with each tree's group chosen by the sampler rather
            than fixed; see \code{?interactions}. An unconstrained fit is
            bit-for-bit identical to previous builds. See \code{?blocks}.
      \item Reduced-precision residual storage for large problems:
            \code{storage = "single"} in \code{dbartsControl} (and
            \code{bart2}) stores the running residual in single precision
            (fp32) while accumulating every reduction in double, roughly
            halving the memory traffic of the dominant per-sweep sufficient
            statistic pass. It is opt-in and currently supported only for
            gaussian responses with constant leaves (an error is raised
            otherwise); the leaf draws stay double precision, so posteriors
            are statistically indistinguishable from the default
            double-precision fit. The speedup is largest where memory
            bandwidth binds - large \eqn{n}{n} and multi-chain runs -
            roughly 10-30\% at \eqn{n = 10^6}{n = 1e6} in testing, growing
            with the number of chains. The default \code{storage = "double"}
            is bit-for-bit identical to previous builds. Separately, the
            internal per-observation index storage was narrowed, lowering
            sampler memory use with no change to results. See
            \code{?dbartsControl}.
      \item Semicontinuous two-part (hurdle) responses:
            \code{family = "hurdle.lognormal"} (alias \code{"twopart"}) in
            \code{dbarts} and \code{bart2} fits a non-negative response with
            exact zeros by composing, R-side, an occupancy probit fit of
            \eqn{z = 1\{y > 0\}}{z = 1{y > 0}} over all n observations with a
            lognormal positive-part fit - an ordinary gaussian fit of
            \eqn{\log y}{log y} - over the \eqn{y > 0} subset, glued at
            report time. The two parts share no parameters (no engine code
            is added) and fit at independently derived seeds. Predictions
            default to the natural scale via posterior-predictive Monte
            Carlo, \eqn{E[y \mid x] = P(y > 0 \mid x) e^{f(x) + \sigma^2/2}}{
            E[y | x] = Pr(y > 0 | x) * exp(f(x) + sigma^2/2)},
            heteroscedasticity-aware when the positive part carries a
            \code{variance = ~x} surface; \code{type = "prob"} gives the
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
            and a heteroscedastic positive part are recorded follow-ups).
            See \code{?bart2}.
      \item Heteroscedastic BART: \code{variance = ~ x1 + x2} (or
            \code{variance = TRUE}) in \code{dbarts} and \code{bart2}
            adds a second forest modeling the residual variance surface
            \eqn{s^2(x)} as a product of scaled-inverse-chi-squared
            leaves (HBART; Pratola, Chipman, George, and McCulloch 2020),
            with \code{n.trees.variance}, \code{power.variance}, and
            \code{base.variance} controlling it. The fit gains
            \code{s.train}/\code{s.test} draws of \eqn{s(x)}, and
            \code{predict} returns \eqn{s(x)} on new data. Gaussian
            responses only.
      \item Monotone constraints: \code{monotone = c(x1 = "+", x3 = "-")}
            in \code{dbarts} and \code{bart2} constrains the fit to be
            monotone increasing or decreasing in named predictors
            (monotone BART; Chipman, George, McCulloch, and Shively
            2022). Directions accept the sign glyphs, the words
            \code{"increasing"}/\code{"decreasing"}, or \eqn{+1}/\eqn{-1};
            an unnamed length-p vector assigns them positionally. Only
            numeric and ordered predictors are eligible (a direction on
            an unordered factor is an error). A constraint forces
            birth/death-only proposals and a fixed \code{k = 2}; linear
            and Gaussian-process leaves are not supported under it.
            Unconstrained fits are unchanged. See \code{?dbarts}.
      \item Ordinal (ordered categorical) responses:
            \code{family = "ordinal"} in \code{dbarts} and \code{bart2}
            fits a cumulative probit - probit's truncated-normal latents
            generalized to K ordered categories split by sampled
            cutpoints (\eqn{\gamma_1 = 0}{gamma_1 = 0} pinned; the free
            cutpoints move by a marginal Metropolis update with the
            latents integrated out, under a log-gap
            \eqn{N(0, 1.5^2)}{N(0, 1.5^2)} prior). At K = 2 the model is
            probit exactly. An ordered factor response with three or
            more levels is auto-detected under \code{family = "auto"},
            announced; the level order defines the category order.
            \code{bart2} returns a \code{"bartOrdinal"} fit carrying the
            levels, the n x K category-probability draws, the latent
            draws, and the cutpoint draws, with \code{fitted},
            \code{extract}, \code{predict}, and \code{residuals}
            methods. Weights are not supported; \code{rbart_vi} and
            \code{xbart} refuse ordinal responses (grouped ordinal
            models and ordinal cross-validation losses are recorded
            follow-ups). See \code{?bart2}.
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
            is refused). \code{bart2} returns a \code{"bartNegbin"} fit
            carrying the mean-count draws, the latent draws, and the
            per-draw dispersion, with \code{fitted}, \code{extract},
            \code{predict}, and \code{residuals} methods. Weights are not
            supported (exposure belongs in the offset); \code{rbart_vi}
            and \code{xbart} refuse count responses (grouped
            negative-binomial models, real dispersion, and a Poisson
            family are recorded follow-ups). See \code{?bart2}.
      \item Discrete-time survival hazard responses:
            \code{family = "hazard"} (probit link) and
            \code{family = "hazard.logistic"} in \code{dbarts} and
            \code{bart2} fit a discrete-time hazard model by
            person-period expansion. The response is the same \code{Surv}
            object or two-column \code{(time, status)} the \code{"aft"}
            family takes; each subject is expanded into its at-risk rows
            with an ordinal period column, and the discrete hazard
            \eqn{h(k \mid x) = g(f(x, k) + o)}{h(k | x) = g(f(x, k) + o)}
            is fit as an ordinary binary model on the expanded rows - so a
            hazard fit is byte-identical, draw for draw, to the binary fit
            on the hand-expanded design. The time grid defaults to the
            distinct observed times; \code{breaks} coarsens it (an integer
            quantile count or an explicit boundary vector) and
            \code{max.rows} guards the expansion size. The fit is an
            ordinary \code{"bart"} object whose \code{$family} records the
            binary link and whose \code{$periods} carries the grid;
            \code{\link{survivalProbabilities}} produces survival-curve
            draws (requiring \code{keepTrees}). The matrix interface is
            required, and weights follow the chosen binary family's
            policy. \code{rbart_vi} and \code{xbart} refuse hazard
            responses (grouped/frailty hazard and a cloglog link are
            recorded follow-ups). See \code{?bart2}.
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
            so an estimated \eqn{\nu}{nu} should be posterior-checked.
      \item Sparse predictor matrices: the x/y interfaces of
            \code{dbarts}, \code{dbartsData}, \code{bart}, \code{bart2},
            and \code{xbart} accept a \code{Matrix::dgCMatrix}. Columns
            are ordinal only; per column, the engine stores codes in a
            rank-bitmap layout when at most 20 percent of rows are
            nonzero (roughly a tenth of dense memory at real sparsity)
            and densifies otherwise, so the sparse object is never
            expanded. Missing values enter as stored \code{NA} entries.
            Test matrices stay dense. A \code{dgCMatrix} design accepts
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
            usual. \code{rbart_vi} and \code{node.prior = linear()} do not
            yet support sparse inputs.
      \item Mixed dense and sparse predictors: a data frame passed to
            the x/y interfaces may hold \code{Matrix::sparseVector} and
            \code{dgCMatrix} columns (assign them into the frame; they
            do not survive \code{data.frame(...)} or \code{I()})
            alongside ordinary numeric and factor columns. Sparse
            columns pass through unexpanded and tier exactly as an
            all-sparse design does, while dense-backed columns keep the
            full feature set: factors still split categorically and
            \code{node.prior = linear()} may designate dense-backed
            columns. A mixed design accepts between-sweep column
            mutation - \code{sampler$setPredictor(x, column = j)} over
            dense-backed and sparse-backed columns alike, including a
            sparse \code{sparseFactor} column, and one call may name both
            kinds - with per-observation updates
            (\code{forceUpdate = "partial"}) open on dense-backed columns.
            Whole-matrix replacement (\code{sampler$setPredictor(x)}) is
            accepted too and leaves \code{data@x} a container; a replaced
            sparse column's storage densifies permanently, since every row
            then differs from the implicit value. A test data frame expands to a
            dense matrix, and a predictor container built from other data
            is re-coded against the training factor levels before it is
            used, so its codes mean what the training design's mean.
            \code{rbart_vi} still requires dense inputs.
      \item Categorical predictors support up to 65535 levels, up from
            53: rules on wider factors pool their direction masks per
            tree. In \code{getTrees} such a rule's \code{value} is \code{NA} (its mask
            is not exactly representable in a double) and the
            \code{directions} column carries the decode as before. Everything
            else composes unchanged, including missing values, saved
            trees with \code{predict}, and state serialization.
      \item GP leaves: \code{node.prior = gp(columns, k, lengthscale,
            max.leaf.size)} in \code{dbarts} fits each bottom node with a
            smooth Gaussian-process function of the designated
            continuous predictor columns under a squared-exponential
            kernel, the prior scale tied to \code{k} exactly as the other
            leaf priors, giving smooth fits where constant-leaf BART
            steps. \code{k} may be fixed or given the \code{chi} hyperprior (the
            binary default), sampled from the drawn functions'
            standardized magnitudes; lengthscales fix at creation, per
            column or by
            a median-distance heuristic; leaves larger than
            \code{max.leaf.size} fall back to constant fits, so the cubic
            kernel cost is confined to small leaves. The sampler surface
            composes: weights, probit and logistic responses,
            \code{keepTrees} with \code{predict} (saved replays reproduce recorded
            test fits exactly), state serialization, and the
            predictor-mutation methods. Function-valued fits ride
            prediction only: \code{getTrees} reports \code{NA} leaf values and no
            coefficient columns, and \code{keepTrees} storage grows with leaf
            sizes. \code{xbart} accepts the same specification through
            \code{node.prior}.
      \item Linear leaves: \code{node.prior = linear(columns, k)} in \code{dbarts}
            fits each bottom node with an intercept plus a linear term
            in the designated continuous predictor columns (standardized
            internally, every coefficient sharing the \code{normal(k)}
            prior), so a forest models smoothly-varying coefficients -
            treed regression in the sense of Chipman, George, and
            McCulloch, or the varying-slope structure of Bayesian Causal
            Forests when the covariate is a treatment indicator. The
            full sampler surface composes: weights, probit and logistic
            responses, \code{keepTrees} with \code{predict}, state serialization,
            and the predictor-mutation methods. \code{getTrees} reports leaf
            intercepts in \code{value} plus one \code{beta.<column>} slope column
            per covariate, and \code{plotTree} labels linear leaves with
            their coefficients. \code{xbart} accepts the same specification
            through a new \code{node.prior} argument, crossvalidating over
            its \code{k} grid as usual; folds standardize the leaf
            covariates with full-data constants, matching the final
            fit's calibration the way folds bin on the full-data cut
            grid.
      \item \code{rbart_vi} with a built-in \code{prior} (\code{cauchy}, \code{gamma}) and
            no \code{callback} now runs its whole Gibbs sampler - random
            intercepts and their standard deviation included - inside
            the engine: one sampler with \code{n.chains} chains on worker
            threads instead of per-chain samplers driven
            iteration-by-iteration from R over a PSOCK/FORK cluster.
            Custom prior functions and callbacks keep the previous
            implementation. Results are statistically equivalent but
            not draw-for-draw identical: the response scaling is fixed
            at creation from the data alone (the old loop re-anchored
            it during warmup with the current random effects
            subtracted), and with case weights the conjugate intercept
            update now weights the group contributions (the old loop's
            update ignored weights). Fitted objects keep their shape;
            with \code{keepSampler} the \code{fit} element holds one multi-chain
            sampler instead of one sampler per chain, and save/load of
            such fits continues to work.
      \item The rewritten engine replaces the pre-1.0 engine: every
            fitting function runs on it. Fits are statistically
            equivalent to the pre-1.0 engine's but not draw-for-draw
            identical, so seeded results change. \code{bart} keeps its
            historical behavior apart from the engine (indicator
            expansion of factors, probit for binary responses).
            \code{dbartsControl} loses the \code{engine}, \code{rngKind}, and
            \code{rngNormalKind} arguments; explicit generator kinds were a
            pre-1.0-engine feature. Saved sampler states from previous
            versions cannot be restored.
      \item Sampling is faster: the benchmark run scenarios take 25 to
            30 percent less time per iteration than the pre-1.0 engine,
            the gap widening with the number of observations. The final
            restructurings - the run loop keeps a running residual
            across each sweep, and change and swap proposals score the
            tree prior on the modified subtree alone - reorder floating
            point for speed, trading bit-for-bit reproducibility with
            earlier development builds of this version for posteriors
            that are statistically identical.
      \item A flat C interface for packages compiling against dbarts,
            \code{dbarts/dbarts.h}: an opaque sampler handle over the rewritten
            engine, created from the R specification objects, with entry
            points for running into caller-owned buffers, the mutation
            surface (response, offset, weights, sigma, predictors, test
            data), latents, prediction, tree storage and introspection,
            and state save and restore. Entry points are reached through
            \code{R_GetCCallable} and versioned by \code{DBARTS_C_API_VERSION}.
            It replaces the removed C++ interface in
            \code{dbarts/R_C_interface.hpp}.
      \item Every chain runs its own generator. \code{rngSeed} makes results
            reproducible for any number of chains without touching R's
            stream, seeded results no longer depend on the thread
            count, and a single-chain run with a given seed reproduces
            chain 0 of a multi-chain run with the same seed. Without a
            seed, chain generators are seeded from R's stream at
            creation, so \code{set.seed} alone suffices for reproducibility
            and sampling itself never advances R's stream (the pre-1.0
            engine ran single chains through R's generator and seeded
            multiple chains from the clock).
      \item \code{resid.prior = fixed(value)} now follows its documentation:
            the residual standard deviation holds at \code{sqrt(value)}, at
            creation and through \code{setModel} alike. The
            pre-1.0 engine ignored the value at creation (holding sigma
            at its starting estimate) and installed it as a standard
            deviation in \code{setModel}.
      \item Factor columns in data frames now default to
            \code{factors = "categorical"} in \code{dbarts}, \code{dbartsData}, \code{bart2},
            \code{xbart}, and \code{rbart_vi} (the latter two gain the argument):
            each unordered factor stays a single predictor with subset
            splits instead of expanding into indicator columns. Pass
            \code{factors = "indicators"} for the previous behavior.
      \item Prior specifications are ordinary objects built by the
            constructors in the new exported list \code{dbartsPriors} (\code{cgm},
            \code{dart}, \code{normal}, \code{chi}, \code{chisq}, \code{fixed}), usable
            programmatically and passable to the
            \code{tree.prior}/\code{node.prior}/\code{resid.prior} arguments. The
            bare-name forms (\code{node.prior = normal(chi(1.25))}) keep working
            and now resolve against dbarts' own vocabulary inside those
            arguments, so attaching packages that export \code{normal} or
            similar cannot change their meaning.
      \item \code{tree.prior = dart(...)} runs DART (Linero 2018), a Dirichlet
            prior over the split-variable probabilities inducing variable
            selection, with optional sampling of the concentration.
            \code{bart2}, \code{rbart_vi}, and \code{xbart} expose it through a
            \code{dart} argument taking \code{TRUE} or a full prior specification,
            along with \code{split.probs} for fixed, non-uniform variable
            weights. Runs return the per-sample probabilities alongside
            the variable counts, as \code{varprobs} in the sampler's results
            and on \code{bart2} and \code{rbart_vi} fits.
      \item \code{xbart} is reimplemented as an R-level driver over the sampler
            interface, gaining every sampler feature. Chains now warm-start
            only across
            parameter settings over the same data split, never across
            folds or replications: the retired implementation carried
            chains between splits whose held-out rows the previous
            training set contained, letting slowly-mixing settings (small
            ensembles especially) score against data they had effectively
            seen. The third element of \code{n.burn} is accordingly ignored.
            Parallelism uses a \code{parallel} cluster over replications;
            results are reproducible for a fixed \code{seed} and \code{n.threads}.
            Also fixed: the weighted RMSE loss was not on the scale of the
            response, and \code{loss = "log"} on a continuous response now
            errors instead of silently computing RMSE. Folds now share the
            full data's cut grid through an internal data handle instead
            of each deriving a grid from its own training rows, so every
            fold bins the predictors identically to a fit on all of the
            data and creation does not re-quantize per fold; numeric
            results shift accordingly.
      \item \code{dbarts}, \code{bart2}, and \code{xbart} accept \code{family}: \code{"auto"} keeps
            the existing dispatch (gaussian for continuous responses, probit
            for 0/1), \code{"gaussian"} forces a continuous fit even for a 0/1
            response, and \code{"logistic"} fits a logistic latent model via
            Polya-Gamma augmentation. \code{xbart}'s built-in binary losses
            transform test predictions through the resolved family's link.
            A logistic
            model accepts positive-integer count weights, a weight-w
            observation contributing its likelihood to the w-th power;
            a probit refuses weights (the weighted probit the pre-1.0
            engine fit was incorrect and is not carried forward), except
            that weights identically 1 are treated as absent, as
            SuperLearner-style callers pass unit weights unconditionally.
            Fits returned by \code{bart}, \code{bart2}, and \code{rbart_vi} record the
            resolved family in a \code{family} element, and \code{predict},
            \code{extract}, \code{fitted}, and \code{plot} transform latent draws
            through that family's link; fits saved by earlier versions
            lack the element and are treated as probit.
      \item \code{fitted} and \code{predict} gain a \code{ci.level} argument: given a
            level in (0, 1) they return a matrix of \code{est}, \code{ci.lower}, and
            \code{ci.upper} per observation rather than the posterior mean. The
            interval kind follows \code{type} - \code{"ev"} gives a credible interval
            for the mean (a probability for binary responses), \code{"ppd"} a
            wider prediction interval that carries the residual noise, and
            \code{"bart"} a credible interval on the latent scale.
      \item \code{dbarts} gains a \code{seed} argument mirroring
            \code{dbartsControl(rngSeed = )}; \code{residuals} gains a \code{type} argument
            that reaches \code{fitted}, so link-scale residuals are available;
            and \code{plotTree} is now a generic that plots a tree from a kept
            \code{bart} or \code{rbart} fit directly rather than through the stored
            sampler.
      \item MCMC runs are now interruptible: pressing \code{Ctrl-C} during a
            \code{bart}, \code{bart2}, \code{rbart_vi}, \code{xbart}, or sampler \code{run} stops it
            between iterations - joining any worker threads first - and
            returns to the prompt with an error, instead of deferring the
            interrupt until the whole run finishes. The check is polled about
            ten times a second and does not affect the draws.
      \item Under \code{factors = "categorical"}, splits send a subset of a
            factor's levels (at most 65535, see above) down each branch;
            ordered factors are coded as ordinal. Level tables are retained
            so test data are coded identically, with unseen levels an error.
      \item \code{getTrees} gains, when the sampler has categorical predictors, a
            \code{directions} column decoding each categorical rule into one
            "L"/"R" character per level; such a rule's \code{value} is \code{NA}, as it
            has no scalar cut point. \code{plotTree} labels categorical rules with
            the levels sent down the left branch.
      \item Missing values in the predictors are supported. Under
            \code{missing = "incorporate"}, the default on \code{dbarts},
            \code{dbartsData}, \code{bart2}, \code{xbart}, and \code{rbart_vi}, every split
            rule learns a direction for missing values of its variable,
            so an NA observation follows that rule's chosen branch in
            training and test data alike (Missingness Incorporated in
            Attributes; Twala et al. 2008). \code{missing = "error"} rejects
            incomplete predictors. The response, weights, and offset must
            be complete, and predictor columns cannot be entirely
            missing; formula inputs previously na.omit-dropped incomplete
            rows silently, and matrix inputs silently mis-coded NAs.
            \code{getTrees} reports each rule's missing route in a \code{missing}
            column and \code{plotTree} appends it to rule labels. \code{bart} keeps
            rejecting missing data. Fits on complete data are unchanged,
            drawing the same random stream as before.
      \item The rewritten engine implements the full \code{dbartsSampler}
            interface: every method works, including \code{keepTrees},
            \code{getTrees}, \code{predict}, \code{plotTree}, \code{setWeights}, test offsets,
            \code{printTrees}, sampling from the prior, \code{setControl},
            \code{setModel}, verbose progress output, and saving/loading
            sampler state (a restored sampler continues bitwise
            identically; the stored state is opaque and incompatible with
            states saved by previous versions). Not supported: weights
            with binary responses, and \code{setControl} changes to settings
            fixed at creation (chain/tree counts, generators, the cut
            grid).
      \item \code{setPredictor} accepts \code{forceUpdate = "partial"}, installing a column
            one observation at a time and rolling back only those that would empty
            a node.
      \item \code{updatePredictorPerObservationJointly} applies a \code{"partial"}-style update across several
            samplers at once, installing an observation only when it keeps every
            sampler's nodes non-empty.
      \item The \code{getTrees} method of a sampler accepts \code{current = TRUE} to return
            the live working trees (the current sampler position) even when
            \code{keepTrees} is \code{TRUE}, rather than the saved samples.
      \item \code{getTrees} and \code{extract(., "trees")} accept a \code{newdata} matrix, routing
            those observations through each tree so the \code{n} column counts them
            instead of the training data.
      \item A warm start (\code{bart2}'s \code{warm.start}, or a sampler's
            \code{installTrees}) accepts a donor fit on a different cut grid
            (differing in \code{n.cuts}, \code{useQuantiles}, or its predictor
            values): the donor's splits are remapped onto the
            destination's grid, collapsing any the grid starves, the same
            way replacing a sampler's data remaps its existing splits. A
            donor with a different tree count or DART setting is still
            refused.
      \item \code{dbarts}, \code{bart2}, and \code{rbart_vi} accept \code{family = "aft"}: an
            accelerated failure time, log-normal survival model (the forest
            models the conditional mean of log survival time, with \code{sigma} the
            log-scale residual standard deviation). The response is a
            \code{survival::Surv} object (right-censoring only) or a two-column
            \code{(time, status)} matrix or data frame; under \code{family = "auto"} a
            \code{Surv} response selects \code{"aft"} on its own. \code{rbart_vi} adds a
            random intercept per group to the same model, through its formula
            interface; \code{dbarts} and \code{bart2} take the response through their
            x/y interface. None of the three support \code{weights} or \code{subset} for
            a survival fit. \code{predict}, \code{extract}, and \code{fitted} return the
            log-time-scale linear predictor, never the time scale; the new
            \code{survivalProbabilities} generic draws the survival probability
            over a user-supplied time grid instead.
      \item \code{bart2} accepts \code{family = "multinomial"}: a K-category softmax
            classifier, K forests coupled through an interleaved Polya-Gamma
            one-vs-rest augmentation. The response is a factor (or character)
            with three or more levels, an n x K matrix of nonnegative integer
            trial counts (one column per category), or the equivalent
            \code{cbind(c1, ..., cK) ~ x} formula; under \code{family = "auto"} a factor
            (or character) response with three or more levels is detected and
            routed here, reporting the choice in a one-line message. The
            fit's class is \code{"bartMultinomial"}, with its own \code{residuals},
            \code{plot}, \code{summary}, \code{print}, and \code{predict} (\code{type = "ev"}/\code{"ppd"})
            methods, separate from \code{bart}'s.
      \item \code{family = "auto"} (\code{dbarts}, \code{bart2}, \code{xbart}, and \code{rbart_vi}
            where applicable) now detects a factor, logical, or character
            response: two levels resolve to probit and three or more to
            multinomial (\code{bart2} only; elsewhere an error names the
            limitation), each reporting the choice in a one-line message.
            Previously such a response was either rejected outright or,
            for three or more levels given through the x/y interface,
            silently fit as a continuous response on the integer level codes.
      \item New \code{dbartsSpec} resolves the \code{control}/\code{model}/\code{data} triple and
            family token that describe a sampler, without constructing one,
            for packages that embed dbarts and supply their own design
            matrix. It shares its resolution with \code{dbarts}, so a sampler
            built from its output draws identically to the equivalent
            \code{dbarts} call, and its returned family token is the one
            \code{dbarts_sampler_create} takes. Together with \code{dbartsData} it
            makes the whole single-forest feature set - response families,
            monotonicity, interaction and block constraints, a variance
            forest, robust residuals - reachable from a \code{LinkingTo: dbarts}
            package without reaching into an unexported function.
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
            \code{dbarts_sampler_setTreatment} swaps the treatment vector, and
            \code{dbarts_sampler_forestFits} and \code{dbarts_sampler_bcfGlue}
            read the per-forest fits and the glue the two combine through.
            \code{dbarts_sampler_setResponse} takes an explicit
            \code{updateScale} argument matching
            \code{dbarts_sampler_setOffset}, in place of always re-deriving the
            response transform, and the mutation entries now carry the same
            per-model rules as their R counterparts: a two-forest sampler takes
            a scale-pinned response, offset or weight swap and refuses its whole
            test surface.
      \item The \code{x_test} and \code{offset_test} parameters of
            \code{dbarts_sampler_setTestPredictors}, \code{dbarts_sampler_setTestOffset},
            and \code{dbarts_sampler_predict} are renamed \code{xTest} and \code{offsetTest}
            for the header's camelCase convention. Parameter names are not part
            of the ABI, so callers are unaffected; the rename is baked into the
            1.0-0 signature hash and does not change the version.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
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
            fitting to model only the observed levels. One consequence for
            states saved by an earlier development build: a factor whose
            declared level count crosses the 63/64 inline/pooled rule-mask
            boundary (64 or more declared, fewer observed) now pools its
            masks, and restoring such a state reports a missing
            \code{'tree.masks'} block; states whose tier is unchanged
            restore correctly.
      \item \code{setPredictor} on a mixed or \code{dgCMatrix} design now
            refuses a whole-matrix replacement argument whose shape
            disagrees with the design instead of reinterpreting it
            column-major; previously only \code{is.matrix} arguments (never
            a \code{Matrix} class) had their dimensions checked, so a
            transposed argument of the same total length was silently
            accepted.
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
      \item \code{copy} now works for samplers with a stored state.
      \item Seeded parallel \code{rbart_vi} now restores R's stream after
            drawing its chain seeds; the restore assigned a local
            variable instead of the global state.
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
            unaffected and usable afterward.
      \item A heteroscedastic sampler (\code{variance = TRUE}) now refuses
            transactional and per-observation predictor mutation rather
            than silently reporting a variance surface routed by the
            previous predictors. The refusals on the transactional
            \code{$setPredictor} and \code{$updatePredictor} and on
            per-observation updates are lifted in a later release, as each
            path is routed correctly.
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
      \item \code{$installTrees} refuses a donor and destination that
            disagree on whether they carry a variance forest.
      \item A warm start carries the donor's variance trees instead of
            cold-starting the variance surface, remapping them onto the
            destination's cut grid when the two differ, and
            \code{n.grow.sweeps} initializes against the variance-weighted
            precisions. A donor variance tree that splits outside a
            restricted destination's \code{variance = ~ x} columns, or that
            leaves a leaf empty on the destination's data, is refused.
      \item \code{samplePriorPredictive(type = "ppd")} refuses a
            heteroscedastic sampler rather than returning a homoscedastic
            prior predictive; \code{type = "ev"} is unaffected.
      \item A BCF forest whose scale multiplier is zero now contributes
            exactly zero weight and exactly zero response for the affected
            rows instead of a near-zero floor, so the treatment surface
            reported at control rows under a fixed glue is exact rather
            than accurate to about seven digits.
      \item A BCF sampler accepts a per-forest, per-observation weight,
            so a caller can exclude rows from one forest's leaf
            conditionals without excluding them from the model.
            \code{dbarts:::}-only; no public surface in this release.
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
      \item The sampler retained as \code{$fit} by
            \code{bart2(family = "multinomial", keepTrees = TRUE)} refuses
            mutation. It is only the host shell of the fit, carrying the
            design and priors the K-forest model was built from but none
            of the model, so \code{$setResponse}, \code{$run} and the rest
            silently changed nothing about the fit; its read surface,
            which \code{predict} uses, is unaffected.
    }
  }
}
\section{CHANGES IN VERSION 0.9-31}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item \code{predict} now accepts a \code{n.threads} argument and will use native
            threads to parallelize across chains.
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
