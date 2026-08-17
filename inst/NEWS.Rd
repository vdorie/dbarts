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
      \item The prior defaults for a multi-forest model change, so seeded
            results move for two classes of fit: a \code{"probit"} or
            \code{"logistic"} model carrying a forest with no
            \code{basis}, and any model of three or more forests under any
            family. Declaring \code{forest(sd = )} on every forest - 2 on
            a basis-free one and 1 on each basis one - reproduces the
            previous model exactly.
      \item A \code{forests} list of length one that declares a
            \code{basis} is refused rather than fit as an ordinary
            single-forest model with the declaration dropped.
    }
  }
  \subsection{NEW FEATURES}{
    \itemize{
      \item A model's leaf-prior calibration can be named at creation, in
            response units, rather than inherited from the range of the
            vector the sampler happened to be constructed on:
            \code{prior.scale} (the prior standard deviation of the forest
            total at \code{k = 1}) as a \code{dbartsModel} slot, as
            \code{node.prior = normal(scale = )} on \code{dbarts} and
            \code{xbart} (also \code{sd = }, the same quantity at the
            resolved \code{k}, refused under a \code{k} hyperprior), and as
            a \code{prior.scale} argument of \code{bart}, \code{bart2}, and
            \code{rbart_vi}. It is re-derived at every \code{setModel}, so a
            round trip through the model object holds it. Unset, nothing
            changes.
      \item A sampler reports the calibration in force, per chain and per
            forest, with \code{$getCalibration}: \code{prior.scale},
            \code{prior.sd}, \code{prior.mean}, \code{k},
            \code{k.has.hyperprior}, the response transform, and a
            \code{leaf.model} tag. It reads the engine rather than the model
            object, so it is the first surface to show a re-anchoring channel
            (a \code{setResponse}, \code{setOffset}, or \code{setData} that
            re-derives the response transform) moving what is in force while
            the named intent on the model stays where it was put. For the
            monotone, linear, and Gaussian-process leaf models the reported
            \code{prior.sd} is a stated bound rather than an equality.
            \code{$setCalibration} rewrites its scale half on every chain,
            mid-chain, touching nothing else - unlike \code{setModel}, which
            also re-pins \code{sigma} and refuses a DART sampler outright.
            Its \code{prior.sd} spelling meets the same sampled-\code{k}
            refusal as at creation.
      \item \code{$getCalibration} also reports the five calibration-map
            quantities a multi-forest model runs under and no surface showed:
            the forest's amplitude prior variance or half-Cauchy scale
            (\code{amplitude.prior.variance} and
            \code{amplitude.prior.scale}, exclusive per forest), its
            \code{node.scale.factor} and \code{node.scale.divisor}, and the
            \code{basis.row.norm} the map divides out. The induced prior on a
            multi-forest model's combined location can now be reconstructed
            from the sampler rather than from the call, including after a
            \code{$setForestBasis} swap re-derives the row norm; they are
            \code{NaN} on a forest whose scale the map does not own. The flat
            C \code{dbarts_forest_calibration} gains the matching fields,
            appended below its documented field boundary, so a caller
            compiled against the earlier struct reads as before and one that
            recompiles picks them up.
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
            \code{sd} is now stated per unit of basis row norm, the median
            of that basis's nonzero row norms dividing it out, so an
            unnormalized basis no longer multiplies the prior and a Gaussian
            fit whose basis rows are not unit norm changes. \code{"aft"},
            \code{"ordinal"}, and \code{"nbinom"} are refused at creation,
            naming what each is missing; every option the multi-forest model
            does not read, and every declaration the engine cannot honour, is
            refused the same way. \code{$setForestBasis} is the sole route by
            which a basis changes after creation, at any forest and any width,
            mirroring the matrix into \code{data@bases} as \code{$setWeights}
            does for weights; the amplitudes are preserved and remapped, so a
            width-preserving install leaves every one of them unchanged and a
            width change carries each forest's block to its new offset.
            \code{$getForestFits}, \code{$getForestAmplitudes}, and
            \code{$getForestVariableCounts} read the forests' fitted values,
            amplitudes, and split counts separately;
            \code{$getForestAmplitudes} takes an optional \code{forest},
            returning that forest's own block, since the amplitude vector is
            ragged. \code{$setControl} no longer discards a sampler's
            \code{bartcore.*} model configuration (the multi-forest structure,
            a variance forest, and similar), which used to make the next
            re-creation refuse loudly. A run reports every surface for every
            draw rather than only the current one: \code{$run} gains a
            \code{forestFits} element holding each forest's fitted values
            (n.obs x n.forests x n.samples) and a \code{glue} element holding
            the amplitudes that draw combines them through, stacked
            forest-major (sum(q) x n.samples), so a caller reads every surface
            and the amplitudes from one call instead of driving the sampler
            one draw at a time. Neither element is present for any other
            model. The \code{varcount} element widens on the same principle:
            on a multi-forest sampler it carries a forest axis between the
            predictors and the samples (n.predictors x n.forests x
            n.samples), forest-major within a draw and the prognostic forest
            first, where it reported that forest alone before. A
            single-forest sampler's array is unchanged, and so is the flat C
            API's, whose \code{dbarts_results} declares no forest count and
            still receives the prognostic forest's counts.
            \code{dbartsData} gains a \code{bases} argument naming the
            same matrices directly, replacing the \code{treatment} argument
            and slot of the development series. A list of length one is
            refused by name rather than dropped: a lone forest has no
            second ensemble for its amplitudes to distinguish it from.
      \item The prior defaults a multi-forest model takes when a forest
            says nothing are stated in the terms a caller can defend. A
            forest with no \code{basis} takes a half-Cauchy median of 1
            rather than 2 under \code{"probit"} and \code{"logistic"}
            (unchanged at 2 under a gaussian response, where a drawn
            \code{sigma} absorbs the difference), which puts the prior on
            the latent index where the single-forest binary default
            already puts it. A forest with a \code{basis} takes
            \eqn{\sqrt{2/K}}{sqrt(2/K)} rather than 1 for \eqn{K}{K}
            forests, so declaring more forests no longer widens the prior
            on the combined location without bound: with every forest
            carrying a basis it is now the same at every \eqn{K}{K}, and
            otherwise it is bounded by that value instead of growing.
            \code{forest(sd = )} overrides both and keeps its per-forest
            reading at every \eqn{K}{K}; two forests is the fixed point of
            the second, so nothing at two forests moves.
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
            \code{R_GetCCallable} and versioned by \code{DBARTS_C_API_MAJOR}
            and \code{DBARTS_C_API_MINOR}, with \code{dbarts_apiHash()} as an
            exact signature check.
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
      \item \code{n.grow.sweeps} (\code{bart2}) and the sampler's
            \code{growFromRoot} method build the initial forest by
            XBART-style grow-from-root (He, Yalov and Hahn 2019): each
            requested sweep rebuilds every tree from the root against the
            current residual, in place of a draw from the tree prior,
            before the ordinary Metropolis-Hastings sweeps take over. It is
            opt-in (the default \code{n.grow.sweeps = 0L} leaves fits
            unchanged), mutually exclusive with \code{warm.start}, and
            leaves the targeted posterior unchanged - only where the chain
            starts.
      \item The \code{n.grow.sweeps} initializer now places real categorical
            split rules instead of skipping categorical predictors, so a
            factor-heavy design no longer pays those columns into the
            split-variable probability while drawing nothing out. Below ten
            categories present at a node the node's choice is exact over all
            category subsets; above it a sorted-prefix family carries the
            same total prior mass, so a variable's realized split mass stays
            continuous in its level count.
      \item With \code{dart = TRUE}, categorical predictors now receive split
            counts during the grow sweeps, instead of the structurally zero
            counts they always drew before, so factor signals inform the
            sampled split-variable probabilities from the first regular
            sweep.
      \item The \code{n.grow.sweeps} initializer no longer under-weights cut
            candidates on ordinal predictors that carry missing values. Such
            a candidate stands for the two rules - missing values routed
            left, missing values routed right - that the initializer's
            post-draw coin picks between, and now carries the prior mass of
            both instead of one, so the pair no longer enters the node's
            choice at half its weight relative to not splitting.
      \item A sampler accepts a per-observation 0/1 active-row mask
            (\code{$setActiveRows}), which removes a row from the likelihood
            and from the latent draws of a gaussian, Student-t, \code{probit}
            or \code{ordinal} sampler between draws while keeping its fitted
            value and its leaf occupancy - the row stays in the forest, just
            not in the model this sweep. \code{NULL} clears; an all-ones mask
            is treated as no mask; fractional or \code{NA} values are refused.
            Independent of \code{setWeights} and composes with it.
      \item The active-row mask extends to \code{logistic}, \code{nbinom} and
            \code{aft} samplers.
      \item The active-row mask extends to multinomial samplers, where it
            applies to every category at once; a per-category mask is
            refused, since the softmax margin reads every category's forest.
            \code{dbarts:::}-only; no public surface in this release.
      \item A sampler reports its combined per-observation location on the
            response scale WITHOUT the installed offset
            (\code{$getFitsWithoutOffset()}), which is what an outer block
            conditioning on \eqn{f(x)}{f(x)} needs: \code{run()$train}
            carries the offset, so differencing it against
            \code{$getLatents()} is biased. On a Bayesian causal forest it is
            the only route from R to the combined fit at all. Refused on a
            multinomial sampler, whose channels are per-category softmax
            probabilities and not one additive location. There is no
            with-offset twin; the caller adds back the offset it installed.
            \code{?dbartsSampler} gains a table of which method reports which
            scale, with or without the offset, and from the current state or
            from stored draws.
      \item A \code{"nbinom"} sampler reports its dispersion \eqn{r}{r} as a
            first-class per-draw channel: \code{run()} gains a
            \code{dispersion} element, shaped exactly as \code{sigma}, and
            \code{$getDispersion()} reads the value in force mid-sweep. Both
            are pure reads of state the sweep already settled on, so neither
            consumes random numbers or changes a draw. A host driving the
            sampler one sweep at a time no longer has to serialize the whole
            state block to recover \eqn{r}{r}, which is how
            \code{bart2(family = "nbinom")} itself now reads it. Both are
            \code{NULL} on every other family - a channel absent rather than
            an error - so a caller tests \code{!is.null()} rather than
            comparing.
      \item A sampler carrying grouped random effects - the one
            \code{\link{rbart_vi}} returns in \code{$fit} - accepts
            \code{$setResponse} and \code{$setOffset} at the pinned scale
            (\code{updateScale = FALSE}, the default) instead of refusing
            every response-side swap, so a grouped BART fits inside a larger
            Gibbs sampler on the same terms as an ungrouped one. The group
            effects and their scale are carried across the swap as the Gibbs
            block they are. \code{updateScale = TRUE} is now refused on a
            gaussian, Student-t or \code{"aft"} grouped sampler, naming the
            effects and their scale: it re-anchored the response transform
            underneath both while converting neither, silently restating them
            in response units - and on \code{$setOffset} it did so without
            any refusal at all. A grouped \code{"probit"} or
            \code{"logistic"} sampler still accepts \code{TRUE}, where it is
            the no-op it has always been. \code{$setData} remains refused.
      \item A \code{"nbinom"} response is capped at \eqn{10^6}{1e6} and a
            larger count refused, at creation and on every conduit that swaps
            the response, on the R and the flat C surfaces alike. The
            dispersion grid's count histogram is sized from the largest
            count, so a bigger one allocated without bound - \code{y = 1e9}
            asked for 8 GB and took the process down where no R error could
            be raised. Counts at or under the bound draw exactly as before.
      \item The \code{$fit} of a \code{bart2(family = "ordinal")} or
            \code{bart2(family = "nbinom")} fit refuses mutation, as a
            multinomial fit's already did. It is only the host shell of the
            model - it carries the design and priors the model was built
            from but none of the model itself - so \code{$setResponse},
            \code{$run} and the rest used to change a sampler nothing reads.
            \code{$getDispersion()} on the count shell is refused for the
            same reason: the shell's own \eqn{r}{r} is not the fit's.
      \item \code{dbarts_sampler_setWeights} refuses a weight change under
            any family but gaussian, the rule the R sampler's
            \code{$setWeights} and the header's own contract already stated.
            A \code{"probit"}, \code{"ordinal"}, \code{"aft"} or
            \code{"nbinom"} sampler installed a vector nothing reads, and a
            \code{"logistic"} one carried a negative weight into a division
            by zero in its working response.
      \item \code{$getLatents} is documented per family. It reports a
            LOCATION for \code{"probit"}, \code{"ordinal"} and \code{"aft"}
            and a PRECISION for \code{"logistic"}, \code{"nbinom"} and a
            Student-t residual distribution; the flat C API's
            \code{dbarts_sampler_getLatents} comment said \dQuote{gaussian}
            where it meant \dQuote{augments nothing}, and a gaussian-family
            sampler with \code{resid.dist = student()} does report latents.
            Documentation only; no value changes.
      \item The pointwise log-likelihood channel - the flat C API's
            \code{logLikelihood} results field, and, for a stored zero-weight
            \code{gaussian} fit, R's \code{pointwiseLogLikelihood}/
            \code{\link{extract}(type = "loglik")} - now reports \code{NaN}
            at any row whose composed weight is zero, rather than a finite
            value or \code{-Inf}: a row not in the model has no likelihood to
            report. Previously a zero-weight \code{gaussian} row reported
            \code{-Inf}, and a zero-weight \code{probit} or \code{ordinal}
            row (reachable only through the active-row mask above) reported a
            finite value. Multinomial and Bayesian causal forest fits report no
            pointwise log-likelihood at all, with or without a mask.
      \item The augmentation draws a non-gaussian family runs inside each
            sweep are callable from R: \code{\link{dbartsDrawLatents}} draws
            the latent at a supplied fit, and
            \code{\link[=dbartsDrawLatents]{dbartsWorkingResponse}} turns it
            into the quantity a host regresses on - which for
            \code{"logistic"} and \code{"nbinom"} is not the latent at all,
            those being precisions. The draws come from R's own stream, so
            \code{\link{set.seed}} reproduces them and they compose with any
            other R draw; a sampler's chain generator is untouched. The
            \code{fit} argument is the location WITHOUT the offset, as
            \code{$getFitsWithoutOffset()} reports it.
      \item \code{\link{dbartsValidateComposition}} runs simulation-based
            calibration over a one-sweep step the caller supplies, for a model
            the caller states as a prior draw and a simulation, and reports per
            scalar functional whether the composed kernel targets the posterior
            it claims. It catches the defect class in which a posterior mean is
            reported where a draw belongs, and needs no sampler: the functionals
            may be derived, the ranks are tie-broken, and the uniformity band is
            corrected over the functionals ranked in the call.
      \item Six named compositions ship as a vignette,
            \code{vignette("dbarts-as-a-component")}: data augmentation written
            in R, an outer block on the offset channel and its calibration
            check, a K-sampler decomposition sharing one prior budget, a latent
            predictor shared across samplers, an outer-owned sigma, and the
            incremental read. Each is a runnable outer sampler and each is also
            a test cell. \code{\link{dbarts-embedding}} is the reference page
            behind them: which channel an outer block writes, which fit it
            reads back, and whose prior is in force after a swap.
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
            \code{dbarts_sampler_forestFits} and the amplitude pair
            (\code{dbarts_sampler_numForestAmplitudes},
            \code{dbarts_sampler_forestAmplitudes}) read the per-forest fits
            and the amplitudes the two combine through.
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
            a consumer an exact signature-lockstep check alongside the
            major/minor handshake; and \code{dbarts_sampler_setForestWeights},
            \code{dbarts_sampler_setActiveRows},
            \code{dbarts_sampler_forestCalibration} and
            \code{dbarts_sampler_setForestPriorScale} reach the per-forest
            weight, active-row mask and leaf-prior calibration channels from C.
      \item A C host reads the negative-binomial dispersion \eqn{r}{r} without
            an R round trip: \code{dbarts_results} gains a \code{dispersion}
            field shaped as \code{sigma}, and \code{dbarts_sampler_dispersion}
            reads the value in force mid-sweep, returning 0 and touching
            nothing off a family carrying none. The field is appended, so a
            consumer compiled against the smaller struct is unaffected.
      \item \code{dbarts_drawLatents} and \code{dbarts_workingResponse} mirror
            \code{\link{dbartsDrawLatents}} and \code{dbartsWorkingResponse}
            over caller-supplied arrays: the same laws, the same family rules
            and response-support validation, drawn from R's own stream.
      \item \code{DBARTS_C_API_HASH} is re-baked for the three entry points
            above. They are append-only, so no existing call site changes, but
            a consumer checking \code{dbarts_apiHash()} against the header
            constant must be rebuilt - the lockstep handshake as designed.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item A leaf containing only observations of weight zero now counts as
            empty for the tree-move vetoes, so no tree can keep a leaf that
            no observation in the likelihood reaches. Previously such a leaf
            was treated as occupied and drew its value from the prior. Fits
            with no weights, or with strictly positive weights, draw exactly
            as before.
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
      \item \code{$setPredictor} and \code{$updatePredictor} accept
            transactional (rollback) updates on BCF and multinomial
            samplers, which previously refused them. A replacement is
            installed only if no leaf of any tree of any forest of any
            chain would empty; otherwise the whole change rolls back and
            the call reports \code{FALSE}, exactly as on a single-forest
            sampler. Multinomial's rejection rate was measured as ungated
            stress rather than covered by a verdict.
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
      \item The sampler retained as \code{$fit} by
            \code{bart2(family = "multinomial", keepTrees = TRUE)} refuses
            mutation. It is only the host shell of the fit, carrying the
            design and priors the K-forest model was built from but none
            of the model, so \code{$setResponse}, \code{$run} and the rest
            silently changed nothing about the fit; its read surface,
            which \code{predict} uses, is unaffected.
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
            retains \code{$bc} and the family's raw draw channels, as
            before.
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
            correction, not a behavior change; \code{man/rbart.Rd} is
            fixed to match.
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
            \code{bart} never emits it.
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
            documented arguments. \code{\dots} is otherwise rejection-only on
            both: any other name it receives is now an error naming the
            argument and, when one is close, suggesting the nearest formal
            (\code{bart2(n.tres = 5)} suggests \code{n.trees}). Previously an
            unknown name that happened to match a \code{dbartsControl}
            formal was silently accepted and forwarded; now only
            \code{storage}, \code{updateState}, and (until the
            \code{rngSeed} -> \code{seed} rename) \code{rngSeed} are.
            Because \code{updateState} is now a formal, the partial match
            \code{u = } no longer resolves to \code{useQuantiles} on either
            function.
      \item \code{dbartsControl}'s \code{rngSeed} argument and slot are
            renamed \code{seed}, matching \code{dbarts}, \code{bart2},
            \code{xbart}, and \code{rbart_vi}, which already spelled it that
            way; \code{dbarts}'s own convenience \code{seed} argument
            (which overrides the control's) is unchanged. A direct
            \code{dbartsControl(rngSeed = )} call now gets an
            unused-argument error. \code{bart2}/\code{rbart_vi}'s
            \code{\dots} still recognizes \code{rngSeed} by name, but only
            to refuse it - the first row of the retired-name table:
            \code{"unknown argument 'rngSeed': 'rngSeed' was renamed to
            'seed'"}. A consumer that forwards its own arguments through a
            \code{formals()}-derived accept filter - stan4bart's
            \code{bart_args}, and CRAN packages such as WeightIt and
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
      \item \code{variance}'s three companion knobs - \code{n.trees.variance},
            \code{power.variance}, and \code{base.variance} - are removed
            from \code{bart2}, \code{dbarts}, and \code{dbartsSpec}, replaced
            by a new \code{\link{varianceForest}(vars = , n.trees = , base =
            , power = )} constructor. \code{variance} now accepts either its
            existing shorthand (\code{NULL}/\code{FALSE} for none,
            \code{TRUE}/a one-sided formula/a character or integer column
            selector for the predictors) or a \code{varianceForest} object,
            whose \code{vars} resolves through the SAME selector the
            shorthand uses and whose \code{n.trees}/\code{base}/\code{power}
            replace the shorthand path's own defaults (\code{40} trees, the
            mean forest's tree-structure prior). \code{variance = TRUE,
            n.trees.variance = 10} becomes \code{variance =
            varianceForest(n.trees = 10)}; a restricted \code{variance = ~
            x1, n.trees.variance = 10} becomes \code{variance =
            varianceForest(vars = ~ x1, n.trees = 10)}. A caller who wrote
            any of the three removed names now gets an unused-argument
            error (\code{bart2}'s \code{\dots} reports \code{"unknown
            argument"} for these three, with no nearest-formal suggestion, since
            none of the surviving formal names is close enough for
            \code{agrep} to propose one). A \code{formals()}-derived accept
            filter - stan4bart's \code{bart_args}, and CRAN packages such as
            WeightIt and MatchIt that expose the dbarts surface the same way
            - drops a user's \code{n.trees.variance}/\code{power.variance}/
            \code{base.variance} silently instead, since the filter no
            longer recognizes the names.
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
            gain these three formals (its model-surface gap is a separate,
            recorded ticket). Two more \code{bart2} abbreviations stop
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
