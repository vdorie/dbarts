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
            (\{1, 2, ..., 20, 30, 50\}) under a renormalized
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
            Test matrices stay dense, and the predictor-mutation surface
            (\code{setPredictor} and relatives, \code{setData}) is fixed
            at creation - make a new sampler instead; state
            serialization and \code{predict} work as usual.
            \code{rbart_vi} and \code{node.prior = linear()} do not yet
            support sparse inputs.
      \item Mixed dense and sparse predictors: a data frame passed to
            the x/y interfaces may hold \code{Matrix::sparseVector} and
            \code{dgCMatrix} columns (assign them into the frame; they
            do not survive \code{data.frame(...)} or \code{I()})
            alongside ordinary numeric and factor columns. Sparse
            columns pass through unexpanded and tier exactly as an
            all-sparse design does, while dense-backed columns keep the
            full feature set: factors still split categorically and
            \code{node.prior = linear()} may designate dense-backed
            columns. Test data frames expand to dense matrices;
            predictor mutation remains fixed at creation and
            \code{rbart_vi} still requires dense inputs.
      \item Categorical predictors support up to 65535 levels, up from
            53: rules on wider factors pool their direction masks per
            tree. In `getTrees` such a rule's `value` is `NA` (its mask
            is not exactly representable in a double) and the
            `directions` column carries the decode as before. Everything
            else composes unchanged, including missing values, saved
            trees with `predict`, and state serialization.
      \item GP leaves: `node.prior = gp(columns, k, lengthscale,
            max.leaf.size)` in `dbarts` fits each bottom node with a
            smooth Gaussian-process function of the designated
            continuous predictor columns under a squared-exponential
            kernel, the prior scale tied to `k` exactly as the other
            leaf priors, giving smooth fits where constant-leaf BART
            steps. `k` may be fixed or given the `chi` hyperprior (the
            binary default), sampled from the drawn functions'
            standardized magnitudes; lengthscales fix at creation, per
            column or by
            a median-distance heuristic; leaves larger than
            `max.leaf.size` fall back to constant fits, so the cubic
            kernel cost is confined to small leaves. The sampler surface
            composes: weights, probit and logistic responses,
            `keepTrees` with `predict` (saved replays reproduce recorded
            test fits exactly), state serialization, and the
            predictor-mutation methods. Function-valued fits ride
            prediction only: `getTrees` reports `NA` leaf values and no
            coefficient columns, and `keepTrees` storage grows with leaf
            sizes. `xbart` accepts the same specification through
            `node.prior`.
      \item Linear leaves: `node.prior = linear(columns, k)` in `dbarts`
            fits each bottom node with an intercept plus a linear term
            in the designated continuous predictor columns (standardized
            internally, every coefficient sharing the `normal(k)`
            prior), so a forest models smoothly-varying coefficients -
            treed regression in the sense of Chipman, George, and
            McCulloch, or the varying-slope structure of Bayesian Causal
            Forests when the covariate is a treatment indicator. The
            full sampler surface composes: weights, probit and logistic
            responses, `keepTrees` with `predict`, state serialization,
            and the predictor-mutation methods. `getTrees` reports leaf
            intercepts in `value` plus one `beta.<column>` slope column
            per covariate, and `plotTree` labels linear leaves with
            their coefficients. `xbart` accepts the same specification
            through a new `node.prior` argument, crossvalidating over
            its `k` grid as usual; folds standardize the leaf
            covariates with full-data constants, matching the final
            fit's calibration the way folds bin on the full-data cut
            grid.
      \item `rbart_vi` with a built-in `prior` (`cauchy`, `gamma`) and
            no `callback` now runs its whole Gibbs sampler - random
            intercepts and their standard deviation included - inside
            the engine: one sampler with `n.chains` chains on worker
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
            with `keepSampler` the `fit` element holds one multi-chain
            sampler instead of one sampler per chain, and save/load of
            such fits continues to work.
      \item The rewritten engine replaces the pre-1.0 engine: every
            fitting function runs on it. Fits are statistically
            equivalent to the pre-1.0 engine's but not draw-for-draw
            identical, so seeded results change. `bart` keeps its
            historical behavior apart from the engine (indicator
            expansion of factors, probit for binary responses).
            `dbartsControl` loses the `engine`, `rngKind`, and
            `rngNormalKind` arguments; explicit generator kinds were a
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
            `dbarts/dbarts.h`: an opaque sampler handle over the rewritten
            engine, created from the R specification objects, with entry
            points for running into caller-owned buffers, the mutation
            surface (response, offset, weights, sigma, predictors, test
            data), latents, prediction, tree storage and introspection,
            and state save and restore. Entry points are reached through
            `R_GetCCallable` and versioned by `DBARTS_C_API_VERSION`.
            It replaces the removed C++ interface in
            `dbarts/R_C_interface.hpp`.
      \item Every chain runs its own generator. `rngSeed` makes results
            reproducible for any number of chains without touching R's
            stream, seeded results no longer depend on the thread
            count, and a single-chain run with a given seed reproduces
            chain 0 of a multi-chain run with the same seed. Without a
            seed, chain generators are seeded from R's stream at
            creation, so `set.seed` alone suffices for reproducibility
            and sampling itself never advances R's stream (the pre-1.0
            engine ran single chains through R's generator and seeded
            multiple chains from the clock).
      \item `resid.prior = fixed(value)` now follows its documentation:
            the residual standard deviation holds at `sqrt(value)`, at
            creation and through `setModel` alike. The
            pre-1.0 engine ignored the value at creation (holding sigma
            at its starting estimate) and installed it as a standard
            deviation in `setModel`.
      \item Factor columns in data frames now default to
            `factors = "categorical"` in `dbarts`, `dbartsData`, `bart2`,
            `xbart`, and `rbart_vi` (the latter two gain the argument):
            each unordered factor stays a single predictor with subset
            splits instead of expanding into indicator columns. Pass
            `factors = "indicators"` for the previous behavior.
      \item Prior specifications are ordinary objects built by the
            constructors in the new exported list `dbartsPriors` (`cgm`,
            `dart`, `normal`, `chi`, `chisq`, `fixed`), usable
            programmatically and passable to the
            `tree.prior`/`node.prior`/`resid.prior` arguments. The
            bare-name forms (`node.prior = normal(chi(1.25))`) keep working
            and now resolve against dbarts' own vocabulary inside those
            arguments, so attaching packages that export `normal` or
            similar cannot change their meaning.
      \item `tree.prior = dart(...)` runs DART (Linero 2018), a Dirichlet
            prior over the split-variable probabilities inducing variable
            selection, with optional sampling of the concentration.
            `bart2`, `rbart_vi`, and `xbart` expose it through a
            `dart` argument taking `TRUE` or a full prior specification,
            along with `split.probs` for fixed, non-uniform variable
            weights. Runs return the per-sample probabilities alongside
            the variable counts, as `varprobs` in the sampler's results
            and on `bart2` and `rbart_vi` fits.
      \item `xbart` is reimplemented as an R-level driver over the sampler
            interface, gaining every sampler feature. Chains now warm-start
            only across
            parameter settings over the same data split, never across
            folds or replications: the retired implementation carried
            chains between splits whose held-out rows the previous
            training set contained, letting slowly-mixing settings (small
            ensembles especially) score against data they had effectively
            seen. The third element of `n.burn` is accordingly ignored.
            Parallelism uses a `parallel` cluster over replications;
            results are reproducible for a fixed `seed` and `n.threads`.
            Also fixed: the weighted RMSE loss was not on the scale of the
            response, and `loss = "log"` on a continuous response now
            errors instead of silently computing RMSE. Folds now share the
            full data's cut grid through an internal data handle instead
            of each deriving a grid from its own training rows, so every
            fold bins the predictors identically to a fit on all of the
            data and creation does not re-quantize per fold; numeric
            results shift accordingly.
      \item `dbarts`, `bart2`, and `xbart` accept `family`: `"auto"` keeps
            the existing dispatch (gaussian for continuous responses, probit
            for 0/1), `"gaussian"` forces a continuous fit even for a 0/1
            response, and `"logistic"` fits a logistic latent model via
            Polya-Gamma augmentation. `xbart`'s built-in binary losses
            transform test predictions through the resolved family's link.
            A logistic
            model accepts positive-integer count weights, a weight-w
            observation contributing its likelihood to the w-th power;
            a probit refuses weights (the weighted probit the pre-1.0
            engine fit was incorrect and is not carried forward), except
            that weights identically 1 are treated as absent, as
            SuperLearner-style callers pass unit weights unconditionally.
            Fits returned by `bart`, `bart2`, and `rbart_vi` record the
            resolved family in a `family` element, and `predict`,
            `extract`, `fitted`, and `plot` transform latent draws
            through that family's link; fits saved by earlier versions
            lack the element and are treated as probit.
      \item `fitted` and `predict` gain a `ci.level` argument: given a
            level in (0, 1) they return a matrix of `est`, `ci.lower`, and
            `ci.upper` per observation rather than the posterior mean. The
            interval kind follows `type` - `"ev"` gives a credible interval
            for the mean (a probability for binary responses), `"ppd"` a
            wider prediction interval that carries the residual noise, and
            `"bart"` a credible interval on the latent scale.
      \item `dbarts` gains a `seed` argument mirroring
            `dbartsControl(rngSeed = )`; `residuals` gains a `type` argument
            that reaches `fitted`, so link-scale residuals are available;
            and `plotTree` is now a generic that plots a tree from a kept
            `bart` or `rbart` fit directly rather than through the stored
            sampler.
      \item MCMC runs are now interruptible: pressing `Ctrl-C` during a
            `bart`, `bart2`, `rbart_vi`, `xbart`, or sampler `run` stops it
            between iterations - joining any worker threads first - and
            returns to the prompt with an error, instead of deferring the
            interrupt until the whole run finishes. The check is polled about
            ten times a second and does not affect the draws.
      \item Under `factors = "categorical"`, splits send a subset of a
            factor's levels (at most 53) down each branch; ordered factors
            are coded as ordinal. Level tables are retained so test data
            are coded identically, with unseen levels an error.
      \item `getTrees` reports a categorical rule's raw direction mask in
            `value` and, when the sampler has categorical predictors, adds
            a `directions` column decoding each rule into one "L"/"R"
            character per level. `plotTree` labels categorical rules with
            the levels sent down the left branch.
      \item Missing values in the predictors are supported. Under
            `missing = "incorporate"`, the default on `dbarts`,
            `dbartsData`, `bart2`, `xbart`, and `rbart_vi`, every split
            rule learns a direction for missing values of its variable,
            so an NA observation follows that rule's chosen branch in
            training and test data alike (Missingness Incorporated in
            Attributes; Twala et al. 2008). `missing = "error"` rejects
            incomplete predictors. The response, weights, and offset must
            be complete, and predictor columns cannot be entirely
            missing; formula inputs previously na.omit-dropped incomplete
            rows silently, and matrix inputs silently mis-coded NAs.
            `getTrees` reports each rule's missing route in a `missing`
            column and `plotTree` appends it to rule labels. `bart` keeps
            rejecting missing data. Fits on complete data are unchanged,
            drawing the same random stream as before.
      \item The rewritten engine implements the full `dbartsSampler`
            interface: every method works, including `keepTrees`,
            `getTrees`, `predict`, `plotTree`, `setWeights`, test offsets,
            `printTrees`, sampling from the prior, `setControl`,
            `setModel`, verbose progress output, and saving/loading
            sampler state (a restored sampler continues bitwise
            identically; the stored state is opaque and incompatible with
            states saved by previous versions). Not supported: weights
            with binary responses, and `setControl` changes to settings
            fixed at creation (chain/tree counts, generators, the cut
            grid).
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
      \item `dbarts`, `bart2`, and `rbart_vi` accept `family = "aft"`: an
            accelerated failure time, log-normal survival model (the forest
            models the conditional mean of log survival time, with `sigma` the
            log-scale residual standard deviation). The response is a
            `survival::Surv` object (right-censoring only) or a two-column
            `(time, status)` matrix or data frame; under `family = "auto"` a
            `Surv` response selects `"aft"` on its own. `rbart_vi` adds a
            random intercept per group to the same model, through its formula
            interface; `dbarts` and `bart2` take the response through their
            x/y interface. None of the three support `weights` or `subset` for
            a survival fit. `predict`, `extract`, and `fitted` return the
            log-time-scale linear predictor, never the time scale; the new
            `survivalProbabilities` generic draws the survival probability
            over a user-supplied time grid instead.
      \item `bart2` accepts `family = "multinomial"`: a K-category softmax
            classifier, K forests coupled through an interleaved Polya-Gamma
            one-vs-rest augmentation. The response is a factor (or character)
            with three or more levels, an n x K matrix of nonnegative integer
            trial counts (one column per category), or the equivalent
            `cbind(c1, ..., cK) ~ x` formula; under `family = "auto"` a factor
            (or character) response with three or more levels is detected and
            routed here, reporting the choice in a one-line message. The
            fit's class is `"bartMultinomial"`, with its own `residuals`,
            `plot`, `summary`, `print`, and `predict` (`type = "ev"`/`"ppd"`)
            methods, separate from `bart`'s.
      \item `family = "auto"` (`dbarts`, `bart2`, `xbart`, and `rbart_vi`
            where applicable) now detects a factor, logical, or character
            response: two levels resolve to probit and three or more to
            multinomial (`bart2` only; elsewhere an error names the
            limitation), each reporting the choice in a one-line message.
            Previously such a response was either rejected outright or,
            for three or more levels given through the x/y interface,
            silently fit as a continuous response on the integer level codes.
    }
  }
  \subsection{C API}{
    \itemize{
      \item The C++ interface in `dbarts/R_C_interface.hpp` is removed with
            the pre-1.0 engine it exposed; packages compiling against
            dbarts use the flat C interface in `dbarts/dbarts.h` instead.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item `copy` now works for samplers with a stored state.
      \item Seeded parallel `rbart_vi` now restores R's stream after
            drawing its chain seeds; the restore assigned a local
            variable instead of the global state.
      \item Fixed reading a saved tree back from binary IO, which recursed with
            the writing routine instead of the reading one.
      \item The sampler's `getSumsOfSquaredResiduals` method now returns
            residual sums of squares on the original response scale as
            documented; it previously descaled by a single factor of the
            data range instead of its square. Binary responses report on
            the latent scale.
      \item A `weights` vector whose length matches neither the response
            nor one is now an error for matrix and sparse inputs; it was
            previously recycled silently. Scalars still broadcast.
      \item `rbart_vi` accepts a custom prior function passed by name;
            any symbol-valued `prior` used to enter the built-in lookup
            and fail.
      \item The `chi(degreesOfFreedom = nu)` hyperprior on `k` now samples
            a chi distribution with `nu` degrees of freedom as documented;
            its Gibbs update used a shape corresponding to `chi(2 * nu - 1)`,
            a dropped Jacobian in the k-squared transform. The binary
            default is relabeled from `chi(1.25, Inf)` to `chi(1.5, Inf)`,
            which draws identically under the corrected shape, so seeded
            default binary fits are unchanged; fits that set `k` to a
            non-default `chi()` now draw the degrees of freedom requested
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
      \item The residual-variance (`sigma`) posterior no longer counts
            zero-weight observations in its degrees of freedom. Such rows
            are documented as ignored and contribute nothing to the
            weighted residual sum of squares or any leaf, but the
            posterior degrees of freedom still counted them, which
            over-concentrated `sigma` and biased it low - severely when
            zero-weight rows were numerous. Fits that pass `weights` of 0
            change accordingly; all-positive-weight fits, including the
            unweighted default, are unaffected.
      \item `rbart_vi`'s `group.by`/`group.by.test` are now looked up by name
            instead of a positional scan that silently grabbed the wrong
            column (typically the response) when the given symbol was absent;
            `predict` also coerces `group.by` to a factor exactly as at fit
            time, so it now recognizes trained groups instead of silently
            re-drawing their random intercepts from the prior.
      \item A factor (or character) response no longer hard-errors
            ("'range' not meaningful for factors") when passed through a
            formula to `dbarts`, `bart2`, or `rbart_vi`; the formula path now
            routes it the same way the x/y interface always has.
      \item Several inputs that previously failed with an unrelated low-level
            error now raise a named one instead: a non-finite response
            (previously "missing value where TRUE/FALSE needed"), zero-row
            data ("subscript out of bounds"), and a non-positive `n.samples`
            ("dim(X) must have a positive length"). A `newdata` variable
            missing from the fit's predictors and sharing its name with a
            builtin function is now named in the error instead of reported as
            "invalid type (builtin)". An unnamed-matrix `newdata` with
            swapped columns now warns naming the assumed column mapping.
    }
  }
}
\section{CHANGES IN VERSION 0.9-31}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item `predict` now accepts a `n.threads` argument and will use native
            threads to parallelize across chains.
    }
  }
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item When running with multiple threads, threads no longer persist after
            `run`ing or `predict`ing and will be re-created on the next call.
            Prior to this, threads were initialized at object creation and
            would run in the background for the life-cycle of the sampler
            object, waiting for work. To prevent frequent stopping/starting of
            threads, new sampler-object functions `startThreads()` and
            `stopThreads()` have been added that can bookend sampling or
            prediction loops. 
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item `rbart_vi` no longer fails when there are multiple out-of-sample
            grouping factors in test data. Bug report thanks to github user
            drcthal.
      \item Should now correctly install on non-Mac OS arm64 systems. Bug
            report thanks to Bob Dröge.
      \item `rbart_vi` now works when `k` is a variable. Bug report thanks
            to Bruno Tancredi.
    }
  }
}
\section{CHANGES IN VERSION 0.9-29}{
  \subsection{BUG FIXES}{
    \itemize{
      \item `dbarts` now works to initialize a sampler for binary outcomes when
            `resid.prior` is left at its default. Bug report thanks to github
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
