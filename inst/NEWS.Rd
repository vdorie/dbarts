\name{NEWS}
\title{dbarts News}
\encoding{UTF-8}
\section{CHANGES IN VERSION 1.0-0}{
  \subsection{NEW FEATURES}{
    \itemize{
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
            their coefficients. Unavailable in `xbart`.
      \item The bartcore engine replaces the classic engine: every
            fitting function runs on it. Fits are statistically
            equivalent to the classic engine's but not draw-for-draw
            identical, so seeded results change. `bart` keeps its
            historical behavior apart from the engine (indicator
            expansion of factors, probit for binary responses).
            `dbartsControl` loses the `engine`, `rngKind`, and
            `rngNormalKind` arguments; explicit generator kinds were a
            classic-engine feature. Saved sampler states from previous
            versions cannot be restored.
      \item A flat C interface for packages compiling against dbarts,
            `dbarts/dbarts.h`: an opaque sampler handle over the bartcore
            engine, created from the R specification objects, with entry
            points for running into caller-owned buffers, the mutation
            surface (response, offset, weights, sigma, predictors, test
            data), latents, prediction, tree storage and introspection,
            and state save and restore. Entry points are reached through
            `R_GetCCallable` and versioned by `DBARTS_C_API_VERSION`.
            It replaces the removed C++ interface in
            `dbarts/R_C_interface.hpp`.
      \item `rngSeed` makes results reproducible for any number of
            chains without touching R's stream, and seeded results no
            longer depend on the thread count. Without a seed,
            multi-chain samplers seed their generators from R's stream,
            so `set.seed` alone suffices for reproducibility (the
            classic engine seeded such chains from the clock).
      \item `resid.prior = fixed(value)` now follows its documentation:
            the residual standard deviation holds at `sqrt(value)`, at
            creation and through `setModel` alike. The
            classic engine ignored the value at creation (holding sigma
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
            `bart2` and `rbart_vi` expose it through a
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
            Binary
            families do not support weights (the weighted probit the
            classic engine fit was incorrect and is not carried forward).
            Fits returned by `bart`, `bart2`, and `rbart_vi` record the
            resolved family in a `family` element, and `predict`,
            `extract`, `fitted`, and `plot` transform latent draws
            through that family's link; fits saved by earlier versions
            lack the element and are treated as probit.
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
      \item The bartcore engine implements the full `dbartsSampler`
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
    }
  }
  \subsection{C API}{
    \itemize{
      \item The C++ interface in `dbarts/R_C_interface.hpp` is removed with
            the classic engine it exposed; packages compiling against
            dbarts use the flat C interface in `dbarts/dbarts.h` instead.
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
      \item The sampler's `getSumsOfSquaredResiduals` method now returns
            residual sums of squares on the original response scale as
            documented; it previously descaled by a single factor of the
            data range instead of its square. Binary responses report on
            the latent scale.
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
