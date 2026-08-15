# Package index

## Overview

- [`dbarts-package`](https://vdorie.github.io/dbarts/reference/dbarts-package.md)
  : Discrete Bayesian Additive Regression Trees Sampler

## Model fitting

Fit BART models to fixed data.

- [`bart()`](https://vdorie.github.io/dbarts/reference/bart.md)
  [`bart2()`](https://vdorie.github.io/dbarts/reference/bart.md)
  [`plot(`*`<bart>`*`)`](https://vdorie.github.io/dbarts/reference/bart.md)
  [`predict(`*`<bart>`*`)`](https://vdorie.github.io/dbarts/reference/bart.md)
  [`extract()`](https://vdorie.github.io/dbarts/reference/bart.md)
  [`fitted(`*`<bart>`*`)`](https://vdorie.github.io/dbarts/reference/bart.md)
  [`residuals(`*`<bart>`*`)`](https://vdorie.github.io/dbarts/reference/bart.md)
  [`fitted(`*`<bartMultinomial>`*`)`](https://vdorie.github.io/dbarts/reference/bart.md)
  [`predict(`*`<bartMultinomial>`*`)`](https://vdorie.github.io/dbarts/reference/bart.md)
  [`print(`*`<bartMultinomial>`*`)`](https://vdorie.github.io/dbarts/reference/bart.md)
  [`residuals(`*`<bartMultinomial>`*`)`](https://vdorie.github.io/dbarts/reference/bart.md)
  [`plot(`*`<bartMultinomial>`*`)`](https://vdorie.github.io/dbarts/reference/bart.md)
  [`summary(`*`<bartMultinomial>`*`)`](https://vdorie.github.io/dbarts/reference/bart.md)
  : Bayesian Additive Regression Trees
- [`rbart_vi()`](https://vdorie.github.io/dbarts/reference/rbart.md)
  [`plot(`*`<rbart>`*`)`](https://vdorie.github.io/dbarts/reference/rbart.md)
  [`fitted(`*`<rbart>`*`)`](https://vdorie.github.io/dbarts/reference/rbart.md)
  [`extract(`*`<rbart>`*`)`](https://vdorie.github.io/dbarts/reference/rbart.md)
  [`predict(`*`<rbart>`*`)`](https://vdorie.github.io/dbarts/reference/rbart.md)
  [`residuals(`*`<rbart>`*`)`](https://vdorie.github.io/dbarts/reference/rbart.md)
  : Bayesian Additive Regression Trees with Random Effects
- [`xbart()`](https://vdorie.github.io/dbarts/reference/xbart.md) :
  Crossvalidation For Bayesian Additive Regression Trees
- [`pdbart()`](https://vdorie.github.io/dbarts/reference/pdbart.md)
  [`plot(`*`<pdbart>`*`)`](https://vdorie.github.io/dbarts/reference/pdbart.md)
  [`pd2bart()`](https://vdorie.github.io/dbarts/reference/pdbart.md)
  [`plot(`*`<pd2bart>`*`)`](https://vdorie.github.io/dbarts/reference/pdbart.md)
  : Partial Dependence Plots for BART

## The mutable sampler

Create and drive a sampler whose predictors, response, offset, and
weights can change between draws, for use inside larger MCMC schemes.

- [`dbarts()`](https://vdorie.github.io/dbarts/reference/dbarts.md) :
  Discrete Bayesian Additive Regression Trees Sampler
- [`run(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`sampleTreesFromPrior(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`sampleNodeParametersFromPrior(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`growFromRoot(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`copy(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`show(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`predict(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`setControl(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`setModel(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`setData(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`setResponse(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`setOffset(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`setWeights(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`setActiveRows(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`setForestWeights(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`setForestBasis(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`setSigma(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`setPredictor(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`setCutPoints(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`setTestPredictor(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`setTestPredictorAndOffset(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`setTestOffset(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`printTrees(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`getTrees(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`getSigmas(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`getDispersion(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`getLatents(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`getSumsOfSquaredResiduals(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`getFitsWithoutOffset(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`getForestFits(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`getForestAmplitudes(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`getForestVariableCounts(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`getCalibration(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`setCalibration(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`installTrees(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`storeState(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`setState(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`plotTree(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  : Class "dbartsSampler" of Discrete Bayesian Additive Regression Trees
  Sampler
- [`dbartsControl()`](https://vdorie.github.io/dbarts/reference/dbartsControl.md)
  : Discrete Bayesian Additive Regression Trees Sampler Control
- [`dbartsData()`](https://vdorie.github.io/dbarts/reference/dbartsData.md)
  : Discrete Bayesian Additive Regression Trees Sampler Data
- [`dbartsSpec()`](https://vdorie.github.io/dbarts/reference/dbartsSpec.md)
  : Discrete Bayesian Additive Regression Trees Sampler Specification
- [`dbartsPriors`](https://vdorie.github.io/dbarts/reference/dbartsPriors.md)
  : Prior Specification Constructors
- [`interactions()`](https://vdorie.github.io/dbarts/reference/interactions.md)
  : Interaction Constraints for BART
- [`blocks()`](https://vdorie.github.io/dbarts/reference/blocks.md) :
  Block-Additive Constraints for BART
- [`forest()`](https://vdorie.github.io/dbarts/reference/forest.md) :
  Forest Specification for Multi-Forest Models
- [`extract(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/extract.dbartsSampler.md)
  : Extract the Predictor Matrix From a Sampler
- [`samplePriorPredictive()`](https://vdorie.github.io/dbarts/reference/samplePriorPredictive.md)
  : Draw From The BART Prior
- [`sparseFactor()`](https://vdorie.github.io/dbarts/reference/sparseFactor.md)
  : Sparse Unordered Factors
- [`updatePredictorPerObservationJointly()`](https://vdorie.github.io/dbarts/reference/updatePredictorPerObservationJointly.md)
  : Jointly Update a Shared Predictor per Observation Across Samplers

## Survival

Accelerated failure time survival curves.

- [`survivalProbabilities()`](https://vdorie.github.io/dbarts/reference/survivalProbabilities.md)
  : Survival Probability Draws from a Survival Fit

## Diagnostics

Convergence summaries and posterior draw extraction.

- [`summary(`*`<bart>`*`)`](https://vdorie.github.io/dbarts/reference/summary.bart.md)
  [`summary(`*`<rbart>`*`)`](https://vdorie.github.io/dbarts/reference/summary.bart.md)
  [`print(`*`<summary.bart>`*`)`](https://vdorie.github.io/dbarts/reference/summary.bart.md)
  [`as_draws_array(`*`<bart>`*`)`](https://vdorie.github.io/dbarts/reference/summary.bart.md)
  [`as_draws_array(`*`<rbart>`*`)`](https://vdorie.github.io/dbarts/reference/summary.bart.md)
  [`as_draws_df(`*`<bart>`*`)`](https://vdorie.github.io/dbarts/reference/summary.bart.md)
  [`as_draws_df(`*`<rbart>`*`)`](https://vdorie.github.io/dbarts/reference/summary.bart.md)
  : Convergence Diagnostics and Posterior-Package Draws for BART Fits

## Utilities

- [`plotTree()`](https://vdorie.github.io/dbarts/reference/plotTree.md)
  : Plot a Single Tree From a Fitted BART Model
- [`makeModelMatrixFromDataFrame()`](https://vdorie.github.io/dbarts/reference/makeind.md)
  [`makeind()`](https://vdorie.github.io/dbarts/reference/makeind.md)
  [`makeTestModelMatrix()`](https://vdorie.github.io/dbarts/reference/makeind.md)
  : Make Model Matrix from Data Frame
- [`guessNumCores()`](https://vdorie.github.io/dbarts/reference/guessNumCores.md)
  : Guess Number of Cores
