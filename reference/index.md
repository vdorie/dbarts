# Package index

## Model fitting

Fit BART models to fixed data.

- [`bart()`](https://vdorie.github.io/dbarts/reference/bart.md)
  [`bart2()`](https://vdorie.github.io/dbarts/reference/bart.md)
  [`plot(`*`<bart>`*`)`](https://vdorie.github.io/dbarts/reference/bart.md)
  [`predict(`*`<bart>`*`)`](https://vdorie.github.io/dbarts/reference/bart.md)
  [`extract()`](https://vdorie.github.io/dbarts/reference/bart.md)
  [`fitted(`*`<bart>`*`)`](https://vdorie.github.io/dbarts/reference/bart.md)
  [`residuals(`*`<bart>`*`)`](https://vdorie.github.io/dbarts/reference/bart.md)
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
  [`copy(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`show(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`predict(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`setControl(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`setModel(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`setData(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`setResponse(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`setOffset(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`setWeights(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`setSigma(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`setPredictor(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`setCutPoints(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`setTestPredictor(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`setTestPredictorAndOffset(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`setTestOffset(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`printTrees(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`getTrees(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`plotTree(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`startThreads(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  [`stopThreads(`*`<dbartsSampler>`*`)`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  : Class "dbartsSampler" of Discrete Bayesian Additive Regression Trees
  Sampler
- [`dbartsControl()`](https://vdorie.github.io/dbarts/reference/dbartsControl.md)
  : Discrete Bayesian Additive Regression Trees Sampler Control
- [`dbartsData()`](https://vdorie.github.io/dbarts/reference/dbartsData.md)
  : Discrete Bayesian Additive Regression Trees Sampler Data
- [`dbartsPriors`](https://vdorie.github.io/dbarts/reference/dbartsPriors.md)
  : Prior Specification Constructors
- [`updatePredictorPerObservationJointly()`](https://vdorie.github.io/dbarts/reference/updatePredictorPerObservationJointly.md)
  : Jointly Update a Shared Predictor per Observation Across Samplers

## Utilities

- [`plotTree()`](https://vdorie.github.io/dbarts/reference/plotTree.md)
  : Plot a Single Tree From a Fitted BART Model
- [`makeModelMatrixFromDataFrame()`](https://vdorie.github.io/dbarts/reference/makeind.md)
  [`makeind()`](https://vdorie.github.io/dbarts/reference/makeind.md)
  [`makeTestModelMatrix()`](https://vdorie.github.io/dbarts/reference/makeind.md)
  : Make Model Matrix from Data Frame
- [`guessNumCores()`](https://vdorie.github.io/dbarts/reference/guessNumCores.md)
  : Guess Number of Cores
