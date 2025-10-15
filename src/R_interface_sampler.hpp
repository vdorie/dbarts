#ifndef R_INTERFACE_SAMPLER_HPP
#define R_INTERFACE_SAMPLER_HPP

// basic sampler creation, running, and accessing

#include <external/Rinternals.h> // SEXP

extern "C" {
  
  SEXP create(SEXP control, SEXP model, SEXP data);
  SEXP run(SEXP fit, SEXP numBurnIn, SEXP numThreads, SEXP numSamples);
  SEXP sampleTreesFromPrior(SEXP fit);
  SEXP sampleNodeParametersFromPrior(SEXP fit);
  
  SEXP setData(SEXP fit, SEXP data);
  SEXP setControl(SEXP fit, SEXP control);
  SEXP setModel(SEXP fit, SEXP model);
  
  SEXP predict(SEXP fit, SEXP x_test, SEXP offset_test, SEXP numThreads);
  SEXP setResponse(SEXP fit, SEXP y);
  SEXP setOffset(SEXP fit, SEXP offset, SEXP updateScale);
  SEXP setWeights(SEXP fit, SEXP weights);
  SEXP setSigma(SEXP fit, SEXP sigma);
  SEXP setPredictor(SEXP fit, SEXP x, SEXP forceUpdate, SEXP updateCutPoints);
  SEXP updatePredictor(SEXP fit, SEXP x, SEXP cols, SEXP forceUpdate, SEXP updateCutPoints);
  SEXP setCutPoints(SEXP fitExpr, SEXP cutPointsExpr, SEXP colsExpr);
  SEXP setTestPredictor(SEXP fit, SEXP x_test);
  SEXP setTestOffset(SEXP fit, SEXP offset_test);
  SEXP setTestPredictorAndOffset(SEXP fit, SEXP x_test, SEXP offset_test);
  SEXP storeLatents(SEXP fit, SEXP result);

  SEXP startThreads(SEXP fit, SEXP numThreads);
  SEXP stopThreads(SEXP fit);
  
  SEXP updateTestPredictor(SEXP fit, SEXP x_test, SEXP cols);

  SEXP getSigmas(SEXP fit);
  SEXP getSumsOfSquaredResiduals(SEXP fit);
   
  SEXP createState(SEXP fit);
  SEXP storeState(SEXP fit, SEXP state);
  SEXP restoreState(SEXP fit, SEXP state);
  
  SEXP getTrees(SEXP fit, SEXP chainIndices, SEXP sampleIndices, SEXP treeIndices);
  SEXP printTrees(SEXP fit, SEXP chainIndices, SEXP sampleIndices, SEXP treeIndices);
}

#endif

