#ifndef R_INTERFACE_SAMPLER_HPP
#define R_INTERFACE_SAMPLER_HPP

// basic sampler creation, running, and accessing

#include <external/Rinternals.h> // SEXP

extern "C" {
  
  SEXP create(SEXP control, SEXP model, SEXP data);
  SEXP run(SEXP fit, SEXP numBurnIn, SEXP numSamples);
  SEXP sampleTreesFromPrior(SEXP fit);
  
  SEXP setData(SEXP fit, SEXP data);
  SEXP setControl(SEXP fit, SEXP control);
  SEXP setModel(SEXP fit, SEXP model);
  
  SEXP predict(SEXP fit, SEXP x_test, SEXP offset_test);
  SEXP setResponse(SEXP fit, SEXP y);
  SEXP setOffset(SEXP fit, SEXP offset);
  SEXP setPredictor(SEXP fit, SEXP x, SEXP forceUpdate, SEXP updateCutPoints);
  SEXP updatePredictor(SEXP fit, SEXP x, SEXP cols, SEXP forceUpdate, SEXP updateCutPoints);
  SEXP setTestPredictor(SEXP fit, SEXP x_test);
  SEXP setTestOffset(SEXP fit, SEXP offset_test);
  SEXP setTestPredictorAndOffset(SEXP fit, SEXP x_test, SEXP offset_test);
  
  SEXP updateTestPredictor(SEXP fit, SEXP x_test, SEXP cols);
   
  SEXP createState(SEXP fit);
  SEXP storeState(SEXP fit, SEXP state);
  SEXP restoreState(SEXP fit, SEXP state);
  
  SEXP getTrees(SEXP fit, SEXP chainIndices, SEXP sampleIndices, SEXP treeIndices);
  SEXP printTrees(SEXP fit, SEXP chainIndices, SEXP sampleIndices, SEXP treeIndices);
  
  SEXP saveToFile(SEXP fit, SEXP fileName);
  SEXP loadFromFile(SEXP fileName);
  
}

#endif

