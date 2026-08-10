#ifndef R_INTERFACE_BARTCORE_HPP
#define R_INTERFACE_BARTCORE_HPP

// internal bridge to the generalized core (bartcore/); see
// docs/design/core-generalization.md

#include <external/Rinternals.h> // SEXP

extern "C" {

SEXP bartcore_create(SEXP control, SEXP model, SEXP data, SEXP family);
SEXP bartcore_createDataHandle(SEXP control, SEXP data,
                               SEXP leafCovariateColumns);
SEXP bartcore_createFromHandle(SEXP control, SEXP model, SEXP data,
                               SEXP handle, SEXP trainRows, SEXP testRows,
                               SEXP family, SEXP columns);
SEXP bartcore_createBCF(SEXP control, SEXP model, SEXP data, SEXP z,
                        SEXP bcfParams, SEXP moderators,
                        SEXP muInteractions, SEXP tauInteractions,
                        SEXP muBlocks, SEXP tauBlocks);
SEXP bartcore_createMultinomial(SEXP control, SEXP model, SEXP data,
                                SEXP labels, SEXP numCategories);
SEXP bartcore_createMultinomialCounts(SEXP control, SEXP model, SEXP data,
                                      SEXP counts, SEXP numCategories);
SEXP bartcore_setTreatment(SEXP ptr, SEXP z);
SEXP bartcore_setForestWeights(SEXP ptr, SEXP forest, SEXP weights);
SEXP bartcore_getBCFGlue(SEXP ptr);
SEXP bartcore_getForestFits(SEXP ptr, SEXP forest);
SEXP bartcore_getForestVariableCounts(SEXP ptr, SEXP forest);
SEXP bartcore_run(SEXP ptr, SEXP numBurnIn, SEXP numSamples);
SEXP bartcore_runWithCallback(SEXP ptr, SEXP numBurnIn, SEXP numSamples,
                              SEXP results, SEXP callback, SEXP rho);
SEXP bartcore_setOffset(SEXP ptr, SEXP offset, SEXP updateScale);
SEXP bartcore_setResponse(SEXP ptr, SEXP y, SEXP updateScale);
SEXP bartcore_setSigma(SEXP ptr, SEXP sigma);
SEXP bartcore_setData(SEXP ptr, SEXP data);
SEXP bartcore_setTestPredictor(SEXP ptr, SEXP x_test);
SEXP bartcore_setTestOffset(SEXP ptr, SEXP offset_test);
SEXP bartcore_setTestPredictorAndOffset(SEXP ptr, SEXP x_test,
                                        SEXP offset_test);
SEXP bartcore_setWeights(SEXP ptr, SEXP weights);
SEXP bartcore_setPredictor(SEXP ptr, SEXP x, SEXP forceUpdate,
                           SEXP updateCutPoints);
SEXP bartcore_updatePredictor(SEXP ptr, SEXP x, SEXP columns, SEXP forceUpdate,
                              SEXP updateCutPoints);
SEXP bartcore_updatePredictorPerObservation(SEXP ptr, SEXP x, SEXP column);
SEXP bartcore_updatePredictorPerObservationJointly(SEXP ptrs, SEXP x,
                                                   SEXP columns);
SEXP bartcore_setCutPoints(SEXP ptr, SEXP cutPoints, SEXP columns,
                           SEXP currentPredictors);
SEXP bartcore_getSigmas(SEXP ptr);
SEXP bartcore_isValidPointer(SEXP ptr);
SEXP bartcore_getLatents(SEXP ptr, SEXP result);
SEXP bartcore_predict(SEXP ptr, SEXP x_test, SEXP offset_test);
SEXP bartcore_getTrees(SEXP ptr, SEXP chainNums, SEXP sampleNums,
                       SEXP treeNums, SEXP current, SEXP newdata,
                       SEXP trainingData, SEXP forest);
SEXP bartcore_storeState(SEXP ptr);
SEXP bartcore_setState(SEXP ptr, SEXP state, SEXP currentPredictors);
SEXP bartcore_installForests(SEXP ptr, SEXP donorState, SEXP samples);
SEXP bartcore_sampleTreesFromPrior(SEXP ptr);
SEXP bartcore_sampleNodeParametersFromPrior(SEXP ptr);
SEXP bartcore_growFromRoot(SEXP ptr, SEXP numSweeps);
SEXP bartcore_printTrees(SEXP ptr, SEXP chainNums, SEXP sampleNums,
                         SEXP treeNums);
SEXP bartcore_setControl(SEXP ptr, SEXP control);
SEXP bartcore_setModel(SEXP ptr, SEXP model, SEXP control, SEXP data);
SEXP bartcore_getSumsOfSquaredResiduals(SEXP ptr);

}

#endif // R_INTERFACE_BARTCORE_HPP
