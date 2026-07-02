#ifndef R_INTERFACE_BARTCORE_HPP
#define R_INTERFACE_BARTCORE_HPP

// internal bridge to the generalized core (bartcore/); see
// docs/design/core-generalization.md, phase 2

#include <external/Rinternals.h> // SEXP

extern "C" {

SEXP bartcore_create(SEXP control, SEXP model, SEXP data, SEXP family);
SEXP bartcore_run(SEXP ptr, SEXP numBurnIn, SEXP numSamples);
SEXP bartcore_setOffset(SEXP ptr, SEXP offset, SEXP updateScale);
SEXP bartcore_setResponse(SEXP ptr, SEXP y);
SEXP bartcore_setSigma(SEXP ptr, SEXP sigma);
SEXP bartcore_setData(SEXP ptr, SEXP data);
SEXP bartcore_setTestPredictor(SEXP ptr, SEXP x_test);
SEXP bartcore_setPredictor(SEXP ptr, SEXP x, SEXP forceUpdate,
                           SEXP updateCutPoints);
SEXP bartcore_updatePredictor(SEXP ptr, SEXP x, SEXP columns, SEXP forceUpdate,
                              SEXP updateCutPoints);
SEXP bartcore_updatePredictorPerObservation(SEXP ptr, SEXP x, SEXP column);
SEXP bartcore_updatePredictorPerObservationJointly(SEXP ptrs, SEXP x,
                                                   SEXP columns);
SEXP bartcore_setCutPoints(SEXP ptr, SEXP cutPoints, SEXP columns);
SEXP bartcore_getSigmas(SEXP ptr);
SEXP bartcore_isValidPointer(SEXP ptr);
SEXP bartcore_getLatents(SEXP ptr);

}

#endif // R_INTERFACE_BARTCORE_HPP
