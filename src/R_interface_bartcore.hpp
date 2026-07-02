#ifndef R_INTERFACE_BARTCORE_HPP
#define R_INTERFACE_BARTCORE_HPP

// internal bridge to the generalized core (bartcore/); see
// docs/design/core-generalization.md, phase 2

#include <external/Rinternals.h> // SEXP

extern "C" {

SEXP bartcore_create(SEXP control, SEXP model, SEXP data);
SEXP bartcore_run(SEXP ptr, SEXP numBurnIn, SEXP numSamples);
SEXP bartcore_setOffset(SEXP ptr, SEXP offset, SEXP updateScale);
SEXP bartcore_setResponse(SEXP ptr, SEXP y);
SEXP bartcore_setSigma(SEXP ptr, SEXP sigma);
SEXP bartcore_setTestPredictor(SEXP ptr, SEXP x_test);
SEXP bartcore_getLatents(SEXP ptr);

}

#endif // R_INTERFACE_BARTCORE_HPP
