#ifndef R_INTERFACE_COMMON_HPP
#define R_INTERFACE_COMMON_HPP

// struct SEXPREC;
// typedef struct SEXPREC* SEXP;
#include <external/Rinternals.h> // SEXP

namespace dbarts {
  struct BARTFit;
  struct Control;
  struct Data;
  struct Model;
  struct State;
  
  void deleteFit(BARTFit* fit);
  
  void initializeControlFromExpression(Control& control, SEXP controlExpr);
  void initializeModelFromExpression(Model& model, SEXP modelExpr, const Control& control);
  void initializeDataFromExpression(Data& data, SEXP dataExpr);
  
  void invalidateControl(Control& control);
  void invalidateModel(Model& model);
  void invalidateData(Data& data);
    
  void initializeStateFromExpression(const BARTFit& fit, State& state, SEXP stateExpr);
  SEXP createStateExpressionFromFit(const BARTFit& fit); // result has a protect count of 1
  void storeStateExpressionFromFit(const BARTFit& fit, SEXP state);
}

#endif

