#include "config.hpp"
#include "R_interface_common.hpp"

#include <cstddef> // size_t
#include <dbarts/cstdint.hpp>
#include <cmath> // fabs
#include <cstring> // memcpy

#include <external/alloca.h>
#include <external/random.h>

#include <Rmath.h> // unif_rand, norm_rand

#include <rc/bounds.h>
#include <rc/util.h>

#include <dbarts/bartFit.hpp>
#include <dbarts/control.hpp>
#include <dbarts/data.hpp>
#include <dbarts/model.hpp>
#include <dbarts/state.hpp>

#include "R_interface.hpp"

using std::size_t;
using std::uint32_t;

namespace dbarts {
  
  void deleteFit(BARTFit* fit) {
#ifdef THREAD_SAFE_UNLOAD
    Rprintf("deleting   %p\n", fit);
#endif
    if (fit == NULL) return;
    
    ext_rng_destroy(fit->control.rng);
    
    delete fit->model.treePrior;
    delete fit->model.muPrior;
    delete fit->model.sigmaSqPrior;
    
    delete [] fit->data.maxNumCuts;
    delete [] fit->data.variableTypes;
    
    delete fit;
  }
  
  void initializeControlFromExpression(Control& control, SEXP controlExpr)
  {
    int i_temp;
    
    SEXP slotExpr = Rf_getAttrib(controlExpr, Rf_install("binary"));
    control.responseIsBinary = rc_getBool(slotExpr, "binary response signifier", RC_LENGTH | RC_GEQ, rc_asRLength(1), RC_END);
    
    slotExpr = Rf_getAttrib(controlExpr, Rf_install("verbose"));
    control.verbose = rc_getBool(slotExpr, "verbose", RC_LENGTH | RC_GEQ, rc_asRLength(1), RC_END);
    
    slotExpr = Rf_getAttrib(controlExpr, Rf_install("keepTrainingFits"));
    control.keepTrainingFits = rc_getBool(slotExpr, "keep training fits", RC_LENGTH | RC_EQ, rc_asRLength(1), RC_END);
    
    slotExpr = Rf_getAttrib(controlExpr, Rf_install("useQuantiles"));
    control.useQuantiles = rc_getBool(slotExpr, "use quantiles", RC_LENGTH | RC_EQ, rc_asRLength(1), RC_END);
    
    slotExpr = Rf_getAttrib(controlExpr, Rf_install("n.samples"));
    i_temp = rc_getInt(slotExpr, "number of samples", RC_LENGTH | RC_EQ, rc_asRLength(1), RC_VALUE | RC_GEQ, 0, RC_END);
    control.numSamples = static_cast<size_t>(i_temp);
  
    slotExpr = Rf_getAttrib(controlExpr, Rf_install("n.burn"));
    i_temp = rc_getInt(slotExpr, "number of burn-in steps", RC_LENGTH | RC_EQ, rc_asRLength(1), RC_VALUE | RC_GEQ, 0, RC_END);
    control.numBurnIn = static_cast<size_t>(i_temp);
            
    slotExpr = Rf_getAttrib(controlExpr, Rf_install("n.trees"));
    i_temp = rc_getInt(slotExpr, "number of trees", RC_LENGTH | RC_EQ, rc_asRLength(1), RC_VALUE | RC_GEQ, 1, RC_END);
    control.numTrees = static_cast<size_t>(i_temp);
        
    slotExpr = Rf_getAttrib(controlExpr, Rf_install("n.threads"));
    i_temp = rc_getInt(slotExpr, "number of threads", RC_LENGTH | RC_EQ, rc_asRLength(1), RC_VALUE | RC_GEQ, 1, RC_END);
    control.numThreads = static_cast<size_t>(i_temp);
            
    slotExpr = Rf_getAttrib(controlExpr, Rf_install("n.thin"));
    i_temp = rc_getInt(slotExpr, "tree thinning rate", RC_LENGTH | RC_EQ, rc_asRLength(1), RC_VALUE | RC_GEQ, 0, RC_END);
    control.treeThinningRate = static_cast<uint32_t>(i_temp);
        
    slotExpr = Rf_getAttrib(controlExpr, Rf_install("printEvery"));
    i_temp = rc_getInt(slotExpr, "print every", RC_LENGTH | RC_EQ, rc_asRLength(1), RC_VALUE | RC_GEQ, 1, RC_NA | RC_YES, RC_END);
    if (i_temp != NA_INTEGER) control.printEvery = static_cast<uint32_t>(i_temp);
    
    slotExpr = Rf_getAttrib(controlExpr, Rf_install("printCutoffs"));
    i_temp = rc_getInt(slotExpr, "print cutoffs", RC_LENGTH | RC_EQ, rc_asRLength(1), RC_VALUE | RC_GEQ, 0, RC_NA | RC_YES, RC_END);
    if (i_temp == NA_INTEGER) i_temp = 0;
    control.printCutoffs = static_cast<uint32_t>(i_temp);
    
    
    if (control.rng == NULL) {
      ext_rng_userFunction uniformFunction;
      uniformFunction.f.stateless = &unif_rand;
      uniformFunction.state = NULL;
      control.rng = ext_rng_create(EXT_RNG_ALGORITHM_USER_UNIFORM, &uniformFunction);
      
      ext_rng_userFunction normalFunction;
      normalFunction.f.stateless = &norm_rand;
      normalFunction.state = NULL;
      ext_rng_setStandardNormalAlgorithm(control.rng, EXT_RNG_STANDARD_NORMAL_USER_NORM, &normalFunction);
    }
  }
  
  void invalidateControl(Control& control)
  {
    if (control.rng != NULL) { ext_rng_destroy(control.rng); control.rng = NULL; }
  }
    
  void initializeModelFromExpression(Model& model, SEXP modelExpr, const Control& control)
  {
    SEXP slotExpr = Rf_getAttrib(modelExpr, Rf_install("p.birth_death"));
    model.birthOrDeathProbability =
      rc_getDouble(slotExpr, "probability of birth/death rule", RC_LENGTH | RC_EQ, rc_asRLength(1),
                   RC_VALUE | RC_GEQ, 0.0, RC_VALUE | RC_LT, 1.0, RC_END);
    
    slotExpr = Rf_getAttrib(modelExpr, Rf_install("p.swap"));
    model.swapProbability =
      rc_getDouble(slotExpr, "probability of swap rule", RC_LENGTH | RC_EQ, rc_asRLength(1),
                   RC_VALUE | RC_GEQ, 0.0, RC_VALUE | RC_LT, 1.0, RC_END);
        
    slotExpr = Rf_getAttrib(modelExpr, Rf_install("p.change"));
    model.changeProbability =
      rc_getDouble(slotExpr, "probability of change rule", RC_LENGTH | RC_EQ, rc_asRLength(1),
                   RC_VALUE | RC_GEQ, 0.0, RC_VALUE | RC_LT, 1.0, RC_END);
        
    if (std::fabs(model.birthOrDeathProbability + model.swapProbability + model.changeProbability - 1.0) >= 1.0e-10)
      Rf_error("rule proposal probabilities must sum to 1.0");
    
    slotExpr = Rf_getAttrib(modelExpr, Rf_install("p.birth"));
    model.birthProbability =
      rc_getDouble(slotExpr, "probability of birth in birth/death rule", RC_LENGTH | RC_EQ, rc_asRLength(1),
                   RC_VALUE | RC_GT, 0.0, RC_VALUE | RC_LT, 1.0, RC_END);
    
    
    SEXP priorExpr = Rf_getAttrib(modelExpr, Rf_install("tree.prior"));
    CGMPrior* treePrior = new CGMPrior;
    model.treePrior = treePrior;
    
    slotExpr = Rf_getAttrib(priorExpr, Rf_install("power"));
    treePrior->power =
      rc_getDouble(slotExpr, "tree prior power", RC_LENGTH | RC_EQ, rc_asRLength(1),
                   RC_VALUE | RC_GT, 0.0, RC_END);
      
    slotExpr = Rf_getAttrib(priorExpr, Rf_install("base"));
    treePrior->base =
      rc_getDouble(slotExpr, "tree prior base", RC_LENGTH | RC_EQ, rc_asRLength(1),
                   RC_VALUE | RC_GT, 0.0, RC_VALUE | RC_LT, 1.0, RC_END);
    
    
    priorExpr = Rf_getAttrib(modelExpr, Rf_install("node.prior"));
      
    slotExpr = Rf_getAttrib(priorExpr, Rf_install("k"));
    double k = rc_getDouble(slotExpr, "k", RC_LENGTH | RC_EQ, rc_asRLength(1), RC_VALUE | RC_GT, 0.0, RC_END);
    model.muPrior = new NormalPrior(control, k);
    
    
    
    priorExpr = Rf_getAttrib(modelExpr, Rf_install("resid.prior"));
      
    slotExpr = Rf_getAttrib(priorExpr, Rf_install("df"));
    double sigmaPriorDf =
      rc_getDouble(slotExpr, "sigma prior degrees of freedom", RC_LENGTH | RC_EQ, rc_asRLength(1),
                   RC_VALUE | RC_GT, 0.0, RC_END);
    
    slotExpr = Rf_getAttrib(priorExpr, Rf_install("quantile"));
    double sigmaPriorQuantile =
      rc_getDouble(slotExpr, "sigma prior quantile", RC_LENGTH | RC_EQ, rc_asRLength(1),
                   RC_VALUE | RC_GT, 0.0, RC_VALUE | RC_LT, 1.0, RC_END);
    
    model.sigmaSqPrior = new ChiSquaredPrior(sigmaPriorDf, sigmaPriorQuantile);
  }
  
  void invalidateModel(Model& model)
  {
    delete static_cast<ChiSquaredPrior*>(model.sigmaSqPrior); model.sigmaSqPrior = NULL;
    delete static_cast<NormalPrior*>(model.muPrior);          model.muPrior = NULL;
    delete static_cast<CGMPrior*>(model.treePrior);           model.treePrior = NULL;
  }
  
  void initializeDataFromExpression(Data& data, SEXP dataExpr)
  {
    // SEXP dimsExpr;
    int* dims;
    
    SEXP slotExpr = Rf_getAttrib(dataExpr, Rf_install("y"));
    if (!Rf_isReal(slotExpr)) Rf_error("y must be of type real");
    if (rc_getLength(slotExpr) <= 0) Rf_error("length of y must be greater than 0");
    data.y = REAL(slotExpr);
    data.numObservations = rc_getLength(slotExpr);
    
    slotExpr = Rf_getAttrib(dataExpr, Rf_install("x"));
    if (!Rf_isReal(slotExpr)) Rf_error("x must be of type real");
    rc_assertDimConstraints(slotExpr, "dimensions of x", RC_LENGTH | RC_EQ, rc_asRLength(2), RC_VALUE | RC_EQ, static_cast<int>(data.numObservations), RC_END);
    dims = INTEGER(Rf_getAttrib(slotExpr, R_DimSymbol));
    
    // rc_assertIntConstraints(dimsExpr = Rf_getAttrib(slotExpr, R_DimSymbol), "dimensions of x", RC_LENGTH | RC_EQ, rc_asRLength(2), RC_END);
    // dims = INTEGER(dimsExpr);
    // if (static_cast<size_t>(dims[0]) != data.numObservations) Rf_error("number of rows of x and length of y must be equal");
    data.x = REAL(slotExpr);
    data.numPredictors = static_cast<size_t>(dims[1]);
    
    slotExpr = Rf_getAttrib(dataExpr, Rf_install("varTypes"));
    rc_assertIntConstraints(slotExpr, "variable types", RC_LENGTH | RC_EQ, rc_asRLength(data.numPredictors), RC_END);
    int* i_variableTypes = INTEGER(slotExpr);
    VariableType* variableTypes = new VariableType[data.numPredictors];
    for (size_t i = 0; i < data.numPredictors; ++i) variableTypes[i] = (i_variableTypes[i] == 0 ? ORDINAL : CATEGORICAL);
    data.variableTypes = variableTypes;
    
    slotExpr = Rf_getAttrib(dataExpr, Rf_install("x.test"));
    if (rc_isS4Null(slotExpr) || Rf_isNull(slotExpr) || rc_getLength(slotExpr) == 0) {
      data.x_test = NULL;
      data.numTestObservations = 0;
    } else {
      if (!Rf_isReal(slotExpr)) Rf_error ("x.test must be of type real");
      rc_assertDimConstraints(slotExpr, "dimensions of x.test", RC_LENGTH | RC_EQ, rc_asRLength(2),
                    RC_NA, RC_VALUE | RC_EQ, static_cast<int>(data.numPredictors), RC_END);
      dims = INTEGER(Rf_getAttrib(slotExpr, R_DimSymbol));
      
      // rc_assertIntConstraints(dimsExpr = Rf_getAttrib(slotExpr, R_DimSymbol), "x.test dimensions", RC_LENGTH | RC_EQ, rc_asRLength(2), RC_END);
      // dims = INTEGER(dimsExpr);
      // if (static_cast<size_t>(dims[1]) != data.numPredictors) Rf_error("number of columns of x.test and x must be equal");
      data.x_test = REAL(slotExpr);
      data.numTestObservations = static_cast<size_t>(dims[0]);
    }
    
    slotExpr = Rf_getAttrib(dataExpr, Rf_install("weights"));
    if (rc_isS4Null(slotExpr) || Rf_isNull(slotExpr) || rc_getLength(slotExpr) == 0) {
      data.weights = NULL;
    } else {
      rc_assertDoubleConstraints(slotExpr, "weights", RC_LENGTH | RC_EQ, rc_asRLength(data.numObservations), RC_END);
      data.weights = REAL(slotExpr);
    }
    
    slotExpr = Rf_getAttrib(dataExpr, Rf_install("offset"));
    if (rc_isS4Null(slotExpr) || Rf_isNull(slotExpr) || rc_getLength(slotExpr) == 0) {
      data.offset = NULL;
    } else {
      rc_assertDoubleConstraints(slotExpr, "offset", RC_LENGTH | RC_EQ, rc_asRLength(data.numObservations), RC_END);
      data.offset = REAL(slotExpr);
    }
    
    slotExpr = Rf_getAttrib(dataExpr, Rf_install("offset.test"));
    if (rc_isS4Null(slotExpr) || Rf_isNull(slotExpr) || rc_getLength(slotExpr) == 0) {
      data.testOffset = NULL;
    } else {
      rc_assertDoubleConstraints(slotExpr, "test offset", RC_LENGTH | RC_EQ, rc_asRLength(data.numTestObservations), RC_END);
      data.testOffset = REAL(slotExpr);
    }
    
    slotExpr = Rf_getAttrib(dataExpr, Rf_install("sigma"));
    data.sigmaEstimate = rc_getDouble(slotExpr, "sigma estimate", RC_LENGTH | RC_EQ, rc_asRLength(1), RC_VALUE | RC_GT, 0.0, RC_END);
    
    slotExpr = Rf_getAttrib(dataExpr, Rf_install("n.cuts"));
    rc_assertIntConstraints(slotExpr, "maximum number of cuts", RC_LENGTH | RC_EQ, rc_asRLength(data.numPredictors), RC_END);
    int* i_maxNumCuts = INTEGER(slotExpr);
    uint32_t* maxNumCuts = new uint32_t[data.numPredictors];
    for (size_t i = 0; i < data.numPredictors; ++i) maxNumCuts[i] = static_cast<uint32_t>(i_maxNumCuts[i]);
    data.maxNumCuts = maxNumCuts;
  }
  
  void invalidateData(Data& data)
  {
    delete [] data.maxNumCuts;    data.maxNumCuts = NULL;
    delete [] data.variableTypes; data.variableTypes = NULL;
  }
  
  SEXP createStateExpressionFromFit(const BARTFit& fit)
  {
    const Control& control(fit.control);
    const Data& data(fit.data);
    const State& state(fit.state);
    
    SEXP result = PROTECT(R_do_new_object(R_do_MAKE_CLASS("dbartsState")));
    
    SEXP slotExpr = rc_allocateInSlot(result, Rf_install("fit.tree"), REALSXP, static_cast<R_xlen_t>(data.numObservations * control.numTrees));
    rc_setDims(slotExpr, static_cast<int>(data.numObservations), static_cast<int>(control.numTrees), -1);
    std::memcpy(REAL(slotExpr), state.treeFits, data.numObservations * control.numTrees * sizeof(double));
    
    slotExpr = rc_allocateInSlot(result, Rf_install("fit.total"), REALSXP, static_cast<R_xlen_t>(data.numObservations));
    std::memcpy(REAL(slotExpr), state.totalFits, data.numObservations * sizeof(double));
    
    if (data.numTestObservations == 0) {
      R_do_slot_assign(result, Rf_install("fit.test"), R_NilValue);
    } else {
      slotExpr = rc_allocateInSlot(result, Rf_install("fit.test"), REALSXP, static_cast<R_xlen_t>(data.numTestObservations));
      std::memcpy(REAL(slotExpr), state.totalTestFits, data.numTestObservations * sizeof(double));
    }
    
    slotExpr = rc_allocateInSlot(result, Rf_install("sigma"), REALSXP, 1);
    REAL(slotExpr)[0] = state.sigma;
    
    slotExpr = rc_allocateInSlot(result, Rf_install("runningTime"), REALSXP, 1);
    REAL(slotExpr)[0] = state.runningTime;
    
    slotExpr = rc_allocateInSlot(result, Rf_install("trees"), STRSXP, static_cast<R_xlen_t>(control.numTrees));
  
    const char** treeStrings = const_cast<const char**>(state.createTreeStrings(fit));
    for (size_t i = 0; i < control.numTrees; ++i) {
      SET_STRING_ELT(slotExpr, static_cast<R_xlen_t>(i), Rf_mkChar(treeStrings[i]));
      delete [] treeStrings[i];
    }
    delete [] treeStrings;
    
    return result;
  }
  
  void storeStateExpressionFromFit(const BARTFit& fit, SEXP stateExpr)
  {
    const Control& control(fit.control);
    const Data& data(fit.data);
    const State& state(fit.state);
    
    SEXP slotExpr = Rf_getAttrib(stateExpr, Rf_install("fit.tree"));
    SEXP dimsExpr = Rf_getAttrib(slotExpr, R_DimSymbol);
    if (rc_getLength(dimsExpr) != 2) Rf_error("dimensions of state@fit.tree indicate that it is not a matrix");
    int* dims = INTEGER(dimsExpr);
    if (static_cast<size_t>(dims[0]) != data.numObservations || static_cast<size_t>(dims[1]) != control.numTrees) {
      // Rf_error("dimensions of state@fit.tree do not match object");
      slotExpr = rc_allocateInSlot(stateExpr, Rf_install("fit.tree"), REALSXP, rc_asRLength(data.numObservations * control.numTrees));
      rc_setDims(slotExpr, static_cast<int>(data.numObservations), static_cast<int>(control.numTrees), -1);
      
      rc_allocateInSlot(stateExpr, Rf_install("fit.total"), REALSXP, rc_asRLength(data.numObservations));
      rc_allocateInSlot(stateExpr, Rf_install("trees"), STRSXP, rc_asRLength(control.numTrees));
    }
    std::memcpy(REAL(slotExpr), state.treeFits, data.numObservations * control.numTrees * sizeof(double));
    
    slotExpr = Rf_getAttrib(stateExpr, Rf_install("fit.total"));
    if (rc_getLength(slotExpr) != data.numObservations) Rf_error("length of state@fit.total does not match object");
    std::memcpy(REAL(slotExpr), state.totalFits, data.numObservations * sizeof(double));
    
    slotExpr = Rf_getAttrib(stateExpr, Rf_install("fit.test"));
    if (data.numTestObservations != 0) {
      // Rf_error("length of state@fit.test does not match object");
      if (!Rf_isNull(slotExpr) && rc_getLength(slotExpr) != data.numTestObservations)
        slotExpr = rc_allocateInSlot(stateExpr, Rf_install("fit.test"), REALSXP, rc_asRLength(data.numTestObservations));
      std::memcpy(REAL(slotExpr), state.totalTestFits, data.numTestObservations * sizeof(double));
    } else {
      R_do_slot_assign(stateExpr, Rf_install("fit.test"), R_NilValue);
    }
    
    slotExpr = Rf_getAttrib(stateExpr, Rf_install("sigma"));
    if (rc_getLength(slotExpr) != 1) Rf_error("length of state@sigma does not match object");
    REAL(slotExpr)[0] = state.sigma;
    
    slotExpr = Rf_getAttrib(stateExpr, Rf_install("runningTime"));
    if (rc_getLength(slotExpr) != 1) Rf_error("length of state@runningTime does not match object");
    REAL(slotExpr)[0] = state.runningTime;
    
    slotExpr = Rf_getAttrib(stateExpr, Rf_install("trees"));
    // if (rc_getLength(slotExpr) != control.numTrees) Rf_error("length of state@trees does not match object");
    
    const char** treeStrings = const_cast<const char**>(state.createTreeStrings(fit));
    for (size_t i = 0; i < control.numTrees; ++i) {
      SET_STRING_ELT(slotExpr, static_cast<R_xlen_t>(i), Rf_mkChar(treeStrings[i]));
      delete [] treeStrings[i];
    }
    delete [] treeStrings;
  }
  
  void initializeStateFromExpression(const BARTFit& fit, State& state, SEXP stateExpr)
  {
    const Control& control(fit.control);
    const Data& data(fit.data);
    
    SEXP slotExpr = Rf_getAttrib(stateExpr, Rf_install("fit.tree"));
    std::memcpy(state.treeFits, const_cast<const double*>(REAL(slotExpr)), data.numObservations * control.numTrees * sizeof(double));
    
    slotExpr = Rf_getAttrib(stateExpr, Rf_install("fit.total"));
    std::memcpy(state.totalFits, const_cast<const double*>(REAL(slotExpr)), data.numObservations * sizeof(double));
    
    if (data.numTestObservations != 0) {
      slotExpr = Rf_getAttrib(stateExpr, Rf_install("fit.test"));
      std::memcpy(state.totalTestFits, const_cast<const double*>(REAL(slotExpr)), data.numTestObservations * sizeof(double));
    }
    
    slotExpr = Rf_getAttrib(stateExpr, Rf_install("sigma"));
    state.sigma = REAL(slotExpr)[0];
    
    slotExpr = Rf_getAttrib(stateExpr, Rf_install("runningTime"));
    state.runningTime = REAL(slotExpr)[0];
    
    slotExpr = Rf_getAttrib(stateExpr, Rf_install("trees"));
    const char** treeStrings = ext_stackAllocate(control.numTrees, const char*);
    for (size_t i = 0; i < control.numTrees; ++i) {
      treeStrings[i] = CHAR(STRING_ELT(slotExpr, rc_asRLength(i)));
    }
    state.recreateTreesFromStrings(fit, treeStrings);
    
    ext_stackFree(treeStrings);
  }
}

