#include "config.hpp"
#include "R_interface_common.hpp"

#include <cstddef> // size_t
#include <dbarts/cstdint.hpp>
#include <cmath> // fabs
#include <cstring> // memcpy

#include <external/alloca.h>
#include <external/random.h>
#include <external/string.h>

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

namespace {
  const char* const rngNames[] = {
    "Wichmann-Hill",
    "Marsaglia-Multicarry",
    "Super-Duper",
    "Mersenne-Twister",
    "Knuth-TAOCP",
    "user-supplied",
    "Knuth-TAOCP-2002",
    "L'Ecuyer-CMRG",
    "default"
  };
  
  const char* const rngNormalNames[] = {
    "Buggy Kinderman-Ramage",
    "Ahrens-Dieter",
    "Box-Muller",
    "user-supplied",
    "Inversion",
    "Kinderman-Ramage",
    "default"
  };
  
  const char* const runModes[] = {
    "sequentialUpdates",
    "fixedSamples"
  };
}

namespace dbarts {
  
  void deleteFit(BARTFit* fit) {
#ifdef THREAD_SAFE_UNLOAD
    Rprintf("deleting   %p\n", fit);
#endif
    if (fit == NULL) return;
    
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
    
    slotExpr = Rf_getAttrib(controlExpr, Rf_install("runMode"));
    if (rc_getLength(slotExpr) != 1) Rf_error("slot 'runMod' must be of length 1");
    const char* runModeName = CHAR(STRING_ELT(slotExpr, 0));
    
    size_t runModeNumber;
    int errorCode = ext_str_matchInArray(runModeName, runModes, static_cast<size_t>(FIXED_SAMPLES - SEQUENTIAL_UPDATES + 1), &runModeNumber);
    if (errorCode != 0) Rf_error("error matching run mode string: %s", std::strerror(errorCode));
    if (runModeNumber == EXT_STR_NO_MATCH) Rf_error("unsupported run mode '%s'", runModeName);
    
    control.runMode = static_cast<RunMode>(runModeNumber);
    
    slotExpr = Rf_getAttrib(controlExpr, Rf_install("n.samples"));
    i_temp = rc_getInt(slotExpr, "number of samples", RC_LENGTH | RC_EQ, rc_asRLength(1), RC_VALUE | RC_GEQ, 0, RC_END);
    control.defaultNumSamples = static_cast<size_t>(i_temp);
  
    slotExpr = Rf_getAttrib(controlExpr, Rf_install("n.burn"));
    i_temp = rc_getInt(slotExpr, "number of burn-in steps", RC_LENGTH | RC_EQ, rc_asRLength(1), RC_VALUE | RC_GEQ, 0, RC_END);
    control.defaultNumBurnIn = static_cast<size_t>(i_temp);
            
    slotExpr = Rf_getAttrib(controlExpr, Rf_install("n.trees"));
    i_temp = rc_getInt(slotExpr, "number of trees", RC_LENGTH | RC_EQ, rc_asRLength(1), RC_VALUE | RC_GEQ, 1, RC_END);
    control.numTrees = static_cast<size_t>(i_temp);
    
    slotExpr = Rf_getAttrib(controlExpr, Rf_install("n.chains"));
    i_temp = rc_getInt(slotExpr, "number of chains", RC_LENGTH | RC_EQ, rc_asRLength(1), RC_VALUE | RC_GEQ, 1, RC_END);
    control.numChains = static_cast<size_t>(i_temp);
    
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
    

    slotExpr = Rf_getAttrib(controlExpr, Rf_install("rngKind"));
    size_t slotLength = rc_getLength(slotExpr);
    if (slotLength != 1) Rf_error("slot 'rngKind' must be of length 1");
    const char* rngKindName = CHAR(STRING_ELT(slotExpr, 0));
    
    size_t rngKindNumber;
    errorCode = ext_str_matchInArray(rngKindName, rngNames, static_cast<size_t>(EXT_RNG_ALGORITHM_INVALID - EXT_RNG_ALGORITHM_WICHMANN_HILL + 1), &rngKindNumber);
    if (errorCode != 0) Rf_error("error matching rng kind string: %s", std::strerror(errorCode));
    if (rngKindNumber == EXT_STR_NO_MATCH) Rf_error("unsupported rng kind '%s'", rngKindName);
    
    control.rng_algorithm = static_cast<ext_rng_algorithm_t>(rngKindNumber);
    
    
    slotExpr = Rf_getAttrib(controlExpr, Rf_install("rngNormalKind"));
    slotLength = rc_getLength(slotExpr);
    if (slotLength != 1) Rf_error("slot 'rngNormalKind' must be of length 1");
    const char* rngNormalKindName = CHAR(STRING_ELT(slotExpr, 0));
    
    size_t rngNormalKindNumber;
    errorCode = ext_str_matchInArray(rngNormalKindName, rngNormalNames, static_cast<size_t>(EXT_RNG_STANDARD_NORMAL_INVALID - EXT_RNG_STANDARD_NORMAL_BUGGY_KINDERMAN_RAMAGE + 1), &rngNormalKindNumber);
    if (errorCode != 0) Rf_error("error matching rng normal kind string: %s", std::strerror(errorCode));
    if (rngNormalKindNumber == EXT_STR_NO_MATCH) Rf_error("unsupported rng normal kind '%s'", rngNormalKindName);
    
    control.rng_standardNormal = static_cast<ext_rng_standardNormal_t>(rngNormalKindNumber);
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
      if (!Rf_isReal(slotExpr)) Rf_error("x.test must be of type real");
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
    data.sigmaEstimate = rc_getDouble(slotExpr, "sigma estimate", RC_LENGTH | RC_EQ, rc_asRLength(1), RC_NA | RC_YES, RC_VALUE | RC_GT, 0.0, RC_END);
    
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
    const State* state(fit.state);
    
    SEXP result = PROTECT(rc_newList(control.numChains));
    size_t numSamples = control.runMode == FIXED_SAMPLES ? fit.currentNumSamples : 1;
    
    SEXP slotExpr = rc_allocateInSlot(result, Rf_install("currentNumSamples"), INTSXP, 1);
    INTEGER(slotExpr)[0] = static_cast<int>(fit.currentNumSamples);
    
    for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum) {
      SEXP result_i = PROTECT(R_do_new_object(R_do_MAKE_CLASS("dbartsState")));
      SET_VECTOR_ELT(result, chainNum, result_i);
      
      slotExpr = rc_allocateInSlot(result_i, Rf_install("fit.tree"), REALSXP, static_cast<R_xlen_t>(data.numObservations * control.numTrees * numSamples));
      rc_setDims(slotExpr, static_cast<int>(data.numObservations), static_cast<int>(control.numTrees), static_cast<int>(numSamples), -1);
      std::memcpy(REAL(slotExpr), state[chainNum].treeFits, data.numObservations * control.numTrees * numSamples * sizeof(double));
      
      slotExpr = rc_allocateInSlot(result_i, Rf_install("sigma"), REALSXP, numSamples);
      std::memcpy(REAL(slotExpr), state[chainNum].sigma, numSamples * sizeof(double));
      
      slotExpr = rc_allocateInSlot(result_i, Rf_install("trees"), STRSXP, static_cast<R_xlen_t>(control.numTrees * numSamples));
      rc_setDims(slotExpr, static_cast<int>(control.numTrees), static_cast<int>(numSamples), -1);
    
      const char** treeStrings = const_cast<const char**>(state[chainNum].createTreeStrings(fit));
      for (size_t sampleNum = 0; sampleNum < numSamples; ++sampleNum) {
        for (size_t treeNum = 0; treeNum < control.numTrees; ++treeNum) {
          size_t treeOffset = treeNum + sampleNum * control.numTrees;
          SET_STRING_ELT(slotExpr, static_cast<R_xlen_t>(treeOffset), Rf_mkChar(treeStrings[treeOffset]));
          delete [] treeStrings[treeOffset];
        }
      }
      delete [] treeStrings;
      
      size_t rngStateLength = ext_rng_getSerializedStateLength(state[chainNum].rng) / sizeof(int);
      slotExpr = rc_allocateInSlot(result_i, Rf_install("rng.state"), INTSXP, rc_asRLength(rngStateLength));
      ext_rng_writeSerializedState(state[chainNum].rng, INTEGER(slotExpr));
    }
    
    slotExpr = rc_allocateInSlot(result, Rf_install("runningTime"), REALSXP, 1);
    REAL(slotExpr)[0] = fit.runningTime;
    
    UNPROTECT(1 + control.numChains);
    
    return result;
  }
  
  void storeStateExpressionFromFit(const BARTFit& fit, SEXP stateExpr)
  {
    const Control& control(fit.control);
    const Data& data(fit.data);
    const State* state(fit.state);
    
    SEXP fitTreeSym  = Rf_install("fit.tree");
    SEXP treeSym     = Rf_install("trees");
    SEXP sigmaSym    = Rf_install("sigma");
    SEXP rngStateSym = Rf_install("rng.state");
    
    // check to see if it is an old-style saved object with only a single state
    SEXP classExpr = rc_getClass(stateExpr);
    if (!Rf_isNull(classExpr) && std::strcmp(CHAR(STRING_ELT(classExpr, 0)), "dbartsState") == 0) 
      Rf_error("object from earlier version detected - model must be refit");
    
    if (rc_getLength(stateExpr) != control.numChains)
      Rf_error("length of state list not equal to number of chains");
    
    size_t numSamples = control.runMode == FIXED_SAMPLES ? fit.currentNumSamples : 1;
    
    SEXP slotExpr = Rf_getAttrib(stateExpr, Rf_install("currentNumSamples"));
    INTEGER(slotExpr)[0] = static_cast<int>(fit.currentNumSamples);
    
    for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum)
    {
      SEXP stateExpr_i = VECTOR_ELT(stateExpr, static_cast<R_xlen_t>(chainNum));
      classExpr = rc_getClass(stateExpr_i);
      if (std::strcmp(CHAR(STRING_ELT(classExpr, 0)), "dbartsState") != 0) Rf_error("'state' not of class 'dbartsState'");
      
      slotExpr = Rf_getAttrib(stateExpr_i, fitTreeSym);
      SEXP dimsExpr = rc_getDims(slotExpr);
      if (rc_getLength(dimsExpr) != 3) Rf_error("dimensions of state@fit.tree indicate that it is not a array");
      int* dims = INTEGER(dimsExpr);
      if (static_cast<size_t>(dims[0]) != data.numObservations ||
          static_cast<size_t>(dims[1]) != control.numTrees ||
          static_cast<size_t>(dims[2]) != numSamples)
      {
        // Rf_error("dimensions of state@fit.tree do not match object");
        slotExpr = rc_allocateInSlot(stateExpr_i, fitTreeSym, REALSXP, rc_asRLength(data.numObservations * control.numTrees * numSamples));
        rc_setDims(slotExpr, static_cast<int>(data.numObservations), static_cast<int>(control.numTrees), static_cast<int>(numSamples), -1);
        
        // rc_allocateInSlot(stateExpr_i, fitTotalSym, REALSXP, rc_asRLength(data.numObservations));
        SEXP treeSlot = rc_allocateInSlot(stateExpr_i, treeSym, STRSXP, rc_asRLength(control.numTrees * numSamples));
        rc_setDims(treeSlot, static_cast<int>(control.numTrees), static_cast<int>(numSamples), -1);
      }
      std::memcpy(REAL(slotExpr), state[chainNum].treeFits, data.numObservations * control.numTrees * numSamples * sizeof(double));
            
      slotExpr = Rf_getAttrib(stateExpr_i, sigmaSym);
      if (rc_getLength(slotExpr) != numSamples) Rf_error("length of state@sigma does not match object");
      std::memcpy(REAL(slotExpr), state[chainNum].sigma, numSamples * sizeof(double));
      
      slotExpr = Rf_getAttrib(stateExpr_i, treeSym);
      
      const char** treeStrings = const_cast<const char**>(state[chainNum].createTreeStrings(fit));
      for (size_t sampleNum = 0; sampleNum < numSamples; ++sampleNum) {
        for (size_t treeNum = 0; treeNum < control.numTrees; ++treeNum) {
          size_t treeOffset = treeNum + sampleNum * control.numTrees;
          SET_STRING_ELT(slotExpr, static_cast<R_xlen_t>(treeOffset), Rf_mkChar(treeStrings[treeOffset]));
          delete [] treeStrings[treeOffset];
        }
      }
      delete [] treeStrings;
      
      size_t rngStateLength = ext_rng_getSerializedStateLength(state[chainNum].rng) / sizeof(int);
      slotExpr = Rf_getAttrib(stateExpr_i, rngStateSym);
      if (rc_getLength(slotExpr) != rngStateLength)
        slotExpr = rc_allocateInSlot(stateExpr_i, rngStateSym, INTSXP, rc_asRLength(rngStateLength));
      ext_rng_writeSerializedState(state[chainNum].rng, INTEGER(slotExpr));
    }
    
    slotExpr = Rf_getAttrib(stateExpr, Rf_install("runningTime"));
    REAL(slotExpr)[0] = fit.runningTime;
  }
  
  void initializeStateFromExpression(BARTFit& fit, State* state, SEXP stateExpr)
  {
    const Control& control(fit.control);
    const Data& data(fit.data);
    
    // check to see if it is an old-style saved object with only a single state
    SEXP classExpr = rc_getClass(stateExpr);
    if (!Rf_isNull(classExpr) && std::strcmp(CHAR(STRING_ELT(classExpr, 0)), "dbartsState") == 0) 
      Rf_error("object from earlier version detected - model must be refit");
    
    size_t currentNumSamples = static_cast<size_t>(INTEGER(Rf_getAttrib(stateExpr, Rf_install("currentNumSamples")))[0]);
    if (fit.currentNumSamples != currentNumSamples && control.runMode == FIXED_SAMPLES)
      for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum)
        state[chainNum].resize(fit, currentNumSamples);
    fit.currentNumSamples = currentNumSamples;
    
    size_t numSamples = control.runMode == FIXED_SAMPLES ? fit.currentNumSamples : 1;
    
    for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum)
    {
      SEXP stateExpr_i = VECTOR_ELT(stateExpr, static_cast<R_xlen_t>(chainNum));
      classExpr = rc_getClass(stateExpr_i);
      if (std::strcmp(CHAR(STRING_ELT(classExpr, 0)), "dbartsState") != 0) Rf_error("'state' not of class 'dbartsState'");
      
      SEXP slotExpr = Rf_getAttrib(stateExpr_i, Rf_install("fit.tree"));
      std::memcpy(state[chainNum].treeFits, const_cast<const double*>(REAL(slotExpr)), data.numObservations * control.numTrees * numSamples * sizeof(double));
            
      slotExpr = Rf_getAttrib(stateExpr_i, Rf_install("sigma"));
      std::memcpy(state[chainNum].sigma, const_cast<const double*>(REAL(slotExpr)), numSamples * sizeof(double));
      
      slotExpr = Rf_getAttrib(stateExpr_i, Rf_install("trees"));
      const char** treeStrings = ext_stackAllocate(control.numTrees * numSamples, const char*);
      for (size_t treeNum = 0; treeNum < control.numTrees * numSamples; ++treeNum) {
        treeStrings[treeNum] = CHAR(STRING_ELT(slotExpr, rc_asRLength(static_cast<R_xlen_t>(treeNum))));
      }
      state[chainNum].recreateTreesFromStrings(fit, treeStrings);
      
      ext_stackFree(treeStrings);
      
      ext_rng_readSerializedState(state[chainNum].rng, INTEGER(Rf_getAttrib(stateExpr_i, Rf_install("rng.state"))));
    }
    
    fit.rebuildScratchFromState();
      
    SEXP slotExpr = Rf_getAttrib(stateExpr, Rf_install("runningTime"));
    fit.runningTime = REAL(slotExpr)[0];
  }
}

