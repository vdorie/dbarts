#include "config.hpp"
#include "R_interface_common.hpp"

#include <cstddef> // size_t
#include <dbarts/cstdint.hpp>
#include <cmath> // fabs
#include <cstring> // memcpy

#include <misc/alloca.h>
#include <misc/string.h>

#include <external/random.h>

#include <Rmath.h> // unif_rand, norm_rand

#include <rc/bounds.h>
#include <rc/util.h>

#include <dbarts/bartFit.hpp>
#include <dbarts/control.hpp>
#include <dbarts/data.hpp>
#include <dbarts/model.hpp>
#include <dbarts/state.hpp>

// #include "R_interface.hpp"

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
}

namespace dbarts {
  
  void deleteFit(BARTFit* fit) {
#ifdef THREAD_SAFE_UNLOAD
    Rprintf("deleting   %p\n", fit);
#endif
    if (fit == NULL) return;
    
    delete fit->model.kPrior;
    delete fit->model.sigmaSqPrior;
    delete fit->model.muPrior;
    delete fit->model.treePrior;
    
    delete [] fit->data.variableTypes;
    delete [] fit->data.maxNumCuts;
    
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
    
    slotExpr = Rf_getAttrib(controlExpr, Rf_install("keepTrees"));
    if (rc_getLength(slotExpr) != 1) Rf_error("slot 'keepTrees' must be of length 1");
    control.keepTrees = rc_getBool(slotExpr, "keep trees", RC_LENGTH | RC_EQ, rc_asRLength(1), RC_END);
    
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
    int errorCode = misc_str_matchInArray(rngKindName, rngNames, static_cast<size_t>(EXT_RNG_ALGORITHM_INVALID - EXT_RNG_ALGORITHM_WICHMANN_HILL + 1), &rngKindNumber);
    if (errorCode != 0) Rf_error("error matching rng kind string: %s", std::strerror(errorCode));
    if (rngKindNumber == MISC_STR_NO_MATCH) Rf_error("unsupported rng kind '%s'", rngKindName);
    
    control.rng_algorithm = static_cast<rng_algorithm_t>(rngKindNumber);
    
    
    slotExpr = Rf_getAttrib(controlExpr, Rf_install("rngNormalKind"));
    slotLength = rc_getLength(slotExpr);
    if (slotLength != 1) Rf_error("slot 'rngNormalKind' must be of length 1");
    const char* rngNormalKindName = CHAR(STRING_ELT(slotExpr, 0));
    
    size_t rngNormalKindNumber;
    errorCode = misc_str_matchInArray(rngNormalKindName, rngNormalNames, static_cast<size_t>(EXT_RNG_STANDARD_NORMAL_INVALID - EXT_RNG_STANDARD_NORMAL_BUGGY_KINDERMAN_RAMAGE + 1), &rngNormalKindNumber);
    if (errorCode != 0) Rf_error("error matching rng normal kind string: %s", std::strerror(errorCode));
    if (rngNormalKindNumber == MISC_STR_NO_MATCH) Rf_error("unsupported rng normal kind '%s'", rngNormalKindName);
    
    control.rng_standardNormal = static_cast<rng_standardNormal_t>(rngNormalKindNumber);
  }
}

namespace {
  struct ModelStackDeconstructor {
    dbarts::TreePrior* treePrior;
    dbarts::EndNodePrior* muPrior;
    dbarts::ResidualVariancePrior* sigmaSqPrior;
    dbarts::EndNodeHyperprior* kPrior;
    
    ModelStackDeconstructor() : treePrior(NULL), muPrior(NULL), sigmaSqPrior(NULL), kPrior(NULL) { }
    ~ModelStackDeconstructor() {
      delete kPrior;
      delete sigmaSqPrior;
      delete muPrior;
      delete treePrior;
    }
  };
}

namespace dbarts {
  void initializeModelFromExpression(Model& model, SEXP modelExpr, const Control& control)
  {
    int errorCode;
    size_t priorType;
    
    ModelStackDeconstructor stackModel;
    
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
    
    slotExpr = Rf_getAttrib(modelExpr, Rf_install("node.scale"));
    model.nodeScale = 
      rc_getDouble(slotExpr, "scale of node prior", RC_LENGTH | RC_EQ, rc_asRLength(1),
                   RC_VALUE | RC_GT, 0.0, RC_END);
    
    
    SEXP priorExpr = Rf_getAttrib(modelExpr, Rf_install("tree.prior"));
    CGMPrior* treePrior = new CGMPrior;
    stackModel.treePrior = treePrior;
    
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
    if (Rf_isReal(slotExpr)) {
      double k = rc_getDouble(slotExpr, "k", RC_LENGTH | RC_EQ, rc_asRLength(1), RC_VALUE | RC_GT, 0.0, RC_END);
      stackModel.muPrior = new NormalPrior(control, model, k);
    } else {
      SEXP classExpr = rc_getClass(slotExpr);
      const char* classStr = CHAR(STRING_ELT(classExpr, 0));
      errorCode = misc_str_matchInVArray(classStr, &priorType, "dbartsChiHyperprior", NULL);
      if (errorCode != 0) Rf_error("error matching node hyper prior: %s", std::strerror(errorCode));
      if (priorType == MISC_STR_NO_MATCH) Rf_error("unsupported hyperprior type '%s'", classStr);
      // some day, we will have a switch statement here
      double degreesOfFreedom =
        rc_getDouble(Rf_getAttrib(slotExpr, Rf_install("degreesOfFreedom")), "degreesOfFreedom",
                     RC_LENGTH | RC_EQ, rc_asRLength(1),
                     RC_VALUE | RC_GT, 0.0, RC_END);
      double scale =
        rc_getDouble(Rf_getAttrib(slotExpr, Rf_install("scale")), "scale",
                     RC_LENGTH | RC_EQ, rc_asRLength(1),
                     RC_VALUE | RC_GT, 0.0, RC_END);
      
      stackModel.muPrior = new NormalPrior(control, model, 2.0);
      stackModel.kPrior = new ChiHyperprior(degreesOfFreedom, scale);
    }
    
    priorExpr = Rf_getAttrib(modelExpr, Rf_install("resid.prior"));
    SEXP classExpr = rc_getClass(priorExpr);
    const char* classStr = CHAR(STRING_ELT(classExpr, 0));
    errorCode = misc_str_matchInVArray(classStr, &priorType, "dbartsChiSqPrior", "dbartsFixedPrior", NULL);
    if (errorCode != 0) Rf_error("error matching residual variance prior: %s", std::strerror(errorCode));
    if (priorType == MISC_STR_NO_MATCH) Rf_error("unsupported residual variance prior type '%s'", classStr);
    
    switch (priorType) {
      case 0:
      {
        slotExpr = Rf_getAttrib(priorExpr, Rf_install("df"));
        double sigmaPriorDf =
          rc_getDouble(slotExpr, "sigma prior degrees of freedom", RC_LENGTH | RC_EQ, rc_asRLength(1),
                     RC_VALUE | RC_GT, 0.0, RC_END);
        
        slotExpr = Rf_getAttrib(priorExpr, Rf_install("quantile"));
        double sigmaPriorQuantile =
          rc_getDouble(slotExpr, "sigma prior quantile", RC_LENGTH | RC_EQ, rc_asRLength(1),
                     RC_VALUE | RC_GT, 0.0, RC_VALUE | RC_LT, 1.0, RC_END);
        
        stackModel.sigmaSqPrior = new ChiSquaredPrior(sigmaPriorDf, sigmaPriorQuantile);
      }
      break;
      default:
      slotExpr = Rf_getAttrib(priorExpr, Rf_install("value"));
      stackModel.sigmaSqPrior = new FixedPrior(
          rc_getDouble(slotExpr, "residual variance prior fixed value", RC_LENGTH | RC_EQ, rc_asRLength(1),
                       RC_VALUE | RC_GT, 0.0, RC_END));
      break;
    }
    
    model.kPrior       = stackModel.kPrior;       stackModel.kPrior = NULL;
    model.sigmaSqPrior = stackModel.sigmaSqPrior; stackModel.sigmaSqPrior = NULL;
    model.muPrior      = stackModel.muPrior;      stackModel.muPrior = NULL;
    model.treePrior    = stackModel.treePrior;    stackModel.treePrior = NULL;
  }
  
  void invalidateModel(Model& model)
  {
    if (model.kPrior != NULL) {
      delete model.kPrior; model.kPrior = NULL;
    }
      
    delete model.sigmaSqPrior; model.sigmaSqPrior = NULL;
    delete model.muPrior;      model.muPrior = NULL;
    delete model.treePrior;    model.treePrior = NULL;
  }
}

// begging to be an auto_ptr or unique, but I don't want to get into
// C++ versioning
namespace {
  struct DataStackDeconstructor {
    dbarts::VariableType* variableTypes;
    
    DataStackDeconstructor() : variableTypes(NULL) { }
    ~DataStackDeconstructor() {
      delete variableTypes;
    }
  };
}

namespace dbarts {
  void initializeDataFromExpression(Data& data, SEXP dataExpr)
  {
    DataStackDeconstructor stackData;
    
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
    
    data.x = REAL(slotExpr);
    data.numPredictors = static_cast<size_t>(dims[1]);
    
    slotExpr = Rf_getAttrib(dataExpr, Rf_install("varTypes"));
    rc_assertIntConstraints(slotExpr, "variable types", RC_LENGTH | RC_EQ, rc_asRLength(data.numPredictors), RC_END);
    int* i_variableTypes = INTEGER(slotExpr);
    stackData.variableTypes = new VariableType[data.numPredictors];
    for (size_t i = 0; i < data.numPredictors; ++i) stackData.variableTypes[i] = (i_variableTypes[i] == 0 ? ORDINAL : CATEGORICAL);
    data.variableTypes = stackData.variableTypes;
    
    slotExpr = Rf_getAttrib(dataExpr, Rf_install("x.test"));
    if (rc_isS4Null(slotExpr) || Rf_isNull(slotExpr) || rc_getLength(slotExpr) == 0) {
      data.x_test = NULL;
      data.numTestObservations = 0;
    } else {
      if (!Rf_isReal(slotExpr)) Rf_error("x.test must be of type real");
      rc_assertDimConstraints(slotExpr, "dimensions of x.test", RC_LENGTH | RC_EQ, rc_asRLength(2),
                    RC_NA, RC_VALUE | RC_EQ, static_cast<int>(data.numPredictors), RC_END);
      dims = INTEGER(Rf_getAttrib(slotExpr, R_DimSymbol));
      
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
    
    stackData.variableTypes = NULL;
  }
  
  void invalidateData(Data& data)
  {
    delete [] data.variableTypes; data.variableTypes = NULL;
    delete [] data.maxNumCuts;    data.maxNumCuts = NULL;
  }
  
  SEXP createStateExpressionFromFit(const BARTFit& fit)
  {
    const Control& control(fit.control);
    const Data& data(fit.data);
    const Model& model(fit.model);
    const State* state(fit.state);
    
    SEXP treesSym         = Rf_install("trees");
    SEXP treeFitsSym      = Rf_install("treeFits");
    SEXP savedTreesSym    = Rf_install("savedTrees");
    SEXP sigmaSym         = Rf_install("sigma");
    SEXP kSym             = Rf_install("k");
    SEXP rngStateSym      = Rf_install("rng.state");
    
    SEXP result = PROTECT(rc_newList(control.numChains));
    
    SEXP slotExpr;
    
    rc_allocateInSlot2(slotExpr, result, Rf_install("runningTime"), REALSXP, 1);
    REAL(slotExpr)[0] = fit.runningTime;
    
    rc_allocateInSlot2(slotExpr, result, Rf_install("currentNumSamples"), INTSXP, 1);
    INTEGER(slotExpr)[0] = static_cast<int>(fit.currentNumSamples);
    
    rc_allocateInSlot2(slotExpr, result, Rf_install("currentSampleNum"), INTSXP, 1);
    INTEGER(slotExpr)[0] = static_cast<int>(fit.currentSampleNum);
    
    rc_allocateInSlot2(slotExpr, result, Rf_install("numCuts"), INTSXP, data.numPredictors);
    int* numCuts = INTEGER(slotExpr);
    rc_allocateInSlot2(slotExpr, result, Rf_install("cutPoints"), VECSXP, data.numPredictors);
    for (size_t j = 0; j < data.numPredictors; ++j) {
      numCuts[j] = static_cast<int>(fit.numCutsPerVariable[j]);
      SEXP cutPointsExpr = PROTECT(rc_newReal(rc_asRLength(fit.numCutsPerVariable[j])));
      std::memcpy(REAL(cutPointsExpr), fit.cutPoints[j], fit.numCutsPerVariable[j] * sizeof(double));
      SET_VECTOR_ELT(slotExpr, j, cutPointsExpr);
      UNPROTECT(1);
    }
    
    SEXP classDef = PROTECT(R_getClassDef("dbartsState"));
    for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum) {
      SEXP result_i = PROTECT(R_do_new_object(classDef));
      SET_VECTOR_ELT(result, chainNum, result_i);
      UNPROTECT(1);
      
      size_t treeStateLength = state[chainNum].getSerializedTreesLength(fit) / sizeof(int);
      rc_allocateInSlot2(slotExpr, result_i, treesSym, INTSXP, rc_asRLength(treeStateLength));
      state[chainNum].serializeTrees(fit, INTEGER(slotExpr));
      
      rc_allocateInSlot2(slotExpr, result_i, treeFitsSym, REALSXP, rc_asRLength(data.numObservations * control.numTrees));
      rc_setDims(slotExpr, static_cast<int>(data.numObservations), static_cast<int>(control.numTrees), -1);
      std::memcpy(REAL(slotExpr), state[chainNum].treeFits, data.numObservations * control.numTrees * sizeof(double));
      
      if (control.keepTrees) {
        treeStateLength = state[chainNum].getSerializedSavedTreesLength(fit) / sizeof(int);
        rc_allocateInSlot2(slotExpr, result_i, savedTreesSym, INTSXP, rc_asRLength(treeStateLength));
        state[chainNum].serializeSavedTrees(fit, INTEGER(slotExpr));
      } else {
        rc_allocateInSlot(result_i, savedTreesSym, INTSXP, 0);
      }
      
      rc_allocateInSlot2(slotExpr, result_i, sigmaSym, REALSXP, 1);
      REAL(slotExpr)[0] = state[chainNum].sigma;
      
      if (model.kPrior != NULL) {
        rc_allocateInSlot2(slotExpr, result_i, kSym, REALSXP, 1);
        REAL(slotExpr)[0] = state[chainNum].k;
      } else {
        Rf_setAttrib(slotExpr, kSym, R_NilValue);
      }
      
      size_t rngStateLength = ext_rng_getSerializedStateLength(state[chainNum].rng) / sizeof(int);
      rc_allocateInSlot2(slotExpr, result_i, rngStateSym, INTSXP, rc_asRLength(rngStateLength));
      ext_rng_writeSerializedState(state[chainNum].rng, INTEGER(slotExpr));
    }
    
    UNPROTECT(2);
    
    return result;
  }

  void storeStateExpressionFromFit(const BARTFit& fit, SEXP stateExpr)
  {
    const Control& control(fit.control);
    const Data& data(fit.data);
    const Model& model(fit.model);
    const State* state(fit.state);
    
    SEXP treesSym         = Rf_install("trees");
    SEXP treeFitsSym      = Rf_install("treeFits");
    SEXP savedTreesSym    = Rf_install("savedTrees");
    SEXP sigmaSym         = Rf_install("sigma");
    SEXP kSym             = Rf_install("k");
    SEXP rngStateSym      = Rf_install("rng.state");
    
    // check to see if it is an old-style saved object with only a single state
    SEXP classExpr = rc_getClass(stateExpr);
    if (!Rf_isNull(classExpr) && std::strcmp(CHAR(STRING_ELT(classExpr, 0)), "dbartsState") == 0) 
      Rf_error("object from earlier version detected - model must be refit");
    
    if (rc_getLength(stateExpr) != control.numChains)
      Rf_error("length of state list not equal to number of chains");
    
    SEXP slotExpr = Rf_getAttrib(stateExpr, Rf_install("runningTime"));
    REAL(slotExpr)[0] = fit.runningTime;
    
    slotExpr = Rf_getAttrib(stateExpr, Rf_install("currentNumSamples"));
    INTEGER(slotExpr)[0] = static_cast<int>(fit.currentNumSamples);
    
    slotExpr = Rf_getAttrib(stateExpr, Rf_install("currentSampleNum"));
    INTEGER(slotExpr)[0] = static_cast<int>(fit.currentSampleNum);
    
    slotExpr = Rf_getAttrib(stateExpr, Rf_install("numCuts"));
    if (rc_getLength(slotExpr) != data.numPredictors) {
      rc_allocateInSlot2(slotExpr, stateExpr, Rf_install("numCuts"), INTSXP, data.numPredictors);
      int* numCuts = INTEGER(slotExpr);
      for (size_t j = 0; j < data.numPredictors; ++j) numCuts[j] = static_cast<int>(fit.numCutsPerVariable[j]);
    }
    slotExpr = Rf_getAttrib(stateExpr, Rf_install("cutPoints"));
    if (rc_getLength(slotExpr) != data.numPredictors) {
      rc_allocateInSlot2(slotExpr, stateExpr, Rf_install("cutPoints"), VECSXP, data.numPredictors);
      for (size_t j = 0; j < data.numPredictors; ++j) {
        SEXP cutPointsExpr = PROTECT(rc_newReal(rc_asRLength(fit.numCutsPerVariable[j])));
        std::memcpy(REAL(cutPointsExpr), fit.cutPoints[j], fit.numCutsPerVariable[j] * sizeof(double));
        SET_VECTOR_ELT(slotExpr, j, cutPointsExpr);
        UNPROTECT(1);
      }
    } else {
      for (size_t j = 0; j < data.numPredictors; ++j) {
        SEXP cutPointsExpr = VECTOR_ELT(slotExpr, j);
        if (rc_getLength(cutPointsExpr) != static_cast<size_t>(fit.numCutsPerVariable[j])) {
          cutPointsExpr = PROTECT(rc_newReal(rc_asRLength(fit.numCutsPerVariable[j])));
          std::memcpy(REAL(cutPointsExpr), fit.cutPoints[j], fit.numCutsPerVariable[j] * sizeof(double));
          SET_VECTOR_ELT(slotExpr, j, cutPointsExpr);
          UNPROTECT(1);
        } else {
          std::memcpy(REAL(cutPointsExpr), fit.cutPoints[j], fit.numCutsPerVariable[j] * sizeof(double));
        }
      }
    }
    
    for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum)
    {
      SEXP stateExpr_i = VECTOR_ELT(stateExpr, rc_asRLength(chainNum));
      classExpr = rc_getClass(stateExpr_i);
      if (std::strcmp(CHAR(STRING_ELT(classExpr, 0)), "dbartsState") != 0)
        Rf_error("'state' not of class 'dbartsState'");
      
      
      slotExpr = Rf_getAttrib(stateExpr_i, treeFitsSym);
      SEXP dimsExpr = rc_getDims(slotExpr);
      if (rc_getLength(dimsExpr) != 2) Rf_error("dimensions of state@treeFits indicate that it is not a matrix");
      int* dims = INTEGER(dimsExpr);
       
      if (static_cast<size_t>(dims[0]) != data.numObservations || static_cast<size_t>(dims[1]) != control.numTrees) {
        rc_allocateInSlot2(slotExpr, stateExpr_i, treeFitsSym, REALSXP, rc_asRLength(data.numObservations * control.numTrees));
        rc_setDims(slotExpr, static_cast<int>(data.numObservations), static_cast<int>(control.numTrees), -1);
      }
      
      size_t treeStateLength = state[chainNum].getSerializedTreesLength(fit) / sizeof(int);
      rc_allocateInSlot2(slotExpr, stateExpr_i, treesSym, INTSXP, rc_asRLength(treeStateLength));
      state[chainNum].serializeTrees(fit, INTEGER(slotExpr));
           
      slotExpr = Rf_getAttrib(stateExpr_i, treeFitsSym);
      std::memcpy(REAL(slotExpr), state[chainNum].treeFits, data.numObservations * control.numTrees * sizeof(double));
      
      if (control.keepTrees) {
        treeStateLength = state[chainNum].getSerializedSavedTreesLength(fit) / sizeof(int);
        rc_allocateInSlot2(slotExpr, stateExpr_i, savedTreesSym, INTSXP, rc_asRLength(treeStateLength));
        state[chainNum].serializeSavedTrees(fit, INTEGER(slotExpr));
      } else {
        rc_allocateInSlot(stateExpr_i, savedTreesSym, INTSXP, 0);
      }
            
      slotExpr = Rf_getAttrib(stateExpr_i, sigmaSym);
      REAL(slotExpr)[0] = state[chainNum].sigma;
      
      if (model.kPrior != NULL) {
        slotExpr = Rf_getAttrib(stateExpr_i, kSym);
        REAL(slotExpr)[0] = state[chainNum].k;
      }
      
      size_t rngStateLength = ext_rng_getSerializedStateLength(state[chainNum].rng) / sizeof(int);
      slotExpr = Rf_getAttrib(stateExpr_i, rngStateSym);
      if (rc_getLength(slotExpr) != rngStateLength)
        rc_allocateInSlot2(slotExpr, stateExpr_i, rngStateSym, INTSXP, rc_asRLength(rngStateLength));
      ext_rng_writeSerializedState(state[chainNum].rng, INTEGER(slotExpr));
    }
  }
  
  void initializeStateFromExpression(BARTFit& fit, State* state, SEXP stateExpr)
  {
    const Control& control(fit.control);
    const Data& data(fit.data);
    const Model& model(fit.model);
    
    // check to see if it is an old-style saved object with only a single state
    SEXP classExpr = rc_getClass(stateExpr);
    if (!Rf_isNull(classExpr) && std::strcmp(CHAR(STRING_ELT(classExpr, 0)), "dbartsState") == 0) 
      Rf_error("object from earlier version detected - model must be refit");
    
    SEXP slotExpr = Rf_getAttrib(stateExpr, Rf_install("runningTime"));
    fit.runningTime = REAL(slotExpr)[0];
    
    fit.currentSampleNum = static_cast<size_t>(INTEGER(Rf_getAttrib(stateExpr, Rf_install("currentSampleNum")))[0]);
    
    size_t currentNumSamples = static_cast<size_t>(INTEGER(Rf_getAttrib(stateExpr, Rf_install("currentNumSamples")))[0]);
    if (fit.currentNumSamples != currentNumSamples && control.keepTrees) {
      for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum)
        state[chainNum].resize(fit, currentNumSamples);
      fit.currentSampleNum = 0;
    }
    fit.currentNumSamples = currentNumSamples;
    
    for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum)
    {
      SEXP stateExpr_i = VECTOR_ELT(stateExpr, rc_asRLength(chainNum));
      classExpr = rc_getClass(stateExpr_i);
      if (std::strcmp(CHAR(STRING_ELT(classExpr, 0)), "dbartsState") != 0)
        Rf_error("'state' not of class 'dbartsState'");
      
      SEXP slotExpr = Rf_getAttrib(stateExpr_i, Rf_install("trees"));
      state[chainNum].deserializeTrees(fit, INTEGER(slotExpr));
      
      slotExpr = Rf_getAttrib(stateExpr_i, Rf_install("treeFits"));
      std::memcpy(state[chainNum].treeFits, const_cast<const double*>(REAL(slotExpr)), data.numObservations * control.numTrees * sizeof(double));
      
      if (control.keepTrees) {
        slotExpr = Rf_getAttrib(stateExpr_i, Rf_install("savedTrees"));
        state[chainNum].deserializeSavedTrees(fit, INTEGER(slotExpr));
      }
            
      slotExpr = Rf_getAttrib(stateExpr_i, Rf_install("sigma"));
      state[chainNum].sigma = REAL(slotExpr)[0];
      
      if (model.kPrior != NULL) {
        slotExpr = Rf_getAttrib(stateExpr_i, Rf_install("k"));
        state[chainNum].k = REAL(slotExpr)[0];
      }
      
      ext_rng_readSerializedState(state[chainNum].rng, INTEGER(Rf_getAttrib(stateExpr_i, Rf_install("rng.state"))));
    }
    
    uint32_t* numCuts        = new uint32_t[data.numPredictors];
    const double** cutPoints = new const double*[data.numPredictors];
    size_t* columns          = new size_t[data.numPredictors];
    
    const int* numCuts_i = INTEGER(Rf_getAttrib(stateExpr, Rf_install("numCuts")));
    
    slotExpr = Rf_getAttrib(stateExpr, Rf_install("cutPoints"));
    for (size_t j = 0; j < data.numPredictors; ++j) {
      numCuts[j] = static_cast<uint32_t>(numCuts_i[j]);
      cutPoints[j] = REAL(VECTOR_ELT(slotExpr, j));
      columns[j] = j;
    }
    fit.setCutPoints(cutPoints, numCuts, columns, data.numPredictors);
    
    delete [] columns;
    delete [] cutPoints;
    delete [] numCuts;
    
    fit.rebuildScratchFromState();
  }
}

