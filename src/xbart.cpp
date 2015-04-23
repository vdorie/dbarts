#include "xbart.hpp"

#include <cstddef> // size_t
#include <cstring> // memcpy
#include <cmath>   // sqrt
#include <algorithm> // sort

#include <external/random.h>
#include <external/stats.h>
#include <rc/bounds.h>
#include <rc/util.h>

#include <dbarts/bartFit.hpp>
#include <dbarts/control.hpp>
#include <dbarts/data.hpp>
#include <dbarts/model.hpp>

#include <Rinternals.h> // R_xlen_t

using std::size_t;
using namespace dbarts;

#define Z_(_X_) static_cast<R_xlen_t>(_X_)

namespace {

  void permuteIndexArray(ext_rng* generator, size_t* indices, size_t length);
}

// #include <dbarts/cstdint.hpp>
#include <dbarts/scratch.hpp>
#include <dbarts/state.hpp>
#include <dbarts/results.hpp>
#include "dbarts/node.hpp"
#include "dbarts/tree.hpp"

namespace {
  SEXP allocateResult(bool dropUnusedDims, size_t numNTrees, size_t numKs, size_t numPowers, size_t numBases, size_t numReps);
  void allocateDataStorage(const Data& origData, Data& repData, size_t inSampleSize, size_t outSampleSize);
  void allocateModelStorage(const Model& origModel, Model& repModel);
  void divideData(const Data& origData, Data& repData, double* targetY, size_t inSampleSize, size_t outSampleSize, const size_t* permutation);
  void freeDataStorage(Data& repData);
  void freeModelStorage(Model& repModel);
  
  /* struct ResultsCalculator {
    size_t (*getNumResults)(void);
    void (*computeResults)(const double* response, size_t numResponse, const double* predictions, size_t numPredictions, double* results);
    void* storage;
  }; */
}

namespace dbarts {
  
  SEXP xbart(SEXP fitExpr, SEXP ntreeExpr, SEXP kExpr, SEXP powerExpr, SEXP baseExpr, SEXP nskipExpr,
             SEXP KExpr, SEXP nrepsExpr, SEXP resultTypeExpr, SEXP dropExpr)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) error("xbart called on NULL external pointer");
    
    Control origControl = fit->control;
    Data    origData    = fit->data;
    Model   origModel   = fit->model;
    
    if (origData.numObservations == 0) error("xbart called on empty data set");
    
    rc_checkInts(ntreeExpr, "num trees", RC_LENGTH | RC_GEQ, Z_(1), RC_VALUE | RC_GT, 0, RC_NA | RC_NO, RC_END);
    rc_checkDoubles(kExpr, "k", RC_LENGTH | RC_GEQ, Z_(1), RC_VALUE | RC_GT, 0.0, RC_NA | RC_NO, RC_END);
    rc_checkDoubles(powerExpr, "power", RC_LENGTH | RC_GEQ, Z_(1), RC_VALUE | RC_GT, 0.0, RC_NA | RC_NO, RC_END);
    rc_checkDoubles(baseExpr, "power", RC_LENGTH | RC_GEQ, Z_(1), RC_VALUE | RC_GT, 0.0, RC_VALUE | RC_LT, 1.0, RC_NA | RC_NO, RC_END);
    rc_checkInts(nskipExpr, "num skip",  RC_LENGTH | RC_GEQ, Z_(1), RC_LENGTH | RC_LEQ, Z_(2), RC_VALUE | RC_GEQ, 0, RC_NA | RC_NO, RC_END);
    
    size_t numFolds = static_cast<size_t>(
      rc_getInt(KExpr, "num folds", RC_LENGTH | RC_EQ, Z_(1), RC_VALUE | RC_GT, 0, RC_VALUE | RC_LEQ, static_cast<int>(origData.numObservations) - 1, RC_NA | RC_NO, RC_END));
    
    size_t numReps = static_cast<size_t>(
      rc_getInt(nrepsExpr, "num reps", RC_LENGTH | RC_GEQ, Z_(1), RC_VALUE | RC_GT, 0, RC_NA | RC_NO, RC_END));
    
    bool dropUnusedDims = rc_getBool(dropExpr, "drop", RC_LENGTH | RC_EQ, Z_(1), RC_NA | RC_NO, RC_END);
    
    size_t outSampleSize = origData.numObservations / numFolds;
    size_t inSampleSize  = origData.numObservations - outSampleSize;
    
    size_t numNTrees = rc_getLength(ntreeExpr);
    size_t numKs     = rc_getLength(kExpr);
    size_t numPowers = rc_getLength(powerExpr);
    size_t numBases  = rc_getLength(baseExpr);
    
    int* ntree = INTEGER(ntreeExpr);    
    double* k     = REAL(kExpr);
    double* power = REAL(powerExpr);
    double* base  = REAL(baseExpr);
    
    SEXP resultExpr = allocateResult(dropUnusedDims, numNTrees, numKs, numPowers, numBases, numReps); // increases protect count by 1
    double* result = REAL(resultExpr);
    
    Data repData;
    Control repControl = fit->control;
    Model repModel;
    
    allocateDataStorage(origData, repData, inSampleSize, outSampleSize);
    allocateModelStorage(origModel, repModel);
    double* targetY = new double[outSampleSize];
    double* probabilities = new double[origControl.numSamples];
    
    size_t* permutation = new size_t[origData.numObservations];
    for (size_t i = 0; i < origData.numObservations; ++i) permutation[i] = i;
    
    size_t numInitialBurnIn    = static_cast<size_t>(INTEGER(nskipExpr)[0]);
    size_t numSubsequentBurnIn = XLENGTH(nskipExpr) == 2 ? static_cast<size_t>(INTEGER(nskipExpr)[1]) : numInitialBurnIn / 5;
    
    size_t resultIndex = 0;
    
    for (size_t nIndex = 0; nIndex < numNTrees; ++nIndex) {
      repControl.numTrees = static_cast<size_t>(ntree[nIndex]);
      fit->setControl(repControl);
      
      for (size_t kIndex = 0; kIndex < numKs; ++kIndex) {
        double endNodeSd = (origControl.responseIsBinary ? 3.0 : 0.5) /  (k[kIndex] * std::sqrt(static_cast<double>(repControl.numTrees)));
        static_cast<NormalPrior*>(repModel.muPrior)->precision = 1.0 / (endNodeSd * endNodeSd);
        
        for (size_t pIndex = 0; pIndex < numPowers; ++pIndex) {
          static_cast<CGMPrior*>(repModel.treePrior)->power = power[pIndex];
          
          for (size_t bIndex = 0; bIndex < numBases; ++bIndex) {
            static_cast<CGMPrior*>(repModel.treePrior)->base = base[bIndex];
            fit->setModel(repModel);
            
            for (size_t repIndex = 0; repIndex < numReps; ++repIndex) {
              permuteIndexArray(origControl.rng, permutation, origData.numObservations);
              std::sort(permutation, permutation + outSampleSize);
              std::sort(permutation + outSampleSize, permutation + origData.numObservations);
              
              divideData(origData, repData, targetY, inSampleSize, outSampleSize, permutation);
              
              fit->setData(repData);
              
              Results* samples = fit->runSampler(resultIndex == 0 ? numInitialBurnIn : numSubsequentBurnIn, origControl.numSamples);
              size_t tp = 0, fp = 0, tn = 0, fn = 0;
              for (size_t i = 0; i < outSampleSize; ++i) {
                for (size_t j = 0; j < origControl.numSamples; ++j) {
                  probabilities[j] = ext_cumulativeProbabilityOfNormal(samples->trainingSamples[j + i * origControl.numSamples], 0.0, 1.0);
                }
                double prediction = ext_computeMean(probabilities, origControl.numSamples) > 0.5 ? 1.0 : 0.0;
                
                if (prediction == 1.0) {
                  if (targetY[i] == 1.0) ++tp; else ++fp;
                } else {
                  if (targetY[i] == 0.0) ++tn; else ++fn;
                }
              }
              result[resultIndex++] = static_cast<double>(fp + fn) / static_cast<double>(outSampleSize);
              
              
              delete samples;
            }
          }
        }
      }
    }
    
    fit->setControl(origControl);
    fit->setModel(origModel);
    fit->setData(origData);
    
    delete [] permutation;
    
    delete [] probabilities;
    delete [] targetY;
    
    freeModelStorage(repModel);
    freeDataStorage(repData);
    
    UNPROTECT(1);
       
    return resultExpr;
  }
}

namespace {
  void permuteIndexArray(ext_rng* generator, size_t* indices, size_t length)
  {
    size_t temp, swapPos;
    for (size_t i = 0; i < length - 1; ++i) {
      swapPos = static_cast<size_t>(ext_rng_simulateUnsignedIntegerUniformInRange(generator, i, length));
      
      temp = indices[i];
      indices[i] = indices[swapPos];
      indices[swapPos] = temp;
    }
  }
  
  void divideData(const Data& origData, Data& repData, double* targetY, size_t inSampleSize, size_t outSampleSize, const size_t* permutation)
  {
    size_t i, j, obsIndex;
    double* y = const_cast<double*>(repData.y);
    double* X = const_cast<double*>(repData.X);
    double* X_test = const_cast<double*>(repData.X_test);
    
    for (i = 0; i < outSampleSize; ++i) {
      obsIndex = permutation[i];
      targetY[i] = origData.y[obsIndex];
      for (j = 0; j < origData.numPredictors; ++j) {
        X_test[j + i * outSampleSize] = origData.X[j + obsIndex * origData.numObservations];
      }
    }
    for (i = 0; i < inSampleSize; ++i) {
      obsIndex = permutation[i + outSampleSize];
      y[i] = origData.y[obsIndex];
      for (j = 0; j < origData.numPredictors; ++j) {
        X[j + i * inSampleSize] = origData.X[j + obsIndex * origData.numObservations];
      }
    }
  }
  
  SEXP allocateResult(bool dropUnusedDims, size_t numNTrees, size_t numKs, size_t numPowers, size_t numBases, size_t numReps)
  {
    SEXP resultExpr = PROTECT(allocVector(REALSXP, numNTrees * numKs * numPowers * numBases * numReps));
    if (dropUnusedDims) {
      size_t numDims = 1 + (numNTrees > 1 ? 1 : 0) + (numKs > 1 ? 1 : 0) + (numPowers > 1 ? 1 : 0) + (numBases > 1 ? 1 : 0);
      if (numDims > 1) {
        SEXP dimsExpr = allocVector(INTSXP, numDims);
        int* dims = INTEGER(dimsExpr);
        numDims = 0;
        if (numNTrees > 0) dims[numDims++] = static_cast<int>(numNTrees);
        if (numKs > 0)     dims[numDims++] = static_cast<int>(numKs);
        if (numPowers > 0) dims[numDims++] = static_cast<int>(numPowers);
        if (numBases > 0)  dims[numDims++] = static_cast<int>(numBases);
        dims[numDims] = static_cast<int>(numReps);
        
        R_do_slot_assign(resultExpr, R_DimSymbol, dimsExpr);
      }
    } else {
      rc_setDims(resultExpr, static_cast<int>(numNTrees), static_cast<int>(numKs), static_cast<int>(numPowers),
                 static_cast<int>(numBases), static_cast<int>(numReps), -1);
    }
    return resultExpr;
  }

  
  void allocateDataStorage(const Data& origData, Data& repData, size_t inSampleSize, size_t outSampleSize)
  {
    repData.y = new double[inSampleSize];
    repData.X = new double[inSampleSize * origData.numPredictors];
    repData.X_test = new double[outSampleSize * origData.numPredictors];
     
    repData.weights    = origData.weights != NULL ? new double[inSampleSize]  : NULL;
    repData.offset     = origData.offset  != NULL ? new double[inSampleSize]  : NULL;
    repData.testOffset = origData.offset  != NULL ? new double[outSampleSize] : NULL;
    
    repData.numObservations     = inSampleSize;
    repData.numPredictors       = origData.numPredictors;
    repData.numTestObservations = outSampleSize;
    repData.sigmaEstimate       = origData.sigmaEstimate;
    
    repData.variableTypes = origData.variableTypes;
    repData.maxNumCuts    = origData.maxNumCuts;
  }
  
  void allocateModelStorage(const Model& origModel, Model& repModel)
  {
    repModel.birthOrDeathProbability = origModel.birthOrDeathProbability;
    repModel.swapProbability = origModel.swapProbability;
    repModel.changeProbability = origModel.changeProbability;
    repModel.birthProbability = origModel.birthProbability;
    
    CGMPrior* repTreePrior = new CGMPrior();
    const CGMPrior* oldTreePrior = static_cast<CGMPrior*>(origModel.treePrior);
    repTreePrior->base = oldTreePrior->base;
    repTreePrior->power = oldTreePrior->power;
    
    repModel.treePrior = repTreePrior;
    
    
    NormalPrior* repNodePrior = new NormalPrior();
    const NormalPrior* oldNodePrior = static_cast<NormalPrior*>(origModel.muPrior);
    repNodePrior->precision = oldNodePrior->precision;
    
    repModel.muPrior = repNodePrior;
    
    
    ChiSquaredPrior* repResidPrior = new ChiSquaredPrior();
    const ChiSquaredPrior* oldResidPrior = static_cast<ChiSquaredPrior*>(origModel.sigmaSqPrior);
    repResidPrior->degreesOfFreedom = oldResidPrior->degreesOfFreedom;
    repResidPrior->scale = oldResidPrior->scale;
    
    repModel.sigmaSqPrior = repResidPrior;
  }
  
  void freeDataStorage(Data& repData)
  {
    delete [] repData.testOffset; repData.testOffset = NULL;
    delete [] repData.offset;     repData.offset     = NULL;
    delete [] repData.weights;    repData.weights    = NULL;
    
    delete [] repData.X_test; repData.X_test = NULL;
    delete [] repData.X;      repData.X      = NULL;
    delete [] repData.y;      repData.y      = NULL;
  }
  
  void freeModelStorage(Model& repModel)
  {
    delete repModel.sigmaSqPrior;
    delete repModel.muPrior;
    delete repModel.treePrior;
  }
}

