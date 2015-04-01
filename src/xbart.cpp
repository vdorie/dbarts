#include "xbart.hpp"

#include <cstddef> // size_t
#include <cstring> // memcpy

#include <external/random.h>
#include <rc/bounds.h>

#include <dbarts/bartFit.hpp>
#include <dbarts/control.hpp>
#include <dbarts/data.hpp>

#include <Rinternals.h> // R_xlen_t

using std::size_t;
using namespace dbarts;

#define Z_(_X_) static_cast<R_xlen_t>(_X_)

namespace {

  void permuteIndexArray(ext_rng* generator, size_t* indices, size_t length);
}

#include <dbarts/state.hpp>
#include <dbarts/scratch.hpp>
#include "dbarts/node.hpp"
#include "dbarts/tree.hpp"
#include <dbarts/cstdint.hpp>

namespace {
  
  using namespace dbarts;
  bool nodeSplitsOnVariableAndCut(Node& node, int32_t varIndex, int32_t splitIndex);
}

namespace dbarts {
  
  SEXP xbart(SEXP fitExpr, SEXP kExpr, SEXP powerExpr, SEXP baseExpr, SEXP ntreeExpr, SEXP nskipExpr,
             SEXP KExpr, SEXP resultTypeExpr, SEXP dropExpr)
  {
    BARTFit* fit = static_cast<BARTFit*>(R_ExternalPtrAddr(fitExpr));
    if (fit == NULL) error("xbart called on NULL external pointer");
    
    rc_checkInts(ntreeExpr, "num trees", RC_LENGTH | RC_GEQ, Z_(1), RC_VALUE | RC_GT, 0, RC_NA | RC_NO, RC_END);
    rc_checkInts(nskipExpr, "num skip", RC_LENGTH | RC_GEQ, Z_(1), RC_LENGTH | RC_LEQ, Z_(2), RC_VALUE | RC_GEQ, 0, RC_NA | RC_NO, RC_END);
    
    size_t numFolds = static_cast<size_t>(
      rc_getInt(KExpr, "num folds", RC_LENGTH | RC_EQ, Z_(1), RC_VALUE | RC_GT, 0, RC_NA | RC_NO, RC_END));
    
    size_t numKs     = (size_t) LENGTH(kExpr);
    size_t numPowers = (size_t) LENGTH(powerExpr);
    size_t numBases  = (size_t) LENGTH(baseExpr);
    size_t numNTrees = (size_t) LENGTH(ntreeExpr);
    
    double* k     = REAL(kExpr);
    double* power = REAL(powerExpr);
    double* base  = REAL(baseExpr);
    int* ntree = INTEGER(ntreeExpr);
    
    Control& control(fit->control);
    Data& data(fit->data);
    
    const double* x_orig = data.X;
    const double* y_orig = data.y;
  
    Data repData;
    
    size_t* permutation = new size_t[data.numObservations];
    for (size_t i = 0; i < data.numObservations; ++i) permutation[i] = i;
    
    permuteIndexArray(control.rng, permutation, data.numObservations);
    
    delete [] permutation;
        
    State& state(fit->state);
    Scratch& scratch(fit->scratch);
    
    control.verbose = false;
    size_t iter = 0;
    while (iter < 1000) {
      fit->runSampler(1, 0);
      ++iter;
      //bool splitsHappy = false;
      //for (size_t j = 0; j < control.numTrees; ++j) {
      //  if (nodeSplitsOnVariableAndCut(state.trees[j].top, 1, 4)) { splitsHappy = true; break; }
      //}
      // if (splitsHappy) break;
      
    }
    control.verbose = true;
    
    double* x_rep = new double[90 * 2];
    double* y_rep = new double[90];
    
    for (size_t i = 0; i < 10; ++i) {
      x_rep[i] = x_orig[i];
      x_rep[i + 90] = x_orig[i + 100];
      y_rep[i] = y_orig[i];
    }
    for (size_t i = 10; i < 90; ++i) {
      x_rep[i] = x_orig[i + 10];
      x_rep[i + 90] = x_orig[i + 110];
      y_rep[i] = y_orig[i + 10];
    }
    
    repData.X = x_rep;
    repData.y = y_rep;
    repData.numObservations = 90;
    repData.numPredictors = 2;
    repData.sigmaEstimate = data.sigmaEstimate;
    repData.variableTypes = data.variableTypes;
    repData.maxNumCuts = data.maxNumCuts;
    
    fit->setData(repData);
    
    fit->setData(data);
    
    delete [] y_rep;
    delete [] x_rep;
    
    return R_NilValue;
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
}

namespace {
  using namespace dbarts;
  bool nodeSplitsOnVariableAndCut(Node& node, int32_t varIndex, int32_t splitIndex) {
    if (node.isBottom()) return false;
    
    if (node.p.rule.variableIndex == varIndex && node.p.rule.splitIndex == splitIndex) return true;

    if (nodeSplitsOnVariableAndCut(*node.getLeftChild(), varIndex, splitIndex)) return true;
    return nodeSplitsOnVariableAndCut(*node.getRightChild(), varIndex, splitIndex);
  }
}

