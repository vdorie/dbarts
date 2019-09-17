#include "config.hpp"
#include <dbarts/bartFit.hpp>

#include <cmath>     // sqrt
#include <cstring>   // memcpy
#include <cstddef>   // size_t
#include <limits>    // nan

#ifdef __INTEL_COMPILER
#  define __need_timespec 1
#endif
#include <time.h>

#if !defined(HAVE_SYS_TIME_H) && defined(HAVE_GETTIMEOFDAY)
#  undef HAVE_GETTIMEOFDAY
#endif
#ifdef HAVE_SYS_TIME_H
#  include <sys/time.h> // gettimeofday
#endif


#include <set>       // used to sort and find 
#include <vector>    //   split points
#include <algorithm> // integer min

#include <misc/alloca.h>
#include <misc/linearAlgebra.h>
#include <misc/stats.h>

#include <external/io.h>
#include <external/linearAlgebra.h>
#include <external/random.h>
#include <external/stats.h>

#include <dbarts/results.hpp>
#include "functions.hpp"
#include "tree.hpp"

using std::size_t;
using std::uint32_t;

namespace {
  using namespace dbarts;

  void allocateMemory(BARTFit& fit);
  void createRNG(BARTFit& fit);
  void destroyRNG(BARTFit& fit);
  void setInitialCutPoints(BARTFit& fit);
  void setXIntegerCutMap(BARTFit& fit);
  void setXIntegerCutMap(BARTFit& fit, const size_t* columns, size_t numColumns);
  void setXTestIntegerCutMap(const BARTFit& fit, const double* x_test, size_t numTestObservations, xint_t* xt_test);
  void setXTestIntegerCutMap(const BARTFit& fit, const double* x_test, size_t numTestObservations,
                             xint_t* xt_test, const size_t* columns, size_t numColumns);
  void setInitialFit(BARTFit& fit);
  
  void setPrior(BARTFit& fit);
  
  void setCutPoints(BARTFit& fit, const size_t* columns, size_t numColumns);
  void setCutPointsFromQuantiles(BARTFit& fit, const double* x, uint32_t maxNumCuts,
                                 uint32_t& numCutsPerVariable, double*& cutPoints,
                                 std::set<double>& uniqueElements, std::vector<double>& sortedElements);
  void setCutPointsUniformly(BARTFit& fit, const double* x, uint32_t maxNumCuts,
                             uint32_t& numCutsPerVariable, double*& cutPoints);
  
  void printInitialSummary(const BARTFit& fit);
  void printTerminalSummary(const BARTFit& fit);
  
  void initializeLatents(BARTFit& fit);
  void initializeLatents(BARTFit& fit, size_t chainNum);
  void rescaleResponse(BARTFit& fit);
  
  // void resampleTreeFits(BARTFit& fit);
  
  void sampleProbitLatentVariables(const BARTFit& fit, State& state, const double* fits, double* yRescaled);
  void storeSamples(const BARTFit& fit, size_t chainNum, Results& results,
                    const double* trainingSample, const double* testSample,
                    double sigma, double k, const uint32_t* variableCounts, size_t simNum);
  void countVariableUses(const BARTFit& fit, const State& state, uint32_t* variableCounts);
  
#ifdef HAVE_SYS_TIME_H
  double subtractTimes(struct timeval end, struct timeval start);
#else
  double subtractTimes(time_t end, time_t start);
#endif
}

namespace dbarts {
  
  void BARTFit::rebuildScratchFromState()
  {
    for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum) {
      
      misc_setVectorToConstant(chainScratch[chainNum].totalFits, data.numObservations, 0.0);
      
      for (size_t treeNum = 0; treeNum < control.numTrees; ++treeNum)
        misc_addVectorsInPlace(const_cast<const double*>(state[chainNum].treeFits + treeNum * data.numObservations),
                               data.numObservations, 1.0,
                               chainScratch[chainNum].totalFits);
      
      
      if (data.numTestObservations > 0) {
        double* testFits = new double[data.numTestObservations];
        
        misc_setVectorToConstant(chainScratch[chainNum].totalTestFits, data.numTestObservations, 0.0);
        
        for (size_t treeNum = 0; treeNum < control.numTrees; ++treeNum) {
          double* treeFits = state[chainNum].treeFits + treeNum * data.numObservations;
        
          // next allocates memory
          double* nodeParams = state[chainNum].trees[treeNum].recoverParametersFromFits(*this, treeFits);
          
          state[chainNum].trees[treeNum].setCurrentFitsFromParameters(*this, nodeParams, treeFits, testFits);
          
          misc_addVectorsInPlace(const_cast<const double*>(testFits), data.numTestObservations, 1.0, chainScratch[chainNum].totalTestFits);
          
          delete [] nodeParams;
        }
        
        delete [] testFits;
      }
    }
  }
  
  void BARTFit::setResponse(const double* newY) {
    
    if (!control.responseIsBinary) {
      double* sigmaUnscaled = misc_stackAllocate(control.numChains, double);
      for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum)
        sigmaUnscaled[chainNum] = state[chainNum].sigma * sharedScratch.dataScale.range;
      
      double priorUnscaled = model.sigmaSqPrior->getScale() * sharedScratch.dataScale.range * sharedScratch.dataScale.range;
      
      data.y = newY;
      
      rescaleResponse(*this);
      
      model.sigmaSqPrior->setScale(priorUnscaled / (sharedScratch.dataScale.range * sharedScratch.dataScale.range));
      
      for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum)
        state[chainNum].sigma = sigmaUnscaled[chainNum] / sharedScratch.dataScale.range;
      
      misc_stackFree(sigmaUnscaled);
     
    } else {
      data.y = newY;
      
      for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum)
        sampleProbitLatentVariables(*this, state[chainNum], const_cast<const double*>(chainScratch[chainNum].totalFits), chainScratch[chainNum].probitLatents);
    }
  }
  
  void BARTFit::setOffset(const double* newOffset) {
    if (!control.responseIsBinary) {
      double* sigmaUnscaled = misc_stackAllocate(control.numChains, double);
      for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum)
        sigmaUnscaled[chainNum] = state[chainNum].sigma * sharedScratch.dataScale.range;
      
      double priorUnscaled = model.sigmaSqPrior->getScale() * sharedScratch.dataScale.range * sharedScratch.dataScale.range;
      
      data.offset = newOffset;
      
      rescaleResponse(*this);
      
      model.sigmaSqPrior->setScale(priorUnscaled / (sharedScratch.dataScale.range * sharedScratch.dataScale.range));
      
      for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum)
        state[chainNum].sigma = sigmaUnscaled[chainNum] / sharedScratch.dataScale.range;
      
      misc_stackFree(sigmaUnscaled);
      
    } else {
      data.offset = newOffset;
      
      // turns out, this isn't necessary for a valid statistical model
      // for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum)
      //   sampleProbitLatentVariables(*this, state[chainNum], const_cast<const double*>(chainScratch[chainNum].totalFits), chainScratch[chainNum].probitLatents);
    }
  }
  
  void BARTFit::setWeights(const double* newWeights) {
    data.weights = newWeights;
  }
  
  void BARTFit::predict(const double* x_test, size_t numTestObservations, const double* testOffset, double* result) const
  {
    double* currTestFits = new double[numTestObservations];
    double* totalTestFits = new double[numTestObservations];
    
    if (control.keepTrees) {
      double* xt_test = new double[numTestObservations * data.numPredictors];
      misc_transposeMatrix(x_test, numTestObservations, data.numPredictors, xt_test);
      
      for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum) {
        for (size_t sampleNum = 0; sampleNum < currentNumSamples; ++sampleNum) {
          
          misc_setVectorToConstant(totalTestFits, numTestObservations, 0.0);
          
          for (size_t treeNum = 0; treeNum < control.numTrees; ++treeNum) {
            size_t treeOffset = treeNum + sampleNum * control.numTrees;
            
            state[chainNum].savedTrees[treeOffset].getPredictions(*this, xt_test, numTestObservations, currTestFits);
            
            misc_addVectorsInPlace(const_cast<const double*>(currTestFits), numTestObservations, 1.0, totalTestFits);
          }
          
          double* result_i = result + (sampleNum + chainNum * currentNumSamples) * numTestObservations;
          if (control.responseIsBinary) {
            std::memcpy(result_i, const_cast<const double*>(totalTestFits), numTestObservations * sizeof(double));
          } else {
            misc_setVectorToConstant(result_i, numTestObservations, sharedScratch.dataScale.range * 0.5 + sharedScratch.dataScale.min);
            misc_addVectorsInPlace(const_cast<const double*>(totalTestFits), numTestObservations, sharedScratch.dataScale.range, result_i);
          }
          
          if (testOffset != NULL) misc_addVectorsInPlace(testOffset, numTestObservations, 1.0, result_i);
        }
      }
      delete [] xt_test;
    } else {
      xint_t* xt_test = new xint_t[numTestObservations * data.numPredictors];
      setXTestIntegerCutMap(*this, x_test, numTestObservations, xt_test);

      for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum) {
      
        misc_setVectorToConstant(totalTestFits, numTestObservations, 0.0);
        
        for (size_t treeNum = 0; treeNum < control.numTrees; ++treeNum) {
          const double* treeFits = state[chainNum].treeFits + treeNum * data.numObservations;
          const double* nodeParams = state[chainNum].trees[treeNum].recoverParametersFromFits(*this, treeFits);
          
          state[chainNum].trees[treeNum].setCurrentFitsFromParameters(*this, nodeParams, xt_test, numTestObservations, currTestFits);
          
          misc_addVectorsInPlace(const_cast<const double*>(currTestFits), numTestObservations, 1.0, totalTestFits);
        }
        
        double* result_i = result + chainNum * numTestObservations;
        if (control.responseIsBinary) {
          std::memcpy(result_i, const_cast<const double*>(totalTestFits), numTestObservations * sizeof(double));
        } else {
          misc_setVectorToConstant(result_i, numTestObservations, sharedScratch.dataScale.range * 0.5 + sharedScratch.dataScale.min);
          misc_addVectorsInPlace(const_cast<const double*>(totalTestFits), numTestObservations, sharedScratch.dataScale.range, result_i);
        }
        
        if (testOffset != NULL) misc_addVectorsInPlace(testOffset, numTestObservations, 1.0, result_i);
      }
      delete [] xt_test;
    }
    
    delete [] totalTestFits;
    delete [] currTestFits;
  }
}

namespace {
  void forceUpdateTrees(BARTFit& fit) {
    const Control& control(fit.control);
    const Data& data(fit.data);
    State* state(fit.state);
    ChainScratch* chainScratch(fit.chainScratch);
    
    for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum) {
      misc_setVectorToConstant(chainScratch[chainNum].totalFits, data.numObservations, 0.0);
      
      for (size_t treeNum = 0; treeNum < control.numTrees; ++treeNum) {
        double* treeFits = const_cast<double*>(state[chainNum].treeFits + treeNum * data.numObservations);
        
        // duplicate old node parameters
        double* nodeParams = state[chainNum].trees[treeNum].recoverParametersFromFits(fit, treeFits);
        
        state[chainNum].trees[treeNum].top.addObservationsToChildren(fit);
        state[chainNum].trees[treeNum].collapseEmptyNodes(fit, nodeParams);
        // this could be combined with collapseEmptyNodes
        for (int32_t j = 0; j < static_cast<int32_t>(data.numPredictors); ++j)
          updateVariablesAvailable(fit, state[chainNum].trees[treeNum].top, j);
                
        state[chainNum].trees[treeNum].setCurrentFitsFromParameters(fit, nodeParams, treeFits, NULL);
        misc_addVectorsInPlace(treeFits, data.numObservations, 1.0, chainScratch[chainNum].totalFits);
        
        delete [] nodeParams;
      }
    }
  }
  
  bool updateTreesWithNewPredictor(BARTFit& fit) {
    const Control& control(fit.control);
    const Data& data(fit.data);
    State* state(fit.state);
    ChainScratch* chainScratch(fit.chainScratch);
    
    size_t totalNumTrees = control.numTrees * control.numChains;
    double** treeParams = new double*[totalNumTrees];
    for (size_t treeNum = 0; treeNum < totalNumTrees; ++treeNum)
      treeParams[treeNum] = NULL;
    
    bool allTreesAreValid = true;
    
    for (size_t chainNum = 0; chainNum < control.numChains && allTreesAreValid; ++chainNum) {
      
      for (size_t treeNum = 0; treeNum < control.numTrees && allTreesAreValid; ++treeNum) {
        const double* treeFits = state[chainNum].treeFits + treeNum * data.numObservations;
        
        // next allocates memory
        treeParams[treeNum + chainNum * control.numTrees] = 
          state[chainNum].trees[treeNum].recoverParametersFromFits(fit, treeFits);
        
        state[chainNum].trees[treeNum].top.addObservationsToChildren(fit);
        
        allTreesAreValid &= state[chainNum].trees[treeNum].isValid();
      }
    }
    
    if (!allTreesAreValid) goto updateTreesWithNewPredictor_cleanup;
    
    // go back across bottoms and set predictions to those mus for obs now in node
    for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum) {
      for (size_t treeNum = 0; treeNum < control.numTrees; ++treeNum) {
        double* treeFits = state[chainNum].treeFits + treeNum * data.numObservations;
        const double* nodeParams = treeParams[treeNum + chainNum * control.numTrees];
        
        misc_addVectorsInPlace(treeFits, data.numObservations, -1.0, chainScratch[chainNum].totalFits);
        
        state[chainNum].trees[treeNum].setCurrentFitsFromParameters(fit, nodeParams, treeFits, NULL);
        for (int32_t j = 0; j < static_cast<int32_t>(data.numPredictors); ++j)
          updateVariablesAvailable(fit, state[chainNum].trees[treeNum].top, j);
        
        misc_addVectorsInPlace(treeFits, data.numObservations, 1.0, chainScratch[chainNum].totalFits);
      }
    }
        
updateTreesWithNewPredictor_cleanup:
    for (size_t treeNum = totalNumTrees; treeNum > 0; --treeNum)
      delete [] treeParams[treeNum - 1];
    
    delete [] treeParams;
    
    return allTreesAreValid;
  }
}

namespace dbarts {
  
  bool BARTFit::setPredictor(const double* newPredictor, bool forceUpdate, bool updateCutPoints)
  {
    const double* oldPredictor = data.x;
    double** oldCutPoints = NULL;
    
    // cache old cutpoints if necessary
    if (!forceUpdate && updateCutPoints) {
      oldCutPoints = new double*[data.numPredictors];
      for (size_t j = 0; j < data.numPredictors; ++j) {
        oldCutPoints[j] = new double[numCutsPerVariable[j]];
        std::memcpy(oldCutPoints[j], cutPoints[j], numCutsPerVariable[j] * sizeof(double));
      }
    }
    
    data.x = newPredictor;
    
    if (updateCutPoints) {
      // find new cut points from default rule
      size_t* columns = misc_stackAllocate(data.numPredictors, size_t);
      for (size_t i = 0; i < data.numPredictors; ++i) columns[i] = i;
      
      ::setCutPoints(*this, columns, data.numPredictors);
      
      misc_stackFree(columns);
    }
    
    setXIntegerCutMap(*this);
    
    if (forceUpdate) {
      forceUpdateTrees(*this);
      if (data.numTestObservations > 0 && updateCutPoints)
        setXTestIntegerCutMap(*this, data.x_test, data.numTestObservations, const_cast<xint_t*>(sharedScratch.xt_test));
      return true;
    }
    
    bool treesAreValid = updateTreesWithNewPredictor(*this);
    
    // rollback
    if (!treesAreValid) {
      data.x = oldPredictor;
      
      if (updateCutPoints) for (size_t j = 0; j < data.numPredictors; ++j)
        std::memcpy(const_cast<double**>(cutPoints)[j], oldCutPoints[j], numCutsPerVariable[j] * sizeof(double));
      
      setXIntegerCutMap(*this);
            
      for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum) {
        for (size_t treeNum = 0; treeNum < control.numTrees; ++treeNum)
          state[chainNum].trees[treeNum].top.addObservationsToChildren(*this);
      }
    } else {
      if (updateCutPoints && data.numTestObservations > 0)
        setXTestIntegerCutMap(*this, data.x_test, data.numTestObservations, const_cast<xint_t*>(sharedScratch.xt_test));
    }
    
    if (updateCutPoints) {
      for (size_t j = data.numPredictors; j > 0; --j) delete [] oldCutPoints[j - 1];
      delete [] oldCutPoints;
    }
    
    return treesAreValid;
  }
  
  bool BARTFit::updatePredictor(const double* newPredictor, const size_t* columns, size_t numColumns, bool forceUpdate, bool updateCutPoints)
  {
    // store current
    double* oldPredictor = NULL;
    double** oldCutPoints = NULL;
    
    if (!forceUpdate) {
      oldPredictor = new double[data.numObservations * numColumns];
      for (size_t j = 0; j < numColumns; ++j)
        std::memcpy(oldPredictor + j * data.numObservations, data.x + columns[j] * data.numObservations,
                    data.numObservations * sizeof(double));
      
      if (updateCutPoints) {
        oldCutPoints = new double*[numColumns];
        for (size_t j = 0; j < numColumns; ++j) {
          oldCutPoints[j] = new double[numCutsPerVariable[columns[j]]];
          std::memcpy(oldCutPoints[j], cutPoints[columns[j]], numCutsPerVariable[columns[j]] * sizeof(double));
        }
      }
    }
    
    double* x = const_cast<double*>(data.x);
    for (size_t j = 0; j < numColumns; ++j) {
      size_t col = columns[j];
      std::memcpy(x + col * data.numObservations, newPredictor + j * data.numObservations, data.numObservations * sizeof(double));
    }
    
    // install new
    if (updateCutPoints)
      ::setCutPoints(*this, columns, numColumns);
    
    setXIntegerCutMap(*this, columns, numColumns);
    
    if (forceUpdate) {
      forceUpdateTrees(*this);
      if (data.numTestObservations > 0 && updateCutPoints)
        setXTestIntegerCutMap(*this, data.x_test, data.numTestObservations,
                              const_cast<xint_t*>(sharedScratch.xt_test), columns, numColumns);
      return true;
    }
    
    bool treesAreValid = updateTreesWithNewPredictor(*this);
    
    if (!treesAreValid) {
      for (size_t j = 0; j < numColumns; ++j) {
        size_t col = columns[j];
        
        std::memcpy(x + col * data.numObservations, oldPredictor + j * data.numObservations, data.numObservations * sizeof(double));
        if (updateCutPoints)
          std::memcpy(const_cast<double**>(cutPoints)[col], oldCutPoints[j], numCutsPerVariable[col] * sizeof(double));
      }
      setXIntegerCutMap(*this, columns, numColumns);
      
      for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum) {
        for (size_t treeNum = 0; treeNum < control.numTrees; ++treeNum)
          state[chainNum].trees[treeNum].top.addObservationsToChildren(*this);
      }
    } else {
      if (updateCutPoints && data.numTestObservations > 0)
        setXTestIntegerCutMap(*this, data.x_test, data.numTestObservations,
                              const_cast<xint_t*>(sharedScratch.xt_test), columns, numColumns);
    }
    
    if (updateCutPoints) {
      for (size_t j = numColumns; j > 0; --j) delete [] oldCutPoints[j - 1];
      delete [] oldCutPoints;
    }
    delete [] oldPredictor;
    
    return treesAreValid;
  }
  
  void BARTFit::setCutPoints(const double* const* newCutPoints, const uint32_t* numCutPoints,
                             const size_t* columns, size_t numColumns)
  {
    for (size_t j = 0; j < numColumns; ++j) {
      if (numCutPoints[j] > static_cast<uint32_t>(((1 >> sizeof(xint_t)) - 1))) {
        ext_throwError("number of cut points for column %lu exceeds the max allowable; "
                        "recompile package with -with-xint-size=32 or 64");
      }
    }
    
    // install the new cutpoints
    for (size_t j = 0; j < numColumns; ++j) {
      size_t col = columns[j];
      if (numCutsPerVariable[col] != numCutPoints[j]) {
        delete [] const_cast<double**>(cutPoints)[col];
        const_cast<double**>(cutPoints)[col] = new double[numCutPoints[j]];
        const_cast<uint32_t*>(numCutsPerVariable)[col] = numCutPoints[j];
        
        if (data.maxNumCuts[col] <= numCutPoints[j])
          const_cast<uint32_t*>(data.maxNumCuts)[col] = numCutPoints[j];
      }
        
      std::memcpy(const_cast<double**>(cutPoints)[col], newCutPoints[j], numCutsPerVariable[col] * sizeof(double));
    }
    setXIntegerCutMap(*this, columns, numColumns);
    if (data.numTestObservations > 0)
      setXTestIntegerCutMap(*this, data.x_test, data.numTestObservations,
                              const_cast<xint_t*>(sharedScratch.xt_test), columns, numColumns);
    
    forceUpdateTrees(*this);
  }
  
#define INVALID_ADDRESS reinterpret_cast<const double*>(this)
  void BARTFit::setTestPredictor(const double* newTestPredictor, size_t numTestObservations) {
    setTestPredictorAndOffset(newTestPredictor, INVALID_ADDRESS, numTestObservations);
  }
  
  void BARTFit::setTestOffset(const double* newTestOffset) {
     data.testOffset = newTestOffset;
  }
}

namespace {
  void updateTestFitsWithNewPredictor(const BARTFit& fit, ChainScratch* chainScratch) {
    const Control& control(fit.control);
    const Data& data(fit.data);
    const State* state(fit.state);
    
    double* currTestFits = new double[data.numTestObservations];
    
    for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum) {
      
      misc_setVectorToConstant(chainScratch[chainNum].totalTestFits, data.numTestObservations, 0.0);
      
      for (size_t treeNum = 0; treeNum < control.numTrees; ++treeNum) {
        const double* treeFits = state[chainNum].treeFits + treeNum * data.numObservations;
   
        const double* nodeParams = state[chainNum].trees[treeNum].recoverParametersFromFits(fit, treeFits);
      
        state[chainNum].trees[treeNum].setCurrentFitsFromParameters(fit, nodeParams, NULL, currTestFits);
      
        misc_addVectorsInPlace(currTestFits, data.numTestObservations, 1.0, chainScratch[chainNum].totalTestFits);
      
        delete [] nodeParams;
      }
    }
    
    delete [] currTestFits;
  }
}

namespace dbarts {
  // setting testOffset to NULL is valid
  // an invalid pointer address for testOffset is the object itself; when invalid, it is not updated
  void BARTFit::setTestPredictorAndOffset(const double* x_test, const double* testOffset, size_t numTestObservations) {
    if (numTestObservations == 0 || x_test == NULL) {
      if (sharedScratch.xt_test != NULL) { delete [] sharedScratch.xt_test; sharedScratch.xt_test = NULL; }
      for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum)
        if (chainScratch[chainNum].totalTestFits != NULL) { delete [] chainScratch[chainNum].totalTestFits; chainScratch[chainNum].totalTestFits = NULL; }
      
      data.x_test = NULL;
      data.numTestObservations = 0;
      data.testOffset = NULL;
    } else {
      data.x_test = x_test;
      
      if (numTestObservations != data.numTestObservations) {
        if (sharedScratch.xt_test != NULL) { delete [] sharedScratch.xt_test; sharedScratch.xt_test = NULL; }
        for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum)
          if (chainScratch[chainNum].totalTestFits != NULL) { delete [] chainScratch[chainNum].totalTestFits; chainScratch[chainNum].totalTestFits = NULL; }
        data.numTestObservations = numTestObservations;
        
        sharedScratch.xt_test = new xint_t[data.numTestObservations * data.numPredictors];
        for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum)
          chainScratch[chainNum].totalTestFits = new double[data.numTestObservations];
      }
      
      setXTestIntegerCutMap(*this, x_test, numTestObservations, const_cast<xint_t*>(sharedScratch.xt_test));
      // misc_transposeMatrix(data.x_test, data.numTestObservations, data.numPredictors, const_cast<double*>(sharedScratch.xt_test));
      
      if (testOffset != INVALID_ADDRESS) data.testOffset = testOffset;
      
      updateTestFitsWithNewPredictor(*this, chainScratch);
    }
  }
#undef INVALID_ADDRESS
  
  void BARTFit::updateTestPredictor(const double* newTestPredictor, size_t column) {
    updateTestPredictors(newTestPredictor, &column, 1);
  }
  
  void BARTFit::updateTestPredictors(const double* newTestPredictor, const size_t* columns, size_t numColumns) {
    double* x_test = const_cast<double*>(data.x_test);
    xint_t* xt_test = const_cast<xint_t*>(sharedScratch.xt_test);
    
    for (size_t j = 0; j < numColumns; ++j) {
      size_t col = columns[j];
      std::memcpy(x_test + col * data.numTestObservations, newTestPredictor + j * data.numTestObservations, data.numTestObservations * sizeof(double));
      
      for (size_t i = 0; i < data.numTestObservations; ++i) {
        xint_t k = 0;
        while (k < numCutsPerVariable[col] &&
               x_test[i + col * data.numTestObservations] > cutPoints[col][k]) ++k;
        xt_test[i * data.numPredictors + col] = k;
      }
    }
    
    updateTestFitsWithNewPredictor(*this, chainScratch);
  }
  
  void BARTFit::storeLatents(double* target) const {
    for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum)
      std::memcpy(target + chainNum * data.numObservations, chainScratch[chainNum].probitLatents, data.numObservations * sizeof(double));
  }
  
  /* to update data, we need to keep the scratch and the state sane
   * that means updating:
   *
   *   sharedScratch.yRescaled    - can just resize and copy in
   *   sharedScratch.xt           - same
   *   sharedScratch.xt_test      - same
   *   chainScratch.treeY         - just resize, is a temp array
   *   chainScratch.probitLatents - resize and initialize to new values
   *   sharedScratch.dataScale    - compute
   *   cutPoints and numCutsPerVariable - compute from data.maxNumCuts and new data
   *
   *   state.trees         - we have to go through these and prune any now-invalid end nodes
   *   state.treeIndices   - resize and assign into trees
   *   state.treeFits      - resize and recompute for new xt
   *   chainScratch.totalFits     - same
   *   chainScratch.totalTestFits - same
   *   state.sigma         - rescale using new scratch.dataScale
   */
  void BARTFit::setData(const Data& newData)
  {
    size_t oldNumObservations     = data.numObservations;
    size_t oldNumTestObservations = data.numTestObservations;
    
    data = newData;
    
    if (oldNumObservations != data.numObservations) {
      // handle resizing arrays
      delete [] sharedScratch.x;
      
      sharedScratch.x = new xint_t[data.numObservations * data.numPredictors];
      
      if (!control.responseIsBinary) {
        delete [] sharedScratch.yRescaled;
        sharedScratch.yRescaled = new double[data.numObservations];
      }
    }
    
    size_t** oldTreeIndices      = misc_stackAllocate(control.numChains, size_t*);
    double** oldTreeFits         = misc_stackAllocate(control.numChains, double*);
    
    double** currTestFits = misc_stackAllocate(control.numChains, double*);
    
    for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum) {
      // extract from old data what we'll need to update
      oldTreeIndices[chainNum]      = state[chainNum].treeIndices;
      oldTreeFits[chainNum]         = state[chainNum].treeFits;
      
      currTestFits[chainNum] = NULL;
      
      if (oldNumObservations != data.numObservations) {
        delete [] chainScratch[chainNum].totalFits;
        delete [] chainScratch[chainNum].treeY;
      
        chainScratch[chainNum].treeY     = new double[data.numObservations];
        chainScratch[chainNum].totalFits = new double[data.numObservations];
        
        if (control.responseIsBinary) {
          delete [] chainScratch[chainNum].probitLatents;
          chainScratch[chainNum].probitLatents = new double[data.numObservations];
        }
        
        state[chainNum].treeIndices = new size_t[data.numObservations * control.numTrees];
        state[chainNum].treeFits    = new double[data.numObservations * control.numTrees];
      }
    }
    
    // update sharedScratch.yRescaled/chainScratch.probitLatents and state.sigma
    if (control.responseIsBinary) {
      initializeLatents(*this);
    } else {
      double* sigmaUnscaled = misc_stackAllocate(control.numChains, double);
      for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum)
        sigmaUnscaled[chainNum] = state[chainNum].sigma * sharedScratch.dataScale.range;
      
      double priorUnscaled = model.sigmaSqPrior->getScale() * sharedScratch.dataScale.range * sharedScratch.dataScale.range;
      
      rescaleResponse(*this);
      
      model.sigmaSqPrior->setScale(priorUnscaled / (sharedScratch.dataScale.range * sharedScratch.dataScale.range));
      
      for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum)
        state[chainNum].sigma = sigmaUnscaled[chainNum] / sharedScratch.dataScale.range;
      
      misc_stackFree(sigmaUnscaled);
    }
           
    // cache old cut points, for use in updating trees
    const double** oldCutPoints = misc_stackAllocate(data.numPredictors, const double*);
    for (size_t j = 0; j < data.numPredictors; ++j) {
      oldCutPoints[j] = cutPoints[j];
      // next assignments 'reset' the variables, so setCutPoints() ignores old values
      const_cast<uint32_t*>(numCutsPerVariable)[j] = static_cast<uint32_t>(-1);
      const_cast<double**>(cutPoints)[j] = NULL;
    }
    // set new cut points
    size_t* columns = misc_stackAllocate(data.numPredictors, size_t);
    for (size_t j = 0; j < data.numPredictors; ++j) columns[j] = j;
    ::setCutPoints(*this, columns, data.numPredictors);
    misc_stackFree(columns);
    
    // now initialize remaining arrays that use numObs
    setXIntegerCutMap(*this);
    for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum)
      misc_setVectorToConstant(chainScratch[chainNum].totalFits, data.numObservations, 0.0);
    
    
    if (data.numTestObservations == 0 || data.x_test == NULL) {
      // no new data in test set
      delete [] sharedScratch.xt_test;
      sharedScratch.xt_test = NULL;
      for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum) {
        delete [] chainScratch[chainNum].totalTestFits;
        chainScratch[chainNum].totalTestFits = NULL;
      }
    } else {
      // handle resizing test arrays and initializing them
      if (oldNumTestObservations != data.numTestObservations) {
        delete [] sharedScratch.xt_test;
        sharedScratch.xt_test = new xint_t[data.numTestObservations * data.numPredictors];
        
        for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum) {
          delete [] chainScratch[chainNum].totalTestFits;
          chainScratch[chainNum].totalTestFits = new double[data.numTestObservations];
        }
      }
      setXTestIntegerCutMap(*this, data.x_test, data.numTestObservations, const_cast<xint_t*>(sharedScratch.xt_test));
      
      for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum) {
        currTestFits[chainNum] = new double[data.numTestObservations];
        misc_setVectorToConstant(chainScratch[chainNum].totalTestFits, data.numTestObservations, 0.0);
      }
    }
    
    // now update the trees, which is a bit messy
    for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum) {
      for (size_t treeNum = 0; treeNum < control.numTrees; ++treeNum) {
        const double* oldTreeFits_i = oldTreeFits[chainNum] + treeNum * oldNumObservations;
        
        // Use the bottom node enumeration to determine which fits to use.
        // The bottom nodes themselves keep an enumeration number, so that when prunned
        // we can still find the right one.
        state[chainNum].trees[treeNum].top.enumerateBottomNodes();
        
        // this allocates memory; predictions are of length equal to the number of bottom nodes
        double* nodeParams = state[chainNum].trees[treeNum].recoverParametersFromFits(*this, oldTreeFits_i);
        
        // the mapping can end up with some end-nodes that no longer exist, handles that internally
        state[chainNum].trees[treeNum].mapOldCutPointsOntoNew(*this, oldCutPoints, nodeParams);
        
        if (oldNumObservations != data.numObservations) {
          state[chainNum].trees[treeNum].top.observationIndices = state[chainNum].treeIndices + treeNum * data.numObservations;
          state[chainNum].trees[treeNum].top.numObservations = data.numObservations;
        }
        
        state[chainNum].trees[treeNum].top.addObservationsToChildren(*this);
        state[chainNum].trees[treeNum].collapseEmptyNodes(*this, nodeParams);
        for (int32_t i = 0; i < static_cast<int32_t>(data.numPredictors); ++i)
          updateVariablesAvailable(*this, state[chainNum].trees[treeNum].top, i);
        
        double* currTreeFits = state[chainNum].treeFits + treeNum * data.numObservations;
        state[chainNum].trees[treeNum].setCurrentFitsFromParameters(*this, nodeParams, currTreeFits, currTestFits[chainNum]);
        misc_addVectorsInPlace(currTreeFits, data.numObservations, 1.0, chainScratch[chainNum].totalFits);
        
        if (data.numTestObservations > 0)
          misc_addVectorsInPlace(currTestFits[chainNum], data.numTestObservations, 1.0, chainScratch[chainNum].totalTestFits);
        
        delete [] nodeParams;
      }
    }
    
    for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum)
      delete [] currTestFits[chainNum];
    
    for (size_t i = 0; i < data.numPredictors; ++i) delete [] oldCutPoints[i];
    misc_stackFree(oldCutPoints);
    
    if (oldNumObservations != data.numObservations) {
      for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum) {
        delete [] oldTreeFits[chainNum];
        delete [] oldTreeIndices[chainNum];
      }
    }
    
    misc_stackFree(currTestFits);
    
    misc_stackFree(oldTreeFits);
    misc_stackFree(oldTreeIndices);
  }
  
  void BARTFit::setControl(const Control& newControl)
  {
    bool stateResized = false;
    if (control.numChains == newControl.numChains) {
      for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum)
        stateResized |= state[chainNum].resize(*this, newControl);
    } else {
    
      size_t resizeEnd = std::min(control.numChains, newControl.numChains);
      
      State* oldState = state;
      state = static_cast<State*>(::operator new (newControl.numChains * sizeof(State)));
      
      for (size_t chainNum = 0; chainNum < resizeEnd; ++chainNum) {
        state[chainNum] = oldState[chainNum];
        stateResized |= state[chainNum].resize(*this, newControl);
      }
      
      for (size_t chainNum = resizeEnd; chainNum < newControl.numChains; ++chainNum) {
        new (state + chainNum) State(newControl, data);
        stateResized = true;
      }
      
      for (size_t chainNum = control.numChains; chainNum > resizeEnd; --chainNum)
        oldState[chainNum - 1].invalidate(control.numTrees, currentNumSamples);
      
      ::operator delete (oldState);
    }
    
    if (control.numTrees != newControl.numTrees) {
      NormalPrior* nodePrior = static_cast<NormalPrior*>(model.muPrior);
      double precisionUnscaled = nodePrior->precision / static_cast<double>(control.numTrees);
      nodePrior->precision = precisionUnscaled * static_cast<double>(newControl.numTrees);
    }
    
    rng_algorithm_t old_rng_algorithm = control.rng_algorithm;
    rng_standardNormal_t old_rng_standardNormal = control.rng_standardNormal;
    
    control = newControl;
    
    if (old_rng_algorithm != control.rng_algorithm || old_rng_standardNormal != control.rng_standardNormal) {
      destroyRNG(*this);
      createRNG(*this);
    }
    
    if (stateResized) {
      rebuildScratchFromState();
      currentSampleNum = 0;
    }
  }
  
  void BARTFit::setModel(const Model& newModel)
  {
    // double priorUnscaled = model.sigmaSqPrior->getScale() * sharedScratch.dataScale.range * sharedScratch.dataScale.range;
    
    model = newModel;
    
    // TODO: currently new model is assumed to be tweaked like model is internally,
    // which won't work for sigmasq priors with different specified quantiles or DoF
    
    // model.sigmaSqPrior->setScale(priorUnscaled / (sharedScratch.dataScale.range * sharedScratch.dataScale.range));
  }
  
  void BARTFit::printTrees(const size_t* chainIndices, size_t numChainIndices,
                           const size_t* sampleIndices, size_t numSampleIndices,
                           const size_t* treeIndices, size_t numTreeIndices) const
  {
    size_t indent = 0;
    
    for (size_t i = 0; i < numChainIndices; ++i) {
      size_t chainNum = chainIndices[i];
      if (numChainIndices > 1) {
        ext_printf("chain %lu\n", chainNum + 1);
        indent += 2;
      }
      if (!control.keepTrees) {
        for (size_t k = 0; k < numTreeIndices; ++k) {
          size_t treeNum = treeIndices[k];
          
          const double* treeFits = state[chainNum].treeFits + treeNum * data.numObservations;
          double* nodeParams = state[chainNum].trees[treeNum].recoverParametersFromFits(*this, treeFits);
          
          NodeVector bottomNodes(const_cast<Tree*>(&state[chainNum].trees[treeNum])->top.getBottomVector());
          size_t numBottomNodes = bottomNodes.size();
          for (size_t k = 0; k < numBottomNodes; ++k) bottomNodes[k]->setAverage(nodeParams[k]);
          delete [] nodeParams;
          
          state[chainNum].trees[treeNum].top.print(*this, indent);
        }
      } else for (size_t j = 0; j < numSampleIndices; ++j) {
        size_t sampleNum = sampleIndices[j];
        if (numSampleIndices > 1) {
          ext_printf("%*ssample %lu\n", indent, "", sampleNum + 1);
          indent += 2;
        }
        for (size_t k = 0; k < numTreeIndices; ++k) {
          size_t treeNum = treeIndices[k];
          
          size_t treeOffset = treeNum + sampleNum * control.numTrees;
          
          state[chainNum].savedTrees[treeOffset].top.print(*this, indent);
        }
        if (numSampleIndices > 1) indent -= 2;
      }
      if (numChainIndices > 1) indent -= 2;
    }
  }
}

namespace {
  using namespace dbarts;
  
  size_t storeFlattenedTree(const BARTFit& fit, const Node& node, size_t* numObservations,
                            int32_t* variable, double* value)
  {
    if (node.isBottom()) {
      *numObservations = node.numObservations;
      *variable = DBARTS_INVALID_RULE_VARIABLE;
      *value = node.m.average;
      return 1;
    }
    
    *numObservations = node.numObservations;
    *variable = node.p.rule.variableIndex;
    *value = fit.cutPoints[node.p.rule.variableIndex][node.p.rule.splitIndex];
    
    size_t numNodes = 1;
    numNodes += storeFlattenedTree(fit, *node.getLeftChild(), numObservations + numNodes,
                                   variable + numNodes, value + numNodes);
    numNodes += storeFlattenedTree(fit, *node.getRightChild(), numObservations + numNodes,
                                   variable + numNodes, value + numNodes);
    
    return numNodes;  
  }
  
  size_t storeFlattenedTree(const BARTFit& fit, const SavedNode& node, std::set<size_t>& indexSet,
                            size_t* numObservations, int32_t* variable, double* value)
  {
    if (node.isBottom()) {
      *numObservations = indexSet.size();
      *variable = DBARTS_INVALID_RULE_VARIABLE;
      *value = node.prediction;
      return 1;
    }
    
    *numObservations = indexSet.size();
    *variable = node.variableIndex;
    *value = node.split;
    
    std::set<size_t> leftIndexSet;
    std::set<size_t> rightIndexSet;
    
    for (std::set<size_t>::iterator it = indexSet.begin(); it != indexSet.end(); ++it) {
      size_t i = *it;
      
      if (fit.data.x[i + fit.data.numObservations * node.variableIndex] <= node.prediction) {
        rightIndexSet.insert(i);
      } else {
        leftIndexSet.insert(i); 
      }
    }
    
    size_t numNodes = 1;
    numNodes += storeFlattenedTree(fit, *node.getLeftChild(), leftIndexSet,
                                   numObservations + numNodes, variable + numNodes, value + numNodes);
    numNodes += storeFlattenedTree(fit, *node.getRightChild(), rightIndexSet,
                                   numObservations + numNodes, variable + numNodes, value + numNodes);
    
    return numNodes;  
  }
}

namespace dbarts {
  FlattenedTrees::FlattenedTrees(size_t totalNumNodes) : totalNumNodes(totalNumNodes), chainNumber(NULL),
    sampleNumber(NULL), treeNumber(NULL), numObservations(NULL), variable(NULL), value(NULL)
  {
    chainNumber     = new size_t[totalNumNodes];
    sampleNumber    = new size_t[totalNumNodes];
    treeNumber      = new size_t[totalNumNodes];
    numObservations = new size_t[totalNumNodes];
    variable = new int32_t[totalNumNodes];
    value = new double[totalNumNodes];
  }
  FlattenedTrees::~FlattenedTrees() {
    delete [] value;
    delete [] variable;
    delete [] numObservations;
    delete [] treeNumber;
    delete [] sampleNumber;
    delete [] chainNumber;
  }
  
  FlattenedTrees* BARTFit::getFlattenedTrees(const size_t* chainIndices, size_t numChainIndices,
                                             const size_t* sampleIndices, size_t numSampleIndices,
                                             const size_t* treeIndices, size_t numTreeIndices) const
  {
    
    size_t totalNumNodes = 0;
    
    // count how many nodes we're getting
    for (size_t i = 0; i < numChainIndices; ++i) {
      size_t chainNum = chainIndices[i];
      
      if (!control.keepTrees) { 
        for (size_t k = 0; k < numTreeIndices; ++k) {
          size_t treeNum = treeIndices[k];
          totalNumNodes += 1 + state[chainNum].trees[treeNum].top.getNumNodesBelow();
        }
      } else for (size_t j = 0; j < numSampleIndices; ++j) {
        size_t sampleNum = sampleIndices[j];
        for (size_t k = 0; k < numTreeIndices; ++k) {
          size_t treeNum = treeIndices[k];
          size_t treeOffset = treeNum + sampleNum * control.numTrees;
          totalNumNodes += 1 + state[chainNum].savedTrees[treeOffset].top.getNumNodesBelow();
        }
      }
    }
    
    FlattenedTrees* resultPtr = new FlattenedTrees(totalNumNodes);
    FlattenedTrees& result(*resultPtr);
    
    size_t offset = 0;
    for (size_t i = 0; i < numChainIndices; ++i) {
      size_t chainNum = chainIndices[i];
      
      if (!control.keepTrees) { 
        for (size_t k = 0; k < numTreeIndices; ++k) {
          size_t treeNum = treeIndices[k];
          
          // get nodes to store averages
          const double* treeFits = state[chainNum].treeFits + treeNum * data.numObservations;
          double* nodeParams = state[chainNum].trees[treeNum].recoverParametersFromFits(*this, treeFits);
          
          NodeVector bottomNodes(const_cast<Tree*>(&state[chainNum].trees[treeNum])->top.getBottomVector());
          size_t numBottomNodes = bottomNodes.size();
          for (size_t k = 0; k < numBottomNodes; ++k) bottomNodes[k]->setAverage(nodeParams[k]);
          delete [] nodeParams;
          
          
          size_t numNodesInTree = storeFlattenedTree(*this, state[chainNum].trees[treeNum].top,
                                                     result.numObservations + offset,
                                                     result.variable + offset,
                                                     result.value + offset);
          
          for (size_t l = 0; l < numNodesInTree; ++l) {
            result.chainNumber [offset + l] = chainNum;
            result.sampleNumber[offset + l] = 0;
            result.treeNumber  [offset + l] = treeNum;
          }
          offset += numNodesInTree;
        }
      } else {
        std::set<size_t> indexSet;
        for (size_t i = 0; i < data.numObservations; ++i) indexSet.insert(i);
        
        for (size_t j = 0; j < numSampleIndices; ++j) {
          size_t sampleNum = sampleIndices[j];
          for (size_t k = 0; k < numTreeIndices; ++k) {
            size_t treeNum = treeIndices[k];
            size_t treeOffset = treeNum + sampleNum * control.numTrees;
            
            size_t numNodesInTree = storeFlattenedTree(*this, state[chainNum].savedTrees[treeOffset].top,
                                                       indexSet,
                                                       result.numObservations + offset,
                                                       result.variable + offset,
                                                       result.value + offset);
            for (size_t l = 0; l < numNodesInTree; ++l) {
              result.chainNumber [offset + l] = chainNum;
              result.sampleNumber[offset + l] = sampleNum;
              result.treeNumber  [offset + l] = treeNum;
            }
            offset += numNodesInTree;
          }
        }
      }
    }
    
    return resultPtr;
  }
  
  BARTFit::BARTFit(Control control, Model model, Data data) :
    control(control), model(model), data(data), sharedScratch(), chainScratch(NULL), state(NULL),
    runningTime(0.0), currentNumSamples(control.defaultNumSamples), currentSampleNum(0), threadManager(NULL)
  {
    allocateMemory(*this);
    
    if (control.responseIsBinary) initializeLatents(*this);
    else rescaleResponse(*this);

    createRNG(*this);
    
    setPrior(*this);
    setInitialCutPoints(*this);
    setXIntegerCutMap(*this);
    if (data.numTestObservations > 0)
      setXTestIntegerCutMap(*this, data.x_test, data.numTestObservations, const_cast<xint_t*>(sharedScratch.xt_test));
    setInitialFit(*this);

    if (this->control.verbose) printInitialSummary(*this);
  }
  
  BARTFit::~BARTFit()
  {
    destroyRNG(*this);
    
    delete [] sharedScratch.yRescaled; sharedScratch.yRescaled = NULL;
    delete [] sharedScratch.x; sharedScratch.x = NULL;
    delete [] sharedScratch.xt_test; sharedScratch.xt_test = NULL;
    for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum) {
      delete [] chainScratch[chainNum].totalTestFits; chainScratch[chainNum].totalTestFits = NULL;
      delete [] chainScratch[chainNum].totalFits; chainScratch[chainNum].totalFits = NULL;
      delete [] chainScratch[chainNum].probitLatents; chainScratch[chainNum].probitLatents = NULL;
      delete [] chainScratch[chainNum].treeY; chainScratch[chainNum].treeY = NULL;
    }
    
    delete [] chainScratch;
    
    delete [] numCutsPerVariable; numCutsPerVariable = NULL;
    if (cutPoints != NULL) {
      for (size_t i = 0; i < data.numPredictors; ++i) delete [] cutPoints[i];
    }
    delete [] cutPoints; cutPoints = NULL;
    
    for (size_t chainNum = control.numChains; chainNum > 0; --chainNum)
      state[chainNum - 1].invalidate(control.numTrees, currentNumSamples);
    
    ::operator delete (state);
    
    misc_htm_destroy(threadManager);
  }
  
  Results* BARTFit::runSampler()
  {
    // ensure at least one sample for state's sake
    Results* resultsPointer = new Results(data.numObservations, data.numPredictors,
                                          data.numTestObservations,
                                          control.defaultNumSamples == 0 ? 1 : control.defaultNumSamples,
                                          control.numChains,
                                          model.kPrior != NULL);
    size_t numBurnIn = control.defaultNumBurnIn - (control.defaultNumSamples == 0 && control.defaultNumBurnIn > 0 ? 1 : 0);
    
    runSampler(numBurnIn, resultsPointer);
    
    if (control.defaultNumSamples == 0) {
      delete resultsPointer;
      return NULL;
    }
    
    return resultsPointer;
  }
  
  Results* BARTFit::runSampler(size_t numBurnIn, size_t numSamples)
  {
    Results* resultsPointer = new Results(data.numObservations, data.numPredictors,
                                          data.numTestObservations,
                                          numSamples == 0 ? 1 : numSamples,
                                          control.numChains,
                                          model.kPrior != NULL);
    numBurnIn -= numSamples == 0 && numBurnIn > 0 ? 1 : 0;
    
    runSampler(numBurnIn, resultsPointer);
    
    if (numSamples == 0) {
      delete resultsPointer;
      return NULL;
    }
    
    return resultsPointer;
  }
  
    
  void BARTFit::sampleTreesFromPrior()
  {
    for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum) {
      for (size_t treeNum = 0; treeNum < control.numTrees; ++treeNum) {
        state[chainNum].trees[treeNum].sampleFromPrior(*this, state[chainNum].rng);
      }
    }
  }
}

extern "C" {
  struct ThreadData {
    BARTFit* fit;
    size_t chainNum;
    size_t numBurnIn;
    Results* results;
  };
  
  void samplerThreadFunction(std::size_t taskId, void* threadDataPtr) {
    ThreadData* threadData(reinterpret_cast<ThreadData*>(threadDataPtr));
    
    BARTFit& fit(*threadData->fit);
    size_t chainNum = threadData->chainNum;
    size_t numBurnIn = threadData->numBurnIn;
    Results& results(*threadData->results);
    
    Control& control(fit.control);
    Model& model(fit.model);
    Data& data(fit.data);
    
    SharedScratch& sharedScratch(fit.sharedScratch);
    ChainScratch& chainScratch(fit.chainScratch[chainNum]);
    State& state(fit.state[chainNum]);
    
    chainScratch.taskId = taskId;
    
    bool stepTaken;
    StepType ignored;
    
    size_t numSamples = results.numSamples;
    
    double* currFits = new double[data.numObservations];
    double* currTestFits = data.numTestObservations > 0 ? new double[data.numTestObservations] : NULL;
    
    uint32_t* variableCounts = misc_stackAllocate(data.numPredictors, uint32_t);
    
    size_t totalNumIterations = (numBurnIn + numSamples) * control.treeThinningRate;
    
    // reserve once at the start if possible
    if (control.numThreads > 1 && control.numChains == 1)
      misc_htm_reserveThreadsForSubTask(fit.threadManager, 0, 0);
    
    // const cast b/c yRescaled doesn't change, but probit latents do
    double* y = control.responseIsBinary ? chainScratch.probitLatents : const_cast<double*>(sharedScratch.yRescaled);
    
    for (size_t k = 0; k < totalNumIterations; ++k) {
      if (control.numThreads > 1 && control.numChains > 1)
        misc_htm_reserveThreadsForSubTask(fit.threadManager, taskId, k);
      
      bool isThinningIteration = ((k + 1) % control.treeThinningRate != 0);
      size_t majorIterationNum = k / control.treeThinningRate;
     
      bool isBurningIn = majorIterationNum < numBurnIn;
      size_t resultSampleNum = !isBurningIn ? majorIterationNum - numBurnIn : 0;
            
      if (control.verbose && !isThinningIteration && (majorIterationNum + 1) % control.printEvery == 0) {
        if (control.numChains > 1)
          misc_htm_printf(fit.threadManager, "[%lu] iteration: %u (of %u)\n", chainNum + 1, k + 1, totalNumIterations);
        else
          ext_printf("iteration: %u (of %u)\n", k + 1, totalNumIterations);
      }
      
      if (!isThinningIteration && data.numTestObservations > 0) misc_setVectorToConstant(chainScratch.totalTestFits, data.numTestObservations, 0.0);
      
      for (size_t treeNum = 0; treeNum < control.numTrees; ++treeNum) {
        double* oldTreeFits = state.treeFits + treeNum * data.numObservations;
        
        // treeY = y - (totalFits - oldTreeFits)
        // is residual from every *other* tree, so what is left for this tree to do
        std::memcpy(chainScratch.treeY, y, data.numObservations * sizeof(double));
        misc_addVectorsInPlace(const_cast<const double*>(chainScratch.totalFits), data.numObservations, -1.0, chainScratch.treeY);
        misc_addVectorsInPlace(const_cast<const double*>(oldTreeFits), data.numObservations, 1.0, chainScratch.treeY);
        
        state.trees[treeNum].setNodeAverages(fit, chainNum, chainScratch.treeY);
        
        metropolisJumpForTree(fit, chainNum, state.trees[treeNum], chainScratch.treeY, state.sigma, &stepTaken, &ignored);
        
        state.trees[treeNum].sampleParametersAndSetFits(fit, chainNum, state.sigma, currFits, isThinningIteration ? NULL : currTestFits);
        
        // totalFits += currFits - treeFits
        misc_addVectorsInPlace(const_cast<const double*>(oldTreeFits), data.numObservations, -1.0, chainScratch.totalFits);
        misc_addVectorsInPlace(const_cast<const double*>(currFits), data.numObservations, 1.0, chainScratch.totalFits);
        
        if (!isThinningIteration && data.numTestObservations > 0)
          misc_addVectorsInPlace(const_cast<const double*>(currTestFits), data.numTestObservations, 1.0, chainScratch.totalTestFits);
        
        std::memcpy(oldTreeFits, const_cast<const double*>(currFits), data.numObservations * sizeof(double));
      }
      
      if (control.keepTrees & !isBurningIn && !isThinningIteration) {
        size_t treeSampleNum = (fit.currentSampleNum + resultSampleNum) % fit.currentNumSamples;
        
        for (size_t treeNum = 0; treeNum < control.numTrees; ++treeNum)
          state.savedTrees[treeNum + treeSampleNum * control.numTrees].copyStructureFrom(fit, state.trees[treeNum], const_cast<const double*>(state.treeFits + treeNum * data.numObservations));
      }
      
      if (control.responseIsBinary)
        sampleProbitLatentVariables(fit, state, chainScratch.totalFits, y);
      
      state.sigma = std::sqrt(model.sigmaSqPrior->drawFromPosterior(fit, chainNum, y, chainScratch.totalFits));
      if (model.kPrior != NULL) {
        state.k = model.kPrior->drawFromPosterior(fit, chainNum);
        model.muPrior->setK(fit, state.k);
      }
      
            
      if (!isThinningIteration) {
        // if not out of burn-in, store result in first result; start
        // overwriting after that
        for (size_t j = 0; j < fit.data.numPredictors; ++j) variableCounts[j] = 0;
        countVariableUses(fit, state, variableCounts);
        
        storeSamples(fit, chainNum, results, chainScratch.totalFits, chainScratch.totalTestFits, state.sigma,
                     state.k, variableCounts, resultSampleNum);
          
        if (control.callback != NULL) {
          size_t chainStride = chainNum * numSamples;
          control.callback(control.callbackData, fit, isBurningIn,
                           results.trainingSamples + (resultSampleNum + chainStride) * data.numObservations,
                           results.testSamples + (resultSampleNum + chainStride) * data.numTestObservations,
                           results.sigmaSamples[resultSampleNum + chainStride],
                           model.kPrior != NULL ? results.kSamples + resultSampleNum + chainStride : NULL);
        }
      }
    }
    
    delete [] currFits;
    if (data.numTestObservations > 0) delete [] currTestFits;
    misc_stackFree(variableCounts);
  }
}

namespace dbarts {
  
  void BARTFit::runSampler(size_t numBurnIn, Results* resultsPointer)
  {
    if (control.verbose) ext_printf("Running mcmc loop:\n");
    
#ifdef HAVE_SYS_TIME_H
    struct timeval startTime;
    struct timeval endTime;
    gettimeofday(&startTime, NULL);
#else
    time_t startTime;
    time_t endTime;
    startTime = time(NULL);
#endif
    
    if (control.keepTrees && currentNumSamples == 0) {
      currentNumSamples = 1;
      for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum)
        state[chainNum].resize(*this, currentNumSamples);
      currentSampleNum = 0;
    }
    
    if (control.numThreads <= 1) {
      // run single threaded, chains in sequence
      ThreadData threadData = { this, 0, numBurnIn, resultsPointer };
      for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum) {
        threadData.chainNum = chainNum;
        samplerThreadFunction(static_cast<size_t>(-1), reinterpret_cast<void*>(&threadData));
      }
    } else {
      ThreadData* threadData = new ThreadData[control.numChains];
      void** threadDataPtr = new void*[control.numChains];
      
      for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum) {
        threadData[chainNum].fit = this;
        threadData[chainNum].chainNum = chainNum;
        threadData[chainNum].numBurnIn = numBurnIn;
        threadData[chainNum].results = resultsPointer;
        threadDataPtr[chainNum] = reinterpret_cast<void*>(&threadData[chainNum]);
      }
      
      if (control.verbose) {
        struct timespec outputDelay;
        outputDelay.tv_sec = 0;
        outputDelay.tv_nsec = 100000000; // every 0.1 seconds
        misc_htm_runTopLevelTasksWithOutput(threadManager, &samplerThreadFunction, threadDataPtr, control.numChains, &outputDelay);
      } else {
        misc_htm_runTopLevelTasks(threadManager, &samplerThreadFunction, threadDataPtr, control.numChains);
      }
      
      delete [] threadDataPtr;
      delete [] threadData;
    }
    
    if (control.keepTrees)
      currentSampleNum = (currentSampleNum + resultsPointer->numSamples) % currentNumSamples;
    
#ifdef HAVE_SYS_TIME_H
    gettimeofday(&endTime, NULL);
#else
    endTime = time(NULL);
#endif
    
    runningTime += subtractTimes(endTime, startTime);
    
    if (control.verbose) printTerminalSummary(*this);
  }
} // namespace dbarts


namespace {
  using namespace dbarts;
  
  void printInitialSummary(const BARTFit& fit) {
    const Control& control(fit.control);
    const Data& data(fit.data);
    const Model& model(fit.model);
    
    const SharedScratch& sharedScratch(fit.sharedScratch);
    
    if (control.responseIsBinary)
      ext_printf("\nRunning BART with binary y\n\n");
    else
      ext_printf("\nRunning BART with numeric y\n\n");
    
    ext_printf("number of trees: %lu\n", control.numTrees);
    ext_printf("number of chains: %lu, number of threads %lu\n", control.numChains, control.numThreads);
    
    ext_printf("Prior:\n");
    // dirty hack... should have priors print themselves
    ext_printf("\tk: %f\n", model.muPrior->getK(fit));
    if (model.kPrior != NULL) model.kPrior->print(fit);
    if (!control.responseIsBinary)
      model.sigmaSqPrior->print(fit);
    
    CGMPrior* treePrior = static_cast<CGMPrior*>(model.treePrior);
    ext_printf("\tpower and base for tree prior: %f %f\n", treePrior->power, treePrior->base);
    ext_printf("\tuse quantiles for rule cut points: %s\n", control.useQuantiles ? "true" : "false");
    ext_printf("data:\n");
    ext_printf("\tnumber of training observations: %u\n", data.numObservations);
    ext_printf("\tnumber of test observations: %u\n", data.numTestObservations);
    ext_printf("\tnumber of explanatory variables: %u\n", data.numPredictors);
    if (!control.responseIsBinary) ext_printf("\tinit sigma: %f, curr sigma: %f\n", data.sigmaEstimate, fit.state[0].sigma * sharedScratch.dataScale.range);
    if (data.weights != NULL) ext_printf("\tusing observation weights\n");
    ext_printf("\n");
    
    
    ext_printf("Cutoff rules c in x<=c vs x>c\n");
    ext_printf("Number of cutoffs: (var: number of possible c):\n");
    for (size_t j = 0; j < data.numPredictors; ++j) {
      ext_printf("(%u: %u) ", j + 1, fit.numCutsPerVariable[j]);
      if ((j + 1) % 5 == 0) ext_printf("\n");
    }
    ext_printf("\n");
    if (control.printCutoffs > 0) {
      ext_printf("cutoffs:\n");
      for (size_t j = 0; j < data.numPredictors; ++j) {
        ext_printf("x(%u) cutoffs: ", j + 1);
        
        size_t k;
        for (k = 0; k < fit.numCutsPerVariable[j] - 1 && k < control.printCutoffs - 1; ++k) {
          ext_printf("%f", fit.cutPoints[j][k]);
          if ((k + 1) % 5 == 0) ext_printf("\n\t");
        }
        if (k > 2 && k == control.printCutoffs && k < fit.numCutsPerVariable[j] - 1)
          ext_printf("...");
        
        ext_printf("%f", fit.cutPoints[j][fit.numCutsPerVariable[j] - 1]);
        ext_printf("\n");
      }
    }
    
    if (data.offset != NULL || (data.numTestObservations > 0 && data.testOffset != NULL)) {
      ext_printf("\noffsets:\n");
      
      if (data.offset != NULL) {
        ext_printf("\treg : %.2f", data.offset[0]);
        for (size_t i = 1; i < (5 < data.numObservations ? 5 : data.numObservations); ++i) ext_printf(" %.2f", data.offset[i]);
        ext_printf("\n");
      }
      if (data.numTestObservations > 0 && data.testOffset != NULL) {
        ext_printf("\ttest: %.2f", data.testOffset[0]);
        for (size_t i = 1; i < (5 < data.numTestObservations ? 5 : data.numTestObservations); ++i) ext_printf(" %.2f", data.testOffset[i]);
        ext_printf("\n");
      }
    }
  }
  
  void printTerminalSummary(const BARTFit& fit) {
    ext_printf("total seconds in loop: %f\n", fit.runningTime);
    
    ext_printf("\nTree sizes, last iteration:\n");
    for (size_t chainNum = 0; chainNum < fit.control.numChains; ++chainNum) {
      size_t linePrintCount = 0;
      if (fit.control.numChains > 0) {
        ext_printf("[%u] ", chainNum + 1);
        linePrintCount += 2;
      }
      for (size_t treeNum = 0; treeNum < fit.control.numTrees; ++treeNum) {
        ext_printf("%u ", fit.state[chainNum].trees[treeNum].getNumBottomNodes());
        if ((linePrintCount++ + 1) % 20 == 0) ext_printf("\n");
      }
      if ((linePrintCount % 20) != 0) ext_printf("\n");
    }
    ext_printf("\n");
    
    uint32_t* variableCounts = misc_stackAllocate(fit.data.numPredictors, uint32_t);
    
    ext_printf("Variable Usage, last iteration (var:count):\n");
    for (size_t j = 0; j < fit.data.numPredictors; ++j) variableCounts[j] = 0;
    for (size_t chainNum = 0; chainNum < fit.control.numChains; ++chainNum)
      countVariableUses(fit, fit.state[chainNum], variableCounts);
    for (size_t j = 0; j < fit.data.numPredictors; ++j) {
      ext_printf("(%lu: %u) ", static_cast<unsigned long int>(j + 1), variableCounts[j]);
      if ((j + 1) % 5 == 0) ext_printf("\n");
    }
    
    misc_stackFree(variableCounts);
    
    ext_printf("\nDONE BART\n\n");
  }
  
  void allocateMemory(BARTFit& fit) {
    Control& control(fit.control);
    Data& data(fit.data);
    SharedScratch& sharedScratch(fit.sharedScratch);
    
    fit.chainScratch = new ChainScratch[control.numChains];
    ChainScratch* chainScratch = fit.chainScratch;
    
    if (!control.responseIsBinary) {
      sharedScratch.yRescaled = new double[data.numObservations];
      for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum)
        chainScratch[chainNum].probitLatents = NULL;
    } else {
      sharedScratch.yRescaled = NULL;
      for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum)
        chainScratch[chainNum].probitLatents = new double[data.numObservations];
    }
    
    // chain scratches
    for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum) {
      chainScratch[chainNum].treeY = new double[data.numObservations];
      double* y = control.responseIsBinary ? chainScratch[chainNum].probitLatents : const_cast<double*>(sharedScratch.yRescaled);
      
      for (size_t i = 0; i < data.numObservations; ++i) chainScratch[chainNum].treeY[i] = y[i];
      
      chainScratch[chainNum].totalFits = new double[data.numObservations];
      chainScratch[chainNum].totalTestFits = data.numTestObservations > 0 ? new double[data.numTestObservations] : NULL;
      
      chainScratch[chainNum].taskId = static_cast<size_t>(-1);
    }
    
    // shared scratch
    sharedScratch.x = new xint_t[data.numObservations * data.numPredictors];
    if (data.numTestObservations > 0)
      sharedScratch.xt_test = new xint_t[data.numTestObservations * data.numPredictors];
    
    fit.numCutsPerVariable = new uint32_t[data.numPredictors];

    fit.cutPoints = new double*[data.numPredictors];
    const double** cutPoints = const_cast<const double**>(fit.cutPoints);
    for (size_t j = 0; j < data.numPredictors; ++j) cutPoints[j] = NULL;
    
    // states
    fit.state = static_cast<State*>(::operator new (control.numChains * sizeof(State)));
    for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum)
      new (fit.state + chainNum) State(control, data);
    
    if (control.numThreads > 1 && misc_htm_create(&fit.threadManager, control.numThreads) != 0) {
      ext_printMessage("Unable to multi-thread, defaulting to single.");
      control.numThreads = 1;
    }
  }
  
  void setPrior(BARTFit& fit) {
    Control& control(fit.control);
    Data& data(fit.data);
    Model& model(fit.model);
    SharedScratch& sharedScratch(fit.sharedScratch);
    State* state(fit.state);
    
    if (control.responseIsBinary) {
      for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum) {
        state[chainNum].sigma = 1.0;
        state[chainNum].k = model.kPrior != NULL ? model.muPrior->getK(fit) : std::numeric_limits<double>::quiet_NaN();
      }
    } else {
      for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum) {
        state[chainNum].sigma = data.sigmaEstimate / sharedScratch.dataScale.range;
        state[chainNum].k = model.kPrior != NULL ? model.muPrior->getK(fit) : std::numeric_limits<double>::quiet_NaN();
      }
      model.sigmaSqPrior->setScale(state[0].sigma * state[0].sigma * model.sigmaSqPrior->getScale());
    }
  }
  
  void setInitialCutPoints(BARTFit& fit) {
    Data& data(fit.data);
    
    uint32_t* numCutsPerVariable = const_cast<uint32_t*>(fit.numCutsPerVariable);
    double** cutPoints = const_cast<double**>(fit.cutPoints);
    for (size_t i = 0; i < data.numPredictors; ++i) {
      numCutsPerVariable[i] = static_cast<uint32_t>(-1);
      cutPoints[i] = NULL;
    }
    
    size_t* columns = misc_stackAllocate(data.numPredictors, size_t);
    for (size_t j = 0; j < data.numPredictors; ++j) columns[j] = j;
    
    ::setCutPoints(fit, columns, data.numPredictors);
    
    misc_stackFree(columns);
  }
  
  void setXIntegerCutMap(BARTFit& fit)
  {
    const Data& data(fit.data);
    
    xint_t* x = const_cast<xint_t*>(fit.sharedScratch.x);
    
    for (size_t j = 0; j < data.numPredictors; ++j) {
      for (size_t i = 0; i < data.numObservations; ++i) {
        
        xint_t k = 0;
        
        // min cut such that x_ij <= c_jk; can possibly be out of range if variable is on the far right
        while (k < fit.numCutsPerVariable[j] &&
               data.x[i + j * data.numObservations] > fit.cutPoints[j][k]) ++k;
        
        x[i + j * data.numObservations] = k;
      }
    }
  }
  
  void setXIntegerCutMap(BARTFit& fit, const size_t* columns, size_t numColumns)
  {
    const Data& data(fit.data);
    
    xint_t* x = const_cast<xint_t*>(fit.sharedScratch.x);
    
    for (size_t j = 0; j < numColumns; ++j) {
      size_t col = columns[j];
      for (size_t i = 0; i < data.numObservations; ++i) {
        xint_t k = 0;
        while (k < fit.numCutsPerVariable[col] &&
               data.x[i + col * data.numObservations] > fit.cutPoints[col][k]) ++k;
        
        x[i + col * data.numObservations] = k;
      }
    }
  }
  
  void setXTestIntegerCutMap(const BARTFit& fit, const double* x_test, size_t numTestObservations, xint_t* xt_test)
  {
    const Data& data(fit.data);
    
    for (size_t j = 0; j < data.numPredictors; ++j) {
      for (size_t i = 0; i < numTestObservations; ++i) {
        xint_t k = 0;
        
        while (k < fit.numCutsPerVariable[j] &&
               x_test[i + j * numTestObservations] > fit.cutPoints[j][k]) ++k;
      
        xt_test[i * data.numPredictors + j] = k;
      }
    }
  }
  
  void setXTestIntegerCutMap(const BARTFit& fit, const double* x_test, size_t numTestObservations,
                             xint_t* xt_test, const size_t* columns, size_t numColumns)
  {
    const Data& data(fit.data);
    
    for (size_t j = 0; j < numColumns; ++j) {
      for (size_t i = 0; i < numTestObservations; ++i) {
        size_t col = columns[j];
        
        xint_t k = 0;
        
        while (k < fit.numCutsPerVariable[col] &&
               x_test[i + col * numTestObservations] > fit.cutPoints[col][k]) ++k;
      
        xt_test[i * data.numPredictors + col] = k;
      }
    }
  }
  
  void setCutPoints(BARTFit& fit, const size_t* columns, size_t numColumns)
  {
    Control& control(fit.control);
    Data& data(fit.data);
    
    uint32_t* numCutsPerVariable = const_cast<uint32_t*>(fit.numCutsPerVariable);
    double** cutPoints = const_cast<double**>(fit.cutPoints);
        
    if (control.useQuantiles) {
      if (data.maxNumCuts == NULL) ext_throwError("Num cuts cannot be NULL if useQuantiles is true.");
      
       // sets are inherently sorted, should be a binary tree back there somewhere
      std::set<double> uniqueElements;
      std::vector<double> sortedElements(data.numObservations);
      
      for (size_t j = 0; j < numColumns; ++j) {
        size_t col = columns[j];
        
        setCutPointsFromQuantiles(fit, data.x + col * data.numObservations, data.maxNumCuts[col],
                                  numCutsPerVariable[col], cutPoints[col],
                                  uniqueElements, sortedElements);
      }
    } else {
      for (size_t j = 0; j < numColumns; ++j) {
        size_t col = columns[j];
        
        setCutPointsUniformly(fit, data.x + col * data.numObservations, data.maxNumCuts[col],
                              numCutsPerVariable[col], cutPoints[col]);
      }
    }
  }
  
  void setCutPointsFromQuantiles(BARTFit& fit, const double* x, uint32_t maxNumCuts,
                                 uint32_t& numCutsPerVariable, double*& cutPoints,
                                 std::set<double>& uniqueElements, std::vector<double>& sortedElements)
  {
    Data& data(fit.data);
    
    // sets are inherently sorted, should be a binary tree back there somewhere
    uniqueElements.clear();
    for (size_t i = 0; i < data.numObservations; ++i) uniqueElements.insert(x[i]);
    
    size_t numUniqueElements = uniqueElements.size();
      
    size_t step, numCuts, offset;
    if (numUniqueElements <= maxNumCuts + 1) {
      step = 1;
      numCuts = numUniqueElements - 1;
      offset = 0;
    } else {
      numCuts = maxNumCuts;
      step = numUniqueElements / numCuts;
      offset = step / 2;
    }
    
    if (numCutsPerVariable != static_cast<uint32_t>(-1)) {
      if (numCuts < numCutsPerVariable) ext_throwError("number of induced cut points in new predictor less than previous: old splits would be invalid");
      if (numCuts > numCutsPerVariable) ext_issueWarning("number of induced cut points in new predictor greater than previous: ignoring extra quantiles");
    } else {
      numCutsPerVariable = static_cast<uint32_t>(numCuts);
      cutPoints = new double[numCuts];
    }
    
    sortedElements.clear();
    sortedElements.assign(uniqueElements.begin(), uniqueElements.end());
      
    for (size_t k = 0; k < numCutsPerVariable; ++k) {
      size_t index = std::min(k * step + offset, numUniqueElements - 2);
      cutPoints[k] = 0.5 * (sortedElements[index] + sortedElements[index + 1]);
    }
  }
  
  void setCutPointsUniformly(BARTFit& fit, const double* x, uint32_t maxNumCuts,
                             uint32_t& numCutsPerVariable, double*& cutPoints)
  {
    Data& data(fit.data);
    
    double xMax, xMin, xIncrement;
    
    xMax = x[0]; xMin = x[0];
    for (size_t i = 1; i < data.numObservations; ++i) {
      double x_i = x[i];
      if (x_i < xMin) xMin = x_i;
      if (x_i > xMax) xMax = x_i;
    }
    
    if (numCutsPerVariable == static_cast<uint32_t>(-1)) {
      numCutsPerVariable = maxNumCuts;
      cutPoints = new double[numCutsPerVariable];
    }
      
    xIncrement = (xMax - xMin) / static_cast<double>(numCutsPerVariable + 1);
      
    for (size_t k = 0; k < numCutsPerVariable; ++k) cutPoints[k] = xMin + (static_cast<double>(k + 1)) * xIncrement;
  }
    
  void createRNG(BARTFit& fit) {
    Control& control(fit.control);
    State* state(fit.state);
    
    // if only one chain or one thread, can use environment's rng since randomization calls will all be
    // serial
    bool useNativeRNG = control.numChains == 1 || control.numThreads == 1;
    
    size_t numSeedResets;
    size_t chainNum;
    const char* errorMessage = NULL;
    
    ext_rng_algorithm_t rng_algorithm = static_cast<ext_rng_algorithm_t>(control.rng_algorithm);
    ext_rng_standardNormal_t rng_standardNormal = static_cast<ext_rng_standardNormal_t>(control.rng_standardNormal);
    
    for (chainNum = 0; chainNum < control.numChains; ++chainNum) {
      if (rng_algorithm == EXT_RNG_ALGORITHM_INVALID) { // use default of some kind
        if ((state[chainNum].rng = ext_rng_createDefault(useNativeRNG)) == NULL) {
          errorMessage = "could not allocate rng";
          goto createRNG_cleanup;
        }
        
        if (rng_standardNormal != EXT_RNG_STANDARD_NORMAL_INVALID &&
            rng_standardNormal != EXT_RNG_STANDARD_NORMAL_USER_NORM &&
            ext_rng_setStandardNormalAlgorithm(state[chainNum].rng, rng_standardNormal, NULL) != 0) {
          errorMessage = "could not set rng standard normal";
          goto createRNG_cleanup;
        }
        // if not using envirnoment's rng, we have to seed
        if (!useNativeRNG) {
          if (ext_rng_setSeedFromClock(state[chainNum].rng) != 0) { errorMessage = "could not seed rng"; goto createRNG_cleanup; }
          
          if (chainNum > 0) {
            // check that seed is unique
            for (numSeedResets = 0; numSeedResets < static_cast<size_t>(-1); ++numSeedResets) {
              if (!ext_rng_seedsAreEqual(state[chainNum].rng, state[chainNum - 1].rng)) break;
              
              if (ext_rng_setSeedFromClock(state[chainNum].rng) != 0) { errorMessage = "could not seed rng"; goto createRNG_cleanup; }
            }
            if (numSeedResets == static_cast<size_t>(-1)) for (numSeedResets = 0; numSeedResets < static_cast<size_t>(-1); ++numSeedResets) {
              if (!ext_rng_seedsAreEqual(state[chainNum].rng, state[chainNum - 1].rng)) break;
              
              if (ext_rng_setSeed(state[chainNum].rng, ext_rng_simulateUnsignedIntegerUniformInRange(state[chainNum - 1].rng, 0, static_cast<uint_least32_t>(-1))) != 0) { errorMessage = "could not seed rng"; goto createRNG_cleanup; }
            }
            if (numSeedResets == static_cast<size_t>(-1)) { errorMessage = "could not obtain unique seed"; goto createRNG_cleanup; }
          }
        }
      } else {
        if ((state[chainNum].rng = ext_rng_create(rng_algorithm, NULL)) == NULL) {
          errorMessage = "could not allocate rng";
          goto createRNG_cleanup;
        }
        
        if (rng_standardNormal != EXT_RNG_STANDARD_NORMAL_INVALID &&
            rng_standardNormal != EXT_RNG_STANDARD_NORMAL_USER_NORM &&
            ext_rng_setStandardNormalAlgorithm(state[chainNum].rng, rng_standardNormal, NULL) != 0) {
          errorMessage = "could not set rng standard normal";
          goto createRNG_cleanup;
        }
        
        if (rng_algorithm != EXT_RNG_ALGORITHM_USER_UNIFORM) {
          if (ext_rng_setSeedFromClock(state[chainNum].rng) != 0) { errorMessage = "could not seed rng"; goto createRNG_cleanup; }
          
          numSeedResets = 0;
          if (chainNum > 0) {
            for (numSeedResets = 0; numSeedResets < static_cast<size_t>(-1); ++numSeedResets) {
              if (!ext_rng_seedsAreEqual(state[chainNum].rng, state[chainNum - 1].rng)) break;
              
              if (ext_rng_setSeedFromClock(state[chainNum].rng) != 0) { errorMessage = "could not seed rng"; goto createRNG_cleanup; }
            }
            if (numSeedResets == static_cast<size_t>(-1)) for (numSeedResets = 0; numSeedResets < static_cast<size_t>(-1); ++numSeedResets) {
              if (!ext_rng_seedsAreEqual(state[chainNum].rng, state[chainNum - 1].rng)) break;
              
              if (ext_rng_setSeed(state[chainNum].rng, ext_rng_simulateUnsignedIntegerUniformInRange(state[chainNum - 1].rng, 0, static_cast<uint_least32_t>(-1))) != 0) { errorMessage = "could not seed rng"; goto createRNG_cleanup; }
            }
            if (numSeedResets == static_cast<size_t>(-1)) { errorMessage = "could not obtain unique seed"; goto createRNG_cleanup; }
          }
        }
      }
    }
    
    return;
    
createRNG_cleanup:
    for ( /* */ ; chainNum > 0; --chainNum) {
      ext_rng_destroy(state[chainNum - 1].rng);
      state[chainNum - 1].rng = NULL;
    }
      
    ext_throwError(errorMessage);
  }
}

namespace dbarts {
  
  void BARTFit::setRNGState(const void* const* uniformState, const void* const* normalState)
  {
    for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum) {
      if (uniformState != NULL && uniformState[chainNum] != NULL)
        ext_rng_setState(state[chainNum].rng, uniformState[chainNum]);
      if (normalState  != NULL && normalState[chainNum]  != NULL)
        ext_rng_setStandardNormalAlgorithm(state[chainNum].rng, static_cast<ext_rng_standardNormal_t>(control.rng_standardNormal), normalState[chainNum]);
    }
  }
  
}

namespace {
  
  void destroyRNG(BARTFit& fit) {
    for (size_t chainNum = fit.control.numChains; chainNum > 0; --chainNum) {
      ext_rng_destroy(fit.state[chainNum - 1].rng);
      fit.state[chainNum - 1].rng = NULL;
    }
  }
  
  void setInitialFit(BARTFit& fit) {
    Control& control(fit.control);
    Data& data(fit.data);
    ChainScratch* chainScratch(fit.chainScratch);
    
    for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum) {
      misc_setVectorToConstant(chainScratch[chainNum].totalFits, data.numObservations, 0.0);
      
      if (data.numTestObservations > 0)
        misc_setVectorToConstant(chainScratch[chainNum].totalTestFits, data.numTestObservations, 0.0);
    }
  }
  
  void initializeLatents(BARTFit& fit) {
    for (size_t chainNum = 0; chainNum < fit.control.numChains; ++chainNum)
      initializeLatents(fit, chainNum);
   
#ifndef MATCH_BAYES_TREE
    // shouldn't be used, but will leave at reasonable values; if anyone cares, should
    // look at offset var for min/max/range
    fit.sharedScratch.dataScale.min = -1.0;
    fit.sharedScratch.dataScale.max =  1.0;
    fit.sharedScratch.dataScale.range = 2.0;
#else
    fit.sharedScratch.dataScale.min = -2.0;
    fit.sharedScratch.dataScale.max =  0.0;
    fit.sharedScratch.dataScale.range = 2.0;
#endif
  }
  
  void initializeLatents(BARTFit& fit, size_t chainNum) {
    const Data& data(fit.data);
    
    double* z = fit.chainScratch[chainNum].probitLatents;
    
    // z = 2.0 * y - 1.0 - offset; so -1 if y == 0 and 1 if y == 1 when offset == 0
#ifndef MATCH_BAYES_TREE
    misc_setVectorToConstant(z, data.numObservations, -1.0);
    if (data.offset != NULL) misc_addVectorsInPlace(data.offset, data.numObservations, -1.0, z);
    misc_addVectorsInPlace(data.y, data.numObservations, 2.0, z);
#else
    // BayesTree initialized the latents to be -2 and 0; was probably a bug
    misc_setVectorToConstant(z, data.numObservations, -2.0);
    if (data.offset != NULL) misc_addVectorsInPlace(data.offset, data.numObservations, -1.0, z);
    misc_addVectorsInPlace(data.y, data.numObservations, 2.0, z);
#endif
  }
  
  void rescaleResponse(BARTFit& fit) {
    const Data& data(fit.data);
    SharedScratch& sharedScratch(fit.sharedScratch);
    
    double* yRescaled = const_cast<double*>(fit.sharedScratch.yRescaled);
    
    if (data.offset != NULL) {
      misc_addVectors(data.offset, data.numObservations, -1.0, data.y, yRescaled);
    } else {
      std::memcpy(yRescaled, data.y, data.numObservations * sizeof(double));
    }
    
    sharedScratch.dataScale.min = yRescaled[0];
    sharedScratch.dataScale.max = yRescaled[0];
    for (size_t i = 1; i < data.numObservations; ++i) {
      if (yRescaled[i] < sharedScratch.dataScale.min) sharedScratch.dataScale.min = yRescaled[i];
      if (yRescaled[i] > sharedScratch.dataScale.max) sharedScratch.dataScale.max = yRescaled[i];
    }
    sharedScratch.dataScale.range = sharedScratch.dataScale.max - sharedScratch.dataScale.min;
    if (sharedScratch.dataScale.max == sharedScratch.dataScale.min) sharedScratch.dataScale.range = 1.0;
    
    // yRescaled = (y - offset - min) / (max - min) - 0.5
    misc_addScalarToVectorInPlace(   yRescaled, data.numObservations, -sharedScratch.dataScale.min);
    misc_scalarMultiplyVectorInPlace(yRescaled, data.numObservations, 1.0 / sharedScratch.dataScale.range);
    misc_addScalarToVectorInPlace(   yRescaled, data.numObservations, -0.5);
  }
  
  // multithread-this!
  // 
  void sampleProbitLatentVariables(const BARTFit& fit, State& state, const double* fits, double* z) {
    for (size_t i = 0; i < fit.data.numObservations; ++i) {      
#ifndef MATCH_BAYES_TREE
      double mean = fits[i];
      double offset = 0.0;
      if (fit.data.offset != NULL) offset = fit.data.offset[i];
      
      if (fit.data.y[i] > 0.0) {
        z[i] = ext_rng_simulateLowerTruncatedNormalScale1(state.rng, mean, -offset);
      } else {
        z[i] = ext_rng_simulateUpperTruncatedNormalScale1(state.rng, mean, -offset);
      }
#else
      double prob;
      
      double mean = fits[i];
      if (fit.data.offset != NULL) mean += fit.data.offset[i];
      
      double u = ext_rng_simulateContinuousUniform(state.rng);
      if (fit.data.y[i] > 0.0) {
        prob = u + (1.0 - u) * ext_cumulativeProbabilityOfNormal(0.0, mean, 1.0);
        z[i] = ext_quantileOfNormal(prob, mean, 1.0);
      } else {
        prob = u + (1.0 - u) * ext_cumulativeProbabilityOfNormal(0.0, -mean, 1.0);
        z[i] = mean - ext_quantileOfNormal(prob, 0.0, 1.0);
      }
#endif
      
    }
  }
  
  void storeSamples(const BARTFit& fit, size_t chainNum, Results& results,
                    const double* trainingSample, const double* testSample,
                    double sigma, double k, const uint32_t* variableCounts, size_t simNum)
  {
    const Data& data(fit.data);
    const Control& control(fit.control);
    const Model& model(fit.model);
    const SharedScratch& sharedScratch(fit.sharedScratch);
    
    size_t chainStride = chainNum * results.numSamples;
    if (control.responseIsBinary) {
      if (control.keepTrainingFits) {
        double* trainingSamples = results.trainingSamples + (simNum + chainStride) * data.numObservations;
        std::memcpy(trainingSamples, trainingSample, data.numObservations * sizeof(double));
        
        if (data.offset != NULL) misc_addVectorsInPlace(data.offset, data.numObservations, 1.0, trainingSamples);
      }
      
      if (data.numTestObservations > 0) {
        double* testSamples = results.testSamples + (simNum + chainStride) * data.numTestObservations;
        std::memcpy(testSamples, testSample, data.numTestObservations * sizeof(double));
        if (data.testOffset != NULL) misc_addVectorsInPlace(data.testOffset, data.numTestObservations, 1.0, testSamples);
      }
      
      results.sigmaSamples[simNum + chainStride] = 1.0;
    } else {
      if (control.keepTrainingFits) {
        double* trainingSamples = results.trainingSamples + (simNum + chainStride) * data.numObservations;
        // set training to dataScale.range * (totalFits + 0.5) + dataScale.min + offset
        misc_setVectorToConstant(trainingSamples, data.numObservations, sharedScratch.dataScale.range * 0.5 + sharedScratch.dataScale.min);
        misc_addVectorsInPlace(trainingSample, data.numObservations, sharedScratch.dataScale.range, trainingSamples);
        if (data.offset != NULL) misc_addVectorsInPlace(data.offset, data.numObservations, 1.0, trainingSamples);
      }
      
      if (data.numTestObservations > 0) {
        double* testSamples = results.testSamples + (simNum + chainStride) * data.numTestObservations;
        misc_setVectorToConstant(testSamples, data.numTestObservations, sharedScratch.dataScale.range * 0.5 + sharedScratch.dataScale.min);
        misc_addVectorsInPlace(testSample, data.numTestObservations, sharedScratch.dataScale.range, testSamples);
        if (data.testOffset != NULL) misc_addVectorsInPlace(data.testOffset, data.numTestObservations, 1.0, testSamples);
      }
       
      results.sigmaSamples[simNum + chainStride] = sigma * sharedScratch.dataScale.range;
    }
    
    if (model.kPrior != NULL)
      results.kSamples[simNum + chainStride] = k;
    
    double* variableCountSamples = results.variableCountSamples + (simNum + chainStride) * data.numPredictors;
    for (size_t j = 0; j < data.numPredictors; ++j) variableCountSamples[j] = static_cast<double>(variableCounts[j]);
  }
  
  
  void countVariableUses(const BARTFit& fit, const State& state, uint32_t* variableCounts)
  {
    for (size_t treeNum = 0; treeNum < fit.control.numTrees; ++treeNum)
      state.trees[treeNum ].countVariableUses(variableCounts);
  }

#ifdef HAVE_GETTIMEOFDAY
  double subtractTimes(struct timeval end, struct timeval start) {
    return (1.0e6 * (static_cast<double>(end.tv_sec - start.tv_sec)) + static_cast<double>(end.tv_usec - start.tv_usec)) / 1.0e6;
  }
#else
  double subtractTimes(time_t end, time_t start) { return static_cast<double>(end - start); }
#endif
}

