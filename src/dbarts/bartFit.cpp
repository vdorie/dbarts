#include "config.hpp"
#include <dbarts/bartFit.hpp>

#include <cmath>     // sqrt
#include <cstring>   // memcpy
#include <cstddef>   // size_t

#if !defined(HAVE_SYS_TIME_H) && defined(HAVE_GETTIMEOFDAY)
#  undef HAVE_GETTIMEOFDAY
#endif
#ifdef HAVE_SYS_TIME_H
#  include <sys/time.h> // gettimeofday
#else
#  include <time.h>
#endif


#include <set>       // used to sort and find 
#include <vector>    //   split points
#include <algorithm> // integer min

#include <external/alloca.h>
#include <external/io.h>
#include <external/random.h>
#include <external/stats.h>
#include <external/stats_mt.h>
#include <external/linearAlgebra.h>

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
  void setInitialFit(BARTFit& fit);
  
  void setPrior(BARTFit& fit);
  
  void setCutPoints(BARTFit& fit, const size_t* columns, size_t numColumns);
  void setCutPointsFromQuantiles(BARTFit& fit, const double* x, uint32_t maxNumCuts,
                                 uint32_t& numCutsPerVariable, double*& cutPoints,
                                 std::set<double>& uniqueElements, std::vector<double>& sortedElements);
  void setCutPointsUniformly(BARTFit& fit, const double* x, uint32_t maxNumCuts,
                             uint32_t& numCutsPerVariable, double*& cutPoints);
  
  void setNewObservationIndices(Node& newNode, size_t* indices, const Node& oldNode);
  
  void printInitialSummary(const BARTFit& fit);
  void printTerminalSummary(const BARTFit& fit);
  
  void initializeLatents(BARTFit& fit);
  void initializeLatents(BARTFit& fit, size_t chainNum);
  void rescaleResponse(BARTFit& fit);
  // void resampleTreeFits(BARTFit& fit);
  
  void runSampler(BARTFit& fit, size_t chainNum, size_t numBurnIn, Results* resultsPointer);
  
  void sampleProbitLatentVariables(const BARTFit& fit, State& state, const double* fits, double* yRescaled);
  void storeSamples(const BARTFit& fit, size_t chainNum, Results& results,
                    const double* trainingSample, const double* testSample,
                    double sigma, const uint32_t* variableCounts, size_t simNum);
  void countVariableUses(const BARTFit& fit, const State& state, uint32_t* variableCounts);
  
#ifdef HAVE_SYS_TIME_H
  double subtractTimes(struct timeval end, struct timeval start);
#else
  double subtractTimes(time_t end, time_t start);
#endif
}

namespace dbarts {
  
  void BARTFit::setResponse(const double* newY) {
    if (!control.responseIsBinary) {
      double* sigmaUnscaled = ext_stackAllocate(control.numChains, double);
      for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum)
        sigmaUnscaled[chainNum] = state[chainNum].sigma * sharedScratch.dataScale.range;
      
      double priorUnscaled = model.sigmaSqPrior->getScale() * sharedScratch.dataScale.range * sharedScratch.dataScale.range;
      
      data.y = newY;
      
      rescaleResponse(*this);
      
      model.sigmaSqPrior->setScale(priorUnscaled / (sharedScratch.dataScale.range * sharedScratch.dataScale.range));
      
      for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum)
        state[chainNum].sigma = sigmaUnscaled[chainNum] / sharedScratch.dataScale.range;
      ext_stackFree(sigmaUnscaled);
     
    } else {
      data.y = newY;
      
      for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum)
        sampleProbitLatentVariables(*this, state[chainNum], const_cast<const double*>(state[chainNum].totalFits), chainScratch[chainNum].probitLatents);
    }
    
    // resampleTreeFits(*this);
  }
  
  void BARTFit::setOffset(const double* newOffset) {
    if (!control.responseIsBinary) {
      double* sigmaUnscaled = ext_stackAllocate(control.numChains, double);
      for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum)
        sigmaUnscaled[chainNum] = state[chainNum].sigma * sharedScratch.dataScale.range;
      
      double priorUnscaled = model.sigmaSqPrior->getScale() * sharedScratch.dataScale.range * sharedScratch.dataScale.range;
      
      data.offset = newOffset;
      
      rescaleResponse(*this);
      
      model.sigmaSqPrior->setScale(priorUnscaled / (sharedScratch.dataScale.range * sharedScratch.dataScale.range));
      
      for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum)
        state[chainNum].sigma = sigmaUnscaled[chainNum] / sharedScratch.dataScale.range;
      ext_stackFree(sigmaUnscaled);
    } else {
      data.offset = newOffset;
      
      for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum)
        sampleProbitLatentVariables(*this, state[chainNum], const_cast<const double*>(state[chainNum].totalFits), chainScratch[chainNum].probitLatents);
    }
  }

  // this can leave the tree structures in an invalid state and doesn't roll-back
  bool BARTFit::setPredictor(const double* newPredictor)
  {
    size_t* columns = ext_stackAllocate(data.numPredictors, size_t);
    for (size_t i = 0; i < data.numPredictors; ++i) columns[i] = i;
    
    setCutPoints(*this, columns, data.numPredictors);
    
    ext_stackFree(columns);
    
    data.x = newPredictor;
    
    ext_transposeMatrix(data.x, data.numObservations, data.numPredictors, const_cast<double*>(sharedScratch.xt));
    
    double** nodePosteriorPredictions = new double*[control.numTrees];
    for (size_t treeNum = 0; treeNum < control.numTrees; ++treeNum) nodePosteriorPredictions[treeNum] = NULL;
    
    bool treesAreValid = true;
    for (size_t chainNum = 0; chainNum < control.numChains && treesAreValid == true; ++chainNum) {
      size_t treeNum = 0;
      for ( ; treeNum < control.numTrees && treesAreValid == true; ++treeNum) {
        const double* treeFits = state[chainNum].treeFits + treeNum * data.numObservations;
        
        // next allocates memory
        nodePosteriorPredictions[treeNum] = state[chainNum].trees[treeNum].recoverAveragesFromFits(*this, treeFits);
        
        state[chainNum].trees[treeNum].top.addObservationsToChildren(*this);
        
        treesAreValid &= state[chainNum].trees[treeNum].isValid();
      }
      
      if (treesAreValid) {
        // go back across bottoms and set predictions to those mus for obs now in node
        for (size_t treeNum = 0; treeNum < control.numTrees; ++treeNum) {
          double* treeFits = state[chainNum].treeFits + treeNum * data.numObservations;
          
          ext_addVectorsInPlace(treeFits, data.numObservations, -1.0, state[chainNum].totalFits);
          
          state[chainNum].trees[treeNum].setCurrentFitsFromAverages(*this, nodePosteriorPredictions[treeNum], treeFits, NULL);
          for (int32_t j = 0; j < static_cast<int32_t>(data.numPredictors); ++j)
            updateVariablesAvailable(*this, state[chainNum].trees[treeNum].top, j);
          
          ext_addVectorsInPlace(treeFits, data.numObservations, 1.0, state[chainNum].totalFits);
        }
      }
      
      for (size_t treeNum = control.numTrees; treeNum > 0; --treeNum) {
        delete [] nodePosteriorPredictions[treeNum - 1];
        nodePosteriorPredictions[treeNum - 1] = NULL;
      }
    }
    
    delete [] nodePosteriorPredictions;
    
    return treesAreValid;
  }
  
  bool BARTFit::updatePredictor(const double* newPredictor, size_t column)
  {
    return updatePredictors(newPredictor, &column, 1);
  }
  
  bool BARTFit::updatePredictors(const double* newPredictor, const size_t* columns, size_t numColumns)
  {
    // store current
    double* oldPredictor = new double[data.numObservations * numColumns];
    double** oldCutPoints = new double*[numColumns];
    
    for (size_t i = 0; i < numColumns; ++i) {
      std::memcpy(oldPredictor + i * data.numObservations, data.x + columns[i] * data.numObservations, data.numObservations * sizeof(double));
      oldCutPoints[i] = new double[sharedScratch.numCutsPerVariable[columns[i]]];
      std::memcpy(oldCutPoints[i], sharedScratch.cutPoints[columns[i]], sharedScratch.numCutsPerVariable[columns[i]] * sizeof(double));
    }
    
    
    // install new
    setCutPoints(*this, columns, numColumns);
    
    double* x  = const_cast<double*>(data.x);
    double* xt = const_cast<double*>(sharedScratch.xt);
    for (size_t j = 0; j < numColumns; ++j) {
      std::memcpy(x + columns[j] * data.numObservations, newPredictor + j * data.numObservations, data.numObservations * sizeof(double));
      for (size_t i = 0; i < data.numObservations; ++i) {
        xt[i * data.numPredictors + columns[j]] = newPredictor[i + j * data.numObservations];
      }
    }
    
    
    // check validity of new columns and recover node posterior samples
    size_t totalNumTrees = control.numChains * control.numTrees;
    double** nodePosteriorPredictions = new double*[totalNumTrees];
    for (size_t treeNum = 0; treeNum < totalNumTrees; ++treeNum) nodePosteriorPredictions[treeNum] = NULL;
    
    bool treesAreValid = true;
    for (size_t chainNum = 0; chainNum < control.numChains && treesAreValid == true; ++chainNum) {
      double** chainNodePosteriorPredictions = nodePosteriorPredictions + chainNum * control.numTrees;
      
      size_t treeNum;
      for (treeNum = 0; treeNum < control.numTrees && treesAreValid == true; ++treeNum) {
        const double* treeFits = state[chainNum].treeFits + treeNum * data.numObservations;
        
        chainNodePosteriorPredictions[treeNum] = state[chainNum].trees[treeNum].recoverAveragesFromFits(*this, treeFits);
        
        state[chainNum].trees[treeNum].top.addObservationsToChildren(*this);
        
        treesAreValid &= state[chainNum].trees[treeNum].isValid();
      }
    }
    
    if (!treesAreValid) {
      // copy back in old
      for (size_t j = 0; j < numColumns; ++j) {
        std::memcpy(x + columns[j] * data.numObservations, oldPredictor + j * data.numObservations, data.numObservations * sizeof(double));
        
        std::memcpy(const_cast<double**>(sharedScratch.cutPoints)[columns[j]], oldCutPoints[j], sharedScratch.numCutsPerVariable[columns[j]] * sizeof(double));
          
        for (size_t i = 0; i < data.numObservations; ++i) {
          xt[i * data.numPredictors + columns[j]] = oldPredictor[i + j * data.numObservations];
        }
      }
      
      for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum)  
        for (size_t treeNum = 0; treeNum < control.numTrees; ++treeNum) state[chainNum].trees[treeNum].top.addObservationsToChildren(*this);
    } else {
      // go back across bottoms and set predictions to those mus for obs now in node
      for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum) {
        double** chainNodePosteriorPredictions = nodePosteriorPredictions + chainNum * control.numTrees;
        
        for (size_t treeNum = 0; treeNum < control.numTrees; ++treeNum) {
          double* treeFits = state[chainNum].treeFits + treeNum * data.numObservations;
          
          ext_addVectorsInPlace(treeFits, data.numObservations, -1.0, state[chainNum].totalFits);
          
          state[chainNum].trees[treeNum].setCurrentFitsFromAverages(*this, chainNodePosteriorPredictions[treeNum], treeFits, NULL);
          for (int32_t j = 0; j < static_cast<int32_t>(data.numPredictors); ++j)
            updateVariablesAvailable(*this, state[chainNum].trees[treeNum].top, j);
          
          ext_addVectorsInPlace(treeFits, data.numObservations, 1.0, state[chainNum].totalFits);
        }
      }
    }
    
    for (size_t treeNum = totalNumTrees; treeNum > 0; --treeNum) delete [] nodePosteriorPredictions[treeNum - 1];
    delete [] nodePosteriorPredictions;
    
    for (size_t j = 0; j < numColumns; ++j) delete [] oldCutPoints[j];
    delete [] oldCutPoints;
    delete [] oldPredictor;
    
    return treesAreValid;
  }
  
#define INVALID_ADDRESS reinterpret_cast<const double*>(this)
  void BARTFit::setTestPredictor(const double* newTestPredictor, size_t numTestObservations) {
    setTestPredictorAndOffset(newTestPredictor, INVALID_ADDRESS, numTestObservations);
  }
  
  void BARTFit::setTestOffset(const double* newTestOffset) {
     data.testOffset = newTestOffset;
  }
  // setting testOffset to NULL is valid
  // an invalid pointer address for testOffset is the object itself; when invalid, it is not updated
  void BARTFit::setTestPredictorAndOffset(const double* x_test, const double* testOffset, size_t numTestObservations) {
    if (numTestObservations == 0 || x_test == NULL) {
      if (sharedScratch.xt_test != NULL) { delete [] sharedScratch.xt_test; sharedScratch.xt_test = NULL; }
      for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum)
        if (state[chainNum].totalTestFits != NULL) { delete [] state[chainNum].totalTestFits; state[chainNum].totalTestFits = NULL; }
      
      data.x_test = NULL;
      data.numTestObservations = 0;
      data.testOffset = NULL;
    } else {
      data.x_test = x_test;
      
      if (numTestObservations != data.numTestObservations) {
        if (sharedScratch.xt_test != NULL) { delete [] sharedScratch.xt_test; sharedScratch.xt_test = NULL; }
        for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum)
          if (state[chainNum].totalTestFits != NULL) { delete [] state[chainNum].totalTestFits; state[chainNum].totalTestFits = NULL; }
        data.numTestObservations = numTestObservations;
        
        sharedScratch.xt_test = new double[data.numTestObservations * data.numPredictors];
        for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum)
          state[chainNum].totalTestFits = new double[data.numTestObservations];
      }
      
      ext_transposeMatrix(data.x_test, data.numTestObservations, data.numPredictors, const_cast<double*>(sharedScratch.xt_test));
      
      if (testOffset != INVALID_ADDRESS) data.testOffset = testOffset;
      
      double* currTestFits = new double[data.numTestObservations];
      
      for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum) {
        ext_setVectorToConstant(state[chainNum].totalTestFits, data.numTestObservations, 0.0);
        
        for (size_t treeNum = 0; treeNum < control.numTrees; ++treeNum) {
          const double* treeFits = state[chainNum].treeFits + treeNum * data.numObservations;
        
          const double* nodePosteriorPredictions = state[chainNum].trees[treeNum].recoverAveragesFromFits(*this, treeFits);
      
          state[chainNum].trees[treeNum].setCurrentFitsFromAverages(*this, nodePosteriorPredictions, NULL, currTestFits);
          
          ext_addVectorsInPlace(currTestFits, data.numTestObservations, 1.0, state[chainNum].totalTestFits);
        
          delete [] nodePosteriorPredictions;
        }
      }
      
      delete [] currTestFits;
    }
  }
#undef INVALID_ADDRESS
  
  void BARTFit::updateTestPredictor(const double* newTestPredictor, size_t column) {
    updateTestPredictors(newTestPredictor, &column, 1);
  }
  
  void BARTFit::updateTestPredictors(const double* newTestPredictor, const size_t* columns, size_t numColumns) {
    double* x_test = const_cast<double*>(data.x_test);
    double* xt_test = const_cast<double*>(sharedScratch.xt_test);
    
    for (size_t j_ind = 0; j_ind < numColumns; ++j_ind) {
      size_t j = columns[j_ind];
      std::memcpy(x_test + j * data.numTestObservations, newTestPredictor + j_ind * data.numTestObservations, data.numTestObservations * sizeof(double));
      
      for (size_t i = 0; i < data.numTestObservations; ++i) {
        xt_test[i * data.numPredictors + j] = newTestPredictor[i + j_ind * data.numTestObservations];
      }
    }
    
    double* currTestFits = new double[data.numTestObservations];
    
    for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum) {
      ext_setVectorToConstant(state[chainNum].totalTestFits, data.numTestObservations, 0.0);
      
      for (size_t treeNum = 0; treeNum < control.numTrees; ++treeNum) {
        const double* treeFits = state[chainNum].treeFits + treeNum * data.numObservations;
        
        const double* nodePosteriorPredictions = state[chainNum].trees[treeNum].recoverAveragesFromFits(*this, treeFits);
      
        state[chainNum].trees[treeNum].setCurrentFitsFromAverages(*this, nodePosteriorPredictions, NULL, currTestFits);
      
        ext_addVectorsInPlace(currTestFits, data.numTestObservations, 1.0, state[chainNum].totalTestFits);
      
        delete [] nodePosteriorPredictions;
      }
    }
    
    delete [] currTestFits;
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
   *   sharedScratch.cutPoints and sharedScratch.numCutsPerVariable - compute from data.maxNumCuts and new data
   *
   *   state.trees         - we have to go through these and prune any now-invalid end nodes
   *   state.treeIndices   - resize and assign into trees
   *   state.treeFits      - resize and recompute for new xt
   *   state.totalFits     - same
   *   state.totalTestFits - same
   *   state.sigma         - rescale using new scratch.dataScale
   */
  void BARTFit::setData(const Data& newData)
  {
    size_t oldNumObservations     = data.numObservations;
    size_t oldNumTestObservations = data.numTestObservations;
    
    data = newData;
    
    if (oldNumObservations != data.numObservations) {
      // handle resizing arrays
      delete [] sharedScratch.xt;
      
      sharedScratch.xt = new double[data.numObservations * data.numPredictors];
      
      if (!control.responseIsBinary) {
        delete [] sharedScratch.yRescaled;
        sharedScratch.yRescaled = new double[data.numObservations];
      }
    }
    
    size_t** oldTreeIndices = ext_stackAllocate(control.numChains, size_t*);
    double** oldTreeFits    = ext_stackAllocate(control.numChains, double*);
    double** currTestFits   = ext_stackAllocate(control.numChains, double*);
    
    for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum) {
      // extract from old data what we'll need to update
      oldTreeIndices[chainNum] = state[chainNum].treeIndices;
      oldTreeFits[chainNum]    = state[chainNum].treeFits;
      currTestFits[chainNum]   = NULL;
      
      if (oldNumObservations != data.numObservations) {
        delete [] state[chainNum].totalFits;
        delete [] chainScratch[chainNum].treeY;
      
        chainScratch[chainNum].treeY = new double[data.numObservations];
        state[chainNum].totalFits    = new double[data.numObservations];
        
        if (control.responseIsBinary) {
          delete [] chainScratch[chainNum].probitLatents;
          chainScratch[chainNum].probitLatents = new double[data.numObservations];
        }
        
        state[chainNum].treeIndices = new size_t[data.numObservations * control.numTrees];
        state[chainNum].treeFits    = new double[data.numObservations * control.numTrees];
      }
    }
    
    // update sharedScratch.yRescaled/chainScratch.probitLatents and state.sigma
    if (control.responseIsBinary)
      initializeLatents(*this);
    else {
      double* sigmaUnscaled = ext_stackAllocate(control.numChains, double);
      for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum)
        sigmaUnscaled[chainNum] = state[chainNum].sigma * sharedScratch.dataScale.range;
      
      double priorUnscaled = model.sigmaSqPrior->getScale() * sharedScratch.dataScale.range * sharedScratch.dataScale.range;
      // ext_printf("set data:\n  sigma est: %f\n", data.sigmaEstimate);
      // ext_printf("  before sigma: %f, scale: %f, range: %f\n", state[0].sigma, model.sigmaSqPrior->getScale(), sharedScratch.dataScale.range);
      
      rescaleResponse(*this);
      
      model.sigmaSqPrior->setScale(priorUnscaled / (sharedScratch.dataScale.range * sharedScratch.dataScale.range));
      
      for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum)
        state[chainNum].sigma = sigmaUnscaled[chainNum] / sharedScratch.dataScale.range;
      
      ext_stackFree(sigmaUnscaled);
      
      // ext_printf("  after sigma: %f, scale: %f, range: %f\n", state[0].sigma, model.sigmaSqPrior->getScale(), sharedScratch.dataScale.range);
    }
        
    // cache old cut points, for use in updating trees
    const double** oldCutPoints = ext_stackAllocate(data.numPredictors, const double*);
    for (size_t j = 0; j < data.numPredictors; ++j) {
      oldCutPoints[j] = sharedScratch.cutPoints[j];
      // next assignments 'reset' the variables, so setCutPoints() ignores old values
      const_cast<uint32_t*>(sharedScratch.numCutsPerVariable)[j] = static_cast<uint32_t>(-1);
      const_cast<double**>(sharedScratch.cutPoints)[j] = NULL;
    }
    // set new cut points
    size_t* columns = ext_stackAllocate(data.numPredictors, size_t);
    for (size_t j = 0; j < data.numPredictors; ++j) columns[j] = j;
    setCutPoints(*this, columns, data.numPredictors);
    ext_stackFree(columns);
    
    // now initialize remaining arrays that use numObs
    ext_transposeMatrix(data.x, data.numObservations, data.numPredictors, const_cast<double*>(sharedScratch.xt)); 
    for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum)
      ext_setVectorToConstant(state[chainNum].totalFits, data.numObservations, 0.0);
    
    
    if (data.numTestObservations == 0 || data.x_test == NULL) {
      delete [] sharedScratch.xt_test;
      sharedScratch.xt_test = NULL;
      for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum) {
        delete [] state[chainNum].totalTestFits;
        state[chainNum].totalTestFits = NULL;
      }
    } else {
      // handle resizing test arrays and initializing them
      if (oldNumTestObservations != data.numTestObservations) {
        delete [] sharedScratch.xt_test;
        sharedScratch.xt_test = new double[data.numTestObservations * data.numPredictors];
        
        for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum) {
          delete [] state[chainNum].totalTestFits;
          state[chainNum].totalTestFits = new double[data.numTestObservations];
        }
      }
      
      ext_transposeMatrix(data.x_test, data.numTestObservations, data.numPredictors, const_cast<double*>(sharedScratch.xt_test));

      for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum) {
        currTestFits[chainNum] = new double[data.numTestObservations];
        ext_setVectorToConstant(state[chainNum].totalTestFits, data.numTestObservations, 0.0);
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
        double* nodePosteriorPredictions = state[chainNum].trees[treeNum].recoverAveragesFromFits(*this, oldTreeFits_i);
      
        // the mapping can end up with some end-nodes that no longer exist, handles that internally
        state[chainNum].trees[treeNum].mapOldCutPointsOntoNew(*this, oldCutPoints, nodePosteriorPredictions);
        
        if (oldNumObservations != data.numObservations) {
          state[chainNum].trees[treeNum].top.observationIndices = state[chainNum].treeIndices + treeNum * data.numObservations;
          state[chainNum].trees[treeNum].top.numObservations = data.numObservations;
        }
        
        state[chainNum].trees[treeNum].top.addObservationsToChildren(*this);
        state[chainNum].trees[treeNum].collapseEmptyNodes(*this, nodePosteriorPredictions);
        for (int32_t i = 0; i < static_cast<int32_t>(data.numPredictors); ++i)
          updateVariablesAvailable(*this, state[chainNum].trees[treeNum].top, i);
        
        double* currTreeFits = state[chainNum].treeFits + treeNum * data.numObservations;
        state[chainNum].trees[treeNum].setCurrentFitsFromAverages(*this, nodePosteriorPredictions, currTreeFits, currTestFits[chainNum]);
        ext_addVectorsInPlace(currTreeFits, data.numObservations, 1.0, state[chainNum].totalFits);
        
        if (data.numTestObservations > 0)
          ext_addVectorsInPlace(currTestFits[chainNum], data.numTestObservations, 1.0, state[chainNum].totalTestFits);
        
        delete [] nodePosteriorPredictions;
      }
    }
    
    for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum)
      delete [] currTestFits[chainNum]; // can be NULL, no big deal
    
    for (size_t i = 0; i < data.numPredictors; ++i) delete [] oldCutPoints[i];
    ext_stackFree(oldCutPoints);
    
    if (oldNumObservations != data.numObservations) {
      for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum) {
        delete [] oldTreeFits[chainNum];
        delete [] oldTreeIndices[chainNum];
      }
    }
    
    ext_stackFree(currTestFits);
    ext_stackFree(oldTreeFits);
    ext_stackFree(oldTreeIndices);
    
    // ext_printf("after setting data\n");
    // printInitialSummary(*this);
  }
  
  void BARTFit::setControl(const Control& newControl)
  {
    if (newControl.numTrees != control.numTrees) {
      for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum) {
        Tree* oldTrees = state[chainNum].trees;
        size_t* oldTreeIndices = state[chainNum].treeIndices;
        double* oldTreeFits = state[chainNum].treeFits;
        
        state[chainNum].trees       = static_cast<Tree*>(::operator new (newControl.numTrees * sizeof(Tree)));
        state[chainNum].treeIndices = new size_t[data.numObservations * newControl.numTrees];
        state[chainNum].treeFits    = new double[data.numObservations * newControl.numTrees];
        
        size_t assignEnd = std::min(control.numTrees, newControl.numTrees);
        for (size_t treeNum = 0; treeNum < assignEnd; ++treeNum) {
          state[chainNum].trees[treeNum] = oldTrees[treeNum];
          setNewObservationIndices(state[chainNum].trees[treeNum].top,
                                   state[chainNum].treeIndices + treeNum * data.numObservations,
                                   oldTrees[treeNum].top);
  
          if (!state[chainNum].trees[treeNum].top.isBottom()) {
            state[chainNum].trees[treeNum].top.getRightChild()->parent = &state[chainNum].trees[treeNum].top;
            state[chainNum].trees[treeNum].top.getLeftChild()->parent  = &state[chainNum].trees[treeNum].top;
            // oldTrees[i].~Tree(); // destructor not needed, as assignment means no memory is lost
          }
        }
        std::memcpy(state[chainNum].treeIndices, oldTreeIndices, assignEnd * data.numObservations * sizeof(size_t));
        std::memcpy(state[chainNum].treeFits,    oldTreeFits,    assignEnd * data.numObservations * sizeof(double));
        
        for (size_t treeNum = assignEnd; treeNum < newControl.numTrees; ++treeNum) {
          new (state[chainNum].trees + treeNum) Tree(state[chainNum].treeIndices + treeNum * data.numObservations, data.numObservations, data.numPredictors);
          for (size_t i = 0; i < data.numObservations; ++i) {
            state[chainNum].treeFits[treeNum * data.numObservations + i] = 0.0;
          }
        }
        
        for (size_t treeNum = control.numTrees; treeNum > assignEnd; /* */) oldTrees[--treeNum].~Tree();
        
        ::operator delete (oldTrees);
        delete [] oldTreeIndices;
        delete [] oldTreeFits;
        
        NormalPrior* nodePrior = static_cast<NormalPrior*>(model.muPrior);
        double precisionUnscaled = nodePrior->precision / static_cast<double>(control.numTrees);
        nodePrior->precision = precisionUnscaled * static_cast<double>(newControl.numTrees);
      }
    }
    
    ext_rng_algorithm_t old_rng_algorithm = control.rng_algorithm;
    ext_rng_standardNormal_t old_rng_standardNormal = control.rng_standardNormal;
    
    control = newControl;
    
    if (old_rng_algorithm != control.rng_algorithm || old_rng_standardNormal != control.rng_standardNormal) {
      destroyRNG(*this);
      createRNG(*this);
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
                           const size_t* treeIndices, size_t numTreeIndices) const {
    for (size_t i = 0; i < numChainIndices; ++i) {
      size_t chainNum = chainIndices[i];
      for (size_t j = 0; j < numTreeIndices; ++j) {
        size_t treeNum = treeIndices[j];
        
        const double* treeFits = state[chainNum].treeFits + treeNum * data.numObservations;
        double* nodePosteriorPredictions = state[chainNum].trees[treeNum].recoverAveragesFromFits(*this, treeFits);
        
        NodeVector bottomNodes(const_cast<Tree*>(state[chainNum].trees + treeNum)->top.getBottomVector());
        size_t numBottomNodes = bottomNodes.size();
        for (size_t k = 0; k < numBottomNodes; ++k) bottomNodes[k]->setAverage(nodePosteriorPredictions[k]);
        delete [] nodePosteriorPredictions;
        
        state[chainNum].trees[treeNum].top.print(*this);
      }
    }
  }
  
  BARTFit::BARTFit(Control control, Model model, Data data) :
    control(control), model(model), data(data), sharedScratch(), chainScratch(NULL), state(NULL),
    runningTime(0.0), threadManager(NULL)
  {
    allocateMemory(*this);

    createRNG(*this);
    
    setPrior(*this);
    setInitialCutPoints(*this);
    setInitialFit(*this);

    if (this->control.verbose) printInitialSummary(*this);
  }
  
  BARTFit::~BARTFit()
  {
    destroyRNG(*this);
    
    delete [] sharedScratch.yRescaled; sharedScratch.yRescaled = NULL;
    delete [] sharedScratch.xt; sharedScratch.xt = NULL;
    delete [] sharedScratch.xt_test; sharedScratch.xt_test = NULL;
    for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum) {
      delete [] chainScratch[chainNum].treeY; chainScratch[chainNum].treeY = NULL;
      delete [] chainScratch[chainNum].probitLatents; chainScratch[chainNum].probitLatents = NULL;
    }
    
    delete [] chainScratch;
    
    delete [] sharedScratch.numCutsPerVariable; sharedScratch.numCutsPerVariable = NULL;
    if (sharedScratch.cutPoints != NULL) {
      for (size_t i = 0; i < data.numPredictors; ++i) delete [] sharedScratch.cutPoints[i];
    }
    delete [] sharedScratch.cutPoints; sharedScratch.cutPoints = NULL;
    
    for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum) {
      if (state[chainNum].trees != NULL) for (size_t treeNum = control.numTrees; treeNum > 0; ) state[chainNum].trees[--treeNum].~Tree();
      ::operator delete (state[chainNum].trees); state[chainNum].trees = NULL;
      delete [] state[chainNum].treeIndices; state[chainNum].treeIndices = NULL;
      
      delete [] state[chainNum].treeFits; state[chainNum].treeFits = NULL;
      delete [] state[chainNum].totalFits; state[chainNum].totalFits = NULL;
      if (data.numTestObservations > 0) {
        delete [] state[chainNum].totalTestFits;
        state[chainNum].totalTestFits = NULL;
      }
    }
    
    delete [] state;
    
    ext_htm_destroy(threadManager);
  }
  
  Results* BARTFit::runSampler()
  {
    // ensure at least one sample for state's sake
    Results* resultsPointer = new Results(data.numObservations, data.numPredictors,
                                          data.numTestObservations, control.numSamples == 0 ? 1 : control.numSamples,
                                          control.numChains);
    
    runSampler(control.numBurnIn, resultsPointer);
    
    if (control.numSamples == 0) {
      delete resultsPointer;
      return NULL;
    }
    
    return resultsPointer;
  }
  
  Results* BARTFit::runSampler(size_t numBurnIn, size_t numSamples)
  {
    Results* resultsPointer = new Results(data.numObservations, data.numPredictors,
                                          data.numTestObservations, numSamples == 0 ? 1 : numSamples,
                                          control.numChains);
    
    runSampler(numBurnIn, resultsPointer);
    
    if (numSamples == 0) {
      delete resultsPointer;
      return NULL;
    }
    
    return resultsPointer;
  }
  
    
  void BARTFit::sampleTreesFromPrior() {
    
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
    
    bool stepTaken, isThinningIteration;
    StepType ignored;
    
    size_t numSamples = results.numSamples;
    
    double* currFits = new double[data.numObservations];
    double* currTestFits = data.numTestObservations > 0 ? new double[data.numTestObservations] : NULL;
    
    uint32_t* variableCounts = ext_stackAllocate(data.numPredictors, uint32_t);
    
    size_t totalNumIterations = (numBurnIn + numSamples) * control.treeThinningRate;
    size_t majorIterationNum = 0;
    
    // reserve once at the start if possible
    if (control.numThreads > 1 && control.numChains == 1)
      ext_htm_reserveThreadsForSubTask(fit.threadManager, 0, 0);
    
    // const cast b/c yRescaled doesn't change, but probit latents do
    double* y = control.responseIsBinary ? chainScratch.probitLatents : const_cast<double*>(sharedScratch.yRescaled);
    
    for (size_t k = 0; k < totalNumIterations; ++k) {
      if (control.numThreads > 1 && control.numChains > 1)
        ext_htm_reserveThreadsForSubTask(fit.threadManager, taskId, k);
      
      isThinningIteration = ((k + 1) % control.treeThinningRate != 0);
            
      majorIterationNum = k / control.treeThinningRate;
      
      if (control.verbose && !isThinningIteration && (majorIterationNum + 1) % control.printEvery == 0) {
        if (control.numChains > 1)
          ext_htm_printf(fit.threadManager, "[%lu] iteration: %u (of %u)\n", chainNum + 1, k + 1, totalNumIterations);
        else
          ext_printf("iteration: %u (of %u)\n", k + 1, totalNumIterations);
     }
      
      if (!isThinningIteration && data.numTestObservations > 0) ext_setVectorToConstant(state.totalTestFits, data.numTestObservations, 0.0);
      
      for (size_t treeNum = 0; treeNum < control.numTrees; ++treeNum) {
        double* oldTreeFits = state.treeFits + treeNum * data.numObservations;
        
        // treeY = y - (totalFits - oldTreeFits)
        // is residual from every *other* tree, so what is left for this tree to do
        std::memcpy(chainScratch.treeY, y, data.numObservations * sizeof(double));
        ext_addVectorsInPlace(const_cast<const double*>(state.totalFits), data.numObservations, -1.0, chainScratch.treeY);
        ext_addVectorsInPlace(const_cast<const double*>(oldTreeFits), data.numObservations, 1.0, chainScratch.treeY);
        
        state.trees[treeNum].setNodeAverages(fit, chainNum, chainScratch.treeY);
        
        metropolisJumpForTree(fit, chainNum, state.trees[treeNum], chainScratch.treeY, &stepTaken, &ignored);
                
        state.trees[treeNum].sampleAveragesAndSetFits(fit, chainNum, currFits, isThinningIteration ? NULL : currTestFits);
        
        // totalFits += currFits - oldTreeFits
        ext_addVectorsInPlace(const_cast<const double*>(oldTreeFits), data.numObservations, -1.0, state.totalFits);
        ext_addVectorsInPlace(const_cast<const double*>(currFits), data.numObservations, 1.0, state.totalFits);
        
        if (!isThinningIteration && data.numTestObservations > 0) {
          ext_addVectorsInPlace(const_cast<const double*>(currTestFits), data.numTestObservations, 1.0, state.totalTestFits);
        }
        
        std::memcpy(oldTreeFits, const_cast<const double*>(currFits), data.numObservations * sizeof(double));
      }
      
      if (control.responseIsBinary) {
        sampleProbitLatentVariables(fit, state, state.totalFits, y);
      } else {
        double sumOfSquaredResiduals;
        if (data.weights != NULL) {
          sumOfSquaredResiduals = ext_htm_computeWeightedSumOfSquaredResiduals(fit.threadManager, taskId, y, data.numObservations, data.weights, state.totalFits);
        } else {
          sumOfSquaredResiduals = ext_htm_computeSumOfSquaredResiduals(fit.threadManager, taskId, y, data.numObservations, state.totalFits);
        }
        state.sigma = std::sqrt(model.sigmaSqPrior->drawFromPosterior(state.rng, static_cast<double>(data.numObservations), sumOfSquaredResiduals));
      }
      
      if (!isThinningIteration) {
        // if not out of burn-in, store result in first result; start
        // overwriting after that
        bool isBurningIn = majorIterationNum < numBurnIn;
        size_t simNum = (!isBurningIn ? majorIterationNum - numBurnIn : 0);
        
        for (size_t j = 0; j < fit.data.numPredictors; ++j) variableCounts[j] = 0;
        countVariableUses(fit, state, variableCounts);
        
        storeSamples(fit, chainNum, results, state.totalFits, state.totalTestFits, state.sigma, variableCounts, simNum);
        
        if (control.callback != NULL) {
          size_t chainStride = chainNum * results.numSamples;
          control.callback(control.callbackData, fit, isBurningIn,
                           results.trainingSamples + (simNum + chainStride) * data.numObservations,
                           results.testSamples + (simNum + chainStride) * data.numTestObservations,
                           results.sigmaSamples[simNum + chainStride]);
        }
      }
    }
    
    delete [] currFits;
    if (data.numTestObservations > 0) delete [] currTestFits;
    ext_stackFree(variableCounts);
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
      
      ext_htm_runTopLevelTasks(threadManager, &samplerThreadFunction, threadDataPtr, control.numChains);
      
      delete [] threadDataPtr;
      delete [] threadData;
    }
        
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
    
    ext_printf("number of trees: %u\n", control.numTrees);
    
    ext_printf("Prior:\n");
    // dirty hack... should have priors print themselves
    double sigma = std::sqrt(1.0 / static_cast<NormalPrior*>(model.muPrior)->precision);
    double k = (control.responseIsBinary ? 3.0 : 0.5) /  (sigma * std::sqrt(static_cast<double>(control.numTrees)));
    ext_printf("\tk: %f\n", k);
    if (!control.responseIsBinary) {
      ChiSquaredPrior* residPrior = static_cast<ChiSquaredPrior*>(model.sigmaSqPrior);
      ext_printf("\tdegrees of freedom in sigma prior: %f\n", residPrior->degreesOfFreedom);
      // double sigma = data.sigmaEstimate / sharedScratch.dataScale.range;
      double quantile = 1.0 - ext_percentileOfChiSquared(residPrior->scale * residPrior->degreesOfFreedom / fit.state[0].sigma / fit.state[0].sigma, residPrior->degreesOfFreedom);
      ext_printf("\tquantile in sigma prior: %f\n", quantile);
      ext_printf("\tscale in sigma prior: %f\n", residPrior->scale);
    }
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
      ext_printf("(%u: %u) ", j + 1, sharedScratch.numCutsPerVariable[j]);
      if ((j + 1) % 5 == 0) ext_printf("\n");
    }
    ext_printf("\n");
    if (control.printCutoffs > 0) {
      ext_printf("cutoffs:\n");
      for (size_t j = 0; j < data.numPredictors; ++j) {
        ext_printf("x(%u) cutoffs: ", j + 1);
        
        size_t k;
        for (k = 0; k < sharedScratch.numCutsPerVariable[j] - 1 && k < control.printCutoffs - 1; ++k) {
          ext_printf("%f", sharedScratch.cutPoints[j][k]);
          if ((k + 1) % 5 == 0) ext_printf("\n\t");
        }
        if (k > 2 && k == control.printCutoffs && k < sharedScratch.numCutsPerVariable[j] - 1)
          ext_printf("...");
        
        ext_printf("%f", sharedScratch.cutPoints[j][sharedScratch.numCutsPerVariable[j] - 1]);
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
    
    uint32_t* variableCounts = ext_stackAllocate(fit.data.numPredictors, uint32_t);
    
    ext_printf("Variable Usage, last iteration (var:count):\n");
    for (size_t j = 0; j < fit.data.numPredictors; ++j) variableCounts[j] = 0;
    for (size_t chainNum = 0; chainNum < fit.control.numChains; ++chainNum)
      countVariableUses(fit, fit.state[chainNum], variableCounts);
    for (size_t j = 0; j < fit.data.numPredictors; ++j) {
      ext_printf("(%lu: %u) ", static_cast<unsigned long int>(j + 1), variableCounts[j]);
      if ((j + 1) % 5 == 0) ext_printf("\n");
    }
    
    ext_stackFree(variableCounts);
    
    ext_printf("\nDONE BART\n\n");
  }
  
  void allocateMemory(BARTFit& fit) {
    Control& control(fit.control);
    Data& data(fit.data);
    SharedScratch& sharedScratch(fit.sharedScratch);
    
    fit.chainScratch = new ChainScratch[control.numChains];
    ChainScratch* chainScratch = fit.chainScratch;
    
    fit.state = new State[control.numChains];
    State* state = fit.state;
    
    if (!control.responseIsBinary) {
      sharedScratch.yRescaled = new double[data.numObservations];
      for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum)
        chainScratch[chainNum].probitLatents = NULL;
    } else {
      sharedScratch.yRescaled = NULL;
      for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum)
        chainScratch[chainNum].probitLatents = new double[data.numObservations];
    }
    
    if (control.responseIsBinary) initializeLatents(fit);
    else rescaleResponse(fit);
    
    sharedScratch.xt = new double[data.numObservations * data.numPredictors];
    ext_transposeMatrix(data.x, data.numObservations, data.numPredictors, const_cast<double*>(sharedScratch.xt));
    
    if (data.numTestObservations > 0) {
      sharedScratch.xt_test = new double[data.numTestObservations * data.numPredictors];
      ext_transposeMatrix(data.x_test, data.numTestObservations, data.numPredictors, const_cast<double*>(sharedScratch.xt_test));
    }

    // chain scratches
    for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum) {
      chainScratch[chainNum].treeY = new double[data.numObservations];
      double* y = control.responseIsBinary ? chainScratch[chainNum].probitLatents : const_cast<double*>(sharedScratch.yRescaled);
      
      for (size_t i = 0; i < data.numObservations; ++i) chainScratch[chainNum].treeY[i] = y[i];
      
      chainScratch[chainNum].taskId = static_cast<size_t>(-1);
    }
    
    // shared scratch
    sharedScratch.numCutsPerVariable = new uint32_t[data.numPredictors];

    sharedScratch.cutPoints = new double*[data.numPredictors];
    const double** cutPoints = const_cast<const double**>(sharedScratch.cutPoints);
    for (size_t j = 0; j < data.numPredictors; ++j) cutPoints[j] = NULL;
    
    // states
    for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum) {
      state[chainNum].trees = static_cast<Tree*>(::operator new (control.numTrees * sizeof(Tree)));
      state[chainNum].treeIndices = new size_t[data.numObservations * control.numTrees];
      
      for (size_t treeNum = 0; treeNum < control.numTrees; ++treeNum) {
        new (state[chainNum].trees + treeNum) Tree(state[chainNum].treeIndices + treeNum * data.numObservations, data.numObservations, data.numPredictors);
      }
    }
    
    if (control.numThreads > 1 && ext_htm_create(&fit.threadManager, control.numThreads) != 0) {
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
    
    for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum)
      state[chainNum].sigma = control.responseIsBinary ? 1.0 : (data.sigmaEstimate / sharedScratch.dataScale.range);
    model.sigmaSqPrior->setScale(state[0].sigma * state[0].sigma * model.sigmaSqPrior->getScale());
  }
  
  void setInitialCutPoints(BARTFit& fit) {
    Data& data(fit.data);
    SharedScratch& sharedScratch(fit.sharedScratch);
    
    uint32_t* numCutsPerVariable = const_cast<uint32_t*>(sharedScratch.numCutsPerVariable);
    double** cutPoints = const_cast<double**>(sharedScratch.cutPoints);
    for (size_t i = 0; i < data.numPredictors; ++i) {
      numCutsPerVariable[i] = static_cast<uint32_t>(-1);
      cutPoints[i] = NULL;
    }
    
    size_t* columns = ext_stackAllocate(data.numPredictors, size_t);
    for (size_t j = 0; j < data.numPredictors; ++j) columns[j] = j;
    
    setCutPoints(fit, columns, data.numPredictors);
    
    ext_stackFree(columns);
  }
  
  void setCutPoints(BARTFit& fit, const size_t* columns, size_t numColumns)
  {
    Control& control(fit.control);
    Data& data(fit.data);
    SharedScratch& sharedScratch(fit.sharedScratch);
    
    uint32_t* numCutsPerVariable = const_cast<uint32_t*>(sharedScratch.numCutsPerVariable);
    double** cutPoints = const_cast<double**>(sharedScratch.cutPoints);
        
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
    
    for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum) {
      if (control.rng_algorithm == EXT_RNG_ALGORITHM_INVALID) {
        // if only one chain, can use environment's rng
        if ((state[chainNum].rng = ext_rng_createDefault(control.numChains == 1)) == NULL) ext_throwError("could not allocate rng");
        
        if (control.rng_standardNormal != EXT_RNG_STANDARD_NORMAL_INVALID &&
            ext_rng_setStandardNormalAlgorithm(state[chainNum].rng, control.rng_standardNormal, NULL) != 0)
          ext_throwError("could not set rng standard normal");
        // if not using envirnoment's rng, we have to seed
        if (control.numChains > 1 && ext_rng_setSeedFromClock(state[chainNum].rng) != 0)
          ext_throwError("could not seed rng");
      } else {
        
        if ((state[chainNum].rng = ext_rng_create(control.rng_algorithm, NULL)) == NULL) ext_throwError("could not allocate rng");
        
        if (control.rng_standardNormal != EXT_RNG_STANDARD_NORMAL_INVALID &&
            ext_rng_setStandardNormalAlgorithm(state[chainNum].rng, control.rng_standardNormal, NULL) != 0)
          ext_throwError("could not set rng standard normal");
        
        if (ext_rng_setSeedFromClock(state[chainNum].rng) != 0)
          ext_throwError("could not seed rng");
      }
    }
  }
  void destroyRNG(BARTFit& fit) {
    for (size_t chainNum = 0; chainNum < fit.control.numChains; ++chainNum)
      ext_rng_destroy(fit.state[chainNum].rng);
  }
  
  void setInitialFit(BARTFit& fit) {
    Control& control(fit.control);
    Data& data(fit.data);
    State* state(fit.state);
    
    for (size_t chainNum = 0; chainNum < control.numChains; ++chainNum) {
      size_t length = data.numObservations * control.numTrees;
      state[chainNum].treeFits = new double[length];
      for (size_t offset = 0; offset < length; ++offset) state[chainNum].treeFits[offset] = 0.0;
      
      state[chainNum].totalFits = new double[data.numObservations];
      for(size_t i = 0; i < data.numObservations; ++i) state[chainNum].totalFits[i] = 0.0;
      
      if (data.numTestObservations > 0) {
        state[chainNum].totalTestFits = new double[data.numTestObservations];
        for (size_t i = 0; i < data.numTestObservations; ++i) state[chainNum].totalTestFits[i] = 0.0;
      }
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
    ext_setVectorToConstant(z, data.numObservations, -1.0);
    if (data.offset != NULL) ext_addVectorsInPlace(data.offset, data.numObservations, -1.0, z);
    ext_addVectorsInPlace(data.y, data.numObservations, 2.0, z);
#else
    // BayesTree initialized the latents to be -2 and 0; was probably a bug
    ext_setVectorToConstant(z, data.numObservations, -2.0);
    if (data.offset != NULL) ext_addVectorsInPlace(data.offset, data.numObservations, -1.0, z);
    ext_addVectorsInPlace(data.y, data.numObservations, 2.0, z);
#endif
  }
  
  void rescaleResponse(BARTFit& fit) {
    const Data& data(fit.data);
    SharedScratch& sharedScratch(fit.sharedScratch);
    
    double* yRescaled = const_cast<double*>(fit.sharedScratch.yRescaled);
    
    if (data.offset != NULL) {
      ext_addVectors(data.offset, data.numObservations, -1.0, data.y, yRescaled);
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
    ext_addScalarToVectorInPlace(   yRescaled, data.numObservations, -sharedScratch.dataScale.min);
    ext_scalarMultiplyVectorInPlace(yRescaled, data.numObservations, 1.0 / sharedScratch.dataScale.range);
    ext_addScalarToVectorInPlace(   yRescaled, data.numObservations, -0.5);
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
                    double sigma, const uint32_t* variableCounts, size_t simNum)
  {
    const Data& data(fit.data);
    const Control& control(fit.control);
    const SharedScratch& sharedScratch(fit.sharedScratch);
    
    size_t chainStride = chainNum * results.numSamples;
    if (control.responseIsBinary) {
      if (control.keepTrainingFits) {
        double* trainingSamples = results.trainingSamples + (simNum + chainStride) * data.numObservations;
        std::memcpy(trainingSamples, trainingSample, data.numObservations * sizeof(double));
        if (data.offset != NULL) ext_addVectorsInPlace(data.offset, data.numObservations, 1.0, trainingSamples);
      }
      
      if (data.numTestObservations > 0) {
        double* testSamples = results.testSamples + (simNum + chainStride) * data.numTestObservations;
        std::memcpy(testSamples, testSample, data.numTestObservations * sizeof(double));
        if (data.testOffset != NULL) ext_addVectorsInPlace(data.testOffset, data.numTestObservations, 1.0, testSamples);
      }
      
      results.sigmaSamples[simNum + chainStride] = 1.0;
      
    } else {
      if (control.keepTrainingFits) {
        double* trainingSamples = results.trainingSamples + (simNum + chainStride) * data.numObservations;
        // set training to dataScale.range * (totalFits + 0.5) + dataScale.min + offset
        ext_setVectorToConstant(trainingSamples, data.numObservations, sharedScratch.dataScale.range * 0.5 + sharedScratch.dataScale.min);
        ext_addVectorsInPlace(trainingSample, data.numObservations, sharedScratch.dataScale.range, trainingSamples);
        if (data.offset != NULL) ext_addVectorsInPlace(data.offset, data.numObservations, 1.0, trainingSamples);
      }
      
      if (data.numTestObservations > 0) {
        double* testSamples = results.testSamples + (simNum + chainStride) * data.numTestObservations;
        ext_setVectorToConstant(testSamples, data.numTestObservations, sharedScratch.dataScale.range * 0.5 + sharedScratch.dataScale.min);
        ext_addVectorsInPlace(testSample, data.numTestObservations, sharedScratch.dataScale.range, testSamples);
        if (data.testOffset != NULL) ext_addVectorsInPlace(data.testOffset, data.numTestObservations, 1.0, testSamples);
      }
       
      results.sigmaSamples[simNum + chainStride] = sigma * sharedScratch.dataScale.range;
    }
    
    double* variableCountSamples = results.variableCountSamples + (simNum + chainStride) * data.numPredictors;
    for (size_t j = 0; j < data.numPredictors; ++j) variableCountSamples[j] = static_cast<double>(variableCounts[j]);
  }
  
  
  void countVariableUses(const BARTFit& fit, const State& state, uint32_t* variableCounts)
  {
    for (size_t treeNum = 0; treeNum < fit.control.numTrees; ++treeNum)
      state.trees[treeNum].countVariableUses(variableCounts);
  }

#ifdef HAVE_GETTIMEOFDAY
  double subtractTimes(struct timeval end, struct timeval start) {
    return (1.0e6 * (static_cast<double>(end.tv_sec - start.tv_sec)) + static_cast<double>(end.tv_usec - start.tv_usec)) / 1.0e6;
  }
#else
  double subtractTimes(time_t end, time_t start) { return static_cast<double>(end - start); }
#endif
  
  void setNewObservationIndices(Node& newNode, size_t* indices, const Node& oldNode)
  {
    newNode.setObservationIndices(indices);
    if (!newNode.isBottom()) {
      setNewObservationIndices(*newNode.getLeftChild(), indices, *oldNode.getLeftChild());
      setNewObservationIndices(*newNode.getRightChild(), indices + oldNode.getLeftChild()->getNumObservations(), *oldNode.getRightChild());
    }
  }
}

#include <external/binaryIO.h>
#include <sys/stat.h> // permissions
#include <fcntl.h>    // open flags
#include <unistd.h>   // unlink
#include "binaryIO.hpp"

#ifndef S_IRGRP
#define S_IRGRP 0
#endif
#ifndef S_IROTH
#define S_IROTH 0
#endif

#define FILE_VERSION_STRING_LENGTH 8
#define FILE_VERSION_STRING "00.09.00"

namespace dbarts {
  
  bool BARTFit::saveToFile(const char* fileName) const
  {
    ext_binaryIO bio;
    int errorCode = ext_bio_initialize(&bio, fileName, O_WRONLY | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
    
    
    if (errorCode != 0) {
      ext_issueWarning("unable to open file: %s", std::strerror(errorCode));
      
      return false;
    }
    
    // because of a peculiarity of how this gets mucked around on creation, this is necessary
    double scaleFactor = control.responseIsBinary ? 1.0 : (data.sigmaEstimate / sharedScratch.dataScale.range);
    double originalScale = model.sigmaSqPrior->getScale();
    model.sigmaSqPrior->setScale(originalScale / (scaleFactor * scaleFactor));
    
    if (ext_bio_writeNChars(&bio, FILE_VERSION_STRING, FILE_VERSION_STRING_LENGTH) != 0) goto save_failed;
    
    if (writeControl(&bio, control) == false) goto save_failed;
    ext_printf("wrote control\n");
    if (writeModel(&bio, model) == false) goto save_failed;
    ext_printf("wrote model\n");
    if (writeData(&bio, data) == false) goto save_failed;
    ext_printf("wrote model\n");
    
    if (writeState(&bio, state, control, data) == false) goto save_failed;
    ext_printf("wrote state\n");
    
    ext_bio_invalidate(&bio);
    
    model.sigmaSqPrior->setScale(originalScale);
    
    printTerminalSummary(*this);
    
    return true;
    
save_failed:
    ext_bio_invalidate(&bio);
    model.sigmaSqPrior->setScale(originalScale);
    unlink(fileName);
    return false; 
  }
  
  
  BARTFit* BARTFit::loadFromFile(const char* fileName) {
    ext_binaryIO bio;
    int errorCode = ext_bio_initialize(&bio, fileName, O_RDONLY, 0);
    if (errorCode != 0) {
      ext_issueWarning("unable to open file: %s", std::strerror(errorCode));
      return NULL;
    }
    
    Version version;
    char* versionString = NULL;
    if ((errorCode = readVersion(&bio, version, &versionString)) != 0) {
      ext_issueWarning("unable to parse version string '%s': %s", versionString, std::strerror(errorCode));
      delete [] versionString;
      return NULL;
    }
    
    Control control;
    Model model;
    Data data;
    BARTFit* result = NULL;
    
    if (readControl(&bio, control, version) == false) goto load_failed;
    ext_printf("read control\n");
    if (readModel(&bio, model) == false) goto load_failed;
    ext_printf("read model\n");
    if (readData(&bio, data) == false) goto load_failed;
    ext_printf("read data\n");
    
    result = new BARTFit(control, model, data);
    
    if (readState(&bio, result->state, result->control, result->data, version) == false) goto load_failed;
    ext_printf("read state\n");
    
    // version 00.08.00 stored running time in state, but it was at the very end so this will work on
    // old objects and new ones
    if ((errorCode = ext_bio_readDouble(&bio, &result->runningTime)) != 0) goto load_failed;
    
    ext_bio_invalidate(&bio);
    
    printTerminalSummary(*result);
    
    return result;
    
load_failed:
    ext_bio_invalidate(&bio);
    
    delete result;
    
    delete [] data.maxNumCuts;
    delete [] data.variableTypes;
    delete [] data.testOffset;
    delete [] data.offset;
    delete [] data.weights;
    delete [] data.x_test;
    delete [] data.x;
    delete [] data.y;
    
    delete model.sigmaSqPrior;
    delete model.muPrior;
    delete model.treePrior;
    
    return NULL;
  }
}
