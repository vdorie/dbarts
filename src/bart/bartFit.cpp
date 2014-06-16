#include "config.hpp"
#include <bart/bartFit.hpp>

#include <cmath>     // sqrt
#include <cstring>   // memcpy
#include <cstddef>   // size_t
#include <bart/cstdint>
#if !defined(HAVE_SYS_TIME_H) && defined(HAVE_GETTIMEOFDAY)
#undef HAVE_GETTIMEOFDAY
#endif
#ifdef HAVE_SYS_TIME_H
#include <sys/time.h> // gettimeofday
#else
#include <time.h>
#endif

#include <set>       // used to sort and find 
#include <vector>    //   split points;
#include <algorithm> // integer min

#include <external/alloca.h>
#include <external/io.h>
#include <external/stats.h>
#include <external/linearAlgebra.h>

#include <bart/results.hpp>
#include "functions.hpp"
#include "tree.hpp"

using std::size_t;
using std::uint32_t;


namespace {
  using namespace bart;

  void allocateMemory(BARTFit& fit);
  void setPrior(BARTFit& fit);
  void setInitialCutPoints(BARTFit& fit);
  void setInitialFit(BARTFit& fit);
  
  void printInitialSummary(const BARTFit& fit);
  void printTerminalSummary(const BARTFit& fit, double runningTime);
  
  void rescaleResponse(BARTFit& fit);
  
  void sampleBinaryOffsets(BARTFit& fit, const double* fits, double* yRescaled);
  void storeSamples(const BARTFit& fit, Results& results, const double* trainingSample, const double* testSample,
                    double sigma, const uint32_t* variableCounts, size_t simNum);
  void countVariableUses(const BARTFit& fit, uint32_t* variableCounts);
  
#ifdef HAVE_GETTIMEOFDAY
  double subtractTimes(struct timeval end, struct timeval start);
#else
  double subtractTimes(time_t end, time_t start);
#endif
}

namespace bart {
  void BARTFit::setResponse(const double* newY) {
    double sigmaUnscaled = state.sigma * scratch.dataScale.range;
    double priorUnscaled = model.sigmaSqPrior->getScale() * scratch.dataScale.range * scratch.dataScale.range;
    
    data.y = newY;
    
    rescaleResponse(*this);
    
    state.sigma = sigmaUnscaled / scratch.dataScale.range;
    model.sigmaSqPrior->setScale(priorUnscaled / (scratch.dataScale.range * scratch.dataScale.range));
    
    
    // rebuild the total fit and tree fits, manually
    ext_setVectorToConstant(state.totalFits, data.numObservations, 0.0);
    for (size_t i = 0; i < control.numTrees; ++i) {
      double* currFits = state.treeFits + i * data.numObservations;
      
      // treeY = y - totalFits
      std::memcpy(scratch.treeY, scratch.yRescaled, data.numObservations * sizeof(double));
      ext_addVectorsInPlace((const double*) state.totalFits, data.numObservations, -1.0, scratch.treeY);
    
      state.trees[i].setNodeAverages(*this, scratch.treeY);
      state.trees[i].getCurrentFits(*this, currFits, NULL);
    
      // totalFits += currFits
      ext_addVectorsInPlace((const double*) currFits, data.numObservations, 1.0, state.totalFits);
    }
  }
  
  BARTFit::BARTFit(Control control, Model model, Data data) :
    control(control), model(model), data(data), scratch(), state(), threadManager(NULL)
  {
    allocateMemory(*this);

    setPrior(*this);
    setInitialCutPoints(*this);
    setInitialFit(*this);

    if (control.verbose) printInitialSummary(*this);
  }
  
  BARTFit::~BARTFit()
  {
    delete [] scratch.yRescaled; scratch.yRescaled = NULL;
    delete [] scratch.Xt; scratch.Xt = NULL;
    delete [] scratch.Xt_test; scratch.Xt_test = NULL;
    delete [] scratch.treeY; scratch.treeY = NULL;
    
    delete [] scratch.numCutsPerVariable; scratch.numCutsPerVariable = NULL;
    if (scratch.cutPoints != NULL) {
      for (size_t i = 0; i < data.numPredictors; ++i) delete [] scratch.cutPoints[i];
    }
    delete [] scratch.cutPoints; scratch.cutPoints = NULL;
    
    if (state.trees != NULL) for (size_t i = 0; i < control.numTrees; ++i) state.trees[i].~Tree();
    ::operator delete(state.trees); state.trees = NULL;
    delete [] state.treeIndices; state.treeIndices = NULL;
    
    delete [] state.treeFits; state.treeFits = NULL;
    delete [] state.totalFits; state.totalFits = NULL;
    if (data.numTestObservations > 0) delete [] state.totalTestFits;
    state.totalTestFits = NULL;
    
    ext_mt_destroy(threadManager);
  }
  
  Results* BARTFit::runSampler()
  {
    return runSampler(control.numBurnIn, control.numSamples);
  }
  
  Results* BARTFit::runSampler(size_t numBurnIn, size_t numSamples)
  {
    bool stepTaken, isThinningIteration;
    StepType ignored;
    
    Results* resultsPointer = new Results(data.numObservations, data.numPredictors,
                                          data.numTestObservations, numSamples);
    Results& results(*resultsPointer);
    
    
    double* currFits = new double[data.numObservations];
    double* currTestFits = NULL;
    if (data.numTestObservations > 0) currTestFits = new double[data.numTestObservations];
    
    uint32_t* variableCounts = ext_stackAllocate(data.numPredictors, uint32_t);
    
    
    size_t totalNumIterations = (numBurnIn + numSamples) * control.treeThinningRate;
    uint32_t majorIterationNum = 0;
    
    if (control.verbose) ext_printf("Running mcmc loop:\n");
    
#ifdef HAVE_GETTIMEOFDAY
    struct timeval startTime;
    struct timeval endTime;
    gettimeofday(&startTime, NULL);
#else
    time_t startTime;
    time_t endTime;
    startTime = time(NULL);
#endif
  
    for (uint32_t k = 0; k < totalNumIterations; ++k) {
      majorIterationNum = k / control.treeThinningRate;
      if (control.verbose && (majorIterationNum + 1) % control.printEvery == 0)
        ext_printf("iteration: %u (of %u)\n", majorIterationNum + 1, totalNumIterations);
      
      isThinningIteration = ((k + 1) % control.treeThinningRate != 0);
      
      if (!isThinningIteration && data.numTestObservations > 0) ext_setVectorToConstant(state.totalTestFits, data.numTestObservations, 0.0);
            
      for (size_t i = 0; i < control.numTrees; ++i) {
        double* oldTreeFits = state.treeFits + i * data.numObservations;
        
        // treeY = y - (totalFits - oldTreeFits)
        // is residual from every *other* tree, so what is left for this tree to do
        std::memcpy(scratch.treeY, scratch.yRescaled, data.numObservations * sizeof(double));
        ext_addVectorsInPlace((const double*) state.totalFits, data.numObservations, -1.0, scratch.treeY);
        ext_addVectorsInPlace((const double*) oldTreeFits, data.numObservations, 1.0, scratch.treeY);
        
        state.trees[i].setNodeAverages(*this, scratch.treeY);
        
        /*if (k == 6 && i == 131) {
          ext_printf("**before:\n");
          trees[i].top.print(*this);
          //ext_printf("  left child obs :\n    ");
          //for (size_t j = 0; j < trees[i].top.leftChild->numObservationsInNode; ++j) ext_printf("%2lu, ", trees[i].top.leftChild->observationIndices[j]);
          //ext_printf("\n  right child obs:\n    ");
          //for (size_t j = 0; j < trees[i].top.rightChild->numObservationsInNode; ++j) ext_printf("%2lu, ", trees[i].top.rightChild->observationIndices[j]);
          //ext_printf("\n");
        }*/
        //ext_printf("iter %lu, tree %lu: ", k + 1, i + 1);
        metropolisJumpForTree(*this, state.trees[i], scratch.treeY, &stepTaken, &ignored);
        /*if (k == 1 && i == 20) {
          ext_printf("**after:\n");
          trees[i].top.print(*this);
          ext_printf("  left child obs :\n    ");
          for (size_t j = 0; j < trees[i].top.leftChild->numObservationsInNode; ++j) ext_printf("%2lu, ", trees[i].top.leftChild->observationIndices[j]);
          ext_printf("\n  right child obs:\n    ");
          for (size_t j = 0; j < trees[i].top.rightChild->numObservationsInNode; ++j) ext_printf("%2lu, ", trees[i].top.rightChild->observationIndices[j]);
          ext_printf("\n");
        } */
        //trees[i].top.print(*this);
        
        state.trees[i].getCurrentFits(*this, currFits, isThinningIteration ? NULL : currTestFits);
        
        // totalFits += currFits - oldTreeFits
        ext_addVectorsInPlace((const double*) oldTreeFits, data.numObservations, -1.0, state.totalFits);
        ext_addVectorsInPlace((const double*) currFits, data.numObservations, 1.0, state.totalFits);
        
        if (!isThinningIteration && data.numTestObservations > 0) {
          ext_addVectorsInPlace((const double*) currTestFits, data.numTestObservations, 1.0, state.totalTestFits);
        }
        
        std::memcpy(oldTreeFits, (const double*) currFits, data.numObservations * sizeof(double));
      }
      
      if (control.responseIsBinary) {
        sampleBinaryOffsets(*this, scratch.yRescaled, state.totalFits);
      } else {
        double sumOfSquaredResiduals = ext_computeAndSumSquaresOfResidualsForVector(scratch.yRescaled, data.numObservations, state.totalFits);
        state.sigma = std::sqrt(model.sigmaSqPrior->drawFromPosterior(data.numObservations, sumOfSquaredResiduals));
      }
      
      if (!isThinningIteration) {
        // if not out of burn-in, store result in first result; start
        // overwriting after that
        bool isBurningIn = majorIterationNum < numBurnIn;
        size_t simNum = (!isBurningIn ? majorIterationNum - numBurnIn : 0);
        
        countVariableUses(*this, variableCounts);
        storeSamples(*this, results, state.totalFits, state.totalTestFits, state.sigma, variableCounts, simNum);
        
        if (control.callback != NULL) {
          control.callback(control.callbackData, *this, isBurningIn,
                           results.trainingSamples + simNum * data.numObservations,
                           results.testSamples + simNum * data.numTestObservations,
                           results.sigmaSamples[simNum]);
        }
      }
    }

#ifdef HAVE_GETTIMEOFDAY
    gettimeofday(&endTime, NULL);
#else
    endTime = time(NULL);
#endif
    
    if (control.verbose) printTerminalSummary(*this, subtractTimes(endTime, startTime));
    
    delete [] currFits;
    if (data.numTestObservations > 0) delete [] currTestFits;
    ext_stackFree(variableCounts);
    
    
    return resultsPointer;
  }
} // namespace bart


namespace {
  using namespace bart;
  
  void printInitialSummary(const BARTFit& fit) {
    const Control& control(fit.control);
    const Data& data(fit.data);
    const Model& model(fit.model);
    const Scratch& scratch(fit.scratch);
    
    if (control.responseIsBinary)
      ext_printf("\nRunning BART with binary y\n\n");
    else
      ext_printf("\nRunning BART with numeric y\n\n");
    
    ext_printf("number of trees: %u\n", control.numTrees);
    
    ext_printf("Prior:\n");
    // dirty hack... should have priors print themselves
    double sigma = std::sqrt(1.0 / static_cast<NormalPrior*>(model.muPrior)->precision);
    double k = (control.responseIsBinary ? 3.0 : 0.5) /  (sigma * std::sqrt((double) control.numTrees));
    ext_printf("\tk: %f\n", k);
    if (!control.responseIsBinary) {
      ChiSquaredPrior* residPrior = static_cast<ChiSquaredPrior*>(model.sigmaSqPrior);
      ext_printf("\tdegrees of freedom in sigma prior: %f\n", residPrior->degreesOfFreedom);
      double quantile = 1.0 - ext_percentileOfChiSquared(residPrior->scale * residPrior->degreesOfFreedom / fit.state.sigma / fit.state.sigma, residPrior->degreesOfFreedom);
      ext_printf("\tquantile in sigma prior: %f\n", quantile);
    }
    CGMPrior* treePrior = static_cast<CGMPrior*>(model.treePrior);
    ext_printf("\tpower and base for tree prior: %f %f\n", treePrior->power, treePrior->base);
    ext_printf("\tuse quantiles for rule cut points: %s\n", control.useQuantiles ? "true" : "false");
    ext_printf("data:\n");
    ext_printf("\tnumber of training observations: %u\n", data.numObservations);
    ext_printf("\tnumber of test observations: %u\n", data.numTestObservations);
    ext_printf("\tnumber of explanatory variables: %u\n\n", data.numPredictors);
    
    
    ext_printf("\nCutoff rules c in x<=c vs x>c\n");
    ext_printf("Number of cutoffs: (var: number of possible c):\n");
    for (size_t i = 0; i < data.numPredictors; ++i ) {
      ext_printf("(%u: %u) ", i + 1, scratch.numCutsPerVariable[i]);
      if ((i + 1) % 5 == 0) ext_printf("\n");
    }
    ext_printf("\n");
    if (control.printCutoffs > 0) {
      ext_printf("cutoffs:\n");
      for (size_t i = 0; i < data.numPredictors; ++i) {
        ext_printf("x(%u) cutoffs: ", i + 1);
        
        size_t j;
        for (j = 0; j < scratch.numCutsPerVariable[i] - 1 && j < control.printCutoffs - 1; ++j) {
          ext_printf("%f", scratch.cutPoints[i][j]);
          if ((j + 1) % 5 == 0) ext_printf("\n\t");
        }
        if (j > 2 && j == control.printCutoffs && j < scratch.numCutsPerVariable[i] - 1)
          ext_printf("...");
        
        ext_printf("%f", scratch.cutPoints[i][scratch.numCutsPerVariable[i] - 1]);
        ext_printf("\n");
      }
    }
  }
  
  void printTerminalSummary(const BARTFit& fit, double runningTime) {
    ext_printf("seconds in loop: %f\n", runningTime);
    
    ext_printf("\nTree sizes, last iteration:\n");
    for (size_t i = 0; i < fit.control.numTrees; ++i) {
      ext_printf("%u ", fit.state.trees[i].getNumBottomNodes());
      if ((i + 1) % 20 == 0) ext_printf("\n");
    }
    ext_printf("\n");
    
    uint32_t* variableCounts = ext_stackAllocate(fit.data.numPredictors, uint32_t);
    
    ext_printf("Variable Usage, last iteration (var:count):\n");
    countVariableUses(fit, variableCounts);
    for (size_t i = 0; i < fit.data.numPredictors; ++i) {
      ext_printf("(%lu: %u) ", (unsigned long int) (i + 1), variableCounts[i]);
      if ((i + 1) % 5 == 0) ext_printf("\n");
    }
    
    ext_stackFree(variableCounts);
    
    
    ext_printf("\nDONE BART\n\n");
  }
  
  void allocateMemory(BARTFit& fit) {
    Control& control(fit.control);
    Data& data(fit.data);
    Scratch& scratch(fit.scratch);
    State& state(fit.state);
        
    scratch.yRescaled = new double[data.numObservations];
    rescaleResponse(fit);
    
    scratch.Xt = new double[data.numObservations * data.numPredictors];
    double* Xt = const_cast<double*>(scratch.Xt);
    for (size_t col = 0; col < data.numPredictors; ++col) {
      for (size_t row = 0; row < data.numObservations; ++row) {
        Xt[row * data.numPredictors + col] = data.X[col * data.numObservations + row];
      }
    }
    
    if (data.numTestObservations > 0) {
      scratch.Xt_test = new double[data.numTestObservations * data.numPredictors];
      double* Xt_test = const_cast<double*>(scratch.Xt_test);
      for (size_t col = 0; col < data.numPredictors; ++col) {
        for (size_t row = 0; row < data.numTestObservations; ++row) {
          Xt_test[row * data.numPredictors + col] = data.X_test[col * data.numTestObservations + row];
        }
      }
    }

    scratch.treeY = new double[data.numObservations];
    for (size_t i = 0; i < data.numObservations; ++i) scratch.treeY[i] = scratch.yRescaled[i];
    
    scratch.numCutsPerVariable = new uint32_t[data.numPredictors];

    scratch.cutPoints = new double*[data.numPredictors];
    const double** cutPoints = const_cast<const double**>(scratch.cutPoints);
    for (size_t i = 0; i < data.numPredictors; ++i) cutPoints[i] = NULL;
    
    state.trees = static_cast<Tree*>(::operator new[] (control.numTrees * sizeof(Tree)));
    state.treeIndices = new size_t[data.numObservations * control.numTrees];
    
    if (control.numThreads > 1 && ext_mt_create(&fit.threadManager, control.numThreads) != 0) {
      ext_printMessage("Unable to multi-thread, defaulting to single.");
    }
  }
  
  void setPrior(BARTFit& fit) {
    Control& control(fit.control);
    Data& data(fit.data);
    Model& model(fit.model);
    Scratch& scratch(fit.scratch);
    State& state(fit.state);
    
    state.sigma = control.responseIsBinary ? 1.0 : (data.sigmaEstimate / scratch.dataScale.range);
    model.sigmaSqPrior->setScale(state.sigma * state.sigma * model.sigmaSqPrior->getScale());
  }
  
  void setInitialCutPoints(BARTFit& fit) {
    Control& control(fit.control);
    Data& data(fit.data);
    Scratch& scratch(fit.scratch);
    
    const double* xCol;
    uint32_t* numCutsPerVariable = const_cast<uint32_t*>(scratch.numCutsPerVariable);
    double** cutPoints = const_cast<double**>(scratch.cutPoints);
    
    if (control.useQuantiles) {
      size_t colNumElements, colFactor, colNumCuts, colIndex, colOffset;
      
      if (data.maxNumCuts == NULL) ext_throwError("num cuts cannot be NULL if useQuantiles is true.");
      
       // sets are inherently sorted, should be a binary tree back there somewhere
      std::set<double> uniqueElements;
      std::vector<double> sortedElements(data.numObservations);
      
      for (size_t i = 0; i < data.numPredictors; ++i) {
        uniqueElements.clear();
        xCol = fit.data.X + i * data.numObservations;
        
        for (size_t j = 0; j < data.numObservations; ++j) uniqueElements.insert(xCol[j]);
        colNumElements = uniqueElements.size();
        
        if (colNumElements <= data.maxNumCuts[i] + 1) {
          colFactor = 1;
          colNumCuts = colNumElements - 1;
          colOffset = 0;
        } else {
          colNumCuts = data.maxNumCuts[i];
          colFactor = colNumCuts / colNumElements;
          colOffset = colFactor / 2;
        }
        
        numCutsPerVariable[i] = colNumCuts;
        cutPoints[i] = new double[colNumCuts];
        
        sortedElements.clear();
        sortedElements.assign(uniqueElements.begin(), uniqueElements.end());
        
        for (size_t j = 0; j < colNumCuts; ++j) {
          colIndex = std::min(j * colFactor + colOffset, colNumElements - 2);
          cutPoints[i][j] = 0.5 * (sortedElements[colIndex] + sortedElements[colIndex + 1]);
        }
      }
    } else {
      double xMax, xMin, xIncrement;
      
      for (size_t i = 0; i < data.numPredictors; ++i) {
        xCol = fit.data.X + i * data.numObservations;
        
        xMax = xCol[0]; xMin = xCol[0];
        for (size_t j = 1; j < data.numObservations; ++j) {
          double x_j = xCol[j];
          if (x_j < xMin) xMin = x_j;
          if (x_j > xMax) xMax = x_j;
        }
        
        numCutsPerVariable[i] = data.maxNumCuts[i];
        cutPoints[i] = new double[numCutsPerVariable[i]];
        
        xIncrement = (xMax - xMin) / (numCutsPerVariable[i] + 1);
        
        for (size_t j = 0; j < numCutsPerVariable[i]; ++j) cutPoints[i][j] = xMin + ((double) (j + 1)) * xIncrement;
      }
    }
  }
  
  void setInitialFit(BARTFit& fit) {
    Control& control(fit.control);
    Data& data(fit.data);
    State& state(fit.state);
    
    for (size_t i = 0; i < control.numTrees; ++i) {
      new (state.trees + i) Tree(state.treeIndices + i * data.numObservations, data.numObservations, data.numPredictors);
    }
    
    size_t length = data.numObservations * control.numTrees;
    state.treeFits = new double[length];
    for (size_t offset = 0; offset < length; ++offset) state.treeFits[offset] = 0.0;
    
    state.totalFits = new double[data.numObservations];
    for(size_t i = 0; i < data.numObservations; ++i) state.totalFits[i] = 0.0;
    
    if (data.numTestObservations > 0) {
      state.totalTestFits = new double[data.numTestObservations];
      for (size_t offset = 0; offset < data.numTestObservations; ++offset) state.totalTestFits[offset] = 0.0;
    }
  }
  
  void rescaleResponse(BARTFit& fit) {
    const Data& data(fit.data);
    const Control& control(fit.control);
    Scratch& scratch(fit.scratch);
    
    double* yRescaled = const_cast<double*>(fit.scratch.yRescaled);
    if (control.responseIsBinary) {
      // yRescaled = 2.0 * y - 1.0 - offset => map to -1 and 1
      ext_setVectorToConstant(yRescaled, data.numObservations, -1.0);
      if (data.offset != NULL) ext_addVectorsInPlace(data.offset, data.numObservations, -1.0, yRescaled);
      ext_addVectorsInPlace((const double*) data.y, data.numObservations, 2.0, yRescaled);
      
      scratch.dataScale.min = 0.0;
      scratch.dataScale.max = 1.0;
      scratch.dataScale.range = 1.0;
    } else {
      scratch.dataScale.min = data.y[0];
      scratch.dataScale.max = data.y[0];
      for (size_t i = 1; i < data.numObservations; ++i) {
        if (data.y[i] < scratch.dataScale.min) scratch.dataScale.min = data.y[i];
        if (data.y[i] > scratch.dataScale.max) scratch.dataScale.max = data.y[i];
      }
      scratch.dataScale.range = scratch.dataScale.max - scratch.dataScale.min;
      
      // yRescaled = -0.5 + (y - min) / (max - min)
      ext_setVectorToConstant(yRescaled, data.numObservations, -0.5 - scratch.dataScale.min / scratch.dataScale.range);
      ext_addVectorsInPlace((const double*) data.y, data.numObservations, 1.0 / scratch.dataScale.range, yRescaled);
    }
  }
  
  void sampleBinaryOffsets(BARTFit& fit, const double* fits, double* yRescaled) {
    double u, zScore;
    
    for (size_t i = 0; i < fit.data.numObservations; ++i) {
      u = ext_simulateContinuousUniform();
      
      double quantile = fits[i];
      if (fit.data.offset != NULL) quantile += fit.data.offset[i];
      
      if (fit.data.y[i] > 0.0) {
        double prob = u + (1.0 - u) * ext_cumulativeProbabilityOfNormal(-quantile, 0.0, 1.0);
        
        zScore =  ext_quantileOfNormal(prob, 0.0, 1.0);
      } else {
        double prob = u + (1.0 - u) * ext_cumulativeProbabilityOfNormal( quantile, 0.0, 1.0);
        
        zScore = -ext_quantileOfNormal(prob, 0.0, 1.0);
      }
      
      yRescaled[i] = fits[i] + zScore;
    }
  }
  
  void storeSamples(const BARTFit& fit, Results& results, const double* trainingSample, const double* testSample,
                    double sigma, const uint32_t* variableCounts, size_t simNum)
  {
    const Data& data(fit.data);
    const Control& control(fit.control);
    const Scratch& scratch(fit.scratch);
    
    size_t offset;
    
    if (control.responseIsBinary) {
      offset = simNum * data.numObservations;
      if (control.keepTrainingFits) {
        double* trainingSamples = results.trainingSamples + offset;
        std::memcpy(trainingSamples, trainingSample, data.numObservations * sizeof(double));
      }
      
      offset = simNum * data.numTestObservations;
      std::memcpy(results.testSamples + offset, testSample, data.numTestObservations * sizeof(double));
      
      results.sigmaSamples[simNum] = 1.0;
      
    } else {
      offset = simNum * data.numObservations;
      if (control.keepTrainingFits) {
        double* trainingSamples = results.trainingSamples + offset;
        // set training to dataScale.range * (totalFits + 0.5) + dataScale.min
        ext_setVectorToConstant(trainingSamples, data.numObservations, scratch.dataScale.range * 0.5 + scratch.dataScale.min);
        ext_addVectorsInPlace(trainingSample, data.numObservations, scratch.dataScale.range, trainingSamples);
      }
      
      offset = simNum * data.numTestObservations;
      double* testSamples = results.testSamples + offset;
      ext_setVectorToConstant(testSamples, data.numTestObservations, scratch.dataScale.range * 0.5 + scratch.dataScale.min);
      ext_addVectorsInPlace(testSample, data.numTestObservations, scratch.dataScale.range, testSamples);
      
      results.sigmaSamples[simNum] = sigma * scratch.dataScale.range;
    }
    
    offset = simNum * data.numPredictors;
    for (size_t i = 0; i < data.numPredictors; ++i) results.variableCountSamples[offset + i] = variableCounts[i];
  }
  
  
  void countVariableUses(const BARTFit& fit, uint32_t* variableCounts)
  {
    for (size_t i = 0; i < fit.data.numPredictors; ++i) variableCounts[i] = 0;
    
    for (size_t i = 0; i < fit.control.numTrees; ++i) {
      fit.state.trees[i].countVariableUses(variableCounts);
    }
  }

#ifdef HAVE_GETTIMEOFDAY
  double subtractTimes(struct timeval end, struct timeval start) {
    return (1.0e6 * ((double) (end.tv_sec - start.tv_sec)) + (double) (end.tv_usec - start.tv_usec)) / 1.0e6;
  }
#else
  double subtractTimes(time_t end, time_t start) { return (double) (end - start); }
#endif
}
