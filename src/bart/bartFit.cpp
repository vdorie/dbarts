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
    double sigmaUnscaled = sigma * dataScale.range;
    double priorUnscaled = model.sigmaSqPrior->getScale() * dataScale.range * dataScale.range;
    
    data.y = newY;
    
    rescaleResponse(*this);
    
    sigma = sigmaUnscaled / dataScale.range;
    model.sigmaSqPrior->setScale(priorUnscaled / (dataScale.range * dataScale.range));
    
    
    // rebuild the total fit and tree fits, manually
    ext_setVectorToConstant(totalFits, data.numObservations, 0.0);
    for (size_t i = 0; i < control.numTrees; ++i) {
      double* currFits = treeFits + i * data.numObservations;
      
      // treeY = y - totalFits
      std::memcpy(treeY, (const double*) yRescaled, data.numObservations * sizeof(double));
      ext_addVectorsInPlace((const double*) totalFits, data.numObservations, -1.0, treeY);
    
      trees[i].setNodeAverages(*this, treeY);
      trees[i].getCurrentFits(*this, currFits, NULL);
    
      // totalFits += currFits
      ext_addVectorsInPlace((const double*) currFits, data.numObservations, 1.0, totalFits);
    }
  }
  
  BARTFit::BARTFit(Control control, Model model, Data data) :
    control(control), model(model), data(data),
    yRescaled(NULL), Xt(NULL), Xt_test(NULL), treeY(NULL), weights(NULL),
    variableTypes(NULL), numCutsPerVariable(NULL), cutPoints(NULL), trees(NULL),
    treeIndices(NULL), treeFits(NULL), totalFits(NULL), totalTestFits(NULL), threadManager(NULL)
  {
    allocateMemory(*this);

    setPrior(*this);
    setInitialCutPoints(*this);
    setInitialFit(*this);

    if (control.verbose) printInitialSummary(*this);
  }
  
  BARTFit::~BARTFit()
  {
    delete [] yRescaled; yRescaled = NULL;
    delete [] Xt; Xt = NULL;
    delete [] Xt_test; Xt_test = NULL;
    delete [] treeY; treeY = NULL;
    delete [] variableTypes; variableTypes = NULL;
    delete [] weights; weights = NULL;
    
    delete [] numCutsPerVariable; numCutsPerVariable = NULL;
    if (cutPoints != NULL) {
      for (size_t i = 0; i < data.numPredictors; ++i) delete [] cutPoints[i];
    }
    delete [] cutPoints; cutPoints = NULL;
    
    if (trees != NULL) for (size_t i = 0; i < control.numTrees; ++i) trees[i].~Tree();
    ::operator delete(trees); trees = NULL;
    delete [] treeIndices; treeIndices = NULL;
    
    delete [] treeFits; treeFits = NULL;
    delete [] totalFits; totalFits = NULL;
    if (data.numTestObservations > 0) delete [] totalTestFits;
    totalTestFits = NULL;
    
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
      
      if (!isThinningIteration && data.numTestObservations > 0) ext_setVectorToConstant(totalTestFits, data.numTestObservations, 0.0);
            
      for (size_t i = 0; i < control.numTrees; ++i) {
        double* oldTreeFits = treeFits + i * data.numObservations;
        
        // treeY = y - (totalFits - oldTreeFits)
        // is residual from every *other* tree, so what is left for this tree to do
        std::memcpy(treeY, (const double*) yRescaled, data.numObservations * sizeof(double));
        ext_addVectorsInPlace((const double*) totalFits, data.numObservations, -1.0, treeY);
        ext_addVectorsInPlace((const double*) oldTreeFits, data.numObservations, 1.0, treeY);
        
        trees[i].setNodeAverages(*this, treeY);
        
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
        metropolisJumpForTree(*this, trees[i], treeY, &stepTaken, &ignored);
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
        
        trees[i].getCurrentFits(*this, currFits, isThinningIteration ? NULL : currTestFits);
        
        // totalFits += currFits - oldTreeFits
        ext_addVectorsInPlace((const double*) oldTreeFits, data.numObservations, -1.0, totalFits);
        ext_addVectorsInPlace((const double*) currFits, data.numObservations, 1.0, totalFits);
        
        if (!isThinningIteration && data.numTestObservations > 0) {
          ext_addVectorsInPlace((const double*) currTestFits, data.numTestObservations, 1.0, totalTestFits);
        }
        
        std::memcpy(oldTreeFits, (const double*) currFits, data.numObservations * sizeof(double));
      }
      
      if (control.responseIsBinary) {
        sampleBinaryOffsets(*this, yRescaled, totalFits);
      } else {
        double sumOfSquaredResiduals = ext_computeAndSumSquaresOfResidualsForVector(yRescaled, data.numObservations, totalFits);
        sigma = std::sqrt(model.sigmaSqPrior->drawFromPosterior(data.numObservations, sumOfSquaredResiduals));
      }
      
      if (!isThinningIteration) {
        // if not out of burn-in, store result in first result; start
        // overwriting after that
        bool isBurningIn = majorIterationNum < numBurnIn;
        size_t simNum = (!isBurningIn ? majorIterationNum - numBurnIn : 0);
        
        countVariableUses(*this, variableCounts);
        storeSamples(*this, results, totalFits, totalTestFits, sigma, variableCounts, simNum);
        
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
    
    if (control.responseIsBinary)
      ext_printf("\n\nRunning BART with binary y\n\n");
    else
      ext_printf("\n\nRunning BART with numeric y\n\n");
    
    ext_printf("number of trees: %u\n", control.numTrees);
    
    ext_printf("Prior:\n");
    ext_printf("\tk: %f\n", control.kFactor);
    if (control.responseIsBinary) {
      ext_printf("\tbinary offset is: %lf\n", control.binaryOffset);
    } else {
      ext_printf("\tdegrees of freedom in sigma prior: %u\n", control.sigmaDf);
      ext_printf("\tquantile in sigma prior: %f\n", control.sigmaQuantile);
    }
    ext_printf("\tpower and base for tree prior: %f %f\n", control.power, control.base);
    ext_printf("\tuse quantiles for rule cut points: %s\n", control.useQuantiles ? "true" : "false");
    ext_printf("data:\n");
    ext_printf("\tnumber of training observations: %u\n", data.numObservations);
    ext_printf("\tnumber of test observations: %u\n", data.numTestObservations);
    ext_printf("\tnumber of explanatory variables: %u\n\n", data.numPredictors);
    
    
    ext_printf("\nCutoff rules c in x<=c vs x>c\n");
    ext_printf("Number of cutoffs: (var: number of possible c):\n");
    for (size_t i = 0; i < data.numPredictors; ++i ) {
      ext_printf("(%u: %u) ", i + 1, fit.numCutsPerVariable[i]);
      if ((i + 1) % 5 == 0) ext_printf("\n");
    }
    ext_printf("\n");
    if (control.printCutoffs > 0) {
      ext_printf("cutoffs:\n");
      for (size_t i = 0; i < data.numPredictors; ++i) {
        ext_printf("x(%u) cutoffs: ", i + 1);
        
        size_t j;
        for (j = 0; j < fit.numCutsPerVariable[i] - 1 && j < control.printCutoffs - 1; ++j) {
          ext_printf("%f", fit.cutPoints[i][j]);
          if ((j + 1) % 5 == 0) ext_printf("\n\t");
        }
        if (j > 2 && j == control.printCutoffs && j < fit.numCutsPerVariable[i] - 1)
          ext_printf("...");
        
        ext_printf("%f", fit.cutPoints[i][fit.numCutsPerVariable[i] - 1]);
        ext_printf("\n");
      }
    }
    ext_printf("\n\n");
  }
  
  void printTerminalSummary(const BARTFit& fit, double runningTime) {
    ext_printf("seconds in loop: %f\n", runningTime);
    
    ext_printf("\nTree sizes, last iteration:\n");
    for (size_t i = 0; i < fit.control.numTrees; ++i) {
      ext_printf("%u ", fit.trees[i].getNumBottomNodes());
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
        
    fit.yRescaled = new double[data.numObservations];
    rescaleResponse(fit);
    
    fit.Xt = new double[data.numObservations * data.numPredictors];
    for (size_t col = 0; col < data.numPredictors; ++col) {
      for (size_t row = 0; row < data.numObservations; ++row) {
        fit.Xt[row * data.numPredictors + col] = data.X[col * data.numObservations + row];
      }
    }
    
    if (data.numTestObservations > 0) {
      fit.Xt_test = new double[data.numTestObservations * data.numPredictors];
      for (size_t col = 0; col < data.numPredictors; ++col) {
        for (size_t row = 0; row < data.numTestObservations; ++row) {
          fit.Xt_test[row * data.numPredictors + col] = data.X_test[col * data.numTestObservations + row];
        }
      }
    }

    fit.treeY = new double[data.numObservations];
    for (size_t i = 0; i < data.numObservations; ++i) fit.treeY[i] = fit.yRescaled[i];
    
    fit.variableTypes = new VariableType[data.numPredictors];
    for (size_t i = 0; i < data.numPredictors; ++i) fit.variableTypes[i] = ORDINAL;
    
    fit.weights = new double[data.numObservations];
    for (size_t i = 0; i < data.numObservations; ++i) fit.weights[i] = 1.0;
    
    fit.numCutsPerVariable = new uint32_t[data.numPredictors];

    fit.cutPoints = new double*[data.numPredictors];
    for (size_t i = 0; i < data.numPredictors; ++i) fit.cutPoints[i] = NULL;
    
    fit.trees = static_cast<Tree*>(::operator new[] (control.numTrees * sizeof(Tree)));
    fit.treeIndices = new size_t[data.numObservations * control.numTrees];
    
    if (control.numThreads > 1 && ext_mt_create(&fit.threadManager, control.numThreads) != 0) {
      ext_printMessage("Unable to multi-thread, defaulting to single.");
    }
  }
  
  void setPrior(BARTFit& fit) {
    Control& control(fit.control);
    Data& data(fit.data);
    
    fit.sigma = control.responseIsBinary ? 1.0 : (data.sigmaEstimate / fit.dataScale.range);
    fit.model.sigmaSqPrior->setScale(fit.sigma * fit.sigma * fit.model.sigmaSqPrior->getScale());
  }
  
  void setInitialCutPoints(BARTFit& fit) {
    Control& control(fit.control);
    Data& data(fit.data);
    
    const double* xCol;
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
        
        fit.numCutsPerVariable[i] = colNumCuts;
        fit.cutPoints[i] = new double[colNumCuts];
        
        sortedElements.clear();
        sortedElements.assign(uniqueElements.begin(), uniqueElements.end());
        
        for (size_t j = 0; j < colNumCuts; ++j) {
          colIndex = std::min(j * colFactor + colOffset, colNumElements - 2);
          fit.cutPoints[i][j] = 0.5 * (sortedElements[colIndex] + sortedElements[colIndex + 1]);
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
        
        fit.numCutsPerVariable[i] = data.maxNumCuts[i];
        fit.cutPoints[i] = new double[fit.numCutsPerVariable[i]];
        
        xIncrement = (xMax - xMin) / (fit.numCutsPerVariable[i] + 1);
        
        for (size_t j = 0; j < fit.numCutsPerVariable[i]; ++j) fit.cutPoints[i][j] = xMin + ((double) (j + 1)) * xIncrement;
      }
    }
  }
  
  void setInitialFit(BARTFit& fit) {
    Control& control(fit.control);
    Data& data(fit.data);
    
    for (size_t i = 0; i < control.numTrees; ++i) {
      new (fit.trees + i) Tree(fit.treeIndices + i * data.numObservations, fit.data.numObservations, fit.data.numPredictors);
    }
    
    size_t length = data.numObservations * control.numTrees;
    fit.treeFits = new double[length];
    for (size_t offset = 0; offset < length; ++offset) fit.treeFits[offset] = 0.0;
    
    fit.totalFits = new double[data.numObservations];
    for(size_t i = 0; i < data.numObservations; ++i) fit.totalFits[i] = 0.0;
    
    if (data.numTestObservations > 0) {
      fit.totalTestFits = new double[data.numTestObservations];
      for (size_t offset = 0; offset < data.numTestObservations; ++offset) fit.totalTestFits[offset] = 0.0;
    }
  }
  
  void rescaleResponse(BARTFit& fit) {
    Data& data(fit.data);
    Control& control(fit.control);
        
    if (control.responseIsBinary) {
      // yRescaled = 2.0 * y - 1.0 - offset => map to -1 and 1
      ext_setVectorToConstant(fit.yRescaled, data.numObservations, -1.0 - control.binaryOffset);
      ext_addVectorsInPlace((const double*) data.y, data.numObservations, 2.0, fit.yRescaled);
      
      fit.dataScale.min = 0.0;
      fit.dataScale.max = 1.0;
      fit.dataScale.range = 1.0;
    } else {
      fit.dataScale.min = data.y[0];
      fit.dataScale.max = data.y[0];
      for (size_t i = 1; i < data.numObservations; ++i) {
        if (data.y[i] < fit.dataScale.min) fit.dataScale.min = data.y[i];
        if (data.y[i] > fit.dataScale.max) fit.dataScale.max = data.y[i];
      }
      fit.dataScale.range = fit.dataScale.max - fit.dataScale.min;
      
      // yRescaled = -0.5 + (y - min) / (max - min)
      ext_setVectorToConstant(fit.yRescaled, data.numObservations, -0.5 - fit.dataScale.min / fit.dataScale.range);
      ext_addVectorsInPlace((const double*) data.y, data.numObservations, 1.0 / fit.dataScale.range, fit.yRescaled);
    }
  }
  
  void sampleBinaryOffsets(BARTFit& fit, const double* fits, double* yRescaled) {
    double u, zScore;
    
    for (size_t i = 0; i < fit.data.numObservations; ++i) {
      u = ext_simulateContinuousUniform();
      
      double quantile = fits[i] + fit.control.binaryOffset;
      
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
      const ScaleFactor& dataScale(fit.dataScale);
      
      offset = simNum * data.numObservations;
      if (control.keepTrainingFits) {
        double* trainingSamples = results.trainingSamples + offset;
        // set training to dataScale.range * (totalFits + 0.5) + dataScale.min
        ext_setVectorToConstant(trainingSamples, data.numObservations, dataScale.range * 0.5 + dataScale.min);
        ext_addVectorsInPlace(trainingSample, data.numObservations, dataScale.range, trainingSamples);
      }
      
      offset = simNum * data.numTestObservations;
      double* testSamples = results.testSamples + offset;
      ext_setVectorToConstant(testSamples, data.numTestObservations, dataScale.range * 0.5 + dataScale.min);
      ext_addVectorsInPlace(testSample, data.numTestObservations, dataScale.range, testSamples);
      
      results.sigmaSamples[simNum] = sigma * dataScale.range;
    }
    
    offset = simNum * data.numPredictors;
    for (size_t i = 0; i < data.numPredictors; ++i) results.variableCountSamples[offset + i] = variableCounts[i];
  }
  
  
  void countVariableUses(const BARTFit& fit, uint32_t* variableCounts)
  {
    for (size_t i = 0; i < fit.data.numPredictors; ++i) variableCounts[i] = 0;
    
    for (size_t i = 0; i < fit.control.numTrees; ++i) {
      fit.trees[i].countVariableUses(variableCounts);
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
