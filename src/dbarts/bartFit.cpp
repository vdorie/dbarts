#include "config.hpp"
#include <dbarts/bartFit.hpp>

#include <cmath>     // sqrt
#include <cstring>   // memcpy
#include <cstddef>   // size_t

#if !defined(HAVE_SYS_TIME_H) && defined(HAVE_GETTIMEOFDAY)
#undef HAVE_GETTIMEOFDAY
#endif
#ifdef HAVE_SYS_TIME_H
#  include <sys/time.h> // gettimeofday
#else
#  include <time.h>
#endif


#include <set>       // used to sort and find 
#include <vector>    //   split points;
#include <algorithm> // integer min

#include <external/alloca.h>
#include <external/io.h>
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
  void setPrior(BARTFit& fit);
  void setInitialCutPoints(BARTFit& fit);
  void setCutPointsFromQuantiles(BARTFit& fit, const double* x, uint32_t maxNumCuts,
                                 uint32_t& numCutsPerVariable, double*& cutPoints,
                                 std::set<double>& uniqueElements, std::vector<double>& sortedElements);
  void setCutPointsUniformly(BARTFit& fit, const double* x, uint32_t maxNumCuts,
                             uint32_t& numCutsPerVariable, double*& cutPoints);
  void setInitialFit(BARTFit& fit);
  
  void printInitialSummary(const BARTFit& fit);
  void printTerminalSummary(const BARTFit& fit);
  
  void initializeLatents(BARTFit& fit);
  void rescaleResponse(BARTFit& fit);
  // void resampleTreeFits(BARTFit& fit);
  
  void sampleProbitLatentVariables(BARTFit& fit, const double* fits, double* yRescaled);
  void storeSamples(const BARTFit& fit, Results& results, const double* trainingSample, const double* testSample,
                    double sigma, const uint32_t* variableCounts, size_t simNum);
  void countVariableUses(const BARTFit& fit, uint32_t* variableCounts);
  
#ifdef HAVE_SYS_TIME_H
  double subtractTimes(struct timeval end, struct timeval start);
#else
  double subtractTimes(time_t end, time_t start);
#endif
}

namespace dbarts {
  
  void BARTFit::setResponse(const double* newY) {
    if (!control.responseIsBinary) {
      double sigmaUnscaled = state.sigma * scratch.dataScale.range;
      double priorUnscaled = model.sigmaSqPrior->getScale() * scratch.dataScale.range * scratch.dataScale.range;
      
      data.y = newY;
      
      rescaleResponse(*this);
      
      state.sigma = sigmaUnscaled / scratch.dataScale.range;
      model.sigmaSqPrior->setScale(priorUnscaled / (scratch.dataScale.range * scratch.dataScale.range));
    } else {
      data.y = newY;
      
      sampleProbitLatentVariables(*this, const_cast<const double*>(state.totalFits), const_cast<double*>(scratch.yRescaled));
    }
    
    // resampleTreeFits(*this);
  }
  
  void BARTFit::setOffset(const double* newOffset) {
    if (!control.responseIsBinary) {
      double sigmaUnscaled = state.sigma * scratch.dataScale.range;
      double priorUnscaled = model.sigmaSqPrior->getScale() * scratch.dataScale.range * scratch.dataScale.range;
      
      data.offset = newOffset;
      
      rescaleResponse(*this);
      
      state.sigma = sigmaUnscaled / scratch.dataScale.range;
      model.sigmaSqPrior->setScale(priorUnscaled / (scratch.dataScale.range * scratch.dataScale.range));
    } else {
      data.offset = newOffset;
      
      sampleProbitLatentVariables(*this, const_cast<const double*>(state.totalFits), const_cast<double*>(scratch.yRescaled));
    }
  }
  
  void BARTFit::setPredictor(const double* newPredictor, size_t predictorColumn) {
    uint32_t& numCutsPerVariable(const_cast<uint32_t*>(scratch.numCutsPerVariable)[predictorColumn]);
    double* cutPoints = const_cast<double**>(scratch.cutPoints)[predictorColumn];
    
    if (control.useQuantiles) {
      if (data.maxNumCuts == NULL) ext_throwError("Num cuts cannot be NULL if useQuantiles is true.");
      
      std::set<double> uniqueElements;
      std::vector<double> sortedElements(data.numObservations);
      
      setCutPointsFromQuantiles(*this, newPredictor, data.maxNumCuts[predictorColumn],
                                numCutsPerVariable, cutPoints,
                                uniqueElements, sortedElements);
    } else {
      setCutPointsUniformly(*this, newPredictor, data.maxNumCuts[predictorColumn],
                            numCutsPerVariable, cutPoints);
    }

    std::memcpy(const_cast<double*>(data.X) + predictorColumn * data.numObservations, newPredictor, data.numObservations * sizeof(double));    
    double* Xt = const_cast<double*>(scratch.Xt);
    for (size_t row = 0; row < data.numObservations; ++row) {
      Xt[row * data.numPredictors + predictorColumn] = newPredictor[row];
    }
  }
  
  void BARTFit::setTestPredictor(const double* newTestPredictor, size_t predictorColumn) {
    std::memcpy(const_cast<double*>(data.X_test) + predictorColumn * data.numTestObservations, newTestPredictor, data.numTestObservations * sizeof(double));    
    double* Xt_test = const_cast<double*>(scratch.Xt_test);
    for (size_t row = 0; row < data.numTestObservations; ++row) {
      Xt_test[row * data.numPredictors + predictorColumn] = newTestPredictor[row];
    }
  }
  
  void BARTFit::setTestPredictors(const double* X_test, size_t numTestObservations) {
    if (numTestObservations == 0 || X_test == NULL) {
      if (scratch.Xt_test != NULL) delete [] scratch.Xt_test; scratch.Xt_test = NULL;
      if (state.totalTestFits != NULL) delete [] state.totalTestFits; state.totalTestFits = NULL;
      
      data.X_test = NULL;
      data.numTestObservations = 0;
    } else {
      data.X_test = X_test;
      
      if (numTestObservations != data.numTestObservations) {
        if (scratch.Xt_test != NULL) delete [] scratch.Xt_test; scratch.Xt_test = NULL;
        if (state.totalTestFits != NULL) delete [] state.totalTestFits; state.totalTestFits = NULL;
        data.numTestObservations = numTestObservations;
        
        scratch.Xt_test = new double[data.numTestObservations * data.numPredictors];
        state.totalTestFits = new double[data.numTestObservations];
      }
      
      double* Xt_test = const_cast<double*>(scratch.Xt_test);
      for (size_t col = 0; col < data.numPredictors; ++col) {
        for (size_t row = 0; row < data.numTestObservations; ++row) {
          Xt_test[row * data.numPredictors + col] = data.X_test[col * data.numTestObservations + row];
        }
      }
    }
  }
  
  BARTFit::BARTFit(Control control, Model model, Data data) :
    control(control), model(model), data(data), scratch(), state(), threadManager(NULL)
  {
    allocateMemory(*this);

    setPrior(*this);
    setInitialCutPoints(*this);
    setInitialFit(*this);

    state.runningTime = 0.0;

    if (this->control.verbose) printInitialSummary(*this);
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
    
    if (state.trees != NULL) for (size_t i = control.numTrees; i > 0; ) state.trees[--i].~Tree();
    ::operator delete (state.trees); state.trees = NULL;
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
    
    double numEffectiveObservations = 
      data.weights == NULL ? (double) data.numObservations : ext_sumVectorElements(data.weights, data.numObservations);
    
    
    double* currFits = new double[data.numObservations];
    double* currTestFits = NULL;
    if (data.numTestObservations > 0) currTestFits = new double[data.numTestObservations];
    
    uint32_t* variableCounts = ext_stackAllocate(data.numPredictors, uint32_t);
    
    
    size_t totalNumIterations = (numBurnIn + numSamples) * control.treeThinningRate;
    uint32_t majorIterationNum = 0;
    
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
    
    for (uint32_t k = 0; k < totalNumIterations; ++k) {
      isThinningIteration = ((k + 1) % control.treeThinningRate != 0);
            
      majorIterationNum = k / control.treeThinningRate;
      
      if (control.verbose && !isThinningIteration && (majorIterationNum + 1) % control.printEvery == 0)
        ext_printf("iteration: %u (of %u)\n", majorIterationNum + 1, totalNumIterations / control.treeThinningRate);
      
      if (!isThinningIteration && data.numTestObservations > 0) ext_setVectorToConstant(state.totalTestFits, data.numTestObservations, 0.0);
      

      for (size_t i = 0; i < control.numTrees; ++i) {
        double* oldTreeFits = state.treeFits + i * data.numObservations;
        
        // treeY = y - (totalFits - oldTreeFits)
        // is residual from every *other* tree, so what is left for this tree to do
        std::memcpy(scratch.treeY, scratch.yRescaled, data.numObservations * sizeof(double));
        ext_addVectorsInPlace((const double*) state.totalFits, data.numObservations, -1.0, scratch.treeY);
        ext_addVectorsInPlace((const double*) oldTreeFits, data.numObservations, 1.0, scratch.treeY);
        
        state.trees[i].setNodeAverages(*this, scratch.treeY);
        
        /*if (k == 1 && i <= 1) {
          ext_printf("**before:\n");
          state.trees[i].top.print(*this);
          if (!state.trees[i].top.isBottom()) {
            ext_printf("  left child obs :\n    ");
            for (size_t j = 0; j < state.trees[i].top.getLeftChild()->getNumObservations(); ++j) ext_printf("%2lu, ", state.trees[i].top.getLeftChild()->observationIndices[j]);
            ext_printf("\n  right child obs:\n    ");
            for (size_t j = 0; j < state.trees[i].top.getRightChild()->getNumObservations(); ++j) ext_printf("%2lu, ", state.trees[i].top.getRightChild()->observationIndices[j]);
            ext_printf("\n");
          }
        } */
        // ext_printf("iter %lu, tree %lu: ", k + 1, i + 1);
        metropolisJumpForTree(*this, state.trees[i], scratch.treeY, &stepTaken, &ignored);
        /* if (k == 1 && i <= 3) {
         ext_printf("**after:\n");
          state.trees[i].top.print(*this);
          if (!state.trees[i].top.isBottom()) {
            ext_printf("  left child obs :\n    ");
            for (size_t j = 0; j < state.trees[i].top.getLeftChild()->getNumObservations(); ++j) ext_printf("%2lu, ", state.trees[i].top.getLeftChild()->observationIndices[j]);
            ext_printf("\n  right child obs:\n    ");
            for (size_t j = 0; j < state.trees[i].top.getRightChild()->getNumObservations(); ++j) ext_printf("%2lu, ", state.trees[i].top.getRightChild()->observationIndices[j]);
          }
          ext_printf("\n");
        } */
        // state.trees[i].top.print(*this);
        
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
        sampleProbitLatentVariables(*this, state.totalFits, const_cast<double*>(scratch.yRescaled));
        // resampleTreeFits(*this);
      } else {
        double sumOfSquaredResiduals;
        if (data.weights != NULL) {
          sumOfSquaredResiduals = ext_mt_computeWeightedSumOfSquaredResiduals(threadManager, scratch.yRescaled, data.numObservations, data.weights, state.totalFits);
        } else {
          sumOfSquaredResiduals = ext_mt_computeSumOfSquaredResiduals(threadManager, scratch.yRescaled, data.numObservations, state.totalFits);
        }
        // double sumOfSquaredResiduals = ext_computeAndSumSquaresOfResidualsForVector(scratch.yRescaled, data.numObservations, state.totalFits);
        state.sigma = std::sqrt(model.sigmaSqPrior->drawFromPosterior(numEffectiveObservations, sumOfSquaredResiduals));
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
    
#ifdef HAVE_SYS_TIME_H
    gettimeofday(&endTime, NULL);
#else
    endTime = time(NULL);
#endif
    
    state.runningTime += subtractTimes(endTime, startTime);
    
    if (control.verbose) printTerminalSummary(*this);
    
    delete [] currFits;
    if (data.numTestObservations > 0) delete [] currTestFits;
    ext_stackFree(variableCounts);
    
    
    return resultsPointer;
  }
} // namespace dbarts


namespace {
  using namespace dbarts;
  
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
    if (data.weights != NULL) ext_printf("\tusing observation weights\n");
    
    
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
  
  void printTerminalSummary(const BARTFit& fit) {
    ext_printf("total seconds in loop: %f\n", fit.state.runningTime);
    
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
    
    if (control.responseIsBinary) initializeLatents(fit);
    else rescaleResponse(fit);
    
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
    
    state.trees = static_cast<Tree*>(::operator new (control.numTrees * sizeof(Tree)));
    state.treeIndices = new size_t[data.numObservations * control.numTrees];
    
    
    for (size_t i = 0; i < control.numTrees; ++i) {
      new (state.trees + i) Tree(state.treeIndices + i * data.numObservations, data.numObservations, data.numPredictors);
    }
    
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
    
    uint32_t* numCutsPerVariable = const_cast<uint32_t*>(scratch.numCutsPerVariable);
    double** cutPoints = const_cast<double**>(scratch.cutPoints);
    for (size_t i = 0; i < data.numPredictors; ++i) {
      numCutsPerVariable[i] = (uint32_t) -1;
      cutPoints[i] = NULL;
    }
    
    if (control.useQuantiles) {
      if (data.maxNumCuts == NULL) ext_throwError("Num cuts cannot be NULL if useQuantiles is true.");
      
       // sets are inherently sorted, should be a binary tree back there somewhere
      std::set<double> uniqueElements;
      std::vector<double> sortedElements(data.numObservations);
      
      for (size_t i = 0; i < data.numPredictors; ++i) {
        setCutPointsFromQuantiles(fit, data.X + i * data.numObservations, data.maxNumCuts[i],
                                  numCutsPerVariable[i], cutPoints[i],
                                  uniqueElements, sortedElements);
      }
    } else {
      for (size_t i = 0; i < data.numPredictors; ++i) {
        setCutPointsUniformly(fit, data.X + i * data.numObservations, data.maxNumCuts[i],
                              numCutsPerVariable[i], cutPoints[i]);
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
      numCuts = numUniqueElements - 1;
      step = 1;
      offset = 0;
    } else {
      numCuts = maxNumCuts;
      step = numCuts / numUniqueElements;
      offset = step / 2;
    }
    
    if (numCutsPerVariable != (uint32_t) -1) {
      if (numCuts < numCutsPerVariable) ext_throwError("Number of induced cut points in new predictor less than previous: old splits would be invalid.");
      if (numCuts > numCutsPerVariable) ext_issueWarning("Number of induced cut points in new predictor greater than previous: ignoring extra quantiles.");
    } else {
      numCutsPerVariable = (uint32_t) numCuts;
      cutPoints = new double[numCuts];
    }
    
    sortedElements.clear();
    sortedElements.assign(uniqueElements.begin(), uniqueElements.end());
      
    for (size_t i = 0; i < numCutsPerVariable; ++i) {
      size_t index = std::min(i * step + offset, numUniqueElements - 2);
      cutPoints[i] = 0.5 * (sortedElements[index] + sortedElements[index + 1]);
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
    
    if (numCutsPerVariable == (uint32_t) -1) {
      numCutsPerVariable = maxNumCuts;
      cutPoints = new double[numCutsPerVariable];
    }
      
    xIncrement = (xMax - xMin) / (double) (numCutsPerVariable + 1);
      
    for (size_t i = 0; i < numCutsPerVariable; ++i) cutPoints[i] = xMin + ((double) (i + 1)) * xIncrement;
  }
  
  void setInitialFit(BARTFit& fit) {
    Control& control(fit.control);
    Data& data(fit.data);
    State& state(fit.state);
    
    size_t length = data.numObservations * control.numTrees;
    state.treeFits = new double[length];
    for (size_t offset = 0; offset < length; ++offset) state.treeFits[offset] = 0.0;
    
    state.totalFits = new double[data.numObservations];
    for(size_t i = 0; i < data.numObservations; ++i) state.totalFits[i] = 0.0;
    
    if (data.numTestObservations > 0) {
      state.totalTestFits = new double[data.numTestObservations];
      for (size_t i = 0; i < data.numTestObservations; ++i) state.totalTestFits[i] = 0.0;
    }
  }
  
  void initializeLatents(BARTFit& fit) {
    const Data& data(fit.data);
    Scratch& scratch(fit.scratch);
    
    double* z = const_cast<double*>(fit.scratch.yRescaled);
    
    // z = 2.0 * y - 1.0 - offset; so -1 if y == 0 and 1 if y == 1 when offset == 0
#ifndef MATCH_BAYES_TREE
    ext_setVectorToConstant(z, data.numObservations, -1.0);
    if (data.offset != NULL) ext_addVectorsInPlace(data.offset, data.numObservations, -1.0, z);
    ext_addVectorsInPlace(data.y, data.numObservations, 2.0, z);
    
    // shouldn't be used, but will leave at reasonable values; if anyone cares, should
    // look at offset var for min/max/range
    scratch.dataScale.min = -1.0;
    scratch.dataScale.max =  1.0;
    scratch.dataScale.range = 2.0;
#else
    // BayesTree initialized the latents to be -2 and 0; was probably a bug
    ext_setVectorToConstant(z, data.numObservations, -2.0);
    if (data.offset != NULL) ext_addVectorsInPlace(data.offset, data.numObservations, -1.0, z);
    ext_addVectorsInPlace(data.y, data.numObservations, 2.0, z);
    
    scratch.dataScale.min = -2.0;
    scratch.dataScale.max =  0.0;
    scratch.dataScale.range = 2.0;
#endif
  }
  
  void rescaleResponse(BARTFit& fit) {
    const Data& data(fit.data);
    Scratch& scratch(fit.scratch);
    
    double* yRescaled = const_cast<double*>(fit.scratch.yRescaled);
    
    if (data.offset != NULL) {
      ext_addVectors(data.offset, data.numObservations, -1.0, data.y, yRescaled);
    } else {
      std::memcpy(yRescaled, data.y, data.numObservations * sizeof(double));
    }
    
    scratch.dataScale.min = yRescaled[0];
    scratch.dataScale.max = yRescaled[0];
    for (size_t i = 1; i < data.numObservations; ++i) {
      if (yRescaled[i] < scratch.dataScale.min) scratch.dataScale.min = yRescaled[i];
      if (yRescaled[i] > scratch.dataScale.max) scratch.dataScale.max = yRescaled[i];
    }
    scratch.dataScale.range = scratch.dataScale.max - scratch.dataScale.min;
    
    // yRescaled = (y - offset - min) / (max - min) - 0.5
    ext_addScalarToVectorInPlace(   yRescaled, data.numObservations, -scratch.dataScale.min);
    ext_scalarMultiplyVectorInPlace(yRescaled, data.numObservations, 1.0 / scratch.dataScale.range);
    ext_addScalarToVectorInPlace(   yRescaled, data.numObservations, -0.5);
  }
  
  /* void resampleTreeFits(BARTFit& fit) {
    const Data& data(fit.data);
    const Control& control(fit.control);
    State& state(fit.state);
    Scratch& scratch(fit.scratch);
    
    // rebuild the total fit and tree fits, manually
    ext_setVectorToConstant(state.totalFits, data.numObservations, 0.0);
    for (size_t i = 0; i < control.numTrees; ++i) {
      double* currFits = state.treeFits + i * data.numObservations;
      
      // treeY = y - totalFits
      std::memcpy(scratch.treeY, scratch.yRescaled, data.numObservations * sizeof(double));
      ext_addVectorsInPlace((const double*) state.totalFits, data.numObservations, -1.0, scratch.treeY);
      
      state.trees[i].setNodeAverages(fit, scratch.treeY);
      state.trees[i].getCurrentFits(fit, currFits, NULL);
      
      // totalFits += currFits
      ext_addVectorsInPlace((const double*) currFits, data.numObservations, 1.0, state.totalFits);
    }
  } */
  
  // multithread-this!
  void sampleProbitLatentVariables(BARTFit& fit, const double* fits, double* z) {
    for (size_t i = 0; i < fit.data.numObservations; ++i) {      
#ifndef MATCH_BAYES_TREE
      double mean = fits[i];
      double offset = 0.0;
      if (fit.data.offset != NULL) offset = fit.data.offset[i];
      
      if (fit.data.y[i] > 0.0) {
        z[i] = ext_simulateLowerTruncatedNormalScale1(mean, offset);
      } else {
        z[i] = ext_simulateUpperTruncatedNormalScale1(mean, offset);
      }
#else
      double prob;
      
      double mean = fits[i];
      if (fit.data.offset != NULL) mean += fit.data.offset[i];
      
      double u = ext_simulateContinuousUniform();
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
  
  void storeSamples(const BARTFit& fit, Results& results, const double* trainingSample, const double* testSample,
                    double sigma, const uint32_t* variableCounts, size_t simNum)
  {
    const Data& data(fit.data);
    const Control& control(fit.control);
    const Scratch& scratch(fit.scratch);
    
    size_t sampleOffset;
    
    if (control.responseIsBinary) {
      sampleOffset = simNum * data.numObservations;
      if (control.keepTrainingFits) {
        double* trainingSamples = results.trainingSamples + sampleOffset;
        std::memcpy(trainingSamples, trainingSample, data.numObservations * sizeof(double));
      }
      
      sampleOffset = simNum * data.numTestObservations;
      std::memcpy(results.testSamples + sampleOffset, testSample, data.numTestObservations * sizeof(double));
      
      results.sigmaSamples[simNum] = 1.0;
      
    } else {
      sampleOffset = simNum * data.numObservations;
      if (control.keepTrainingFits) {
        double* trainingSamples = results.trainingSamples + sampleOffset;
        // set training to dataScale.range * (totalFits + 0.5) + dataScale.min + offset
        ext_setVectorToConstant(trainingSamples, data.numObservations, scratch.dataScale.range * 0.5 + scratch.dataScale.min);
        ext_addVectorsInPlace(trainingSample, data.numObservations, scratch.dataScale.range, trainingSamples);
        // if (data.offset != NULL) ext_addVectorsInPlace(data.offset, data.numObservations, 1.0, trainingSamples);
      }
      
      sampleOffset = simNum * data.numTestObservations;
      double* testSamples = results.testSamples + sampleOffset;
      ext_setVectorToConstant(testSamples, data.numTestObservations, scratch.dataScale.range * 0.5 + scratch.dataScale.min);
      ext_addVectorsInPlace(testSample, data.numTestObservations, scratch.dataScale.range, testSamples);
      
      results.sigmaSamples[simNum] = sigma * scratch.dataScale.range;
    }
    
    sampleOffset = simNum * data.numPredictors;
    for (size_t i = 0; i < data.numPredictors; ++i) results.variableCountSamples[sampleOffset + i] = variableCounts[i];
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
