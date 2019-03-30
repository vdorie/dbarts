#ifndef DBARTS_RESULTS_HPP
#define DBARTS_RESULTS_HPP

#include <cstddef> // size_t

namespace dbarts {
  struct Results {
    double* sigmaSamples;         // 1 x numSamples x numChains
    double* trainingSamples;      // numObservations x numSamples x numChains
    double* testSamples;          // numTestObservations x numSamples x numChains
    double* variableCountSamples; // numPredictors x numSamples x numChains
    double* kSamples;             // 1 x numSamples x numChains, if model over k
    
    std::size_t numObservations;
    std::size_t numPredictors;
    std::size_t numTestObservations;
    std::size_t numSamples;
    std::size_t numChains;
  
    Results(std::size_t numObservations, std::size_t numPredictors,
            std::size_t numTestObservations, std::size_t numSamples, std::size_t numChains,
            bool kIsModeled) :
      sigmaSamples(NULL), trainingSamples(NULL), testSamples(NULL), variableCountSamples(NULL), kSamples(NULL),
      numObservations(numObservations), numPredictors(numPredictors), numTestObservations(numTestObservations),
      numSamples(numSamples), numChains(numChains)
    {
      sigmaSamples = new double[getNumSigmaSamples()];
      trainingSamples = new double[getNumTrainingSamples()];
      if (this->numTestObservations > 0) testSamples = new double[getNumTestSamples()];
      variableCountSamples = new double[getNumVariableCountSamples()];
      if (kIsModeled) kSamples = new double[getNumSigmaSamples()];
    }
    
    // note how dangerous this constructor is, as it accepts pointers
    // but the destructor deletes them; set to NULL before deleting
    // if you use
    Results(std::size_t numObservations, std::size_t numPredictors, std::size_t numTestObservations,
            std::size_t numSamples, std::size_t numChains,
            double* sigmaSamples, double* trainingSamples,
            double* testSamples, double* variableCountSamples, double* kSamples) :
      sigmaSamples(sigmaSamples), trainingSamples(trainingSamples), testSamples(testSamples),
      variableCountSamples(variableCountSamples), kSamples(kSamples), numObservations(numObservations),
      numPredictors(numPredictors), numTestObservations(numTestObservations), numSamples(numSamples),
      numChains(numChains)
    {
    }
    
    ~Results() {
      delete [] kSamples; kSamples = NULL;
      delete [] variableCountSamples; variableCountSamples = NULL;
      delete [] testSamples; testSamples = NULL;
      delete [] trainingSamples; trainingSamples = NULL;
      delete [] sigmaSamples; sigmaSamples = NULL;
    }
    
    std::size_t getNumSigmaSamples() { return numSamples * numChains; }
    std::size_t getNumTrainingSamples() { return numObservations * numSamples * numChains; }
    std::size_t getNumTestSamples() { return numTestObservations * numSamples * numChains; }
    std::size_t getNumVariableCountSamples() { return numPredictors * numSamples * numChains; }
  };
} // namespace dbarts

#endif // DBARTS_RESULTS_HPP

