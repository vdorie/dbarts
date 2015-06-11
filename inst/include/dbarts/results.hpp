#ifndef DBARTS_RESULTS_HPP
#define DBARTS_RESULTS_HPP

#include <cstddef> // size_t

namespace dbarts {
  struct Results {
    double* sigmaSamples;         // 1 x numSamples
    double* trainingSamples;      // numObservations x numSamples
    double* testSamples;          // numTestObservations x numSamples
    double* variableCountSamples; // numPredictors x numSamples
    
    std::size_t numObservations;
    std::size_t numPredictors;
    std::size_t numTestObservations;
    std::size_t numSamples;
  
    Results(std::size_t numObservations, std::size_t numPredictors,
            std::size_t numTestObservations, std::size_t numSamples) :
      sigmaSamples(NULL), trainingSamples(NULL), testSamples(NULL), variableCountSamples(NULL),
      numObservations(numObservations), numPredictors(numPredictors), numTestObservations(numTestObservations), numSamples(numSamples)
    {
      sigmaSamples = new double[getNumSigmaSamples()];
      trainingSamples = new double[getNumTrainingSamples()];
      if (this->numTestObservations > 0) testSamples = new double[getNumTestSamples()];
      variableCountSamples = new double[getNumVariableCountSamples()];
    }
    
    // note how dangerous this constructor is, as it accepts pointers
    // but the destructor deletes them; set to NULL before deleting
    // if you use
    Results(std::size_t numObservations, std::size_t numPredictors,
            std::size_t numTestObservations, std::size_t numSamples,
            double* sigmaSamples, double* trainingSamples,
            double* testSamples, double* variableCountSamples) :
      sigmaSamples(sigmaSamples), trainingSamples(trainingSamples), testSamples(testSamples),
      variableCountSamples(variableCountSamples), numObservations(numObservations),
      numPredictors(numPredictors), numTestObservations(numTestObservations), numSamples(numSamples)
    {
    }
    
    ~Results() {
      delete [] sigmaSamples; sigmaSamples = NULL;
      delete [] trainingSamples; trainingSamples = NULL;
      delete [] testSamples; testSamples = NULL;
      delete [] variableCountSamples; variableCountSamples = NULL;
    }
    
    std::size_t getNumSigmaSamples() { return numSamples; }
    std::size_t getNumTrainingSamples() { return numObservations * numSamples; }
    std::size_t getNumTestSamples() { return numTestObservations * numSamples; }
    std::size_t getNumVariableCountSamples() { return numPredictors * numSamples; }
  };
} // namespace dbarts

#endif // DBARTS_RESULTS_HPP
