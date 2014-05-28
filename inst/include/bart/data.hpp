#ifndef BART_DATA_HPP
#define BART_DATA_HPP

#include <cstddef> // size_t
#include <stdint.h>

namespace bart {
  struct Data {
    const double* y;
    const double* X;
    const double* X_test;
    
    std::size_t numObservations;
    std::size_t numPredictors;
    std::size_t numTestObservations;
    double sigmaEstimate;
        
    const std::uint32_t* maxNumCuts; // length = numPredictors if control.useQuantiles is true
    
    Data() :
      y(NULL), X(NULL), X_test(NULL),
      numObservations(0), numPredictors(0), numTestObservations(0),
      sigmaEstimate(1.0), maxNumCuts(NULL)
    { }
    
    Data(const double* y, const double* X, const double* X_test,
         std::size_t numObservations, std::size_t numPredictors, std::size_t numTestObservations,
         double sigmaEstimate, const std::uint32_t* maxNumCuts) :
      y(y), X(X), X_test(X_test),
      numObservations(numObservations), numPredictors(numPredictors), numTestObservations(numTestObservations),
      sigmaEstimate(sigmaEstimate), maxNumCuts(maxNumCuts)
    {
    }
  };
} // namespace bart

#endif // BART_DATA_HPP
