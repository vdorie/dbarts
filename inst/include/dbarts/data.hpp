#ifndef DBARTS_DATA_HPP
#define DBARTS_DATA_HPP

#include <cstddef> // size_t
#include "cstdint.hpp"

#include "types.hpp"

namespace dbarts {
  struct Data {
    const double* y;
    const double* X;
    const double* X_test;
    
    const double* weights;
    const double* offset;
    const double* testOffset;
    
    std::size_t numObservations;
    std::size_t numPredictors;
    std::size_t numTestObservations;
    double sigmaEstimate;
    
    const VariableType* variableTypes;
    const std::uint32_t* maxNumCuts; // length = numPredictors if control.useQuantiles is true
    
    Data() :
      y(NULL), X(NULL), X_test(NULL), weights(NULL), offset(NULL), testOffset(NULL),
      numObservations(0), numPredictors(0), numTestObservations(0),
      sigmaEstimate(1.0), variableTypes(NULL), maxNumCuts(NULL)
    { }
    
    Data(const double* y,
         const double* X,
         const double* X_test,
         const double* weights,
         const double* offset,
         const double* testOffset,
         std::size_t numObservations,
         std::size_t numPredictors,
         std::size_t numTestObservations,
         double sigmaEstimate,
         const VariableType* variableTypes,
         const std::uint32_t* maxNumCuts) :
      y(y), X(X), X_test(X_test), weights(weights), offset(offset), testOffset(testOffset),
      numObservations(numObservations), numPredictors(numPredictors), numTestObservations(numTestObservations),
      sigmaEstimate(sigmaEstimate), variableTypes(variableTypes), maxNumCuts(maxNumCuts)
    {
    }
  };
} // namespace dbarts

#endif // DBARTS_DATA_HPP
