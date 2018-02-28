#ifndef DBARTS_XVAL_CROSSVALIDATE_HPP
#define DBARTS_XVAL_CROSSVALIDATE_HPP

#include <cstddef>

#include <dbarts/control.hpp>
#include <dbarts/data.hpp>
#include <dbarts/model.hpp>

namespace dbarts {
  namespace xval {
    // overload any structs to provide needed members
    struct LossFunctor {
    };
    
    typedef void(*LossFunction)(LossFunctor& restrict instance,
                                const double* restrict y_test, std::size_t numTestObservations,
                                const double* restrict testSamples, std::size_t numSamples, // numTestObservations x numSamples
                                double* restrict results);
    
    struct LossFunctorDefinition {
      std::ptrdiff_t y_testOffset;         // offset into functor that points to y_test, or negative if not supplied
      std::ptrdiff_t testSamplesOffset;    // offset into functor that points to testSamples, or negative if not supplied
      std::size_t numResults;              // how much to increment results pointer after calculation
      const char* displayString;
      
      // member functions
      LossFunction calculateLoss;
      LossFunctor* (*createFunctor)(const LossFunctorDefinition& def, std::size_t numTestObservations, std::size_t numSamples);
      void (*deleteFunctor)(LossFunctor* instance);
      
      virtual ~LossFunctorDefinition() { }
    };
    
    
    void crossvalidate(const Control& control, const Model& model, const Data& data,
                       double testSampleProp, std::size_t numReps,
                       std::size_t numInitialBurnIn, std::size_t numContextShiftBurnIn, std::size_t numRepBurnIn,
                       const LossFunctorDefinition& lossFunctorDef, std::size_t numThreads,
                       const std::size_t* nTrees, std::size_t numNTrees, const double* k, std::size_t numKs,
                       const double* power, std::size_t numPowers, const double* base, std::size_t numBases,
                       double* results);
  }                       
}

#endif // DBARTS_XVAL_CROSSVALIDATE_HPP

