#ifndef DBARTS_CONTROL_HPP
#define DBARTS_CONTROL_HPP

#include <cstddef> // size_t
#include "cstdint.hpp" // int types

namespace dbarts {
  struct BARTFit;
  
  typedef void (*CallbackFunction)(void* data, BARTFit& fit, bool isBurningIn,
                                   const double* trainingDraw,
                                   const double* testDraw,
                                   double sigma);
  
  struct Control {
    bool responseIsBinary;
    bool verbose;
    bool keepTrainingFits;
    bool useQuantiles;
    
    std::size_t numSamples;
    std::size_t numBurnIn;
    std::size_t numTrees;
    std::size_t numThreads;
    std::uint32_t treeThinningRate;
    std::uint32_t printEvery;
    std::uint32_t printCutoffs;
    
    CallbackFunction callback;
    void* callbackData;
    
    // I think these defaults are ridiculous, but they're what BART in R uses
    Control() :
      responseIsBinary(false), verbose(true), keepTrainingFits(true), useQuantiles(false),
      numSamples(1000), numBurnIn(100), numTrees(200), numThreads(1), treeThinningRate(1), printEvery(100),
      printCutoffs(0), callback(NULL), callbackData(NULL)
    { }
    Control(std::size_t numSamples,
            std::size_t numBurnIn,
            std::size_t numTrees,
            std::size_t numThreads,
            std::uint32_t treeThinningRate,
            bool keepTrainingFits,
            bool verbose,
            uint32_t printEvery,
            bool responseIsBinary,
            bool useQuantiles,
            uint32_t printCutoffs,
            CallbackFunction callback,
            void* callbackData) :
      responseIsBinary(responseIsBinary), verbose(verbose), keepTrainingFits(keepTrainingFits), useQuantiles(useQuantiles),
      numSamples(numSamples), numBurnIn(numBurnIn), numTrees(numTrees), numThreads(numThreads), treeThinningRate(treeThinningRate), printEvery(printEvery),
      printCutoffs(printCutoffs), callback(callback), callbackData(callbackData)
    { }
  };
} // namespace dbarts

#endif // BART_CONTROL_HPP
