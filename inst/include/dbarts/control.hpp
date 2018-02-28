#ifndef DBARTS_CONTROL_HPP
#define DBARTS_CONTROL_HPP

#include <cstddef> // size_t
#include "cstdint.hpp" // int types

#include <external/random.h>

namespace dbarts {
  struct BARTFit;
  
  typedef void (*CallbackFunction)(void* data, BARTFit& fit, bool isBurningIn,
                                   const double* trainingDraw,
                                   const double* testDraw,
                                   double sigma);
  
  enum RunMode {
    SEQUENTIAL_UPDATES, FIXED_SAMPLES
  };
  
  struct Control {
    bool responseIsBinary;
    bool verbose;
    bool keepTrainingFits;
    bool useQuantiles;
    RunMode runMode;
    
    std::size_t defaultNumSamples;
    std::size_t defaultNumBurnIn;
    std::size_t numTrees;
    std::size_t numChains;
    std::size_t numThreads;
    std::uint32_t treeThinningRate;
    std::uint32_t printEvery;
    std::uint32_t printCutoffs;
    
    // these should be from external/random.h with the exception that we catch "INVALID" codes
    // and use them to construct a default RNG (e.g., one that matches the environment's)
    ext_rng_algorithm_t rng_algorithm;
    ext_rng_standardNormal_t rng_standardNormal;
    
    CallbackFunction callback;
    void* callbackData;
    
    Control() :
      responseIsBinary(false), verbose(true), keepTrainingFits(true), useQuantiles(false), runMode(SEQUENTIAL_UPDATES),
      defaultNumSamples(800), defaultNumBurnIn(200), numTrees(75), numChains(1), numThreads(1), treeThinningRate(1),
      printEvery(100), printCutoffs(0), rng_algorithm(EXT_RNG_ALGORITHM_MERSENNE_TWISTER),
      rng_standardNormal(EXT_RNG_STANDARD_NORMAL_INVERSION), callback(NULL), callbackData(NULL)
    { }
    Control(std::size_t defaultNumSamples,
            std::size_t defaultNumBurnIn,
            std::size_t numTrees,
            std::size_t numChains,
            std::size_t numThreads,
            std::uint32_t treeThinningRate,
            bool keepTrainingFits,
            bool verbose,
            uint32_t printEvery,
            bool responseIsBinary,
            bool useQuantiles,
            RunMode runMode,
            uint32_t printCutoffs,
            ext_rng_algorithm_t rng_algorithm,
            ext_rng_standardNormal_t rng_standardNormal,
            CallbackFunction callback,
            void* callbackData) :
      responseIsBinary(responseIsBinary), verbose(verbose), keepTrainingFits(keepTrainingFits), useQuantiles(useQuantiles),
      runMode(runMode), defaultNumSamples(defaultNumSamples), defaultNumBurnIn(defaultNumBurnIn), numTrees(numTrees),
      numChains(numChains), numThreads(numThreads), treeThinningRate(treeThinningRate), printEvery(printEvery),
      printCutoffs(printCutoffs), rng_algorithm(rng_algorithm), rng_standardNormal(rng_standardNormal),
      callback(callback), callbackData(callbackData)
    { }
  };
} // namespace dbarts

#endif // BART_CONTROL_HPP
