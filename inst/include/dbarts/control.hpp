#ifndef DBARTS_CONTROL_HPP
#define DBARTS_CONTROL_HPP

#include <cstddef> // size_t
#include <dbarts/cstdint.hpp> // int types
#include <limits> // numeric_limits

#include <dbarts/random.hpp>

#define DBARTS_CONTROL_INVALID_SEED static_cast<std::uint_least32_t>(std::numeric_limits<int>::min())

namespace dbarts {
  struct BARTFit;
  
  typedef void (*CallbackFunction)(
    void* data,
    BARTFit& fit,
    bool isBurningIn,
    const double* trainingDraw,
    const double* testDraw,
    double sigma,
    const double* k
  );
  
  struct Control {
    bool responseIsBinary;
    bool verbose;
    bool keepTrainingFits;
    bool useQuantiles;
    bool keepTrees;
    
    std::size_t defaultNumSamples;
    std::size_t defaultNumBurnIn;
    std::size_t numTrees;
    std::size_t numChains;
    std::size_t numThreads;
    std::uint32_t treeThinningRate;
    std::uint32_t printEvery;
    std::uint32_t printCutoffs;
    
    rng_algorithm_t rng_algorithm;
    rng_standardNormal_t rng_standardNormal;
    std::uint_least32_t rng_seed;
    
    CallbackFunction callback;
    void* callbackData;
    
    Control() :
      responseIsBinary(false),
      verbose(true),
      keepTrainingFits(true),
      useQuantiles(false),
      keepTrees(false),
      defaultNumSamples(800),
      defaultNumBurnIn(200),
      numTrees(75),
      numChains(1),
      numThreads(1),
      treeThinningRate(1),
      printEvery(100),
      printCutoffs(0),
      rng_algorithm(RNG_ALGORITHM_MERSENNE_TWISTER),
      rng_standardNormal(RNG_STANDARD_NORMAL_INVERSION),
      rng_seed(DBARTS_CONTROL_INVALID_SEED),
      callback(NULL),
      callbackData(NULL)
    { }

    Control(
      std::size_t defaultNumSamples,
      std::size_t defaultNumBurnIn,
      std::size_t numTrees,
      std::size_t numChains,
      std::size_t numThreads,
      std::uint32_t treeThinningRate,
      bool keepTrainingFits,
      bool verbose,
      std::uint32_t printEvery,
      bool responseIsBinary,
      bool useQuantiles,
      bool keepTrees,
      std::uint32_t printCutoffs,
      rng_algorithm_t rng_algorithm,
      rng_standardNormal_t rng_standardNormal,
      std::uint_least32_t rng_seed,
      CallbackFunction callback,
      void* callbackData
    ) :
      responseIsBinary(responseIsBinary),
      verbose(verbose),
      keepTrainingFits(keepTrainingFits),
      useQuantiles(useQuantiles),
      keepTrees(keepTrees),
      defaultNumSamples(defaultNumSamples),
      defaultNumBurnIn(defaultNumBurnIn),
      numTrees(numTrees),
      numChains(numChains),
      numThreads(numThreads),
      treeThinningRate(treeThinningRate),
      printEvery(printEvery),
      printCutoffs(printCutoffs),
      rng_algorithm(rng_algorithm),
      rng_standardNormal(rng_standardNormal),
      rng_seed(rng_seed),
      callback(callback),
      callbackData(callbackData)
    { }
  };
} // namespace dbarts

#endif // BART_CONTROL_HPP

