#ifndef BART_CONTROL_HPP
#define BART_CONTROL_HPP

#include <cstddef> // size_t
#include "cstdint" // ints types

namespace bart {
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
    double binaryOffset;

    double sigmaQuantile;
    std::uint32_t sigmaDf;

    double kFactor;
    double base, power;
    
    std::size_t numSamples;
    std::size_t numTrees;
    std::size_t numBurnIn;
    std::size_t numThreads;
    std::uint32_t treeThinningRate;
    std::uint32_t printEvery;
    std::uint32_t printCutoffs;
    
    CallbackFunction callback;
    void* callbackData;
    
    // I think these defaults are ridiculous, but they're what BART in R uses
    Control() :
      responseIsBinary(false), verbose(true), keepTrainingFits(true), useQuantiles(false),
      binaryOffset(0.0), sigmaQuantile(0.9), sigmaDf(3), kFactor(2.0), base(0.95), power(2.0),
      numSamples(1000), numTrees(200), numBurnIn(100), numThreads(1), treeThinningRate(1), printEvery(100),
      printCutoffs(0), callback(NULL), callbackData(NULL)
    { }
    Control(std::size_t numSamples, std::size_t numTrees, std::size_t numBurnIn, std::size_t numThreads, std::uint32_t treeThinningRate, bool keepTrainingFits,
            bool verbose, uint32_t printEvery,
            bool responseIsBinary, double binaryOffset,
            uint32_t sigmaDf, double sigmaQuantile, double kFactor, double base, double power,
            bool useQuantiles, uint32_t printCutoffs,
            CallbackFunction callback, void* callbackData) :
      responseIsBinary(responseIsBinary), verbose(verbose), keepTrainingFits(keepTrainingFits), useQuantiles(useQuantiles),
      binaryOffset(binaryOffset), sigmaQuantile(sigmaQuantile), sigmaDf(sigmaDf), kFactor(kFactor), base(base), power(power),
      numSamples(numSamples), numTrees(numTrees), numBurnIn(numBurnIn), numThreads(numThreads), treeThinningRate(treeThinningRate), printEvery(printEvery),
      printCutoffs(printCutoffs), callback(callback), callbackData(callbackData)
    { }
  };
} // namespace bart

#endif // BART_CONTROL_HPP
