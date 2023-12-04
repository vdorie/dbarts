#include "config.hpp"
#include "crossvalidate.hpp"

#include <algorithm> // sort
#include <cstddef> // size_t
#include <dbarts/cstdint.hpp> // uint_least32_t
#include <cmath>   // sqrt, floor
#include <errno.h>

#include <misc/alloca.h>
#include <misc/thread.h>

#include <external/io.h>
#include <external/random.h>

#include <dbarts/bartFit.hpp>
#include <dbarts/results.hpp>

#define INVALID_INDEX static_cast<size_t>(-1)

#if __cplusplus < 201112L
#  if defined(_WIN64) || SIZEOF_SIZE_T == 8
#    define SIZE_T_SPECIFIER "%lu"
#  else
#    define SIZE_T_SPECIFIER "%u"
#  endif
#else
#  define SIZE_T_SPECIFIER "%zu"
#endif

using std::size_t;

namespace {
  using namespace dbarts;
  using namespace xval;
  
  void permuteIndexArray(ext_rng* restrict generator, size_t* restrict indices, size_t length);
  
  void allocateDataStorage(const Data& origData, Data& repData, size_t numTrainingSamples, size_t numTestSamples);
  void allocateModelStorage(const Model& origModel, Model& repModel);
  
  void freeDataStorage(Data& repData);
  void freeModelStorage(Model& repModel);

  struct CellParameters {
    size_t numTrees;
    double k;
    double power;
    double base;
    std::uint8_t status;
  };
  
  void updateFitForCell(BARTFit& fit, Control& repControl, Model& repModel, const CellParameters& parameters,
                        size_t threadId, size_t cellIndex, misc_btm_manager_t manager, bool verbose);
  
  struct SharedData {
    misc_btm_manager_t threadManager;
    Method method;
    
    const Control& control;
    const Model& model;
    const Data& data;
    
    size_t numInitialBurnIn;
    size_t numContextShiftBurnIn;
    size_t numRepBurnIn;
    
    sizetOrDouble testSampleSize;
    
    const LossFunctorDefinition& lossFunctorDef;
    
    size_t numReps;
    const CellParameters* parameters;
  };
  
  struct ThreadData {
    SharedData* shared;
    ext_rng* rng;
    
    size_t repCellOffset;
    size_t numRepCells;
    size_t threadId;
    
    double* results;
  };
  
  ext_rng* createSingleThreadedRNG(
    ext_rng_algorithm_t rng_algorithm,
    ext_rng_standardNormal_t rng_standardNormal,
    std::uint_least32_t rng_seed,
    const char*& errorMessage);
  
  ext_rng* createMultiThreadedRNG(
    ext_rng_algorithm_t rng_algorithm,
    ext_rng_standardNormal_t rng_standardNormal,
    ext_rng* seedGenerator,
    const char*& errorMessage);
  
  ext_rng* createSeedRNG(
    ext_rng_algorithm_t rng_algorithm,
    std::uint_least32_t rng_seed,
    const char*& errorMessage);
  
  bool ensureRNGSeedsAreUnique(
    const ext_rng* rng_1,
    ext_rng* rng_2,
    ext_rng* seedGenerator);
}

extern "C" { static void crossvalidationTask(void* data); }
// extern "C" { static void printInfo(void** data, size_t numThreads); }

namespace dbarts { namespace xval {
    void crossvalidate(const Control& origControl, const Model& origModel, const Data& origData,
                       Method method, sizetOrDouble testSampleSize, size_t numReps,
                       size_t numInitialBurnIn, size_t numContextShiftBurnIn, size_t numRepBurnIn,
                       const LossFunctorDefinition& lossFunctorDef, size_t numThreads,
                       const size_t* nTrees, size_t numNTrees, const double* k, size_t numKs,
                       const double* power, size_t numPowers, const double* base, size_t numBases,
                       double* results)

  {
    if (origControl.verbose) {
      ext_printf("starting %s crossvalidation with ", method == RANDOM_SUBSAMPLE ? "random subsample" : "k-fold");
      if (method == RANDOM_SUBSAMPLE)
        ext_printf("%.2f%% test obs", 100.0 * testSampleSize.p);
      else
        ext_printf(SIZE_T_SPECIFIER " folds", testSampleSize.n);
      ext_printf(", " SIZE_T_SPECIFIER " replications\n", numReps);
      ext_printf("  " SIZE_T_SPECIFIER " tree par(s), ", numNTrees);
      if (k != NULL)
        ext_printf(SIZE_T_SPECIFIER " k par(s), ", numKs);
      else 
        ext_printf("k w/hyperprior, ");
      ext_printf(" " SIZE_T_SPECIFIER " power par(s), " SIZE_T_SPECIFIER " base par(s)\n", numPowers, numBases);
      ext_printf("  results of type: %s, " SIZE_T_SPECIFIER " thread(s)\n", lossFunctorDef.displayString, numThreads);
      ext_printf("  num samp: " SIZE_T_SPECIFIER ", num reps: " SIZE_T_SPECIFIER "\n", origControl.defaultNumSamples, numReps);
      ext_printf("  burn in: " SIZE_T_SPECIFIER " first, " SIZE_T_SPECIFIER " shift, " SIZE_T_SPECIFIER " rep\n\n", numInitialBurnIn, numContextShiftBurnIn, numRepBurnIn);

      if (numThreads > 1)
        ext_printf("  parameters for [thread num, cell number]; some cells may be split across multiple threads:\n\n");
      ext_fflush_stdout();
    }
    
    size_t numRepCells = numNTrees * numKs * numPowers * numBases * numReps;
    if (numRepCells < numThreads) numThreads = numRepCells;
    
    size_t numCells = numRepCells / numReps;
    CellParameters* cellParameters = new CellParameters[numCells];
    
    size_t cellNumber = 0;
    for (size_t nIndex = 0; nIndex < numNTrees; ++nIndex) {
      for (size_t kIndex = 0; kIndex < numKs; ++kIndex) {
        for (size_t pIndex = 0; pIndex < numPowers; ++pIndex) {
          for (size_t bIndex = 0; bIndex < numBases; ++bIndex) {
            
            cellParameters[cellNumber].numTrees = nTrees[nIndex];
            cellParameters[cellNumber].k        = k != NULL ? k[kIndex] : -1.0;
            cellParameters[cellNumber].power    = power[pIndex];
            cellParameters[cellNumber].base     = base[bIndex];
            cellNumber++;
          }
        }
      }
    }
    
    Control threadControl = origControl;
    
    SharedData sharedData = { 0, method, threadControl, origModel, origData,
                              numInitialBurnIn, numContextShiftBurnIn, numRepBurnIn,
                              testSampleSize,
                              lossFunctorDef, numReps, cellParameters };
    
    ext_rng_algorithm_t rng_algorithm = static_cast<ext_rng_algorithm_t>(threadControl.rng_algorithm);
    ext_rng_standardNormal_t rng_standardNormal = static_cast<ext_rng_standardNormal_t>(threadControl.rng_standardNormal);
    std::uint_least32_t rng_seed = threadControl.rng_seed;
    threadControl.rng_algorithm = RNG_ALGORITHM_USER_POINTER;
    
    if (numThreads <= 1) {
      // Compiler complains that errorMessage may be unitialized unless I set it
      // explicitly. I have no idea how that could be the case.
      const char* errorMessage = "";
      ext_rng* rng = createSingleThreadedRNG(
        rng_algorithm, rng_standardNormal, rng_seed, errorMessage);
      
      if (rng == NULL) {
        delete [] cellParameters;
        ext_throwError(errorMessage);
      }

      ThreadData threadData = { &sharedData, rng, 0, numRepCells, 0, results };
      
      crossvalidationTask(&threadData);
      
      ext_rng_destroy(rng);
    } else {
      ext_rng* seedGenerator = NULL;
      const char* errorMessage;
      if (rng_seed != DBARTS_CONTROL_INVALID_SEED &&
          (seedGenerator = createSeedRNG(rng_algorithm, rng_seed, errorMessage)) == NULL)
      {
        delete [] cellParameters;
        ext_throwError(errorMessage);
      }
      
      misc_btm_manager_t threadManager;
      misc_btm_create(&threadManager, numThreads);
      sharedData.threadManager = threadManager;
      
      size_t numRepCellsPerThread;
      size_t offByOneIndex;
      
      misc_btm_getNumThreadsForJob(threadManager, numRepCells, 0, NULL, &numRepCellsPerThread, &offByOneIndex);

      
      ThreadData* threadData = misc_stackAllocate(numThreads, ThreadData);
      void** threadDataPtrs  = misc_stackAllocate(numThreads, void*);
      for (size_t i = 0; i < numThreads; ++i) {
        threadData[i].shared = &sharedData;
        threadData[i].rng = createMultiThreadedRNG(
          rng_algorithm,
          rng_standardNormal,
          seedGenerator,
          errorMessage);
        if (threadData[i].rng == NULL) {
          for (size_t j = i; j > 0; --j) ext_rng_destroy(threadData[j - 1].rng);
          if (seedGenerator != NULL) ext_rng_destroy(seedGenerator);
          delete [] cellParameters;
          ext_throwError(errorMessage);
        }
        if (i > 0 && !ensureRNGSeedsAreUnique(threadData[i - 1].rng, threadData[i].rng, seedGenerator)) {
          ext_rng_destroy(threadData[i].rng);
          for (size_t j = i; j > 0; --j) ext_rng_destroy(threadData[j - 1].rng);
          if (seedGenerator != NULL) ext_rng_destroy(seedGenerator);
          delete [] cellParameters;
          ext_throwError("unable to ensure seed uniqueness");
        }
        threadData[i].threadId = i;
        if (i < offByOneIndex) {
          threadData[i].numRepCells   = numRepCellsPerThread;
          threadData[i].repCellOffset = i * numRepCellsPerThread;
          threadData[i].results = results + i * numRepCellsPerThread * lossFunctorDef.numResults;
        } else {
          threadData[i].numRepCells   = numRepCellsPerThread - 1;
          threadData[i].repCellOffset = offByOneIndex * numRepCellsPerThread + (i - offByOneIndex) * (numRepCellsPerThread - 1);
          threadData[i].results = results + threadData[i].repCellOffset * lossFunctorDef.numResults;
        }
        threadDataPtrs[i] = threadData + i;
      }
      
      if (seedGenerator != NULL) ext_rng_destroy(seedGenerator);
      
      misc_btm_runTasks(threadManager, crossvalidationTask, threadDataPtrs, numThreads);
      
      for (size_t i = 0; i < numThreads; ++i)
        ext_rng_destroy(threadData[i].rng);
         
      misc_stackFree(threadDataPtrs);
      misc_stackFree(threadData);
        
      misc_btm_destroy(threadManager);
      
    }
    
    delete [] cellParameters;
    
  }
} }

namespace {
  struct CrossvalidationData {
    BARTFit& fit;
    const Data& origData;
    Data& repData;
    size_t numBurnIn;
  };
  struct ThreadScratch {
    size_t maxNumTrainingObservations;
    size_t maxNumTestObservations;
    double* y_test;
    Results* samples;
    double* weights;
    LossFunctor* lf;
    ext_rng* generator;
    size_t* permutation;
  };
  struct RandomSubsampleThreadScratch : ThreadScratch {
    
  };
  struct KFoldThreadScratch : ThreadScratch {
    size_t numFolds;
    size_t numFullSizedFolds;
    size_t numResults;
    size_t numRepBurnIn;
  };
  
  void randomSubsampleCrossvalidate(CrossvalidationData& xvalData,
                                    Results* restrict samples, size_t numSamples, double* restrict results,
                                    LossFunction calculateLoss,
                                    misc_btm_manager_t manager, size_t threadId, bool lossRequiresMutex,
                                    ThreadScratch* v_scratch);
  void kFoldCrossvalidate(CrossvalidationData& data,
                          Results* restrict samples, size_t numSamples, double* restrict results,
                          LossFunction calculateLoss,
                          misc_btm_manager_t manager, size_t threadId, bool lossRequiresMutex,
                          ThreadScratch* v_scratch);
  
  void randomSubsampleDivideData(const Data& restrict origData, Data& restrict repData,
                                 double* restrict y_test, double* restrict weights,
                                 ext_rng* restrict generator, size_t* restrict permutation);
  void kFoldDivideData(const Data& restrict origData, Data& restrict repData,
                       double* restrict y_test, double* restrict weights,
                       size_t k, size_t maxNumFoldObservations, size_t numFullSizedFolds,
                       const size_t* restrict permutation);
}

namespace {
  struct LossFunctorCreatorData {
    const LossFunctorDefinition& lfDef;
    Method method;
    size_t maxNumTestObservations;
    size_t numSamples;
    bool hasWeights;
    LossFunctor** result;
  };
  struct LossFunctorDestructorData {
    const LossFunctorDefinition& lfDef;
    LossFunctor* lf;
  };
}

extern "C" void lossFunctorCreatorTask(void* data) {
  LossFunctorCreatorData& lfcd(*static_cast<LossFunctorCreatorData*>(data));
  *lfcd.result = lfcd.lfDef.createFunctor(lfcd.lfDef, lfcd.method, lfcd.maxNumTestObservations, lfcd.numSamples, lfcd.hasWeights);
}
extern "C" void lossFunctorDestructorTask(void* data) {
  LossFunctorDestructorData& lfdd(*static_cast<LossFunctorDestructorData*>(data));
  lfdd.lfDef.deleteFunctor(lfdd.lf);
}

extern "C" {
  using namespace dbarts;
  using namespace dbarts::xval;
  
  static void crossvalidationTask(void* v_data)
  {
    ThreadData& threadData(*static_cast<ThreadData*>(v_data));
    const SharedData& sharedData(*threadData.shared);
    
    const Control& origControl(sharedData.control);
    const Model& origModel(sharedData.model);
    const Data& origData(sharedData.data);
    
    size_t maxNumTestObservations, maxNumTrainingObservations;
    if (sharedData.method == K_FOLD) {
      maxNumTestObservations = origData.numObservations / sharedData.testSampleSize.n;
      if (origData.numObservations % sharedData.testSampleSize.n == 0)
        maxNumTrainingObservations = origData.numObservations - maxNumTestObservations;
      else
        maxNumTrainingObservations = origData.numObservations - maxNumTestObservations++;
    } else {
      maxNumTestObservations =
        static_cast<size_t>(std::floor(static_cast<double>(origData.numObservations) * sharedData.testSampleSize.p + 0.5));
      maxNumTrainingObservations = origData.numObservations - maxNumTestObservations;
    }
    
    size_t numSamples = origControl.defaultNumSamples;
        
    const LossFunctorDefinition& lfDef(sharedData.lossFunctorDef);
    bool lossRequiresMutex = lfDef.requiresMutex && !misc_btm_isNull(sharedData.threadManager);
    bool hasWeights = origData.weights != NULL;
    
    LossFunctor* lf = NULL;
    if (lossRequiresMutex) {
      LossFunctorCreatorData lfcd = { lfDef, sharedData.method, maxNumTestObservations, numSamples, hasWeights, &lf };
      misc_btm_runTaskInParentThread(sharedData.threadManager, threadData.threadId, &lossFunctorCreatorTask, &lfcd);
    } else {
      lf = lfDef.createFunctor(lfDef, sharedData.method, maxNumTestObservations, numSamples, hasWeights);
    }
    
    double* suppliedY_test      = lfDef.y_testOffset      >= 0 ?
                                  *reinterpret_cast<double**>(reinterpret_cast<char*>(lf) + lfDef.y_testOffset) :
                                  NULL;
    double* suppliedTestSamples = lfDef.testSamplesOffset >= 0 ?
                                  *reinterpret_cast<double**>(reinterpret_cast<char*>(lf) + lfDef.testSamplesOffset) :
                                  NULL;
    double* suppliedWeights = NULL;
    if (hasWeights) {
      suppliedWeights = lfDef.weightsOffset >= 0 ? 
                        *reinterpret_cast<double**>(reinterpret_cast<char*>(lf) + lfDef.weightsOffset) :
                        NULL;
    }
    double* y_test = (suppliedY_test == NULL ? new double[maxNumTestObservations] : suppliedY_test);
    double* weights = NULL;
    if (hasWeights)
      weights = (suppliedWeights == NULL ? new double[maxNumTestObservations] : suppliedWeights);
    
    Results* samples =
      suppliedTestSamples == NULL ?
        new Results(maxNumTrainingObservations, origData.numPredictors, maxNumTestObservations, numSamples, 1, origModel.kPrior != NULL) :
        new Results(maxNumTrainingObservations, origData.numPredictors, maxNumTestObservations, numSamples, 1,
                    new double[numSamples],
                    new double[maxNumTrainingObservations * numSamples],
                    suppliedTestSamples,
                    new std::uint32_t[origData.numPredictors * numSamples],
                    origModel.kPrior != NULL ? new double[numSamples] : NULL);
    
    Control repControl = origControl;
    repControl.numThreads = 1;
    repControl.verbose = false;
    bool verbose = origControl.verbose;

    Data repData;
    Model repModel;
    allocateDataStorage(origData, repData, maxNumTrainingObservations, maxNumTestObservations);
    allocateModelStorage(origModel, repModel);
    
    BARTFit* fit = new BARTFit(repControl, repModel, repData);
    fit->state[0].rng = threadData.rng;
    
    CrossvalidationData xvalData = { *fit, origData, repData, 0 };
    
    void (*crossvalidate)(CrossvalidationData& data,
                          Results* restrict samples, size_t numSamples, double* restrict results,
                          LossFunction calculateLoss, misc_btm_manager_t manager, size_t threadId, bool lossRequiresMutex,
                          ThreadScratch* v_scratch);

    
    ThreadScratch* v_threadScratch;
    if (sharedData.method == K_FOLD) {
      const size_t& numFolds(sharedData.testSampleSize.n);
      size_t numFullSizedFolds = maxNumTestObservations * numFolds == origData.numObservations ? numFolds : origData.numObservations % numFolds;
      
      KFoldThreadScratch* threadScratch = new KFoldThreadScratch;
      threadScratch->numFolds = numFolds;
      threadScratch->numFullSizedFolds = numFullSizedFolds;
      threadScratch->numResults = lfDef.numResults;
      threadScratch->numRepBurnIn = sharedData.numRepBurnIn;
      
      v_threadScratch = threadScratch;

      crossvalidate = &kFoldCrossvalidate;
    } else {
      RandomSubsampleThreadScratch* threadScratch = new RandomSubsampleThreadScratch;
      v_threadScratch = threadScratch;

      crossvalidate = &randomSubsampleCrossvalidate;
    }
    v_threadScratch->maxNumTrainingObservations = maxNumTrainingObservations;
    v_threadScratch->maxNumTestObservations = maxNumTestObservations;
    v_threadScratch->y_test = y_test;
    v_threadScratch->weights = weights;
    v_threadScratch->samples = samples;
    v_threadScratch->lf = lf;
    v_threadScratch->generator = threadData.rng;
    v_threadScratch->permutation = new size_t[origData.numObservations];
    for (size_t i = 0; i < origData.numObservations; ++i) v_threadScratch->permutation[i] = i;
    
        
    // each cell is handled completely by this thread
    // cell-reps are partially handled by other threads; first is what we have to clean up from
    // someone else, last is what we partially get into ourselves
    //
    // "Rows" are "cells", that is specific sets of hyperparameters
    // "Cols" are the reps of each set of hyperpars
    size_t firstCell = threadData.repCellOffset / sharedData.numReps; // start row
    size_t firstCellRep = threadData.repCellOffset % sharedData.numReps; // start col
    size_t lastCell = (threadData.repCellOffset + threadData.numRepCells - 1) / sharedData.numReps; // end row
    size_t lastCellRep = (threadData.repCellOffset + threadData.numRepCells - 1) % sharedData.numReps + 1; // end col

    size_t resultIndex = 0;
    xvalData.numBurnIn = sharedData.numInitialBurnIn;
    
    // first and last cells are a bit of a mess, since there can be a lot of off-by-one stuff
    // process the first incomplete row
    if (firstCellRep != 0) {
      updateFitForCell(*fit, repControl, repModel, sharedData.parameters[firstCell],
                       threadData.threadId, firstCell, sharedData.threadManager, verbose);
      size_t repEnd;
      if (firstCell == lastCell) {
         // If all one row, handle it now and reset the end column
         repEnd = lastCellRep;
         lastCellRep = 0;
      } else {
         // Otherwise, handle next full row and whatever in next is necessary
         repEnd = sharedData.numReps;
         ++firstCell;
      }
      for (size_t repIndex = firstCellRep; repIndex < repEnd; ++repIndex)
      {
        crossvalidate(xvalData, samples, numSamples, threadData.results + resultIndex,
                      lfDef.calculateLoss, sharedData.threadManager, threadData.threadId, lossRequiresMutex, v_threadScratch);

        resultIndex += lfDef.numResults;
        
        xvalData.numBurnIn = sharedData.numRepBurnIn;
      }
      
      xvalData.numBurnIn = sharedData.numContextShiftBurnIn;
    }
    
    for (size_t cellIndex = firstCell; cellIndex < lastCell; ++cellIndex) {
      updateFitForCell(*fit, repControl, repModel, sharedData.parameters[cellIndex],
                       threadData.threadId, cellIndex, sharedData.threadManager, verbose);
      
      for (size_t repIndex = 0; repIndex < sharedData.numReps; ++repIndex)
      {
        crossvalidate(xvalData, samples, numSamples, threadData.results + resultIndex,
                      lfDef.calculateLoss, sharedData.threadManager, threadData.threadId, lossRequiresMutex, v_threadScratch);
        resultIndex += lfDef.numResults;
        
        xvalData.numBurnIn = sharedData.numRepBurnIn;
      }
      
      xvalData.numBurnIn = sharedData.numContextShiftBurnIn;
    }
    
    if (lastCellRep != 0) {
      updateFitForCell(*fit, repControl, repModel, sharedData.parameters[lastCell],
                       threadData.threadId, lastCell, sharedData.threadManager, verbose);
      
      for (size_t repIndex = 0; repIndex < lastCellRep; ++repIndex)
      {
        crossvalidate(xvalData, samples, numSamples, threadData.results + resultIndex,
                      lfDef.calculateLoss, sharedData.threadManager, threadData.threadId, lossRequiresMutex, v_threadScratch);
        
        resultIndex += lfDef.numResults;
        
        xvalData.numBurnIn = sharedData.numRepBurnIn;
      }
    }
    
    
    delete [] v_threadScratch->permutation;
    if (sharedData.method == K_FOLD) {
      delete reinterpret_cast<KFoldThreadScratch*>(v_threadScratch);
    } else {
      delete reinterpret_cast<RandomSubsampleThreadScratch*>(v_threadScratch);
    }
    
    if (suppliedY_test == NULL) delete [] y_test;
    if (hasWeights && suppliedWeights == NULL) delete [] weights;
    
    delete fit;
    
    freeModelStorage(repModel);
    freeDataStorage(repData);
    
    if (suppliedTestSamples != NULL) samples->testSamples = NULL;
    delete samples;
    
    if (lossRequiresMutex) {
      LossFunctorDestructorData lfdd = { lfDef, lf };
      misc_btm_runTaskInParentThread(sharedData.threadManager, threadData.threadId, &lossFunctorDestructorTask, &lfdd);
    } else {
      lfDef.deleteFunctor(lf);
    }
  }
}


namespace {
  struct LossFunctorData {
    LossFunction calculateLoss;
    LossFunctor& lf;
    const double* y_test;
    size_t numTestObservations;
    const double* testSamples;
    const double* weights;
    size_t numSamples;
    double* results;
  };
}

extern "C" void lossFunctorTask(void* data) {
  LossFunctorData& lfd(*static_cast<LossFunctorData*>(data));
  lfd.calculateLoss(lfd.lf, lfd.y_test, lfd.numTestObservations, lfd.testSamples, lfd.numSamples, lfd.weights, lfd.results);
}

namespace {
  void randomSubsampleCrossvalidate(CrossvalidationData& xvalData,
                                    Results* restrict samples, size_t numSamples, double* restrict results,
                                    LossFunction calculateLoss, misc_btm_manager_t manager, size_t threadId, bool lossRequiresMutex,
                                    ThreadScratch* v_scratch)
  {
    RandomSubsampleThreadScratch& threadScratch(*reinterpret_cast<RandomSubsampleThreadScratch *>(v_scratch));
    
    randomSubsampleDivideData(xvalData.origData, xvalData.repData, threadScratch.y_test, threadScratch.weights,
                              threadScratch.generator, threadScratch.permutation);
    xvalData.fit.setData(xvalData.repData);

    xvalData.fit.runSampler(xvalData.numBurnIn, samples);
    
    if (lossRequiresMutex) {
      LossFunctorData ldf = { calculateLoss, *threadScratch.lf, threadScratch.y_test, threadScratch.maxNumTestObservations, threadScratch.weights, samples->testSamples, numSamples, results };
      misc_btm_runTaskInParentThread(manager, threadId, &lossFunctorTask, &ldf);
    } else {
      calculateLoss(*threadScratch.lf, threadScratch.y_test, threadScratch.maxNumTestObservations,
                    samples->testSamples, numSamples, threadScratch.weights, results);
    }
  }
  
  void kFoldCrossvalidate(CrossvalidationData& xvalData,
                          Results* restrict samples, size_t numSamples, double* restrict results,
                          LossFunction calculateLoss, misc_btm_manager_t manager, size_t threadId, bool lossRequiresMutex,
                          ThreadScratch* v_scratch)
  {
    KFoldThreadScratch& threadScratch(*reinterpret_cast<KFoldThreadScratch *>(v_scratch));
    
    permuteIndexArray(threadScratch.generator, threadScratch.permutation, xvalData.origData.numObservations);
    
    for (size_t k = 0; k < threadScratch.numFolds; ++k) {
      size_t numTestObservations, foldStartIndex;
      if (k < threadScratch.numFullSizedFolds) {
         numTestObservations = threadScratch.maxNumTestObservations;
         foldStartIndex = k * threadScratch.maxNumTestObservations;
      } else {
         numTestObservations = threadScratch.maxNumTestObservations - 1;
         foldStartIndex = threadScratch.numFullSizedFolds * threadScratch.maxNumTestObservations + (k - threadScratch.numFullSizedFolds) * (threadScratch.maxNumTestObservations - 1);
      }
      
      std::sort(threadScratch.permutation + foldStartIndex, threadScratch.permutation + foldStartIndex + numTestObservations);
    }
    
    double* foldResults = misc_stackAllocate(threadScratch.numResults, double);
    
    for (size_t i = 0; i < threadScratch.numResults; ++i) results[i] = 0.0;
    
    for (size_t k = 0; k < threadScratch.numFolds; ++k) {
      size_t numTestObservations = k < threadScratch.numFullSizedFolds ? threadScratch.maxNumTestObservations : threadScratch.maxNumTestObservations - 1;
      size_t numTrainingObservations = xvalData.origData.numObservations - numTestObservations;
      
      xvalData.repData.numObservations = numTrainingObservations;
      xvalData.repData.numTestObservations = numTestObservations;
      
      kFoldDivideData(xvalData.origData, xvalData.repData, threadScratch.y_test, threadScratch.weights,
                      k, threadScratch.maxNumTestObservations, threadScratch.numFullSizedFolds, threadScratch.permutation);
      xvalData.fit.setData(xvalData.repData);
      
      xvalData.fit.runSampler(xvalData.numBurnIn, samples);
    
      if (lossRequiresMutex) {
        LossFunctorData ldf = { calculateLoss, *threadScratch.lf, threadScratch.y_test, numTestObservations, samples->testSamples, threadScratch.weights, numSamples, foldResults };
        misc_btm_runTaskInParentThread(manager, threadId, &lossFunctorTask, &ldf);
      } else {
        calculateLoss(*threadScratch.lf, threadScratch.y_test, numTestObservations, samples->testSamples, numSamples, threadScratch.weights, foldResults);
      }
      
      for (size_t i = 0; i < threadScratch.numResults; ++i) results[i] += foldResults[i];
      
      if (k > 0) xvalData.numBurnIn = threadScratch.numRepBurnIn;
    }
    
    for (size_t i = 0; i < threadScratch.numResults; ++i) results[i] /= static_cast<double>(threadScratch.numFolds);
    
    misc_stackFree(foldResults);
  }
  
  void permuteIndexArray(ext_rng* restrict generator, size_t* restrict indices, size_t length)
  {
    size_t temp, swapPos;
    for (size_t i = 0; i < length - 1; ++i) {
      swapPos = static_cast<size_t>(ext_rng_simulateUnsignedIntegerUniformInRange(generator, i, length));
      
      temp = indices[i];
      indices[i] = indices[swapPos];
      indices[swapPos] = temp;
    }
  }
  
  void randomSubsampleDivideData(const Data& restrict origData, Data& restrict repData,
                                 double* restrict y_test, double* restrict weights,
                                 ext_rng* restrict generator, size_t* restrict permutation)
  {
    double* restrict y = const_cast<double* restrict>(repData.y);
    double* restrict x = const_cast<double* restrict>(repData.x);
    double* restrict x_test = const_cast<double* restrict>(repData.x_test);
    double* restrict repWeights = const_cast<double* restrict>(repData.weights);
    
    size_t numTrainingObservations = repData.numObservations;
    size_t numTestObservations     = repData.numTestObservations;
    
    permuteIndexArray(generator, permutation, origData.numObservations);
    std::sort(permutation, permutation + numTestObservations);
    std::sort(permutation + numTestObservations, permutation + origData.numObservations);
    
    for (size_t i = 0; i < numTestObservations; ++i) {
      size_t obsIndex = *permutation++;
      y_test[i] = origData.y[obsIndex];
      for (size_t j = 0; j < origData.numPredictors; ++j) {
        x_test[i + j * numTestObservations] = origData.x[obsIndex + j * origData.numObservations];
      }
      if (weights != NULL) weights[i] = origData.weights[obsIndex];
    }
    for (size_t i = 0; i < numTrainingObservations; ++i) {
      size_t obsIndex = *permutation++;
      y[i] = origData.y[obsIndex];
      for (size_t j = 0; j < origData.numPredictors; ++j) {
        x[i + j * numTrainingObservations] = origData.x[obsIndex + j * origData.numObservations];
      }
      if (repWeights != NULL) repWeights[i] = origData.weights[obsIndex];
    }
    
    /* ext_printf("training data:\n");
    for (size_t i = 0; i < numTrainingObservations; ++i) {
      ext_printf("%.4f %.4f", repData.y[i], repData.x[i]);
      for (size_t j = 1; j < origData.numPredictors; ++j)
        ext_printf(" %.4f", repData.x[i + j * numTrainingObservations]);
      ext_printf("\n");
    }
    ext_printf("\ntest data:\n");
    for (size_t i = 0; i < numTestObservations; ++i) {
      ext_printf("%.4f %.4f", y_test[i], repData.x_test[i]);
      for (size_t j = 1; j < origData.numPredictors; ++j)
        ext_printf(" %.4f", repData.x_test[i + j * numTestObservations]);
      ext_printf("\n");
    }
    ext_printf("\n"); */
  }
  
  void kFoldDivideData(const Data& restrict origData, Data& restrict repData,
                       double* restrict y_test, double* restrict weights,
                       size_t k, size_t maxNumFoldObservations, size_t numFullSizedFolds,
                       const size_t* restrict permutation)
  {
    size_t i, j, obsIndex, foldStartIndex;
    double* restrict y = const_cast<double* restrict>(repData.y);
    double* restrict x = const_cast<double* restrict>(repData.x);
    double* restrict x_test = const_cast<double* restrict>(repData.x_test);
    double* restrict repWeights = const_cast<double* restrict>(repData.weights);
    
    size_t numTrainingObservations = repData.numObservations;
    size_t numTestObservations     = repData.numTestObservations;
    
    if (k < numFullSizedFolds) {
      foldStartIndex = k * maxNumFoldObservations;
    } else {
      foldStartIndex = numFullSizedFolds * maxNumFoldObservations + (k - numFullSizedFolds) * (maxNumFoldObservations - 1);
    }
    
    // i is always the target observation number, so adjust the source index
    for (i = 0; i < numTestObservations; ++i) {
      obsIndex = permutation[i + foldStartIndex];
      y_test[i] = origData.y[obsIndex];
      for (j = 0; j < origData.numPredictors; ++j)
        x_test[i + j * numTestObservations] = origData.x[obsIndex + j * origData.numObservations];
      if (weights != NULL) weights[i] = origData.weights[obsIndex];

    }
    for (i = 0; i < foldStartIndex; ++i) {
      obsIndex = permutation[i];
      y[i] = origData.y[obsIndex];
      for (j = 0; j < origData.numPredictors; ++j) {
        x[i + j * numTrainingObservations] = origData.x[obsIndex + j * origData.numObservations];
      }
      if (repWeights != NULL) repWeights[i] = origData.weights[obsIndex];
    }
    for ( /* */; i < numTrainingObservations; ++i) {
      obsIndex = permutation[i + numTestObservations];
      y[i] = origData.y[obsIndex];
      for (j = 0; j < origData.numPredictors; ++j) {
        x[i + j * numTrainingObservations] = origData.x[obsIndex + j * origData.numObservations];
      }
      if (repWeights != NULL) repWeights[i] = origData.weights[obsIndex];
    }
  }
  
  void allocateDataStorage(const Data& origData, Data& repData, size_t numTrainingObservations, size_t numTestObservations)
  {
    repData.y = new double[numTrainingObservations];
    repData.x = new double[numTrainingObservations * origData.numPredictors];
    repData.x_test = new double[numTestObservations * origData.numPredictors];
    
    // copy in some default values for use with the BART constructor
    // these could be skipped by creating a constructor that waits for
    // data to be assigned or through tricky use of placement new
    std::memcpy(const_cast<double*>(repData.y), origData.y, numTrainingObservations * sizeof(double));
    std::memcpy(const_cast<double*>(repData.x), origData.x, numTrainingObservations * origData.numPredictors * sizeof(double));
    std::memcpy(const_cast<double*>(repData.x_test), origData.x, numTestObservations * origData.numPredictors * sizeof(double));
     
    // should be OK to leave weights and test offset garbage until run() is called
    repData.weights = origData.weights != NULL ? new double[numTrainingObservations]  : NULL;
    
    if (origData.offset == NULL) {
      repData.offset = NULL;
      repData.testOffset = NULL;
    } else {
      repData.offset     = new double[numTrainingObservations];
      repData.testOffset = new double[numTestObservations];
     
      std::memcpy(const_cast<double*>(repData.offset), origData.offset, numTrainingObservations * sizeof(double));
    }
    
    repData.numObservations     = numTrainingObservations;
    repData.numPredictors       = origData.numPredictors;
    repData.numTestObservations = numTestObservations;
    repData.sigmaEstimate       = origData.sigmaEstimate;
    
    repData.variableTypes = origData.variableTypes;
    repData.maxNumCuts    = origData.maxNumCuts;
  }
  
  void allocateModelStorage(const Model& origModel, Model& repModel)
  {
    repModel.birthOrDeathProbability = origModel.birthOrDeathProbability;
    repModel.swapProbability = origModel.swapProbability;
    repModel.changeProbability = origModel.changeProbability;
    repModel.birthProbability = origModel.birthProbability;
    repModel.nodeScale = origModel.nodeScale;
    
    CGMPrior* repTreePrior = new CGMPrior();
    const CGMPrior* oldTreePrior = static_cast<CGMPrior*>(origModel.treePrior);
    repTreePrior->base = oldTreePrior->base;
    repTreePrior->power = oldTreePrior->power;
    repTreePrior->splitProbabilities = oldTreePrior->splitProbabilities;
    
    repModel.treePrior = repTreePrior;
    
    
    NormalPrior* repNodePrior = new NormalPrior();
    const NormalPrior* oldNodePrior = static_cast<NormalPrior*>(origModel.muPrior);
    repNodePrior->scale = oldNodePrior->scale;
    
    repModel.muPrior = repNodePrior;
    
    repModel.sigmaSqPrior = origModel.sigmaSqPrior->duplicate();
    
    if (!origModel.kPrior->isFixed) {
      const ChiHyperprior* oldKPrior = static_cast<ChiHyperprior*>(origModel.kPrior);
      repModel.kPrior = new ChiHyperprior(oldKPrior->degreesOfFreedom, oldKPrior->scale);
    } else {
      const FixedHyperprior* oldKPrior = static_cast<FixedHyperprior*>(origModel.kPrior);
      repModel.kPrior = new FixedHyperprior(oldKPrior->getK());
    }
  }
  
  void freeDataStorage(Data& repData)
  {
    delete [] repData.testOffset; repData.testOffset = NULL;
    delete [] repData.offset;     repData.offset     = NULL;
    delete [] repData.weights;    repData.weights    = NULL;
    
    delete [] repData.x_test; repData.x_test = NULL;
    delete [] repData.x;      repData.x      = NULL;
    delete [] repData.y;      repData.y      = NULL;
  }
  
  void freeModelStorage(Model& repModel)
  {
    delete repModel.kPrior;
    delete repModel.sigmaSqPrior;
    delete repModel.muPrior;
    delete repModel.treePrior;
  }
}

/* extern "C" {
  void printInfo(void** vt_data, size_t numThreads)
  {
    SharedData& sharedData(*static_cast<ThreadData*>(vt_data[0])->shared);
    
    for (size_t i = 0; i < numThreads; ++i) {
      const ThreadData& threadData(*static_cast<const ThreadData*>(vt_data[i]));
      
      size_t fittingCell = threadData.fittingCell;
      size_t printedCell = sharedData.printedCells[i];
      
      
      if (fittingCell == INVALID_INDEX || (printedCell != INVALID_INDEX && printedCell >= fittingCell)) continue;
      
      size_t firstCell = threadData.repCellOffset / sharedData.numReps;
      size_t lastCell  = (threadData.repCellOffset + threadData.numRepCells) / sharedData.numReps;
      size_t lastCellRep  = (threadData.repCellOffset + threadData.numRepCells) % sharedData.numReps;
      if (lastCellRep == 0) --lastCell;
      
      size_t firstPrintCell = printedCell == INVALID_INDEX ? firstCell : printedCell + 1;
      size_t lastPrintCell = fittingCell < lastCell ? fittingCell : lastCell;
      
      if (firstPrintCell <= lastPrintCell) {
        for (size_t j = firstPrintCell; j <= lastPrintCell; ++j) {
          ext_printf("    [" SIZE_T_SPECIFIER ", " SIZE_T_SPECIFIER "] n.tree: " SIZE_T_SPECIFIER ", k: %f, power: %f, base: %f\n",
                     i + 1, j + 1,
                     sharedData.parameters[j].numTrees, sharedData.parameters[j].k,
                     sharedData.parameters[j].power, sharedData.parameters[j].base);
        }
        sharedData.printedCells[i] = fittingCell;
      }
    }
  }
} */

namespace {
  struct PrintData {
    size_t threadId;
    size_t cellIndex;
    size_t numTrees;
    double k;
    double power;
    double base;
  };
}
extern "C" void printTask(void* v_data) {
  PrintData& data(*static_cast<PrintData*>(v_data));
  ext_printf("    [" SIZE_T_SPECIFIER ", " SIZE_T_SPECIFIER "] n.trees: " SIZE_T_SPECIFIER ", ", data.threadId + 1, data.cellIndex + 1,
             data.numTrees);
  if (data.k > 0.0) ext_printf("k: %.2f, ", data.k);
  ext_printf("power: %.2f, base: %.2f\n", data.power, data.base);
}

namespace {
  void updateFitForCell(BARTFit& fit, Control& repControl, Model& repModel, const CellParameters& parameters,
                        size_t threadId, size_t cellIndex, misc_btm_manager_t manager, bool verbose)
  {
    size_t numTrees = parameters.numTrees;
    double k        = parameters.k;
    double power    = parameters.power;
    double base     = parameters.base;
      
    if (verbose) {
      if (misc_btm_isNull(manager)) {
        ext_printf("    [" SIZE_T_SPECIFIER "] n.trees: " SIZE_T_SPECIFIER ", ", cellIndex, numTrees);
        if (k > 0.0) ext_printf("k: %.2f, ", k);
        ext_printf("power: %.2f, base: %.2f\n", power, base);
      } else {
        PrintData printData = { threadId, cellIndex, numTrees, k, power, base };
        misc_btm_runTaskInParentThread(manager, threadId, &printTask, &printData);
      }
    }
    
    repControl.numTrees = numTrees;
    
    if (k > 0.0 && repModel.kPrior->isFixed)
      static_cast<FixedHyperprior*>(repModel.kPrior)->setK(k);
    
    static_cast<CGMPrior*>(repModel.treePrior)->power = power;
    static_cast<CGMPrior*>(repModel.treePrior)->base  = base;
    
    fit.setControl(repControl);
    fit.setModel(repModel);
  }
}

namespace {
  ext_rng* createSingleThreadedRNG(
    ext_rng_algorithm_t rng_algorithm,
    ext_rng_standardNormal_t rng_standardNormal,
    std::uint_least32_t rng_seed,
    const char*& errorMessage)
  {
    ext_rng* result;
    if (rng_algorithm == EXT_RNG_ALGORITHM_INVALID &&
        rng_standardNormal == EXT_RNG_STANDARD_NORMAL_INVALID)
    {
      if ((result = ext_rng_createDefault(true)) == NULL) {
        // Use built-in
        errorMessage = "could not allocate rng";
        return NULL;
      }
      if (rng_seed != DBARTS_CONTROL_INVALID_SEED && ext_rng_setSeed(result, rng_seed) != 0) {
        // Set seed if specified
        errorMessage = "could not seed rng";
        ext_rng_destroy(result);
        return NULL;
      }
      return result;
    }
    // Use custom
    if (rng_algorithm == EXT_RNG_ALGORITHM_INVALID)
      rng_algorithm = ext_rng_getDefaultAlgorithmType();
    if (rng_standardNormal == EXT_RNG_STANDARD_NORMAL_INVALID)
      rng_standardNormal = ext_rng_getDefaultStandardNormalType();
    
    if (rng_algorithm == EXT_RNG_ALGORITHM_USER_UNIFORM) {
      errorMessage = "cannot create rng with user-specified uniform distribution function";
      return NULL;
    }
    
    if (rng_standardNormal == EXT_RNG_STANDARD_NORMAL_USER_NORM) {
      errorMessage = "cannot create rng with user-specified normal distribution function";
      return NULL;
    }
    
    if (rng_seed == DBARTS_CONTROL_INVALID_SEED) {
      // No need to use a custom seed
      int errorCode = ext_rng_createAndSeed(&result, rng_algorithm, rng_standardNormal);
      if (errorCode == ENOMEM) {
        errorMessage = "could not allocate rng";
        return NULL;
      }
      if (errorCode == EINVAL) {
        errorMessage = "failure to set standard normal type or possibly seed";
        return NULL;
      }
      if (errorCode != 0) {
        errorMessage = "unknown error trying to create rng";
        return NULL;
      }
    } else {
      // Create first, then set seed
      if ((result = ext_rng_create(rng_algorithm, NULL)) == NULL) {
        errorMessage = "could not allocate rng";
        return NULL;
      }
      
      if (ext_rng_setStandardNormalAlgorithm(result, rng_standardNormal, NULL) != 0) {
        errorMessage = "could not set rng standard normal";
        ext_rng_destroy(result);
        return NULL;
      }
      if (ext_rng_setSeed(result, rng_seed) != 0) {
        errorMessage = "could not seed rng";
        ext_rng_destroy(result);
        return NULL;
      }
    }
    
    return result;
  }

  ext_rng* createMultiThreadedRNG(
    ext_rng_algorithm_t rng_algorithm,
    ext_rng_standardNormal_t rng_standardNormal,
    ext_rng* seedGenerator,
    const char*& errorMessage)
  {
    ext_rng* result = NULL;
    if (rng_algorithm == EXT_RNG_ALGORITHM_INVALID &&
        rng_standardNormal == EXT_RNG_STANDARD_NORMAL_INVALID) {
      if ((result = ext_rng_createDefault(false)) == NULL) {
        errorMessage = "could not allocate rng";
        return NULL;
      }
    } else {
      if (rng_algorithm == EXT_RNG_ALGORITHM_INVALID)
        rng_algorithm = ext_rng_getDefaultAlgorithmType();
      if (rng_standardNormal == EXT_RNG_STANDARD_NORMAL_INVALID)
        rng_standardNormal = ext_rng_getDefaultStandardNormalType();
      
      if (rng_algorithm == EXT_RNG_ALGORITHM_USER_UNIFORM) {
        errorMessage = "cannot create rng with user-specified uniform distribution function";
        return NULL;
      }

      if (rng_standardNormal == EXT_RNG_STANDARD_NORMAL_USER_NORM) {
        errorMessage = "cannot create rng with user-specified normal distribution function";
        return NULL;
      }
      
      if ((result = ext_rng_create(rng_algorithm, NULL)) == NULL) {
        errorMessage = "could not allocate rng";
        return NULL;
      }

      if (ext_rng_setStandardNormalAlgorithm(result, rng_standardNormal, NULL) != 0) {
        errorMessage = "could not set rng standard normal";
        ext_rng_destroy(result);
        return NULL;
      }
    }
    
    if (seedGenerator != NULL) {
      std::uint_least32_t seed = static_cast<std::uint_least32_t>(ext_rng_simulateUnsignedIntegerUniformInRange(seedGenerator, 0, static_cast<std::uint_least32_t>(-1)));
      if (ext_rng_setSeed(result, seed) != 0) {
        errorMessage = "could not seed rng";
        ext_rng_destroy(result);
        return NULL;
      }
    } else {
      if (ext_rng_setSeedFromClock(result) != 0) {
        errorMessage = "could not seed rng";
        ext_rng_destroy(result);
        return NULL;
      }
    }

    return result;
  }
  
  ext_rng* createSeedRNG(
    ext_rng_algorithm_t rng_algorithm,
    std::uint_least32_t rng_seed,
    const char*& errorMessage)
  {
    ext_rng* rng = NULL;
    if (rng_algorithm == EXT_RNG_ALGORITHM_INVALID) {
      rng = ext_rng_createDefault(false);
    } else {
      if (rng_algorithm == EXT_RNG_ALGORITHM_USER_UNIFORM) {
        errorMessage = "cannot create rng with user-specified uniform distribution function";
        return NULL;
      }
      
      rng = ext_rng_create(rng_algorithm, NULL);
    }
    if (rng == NULL) {
      errorMessage = "could not allocate rng";
      return NULL;
    }
    if (ext_rng_setSeed(rng, rng_seed) != 0) {
      errorMessage = "could not seed rng";
      ext_rng_destroy(rng);
      return NULL;
    }
    return rng;
  }

  bool ensureRNGSeedsAreUnique(const ext_rng* rng_1, ext_rng* rng_2, ext_rng* seedGenerator)
  {
    size_t numSeedResets;
    for (numSeedResets = 0; numSeedResets < static_cast<size_t>(-1); ++numSeedResets) {
      if (!ext_rng_seedsAreEqual(rng_1, rng_2)) return true;
      if (seedGenerator != NULL) {
        if (ext_rng_setSeed(rng_2, static_cast<std::uint_least32_t>(ext_rng_simulateUnsignedIntegerUniformInRange(seedGenerator, 0, static_cast<std::uint_least32_t>(-1)))) != 0)
          return false;
      } else {
        if (ext_rng_setSeedFromClock(rng_2) != 0)
          return false;
      }
    }
    
    if (!ext_rng_seedsAreEqual(rng_1, rng_2)) return true;

    return false;
  }
}
