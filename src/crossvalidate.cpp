#include "config.hpp"
#include "crossvalidate.hpp"

#include <algorithm> // sort
#include <cstddef> // size_t
#include <cmath>   // sqrt, floor

#include <misc/alloca.h>
#include <misc/thread.h>

#include <external/io.h>
#include <external/random.h>

#include <dbarts/bartFit.hpp>
#include <dbarts/results.hpp>

#define INVALID_INDEX static_cast<size_t>(-1)

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
    uint8_t status;
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
}

extern "C" { static void crossvalidationTask(void* data); }
// extern "C" { static void printInfo(void** data, size_t numThreads); }

namespace dbarts { namespace xval {
    void crossvalidate(const Control& origControl, const Model& origModel, const Data& origData,
                       Method method, sizetOrDouble testSampleSize, size_t numReps,
                       size_t numInitialBurnIn, size_t numContextShiftBurnIn, size_t numRepBurnIn,
                       const LossFunctorDefinition& lossFunctorDef, size_t numThreads,
                       const std::size_t* nTrees, size_t numNTrees, const double* k, size_t numKs,
                       const double* power, size_t numPowers, const double* base, size_t numBases,
                       double* results)

  {
    if (origControl.verbose) {
      ext_printf("starting %s crossvalidation with ", method == RANDOM_SUBSAMPLE ? "random subsample" : "k-fold");
      if (method == RANDOM_SUBSAMPLE)
        ext_printf("%.2f%% test obs", 100.0 * testSampleSize.p);
      else
        ext_printf("%lu folds", testSampleSize.n);
      ext_printf(", %lu replications\n", numReps);
      ext_printf("  %lu tree par(s), %lu k par(s), %lu power par(s), %lu base par(s)\n",
                 numNTrees, numKs, numPowers, numBases);
      ext_printf("  results of type: %s\n", lossFunctorDef.displayString);
      ext_printf("  num samp: %lu, num reps: %lu\n", origControl.defaultNumSamples, numReps);
      ext_printf("  burn in: %lu first, %lu shift, %lu rep\n\n", numInitialBurnIn, numContextShiftBurnIn, numRepBurnIn);

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
            cellParameters[cellNumber].k        = k[kIndex];
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
    
    if (numThreads <= 1) {
      ext_rng* rng;
      
      // both BART sampler and xval algorithm can use native sampler
      if (rng_algorithm == EXT_RNG_ALGORITHM_INVALID &&
          rng_standardNormal == EXT_RNG_STANDARD_NORMAL_INVALID) {
        
        if ((rng = ext_rng_createDefault(true)) == NULL)
          ext_throwError("could not allocate rng");
      
      } else {
        if (ext_rng_createAndSeed(&rng, rng_algorithm, rng_standardNormal) != 0)
          ext_throwError("could not allocate rng");
      }
      
      ThreadData threadData = { &sharedData, rng, 0, numRepCells, 0, results };
      
      crossvalidationTask(&threadData);
      
      ext_rng_destroy(rng);
    } else {
      if (rng_algorithm == EXT_RNG_ALGORITHM_INVALID)
        rng_algorithm = ext_rng_getDefaultAlgorithmType();
      if (rng_standardNormal == EXT_RNG_STANDARD_NORMAL_INVALID)
        rng_standardNormal = ext_rng_getDefaultStandardNormalType();
      
      misc_btm_manager_t threadManager;
      misc_btm_create(&threadManager, numThreads);
      sharedData.threadManager = threadManager;
      
      size_t numRepCellsPerThread;
      size_t offByOneIndex;
      
      misc_btm_getNumThreadsForJob(threadManager, numRepCells, 0, NULL, &numRepCellsPerThread, &offByOneIndex);
      
      
      ThreadData* threadData = misc_stackAllocate(numThreads, ThreadData);
      void** threadDataPtrs  = misc_stackAllocate(numThreads, void*);
      for (size_t i = 0; i < offByOneIndex; ++i) {
        threadData[i].shared = &sharedData;
        if (ext_rng_createAndSeed(&threadData[i].rng, rng_algorithm, rng_standardNormal) != 0)
          ext_throwError("could not allocate rng");
        threadData[i].repCellOffset = i * numRepCellsPerThread;
        threadData[i].numRepCells   = numRepCellsPerThread;
        threadData[i].threadId = i;
        threadData[i].results = results + i * numRepCellsPerThread * lossFunctorDef.numResults;
        threadDataPtrs[i] = threadData + i;
      }
      
      for (size_t i = offByOneIndex; i < numThreads; ++i) {
        threadData[i].shared = &sharedData;
        if (ext_rng_createAndSeed(&threadData[i].rng, rng_algorithm, rng_standardNormal) != 0)
          ext_throwError("could not allocate rng");
        threadData[i].repCellOffset = offByOneIndex * numRepCellsPerThread + (i - offByOneIndex) * (numRepCellsPerThread - 1);
        threadData[i].numRepCells   = numRepCellsPerThread - 1;
        threadData[i].threadId = i;
        threadData[i].results = results + threadData[i].repCellOffset * lossFunctorDef.numResults;
        threadDataPtrs[i] = threadData + i;
      }
      
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
  
  void randomSubsampleDivideData(const Data& restrict origData, Data& restrict repData, double* restrict y_test,
                                 ext_rng* restrict generator, size_t* restrict permutation);
  void kFoldDivideData(const Data& restrict origData, Data& restrict repData, double* restrict y_test,
                       size_t k, size_t maxNumFoldObservations, size_t numFullSizedFolds,
                       const size_t* restrict permutation);
}

namespace {
  struct LossFunctorCreatorData {
    const LossFunctorDefinition& lfDef;
    Method method;
    size_t maxNumTestObservations;
    size_t numSamples;
    LossFunctor** result;
  };
  struct LossFunctorDestructorData {
    const LossFunctorDefinition& lfDef;
    LossFunctor* lf;
  };
}

extern "C" void lossFunctorCreatorTask(void* data) {
  LossFunctorCreatorData& lfcd(*static_cast<LossFunctorCreatorData*>(data));
  *lfcd.result = lfcd.lfDef.createFunctor(lfcd.lfDef, lfcd.method, lfcd.maxNumTestObservations, lfcd.numSamples);
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
      if (origData.numObservations % maxNumTestObservations == 0)
        maxNumTrainingObservations = origData.numObservations - maxNumTestObservations;
      else
        maxNumTrainingObservations = origData.numObservations - maxNumTestObservations++;
    } else {
      maxNumTestObservations = origData.numObservations -
        static_cast<size_t>(std::floor(static_cast<double>(origData.numObservations) * sharedData.testSampleSize.p + 0.5));
      maxNumTrainingObservations = origData.numObservations - maxNumTestObservations;
    }
    
    size_t numSamples              = origControl.defaultNumSamples;
        
    const LossFunctorDefinition& lfDef(sharedData.lossFunctorDef);
    bool lossRequiresMutex = lfDef.requiresMutex && !misc_btm_isNull(sharedData.threadManager);
    
    LossFunctor* lf = NULL;
    if (lossRequiresMutex) {
      LossFunctorCreatorData lfcd = { lfDef, sharedData.method, maxNumTestObservations, numSamples, &lf };
      misc_btm_runTaskInParentThread(sharedData.threadManager, threadData.threadId, &lossFunctorCreatorTask, &lfcd);
    } else {
      lf = lfDef.createFunctor(lfDef, sharedData.method, maxNumTestObservations, numSamples);
    }
    
    double* suppliedY_test      = lfDef.y_testOffset      >= 0 ?
                                  *reinterpret_cast<double**>(reinterpret_cast<char*>(lf) + lfDef.y_testOffset) :
                                  NULL;
    double* suppliedTestSamples = sharedData.lossFunctorDef.testSamplesOffset >= 0 ?
                                  *reinterpret_cast<double**>(reinterpret_cast<char*>(lf) + lfDef.testSamplesOffset) :
                                  NULL;
    double* y_test = (suppliedY_test == NULL ? new double[maxNumTestObservations] : suppliedY_test);
    
    Results* samples =
      suppliedTestSamples == NULL ?
        new Results(maxNumTrainingObservations, origData.numPredictors, maxNumTestObservations, numSamples, 1) :
        new Results(maxNumTrainingObservations, origData.numPredictors, maxNumTestObservations, numSamples, 1,
                    new double[numSamples],
                    new double[maxNumTrainingObservations * numSamples],
                    suppliedTestSamples,
                    new double[origData.numPredictors * numSamples]);
    
    Control repControl = origControl;
    repControl.numThreads = 1;
    repControl.verbose = false;
    bool verbose = origControl.verbose;
    
    BARTFit* fit = new BARTFit(repControl, origModel, origData);
    
    Data repData;
    Model repModel;
    
    allocateDataStorage(fit->data, repData, maxNumTrainingObservations, maxNumTestObservations);
    allocateModelStorage(fit->model, repModel);
    
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
    v_threadScratch->samples = samples;
    v_threadScratch->lf = lf;
    v_threadScratch->generator = threadData.rng;
    v_threadScratch->permutation = new size_t[origData.numObservations];
    for (size_t i = 0; i < origData.numObservations; ++i) v_threadScratch->permutation[i] = i;
    
        
    size_t firstCell    = threadData.repCellOffset / sharedData.numReps;
    size_t firstCellRep = threadData.repCellOffset % sharedData.numReps;
    size_t lastCell     = (threadData.repCellOffset + threadData.numRepCells) / sharedData.numReps;
    size_t lastCellRep  = (threadData.repCellOffset + threadData.numRepCells) % sharedData.numReps;
    
    size_t resultIndex = 0;
    xvalData.numBurnIn = sharedData.numInitialBurnIn;
    
    // first and last cells are a bit of a mess, since there can be a lot of off-by-one stuff
    if (firstCellRep != 0) {
      updateFitForCell(*fit, repControl, repModel, sharedData.parameters[firstCell],
                       threadData.threadId, firstCell, sharedData.threadManager, verbose);
      
      for (size_t repIndex = firstCellRep; repIndex < sharedData.numReps; ++repIndex)
      {
        crossvalidate(xvalData, samples, numSamples, threadData.results + resultIndex,
                      lfDef.calculateLoss, sharedData.threadManager, threadData.threadId, lossRequiresMutex, v_threadScratch);

        resultIndex += lfDef.numResults;
        
        xvalData.numBurnIn = sharedData.numRepBurnIn;
      }
      
      ++firstCell;
      firstCellRep = 0;
      
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
    size_t numSamples;
    double* results;
  };
}

extern "C" void lossFunctorTask(void* data) {
  LossFunctorData& lfd(*static_cast<LossFunctorData*>(data));
  lfd.calculateLoss(lfd.lf, lfd.y_test, lfd.numTestObservations, lfd.testSamples, lfd.numSamples, lfd.results);
}

namespace {
  void randomSubsampleCrossvalidate(CrossvalidationData& xvalData,
                                    Results* restrict samples, size_t numSamples, double* restrict results,
                                    LossFunction calculateLoss, misc_btm_manager_t manager, size_t threadId, bool lossRequiresMutex,
                                    ThreadScratch* v_scratch)
  {
    RandomSubsampleThreadScratch& threadScratch(*reinterpret_cast<RandomSubsampleThreadScratch *>(v_scratch));
    
    randomSubsampleDivideData(xvalData.origData, xvalData.repData, threadScratch.y_test,
                              threadScratch.generator, threadScratch.permutation);
    xvalData.fit.setData(xvalData.repData);
    
    xvalData.fit.runSampler(xvalData.numBurnIn, samples);
    
    if (lossRequiresMutex) {
      LossFunctorData ldf = { calculateLoss, *threadScratch.lf, threadScratch.y_test, threadScratch.maxNumTestObservations, samples->testSamples, numSamples, results };
      misc_btm_runTaskInParentThread(manager, threadId, &lossFunctorTask, &ldf);
    } else {
      calculateLoss(*threadScratch.lf, threadScratch.y_test, threadScratch.maxNumTestObservations, samples->testSamples, numSamples, results);
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
      
      kFoldDivideData(xvalData.origData, xvalData.repData, threadScratch.y_test,
                      k, threadScratch.maxNumTestObservations, threadScratch.numFullSizedFolds, threadScratch.permutation);
      xvalData.fit.setData(xvalData.repData);
      
      xvalData.fit.runSampler(xvalData.numBurnIn, samples);
    
      if (lossRequiresMutex) {
        LossFunctorData ldf = { calculateLoss, *threadScratch.lf, threadScratch.y_test, numTestObservations, samples->testSamples, numSamples, foldResults };
        misc_btm_runTaskInParentThread(manager, threadId, &lossFunctorTask, &ldf);
      } else {
        calculateLoss(*threadScratch.lf, threadScratch.y_test, numTestObservations, samples->testSamples, numSamples, foldResults);
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
  
  void randomSubsampleDivideData(const Data& restrict origData, Data& restrict repData, double* restrict y_test,
                                 ext_rng* restrict generator, size_t* restrict permutation)
  {
    size_t i, j, obsIndex;
    double* restrict y = const_cast<double* restrict>(repData.y);
    double* restrict x = const_cast<double* restrict>(repData.x);
    double* restrict x_test = const_cast<double* restrict>(repData.x_test);
    
    size_t numTrainingObservations = repData.numObservations;
    size_t numTestObservations     = repData.numTestObservations;
    
    permuteIndexArray(generator, permutation, origData.numObservations);
    std::sort(permutation, permutation + numTestObservations);
    std::sort(permutation + numTestObservations, permutation + origData.numObservations);
    
    for (i = 0; i < numTestObservations; ++i) {
      obsIndex = *permutation++;
      y_test[i] = origData.y[obsIndex];
      for (j = 0; j < origData.numPredictors; ++j) {
        x_test[i + j * numTestObservations] = origData.x[obsIndex + j * origData.numObservations];
      }
    }
    for (i = 0; i < numTrainingObservations; ++i) {
      obsIndex = *permutation++;
      y[i] = origData.y[obsIndex];
      for (j = 0; j < origData.numPredictors; ++j) {
        x[i + j * numTrainingObservations] = origData.x[obsIndex + j * origData.numObservations];
      }
    }
  }
  
  void kFoldDivideData(const Data& restrict origData, Data& restrict repData, double* restrict y_test,
                       size_t k, size_t maxNumFoldObservations, size_t numFullSizedFolds,
                       const size_t* restrict permutation)
  {
    size_t i, j, obsIndex, foldStartIndex;
    double* restrict y = const_cast<double* restrict>(repData.y);
    double* restrict x = const_cast<double* restrict>(repData.x);
    double* restrict x_test = const_cast<double* restrict>(repData.x_test);
    
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
    }
    for (i = 0; i < foldStartIndex; ++i) {
      obsIndex = permutation[i];
      y[i] = origData.y[obsIndex];
      for (j = 0; j < origData.numPredictors; ++j) {
        x[i + j * numTrainingObservations] = origData.x[obsIndex + j * origData.numObservations];
      }
    }
    for ( /* */; i < numTrainingObservations; ++i) {
      obsIndex = permutation[i + numTestObservations];
      y[i] = origData.y[obsIndex];
      for (j = 0; j < origData.numPredictors; ++j) {
        x[i + j * numTrainingObservations] = origData.x[obsIndex + j * origData.numObservations];
      }
    }
  }
  
  void allocateDataStorage(const Data& origData, Data& repData, size_t numTrainingObservations, size_t numTestObservations)
  {
    repData.y = new double[numTrainingObservations];
    repData.x = new double[numTrainingObservations * origData.numPredictors];
    repData.x_test = new double[numTestObservations * origData.numPredictors];
     
    repData.weights    = origData.weights != NULL ? new double[numTrainingObservations]  : NULL;
    repData.offset     = origData.offset  != NULL ? new double[numTrainingObservations]  : NULL;
    repData.testOffset = origData.offset  != NULL ? new double[numTestObservations] : NULL;
    
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
    
    CGMPrior* repTreePrior = new CGMPrior();
    const CGMPrior* oldTreePrior = static_cast<CGMPrior*>(origModel.treePrior);
    repTreePrior->base = oldTreePrior->base;
    repTreePrior->power = oldTreePrior->power;
    
    repModel.treePrior = repTreePrior;
    
    
    NormalPrior* repNodePrior = new NormalPrior();
    const NormalPrior* oldNodePrior = static_cast<NormalPrior*>(origModel.muPrior);
    repNodePrior->precision = oldNodePrior->precision;
    
    repModel.muPrior = repNodePrior;
    
    
    ChiSquaredPrior* repResidPrior = new ChiSquaredPrior();
    const ChiSquaredPrior* oldResidPrior = static_cast<ChiSquaredPrior*>(origModel.sigmaSqPrior);
    repResidPrior->degreesOfFreedom = oldResidPrior->degreesOfFreedom;
    repResidPrior->scale = oldResidPrior->scale;
    
    repModel.sigmaSqPrior = repResidPrior;
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
          ext_printf("    [%lu, %lu] n.tree: %lu, k: %f, power: %f, base: %f\n",
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
  ext_printf("    [%lu, %lu] n.trees: %lu, k: %.2f, power: %.2f, base: %.2f\n",
             data.threadId + 1, data.cellIndex + 1,
             data.numTrees, data.k, data.power, data.base);
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
        ext_printf("    [%lu] n.trees: %lu, k: %.2f, power: %.2f, base: %.2f\n", cellIndex, numTrees, k, power, base);
      } else {
        PrintData printData = { threadId, cellIndex, numTrees, k, power, base };
        misc_btm_runTaskInParentThread(manager, threadId, &printTask, &printData);
      }
    }
    
    repControl.numTrees = numTrees;
    
    double endNodeSd = (repControl.responseIsBinary ? 3.0 : 0.5) / (k * std::sqrt(static_cast<double>(repControl.numTrees)));
    static_cast<NormalPrior*>(repModel.muPrior)->precision = 1.0 / (endNodeSd * endNodeSd);
    
    static_cast<CGMPrior*>(repModel.treePrior)->power = power;
    static_cast<CGMPrior*>(repModel.treePrior)->base  = base;
    
    fit.setControl(repControl);
    fit.setModel(repModel);
  }
}

