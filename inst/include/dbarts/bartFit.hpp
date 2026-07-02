#ifndef DBARTS_BART_FIT_HPP
#define DBARTS_BART_FIT_HPP

#include <cstddef>            // size_t
#include <dbarts/cstdint.hpp> // uint32_t

#include <dbarts/control.hpp>
#include <dbarts/data.hpp>
#include <dbarts/model.hpp>
#include <dbarts/scratch.hpp>
#include <dbarts/state.hpp>
#include <dbarts/types.hpp>

extern "C" {
struct _misc_htm_manager_t;
typedef struct _misc_htm_manager_t* misc_htm_manager_t;
}

namespace dbarts {
  struct Results;
  struct SharedScratch;
  struct FlattenedTrees;

  struct BARTFit {
    // control, model, and data are supplied by the caller at construction.
    Control control;
    Model model;
    Data data;

    SharedScratch sharedScratch;
    ChainScratch* chainScratch;
    State* state;

    double runningTime;
    std::size_t currentNumSamples;
    std::size_t currentSampleNum;

    misc_htm_manager_t threadManager;

    const std::uint32_t* numCutsPerVariable;
    const double* const* cutPoints;

    BARTFit(Control control, Model model, Data data);
    ~BARTFit();

    void setRNGState(
      const void* const* uniformState,
      const void* const* normalState
    );

    Results* runSampler();
    Results* runSampler(
      std::size_t numBurnIn,
      std::size_t numThreads,
      std::size_t numSamples
    );
    void
    runSampler(std::size_t numBurnIn, std::size_t numThreads, Results* results);

    void predict(
      const double* x_test,
      std::size_t numTestObservations,
      const double* testOffset,
      double* result
    ) const;
    void predict(
      const double* x_test,
      std::size_t numTestObservations,
      const double* testOffset,
      size_t numThreads,
      double* result
    ) const;
    // These replace the model's copy with the given quantity; dimensions must
    // match.
    void setResponse(const double* newResponse);
    void setOffset(const double* newOffset, bool updateScale);
    void setWeights(const double* newWeights);

    // These set the value in every chain and ignore any fixed prior; installing
    // a new fixed prior also sets the parameters. The pointer overloads take
    // one value per chain.
    void setSigma(double newSigma);
    void setK(double newK);
    void setSigma(const double* newSigma);
    void setK(const double* newK);

    // These change the predictors, or the cut points derived from them, and so
    // can change the trees; their semantics differ by intended use.
    // set/updatePredictor are for a Gibbs sampler (set swaps 'x' in memory,
    // update modifies it in place); setCutPoints is for the start of a run, so
    // the cut points aren't taken from a random sample; setData is for fitting
    // 'near' another model, preserving as much tree structure as possible.
    //
    // A change can leave a node empty, giving the tree 0 prior probability, so
    // set/updatePredictor can roll back and be used for rejection sampling. The
    // per-observation variants assume observation-level independence, so they
    // install what they can and roll back what they can't, reporting the
    // outcome for each observation in the supplied boolean array.
    //
    // Changing the number of cut points can also prematurely exhaust a tree's
    // splits. set/updatePredictor and setCutPoints treat cut points as reals
    // and don't preserve tree structure; setData maps old cut points onto new
    // ones.

    /// Replaces \c x in memory with a new matrix, rolling back the change if it
    /// would leave a leaf empty.
    /// \param forceUpdate If \c true, apply the change even when it empties a
    ///   leaf (the empty leaf is then pruned); if \c false, an invalid change
    ///   is rolled back, so this can drive a rejection sampler.
    /// \param updateCutPoints If \c true, derive fresh cut points from the new
    ///   predictor with the default rule (max count set elsewhere) and
    ///   rebalance observations.
    /// \returns \c true if the change was kept, \c false if it was rolled back
    ///   (only possible when \p forceUpdate is \c false).
    bool setPredictor(
      const double* newPredictor,
      bool forceUpdate,
      bool updateCutPoints
    );
    /// Modifies the given columns of \c x in place, rolling back the change if
    /// it would leave a leaf empty.
    /// \param columns The \p numColumns column indices to replace.
    /// \param forceUpdate As in \c setPredictor.
    /// \param updateCutPoints As in \c setPredictor.
    /// \returns \c true if the change was kept, \c false if it was rolled back
    ///   (only possible when \p forceUpdate is \c false).
    bool updatePredictor(
      const double* newPredictor,
      const std::size_t* columns,
      std::size_t numColumns,
      bool forceUpdate,
      bool updateCutPoints
    );
    /// Per-observation rollback for a single column: installs each
    /// observation's new value on its own, rejecting only those that would
    /// empty a leaf. Observations are visited in a random order drawn from the
    /// first chain's generator.
    /// \param[out] installed One entry per observation, \c true where the new
    ///   value was kept and \c false where it was rolled back.
    void updatePredictorPerObservation(
      const double* newColumn,
      std::size_t column,
      bool* installed
    );
    /// Joint per-observation rollback across several index-aligned samplers: a
    /// single sweep installs each observation in every sampler or in none, so
    /// no sampler is left with an empty leaf and the samplers never disagree.
    /// \param samplers The \p numFits samplers, which share index-aligned
    ///   observations.
    /// \param columns \p columns[j] is the column to replace in sampler \c j;
    ///   these may differ across samplers.
    /// \param[out] installed One entry per observation, identical across all
    ///   fits; \c true where the new value was kept.
    static void updatePredictorPerObservationJointly(
      BARTFit** fits,
      std::size_t numFits,
      const double* newColumn,
      const std::size_t* columns,
      bool* installed
    );
    /// Installs explicit cut points for the given columns instead of deriving
    /// them from the data.
    /// \param cutPoints \p cutPoints[j] holds \p numCutPoints[j] ascending cut
    ///   points for column \p columns[j].
    /// \param columns The \p numColumns column indices to set.
    void setCutPoints(
      const double* const* cutPoints,
      const std::uint32_t* numCutPoints,
      const std::size_t* columns,
      std::size_t numColumns
    );
    /// Swaps in new data, mapping the old cut points onto the new ones to
    /// preserve as much tree structure as possible.
    void setData(const Data& data);

    void setTestPredictor(
      const double* newTestPredictor,
      std::size_t numTestObservations
    );
    void setTestOffset(const double* newTestOffset);
    void setTestPredictorAndOffset(
      const double* newTestPredictor,
      const double* newTestOffset,
      std::size_t numTestObservations
    );

    void
    updateTestPredictor(const double* newTestPredictor, std::size_t column);
    void updateTestPredictors(
      const double* newTestPredictor,
      const std::size_t* columns,
      std::size_t numColumns
    );
    void storeLatents(double* target) const;

    /// \returns The number of threads actually started.
    std::size_t startThreads();
    /// \returns The number of threads actually started.
    std::size_t startThreads(std::size_t numThreads);
    void stopThreads();

    void sampleTreesFromPrior();
    void sampleNodeParametersFromPrior();

    void rebuildScratchFromState();

    void printTrees(
      const std::size_t* chainIndices,
      std::size_t numChainIndices,
      const std::size_t* sampleIndices,
      std::size_t numSampleIndices,
      const std::size_t* treeIndices,
      std::size_t numTreeIndices
    ) const;
    // When useLiveTrees is true, returns the current working trees even for a
    // keepTrees sampler; sampleIndices are then ignored.
    FlattenedTrees* getFlattenedTrees(
      const std::size_t* chainIndices,
      std::size_t numChainIndices,
      const std::size_t* sampleIndices,
      std::size_t numSampleIndices,
      const std::size_t* treeIndices,
      std::size_t numTreeIndices,
      bool useLiveTrees
    ) const;
    // Similar to getFlattenedTrees, but routes the supplied x_test through each
    // tree so numObservations counts that data instead of the training
    // predictors.
    FlattenedTrees* getFlattenedTreesCountingData(
      const std::size_t* chainIndices,
      std::size_t numChainIndices,
      const std::size_t* sampleIndices,
      std::size_t numSampleIndices,
      const std::size_t* treeIndices,
      std::size_t numTreeIndices,
      bool useLiveTrees,
      const double* x_test,
      std::size_t numTestObservations
    ) const;

    /// \warning The new control must have the same number of chains as the old,
    ///   or it may segfault.
    void setControl(const Control& control);
    void setModel(const Model& model);

    void printInitialSummary() const;
  };

  struct FlattenedTrees {
    std::size_t totalNumNodes;
    std::size_t* chainNumber;
    std::size_t* sampleNumber;
    std::size_t* treeNumber;
    std::size_t* numObservations;
    std::int32_t* variable;
    double* value;

    FlattenedTrees(std::size_t totalNumNodes);
    ~FlattenedTrees();
  };
} // namespace dbarts

#endif // DBARTS_BART_FIT_HPP
