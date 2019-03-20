#ifndef DBARTS_BART_FIT_HPP
#define DBARTS_BART_FIT_HPP

#include <cstddef> // size_t
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
    Control control; // top three are passed in from elsewhere
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
    
    void setRNGState(const void* const* uniformState, const void* const* normalState);
    
    Results* runSampler();
    Results* runSampler(std::size_t numBurnIn, std::size_t numSamples);
    void runSampler(std::size_t numBurnIn, Results* results);
    
    
    void predict(const double* x_test, std::size_t numTestObservations, const double* testOffset, double* result) const;
    // settors simply replace local pointers to variables. dimensions much match
    // update modifies the local copy (which may belong to someone else)
    void setResponse(const double* newResponse); 
    void setOffset(const double* newOffset);
    
    /* These functions change the predictors or the cut points derived from them, and thus can change the
     * trees. Consequently, they have differring semantics depending on intended use.
     *
     *   setPredictor/updatePredictor - for use in a sampler; set replaces 'x' in memory with
     *                                  a different matrix, update modifies existing
     *   setCutPoints - for use at the beginning of a process so that the cut points aren't
     *                  derived from a random sample
     *   setData - for use when fitting a model 'near' another model, so as to preserve 
     *             the tree structure as a starting point
     *
     * Changing a predictor can mean that some nodes are empty and the tree would have 0 prior
     * probability. Consequently, setPredictor and updatePredictor can roll back the change
     * and be used in a rejection-sampling scheme.
     * 
     * Changing cut points can invalidate trees by leaving a node empty as before, but if the number of
     * cut points change it is also possible to cause a tree to prematurely exhaust all available splits.
     *
     * setPredictor, updatePredictor, and setCutPoints - as they are meant to be used within a single
     * context - view the cutpoints as real numbers and don't attempt to preserve the tree structure.
     * setData will attempt to map the old cut points onto new ones so as to preserve as much of
     * the tree as possible.
     */
    // updateCutPoints == true uses the default rule and max number of cut points set elsewhere
    bool setPredictor(const double* newPredictor, bool forceUpdate, bool updateCutPoints);
    bool updatePredictor(const double* newPredictor, const std::size_t* columns, std::size_t numColumns, bool forceUpdate, bool updateCutPoints); 
    void setCutPoints(const double* const* cutPoints, const std::uint32_t* numCutPoints, const std::size_t* columns, std::size_t numColumns);
    void setData(const Data& data);
    
    void setTestPredictor(const double* newTestPredictor, std::size_t numTestObservations);
    void setTestOffset(const double* newTestOffset);
    void setTestPredictorAndOffset(const double* newTestPredictor, const double* newTestOffset, std::size_t numTestObservations);
    
    void updateTestPredictor(const double* newTestPredictor, std::size_t column);
    void updateTestPredictors(const double* newTestPredictor, const std::size_t* columns, std::size_t numColumns);
    
    void sampleTreesFromPrior();
    
    void rebuildScratchFromState();
    
    void printTrees(const std::size_t* chains, std::size_t numChains,
                    const std::size_t* samples, std::size_t numSamples,
                    const std::size_t* indices, std::size_t numIndices) const;
    FlattenedTrees* getFlattenedTrees(const size_t* chainIndices, size_t numChainIndices,
                                      const size_t* sampleIndices, size_t numSampleIndices,
                                      const size_t* treeIndices, size_t numTreeIndices) const;
    
    
    // the new control must have the same number of chains as the previous or else prob seg fault
    void setControl(const Control& control);
    void setModel(const Model& model);
    
    bool saveToFile(const char* fileName) const;
    static BARTFit* loadFromFile(const char* fileName);
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

