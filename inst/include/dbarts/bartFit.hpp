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
    
    // returns true/false if update is successful; if forceUpdate is true, will prune empty leaves
    bool setPredictor(const double* newPredictor, bool forceUpdate, bool updateCutPoints);
    bool updatePredictor(const double* newPredictor, const std::size_t* columns, std::size_t numColumns, bool forceUpdate, bool updateCutPoints); 
    
    void setCutPoints(const double* const* cutPoints, const std::uint32_t* numCutPoints, const std::size_t* columns, std::size_t numColumns);
    
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
    
    // this assumes that the new data has as many predictors as the old, and that they correspond to each other;
    // it'll attempt to map cut points from the old to the new, and prune any trees that may have been left in an
    // invalid state
    void setData(const Data& data);
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

