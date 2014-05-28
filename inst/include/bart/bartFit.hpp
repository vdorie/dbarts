#ifndef BART_BART_FIT_HPP
#define BART_BART_FIT_HPP

#include <cstddef> // size_t
#include "cstdint" // uint32_t

#include <external/thread.h>

#include "types.hpp"

#include "control.hpp"
#include "data.hpp"
#include "model.hpp"

namespace bart {
  struct Results;
  
  struct Tree;
  
  struct ScaleFactor {
    double min, max, range;
  };
  
  struct BARTFit {
    Control control;
    Model model;
    Data data;
    
    double* yRescaled;
    double* Xt; // used for pulling out rows of X
    double* Xt_test;
    double* treeY;
    double* weights;
    
    VariableType* variableTypes;
    uint32_t* numCutsPerVariable;
    double** cutPoints; // used to be RuleMat
    
    Tree* trees;
    size_t* treeIndices;
    
    double* treeFits;      // numObs x numTrees;     vals for tree <=> x + i * numObs
    double* totalFits;
    double* totalTestFits; // numTestObs x 1
    
    
    ScaleFactor dataScale;
    double sigma;
    
    ext_mt_manager_t threadManager;
    
    BARTFit(Control control, Model model, Data data);
    ~BARTFit();
    
    
    // creates a results object
    Results* runSampler();
    Results* runSampler(std::size_t numBurnIn, std::size_t numSamples);
    
    // general work flow:
    //   sample treeThinningRate trees
    //   store current test/training predictions, sigma
    //   pass to call back function
    //
    // call back can change Y, but must call the following. can even
    // be changed in place
    void setResponse(const double* newResponse);
  };
} // namespace bart

#endif // BART_BART_FIT_HPP
