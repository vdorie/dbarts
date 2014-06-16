#ifndef BART_BART_FIT_HPP
#define BART_BART_FIT_HPP

#include <cstddef> // size_t
#include "cstdint" // uint32_t

#include <external/thread.h>

#include "types.hpp"

#include "control.hpp"
#include "data.hpp"
#include "model.hpp"
#include "scratch.hpp"
#include "state.hpp"

namespace bart {
  struct Results;
  
  struct BARTFit {
    Control control; // top three are passed in from elsewhere
    Model model;
    Data data;
    
    Scratch scratch;
    State state;
    
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
