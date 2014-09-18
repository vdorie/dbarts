#ifndef DBARTS_BART_FIT_HPP
#define DBARTS_BART_FIT_HPP

#include <cstddef> // size_t
#include "cstdint.hpp" // uint32_t

#include <external/thread.h>

#include "types.hpp"

#include "control.hpp"
#include "data.hpp"
#include "model.hpp"
#include "scratch.hpp"
#include "state.hpp"

namespace dbarts {
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
    
    Results* runSampler();
    Results* runSampler(std::size_t numBurnIn, std::size_t numSamples);
    
    
    // setResponse and setOffset replace the local copy (or can even be modified 
    // in place), however because setPredictor only works on a single column it
    // necessarily modifies in place
    void setResponse(const double* newResponse); 
    void setOffset(const double* newOffset);
    void setPredictor(const double* newPredictor, std::size_t predictorColumn);
    void setTestPredictor(const double* newTestPredictor, std::size_t predictorColumn);
    
    void setTestPredictors(const double* newTestPredictor, std::size_t numTestObservations);
    void setTestOffset(const double* newTestOffset);
    void setTestPredictors(const double* newTestPredictor, const double* newTestOffset, std::size_t numTestObservations);
    
    bool saveToFile(const char* fileName) const;
    BARTFit* BARTFit::loadFromFile(const char* fileName);
  };
} // namespace dbarts

#endif // DBARTS_BART_FIT_HPP
