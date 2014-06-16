// stuff that can be figured out from Control, Model, and Data
#ifndef BART_SCRATCH_HPP
#define BART_SCRATCH_HPP

#include "cstdint" // int types

namespace bart {
  struct ScaleFactor { double min, max, range; };
  
  struct Scratch {
    const double* yRescaled;
    const double* Xt;
    const double* Xt_test;
    double* treeY;
    
    ScaleFactor dataScale;
    
    const uint32_t* numCutsPerVariable;
    const double* const* cutPoints;
  };
} // namespace bart

#endif // BART_SCRATCH_HPP
