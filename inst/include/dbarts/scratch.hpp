// stuff that can be figured out from Control, Model, and Data
#ifndef DBARTS_SCRATCH_HPP
#define DBARTS_SCRATCH_HPP

#include "cstdint" // int types

namespace dbarts {
  struct ScaleFactor { double min, max, range; };
  
  struct Scratch {
    const double* yRescaled;
    const double* Xt;
    const double* Xt_test;
    double* treeY;
    
    ScaleFactor dataScale;
    
    const std::uint32_t* numCutsPerVariable;
    const double* const* cutPoints;
  };
} // namespace dbarts

#endif // DBARTS_SCRATCH_HPP
