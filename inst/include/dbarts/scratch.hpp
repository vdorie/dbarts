// stuff that can be figured out from Control, Model, and Data
#ifndef DBARTS_SCRATCH_HPP
#define DBARTS_SCRATCH_HPP

#include <cstddef>
#include "cstdint.hpp" // int types
#include "types.hpp"

namespace dbarts {
  struct ScaleFactor { double min, max, range; };
  
  struct SharedScratch {
    const double* yRescaled;
    const xint_t* xt; // x transpose
    const xint_t* xt_test;
    // const double* xt; // x transpose
    //const double* xt_test;
    
    ScaleFactor dataScale;
    
    const std::uint32_t* numCutsPerVariable;
    const double* const* cutPoints;
  };
  struct ChainScratch {
    double* treeY;
    double* probitLatents;
    
    double* totalFits;     // numObs
    double* totalTestFits; // numTestObs
    
    std::size_t taskId;
  };
} // namespace dbarts

#endif // DBARTS_SCRATCH_HPP

