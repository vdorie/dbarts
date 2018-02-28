#ifndef DBARTS_SWAP_RULE_HPP
#define DBARTS_SWAP_RULE_HPP

#include <cstddef>

namespace dbarts {
  struct BARTFit;
  struct Tree;
  
  double swapRule(const BARTFit& fit, std::size_t chainNum, Tree& tree, const double* y, double sigma, bool* stepTaken);
}

#endif
