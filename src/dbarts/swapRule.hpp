#ifndef DBARTS_SWAP_RULE_HPP
#define DBARTS_SWAP_RULE_HPP

namespace dbarts {
  struct BARTFit;
  struct Tree;
  
  double swapRule(const BARTFit& fit, Tree& tree, const double* y, bool* stepTaken);
}

#endif
