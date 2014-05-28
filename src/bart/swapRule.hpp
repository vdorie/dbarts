#ifndef BART_SWAP_RULE_HPP
#define BART_SWAP_RULE_HPP

namespace bart {
  struct BARTFit;
  struct Tree;
  
  double swapRule(const BARTFit& fit, Tree& tree, const double* y, bool* stepTaken);
}

#endif
