#ifndef BART_CHANGE_RULE_HPP
#define BART_CHANGE_RULE_HPP

namespace bart {
  struct BARTFit;
  struct Tree;
  
  double changeRule(const BARTFit& fit, Tree& tree, const double* y, bool* stepTaken);
}  

#endif
