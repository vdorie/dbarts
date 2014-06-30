#ifndef DBARTS_CHANGE_RULE_HPP
#define DBARTS_CHANGE_RULE_HPP

namespace dbarts {
  struct BARTFit;
  struct Tree;
  
  double changeRule(const BARTFit& fit, Tree& tree, const double* y, bool* stepTaken);
}  

#endif
