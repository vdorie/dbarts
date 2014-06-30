#ifndef DBARTS_BIRTH_DEATH_RULE_HPP
#define DBARTS_BIRTH_DEATH_RULE_HPP

namespace dbarts {
  struct BARTFit;
  struct Tree;
  
  double birthOrDeathNode(const BARTFit& fit, Tree& tree, const double* y, bool* stepWasTaken, bool* birthedNode);
}
  

#endif
