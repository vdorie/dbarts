#ifndef BART_BIRTH_DEATH_RULE_HPP
#define BART_BIRTH_DEATH_RULE_HPP

namespace bart {
  struct BARTFit;
  struct Tree;
  
  double birthOrDeathNode(const BARTFit& fit, Tree& tree, const double* y, bool* stepWasTaken, bool* birthedNode);
}
  

#endif
