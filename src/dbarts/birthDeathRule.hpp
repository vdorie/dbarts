#ifndef DBARTS_BIRTH_DEATH_RULE_HPP
#define DBARTS_BIRTH_DEATH_RULE_HPP

#include <cstddef>

namespace dbarts {
  struct BARTFit;
  struct Tree;
  
  double birthOrDeathNode(const BARTFit& fit, std::size_t chainNum, Tree& tree, const double* y, double sigma, bool* stepWasTaken, bool* birthedNode);
}
  

#endif
