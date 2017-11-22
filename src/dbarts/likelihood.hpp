#ifndef DBARTS_LIKELIHOOD_HPP
#define DBARTS_LIKELIHOOD_HPP

#include <cstddef>

namespace dbarts {
  struct BARTFit;
  struct Node;
  
  double computeLogLikelihoodForBranch(const BARTFit& fit, std::size_t chainNum, const Node& branch, const double* y, double sigma);
}

#endif
