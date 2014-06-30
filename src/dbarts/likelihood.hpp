#ifndef DBARTS_LIKELIHOOD_HPP
#define DBARTS_LIKELIHOOD_HPP

namespace dbarts {
  struct BARTFit;
  struct Node;
  
  double computeLogLikelihoodForBranch(const BARTFit& fit, const Node& branch, const double* y);
}

#endif
