#ifndef BART_LIKELIHOOD_HPP
#define BART_LIKELIHOOD_HPP

namespace bart {
  struct BARTFit;
  struct Node;
  
  double computeLogLikelihoodForBranch(const BARTFit& fit, const Node& branch, const double* y);
}

#endif
