#include "config.hpp"
#include "likelihood.hpp"

#include <cstddef>

#include <bart/bartFit.hpp>
#include <bart/model.hpp>
#include <bart/state.hpp>
#include "node.hpp"

using std::size_t;

namespace bart {
  double computeLogLikelihoodForBranch(const BARTFit& fit, const Node& branch, const double* y)
  {
    NodeVector bottomVector(branch.getBottomVector());
    size_t numBottomNodes = bottomVector.size();
    
    double logProbability = 0.0;
    for (size_t i = 0; i < numBottomNodes; ++i) {
      const Node& bottomNode(*bottomVector[i]);
      
      if (bottomNode.getNumObservationsInNode() == 0) return -10000000.0;
      
      logProbability += fit.model.muPrior->computeLogIntegratedLikelihood(fit, bottomNode, y, fit.state.sigma * fit.state.sigma);
    }
    
    return logProbability;
  }
}
