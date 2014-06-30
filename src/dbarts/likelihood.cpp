#include "config.hpp"
#include "likelihood.hpp"

#include <cstddef>

#include <dbarts/bartFit.hpp>
#include <dbarts/model.hpp>
#include <dbarts/state.hpp>
#include "node.hpp"

using std::size_t;

namespace dbarts {
  double computeLogLikelihoodForBranch(const BARTFit& fit, const Node& branch, const double* y)
  {
    NodeVector bottomVector(branch.getBottomVector());
    size_t numBottomNodes = bottomVector.size();
    
    double logProbability = 0.0;
    for (size_t i = 0; i < numBottomNodes; ++i) {
      const Node& bottomNode(*bottomVector[i]);
      
      if (bottomNode.getNumObservations() == 0) return -10000000.0;
      
      logProbability += fit.model.muPrior->computeLogIntegratedLikelihood(fit, bottomNode, y, fit.state.sigma * fit.state.sigma);
    }
    
    return logProbability;
  }
}
