#include "config.hpp"
#include <bart/model.hpp>

#include <cmath>

#include <external/stats.h>

#include <bart/control.hpp>
#include "node.hpp"

using std::size_t;
using std::uint32_t;

namespace bart {
  NormalPrior::NormalPrior(const Control& control) :
    precision(1.0)
  {
    double sigma = (control.responseIsBinary ? 3.0 : 0.5) /  (control.kFactor * std::sqrt((double) control.numTrees));
    precision = 1.0 / (sigma * sigma);
  }
  
  double NormalPrior::drawFromPosterior(double ybar, size_t numObservations, double residualVariance) const {
    double dataWeight = ((double) numObservations) / residualVariance;
    
    double posteriorMean = dataWeight * ybar / (precision + dataWeight);
    double posteriorSd   = 1.0 / std::sqrt(precision + dataWeight);
    
    return posteriorMean + posteriorSd * ext_simulateStandardNormal();
  }
  
  double NormalPrior::computeLogIntegratedLikelihood(const BARTFit& fit, const Node& node, const double* y, double residualVariance) const {
    size_t numObservationsInNode = node.getNumObservationsInNode();
    
    double dataWeight = 0.0, varianceInNode = 0.0, ybar = 0.0;
    
    if (numObservationsInNode > 0) {
      ybar = node.getAverage();
      varianceInNode = node.computeVariance(fit, y);
      
      dataWeight = (double) numObservationsInNode / residualVariance;
    }
    
    double result = 0.5 * std::log(precision / (precision + dataWeight));
    result -= 0.5 * (varianceInNode / residualVariance) * (double) (numObservationsInNode - 1);
    result -= 0.5 * (precision * dataWeight * ybar * ybar) / (precision + dataWeight);
    return result;
  }
  
  ChiSquaredPrior::ChiSquaredPrior(const Control& control) :
    degreesOfFreedom((double) control.sigmaDf),
    scale(ext_quantileOfChiSquared(1.0 - control.sigmaQuantile, (double) control.sigmaDf) / (double) control.sigmaDf)
  {
  }
  
  double ChiSquaredPrior::drawFromPosterior(size_t numObservations, double sumOfSquaredResiduals) const {
    uint32_t posteriorDegreesOfFreedom = degreesOfFreedom + numObservations;
    
    double posteriorScale = ((double) degreesOfFreedom * scale + sumOfSquaredResiduals);
    
    return posteriorScale / ext_simulateChiSquared((double) posteriorDegreesOfFreedom);
  }
}
