#include "config.hpp"
#include <dbarts/model.hpp>

#include <cmath>

#include <external/random.h>
#include <external/stats.h>

#include <dbarts/control.hpp>
#include "node.hpp"

using std::size_t;
using std::uint32_t;

namespace dbarts {
  NormalPrior::NormalPrior(const Control& control, double k)
  {
    double sigma = (control.responseIsBinary ? 3.0 : 0.5) /  (k * std::sqrt(static_cast<double>(control.numTrees)));
    precision = 1.0 / (sigma * sigma);
  }
  
  double NormalPrior::drawFromPosterior(ext_rng* rng, double ybar, double numEffectiveObservations, double residualVariance) const {
    double posteriorPrecision = numEffectiveObservations / residualVariance;
    
    double posteriorMean = posteriorPrecision * ybar / (this->precision + posteriorPrecision);
    double posteriorSd   = 1.0 / std::sqrt(this->precision + posteriorPrecision);
    
    return posteriorMean + posteriorSd * ext_rng_simulateStandardNormal(rng);
  }
  
  double NormalPrior::computeLogIntegratedLikelihood(const BARTFit& fit, size_t chainNum, const Node& node, const double* y, double residualVariance) const
  {
    size_t numObservationsInNode = node.getNumObservations();
    if (numObservationsInNode == 0) return 0.0;
    
    double y_bar = node.getAverage();
    double var_y = node.computeVariance(fit, chainNum, y);
      
    double posteriorPrecision = node.getNumEffectiveObservations() / residualVariance;
    
    double result;
    // we divide out by the actual n - 1 instead of a weighted sum because the variance
    // divides that number regardless; the calculations are correct in both cases
    result  = 0.5 * std::log(this->precision / (this->precision + posteriorPrecision));
    result -= 0.5 * (var_y / residualVariance) * static_cast<double>(numObservationsInNode - 1);
    result -= 0.5 * ((this->precision * y_bar) * (posteriorPrecision * y_bar)) / (this->precision + posteriorPrecision);
    
    return result;
  }
  
  ChiSquaredPrior::ChiSquaredPrior(double degreesOfFreedom, double quantile) :
    degreesOfFreedom(degreesOfFreedom),
    scale(ext_quantileOfChiSquared(1.0 - quantile, degreesOfFreedom) / degreesOfFreedom)
  {
  }
  
  double ChiSquaredPrior::drawFromPosterior(ext_rng* rng, double numObservations, double sumOfSquaredResiduals) const {
    double posteriorDegreesOfFreedom = degreesOfFreedom + numObservations;
    
    double posteriorScale = degreesOfFreedom * scale + sumOfSquaredResiduals;
    
    return posteriorScale / ext_rng_simulateChiSquared(rng, posteriorDegreesOfFreedom);
  }
}
