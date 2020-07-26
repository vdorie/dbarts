#include "config.hpp"
#include <dbarts/model.hpp>

#include <cmath>
#include <limits>

#include <external/io.h>
#include <external/random.h>
#include <external/stats.h>
#include <misc/stats.h>

#include <dbarts/bartFit.hpp>
#include <dbarts/control.hpp>

#include "node.hpp"
#include "tree.hpp"

using std::size_t;
using std::uint32_t;

namespace dbarts {
  NormalPrior::NormalPrior(const Control& control, const Model& model, double k)
  {
    double sigma = model.nodeScale / (k * std::sqrt(static_cast<double>(control.numTrees)));
    precision = 1.0 / (sigma * sigma);
  }
  
  double NormalPrior::drawFromPosterior(ext_rng* rng, double ybar, double numEffectiveObservations, double residualVariance) const {
    double posteriorPrecision = numEffectiveObservations / residualVariance;
    
    double posteriorMean = posteriorPrecision * ybar / (this->precision + posteriorPrecision);
    double posteriorSd   = 1.0 / std::sqrt(this->precision + posteriorPrecision);
    
    return posteriorMean + posteriorSd * ext_rng_simulateStandardNormal(rng);
  }
  
  double NormalPrior::drawFromPrior(ext_rng* rng) const {
    return ext_rng_simulateStandardNormal(rng) / std::sqrt(this->precision);
  }
  
  double NormalPrior::computeLogIntegratedLikelihood(const BARTFit& fit, size_t chainNum, const Node& node, const double* y, double residualVariance) const
  {
    size_t numObservationsInNode = node.getNumObservations();
    if (numObservationsInNode == 0) return 0.0;
    
    double y_bar = node.getAverage();
    double var_y = node.computeVariance(fit, chainNum, y);
      
    double posteriorPrecision = node.getNumEffectiveObservations() / residualVariance;
    
    
    // symbolically:
    //   tau.mu = 1 / sigma.sq.mu
    //   tau.y  = n / sigma.sq.y
    //   0.5 * log ( tau.mu / (tau.mu + tau.y) ) - 0.5 * SSR / sigma.sq.y -
    //     y.bar^2 * tau.mu * tau.y / (tau.mu + tau.y)
    double result;
    // we divide out by the actual n - 1 instead of a weighted sum because the variance
    // divides that number regardless; the calculations are correct in both cases
    result  = 0.5 * std::log(this->precision / (this->precision + posteriorPrecision));
    result -= 0.5 * (var_y / residualVariance) * static_cast<double>(numObservationsInNode - 1);
    result -= 0.5 * ((this->precision * y_bar) * (posteriorPrecision * y_bar)) / (this->precision + posteriorPrecision);
    
    return result;
  }
  
  void NormalPrior::setK(const BARTFit& fit, double k) {
    double sigma = fit.model.nodeScale / (k * std::sqrt(static_cast<double>(fit.control.numTrees)));
    precision = 1.0 / (sigma * sigma);
  }
  
  double NormalPrior::getK(const BARTFit& fit) const {
    double sigma = std::sqrt(1.0 / precision);
    return fit.model.nodeScale / (sigma * std::sqrt(static_cast<double>(fit.control.numTrees)));
  }
  
  ChiSquaredPrior::ChiSquaredPrior(double degreesOfFreedom, double quantile) :
    degreesOfFreedom(degreesOfFreedom),
    scale(ext_quantileOfChiSquared(1.0 - quantile, degreesOfFreedom) / degreesOfFreedom)
  {
  }
  
  ResidualVariancePrior* ChiSquaredPrior::duplicate() const {
    ChiSquaredPrior* result = new ChiSquaredPrior();
    result->degreesOfFreedom = degreesOfFreedom;
    result->scale = scale;
    
    return result;
  }
  
  double ChiSquaredPrior::drawFromPosterior(const BARTFit& fit, size_t chainNum,
                                            const double* y,
                                            const double* y_hat) const
  {
    const Data& data(fit.data);
    
    double sumOfSquaredResiduals;
    if (fit.data.weights != NULL) {
      sumOfSquaredResiduals =
        misc_htm_computeWeightedSumOfSquaredResiduals(fit.threadManager, fit.chainScratch[chainNum].taskId,
                                                      y, data.numObservations, data.weights, y_hat);
    } else {
      sumOfSquaredResiduals =
        misc_htm_computeSumOfSquaredResiduals(fit.threadManager, fit.chainScratch[chainNum].taskId,
                                              y, data.numObservations, y_hat);
    }
    
    double posteriorDegreesOfFreedom = degreesOfFreedom + static_cast<double>(data.numObservations);
    
    double posteriorScale = degreesOfFreedom * scale + sumOfSquaredResiduals;
    
    return posteriorScale / ext_rng_simulateChiSquared(fit.state[chainNum].rng, posteriorDegreesOfFreedom);
  }
  
  void ChiSquaredPrior::print(const BARTFit& fit) const
  {
    ext_printf("\tdegrees of freedom in sigma prior: %f\n", degreesOfFreedom);
    double quantile = 1.0 - ext_percentileOfChiSquared(scale * degreesOfFreedom / fit.state[0].sigma / fit.state[0].sigma, degreesOfFreedom);
    ext_printf("\tquantile in sigma prior: %f\n", quantile);
    ext_printf("\tscale in sigma prior: %f\n", scale);
  }
  
  ResidualVariancePrior* FixedPrior::duplicate() const {
    return new FixedPrior(value);
  }

  
  void FixedPrior::print(const BARTFit&) const
  {
    ext_printf("\tresidual variance prior fixed to %f\n", value);
  }
  
  void ChiHyperprior::print(const BARTFit&) const
  {
    ext_printf("\tprior on k: chi with %f degrees of freedom and %f scale\n", degreesOfFreedom, scale);
  }
  
  double ChiHyperprior::drawFromPosterior(const BARTFit& fit, size_t chainNum) const
  {
    const Control& control(fit.control);
    const Data& data(fit.data);
    const Model& model(fit.model);
    
    const State& state(fit.state[chainNum]);
    
    double s_sq = 0.0;
    double totalNumBottomNodes = 0.0;
    for (size_t treeNum = 0; treeNum < control.numTrees; ++treeNum) {
      const double* treeFits = state.treeFits + treeNum * data.numObservations;
      
      size_t numBottomNodes;
      double* nodeParameters = state.trees[treeNum].recoverParametersFromFits(fit, treeFits, &numBottomNodes);
      totalNumBottomNodes += static_cast<double>(numBottomNodes);
      for (size_t i = 0; i < numBottomNodes; ++i)
        s_sq += nodeParameters[i] * nodeParameters[i];
      delete [] nodeParameters;
    }
    
    double numTrees = static_cast<double>(control.numTrees);
    
    double shape = 0.5 * (totalNumBottomNodes + 2.0 * degreesOfFreedom - 1.0);
    double rate = 0.5 * (numTrees * s_sq / (model.nodeScale * model.nodeScale));
    if (std::fabs(scale) <= std::numeric_limits<double>::max())
      rate += 0.5 / (scale * scale);
    
    return std::sqrt(ext_rng_simulateGamma(state.rng, shape, 1.0 / rate));
  }
}
