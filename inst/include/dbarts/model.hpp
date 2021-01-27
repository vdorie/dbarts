#ifndef DBARTS_MODEL_HPP
#define DBARTS_MODEL_HPP

#include <cstddef>
#include <cmath> // sqrt
#include <dbarts/cstdint.hpp>

// can make these kinds of adjustments to trees during MCMC
#define DBARTS_BIRTH_OR_DEATH_PROBABILITY 0.5
#define DBARTS_SWAP_PROBABILITY           0.1
#define DBARTS_CHANGE_PROBABILITY         0.4
// conditional on being inside a birthOrDeath step
#define DBARTS_BIRTH_PROBABILITY          0.5

#define DBARTS_DEFAULT_NORMAL_PRIOR_K       2.0
#define DBARTS_DEFAULT_CHISQ_PRIOR_DF       3.0
#define DBARTS_DEFAULT_CHISQ_PRIOR_QUANTILE 0.9

#define DBARTS_DEFAULT_TREE_PRIOR_POWER 2.0
#define DBARTS_DEFAULT_TREE_PRIOR_BASE  0.95

struct ext_rng;

namespace dbarts {
  struct TreePrior;
  struct EndNodePrior;
  struct ResidualVariancePrior;
  struct EndNodeHyperprior;
  
  struct Model {
    double birthOrDeathProbability;
    double swapProbability;
    double changeProbability;

    double birthProbability;
    
    double nodeScale; // originally 3 for binary response, 0.5 for continuous
    
    TreePrior* treePrior;
    EndNodePrior* muPrior;
    ResidualVariancePrior* sigmaSqPrior;
    EndNodeHyperprior* kPrior;
    
    Model() : 
      birthOrDeathProbability(DBARTS_BIRTH_OR_DEATH_PROBABILITY),
      swapProbability(DBARTS_SWAP_PROBABILITY),
      changeProbability(DBARTS_CHANGE_PROBABILITY),
      birthProbability(DBARTS_BIRTH_PROBABILITY),
      nodeScale(-1.0),
      treePrior(NULL), muPrior(NULL), sigmaSqPrior(NULL), kPrior(NULL)
    {
    }
    Model(bool responseIsBinary) : 
      birthOrDeathProbability(DBARTS_BIRTH_OR_DEATH_PROBABILITY),
      swapProbability(DBARTS_SWAP_PROBABILITY),
      changeProbability(DBARTS_CHANGE_PROBABILITY),
      birthProbability(DBARTS_BIRTH_PROBABILITY),
      nodeScale(responseIsBinary ? 3.0 : 0.5),
      treePrior(NULL), muPrior(NULL), sigmaSqPrior(NULL), kPrior(NULL)
    {
    }
  };
  
  struct BARTFit;
  struct Control;
  struct Node;
  struct Tree;
  struct Rule;
  
  struct TreePrior {
    virtual double computeGrowthProbability(const BARTFit& fit, const Node& node) const = 0;
    virtual double computeTreeLogProbability(const BARTFit& fit, const Tree& tree) const = 0;
    
    virtual double computeSplitVariableLogProbability(const BARTFit& fit, const Node& node) const = 0;
    virtual double computeRuleForVariableLogProbability(const BARTFit& fit, const Node& node) const = 0;
    

    virtual Rule drawRuleAndVariable(const BARTFit& fit, ext_rng* rng, const Node& node, bool* exhaustedLeftSplits, bool* exhaustedRightSplits) const = 0;
    virtual std::int32_t drawSplitVariable(const BARTFit& fit, ext_rng* rng, const Node& node) const = 0;
    virtual Rule drawRuleForVariable(const BARTFit& fit, ext_rng* rng, const Node& node, std::int32_t variableIndex, bool* exhaustedLeftSplits, bool* exhaustedRightSplits) const = 0;
    
    virtual ~TreePrior() { }
  };
  
  struct EndNodePrior {
    virtual double computeLogIntegratedLikelihood(const BARTFit& fit, std::size_t chainNum, double k, const Node& node, const double* y, double residualVariance) const = 0;
    virtual double drawFromPosterior(ext_rng* rng, double k, double ybar, double numEffectiveObservations, double residualVariance) const = 0;
    virtual double drawFromPrior(ext_rng* rng, double k) const = 0;
    
    virtual ~EndNodePrior() { }
  };
  
  struct EndNodeHyperprior {
    bool isFixed;
    
    explicit EndNodeHyperprior(bool isFixed) : isFixed(isFixed) {}
    virtual ~EndNodeHyperprior() { }
    
    virtual double drawFromPosterior(const BARTFit& fit, std::size_t chainNum) const = 0;
    virtual void print(const BARTFit& fit) const = 0;
    
  };
  
  // the virtual scale accessors are for the conditional bart, which can have its data rescaled
  // if your prior doesn't use them, ignore them
  struct ResidualVariancePrior {
    bool isFixed;
    
    explicit ResidualVariancePrior(bool isFixed) : isFixed(isFixed) {}
    virtual ~ResidualVariancePrior() { }
    
    virtual double drawFromPosterior(const BARTFit& fit, std::size_t chainNum,
                                     const double* y,
                                     const double* y_hat) const = 0;
    
    virtual double getScale() const = 0;
    virtual void setScale(double scale) = 0;
    virtual ResidualVariancePrior* duplicate() const = 0;
    
    virtual void print(const BARTFit& fit) const = 0;
  };
  
  // for lack of a better name, calling it the Chipman, George, and McCullough prior
  // Pr(node splits) = base / (1 + depth)^power
  
  struct CGMPrior : TreePrior {
    double base;
    double power;
    
    CGMPrior() { }
    CGMPrior(double base, double power) : base(base), power(power) { }
    virtual ~CGMPrior() { }
    
    virtual double computeGrowthProbability(const BARTFit& fit, const Node& node) const;
    virtual double computeTreeLogProbability(const BARTFit& fit, const Tree& tree) const;
    
    virtual double computeSplitVariableLogProbability(const BARTFit& fit, const Node& node) const;
    virtual double computeRuleForVariableLogProbability(const BARTFit& fit, const Node& node) const;
    
    virtual Rule drawRuleAndVariable(const BARTFit& fit, ext_rng* rng, const Node& node, bool* exhaustedLeftSplits, bool* exhaustedRightSplits) const;
    virtual std::int32_t drawSplitVariable(const BARTFit& fit, ext_rng* rng, const Node& node) const;
    virtual Rule drawRuleForVariable(const BARTFit& fit, ext_rng* rng, const Node& node, std::int32_t variableIndex, bool* exhaustedLeftSplits, bool* exhaustedRightSplits) const;
  };
  
  // nodeMu ~ normal(0, 1 / precision)
  struct NormalPrior : EndNodePrior {
    double scale;
    
    NormalPrior() : scale(1.0) { }
    NormalPrior(const Control& control, const Model& model);
    virtual ~NormalPrior() { }
    
    virtual double computeLogIntegratedLikelihood(const BARTFit& fit, std::size_t chainNum, double k, const Node& node, const double* y, double residualVariance) const;
    virtual double drawFromPosterior(ext_rng* rng, double k, double ybar, double numEffectiveObservations, double residualVariance) const;
    virtual double drawFromPrior(ext_rng* rng, double k) const;
    virtual void setScale(double newScale) { scale = newScale; }
    virtual void setScale(const Control& control, const Model& model);
    virtual double getScale() const { return scale; }
  };
  
  /* struct CauchyHyperprior : EndNodeHyperprior {
    double scale;
    
    CauchyHyperprior() : scale(1.0) { }
    CauchyHyperprior(double scale) : location(scale) { }
    
    virtual double drawFromPosterior(const BARTFit& fit, std::size_t chainNum) const;
  }; */
  
  struct ChiHyperprior : EndNodeHyperprior {
    double degreesOfFreedom;
    double scale;
    
    ChiHyperprior() : EndNodeHyperprior(false), degreesOfFreedom(1.25), scale(1.0) { }
    ChiHyperprior(double degreesOfFreedom, double scale) : EndNodeHyperprior(false), degreesOfFreedom(degreesOfFreedom), scale(scale) { }
    virtual ~ChiHyperprior() { }
    
    virtual void print(const BARTFit& fit) const;
    virtual double drawFromPosterior(const BARTFit& fit, std::size_t chainNum) const;
  };
  
  struct FixedHyperprior : EndNodeHyperprior {
    double k;
    
    FixedHyperprior() : EndNodeHyperprior(true), k(2.0) { }
    FixedHyperprior(double k) : EndNodeHyperprior(true), k(k) { }
    virtual ~FixedHyperprior() { }
    
    double getK() const { return k; }
    void setK(double newK) { k = newK; }
    
    virtual void print(const BARTFit& fit) const;
    virtual double drawFromPosterior(const BARTFit&, std::size_t) const { return k; }
  };
  
  // sigmaSq ~ chisq(df, scale)
  struct ChiSquaredPrior : ResidualVariancePrior {
    double degreesOfFreedom;
    double scale;
    
    ChiSquaredPrior() : ResidualVariancePrior(false),
      degreesOfFreedom(DBARTS_DEFAULT_CHISQ_PRIOR_DF),
      scale(1.0) { }
    ChiSquaredPrior(double degreesOfFreedom, double quantile);
    virtual ~ChiSquaredPrior() { }
    virtual ResidualVariancePrior* duplicate() const;
    
    virtual double getScale() const { return scale; }
    virtual void setScale(double newScale) { scale = newScale; }
    
    virtual void print(const BARTFit& fit) const;
    
    virtual double drawFromPosterior(const BARTFit& fit, std::size_t chainNum,
                                     const double* y,
                                     const double* y_hat) const;
  };
  struct FixedPrior : ResidualVariancePrior {
    double value;
    
    FixedPrior() : ResidualVariancePrior(true), value(1.0) { }
    FixedPrior(double value) : ResidualVariancePrior(true), value(value) { }
    virtual ResidualVariancePrior* duplicate() const;
    
    virtual double drawFromPosterior(const BARTFit&, std::size_t,
                                     const double*,
                                     const double*) const {
      return value;
    }
    
    virtual double getScale() const { return value; }
    virtual void setScale(double scale) { value = scale; }
    
    virtual void print(const BARTFit& fit) const;
    
    virtual ~FixedPrior() { }
  };
} // namespace dbarts

#endif // DBARTS_MODEL_HPP

