#include "common.hpp"

#include <cstdarg>

namespace {

std::string* printSink = nullptr;

}  // namespace

extern "C" void Rprintf(const char* format, ...) {
  if (printSink == nullptr) return;
  char line[4096];
  va_list args;
  va_start(args, format);
  int written = std::vsnprintf(line, sizeof(line), format, args);
  va_end(args);
  if (written <= 0) return;
  size_t length = static_cast<size_t>(written);
  printSink->append(line, length < sizeof(line) ? length : sizeof(line) - 1);
}

void beginPrintCapture(std::string& sink) {
  sink.clear();
  printSink = &sink;
}

void endPrintCapture() { printSink = nullptr; }

bool sameFlatTrees(const std::vector<std::vector<FlatNode>>& a,
                   const std::vector<std::vector<FlatNode>>& b) {
  if (a.size() != b.size()) return false;
  for (size_t t = 0; t < a.size(); ++t) {
    if (a[t].size() != b[t].size()) return false;
    for (size_t i = 0; i < a[t].size(); ++i) {
      const FlatNode& x(a[t][i]);
      const FlatNode& y(b[t][i]);
      if (x.variable != y.variable || x.numMaskWords != y.numMaskWords ||
          x.flags != y.flags)
        return false;
      if (x.variable == invalidVariable) {
        if (std::fabs(x.value - y.value) > 1e-9 * (1.0 + std::fabs(x.value)))
          return false;
      } else if (x.mask != y.mask) {
        return false;
      }
    }
  }
  return true;
}

bool sendFirstOrdinalRuleMissingRight(
    std::vector<std::vector<FlatNode>>& trees) {
  for (std::vector<FlatNode>& tree : trees)
    for (FlatNode& node : tree)
      if (node.variable != invalidVariable &&
          flatKindOf(node) == FlatKind::ordinal) {
        node.flags |= flatMissingGoesRight;
        return true;
      }
  return false;
}

// Tripwire for the comparison below and for the fuzz snapshot built on it: a
// new PERSISTED field must gain a comparison here, and a state field that
// nothing compares is a state field a rollback or restore gate cannot see. The
// size is the LP64 layout (std::vector 24 bytes, size_t 8); other data models
// are let through rather than guessed at. Honest, not airtight - a small field
// can hide in existing padding - which is why the table-driven coverage test
// beside the fuzz snapshot exists as well.
static_assert(sizeof(void*) != 8 || sizeof(ChainStateData) == 384,
              "ChainStateData gained or lost a field; add its comparison to "
              "statesAgree below and update this size");
static_assert(sizeof(void*) != 8 || sizeof(ForestStateData) == 224,
              "ForestStateData gained or lost a field; add its comparison to "
              "statesAgree below and update this size");

// a scalar a state holds only where it is drawn is NaN, absent, otherwise;
// two absent values agree
static bool sameScalar(double x, double y) {
  return (std::isnan(x) && std::isnan(y)) || x == y;
}

// a leaf covariate's scale is absent where its column had no spread
static bool sameScales(const std::vector<double>& x,
                       const std::vector<double>& y) {
  if (x.size() != y.size()) return false;
  for (size_t j = 0; j < x.size(); ++j)
    if (!sameScalar(x[j], y[j])) return false;
  return true;
}

bool statesAgree(const SamplerStateData& a, const SamplerStateData& b) {
  if (a.chains.size() != b.chains.size()) return false;
  for (size_t c = 0; c < a.chains.size(); ++c) {
    const ChainStateData& x(a.chains[c]);
    const ChainStateData& y(b.chains[c]);
    if (x.forests.size() != y.forests.size()) return false;
    for (size_t f = 0; f < x.forests.size(); ++f) {
      const ForestStateData& xf(x.forests[f]);
      const ForestStateData& yf(y.forests[f]);
      if (!sameFlatTrees(xf.trees, yf.trees) ||
          !sameFlatTrees(xf.savedTrees, yf.savedTrees))
        return false;
      if (xf.treeParams != yf.treeParams ||
          xf.savedTreeParams != yf.savedTreeParams ||
          xf.treeMasks != yf.treeMasks ||
          xf.savedTreeMasks != yf.savedTreeMasks || !sameScalar(xf.k, yf.k) ||
          xf.leafCovariateCenters != yf.leafCovariateCenters ||
          !sameScales(xf.leafCovariateScales, yf.leafCovariateScales) ||
          xf.leafLengthscales != yf.leafLengthscales)
        return false;
    }
    // the variance forest sits outside forests_, so its flat trees are
    // sibling fields rather than ForestStateData members; a homoscedastic state
    // carries none on either side and agrees vacuously. The SAVED buffer is
    // compared too: it is the only record of the kept samples' scale surface,
    // and a restore that rebuilt the live trees alone would predict off the
    // destination's identity fill unnoticed.
    if (!sameFlatTrees(x.varianceTrees, y.varianceTrees) ||
        !sameFlatTrees(x.savedVarianceTrees, y.savedVarianceTrees))
      return false;
    if (x.varianceTreeMasks != y.varianceTreeMasks ||
        x.savedVarianceTreeMasks != y.savedVarianceTreeMasks)
      return false;
    if (x.latents != y.latents ||
        x.dartProbabilities != y.dartProbabilities ||
        x.rngState != y.rngState)
      return false;
    // nu round-trips bitwise where drawn and is absent on both sides otherwise
    if (!sameScalar(x.residualDf, y.residualDf)) return false;
    // the ordinal threshold vector round-trips bitwise
    // (restoreOrdinalThresholds is a copy); a non-ordinal state carries an
    // empty vector on both sides
    if (x.ordinalThresholds != y.ordinalThresholds) return false;
    // the nbinom shape r round-trips bitwise (restoreShape is a copy) where
    // drawn
    if (!sameScalar(x.shape, y.shape)) return false;
    if (x.fitMin != y.fitMin || x.fitMax != y.fitMax ||
        !sameScalar(x.dartAlpha, y.dartAlpha) ||
        x.dartNumUpdatesSkipped != y.dartNumUpdatesSkipped)
      return false;
    if (x.hasAmplitudes != y.hasAmplitudes ||
        x.amplitudeWidths != y.amplitudeWidths ||
        x.amplitudes != y.amplitudes ||
        x.amplitudeVariances != y.amplitudeVariances)
      return false;
    if (std::isnan(x.sigma) != std::isnan(y.sigma) ||
        std::fabs(x.sigma - y.sigma) > 1e-9 * (1.0 + std::fabs(x.sigma)))
      return false;
  }
  return true;
}

// A burned-in sampler for mutation tests: strong signal in both columns so
// trees certainly split.
std::unique_ptr<ConstantLeafSampler> makeBurnedInSampler(
  std::vector<double>& x, std::vector<double>& y, size_t n, ext_rng* rng) {
  SamplerOptions options;
  options.numTrees = 25;
  auto sampler = std::make_unique<ConstantLeafSampler>(
    x.data(), y.data(), n, size_t(2), nullptr, nullptr, ResponseFamily::gaussian, 1.0, 3.0,
    0.37804942330213542, options, &rng);
  Results empty;
  sampler->run(100, 0, empty);
  return sampler;
}

void makeMutationData(std::vector<double>& x, std::vector<double>& y,
                       size_t n) {
  x.resize(n * 2);
  y.resize(n);
  for (double& v : x) v = runif01();
  for (size_t i = 0; i < n; ++i) {
    double u1 = runif01(), u2 = runif01();
    double normal = std::sqrt(-2.0 * std::log(u1)) *
                    std::cos(6.283185307179586 * u2);
    y[i] = 4.0 * (x[i] - 0.5) + 2.0 * x[i + n] + 0.2 * normal;
  }
}

std::vector<double> normalizedFromLogWeights(
    const std::vector<double>& logWeights) {
  double maxLogWeight = *std::max_element(logWeights.begin(), logWeights.end());
  std::vector<double> probabilities(logWeights.size());
  double sum = 0.0;
  for (size_t i = 0; i < logWeights.size(); ++i) {
    probabilities[i] = std::exp(logWeights[i] - maxLogWeight);
    sum += probabilities[i];
  }
  for (double& probability : probabilities) probability /= sum;
  return probabilities;
}

// Pearson goodness of fit against a fully specified law: df = cells - 1
double chiSquareStatistic(const std::vector<double>& counts,
                          const std::vector<double>& probabilities,
                          double numDraws) {
  double statistic = 0.0;
  for (size_t i = 0; i < counts.size(); ++i) {
    double expected = numDraws * probabilities[i];
    double deviation = counts[i] - expected;
    statistic += deviation * deviation / expected;
  }
  return statistic;
}

// Regularized upper incomplete gamma Q(a, x): the series for P below the
// crossover, the Lentz continued fraction for Q above it. Coded here because
// libR's own pchisq silently returns zero without an initialized R runtime,
// which this standalone host is not; agrees with R's to 6 figures.
static double upperIncompleteGamma(double a, double x) {
  double logGammaA = std::lgamma(a);
  if (x < a + 1.0) {
    double term = 1.0 / a, sum = term;
    for (int i = 1; i < 1000; ++i) {
      term *= x / (a + i);
      sum += term;
      if (std::fabs(term) < std::fabs(sum) * 1e-16) break;
    }
    return 1.0 - sum * std::exp(-x + a * std::log(x) - logGammaA);
  }
  const double tiny = 1e-300;
  double b = x + 1.0 - a, c = 1.0 / tiny, d = 1.0 / b, h = d;
  for (int i = 1; i < 1000; ++i) {
    double an = -i * (i - a);
    b += 2.0;
    d = an * d + b;
    if (std::fabs(d) < tiny) d = tiny;
    c = b + an / c;
    if (std::fabs(c) < tiny) c = tiny;
    d = 1.0 / d;
    double delta = d * c;
    h *= delta;
    if (std::fabs(delta - 1.0) < 1e-16) break;
  }
  return h * std::exp(-x + a * std::log(x) - logGammaA);
}

double chiSquareUpperTail(double statistic, double df) {
  return upperIncompleteGamma(0.5 * df, 0.5 * statistic);
}
