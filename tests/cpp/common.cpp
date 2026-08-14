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

// Tripwire for the comparison below and for the fuzz snapshot built on it: a
// new PERSISTED field must gain a comparison here, and a state field that
// nothing compares is a state field a rollback or restore gate cannot see. The
// size is the LP64 layout (std::vector 24 bytes, size_t 8); other data models
// are let through rather than guessed at. Honest, not airtight - a small field
// can hide in existing padding - which is why the table-driven coverage test
// beside the fuzz snapshot exists as well.
static_assert(sizeof(void*) != 8 || sizeof(ChainStateData) == 344,
              "ChainStateData gained or lost a field; add its comparison to "
              "statesAgree below and update this size");
static_assert(sizeof(void*) != 8 || sizeof(ForestStateData) == 160,
              "ForestStateData gained or lost a field; add its comparison to "
              "statesAgree below and update this size");

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
          xf.savedTreeMasks != yf.savedTreeMasks || xf.k != yf.k ||
          xf.leafScale != yf.leafScale)
        return false;
    }
    // the variance forest sits outside forests_, so its flat trees are a
    // sibling field rather than a ForestStateData member; a homoscedastic state
    // carries none on either side and agrees vacuously
    if (!sameFlatTrees(x.varianceTrees, y.varianceTrees)) return false;
    if (x.latents != y.latents || x.groupEffects != y.groupEffects ||
        x.dartProbabilities != y.dartProbabilities ||
        x.rngState != y.rngState)
      return false;
    // nu round-trips bitwise; a gaussian state carries NaN on both sides, which
    // an == would reject, so treat both-absent as agreement
    bool bothNu = std::isnan(x.residualDf) && std::isnan(y.residualDf);
    if (!bothNu && x.residualDf != y.residualDf) return false;
    // the ordinal cutpoint vector round-trips bitwise (restoreCutpoints is a
    // copy); a non-ordinal state carries an empty vector on both sides
    if (x.cutpoints != y.cutpoints) return false;
    // the nbinom dispersion r round-trips bitwise (restoreDispersion is a copy);
    // a non-count state carries NaN on both sides, which an == would reject
    bool bothDispersion =
      std::isnan(x.dispersion) && std::isnan(y.dispersion);
    if (!bothDispersion && x.dispersion != y.dispersion) return false;
    if (x.fitMin != y.fitMin || x.fitMax != y.fitMax ||
        x.groupTau != y.groupTau || x.dartAlpha != y.dartAlpha ||
        x.dartNumUpdatesSkipped != y.dartNumUpdatesSkipped)
      return false;
    if (x.hasBCF != y.hasBCF || x.amplitudeWidths != y.amplitudeWidths ||
        x.amplitudes != y.amplitudes ||
        x.amplitudeVariances != y.amplitudeVariances)
      return false;
    if (std::fabs(x.sigma - y.sigma) > 1e-9 * (1.0 + std::fabs(x.sigma)))
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
