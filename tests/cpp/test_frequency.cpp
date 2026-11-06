// Frequency weights on the gaussian-type families: a row of weight w read as
// w identical rows. The oracle throughout is the fit to the data with every
// row written out w times as adjacent copies, run under unit weights. A
// gaussian sweep's random numbers do not depend on how many times a row
// appears, so the gaussian, heteroscedastic, linear-leaf and two-forest fits
// must agree with that expansion draw for draw, to rounding; Student-t's
// latent is a sum over the copies and is held to it in distribution and, for
// its df step, to the exact grid probabilities. The suite takes its own
// generators and restores rngState.

#include "common.hpp"

#include <array>
#include <numbers>

namespace {

constexpr double sigmaDf = 3.0;
constexpr double sigmaRawScale = 0.37804942330213542;

ext_rng* seededRng(unsigned int seed) {
  ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  ext_rng_setSeed(rng, seed);
  return rng;
}

/// n rows over two columns, a smooth signal plus noise, and whole weights 1
/// to 4; column-major x.
struct Fixture {
  std::size_t n = 0, p = 2;
  std::vector<double> x, y, w;

  explicit Fixture(std::size_t n_) : n(n_) {
    x.resize(n * p);
    y.resize(n);
    w.resize(n);
    for (double& v : x) v = runif01();
    for (std::size_t i = 0; i < n; ++i) {
      y[i] = std::sin(6.0 * x[i]) + 2.0 * x[n + i] + 0.4 * (runif01() - 0.5);
      w[i] = 1.0 + std::floor(4.0 * runif01());
    }
  }
};

/// The rows written out: row i of \p columns (column-major over \p n rows,
/// \p p columns) repeated \p w[i] times, adjacently. rowOf maps each written
/// row to its source row.
struct Expansion {
  std::size_t n = 0;
  std::vector<double> x, y, extra;
  std::vector<std::size_t> rowOf;

  Expansion(const Fixture& f, const std::vector<double>& weights,
            const std::vector<double>* column = nullptr) {
    for (std::size_t i = 0; i < f.n; ++i)
      for (int c = 0; c < static_cast<int>(weights[i]); ++c) rowOf.push_back(i);
    n = rowOf.size();
    x.resize(n * f.p);
    y.resize(n);
    for (std::size_t r = 0; r < n; ++r) {
      for (std::size_t j = 0; j < f.p; ++j)
        x[j * n + r] = f.x[j * f.n + rowOf[r]];
      y[r] = f.y[rowOf[r]];
    }
    if (column != nullptr) {
      extra.resize(n);
      for (std::size_t r = 0; r < n; ++r) extra[r] = (*column)[rowOf[r]];
    }
  }
};

/// A run's reported draws.
struct Draws {
  std::vector<double> sigma, fits, logLikelihood, variance, df;

  Draws(std::size_t n, std::size_t numSamples)
    : sigma(numSamples), fits(n * numSamples), logLikelihood(n * numSamples),
      variance(n * numSamples), df(numSamples) {}

  Results results(bool heteroscedastic = false, bool student = false) {
    Results r;
    r.sigma = sigma.data();
    r.trainingFits = fits.data();
    r.logLikelihood = logLikelihood.data();
    if (heteroscedastic) r.varianceFits = variance.data();
    if (student) r.residualDf = df.data();
    return r;
  }
};

double relativeGap(double a, double b) {
  return std::fabs(a - b) / std::max(1.0, std::fabs(b));
}

/// The largest gap between a weighted run and its expansion over sigma and
/// every training fit, each expanded row read against its source row.
double drawGap(const Draws& weighted, const Draws& expanded,
               const Expansion& e, std::size_t n, std::size_t numSamples,
               bool variance = false) {
  double gap = 0.0;
  for (std::size_t s = 0; s < numSamples; ++s) {
    gap = std::max(gap, relativeGap(weighted.sigma[s], expanded.sigma[s]));
    for (std::size_t r = 0; r < e.n; ++r) {
      std::size_t i = e.rowOf[r];
      gap = std::max(gap, relativeGap(weighted.fits[s * n + i],
                                      expanded.fits[s * e.n + r]));
      if (variance)
        gap = std::max(gap, relativeGap(weighted.variance[s * n + i],
                                        expanded.variance[s * e.n + r]));
    }
  }
  return gap;
}

bool allFinite(const std::vector<double>& values) {
  for (double v : values)
    if (!std::isfinite(v)) return false;
  return true;
}

template <typename S>
std::unique_ptr<S> makeSampler(const double* x, const double* y,
                               std::size_t n, std::size_t p,
                               const double* weights,
                               const SamplerOptions& options, ext_rng** rng) {
  return std::make_unique<S>(x, y, n, p, weights, nullptr,
                             ResponseFamily::gaussian, 1.0, sigmaDf,
                             sigmaRawScale, options, rng);
}

SamplerOptions frequencyOptions(SamplerOptions options) {
  options.weightsReading = WeightsReading::frequency;
  return options;
}

/// A frequency run against its expansion, draw for draw, for the leaf model
/// S and the options given.
template <typename S>
void checkAgainstExpansion(const SamplerOptions& options, unsigned int seed,
                           const char* what, bool heteroscedastic = false) {
  const std::size_t n = 40, numSamples = 300;
  Fixture f(n);
  Expansion e(f, f.w);
  ext_rng* rngW = seededRng(seed);
  ext_rng* rngE = seededRng(seed);
  auto weighted = makeSampler<S>(f.x.data(), f.y.data(), n, f.p, f.w.data(),
                                 frequencyOptions(options), &rngW);
  auto expanded = makeSampler<S>(e.x.data(), e.y.data(), e.n, f.p, nullptr,
                                 options, &rngE);
  Draws dw(n, numSamples), de(e.n, numSamples);
  Results rw = dw.results(heteroscedastic), re = de.results(heteroscedastic);
  bool ran = !weighted->run(0, numSamples, rw) && !expanded->run(0, numSamples, re);
  double gap = drawGap(dw, de, e, n, numSamples, heteroscedastic);
  printf("  %s: largest gap %.3g over %zu sweeps\n", what, gap, numSamples);
  check(ran && gap < 1e-12, what);
  ext_rng_destroy(rngW);
  ext_rng_destroy(rngE);
}

void testGaussianAgainstExpansion() {
  SamplerOptions options;
  options.numTrees = 20;
  checkAgainstExpansion<ConstantLeafSampler>(
    options, 20261106u,
    "frequency gaussian equals its expansion draw for draw (sigma, fits)");
  options.useQuantiles = true;
  options.maxNumCuts = 8;
  checkAgainstExpansion<ConstantLeafSampler>(
    options, 20261107u,
    "frequency gaussian under quantile cuts equals its expansion");
}

void testHeteroscedasticAgainstExpansion() {
  SamplerOptions options;
  options.numTrees = 20;
  options.numVarianceTrees = 5;
  checkAgainstExpansion<ConstantLeafSampler>(
    options, 20261108u,
    "frequency heteroscedastic equals its expansion draw for draw "
    "(fits, variance surface)",
    true);
}

void testLinearLeafAgainstExpansion() {
  const std::size_t columns[] = {0};
  SamplerOptions options;
  options.numTrees = 10;
  options.leafCovariateColumns = columns;
  options.numLeafCovariates = 1;
  checkAgainstExpansion<Sampler<LinearGaussianLeaf>>(
    options, 20261109u,
    "frequency linear leaf equals its expansion draw for draw (the "
    "covariate's weighted center and scale)");
}

/// Two forests, the second multiplying a numeric column that states no sd, so
/// its default sd is the column's scale over the rows, and the calibration
/// anchor the response's sd: both read off the rows written out.
void testTwoForestAgainstExpansion() {
  const std::size_t n = 40, numSamples = 300;
  Fixture f(n);
  std::vector<double> dose(n);
  for (std::size_t i = 0; i < n; ++i) dose[i] = 0.5 + runif01();
  Expansion e(f, f.w, &dose);
  auto spec = [](const double* basis) {
    AmplitudeSpec s;
    s.forests.resize(2);
    s.forests[0].forest.numTrees = 15;
    s.forests[1].forest.numTrees = 10;
    s.forests[1].kind = ForestKind::numeric;
    s.forests[1].basis = basis;
    return s;
  };
  SamplerOptions options;
  ext_rng* rngW = seededRng(20261110u);
  ext_rng* rngE = seededRng(20261110u);
  ConstantLeafSampler weighted(f.x.data(), f.y.data(), n, f.p, f.w.data(),
                               nullptr, 1.0, sigmaDf, sigmaRawScale,
                               frequencyOptions(options), spec(dose.data()),
                               &rngW);
  ConstantLeafSampler expanded(e.x.data(), e.y.data(), e.n, f.p, nullptr,
                               nullptr, 1.0, sigmaDf, sigmaRawScale, options,
                               spec(e.extra.data()), &rngE);
  ForestCalibration cw = weighted.forestCalibration(0, 1),
                    ce = expanded.forestCalibration(0, 1);
  check(relativeGap(cw.basisScale, ce.basisScale) < 1e-14 &&
          relativeGap(cw.mapAnchor, ce.mapAnchor) < 1e-14,
        "frequency two-forest: the basis scale and the anchor are the "
        "expansion's");
  Draws dw(n, numSamples), de(e.n, numSamples);
  Results rw = dw.results(), re = de.results();
  rw.logLikelihood = re.logLikelihood = nullptr;
  bool ran = !weighted.run(0, numSamples, rw) && !expanded.run(0, numSamples, re);
  double gap = drawGap(dw, de, e, n, numSamples);
  printf("  two-forest: largest gap %.3g over %zu sweeps\n", gap, numSamples);
  check(ran && gap < 1e-12,
        "frequency two-forest equals its expansion draw for draw");
  ext_rng_destroy(rngW);
  ext_rng_destroy(rngE);
}

/// At all-ones weights the two readings are one model and one expression:
/// bitwise, on gaussian, Student-t and heteroscedastic fits, and on a
/// gaussian fit at 0/1 weights.
void testUnitWeightIdentity() {
  const std::size_t n = 40, numSamples = 60;
  Fixture f(n);
  std::vector<double> ones(n, 1.0), zeroOne(n, 1.0);
  for (std::size_t i = 0; i < n; i += 5) zeroOne[i] = 0.0;
  struct Arm {
    const char* what;
    const std::vector<double>* weights;
    double residualDf;
    std::size_t varianceTrees;
  };
  const Arm arms[] = {
    {"gaussian, all ones", &ones, std::nan(""), 0},
    {"gaussian, 0/1", &zeroOne, std::nan(""), 0},
    {"Student-t, all ones", &ones, 0.0, 0},
    {"heteroscedastic, all ones", &ones, std::nan(""), 4}};
  for (const Arm& arm : arms) {
    SamplerOptions options;
    options.numTrees = 15;
    options.residualDf = arm.residualDf;
    options.numVarianceTrees = arm.varianceTrees;
    bool student = std::isfinite(arm.residualDf);
    bool hetero = arm.varianceTrees > 0;
    ext_rng* rngP = seededRng(4401u);
    ext_rng* rngF = seededRng(4401u);
    auto precision = makeSampler<ConstantLeafSampler>(
      f.x.data(), f.y.data(), n, f.p, arm.weights->data(), options, &rngP);
    auto frequency = makeSampler<ConstantLeafSampler>(
      f.x.data(), f.y.data(), n, f.p, arm.weights->data(),
      frequencyOptions(options), &rngF);
    Draws dp(n, numSamples), df(n, numSamples);
    Results rp = dp.results(hetero, student), rf = df.results(hetero, student);
    rp.logLikelihood = rf.logLikelihood = nullptr;
    bool ran = !precision->run(10, numSamples, rp) &&
               !frequency->run(10, numSamples, rf);
    bool same = dp.sigma == df.sigma && dp.fits == df.fits &&
                (!student || dp.df == df.df) &&
                (!hetero || dp.variance == df.variance);
    std::string label = std::string("frequency at unit weights is precision "
                                    "bitwise: ") + arm.what;
    check(ran && same, label.c_str());
    ext_rng_destroy(rngP);
    ext_rng_destroy(rngF);
  }
}

/// A row of frequency weight 0 is out of the likelihood and keeps the run
/// finite on each family, its log-likelihood entry NaN and every other one
/// finite.
void testZeroWeights() {
  const std::size_t n = 40, numSamples = 40;
  Fixture f(n);
  std::vector<double> w(f.w);
  for (std::size_t i = 3; i < n; i += 7) w[i] = 0.0;
  w[10] = 2.5;  // a real count is taken
  struct Arm { const char* what; double residualDf; std::size_t varianceTrees; };
  const Arm arms[] = {{"gaussian", std::nan(""), 0},
                      {"Student-t, drawn df", 0.0, 0},
                      {"Student-t, fixed df", 5.0, 0},
                      {"heteroscedastic", std::nan(""), 4}};
  for (const Arm& arm : arms) {
    SamplerOptions options;
    options.numTrees = 15;
    options.residualDf = arm.residualDf;
    options.numVarianceTrees = arm.varianceTrees;
    bool student = std::isfinite(arm.residualDf);
    bool hetero = arm.varianceTrees > 0;
    ext_rng* rng = seededRng(7301u);
    auto sampler = makeSampler<ConstantLeafSampler>(
      f.x.data(), f.y.data(), n, f.p, w.data(), frequencyOptions(options),
      &rng);
    Draws d(n, numSamples);
    Results r = d.results(hetero, student);
    bool ran = !sampler->run(20, numSamples, r);
    bool finite = allFinite(d.sigma) && allFinite(d.fits) &&
                  (!student || allFinite(d.df)) &&
                  (!hetero || allFinite(d.variance));
    for (double s : d.sigma) finite = finite && s > 0.0;
    if (student)
      for (double nu : d.df)
        finite = finite &&
                 (nu == 5.0 || std::find(std::begin(ResidualDfPrior::grid),
                                         std::end(ResidualDfPrior::grid),
                                         nu) != std::end(ResidualDfPrior::grid));
    bool logLikelihood = true;
    for (std::size_t s = 0; s < numSamples; ++s)
      for (std::size_t i = 0; i < n; ++i) {
        double v = d.logLikelihood[s * n + i];
        logLikelihood = logLikelihood &&
                        (w[i] == 0.0 ? std::isnan(v) : std::isfinite(v));
      }
    std::string label =
      std::string("frequency weight 0 runs finite, its log-likelihood NaN: ") +
      arm.what;
    check(ran && finite && logLikelihood, label.c_str());
    ext_rng_destroy(rng);
  }
}

/// A frequency sampler whose rows a mask takes out is the one created with
/// those rows' weights at 0, draw for draw: the count, the composite and, on
/// Student-t, the df step's constant are recomputed at the mask change. The
/// rows masked hold neither end of the response, so both read one range.
void testMaskAgainstZeroWeights() {
  const std::size_t n = 40, numSamples = 80;
  Fixture f(n);
  std::size_t lowest = 0, highest = 0;
  for (std::size_t i = 1; i < n; ++i) {
    if (f.y[i] < f.y[lowest]) lowest = i;
    if (f.y[i] > f.y[highest]) highest = i;
  }
  std::vector<double> mask(n, 1.0), zeroed(f.w);
  for (std::size_t i = 2; i < n; i += 6) {
    if (i == lowest || i == highest) continue;
    mask[i] = 0.0;
    zeroed[i] = 0.0;
  }
  for (double residualDf : {std::nan(""), 0.0}) {
    bool student = std::isfinite(residualDf);
    SamplerOptions options = frequencyOptions(SamplerOptions());
    options.numTrees = 15;
    options.residualDf = residualDf;
    ext_rng* rngM = seededRng(3307u);
    ext_rng* rngZ = seededRng(3307u);
    auto masked = makeSampler<ConstantLeafSampler>(
      f.x.data(), f.y.data(), n, f.p, f.w.data(), options, &rngM);
    auto zero = makeSampler<ConstantLeafSampler>(
      f.x.data(), f.y.data(), n, f.p, zeroed.data(), options, &rngZ);
    bool installed = masked->setActiveRows(mask.data());
    Draws dm(n, numSamples), dz(n, numSamples);
    Results rm = dm.results(false, student), rz = dz.results(false, student);
    bool ran = !masked->run(10, numSamples, rm) && !zero->run(10, numSamples, rz);
    check(installed && ran && dm.sigma == dz.sigma && dm.fits == dz.fits &&
            (!student || dm.df == dz.df),
          student ? "frequency Student-t masked by setActiveRows is the one "
                    "created at weight 0, draw for draw"
                  : "frequency gaussian masked by setActiveRows is the one "
                    "created at weight 0, draw for draw");
    ext_rng_destroy(rngM);
    ext_rng_destroy(rngZ);
  }
}

/// log t_nu(x), from its definition.
double logStudentT(double x, double nu) {
  return std::lgamma(0.5 * (nu + 1.0)) - std::lgamma(0.5 * nu) -
         0.5 * std::log(nu * std::numbers::pi) -
         0.5 * (nu + 1.0) * std::log1p(x * x / nu);
}

/// The df step on fixed residuals: alternating the latent draw and the grid
/// draw leaves nu at its exact conditional given the residuals, a row of
/// count m contributing m t densities, prior on the grid nu exp(-nu / 10).
/// Run through TResponse's own sweep at fixed fits and sigma, at whole and
/// real counts, a zero count included, and under precision at unit weights.
void testResidualDfGrid() {
  const std::size_t n = 30, numDraws = 60000;
  std::vector<double> y(n), counts(n), ones(n, 1.0), fits(n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    // heavy-tailed residuals, so the df is informed
    double u = runif01() - 0.5;
    y[i] = u * (1.0 + 4.0 * u * u * (i % 3 == 0 ? 6.0 : 1.0));
    counts[i] = i % 5 == 0 ? 0.4 + i / 10.0 : 1.0 + (i % 4);
  }
  counts[7] = 0.0;
  const double sigma = 0.15;
  struct Arm {
    const char* what;
    const double* weights;
    WeightsReading reading;
  };
  const Arm arms[] = {
    {"frequency, whole and real counts", counts.data(),
     WeightsReading::frequency},
    {"precision, unit weights", ones.data(), WeightsReading::precision}};
  for (const Arm& arm : arms) {
    TResponse response(y.data(), nullptr, arm.weights, n, 1.0, sigmaDf,
                       sigmaRawScale, 0.0, arm.reading);
    const double* z = response.workingResponse();
    std::array<double, ResidualDfPrior::gridSize> logExact{};
    double total = 0.0, largest = -HUGE_VAL;
    for (std::size_t k = 0; k < ResidualDfPrior::gridSize; ++k) {
      double nu = ResidualDfPrior::grid[k];
      double value = std::log(nu) - 0.1 * nu;
      for (std::size_t i = 0; i < n; ++i)
        value += arm.weights[i] *
                 (logStudentT(z[i] / sigma, nu) - std::log(sigma));
      logExact[k] = value;
      largest = std::max(largest, value);
    }
    std::array<double, ResidualDfPrior::gridSize> exact{}, observed{};
    for (std::size_t k = 0; k < ResidualDfPrior::gridSize; ++k)
      total += exact[k] = std::exp(logExact[k] - largest);
    for (double& v : exact) v /= total;
    ext_rng* rng = seededRng(9091u);
    for (std::size_t s = 0; s < 500; ++s)
      response.refreshLatents(rng, fits.data(), sigma);
    for (std::size_t s = 0; s < numDraws; ++s) {
      response.refreshLatents(rng, fits.data(), sigma);
      for (std::size_t k = 0; k < ResidualDfPrior::gridSize; ++k)
        if (response.residualDf() == ResidualDfPrior::grid[k])
          observed[k] += 1.0 / static_cast<double>(numDraws);
    }
    double gap = 0.0, spread = 0.0;
    for (std::size_t k = 0; k < ResidualDfPrior::gridSize; ++k) {
      gap = std::max(gap, std::fabs(observed[k] - exact[k]));
      spread = std::max(spread, exact[k]);
    }
    printf("  df grid, %s: largest gap %.4f, largest probability %.3f\n",
           arm.what, gap, spread);
    std::string label =
      std::string("Student-t df step on fixed residuals reaches the exact "
                  "grid probabilities: ") + arm.what;
    check(gap < 0.01 && spread < 0.5, label.c_str());
    ext_rng_destroy(rng);
  }
}

/// Batch means of \p values: the mean and its standard error over 20
/// batches.
std::pair<double, double> batchMean(const std::vector<double>& values) {
  const std::size_t numBatches = 20, size = values.size() / numBatches;
  double grand = 0.0, sumSquares = 0.0;
  std::vector<double> means(numBatches, 0.0);
  for (std::size_t b = 0; b < numBatches; ++b) {
    for (std::size_t j = 0; j < size; ++j) means[b] += values[b * size + j];
    means[b] /= static_cast<double>(size);
    grand += means[b];
  }
  grand /= static_cast<double>(numBatches);
  for (double m : means) sumSquares += (m - grand) * (m - grand);
  return {grand, std::sqrt(sumSquares / (numBatches - 1.0) / numBatches)};
}

/// Student-t against its expansion in distribution: posterior means of
/// sigma, nu and two rows' fits within four standard errors.
void testStudentAgainstExpansion() {
  const std::size_t n = 40, numSamples = 4000;
  Fixture f(n);
  for (std::size_t i = 0; i < n; i += 9) f.y[i] += 1.5;  // outliers
  Expansion e(f, f.w);
  SamplerOptions options;
  options.numTrees = 15;
  options.residualDf = 0.0;
  ext_rng* rngW = seededRng(6100u);
  ext_rng* rngE = seededRng(6200u);
  auto weighted = makeSampler<ConstantLeafSampler>(
    f.x.data(), f.y.data(), n, f.p, f.w.data(), frequencyOptions(options),
    &rngW);
  auto expanded = makeSampler<ConstantLeafSampler>(
    e.x.data(), e.y.data(), e.n, f.p, nullptr, options, &rngE);
  Draws dw(n, numSamples), de(e.n, numSamples);
  Results rw = dw.results(false, true), re = de.results(false, true);
  rw.logLikelihood = re.logLikelihood = nullptr;
  bool ran = !weighted->run(500, numSamples, rw) &&
             !expanded->run(500, numSamples, re);
  auto compare = [](const std::vector<double>& a, const std::vector<double>& b,
                    const char* name) {
    auto [ma, sa] = batchMean(a);
    auto [mb, sb] = batchMean(b);
    double z = (ma - mb) / std::sqrt(sa * sa + sb * sb);
    printf("  Student-t %s: %.4f against %.4f, z %.2f\n", name, ma, mb, z);
    return std::fabs(z) < 4.0;
  };
  std::vector<double> fitW(numSamples), fitE(numSamples);
  std::size_t row = 9, copy = 0;
  while (e.rowOf[copy] != row) ++copy;
  for (std::size_t s = 0; s < numSamples; ++s) {
    fitW[s] = dw.fits[s * n + row];
    fitE[s] = de.fits[s * e.n + copy];
  }
  // every comparison prints, so none short-circuits
  bool agreeSigma = compare(dw.sigma, de.sigma, "sigma");
  bool agreeDf = compare(dw.df, de.df, "nu");
  bool agreeFit = compare(fitW, fitE, "an outlier's fit");
  bool agree = agreeSigma && agreeDf && agreeFit;
  check(ran && agree,
        "frequency Student-t agrees with its expansion in distribution");
  ext_rng_destroy(rngW);
  ext_rng_destroy(rngE);
}

/// The log-likelihood channel under frequency: w log f(y | one observation),
/// read against the reported fits, sigma, variance surface and df.
void testLogLikelihood() {
  const std::size_t n = 30, numSamples = 5;
  Fixture f(n);
  f.w[4] = 0.0;
  f.w[5] = 1.75;
  struct Arm { const char* what; double residualDf; std::size_t varianceTrees; };
  const Arm arms[] = {{"gaussian", std::nan(""), 0},
                      {"Student-t", 0.0, 0},
                      {"heteroscedastic", std::nan(""), 3}};
  for (const Arm& arm : arms) {
    SamplerOptions options = frequencyOptions(SamplerOptions());
    options.numTrees = 10;
    options.residualDf = arm.residualDf;
    options.numVarianceTrees = arm.varianceTrees;
    bool student = std::isfinite(arm.residualDf);
    bool hetero = arm.varianceTrees > 0;
    ext_rng* rng = seededRng(1203u);
    auto sampler = makeSampler<ConstantLeafSampler>(
      f.x.data(), f.y.data(), n, f.p, f.w.data(), options, &rng);
    Draws d(n, numSamples);
    Results r = d.results(hetero, student);
    bool ran = !sampler->run(10, numSamples, r);
    double gap = 0.0;
    bool zeroIsNaN = true;
    for (std::size_t s = 0; s < numSamples; ++s)
      for (std::size_t i = 0; i < n; ++i) {
        double reported = d.logLikelihood[s * n + i];
        if (f.w[i] == 0.0) {
          zeroIsNaN = zeroIsNaN && std::isnan(reported);
          continue;
        }
        double sd = hetero ? std::sqrt(d.variance[s * n + i]) : d.sigma[s];
        double residual = f.y[i] - d.fits[s * n + i];
        double one = student
          ? logStudentT(residual / sd, d.df[s]) - std::log(sd)
          : -0.5 * std::log(2.0 * std::numbers::pi) - std::log(sd) -
              0.5 * residual * residual / (sd * sd);
        gap = std::max(gap, relativeGap(reported, f.w[i] * one));
      }
    std::string label =
      std::string("frequency log-likelihood is w times one observation's: ") +
      arm.what;
    check(ran && zeroIsNaN && gap < 1e-9, label.c_str());
    ext_rng_destroy(rng);
  }
}

/// A frequency sampler's state reinstalled into a sampler created on the same
/// data continues as the source does; a Student-t frequency state installs
/// into a precision sampler of the same data without consuming a variate.
void testStateRoundTrip() {
  const std::size_t n = 40, numSamples = 30;
  Fixture f(n);
  f.w[6] = 0.0;
  for (double residualDf : {std::nan(""), 0.0}) {
    bool student = std::isfinite(residualDf);
    SamplerOptions options;
    options.numTrees = 15;
    options.residualDf = residualDf;
    ext_rng* rngA = seededRng(5501u);
    ext_rng* rngB = seededRng(5502u);
    auto source = makeSampler<ConstantLeafSampler>(
      f.x.data(), f.y.data(), n, f.p, f.w.data(), frequencyOptions(options),
      &rngA);
    auto copy = makeSampler<ConstantLeafSampler>(
      f.x.data(), f.y.data(), n, f.p, f.w.data(), frequencyOptions(options),
      &rngB);
    Results none;
    bool ran = !source->run(40, 0, none);
    SamplerStateData state;
    source->getState(state);
    bool installed = restoresExactly(*copy, state, f.x.data());
    checkStructuralRoundTrip(state, *copy,
                             student ? "frequency Student-t state round trip"
                                     : "frequency gaussian state round trip");
    Draws da(n, numSamples), db(n, numSamples);
    Results ra = da.results(false, student), rb = db.results(false, student);
    ran = ran && !source->run(0, numSamples, ra) && !copy->run(0, numSamples, rb);
    double gap = 0.0;
    for (std::size_t s = 0; s < numSamples; ++s) {
      gap = std::max(gap, relativeGap(da.sigma[s], db.sigma[s]));
      for (std::size_t i = 0; i < n; ++i)
        gap = std::max(gap,
                       relativeGap(da.fits[s * n + i], db.fits[s * n + i]));
    }
    check(ran && installed && gap < 1e-9 && (!student || da.df == db.df),
          student ? "frequency Student-t continues from its state"
                  : "frequency gaussian continues from its state");
    if (student) {
      ext_rng* rngC = seededRng(5503u);
      auto precision = makeSampler<ConstantLeafSampler>(
        f.x.data(), f.y.data(), n, f.p, f.w.data(), options, &rngC);
      bool clean = installsClean(*precision, state, f.x.data());
      SamplerStateData after;
      precision->getState(after);
      check(clean && after.chains[0].rngState == state.chains[0].rngState &&
              after.chains[0].latents == state.chains[0].latents,
            "a frequency Student-t state installs into a precision sampler "
            "consuming no variate");
      Draws dc(n, 3);
      Results rc = dc.results(false, true);
      check(!precision->run(0, 3, rc) && allFinite(dc.sigma),
            "the precision sampler runs on from the frequency state");
      ext_rng_destroy(rngC);
    }
    ext_rng_destroy(rngA);
    ext_rng_destroy(rngB);
  }
}

}  // namespace

void runFrequencyTests() {
  std::uint64_t savedState = rngState;
  rngState = 0x5eedf00du;
  printf("frequency weights:\n");
  testGaussianAgainstExpansion();
  testHeteroscedasticAgainstExpansion();
  testLinearLeafAgainstExpansion();
  testTwoForestAgainstExpansion();
  testUnitWeightIdentity();
  testZeroWeights();
  testMaskAgainstZeroWeights();
  testResidualDfGrid();
  testStudentAgainstExpansion();
  testLogLikelihood();
  testStateRoundTrip();
  rngState = savedState;
}
