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
#include <string>

extern "C" double Rf_dbinom_raw(double, double, double, double, int);

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
/// finite. A drawn df spreads over the grid here, so a df step whose every
/// grid weight is -inf (a 0 count's term) shows as a df that never moves.
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
    if (student && arm.residualDf == 0.0)
      finite = finite && std::any_of(d.df.begin(), d.df.end(), [&](double nu) {
                 return nu != d.df.front();
               });
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

/// \p values over the source rows written out as \p e writes the rows.
std::vector<double> expanded(const Expansion& e,
                             const std::vector<double>& values) {
  std::vector<double> result(e.n);
  for (std::size_t r = 0; r < e.n; ++r) result[r] = values[e.rowOf[r]];
  return result;
}

/// The doors that re-read a scale off the rows, under frequency against the
/// same door on the expansion, draw for draw: a leaf covariate's
/// standardization re-derived after a mask change, the variance forest
/// re-anchored by a response replaced with its scale, and a numeric basis
/// replaced with its default sd re-read.
void testDoorsAgainstExpansion() {
  const std::size_t n = 40, numSamples = 150;
  Fixture f(n);
  Expansion e(f, f.w);
  std::vector<double> mask(n, 1.0);
  for (std::size_t i = 4; i < n; i += 7) mask[i] = 0.0;
  std::vector<double> maskE = expanded(e, mask);
  auto compare = [&](auto& weighted, auto& expansion, ext_rng* rngW,
                     ext_rng* rngE, bool variance, const char* what) {
    Draws dw(n, numSamples), de(e.n, numSamples);
    Results rw = dw.results(variance), re = de.results(variance);
    rw.logLikelihood = re.logLikelihood = nullptr;
    bool ran = !weighted.run(0, numSamples, rw) &&
               !expansion.run(0, numSamples, re);
    double gap = drawGap(dw, de, e, n, numSamples, variance);
    printf("  %s: largest gap %.3g over %zu sweeps\n", what, gap, numSamples);
    check(ran && gap < 1e-12, what);
    ext_rng_destroy(rngW);
    ext_rng_destroy(rngE);
  };

  {
    const std::size_t columns[] = {0};
    SamplerOptions options;
    options.numTrees = 10;
    options.leafCovariateColumns = columns;
    options.numLeafCovariates = 1;
    ext_rng* rngW = seededRng(4401u);
    ext_rng* rngE = seededRng(4401u);
    auto weighted = makeSampler<Sampler<LinearGaussianLeaf>>(
      f.x.data(), f.y.data(), n, f.p, f.w.data(), frequencyOptions(options),
      &rngW);
    auto expansion = makeSampler<Sampler<LinearGaussianLeaf>>(
      e.x.data(), e.y.data(), e.n, f.p, nullptr, options, &rngE);
    weighted->setActiveRows(mask.data());
    expansion->setActiveRows(maskE.data());
    weighted->rederiveStandardization();
    expansion->rederiveStandardization();
    compare(*weighted, *expansion, rngW, rngE, false,
            "frequency linear leaf re-derived after a mask change equals "
            "its expansion's");
  }
  {
    SamplerOptions options;
    options.numTrees = 15;
    options.numVarianceTrees = 4;
    ext_rng* rngW = seededRng(4402u);
    ext_rng* rngE = seededRng(4402u);
    auto weighted = makeSampler<ConstantLeafSampler>(
      f.x.data(), f.y.data(), n, f.p, f.w.data(), frequencyOptions(options),
      &rngW);
    auto expansion = makeSampler<ConstantLeafSampler>(
      e.x.data(), e.y.data(), e.n, f.p, nullptr, options, &rngE);
    std::vector<double> y2(f.y);
    for (std::size_t i = 0; i < n; ++i) y2[i] = 3.0 * y2[i] + f.x[n + i];
    std::vector<double> y2E = expanded(e, y2);
    weighted->setResponse(y2.data(), true);
    expansion->setResponse(y2E.data(), true);
    compare(*weighted, *expansion, rngW, rngE, true,
            "frequency heteroscedastic re-anchored by a rescaled response "
            "equals its expansion's");
  }
  {
    std::vector<double> dose(n), dose2(n);
    for (std::size_t i = 0; i < n; ++i) {
      dose[i] = 0.5 + runif01();
      dose2[i] = 2.0 * runif01() - 0.3;
    }
    std::vector<double> doseE = expanded(e, dose), dose2E = expanded(e, dose2);
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
    ext_rng* rngW = seededRng(4403u);
    ext_rng* rngE = seededRng(4403u);
    ConstantLeafSampler weighted(f.x.data(), f.y.data(), n, f.p, f.w.data(),
                                 nullptr, 1.0, sigmaDf, sigmaRawScale,
                                 frequencyOptions(options), spec(dose.data()),
                                 &rngW);
    ConstantLeafSampler expansion(e.x.data(), e.y.data(), e.n, f.p, nullptr,
                                  nullptr, 1.0, sigmaDf, sigmaRawScale,
                                  options, spec(doseE.data()), &rngE);
    bool rescaled =
      weighted.setForestBasisRescaled(1, dose2.data(), 1) ==
        BasisRescaleResult::rescaled &&
      expansion.setForestBasisRescaled(1, dose2E.data(), 1) ==
        BasisRescaleResult::rescaled;
    double scaleGap = relativeGap(weighted.forestCalibration(0, 1).basisScale,
                                  expansion.forestCalibration(0, 1).basisScale);
    check(rescaled && scaleGap < 1e-14,
          "frequency basis rescaled through the door: the expansion's scale");
    compare(weighted, expansion, rngW, rngE, false,
            "frequency two-forest after a basis rescale equals its "
            "expansion's");
  }
}

/// A logistic numeric basis's default sd under frequency is the expansion's,
/// with and without a mask, and where the mask keeps one row in it reads
/// every row at its count, as the expansion's every row.
void testLogisticBasisScale() {
  const std::size_t n = 40;
  Fixture f(n);
  std::vector<double> yBinary(n), dose(n), mask(n, 1.0), oneIn(n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    yBinary[i] = runif01() < 0.4 ? 1.0 : 0.0;
    dose[i] = 0.5 + runif01() * (i % 3 == 0 ? 3.0 : 1.0);
  }
  for (std::size_t i = 2; i < n; i += 5) mask[i] = 0.0;
  oneIn[7] = 1.0;
  Expansion e(f, f.w, &dose);
  std::vector<double> yE = expanded(e, yBinary), maskE = expanded(e, mask);
  auto scaleOf = [&](const double* x, const double* y, std::size_t rows,
                     const double* w, const double* active,
                     const double* basis, bool frequency) {
    AmplitudeSpec spec;
    spec.family = ResponseFamily::logistic;
    spec.activeRows = active;
    spec.forests.resize(2);
    spec.forests[1].kind = ForestKind::numeric;
    spec.forests[1].basis = basis;
    SamplerOptions options;
    if (frequency) options = frequencyOptions(options);
    ext_rng* rng = seededRng(7701u);
    ConstantLeafSampler sampler(x, y, rows, f.p, w, nullptr, 1.0, sigmaDf,
                                sigmaRawScale, options, spec, &rng);
    double scale = sampler.forestCalibration(0, 1).basisScale;
    ext_rng_destroy(rng);
    return scale;
  };
  double plain = scaleOf(f.x.data(), yBinary.data(), n, f.w.data(), nullptr,
                         dose.data(), true);
  double masked = scaleOf(f.x.data(), yBinary.data(), n, f.w.data(),
                          mask.data(), dose.data(), true);
  double fallback = scaleOf(f.x.data(), yBinary.data(), n, f.w.data(),
                            oneIn.data(), dose.data(), true);
  double plainE = scaleOf(e.x.data(), yE.data(), e.n, nullptr, nullptr,
                          e.extra.data(), false);
  double maskedE = scaleOf(e.x.data(), yE.data(), e.n, nullptr, maskE.data(),
                           e.extra.data(), false);
  check(relativeGap(plain, plainE) < 1e-14 &&
          relativeGap(masked, maskedE) < 1e-14 &&
          relativeGap(fallback, plainE) < 1e-14,
        "frequency logistic basis scale is the expansion's, with and without "
        "a mask and where a mask keeps one row");
}

/// Where a creation mask keeps one row in, a linear leaf's covariate center
/// and scale are read off every row at its count: the expansion's every row.
void testLeafFallbackWeighted() {
  const std::size_t n = 40;
  Fixture f(n);
  Expansion e(f, f.w);
  std::vector<double> oneIn(n, 0.0);
  oneIn[11] = 1.0;
  const std::size_t columns[] = {0};
  SamplerOptions options;
  options.numTrees = 5;
  options.leafCovariateColumns = columns;
  options.numLeafCovariates = 1;
  SamplerOptions masked = frequencyOptions(options);
  masked.activeRows = oneIn.data();
  ext_rng* rngW = seededRng(8801u);
  ext_rng* rngE = seededRng(8801u);
  auto weighted = makeSampler<Sampler<LinearGaussianLeaf>>(
    f.x.data(), f.y.data(), n, f.p, f.w.data(), masked, &rngW);
  auto expansion = makeSampler<Sampler<LinearGaussianLeaf>>(
    e.x.data(), e.y.data(), e.n, f.p, nullptr, options, &rngE);
  const auto& lw = weighted->chain(0).leaf();
  const auto& le = expansion->chain(0).leaf();
  check(relativeGap(lw.covariateMeans()[0], le.covariateMeans()[0]) < 1e-14 &&
          relativeGap(lw.covariateSds()[0], le.covariateSds()[0]) < 1e-14,
        "frequency linear leaf read off every row at its count where one row "
        "is in");
  ext_rng_destroy(rngW);
  ext_rng_destroy(rngE);
}

/// A Student-t frequency sampler handed new weights through setWeights, a
/// new pointer, is the one created with them: every count the response
/// holds follows the pointer.
void testStudentSetWeights() {
  const std::size_t n = 40, numSamples = 80;
  Fixture f(n);
  std::vector<double> w2(n);
  for (std::size_t i = 0; i < n; ++i) w2[i] = 1.0 + std::floor(5.0 * runif01());
  w2[6] = 2.5;
  SamplerOptions options = frequencyOptions(SamplerOptions());
  options.numTrees = 15;
  options.residualDf = 0.0;
  ext_rng* rngA = seededRng(5501u);
  ext_rng* rngB = seededRng(5501u);
  auto moved = makeSampler<ConstantLeafSampler>(
    f.x.data(), f.y.data(), n, f.p, f.w.data(), options, &rngA);
  auto created = makeSampler<ConstantLeafSampler>(
    f.x.data(), f.y.data(), n, f.p, w2.data(), options, &rngB);
  moved->setWeights(w2.data());
  Draws da(n, numSamples), db(n, numSamples);
  Results ra = da.results(false, true), rb = db.results(false, true);
  bool ran = !moved->run(5, numSamples, ra) && !created->run(5, numSamples, rb);
  check(ran && da.sigma == db.sigma && da.fits == db.fits && da.df == db.df &&
          da.logLikelihood == db.logLikelihood,
        "frequency Student-t given new weights by setWeights is the one "
        "created with them, draw for draw");
  ext_rng_destroy(rngA);
  ext_rng_destroy(rngB);
}

/// A row out of the likelihood under frequency, at weight 0 or masked, is
/// inert: moving its response inside the range leaves every draw and every
/// row in's fit bitwise, on gaussian and on Student-t with a drawn and a
/// fixed df.
void testRowOutIsInert() {
  const std::size_t n = 40, numSamples = 60;
  Fixture f(n);
  std::size_t lowest = 0, highest = 0;
  for (std::size_t i = 1; i < n; ++i) {
    if (f.y[i] < f.y[lowest]) lowest = i;
    if (f.y[i] > f.y[highest]) highest = i;
  }
  double middle = 0.5 * (f.y[lowest] + f.y[highest]);
  std::vector<double> mask(n, 1.0), wZero(f.w), yMoved(f.y);
  for (std::size_t i = 3; i < n; i += 6) {
    if (i == lowest || i == highest) continue;
    mask[i] = 0.0;
    wZero[i] = 0.0;
    yMoved[i] = middle + 0.1 * (f.y[i] - middle);
  }
  for (double residualDf : {std::nan(""), 0.0, 5.0}) {
    for (bool byMask : {false, true}) {
      SamplerOptions options = frequencyOptions(SamplerOptions());
      options.numTrees = 15;
      options.residualDf = residualDf;
      const double* weights = byMask ? f.w.data() : wZero.data();
      ext_rng* rngA = seededRng(811u);
      ext_rng* rngB = seededRng(811u);
      auto a = makeSampler<ConstantLeafSampler>(
        f.x.data(), f.y.data(), n, f.p, weights, options, &rngA);
      auto b = makeSampler<ConstantLeafSampler>(
        f.x.data(), yMoved.data(), n, f.p, weights, options, &rngB);
      if (byMask) {
        a->setActiveRows(mask.data());
        b->setActiveRows(mask.data());
      }
      bool student = std::isfinite(residualDf);
      Draws da(n, numSamples), db(n, numSamples);
      Results ra = da.results(false, student), rb = db.results(false, student);
      ra.logLikelihood = rb.logLikelihood = nullptr;
      bool ran = !a->run(10, numSamples, ra) && !b->run(10, numSamples, rb);
      bool inRowsSame = true;
      for (std::size_t s = 0; s < numSamples; ++s)
        for (std::size_t i = 0; i < n; ++i)
          if (wZero[i] != 0.0)
            inRowsSame = inRowsSame && da.fits[s * n + i] == db.fits[s * n + i];
      std::string label =
        std::string("frequency row out is inert: ") +
        (!student ? "gaussian" : residualDf == 0.0 ? "Student-t, drawn df"
                                                   : "Student-t, fixed df") +
        (byMask ? ", masked" : ", weight 0");
      check(ran && da.sigma == db.sigma && da.df == db.df && inRowsSame,
            label.c_str());
      ext_rng_destroy(rngA);
      ext_rng_destroy(rngB);
    }
  }
}

/// A gp leaf's default lengthscale under frequency is the median pairwise
/// distance over the rows written out: the expansion's at whole counts, on
/// both sides of the 512-row subsample; precision's bitwise at all-ones
/// counts; and a row of count 2.5 reads as two rows of 1.25 at one point.
void testGPLengthscale() {
  const std::size_t columns[] = {0};
  SamplerOptions options;
  options.numTrees = 5;
  options.leafCovariateColumns = columns;
  options.numLeafCovariates = 1;
  auto lengthscale = [&](const Fixture& f, std::size_t n, const double* x,
                         const double* y, const double* w,
                         const SamplerOptions& o) {
    ext_rng* rng = seededRng(6601u);
    auto sampler = makeSampler<Sampler<GPGaussianLeaf>>(x, y, n, f.p, w, o,
                                                        &rng);
    double value = sampler->chain(0).leaf().lengthscales()[0];
    ext_rng_destroy(rng);
    return value;
  };
  for (std::size_t n : {std::size_t(40), std::size_t(300)}) {
    Fixture f(n);
    Expansion e(f, f.w);
    double weighted = lengthscale(f, n, f.x.data(), f.y.data(), f.w.data(),
                                  frequencyOptions(options));
    double expansion =
      lengthscale(f, e.n, e.x.data(), e.y.data(), nullptr, options);
    std::vector<double> ones(n, 1.0);
    double unitFrequency = lengthscale(f, n, f.x.data(), f.y.data(),
                                       ones.data(), frequencyOptions(options));
    double unitPrecision =
      lengthscale(f, n, f.x.data(), f.y.data(), ones.data(), options);
    double precision =
      lengthscale(f, n, f.x.data(), f.y.data(), f.w.data(), options);
    printf("  gp lengthscale, %zu rows (%zu written out): %.6f against %.6f, "
           "precision %.6f\n",
           n, e.n, weighted, expansion, precision);
    std::string label = "frequency gp lengthscale is the expansion's, " +
                        std::to_string(n) + " rows";
    check(relativeGap(weighted, expansion) < 1e-12 &&
            unitFrequency == unitPrecision,
          label.c_str());
  }
  {
    const std::size_t n = 30;
    Fixture f(n);
    std::vector<double> w(n, 1.0);
    w[5] = 2.5;
    // row 5 split into two adjacent rows of 1.25
    Fixture s(n + 1);
    s.x.assign((n + 1) * f.p, 0.0);
    s.y.assign(n + 1, 0.0);
    s.w.assign(n + 1, 1.0);
    for (std::size_t r = 0, i = 0; r <= n; ++r) {
      for (std::size_t j = 0; j < f.p; ++j) s.x[j * (n + 1) + r] = f.x[j * n + i];
      s.y[r] = f.y[i];
      if (i == 5) s.w[r] = 1.25;
      if (!(i == 5 && r == 5)) ++i;
    }
    double whole = lengthscale(f, n, f.x.data(), f.y.data(), w.data(),
                               frequencyOptions(options));
    double split = lengthscale(s, n + 1, s.x.data(), s.y.data(), s.w.data(),
                               frequencyOptions(options));
    check(relativeGap(whole, split) < 1e-12,
          "frequency gp lengthscale: a count of 2.5 reads as two of 1.25");
  }
}

/// The counted median pairwise distance computed independently of the
/// leaf: the rows read (mask kept, count positive) laid out in row order on
/// a line of length sum(counts), the stride over its ceiling, each row
/// weighted by its overlap with the unit cells starting at every stride-th
/// point, thinned where more than 512 carry a weight, pairs v_a v_b and self
/// pairs v (v - 1) / 2 where positive, the median at the unweighted median's
/// ranks.
double referenceCountedMedian(const std::vector<double>& u,
                              const std::vector<double>& counts,
                              const std::vector<double>& mask) {
  std::vector<double> uu, ww;
  for (std::size_t i = 0; i < u.size(); ++i)
    if (mask[i] != 0.0 && counts[i] > 0.0) {
      uu.push_back(u[i]);
      ww.push_back(counts[i]);
    }
  double total = 0.0;
  for (double w : ww) total += w;
  double numWritten = std::ceil(total);
  double s =
    numWritten < 1.0 ? 1.0 : 1.0 + std::floor((numWritten - 1.0) / 512.0);
  std::vector<double> vu, vv;
  double start = 0.0;
  for (std::size_t i = 0; i < ww.size(); ++i) {
    double end = start + ww[i], v = 0.0;
    if (s == 1.0) {
      v = ww[i];
    } else {
      for (double k = 0.0; k * s < end; k += 1.0) {
        double lo = std::max(start, k * s), hi = std::min(end, k * s + 1.0);
        if (hi > lo) v += hi - lo;
      }
    }
    if (v > 0.0) {
      vu.push_back(uu[i]);
      vv.push_back(v);
    }
    start = end;
  }
  if (vu.size() > 512) {
    std::size_t thin = 1 + (vu.size() - 1) / 512;
    std::vector<double> a, b;
    for (std::size_t k = 0; k < vu.size(); k += thin) {
      a.push_back(vu[k]);
      b.push_back(vv[k]);
    }
    vu.swap(a);
    vv.swap(b);
  }
  double sum = 0.0;
  for (double v : vv) sum += v;
  if (sum < 2.0) return 1.0;
  std::vector<std::pair<double, double>> d;
  for (std::size_t a = 0; a < vu.size(); ++a) {
    if (vv[a] * (vv[a] - 1.0) > 0.0)
      d.emplace_back(0.0, 0.5 * vv[a] * (vv[a] - 1.0));
    for (std::size_t b = a + 1; b < vu.size(); ++b)
      d.emplace_back(std::fabs(vu[a] - vu[b]), vv[a] * vv[b]);
  }
  std::sort(d.begin(), d.end());
  double numPairs = 0.0;
  for (const auto& e : d) numPairs += e.second;
  auto at = [&](double r) {
    double c = 0.0;
    for (const auto& e : d) {
      c += e.second;
      if (c > r) return e.first;
    }
    return d.back().first;
  };
  double lo = std::floor(0.5 * (numPairs - 1.0)), hi = std::floor(0.5 * numPairs);
  double m = lo == hi ? at(lo) : 0.5 * (at(lo) + at(hi));
  return m > 0.0 ? m : 1.0;
}

/// A gp leaf's default lengthscale at a creation mask \p mask (null for
/// none), and in \p u (when non-null) the standardized covariate it read.
double gpLengthscale(const double* x, const double* y, std::size_t n,
                     std::size_t p, const double* w, const double* mask,
                     bool frequency, std::vector<double>* u) {
  const std::size_t columns[] = {0};
  SamplerOptions options;
  options.numTrees = 5;
  options.leafCovariateColumns = columns;
  options.numLeafCovariates = 1;
  options.activeRows = mask;
  if (frequency) options = frequencyOptions(options);
  ext_rng* rng = seededRng(6601u);
  auto sampler =
    makeSampler<Sampler<GPGaussianLeaf>>(x, y, n, p, w, options, &rng);
  const auto& leaf = sampler->chain(0).leaf();
  if (u != nullptr) {
    u->resize(n);
    for (std::size_t i = 0; i < n; ++i)
      (*u)[i] = (x[i] - leaf.covariateMeans()[0]) / leaf.covariateSds()[0];
  }
  double value = leaf.lengthscales()[0];
  ext_rng_destroy(rng);
  return value;
}

/// The gp lengthscale under frequency with zero weights and a creation mask
/// against the expansion's masked copies (whole counts, past the subsample),
/// and at real counts on both sides of the subsample against the counted
/// median computed independently over the fit's own standardized covariate.
void testGPLengthscaleMaskedAndReal() {
  {
    const std::size_t n = 200;
    Fixture f(n);
    std::vector<double> mask(n, 1.0);
    for (std::size_t i = 0; i < n; ++i) {
      f.w[i] = 1.0 + std::floor(5.0 * runif01());
      if (i % 9 == 4) f.w[i] = 0.0;
      if (i % 4 == 1) mask[i] = 0.0;
    }
    Expansion e(f, f.w);
    std::vector<double> maskE = expanded(e, mask);
    double weighted = gpLengthscale(f.x.data(), f.y.data(), n, f.p, f.w.data(),
                                    mask.data(), true, nullptr);
    double expansion = gpLengthscale(e.x.data(), e.y.data(), e.n, f.p, nullptr,
                                     maskE.data(), false, nullptr);
    printf("  gp lengthscale, zero weights and a mask (%zu written out): %.6f "
           "against %.6f\n",
           e.n, weighted, expansion);
    check(relativeGap(weighted, expansion) < 1e-12,
          "frequency gp lengthscale with zero weights and a mask is the "
          "expansion's masked copies'");
  }
  for (std::size_t n : {std::size_t(80), std::size_t(300)}) {
    Fixture f(n);
    std::vector<double> mask(n, 1.0);
    for (std::size_t i = 0; i < n; ++i) {
      f.w[i] = 0.3 + 6.7 * runif01();
      if (i % 5 == 2) mask[i] = 0.0;
    }
    std::vector<double> u;
    double value = gpLengthscale(f.x.data(), f.y.data(), n, f.p, f.w.data(),
                                 mask.data(), true, &u);
    double reference = referenceCountedMedian(u, f.w, mask);
    std::string label = std::string("frequency gp lengthscale at real counts "
                                    "is the counted median, ") +
                        (n == 80 ? "stride 1" : "stride above 1");
    check(relativeGap(value, reference) < 1e-12, label.c_str());
  }
}

/// A gaussian two-forest fit whose rows in hold one response value while
/// the rows of weight 0 vary: the response's range falls back to every row,
/// and under frequency the basis scale and the anchor read every row at its
/// count, the expansion's every row.
void testGaussianFallbackWeighted() {
  const std::size_t n = 40;
  Fixture f(n);
  std::vector<double> dose(n);
  for (std::size_t i = 0; i < n; ++i) {
    f.w[i] = i % 4 == 1 ? 0.0 : 1.0 + std::floor(4.0 * runif01());
    f.y[i] = f.w[i] == 0.0 ? 3.0 * runif01() : 0.7;
    dose[i] = 0.5 + runif01() * (i % 3 == 0 ? 3.0 : 1.0);
  }
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
  ext_rng* rngW = seededRng(9901u);
  ext_rng* rngE = seededRng(9901u);
  ConstantLeafSampler weighted(f.x.data(), f.y.data(), n, f.p, f.w.data(),
                               nullptr, 1.0, sigmaDf, sigmaRawScale,
                               frequencyOptions(options), spec(dose.data()),
                               &rngW);
  ConstantLeafSampler expansion(e.x.data(), e.y.data(), e.n, f.p, nullptr,
                                nullptr, 1.0, sigmaDf, sigmaRawScale, options,
                                spec(e.extra.data()), &rngE);
  ForestCalibration cw = weighted.forestCalibration(0, 1),
                    ce = expansion.forestCalibration(0, 1);
  check(relativeGap(cw.basisScale, ce.basisScale) < 1e-14,
        "frequency gaussian basis scale where the rows in hold one value is "
        "the expansion's");
  check(relativeGap(cw.mapAnchor, ce.mapAnchor) < 1e-14,
        "frequency gaussian anchor where the rows in hold one value is the "
        "expansion's");
  ext_rng_destroy(rngW);
  ext_rng_destroy(rngE);
}

/// The gaussian anchor where the counts over the rows read total no more
/// than 1, the rows varying: the sd has no positive divisor, so the anchor
/// is the one a response with no spread gets. At 32 rows of weight 1/32,
/// a total of exactly 1, and at whole weights rescaled to sum to 1.
void testAnchorTotalAtMostOne() {
  auto spec = []() {
    AmplitudeSpec s;
    s.forests.resize(2);
    s.forests[0].forest.numTrees = 15;
    s.forests[1].forest.numTrees = 10;
    return s;
  };
  auto anchorAt = [&](const Fixture& f, const std::vector<double>& w) {
    ext_rng* rng = seededRng(4242u);
    ConstantLeafSampler sampler(f.x.data(), f.y.data(), f.n, f.p, w.data(),
                                nullptr, 1.0, sigmaDf, sigmaRawScale,
                                frequencyOptions(SamplerOptions()), spec(),
                                &rng);
    double anchor = sampler.forestCalibration(0, 1).mapAnchor;
    ext_rng_destroy(rng);
    return anchor;
  };
  Fixture exact(32);
  std::vector<double> w(32, 1.0 / 32.0);
  double anchor = anchorAt(exact, w);
  check(anchor == singleForestDefaultSize(),
        "frequency gaussian anchor at a total count of exactly 1 is the "
        "no-spread one");
  Fixture rescaled(40);
  double sum = 0.0;
  for (double v : rescaled.w) sum += v;
  std::vector<double> u(40);
  double total = 0.0;
  for (std::size_t i = 0; i < 40; ++i) {
    u[i] = rescaled.w[i] / sum;
    total += u[i];
  }
  anchor = anchorAt(rescaled, u);
  printf("  anchor at weights rescaled to sum to 1: total - 1 = %.3g, "
         "anchor %.6g\n", total - 1.0, anchor);
  check(anchor == singleForestDefaultSize(),
        "frequency gaussian anchor at whole weights rescaled to sum to 1 is "
        "the no-spread one");
  // the two siblings read the same total: a leaf covariate's sd is its
  // placeholder and a basis column has no scale
  double mean = 0.0, sd = 0.0;
  standardizationMomentsForColumn(rescaled.y.data(), 40, &mean, &sd, nullptr,
                                  u.data());
  bool refused = false;
  try {
    countedBasisColumnScale(rescaled.y.data(), 40, nullptr, u.data());
  } catch (const BasisScaleError& error) {
    // the column varies: the refusal names the total, not a constant
    refused = error.totalAtMostOne && std::isnan(error.value) &&
              std::string(error.what()).find("total count of at most 1") !=
                std::string::npos;
  }
  check(sd == 1.0 && refused,
        "frequency leaf moments and basis scale at weights rescaled to sum "
        "to 1 take no scale, the basis refused by its total");
}

// --- the Polya-Gamma families: nbinom and multinomial weights are counts. A
// row of count w draws the sum of its copies' unit draws, consumed in row
// order, so a fit equals the expansion's draw for draw, to rounding. ---

/// nbinom counts, row-major in nothing: y whole, x as Fixture's.
struct CountFixture {
  std::size_t n = 0, p = 2;
  std::vector<double> x, y, w;

  explicit CountFixture(std::size_t n_) : n(n_) {
    Fixture f(n_);
    x = f.x;
    w = f.w;
    y.resize(n);
    for (std::size_t i = 0; i < n; ++i)
      y[i] = std::floor(6.0 * std::exp(std::sin(6.0 * x[i])) * runif01());
  }
};

std::unique_ptr<ConstantLeafSampler> makeCountSampler(
    const double* x, const double* y, std::size_t n, std::size_t p,
    const double* w, double shape, ext_rng** rng) {
  SamplerOptions options;
  options.numTrees = 20;
  options.shape = shape;
  return std::make_unique<ConstantLeafSampler>(x, y, n, p, w, nullptr,
                                               ResponseFamily::nbinom, 1.0,
                                               sigmaDf, 1.0, options, rng);
}

/// A count run's reported draws.
struct CountDraws {
  std::vector<double> fits, logLikelihood, shape;

  CountDraws(std::size_t n, std::size_t numSamples)
    : fits(n * numSamples), logLikelihood(n * numSamples), shape(numSamples) {}

  Results results() {
    Results r;
    r.trainingFits = fits.data();
    r.logLikelihood = logLikelihood.data();
    r.shape = shape.data();
    return r;
  }
};

/// The expansion's rows of \p c, written out \p weights[i] times.
Expansion countExpansion(const CountFixture& c,
                         const std::vector<double>& weights) {
  Fixture f(1);
  f.n = c.n;
  f.p = c.p;
  f.x = c.x;
  f.y = c.y;
  return Expansion(f, weights);
}

/// An nbinom row of count w against w adjacent copies, at a fixed and at a
/// drawn r: the fits, r, and a row's log-likelihood against its copies' sum.
void testNBinomAgainstExpansion() {
  const std::size_t n = 40, numSamples = 200;
  CountFixture c(n);
  Expansion e = countExpansion(c, c.w);
  for (double shape : {3.0, -1.0}) {
    ext_rng* rngW = seededRng(7301u);
    ext_rng* rngE = seededRng(7301u);
    auto weighted = makeCountSampler(c.x.data(), c.y.data(), n, c.p,
                                     c.w.data(), shape, &rngW);
    auto expansion = makeCountSampler(e.x.data(), e.y.data(), e.n, c.p,
                                      nullptr, shape, &rngE);
    CountDraws dw(n, numSamples), de(e.n, numSamples);
    Results rw = dw.results(), re = de.results();
    bool ran = !weighted->run(0, numSamples, rw) &&
               !expansion->run(0, numSamples, re);
    double gap = 0.0;
    bool shapeSame = true;
    for (std::size_t s = 0; s < numSamples; ++s) {
      shapeSame = shapeSame && dw.shape[s] == de.shape[s];
      std::vector<double> summed(n, 0.0);
      for (std::size_t r = 0; r < e.n; ++r) {
        std::size_t i = e.rowOf[r];
        gap = std::max(gap, relativeGap(dw.fits[s * n + i],
                                        de.fits[s * e.n + r]));
        summed[i] += de.logLikelihood[s * e.n + r];
      }
      for (std::size_t i = 0; i < n; ++i)
        gap = std::max(gap,
                       relativeGap(dw.logLikelihood[s * n + i], summed[i]));
    }
    const char* what = shape > 0.0 ? "nbinom at a fixed r" : "nbinom at a drawn r";
    printf("  %s: largest gap %.3g over %zu sweeps\n", what, gap, numSamples);
    check(ran && shapeSame && gap < 1e-10,
          shape > 0.0 ? "frequency nbinom at a fixed r equals its expansion "
                        "draw for draw"
                      : "frequency nbinom at a drawn r equals its expansion "
                        "draw for draw");
    ext_rng_destroy(rngW);
    ext_rng_destroy(rngE);
  }
}

/// Whether two runs' fits over \p n rows and \p locations per row agree:
/// bitwise at the rows \p out leaves in, and to rounding at the rows out,
/// whose fit the trees carry through a working response nothing weights.
bool agreeOnRowsIn(const std::vector<double>& a, const std::vector<double>& b,
                   const std::vector<double>& out, std::size_t n) {
  if (a.size() != b.size()) return false;
  for (std::size_t j = 0; j < a.size(); ++j) {
    bool rowOut = out[j % n] == 0.0;
    if (rowOut ? !(relativeGap(a[j], b[j]) < 1e-12) : a[j] != b[j])
      return false;
  }
  return true;
}

/// The largest gap between two runs' draws.
double largestGap(const std::vector<double>& a, const std::vector<double>& b) {
  double gap = a.size() == b.size() ? 0.0 : HUGE_VAL;
  for (std::size_t j = 0; j < a.size() && j < b.size(); ++j)
    gap = std::max(gap, relativeGap(a[j], b[j]));
  return gap;
}

/// An nbinom row of count 0 is out of the likelihood: the fit equals one
/// whose row is masked, its count 1, with the center re-derived over the
/// rows kept, bitwise at the rows in; every fit, latent and kept row's
/// log-likelihood is finite and the row's log-likelihood NaN.
void testNBinomZeroWeights() {
  const std::size_t n = 40, numSamples = 60;
  CountFixture c(n);
  std::vector<double> wZero(c.w), wOne(c.w), mask(n, 1.0);
  for (std::size_t i = 2; i < n; i += 7) {
    wZero[i] = 0.0;
    wOne[i] = 1.0;
    mask[i] = 0.0;
  }
  for (double shape : {3.0, -1.0}) {
    ext_rng* rngA = seededRng(7401u);
    ext_rng* rngB = seededRng(7401u);
    auto zero = makeCountSampler(c.x.data(), c.y.data(), n, c.p, wZero.data(),
                                 shape, &rngA);
    auto masked = makeCountSampler(c.x.data(), c.y.data(), n, c.p, wOne.data(),
                                   shape, &rngB);
    masked->setActiveRows(mask.data());
    masked->rederiveScale();
    CountDraws da(n, numSamples), db(n, numSamples);
    Results ra = da.results(), rb = db.results();
    bool ran = !zero->run(5, numSamples, ra) && !masked->run(5, numSamples, rb);
    bool finite = allFinite(da.fits);
    for (std::size_t i = 0; i < n; ++i)
      finite = finite && std::isfinite(zero->latents(0)[i]) &&
               zero->latents(0)[i] > 0.0;
    bool logLikelihood = true;
    for (std::size_t s = 0; s < numSamples; ++s)
      for (std::size_t i = 0; i < n; ++i)
        logLikelihood = logLikelihood &&
          (wZero[i] == 0.0 ? std::isnan(da.logLikelihood[s * n + i])
                           : std::isfinite(da.logLikelihood[s * n + i]));
    check(ran && finite && logLikelihood,
          shape > 0.0 ? "frequency nbinom at weight 0, a fixed r: finite "
                        "throughout, the row's log-likelihood NaN"
                      : "frequency nbinom at weight 0, a drawn r: finite "
                        "throughout, the row's log-likelihood NaN");
    check(agreeOnRowsIn(da.fits, db.fits, mask, n) && da.shape == db.shape,
          shape > 0.0 ? "frequency nbinom at weight 0 equals the masked row, "
                        "a fixed r"
                      : "frequency nbinom at weight 0 equals the masked row, "
                        "a drawn r");
    ext_rng_destroy(rngA);
    ext_rng_destroy(rngB);
  }
}

/// nbinom setWeights draws nothing; with the center re-derived it equals the
/// sampler created with those counts, draw for draw; and a state round trip
/// continues to rounding.
void testNBinomSetWeightsAndState() {
  const std::size_t n = 40, numSamples = 40;
  CountFixture c(n);
  std::vector<double> w2(n);
  for (std::size_t i = 0; i < n; ++i) w2[i] = i % 5 == 0 ? 0.0 : 5.0 - c.w[i];
  ext_rng* rngA = seededRng(7501u);
  ext_rng* rngB = seededRng(7501u);
  ext_rng* rngC = seededRng(7501u);
  auto changed = makeCountSampler(c.x.data(), c.y.data(), n, c.p, c.w.data(),
                                  -1.0, &rngA);
  auto twin = makeCountSampler(c.x.data(), c.y.data(), n, c.p, c.w.data(),
                               -1.0, &rngB);
  auto created = makeCountSampler(c.x.data(), c.y.data(), n, c.p, w2.data(),
                                  -1.0, &rngC);
  changed->setWeights(w2.data());
  check(ext_rng_simulateContinuousUniform(rngA) ==
          ext_rng_simulateContinuousUniform(rngB),
        "frequency nbinom setWeights consumes no variate");
  ext_rng_setSeed(rngA, 7501u);
  changed->rederiveScale();
  CountDraws da(n, numSamples), dc(n, numSamples);
  Results ra = da.results(), rc = dc.results();
  bool ran = !changed->run(0, numSamples, ra) &&
             !created->run(0, numSamples, rc);
  check(ran && da.fits == dc.fits && da.shape == dc.shape,
        "frequency nbinom after setWeights and a re-derived center equals the "
        "sampler created with those counts");

  ext_rng* rngD = seededRng(7502u);
  auto copy = makeCountSampler(c.x.data(), c.y.data(), n, c.p, w2.data(),
                               -1.0, &rngD);
  SamplerStateData state;
  changed->getState(state);
  bool installed = restoresExactly(*copy, state, c.x.data());
  checkStructuralRoundTrip(state, *copy, "frequency nbinom state round trip");
  CountDraws d1(n, numSamples), d2(n, numSamples);
  Results r1 = d1.results(), r2 = d2.results();
  ran = !changed->run(0, numSamples, r1) && !copy->run(0, numSamples, r2);
  check(ran && installed && largestGap(d1.fits, d2.fits) < 1e-9 &&
          d1.shape == d2.shape,
        "frequency nbinom continues from its state");
  ext_rng_destroy(rngA);
  ext_rng_destroy(rngB);
  ext_rng_destroy(rngC);
  ext_rng_destroy(rngD);
}

/// Multinomial rows: n x K category-major counts of 0 to 2 per category, and
/// whole weights 1 to 4.
struct MultinomialFixture {
  static constexpr std::size_t K = 3;
  std::size_t n = 0, p = 2;
  std::vector<double> x, w;
  std::vector<int> counts, trials;

  explicit MultinomialFixture(std::size_t n_) : n(n_) {
    Fixture f(n_);
    x = f.x;
    w = f.w;
    counts.assign(n * K, 0);
    trials.assign(n, 0);
    for (std::size_t i = 0; i < n; ++i) {
      double xi = x[i];
      std::size_t favored = xi < 0.34 ? 0 : (xi < 0.67 ? 1 : 2);
      for (std::size_t k = 0; k < K; ++k) {
        int count = static_cast<int>(std::floor(3.0 * runif01()));
        if (k == favored) ++count;
        counts[k * n + i] = count;
        trials[i] += count;
      }
    }
  }

  /// The rows written out weights[i] times, at unit weight.
  MultinomialFixture expanded(const std::vector<double>& weights,
                              std::vector<std::size_t>& rowOf) const {
    MultinomialFixture e(0);
    rowOf.clear();
    for (std::size_t i = 0; i < n; ++i)
      for (int c = 0; c < static_cast<int>(weights[i]); ++c) rowOf.push_back(i);
    e.n = rowOf.size();
    e.x.resize(e.n * p);
    e.counts.resize(e.n * K);
    e.trials.resize(e.n);
    for (std::size_t r = 0; r < e.n; ++r) {
      for (std::size_t j = 0; j < p; ++j) e.x[j * e.n + r] = x[j * n + rowOf[r]];
      for (std::size_t k = 0; k < K; ++k)
        e.counts[k * e.n + r] = counts[k * n + rowOf[r]];
      e.trials[r] = trials[rowOf[r]];
    }
    return e;
  }
};

std::unique_ptr<Sampler<ConstantGaussianLeaf>> makeMultinomialSampler(
    const MultinomialFixture& m, const int* counts, const int* trials,
    const double* weights, ext_rng** rng, const double* offset = nullptr) {
  SamplerOptions options;
  options.numTrees = 15;
  MultinomialSpec spec;
  spec.numCategories = MultinomialFixture::K;
  spec.counts = counts;
  spec.trials = trials;
  spec.weights = weights;
  spec.offset = offset;
  spec.forest.numTrees = 15;
  return std::make_unique<Sampler<ConstantGaussianLeaf>>(m.x.data(), m.n, m.p,
                                                         options, spec, rng);
}

/// \p numSamples recorded draws of the K softmax probabilities per row.
std::vector<double> multinomialDraws(Sampler<ConstantGaussianLeaf>& sampler,
                                     std::size_t n, std::size_t numSamples,
                                     bool* ran) {
  std::vector<double> fits(n * MultinomialFixture::K * numSamples);
  Results r;
  r.trainingFits = fits.data();
  r.numReportedLocations = MultinomialFixture::K;
  *ran = !sampler.run(0, numSamples, r) && *ran;
  return fits;
}

/// A multinomial row of count w against w adjacent copies, draw for draw,
/// and against the counts and trials times w, bitwise; with a category offset
/// the working response takes its own branch, so both are run.
void testMultinomialAgainstExpansion(bool withOffset) {
  const std::size_t n = 40, numSamples = 150, K = MultinomialFixture::K;
  MultinomialFixture m(n);
  std::vector<std::size_t> rowOf;
  MultinomialFixture e = m.expanded(m.w, rowOf);
  std::vector<int> countsW(m.counts), trialsW(m.trials);
  for (std::size_t i = 0; i < n; ++i) {
    int w = static_cast<int>(m.w[i]);
    trialsW[i] *= w;
    for (std::size_t k = 0; k < K; ++k) countsW[k * n + i] *= w;
  }
  std::vector<double> offset(n * K), offsetE(e.n * K);
  for (std::size_t k = 0; k < K; ++k) {
    for (std::size_t i = 0; i < n; ++i)
      offset[k * n + i] = 0.4 * (static_cast<double>(k) - 1.0) + 0.3 * m.x[i];
    for (std::size_t r = 0; r < e.n; ++r)
      offsetE[k * e.n + r] = offset[k * n + rowOf[r]];
  }
  const double* o = withOffset ? offset.data() : nullptr;
  const double* oE = withOffset ? offsetE.data() : nullptr;
  ext_rng* rngW = seededRng(7601u);
  ext_rng* rngE = seededRng(7601u);
  ext_rng* rngT = seededRng(7601u);
  auto weighted = makeMultinomialSampler(m, m.counts.data(), m.trials.data(),
                                         m.w.data(), &rngW, o);
  auto expansion = makeMultinomialSampler(e, e.counts.data(), e.trials.data(),
                                          nullptr, &rngE, oE);
  auto times = makeMultinomialSampler(m, countsW.data(), trialsW.data(),
                                      nullptr, &rngT, o);
  bool ran = true;
  std::vector<double> fw = multinomialDraws(*weighted, n, numSamples, &ran);
  std::vector<double> fe = multinomialDraws(*expansion, e.n, numSamples, &ran);
  std::vector<double> ft = multinomialDraws(*times, n, numSamples, &ran);
  double gap = 0.0;
  for (std::size_t s = 0; s < numSamples; ++s)
    for (std::size_t k = 0; k < K; ++k)
      for (std::size_t r = 0; r < e.n; ++r)
        gap = std::max(gap,
                       relativeGap(fw[(s * K + k) * n + rowOf[r]],
                                   fe[(s * K + k) * e.n + r]));
  const char* label = withOffset ? " under a category offset" : "";
  printf("  multinomial%s: largest gap %.3g over %zu sweeps\n", label, gap,
         numSamples);
  check(ran && gap < 1e-10,
        withOffset ? "frequency multinomial under a category offset equals its "
                     "expansion draw for draw"
                   : "frequency multinomial equals its expansion draw for draw");
  check(fw == ft,
        withOffset ? "frequency multinomial under a category offset equals the "
                     "counts and trials times w, bitwise"
                   : "frequency multinomial equals the counts and trials times "
                     "w, bitwise");
  ext_rng_destroy(rngW);
  ext_rng_destroy(rngE);
  ext_rng_destroy(rngT);
}

/// A multinomial row of count 0 is out as a masked row is, bitwise at the
/// rows in, and finite; setWeights draws nothing and equals the sampler
/// created with the counts; a state round trip continues to rounding.
void testMultinomialZeroWeightsAndState() {
  const std::size_t n = 40, numSamples = 40, K = MultinomialFixture::K;
  MultinomialFixture m(n);
  std::vector<double> wZero(m.w), wOne(m.w), mask(n, 1.0);
  for (std::size_t i = 1; i < n; i += 6) {
    wZero[i] = 0.0;
    wOne[i] = 1.0;
    mask[i] = 0.0;
  }
  ext_rng* rngA = seededRng(7701u);
  ext_rng* rngB = seededRng(7701u);
  auto zero = makeMultinomialSampler(m, m.counts.data(), m.trials.data(),
                                     wZero.data(), &rngA);
  auto masked = makeMultinomialSampler(m, m.counts.data(), m.trials.data(),
                                       wOne.data(), &rngB);
  masked->setActiveRows(mask.data());
  bool ran = true;
  std::vector<double> fz = multinomialDraws(*zero, n, numSamples, &ran);
  std::vector<double> fm = multinomialDraws(*masked, n, numSamples, &ran);
  check(ran && allFinite(fz) && agreeOnRowsIn(fz, fm, mask, n),
        "frequency multinomial at weight 0 is finite and equals the masked "
        "row, bitwise at the rows in");

  ext_rng* rngC = seededRng(7702u);
  ext_rng* rngD = seededRng(7702u);
  ext_rng* rngE = seededRng(7702u);
  auto changed = makeMultinomialSampler(m, m.counts.data(), m.trials.data(),
                                        m.w.data(), &rngC);
  auto twin = makeMultinomialSampler(m, m.counts.data(), m.trials.data(),
                                     m.w.data(), &rngD);
  auto created = makeMultinomialSampler(m, m.counts.data(), m.trials.data(),
                                        wZero.data(), &rngE);
  changed->setWeights(wZero.data());
  check(ext_rng_simulateContinuousUniform(rngC) ==
          ext_rng_simulateContinuousUniform(rngD),
        "frequency multinomial setWeights consumes no variate");
  ext_rng_setSeed(rngC, 7702u);
  ran = true;
  std::vector<double> fc = multinomialDraws(*changed, n, numSamples, &ran);
  std::vector<double> fr = multinomialDraws(*created, n, numSamples, &ran);
  check(ran && fc == fr,
        "frequency multinomial after setWeights equals the sampler created "
        "with those counts");

  ext_rng* rngF = seededRng(7703u);
  auto copy = makeMultinomialSampler(m, m.counts.data(), m.trials.data(),
                                     wZero.data(), &rngF);
  SamplerStateData state;
  changed->getState(state);
  bool installed = restoresExactly(*copy, state, m.x.data());
  ran = true;
  std::vector<double> f1 = multinomialDraws(*changed, n, numSamples, &ran);
  std::vector<double> f2 = multinomialDraws(*copy, n, numSamples, &ran);
  check(ran && installed && largestGap(f1, f2) < 1e-9 &&
          f1.size() == n * K * numSamples,
        "frequency multinomial continues from its state");
  ext_rng_destroy(rngA);
  ext_rng_destroy(rngB);
  ext_rng_destroy(rngC);
  ext_rng_destroy(rngD);
  ext_rng_destroy(rngE);
  ext_rng_destroy(rngF);
}

// --- the binomial row: a logistic row of k successes in n trials stands for
// n binary rows of the same x, k of them successes, and draws its n unit
// Polya-Gamma variates at the one psi those rows share. ---

/// Binomial rows: trials 1 to 6 with 0 < k < n wherever n > 1 (\p interior),
/// or any k in 0..n, and whole weights 1 to 4; x as Fixture's.
struct BinomialFixture {
  std::size_t n = 0, p = 2;
  std::vector<double> x, w, k, trials;

  BinomialFixture(std::size_t n_, bool interior) : n(n_) {
    Fixture f(n_);
    x = f.x;
    w = f.w;
    k.resize(n);
    trials.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
      double size = interior ? 2.0 + std::floor(5.0 * runif01())
                             : 1.0 + std::floor(6.0 * runif01());
      double prob = 1.0 / (1.0 + std::exp(-3.0 * std::sin(6.0 * x[i])));
      double successes = 0.0;
      for (int t = 0; t < static_cast<int>(size); ++t)
        successes += runif01() < prob ? 1.0 : 0.0;
      if (interior)
        successes = std::min(std::max(successes, 1.0), size - 1.0);
      trials[i] = size;
      k[i] = successes;
    }
  }
};

std::unique_ptr<ConstantLeafSampler> makeBinomialSampler(
    const double* x, const double* y, std::size_t n, std::size_t p,
    const double* w, const double* trials, ext_rng** rng,
    const double* offset = nullptr) {
  SamplerOptions options;
  options.numTrees = 20;
  options.trials = trials;
  return std::make_unique<ConstantLeafSampler>(x, y, n, p, w, offset,
                                               ResponseFamily::logistic, 1.0,
                                               sigmaDf, 1.0, options, rng);
}

/// A binomial row of k successes in n trials against two adjacent binary
/// rows of weights k (a success) and n - k (a failure), 0 < k < n: the same
/// unit draws at the same psi, so the fits agree to the association of the
/// omega sums over 300 sweeps and both runs consume the same variates.
void testBinomialAgainstSplitRows() {
  const std::size_t n = 40, numSamples = 300;
  BinomialFixture b(n, true);
  std::size_t m = 2 * n;
  std::vector<double> xs(m * b.p), ys(m), ws(m);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < b.p; ++j)
      xs[j * m + 2 * i] = xs[j * m + 2 * i + 1] = b.x[j * n + i];
    ys[2 * i] = 1.0;
    ys[2 * i + 1] = 0.0;
    ws[2 * i] = b.k[i];
    ws[2 * i + 1] = b.trials[i] - b.k[i];
  }
  ext_rng* rngB = seededRng(7601u);
  ext_rng* rngS = seededRng(7601u);
  auto binomial = makeBinomialSampler(b.x.data(), b.k.data(), n, b.p, nullptr,
                                      b.trials.data(), &rngB);
  auto split = makeBinomialSampler(xs.data(), ys.data(), m, b.p, ws.data(),
                                   nullptr, &rngS);
  CountDraws db(n, numSamples), ds(m, numSamples);
  Results rb = db.results(), rs = ds.results();
  rb.shape = rs.shape = nullptr;
  bool ran = !binomial->run(0, numSamples, rb) && !split->run(0, numSamples, rs);
  double gap = 0.0;
  for (std::size_t s = 0; s < numSamples; ++s)
    for (std::size_t r = 0; r < m; ++r)
      gap = std::max(gap, relativeGap(db.fits[s * n + r / 2],
                                      ds.fits[s * m + r]));
  bool sameVariates = ext_rng_simulateContinuousUniform(rngB) ==
                      ext_rng_simulateContinuousUniform(rngS);
  printf("  binomial row against its split rows: largest gap %.3g over %zu "
         "sweeps\n", gap, numSamples);
  check(ran && gap < 1e-10 && sameVariates,
        "a binomial row equals two rows of weights k and n - k to 1e-10, "
        "consuming the same variates");
  ext_rng_destroy(rngB);
  ext_rng_destroy(rngS);
}

/// A binomial row at frequency weight w against the row of k w successes in
/// n w trials, and a binomial response of one trial a row against the binary
/// response, weighted and not: fits bitwise, and the one-trial
/// log-likelihood bitwise the binary one's.
void testBinomialScaledAndBinary() {
  const std::size_t n = 40, numSamples = 100;
  BinomialFixture b(n, false);
  std::vector<double> kw(n), nw(n), ones(n, 1.0), yBinary(n);
  for (std::size_t i = 0; i < n; ++i) {
    kw[i] = b.k[i] * b.w[i];
    nw[i] = b.trials[i] * b.w[i];
    yBinary[i] = b.k[i] > 0.5 * b.trials[i] ? 1.0 : 0.0;
  }
  ext_rng* rngA = seededRng(7701u);
  ext_rng* rngB = seededRng(7701u);
  auto weighted = makeBinomialSampler(b.x.data(), b.k.data(), n, b.p,
                                      b.w.data(), b.trials.data(), &rngA);
  auto scaled = makeBinomialSampler(b.x.data(), kw.data(), n, b.p, nullptr,
                                    nw.data(), &rngB);
  CountDraws da(n, numSamples), db(n, numSamples);
  Results ra = da.results(), rb = db.results();
  ra.shape = rb.shape = nullptr;
  bool ran = !weighted->run(0, numSamples, ra) &&
             !scaled->run(0, numSamples, rb);
  check(ran && da.fits == db.fits,
        "a binomial row at weight w equals the row of k w successes in n w "
        "trials, bitwise");
  ext_rng_destroy(rngA);
  ext_rng_destroy(rngB);

  for (bool withWeights : {false, true}) {
    const double* w = withWeights ? b.w.data() : nullptr;
    ext_rng* rngC = seededRng(7702u);
    ext_rng* rngD = seededRng(7702u);
    auto oneTrial = makeBinomialSampler(b.x.data(), yBinary.data(), n, b.p, w,
                                        ones.data(), &rngC);
    auto binary = makeBinomialSampler(b.x.data(), yBinary.data(), n, b.p, w,
                                      nullptr, &rngD);
    CountDraws dc(n, numSamples), dd(n, numSamples);
    Results rc = dc.results(), rd = dd.results();
    rc.shape = rd.shape = nullptr;
    ran = !oneTrial->run(0, numSamples, rc) && !binary->run(0, numSamples, rd);
    check(ran && dc.fits == dd.fits && dc.logLikelihood == dd.logLikelihood,
          withWeights ? "a weighted binomial response of one trial a row is "
                        "the binary one, bitwise"
                      : "a binomial response of one trial a row is the "
                        "binary one, bitwise");
    ext_rng_destroy(rngC);
    ext_rng_destroy(rngD);
  }
}

/// The binomial row's log-likelihood against R's dbinom times the weight, to
/// 1e-10 relative, at trials up to one million: among them k = n / 2 and
/// k = n / 1000, where a log choose constant from std::lgamma is 1e-10 off.
/// A row of no trials reports 0 and a row the mask leaves out NaN.
void testBinomialLogLikelihood() {
  std::vector<double> k, trials, eta, w;
  auto add = [&](double size, double successes) {
    double rate = successes / size;
    std::vector<double> etas = {0.0, 0.7, -1.3};
    if (successes > 0.0 && successes < size)
      etas.push_back(std::log(rate / (1.0 - rate)));
    for (double e : etas)
      for (double weight : {1.0, 3.0}) {
        k.push_back(successes);
        trials.push_back(size);
        eta.push_back(e);
        w.push_back(weight);
      }
  };
  for (double size : {1.0, 7.0, 100.0, 1e4, 999983.0, 1e6})
    for (double successes : {0.0, 1.0, std::floor(size / 1000.0),
                             std::floor(size / 2.0), size - 1.0, size})
      if (successes >= 0.0 && successes <= size) add(size, successes);
  std::size_t n = k.size();
  LogisticResponse response(k.data(), nullptr, w.data(), n, trials.data());
  std::vector<double> out(n);
  response.computeLogLikelihood(eta.data(), 1.0, n, out.data());
  double worst = 0.0, worstHalf = 0.0, worstThousandth = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    double p = 1.0 / (1.0 + std::exp(-eta[i]));
    double q = 1.0 / (1.0 + std::exp(eta[i]));
    double expected = w[i] * Rf_dbinom_raw(k[i], trials[i], p, q, 1);
    double gap = std::fabs(out[i] - expected) /
                 std::max(std::fabs(expected), 1e-300);
    if (expected == 0.0) gap = std::fabs(out[i]);
    worst = std::max(worst, gap);
    if (trials[i] == 1e6 && k[i] == 5e5) worstHalf = std::max(worstHalf, gap);
    if (trials[i] == 1e6 && k[i] == 1e3)
      worstThousandth = std::max(worstThousandth, gap);
  }
  printf("  binomial log-likelihood against dbinom: %zu rows, worst %.3g, "
         "n = 1e6 at k = n/2 %.3g and k = n/1000 %.3g\n", n, worst, worstHalf,
         worstThousandth);
  check(worst < 1e-10 && worstHalf < 1e-10 && worstThousandth < 1e-10,
        "the binomial log-likelihood is w log dbinom(k; n, p) to 1e-10 relative");

  std::vector<double> k2 = {0.0, 2.0, 1.0}, n2 = {0.0, 5.0, 3.0},
                      eta2 = {0.4, 0.4, 0.4}, mask = {1.0, 1.0, 0.0}, out2(3);
  LogisticResponse small(k2.data(), nullptr, nullptr, 3, n2.data());
  small.setActiveRows(mask.data());
  small.computeLogLikelihood(eta2.data(), 1.0, 3, out2.data());
  check(out2[0] == 0.0 && std::isfinite(out2[1]) && std::isnan(out2[2]),
        "a binomial row of no trials reports 0, a masked one NaN");
}

/// A row of no trials and a binomial row of weight 0: finite throughout,
/// every fit, omega and working response, and the zero-trial row consumes
/// nothing and reads as a masked row of one trial, bitwise at the rows in.
void testBinomialZeroRows() {
  const std::size_t n = 40, numSamples = 300;
  BinomialFixture b(n, false);
  std::vector<double> kZero(b.k), nZero(b.trials), nOne(b.trials), mask(n, 1.0),
    wZero(b.w);
  for (std::size_t i = 3; i < n; i += 9) {
    kZero[i] = 0.0;
    nZero[i] = 0.0;
    nOne[i] = 1.0;
    mask[i] = 0.0;
  }
  wZero[5] = 0.0;
  ext_rng* rngA = seededRng(7801u);
  ext_rng* rngB = seededRng(7801u);
  auto zero = makeBinomialSampler(b.x.data(), kZero.data(), n, b.p,
                                  wZero.data(), nZero.data(), &rngA);
  auto masked = makeBinomialSampler(b.x.data(), kZero.data(), n, b.p,
                                    wZero.data(), nOne.data(), &rngB);
  masked->setActiveRows(mask.data());
  CountDraws da(n, numSamples), db(n, numSamples);
  Results ra = da.results(), rb = db.results();
  ra.shape = rb.shape = nullptr;
  bool ran = !zero->run(0, numSamples, ra) && !masked->run(0, numSamples, rb);
  bool finite = allFinite(da.fits) && allFinite(da.logLikelihood);
  const double* omega = zero->latents(0);
  const double* working = TestPeer::workingResponse(zero->chain(0));
  for (std::size_t i = 0; i < n; ++i)
    finite = finite && std::isfinite(omega[i]) && omega[i] > 0.0 &&
             std::isfinite(working[i]);
  bool zeroReportsZero = true;
  for (std::size_t s = 0; s < numSamples; ++s)
    for (std::size_t i = 0; i < n; ++i)
      if (nZero[i] == 0.0 || wZero[i] == 0.0)
        zeroReportsZero = zeroReportsZero && da.logLikelihood[s * n + i] == 0.0;
  check(ran && finite && zeroReportsZero,
        "a binomial row of no trials and one of weight 0 are finite "
        "throughout, each reporting a log-likelihood of 0");
  check(agreeOnRowsIn(da.fits, db.fits, mask, n) &&
          ext_rng_simulateContinuousUniform(rngA) ==
            ext_rng_simulateContinuousUniform(rngB),
        "a binomial row of no trials is a masked row of one, bitwise at the "
        "rows in, consuming the same variates");
  ext_rng_destroy(rngA);
  ext_rng_destroy(rngB);
}

/// setBinomialCounts draws nothing and leaves the working response
/// w (k - n / 2) / omega - offset under the new counts at every row, a row
/// of no trials and a masked one included; a refusal on a binary response;
/// and a state round trip of a binomial sampler, copied into a sampler
/// built on the same rows and reloaded into a fresh one, continues to
/// rounding.
void testBinomialSetCountsAndState() {
  const std::size_t n = 40, numSamples = 40;
  BinomialFixture b(n, false);
  std::vector<double> offset(n), k2(n), n2(n), mask(n, 1.0);
  for (std::size_t i = 0; i < n; ++i) {
    offset[i] = 0.2 * std::sin(static_cast<double>(i));
    n2[i] = i % 7 == 0 ? 0.0 : b.trials[i] + static_cast<double>(i % 3);
    k2[i] = std::min(n2[i], b.k[i] + 1.0);
  }
  mask[4] = 0.0;
  ext_rng* rngA = seededRng(7901u);
  auto sampler = makeBinomialSampler(b.x.data(), b.k.data(), n, b.p,
                                     b.w.data(), b.trials.data(), &rngA,
                                     offset.data());
  sampler->setActiveRows(mask.data());
  Results none;
  bool ran = !sampler->run(10, 0, none);
  SamplerStateData before;
  sampler->getState(before);
  bool installed = sampler->setBinomialCounts(k2.data(), n2.data());
  SamplerStateData after;
  sampler->getState(after);
  const double* omega = sampler->latents(0);
  const double* working = TestPeer::workingResponse(sampler->chain(0));
  bool rebuilt = true;
  for (std::size_t i = 0; i < n; ++i)
    rebuilt = rebuilt &&
      working[i] == b.w[i] * (k2[i] - 0.5 * n2[i]) / omega[i] - offset[i];
  check(ran && installed && rebuilt &&
          after.chains[0].rngState == before.chains[0].rngState &&
          after.chains[0].latents == before.chains[0].latents,
        "binomial setBinomialCounts draws nothing and rebuilds the working "
        "response from omega, the counts and the offset at every row");

  ext_rng* rngC = seededRng(7902u);
  auto binary = makeBinomialSampler(b.x.data(), mask.data(), n, b.p, nullptr,
                                    nullptr, &rngC);
  check(!binary->setBinomialCounts(k2.data(), n2.data()) &&
          !binary->holdsTrials() && sampler->holdsTrials(),
        "setBinomialCounts is refused on a binary response");
  check(!sampler->setData(b.x.data(), k2.data(), n, b.w.data(), offset.data(),
                          nullptr, 0),
        "setData is refused on a binomial response");

  ext_rng* rngD = seededRng(7903u);
  ext_rng* rngE = seededRng(7904u);
  auto copy = makeBinomialSampler(b.x.data(), k2.data(), n, b.p, b.w.data(),
                                  n2.data(), &rngD, offset.data());
  copy->setActiveRows(mask.data());
  auto reloaded = makeBinomialSampler(b.x.data(), k2.data(), n, b.p,
                                      b.w.data(), n2.data(), &rngE,
                                      offset.data());
  reloaded->setActiveRows(mask.data());
  ran = !sampler->run(5, 0, none);
  SamplerStateData state;
  sampler->getState(state);
  bool copied = restoresExactly(*copy, state, b.x.data());
  bool reloadedOk = restoresExactly(*reloaded, state, b.x.data());
  checkStructuralRoundTrip(state, *copy, "binomial state round trip");
  CountDraws d1(n, numSamples), d2(n, numSamples), d3(n, numSamples);
  Results r1 = d1.results(), r2 = d2.results(), r3 = d3.results();
  r1.shape = r2.shape = r3.shape = nullptr;
  ran = ran && !sampler->run(0, numSamples, r1) &&
        !copy->run(0, numSamples, r2) && !reloaded->run(0, numSamples, r3);
  check(ran && copied && reloadedOk && largestGap(d1.fits, d2.fits) < 1e-9 &&
          largestGap(d1.fits, d3.fits) < 1e-9 &&
          largestGap(d1.logLikelihood, d2.logLikelihood) < 1e-9,
        "a binomial sampler continues from its state, copied and reloaded");
  ext_rng_destroy(rngA);
  ext_rng_destroy(rngC);
  ext_rng_destroy(rngD);
  ext_rng_destroy(rngE);
}

/// The read-offs of a binomial row are its n binary rows' of its weight,
/// under either reading: a numeric basis's default scale, at creation and
/// re-read by a rescaling basis change, and a linear leaf's standardization
/// equal the expansion's to rounding, a row of no trials out of each.
void testBinomialReadOffs() {
  const std::size_t n = 40;
  BinomialFixture b(n, false);
  for (std::size_t i = 1; i < n; i += 8) {
    b.trials[i] = 0.0;
    b.k[i] = 0.0;
  }
  Fixture f(n);
  f.x = b.x;
  std::vector<double> counts(n), dose(n), dose2(n);
  for (std::size_t i = 0; i < n; ++i) {
    counts[i] = b.trials[i] * b.w[i];
    dose[i] = 0.5 + runif01() * (i % 3 == 0 ? 3.0 : 1.0);
    dose2[i] = runif01() * (i % 4 == 0 ? 5.0 : 1.0);
  }
  // a row of no trials far out on both columns, which only its absence
  // leaves the scale untouched by
  dose[1] = dose2[1] = 1.0e3;
  Expansion e(f, counts, &dose);
  std::vector<double> yE(e.n), dose2E = expanded(e, dose2);
  for (std::size_t r = 0; r < e.n; ++r) yE[r] = r % 2 == 0 ? 1.0 : 0.0;
  auto scalesOf = [&](const double* x, const double* y, std::size_t rows,
                      const double* w, const double* trials,
                      const double* basis, const double* basis2,
                      bool frequency) {
    AmplitudeSpec spec;
    spec.family = ResponseFamily::logistic;
    spec.forests.resize(2);
    spec.forests[1].kind = ForestKind::numeric;
    spec.forests[1].basis = basis;
    SamplerOptions options;
    if (frequency) options = frequencyOptions(options);
    options.trials = trials;
    ext_rng* rng = seededRng(7711u);
    ConstantLeafSampler sampler(x, y, rows, f.p, w, nullptr, 1.0, sigmaDf,
                                sigmaRawScale, options, spec, &rng);
    double created = sampler.forestCalibration(0, 1).basisScale;
    bool rescaled = sampler.setForestBasisRescaled(1, basis2, 1) ==
                    BasisRescaleResult::rescaled;
    double changed = sampler.forestCalibration(0, 1).basisScale;
    ext_rng_destroy(rng);
    return std::array<double, 3>{created, changed, rescaled ? 1.0 : 0.0};
  };
  auto expansion = scalesOf(e.x.data(), yE.data(), e.n, nullptr, nullptr,
                            e.extra.data(), dose2E.data(), false);
  bool basisAgrees = expansion[2] == 1.0;
  for (bool frequency : {false, true}) {
    auto binomial = scalesOf(b.x.data(), b.k.data(), n, b.w.data(),
                             b.trials.data(), dose.data(), dose2.data(),
                             frequency);
    basisAgrees = basisAgrees && binomial[2] == 1.0 &&
                  relativeGap(binomial[0], expansion[0]) < 1e-14 &&
                  relativeGap(binomial[1], expansion[1]) < 1e-14;
  }
  check(basisAgrees,
        "a binomial basis scale is its binary rows', at creation and at a "
        "rescaling basis change, under either reading");

  const std::size_t columns[] = {0};
  SamplerOptions options;
  options.numTrees = 5;
  options.leafCovariateColumns = columns;
  options.numLeafCovariates = 1;
  std::vector<double> xFar(b.x);
  xFar[1] = 50.0;  // the row of no trials, out
  Fixture g(n);
  g.x = xFar;
  Expansion eLeaf(g, counts);
  ext_rng* rngB = seededRng(7712u);
  ext_rng* rngE = seededRng(7712u);
  SamplerOptions binomialOptions = options;
  binomialOptions.trials = b.trials.data();
  Sampler<LinearGaussianLeaf> binomial(xFar.data(), b.k.data(), n, b.p,
                                       b.w.data(), nullptr,
                                       ResponseFamily::logistic, 1.0, sigmaDf,
                                       1.0, binomialOptions, &rngB);
  Sampler<LinearGaussianLeaf> binary(eLeaf.x.data(), yE.data(), eLeaf.n, b.p,
                                     nullptr, nullptr,
                                     ResponseFamily::logistic, 1.0, sigmaDf,
                                     1.0, options, &rngE);
  const auto& lb = binomial.chain(0).leaf();
  const auto& le = binary.chain(0).leaf();
  check(relativeGap(lb.covariateMeans()[0], le.covariateMeans()[0]) < 1e-14 &&
          relativeGap(lb.covariateSds()[0], le.covariateSds()[0]) < 1e-14,
        "a binomial linear leaf's standardization is its binary rows'");
  ext_rng_destroy(rngB);
  ext_rng_destroy(rngE);
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
  testDoorsAgainstExpansion();
  testLogisticBasisScale();
  testLeafFallbackWeighted();
  testStudentSetWeights();
  testRowOutIsInert();
  testGPLengthscale();
  testGPLengthscaleMaskedAndReal();
  testGaussianFallbackWeighted();
  testAnchorTotalAtMostOne();
  testNBinomAgainstExpansion();
  testNBinomZeroWeights();
  testNBinomSetWeightsAndState();
  testMultinomialAgainstExpansion(false);
  testMultinomialAgainstExpansion(true);
  testMultinomialZeroWeightsAndState();
  testBinomialAgainstSplitRows();
  testBinomialScaledAndBinary();
  testBinomialLogLikelihood();
  testBinomialZeroRows();
  testBinomialSetCountsAndState();
  testBinomialReadOffs();
  rngState = savedState;
}
