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
  rngState = savedState;
}
