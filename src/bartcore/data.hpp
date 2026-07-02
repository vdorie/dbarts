#ifndef BARTCORE_DATA_HPP
#define BARTCORE_DATA_HPP

#include <cstddef>
#include <cstdint>
#include <vector>

namespace bartcore {

using std::size_t;

/// Predictor code type; per-column widths are planned (see
/// docs/design/core-generalization.md), uint16_t matches the classic engine.
using xint_t = std::uint16_t;

/// Classic dense column store: borrowed column-major doubles quantized once
/// into per-column integer codes against uniformly spaced cut points.
struct ColumnStore {
  size_t numObservations = 0;
  size_t numPredictors = 0;
  const double* x = nullptr;  // borrowed, column-major

  std::vector<xint_t> codes;  // column-major, numObservations x numPredictors
  std::vector<std::vector<double>> cutPoints;
  std::vector<std::uint32_t> numCuts;

  size_t numTestObservations = 0;
  std::vector<xint_t> testCodes;  // row-major, numTestObservations x numPredictors

  // Codes are k such that cutPoints[k - 1] < value <= cutPoints[k], with
  // value > all cuts mapping to numCuts (always right of any split).
  xint_t codeFor(size_t variable, double value) const {
    const std::vector<double>& cuts = cutPoints[variable];
    std::uint32_t k = 0;
    while (k < numCuts[variable] && value > cuts[k]) ++k;
    return static_cast<xint_t>(k);
  }

  void build(const double* x_, size_t n, size_t p, std::uint32_t maxNumCuts) {
    x = x_;
    numObservations = n;
    numPredictors = p;
    cutPoints.resize(p);
    numCuts.resize(p);
    codes.resize(n * p);

    for (size_t j = 0; j < p; ++j) {
      const double* column = x + j * n;
      double xMin = column[0], xMax = column[0];
      for (size_t i = 1; i < n; ++i) {
        if (column[i] < xMin) xMin = column[i];
        if (column[i] > xMax) xMax = column[i];
      }

      numCuts[j] = maxNumCuts;
      cutPoints[j].resize(maxNumCuts);
      double increment = (xMax - xMin) / static_cast<double>(maxNumCuts + 1);
      for (std::uint32_t k = 0; k < maxNumCuts; ++k)
        cutPoints[j][k] = xMin + static_cast<double>(k + 1) * increment;

      for (size_t i = 0; i < n; ++i)
        codes[i + j * n] = codeFor(j, column[i]);
    }
  }

  void buildTest(const double* x_test, size_t numTest) {
    numTestObservations = numTest;
    testCodes.resize(numTest * numPredictors);
    for (size_t i = 0; i < numTest; ++i)
      for (size_t j = 0; j < numPredictors; ++j)
        testCodes[i * numPredictors + j] =
          codeFor(j, x_test[i + j * numTest]);
  }

  const xint_t* column(size_t variable) const {
    return codes.data() + variable * numObservations;
  }
  const xint_t* testRow(size_t testObservation) const {
    return testCodes.data() + testObservation * numPredictors;
  }
};

}  // namespace bartcore

#endif  // BARTCORE_DATA_HPP
