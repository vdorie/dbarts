#ifndef BARTCORE_DATA_HPP
#define BARTCORE_DATA_HPP

#include <cstddef>
#include <cstdint>
#include <cstring>
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
  const double* x_test = nullptr;  // borrowed, column-major
  std::vector<xint_t> testCodes;  // row-major, numTestObservations x numPredictors

  // Codes are k such that cutPoints[k - 1] < value <= cutPoints[k], with
  // value > all cuts mapping to numCuts (always right of any split).
  xint_t codeFor(size_t variable, double value) const {
    const std::vector<double>& cuts = cutPoints[variable];
    std::uint32_t k = 0;
    while (k < numCuts[variable] && value > cuts[k]) ++k;
    return static_cast<xint_t>(k);
  }

  void updateCutsForColumn(size_t j) {
    const double* column = x + j * numObservations;
    double xMin = column[0], xMax = column[0];
    for (size_t i = 1; i < numObservations; ++i) {
      if (column[i] < xMin) xMin = column[i];
      if (column[i] > xMax) xMax = column[i];
    }
    std::uint32_t maxNumCuts = numCuts[j];
    cutPoints[j].resize(maxNumCuts);
    double increment = (xMax - xMin) / static_cast<double>(maxNumCuts + 1);
    for (std::uint32_t k = 0; k < maxNumCuts; ++k)
      cutPoints[j][k] = xMin + static_cast<double>(k + 1) * increment;
  }

  void quantizeColumn(size_t j) {
    const double* column = x + j * numObservations;
    for (size_t i = 0; i < numObservations; ++i)
      codes[i + j * numObservations] = codeFor(j, column[i]);
  }

  void quantizeTestColumn(size_t j) {
    for (size_t i = 0; i < numTestObservations; ++i)
      testCodes[i * numPredictors + j] =
        codeFor(j, x_test[i + j * numTestObservations]);
  }

  void build(const double* x_, size_t n, size_t p, std::uint32_t maxNumCuts) {
    x = x_;
    numObservations = n;
    numPredictors = p;
    cutPoints.resize(p);
    numCuts.assign(p, maxNumCuts);
    codes.resize(n * p);

    for (size_t j = 0; j < p; ++j) {
      updateCutsForColumn(j);
      quantizeColumn(j);
    }
  }

  void buildTest(const double* x_test_, size_t numTest) {
    x_test = x_test_;
    numTestObservations = numTest;
    testCodes.resize(numTest * numPredictors);
    for (size_t j = 0; j < numPredictors; ++j) quantizeTestColumn(j);
  }

  // Mutation. Snapshot/rollback of x, cutPoints, and codes is the caller's
  // (the sampler's) responsibility; these only install new values.

  /// Replace the whole predictor matrix by pointer swap.
  void setPredictors(const double* newX, bool updateCuts) {
    x = newX;
    for (size_t j = 0; j < numPredictors; ++j) {
      if (updateCuts) updateCutsForColumn(j);
      quantizeColumn(j);
    }
  }

  /// Overwrite a subset of columns in place; newColumns is column-major,
  /// numObservations x numColumns.
  void setColumns(const double* newColumns, const size_t* columns,
                  size_t numColumns, bool updateCuts) {
    double* x_mutable = const_cast<double*>(x);
    for (size_t k = 0; k < numColumns; ++k) {
      size_t j = columns[k];
      std::memcpy(x_mutable + j * numObservations,
                  newColumns + k * numObservations,
                  numObservations * sizeof(double));
      if (updateCuts) updateCutsForColumn(j);
      quantizeColumn(j);
    }
  }

  /// Overwrite a single cell in place, re-quantizing against existing cuts.
  void setCell(size_t i, size_t j, double value) {
    const_cast<double*>(x)[i + j * numObservations] = value;
    codes[i + j * numObservations] = codeFor(j, value);
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
