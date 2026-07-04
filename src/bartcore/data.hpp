#ifndef BARTCORE_DATA_HPP
#define BARTCORE_DATA_HPP

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <vector>

namespace bartcore {

using std::size_t;

/// Predictor code type; per-column widths are planned (see
/// docs/design/core-generalization.md), uint16_t matches the classic engine.
using xint_t = std::uint16_t;

/// Reserved code for a missing ordinal value; cut counts cap below it so
/// real codes (which reach numCuts) never collide.
constexpr xint_t naCode = 0xFFFFu;
constexpr std::uint32_t maxNumCutsRepresentable = 0xFFFEu - 1u;

/// Missing categorical values take a category position above the real ones
/// so the reachable-mask machinery routes them like any other category: the
/// fixed position 63 (the top of the rule word) for columns whose mask is
/// inline, position K for pooled columns of K > 63 categories.
constexpr std::uint32_t naCategory = 63;

/// Words a categorical mask occupies: category bits 0..K-1 plus the
/// missing-direction position (63 inline, K pooled), ceil((K + 1) / 64).
constexpr size_t maskWordsForCount(std::uint32_t numCategories) {
  return (static_cast<size_t>(numCategories) + 64) / 64;
}

/// The reserved code a missing value takes in a categorical column.
constexpr xint_t missingCategoryCode(std::uint32_t numCategories) {
  return maskWordsForCount(numCategories) > 1
    ? static_cast<xint_t>(numCategories)
    : static_cast<xint_t>(naCategory);
}

inline bool isNA(double value) { return value != value; }

/// Column types the store distinguishes. Ordinal columns quantize against
/// cut points and split by threshold; categorical columns hold integer
/// category codes 0..numCategories-1 directly and split by subset. Masks of
/// up to 63 categories live inline in the rule word; wider columns pool
/// their mask words per tree, and the flattened format moves masks past 53
/// categories (the double-exactness bound) into a side channel (see
/// docs/design/pooled-masks.md). Codes must fit xint_t, including the
/// reserved missing code K of a pooled column.
enum class ColumnType : std::uint8_t { ordinal, categorical };

constexpr std::uint32_t maxCategories = 0xFFFFu;

/// The widest category count whose full direction mask a double represents
/// exactly; the flattened format's value-encoding boundary.
constexpr std::uint32_t maxValueEncodableCategories = 53;

/// Classic dense column store: borrowed column-major doubles quantized once
/// into per-column integer codes against per-column cut points, either
/// uniformly spaced over the column's range or at unique-value midpoints
/// (quantile mode). numCuts is fixed once built; recomputing cuts for new
/// values keeps the count so existing split indices stay in range. For
/// categorical columns numCuts holds the (fixed) category count, cutPoints
/// stays empty, and codes are the values themselves.
struct ColumnStore {
  size_t numObservations = 0;
  size_t numPredictors = 0;
  const double* x = nullptr;  // borrowed, column-major
  bool useQuantiles = false;

  std::vector<ColumnType> types;
  std::vector<xint_t> codes;  // column-major, numObservations x numPredictors
  std::vector<std::vector<double>> cutPoints;
  std::vector<std::uint32_t> numCuts;
  std::vector<std::uint32_t> maxNumCuts;  // cap on quantile-induced counts
  // per column, whether any training value is missing; gates the extra
  // missing-direction draw in rules and the NA-aware partition kernel, so
  // NA-free columns keep today's draws and code paths exactly
  std::vector<std::uint8_t> hasMissing;
  // categorical tier flags, fixed once category counts are (they never
  // change after build): wide columns (> 53 levels) need the flattened
  // format's mask side channel, pooled ones (> 63) the per-tree mask pool
  bool hasWideCategorical = false;
  bool hasPooledCategorical = false;

  /// Whether column j's flattened masks move to the side channel.
  bool columnHasWideMask(size_t j) const {
    return types[j] == ColumnType::categorical &&
           numCuts[j] > maxValueEncodableCategories;
  }
  /// Whether column j's rules store pool offsets instead of inline masks.
  bool columnIsPooled(size_t j) const {
    return types[j] == ColumnType::categorical &&
           maskWordsForCount(numCuts[j]) > 1;
  }

  void refreshCategoricalTiers() {
    hasWideCategorical = hasPooledCategorical = false;
    for (size_t j = 0; j < numPredictors; ++j) {
      if (columnHasWideMask(j)) hasWideCategorical = true;
      if (columnIsPooled(j)) hasPooledCategorical = true;
    }
  }

  size_t numTestObservations = 0;
  const double* x_test = nullptr;  // borrowed, column-major
  std::vector<xint_t> testCodes;  // row-major, numTestObservations x numPredictors
  // borrowed; added to recorded test fits. buildTest leaves it alone (the
  // caller keeps the lengths consistent), clearTest clears it.
  const double* testOffset = nullptr;

  // Ordinal codes are k such that cutPoints[k - 1] < value <= cutPoints[k],
  // with value > all cuts mapping to numCuts (always right of any split);
  // categorical codes are the values themselves. Missing values take the
  // reserved codes.
  xint_t codeFor(size_t variable, double value) const {
    if (types[variable] == ColumnType::categorical)
      return isNA(value) ? missingCategoryCode(numCuts[variable])
                         : static_cast<xint_t>(value);
    if (isNA(value)) return naCode;
    const std::vector<double>& cuts = cutPoints[variable];
    std::uint32_t k = 0;
    while (k < numCuts[variable] && value > cuts[k]) ++k;
    return static_cast<xint_t>(k);
  }

  /// A categorical value is representable when it is an integral code of an
  /// existing category or missing; the category count is fixed once built.
  bool categoricalValueIsValid(size_t variable, double value) const {
    if (isNA(value)) return true;
    return value >= 0.0 && value < static_cast<double>(numCuts[variable]) &&
           value == static_cast<double>(static_cast<xint_t>(value));
  }

  /// Quantile-mode support: sorted unique values of a column and the
  /// stepping that thins their midpoints to at most maxNumCuts[j] cuts.
  struct QuantileGrid {
    std::vector<double> sortedUnique;
    std::uint32_t inducedNumCuts = 0;
    size_t step = 1, offset = 0;
  };

  QuantileGrid quantileGridForColumn(size_t j, const double* values) const {
    QuantileGrid grid;
    grid.sortedUnique.reserve(numObservations);
    // NaN would break the sort's ordering; the grid is over observed values
    for (size_t i = 0; i < numObservations; ++i)
      if (!isNA(values[i])) grid.sortedUnique.push_back(values[i]);
    std::sort(grid.sortedUnique.begin(), grid.sortedUnique.end());
    grid.sortedUnique.erase(
      std::unique(grid.sortedUnique.begin(), grid.sortedUnique.end()),
      grid.sortedUnique.end());

    size_t numUnique = grid.sortedUnique.size();
    if (numUnique <= 1) {  // constant or fully missing column
      grid.inducedNumCuts = 0;
    } else if (numUnique <= static_cast<size_t>(maxNumCuts[j]) + 1) {
      grid.inducedNumCuts = static_cast<std::uint32_t>(numUnique - 1);
    } else {
      grid.inducedNumCuts = maxNumCuts[j];
      grid.step = numUnique / grid.inducedNumCuts;
      grid.offset = grid.step / 2;
    }
    return grid;
  }

  void fillCutsFromQuantileGrid(size_t j, const QuantileGrid& grid) {
    cutPoints[j].resize(numCuts[j]);
    for (std::uint32_t k = 0; k < numCuts[j]; ++k) {
      size_t index = std::min(static_cast<size_t>(k) * grid.step + grid.offset,
                              grid.sortedUnique.size() - 2);
      cutPoints[j][k] =
        0.5 * (grid.sortedUnique[index] + grid.sortedUnique[index + 1]);
    }
  }

  void fillCutsUniformly(size_t j) {
    const double* column = x + j * numObservations;
    // the range is over observed values; NaN never satisfies a comparison,
    // so only the seed needs the explicit skip
    double xMin = 0.0, xMax = 0.0;
    size_t i = 0;
    while (i < numObservations && isNA(column[i])) ++i;
    if (i < numObservations) {
      xMin = xMax = column[i];
      for (++i; i < numObservations; ++i) {
        if (column[i] < xMin) xMin = column[i];
        if (column[i] > xMax) xMax = column[i];
      }
    }
    cutPoints[j].resize(numCuts[j]);
    double increment = (xMax - xMin) / static_cast<double>(numCuts[j] + 1);
    for (std::uint32_t k = 0; k < numCuts[j]; ++k)
      cutPoints[j][k] = xMin + static_cast<double>(k + 1) * increment;
  }

  /// Initial cut construction; sets numCuts[j]. A categorical column takes
  /// its (fixed) category count from its largest value and keeps no cuts.
  void buildCutsForColumn(size_t j) {
    if (types[j] == ColumnType::categorical) {
      const double* column = x + j * numObservations;
      double maxValue = -1.0;
      for (size_t i = 0; i < numObservations; ++i)
        if (!isNA(column[i]) && column[i] > maxValue) maxValue = column[i];
      numCuts[j] = maxValue < 0.0
        ? 0 : static_cast<std::uint32_t>(maxValue) + 1;
      cutPoints[j].clear();
    } else if (useQuantiles) {
      QuantileGrid grid = quantileGridForColumn(j, x + j * numObservations);
      numCuts[j] = grid.inducedNumCuts;
      fillCutsFromQuantileGrid(j, grid);
    } else {
      numCuts[j] = maxNumCuts[j];
      fillCutsUniformly(j);
    }
  }

  /// Recompute cuts for a column's current values, keeping numCuts[j] fixed.
  /// In quantile mode fewer induced cuts than existing would leave splits
  /// invalid: returns false, having changed nothing (extra induced cuts are
  /// silently thinned, as in the reference engine). Categorical columns have
  /// nothing to refresh; the caller pre-checked value validity.
  bool refreshCutsForColumn(size_t j) {
    if (types[j] == ColumnType::categorical) return true;
    if (useQuantiles) {
      QuantileGrid grid = quantileGridForColumn(j, x + j * numObservations);
      if (grid.inducedNumCuts < numCuts[j]) return false;
      fillCutsFromQuantileGrid(j, grid);
    } else {
      fillCutsUniformly(j);
    }
    return true;
  }

  /// Non-mutating feasibility check for values not yet installed: quantile
  /// refresh feasibility for ordinal columns, category-code validity for
  /// categorical ones.
  bool cutsWouldRemainValid(size_t j, const double* values) const {
    if (types[j] == ColumnType::categorical) {
      for (size_t i = 0; i < numObservations; ++i)
        if (!categoricalValueIsValid(j, values[i])) return false;
      return true;
    }
    if (!useQuantiles) return true;
    return quantileGridForColumn(j, values).inducedNumCuts >= numCuts[j];
  }

  /// Install externally chosen cut points (ascending) for a column; the cut
  /// count may shrink or grow, and existing splits beyond the new range are
  /// the caller's problem (the sampler collapses them).
  void setCutPointsForColumn(size_t j, const double* cuts,
                             std::uint32_t numCutPoints) {
    cutPoints[j].assign(cuts, cuts + numCutPoints);
    numCuts[j] = numCutPoints;
    if (maxNumCuts[j] < numCutPoints) maxNumCuts[j] = numCutPoints;
    quantizeColumn(j);
    if (numTestObservations > 0) quantizeTestColumn(j);
  }

  void quantizeColumn(size_t j) {
    const double* column = x + j * numObservations;
    std::uint8_t anyMissing = 0;
    for (size_t i = 0; i < numObservations; ++i) {
      xint_t code = codeFor(j, column[i]);
      if (isNA(column[i])) anyMissing = 1;
      codes[i + j * numObservations] = code;
    }
    hasMissing[j] = anyMissing;
  }

  void quantizeTestColumn(size_t j) {
    for (size_t i = 0; i < numTestObservations; ++i)
      testCodes[i * numPredictors + j] =
        codeFor(j, x_test[i + j * numTestObservations]);
  }

  /// columnTypes may be null for all-ordinal. Categorical columns must hold
  /// integral values 0..K-1 with K <= maxCategories; the caller validates.
  void build(const double* x_, size_t n, size_t p,
             const std::uint32_t* maxNumCuts_, bool useQuantiles_,
             const ColumnType* columnTypes = nullptr) {
    x = x_;
    numObservations = n;
    numPredictors = p;
    useQuantiles = useQuantiles_;
    if (columnTypes != nullptr) {
      types.assign(columnTypes, columnTypes + p);
    } else {
      types.assign(p, ColumnType::ordinal);
    }
    cutPoints.resize(p);
    numCuts.resize(p);
    maxNumCuts.assign(maxNumCuts_, maxNumCuts_ + p);
    // keep the reserved missing code out of the real code range
    for (size_t j = 0; j < p; ++j)
      if (maxNumCuts[j] > maxNumCutsRepresentable)
        maxNumCuts[j] = maxNumCutsRepresentable;
    codes.resize(n * p);
    hasMissing.assign(p, 0);

    for (size_t j = 0; j < p; ++j) {
      buildCutsForColumn(j);
      quantizeColumn(j);
    }
    refreshCategoricalTiers();
  }

  void build(const double* x_, size_t n, size_t p, std::uint32_t maxNumCuts_,
             bool useQuantiles_ = false,
             const ColumnType* columnTypes = nullptr) {
    std::vector<std::uint32_t> maxPerColumn(p, maxNumCuts_);
    build(x_, n, p, maxPerColumn.data(), useQuantiles_, columnTypes);
  }

  void buildTest(const double* x_test_, size_t numTest) {
    x_test = x_test_;
    numTestObservations = numTest;
    testCodes.resize(numTest * numPredictors);
    for (size_t j = 0; j < numPredictors; ++j) quantizeTestColumn(j);
  }

  /// A row-subset view of a built parent store: copies the parent's cut
  /// structure and gathers the subset's codes, so the view bins identically
  /// to the parent by construction; testRows also index the parent's
  /// observations. Views hold no raw values (x and x_test stay null), which
  /// rules out every raw-x path here (quantizeColumn and the mutation
  /// surface); callers enforce that. The view is self-contained: nothing
  /// references the parent after this returns.
  void buildFromParent(const ColumnStore& parent, const size_t* rows,
                       size_t numRows, const size_t* testRows,
                       size_t numTestRows) {
    x = nullptr;
    numObservations = numRows;
    numPredictors = parent.numPredictors;
    useQuantiles = parent.useQuantiles;
    types = parent.types;
    cutPoints = parent.cutPoints;
    numCuts = parent.numCuts;
    maxNumCuts = parent.maxNumCuts;
    codes.resize(numRows * numPredictors);
    hasMissing.assign(numPredictors, 0);
    refreshCategoricalTiers();
    for (size_t j = 0; j < numPredictors; ++j) {
      const xint_t* parentColumn =
        parent.codes.data() + j * parent.numObservations;
      xint_t missingCode = types[j] == ColumnType::categorical
        ? missingCategoryCode(numCuts[j]) : naCode;
      xint_t* column = codes.data() + j * numRows;
      for (size_t i = 0; i < numRows; ++i) {
        column[i] = parentColumn[rows[i]];
        if (column[i] == missingCode) hasMissing[j] = 1;
      }
    }
    x_test = nullptr;
    numTestObservations = numTestRows;
    testCodes.resize(numTestRows * numPredictors);
    for (size_t i = 0; i < numTestRows; ++i)
      for (size_t j = 0; j < numPredictors; ++j)
        testCodes[i * numPredictors + j] =
          parent.codes[testRows[i] + j * parent.numObservations];
    testOffset = nullptr;
  }

  // Mutation. Snapshot/rollback of x, cutPoints, and codes is the caller's
  // (the sampler's) responsibility; these only install new values. Cut
  // refreshes assume the caller pre-checked quantile feasibility with
  // cutsWouldRemainValid.

  /// Replace the whole predictor matrix by pointer swap.
  void setPredictors(const double* newX, bool updateCuts) {
    x = newX;
    for (size_t j = 0; j < numPredictors; ++j) {
      if (updateCuts) refreshCutsForColumn(j);
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
      if (updateCuts) refreshCutsForColumn(j);
      quantizeColumn(j);
    }
  }

  /// Overwrite a single cell in place, re-quantizing against existing cuts.
  /// A missing value marks the column; the flag only clears on a full
  /// column re-quantize (conservative but never wrong - the NA-aware
  /// partition handles columns without missing values too).
  void setCell(size_t i, size_t j, double value) {
    const_cast<double*>(x)[i + j * numObservations] = value;
    codes[i + j * numObservations] = codeFor(j, value);
    if (isNA(value)) hasMissing[j] = 1;
  }

  /// Whole-data replacement: new values for the same predictors, possibly a
  /// new number of observations. Ordinal cuts are rebuilt from scratch, so
  /// unlike refreshCutsForColumn a quantile-mode count may shrink and the
  /// caller remaps existing splits onto the new grid; categorical category
  /// counts stay fixed (the caller validates the new values).
  void setData(const double* x_, size_t n) {
    x = x_;
    numObservations = n;
    codes.resize(n * numPredictors);
    for (size_t j = 0; j < numPredictors; ++j) {
      if (types[j] != ColumnType::categorical) buildCutsForColumn(j);
      quantizeColumn(j);
    }
  }

  void clearTest() {
    x_test = nullptr;
    numTestObservations = 0;
    testCodes.clear();
    testOffset = nullptr;
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
