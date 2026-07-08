#ifndef BARTCORE_DATA_HPP
#define BARTCORE_DATA_HPP

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <vector>

namespace bartcore {

using std::size_t;

/// Predictor code type: the per-cell quantized ordinal rank or categorical
/// code. uint16_t caps a column at maxNumCutsRepresentable cuts, keeping every
/// real code below naCode (0xFFFF), the reserved missing marker.
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

/// Mean and sample sd of the observed (non-missing) values of a raw column:
/// the leaf-covariate standardization constants. A constant (or all-missing)
/// column keeps sd 1. One definition shared by LinearGaussianLeaf and the
/// view gather in buildFromParent, so a full-rows view standardizes
/// bit-identically to a sampler over the raw data.
inline void standardizationMomentsForColumn(const double* column, size_t n,
                                            double* mean_, double* sd_) {
  double total = 0.0;
  size_t numObserved = 0;
  for (size_t i = 0; i < n; ++i) {
    if (isNA(column[i])) continue;
    total += column[i];
    ++numObserved;
  }
  double mean = numObserved > 0 ? total / static_cast<double>(numObserved)
                                : 0.0;
  double sumOfSquares = 0.0;
  for (size_t i = 0; i < n; ++i)
    if (!isNA(column[i]))
      sumOfSquares += (column[i] - mean) * (column[i] - mean);
  double sd = numObserved > 1
    ? std::sqrt(sumOfSquares / static_cast<double>(numObserved - 1))
    : 0.0;
  if (!(sd > 0.0)) sd = 1.0;
  *mean_ = mean;
  *sd_ = sd;
}

/// Column types the store distinguishes. Ordinal columns quantize against
/// cut points and split by threshold; categorical columns hold integer
/// category codes 0..numCategories-1 directly and split by subset. Masks of
/// up to 63 categories live inline in the rule word (and inline in the
/// flattened node); wider columns pool their mask words per tree and the
/// flattened format references them through a side channel (see
/// docs/design/pooled-masks.md). Codes must fit xint_t, including the
/// reserved missing code K of a pooled column.
enum class ColumnType : std::uint8_t { ordinal, categorical };

constexpr std::uint32_t maxCategories = 0xFFFFu;

/// CSC-built columns at or below this nonzero fraction take rank-bitmap
/// storage; denser ones densify their codes (docs/design/sparse-columns.md).
constexpr double sparseDensityThreshold = 0.2;

/// Rank-bitmap hot storage of a sparse ordinal column: code(i) is zeroCode
/// when bit i is clear, otherwise the packed nonzero code at rank(i), an
/// O(1) word rank plus masked popcount. The pattern (bits, wordRanks) is
/// fixed at build; re-quantization rewrites nzCodes and zeroCode only.
struct SparseColumnData {
  std::vector<std::uint64_t> bits;       // ceil(n / 64) words
  std::vector<std::uint32_t> wordRanks;  // nonzeros before each word
  std::vector<xint_t> nzCodes;           // ascending-row order
  xint_t zeroCode = 0;

  xint_t at(size_t i) const {
    std::uint64_t word = bits[i >> 6];
    std::uint64_t bit = std::uint64_t{1} << (i & 63u);
    if ((word & bit) == 0) return zeroCode;
    return nzCodes[wordRanks[i >> 6] +
                   static_cast<size_t>(std::popcount(word & (bit - 1u)))];
  }
};

/// A CSC column's borrowed slices, kept for re-quantization; rows within a
/// column are unique and in range (the host validates before building).
struct CscColumnSlice {
  const double* values = nullptr;
  const int* rows = nullptr;
  size_t numNonzero = 0;
};

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
  // packed dense codes; per-column starts in codeOffsets (j * numObservations
  // for dense-matrix builds, packed among densified columns for CSC builds,
  // unused for rank-stored columns)
  std::vector<xint_t> codes;
  std::vector<size_t> codeOffsets;
  // per column, the rank-storage slot in sparseColumns or -1 for dense
  std::vector<std::int32_t> sparseSlot;
  std::vector<SparseColumnData> sparseColumns;
  // borrowed CSC slices per column when builtFromCsc, for re-quantization
  std::vector<CscColumnSlice> cscSlices;
  // per column, borrowed raw dense values on mixed builds (null for
  // CSC-backed columns); empty for every other build kind
  std::vector<const double*> mixedRawColumns;
  bool builtFromCsc = false;
  bool hasSparse = false;
  std::vector<std::vector<double>> cutPoints;
  std::vector<std::uint32_t> numCuts;
  std::vector<std::uint32_t> maxNumCuts;  // cap on quantile-induced counts
  // per column, whether any training value is missing; gates the extra
  // missing-direction draw in rules and the NA-aware partition kernel, so
  // NA-free columns keep today's draws and code paths exactly
  std::vector<std::uint8_t> hasMissing;
  // categorical tier flag, fixed once category counts are (they never change
  // after build): pooled columns (> 63 levels) need the per-tree mask pool
  // and the flattened format's mask side channel; narrower masks are inline
  bool hasPooledCategorical = false;

  /// Whether column j's rules store pool offsets instead of inline masks
  /// (more than 63 categories, so the mask spans more than one word).
  bool columnIsPooled(size_t j) const {
    return types[j] == ColumnType::categorical &&
           maskWordsForCount(numCuts[j]) > 1;
  }

  void refreshCategoricalTiers() {
    hasPooledCategorical = false;
    for (size_t j = 0; j < numPredictors; ++j)
      if (columnIsPooled(j)) hasPooledCategorical = true;
  }

  size_t numTestObservations = 0;
  const double* x_test = nullptr;  // borrowed, column-major
  std::vector<xint_t> testCodes;  // row-major, numTestObservations x numPredictors
  // borrowed; added to recorded test fits. buildTest leaves it alone (the
  // caller keeps the lengths consistent), clearTest clears it.
  const double* testOffset = nullptr;

  // Raw values of designated columns, gathered by buildFromParent so a view
  // can serve a leaf model that reads raw predictors despite holding no x.
  // The standardization constants come from the parent's full training data:
  // the same calibration inheritance as the copied cut grid, so every fold
  // runs the prior a full-data fit would. Owned, so views stay
  // self-contained.
  std::vector<size_t> gatheredRawColumns;
  std::vector<double> gatheredRawValues;      // column-major, numObservations x q
  std::vector<double> gatheredRawTestValues;  // column-major, numTestObservations x q
  std::vector<double> gatheredMeans;
  std::vector<double> gatheredSds;

  /// Borrowed raw values of column j when a dense source backs it: the
  /// matrix column of a dense build, the dense-slice column of a mixed
  /// build, null for CSC-backed columns and views.
  const double* denseSourceColumn(size_t j) const {
    if (x != nullptr) return x + j * numObservations;
    if (!mixedRawColumns.empty()) return mixedRawColumns[j];
    return nullptr;
  }

  /// Whether column j quantizes from a borrowed CSC slice (every column of
  /// a buildFromCsc store, the sparse-source columns of a mixed build).
  bool columnIsCscBacked(size_t j) const {
    return builtFromCsc &&
           (mixedRawColumns.empty() || mixedRawColumns[j] == nullptr);
  }

  /// Raw training values of column j: the dense source when one backs it,
  /// a gathered copy on views, null when neither serves it.
  const double* rawColumn(size_t j) const {
    const double* raw = denseSourceColumn(j);
    if (raw != nullptr) return raw;
    for (size_t k = 0; k < gatheredRawColumns.size(); ++k)
      if (gatheredRawColumns[k] == j)
        return gatheredRawValues.data() + k * numObservations;
    return nullptr;
  }

  /// Raw test values of column j; a borrowed x_test (setTestPredictors)
  /// supersedes the values gathered at view construction.
  const double* rawTestColumn(size_t j) const {
    if (x_test != nullptr) return x_test + j * numTestObservations;
    for (size_t k = 0; k < gatheredRawColumns.size(); ++k)
      if (gatheredRawColumns[k] == j)
        return gatheredRawTestValues.data() + k * numTestObservations;
    return nullptr;
  }

  /// Parent-derived standardization constants for column j, when the store
  /// is a view that gathered them; false tells the caller to compute its own.
  bool suppliedStandardization(size_t j, double* mean, double* sd) const {
    for (size_t k = 0; k < gatheredRawColumns.size(); ++k) {
      if (gatheredRawColumns[k] != j) continue;
      *mean = gatheredMeans[k];
      *sd = gatheredSds[k];
      return true;
    }
    return false;
  }

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

  void finishQuantileGrid(QuantileGrid& grid, size_t j) const {
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
  }

  QuantileGrid quantileGridForColumn(size_t j, const double* values) const {
    QuantileGrid grid;
    grid.sortedUnique.reserve(numObservations);
    // NaN would break the sort's ordering; the grid is over observed values
    for (size_t i = 0; i < numObservations; ++i)
      if (!isNA(values[i])) grid.sortedUnique.push_back(values[i]);
    finishQuantileGrid(grid, j);
    return grid;
  }

  /// The quantile grid of a CSC column's logical values: the observed stored
  /// entries plus a single zero when implicit zeros exist, the same value
  /// set the dense collector sees, so the induced grid is identical.
  QuantileGrid quantileGridForCscColumn(size_t j) const {
    const CscColumnSlice& slice = cscSlices[j];
    QuantileGrid grid;
    grid.sortedUnique.reserve(slice.numNonzero + 1);
    if (slice.numNonzero < numObservations) grid.sortedUnique.push_back(0.0);
    for (size_t k = 0; k < slice.numNonzero; ++k)
      if (!isNA(slice.values[k])) grid.sortedUnique.push_back(slice.values[k]);
    finishQuantileGrid(grid, j);
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

  void fillCutsOverRange(size_t j, double xMin, double xMax) {
    cutPoints[j].resize(numCuts[j]);
    double increment = (xMax - xMin) / static_cast<double>(numCuts[j] + 1);
    for (std::uint32_t k = 0; k < numCuts[j]; ++k)
      cutPoints[j][k] = xMin + static_cast<double>(k + 1) * increment;
  }

  void fillCutsUniformly(size_t j) {
    const double* column = denseSourceColumn(j);
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
    fillCutsOverRange(j, xMin, xMax);
  }

  /// The dense range scan over a CSC column's logical values: implicit
  /// zeros seed the range at 0, observed stored entries fold in.
  void fillCutsUniformlyCsc(size_t j) {
    const CscColumnSlice& slice = cscSlices[j];
    double xMin = 0.0, xMax = 0.0;
    size_t k = 0;
    if (slice.numNonzero == numObservations) {  // no implicit zeros
      while (k < slice.numNonzero && isNA(slice.values[k])) ++k;
      if (k < slice.numNonzero) {
        xMin = xMax = slice.values[k];
        ++k;
      }
    }
    for (; k < slice.numNonzero; ++k) {
      if (slice.values[k] < xMin) xMin = slice.values[k];
      if (slice.values[k] > xMax) xMax = slice.values[k];
    }
    fillCutsOverRange(j, xMin, xMax);
  }

  /// Initial cut construction; sets numCuts[j]. A categorical column takes
  /// its (fixed) category count from its largest value and keeps no cuts.
  /// CSC columns are always ordinal.
  void buildCutsForColumn(size_t j) {
    if (types[j] == ColumnType::categorical) {
      const double* column = denseSourceColumn(j);
      double maxValue = -1.0;
      for (size_t i = 0; i < numObservations; ++i)
        if (!isNA(column[i]) && column[i] > maxValue) maxValue = column[i];
      numCuts[j] = maxValue < 0.0
        ? 0 : static_cast<std::uint32_t>(maxValue) + 1;
      cutPoints[j].clear();
    } else if (useQuantiles) {
      QuantileGrid grid = columnIsCscBacked(j)
        ? quantileGridForCscColumn(j)
        : quantileGridForColumn(j, denseSourceColumn(j));
      numCuts[j] = grid.inducedNumCuts;
      fillCutsFromQuantileGrid(j, grid);
    } else {
      numCuts[j] = maxNumCuts[j];
      if (columnIsCscBacked(j)) fillCutsUniformlyCsc(j);
      else fillCutsUniformly(j);
    }
  }

  /// No two observed values differ, so a fixed-count uniform grid of two or
  /// more cuts would repeat a value rather than strictly ascend.
  bool valuesAreDegenerate(const double* values) const {
    size_t i = 0;
    while (i < numObservations && isNA(values[i])) ++i;
    if (i >= numObservations) return true;
    double first = values[i];
    for (++i; i < numObservations; ++i)
      if (!isNA(values[i]) && values[i] != first) return false;
    return true;
  }

  /// Recompute cuts for a column's current values, keeping numCuts[j] fixed.
  /// Refuses (returns false, keeping the old grid) when the fixed count cannot
  /// yield a strictly ascending grid: quantile mode with fewer induced cuts
  /// than existing (extra induced cuts are silently thinned), or a degenerate
  /// range under two or more uniform cuts
  /// (a re-cut there would repeat a value). A forced update then routes the
  /// new values through the retained grid and collapses what empties.
  /// Categorical columns have nothing to refresh; the caller pre-checked.
  bool refreshCutsForColumn(size_t j) {
    if (types[j] == ColumnType::categorical) return true;
    if (useQuantiles) {
      QuantileGrid grid = quantileGridForColumn(j, denseSourceColumn(j));
      if (grid.inducedNumCuts < numCuts[j]) return false;
      fillCutsFromQuantileGrid(j, grid);
    } else {
      if (numCuts[j] >= 2 && valuesAreDegenerate(denseSourceColumn(j)))
        return false;
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
    if (columnIsCscBacked(j)) {
      quantizeCscColumn(j);
      return;
    }
    const double* column = denseSourceColumn(j);
    std::uint8_t anyMissing = 0;
    for (size_t i = 0; i < numObservations; ++i) {
      xint_t code = codeFor(j, column[i]);
      if (isNA(column[i])) anyMissing = 1;
      codes[codeOffsets[j] + i] = code;
    }
    hasMissing[j] = anyMissing;
  }

  /// Quantize a CSC column against its current cuts: rank columns rewrite
  /// their packed codes and zero code in place (the pattern is fixed),
  /// densified ones fill with the zero code and scatter the stored entries.
  /// Missing values are stored NaN entries and take the reserved code.
  void quantizeCscColumn(size_t j) {
    const CscColumnSlice& slice = cscSlices[j];
    xint_t zeroCode = codeFor(j, 0.0);
    std::uint8_t anyMissing = 0;
    if (columnIsSparse(j)) {
      SparseColumnData& sparse =
        sparseColumns[static_cast<size_t>(sparseSlot[j])];
      sparse.zeroCode = zeroCode;
      for (size_t k = 0; k < slice.numNonzero; ++k) {
        size_t i = static_cast<size_t>(slice.rows[k]);
        std::uint64_t word = sparse.bits[i >> 6];
        size_t rank = sparse.wordRanks[i >> 6] +
          static_cast<size_t>(std::popcount(
            word & ((std::uint64_t{1} << (i & 63u)) - 1u)));
        sparse.nzCodes[rank] = codeFor(j, slice.values[k]);
        if (isNA(slice.values[k])) anyMissing = 1;
      }
    } else {
      xint_t* column = codes.data() + codeOffsets[j];
      for (size_t i = 0; i < numObservations; ++i) column[i] = zeroCode;
      for (size_t k = 0; k < slice.numNonzero; ++k) {
        column[slice.rows[k]] = codeFor(j, slice.values[k]);
        if (isNA(slice.values[k])) anyMissing = 1;
      }
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
    codeOffsets.resize(p);
    for (size_t j = 0; j < p; ++j) codeOffsets[j] = j * n;
    sparseSlot.assign(p, -1);
    sparseColumns.clear();
    cscSlices.clear();
    mixedRawColumns.clear();
    builtFromCsc = false;
    hasSparse = false;
    hasMissing.assign(p, 0);
    gatheredRawColumns.clear();
    gatheredRawValues.clear();
    gatheredRawTestValues.clear();
    gatheredMeans.clear();
    gatheredSds.clear();

    for (size_t j = 0; j < p; ++j) {
      buildCutsForColumn(j);
      quantizeColumn(j);
    }
    refreshCategoricalTiers();
  }

  /// Build from mixed dense and CSC sources: column j reads dense column
  /// columnSources[j] of denseValues (column-major) when the source is
  /// nonnegative - quantized with the dense arithmetic exactly as build(),
  /// categorical allowed, rawColumn served - or CSC column ~columnSources[j]
  /// of the triple otherwise: ordinal only, rank-bitmap storage at or below
  /// sparseDensityThreshold nonzero fraction, densified codes above, the
  /// borrowed slices serving re-quantization either way. The host validates
  /// structure (row indices unique and in range per column) and that
  /// categorical columns are dense-backed. All pointers are borrowed for
  /// the store's lifetime.
  void buildMixed(const double* denseValues, const int* columnPointers,
                  const int* rowIndices, const double* values,
                  const std::int32_t* columnSources, size_t n, size_t p,
                  const std::uint32_t* maxNumCutsPerColumn,
                  std::uint32_t maxNumCutsScalar, bool useQuantiles_,
                  const ColumnType* columnTypes = nullptr) {
    x = nullptr;
    numObservations = n;
    numPredictors = p;
    useQuantiles = useQuantiles_;
    builtFromCsc = true;
    if (columnTypes != nullptr) {
      types.assign(columnTypes, columnTypes + p);
    } else {
      types.assign(p, ColumnType::ordinal);
    }
    cutPoints.resize(p);
    numCuts.resize(p);
    if (maxNumCutsPerColumn != nullptr) {
      maxNumCuts.assign(maxNumCutsPerColumn, maxNumCutsPerColumn + p);
    } else {
      maxNumCuts.assign(p, maxNumCutsScalar);
    }
    for (size_t j = 0; j < p; ++j)
      if (maxNumCuts[j] > maxNumCutsRepresentable)
        maxNumCuts[j] = maxNumCutsRepresentable;
    hasMissing.assign(p, 0);
    gatheredRawColumns.clear();
    gatheredRawValues.clear();
    gatheredRawTestValues.clear();
    gatheredMeans.clear();
    gatheredSds.clear();

    cscSlices.assign(p, CscColumnSlice());
    mixedRawColumns.assign(p, nullptr);
    sparseSlot.assign(p, -1);
    sparseColumns.clear();
    codeOffsets.assign(p, 0);
    size_t numDenseCodes = 0;
    for (size_t j = 0; j < p; ++j) {
      if (columnSources[j] >= 0) {
        mixedRawColumns[j] =
          denseValues + static_cast<size_t>(columnSources[j]) * n;
        codeOffsets[j] = numDenseCodes;
        numDenseCodes += n;
        continue;
      }
      std::int32_t sourceIndex = ~columnSources[j];
      size_t source = static_cast<size_t>(sourceIndex);
      size_t begin = static_cast<size_t>(columnPointers[source]);
      size_t end = static_cast<size_t>(columnPointers[source + 1]);
      cscSlices[j] = { values + begin, rowIndices + begin, end - begin };
      bool sparse = static_cast<double>(end - begin) <=
        sparseDensityThreshold * static_cast<double>(n);
      if (sparse) {
        sparseSlot[j] = static_cast<std::int32_t>(sparseColumns.size());
        sparseColumns.emplace_back();
      } else {
        codeOffsets[j] = numDenseCodes;
        numDenseCodes += n;
      }
    }
    codes.resize(numDenseCodes);
    hasSparse = !sparseColumns.empty();

    size_t numWords = (n + 63) / 64;
    for (size_t j = 0; j < p; ++j) {
      if (columnIsSparse(j)) {
        SparseColumnData& sparse =
          sparseColumns[static_cast<size_t>(sparseSlot[j])];
        const CscColumnSlice& slice = cscSlices[j];
        sparse.bits.assign(numWords, 0);
        sparse.wordRanks.assign(numWords, 0);
        sparse.nzCodes.assign(slice.numNonzero, 0);
        for (size_t k = 0; k < slice.numNonzero; ++k) {
          size_t i = static_cast<size_t>(slice.rows[k]);
          sparse.bits[i >> 6] |= std::uint64_t{1} << (i & 63u);
        }
        std::uint32_t runningRank = 0;
        for (size_t w = 0; w < numWords; ++w) {
          sparse.wordRanks[w] = runningRank;
          runningRank +=
            static_cast<std::uint32_t>(std::popcount(sparse.bits[w]));
        }
      }
      buildCutsForColumn(j);
      quantizeColumn(j);
    }
    refreshCategoricalTiers();

    numTestObservations = 0;
    x_test = nullptr;
    testCodes.clear();
    testOffset = nullptr;
  }

  /// Build from a CSC (dgCMatrix-layout) predictor matrix: buildMixed with
  /// every column CSC-backed, so all columns ordinal and no raw dense
  /// matrix (x stays null).
  void buildFromCsc(const int* columnPointers, const int* rowIndices,
                    const double* values, size_t n, size_t p,
                    const std::uint32_t* maxNumCutsPerColumn,
                    std::uint32_t maxNumCutsScalar, bool useQuantiles_) {
    std::vector<std::int32_t> allCsc(p);
    for (size_t j = 0; j < p; ++j) allCsc[j] = ~static_cast<std::int32_t>(j);
    buildMixed(nullptr, columnPointers, rowIndices, values, allCsc.data(),
               n, p, maxNumCutsPerColumn, maxNumCutsScalar, useQuantiles_);
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
  /// surface); callers enforce that. rawColumnsToGather names columns whose
  /// raw values a leaf model reads (linear leaves): their subset values and
  /// parent-derived standardization constants are copied so rawColumn and
  /// suppliedStandardization serve them. The view is self-contained: nothing
  /// references the parent after this returns.
  void buildFromParent(const ColumnStore& parent, const size_t* rows,
                       size_t numRows, const size_t* testRows,
                       size_t numTestRows,
                       const size_t* rawColumnsToGather = nullptr,
                       size_t numRawColumnsToGather = 0) {
    x = nullptr;
    numObservations = numRows;
    numPredictors = parent.numPredictors;
    useQuantiles = parent.useQuantiles;
    types = parent.types;
    cutPoints = parent.cutPoints;
    numCuts = parent.numCuts;
    maxNumCuts = parent.maxNumCuts;
    // views densify: gathered codes are fully dense whatever the parent's
    // per-column storage (docs/design/sparse-columns.md)
    codes.resize(numRows * numPredictors);
    codeOffsets.resize(numPredictors);
    for (size_t j = 0; j < numPredictors; ++j) codeOffsets[j] = j * numRows;
    sparseSlot.assign(numPredictors, -1);
    sparseColumns.clear();
    cscSlices.clear();
    mixedRawColumns.clear();
    builtFromCsc = false;
    hasSparse = false;
    hasMissing.assign(numPredictors, 0);
    refreshCategoricalTiers();

    gatheredRawColumns.clear();
    gatheredRawValues.clear();
    gatheredRawTestValues.clear();
    gatheredMeans.clear();
    gatheredSds.clear();
    // gather what the parent can serve raw (dense-backed columns, or its
    // own gathered copies); columns it cannot are left ungathered, so the
    // view's rawColumn returns null and the facade refuses the designation
    for (size_t k = 0; k < numRawColumnsToGather; ++k) {
      const double* parentColumn = parent.rawColumn(rawColumnsToGather[k]);
      if (parentColumn == nullptr) continue;
      size_t slot = gatheredRawColumns.size();
      gatheredRawColumns.push_back(rawColumnsToGather[k]);
      gatheredRawValues.resize(gatheredRawValues.size() + numRows);
      gatheredRawTestValues.resize(gatheredRawTestValues.size() + numTestRows);
      gatheredMeans.resize(slot + 1);
      gatheredSds.resize(slot + 1);
      double* values = gatheredRawValues.data() + slot * numRows;
      for (size_t i = 0; i < numRows; ++i)
        values[i] = parentColumn[rows[i]];
      double* testValues = gatheredRawTestValues.data() + slot * numTestRows;
      for (size_t i = 0; i < numTestRows; ++i)
        testValues[i] = parentColumn[testRows[i]];
      double mean, sd;
      if (!parent.suppliedStandardization(rawColumnsToGather[k], &mean, &sd))
        standardizationMomentsForColumn(parentColumn, parent.numObservations,
                                        &mean, &sd);
      gatheredMeans[slot] = mean;
      gatheredSds[slot] = sd;
    }
    for (size_t j = 0; j < numPredictors; ++j) {
      xint_t missingCode = types[j] == ColumnType::categorical
        ? missingCategoryCode(numCuts[j]) : naCode;
      xint_t* column = codes.data() + j * numRows;
      for (size_t i = 0; i < numRows; ++i) {
        column[i] = parent.codeAt(j, rows[i]);
        if (column[i] == missingCode) hasMissing[j] = 1;
      }
    }
    x_test = nullptr;
    numTestObservations = numTestRows;
    testCodes.resize(numTestRows * numPredictors);
    for (size_t i = 0; i < numTestRows; ++i)
      for (size_t j = 0; j < numPredictors; ++j)
        testCodes[i * numPredictors + j] = parent.codeAt(j, testRows[i]);
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
    codes[codeOffsets[j] + i] = codeFor(j, value);
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
    for (size_t j = 0; j < numPredictors; ++j) codeOffsets[j] = j * n;
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

  /// Dense-stored columns only; rank columns have no contiguous codes.
  const xint_t* column(size_t variable) const {
    return codes.data() + codeOffsets[variable];
  }
  const xint_t* testRow(size_t testObservation) const {
    return testCodes.data() + testObservation * numPredictors;
  }

  bool columnIsSparse(size_t variable) const {
    return sparseSlot[variable] >= 0;
  }
  const SparseColumnData& sparseColumn(size_t variable) const {
    return sparseColumns[static_cast<size_t>(sparseSlot[variable])];
  }
  /// Storage-aware single-code access (tree descent, restore validation).
  xint_t codeAt(size_t variable, size_t i) const {
    return columnIsSparse(variable) ? sparseColumn(variable).at(i)
                                    : codes[codeOffsets[variable] + i];
  }
};

}  // namespace bartcore

#endif  // BARTCORE_DATA_HPP
