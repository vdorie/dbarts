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

/// Per-observation GATHER INDEX type for the hot index buffers (Tree::indices,
/// Forest/VarianceForest::indexBuffer): a subscript into the n-length columns,
/// narrowed from size_t to halve the single biggest hot array
/// (docs/design/reduced-precision-storage.md sec 3a). Must match the C kernel
/// index type misc_index_t (pinned by a static_assert in tree.hpp); DISTINCT
/// from length/count storage, which stays size_t. Capped by a guard that
/// refuses numObservations > UINT32_MAX at ingestion.
using index_t = std::uint32_t;

/// Reserved code for a missing ordinal value; cut counts cap below it so
/// real codes (which reach numCuts) never collide.
constexpr xint_t naCode = 0xFFFFu;
constexpr std::uint32_t maxNumCutsRepresentable = 0xFFFDu;

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

/// Where one column's codes live and what re-quantization reads.
enum class ColumnSourceKind : std::uint8_t {
  denseOwned,     // dense codes in codes[]; re-quantize from the side's owned
                  // dense source (train: the call-time x; test: ownedTestValues)
  denseBorrowed,  // dense codes in codes[]; re-quantize from denseRaw
  cscRank,        // rank-bitmap in the side's sparseColumns[rankSlot]; from slice
  cscDensified    // dense codes in codes[]; re-quantize from slice
};

/// A per-column source descriptor: one instance per column of a row set,
/// carried in a vector sized to numPredictors on every built store. The
/// discriminated fields are read only for the kinds that own them (rankSlot for
/// cscRank; denseRaw for denseBorrowed; slice for the two CSC kinds;
/// declaredCategoryCount for any categorical column whose host declared a level
/// table, train-side only, and refCode for CSC-backed categorical columns on
/// both sides). denseOwned has a side-specific raw source the accessor supplies.
struct ColumnSource {
  ColumnSourceKind kind = ColumnSourceKind::denseOwned;
  std::int32_t rankSlot = -1;          // cscRank: slot into the side's sparseColumns
  const double* denseRaw = nullptr;    // denseBorrowed: borrowed/owned dense column
  CscColumnSlice slice;                // cscRank/cscDensified: retained nonzeros
  std::uint32_t declaredCategoryCount = 0;  // categorical: host's level count K (train only)
  xint_t refCode = 0;                  // CSC categorical: reference-level code

  bool isCscBacked() const {
    return kind == ColumnSourceKind::cscRank ||
           kind == ColumnSourceKind::cscDensified;
  }
  bool isRank() const { return kind == ColumnSourceKind::cscRank; }
};

/// One row set's code storage over the store's shared cut grid: the packed
/// dense codes and their per-column starts, the rank storage of the sparse
/// columns, and the per-column source descriptors. Instantiated train and test
/// on ColumnStore; the grid-dependent operations that fill a block are
/// ColumnStore methods parameterized by the block and its row count.
struct CodeBlock {
  // packed dense codes; per-column starts in codeOffsets (j * numRows for
  // dense-matrix builds, packed among densified columns for CSC builds, unused
  // for rank-stored columns)
  std::vector<xint_t> codes;
  std::vector<size_t> codeOffsets;
  std::vector<SparseColumnData> sparseColumns;
  // per column, where its codes live and what re-quantization reads: the
  // storage kind, the rank slot into sparseColumns, the borrowed dense raw of a
  // mixed build, the retained CSC slice, a categorical column's host-declared
  // level count K, and a CSC-backed categorical column's reference-level code.
  // The declared count carries the levels no training row happens to observe,
  // which no scan of the codes can recover; a CSC categorical column stores
  // only its non-reference entries, so its implicit rows' code is the reference
  // level's own level-order index, not the storage's structural zero.
  std::vector<ColumnSource> sources;

  /// Whether column j quantizes from a borrowed CSC slice (every column of a
  /// buildFromCsc store, the sparse-source columns of a mixed build).
  bool columnIsCscBacked(size_t j) const { return sources[j].isCscBacked(); }

  bool columnIsSparse(size_t j) const { return sources[j].isRank(); }
  const SparseColumnData& sparseColumn(size_t j) const {
    return sparseColumns[static_cast<size_t>(sources[j].rankSlot)];
  }

  /// Dense-stored columns only; rank columns have no contiguous codes.
  const xint_t* column(size_t j) const { return codes.data() + codeOffsets[j]; }

  /// Storage-aware single-code access (tree descent, restore validation),
  /// reading only the columns a rule visits rather than materializing a row.
  xint_t codeAt(size_t j, size_t i) const {
    return columnIsSparse(j) ? sparseColumn(j).at(i)
                             : codes[codeOffsets[j] + i];
  }

  /// A dense-backed column's raw source for re-quantization: the borrowed dense
  /// raw when denseBorrowed, else the side's owned dense fallback at column j
  /// (the call-time x for train, ownedTestValues for test).
  const double* denseRawForColumn(size_t j, size_t numRows,
                                  const double* ownedFallback) const {
    return sources[j].kind == ColumnSourceKind::denseBorrowed
      ? sources[j].denseRaw : ownedFallback + j * numRows;
  }
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
  bool useQuantiles = false;
  // Provenance: built from a parent store (buildFromParent). Gates only the
  // parent-derived standardization constants (suppliedStandardization);
  // refusal decisions read the capability predicates hasRequantizeSource and
  // acceptsNewRawPredictors instead.
  bool isView = false;

  std::vector<ColumnType> types;
  // the training rows' codes over this store's cut grid
  CodeBlock train;
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
  // Owned copy of the test predictors, column-major numTestObservations x
  // numPredictors, taken at buildTest. The engine borrows no test matrix: cut
  // changes re-quantize the test codes from this copy, and rawTestColumn serves
  // it to leaf models. Views hold none (they gather test-subset columns below).
  std::vector<double> ownedTestValues;
  // The test rows' codes over this store's cut grid (types, numCuts,
  // cutPoints), the training-side layout mirrored for the test rows. testCodeAt
  // reads whichever storage a column takes, so descent routes a test row
  // without materializing it. Views hold dense codes only (every test.sources
  // entry denseOwned).
  CodeBlock test;
  // Owned raw of a mixed/CSC test build (the test store owns its raw, so no
  // borrowed pointer survives the build): ownedTestCscValues/ownedTestCscRows
  // pack every CSC-backed test column's nonzeros, which each test.sources entry's
  // slice points into. A dense-backed column's entry instead borrows into
  // ownedTestValues via denseRaw, and a CSC-backed categorical column's carries
  // the reference level's code. Both owned CSC buffers are empty on the dense
  // buildTest path and on views.
  std::vector<double> ownedTestCscValues;
  std::vector<int> ownedTestCscRows;
  // borrowed; added to recorded test fits. buildTest leaves it alone (the
  // caller keeps the lengths consistent), clearTest clears it.
  const double* testOffset = nullptr;

  // Owned raw values of designated columns. build gathers a sampler's leaf
  // covariates (or a data handle's declared columns) so rawColumn serves owned
  // memory after the build borrow releases; buildFromParent gathers a view's
  // leaf covariates from its parent, with standardization constants from the
  // parent's full training data - the same calibration inheritance as the
  // copied cut grid, so every fold runs the prior a full-data fit would.
  std::vector<size_t> gatheredRawColumns;
  std::vector<double> gatheredRawValues;      // column-major, numObservations x q
  std::vector<double> gatheredRawTestValues;  // column-major, numTestObservations x q
  std::vector<double> gatheredMeans;
  std::vector<double> gatheredSds;

  // Owned re-quantize sources for CSC-backed columns after mutation. A column
  // keeps its build-time borrowed slice (R's dgCMatrix i/x slots) until first
  // mutated; setPredictor/updatePredictor writes the new nonzeros here and
  // repoints train.sources[j].slice at them, since the borrow no longer
  // reflects the live values (setCutPoints and state restore re-quantize from
  // the slice). Sized numPredictors on a CSC/mixed build, empty on dense
  // builds and views; ownedCsc*[j] stay empty and cscColumnOwned[j] false
  // until column j is first mutated (docs/design/sparse-columns.md).
  std::vector<std::vector<int>> ownedCscRows;
  std::vector<std::vector<double>> ownedCscValues;
  std::vector<std::uint8_t> cscColumnOwned;

  /// The slot column j occupies in the gathered-raw buffers, or -1 when it is
  /// not gathered. Few columns are gathered (leaf covariates), so a scan.
  std::int32_t gatheredSlotForColumn(size_t j) const {
    for (size_t k = 0; k < gatheredRawColumns.size(); ++k)
      if (gatheredRawColumns[k] == j) return static_cast<std::int32_t>(k);
    return -1;
  }

  /// Whether column j quantizes from a borrowed CSC slice (every column of
  /// a buildFromCsc store, the sparse-source columns of a mixed build).
  bool columnIsCscBacked(size_t j) const {
    return train.columnIsCscBacked(j);
  }

  /// Whether any re-quantize source survives the build: retained CSC slices,
  /// borrowed dense slices, or a caller-supplied x. A view gathers codes from
  /// its parent and retains none, so cut installation and state restore are
  /// refused on it.
  bool hasRequantizeSource() const { return !isView; }

  /// Whether the raw-x mutation surface (setData, setPredictor, the
  /// updatePredictor family) can run: a CSC/mixed build fixes the design at
  /// creation and a view retains no raw, so only a dense-built top-level
  /// store accepts new raw predictors.
  bool acceptsNewRawPredictors() const {
    return hasRequantizeSource() && !builtFromCsc;
  }

  /// Column j's raw values for a re-quantize, given the call-time predictor
  /// matrix x (which the caller supplies for the build's duration): the mixed
  /// build's retained dense slice, x's column for a dense build, null for
  /// CSC-backed columns (which re-quantize from their retained slices instead).
  const double* rawColumnForRequantize(size_t j, const double* x) const {
    if (train.sources[j].isCscBacked()) return nullptr;
    if (train.sources[j].kind == ColumnSourceKind::denseBorrowed)
      return train.sources[j].denseRaw;
    return x != nullptr ? x + j * numObservations : nullptr;
  }

  /// Owned raw training values of column j: a gathered copy (leaf covariates,
  /// or a data handle's declared columns), the mixed build's retained dense
  /// slice, null when neither serves it.
  const double* rawColumn(size_t j) const {
    std::int32_t slot = gatheredSlotForColumn(j);
    if (slot >= 0)
      return gatheredRawValues.data() +
             static_cast<size_t>(slot) * numObservations;
    return train.sources[j].kind == ColumnSourceKind::denseBorrowed
      ? train.sources[j].denseRaw : nullptr;
  }

  /// Raw test values of column j: the mixed build's owned dense slice (null for
  /// its CSC-backed columns, which serve no leaf raw), the dense buildTest copy,
  /// or the values gathered at view construction.
  const double* rawTestColumn(size_t j) const {
    if (!test.sources.empty()) {
      const ColumnSource& source = test.sources[j];
      if (source.kind == ColumnSourceKind::denseBorrowed) return source.denseRaw;
      if (source.isCscBacked()) return nullptr;
    }
    if (!ownedTestValues.empty())
      return ownedTestValues.data() + j * numTestObservations;
    std::int32_t slot = gatheredSlotForColumn(j);
    if (slot >= 0)
      return gatheredRawTestValues.data() +
             static_cast<size_t>(slot) * numTestObservations;
    return nullptr;
  }

  /// Parent-derived standardization constants for column j, when the store
  /// is a view that gathered them; false tells the caller to compute its own
  /// (the top-level store standardizes from its own gathered values).
  bool suppliedStandardization(size_t j, double* mean, double* sd) const {
    if (!isView) return false;
    std::int32_t slot = gatheredSlotForColumn(j);
    if (slot < 0) return false;
    *mean = gatheredMeans[static_cast<size_t>(slot)];
    *sd = gatheredSds[static_cast<size_t>(slot)];
    return true;
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
    const double* first = cuts.data();
    return static_cast<xint_t>(
      std::lower_bound(first, first + numCuts[variable], value) - first);
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
    const CscColumnSlice& slice = train.sources[j].slice;
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

  void fillCutsUniformly(size_t j, const double* column) {
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
    const CscColumnSlice& slice = train.sources[j].slice;
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

  /// The category count a dense column's own codes imply: the largest observed
  /// code plus one, and 0 when every value is missing.
  std::uint32_t inferredCategoryCount(const double* column) const {
    double maxValue = -1.0;
    for (size_t i = 0; i < numObservations; ++i)
      if (!isNA(column[i]) && column[i] > maxValue) maxValue = column[i];
    return maxValue < 0.0 ? 0u : static_cast<std::uint32_t>(maxValue) + 1u;
  }

  /// The same over a CSC-backed column's logical values: its retained
  /// nonzeros, plus the reference code the implicit rows read when it has any.
  std::uint32_t inferredCategoryCountCsc(size_t j) const {
    const ColumnSource& source = train.sources[j];
    double maxValue = source.slice.numNonzero < numObservations
      ? static_cast<double>(source.refCode) : -1.0;
    for (size_t k = 0; k < source.slice.numNonzero; ++k)
      if (!isNA(source.slice.values[k]) && source.slice.values[k] > maxValue)
        maxValue = source.slice.values[k];
    return maxValue < 0.0 ? 0u : static_cast<std::uint32_t>(maxValue) + 1u;
  }

  /// Initial cut construction; sets numCuts[j]. column supplies the dense raw
  /// values (null for CSC-backed columns, which read their retained slice). A
  /// categorical column keeps no cuts and takes its (fixed) category count from
  /// the host's declared level table when it supplied one, but never below what
  /// its own codes reach - a declared count short of an observed code would
  /// strand that code past its own grid.
  void buildCutsForColumn(size_t j, const double* column) {
    if (types[j] == ColumnType::categorical) {
      std::uint32_t inferred = columnIsCscBacked(j)
        ? inferredCategoryCountCsc(j) : inferredCategoryCount(column);
      std::uint32_t declared = train.sources[j].declaredCategoryCount;
      numCuts[j] = declared > inferred ? declared : inferred;
      cutPoints[j].clear();
    } else if (useQuantiles) {
      QuantileGrid grid = columnIsCscBacked(j)
        ? quantileGridForCscColumn(j)
        : quantileGridForColumn(j, column);
      numCuts[j] = grid.inducedNumCuts;
      fillCutsFromQuantileGrid(j, grid);
    } else {
      numCuts[j] = maxNumCuts[j];
      if (columnIsCscBacked(j)) fillCutsUniformlyCsc(j);
      else fillCutsUniformly(j, column);
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
  bool refreshCutsForColumn(size_t j, const double* column) {
    if (types[j] == ColumnType::categorical) return true;
    if (useQuantiles) {
      QuantileGrid grid = quantileGridForColumn(j, column);
      if (grid.inducedNumCuts < numCuts[j]) return false;
      fillCutsFromQuantileGrid(j, grid);
    } else {
      if (numCuts[j] >= 2 && valuesAreDegenerate(column)) return false;
      fillCutsUniformly(j, column);
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

  /// Install externally chosen cut points (ascending) for a column and
  /// re-quantize its codes against them; x is the call-time predictor matrix
  /// the raw is read from (ignored for CSC-backed columns, which use their
  /// retained slice). The cut count may shrink or grow, and existing splits
  /// beyond the new range are the caller's problem (the sampler collapses them).
  void setCutPointsForColumn(size_t j, const double* cuts,
                             std::uint32_t numCutPoints, const double* x) {
    cutPoints[j].assign(cuts, cuts + numCutPoints);
    numCuts[j] = numCutPoints;
    if (maxNumCuts[j] < numCutPoints) maxNumCuts[j] = numCutPoints;
    quantizeColumn(j, rawColumnForRequantize(j, x));
    if (numTestObservations > 0) quantizeTestColumn(j);
  }

  /// Quantize column j's dense raw values into block against the current cuts,
  /// recording any missingness through hasMissingOut[j] when non-null. observe
  /// is invoked as observe(i, column, code) for each cell just before column[i]
  /// is overwritten, giving a caller the old code without a second pass; it is a
  /// compile-time functor, inlined away for the no-op case.
  template <typename Observer>
  void quantizeDenseObserved(CodeBlock& block, size_t numRows, size_t j,
                             const double* raw, std::uint8_t* hasMissingOut,
                             Observer observe) {
    xint_t* column = block.codes.data() + block.codeOffsets[j];
    std::uint8_t anyMissing = 0;
    for (size_t i = 0; i < numRows; ++i) {
      xint_t code = codeFor(j, raw[i]);
      if (isNA(raw[i])) anyMissing = 1;
      observe(i, column, code);
      column[i] = code;
    }
    if (hasMissingOut != nullptr) hasMissingOut[j] = anyMissing;
  }

  /// train tracks missingness; test tracks none and passes null.
  void quantizeDenseInto(CodeBlock& block, size_t numRows, size_t j,
                         const double* raw, std::uint8_t* hasMissingOut) {
    quantizeDenseObserved(block, numRows, j, raw, hasMissingOut,
                          [](size_t, const xint_t*, xint_t) {});
  }

  /// Quantize a CSC-backed column j into block against its current cuts: rank
  /// columns rewrite their packed codes and zero code in place (the pattern is
  /// fixed at build), densified ones fill with the zero code and scatter the
  /// stored entries. A categorical column's implicit rows carry the reference
  /// level's level-order code; an ordinal column's carry the quantized zero.
  /// Missing values are stored NaN entries and take the reserved code, recorded
  /// through hasMissingOut[j] when non-null.
  void quantizeCscColumnInto(CodeBlock& block, size_t numRows, size_t j,
                             std::uint8_t* hasMissingOut) {
    const CscColumnSlice& slice = block.sources[j].slice;
    xint_t zeroCode = types[j] == ColumnType::categorical
      ? block.sources[j].refCode : codeFor(j, 0.0);
    std::uint8_t anyMissing = 0;
    if (block.columnIsSparse(j)) {
      SparseColumnData& sparse =
        block.sparseColumns[static_cast<size_t>(block.sources[j].rankSlot)];
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
      xint_t* column = block.codes.data() + block.codeOffsets[j];
      for (size_t i = 0; i < numRows; ++i) column[i] = zeroCode;
      for (size_t k = 0; k < slice.numNonzero; ++k) {
        column[static_cast<size_t>(slice.rows[k])] =
          codeFor(j, slice.values[k]);
        if (isNA(slice.values[k])) anyMissing = 1;
      }
    }
    if (hasMissingOut != nullptr) hasMissingOut[j] = anyMissing;
  }

  /// Build column j's rank-bitmap pattern in block from its retained CSC slice:
  /// the fixed bits/wordRanks and the nzCodes capacity. A non-sparse column has
  /// none. Must precede the quantize that writes nzCodes by rank.
  void buildRankStorageInto(CodeBlock& block, size_t numRows, size_t j) {
    if (!block.columnIsSparse(j)) return;
    SparseColumnData& sparse =
      block.sparseColumns[static_cast<size_t>(block.sources[j].rankSlot)];
    const CscColumnSlice& slice = block.sources[j].slice;
    size_t numWords = (numRows + 63) / 64;
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
      runningRank += static_cast<std::uint32_t>(std::popcount(sparse.bits[w]));
    }
  }

  /// Re-quantize column j's codes from column (the dense raw, or null for a
  /// CSC-backed column, which reads its retained slice), refreshing the
  /// gathered raw copy of a leaf-covariate column in the same pass.
  void quantizeColumn(size_t j, const double* column) {
    if (columnIsCscBacked(j)) {
      quantizeCscColumnInto(train, numObservations, j, hasMissing.data());
      return;
    }
    quantizeDenseInto(train, numObservations, j, column, hasMissing.data());
    std::int32_t slot = gatheredSlotForColumn(j);
    if (slot >= 0)
      std::memcpy(gatheredRawValues.data() +
                    static_cast<size_t>(slot) * numObservations,
                  column, numObservations * sizeof(double));
  }

  /// Whether test column j re-quantizes from a retained CSC test slice (the
  /// CSC-backed columns of a mixed test build), test analog of columnIsCscBacked.
  bool testColumnIsCscBacked(size_t j) const {
    return test.columnIsCscBacked(j);
  }

  /// Re-quantize test column j against the current cuts. A CSC-backed column
  /// reads its retained owned slice; a dense-backed one reads the owned dense
  /// raw (the mixed build's slice or the buildTest per-column copy). The test
  /// side gates no draws, so it tracks no missingness.
  void quantizeTestColumn(size_t j) {
    if (testColumnIsCscBacked(j)) {
      quantizeCscColumnInto(test, numTestObservations, j, nullptr);
      return;
    }
    const double* raw =
      test.denseRawForColumn(j, numTestObservations, ownedTestValues.data());
    quantizeDenseInto(test, numTestObservations, j, raw, nullptr);
  }

  /// Designate the columns build/setData own an owned raw copy of, so
  /// rawColumn serves them after the build borrow releases. gatherColumns are
  /// a sampler's leaf covariates, or a data handle's declared columns; the
  /// buffers fill as each column quantizes. Cleared when nothing is gathered.
  void setupGatheredColumns(const size_t* gatherColumns,
                            size_t numGatherColumns) {
    gatheredRawColumns.assign(gatherColumns, gatherColumns + numGatherColumns);
    gatheredRawValues.assign(numGatherColumns * numObservations, 0.0);
    gatheredRawTestValues.clear();
    gatheredMeans.clear();
    gatheredSds.clear();
  }

  /// Reset the per-column source storage to the dense-empty baseline: every
  /// column denseOwned, no CSC or rank backing, no gathered raw, no recorded
  /// missingness. Every train builder resets through here, then overwrites the
  /// per-column source descriptors its storage kind owns (buildMixed the CSC
  /// slices, counts, and reference codes; build the gathered raw via
  /// setupGatheredColumns). numPredictors must already hold the new column
  /// count; codes/codeOffsets and the cut grid are sized by each builder.
  void resetTrainStorage() {
    train.sources.assign(numPredictors, ColumnSource{});
    train.sparseColumns.clear();
    builtFromCsc = false;
    hasSparse = false;
    hasMissing.assign(numPredictors, 0);
    gatheredRawColumns.clear();
    gatheredRawValues.clear();
    gatheredRawTestValues.clear();
    gatheredMeans.clear();
    gatheredSds.clear();
    // no CSC-backed columns until a CSC/mixed build re-sizes these; a dense
    // build leaves them empty and never reaches the mutation-owned path
    ownedCscRows.clear();
    ownedCscValues.clear();
    cscColumnOwned.clear();
  }

  /// columnTypes may be null for all-ordinal. Categorical columns must hold
  /// integral values 0..K-1 with K <= maxCategories; the caller validates.
  /// gatherColumns names the columns whose raw values a leaf model reads (or a
  /// data handle's declared columns): build owns a copy so the raw source
  /// borrow x_ need not outlive the call. categoryCounts_, when supplied, is
  /// the host's declared level count per column (0 = infer from the codes).
  void build(const double* x_, size_t n, size_t p,
             const std::uint32_t* maxNumCuts_, bool useQuantiles_,
             const ColumnType* columnTypes = nullptr,
             const size_t* gatherColumns = nullptr,
             size_t numGatherColumns = 0,
             const std::uint32_t* categoryCounts_ = nullptr) {
    isView = false;
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
    train.codes.resize(n * p);
    train.codeOffsets.resize(p);
    for (size_t j = 0; j < p; ++j) train.codeOffsets[j] = j * n;
    resetTrainStorage();
    setupGatheredColumns(gatherColumns, numGatherColumns);
    // the declared level tables, consumed here by buildCutsForColumn; a null
    // array leaves every count inferred
    if (categoryCounts_ != nullptr)
      for (size_t j = 0; j < p; ++j)
        train.sources[j].declaredCategoryCount = categoryCounts_[j];

    for (size_t j = 0; j < p; ++j) {
      const double* column = x_ + j * n;
      buildCutsForColumn(j, column);
      quantizeColumn(j, column);
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
                  const ColumnType* columnTypes = nullptr,
                  const std::uint32_t* categoryCounts_ = nullptr,
                  const xint_t* cscReferenceCodes_ = nullptr) {
    isView = false;
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
    if (maxNumCutsPerColumn != nullptr) {
      maxNumCuts.assign(maxNumCutsPerColumn, maxNumCutsPerColumn + p);
    } else {
      maxNumCuts.assign(p, maxNumCutsScalar);
    }
    for (size_t j = 0; j < p; ++j)
      if (maxNumCuts[j] > maxNumCutsRepresentable)
        maxNumCuts[j] = maxNumCutsRepresentable;
    resetTrainStorage();
    // this builder owns the per-column CSC sources: borrowed slices,
    // densified/rank backing, and the reference codes of CSC-backed
    // categoricals (declared level counts ride both storage kinds, below)
    builtFromCsc = true;
    // per-column owned-mutation storage, empty until a column is first mutated
    ownedCscRows.assign(p, {});
    ownedCscValues.assign(p, {});
    cscColumnOwned.assign(p, 0);
    train.codeOffsets.assign(p, 0);
    size_t numDenseCodes = 0;
    for (size_t j = 0; j < p; ++j) {
      ColumnSource& desc = train.sources[j];
      // the declared level table rides either storage kind: a dense-backed
      // factor inside the container declares one exactly as a CSC-backed one
      if (categoryCounts_ != nullptr)
        desc.declaredCategoryCount = categoryCounts_[j];
      if (columnSources[j] >= 0) {
        desc.kind = ColumnSourceKind::denseBorrowed;
        desc.denseRaw = denseValues + static_cast<size_t>(columnSources[j]) * n;
        train.codeOffsets[j] = numDenseCodes;
        numDenseCodes += n;
        continue;
      }
      std::int32_t sourceIndex = ~columnSources[j];
      size_t source = static_cast<size_t>(sourceIndex);
      size_t begin = static_cast<size_t>(columnPointers[source]);
      size_t end = static_cast<size_t>(columnPointers[source + 1]);
      desc.slice = { values + begin, rowIndices + begin, end - begin };
      if (cscReferenceCodes_ != nullptr) desc.refCode = cscReferenceCodes_[j];
      bool sparse = static_cast<double>(end - begin) <=
        sparseDensityThreshold * static_cast<double>(n);
      if (sparse) {
        desc.kind = ColumnSourceKind::cscRank;
        desc.rankSlot = static_cast<std::int32_t>(train.sparseColumns.size());
        train.sparseColumns.emplace_back();
      } else {
        desc.kind = ColumnSourceKind::cscDensified;
        train.codeOffsets[j] = numDenseCodes;
        numDenseCodes += n;
      }
    }
    train.codes.resize(numDenseCodes);
    hasSparse = !train.sparseColumns.empty();

    for (size_t j = 0; j < p; ++j) {
      buildRankStorageInto(train, n, j);
      // train.sources[j].denseRaw is the dense slice for a dense-backed column
      // and null for a CSC-backed one (which quantizes from its retained slice)
      buildCutsForColumn(j, train.sources[j].denseRaw);
      quantizeColumn(j, train.sources[j].denseRaw);
    }
    refreshCategoricalTiers();

    resetTestStorage();
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
             const ColumnType* columnTypes = nullptr,
             const size_t* gatherColumns = nullptr,
             size_t numGatherColumns = 0,
             const std::uint32_t* categoryCounts_ = nullptr) {
    std::vector<std::uint32_t> maxPerColumn(p, maxNumCuts_);
    build(x_, n, p, maxPerColumn.data(), useQuantiles_, columnTypes,
          gatherColumns, numGatherColumns, categoryCounts_);
  }

  /// Clear the owned-CSC test payload (the packed nonzeros); the per-column CSC
  /// state now lives in test.sources, which callers reset directly.
  void clearTestCscSources() {
    ownedTestCscValues.clear();
    ownedTestCscRows.clear();
  }

  /// Reset the whole test store to empty: no test rows or codes, no CSC or rank
  /// test backing, no test offset. The reset-to-empty test sites route through
  /// here; active test builders size the fields they populate and keep any
  /// caller-owned test offset.
  void resetTestStorage() {
    numTestObservations = 0;
    ownedTestValues.clear();
    test.codes.clear();
    test.codeOffsets.clear();
    test.sources.clear();
    test.sparseColumns.clear();
    clearTestCscSources();
    testOffset = nullptr;
  }

  /// Own a copy of the test predictors and quantize them against the current
  /// cuts; the raw borrow x_test_ need not outlive the call. Re-quantization
  /// on a later cut change reads the owned copy.
  void buildTest(const double* x_test_, size_t numTest) {
    numTestObservations = numTest;
    ownedTestValues.assign(x_test_, x_test_ + numTest * numPredictors);
    test.codes.resize(numTest * numPredictors);
    test.codeOffsets.resize(numPredictors);
    for (size_t j = 0; j < numPredictors; ++j) test.codeOffsets[j] = j * numTest;
    test.sources.assign(numPredictors, ColumnSource{});
    test.sparseColumns.clear();
    clearTestCscSources();
    for (size_t j = 0; j < numPredictors; ++j) quantizeTestColumn(j);
  }

  /// Build the test store from mixed dense and CSC sources against the training
  /// cut grid (already built, shared by identity): column j reads dense column
  /// columnSources[j] of denseValues when nonnegative - quantized like the dense
  /// buildTest, raw served through rawTestColumn - or CSC column
  /// ~columnSources[j] of the triple otherwise, rank-bitmap storage at or below
  /// sparseDensityThreshold nonzero fraction and densified codes above, the
  /// training tier rule per column. cscReferenceCodes gives each CSC-backed
  /// categorical test column its reference level code (the code the implicit
  /// rows take), null when none is categorical. The test store owns its raw:
  /// the dense block and the CSC nonzero values+rows are copied, so no borrowed
  /// pointer survives the call. numCuts and cutPoints are not rebuilt.
  void buildTestMixed(const double* denseValues, const int* columnPointers,
                      const int* rowIndices, const double* values,
                      const std::int32_t* columnSources, size_t numTest,
                      const xint_t* cscReferenceCodes = nullptr) {
    size_t p = numPredictors;
    numTestObservations = numTest;

    // own the dense block, sized by the largest dense source it is indexed by
    std::int32_t maxDenseSource = -1;
    for (size_t j = 0; j < p; ++j)
      if (columnSources[j] > maxDenseSource) maxDenseSource = columnSources[j];
    size_t numDenseColumns =
      maxDenseSource >= 0 ? static_cast<size_t>(maxDenseSource) + 1 : 0;
    ownedTestValues.assign(denseValues, denseValues + numDenseColumns * numTest);

    test.sources.assign(p, ColumnSource{});
    test.sparseColumns.clear();
    test.codeOffsets.assign(p, 0);

    // own the CSC nonzeros: pack them so each column's slice points into
    // storage that never reallocates for the store's lifetime
    size_t totalNonzero = 0;
    for (size_t j = 0; j < p; ++j)
      if (columnSources[j] < 0) {
        size_t source = static_cast<size_t>(~columnSources[j]);
        totalNonzero += static_cast<size_t>(columnPointers[source + 1] -
                                            columnPointers[source]);
      }
    ownedTestCscValues.resize(totalNonzero);
    ownedTestCscRows.resize(totalNonzero);

    size_t numDenseCodes = 0, cursor = 0;
    for (size_t j = 0; j < p; ++j) {
      ColumnSource& desc = test.sources[j];
      if (columnSources[j] >= 0) {
        desc.kind = ColumnSourceKind::denseBorrowed;
        desc.denseRaw = ownedTestValues.data() +
          static_cast<size_t>(columnSources[j]) * numTest;
        test.codeOffsets[j] = numDenseCodes;
        numDenseCodes += numTest;
        continue;
      }
      size_t source = static_cast<size_t>(~columnSources[j]);
      size_t begin = static_cast<size_t>(columnPointers[source]);
      size_t end = static_cast<size_t>(columnPointers[source + 1]);
      size_t numNonzero = end - begin;
      std::memcpy(ownedTestCscValues.data() + cursor, values + begin,
                  numNonzero * sizeof(double));
      std::memcpy(ownedTestCscRows.data() + cursor, rowIndices + begin,
                  numNonzero * sizeof(int));
      desc.slice = { ownedTestCscValues.data() + cursor,
                     ownedTestCscRows.data() + cursor, numNonzero };
      if (cscReferenceCodes != nullptr) desc.refCode = cscReferenceCodes[j];
      cursor += numNonzero;
      bool sparse = static_cast<double>(numNonzero) <=
        sparseDensityThreshold * static_cast<double>(numTest);
      if (sparse) {
        desc.kind = ColumnSourceKind::cscRank;
        desc.rankSlot = static_cast<std::int32_t>(test.sparseColumns.size());
        test.sparseColumns.emplace_back();
      } else {
        desc.kind = ColumnSourceKind::cscDensified;
        test.codeOffsets[j] = numDenseCodes;
        numDenseCodes += numTest;
      }
    }
    test.codes.resize(numDenseCodes);

    for (size_t j = 0; j < p; ++j) {
      buildRankStorageInto(test, numTest, j);
      quantizeTestColumn(j);
    }
  }

  /// A row- and column-subset view of a built parent store: copies the
  /// parent's cut structure and gathers the subset's codes, so the view bins
  /// identically to the parent by construction; testRows also index the
  /// parent's observations. columns, when given, selects the parent columns
  /// the view spans (view-local column j reads parent column columns[j]); null
  /// spans every parent column, the full-span view unchanged. Views hold no
  /// re-quantizable raw (hasRequantizeSource is false), which rules out every
  /// mutation and re-quantize path here; callers enforce that; isView records
  /// the provenance. rawColumnsToGather names the
  /// view's own columns whose raw values a leaf model reads (linear leaves):
  /// their subset values and parent-derived standardization constants are
  /// copied so rawColumn and suppliedStandardization serve them. The view is
  /// self-contained: nothing references the parent after this returns.
  void buildFromParent(const ColumnStore& parent, const size_t* rows,
                       size_t numRows, const size_t* testRows,
                       size_t numTestRows,
                       const size_t* rawColumnsToGather = nullptr,
                       size_t numRawColumnsToGather = 0,
                       const size_t* columns = nullptr, size_t numColumns = 0) {
    isView = true;
    numObservations = numRows;
    numPredictors = columns != nullptr ? numColumns : parent.numPredictors;
    useQuantiles = parent.useQuantiles;
    // view-local column j reads parent column parentColumns[j]; the absent
    // list is the identity, so the full-span view stays byte-for-byte
    std::vector<size_t> parentColumns(numPredictors);
    for (size_t j = 0; j < numPredictors; ++j)
      parentColumns[j] = columns != nullptr ? columns[j] : j;
    if (columns != nullptr) {
      types.resize(numPredictors);
      cutPoints.resize(numPredictors);
      numCuts.resize(numPredictors);
      maxNumCuts.resize(numPredictors);
      for (size_t j = 0; j < numPredictors; ++j) {
        types[j] = parent.types[parentColumns[j]];
        cutPoints[j] = parent.cutPoints[parentColumns[j]];
        numCuts[j] = parent.numCuts[parentColumns[j]];
        maxNumCuts[j] = parent.maxNumCuts[parentColumns[j]];
      }
    } else {
      types = parent.types;
      cutPoints = parent.cutPoints;
      numCuts = parent.numCuts;
      maxNumCuts = parent.maxNumCuts;
    }
    // views densify: gathered codes are fully dense whatever the parent's
    // per-column storage (docs/design/sparse-columns.md)
    train.codes.resize(numRows * numPredictors);
    train.codeOffsets.resize(numPredictors);
    for (size_t j = 0; j < numPredictors; ++j) train.codeOffsets[j] = j * numRows;
    resetTrainStorage();
    refreshCategoricalTiers();
    ownedTestValues.clear();
    // gather what the parent can serve raw (dense-backed columns, or its
    // own gathered copies); columns it cannot are left ungathered, so the
    // view's rawColumn returns null and the facade refuses the designation.
    // rawColumnsToGather index the view's columns, read from their parent map
    for (size_t k = 0; k < numRawColumnsToGather; ++k) {
      size_t parentColumnIndex = parentColumns[rawColumnsToGather[k]];
      const double* parentColumn = parent.rawColumn(parentColumnIndex);
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
      if (!parent.suppliedStandardization(parentColumnIndex, &mean, &sd))
        standardizationMomentsForColumn(parentColumn, parent.numObservations,
                                        &mean, &sd);
      gatheredMeans[slot] = mean;
      gatheredSds[slot] = sd;
    }
    for (size_t j = 0; j < numPredictors; ++j) {
      xint_t missingCode = types[j] == ColumnType::categorical
        ? missingCategoryCode(numCuts[j]) : naCode;
      xint_t* column = train.codes.data() + j * numRows;
      for (size_t i = 0; i < numRows; ++i) {
        column[i] = parent.codeAt(parentColumns[j], rows[i]);
        if (column[i] == missingCode) hasMissing[j] = 1;
      }
    }
    // views densify their test codes too: dense per-column codes gathered
    // from the parent's storage-aware codeAt, no sparse test slots
    numTestObservations = numTestRows;
    test.codes.resize(numTestRows * numPredictors);
    test.codeOffsets.resize(numPredictors);
    for (size_t j = 0; j < numPredictors; ++j)
      test.codeOffsets[j] = j * numTestRows;
    test.sources.assign(numPredictors, ColumnSource{});
    test.sparseColumns.clear();
    clearTestCscSources();
    for (size_t j = 0; j < numPredictors; ++j) {
      xint_t* column = test.codes.data() + j * numTestRows;
      for (size_t i = 0; i < numTestRows; ++i)
        column[i] = parent.codeAt(parentColumns[j], testRows[i]);
    }
    testOffset = nullptr;
  }

  // Mutation. The new raw values arrive as call arguments (the engine keeps no
  // predictor matrix to write through); snapshot/rollback of cutPoints, codes,
  // and the gathered leaf raw is the caller's (the sampler's) responsibility.
  // Cut refreshes assume the caller pre-checked quantile feasibility with
  // cutsWouldRemainValid. These paths run only on dense-built stores (the
  // bridge refuses mutation on CSC/mixed and view stores).

  /// Replace the whole predictor matrix; newX is column-major and read for
  /// the call only, quantized into the owned codes. CSC-backed columns route
  /// through the sparse mutation path (their new dense column rebuilds the
  /// owned slice and rank/densified storage).
  void setPredictors(const double* newX, bool updateCuts) {
    for (size_t j = 0; j < numPredictors; ++j) {
      const double* column = newX + j * numObservations;
      if (columnIsCscBacked(j)) {
        mutateCscColumnFromDense(j, column, updateCuts);
        continue;
      }
      if (updateCuts) refreshCutsForColumn(j, column);
      quantizeColumn(j, column);
    }
  }

  /// Overwrite a subset of columns; newColumns is column-major,
  /// numObservations x numColumns, read for the call only. CSC-backed columns
  /// route through the sparse mutation path.
  void setColumns(const double* newColumns, const size_t* columns,
                  size_t numColumns, bool updateCuts) {
    for (size_t k = 0; k < numColumns; ++k) {
      size_t j = columns[k];
      const double* column = newColumns + k * numObservations;
      if (columnIsCscBacked(j)) {
        mutateCscColumnFromDense(j, column, updateCuts);
        continue;
      }
      if (updateCuts) refreshCutsForColumn(j, column);
      quantizeColumn(j, column);
    }
  }

  /// Repoint CSC-backed column j's slice at its owned nonzero buffers, so
  /// re-quantization no longer reads the build-time borrow (R's dgCMatrix).
  void repointOwnedSlice(size_t j) {
    train.sources[j].slice = { ownedCscValues[j].data(),
                               ownedCscRows[j].data(),
                               ownedCscValues[j].size() };
  }

  /// Mutate CSC-backed column j from a new dense column of numObservations
  /// values (the mutation surface hands a dense column even for sparse
  /// storage). The nonzero pattern is keyed on the column's KIND, since that
  /// fixes what an implicit row reads: an ordinal column's implicit rows read a
  /// structural zero, so the pattern is {i : value != 0}; a categorical
  /// column's read the reference level's own level-order code, so the pattern
  /// is {i : code != refCode}. Either way a stored NaN (missing) stays stored
  /// (NaN compares unequal to both) and the minimal pattern a dense equivalent
  /// would carry is produced - codes stay bitwise identical to a dense build of
  /// the same values. When the pattern is unchanged the nonzero values
  /// re-quantize IN PLACE (rank: nzCodes and zeroCode; densified: the codes
  /// segment); when it changes the rank bitmap and index REBUILD
  /// (O(n / 64 + nnz)). The owned slice repoints at the new nonzeros so later
  /// re-quantizes (setCutPoints, state restore) read them. The storage tier
  /// (rank vs densified) is fixed at build and never flips. updateCuts
  /// refreshes the ordinal cut grid from the dense column exactly as the dense
  /// path does (the CSC cut builders fold the same implicit zeros the dense
  /// column carries explicitly, so the grids match); a categorical column has
  /// no grid to refresh and keeps its creation-pinned category count. The
  /// caller snapshots for rollback first.
  void mutateCscColumnFromDense(size_t j, const double* column,
                                bool updateCuts) {
    if (updateCuts && types[j] != ColumnType::categorical)
      refreshCutsForColumn(j, column);

    const double implicitValue = types[j] == ColumnType::categorical
      ? static_cast<double>(train.sources[j].refCode) : 0.0;
    std::vector<int> newRows;
    std::vector<double> newValues;
    for (size_t i = 0; i < numObservations; ++i)
      if (column[i] != implicitValue) {
        newRows.push_back(static_cast<int>(i));
        newValues.push_back(column[i]);
      }

    const CscColumnSlice& oldSlice = train.sources[j].slice;
    bool patternChanged = newRows.size() != oldSlice.numNonzero;
    for (size_t k = 0; !patternChanged && k < newRows.size(); ++k)
      if (newRows[k] != oldSlice.rows[k]) patternChanged = true;

    ownedCscRows[j] = std::move(newRows);
    ownedCscValues[j] = std::move(newValues);
    cscColumnOwned[j] = 1;
    repointOwnedSlice(j);

    if (train.columnIsSparse(j) && patternChanged)
      buildRankStorageInto(train, numObservations, j);
    quantizeCscColumnInto(train, numObservations, j, hasMissing.data());
  }

  /// A CSC-backed column's rollback record for a transactional mutation: the
  /// pre-change source descriptor (slice pointers, kind, rank slot), the
  /// pre-change rank storage (rank tier) or codes segment (densified tier),
  /// the pre-change owned nonzero buffers, and whether the slice was already
  /// owned. Drives the subset (updatePredictor) transaction; the whole-matrix
  /// transaction snapshots the same state in bulk.
  struct CscColumnRollback {
    ColumnSource source;
    SparseColumnData sparse;
    std::vector<xint_t> denseCodes;
    std::vector<int> ownedRows;
    std::vector<double> ownedValues;
    bool wasOwned = false;
  };

  /// Snapshot CSC-backed column j before a mutation, into rollback.
  void snapshotCscColumn(size_t j, CscColumnRollback& rollback) {
    rollback.source = train.sources[j];
    rollback.wasOwned = cscColumnOwned[j] != 0;
    rollback.ownedRows = ownedCscRows[j];
    rollback.ownedValues = ownedCscValues[j];
    if (train.columnIsSparse(j)) {
      rollback.sparse = train.sparseColumns[
        static_cast<size_t>(train.sources[j].rankSlot)];
    } else {
      const xint_t* col = train.codes.data() + train.codeOffsets[j];
      rollback.denseCodes.assign(col, col + numObservations);
    }
  }

  /// Undo a CSC-backed column's mutation from its rollback record, restoring
  /// the store byte-for-byte (owned slices repoint at the restored buffers,
  /// since a moved-then-restored buffer need not sit at the pre-change
  /// address; a slice that was borrowed keeps the restored R pointer).
  void restoreCscColumn(size_t j, const CscColumnRollback& rollback) {
    train.sources[j] = rollback.source;
    cscColumnOwned[j] = rollback.wasOwned ? 1 : 0;
    ownedCscRows[j] = rollback.ownedRows;
    ownedCscValues[j] = rollback.ownedValues;
    if (train.columnIsSparse(j)) {
      train.sparseColumns[static_cast<size_t>(train.sources[j].rankSlot)] =
        rollback.sparse;
    } else {
      std::memcpy(train.codes.data() + train.codeOffsets[j],
                  rollback.denseCodes.data(),
                  numObservations * sizeof(xint_t));
    }
    if (cscColumnOwned[j]) repointOwnedSlice(j);
  }

  /// A column's rollback record for a journaled re-quantize: either the
  /// changed cells' pre-change codes, or, once too many cells change, the whole
  /// pre-change column. Exactly one is populated.
  struct ColumnCodeRollback {
    struct Cell { std::uint32_t index; xint_t oldCode; };
    std::vector<Cell> journal;
    std::vector<xint_t> fullColumn;
    bool full = false;
  };

  /// Re-quantize column j from newColumn (refreshing its cuts first when
  /// updateCuts) exactly as setColumns does, but record into rollback only what
  /// a reject must undo: each changed cell's old code, or the whole pre-change
  /// column once more than maxJournal cells change (journaling then stops).
  /// Dense stores only: the bridge's refusePredictorMutation blocks this
  /// mutation surface on CSC-built samplers, so no CSC-backed column reaches
  /// here.
  void setColumnJournaled(size_t j, const double* newColumn, bool updateCuts,
                          size_t maxJournal, ColumnCodeRollback& rollback) {
    if (updateCuts) refreshCutsForColumn(j, newColumn);
    quantizeDenseObserved(train, numObservations, j, newColumn,
                          hasMissing.data(),
      [&](size_t i, const xint_t* column, xint_t code) {
        if (rollback.full || code == column[i]) return;
        rollback.journal.push_back({static_cast<std::uint32_t>(i), column[i]});
        if (rollback.journal.size() > maxJournal) {
          // reconstruct the pre-change column (current cells with the journaled
          // old codes restored) so the journal can be dropped from here on
          rollback.fullColumn.assign(column, column + numObservations);
          for (const ColumnCodeRollback::Cell& cell : rollback.journal)
            rollback.fullColumn[cell.index] = cell.oldCode;
          rollback.journal.clear();
          rollback.full = true;
        }
      });
    std::int32_t slot = gatheredSlotForColumn(j);
    if (slot >= 0)
      std::memcpy(gatheredRawValues.data() +
                    static_cast<size_t>(slot) * numObservations,
                  newColumn, numObservations * sizeof(double));
  }

  /// Undo a journaled re-quantize, restoring column j's codes from rollback.
  void restoreColumn(size_t j, const ColumnCodeRollback& rollback) {
    xint_t* column = train.codes.data() + train.codeOffsets[j];
    if (rollback.full)
      std::memcpy(column, rollback.fullColumn.data(),
                  numObservations * sizeof(xint_t));
    else
      for (const ColumnCodeRollback::Cell& cell : rollback.journal)
        column[cell.index] = cell.oldCode;
  }

  /// Overwrite a single cell's code against existing cuts, refreshing the
  /// gathered raw copy of a leaf-covariate column. A missing value marks the
  /// column; the flag only clears on a full column re-quantize (conservative
  /// but never wrong - the NA-aware partition handles NA-free columns too).
  void setCell(size_t i, size_t j, double value) {
    train.codes[train.codeOffsets[j] + i] = codeFor(j, value);
    if (isNA(value)) hasMissing[j] = 1;
    std::int32_t slot = gatheredSlotForColumn(j);
    if (slot >= 0)
      gatheredRawValues[static_cast<size_t>(slot) * numObservations + i] =
        value;
  }

  /// Whole-data replacement: new values for the same predictors, possibly a
  /// new number of observations, read for the call only. Ordinal cuts are
  /// rebuilt from scratch, so unlike refreshCutsForColumn a quantile-mode count
  /// may shrink and the caller remaps existing splits onto the new grid;
  /// categorical category counts stay fixed (the caller validates the new
  /// values). Dense stores only (setData is refused on CSC/mixed).
  void setData(const double* x_, size_t n) {
    numObservations = n;
    train.codes.resize(n * numPredictors);
    for (size_t j = 0; j < numPredictors; ++j) train.codeOffsets[j] = j * n;
    // resize the gathered leaf-covariate copies for the new observation count
    gatheredRawValues.assign(gatheredRawColumns.size() * n, 0.0);
    for (size_t j = 0; j < numPredictors; ++j) {
      const double* column = x_ + j * n;
      if (types[j] != ColumnType::categorical) buildCutsForColumn(j, column);
      quantizeColumn(j, column);
    }
  }

  void clearTest() {
    resetTestStorage();
  }

  /// Dense-stored columns only; rank columns have no contiguous codes.
  const xint_t* column(size_t variable) const { return train.column(variable); }

  bool columnIsSparse(size_t variable) const {
    return train.columnIsSparse(variable);
  }
  const SparseColumnData& sparseColumn(size_t variable) const {
    return train.sparseColumn(variable);
  }
  /// Storage-aware single-code access (tree descent, restore validation).
  xint_t codeAt(size_t variable, size_t i) const {
    return train.codeAt(variable, i);
  }

  bool testColumnIsSparse(size_t variable) const {
    return test.columnIsSparse(variable);
  }
  const SparseColumnData& testSparseColumn(size_t variable) const {
    return test.sparseColumn(variable);
  }
  /// Storage-aware single test-code access (test-row descent), reading only
  /// the columns a rule visits rather than materializing a row.
  xint_t testCodeAt(size_t variable, size_t i) const {
    return test.codeAt(variable, i);
  }
};

/// Temporarily install a donor cut grid over a store's ordinal columns for a
/// structural tree rebuild, restoring the live grid on scope exit. Only
/// cutPoints and the per-column cut counts move; the quantized observation
/// codes are the live data's and are never touched, so buildFromFlat resolves a
/// donor tree's split values against the donor grid while a later repartition
/// still routes the live observations. Categorical columns (empty donor cuts,
/// fixed category counts) pass through unchanged. Cross-grid warm start's
/// build-then-remap uses this (Chain::rebuildLiveForestRemapped); the caller
/// guarantees donorCutPoints has one entry per predictor, matching the store's
/// categorical/ordinal split column for column.
class ScopedCutGrid {
public:
  ScopedCutGrid(ColumnStore& store,
                const std::vector<std::vector<double>>& donorCutPoints)
      : store_(store), savedCutPoints_(store.cutPoints),
        savedNumCuts_(store.numCuts) {
    store_.cutPoints = donorCutPoints;
    for (size_t j = 0; j < store_.numPredictors; ++j)
      if (store_.types[j] != ColumnType::categorical)
        store_.numCuts[j] =
          static_cast<std::uint32_t>(donorCutPoints[j].size());
  }
  ~ScopedCutGrid() {
    store_.cutPoints = std::move(savedCutPoints_);
    store_.numCuts = std::move(savedNumCuts_);
  }
  ScopedCutGrid(const ScopedCutGrid&) = delete;
  ScopedCutGrid& operator=(const ScopedCutGrid&) = delete;

private:
  ColumnStore& store_;
  std::vector<std::vector<double>> savedCutPoints_;
  std::vector<std::uint32_t> savedNumCuts_;
};

}  // namespace bartcore

#endif  // BARTCORE_DATA_HPP
