#ifndef BARTCORE_DATA_HPP
#define BARTCORE_DATA_HPP

#include <algorithm>
#include <bit>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <type_traits>
#include <vector>

namespace bartcore {

using std::size_t;

/// Predictor code type: the per-cell quantized ordinal rank or categorical
/// code. uint16_t caps a column at maxNumCutsRepresentable cuts, keeping every
/// real code below naCode (0xFFFF), the reserved missing marker.
using xint_t = std::uint16_t;

/// Per-observation GATHER INDEX type for the hot index buffers (Tree::indices,
/// Forest/VarianceForest::indexBuffer): a subscript into the n-length columns,
/// narrowed from size_t to halve the single biggest hot array. Must match
/// the C kernel index type misc_index_t (pinned by a static_assert in
/// tree.hpp); DISTINCT
/// from length/count storage, which stays size_t. Capped by a guard that
/// refuses numObservations > UINT32_MAX at ingestion.
using index_t = std::uint32_t;

/// Reserved code for a missing ordinal value; cut counts cap below it so
/// real codes (which reach numCuts) never collide.
constexpr xint_t naCode = 0xFFFFu;
constexpr std::uint32_t maxNumCutsRepresentable = 0xFFFDu;

/// The missing marker of a HOST's int32 code channel (PredictorSource's
/// denseCodes), which is a different reservation from naCode: the minimum
/// int32, so every other value is a candidate level code and no real code is
/// spent. Distinct again from the code a missing value takes in the store,
/// missingCategoryCode below - the channel is what a host hands over, the
/// store's code is what quantization writes.
constexpr std::int32_t naDenseCode =
  std::numeric_limits<std::int32_t>::min();

/// The largest int32 code c for which (double)c <= \p value, widened to
/// int64 so both ends of the double range clamp rather than wrap. Exact:
/// every int32 converts to a double without rounding, so (double)c <= value
/// is c <= floor(value). A NaN value admits no code, which is what the
/// double comparison it stands in for does. The flat replay reads it once
/// per node so a coded column's threshold rule runs without converting a
/// cell.
inline std::int64_t codeThresholdBelow(double value) {
  // NaN, or under every int32: nothing is at or below it
  if (!(value >= -2147483649.0))
    return std::numeric_limits<std::int64_t>::min();
  if (value >= 2147483647.0) return 2147483647;
  return static_cast<std::int64_t>(std::floor(value));
}

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

/// One dense column's raw as ingestion reads it: the host's doubles, or its
/// int32 level codes. Exactly one channel is present, so isCoded discriminates
/// and at() serves either - a code widens, and naDenseCode becomes the NaN
/// every double-typed consumer already tests for. The code channel exists so a
/// host whose factor columns are integers hands them over as integers: the
/// widen-then-narrow round trip through a transient double block is what it
/// removes, and the codes it carries are the same integers that block held.
struct DenseColumnValues {
  const double* values = nullptr;
  const std::int32_t* codes = nullptr;

  DenseColumnValues() = default;
  // implicit from either channel, so a caller that has only one spells only
  // that one
  DenseColumnValues(const double* values_) : values(values_) {}
  DenseColumnValues(const std::int32_t* codes_) : codes(codes_) {}

  bool isCoded() const { return codes != nullptr; }
  bool isPresent() const { return values != nullptr || codes != nullptr; }
  double at(size_t i) const {
    if (codes == nullptr) return values[i];
    return codes[i] == naDenseCode
      ? std::numeric_limits<double>::quiet_NaN()
      : static_cast<double>(codes[i]);
  }
};

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

/// A column's semantic kind, which fixes how its values are read. A numeric
/// column quantizes real values against cut points; a categorical column
/// holds integer category codes 0..numCategories-1 directly; an ordered
/// factor holds level codes whose order is meaningful, so it quantizes
/// against cut points over those codes. The values are pinned: 0 and 1 are
/// the two the shipped C API has always carried, and the third is appended.
///
/// Splitting is the SEPARATE, derived axis - splitsBySubset below - because
/// only the two ends of the pipeline (grid construction, ingestion
/// validation, reporting) care which kind a column is, while every rule,
/// scan, mask and replay site cares only whether it partitions by subset
/// mask or by threshold. Masks of up to 63 categories live inline in the
/// rule word (and inline in the flattened node); wider columns pool their
/// mask words per tree and the flattened format references them through a
/// side channel. Codes must fit xint_t, including the reserved missing code
/// K of a pooled column.
enum class ColumnKind : std::uint8_t {
  numeric = 0,
  categorical = 1,
  orderedFactor = 2
};

/// The mechanic axis: whether a column's rules partition its codes by a
/// subset mask rather than by a threshold on their order. One definition,
/// read through ColumnStore::splitsBySubset, PredictorSource::splitsBySubset,
/// and directly wherever only a bare type channel is in hand.
constexpr bool kindSplitsBySubset(ColumnKind kind) {
  return kind == ColumnKind::categorical;
}

constexpr std::uint32_t maxCategories = 0xFFFFu;

/// The most levels a factor column of a kind can carry. A categorical column
/// spends codes 0..K-1 plus a missing position above them, which xint_t holds
/// while K reaches maxCategories. An ordered factor spends one more: the
/// upper bin of its K - 1 midpoint grid, so that grid must fit
/// maxNumCutsRepresentable and the kind's ceiling stops one lower. Read only
/// for a factor column.
constexpr std::uint32_t maxLevelsForKind(ColumnKind kind) {
  return kindSplitsBySubset(kind) ? maxCategories
                                  : maxNumCutsRepresentable + 1u;
}

/// The count a factor column's level-code sweep reports when some cell is not
/// a representable level code. Above every real count, which maxCategories
/// bounds.
constexpr std::uint32_t invalidCategoryCount = 0xFFFFFFFFu;

/// Whether a factor cell is a level code the store can hold: a whole number
/// in [0, maxCategories), or missing. The two kinds need it for different
/// reasons. A CATEGORICAL cell is narrowed to xint_t by a bare cast, which is
/// undefined outside that range, and a code past its column's count shifts
/// past an inline category mask or over-reads a pooled one. An ORDERED FACTOR
/// cell takes the ordinal arm instead, where the cast is safe and the reason
/// is the model: a value that is not one of the declared levels has no
/// position on a grid built from the level table. The R bridge refuses both
/// before the store is touched, so this exists for the header-only engine -
/// bartcore driven without that bridge - which has no other check between a
/// host's values and the store.
inline bool levelCodeIsRepresentable(double value) {
  return isNA(value) ||
         (value >= 0.0 && value < static_cast<double>(maxCategories) &&
          value == static_cast<double>(static_cast<std::uint32_t>(value)));
}

/// The same over the int32 code channel, where integrality is free and the
/// missing marker is naDenseCode rather than a NaN. Signed throughout: reading
/// naDenseCode as unsigned would make it the largest count in the column and
/// inflate K past the tier the mask machinery sizes from it.
inline bool levelCodeIsRepresentable(std::int32_t code) {
  return code == naDenseCode ||
         (code >= 0 && code < static_cast<std::int32_t>(maxCategories));
}

/// CSC-built columns at or below this nonzero fraction take rank-bitmap
/// storage; denser ones densify their codes.
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

/// Where one column's codes live and what re-quantization reads. The
/// discriminator is where the RAW a re-quantize reads lives, not who owns the
/// codes: every kind's codes are the store's.
enum class ColumnSourceKind : std::uint8_t {
  denseCallSupplied,  // dense codes in codes[]; the raw arrives with the call
                      // (train: the call-time x; test: none - the test store
                      // owns every raw it keeps)
  denseResident,      // dense codes in codes[]; the raw lives in the store, at
                      // residentRaw, a slice of the side's owned dense block
  denseCodesOnly,     // dense codes in codes[] and NO raw anywhere: a factor
                      // column's cells are level codes, which the codes
                      // already carry, and its grid follows the level table
                      // rather than its values, so it never re-quantizes
  cscRank,            // rank-bitmap in the side's sparseColumns[rankSlot];
                      // re-quantize from slice
  cscDensified        // dense codes in codes[]; re-quantize from slice
};

/// A per-column source descriptor: one instance per column of a row set,
/// carried in a vector sized to numPredictors on every built store. The kind
/// names WHICH POOL the column's raw lives in - the store's own double block
/// (denseResident), the caller's matrix (denseCallSupplied), the retained CSC
/// nonzeros (the two CSC kinds), or none at all (denseCodesOnly) - so "this
/// column has no double" is a state the descriptor states rather than a null a
/// reader must know to expect. The discriminated fields are read only for the
/// kinds that own them (rankSlot for cscRank; residentRaw for denseResident;
/// slice for the two CSC kinds; declaredCategoryCount for any factor column
/// whose host declared a level table, train-side only, and refCode for
/// CSC-backed categorical columns on both sides). denseCodesOnly reads none of
/// them.
struct ColumnSource {
  ColumnSourceKind kind = ColumnSourceKind::denseCallSupplied;
  std::int32_t rankSlot = -1;          // cscRank: slot into the side's sparseColumns
  // denseResident: the store's own dense column (train: ownedDenseValues,
  // test: ownedTestValues), writable so a mutation keeps the raw current
  double* residentRaw = nullptr;
  CscColumnSlice slice;                // cscRank/cscDensified: retained nonzeros
  std::uint32_t declaredCategoryCount = 0;  // factor: declared K (train only)
  xint_t refCode = 0;                  // CSC categorical: reference-level code

  bool isCscBacked() const {
    return kind == ColumnSourceKind::cscRank ||
           kind == ColumnSourceKind::cscDensified;
  }
};

/// One borrowed view of a host's predictor values, taken by creation, test
/// ingestion, and mutation alike: the storage (a column-major dense block, a
/// CSC triple, or any mix of the two) plus the per-column typing channel the
/// grid needs. The single spelling every entrance reads, so a validation or
/// ownership rule is stated once rather than per storage kind.
///
/// Per column: columnSources == nullptr means column j is dense column j of
/// denseValues, the call's transient block; a nonnegative entry names a dense
/// column of a host-pinned denseValues; a negative one names CSC column
/// ~columnSources[j] of the triple. sourceOf(j) answers uniformly, so no
/// consumer dereferences an absent map.
///
/// Ownership: everything here is borrowed for the consuming call only. A train
/// build over a MAPPED source retains what it must - the CSC slices stay
/// borrowed until a column is first mutated, and the REAL-VALUED dense columns
/// are copied into the store (ColumnStore::ownedDenseValues), so no host pins
/// a mixed build's block. A train build over an UNMAPPED source retains
/// nothing: the store re-quantizes from whatever matrix a later call hands it.
/// A test build copies every raw it KEEPS, so no borrowed pointer survives it.
/// Neither side keeps a FACTOR column's values: those are level codes, which
/// the store's own codes carry, so what a designated leaf covariate reads is a
/// gathered copy rather than a slice of either block.
///
/// The grid spec (cut counts, quantile mode) is deliberately absent: it is a
/// prior, not data, and a view built from a parent store carries the parent's.
struct PredictorSource {
  size_t numRows = 0;
  size_t numColumns = 0;

  const double* denseValues = nullptr;  // column-major, numRows x its columns
  // the code channel: a second column-major block holding int32 level codes,
  // naDenseCode for missing, for the dense columns a host stores as integers.
  // Read only for a FACTOR column - a numeric column's values are real numbers
  // whatever the host holds them in, and build widens a coded one back before
  // it cuts a grid over them.
  const std::int32_t* denseCodes = nullptr;
  // per column, where a dense-backed column's slice lives when the host splits
  // its dense storage across the two channels: k >= 0 names double column k of
  // denseValues, k < 0 code column ~k of denseCodes. The two channels index
  // independently, so each is packed over the columns it serves. Null leaves
  // every dense-backed column at column sourceOf(j) of denseValues, the single
  // block every entrance took before the code channel existed.
  const std::int32_t* denseChannels = nullptr;
  const int* cscColumnPointers = nullptr;  // length of the CSC source + 1
  const int* cscRowIndices = nullptr;
  const double* cscValues = nullptr;
  const std::int32_t* columnSources = nullptr;  // length numColumns

  // the typing channel: per-column types (null for all-ordinal), the level
  // count K each categorical column's host declares (null or 0 leaves it
  // inferred from the observed codes), and the reference-level code of each
  // CSC-backed categorical column (null when none is categorical)
  const ColumnKind* columnTypes = nullptr;
  const std::uint32_t* categoryCounts = nullptr;
  const xint_t* referenceCodes = nullptr;

  /// Whether the host supplied a per-column source map; without one the
  /// columns are dense column for column.
  bool isMapped() const { return columnSources != nullptr; }

  /// Column j's source in engine convention: a nonnegative dense column, or
  /// the complement of a CSC column. An absent map is the identity.
  std::int32_t sourceOf(size_t j) const {
    return columnSources != nullptr ? columnSources[j]
                                    : static_cast<std::int32_t>(j);
  }

  ColumnKind typeOf(size_t j) const {
    return columnTypes != nullptr ? columnTypes[j] : ColumnKind::numeric;
  }

  /// The mechanic axis over the view's own typing channel; see
  /// kindSplitsBySubset.
  bool splitsBySubset(size_t j) const { return kindSplitsBySubset(typeOf(j)); }

  std::uint32_t categoryCountOf(size_t j) const {
    return categoryCounts != nullptr ? categoryCounts[j] : 0u;
  }

  xint_t referenceCodeOf(size_t j) const {
    return referenceCodes != nullptr ? referenceCodes[j] : xint_t{0};
  }

  /// Whether the host split its dense storage across the two channels, so a
  /// dense-backed column's slice is addressed through denseChannels rather
  /// than as column sourceOf(j) of denseValues.
  bool hasSplitDenseChannels() const { return denseChannels != nullptr; }

  /// Whether column j's dense values arrive as int32 codes.
  bool denseColumnIsCoded(size_t j) const {
    return denseChannels != nullptr && sourceOf(j) >= 0 && denseChannels[j] < 0;
  }

  /// Column j's dense raw, in whichever channel holds it; both pointers null
  /// for a CSC-backed column and for a view carrying no values at all.
  DenseColumnValues denseColumn(size_t j) const {
    if (sourceOf(j) < 0) return {};
    if (denseChannels != nullptr) {
      std::int32_t which = denseChannels[j];
      if (which < 0)
        return denseCodes == nullptr
          ? DenseColumnValues{}
          : DenseColumnValues(denseCodes +
                              static_cast<size_t>(~which) * numRows);
      return denseValues == nullptr
        ? DenseColumnValues{}
        : DenseColumnValues(denseValues +
                            static_cast<size_t>(which) * numRows);
    }
    if (denseValues == nullptr) return {};
    return DenseColumnValues(denseValues +
                             static_cast<size_t>(sourceOf(j)) * numRows);
  }

  /// Whether the view is a plain column-major dense block every consumer can
  /// index as denseValues + j * numRows: values present, no CSC storage, and
  /// either no map or the identity one. A MAPPED dense view fails this - its
  /// columns need not sit at their own index - which is exactly why the
  /// mutation entrances test the block rather than merely the absence of CSC
  /// values.
  bool isDenseBlock() const {
    if (denseValues == nullptr || denseChannels != nullptr) return false;
    return isDenseColumnar();
  }

  /// The weaker question the entrances that consume a view COLUMN BY COLUMN
  /// ask: every column dense-backed, no CSC storage, no reordering map -
  /// whichever channel each column's values sit in. isDenseBlock is the
  /// stricter one, and only the kernels that index denseValues as a single
  /// block may ask it.
  bool isDenseColumnar() const {
    if (cscColumnPointers != nullptr || cscRowIndices != nullptr ||
        cscValues != nullptr)
      return false;
    for (size_t j = 0; columnSources != nullptr && j < numColumns; ++j)
      if (columnSources[j] != static_cast<std::int32_t>(j)) return false;
    return true;
  }
};

// Bridge entry points longjmp out of Rf_error past every destructor, so no
// view field may own storage.
static_assert(std::is_trivially_destructible_v<PredictorSource>,
              "PredictorSource must not own storage: Rf_error skips its dtor");

/// The dense spelling of a view: a plain column-major block, no map.
inline PredictorSource densePredictorSource(
    const double* values, size_t numRows, size_t numColumns,
    const ColumnKind* columnTypes = nullptr,
    const std::uint32_t* categoryCounts = nullptr) {
  PredictorSource source;
  source.numRows = numRows;
  source.numColumns = numColumns;
  source.denseValues = values;
  source.columnTypes = columnTypes;
  source.categoryCounts = categoryCounts;
  return source;
}

/// The view a sampler constructor builds from: the host's borrowed source
/// carrying the constructor's own row and column counts, plus - when the host
/// supplied no map - the constructor's transient dense block as the values a
/// null map reads.
inline PredictorSource creationPredictorSource(const PredictorSource& source,
                                               const double* x, size_t numRows,
                                               size_t numColumns) {
  PredictorSource view = source;
  view.numRows = numRows;
  view.numColumns = numColumns;
  // a host that split its dense storage across the two channels addresses
  // both through denseChannels, so the constructor's block is already the
  // double half of it rather than a replacement for the whole
  if (!view.isMapped() && !view.hasSplitDenseChannels()) view.denseValues = x;
  return view;
}

/// The dense-only creation view: the host's typing channel over the
/// constructor's dense block, any mapped or CSC storage dropped. The BCF and
/// multinomial samplers offer no sparse ingestion path. \p x is the block
/// only while the host kept ONE: a host that split its dense storage across
/// the two channels addresses both through denseChannels, and both ride from
/// \p source instead.
inline PredictorSource denseCreationPredictorSource(
    const PredictorSource& source, const double* x, size_t numRows,
    size_t numColumns) {
  PredictorSource view = densePredictorSource(
    x, numRows, numColumns, source.columnTypes, source.categoryCounts);
  // the code channel is dense storage, so it rides: dropping it would leave
  // the double half addressed as though it held every column
  view.denseCodes = source.denseCodes;
  view.denseChannels = source.denseChannels;
  if (view.hasSplitDenseChannels()) view.denseValues = source.denseValues;
  return view;
}

/// Whether every column of a source is CSC-backed - the bare-dgCMatrix design,
/// which serves no column's contiguous raw values and so admits no leaf
/// covariate at all. An unmapped source is dense column for column.
inline bool predictorSourceIsAllCsc(const PredictorSource& source,
                                    size_t numColumns) {
  if (!source.isMapped() || numColumns == 0) return false;
  for (size_t j = 0; j < numColumns; ++j)
    if (source.sourceOf(j) >= 0) return false;
  return true;
}

/// Materialize rows [rowBegin, rowEnd) of a borrowed view as a column-major
/// (rowEnd - rowBegin) x source.numColumns block of raw doubles: each column is
/// filled with its implicit value, then its stored entries scatter over it.
///
/// A CSC column's implicit rows read \p storeTypes[j] == categorical ?
/// referenceCodeOf(j) : 0 - the rule quantizeCscColumnInto applies, keyed on
/// the STORE's type rather than on whether the source declared a reference. A
/// null \p storeTypes means all-ordinal. \p out holds one block's worth of
/// doubles; the row range exists so a caller may stream a wide source in
/// chunks without materializing it whole.
inline void materializePredictorSource(const PredictorSource& source,
                                       const ColumnKind* storeTypes,
                                       size_t rowBegin, size_t rowEnd,
                                       double* out) {
  size_t numRows = rowEnd - rowBegin;
  for (size_t j = 0; j < source.numColumns; ++j) {
    double* target = out + j * numRows;
    std::int32_t which = source.sourceOf(j);
    if (which >= 0) {
      DenseColumnValues column = source.denseColumn(j);
      if (column.isCoded()) {
        // the one place the code channel meets a double-typed consumer: the
        // widen is total, naDenseCode becoming the NaN the reader tests for
        for (size_t i = 0; i < numRows; ++i)
          target[i] = column.at(rowBegin + i);
        continue;
      }
      std::memcpy(target, column.values + rowBegin,
                  numRows * sizeof(double));
      continue;
    }
    bool categorical =
      storeTypes != nullptr && kindSplitsBySubset(storeTypes[j]);
    double implicitValue =
      categorical ? static_cast<double>(source.referenceCodeOf(j)) : 0.0;
    for (size_t i = 0; i < numRows; ++i) target[i] = implicitValue;
    size_t column = static_cast<size_t>(~which);
    for (int k = source.cscColumnPointers[column];
         k < source.cscColumnPointers[column + 1]; ++k) {
      size_t row = static_cast<size_t>(source.cscRowIndices[k]);
      if (row >= rowBegin && row < rowEnd)
        target[row - rowBegin] = source.cscValues[k];
    }
  }
}

/// SparseColumnData's shape over RAW values: the same fixed pattern (bits,
/// wordRanks) and O(1) at(i), holding the borrowed nonzero doubles of one CSC
/// column and the double its implicit rows read. No cut grid and no
/// quantization, so nothing here touches a draw law - it is the read side of a
/// replay over a caller's sparse argument.
///
/// Built per call and never cached: a consumer reads the source it was handed,
/// so no store owns one and no mutation can leave one stale. nzValues borrows
/// the CSC values of its column, which are in ascending row order, so the rank
/// indexes the pattern rather than the values and duplicate values are no
/// different from distinct ones.
struct SparseRawColumn {
  std::vector<std::uint64_t> bits;       // ceil(n / 64) words
  std::vector<std::uint32_t> wordRanks;  // nonzeros before each word
  const double* nzValues = nullptr;      // borrowed, ascending-row order
  double implicitValue = 0.0;

  double at(size_t i) const {
    std::uint64_t word = bits[i >> 6];
    std::uint64_t bit = std::uint64_t{1} << (i & 63u);
    if ((word & bit) == 0) return implicitValue;
    return nzValues[wordRanks[i >> 6] +
                    static_cast<size_t>(std::popcount(word & (bit - 1u)))];
  }

  /// Lay the pattern of one CSC column's \p numNonzero ascending \p rows over
  /// \p numRows rows; \p values is borrowed for the reader's lifetime.
  void build(const int* rows, const double* values, size_t numNonzero,
             size_t numRows, double implicit) {
    size_t numWords = (numRows + 63) / 64;
    bits.assign(numWords, 0);
    wordRanks.assign(numWords, 0);
    nzValues = values;
    implicitValue = implicit;
    for (size_t k = 0; k < numNonzero; ++k) {
      size_t row = static_cast<size_t>(rows[k]);
      bits[row >> 6] |= std::uint64_t{1} << (row & 63u);
    }
    std::uint32_t runningRank = 0;
    for (size_t w = 0; w < numWords; ++w) {
      wordRanks[w] = runningRank;
      runningRank += static_cast<std::uint32_t>(std::popcount(bits[w]));
    }
  }
};

/// The reader a borrowed view hands the flat replay: one indexed load off a
/// dense-backed column, the rank lookup off a CSC-backed one. A dense column
/// is read through whichever channel holds it, so a host whose factor columns
/// are int32 level codes routes its rows straight off them. The code arm is
/// DenseColumnValues::at, the one widen in the header - naDenseCode becomes
/// the NaN the missing arms below test for - so a coded column routes exactly
/// as the same values laid out as doubles would, cell for cell.
struct PredictorSourceColumnReader {
  DenseColumnValues dense;
  const SparseRawColumn* sparse;
  double at(size_t row) const {
    if (dense.values != nullptr) return dense.values[row];
    if (dense.codes != nullptr) return dense.at(row);
    return sparse->at(row);
  }
  /// The column's own int32 codes, or null for one that is not stored in the
  /// code channel. The flat replay asks once per node and routes off them
  /// where they exist, so a coded column pays no conversion per row; every
  /// reader the replay accepts answers this.
  const std::int32_t* codes() const { return dense.codes; }
};

/// A borrowed PredictorSource in the Columns shape the flat replay reads
/// predictors through, so a sparse or mixed test set routes rows without a
/// dense n x p materialization. column(j) answers for STORE column j, which is
/// what the routing hoist and the leaf-covariate loads both index by.
///
/// A CSC-backed column's implicit rows read \p storeTypes[j] == categorical ?
/// referenceCodeOf(j) : 0, the same rule materializePredictorSource applies;
/// a null \p storeTypes means all-ordinal. A dense-backed column's slice is
/// whichever channel holds it, so a view whose host split its dense storage
/// across the two replays off the codes rather than off a block built to hold
/// them widened. Everything is built in the constructor and freed with the
/// object: the rank bitmaps cost O(numRows / 64) words per CSC-backed column
/// against the O(numRows) doubles a densification would cost, and the values
/// themselves stay borrowed.
struct PredictorSourceColumns {
  PredictorSourceColumns(const PredictorSource& source,
                         const ColumnKind* storeTypes) {
    size_t numRows = source.numRows;
    size_t numColumns = source.numColumns;
    // sized once, before any reader points into it, so no build reallocates;
    // an all-dense source allocates no column storage at all
    size_t numSparse = 0;
    for (size_t j = 0; j < numColumns; ++j)
      if (source.sourceOf(j) < 0) ++numSparse;
    sparseColumns_.resize(numSparse);
    readers_.resize(numColumns);
    numSparse = 0;
    for (size_t j = 0; j < numColumns; ++j) {
      std::int32_t which = source.sourceOf(j);
      if (which >= 0) {
        readers_[j] = { source.denseColumn(j), nullptr };
        continue;
      }
      size_t column = static_cast<size_t>(~which);
      int begin = source.cscColumnPointers[column];
      int end = source.cscColumnPointers[column + 1];
      bool categorical =
        storeTypes != nullptr && kindSplitsBySubset(storeTypes[j]);
      sparseColumns_[numSparse].build(
        source.cscRowIndices + begin, source.cscValues + begin,
        static_cast<size_t>(end - begin), numRows,
        categorical ? static_cast<double>(source.referenceCodeOf(j)) : 0.0);
      readers_[j] = { DenseColumnValues{}, &sparseColumns_[numSparse++] };
    }
  }

  PredictorSourceColumnReader column(size_t j) const { return readers_[j]; }

private:
  // one entry per CSC-backed column, packed in column order; the readers point
  // into it, so it is filled to its final size before any reader is published
  std::vector<SparseRawColumn> sparseColumns_;
  std::vector<PredictorSourceColumnReader> readers_;
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
  // storage kind, the rank slot into sparseColumns, the owned dense raw of a
  // mixed build, the retained CSC slice, a categorical column's host-declared
  // level count K, and a CSC-backed categorical column's reference-level code.
  // The declared count carries the levels no training row happens to observe,
  // which no scan of the codes can recover; a CSC categorical column stores
  // only its non-reference entries, so its implicit rows' code is the reference
  // level's own level-order index, not the storage's structural zero.
  std::vector<ColumnSource> sources;

  /// Whether column j quantizes from a borrowed CSC slice (every column of an
  /// all-CSC build, the sparse-source columns of a mixed one).
  bool columnIsCscBacked(size_t j) const { return sources[j].isCscBacked(); }

  bool columnIsSparse(size_t j) const {
    return sources[j].kind == ColumnSourceKind::cscRank;
  }
  const SparseColumnData& sparseColumn(size_t j) const {
    return sparseColumns[static_cast<size_t>(sources[j].rankSlot)];
  }

  /// Dense-stored columns only; rank columns have no contiguous codes. A rank
  /// column's codeOffsets entry is never assigned past the build's zero fill,
  /// so this would silently hand back whichever dense column starts at zero -
  /// or null, an all-rank block holding no codes at all. Every caller must
  /// have tested columnIsSparse; the assert is debug-only (R's build defines
  /// NDEBUG), so the descent and partition paths pay nothing for it.
  const xint_t* column(size_t j) const {
    assert(!columnIsSparse(j));
    return codes.data() + codeOffsets[j];
  }

  /// Storage-aware single-code access (tree descent, restore validation),
  /// reading only the columns a rule visits rather than materializing a row.
  xint_t codeAt(size_t j, size_t i) const {
    return columnIsSparse(j) ? sparseColumn(j).at(i)
                             : codes[codeOffsets[j] + i];
  }
};

/// Classic dense column store: borrowed column-major doubles quantized once
/// into per-column integer codes against per-column cut points, either
/// uniformly spaced over the column's range or at unique-value midpoints
/// (quantile mode). numCuts is fixed once built; recomputing cuts for new
/// values keeps the count so existing split indices stay in range. A
/// categorical column has no cut grid at all - numCuts 0, cutPoints empty,
/// its (fixed) category count in categoryCounts, and codes the values
/// themselves.
struct ColumnStore {
  // Move-only: a denseResident source points into this store's own
  // ownedDenseValues/ownedTestValues, so a copy would leave the duplicate's
  // sources aliasing the original's buffers. Moves keep the heap buffers at
  // their addresses, which is all a cached slice needs, and are all production
  // performs (the builders hand a finished store to the sampler).
  ColumnStore() = default;
  ColumnStore(const ColumnStore&) = delete;
  ColumnStore& operator=(const ColumnStore&) = delete;
  ColumnStore(ColumnStore&&) = default;
  ColumnStore& operator=(ColumnStore&&) = default;

  size_t numObservations = 0;
  size_t numPredictors = 0;
  bool useQuantiles = false;
  // Provenance: built from a parent store (buildFromParent). Gates only the
  // parent-derived standardization constants (suppliedStandardization);
  // refusal decisions read the capability predicates hasRequantizeSource and
  // acceptsNewRawPredictors instead.
  bool isView = false;

  std::vector<ColumnKind> types;
  // the training rows' codes over this store's cut grid
  CodeBlock train;
  bool builtFromCsc = false;
  bool hasSparse = false;
  std::vector<std::vector<double>> cutPoints;
  std::vector<std::uint32_t> numCuts;
  // per column, the level count K of a FACTOR column (either kind), 0 for a
  // numeric one. Fixed once built: every mask tier, every reserved missing
  // code and every category histogram width derives from it, so it is the one
  // channel that says how many levels a column has - numCuts answers only how
  // many cut points it has.
  std::vector<std::uint32_t> categoryCounts;
  std::vector<std::uint32_t> maxNumCuts;  // cap on quantile-induced counts
  // per column, whether any training value is missing; gates the extra
  // missing-direction draw in rules and the NA-aware partition kernel.
  // INVARIANT: a column flagged 0 consumes no missing-direction draw from the
  // rng and takes the plain partition path, so missingness support costs an
  // NA-free column neither a draw nor a branch of its own.
  std::vector<std::uint8_t> hasMissing;
  // categorical tier flag, fixed once category counts are (they never change
  // after build): pooled columns (> 63 levels) need the per-tree mask pool
  // and the flattened format's mask side channel; narrower masks are inline
  bool hasPooledCategorical = false;

  /// The mechanic axis: whether column j's rules partition by a subset mask
  /// over its category codes rather than by a threshold on their order. Every
  /// rule, scan, mask and replay site branches on this rather than on the
  /// kind, so a kind that splits by threshold needs no site of its own.
  bool splitsBySubset(size_t j) const { return kindSplitsBySubset(types[j]); }
  bool splitsByThreshold(size_t j) const { return !splitsBySubset(j); }

  /// The semantic axis, where the two factor kinds behave alike: their cells
  /// are level codes of a fixed table rather than measurements, so their grid
  /// follows the level count rather than the values and a value outside the
  /// table is malformed rather than merely extreme.
  bool isFactor(size_t j) const { return types[j] != ColumnKind::numeric; }

  /// Whether column j's rules store pool offsets instead of inline masks
  /// (more than 63 categories, so the mask spans more than one word).
  bool columnIsPooled(size_t j) const {
    return splitsBySubset(j) && maskWordsForCount(categoryCounts[j]) > 1;
  }

  void refreshCategoricalTiers() {
    hasPooledCategorical = false;
    for (size_t j = 0; j < numPredictors; ++j)
      if (columnIsPooled(j)) hasPooledCategorical = true;
  }

  size_t numTestObservations = 0;
  // Owned copy of the REAL-VALUED dense-backed test columns, column-major and
  // packed per predictor over the columns it serves, taken at buildTest; each
  // denseResident test source points into it. A FACTOR test column keeps only
  // its codes: its cells are level codes, its grid follows the training level
  // table rather than its values, so it never re-quantizes, and a designated
  // leaf covariate among them is served by the test gather below instead. The
  // engine borrows no test matrix either way: cut changes re-quantize a
  // numeric test column from this copy, and rawTestColumn serves it to leaf
  // models. Views hold none (they gather test-subset columns below).
  std::vector<double> ownedTestValues;
  // The test rows' codes over this store's cut grid (types, numCuts,
  // cutPoints), the training-side layout mirrored for the test rows. testCodeAt
  // reads whichever storage a column takes, so descent routes a test row
  // without materializing it. Views hold dense codes only (every test.sources
  // entry denseCallSupplied).
  CodeBlock test;
  // Owned raw of a mixed/CSC test build (the test store owns its raw, so no
  // borrowed pointer survives the build): ownedTestCscValues/ownedTestCscRows
  // pack every CSC-backed test column's nonzeros, which each test.sources entry's
  // slice points into. A dense-backed column's entry instead borrows into
  // ownedTestValues via residentRaw, and a CSC-backed categorical column's
  // entry carries the reference level's code. Both owned CSC buffers are empty
  // on the dense buildTest path and on views.
  std::vector<double> ownedTestCscValues;
  std::vector<int> ownedTestCscRows;
  // borrowed; added to recorded test fits. buildTest leaves it alone (the
  // caller keeps the lengths consistent), resetTestStorage clears it.
  const double* testOffset = nullptr;

  // Owned raw values of designated columns. build gathers a sampler's leaf
  // covariates (or a data handle's declared columns) so rawColumn serves owned
  // memory after the build borrow releases; buildFromParent gathers a view's
  // leaf covariates from its parent, with standardization constants from the
  // parent's full training data - the same calibration inheritance as the
  // copied cut grid, so every fold runs the prior a full-data fit would.
  std::vector<size_t> gatheredRawColumns;
  std::vector<double> gatheredRawValues;      // column-major, numObservations x q
  // the test-side twin, with its OWN column list because the two sides gather
  // different subsets: a top-level test store owns every real-valued dense
  // test column already and gathers only the designated columns it cannot
  // serve from that block, while a view gathers each side of the same
  // designation. Slots index gatheredRawTestColumns, never gatheredRawColumns.
  std::vector<size_t> gatheredRawTestColumns;
  std::vector<double> gatheredRawTestValues;  // column-major, numTestObservations x its columns
  std::vector<double> gatheredMeans;
  std::vector<double> gatheredSds;

  // Owned dense block of a mixed build, holding its REAL-VALUED dense-backed
  // columns and packed per predictor over the columns it serves; every
  // denseResident train source points into it. The store owns its raw on both
  // sides: the host assembles the block transiently and the build copies it,
  // so a mutation writes the new values through residentRaw and every later
  // reader - setCutPoints, state restore, a linear/GP leaf's regather - sees
  // the live column rather than the creation-time one. A FACTOR column keeps
  // no slice: its cells are level codes the codes already carry, its grid
  // follows the level table rather than its values so it never re-quantizes,
  // and a designated leaf covariate among them is gathered instead. Empty on
  // dense builds and views.
  std::vector<double> ownedDenseValues;

  // Owned re-quantize sources for CSC-backed columns after mutation. A column
  // keeps its build-time borrowed slice (R's dgCMatrix i/x slots) until first
  // mutated; setPredictor/updatePredictor writes the new nonzeros here and
  // repoints train.sources[j].slice at them, since the borrow no longer
  // reflects the live values (setCutPoints and state restore re-quantize from
  // the slice). Sized numPredictors on a CSC/mixed build, empty on dense
  // builds and views; ownedCsc*[j] stay empty and cscColumnOwned[j] false
  // until column j is first mutated.
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

  /// The same over the test-side gather, whose column list is its own.
  std::int32_t gatheredTestSlotForColumn(size_t j) const {
    for (size_t k = 0; k < gatheredRawTestColumns.size(); ++k)
      if (gatheredRawTestColumns[k] == j) return static_cast<std::int32_t>(k);
    return -1;
  }

  /// Whether column j quantizes from a borrowed CSC slice (every column of
  /// an all-CSC build, the sparse-source columns of a mixed one).
  bool columnIsCscBacked(size_t j) const {
    return train.columnIsCscBacked(j);
  }

  /// Whether any re-quantize source survives the build: retained CSC slices,
  /// the owned dense block, or a caller-supplied x. A view gathers codes from
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
  /// build's owned dense slice, x's column for a dense build, null for
  /// CSC-backed columns (which re-quantize from their retained slices instead)
  /// and for factor columns, which never re-quantize at all - their grid
  /// follows the level table, fixed at build.
  const double* rawColumnForRequantize(size_t j, const double* x) const {
    if (isFactor(j)) return nullptr;
    if (train.sources[j].isCscBacked()) return nullptr;
    if (train.sources[j].kind == ColumnSourceKind::denseResident)
      return train.sources[j].residentRaw;
    return x != nullptr ? x + j * numObservations : nullptr;
  }

  /// Owned raw training values of column j: a gathered copy (leaf covariates,
  /// or a data handle's declared columns), the mixed build's owned dense
  /// slice, null when neither serves it - a CSC-backed column, and a factor
  /// column the build was not asked to gather, which keeps only its codes.
  /// The designation is refused upstream where this is null, so a leaf model
  /// never reads one.
  const double* rawColumn(size_t j) const {
    std::int32_t slot = gatheredSlotForColumn(j);
    if (slot >= 0)
      return gatheredRawValues.data() +
             static_cast<size_t>(slot) * numObservations;
    return train.sources[j].kind == ColumnSourceKind::denseResident
      ? train.sources[j].residentRaw : nullptr;
  }

  /// Raw test values of column j: the test store's own dense slice, a gathered
  /// copy (a designated factor column, which retains no slice of its own, or a
  /// view's subset values), null when neither serves it - a CSC-backed column,
  /// which serves no leaf raw, and an undesignated factor column, which keeps
  /// only its codes. The designation is refused upstream where this is null,
  /// so a leaf model never reads one.
  const double* rawTestColumn(size_t j) const {
    if (!test.sources.empty()) {
      const ColumnSource& source = test.sources[j];
      if (source.kind == ColumnSourceKind::denseResident)
        return source.residentRaw;
      if (source.isCscBacked()) return nullptr;
    }
    std::int32_t slot = gatheredTestSlotForColumn(j);
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
    if (splitsBySubset(variable))
      return isNA(value) ? missingCategoryCode(categoryCounts[variable])
                         : static_cast<xint_t>(value);
    if (isNA(value)) return naCode;
    const std::vector<double>& cuts = cutPoints[variable];
    const double* first = cuts.data();
    return static_cast<xint_t>(
      std::lower_bound(first, first + numCuts[variable], value) - first);
  }

  /// A factor value of either kind is representable when it is an integral
  /// code of an existing level or missing; the level count is fixed once
  /// built.
  bool categoricalValueIsValid(size_t variable, double value) const {
    if (isNA(value)) return true;
    return value >= 0.0 &&
           value < static_cast<double>(categoryCounts[variable]) &&
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
    if (numUnique <= 1) {
      // A constant or fully missing column induces no interior split, but a
      // column with no cut at all is a state neither the store's own validator
      // nor its consumers admit; one degenerate cut at the observed value is
      // what a uniform grid places over the same column, and it splits
      // nothing for the same reason.
      grid.inducedNumCuts = 1;
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
    if (grid.sortedUnique.size() < 2) {  // degenerate: no midpoint to take
      double value = grid.sortedUnique.empty() ? 0.0 : grid.sortedUnique[0];
      std::fill(cutPoints[j].begin(), cutPoints[j].end(), value);
      return;
    }
    for (std::uint32_t k = 0; k < numCuts[j]; ++k) {
      size_t index = std::min(static_cast<size_t>(k) * grid.step + grid.offset,
                              grid.sortedUnique.size() - 2);
      cutPoints[j][k] =
        0.5 * (grid.sortedUnique[index] + grid.sortedUnique[index + 1]);
    }
  }

  /// An ordered factor's grid: K - 1 cuts at the midpoints between
  /// consecutive DECLARED level codes, so every adjacent level pair is
  /// separable, each distinct threshold occupies exactly one index slot (the
  /// grow prior's logCut then normalizes over the real partition set), and a
  /// level's code is its own index. The cap is RAISED to fit rather than the
  /// thinning bypassed, keeping numCuts[j] <= maxNumCuts[j] as every other
  /// externally determined grid does; n.cuts therefore does not apply to the
  /// column, exactly as it does not to a categorical one. Fewer than two
  /// levels admit no interior split and take the single degenerate cut a
  /// column with no cut at all cannot be. The raise needs no representability
  /// clamp of its own: the level count is bounded by maxLevelsForKind before
  /// the grid is built, so K - 1 never passes maxNumCutsRepresentable.
  void fillCutsAtLevelMidpoints(size_t j) {
    std::uint32_t numLevels = categoryCounts[j];
    numCuts[j] = numLevels >= 2u ? numLevels - 1u : 1u;
    if (maxNumCuts[j] < numCuts[j]) maxNumCuts[j] = numCuts[j];
    cutPoints[j].resize(numCuts[j]);
    for (std::uint32_t k = 0; k < numCuts[j]; ++k)
      cutPoints[j][k] = static_cast<double>(k) + 0.5;
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
  /// code plus one, and 0 when every value is missing. invalidCategoryCount
  /// when a cell is not a representable level code - the ingestion check rides
  /// this sweep, which already reads every cell of a factor column, so no cell
  /// is read twice and no column that is not a factor pays for it.
  std::uint32_t inferredCategoryCount(const double* column) const {
    double maxValue = -1.0;
    for (size_t i = 0; i < numObservations; ++i) {
      double value = column[i];
      if (!levelCodeIsRepresentable(value)) return invalidCategoryCount;
      if (!isNA(value) && value > maxValue) maxValue = value;
    }
    return maxValue < 0.0 ? 0u : static_cast<std::uint32_t>(maxValue) + 1u;
  }

  /// The same over the int32 code channel: naDenseCode is the missing marker
  /// rather than a NaN, and the running maximum stays SIGNED so the marker
  /// cannot become the largest code and inflate K past the mask tier sized
  /// from it.
  std::uint32_t inferredCategoryCount(const std::int32_t* column) const {
    std::int32_t maxCode = -1;
    for (size_t i = 0; i < numObservations; ++i) {
      std::int32_t code = column[i];
      if (!levelCodeIsRepresentable(code)) return invalidCategoryCount;
      if (code != naDenseCode && code > maxCode) maxCode = code;
    }
    return maxCode < 0 ? 0u : static_cast<std::uint32_t>(maxCode) + 1u;
  }

  /// The same over a CSC-backed column's logical values: its retained
  /// nonzeros, plus - for a subset-splitting column only - the reference code
  /// its implicit rows read. An ordered factor's implicit rows quantize a
  /// structural zero rather than the reference, so folding the reference into
  /// its count would inflate the grid and, at the ceiling, refuse a column
  /// nothing is wrong with. The reference is an xint_t and so needs no code
  /// check of its own; a reference at the top of the type is caught by the
  /// count ceiling instead, since the count it induces is one higher.
  std::uint32_t inferredCategoryCountCsc(size_t j) const {
    const ColumnSource& source = train.sources[j];
    double maxValue =
      splitsBySubset(j) && source.slice.numNonzero < numObservations
        ? static_cast<double>(source.refCode) : -1.0;
    for (size_t k = 0; k < source.slice.numNonzero; ++k) {
      double value = source.slice.values[k];
      if (!levelCodeIsRepresentable(value)) return invalidCategoryCount;
      if (!isNA(value) && value > maxValue) maxValue = value;
    }
    return maxValue < 0.0 ? 0u : static_cast<std::uint32_t>(maxValue) + 1u;
  }

  /// Initial cut construction; sets numCuts[j] and, for a factor column,
  /// categoryCounts[j]. column supplies the dense raw values (null for
  /// CSC-backed columns, which read their retained slice). A factor column of
  /// either kind takes its (fixed) level count from the host's declared level
  /// table when it supplied one, but never below what its own codes reach - a
  /// declared count short of an observed code would strand that code past its
  /// own grid. A categorical column then keeps no cuts at all, and an ordered
  /// factor takes the midpoint grid of that level table rather than a uniform
  /// or quantile grid over its observed values.
  ///
  /// keepCategoryCount holds a factor column's creation-time level count
  /// against a whole-data replacement, whose new values are a new sample of
  /// the SAME factor: the level table is a property of the factor rather than
  /// of the sample, so re-deriving it from the replacement would make K
  /// depend on call history.
  ///
  /// False refuses the column: some cell of a factor column is not a
  /// representable level code. keepCategoryCount takes the count as already
  /// checked, since it was checked when the column was built.
  bool buildCutsForColumn(size_t j, DenseColumnValues column,
                          bool keepCategoryCount = false) {
    if (isFactor(j) && !keepCategoryCount) {
      std::uint32_t inferred = columnIsCscBacked(j)
        ? inferredCategoryCountCsc(j)
        : (column.isCoded() ? inferredCategoryCount(column.codes)
                            : inferredCategoryCount(column.values));
      if (inferred == invalidCategoryCount) return false;
      std::uint32_t declared = train.sources[j].declaredCategoryCount;
      categoryCounts[j] = declared > inferred ? declared : inferred;
      // the count itself must be representable, not merely each code: a
      // declared count reaches the store without passing the sweep, and the
      // midpoint grid is sized from it
      if (categoryCounts[j] > maxLevelsForKind(types[j])) return false;
    }
    if (splitsBySubset(j)) {
      numCuts[j] = 0;
      cutPoints[j].clear();
    } else if (types[j] == ColumnKind::orderedFactor) {
      fillCutsAtLevelMidpoints(j);
    } else if (useQuantiles) {
      // a numeric column's grid is over its real values, which build widens a
      // coded column back to before it gets here
      QuantileGrid grid = columnIsCscBacked(j)
        ? quantileGridForCscColumn(j)
        : quantileGridForColumn(j, column.values);
      numCuts[j] = grid.inducedNumCuts;
      fillCutsFromQuantileGrid(j, grid);
    } else {
      numCuts[j] = maxNumCuts[j];
      if (columnIsCscBacked(j)) fillCutsUniformlyCsc(j);
      else fillCutsUniformly(j, column.values);
    }
    return true;
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
  /// A factor column of either kind has nothing to refresh: its grid follows
  /// its level table rather than its values, so re-deriving one from a
  /// replacement column would make the grid depend on call history - and for
  /// an ordered factor would replace the level midpoints with a uniform or
  /// quantile grid over the replacement's observed range, merging every level
  /// that range does not separate. The caller pre-checked the values against
  /// the table.
  bool refreshCutsForColumn(size_t j, const double* column) {
    if (isFactor(j)) return true;
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
  /// refresh feasibility for numeric columns, level-code validity for factor
  /// columns of either kind - whose grid is fixed by the level table, so the
  /// question is membership in it rather than an induced-cut-count comparison.
  bool cutsWouldRemainValid(size_t j, const double* values) const {
    if (isFactor(j)) {
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
    // a factor column's grid is the level table's, fixed at build: an
    // externally chosen one would strand its codes off their own levels, and
    // it retains no raw to re-quantize against one. Every reachable caller
    // refuses one first; this is the backstop
    if (isFactor(j)) return;
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

  /// The code channel's sibling: the same codes the double arm would write,
  /// reached without the widen. Two things must match it cell for cell. The
  /// missing marker is naDenseCode rather than a NaN, and it takes the SAME
  /// reserved code a NaN takes - narrowing it as though it were a level would
  /// make it a legal-looking category. And it must set hasMissingOut, since a
  /// column that goes from flagged to unflagged consumes no missing-direction
  /// draw and shifts every draw after it.
  void quantizeDenseCodesInto(CodeBlock& block, size_t numRows, size_t j,
                              const std::int32_t* raw,
                              std::uint8_t* hasMissingOut) {
    xint_t* column = block.codes.data() + block.codeOffsets[j];
    xint_t missingCode = splitsBySubset(j)
      ? missingCategoryCode(categoryCounts[j]) : naCode;
    std::uint8_t anyMissing = 0;
    for (size_t i = 0; i < numRows; ++i) {
      if (raw[i] == naDenseCode) {
        column[i] = missingCode;
        anyMissing = 1;
        continue;
      }
      column[i] = codeFor(j, static_cast<double>(raw[i]));
    }
    if (hasMissingOut != nullptr) hasMissingOut[j] = anyMissing;
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
    xint_t zeroCode =
      splitsBySubset(j) ? block.sources[j].refCode : codeFor(j, 0.0);
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
  void quantizeColumn(size_t j, DenseColumnValues column) {
    if (columnIsCscBacked(j)) {
      quantizeCscColumnInto(train, numObservations, j, hasMissing.data());
      return;
    }
    if (column.isCoded())
      quantizeDenseCodesInto(train, numObservations, j, column.codes,
                             hasMissing.data());
    else
      quantizeDenseInto(train, numObservations, j, column.values,
                        hasMissing.data());
    std::int32_t slot = gatheredSlotForColumn(j);
    if (slot < 0) return;
    // a leaf covariate reads raw doubles whatever channel the column arrived
    // in, so a coded one widens here rather than losing its gathered copy
    double* gathered =
      gatheredRawValues.data() + static_cast<size_t>(slot) * numObservations;
    if (!column.isCoded()) {
      std::memcpy(gathered, column.values, numObservations * sizeof(double));
      return;
    }
    for (size_t i = 0; i < numObservations; ++i) gathered[i] = column.at(i);
  }

  /// Re-quantize test column j against the current cuts. A CSC-backed column
  /// (the CSC-backed columns of a mixed test build) reads its retained owned
  /// slice; a real-valued dense-backed one reads its slice of ownedTestValues.
  /// A FACTOR test column has nothing to re-quantize from and needs nothing:
  /// its grid follows the training level table, which is fixed at build, so
  /// the codes buildTest wrote against it stay current. The test side gates no
  /// draws, so it tracks no missingness.
  void quantizeTestColumn(size_t j) {
    if (test.columnIsCscBacked(j)) {
      quantizeCscColumnInto(test, numTestObservations, j, nullptr);
      return;
    }
    const double* raw = test.sources[j].residentRaw;
    if (raw == nullptr) return;
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
    gatheredRawTestColumns.clear();
    gatheredRawTestValues.clear();
    gatheredMeans.clear();
    gatheredSds.clear();
  }

  /// Reset the per-column source storage to the dense-empty baseline: every
  /// column denseCallSupplied, no CSC or rank backing, no gathered raw, no
  /// recorded missingness. Every train builder resets through here, then
  /// overwrites the per-column source descriptors the view's storage kinds own
  /// (a mapped build the CSC slices, counts, and reference codes; an unmapped
  /// one the gathered raw via setupGatheredColumns). numPredictors must
  /// already hold the new column count; codes/codeOffsets and the cut grid are
  /// sized by each builder.
  void resetTrainStorage() {
    train.sources.assign(numPredictors, ColumnSource{});
    train.sparseColumns.clear();
    builtFromCsc = false;
    hasSparse = false;
    hasMissing.assign(numPredictors, 0);
    gatheredRawColumns.clear();
    gatheredRawValues.clear();
    gatheredRawTestColumns.clear();
    gatheredRawTestValues.clear();
    gatheredMeans.clear();
    gatheredSds.clear();
    // no CSC-backed columns until a CSC/mixed build re-sizes these; a dense
    // build leaves them empty and never reaches the mutation-owned path
    ownedCscRows.clear();
    ownedCscValues.clear();
    cscColumnOwned.clear();
    // likewise no owned dense block until a mapped build copies one
    ownedDenseValues.clear();
  }

  /// Build the training store from a borrowed predictor view (PredictorSource)
  /// against a fresh cut grid. maxNumCutsPerColumn, when non-null, overrides
  /// maxNumCutsScalar per column. The host validates STRUCTURE - CSC row
  /// indices unique and in range per column - since the quantize trusts it;
  /// factor level codes the build checks itself, below.
  ///
  /// An UNMAPPED view is the dense build: column j quantizes from dense column
  /// j of source.denseValues, read for the call only, and the store retains no
  /// raw - a later re-quantize reads whatever matrix its caller hands over.
  /// gatherColumns then names the columns whose raw values a leaf model (or a
  /// data handle's declared designation) must still see after the borrow
  /// releases, and the build owns a copy of each.
  ///
  /// A MAPPED view mixes storage. A nonnegative source names a dense column of
  /// the view, read through whichever value channel holds it and quantized
  /// with the same dense arithmetic, categorical allowed; a negative one names
  /// CSC column ~sourceOf(j) of the triple, which takes rank-bitmap storage at
  /// or below sparseDensityThreshold nonzero fraction and densified codes
  /// above, its borrowed slice serving re-quantization either way. The CSC
  /// triple stays borrowed for the store's lifetime (until a column's first
  /// mutation repoints it at owned nonzeros).
  ///
  /// Of the dense-backed columns, the REAL-VALUED ones are COPIED into
  /// ownedDenseValues, so the host's block need not outlive the call and a
  /// mutation writes the new values through the store's own copy; those serve
  /// rawColumn from that copy, so gatherColumns does not apply to them. A
  /// FACTOR column is copied nowhere - its cells are level codes the codes
  /// already carry, and its grid follows the level table, so it never
  /// re-quantizes - and serves rawColumn only when gatherColumns names it,
  /// which is how an ordered factor stays admissible as a leaf covariate here.
  /// A CSC-backed column is never gathered: sparse storage serves no dense raw
  /// at all.
  ///
  /// False REFUSES the build: some cell of a factor column is not a level code
  /// the store can represent. The store is left partly built and the caller
  /// discards it - a creation build has nothing to preserve - and the refusal
  /// travels out as a status rather than an exception, since the hosts that
  /// raise on it cross a C boundary.
  [[nodiscard]] bool build(const PredictorSource& source,
                           const std::uint32_t* maxNumCutsPerColumn,
                           std::uint32_t maxNumCutsScalar, bool useQuantiles_,
                           const size_t* gatherColumns = nullptr,
                           size_t numGatherColumns = 0) {
    size_t n = source.numRows, p = source.numColumns;
    bool mapped = source.isMapped();
    isView = false;
    numObservations = n;
    numPredictors = p;
    useQuantiles = useQuantiles_;
    if (source.columnTypes != nullptr) {
      types.assign(source.columnTypes, source.columnTypes + p);
    } else {
      types.assign(p, ColumnKind::numeric);
    }
    cutPoints.resize(p);
    numCuts.resize(p);
    categoryCounts.assign(p, 0);
    if (maxNumCutsPerColumn != nullptr) {
      maxNumCuts.assign(maxNumCutsPerColumn, maxNumCutsPerColumn + p);
    } else {
      maxNumCuts.assign(p, maxNumCutsScalar);
    }
    // keep the reserved missing code out of the real code range
    for (size_t j = 0; j < p; ++j)
      if (maxNumCuts[j] > maxNumCutsRepresentable)
        maxNumCuts[j] = maxNumCutsRepresentable;
    train.codeOffsets.assign(p, 0);
    resetTrainStorage();
    if (mapped) {
      // this branch owns the per-column CSC sources: borrowed slices,
      // densified/rank backing, and the reference codes of CSC-backed
      // categoricals (declared level counts ride both storage kinds, below)
      builtFromCsc = true;
      // per-column owned-mutation storage, empty until a column is first
      // mutated
      ownedCscRows.assign(p, {});
      ownedCscValues.assign(p, {});
      cscColumnOwned.assign(p, 0);
    }
    // the mapped arm's own dense block: one slice per predictor whose values
    // are real numbers, so a factor column costs it nothing. Creation peaks at
    // the host's transient assembly beside this copy of the columns it keeps
    std::vector<size_t> residentSlot(p, 0);
    size_t numResidentColumns = 0;
    if (mapped)
      for (size_t j = 0; j < p; ++j) {
        if (source.sourceOf(j) < 0 || isFactor(j)) continue;
        residentSlot[j] = numResidentColumns++;
      }
    ownedDenseValues.assign(numResidentColumns * n, 0.0);

    size_t numDenseCodes = 0;
    for (size_t j = 0; j < p; ++j) {
      ColumnSource& desc = train.sources[j];
      // the declared level table rides every storage kind: a dense-backed
      // factor inside a container declares one exactly as a CSC-backed one
      desc.declaredCategoryCount = source.categoryCountOf(j);
      std::int32_t columnSource = source.sourceOf(j);
      if (!mapped) {
        // dense build: codes packed column for column, re-quantized from
        // whatever matrix a later caller supplies
        train.codeOffsets[j] = numDenseCodes;
        numDenseCodes += n;
        continue;
      }
      if (columnSource >= 0) {
        if (isFactor(j)) {
          desc.kind = ColumnSourceKind::denseCodesOnly;
        } else {
          desc.kind = ColumnSourceKind::denseResident;
          // the block is sized above, so no later assign moves these
          desc.residentRaw = ownedDenseValues.data() + residentSlot[j] * n;
          // the copy widens a coded real-valued column exactly as a per-cell
          // read of it would, so the grid below is over real values whichever
          // channel the host held them in
          DenseColumnValues column = source.denseColumn(j);
          if (column.isCoded()) {
            for (size_t i = 0; i < n; ++i) desc.residentRaw[i] = column.at(i);
          } else {
            std::memcpy(desc.residentRaw, column.values, n * sizeof(double));
          }
        }
        train.codeOffsets[j] = numDenseCodes;
        numDenseCodes += n;
        continue;
      }
      size_t cscColumn = static_cast<size_t>(~columnSource);
      size_t begin = static_cast<size_t>(source.cscColumnPointers[cscColumn]);
      size_t end = static_cast<size_t>(source.cscColumnPointers[cscColumn + 1]);
      desc.slice = { source.cscValues + begin, source.cscRowIndices + begin,
                     end - begin };
      desc.refCode = source.referenceCodeOf(j);
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

    // the leaf-covariate gather, over the designated columns THIS build cannot
    // otherwise serve raw. A dense build retains no raw at all, so every
    // designated column is copied. A mapped build already serves its
    // real-valued dense-backed columns from its own block, so only a FACTOR
    // column among them is gathered - its cells are level codes, which the
    // block does not hold - and never a CSC-backed one, whose sparse storage
    // serves no dense raw and which the quantize would leave a slot of zeros.
    if (mapped) {
      std::vector<size_t> mappedGather;
      for (size_t k = 0; k < numGatherColumns; ++k) {
        size_t j = gatherColumns[k];
        if (j < p && isFactor(j) && source.sourceOf(j) >= 0)
          mappedGather.push_back(j);
      }
      setupGatheredColumns(mappedGather.data(), mappedGather.size());
    } else {
      setupGatheredColumns(gatherColumns, numGatherColumns);
    }

    // scratch for a coded column of NUMERIC kind, whose grid is over real
    // values: one column at a time, and only for a host that codes a column
    // the store does not read as a factor
    std::vector<double> widened;
    for (size_t j = 0; j < p; ++j) {
      buildRankStorageInto(train, n, j);
      // the raw a column quantizes from: the store's own dense slice where it
      // keeps one, the call's own channel otherwise - an unmapped column, and
      // a mapped FACTOR column, whose codes are its values - and neither for a
      // CSC-backed column (which reads its retained slice)
      DenseColumnValues column =
        train.sources[j].kind == ColumnSourceKind::denseResident
          ? DenseColumnValues(train.sources[j].residentRaw)
          : source.denseColumn(j);
      if (column.isCoded() && !isFactor(j)) {
        widened.resize(n);
        for (size_t i = 0; i < n; ++i) widened[i] = column.at(i);
        column = DenseColumnValues(widened.data());
      }
      if (!buildCutsForColumn(j, column)) return false;
      quantizeColumn(j, column);
    }
    refreshCategoricalTiers();

    resetTestStorage();
    return true;
  }

  /// Dense convenience spelling: a plain column-major matrix, no map.
  /// columnTypes may be null for all-ordinal; categoryCounts_, when supplied,
  /// is the host's declared level count per column (0 = infer from the codes).
  [[nodiscard]] bool build(const double* x_, size_t n, size_t p,
                           const std::uint32_t* maxNumCuts_, bool useQuantiles_,
                           const ColumnKind* columnTypes = nullptr,
                           const size_t* gatherColumns = nullptr,
                           size_t numGatherColumns = 0,
                           const std::uint32_t* categoryCounts_ = nullptr) {
    return build(densePredictorSource(x_, n, p, columnTypes, categoryCounts_),
                 maxNumCuts_, 0, useQuantiles_, gatherColumns,
                 numGatherColumns);
  }

  [[nodiscard]] bool build(const double* x_, size_t n, size_t p,
                           std::uint32_t maxNumCuts_,
                           bool useQuantiles_ = false,
                           const ColumnKind* columnTypes = nullptr,
                           const size_t* gatherColumns = nullptr,
                           size_t numGatherColumns = 0,
                           const std::uint32_t* categoryCounts_ = nullptr) {
    return build(densePredictorSource(x_, n, p, columnTypes, categoryCounts_),
                 nullptr, maxNumCuts_, useQuantiles_, gatherColumns,
                 numGatherColumns);
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
    gatheredRawTestColumns.clear();
    gatheredRawTestValues.clear();
    test.codes.clear();
    test.codeOffsets.clear();
    test.sources.clear();
    test.sparseColumns.clear();
    clearTestCscSources();
    testOffset = nullptr;
  }

  /// Whether every factor cell a test view would ingest is a level code of its
  /// column's table: each dense-backed column's slice, each CSC-backed
  /// column's stored nonzeros, and - for a subset-splitting column, whose
  /// implicit rows read it rather than a quantized zero - its reference code.
  /// A column the host declared a factor is a factor on both sides of the fit,
  /// so a value between two levels is refused rather than quantized onto the
  /// nearest boundary.
  bool testSourceLevelCodesAreValid(const PredictorSource& source) const {
    for (size_t j = 0; j < numPredictors; ++j) {
      if (!isFactor(j)) continue;
      std::int32_t columnSource = source.sourceOf(j);
      if (columnSource >= 0) {
        DenseColumnValues column = source.denseColumn(j);
        // buildTest would read through the same absence to copy the block
        if (!column.isPresent()) return false;
        for (size_t i = 0; i < source.numRows; ++i)
          if (!categoricalValueIsValid(j, column.at(i))) return false;
        continue;
      }
      size_t cscColumn = static_cast<size_t>(~columnSource);
      size_t end = static_cast<size_t>(source.cscColumnPointers[cscColumn + 1]);
      for (size_t k = static_cast<size_t>(source.cscColumnPointers[cscColumn]);
           k < end; ++k)
        if (!categoricalValueIsValid(j, source.cscValues[k])) return false;
      if (splitsBySubset(j) &&
          !categoricalValueIsValid(
            j, static_cast<double>(source.referenceCodeOf(j))))
        return false;
    }
    return true;
  }

  /// Whether buildTest would ingest this view, asked without touching the test
  /// store, so an entrance that has training values to install can decide its
  /// whole answer before it installs any of them. Keep this in step with
  /// buildTest's own opening test.
  bool testSourceIsIngestible(const PredictorSource& source) const {
    return testSourceLevelCodesAreValid(source);
  }

  /// Dense spelling of the same question, for the entrances that hand over a
  /// plain column-major test matrix.
  bool testSourceIsIngestible(const double* x_test, size_t numTest) const {
    return testSourceIsIngestible(
      densePredictorSource(x_test, numTest, numPredictors));
  }

  /// Whether every factor cell of a WHOLE-DATA replacement is a level code of
  /// its column's table. The table is the creation-time one, pinned rather
  /// than re-derived, so the question is membership in it rather than the
  /// representability the build's count sweep asks - and it is asked over the
  /// replacement's own row count, which need not be the store's.
  bool replacementLevelCodesAreValid(const double* x_, size_t n) const {
    for (size_t j = 0; j < numPredictors; ++j) {
      if (!isFactor(j)) continue;
      const double* column = x_ + j * n;
      for (size_t i = 0; i < n; ++i)
        if (!categoricalValueIsValid(j, column[i])) return false;
    }
    return true;
  }

  /// Build the test store from a borrowed predictor view against the training
  /// cut grid (already built, shared by identity; numCuts and cutPoints are not
  /// rebuilt). Column j reads dense column sourceOf(j) of the view - in
  /// whichever value channel holds it, so a host whose factor columns are
  /// int32 level codes hands them over as integers - or CSC column ~sourceOf(j)
  /// of the triple otherwise, which takes rank-bitmap storage at or below
  /// sparseDensityThreshold nonzero fraction and densified codes above, the
  /// training tier rule per column. An unmapped view is the plain test matrix,
  /// dense column for dense column. source.referenceCodes gives each
  /// CSC-backed categorical test column its reference level code (the code the
  /// implicit rows take).
  ///
  /// The test store owns every raw it KEEPS, whatever the view's shape, so no
  /// borrowed pointer survives the call: the real-valued dense columns are
  /// copied into ownedTestValues (packed per predictor over the columns it
  /// serves) and the CSC nonzero values+rows into their own buffers, and a
  /// later cut change re-quantizes from those copies. A FACTOR test column
  /// keeps only its codes - it never re-quantizes - and a designated leaf
  /// covariate among them is gathered instead, the test-side twin of the
  /// training gather.
  ///
  /// False REFUSES the view, leaving the test store untouched: some factor
  /// cell it would ingest is not a level code of that column's table. The
  /// table is the training one, fixed at build, so the check runs before
  /// anything is written rather than riding the quantize as the training
  /// side's does.
  [[nodiscard]] bool buildTest(const PredictorSource& source) {
    if (!testSourceLevelCodesAreValid(source)) return false;
    size_t p = numPredictors, numTest = source.numRows;
    numTestObservations = numTest;

    test.sources.assign(p, ColumnSource{});
    test.sparseColumns.clear();
    test.codeOffsets.assign(p, 0);

    // the dense block holds the real-valued columns only, one slice per
    // predictor that takes one; a factor column's slot stays absent
    std::vector<size_t> residentSlot(p, 0);
    size_t numResidentColumns = 0;
    for (size_t j = 0; j < p; ++j) {
      if (source.sourceOf(j) < 0 || isFactor(j)) continue;
      residentSlot[j] = numResidentColumns++;
    }
    ownedTestValues.assign(numResidentColumns * numTest, 0.0);

    // own the CSC nonzeros: pack them so each column's slice points into
    // storage that never reallocates for the store's lifetime
    size_t totalNonzero = 0;
    for (size_t j = 0; j < p; ++j)
      if (source.sourceOf(j) < 0) {
        size_t cscColumn = static_cast<size_t>(~source.sourceOf(j));
        totalNonzero +=
          static_cast<size_t>(source.cscColumnPointers[cscColumn + 1] -
                              source.cscColumnPointers[cscColumn]);
      }
    ownedTestCscValues.resize(totalNonzero);
    ownedTestCscRows.resize(totalNonzero);

    size_t numDenseCodes = 0, cursor = 0;
    for (size_t j = 0; j < p; ++j) {
      ColumnSource& desc = test.sources[j];
      std::int32_t columnSource = source.sourceOf(j);
      if (columnSource >= 0) {
        if (isFactor(j)) {
          desc.kind = ColumnSourceKind::denseCodesOnly;
        } else {
          desc.kind = ColumnSourceKind::denseResident;
          // the block is sized above, so no later assign moves these
          desc.residentRaw =
            ownedTestValues.data() + residentSlot[j] * numTest;
          // the copy widens a coded real-valued column exactly as a per-cell
          // read of it would
          DenseColumnValues column = source.denseColumn(j);
          if (column.isCoded()) {
            for (size_t i = 0; i < numTest; ++i)
              desc.residentRaw[i] = column.at(i);
          } else {
            std::memcpy(desc.residentRaw, column.values,
                        numTest * sizeof(double));
          }
        }
        test.codeOffsets[j] = numDenseCodes;
        numDenseCodes += numTest;
        continue;
      }
      size_t cscColumn = static_cast<size_t>(~columnSource);
      size_t begin = static_cast<size_t>(source.cscColumnPointers[cscColumn]);
      size_t end = static_cast<size_t>(source.cscColumnPointers[cscColumn + 1]);
      size_t numNonzero = end - begin;
      std::memcpy(ownedTestCscValues.data() + cursor,
                  source.cscValues + begin, numNonzero * sizeof(double));
      std::memcpy(ownedTestCscRows.data() + cursor,
                  source.cscRowIndices + begin, numNonzero * sizeof(int));
      desc.slice = { ownedTestCscValues.data() + cursor,
                     ownedTestCscRows.data() + cursor, numNonzero };
      desc.refCode = source.referenceCodeOf(j);
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
      if (test.columnIsCscBacked(j) ||
          test.sources[j].kind == ColumnSourceKind::denseResident) {
        quantizeTestColumn(j);
        continue;
      }
      // a factor test column retains no raw, so it quantizes here from the
      // borrowed channel that holds it - the one pass its codes ever need
      DenseColumnValues column = source.denseColumn(j);
      if (column.isCoded())
        quantizeDenseCodesInto(test, numTest, j, column.codes, nullptr);
      else
        quantizeDenseInto(test, numTest, j, column.values, nullptr);
    }

    // the test-side gather: a designated column the block above no longer
    // holds is copied out of the borrowed channel, so rawTestColumn serves a
    // leaf model owned memory once the borrow releases
    gatheredRawTestColumns.clear();
    gatheredRawTestValues.clear();
    for (size_t k = 0; k < gatheredRawColumns.size(); ++k) {
      size_t j = gatheredRawColumns[k];
      if (j >= p || source.sourceOf(j) < 0 ||
          test.sources[j].kind != ColumnSourceKind::denseCodesOnly)
        continue;
      size_t slot = gatheredRawTestColumns.size();
      gatheredRawTestColumns.push_back(j);
      gatheredRawTestValues.resize((slot + 1) * numTest);
      double* values = gatheredRawTestValues.data() + slot * numTest;
      DenseColumnValues column = source.denseColumn(j);
      for (size_t i = 0; i < numTest; ++i) values[i] = column.at(i);
    }
    return true;
  }

  /// Dense convenience spelling: own a copy of a plain column-major test
  /// matrix and quantize it against the current cuts.
  [[nodiscard]] bool buildTest(const double* x_test_, size_t numTest) {
    return buildTest(densePredictorSource(x_test_, numTest, numPredictors));
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
      categoryCounts.resize(numPredictors);
      maxNumCuts.resize(numPredictors);
      for (size_t j = 0; j < numPredictors; ++j) {
        types[j] = parent.types[parentColumns[j]];
        cutPoints[j] = parent.cutPoints[parentColumns[j]];
        numCuts[j] = parent.numCuts[parentColumns[j]];
        categoryCounts[j] = parent.categoryCounts[parentColumns[j]];
        maxNumCuts[j] = parent.maxNumCuts[parentColumns[j]];
      }
    } else {
      types = parent.types;
      cutPoints = parent.cutPoints;
      numCuts = parent.numCuts;
      categoryCounts = parent.categoryCounts;
      maxNumCuts = parent.maxNumCuts;
    }
    // views densify: gathered codes are fully dense whatever the parent's
    // per-column storage
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
      // a view gathers both sides of the same designation, so the two lists
      // agree slot for slot here even though a top-level store's do not
      gatheredRawTestColumns.push_back(rawColumnsToGather[k]);
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
      xint_t missingCode = splitsBySubset(j)
        ? missingCategoryCode(categoryCounts[j]) : naCode;
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

  // Mutation. The new raw values arrive as call arguments (a dense build keeps
  // no predictor matrix to write through; a mixed one writes them into its own
  // dense block and a CSC-backed column into its own nonzeros);
  // snapshot/rollback of cutPoints, codes, the gathered leaf raw, and the owned
  // dense raw is the caller's (the sampler's) responsibility. Cut refreshes
  // assume the caller pre-checked quantile feasibility with
  // cutsWouldRemainValid. Views hold no raw source and are refused upstream;
  // dense, CSC, and mixed builds all reach here.

  /// Keep a column's owned raw current with the values a mutation installs, so
  /// the re-quantize sources (setCutPoints, state restore) and the
  /// leaf-covariate regather read the live column. denseResident is the one
  /// kind that keeps a raw slice, and every other keeps none by construction
  /// rather than by accident: a dense build re-reads the caller's matrix, a
  /// FACTOR column's cells are the codes themselves, and a CSC-backed column's
  /// slice is rebuilt by its own mutation path. The CODES are never at stake
  /// here - every caller pairs this with the quantize that writes them - and a
  /// gathered leaf covariate is refreshed by that same quantize. The
  /// self-write guard covers the builders and re-quantizes that pass the owned
  /// column back in.
  void writeOwnedDenseColumn(size_t j, const double* column) {
    if (train.sources[j].kind != ColumnSourceKind::denseResident) return;
    double* raw = train.sources[j].residentRaw;
    if (raw == column) return;
    std::memcpy(raw, column, numObservations * sizeof(double));
  }

  /// The one-cell analogue, for the per-observation update session.
  void writeOwnedDenseCell(size_t i, size_t j, double value) {
    if (train.sources[j].kind != ColumnSourceKind::denseResident) return;
    train.sources[j].residentRaw[i] = value;
  }

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
      writeOwnedDenseColumn(j, column);
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
      writeOwnedDenseColumn(j, column);
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
  /// refreshes the numeric cut grid from the dense column exactly as the dense
  /// path does (the CSC cut builders fold the same implicit zeros the dense
  /// column carries explicitly, so the grids match); a factor column of either
  /// kind has no grid to refresh and keeps its creation-pinned level count.
  /// The caller snapshots for rollback first.
  void mutateCscColumnFromDense(size_t j, const double* column,
                                bool updateCuts) {
    if (updateCuts) refreshCutsForColumn(j, column);

    const double implicitValue = splitsBySubset(j)
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

  /// A transaction's snapshot of the owned dense raw of the columns it touches:
  /// the dense-backed ones among them and their pre-change values, packed
  /// column-major. Per column rather than per block so a per-sweep single-column
  /// update stays O(numObservations) rather than O(numObservations * p).
  struct OwnedDenseRollback {
    std::vector<size_t> columns;
    std::vector<double> values;
  };

  /// Snapshot the owned dense raw of the columns a transaction is about to
  /// write. columns names them (null means every predictor in order, the
  /// whole-matrix convention); non-dense-backed entries and stores with no
  /// owned block record nothing.
  void snapshotOwnedDenseColumns(const size_t* columns, size_t numColumns,
                                 OwnedDenseRollback& rollback) const {
    rollback.columns.clear();
    rollback.values.clear();
    if (ownedDenseValues.empty()) return;
    for (size_t k = 0; k < numColumns; ++k) {
      size_t j = columns != nullptr ? columns[k] : k;
      // only the columns that keep a slice; a factor column keeps none and
      // has no raw write to undo
      if (train.sources[j].kind != ColumnSourceKind::denseResident) continue;
      const double* raw = train.sources[j].residentRaw;
      rollback.columns.push_back(j);
      rollback.values.insert(rollback.values.end(), raw, raw + numObservations);
    }
  }

  /// Undo the raw writes of a rejected transaction, restoring each snapshotted
  /// column in place. The copy is a memcpy, never a buffer swap: every
  /// denseResident source caches a pointer into ownedDenseValues, so relocating
  /// it would dangle them all.
  void restoreOwnedDenseColumns(const OwnedDenseRollback& rollback) {
    for (size_t k = 0; k < rollback.columns.size(); ++k)
      std::memcpy(train.sources[rollback.columns[k]].residentRaw,
                  rollback.values.data() + k * numObservations,
                  numObservations * sizeof(double));
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
  /// Dense-STORED columns only - the subset transaction routes a CSC-backed one
  /// to the sparse mutation path instead - which includes the dense-backed
  /// columns of a mixed store, whose owned raw this writes through.
  void setColumnJournaled(size_t j, const double* newColumn, bool updateCuts,
                          size_t maxJournal, ColumnCodeRollback& rollback) {
    if (updateCuts) refreshCutsForColumn(j, newColumn);
    writeOwnedDenseColumn(j, newColumn);
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
    assert(!train.columnIsSparse(j));
    xint_t* column = train.codes.data() + train.codeOffsets[j];
    if (rollback.full)
      std::memcpy(column, rollback.fullColumn.data(),
                  numObservations * sizeof(xint_t));
    else
      for (const ColumnCodeRollback::Cell& cell : rollback.journal)
        column[cell.index] = cell.oldCode;
  }

  /// Overwrite a single cell's code against existing cuts, refreshing the owned
  /// dense raw and the gathered raw copy of a leaf-covariate column. A missing
  /// value marks the column; the flag only clears on a full column re-quantize
  /// (conservative but never wrong - the NA-aware partition handles NA-free
  /// columns too).
  void setCell(size_t i, size_t j, double value) {
    assert(!train.columnIsSparse(j));
    writeOwnedDenseCell(i, j, value);
    train.codes[train.codeOffsets[j] + i] = codeFor(j, value);
    if (isNA(value)) hasMissing[j] = 1;
    std::int32_t slot = gatheredSlotForColumn(j);
    if (slot >= 0)
      gatheredRawValues[static_cast<size_t>(slot) * numObservations + i] =
        value;
  }

  /// Whole-data replacement: new values for the same predictors, possibly a
  /// new number of observations, read for the call only. Threshold cuts are
  /// rebuilt from scratch, so unlike refreshCutsForColumn a quantile-mode count
  /// may shrink and the caller remaps existing splits onto the new grid; a
  /// factor column's level count stays fixed whichever kind it is. Dense
  /// stores only (setData is refused on CSC/mixed).
  ///
  /// False REFUSES the replacement, leaving the store untouched: some cell of
  /// a factor column is not a level code of that column's table. The count is
  /// pinned, so the count sweep build runs cannot arise here and the check is
  /// membership instead - and it runs before a single code is written, since
  /// unlike a creation build this store has values to preserve.
  [[nodiscard]] bool setData(const double* x_, size_t n) {
    if (!replacementLevelCodesAreValid(x_, n)) return false;
    numObservations = n;
    train.codes.resize(n * numPredictors);
    for (size_t j = 0; j < numPredictors; ++j) train.codeOffsets[j] = j * n;
    // resize the gathered leaf-covariate copies for the new observation count
    gatheredRawValues.assign(gatheredRawColumns.size() * n, 0.0);
    for (size_t j = 0; j < numPredictors; ++j) {
      const double* column = x_ + j * n;
      // the level count is pinned, so the count sweep does not run here; the
      // membership check above is what stands in for it
      if (splitsByThreshold(j)) (void)buildCutsForColumn(j, column, true);
      quantizeColumn(j, column);
    }
    return true;
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

  bool testColumnIsSparseForTesting(size_t variable) const {
    return test.columnIsSparse(variable);
  }
  const SparseColumnData& testSparseColumnForTesting(size_t variable) const {
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
      if (store_.splitsByThreshold(j))
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
