// Implementation of the flat C API (inst/include/dbarts/dbarts.h) over the
// bartcore engine; entry points are registered with R_RegisterCCallable in
// R_interface.cpp. Argument validation is minimal - consumers are compiled
// packages - except where an engine invariant is at stake (categorical
// category codes, column ranges).

#include <dbarts/dbarts.h>

#include <cstddef> // size_t
#include <cstdint> // uint64_t
#include <cstring> // memcpy
#include <type_traits> // is_same
#include <vector>

#include <external/Rinternals.h>

#include <misc/linearAlgebra.h> // misc_addVectorsInPlace

#include "R_interface_bartcore_common.hpp"

using std::size_t;
using bartcore_bridge::adoptVector;
using bartcore_bridge::callConvertingExceptions;
using bartcore_bridge::DrawCallbackProtection;
using bartcore_bridge::refuseCscReferenceAgainstStore;
using bartcore_bridge::refuseEmptyTreeStore;
using bartcore_bridge::refuseMultiForestResponseMutation;
using bartcore_bridge::refuseSparseLeafCovariate;
using bartcore_bridge::responseConduitIsFixed;
using bartcore_bridge::ResponseConduit;
using bartcore_bridge::sigmaIsPinned;
using bartcore_bridge::testFitsAreUndefined;
using bartcore_bridge::UnwindJump;
using bartcore_bridge::refuseNonFinite;
using bartcore_bridge::validateResponseSupport;
using bartcore_bridge::validateTestContainerAgainstStore;

namespace {

// dbarts_sampler_t IS the bridge's holder (R_interface_bartcore_common.hpp), so
// the address in an R sampler object's external pointer is the handle dbarts.h
// declares. Nothing is wrapped and nothing is owned here.
inline bartcore::SamplerBase& samplerOf(dbarts_sampler* sampler) {
  return *sampler->sampler;
}
inline const bartcore::SamplerBase& samplerOf(const dbarts_sampler* sampler) {
  return *sampler->sampler;
}


// Present-by-size read of a caller-filled struct's member: the input-side twin
// of the dbarts_results write guard, so a caller compiled against an older
// (smaller) layout is never READ past its own buffer either. Every optional
// member of dbarts_predictor_source below arrives through it.
#define SOURCE_PTR(source, field) \
  (DBARTS_HAS_FIELD(dbarts_predictor_source, (source), field) \
     ? (source)->field : NULL)
#define SOURCE_NUM(source, field) \
  (DBARTS_HAS_FIELD(dbarts_predictor_source, (source), field) \
     ? (source)->field : 0)

// A caller's predictor source, validated and translated: the engine's borrowed
// view, plus the STORE's column types indexed by the VIEW's own columns, which
// is what the CSC implicit-value rule and the reference refusal both key on
// (a view that names store columns of its own would index the two
// differently; the flat surface has no such entry).
struct TranslatedSource {
  bartcore::PredictorSource view;
  const bartcore::ColumnKind* storeTypes;
};

// Validate a caller's source against the sampler and translate it into the
// engine's view. \p columns maps view column j onto store column columns[j]
// (null for the identity, and range-checked by the caller before it gets
// here), \p numColumns is the width the entry requires, and \p numRowsRequired
// the row count it requires (0 leaves the rows the source's own, which is what
// the test-side entries take). Scratch rides R's transient stack, so the
// refusals below cost no cleanup on the way out.
//
// Every refusal here exists because the source declares its own shape: an
// argument that does not describe itself is refused rather than read to the
// sampler's width, which is the whole reason these entries take a struct.
TranslatedSource translateSource(const bartcore::ColumnStore& store,
                                 const dbarts_predictor_source* source,
                                 const size_t* columns, size_t numColumns,
                                 size_t numRowsRequired, const char* caller) {
  if (source == NULL)
    Rf_error("%s: the predictor source cannot be NULL", caller);
  // as dbarts_sampler_run's results do: a zero structSize means the caller
  // forgot to set it, and reading every member as absent would silently make a
  // null source of a populated one
  if (source->structSize == 0)
    Rf_error("%s: source.structSize is 0 - set it to "
             "sizeof(dbarts_predictor_source) (e.g. dbarts_predictor_source x "
             "= DBARTS_PREDICTOR_SOURCE_INIT)", caller);

  size_t numRows = SOURCE_NUM(source, numRows);
  if (SOURCE_NUM(source, numColumns) != numColumns)
    Rf_error("%s: the source declares %lu columns; %lu are required", caller,
             static_cast<unsigned long>(SOURCE_NUM(source, numColumns)),
             static_cast<unsigned long>(numColumns));
  if (numRowsRequired != 0 && numRows != numRowsRequired)
    Rf_error("%s: the source declares %lu rows; %lu are required", caller,
             static_cast<unsigned long>(numRows),
             static_cast<unsigned long>(numRowsRequired));

  const double* denseValues = SOURCE_PTR(source, denseValues);
  size_t numCscColumns = SOURCE_NUM(source, numCscColumns);
  const std::int32_t* cscColumnPointers = SOURCE_PTR(source, cscColumnPointers);
  const std::int32_t* cscRowIndices = SOURCE_PTR(source, cscRowIndices);
  const double* cscValues = SOURCE_PTR(source, cscValues);
  const std::int32_t* columnSources = SOURCE_PTR(source, columnSources);
  const std::int32_t* columnTypes = SOURCE_PTR(source, columnTypes);
  const std::uint32_t* categoryCounts = SOURCE_PTR(source, categoryCounts);
  const std::int32_t* referenceCodes = SOURCE_PTR(source, referenceCodes);
  const std::int32_t* denseCodes = SOURCE_PTR(source, denseCodes);
  size_t numDenseCodeColumns = SOURCE_NUM(source, numDenseCodeColumns);

  // the STORE's types, gathered onto the view's own columns; read below to
  // decide which channel each dense-backed column's index bounds against, so
  // it must precede the sweep rather than follow it
  bartcore::ColumnKind* storeTypes = reinterpret_cast<bartcore::ColumnKind*>(
    R_alloc(numColumns > 0 ? numColumns : 1, sizeof(bartcore::ColumnKind)));
  for (size_t j = 0; j < numColumns; ++j)
    storeTypes[j] = store.types[columns != NULL ? columns[j] : j];

  // per channel, whether any dense-backed column reads it: a channel a column
  // needs must be present, which is what the sweep below refuses on
  bool anyDoubleDense = false, anyCodedDense = false, anyCsc = false;
  for (size_t j = 0; j < numColumns; ++j) {
    if (columnTypes != NULL && columnTypes[j] != DBARTS_COLUMN_ORDINAL &&
        columnTypes[j] != DBARTS_COLUMN_CATEGORICAL &&
        columnTypes[j] != DBARTS_COLUMN_ORDERED_FACTOR)
      Rf_error("%s: source.columnTypes[%lu] is not one of "
               "DBARTS_COLUMN_ORDINAL, DBARTS_COLUMN_CATEGORICAL, "
               "DBARTS_COLUMN_ORDERED_FACTOR", caller,
               static_cast<unsigned long>(j));
    if (categoryCounts != NULL && categoryCounts[j] > bartcore::maxCategories)
      Rf_error("%s: source.categoryCounts[%lu] exceeds the %lu category limit",
               caller, static_cast<unsigned long>(j),
               static_cast<unsigned long>(bartcore::maxCategories));
    if (referenceCodes != NULL &&
        static_cast<std::int64_t>(referenceCodes[j]) >
          static_cast<std::int64_t>(bartcore::maxCategories))
      Rf_error("%s: source.referenceCodes[%lu] exceeds the %lu category limit",
               caller, static_cast<unsigned long>(j),
               static_cast<unsigned long>(bartcore::maxCategories));
    std::int32_t which =
      columnSources != NULL ? columnSources[j] : static_cast<std::int32_t>(j);
    if (which >= 0) {
      // one index, bounded against the channel the column's kind selects
      if (denseCodes != NULL &&
          storeTypes[j] != bartcore::ColumnKind::numeric) {
        if (static_cast<size_t>(which) >= numDenseCodeColumns)
          Rf_error("%s: source.columnSources[%lu] names code column %lu, but "
                   "the source declares %lu", caller,
                   static_cast<unsigned long>(j),
                   static_cast<unsigned long>(which),
                   static_cast<unsigned long>(numDenseCodeColumns));
        anyCodedDense = true;
      } else {
        if (static_cast<size_t>(which) >= numColumns)
          Rf_error("%s: source.columnSources[%lu] names dense column %lu, past "
                   "the source's own width", caller,
                   static_cast<unsigned long>(j),
                   static_cast<unsigned long>(which));
        anyDoubleDense = true;
      }
    } else {
      if (static_cast<size_t>(~which) >= numCscColumns)
        Rf_error("%s: source.columnSources[%lu] names CSC column %lu, but the "
                 "source declares %lu", caller, static_cast<unsigned long>(j),
                 static_cast<unsigned long>(~which),
                 static_cast<unsigned long>(numCscColumns));
      anyCsc = true;
    }
  }
  // a channel is required by the columns that read it, not by the other one
  // being present: a code channel beside a numeric dense column does not
  // supply that column's values
  if (anyDoubleDense && denseValues == NULL)
    Rf_error("%s: a dense-backed column names no denseValues", caller);
  if (anyCodedDense && denseCodes == NULL)
    Rf_error("%s: a factor column names no denseCodes", caller);
  if (anyCsc && (cscColumnPointers == NULL || cscRowIndices == NULL ||
                 cscValues == NULL))
    Rf_error("%s: a CSC-backed column names an incomplete CSC triple", caller);

  // one rule, one implementation: the bridge's own refusal, over the
  // per-CSC-column NA_INTEGER encoding it keys on (< 0 here is "declared
  // none", which is the absence a uint code cannot express)
  if (anyCsc && referenceCodes != NULL) {
    int* referenceMeta =
      reinterpret_cast<int*>(R_alloc(numCscColumns, sizeof(int)));
    for (size_t s = 0; s < numCscColumns; ++s) referenceMeta[s] = NA_INTEGER;
    for (size_t j = 0; j < numColumns; ++j) {
      if (columnSources[j] >= 0 || referenceCodes[j] < 0) continue;
      referenceMeta[static_cast<size_t>(~columnSources[j])] = referenceCodes[j];
    }
    refuseCscReferenceAgainstStore(storeTypes, columnSources, numColumns,
                                   referenceMeta, numCscColumns);
  }

  TranslatedSource translated;
  translated.storeTypes = storeTypes;
  translated.view.numRows = numRows;
  translated.view.numColumns = numColumns;
  translated.view.denseValues = denseValues;
  // published only when a column actually reads it, so a caller that leaves
  // stray CSC pointers beside an all-dense map keeps the dense fast path
  if (anyCsc) {
    translated.view.cscColumnPointers = cscColumnPointers;
    translated.view.cscRowIndices = cscRowIndices;
    translated.view.cscValues = cscValues;
  }
  translated.view.columnSources = columnSources;
  translated.view.categoryCounts = categoryCounts;
  if (referenceCodes != NULL) {
    bartcore::xint_t* codes = reinterpret_cast<bartcore::xint_t*>(R_alloc(
      numColumns > 0 ? numColumns : 1, sizeof(bartcore::xint_t)));
    for (size_t j = 0; j < numColumns; ++j)
      codes[j] = referenceCodes[j] >= 0
        ? static_cast<bartcore::xint_t>(referenceCodes[j]) : bartcore::xint_t{0};
    translated.view.referenceCodes = codes;
  }

  // The code channel, resolved against the STORE's kinds: a dense-backed
  // factor column reads its codes, everything else the double block, both at
  // the same columnSources[j] within the channel the kind selects. Only the
  // dense columns: any CSC storage stays sparse, so a coded source keeps
  // every rule an uncoded one has, the sparse leaf-covariate refusal
  // included. The channel rides the view from here: every consumer reads a
  // view a column at a time - the replay, the refusal sweeps, the test-store
  // build - so the codes stay where they lie.
  if (anyCodedDense) {
    std::int32_t* channels = reinterpret_cast<std::int32_t*>(
      R_alloc(numColumns > 0 ? numColumns : 1, sizeof(std::int32_t)));
    for (size_t j = 0; j < numColumns; ++j) {
      std::int32_t which = translated.view.sourceOf(j);
      channels[j] =
        which >= 0 && storeTypes[j] != bartcore::ColumnKind::numeric
          ? ~which : which;
    }
    translated.view.denseCodes = denseCodes;
    translated.view.denseChannels = channels;
  }
  return translated;
}

#undef SOURCE_PTR
#undef SOURCE_NUM

// A test NA takes a rule's learned missing direction, and a rule learns one
// only where the training column had NAs (ColumnStore::hasMissing gates the
// draw), so on a complete column it would take one fixed branch at every
// split. The R surface refuses first and names the column; this is the flat
// entrances' backstop, which has only the index.
void refuseTestMissingness(const bartcore::ColumnStore& store,
                           const bartcore::PredictorSource& source,
                           const char* caller) {
  // The reader owns heap storage (the CSC columns' rank bitmaps) and Rf_error
  // longjmps past its destructor, so the scan runs in a scope that closes
  // before the raise and reports the offending column by value.
  size_t offending = [&]() -> size_t {
    bartcore::PredictorSourceColumns columns(source, store.types.data());
    for (size_t j = 0; j < source.numColumns; ++j) {
      if (store.hasMissing[j]) continue;
      bartcore::PredictorSourceColumnReader column = columns.column(j);
      for (size_t i = 0; i < source.numRows; ++i)
        if (bartcore::isNA(column.at(i))) return j;
    }
    return source.numColumns;
  }();
  if (offending < source.numColumns)
    Rf_error("%s: test column %zu has missing values but the training "
             "column had none, so no rule routes them", caller, offending + 1);
}

// The test-side entries' shared refusals, in the order the R bridge runs them:
// a designated leaf covariate must be dense (CSC serves no contiguous raw),
// every categorical code is bounded against the STORE's counts, which the
// view's author cannot see, and a test NA is refused wherever the training
// column carried none.
void validateTestSource(const bartcore::SamplerBase& engine,
                        const TranslatedSource& source, const char* caller) {
  refuseSparseLeafCovariate(engine.shape(), source.view);
  validateTestContainerAgainstStore(engine.data(), source.view);
  refuseTestMissingness(engine.data(), source.view, caller);
}

} // namespace

// Layout lock for dbarts_results (dbarts.h): the growable ABI. Fields append
// monotonically and never reorder, so the library's offsetof matches every
// caller's. Every field's exact offset is pinned (all trailing members are
// pointer-sized), so a mid-struct insertion shifts a downstream offset and
// fails here, a reorder fails at the swapped pair, and the size assert forces
// an author who appends a field to update it (and, once 1.0-0 has shipped,
// bump DBARTS_C_API_MINOR).
static_assert(offsetof(dbarts_results, structSize) == 0);
static_assert(offsetof(dbarts_results, sigma) == sizeof(size_t) + 0 * sizeof(double*));
static_assert(offsetof(dbarts_results, train) == sizeof(size_t) + 1 * sizeof(double*));
static_assert(offsetof(dbarts_results, test) == sizeof(size_t) + 2 * sizeof(double*));
static_assert(offsetof(dbarts_results, varcount) == sizeof(size_t) + 3 * sizeof(double*));
static_assert(offsetof(dbarts_results, k) == sizeof(size_t) + 4 * sizeof(double*));
static_assert(offsetof(dbarts_results, varprobs) == sizeof(size_t) + 5 * sizeof(double*));
static_assert(offsetof(dbarts_results, logLikelihood) == sizeof(size_t) + 6 * sizeof(double*));
static_assert(offsetof(dbarts_results, shape) == sizeof(size_t) + 7 * sizeof(double*));
static_assert(offsetof(dbarts_results, residualDf) == sizeof(size_t) + 8 * sizeof(double*));
static_assert(sizeof(dbarts_results) == sizeof(size_t) + 9 * sizeof(double*),
              "dbarts_results layout changed; update these offsets, and bump "
              "DBARTS_C_API_MINOR if a field was appended after 1.0-0");

// The same lock for the two structs a CALLER fills, whose structSize is read
// against the library's offsets rather than written to the caller's. These
// asserts state the layout a reader can check by eye; the token fold below
// carries the same offsets to the consumer's runtime handshake.
static_assert(offsetof(dbarts_predictor_source, structSize) == 0);
static_assert(offsetof(dbarts_predictor_source, numRows) == 1 * sizeof(size_t));
static_assert(offsetof(dbarts_predictor_source, numColumns) == 2 * sizeof(size_t));
static_assert(offsetof(dbarts_predictor_source, denseValues) == 3 * sizeof(size_t));
static_assert(offsetof(dbarts_predictor_source, numCscColumns) ==
              3 * sizeof(size_t) + 1 * sizeof(double*));
static_assert(offsetof(dbarts_predictor_source, cscColumnPointers) ==
              4 * sizeof(size_t) + 1 * sizeof(double*));
static_assert(offsetof(dbarts_predictor_source, cscRowIndices) ==
              4 * sizeof(size_t) + 2 * sizeof(double*));
static_assert(offsetof(dbarts_predictor_source, cscValues) ==
              4 * sizeof(size_t) + 3 * sizeof(double*));
static_assert(offsetof(dbarts_predictor_source, columnSources) ==
              4 * sizeof(size_t) + 4 * sizeof(double*));
static_assert(offsetof(dbarts_predictor_source, columnTypes) ==
              4 * sizeof(size_t) + 5 * sizeof(double*));
static_assert(offsetof(dbarts_predictor_source, categoryCounts) ==
              4 * sizeof(size_t) + 6 * sizeof(double*));
static_assert(offsetof(dbarts_predictor_source, referenceCodes) ==
              4 * sizeof(size_t) + 7 * sizeof(double*));
static_assert(offsetof(dbarts_predictor_source, denseCodes) ==
              4 * sizeof(size_t) + 8 * sizeof(double*));
static_assert(offsetof(dbarts_predictor_source, numDenseCodeColumns) ==
              4 * sizeof(size_t) + 9 * sizeof(double*));
static_assert(sizeof(dbarts_predictor_source) ==
                5 * sizeof(size_t) + 9 * sizeof(double*),
              "dbarts_predictor_source layout changed; update these offsets");

// The same lock for the struct the library fills PER DRAW and the callback
// only reads. Only the boundaries are pinned here, since the token fold below
// carries every field's own name and offset and so fails on a reorder or an
// insertion by itself: what these add is the PORTABILITY of the layout. The
// eleven counts and ten channel pointers are one width and pack without
// padding; the four trailing scalars are pinned against sigma rather than
// against the prefix, because a host that aligns a double more strictly than a
// pointer pads before it - the one thing about this layout that is not the
// same everywhere, and the reason the fold does not carry their offsets.
static_assert(offsetof(dbarts_draw, structSize) == 0);
static_assert(offsetof(dbarts_draw, numOrdinalThresholds) ==
              10 * sizeof(size_t));
static_assert(offsetof(dbarts_draw, train) == 11 * sizeof(size_t));
static_assert(offsetof(dbarts_draw, varcount) ==
              11 * sizeof(size_t) + 9 * sizeof(double*));
static_assert(offsetof(dbarts_draw, k) ==
              offsetof(dbarts_draw, sigma) + 1 * sizeof(double));
static_assert(offsetof(dbarts_draw, shape) ==
              offsetof(dbarts_draw, sigma) + 2 * sizeof(double));
static_assert(offsetof(dbarts_draw, residualDf) ==
              offsetof(dbarts_draw, sigma) + 3 * sizeof(double));
static_assert(sizeof(dbarts_draw) ==
                offsetof(dbarts_draw, sigma) + 4 * sizeof(double),
              "dbarts_draw layout changed; update these offsets, and bump "
              "DBARTS_C_API_MINOR if a field was appended after 1.0-0");

// Compile-time ABI token, checked against the baked DBARTS_C_API_HASH: FNV-1a
// over the stringized DBARTS_C_API_LIST signatures, then the ABI enums'
// enumerator lists, then the layout the compiler gives the two structs that
// cross the ABI. Anything it
// covers moving fails this assert until DBARTS_C_API_HASH is re-baked - the
// mechanical acknowledgment that the ABI changed - and the consumer stubs
// raise on the mismatch at runtime until the consumer is rebuilt.
//
// To re-check that the assert still bites: flip one digit of
// DBARTS_C_API_HASH in inst/include/dbarts/dbarts.h and rebuild; the build
// must stop here.
namespace {
constexpr std::uint64_t dbarts_fnv1aPrime = 0x100000001b3ULL;
constexpr std::uint64_t dbarts_fnv1aBasis = 0xcbf29ce484222325ULL;

constexpr std::uint64_t dbarts_fnv1a(std::uint64_t hash, const char* text) {
  while (*text != '\0') {
    hash ^= static_cast<std::uint64_t>(static_cast<unsigned char>(*text));
    hash *= dbarts_fnv1aPrime;
    ++text;
  }
  return hash;
}
constexpr std::uint64_t dbarts_fnv1a(const char* text) {
  return dbarts_fnv1a(dbarts_fnv1aBasis, text);
}
// Integers enter the state a byte at a time, low byte first, extracted by
// SHIFTING and never by reading the object representation: the token is one
// number across hosts, so a big-endian one must fold the same bytes in the
// same order as a little-endian one.
constexpr std::uint64_t dbarts_fnv1aValue(std::uint64_t hash,
                                          std::uint64_t value) {
  for (int i = 0; i != 8; ++i) {
    hash ^= (value >> (8 * i)) & 0xffULL;
    hash *= dbarts_fnv1aPrime;
  }
  return hash;
}

// The ABI structs' fields, in declaration order, for the layout fold. Each
// field's NAME and OFFSET are folded from one token, so the two cannot drift
// apart and a rename moves the token as surely as a reorder does; offsets and
// sizes fold in POINTER UNITS, identical on ILP32, LP64 and LLP64, which is
// what keeps a platform out of the token. Every member of both is
// pointer-width, which the alignment asserts hold a future author to.
#define DBARTS_RESULTS_FIELDS(X) \
  X(structSize) X(sigma) X(train) X(test) X(varcount) X(k) X(varprobs) \
  X(logLikelihood) X(shape) X(residualDf)
#define DBARTS_PREDICTOR_SOURCE_FIELDS(X) \
  X(structSize) X(numRows) X(numColumns) X(denseValues) X(numCscColumns) \
  X(cscColumnPointers) X(cscRowIndices) X(cscValues) X(columnSources) \
  X(columnTypes) X(categoryCounts) X(referenceCodes) X(denseCodes) \
  X(numDenseCodeColumns)
// dbarts_draw splits in two. Its counts and channel pointers fold like the
// other structs' fields, being one width; its four by-value doubles fold by
// NAME, DECLARATION POSITION and sizeof instead, because a double is two
// pointer units on ILP32 and one on LP64 - an offset fold would make the token
// a different number per platform, and the header bakes ONE literal. Position
// plus sizeof still moves the token on a rename, a reorder, an insertion and a
// retype, which is what the offsets buy for the fields above.
#define DBARTS_DRAW_POINTER_FIELDS(X) \
  X(structSize) X(chainIndex) X(drawIndex) X(numObservations) \
  X(numTestObservations) X(numPredictors) X(numFittedValuesPerObservation) \
  X(numVariableCountForests) X(numForests) X(numAmplitudes) \
  X(numOrdinalThresholds) X(train) X(test) X(varianceFits) \
  X(varianceTestFits) X(forestFits) X(glue) X(splitProbabilities) \
  X(logLikelihood) X(ordinalThresholds) X(varcount)
#define DBARTS_DRAW_SCALAR_FIELDS(X) \
  X(sigma, 0) X(k, 1) X(shape, 2) X(residualDf, 3)
#define DBARTS_ALIGN_ASSERT(type, field) \
  static_assert(offsetof(type, field) % sizeof(void*) == 0, \
                "flat C API field is not pointer-aligned; the token folds " \
                "offsets in pointer units and would divide it away");
#define X(field) DBARTS_ALIGN_ASSERT(dbarts_results, field)
DBARTS_RESULTS_FIELDS(X)
#undef X
#define X(field) DBARTS_ALIGN_ASSERT(dbarts_predictor_source, field)
DBARTS_PREDICTOR_SOURCE_FIELDS(X)
#undef X
#define X(field) DBARTS_ALIGN_ASSERT(dbarts_draw, field)
DBARTS_DRAW_POINTER_FIELDS(X)
#undef X
#define DBARTS_FOLD_STRUCT(hash, type) \
  dbarts_fnv1aValue(dbarts_fnv1a(hash, #type), sizeof(type) / sizeof(void*))
#define DBARTS_FOLD_FIELD(hash, type, field) \
  dbarts_fnv1aValue(dbarts_fnv1a(hash, #field), \
                    offsetof(type, field) / sizeof(void*))
// name, declaration position, width in bytes: the platform-free fold, for a
// member whose offset is not the same number of pointer units everywhere.
#define DBARTS_FOLD_SCALAR(hash, type, field, position) \
  dbarts_fnv1aValue( \
    dbarts_fnv1aValue(dbarts_fnv1a(hash, #field), position), \
    sizeof(((type*) nullptr)->field))

constexpr std::uint64_t dbarts_drawNumFields = 0
#define X(field) + 1
  DBARTS_DRAW_POINTER_FIELDS(X)
#undef X
#define X(field, position) + 1
  DBARTS_DRAW_SCALAR_FIELDS(X)
#undef X
  ;

constexpr std::uint64_t dbarts_foldLayout(std::uint64_t hash) {
  hash = DBARTS_FOLD_STRUCT(hash, dbarts_results);
#define X(field) hash = DBARTS_FOLD_FIELD(hash, dbarts_results, field);
  DBARTS_RESULTS_FIELDS(X)
#undef X
  hash = DBARTS_FOLD_STRUCT(hash, dbarts_predictor_source);
#define X(field) hash = DBARTS_FOLD_FIELD(hash, dbarts_predictor_source, field);
  DBARTS_PREDICTOR_SOURCE_FIELDS(X)
#undef X
  // the per-draw struct: its name and FIELD COUNT rather than its sizeof, for
  // the same reason its scalars fold by position - the size is not one number
  // of pointer units across platforms once a by-value double is in it
  hash = dbarts_fnv1aValue(dbarts_fnv1a(hash, "dbarts_draw"),
                           dbarts_drawNumFields);
#define X(field) hash = DBARTS_FOLD_FIELD(hash, dbarts_draw, field);
  DBARTS_DRAW_POINTER_FIELDS(X)
#undef X
#define X(field, position) \
  hash = DBARTS_FOLD_SCALAR(hash, dbarts_draw, field, position);
  DBARTS_DRAW_SCALAR_FIELDS(X)
#undef X
  return hash;
}

// "NAME=VALUE;" per enumerator.
#define DBARTS_ENUMERATOR_TEXT(name, value) #name "=" #value ";"

// The signature half alone, baked privately so a failed build says WHICH half
// moved: both asserts firing means the entry-point list changed, the combined
// one alone means the ABI moved underneath unchanged signatures (a struct's
// layout, an enumerator).
constexpr std::uint64_t dbarts_apiSignatureToken =
  dbarts_fnv1a(DBARTS_C_API_DECLS);

constexpr std::uint64_t dbarts_apiToken() {
  std::uint64_t hash = dbarts_apiSignatureToken;
  hash = dbarts_fnv1a(hash, DBARTS_COLUMN_TYPE_LIST(DBARTS_ENUMERATOR_TEXT));
  hash = dbarts_fnv1a(hash, DBARTS_LEAF_MODEL_LIST(DBARTS_ENUMERATOR_TEXT));
  hash = dbarts_fnv1a(hash, DBARTS_FAMILY_LIST(DBARTS_ENUMERATOR_TEXT));
  return dbarts_foldLayout(hash);
}
} // namespace
static_assert(dbarts_apiSignatureToken == 0xfbf29fc67c22558bULL,
              "dbarts.h C API signatures moved (the entry-point list, not the "
              "layout fold); re-bake this literal here and DBARTS_C_API_HASH "
              "with it");
static_assert(dbarts_apiToken() == DBARTS_C_API_HASH,
              "dbarts.h C ABI changed - a signature, a struct's layout, or an "
              "ABI enumerator; re-bake DBARTS_C_API_HASH in "
              "inst/include/dbarts/dbarts.h (and bump DBARTS_C_API_MAJOR or "
              "DBARTS_C_API_MINOR as the change warrants)");

namespace {
// R's own dbartsSlowCountWarning, raised by evaluating warnOnSlowCount from the
// dbarts namespace on the carrier the R run builds, so the class, sentence and
// tally are the R route's.
void raiseSlowCountWarning(const bartcore::SamplerBase& sampler) {
  SEXP carrier = PROTECT(Rf_allocVector(VECSXP, 0));
  bartcore_bridge::attachSlowCountTally(carrier, sampler);
  SEXP name = PROTECT(Rf_mkString("dbarts"));
  SEXP ns = PROTECT(R_FindNamespace(name));
  SEXP call = PROTECT(Rf_lang2(Rf_install("warnOnSlowCount"), carrier));
  Rf_eval(call, ns);
  UNPROTECT(4);
}
} // namespace

extern "C" {

int dbarts_apiMajorVersion(void) { return DBARTS_C_API_MAJOR; }
int dbarts_apiMinorVersion(void) { return DBARTS_C_API_MINOR; }
uint64_t dbarts_apiHash(void) { return DBARTS_C_API_HASH; }

void dbarts_sampler_destroy(dbarts_sampler* sampler) {
  // The handle belongs to an R object, which frees the holder from its own
  // finalizer, so this releases the ENGINE and nothing else: the holder stays
  // addressable and reads as dead (bartcore_isValidPointer), which is the
  // state the R object's own methods re-create from a stored state or refuse
  // in. Idempotent by construction - resetting a null unique_ptr is a no-op -
  // which is what makes a second destroy safe where every other entry here
  // would dereference.
  if (sampler == NULL) return;
  sampler->sampler.reset();
}

void dbarts_sampler_run(dbarts_sampler* sampler, size_t numBurnIn,
                        size_t numSamples, dbarts_results* results) {
  // A zero structSize means the caller forgot to set it (see
  // DBARTS_RESULTS_INIT): reject loudly instead of silently skipping every
  // field and handing back an uninitialized buffer - the flat-API footgun that
  // fed garbage draws to a consumer's Gibbs loop. A nonzero older/smaller
  // structSize stays valid.
  if (results != NULL && results->structSize == 0)
    Rf_error("dbarts_sampler_run: results.structSize is 0 - set it to "
             "sizeof(dbarts_results) (e.g. dbarts_results r = "
             "DBARTS_RESULTS_INIT)");

  // The run's two error paths meet here. A registered draw callback that
  // raises jumps at the LEAF - the one call into the host - which converts it
  // into an UnwindJump, so the engine's whole run unwinds with its destructors
  // before this frame catches it and hands the jump back to R. An engine
  // failure arrives as an ordinary exception and becomes an R error, raised
  // below where nothing is in flight.
  bartcore_bridge::CapturedError error;
  SEXP continuation = NULL;
  bool cancelled = false, stoppedByCallback = false;
  try {
    DrawCallbackProtection armed(sampler->drawHook);
    bartcore_bridge::captureExceptions(error, [&]() {
      bartcore::SamplerShape shape = samplerOf(sampler).shape();
      bartcore::Results engineResults;
      // the two widths the header states through its accessors: F fitted
      // values per observation (K on multinomial, 1 elsewhere) for train and
      // test, and V split-count sets per draw for varcount, exactly as the R
      // run sets them, so a caller sizing by the accessors gets R's layout
      engineResults.numReportedLocations = shape.numReportedLocations;
      engineResults.numVariableCountForests = shape.numVariableCountForests;

      if (results != NULL && numSamples > 0) {
        // A field is filled only when present-by-size AND non-null. offsetof
        // is against the library's (newest) layout; fields only append, so it
        // equals the caller's offset and structSize bounds the buffer.
#define FILL(field, member) \
  engineResults.member = DBARTS_RESULTS_HAS(results, field) ? results->field : NULL
        FILL(sigma, sigma);
        FILL(train, trainingFits);
        FILL(test, testFits);
        FILL(varcount, variableCounts);
        FILL(k, k);
        FILL(varprobs, splitProbabilities);
        FILL(logLikelihood, logLikelihood);
        FILL(shape, shape);
        FILL(residualDf, residualDf);
#undef FILL
        // the header leaves k untouched without a k hyperprior, as the R run
        // does; the engine writes it wherever the pointer is non-null
        if (!shape.kIsSampled) engineResults.k = NULL;
      }

      // The engine samples only from each chain's own Mersenne Twister (seeded
      // from R's stream once at creation), never from R's stream during a run,
      // so no GetRNGstate/PutRNGstate bracket is needed here - and none is
      // left unbalanced by a longjmp out of the engine.
      // the registered observer, adapted to the shipped draw struct one draw
      // at a time; an empty hook when nothing is registered. The interrupt
      // poll is the R route's, so a pending interrupt cancels the run at the
      // next sweep boundary or inside a leaf-order count, and the engine says
      // whether a cancel was the callback's stop rather than an interrupt.
      cancelled = samplerOf(sampler).run(
        numBurnIn, numSamples, engineResults, bartcore_bridge::userInterrupted,
        {}, sampler->drawHook.engineHook(), &stoppedByCallback);
    });
  } catch (const UnwindJump& jump) {
    // the protection is disarmed by now: the throw ran every destructor
    // between the callback and here, this frame's guard included
    continuation = jump.continuation;
  }
  // the handler is left before the jump resumes: a longjmp out of a live catch
  // block strands the exception on this thread's caught-exception stack, the
  // same reason captureExceptions copies its message out before raising
  if (continuation != NULL) R_ContinueUnwind(continuation); // does not return
  if (error.failed) Rf_error("dbarts_sampler_run: %s", error.message);
  // a callback's stop keeps its contract and returns normally; a real cancel
  // leaves the sampler as that stop does and raises
  if (cancelled && !stoppedByCallback)
    Rf_error("dbarts_sampler_run: sampler run interrupted");
  // once per holder, last, with nothing of the library's live: any handler
  // that exits on the warning jumps out of a complete run
  if (!sampler->slowCountWarned &&
      samplerOf(sampler).slowCountTally().slowCounts > 0) {
    sampler->slowCountWarned = true;
    raiseSlowCountWarning(samplerOf(sampler));
  }
}

/// The setter copies the pair into the sampler and nothing else: no call is
/// made through fn here, and the sampler neither reads nor frees the context.
/// A null fn clears, dropping the context with it.
void dbarts_sampler_setDrawCallback(dbarts_sampler* sampler,
                                    dbarts_draw_callback fn, void* context) {
  sampler->drawHook.set(fn, context);
}

void dbarts_sampler_sampleTreesFromPrior(dbarts_sampler* sampler) {
  // The prior draw refuses from the bottom of a deep call stack (a tree whose
  // every draw left an empty leaf), and growing a whole forest allocates. Both
  // reach here as C++ exceptions, so the engine's frames unwind first and only
  // this one raises.
  //
  // draws from the chain RNG only, not R's stream (see dbarts_sampler_run)
  callConvertingExceptions("dbarts_sampler_sampleTreesFromPrior", [&]() {
    samplerOf(sampler).sampleTreesFromPrior();
  });
}

int dbarts_sampler_setResponse(dbarts_sampler* sampler, const double* y,
                               int updateScale) {
  // the capability answer, which no argument would have changed: a coupling
  // that fixes its response conduit at creation
  if (responseConduitIsFixed(samplerOf(sampler).shape())) return 0;
  // the shared conduit guard, not the whole-data refusal: a two-forest
  // sampler's response swap is opt-in and scale-pinned rather than refused -
  // the same rule bartcore_setResponse applies
  refuseMultiForestResponseMutation(samplerOf(sampler),
                                    "dbarts_sampler_setResponse",
                                    ResponseConduit::response, updateScale);
  // the one place minimal validation is not enough: an out-of-support y is a
  // silently garbage latent draw for probit/ordinal and, for nbinom, an
  // uncatchable crash inside the count histogram (see validateResponseSupport)
  bartcore::SamplerShape shape = samplerOf(sampler).shape();
  refuseNonFinite(y, shape.numObservations, "dbarts_sampler_setResponse",
                  "response");
  validateResponseSupport(shape.family, shape.numOrdinalThresholds + 1, y,
                          shape.numObservations, "dbarts_sampler_setResponse");
  // the probit latent redraw draws from the chain RNG, not R's stream
  samplerOf(sampler).setResponse(
    adoptVector(sampler->ownedResponse, y, shape.numObservations),
    updateScale != 0);
  return 1;
}

int dbarts_sampler_setOffset(dbarts_sampler* sampler, const double* offset,
                             int updateScale) {
  // the offset is the response-side swap under a different pointer, so it
  // carries the same conditions; see dbarts_sampler_setResponse
  if (responseConduitIsFixed(samplerOf(sampler).shape())) return 0;
  refuseMultiForestResponseMutation(samplerOf(sampler),
                                    "dbarts_sampler_setOffset",
                                    ResponseConduit::offset, updateScale);
  refuseNonFinite(offset, samplerOf(sampler).shape().numObservations,
                  "dbarts_sampler_setOffset", "offset");
  samplerOf(sampler).setOffset(
    adoptVector(sampler->ownedOffset, offset,
                samplerOf(sampler).shape().numObservations),
    updateScale != 0);
  return 1;
}

int dbarts_sampler_setSigma(dbarts_sampler* sampler, double sigma) {
  // reachable here: a handle can name any sampler R builds, the pinned
  // families (probit, logistic, ordinal, nbinom, multinomial) and the
  // heteroscedastic gaussian dbartsSpec(variance = ) builds included
  if (sigmaIsPinned(samplerOf(sampler))) return 0;
  // a non-finite or non-positive sd poisons every later draw, as a
  // non-finite response does
  if (!R_finite(sigma) || sigma <= 0.0)
    Rf_error("dbarts_sampler_setSigma: sigma must be finite and positive");
  samplerOf(sampler).setSigma(sigma);
  return 1;
}

int dbarts_sampler_getLatents(const dbarts_sampler* sampler, double* out) {
  const bartcore::SamplerBase& engine(samplerOf(sampler));
  if (engine.latents(0) == NULL) return 0;

  bartcore::SamplerShape shape = engine.shape();
  size_t numObservations = shape.numObservations;
  for (size_t c = 0; c < shape.numChains; ++c)
    std::memcpy(out + c * numObservations, engine.latents(c),
                numObservations * sizeof(double));
  return 1;
}

int dbarts_sampler_predict(dbarts_sampler* sampler,
                           const dbarts_predictor_source* xTest,
                           const double* offsetTest, size_t numThreads,
                           double* out) {
  bartcore::SamplerBase& engine(samplerOf(sampler));
  // predictColumns opens forests_[0] alone, so a caller would receive the
  // first forest's fit labelled as the whole; see
  // dbarts_sampler_setTestPredictors
  if (testFitsAreUndefined(engine)) return 0;
  // A multi-location (multinomial) surface reads its offset as numRows x F,
  // before the softmax, under R's two refusals. Both read only the shape, the
  // holder's own category offsets and the caller's buffer, so they raise here,
  // before anything is captured or allocated; a malformed source is left to
  // translateSource below.
  size_t numLocations = engine.shape().numReportedLocations;
  if (numLocations > 1) {
    if (offsetTest == NULL) {
      if (!sampler->ownedCategoryOffset.empty() ||
          !sampler->ownedCategoryTestOffset.empty())
        Rf_error("dbarts_sampler_predict: %s",
                 bartcore_bridge::categoryOffsetRequiredMessage);
    } else if (xTest != NULL && xTest->structSize != 0) {
      size_t numRows =
        DBARTS_HAS_FIELD(dbarts_predictor_source, xTest, numRows)
          ? xTest->numRows : 0;
      size_t numOffsets = numRows * numLocations;
      for (size_t i = 0; i < numOffsets; ++i)
        if (!R_finite(offsetTest[i]))
          Rf_error("dbarts_sampler_predict: offset %s",
                   bartcore_bridge::categoryOffsetNotFiniteMessage);
    }
  }
  // The replay builds the CSC rank bitmaps a sparse view reads through and
  // fans across threads, so a worker failure arrives as a C++ exception whose
  // unwind frees both before this frame reports it to R.
  int filled = 0;
  callConvertingExceptions("dbarts_sampler_predict", [&]() {
    refuseEmptyTreeStore(engine, "dbarts_sampler_predict");
    bartcore::SamplerShape shape = engine.shape();
    void* scratch = vmaxget();
    TranslatedSource source = translateSource(
      engine.data(), xTest, NULL, shape.numPredictors, 0,
      "dbarts_sampler_predict");
    // a read-only replay builds no store, so the leaf-covariate rule is
    // checked on the view itself rather than answered by a store build
    validateTestSource(engine, source, "dbarts_sampler_predict");
    size_t numTestObservations = source.view.numRows;

    // F > 1: the offset is the per-category matrix the replay adds to each
    // category's raw fits; F = 1: a vector added to the fits after
    engine.predict(source.view, numTestObservations,
                   numLocations > 1 ? offsetTest : NULL, numThreads, out);
    vmaxset(scratch);

    if (offsetTest != NULL && numLocations == 1) {
      size_t capacity = shape.savedTreeCapacity;
      size_t numSamples = capacity > 0 ? shape.numSavedDraws : 1;
      for (size_t slab = 0; slab < numSamples * shape.numChains; ++slab)
        misc_addVectorsInPlace(offsetTest, numTestObservations,
                               out + slab * numTestObservations);
    }
    filled = 1;
  });
  return filled;
}

void dbarts_sampler_setTreeStorage(dbarts_sampler* sampler, int keepTrees,
                                   size_t numSamplesToStore) {
  // the store resize allocates; a capacity too large to hold raises with the
  // store and its draws as they were
  callConvertingExceptions("dbarts_sampler_setTreeStorage", [&]() {
    samplerOf(sampler).setTreeStorage(keepTrees != 0, numSamplesToStore);
  });
}

void dbarts_sampler_printTrees(dbarts_sampler* sampler, size_t forest,
                               const size_t* chainIndices,
                               size_t numChainIndices,
                               const size_t* sampleIndices,
                               size_t numSampleIndices,
                               const size_t* treeIndices,
                               size_t numTreeIndices, int useLiveTrees) {
  bartcore::SamplerBase& engine(samplerOf(sampler));
  bartcore::SamplerShape shape = engine.shape();
  // the engine's printers index forests_[forest] unchecked, by design (fast
  // over safe), so this is the only thing between a caller's index and a read
  // past the last forest
  if (forest >= shape.numForests)
    Rf_error("dbarts_sampler_printTrees: forest index out of range");
  // mirrors bartcore_bridge::getTrees: printing live trees needs no saved
  // store, so the empty-store refusal and the saved-draw range check below
  // apply only when the saved store actually serves this call
  bool useSaved = shape.savedTreeCapacity > 0 && useLiveTrees == 0;
  if (useSaved) refuseEmptyTreeStore(engine, "dbarts_sampler_printTrees");
  for (size_t i = 0; i < numChainIndices; ++i) {
    if (chainIndices[i] >= shape.numChains)
      Rf_error("dbarts_sampler_printTrees: chain number out of range");
  }
  // sample numbers address RECORDED DRAWS, oldest first, as getTrees does
  if (useSaved) {
    for (size_t i = 0; i < numSampleIndices; ++i) {
      if (sampleIndices[i] >= shape.numSavedDraws)
        Rf_error("dbarts_sampler_printTrees: sample number out of range");
    }
  }
  // against the NAMED forest's own count, which a multi-forest sampler states
  // per forest (shape.numTrees is forest 0's)
  for (size_t i = 0; i < numTreeIndices; ++i) {
    if (treeIndices[i] >= engine.numTreesInForest(forest))
      Rf_error("dbarts_sampler_printTrees: tree number out of range");
  }
  engine.printTrees(chainIndices, numChainIndices, sampleIndices,
                    numSampleIndices, treeIndices, numTreeIndices, forest,
                    useLiveTrees != 0);
}

void dbarts_sampler_setNumThreads(dbarts_sampler* sampler,
                                  size_t numThreads) {
  samplerOf(sampler).setNumThreads(numThreads);
}

void dbarts_sampler_setVerbose(dbarts_sampler* sampler, int verbose,
                               size_t printEvery) {
  // the print condition is a modulo by printEvery, so 0 is a division by zero
  // rather than "never print"; refuse it here as the R bridge does, since a
  // flat-C caller has no R layer ahead of it
  if (printEvery == 0)
    Rf_error("dbarts_sampler_setVerbose: printEvery must be at least 1");
  samplerOf(sampler).setVerbose(verbose != 0, printEvery);
}

size_t dbarts_sampler_numObservations(const dbarts_sampler* sampler) {
  return samplerOf(sampler).shape().numObservations;
}

size_t dbarts_sampler_numPredictors(const dbarts_sampler* sampler) {
  return samplerOf(sampler).shape().numPredictors;
}

size_t dbarts_sampler_numTestObservations(const dbarts_sampler* sampler) {
  return samplerOf(sampler).shape().numTestObservations;
}

size_t dbarts_sampler_numChains(const dbarts_sampler* sampler) {
  return samplerOf(sampler).shape().numChains;
}

size_t dbarts_sampler_numTrees(const dbarts_sampler* sampler, size_t forest) {
  const bartcore::SamplerBase& engine(samplerOf(sampler));
  // a size_t probe carries no refusal channel, so an out-of-range forest
  // errors: a 0 tree count is indistinguishable from a legitimate answer
  if (forest >= engine.shape().numForests)
    Rf_error("dbarts_sampler_numTrees: forest index out of range");
  return engine.numTreesInForest(forest);
}

size_t dbarts_sampler_numSavedSamples(const dbarts_sampler* sampler) {
  return samplerOf(sampler).shape().numSavedDraws;
}

int dbarts_sampler_kIsSampled(const dbarts_sampler* sampler) {
  return samplerOf(sampler).shape().kIsSampled ? 1 : 0;
}

int dbarts_sampler_usesDart(const dbarts_sampler* sampler) {
  return samplerOf(sampler).shape().usesDart ? 1 : 0;
}

int dbarts_sampler_family(const dbarts_sampler* sampler) {
  bartcore::SamplerShape shape = samplerOf(sampler).shape();
  // the counts-mutation capability is the multinomial coupling's own
  // fingerprint, and every multinomial handle is an R-built sampler's, since
  // no entry here builds one. Every other
  // family maps one to one off shape().family, and AUTO/STUDENT are never
  // reported (creation resolves AUTO, and a Student-t sampler's family IS
  // gaussian)
  if (shape.supportsCountsMutation) return DBARTS_FAMILY_MULTINOMIAL;
  using RF = bartcore::ResponseFamily;
  switch (shape.family) {
  case RF::gaussian: return DBARTS_FAMILY_GAUSSIAN;
  case RF::probit: return DBARTS_FAMILY_PROBIT;
  case RF::logistic: return DBARTS_FAMILY_LOGISTIC;
  case RF::aft: return DBARTS_FAMILY_AFT;
  case RF::ordinal: return DBARTS_FAMILY_ORDINAL;
  case RF::nbinom: return DBARTS_FAMILY_NBINOM;
  }
  return DBARTS_FAMILY_GAUSSIAN; // unreached: ResponseFamily is exhausted above
}

size_t dbarts_sampler_numFittedValuesPerObservation(
  const dbarts_sampler* sampler) {
  return samplerOf(sampler).shape().numReportedLocations;
}

size_t dbarts_sampler_numVariableCountForests(const dbarts_sampler* sampler) {
  return samplerOf(sampler).shape().numVariableCountForests;
}

// Provider-side binding: each real function's address must
// have exactly the type the single-source DBARTS_C_API_LIST ascribes to it, so
// any drift between the readable prototypes above and the list fails dbarts's
// own compile. Placed inside extern "C" so the list-formed pointer type carries
// the same C language linkage as the resolved &function.
#define DBARTS_BIND_ASSERT(ret, name, params, args) \
  static_assert(std::is_same<decltype(&name), ret (*) params>::value, \
                #name " signature drifted from DBARTS_C_API_LIST");
DBARTS_C_API_LIST(DBARTS_BIND_ASSERT)
#undef DBARTS_BIND_ASSERT

} // extern "C"
