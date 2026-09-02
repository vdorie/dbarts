// Implementation of the flat C API (inst/include/dbarts/dbarts.h) over the
// bartcore engine; entry points are registered with R_RegisterCCallable in
// R_interface.cpp. Argument validation is minimal - consumers are compiled
// packages - except where an engine invariant is at stake (categorical
// category codes, column ranges).

#include <dbarts/dbarts.h>

#include <cmath> // isfinite
#include <cstddef> // size_t
#include <cstdint> // uint64_t
#include <cstring> // memcpy
#include <type_traits> // is_same
#include <vector>

#include <external/Rinternals.h>

#include <misc/linearAlgebra.h> // misc_addVectorsInPlace

#include "R_interface_bartcore_common.hpp"

using std::size_t;
using bartcore_bridge::AugmentationInputs;
using bartcore_bridge::AugmentationLaw;
using bartcore_bridge::BartcoreHolder;
using bartcore_bridge::computeWorkingResponse;
using bartcore_bridge::enforceBinaryWeightPolicy;
using bartcore_bridge::drawAugmentation;
using bartcore_bridge::familyCarriesNoWeights;
using bartcore_bridge::isMultiForest;
using bartcore_bridge::refuseCscReferenceAgainstStore;
using bartcore_bridge::refuseEmptyTreeStore;
using bartcore_bridge::refuseGroupedScaleUpdate;
using bartcore_bridge::refuseMultiForestResponseMutation;
using bartcore_bridge::refuseNonBinaryMask;
using bartcore_bridge::refuseSparseLeafCovariate;
using bartcore_bridge::refuseVarianceForestScaleUpdate;
using bartcore_bridge::responseConduitIsFixed;
using bartcore_bridge::ResponseConduit;
using bartcore_bridge::sigmaIsPinned;
using bartcore_bridge::supportFamily;
using bartcore_bridge::testFitsAreUndefined;
using bartcore_bridge::validateColumnValues;
using bartcore_bridge::validateResponseSupport;
using bartcore_bridge::validateTestContainerAgainstStore;

struct dbarts_sampler_t {
  BartcoreHolder* holder;
  SEXP data; // preserved against collection for the columns it lends
  dbarts_sampler_callback callback = NULL; // per-sweep conditioning hook
  void* callbackData = NULL;
};

namespace {

inline bartcore::SamplerBase& samplerOf(dbarts_sampler* sampler) {
  return *sampler->holder->sampler;
}
inline const bartcore::SamplerBase& samplerOf(const dbarts_sampler* sampler) {
  return *sampler->holder->sampler;
}

// The dense predictor matrix of a retained dbartsData spec (@x), or null when
// it is absent or non-dense (a sparse creation spec): the call-time raw source
// for saved-tree replay and cross-grid state restore. The spec is preserved
// against collection for the sampler's lifetime, so the borrow is valid.
const double* predictorsFromDataExpr(SEXP dataExpr) {
  if (dataExpr == NULL || Rf_isNull(dataExpr)) return NULL;
  SEXP xExpr = Rf_getAttrib(dataExpr, Rf_install("x"));
  if (!Rf_isReal(xExpr)) return NULL;
  // an attribute of the preserved spec cannot be collected, so the PROTECT is
  // redundant to that rooting and is what the PROTECT-balance analyzer reads
  PROTECT(xExpr);
  SEXP dims = Rf_getAttrib(xExpr, R_DimSymbol);
  UNPROTECT(1);
  if (Rf_isNull(dims) || Rf_xlength(dims) != 2) return NULL;
  return REAL(xExpr);
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
// is what the CSC implicit-value rule and the reference refusal both key on (a
// subset mutation names its own columns, so the two indexings differ).
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
  // build - so the codes stay where they lie, and the one entrance that
  // indexes the dense raw as a block, mutation, materializes it instead.
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

// The dense block a MUTATION reads. Every mutation kernel indexes values
// column-major, so a non-dense source is materialized here exactly as the R
// bridge materializes its own - a split-channel view included, whose codes
// widen in the materializer - and a plain dense block is passed straight
// through, never copied. The engine's transaction has no arm for a source it
// cannot index: what it reports is acceptance or rollback, so a non-dense view
// must be resolved to a block before it reaches one.
const double* mutationValues(const TranslatedSource& source) {
  if (source.view.isDenseBlock()) return source.view.denseValues;
  double* block = reinterpret_cast<double*>(
    R_alloc(source.view.numRows * source.view.numColumns > 0
              ? source.view.numRows * source.view.numColumns : 1,
            sizeof(double)));
  bartcore::materializePredictorSource(source.view, source.storeTypes, 0,
                                       source.view.numRows, block);
  return block;
}

// A test NA takes a rule's learned missing direction, and a rule learns one
// only where the training column had NAs (ColumnStore::hasMissing gates the
// draw), so on a complete column it would take one fixed branch at every
// split. The R surface refuses first and names the column; this is the flat
// entrances' backstop, which has only the index.
void refuseTestMissingness(const bartcore::ColumnStore& store,
                           const bartcore::PredictorSource& source,
                           const char* caller) {
  bartcore::PredictorSourceColumns columns(source, store.types.data());
  for (size_t j = 0; j < source.numColumns; ++j) {
    if (store.hasMissing[j]) continue;
    bartcore::PredictorSourceColumnReader column = columns.column(j);
    for (size_t i = 0; i < source.numRows; ++i)
      if (bartcore::isNA(column.at(i)))
        Rf_error("%s: test column %zu has missing values but the training "
                 "column had none, so no rule routes them", caller, j + 1);
  }
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

// dbarts_leaf_model for the engine's tag; the two enumerations are separate on
// purpose, since one is an ABI constant and the other an engine detail.
int leafModelTag(bartcore::LeafModelKind kind) {
  switch (kind) {
  case bartcore::LeafModelKind::monotone: return DBARTS_LEAF_MONOTONE;
  case bartcore::LeafModelKind::linear: return DBARTS_LEAF_LINEAR;
  case bartcore::LeafModelKind::gp: return DBARTS_LEAF_GP;
  case bartcore::LeafModelKind::constant: break;
  }
  return DBARTS_LEAF_CONSTANT;
}

// The enumerator text for an admission refusal, over DBARTS_FAMILY_LIST so
// every dbarts_family value names itself; a value outside the enum falls
// through to NULL, and the caller reports the bare number instead.
const char* familyEnumeratorName(int family) {
  switch (family) {
#define DBARTS_FAMILY_NAME_CASE(name, value) case name: return #name;
  DBARTS_FAMILY_LIST(DBARTS_FAMILY_NAME_CASE)
#undef DBARTS_FAMILY_NAME_CASE
  }
  return NULL;
}

// dbarts_sampler_create's family admission: AUTO plus the six flat-creatable
// families, mapped onto resolveFamily's string vocabulary so nothing below
// this function changes. Refuses the two dbarts_family carries for a sampler
// this entry cannot build (STUDENT, MULTINOMIAL) and anything outside the
// enum, naming both the entry and the family.
const char* creationFamilyName(int family) {
  switch (family) {
  case DBARTS_FAMILY_AUTO: return "";
  case DBARTS_FAMILY_GAUSSIAN: return "gaussian";
  case DBARTS_FAMILY_PROBIT: return "probit";
  case DBARTS_FAMILY_LOGISTIC: return "logistic";
  case DBARTS_FAMILY_AFT: return "aft";
  case DBARTS_FAMILY_ORDINAL: return "ordinal";
  case DBARTS_FAMILY_NBINOM: return "nbinom";
  }
  {
    const char* name = familyEnumeratorName(family);
    if (name != NULL)
      Rf_error("dbarts_sampler_create: family %s is refused (accepts "
               "DBARTS_FAMILY_AUTO, GAUSSIAN, PROBIT, LOGISTIC, AFT, ORDINAL "
               "and NBINOM)", name);
  }
  Rf_error("dbarts_sampler_create: family %d is refused (accepts "
           "DBARTS_FAMILY_AUTO, GAUSSIAN, PROBIT, LOGISTIC, AFT, ORDINAL and "
           "NBINOM)", family);
  return NULL; // unreached: Rf_error longjmps
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
static_assert(offsetof(dbarts_results, tau) == sizeof(size_t) + 6 * sizeof(double*));
static_assert(offsetof(dbarts_results, groupEffects) == sizeof(size_t) + 7 * sizeof(double*));
static_assert(offsetof(dbarts_results, logLikelihood) == sizeof(size_t) + 8 * sizeof(double*));
static_assert(offsetof(dbarts_results, dispersion) == sizeof(size_t) + 9 * sizeof(double*));
static_assert(offsetof(dbarts_results, residualDf) == sizeof(size_t) + 10 * sizeof(double*));
static_assert(sizeof(dbarts_results) == sizeof(size_t) + 11 * sizeof(double*),
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

static_assert(offsetof(dbarts_forest_calibration, structSize) == 0);
static_assert(offsetof(dbarts_forest_calibration, priorScale) == sizeof(size_t) + 0 * sizeof(double*));
static_assert(offsetof(dbarts_forest_calibration, priorSd) == sizeof(size_t) + 1 * sizeof(double*));
static_assert(offsetof(dbarts_forest_calibration, priorMean) == sizeof(size_t) + 2 * sizeof(double*));
static_assert(offsetof(dbarts_forest_calibration, k) == sizeof(size_t) + 3 * sizeof(double*));
static_assert(offsetof(dbarts_forest_calibration, responseScale) == sizeof(size_t) + 4 * sizeof(double*));
static_assert(offsetof(dbarts_forest_calibration, responseShift) == sizeof(size_t) + 5 * sizeof(double*));
static_assert(offsetof(dbarts_forest_calibration, kHasHyperprior) == sizeof(size_t) + 6 * sizeof(double*));
static_assert(offsetof(dbarts_forest_calibration, leafModel) == sizeof(size_t) + 7 * sizeof(double*));
static_assert(offsetof(dbarts_forest_calibration, amplitudePriorVariance) == sizeof(size_t) + 8 * sizeof(double*));
static_assert(offsetof(dbarts_forest_calibration, amplitudePriorScale) == sizeof(size_t) + 9 * sizeof(double*));
static_assert(offsetof(dbarts_forest_calibration, nodeScaleFactor) == sizeof(size_t) + 10 * sizeof(double*));
static_assert(offsetof(dbarts_forest_calibration, nodeScaleDivisor) == sizeof(size_t) + 11 * sizeof(double*));
static_assert(offsetof(dbarts_forest_calibration, basisRowNorm) == sizeof(size_t) + 12 * sizeof(double*));
// every member is a pointer so that an omitting caller's structSize cannot
// land inside tail padding and make DBARTS_HAS_FIELD claim a field it does not
// carry; the size assert is what forces an appending author back to the list
static_assert(sizeof(dbarts_forest_calibration) ==
                sizeof(size_t) + 13 * sizeof(double*),
              "dbarts_forest_calibration layout changed; update these offsets, "
              "and bump DBARTS_C_API_MINOR if a field was appended after 1.0-0");

// Compile-time ABI token, checked against the baked DBARTS_C_API_HASH: FNV-1a
// over the stringized DBARTS_C_API_LIST signatures, then the two ABI enums'
// enumerator lists, then dbarts_sampler_callback's parameter list, then the
// layout the compiler gives the three structs that cross the ABI. Anything it
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
// what keeps a platform out of the token. Every member of all three is
// pointer-width, which the alignment asserts hold a future author to.
#define DBARTS_RESULTS_FIELDS(X) \
  X(structSize) X(sigma) X(train) X(test) X(varcount) X(k) X(varprobs) \
  X(tau) X(groupEffects) X(logLikelihood) X(dispersion) X(residualDf)
#define DBARTS_PREDICTOR_SOURCE_FIELDS(X) \
  X(structSize) X(numRows) X(numColumns) X(denseValues) X(numCscColumns) \
  X(cscColumnPointers) X(cscRowIndices) X(cscValues) X(columnSources) \
  X(columnTypes) X(categoryCounts) X(referenceCodes) X(denseCodes) \
  X(numDenseCodeColumns)
#define DBARTS_FOREST_CALIBRATION_FIELDS(X) \
  X(structSize) X(priorScale) X(priorSd) X(priorMean) X(k) X(responseScale) \
  X(responseShift) X(kHasHyperprior) X(leafModel) X(amplitudePriorVariance) \
  X(amplitudePriorScale) X(nodeScaleFactor) X(nodeScaleDivisor) X(basisRowNorm)

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
#define X(field) DBARTS_ALIGN_ASSERT(dbarts_forest_calibration, field)
DBARTS_FOREST_CALIBRATION_FIELDS(X)
#undef X

#define DBARTS_FOLD_STRUCT(hash, type) \
  dbarts_fnv1aValue(dbarts_fnv1a(hash, #type), sizeof(type) / sizeof(void*))
#define DBARTS_FOLD_FIELD(hash, type, field) \
  dbarts_fnv1aValue(dbarts_fnv1a(hash, #field), \
                    offsetof(type, field) / sizeof(void*))

constexpr std::uint64_t dbarts_foldLayout(std::uint64_t hash) {
  hash = DBARTS_FOLD_STRUCT(hash, dbarts_results);
#define X(field) hash = DBARTS_FOLD_FIELD(hash, dbarts_results, field);
  DBARTS_RESULTS_FIELDS(X)
#undef X
  hash = DBARTS_FOLD_STRUCT(hash, dbarts_predictor_source);
#define X(field) hash = DBARTS_FOLD_FIELD(hash, dbarts_predictor_source, field);
  DBARTS_PREDICTOR_SOURCE_FIELDS(X)
#undef X
  hash = DBARTS_FOLD_STRUCT(hash, dbarts_forest_calibration);
#define X(field) \
  hash = DBARTS_FOLD_FIELD(hash, dbarts_forest_calibration, field);
  DBARTS_FOREST_CALIBRATION_FIELDS(X)
#undef X
  return hash;
}

// "NAME=VALUE;" per enumerator, and the callback's parameter list as the
// header spells it (stringizing an expanded macro takes the usual two steps).
#define DBARTS_ENUMERATOR_TEXT(name, value) #name "=" #value ";"
#define DBARTS_STRINGIZE_(text) #text
#define DBARTS_STRINGIZE(text) DBARTS_STRINGIZE_(text)

// The signature half alone, baked privately so a failed build says WHICH half
// moved: both asserts firing means the entry-point list changed, the combined
// one alone means the ABI moved underneath unchanged signatures (a struct's
// layout, an enumerator, the callback's parameters).
constexpr std::uint64_t dbarts_apiSignatureToken =
  dbarts_fnv1a(DBARTS_C_API_DECLS);

constexpr std::uint64_t dbarts_apiToken() {
  std::uint64_t hash = dbarts_apiSignatureToken;
  hash = dbarts_fnv1a(hash, DBARTS_COLUMN_TYPE_LIST(DBARTS_ENUMERATOR_TEXT));
  hash = dbarts_fnv1a(hash, DBARTS_LEAF_MODEL_LIST(DBARTS_ENUMERATOR_TEXT));
  hash = dbarts_fnv1a(hash, DBARTS_FAMILY_LIST(DBARTS_ENUMERATOR_TEXT));
  hash = dbarts_fnv1a(hash, DBARTS_STRINGIZE(DBARTS_SAMPLER_CALLBACK_PARAMS));
  return dbarts_foldLayout(hash);
}
} // namespace
static_assert(dbarts_apiSignatureToken == 0x0b33edcf638a3cd3ULL,
              "dbarts.h C API signatures moved (the entry-point list, not the "
              "layout fold); re-bake this literal here and DBARTS_C_API_HASH "
              "with it");
static_assert(dbarts_apiToken() == DBARTS_C_API_HASH,
              "dbarts.h C ABI changed - a signature, a struct's layout, an ABI "
              "enumerator, or the callback's parameters; re-bake "
              "DBARTS_C_API_HASH in inst/include/dbarts/dbarts.h (and bump "
              "DBARTS_C_API_MAJOR or DBARTS_C_API_MINOR as the change "
              "warrants)");

extern "C" {

int dbarts_apiMajorVersion(void) { return DBARTS_C_API_MAJOR; }
int dbarts_apiMinorVersion(void) { return DBARTS_C_API_MINOR; }
uint64_t dbarts_apiHash(void) { return DBARTS_C_API_HASH; }

dbarts_sampler* dbarts_sampler_create(SEXP control, SEXP model, SEXP data,
                                      int family) {
  BartcoreHolder* holder = bartcore_bridge::createHolder(
    control, model, data, creationFamilyName(family));
  R_PreserveObject(data);
  return new dbarts_sampler_t{holder, data};
}

void dbarts_sampler_destroy(dbarts_sampler* sampler) {
  if (sampler == NULL) return;
  delete sampler->holder;
  R_ReleaseObject(sampler->data);
  delete sampler;
}

void dbarts_sampler_run(dbarts_sampler* sampler, size_t numBurnIn,
                        size_t numSamples, dbarts_results* results) {
  bartcore::SamplerShape shape = samplerOf(sampler).shape();
  bartcore::Results engineResults;
  // the internal location stride (invisible to the frozen dbarts_results ABI):
  // 1 for every dbarts.h-created sampler, since the flat C API builds no
  // multi-location model, so the caller's n x numSamples train/test hold
  //
  // numVariableCountForests is deliberately NOT set: dbarts_results declares no
  // forest count, so the field stays at its default 1 and the engine writes the
  // single numPredictors slab per sample this struct documents - the reported
  // (prognostic) forest - even on the BCF samplers this entry point can create
  engineResults.numReportedLocations = shape.numReportedLocations;
  // A zero structSize means the caller forgot to set it (see DBARTS_RESULTS_INIT):
  // reject loudly instead of silently skipping every field and handing back an
  // uninitialized buffer - the flat-API footgun that fed garbage draws to a
  // consumer's Gibbs loop. A nonzero older/smaller structSize stays valid.
  if (results != NULL && results->structSize == 0)
    Rf_error("dbarts_sampler_run: results.structSize is 0 - set it to "
             "sizeof(dbarts_results) (e.g. dbarts_results r = DBARTS_RESULTS_INIT)");

  if (results != NULL && numSamples > 0) {
    // A field is filled only when present-by-size AND non-null. offsetof is
    // against the library's (newest) layout; fields only append, so it
    // equals the caller's offset and structSize bounds the buffer.
#define FILL(field, member) \
  engineResults.member = DBARTS_RESULTS_HAS(results, field) ? results->field : NULL
    FILL(sigma, sigma);
    FILL(train, trainingFits);
    FILL(test, testFits);
    FILL(varcount, variableCounts);
    FILL(k, k);
    FILL(varprobs, splitProbabilities);
    FILL(tau, tau);
    FILL(groupEffects, groupEffects);
    FILL(logLikelihood, logLikelihood);
    FILL(dispersion, dispersion);
    FILL(residualDf, residualDf);
#undef FILL
  }

  bartcore::SweepCallback onSweep;
  if (sampler->callback != NULL) {
    if (shape.numThreads > 1 && shape.numChains > 1)
      Rf_error("dbarts_sampler_run: a per-sweep callback cannot run while "
               "chains execute on worker threads");
    dbarts_sampler_callback fn = sampler->callback;
    void* userData = sampler->callbackData;
    onSweep = [sampler, fn, userData](size_t chainIndex, size_t sweepIndex,
                                      bool isBurnIn) -> bool {
      return fn(userData, sampler, chainIndex, sweepIndex,
                isBurnIn ? 1 : 0) == 0;
    };
  }

  // The engine samples only from each chain's own Mersenne Twister (seeded
  // from R's stream once at creation), never from R's stream during a run, so
  // no GetRNGstate/PutRNGstate bracket is needed here - and none is left
  // unbalanced by a longjmp out of the engine.
  samplerOf(sampler).run(numBurnIn, numSamples, engineResults, {}, onSweep);
}

void dbarts_sampler_setCallback(dbarts_sampler* sampler,
                                dbarts_sampler_callback callback,
                                void* userData) {
  if (callback != NULL) {
    bartcore::SamplerShape shape = samplerOf(sampler).shape();
    if (shape.numThreads > 1 && shape.numChains > 1)
      Rf_error("dbarts_sampler_setCallback: a per-sweep callback requires "
               "chains to run inline (numThreads == 1 or numChains == 1)");
  }
  sampler->callback = callback;
  sampler->callbackData = userData;
}

void dbarts_sampler_sampleTreesFromPrior(dbarts_sampler* sampler) {
  // draws from the chain RNG only, not R's stream (see dbarts_sampler_run)
  samplerOf(sampler).sampleTreesFromPrior();
}

void dbarts_sampler_sampleNodeParametersFromPrior(dbarts_sampler* sampler) {
  // draws from the chain RNG only, not R's stream (see dbarts_sampler_run)
  samplerOf(sampler).sampleNodeParametersFromPrior();
}

int dbarts_sampler_setResponse(dbarts_sampler* sampler, const double* y,
                               int updateScale) {
  // the capability answer, which no argument would have changed: a coupling
  // that fixes its response conduit at creation
  if (responseConduitIsFixed(samplerOf(sampler).shape())) return 0;
  // the shared conduit guard, not the whole-data refusal: a two-forest sampler
  // is flat-creatable, and its response swap is opt-in and scale-pinned rather
  // than refused - the same rule bartcore_setResponse applies
  refuseMultiForestResponseMutation(samplerOf(sampler),
                                    "dbarts_sampler_setResponse",
                                    ResponseConduit::response, updateScale);
  refuseVarianceForestScaleUpdate(samplerOf(sampler),
                                  "dbarts_sampler_setResponse",
                                  ResponseConduit::response, updateScale);
  // reachable here: this entry's control carries whatever bartcore.groups
  // attribute the consumer put on it, so a flat-API sampler can be grouped
  refuseGroupedScaleUpdate(samplerOf(sampler), "dbarts_sampler_setResponse",
                           ResponseConduit::response, updateScale);
  // the one place minimal validation is not enough: an out-of-support y is a
  // silently garbage latent draw for probit/ordinal and, for nbinom, an
  // uncatchable crash inside the count histogram (see validateResponseSupport)
  bartcore::SamplerShape shape = samplerOf(sampler).shape();
  validateResponseSupport(shape.family, shape.numOrdinalThresholds + 1, y,
                          shape.numObservations, "dbarts_sampler_setResponse");
  // the probit latent redraw draws from the chain RNG, not R's stream
  samplerOf(sampler).setResponse(y, updateScale != 0);
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
  refuseVarianceForestScaleUpdate(samplerOf(sampler),
                                  "dbarts_sampler_setOffset",
                                  ResponseConduit::offset, updateScale);
  refuseGroupedScaleUpdate(samplerOf(sampler), "dbarts_sampler_setOffset",
                           ResponseConduit::offset, updateScale);
  samplerOf(sampler).setOffset(offset, updateScale != 0);
  return 1;
}

int dbarts_sampler_setWeights(dbarts_sampler* sampler,
                              const double* weights) {
  // both capability answers: a coupling that fixes the weight conduit at
  // creation (the weight conduit has no scale to pin, so that is its whole
  // condition), and the family rule bartcore_setWeights states, which this
  // entry used to drop on the floor - probit/ordinal/aft/nbinom would install
  // a vector nothing reads
  bartcore::SamplerShape weightShape = samplerOf(sampler).shape();
  if (responseConduitIsFixed(weightShape) ||
      familyCarriesNoWeights(samplerOf(sampler)))
    return 0;
  // the value half stays a raise: a logistic count that is not a positive
  // integer leaves a row carrying a precision no observation of it justifies,
  // a gaussian weight that is negative or not finite corrupts the leaf
  // sufficient statistics outright, and a different vector would have worked
  enforceBinaryWeightPolicy(weightShape.family, weights,
                            weightShape.numObservations);
  samplerOf(sampler).setWeights(weights);
  return 1;
}

int dbarts_sampler_setSigma(dbarts_sampler* sampler, double sigma) {
  // reachable here: this entry creates probit/logistic/ordinal/nbinom samplers
  // by family name, and dbartsSpec(variance = ) hands a consumer a
  // heteroscedastic control (multinomial has no flat creation path)
  if (sigmaIsPinned(samplerOf(sampler))) return 0;
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

int dbarts_sampler_getDispersion(const dbarts_sampler* sampler, double* out) {
  const bartcore::SamplerBase& engine(samplerOf(sampler));
  bartcore::SamplerShape shape = engine.shape();
  // the capability answer is the return value, as it is for the R bridge's
  // NULL: the value itself would mean nothing off a family carrying one
  if (!shape.carriesDispersion) return 0;
  for (size_t c = 0; c < shape.numChains; ++c) out[c] = engine.dispersion(c);
  return 1;
}

int dbarts_sampler_setPredictor(dbarts_sampler* sampler,
                                const dbarts_predictor_source* x,
                                int forceUpdate, int updateCutPoints) {
  // Unguarded again, as it was before the stop-loss, but for the opposite
  // reason: the two-phase transaction now covers every forest and the variance
  // forest, so an unforced call here vetoes or rolls back rather than
  // misrouting. Mirrors the R bridge's bartcore_setPredictor, which carries no
  // guard either.
  bartcore::SamplerBase& engine(samplerOf(sampler));
  bartcore::SamplerShape shape = engine.shape();
  size_t numObservations = shape.numObservations;
  void* scratch = vmaxget();
  TranslatedSource source =
    translateSource(engine.data(), x, NULL, shape.numPredictors,
                    numObservations, "dbarts_sampler_setPredictor");
  const double* values = mutationValues(source);
  for (size_t j = 0; j < shape.numPredictors; ++j)
    validateColumnValues(engine.data(), j, values + j * numObservations,
                         numObservations);

  bartcore::PredictorUpdateResult result =
    engine.setPredictor(values, forceUpdate != 0, updateCutPoints != 0);
  vmaxset(scratch);
  if (result == bartcore::PredictorUpdateResult::invalidCutPoints)
    Rf_error("number of induced cut points in new predictor less than "
             "previous: old splits would be invalid");
  return result == bartcore::PredictorUpdateResult::accepted ? 1 : 0;
}

int dbarts_sampler_updatePredictor(dbarts_sampler* sampler,
                                   const dbarts_predictor_source* x,
                                   const size_t* columns, size_t numColumns,
                                   int forceUpdate, int updateCutPoints) {
  // unguarded, as dbarts_sampler_setPredictor above
  bartcore::SamplerBase& engine(samplerOf(sampler));
  bartcore::SamplerShape shape = engine.shape();
  size_t numObservations = shape.numObservations;
  for (size_t k = 0; k < numColumns; ++k)
    if (columns[k] >= shape.numPredictors)
      Rf_error("dbarts_sampler_updatePredictor: column out of range");

  void* scratch = vmaxget();
  // the source's columns are in ARGUMENT order, so column k of the source is
  // store column columns[k] - what the type gather and the validation below
  // both index by
  TranslatedSource source =
    translateSource(engine.data(), x, columns, numColumns, numObservations,
                    "dbarts_sampler_updatePredictor");
  const double* values = mutationValues(source);
  for (size_t k = 0; k < numColumns; ++k)
    validateColumnValues(engine.data(), columns[k],
                         values + k * numObservations, numObservations);

  bartcore::PredictorUpdateResult result = engine.updatePredictor(
    values, columns, numColumns, forceUpdate != 0, updateCutPoints != 0);
  vmaxset(scratch);
  if (result == bartcore::PredictorUpdateResult::invalidCutPoints)
    Rf_error("number of induced cut points in new predictor less than "
             "previous: old splits would be invalid");
  return result == bartcore::PredictorUpdateResult::accepted ? 1 : 0;
}

int dbarts_sampler_setTestPredictors(dbarts_sampler* sampler,
                                     const dbarts_predictor_source* xTest) {
  bartcore::SamplerBase& engine(samplerOf(sampler));
  // an amplitude coupling's test blend is undefined without an off-sample
  // basis, and its whole test surface is refused; a multi-forest model whose
  // blend IS defined passes (the predicate is the blend, not the forest count)
  if (testFitsAreUndefined(engine)) return 0;
  if (xTest == NULL) {
    // removal is the whole no-test-data state, the offset included: the engine
    // preserves a test offset across a test-store REBUILD (the caller keeps the
    // two lengths consistent) but a removal leaves it describing rows that no
    // longer exist, and the next install would silently re-adopt it
    engine.setTestPredictors(NULL, 0);
    engine.setTestOffset(NULL);
    return 1;
  }
  void* scratch = vmaxget();
  TranslatedSource source =
    translateSource(engine.data(), xTest, NULL, engine.shape().numPredictors, 0,
                    "dbarts_sampler_setTestPredictors");
  validateTestSource(engine, source, "dbarts_sampler_setTestPredictors");
  // and a rebuild at a different row count would read the caller's offset past
  // its end on every recorded test fit, since nothing downstream re-checks the
  // two against each other; the pair moves together
  if (engine.data().testOffset != NULL &&
      source.view.numRows != engine.shape().numTestObservations)
    Rf_error("dbarts_sampler_setTestPredictors: test offset length would no "
             "longer match; set the predictors and offset together");
  // the store build answers the leaf-covariate refusal with a false return;
  // defense in depth, since validateTestSource has already raised it - a
  // discarded false would leave the store holding its PREVIOUS rows and report
  // them as the new test set. The build's own level-code refusal shares the
  // return and not this text, and cannot reach it: validateTestSource bounds
  // every code against the store first
  if (!engine.setTestData(source.view))
    Rf_error("a leaf covariate column cannot be a sparse test column; "
             "supply it as a dense test column");
  vmaxset(scratch);
  return 1;
}

int dbarts_sampler_setTestOffset(dbarts_sampler* sampler,
                                 const double* offsetTest) {
  // a multi-forest test offset lands after the forests are blended, which is a
  // fixed property of the sampler; see dbarts_sampler_setResponse
  if (isMultiForest(samplerOf(sampler))) return 0;
  samplerOf(sampler).setTestOffset(offsetTest);
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
  refuseEmptyTreeStore(engine, "dbarts_sampler_predict");
  bartcore::SamplerShape shape = engine.shape();
  void* scratch = vmaxget();
  TranslatedSource source = translateSource(
    engine.data(), xTest, NULL, shape.numPredictors, 0,
    "dbarts_sampler_predict");
  // a read-only replay builds no store, so the leaf-covariate rule is checked
  // on the view itself rather than answered by a store build
  validateTestSource(engine, source, "dbarts_sampler_predict");
  size_t numTestObservations = source.view.numRows;

  engine.predict(source.view, numTestObservations, NULL, numThreads, out);
  vmaxset(scratch);

  if (offsetTest != NULL) {
    size_t capacity = shape.savedTreeCapacity;
    size_t numSamples = capacity > 0 ? shape.numSavedDraws : 1;
    for (size_t slab = 0; slab < numSamples * shape.numChains; ++slab)
      misc_addVectorsInPlace(offsetTest, numTestObservations,
                             out + slab * numTestObservations);
  }
  return 1;
}

void dbarts_sampler_setTreeStorage(dbarts_sampler* sampler, int keepTrees,
                                   size_t numSamplesToStore) {
  samplerOf(sampler).setTreeStorage(keepTrees != 0, numSamplesToStore);
}

SEXP dbarts_sampler_getTrees(dbarts_sampler* sampler, size_t forest,
                             const size_t* chainIndices,
                             size_t numChainIndices,
                             const size_t* sampleIndices,
                             size_t numSampleIndices,
                             const size_t* treeIndices, size_t numTreeIndices,
                             int useLiveTrees) {
  if (forest >= samplerOf(sampler).shape().numForests)
    Rf_error("dbarts_sampler_getTrees: forest index out of range");
  // the n column replays the retained creation spec's predictors through each
  // saved tree; the engine keeps no matrix, and a caller that mutated
  // predictors since creation sees the pre-mutation spec
  const double* replay = predictorsFromDataExpr(sampler->data);
  return bartcore_bridge::getTrees(
    samplerOf(sampler), chainIndices, numChainIndices, sampleIndices,
    numSampleIndices, treeIndices, numTreeIndices, useLiveTrees != 0, NULL, 0,
    replay, samplerOf(sampler).shape().numObservations, forest,
    "dbarts_sampler_getTrees");
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

SEXP dbarts_sampler_storeState(dbarts_sampler* sampler) {
  return bartcore_bridge::storeState(samplerOf(sampler));
}

void dbarts_sampler_setState(dbarts_sampler* sampler, SEXP state) {
  // a cross-grid restore re-quantizes dense columns from the retained creation
  // spec; a same-spec continuation (the contract) skips per column and reads no
  // raw. A flat-C caller that mutated predictors since sees the creation spec.
  bartcore_bridge::setState(samplerOf(sampler), state,
                            predictorsFromDataExpr(sampler->data));
}

void dbarts_sampler_setNumThreads(dbarts_sampler* sampler,
                                  size_t numThreads) {
  samplerOf(sampler).setNumThreads(numThreads);
}

void dbarts_sampler_setNumThin(dbarts_sampler* sampler, size_t numThin) {
  samplerOf(sampler).setNumThin(numThin);
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
  // fingerprint; DBARTS_FAMILY_MULTINOMIAL is reserved for multinomial
  // creation opening, though no entry here builds one yet. Every other
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

size_t dbarts_sampler_numForests(const dbarts_sampler* sampler) {
  return samplerOf(sampler).shape().numForests;
}

int dbarts_sampler_setForestBasis(dbarts_sampler* sampler, size_t forest,
                                  const double* basisRowMajor,
                                  size_t numColumns) {
  bartcore::SamplerBase& engine(samplerOf(sampler));
  bartcore::SamplerShape shape = engine.shape();
  // a capability probe rather than a forest count, matching the bridge: a
  // basis is defined only as what the amplitudes multiply, and a K-forest
  // multinomial would defeat a numForests test
  if (engine.totalAmplitudes() == 0) return 0;
  if (forest >= shape.numForests) return 0;
  if (basisRowMajor == NULL)
    Rf_error("dbarts_sampler_setForestBasis: 'basisRowMajor' cannot be NULL");
  if (numColumns == 0)
    Rf_error("dbarts_sampler_setForestBasis: a basis needs at least one "
             "column");
  // ROW-major, row i at basisRowMajor + i * numColumns, as the parameter name
  // states and the engine's own contraction reads it; the values are copied
  // through
  size_t numValues = shape.numObservations * numColumns;
  for (size_t i = 0; i < numValues; ++i)
    if (!std::isfinite(basisRowMajor[i]))
      Rf_error("dbarts_sampler_setForestBasis: a basis value is not finite");
  return engine.setForestBasis(forest, basisRowMajor, numColumns) ? 1 : 0;
}

int dbarts_sampler_getForestFits(const dbarts_sampler* sampler, size_t forest,
                                 double* out) {
  const bartcore::SamplerBase& engine(samplerOf(sampler));
  bartcore::SamplerShape shape = engine.shape();
  if (forest >= shape.numForests) return 0;
  for (size_t c = 0; c < shape.numChains; ++c)
    engine.forestTotalFits(c, forest, out + c * shape.numObservations);
  return 1;
}

size_t dbarts_sampler_numForestAmplitudes(const dbarts_sampler* sampler,
                                          size_t forest) {
  const bartcore::SamplerBase& engine(samplerOf(sampler));
  // a size_t probe carries no refusal channel; see dbarts_sampler_numTrees
  if (forest >= engine.shape().numForests)
    Rf_error("dbarts_sampler_numForestAmplitudes: forest index out of range");
  // the amplitude vector is ragged by construction: as wide as the forest's
  // own basis, which every forest carries independently
  return engine.numForestAmplitudes(forest);
}

int dbarts_sampler_getForestAmplitudes(const dbarts_sampler* sampler,
                                       size_t forest, double* out) {
  const bartcore::SamplerBase& engine(samplerOf(sampler));
  bartcore::SamplerShape shape = engine.shape();
  if (forest >= shape.numForests) return 0;
  // the amplitudes are a model property, so the length answers whether there
  // are any at all; a coupling carrying none writes nothing
  size_t total = engine.totalAmplitudes();
  if (total == 0) return 0;
  size_t offset = 0;
  for (size_t f = 0; f < forest; ++f) offset += engine.numForestAmplitudes(f);
  size_t numAmplitudes = engine.numForestAmplitudes(forest);
  std::vector<double> amplitudes(total);
  for (size_t c = 0; c < shape.numChains; ++c) {
    engine.amplitudes(c, amplitudes.data());
    for (size_t j = 0; j < numAmplitudes; ++j)
      out[c * numAmplitudes + j] = amplitudes[offset + j];
  }
  return 1;
}

int dbarts_sampler_setForestWeights(dbarts_sampler* sampler, size_t forest,
                                    const double* weights) {
  bartcore::SamplerBase& engine(samplerOf(sampler));
  bartcore::SamplerShape shape = engine.shape();
  // the capability probe comes FIRST, and it is not a forest count: a K-forest
  // multinomial carries several forests and admits no such weight
  if (!shape.supportsForestWeights) return 0;
  if (forest >= shape.numForests) return 0;
  if (weights != NULL)
    for (size_t i = 0; i < shape.numObservations; ++i)
      if (!R_FINITE(weights[i]) || weights[i] < 0.0)
        Rf_error("dbarts_sampler_setForestWeights: weights must be finite and "
                 "non-negative");
  // BORROWED, unlike the bridge's copy into a holder-owned buffer: the chains
  // hold this pointer until it is replaced, so the caller owns the array for
  // the sampler's life
  return engine.setForestWeights(forest, weights) ? 1 : 0;
}

int dbarts_sampler_getForestCalibration(const dbarts_sampler* sampler,
                                        size_t forest,
                                        dbarts_forest_calibration* out) {
  if (out == NULL || out->structSize == 0)
    Rf_error("dbarts_sampler_getForestCalibration: out.structSize is 0 - set it "
             "to sizeof(dbarts_forest_calibration) (e.g. "
             "dbarts_forest_calibration c = DBARTS_FOREST_CALIBRATION_INIT)");
  const bartcore::SamplerBase& engine(samplerOf(sampler));
  bartcore::SamplerShape shape = engine.shape();
  if (forest >= shape.numForests) return 0;
  int leafModel = leafModelTag(shape.leafModel);
  for (size_t c = 0; c < shape.numChains; ++c) {
    bartcore::ForestCalibration calibration =
      engine.forestCalibration(c, forest);
    // a member is filled only when both present-by-size and non-null, the
    // dbarts_results contract read in the same direction
#define FILL(field, value) \
  if (DBARTS_HAS_FIELD(dbarts_forest_calibration, out, field) && \
      out->field != NULL) \
    out->field[c] = (value)
    FILL(priorScale, calibration.priorScale);
    FILL(priorSd, calibration.priorSd);
    FILL(priorMean, calibration.priorMean);
    FILL(k, calibration.k);
    FILL(responseScale, calibration.responseScale);
    FILL(responseShift, calibration.responseShift);
    FILL(kHasHyperprior, calibration.kHasHyperprior ? 1 : 0);
    FILL(leafModel, leafModel);
    // the calibration map's five, appended below the 1.0-0 boundary: a
    // pre-append caller's structSize stops here and its buffers are untouched
    FILL(amplitudePriorVariance, calibration.amplitudePriorVariance);
    FILL(amplitudePriorScale, calibration.amplitudePriorScale);
    FILL(nodeScaleFactor, calibration.nodeScaleFactor);
    FILL(nodeScaleDivisor, calibration.nodeScaleDivisor);
    FILL(basisRowNorm, calibration.basisRowNorm);
#undef FILL
  }
  return 1;
}

int dbarts_sampler_setForestPriorScale(dbarts_sampler* sampler, size_t forest,
                                       double priorScale) {
  bartcore::SamplerBase& engine(samplerOf(sampler));
  // two channels, as the getter has: a capability answer returns, a malformed
  // value raises
  if (forest >= engine.shape().numForests) return 0;
  if (!R_FINITE(priorScale) || priorScale <= 0.0)
    Rf_error("dbarts_sampler_setForestPriorScale: 'priorScale' must be a "
             "positive finite number");
  return engine.setForestPriorScale(forest, priorScale) ? 1 : 0;
}

int dbarts_sampler_setActiveRows(dbarts_sampler* sampler,
                                 const double* active) {
  bartcore::SamplerBase& engine(samplerOf(sampler));
  // the capability probe comes FIRST and never switches on the family -
  // Student-t reports as gaussian - and it is the only thing this return value
  // reports; the all-ones normalization and the copy are the engine's
  if (!engine.shape().supportsActiveRows) return 0;
  // the value refusal is the other channel, shared with the R bridge: a
  // fractional element is recoverable, so it raises rather than answering the
  // capability question with a no. The engine scans again on its own contract,
  // which leaves the ? 1 : 0 below as defense in depth
  refuseNonBinaryMask(active, engine.shape().numObservations);
  return engine.setActiveRows(active) ? 1 : 0;
}

// dbarts_drawLatents/dbarts_workingResponse's family admission, and the one
// place a header family token becomes an augmentation law. The map is
// deliberately not one-to-one in either direction: DBARTS_FAMILY_STUDENT names
// no sampler family (a Student-t sampler's family IS gaussian) yet is the only
// token carrying the scale-mixture law, while AUTO, GAUSSIAN and MULTINOMIAL
// name families that carry no augmentation law at all and are refused here,
// as is anything outside the enum. Every refusal names both the entry and the
// family.
static AugmentationLaw augmentationLawOf(int family, const char* caller) {
  switch (family) {
  case DBARTS_FAMILY_PROBIT: return AugmentationLaw::probit;
  case DBARTS_FAMILY_LOGISTIC: return AugmentationLaw::logistic;
  case DBARTS_FAMILY_ORDINAL: return AugmentationLaw::ordinal;
  case DBARTS_FAMILY_AFT: return AugmentationLaw::aft;
  case DBARTS_FAMILY_NBINOM: return AugmentationLaw::nbinom;
  case DBARTS_FAMILY_STUDENT: return AugmentationLaw::studentT;
  }
  {
    const char* name = familyEnumeratorName(family);
    if (name != NULL)
      Rf_error("%s: family %s is refused (accepts DBARTS_FAMILY_PROBIT, "
               "LOGISTIC, ORDINAL, AFT, NBINOM and STUDENT)", caller, name);
  }
  Rf_error("%s: family %d is refused (accepts DBARTS_FAMILY_PROBIT, "
           "LOGISTIC, ORDINAL, AFT, NBINOM and STUDENT)", caller, family);
  return AugmentationLaw::probit; // unreached: Rf_error longjmps
}

// This file reads the laws through if-chains rather than a switch - the token
// map above and the per-law argument rules below - so -Wswitch cannot see a law
// either one omits; this assertion is the only tripwire. Update it only
// together with both.
static_assert(bartcore_bridge::numAugmentationLaws == 6,
              "a law was added to AugmentationLaw: augmentationLawOf's token "
              "map and the argument rules below each need an arm for it");

// Applies the rules R/augmentation.R applies ahead of the R helpers, which the
// wrapped forms below have no R layer to inherit. Only the half a wrong
// argument cannot survive: the parameter each law REQUIRES, there being no
// default to fall back on (the ordinal arm indexes its cut points
// unconditionally, and the Polya-Gamma working response divides by its
// latent), and the logistic counts. A parameter no law reads is IGNORED rather
// than refused by name as R refuses it, a C caller having no way to leave one
// out. drawing selects the draw's own scalars.
static void augmentationArguments(AugmentationLaw law,
                                  const AugmentationInputs& in, bool drawing,
                                  const char* caller) {
  using AL = AugmentationLaw;
  if (in.fit == NULL || in.y == NULL)
    Rf_error("%s: the %s and response vectors are required", caller,
             drawing ? "fit" : "latent");
  if (drawing && law == AL::ordinal &&
      (in.ordinalThresholds == NULL || in.numOrdinalThresholds == 0))
    Rf_error("%s: family \"ordinal\" requires cut points", caller);
  if (drawing && (law == AL::aft || law == AL::studentT) &&
      !(R_FINITE(in.sigma) && in.sigma > 0.0))
    Rf_error("%s: the aft and student laws require a positive 'sigma'", caller);
  if (drawing && law == AL::studentT && !(R_FINITE(in.df) && in.df > 0.0))
    Rf_error("%s: family \"student\" requires positive 'df'", caller);
  if (law == AL::nbinom &&
      !(R_FINITE(in.dispersion) && in.dispersion > 0.0))
    Rf_error("%s: family \"nbinom\" requires a positive 'dispersion'", caller);

  bool precisionLatent =
    !drawing && (law == AL::logistic || law == AL::nbinom);
  // the counts rule belongs to the logistic law alone, which is the only one
  // that reads a weight (as its Polya-Gamma copy count); a law that does not
  // read them ignores them, as this entry's contract states
  bool readsWeights = in.weights != NULL && law == AL::logistic;
  if (!readsWeights && !precisionLatent) return;
  for (size_t i = 0; i < in.numObservations; ++i) {
    if (readsWeights && (!(in.weights[i] > 0.0) ||
                         in.weights[i] != std::floor(in.weights[i])))
      Rf_error("%s: 'weights' are observation counts: positive whole numbers",
               caller);
    if (precisionLatent && !(in.fit[i] > 0.0))
      Rf_error("%s: 'latent' is the precision the working response divides by: "
               "it must be positive", caller);
  }
}

void dbarts_drawLatents(int family, size_t numObservations,
                        const double* fit, const double* y,
                        const double* weights, const double* offset,
                        double sigma, double dispersion,
                        const double* ordinalThresholds,
                        size_t numOrdinalThresholds, double df, double* out) {
  AugmentationLaw law = augmentationLawOf(family, "dbarts_drawLatents");
  AugmentationInputs in{.numObservations = numObservations, .fit = fit, .y = y,
                        .weights = weights, .offset = offset,
                        .ordinalThresholds = ordinalThresholds,
                        .numOrdinalThresholds = numOrdinalThresholds,
                        .sigma = sigma, .dispersion = dispersion, .df = df};
  augmentationArguments(law, in, true, "dbarts_drawLatents");
  if (out == NULL) Rf_error("dbarts_drawLatents: 'out' cannot be NULL");
  // the same support rule every conduit that swaps a y states
  validateResponseSupport(supportFamily(law), in.numOrdinalThresholds + 1, in.y,
                          in.numObservations, "dbarts_drawLatents");
  drawAugmentation(law, in, out, "dbarts_drawLatents");
}

void dbarts_workingResponse(int family, size_t numObservations,
                            const double* latent, const double* y,
                            const double* weights, const double* offset,
                            double dispersion, double* out) {
  AugmentationLaw law = augmentationLawOf(family, "dbarts_workingResponse");
  AugmentationInputs in{.numObservations = numObservations, .fit = latent,
                        .y = y, .weights = weights, .offset = offset,
                        .dispersion = dispersion};
  augmentationArguments(law, in, false, "dbarts_workingResponse");
  if (out == NULL) Rf_error("dbarts_workingResponse: 'out' cannot be NULL");
  // the ordinal working response is the latent less the offset, so y never
  // enters it and there is no category count here to state its support against
  if (law != AugmentationLaw::ordinal)
    validateResponseSupport(supportFamily(law), 0, in.y, in.numObservations,
                            "dbarts_workingResponse");
  computeWorkingResponse(law, in, latent, out);
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
