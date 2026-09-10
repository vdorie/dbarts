#ifndef DBARTS_DBARTS_H
#define DBARTS_DBARTS_H

/// \file dbarts.h
/// Flat C interface to the dbarts sampler, for packages that drive BART as
/// a conditional model inside a larger sampler (LinkingTo: dbarts).
///
/// LinkingTo compiles a consumer against this header, but does not load
/// dbarts at runtime: R loads a package's declared dependencies from its
/// NAMESPACE file's import directives, not from DESCRIPTION's Imports field,
/// and R_GetCCallable (what the stubs below call on first use) can only find
/// an entry point in a package whose shared object is already loaded. A
/// consumer must therefore also carry an importFrom(dbarts, ...) or
/// import(dbarts) directive in its own NAMESPACE - Imports: dbarts alone
/// compiles clean and then fails at load time with "function
/// 'dbarts_apiHash' not provided by package 'dbarts'" (or any other stub, on
/// first call), naming neither the missing directive nor the cure.
///
/// Every function is registered with R_RegisterCCallable under its own name.
/// A consumer has two ways to reach an entry point. It can look each one up by
/// hand with R_GetCCallable("dbarts", "<name>") and cast the DL_FUNC to the
/// matching signature, or - the supported path - it can define DBARTS_USE_STUBS
/// before including this header, which replaces the prototypes below with
/// same-name static inline stubs that resolve and cache the pointer on first
/// call. The stubs are generated from DBARTS_C_API_LIST, the single source of
/// truth for the surface, so a rebuild always re-derives the consumer's call
/// types from this header and a stale hand-rolled signature cannot drift.
///
/// The version is two components: check dbarts_apiMajorVersion() ==
/// DBARTS_C_API_MAJOR && dbarts_apiMinorVersion() >= DBARTS_C_API_MINOR at load
/// time (major = incompatible change, minor = additive). The DBARTS_USE_STUBS
/// stubs enforce exactly this handshake on first use of any entry point.
/// dbarts_apiHash() adds an exact ABI token beside it, for a consumer that
/// wants lockstep rather than compatibility; the stubs check it too, but only
/// when the consumer defines DBARTS_REQUIRE_EXACT_ABI before including this
/// header. The interface only ever grows: names and
/// signatures below are stable and function additions arrive under new names (a
/// minor bump), while dbarts_results grows in place by appending fields - its
/// leading structSize keeps callers compiled against an older layout safe, and
/// the caller-filled struct below (dbarts_predictor_source) carries the same
/// leading member for the same reason, in the other direction.
///
/// Contracts common to all entry points:
/// - Validation is deliberately partial: consumers are compiled packages, so
///   what an entry checks is what the engine's invariants or the caller's own
///   buffers depend on - a struct's structSize, a source's declared shape, a
///   response value outside its family's support, a categorical code the
///   sampler does not hold, a capability the model does not carry. What is NOT
///   checked is the plain pointer: the sampler handle, an output buffer, and a
///   required input vector are dereferenced as handed, so a null (or
///   destroyed, or short) one crashes rather than raising. The one exception
///   is dbarts_sampler_destroy itself, which is idempotent.
/// - A non-void return is one of three things, and each entry's own doc says
///   which. A VALUE: the number IS the answer and carries no refusal, which is
///   what the counts, dbarts_sampler_kIsSampled, dbarts_sampler_usesDart,
///   dbarts_sampler_family and the version accessors report - an int here is
///   not a status. Or a CAPABILITY STATUS, whose rule
///   is that 0 means the SAMPLER cannot do this at all - no argument would
///   have worked - and nothing was touched, so a host driving a sampler it did
///   not build probes the channel instead of unwinding through its own frames.
///   An Rf_error means THIS call is wrong - a value outside a family's
///   support, a source whose declared shape disagrees, a malformed struct, a
///   scale update against a transform the sampler pins - and longjmps through
///   the caller's frames, so call from contexts that are safe to unwind. A
///   DISCARDED capability 0 leaves the sampler unchanged and the run
///   conditioned on what it held before; that answer is a fixed property of
///   the sampler, so test it once at setup rather than every sweep.
/// - Where an entry names a forest, forest is the argument after the sampler:
///   a 0-based index over the sampler's own forests, qualifying every argument
///   after it - a tree index list is read against THAT forest's tree count.
///   0 is the only index a single-forest sampler has, and a sampler this
///   header can drive but not describe (a multi-forest one, built from R)
///   states its count nowhere here. An index past the last forest RAISES on
///   both entries that take one: dbarts_sampler_numTrees, whose size_t value
///   carries no refusal a caller could tell from a legitimate answer, and
///   dbarts_sampler_printTrees, which carries no status channel at all.
/// - The functions that draw (dbarts_sampler_run,
///   dbarts_sampler_sampleTreesFromPrior) manage R's RNG state internally and
///   must be called from the main R thread. Do not wrap them in a
///   GetRNGstate/PutRNGstate bracket that spans your own draws through R's
///   API. dbarts_sampler_predict is main-R-thread-only for a separate reason:
///   it is R_alloc-backed internally, and R_alloc is unsafe off that thread.
/// - THE HANDLE. There is no creation entry here: a dbarts_sampler* is the
///   address stored in an R dbartsSampler object's external pointer, which a
///   consumer builds from R and reads with R_ExternalPtrAddr in its own
///   code. dbarts() takes a formula or a matrix, not a model object, so a
///   consumer starting from its own already-built (control, model, data)
///   triple (dbartsSpec() resolves one without constructing a sampler)
///   builds the sampler with methods::new("dbartsSampler", control, model,
///   data) and reads the handle back with that object's $getPointer(). The
///   handle is valid until that R object is
///   garbage collected or REPLACES its pointer, which the object does
///   whenever it re-creates its engine from a stored state, so re-read the
///   handle after any R-side restore and keep the R object reachable for as
///   long as you hold one. Everything this header does not reach - creation,
///   the predictor and weight conduits, the state round trip, tree
///   extraction, the multi-forest surface - is an R method on that same
///   object, and the two views drive one engine.
/// - COPY-ON-SET. dbarts_sampler_setResponse and dbarts_sampler_setOffset
///   COPY into buffers the sampler owns and allocated at creation: the
///   caller's array is free the moment the call returns, and no set here
///   allocates. (Changing the observation count or the test row count is an
///   R-side conduit, and it is the one thing that resizes those buffers.) The
///   copy is what the caller conditions on, so a value changed by writing
///   through the array afterwards is not seen at all - call the setter again.
///   Nothing else here retains a pointer either: dbarts_sampler_predict reads
///   its source and its offset for the call alone.
/// - Matrices are column-major. Result and prediction layouts put samples and
///   then chains in trailing dimensions.

// The only headers the prototype view needs: no R header is included here, so
// this file's declarations compile as plain C (gcc -std=c99 -pedantic) with no
// R include path. Nothing below names an R type. The stub view alone reaches
// into R, and includes what it needs where it needs it.
#include <stddef.h> // size_t
#include <stdint.h> // uint32_t

/// The C ABI version, two components. major changes are
/// incompatible; minor changes are additive-only. The safe consumer handshake
/// is major-equality with a minor floor:
///   dbarts_apiMajorVersion() == DBARTS_C_API_MAJOR &&
///   dbarts_apiMinorVersion() >= DBARTS_C_API_MINOR
/// The constants become a compatibility contract at the first release: no
/// version of this API has shipped yet, so whatever they read then simply IS
/// the initial contract, and they do not move before it. A consumer may
/// pre-define either to force a mismatch; nothing but a test of the handshake
/// itself has reason to.
#ifndef DBARTS_C_API_MAJOR
#  define DBARTS_C_API_MAJOR 1
#endif
#ifndef DBARTS_C_API_MINOR
#  define DBARTS_C_API_MINOR 0
#endif

/// FNV-1a token over the ABI this header declares, baked here and
/// static_assert'd against a recomputation in dbarts's own C++ build. Any
/// change it covers fails dbarts's compile until this literal is re-baked;
/// that re-bake is the mechanical acknowledgment of an ABI change. Plain hex
/// literal, usable from C.
///
/// It folds, in fixed order: the stringized DBARTS_C_API_LIST signatures
/// (return types, names, parameter lists), the three ABI enums' enumerator
/// names and values, and the LAYOUT of every struct that crosses the ABI -
/// dbarts_results, dbarts_predictor_source, dbarts_draw - as the compiler
/// reports it: each struct's size, and each field's name paired with its
/// offset, both in pointer units so the token is one number on every supported
/// platform. dbarts_draw is the exception that rule needs: it is the one
/// struct carrying by-value doubles, whose offset in pointer units is not the
/// same number everywhere, so it folds its FIELD COUNT in place of its size
/// and folds those doubles by name, declaration position and width instead of
/// by offset. So a field appended, removed, reordered or retyped to a different
/// width, a renamed field, and a renumbered or added enumerator all move it.
/// What it still does not see is an in-place
/// type swap of the SAME width (double* -> int64_t* under an unchanged name),
/// which moves no offset and no name; that residue is announced to consumers
/// by hand.
///
/// dbarts_apiHash() == DBARTS_C_API_HASH is an OPT-IN exact-ABI check beside
/// the major/minor handshake the stubs always enforce: while the version
/// constants do not move, this token is the only runtime signal that
/// distinguishes a consumer binary built against a different header from one
/// built against this one. A DBARTS_USE_STUBS consumer gets it enforced only
/// when it defines DBARTS_REQUIRE_EXACT_ABI before including this header - the
/// first stub resolution then raises on a mismatch - and one that resolves
/// symbols by hand checks it alongside the major/minor handshake by the same
/// choice.
///
/// It moves on ADDITIVE changes too, since an append moves a struct's size -
/// its field count, for dbarts_draw:
/// after 1.0-0 such an append bumps DBARTS_C_API_MINOR and re-bakes this token
/// together, so a consumer that wants append-COMPATIBILITY gates on
/// major-equality plus a minor floor (plus each struct's structSize contract)
/// and leaves DBARTS_REQUIRE_EXACT_ABI undefined, reserving hash equality for
/// the lockstep EXACTNESS check.
///
/// A consumer may pre-define DBARTS_C_API_HASH to force a mismatch; nothing
/// but a test of the handshake itself has reason to.
#ifndef DBARTS_C_API_HASH
#  define DBARTS_C_API_HASH 0x6380bf095d5cae3fULL
#endif

#ifdef __cplusplus
extern "C" {
#endif

// ---------------------------------------------------------------------------
// ABI types, shared by both the prototype view and the stub view below. Every
// type that crosses the ABI is defined here, so that the single-source list and
// the compile-time token see the whole surface. Nothing here is an R type: an R
// object reaches the sampler through the R methods on the object the handle
// came from, never through this file.
// ---------------------------------------------------------------------------

/// Opaque sampler handle: the address in an R dbartsSampler object's external
/// pointer, read with R_ExternalPtrAddr (see THE HANDLE above). This header
/// creates none and frees none - dbarts_sampler_destroy releases the engine
/// early, and the R object owns what is left.
typedef struct dbarts_sampler_t dbarts_sampler;

/// Caller-owned, growable output buffers for dbarts_sampler_run. The caller
/// MUST set structSize to sizeof(dbarts_results) as the caller compiled it;
/// the library fills only fields whose end offset falls within structSize,
/// so a caller built against an older (smaller) header is never written
/// past. Fields append monotonically below the marked boundary and never
/// reorder across releases; an append after 1.0-0 bumps DBARTS_C_API_MINOR,
/// while a pre-1.0-0 append extends the initial field set and moves no version
/// constant. REMOVING a field is a pre-1.0-0 action only: the initial field
/// set is still being fixed, so a removal shifts every field below it, shrinks
/// sizeof, re-bakes DBARTS_C_API_HASH and moves no version constant - and the
/// structSize contract does NOT cover it, since a stale caller passes a LARGER
/// structSize and the library then fills the shifted fields at the caller's
/// old offsets. Rebuild every consumer against the new header. After 1.0-0 no
/// field is ever removed. A field is
/// filled only when both present-by-size and non-null: a null member skips
/// that quantity, and a zero or unset structSize makes dbarts_sampler_run error
/// rather than silently produce no output. k requires a k
/// hyperprior (dbarts_sampler_kIsSampled), varprobs a DART tree prior
/// (dbarts_sampler_usesDart), dispersion a count (nbinom) response, and
/// residualDf a Student-t residual law; each is
/// left untouched otherwise. logLikelihood
/// carries the per-draw training-data log-likelihood for the gaussian,
/// binary, and aft families; aft reports the log density for events and the
/// log survival tail for right-censored observations. It is
/// NaN-filled wherever the combined per-observation location is not visible to
/// the response model to score - any sampler whose forests combine through
/// amplitudes, at any forest count, and the multinomial softmax - and skipping
/// it (null or absent-by-size) elides all of its computation.
/// varcount reports the sampler's REPORTED forest, which on a multi-forest
/// model is the prognostic forest (the first): this struct declares no forest
/// count, so the engine writes exactly the numPredictors x numSamples x
/// numChains slab documented below whatever the sampler's forest count is. A
/// caller wanting every forest's split counts drives the sampler from R, whose
/// run channel carries a forest axis.
/// Value-initialize with DBARTS_RESULTS_INIT (sets structSize, zeroes the rest):
///   dbarts_results results = DBARTS_RESULTS_INIT;
typedef struct dbarts_results_t {
  size_t structSize;  ///< caller sets to sizeof(dbarts_results)
  double* sigma;      ///< numSamples x numChains
  double* train;      ///< numObservations x numSamples x numChains
  double* test;       ///< numTestObservations x numSamples x numChains
  uint32_t* varcount; ///< numPredictors x numSamples x numChains
  double* k;          ///< numSamples x numChains
  double* varprobs;   ///< numPredictors x numSamples x numChains
  double* logLikelihood; ///< numObservations x numSamples x numChains
  double* dispersion;    ///< numSamples x numChains, the nbinom r per draw
  double* residualDf;    ///< numSamples x numChains, the Student-t nu per draw
  /* 1.0-0 field boundary: every future append goes below this line, never
     above. An append after 1.0-0 bumps DBARTS_C_API_MINOR; a pre-1.0-0 one
     extends the initial field set above it and moves no version constant. */
} dbarts_results;

/// True when the caller's struct (per structSize) actually carries `field`.
/// The sizeof operand is unevaluated, so this never dereferences past the
/// caller's buffer. Every size-first struct here reads its optional members
/// through it, whether the library writes them (dbarts_results) or reads them
/// (dbarts_predictor_source).
#define DBARTS_HAS_FIELD(type, ptr, field) \
  ((ptr)->structSize >= offsetof(type, field) + sizeof((ptr)->field))

/// The dbarts_results spelling of DBARTS_HAS_FIELD.
#define DBARTS_RESULTS_HAS(r, field) DBARTS_HAS_FIELD(dbarts_results, r, field)

/// Value-initializer: sets structSize (the leading member, offset 0) and zeroes
/// the field pointers. Prefer it to hand-setting structSize - dbarts_sampler_run
/// rejects a zero structSize rather than silently producing no output. Every
/// member is spelled out so a consumer building under -Wextra sees no
/// missing-field-initializer warning pointing into this header; an appended
/// field extends this list too.
///   dbarts_results results = DBARTS_RESULTS_INIT;
#define DBARTS_RESULTS_INIT \
  { sizeof(dbarts_results), NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, \
    NULL }

/// One SAVED draw of one chain, handed to a dbarts_draw_callback while the
/// engine holds it. The LIBRARY fills structSize with its own sizeof here (the
/// other direction from dbarts_results, which the caller fills), so a consumer
/// compiled against a NEWER header reads a field the installed library may not
/// have through DBARTS_DRAW_HAS, and one compiled against an older header
/// simply reads the prefix it knows. Fields append monotonically below the
/// marked boundary and never reorder.
///
/// ABSENCE. A channel pointer is NULL wherever the fit does not carry that
/// channel (the variance pair off a heteroscedastic model, the forest pair and
/// glue off a multi-forest coupling, ordinalThresholds off ordinal,
/// splitProbabilities off DART) and wherever the model cannot define it for
/// this draw - the test blend of a coupling with no test treatment vector, the
/// log-likelihood of one whose blended location the response cannot score. A
/// scalar that does not apply is NaN, since a double has no absent value. So a
/// callback tests the channel, never the family.
///
/// LAYOUT, observation fastest throughout: train is numObservations x
/// numReportedLocations (L is 1 on every model but a multi-location one, which
/// folds the offset in at L = 1), test the same over numTestObservations,
/// varianceFits and logLikelihood numObservations, varianceTestFits
/// numTestObservations, forestFits numObservations x numForests forest-major,
/// glue the ragged per-forest amplitude vector numAmplitudes long and
/// forest-major, splitProbabilities numPredictors, ordinalThresholds
/// numOrdinalThresholds, and varcount numPredictors x numVariableCountForests
/// forest-major within the draw (one slab for a single-forest model, K for a
/// multinomial or multi-forest one). Both indices are 0-based, and drawIndex
/// counts the saved draws of THIS run call rather than the sampler's life.
///
/// VALIDITY IS THE CALL AND NO LONGER. A stored channel's pointer is into the
/// caller's own result buffer and outlives the call, but a channel the caller
/// opted out of storing is per-chain SCRATCH that the chain's next draw
/// overwrites, and nothing in this struct tells the two apart. Copy or reduce
/// inside the call.
typedef struct dbarts_draw_t {
  size_t structSize;    ///< the library sets it; read fields through DBARTS_DRAW_HAS
  size_t chainIndex;
  size_t drawIndex;     ///< 0-based over the saved draws of THIS run call
  size_t numObservations;
  size_t numTestObservations;
  size_t numPredictors;
  size_t numReportedLocations;    ///< L: 1, or K for a multi-location model
  size_t numVariableCountForests; ///< varcount slabs in this draw
  size_t numForests;
  size_t numAmplitudes;           ///< glue entries in this draw
  size_t numOrdinalThresholds;
  const double* train;
  const double* test;
  const double* varianceFits;
  const double* varianceTestFits;
  const double* forestFits;
  const double* glue;
  const double* splitProbabilities;
  const double* logLikelihood;
  const double* ordinalThresholds;
  const uint32_t* varcount;
  double sigma;
  double k;
  double dispersion;   ///< NaN off a count (nbinom) response
  double residualDf;   ///< NaN off a Student-t residual law
  /* 1.0-0 field boundary: every future append goes below this line, never
     above, and bumps DBARTS_C_API_MINOR after 1.0-0. */
} dbarts_draw;

/// The dbarts_draw spelling of DBARTS_HAS_FIELD: true when the LIBRARY that
/// filled this draw carries `field`. A consumer reads it only for a field
/// appended after the header it was built against.
#define DBARTS_DRAW_HAS(d, field) DBARTS_HAS_FIELD(dbarts_draw, d, field)

/// A per-draw observer, registered with dbarts_sampler_setDrawCallback and
/// called once per SAVED draw per chain, immediately after the engine settles
/// that draw, with the context registered beside it. A sweep discarded as
/// burn-in never reaches it. It observes: nothing in the draw may be written
/// through, and no sampler state may be mutated from inside it.
///
/// RETURN. 0 continues. NONZERO ABORTS the run: every chain stops at its next
/// sweep boundary, and the sampler is then INCONSISTENT with the results - the
/// sample cursors have not advanced past draws already written into the slots
/// they count - so the caller discards both the results and any saved trees,
/// exactly as on the interrupt path. dbarts_sampler_run reports no status, so
/// a callback that wants to say WHY records it in its own context and the
/// caller reads that after the run; a callback that merely disagrees with a
/// draw records and returns 0.
///
/// CONCURRENCY. Calls for different chains may run CONCURRENTLY, on whichever
/// worker thread owns each chain, and the engine takes NO lock around the
/// call: one would make every chain wait on the slowest callback, which is the
/// cost this mechanism exists to avoid. A callback touching shared state owns
/// its own synchronization; the discipline that needs none is a write
/// addressed by (chainIndex, drawIndex), disjoint by construction. Calls
/// within one chain are ordered by drawIndex.
///
/// NO R API INSIDE THE CALLBACK, EVER - not Rf_allocVector, not PROTECT, not
/// Rf_error, and in C++ not the CONSTRUCTION OR DESTRUCTION of an Rcpp proxy
/// type (Rcpp::NumericVector and its siblings touch the protection stack on
/// both, allocation visible or not; take the raw double* out before the run).
/// R's evaluator, allocator and protection stack are single-threaded, any R
/// allocation may collect objects nothing protected on the worker's behalf,
/// and the callback MUST NOT LONGJMP: Rf_error unwinds a context the main
/// thread established, skipping every C++ destructor between raise and catch
/// even when it is reached from the main thread. An interrupt cannot land
/// while a call is running either, so a callback that blocks hangs the session
/// with no Ctrl-C.
typedef int (*dbarts_draw_callback)(void* context, const dbarts_draw* draw);

/// A predictor column's type. Ordinal columns are cut on their values;
/// categorical ones carry 0-based category codes and split by subset mask;
/// ordered-factor ones carry 0-based level codes whose order is meaningful,
/// so they are cut on those codes rather than masked.
///
/// All three ABI enums are generated from a list macro, which is also what
/// the compile-time token folds: an enumerator ADDED to the enum body
/// directly would be invisible to a hand-kept fold, while one added to the
/// list cannot be. DBARTS_ENUMERATOR is the shared body emitter, #undef'd
/// after the last list below.
#define DBARTS_COLUMN_TYPE_LIST(X) \
  X(DBARTS_COLUMN_ORDINAL, 0) \
  X(DBARTS_COLUMN_CATEGORICAL, 1) \
  X(DBARTS_COLUMN_ORDERED_FACTOR, 2)
#define DBARTS_ENUMERATOR(name, value) name = value,
typedef enum { DBARTS_COLUMN_TYPE_LIST(DBARTS_ENUMERATOR) } dbarts_column_type;

/// A borrowed, self-describing view of predictor values: what
/// dbarts_sampler_predict takes instead of a bare pointer. The caller MUST set
/// structSize to sizeof(dbarts_predictor_source) as it compiled it; the library
/// reads only fields whose end offset falls within structSize, so a caller
/// built against an older (smaller) header is never read past, and a zero
/// structSize is an error rather than a source of silent nulls. Fields append
/// monotonically below the marked boundary and never reorder. Nothing here is
/// retained: the values are replayed during the call and nothing outlives it.
///
/// numRows x numColumns is the shape the argument declares of ITSELF, which is
/// what a caller's own width must agree with - the entry refuses a numColumns
/// that disagrees with the sampler rather than reading whatever lies past the
/// caller's matrix.
///
/// Storage is dense, compressed-column (CSC), or a mix. Without a map
/// (columnSources null) a dense-backed column sits at its own index in
/// whichever channel its kind selects. With one, column j reads dense column
/// columnSources[j] when that is >= 0 and CSC column ~columnSources[j] when
/// it is < 0; the CSC triple is the usual (column pointers of length
/// numCscColumns + 1, row indices, values), and a CSC column's absent rows
/// read its declared reference code when the sampler holds that column
/// categorical, 0 otherwise.
///
/// columnTypes and categoryCounts describe the view's own typing, and the
/// entry that takes this struct is a replay path, where the SAMPLER's store
/// already fixes both for each column it names. They are
/// therefore checked for well-formedness and otherwise IGNORED: no column's
/// type or category count is ever taken from the argument, and declaring them
/// against a sampler that holds otherwise changes nothing. referenceCodes, by
/// contrast, IS read - it is what a CSC column's absent rows decode to.
///
/// Both stay ABOVE the boundary rather than being dropped. They are the only
/// fields a source needs to describe a typing the SAMPLER does not already fix,
/// which is the one thing a predictor view cannot say about itself once they
/// are gone, and this is the last release in which they can be placed rather
/// than appended. They cost a filling caller nothing - null is "nothing
/// declared" - and what a caller does fill is checked, so a malformed
/// declaration is refused here rather than the first time an entry reads one.
///
/// Value-initialize with DBARTS_PREDICTOR_SOURCE_INIT (sets structSize, zeroes
/// the rest), or build the dense case with dbarts_dense_predictor_source():
///   dbarts_predictor_source x = DBARTS_PREDICTOR_SOURCE_INIT;
typedef struct dbarts_predictor_source_t {
  size_t structSize;              ///< caller sets to sizeof(dbarts_predictor_source)
  size_t numRows;
  size_t numColumns;
  const double* denseValues;      ///< column-major, numRows x numColumns
  size_t numCscColumns;           ///< CSC columns; 0 when there is no CSC part
  const int32_t* cscColumnPointers; ///< length numCscColumns + 1
  const int32_t* cscRowIndices;
  const double* cscValues;
  const int32_t* columnSources;   ///< NULL = identity; >= 0 dense col; < 0 CSC col ~v
  const int32_t* columnTypes;     ///< dbarts_column_type per column; declared
                                  ///< and checked; no entry reads it
  const uint32_t* categoryCounts; ///< declared and checked; no entry reads it
  const int32_t* referenceCodes;  ///< per column; < 0 = declared none
  /* 1.0-0 field boundary: appends go below, never above. */
  /// The code channel: a second column-major block, numDenseCodeColumns
  /// columns by numRows rows, holding the int32 level codes of the dense
  /// columns the SAMPLER holds as factors, with INT_MIN (R's NA_INTEGER)
  /// marking a missing code. It saves a caller whose factor columns are
  /// already integers the widening loop and the double block that loop needs.
  /// The replay saves the block as well: it reads the codes where they lie
  /// rather than widening one, dbarts_sampler_predict routing each row by
  /// comparing its code against an integer threshold.
  /// A missing code is admitted only where the sampler's own training column
  /// had missing values; a test source carrying one on a column complete in
  /// training is refused, exactly as a NaN in the double channel is.
  ///
  /// ONE rule decides which channel a column reads, and it is the SAME index
  /// in either. A dense-backed column j (columnSources[j] >= 0, the identity
  /// when columnSources is NULL) that the sampler holds as a factor reads
  /// code column columnSources[j] of this block, which must be less than
  /// numDenseCodeColumns; every other dense-backed column reads double column
  /// columnSources[j] of denseValues, which must be less than numColumns. An
  /// entry refuses a source that breaks either bound, or that leaves NULL a
  /// channel one of its columns needs.
  ///
  /// So under the identity map the code block is indexed by predictor
  /// position and must be wide enough for the last factor column's position,
  /// its other columns unread. A caller that wants both blocks packed over
  /// the columns they serve supplies an explicit columnSources naming each
  /// column's index WITHIN the channel its kind selects. A caller that does
  /// not know which columns the sampler holds as factors leaves denseCodes
  /// NULL, which is exactly what every caller written before this field does.
  const int32_t* denseCodes;
  size_t numDenseCodeColumns;     ///< columns in denseCodes; 0 when there is
                                  ///< no code channel
} dbarts_predictor_source;

/// Value-initializer: sets structSize (the leading member, offset 0) and zeroes
/// everything else, which reads as "dense block, no map, nothing declared".
#define DBARTS_PREDICTOR_SOURCE_INIT \
  { sizeof(dbarts_predictor_source), 0, 0, NULL, 0, NULL, NULL, NULL, NULL, \
    NULL, NULL, NULL, NULL, 0 }

/// The dense spelling: a plain column-major numRows x numColumns block, which
/// is the shape the predictor argument took before this struct existed.
static inline dbarts_predictor_source
dbarts_dense_predictor_source(const double* values, size_t numRows,
                              size_t numColumns) {
  dbarts_predictor_source source = DBARTS_PREDICTOR_SOURCE_INIT;
  source.numRows = numRows;
  source.numColumns = numColumns;
  source.denseValues = values;
  return source;
}

/// The leaf model a sampler's forests carry. No entry here reports one: the
/// enumeration is the ABI's name for what an R-side calibration read answers,
/// kept so the two surfaces number the leaf models alike.
#define DBARTS_LEAF_MODEL_LIST(X) \
  X(DBARTS_LEAF_CONSTANT, 0) \
  X(DBARTS_LEAF_MONOTONE, 1) \
  X(DBARTS_LEAF_LINEAR, 2) \
  X(DBARTS_LEAF_GP, 3)
typedef enum { DBARTS_LEAF_MODEL_LIST(DBARTS_ENUMERATOR) } dbarts_leaf_model;

/// The response family a sampler is built with, or (the R side's default)
/// dispatch on the response shape. 1-6 mirror the engine's own family order;
/// 7 and 8 name two cases that order alone cannot carry: STUDENT is a
/// Student-t residual law riding a gaussian family (dbarts_sampler_family
/// never reports it - a Student-t sampler's family IS gaussian), and
/// MULTINOMIAL is the K-forest softmax. AUTO likewise never comes back: it is
/// what an R-side specification may ask for, and creation resolves it. So
/// dbarts_sampler_family is the only entry here that speaks this enumeration,
/// and it reports 1-6 or 8.
#define DBARTS_FAMILY_LIST(X) \
  X(DBARTS_FAMILY_AUTO, 0) \
  X(DBARTS_FAMILY_GAUSSIAN, 1) \
  X(DBARTS_FAMILY_PROBIT, 2) \
  X(DBARTS_FAMILY_LOGISTIC, 3) \
  X(DBARTS_FAMILY_AFT, 4) \
  X(DBARTS_FAMILY_ORDINAL, 5) \
  X(DBARTS_FAMILY_NBINOM, 6) \
  X(DBARTS_FAMILY_STUDENT, 7) \
  X(DBARTS_FAMILY_MULTINOMIAL, 8)
typedef enum { DBARTS_FAMILY_LIST(DBARTS_ENUMERATOR) } dbarts_family;
#undef DBARTS_ENUMERATOR

// ---------------------------------------------------------------------------
// The single source of truth for the entry-point surface. Each
// entry is X(returnType, name, (parameterList), (argumentList)): the parameter
// list carries names so it also spells the forwarding stub's signature, and the
// argument list forwards those names. Registration (R_interface.cpp), the
// consumer stubs below, the provider-side binding asserts, and the compile-time
// token (C_interface.cpp) are all expansions of this one list, so a signature
// stated here is the only place it is stated. The readable Doxygen prototypes
// kept in the #else branch below are compile-time bound to this list in
// dbarts's own build, so any drift between them fails dbarts's compile.
// ---------------------------------------------------------------------------
#define DBARTS_C_API_LIST(X) \
  X(int, dbarts_apiMajorVersion, (void), ()) \
  X(int, dbarts_apiMinorVersion, (void), ()) \
  X(uint64_t, dbarts_apiHash, (void), ()) \
  X(void, dbarts_sampler_destroy, (dbarts_sampler* sampler), (sampler)) \
  X(void, dbarts_sampler_run, \
    (dbarts_sampler* sampler, size_t numBurnIn, size_t numSamples, \
     dbarts_results* results), \
    (sampler, numBurnIn, numSamples, results)) \
  X(void, dbarts_sampler_sampleTreesFromPrior, (dbarts_sampler* sampler), \
    (sampler)) \
  X(void, dbarts_sampler_setDrawCallback, \
    (dbarts_sampler* sampler, dbarts_draw_callback fn, void* context), \
    (sampler, fn, context)) \
  X(int, dbarts_sampler_setResponse, \
    (dbarts_sampler* sampler, const double* y, int updateScale), \
    (sampler, y, updateScale)) \
  X(int, dbarts_sampler_setOffset, \
    (dbarts_sampler* sampler, const double* offset, int updateScale), \
    (sampler, offset, updateScale)) \
  X(int, dbarts_sampler_setSigma, \
    (dbarts_sampler* sampler, double sigma), (sampler, sigma)) \
  X(int, dbarts_sampler_getLatents, \
    (const dbarts_sampler* sampler, double* out), (sampler, out)) \
  X(int, dbarts_sampler_predict, \
    (dbarts_sampler* sampler, const dbarts_predictor_source* xTest, \
     const double* offsetTest, size_t numThreads, double* out), \
    (sampler, xTest, offsetTest, numThreads, out)) \
  X(void, dbarts_sampler_setTreeStorage, \
    (dbarts_sampler* sampler, int keepTrees, size_t numSamplesToStore), \
    (sampler, keepTrees, numSamplesToStore)) \
  X(void, dbarts_sampler_printTrees, \
    (dbarts_sampler* sampler, size_t forest, const size_t* chainIndices, \
     size_t numChainIndices, const size_t* sampleIndices, \
     size_t numSampleIndices, const size_t* treeIndices, \
     size_t numTreeIndices, int useLiveTrees), \
    (sampler, forest, chainIndices, numChainIndices, sampleIndices, \
     numSampleIndices, treeIndices, numTreeIndices, useLiveTrees)) \
  X(void, dbarts_sampler_setNumThreads, \
    (dbarts_sampler* sampler, size_t numThreads), (sampler, numThreads)) \
  X(void, dbarts_sampler_setVerbose, \
    (dbarts_sampler* sampler, int verbose, size_t printEvery), \
    (sampler, verbose, printEvery)) \
  X(size_t, dbarts_sampler_numObservations, (const dbarts_sampler* sampler), \
    (sampler)) \
  X(size_t, dbarts_sampler_numPredictors, (const dbarts_sampler* sampler), \
    (sampler)) \
  X(size_t, dbarts_sampler_numTestObservations, \
    (const dbarts_sampler* sampler), (sampler)) \
  X(size_t, dbarts_sampler_numChains, (const dbarts_sampler* sampler), \
    (sampler)) \
  X(size_t, dbarts_sampler_numTrees, \
    (const dbarts_sampler* sampler, size_t forest), (sampler, forest)) \
  X(size_t, dbarts_sampler_numSavedSamples, (const dbarts_sampler* sampler), \
    (sampler)) \
  X(int, dbarts_sampler_kIsSampled, (const dbarts_sampler* sampler), \
    (sampler)) \
  X(int, dbarts_sampler_usesDart, (const dbarts_sampler* sampler), (sampler)) \
  X(int, dbarts_sampler_family, (const dbarts_sampler* sampler), (sampler))

/// One stringized "returnType name(parameterList);" per list entry, adjacent
/// string literals that concatenate into the full declaration text the
/// compile-time token hashes.
#define DBARTS_API_STRINGIZE(ret, name, params, args) #ret " " #name #params ";"
#define DBARTS_C_API_DECLS DBARTS_C_API_LIST(DBARTS_API_STRINGIZE)

#ifdef DBARTS_USE_STUBS

// Same-name cached-pointer forwarders generated from the list, one per entry,
// in place of the extern prototypes (the xts inst/include/xtsAPI.h idiom). The
// first call resolves the symbol through R_GetCCallable and caches it in a
// block-scope static; later calls forward directly. A consumer defines
// DBARTS_USE_STUBS to opt in, and then never restates a signature.
#include <R_ext/Error.h>    // Rf_error, raised by the handshake below
#include <R_ext/Rdynload.h> // R_GetCCallable, DL_FUNC

/// The ABI handshake, enforced. Every stub runs this once, on the resolution
/// branch it already has, so it covers every path a stubs consumer can take
/// into the library and costs nothing after the first call of each entry
/// point. The default check is major-equality with a minor floor, resolving
/// dbarts_apiMajorVersion and dbarts_apiMinorVersion raw (calling the stubs
/// would re-enter this function). A consumer that also defines
/// DBARTS_REQUIRE_EXACT_ABI before including this header gets the exact-ABI
/// token checked too, resolved the same raw way: a header the installed
/// dbarts does not carry has the wrong struct layouts, enumerator values and
/// signatures inlined here, so it raises through the same R error path a
/// failed R_GetCCallable resolution takes rather than reading fields the
/// library never wrote. A consumer that wants to negotiate on the
/// major/minor handshake instead of failing resolves those two entries by
/// hand.
static inline void dbarts_stub_checkApi(void) {
  static int dbarts_stub_apiChecked = 0;
  int (*majorFn)(void);
  int (*minorFn)(void);
  int major, minor;
  if (dbarts_stub_apiChecked) return;
  majorFn = (int (*)(void)) (void (*)(void))
    R_GetCCallable("dbarts", "dbarts_apiMajorVersion");
  minorFn = (int (*)(void)) (void (*)(void))
    R_GetCCallable("dbarts", "dbarts_apiMinorVersion");
  major = majorFn();
  minor = minorFn();
  if (major != DBARTS_C_API_MAJOR || minor < DBARTS_C_API_MINOR)
    Rf_error("dbarts C API version mismatch: this package was built against "
             "%d.%d, the installed dbarts provides %d.%d; rebuild this "
             "package against the installed dbarts",
             DBARTS_C_API_MAJOR, DBARTS_C_API_MINOR, major, minor);
#ifdef DBARTS_REQUIRE_EXACT_ABI
  {
    uint64_t (*apiHash)(void) = (uint64_t (*)(void)) (void (*)(void))
      R_GetCCallable("dbarts", "dbarts_apiHash");
    uint64_t installed = apiHash();
    if (installed != (uint64_t) DBARTS_C_API_HASH)
      Rf_error("dbarts C ABI mismatch: this package was built against token "
               "0x%016llx, the installed dbarts reports 0x%016llx; rebuild "
               "this package against the installed dbarts",
               (unsigned long long) DBARTS_C_API_HASH,
               (unsigned long long) installed);
  }
#endif
  dbarts_stub_apiChecked = 1;
}

// void-return detection: ISO C forbids `return <expr>;` in a void function, so
// a void stub must forward without `return`. These helpers pick the right body
// from the entry's return type; all are #undef'd right after the expansion so
// none leak into the consumer's macro namespace.
#define DBARTS_CAT(a, b) DBARTS_CAT_(a, b)
#define DBARTS_CAT_(a, b) a##b
#define DBARTS_CHECK_N(x, n, ...) n
#define DBARTS_CHECK(...) DBARTS_CHECK_N(__VA_ARGS__, 0,)
#define DBARTS_PROBE(x) x, 1,
#define DBARTS_VOID_void DBARTS_PROBE(~)
#define DBARTS_IS_VOID(ret) DBARTS_CHECK(DBARTS_CAT(DBARTS_VOID_, ret))
#define DBARTS_IIF(c) DBARTS_CAT(DBARTS_IIF_, c)
#define DBARTS_IIF_0(t, f) f
#define DBARTS_IIF_1(t, f) t

// The resolution retypes DL_FUNC (void *(*)(void)) to the entry's own
// signature through void (*)(void), the universal function-pointer type the
// -Wextra cast diagnostics exempt (gcc -Wcast-function-type, clang
// -Wcast-function-type-mismatch): cast directly, and every consumer
// translation unit warns once per list entry. The round trip is a retype
// only and generates no code.
#define DBARTS_API_STUB(ret, name, params, args) \
  static inline ret name params { \
    static ret (*dbarts_stub_fn) params = NULL; \
    if (dbarts_stub_fn == NULL) { \
      dbarts_stub_checkApi(); \
      dbarts_stub_fn = \
        (ret (*) params) (void (*)(void)) R_GetCCallable("dbarts", #name); \
    } \
    DBARTS_IIF(DBARTS_IS_VOID(ret))(dbarts_stub_fn args;, \
                                    return dbarts_stub_fn args;) \
  }
DBARTS_C_API_LIST(DBARTS_API_STUB)

#undef DBARTS_API_STUB
#undef DBARTS_IIF_1
#undef DBARTS_IIF_0
#undef DBARTS_IIF
#undef DBARTS_IS_VOID
#undef DBARTS_VOID_void
#undef DBARTS_PROBE
#undef DBARTS_CHECK
#undef DBARTS_CHECK_N
#undef DBARTS_CAT_
#undef DBARTS_CAT

#else // !DBARTS_USE_STUBS: the readable, Doxygen-documented prototypes.

/// Returns DBARTS_C_API_MAJOR of the installed package (incompatible-change
/// component of the version; equal to the caller's for a usable library).
int dbarts_apiMajorVersion(void);
/// Returns DBARTS_C_API_MINOR of the installed package (additive component;
/// at least the caller's for a usable library).
int dbarts_apiMinorVersion(void);
/// Returns the installed package's DBARTS_C_API_HASH: the FNV-1a token over
/// the entry-point signatures, the ABI enums' enumerators, and the ABI
/// structs' compiler-reported layout. Equality with
/// the caller's DBARTS_C_API_HASH says the two were built from the same ABI -
/// the lockstep check, and the only runtime signal that moves while the
/// version constants do not. It moves on additive appends as well, so
/// append-compatibility gates on major/minor and structSize instead (see
/// DBARTS_C_API_HASH). A DBARTS_USE_STUBS consumer does not call it directly:
/// the stubs check it on their first resolution only when the consumer
/// defines DBARTS_REQUIRE_EXACT_ABI before including this header.
uint64_t dbarts_apiHash(void);

/// Releases the engine behind the handle EARLY, before the R object that owns
/// it is collected: for a host that wants a bounded lifetime rather than the
/// collector's. The R object is left in its dead-pointer state, from which its
/// own methods re-create the engine from a stored state or refuse for want of
/// one, so this is a release and not a corruption - but the handle names no
/// sampler afterwards, and passing it to anything else here is the crash a
/// null one is. A second destroy is the one call that is safe: it is a no-op.
void dbarts_sampler_destroy(dbarts_sampler* sampler);

/// Runs numBurnIn discarded then numSamples recorded iterations per chain,
/// recording into results. Set results->structSize before calling; fields whose
/// end offset exceeds it are skipped (never read, never written). Thinning
/// applies within recorded iterations at the rate the control set.
///
/// A null results is legal at any numSamples and records NOTHING: the sweeps
/// still run and the chains advance, which is how a host advances a sampler it
/// is not reading this round. It is not an error, so a null passed by mistake
/// under a positive numSamples returns cleanly with the caller's buffers
/// untouched.
///
/// A callback registered with dbarts_sampler_setDrawCallback fires once per
/// recorded draw of each chain, whether or not results is null - which is how
/// a host reduces draws as they are produced instead of materializing them.
/// One that returns nonzero stops the run early and this entry STILL RETURNS
/// NORMALLY, so a caller that registered one discards these buffers and any
/// saved trees on the status its own context carries.
void dbarts_sampler_run(dbarts_sampler* sampler, size_t numBurnIn,
                        size_t numSamples, dbarts_results* results);
void dbarts_sampler_sampleTreesFromPrior(dbarts_sampler* sampler);

/// Registers a per-draw observer, or clears one with a null fn. fn and context
/// are COPIED into the sampler on the copy-on-set rule above and stay in force
/// for every later dbarts_sampler_run until the next call here; a second
/// registration REPLACES both, and clearing drops the context with the
/// function. Nothing is called through the pointer at registration time, so a
/// wrong one crashes at the first draw of the next run rather than here.
/// context is handed back untouched and may be null; the sampler neither reads
/// nor frees what it points at.
///
/// The callback fires once per SAVED draw per chain (see dbarts_draw_callback
/// for the contract the callback must keep, which is the whole of what makes
/// this safe from a worker thread). A run it aborts returns NORMALLY, having
/// stopped early: this entry and dbarts_sampler_run both report no status, so
/// the caller learns of an abort from its own context and discards the
/// results and any saved trees it holds.
void dbarts_sampler_setDrawCallback(dbarts_sampler* sampler,
                                    dbarts_draw_callback fn, void* context);

/// y has numObservations values, which must lie in the family's support: 0/1
/// for probit and logistic, an integer category index in [1, K] for ordinal, a
/// finite non-negative integer count no larger than 1e6 for nbinom (the
/// dispersion grid's count histogram is sized from the largest count, so a
/// larger one allocates without bound). Out-of-support values are an
/// error, as they are at creation; gaussian and aft (log survival times)
/// constrain nothing. updateScale re-derives the internal response transform
/// from the new response, as dbarts_sampler_setOffset's argument does (gaussian
/// only); pass false once burnt in so fits stay comparable. true is refused on
/// any multi-forest sampler, at any forest count, whose per-forest leaf
/// calibrations are stated against the transform it was built with, and on a
/// heteroscedastic one, whose variance forest is calibrated the same way. The
/// swap itself is refused outright on a coupling that caches per-forest state
/// across sweeps rather than re-deriving it. COPIED, on the copy-on-set rule
/// above: the caller's y is free on return.
///
/// A CAPABILITY STATUS: 1 on a swap, or 0 touching nothing where the coupling
/// admits no response conduit at all. The updateScale refusals above are the
/// other channel and raise.
int dbarts_sampler_setResponse(dbarts_sampler* sampler, const double* y,
                               int updateScale);
/// offset has numObservations values or is null to remove. updateScale
/// rescales the internal response transform to the offset-adjusted range
/// (gaussian only); pass false once burnt in so fits stay comparable. A
/// multi-forest sampler, at any forest count, or a heteroscedastic one refuses
/// true (see setResponse). COPIED, on the copy-on-set rule above: the caller's
/// offset is free on return.
///
/// A CAPABILITY STATUS on dbarts_sampler_setResponse's rule: 1 on a swap, 0
/// where the coupling carries no offset at all.
int dbarts_sampler_setOffset(dbarts_sampler* sampler, const double* offset,
                             int updateScale);
/// Holds the residual standard deviation at sigma (original response scale)
/// until the next call or gaussian draw; the Gibbs conditioning hook. Only a
/// sampler that HAS a residual sd to set takes it: gaussian (Student-t
/// included) and aft. The fixed-unit-scale families (probit, logistic, ordinal,
/// nbinom) pin it at 1 by their own definition and a heteroscedastic sampler's
/// variance forest owns it row by row, so both refuse rather than accept a
/// value nothing would read - read dbarts_sampler_family before calling on a
/// sampler whose family the caller did not choose. A heteroscedastic sampler
/// answers DBARTS_FAMILY_GAUSSIAN and is still refused here, so the accessor
/// does not by itself predict this refusal - but the return value does,
/// without unwinding.
///
/// A CAPABILITY STATUS: 1 on a write, 0 touching nothing on either pinned
/// case.
int dbarts_sampler_setSigma(dbarts_sampler* sampler, double sigma);
/// Copies the current draw of the augmentation variable (numObservations x
/// numChains) into out. Returns 1, or 0 without touching out when the family
/// augments nothing - a plain gaussian response, and a multinomial one.
///
/// WHAT the variable is depends on the family and is not uniform. A LOCATION,
/// on the sampler's own latent scale, for probit (the truncated normal z),
/// ordinal (the same z under the cut points) and aft (the imputed log survival
/// time): a host regresses on these directly. A PRECISION, one per
/// observation, for logistic and nbinom (the Polya-Gamma omega) and for a
/// Student-t residual distribution (the scale-mixing lambda): these WEIGHT a
/// working response and are not on the response scale at all. Note the last
/// case - a sampler whose family is gaussian but whose residual distribution
/// is Student-t DOES report latents, and they are precisions.
///
/// out is written, never read; nothing is retained.
int dbarts_sampler_getLatents(const dbarts_sampler* sampler, double* out);

/// Fits for new data on the original response scale (binary families give
/// the latent scale), from a borrowed source declaring numPredictors columns
/// over the rows to predict. With tree storage out is xTest->numRows x
/// numSavedSamples x numChains from the saved trees; without, one set per
/// chain from the live trees. The saved draws come out OLDEST FIRST - the
/// numSavedSamples most recent recorded draws, however many runs recorded
/// them - so the draw axis is chronological and pairs with the run channels
/// that recorded it, and a store the sampler has recorded nothing into is
/// refused rather than answered from its unwritten slots. offsetTest, when
/// non-null, is added to every sample's fits. A CSC-backed source routes its rows resident, without
/// a dense materialization. Refused on any sampler whose blend is undefined -
/// the predicate is the blend, not the forest count - which is a fixed
/// property of how the sampler was built, so test it once at setup; a host
/// driving such a model reads its forests through the R methods on the object
/// the handle came from.
/// numThreads is a PER-CALL override that does not persist: 0 means the
/// sampler's own count (dbarts_sampler_setNumThreads), and a resolved count
/// below 1 - including a sampler whose own count was set to 0 - is treated as
/// 1. The replay is bitwise identical at every value: the work is partitioned
/// by (chain, draw), each partition owning its output range whole, and nothing
/// is reduced across threads.
///
/// A CAPABILITY STATUS: 1 with out written, or 0 leaving out untouched where
/// the sampler's blend is undefined. An empty tree store and every source
/// refusal are the other channel and raise.
int dbarts_sampler_predict(dbarts_sampler* sampler,
                           const dbarts_predictor_source* xTest,
                           const double* offsetTest, size_t numThreads,
                           double* out);

/// Turns saved-tree storage on or off; numSamplesToStore sizes the buffer
/// when on. Turn on for recorded iterations to predict from them later.
/// Changing either discards what the store held: its recorded-draw count
/// returns to 0.
void dbarts_sampler_setTreeStorage(dbarts_sampler* sampler, int keepTrees,
                                   size_t numSamplesToStore);
/// Prints forest number forest's trees to R's console: pre-order, var 1-based
/// with -1 marking leaves, leaf values on the engine's internal response
/// scale. Every index here is 0-based. Saved trees are read unless
/// useLiveTrees, which ignores sampleIndices; sampleIndices otherwise address
/// RECORDED DRAWS on dbarts_sampler_predict's oldest-first axis, in
/// [0, dbarts_sampler_numSavedSamples), not store slots. Carries no status: an
/// out-of-range forest, chain, sample or tree number is recoverable and
/// raises.
void dbarts_sampler_printTrees(dbarts_sampler* sampler, size_t forest,
                               const size_t* chainIndices,
                               size_t numChainIndices,
                               const size_t* sampleIndices,
                               size_t numSampleIndices,
                               const size_t* treeIndices,
                               size_t numTreeIndices, int useLiveTrees);

void dbarts_sampler_setNumThreads(dbarts_sampler* sampler, size_t numThreads);
/// printEvery counts kept iterations between progress lines and must be at
/// least 1, verbose or not; 0 is an error rather than "never print".
void dbarts_sampler_setVerbose(dbarts_sampler* sampler, int verbose,
                               size_t printEvery);

size_t dbarts_sampler_numObservations(const dbarts_sampler* sampler);
size_t dbarts_sampler_numPredictors(const dbarts_sampler* sampler);
size_t dbarts_sampler_numTestObservations(const dbarts_sampler* sampler);
size_t dbarts_sampler_numChains(const dbarts_sampler* sampler);
/// Forest number forest's tree count; an index past the sampler's last forest
/// is an error, since a size_t probe carries no refusal a caller could tell
/// from a legitimate answer.
size_t dbarts_sampler_numTrees(const dbarts_sampler* sampler, size_t forest);
/// The recorded draws the saved-tree store holds - 0 without tree storage,
/// and 0 until a run records into it. The sample count predict produces and
/// the bound on printTrees' sample indices. It stops at the
/// buffer size: the store is circular, so a longer run keeps the most recent
/// draws.
size_t dbarts_sampler_numSavedSamples(const dbarts_sampler* sampler);
int dbarts_sampler_kIsSampled(const dbarts_sampler* sampler);
int dbarts_sampler_usesDart(const dbarts_sampler* sampler);
/// The dbarts_family this sampler was built with, or the one creation
/// resolved DBARTS_FAMILY_AUTO to: never AUTO (creation always resolves it)
/// and never STUDENT (a Student-t residual sampler's family IS gaussian).
/// Total over every sampler any construction path can build:
/// DBARTS_FAMILY_MULTINOMIAL on a K-forest softmax coupling, the matching
/// enumerator otherwise.
int dbarts_sampler_family(const dbarts_sampler* sampler);

#endif // DBARTS_USE_STUBS

#ifdef __cplusplus
}
#endif

#endif // DBARTS_DBARTS_H
