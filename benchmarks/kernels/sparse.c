// Sparse-column partition prototype (docs/design/sparse-columns.md).
//
// BART's hot per-node operation partitions an arbitrary, scrambled index
// set by a column's codes. This prototypes candidate storage layouts for
// sparse columns and times their partition kernels against the dense u16
// baseline, to settle the representation before any interface freezes
// (docs/design/kernel-vocabulary.md, planned addition 4).
//
// Candidates:
//   dense   - the status quo: n u16 codes, misc_partitionIndices.
//   rank    - bitmap of nonzero rows + per-word cumulative popcounts +
//             packed nonzero codes; code(i) is a bit test, plus a rank
//             lookup only when set. Keeps the unstable two-pointer
//             partition, so nothing about the engine's index handling
//             changes.
//   bsearch - CSC (sorted nonzero rows + codes); code(i) by binary
//             search. The naive layout.
//   merge   - CSC walked in tandem with a SORTED index segment (stable
//             partition into scratch). Hypothetical: every partition in
//             the engine would have to become order-preserving for
//             segments to stay sorted.
//
// Build the package first (src/misc.a); then `make sparse && ./sparse`.
// Output is CSV: kernel,f_nonzero,segment,n,reps,ns_per_elem

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <stdarg.h>

#include <misc/stddef.h>
#include <misc/types.h>
#include <misc/simd.h>
#include <misc/partition.h>

void Rprintf(const char* format, ...) {
  va_list args;
  va_start(args, format);
  vprintf(format, args);
  va_end(args);
}
void R_FlushConsole(void) { fflush(stdout); }

#define N ((size_t) 262144)
#define TARGET_ELEMS ((size_t) 1 << 27)
#define MAX_CODE 250

static volatile size_t sink_s;

static uint64_t rngState = 0x243F6A8885A308D3ull;
static uint64_t nextRand(void) {
  rngState ^= rngState << 13;
  rngState ^= rngState >> 7;
  rngState ^= rngState << 17;
  return rngState;
}

static uint64_t nsecNow(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (uint64_t) ts.tv_sec * 1000000000ull + (uint64_t) ts.tv_nsec;
}

// ----- rank-bitmap layout: 1 bit/row, a cumulative nonzero count per
// 64-row word (4 bytes/64 rows), and 2 bytes per nonzero code
typedef struct {
  uint64_t* bits;      // N / 64 words
  uint32_t* wordRanks; // nonzeros before each word
  misc_xint_t* nzCodes;
  misc_xint_t zeroCode;
} RankColumn;

static inline misc_xint_t rankCodeFor(const RankColumn* col, size_t i) {
  uint64_t word = col->bits[i >> 6];
  uint64_t bit = 1ull << (i & 63u);
  if ((word & bit) == 0) return col->zeroCode;
  size_t rank = col->wordRanks[i >> 6] +
                (size_t) __builtin_popcountll(word & (bit - 1));
  return col->nzCodes[rank];
}

static size_t partitionRank(const RankColumn* col, misc_xint_t cut,
                            size_t* indices, size_t length) {
  size_t lo = 0, hi = length;
  while (1) {
    while (lo < hi && rankCodeFor(col, indices[lo]) <= cut) ++lo;
    while (lo < hi && rankCodeFor(col, indices[hi - 1]) > cut) --hi;
    if (hi - lo < 2) break;
    size_t temp = indices[lo];
    indices[lo] = indices[hi - 1];
    indices[hi - 1] = temp;
    ++lo;
    --hi;
  }
  return lo;
}

// ----- CSC layout: sorted nonzero row ids + codes
typedef struct {
  uint32_t* rows;
  misc_xint_t* codes;
  size_t numNonzero;
  misc_xint_t zeroCode;
} CscColumn;

static inline misc_xint_t cscCodeFor(const CscColumn* col, size_t i) {
  size_t lo = 0, hi = col->numNonzero;
  while (lo < hi) {
    size_t mid = (lo + hi) / 2;
    if (col->rows[mid] < i) lo = mid + 1;
    else hi = mid;
  }
  return lo < col->numNonzero && col->rows[lo] == i ? col->codes[lo]
                                                    : col->zeroCode;
}

static size_t partitionBsearch(const CscColumn* col, misc_xint_t cut,
                               size_t* indices, size_t length) {
  size_t lo = 0, hi = length;
  while (1) {
    while (lo < hi && cscCodeFor(col, indices[lo]) <= cut) ++lo;
    while (lo < hi && cscCodeFor(col, indices[hi - 1]) > cut) --hi;
    if (hi - lo < 2) break;
    size_t temp = indices[lo];
    indices[lo] = indices[hi - 1];
    indices[hi - 1] = temp;
    ++lo;
    --hi;
  }
  return lo;
}

// stable partition of a SORTED segment against the CSC nonzeros in tandem;
// left-bound rows to scratch front, right-bound to scratch back (reversed),
// then copied back preserving order within each side
static size_t partitionMerge(const CscColumn* col, misc_xint_t cut,
                             size_t* indices, size_t length,
                             size_t* scratch) {
  size_t numLeft = 0, numRight = 0;
  size_t nz = 0;
  int zeroGoesLeft = col->zeroCode <= cut;
  for (size_t k = 0; k < length; ++k) {
    size_t i = indices[k];
    while (nz < col->numNonzero && col->rows[nz] < i) ++nz;
    int goesLeft;
    if (nz < col->numNonzero && col->rows[nz] == i)
      goesLeft = col->codes[nz] <= cut;
    else
      goesLeft = zeroGoesLeft;
    if (goesLeft) scratch[numLeft++] = i;
    else scratch[length - ++numRight] = i;
  }
  memcpy(indices, scratch, numLeft * sizeof(size_t));
  // keep ascending order on the right side too
  for (size_t k = 0; k < numRight; ++k)
    indices[numLeft + k] = scratch[length - 1 - k];
  return numLeft;
}

static int compareIndices(const void* a, const void* b) {
  size_t x = *(const size_t*) a, y = *(const size_t*) b;
  return x < y ? -1 : (x > y ? 1 : 0);
}

static misc_xint_t denseCodes[N];
static size_t masterIndices[N];
static size_t sortedIndices[N];
static size_t workIndices[N];
static size_t scratchIndices[N];

int main(void) {
  misc_simd_init();

  const double fractions[] = { 0.5, 0.1, 0.05, 0.01 };
  const size_t segments[] = { N, N / 16, N / 256 };

  printf("# sparse partition prototype, n=%zu, codes 1..%d, cut at median "
         "of nonzeros, zero code 0\n", N, MAX_CODE);
  printf("kernel,f_nonzero,segment,n,reps,ns_per_elem\n");

  for (size_t fi = 0; fi < sizeof(fractions) / sizeof(fractions[0]); ++fi) {
    double f = fractions[fi];

    // dense codes with an approximate nonzero fraction f; nonzero codes
    // uniform in 1..MAX_CODE, zero rows code 0
    size_t numNonzero = 0;
    for (size_t i = 0; i < N; ++i) {
      if ((double) (nextRand() >> 11) * 0x1.0p-53 < f) {
        denseCodes[i] = (misc_xint_t) (1 + nextRand() % MAX_CODE);
        ++numNonzero;
      } else {
        denseCodes[i] = 0;
      }
    }

    RankColumn rank;
    rank.bits = calloc(N / 64, sizeof(uint64_t));
    rank.wordRanks = malloc((N / 64) * sizeof(uint32_t));
    rank.nzCodes = malloc(numNonzero * sizeof(misc_xint_t));
    rank.zeroCode = 0;
    CscColumn csc;
    csc.rows = malloc(numNonzero * sizeof(uint32_t));
    csc.codes = malloc(numNonzero * sizeof(misc_xint_t));
    csc.numNonzero = numNonzero;
    csc.zeroCode = 0;
    size_t nz = 0;
    for (size_t i = 0; i < N; ++i) {
      if (i % 64 == 0) rank.wordRanks[i / 64] = (uint32_t) nz;
      if (denseCodes[i] != 0) {
        rank.bits[i >> 6] |= 1ull << (i & 63u);
        rank.nzCodes[nz] = denseCodes[i];
        csc.rows[nz] = (uint32_t) i;
        csc.codes[nz] = denseCodes[i];
        ++nz;
      }
    }

    // the cut splits the nonzero codes roughly in half; zeros go left
    misc_xint_t cut = MAX_CODE / 2;

    for (size_t si = 0; si < sizeof(segments) / sizeof(segments[0]); ++si) {
      size_t k = segments[si];

      // a scrambled index subset, the shape partitions actually see; the
      // sorted copy serves the hypothetical merge kernel
      for (size_t i = 0; i < N; ++i) masterIndices[i] = i;
      for (size_t i = N - 1; i > 0; --i) {
        size_t j = (size_t) (nextRand() % (i + 1));
        size_t temp = masterIndices[i];
        masterIndices[i] = masterIndices[j];
        masterIndices[j] = temp;
      }
      memcpy(sortedIndices, masterIndices, k * sizeof(size_t));
      qsort(sortedIndices, k, sizeof(size_t), compareIndices);

      size_t reps = TARGET_ELEMS / k;
      if (reps < 8) reps = 8;

      // correctness: every kernel agrees with dense on the left count and
      // the left-side membership sum
      memcpy(workIndices, masterIndices, k * sizeof(size_t));
      size_t leftDense = misc_partitionIndices(denseCodes, cut, workIndices, k);
      uint64_t sumDense = 0;
      for (size_t i = 0; i < leftDense; ++i) sumDense += workIndices[i];
      memcpy(workIndices, masterIndices, k * sizeof(size_t));
      size_t leftRank = partitionRank(&rank, cut, workIndices, k);
      uint64_t sumRank = 0;
      for (size_t i = 0; i < leftRank; ++i) sumRank += workIndices[i];
      memcpy(workIndices, sortedIndices, k * sizeof(size_t));
      size_t leftMerge = partitionMerge(&csc, cut, workIndices, k,
                                        scratchIndices);
      uint64_t sumMerge = 0;
      for (size_t i = 0; i < leftMerge; ++i) sumMerge += workIndices[i];
      if (leftRank != leftDense || sumRank != sumDense ||
          leftMerge != leftDense || sumMerge != sumDense) {
        fprintf(stderr, "kernel disagreement at f=%g k=%zu\n", f, k);
        return 1;
      }

      uint64_t start, elapsed;

      start = nsecNow();
      for (size_t r = 0; r < reps; ++r) {
        memcpy(workIndices, masterIndices, k * sizeof(size_t));
        sink_s = misc_partitionIndices(denseCodes, cut, workIndices, k);
      }
      elapsed = nsecNow() - start;
      printf("dense,%.2f,%zu,%zu,%zu,%.3f\n", f, k, N, reps,
             (double) elapsed / ((double) reps * (double) k));

      start = nsecNow();
      for (size_t r = 0; r < reps; ++r) {
        memcpy(workIndices, masterIndices, k * sizeof(size_t));
        sink_s = partitionRank(&rank, cut, workIndices, k);
      }
      elapsed = nsecNow() - start;
      printf("rank,%.2f,%zu,%zu,%zu,%.3f\n", f, k, N, reps,
             (double) elapsed / ((double) reps * (double) k));

      // binary search is O(k log nnz); keep its reps manageable
      size_t bsReps = reps / 8 > 0 ? reps / 8 : 1;
      start = nsecNow();
      for (size_t r = 0; r < bsReps; ++r) {
        memcpy(workIndices, masterIndices, k * sizeof(size_t));
        sink_s = partitionBsearch(&csc, cut, workIndices, k);
      }
      elapsed = nsecNow() - start;
      printf("bsearch,%.2f,%zu,%zu,%zu,%.3f\n", f, k, N, bsReps,
             (double) elapsed / ((double) bsReps * (double) k));

      start = nsecNow();
      for (size_t r = 0; r < reps; ++r) {
        memcpy(workIndices, sortedIndices, k * sizeof(size_t));
        sink_s = partitionMerge(&csc, cut, workIndices, k, scratchIndices);
      }
      elapsed = nsecNow() - start;
      printf("merge,%.2f,%zu,%zu,%zu,%.3f\n", f, k, N, reps,
             (double) elapsed / ((double) reps * (double) k));
    }

    free(rank.bits);
    free(rank.wordRanks);
    free(rank.nzCodes);
    free(csc.rows);
    free(csc.codes);
  }

  return 0;
}
