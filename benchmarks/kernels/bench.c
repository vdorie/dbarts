// Microbenchmarks for the misc kernel vocabulary (docs/design/kernel-vocabulary.md).
//
// Links against the in-tree static library; build the package first so that
// src/misc.a exists (R CMD INSTALL . leaves it behind), then `make && ./bench`.
//
// Output is CSV on stdout: kernel,variant,inst,n,reps,ns_per_elem
// Metadata is emitted as leading '#' comment lines.

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
#include <misc/stats.h>
#include <misc/linearAlgebra.h>

// misc.a's thread manager prints through the R runtime; stub it out here.
void Rprintf(const char* format, ...) {
  va_list args;
  va_start(args, format);
  vprintf(format, args);
  va_end(args);
}
void R_FlushConsole(void) {
  fflush(stdout);
}

static const size_t sizes[] = { 1024, 16384, 262144 };
#define NUM_SIZES (sizeof(sizes) / sizeof(sizes[0]))
#define MAX_N 262144
#define TARGET_ELEMS ((size_t) 1 << 25)
#define MAX_CODE 250

static volatile double sink_d;
static volatile size_t sink_s;

static uint64_t rngState = 0x243F6A8885A308D3ull;
static uint64_t nextRand(void) {
  rngState ^= rngState << 13;
  rngState ^= rngState >> 7;
  rngState ^= rngState << 17;
  return rngState;
}
static double nextUniform(void) {
  return (double) (nextRand() >> 11) * 0x1.0p-53;
}

static uint64_t nsecNow(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (uint64_t) ts.tv_sec * 1000000000ull + (uint64_t) ts.tv_nsec;
}

static const char* instName(misc_simd_instructionSet i) {
#if defined(__i386) || defined(_X86_) || defined(__x86_64__) || defined(_M_AMD64) || defined(_M_X64)
  static const char* const names[] = {
    "C", "SSE", "SSE2", "SSE3", "SSSE3", "SSE4_1", "SSE4_2",
    "AVX", "AVX2", "AVX512F", "AVX512VL", "AVX512BW"
  };
#elif defined(__arm__) || defined(__aarch64__) || defined(_ARM) || defined(_M_ARM)
  static const char* const names[] = { "C", "NEON", "SVE", "SVE2" };
#else
  static const char* const names[] = { "C" };
#endif
  return i < MISC_INST_INVALID ? names[i] : "?";
}

static size_t repsFor(size_t n) {
  size_t reps = TARGET_ELEMS / n;
  return reps < 4 ? 4 : reps;
}

static void report(const char* kernel, const char* variant, const char* inst,
                   size_t n, size_t reps, uint64_t elapsed) {
  printf("%s,%s,%s,%zu,%zu,%.4f\n",
         kernel, variant, inst, n, reps,
         (double) elapsed / ((double) reps * (double) n));
}

// Shared buffers, sized once.
static misc_xint_t codes[MAX_N];
static size_t indices[MAX_N];
static size_t pristineIndices[MAX_N];
static double z[MAX_N], w[MAX_N], y2[MAX_N];

static void fillInputs(void) {
  for (size_t i = 0; i < MAX_N; ++i) {
    codes[i] = (misc_xint_t) (nextRand() % (MAX_CODE + 1));
    z[i] = nextUniform() - 0.5;
    w[i] = nextUniform() + 0.5;
    y2[i] = nextUniform();
    pristineIndices[i] = i;
  }
  // Fisher-Yates so indexed kernels see scattered access.
  for (size_t i = MAX_N - 1; i > 0; --i) {
    size_t j = nextRand() % (i + 1);
    size_t temp = pristineIndices[i];
    pristineIndices[i] = pristineIndices[j];
    pristineIndices[j] = temp;
  }
}

static void benchPartition(const char* inst) {
  static const struct { const char* name; misc_xint_t cut; } cuts[] = {
    { "balanced", (misc_xint_t) (MAX_CODE / 2) },
    { "skewed",   (misc_xint_t) (MAX_CODE / 10) }
  };

  for (size_t c = 0; c < 2; ++c) {
    for (size_t s = 0; s < NUM_SIZES; ++s) {
      size_t n = sizes[s], reps = repsFor(n);

      // partitionRange writes indices from scratch; x is untouched.
      misc_partitionRange(codes, cuts[c].cut, indices, n);
      uint64_t start = nsecNow();
      for (size_t r = 0; r < reps; ++r)
        sink_s = misc_partitionRange(codes, cuts[c].cut, indices, n);
      report("partitionRange", cuts[c].name, inst, n, reps, nsecNow() - start);

      // partitionIndices permutes its input, so restore a shuffled index set
      // each rep; report the memcpy alone so it can be subtracted.
      memcpy(indices, pristineIndices, n * sizeof(size_t));
      misc_partitionIndices(codes, cuts[c].cut, indices, n);
      start = nsecNow();
      for (size_t r = 0; r < reps; ++r) {
        memcpy(indices, pristineIndices, n * sizeof(size_t));
        sink_s = misc_partitionIndices(codes, cuts[c].cut, indices, n);
      }
      report("partitionIndices", cuts[c].name, inst, n, reps, nsecNow() - start);

      start = nsecNow();
      for (size_t r = 0; r < reps; ++r) {
        memcpy(indices, pristineIndices, n * sizeof(size_t));
        sink_s = indices[r % n];
      }
      report("memcpyBaseline", cuts[c].name, inst, n, reps, nsecNow() - start);
    }
  }
}

static void benchMoments(const char* inst) {
  for (size_t s = 0; s < NUM_SIZES; ++s) {
    size_t n = sizes[s], reps = repsFor(n);
    double nEff;

    uint64_t start = nsecNow();
    for (size_t r = 0; r < reps; ++r)
      sink_d = misc_computeMean(z, n);
    report("computeMean", "full", inst, n, reps, nsecNow() - start);

    start = nsecNow();
    for (size_t r = 0; r < reps; ++r)
      sink_d = misc_computeIndexedMean(z, pristineIndices, n);
    report("computeIndexedMean", "shuffled", inst, n, reps, nsecNow() - start);

    start = nsecNow();
    for (size_t r = 0; r < reps; ++r)
      sink_d = misc_computeWeightedMean(z, n, w, &nEff);
    report("computeWeightedMean", "full", inst, n, reps, nsecNow() - start);

    start = nsecNow();
    for (size_t r = 0; r < reps; ++r)
      sink_d = misc_computeVarianceForKnownMean(z, n, 0.0);
    report("computeVarianceForKnownMean", "full", inst, n, reps, nsecNow() - start);

    start = nsecNow();
    for (size_t r = 0; r < reps; ++r)
      sink_d = misc_computeIndexedWeightedVarianceForKnownMean(z, pristineIndices, n, w, 0.0);
    report("computeIndexedWeightedVarianceForKnownMean", "shuffled", inst, n, reps, nsecNow() - start);
  }
}

static void benchVectorOps(const char* inst) {
  for (size_t s = 0; s < NUM_SIZES; ++s) {
    size_t n = sizes[s], reps = repsFor(n);

    uint64_t start = nsecNow();
    for (size_t r = 0; r < reps; ++r)
      misc_addVectorsInPlace(z, n, y2);
    sink_d = y2[n - 1];
    report("addVectorsInPlace", "full", inst, n, reps, nsecNow() - start);

    start = nsecNow();
    for (size_t r = 0; r < reps; ++r)
      misc_subtractVectorsInPlace(z, n, y2);
    sink_d = y2[n - 1];
    report("subtractVectorsInPlace", "full", inst, n, reps, nsecNow() - start);
  }
}

int main(void) {
  misc_simd_init();
  fillInputs();

  misc_simd_instructionSet maxSet = misc_simd_getMaxSIMDInstructionSet();

  misc_simd_instructionSet candidates[8];
  size_t numCandidates = 0;
  candidates[numCandidates++] = MISC_INST_C;
#if defined(__i386) || defined(_X86_) || defined(__x86_64__) || defined(_M_AMD64) || defined(_M_X64)
  misc_simd_instructionSet intermediates[] = { MISC_INST_SSE2, MISC_INST_SSE4_1, MISC_INST_AVX2 };
  for (size_t i = 0; i < 3; ++i)
    if (intermediates[i] < maxSet) candidates[numCandidates++] = intermediates[i];
#endif
  if (maxSet > MISC_INST_C) candidates[numCandidates++] = maxSet;

  printf("# xint bytes: %zu, max instruction set: %s\n",
         sizeof(misc_xint_t), instName(maxSet));
  printf("kernel,variant,inst,n,reps,ns_per_elem\n");

  for (size_t c = 0; c < numCandidates; ++c) {
    misc_simd_setSIMDInstructionSet(candidates[c]);
    const char* inst = instName(candidates[c]);
    benchPartition(inst);
    benchMoments(inst);
    benchVectorOps(inst);
  }

  (void) sink_d; (void) sink_s;
  return 0;
}
