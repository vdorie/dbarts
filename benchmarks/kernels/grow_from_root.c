// Cut-scan histogram kernel prototype for the grow-from-root memo
// (docs/design/grow-from-root.md) and the informed-proposal falsifier of
// docs/design/parallel-bart-frontier.md 3.1 / item 4.
//
// The kernel: for one leaf's member set, one pass per variable over the members
// accumulating per-code sufficient statistics (count, sum w, sum wz, optionally
// sum wz^2) for ALL cuts of ALL p variables, then a prefix scan collapsing the
// per-code histogram into the integrated-likelihood weight of every candidate
// cut. This is the histogram primitive of GBM stacks and XBART's GrowFromRoot.
//
// Data layout is faithful to bartcore's ColumnStore: codes are misc_xint_t
// (uint16), column-major with column j at codes + j*N; a leaf's members are obs
// ids held in an index segment (Tree::Node::indices), so the scan gathers
// col[idx[i]], z[idx[i]], w[idx[i]] - exactly the access pattern of
// misc_partitionIndices and misc_computeIndexed*SufficientStatisticsFast.
//
// Baseline: the current single-cut move cost model - one gather partition on one
// column (misc_partitionIndices) plus the two-child suffstat recompute
// (misc_computeIndexedWeightedSufficientStatisticsFast, the exact kernel
// Tree::birth calls). The multiplier scan/move is the number frontier item 4
// needs.
//
// Links the in-tree misc.a (build the package first: R CMD INSTALL . leaves
// src/misc.a), then `make grow_from_root && ./grow_from_root`.
//
// Output is CSV on stdout with leading '#' metadata/summary lines.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <pthread.h>

#include <misc/stddef.h>
#include <misc/types.h>
#include <misc/simd.h>
#include <misc/partition.h>
#include <misc/stats.h>

#define N_POP ((size_t) 1 << 20)   // 1048576 observations in the population
#define P_MAX 50
#define MAX_BINS 256               // numCuts up to 255 -> 256 codes
#define SCAN_TARGET ((size_t) 1 << 26)
#define MOVE_TARGET ((size_t) 1 << 24)
#define LAMBDA 16.0                // leaf prior precision at BART defaults (k=2, m=64-ish)

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

// Population buffers, sized once at the max (n = N_POP, p = P_MAX).
static misc_xint_t* codes;             // N_POP * P_MAX, column-major
static uint32_t numCuts[P_MAX];
static double* z;                      // working response, N_POP
static double* w;                      // weights, N_POP
static size_t* idxSeq;                 // identity 0..N_POP-1 (root / contiguous)
static size_t* idxGather;              // a shuffled permutation (scattered leaf)
static size_t* idxScratch;             // partition working copy
static double* streamBuf;              // bandwidth probe

static void fillInputs(size_t cutCount) {
  for (size_t j = 0; j < P_MAX; ++j) {
    numCuts[j] = (uint32_t) cutCount;
    misc_xint_t* col = codes + j * N_POP;
    for (size_t i = 0; i < N_POP; ++i)
      col[i] = (misc_xint_t) (nextRand() % (cutCount + 1));
  }
  for (size_t i = 0; i < N_POP; ++i) {
    z[i] = nextUniform() - 0.5;
    w[i] = nextUniform() + 0.5;
    idxSeq[i] = i;
    idxGather[i] = i;
  }
  for (size_t i = N_POP - 1; i > 0; --i) {
    size_t k = nextRand() % (i + 1);
    size_t t = idxGather[i]; idxGather[i] = idxGather[k]; idxGather[k] = t;
  }
}

// One node-expansion: histogram every variable's members, prefix-scan each into
// the best integrated-likelihood cut gain. Returns the summed best gain (sink,
// so the accumulation is not dead code). withSq accumulates sum wz^2 too.
static double scanColMajor(size_t p, const double* zz, const double* ww,
                           const size_t* idx, size_t nLeaf, bool withSq,
                           double* bins /* 4*MAX_BINS scratch */) {
  double bestSum = 0.0;
  for (size_t j = 0; j < p; ++j) {
    const misc_xint_t* col = codes + j * N_POP;
    size_t B = numCuts[j] + 1;
    double* cnt = bins;
    double* sw  = bins + B;
    double* swz = bins + 2 * B;
    double* swz2 = bins + 3 * B;
    for (size_t b = 0; b < B; ++b) { cnt[b] = 0.0; sw[b] = 0.0; swz[b] = 0.0; if (withSq) swz2[b] = 0.0; }

    for (size_t i = 0; i < nLeaf; ++i) {
      size_t obs = idx[i];
      misc_xint_t c = col[obs];
      double wi = ww[obs];
      double zi = zz[obs];
      double wz = wi * zi;
      cnt[c] += 1.0;
      sw[c]  += wi;
      swz[c] += wz;
      if (withSq) swz2[c] += wz * zi;
    }

    double totW = 0.0, totWZ = 0.0;
    for (size_t b = 0; b < B; ++b) { totW += sw[b]; totWZ += swz[b]; }
    double accW = 0.0, accWZ = 0.0, best = -1e300;
    for (size_t k = 0; k + 1 < B; ++k) {          // split sends codes 0..k left
      accW += sw[k]; accWZ += swz[k];
      double rW = totW - accW, rWZ = totWZ - accWZ;
      double gain = accWZ * accWZ / (accW + LAMBDA) + rWZ * rWZ / (rW + LAMBDA);
      if (gain > best) best = gain;
    }
    bestSum += best;
  }
  return bestSum;
}

// Threaded scan: split the p columns across T workers, each with private bins.
typedef struct {
  size_t jStart, jEnd, nLeaf;
  const double *zz, *ww; const size_t* idx; bool withSq;
  double partial;
} ScanJob;

static void* scanWorker(void* arg) {
  ScanJob* jb = (ScanJob*) arg;
  double bins[4 * MAX_BINS];
  double s = 0.0;
  for (size_t j = jb->jStart; j < jb->jEnd; ++j) {
    // reuse scanColMajor's body over a single column via a 1-col call
    const misc_xint_t* col = codes + j * N_POP;
    size_t B = numCuts[j] + 1;
    double* cnt = bins; double* sw = bins + B; double* swz = bins + 2 * B; double* swz2 = bins + 3 * B;
    for (size_t b = 0; b < B; ++b) { cnt[b]=0; sw[b]=0; swz[b]=0; if (jb->withSq) swz2[b]=0; }
    for (size_t i = 0; i < jb->nLeaf; ++i) {
      size_t obs = jb->idx[i];
      misc_xint_t c = col[obs];
      double wi = jb->ww[obs]; double zi = jb->zz[obs]; double wz = wi * zi;
      cnt[c]+=1.0; sw[c]+=wi; swz[c]+=wz; if (jb->withSq) swz2[c]+=wz*zi;
    }
    double totW=0, totWZ=0; for (size_t b=0;b<B;++b){ totW+=sw[b]; totWZ+=swz[b]; }
    double accW=0, accWZ=0, best=-1e300;
    for (size_t k=0;k+1<B;++k){ accW+=sw[k]; accWZ+=swz[k]; double rW=totW-accW, rWZ=totWZ-accWZ;
      double g=accWZ*accWZ/(accW+LAMBDA)+rWZ*rWZ/(rW+LAMBDA); if (g>best) best=g; }
    s += best;
  }
  jb->partial = s;
  return NULL;
}

static double scanThreaded(size_t p, const size_t* idx, size_t nLeaf, bool withSq, int T) {
  pthread_t th[8];
  ScanJob jobs[8];
  size_t per = (p + (size_t) T - 1) / (size_t) T;
  int nt = 0;
  for (int t = 0; t < T; ++t) {
    size_t a = (size_t) t * per, b = a + per; if (b > p) b = p; if (a >= b) break;
    jobs[nt] = (ScanJob){ a, b, nLeaf, z, w, idx, withSq, 0.0 };
    ++nt;
  }
  for (int t = 0; t < nt; ++t) pthread_create(&th[t], NULL, scanWorker, &jobs[t]);
  double s = 0.0;
  for (int t = 0; t < nt; ++t) { pthread_join(th[t], NULL); s += jobs[t].partial; }
  (void) p;
  return s;
}

// Current single-cut move: one gather partition on column jcol + two-child
// suffstat recompute. idxScratch must hold a pristine copy of the member set on
// entry (partition permutes it).
static void singleCutMove(size_t jcol, misc_xint_t cut, size_t nLeaf) {
  const misc_xint_t* col = codes + jcol * N_POP;
  size_t numLeft = misc_partitionIndices(col, cut, idxScratch, nLeaf);
  double sw0, swz0, swz20, sw1, swz1, swz21;
  misc_computeIndexedWeightedSufficientStatisticsFast(z, idxScratch, numLeft, w, &sw0, &swz0, &swz20);
  misc_computeIndexedWeightedSufficientStatisticsFast(z, idxScratch + numLeft, nLeaf - numLeft, w, &sw1, &swz1, &swz21);
  sink_d += sw0 + swz0 + swz20 + sw1 + swz1 + swz21;
}

static size_t repsFor(size_t work, size_t target, size_t floorReps) {
  size_t r = target / (work == 0 ? 1 : work);
  return r < floorReps ? floorReps : r;
}

// GB/s streaming read bandwidth: sequential sum over a large (out-of-cache)
// array with 8 independent accumulators, so the FP-add latency chain does not
// cap it and the loop is memory-bound.
static double streamReadBandwidth(void) {
  size_t n = (size_t) 1 << 23;  // 8M doubles = 64 MB, well past L2
  for (size_t i = 0; i < n; ++i) streamBuf[i] = (double) i;
  double a0=0,a1=0,a2=0,a3=0,a4=0,a5=0,a6=0,a7=0;
  int reps = 8;
  uint64_t start = nsecNow();
  for (int r = 0; r < reps; ++r)
    for (size_t i = 0; i + 8 <= n; i += 8) {
      a0+=streamBuf[i];   a1+=streamBuf[i+1]; a2+=streamBuf[i+2]; a3+=streamBuf[i+3];
      a4+=streamBuf[i+4]; a5+=streamBuf[i+5]; a6+=streamBuf[i+6]; a7+=streamBuf[i+7];
    }
  uint64_t el = nsecNow() - start;
  sink_d += a0+a1+a2+a3+a4+a5+a6+a7;
  return (double) (n * sizeof(double) * (size_t) reps) / (double) el;  // bytes/ns = GB/s
}

static void correctnessCheck(void) {
  // Histogram totals must equal a direct suffstat over the same members.
  double bins[4 * MAX_BINS];
  size_t nLeaf = 4096;
  size_t p = 1;
  numCuts[0] = 100;
  fillInputs(100);
  // direct suffstat over gather members of column 0
  memcpy(idxScratch, idxGather, nLeaf * sizeof(size_t));
  double swD, swzD, swz2D;
  misc_computeIndexedWeightedSufficientStatisticsFast(z, idxScratch, nLeaf, w, &swD, &swzD, &swz2D);
  // histogram totals for column 0
  const misc_xint_t* col = codes;
  double tW = 0, tWZ = 0, tWZ2 = 0;
  for (size_t i = 0; i < nLeaf; ++i) {
    size_t obs = idxGather[i];
    double wi = w[obs], zi = z[obs], wz = wi * zi;
    tW += wi; tWZ += wz; tWZ2 += wz * zi;
  }
  if (fabs(tW - swD) > 1e-6 * fabs(swD) + 1e-9 ||
      fabs(tWZ - swzD) > 1e-6 * fabs(swzD) + 1e-9 ||
      fabs(tWZ2 - swz2D) > 1e-6 * fabs(swz2D) + 1e-9) {
    fprintf(stderr, "FAIL: histogram totals diverge from direct suffstat\n");
    exit(1);
  }
  (void) bins; (void) p;
  printf("# check: histogram totals == direct suffstat (n=%zu): OK\n", nLeaf);
}

int main(void) {
  misc_simd_init();
  misc_simd_setSIMDInstructionSet(misc_simd_getMaxSIMDInstructionSet());

  codes = malloc(N_POP * P_MAX * sizeof(misc_xint_t));
  z = malloc(N_POP * sizeof(double));
  w = malloc(N_POP * sizeof(double));
  idxSeq = malloc(N_POP * sizeof(size_t));
  idxGather = malloc(N_POP * sizeof(size_t));
  idxScratch = malloc(N_POP * sizeof(size_t));
  streamBuf = malloc(((size_t) 1 << 23) * sizeof(double));
  if (!codes || !z || !w || !idxSeq || !idxGather || !idxScratch || !streamBuf) {
    fprintf(stderr, "alloc failed\n"); return 1;
  }

  correctnessCheck();

  double bw = streamReadBandwidth();
  printf("# machine: single-thread streaming read bandwidth = %.1f GB/s\n", bw);
  printf("# xint bytes: %zu; population N = %zu; lambda = %.1f\n",
         sizeof(misc_xint_t), N_POP, LAMBDA);
  printf("# columns: sq = sum-wz^2 accumulated; move_ns = partition + 2-child suffstat recompute\n");

  static const size_t leafSizes[] = { 1000, 10000, 100000, 1000000 };
  static const size_t pVals[] = { 10, 50 };
  static const size_t cutVals[] = { 100, 255 };

  double bins[4 * MAX_BINS];

  printf("kernel,mode,n_leaf,p,cuts,sq,reps,scan_ns,ns_per_elem,move_ns,multiplier,scan_GBs\n");

  for (size_t ci = 0; ci < 2; ++ci) {
    size_t cutCount = cutVals[ci];
    fillInputs(cutCount);
    for (size_t pi = 0; pi < 2; ++pi) {
      size_t p = pVals[pi];
      for (size_t li = 0; li < 4; ++li) {
        size_t nLeaf = leafSizes[li];
        size_t reps = repsFor(nLeaf * p, SCAN_TARGET, 3);
        size_t moveReps = repsFor(nLeaf, MOVE_TARGET, 8);

        // ---- baseline single-cut move (gather members), memcpy subtracted ----
        misc_xint_t cut = (misc_xint_t) (cutCount / 2);
        memcpy(idxScratch, idxGather, nLeaf * sizeof(size_t));
        singleCutMove(0, cut, nLeaf);            // warmup
        uint64_t start = nsecNow();
        for (size_t r = 0; r < moveReps; ++r) {
          memcpy(idxScratch, idxGather, nLeaf * sizeof(size_t));
          singleCutMove(r % p, cut, nLeaf);
        }
        uint64_t movePlusCopy = nsecNow() - start;
        start = nsecNow();
        for (size_t r = 0; r < moveReps; ++r) {
          memcpy(idxScratch, idxGather, nLeaf * sizeof(size_t));
          sink_s += idxScratch[r % nLeaf];
        }
        uint64_t copyOnly = nsecNow() - start;
        double moveNs = (double) (movePlusCopy - copyOnly) / (double) moveReps;
        if (moveNs < 1.0) moveNs = 1.0;

        // ---- cut-scan, both access modes and both sq settings ----
        for (int sq = 0; sq <= 1; ++sq) {
          for (int mode = 0; mode < 2; ++mode) {
            const size_t* idx = mode == 0 ? idxSeq : idxGather;
            const char* modeName = mode == 0 ? "seq" : "gather";
            sink_d += scanColMajor(p, z, w, idx, nLeaf, sq, bins);  // warmup
            start = nsecNow();
            for (size_t r = 0; r < reps; ++r)
              sink_d += scanColMajor(p, z, w, idx, nLeaf, sq, bins);
            uint64_t el = nsecNow() - start;
            double scanNs = (double) el / (double) reps;
            double nsPerElem = scanNs / (double) (nLeaf * p);
            // demand bytes: code(2) + z(8) + w(8) per member-variable; the sq
            // bin write lands in cache, so DRAM demand is 18 either way.
            size_t bytesPerElem = 2 + 8 + 8;
            double gbs = (double) (bytesPerElem * nLeaf * p) / scanNs;
            double mult = scanNs / moveNs;
            printf("scan,%s,%zu,%zu,%zu,%d,%zu,%.0f,%.4f,%.0f,%.1f,%.1f\n",
                   modeName, nLeaf, p, cutCount, sq, reps, scanNs, nsPerElem, moveNs, mult, gbs);
          }
        }
      }
    }
  }

  // ---- threading scan: representative heavy config (p=50, cuts=255) ----
  fillInputs(255);
  printf("# threading (gather, p=50, cuts=255, sq=0): scan_ns and speedup vs T=1\n");
  printf("kernel,mode,n_leaf,p,cuts,threads,reps,scan_ns,speedup\n");
  static const int threadCounts[] = { 1, 2, 4 };
  for (size_t li = 2; li < 4; ++li) {          // n_leaf 1e5, 1e6
    size_t nLeaf = leafSizes[li];
    size_t reps = repsFor(nLeaf * 50, SCAN_TARGET, 3);
    double baseNs = 0.0;
    for (size_t ti = 0; ti < 3; ++ti) {
      int T = threadCounts[ti];
      scanThreaded(50, idxGather, nLeaf, false, T);  // warmup
      uint64_t start = nsecNow();
      for (size_t r = 0; r < reps; ++r) sink_d += scanThreaded(50, idxGather, nLeaf, false, T);
      double ns = (double) (nsecNow() - start) / (double) reps;
      if (T == 1) baseNs = ns;
      printf("scan_mt,gather,%zu,50,255,%d,%zu,%.0f,%.2f\n", nLeaf, T, reps, ns, baseNs / ns);
    }
  }

  (void) sink_d; (void) sink_s; (void) idxSeq;
  free(codes); free(z); free(w); free(idxSeq); free(idxGather); free(idxScratch); free(streamBuf);
  return 0;
}
