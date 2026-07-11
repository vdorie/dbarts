# dbarts x86 AVX2 SIMD implementation plan (READ-ONLY investigation)

Job 120e5d72. Box: dbarts-bench, x86-64 Ubuntu, 16 core, gcc 13.3, R 4.3.3.
CPU: fma sse4_1 sse4_2 avx f16c bmi1 avx2 bmi2 (NO avx512). perf present,
perf_event_paranoid=4 (no perf record without sudo; perf stat user counters
may still work partially). Engine at bartcore HEAD 2e2b1c9.

This memo is the x86-MEASURED follow-up to simd-survey.md (arm64 theory pass).
No repo files modified. Microbenches are throwaway scratch on the box.

## Dispatch reality on THIS box (confirmed from source)
- misc_simd_getMaxSIMDInstructionSet -> AVX2 (has avx2 bit). So:
  - partition family: dispatched to *_avx2 (compiled -mavx2). 8/16-lane.
  - linalg family (add/sub/AXPY/addScalar/setConstant/transpose + Aligned):
    dispatched to *_avx (compiled -mavx). BUT the *_avx bodies are SCALAR
    unroll-by-8 loops with the _mm256 intrinsics COMMENTED OUT
    (linearAlgebra_avx.c:158-526) -> rely on gcc autovec at -mavx = 4-wide,
    and they DROPPED the non-temporal _mm256_stream_pd stores.
  - moments (Mean/Variance): dispatched to *_sse2 only (no avx/avx2/neon variant
    exists). moments_sse2.c.
- The SSE2 linalg bodies (linearAlgebra_sse2.c:158-503) are ALSO scalar
  unroll-by-4, intrinsics + _mm_stream_pd commented out -> autovec at -msse2
  = 2-wide, no NT stores.
- CRITICAL: the bartcore engine headers (chain.hpp etc.) compile into the
  R_interface .o at the PACKAGE Makevars flags = CRAN x86-64 baseline (SSE2,
  no -mavx). So every INLINE elementwise loop in the headers is 2-wide autovec
  on x86, while the dispatched misc kernels it sits next to are 4-wide AVX.
  THIS asymmetry is the whole x86 opportunity. [confirm flags below]

## FINDING 0 (side): 2e2b1c9 does not compile on Linux/gcc (latent, macOS-only-passing)
src/include/external/stats.h sets R_NO_REMAP_RMATH before <Rmath.h>, then
hand-declares only Rf_pnorm5/Rf_pchisq/Rf_qchisq (lines 34-36). It OMITS
Rf_dnorm4 and Rf_qnorm5, which the ext_densityOfNormal/ext_quantileOfNormal
macros (lines 30,32) and model.hpp:2022,2470 use directly. Under R 4.3.3's
Rmath.h on Linux/gcc those symbols are then undeclared -> C_interface.cpp (first
TU) fails: "Rf_dnorm4 was not declared in this scope". macOS clang/R headers
still expose them so the maintainer never sees it. Confirmed both box (2e2b1c9)
and the local worktree copy omit the decls. FIX (2 lines, perf-neutral):
  extern double Rf_dnorm4(double, double, double, int);
  extern double Rf_qnorm5(double, double, double, int, int);
Applied to the BOX CHECKOUT ONLY (declarations; defs come from libR; zero
codegen/perf impact) so profiling can proceed. Repo left untouched. Report to
maintainer as a portability bug independent of the SIMD work.

## Package flag regime (confirmed)
Makevars.in line 60 passes SSE2/SSE4_1/AVX/AVX2/NEON flags ONLY into the
src/{external,misc,rc} sub-make. The main package .o (R_interface.cpp,
C_interface.cpp -> which is where bartcore/*.hpp compiles) get R's default
CXXFLAGS = -g -O2 with NO -mavx/-march (box log line 63 confirms: baseline
x86-64 = SSE2). So bartcore inline loops are 2-wide; misc kernels are 4-wide AVX.

## PIVOT (coordinator correction): suffstat kernel is the PRIORITY target
misc_compute*SufficientStatisticsFast (moments.c:311/334/358/387; 4 variants
contig/indexed x unweighted/weighted) accumulate (Sw, Swz, Swz2) per constant-
leaf node EVERY sweep (tree.hpp:497-513). SCALAR-ONLY, unroll-by-5, single
accumulator, NOT in simd.c dispatch. This is the live hot reduction (NOT the
Mean/Variance family, which is genuinely uncalled - both statements true).
PIVOTAL DELIVERABLE:
 (1) profile-confirm suffstat is a top x86 hotspot;
 (2) microbench a FIXED-LANE bit-identical SIMD (Sw,Swz,Swz2) fused accumulation
     vs current scalar, for BOTH contiguous and indexed/gather, node-n swept
     full->tiny (most nodes small);
 (3) compute-bound (SIMD wins) vs memory/GATHER-bound (low ceiling), reported
     SEPARATELY for contig vs indexed.
BIT-IDENTITY RULE for the microbenched kernel: fixed accumulator lane structure
(identical partial-sum tree + element->lane map on scalar AND every ISA,
independent of width); NO FMA (single-rounding unmatchable in scalar fallback ->
CPU-dependent draws). A fixed-lane scalar != current scalar's order -> a one-time
equivalence re-record is REQUIRED to adopt (maintainer decision). Measure the
fixed-lane SIMD throughput vs current scalar anyway.

## IN-SITU x86 PROFILE (gprof, -pg -O2, THE arm-survey-impossible measurement)
Driver: ConstantLeafSampler, n=10000, p=20, ntree=200, 300 burn + 300 sample =
600 sweeps, gaussian. 3.82s under -pg. Flat profile (self %):
  32.0%  misc_computeIndexedSufficientStatisticsFast   <- suffstat GATHER reduce
  28.7%  Chain::run  (self = the inlined residual-roll loops + orchestration)
  14.9%  misc_setIndexedVectorToConstant               <- fit SCATTER
  14.0%  misc_partitionRange_sse2  }
   8.4%  misc_partitionIndices_sse2} = 22.5% partition (SSE2! see FINDING 1)
   0.6%  misc_setVectorToConstant_avx (root fit fill; AVX ok)
   0.3%  misc_computeSufficientStatisticsFast (contiguous suffstat; root only)
  rest: tree walk / prior / vectors, all <0.6% each.
Top-3 x86 hotspots: (1) indexed suffstat 32%, (2) residual-roll+run 29%,
(3) fit scatter 15%; partition 22.5% is #2.5 but runs at the WRONG ISA.
NOTE the shares are for the AS-SHIPPED regime (buggy detection -> SSE2 partition)
= exactly what CRAN users get. Suffstat stays #1 even after any partition fix.

## FINDING 1 (BIG): AVX2 is MISDETECTED as AVX -> AVX2 partition never runs
misc_simd_getMaxSIMDInstructionSet() returns MISC_INST_AVX(7), NOT AVX2(8), on a
CPU that HAS avx2. Root cause reproduced directly on the box:
  __get_cpuid(7,&a,&b,&c,&d)       -> ebx=0x00000000  (AVX2 bit5 = 0)  WRONG
  __get_cpuid_count(7,0,&a,&b,&c,&d)-> ebx=0x219c91a9  (AVX2 bit5 = 1)  RIGHT
simd.c's cpuid() wrapper (line 101-105) calls __get_cpuid(7,...) for leaf 7, but
leaf 7 needs ecx=subleaf=0 (__get_cpuid does NOT zero ecx). Indeterminate ecx ->
reads a nonzero subleaf -> ebx=0 -> AVX2/AVX512/BMI bits all missed. So dispatch
caps at AVX: linalg gets _avx (fine, that's 256-bit) but PARTITION, whose avx2
branch needs MISC_INST_AVX2, falls to _sse2. The 8-lane u16 AVX2 partition
(partition_avx2.o, built, COMPILER_SUPPORTS_AVX2 defined) is DEAD on this box.
- Bitwise: partition is an exact permutation -> avx2 == sse2 == scalar bit-for-bit
  (that's the equivalence-anchor invariant). Fixing detection needs ZERO re-record.
- This is machine-dependent/latent: works on maintainer macOS (clang zeroes ecx or
  ecx happens 0), fails on this Linux/gcc build. Likely hits many CRAN x86 users
  on AVX2 boxes silently. Fix: give cpuid() a subleaf and use
  __get_cpuid_count(leaf,0,...) (or set ecx=0 before the asm) for leaf>=7.
- Effort: ~5 lines in simd.c, no configure change. Unlocks the ALREADY-WRITTEN
  AVX2 partition for the 22.5% partition share. Quantify via microbench below.

## FINDING 1b (minor): COMPILER_SUPPORTS_SSE4_1 undef'd though flag set + obj built
config.h has /* #undef COMPILER_SUPPORTS_SSE4_1 */ while src/Makevars has
SSE4_1_FLAG=-msse4.1 and partition_sse4_1.o is 22648 bytes (compiled). So the
sse4_1 dispatch branch is #ifdef'd out even though the object exists. Moot once
AVX2 detection works (avx2 supersedes), but it is a configure.ac inconsistency
(the AVX2/AVX/SSE2 macros are defined; SSE4_1 is not). Low priority.

## VD DECISION UPDATE (supersedes fixed-lane): bit-identity is a TOGGLE
Fast path (SIMD/FMA) default; a switch forces the current SCALAR kernel as the
bit-identical reference (tests / exact reproducibility). Fast path must match the
reference only DISTRIBUTIONALLY -> NO cross-ISA fixed-lane constraint. So the
suffstat microbench ladder per variant (contig AND indexed) x node-n (full->tiny):
  tier1 = current scalar (reference)
  tier2 = straightforward max-speed SIMD, natural lane order, NO fma
  tier3 = SIMD + FMA (compile microbench -mfma directly; measure only)
Report each tier's speedup over scalar; CALL OUT tier3-vs-tier2 (the FMA
increment) - VD wants to know if FMA is SIGNIFICANT (decides if FMA ships as a
toggle). Plus compute-bound vs memory/gather-bound per variant.

## SUFFSTAT MICROBENCH (tier1 scalar / tier2 SIMD noFMA / tier3 SIMD+FMA)
Realistic model: full scattered pass over working set N, split into disjoint
nodes of size m (chunks of a random permutation) = one setNodeAverages pass.
ns/element, best of many reps. -O2 -mavx2 -mfma -ffp-contract=off.

INDEXED (the hot 32% variant) -- ns/elem:
 N=2000   m=32 : scalar .351 simd .333 fma .385 hwgather .660
 N=50000  m=32 : scalar .667 simd .607 fma .606 hwgather .943
 N=50000  m=128: scalar .670 simd .586 fma .585 hwgather .770
 N=500000 m=32 : scalar .743 simd .677 fma .662 hwgather .987
 N=500000 m=512: scalar .666 simd .645 fma .643 hwgather .919
 => tier2 SIMD 1.05-1.15x over scalar. tier3 FMA increment ~0 (often negative).
 => HARDWARE AVX2 GATHER (_mm256_i64gather_pd) is WORSE THAN SCALAR everywhere.
 => per-elem cost RISES with N (.35 -> .61 -> .68) = cache-miss/gather signature.
 VERDICT: INDEXED suffstat is MEMORY/GATHER-BOUND. SIMD ceiling ~1.1x; FMA adds
 nothing. Even scalar-indexed is ~2.6x slower/elem than scalar-contiguous purely
 from the gather -> that gap is why it's 32% of runtime, and SIMD can't close it.

CONTIGUOUS (only 0.3% of runtime, root nodes) -- ns/elem:
 N=50000 m=32 : scalar .284 simd .180 fma .159
 N=50000 m=128: scalar .258 simd .153 fma .187
 => tier2 SIMD ~1.6x; tier3 FMA mixed (my SINGLE-accumulator loop is latency-
 bound so FMA's longer dep chain can lose at large m). Compute-bound, SIMD-
 friendly, but negligible share.

FMA VERDICT (VD's question): for the suffstat kernel FMA is NOT significant --
the hot (indexed) path is gather-bound where FMA ~ 0; the compute-bound
(contiguous) path is 0.3% of runtime. FMA does not justify a toggle on suffstat.
NET: SIMD-izing indexed suffstat saves ~10-15% of a 32% kernel = ~3-4% total,
gather-bound. This is the ceiling; it is MODEST, not the big lever people expect.

## PARTITION MICROBENCH (sse2 vs avx2; quantifies the FINDING 1 payoff)
All three variants return the SAME split length L (scalar==sse2==avx2) at every
size -> partition is a bitwise-identical permutation across ISAs; enabling AVX2
needs NO re-record. ns/elem:
 partitionRange (contiguous compare), ncut=100:
   m=512 : scalar 1.065 sse2 .830 avx2 .714  (avx2/sse2 1.16x, sse2/scalar 1.28x)
   m=2048: scalar 1.075 sse2 .790 avx2 .713  (avx2/sse2 1.11x, sse2/scalar 1.36x)
   m=10000:scalar 1.577 sse2 .792 avx2 .718  (avx2/sse2 1.10x, sse2/scalar 1.99x)
 partitionIndices (gather), ncut=100:
   m=128 N=50000 : sse2 .867 avx2 .926 (avx2/sse2 0.94x = WORSE)
   m=512 N=50000 : sse2 1.022 avx2 1.153 (0.89x WORSE)
   m=2048 N=500000: sse2 1.279 avx2 1.320 (0.97x)
VERDICT: the big partition win is scalar->SSE2 (~2x) and is ALREADY captured
(ships at SSE2). AVX2-over-SSE2 is only ~1.10x for Range and <=1.0x (gather-
bound) for Indices. So fixing FINDING 1 buys ~1.1x on the 14% Range share + ~0
on the 8.4% Indices share = ~1-1.5% total. The "22.5% at the wrong ISA" framing
OVERSTATES it: the two-pointer swaps/branches and the gather cap the SIMD, not
compare width. Still worth fixing for CORRECTNESS (+ future AVX2/BMI kernels),
but it is a small speed lever, not a large one.

## RESIDUAL-ROLL MICROBENCH (SSE2 build == CRAN baseline; AVX2 == routed to AVX)
ns/elem, fused resid[i]+=oldF[i]-prevF[i] and the first/2pass forms:
   n=10000 : SSE2 fused .428  AVX2 fused .504  (AVX2 NOT faster - slightly worse)
   n=100000: SSE2 fused .395  AVX2 fused .486
   n=1e6   : SSE2 fused 1.156 AVX2 fused 1.151  (both DRAM-bound)
   fusion: fused .428 vs 2pass .744 at n=1e4 -> 1.7x; DO NOT split into 2 misc calls
VERDICT: the residual roll is MEMORY-BANDWIDTH-BOUND (streams 3-4 n-length
arrays). SSE2==AVX (width irrelevant); at large n it is pure DRAM. The existing
HAND-FUSED single-pass loop is already near-optimal; routing it through the misc
add/sub vocabulary (2 passes) would REGRESS ~1.7x. A fused misc AVX kernel would
match, not beat, the current SSE2 autovec. => residual-roll routing is NOT worth
it. (Chain::run self 28.7% = ~14-15% roll [bandwidth-bound, no headroom] + ~14%
MH-move/orchestration [branchy, not SIMD-able].)

## FIT-SCATTER (misc_setIndexedVectorToConstant, 14.9%)
Scattered stores x[indices[i]]=const. No SSE2/AVX scatter instruction (AVX512
only). Not dispatched, but a SIMD variant can't help pre-AVX512. Memory/scatter-
bound. Only lever is algorithmic (contiguous per-node fit layout). No SIMD win.

## THE HONEST BIG PICTURE
All four top hotspots are MEMORY-bound on the shuffled index buffer:
 indexed suffstat 32% (gather), residual roll ~15% (bandwidth), fit scatter 15%
 (scatter), partitionIndices 8% (gather). SIMD width gives them ~0-1.15x. The
 clean elementwise/permutation SIMD wins are ALREADY captured at SSE2. What's
 left hot is gather/scatter/bandwidth that AVX/AVX2 does not fix. The single
 highest-value lever is NON-SIMD and structural: store per-node y/fits in
 index-buffer (contiguous) order so the gathers/scatters become sequential
 (~2.6x on the 32% suffstat alone, then SIMD-able) -- but that is a sampler
 redesign trading gather-on-read for scatter-on-refresh, out of scope here.

## MULTI-ACCUMULATOR contiguous suffstat (clean compute-bound FMA ceiling)
4 independent accumulators (hides latency), cache-resident, ns/elem:
   m=256 : scalar5 .268  simd4x .076  fma4x .076
   m=1024: scalar5 .264  simd4x .066  fma4x .077
   m=4096: scalar5 .264  simd4x .065  fma4x .079
=> proper-ILP SIMD (no FMA) is ~4.0x over scalar on the compute-bound path.
=> FMA gives NOTHING even here (fma4x == or > simd4x): the loop is LOAD/PORT-
   bound (4x 256-bit loads/16 elem), not FP-ALU-bound, so folding mul+add into
   fma frees ALU slots that were never the bottleneck.
FMA FINAL ANSWER (VD): FMA is NOT a significant lever for suffstat -- 0 on the
gather-bound indexed (hot) path AND 0 on the compute-bound contiguous path with
full ILP. Do not ship FMA as a perf toggle for this kernel.

## MOMENTS VERDICT (reconfirmed on box)
- Mean/Variance/*Unrolled* family (moments_sse2.c, SSE2-only, dispatched): ZERO
  callers in src/ or R/ (grep empty). Completing AVX2/NEON optimizes code the
  sampler never runs. IGNORE (matches my prior; coordinator's correction was
  about the SEPARATE SufficientStatisticsFast family, which IS hot).
- SufficientStatisticsFast family: the live hot reduction (indexed 32%). Verdict
  above: gather-bound, SIMD ~1.1x, FMA 0. Contiguous variant 4x but 0.3%.
- SumOfSquaredResiduals (SSR): DOES have callers (model.hpp:1763 sigma draw via
  drawSigmaSqFromPosterior; chain.hpp:908). Contiguous, compute-bound, SIMD-able
  ~4x, BUT once/sweep (1/ntree the suffstat frequency) -> below profile threshold
  (<0.3%). Negligible alone; trivially reuse a suffstat SIMD kernel if one is built.

================================================================================
## RANKED IMPLEMENTATION PLAN (x86-measured)
Legend: SHARE = est. share of total sampler runtime affected; WIN = measured
kernel speedup; NET = SHARE x realized win; BIT = bit-safety; EFF = effort.

R0. FIX THE BUILD BREAK (external/stats.h) -- SHIP BLOCKER, not perf.
    Add: extern Rf_dnorm4(double,double,double,int);
         extern Rf_qnorm5(double,double,double,int,int);
    2e2b1c9 does NOT compile on Linux/gcc without this (macOS masks it). BIT:
    n/a (decls). EFF: trivial. Independent of all SIMD work. DO FIRST.

R1. FIX AVX2 DETECTION (simd.c cpuid leaf-7 subleaf-0). [FINDING 1]
    Give cpuid() a subleaf param; use __get_cpuid_count(leaf,0,...) for leaf>=7
    (or zero ecx before the asm). Unlocks the already-written AVX2 partition +
    all future AVX2/BMI dispatch. WIN: partitionRange avx2/sse2 ~1.10x; Indices
    ~1.0x (gather). SHARE: 14% Range -> NET ~1-1.5% total. BIT: SAFE (partition
    is an exact permutation; sse2==avx2==scalar split verified). EFF: ~5 lines,
    no configure change, no re-record. RECOMMEND: DO IT -- primarily a
    correctness fix (CRAN x86/AVX2 users silently run SSE2 partition); speed is a
    small bonus. (Also fix R1b: config.h COMPILER_SUPPORTS_SSE4_1 undef'd though
    flag set + object built -- configure.ac inconsistency; low priority.)

R2. SIMD indexed suffstat (tier2 natural-lane, NO FMA), under VD's fast/reference
    TOGGLE. WIN: ~1.05-1.15x (GATHER-BOUND). SHARE: 32%. NET: ~3-4% total.
    BIT: fast path sums in lane order != scalar -> distributional-match only;
    needs the toggle + a ONE-TIME reference re-record. EFF: MED (kernel x ISAs +
    dispatch entry + toggle). Bundle the contiguous variant (4x, compute-bound,
    trivial) and reuse for SSR -- per constant-leaf-suffstat.md. RECOMMEND: only
    if the toggle infra is being built anyway; the payoff is MODEST because the
    hot variant is gather-bound. Do NOT use hardware gather (measured slower).

R3. DO NOT bother (measured dead ends):
    - residual-roll routing: bandwidth-bound, SSE2==AVX; fusion already optimal;
      2-pass misc split would REGRESS 1.7x.
    - fit-scatter SIMD: no pre-AVX512 scatter; memory-bound.
    - moments Mean/Variance AVX/NEON: no callers.
    - FMA anywhere: 0 benefit even compute-bound.

MAINTAINER DECISION TO SURFACE (bigger than any SIMD item): the four top
hotspots are all memory-bound on the shuffled index buffer. The real ~2.6x lever
on the 32% suffstat is STRUCTURAL -- store per-node y/fits in index-buffer
(contiguous) order to turn the gather into a sequential read (then it SIMDs 4x
too) -- at the cost of a scatter on each residual refresh. This is a sampler
redesign, not a kernel, and is the highest-value direction if perf is the goal.
Recommend scoping it separately; the SIMD kernels above are incremental.

TOP 2-3 TO IMPLEMENT: (1) R0 build fix [required]; (2) R1 AVX2 detection fix
[correctness + ~1%]; (3) R2 SIMD indexed suffstat under the toggle [~3-4%, only
if toggle exists]. Everything else is measured not worth it.
================================================================================
- [ ] build on box to private lib; log runtime ISA
- [ ] profile representative sampler run (perf stat / gprof / Rprof)
- [ ] enumerate inline elementwise hot loops
- [ ] microbench residual-roll scalar vs sse2 vs avx (1 pass vs 2)
- [ ] moments caller search + verdict
- [ ] ranked plan
