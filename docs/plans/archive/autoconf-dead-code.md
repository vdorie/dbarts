# autoconf-dead-code

agent: sonnet
rng: neutral
budget: ~-800 lines (mostly deletions), no behavior change

## Goal

configure advertises only options that build; misc.a ships only code
with callers. The --with-xint-size knob is removed and the u16 width is
pinned and asserted at the bartcore/misc boundary.

## Context

- Knob: [[configure.ac:22-45@418cc6cc]]. bartcore pins the width independently
  ([[src/bartcore/data.hpp:19@9351574f]] `using xint_t = std::uint16_t;`) and passes
  `const xint_t*` uncast into `misc_partitionRange`/`Indices`
  ([[src/bartcore/tree.hpp:658-660@9351574f]]), so non-16 widths fail to compile;
  the bridge hardcodes 65535 ([[src/R_interface_bartcore.cpp:1877@9351574f]]) and
  partition_body.c uses epi16/vdupq_n_u16 intrinsics regardless.
- Dead misc.a members: binaryIO.c (sole WORDS_BIGENDIAN and direct
  posix_memalign consumer; zero callers), adaptiveRadixTree.c (sole
  SIZEOF_SIZE_T and memalign.h consumer; only user is string.c),
  string.c (zero callers). Dead since before the rewrite.
- Unreachable/orphaned configure machinery: the solaris* block
  ([[configure.ac:78-132@418cc6cc]]; Oracle Studio has no C++20, CRAN dropped
  Solaris) with ax_ext_solaris.m4, ax_gcc_x86_cpuid.m4,
  ax_gcc_x86_avx_xgetbv.m4; ax_ext.m4 is never m4_included at all.
- Consumerless checks: HAVE_CSTDINT + stdint fallback, HAVE_SYS_TIME_H,
  AC_TYPE_INT64_T/UINT64_T, ALIGNOF_VOID_P, AC_C_BIGENDIAN,
  snprintf/log1p/namespace-std five (AX_CXX_*, AC_CHECK_FUNCS,
  AC_FUNC_STRERROR_R is LIVE - [[guessNumCores.cpp:211@9351574f]], [[guessNumCores.cpp:247@9351574f]] - keep it).
  After binaryIO/ART/string deletion also: AX_FUNC_POSIX_MEMALIGN,
  HAVE_MALLOC_H, src/include/misc/memalign.h.
- KEEP: AX_COMPILER_EXT/VENDOR + *_FLAG substitutions, AX_PTHREAD
  (std::thread linking + thread managers), AC_C_RESTRICT, alloca,
  unistd, gettimeofday/clock_gettime, ffs, strerror_r. Thread managers
  (thread.c, blockingThreadManager.c, hierarchicalThreadManager.c)
  stay: within-chain-threading names them its substrate.

## Constraints

Out of scope: any change to live kernels or dispatch; Makevars
restructuring beyond removing sun/mapfile branches; the .win files
must stay in sync with the pruned set.

## Steps

1. Delete src/misc/{binaryIO.c,adaptiveRadixTree.c,string.c}, their
   headers under src/include/misc/, and memalign.h; prune the misc
   Makefile object lists (unix and win).
2. Prune configure.ac: xint knob (keep `AC_DEFINE(XINT_TYPE, uint16_t)`
   or hardcode in types.h.in), solaris block, MAPFILE_FLAG, the
   consumerless checks above; delete the four orphaned m4 files.
3. Pin: static_assert(std::is_same_v<misc_xint_t, std::uint16_t>) where
   tree.hpp first uses the kernels; same pin in types.h.win path.
4. Prune src/Makevars.in COMPILER.sun/PTHREAD_*.sun branches and
   @MAPFILE_FLAG@; prune config.h(.in/.win) templates of removed
   defines; configure.win parity.
5. autoreconf -i; ./cleanup; fresh ./configure.
6. Update CLAUDE.local.md's misc.a description (drop "radix tree").

## Verification

- R CMD INSTALL . --preclean succeeds; grep confirms no orphaned
  HAVE_*/deleted symbols remain outside generated files.
- cd tests/cpp && make && ./test_bartcore passes.
- tinytest::test_package("dbarts") passes.
- Equivalence compare vs current baseline reports exact (neutral).

## Landing note (2026-07-07, 091ffb9)

Landed: binaryIO/adaptiveRadixTree/string (and memalign.h) deleted;
the --with-xint-size knob removed with the width pinned by
AC_DEFINE(XINT_TYPE, uint16_t) plus a static_assert at tree.hpp's
first kernel use (mirrored on the .win types path); Solaris block,
MAPFILE_FLAG, and the consumerless checks pruned from configure.ac
and the config templates (.win kept in sync); seven m4 files deleted
(the plan's four Solaris ones by the implementer; three more left
orphaned by the pruning - posix_memalign and the two namespace-std
probes - removed by the reviewer after a zero-reference grep).
autoreconf regenerated configure (-2500 lines net). Net deletion
~-3700 non-generated lines vs the ~-800 estimate: an estimate miss,
not scope creep (the radix tree and binaryIO alone are 2200). Gates:
fresh cleanup/autoreconf/configure/preclean install; zero-hit greps
over the deleted symbols; component tests + fuzzer; tinytest 2468 ok;
equivalence exact 18/18. Windows parity is template-synced but only
compile-verified by CI after the next push. CLAUDE.local.md's misc.a
description updated at landing (step 6, orchestrator-applied).
