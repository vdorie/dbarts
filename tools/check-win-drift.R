#!/usr/bin/env Rscript

# Windows has no configure: src/*.win files are checked in by hand and
# nothing regenerates them, so they drift silently from the autoconf
# templates (*.in) and from DESCRIPTION. Checks (a) every PACKAGE_STRING /
# PACKAGE_VERSION literal in the config *.win headers against DESCRIPTION's
# Version, and (b) the macro nameset of each config *.win header against its
# *.in template, allowing only macros on the RECORDED expected-absent table
# below (macros the Windows build deliberately does not define). A macro
# absent from a .win header, present in its .in template, and not on the
# table is silent drift (a new configure check missing on Windows); a listed
# macro that has since appeared in its .win header, or dropped from the .in
# template, is a stale table entry.
#
# Usage: Rscript tools/check-win-drift.R [pkg-root]

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) args[1L] else "."
p <- function(...) file.path(root, ...)

# Active (uncommented) "#define NAME" macros, or "#undef NAME" placeholders.
macroNames <- function(file, directive) {
  lines <- trimws(readLines(file))
  lines <- lines[!startsWith(lines, "/*")]
  re <- sprintf("^#[ ]*%s[ ]+([A-Za-z_][A-Za-z0-9_]*)", directive)
  hit <- regmatches(lines, regexec(re, lines))
  unique(vapply(hit[lengths(hit) == 2L], `[[`, character(1L), 2L))
}

# file|comma-separated macros|reason. Populated from the current files (see
# each reason for the forensics behind it); UNKNOWN is an acceptable reason
# but none are needed here.
absentRaw <- c(
  "src/config.hpp.win|COMPILER_SUPPORTS_AVX,COMPILER_SUPPORTS_AVX2,COMPILER_SUPPORTS_NEON,COMPILER_SUPPORTS_SSE2,COMPILER_SUPPORTS_SSE4_1|SIMD dispatch reads misc/config.h, not config.hpp (src/misc/simd.c, src/include/misc/intrinsic.h)",
  "src/config.hpp.win|XINT_TYPE|provided instead by types.h.win's own ifndef-guarded default",
  "src/config.hpp.win|STRERROR_R_CHAR_P|read only in guessNumCores.cpp's #elif defined(__linux__) branch, unreachable under #ifdef _WIN32",
  "src/config.hpp.win|C_ALLOCA,HAVE_ALLOCA,HAVE_ALLOCA_H,HAVE_CLOCK_GETTIME,HAVE_DECL_STRERROR_R,HAVE_FFS,HAVE_GETTIMEOFDAY,HAVE_INTTYPES_H,HAVE_PTHREAD,HAVE_PTHREAD_PRIO_INHERIT,HAVE_STDIO_H,HAVE_STDLIB_H,HAVE_STRERROR_R,HAVE_STRINGS_H,HAVE_STRING_H,HAVE_SYS_STAT_H,HAVE_SYS_TYPES_H,HAVE_UNISTD_H,PACKAGE_URL,PTHREAD_CREATE_JOINABLE,STACK_DIRECTION,restrict,size_t|unused via config.hpp: no dbarts source reads this macro from config.hpp",
  "src/external/config.h.win|HAVE_UNISTD_H|no unistd.h on Windows (POSIX-only header)",
  "src/external/config.h.win|HAVE_CLOCK_GETTIME|no POSIX clock_gettime; falls back to HAVE_GETTIMEOFDAY (src/external/randomBase.c)",
  "src/external/config.h.win|HAVE_ALLOCA_H|no alloca.h on Windows; falls back to _alloca/malloc.h or __builtin_alloca",
  "src/external/config.h.win|HAVE_FFS|no BSD ffs() on Windows; unreached in intrinsic.h's #elif chain (see misc/config.h.win entry)",
  "src/misc/config.h.win|HAVE_UNISTD_H|no unistd.h on Windows (POSIX-only header)",
  "src/misc/config.h.win|HAVE_CLOCK_GETTIME|no POSIX clock_gettime; falls back to HAVE_GETTIMEOFDAY (src/misc/thread.c, hierarchicalThreadManager.c)",
  "src/misc/config.h.win|HAVE_ALLOCA_H|no alloca.h on Windows; src/include/misc/alloca.h falls back to _alloca/malloc.h or __builtin_alloca",
  "src/misc/config.h.win|HAVE_FFS|no BSD ffs() on Windows; src/include/misc/intrinsic.h's #elif chain resolves via _MSC_VER/__GNUC__ first"
)
parts <- strsplit(absentRaw, "|", fixed = TRUE)
field <- function(i) vapply(parts, `[[`, character(1L), i)
macroLists <- strsplit(field(2L), ",", fixed = TRUE)
EXPECTED_ABSENT <- data.frame(
  file = rep(field(1L), lengths(macroLists)),
  macro = unlist(macroLists),
  reason = rep(field(3L), lengths(macroLists)),
  stringsAsFactors = FALSE
)

find <- character(0L)
report <- function(msg) find <<- c(find, msg)

# --- (a) version literals ---
descVersion <- read.dcf(p("DESCRIPTION"), fields = "Version")[1L, 1L]
checkVersionMacro <- function(file, macro, lines, wantSuffix) {
  hit <- grep(sprintf("^#define %s ", macro), lines, value = TRUE)
  if (length(hit) != 1L || !grepl(paste0(wantSuffix, "$"), hit)) {
    report(sprintf(
      "%s: %s does not match DESCRIPTION Version %s",
      file,
      macro,
      descVersion
    ))
  }
}
versionFiles <- c(
  "src/config.hpp.win",
  "src/external/config.h.win",
  "src/misc/config.h.win"
)
for (f in versionFiles) {
  lines <- readLines(p(f))
  checkVersionMacro(f, "PACKAGE_STRING", lines, paste0(" ", descVersion, '"'))
  checkVersionMacro(f, "PACKAGE_VERSION", lines, paste0('"', descVersion, '"'))
}

# --- (b) config-header macro nameset ---
headerPairs <- list(
  c("src/config.hpp.win", "src/config.hpp.in"),
  c("src/external/config.h.win", "src/external/config.h.in"),
  c("src/misc/config.h.win", "src/misc/config.h.in"),
  c("src/include/misc/types.h.win", "src/include/misc/types.h.in")
)
for (pair in headerPairs) {
  win <- pair[1L]
  inSet <- macroNames(p(pair[2L]), "undef")
  winSet <- macroNames(p(win), "define")
  missing <- setdiff(inSet, winSet)
  listed <- EXPECTED_ABSENT[EXPECTED_ABSENT$file == win, , drop = FALSE]
  for (m in setdiff(missing, listed$macro)) {
    report(sprintf(
      "%s: %s is in %s but not listed as expected-absent",
      win,
      m,
      pair[2L]
    ))
  }
  for (m in listed$macro[!(listed$macro %in% missing)]) {
    why <- if (m %in% inSet) {
      "now defined in .win"
    } else {
      "dropped from the .in template"
    }
    report(sprintf("%s: expected-absent entry '%s' is stale (%s)", win, m, why))
  }
}

# --- report ---
if (length(find) > 0L) {
  cat(paste0("DRIFT: ", find), sep = "\n")
  quit(status = 1L)
}
nUnknown <- sum(EXPECTED_ABSENT$reason == "UNKNOWN")
cat(sprintf(
  "check-win-drift: OK (%d version literals, %d config headers, %d expected-absent entries, %d UNKNOWN)\n",
  length(versionFiles) * 2L,
  length(headerPairs),
  nrow(EXPECTED_ABSENT),
  nUnknown
))
quit(status = 0L)
