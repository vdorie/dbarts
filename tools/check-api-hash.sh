#!/bin/sh
# Fails CI when the C API header's hash (DBARTS_C_API_HASH) changed since
# the newest release tag but the version pair (DBARTS_C_API_MAJOR,
# DBARTS_C_API_MINOR) did not: an ABI change with no minor bump, the
# discipline dec-B111's version-pair guard depends on. Armed only once a
# release tag exists (dbarts.h's own rule is that the hash moves freely
# and the constants stay put pre-release, so there is nothing to compare
# against before the first one). A release tag is v1.* or 1.* - this
# repo's own tags predate the v-prefixed spelling (0.8-7), so both arm
# the check; the newest by version, across either spelling, wins.
#
# POSIX sh; no bashisms. Run from the repository root with full tag
# history available (a shallow checkout sees no tags and this prints the
# skip line instead of silently passing).

set -e

header="inst/include/dbarts/dbarts.h"

# Extracts one #define'd integer/hex literal from a header on stdin.
#   $1: the macro name (DBARTS_C_API_HASH, _MAJOR or _MINOR)
extract_macro() {
  sed -n "s/^#[ 	]*define[ 	]*$1[ 	]*\\([0-9A-Za-z]*\\).*/\\1/p" | head -n 1
}

# Newest tag matching v1.* or 1.*, compared numerically on
# major.minor.patch after stripping an optional leading v: git's own
# --sort=-v:refname groups a v-prefixed tag ahead of every unprefixed one
# rather than interleaving them by version, which is wrong here since
# both spellings name the same tag space.
tag=`git tag -l 'v1.*' '1.*' | awk '
  {
    v = $0
    sub(/^v/, "", v)
    n = split(v, part, /[.-]/)
    if (n < 3) next
    key = sprintf("%010d.%010d.%010d", part[1], part[2], part[3])
    if (key > best) { best = key; besttag = $0 }
  }
  END { print besttag }
'`
if [ -z "$tag" ]; then
  echo "no release tag, skipped"
  exit 0
fi

working_hash=`extract_macro DBARTS_C_API_HASH < "$header"`
working_major=`extract_macro DBARTS_C_API_MAJOR < "$header"`
working_minor=`extract_macro DBARTS_C_API_MINOR < "$header"`

tagged_hash=`git show "$tag:$header" | extract_macro DBARTS_C_API_HASH`
tagged_major=`git show "$tag:$header" | extract_macro DBARTS_C_API_MAJOR`
tagged_minor=`git show "$tag:$header" | extract_macro DBARTS_C_API_MINOR`

if [ -z "$working_hash" ] || [ -z "$tagged_hash" ]; then
  echo "check-api-hash: could not read DBARTS_C_API_HASH from $header or from $tag:$header" >&2
  exit 1
fi

if [ "$working_hash" != "$tagged_hash" ]; then
  if [ "$working_major" = "$tagged_major" ] && [ "$working_minor" = "$tagged_minor" ]; then
    echo "check-api-hash: $header's DBARTS_C_API_HASH changed since $tag ($tagged_hash -> $working_hash) with no DBARTS_C_API_MAJOR/DBARTS_C_API_MINOR bump ($tagged_major.$tagged_minor -> $working_major.$working_minor); an ABI change needs a minor bump" >&2
    exit 1
  fi
  echo "hash changed since $tag, version bumped ($tagged_major.$tagged_minor -> $working_major.$working_minor): ok"
  exit 0
fi

echo "hash unchanged since $tag: ok"
exit 0
