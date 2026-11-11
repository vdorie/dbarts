#!/bin/sh
# Builds one pricing library from a commit and stamps it with its identity.
# The tree is exported with git archive into a fresh source directory, the
# build's patch applied (shipped: the move off; mutant1, mutant3: step 1's
# power mutants; move: none), installed into the library directory, and
# kmp-build.dcf written into the installed package: the build's name, the
# commit, the patch's md5. The harness's libId() reads it beside the md5 of
# the installed shared object.
#
# Usage: k-move-pricing-install.sh BUILD COMMIT REPO SRCDIR LIBDIR
# SRCDIR must not exist (a fresh directory each build); LIBDIR may.
set -eu
build=$1
commit=$2
repo=$3
src=$4
lib=$5
here=$(dirname "$0")
case "$build" in
  move) patchfile= ;;
  shipped) patchfile=$here/k-move-pricing-shipped.patch ;;
  mutant1 | mutant3) patchfile=$here/k-move-pricing-$build.patch ;;
  kstudy) patchfile= ;;
  *)
    echo "unknown build $build" >&2
    exit 2
    ;;
esac
if [ -e "$src" ]; then
  echo "$src exists; give a fresh source directory" >&2
  exit 2
fi
full=$(git -C "$repo" rev-parse "$commit^{commit}")
mkdir -p "$src" "$lib"
git -C "$repo" archive "$full" | tar -x -C "$src"
pmd5=none
if [ -n "$patchfile" ]; then
  patch -d "$src" -p1 < "$patchfile"
  pmd5=$(openssl md5 -r < "$patchfile" | cut -d' ' -f1)
fi
R CMD INSTALL -l "$lib" "$src"
printf 'Build: %s\nCommit: %s\nPatch: %s\n' "$build" "$full" "$pmd5" \
  > "$lib/dbarts/kmp-build.dcf"
cat "$lib/dbarts/kmp-build.dcf"
