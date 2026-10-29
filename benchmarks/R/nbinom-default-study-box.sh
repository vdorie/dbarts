#!/bin/sh
# The negative binomial default-prior study's main run on a run machine, as
# docs/plans/nbinom-default-study.md's Box readiness states it. Run from a
# checkout of the study's branch, in this order, before any job:
#
#   nbinom-default-study-box.sh build <library root>
#   nbinom-default-study-box.sh identity <library root> <run directory>
#   nbinom-default-study-box.sh calibrate <library root> <run directory>
#
# then the plan (look 1's queue) and one worker a core on each queue:
#
#   NBS_RUN=<run directory> NBS_BUILD=today R_LIBS=<library root>/today \
#     Rscript nbinom-default-study.R plan
#   nbinom-default-study-box.sh worker <library root> <run directory> <queue>
#
# build installs three libraries under the library root, each from a fresh
# tree of main at 8a78507c ("today", the pinned base): today/ as is, move/
# with the move patch, kprior/ with the kprior patch on top. identity checks
# there that the kprior build draws as the move build with its prior unset
# and as today's with the move off, and stops if not. calibrate times three
# fixed runs against the laptop's and rescales the sizing table. The run
# directory holds data/, copied from the laptop's export with its truths,
# and every output.
set -eu
base=8a78507c
here=$(cd "$(dirname "$0")" && pwd)
repo=$(git -C "$here" rev-parse --show-toplevel)
harness=$here/nbinom-default-study.R

build() {
  root=$1
  mkdir -p "$root/today" "$root/move" "$root/kprior"
  tree=$(mktemp -d)
  git -C "$repo" archive "$base" | tar -x -C "$tree"
  R CMD INSTALL --preclean -l "$root/today" "$tree"
  patch -d "$tree" -p1 < "$here/nbinom-default-study-move.patch"
  R CMD INSTALL --preclean -l "$root/move" "$tree"
  patch -d "$tree" -p1 < "$here/nbinom-default-study-kprior.patch"
  R CMD INSTALL --preclean -l "$root/kprior" "$tree"
  echo "built; the source tree $tree is no longer needed"
}

identity() {
  root=$1
  run=$2
  for seed in 1 7; do
    NBS_RUN=$run NBS_BUILD=today R_LIBS=$root/today \
      Rscript "$harness" identity $seed
    NBS_RUN=$run NBS_BUILD=move R_LIBS=$root/move \
      Rscript "$harness" identity $seed
    NBS_RUN=$run NBS_BUILD=kprior R_LIBS=$root/kprior \
      Rscript "$harness" identity $seed
    NBS_RUN=$run NBS_BUILD=kprior NBS_KPRIOR=floor:1e-6 R_LIBS=$root/kprior \
      Rscript "$harness" identity $seed
    NBS_RUN=$run NBS_BUILD=kprior NBS_MOVE=off R_LIBS=$root/kprior \
      Rscript "$harness" identity $seed
    NBS_RUN=$run NBS_BUILD=kprior NBS_KPRIOR=floor:1e-6 NBS_MOVE=off \
      R_LIBS=$root/kprior Rscript "$harness" identity $seed
  done
  summary=$(NBS_RUN=$run NBS_BUILD=today R_LIBS=$root/today \
    Rscript "$harness" summarize identity)
  echo "$summary"
  if echo "$summary" | grep -q DIFFERENT; then
    echo "the identity fails here: run no job" >&2
    exit 1
  fi
  if [ "$(echo "$summary" | grep -c identical)" -ne 8 ]; then
    echo "the identity is incomplete: run no job" >&2
    exit 1
  fi
}

calibrate() {
  root=$1
  run=$2
  for cell in roaches-w3 crimtab-w3 mite; do
    NBS_RUN=$run NBS_BUILD=kprior R_LIBS=$root/kprior \
      Rscript "$harness" calibrate $cell &
  done
  wait
  NBS_RUN=$run NBS_BUILD=kprior R_LIBS=$root/kprior \
    Rscript "$harness" summarize calibration
}

# Each queue line "<build> <arguments>" is claimed by an atomic mkdir, so
# several workers split one queue without running a job twice. A build of
# kprior, kprior-off (the move off) and either with @<prior> runs the
# kprior library with NBS_MOVE=off and NBS_KPRIOR=<prior> as named.
worker() {
  root=$1
  run=$2
  queue=$3
  mkdir -p "$run/logs" "$run/claims"
  while read -r b rest; do
    name=$(echo "$b $rest" | tr ' /:' '---')
    mkdir "$run/claims/$name" 2>/dev/null || continue
    build=$b
    kprior=
    move=
    case "$b" in
      today) lib=$root/today ;;
      move) lib=$root/move ;;
      kprior*)
        lib=$root/kprior
        build=kprior
        case "$b" in *@*) kprior=${b#*@} ;; esac
        case "$b" in kprior-off*) move=off ;; esac
        ;;
      *)
        echo "bad build $b"
        continue
        ;;
    esac
    start=$(date +%s)
    status=0
    NBS_RUN=$run NBS_BUILD=$build NBS_KPRIOR=$kprior NBS_MOVE=$move \
      R_LIBS=$lib Rscript "$harness" $rest > "$run/logs/$name.log" 2>&1 ||
      status=$?
    echo "$(date +%H:%M:%S) $name exit $status $(($(date +%s) - start)) s"
  done < "$queue"
  echo "worker done"
}

command=$1
shift
case "$command" in
  build) build "$@" ;;
  identity) identity "$@" ;;
  calibrate) calibrate "$@" ;;
  worker) worker "$@" ;;
  *)
    echo "unknown command $command" >&2
    exit 2
    ;;
esac
