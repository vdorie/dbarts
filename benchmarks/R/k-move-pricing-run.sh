#!/bin/sh
# Runs a pricing job file in order, one job at a time (the laptop smoke; the
# box runs the same file through the k study's queue, which runs jobs in
# parallel). A job runs once every job it waits for is done; a job with a
# missing gate file (relative to KMP_RUN) is recorded skipped; a failed job
# blocks the jobs that wait on it, a dependency written ~id excepted.
# Resumable: a job recorded done or skipped in the status file is not rerun.
#
# Usage: k-move-pricing-run.sh JOBFILE STATUSFILE LOGDIR
# Environment: LIB_SHIPPED, LIB_MOVE (and LIB_MUTANT1, LIB_MUTANT3 for step
# 1's power runs, LIB_KSTUDY for Reading 3's model check), LIB_DEPS, plus
# what the harness reads (KMP_RUN, KMP_QUICK, KMP_SMOKE, DBARTS_BIS_DATA,
# ...). Library paths must be absolute; each is built by
# k-move-pricing-install.sh.
set -u
jobs=$1
status=$2
logs=$3
harness="${HARNESS:-$(dirname "$0")/k-move-pricing.R}"
mkdir -p "$logs"
touch "$status"

state_of() {
  awk -F'\t' -v id="$1" '$1 == id { s = $2 } END { print s }' "$status"
}

tab=$(printf '\t')
while IFS="$tab" read -r id deps mem cell est cmd gate; do
  [ -z "$id" ] && continue
  case "$(state_of "$id")" in
    done | skipped) continue ;;
  esac
  ok=1
  if [ "$deps" != "-" ]; then
    for d in $(echo "$deps" | tr ',' ' '); do
      case "$d" in
        "~"*) ;;
        *)
          s=$(state_of "$d")
          [ "$s" = done ] || [ "$s" = skipped ] || ok=0
          ;;
      esac
    done
  fi
  if [ $ok -eq 0 ]; then
    printf '%s\tblocked\n' "$id" >> "$status"
    continue
  fi
  if [ "$gate" != "-" ]; then
    for g in $(echo "$gate" | tr ',' ' '); do
      [ -e "$KMP_RUN/$g" ] || ok=0
    done
    if [ $ok -eq 0 ]; then
      printf '%s\tskipped\n' "$id" >> "$status"
      continue
    fi
  fi
  build=$(echo "$cmd" | sed 's/^@\([a-z0-9]*\) .*/\1/')
  rest=$(echo "$cmd" | sed 's/^@[a-z0-9]* //')
  case "$build" in
    shipped) lib=$LIB_SHIPPED ;;
    move) lib=$LIB_MOVE ;;
    mutant1) lib=${LIB_MUTANT1:-} ;;
    mutant3) lib=${LIB_MUTANT3:-} ;;
    kstudy) lib=${LIB_KSTUDY:-} ;;
    *) lib= ;;
  esac
  if [ -z "$lib" ]; then
    printf '%s\tfailed\n' "$id" >> "$status"
    echo "no library for build $build ($id)" >&2
    continue
  fi
  start=$(date +%s)
  # shellcheck disable=SC2086
  if R_LIBS="$lib:$LIB_DEPS" KMP_BUILD="$build" Rscript "$harness" $rest \
    > "$logs/$id.log" 2>&1 < /dev/null; then
    s=done
  else
    s=failed
  fi
  printf '%s\t%s\t%s\n' "$id" "$s" "$(($(date +%s) - start))" >> "$status"
  echo "$id $s"
done < "$jobs"
