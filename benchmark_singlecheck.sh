#!/bin/bash
###############################################################################
# benchmark_singlecheck.sh
#
# Compare multiple versions/copies of the SingleCheck pipeline (each living in
# its own folder, e.g. SingleCheckNoChanges/, SingleCheckChangeA/, ...) by:
#   (a) VERIFYING the final output is consistent across versions, and
#   (b) TIMING each version so you know which is faster.
#
# Each run is fully isolated: the input file is SYMLINKED (never copied -- safe
# for 300 GB-TB files) into a private working directory, so versions never
# overwrite each other's outputs and never interfere with timing.
#
# USAGE
#   benchmark_singlecheck.sh -i <input> [-I <mate2>] \
#       (-V dir1,dir2,... | -g 'GLOB') \
#       [-r REPS] [-b BASELINE] [-T TOL] [-o OUTDIR] [--warmup] [--srun] \
#       [-- <extra SingleCheck options, e.g. -w 10000000 -t 8 -r ref.fa>]
#
# OPTIONS
#   -i, --input   FILE   primary input: .bam | .cram | _1.fastq.gz | .fastq.gz  (required)
#   -I, --input2  FILE   second mate for paired-end FASTQ (_2.fastq.gz)
#   -V, --versions LIST  comma-separated list of version directories
#                        (each must contain an executable `SingleCheck`)
#   -g, --glob   PATTERN glob to auto-discover version dirs, e.g. '/home/me/SingleCheck*'
#   -r, --reps   N       repetitions per version (default 1; >1 gives stable timings)
#   -b, --baseline NAME  version (folder basename) used as the reference for
#                        consistency + speedup (default: the first version)
#   -T, --tol    VAL     relative tolerance for numeric output fields (default 1e-6)
#   -o, --outdir DIR     results directory (default ./singlecheck_bench_<timestamp>)
#       --warmup         run one unmeasured warmup per version (warms page cache)
#       --srun           prefix each run with `srun` (use inside a SLURM allocation)
#   -h, --help           show this help
#   --                   everything after is passed verbatim to every SingleCheck
#
# EXAMPLE
#   salloc -c8 --mem 80G                 # get an interactive allocation, then:
#   ./benchmark_singlecheck.sh \
#       -i /data/HG002.bam -r 3 \
#       -g '/home/tretorn/SingleCheck*' \
#       -b SingleCheckNoChanges -T 1e-4 \
#       -- -w 10000000 -t 8
###############################################################################

set -uo pipefail

usage() { sed -n '2,/^###########/p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//;s/^#//'; }

# ---- defaults ---------------------------------------------------------------
INPUT="" ; INPUT2="" ; VERSIONS_CSV="" ; GLOB="" ; REPS=1 ; BASELINE=""
TOL="1e-6" ; OUTDIR="" ; WARMUP=0 ; USE_SRUN=0
EXTRA=()

# ---- parse wrapper options --------------------------------------------------
while [ "$#" -gt 0 ]; do
  case "$1" in
    -i|--input)    INPUT="$2"; shift 2 ;;
    -I|--input2)   INPUT2="$2"; shift 2 ;;
    -V|--versions) VERSIONS_CSV="$2"; shift 2 ;;
    -g|--glob)     GLOB="$2"; shift 2 ;;
    -r|--reps)     REPS="$2"; shift 2 ;;
    -b|--baseline) BASELINE="$2"; shift 2 ;;
    -T|--tol)      TOL="$2"; shift 2 ;;
    -o|--outdir)   OUTDIR="$2"; shift 2 ;;
    --warmup)      WARMUP=1; shift ;;
    --srun)        USE_SRUN=1; shift ;;
    -h|--help)     usage; exit 0 ;;
    --)            shift; EXTRA=("$@"); break ;;
    -*)            printf 'Unknown option: %s\n\n' "$1" >&2; usage >&2; exit 1 ;;
    *)             printf 'Unexpected argument: %s (did you forget "--"?)\n\n' "$1" >&2; usage >&2; exit 1 ;;
  esac
done

# ---- validate ---------------------------------------------------------------
[ -n "$INPUT" ]    || { echo "Error: --input is required" >&2; usage >&2; exit 1; }
[ -f "$INPUT" ]    || { echo "Error: input not found: $INPUT" >&2; exit 1; }
[ -z "$INPUT2" ] || [ -f "$INPUT2" ] || { echo "Error: mate not found: $INPUT2" >&2; exit 1; }

# Resolve the list of version directories
declare -a VERSIONS=()
if [ -n "$VERSIONS_CSV" ]; then
  IFS=',' read -r -a VERSIONS <<< "$VERSIONS_CSV"
elif [ -n "$GLOB" ]; then
  for d in $GLOB; do [ -d "$d" ] && [ -x "$d/SingleCheck" ] && VERSIONS+=("$d"); done
else
  echo "Error: specify versions with -V or -g" >&2; usage >&2; exit 1
fi
[ "${#VERSIONS[@]}" -ge 1 ] || { echo "Error: no version directories resolved" >&2; exit 1; }

for d in "${VERSIONS[@]}"; do
  [ -x "$d/SingleCheck" ] || { echo "Error: $d/SingleCheck missing or not executable" >&2; exit 1; }
done

# Baseline defaults to the first version (compare by folder basename)
declare -A VDIR=()
declare -a VNAMES=()
for d in "${VERSIONS[@]}"; do
  n="$(basename "$(cd "$d" && pwd -P)")"
  VDIR["$n"]="$(cd "$d" && pwd -P)"
  VNAMES+=("$n")
done
[ -n "$BASELINE" ] || BASELINE="${VNAMES[0]}"
[ -n "${VDIR[$BASELINE]:-}" ] || { echo "Error: baseline '$BASELINE' is not among the versions" >&2; exit 1; }

OUTDIR="${OUTDIR:-$PWD/singlecheck_bench_$(date +%Y%m%d_%H%M%S)}"
mkdir -p "$OUTDIR"
SUMMARY="$OUTDIR/summary.tsv"
TIMINGS="$OUTDIR/timings.tsv"
printf "version\trep\twall_s\texit_code\tmax_rss_kb\toutput\n" > "$TIMINGS"

# Optional external `time` for max-RSS
TIMEBIN=""
for cand in /usr/bin/time /bin/time; do [ -x "$cand" ] && TIMEBIN="$cand" && break; done

ABS_INPUT="$(readlink -f "$INPUT")"
[ -n "$INPUT2" ] && ABS_INPUT2="$(readlink -f "$INPUT2")"
BASE_IN="$(basename "$INPUT")"
[ -n "$INPUT2" ] && BASE_IN2="$(basename "$INPUT2")"

echo "Benchmark settings"
echo "  input      : $INPUT ${INPUT2:+(+ $INPUT2)}"
echo "  versions   : ${VNAMES[*]}"
echo "  baseline   : $BASELINE"
echo "  reps       : $REPS   warmup: $WARMUP   srun: $USE_SRUN"
echo "  tolerance  : $TOL"
echo "  extra args : ${EXTRA[*]:-<none>}"
echo "  results    : $OUTDIR"
echo

# ---- helper: symlink the input(s) (+index) into a run dir -------------------
link_inputs() {
  local rundir="$1"
  ln -sf "$ABS_INPUT" "$rundir/$BASE_IN"
  [ -f "${ABS_INPUT}.bai" ]  && ln -sf "${ABS_INPUT}.bai"  "$rundir/${BASE_IN}.bai"
  [ -f "${ABS_INPUT}.crai" ] && ln -sf "${ABS_INPUT}.crai" "$rundir/${BASE_IN}.crai"
  if [ -n "$INPUT2" ]; then ln -sf "$ABS_INPUT2" "$rundir/$BASE_IN2"; fi
}

# ---- helper: run one version+rep, record timing, echo the output file -------
run_once() {
  local vname="$1" rep="$2" rundir="$3"
  rm -rf "$rundir"; mkdir -p "$rundir"
  link_inputs "$rundir"

  local -a cmd=("${VDIR[$vname]}/SingleCheck")
  [ "${#EXTRA[@]}" -gt 0 ] && cmd+=("${EXTRA[@]}")
  cmd+=("$rundir/$BASE_IN")
  [ -n "$INPUT2" ] && cmd+=("$rundir/$BASE_IN2")
  [ "$USE_SRUN" -eq 1 ] && cmd=(srun --exclusive -N1 -n1 "${cmd[@]}")

  local start end wall rc rss
  start=$(date +%s.%N)
  if [ -n "$TIMEBIN" ]; then
    "$TIMEBIN" -v -o "$rundir/time.txt" "${cmd[@]}" > "$rundir/run.log" 2>&1
    rc=$?
  else
    "${cmd[@]}" > "$rundir/run.log" 2>&1
    rc=$?
  fi
  end=$(date +%s.%N)
  wall=$(awk -v a="$start" -v b="$end" 'BEGIN{printf "%.3f", b-a}')
  rss="NA"
  [ -f "$rundir/time.txt" ] && rss=$(awk '/Maximum resident set size/{print $NF}' "$rundir/time.txt")

  local out=""
  out=$(ls "$rundir"/*.SingleCheck.txt 2>/dev/null | head -n1)
  printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$vname" "$rep" "$wall" "$rc" "${rss:-NA}" "${out:-NONE}" >> "$TIMINGS"
  RUN_RC=$rc; RUN_WALL=$wall; RUN_OUT="$out"
}

# ---- helper: numeric-aware comparison of two single-line output files -------
# returns 0 if consistent within $3 tolerance, 1 otherwise; prints field diffs
compare_outputs() {
  awk -v tol="$3" '
    function isnum(x){ return (x ~ /^[+-]?([0-9]+\.?[0-9]*|\.[0-9]+)([eE][+-]?[0-9]+)?$/) }
    NR==FNR { n=split($0,a,"\t"); next }
            { m=split($0,b,"\t") }
    END{
      fail=0
      if (n!=m){ printf("    field-count differs: %d vs %d\n", n, m); fail=1 }
      cols=(n<m?n:m)
      for(i=1;i<=cols;i++){
        if(a[i]==b[i]) continue
        if(isnum(a[i]) && isnum(b[i])){
          d=a[i]-b[i]; if(d<0)d=-d
          den=(a[i]<0?-a[i]:a[i]); if(den<1e-12)den=1
          rel=d/den
          if(rel>tol+0){ printf("    col %2d: %s  vs  %s   (reldiff %.3g)\n", i, a[i], b[i], rel); fail=1 }
        } else {
          printf("    col %2d: \"%s\"  vs  \"%s\"\n", i, a[i], b[i]); fail=1
        }
      }
      exit fail
    }' "$1" "$2"
}

# ---- execute all runs -------------------------------------------------------
declare -A WALLS=()        # space-separated wall times per version
declare -A FIRSTOUT=()     # path to rep1 output per version
declare -A STATUS=()       # OK | FAILED | NONDETERMINISTIC
declare -A DETERMIN=()     # yes | no

for vname in "${VNAMES[@]}"; do
  echo ">>> Version: $vname"
  STATUS["$vname"]="OK"; DETERMIN["$vname"]="yes"; WALLS["$vname"]=""

  if [ "$WARMUP" -eq 1 ]; then
    echo "    warmup..."
    run_once "$vname" "warmup" "$OUTDIR/$vname/warmup"
  fi

  for rep in $(seq 1 "$REPS"); do
    run_once "$vname" "$rep" "$OUTDIR/$vname/rep$rep"
    echo "    rep $rep: wall=${RUN_WALL}s exit=${RUN_RC}"
    if [ "$RUN_RC" -ne 0 ] || [ -z "$RUN_OUT" ]; then
      STATUS["$vname"]="FAILED"
      echo "    !! run failed (exit $RUN_RC) -- see $OUTDIR/$vname/rep$rep/run.log" >&2
      continue
    fi
    WALLS["$vname"]+=" $RUN_WALL"
    if [ "$rep" -eq 1 ]; then
      FIRSTOUT["$vname"]="$RUN_OUT"
    else
      # determinism: rep must match rep1 exactly (seed is fixed)
      if ! compare_outputs "${FIRSTOUT[$vname]}" "$RUN_OUT" 0 >/dev/null 2>&1; then
        DETERMIN["$vname"]="no"
      fi
    fi
  done
  # Save a copy of this version's canonical output
  if [ -n "${FIRSTOUT[$vname]:-}" ]; then
    cp "${FIRSTOUT[$vname]}" "$OUTDIR/$vname.SingleCheck.txt"
  fi
  echo
done

# ---- consistency vs baseline + summary --------------------------------------
echo "==================== RESULTS ===================="
base_out="${FIRSTOUT[$BASELINE]:-}"
base_mean=$(awk -v s="${WALLS[$BASELINE]:-}" 'BEGIN{n=split(s,t," ");for(i=1;i<=n;i++)x+=t[i];if(n)printf "%.3f",x/n;else print "NA"}')

printf "%-28s %5s %12s %12s %9s %12s %14s %s\n" \
  "VERSION" "REPS" "WALL_min(s)" "WALL_mean(s)" "SPEEDUP" "MAXRSS(KB)" "DETERMINISTIC" "VS_BASELINE"
printf "%-28s %5s %12s %12s %9s %12s %14s %s\n" \
  "$BASELINE" "$REPS" \
  "$(awk -v s="${WALLS[$BASELINE]:-}" 'BEGIN{n=split(s,t," ");m=1e99;for(i=1;i<=n;i++)if(t[i]<m)m=t[i];if(n)printf "%.3f",m;else print "NA"}')" \
  "${base_mean}" "1.00x" \
  "$(awk -F'\t' -v v="$BASELINE" '$1==v && $5!="NA"{print $5; exit}' "$TIMINGS" 2>/dev/null || echo NA)" \
  "${DETERMIN[$BASELINE]}" "baseline (${STATUS[$BASELINE]})"

printf "%-28s %5s %12s %12s %9s %12s %14s %s\n" "$BASELINE" "$REPS" "" "" "" "" "" "" >/dev/null

{
  printf "version\twall_min_s\twall_mean_s\tspeedup_vs_baseline\tmax_rss_kb\tdeterministic\tconsistency\tstatus\n"
} > "$SUMMARY"

for vname in "${VNAMES[@]}"; do
  wmin=$(awk -v s="${WALLS[$vname]:-}" 'BEGIN{n=split(s,t," ");m=1e99;for(i=1;i<=n;i++)if(t[i]<m)m=t[i];if(n)printf "%.3f",m;else print "NA"}')
  wmean=$(awk -v s="${WALLS[$vname]:-}" 'BEGIN{n=split(s,t," ");for(i=1;i<=n;i++)x+=t[i];if(n)printf "%.3f",x/n;else print "NA"}')
  rss=$(awk -F'\t' -v v="$vname" '$1==v && $5!="NA"{print $5; exit}' "$TIMINGS" 2>/dev/null); rss="${rss:-NA}"
  speed="NA"
  if [ "$wmean" != "NA" ] && [ "$base_mean" != "NA" ]; then
    speed=$(awk -v b="$base_mean" -v v="$wmean" 'BEGIN{if(v>0)printf "%.2fx", b/v; else print "NA"}')
  fi

  consistency="-"
  if [ "$vname" = "$BASELINE" ]; then
    consistency="baseline"
  elif [ "${STATUS[$vname]}" = "FAILED" ]; then
    consistency="N/A(failed)"
  elif [ -n "$base_out" ] && [ -n "${FIRSTOUT[$vname]:-}" ]; then
    diffout=$(compare_outputs "$base_out" "${FIRSTOUT[$vname]}" "$TOL"); drc=$?
    if [ "$drc" -eq 0 ]; then
      consistency="CONSISTENT"
    else
      ndiff=$(printf "%s\n" "$diffout" | grep -c 'col ')
      consistency="DIFFERS($ndiff)"
    fi
  fi

  [ "$vname" = "$BASELINE" ] && continue   # baseline already printed above
  printf "%-28s %5s %12s %12s %9s %12s %14s %s\n" \
    "$vname" "$REPS" "$wmin" "$wmean" "$speed" "$rss" "${DETERMIN[$vname]}" "$consistency (${STATUS[$vname]})"

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$vname" "$wmin" "$wmean" "$speed" "$rss" "${DETERMIN[$vname]}" "$consistency" "${STATUS[$vname]}" >> "$SUMMARY"

  # show field-level diffs when a version differs from baseline
  if [[ "$consistency" == DIFFERS* ]]; then
    echo "    --- $vname differs from baseline $BASELINE (tol=$TOL) ---"
    compare_outputs "$base_out" "${FIRSTOUT[$vname]}" "$TOL"
  fi
done

# write baseline row to summary too
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
  "$BASELINE" \
  "$(awk -v s="${WALLS[$BASELINE]:-}" 'BEGIN{n=split(s,t," ");m=1e99;for(i=1;i<=n;i++)if(t[i]<m)m=t[i];if(n)printf "%.3f",m;else print "NA"}')" \
  "$base_mean" "1.00x" \
  "$(awk -F'\t' -v v="$BASELINE" '$1==v && $5!="NA"{print $5; exit}' "$TIMINGS")" \
  "${DETERMIN[$BASELINE]}" "baseline" "${STATUS[$BASELINE]}" >> "$SUMMARY"

echo "================================================="
echo "Per-run timings : $TIMINGS"
echo "Summary table   : $SUMMARY"
echo "Saved outputs   : $OUTDIR/<version>.SingleCheck.txt"
echo
echo "Notes:"
echo " * WALL_min is the fairest single-number comparison (page cache warms across reps)."
echo " * 'DIFFERS' is expected between Picard- and samtools-downsampling versions:"
echo "   the sampled read set changes, so coverage metrics move. Inspect the field"
echo "   diffs above and raise --tol if the differences are within acceptable noise."
echo " * Run this inside a SLURM allocation (salloc/srun) for representative timings."
