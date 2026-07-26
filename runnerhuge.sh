#!/bin/bash
###############################################################################
# runnerhuge.sh -- single production run of SingleCheck on the 300 GB sample.
#
# Not a benchmark: runner.sh runs the pipeline 6 times (2 versions x 3 reps),
# which at this size is a day of compute to measure a few percent. This script
# runs it ONCE, with the settings that matter at 300 GB, after checking
# everything that could waste those hours.
#
# The run benchmarks itself: SingleCheck times every phase and writes a
# .SingleCheck.timings.tsv next to the result, which costs nothing. --stages
# additionally submits the per-chunk profiler as a second job (~one more
# pipeline run at this size).
#
# Usage:
#   ./runnerhuge.sh [/path/to/Wang300.bam]      # submit
#   ./runnerhuge.sh --stages                    # + per-chunk profile job
#   ./runnerhuge.sh --dry-run                   # show the sbatch lines only
###############################################################################
set -uo pipefail

# Resolve from the script's own location, so this works in any checkout:
# on the cluster REPO=$CURSO/SingleCheck and CURSO=.../curso385 as before.
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd -P)"
CURSO="$(dirname "$REPO")"

CPUS=32
MEM=128G
TIME=48:00:00          # lower it if your QoS rejects it (sbatch will say so)

DRYRUN=""
BAM=""
DO_STAGES=0            # --stages: also submit the per-chunk profiler
for a in "$@"; do
    case "$a" in
        --dry-run) DRYRUN="--dry-run" ;;
        --stages)  DO_STAGES=1 ;;
        -h|--help) sed -n '2,/^###########/p' "$0" | sed 's/^# \{0,1\}//;s/^#//'; exit 0 ;;
        *)         BAM="$a" ;;
    esac
done
BAM="${BAM:-$CURSO/Wang300.bam}"

fail() { printf '\nSTOP: %s\n' "$1" >&2; exit 1; }
cd "$REPO" || fail "cannot cd to the repository: $REPO"
[ -x "$REPO/SingleCheck" ] || fail "no executable SingleCheck in $REPO"
[ -x "$REPO/submit" ]      || fail "no executable submit in $REPO"

echo "=============================================================="
echo " SingleCheck -- huge-input run"
echo " input : $BAM"
echo " repo  : $REPO"
echo "=============================================================="

# ---------------------------------------------------------------------------
# 1. Per-cluster environment. Resolve the tools exactly the way the pipeline
#    will (by sourcing the same file), then verify them HERE -- a broken
#    $MOSDEPTH must not be discovered 30 seconds into a 24 h allocation.
# ---------------------------------------------------------------------------
export SINGLECHECK_ENV=$REPO/singlecheck_env.sh
[ -f "$SINGLECHECK_ENV" ] || fail "no env file at $SINGLECHECK_ENV"

MOSDEPTH_RESOLVED=$( set +u; source "$SINGLECHECK_ENV" >/dev/null 2>&1; echo "${MOSDEPTH:-mosdepth}" )
if [ -d "$MOSDEPTH_RESOLVED" ] || [ ! -f "$MOSDEPTH_RESOLVED" ] || [ ! -x "$MOSDEPTH_RESOLVED" ]; then
    if command -v mosdepth >/dev/null 2>&1; then
        MOSDEPTH_RESOLVED=$(command -v mosdepth)
    else
        printf 'mosdepth is not usable: %s\n' "$MOSDEPTH_RESOLVED" >&2
        printf 'find the real binary with:\n    find %s/mosdepth -maxdepth 3 -name "mosdepth*" -type f\n' "$CURSO" >&2
        fail "fix MOSDEPTH in $SINGLECHECK_ENV"
    fi
fi
export MOSDEPTH="$MOSDEPTH_RESOLVED"
"$MOSDEPTH" --version >/dev/null 2>&1 || fail "$MOSDEPTH exists but does not run"
echo "[ok]   mosdepth : $MOSDEPTH ($("$MOSDEPTH" --version 2>&1 | head -n1))"

command -v samtools >/dev/null 2>&1 || module load samtools 2>/dev/null
command -v samtools >/dev/null 2>&1 || fail "samtools not on \$PATH (module load samtools)"
echo "[ok]   samtools : $(command -v samtools) ($(samtools --version 2>&1 | head -n1))"

# ---------------------------------------------------------------------------
# 2. The input and its index
# ---------------------------------------------------------------------------
[ -r "$BAM" ] || fail "input not readable: $BAM"
[ -s "$BAM" ] || fail "input is empty: $BAM"
echo "[ok]   size     : $(du -Lh "$BAM" 2>/dev/null | awk '{print $1}')"

if [ -f "${BAM}.bai" ]; then
    echo "[ok]   index    : ${BAM}.bai"
    [ "${BAM}.bai" -ot "$BAM" ] && echo "[WARN] the index is older than the BAM -- re-run samtools index"
elif [ -f "${BAM}.csi" ]; then
    fail "only a .csi index exists. SingleCheck looks for .bai (ALNIDX in the
      script). Either build one:   samtools index -@ $CPUS $BAM
      or tell me and I will make ALNIDX accept .csi (2-line change)."
else
    fail "no index next to $BAM -- run: samtools index -@ $CPUS $BAM"
fi

samtools quickcheck "$BAM" 2>/dev/null || fail "samtools quickcheck failed: truncated BAM?"
echo "[ok]   quickcheck passed"

# ---------------------------------------------------------------------------
# 3. Why MetaPhyler is skipped, quantified. idxstats reads only the index, so
#    this is instant even on 300 GB. The branch would write EVERY unmapped read
#    as uncompressed FASTA and re-scan the whole file to do it -- for a column
#    that is NA anyway while the marker database is missing.
# ---------------------------------------------------------------------------
read unmapped_reads fasta_gb < <(samtools idxstats "$BAM" 2>/dev/null | \
    awk '{u+=$4} END{printf "%d %.1f\n", u, u*250/1024/1024/1024}')
echo "[info] unmapped reads: ${unmapped_reads:-?}  (~${fasta_gb:-?} GB of FASTA if MetaPhyler ran)"

# ---------------------------------------------------------------------------
# 3b. The optimal flag set for a huge input. Together these take the pipeline
#     from THREE full passes over the input to ONE (only the downsampling).
#     Each is individually reversible; see OPTIMIZATION_REPORT.md.
# ---------------------------------------------------------------------------
export SINGLECHECK_SKIP_METAPHYLER=1   # no unmapped scan at all (column = NA)
export SINGLECHECK_INDEX_DEPTH=1       # read count from the index, not a full pass
export SINGLECHECK_FAST_UNMAPPED=1     # if MetaPhyler is ever re-enabled: '*' block only
echo "[ok]   MetaPhyler branch disabled     -> one less full pass"
echo "[ok]   depth from the index           -> one less full pass (ESTIMATE, see below)"
echo "[ok]   uncompressed transient BAM + piped shifted track (automatic)"

# ---------------------------------------------------------------------------
# 4. I/O environment: striping on the input, and node-local scratch capacity.
#    NOTE: $SLURM_TMPDIR only exists inside the job, so the space check here is
#    indicative only -- the pipeline decides at run time and prints which
#    location it chose in its "> HPC settings" block.
# ---------------------------------------------------------------------------
if command -v lfs >/dev/null 2>&1; then
    stripe=$(lfs getstripe -c "$BAM" 2>/dev/null | tail -n1)
    if [ "${stripe:-0}" = "1" ]; then
        echo "[WARN] Lustre stripe count is 1: every full-file pass is limited to one OST."
        echo "       Cannot be changed in place; a re-striped copy would read faster:"
        echo "         lfs setstripe -c 8 -S 4M <newdir> && cp $BAM <newdir>/"
    else
        echo "[ok]   Lustre stripe count: ${stripe:-unknown}"
    fi
fi
if [ -d /scratch ]; then
    echo "[info] login-node /scratch free: $(df -Ph /scratch 2>/dev/null | awk 'NR==2{print $4}')"
    echo "       (the job gets its own /scratch/\$SLURM_JOB_ID; the pipeline falls"
    echo "        back to Lustre automatically if it has less than 20 GB free)"
fi

# ---------------------------------------------------------------------------
# 5. Submit. -t is not passed: SingleCheck takes its thread count from
#    $SLURM_CPUS_PER_TASK. The variables are named explicitly on --export so
#    the run does not depend on the site's SBATCH_EXPORT default.
# ---------------------------------------------------------------------------
echo "--------------------------------------------------------------"
echo " cpus=$CPUS  mem=$MEM  time=$TIME"
echo "--------------------------------------------------------------"

EXPORTS="ALL,SINGLECHECK_ENV=$SINGLECHECK_ENV,MOSDEPTH=$MOSDEPTH,SINGLECHECK_SKIP_METAPHYLER=1,SINGLECHECK_INDEX_DEPTH=1,SINGLECHECK_FAST_UNMAPPED=1"

# The run measures itself: SingleCheck times every phase and writes
# <sample>.SingleCheck.timings.tsv next to the result. At this size that is the
# only affordable benchmark -- bench/bench_stages.sh re-runs each chunk, which
# means paying for the full-file passes several times over.
./submit $DRYRUN \
    --cpus "$CPUS" --mem "$MEM" --time "$TIME" \
    --name sc_huge \
    --sbatch "--export=$EXPORTS" \
    SingleCheck "$BAM"

# --stages: per-chunk detail as a SEPARATE job. One repetition only, no legacy
# variants, intermediates on node-local scratch: roughly the cost of one more
# pipeline run. Do not raise -r here without doing the arithmetic first.
if [ "$DO_STAGES" -eq 1 ]; then
    echo
    echo ">>> per-chunk profile (separate job, -r 1, no legacy variants)"
    ./submit $DRYRUN \
        --cpus "$CPUS" --mem "$MEM" --time "$TIME" \
        --name sc_huge_stages \
        --sbatch "--export=$EXPORTS" \
        bench/bench_stages.sh \
        -i "$BAM" \
        -r 1 \
        --threads "$CPUS" \
        --no-compare \
        --work-scratch
fi

[ -n "$DRYRUN" ] && exit 0

cat <<EOF

Submitted. In sc_huge-<jobid>.out check the "> HPC settings" block first:

  intermediates      : /scratch/<jobid>/...    <- must NOT say "shared filesystem"
  MetaPhyler branch  : SKIPPED
  sequencing depth   : from the index (no counting pass; ESTIMATE)
  autocorrelation    : sort-free aggregation, shifted track piped (no temp file)
  downsampled BAM    : uncompressed (transient)
  threads            : $CPUS (samtools/sort/aligner), $CPUS (mosdepth)

The run benchmarks itself. At the end of the log:

  > Timing summary (wall clock)
    PHASE                             SECONDS    SHARE
    downsample                          ...      ...%
    mosdepth                            ...      ...%
    ...

Then:
  squeue -u \$USER
  sacct -j <jobid> --format=JobID,Elapsed,MaxRSS,State   # MaxRSS sizes the next --mem

Results:
  ${BAM%.bam}.SingleCheck.txt            13 tab-separated columns (contaminants = NA)
  ${BAM%.bam}.SingleCheck.timings.tsv    per-phase wall clock, machine readable
$([ "$DO_STAGES" -eq 1 ] && echo "  bench/stages_<host>_<ts>/report.md     per-chunk profile (second job)")

NOTE on SINGLECHECK_INDEX_DEPTH=1: columns 2-3 ("Sequenced bases", "Analysis
depth") become an estimate -- the index gives every alignment record, and the
secondary/supplementary share is measured on the first million reads. Drop that
one variable if you need those two columns exact; it costs one full pass.
bench/bench_stages.sh reports the exact-vs-estimate error on your own data.
EOF
