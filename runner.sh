#!/bin/bash
###############################################################################
# runner.sh -- A/B benchmark: current SingleCheck vs the commit before the HPC
#              optimizations, in ONE job so both arms share node and filesystem.
#
# Usage:
#   ./runner.sh              # exact comparison (default)
#   ./runner.sh --fast       # also enable the optimizations that ESTIMATE
#   ./runner.sh --dry-run
#
# DEFAULT MODE enables every optimization that provably cannot change the
# output, so the final `CONSISTENT` verdict stays meaningful at -T 1e-6:
#   * SINGLECHECK_SKIP_METAPHYLER=1 -- the baseline ignores the variable and runs
#     the full unmapped scan, the current version skips it; both end with
#     class=NA (marker DB missing), so the output line is unaffected while the
#     time saved by dropping that pass IS measured.
#   * uncompressed transient BAM, piped shifted track, sort-free aggregation,
#     streaming mapping stats, node-local intermediates -- all in the code, all
#     bit-identical, all measured automatically.
#
# --fast additionally sets SINGLECHECK_INDEX_DEPTH=1 (read count from the index)
# and SINGLECHECK_FAST_UNMAPPED=1. Those are an ESTIMATE and a read subset, so
# columns 2-3 will move: the tolerance is raised to 1e-2 and a DIFFERS there is
# information, not a bug. bench/bench_stages.sh reports the exact error.
###############################################################################
set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd -P)"
CURSO="$(dirname "$REPO")"
BASELINE_COMMIT=35ec840          # last commit before the round-2 optimizations

CPUS=32
MEM=128G
REPS=3
TOL=1e-6
FAST=0
DRYRUN=""
BAM=""
DO_AB=1        # end-to-end A/B  (benchmark_singlecheck.sh)
DO_STAGES=1    # per-chunk profile (bench/bench_stages.sh)

for a in "$@"; do
    case "$a" in
        --fast)        FAST=1 ;;
        --ab-only)     DO_STAGES=0 ;;
        --stages-only) DO_AB=0 ;;
        --dry-run)     DRYRUN="--dry-run" ;;
        -h|--help)     sed -n '2,/^###########/p' "$0" | sed 's/^# \{0,1\}//;s/^#//'; exit 0 ;;
        *)             BAM="$a" ;;
    esac
done
BAM="${BAM:-$CURSO/Wang5.bam}"

fail() { printf '\nSTOP: %s\n' "$1" >&2; exit 1; }
cd "$REPO" || fail "cannot cd to $REPO"
[ -r "$BAM" ] || fail "input not readable: $BAM"

# The CURRENT arm must understand the flags exported below and carry the phase
# timing; unknown variables are silently ignored, which would look like a
# successful but meaningless comparison. (The baseline arm is expected to lack
# them -- that is the point of the A/B.)
missing_feats=""
for feat in print_timings SINGLECHECK_SKIP_METAPHYLER; do
    grep -q "$feat" "$REPO/SingleCheck" || missing_feats="$missing_feats $feat"
done
[ -n "$missing_feats" ] && fail "this checkout of SingleCheck does not support:$missing_feats -- git pull in $REPO first"

# ---------------------------------------------------------------------------
# 1. baseline worktree -- on Lustre, not /tmp: the compute node must see it
# ---------------------------------------------------------------------------
if [ ! -x "$CURSO/SingleCheckBaseline/SingleCheck" ]; then
    git worktree add "$CURSO/SingleCheckBaseline" "$BASELINE_COMMIT" || exit 1
fi

# ---------------------------------------------------------------------------
# 2. Environment. benchmark_singlecheck.sh has no --env option: it inherits the
#    exported variable and every SingleCheck it launches sources it. Both
#    versions source the SAME file, so the mosdepth fix applies to the baseline
#    worktree too. Verify the tools HERE, before queuing anything.
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
echo "[ok] mosdepth : $MOSDEPTH ($("$MOSDEPTH" --version 2>&1 | head -n1))"

# ---------------------------------------------------------------------------
# 3. flag set (see the header for why the default is what it is)
# ---------------------------------------------------------------------------
EXPORTS="ALL,SINGLECHECK_ENV=$SINGLECHECK_ENV,MOSDEPTH=$MOSDEPTH,SINGLECHECK_SKIP_METAPHYLER=1"
if [ "$FAST" -eq 1 ]; then
    EXPORTS="$EXPORTS,SINGLECHECK_INDEX_DEPTH=1,SINGLECHECK_FAST_UNMAPPED=1"
    TOL=1e-2
    echo "[fast]  index-based depth + '*'-only unmapped: columns 2-3 become estimates, tolerance $TOL"
else
    echo "[exact] only bit-identical optimizations are enabled; tolerance $TOL"
fi

# ---------------------------------------------------------------------------
# 4. Submit.
#    No `--` after the script name: submit stops parsing at the first
#    non-option, so that `--` would reach benchmark_singlecheck.sh as its first
#    argument and swallow -i. The final `--` IS for benchmark_singlecheck.sh:
#    everything after it goes to every SingleCheck invocation.
#    -t $CPUS explicitly: the baseline hard-codes THREADS=3 while the current
#    version reads $SLURM_CPUS_PER_TASK; without it that difference would be
#    mixed into the result.
# ---------------------------------------------------------------------------
if [ "$DO_AB" -eq 1 ]; then
    echo
    echo ">>> end-to-end A/B (both versions, $REPS reps, one job)"
    ./submit $DRYRUN \
        --cpus "$CPUS" --mem "$MEM" --name bench_ab \
        --sbatch "--export=$EXPORTS" \
        benchmark_singlecheck.sh \
        -i "$BAM" \
        -r "$REPS" \
        -V "$CURSO/SingleCheckBaseline,$REPO" \
        -b SingleCheckBaseline \
        -T "$TOL" \
        -- -t "$CPUS"
fi

# Per-chunk profile of the CURRENT pipeline: which chunk costs what, plus the
# optimized-vs-legacy timings and the equivalence verdicts. --work-scratch so
# the text chunks are timed on node-local disk, like the real pipeline.
if [ "$DO_STAGES" -eq 1 ]; then
    echo
    echo ">>> per-chunk profile + optimized-vs-legacy equivalence"
    ./submit $DRYRUN \
        --cpus "$CPUS" --mem "$MEM" --name bench_stages \
        --sbatch "--export=$EXPORTS" \
        bench/bench_stages.sh \
        -i "$BAM" \
        -r "$REPS" \
        --threads "$CPUS" \
        --warmup \
        --work-scratch
fi

[ -n "$DRYRUN" ] && exit 0

cat <<EOF

Two jobs submitted (use --ab-only / --stages-only to run just one).

1. bench_ab     -> singlecheck_bench_<timestamp>/summary.tsv
     WALL_min      fairest single number (the page cache warms across reps)
     SPEEDUP       vs SingleCheckBaseline
     VS_BASELINE   CONSISTENT = the 13 output columns agree within $TOL

2. bench_stages -> bench/stages_<host>_<timestamp>/report.md
     ranked chunk table, "Optimized vs legacy: do they agree?",
     and the exact-vs-estimate error for the flags that approximate

Each SingleCheck run also prints its own "> Timing summary" phase table and
writes <sample>.SingleCheck.timings.tsv next to the result.
EOF
