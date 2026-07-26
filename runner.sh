#!/bin/bash
# A/B benchmark: current SingleCheck vs the commit before the HPC optimizations.
set -uo pipefail

CURSO=/mnt/lustre/scratch/nlsas/home/ulc/cursos/curso385
cd $CURSO/SingleCheck || exit 1

# ---------------------------------------------------------------------------
# 1. baseline = the commit right before the optimizations (they landed in 9730eb6)
# ---------------------------------------------------------------------------
if [ ! -x "$CURSO/SingleCheckBaseline/SingleCheck" ]; then
    git worktree add $CURSO/SingleCheckBaseline 35ec840 || exit 1
fi

# ---------------------------------------------------------------------------
# 2. mosdepth: $CURSO/mosdepth is a DIRECTORY (conda env / release dir), not the
#    executable -- that is what made the run abort with
#      [MISSING] mosdepth found at .../mosdepth but does not run
#      .../SingleCheck: line 368: .../curso385/mosdepth: Is a directory
#    Locate the real binary and export it BEFORE the pipeline runs.
# ---------------------------------------------------------------------------
MOSDEPTH_BIN=""
for cand in "$CURSO"/mosdepth/bin/mosdepth "$CURSO"/mosdepth/mosdepth "$CURSO"/mosdepth.bin; do
    if [ -f "$cand" ] && [ -x "$cand" ]; then MOSDEPTH_BIN="$cand"; break; fi
done
# last resort: search a couple of levels down
if [ -z "$MOSDEPTH_BIN" ] && [ -d "$CURSO/mosdepth" ]; then
    MOSDEPTH_BIN=$(find "$CURSO/mosdepth" -maxdepth 3 -type f -name 'mosdepth' -perm -u+x 2>/dev/null | head -n1)
fi
# or maybe it is simply on $PATH already
if [ -z "$MOSDEPTH_BIN" ] && command -v mosdepth >/dev/null 2>&1; then
    MOSDEPTH_BIN=$(command -v mosdepth)
fi

if [ -z "$MOSDEPTH_BIN" ]; then
    cat >&2 <<EOF
ERROR: could not find a mosdepth executable.
       \$CURSO/mosdepth is a directory; find the binary inside it with:
           find $CURSO/mosdepth -maxdepth 3 -name 'mosdepth*' -type f
       then set it in singlecheck_env.sh:
           export MOSDEPTH=/full/path/to/mosdepth
EOF
    exit 1
fi

export MOSDEPTH="$MOSDEPTH_BIN"
if ! "$MOSDEPTH" --version >/dev/null 2>&1; then
    echo "ERROR: $MOSDEPTH exists but does not run:" >&2
    "$MOSDEPTH" --version 2>&1 | head -n 3 >&2
    exit 1
fi
echo "mosdepth: $MOSDEPTH ($("$MOSDEPTH" --version 2>&1 | head -n1))"

# ---------------------------------------------------------------------------
# 3. the env file must be exported: benchmark_singlecheck.sh has no --env option,
#    it inherits this and every SingleCheck it launches sources it.
#    Both versions source the SAME file, so the mosdepth fix applies to the
#    baseline worktree too.
# ---------------------------------------------------------------------------
export SINGLECHECK_ENV=$CURSO/SingleCheck/singlecheck_env.sh

# ---------------------------------------------------------------------------
# 4. submit
#    No `--` after the script name: submit stops parsing at the first
#    non-option, so that `--` would reach benchmark_singlecheck.sh as its first
#    argument and swallow -i. The final `--` IS for benchmark_singlecheck.sh:
#    everything after it goes to every SingleCheck invocation.
#    -t 32 explicitly, because the baseline hard-codes THREADS=3 while the
#    current version reads $SLURM_CPUS_PER_TASK -- without it you would be
#    measuring the thread default and everything else at the same time.
# ---------------------------------------------------------------------------
./submit --cpus 32 --mem 128G --name bench_ab benchmark_singlecheck.sh \
    -i $CURSO/Wang5.bam \
    -r 3 \
    -V $CURSO/SingleCheckBaseline,$CURSO/SingleCheck \
    -b SingleCheckBaseline \
    -T 1e-6 \
    -- -t 32
