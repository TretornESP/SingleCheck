# singlecheck_env.sh -- per-cluster environment for SingleCheck
#
# This file is SOURCED by SingleCheck (and the benchmark) when $SINGLECHECK_ENV
# points at it. Put here everything that is specific to THIS cluster/account but
# common to ALL SingleCheck code variations: module names and tool locations.
# Keeping it separate means you never have to edit each SingleCheck* folder.
#
# Usage:
#   export SINGLECHECK_ENV=/mnt/lustre/scratch/nlsas/home/ulc/cursos/curso385/singlecheck_env.sh
#   ./SingleCheck sample.bam
# or with the benchmarks (they have no --env option: they inherit the exported
# variable, and every SingleCheck they launch sources it):
#   export SINGLECHECK_ENV=$PWD/singlecheck_env.sh
#   ./benchmark_singlecheck.sh -i sample.bam -g 'SingleCheck*' -- -t 8
#   ./bench/bench_stages.sh -i sample.bam -r 5 --threads 32

# --- modules available on this cluster (CESGA / FinisTerrae) -----------------
module load bedtools/2.31.0 2> /dev/null

# BLAST for MetaPhyler. Find the right name with `module spider blast`.
# MetaPhylerV1.13 drives LEGACY blast (blastall); BLAST+ (blastn/blastx) needs
# the marker database to be re-formatted with makeblastdb.
#module load blast-plus/2.11.0 2> /dev/null

# --- tools not provided as modules: give absolute paths ----------------------
# NOTE: this must be the mosdepth EXECUTABLE, not the directory that contains it
# ("...: Is a directory" at run time). Anything already exported (e.g. by
# runner.sh) wins; otherwise try the usual layouts of that directory.
SC_CURSO=/mnt/lustre/scratch/nlsas/home/ulc/cursos/curso385
if [ -z "${MOSDEPTH:-}" ] || [ -d "${MOSDEPTH:-}" ]; then
    for sc_cand in "$SC_CURSO"/mosdepth/bin/mosdepth \
                   "$SC_CURSO"/mosdepth/mosdepth \
                   "$SC_CURSO"/mosdepth; do
        if [ -f "$sc_cand" ] && [ -x "$sc_cand" ]; then
            export MOSDEPTH="$sc_cand"
            break
        fi
    done
    unset sc_cand
fi
export METAPHYLER=$SC_CURSO/MetaPhylerV1.13/metaphyler.pl

# BLAST is looked up on $PATH by MetaPhyler itself, so an absolute path to the
# binary would not help: give the DIRECTORY that holds blastall / blastn and
# SingleCheck prepends it to $PATH. Create one with, for example:
#   conda create -y -p /mnt/lustre/scratch/nlsas/home/ulc/cursos/curso385/blast-legacy \
#       -c bioconda blast-legacy          # blastall + formatdb (MetaPhyler 1.13)
#   conda create -y -p .../blast -c bioconda blast    # BLAST+: blastn, makeblastdb
#export BLAST_DIR=/mnt/lustre/scratch/nlsas/home/ulc/cursos/curso385/blast-legacy/bin
#cuak
# If samtools / R / bwa are NOT auto-loaded on this cluster, add their correct
# module names here as well, e.g.:
#   module load samtools/1.10 2> /dev/null
#   module load gcc R 2> /dev/null
#   module load bwa-mem2 2> /dev/null
