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
# or with the benchmark:
#   ./benchmark_singlecheck.sh --env "$SINGLECHECK_ENV" -i sample.bam -g 'SingleCheck*' -- -t 8

# --- modules available on this cluster (CESGA / FinisTerrae) -----------------
module load bedtools/2.31.0 2> /dev/null

# --- tools not provided as modules: give absolute paths ----------------------
export MOSDEPTH=/mnt/lustre/scratch/nlsas/home/ulc/cursos/curso385/mosdepth
export METAPHYLER=/mnt/lustre/scratch/nlsas/home/ulc/cursos/curso385/MetaPhylerV1.13/metaphyler.pl
#cuak
# If samtools / R / bwa are NOT auto-loaded on this cluster, add their correct
# module names here as well, e.g.:
#   module load samtools/1.10 2> /dev/null
#   module load gcc R 2> /dev/null
#   module load bwa-mem2 2> /dev/null
