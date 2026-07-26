#!/bin/bash
###############################################################################
# bench_stages.sh -- per-stage ("chunk") profiling of the SingleCheck pipeline
#
# SingleCheck is an orchestrator: it chains samtools, mosdepth, bedtools, awk /
# sort / zcat glue, four R scripts and Metaphyler. This script splits that chain
# into the individual chunks, runs EACH chunk N times, times every repetition
# and writes a benchmark report telling you which tool eats the wall-clock.
#
# The chunks are copies of the exact command lines used by ../SingleCheck, so
# the timings are representative of a real run. Chunks execute in pipeline
# order inside one private working directory: the artifacts produced by chunk i
# are the input of chunk i+1, exactly as in the real pipeline.
#
# USAGE
#   bench/bench_stages.sh -i <input> [-I <mate2>] [-r REPS] [options]
#
# WRAPPER OPTIONS
#   -i, --input   FILE   input: .bam | .cram | .fastq.gz | _1.fastq.gz  (required)
#   -I, --input2  FILE   second mate for paired-end FASTQ (_2.fastq.gz)
#   -r, --reps    N      repetitions per chunk (default 3)
#   -o, --outdir  DIR    results directory (default bench/stages_<host>_<timestamp>)
#   -s, --stages  LIST   comma-separated subset of chunks to run (default: all)
#   -l, --list           list the chunk keys and exit
#       --warmup         one unmeasured warmup repetition per chunk (warms page cache)
#       --full           ALSO time the whole ../SingleCheck end-to-end, REPS times
#                        (lets you compare sum-of-chunks against the real total)
#       --strict         abort when the dependency check finds a problem, exactly
#                        like ../SingleCheck. Default: report it and measure the
#                        chunks that CAN run, skipping the rest.
#       --no-compare     do not run the LEGACY variants of the optimized chunks.
#                        By default the benchmark runs both the old and the new
#                        implementation of the two rewritten steps, times them
#                        side by side AND verifies they produce the same numbers.
#       --env  FILE      per-cluster environment file to source (see
#                        ../singlecheck_env.sh); defaults to $SINGLECHECK_ENV
#       --keep           keep every intermediate artifact (default)
#       --clean          delete the working directory when finished
#   -h, --help           show this help
#
# PIPELINE OPTIONS (same meaning as in ../SingleCheck)
#       --threads N      -t   number of threads                   (default 3)
#       --window  N      -w   window size                         (default 10000000)
#       --delta   N      -i   autocorrelation increment           (default 1000)
#       --depth   F      -d   downsampling depth                  (default 0.1)
#       --ref     FILE   -r   reference fasta (needed for FASTQ / CRAM)
#       --flag    N      -f   flag to filter out                  (default 772)
#       --mapq    N      -q   mapping quality                     (default 20)
#       --chroms  STR    -c   diploid contigs regex, e.g. '[1-9]'
#       --mt      STR    -m   mitochondrial contig name
#       --include-x      -X   include chromosome X
#       --no-downsample  -N   analyse the original file
#
# OUTPUT (inside --outdir)
#   report.md      human-readable benchmark report (ranked chunks, shares, bars)
#   timings.tsv    one row per chunk and repetition (machine readable)
#   logs/          stdout+stderr and /usr/bin/time -v output of every repetition
#   work/          the private working directory with all artifacts
#
# EXAMPLES
#   # quick local run on a small BAM, 5 repetitions per chunk
#   bench/bench_stages.sh -i test/R1.T15.bam -r 5
#
#   # on the cluster, inside an allocation, 8 threads, also time the full run
#   salloc -c8 --mem 80G
#   export SINGLECHECK_ENV=$PWD/singlecheck_env.sh
#   bench/bench_stages.sh -i /data/HG002.bam -r 3 --threads 8 --warmup --full
#
#   # or submit it
#   ./submit --cpus 8 --mem 80G --time 10:00:00 bench/bench_stages.sh -i /data/HG002.bam -r 3 --threads 8
###############################################################################

usage() { sed -n '2,/^###########/p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//;s/^#//'; }

SELF="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd -P)/$(basename "${BASH_SOURCE[0]}")"
BENCHDIR="$(dirname "$SELF")"
ROOTDIR="$(dirname "$BENCHDIR")"      # repo root: holds SingleCheck and src/

###############################################################################
#                        1. THE CHUNKS                                        #
#                                                                             #
# One function per chunk, holding the very same commands as ../SingleCheck.    #
# They run in a child process (`bench_stages.sh __stage <key>`) that sources   #
# the state file, so a chunk sees exactly the variables the pipeline would     #
# have at that point. Values a later chunk needs are handed back with `emit`.  #
# No `set -e` / `pipefail` here either -- the real pipeline does not use them  #
# (e.g. the `| head -n 1000000` chunk relies on SIGPIPE).                      #
###############################################################################

emit() { printf '%s=%q\n' "$1" "$2" >> "$BENCH_VALS"; }

# --- FASTQ input only: alignment ---------------------------------------------
# $ALIGNER was picked by the dependency check, exactly as in ../SingleCheck.
stage_align() {
    $ALIGNER \
        -t $THREADS $REFERENCE \
        $FASTQ1 $FASTQ2 | \
        samtools sort -@$THREADS -m $SORTMEM -T "${SCRATCH}/sort.$$" \
        -o ${NAME}.bam -
    samtools index -@ $THREADS ${NAME}.bam
}

# GPU alternative to the chunk above (item 8.1). Skipped when pbrun is absent.
stage_align_gpu() {
    pbrun fq2bam \
        --ref "$REFERENCE" \
        --in-fq $FASTQ1 $FASTQ2 \
        --out-bam ${NAME}.gpu.bam \
        --tmp-dir "$SCRATCH" \
        ${SLURM_GPUS:+--num-gpus $SLURM_GPUS}
    samtools index -@ $THREADS ${NAME}.gpu.bam
}

# --- raw sequencing depth ----------------------------------------------------
stage_genome_length() {
    genome_length=$(samtools view -H $SAMREF ${ALN} | grep "^@SQ" | awk '{sum+=substr($3,4,length($3)-1)}END{print sum}')
    emit genome_length "$genome_length"
    echo "genome_length=$genome_length"
}

stage_count_reads() {
    raw_reads=$(samtools view -c -F 2304 -@ $THREADS $SAMREF ${ALN})
    emit raw_reads "$raw_reads"
    echo "raw_reads=$raw_reads"
}

stage_mean_readlen() {
    mean_readlength=$(samtools view -F 2304 $SAMREF ${ALN} | head -n 1000000 | cut -f 10 | awk '{ print length }'| sort | uniq -c | awk '{sum+=$1*$2;num+=$1}END{print sum/num}')
    emit mean_readlength "$mean_readlength"
    echo "mean_readlength=$mean_readlength"
}

# Pure arithmetic (awk + bc): raw bases, depth, downsampling probability.
stage_seq_depth() {
    raw_bases=$(awk -v meanl=$mean_readlength -v rawr=$raw_reads 'BEGIN{print rawr*meanl}')
    sequencing_depth=$(awk -v gl=$genome_length -v rawb=$raw_bases 'BEGIN{print rawb/gl}')
    probability=$(bc -l <<< "scale=10; $downsampling_depth / $sequencing_depth")

    if [[ ! -z $(awk -v prob=$probability 'BEGIN{if (prob >= 1) print "Lowest sequencing than downsampling selected"}') ]] || [[ "$DOWNSAMPLE" = "NO" ]]; then
        DS_MODE="full"    ; DEPTH=$sequencing_depth
    else
        DS_MODE="sample"  ; DEPTH=$downsampling_depth
    fi
    emit raw_bases "$raw_bases"
    emit sequencing_depth "$sequencing_depth"
    emit probability "$probability"
    emit DS_MODE "$DS_MODE"
    emit DEPTH "$DEPTH"
    echo "raw_bases=$raw_bases sequencing_depth=$sequencing_depth probability=$probability mode=$DS_MODE"
}

# --- downsampling ------------------------------------------------------------
stage_downsample() {
    if [ -z "${DS_MODE:-}" ]; then
        echo "downsample: downsampling mode unknown -- the 'seq_depth' chunk must run first" >&2
        exit 2
    fi
    if [ "$DS_MODE" = "full" ]; then
        # Target depth >= actual (or -N): analyse the whole file.
        if [[ "$ALN" == *.bam ]]; then
            ln -sf "$(readlink -f ${ALN})" ${NAME}.${downsampling_depth}X.bam
        else
            samtools view -b -@ $THREADS $SAMREF ${ALN} -o ${NAME}.${downsampling_depth}X.bam
        fi
    else
        frac="${probability#0}"          # bc may emit ".0123" or "0.0123" -> ".0123"
        subsample_arg="1${frac}"         # -> "1.0123" (seed.fraction)
        echo "Downsampling from ${sequencing_depth}X to ${downsampling_depth}X. Fraction kept: $probability"
        samtools view -b -@ $THREADS $SAMREF -s ${subsample_arg} ${ALN} -o ${NAME}.${downsampling_depth}X.bam
    fi
}

stage_index_ds() {
    samtools index -@ $THREADS ${NAME}.${downsampling_depth}X.bam
}

# --- per-base / per-window coverage -----------------------------------------
stage_mosdepth() {
    precision=$(echo $WSIZE | wc -c | awk '{print $1-2}')
    MOSDEPTH_PRECISION=${precision} "$MOSDEPTH" \
        -t $MD_THREADS \
        --fast-mode \
        --by $WSIZE \
        --flag $FLAGTOFILTEROUT \
        --mapq $MAPQUAL \
        ${NAME}.${WSIZE} \
        ${NAME}.${downsampling_depth}X.bam
    # keep mosdepth's status: the `rm` below always succeeds and would otherwise
    # report a failed chunk as OK
    rc=$?
    rm -f ${NAME}.${WSIZE}.mosdepth*
    return $rc
}

# --- autocorrelation input ---------------------------------------------------
stage_shift_track() {
    zcat ${NAME}.${WSIZE}.per-base.bed.gz | \
    awk -v alpha=$DELTA -v print_switch=0 \
    '{if (print_switch==1) {start=$2-alpha; if (start<0) {print $1"\t"0"\t"$3-alpha"\t"$4} else{print $1"\t"start"\t"$3-alpha"\t"$4}}
    else if (alpha >= $2 && alpha < $3) {start=$2-alpha;if (start<0) {print $1"\t"0"\t"$3-alpha"\t"$4} else {print $1"\t"start"\t"$3-alpha"\t"$4}; print_switch=1}}' | \
    $GZIP_CMD > ${NAME}.${WSIZE}.${DELTA}.bed.gz
}

# LEGACY autocorrelation aggregation: whole-genome external sort + RLE.
# Kept so the benchmark can measure -- and validate -- what replaced it.
stage_unionbedg_sort() {
    "$BEDTOOLS" unionbedg -filler NA  \
    -i ${NAME}.${WSIZE}.per-base.bed.gz \
    ${NAME}.${WSIZE}.${DELTA}.bed.gz | \
    grep -E "$DIPLOID_REGEX" | \
    sort --parallel=$THREADS -S 80% --compress-program=$SORTCOMP -T "$SCRATCH" --version-sort -k4 -k5 | \
    awk  \
    '{if (FNR==1){diff=$3-$2;value1=$4;value2=$5}
    else if (value1!=$4 || value2!=$5){print value1"\t"value2"\t"diff; diff=$3-$2;value1=$4;value2=$5}
    else{diff=diff+($3-$2);value1=$4;value2=$5}
    }END{print value1"\t"value2"\t"diff}' \
    > ${NAME}.${DELTA}.shiftedcov.legacy.txt
    # the pipeline's canonical artifact, for the chunks that come after
    cp ${NAME}.${DELTA}.shiftedcov.legacy.txt ${NAME}.${DELTA}.shiftedcov.txt
}

# CURRENT autocorrelation aggregation: hash accumulator, no external sort.
# Writes to a separate file so the runner can prove the two agree.
stage_unionbedg_hash() {
    "$BEDTOOLS" unionbedg -filler NA  \
    -i ${NAME}.${WSIZE}.per-base.bed.gz \
    ${NAME}.${WSIZE}.${DELTA}.bed.gz | \
    grep -E "$DIPLOID_REGEX" | \
    awk '{pair[$4"\t"$5] += $3-$2}
         END{for (p in pair) print p"\t"pair[p]}' \
    > ${NAME}.${DELTA}.shiftedcov.fast.txt
    # the pipeline's canonical artifact, for the chunks that come after
    cp ${NAME}.${DELTA}.shiftedcov.fast.txt ${NAME}.${DELTA}.shiftedcov.txt
}

# --- gini / cv input + R -----------------------------------------------------
stage_freq_table() {
    zcat ${NAME}.${WSIZE}.regions.bed.gz | \
    grep -E "$DIPLOID_REGEX" | \
    awk '{print $4}' | \
    sort -n -T "$SCRATCH" | \
    uniq -c | \
    awk '{print $2"\t"$1}' | sort -n -k1,1 -T "$SCRATCH" \
    > ${NAME}.${WSIZE}.freqs.txt
}

stage_gini_R() {
    Rscript ${SRCDIR}/GiniIndex.R ${NAME}.${WSIZE} $WSIZE
}

stage_cv_R() {
    Rscript ${SRCDIR}/CoefficientOfVariation.R ${NAME}.${WSIZE}
}

stage_autocorr_R() {
    Rscript ${SRCDIR}/Autocorrelation.R ${NAME}.${DELTA} $DELTA
}

# --- MAD input + R -----------------------------------------------------------
stage_contiguous() {
    zcat ${NAME}.${WSIZE}.regions.bed.gz | \
        tail -n +2 | \
        paste  <(zcat ${NAME}.${WSIZE}.regions.bed.gz) - | \
        grep -E "$DIPLOID_REGEX" | \
        awk '{if (NF==8 && $1==$5) {print $4"\t"$8} else if (NF==8 && $1!=$5){print $4}}' | \
        sort -T "$SCRATCH" | uniq -c | \
        awk '{if (NF==3) {print $2"\t"$3"\t"$1}else {print $2"\tNA\t"$1}}' \
        > ${NAME}.${WSIZE}.contiguous.txt
}

stage_mad_R() {
    Rscript ${SRCDIR}/MAD.R ${NAME}.${WSIZE}
}

# --- contamination (full original file) --------------------------------------
stage_unmapped_fasta() {
    samtools view -f 0x4 -@ $THREADS $SAMREF ${ALN} | awk '{OFS="\t"; print ">"$1"\n"$10}' > ${NAME}.unmapped.fasta
}

# ../SingleCheck only warns here and carries on with class=NA; the benchmark
# additionally returns non-zero so the repetition is recorded as FAILED -- the
# time of a MetaPhyler that died early is not a measurement worth reporting.
stage_metaphyler() {
    if ! "$METAPHYLER" 2 ${NAME}.unmapped.fasta ${NAME} || [ ! -s ${NAME}.genus.tab ]; then
        printf "\nWarning: MetaPhyler did not produce %s.genus.tab -- the contaminants\n" "${NAME}" >&2
        printf "         column will be NA. Every other metric is unaffected.\n" >&2
        printf "         Check the installation at %s\n" "$METAPHYLER" >&2
        return 1
    fi
}

# --- final statistics --------------------------------------------------------
# LEGACY: write a primary-only copy of the analysis BAM, index it, then read it
# back three times with idxstats. Kept for comparison.
stage_primary_bam_legacy() {
    samtools view -bF 2304 -@ $THREADS ${NAME}.${downsampling_depth}X.bam > ${NAME}.${downsampling_depth}X.primary.bam
    samtools index -@ $THREADS ${NAME}.${downsampling_depth}X.primary.bam
}

stage_idxstats_legacy() {
    MT_mappedreads=$(samtools idxstats ${NAME}.${downsampling_depth}X.primary.bam | grep -E "$MT_REGEX" | awk '{print $3}')
    mt=$(samtools idxstats ${NAME}.${downsampling_depth}X.primary.bam | awk -v mt=$MT_mappedreads '{sum+=($3+$4)}END{print mt/sum*100}')
    un=$(samtools idxstats ${NAME}.${downsampling_depth}X.primary.bam | \
        awk '{mapped+=$3;unmapped+=$4}END{print unmapped/(unmapped+mapped)*100}')
    printf "%s %s\n" "$mt" "$un" > ${NAME}.mapstats.legacy.txt
    cat ${NAME}.mapstats.legacy.txt
}

# CURRENT: one streaming pass, no copy, no index, no idxstats.
stage_mapstats() {
    samtools view -F 2304 -@ $THREADS $SAMREF ${NAME}.${downsampling_depth}X.bam | \
    awk -v mtre="$MT_REGEX" '
        { total++
          if (int($2/4)%2 == 1) { unmapped++ }
          else if ($3 ~ mtre)   { mt++ } }
        END{ if (total>0) printf "%s %s\n", mt/total*100, unmapped/total*100
             else printf "NA NA\n" }' > ${NAME}.mapstats.txt
    cat ${NAME}.mapstats.txt
}

stage_final_stats() {
    read mt_perc_totalreads unmapped_perc_totalreads < ${NAME}.mapstats.txt
    breadth=$(awk '{if ($1==0){sum+=$3}else if ($1!="NA"){rest+=$3}}END{print 100 - ((sum/(rest+sum))*100)}' ${NAME}.${DELTA}.shiftedcov.txt)
    if [ -s ${NAME}.genus.tab ]; then
        class=$(awk '{if ($1 !~ "{") print $0}' ${NAME}.genus.tab | grep -v "^@" | awk '{print $1"-"$2"-"$3"-"$4"-"$5}' | tr -s '\n' ',' | sed 's/,$/\n/')
    else
        class="NA"
    fi
    SAMPLE=$(basename $NAME)
    WORKDIR=$(dirname "$NAME")
    autocorrelation=$(awk '{print $2}' ${WORKDIR}/Autocorrelation.${SAMPLE}.${DELTA}.txt 2>/dev/null)
    gini=$(awk '{print $2}' ${WORKDIR}/Gini.${SAMPLE}.${WSIZE}.txt 2>/dev/null)
    CV=$(awk '{print $2}' ${WORKDIR}/CV.${SAMPLE}.${WSIZE}.txt 2>/dev/null)
    MAD=$(awk '{print $2}' ${WORKDIR}/MAD.${SAMPLE}.${WSIZE}.txt 2>/dev/null)

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$SAMPLE" "$raw_bases" "$DEPTH" "$WSIZE" "$DELTA" "$unmapped_perc_totalreads" \
        "$mt_perc_totalreads" "$breadth" "$autocorrelation" "$CV" "$gini" "$MAD" "$class" \
        > ${NAME}.SingleCheck.txt
    cat ${NAME}.SingleCheck.txt
}

###############################################################################
#                        2. CHILD MODE                                        #
#            bench_stages.sh __stage <key>   -- runs exactly one chunk        #
###############################################################################
if [ "${1:-}" = "__stage" ]; then
    [ -n "${BENCH_STATE:-}" ] && [ -f "$BENCH_STATE" ] || { echo "child: BENCH_STATE not set" >&2; exit 90; }
    source "$BENCH_STATE"
    : > "$BENCH_VALS"
    cd "$WORK" || exit 91
    "stage_$2"
    exit $?
fi

set -uo pipefail

###############################################################################
#                        3. OPTIONS                                           #
###############################################################################
INPUT="" ; INPUT2="" ; REPS=3 ; OUTDIR="" ; STAGES_CSV="" ; LIST=0
WARMUP=0 ; FULL=0 ; CLEAN=0 ; STRICT=0 ; COMPARE=1 ; ENVFILE="${SINGLECHECK_ENV:-}"

# pipeline defaults -- kept identical to ../SingleCheck
WSIZE=10000000
DELTA=1000
FLAGTOFILTEROUT=772
THREADS=3
MAPQUAL=20
DIPLOID_REGEX="^(chr)*[1-9]"
MT_REGEX='^(MT|chrM)'
DOWNSAMPLE="YES"
downsampling_depth=0.1
ds_strategy="ConstantMemory"
REFERENCE=""

while [ "$#" -gt 0 ]; do
  case "$1" in
    -i|--input)      INPUT="$2"; shift 2 ;;
    -I|--input2)     INPUT2="$2"; shift 2 ;;
    -r|--reps)       REPS="$2"; shift 2 ;;
    -o|--outdir)     OUTDIR="$2"; shift 2 ;;
    -s|--stages)     STAGES_CSV="$2"; shift 2 ;;
    -l|--list)       LIST=1; shift ;;
    --warmup)        WARMUP=1; shift ;;
    --full)          FULL=1; shift ;;
    --strict)        STRICT=1; shift ;;
    --no-compare)    COMPARE=0; shift ;;
    --env)           ENVFILE="$2"; shift 2 ;;
    --keep)          CLEAN=0; shift ;;
    --clean)         CLEAN=1; shift ;;
    --threads)       THREADS="$2"; shift 2 ;;
    --window)        WSIZE="$2"; shift 2 ;;
    --delta)         DELTA="$2"; shift 2 ;;
    --depth)         downsampling_depth="$2"; shift 2 ;;
    --ref)           REFERENCE="$2"; shift 2 ;;
    --flag)          FLAGTOFILTEROUT="$2"; shift 2 ;;
    --mapq)          MAPQUAL="$2"; shift 2 ;;
    --chroms)        DIPLOID_REGEX='^('"$2"')'; shift 2 ;;
    --mt)            MT_REGEX='^'"$2"; shift 2 ;;
    --include-x)     DIPLOID_REGEX="^(chr)*[1-9|X]"; shift ;;
    --no-downsample) DOWNSAMPLE="NO"; shift ;;
    -h|--help)       usage; exit 0 ;;
    *) printf 'Unknown option: %s\n\n' "$1" >&2; usage >&2; exit 1 ;;
  esac
done

if [ "$LIST" -eq 1 ]; then
  cat <<'EOF'
Chunks, in pipeline order (key -- what it runs):

  align            bwa-mem2/bwa mem | samtools sort ; samtools index   [FASTQ input only]
  align_gpu        pbrun fq2bam (NVIDIA Parabricks)                    [FASTQ + GPU only]
  genome_length    samtools view -H | grep @SQ | awk
  count_reads      samtools view -c -F 2304                            [full input file]
  mean_readlen     samtools view | head -1000000 | cut | sort | uniq   [full input file]
  seq_depth        awk + bc arithmetic (bases, depth, probability)
  downsample       samtools view -b -s SEED.FRACTION                   [full input file]
  index_ds         samtools index (downsampled BAM)
  mosdepth         mosdepth --fast-mode --by WSIZE
  shift_track      zcat per-base | awk shift | pigz/bgzip/gzip
  unionbedg_sort   LEGACY: unionbedg | grep | sort --version-sort | awk RLE
  unionbedg_hash   unionbedg | grep | awk hash accumulator (no external sort)
  freq_table       zcat regions | grep | awk | sort | uniq -c | awk
  gini_R           Rscript src/GiniIndex.R
  cv_R             Rscript src/CoefficientOfVariation.R
  autocorr_R       Rscript src/Autocorrelation.R
  contiguous       zcat regions | paste | grep | awk | sort | uniq -c
  mad_R            Rscript src/MAD.R
  unmapped_fasta   samtools view -f 0x4 | awk                          [full input file]
  metaphyler       metaphyler.pl (contamination)
  primary_bam_legacy  LEGACY: samtools view -bF 2304 copy + samtools index
  idxstats_legacy     LEGACY: samtools idxstats x3 + awk
  mapstats            one streaming pass over the analysis BAM (MT%, unmapped%)
  final_stats         assemble the result line (breadth + the metric files)
EOF
  exit 0
fi

[ -n "$INPUT" ] || { echo "Error: --input is required" >&2; usage >&2; exit 1; }
[ -f "$INPUT" ] || { echo "Error: input not found: $INPUT" >&2; exit 1; }
[ -z "$INPUT2" ] || [ -f "$INPUT2" ] || { echo "Error: mate not found: $INPUT2" >&2; exit 1; }
[[ "$REPS" =~ ^[0-9]+$ ]] && [ "$REPS" -ge 1 ] || { echo "Error: --reps must be >= 1" >&2; exit 1; }

case "$INPUT" in
  *.bam|*.cram)              METHOD="Aligned" ;;
  *_1.fastq.gz)              METHOD="Paired-end" ;;
  *.fastq.gz)                METHOD="Single-end" ;;
  *) echo "Error: input must be .bam, .cram or .fastq.gz" >&2; exit 1 ;;
esac
if [ "$METHOD" = "Paired-end" ] && [ -z "$INPUT2" ]; then
  echo "Note: _1.fastq.gz given without -I mate2 -- treating as single-end" >&2
  METHOD="Single-end"
fi

###############################################################################
#                        4. ENVIRONMENT (same as ../SingleCheck)              #
#                                                                             #
# `set +u` while sourcing: module init files and conda's activate reference    #
# unset variables, which would abort a `set -u` shell.                        #
###############################################################################
set +u
MODULES_LOADED=""
MODULES_FAILED=""
HAVE_MODULE=0
type module >/dev/null 2>&1 && HAVE_MODULE=1

# try_module <spec> [alternative spec ...] -- same helper as ../SingleCheck
try_module() {
    local spec
    [ "$HAVE_MODULE" -eq 1 ] || return 1
    for spec in "$@"; do
        if module load $spec 2>/dev/null; then
            MODULES_LOADED="${MODULES_LOADED} ${spec}"
            return 0
        fi
    done
    MODULES_FAILED="${MODULES_FAILED} $1"
    return 1
}

[ "$HAVE_MODULE" -eq 1 ] && module purge 2> /dev/null
try_module samtools/1.10 samtools
try_module "gcc/6.4.0 R/3.6.3" R
try_module miniconda3/4.8.2 miniconda3
try_module "gcccore/6.4.0 bedtools/2.28.0" bedtools
try_module "gcccore/6.4.0 bwa-mem2/2.0" bwa-mem2
try_module "gcc/6.4.0 bwa/0.7.17" bwa
try_module blast-plus blast blast+ ncbi-blast blast-legacy
source activate /mnt/netapp1/posadalab/APPS/CommonCondaEnvironments/mosdepth 2> /dev/null

if [ -n "$ENVFILE" ] && [ -f "$ENVFILE" ]; then source "$ENVFILE"; fi
set -u

# MetaPhyler looks BLAST up on $PATH by name, so $BLAST_DIR is a DIRECTORY that
# gets prepended to $PATH -- same handling as ../SingleCheck.
if [ -n "${BLAST_DIR:-}" ]; then
  if [ -d "$BLAST_DIR" ]; then
    case ":$PATH:" in *":$BLAST_DIR:"*) ;; *) PATH="$BLAST_DIR:$PATH"; export PATH ;; esac
  else
    printf "\nWarning: \$BLAST_DIR=%s is not a directory -- ignored\n" "$BLAST_DIR" >&2
  fi
fi

MOSDEPTH="${MOSDEPTH:-mosdepth}"
BEDTOOLS="${BEDTOOLS:-bedtools}"
METAPHYLER="${METAPHYLER:-$HOME/apps/Metaphyler/MetaPhylerSRV0.115/metaphyler.pl}"
SRCDIR="${ROOTDIR}/src"

export LC_ALL=C
SCRATCH="${SLURM_TMPDIR:-${TMPDIR:-/tmp}}"
export TMPDIR="$SCRATCH"

if command -v zstd >/dev/null 2>&1; then SORTCOMP=zstd
elif command -v pigz >/dev/null 2>&1; then SORTCOMP=pigz
else SORTCOMP=gzip; fi

if command -v pigz >/dev/null 2>&1; then GZIP_CMD="pigz -c -p $THREADS"
elif command -v bgzip >/dev/null 2>&1; then GZIP_CMD="bgzip -@ $THREADS -c"
else GZIP_CMD="gzip -c"; fi

# Same thread split and sort sizing as ../SingleCheck (items 5.1, 6.6, 11).
# No mosdepth cap by default: measured slower at -t 4 than at -t 32.
MD_THREADS="${SINGLECHECK_MOSDEPTH_THREADS:-$THREADS}"
case "$MD_THREADS" in ''|*[!0-9]*) MD_THREADS=$THREADS ;; esac
[ "$MD_THREADS" -lt 1 ] && MD_THREADS=1
[ "$MD_THREADS" -gt "$THREADS" ] && MD_THREADS=$THREADS

sort_mem_mb=0
if [ -n "${SLURM_MEM_PER_NODE:-}" ]; then
    sort_mem_mb=$(( ${SLURM_MEM_PER_NODE%%[!0-9]*} * 70 / 100 / THREADS ))
elif [ -n "${SLURM_MEM_PER_CPU:-}" ]; then
    sort_mem_mb=$(( ${SLURM_MEM_PER_CPU%%[!0-9]*} * 70 / 100 ))
fi
[ "$sort_mem_mb" -lt 256 ] && sort_mem_mb=768
SORTMEM="${sort_mem_mb}M"

SAMREF=""
[ -n "$REFERENCE" ] && SAMREF="--reference $REFERENCE"

###############################################################################
#                        5. WORKING DIRECTORY                                 #
###############################################################################
STAMP="$(date +%Y%m%d_%H%M%S)"
OUTDIR="${OUTDIR:-$BENCHDIR/stages_$(hostname -s)_${STAMP}}"
mkdir -p "$OUTDIR/logs" "$OUTDIR/work" || exit 1
OUTDIR="$(cd "$OUTDIR" && pwd -P)"
WORK="$OUTDIR/work"
STATE="$OUTDIR/state.env"
BENCH_VALS="$OUTDIR/vals.env"
TIMINGS="$OUTDIR/timings.tsv"
REPORT="$OUTDIR/report.md"
export BENCH_STATE="$STATE" BENCH_VALS

# Inputs are SYMLINKED, never copied (safe for 300 GB - TB files).
ABS_INPUT="$(readlink -f "$INPUT")"
BASE_IN="$(basename "$INPUT")"
ln -sf "$ABS_INPUT" "$WORK/$BASE_IN"
[ -f "${ABS_INPUT}.bai" ]  && ln -sf "${ABS_INPUT}.bai"  "$WORK/${BASE_IN}.bai"
[ -f "${ABS_INPUT}.crai" ] && ln -sf "${ABS_INPUT}.crai" "$WORK/${BASE_IN}.crai"
if [ -n "$INPUT2" ]; then
  ABS_INPUT2="$(readlink -f "$INPUT2")"; BASE_IN2="$(basename "$INPUT2")"
  ln -sf "$ABS_INPUT2" "$WORK/$BASE_IN2"
fi

FILE="$WORK/$BASE_IN"
FASTQ1="" ; FASTQ2=""
if [ "$METHOD" != "Aligned" ]; then
  FASTQ1="$FILE"
  [ -n "$INPUT2" ] && FASTQ2="$WORK/$BASE_IN2"
fi
NAME=$(echo $FILE | sed 's/.bam$//' | sed 's/.cram$//' | sed 's/_1.fastq.gz$//' | sed 's/.fastq.gz$//')

if [ "$METHOD" = "Aligned" ]; then ALN="$FILE"; else ALN="${NAME}.bam"; fi
case "$ALN" in
  *.cram) ALNIDX="${ALN}.crai" ;;
  *)      ALNIDX="${ALN}.bai"  ;;
esac

###############################################################################
#              5b. DEPENDENCY CHECK -- mirrors ../SingleCheck                  #
#                                                                             #
# Same checks, same messages, same order as the CHECK DEPENDENCIES section of  #
# ../SingleCheck, so the benchmark runs in exactly the environment the         #
# pipeline would. KEEP THE TWO IN SYNC when either changes.                    #
#                                                                             #
# Difference: a missing tool does not abort here by default -- the chunks that #
# need it are reported as SKIPPED and everything else is still measured, which #
# is the point of a benchmark. Use --strict to abort like ../SingleCheck does. #
###############################################################################
echo "> Checking dependencies"
missing=0
warnings=0

err()   { printf "  [MISSING] %s\n" "$1" >&2; missing=$((missing+1)); }
warn()  { printf "  [WARN]    %s\n" "$1" >&2; warnings=$((warnings+1)); }
found() { printf "  [ok]      %s\n" "$1"; }

have()    { command -v "$1" >/dev/null 2>&1 || [ -x "$1" ]; }
resolve() { command -v "$1" 2>/dev/null || printf '%s\n' "$1"; }

version_ge() {
    awk -v a="$1" -v b="$2" 'BEGIN{
        na=split(a,A,"."); nb=split(b,B,"."); n=(na>nb?na:nb)
        for(i=1;i<=n;i++){ x=A[i]+0; y=B[i]+0; if(x>y) exit 0; if(x<y) exit 1 }
        exit 0}'
}

require_tool() {
    local handle="$1" pretty="$2" minv="$3" hint="$4"; shift 4
    local path out ver
    if ! have "$handle"; then
        err "$pretty not found (looked for '$handle')${hint:+ -- $hint}"
        return 1
    fi
    path="$(resolve "$handle")"
    if ! out="$("$@" 2>&1)"; then
        err "$pretty found at $path but does not run${hint:+ -- $hint}"
        printf "%s\n" "$out" | head -n 3 | sed 's/^/            /' >&2
        return 1
    fi
    if [ "$minv" = "-" ]; then
        found "$pretty ($path)"
        return 0
    fi
    ver=$(printf "%s\n" "$out" | grep -oE '[0-9]+(\.[0-9]+)+' | head -n1)
    if [ -z "$ver" ]; then
        warn "could not determine the version of $pretty at $path (need >= $minv)"
        found "$pretty (unknown version, $path)"
    elif version_ge "$ver" "$minv"; then
        found "$pretty $ver ($path)"
    else
        err "$pretty $ver at $path is too old, need >= $minv${hint:+ -- $hint}"
    fi
}

# --- environment modules -----------------------------------------------------
if [ "$HAVE_MODULE" -eq 1 ]; then
    [ -n "$MODULES_LOADED" ] && found "modules loaded:$MODULES_LOADED"
    if [ -n "$MODULES_FAILED" ]; then
        if [ -n "${SINGLECHECK_STRICT_MODULES:-}" ]; then
            err "module(s) not available on this cluster:$MODULES_FAILED (SINGLECHECK_STRICT_MODULES=1)"
        else
            warn "module(s) not available on this cluster:$MODULES_FAILED -- only a problem if a tool below is missing"
        fi
    fi
else
    warn "no 'module' command available -- every tool must be on \$PATH or set in \$SINGLECHECK_ENV"
fi

# --- shell / coreutils helpers the chunks pipe data through ------------------
core_missing=""
for c in grep sed awk sort uniq cut head tail paste tr wc zcat gzip ln readlink basename dirname df; do
    have "$c" || core_missing="$core_missing $c"
done
if [ -n "$core_missing" ]; then
    err "basic shell utilities not found on \$PATH:$core_missing"
else
    found "shell utilities (grep sed awk sort uniq cut head tail paste tr wc zcat gzip)"
fi

if have sort; then
    if printf '' | sort --parallel=2 --compress-program="$SORTCOMP" -T "$SCRATCH" --version-sort >/dev/null 2>&1; then
        found "GNU sort (--parallel, --compress-program=$SORTCOMP, --version-sort)"
    else
        err "'sort' at $(resolve sort) does not support --parallel / --compress-program=$SORTCOMP / --version-sort (GNU coreutils sort required)"
    fi
fi

if ! have bc; then
    err "command 'bc' not found on \$PATH (needed for the downsampling probability)"
elif [ -z "$(printf 'scale=10; 1/3\n' | bc -l 2>/dev/null)" ]; then
    err "'bc' at $(resolve bc) is present but produced no output for a test computation"
else
    found "bc ($(resolve bc))"
fi

# --- the tools the chunks actually measure -----------------------------------
require_tool samtools "samtools" 1.9 \
    "module load samtools/1.10" \
    samtools --version
require_tool "$BEDTOOLS" "bedtools" 2.25 \
    "module load bedtools, or set BEDTOOLS=/path/to/bedtools in \$SINGLECHECK_ENV" \
    "$BEDTOOLS" --version
require_tool "$MOSDEPTH" "mosdepth" 0.2.5 \
    "activate the mosdepth conda env, or set MOSDEPTH=/path/to/mosdepth in \$SINGLECHECK_ENV" \
    "$MOSDEPTH" --version
require_tool Rscript "R (Rscript)" 3.5 \
    "module load R" \
    Rscript --version

# The compressor matters for the timings: say which one the chunks will use.
have pigz  || have bgzip || warn "neither pigz nor bgzip found -- shift_track will be timed with single-threaded gzip"
have zstd  || warn "zstd not found -- sort spill files will be compressed with $SORTCOMP"

# --- R packages used by the four metric scripts ------------------------------
if have Rscript; then
    rmissing=$(Rscript -e 'p <- c("tidyr","dplyr","matrixStats"); m <- p[!sapply(p, requireNamespace, quietly=TRUE)]; cat(paste(m, collapse=" "))' 2>/dev/null)
    rrc=$?
    if [ "$rrc" -ne 0 ]; then
        err "could not query the R library (Rscript exited $rrc) -- check the R module/installation"
    elif [ -n "$rmissing" ]; then
        err "R package(s) not installed: $rmissing -- install.packages(c(\"tidyr\",\"dplyr\",\"matrixStats\"))"
    else
        found "R packages tidyr, dplyr, matrixStats"
    fi
fi

# --- helper scripts shipped with the pipeline --------------------------------
if [ ! -d "$SRCDIR" ]; then
    err "helper directory $SRCDIR not found -- bench/ must sit next to the pipeline's src/ folder"
else
    rs_missing=""
    for rs in GiniIndex.R CoefficientOfVariation.R Autocorrelation.R MAD.R; do
        [ -r "${SRCDIR}/${rs}" ] || rs_missing="$rs_missing $rs"
    done
    if [ -n "$rs_missing" ]; then
        err "helper script(s) missing or unreadable in ${SRCDIR}:$rs_missing"
    else
        found "helper R scripts ($SRCDIR)"
    fi
fi

# --- Metaphyler, the BLAST it shells out to, and its own internal paths ------
if [ ! -f "$METAPHYLER" ]; then
    err "Metaphyler not found at $METAPHYLER (set \$METAPHYLER to the full path of metaphyler.pl)"
elif [ ! -x "$METAPHYLER" ]; then
    err "Metaphyler at $METAPHYLER is not executable (chmod +x it)"
else
    have perl || err "perl not found on \$PATH (needed to run $METAPHYLER)"
    found "Metaphyler ($METAPHYLER)"

    blast_wanted=$(grep -ohE '\b(blastall|blastn|blastx|blastp|tblastn|megablast)\b' "$METAPHYLER" 2>/dev/null | sort -u)
    blast_guessed=0
    if [ -z "$blast_wanted" ]; then
        blast_wanted="blastall blastn"
        blast_guessed=1
    fi
    blast_have="" ; blast_first="" ; blast_miss=""
    for b in $blast_wanted; do
        if have "$b"; then
            blast_have="$blast_have $b"
            [ -n "$blast_first" ] || blast_first="$(resolve "$b")"
        else
            blast_miss="$blast_miss $b"
        fi
    done
    if [ -n "$blast_have" ]; then
        found "BLAST:$blast_have ($blast_first)"
        if [ -n "$blast_miss" ] && [ "$blast_guessed" -eq 0 ]; then
            warn "$(basename "$METAPHYLER") also refers to:$blast_miss -- not on \$PATH"
        fi
    elif [ "$blast_guessed" -eq 1 ]; then
        err "no BLAST on \$PATH (MetaPhyler needs 'blastall' for legacy BLAST or 'blastn' for BLAST+) -- module load blast, or set BLAST_DIR=/path/to/blast/bin in \$SINGLECHECK_ENV"
    else
        err "no BLAST on \$PATH, but $(basename "$METAPHYLER") calls:$blast_miss -- module load blast (conda: 'blast-legacy' provides blastall, 'blast' provides blastn), or set BLAST_DIR=/path/to/blast/bin in \$SINGLECHECK_ENV"
    fi

    METADIR=$(dirname "$METAPHYLER")
    mp_missing=""
    for helper in $(grep -ohE '/[A-Za-z0-9_./-]+\.pl' "$METAPHYLER" 2>/dev/null | sort -u); do
        [ "$helper" = "$METAPHYLER" ] && continue
        root="/$(printf '%s' "$helper" | cut -d/ -f2)"
        [ -d "$root" ] || continue
        [ -e "$helper" ] || mp_missing="$mp_missing $helper"
    done
    if [ -n "$mp_missing" ]; then
        err "MetaPhyler refers to helper script(s) that do not exist:$mp_missing"
        printf "            metaphyler.pl has its install path baked in, so the installation was\n" >&2
        printf "            moved or copied. Re-install it where it lives now:\n" >&2
        printf "                cd %s && perl installMetaphyler.pl\n" "$METADIR" >&2
    fi
    if [ ! -d "${METADIR}/markers" ]; then
        warn "MetaPhyler marker directory ${METADIR}/markers not found -- classification may fail"
    elif [ -z "$(ls "${METADIR}"/markers/*hr 2>/dev/null | head -n1)" ]; then
        warn "no formatted BLAST database in ${METADIR}/markers -- run: perl ${METADIR}/installMetaphyler.pl"
    else
        found "MetaPhyler marker database (${METADIR}/markers)"
    fi
fi

# --- aligner + its reference index (FASTQ input only) ------------------------
# $ALIGNER is chosen here and used by the align chunk, exactly as ../SingleCheck
# chooses it in its dependency check and uses it in the mapping step.
ALIGNER=""
if [ "$METHOD" = "Paired-end" ] || [ "$METHOD" = "Single-end" ]; then
    [ -n "$REFERENCE" ] || err "FASTQ input requires a reference genome (--ref ref.fa)"
    if command -v bwa-mem2 >/dev/null 2>&1 && [ -n "$REFERENCE" ] && [ -f "${REFERENCE}.bwt.2bit.64" ]; then
        ALIGNER="bwa-mem2 mem"
        idx_missing=""
        for ext in 0123 amb ann bwt.2bit.64 pac; do
            [ -f "${REFERENCE}.${ext}" ] || idx_missing="$idx_missing ${REFERENCE}.${ext}"
        done
        if [ -n "$idx_missing" ]; then
            err "incomplete bwa-mem2 index (missing:$idx_missing) -- run: bwa-mem2 index $REFERENCE"
        else
            found "aligner bwa-mem2 ($(resolve bwa-mem2)) + its index"
        fi
    elif command -v bwa >/dev/null 2>&1; then
        ALIGNER="bwa mem"
        if [ -n "$REFERENCE" ]; then
            idx_missing=""
            for ext in amb ann bwt pac sa; do
                [ -f "${REFERENCE}.${ext}" ] || idx_missing="$idx_missing ${REFERENCE}.${ext}"
            done
            if [ -n "$idx_missing" ]; then
                err "incomplete bwa index (missing:$idx_missing) -- run: bwa index $REFERENCE"
            else
                found "aligner bwa ($(resolve bwa)) + its index"
            fi
        fi
    else
        err "no aligner on \$PATH (need bwa-mem2 or bwa) -- module load bwa-mem2 / bwa"
    fi
fi

# --- reference genome --------------------------------------------------------
if [ -n "$REFERENCE" ]; then
    if [ ! -f "$REFERENCE" ]; then
        err "reference fasta not found: $REFERENCE"
    elif [ ! -r "$REFERENCE" ]; then
        err "reference fasta not readable: $REFERENCE"
    elif [ ! -f "${REFERENCE}.fai" ]; then
        err "samtools faidx index missing: ${REFERENCE}.fai -- run: samtools faidx $REFERENCE"
    else
        found "reference $REFERENCE (+ .fai)"
    fi
fi

# --- input file(s) and index -------------------------------------------------
if [ "$METHOD" = "Aligned" ]; then
    if [ ! -r "$ALN" ]; then
        err "input file not readable: $ALN"
    elif [ ! -s "$ALN" ]; then
        err "input file is empty: $ALN"
    else
        if have samtools && ! samtools quickcheck "$ALN" >/dev/null 2>&1; then
            if [[ "$ALN" == *.cram ]]; then
                warn "samtools quickcheck is unhappy about $ALN -- check the file and its reference"
            else
                err "$ALN is not a valid/complete BAM (samtools quickcheck failed: truncated or not a BAM?)"
            fi
        fi
        if [ ! -f "$ALNIDX" ]; then
            err "index missing: $ALNIDX -- run: samtools index $ABS_INPUT"
        else
            [ "$ALNIDX" -ot "$ALN" ] && warn "index $ALNIDX is older than $ALN -- consider re-running samtools index"
            found "input $ABS_INPUT (+ index, symlinked into work/)"
        fi
    fi
else
    for fq in "$FASTQ1" "$FASTQ2"; do
        [ -n "$fq" ] || continue
        if [ ! -r "$fq" ]; then
            err "FASTQ not readable: $fq"
        elif [ ! -s "$fq" ]; then
            err "FASTQ is empty: $fq"
        else
            found "input $fq"
        fi
    done
fi

# --- directories ------------------------------------------------------------
if [ -d "$WORK" ] && [ -w "$WORK" ]; then
    found "working directory $WORK is writable"
else
    err "working directory $WORK is not writable"
fi
if [ -d "$SCRATCH" ] && [ -w "$SCRATCH" ]; then
    found "scratch directory $SCRATCH is writable"
else
    err "scratch directory $SCRATCH does not exist or is not writable (set \$TMPDIR)"
fi

if have df; then
    free_kb=$(df -Pk "$SCRATCH" 2>/dev/null | awk 'NR==2{print $4}')
    case "$free_kb" in
        ''|*[!0-9]*) : ;;
        *) [ "$free_kb" -lt 5000000 ] && warn "only $((free_kb/1024)) MB free in $SCRATCH -- the sort step may run out of space" ;;
    esac
fi

# Oversubscribed threads would make every timing meaningless, so this one is
# worth more here than in the pipeline itself.
ncpus="${SLURM_CPUS_PER_TASK:-$(nproc 2>/dev/null || echo 0)}"
case "$ncpus" in ''|*[!0-9]*) ncpus=0 ;; esac
if [ "$ncpus" -gt 0 ] && [ "$THREADS" -gt "$ncpus" ]; then
    warn "--threads $THREADS is larger than the $ncpus CPU(s) available -- timings will be distorted"
fi

# --- verdict ----------------------------------------------------------------
if [ "$missing" -ne 0 ]; then
    if [ "$STRICT" -eq 1 ]; then
        printf "\nError: %d dependency problem(s) found (see the [MISSING] lines above).\nAborting (--strict).\n" "$missing" >&2
        exit 1
    fi
    printf "  %d dependency problem(s): the chunks that need them will be reported as SKIPPED\n" "$missing"
    printf "  (run with --strict to abort instead, like ../SingleCheck does)\n"
fi
[ "$warnings" -ne 0 ] && printf "  %d warning(s) -- continuing\n" "$warnings"
[ "$missing" -eq 0 ] && echo "  all dependencies OK"
echo

# values filled in by the chunks themselves
genome_length="" ; raw_reads="" ; mean_readlength="" ; raw_bases=""
sequencing_depth="" ; probability="" ; DEPTH="" ; DS_MODE=""

STATE_VARS=(FILE FASTQ1 FASTQ2 METHOD NAME ALN ALNIDX SAMREF REFERENCE THREADS
            WSIZE DELTA FLAGTOFILTEROUT MAPQUAL DIPLOID_REGEX MT_REGEX
            downsampling_depth ds_strategy DOWNSAMPLE SCRATCH SORTCOMP GZIP_CMD
            MOSDEPTH BEDTOOLS METAPHYLER SRCDIR WORK BENCH_VALS ALIGNER
            MD_THREADS SORTMEM
            genome_length raw_reads mean_readlength raw_bases sequencing_depth
            probability DEPTH DS_MODE)

save_state() {
  local v
  : > "$STATE"
  for v in "${STATE_VARS[@]}"; do printf '%s=%q\n' "$v" "${!v-}" >> "$STATE"; done
}
save_state

###############################################################################
#                        6. CHUNK TABLE                                       #
#   NEEDS_CMD  : commands/executables required (skip the chunk if absent)     #
#   NEEDS_FILE : artifacts required as input   (skip the chunk if absent)     #
#   RESET      : artifacts deleted before EVERY repetition, so each rep does  #
#                the full amount of work instead of overwriting a warm file   #
###############################################################################
DS="${NAME}.${downsampling_depth}X"
SAMPLE="$(basename "$NAME")"

STAGE_KEYS=()
declare -A S_GROUP=() S_DESC=() S_CMD=() S_FILE=() S_RESET=()

add_stage() { # key group desc needs_cmd needs_file reset
  STAGE_KEYS+=("$1")
  S_GROUP["$1"]="$2"; S_DESC["$1"]="$3"
  S_CMD["$1"]="$4"; S_FILE["$1"]="$5"; S_RESET["$1"]="$6"
}

if [ "$METHOD" != "Aligned" ]; then
  add_stage align "Alignment" "bwa mem | samtools sort -m $SORTMEM ; samtools index" \
    "samtools" "$FASTQ1" "${NAME}.bam ${NAME}.bam.bai"
  add_stage align_gpu "Alignment" "GPU: pbrun fq2bam (Parabricks)" \
    "pbrun" "$FASTQ1" "${NAME}.gpu.bam ${NAME}.gpu.bam.bai"
fi

add_stage genome_length "Raw depth" "samtools view -H | awk (genome length)" \
  "samtools" "$ALN" ""
add_stage count_reads "Raw depth" "samtools view -c -F 2304 (full file)" \
  "samtools" "$ALN" ""
add_stage mean_readlen "Raw depth" "samtools view | head -1e6 | sort | uniq -c" \
  "samtools" "$ALN" ""
add_stage seq_depth "Raw depth" "awk + bc (bases, depth, probability)" \
  "bc" "" ""
add_stage downsample "Downsampling" "samtools view -b -s (full file)" \
  "samtools" "$ALN" "${DS}.bam"
add_stage index_ds "Downsampling" "samtools index (downsampled BAM)" \
  "samtools" "${DS}.bam" "${DS}.bam.bai"
add_stage mosdepth "Coverage" "mosdepth --fast-mode --by $WSIZE" \
  "$MOSDEPTH" "${DS}.bam ${DS}.bam.bai" \
  "${NAME}.${WSIZE}.per-base.bed.gz ${NAME}.${WSIZE}.per-base.bed.gz.csi ${NAME}.${WSIZE}.regions.bed.gz ${NAME}.${WSIZE}.regions.bed.gz.csi"
add_stage shift_track "Autocorrelation" "zcat per-base | awk shift | gzip" \
  "zcat" "${NAME}.${WSIZE}.per-base.bed.gz" "${NAME}.${WSIZE}.${DELTA}.bed.gz"
if [ "$COMPARE" -eq 1 ]; then
  add_stage unionbedg_sort "Autocorrelation" "LEGACY: unionbedg | grep | sort --version-sort | awk RLE" \
    "$BEDTOOLS" "${NAME}.${WSIZE}.per-base.bed.gz ${NAME}.${WSIZE}.${DELTA}.bed.gz" \
    "${NAME}.${DELTA}.shiftedcov.legacy.txt ${NAME}.${DELTA}.shiftedcov.txt"
fi
add_stage unionbedg_hash "Autocorrelation" "bedtools unionbedg | grep | awk hash (sort-free)" \
  "$BEDTOOLS" "${NAME}.${WSIZE}.per-base.bed.gz ${NAME}.${WSIZE}.${DELTA}.bed.gz" \
  "${NAME}.${DELTA}.shiftedcov.fast.txt ${NAME}.${DELTA}.shiftedcov.txt"
add_stage freq_table "Gini/CV" "zcat regions | sort | uniq -c (freq table)" \
  "zcat" "${NAME}.${WSIZE}.regions.bed.gz" "${NAME}.${WSIZE}.freqs.txt"
add_stage gini_R "Gini/CV" "Rscript GiniIndex.R" \
  "Rscript" "${NAME}.${WSIZE}.freqs.txt" "${WORK}/Gini.${SAMPLE}.${WSIZE}.txt"
add_stage cv_R "Gini/CV" "Rscript CoefficientOfVariation.R" \
  "Rscript" "${NAME}.${WSIZE}.freqs.txt" "${WORK}/CV.${SAMPLE}.${WSIZE}.txt"
add_stage autocorr_R "Autocorrelation" "Rscript Autocorrelation.R" \
  "Rscript" "${NAME}.${DELTA}.shiftedcov.txt" "${WORK}/Autocorrelation.${SAMPLE}.${DELTA}.txt"
add_stage contiguous "MAD" "zcat regions | paste | sort | uniq -c" \
  "zcat" "${NAME}.${WSIZE}.regions.bed.gz" "${NAME}.${WSIZE}.contiguous.txt"
add_stage mad_R "MAD" "Rscript MAD.R" \
  "Rscript" "${NAME}.${WSIZE}.contiguous.txt" "${WORK}/MAD.${SAMPLE}.${WSIZE}.txt"
add_stage unmapped_fasta "Contamination" "samtools view -f 0x4 | awk (full file)" \
  "samtools" "$ALN" "${NAME}.unmapped.fasta"
add_stage metaphyler "Contamination" "metaphyler.pl" \
  "$METAPHYLER" "${NAME}.unmapped.fasta" ""
if [ "$COMPARE" -eq 1 ]; then
  add_stage primary_bam_legacy "Final stats" "LEGACY: samtools view -bF 2304 copy + index" \
    "samtools" "${DS}.bam" "${DS}.primary.bam ${DS}.primary.bam.bai"
  add_stage idxstats_legacy "Final stats" "LEGACY: samtools idxstats x3 + awk" \
    "samtools" "${DS}.primary.bam" "${NAME}.mapstats.legacy.txt"
fi
add_stage mapstats "Final stats" "one streaming pass (MT %, unmapped %)" \
  "samtools" "${DS}.bam" "${NAME}.mapstats.txt"
add_stage final_stats "Final stats" "assemble the result line" \
  "samtools" "${NAME}.mapstats.txt ${NAME}.${DELTA}.shiftedcov.txt" "${NAME}.SingleCheck.txt"

# restrict to --stages
declare -a RUN_KEYS=()
if [ -n "$STAGES_CSV" ]; then
  IFS=',' read -r -a want <<< "$STAGES_CSV"
  for k in "${want[@]}"; do
    stage_found=0
    for kk in "${STAGE_KEYS[@]}"; do [ "$k" = "$kk" ] && stage_found=1 && break; done
    [ "$stage_found" -eq 1 ] || { echo "Error: unknown chunk '$k' (see --list)" >&2; exit 1; }
    RUN_KEYS+=("$k")
  done
else
  RUN_KEYS=("${STAGE_KEYS[@]}")
fi

###############################################################################
#                        7. RUN                                               #
###############################################################################
TIMEBIN=""
for cand in /usr/bin/time /bin/time; do [ -x "$cand" ] && TIMEBIN="$cand" && break; done
# only GNU time understands -v (used for CPU% and peak RSS)
[ -n "$TIMEBIN" ] && ! "$TIMEBIN" -v -o /dev/null true 2>/dev/null && TIMEBIN=""

printf "stage\tgroup\trep\tstatus\twall_s\tuser_s\tsys_s\tcpu_pct\tmax_rss_kb\texit\tnote\n" > "$TIMINGS"

INPUT_SIZE=$(du -Lh "$ABS_INPUT" 2>/dev/null | awk '{print $1}')
echo "SingleCheck per-chunk benchmark"
echo "  input     : $INPUT (${INPUT_SIZE:-?}) [$METHOD]"
echo "  reps      : $REPS   warmup: $WARMUP   full-pipeline run: $FULL"
echo "  threads   : $THREADS   window: $WSIZE   delta: $DELTA   depth: $downsampling_depth"
echo "  scratch   : $SCRATCH   sort-compress: $SORTCOMP   gzip: $GZIP_CMD"
echo "  results   : $OUTDIR"
echo

# run_one <key> <rep>  -> sets R_WALL R_RC ; appends a row to $TIMINGS
run_one() {
  local key="$1" rep="$2" f
  local tf="$OUTDIR/logs/${key}.rep${rep}.time" lf="$OUTDIR/logs/${key}.rep${rep}.log"

  # every RESET path is absolute and inside $WORK -- refuse anything else
  for f in ${S_RESET[$key]}; do
    case "$f" in "$WORK"/*) rm -f "$f" ;; *) echo "internal: refusing to reset '$f'" >&2 ;; esac
  done

  local start end rc
  start=$(date +%s.%N)
  if [ -n "$TIMEBIN" ]; then
    "$TIMEBIN" -v -o "$tf" bash "$SELF" __stage "$key" > "$lf" 2>&1; rc=$?
  else
    bash "$SELF" __stage "$key" > "$lf" 2>&1; rc=$?
  fi
  end=$(date +%s.%N)

  R_WALL=$(awk -v a="$start" -v b="$end" 'BEGIN{printf "%.3f", b-a}')
  R_RC=$rc
  local user=NA sys=NA cpu=NA rss=NA
  if [ -f "$tf" ]; then
    user=$(awk -F': ' '/User time/{print $2}' "$tf"); user="${user:-NA}"
    sys=$(awk -F': ' '/System time/{print $2}' "$tf"); sys="${sys:-NA}"
    cpu=$(awk -F': ' '/Percent of CPU/{gsub("%","",$2); print $2}' "$tf"); cpu="${cpu:-NA}"
    rss=$(awk '/Maximum resident set size/{print $NF}' "$tf"); rss="${rss:-NA}"
  fi
  # A chunk that exits 0 but produced none of the artifacts it declares is a
  # silent failure (a trailing `rm`/`cp` masking the real status, a tool that
  # printed usage and quit, ...). Do not report a time for it.
  local note="-"
  if [ "$rc" -eq 0 ] && [ -n "${S_RESET[$key]}" ]; then
    local produced=0
    for f in ${S_RESET[$key]}; do [ -e "$f" ] && produced=1 && break; done
    if [ "$produced" -eq 0 ]; then
      rc=127
      note="exited 0 but produced none of its outputs -- see logs/${key}.rep${rep}.log"
      echo "      !! $key produced no output despite exit 0 -- see $lf" >&2
    fi
  fi

  local status="OK"; [ "$rc" -ne 0 ] && status="FAILED"
  [ "$rep" = "warmup" ] && status="WARMUP"
  R_RC=$rc
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$key" "${S_GROUP[$key]}" "$rep" "$status" "$R_WALL" "$user" "$sys" "$cpu" "$rss" "$rc" "$note" >> "$TIMINGS"
}

skip_stage() { # key reason
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$1" "${S_GROUP[$1]}" "0" "SKIPPED" "NA" "NA" "NA" "NA" "NA" "NA" "$2" >> "$TIMINGS"
  printf "  [SKIP] %-16s %s\n" "$1" "$2"
}

for key in "${RUN_KEYS[@]}"; do
  # --- prerequisites: tools ---
  miss=""
  for c in ${S_CMD[$key]}; do
    command -v "$c" >/dev/null 2>&1 || [ -x "$c" ] || miss="$miss $c"
  done
  if [ "$key" = "align" ]; then
    command -v bwa-mem2 >/dev/null 2>&1 || command -v bwa >/dev/null 2>&1 || miss="$miss bwa/bwa-mem2"
    [ -n "$REFERENCE" ] || miss="$miss --ref"
  fi
  if [ "$key" = "metaphyler" ]; then
    # MetaPhyler shells out to BLAST (legacy blastall or BLAST+ blastn)
    command -v perl >/dev/null 2>&1 || miss="$miss perl"
    command -v blastall >/dev/null 2>&1 || command -v blastn >/dev/null 2>&1 || miss="$miss blastall/blastn"
  fi
  if [ -n "$miss" ]; then skip_stage "$key" "missing tool(s):$miss"; continue; fi

  # --- prerequisites: input artifacts ---
  missf=""
  for f in ${S_FILE[$key]}; do [ -e "$f" ] || missf="$missf $(basename "$f")"; done
  if [ -n "$missf" ]; then skip_stage "$key" "missing input(s):$missf"; continue; fi

  printf "  %-16s %s\n" "$key" "${S_DESC[$key]}"
  if [ "$WARMUP" -eq 1 ]; then
    run_one "$key" warmup
    printf "      warmup  %8ss\n" "$R_WALL"
  fi
  fails=0
  for rep in $(seq 1 "$REPS"); do
    run_one "$key" "$rep"
    printf "      rep %-3s %8ss  exit=%s\n" "$rep" "$R_WALL" "$R_RC"
    [ "$R_RC" -ne 0 ] && fails=$((fails+1))
  done
  if [ "$fails" -gt 0 ]; then
    echo "      !! $fails/$REPS repetitions failed -- see $OUTDIR/logs/${key}.rep1.log" >&2
  fi
  # pick up the values this chunk produced (genome_length, probability, ...)
  if [ -s "$BENCH_VALS" ]; then source "$BENCH_VALS"; save_state; fi
done

###############################################################################
#        7b. EQUIVALENCE OF THE OPTIMIZED CHUNKS vs THEIR LEGACY FORM         #
#                                                                             #
# Timing an optimization is only half the job: the report must also say the   #
# rewritten step still produces the same numbers, ON THIS DATA.               #
###############################################################################
EQUIV_LINES=()
if [ "$COMPARE" -eq 1 ]; then
  echo
  echo "  verifying the optimized chunks against the legacy ones"

  # autocorrelation input: same set of (depth, shifted depth, length) rows.
  # Row ORDER differs by construction (a hash has no order), so sort both.
  legacy="${NAME}.${DELTA}.shiftedcov.legacy.txt"
  fast="${NAME}.${DELTA}.shiftedcov.fast.txt"
  if [ -s "$legacy" ] && [ -s "$fast" ]; then
    if diff -q <(sort "$legacy") <(sort "$fast") >/dev/null 2>&1; then
      EQUIV_LINES+=("| \`unionbedg_hash\` vs \`unionbedg_sort\` | **identical** | $(wc -l < "$fast") rows |")
      echo "      shiftedcov: IDENTICAL"
    else
      nd=$(diff <(sort "$legacy") <(sort "$fast") | grep -c '^[<>]')
      EQUIV_LINES+=("| \`unionbedg_hash\` vs \`unionbedg_sort\` | **DIFFERS** | $nd differing rows -- see work/ |")
      echo "      shiftedcov: DIFFERS ($nd rows)" >&2
    fi
  fi

  # mapping statistics: two floats, compare with a relative tolerance
  if [ -s "${NAME}.mapstats.txt" ] && [ -s "${NAME}.mapstats.legacy.txt" ]; then
    if awk 'NR==FNR{a1=$1;a2=$2;next}
            {d1=(a1-$1); if(d1<0)d1=-d1; d2=(a2-$2); if(d2<0)d2=-d2
             r1=(a1>0?d1/a1:d1); r2=(a2>0?d2/a2:d2)
             exit (r1<1e-9 && r2<1e-9) ? 0 : 1}' \
           "${NAME}.mapstats.legacy.txt" "${NAME}.mapstats.txt"; then
      EQUIV_LINES+=("| \`mapstats\` vs \`primary_bam_legacy\`+\`idxstats_legacy\` | **identical** | MT % and unmapped % agree to 1e-9 |")
      echo "      mapstats:   IDENTICAL"
    else
      EQUIV_LINES+=("| \`mapstats\` vs \`primary_bam_legacy\`+\`idxstats_legacy\` | **DIFFERS** | legacy: $(cat "${NAME}.mapstats.legacy.txt"), new: $(cat "${NAME}.mapstats.txt") |")
      echo "      mapstats:   DIFFERS" >&2
    fi
  fi
fi

# --- optional end-to-end reference run --------------------------------------
if [ "$FULL" -eq 1 ]; then
  echo
  if [ ! -x "$ROOTDIR/SingleCheck" ]; then
    echo "  [SKIP] full pipeline: $ROOTDIR/SingleCheck not executable"
  else
    echo "  full pipeline (../SingleCheck end-to-end)"
    SC_ARGS=(-w "$WSIZE" -i "$DELTA" -t "$THREADS" -f "$FLAGTOFILTEROUT" -q "$MAPQUAL" -d "$downsampling_depth")
    [ -n "$REFERENCE" ] && SC_ARGS+=(-r "$REFERENCE")
    [ "$DOWNSAMPLE" = "NO" ] && SC_ARGS+=(-N)
    for rep in $(seq 1 "$REPS"); do
      fdir="$OUTDIR/full/rep$rep"; rm -rf "$fdir"; mkdir -p "$fdir"
      ln -sf "$ABS_INPUT" "$fdir/$BASE_IN"
      [ -f "${ABS_INPUT}.bai" ]  && ln -sf "${ABS_INPUT}.bai"  "$fdir/${BASE_IN}.bai"
      [ -f "${ABS_INPUT}.crai" ] && ln -sf "${ABS_INPUT}.crai" "$fdir/${BASE_IN}.crai"
      [ -n "$INPUT2" ] && ln -sf "$ABS_INPUT2" "$fdir/$BASE_IN2"
      cmd=("$ROOTDIR/SingleCheck" "${SC_ARGS[@]}" "$fdir/$BASE_IN")
      [ -n "$INPUT2" ] && cmd+=("$fdir/$BASE_IN2")
      tf="$OUTDIR/logs/FULL.rep${rep}.time" lf="$OUTDIR/logs/FULL.rep${rep}.log"
      start=$(date +%s.%N)
      if [ -n "$TIMEBIN" ]; then "$TIMEBIN" -v -o "$tf" "${cmd[@]}" > "$lf" 2>&1; rc=$?
      else "${cmd[@]}" > "$lf" 2>&1; rc=$?; fi
      end=$(date +%s.%N)
      wall=$(awk -v a="$start" -v b="$end" 'BEGIN{printf "%.3f", b-a}')
      user=NA; sys=NA; cpu=NA; rss=NA
      if [ -f "$tf" ]; then
        user=$(awk -F': ' '/User time/{print $2}' "$tf"); sys=$(awk -F': ' '/System time/{print $2}' "$tf")
        cpu=$(awk -F': ' '/Percent of CPU/{gsub("%","",$2); print $2}' "$tf")
        rss=$(awk '/Maximum resident set size/{print $NF}' "$tf")
      fi
      st="OK"; [ "$rc" -ne 0 ] && st="FAILED"
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "FULL_PIPELINE" "END-TO-END" "$rep" "$st" "$wall" "${user:-NA}" "${sys:-NA}" "${cpu:-NA}" "${rss:-NA}" "$rc" "-" >> "$TIMINGS"
      printf "      rep %-3s %8ss  exit=%s\n" "$rep" "$wall" "$rc"
    done
  fi
fi

###############################################################################
#                        8. REPORT                                            #
###############################################################################
{
  echo "# SingleCheck per-chunk benchmark"
  echo
  echo "| | |"
  echo "|---|---|"
  echo "| date | $(date '+%Y-%m-%d %H:%M:%S %Z') |"
  echo "| host | $(hostname) |"
  echo "| cpus available | $(nproc 2>/dev/null || echo '?') |"
  echo "| input | \`$ABS_INPUT\` (${INPUT_SIZE:-?}, $METHOD) |"
  echo "| repetitions | $REPS $([ "$WARMUP" -eq 1 ] && echo '(+1 warmup, not measured)') |"
  echo "| dependency check | $missing problem(s), $warnings warning(s) |"
  echo "| threads (-t) | $THREADS (mosdepth capped at $MD_THREADS) |"
  echo "| samtools sort memory | $SORTMEM per thread |"
  echo "| window (-w) | $WSIZE |"
  echo "| delta (-i) | $DELTA |"
  echo "| downsampling depth (-d) | $downsampling_depth |"
  echo "| diploid contigs | \`$DIPLOID_REGEX\` |"
  echo "| scratch (\$TMPDIR) | \`$SCRATCH\` |"
  echo "| sort compressor | $SORTCOMP |"
  echo "| track compressor | \`$GZIP_CMD\` |"
  [ -n "${genome_length:-}" ]    && echo "| genome length | $genome_length |"
  [ -n "${raw_reads:-}" ]        && echo "| primary reads | $raw_reads |"
  [ -n "${mean_readlength:-}" ]  && echo "| mean read length | $mean_readlength |"
  [ -n "${sequencing_depth:-}" ] && echo "| raw sequencing depth | ${sequencing_depth}X |"
  [ -n "${probability:-}" ]      && echo "| downsampling fraction | $probability (mode: ${DS_MODE:-?}) |"
  echo
  echo "Tool versions"
  echo
  echo '```'
  command -v samtools >/dev/null 2>&1 && samtools --version 2>/dev/null | head -n1
  command -v "$MOSDEPTH" >/dev/null 2>&1 || [ -x "$MOSDEPTH" ] && echo "mosdepth $("$MOSDEPTH" --version 2>&1 | head -n1)"
  command -v "$BEDTOOLS" >/dev/null 2>&1 || [ -x "$BEDTOOLS" ] && "$BEDTOOLS" --version 2>/dev/null | head -n1
  command -v Rscript >/dev/null 2>&1 && Rscript --version 2>&1 | head -n1
  command -v bwa-mem2 >/dev/null 2>&1 && bwa-mem2 version 2>&1 | tail -n1
  echo '```'
  echo
} > "$REPORT"

awk -F'\t' -v OFS='\t' '
  NR==1 { next }
  {
    key=$1; grp=$2; rep=$3; st=$4; wall=$5; user=$6; sys=$7; cpu=$8; rss=$9
    if (!(key in seen)) { seen[key]=1; order[++nkeys]=key; group[key]=grp }
    if (st=="SKIPPED") { skipped[key]=$11; next }
    if (st=="WARMUP")  { warm[key]=wall; next }
    if (st=="FAILED")  { failed[key]++ ; next }
    n[key]++
    v[key,n[key]]=wall+0
    sum[key]+=wall+0
    if (n[key]==1 || wall+0<min[key]) min[key]=wall+0
    if (n[key]==1 || wall+0>max[key]) max[key]=wall+0
    if (user!="NA") { usum[key]+=user+0; ucnt[key]++ }
    if (sys!="NA")  { ssum[key]+=sys+0 }
    if (cpu!="NA" && cpu+0>cpumax[key]) cpumax[key]=cpu+0
    if (rss!="NA" && rss+0>rssmax[key]) rssmax[key]=rss+0
  }
  function median(k,  i,j,t,a,c) {
    c=n[k]; for(i=1;i<=c;i++) a[i]=v[k,i]
    for(i=2;i<=c;i++){ t=a[i]; for(j=i-1;j>=1 && a[j]>t;j--) a[j+1]=a[j]; a[j+1]=t }
    return (c%2) ? a[(c+1)/2] : (a[c/2]+a[c/2+1])/2
  }
  function sd(k,  i,m,s) {
    if (n[k]<2) return 0
    m=sum[k]/n[k]; for(i=1;i<=n[k];i++) s+=(v[k,i]-m)^2
    return sqrt(s/(n[k]-1))
  }
  function bar(p,  i,s,w) { w=int(p/2+0.5); s=""; for(i=0;i<w;i++) s=s "#"; return s }
  function fmt(x) { return (x>=100) ? sprintf("%.1f",x) : sprintf("%.3f",x) }
  END {
    total=0
    for(i=1;i<=nkeys;i++){ k=order[i]; if(k!="FULL_PIPELINE" && n[k]>0) total+=sum[k]/n[k] }

    print "## Per-chunk timings (pipeline order)"
    print ""
    print "| # | group | chunk | reps | mean s | min s | median s | max s | sd s | share | peak CPU % | peak RSS MB |"
    print "|---|---|---|---|---|---|---|---|---|---|---|---|"
    idx=0
    for(i=1;i<=nkeys;i++){
      k=order[i]; if(k=="FULL_PIPELINE") continue
      idx++
      if (n[k]==0) {
        printf("| %d | %s | `%s` | 0 | - | - | - | - | - | - | - | **%s** |\n", idx, group[k], k,
               (k in skipped) ? "SKIPPED: " skipped[k] : "FAILED")
        continue
      }
      mean=sum[k]/n[k]
      printf("| %d | %s | `%s` | %d | %s | %s | %s | %s | %s | %.1f%% | %d | %.1f |\n",
             idx, group[k], k, n[k], fmt(mean), fmt(min[k]), fmt(median(k)), fmt(max[k]), fmt(sd(k)),
             (total>0?100*mean/total:0), cpumax[k], rssmax[k]/1024)
    }
    print ""
    printf("**Sum of chunk means: %s s** (%.1f min)\n", fmt(total), total/60)
    print ""
    if ("FULL_PIPELINE" in seen && n["FULL_PIPELINE"]>0) {
      fm=sum["FULL_PIPELINE"]/n["FULL_PIPELINE"]
      print "## End-to-end reference"
      print ""
      printf("| metric | value |\n|---|---|\n")
      printf("| full `SingleCheck` mean | %s s (%.1f min, %d reps) |\n", fmt(fm), fm/60, n["FULL_PIPELINE"])
      printf("| full `SingleCheck` min | %s s |\n", fmt(min["FULL_PIPELINE"]))
      printf("| sum of chunks | %s s |\n", fmt(total))
      printf("| chunks / end-to-end | %.1f%% |\n", (fm>0?100*total/fm:0))
      printf("| unaccounted (deps check, cleanup, startup) | %s s |\n", fmt(fm-total))
      print ""
    }

    print "## Ranked by mean wall-clock (where the time goes)"
    print ""
    print "```"
    # insertion sort keys by mean desc
    m=0
    for(i=1;i<=nkeys;i++){ k=order[i]; if(k!="FULL_PIPELINE" && n[k]>0){ rk[++m]=k } }
    for(i=2;i<=m;i++){ t=rk[i]; tm=sum[t]/n[t]
      for(j=i-1;j>=1 && (sum[rk[j]]/n[rk[j]])<tm;j--) rk[j+1]=rk[j]
      rk[j+1]=t }
    cum=0
    printf("%-16s %10s %7s %7s  %s\n", "CHUNK", "MEAN(s)", "SHARE", "CUMUL", "")
    for(i=1;i<=m;i++){
      k=rk[i]; mean=sum[k]/n[k]; pc=(total>0?100*mean/total:0); cum+=pc
      printf("%-16s %10s %6.1f%% %6.1f%%  %s\n", k, fmt(mean), pc, cum, bar(pc))
    }
    print "```"
    print ""

    print "## By group"
    print ""
    print "| group | mean s | share |"
    print "|---|---|---|"
    for(i=1;i<=nkeys;i++){ k=order[i]; if(k=="FULL_PIPELINE"||n[k]==0) continue
      g=group[k]; if(!(g in gseen)){ gseen[g]=1; gorder[++ng]=g }
      gsum[g]+=sum[k]/n[k] }
    for(i=1;i<=ng;i++){ g=gorder[i]
      printf("| %s | %s | %.1f%% |\n", g, fmt(gsum[g]), (total>0?100*gsum[g]/total:0)) }
    print ""

    nsk=0; for(k in skipped) nsk++
    nfa=0; for(k in failed) nfa++
    if (nsk>0 || nfa>0) {
      print "## Skipped / failed chunks"
      print ""
      print "| chunk | reason |"
      print "|---|---|"
      for(i=1;i<=nkeys;i++){ k=order[i]
        if (k in skipped) printf("| `%s` | skipped -- %s |\n", k, skipped[k])
        else if (k in failed) printf("| `%s` | %d repetition(s) exited non-zero -- see logs/%s.rep1.log |\n", k, failed[k], k) }
      print ""
    }
  }
' "$TIMINGS" >> "$REPORT"

{
  if [ "${#EQUIV_LINES[@]}" -gt 0 ]; then
    echo "## Optimized vs legacy: do they agree?"
    echo
    echo "The chunks rewritten for HPC_OPTIMIZATION.md were run in BOTH forms on this input."
    echo "Compare their rows in the timing table above for the speedup, and this table for"
    echo "correctness:"
    echo
    echo "| comparison | verdict | detail |"
    echo "|---|---|---|"
    printf '%s\n' "${EQUIV_LINES[@]}"
    echo
  fi
  echo "## How to read this"
  echo
  echo "* Chunks are the command chains of \`../SingleCheck\`, executed in pipeline order in"
  echo "  \`work/\`; the artifacts of one chunk feed the next, so the timings are those of a"
  echo "  real run and not of an artificial micro-benchmark."
  echo "* **min** is the fairest single number for comparing chunks: the page cache warms up"
  echo "  across repetitions, so rep 1 is usually the slowest. **sd** tells you how noisy the"
  echo "  measurement was -- if sd is close to mean, raise \`--reps\`."
  echo "* Every repetition of a chunk deletes its own outputs first, so no repetition gets a"
  echo "  free ride on a file another one already produced."
  echo "* **peak CPU %** > 100 means the chunk really used the \`-t $THREADS\` threads; a chunk"
  echo "  stuck near 100 % is single-threaded and a candidate for parallelisation."
  echo "* Chunks that read the *original* input (\`count_reads\`, \`mean_readlen\`, \`downsample\`,"
  echo "  \`unmapped_fasta\`) are the ones that scale with input size; everything after"
  echo "  \`mosdepth\` works on the small downsampled BAM. See \`../HPC_OPTIMIZATION.md\` §1."
  echo "* Per-repetition raw numbers: \`timings.tsv\`. Logs and \`/usr/bin/time -v\` output per"
  echo "  repetition: \`logs/\`."
  [ -z "$TIMEBIN" ] && echo "* \`/usr/bin/time\` was not available: CPU % and RSS columns are empty."
  echo
} >> "$REPORT"

if [ "$CLEAN" -eq 1 ]; then
  echo "Removing working directory $WORK"
  rm -rf "$WORK"
fi

echo
echo "================================================================"
awk '/^## Ranked by mean/{f=1} f && /^```$/{c++; next} f && c==1 {print} c==2{exit}' "$REPORT"
echo "================================================================"
echo "Report        : $REPORT"
echo "Per-rep TSV   : $TIMINGS"
echo "Logs          : $OUTDIR/logs"
[ "$CLEAN" -eq 0 ] && echo "Artifacts     : $WORK"
exit 0
