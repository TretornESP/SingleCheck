# SingleCheck Pipeline — Detailed Step-by-Step Documentation

SingleCheck quantifies **coverage dispersion** in single-cell DNA sequencing libraries. It takes aligned reads (BAM) or raw reads (FASTQ) as input and produces four complementary statistics — Gini coefficient, Coefficient of Variation (CV), Autocorrelation, and Median Absolute Deviation (MAD) — that together characterise how uniformly the genome is covered. Non-uniform coverage is the hallmark of whole-genome amplification (WGA) artefacts present in single-cell DNA-seq data.

---

## Table of Contents

1. [Dependencies](#1-dependencies)
2. [Repository layout](#2-repository-layout)
3. [Command-line interface](#3-command-line-interface)
4. [Full pipeline walkthrough](#4-full-pipeline-walkthrough)
   - [Phase 0 — Argument parsing and variable initialisation](#phase-0--argument-parsing-and-variable-initialisation)
   - [Phase 1 — Module loading](#phase-1--module-loading)
   - [Phase 2 — Read alignment (FASTQ input only)](#phase-2--read-alignment-fastq-input-only)
   - [Phase 3 — Raw sequencing depth calculation](#phase-3--raw-sequencing-depth-calculation)
   - [Phase 4 — Downsampling](#phase-4--downsampling)
   - [Phase 5 — Per-window depth with mosdepth](#phase-5--per-window-depth-with-mosdepth)
   - [Phase 6 — Autocorrelation input preparation](#phase-6--autocorrelation-input-preparation)
   - [Phase 7 — Gini coefficient and CV](#phase-7--gini-coefficient-and-cv)
   - [Phase 8 — Autocorrelation](#phase-8--autocorrelation)
   - [Phase 9 — MAD (Median Absolute Deviation)](#phase-9--mad-median-absolute-deviation)
   - [Phase 10 — Contamination detection with Metaphyler](#phase-10--contamination-detection-with-metaphyler)
   - [Phase 11 — Quality metrics from downsampled BAM](#phase-11--quality-metrics-from-downsampled-bam)
   - [Phase 12 — Output generation and cleanup](#phase-12--output-generation-and-cleanup)
5. [R scripts in detail](#5-r-scripts-in-detail)
   - [GiniIndex.R](#giniindexr)
   - [CoefficientOfVariation.R](#coefficientofvariationr)
   - [Autocorrelation.R](#autocorrelationr)
   - [MAD.R](#madr)
6. [Aggregating results — CreateInputApp](#6-aggregating-results--createinputapp)
7. [Output format and interpretation](#7-output-format-and-interpretation)
8. [Helper and batch scripts](#8-helper-and-batch-scripts)
9. [Shiny app](#9-shiny-app)
10. [End-to-end usage examples](#10-end-to-end-usage-examples)

---

## 1. Dependencies

| Tool | Version tested | Purpose |
|------|---------------|---------|
| BWA-MEM | 0.7.17 | Align FASTQ reads to a reference genome |
| samtools | 1.10 | BAM manipulation, indexing, statistics |
| Picard | 2.18.14 | Downsampling (`DownsampleSam`) |
| mosdepth | any recent | Fast per-base and per-window depth calculation |
| bedtools | 2.28.0 | Merging shifted coverage tracks |
| Metaphyler | SRV0.115 | Microbial contamination detection from unmapped reads |
| R | 3.6.3+ | Statistical metric computation |
| R packages (main) | — | `tidyr`, `dplyr`, `matrixStats` |
| R packages (app) | — | `shiny`, `shinydashboard`, `shinycssloaders`, `DT`, `data.table`, `ggplot2`, `plotly` |

The script expects a SLURM HPC environment with `module load` commands, but it also works outside SLURM — in that case `module load` calls silently fail and the tools must be on `$PATH`.

---

## 2. Repository layout

```
SingleCheck/
├── SingleCheck              ← Main pipeline script (SLURM-ready bash)
├── CreateInputApp           ← Aggregates per-sample results into one table
├── README.md
├── SOURCES.md               ← Internal development notes (Spanish)
├── Workflow-SingleCheck.png ← Visual pipeline diagram
├── src/                     ← R scripts called by the main pipeline
│   ├── GiniIndex.R
│   ├── CoefficientOfVariation.R
│   ├── Autocorrelation.R
│   └── MAD.R
├── analysis/                ← Batch / helper scripts
│   ├── SingleCheckArray     ← SLURM array job wrapper
│   ├── SortSam.sh
│   ├── DownsampleSam.sh
│   ├── GetSeqDepth.sh
│   ├── GetSeqDepthFixingReadLength.sh
│   ├── CheckSeqDepth.sh
│   ├── CollectInsertSizeMetrics.sh
│   ├── SimulateSingleCellReads.sh
│   ├── RunExample.sh
│   ├── RunWANG.sh / RunHUANG.sh
│   └── CreatePlotsWANG.sh / CreatePlotsHUANG.sh
├── simulations/             ← R scripts for testing metric behaviour
│   ├── SmallSimulations.R
│   ├── CoefficientOfVariation.R
│   ├── GiniIndex.R
│   ├── MAD.R
│   └── Autocorrelation.R
└── ShinyApp/                ← Interactive results dashboard
    ├── ui.R
    └── server.R
```

---

## 3. Command-line interface

```
SingleCheck [options] <in.bam|in.fastq.gz> [in2.fastq.gz]
```

### Positional arguments

| Position | Value | Description |
|----------|-------|-------------|
| 1 | `<in.bam>` | Pre-aligned BAM (must be indexed with `.bam.bai`) |
| 1 | `<in_1.fastq.gz>` | First mate of paired-end reads (must end `_1.fastq.gz`) |
| 2 | `<in_2.fastq.gz>` | Second mate (must end `_2.fastq.gz`) |

### Options

| Flag | Default | Meaning |
|------|---------|---------|
| `-h` | — | Print help and exit |
| `-w <int>` | `10000000` | Window size in bp used for Gini, CV, and MAD |
| `-i <int>` | `1000` | Lag (Delta) in bp used for autocorrelation |
| `-t <int>` | `3` | Number of threads for parallelisable steps |
| `-f <int>` | `772` | SAM flag — reads matching this flag are **excluded** from mosdepth (772 = unmapped + non-primary + QC-failed) |
| `-q <int>` | `20` | Minimum mapping quality for mosdepth |
| `-r <ref.fa>` | — | Reference FASTA (required when input is FASTQ; must be BWA- and samtools-indexed) |
| `-X` | off | Include chromosome X in the diploid set (human samples only) |
| `-N` | off | Skip downsampling; compute statistics from the original file |
| `-d <float>` | `0.1` | Target downsampling depth in X coverage |
| `-s <strategy>` | `ConstantMemory` | Picard DownsampleSam strategy: `ConstantMemory`, `HighAccuracy`, or `Chained` |
| `-c <REGEX>` | `^(chr)*[1-9]` | Custom regex to select chromosomes to include (overrides `-X`) |
| `-m <STRING>` | `^(MT\|chrM)` | Name of the mitochondrial contig if non-standard |

---

## 4. Full pipeline walkthrough

### Phase 0 — Argument parsing and variable initialisation

**File:** `SingleCheck`, lines 47–131

Default values are set before `getopts` processes the command line:

```bash
WSIZE=10000000          # Window size for Gini/CV/MAD
DELTA=1000              # Lag for autocorrelation
FLAGTOFILTEROUT=772     # SAM flag exclusion mask
THREADS=3
MAPQUAL=20
DIPLOID_REGEX="^(chr)*[1-9]"   # Regex for autosomal chromosomes
MT_REGEX='^(MT|chrM)'          # Regex for mitochondrial contig
DOWNSAMPLE="YES"               # Perform downsampling by default
downsampling_depth=0.1         # Target depth
ds_strategy="ConstantMemory"   # Picard strategy
DOWNSAMPLING="ImpreciseSeqDepthCalc"   # Which depth-calc method to use
```

`getopts` parses each flag and validates its value (numeric checks, strategy membership, etc.). After parsing, the remaining positional arguments determine the **input method**:

- **1 argument ending in `.bam`** → `METHOD="Aligned"`. The sample name is the filename without `.bam`.
- **1 argument ending in `.fastq.gz`** → `METHOD="Single-end"`. The sample name is the filename without `.fastq.gz`.
- **2 arguments** (both `.fastq.gz`) → `METHOD="Paired-end"`. Files must match `*_1.fastq.gz` / `*_2.fastq.gz`. The sample name is stripped of `_1.fastq.gz`.

The script directory (`SCRIPTDIR`) is resolved via `scontrol` under SLURM, or `dirname $0` otherwise.

---

### Phase 1 — Module loading

**File:** `SingleCheck`, lines 210–218

```bash
module purge
module load samtools/1.10
module load gcc/6.4.0 R/3.6.3
module load miniconda3/4.8.2
module load gcccore/6.4.0 bedtools/2.28.0
module load gcc/6.4.0 bwa/0.7.17
module load picard/2.18.14
source activate /mnt/netapp1/posadalab/APPS/CommonCondaEnvironments/mosdepth
```

All `2> /dev/null` redirects suppress errors so the script works outside SLURM (tools must already be on `$PATH`). The `mosdepth` conda environment is activated via `source activate`.

---

### Phase 2 — Read alignment (FASTQ input only)

**File:** `SingleCheck`, lines 226–240

Runs only when `METHOD` is `Paired-end` or `Single-end`. First, the reference is validated: both the `.fai` samtools index and the `.fa.ann` BWA index must exist.

```bash
bwa mem \
    -t $THREADS \
    $REFERENCE \
    $FASTQ1 $FASTQ2 | \
    samtools sort -@$THREADS -o ${NAME}.bam -

samtools index ${NAME}.bam
```

- `bwa mem -t $THREADS $REFERENCE $FASTQ1 $FASTQ2` — Runs BWA-MEM aligner. `$FASTQ2` is empty for single-end input, so BWA runs in single-end mode.
- The SAM output is piped directly into `samtools sort`, which writes a coordinate-sorted BAM to `${NAME}.bam`. Sorting is done in-memory using `$THREADS` threads.
- `samtools index ${NAME}.bam` creates `${NAME}.bam.bai`, required by all downstream tools.

---

### Phase 3 — Raw sequencing depth calculation

**File:** `SingleCheck`, lines 247–291

This phase computes the total number of sequenced bases and the resulting mean genome coverage (sequencing depth) before any downsampling.

#### Step 3a — Validate BAM index

```bash
if [ ! -f ${NAME}.bam.bai ]; then
    printf "\nError: your BAM file must be indexed\n" >&2
    exit 1
fi
```

#### Step 3b — Genome length from BAM header

```bash
genome_length=$(samtools view -H ${NAME}.bam | \
    grep "^@SQ" | \
    awk '{sum+=substr($3,4,length($3)-1)}END{print sum}')
```

- `samtools view -H` prints only the header lines.
- `grep "^@SQ"` isolates sequence dictionary lines (one per reference contig).
- The `awk` expression extracts the `LN:` field (e.g. `LN:249250621`), strips the 4-character `LN:` prefix, and accumulates the total genome length in bp.

#### Step 3c — Count primary reads and calculate raw bases (fast method)

The default method (`ImpreciseSeqDepthCalc`) cannot count soft/hard-clipped bases but is much faster than the precise alternative.

```bash
# Remove supplementary (flag 2048) and secondary (flag 256) alignments
# to avoid counting the same read twice
samtools view -bF 2304 ${NAME}.bam > ${NAME}.primary.bam
samtools index ${NAME}.primary.bam
```

Flag `2304 = 2048 + 256` (supplementary + secondary). Filtering leaves only primary alignments.

```bash
raw_reads=$(samtools idxstats ${NAME}.primary.bam | \
    awk '{sum+=($3+$4)}END{print sum}')
```

`samtools idxstats` outputs one line per contig with: name, length, mapped reads, unmapped reads. Summing columns 3 and 4 gives the total number of read segments.

```bash
mean_readlength=$(samtools view ${NAME}.primary.bam | \
    head -n 1000000 | \
    cut -f10 | \
    awk '{print length}' | \
    sort | uniq -c | \
    awk '{sum+=$1*$2; num+=$1}END{print sum/num}')
```

Reads the SEQ field (column 10) from up to 1 million primary alignments, computes a frequency table of sequence lengths, and calculates the weighted mean read length.

```bash
raw_bases=$(awk -v meanl=$mean_readlength -v rawr=$raw_reads \
    'BEGIN{print rawr*meanl}')
```

Total sequenced bases ≈ total reads × mean read length.

#### Step 3d — Precise method (optional, not default)

When `DOWNSAMPLING` is not `ImpreciseSeqDepthCalc`, the script counts bases from the CIGAR string:

```bash
# Count bases from SEQ field (aligned + soft-clipped + unmapped)
aligned_soft_bases=$(samtools view ${NAME}.primary.bam | \
    cut -f10 | awk '{total+=length}END{print total}')

# Count hard-clipped bases from CIGAR (column 6)
hard_bases=$(samtools view ${NAME}.primary.bam | \
    cut -f6 | grep H | \
    sed 's/\([0-9]*\)\([A-Z]\)/\1\2\n/g' | \
    grep -v "^$" | grep H | sed 's/H//' | \
    awk '{sum+=$1}END{print sum}')

raw_bases=$((aligned_soft_bases + hard_bases))
```

Hard-clipped bases are not present in the SEQ field, so they must be extracted from the CIGAR string separately.

#### Step 3e — Sequencing depth

```bash
sequencing_depth=$(awk -v gl=$genome_length -v rawb=$raw_bases \
    'BEGIN{print rawb/gl}')
```

Formula: **depth = raw_bases / genome_length**

The temporary primary BAM and its index are then removed.

---

### Phase 4 — Downsampling

**File:** `SingleCheck`, lines 298–342

#### Step 4a — Calculate downsampling probability

```bash
probability=`bc -l <<< "scale=10; $downsampling_depth / $sequencing_depth"`
```

`bc -l` with `scale=10` gives 10 decimal places of precision. The probability is the fraction of reads to retain: **P = target_depth / actual_depth**.

#### Step 4b — Decide whether to downsample

```bash
if [[ ! -z $(awk -v prob=$probability \
    'BEGIN{if (prob > 1) print "Lowest sequencing than downsampling selected"}') ]] \
    || [[ "$DOWNSAMPLE" = "NO" ]]
then
    ln -s ${NAME}.bam ${NAME}.${downsampling_depth}X.bam
    samtools index ${NAME}.${downsampling_depth}X.bam
    DEPTH=$sequencing_depth
```

If the sample is already at or below the target depth, or `-N` was passed, a symlink is created instead of actually downsampling. `DEPTH` is set to the actual sequencing depth in this case.

#### Step 4c — Picard DownsampleSam

When downsampling is needed:

```bash
# Calculate accuracy from the exponent of the probability
num=$(printf %e ${probability} | fold -w1 | tail -n 1)
accuracy=$(awk -v num=$num 'BEGIN{print 0.01*10^-num}')
```

- `printf %e` converts the probability to scientific notation (e.g. `1.23456789e-02`).
- `fold -w1 | tail -n 1` extracts the last character, which is the exponent digit.
- `accuracy = 0.01 × 10^(-exponent)` — a tighter accuracy is used for very small probabilities to ensure the correct number of reads is retained.

```bash
java -jar $EBROOTPICARD/picard.jar DownsampleSam \
    INPUT=${NAME}.bam \
    OUTPUT=${NAME}.${downsampling_depth}X.bam \
    RANDOM_SEED=1 \
    PROBABILITY=${probability} \
    STRATEGY=$ds_strategy \
    CREATE_INDEX=true \
    ACCURACY=$accuracy
```

| Parameter | Description |
|-----------|-------------|
| `RANDOM_SEED=1` | Reproducible downsampling |
| `PROBABILITY` | Fraction of read-templates to keep |
| `STRATEGY` | Algorithm: `ConstantMemory` (fast, low memory), `HighAccuracy` (slow, high memory), or `Chained` (intermediate) |
| `CREATE_INDEX=true` | Produces `.bai` index automatically |
| `ACCURACY` | Acceptable deviation from the requested probability |

Read-pairs (and their supplementary/secondary alignments) are kept or discarded as a unit, preserving pairing integrity.

`DEPTH` is set to `$downsampling_depth` after this step.

---

### Phase 5 — Per-window depth with mosdepth

**File:** `SingleCheck`, lines 325–343

#### Step 5a — Calculate mosdepth precision

```bash
precision=$(echo $WSIZE | wc -c | awk '{print $1-2}')
```

`wc -c` counts characters including the newline. Subtracting 2 gives the number of digits minus 1, i.e. the number of decimal places needed to represent the minimum possible coverage for a window of that size. For example, a window of 100 bp (3 chars) gives precision = 1; for 10 000 000 bp (8 chars) gives precision = 6. This sets `MOSDEPTH_PRECISION` to avoid rounding artefacts.

#### Step 5b — Run mosdepth

```bash
MOSDEPTH_PRECISION=${precision} mosdepth \
    -t $THREADS \
    --fast-mode \
    --by $WSIZE \
    --flag $FLAGTOFILTEROUT \
    --mapq $MAPQUAL \
    ${NAME}.${WSIZE} \
    ${NAME}.${downsampling_depth}X.bam
```

| Flag | Description |
|------|-------------|
| `MOSDEPTH_PRECISION` | Number of decimal places in output depth values |
| `-t $THREADS` | Parallelise BAM decompression |
| `--fast-mode` | Skip CIGAR-operation correction and mate-overlap correction (suitable for most use cases; significantly faster) |
| `--by $WSIZE` | Tile the genome into non-overlapping windows of `$WSIZE` bp and report mean depth per window |
| `--flag $FLAGTOFILTEROUT` | Exclude reads with any of these SAM flags set (default 772 = unmapped + non-primary + QC-failed) |
| `--mapq $MAPQUAL` | Exclude reads with MAPQ below this threshold |
| `${NAME}.${WSIZE}` | Output prefix |
| input BAM | The downsampled (or symlinked) BAM |

**Output files produced by mosdepth:**

| File | Content |
|------|---------|
| `${NAME}.${WSIZE}.regions.bed.gz` | Mean depth for each `$WSIZE`-bp window |
| `${NAME}.${WSIZE}.per-base.bed.gz` | Per-base depth across the whole genome |
| `${NAME}.${WSIZE}.mosdepth.*` | Summary files — removed immediately after |

The `.regions.bed.gz` file has columns: `chrom start end mean_depth`.
The `.per-base.bed.gz` file has columns: `chrom start end depth`.

After mosdepth finishes, the mosdepth summary files are removed and the conda environment is deactivated:

```bash
rm ${NAME}.${WSIZE}.mosdepth*
conda deactivate
```

---

### Phase 6 — Autocorrelation input preparation

**File:** `SingleCheck`, lines 345–362

Autocorrelation measures whether coverage at position *x* predicts coverage at position *x + Δ*. To compute it, we need paired (depth, depth_shifted_by_Δ) values for every base pair.

#### Step 6a — Shift the per-base coverage track by Δ

```bash
zcat ${NAME}.${WSIZE}.per-base.bed.gz | \
awk -v alpha=$DELTA -v print_switch=0 \
'{
    if (print_switch==1) {
        start=$2-alpha;
        if (start<0) {print $1"\t"0"\t"$3-alpha"\t"$4}
        else         {print $1"\t"start"\t"$3-alpha"\t"$4}
    }
    else if (alpha >= $2 && alpha < $3) {
        start=$2-alpha;
        if (start<0) {print $1"\t"0"\t"$3-alpha"\t"$4}
        else         {print $1"\t"start"\t"$3-alpha"\t"$4};
        print_switch=1
    }
}' | \
gzip -c > ${NAME}.${WSIZE}.${DELTA}.bed.gz
```

This AWK script reads the per-base BED and outputs a new BED where each interval is shifted **backwards** by `$DELTA` bp, effectively making the depth at position *p* appear at position *p − Δ*. The result is saved as `${NAME}.${WSIZE}.${DELTA}.bed.gz`.

- `print_switch` starts at 0 and flips to 1 once the script reaches the position where the lag `alpha` begins.
- Positions that would go negative are clamped to 0.

#### Step 6b — Merge original and shifted tracks

```bash
bedtools unionbedg \
    -filler NA \
    -i ${NAME}.${WSIZE}.per-base.bed.gz \
       ${NAME}.${WSIZE}.${DELTA}.bed.gz | \
grep -E "$DIPLOID_REGEX" | \
sort --parallel=$THREADS --version-sort -k4 -k5 | \
awk '
{
    if (FNR==1) {diff=$3-$2; value1=$4; value2=$5}
    else if (value1!=$4 || value2!=$5) {
        print value1"\t"value2"\t"diff;
        diff=$3-$2; value1=$4; value2=$5
    }
    else {diff=diff+($3-$2); value1=$4; value2=$5}
}
END{print value1"\t"value2"\t"diff}' \
> ${NAME}.${DELTA}.shiftedcov.txt
```

- `bedtools unionbedg -filler NA` merges the two BED files column-wise. Where one track has no data, the value is `NA`.
- `grep -E "$DIPLOID_REGEX"` retains only chromosomes matching the diploid regex (autosomal by default).
- `sort --version-sort -k4 -k5` sorts by the depth pair (columns 4 and 5) for run-length encoding.
- The `awk` block performs **run-length encoding**: consecutive rows with identical (depth, depth_shifted) pairs are merged into a single row with a summed base-pair count.

**Output:** `${NAME}.${DELTA}.shiftedcov.txt` — three tab-separated columns:

```
depth_original   depth_shifted   count_bp
5                6               1000
NA               10              500
7                NA              300
...
```

`NA` appears in `depth_shifted` for the last `$DELTA` bp of each chromosome (no future position to pair with), and in `depth_original` for the first `$DELTA` bp (no past position to shift from).

This file is also used later for the **breadth of coverage** calculation.

---

### Phase 7 — Gini coefficient and CV

**File:** `SingleCheck`, lines 364–375

#### Step 7a — Build the depth frequency table

```bash
zcat ${NAME}.${WSIZE}.regions.bed.gz | \
grep -E "$DIPLOID_REGEX" | \
awk '{print $4}' | \
sort -n | \
uniq -c | \
awk '{print $2"\t"$1}' | \
sort -n -k1,1 \
> ${NAME}.${WSIZE}.freqs.txt
```

- Reads the window-level mean depths from mosdepth output.
- Keeps only diploid chromosomes.
- Extracts the depth column (`$4`).
- `sort -n | uniq -c` produces a frequency table.
- The two `awk` commands reformat to `depth<TAB>count` and sort numerically by depth.

**Output:** `${NAME}.${WSIZE}.freqs.txt`

```
0.12    150
0.15    230
0.17    410
...
```

Each row: mean depth value observed in some windows, count of windows with that exact depth.

#### Step 7b — Call R scripts

```bash
Rscript ${SCRIPTDIR}/src/GiniIndex.R ${NAME}.${WSIZE} $WSIZE
Rscript ${SCRIPTDIR}/src/CoefficientOfVariation.R ${NAME}.${WSIZE}
```

Both receive the file basename (without `.freqs.txt`) and read the frequency table from it. See [Section 5](#5-r-scripts-in-detail) for the mathematical details.

---

### Phase 8 — Autocorrelation

**File:** `SingleCheck`, line 379

```bash
Rscript ${SCRIPTDIR}/src/Autocorrelation.R ${NAME}.${DELTA} $DELTA
```

Receives the basename for the shifted coverage file (without `.shiftedcov.txt`) and the lag value `$DELTA`. See [Section 5](#5-r-scripts-in-detail).

---

### Phase 9 — MAD (Median Absolute Deviation)

**File:** `SingleCheck`, lines 382–397

#### Step 9a — Build the contiguous-window comparison table

```bash
zcat ${NAME}.${WSIZE}.regions.bed.gz | \
    tail -n +2 | \
    paste <(zcat ${NAME}.${WSIZE}.regions.bed.gz) - | \
    grep -E "$DIPLOID_REGEX" | \
    awk '{
        if (NF==8 && $1==$5) {print $4"\t"$8}
        else if (NF==8 && $1!=$5) {print $4}
    }' | \
    sort | uniq -c | \
    awk '{
        if (NF==3) {print $2"\t"$3"\t"$1}
        else       {print $2"\tNA\t"$1}
    }' \
    > ${NAME}.${WSIZE}.contiguous.txt
```

- `tail -n +2 | paste ... -` creates a two-column view by pairing each line with the **next** line. The result has 8 fields: chr1, start1, end1, depth1, chr2, start2, end2, depth2.
- `grep -E "$DIPLOID_REGEX"` keeps only autosomal windows.
- The `awk` block: if both windows are on the **same chromosome** (`$1==$5`), prints `depth1<TAB>depth2`; if they are on **different chromosomes** (window boundary), prints only `depth1` (the cross-chromosome transition is discarded).
- `sort | uniq -c` run-length encodes the depth pairs.
- Final `awk` reformats to `depth1<TAB>depth2<TAB>count`, using `NA` for unpaired values.

**Output:** `${NAME}.${WSIZE}.contiguous.txt`

```
0.15    0.17    200
0.20    NA      50
0.22    0.21    300
...
```

`NA` in column 2 occurs for the last window on each chromosome (no next window on the same chromosome).

#### Step 9b — Call R script

```bash
Rscript ${SCRIPTDIR}/src/MAD.R ${NAME}.${WSIZE}
```

See [Section 5](#5-r-scripts-in-detail) for the calculation.

---

### Phase 10 — Contamination detection with Metaphyler

**File:** `SingleCheck`, lines 400–402

Unmapped reads are extracted from the **original** (non-downsampled) BAM, maximising sensitivity for contamination detection.

```bash
# Extract unmapped reads and convert to FASTA
samtools view -f 0x4 ${NAME}.bam | \
    awk '{OFS="\t"; print ">"$1"\n"$10}' \
    > ${NAME}.unmapped.fasta
```

- `samtools view -f 0x4` selects reads where the **read unmapped** flag (0x4) is **set**.
- The AWK command reformats from SAM to FASTA: `>readname` on one line, then the SEQ field on the next.

```bash
~/apps/Metaphyler/MetaPhylerSRV0.115/metaphyler.pl 2 \
    ${NAME}.unmapped.fasta \
    ${NAME}
```

- Metaphyler classifies reads against a database of phylogenetic marker genes.
- The `2` argument selects classification mode 2.
- Output: `${NAME}.genus.tab` — a table with genus-level taxonomy assignments.

**Note:** The Metaphyler path is hardcoded to `~/apps/Metaphyler/`. This must be adjusted for each installation.

---

### Phase 11 — Quality metrics from downsampled BAM

**File:** `SingleCheck`, lines 405–413

A temporary primary-only BAM is created from the downsampled file:

```bash
samtools view -bF 2304 ${NAME}.${downsampling_depth}X.bam \
    > ${NAME}.${downsampling_depth}X.primary.bam
samtools index ${NAME}.${downsampling_depth}X.primary.bam
```

#### Mitochondrial read percentage

```bash
MT_mappedreads=$(samtools idxstats ${NAME}.${downsampling_depth}X.primary.bam | \
    grep -E "$MT_REGEX" | \
    awk '{print $3}')

mt_perc_totalreads=$(samtools idxstats ${NAME}.${downsampling_depth}X.primary.bam | \
    awk -v mt=$MT_mappedreads \
    '{sum+=($3+$4)}END{print mt/sum*100}')
```

`idxstats` gives mapped read counts per contig. The mitochondrial contig is identified with `$MT_REGEX` (default: `^(MT|chrM)`). The percentage is computed as `MT_reads / total_reads × 100`.

#### Unmapped read percentage

```bash
unmapped_perc_totalreads=$(samtools idxstats \
    ${NAME}.${downsampling_depth}X.primary.bam | \
    awk '{mapped+=$3; unmapped+=$4}END{print unmapped/(unmapped+mapped)*100}')
```

#### Breadth of coverage

```bash
breadth=$(awk \
    '{if ($1==0) {sum+=$3}
      else if ($1!="NA") {rest+=$3}}
    END{print 100 - ((sum/(rest+sum))*100)}' \
    ${NAME}.${DELTA}.shiftedcov.txt)
```

The shifted coverage file (from Phase 6) is reused here. Rows where `depth_original == 0` represent uncovered bases. Breadth = 100% − (uncovered_bp / total_bp × 100%).

#### Contaminant classification summary

```bash
class=$(awk '{if ($1 !~ "{") print $0}' ${NAME}.genus.tab | \
    grep -v "^@" | \
    awk '{print $1"-"$2"-"$3"-"$4"-"$5}' | \
    tr -s '\n' ',' | \
    sed 's/,$/\n/')
```

Reads Metaphyler's genus-level output, filters header lines, formats each hit as `field1-field2-field3-field4-field5`, and joins all entries with commas into a single-line string.

#### Extract metric values

```bash
SAMPLE=$(basename $NAME)
WORKDIR=$(dirname "$NAME")

autocorrelation=$(awk '{print $2}' ${WORKDIR}/Autocorrelation.${SAMPLE}.${DELTA}.txt)
gini=$(awk '{print $2}' ${WORKDIR}/Gini.${SAMPLE}.${WSIZE}.txt)
CV=$(awk '{print $2}' ${WORKDIR}/CV.${SAMPLE}.${WSIZE}.txt)
MAD=$(awk '{print $2}' ${WORKDIR}/MAD.${SAMPLE}.${WSIZE}.txt)
```

Each R script writes its result to a two-column text file: `sample_name<TAB>metric_value`. Column 2 is extracted here.

---

### Phase 12 — Output generation and cleanup

**File:** `SingleCheck`, lines 423–445

#### Cleanup

```bash
rm ${NAME}.${downsampling_depth}X.ba*          # Downsampled BAM and index
rm ${NAME}.${downsampling_depth}X.primary.bam  # Temporary primary BAM
rm ${NAME}.${downsampling_depth}X.primary.bam.bai
rm ${NAME}.${DELTA}.shiftedcov.txt             # Autocorrelation input
rm ${NAME}.${WSIZE}.${DELTA}.bed.gz            # Shifted coverage BED
rm ${NAME}.${WSIZE}.contiguous.txt             # MAD input
rm ${NAME}.${WSIZE}.freqs.txt                  # Gini/CV input
rm ${NAME}.${WSIZE}.regions.bed.gz             # mosdepth window output
rm ${NAME}.${WSIZE}.regions.bed.gz.csi
rm ${NAME}.${WSIZE}.per-base.bed.gz            # mosdepth per-base output
rm ${NAME}.${WSIZE}.per-base.bed.gz.csi
rm ${NAME}.unmapped.fasta                      # Metaphyler input
rm ${NAME}*.tab                                # Metaphyler output
rm ${NAME}.map                                 # Metaphyler intermediate
rm ${WORKDIR}/CV.${SAMPLE}.${WSIZE}.txt        # Metric intermediate files
rm ${WORKDIR}/Gini.${SAMPLE}.${WSIZE}.txt
rm ${WORKDIR}/Autocorrelation.${SAMPLE}.${DELTA}.txt
rm ${WORKDIR}/MAD.${SAMPLE}.${WSIZE}.txt
```

Only three files are kept: `${NAME}.bam`, `${NAME}.bam.bai`, and the final results file.

#### Final output

```bash
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$SAMPLE" "$raw_bases" "$DEPTH" "$WSIZE" "$DELTA" \
    "$unmapped_perc_totalreads" "$mt_perc_totalreads" \
    "$breadth" "$autocorrelation" "$CV" "$gini" "$MAD" "$class" \
    > ${NAME}.SingleCheck.txt
```

A single tab-delimited line is written to `${NAME}.SingleCheck.txt`.

---

## 5. R scripts in detail

### GiniIndex.R

**File:** `src/GiniIndex.R`

**Input:** `${basename}.freqs.txt` (depth, count pairs from Phase 7a)

**Arguments:**
- `args[1]` — file path prefix (without `.freqs.txt`)
- `args[2]` — window size in bp (used to scale mosdepth's floating-point depth values back to integer base counts)

**Algorithm:**

```r
genomecov <- read.table(paste(args[1], ".freqs.txt", sep=""))
colnames(genomecov) <- c("depth", "count")

# Scale depths: mosdepth reports mean depth per window as a float;
# multiplying by window_size recovers approximate base counts.
genomecov <- genomecov %>%
    mutate(sumcount = cumsum(as.numeric(count))) %>%
    mutate(depth = round(depth * window_size))

lengthSeq <- sum(genomecov$count)          # total number of windows
meanDepth <- sum(genomecov$depth * genomecov$count) / lengthSeq

# Gini formula (area-based definition):
# G = (2 / (N² × μ)) × Σ_i [ (Σ_{j=0}^{rank_i} j) - (Σ_{j=0}^{rank_i - count_i} j) ] × (depth_i - μ)
first_formula_term <- 2 / (lengthSeq^2 * meanDepth)

genomecov <- genomecov %>% rowwise() %>%
    mutate(first = (sum(0:sumcount) - sum(0:(sumcount - count)))) %>%
    mutate(tosum = first * (depth - meanDepth))

gini <- first_formula_term * sum(genomecov$tosum)
```

**Output:** `Gini.${basename}.txt` — two columns: sample name, Gini value.

**Interpretation:** A Gini coefficient of 0 means perfectly uniform coverage; 1 means all sequencing mapped to a single window. Values for single-cell libraries are typically much higher than for bulk DNA.

---

### CoefficientOfVariation.R

**File:** `src/CoefficientOfVariation.R`

**Input:** Same `${basename}.freqs.txt`

**Algorithm:**

```r
genomecov <- read.table(paste(args[1], ".freqs.txt", sep=""))
colnames(genomecov) <- c("depth", "count")

lengthSeq <- sum(genomecov$count)
meanDepth <- sum(genomecov$depth * genomecov$count) / lengthSeq

# Weighted standard deviation
genomecov <- genomecov %>%
    mutate(tosum = (depth - meanDepth)^2 * count)

upper_term <- sqrt(sum(genomecov$tosum) / (sum(genomecov$count) - 1))
cv <- upper_term / meanDepth    # CV = SD / mean
```

**Output:** `CV.${basename}.txt` — two columns: sample name, CV value.

**Interpretation:** CV is the standard deviation expressed as a fraction of the mean. Higher CV = more variable coverage = more WGA amplification bias.

---

### Autocorrelation.R

**File:** `src/Autocorrelation.R`

**Input:** `${basename}.shiftedcov.txt` (depth pairs from Phase 6b)

**Arguments:**
- `args[1]` — file path prefix (without `.shiftedcov.txt`)
- `args[2]` — lag Delta in bp

**Algorithm:**

```r
genomecov <- read.table(paste(args[1], ".shiftedcov.txt", sep=""))
colnames(genomecov) <- c("depth", "depth_fwd", "count")
alpha <- as.numeric(args[2])

lengthSeq <- sum(genomecov$count)

# Mean depth computed only from rows where depth is not NA
meanDepth <- genomecov %>%
    drop_na("depth") %>%
    mutate(c = depth * count) %>%
    summarize(sum(c) / lengthSeq) %>%
    as.numeric()

# Keep only rows where both depth values are present
genomecov <- genomecov %>%
    drop_na("depth_fwd") %>%
    drop_na("depth") %>%
    mutate(c = depth * depth_fwd * count)

# Normalised autocovariance at lag alpha:
# Autocorr = (E[X_t × X_{t+α}] - μ²) / μ²
first_term  <- sum(genomecov$c) / (lengthSeq - alpha)
second_term <- meanDepth^2
autocorrelation <- (first_term - second_term) / second_term
```

**Output:** `Autocorrelation.${basename}.txt` — two columns: sample name, autocorrelation value.

**Interpretation:** A high positive autocorrelation means that regions with high coverage tend to be near other high-coverage regions (clustered amplification). Uniform coverage → autocorrelation near 0. WGA libraries typically show high positive autocorrelation.

---

### MAD.R

**File:** `src/MAD.R`

**Input:** `${basename}.contiguous.txt` (paired consecutive window depths from Phase 9a)

**Algorithm:**

```r
genomecov <- read.table(paste(args[1], ".contiguous.txt", sep=""))
colnames(genomecov) <- c("depth", "depth_fwd", "count")

lengthSeq <- sum(genomecov$count)
meanDepth <- sum(genomecov$depth * genomecov$count) / lengthSeq

# Remove windows at chromosome boundaries (NA in depth_fwd)
genomecov <- genomecov %>%
    drop_na("depth_fwd") %>%
    mutate(c = (depth - depth_fwd) / meanDepth)   # normalised difference

# Weighted median of differences
median_diffs <- weightedMedian(genomecov$c, genomecov$count, ties="weighted")

# Weighted median of absolute deviations from the median
genomecov <- genomecov %>%
    mutate(d = abs(c - median_diffs))

mad <- weightedMedian(genomecov$d, genomecov$count, ties="weighted")
```

**Output:** `MAD.${basename}.txt` — two columns: sample name, MAD value.

**Interpretation:** MAD measures the typical magnitude of depth change between consecutive windows, normalised by the mean depth. In uniform coverage, consecutive windows should have similar depths → low MAD. In single-cell data with amplification artefacts, large jumps between windows → high MAD.

---

## 6. Aggregating results — CreateInputApp

**File:** `CreateInputApp`

Once multiple samples have been processed, this script merges all individual `${sample}.SingleCheck.txt` files into one table for the Shiny app.

```bash
CreateInputApp <Samples.txt>
```

**`Samples.txt` format** — one sample path per line (without extension):

```
path/to/R1.T15
path/to/R20.S5
path/to/R9.S1
```

**Script logic:**

```bash
SamplesFile=$1
path=$(dirname "$SamplesFile")

# Write header
printf "Sample\tSequenced bases\tAnalysis depth\tBin size\tDelta\t\
% of unmapped reads\t% of reads mapped to the mitochondria\t\
Breadth\tAutocorrelation\tCoefficient of variation\t\
Gini coefficient\tMAD\tPotential contaminants\n" \
    > ${path}/SingleCheck.txt

# Append each sample's results
while read sample; do
    name=$(basename $sample)
    if [ ! -f "${sample}.SingleCheck.txt" ]; then
        rm -f ${path}/SingleCheck.txt
        echo "${sample}.SingleCheck.txt does not exist"
        exit 1
    fi
    cat ${sample}.SingleCheck.txt >> ${path}/SingleCheck.txt
done < $SamplesFile
```

**Output:** `${path}/SingleCheck.txt` — a headed, tab-delimited table ready to load into the Shiny app.

---

## 7. Output format and interpretation

### Per-sample output: `${NAME}.SingleCheck.txt`

Tab-delimited, one line, no header:

| Column | Name | Description |
|--------|------|-------------|
| 1 | Sample | Sample name (filename without extension) |
| 2 | Sequenced bases | Total bases in the original BAM (`raw_reads × mean_readlength`) |
| 3 | Analysis depth | Coverage used for metric calculation (downsampled or original) |
| 4 | Bin size | Window size used (`-w`, default 10 Mbp) |
| 5 | Delta | Lag used for autocorrelation (`-i`, default 1000 bp) |
| 6 | % unmapped reads | Fraction of primary reads that are unmapped in the downsampled BAM |
| 7 | % reads mapped to mitochondria | Mitochondrial read fraction — high values may indicate poor nuclear yield |
| 8 | Breadth | % of the genome covered by at least one read |
| 9 | Autocorrelation | Normalised autocovariance at lag Delta |
| 10 | Coefficient of variation | SD / mean of per-window depths |
| 11 | Gini coefficient | Lorenz-curve area metric of depth inequality |
| 12 | MAD | Weighted median of normalised absolute differences between consecutive windows |
| 13 | Potential contaminants | Comma-separated Metaphyler genus classifications |

### Metric interpretation

| Metric | Single-cell WGA | Bulk DNA |
|--------|----------------|----------|
| Gini coefficient | High (e.g. > 0.3) | Low (e.g. < 0.1) |
| CV | High | Low |
| Autocorrelation | High positive | Near zero |
| MAD | High | Low |
| Breadth | Variable | High |

All four metrics should agree directionally. Discordant values can indicate library-specific artefacts.

---

## 8. Helper and batch scripts

### `analysis/SingleCheckArray`

A SLURM array job that processes a list of samples. Uses `$SLURM_ARRAY_TASK_ID` to select a sample from the list and runs the same metric calculations as the main script. Useful for large cohorts.

### `analysis/SortSam.sh`

Wrapper around Picard `SortSam`:

```bash
java -jar picard.jar SortSam \
    INPUT=<in.bam> \
    OUTPUT=<out.bam> \
    SORT_ORDER=coordinate \
    VALIDATION_STRINGENCY=LENIENT \
    CREATE_INDEX=true
```

### `analysis/DownsampleSam.sh`

Standalone downsampling wrapper. Accepts a target depth and read length, computes the probability, and calls Picard `DownsampleSam` (same logic as Phase 4 of the main script).

### `analysis/GetSeqDepth.sh` / `GetSeqDepthFixingReadLength.sh`

Calculate sequencing depth using Picard `CollectAlignmentSummaryMetrics` and samtools. `GetSeqDepthFixingReadLength.sh` accepts a fixed read length instead of inferring it.

### `analysis/CheckSeqDepth.sh`

Verify the actual depth of a downsampled BAM file by running `GetSeqDepth` on it.

### `analysis/CollectInsertSizeMetrics.sh`

Wrapper for Picard `CollectInsertSizeMetrics`, generating insert-size distribution statistics and a histogram PDF.

### `analysis/SimulateSingleCellReads.sh`

Generates synthetic single-cell reads with coverage bias using ART (Artificial Read Tool). Parameters:

| Variable | Meaning |
|----------|---------|
| `AMPLICON_LENGTH` | Length of each simulated amplified region |
| `INTERAMPLICON_LENGTH` | Gap between amplicons (not amplified) |
| `AMPLICON_DEPTH` | Coverage within amplicons (default 20×) |
| `INTERAMPLICON_DEPTH` | Coverage in gaps (default 0×) |

This is used to generate ground-truth datasets for evaluating whether SingleCheck correctly identifies non-uniform libraries.

### `analysis/RunWANG.sh` / `RunHUANG.sh`

Batch submission scripts that loop over multiple window sizes (1 bp to 10 Mbp) and submit `SingleCheckArray` jobs for the Wang and Huang benchmark datasets respectively. Job IDs are saved to `.out` files for monitoring.

### `analysis/CreatePlotsWANG.sh` / `CreatePlotsHUANG.sh`

Aggregate metric outputs across window sizes into combined tables for cross-window-size visualisation.

---

## 9. Shiny app

Located in `ShinyApp/`. Launch from RStudio:

```r
shiny::runApp("ShinyApp/")
```

### Tabs

| Tab | Content |
|-----|---------|
| Documentation | Link to GitHub and rendered README |
| Input | Upload `SingleCheck.txt` (created by `CreateInputApp`) |
| Results | Summary statistics table + interactive ggplot2/plotly charts |

### Features

- Filter samples interactively in the results table.
- Select any metric column for the y-axis.
- Toggle log-scale y-axis.
- Download filtered results as a file.

---

## 10. End-to-end usage examples

### Example 1 — Single sample, pre-aligned BAM

```bash
SingleCheck -w 10000000 -i 1000 -t 3 -q 20 sample.bam
```

Produces: `sample.SingleCheck.txt`

### Example 2 — Single sample, paired-end FASTQ

```bash
SingleCheck \
    -r /path/to/hs37d5.fa \
    -w 10000000 -i 1000 -t 3 \
    sample_1.fastq.gz sample_2.fastq.gz
```

Produces: `sample.bam`, `sample.bam.bai`, `sample.SingleCheck.txt`

### Example 3 — Include chromosome X (human)

```bash
SingleCheck -X -w 10000000 sample.bam
```

### Example 4 — Custom chromosome set (non-human organism)

```bash
SingleCheck -c 'chr[1-9]|chrZ' -m chrMT sample.bam
```

### Example 5 — Skip downsampling (use original depth)

```bash
SingleCheck -N -w 5000000 sample.bam
```

### Example 6 — Batch processing and aggregation

```bash
# 1. Create sample list (paths without extension)
cat > Samples.txt <<EOF
/data/samples/R1.T15
/data/samples/R20.S5
/data/samples/R9.S1
EOF

# 2. Submit each sample as a SLURM job
while read sample; do
    sbatch SingleCheck -w 10000000 -i 1000 ${sample}.bam
done < Samples.txt

# 3. After all jobs finish, aggregate results
CreateInputApp Samples.txt

# 4. Launch Shiny app and load SingleCheck.txt from the Input tab
Rscript -e "shiny::runApp('ShinyApp/')"
```

### Example 7 — High-accuracy downsampling with a different target depth

```bash
SingleCheck \
    -d 0.5 \
    -s HighAccuracy \
    -w 10000000 -i 1000 \
    sample.bam
```

---

## Appendix — File naming conventions

| Pattern | Meaning |
|---------|---------|
| `${NAME}.bam` | Input or aligned BAM |
| `${NAME}.bam.bai` | BAM index |
| `${NAME}.${depth}X.bam` | Downsampled BAM (e.g. `sample.0.1X.bam`) |
| `${NAME}.${WSIZE}.regions.bed.gz` | mosdepth window depths (e.g. `sample.10000000.regions.bed.gz`) |
| `${NAME}.${WSIZE}.per-base.bed.gz` | mosdepth per-base depths |
| `${NAME}.${WSIZE}.${DELTA}.bed.gz` | Shifted per-base coverage |
| `${NAME}.${DELTA}.shiftedcov.txt` | Run-length-encoded depth pairs for autocorrelation |
| `${NAME}.${WSIZE}.freqs.txt` | Depth frequency table (input to Gini and CV) |
| `${NAME}.${WSIZE}.contiguous.txt` | Consecutive window depth pairs (input to MAD) |
| `${NAME}.unmapped.fasta` | Unmapped reads for Metaphyler |
| `Gini.${SAMPLE}.${WSIZE}.txt` | Gini coefficient result |
| `CV.${SAMPLE}.${WSIZE}.txt` | CV result |
| `Autocorrelation.${SAMPLE}.${DELTA}.txt` | Autocorrelation result |
| `MAD.${SAMPLE}.${WSIZE}.txt` | MAD result |
| `${NAME}.SingleCheck.txt` | **Final output** |
| `SingleCheck.txt` | Aggregated multi-sample table (CreateInputApp output) |
