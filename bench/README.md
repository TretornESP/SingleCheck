# bench/ — performance analysis of SingleCheck

Two different questions, two different scripts:

| Question | Script |
|---|---|
| *Which **chunk of the pipeline** is slow?* (samtools vs mosdepth vs bedtools vs R …) | `bench/bench_stages.sh` (here) |
| *Is **version B** of the whole pipeline faster than version A, and does it give the same numbers?* | `../benchmark_singlecheck.sh` |

## bench_stages.sh — per-chunk profiling

`SingleCheck` is an orchestrator: it chains external tools together. This script splits that
chain into 19–20 chunks, runs **each chunk N times**, times every repetition and writes a
report ranking the chunks by wall-clock, so you can see which tool actually costs the time.

The chunks are copies of the exact command lines in [../SingleCheck](../SingleCheck), and they
run **in pipeline order inside one private working directory**: the artifacts produced by
chunk *i* are the input of chunk *i+1*, exactly as in a real run. So the timings are real
pipeline timings, not synthetic micro-benchmarks.

### Run it

```bash
# minimum: an input file and how many repetitions per chunk
bench/bench_stages.sh -i sample.bam -r 5

# on the cluster, inside an allocation
salloc -c 8 --mem 80G
export SINGLECHECK_ENV=$PWD/singlecheck_env.sh
bench/bench_stages.sh -i /data/HG002.bam -r 3 --threads 8 --warmup

# or via the submit wrapper (unknown scripts default to 1 cpu / 8G -- always pass these)
./submit --cpus 8 --mem 80G --time 12:00:00 bench/bench_stages.sh -i /data/HG002.bam -r 3 --threads 8

# also time the whole pipeline end-to-end, to compare against the sum of the chunks
bench/bench_stages.sh -i sample.bam -r 3 --full

# only re-measure a couple of chunks (their inputs must already exist in work/)
bench/bench_stages.sh -i sample.bam -r 10 -s unionbedg,freq_table -o bench/stages_run1

bench/bench_stages.sh --list      # show the chunk keys
bench/bench_stages.sh --help      # all options
```

Wrapper options: `-i/--input`, `-I/--input2`, `-r/--reps` (default 3), `-o/--outdir`,
`-s/--stages`, `-l/--list`, `--warmup`, `--full`, `--strict`, `--env`, `--keep` (default),
`--clean`.
Pipeline options mirror `SingleCheck`'s flags with long names: `--threads`, `--window`,
`--delta`, `--depth`, `--ref`, `--flag`, `--mapq`, `--chroms`, `--mt`, `--include-x`,
`--no-downsample`.

### The chunks

| group | chunk | what it runs | reads the *original* file |
|---|---|---|---|
| Alignment | `align` | `bwa-mem2/bwa mem \| samtools sort`, `samtools index` (FASTQ input only) | yes |
| Raw depth | `genome_length` | `samtools view -H \| grep @SQ \| awk` | header only |
| Raw depth | `count_reads` | `samtools view -c -F 2304` | **yes** |
| Raw depth | `mean_readlen` | `samtools view \| head -1e6 \| cut \| sort \| uniq -c` | first 1e6 reads |
| Raw depth | `seq_depth` | `awk` + `bc` arithmetic | no |
| Downsampling | `downsample` | `samtools view -b -s SEED.FRACTION` | **yes** |
| Downsampling | `index_ds` | `samtools index` | no |
| Coverage | `mosdepth` | `mosdepth --fast-mode --by WSIZE` | no |
| Autocorrelation | `shift_track` | `zcat \| awk` shift `\| pigz/bgzip/gzip` | no |
| Autocorrelation | `unionbedg_sort` | **legacy**: `unionbedg \| grep \| sort --version-sort \| awk` | no |
| Autocorrelation | `unionbedg_hash` | `unionbedg \| grep \| awk` hash accumulator, no external sort | no |
| Gini/CV | `freq_table` | `zcat \| grep \| awk \| sort \| uniq -c` | no |
| Gini/CV | `gini_R` | `Rscript src/GiniIndex.R` | no |
| Gini/CV | `cv_R` | `Rscript src/CoefficientOfVariation.R` | no |
| Autocorrelation | `autocorr_R` | `Rscript src/Autocorrelation.R` | no |
| MAD | `contiguous` | `zcat \| paste \| grep \| awk \| sort \| uniq -c` | no |
| MAD | `mad_R` | `Rscript src/MAD.R` | no |
| Contamination | `unmapped_fasta` | `samtools view -f 0x4 \| awk` | **yes** |
| Contamination | `metaphyler` | `metaphyler.pl` | no |
| Final stats | `primary_bam_legacy` | **legacy**: `samtools view -bF 2304` copy + `samtools index` | no |
| Final stats | `idxstats_legacy` | **legacy**: `samtools idxstats` ×3 + `awk` | no |
| Final stats | `mapstats` | one streaming pass (MT %, unmapped %) | no |
| Final stats | `final_stats` | assemble the result line | no |

The chunks marked **yes** are the ones whose cost grows with the input size; everything from
`mosdepth` onwards works on the small downsampled BAM. That is the split
[../HPC_OPTIMIZATION.md](../HPC_OPTIMIZATION.md) §1 predicts — this script is how you confirm
it with numbers on your own data and cluster.

### Output

Inside the results directory (default `bench/stages_<host>_<timestamp>/`):

| file | content |
|---|---|
| `report.md` | the benchmark report: run metadata, per-chunk table (mean / min / median / max / sd / share / peak CPU % / peak RSS), chunks ranked by wall-clock with bars, totals per group, and skipped/failed chunks |
| `timings.tsv` | one row per chunk *and repetition* — feed it to R/pandas for your own plots |
| `logs/` | stdout+stderr and `/usr/bin/time -v` output of every single repetition |
| `work/` | the private working directory with every artifact (deleted with `--clean`) |
| `state.env` | the variables handed from chunk to chunk (genome length, depth, downsampling fraction …) |

### Dependency check

Before timing anything the script runs **the same dependency check as
[../SingleCheck](../SingleCheck)** — same module loading (`try_module`, `$BLAST_DIR`), same
tools and minimum versions, same R-package, helper-script, MetaPhyler/BLAST, reference-index,
input and scratch checks, same `[ok]`/`[WARN]`/`[MISSING]` output. So the benchmark measures
the pipeline in the environment the pipeline would actually get, and `$ALIGNER` is picked
here and used by the `align` chunk exactly as the pipeline does it.

One deliberate difference: a `[MISSING]` does **not** abort. The chunks that need the missing
tool are reported as `SKIPPED` and everything else is still measured — a partial environment
should still produce a useful report. Pass `--strict` to abort like `SingleCheck` does.

> The check is a copy of the pipeline's, not a shared library: **if you change one, change the
> other**. Say the word and I'll factor it into a single `src/check_dependencies.sh` that both
> source.

### Optimized vs legacy

Two steps were rewritten for [../HPC_OPTIMIZATION.md](../HPC_OPTIMIZATION.md) (see
[../OPTIMIZATION_REPORT.md](../OPTIMIZATION_REPORT.md)). The benchmark runs **both** forms by
default, so you get the speedup *and* the proof they agree, measured on your own data:

| optimized chunk | legacy chunk(s) it replaces |
|---|---|
| `unionbedg_hash` (no external sort) | `unionbedg_sort` |
| `mapstats` (one streaming pass) | `primary_bam_legacy` + `idxstats_legacy` |

`report.md` gets an *Optimized vs legacy: do they agree?* table — the shifted-coverage tables are
diffed row by row (after sorting, since the hash version is unordered) and the mapping
percentages are compared to a relative tolerance of 1e-9. `--no-compare` skips the legacy
variants once you trust them.

### How it works / what to watch out for

* Each chunk runs in a child process (`bench_stages.sh __stage <key>`) that sources
  `state.env`, so it sees exactly the variables the pipeline would have at that point.
  Values later chunks need (`genome_length`, `probability`, …) are handed back through
  `vals.env`.
* Before **every** repetition a chunk deletes its own outputs, so no repetition gets a free
  ride on a file a previous one already wrote.
* The input file is **symlinked**, never copied — safe for 300 GB–TB inputs.
* **min** is the fairest number for comparing chunks: the page cache warms up across
  repetitions, so rep 1 is usually the slowest. Use `--warmup` to discard that first cold
  repetition. If `sd` is large compared with `mean`, raise `-r`.
* A chunk whose peak CPU % sits near 100 is single-threaded — that is where `-t` buys you
  nothing and a different tool/approach is needed.
* Missing tools or missing inputs make a chunk `SKIPPED` (with the reason in the report)
  instead of aborting the run, so a partial environment still produces a useful report.
* Nothing outside the results directory is ever written or deleted.
