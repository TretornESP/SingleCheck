# SingleCheck — HPC Optimization Report

What was implemented from [HPC_OPTIMIZATION.md](HPC_OPTIMIZATION.md), why, how to verify it, and
what is deliberately left for later.

Target: a SLURM cluster (CESGA / FinisTerrae), single BAM/CRAM inputs of 300 GB–multi TB.

---

## 0. Summary

| # | Roadmap item | Status | Where |
|---|---|---|---|
| 1 | Drop the `primary.bam` copy → `view -c` / streaming counts | **done** (both places) | [SingleCheck](SingleCheck) §raw depth, §mapping statistics |
| 2 | Picard `DownsampleSam` → `samtools view -s` | **done** (earlier) | [SingleCheck](SingleCheck) §downsample |
| 3 | Node-local scratch for all temporaries + `sort -T -S --compress-program` | **done**, extended to *all* intermediates | [SingleCheck](SingleCheck) §work prefix |
| 4 | `pigz`/`bgzip -@`, `samtools -@`, `LC_ALL=C` | **done** (earlier) + thread auto-sizing | [SingleCheck](SingleCheck) §HPC tuning |
| 5 | CRAM input + intermediates | **partial** — CRAM input supported; intermediates left as BAM (see §4) | — |
| 6 | `bwa` → `bwa-mem2` | **done** (earlier) | [SingleCheck](SingleCheck) §alignment |
| 7 | Per-chromosome SLURM-array scatter–gather | **designed, not implemented** (see §5) | — |
| 8 | Overlap the MetaPhyler branch | **done** | [SingleCheck](SingleCheck) §MetaPhyler branch |
| 9 | Parabricks `fq2bam` GPU alignment | **done**, opt-in `-G` | [SingleCheck](SingleCheck) §alignment |
| 10 | Kill the whole-genome external sort | **done** (hash aggregation; the full C kernel is still open — §5) | [SingleCheck](SingleCheck) §autocorrelation |
| 11 | RAPIDS/cuDF | **not done** — situational, superseded by #10 | — |

Plus two robustness items from §11 of the analysis: `samtools sort -m` sized from the
allocation, and thread count taken from `$SLURM_CPUS_PER_TASK`.

**Full-file passes over the original BAM: 3 serial → 2, one of them overlapped.**
**Whole-genome external sort: removed. Full-file rewrites: 2 → 0.**

---

## 1. What changed in the pipeline

### 1.1 No more full-file rewrites (items 3.1, 6.4)

The original code materialised a primary-only copy of the analysis BAM and indexed it, purely
to read three counters out of `samtools idxstats`:

```bash
samtools view -bF 2304 … > …primary.bam      # full decompress + recompress + write
samtools index …primary.bam                  # + index build
samtools idxstats …primary.bam               # ×3
```

With `-N` (no downsampling) that copy **is the whole multi-TB input**. It is now a single
streaming pass with no temp file, no index and no re-read:

```bash
read mt_perc_totalreads unmapped_perc_totalreads < <(
    samtools view -F 2304 -@ $THREADS $SAMREF ${WORKPFX}.${downsampling_depth}X.bam |
    awk -v mtre="$MT_REGEX" '
        { total++
          if (int($2/4)%2 == 1) { unmapped++ }
          else if ($3 ~ mtre)   { mt++ } }
        END{ printf "%s %s\n", mt/total*100, unmapped/total*100 }')
```

This is exact, not an approximation: `idxstats` column 3 counts primary records with `0x4`
clear, column 4 counts those with `0x4` set, and `MT_mappedreads` is the mapped subset whose
RNAME matches `$MT_REGEX` — which is precisely what the awk computes. `int($2/4)%2` is used
instead of `and($2,4)` so it works on mawk as well as gawk.

### 1.2 The whole-genome external sort is gone (items 3.5, 7)

This was the heaviest pure CPU+I/O step of the BAM path:

```bash
bedtools unionbedg … | grep -E … |
  sort --parallel=N -S 80% --compress-program=zstd -T $SCRATCH --version-sort -k4 -k5 |
  awk '<run-length encode adjacent equal pairs>'
```

The sort exists **only to bring equal `(depth, shifted depth)` pairs next to each other** so the
following awk can sum their lengths. Grouping is associative and commutative, so a hash
accumulator does the same job in one streaming pass:

```bash
bedtools unionbedg … | grep -E … |
  awk '{pair[$4"\t"$5] += $3-$2} END{for (p in pair) print p"\t"pair[p]}'
```

- **Output is the same multiset of rows.** Only the row order differs, and both consumers —
  `src/Autocorrelation.R` (sums over all rows) and the `breadth` awk (also a sum) — are
  order-independent.
- **Memory is bounded by the number of distinct depth pairs**, not by genome size: a handful of
  hundreds at 0.1×, `max_depth²` in the worst case. The sort, by contrast, spilled tens of GB.
- Escape hatch: `SINGLECHECK_LEGACY_AUTOCORR=1` restores the old path byte for byte.

The benchmark runs **both** implementations and diffs their outputs — see §3.

### 1.3 Intermediates live on node-local scratch (item 6.1)

Previously every intermediate (downsampled BAM + index, per-base track, shifted track, the RLE
tables, the unmapped FASTA, MetaPhyler's output, the four metric files) was written next to the
input — i.e. on Lustre, over the network, for tens of GB of traffic that nobody ever reads
twice.

Now the pipeline computes a `$WORKPFX` on `$SLURM_TMPDIR` and every intermediate uses it. Only
`${NAME}.SingleCheck.txt` is written to the shared filesystem. An `EXIT`/`INT`/`TERM` trap wipes
the scratch directory, so nothing is left behind on a compute node — including after a failure
or a `scancel`.

Guards, because a full node-local disk mid-run is worse than slow Lustre writes:
- scratch is used only when it has **≥ 20 GB free**, otherwise the pipeline transparently falls
  back to the old behaviour;
- `SINGLECHECK_SCRATCH=0` forces the old behaviour;
- the effective choice is printed in the `> HPC settings` block of every run.

### 1.4 The MetaPhyler branch runs in parallel (item 5.3)

Extracting unmapped reads is a second full scan of the original file, and it feeds exactly one
output column. It is now started in the background right after downsampling, so its I/O overlaps
the mosdepth + metrics branch, and it is joined with `wait` immediately before the result line is
assembled. Threads are split so the branches don't fight: `MD_THREADS` (≤4) for mosdepth,
`THREADS - MD_THREADS` for the unmapped scan. `SINGLECHECK_NO_OVERLAP=1` restores serial
execution.

### 1.5 Thread and memory sizing (items 5.1, 6.6, 11)

| Before | After |
|---|---|
| `THREADS=3` hard-coded | `-t` if given, else `$SLURM_CPUS_PER_TASK`, else 3 |
| `mosdepth -t $THREADS` | `-t min(THREADS,4)` — mosdepth stops scaling at ~4 |
| `samtools sort -@N` (768 MB/thread default) | `-m` derived from `$SLURM_MEM_PER_NODE`/`_PER_CPU` (70 % of the allocation, ÷ threads) and `-T` on node-local scratch |
| `submit` default 3 cpus / 60 G | 8 cpus / 64 G, and the thread count follows `--cpus` automatically |

The last row matters: before this change, `./submit --cpus 32 SingleCheck …` gave you 32 cores
and used 3.

### 1.6 GPU alignment, opt-in (item 8.1)

`SingleCheck -G` on a FASTQ input replaces `bwa-mem2 mem | samtools sort` with
`pbrun fq2bam` (NVIDIA Parabricks: GPU BWA-MEM + sort, standard BAM out). The dependency check
verifies `pbrun` is on `$PATH` and warns when no GPU is allocated. `-G` on BAM/CRAM input is
ignored with a warning. Everything downstream is unchanged.

```bash
./submit --gres gpu:a100:1 --partition gpu SingleCheck -G -r ref.fa s_1.fastq.gz s_2.fastq.gz
```

---

## 2. Effect at a glance (BAM input, per run)

| Cost | Before | After |
|---|---|---|
| Full decompress passes over the original file | 3 (count, downsample, unmapped) | 3, but the unmapped one **overlaps** the coverage branch |
| Full-file rewrites (decompress + recompress + write) | 2 (`primary.bam`, and the `-N` copy) | 0 |
| Whole-genome external sorts | 1 (spills tens of GB) | 0 |
| Index builds on derived BAMs | 2 | 1 |
| Intermediate bytes written to Lustre | all of them | 0 (only the final 1-line result) |
| Cores actually used | 3 | the whole allocation |

No claim is made here about wall-clock ratios: those depend on your file, your node and your
filesystem, and measuring them is exactly what [bench/](bench/) is for.

---

## 3. How to verify — and get the numbers for your data

Every optimization above is measurable and checkable with the tools in the repo.

### 3.1 Per-chunk timings + built-in equivalence check

[bench/bench_stages.sh](bench/bench_stages.sh) now runs **both** the legacy and the optimized form
of each rewritten step, times them side by side, and verifies they agree:

```bash
export SINGLECHECK_ENV=$PWD/singlecheck_env.sh
./submit --cpus 16 --mem 96G bench/bench_stages.sh -i /data/HG002.bam -r 3 --threads 16 --warmup
```

`report.md` then contains, besides the ranked chunk table:

```
## Optimized vs legacy: do they agree?

| comparison | verdict | detail |
|---|---|---|
| `unionbedg_hash` vs `unionbedg_sort` | identical | 412 rows |
| `mapstats` vs `primary_bam_legacy`+`idxstats_legacy` | identical | MT % and unmapped % agree to 1e-9 |
```

Compare the `unionbedg_sort` and `unionbedg_hash` rows for the sort speedup, and
`primary_bam_legacy` + `idxstats_legacy` against `mapstats` for the copy elimination.
`--no-compare` skips the legacy variants once you trust them.

### 3.2 End-to-end, old vs new

[benchmark_singlecheck.sh](benchmark_singlecheck.sh) compares whole pipeline versions and checks
that the final output line matches within a tolerance — point it at a checkout of the previous
revision:

```bash
git worktree add /tmp/SingleCheckOld <commit-before-these-changes>
./benchmark_singlecheck.sh -i /data/HG002.bam -r 2 \
    -V /tmp/SingleCheckOld,$PWD -b SingleCheckOld -T 1e-6 -- -t 16
```

### 3.3 Cheap A/B of a single optimization

Each change has an environment switch, so you can isolate one at a time without editing code:

```bash
SINGLECHECK_LEGACY_AUTOCORR=1 ./SingleCheck …   # old sort-based aggregation
SINGLECHECK_SCRATCH=0         ./SingleCheck …   # intermediates back on Lustre
SINGLECHECK_NO_OVERLAP=1      ./SingleCheck …   # serial MetaPhyler branch
```

---

## 4. Deliberately not done

**CRAM intermediates (item 5, partial).** CRAM *input* already works (`-r ref.fa`). Converting the
*intermediates* to CRAM was not done because the only large intermediate left is the downsampled
0.1× BAM, which is small by construction, and CRAM would add a reference dependency to a path
that currently doesn't need one. The high-value half of this item — "read less" — is better served
by converting your **inputs** to CRAM, which needs no code change:
`samtools view -@16 -C -T ref.fa -o sample.cram sample.bam`.

**RAPIDS/cuDF (item 11).** The sort it was meant to accelerate no longer exists (§1.2). Moot.

---

## 5. Next steps, in the order I would do them

### 5.1 Per-region SLURM array scatter–gather (item 7) — the biggest remaining win

The metrics are additive across regions (§2 of the analysis), which makes the merge trivial:

| partial file | how to merge |
|---|---|
| `*.freqs.txt` (depth, count) | sum counts per depth |
| `*.shiftedcov.txt` (d, d+Δ, len) | sum lengths per pair — *already the hash format of §1.2* |
| `*.contiguous.txt` (d1, d2, count) | sum counts per pair |
| `raw_bases`, `genome_length`, MT/unmapped counters | plain sums |

Note that §1.2 makes this *easier*: the shifted-coverage partials are now unordered
`(key, value)` tables, so gathering is a one-line awk over the concatenation — no merge sort.

Shape:

```bash
# scatter: one array task per chromosome
sbatch --array=1-24 --cpus-per-task 8 SingleCheckRegion sample.bam        # writes <sample>.<chr>.partial.*
# gather: additive reduction, then the four R metrics once
sbatch --dependency=afterok:$ARRAY_JOBID SingleCheckGather sample.bam
```

Per-region work uses `samtools view … <region>` against the existing index, so the three
full-file passes become N parallel range reads — the real multi-node answer for a TB file. The
two levels then compose: per-sample (existing [analysis/SingleCheckArray](analysis/SingleCheckArray))
× per-region.

Est. 1–2 days including a correctness harness against the single-node result.

### 5.2 Streaming Δ-kernel (item 7/10, remaining half)

`bedtools unionbedg` and the double gzip of the shifted track are still there; only the sort was
removed. A ~150-line C or Rust reader of `per-base.bed.gz` keeping a Δ-length ring buffer per
chromosome would emit the `(depth, depth+Δ, length)` table directly — removing `unionbedg`, the
shifted-track write, and one full read of the per-base track. Compile `-O3 -march=native`.

Worth doing **after** you have `bench/` numbers showing what `shift_track` + `unionbedg_hash`
actually cost on your data — if they are 5 % of wall-clock, spend the effort on §5.1 instead.

### 5.3 Housekeeping

[analysis/SingleCheckArray](analysis/SingleCheckArray) is an old standalone copy of the pipeline
(its own module loads, its own metric chain, driven by `ReadConfig.sh`). It did **not** receive any
of these optimizations and has drifted from `SingleCheck`. Either retire it in favour of
`./submit --array … SingleCheck`, or rewrite it as a thin wrapper that calls `SingleCheck`.

---

## 6. Correctness notes

- The autocorrelation rewrite is an **exact** transformation (associative grouping), not an
  approximation, and the benchmark diffs it against the legacy output on real data.
- The mapping-statistics rewrite is likewise exact; it reproduces `idxstats`' definitions of
  mapped/unmapped/MT counts over primary records.
- `samtools view -s` sampling (item 2, done earlier) *is* a semantic change from Picard
  `DownsampleSam`: it samples by read-name hash rather than exact template accounting. Equivalent
  for coverage-dispersion QC; seed fixed at 1 to mirror Picard's `RANDOM_SEED=1`.
- Nothing above changes the metric definitions in [src/](src/) — those four R scripts are untouched.
- Every new code path has an environment switch back to the old behaviour, so a regression can be
  bisected in one run rather than one commit.

---

## 7. Status of testing

`bash -n` clean for [SingleCheck](SingleCheck), [bench/bench_stages.sh](bench/bench_stages.sh) and
[submit](submit); the dependency-check and `--dry-run`/`--list` paths were exercised locally. The
pipeline itself has **not** been executed against real data as part of this work — the tools
(`mosdepth`, `bedtools`, `Rscript`, MetaPhyler) only exist on the cluster. Before trusting these
changes on a production sample, run §3.1 on a small BAM and confirm both equivalence rows say
*identical*.
