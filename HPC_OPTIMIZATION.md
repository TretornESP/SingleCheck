# SingleCheck — HPC Optimization Analysis for Very Large Inputs (300 GB – multi-TB)

Target environment: Linux SLURM cluster, Intel CPUs, NVIDIA GPUs.

---

## 0. The single most important framing point

SingleCheck is **not** a numerical program. It is a Bash script that *orchestrates* external
tools (`bwa`, `samtools`, `picard`, `mosdepth`, `bedtools`) and glues them together with
`awk` / `sort` / `grep` / `zcat` / `gzip` plus four tiny R scripts.

That changes which HPC techniques actually apply:

| Technique you asked about | Applicability here | Why |
|---|---|---|
| **OpenMP** | Indirect only | You don't own any hot loop. OpenMP already lives *inside* `bwa`, `samtools`, `mosdepth`, `sort --parallel`. You exploit it by passing more threads, not by writing `#pragma omp`. |
| **MPI** | Yes, but as *scatter–gather across chromosomes/nodes*, not as a linked MPI binary | The metrics are mathematically decomposable per chromosome (see §2). Multi-node parallelism = region splitting + SLURM, not `MPI_Send`. |
| **GPU (CUDA)** | Yes, in two specific places | (1) Alignment via NVIDIA Parabricks `fq2bam`. (2) Optionally the giant sort+aggregate via RAPIDS/cuDF. The statistics themselves are too small to bother. |
| **Cache / SIMD optimization** | Only if you rewrite the awk/sort glue as a C kernel | Until then there is no source-level loop to tune. Covered in §7. |
| **I/O optimization** | **The biggest lever** | At 300 GB–TB scale this pipeline is dominated by reading/writing/recompressing the BAM and by sorting a whole-genome per-base track. This is where most wall-clock is spent. |

**Bottom line:** for this pipeline the win ordering is
**Algorithmic redundancy elimination ≈ I/O ≫ tool replacement > multi-node scatter-gather > GPU alignment ≫ everything else.**
Do the cheap structural fixes *first* — several of them will beat a GPU.

---

## 1. Where the wall-clock actually goes on a multi-TB file

Walking the script ([SingleCheck](SingleCheck)) phase by phase and classifying the cost for a
*large original BAM* (the downsampled 0.1X BAM is tiny, so anything operating on it is cheap):

| Phase | Line(s) | Operation | Cost on a 300 GB–TB BAM | Bound by |
|---|---|---|---|---|
| 2 | [234-239](SingleCheck#L234-L239) | `bwa mem … \| samtools sort` (FASTQ path) | **Hours** | CPU (align) + I/O (sort spill) |
| 3c | [262-263](SingleCheck#L262-L263) | `samtools view -bF 2304 > primary.bam` then index | **Very expensive** — full decompress + **recompress + rewrite of the entire file** just to count reads | I/O + CPU (bgzip) |
| 3c | [268](SingleCheck#L268) | mean read length from first 1e6 reads | Cheap (SIGPIPE stops early) | — |
| 4 | [314-321](SingleCheck#L314-L321) | **Picard `DownsampleSam`** on the full file | **Very expensive** — single-threaded JVM streaming the whole BAM | CPU (1 core) + I/O |
| 5 | [332-339](SingleCheck#L332-L339) | `mosdepth` on the *downsampled* (tiny) BAM | Cheap | — |
| 6 | [346-362](SingleCheck#L346-L362) | shift per-base track + `bedtools unionbedg \| grep \| sort --version-sort \| awk` | **Expensive** — whole-genome per-base RLE, then a full external **sort** | CPU + I/O (sort temp) |
| 7/9 | [366-394](SingleCheck#L366-L394) | freq table + contiguous table (`zcat \| awk \| sort \| uniq`) | Moderate | CPU + I/O |
| 8/7/9 | R scripts | Gini/CV/Autocorr/MAD on small RLE tables | Cheap–moderate (single-threaded R, but small input) | CPU (1 core) |
| 10 | [401-402](SingleCheck#L401-L402) | `samtools view -f 0x4` over the **full original BAM** → FASTA → Metaphyler | **Expensive** — another full scan of the big file | I/O + CPU |

**The three killers are all on the original big file:** (3c) writing a full primary-only copy,
(4) Picard downsampling, and (10) the full unmapped scan. Plus (6) the whole-genome sort.
Notice the original BAM is **fully decompressed at least 3 separate times**. Fixing that
redundancy is free performance.

---

## 2. The structural insight that unlocks multi-node parallelism

Every metric is **decomposable per chromosome**, and the script already throws away
cross-chromosome information:

- **Gini, CV** ([GiniIndex.R](src/GiniIndex.R), [CoefficientOfVariation.R](src/CoefficientOfVariation.R)):
  computed from a **depth-frequency table** of windows. Frequency tables are *additive* — concatenate
  per-chromosome freq tables and sum counts.
- **MAD** ([MAD.R](src/MAD.R)): uses only *consecutive same-chromosome* window pairs; cross-chromosome
  boundaries are explicitly dropped at [SingleCheck:391](SingleCheck#L391).
- **Autocorrelation** ([Autocorrelation.R](src/Autocorrelation.R)): pairs base *p* with *p+Δ*; the
  `unionbedg` + RLE only ever pairs positions **within the same chromosome**.

Therefore the entire mosdepth → text-processing → metric chain can run **independently per chromosome
(or per region)** and be merged at the end with trivial additive reductions. This is the
"embarrassingly parallel" decomposition that maps cleanly onto:

- **GNU `parallel` / xargs -P** across cores on one node, or
- a **SLURM job array / `srun --multi-prog`** across many nodes (your "MPI" answer), or
- `mosdepth`'s own threads for the within-region part.

This is far more effective than trying to thread the awk/sort glue, because it parallelizes the
*whole* dependency chain, not one stage.

---

## 3. Tier 1 — Algorithmic / redundancy fixes (do these first; near-zero risk, large payoff)

These need no GPU and no new infrastructure. They attack the §1 killers directly.

### 3.1 Stop materializing `primary.bam` just to count reads — [262-291](SingleCheck#L262-L291)
You rewrite the entire BAM to disk only to run `idxstats`. Replace with a **header/index-only** count:

```bash
# raw primary read count without writing a copy:
raw_reads=$(samtools idxstats ${NAME}.bam | awk '{m+=$3;u+=$4} END{print m+u}')
# (idxstats reads the .bai — effectively instant; counts include 2ndary/supp,
#  so for primary-only use: samtools view -c -F 2304 -@ $THREADS ${NAME}.bam )
```

`samtools view -c -F 2304` streams once and never writes a multi-hundred-GB temp file.
**This alone removes one full decompress+recompress+write of the entire input.**

### 3.2 Replace Picard `DownsampleSam` with `samtools view -s` — [314-321](SingleCheck#L314-L321)
Picard `DownsampleSam` is a **single-threaded JVM** process — the worst possible thing to run on a
multi-TB file. `samtools view -s SEED.FRACTION` is multithreaded, streams, and keeps mate pairs
together by read-name hashing (sufficient for coverage QC):

```bash
# probability already computed; samtools wants SEED.FRACTION as one number
frac=$(printf '%.6f' "$probability")
samtools view -@ $THREADS -b -s 1${frac#0} ${NAME}.bam -o ${NAME}.${downsampling_depth}X.bam
samtools index -@ $THREADS ${NAME}.${downsampling_depth}X.bam
```
Expect this step to go from a single-core JVM crawl to near-I/O-bound. (Keep Picard available behind
a flag if exact template-level semantics are ever required.)

### 3.3 Fuse the passes over the original BAM
Currently the big file is scanned for: (a) read count, (b) downsampling, (c) unmapped extraction.
- (a) becomes index-only (§3.1).
- (c) unmapped extraction for Metaphyler ([401](SingleCheck#L401)) can be produced **as a by-product
  of the downsampling pass** with `samtools view`'s output filters, or at minimum run concurrently
  (it's independent of everything except the final `class` field).

### 3.4 Use multithreaded compressors everywhere
Every `gzip -c` ([350](SingleCheck#L350)) and the implicit bgzip in `samtools` are compression
hotspots. Replace `gzip` with **`pigz -p $THREADS`** or **`bgzip -@ $THREADS`**, and pass `-@` to all
`samtools view/sort/index`. Consider **zstd** (`--compress-program=zstd` for sort, `bgzip` w/ zstd) —
markedly faster at similar ratios.

### 3.5 Tame the whole-genome `sort` — [356](SingleCheck#L356)
`sort --parallel=$THREADS --version-sort -k4 -k5` over a per-base RLE track is the heaviest pure-CPU+I/O
step in the BAM-input path. Make it:
```bash
sort --parallel=$THREADS -S 80% \
     --compress-program=zstd \
     -T $SLURM_TMPDIR \           # node-local NVMe, NOT Lustre/NFS
     --version-sort -k4 -k5
```
- `-S` gives it real memory (default buffer is tiny → excessive spilling).
- `--compress-program` shrinks spill I/O.
- `-T` on **node-local scratch** is critical (see §6). The default `/tmp` may be small or networked.
- Where byte-order sorting is acceptable (the freq/contiguous tables, [369](SingleCheck#L369),
  [392](SingleCheck#L392)), prefix with `LC_ALL=C` for a large speedup.

> Even better: the whole §6 autocorrelation construction (shift + `unionbedg` + sort + RLE) can be
> replaced by a streaming, **sort-free** single pass that keeps a Δ-length ring buffer of depths per
> chromosome and emits `(depth, depth_{+Δ})` pairs directly. That eliminates the external sort
> entirely. This is the one place a small custom C/Rust kernel pays off (see §7).

---

## 4. Tier 2 — Drop-in faster tools

| Current | Replace with | Win | Notes |
|---|---|---|---|
| `bwa mem` [234](SingleCheck#L234) | **`bwa-mem2`** (already commented out at [215](SingleCheck#L215)!) | ~1.3–1.8× CPU | Bit-identical-enough output; just uncomment + reindex. |
| `picard DownsampleSam` | `samtools view -s` | single-core → multi-core | §3.2 |
| `gzip` | `pigz` / `bgzip -@` / `zstd` | scales with cores | §3.4 |
| BAM intermediates | **CRAM** (`--reference`) | 30–60% less I/O | Huge at TB scale; reference-based compression. |
| Picard `SortSam` ([analysis/SortSam.sh](analysis/SortSam.sh)) | `samtools sort -@ -m` | multi-core | If used in your batch path. |

---

## 5. Tier 3 — Parallelism (OpenMP / MPI / SLURM)

### 5.1 OpenMP / intra-tool threading
You don't write OpenMP; you *feed* it. Today `THREADS` defaults to 3 ([50](SingleCheck#L50)) and
`#SBATCH --cpus-per-task 3` ([6](SingleCheck#L6)). On a big node, request many more cores and pass
them through. Caveats:
- `mosdepth` decompression scales only to ~4 threads — more won't help *it*.
- `bwa-mem2` and `samtools sort` scale well to dozens of cores.
- `sort --parallel` scales to a point then becomes I/O-bound (→ §6).

So "more threads" is necessary but **sub-linear**; it does not solve a multi-TB single file alone.

### 5.2 MPI-style multi-node scatter–gather (the real multi-node answer)
Using the §2 decomposition, split the genome into regions and fan out:

```
                 ┌─ region chr1  → mosdepth → freq/contig/shift partials ─┐
 big BAM (CRAM) ─┼─ region chr2  → mosdepth → partials ───────────────────┼─ gather → sum freq tables
                 ├─ …                                                      │          concat contig/shift
                 └─ region chrN  → mosdepth → partials ───────────────────┘          → 4 R metrics once
```

- **One node:** `ls chroms | parallel -j$N 'process_region {}'`.
- **Many nodes:** SLURM **job array** (`--array=1-N`), one chromosome per task, then a small dependent
  gather job (`sbatch --dependency=afterok:<arrayjobid>`). You already have the array pattern in
  [analysis/SingleCheckArray](analysis/SingleCheckArray) — currently used for *per-sample* fan-out;
  reuse the same idea for *per-region* fan-out within one giant sample.
- True `MPI` (linked binary) buys you nothing over SLURM arrays here because the tasks are independent
  and communicate only a tiny additive reduction at the end. Use the array; skip MPI.

**Two levels of parallelism then coexist:** per-sample (existing, [analysis/SingleCheckArray](analysis/SingleCheckArray))
× per-region (new). For a cohort of huge files this is the dominant speedup.

### 5.3 Pipeline-stage overlap
Metaphyler ([401-402](SingleCheck#L401-L402)) is independent of the Gini/CV/MAD/autocorr branch.
Run it as a backgrounded `&` (or a separate array task) so the full-file unmapped scan overlaps the
mosdepth+stats branch instead of running serially.

---

## 6. Tier 1 — I/O optimization (co-equal #1 priority with §3)

At TB scale, I/O *is* the program. Concrete actions:

1. **Stage to node-local NVMe / burst buffer.** Read the input from Lustre once, copy (or
   `samtools view`-stream) to `$SLURM_TMPDIR` (or `/dev/shm` for the small intermediates), do all
   intermediate writes there, copy only the tiny final `.SingleCheck.txt` back. Set
   `TMPDIR`/`sort -T`/Picard `TMP_DIR` all to node-local scratch. The whole-genome `sort` spilling to
   networked `/tmp` will otherwise dominate wall-clock.
2. **Switch BAM → CRAM** for the input and any large intermediate. 30–60% smaller ⇒ proportionally less
   read bandwidth and decompression. Requires the reference (you already pass `-r`).
3. **Lustre striping** for the big input: `lfs setstripe -c 8 -S 4M <dir>` so a single huge file is read
   in parallel across OSTs. Match stripe count to your read concurrency.
4. **Eliminate gratuitous writes:** §3.1 (no `primary.bam` copy), CRAM intermediates, and prefer
   *streaming pipes* over temp files where a stage's output feeds exactly one consumer.
5. **Single decompression pass:** the original BAM is decompressed ≥3× today (§1). Fuse to ≤1–2.
6. **Right-size `samtools sort` memory** (`-m <per-thread>`), so the FASTQ-path sort spills to
   node-local scratch, never NFS.

---

## 7. Tier 4 — Cache / SIMD (only worthwhile if you write a kernel)

There is no source-level loop to cache-optimize today — it's all awk/sort/external tools. The payoff
appears only if you replace the §6 autocorrelation/MAD glue with a small compiled kernel:

- A streaming **C/C++/Rust** reader of `mosdepth`'s `per-base.bed.gz` that, per chromosome, keeps a
  **ring buffer of Δ depths** and emits `(depth, depth_{+Δ})` pairs — **sequential, cache-friendly,
  branch-light, SIMD-friendly** (auto-vectorizes on Intel AVX2/AVX-512 with `-O3 -march=native`).
  This removes the external sort *and* the `bedtools unionbedg` *and* the double-gzip — likely the
  largest single speedup in the BAM-input path.
- Compile with `-O3 -march=native -funroll-loops`; for the additive freq/RLE reductions, a tight loop
  over `int64` arrays vectorizes cleanly.
- Cache "optimization" at the orchestration level really means **reduce data volume** (CRAM, fewer
  passes, RLE) so the working set fits in node RAM / `/dev/shm` — already covered in §3/§6.

This is a medium-effort, high-reward item; scope it after Tier 1–2 land.

---

## 8. Tier 3/4 — GPU acceleration (Intel + NVIDIA)

Be selective — GPUs help in exactly two places:

### 8.1 Alignment (FASTQ path) — **the big GPU win** — [234-239](SingleCheck#L234-L239)
Replace `bwa mem | samtools sort` with **NVIDIA Parabricks `fq2bam`** (GPU BWA-MEM + sort + optional
mark-duplicates). On multiple A100/H100 GPUs a 30× WGS aligns in tens of minutes vs many CPU-hours;
output is a standard coordinate-sorted BAM/CRAM that drops straight into the rest of the pipeline.
Parabricks is now free for academic/research use. This only matters when input is FASTQ, not
pre-aligned BAM.

```bash
pbrun fq2bam --ref $REFERENCE --in-fq $FASTQ1 $FASTQ2 \
             --out-bam ${NAME}.bam --num-gpus $SLURM_GPUS
```
(Allocate GPUs via `#SBATCH --gres=gpu:a100:N`.)

### 8.2 The sort + aggregate — optional GPU path via RAPIDS/cuDF
The whole-genome sort+groupby in §6 can be done on GPU with **RAPIDS cuDF**: read mosdepth's bed into a
cuDF DataFrame, `groupby`/shift/aggregate on the GPU, emit the freq/contiguous/shifted tables. This can
replace the `awk | sort | uniq | bedtools` chain for the per-base track. Worth it only if (a) you keep
the per-base approach rather than the §7 streaming kernel, and (b) the table fits in GPU memory after
RLE. For most cases the §7 CPU streaming kernel is simpler and competitive — treat cuDF as the
"already have idle GPUs, want it off the CPU" option.

### 8.3 What is **not** worth a GPU
The four R metric computations operate on small RLE tables — kernel-launch overhead would dominate.
Leave them on CPU (or fold them into the §7 kernel). mosdepth has no GPU version and is already fast on
the *downsampled* BAM.

---

## 9. Recommended target architecture (per huge sample)

```
SLURM array task per chromosome/region  ─────────────────────────────┐
  stage region (CRAM) → node-local NVMe                               │
  ├─ branch A: samtools view -s  (downsample, multi-core)             │
  │            → mosdepth (--by W, per-base)                          │
  │            → streaming Δ-kernel (§7)  → shift/contig/freq partials│
  └─ branch B: samtools view -f4 → Metaphyler        (overlapped, &)  │
                                                                      │
gather job (afterok) ── sum freq tables, concat partials             │
                     ── run 4 R metrics once  → ${NAME}.SingleCheck.txt
FASTQ path only: Parabricks fq2bam on GPU node → CRAM, then as above ─┘
```

---

## 10. Prioritized roadmap (impact × effort)

| # | Action | Effort | Impact at TB scale | Tier |
|---|---|---|---|---|
| 1 | Drop `primary.bam` copy → `view -c`/`idxstats` | trivial | High (kills 1 full rewrite) | §3.1 |
| 2 | Picard → `samtools view -s` downsample | low | High (1 core → many) | §3.2 |
| 3 | Node-local scratch for all temp + `sort -T -S --compress-program` | low | High (sort/IO) | §5.1/§6 |
| 4 | `pigz`/`bgzip -@`, `samtools -@`, `LC_ALL=C` | trivial | Medium–High | §3.4/§3.5 |
| 5 | CRAM input + intermediates | low–med | High (I/O) | §6 |
| 6 | `bwa` → `bwa-mem2` (uncomment [215](SingleCheck#L215)) | trivial | Medium (FASTQ) | §4 |
| 7 | Per-chromosome SLURM-array scatter–gather | medium | High (multi-node) | §5.2 |
| 8 | Overlap Metaphyler branch | trivial | Medium | §5.3 |
| 9 | Parabricks `fq2bam` GPU alignment | medium | High (FASTQ only) | §8.1 |
| 10 | Streaming C/Rust Δ-kernel (replace sort+unionbedg) | medium–high | High (BAM path) | §7 |
| 11 | RAPIDS/cuDF sort+aggregate | high | Situational | §8.2 |

**Do 1–6 first.** They are hours of work, carry almost no correctness risk, and on a 300 GB+ file will
likely cut wall-clock more than the GPU items — because today's bottlenecks are redundant full-file
rewrites, a single-threaded JVM downsampler, and an unbounded external sort spilling to shared storage,
not a lack of FLOPs.

---

## 11. Correctness/robustness notes that bite specifically at scale

- **`awk` arithmetic is double-precision.** `raw_bases` for multi-TB ≈ 1e12–1e13, safely < 2^53, so OK —
  but anything multiplying counts by depths in the R scripts should stay in `as.numeric`/double (it does).
- **`sort` temp exhaustion:** without `-T`/`--compress-program`, the §6 sort can fill `/tmp` and fail
  silently mid-pipeline. Set them explicitly.
- **`samtools sort` RAM** on the FASTQ path: set `-m` so total ≤ allocated `--mem`, else OOM-kill.
- **Validate the `samtools view -s` semantics** against your downstream interpretation: it samples by
  read-name hash (template-consistent) rather than Picard's exact template accounting. For *coverage
  dispersion* QC this is equivalent; document the change.
- **Reproducibility:** `samtools view -s` seed is the integer part of `SEED.FRACTION`; keep it fixed (mirror Picard's `RANDOM_SEED=1`).
```
