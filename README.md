
#### SingleCheck is a program for assessing coverage dispersion of single cell DNA-seq libraries 

<!-- background / introduction -->

<!-- ## Installation -->

### Dependencies

* [samtools](http://www.htslib.org/) &#8805; 1.9
* [mosdepth](https://github.com/brentp/mosdepth) &#8805; 0.2.5
* [bedtools](https://bedtools.readthedocs.io/en/latest/) &#8805; 2.25
* [BWA-MEM](https://github.com/lh3/bwa) or [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) (only for FASTQ input)
* [MetaPhyler](http://metaphyler.cbcb.umd.edu/) plus `perl`
  (set `$METAPHYLER` to the full path of its `metaphyler.pl`)
* [BLAST](https://blast.ncbi.nlm.nih.gov/) — MetaPhyler shells out to it. MetaPhyler V1.13
  drives **legacy** BLAST (`blastall`, conda: `blast-legacy`); BLAST+ (`blastn`) works only
  if the marker database was re-formatted with `makeblastdb`. Because MetaPhyler looks BLAST
  up on `$PATH` by name, point `$BLAST_DIR` at the directory holding the binary rather than
  at the binary itself
* GNU coreutils (`sort` with `--parallel` / `--compress-program` / `--version-sort`) and `bc`
* optional, for speed: `pigz` or `bgzip`, and `zstd`
* R &#8805; 3.5 with the packages for the main program:
	* tidyr
	* dplyr
	* matrixStats
* R packages for the shiny app (optional):
  * shiny
  * shinydashboard
  * shinycssloaders
  * DT
  * data.table
  * ggplot2
  * plotly

`SingleCheck` checks all of this before it starts: every module, tool, R package, helper
script, reference index, input file and output/scratch directory is verified — and every
external tool is actually executed once, so a module that loaded but is broken is caught
too. If anything is missing the run aborts with a `[MISSING]` list and does no work.
On a cluster whose modules are named differently, point `$SINGLECHECK_ENV` at a file that
loads them (see `singlecheck_env.sh`) or export `MOSDEPTH`, `BEDTOOLS` and `METAPHYLER`
with absolute paths. Set `SINGLECHECK_STRICT_MODULES=1` to also abort when a `module load`
fails, even if the tool it provides was found elsewhere on `$PATH`.


## Usage

```
SingleCheck [options] <in.bam|in.fastq.gz> <in2.fastq.gz>

options:
  -h  show this help text
  -w <int> set the window size for gini coefficient, coefficient of variation and MAD efault: 10000000)
  -i <int> set increment/Delta for autocorrelation (default: 1000)
  -t <int> number of threads
  -f <int> flag of the reads to filter out
  -q <int>  mapping quality of the reads to analyze
  -r <ref.fa> reference genome
  -X include chromosome X in the analysis (only for human samples)
  -N do not perform downsampling. Extract statistics from original file
  -d downsampling sequencing depth. Ignore if -N (default: 0.1)
  -s <ConstantMemory|HighAccuracy|Chained> downsampling strategy method.
     Read Picard DownsampleSam documentation for details:
     https://broadinstitute.github.io/picard/command-line-overview.html#DownsampleSam
     Ignore if -N (default: ConstantMemory)
  -c <STRING> string containing all the names of the chromosomes you want to analyze parated by vertical bars enclosed in single or double quotation marks. Ignore if your ganism of study is human. Sintax of regular expression allowed
    	examples:
    		'2|3|4|X'
    		'[1-9]'
    		'chr[1-9|Z]'
    		'NC87126|NC78623'
    		'NW_[0-9]*'
  -m <STRING> mitochondrial contig name if different from MT or chrM. 
     Quotes are not required in this case
  -G use GPU alignment (NVIDIA Parabricks 'pbrun fq2bam') when available.
     FASTQ input only
```

`-t` defaults to `$SLURM_CPUS_PER_TASK`, so on a cluster just ask for the cores you want
(`./submit --cpus 32 SingleCheck sample.bam`). Intermediates are written to node-local scratch
(`$SLURM_TMPDIR`) and only the result file goes to the shared filesystem. See
[OPTIMIZATION_REPORT.md](OPTIMIZATION_REPORT.md) for the HPC behaviour and the environment
variables that switch each optimization off (`SINGLECHECK_SCRATCH`,
`SINGLECHECK_LEGACY_AUTOCORR`, `SINGLECHECK_NO_OVERLAP`).

Example
```bash
SingleCheck test/R1.T15.bam
```

##  Output and interpretation

The program generates a text file with the following columns:

* Name of the sample: single cell or unamplified control
* Sequenced bases: original number of bases present in the input BAM file
* Analysis depth
* Window size
* Delta
* % of unmapped reads	
* % of reads mapped to the mitochondria	
* Breadth: % of genome covered by &#8805;1 read	
* Autocorrelation	
* Coefficient of variation	
* Gini coefficient	
* MAD	
* Potential contaminants (Metaphyler genus file information condensed)


For comparing results from the different single cells simultaneously, we created a shiny app which can be run from Rstudio.

In order to create the input for the app you must run the following line:

```bash
CreateInputApp.sh <Samples.txt>
```

Example of Samples.txt
```
R1.T15
R20.S5
R9.S1
```

Then, you can load the output on the Input tab in the menu panel of the app. 

## Current issues

<!--
## FAQ

1. Can I use SingleCheck for performing quality control of single-cell RNA-seq data?

We have not tested the program for this aim so we do not provide support for this. 

2. Should I remove duplicates from my data?

We think is better to keep the duplicates for quality control the single cells.
-->

<!-- ## How to cite -->
