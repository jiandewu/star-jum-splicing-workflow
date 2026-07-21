# STAR–JUM Splicing Workflow

A configurable two-pass STAR alignment and JUM workflow for comparing alternative splicing between two RNA-seq conditions.

This public repository contains workflow code and synthetic tests only. It does not contain sequencing reads, sample metadata, study identifiers, reference indexes, or JUM itself.

## What it does

1. Runs first-pass STAR alignment for every paired-end sample.
2. Combines the discovered splice junction files for a second STAR pass.
3. Sorts and indexes the second-pass BAM files with samtools.
4. Runs JUM stages A, B, and C for a two-condition comparison.

The script defaults to a dry run: it validates inputs and prints a shell-escaped execution plan without creating output. Add `--execute` only after reviewing that plan.

## Requirements

- Bash 4 or newer
- STAR
- samtools
- R and the packages required by your JUM release
- A separately obtained JUM installation containing `JUM_A.sh`, `JUM_B.sh`, `JUM_C.sh`, and `R_script_JUM.R`
- A compatible STAR genome index and `refFlat.txt`

The Conda file installs the core command-line dependencies, but JUM and its release-specific R packages must be installed separately according to its documentation.

```bash
conda env create -f environment.yml
conda activate star-jum-splicing
```

## Sample sheet

Create a tab-separated file with exactly four columns and exactly two condition labels:

```text
sample	condition	read1	read2
control_a	control	/secure/reads/control_a_R1.fastq.gz	/secure/reads/control_a_R2.fastq.gz
treated_a	treated	/secure/reads/treated_a_R1.fastq.gz	/secure/reads/treated_a_R2.fastq.gz
```

Sample and condition names may contain letters, digits, dots, underscores, and hyphens. Sample names must be unique. Keep real sample sheets outside Git because they may disclose identifiers and storage paths.

## Usage

Review a dry-run plan:

```bash
./run_workflow.sh \
  --samples /secure/config/samples.tsv \
  --genome-dir /references/star-index \
  --ref-flat /references/refFlat.txt \
  --jum-dir /opt/JUM \
  --design /secure/config/experiment_design.txt \
  --output-dir /secure/results/splicing \
  --read-files-command zcat \
  --threads 12
```

Execute the same plan:

```bash
./run_workflow.sh \
  --samples /secure/config/samples.tsv \
  --genome-dir /references/star-index \
  --ref-flat /references/refFlat.txt \
  --jum-dir /opt/JUM \
  --design /secure/config/experiment_design.txt \
  --output-dir /secure/results/splicing \
  --read-files-command zcat \
  --threads 12 \
  --execute
```

Run `./run_workflow.sh --help` for all options.

## Important assumptions

- Reads are paired-end and already decompressed by the filesystem or another mechanism understood by STAR. For gzip-compressed reads, supply `--read-files-command zcat` or the equivalent command for your environment.
- The sample sheet contains exactly two biological conditions.
- Condition names are sorted to assign JUM condition 1 and condition 2; inspect the dry-run plan before execution.
- JUM output conventions vary by release. This workflow follows the stage layout used by JUM 2.x and should be validated against the installed version.
- Reference files, genome assembly, read length, strandedness, replication thresholds, and experimental design remain study-specific responsibilities.

## Tests

The tests create temporary empty files and inspect the dry-run plan; no bioinformatics software or real data are used.

```bash
python -m unittest discover -s tests -v
bash -n run_workflow.sh
```
