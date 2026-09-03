# BaseJumper BASEJ-WGS-QC

The BioSkryb BASEJ-WGS-QC pipeline evaluates the quality of single-cell
whole-genome (WGS) and whole-exome (WES) sequencing libraries. It aligns reads,
removes duplicates, and collects a suite of alignment, coverage, insert-size and
duplication QC metrics to assess whether a library is suitable for downstream
analysis.

The pipeline supports two modes via `--mode`:
- `wgs` (default) — whole-genome coverage metrics
- `exome` — target/hybrid-selection metrics for a chosen exome panel

This README documents the **fully open-source** run of the pipeline, which is the
**default** for the public release (`pipeline_tool = opensource`). It uses only
open-source tools and public/locally-built container images — no proprietary
Sentieon license is required. (To use the proprietary Sentieon path instead, pass
`--pipeline_tool sentieon`.)

# Pipeline Overview

Steps and tools used by the open-source run:

- (Optional) subsample reads using **SeqKit** (FASTQ) or **samtools** (Ultima CRAM); subsampling is off by default for WGS
- Merge multi-lane FASTQs for the same biosample (`cat`)
- Map reads to the reference genome using **BWA-MEM2**, sort/index with **samtools**
- Remove duplicate reads using **samtools markdup** (`-r`, matching Sentieon's `Dedup --rmdup`, so duplicates are dropped from the BAM rather than only flagged)
- Collect alignment, GC bias, insert size, mean-quality and quality-yield metrics using **Picard** (`CollectMultipleMetrics`)
- WGS mode: whole-genome coverage via **Picard** `CollectWgsMetrics`
- Exome mode: hybrid-selection coverage via **Picard** `CollectHsMetrics` for the selected `--exome_panel`
- Parse metrics into per-biosample Parquet files (custom Python)
- Generate QC plots, consensus scores, and a QC-status table (custom R)
- Aggregate metrics across biosamples and tools into an overall report using **MultiQC**

The custom container images for the open-source run are **built locally** from the
Dockerfiles in [`container/`](container/README.md) (see that README for the full
image list and build instructions). Public biocontainers (picard 3.0.0 for exome
HsMetrics, multiqc) are pulled directly from `quay.io`.

# Running Locally

Instructions for running BASEJ-WGS-QC on a local Ubuntu server.

## Install Java

```
sudo apt-get install default-jdk
java -version
```

## Install AWS CLI

```
curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
unzip awscliv2.zip
sudo ./aws/install
```

## Install Nextflow

```
wget -qO- https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
```

## Install Docker

```
# Add Docker's official GPG key:
sudo apt-get update
sudo apt-get install ca-certificates curl
sudo install -m 0755 -d /etc/apt/keyrings
sudo curl -fsSL https://download.docker.com/linux/ubuntu/gpg -o /etc/apt/keyrings/docker.asc
sudo chmod a+r /etc/apt/keyrings/docker.asc

# Add the repository to Apt sources:
echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.asc] https://download.docker.com/linux/ubuntu \
  $(. /etc/os-release && echo "$VERSION_CODENAME") stable" | \
  sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
sudo apt-get update

sudo apt-get install docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin
```

> **No license required.** The open-source run does not use Sentieon, so there is
> no license setup step. Build the custom images locally from [`container/`](container/README.md);
> public biocontainers (picard, multiqc) are pulled from `quay.io`.

## Resources Required

WGS alignment is compute-heavy. The open-source `BWAMEM2_ALIGN_DEDUP_METRICS`
process is sized for a full WGS sample; for local testing use subsampled inputs
and constrain resources, e.g.:

```
--max_cpus 8 --max_memory 30.GB
```

# Reference Data

The pipeline needs a reference-genome bundle (FASTA, BWA-MEM2 index, intervals).
These are read from the location given by `--genomes_base`, which defaults to the
BioSkryb shared S3 path (`s3://bioskryb-shared-data`).

For a public/local run, download the reference bundle for your genome once and
point `--genomes_base` at your local copy. The directory must contain the expected
`genomes/...` layout (e.g. `<genomes_base>/genomes/Homo_sapiens/NCBI/GRCh38/...`).
Then pass it on every run:

```
--genomes_base /path/to/local/genomes
```

# Test Pipeline Execution

The example **input** datasets referenced below are publicly available under
`s3://bioskryb-public-data/pipeline_resources` and are fetched automatically
during the run. The **reference genome** is resolved from `--genomes_base` (see
[Reference Data](#reference-data) above).

**Command (WGS, Illumina FASTQ, open-source run)**

```
nextflow run main.nf \
  --input_csv $PWD/tests/data/inputs/input_wgs_subsampled.csv \
  --genomes_base /path/to/local/genomes \
  --architecture x86 \
  --genome GRCh38 \
  --mode wgs \
  --outputDir test \
  --max_cpus 8 --max_memory 30.GB
```

**Command (WES/exome, Illumina FASTQ, open-source run)**

```
nextflow run main.nf \
  --input_csv $PWD/tests/data/inputs/input_exome_subsampled.csv \
  --genomes_base /path/to/local/genomes \
  --architecture x86 \
  --genome GRCh38 \
  --mode exome \
  --exome_panel "xGen Exome Hyb Panel v2" \
  --outputDir test \
  --max_cpus 8 --max_memory 30.GB
```

**Command (Ultima CRAM, open-source run)**

```
nextflow run main.nf \
  --input_csv $PWD/tests/data/inputs/input_ultima.csv \
  --genomes_base /path/to/local/genomes \
  --architecture x86 \
  --genome GRCh38 \
  --mode wgs \
  --outputDir test \
  --max_cpus 8 --max_memory 30.GB
```

> The open-source run is the default, so `--pipeline_tool` is not required. The
> public release is **x86 (amd64) only**.

## Input Options

The input is passed via an `input.csv` file. The platform is auto-detected from
the columns present.

- **Illumina / Element (FASTQ)** — columns `biosampleName`, `read1`, `read2`
  (optional `sampleId`, `reads`, `readLength`):

```
biosampleName,read1,read2
WGS-test1,s3://.../WGS-test1_R1.fastq.gz,s3://.../WGS-test1_R2.fastq.gz
```

  Multi-lane inputs are supported by pipe-delimiting (`|`) multiple paths within
  `read1`/`read2`.

- **Ultima (CRAM)** — columns `biosampleName`, `cram` (the `.cram.crai` index is
  auto-discovered):

```
biosampleName,cram
HG001-UltimaWGS_Z0022,s3://.../HG001-UltimaWGS_Z0022.cram
```

## Command Options

```
    Usage:
        nextflow run main.nf [options]

    Script Options: see nextflow.config

        [required]
        --input_csv         FILE    Path to input csv file

        [pipeline tool]
        --pipeline_tool     STR     Tool selection: 'opensource' (BWA-MEM2 + samtools markdup -r + Picard)
                                    or 'sentieon' (proprietary).
                                    DEFAULT: opensource

        [optional]
        --mode              STR     Analysis mode: 'wgs' (whole genome) or 'exome' (target panel).
                                    DEFAULT: wgs

        --exome_panel       STR     Exome panel name (used in exome mode), e.g.
                                    "xGen Exome Hyb Panel v2", "TWIST", "Twist Exome 2.0".
                                    DEFAULT: "xGen Exome Hyb Panel v2"

        --genomes_base      STR     Base path to the reference-genome bundle. Set this to your
                                    local download directory for public/offline runs.
                                    DEFAULT: s3://bioskryb-shared-data

        --genome            STR     Reference genome. Options: GRCh38, GRCm39
                                    DEFAULT: GRCh38

        --outputDir         DIR     Path to run output directory
                                    DEFAULT: results

        --skip_subsampling  BOOL    Skip read subsampling (subsampling is effectively off by
                                    default; max_total_reads is very large).
                                    DEFAULT: false

        --max_total_reads   VAL     Subsample target (total reads across mates) when enabled.

        --analysis_mode     STR     Sample type: 'single_cell' or 'bulk'. Bulk mode is for deep
                                    bulk libraries (e.g. 200x WGS) and requires
                                    --pipeline_tool sentieon. See "Bulk mode" below.
                                    DEFAULT: single_cell

        --bulk_split_reads  VAL     Reads per split part in bulk mode, counted across both mates.
                                    DEFAULT: 1000000000

        --bulk_coverage_cap VAL     Depth histogram cap passed to Sentieon WgsMetricsAlgo in
                                    bulk mode (wgs mode only).
                                    DEFAULT: 1000

        --help              BOOL    Display help message
```

## Bulk mode (deep libraries)

`--analysis_mode bulk` handles libraries that are too deep to align in a single
task, such as 200x WGS (~4.1 billion reads at 2x150). It is **Sentieon-only** for
now (`--pipeline_tool sentieon`); the open-source equivalent is not implemented yet.

What changes compared with the default `single_cell` mode:

- **Nothing is subsampled.** Instead of capping at `max_total_reads`, each FASTQ
  unit is split into parts of `--bulk_split_reads` reads (**SeqKit** `split2`,
  mates kept in lockstep), so a 200x library becomes 4-5 parts.
- **Parts are aligned in parallel** (`SENTIEON_ALIGN_SHARD`), one task each. All
  parts of a biosample share an identical `@RG` (same `ID`, `SM` and `LB`) so the
  merged library is treated as one read group and one library.
- **Merge, duplicate marking and metrics run in a single task**
  (`SENTIEON_MERGE_DEDUP_METRICS`). `sentieon driver` takes one `-i` per part and
  merges them in place, so a several-hundred-GB BAM is never staged in and out
  between steps. LocusCollector + Dedup see every part at once, which is what
  makes cross-part duplicates detectable.
- **Duplicates are marked, not removed.** Dedup runs without `--rmdup`, so the
  published BAM keeps duplicate reads flagged `0x400`. Coverage and duplication
  metrics still exclude them.
- **Multi-lane input is not merged.** Each lane is aligned as its own unit and the
  lanes regroup at the duplicate-marking step, so there is no full-library `cat`
  pass. Cross-lane duplicates are still marked.
- **`WgsMetricsAlgo` gets a raised `--coverage_cap`** (`--bulk_coverage_cap`,
  default 1000). The Picard-derived default of 250 truncates the depth histogram
  of a 200x library, which would depress `mean_coverage` and inflate
  `pct_exc_capped`. Single-cell metrics are unchanged.

Everything downstream (Parquet metrics, QC scores/plots, MultiQC, published BAM)
is identical to a single-cell run, and the output is alignment + QC only — no CNV
and no variant calling.

Constraints: FASTQ input only (CRAM/Ultima rows are rejected), mutually exclusive
with `--coverage_target`.

```bash
nextflow run main.nf \
  --input_csv bulk_samples.csv \
  --analysis_mode bulk \
  --pipeline_tool sentieon \
  --genome GRCh38 --mode wgs \
  --outputDir results_bulk
```

## Tool versions (open-source run)

- `SeqKit: 2.13.0`
- `BWA-MEM2: 2.2.1`
- `samtools: 1.21` (alignment/markdup image) / `1.23.1` (subsampling image)
- `Picard: 3.1.1` (metrics) / `3.0.0` (exome HsMetrics, quay biocontainer)
- `MultiQC: 1.33`

## Outputs

Outputs are written to the directory specified by `--outputDir`, including
per-biosample QC Parquet files, merged metric tables, QC plots, consensus scores,
and a MultiQC report.

# Testing

## Test Data Access

Test data is stored on Wasabi-backed S3 at `s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/`.

To access the test data:

**Step 1 — Get your access keys**

Retrieve your AWS credentials from BioSkryb support (contact basejumper-support for the access link).

**Step 2 — Set environment variables**

```bash
export AWS_ACCESS_KEY_ID=<provided_access_key>
export AWS_SECRET_ACCESS_KEY=<provided_secret_key>
export AWS_DEFAULT_REGION=us-east-1
```

## Running a Test

Run the pipeline with the provided test input CSV (2 WGS samples, ~35 min):

```bash
nextflow run main.nf \
  --input_csv tests/data/inputs/nftest_input.csv \
  --max_cpus 8 --max_memory 24.GB --architecture x86 \
  --genome GRCh38 --mode wgs --skip_subsampling true \
  --pipeline_tool opensource \
  --outputDir results_test
```

## nf-test (Automated Testing)

Install nf-test (requires Java 11+):

```bash
curl -fsSL https://code.askimed.com/install/nf-test | bash
mv nf-test /usr/local/bin/
```

Run the automated tests:

```bash
# Run all tests
nf-test test

# Run only the GRCh38 opensource test
nf-test test tests/main.nf.test --tag GRCh38
```


# Need Help?

If you need any help, please [submit a helpdesk ticket](https://bioskryb.atlassian.net/servicedesk/customer/portal/3/group/14/create/156).

# References

NOTE: Several studies have utilized BaseJumper pipelines as part of standard
quality control processes implemented through ResolveServices<sup>SM</sup>. While
these pipelines may not be explicitly cited, they are integral to the methodologies described.
