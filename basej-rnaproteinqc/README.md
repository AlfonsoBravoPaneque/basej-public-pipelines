# basej-rnaproteinqc

Single-cell **RNA + Protein (ADT) QC** pipeline for the **ResolveDOGMA** product.

It runs the standard RNA QC backbone (identical to `basej-rnaqc`) and, on the same
per-biosample reads, quantifies antibody-derived tags (ADT) and computes CPM-normalized
protein values.

## Pipeline stages

```
RNA QC (reused from basej-rnaqc):
  MERGE_MULTILANE_FASTQ → SEQKIT_SAMPLE → FASTP_TRIM → STAR_ALIGN →
  SAMTOOLS_INDEX_FILTER → HTSEQ_COUNTS / QUALIMAP_BAMRNA → (summaries) → RNA_QC_PLOTS

Protein / ADT (reused from subworkflows/protein_quant_wf):
  SEQTK_SEQ → CUSTOM_CAT_FASTA → BOWTIE2_BUILD_S → BOWTIE2 → CUSTOM_COUNT

Integration (new):
  PROTEIN_CPM  — merge ADT counts + RNA total reads + DNA-background subtraction → CPM
```

## CPM normalization

Per antibody, per sample:

```
CPM = max(antibody_count − matched_DNA_background, 0) ÷ (total_reads ÷ 1,000,000)
```

- **background** = the paired DNA-amplified sample's count for that antibody. Partners are
  matched by name: `_R_`→`_D_`, or `-R-`→`-D-` matched on well id (e.g. `01A`). No partner ⇒ background 0.
- **total_reads** = the sample's RNA total reads, taken from `rnaqc_all_metrics.tsv`.

This reproduces the manual Excel calculation (`MAX(AY−BS,0) / (DG/1e6)`); reference
implementation is `tmp/resolvedogma_assets/cpm_calculation.py`.

## Inputs

Same samplesheet as `basej-rnaqc` (`--input_csv`): `biosampleName,read1,read2[,groups]`.
Illumina FASTQ only for the protein side (CRAM/Ultima RNA path is preserved but ADT quant
is skipped for Ultima).

Key params:
- `barcode` — antibody barcode reference FASTA (default: TotalSeqA cocktail on shared data).
- `skip_protein_quant` — set true to run RNA QC only.

## Outputs

- `tables/rnaqc_summary/...` — per-biosample RNA QC Parquet (Athena).
- `tables/protein_cpm_summary/...` — per-biosample ADT CPM Parquet, long format (**new Iceberg table**, see below).
- `workflow_outputs/.../protein_cpm/protein_adt_cpm.tsv` + `protein_adt_raw_counts.tsv` — wide tables for the R heatmap scripts.
- `multiqc_report.html` — includes an ADT CPM summary table.

## Containers

All images are in ECR. The one custom image for this pipeline (the `PROTEIN_CPM` step) is
defined in-repo under `containers/protein_cpm/`; everything else is a mirrored biocontainer
or a custom image reused from the RNA QC backbone. See `containers/README.md` for the full
image list and build/push/mirror commands.

## Design assumptions (confirm / revisit)

1. **ADT reads = RNA reads.** Protein quant runs on the same FASTQs as RNA (matching the
   standalone protein-quant pipeline, which aligns the barcode FASTA against the sample's
   own reads). If ADT is a separate sub-library in production, the samplesheet/wiring needs
   dedicated ADT read columns.
2. **DNA background comes from `_D_`/`-D-` samples present in the same run.** For background
   subtraction to apply, the DNA-amplified partners must be rows in `--input_csv`; otherwise
   background is 0 (CPM = count / (reads/1e6)).
3. **`total_reads`** (RNA, R1+R2 combined) is used as the CPM denominator, standing in for the
   Excel "BJ calculated raw read pairs". Adjust if a different denominator is required.

## Validation status

- `nextflow config` parses; `nextflow run -stub-run` builds the full workflow graph
  (RNA + protein processes wired correctly).
- Full `nf-test` (`cd pipelines/basej-rnaproteinqc && nf-test test`) needs S3 input access +
  containers + compute; not run in the dev sandbox.
- **TODO before production:** validate `PROTEIN_CPM` output against a known-good previous run
  (a prior run's inputs + delivered `*_cpm.txt`), and register the `protein_cpm_summary`
  Iceberg table so the new Parquet is queryable in Athena.

## Downstream Iceberg note

`protein_cpm_summary` is a **new** metrics table. Its Parquet schema (columns: `biosample`,
`antibody`, `raw_count`, `background`, `cpm`, `total_reads`, `dna_partner`, plus run
metadata) must be registered/evolved in the platform Iceberg catalog before Athena queries
will see it. `rnaqc_summary` and `gene_count_summary` rows are written with
`pipeline = "basej-rnaproteinqc"` (additive string value, no schema change).
