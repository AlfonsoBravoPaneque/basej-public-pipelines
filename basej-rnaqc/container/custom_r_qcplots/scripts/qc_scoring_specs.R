# Per-assay tier specifications for the QC composition plots.
#
# Kept separate from qc_scoring.R so the thresholds are readable in one place and can be
# retuned without touching the scoring machinery. These are BioSkryb-internal values and
# depend on cell type and biological context - the appropriate cutoff for a tumour biopsy is
# not the appropriate cutoff for a cell line.
#
# Tier 4 mirrors each assay's established gates so `tier >= 4` matches the historical
# "all gates pass" verdict and published counts reconcile. Tier 5 tightens the CNV metrics,
# which are what actually discriminate quality once the coverage gates are met. Tier 3 relaxes
# both, giving a usable-if-capacity-allows band.
#
# Depth handling differs by assay:
#   DNA / methyl : absolute read floor, since the low-pass target is fixed and characterised
#   HiFi / RNA   : median-relative floor, because no absolute number is established. A fixed
#                  bottom percentile was considered and rejected - it condemns a constant
#                  fraction of cells even when the whole run is uniformly good.

source_qc_scoring <- function() {
  for (p in c("/usr/local/bin/qc_scoring.R", "qc_scoring.R",
              file.path(dirname(sys.frame(1)$ofile %||% "."), "qc_scoring.R"))) {
    if (file.exists(p)) { source(p); return(invisible(p)) }
  }
  stop("qc_scoring.R not found alongside the plot script or in /usr/local/bin")
}
`%||%` <- function(a, b) if (is.null(a)) b else a


# ---------------------------------------------------------------- DNA (basej-dnaqc, basej-wgs)
# Gates: PreSeq library complexity, PCT_CHIMERAS, chrM, CNV MAPD, CNV skew.
qc_spec_dna <- function(cutoff_preseq = 3.5e9, cutoff_chimeras = 0.20, cutoff_chrM = 0.20,
                        cutoff_cnv_mapd = 0.25, cutoff_cnv_sk = 0.25,
                        tight_cnv_mapd = 0.20, tight_cnv_sk = 0.20,
                        borderline_preseq = 2.5e9, borderline_cnv = 0.35) {
  list(
    "5" = list(preseq_count = c(gt = cutoff_preseq),
               PCT_CHIMERAS = c(lt = cutoff_chimeras),
               chrM         = c(lt = cutoff_chrM),
               MAPD_CNV_Log2 = c(lt = tight_cnv_mapd),
               SKEW_CNV      = c(lt = tight_cnv_sk)),
    "4" = list(preseq_count = c(gt = cutoff_preseq),
               PCT_CHIMERAS = c(lt = cutoff_chimeras),
               chrM         = c(lt = cutoff_chrM),
               MAPD_CNV_Log2 = c(lt = cutoff_cnv_mapd),
               SKEW_CNV      = c(lt = cutoff_cnv_sk)),
    "3" = list(preseq_count = c(gt = borderline_preseq),
               MAPD_CNV_Log2 = c(lt = borderline_cnv),
               SKEW_CNV      = c(lt = borderline_cnv))
  )
}

# ------------------------------------------------- methylation (basej-methylqc, -methylqc-rs)
# Same amplification-quality metrics as DNA, but PreSeq is not produced by the methyl route:
# depth comes from total_reads, and Methylated_CpG_pct is added as the assay-specific metric.
# A methylation percentage far outside the expected band means conversion or the library
# failed, independent of how evenly the genome amplified.
qc_spec_methyl <- function(cutoff_reads = 1e6, cutoff_chimeras = 0.20, cutoff_chrM = 0.20,
                           cutoff_cnv_mapd = 0.25, cutoff_cnv_sk = 0.25,
                           tight_cnv_mapd = 0.20, tight_cnv_sk = 0.20,
                           methyl_min = 40, methyl_max = 90,
                           tight_methyl_min = 50, tight_methyl_max = 85,
                           borderline_reads = 5e5, borderline_cnv = 0.35) {
  list(
    "5" = list(total_reads   = c(gt = cutoff_reads),
               PCT_CHIMERAS  = c(lt = cutoff_chimeras),
               chrM          = c(lt = cutoff_chrM),
               MAPD_CNV_Log2 = c(lt = tight_cnv_mapd),
               SKEW_CNV      = c(lt = tight_cnv_sk),
               # two-sided: too little and too much conversion both mean a failed library
               Methylated_CpG_pct = c(gt = tight_methyl_min, lt = tight_methyl_max)),
    "4" = list(total_reads   = c(gt = cutoff_reads),
               PCT_CHIMERAS  = c(lt = cutoff_chimeras),
               chrM          = c(lt = cutoff_chrM),
               MAPD_CNV_Log2 = c(lt = cutoff_cnv_mapd),
               SKEW_CNV      = c(lt = cutoff_cnv_sk),
               Methylated_CpG_pct = c(gt = methyl_min, lt = methyl_max)),
    "3" = list(total_reads   = c(gt = borderline_reads),
               MAPD_CNV_Log2 = c(lt = borderline_cnv),
               SKEW_CNV      = c(lt = borderline_cnv))
  )
}

# ------------------------------------------------------------------- HiFi (basej-hifiqc)
# No read gate: HiFi yield per cell is not comparable to short-read counts and no absolute
# cutoff is established, so depth is handled by the median-relative floor instead.
qc_spec_hifi <- function(cutoff_pct5x = 80, cutoff_chimeras = 15,
                         cutoff_cnv_mapd = 0.25, cutoff_cnv_sk = 0.25,
                         tight_cnv_mapd = 0.20, tight_cnv_sk = 0.20,
                         borderline_pct5x = 50, borderline_cnv = 0.35) {
  list(
    "5" = list(pct_bases_5x       = c(gt = cutoff_pct5x),
               pct_chimeric_reads = c(lt = cutoff_chimeras),
               MAPD_CNV_Log2      = c(lt = tight_cnv_mapd),
               SKEW_CNV           = c(lt = tight_cnv_sk)),
    "4" = list(pct_bases_5x       = c(gt = cutoff_pct5x),
               pct_chimeric_reads = c(lt = cutoff_chimeras),
               MAPD_CNV_Log2      = c(lt = cutoff_cnv_mapd),
               SKEW_CNV           = c(lt = cutoff_cnv_sk)),
    "3" = list(pct_bases_5x  = c(gt = borderline_pct5x),
               MAPD_CNV_Log2 = c(lt = borderline_cnv),
               SKEW_CNV      = c(lt = borderline_cnv))
  )
}

# ------------------------------------------- RNA (basej-rnaqc, basej-rnaproteinqc)
# No read gate, per the same reasoning as HiFi; depth uses the median-relative floor on
# aligned reads. RNA has no CNV, so its tiers rest on transcriptome composition.
#
# The tier-4 defaults are basej-rnaqc's established gate values (mappability/exonic > 0.7,
# intergenic/mito < 0.1, >= 500 protein-coding genes), so `tier >= 4` reproduces the old
# "all five pass" and the PASS/Borderline/FAIL mapping already applied in main.nf holds.
qc_spec_rna <- function(cutoff_mappability = 0.7, cutoff_exonic = 0.7,
                        cutoff_intergenic = 0.1, cutoff_mito = 0.1,
                        cutoff_genes = 500,
                        tight_exonic = 0.8, tight_mito = 0.05, tight_genes = 2000,
                        borderline_genes = 250, borderline_mito = 0.2) {
  list(
    "5" = list(PropMappability = c(gt = cutoff_mappability),
               PropExonic      = c(gt = tight_exonic),
               PropIntergenic  = c(lt = cutoff_intergenic),
               ProportionCountsMitochondrialGenes = c(lt = tight_mito),
               ProteinCodingGenesDetected = c(gt = tight_genes)),
    "4" = list(PropMappability = c(gt = cutoff_mappability),
               PropExonic      = c(gt = cutoff_exonic),
               PropIntergenic  = c(lt = cutoff_intergenic),
               ProportionCountsMitochondrialGenes = c(lt = cutoff_mito),
               ProteinCodingGenesDetected = c(gt = cutoff_genes)),
    "3" = list(ProteinCodingGenesDetected = c(gt = borderline_genes),
               ProportionCountsMitochondrialGenes = c(lt = borderline_mito))
  )
}

# ---------------------------------------------------------------- WGS (basej-wgs, wgs mode)
# Full-depth short-read QC, so unlike DNA low-pass / HiFi / RNA the read count IS a
# meaningful gate here: the assay has an established target (50M) and coverage breadth is
# read against it. The original gates used `>=` for reads and coverage and `<` for the
# duplication/chimera rates; both forms are preserved exactly (see qc_cmp).
#
# Tier 5 tightens coverage breadth and duplication, which is what separates a good WGS
# library from an excellent one once the read target is met. Chimeras stay at the gate value:
# on real data the chimera gate is nearly never the binding constraint.
qc_spec_wgs <- function(cutoff_num_reads = 50e6, cutoff_pct_dup = 0.25, cutoff_pct_chim = 0.15,
                        cutoff_1x = 0.9, cutoff_5x = 0.7,
                        tight_pct_dup = 0.15, tight_1x = 0.95, tight_5x = 0.85,
                        borderline_reads = 25e6, borderline_1x = 0.8, borderline_5x = 0.5) {
  list(
    "5" = list(total_reads     = c(ge = cutoff_num_reads),
               pct_duplication = c(lt = tight_pct_dup),
               pct_chimeras    = c(lt = cutoff_pct_chim),
               pct_1x          = c(ge = tight_1x),
               pct_5x          = c(ge = tight_5x)),
    "4" = list(total_reads     = c(ge = cutoff_num_reads),
               pct_duplication = c(lt = cutoff_pct_dup),
               pct_chimeras    = c(lt = cutoff_pct_chim),
               pct_1x          = c(ge = cutoff_1x),
               pct_5x          = c(ge = cutoff_5x)),
    "3" = list(total_reads = c(ge = borderline_reads),
               pct_1x      = c(ge = borderline_1x),
               pct_5x      = c(ge = borderline_5x))
  )
}

# ---------------------------------------------------------------- WES (basej-wgs, exome mode)
# Four gates rather than five, so the old maximum score was 4 - `tier >= 4` therefore still
# means "all gates pass" and reconciles with published exome counts.
#
# fold_80_base_penalty measures how much extra sequencing would be needed to bring 80% of
# targets to the mean depth, so it is capture uniformity: lower is better and it is the metric
# that most often distinguishes a usable exome from a good one.
qc_spec_wes <- function(cutoff_num_reads = 5e5, cutoff_10x = 0.75,
                        cutoff_zero_cov = 0.05, cutoff_fold_80 = 5,
                        tight_10x = 0.90, tight_zero_cov = 0.02, tight_fold_80 = 3,
                        borderline_reads = 2.5e5, borderline_10x = 0.5, borderline_fold_80 = 8) {
  list(
    "5" = list(total_reads           = c(ge = cutoff_num_reads),
               pct_target_bases_10x  = c(ge = tight_10x),
               zero_cvg_targets_pct  = c(lt = tight_zero_cov),
               fold_80_base_penalty  = c(lt = tight_fold_80)),
    "4" = list(total_reads           = c(ge = cutoff_num_reads),
               pct_target_bases_10x  = c(ge = cutoff_10x),
               zero_cvg_targets_pct  = c(lt = cutoff_zero_cov),
               fold_80_base_penalty  = c(lt = cutoff_fold_80)),
    "3" = list(total_reads           = c(ge = borderline_reads),
               pct_target_bases_10x  = c(ge = borderline_10x),
               fold_80_base_penalty  = c(lt = borderline_fold_80))
  )
}
