#!/usr/bin/env Rscript
# Correct the 5'/3'/5'-3' bias columns in qualimap_stats_mqc.csv.
#
# parse_qualimap.R matches each Qualimap metric with an unanchored grep. Because
# "3' bias" and "5'-3' bias" are substrings of the earlier Qualimap header line
# "5'-3' bias region size = 100", both of those columns get populated with the
# constant 100 instead of the real per-sample values (only "5' bias" is
# unambiguous and therefore correct). This helper re-reads the real values
# straight from each sample's rnaseq_qc_results.txt using start-anchored matches
# and overwrites only those three columns.
#
# Everything else is preserved: all other columns are read and written back
# verbatim (colClasses="character", na.strings none), and the output keeps the
# unquoted format that write.csv(quote=FALSE) in parse_qualimap.R produces, so
# the downstream RNA_QC_PLOTS parse is unaffected.
#
# The root cause lives in parse_qualimap.R. It is corrected here, in pipeline
# code, to avoid a container rebuild and to keep that shared script
# byte-identical to its nf-scrnaseq-pipeline source copy.

csv <- "qualimap_stats_mqc.csv"
if (!file.exists(csv)) quit(save = "no", status = 0)

df <- read.csv(csv, check.names = FALSE, colClasses = "character",
               na.strings = character(0))

# Mirror the biosampleName derivation used in main.nf's RNA_QC_PLOTS block.
biosample <- function(x) {
  sub("_Aligned\\.sortedByCoord\\.out$", "",
      sub("\\.(bam|cram)$", "", basename(x)))
}

# Start-anchored extraction. "^\\s*<label>\\s*=" cannot match the
# "5'-3' bias region size = 100" / "... number of top transcripts = 1000"
# header lines (those have text, not "=", immediately after the label), and the
# "3' bias" pattern cannot match a "5'-3' bias" line because the anchored label
# does not start with "5".
getval <- function(lines, label) {
  pat <- paste0("^[[:space:]]*", label, "[[:space:]]*=[[:space:]]*")
  hit <- grep(pat, lines, value = TRUE)
  if (!length(hit)) return(NA_character_)
  trimws(sub(pat, "", hit[1]))
}

labels <- c("5' bias", "3' bias", "5'-3' bias")

if (!"bam file" %in% colnames(df)) {
  message("WARNING: 'bam file' column absent from qualimap_stats_mqc.csv; ",
          "leaving bias columns unchanged")
  quit(save = "no", status = 0)
}

for (i in seq_len(nrow(df))) {
  s <- biosample(df[["bam file"]][i])
  f <- file.path(paste0("qualimap_outdir_", s), "rnaseq_qc_results.txt")
  if (!file.exists(f)) {
    message("WARNING: ", f, " not found; leaving bias values for '", s, "' unchanged")
    next
  }
  lines <- readLines(f, warn = FALSE)
  for (lab in labels) {
    if (!lab %in% colnames(df)) next
    v <- getval(lines, lab)
    if (!is.na(v)) df[i, lab] <- v
  }
}

write.csv(df, csv, row.names = FALSE, quote = FALSE)
cat("Corrected 5'/3'/5'-3' bias columns in qualimap_stats_mqc.csv from rnaseq_qc_results.txt\n")
