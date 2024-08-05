library(SummarizedExperiment)
library(biomaRt)
library(data.table)
library(dplyr)
library(edgeR)
library(M3C)

# These functions may be found on GitHub: https://github.com/jackgisby/covid-longitudinal-multi-omics/tree/main/R
# Repository archived at Zenodo: https://zenodo.org/records/7333789
source("R/preprocessing.R")
source("R/de.R")

# For reproducibility
set.seed(1)

# These datasets may be downloaded from Zenodo: https://zenodo.org/records/6497251
htseq_counts <- data.frame(data.table::fread("htseq_counts.csv"))

mart <- biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
wave1_counts <- get_summarized_experiment(htseq_counts, "w1_metadata.csv", mart = mart)
wave2_counts <- get_summarized_experiment(htseq_counts, "w2_metadata.csv", mart = mart)

# Functions found in aformentioned repository
# Takes counts and outputs TMM-normalised logCPM
# Also filters by counts and variance
wave1_se <- normalize_se(
    wave1_counts,
    normalisation_level = "variance_filtered",
    variance_filter_cutoff = 3000,
    group = "case_control"
)

wave2_se <- normalize_se(
    wave2_counts,
    normalisation_level = "variance_filtered",
    variance_filter_cutoff = 3000,
    group = "case_control"
)

# Keep genes that are shared in both SEs
wave1_se <- wave1_se[rownames(wave1_se) %in% rownames(wave2_se), ]
wave2_se <- wave2_se[rownames(wave2_se) %in% rownames(wave1_se), ]
wave2_se <- wave2_se[match(rownames(wave1_se), rownames(wave2_se)), ]
stopifnot(all(rownames(wave1_se) == rownames(wave2_se)))

# Do an additional step to randomly filter out genes
genes_to_keep <- sample(rownames(wave1_se), 500)
wave1_se <- wave1_se[which(rownames(wave1_se) %in% genes_to_keep), ]
wave2_se <- wave2_se[which(rownames(wave2_se) %in% genes_to_keep), ]
stopifnot(all(rownames(wave1_se) == rownames(wave2_se)))

# Just keep the assay we need
assay(wave1_se, "normal") <- assay(wave1_se, "cpm")
assay(wave2_se, "normal") <- assay(wave2_se, "cpm")
assay(wave1_se, "counts") <- assay(wave2_se, "counts") <- NULL
assay(wave1_se, "cpm") <- assay(wave2_se, "cpm") <- NULL

# Remove late samples
wave1_se <- wave1_se[, which(!wave1_se$sample_id %in% wave1_se$sample_id[which(wave1_se$time_from_first_x > 15)])]
wave2_se <- wave2_se[, which(!wave2_se$sample_id %in% wave2_se$sample_id[which(wave2_se$time_from_first_x > 15)])]

# Order samples by severity
wave1_se$WHO_temp_severity <- factor(wave1_se$WHO_temp_severity, c("mild", "moderate", "severe", "critical"))
wave2_se$WHO_temp_severity <- factor(wave2_se$WHO_temp_severity, c("mild", "moderate", "severe", "critical"))
wave1_se <- wave1_se[, order(wave1_se$WHO_temp_severity, decreasing = TRUE)]
wave2_se <- wave2_se[, order(wave2_se$WHO_temp_severity, decreasing = TRUE)]

# Keep the most severe sample from each individual in the Wave 1 cohort
wave1_se <- wave1_se[, which(!duplicated(wave1_se$individual_id))]

# Keep the three most severe samples from each individual in the Wave 2 cohort
samples_to_keep <- data.frame(colData(wave2_se)) %>%
    group_by(individual_id) %>%
    slice_head(n = 3) %>%
    ungroup()

wave2_se <- wave2_se[, which(wave2_se$case_control == "NEGATIVE" | wave2_se$sample_id %in% samples_to_keep$sample_id)]

# Remove duplicate controls from the wave 1 cohort
wave1_se <- wave1_se[, !wave1_se$individual_id %in% wave2_se$individual_id]

# Remove the recovery samples from the wave 2 cohort
wave2_se <- wave2_se[, which(wave2_se$case_control != "RECOVERY")]

# Remove unnecessary columns from row and colData
colData(wave1_se) <- colData(wave1_se)[, c("sample_id", "individual_id", "sex", "calc_age", "WHO_severity", "WHO_temp_severity", "case_control", "time_from_first_x")]
colData(wave2_se) <- colData(wave2_se)[, c("sample_id", "individual_id", "sex", "calc_age", "WHO_severity", "WHO_temp_severity", "case_control", "time_from_first_x")]
rowData(wave1_se)$gencode_id <- NULL
rowData(wave2_se)$gencode_id <- NULL

# Save the objects
saveRDS(wave1_se, "scripts/20240305_se_from_public/wave1.rds")
wave1_se
saveRDS(wave2_se, "scripts/20240305_se_from_public/wave2.rds")
wave2_se
