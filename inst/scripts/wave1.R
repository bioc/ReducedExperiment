library(SummarizedExperiment)
library(biomaRt)
library(data.table)
library(edgeR)
library(M3C)

# These functions may be found on GitHub: https://github.com/jackgisby/covid-longitudinal-multi-omics/tree/main/R
# Repository archived at Zenodo: https://zenodo.org/records/7333789
source("R/preprocessing.R")
source("R/de.R")

# These datasets may be downloaded from Zenodo: https://zenodo.org/records/6497251
htseq_counts <- data.frame(data.table::fread("htseq_counts.csv"))
wave1_counts <- get_summarized_experiment(htseq_counts, "w1_metadata.csv")
wave2_counts <- get_summarized_experiment(htseq_counts, "w2_metadata.csv")

# Functions found in aformentioned repository - takes counts and outputs TMM-normalised logCPM
wave1_se <- normalize_se(
    wave1_counts,
    normalisation_level = "logcpm",
    group = "case_control"
)

wave2_se <- normalize_se(
    wave2_counts,
    normalisation_level = "logcpm",
    group = "case_control"
)

saveRDS(wave1_se, "inst/extdata/wave1.rds", compress = "gzip")
saveRDS(wave2_se, "inst/extdata/wave2.rds", compress = "gzip")
