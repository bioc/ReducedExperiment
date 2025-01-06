library(ReducedExperiment)

set.seed(2)
airway <- ReducedExperiment:::.getAirwayData(n_features = 2000)
airway_fe <- estimateFactors(
    airway,
    nc = 2,
    use_stability = FALSE,
    method = "imax"
)

overrep_res <- runEnrich(
    airway_fe,
    method = "overrepresentation",
    feature_id_col = "rownames",
    as_dataframe = TRUE,
    p_cutoff = 0.1,
    universe = rownames(airway_fe)
)

t2g <- ReducedExperiment::getMsigdbT2G()
t2g <- t2g[t2g$gs_name %in% overrep_res$ID, ]
t2g <- rbind(t2g, data.frame(gs_name = "universe", ensembl_gene = rownames(airway_fe)))
write.csv(t2g, "inst/extdata/msigdb_t2g_filtered.csv", row.names = FALSE)
