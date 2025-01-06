context("samples")

test_that("Associate components", {
    set.seed(2)
    airway <- .getAirwayData(n_features = 500)

    set.seed(1)
    airway_fe <- estimateFactors(airway, nc = 2)

    fe_res <- associateComponents(airway_fe, "~ cell + dex")

    expect_true(all(paste0("factor_", 1:2) %in% fe_res$anovas$component))
    expect_true(all(paste0("factor_", 1:2) %in% fe_res$summaries$component))
    expect_true(all(paste0("factor_", 1:2) %in% names(fe_res$models)))
    expect_true(inherits(fe_res$models[[1]], "lm"))

    airway_me <- identifyModules(airway, verbose = 0, power = 21)
    me_res <- associateComponents(airway_me, "~ dex + (1|cell)", method = "lmer")

    expect_true(all(paste0("module_", 0:5) %in% me_res$anovas$component))
    expect_true(all(paste0("module_", 0:5) %in% me_res$summaries$component))
    expect_true(all(paste0("module_", 0:5) %in% names(me_res$models)))
    expect_true(inherits(me_res$models[[1]], "lmerModLmerTest"))
})
