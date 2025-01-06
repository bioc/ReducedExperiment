#' Helper function to run a linear mixed model
#'
#' @importFrom lmerTest lmer
#' @importFrom lme4 lmerControl
#'
#' @author Jack Gisby
#'
#' @noRd
#' @keywords internal
.singleLmer <- function(data, formula_string, REML = TRUE, ...) {
    out_model <- tryCatch(
        lmerTest::lmer(
            stats::as.formula(formula_string),
            data = data,
            REML = REML,
            control = lme4::lmerControl(check.conv.singular = "ignore"),
            ...
        ),
        warning = function(w) {
            return(lmerTest::lmer(
                stats::as.formula(formula_string),
                data = data,
                REML = REML,
                control = lme4::lmerControl(
                    optimizer = "Nelder_Mead",
                    check.conv.singular = "ignore"
                ),
                ...
            ))
        }
    )

    if (is(out_model, "lmerModLmerTest")) {
        return(out_model)
    } else {
        stop("Convergence issue not caught by single_lmer")
    }
}

#' Helper function to run a linear model
#'
#' @author Jack Gisby
#'
#' @noRd
#' @keywords internal
.runLinearModel <- function(
    X,
    pheno,
    formula,
    method = "lm",
    type = "II",
    ...
) {
    pheno$component <- X

    formula <- stats::as.formula(paste0("component", formula))

    if (method == "lmer") {
        linear_model <- .singleLmer(pheno, formula, ...)

        anova_res <- data.frame(stats::anova(linear_model, type = type))
        summary_res <- data.frame(summary(linear_model)$coefficients)
    } else if (method == "lm") {
        linear_model <- stats::lm(formula, pheno)

        anova_res <- data.frame(car::Anova(linear_model, type = type))
        summary_res <- data.frame(summary(linear_model)$coefficients)
    }

    anova_res$term <- rownames(anova_res)
    summary_res$term <- rownames(summary_res)

    return(list(
        "model" = linear_model,
        "anova" = anova_res,
        "summary" = summary_res
    ))
}

#' Runs linear models for components and sample-level data
#'
#' Runs either standard linear or linear mixed models, with reduced components
#' (e.g., factors or modules) as the outcomes  and sample-level information
#' (e.g., treatment, disease status) as predictors.
#'
#' @param re An object inheriting from
#' \link[ReducedExperiment]{ReducedExperiment}.
#'
#' @param formula The model formula to apply. Only the right hand side of the
#' model need be specified (e.g., "~ x + y"). The left hand side (outcome)
#' represents the
#' components themselves. The variables in this formula should be present in the
#' `colData` of `re`.
#'
#' @param method If "lm", then the \link[stats]{lm} function is used to run
#' linear models (in tandem with \link[car]{Anova} for running anovas on the
#' model terms). If "lmer", then linear mixed models are run through
#' \link[lmerTest]{lmer}.
#'
#' @param scale_reduced If `TRUE`, the reduced data are scaled (to have a standard
#' deviation of 1) before modelling.
#'
#' @param center_reduced If `TRUE`, the reduced data are centered (to have a mean
#' of 0) before modelling.
#'
#' @param type The type of anova to be applied to the terms of the linear model.
#'
#' @param adj_method The method for adjusting for multiple testing. Passed to
#' the \link[stats]{p.adjust} `method` parameter.
#'
#' @param ... Additional arguments passed to \link[lmerTest]{lmer}, given that
#' `method` is set to "lmer".
#'
#' @returns Returns a list with the entry "models" including a list of the
#' model objects, "anovas" containing the output of anova-based testing,
#' and "summaries" containing the results of running `summary` on the models.
#'
#' @details
#' Multiple testing adjustment is performed separately for each term in the
#' model across all factors.
#' In other words, p-values are adjusted for the number of factors, but not
#' the number of model terms. If you are testing a large number of terms,
#' you could consider applying a custom adjustment method or using
#' penalised regression.
#'
#' @seealso [stats::lm()], [car::Anova()], [lmerTest::lmer()]
#'
#' @author Jack Gisby
#'
#' @examples
#' # Create FactorisedExperiment with random data (100 features, 50 samples,
#' # 10 factors)
#' set.seed(1)
#' fe <- ReducedExperiment:::.createRandomisedFactorisedExperiment(100, 50, 10)
#' fe
#'
#' # Create a sample-level variable describing some sort of treatment
#' colData(fe)$treated <- c(rep("control", 25), rep("treatment", 25))
#' colData(fe)$treated <- factor(colData(fe)$treated, c("control", "treatment"))
#'
#' # Increase the value of factor 1 for the treated samples, simulating some
#' # kind of treatment-related effect picked up by factor analysis
#' reduced(fe)[, 1][colData(fe)$treated == "treatment"] <-
#'     reduced(fe)[, 1][colData(fe)$treated == "treatment"] +
#'     rnorm(25, mean = 1.5, sd = 0.1)
#'
#' # Create a sample-level variable describing a covariate we want to adjust for
#' # We will make the treated patients slightly older on average
#' colData(fe)$age <- 0
#' colData(fe)$age[colData(fe)$treated == "control"] <- rnorm(25, mean = 35, sd = 8)
#' colData(fe)$age[colData(fe)$treated == "treatment"] <- rnorm(25, mean = 40, sd = 8)
#'
#' # Associate the factors with sample-level variable in the colData
#' lm_res <- associateComponents(
#'     fe,
#'     formula = "~ treated + age", # Our model formula
#'     method = "lm", # Use a linear model
#'     adj_method = "BH" # Adjust our p-values with Benjamini-Hochberg
#' )
#'
#' # We see that treatment is significantly associated with factor 1 (adjusted
#' # p-value < 0.05) and is higher in the treated patients. Age is not
#' # significantly associated with factor 1, but there is a slight positive
#' # relationship
#' print(head(lm_res$summaries[
#'     ,
#'     c("term", "component", "estimate", "stderr", "pvalue", "adj_pvalue")
#' ]))
#'
#' # But what if these aren't 50 independent patients, but rather 25 patients
#' # sampled before and after treatment? We can account for this using a
#' # linear mixed model, which can account for repeated measures and paired
#' # designs
#'
#' # First we add in this information
#' colData(fe)$patient_id <- c(paste0("patient_", 1:25), paste0("patient_", 1:25))
#'
#' # Then we run the linear mixed model with a random intercept for patient
#' lmm_res <- associateComponents(
#'     fe,
#'     formula = "~ treated + age + (1 | patient_id)", # Add a random intercept
#'     method = "lmer", # Use a linear mixed model
#'     adj_method = "BH"
#' )
#'
#' # We used a different method, but can obtain a similar summary output
#' print(head(lmm_res$summaries[
#'     ,
#'     c("term", "component", "estimate", "stderr", "pvalue", "adj_pvalue")
#' ]))
#'
#' @export
associateComponents <- function(
    re,
    formula,
    method = "lm",
    scale_reduced = TRUE,
    center_reduced = TRUE,
    type = "II",
    adj_method = "BH",
    ...
) {
    models <- list()
    summaries <- anovas <- data.frame()

    red <- reduced(re,
        scale_reduced = scale_reduced,
        center_reduced = center_reduced
    )

    for (comp in componentNames(re)) {
        linear_model <- .runLinearModel(
            X = red[, comp],
            pheno = data.frame(colData(re)),
            formula = formula,
            method = method,
            type = type,
            ...
        )

        linear_model$anova$component <- linear_model$summary$component <- comp

        models[[comp]] <- linear_model$model
        anovas <- rbind(anovas, linear_model$anova)
        summaries <- rbind(summaries, linear_model$summary)
    }

    colnames(anovas) <- .renameResultsTable(colnames(anovas))
    colnames(summaries) <- .renameResultsTable(colnames(summaries))

    rownames(anovas) <- paste(anovas$component, anovas$term, sep = "_")
    rownames(summaries) <- paste(summaries$component, summaries$term, sep = "_")

    anovas$adj_pvalue <- .adjustByTerm(anovas, method = adj_method)
    summaries$adj_pvalue <- .adjustByTerm(summaries, method = adj_method)

    return(list("models" = models, "anovas" = anovas, "summaries" = summaries))
}

#' Adjusts the p-values of model results for multiple testing per-term
#'
#' @author Jack Gisby
#'
#' @noRd
#' @keywords internal
.adjustByTerm <- function(res, method = "BH") {
    res$adj_pvalue <- NA

    for (term in unique(res$term)) {
        res$adj_pvalue[which(res$term == term)] <-
            stats::p.adjust(res$pvalue[which(res$term == term)],
                method = method
            )
    }

    return(res$adj_pvalue)
}

#' Renames linear model result table column names
#'
#' @author Jack Gisby
#'
#' @noRd
#' @keywords internal
.renameResultsTable <- function(cnames) {
    cname_conversions <- list(
        "Sum.Sq" = "sum_sq",
        "Mean.Sq" = "mean_sq",
        "NumDF" = "num_df",
        "DenDF" = "den_df",
        "F.value" = "fvalue",
        "Pr..F." = "pvalue",
        "Estimate" = "estimate",
        "Std..Error" = "stderr",
        "t.value" = "tvalue",
        "Pr...t.." = "pvalue"
    )

    for (i in seq_along(cnames)) {
        if (cnames[i] %in% names(cname_conversions)) {
            cnames[i] <- cname_conversions[[cnames[i]]]
        }
    }

    return(cnames)
}
