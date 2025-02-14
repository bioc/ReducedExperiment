#' FactorisedExperiment: A container for the results of factor analysis
#'
#' @description
#' A container inheriting from the \link[ReducedExperiment]{ReducedExperiment}
#' class, that
#' contains one or more data matrices, to which factor analysis has been applied
#' to identify a reduced set of features. A
#' \link[ReducedExperiment]{FactorisedExperiment} can be created directly in
#' a similar manner to a \link[SummarizedExperiment]{SummarizedExperiment}.
#' Alternatively, the \link[ReducedExperiment]{estimateFactors} function
#' can be used to both apply factor analysis and generate a
#' \link[ReducedExperiment]{FactorisedExperiment} from the results.
#'
#' @param reduced A `matrix`, produced by factor analysis, with rows
#' representing samples and columns representing factors.
#'
#' @param loadings A `matrix`, produced by factor analysis, with rows
#' representing features and columns representing factors.
#'
#' @param stability A vector containing some measure of stability or variance
#' explained for each factor. If factor analysis was performed using
#' \link[ReducedExperiment]{estimateFactors} and `use_stability = TRUE`, this
#' slot will indicate the stability of the factors across multiple runs of ICA.
#'
#' @param ... Additional arguments to be passed to
#' \link[ReducedExperiment]{ReducedExperiment}.
#'
#' @inheritParams ReducedExperiment
#'
#' @returns Constructor method returns a
#' \link[ReducedExperiment]{FactorisedExperiment} object.
#'
#' @seealso [ReducedExperiment::ReducedExperiment()],
#' [ReducedExperiment::ModularExperiment()],
#' [ReducedExperiment::estimateFactors()]
#'
#' @author Jack Gisby
#'
#' @examples
#' # Create randomised data with the following dimensions
#' i <- 300 # Number of features
#' j <- 100 # Number of samples
#' k <- 10 # Number of factors
#'
#' # In this case we use random assay, reduced and loadings data, but in
#' # practice these will likely be the result of applying some kind of factor
#' # analysis to the assay data (e.g., gene expression data) from some study.
#' rand_assay_data <- ReducedExperiment:::.makeRandomData(i, j, "gene", "sample")
#' rand_reduced_data <- ReducedExperiment:::.makeRandomData(j, k, "sample", "factor")
#' rand_loadings <- ReducedExperiment:::.makeRandomData(i, k, "gene", "factor")
#'
#' fe <- FactorisedExperiment(
#'     assays = list("normal" = rand_assay_data),
#'     reduced = rand_reduced_data,
#'     loadings = rand_loadings
#' )
#'
#' fe
#'
#' @import SummarizedExperiment
#'
#' @rdname factorised_experiment
#' @export
FactorisedExperiment <- function(
    reduced = new("matrix"),
    scale = TRUE,
    center = TRUE,
    loadings = new("matrix"),
    stability = NULL,
    ...
) {
    re <- ReducedExperiment(
        reduced = reduced,
        scale = scale,
        center = center,
        ...
    )

    return(.FactorisedExperiment(
        re,
        loadings = loadings,
        stability = stability
    ))
}

S4Vectors::setValidity2("FactorisedExperiment", function(object) {
    msg <- NULL

    obj_dims <- dim(object)

    # Check feature names/numbers
    if (obj_dims[1] != dim(loadings(object))[1]) {
        msg <- c(msg, "Loadings have invalid row dimensions")
    }

    if (!identical(featureNames(object), rownames(loadings(object)))) {
        msg <- c(msg, "Loadings have incorrect column names (feature labels)")
    }

    # Factors
    if (dim(loadings(object))[2] != dim(reduced(object))[2]) {
        msg <- c(
            msg,
            "Reduced data and loadings have incompatible column dimensions"
        )
    }

    if (!identical(colnames(loadings(object)), colnames(reduced(object)))) {
        msg <- c(
            msg,
            paste0(
                "Reduced data and loadings have incompatible column names ",
                "(factor names)"
            )
        )
    }

    # Stability - check names/length matches
    if (!is.null(stability(object)) & length(stability(object)) > 0) {
        if (length(stability(object)) != nComponents(object)) {
            msg <- c(
                msg,
                "Number of components do not match with component stability"
            )
        }

        # If stability vector has names, check they are correct
        if (!is.null(names(stability(object)))) {
            if (!identical(names(stability(object)), componentNames(object))) {
                msg <- c(
                    msg,
                    "Component names do not match with component stability"
                )
            }
        }
    }

    return(if (is.null(msg)) TRUE else msg)
})

#' @rdname loadings
#' @export
setMethod("loadings", "FactorisedExperiment", function(
    object,
    scale_loadings = FALSE,
    center_loadings = FALSE,
    abs_loadings = FALSE
) {
    l <- scale(
        object@loadings,
        scale = scale_loadings,
        center = center_loadings
    )
    if (abs_loadings) l <- abs(l)
    return(l)
})

#' @rdname loadings
#' @export
setReplaceMethod("loadings", "FactorisedExperiment", function(object, value) {
    object@loadings <- value
    validObject(object)
    return(object)
})

#' @rdname feature_names
#' @export
setReplaceMethod("names", "FactorisedExperiment", function(x, value) {
    rownames(x@loadings) <- value
    x <- callNextMethod(x, value)
    validObject(x)
    return(x)
})

#' @rdname feature_names
#' @export
setReplaceMethod("featureNames", "FactorisedExperiment", function(x, value) {
    names(x) <- value
    return(x)
})

#' @rdname feature_names
#' @export
setReplaceMethod("rownames", "FactorisedExperiment", function(x, value) {
    names(x) <- value
    return(x)
})

#' Get and setting the stability values for factors
#'
#' @param object \link[ReducedExperiment]{FactorisedExperiment} object.
#'
#' @param value New value to replace existing stability vector.
#'
#' @returns A vector with a value for each factor indicating the factor
#' stability. More details are available from the
#' \link[ReducedExperiment]{estimateStability} help page.
#'
#' @seealso [ReducedExperiment::estimateStability()]
#'
#' @author Jack Gisby
#'
#' @examples
#' # Get a random matrix with rnorm, with 100 rows (features)
#' # and 20 columns (observations)
#' X <- ReducedExperiment:::.makeRandomData(100, 20, "feature", "obs")
#'
#' # Run stabilised ICA on the data with 5 components
#' fe <- estimateFactors(X, nc = 5, use_stability = TRUE)
#'
#' stability(fe)
#'
#' stability(fe)[2] <- 10
#' stability(fe)
#'
#' @rdname stability
#' @name stability
#' @aliases stability<-
#' @export stability
NULL

#' @rdname stability
#' @export
setMethod("stability", "FactorisedExperiment", function(object) {
    return(object@stability)
})

#' @rdname stability
#' @export
setReplaceMethod("stability", "FactorisedExperiment", function(object, value) {
    object@stability <- value
    validObject(object)
    return(object)
})

#' @rdname component_names
#' @export
setReplaceMethod("componentNames", "FactorisedExperiment", function(
        object,
        value) {
    colnames(object@loadings) <- value
    if (!is.null(object@stability)) names(object@stability) <- value
    object <- callNextMethod(object, value)
    validObject(object)
    return(object)
})

#' @rdname dollar_names
#' @export
.DollarNames.FactorisedExperiment <- function(x, pattern = "") {
    grep(pattern, names(colData(x)), value = TRUE)
}

#' @rdname slice
#' @export
setMethod(
    "[", c("FactorisedExperiment", "ANY", "ANY", "ANY"),
    function(x, i, j, k, ..., drop = FALSE)
{
    object <- x

    if (1L != length(drop) || (!missing(drop) && drop)) {
        warning("'drop' ignored '[,", class(object), ",ANY,ANY-method'")
    }

    lod <- object@loadings
    stab <- object@stability

    if (!missing(i)) {
        i <- .process_char_index(class(object), rownames(object), i, "i")
        lod <- lod[i, , drop = FALSE]
    }

    if (!missing(k)) {
        k <- .process_char_index(class(object), componentNames(object), k, "k")
        lod <- lod[, k, drop = FALSE]
        stab <- stab[k, drop = FALSE]
    }

    out <- callNextMethod(object, i, j, k, ...)
    BiocGenerics:::replaceSlots(out,
        loadings = lod, stability = stab,
        check = FALSE
    )
})

#' @rdname slice
#' @export
setReplaceMethod(
    "[",
    signature(x = "FactorisedExperiment", value = "FactorisedExperiment"),
    function(x, i, j, k, ..., value)
{
    if (missing(i) & missing(j) & missing(k)) {
        return(value)
    }

    object <- x
    lod <- object@loadings
    stab <- object@stability

    if (!missing(i)) {
        i <- .process_char_index(class(object), rownames(object), i, "i")
    } else {
        i <- seq_len(nrow(object))
    }

    if (!missing(k)) {
        k <- .process_char_index(class(object), componentNames(object), k, "k")
    } else {
        k <- seq_len(nComponents(object))
    }

    stab[k] <- value@stability
    lod[i, k] <- value@loadings

    out <- callNextMethod(object, i, j, k, ..., value = value)
    BiocGenerics:::replaceSlots(
        out,
        loadings = lod,
        stability = stab,
        check = FALSE
    )
})

# Same features/compnames, different samples
#' @rdname cbind_rbind
#' @export
setMethod("cbind", "FactorisedExperiment", function(..., deparse.level = 1) {
    args <- list(...)

    loadings_stability_equal <- vapply(args, function(re) {
        return(identical(re@loadings, args[[1]]@loadings) &
            identical(re@stability, args[[1]]@stability))
    },
    FUN.VALUE = FALSE
    )

    if (!all(loadings_stability_equal)) {
        stop("Column bind expects loadings and stability slots are equal")
    } else {
        args[["deparse.level"]] <- deparse.level
        return(do.call(callNextMethod, args))
    }
})

# Same samples/compnames, different features
#' @rdname cbind_rbind
#' @export
setMethod("rbind", "FactorisedExperiment", function(..., deparse.level = 1) {
    args <- list(...)

    stability_equal <- vapply(args, function(re) {
        return(identical(re@stability, args[[1]]@stability))
    },
    FUN.VALUE = FALSE
    )

    if (!all(stability_equal)) {
        stop("Row bind expects stability slots are equal")
    } else {
        loadings_slot <- do.call(rbind, lapply(args, loadings))

        args[[1]] <- BiocGenerics:::replaceSlots(
            args[[1]],
            loadings = loadings_slot,
            check = FALSE
        )

        args[["deparse.level"]] <- deparse.level
        return(do.call(callNextMethod, args))
    }
})

#' Project new data using pre-defined factors
#'
#' @description
#' Uses a projection approach to calculate factors in new data. Functions in a
#' similar fashion to the `predict` method of \link[stats]{prcomp}. The
#' transposed `newdata` are multiplied by the original loadings matrix.
#'
#' @param object A \link[ReducedExperiment]{FactorisedExperiment} object. The
#' `loadings` slot of this class will be used for projection. Additionally, by
#' default, the `scale` and `center` slots are used to apply the original
#' transformation to the new data.
#'
#' @param newdata New data for projection. Must be a `data.frame` or `matrix`
#' with features as rows and samples as columns, or a
#' \link[SummarizedExperiment]{SummarizedExperiment} object. Assumes that the
#' rows of `newdata` match those of the
#' \link[ReducedExperiment]{FactorisedExperiment} object.
#'
#' @param standardise_reduced Whether or not the reduced data should be standardised
#' (i.e., transformed to have a mean of 0 and standard deviation of 1)
#' after calculation.
#'
#' @param scale_newdata Controls whether the `newdata` are scaled. If `NULL`,
#' performs scaling based on the `FactorisedExperiment`
#' object's `scale` slot. The value of this argument will be passed to the
#' `scale` argument of \link[base]{scale}.
#'
#' @param center_newdata Controls whether the `newdata` are centered If `NULL`,
#' performs centering based on the
#' `FactorisedExperiment` object's `center` slot. The
#' value of this argument will be passed to the `center` argument of
#' \link[base]{scale}.
#'
#' @param assay_name If a \link[SummarizedExperiment]{SummarizedExperiment}
#' object is passed as new data, this argument indicates which assay should be
#' used for projection.
#'
#' @param ... Additional arguments to be passed to `projectData.`
#'
#' @returns Calculates a matrix with samples as rows and factors as columns. If
#' `newdata` was a `matrix` or `data.frame`, this will be returned as a
#' `matrix`.
#' If a \link[SummarizedExperiment]{SummarizedExperiment} object was passed
#' instead, then a `FactorisedExperiment`
#' object will be created containing this `matrix` in its `reduced` slot.
#'
#' @details
#' If `scale_newdata` and `center_newdata` are left as `NULL`, then the
#' projection method assumes that the `newdata` are on the same scale as the
#' original data of the `object`. It will therefore use the values of the
#' `center` and `scale` slots of the `object`. For instance, if the `scale` slot
#' is `TRUE`, the `newdata` will be scaled. If the `scale` slot is a vector,
#' the values of this vector will be applied to scale the `newdata`.
#'
#' @seealso \code{\link[ReducedExperiment]{calcEigengenes}}, [stats::prcomp]
#'
#' @author Jack Gisby
#'
#' @examples
#' # Get two random matrices with rnorm
#' # 1: 100 rows (features) and 20 columns (observations)
#' X_1 <- ReducedExperiment:::.makeRandomData(100, 20, "feature", "obs")
#'
#' # Both matrices must have the same features, but they may have different obs
#' # 2: 100 rows (features) and 30 columns (observations)
#' X_2 <- ReducedExperiment:::.makeRandomData(100, 30, "feature", "obs")
#'
#' # Estimate 5 factors based on the data matrix
#' fe_1 <- estimateFactors(X_1, nc = 5)
#' fe_1
#'
#' # Project the fe_1 factors for the samples in X_2
#' projected_data <- projectData(fe_1, X_2)
#' projected_data
#'
#' @rdname projectData
#' @name projectData
#' @export projectData
NULL

#' @rdname projectData
#' @export
setMethod("projectData", c("FactorisedExperiment", "matrix"),  function(
    object,
    newdata,
    standardise_reduced = TRUE,
    scale_newdata = NULL,
    center_newdata = NULL
) {
    if (!identical(rownames(object), rownames(newdata))) {
        stop("Rownames of x do not match those of newdata")
    }

    # apply known vectors for scaling and centering (returned as attributes
    # by `scale`)
    if (is.null(scale_newdata)) scale_newdata <- object@scale
    if (is.null(center_newdata)) center_newdata <- object@center

    newdata <- t(scale(t(newdata),
        scale = scale_newdata,
        center = center_newdata
    ))
    red <- .projectICA(newdata, loadings(object))

    if (standardise_reduced) red <- scale(red)

    return(red)
})

#' @rdname projectData
#' @export
setMethod("projectData", c("FactorisedExperiment", "data.frame"), function(
    object,
    newdata,
    standardise_reduced = TRUE,
    scale_newdata = NULL,
    center_newdata = NULL
) {
    return(projectData(
        object,
        as.matrix(newdata),
        standardise_reduced = standardise_reduced,
        scale_newdata = scale_newdata,
        center_newdata = center_newdata
    ))
})

#' @rdname projectData
#' @export
setMethod(
    "projectData",
    c("FactorisedExperiment", "SummarizedExperiment"),
    function(
        object,
        newdata,
        standardise_reduced = TRUE,
        scale_newdata = NULL,
        center_newdata = NULL,
        assay_name = "normal"
) {
    projected_data <- projectData(
        object,
        assay(newdata, assay_name),
        standardise_reduced = standardise_reduced,
        scale_newdata = scale_newdata,
        center_newdata = center_newdata
    )

    return(.seToFe(
        newdata,
        reduced = projected_data,
        loadings = loadings(object),
        stability = stability(object),
        center_X = object@center,
        scale_X = object@scale
    ))
})

#' @rdname projectData
#' @export
setMethod("predict", c("FactorisedExperiment"), function(object, newdata, ...) {
    return(projectData(object, newdata, ...))
})

#' Get feature alignments with factors
#'
#' Retrieves features (usually genes) and their alignment (`loadings`) with the
#' factors. Allows for the selection of features whose alignments are high
#' relative to other features. Useful for functional interpretation of factors.
#'
#' @param object A \link[ReducedExperiment]{FactorisedExperiment} object.
#'
#' @param loading_threshold A value between 0 and 1 indicating the proportion of
#' the maximal loading to be used as a threshold. A value of 0.5 (default) means
#' that genes will be selected if their factor alignment
#' (derived from the `loadings` slot) exceeds or equals 50% of the maximally
#' aligned feature.
#'
#' @param proportional_threshold A value between 0 and 1 indicating the maximal
#' proportion of features to be returned. A value of 0.01 (default) means that a
#' maximum of 1% of the input features (usually genes) will be returned for
#' each factor. These will be the genes in the top percentile with respect to
#' the `loadings`
#'
#' @param feature_id_col The column in `rowData(object)` that will be used as a
#' feature ID. Setting this to "rownames" (default) instead uses
#' `rownames(object)`.
#'
#' @param format A string specifying the format in which to return the results.
#' See the `value` section below.
#'
#' @param center_loadings If `TRUE`, loadings will be centered column-wise to
#' have a mean of 0.
#'
#' @returns If the `format` argument is "list", then a
#' list will be returned with an entry for each factor, each containing a vector
#' of input features. Otherwise, if `format` is `"data.frame"`, a data.frame is
#' returned with a row for each gene-factor combination. The `format` argument
#' can also be a function to be applied to the output data.frame before
#' returning the results.
#'
#' @seealso [ReducedExperiment::getCommonFeatures()]
#'
#' @author Jack Gisby
#'
#' @examples
#' # Get a random matrix with rnorm, with 100 rows (features)
#' # and 20 columns (observations)
#' X <- ReducedExperiment:::.makeRandomData(100, 20, "feature", "obs")
#'
#' # Estimate 5 factors based on the data matrix
#' fe <- estimateFactors(X, nc = 5)
#'
#' # Get the genes highly aligned with each factor as a list
#' aligned_features <- getAlignedFeatures(fe, proportional_threshold = 0.03)
#' aligned_features
#'
#' # Can also view as a data.frame
#' head(getAlignedFeatures(fe, format = "data.frame", proportional_threshold = 0.03))
#'
#' @rdname getAlignedFeatures
#' @name getAlignedFeatures
#' @export getAlignedFeatures
NULL

#' @rdname getAlignedFeatures
#' @export
setMethod("getAlignedFeatures", c("FactorisedExperiment"), function(
    object,
    loading_threshold = 0.5,
    proportional_threshold = 0.01,
    feature_id_col = "rownames",
    format = "list",
    center_loadings = FALSE
) {
    S <- loadings(
        object, scale_loadings = TRUE, center_loadings = center_loadings
    )

    if (feature_id_col != "rownames")
        rownames(S) <- rowData(object)[[feature_id_col]]

    abs_thresholds <- apply(S, 2, function(l) {
        max(abs(l)) * loading_threshold
    })
    perc_thresholds <- apply(S, 2, function(l) {
        stats::quantile(abs(l), probs = 1 - proportional_threshold)
    })

    factor_features <- data.frame()
    for (f in componentNames(object)) {
        abs_loadings <- abs(S[, f])

        which_features <- which(abs_loadings >= abs_thresholds[f] &
                                abs_loadings >= perc_thresholds[f])

        if (length(which_features) > 0) {
            factor_features <- rbind(factor_features, data.frame(
                component = f,
                feature = rownames(S)[which_features],
                value = S[, f][which_features],
                loadings_centered = center_loadings
            ))
        }
    }

    factor_features$loading_threshold <- loading_threshold
    factor_features$proportional_threshold <- proportional_threshold

    factor_features <- factor_features[
        order(abs(factor_features$value), decreasing = TRUE), ]
    factor_features <- factor_features[order(factor_features$component), ]

    if (is.function(format)) {
        return(format(factor_features))
    } else if (format == "data.frame") {
        return(factor_features)
    } else if (format == "list") {
        factor_list <- list()

        for (f in unique(factor_features$component)) {
            factor_list[[f]] <-
                factor_features$feature[which(factor_features$component == f)]
        }

        return(factor_list)
    }
})


#' Functional enrichment analyses for dimensionally-reduced data
#'
#' Method for applying pathway enrichment analysis to components identified
#' through dimensionality reduction (e.g., factors or modules).
#' Enrichment analyses are applied to each component
#' separately.
#'
#' @param object \link[ReducedExperiment]{FactorisedExperiment}  or
#' \link[ReducedExperiment]{ModularExperiment} object.
#'
#' @param method The method to use for identifying enriched pathways. One of
#' "overrepresentation" or "gsea". The "overrepresentation" method calls
#' \link[clusterProfiler]{enricher} whereas the "gsea" method calls
#' \link[clusterProfiler]{GSEA}. Note that "gsea" is not available for modules.
#'
#' @param feature_id_col The column in `rowData(object)` that will be used as a
#' feature ID. Setting this to `"rownames"` (default) instead uses
#' `rownames(object)`.
#'
#' @param as_dataframe If `TRUE`, the results will be returned as a data.frame.
#' Otherwise, the results will be returned as a list of objects created by
#' either \link[clusterProfiler]{enricher}, in the case of overrepresentation
#' analysis, or \link[clusterProfiler]{GSEA}, in the case of GSEA.
#'
#' @param center_loadings If `TRUE`, loadings will be centered
#' column-wise to have a mean of 0.
#'
#' @param abs_loadings If `TRUE`, the absolute values of the
#' loadings will be used for enrichment analysis. If `FALSE`, the signed
#' loadings will be used for GSEA enrichment. Note that, regardless of the
#' value of this term, the process used to select genes for overrepresentation
#' analysis will be based on absolute loadings.
#'
#' @param loading_threshold See
#' \link[ReducedExperiment]{getAlignedFeatures}. Only relevant for
#' overrepresentation analysis.
#'
#' @param proportional_threshold See
#' \link[ReducedExperiment]{getAlignedFeatures}. Only relevant for
#' overrepresentation analysis.
#'
#' @param ... Additional arguments to be passed to
#' \link[clusterProfiler]{GSEA} (if `method == "gsea"`) or
#' \link[clusterProfiler]{enricher} (if `method == "overrepresentation"`).
#'
#' @details
#' When running module analysis, the overrepresentation method identifies
#' pathways that are overrepresented in each module.
#'
#' For factor analysis, the overrepresentation method first identifies the genes
#' most highly aligned with each factor
#' (using \link[ReducedExperiment]{getAlignedFeatures}), then uses
#' the resulting gene lists to perform overrepresentation analysis. The GSEA
#' method instead uses the entire set of factor loadings, and identifies
#' pathways that are overrepresented in the tails of this distribution.
#'
#' @returns If `as_dataframe` is `TRUE`, the results will be returned as a
#' data.frame. Otherwise, the results will be returned as a list of objects
#' created by either \link[clusterProfiler]{enricher}, in the case of
#' overrepresentation analysis, or \link[clusterProfiler]{GSEA}, in the case of
#' GSEA.
#'
#' @seealso [ReducedExperiment::getMsigdbT2G()]
#'
#' @author Jack Gisby
#'
#' @examples
#' set.seed(2)
#' airway <- ReducedExperiment:::.getAirwayData(n_features = 2000)
#' airway_fe <- estimateFactors(
#'     airway,
#'     nc = 2,
#'     use_stability = FALSE,
#'     method = "imax"
#' )
#'
#' # Get pathways (e.g., by using ReducedExperiment::getMsigdbT2G())
#' t2g <- read.csv(system.file(
#'     "extdata",
#'     "msigdb_t2g_filtered.csv",
#'     package = "ReducedExperiment"
#' ))
#'
#' # Run overrepresentation analysis
#' overrep_res <- runEnrich(
#'     airway_fe,
#'     method = "overrepresentation",
#'     feature_id_col = "rownames",
#'     as_dataframe = TRUE,
#'     p_cutoff = 0.1,
#'     TERM2GENE = t2g,
#'     universe = rownames(airway_fe)
#' )
#'
#' head(overrep_res)
#'
#' @rdname enrichment
#' @name runEnrich
#' @export runEnrich
NULL

#' @rdname enrichment
#' @export
setMethod("runEnrich", c("FactorisedExperiment"), function(
    object,
    method = "gsea",
    feature_id_col = "rownames",
    center_loadings = FALSE,
    abs_loadings = FALSE,
    loading_threshold = 0.5,
    proportional_threshold = 0.01,
    as_dataframe = FALSE,
    ...
) {
    if (method == "gsea") {
        proportional_threshold <- loading_threshold <- NULL

        S <- loadings(
            object, scale_loadings = TRUE,
            center_loadings = center_loadings,
            abs_loadings = abs_loadings
        )

        if (feature_id_col != "rownames")
            rownames(S) <- rowData(object)[[feature_id_col]]

        enrich_res <- reducedGSEA(S, ...)

    } else if (method == "overrepresentation") {
        factor_features <- getAlignedFeatures(object,
            feature_id_col = feature_id_col,
            center_loadings = center_loadings,
            loading_threshold = loading_threshold,
            proportional_threshold = proportional_threshold
        )

        enrich_res <- reducedOA(factor_features, ...)
    } else {
        stop("Enrichment method not recognised")
    }

    for (comp in names(enrich_res)) {
        if (is.null(enrich_res[[comp]]@result)) next
        if (nrow(enrich_res[[comp]]@result) == 0) next

        enrich_res[[comp]]@result$loading_threshold <- loading_threshold
        enrich_res[[comp]]@result$proportional_threshold <-
            proportional_threshold
        enrich_res[[comp]]@result$loadings_centered <- center_loadings
        enrich_res[[comp]]@result$loadings_scaled <- TRUE
        enrich_res[[comp]]@result$abs_loadings <- abs_loadings
    }

    if (as_dataframe) {
        enrich_res <- lapply(enrich_res, function(object) {object@result})
        enrich_res <- do.call("rbind", enrich_res)
    }

    return(enrich_res)
})
