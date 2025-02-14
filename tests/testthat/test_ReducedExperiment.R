context("ReducedExperiment")

test_that("Build and subset", {
    i <- 300
    j <- 100
    k <- 10

    set.seed(1)
    rrs <- ReducedExperiment:::.createRandomisedReducedExperiment(i = i, j = j, k = k)

    # Test reduced experiment
    expect_equal(dim(rrs), c("Features" = i, "Samples" = j, "Components" = k))
    expect_equal(dim(rrs), c(nFeatures(rrs), nSamples(rrs), nComponents(rrs)))

    expect_equal(colnames(assay(rrs, "normal")), sampleNames(rrs))
    expect_equal(rownames(reduced(rrs)), sampleNames(rrs))
    expect_equal(rownames(rowData(rrs)), rownames(rrs))
    show(rrs)

    rrs@scale <- setNames(1:i, featureNames(rrs))
    rrs@center <- setNames(1:i, featureNames(rrs))

    # Subset and re-test
    rrs_subset <- rrs[5:10, 50:90, 1:2]
    expect_equal(dim(rrs_subset), c("Features" = 6, "Samples" = 41, "Components" = 2))
    expect_equal(nComponents(rrs_subset), c("Components" = 2))
    expect_equal(nSamples(rrs_subset), c("Samples" = 41))
    expect_equal(nFeatures(rrs_subset), c("Features" = 6))
    expect_true(ncol(rrs_subset) == nSamples(rrs_subset))
    expect_true(nrow(rrs_subset) == nFeatures(rrs_subset))
    expect_equal(rownames(reduced(rrs_subset)), sampleNames(rrs_subset))
    expect_equal(paste0("sample_", 50:90), sampleNames(rrs_subset))
    expect_equal(names(rrs_subset@scale), featureNames(rrs_subset))
    expect_equal(names(rrs_subset@center), featureNames(rrs_subset))

    names(rrs_subset) <- paste0("123_", rownames(rrs_subset))
    expect_equal(rownames(rrs_subset)[1], "123_gene_5")
    expect_true(validObject(rrs_subset))

    # Now test an empty object
    rrs_empy <- ReducedExperiment()
    expect_equal(dim(rrs_empy), c("Features" = 0, "Samples" = 0, "Components" = 0))
    expect_equal(reduced(rrs_empy), matrix(0, 0, 0), check.attributes = FALSE)
    expect_equal(rrs_empy@scale, TRUE)
    expect_equal(rrs_empy@center, TRUE)
    expect_true(validObject(rrs_subset))

    # Subset with characters
    expect_equal(
        dimnames(rrs[5:10, 50:90, 1:2]),
        dimnames(rrs[paste0("gene_", 5:10), paste0("sample_", 50:90), paste0("factor_", 1:2)])
    )

    # Test subset replacement
    rrs[5:10, 50:90, 1:2] <- rrs[paste0("gene_", 5:10), paste0("sample_", 50:90), paste0("factor_", 1:2)]
    expect_true(validObject(rrs))

    rrs[paste0("gene_", 5:10), paste0("sample_", 50:90), paste0("factor_", 1:2)] <- rrs[5:10, 50:90, 1:2]
    expect_true(validObject(rrs))

    # Check changes to subsetted object are reflected after allocation
    rrs_subset <- rrs[5:10, 50:90, 1:2]
    rownames(rrs_subset)[3] <- "renamed_feature"
    colnames(rrs_subset)[2] <- "renamed_sample"
    componentNames(rrs_subset)[1] <- "renamed_comp"
    rownames(rrs)[7] <- "renamed_feature"
    colnames(rrs)[51] <- "renamed_sample"
    componentNames(rrs)[1] <- "renamed_comp"

    reduced(rrs_subset)["renamed_sample", "renamed_comp"] <- 10
    rrs_subset@scale["renamed_feature"] <- 8

    rrs[5:10, 50:90, 1:2] <- rrs_subset
    expect_true(validObject(rrs))
    expect_true(validObject(rrs_subset))
    expect_equal(reduced(rrs)["renamed_sample", "renamed_comp"], 10)
    expect_true(rrs_subset@scale["renamed_feature"] == 8)
})

test_that("Access and replace reduced data", {
    i <- 300
    j <- 100
    k <- 10

    set.seed(1)
    reduced_data <- .makeRandomData(j, k, "sample", "factor")

    # Make a ReducedExperiment and save original reduced data to an object
    set.seed(1)
    rrs <- ReducedExperiment(
        assays = list("normal" = .makeRandomData(i, j, "gene", "sample")),
        reduced = reduced_data
    )

    # Expect that the reduced method returns the original data
    expect_equal(reduced(rrs, scale_reduced = TRUE, center_reduced = TRUE), scale(reduced_data))
    expect_equal(reduced(rrs), reduced_data)

    # Replacing a value should work
    reduced(rrs)[3, 5] <- 5
    expect_equal(reduced(rrs)[3, 5], 5)

    reduced(rrs)[3, 5] <- reduced_data[3, 5]
    expect_equal(reduced(rrs), reduced_data, check.attributes = FALSE)

    # This should work
    reduced(rrs) <- reduced_data[, 1:5]

    # This should not work (validity should fail because different number of samples)
    expect_error((function() {
        reduced(rrs) <- reduced_data[1:5, ]
    })())
})

test_that("Access and replace component names", {
    set.seed(1)
    rrs <- .createRandomisedReducedExperiment(i = 300, j = 100, k = 10)

    expect_equal(componentNames(rrs), paste0("factor_", 1:10))
    expect_equal(colnames(reduced(rrs)), paste0("factor_", 1:10))

    componentNames(rrs)[5] <- "new_name"
    expect_equal(componentNames(rrs)[5], "new_name")
    expect_equal(colnames(reduced(rrs))[5], "new_name")
})

test_that("Access and replace sample names", {
    set.seed(1)
    rrs <- ReducedExperiment:::.createRandomisedReducedExperiment(i = 300, j = 100, k = 10)

    expect_equal(sampleNames(rrs), paste0("sample_", 1:100))
    expect_equal(rownames(colData(rrs)), paste0("sample_", 1:100))
    expect_equal(colnames(assay(rrs, 1)), paste0("sample_", 1:100))
    expect_equal(rownames(reduced(rrs)), paste0("sample_", 1:100))

    sampleNames(rrs)[5] <- "new_name"
    expect_equal(sampleNames(rrs)[5], "new_name")
    expect_equal(rownames(colData(rrs))[5], "new_name")
    expect_equal(colnames(assay(rrs, 1))[5], "new_name")
    expect_equal(rownames(reduced(rrs))[5], "new_name")
})

test_that("Access and replace feature names", {
    set.seed(1)
    rrs <- ReducedExperiment:::.createRandomisedReducedExperiment(i = 300, j = 100, k = 10)

    rrs@scale <- setNames(1:300, featureNames(rrs))
    rrs@center <- setNames(1:300, featureNames(rrs))

    expect_equal(featureNames(rrs), paste0("gene_", 1:300))
    expect_equal(rownames(assay(rrs, 1)), paste0("gene_", 1:300))
    expect_equal(names(rrs@center), paste0("gene_", 1:300))
    expect_equal(names(rrs@scale), paste0("gene_", 1:300))

    featureNames(rrs)[5] <- "new_name"
    expect_equal(featureNames(rrs)[5], "new_name")
    expect_equal(rownames(assay(rrs, 1))[5], "new_name")
    expect_equal(names(rrs@center)[5], "new_name")
    expect_equal(names(rrs@scale)[5], "new_name")
})

test_that("Access and replace scale/center", {
    set.seed(1)
    rrs <- .createRandomisedReducedExperiment(i = 300, j = 100, k = 10)

    # This should work
    rrs@scale <- setNames(1:300, featureNames(rrs))
    rrs@center <- setNames(1:300, featureNames(rrs))

    # This should not work (mismatch with feature names/length)
    expect_error((function() {
        rrs@scale <- setNames(1:10, paste0("module_", 1:10))
        validObject(rrs)
    })())
    expect_error((function() {
        rrs@scale <- setNames(1:5, paste0("factor_", 1:5))
        validObject(rrs)
    })())
    expect_error((function() {
        rrs@scale <- 1:5
        validObject(rrs)
    })())
    expect_error((function() {
        rrs@center <- setNames(1:10, paste0("module_", 1:10))
        validObject(rrs)
    })())
    expect_error((function() {
        rrs@center <- setNames(1:5, paste0("factor_", 1:5))
        validObject(rrs)
    })())
    expect_error((function() {
        rrs@center <- 1:5
        validObject(rrs)
    })())
})

test_that("Combine ReducedExperiments with cbind and rbind", {
    set.seed(1)
    rrs_a <- ReducedExperiment:::.createRandomisedReducedExperiment(i = 300, j = 100, k = 10)

    set.seed(2)
    rrs_b <- ReducedExperiment:::.createRandomisedReducedExperiment(i = 300, j = 100, k = 10)

    # Objects should be cbind-able due to matching names
    rrs_a_a <- cbind(rrs_a, rrs_a)
    expect_true(validObject(rrs_a_a))
    rrs_a_b <- cbind(rrs_a, rrs_b)
    expect_true(validObject(rrs_a_b))

    expect_equal(dim(rrs_a_b), c("Features" = 300, "Samples" = 200, "Components" = 10))

    # Not rbindable when non-matching reduced
    rrs_a_a <- rbind(rrs_a, rrs_a)
    expect_true(validObject(rrs_a_a))
    expect_error(rbind(rrs_a, rrs_b))

    reduced(rrs_b) <- reduced(rrs_a)
    rrs_a_b <- rbind(rrs_a, rrs_b)
    expect_true(validObject(rrs_a_b))

    expect_equal(dim(rrs_a_b), c("Features" = 600, "Samples" = 100, "Components" = 10))

    # Add scaling information to rrs_b but not a
    rrs_b@scale <- 1:300
    names(rrs_b@scale) <- rownames(rrs_b)
    expect_true(validObject(rrs_b))

    # Both should fail due to non-matching scaling slot types
    expect_error(cbind(rrs_a, rrs_b))
    expect_error(rbind(rrs_a, rrs_b))

    # Add scaling information to rrs_a and it should work again
    rrs_a@scale <- rrs_b@scale
    expect_no_error(cbind(rrs_a, rrs_b))
    expect_no_error(rbind(rrs_a, rrs_b))

    # rbind works with non-matching information if matching types
    rrs_b@scale <- 1:300 * 2
    names(rrs_b@scale) <- rownames(rrs_b)
    expect_true(validObject(rrs_b))
    expect_error(cbind(rrs_a, rrs_b))
    expect_no_error(rbind(rrs_a, rrs_b))

    rrs_a@scale <- rrs_b@scale <- TRUE

    # Same with centering
    rrs_b@center <- 1:300
    names(rrs_b@center) <- rownames(rrs_b)
    expect_true(validObject(rrs_b))
    expect_error(cbind(rrs_a, rrs_b))
    expect_error(rbind(rrs_a, rrs_b))

    rrs_a@center <- rrs_b@center
    expect_no_error(cbind(rrs_a, rrs_b))
    expect_no_error(rbind(rrs_a, rrs_b))

    rrs_b@center <- 1:300 * 2
    names(rrs_b@center) <- rownames(rrs_b)
    expect_true(validObject(rrs_b))
    expect_error(cbind(rrs_a, rrs_b))
    expect_no_error(rbind(rrs_a, rrs_b))

    # Both should fail with non-matching comp-names
    componentNames(rrs_b)[5] <- "new_name"
    expect_error(cbind(rrs_a, rrs_b))
    expect_error(rbind(rrs_a, rrs_b))
})
