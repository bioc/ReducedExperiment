
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ReducedExperiment

<!-- <img src="inst/ReducedExperiment_hex.png" align="right" height="174" width="150" /> -->
<!-- badges: start -->

[![build-test](https://github.com/jackgisby/ReducedExperiment/actions/workflows/build-test-deploy.yml/badge.svg)](https://github.com/jackgisby/ReducedExperiment/actions/workflows/build-test-deploy.yml)
[![codecov](https://codecov.io/gh/jackgisby/ReducedExperiment/graph/badge.svg?token=FHNH7AA6S3)](https://codecov.io/gh/jackgisby/ReducedExperiment)
[![GitHub
issues](https://img.shields.io/github/issues/jackgisby/ReducedExperiment)](https://github.com/jackgisby/ReducedExperiment/issues)
[![GitHub
pulls](https://img.shields.io/github/issues-pr/jackgisby/ReducedExperiment)](https://github.com/jackgisby/ReducedExperiment/pulls)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- [![check-bioc](https://github.com/jackgisby/ReducedExperiment/actions/workflows/check-bioc.yml/badge.svg)](https://github.com/jackgisby/ReducedExperiment/actions/workflows/check-bioc.yml) -->
<!-- [![Bioc release status](http://www.bioconductor.org/shields/build/release/bioc/ReducedExperiment.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/ReducedExperiment) -->
<!-- [![Bioc devel status](http://www.bioconductor.org/shields/build/devel/bioc/ReducedExperiment.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/ReducedExperiment) -->
<!-- [![Bioc downloads rank](https://bioconductor.org/shields/downloads/release/ReducedExperiment.svg)](http://bioconductor.org/packages/stats/bioc/ReducedExperiment/) -->
<!-- [![Bioc support](https://bioconductor.org/shields/posts/ReducedExperiment.svg)](https://support.bioconductor.org/tag/ReducedExperiment) -->
<!-- [![Bioc history](https://bioconductor.org/shields/years-in-bioc/ReducedExperiment.svg)](https://bioconductor.org/packages/release/bioc/html/ReducedExperiment.html#since) -->
<!-- [![Bioc last commit](https://bioconductor.org/shields/lastcommit/devel/bioc/ReducedExperiment.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/ReducedExperiment/) -->
<!-- [![Bioc dependencies](https://bioconductor.org/shields/dependencies/release/ReducedExperiment.svg)](https://bioconductor.org/packages/release/bioc/html/ReducedExperiment.html#since) -->

<!-- badges: end -->

ReducedExperiment provides containers for storing and manipulating
dimensionally-reduced assay data. The ReducedExperiment classes allow
users to simultaneously manipulate their original dataset and their
decomposed data, in addition to other method-specific outputs like
pathway analysis. Implements utilities and specialised classes for the
application of stabilised independent component analysis (sICA) and
weighted gene correlation network analysis (WGCNA).

## Installation

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then install `ReducedExperiment`
from
[Bioconductor](https://bioconductor.org/packages/release/bioc/html/ReducedExperiment.html)
using the following code:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("ReducedExperiment")
```

Alternatively, the development version of `ReducedExperiment` is
available from
[Bioconductor](https://bioconductor.org/packages/devel/bioc/html/ReducedExperiment.html)
or [GitHub](https://github.com/jackgisby/ReducedExperiment) with:

``` r
BiocManager::install("ReducedExperiment", version = "devel")

devtools::install_github("jackgisby/ReducedExperiment")
```

The development version of the package is also available as a container
on
[DockerHub](https://hub.docker.com/repository/docker/jackgisby/reducedexperiment/).

## Usage

`ReducedExperiment` objects are derived from `SummarizedExperiment`
objects, with additional slots designed to store and manipulate the
outputs of common dimensionality reduction techniques.

As an example, the `SummarizedExperiment` described below contains gene
expression data from individuals with COVID-19. It contains the
following slots:

- `assays` - A features by samples matrix containing the expression
  data.
- `colData` - Contains a row for each sample containing phenotype data.
- `rowData` - Contains a row for each feature containing gene IDs.

The `SummarizedExperiment` objects are convenient because, when we slice
the rows or columns of the expression matrix, the metadata for the rows
and columns are sliced accordingly.

``` r
library("SummarizedExperiment")

se <- readRDS(system.file(
    "extdata",
    "wave1.rds",
    package = "ReducedExperiment"
))

se
#> class: SummarizedExperiment 
#> dim: 500 83 
#> metadata(0):
#> assays(1): normal
#> rownames(500): ENSG00000004799 ENSG00000007038 ... ENSG00000287935
#>   ENSG00000288049
#> rowData names(2): ensembl_id gene_id
#> colnames(83): C37_positive_9 C48_positive_4 ... C85_negative
#>   C89_negative
#> colData names(8): sample_id individual_id ... case_control
#>   time_from_first_x
```

The `SummarizedExperiment` has two dimensions, representing the features
(2,184) and samples (234).

We can perform a factor analysis on these data, the result of which is a
set of reduced components and feature loadings.

``` r
library("ReducedExperiment")

fe <- estimate_factors(se, nc = 35)
fe
#> class: FactorisedExperiment 
#> dim: 500 83 35 
#> metadata(0):
#> assays(2): normal transformed
#> rownames(500): ENSG00000004799 ENSG00000007038 ... ENSG00000287935
#>   ENSG00000288049
#> rowData names(2): ensembl_id gene_id
#> colnames(83): C37_positive_9 C48_positive_4 ... C85_negative
#>   C89_negative
#> colData names(8): sample_id individual_id ... case_control
#>   time_from_first_x
#> 35 components
```

This `FactorisedExperiment` object has an additional dimension
representing the 35 factors. It also has additional slots, including:

- `reduced` - A samples by factors matrix containing the
  dimensionally-reduced data.
- `loadings` - Contains a features by factors matrix containing the
  loadings.

The `ReducedExperiment` objects allow users to simultaneously slice and
modify the `assays`, `rowData`, `colData`, `reduced` and `loadings`
matrices. Here, we provided a `SummarizedExperiment` object to
`estimate_factors`, but we could just have easily provided a simple
expression matrix.

Alternatively, you may have already applied dimensionality reduction to
your data and simply wish to package it into a `ReducedExperiment`
container. For instance, below we apply principal components analysis,
and construct a `FactorisedExperiment` object from the results.

``` r
prcomp_res <- stats::prcomp(t(assay(se)), center = TRUE, scale. = TRUE)

fe_prcomp <- FactorisedExperiment(
    se,
    reduced = prcomp_res$x,
    loadings = prcomp_res$rotation,
    stability = prcomp_res$sdev,
    center = prcomp_res$center,
    scale = prcomp_res$scale
)

fe_prcomp
#> class: FactorisedExperiment 
#> dim: 500 83 83 
#> metadata(0):
#> assays(1): ''
#> rownames(500): ENSG00000004799 ENSG00000007038 ... ENSG00000287935
#>   ENSG00000288049
#> rowData names(0):
#> colnames(83): C37_positive_9 C48_positive_4 ... C85_negative
#>   C89_negative
#> colData names(0):
#> 83 components
```

## Functionality

The package currently provides three types of container:

- `ReducedExperiment` - A basic container that can store
  dimensionally-reduced components.
- `FactorisedExperiment` - A container based on `ReducedExperiment`
  designed for working with the results of factor analysis. It can
  contain feature loadings and factor stability values.
- `ModularExperiment` - A container based on `ReducedExperiment`
  designed for working with modules of features (usually genes), as is
  produced by the popular Weighted Gene Correlation Network Analysis
  (WGCNA) approach. It contains the mapping of features to modules

Various tools are provided by the package for applying dimensionality
reduction and manipulating their results. These include:

- Workflows for applying independent component analysis (ICA) and WGCNA.
  We additionally developed an R implementation of the stabilised ICA
  algorithm.
- Methods for applying pathway enrichment analysis to factors and
  modules.
- Functions for identifying associations between factors/modules and
  sample-level variables.
- Methods for applying identified factors or modules to new datasets.
- Other approach-specific plots and utilities, such as factor stability
  plots and module preservation plots.

Many of these are demonstrated in more detail in the [package’s
vignette](https://jackgisby.github.io/ReducedExperiment/articles/ReducedExperiment.html).

The containers implemented in `ReducedExperiment` are designed to be
extensible. We encourage the development of children classes with
additional, or alternative, slots and methods.

## Citation

Below is the citation output from using `citation('ReducedExperiment')`
in R.

``` r
print(citation("ReducedExperiment"), bibtex = TRUE)
#> To cite package 'ReducedExperiment' in publications use:
#> 
#>   Gisby JS, Barnes MR (2025). _ReducedExperiment: Containers and tools
#>   for dimensionally-reduced -omics data_.
#>   doi:10.18129/B9.bioc.ReducedExperiment
#>   <https://doi.org/10.18129/B9.bioc.ReducedExperiment>, v0.99.3,
#>   <http://www.bioconductor.org/packages/ReducedExperiment>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {ReducedExperiment: Containers and tools for dimensionally-reduced -omics data},
#>     author = {Jack S. Gisby and Michael R. Barnes},
#>     year = {2025},
#>     url = {http://www.bioconductor.org/packages/ReducedExperiment},
#>     note = {v0.99.3},
#>     doi = {10.18129/B9.bioc.ReducedExperiment},
#>   }
```

The `ReducedExperiment` package relies on many software packages and
development tools. The packages used are listed in the vignette and
relevant papers are cited.
