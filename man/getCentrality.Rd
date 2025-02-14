% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods-ModularExperiment.R
\name{getCentrality}
\alias{getCentrality}
\alias{getCentrality,ModularExperiment-method}
\title{Get correlation of features with module eigengenes}
\usage{
\S4method{getCentrality}{ModularExperiment}(object, assay_name = "normal", feature_id_col = "rownames")
}
\arguments{
\item{object}{A \link[ReducedExperiment]{ModularExperiment} object.}

\item{assay_name}{The name of the assay to be used for calculation of
module centrality.}

\item{feature_id_col}{The column in \code{rowData(object)} that will be used as a
feature ID. Setting this to "rownames" (default) instead uses
\code{rownames(object)}.}
}
\value{
Returns a \code{data.frame} with columns for \code{feature},
\code{r} (signed correlation with the eigengene),
\code{rsq} (squared correlation with the eigengene),
\code{rank_r} (feature rank based on \code{r}) and
\code{rank_rsq} (feature rank based on rsq).
}
\description{
Provides a wrapper around \link[WGCNA]{signedKME}. Provides a measure
of module centrality/connectivity of each feature. Calculates
correlation (Pearson's r) of each feature with the module
eigengene (i.e., the column of \code{reduced} to which the feature belongs).
}
\examples{
# Create ModularExperiment with random data (100 features, 50 samples,
# 10 modules)
me <- ReducedExperiment:::.createRandomisedModularExperiment(100, 50, 10)
me

# Calculate centrality of each feature for the corresponding module
head(getCentrality(me))

}
\seealso{
\code{\link[WGCNA:signedKME]{WGCNA::signedKME()}}
}
\author{
Jack Gisby
}
