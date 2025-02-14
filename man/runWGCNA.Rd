% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modules.R
\name{runWGCNA}
\alias{runWGCNA}
\title{Run WGCNA for a data matrix}
\usage{
runWGCNA(
  X,
  power,
  cor_type = "pearson",
  networkType = "signed",
  module_labels = "numbers",
  maxBlockSize = 30000,
  verbose = 0,
  standardise_reduced = TRUE,
  ...
)
}
\arguments{
\item{X}{Either a \link[SummarizedExperiment]{SummarizedExperiment} object
or a matrix containing data to be subject to WGCNA. \code{X} should have rows as
features and columns as samples.}

\item{power}{ soft-thresholding power for network construction. }

\item{cor_type}{The type of correlation to be used to generate a correlation
matrix during network formation. One of "pearson" (\link[WGCNA]{cor}) and
"bicor" (\link[WGCNA]{bicor}).}

\item{networkType}{ network type. Allowed values are (unique abbreviations of) \code{"unsigned"},
\code{"signed"}, \code{"signed hybrid"}. See \code{\link[WGCNA]{adjacency}}. }

\item{module_labels}{Specifies whether the modules should be named based on
"numbers" or "colours. If \code{module_labels} is set to "numbers", then
"module_0" represents unclustered genes, whereas if it is set to "colours"
then "grey" represents unclustered genes.}

\item{maxBlockSize}{The chunk size (in terms of the number of features/genes)
to process the data. See \link[WGCNA]{blockwiseModules} for
more details. The default (30000) should process standard transcriptomic
datasets in a single chunk. Results may differ if the number of features
exceeds the chunk size. Lower values of this parameter may use less memory
to calculate networks.}

\item{verbose}{ integer level of verbosity. Zero means silent, higher values make the output
progressively more and more verbose. }

\item{standardise_reduced}{If \code{TRUE}, the reduced data (eigengenes) are
standardised to have a mean of 0 and a standard deviation of 1.}

\item{...}{Additional arguments to be passed to
\link[WGCNA]{blockwiseModules}.}
}
\value{
Returns a list containing:
\describe{
\item{"E"}{The reduced data (eigengenes).}
\item{"L"}{The module loadings. This represents the values of the PCA
rotation matrix for the first principal component generated for each
module.}
\item{"assignments"}{A named vector representing the assignments of
genes to modules.}
}
}
\description{
Runs WGCNA. Largely a wrapper for the \link[WGCNA]{blockwiseModules}
function that reformats data
into a format convenient for creating a
\link[ReducedExperiment]{ModularExperiment} object and changes module names
from colours to numbers by default.
}
\details{
Note that if \code{module_labels} is set to "numbers", then
"module_0" represents unclustered genes, whereas if it is set to "colours"
then "grey" represents unclustered genes.

The function also stores the loadings matrices generated when PCA is
performed for each module to calculate eigengenes. These loadings can be
used to recalculate the reduced data matrix (eigengenes).
}
\examples{
# Get the airway data as a SummarizedExperiment (with a subset of features)
set.seed(2)
airway_se <- ReducedExperiment:::.getAirwayData(n_features = 500)

# Choose an appropriate soft-thresholding power
WGCNA::disableWGCNAThreads()
fit_indices <- assessSoftThreshold(airway_se)
estimated_power <- fit_indices$Power[fit_indices$estimated_power]

# Identify modules using the airway expression matrix
wgcna_res <- runWGCNA(
    assay(airway_se, "normal"),
    verbose = 0,
    power = estimated_power
)

# We find just one module for this small dataset (module_0 indicates unclustered genes)
table(names(wgcna_res$assignments))

}
\seealso{
\code{\link[WGCNA:blockwiseModules]{WGCNA::blockwiseModules()}},
\code{\link[=assessSoftThreshold]{assessSoftThreshold()}},
\code{\link[WGCNA:pickSoftThreshold]{WGCNA::pickSoftThreshold()}},
}
\author{
Jack Gisby
}
