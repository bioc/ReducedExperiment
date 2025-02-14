% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/features.R
\name{getMsigdbT2G}
\alias{getMsigdbT2G}
\title{Get TERM2GENE dataframe from MSigDB}
\usage{
getMsigdbT2G(
  species = "Homo sapiens",
  category = "C2",
  subcategory = NULL,
  subcategory_to_remove = "CGP",
  gene_id = "ensembl_gene"
)
}
\arguments{
\item{species}{The species for which to obtain MSigDB pathways. See
\link[msigdbr]{msigdbr} for more details.}

\item{category}{The MSigDB category to retrieve pathways for. See
\link[msigdbr]{msigdbr} for more details.}

\item{subcategory}{The MSigDB subcategory to retrieve pathways for. See
\link[msigdbr]{msigdbr} for more details.}

\item{subcategory_to_remove}{If not NULL, this is a character string
indicating a subcategory to be removed from the results of
\link[msigdbr]{msigdbr}.}

\item{gene_id}{The name to be given to the gene_id column of the resulting
data.frame.}
}
\value{
Returns a data.frame, where the \code{gs_name} column indicates the name
of a pathway, and the \code{gene_id} column indicates genes that belong to
said pathway.
}
\description{
Gets pathways from the MSigDB database in the format required by
\code{clusterProfiler} enrichment functions, such as
\link[clusterProfiler]{enricher} and \link[clusterProfiler]{GSEA}.
May be used as input to \link[ReducedExperiment]{runEnrich}. By default,
retrieves the C2 canonical pathways.
}
\examples{
pathways <- getMsigdbT2G(
    species = "Homo sapiens",
    category = "C2",
    subcategory_to_remove = "CGP",
    gene_id = "ensembl_gene"
)

# A data.frame indicating gene-pathway mappings for use in pathway analysis
head(pathways)

}
\author{
Jack Gisby
}
