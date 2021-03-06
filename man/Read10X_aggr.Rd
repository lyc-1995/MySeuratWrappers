% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{Read10X_aggr}
\alias{Read10X_aggr}
\title{Load in data from 10X Cellranger aggr}
\usage{
Read10X_aggr(
  data.dir,
  sample.names,
  split = FALSE,
  gene.column = 2,
  unique.features = TRUE
)
}
\arguments{
\item{data.dir}{Directory containing the matrix.mtx, genes.tsv (or features.tsv), and barcodes.tsv
files provided by by 10X Cellranger aggr}

\item{sample.names}{A vector of characters indicating the sample names. The length of sample.names must be equal
to the number of libraries input to aggr pipeline}

\item{split}{Whether or not to split the matrix by each sample}

\item{gene.column}{Specify which column of genes.tsv or features.tsv to use for gene names; default is 2}

\item{unique.features}{Make feature names unique (default TRUE)}
}
\value{
If split=FALSE, will return a sparse matrix with cell barcodes attched to each sample name. If split=TRUE,
will return a list of splited matrices
}
\description{
Enables easy loading of sparse data matrices provided by 10X Cellranger aggr.
}
