% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quality_control.R
\name{CalculateSeuratQC}
\alias{CalculateSeuratQC}
\title{Calculate some QC metrics}
\usage{
CalculateSeuratQC(
  object,
  assay = NULL,
  calculate.mapped.reads = FALSE,
  mol.info = NULL,
  total.reads = NULL,
  control.feature.list = NULL,
  use.names = FALSE,
  symbol.to.id = NULL,
  log10 = TRUE
)
}
\arguments{
\item{object}{Seurat object}

\item{assay}{Name of assay to use}

\item{calculate.mapped.reads}{Whether or not to calculate mapped reads metrics. Default is FALSE. If TRUE, will need to set molecule information}

\item{mol.info}{Molecule information for mapped reads calculation (molecule_info.h5 files of 10x)}

\item{total.reads}{a data.frame with column named 'total.reads', collecting total count of reads per barcode.This can be extracted from .BAM file that CellRanger output and is used to calculate mapping rate per cell.}

\item{control.feature.list}{A list containing one or more vectors (a character vector of ensembl IDs named by feature symbols), used to identify feature controls}

\item{use.names}{Whether or not to use feature names instead of ID. If FALSE, shall specify symbol.to.id data.frame for transferation}

\item{symbol.to.id}{Data.frame for convertion from gene symbols to such as ENSEMBL ID, must contain at least two columns (use 'ID' and 'Symbol' as colnames)}

\item{log10}{Whether or not to calculate log10 transformed QC metrics}
}
\value{
Returns Seurat object with new meta.data storing QC metrics
}
\description{
Compute quality control (QC) metrics for each control feature and cell in a Seurat object, accounting for specified control sets.
}
