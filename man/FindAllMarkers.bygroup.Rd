% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{FindAllMarkers.bygroup}
\alias{FindAllMarkers.bygroup}
\title{Gene expression markers shared by multiple groups of cells}
\usage{
FindAllMarkers.bygroup(
  object,
  split.by = NULL,
  only.pos = TRUE,
  features = NULL,
  min.cells.group = 10,
  min.freq = 0.5,
  ...
)
}
\arguments{
\item{object}{A Seurat object}

\item{split.by}{A metadata column name - the data will be split by this column to calculate \link[Seurat]{FindAllMarkers}
separately for each data split}

\item{only.pos}{Only return positive markers (TRUE by default)}

\item{features}{Genes to test. Default is to use all genes}

\item{min.cells.group}{Minimum number of cells in the group - if lower the group is skipped}

\item{min.freq}{Only return markers which are differentially expressed in at least this fraction of datasets.}

\item{...}{Additional paramters to \link[Seurat]{FindAllMarkers}}
}
\value{
A list of marker genes for each identity class (typically clusters), with two associated numerical values:
    i) the fraction of datasets for which the marker was found to be differentially expressed; ii) the
    average log-fold change for the genes across datasets
}
\description{
This function expands \link[Seurat]{FindAllMarkers} to find markers that are differentially expressed across multiple
datasets or samples. Given a Seurat object with identity classes (for example annotated clusters) and a grouping
variable (for example a Sample ID), it calculate differentially expressed genes (DEGs) individually for each sample.
Then it determines the fraction of samples for which the gene was found to be differentially expressed.
}
\details{
This function can be useful to find marker genes that are specific for individual cell types, and that are found
to be so consistently across multiple samples.
}
\examples{
\dontshow{if (interactive()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
library(Seurat)
ref <- load.reference.map(ref = "https://figshare.com/ndownloader/files/38921366")
Idents(ref) <- "functional.cluster"
FindAllMarkers.bygroup(ref, split.by = "Sample", min.cells.group=30, min.freq=0.8)
\dontshow{\}) # examplesIf}
}
