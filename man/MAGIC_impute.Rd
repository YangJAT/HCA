% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anno_function.R
\name{MAGIC_impute}
\alias{MAGIC_impute}
\title{MAGIC_impute}
\usage{
MAGIC_impute(
  datafilt,
  knn = 30,
  t = 3,
  npca = 50,
  only_marker = TRUE,
  organism = "human",
  ncore = 10
)
}
\arguments{
\item{datafilt}{your Seurat object(support V4)}

\item{only_marker}{impute marker genes}

\item{organism}{human or mouse}

\item{ncore}{number of cores}
}
\value{
Seurat object
}
\description{
MAGIC_impute
}
