% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/delta.R
\name{delta_score}
\alias{delta_score}
\title{Delta Scores}
\usage{
delta_score(data, taxa, method = "hamming")
}
\arguments{
\item{data}{The data matrix used to calculate distance between taxa. Each taxa should be a column in the data.}

\item{taxa}{A vector of the taxa analysed. Should relate to the columns in data}
}
\value{
The delta score for the entire dataset, and the average delta score for each taxon.
}
\description{
Delta Scores
}
