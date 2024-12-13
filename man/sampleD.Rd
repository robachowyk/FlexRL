% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{sampleD}
\alias{sampleD}
\title{sampleD}
\usage{
sampleD(S, LLA, LLB, LLL, gamma, loglik, nlinkrec, sumRowD, sumColD)
}
\arguments{
\item{S}{IntegerMatrix where each row correspond to the indices (from source A and source B) of records for which the true values matches (representing the potential links)}

\item{LLA}{NumericVector gives the likelihood contribution of each non linked record from A}

\item{LLB}{NumericVector gives the likelihood contribution of each non linked record from B}

\item{LLL}{NumericVector gives the likelihood contribution of each potential linked records (from select)}

\item{gamma}{NumericVector repeats the value of the parameter gamma (proportion of linked records) number of potential linked records (nrow of S) times}

\item{loglik}{double for the value of the current complete log likelihood of the model}

\item{nlinkrec}{integer for the current number of linked records}

\item{sumRowD}{A LogicalVector vector indicating, for each row of the linkage matrix, i.e. for each record in the smallest file A, whether the record has a link in B or not.}

\item{sumColD}{A LogicalVector vector indicating, for each column of the linkage matrix, i.e. for each record in the largest file B, whether the record has a link in A or not.}
}
\value{
List:
\itemize{
\item new set of links
\item new sumRowD
\item new sumColD
\item new value of the complete log likelihood
\item new number fo linked records
}
}
\description{
sampleD
}
