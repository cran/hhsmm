% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/homogeneity.R
\name{homogeneity}
\alias{homogeneity}
\title{Computing maximum homogeneity of two state sequences}
\usage{
homogeneity(state.seq1, state.seq2)
}
\arguments{
\item{state.seq1}{first state sequence}

\item{state.seq2}{second state sequence}
}
\value{
a vector of a length equal to the maximum number of states giving
the maximum homogeneity ratios
}
\description{
A function to compute the maximum homogeneity of two state sequences.
}
\examples{
state.seq1 = c(3, 3, 3, 1, 1, 2, 2, 2, 2)
state.seq2 = c(2, 2, 2, 3, 3, 1, 1, 1, 1)
homogeneity(state.seq1, state.seq2)

}
\author{
Morteza Amini, \email{morteza.amini@ut.ac.ir}
}
