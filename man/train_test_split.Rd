% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/train-test-split.R
\name{train_test_split}
\alias{train_test_split}
\title{Splitting the data sets to train and test}
\usage{
train_test_split(train, train.ratio = 0.7, trim = FALSE, trim.ratio = NULL)
}
\arguments{
\item{train}{the train data of class \code{"hhsmmdata"}}

\item{train.ratio}{a number in (0,1] which determines the ratio of the train subset. It can be
equal to 1, if we need the test set to be equal to the train set and we only
need to right trim the sequences}

\item{trim}{logical. if TRUE the sequences will be right trimmed with random lengths}

\item{trim.ratio}{a vector of trim ratios with a length equal to that of \code{train$N},
or a single trim ratio for all sequences. If it is \code{NULL}, then random trim
ratios will be used}
}
\value{
a list containing:
\itemize{
\item\code{train} the randomly selected subset of train data of class \code{"hhsmmdata"}
\item\code{test} the randomly selected subset of test data of class \code{"hhsmmdata"}
\item\code{trimmed} right trimmed test subset, if \code{trim}=TRUE, with trim ratios equal to \code{trim.ratio}
\item\code{trimmed.count} the number of right trimmed individuals in each sequence of the test subset, if \code{trim}=TRUE
}
}
\description{
A function to split the train data of class \code{"hhsmmdata"}
to train and test subsets with an option to right trim the sequences
}
\details{
This function splits the sample to train and test samples and
trims the test sample from right, in order to provide a sample for examination
of the prediction tools.
In reliability applications, the hhsmm models are often left-to-right
and the modeling aims to predict the future states. In such cases, the
test sets are right trimmed and the prediction aims to predict the
residual useful lifetime (RUL) of a new sequence.
}
\examples{
\donttest{
data(CMAPSS)
tt = train_test_split(CMAPSS$train, train.ratio = 0.7, trim = TRUE)
}

}
\author{
Morteza Amini, \email{morteza.amini@ut.ac.ir}
}
