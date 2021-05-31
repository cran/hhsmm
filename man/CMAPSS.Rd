\name{CMAPSS}
\alias{CMAPSS}
\docType{data}
\title{
CMAPSS data set
}
\description{
Commercial Modular Aero-Propulsion System Simulation (C-MAPSS) data set provided by NASA.
}
\usage{data("CMAPSS")}
\format{
  A list of the following 2 objects:
  \describe{
    \item{\code{train}}{a list of class "hhsmm.data" as the train dataset}
    \item{\code{test}}{a list of class "hhsmm.data" as the test dataset}
  }
}
\details{
The turbofan engine data is from the Prognostic Center of Excellence
(PCoE) of NASA Ames Research Center,
which is simulated by the Commercial Modular Aero-Propulsion
System Simulation (C-MAPSS). Only 14 out of 21 variables, by a method mentioned by Li, et al. (2019) are selected. The \code{train} and \code{test} lists are of class \code{"hhsmm.data"}.
}
\references{
Frederick, D. K., DeCastro, J. A., & Litt, J. S. (2007). User's guide for the commercial modular aero-propulsion system simulation (C-MAPSS).

Saxena, A., Goebel, K., Simon, D., & Eklund, N. (2008, October). Damage propagation modeling for aircraft engine run-to-failure simulation. In \emph{2008 international conference on prognostics and health management} (pp. 1-9). IEEE.

Li, J., Li, X., & He, D. (2019). A directed acyclic graph network combined with CNN and LSTM for remaining useful life prediction. \emph{IEEE Access}, 7, 75464-75475.
}
\examples{
data(CMAPSS)
str(CMAPSS$train)
str(CMAPSS$test)
}
\keyword{datasets}
