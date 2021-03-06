% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modules_pESCA_algorithm.R
\name{penalty_concave_composite}
\alias{penalty_concave_composite}
\title{Composition of group-wise and element-wise conave penalty}
\usage{
penalty_concave_composite(B_i, fun_concave, gamma, R)
}
\arguments{
\item{B_i}{The loading matrix for the ith data set}

\item{fun_concave}{A string indicates the used concave function}

\item{gamma}{The hyper-parameter of the concave penalty}

\item{R}{The number of PCs}
}
\value{
This function returns the value of the composition
of group-wise and element-wise conave penalty for the pESCA model
}
\description{
This is an intermediate step of the algorithm for fitting pESCA model. The
details of this function can be found in ref thesis.
}
\examples{
\dontrun{
penalty_concave_composite(B_i, fun_concave, gamma, R)
}
}
