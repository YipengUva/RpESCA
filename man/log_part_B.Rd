% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/log_partitions_functions.R
\name{log_part_B}
\alias{log_part_B}
\title{The log-partation function of Bernoulli distribution}
\usage{
log_part_B(Theta)
}
\arguments{
\item{Theta}{a matrix of natural parameter}
}
\value{
The value of the corresponding log-partation function
}
\description{
This function define the log-partation function
of Bernoulli distribution. The formula is
\code{log_part_B(Theta) = log(1+exp(Theta))}.
Details can be found in \url{https://arxiv.org/abs/1902.06241}.
}
\examples{
\dontrun{log_part_B(matrix(0:9,2,5))}
}
