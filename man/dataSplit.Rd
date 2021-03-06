% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/supplementary_functions.R
\name{dataSplit}
\alias{dataSplit}
\title{Split multiple data sets into training and test sets}
\usage{
dataSplit(dataSets, dataTypes, ratio_mis = 0.1)
}
\arguments{
\item{dataSets}{a list contains multiple matrices with same number of rows.
Each matrix (\code{samples * variables}) indicates a data set.}

\item{dataTypes}{a string indicates the data types of the multiple data sets.}

\item{ratio_mis}{how many percent of test set could be? default: 0.1}
}
\value{
This function returns a list contains \itemize{
\item trainSets: a list contains the training sets;
\item testSets: a list contains the test sets;
\item indexSets: a list contains the index sets.
}
}
\description{
This function will split multiple data sets into training and test
sets. Nonmissing elements are randomly selected as the test sets.
Then the selected elements are taken as missing, and regarded as
training sets. The details can be found in \url{https://arxiv.org/abs/1902.06241}.
}
\examples{
\dontrun{dataSplit(dataSets,dataTypes,ratio_mis=0.1)}
}
