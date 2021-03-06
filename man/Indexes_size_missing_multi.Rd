% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/4_2Indexes_size_missing.R
\name{Indexes_size_missing_multi}
\alias{Indexes_size_missing_multi}
\title{Indexing gaps size}
\usage{
Indexes_size_missing_multi(data)
}
\arguments{
\item{data}{multivariate signal}
}
\value{
returns a list with one element per signal.
Within each element of this list, the first column gives the position of the begining of each gap and the second column its size.
}
\description{
Stores the position of the begining of each gap and their respective size within a multivariate signal.
}
\examples{
data(dataDTWUMI)
id_NA <- Indexes_size_missing_multi(dataDTWUMI$incomplete_signal)
}
\author{
DEZECACHE Camille, PHAN Thi Thu Hong, POISSON-CAILLAULT Emilie
}
