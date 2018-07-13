\name{DTWUMI-package}
\alias{DTWUMI-package}
\alias{DTWUMI}
\docType{package}
\title{
\packageTitle{DTWUMI}
}
\description{
\packageDescription{DTWUMI}
}
\details{
%\packageDESCRIPTION{DTWUMI}
\packageIndices{DTWUMI}
}
\author{
\packageAuthor{DTWUMI}

Maintainer: \packageMaintainer{DTWUMI}
}
\references{
Thi-Thu-Hong Phan, Emilie Poisson-Caillault, Alain Lefebvre, Andre Bigand. Dynamic time warping-based imputation for univariate time series data.  Pattern Recognition Letters, Elsevier, 2017, <DOI:10.1016/j.patrec.2017.08.019>. <hal-01609256>
}
\keyword{package}
\keyword{DTW}
\keyword{times series}
\keyword{imputation}
\examples{
data(dataDTWUMI)
dataDTWUMI_gap <- dataDTWUMI[["incomplete_signal"]]
imputation <- DTWUMI_imputation(dataDTWUMI_gap, gap_size_threshold = 10, DTW_method = "DTW")
plot(dataDTWUMI_gap[, 1], type = "l", lwd = 2)
lines(imputation$output[, 1], col = "red")
plot(dataDTWUMI_gap[, 2], type = "l", lwd = 2)
lines(imputation$output[, 2], col = "red")
plot(dataDTWUMI_gap[, 3], type = "l", lwd = 2)
lines(imputation$output[, 3], col = "red")
}
