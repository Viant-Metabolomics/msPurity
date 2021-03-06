% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iw-norm.R
\name{iwNormGauss}
\alias{iwNormGauss}
\title{Gaussian normalisation for isolation window efficiency}
\usage{
iwNormGauss(sdlim = 3, minOff = -0.5, maxOff = +0.5)
}
\arguments{
\item{sdlim}{numerical; Standard deviation limit for gaussian curve}

\item{minOff}{numerical; Offset to the 'left' for the precursor range. (Should be negative)}

\item{maxOff}{character; Offset to the 'left' for the precursor range. (Should be positive)}
}
\value{
normalisation function for selected range.
}
\description{
Creates a function based on a gaussian curve shape that will normalise any intensity values within
a defined isolation window.

The function that is created will output a value between 0 to 1 based on the position between
the minOff and maxOff params. (The value 1.0 being equivalent to 100\% efficient)
}
\examples{

iwNormFun <- iwNormGauss(minOff=-0.5, maxOff=0.5)
pm <- data.frame(mz=c(99.5, 99.9, 100, 100.1, 100.5),i=c(1000, 1000, 1000, 1000, 1000))
mzmax = 100.5
mzmin = 99.5
middle <- mzmax-(mzmax-mzmin)/2
adjustmz = pm$mz-middle

# normalise the intensities
pm$normi = pm$i*iwNormFun(adjustmz)


}
