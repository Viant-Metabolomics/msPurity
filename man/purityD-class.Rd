% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/purityD-class.R
\docType{class}
\name{purityD-class}
\alias{purityD-class}
\alias{purityD}
\title{An S4 class to represent a DI-MS purityD}
\arguments{
\item{.Object}{object; purityD object}

\item{fileList}{data.frame; Created using GetFiles, data.frame with filepaths and sample class information}

\item{cores}{numeric; Number of cores used to perform Hierarchical clustering WARNING: memory intensive, default 1}

\item{mzML}{boolean; TRUE if mzML to be used FALSE if .csv file to be used}
}
\value{
purityD object
}
\description{
The class used to assess anticipated purity from a DI-MS run
}
\examples{
datapth <- system.file("extdata", "dims", "mzML", package="msPurityData")
inDF <- Getfiles(datapth, pattern=".mzML", check = FALSE, cStrt = FALSE)
ppDIMS <- purityD(fileList=inDF, cores=1, mzML=TRUE)
}
