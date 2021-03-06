% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pcalc.R
\name{pcalc}
\alias{pcalc}
\title{Perform purity calculation on a peak matrix}
\usage{
pcalc(
  peaks,
  mzmin,
  mzmax,
  mztarget,
  ppm = NA,
  iwNorm = FALSE,
  iwNormFun = NULL,
  ilim = 0,
  targetMinMZ = NA,
  targetMaxMZ = NA,
  isotopes = FALSE,
  im = NULL
)
}
\arguments{
\item{peaks}{matrix; Matrix of peaks consisting of 2 columns: mz and i}

\item{mzmin}{numeric; Isolation window (min)}

\item{mzmax}{numeric; Isolation window (max)}

\item{mztarget}{numeric; The mz window to target in the isolation window}

\item{ppm}{numeric; PPM tolerance for the target mz value. If NA will presume targetMinMZ and targetMaxMZ will be used}

\item{iwNorm}{boolean; If TRUE then the intensity of the isolation window will be normalised based on the iwNormFun function}

\item{iwNormFun}{function; A function to normalise the isolation window intensity. The default function is very generalised and just accounts for edge effects}

\item{ilim}{numeric; All peaks less than this percentage of the target peak will be removed from the purity calculation, default is 5\% (0.05)}

\item{targetMinMZ}{numeric; Range to look for the mztarget (min)}

\item{targetMaxMZ}{numeric; Range to look for the mztarget (max)}

\item{isotopes}{boolean; TRUE if isotopes are to be removed}

\item{im}{matrix; Isotope matrix, default removes C13 isotopes (single, double and triple bonds)}
}
\value{
a vector of the purity score and the number of peaks in the window e.g c(purity, pknm)
}
\description{
This is the main purity calculation that is performed in purityX,
purityD and purityA.
\itemize{
\item{Takes in a matrix of peaks}
\item{gets isolation window based on mzmin mzmax}
\item{locates the mz target in the peak matrix}
\item{removes isotopic peaks}
\item{removes any peaks below limit (percentage of target peak intensity)}
\item{normalises}
\item{Calculates purity: Divides the target peak intensity by the total peak intensity for
the isolation window}
}
}
\examples{
pm <- rbind(c(100, 1000),c(101.003, 10))
pcalc(pm, mzmin = 98, mzmax = 102, mztarget=100, ppm=5)
pcalc(pm, mzmin = 98, mzmax = 102, mztarget=100, ppm=5, isotopes = TRUE)

}
