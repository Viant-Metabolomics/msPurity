% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/flag-filter-remove.R
\name{flag_remove}
\alias{flag_remove}
\title{Flag and remove unwanted peaks}
\usage{
flag_remove(
  xset,
  pol = NA,
  rsd_i_blank = NA,
  minfrac_blank = 0.5,
  rsd_rt_blank = NA,
  ithres_blank = NA,
  s2b = 10,
  ref.class = "blank",
  egauss_thr = NA,
  rsd_i_sample = NA,
  minfrac_sample = 0.7,
  rsd_rt_sample = NA,
  ithres_sample = NA,
  minfrac_xcms = 0.7,
  mzwid = 0.025,
  bw = 5,
  out_dir = ".",
  temp_save = FALSE,
  remove_spectra = TRUE,
  grp_rm_ids = NA
)
}
\arguments{
\item{xset}{object; xcmsSet object}

\item{pol}{str; polarity (just used for naming purpose for files being saved) [positive, negative, NA]}

\item{rsd_i_blank}{numeric; RSD threshold for the blank}

\item{minfrac_blank}{numeric; minimum fraction of files for features needed for the blank}

\item{rsd_rt_blank}{numeric; RSD threshold for the RT of the blank}

\item{ithres_blank}{numeric; Intensity threshold for the blank}

\item{s2b}{numeric; fold change (sample/blank) needed for sample peak to be allowed. e.g.
if s2b set to 10 and the recorded sample 'intensity' value was 100 and blank = 10.
1000/10 = 100 so sample has fold change higher than the threshold and the peak is not considered a blank}

\item{ref.class}{str; A string representing the class that will be used for the blank.}

\item{egauss_thr}{numeric; Threshold for filtering out non gaussian shaped peaks. Note this only works
if the verbose option was set for XCMS;}

\item{rsd_i_sample}{numeric; RSD threshold for the sample}

\item{minfrac_sample}{numeric; minimum fraction of files for features needed for the sample}

\item{rsd_rt_sample}{numeric; RSD threshold for the RT of the sample}

\item{ithres_sample}{numeric; Intensity threshold for the sample}

\item{minfrac_xcms}{numeric; minfrac for xcms  grouping}

\item{mzwid}{numeric; xcms grouping parameter}

\item{bw}{numeric; xcms grouping parameter}

\item{out_dir}{str; out directory}

\item{temp_save}{boolean; Assign True if files for each step saved (for testing purpsoses)}

\item{remove_spectra}{bool; TRUE if flagged spectra is to be removed}

\item{grp_rm_ids}{vector; vector of grouped_xcms peaks to remove (coresponds to the row from xcms::group output)}
}
\value{
list(xset, grp_peaklist, removed_peaks)
}
\description{
On an xcmsSet object, filter flag and remove unwanted peaks. When the peaks are removed, the
the xcmsSet object can be regrouped using xcms::group. The function then checks if any blank
peaks are still present and the process is repeated.

The output is a list of the updated xcmsSet object, grouped peaklist and the blank removed peaks
}
\examples{

msPths <- list.files(system.file("extdata", "lcms", "mzML",
                     package="msPurityData"), full.names = TRUE)
xset <- xcms::xcmsSet(msPths)
xset@phenoData[,1] <- c('blank', 'blank', 'sample', 'sample')
xset <- xcms::group(xset)
fr = flag_remove(xset)
}
