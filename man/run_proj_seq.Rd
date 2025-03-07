% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_proj_seq.r
\name{run_proj_seq}
\alias{run_proj_seq}
\title{Sequence of projection and estimation}
\usage{
run_proj_seq(
  CommonName = NULL,
  fileName = "bam",
  dat_obj = NULL,
  tpl_obj = NULL,
  cxx_obj = NULL,
  rdat_obj = NULL,
  run_proj_args = list(),
  nyps = NULL,
  plot_type = "final",
  write_bam_files = "final",
  parallel = TRUE,
  ncores = NULL
)
}
\arguments{
\item{CommonName}{Common name of species associated with dat, tpl, and cxx files}

\item{fileName}{Name given to BAM files, not including file extensions.}

\item{dat_obj}{dat file read in as a character vector with readLines(con=dat_file)}

\item{tpl_obj}{tpl file read in as a character vector with readLines(con=tpl_file)}

\item{cxx_obj}{cxx file read in as a character vector with readLines(con=cxx_file)}

\item{run_proj_args}{List of arguments to pass to \code{run_proj}}

\item{nyps}{Total number of years to project which may include multiple assessment cycles. nyps = Number of Years to Project in Sequence}

\item{plot_type}{Indicate which bam results should be plotted by calling \code{plot_bam}. Current options
include "all", "final", or "none", to plot results from all stock assessments in the projection
sequence, the final assessment, or none of the assessments.}

\item{write_bam_files}{indicate which bam files should be written. Current options
include "all", "final", or "none", to write files from all stock assessments in the projection
sequence, the final assessment, or none of the assessments.}

\item{parallel}{logical. Should parallel processing be used? This setting is
only used when \code{dat_obj} is a list of multiple dat objects. The function
will not use parallel processing for a single dat object (it would actually
be somewhat slower).}

\item{ncores}{number of cores to use for parallel processing}

\item{rdat}{(list) object read in with dget(). Specifically, the BAM output file produced by the cxx file.}
}
\value{
By default, returns the bam model object from the end of the projection
sequence.
}
\description{
This function runs a sequence of projections with \code{run_proj}, and estimation
with \code{run_bam}, to build long-term closed-loop simulations. If a single
\code{dat} object is provided, it will run a single deterministic projection.
If a list of \code{dat} objects is provided, for example as output from
\code{run_MCBE(run_est=FALSE)}, the function will run
multiple projections in parallel.
}
\author{
Nikolai Klibansky
}
\keyword{assessment}
\keyword{bam}
\keyword{dynamics}
\keyword{fisheries}
\keyword{population}
\keyword{stock}
