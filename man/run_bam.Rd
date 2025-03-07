% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_bam.R
\name{run_bam}
\alias{run_bam}
\title{Run bam model}
\usage{
run_bam(
  CommonName = NULL,
  fileName = "bam",
  dir_bam = NULL,
  bam = NULL,
  dat_file = NULL,
  tpl_file = NULL,
  cxx_file = NULL,
  exe_file = NULL,
  dat_obj = NULL,
  tpl_obj = NULL,
  cxx_obj = NULL,
  standardize = TRUE,
  subset_rdat = list(eq.series = 101, pr.series = 101, t.series = "styr:endyr"),
  unlink_dir_bam = TRUE,
  admb_options = "-nox",
  admb2r_obj = admb2r.cpp,
  cleanup = list(del = c("*.r0*", "*.p0*", "*.b0*", "*.log", "*.rpt", "*.obj", "*.htp",
    "*.eva", "*.bar", "*.tds", "*.o", "tmp_admb", "variance", "*.dep", "*.hes", "*.tmp")),
  run_cleanup = TRUE,
  shell_args_compile = list(intern = FALSE),
  shell_args_run = list(intern = FALSE),
  return_obj = "all"
)
}
\arguments{
\item{CommonName}{Common name of species associated with dat, tpl, and cxx files}

\item{fileName}{Name given to BAM files, not including file extensions.}

\item{dir_bam}{Name of directory to write BAM files to, relative to the working directory.}

\item{bam}{Output of \code{bam2r}.}

\item{dat_file}{dat file path}

\item{tpl_file}{tpl file path}

\item{cxx_file}{cxx file path}

\item{exe_file}{exe file path. If you specify only exe_file and dat_file or dat_obj,
the function will copy the executable and will not try to compile bam. It will just
run teh executable with the dat file you give it. This option is used with
run_MCBE.}

\item{dat_obj}{dat file read in as a character vector with readLines(con=dat_file)}

\item{tpl_obj}{tpl file read in as a character vector with readLines(con=tpl_file)}

\item{cxx_obj}{cxx file read in as a character vector with readLines(con=cxx_file)}

\item{standardize}{Should \code{\link[bamExtras]{standardize_bam}} be run by the function before running the BAM}

\item{subset_rdat}{list of rdat objects to subset in different ways. For eq.series and pr.series
specify number of evenly spaced values to retain. For t.series specify a series
of years to retain, either as a numeric vector, or as a character string that will be
evaluated as code internally with \code{eval(parse(text=my_string))}.
When subsetting eq.series or pr.series, this option can
substantially decrease rdat file size, without affecting precision of reference
point calculations. This is particularly helpful for the MCBE runs.}

\item{unlink_dir_bam}{Should \code{dir_bam} be deleted after this function is run?}

\item{admb_options}{Character string pasted to fileName to build \code{run_command} when running BAM with \code{shell(run_command)}. See \href{https://www.admb-project.org/docs/refcards/admb-additional-reference-card.pdf}{ADMB reference card} for more options.
(i.e. \code{run_command <- paste(fileName, admb_options)})}

\item{admb2r_obj}{Character string containing admb2r C++ code, which is written with \code{base::writeLines} to \code{dir_bam}}

\item{cleanup}{List object written to \code{cleanup.bat} file in \code{dir_bam}.}

\item{run_cleanup}{Should the \code{cleanup.bat} be run?}

\item{return_obj}{names of objects to return from the function. May include one or more
of the following values which refer to ADMB file extensions:
\code{dat}, \code{tpl}, \code{cxx}, \code{rdat} or \code{par}. These files will be
read with \code{readLines} and passed to the result as character vectors. May
also include \code{admb} which returns a list of basic run results including
\code{lk_total} (total likelihood), \code{grad_max} (maximum gradient), and
\code{file_names} (list of file names in \code{dir_bam} after running
\code{\link[base]{shell}}).
Character vector. If NULL, all objects are returned.}
}
\value{
See \code{return_obj}
}
\description{
Creates a folder in a directory \code{dir_bam} to write BAM files including model-specific files
(.dat, .tpl, and .cxx) a standard \code{admb2r.cpp} file stored in \code{bamExtras}, and a standard
but customizable \code{cleanup.bat} file. Model-specific files are supplied by the user in any one of several different ways
(see Arguments). It temporarily changes the working directory to \code{dir_bam}
then calls a shell script to run the BAM model, and changes the working directory back to the previous path.
The function invisibly returns the bam file objects (indicated by return_obj) in a list.
}
\examples{
\dontrun{
# Run a bam model and assign rdat output to object
rdat_AtMe <- run_bam("AtlanticMenhaden")$rdat
rdat_BlSB <- run_bam("BlackSeaBass")$rdat
rdat_BlTi <- run_bam("BluelineTilefish")$rdat
rdat_Cobi <- run_bam("Cobia")$rdat
rdat_GagG <- run_bam("GagGrouper")$rdat
rdat_GrTr <- run_bam("GrayTriggerfish")$rdat
rdat_GrAm <- run_bam("GreaterAmberjack")$rdat
rdat_ReGr <- run_bam("RedGrouper")$rdat
rdat_RePo <- run_bam("RedPorgy")$rdat
rdat_ReSn <- run_bam("RedSnapper")$rdat
rdat_ScGr <- run_bam("ScampGrouper")$rdat
rdat_SnGr <- run_bam("SnowyGrouper")$rdat
rdat_Tile <- run_bam("Tilefish")$rdat
rdat_VeSn <- run_bam("VermilionSnapper")$rdat

### Change a value in the init object, and run bam
bam_RePo <- bam_RePo_raw <- bam2r("RedPorgy")
bam_RePo$init$obs_cpue_sCT[1] <- paste(2*as.numeric(bam_RePo$init$obs_cpue_sCT[1]))
bam_RePo <- bam2r(init=bam_RePo$init,dat_obj = bam_RePo$dat,tpl_obj = bam_RePo$tpl,cxx_obj = bam_RePo$cxx)
bamout_RePo <- run_bam(bam=bam_RePo,return_obj = NULL,
                       shell_args_compile = list("intern"=TRUE),
                       shell_args_run = list("intern"=TRUE)
)

##### Create and run examples where errors occur #####
### silly value in dat file. Data file NOT READ CORRECTLY.
# Note that although this run fails, it is not a shell error (i.e. shell error code = 0)
bam_RePo_err1 <- bam_RePo_raw
bam_RePo_err1$init$endyr <- "pickles"
bam_RePo_err1 <- bam2r(init=bam_RePo_err1$init,dat_obj = bam_RePo_err1$dat,tpl_obj = bam_RePo_err1$tpl,cxx_obj = bam_RePo_err1$cxx)
bamout_RePo_err1 <- run_bam(bam=bam_RePo_err1,return_obj = NULL,
                            shell_args_compile = list("intern"=TRUE),
                            shell_args_run = list("intern"=TRUE)
)
# look at results to see the Data File reading error
bamout_RePo_err1$admb$run_out[grepl("Data File",bamout_RePo_err1$admb$run_out)]

### Negative value in cpue. Data file read correctly, but shell(run_command) fails with error code 1
bam_RePo_err2 <- bam_RePo_raw
bam_RePo_err2$init$obs_cpue_sCT[1] <- paste(-10)
bam_RePo_err2 <- bam2r(init=bam_RePo_err2$init,dat_obj = bam_RePo_err2$dat,tpl_obj = bam_RePo_err2$tpl,cxx_obj = bam_RePo_err2$cxx)
bamout_RePo_err2 <- run_bam(bam=bam_RePo_err2,return_obj = NULL,
                            shell_args_compile = list("intern"=TRUE),
                            shell_args_run = list("intern"=TRUE)
)
# look at results to see that the Data File was read correctly
bamout_RePo_err2$admb$run_out[grepl("Data File",bamout_RePo_err2$admb$run_out)]
# get run error code
bamout_RePo_err2$admb$error_code_run
# also note that lk_total is NaN
bamout_RePo_err2$admb$lk_total

### Crazy values of M. Data file read correctly, no shell error, but Hessian doesn't converge
bam_RePo_err3 <- bam_RePo_raw
bam_RePo_err3$init$set_M <- setNames(paste(as.numeric(bam_RePo_err3$init$set_M)*10),names(bam_RePo_err3$init$set_M))
bam_RePo_err3 <- bam2r(init=bam_RePo_err3$init,dat_obj = bam_RePo_err3$dat,tpl_obj = bam_RePo_err3$tpl,cxx_obj = bam_RePo_err3$cxx)
bamout_RePo_err3 <- run_bam(bam=bam_RePo_err3,return_obj = NULL,
                            shell_args_compile = list("intern"=TRUE),
                            shell_args_run = list("intern"=TRUE)
)
# look at results to see the indefinite Hessian error
bamout_RePo_err3$admb$run_out[grepl("Hessian",bamout_RePo_err3$admb$run_out)]
}
}
\keyword{assessment}
\keyword{bam}
\keyword{fisheries}
\keyword{stock}
