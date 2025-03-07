% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_MCBE.R
\name{run_MCBE}
\alias{run_MCBE}
\title{Run MCBE uncertainty analysis}
\usage{
run_MCBE(
  CommonName = NULL,
  fileName = "bam",
  dir_bam_sim = "sim",
  dir_bam_base = "base",
  bam = NULL,
  dat_file = NULL,
  tpl_file = NULL,
  cxx_file = NULL,
  dat_obj = NULL,
  tpl_obj = NULL,
  cxx_obj = NULL,
  data_sim = list(cv_U = NULL, cv_L = NULL, cv_D = NULL),
  par_default = list(cv_U = 0.2, cv_L = 0.2, cv_D = 0.2),
  standardize = TRUE,
  nsim = 11,
  sclim = list(),
  sclim_gen = c(0.9, 1.1),
  data_type_resamp = c("U", "L", "D", "age", "len"),
  fn_par = list(M_constant = "runif(nsim,M_min,M_max)", K = "runif(nsim,K_min,K_max)",
    Linf = "runif(nsim,Linf_min,Linf_max)", t0 = "runif(nsim,t0_min,t0_max)", steep =
    "runif(nsim,steep_min,steep_max)", rec_sigma =
    "msm::rtnorm(n=nsim,mean=0.6,sd=0.15,lower=0.3,upper=1.0)", Dmort =
    "apply(Dmort_lim,2,function(x){runif(nsim,min(x),max(x))})", Pfa =
    "runif(nsim,Pfa_min,Pfa_max)", Pfb = "runif(nsim,Pfb_min,Pfb_max)", Pfma =
    "runif(nsim,Pfma_min,Pfma_max)", Pfmb = "runif(nsim,Pfmb_min,Pfmb_max)", fecpar_a =
    "runif(nsim,fecpar_a_min,fecpar_a_max)", 
     fecpar_b =
    "runif(nsim,fecpar_b_min,fecpar_b_max)", fecpar_c =
    "runif(nsim,fecpar_c_min,fecpar_c_max)", fecpar_batches_sc =
    "runif(nsim,sclim_gen[1],sclim_gen[1])"),
  repro_sim = list(P_repro = 0.1, nagecfishage = 10, plot_fits = FALSE),
  fix_par = c(),
  ncores = NULL,
  ndigits = 4,
  unlink_dir_bam_base = FALSE,
  unlink_dir_bam_sim = FALSE,
  run_bam_base = TRUE,
  overwrite_bam_base = TRUE,
  admb_options_base = "-nox",
  run_est = TRUE,
  sim_type_return = "dat",
  admb_options_sim = "-est -nox",
  prompt_me = FALSE,
  subset_rdat = list(eq.series = 101, pr.series = 101, t.series = "styr:endyr"),
  random_seed = 12345
)
}
\arguments{
\item{CommonName}{Common name of species associated with dat, tpl, and cxx files}

\item{fileName}{Name given to BAM files, not including file extensions.}

\item{dir_bam_sim}{Name of directory to write MCBE files to, relative to the working directory.}

\item{dir_bam_base}{Name of directory to write bam base model files to, relative to the working directory.}

\item{bam}{Output of \code{bam2r}.}

\item{dat_file}{dat file path}

\item{tpl_file}{tpl file path}

\item{cxx_file}{cxx file path}

\item{dat_obj}{dat file read in as a character vector with readLines(con=dat_file)}

\item{tpl_obj}{tpl file read in as a character vector with readLines(con=tpl_file)}

\item{cxx_obj}{cxx file read in as a character vector with readLines(con=cxx_file)}

\item{data_sim}{list for supplying optional data not available in input data or rdat (e.g. cvs for simulating time series of landings, discards, and cpue)}

\item{par_default}{list for supplying default values for particular parameters to use if values cannot be found by in the usual locations (e.g. if a time series of landings does not have a corresponding time series of cvs, the default cv_L will be used)}

\item{standardize}{Should \code{\link[bamExtras]{standardize_bam}} be run by the function before running the BAM?}

\item{nsim}{number of simulations to run}

\item{sclim}{Optional list of scalars for computing point parameter limits. By default, limits are (generically) computed using sclim_gen. Numeric vectors of length 2, usually centered around 1. e.g. sclim = list(M_constant=c(0.9,1.1), steep=(c(0.8,1.2))). Note that if M_constant is not available in the base model output, sclim$M_constant values will be used to scale M at age.}

\item{sclim_gen}{Scalar (multipliers) for computing upper and lower bounds of random uniform distribution from mean value from base run output}

\item{data_type_resamp}{character vector of abbreviations for types of data sets that should be resampled in the MCBE simulations. L = landings, D = discards, U = cpue indices, age = age compositions, len = length compositions. If you don't want to resample any of the data sets, set data_type_resamp = c().}

\item{fn_par}{List of character strings used to simulate values of fixed parameters (see \code{Details}). Strings are internally passed to \code{eval(parse(text=mystring))}. Functions should produce vectors of length nsim, or in some cases matrices with nrow nsim.}

\item{repro_sim}{List of settings to be used to simulate reproductive data sets. Set to NULL to shut it off. \code{P_repro} defines the proportion of age comps for which reproductive data are available.  \code{nagecfishage} is the assumed number of age samples per age, used only if age comp data is not found in BAM.}

\item{fix_par}{Optional character vector of parameter names to fix in the simulations using tpl init object names with a phase setting (e.g. "set_M_constant", "set_steep"). This is mostly used for running sensitivities and parameter profiles. Note that this has no effect on the base model.}

\item{ncores}{number of cores to use for parallel processing}

\item{ndigits}{number of digits to round simulated values to}

\item{unlink_dir_bam_base}{Should dir_bam_base be deleted after this function is run?}

\item{unlink_dir_bam_sim}{Should dir_bam_sim be deleted after this function is run?}

\item{run_bam_base}{If FALSE, the function will look for an executable named fileName.exe in dir_bam_base and use it as the base model.}

\item{overwrite_bam_base}{If FALSE, the files in dir_bam_base will not be overwritten if run_bam_base=TRUE}

\item{admb_options_base}{Character string pasted to fileName to build \code{run_command} for the base model when running BAM with \code{shell(run_command)}
(i.e. \code{run_command <- paste(fileName, admb_options)})}

\item{run_est}{Do you want to run the estimation model (BAM)? If FALSE, the simulated data will be generated but won't be used in new BAM runs}

\item{sim_type_return}{Type of simulated data to return: "init" or "dat". Has no effect in run_est=TRUE.}

\item{admb_options_sim}{ADMB code snippet used in shell script when running bam}

\item{subset_rdat}{list of rdat objects to subset in different ways. For eq.series and pr.series
specify number of evenly spaced values to retain. For t.series specify a series
of years to retain. When subsetting eq.series or pr.series, this option can
substantially decrease rdat file size, without affecting precision of reference
point calculations. This is particularly helpful for the MCBE runs.}

\item{random_seed}{random seed value. If NULL, random seed is not set within the function.}
}
\value{
Invisibly returns a data frame containing basic results of sim runs, including the following objects
typically stored in the rdat output from BAM: \code{parms}, \code{parm.cons} (estimated values only), \code{like}, and \code{sdnr}. This data frame is also written to a csv file in \code{dir_bam_sim}.
}
\description{
This function reruns the base model, from user supplied input files or objects, creating all standard BAM ADMB files
most notably the executable (\code{exe}) file used in the MCBE process, stored in \code{dir_bam_base}.
Simulated input data sets are generated for specified fixed parameters and data sets. New BAM data input (\code{dat}) files
are created and the BAM base executable is rerun with each \code{dat} file, using parallel processing. After each simulation
is run, \code{run_bam} checks if the objective function value (i.e. total likelihood) for a numeric value; if the value is non numeric
(e.g. nan), then the run is considered to have failed. The \code{dat} and results (\code{rdat}) files are for successful
runs are moved to \code{dir_bam_sim}. If any runs fail (which is not common) a folder is created (named \code{dir_bam_sim}
with suffix '_fail') and failed runs are moved there. \code{run_bam}
and results file (\code{rdat}) are copied to \code{dir_bam_sim}.
}
\details{
\strong{Estimating reproductive (sex and maturity) functions}
\cr
Sex ratio (proportion of females; \code{obs_prop_f}) and female maturity (proportion of females mature;
\code{obs_maturity_f}) at age are typically provided to BAM data as vectors of proportions at age. In
most cases, those proportions have been computed with a logistic model. The \code{run_MCBE} function
characterizes uncertainty in those relationships by first building a composite age-frequency distribution
of all age samples used in to the assessment (in numbers) multiplied by an assumed proportion
(\code{repro_sim$P_repro}) of age samples which also have sex and maturity estimates to approximate the
number of fish per age for which sex was determined (\code{nrepro}. \code{nrepro} is then used as the number
of trials-at-age while successes-at-age approximated by \code{nrepro} times \code{obs_prop_f}, in logistic
regression of the form \eqn{1/(1+exp(-a*(x-b)))}. The computed number of females-at-age (\code{nf}) is then used as the
numbers of trials in a logistic regression to estimate a continuous function predicting proportion of females
mature. Uncertainty or alternate values of parameters \code{a} (slope) and \code{b} (a50) of either relationship
can then be specified by \code{fn_par$Pfa}, \code{fn_par$Pfb}, \code{fn_par$Pfma}, or \code{fn_par$Pfmb} when running the MCBE analysis
to simulate uncertainty in sex- or female maturity-at-age.
\cr \cr
There are two special cases in which the function will not try to estimate a logistic fit: \code{obs_prop_f} or
\code{obs_maturity_f} is either constant (e.g. all equal to 0.5) or completely binary (i.e. only 0 or 1). In the constant
case the \code{a} parameter can be used to simulate variation in constant sex or maturity. In the binary case, the \code{b} parameter
can be used to simulate variation in \code{a50}.
\cr \cr
Note that this process doesn't currently fit reproductive data in a time varying manner,
so if \code{obs_prop_f} or \code{obs_maturity_f} is a time-varying matrix,
a logistic model will only be fit to the last year of data, and the results
will be repeated for all years of the matrix.

\strong{Values currently supported by \code{fn_par}}
\tabular{ll}{
\code{M_constant} \tab Natural mortality constant. Corresponds to the \code{set_M_constant} parameter in the BAM tpl file. \cr
\code{K} \tab Von-Bertalanffy growth model \code{K} K parameter. Corresponds to the \code{set_K} parameter in the BAM tpl file. \cr
\code{Linf} \tab Von-Bertalanffy growth model \code{Linf} parameter. Corresponds to the \code{set_Linf} parameter in the BAM tpl file. \cr
\code{steep} \tab Beverton-Holt stock recruitment function \code{h} (steepness) parameter. Corresponds to the \code{set_steep} parameter in the BAM tpl file. \cr
\code{R0} \tab  Beverton-Holt stock recruitment function \code{R0} (unfished recruitment) parameter. Corresponds to the \code{set_rec_sigma} parameter in the BAM tpl file. \cr
\code{Dmort} \tab Discard mortality rate parameter(s). Corresponds to the parameters beginning with \code{set_Dmort} in the BAM tpl file. \cr
\code{Pfa, Pfb} \tab Scalar for sex ratio model parameter \code{a} or \code{b}. \code{run_MCBE} fits a logistic function of the form \eqn{1/(1+exp(-a*(x-b)))} to \code{obs_prop_f} (see Details).\cr
\code{Pfma, Pfmb} \tab Scalar for female maturity model parameter \code{a} or \code{b}. \code{run_MCBE} fits a logistic function of the form \eqn{1/(1+exp(-a*(x-b)))} to \code{obs_maturity_f}  (see Details).\cr

}
}
\examples{
\dontrun{
Run MCBE, writing files to dir_bam_sim
run_MCBE("AtlanticMenhaden", dir_bam_base="AtMe_base", dir_bam_sim="AtMe_sim")
run_MCBE("BlackSeaBass", dir_bam_base="BlSB_base", dir_bam_sim="BlSB_sim")
run_MCBE("BluelineTilefish", dir_bam_base="BlTi_base", dir_bam_sim="BlTi_sim")
run_MCBE("Cobia", dir_bam_base="Cobi_base", dir_bam_sim="Cobi_sim")
run_MCBE("GagGrouper", dir_bam_base="GaGr_base", dir_bam_sim="GaGr_sim")
run_MCBE("GrayTriggerfish", dir_bam_base="GrTr_base", dir_bam_sim="GrTr_sim")
run_MCBE("GreaterAmberjack", dir_bam_base="GrAm_base", dir_bam_sim="GrAm_sim")
run_MCBE("RedGrouper", dir_bam_base="ReGr_base", dir_bam_sim="ReGr_sim")
run_MCBE("RedPorgy", dir_bam_base="RePo_base", dir_bam_sim="RePo_sim")
run_MCBE("RedSnapper", dir_bam_base="ReSn_base", dir_bam_sim="ReSn_sim")
run_MCBE("ScampGrouper", dir_bam_base="ScGr_base", dir_bam_sim="ScGr_sim")
run_MCBE("SnowyGrouper", dir_bam_base="SnGr_base", dir_bam_sim="SnGr_sim")
run_MCBE("Tilefish", dir_bam_base="Tile_base", dir_bam_sim="Tile_sim")
run_MCBE("VermilionSnapper", dir_bam_base="VeSn_base", dir_bam_sim="VeSn_sim")

MCBE_ReSn <- run_MCBE("RedSnapper",steep_sclim = c(1,1))
}

}
\keyword{MCBE}
\keyword{assessment}
\keyword{bam}
\keyword{fisheries}
\keyword{stock}
