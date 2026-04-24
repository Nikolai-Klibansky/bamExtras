#' Run MCBE uncertainty analysis
#'
#' This function reruns the base model, from user supplied input files or objects, creating all standard BAM ADMB files
#' most notably the executable (\code{exe}) file used in the MCBE process, stored in \code{dir_bam_base}.
#' Simulated input data sets are generated for specified fixed parameters and data sets. New BAM data input (\code{dat}) files
#' are created and the BAM base executable is rerun with each \code{dat} file, using parallel processing. After each simulation
#' is run, \code{run_bam} checks if the objective function value (i.e. total likelihood) for a numeric value; if the value is non numeric
#' (e.g. nan), then the run is considered to have failed. The \code{dat} and results (\code{rdat}) files for successful
#' runs are moved to \code{dir_bam_sim}. If any runs fail (which is not common) a folder is created (named \code{dir_bam_sim}
#' with suffix '_fail') and failed runs are moved there. \code{run_bam}
#' and results file (\code{rdat}) are copied to \code{dir_bam_sim}.
#' @param CommonName Common name of species associated with dat, tpl, and cxx files
#' @param fileName Name given to BAM files, not including file extensions.
#' @param dir_bam_sim Name of directory to write MCBE files to, relative to the working directory.
#' @param dir_bam_base Name of directory to write bam base model files to, relative to the working directory.
#' @param bam Output of \code{bam2r}.
#' @param dat_file dat file path
#' @param tpl_file tpl file path
#' @param cxx_file cxx file path
#' @param rdat_base_file rdat_base file path
#' @param dat_obj dat object. (e.g. file read in as a character vector with readLines(con=dat_file))
#' @param tpl_obj tpl object. (e.g. file read in as a character vector with readLines(con=tpl_file))
#' @param cxx_obj cxx object. (e.g. file read in as a character vector with readLines(con=cxx_file))
#' @param rdat_base_obj rdat_base object. (e.g. file read in as list with dget(file=rdat_base_file))
#' @param nsim number of simulations to run
#' @param data_sim list for supplying optional data not available in input data or rdat (e.g. cvs for simulating time series of landings, discards, and cpue). Can be provided either as a function that will be applied internally to cvs used for fitting in BAM and a data frame where rows represent model years (styr:endyr) and columns represent fleets. Column names need to be fleet abbreviations matching the fleet abbreviations used in the tpl (e.g. cHL, rHB)
#' @param par_default list for supplying default values for particular parameters to use if values cannot be found by in the usual locations (e.g. if a time series of landings does not have a corresponding time series of cvs, the default cv_L will be used)
#' @param standardize Should \code{\link[bamExtras]{standardize_bam}} be run by the function before running the BAM?
#' @param sclim_gen Scalar (multipliers) for computing upper and lower bounds of random uniform distribution from mean value from base run output
#' @param sclim Optional list of scalars for computing point parameter limits. By default, limits are (generically) computed using sclim_gen. Numeric vectors of length 2, usually centered around 1. e.g. sclim = list(M_constant=c(0.9,1.1), steep=(c(0.8,1.2))). Note that if M_constant is not available in the base model output, sclim$M_constant values will be used to scale M at age.
#' @param data_type_resamp character vector of abbreviations for types of data sets that should be resampled in the MCBE simulations. L = landings, D = discards, U = cpue indices, age = age compositions, len = length compositions. If you don't want to resample any of the data sets, set data_type_resamp = c().
#' @param fn_par List of character strings used to simulate values of fixed parameters (see \code{Details}). Strings are internally passed to \code{eval(parse(text=mystring))}. Functions should produce vectors of length nsim, or in some cases matrices with nrow nsim.
#' @param repro_sim List of settings to be used to simulate reproductive data sets. Set to NULL to shut it off. \code{P_repro} defines the proportion of age comps for which reproductive data are available.  \code{nagecfishage} is the assumed number of age samples per age, used only if age comp data is not found in BAM.
#' @param fix_par Optional character vector of parameter names to fix in the simulations using tpl init object names with a phase setting (e.g. "set_M_constant", "set_steep"). This is mostly used for running sensitivities and parameter profiles. Note that this has no effect on the base model.
#' @param subset_rdat  Subset objects in sim rdat files. Value should be a list of the object names and either a numeric value indicating how many evenly spaced rows to include in the subset of a matrix, or NULL to set the object to NULL
#' @param ncores number of cores to use for parallel processing
#' @param ndigits number of digits to round simulated values to
#' @param unlink_dir_bam_base Should dir_bam_base be deleted after this function is run?
#' @param unlink_dir_bam_sim  Should dir_bam_sim be deleted after this function is run?
#' @param run_bam_base If FALSE, the function will look for an executable named fileName.exe in dir_bam_base and use it as the base model.
# If TRUE and overwrite_bam_base=TRUE, the function will call run_bam.
#' @param overwrite_bam_base If FALSE, the files in dir_bam_base will not be overwritten if run_bam_base=TRUE
#' @param admb_options_base Character string pasted to fileName to build \code{run_command} for the base model when running BAM with \code{shell(run_command)}
#' (i.e. \code{run_command <- paste(fileName, admb_options)})
#' @param run_est Do you want to run the estimation model (BAM)? If FALSE, the simulated data will be generated but won't be used in new BAM runs
#' @param run_mc Logical switch for turning on/off Monte Carlo (MC) draws. If FALSE, the MC process will not run. Note that this includes reproductive inputs, thus it effectively shuts of repro_sim
#' @param run_boot Logical switch for turning on/off bootstrapping (resampling) of data. If FALSE, the bootstrapping process will not run.
#' @param sim_type_return Type of simulated data to return: "init" or "dat". Has no effect if run_est=TRUE.
#' @param admb_options_sim ADMB code snippet used in shell script when running bam
#' @param subset_rdat list of rdat objects to subset in different ways. For eq.series and pr.series
#' specify number of evenly spaced values to retain. For t.series specify a series
#' of years to retain. When subsetting eq.series or pr.series, this option can
#' substantially decrease rdat file size, without affecting precision of reference
#' point calculations. This is particularly helpful for the MCBE runs.
#' @param random_seed random seed value. If NULL, random seed is not set within the function.
#' @details
#' \strong{Estimating reproductive (sex and maturity) functions}
#' \cr
#' Sex ratio (proportion of females; \code{obs_prop_f}) and female maturity (proportion of females mature;
#' \code{obs_maturity_f}) at age are typically provided to BAM data as vectors of proportions at age. In
#' most cases, those proportions have been computed with a logistic model. The \code{run_MCBE} function
#' characterizes uncertainty in those relationships by first building a composite age-frequency distribution
#' of all age samples used in to the assessment (in numbers) multiplied by an assumed proportion
#' (\code{repro_sim$P_repro}) of age samples which also have sex and maturity estimates to approximate the
#' number of fish per age for which sex was determined (\code{nrepro}. \code{nrepro} is then used as the number
#' of trials-at-age while successes-at-age approximated by \code{nrepro} times \code{obs_prop_f}, in logistic
#' regression of the form \eqn{1/(1+exp(-a*(x-b)))}. The computed number of females-at-age (\code{nf}) is then used as the
#' numbers of trials in a logistic regression to estimate a continuous function predicting proportion of females
#' mature. Uncertainty or alternate values of parameters \code{a} (slope) and \code{b} (a50) of either relationship
#' can then be specified by \code{fn_par$Pfa}, \code{fn_par$Pfb}, \code{fn_par$Pfma}, or \code{fn_par$Pfmb} when running the MCBE analysis
#' to simulate uncertainty in sex- or female maturity-at-age.
#' \cr \cr
#' There are two special cases in which the function will not try to estimate a logistic fit: \code{obs_prop_f} or
#' \code{obs_maturity_f} is either constant (e.g. all equal to 0.5) or completely binary (i.e. only 0 or 1). In the constant
#' case the \code{a} parameter can be used to simulate variation in constant sex or maturity. In the binary case, the \code{b} parameter
#' can be used to simulate variation in \code{a50}.
#' \cr \cr
#' Note that this process doesn't currently fit reproductive data in a time varying manner,
#' so if \code{obs_prop_f} or \code{obs_maturity_f} is a time-varying matrix,
#' a logistic model will only be fit to the last year of data, and the results
#' will be repeated for all years of the matrix.
#'
#' \strong{Values currently supported by \code{fn_par}}
#' \tabular{ll}{
#'  \code{M_constant} \tab Natural mortality constant. Corresponds to the \code{set_M_constant} parameter in the BAM tpl file. \cr
#'  \code{K} \tab Von-Bertalanffy growth model \code{K} K parameter. Corresponds to the \code{set_K} parameter in the BAM tpl file. \cr
#'  \code{Linf} \tab Von-Bertalanffy growth model \code{Linf} parameter. Corresponds to the \code{set_Linf} parameter in the BAM tpl file. \cr
#'  \code{steep} \tab Beverton-Holt stock recruitment function \code{h} (steepness) parameter. Corresponds to the \code{set_steep} parameter in the BAM tpl file. \cr
#'  \code{R0} \tab  Beverton-Holt stock recruitment function \code{R0} (unfished recruitment) parameter. Corresponds to the \code{set_rec_sigma} parameter in the BAM tpl file. \cr
#'  \code{Dmort} \tab Discard mortality rate parameter(s). Corresponds to the parameters beginning with \code{set_Dmort} in the BAM tpl file. \cr
#'  \code{Pfa, Pfb} \tab Scalar for sex ratio model parameter \code{a} or \code{b}. \code{run_MCBE} fits a logistic function of the form \eqn{1/(1+exp(-a*(x-b)))} to \code{obs_prop_f} (see Details).\cr
#'  \code{Pfma, Pfmb} \tab Scalar for female maturity model parameter \code{a} or \code{b}. \code{run_MCBE} fits a logistic function of the form \eqn{1/(1+exp(-a*(x-b)))} to \code{obs_maturity_f}  (see Details).\cr
#'
#'
#' }
#' @return Invisibly returns a data frame containing basic results of sim runs, including the following objects
#' typically stored in the rdat output from BAM: \code{parms}, \code{parm.cons} (estimated values only), \code{like}, and \code{sdnr}. This data frame is also written to a csv file in \code{dir_bam_sim}.
#' @keywords bam MCBE stock assessment fisheries
#' @export
#' @examples
#' \dontrun{
#' Run MCBE, writing files to dir_bam_sim
#' run_MCBE("AtlanticMenhaden", dir_bam_base="AtMe_base", dir_bam_sim="AtMe_sim")
#' run_MCBE("BlackSeaBass", dir_bam_base="BlSB_base", dir_bam_sim="BlSB_sim")
#' run_MCBE("BluelineTilefish", dir_bam_base="BlTi_base", dir_bam_sim="BlTi_sim")
#' run_MCBE("Cobia", dir_bam_base="Cobi_base", dir_bam_sim="Cobi_sim")
#' run_MCBE("GagGrouper", dir_bam_base="GaGr_base", dir_bam_sim="GaGr_sim")
#' run_MCBE("GrayTriggerfish", dir_bam_base="GrTr_base", dir_bam_sim="GrTr_sim")
#' run_MCBE("GreaterAmberjack", dir_bam_base="GrAm_base", dir_bam_sim="GrAm_sim")
#' run_MCBE("RedGrouper", dir_bam_base="ReGr_base", dir_bam_sim="ReGr_sim")
#' run_MCBE("RedPorgy", dir_bam_base="RePo_base", dir_bam_sim="RePo_sim")
#' run_MCBE("RedSnapper", dir_bam_base="ReSn_base", dir_bam_sim="ReSn_sim")
#' run_MCBE("ScampGrouper", dir_bam_base="ScGr_base", dir_bam_sim="ScGr_sim")
#' run_MCBE("SnowyGrouper", dir_bam_base="SnGr_base", dir_bam_sim="SnGr_sim")
#' run_MCBE("Tilefish", dir_bam_base="Tile_base", dir_bam_sim="Tile_sim")
#' run_MCBE("VermilionSnapper", dir_bam_base="VeSn_base", dir_bam_sim="VeSn_sim")
#'
#' MCBE_ReSn <- run_MCBE("RedSnapper",steep_sclim = c(1,1))
#'
#' # Add more observation error to landings
#' ts <- rdat_GagGrouper$t.series
#' cv_L <- ts[,names(ts)[grepl("^cv.L",names(ts))],drop=FALSE]*4
#' names(cv_L) <- gsub("^cv.L.","",names(cv_L))
#' data_sim <- list(cv_U = NULL, cv_L = cv_L, cv_D = NULL)
#' run_MCBE("GagGrouper", dir_bam_base="GaGr_base", dir_bam_sim="GaGr_sim", data_sim=data_sim)
#' simsum_GaGr <- summarize_MCBE(dir_bam_sim="GaGr_sim",dir_bam_base = "GaGr_base") # summarize results
#' plot_MCBE(sim_summary = simsum_GaGr,filter_MCBE_args = NULL,dir_figs="figs_GaGr") # plot results
#' }
#'
run_MCBE <- function(CommonName = NULL,
                     fileName     = "bam",
                     dir_bam_sim  = "sim",
                     dir_bam_base = "base",
                     bam=NULL,
                     dat_file=NULL,tpl_file=NULL,cxx_file=NULL,rdat_base_file=NULL,
                     dat_obj=NULL, tpl_obj=NULL,cxx_obj=NULL,rdat_base_obj=NULL,
                     data_sim =  list(cv_U=function(x,a=1){x*a},
                                      cv_L=function(x,a=4){x*a},
                                      cv_D=function(x,a=4){x*a}),
                     par_default = list(cv_U=0.2,
                                        cv_L=0.05,
                                        cv_D=0.05),
                     standardize=TRUE,
                     nsim=11,
                     sclim = list(),
                     sclim_gen = c(0.9,1.1),
                     data_type_resamp = c("U","L","D","age","len"),
                     fn_par = list(
                       M_constant = "runif(nsim,M_min,M_max)",
                       K = "runif(nsim,K_min,K_max)",
                       Linf = "runif(nsim,Linf_min,Linf_max)",
                       t0 = "runif(nsim,t0_min,t0_max)",
                       steep = "runif(nsim,steep_min,steep_max)",
                       rec_sigma = "msm::rtnorm(n=nsim,mean=0.6,sd=0.15,lower=0.3,upper=1.0)",
                       Dmort = "apply(Dmort_lim,2,function(x){runif(nsim,min(x),max(x))})",
                       Pfa=    "runif(nsim,Pfa_min,Pfa_max)",
                       Pfb=    "runif(nsim,Pfb_min,Pfb_max)",
                       Pfma=   "runif(nsim,Pfma_min,Pfma_max)",
                       Pfmb=   "runif(nsim,Pfmb_min,Pfmb_max)",
                       fecpar_a= "runif(nsim,fecpar_a_min,fecpar_a_max)",
                       fecpar_b= "runif(nsim,fecpar_b_min,fecpar_b_max)",
                       fecpar_c= "runif(nsim,fecpar_c_min,fecpar_c_max)",
                       fecpar_batches_sc="runif(nsim,sclim_gen[1],sclim_gen[1])"    # scalars, not actual parameter values
                     ),
                     repro_sim= list(P_repro = 0.10,
                                     nagecfishage = 10,
                                     plot_fits = FALSE
                                     ),
                     fix_par = c(),
                     # parallel=TRUE, # Right now it has to be in parallel
                     ncores=NULL,
                     ndigits=4,
                     unlink_dir_bam_base=FALSE,
                     unlink_dir_bam_sim=FALSE,
                     run_bam_base=TRUE,
                     overwrite_bam_base=TRUE,
                     admb_options_base = '-nox',
                     run_est = TRUE,
                     run_mc = TRUE,
                     run_boot = TRUE,
                     sim_type_return = "dat",
                     admb_options_sim = '-est -nox',
                     prompt_me=FALSE,
                     subset_rdat=list("eq.series"=101,"pr.series"=101,"t.series"="styr:endyr"),
                     random_seed=NULL


                     # admb2r_obj = admb2r.cpp,
                     # cleanup = list(del=c("*.r0*","*.p0*","*.b0*","*.log","*.rpt","*.obj",
                     #           "*.htp","*.eva","*.bar","*.tds","*.o","tmp_admb",
                     #           "variance","*.dep","*.hes","*.tmp"))
){
  #######################
  # library(doParallel)
  # library(foreach)
  # library(msm)

  dir_bam_sim_fail <- paste0(dir_bam_sim,"_fail")

  if(is.numeric(random_seed)){
    set.seed(random_seed)
  }

  # parallel setup
  if(is.null(ncores)){
    coresAvail <- parallel::detectCores()
    ncores <- coresAvail-1
  }
  ncores  <- min(c(ncores,coresAvail))
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)

  nm_sim <- sprintf(paste("%0",nchar(nsim),".0f",sep=""),1:nsim)

  ## Get base model stuff

  # Identify working directory
  wd <- getwd()
  message(paste("working directory:",wd))

  # Define bam base model objects
  if(!is.null(CommonName)){
    dat <- get(paste0("dat_",CommonName))
    tpl <- get(paste0("tpl_",CommonName))
    cxx <- get(paste0("cxx_",CommonName))
  }

  if(!is.null(bam)){
    dat <- bam$dat
    tpl <- bam$tpl
    cxx <- bam$cxx
  }

  if(!is.null(dat_obj)&!is.null(tpl_obj)&!is.null(cxx_obj)){
    dat <- dat_obj
    tpl <- tpl_obj
    cxx <- cxx_obj
  }

  # Read in dat, tpl, and cxx files
  if(!is.null(dat_file)&!is.null(tpl_file)&!is.null(cxx_file)){
    dat <- readLines(con=dat_file)
    tpl <- readLines(con=tpl_file)
    cxx <- readLines(con=cxx_file)
  }

  if(standardize){
    message("Running bamExtras::standardize_bam()")
    bam <- standardize_bam(dat_obj=dat, tpl_obj=tpl,cxx_obj=cxx)
    dat <- bam$dat
    tpl <- bam$tpl
    cxx <- bam$cxx
  }else{
    message("Running bamExtras::bam2r()")
    bam <- bam2r(dat_obj=dat, tpl_obj=tpl,cxx_obj=cxx)
  }
  init <- bam$init
  nm_rdat_base <- paste0(fileName,".rdat")

  # Run bam base
  if(run_bam_base){
    if(dir.exists(dir_bam_base)){
      message(paste("The folder",paste0("'",dir_bam_base,"'"),"already exists."))
      if(overwrite_bam_base){
        message(paste("Since overwrite_bam_base = TRUE, files in ",dir_bam_base,"will be overwritten when run_bam is called"))
        rdat_base <- run_bam(bam=bam, dir_bam = dir_bam_base, unlink_dir_bam=unlink_dir_bam_base,admb_options=admb_options_base)$rdat
      }else{
        message(paste("Since overwrite_bam_base = FALSE, run_bam will not be called to rerun the base run"))
        rdat_base <- dget(file.path(dir_bam_base,nm_rdat_base))
      }
    }else{
      dir.create(dir_bam_base)
      rdat_base <- run_bam(bam=bam, dir_bam = dir_bam_base, unlink_dir_bam=unlink_dir_bam_base,admb_options=admb_options_base)$rdat

    }
  }else if(!is.null(CommonName)){
    rdat_base <- get(paste0("rdat_",CommonName))
  }else{
    if(!is.null(rdat_base_file)){
      rdat_base <- dget(rdat_base_file)
      message("got user supplied rdat_base_file, to define rdat_base")
    }else if(!is.null(rdat_base_obj)){
      message("using user supplied rdat_base_obj, to define rdat_base")
      rdat_base <- rdat_base_obj
    }else if(nm_rdat_base%in%list.files(file.path(dir_bam_base))){
      message(paste("Found",nm_rdat_base,"in directory:",dir_bam_base))
      rdat_base <- dget(file.path(dir_bam_base,nm_rdat_base))
    }else{
      message(paste("Did not find",nm_rdat_base,"in directory:",dir_bam_base))
    }
  }
  ## Identify objects from bam base output
  comp.mats <- rdat_base$comp.mats
  t.series <- rdat_base$t.series
  styr <- rdat_base$parms$styr
  endyr <- rdat_base$parms$endyr
  yr <- styr:endyr

  if(!is.null(repro_sim)){
    P_repro <- repro_sim$P_repro
    nagecfishage <- repro_sim$nagecfishage
    plot_fits_repro <- repro_sim$plot_fits
  }

  ###############################################################
  ## Conduct Monte Carlo draws of parameters and bootstrap data #
  ###############################################################
  agec_nm <- names(init)[grepl(pattern="^obs_agec",names(init))]

    ###  Monte Carlo ###
  if(run_mc){
    # sclim
    sclim_user <- sclim
    sclim_default <- list(M_constant=sclim_gen,
                          steep=sclim_gen,
                          rec_sigma=sclim_gen,
                          Linf=sclim_gen,
                          K=sclim_gen,
                          t0=sclim_gen,
                          Dmort=sclim_gen,
                          Pfa=sclim_gen,
                          Pfb=sclim_gen,
                          Pfma=sclim_gen,
                          Pfmb=sclim_gen,
                          fecpar_a= sclim_gen,
                          fecpar_b= sclim_gen,
                          fecpar_c= sclim_gen,
                          fecpar_batches_sc=sclim_gen
    )

    sclim <- modifyList(sclim_default,sclim_user)

    # fn_par
    fn_par_user <- fn_par
    # Set some defaults which will effectively repeat the base values for parameters
    # if nothing else is supplied by the user in the arguments.
    fn_par_default = list(
      M_constant = "rep(M_constant,nsim)",
      K = "rep(K,nsim)",
      Linf = "rep(Linf,nsim)",
      t0 = "rep(t0,nsim)",
      steep = "rep(steep,nsim)",
      rec_sigma = "rep(rec_sigma,nsim)",
      SR_switch = "rep(SR_switch,nsim)",
      Dmort = "apply(Dmort,2,function(x){rep(x,nsim)})",
      Pfa=     "rep(Pfa,nsim)",
      Pfb=     "rep(Pfb,nsim)",
      Pfma=     "rep(Pfma,nsim)",
      Pfmb=     "rep(Pfmb,nsim)",
      fecpar_a= "rep(fecpar_a,nsim)",
      fecpar_b= "rep(fecpar_b,nsim)",
      fecpar_c= "rep(fecpar_b,nsim)",
      fecpar_batches_sc= "rep(1,nsim)"
    )
    fn_par <- modifyList(fn_par_default, fn_par_user)

    # Check for values
    M_constant_is <- "set_M_constant"%in%names(init)
    Dmort_is <- length(grep("^set_Dmort",names(init),value=TRUE))>0
    fecpar_a_is <- "fecpar_a"%in%names(init)
    fecpar_b_is <- "fecpar_b"%in%names(init)
    fecpar_c_is <- "fecpar_c"%in%names(init)
    fecpar_batches_is <- "fecpar_batches"%in%names(init)

    # get base parameter values
    agebins <- as.numeric(init$agebins)

    Linf <- as.numeric(init$set_Linf[1])
    K <- as.numeric(init$set_K[1])
    t0 <- as.numeric(init$set_t0[1])

    # set defaults
    Pfa <- Pfb <- Pfma <- Pfmb <- NA
    fecpar_a <- if(fecpar_a_is){as.numeric(init$fecpar_a)}else{NA}
    fecpar_b <- if(fecpar_b_is){as.numeric(init$fecpar_b)}else{NA}
    fecpar_c <- if(fecpar_c_is){as.numeric(init$fecpar_c)}else{NA}
    fecpar_batches <- if(fecpar_batches_is){
      a <- init$fecpar_batches
      class(a) <- "numeric"
      a
    }else{NA}


    if(M_constant_is){
      M_constant <- as.numeric(init$set_M_constant[1])
    }else{
      M_constant <- NA
      message("M_constant not found in names(init)")
    }
    if("set_M"%in%names(init)){
      Ma <- as.numeric(init[["set_M"]])
    }else{
      message("set_M not found in init object. getting rdat_base$a.series$M")
      Ma <- setNames(rdat_base$a.series$M,rdat_base$a.series$age)
    }
    steep <- as.numeric(init$set_steep[1])
    rec_sigma <- as.numeric(init$set_rec_sigma[1])
    SR_switch <- init$SR_switch


    if(Dmort_is){
      Dmort <- local({
        a <- init[grep("^set_Dmort",names(init),value=TRUE)] # list
        b <- lapply(a,as.numeric)
        array(unlist(b),dim=c(1,length(b)),dimnames = list(NULL,names(b)))
      })
    }else{
      message("Discard mortality rates not found in the base model (i.e. no names(init) beginning with set_Dmort).")
    }

    ### Simulate reproductive input (repro_sim) ###

    # compute age comps in numbers of fish for repro_sim
    nfish_agec_nm <- names(init)[grepl("^nfish_agec",names(init))]
    agec <-  init[agec_nm] # agec in proportions
    nfish_agec <- init[nfish_agec_nm]
    nagec <- lapply(seq_along(agec),function(j){
      nmj <- names(agec)[j]
      abbj <- gsub("^obs_agec_","",nmj)
      agecj <- agec[[j]]
      # convert to numeric and retain attributes
      att <- attributes(agecj)
      agecj <- apply(agecj,2,as.numeric)
      attributes(agecj) <- att
      nfish_agecj <- nfish_agec[[paste0("nfish_agec_",abbj)]]
      class(nfish_agecj) <- "numeric"
      if(!all(rownames(agecj)==names(nfish_agecj))){
        warning(paste("For",nmi,"not all rownames (years) in agecj match names of nfish_agecj"))
      }
      round(agecj*nfish_agecj)
    })
    names(nagec) <- names(agec)

    # pool age comps in number to a single distribution
    nagec_pool <-
      if(length(nagec)>0){
        colSums(comp_combine(nagec,scale_rows = FALSE))
      }else{
        warning(paste("No age comps found in model. Assuming",nagecfishage, "age samples per age class for estimating uncertainty in reproductive estimates."))
        setNames(rep(nagecfishage,length(init$obs_prop_f)),names(init$obs_prop_f))
      }

    # General repro stuff (should work whether or not repro_sim is TRUE)
    # Proportion female
    Pf <- local({
      a <- init$obs_prop_f
      if(is.matrix(a)){
        warning(paste("obs_prop_f is a matrix. The last row will be used and repeated to fill the matrix in simulations."))
        a <- setNames(as.numeric(tail(a,1)),colnames(tail(a,1)))
      }else{
        a
      }
      class(a) <- "numeric"
      a
    })
    age_Pf <- as.numeric(names(Pf))

    Pfm <- local({
      a <- init$obs_maturity_f
      if(is.matrix(a)){
        warning(paste("obs_maturity_f is a matrix. The last row will be used and repeated to fill the matrix in simulations."))
        a <- setNames(as.numeric(tail(a,1)),colnames(tail(a,1)))
      }else{
        a
      }
      class(a) <- "numeric"
      a
    })

    age_Pfm <- as.numeric(names(Pfm))

    ### more repro_sim stuff
    if(!is.null(repro_sim)){
      nrepro <- local({
        agec_agemax <- max(as.numeric(names(nagec_pool)))
        a <- round(nagec_pool*P_repro) # estimated number of repro samples

        # If any ages in Pf are older than in agec, spread out the n fish from the
        # oldest age of agec among that age and the older ages
        if(any(as.numeric(names(Pf))>agec_agemax)){
          ages2add <- names(Pf)[which(as.numeric(names(Pf))>=agec_agemax)]
          Nage_est <- setNames(bamExtras::exp_decay(age=as.numeric(names(Pf)),Z=as.numeric(init[["set_M"]]),N1=1),names(Pf))
          nagec_plus <- tail(a,1)
          Page_plus <- Nage_est[ages2add]*nagec_plus
          nrepro2add <- round(nagec_plus*(Page_plus/sum(Page_plus)))
          a <- c(a[1:(length(a)-1)],nrepro2add)
        }
        a
      })

      # Check for atypical situations
      Pf_constant <- ifelse(length(unique(Pf))==1,TRUE,FALSE)
      Pf_binary <-   ifelse((all(Pf%in%c(0,1))),TRUE,FALSE)
      if(Pf_constant){ # If all values are the same
        Pfa <- unique(Pf)
        Pfb <- NA
        message(paste("all values of Pf are equal to",Pfa,". Will not attempt to estimate logistic function."))
      }else if(Pf_binary){ # If all values are either 0 or 1
        Pfa <- NA
        Pfb <- mean(c(max(age_Pf[which(Pf==0)]),min(age_Pf[which(Pf==1)])))
        message(paste("all values of Pf are either 0 or 1. Will not attempt to estimate logistic function. a50 is approximately",Pfb))
      }

      # Estimate logistic relationship if possible
      if(!any(c(Pf_constant,Pf_binary))){

        data_Pf <- local({
          y1 <- round(nrepro*Pf)
          y0 <- nrepro-y1
          x <- age_Pf
          data.frame(x=c(rep(x,y0),rep(x,y1)), y=c(rep(0,sum(y0)),rep(1,sum(y1))))
        })
        # logistic fit to prop_f
        fit_Pf <- glm(y~x,data=data_Pf,family="binomial")
        coef_Pf <- setNames(coef(fit_Pf),c("intercept","slope"))
        intercept_Pf <- coef_Pf[[1]]
        Pfa <- slope_Pf <- coef_Pf[[2]]
        Pfb <- a50_Pf <- max(-1e6,as.numeric(coef_Pf[[1]]/-coef_Pf[[2]])) # don't let it be -Inf, as when Pf is a constant 0.5
        Pf_pr <- lgs(x=age_Pf,a=Pfa,b=Pfb)

        # Plot observed values and prediction
        if(plot_fits_repro){
          plot(age_Pf,Pf,xlab="age",ylab="proportion female",ylim=c(0,1))
          lines(age_Pf,Pf_pr)
        }
      }

      # NOTE: Code to make random draws of parameters from multivariate normal distribution
      # mvtnorm::rmvnorm(n=nsim, mean=coef_Pf, sigma=vcov(fit_Pf), method="chol")

      # Female maturity
      nf <- round(nrepro*Pf)  # estimated number of females

      # Check for atypical situations
      Pfm_constant <- ifelse(length(unique(Pfm))==1,TRUE,FALSE)
      Pfm_binary <-   ifelse((all(Pfm%in%c(0,1))),TRUE,FALSE)
      if(Pfm_constant){ # If all values are the same
        Pfma <- unique(Pfm)
        Pfmb <- NA
        warning(paste("all values of Pfm are equal to",Pfma,". Will not attempt to estimate logistic function."))
      }else if(Pfm_binary){ # If all values are either 0 or 1
        Pfma <- NA
        Pfmb <- mean(c(max(age_Pfm[which(Pfm==0)]),min(age_Pfm[which(Pfm==1)])))
        warning(paste("all values of Pfm are either 0 or 1. Will not attempt to estimate logistic function. a50 is approximately",Pfmb))
      }

      # Estimate logistic relationship if possible
      if(!any(c(Pfm_constant,Pfm_binary))){
        data_Pfm <- local({
          y1 <- round(nf*Pfm)
          y0 <- nf-y1
          x <- age_Pfm
          data.frame(x=c(rep(x,y0),rep(x,y1)), y=c(rep(0,sum(y0)),rep(1,sum(y1))))
        })
        # logistic fit to prop_f
        fit_Pfm <- glm(y~x,data=data_Pfm,family="binomial")
        coef_Pfm <- setNames(coef(fit_Pfm),c("intercept","slope"))
        intercept_Pfm <- coef_Pfm[[1]]
        Pfma <- coef_Pfm[[2]]
        Pfmb <- as.numeric(coef_Pfm[[1]]/-coef_Pfm[[2]])
        Pfm_pr <- lgs(x=age_Pfm,a=Pfma,b=Pfmb)

        # Plot observed values and prediction
        if(plot_fits_repro){
          plot(age_Pfm,Pfm,xlab="age",ylab="proportion female mature",ylim=c(0,1))
          lines(age_Pfm,Pfm_pr)
        }
      }
    }else{ # end of repro_sim stuff
      message("Reproductive traits are not being simulated since repro_sim = NULL.")
    }

    # compute limits of parameter values
    Linf_lim  <- Linf*sclim$Linf
    K_lim     <- K*sclim$K
    t0_lim    <- t0*sclim$t0
    M_lim     <- M_constant*sclim$M_constant
    steep_lim <- steep*sclim$steep
    Pfa_lim   <- Pfa*sclim$Pfa
    Pfb_lim   <- Pfb*sclim$Pfb
    Pfma_lim   <- Pfma*sclim$Pfma
    Pfmb_lim   <- Pfmb*sclim$Pfmb
    fecpar_a_lim <- fecpar_a*sclim$fecpar_a
    fecpar_b_lim <- fecpar_b*sclim$fecpar_b
    fecpar_c_lim <- fecpar_c*sclim$fecpar_c
    fecpar_batches_sc_lim <- sclim$fecpar_batches_sc

    Linf_min  <- min(Linf_lim)
    K_min     <- min(K_lim)
    t0_min    <- min(t0_lim)
    M_min     <- min(M_lim)
    steep_min <- min(steep_lim)
    Pfa_min   <- min(Pfa_lim)
    Pfb_min   <- min(Pfb_lim)
    Pfma_min  <- min(Pfma_lim)
    Pfmb_min  <- min(Pfmb_lim)
    fecpar_a_min  <- min(fecpar_a_lim)
    fecpar_b_min  <- min(fecpar_b_lim)
    fecpar_c_min  <- min(fecpar_c_lim)
    fecpar_batches_sc_min  <- min(fecpar_batches_sc_lim)

    Linf_max  <- max(Linf_lim)
    K_max     <- max(K_lim)
    t0_max    <- max(t0_lim)
    M_max     <- max(M_lim)
    steep_max <- max(steep_lim)
    Pfa_max   <- max(Pfa_lim)
    Pfb_max   <- max(Pfb_lim)
    Pfma_max  <- max(Pfma_lim)
    Pfmb_max  <- max(Pfmb_lim)
    fecpar_a_max  <- max(fecpar_a_lim)
    fecpar_b_max  <- max(fecpar_b_lim)
    fecpar_c_max  <- max(fecpar_c_lim)
    fecpar_batches_sc_max  <- max(fecpar_batches_sc_lim)

    # NOTE column names of Dmort_lim must match parameter names in init object
    if(Dmort_is){Dmort_lim <- Dmort[c(1,1),,drop=FALSE] * sclim$Dmort}
    ## Conduct Monte Carlo draws

    # Vectors (sim)
    Linf_sim <- eval(parse(text=fn_par$Linf))
    K_sim <- eval(parse(text=fn_par$K))
    t0_sim <- eval(parse(text=fn_par$t0))
    if(M_constant_is){
      M_sim <- eval(parse(text=fn_par$M_constant))
      Masc_sim <- M_sim/M_constant # Multiply by M-at-age to rescale appropriately with sim M_constant
    }else{
      M_sim <- rep(NA,nsim)
      Masc_sim <- runif(nsim,sclim$M_constant[1],sclim$M_constant[2])
    }
    rec_sigma_sim <- eval(parse(text=fn_par$rec_sigma))
    steep_sim <- eval(parse(text=fn_par$steep))
    SR_switch_sim <- eval(parse(text=fn_par$SR_switch))

    if(Dmort_is){
      Dmort_sim <- eval(parse(text=fn_par$Dmort))
      # Assure matrix structure
      Dmort_sim <- matrix(Dmort_sim,nrow=nsim,dimnames=list(NULL,dimnames(Dmort_sim)[[2]]))
    }

    # Initialize
    Pfa_sim <- Pfb_sim <- Pfma_sim <- Pfmb_sim <- rep(NA,nsim)

    # These should have actual values if they are in the assessment
    fecpar_a_sim <- rep(fecpar_a,nsim)
    fecpar_b_sim <- rep(fecpar_b,nsim)
    fecpar_c_sim <- rep(fecpar_c,nsim)
    Pf_sim <-  matrix(Pf,nrow=nsim,ncol=length(age_Pf),dimnames=list(sim=nm_sim,age=age_Pf),byrow = TRUE)
    Pfm_sim <- matrix(Pfm,nrow=nsim,ncol=length(age_Pfm),dimnames=list(sim=nm_sim,age=age_Pfm),byrow = TRUE)

    if(!is.null(repro_sim)){
      # Proportion female (sex ratio)
      Pfa_sim <- if(!Pf_binary){eval(parse(text=fn_par$Pfa))}else{
        rep(NA,nsim)
      }
      Pfb_sim <- if(!Pf_constant){eval(parse(text=fn_par$Pfb))}else{
        rep(NA,nsim)
      }

      Pf_sim <- local({
        if(Pf_constant){
          fn <- expression(rep(Pfa_sim[i],length(age_Pf)))
        }else if(Pf_binary){
          fn <- expression((age_Pf>=Pfb_sim[i])*1)
        }else{
          fn <- expression(bamExtras::lgs(x=age_Pf,a=Pfa_sim[i],b=Pfb_sim[i],type = "slope_inflection"))
        }
        a <- lapply(1:nsim,function(i){
          round(pmin(1,pmax(0,eval(fn))),ndigits)
        })
        b <- do.call(rbind,a)
        dimnames(b) <- list(sim=nm_sim,age=age_Pf)
        return(b)
      })

      # Proportion of females mature (maturity)
      Pfma_sim <- if(!Pfm_binary){eval(parse(text=fn_par$Pfma))}else{
        rep(NA,nsim)
      }
      Pfmb_sim <- if(!Pfm_constant){eval(parse(text=fn_par$Pfmb))}else{
        rep(NA,nsim)
      }

      Pfm_sim <- local({
        if(Pfm_constant){
          fn <- expression(rep(Pfma_sim[i],length(age_Pfm)))
        }else if(Pfm_binary){
          fn <- expression((age_Pfm>=Pfmb_sim[i])*1)
        }else{
          fn <- expression(bamExtras::lgs(x=age_Pfm,a=Pfma_sim[i],b=Pfmb_sim[i],type = "slope_inflection"))
        }
        a <- lapply(1:nsim,function(i){
          round(pmin(1,pmax(0,eval(fn))),ndigits)
        })
        b <- do.call(rbind,a)
        dimnames(b) <- list(sim=nm_sim,age=age_Pfm)
        return(b)
      })

      fecpar_a_sim <- if(fecpar_a_is){
        eval(parse(text=fn_par$fecpar_a))
      }else{
        rep(NA,nsim)
      }
      fecpar_b_sim <- if(fecpar_b_is){
        eval(parse(text=fn_par$fecpar_b))
      }else{
        rep(NA,nsim)
      }
      fecpar_c_sim <- if(fecpar_c_is){
        eval(parse(text=fn_par$fecpar_c))
      }else{
        rep(NA,nsim)
      }
    } # end if(!is.null(repro_sim)){

    # This should always work
    if(fecpar_batches_is){
      fecpar_batches_sc_sim <- eval(parse(text=fn_par$fecpar_batches_sc))
      if(length(fecpar_batches)==1){
        fecpar_batches_sim <- fecpar_batches_sc_sim*fecpar_batches
      }else{
        fecpar_batches_sim <- round(t(matrix(rep(fecpar_batches,nsim),ncol=nsim,dimnames=list(age=agebins,sim=nm_sim)))*fecpar_batches_sc_sim,ndigits)
      }
    }else{
      fecpar_batches_sc_sim <- rep(NA,nsim)
    }

    # Custom fn_par (2026-01-15 NPK I don't think this is really in use yet. This should work and will be necessary for atypical parameters such as Linf_L and K_L)
    fn_par_custom <- fn_par_user[!names(fn_par_user)%in%names(fn_par_default)]
    par_sim_custom <- list()
    if(length(fn_par_custom)>0){
      for(i in seq_along(fn_par_custom)){
        par_sim_custom[[i]] <- setNames(eval(parse(text=fn_par_custom[[i]])),nm_sim)
      }
      names(par_sim_custom) <- names(fn_par_custom)
    }

    ## Matrices (sim,age)
    # collected par values
    par_sim <- local({
      a <- data.frame(Linf=Linf_sim,
                      K=K_sim,
                      t0=t0_sim,
                      M_constant=M_sim,
                      Masc=Masc_sim,
                      steep=steep_sim,
                      rec_sigma=rec_sigma_sim,
                      SR_switch=SR_switch_sim,
                      Pfa=Pfa_sim,
                      Pfb=Pfb_sim,
                      Pfma=Pfma_sim,
                      Pfmb=Pfmb_sim,
                      fecpar_a=fecpar_a_sim,
                      fecpar_b=fecpar_b_sim,
                      fecpar_c=fecpar_c_sim,
                      fecpar_batches_sc=fecpar_batches_sc_sim
      )
      if(Dmort_is) {b <- cbind(a,Dmort_sim)
      }else{
        b <- a
      }
      c <- as.data.frame(lapply(b, function(x){if(is.numeric(x)){round(x,ndigits)}else{x}}))
      rownames(c) <- nm_sim
      return(c)
    })
    # Natural mortality
    Ma_sim <- round(t(matrix(rep(Ma,nsim),ncol=nsim,dimnames=list(age=agebins,sim=nm_sim)))*Masc_sim,ndigits)

    # Check SR_switch settings
    SR_switch_check <- sum(!paste(par_sim$SR_switch)%in%c("1","2")) # Number of SR_switch settings not 1 or 2 (Beverton Holt or Ricker)
    if(SR_switch_check>0){message(paste(SR_switch_check,"simulated values of SR_switch were not 1 or 2 (Beverton-Holt or Ricker). Steepness is probably not being used in model runs."))}


  } # end if(run_mc)

  # Bootstrap data
  if(run_boot){
  #%%%% bootstrap setup %%%%#
  #%% Indices %%#
  obs_cpue_nm <- names(init)[grepl(pattern="obs_cpue_",names(init))]
  cv_cpue_nm <- names(init)[grepl(pattern="cv_cpue_",names(init))]
  cpue_root <- gsub("obs_cpue_","",obs_cpue_nm)

  # if(!is.null(data_sim$cv_cpue)){
  #   cv_cpue <-
  # }

  sim_cpue <- list()
  for(nm_i in obs_cpue_nm){
    cpue_root_i <- gsub("obs_cpue_","",nm_i)
    cv_cpue_nm_i <- gsub("obs_cpue_","obs_cv_cpue_",nm_i)

    yrs_cpue_nm_i <- paste0("yrs_cpue_",cpue_root_i)
    styr_cpue_nm_i <- gsub("obs","styr",nm_i)
    endyr_cpue_nm_i <- gsub("obs","endyr",nm_i)
    if(yrs_cpue_nm_i%in%names(init)){ # If the specific years of the index are given, use them. Otherwise look for range of years given and use those
      yrs_cpue_i <- init[[yrs_cpue_nm_i]]
    }else{
      styr_cpue_i  <- init[[styr_cpue_nm_i]]
      endyr_cpue_i <- init[[endyr_cpue_nm_i]]
      yrs_cpue_i <- paste(styr_cpue_i:endyr_cpue_i)
    }

    cpue_i <- setNames(as.numeric(init[[nm_i]]),yrs_cpue_i)

    # Set default cvs
    cv_cpue_i <- local({
      a <- if(cv_cpue_nm_i%in%names(init)){
      as.numeric(init[[cv_cpue_nm_i]])
      }else{
        message(paste0(cv_cpue_nm_i," not found in names(init). Using cvs of par_default$cv_U = ",par_default$cv_U, " to simulate ",nm_i))
        rep(par_default$cv_U,length(cpue_i))
      }
      setNames(a,yrs_cpue_i)
    })

    # Override and replace cvs if function or values are provided in data_sim
    # get data_sim$cv_U
    xi <- data_sim$cv_U
    if(is.function(xi)){ # If xi is a function, apply it
      cv_cpue_i <- xi(cv_cpue_i)
    }else if(is.data.frame(xi)){
      if(cpue_root_i%in%dimnames(xi)[[2]]){
        if(any(!is.na(xi[,cpue_root_i]))){
          message(paste("Using values from data_sim for",nm_i))
          ci <- local({
            ai <- xi[,cpue_root_i,drop=FALSE]
            ai[complete.cases(ai),,drop=FALSE]
          })
          if(!all(rownames(ci)%in%yrs_cpue_i)){
            message(paste("Some rownames (years) in data_sim$cv_U for",cpue_root_i,"do not match values in",nm_i,". These values will not be used."))
          }
          cv_cpue_i <- setNames(ci[yrs_cpue_i,cpue_root_i],yrs_cpue_i)
          if(any(yrs_cpue_i != names(cv_cpue_i))){warning(paste("years of",paste0("data_sim$cv_U$",cpue_root_i),"do not match years of",nm_i,"in the base model"))}
        }
      }else{
        message(paste(nm_i,"found in model but",cpue_root_i,"not found in data_sim$cv_U"))
      }
    }else{
      message("data_sim$cv_U is neither a function nor a data frame and will not be applied")
    }
    rm(xi)

    if(any(c("U","cpue")%in%data_type_resamp)){

      message(paste("bootstrapping",nm_i))
      # look for negative values in observed data (e.g. missing values in bam dat coded as -999)
      # and replace them with NA for bootstrapping
      cpue_i2 <- cpue_i
      cv_cpue_i2 <- cv_cpue_i
      neg_cpue_ij <- which(cpue_i<0)
      if(length(neg_cpue_ij)>0){
        message(paste("Negative values found in",nm_i,"will be ignored when bootstrapping."))
      }

      neg_cv_cpue_ij <- which(cv_cpue_i<0)
      if(length(neg_cv_cpue_ij)>0){
        message(paste("Negative values found in",cv_cpue_nm_i,"will be ignored when bootstrapping."))
      }

      cpue_i2[neg_cpue_ij] <- NA
      cv_cpue_i2[neg_cv_cpue_ij] <- NA

      sim_cpue_i <- lnorm_vector_boot(cpue_i2,cv_cpue_i2,nsim,
                                      standardize=TRUE,digits=ndigits)
      sim_cpue_i[neg_cpue_ij,] <- cpue_i[neg_cpue_ij] # replace NA with original negative values to run bam

    }else{
      sim_cpue_i <- matrix(cpue_i, nrow=length(yrs_cpue_i), ncol=nsim)
    }
    dimnames(sim_cpue_i) <- list("year"=yrs_cpue_i,"sim"=nm_sim)
    sim_cpue[[cpue_root_i]] <- sim_cpue_i
  }

  #%% Landings %%#
  obs_L_nm <- names(init)[grepl(pattern="obs_L_",names(init))]
  sim_L <- list()
  cv_L_nm <- names(init)[grepl(pattern="cv_L_",names(init))]

  for(nm_i in obs_L_nm){
    L_root_i <- gsub("obs_L_","",nm_i)
    cv_L_nm_i <- gsub("obs_L_","obs_cv_L_",nm_i)

    yrs_L_nm_i <- paste0("yrs_L_",L_root_i)
    styr_L_nm_i <- gsub("obs","styr",nm_i)
    endyr_L_nm_i <- gsub("obs","endyr",nm_i)
    if(yrs_L_nm_i%in%names(init)){ # If the specific years are given, use them. Otherwise look for range of years given and use those.
      yrs_L_i <- init[[yrs_L_nm_i]]
    }else{
      styr_L_i  <- init[[styr_L_nm_i]]
      endyr_L_i <- init[[endyr_L_nm_i]]
      yrs_L_i <- paste(styr_L_i:endyr_L_i)
    }

    L_i <- as.numeric(init[[nm_i]])

    # Set default cvs
    cv_L_i <- if(cv_L_nm_i%in%names(init)){
      as.numeric(init[[cv_L_nm_i]])}else{
        message(paste0(cv_L_nm_i," not found in names(init). Using cvs of par_default$cv_L = ",par_default$cv_L, " to simulate ",nm_i))
        rep(par_default$cv_L,length(L_i))
      }
    # Override and replace cvs if function or values are provided in data_sim
    # get data_sim$cv_L
    xi <- data_sim$cv_L
    if(is.function(xi)){ # If xi is a function, apply it
      cv_L_i <- xi(cv_L_i)
    }else if(is.data.frame(xi)){
      if(L_root_i%in%dimnames(xi)[[2]]){
        if(any(!is.na(xi[,L_root_i]))){
          message(paste("Using values from data_sim for",nm_i))
          ci <- local({
            ai <- xi[,L_root_i,drop=FALSE]
            ai[complete.cases(ai),,drop=FALSE]
          })
          if(!all(rownames(ci)%in%yrs_L_i)){
            message(paste("Some rownames (years) in data_sim$cv_L for",L_root_i,"do not match values in",nm_i,". These values will not be used."))
          }
          cv_L_i <- setNames(ci[yrs_L_i,L_root_i],yrs_L_i)
          if(any(yrs_L_i != names(cv_L_i))){warning(paste("years of",paste0("data_sim$cv_L$",L_root_i),"do not match years of",nm_i,"in the base model"))}
        }
      }else{
        message(paste(nm_i,"found in model but",L_root_i,"not found in data_sim$cv_L"))
      }
    }else{
      message("data_sim$cv_L is neither a function nor a data frame and will not be applied")
    }
    rm(xi)


    if("L"%in%data_type_resamp){
      message(paste("bootstrapping",nm_i))
      sim_L_i <- lnorm_vector_boot(L_i,cv_L_i,nsim,
                                   standardize=FALSE,digits=ndigits)
    }else{
      sim_L_i <- matrix(L_i,nrow=length(yrs_L_i),ncol=nsim)
    }
    dimnames(sim_L_i) <- list("year"=yrs_L_i,"sim"=nm_sim)
    sim_L[[L_root_i]] <- sim_L_i
  }

  #%% Discards (released) %%#
  obs_released_nm <- names(init)[grepl(pattern="obs_released_",names(init))]
  sim_D <- list()

  if(length(obs_released_nm)>0){
    cv_D_nm <- names(init)[grepl(pattern="obs_cv_D_",names(init))]

    for(nm_i in obs_released_nm){
      D_root_i <- gsub("obs_released_","",nm_i)
      cv_D_nm_i <- gsub("obs_released_","obs_cv_D_",nm_i)

      yrs_D_nm_i <- paste0("yrs_D_",D_root_i)
      styr_D_nm_i <- gsub("obs_released","styr_D",nm_i)
      endyr_D_nm_i <- gsub("obs_released","endyr_D",nm_i)
      if(yrs_D_nm_i%in%names(init)){ # If the specific years are given, use them. Otherwise look for range of years given and use those
        yrs_D_i <- init[[yrs_D_nm_i]]
      }else{
        styr_D_i  <- init[[styr_D_nm_i]]
        endyr_D_i <- init[[endyr_D_nm_i]]
        yrs_D_i <- paste(styr_D_i:endyr_D_i)
      }

      D_i <- as.numeric(init[[nm_i]])

      # Set default cvs
      cv_D_i <- if(cv_D_nm_i%in%names(init)){
        as.numeric(init[[cv_D_nm_i]])}else{
          message(paste0(cv_D_nm_i," not found in names(init). Using cvs of par_default$cv_D = ",par_default$cv_D, " to simulate ",nm_i))
          rep(par_default$cv_D,length(D_i))
        }

      # Override and replace cvs if function or values are provided in data_sim
      # get data_sim$cv_D
      xi <- data_sim$cv_D
      if(is.function(xi)){ # If xi is a function, apply it
        cv_D_i <- xi(cv_D_i)
      }else if(is.data.frame(xi)){
        if(D_root_i%in%dimnames(xi)[[2]]){
          if(any(!is.na(xi[,D_root_i]))){
            message(paste("Using values from data_sim for",nm_i))
            ci <- local({
              ai <- xi[,D_root_i,drop=FALSE]
              ai[complete.cases(ai),,drop=FALSE]
            })
            if(!all(rownames(ci)%in%yrs_D_i)){
              message(paste("Some rownames (years) in data_sim$cv_D for",D_root_i,"do not match values in",nm_i,". These values will not be used."))
            }
            cv_D_i <- setNames(ci[yrs_D_i,D_root_i],yrs_D_i)
            if(any(yrs_D_i != names(cv_D_i))){warning(paste("years of",paste0("data_sim$cv_D$",D_root_i),"do not match years of",nm_i,"in the base model"))}
          }
        }else{
          message(paste(nm_i,"found in model but",D_root_i,"not found in data_sim$cv_D"))
        }
      }else{
        message("data_sim$cv_D is neither a function nor a data frame and will not be applied")
      }
      rm(xi)


      if("D"%in%data_type_resamp){
        message(paste("bootstrapping",nm_i))
        sim_D_i <- lnorm_vector_boot(D_i,cv_D_i,nsim,
                                     standardize=FALSE,digits=ndigits)
      }else{
        sim_D_i <- matrix(D_i,nrow=length(yrs_D_i),ncol=nsim)
      }
      dimnames(sim_D_i) <- list("year"=yrs_D_i,"sim"=nm_sim)
      sim_D[[D_root_i]] <- sim_D_i
    }
  }else{
    message("No discard time series found in the base model (i.e. no names(init) beginning with obs_released).")
  }

  #%% Age composition %%#
  # agec_nm <- names(init)[grepl(pattern="obs_agec",names(init))]
  sim_acomp <- list()

  if(length(agec_nm)>0){
    agec_root <- gsub("^obs_agec_","",agec_nm)

    for(agec_root_i in agec_root){
      acomp_ob_nm_i <- gsub("\\_","\\.",paste0("acomp.",agec_root_i,".ob"))
      x_i <- comp.mats_i <- comp.mats[[acomp_ob_nm_i]]

      agebins_acomp_i <- factor(colnames(x_i),levels=colnames(x_i))

      yrs_i <- as.numeric(rownames(x_i))

      # Get nfish from init instead of rdat_base. If years of comps are excluded from
      # fitting by bam based on minSS then values of -99999 will appear for
      # values of nfish in the rdat. This won't affect results, but the -99999
      # values cause an error when trying to simulate comps.
      nfish_i <- setNames(as.numeric(init[[paste0("nfish_agec_",agec_root_i)]]),paste(yrs_i))
      dimnames(x_i) <- list("year"=yrs_i,"age"=agebins_acomp_i)

      x_i <- as.data.frame.matrix(x_i)
      x_i$nfish <- nfish_i

      if(any(c("age","agec")%in%data_type_resamp)){
        message(paste("bootstrapping",paste0("obs_agec_",agec_root_i)))
        sim_acomp_i <- foreach::foreach(i=1:nsim) %dopar% {
          out <- apply(x_i,1,function(y){
            if(sum(y[names(y)!="nfish"])>0&y["nfish"]>0){
              ages_y <- sample(x=agebins_acomp_i, size=y["nfish"], replace=TRUE, prob=y[paste(agebins_acomp_i)])
              round(table(factor(ages_y,agebins_acomp_i))/y["nfish"],4)
            }else{
              ages_y <- rep(0,length(agebins_acomp_i))
            }
          })
          t(out)
        }
        sim_acomp_i <- array(unlist(sim_acomp_i),dim=c(dim(comp.mats_i),length(sim_acomp_i)),
                             dimnames=list(year=yrs_i,age=levels(agebins_acomp_i),sim=nm_sim))

        sim_acomp[[agec_root_i]] <- sim_acomp_i
      }else{
        sim_acomp[[agec_root_i]] <- array(unlist(comp.mats_i),dim=c(dim(comp.mats_i),nsim),
                                          dimnames=list(year=yrs_i,age=levels(agebins_acomp_i),sim=nm_sim))
      }
    }
  }else{
    message("No age compositions found in the base model (i.e. no names(init) beginning with obs_agec).")
  }

  #%% Length composition %%#
  lenc_nm <- names(init)[grepl(pattern="obs_lenc",names(init))]
  sim_lcomp <- list()

  if(length(lenc_nm)>0){
    lenc_root <- gsub("^obs_lenc_","",lenc_nm)

    for(lenc_root_i in lenc_root){
      lcomp_ob_nm_i <- gsub("\\_","\\.",paste0("lcomp.",lenc_root_i,".ob"))
      x_i <- comp.mats_i <- comp.mats[[lcomp_ob_nm_i]]

      lenbins_lcomp_i <- factor(colnames(x_i),levels=colnames(x_i))

      yrs_i <- as.numeric(rownames(x_i))

      # Get nfish from init instead of rdat_base. If years of comps are excluded from
      # fitting by bam based on minSS then values of -99999 will appear for
      # values of nfish in the rdat. This won't affect results, but the -99999
      # values cause an error when trying to simulate comps.
      nfish_i <- setNames(as.numeric(init[[paste0("nfish_lenc_",lenc_root_i)]]),paste(yrs_i))

      dimnames(x_i) <- list("year"=yrs_i,"len"=lenbins_lcomp_i)

      x_i <- as.data.frame.matrix(x_i)
      x_i$nfish <- nfish_i

      if(any(c("len","lenc")%in%data_type_resamp)){
        message(paste("bootstrapping",paste0("obs_lenc_",lenc_root_i)))
        sim_lcomp_i <- foreach::foreach(i=1:nsim) %dopar% {
          out <- apply(x_i,1,function(y){
            if(sum(y[names(y)!="nfish"])>0&y["nfish"]>0){
              lens_y <- sample(x=lenbins_lcomp_i, size=y["nfish"], replace=TRUE, prob=y[paste(lenbins_lcomp_i)])
              round(table(factor(lens_y,lenbins_lcomp_i))/y["nfish"],4)
            }else{
              lens_y <- rep(0,length(lenbins_lcomp_i))
            }
          })
          t(out)
        }
        sim_lcomp_i <- array(unlist(sim_lcomp_i),dim=c(dim(comp.mats_i),length(sim_lcomp_i)),
                             dimnames=list(year=yrs_i,len=levels(lenbins_lcomp_i),sim=nm_sim))
        sim_lcomp[[lenc_root_i]] <- sim_lcomp_i
      }else{
        sim_lcomp[[lenc_root_i]] <- array(unlist(comp.mats_i),dim=c(dim(comp.mats_i),nsim),
                                          dimnames=list(year=yrs_i,len=levels(lenbins_lcomp_i),sim=nm_sim))
      }
    }
  }else{
    message("No length compositions found in the base model (i.e. no names(init) beginning with obs_lenc).")
  }

} #end  if(run_boot)
  ###########################
  ## Build list of init objects for all sim

  sim_out <- foreach::foreach(i=1:nsim,
                     #, .export = ls()[grepl("xboot.obs",ls())] # Pass additional objects to foreach
                     .packages=c("bamExtras")
  ) %dopar% {
    nm_sim_i <- nm_sim[i]

    init_i <- init # Copy init from base for each bootstrap run, to initialize

    #%%%%  Modify init_i %%%%#

    #%% Create Monte Carlo/bootstrap data set %%#
    # Fix parameters
    init_i[fix_par] <- lapply(init_i[fix_par],function(x){x["phase"] <- paste(-abs(as.numeric(x["phase"]))); x})


    ##### Monte Carlo #####
    if(run_mc){
      # natural mortality
      if(M_constant_is){
        init_i$set_M_constant[c(1,5)] <- paste(par_sim$M_constant[i])
      }
      if("set_M"%in%names(init)){
        init_i$set_M[] <- paste(Ma_sim[i,])
      }

      # Dmort
      if(Dmort_is){
        for(j in colnames(Dmort)){
          set_Dmort_ij <- par_sim[i,j]
          init_i[[j]] <- paste(set_Dmort_ij)
        }
      }

      # K
      init_i$set_K[c(1,5)] <- paste(par_sim$K[i])

      # Linf
      init_i$set_Linf[c(1,5)] <- paste(par_sim$Linf[i])

      # t0
      init_i$set_t0[c(1,5)] <- paste(par_sim$t0[i])

      # steepness
      init_i$set_steep[c(1,5)] <- paste(par_sim$steep[i])

      # rec_sigma
      init_i$set_rec_sigma[c(1,5)] <- paste(par_sim$rec_sigma[i])

      # SR_switch
      init_i$SR_switch <- paste(par_sim$SR_switch[i])

      # obs_prop_f
      init_i$obs_prop_f <- Pf_sim[i,]

      # # obs_maturity_f
      # init_i$obs_maturity_f <- Pfm_sim[i,]
      # obs_maturity_f
      if(is.matrix(init_i$obs_maturity_f)){
        a <- init_i$obs_maturity_f
        init_i$obs_maturity_f <- matrix(Pfm_sim[i,],byrow=TRUE,nrow=nrow(a),ncol=ncol(a),dimnames=dimnames(a))
      }else{
        init_i$obs_maturity_f <- Pfm_sim[i,]
      }

      if(fecpar_a_is){init_i$fecpar_a <- fecpar_a_sim[i]}
      if(fecpar_b_is){init_i$fecpar_b <- fecpar_b_sim[i]}
      if(fecpar_c_is){init_i$fecpar_c <- fecpar_c_sim[i]}
      if(fecpar_batches_is){
        init_i$fecpar_batches <- if(is.matrix(fecpar_batches_sim)){
          fecpar_batches_sim[i,]
        }else{
          fecpar_batches_sim[i]
        }
      }

      # Add values from custom parameters
      for(j in seq_along(par_sim_custom)){
        nm_j <- names(par_sim_custom)[j]
        x_j <- init_i[[nm_j]]
        par_ji <- par_sim_custom[[nm_j]][[i]] # New par value
        if(length(x_j)==1){
          init_i[[nm_j]] <- par_ji
        }else if(length(x_j)==7){
          init_i[[nm_j]][c(1,5)] <- par_ji
        }
      }
    } # end if(run_mc)

    ##### Bootstrap #####
    if(run_boot){
      spf <- paste0("%01.",ndigits,"f")
      # Add indices (cpue)
      for(nm_j in obs_cpue_nm){
        cpue_root_j <- gsub("obs_cpue_","",nm_j)
        init_i[[nm_j]][] <- sprintf(spf,sim_cpue[[cpue_root_j]][,i])
      }

      # Add landings (L)
      for(nm_j in obs_L_nm){
        L_root_j <- gsub("obs_L_","",nm_j)
        init_i[[nm_j]][] <- sprintf(spf,sim_L[[L_root_j]][,i])
      }

      # Add discards (D)
      if(length(obs_released_nm)>0){
        for(nm_j in obs_released_nm){
          released_root_j <- gsub("obs_released_","",nm_j)
          init_i[[nm_j]][] <- sprintf(spf,sim_D[[released_root_j]][,i])
        }
      }

      # Add age comps
      if(length(agec_nm)>0){
        for(agec_root_j in agec_root){
          # acomp_ob_nm_j <- paste0("acomp.",agec_root_j,".ob") # rdat naming convention
          agec_nm_j <- gsub("\\.","\\_",paste0("obs_agec_",agec_root_j))    # tpl naming convention
          agec_ij <- sim_acomp[[agec_root_j]][,,i,drop=FALSE]
          dim_agec_ij <- dim(agec_ij)
          agec_ij <- matrix(apply(agec_ij,2,function(x){sprintf(spf,x)}),nrow=dim_agec_ij[1])
          dimnames(agec_ij) <- NULL
          init_i[[agec_nm_j]][] <- agec_ij
        }
      }

      if(length(lenc_nm)>0){
        for(lenc_root_j in lenc_root){
          # lcomp_ob_nm_j <- paste0("lcomp.",lenc_root_j,".ob") # rdat naming convention
          lenc_nm_j <- gsub("\\.","\\_",paste0("obs_lenc_",lenc_root_j))    # tpl naming convention
          lenc_ij <- sim_lcomp[[lenc_root_j]][,,i,drop=FALSE]
          dim_lenc_ij <- dim(lenc_ij)
          lenc_ij <- matrix(apply(lenc_ij,2,function(x){sprintf(spf,x)}),nrow=dim_lenc_ij[1])
          dimnames(lenc_ij) <- NULL
          init_i[[lenc_nm_j]][] <- lenc_ij
        }
      }
    } # end if(run_boot)

    #%% Incorporate changes to init_i back into BAM dat %%#
    if(!run_est&sim_type_return!="init"){
    bam_i <- bam2r(
      dat_obj = dat,
      tpl_obj = tpl,
      cxx_obj = cxx,
      init = init_i)

    out_i <- bam_i[[sim_type_return]]
    }else{
      out_i <- init_i
    }
  invisible(out_i)
}
names(sim_out) <- nm_sim

  ############################
  ## Run the bam base model to build the executable and establish a reference?
  #%%%%%%%%%%%%%%%%%%%%%
  #### Run MCBE (maybe this should be a separate function?) ####
  #%%%%%%%%%%%%%%%%%%%%%
  if(run_est){
    if(dir.exists(dir_bam_sim)){
      if(prompt_me){
        delete_files <- readline(prompt=paste0(paste("The folder",paste0("'",dir_bam_sim,"'"),"already exists. "),"Are you sure you want to delete all files in ",dir_bam_sim,"? (Enter TRUE or FALSE)"))
      }else{
        delete_files <- TRUE
      }
      if(delete_files){
        res <- unlink(paste0(dir_bam_sim,"/*"),recursive=TRUE,force=TRUE)
        txt_res <- ifelse(res!=0,"could not be","were")
        message(paste("The contents of the folder",paste0("'",dir_bam_sim,"'"),txt_res,"deleted."))
        if(dir.exists(dir_bam_sim_fail)){
          unlink(paste0(dir_bam_sim_fail,"/*"),force=TRUE)
          message(paste("The contents of the folder",paste0("'",dir_bam_sim_fail,"'"),"have been deleted."))
        }
        if(length(list.files(dir_bam_sim))>0){
          stop(paste("The contents of the folder",paste0("'",dir_bam_sim,"'"),"were not actually deleted."))
        }

      }else{
        message(paste("If you do not want to delete the contents of ",dir_bam_sim,", please specify a new value of dir_bam_sim and rerun."))
      }
    }else{
      dir.create(dir_bam_sim)
      message(paste("Created the folder",paste0("'",dir_bam_sim,"'.")))
    }

    message(paste("Running MCBE for", nsim,"sims in parallel on",ncores,"cores at",
                  format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

    est_out <- foreach::foreach(i=1:nsim,
                        #, .export = ls()[grepl("xboot.obs",ls())] # Pass additional objects to foreach
                        .packages=c("bamExtras")
    ) %dopar% {
      nm_sim_i <- nm_sim[i]

      init_i <- sim_out[[i]] # Copy init from base for each bootstrap run, to initialize

      #%% Incorporate changes to init_i back into BAM dat %%#
      bam_i <- bam2r(
        dat_obj = dat,
        tpl_obj = tpl,
        cxx_obj = cxx,
        init = init_i)

      #%%  File management stuff
      sim_dir_i <- paste0("sim_",nm_sim_i)
      dir.create(sim_dir_i)
      setwd(sim_dir_i)

      fileName_exe_base <- paste0(fileName,".exe")
      fileName_sim_i <- paste0(nm_sim_i,"-",fileName)

      fileName_dat_i  <- paste0(fileName_sim_i,".dat")  # Name of dat file for i
      fileName_rdat_i <- paste0(fileName_sim_i,".rdat") # Name of rdat file for i
      fileName_exe_i  <- paste0(fileName_sim_i,".exe")  # Name of exe file for i
      fileName_par_i  <- paste0(fileName_sim_i,".par")  # Name of par file for i

      # # Copy executable to directory for sim i
      # file.copy(file.path("..",dir_bam_base,fileName_exe_base), fileName_exe_i, overwrite=TRUE)
      #
      # # Rewrite dat file incorporating modified data
      # writeLines(text=bam_i$dat, con=fileName_dat_i)
      #
      # #######Run bam for sim_i
      # code_i <- shell(paste(fileName_exe_i, admb_options_sim, fileName_dat_i, sep=" "))
      dir_bam_i <- "bam"
      bamout_i <- run_bam(
        fileName=fileName_sim_i,
        dir_bam = dir_bam_i,
        admb_options = admb_options_sim,
        dat_obj=bam_i$dat,
        exe_file=file.path("..",dir_bam_base,fileName_exe_base),
        return_obj = NULL,
        subset_rdat=subset_rdat,
        unlink_dir_bam=FALSE
        )


      par_i <- bamout_i$par
      lk_total_i <- bamout_i$admb$lk_total
      grad_max_i <- bamout_i$admb$grad_max
      error_code_run_i <- bamout_i$admb$error_code_run

      nm_files <- list.files(dir_bam_i)
      is_rdat_file <- fileName_rdat_i%in%nm_files

      if(is_rdat_file&is.finite(as.numeric(lk_total_i))&is.finite(as.numeric(grad_max_i))&error_code_run_i==0){
        rdat_i <- dget(file=file.path(dir_bam_i,fileName_rdat_i))
        # Collect the most important values from rdat
        parms_i <- unlist(rdat_i$parms)
        # Add values from par.sim to parms, unless there is already a parameter with the same name in parms
        if(run_mc){
          rdat_i$parms <- c(parms_i,par_sim[i,which(!names(par_sim)%in%names(parms_i))])
        }

        parm.cons.result_i <- unlist(lapply(rdat_i$parm.cons,function(x){x[8]}))
        like_i <- rdat_i$like
        sdnr_i <- rdat_i$sdnr
        vals_i <- c(parms_i[!names(parms_i)%in%names(parm.cons.result_i)],
                    parm.cons.result_i,
                    like_i,
                    sdnr_i)

        # # Subset rdat to decrease size on disk
        # for(nm_k in names(subset_rdat)){
        #   if(!is.null(subset_rdat[[nm_k]])){ # won't subset anything if the value of nm_k is NULL
        #     k <- rdat_i[[nm_k]]
        #     if(nm_k%in%c("eq.series","pr.series")){
        #       rdat[[nm_k]] <- k[seq(1,nrow(k),length=subset_rdat[[nm_k]]),,drop=FALSE]
        #     }else if(nm_k=="t.series"){
        #       rdat[[nm_k]] <- k[paste(subset_rdat[[nm_k]]),]
        #     }else{
        #       message(paste("no subsetting method is currently specified for",nm_k,"so no change was made"))
        #     }
        #   }
        # }
        # Save file over previous
        # dput(x=rdat_i,file=fileName_rdat_i)
        file_to <- file.path("..",dir_bam_sim)
      }else{
        rdat_i <- rdat_base
        # Collect the most important values from rdat
        parms_i <- unlist(rdat_i$parms)
        parm.cons.result_i <- unlist(lapply(rdat_i$parm.cons,function(x){x[8]}))
        like_i <- rdat_i$like
        sdnr_i <- rdat_i$sdnr
        vals_i <- c(parms_i[!names(parms_i)%in%names(parm.cons.result_i)],
                    parm.cons.result_i,
                    like_i,
                    sdnr_i)*NA
        file_to <- file.path("..",dir_bam_sim_fail)
        if(!dir.exists(file_to)){
          dir.create(file_to)
        }
      }

      ########## Copy data files to folder
      file.copy(from=file.path(dir_bam_i,c(fileName_dat_i,fileName_rdat_i)),
                to=file_to,
                overwrite = TRUE)

      #######Remove individual processing folders
      setwd(wd)
      unlink(sim_dir_i, recursive=T,force=TRUE)

      return(c(setNames(c(lk_total_i,grad_max_i),c("lk_total","grad_max")),vals_i))
    } #end MCBE foreach loop


    message(paste("Finished running MCBE for", nsim,"sims in parallel on",ncores,"cores at",
                  format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

    est_out <- apply(do.call(rbind,est_out),2,as.numeric)
    est_out <- cbind("sim"=nm_sim,as.data.frame(est_out)
                      )

    write.csv(x=est_out,file=file.path(dir_bam_sim,"MCBE_results.csv"),row.names = FALSE)

    if(unlink_dir_bam_sim){
      # delete a directory -- must add recursive = TRUE
      unlink(dir_bam_sim, recursive = TRUE,force=TRUE)
    }

    ### Return stuff
    out <- est_out

  }else{
    out <- sim_out
  } # end if(run_est)

  ## parallel shut down
  parallel::stopCluster(cl)

  invisible(out)
}
