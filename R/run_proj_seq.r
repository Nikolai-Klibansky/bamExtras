#' Sequence of projection and estimation
#'
#' This function runs a sequence of projections with \code{run_proj}, and estimation
#' with \code{run_bam}, to build long-term closed-loop simulations. If a single
#' \code{dat} object is provided, it will run a single deterministic projection.
#' If a list of \code{dat} objects is provided, for example as output from
#' \code{run_MCBE(run_est=FALSE)}, the function will run
#' multiple projections in parallel.
#' @param CommonName Common name of species associated with dat, tpl, and cxx files
#' @param fileName Name given to BAM files, not including file extensions.
#' @param dat_obj dat file read in as a character vector with readLines(con=dat_file)
#' @param tpl_obj tpl file read in as a character vector with readLines(con=tpl_file)
#' @param cxx_obj cxx file read in as a character vector with readLines(con=cxx_file)
#' @param rdat (list) object read in with dget(). Specifically, the BAM output file produced by the cxx file.
#' @param args_run_proj List of arguments to pass to \code{run_proj}
#' @param nyps Total number of years to project which may include multiple assessment cycles. nyps = Number of Years to Project in Sequence
#' @param plot_type Indicate which bam results should be plotted by calling \code{plot_bam}. Current options
#' include "all", "final", or "none", to plot results from all stock assessments in the projection
#' sequence, the final assessment, or none of the assessments.
#' @param write_bam_files indicate which bam files should be written. Current options
#' include "all", "final", or "none", to write files from all stock assessments in the projection
#' sequence, the final assessment, or none of the assessments.
#' @param parallel logical. Should parallel processing be used? This setting is
#' only used when \code{dat_obj} is a list of multiple dat objects. The function
#' will not use parallel processing for a single dat object (it would actually
#' be somewhat slower).
#' @param ncores number of cores to use for parallel processing
#' @returns By default, returns the bam model object from the end of the projection
#' sequence.
#' @keywords bam stock assessment fisheries population dynamics
#' @author Nikolai Klibansky
#' @export

run_proj_seq <- function(CommonName = NULL,
                         fileName = "bam",
                         dat_obj = NULL,
                         tpl_obj = NULL,
                         cxx_obj = NULL,
                         rdat_obj = NULL,
                         args_run_proj = list(),
                         args_run_MCBE=list(run_mc=FALSE,sim_type_return = "dat"),
                         nyps = NULL,
                         plot_type="final",
                         write_bam_files = "final",
                         parallel = TRUE,
                         ncores=NULL


){
  draw_plot <- ifelse(plot_type=="all",TRUE,FALSE)
  wdt <- getwd()
  message(paste("working directory:",wdt))

  args_run_proj_default <- as.list(formals(run_proj))

  args_run_proj <- modifyList(args_run_proj_default,args_run_proj)

  args_run_MCBE_default <- modifyList(as.list(formals(run_MCBE)),
                                      list(run_est=FALSE,run_bam_base = FALSE,nsim=1
                                      )
                                      )
  args_run_MCBE <- modifyList(args_run_MCBE_default,args_run_MCBE)


  nyp <- args_run_proj$nyp
  # if(!"nyp"%in%names(args_run_proj)){
  #   nyp <- formals(run_proj)$nyp
  # }else{
  #   nyp <- args_run_proj$nyp
  # }
  if(is.null(nyps)){
    n_cyc <- 2 # Number of projection cycles default to 2
    nyps <- nyp*n_cyc # get default length of projections
  }else{
    n_cyc <- ceiling(nyps/nyp)
  }

  # Determine type
  if(is.list(dat_obj)){
    nsim <- length(dat_obj)
  }else{
    nsim <- 1
    parallel <- FALSE
  }

  `%dowhat%` <- ifelse(parallel, `%dopar%`, `%do%`)

  nm_sim <- paste0("sim_",sprintf(paste("%0",nchar(nsim),".0f",sep=""),1:nsim))

  #### parallel setup
  if(parallel){
    if(is.null(ncores)){
      coresAvail <- parallel::detectCores()
      ncores <- coresAvail-1
    }
    ncores  <- min(c(ncores,coresAvail))
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
  }

  # Run multiple projection
  proj_seq_result <- foreach::foreach(i=1:nsim,
                                 .packages="bamExtras"
  ) %dowhat% {
    nm_sim_i <- nm_sim[i]
    fileName_i <- if(nsim==1){
      fileName
    }else{
      file.path("sim",nm_sim_i,fileName)
    }

    dat <-  if(is.null(dat_obj)){
      get(paste0("dat_",CommonName))
    }else{
      if(!is.list(dat_obj)){
        dat_obj
      }else{
        dat_obj[[i]]
      }
    }
    tpl <- if(is.null(tpl_obj)){
      get(paste0("tpl_",CommonName))
    }else{
      tpl_obj
    }
    cxx <- if(is.null(cxx_obj)){
      get(paste0("cxx_",CommonName))
    }else{
      cxx_obj
    }
    rdat <- if(is.null(rdat_obj)){
      get(paste0("rdat_",CommonName))
    }else{
      rdat_obj
    }

  ## Run a single projection sequence
  # For each cycle of the projection sequence
  for(j in 1:n_cyc){
  endyr <- rdat$parms$endyr
  endyr_proj <- as.numeric(endyr) + nyp
  dir_curr <- file.path(fileName_i,paste0("endyr_",endyr))
  dir_proj <- file.path(fileName_i,paste0("endyr_",endyr_proj))
  dir_curr_figs <- "figs" #file.path(dir_curr,"figs")
  dir_proj_figs <- "figs" #file.path(dir_proj,"figs")

  # if(write_bam_files!="none"){
    if(write_bam_files=="all"){
      if(!dir.exists(dir_curr)){
        dir.create(dir_curr,recursive = TRUE)
      }
      if(!dir.exists(dir_proj)){
        dir.create(dir_proj,recursive=TRUE)
      }
      fp_dat <- file.path(dir_curr, "bam.dat")
      fp_tpl <- file.path(dir_curr, "bam.tpl")
      fp_cxx <- file.path(dir_curr, "bam.cxx")
      fp_rdat <- file.path(dir_curr,"bam.rdat")

      if(!exists(fp_dat)){
        writeLines(dat,fp_dat)
      }
      if(!exists(fp_tpl)){
        writeLines(tpl,fp_tpl)
      }
      if(!exists(fp_cxx)){
        writeLines(cxx,fp_cxx)
      }
      if(!exists(fp_rdat)){
        dput(rdat,fp_rdat)
      }

      ########## Copy data files to folder
      if(j==1&draw_plot){
        # plot initial (base) run
        setwd(dir_curr)
        message(paste("Working directory temporarily changed to:\n",dir_curr))
        plot_bam(rdat)

        if(!dir.exists(dir_curr_figs)){
          dir.create(dir_curr_figs,recursive = TRUE)
        }
        file.copy(from=file.path("spp-figs",list.files("spp-figs")),
                  to=dir_curr_figs,recursive = TRUE)
        unlink("spp-figs", recursive=T)
        setwd(wdt)
        message(paste("Working directory changed back to:\n",wdt))
      }

      unlink_dir_bam_update <- FALSE
  }else{
    dir_proj <- fileName_i
    dir_proj_figs <- "figs"#file.path(dir_proj,"figs")
    if(!dir.exists(dir_proj)){
      dir.create(dir_proj,recursive=TRUE)
    }
    unlink_dir_bam_update <- ifelse(write_bam_files=="none"|(write_bam_files=="final"&j<n_cyc),TRUE,FALSE)
  }

  # get unprojected bam
  bam <- do.call(bam2r,list(dat_obj=dat,tpl_obj=tpl,cxx_obj=cxx))
  init <- bam$init
  styr_bam <- init$styr
  endyr_bam <- init$endyr
  yr_bam <- styr_bam:endyr_bam

  # project BAM
  args_run_proj$rdat <- rdat
  args_run_proj <- modifyList(args_run_proj,
                              list(args_bam2r = list(dat_obj=dat,tpl_obj=tpl,cxx_obj=cxx),
                                   plot=draw_plot)
                              )


  bam_proj <- do.call(run_proj,args_run_proj)
  message(paste("run_proj ran successfully to",endyr_proj))

  # get projected bam
  bam_p <- bam_proj$bam_p
  rdat_proj <- bam_proj$rdat_proj

  # Completely replace values in args_run_MCBE (modifyList will try to replace elements of sublists like rdat_base_obj, which can cause errors here)
  args_run_MCBE$bam <- bam_p
  args_run_MCBE$rdat_base_obj <- rdat_proj

  # Add observation error to projected bam (NOTE: adds error to all years)
  dat_sim <- do.call(run_MCBE,args_run_MCBE)[[1]]
  bam_p <- bam2r(dat_obj = dat_sim,tpl_obj = bam_p$tpl,cxx_obj = bam_p$cxx)
  init_p <- bam_p$init

  if(j>1){ # After the first assessment cycle, overwrite values of init_p for with values from base years. That is, keep the data from the previous assessment for past years.
    # landings, discards, indices, and cvs
    nm_obs_ts <- names(init)[grepl("^obs_L|^obs_cv_L|^obs_released|^obs_cv_D|^obs_cpue|^obs_cv_cpue",names(init))]
    for(k in seq_along(nm_obs_ts)){
      nmk <- nm_obs_ts[k]
      xk <- init[[nmk]]
      init_p[[nmk]][names(xk)] <- xk
    }

    # age and length compositions
    nm_obs_comp <- names(init)[grepl("^obs_agec|^obs_lenc",names(init))]
    for(k in seq_along(nm_obs_comp)){
      nmk <- nm_obs_comp[k]
      xk <- init[[nmk]]
      init_p[[nmk]][rownames(xk),] <- xk
    }
  }

  # add modified init_p back intor bam_p
  bam_p <- bam2r(dat_obj = bam_p$dat,tpl_obj = bam_p$tpl,cxx_obj = bam_p$cxx, init=init_p)

  bam_update <- run_bam(bam=bam_p,dir_bam=dir_proj,unlink_dir_bam = unlink_dir_bam_update,
                        return_obj=c("dat","tpl","cxx","rdat")
                        )

  # Plot bam update results
  if(plot_type=="final"&j==n_cyc){draw_plot <- TRUE}

  if(draw_plot){
    setwd(dir_proj)
    message(paste("Working directory temporarily changed to:\n",dir_proj))
    plot_bam(bam_update$rdat)

    ########## Copy data files to folder
    if(!dir.exists(dir_proj_figs)){
      dir.create(dir_proj_figs,recursive=TRUE)
    }
    file.copy(from=file.path("spp-figs",list.files("spp-figs")),
              to=dir_proj_figs,recursive = TRUE)
    unlink("spp-figs", recursive=T)
    setwd(wdt)
    message(paste("Working directory changed back to:\n",wdt))
  }

  # update object values
    dat <-  bam_update$dat
    tpl <-  bam_update$tpl
    cxx <-  bam_update$cxx
    rdat <- rdat_proj

  # update catch advice based on estimation model
    rdat_sim <- bam_update$rdat
    parms_sim <- rdat_sim$parms
    ts_sim <- rdat_sim$t.series
    nyb_rcn <- args_run_proj$nyb_rcn
    yrs_L_b <- paste((parms_sim$endyr-(nyb_rcn$L-1)):parms_sim$endyr) # years of recent landings (i.e. from the base model)
    F_cur_sim <- bamExtras::geomean2(ts_sim[yrs_L_b,"F.full"],na.rm=TRUE)
    F_proj_sim <- parms_sim$Fmsy # NOTE: NEED TO UPDATE TO ACCOMODATE MODELS THAT DON'T USE Fmsy!!!! 2026-01-30 NPK

    args_run_proj$F_cur <- F_cur_sim
    args_run_proj$F_proj <- F_proj_sim
  } # end for(j in seq_along(n_cyc))

    return(list("bam_proj"=bam_proj,
                "bam_update"=bam_update
                )
    )
  }##### END FOREACH LOOP ######
  names(proj_seq_result) <- nm_sim

  if(parallel){
    ## parallel shut down
    parallel::stopCluster(cl)
  }

  # Return bam from final cycle in sequence
  invisible(proj_seq_result)
}
