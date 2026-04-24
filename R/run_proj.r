#' run_proj
#'
#' Project BAM forward in time to conduct deterministic or stochastic projections.
#' @param rdat BAM output rdat (list) object read in with dget().
#' @param args_run_bam list of arguments passed internally to run_bam if any are provided
#' @param args_bam2r list of arguments passed internally to bam2r if any are provided
#' @param stochastic logical. Is this a stochastic projection? Mostly runs the same code
#'  either way, but stochastic projections add normal errors to recruitment. Stochastic
#'  setting is intended to be run with rdat input from MCBE.
#' @param pstar Value to use for applying pstar scaling of Fmsy. If set to NULL, it is not used and pstar code will not be run.
#' @param styr_proj start year (i.e. first calendar year) of projection
#' @param nyb_rcn number of years in calculations that include recent landings
#' (L) or discards and their cvs, index cvs (U), recruitment (R), or numbers of
#' samples in age or length compositions (comp). Fleets without L or D within the
#' last nyb_rcn years of the base model will not be projected forward. named list
#' @param nyp_cur  number of years to maintain current conditions before implementing management
#' @param nyp  number of years of projection. Must be >nyp_cur
#' @param nm_yr_p names of \code{init} objects to project forward by \code{nyp}
#' @param F_cur current fully selected fishing mortality rate
#' @param F_proj fully selected fishing mortality rate
#' @param L_cur current level of landings. (2026-01-30 NPK: Not yet implemented in run_proj but will be used to set catch in projections instead of F).
#' @param ages   ages. numeric vector
#' @param N_styr_proj abundance at age in styr_proj. numeric vector
#' @param S_styr_proj spawning stock size at age in styr_proj. Often biomass in mt, sometimes eggs in n. numeric vector
#' @param M natural mortality rate, at age
#' @param wgt_mt weight at age of population in metric tons (mt). numeric vector
#' @param wgt_L_klb  weight at age of landings in thousand pounds (klb). numeric vector
#' @param wgt_D_klb  weight at age of discards in thousand pounds (klb). numeric vector
#' @param wgt_F_flt_klb list of weight vectors (klb by age) or matrices (by year,age) for each fleet associated with landings or discards. numeric vector or matrix
#' @param len_F_flt_mm list of length vectors (mm by age) or matrices (by year,age) for each fleet associated with landings or discards. numeric vector or matrix
#' @param reprod reproductive contribution at age to SSB. numeric vector
#' @param sel_L  selectivity at age to compute landings. numeric vector
#' @param sel_D  selectivity at age to compute dead discards. numeric vector
#' @param sel_tot  selectivity at age to compute Z (includes landings and discards). numeric vector
#' @param sel_F_flt list of selectivity vectors (by age) or matrices (by year,age) for each fleet associated with landings or discards.
#' @param spawn_time time of year for peak spawning
#' @param SR_par list of parameters associated with the stock-recruit (SR) relationship.
#' Currently only works with a Beverton-Holt SR relationship for which the parameters are:
#' h = steepness of spawner-recruit function, R0 = virgin recruitment of spawner-recruit
#' function, Phi0 = virgin spawners per recruit, biascorr = bias correction
#' @param SR_method Spawner-recruit function BH = Beverton-Holt, R0 = model estimated
#' mean recruitment (R0), GM = constant geometric mean of nyb_rcn$R recent years of
#' t.series$recruits. By default, run_proj will try to use BH, or R0, but if it can't
#' find all the parameters, it will resort to GM. Note GM method has not been reviewed
#' and should be used with caution.
#' @param age_error age error matrix to use for the projection years. By default, the function uses the same age_error matrix from the base model
#' @param project_bam logical. After the projection is run, should the function build a new set of projected bam files?
#' @param args_ia list of arguments to pass to interim adjustment function. See Details section below.
#' @param model_obs list of optional arguments for incorporating observation error into simulated sampling when building a dat file when project_bam=TRUE. See Details section below. NOT CURRENTLY USED 2025-09-17.
#' @param plot logical. Produce plots of extended base model including projected data
#' @param nm_comp_sfx Character vector. Possible suffix added to fleet abbreviation in comps
#' (e.g. c("a","b") in SEDAR 82 gray triggerfish). Used to match comps with catch when projecting
#' comp data when \code{project_bam=TRUE}
#' @param t_series_na_vals Numeric vector of values to consider as NA in t.series object.
#' @param key_lenprob Optional list for specifying cv values for length-at-age,
#' for each set of length compositions. (e.g. key_lenprob = c("D.rHD"="len.cv.L", "L.rGN"="len.cv.L")). Use keyword 'all', as in default to set all
#' values to a particular value in parm.cons.
#' @details
#' \code{args_ia}:
#' \itemize{
#' \item{U_nm}{The fleet abbreviation for an index of abundance that you want to use to adjust F (e.g. sTV, sCT). If set to \code{NULL}, no interim adjustment will be used.}
#' \item{type}{Method used for computing F adjustment. "Uprop" computes the mean U from recent years, divided by U from
#' reference year. The reference year is the terminal year of the last full stock assessment. The set of recent years is computed
#' as the current projection year plus yr_U_lag. By default, type = "Uprop"}
#' \item{yr_U_lag}{A vector of negative integers used to compute the set of years used for computing an interim adjustment
#' from the reference index (comma separated example values: -1, -1:-3, -2:4). By default, yr_U_lag = -(2:4).}
#' \item{yr_ia_by}{Frequency of interim adjustments, in year units. Used to determine the years between stock
#' assessments when interim adjustments will be computed and applied to F. The interim adjustment years are computed as
#' \code{yrlim_p <- endyr + c(1,nyp); yr_ia <- seq(yrlim_p[1], yrlim_p[2], by = yr_ia_by)}. By default, yr_ia_by = 5}
#' \item{include_yrlim}{Logical. Should interim analysis be conducted in the first and last years of the projection period?
#' By default, include_yrlim_p = FALSE}
#' }
#' \code{model_obs}:
#' Values supplied to model_obs should be provided as a vector of function class objects of length 1 or length equal to the number of fleets associated with the data type. For cv_U, ntrip, and nsamp, these function can use the argument x, where x is the geometric mean of values for a fleet for the nyb_rcn years of the assessment.
#' \itemize{
#' \item{cv_L}{Function for computing cvs of landings during projection period.}
#' \item{cv_D}{Function for computing cvs of discards during projection period.}
#' \item{cv_U}{Function for computing cvs of indices during projection period.}
#' \item{ntrip}{Function for computing the number of fish sampled for comps during projection period.}
#' \item{nsamp}{Function for computing the number of trips sampled for comps during projection period.}
#' }
#' @keywords bam stock assessment fisheries population dynamics
#' @author Kyle Shertzer, Erik Williams, and Nikolai Klibansky
#' @export
#' @examples
#' \dontrun{
#' # Run deterministic projection with defaults
#' proj_VeSn <- run_proj(rdat_VermilionSnapper,plot=TRUE)
#'
#' # Run deterministic projection with defaults
#' # and project bam inputs, by adding the bam files with args_bam2r
#' proj_VeSn <- run_proj(rdat_VermilionSnapper,
#' project_bam=TRUE,
#' args_bam2r = list(
#' dat_obj=dat_VermilionSnapper,
#' tpl_obj=tpl_VermilionSnapper,
#' cxx_obj=cxx_VermilionSnapper
#' ))
#' # The projections work for most of the assessments
#' proj_BlSb <- run_proj(rdat_BlackSeaBass,plot=TRUE)
#' proj_GaGr <- run_proj(rdat_GagGrouper,plot=TRUE)
#' proj_RdGr <- run_proj(rdat_RedGrouper,plot=TRUE)
#' proj_RdSn <- run_proj(rdat_RedSnapper,plot=TRUE)
#' proj_VeSn <- run_proj(rdat_VermilionSnapper,plot=TRUE)
#' # Some need additional arguments to run correctly
#' proj_GrTr <- run_proj(rdat_GrayTriggerfish,
#' nm_comp_sfx = c("a","b"),
#' plot=TRUE,
#' key_lenprob = c("D.rHDs"="len.cv.L", "L.rGNs"="len.cv.L")
#' )
#' }

run_proj <- function(rdat = NULL,
                     args_run_bam = NULL,
                     args_bam2r = NULL,
                     stochastic = FALSE,
                     pstar = NULL,
                     styr_proj = NULL,
                     nyb_rcn = list(L=3, U=3, R=5, comp=3),
                     nyp_cur = 3,
                     nyp = 5,
                     nm_yr_p=c("endyr","endyr_dev_rec","endyr_rec_dev","endyr_rec_phase2","endyr_rec_spr","styr_regs","endyr_proj"),
                     F_cur = NULL,
                     F_proj = NULL,
                     L_cur = NULL,
                     ages = NULL,
                     N_styr_proj = NULL,
                     S_styr_proj = NULL,
                     sel_L = NULL,
                     sel_D=NULL,
                     sel_tot = NULL,
                     sel_F_flt = NULL,
                     wgt_mt = NULL,
                     wgt_L_klb = NULL,
                     wgt_D_klb = NULL,
                     wgt_F_flt_klb = NULL,
                     len_F_flt_mm = NULL,
                     reprod = NULL,
                     M = NULL,
                     spawn_time = NULL,
                     SR_par = NULL,
                     SR_method = "BH",
                     age_error = NULL,
                     project_bam = FALSE,
                     args_ia = list(U_nm = NULL),
                     model_obs = list(
                       cv_L=  function(x=0.2,n=nyp){rep(x,n)}, # cvs of landings
                       cv_D=  function(x=0.2,n=nyp){rep(x,n)}, # cvs of discards
                       cv_U=  function(x,    n=nyp){rep(x,n)}, # cvs of indices
                       nfish= function(x,    n=nyp){rep(x,n)}, # number of fish for comps
                       ntrip= function(x,    n=nyp){rep(x,n)}  # number of fish for length comps
                     ),
                     plot = FALSE,
                     nm_comp_sfx = c(""),
                     t_series_na_vals=c(-9999,-99999,-999999),
                     key_lenprob = c("all"="^len.cv")

) {
library(ggplot2)
library(tidyr)
mt2klb <- 2.20462              # conversion of metric tons to 1000 lb

# Compute adjustment to target fishing mortality (F_target)
# An internal function for now

# data: data frame of index values by year, possibly for multiple indices
# U_nm: name of index you want to use to adjust F between stock assessments (must match column name of index in data)
# yr_ref: reference year to compare current index value with
# yr: year(s) to use for computing the current index value
F_adjust <- function(data, U_nm, yr_ref, yr,
                     type = "Uprop"){
  Ux <- data[,U_nm,drop=FALSE]
  Uyr <- Ux[paste(yr),]
  Uref <- Ux[paste(yr_ref),]
  switch(type,
         Uprop = mean(Uyr,na.rm=TRUE)/Uref
  )
}

args_ia_default = list(
  type = "Uprop", yr_U_lag = -(2:4), yr_ia_by = 5, include_yrlim_p=FALSE
  )
args_ia <- modifyList(args_ia_default,args_ia)

if(!is.null(args_run_bam)){
  if(!"return_obj"%in%names(args_run_bam)){
    args_run_bam <- c(args_run_bam,list(return_obj=c("dat","tpl","cxx","rdat")))
  }
  run_bam_out <- do.call(run_bam,args_run_bam)
  rdat <- run_bam_out$rdat
}

### Get parts of the rdat used later in the function and compute stuff
  # if(!is.null(rdat)){
   parms <- rdat$parms
   a.series <- rdat$a.series
   t.series <- rdat$t.series
   parm.cons <- rdat$parm.cons
   sel.age <- rdat$sel.age
   CLD.est.mats <- rdat$CLD.est.mats
   Z.age <- rdat$Z.age
   N.age <- rdat$N.age
   N.age.mdyr <- rdat$N.age.mdyr
   size.age.fishery <- rdat$size.age.fishery
   comp.mats <- rdat$comp.mats

   nm_parms <- names(parms)

    styr <- parms$styr
    endyr <- parms$endyr
    yb <- styr:endyr # years of the base model
    nyb <- length(yb)

    if(any(unlist(t.series)%in%t_series_na_vals)){
      t.series <- apply(t.series,2,
                        function(x){x[x%in%t_series_na_vals] <- NA
                        x}
                        )
      t.series <- as.data.frame(t.series)
      message(paste0("values in t.series in the set '", paste(t_series_na_vals,collapse=", "), "' found and changed to NA"))
    }
    # t.series[t.series==-99999] <- NA

    sel_age_1 <- sel.age[names(sel.age)%in%c("sel.v.wgted.L","sel.v.wgted.D","sel.v.wgted.tot")]
    sel_age_2 <- sel.age[!names(sel.age)%in%names(sel_age_1)]

    yrs_L_b <- paste((endyr-(nyb_rcn$L-1)):endyr) # years of recent landings (i.e. from the base model)
    yrs_R_b <- paste((endyr-(nyb_rcn$R-1)):endyr)
    R_b <- t.series[yrs_R_b,"recruits"]
    R_b_gm <- bamExtras::geomean2(R_b)

    # Identify F-at-age for each fleet in endyr
    # L = (F/Z)*N*(1-exp(-Z))
    # L/(N*(1-exp(-Z))) = F/Z
    # L*Z/(N*(1-exp(-Z))) = F
    # F = L*Z/(N*(1-exp(-Z)))
    nm_Cn_CLD <- names(CLD.est.mats)[grepl("^(Ln|Dn)((?!total).)*$",names(CLD.est.mats),perl=TRUE)]
    Cn <- CLD.est.mats[nm_Cn_CLD ]
    names(Cn) <- gsub("^([LD])(n)(.*)","Cn.\\1\\3",names(Cn))
    key_nm_Cn <- setNames(names(Cn),nm_Cn_CLD)

    nm_Cw_CLD <- names(CLD.est.mats)[grepl("^(Lw|Dw)((?!total).)*$",names(CLD.est.mats),perl=TRUE)]
    Cw <- CLD.est.mats[nm_Cw_CLD]
    names(Cw) <- gsub("^([LD])(w)(.*)","Cw.\\1\\3",names(Cw))
    key_nm_Cw <- setNames(names(Cw),nm_Cw_CLD)

    cv_C <- t.series[paste(yb),grepl("^cv.[DL]",names(t.series))] # cvs associated with components of the catch (removals) during the base years
    Z <- Z.age[paste(yb),]
    N <- N.age[paste(yb),]
    Nmdyr <- N.age.mdyr[paste(yb),]
    logRdev <- rdat$t.series[paste(yb),"logR.dev"]
    if("SSB"%in%names(t.series)){
    S <- t.series[paste(yb),"SSB"]
    }else{
      warning("SSB not found in t.series. run_proj is not sure what to use to characterize stock size.")
    }
    R <- N[,1] # Recruits
    # Like F_fleet in bam tpl (e.g. F_cHL). Same computation as bam
    F_flt <- lapply(Cn,function(x){x*Z/(N*(1-exp(-Z)))}) # Computes a list of F (year,age) for each fleet, which are not in the rdat
    names(F_flt) <- gsub("^Cn","F",names(F_flt))
    Fage <- Reduce("+",F_flt)
    Ffull <- apply(Fage,1,max)

    if(is.null(ages)){ages <- a.series$age}
    nages <- length(ages)
    len <- rdat$a.series$length
    if(is.null(age_error)){age_error <- rdat$age.error$error.mat
    if(is.null(age_error)){ # If there is no age error matrix in the rdat, just use an identity matrix
      age_error <- diag(length(ages))
      dimnames(age_error) <- list("age"=ages,"age"=ages)
    }
    }
    if(is.null(spawn_time)){spawn_time <- parms$spawn.time}
    if(is.null(reprod)){reprod <- a.series$reprod}
    if(is.null(styr_proj)){styr_proj <- endyr+1}
    yp <- styr_proj:(styr_proj+nyp-1) # years of the projection period
    if(is.null(F_cur)){F_cur <- bamExtras::geomean2(t.series[yrs_L_b,"F.full"],na.rm=TRUE)}
    if(is.null(F_proj)){F_proj <- parms$Fmsy}
    if(is.null(L_cur)){
      is.total.L.klb <- "total.L.klb"%in%names(t.series)
      if(!is.total.L.klb){
        total.L.name <- names(t.series)[grepl("total.L",names(t.series))][1]
        message(paste("total.L.klb not found in t.series. L_cur is computed from",total.L.name,"instead\n"))
        total.L <- t.series[,total.L.name,drop=FALSE]
      }else{
        total.L <- t.series[,"total.L.klb",drop=FALSE]
        message(paste("L_cur is computed from t.series$total.L.klb\n"))

      }
      L_cur <- mean(total.L[yrs_L_b,])
      }


    if(is.null(sel_L)){sel_L <- sel_age_1$sel.v.wgted.L}

    if(is.null(sel_D)){
      if(!is.null(sel_age_1$sel.v.wgted.D)){
        sel_D <- sel_age_1$sel.v.wgted.D
      }else{
        message("sel_age_1$sel.v.wgted.D not found. Assessment may not model discards.\n")
      }
    }

    if(is.null(sel_tot)){sel_tot <- sel_age_1$sel.v.wgted.tot}

    if(is.null(sel_F_flt)){
      abb_LD <- gsub(".pr$","",names(t.series)[grepl("^[LD].*.pr$",names(t.series),perl=TRUE)]) # landings abbreviations
      abb_F_1 <- paste0("F.",gsub("^(D.)(.*)","\\2.D",gsub("^L.","",abb_LD))) # pattern 1 (typical)
      abb_F_2 <- paste0("F.",abb_LD) # pattern 2 (less typical but preferred)
      if(all(abb_F_1%in%names(t.series))){
        abb_F <- abb_F_1
      }else if(all(abb_F_2%in%names(t.series))){
        abb_F <- abb_F_2
      }else {
        warning(paste0("Can't match Fs. Neither set 1 (",paste(abb_F_1,collapse=", "),
                       ") nor set 2 (",paste(abb_F_2,collapse=", "),
                       ") column names can all be found in t.series"))
      }
      Fsum_flt <- t.series[paste(yb),abb_F] # F time series for each fleet
      names(Fsum_flt) <- paste0("F.",abb_LD)
      key_abb_F <- setNames(abb_F,names(Fsum_flt)) # used to match up column names in different structures. Should eventually change the t.series names to the preferred format.
      Fsum <- rowSums(Fsum_flt) # Should be equal to t.series$Fsum

      # Compute selectivity for each fleet associated with an F value.
      # NOTE: Computing selectivities is somewhat more reliable than trying to find
      # them in the rdat, in instances when selectivities from one fleet are used
      # for multiple fleets (e.g. when the headboat selectivity is used for the MRIP landings)
      # 2025-09-08 NPK should probably make sure selectivities for all fleets are in the rdat, even when duplicating
      sel_F_flt <- lapply(names(F_flt),function(x){
        a <- pmax(pmin(F_flt[[x]]/Fsum_flt[,x],1),0)
        a[is.na(a)] <- 0 # Replace NA with zero (NaN will occur when Fsum_flt values are zero, as in years when fleets do not have any landings)
        a
      })
      names(sel_F_flt) <- gsub("^F.","sel.",names(F_flt))

      F_ybgm_flt <- apply(tail(Fsum_flt,nyb_rcn$L),2,bamExtras::geomean2) # geomean Fsum by fleet during the last nyb_rcn years of the base model
      F_prop_flt <- F_ybgm_flt/sum(F_ybgm_flt) # Proportion of Fsum attributed to each fleet at the end of the assessment

      # This is just F_ybgm_flt multiplied by the selectivity in the endyr of the base model for each fleet
      # In the bam tpl, these vectors are named F_end although they are summed by landings (F_end_L), discards (F_end_D), or both (F_end)
      # It's basically Fsum-at-age at the end of the assessment.
      F_ybgm_flt_a <- as.data.frame(
        lapply(1:length(sel_F_flt),function(i){
          nm_x <- names(sel_F_flt)[i]
          sel_endyr_x <- sel_F_flt[[nm_x]][paste(endyr),]
          F_ybgm_flt_x <- F_ybgm_flt[[gsub("^sel.","F.",nm_x)]]
          F_ybgm_flt_x*sel_endyr_x
        }
        )
      )
      names(F_ybgm_flt_a) <- names(F_ybgm_flt)

      # Sum across fleets (these objects are named as in bam tpl)
      F_end_L <- rowSums(F_ybgm_flt_a[grepl("^F.L",names(F_ybgm_flt_a))])
      F_end_D <- rowSums(F_ybgm_flt_a[grepl("^F.D",names(F_ybgm_flt_a))])
      F_end <- rowSums(as.data.frame(F_ybgm_flt_a))
      F_end_apex <- max(F_end)
      sel_wgted_tot <- F_end/F_end_apex # Should equal to rdat$sel.age$sel.v.wgted.tot
      sel_wgted_L <-   sel_wgted_tot*(F_end_L/F_end) # Should equal rdat$sel.age$sel.v.wgted.L # F_end_L/F_end_apex # as in tpl
      sel_wgted_D <-   sel_wgted_tot*(F_end_D/F_end) # Should equal rdat$sel.age$sel.v.wgted.D # F_end_D/F_end_apex # as in tpl
      sel_wgted_F_flt <- sel_wgted_tot*F_ybgm_flt_a/F_end # by extension, scale selectivity at age for each fleet
      #sel_wgted_F_flt <- F_ybgm_flt_a/F_end_apex     # by extension, scale selectivity at age for each fleet # same calculation as above line
      # Note that: (sel_wgted_L+sel_wgted_D) == rowSums(sel_wgted_F_flt)
    }

    ## Weights of fish (year, age)
    if(is.null(wgt_mt)){wgt_mt <- a.series$wgt.mt}
    wgt_klb <- wgt_mt*mt2klb
    wgt_b_klb <- matrix(wgt_klb,nrow=nyb,ncol=nages,byrow=TRUE,dimnames=list(yb,ages))
    if(is.null(wgt_L_klb)){
      wgt.wgted.L.klb_nm <- names(a.series)[grepl("^[A-Za-z]*wgt.wgted.L.klb",names(a.series))]
      if(length(wgt.wgted.L.klb_nm)>0){
        message(paste0(wgt.wgted.L.klb_nm, " found in names(a.series) and used to set wgt_L_klb\n"))
        wgt_L_klb <- a.series[,wgt.wgted.L.klb_nm[1]]

      }else{
        warning("no wgt.wgted.L.klb found in names(a.series).\n")
      }
    }
    if(is.null(wgt_D_klb)){
      wgt.wgted.D.klb_nm <- names(a.series)[grepl("^[A-Za-z]*wgt.wgted.D.klb",names(a.series))]
      if(length(wgt.wgted.D.klb_nm)>0){
        message(paste0(wgt.wgted.D.klb_nm, " found in names(a.series) and used to set wgt_D_klb\n"))
        wgt_D_klb <- a.series[,wgt.wgted.D.klb_nm[1]]

      }else{
        warning("no wgt.wgted.D.klb found in names(a.series). Does this assessment model discards?\n")
      }
    }

    ## Weights of fish (year, age) by fleet.
    if(is.null(wgt_F_flt_klb)){
      # initialize
      wgt_F_flt_klb <- size.age.fishery[grepl("^[a-zA_Z]*.*wgt",names(size.age.fishery))]

      nm_unit_wgt <- c("lb","klb")
      wgt_F_flt_klb_type <- gsub(".*\\<([a-zA-Z]*wgt[a-zA-Z]*)\\>.*","\\1",names(wgt_F_flt_klb))
      wgt_F_flt_klb_units <- gsub(paste0(".*\\<([:alpha:]*",paste(nm_unit_wgt,collapse = "|"),"[:alpha:]*)\\>.*"),
                              "\\1",names(wgt_F_flt_klb))

      # Reformat names to wgt.fleetType.fleet which should be a combination of
      # type (e.g. wgt or wholewgt), unit (e.g. lb or klb) and fleet abbreviation (e.g. rHB or cHL)
      # (e.g. wgt.L.cHL for weight of landings in the commercial hook and line fleet)
      wgt_F_flt_abb <- local({
        a <- gsub(paste0(c(unique(wgt_F_flt_klb_type),unique(wgt_F_flt_klb_units)),collapse="|"),"",names(wgt_F_flt_klb))
        b <- gsub("^\\.*|\\.*$","",a)

        #a <- gsub("^([a-zA-Z]+)(.)(.*?)(.)([a-zA-Z]+)$","\\3",names(wgt_F_flt_klb),perl=TRUE)
        c <- gsub("(.*?)(?<=.)(\\.D)$","D.\\1",b,perl=TRUE) # Change .D suffix to D. prefix
        d <- gsub("^([a-zA-Z]{2})(D)","D.\\1\\2",c)         # Add D. prefix if third letter is D
        e <- gsub("^([cr])(.*)","L.\\1\\2",d) # Add L. prefix to comps starting with c or r (landings not discards. won't modify fishery independent survey comps)
        e
      })

      # wgt_F_flt_klb_units <- gsub("^([a-zA-Z]+)(.)(.*?)(.)([a-zA-Z]+)$","\\1.\\5",names(wgt_F_flt_klb),perl=TRUE)
      names(wgt_F_flt_klb) <- names(wgt_F_flt_klb_units) <- paste("wgt", #wgt_F_flt_klb_units,
                                                          wgt_F_flt_abb,sep=".")

      # Convert weights in lb to klb
      for(nm_i in names(wgt_F_flt_klb)){
        xi <- wgt_F_flt_klb[[nm_i]]
        if(wgt_F_flt_klb_units[[nm_i]]=="klb"){
          # do nothing
        } else if(wgt_F_flt_klb_units[[nm_i]]=="lb"){
          wgt_F_flt_klb[[nm_i]] <- xi/1000
          wgt_F_flt_klb_units[[nm_i]] <- "klb"
        } else {
          warning(paste("units of",nm_i,"are not lb or klb"))
        }
      }
    }

    # Lengths of fish (year, age) by fleet.
    if(is.null(len_F_flt_mm)){
      len_F_flt_mm <- size.age.fishery[grepl("^[a-zA_Z]*.*len",names(size.age.fishery))]
      nm_unit_len <- c("mm")
      len_F_flt_mm_type <- gsub(".*\\<([a-zA-Z]*len[a-zA-Z]*)\\>.*","\\1",names(len_F_flt_mm))
      len_F_flt_mm_units <- gsub(paste0(".*\\<([:alpha:]*",paste(nm_unit_len,collapse = "|"),"[:alpha:]*)\\>.*"),
                              "\\1",names(len_F_flt_mm))

      # Reformat names to len.fleetType.fleet
      # (e.g. len.L.cHL for length of landings in the commercial hook and line fleet)
      len_F_flt_abb <- local({
        a <- gsub(paste0(c(unique(len_F_flt_mm_type),unique(len_F_flt_mm_units)),collapse="|"),"",names(len_F_flt_mm))
        b <- gsub("^\\.*|\\.*$","",a)
        c <- gsub("(.*?)(?<=.)(\\.D)$","D.\\1",b,perl=TRUE) # Change .D suffix to D. prefix
        d <- gsub("^([a-zA-Z]{2})(D)","D.\\1\\2",c)         # Add D. prefix if third letter is D
        e <- gsub("^([cr])(.*)","L.\\1\\2",d) # Add L. prefix to comps starting with c or r (landings not discards. won't modify fishery independent survey comps)
        e

      })


      len_F_flt_mm_units <- gsub("^([a-zA-Z]+)(.)(.*?)(.)([a-zA-Z]+)$","\\1.\\5",names(len_F_flt_mm),perl=TRUE)
      names(len_F_flt_mm) <- names(len_F_flt_mm_units) <- paste("len", #len_F_flt_mm_units,
                                                          len_F_flt_abb,sep=".")
    }


    # # Recompute landings by fleet to compare with calculations used in projections
    # Cn2 <- lapply(Cn,function(x){x*NA})
    # Cw2 <- lapply(Cw,function(x){x*NA})
    # for(i in 1:length(yb)){
    # for(j in 1:ncol(Fsum_flt)){
    #   Cn2[[j]][i,] <- L_calc(F_flt[[j]][i,], Z[i,], N[i,])
    #   Cw2[[j]][i,] <- L_calc(F_flt[[j]][i,], Z[i,], N[i,], wgt_F_flt_klb[[j]][i])
    # }
    # }

    # CPUE
    #   From bam tpl for SEDAR53:
    #     N_cHL(iyear)=elem_prod(elem_prod(Nmdyr(iyear),sel_cHL(iyear)),wholewgt_cHL_klb(iyear));
    #     pred_cHL_cpue(iyear)=q_cHL(iyear)*q_rate_fcn_cHL(iyear)*q_DD_fcn(iyear)*sum(N_cHL(iyear));
    U_a <- list()
    NU <- list()
    unit_U <- list()
    U_pr <- t.series[paste(yb),grepl("^U.*pr$",names(t.series)),drop=FALSE]
    U_ob <- t.series[paste(yb),grepl("^U.*ob$",names(t.series)),drop=FALSE]
    cv_U <- t.series[paste(yb),grepl("^cv.U",names(t.series)),drop=FALSE] # cvs associated with U during the base years
    # names(U_pr) <- gsub("^(U.)(.*)(.)(pr)$","\\1\\4.\\2",names(U_pr))
    U_abb <- gsub("^U.|.pr$","",names(U_pr))

    names(U_pr) <- names(U_ob) <- U_abb
    # Find q values in rdat
    # Are the q values provided in t.series? (2025-09-08 NPK This should be standard in rdat. All q should just be in t.series.)
    #nm_q <- paste0("q.",U_abb)
    nm_q_tsY <- names(rdat$t.series)[grepl(paste0("^q..*(",paste(U_abb,collapse="|"),")$"),names(rdat$t.series))] # names of q in t.series

    U_abb_tsY <- U_abb[unlist(lapply(U_abb,function(x){any(grepl(x,nm_q_tsY))}))]
    U_abb_tsN <- U_abb[!U_abb%in%U_abb_tsY]
    # nm_q_tsN <- nm_q[!nm_q%in%names(t.series)]
    #nm_q_tsNpcY <- names(parm.cons)[grepl(paste0("^log.q..*(",paste(U_abb,collapse="|"),")$"),names(parm.cons))]#paste0("log.",nm_q_tsN)[paste0("log.",nm_q_tsN)%in%names(parm.cons)]
    # nm_q_tsNpcN <- nm_q[which((!nm_q%in%names(t.series))&(!paste0("log.",nm_q)%in%names(parm.cons)))]
    if(length(nm_q_tsY)>0){
      q_mn <- t.series[paste(yb),nm_q_tsY,drop=FALSE]  # This is really a kind of mean q, since q_ratemult and q_DD_mult might scale it to compute U
    }
    # If any of the q values are not in t.series, look in parms.cons
    if(length(U_abb_tsN)>0){
      message(paste0("q values for ",paste(U_abb_tsN,collapse=", "), " index not found in t.series"))
      nm_q_tsNpcY <- names(parm.cons)[grepl(paste0("^log.q..*(",paste(U_abb_tsN,collapse="|"),")$"),names(parm.cons))]
      if(length(nm_q_tsNpcY)>0){
        message(paste0("time invariant values ",paste(nm_q_tsNpcY,collapse=", "), " found in parm.cons will be exp transformed and used instead."))
        U_abb_tsNpcY <- gsub(paste0("^(.*)(",paste(U_abb,collapse="|"),")"),"\\2",nm_q_tsNpcY)
        q_tsNpcY <- setNames(exp(parm.cons[nm_q_tsNpcY][8,]),U_abb_tsNpcY)
        # names(q_tsNpcY) <- gsub("^log.q.","",names(q_tsNpcY))
        for(nm_i in names(q_tsNpcY)){
          a <- U_pr[,nm_i]
          a[!is.na(a)] <- q_tsNpcY[[nm_i]]
          q_mn[,paste0("q.",nm_i)] <- a
          }
      }else{
        U_abb_tsNpcY <- NULL
      }
      # If any q not in t.series or parm.cons
      U_abb_tsNpcN <- U_abb[!U_abb%in%c(U_abb_tsY,U_abb_tsNpcY)]
      if(length(U_abb_tsNpcN)>0){
        warning(paste0("q values for ",paste(U_abb_tsNpcN,collapse=", "), " not found in t.series or parm.cons. Can't compute these cpue indices in the projections."))
      }
    }

    q_DD_mult <- t.series[paste(yb),"q.DD.mult",drop=FALSE] # density dependent function as a multiple of q (scaled a la Katsukawa and Matsuda. 2003)
    q_ratemult <- q_mn*0+1
    q <- q_mn*NA # Initialize values for final q value (q_mn * q_ratemult * q_DD_mult) multiplied by N to compute CPUE

    sel_U <- sel_age_2[which(gsub("sel.[vm].","",names(sel_age_2))%in%U_abb)]
    names(sel_U) <- paste0("sel.U.",U_abb)

    for(i in 1:length(sel_U)){
      if(is.vector(sel_U[[i]])){
        sel_U[[i]] <- matrix(sel_U[[i]], nrow=length(yb),ncol=length(sel_U[[i]]),byrow=TRUE,dimnames=list(yb,names(sel_U[[i]])))
      }
      sel_U_i <- sel_U[[i]]
      sel_nm_i <- names(sel_U)[i]
      sel_U_abb_i <- gsub("sel.U.","",sel_nm_i)

      nm_q_i <- paste0("q.",sel_U_abb_i)
      q_mn_i <- q_mn[paste(yb),nm_q_i,drop=FALSE]
      q_ratemult_nm_i <- paste0("q.",sel_U_abb_i,".rate.mult")
      if(q_ratemult_nm_i%in%names(t.series)){
        q_ratemult_i <- t.series[paste(yb),q_ratemult_nm_i,drop=FALSE]
      }else{
        q_ratemult_i <- q_mn_i*0+1
        message(paste("message:",q_ratemult_nm_i,"not found in t.series. Setting q rate multiplier values to 1.\n"))
      }
      q_ratemult[,nm_q_i] <- q_ratemult_i

      q_i <- q_mn_i*q_ratemult_i*q_DD_mult
      q[,nm_q_i] <- q_i

      unit_U_n_i <- local({ # Set default unit = 1 for keeping the index in numbers
        a <- sel_U_i
        a*0+1
      })
      # Since it can be hard to know if the index was in units of weight or numbers
      # which affects q, compute it both ways and then see which matches the
      # values in the base model better, then change the appropriate values.

      # If the U type is commercial or recreational (not a fishery independent survey), get unit_U_w_i from wgt_F_flt_klb..
      # WARNING: IN SOME CASES OF DISCARD INDICES THIS SHOULD BE LOOKING FOR WEIGHTS IN wgt.D!!! CODE IT NIKOLAI!!
      if(gsub("^(.{1})(.*)","\\1",sel_U_abb_i)%in%c("r","c")){
        unit_U_nm_i <- paste0("wgt.L.",sel_U_abb_i)
        if(unit_U_nm_i%in%names(wgt_F_flt_klb)){
          unit_U_w_i <- wgt_F_flt_klb[[paste0("wgt.L.",sel_U_abb_i)]]
        }else{
          message(paste("I could not find",unit_U_nm_i, "to multiply by", sel_nm_i,"when computing", paste0("N_",sel_U_abb_i,"."), paste0("U_pr_",sel_U_abb_i), "will be computed in numbers instead of weight.\n"))
        }
      }else{
        # ..otherwise use population weights
        unit_U_w_i <- wgt_b_klb
      }

      NUn_i <- Nmdyr*sel_U_i*unit_U_n_i
      NUw_i <- Nmdyr*sel_U_i*unit_U_w_i

      # U_a_i <- unlist(q_i)*NU_i
      U_a_n_i <- unlist(q_i)*NUn_i
      U_a_w_i <- unlist(q_i)*NUw_i
      U_n_i <- rowSums(U_a_n_i)
      U_w_i <- rowSums(U_a_w_i)
      U_ob_i <- U_ob[,sel_U_abb_i]
      U_pr_i <- U_pr[,sel_U_abb_i]

      # Compute the sums of squared deviations to determine whether the indices are
      # in numbers or in weight. Indices in the appropriate units are then added
      # to the the appropriate objects
      SSUi <- unlist(lapply(list(U_n_i=U_n_i,U_w_i=U_w_i),function(x){sum((U_pr_i-x)^2,na.rm=TRUE)}))
      if(names(SSUi)[which.min(SSUi)]=="U_n_i"){
        message(paste("The",sel_U_abb_i,"index appears to be in numbers in the base years and will therefore be projected in numbers."))
        NU_i <- NUn_i
        U_a_i <- U_a_n_i
        unit_U_i <- unit_U_n_i
      }else{
        message(paste("The",sel_U_abb_i,"index appears to be in weight in the base years and will therefore be projected in weight."))
        NU_i <- NUw_i
        U_a_i <- U_a_w_i
        unit_U_i <- unit_U_w_i
      }

      NU[[i]] <- NU_i
      U_a[[i]] <- U_a_i
      unit_U[[i]] <- unit_U_i
    }
    names(U_a) <- names(unit_U) <- names(NU) <- U_abb

    U <- as.data.frame(lapply(U_a,rowSums))

    # Identify other selectivities (e.g. associated with comps or aggregate landings or discards)
    # and compute associated numbers at age matrices
    Nmisc <- list()
    unit_misc <- list()
    sel_misc <- sel_age_2[!gsub("sel.[vm].","",names(sel_age_2))%in%unique(c(U_abb,gsub("sel.[LD].","",names(sel_F_flt))))]
    if(length(sel_misc)>0){
      message(paste0("selectivities found in sel_age_2 that don't clearly match a source of F or an index: ", paste(names(sel_misc),collapse=", ")))

    sel_misc_abb <- gsub("sel.[mv].","",names(sel_misc))
    names(sel_misc) <- paste0("sel.misc.",sel_misc_abb)

    for(i in 1:length(sel_misc)){
      if(is.vector(sel_misc[[i]])){
        sel_misc[[i]] <- matrix(sel_misc[[i]], nrow=length(yb),ncol=length(sel_misc[[i]]),byrow=TRUE,dimnames=list(yb,names(sel_misc[[i]])))
      }
      sel_misc_i <- sel_misc[[i]]
      sel_misc_nm_i <- names(sel_misc)[i]
      sel_misc_abb_i <- gsub("sel.[mv].","",sel_misc_nm_i)

      unit_misc_n_i <- local({ # Set default unit = 1 for computing numbers-at-age
        a <- sel_misc_i
        a*0+1
      })
      unit_misc_i <- unit_misc_n_i

      Nmisc_n_i <- Nmdyr*sel_misc_i*unit_misc_n_i
      Nmisc_i <- Nmisc_n_i

      Nmisc[[i]] <- Nmisc_i
      unit_misc[[i]] <- unit_misc_i
    }
    names(unit_misc) <- names(Nmisc) <- sel_misc_abb
}

    ## age compositions from the base years
    # observed
    nm_acomp_ob_comp_mats <- names(comp.mats)[grepl("^acomp.*.ob$",names(comp.mats))]
    acomp_ob_b <- comp.mats[nm_acomp_ob_comp_mats]
    names(acomp_ob_b) <- local({
      a <- gsub("^(acomp.)(.*)(.ob)$","\\2",names(acomp_ob_b))
      b <- gsub("(.*)(.D)$","D.\\1",a) # Convert .D suffix to D. prefix
      c <- gsub("^([cr])(.*)","L.\\1\\2",b) # Add L. prefix to comps starting with c or r (landings not discards. won't modify fishery independent survey comps)
      c
    })
    key_nm_acomp_ob <- setNames(names(acomp_ob_b),nm_acomp_ob_comp_mats)

    # predicted
    nm_acomp_pr_comp_mats <- names(comp.mats)[grepl("^acomp.*.pr$",names(comp.mats))]
    acomp_pr_b <- comp.mats[nm_acomp_pr_comp_mats]
    names(acomp_pr_b) <- local({
      a <- gsub("^(acomp.)(.*)(.pr)$","\\2",names(acomp_pr_b))
      b <- gsub("(.*)(.D)$","D.\\1",a) # Convert .D suffix to D. prefix
      c <- gsub("^([cr])(.*)","L.\\1\\2",b) # Add L. prefix to comps starting with c or r (landings not discards. won't modify fishery independent survey comps)
      c
    })
    key_nm_acomp_pr <- setNames(names(acomp_pr_b),nm_acomp_pr_comp_mats)

    ## length compositions from the base years
    # observed
    nm_lcomp_ob_comp_mats <- names(comp.mats)[grepl("^lcomp.*.ob$",names(comp.mats))]
    lcomp_ob_b <- comp.mats[nm_lcomp_ob_comp_mats]
    names(lcomp_ob_b) <- local({
      a <- gsub("^(lcomp.)(.*)(.ob)$","\\2",names(lcomp_ob_b))
      b <- gsub("(.*)(\\.D)$","D.\\1",a) # Convert .D suffix to D. prefix
      c <- gsub("^([a-zA-Z]{2})(D)","D.\\1\\2",b) # Add D. prefix if third letter is D
      d <- gsub("^([cr])(.*)","L.\\1\\2",c) # Add L. prefix to comps starting with c or r (landings not discards. won't modify fishery independent survey comps)
      d
    })
    key_nm_lcomp_ob <- setNames(names(lcomp_ob_b),nm_lcomp_ob_comp_mats)

    # predicted
    nm_lcomp_pr_comp_mats <- names(comp.mats)[grepl("^lcomp.*.pr$",names(comp.mats))]
    lcomp_pr_b <- comp.mats[nm_lcomp_pr_comp_mats]
    names(lcomp_pr_b) <- local({
      a <- gsub("^(lcomp.)(.*)(.pr)$","\\2",names(lcomp_pr_b))
      b <- gsub("(.*)(\\.D)$","D.\\1",a) # Convert .D suffix to D. prefix
      c <- gsub("^([a-zA-Z]{2})(D)","D.\\1\\2",b) # Add D. prefix if third letter is D
      d <- gsub("^([cr])(.*)","L.\\1\\2",c) # Add L. prefix to comps starting with c or r (landings not discards. won't modify fishery independent survey comps)
      d
    })
    key_nm_lcomp_pr <- setNames(names(lcomp_pr_b),nm_lcomp_pr_comp_mats)

    lenprob <- list()

    ## Build age-length conversion matrix associated with each set of length comps
    ## (will often be the same for all comps)
    if("all"%in%names(key_lenprob)){
      rgx_len_cv <- key_lenprob[["all"]] # pattern to use when searching names of parm.cons objects
      message(paste("keyword 'all' found in names of key_lenprob. ",rgx_len_cv, "will be used to search names of parm.cons to identify the cv of length-at-age to use when generating age-length converstion matrices."))
    }else{
      rgx_len_cv <- NULL
    }
      avail_len_cv <- names(parm.cons)[grepl("^len.cv",names(parm.cons))]
      message(paste("The following len.cv parameters were found in parm.cons:",paste(avail_len_cv,collapse=", ")))

    for(i in 1:length(lcomp_pr_b)){
      nmi <- names(lcomp_pr_b)[i]
      if("all"%in%names(key_lenprob)){
        rgx_len_cv_i <- rgx_len_cv
        nm_len_cv_i <- names(parm.cons)[grepl(rgx_len_cv_i,names(parm.cons))][1]
      }else if(nmi%in%names(key_lenprob)){
        rgx_len_cv_i <- key_lenprob[[nmi]]
        nm_len_cv_i <- names(parm.cons)[grepl(rgx_len_cv_i,names(parm.cons))][1]
      }else{
        rgx_len_cv_i <- "^len.cv"
        nm_len_cv_i <- names(parm.cons)[grepl(rgx_len_cv_i,names(parm.cons))][1]
        warning(paste(nmi, "not found in names(key_lenprob). Will use estimated value of",nm_len_cv_i,"from parm.cons.") )
      }

      if(is.na(nm_len_cv_i)){
        nm_len_cv_i <- avail_len_cv[1]
        warning(paste("Couldn't match",rgx_len_cv_i, "in names(parm.cons). Using",avail_len_cv[1],"instead for.",nmi,"length comps."))
      }

      # if(length(avail_len_cv)>1){
      #   message(paste("found multiple len.cv.val in rdat:",paste(avail_len_cv,collapse=", "),"Currently only using len.cv.val"))
      # }

      len_cv_i <- parm.cons[8,nm_len_cv_i]
      len_sd_i <- len*len_cv_i

      lc_i <- lcomp_pr_b[[i]]
      lc_bm_i <- as.numeric(colnames(lc_i))  # bin mid point
      lc_bw_i <- median(diff(lc_bm_i)) # bin width
      lc_bl_i <- lc_bm_i-lc_bw_i/2  # bin lo (minimum)
      lenprob_i <- local({
        mat_i <- matrix(NA,nrow=length(ages),ncol=length(lc_bm_i),dimnames = list("age"=ages,"lenbin"=lc_bm_i))
        mat_i[,1] <- pnorm((lc_bl_i[2]-len)/len_sd_i)
        for(i in 2:(ncol(mat_i)-1)){
          mat_i[,i] <- pnorm((lc_bl_i[i+1]-len)/len_sd_i)-pnorm((lc_bl_i[i]-len)/len_sd_i)
        }
        mat_i[,ncol(mat_i)] <- 1-rowSums(mat_i[,1:(ncol(mat_i)-1)])
        mat_i
      })
      lenprob[[i]] <- lenprob_i
    }

    names(lenprob) <- names(lcomp_pr_b)


    ## Recompute predicted age and length comps
    # Note: These should be the same or very similar to age and length comps computed by BAM.
    #  But these are recomputed for the base years as a check since the computation
    #  will be used in the projection years
    # acomp_b_pr_2 <-  list()
    # for(nm_i in names(acomp_pr_b)){
    #   yrs_ac_i <- rownames(acomp_pr_b[[nm_i]])
    #   if(grepl("^[LD]",nm_i)){ # If length comps are associated with catch (landings or discards), used numbers from landings
    #     n_i <- Cn[[paste0("Cn.",nm_i)]][yrs_ac_i,]
    #   }else{
    #     n_i <- NU[[nm_i]][yrs_ac_i,] # If not (i.e. it's a fishery independent survey) use numbers associated with the index
    #   }
    #   P_n_i <- t(apply(n_i,1,function(x){x/sum(x)}))
    #   acomp_b_pr_2_i <- t(apply(P_n_i,1,function(x){colSums(x*age_error)}))
    #   acomp_b_pr_2[[nm_i]] <- acomp_b_pr_2_i
    # }

    # lcomp_b_pr_2 <-  list()
    # for(nm_i in names(lcomp_pr_b)){
    #   yrs_lc_i <- rownames(lcomp_pr_b[[nm_i]])
    #   lenprob_i <- lenprob[[nm_i]]
    #   if(grepl("^[LD]",nm_i)){ # If length comps are associated with catch (landings or discards), used numbers from landings
    #     n_i <- Cn[[paste0("Cn.",nm_i)]][yrs_lc_i,]
    #   }else{
    #     n_i <- NU[[nm_i]][yrs_lc_i,] # If not (i.e. it's a fishery independent survey) use numbers associated with the index
    #   }
    #   P_n_i <- t(apply(n_i,1,function(x){x/sum(x)}))
    #   lcomp_b_pr_2_i <- t(apply(P_n_i,1,function(x){colSums(x*lenprob_i)}))
    #   lcomp_b_pr_2[[nm_i]] <- lcomp_2_i
    # }

    ##
    # number of trips in compositions
    nm_ntrip_ts <- names(t.series)[grepl("^([la]comp).*(.n)$",names(t.series))] # Names of ntrip columns in t.series
    ntrip <- t.series[paste(yb),nm_ntrip_ts]

    # ntrip[ntrip==-99999] <- NA
    names(ntrip) <- local({
      a <- gsub("(.n)$","",names(ntrip))
      b <- gsub("^([al]comp.)(.*)(.)(D)$","\\1D.\\2",a) # Convert .D suffix to D. prefix
      c <- gsub("^([al]comp.)([cr])(.*)","\\1L.\\2\\3",b) # Add L. prefix to comps starting with c or r (landings not discards. won't modify fishery independent survey comps)
      c
    })
    key_nm_ntrip <- setNames(names(ntrip),nm_ntrip_ts)

    # effective number of trips in compositions
    nm_neff_ts <- names(t.series)[grepl("^([la]comp).*(.neff)$",names(t.series))] # Names of neff columns in t.series
    neff <- t.series[paste(yb),nm_neff_ts]
    names(neff) <- names(ntrip)
    key_nm_neff <- setNames(names(neff),nm_neff_ts)

    # number of fish in compositions
    nm_nfish_ts <- names(t.series)[grepl("^([la]comp).*(.nfish)$",names(t.series))]  # Names of nfish columns in t.series
    nfish <- t.series[paste(yb),nm_nfish_ts]
    # nfish[nfish==-99999] <- NA
    names(nfish) <- local({
      a <- gsub("(.nfish)$","",names(nfish))
      b <- gsub("^([al]comp.)(.*)(.)(D)$","\\1D.\\2",a) # Convert .D suffix to D. prefix
      c <- gsub("^([al]comp.)([cr])(.*)","\\1L.\\2\\3",b) # Add L. prefix to comps starting with c or r (landings not discards. won't modify fishery independent survey comps)
      c
    })
    key_nm_nfish <- setNames(names(nfish),nm_nfish_ts)

    #Initial conditions
    #Recruits in yr 1 constrained to S-R curve, or else varies in stochastic runs
    if(is.null(N_styr_proj)){N_styr_proj <- tail(N.age,1)}
    N_endyr <- N.age[paste(endyr),]
    if(is.null(S_styr_proj)){
      S_styr_proj <- local({
        Z_styr_proj <- Z.age[paste(endyr),]
        N_endyr_spn <- N_endyr*exp(-1.0*Z_styr_proj*spawn_time)
        sum(N_endyr_spn*reprod)
      })
    }

    ## SR stuff. Kind of long but used to accommodate varying naming conventions.
    # Note NK 2023-10-15: I think that biascorr and R.autocorr should be available
    # in all models, but maybe with a different name.

    biascorr <-  if("BH.biascorr"%in%nm_parms){parms$BH.biascorr
    }else if("biascorr"%in%nm_parms){
      nm_biascorr <- "biascorr"
      message(paste("BH.biascorr not found in parms.",nm_biascorr,"found and will be used instead."))
      parms[[nm_biascorr]]
    }else if(any(grepl("biascorr",nm_parms))){
      nm_biascorr <- nm_parms[which(grepl("biascorr",nm_parms))][1]
      message(paste("BH.biascorr not found in parms.",nm_biascorr,"found and will be used instead."))
      parms[[nm_biascorr]]
    }else{
      message(paste("Can't match biascorr in parms. Setting biascorr = 1."))
      1
    }

    R.autocorr <-  if("R.autocorr"%in%nm_parms){parms$R.autocorr
    }else{
      message(paste("R.autocorr not found in parms. Setting R.autocorr = 0."))
      0
    }

    if(SR_method=="BH"){
      if(is.null(SR_par)){
        h <-  if("BH.steep"%in%nm_parms){parms$BH.steep
        }else if("steep"%in%nm_parms){
          nm_steep <- "steep"
          message(paste("BH.steep not found in parms.",nm_steep,"found and will be used instead."))
          parms[[nm_steep]]
        }else if(any(grepl("steep",nm_parms))){
          nm_steep <- nm_parms[which(grepl("steep",nm_parms))][1]
          message(paste("BH.steep not found in parms.",nm_steep,"found and will be used instead."))
          parms[[nm_steep]]
        }else{
          message(paste("Can't match steep in parms."))
          NULL
        }
        R0 <-  if("BH.R0"%in%nm_parms){parms$BH.R0
        }else if("R0"%in%nm_parms){
          nm_R0 <- "R0"
          message(paste("BH.R0 not found in parms.",nm_R0,"found and will be used instead."))
          parms[[nm_R0]]
        }else if(any(grepl("R0",nm_parms))){
          nm_R0 <- nm_parms[which(grepl("R0",nm_parms))][1]
          message(paste("BH.R0 not found in parms.",nm_R0,"found and will be used instead."))
          parms[[nm_R0]]
        }else{
          message(paste("Can't match R0 in parms."))
          NULL
        }
        Phi0 <-  if("BH.Phi0"%in%nm_parms){parms$BH.Phi0
        }else if("Phi0"%in%nm_parms){
          nm_Phi0 <- "Phi0"
          message(paste("BH.Phi0 not found in parms.",nm_Phi0,"found and will be used instead."))
          parms[[nm_Phi0]]
        }else if(any(grepl("Phi0",nm_parms))){
          nm_Phi0 <- nm_parms[which(grepl("Phi0",nm_parms))][1]
          message(paste("BH.Phi0 not found in parms.",nm_Phi0,"found and will be used instead."))
          parms[[nm_Phi0]]
        }else{
          message(paste("Can't match Phi0 in parms."))
          NULL
        }
      }else{
        for(i in 1:length(SR_par)){xi <- SR_par[i]; assign(x=names(xi),value=xi)}
      }
      SR_par_null <- c("h","R0","Phi0")[c(is.null(h),is.null(R0),is.null(Phi0))]
      # If not all of the BH SR parameters are found, resort to using the GM method
      if(length(SR_par_null)>0){
        message(paste("BH SR parameters not found:",paste(SR_par_null,collapse=", ")," Can't use SR_method = BH.\n"))
        if(!is.null(R0)){
        message(paste("R0 parameter found. Using SR_method = R0.\n"))
        SR_method <- "R0"
        }else{
          message(paste("Using SR_method = GM.\n"))
          SR_method <- "GM"
        }
      }else{
        message("all BH SR parameters found. Using SR_method = BH.\n")
      }
    }else if(SR_method=="R0"){
      R0 <-  if("BH.R0"%in%nm_parms){parms$BH.R0
      }else if("R0"%in%nm_parms){
        nm_R0 <- "R0"
        message(paste("BH.R0 not found in parms.",nm_R0,"found and will be used instead."))
        parms[[nm_R0]]
      }else if(any(grepl("R0",nm_parms))){
        nm_R0 <- nm_parms[which(grepl("R0",nm_parms))][1]
        message(paste("BH.R0 not found in parms.",nm_R0,"found and will be used instead."))
        parms[[nm_R0]]
      }else{
        message(paste("Can't match R0 in parms."))
        NULL
      }
      if(!any(is.null(c(R0,biascorr)))){
        message(paste("Found values of R0 and biascorr. Using SR_method = R0.\n"))
      }

    }else if(SR_method=="GM"){
      message(paste("Using SR_method = GM.\n"))
    }else{
      message(paste("SR_method value not valid. Using default SR_method = GM.\n"))
      SR_method <- "GM"
    }

    if(is.null(M)){
      M <- setNames(a.series$M,rownames(a.series))
    }

  # } # end if(!is.null(rdat))


  ##############################################################################
  ###### Setup projections #####################################################
  ##############################################################################
    if(nyp>0){ # Only run this if there is actually a projection
  ## Build empty objects
  # generic empty objects
    ema <-  setNames(rep(NA,nages),ages) # age vector; a = ages
    emyp <- setNames(rep(NA,nyp),yp)     # year vector; yp = years of the projection period
    emypa <- matrix(NA,nrow=nyp,ncol=nages,dimnames=list("year"=yp,"age"=ages)) # matrix (year,age) during the projection period
    emfyp <- matrix(NA,nrow=nyp,ncol=ncol(Fsum_flt),dimnames=list(yp,colnames(Fsum_flt)))
    emuyp <- matrix(NA,nrow=nyp,ncol=length(U_a),dimnames=list("year"=yp,"U"=names(U_a))) # matrix (year, U) during the projection period

    emfypa <- local({
      a <- lapply(1:length(Fsum_flt),function(x){emypa})
      names(a) <- names(Fsum_flt)
      a
    })
    emuypa <- local({
      a <- lapply(1:length(U_a),function(x){emypa})
      names(a) <- names(U_a)
      a
    })

  Z_p   <-      emypa  # total mortality by age
  F_L_p <-      emypa  # fishing mortality of landings by age
  F_D_p <-      emypa  # fishing mortality of discards by age
  Nspwn_p <- emypa  # numbers at age at spawn_time
  Nmdyr_p <- emypa  # numbers at age at midyear
  N_p <-      emypa  # numbers at age by year
  NU_p <-    emuypa # numbers (or biomass) at age for each index of abundance, used in U_a calculations
  if(length(sel_misc)>0){
    Nmisc_p <- lapply(Nmisc,function(x){matrix(NA,nrow=nyp,ncol=ncol(x),dimnames=list("year"=yp,"age"=colnames(x)))})
  }else{
    Nmisc_p <- NULL
  }

  S_p <-    emyp # spawning stock (often biomass in mt, sometimes eggs in n)
  B_p <-    emyp # population biomass (mt)
  R_p <-    emyp # recruits (n)
  logRdev_p <- emyp # log recruitment deviation
  Ffull_p <- emyp # Max F across by year, across all fleets (/yr) (comment from bam: Max across ages, fishing mortality rate by year (may differ from Fsum bc of dome-shaped sel)

  Ln_p <- emyp	# total landings (n)
  Lw_p <- emyp	# total landings (weight)
  Dn_p <- emyp	# total dead discards (n)
  Dw_p <- emyp	# total dead discards (weight)

  ## New objects
  # sel_wgted_F_flt_p <- lapply(1:length(Fsum_flt),function(x){emypa})
  # names(sel_wgted_F_flt_p) <- gsub("^F.","sel.",names(Fsum_flt))

  U_a_p <-  emuypa       # predicted cpue at age by fleet
  unit_U_p <- emuypa     # weights associated with cpue at age by fleet
  q_p <- emuyp           # catchability (q; final value used to scale U; q = q_mn*q_ratemult*q_DD)
  q_mn_p <- emuyp        # catchability (mean value usually presented in t.series)
  q_ratemult_p <- emuyp  # catchability rate multiplier
  q_DD_mult_p <- emyp

  F_flt_p <- emfypa
  Fsum_flt_p <-  emfyp

  wgt_F_flt_p <- emfypa
  names(wgt_F_flt_klb) <- gsub("^F","wgt",names(emfypa))

  len_F_flt_p <- emfypa
  names(len_F_flt_mm) <- gsub("^F","len",names(emfypa))

  Cn_p <- emfypa # Catch in numbers by fleet (landings or discards)
  names(Cn_p) <- gsub("^F","Cn",names(emfypa))

  Cw_p <- emfypa # Catch in weight by fleet (landings or discards)
  names(Cw_p) <- gsub("^F","Cw",names(emfypa))

  ## Initialization
  S_p[1] <- S_styr_proj

  N_p[1,] <- N_styr_proj
  if(stochastic){
    # Compute rec devs (independent of SR_method)
    logRdev_endyr <- rdat$t.series[paste(rdat$parms$endyr),"logR.dev"]
      logRdev_p[1] <- rnorm(n=1 ,mean=logRdev_endyr*R.autocorr,    sd=sqrt(2.0*log(biascorr)))
    for(i in 2:nyp){
      logRdev_p[i] <- rnorm(n=1 ,mean=logRdev_p[i-1]*R.autocorr, sd=sqrt(2.0*log(biascorr)))
    }
    if(SR_method=="BH"){
      N_p[1,1] = exp(logRdev_p[1]) * (0.8*R0*h*S_p[1])/(0.2*R0*Phi0*(1.0-h)+ (h-0.2)*S_p[1])
    }else if(SR_method=="R0"){
      N_p[1,1] = exp(logRdev_p[1]) * R0
    }else{
      # Use default GM method
      N_p[1,1] = exp(logRdev_p[1]) * R_b_gm
    }
  }else{
    logRdev_p[paste(yp)] <- 0
  }
  B_p[1] <- sum(N_p[1,]*wgt_mt)

  # Project cvs for catch and cpue
  cv_C_p <- local({
    a <- matrix(apply(cv_C,2,
                      function(x){
                        y <- rep(bamExtras::geomean2(tail(x,nyb_rcn$L)),nyp)
                        if(is.na(tail(x,1))){
                          y <- y*NA
                        }
                        y
                      }
    ),
    nrow=nyp,dimnames=list(paste(yp),names(cv_C)))
    # rownames(a) <- paste(yp)
    a
  })

  cv_U_p <- local({
    a <- matrix(apply(cv_U,2,
                      function(x){
                        y <- rep(bamExtras::geomean2(tail(x,nyb_rcn$U)),nyp)
                        if(is.na(tail(x,1))){
                          y <- y*NA
                        }
                        y
                      }
                      ),
                nrow=nyp,dimnames=list(paste(yp),names(cv_U)))
    # rownames(a) <- paste(yp)
    a
  })

  ## sel during projection years (by py, age, fleet)
  # for any source of F (landings or discards)
  # Fill with weighted selectivity from last year of base (endyr)
  # NOTE: 2025-08-29 This isn't quite the selectivity, but actually fleet-specific selectivity multiplied by fleet-specific Fsum
  # divided by overall Ffull (i.e. (F_ybgm_flt_x*sel_endyr_x)/F_end_apex)
  sel_wgted_F_flt_p <- lapply(1:length(sel_wgted_F_flt),function(i){
    nm_i <- names(sel_wgted_F_flt)[i]
    sel_wgted_F_endyr_i <- sel_wgted_F_flt[,gsub("^sel.","F.",nm_i)]
    matrix(sel_wgted_F_endyr_i, nrow=nyp, ncol=nages, byrow=TRUE, dimnames = dimnames(emypa))
  })
  names(sel_wgted_F_flt_p) <- gsub("^F.","sel.",names(Fsum_flt))

  sel_F_flt_p <- lapply(1:length(sel_F_flt),function(i){
    nm_i <- names(sel_F_flt)[i]
    sel_F_endyr_i <- sel_F_flt[[nm_i]][paste(endyr),] # Gets unweighted selectivities
    matrix(sel_F_endyr_i, nrow=nyp, ncol=nages, byrow=TRUE, dimnames = dimnames(emypa))
  })
  names(sel_F_flt_p) <- gsub("^F.","sel.",names(Fsum_flt))

  # sel during projection years (by py, age, fleet) for any source of indices of abundance (cpue)
  # Fill with selectivity from last year of base (endyr)
  sel_U_p <- lapply(1:length(sel_U),function(i){
    nm_i <- names(sel_U)[i]
    sel_U_endyr_i <- sel_U[[nm_i]][paste(endyr),]
    matrix(sel_U_endyr_i, nrow=nyp, ncol=nages, byrow=TRUE, dimnames = dimnames(emypa))
  })
  names(sel_U_p) <- names(sel_U)

  # sel during projection years (by py, age, fleet) for any miscellaneous data sources
  # Fill with selectivity from last year of base (endyr)
  if(length(sel_misc)>0){
  sel_misc_p <- lapply(1:length(sel_misc),function(i){
    nm_i <- names(sel_misc)[i]
    sel_misc_endyr_i <- sel_misc[[nm_i]][paste(endyr),]
    matrix(sel_misc_endyr_i, nrow=nyp, ncol=nages, byrow=TRUE, dimnames = dimnames(emypa))
  })
  names(sel_misc_p) <- names(sel_misc)
  }

  # weight associated with cpue during the projection years
  # (mostly used to convert commercial cpue to weight. Otherwise set to 1 which leaves cpue in numbers)
  unit_U_p <- lapply(1:length(unit_U),function(i){
    nm_i <- names(unit_U)[i]
    x_i <- unit_U[[nm_i]]
    x_endyr_i <- x_i[paste(endyr),]
    matrix(x_endyr_i, nrow=nyp, ncol=nages, byrow=TRUE, dimnames = dimnames(emypa))
  })
  names(unit_U_p) <- names(unit_U)

  # catchability associated with cpue during the projection years
  q_mn_p[paste(yp),] <- rep(unlist(q_mn[paste(endyr),paste0("q.",colnames(q_mn_p))]),each=nyp)
  # Note: q_ratemult and q_DD_mult currently do nothing in the projections.
  q_ratemult_p[paste(yp),] <- rep(unlist(q_ratemult[paste(endyr),paste0("q.",colnames(q_ratemult_p))]),each=nyp)
  q_DD_mult_p[paste(yp)] <- rep(q_DD_mult[paste(endyr),],nyp)
  q_p[paste(yp),] <- q_mn_p*q_ratemult_p*q_DD_mult_p


  # compositions
  acomp_ob_p <- list()
  for(i in names(acomp_ob_b)){
  ac_i <- acomp_ob_b[[i]]
  acomp_ob_p[[i]] <- matrix(NA,nrow=nyp,ncol=ncol(ac_i),dimnames=list("year"=yp,"age"=colnames(ac_i)))
  }
  acomp_pr_p <- acomp_ob_p

  lcomp_ob_p <- list()
  for(i in names(lcomp_ob_b)){
    lc_i <- lcomp_ob_b[[i]]
    lcomp_ob_p[[i]] <- matrix(NA,nrow=nyp,ncol=ncol(lc_i),dimnames=list("year"=yp,"lenbin"=colnames(lc_i)))
  }
  lcomp_pr_p <- lcomp_ob_p

  ntrip_p <- local({
    a <- ceiling(apply(tail(ntrip,nyb_rcn$comp),2,function(x){bamExtras::geomean2(x)}))
    matrix(a,nrow=nyp,ncol=length(a),dimnames=list(year=paste(yp),fleet=names(a)),byrow=TRUE)
  })

  # Compute neff_p (kind of long, but I think this is the way to do it)
  #neff
  # 1. Get ntrip
  # 2. Get log.dm paramater that matches each column in ntrip
  # 3. Compute neff

  key_ntrip_log_dm <- local({
    #nm_neff_ts <- names(t.series)[grepl("neff$",names(t.series))]
    nm_log_dm <- names(parm.cons[grepl("log.dm",names(parm.cons))])

    # names of log_dm parameters that match comp neff names in t.series
    nm_ntrip_log_dm <- gsub(".n$","",nm_ntrip_ts) %>%
      gsub("lcomp","lenc",.) %>%
      gsub("acomp","agec",.) %>%
      paste0("log.dm.",.)

    setNames(nm_ntrip_log_dm,nm_ntrip_ts)
  })

  # log_dm values that match ntrip columns in t.series
  compute_neff <- function(n,log_dm){
    (1+n*exp(log_dm))/(1+exp(log_dm))
  }
  log_dm_ntrip <- setNames(unlist(parm.cons[8,key_ntrip_log_dm]),names(key_ntrip_log_dm))

  neff_p <- ntrip_p*NA
  for(i in 1:ncol(ntrip_p)){
    nm_i <- dimnames(ntrip_p)[[2]][i]
    ntrip_p_i <- ntrip_p[,i]
    log_dm_i <- log_dm_ntrip[names(key_nm_ntrip)[which(key_nm_ntrip==nm_i)]]
    neff_p[,i] <- compute_neff(ntrip_p_i,log_dm_i)
  }

  nfish_p <- local({
    a <- ceiling(apply(tail(nfish,nyb_rcn$comp),2,function(x){bamExtras::geomean2(x)}))
    matrix(a,nrow=nyp,ncol=length(a),dimnames=list(year=paste(yp),fleet=names(a)),byrow=TRUE)
  })

  ##############################################################################
  #### Projection loop #########################################################
  ##############################################################################
  ### years 1 to nyp-1 (i.e. not the last year)
  if(nyp>1){
    for (i in 1:nyp) {
      yp_i <- paste(yp[i])

      # Set fishing mortality
      if(i%in%1:nyp_cur){
        Ffull_p_i <- F_cur
      }else if(i%in%((nyp_cur+1):nyp)){
        Ffull_p_i <- F_proj
      }

      # Interim adjustment of F (optional)
      yrlim_p <- paste(endyr+c(1,nyp))
      yr_ia <- paste(seq(yrlim_p[1],yrlim_p[2],by=args_ia$yr_ia_by))
      if(!args_ia$include_yrlim_p){
        yr_ia <- yr_ia[!yr_ia%in%yrlim_p]
      }

      if(!is.null(args_ia$U_nm)&yp_i%in%yr_ia){
        F_adj <- F_adjust(data=U,
                          U_nm=args_ia$U_nm,
                          yr_ref = paste(endyr),
                          yr = paste(sort(as.numeric(yp_i)+args_ia$yr_U_lag)),
                          type = args_ia$type
                          )
        message(paste("For year =",yp_i,"full F was adjusted by", F_adj,"based on the",args_ia$U_nm, "index"))

      }else{
        F_adj <- 1
      }

      Ffull_p[i] <- Ffull_p_i*F_adj

      # Compute total Z and F at-age for year i
      # Z_p <-   outer(Ffull_p,M+sel_tot,FUN="*") # Z for population
      Z_p[i,] <-   M+Ffull_p[i]*sel_tot   # Z for population  #t(M+t(outer(Ffull_p,sel_tot,FUN="*")))
      F_L_p[i,] <- Ffull_p[i]*sel_L       # F for landings outer(Ffull_p,sel_L,  FUN="*")
      if(!is.null(sel_D)){
        F_D_p[i,] <-  Ffull_p[i]*sel_D # F for discards outer(Ffull_p,sel_D,  FUN="*")
      }else{
        F_D_p[i,] <- F_L_p[i,]*0
      }

      # F during projection years (by py, age, fleet) for any source of F (landings or discards)
      # (analogous to F_L_p or F_D_p but by fleet)
      # Note that:
      #      F.L.tmp <- Reduce("+",F_flt_p[grepl("^F.L.",names(F_flt_p))])
      #      F.L.tmp == F_L_p
      F_flt_p <- lapply(1:length(F_flt_p),function(j){
        sel_wgted_F_flt_p[[j]]*Ffull_p
      })
      names(F_flt_p) <- names(Fsum_flt)

      Fsum_flt_p_i <- local({
        a <- lapply(1:length(F_flt_p),function(j){
          #out <- (F_flt_p[[j]]/sel_F_flt_p[[j]])
          out <- (F_flt_p[[j]][yp_i,,drop=FALSE]/sel_F_flt_p[[j]][yp_i,,drop=FALSE])
          out[is.na(out)] <- 0 # Replace NA with zero. NA will occur when values in sel_F_flt_p are zero

          # Run this check in the last projection year
#          if(i==nyp){
            if(length(unique(round(as.numeric(out),6)))>1){
              warning("When dividing Ffull/selectivity matrices (age,year) to compute Fsum vectors (year), in projection years, results are not identical for all ages. This is likely a rounding error but may be something more concerning.")
            }
#          }
          out[,1] # The results of dividing an Ffull matrix by a selectivity matrix should yield a matrix out where all columns are identical
                  # since the Fsum vector (year) is multiplied by the selectivity matrix (age,year) to get the F
        })
        b <- unlist(a)
        names(b) <- names(Fsum_flt)
        b
      })
      Fsum_flt_p[yp_i,names(Fsum_flt_p_i)] <- Fsum_flt_p_i


      Fage_p <- Reduce("+",F_flt_p)
      Fsum_p <- rowSums(Fsum_flt_p)

      # wgt during projection years (by py, age, fleet) for any source of F (landings or discards)
      # Fill with fish weights from last year of base (endyr)
      wgt_F_flt_p <- lapply(1:length(F_flt_p),function(j){
        nm_i <- names(F_flt_p)[j]
        wgt_i <- wgt_F_flt_klb[[gsub("^F.","wgt.",nm_i)]]
        wgt_endyr_i <- wgt_i[paste(endyr),]
        matrix(wgt_endyr_i, nrow=nyp, ncol=nages, byrow=TRUE, dimnames = dimnames(emypa))
      })
      names(wgt_F_flt_p) <- gsub("^F.","wgt.",names(F_flt_p))

      # len during projection years (by py, age, fleet) for any source of F (landings or discards)
      # Fill with fish weights from last year of base (endyr)
      len_F_flt_p <- lapply(1:length(F_flt_p),function(j){
        nm_i <- names(F_flt_p)[j]
        len_i <- len_F_flt_mm[[gsub("^F.","len.",nm_i)]]
        len_endyr_i <- len_i[paste(endyr),]
        matrix(len_endyr_i, nrow=nyp, ncol=nages, byrow=TRUE, dimnames = dimnames(emypa))
      })
      names(len_F_flt_p) <- gsub("^F.","len.",names(F_flt_p))

      ## Population (year i)
      # N_p = number of fish at-age at the start of the year
      B_p[i] <-  sum(N_p[i,]*wgt_mt) # biomass of fish at-age at that start of the year
      Nmdyr_p[i,] <- N_p[i,]*(exp(-1.0*Z_p[i,]*0.5))        # number of fish at-age at mid-year
      Nspwn_p[i,] <- N_p[i,]*(exp(-1.0*Z_p[i,]*spawn_time)) # number of fish at-age at spawning time
      S_p[i] <- sum(Nspwn_p[i,]*reprod)                     # spawning stock at-age at spawning time

      ## cpue
      # By fleet
      for(j in 1:length(U_a_p)){
        nm_j <- names(U_a_p)[j]
        NU_p_ij <- Nmdyr_p[yp_i,]*sel_U_p[[paste0("sel.U.",nm_j)]][yp_i,]*unit_U_p[[nm_j]][yp_i,]
        NU_p[[nm_j]][yp_i,] <- NU_p_ij
        U_a_p_ij <- NU_p_ij*q_p[yp_i,nm_j]
        # # Add observation error to index (Note: should only be included when generating BAM dat file 2025-09-15)
        # U_p_ij <- local({
        #   a <- sum(U_a_p_ij)
        #   cv_nm_j <- paste0("cv.U.",nm_j)
        #   cv_U_p_ij <- cv_U_p[yp_i,cv_nm_j]*model_obs$cv_U_sc
        #   as.numeric(lnorm_vector_boot(a,cv_U_p_ij))
        # })
        # U_a_p[[nm_j]][yp_i,] <- (U_a_p_ij/sum(U_a_p_ij))*U_p_ij
        U_a_p[[nm_j]][yp_i,] <- U_a_p_ij
      }

      ## Fishery (year i)
      # Aggregate (Standard calculations)
      Ln_p[i] <- sum(L_calc(F_L_p[i,], Z_p[i,], N_p[i,]))
      Lw_p[i] <- sum(L_calc(F_L_p[i,], Z_p[i,], N_p[i,], wgt_L_klb))
      if(any(!is.null(c(sel_D,wgt_D_klb)))){
        Dn_p[i] <- sum(L_calc(F_D_p[i,], Z_p[i,], N_p[i,]))
        Dw_p[i] <- sum(L_calc(F_D_p[i,], Z_p[i,], N_p[i,], wgt_D_klb))
      }
      # By fleet
      for(j in 1:ncol(Fsum_flt)){
        Cn_p[[j]][i,] <- L_calc(F_flt_p[[j]][i,], Z_p[i,], N_p[i,])
        Cw_p[[j]][i,] <- L_calc(F_flt_p[[j]][i,], Z_p[i,], N_p[i,], wgt_F_flt_p[[j]][i,])
      }

      ## Nmisc_p
      if(length(sel_misc)>0){
        for(j in 1:length(Nmisc_p)){
          nm_j <- names(Nmisc_p)[j]
          sel_misc_ij <- sel_misc_p[[paste0("sel.misc.",nm_j)]][yp_i,]
          unit_misc_n_ij <- local({ # Set default unit = 1 for computing numbers-at-age
            a <- sel_misc_ij
            a*0+1
          })
          unit_misc_n_ij <- unit_misc_n_ij

          Nmisc_p_n_ij <- Nmdyr_p[yp_i,]*sel_misc_ij*unit_misc_n_ij
          Nmisc_p_ij <- Nmisc_p_n_ij
          Nmisc_p[[nm_j]][yp_i,] <- Nmisc_p_ij
        }
      }

      # Update data structures
      U_p <- as.data.frame(lapply(U_a_p,rowSums)) # Compute cpue time series (summing across ages)
      # 2025-09-05 NPK not standardized here so that U=N*q. Should be standardized when projecting dat file.
      U <- rbind(U,U_p[yp_i,,drop=FALSE])

      if(i<nyp){
        ## Population (year i+1)
        # Calculate number of offspring at midyear in year i which will be recruits in first age class in year i+1
        if(SR_method=="BH"){ # Beverton-Holt stock-recruit relationship
          if(stochastic){
            N_p[i+1,1] <- exp(logRdev_p[i+1]) * (0.8*R0*h*S_p[i])/(0.2*R0*Phi0*(1.0-h)+ (h-0.2)*S_p[i])
          }else{
            N_p[i+1,1] <- biascorr            * (0.8*R0*h*S_p[i])/(0.2*R0*Phi0*(1.0-h)+ (h-0.2)*S_p[i])
          }
        }else if(SR_method=="R0"){ # model estimated mean recruitment
          if(stochastic){
            N_p[i+1,1] <- exp(logRdev_p[i+1]) * R0
          }else{
            N_p[i+1,1] <- biascorr            * R0
          }
        }else if(SR_method=="GM"){ # Geometric mean recruitment
          if(stochastic){
            N_p[i+1,1] <- exp(logRdev_p[i+1]) * R_b_gm
          }else{
            N_p[i+1,1] <- biascorr            * R_b_gm
          }
          N_p[i+1,1] <- R_b_gm
        }

        # Fish from year i either die or become 1 year older in year i+1
        N_p[(i+1),2:nages] <- N_p[i,1:(nages-1)]*(exp(-1.0*Z_p[i,1:(nages-1)]))
        N_p[(i+1),nages]   <- N_p[(i+1),nages] +
          N_p[i,nages]*(exp(-1.0*Z_p[i,nages])) #plus group
      }
    }  # end i
  }  # end if(nyp>1)

  ##############################################################################
  #### Post-projection computations ############################################
  ##############################################################################

  R_p <- N_p[,1] # Get recruitment vector from N-at-age matrix

  ## Project observed age and length comps during projection years.
  #  Note: Comps are projected without error based on BAM model predicted numbers-at-age for each fleet
  for(nm_i in names(acomp_ob_p)){
    yrs_acomp_ob_b_i <- rownames(acomp_ob_b[[nm_i]])
    yrs_acomp_ob_p_i <- rownames(acomp_ob_p[[nm_i]])
    # If the last year of comps is within the last nyb_rcn$comp of the assessment, then project them.
    if(max(as.numeric(yrs_acomp_ob_b_i))%in%tail(as.numeric(yb),nyb_rcn$comp)){
    if(grepl("^[LD]",nm_i)){ # If length comps are associated with catch (landings or discards), use numbers from landings
      # Remove nm_comp_sfx from nm_i, if specified
      nm_Cn_i <- nm_i
      if(nm_comp_sfx[1]!=""){
        pat_i <- paste0(nm_comp_sfx,collapse="|")
        nm_Cn_i <- gsub(paste0("(.*)(",pat_i,")$"),"\\1",nm_i)
      }
      n_i <- Cn_p[[paste0("Cn.",nm_Cn_i)]][yrs_acomp_ob_p_i,,drop=FALSE]
    }else if(nm_i%in%names(NU_p)){
      n_i <- NU_p[[nm_i]][yrs_acomp_ob_p_i,,drop=FALSE] # If not, see if it's associated with an index (e.g. fishery independent) and use numbers associated with the index
    }else if(nm_i%in%names(Nmisc_p)){
      n_i <- Nmisc_p[[nm_i]][yrs_acomp_ob_p_i,,drop=FALSE] # If not, there should be an N-at-age matrix in Nmisc that matches
    }else{
      warning(paste0(nm_i, " from acomp_ob_p doesn't appear to match any sources of F, indices, or other data sources associated with a selectivity."))
    }
    P_n_i <- t(apply(n_i,1,function(x){x/sum(x)}))
    agebins_agec_i <- dimnames(acomp_ob_b[[nm_i]])[[2]]
    # Make sure the agebins match the observed age comps and apply a plus group if necessary
    acomp_ob_p_i <- local({
      a <- t(apply(P_n_i,1,function(x){colSums(x*age_error)}))
      b <- a[,agebins_agec_i]
      bplusgroup <- tail(agebins_agec_i,1)
      b[,bplusgroup] <- b[,bplusgroup]+rowSums(a[,as.numeric(dimnames(a)[[2]])>as.numeric(bplusgroup)])
      b
    })
    acomp_ob_p[[nm_i]] <- acomp_ob_p_i
    }
  }

  for(nm_i in names(lcomp_ob_p)){
    yrs_lcomp_ob_b_i <- rownames(lcomp_ob_b[[nm_i]])
    yrs_lcomp_ob_p_i <- rownames(lcomp_ob_p[[nm_i]])
    # If the last year of comps is within the last nyb_rcn$comp of the assessment, then project them.
    if(max(as.numeric(yrs_lcomp_ob_b_i))%in%tail(as.numeric(yb),nyb_rcn$comp)){
      lenprob_i <- lenprob[[nm_i]]
      if(grepl("^[LD]",nm_i)){ # If length comps are associated with catch (landings or discards), used numbers from landings
        n_i <- Cn_p[[paste0("Cn.",nm_i)]][yrs_lcomp_ob_p_i,,drop=FALSE]
      }else if(nm_i%in%names(NU_p)){
        n_i <- NU_p[[nm_i]][yrs_lcomp_ob_p_i,,drop=FALSE] # If not, see if it's associated with an index (e.g. fishery independent) and use numbers associated with the index
      }else if(nm_i%in%names(Nmisc_p)){
        n_i <- Nmisc_p[[nm_i]][yrs_lcomp_ob_p_i,,drop=FALSE] # If not, their should be an N-at-age matrix in Nmisc that matches
      }else{
        warning(paste0(nm_i, "from lcomp_ob_p doesn't appear to match any sources of F, indices, or other data sources associated with a selectivity."))
      }
      P_n_i <- t(apply(n_i,1,function(x){x/sum(x)}))
      lenbins_lenc_i <- dimnames(lcomp_ob_b[[nm_i]])[[2]]
      # Make sure the lenbins match the observed length comps and apply a plus group if necessary
      lcomp_ob_p_i <- local({
        a <- t(apply(P_n_i,1,function(x){colSums(x*lenprob_i)}))
        b <- a[,lenbins_lenc_i]
        bplusgroup <- tail(lenbins_lenc_i,1)
        b[,bplusgroup] <- b[,bplusgroup]+rowSums(a[,as.numeric(dimnames(a)[[2]])>as.numeric(bplusgroup)])
        b
      })
      lcomp_ob_p[[nm_i]] <- lcomp_ob_p_i
    }
  }

# Add projected values to base values
  ybp <- c(yb,yp)
  nybp <- length(ybp)

  N <- rbind(N,N_p)
  Nmdyr <- rbind(Nmdyr,Nmdyr_p)

  logRdev <- c(logRdev,logRdev_p)
  R <- c(R,R_p)

  B0 <- ifelse("B0"%in%names(parms),parms$B0,NA)
  SSB <- S <- c(S,S_p)
  SSB.SSBmsy <- if("SSBmsy"%in%names(parms)){
    SSB/parms$SSBmsy}else{
    SSB*NA
    }
  SSB.SSBF30 <- if("SSB.F30"%in%names(parms)){
    SSB/parms$SSB.F30}else{
      SSB*NA
    }
  SSB.msst <- if("msst"%in%names(parms)){
    SSB/parms$msst}else{
      SSB*NA
    }
  SSB.msst.F30 <- if("msst.F30"%in%names(parms)){
    SSB/parms$msst.F30}else{
      SSB*NA
    }

  # U <- rbind(U,U_p) # U is updated above
  # Round indices to number of digits used in observed data
  U <- local({
    a <- U
    for(i in 1:ncol(a)){
      nmi <- names(a)[i]
      ai <- a[,i]
      ndigi <- round(median(nchar(gsub("([0-9]*.)([0-9]*)","\\2",na.omit(U_ob[,nmi])))))
      a[,i] <- round(ai,ndigi)
    }
    a
  })
  U_p <- U[paste(yp),] # get the rounded version of U_p

  U_a <- local({
    a <- lapply(seq_along(U_a),function(i){rbind(U_a[[i]],U_a_p[[i]])})
    names(a) <- names(U_a)
    a
  })
  cv_U <- rbind(cv_U,cv_U_p)
  dimnames(q_p)[[2]] <- paste0("q.",dimnames(q_p)[[2]])
  q <- rbind(q,q_p)

  cv_C <- rbind(cv_C,cv_C_p)
  nfish <- rbind(nfish,nfish_p)
  ntrip <- rbind(ntrip,ntrip_p)
  neff <- rbind(neff,neff_p)

  # Fishgraph gets mad if you have NA values in nfish, ntrip, or neff
  nfish[is.na(nfish)] <- -99999
  ntrip[is.na(ntrip)] <- -99999
  neff[is.na(neff)] <- -99999


  # Add projected acomp_ob
  acomp_ob <- local({
    a <- lapply(seq_along(names(acomp_ob_b)),function(i){
      nm_i <- names(acomp_ob_b)[i]
      ai <- rbind(acomp_ob_b[[i]],acomp_ob_p[[i]])
      ai[complete.cases(ai),,drop=FALSE]
    })
    names(a) <- names(acomp_ob_b)
    a
  })

  # Add NA cells to match structure of projected acomp_ob
  acomp_pr <- local({
    a <- lapply(seq_along(names(acomp_pr_b)),function(i){
      nm_i <- names(acomp_pr_b)[i]
      ai <- acomp_ob[[nm_i]]*NA # Initialize with structure of acomp_ob
      ai[rownames(acomp_pr_b[[i]]),] <- acomp_pr_b[[i]] # Fill with base values
      # This will satisfy FishGraph, but it seems inappropriate to fill in bam predicted comps in the projection
      # ai <- rbind(acomp_pr_b[[i]],acomp_ob_p[[i]])
      # ai[complete.cases(ai),,drop=FALSE]
      ai
    })
    names(a) <- names(acomp_pr_b)
    a
  })

  lcomp_ob <- local({
    a <- lapply(seq_along(names(lcomp_ob_b)),function(i){
      nm_i <- names(lcomp_ob_b)[i]
      ai <- rbind(lcomp_ob_b[[i]],lcomp_ob_p[[i]])
      ai[complete.cases(ai),,drop=FALSE]
    })
    names(a) <- names(lcomp_ob_b)
    a
  })

  # Add NA cells to match structure of projected lcomp_ob
  lcomp_pr <- local({
    a <- lapply(seq_along(names(lcomp_pr_b)),function(i){
      nm_i <- names(lcomp_pr_b)[i]
      ai <- lcomp_ob[[nm_i]]*NA # Initialize with structure of lcomp_ob
      ai[rownames(lcomp_pr_b[[i]]),] <- lcomp_pr_b[[i]] # Fill with base values
      # This will satisfy FishGraph, but it seems inappropriate to fill in bam predicted comps in the projection
      # ai <- rbind(lcomp_pr_b[[i]],lcomp_ob_p[[i]])
      # ai[complete.cases(ai),,drop=FALSE]
      ai
    })
    names(a) <- names(lcomp_pr_b)
    a
  })

  Z <- rbind(Z,Z_p)
  Fsum <-  c(Fsum, Fsum_p)
  Ffull <- c(Ffull, Ffull_p)
  Fsum_flt <- rbind(Fsum_flt,Fsum_flt_p)
  F.Fmsy <- Ffull/parms$Fmsy

  for(nm_i in names(Cn)){
    Cn[[nm_i]] <- rbind(Cn[[nm_i]],Cn_p[[nm_i]])
  }
  for(nm_i in names(Cw)){
    Cw[[nm_i]] <- rbind(Cw[[nm_i]],Cw_p[[nm_i]])
  }
    #######################################
  } ########## end if nyp>0 ###############
    #######################################

    B <- wgt_mt*N
    Nsum <- rowSums(N)
    Bsum <- rowSums(B)


  # Values by fleet  in long format (for ggplot)
  # Landings (n)
  Cn.L <- local({
    a <- as.data.frame(lapply(Cn[grepl("^Cn.L.",names(Cn))],rowSums))
    b <- cbind("year"=as.numeric(rownames(a)),a)
    colnames(b) <- gsub("Cn.L.","",colnames(b))
    b
  })
  Cn.L.long <- gather(data=Cn.L,key=fleet,value=Ln,names(Cn.L)[names(Cn.L)!="year"],factor_key = TRUE)

  # Discards (n)
  Cn.D <- local({
    a <- as.data.frame(lapply(Cn[grepl("^Cn.D.",names(Cn))],rowSums))
    b <- cbind("year"=as.numeric(rownames(a)),a)
    colnames(b) <- gsub("Cn.D.","",colnames(b))
    b
  })
  Cn.D.long <- gather(data=Cn.D,key=fleet,value=Dn,names(Cn.D)[names(Cn.D)!="year"],factor_key = TRUE)

  # Landings (w)
  Cw.L <- local({
    a <- as.data.frame(lapply(Cw[grepl("^Cw.L.",names(Cw))],rowSums))
    b <- cbind("year"=as.numeric(rownames(a)),a)
    colnames(b) <- gsub("Cw.L.","",colnames(b))
    b
  })
  Cw.L.long <- gather(data=Cw.L,key=fleet,value=Lw,names(Cw.L)[names(Cw.L)!="year"],factor_key = TRUE)

  # Discards (w)
  Cw.D <- local({
    a <- as.data.frame(lapply(Cw[grepl("^Cw.D.",names(Cw))],rowSums))
    b <- cbind("year"=as.numeric(rownames(a)),a)
    colnames(b) <- gsub("Cw.D.","",colnames(b))
    b
    })

  Cw.D.long <- gather(data=Cw.D,key=fleet,value=Dw,names(Cw.D)[names(Cw.D)!="year"],factor_key = TRUE)

  Cn.L.tot <- tapply(Cn.L.long[,"Ln"],Cn.L.long[,"year"],sum)
  Cw.L.tot <- tapply(Cw.L.long[,"Lw"],Cw.L.long[,"year"],sum)

  if(nrow(Cn.D.long)>0){
  Cn.D.tot <- tapply(Cn.D.long[,"Dn"],Cn.D.long[,"year"],sum)
  Cw.D.tot <- tapply(Cw.D.long[,"Dw"],Cw.D.long[,"year"],sum)
  }else{
    Cn.D.tot <- setNames(rep(NA,nybp),paste(ybp))
    Cw.D.tot <- setNames(rep(NA,nybp),paste(ybp))
  }


################################################################################
### Build projection output ####################################################
################################################################################
out_proj <- list(
  results = list(
    t_series=cbind(data.frame(years=ybp,
                              # Fsum=Fsum,
                              Ffull=Ffull,
                              L_wgt=Cw.L.tot,
                              L_num=Cn.L.tot,
                              D_wgt=Cw.D.tot,
                              D_num=Cn.D.tot,
                              N=Nsum,
                              S=S,
                              B=Bsum,
                              R=R,
                              Rlogdev=logRdev
    ),
    U),
    Nage = N,
    Uage = U_a,
    Cn=Cn,
    Cw=Cw,
    Cn.L.long=Cn.L.long,
    Cw.L.long=Cw.L.long,
    Cn.D.long=Cn.D.long,
    Cw.D.long=Cw.D.long
  ),
  results_p = list(
    t_series=cbind(data.frame(years_p=yp,
                              # Fsum=Fsum_p,
                              Ffull=Ffull_p,
                              L_wgt=Lw_p,
                              L_num=Ln_p,
                              D_wgt=Dw_p,
                              D_num=Dn_p,
                              S=S_p,
                              B=B_p,
                              R=R_p,
                              Rlogdev=logRdev_p
    ),
    U_p),
    Nage = N_p,
    Uage = U_a_p
  ),
  parms_bam = parms#,
  # Fage=Fage,
  # Fage_p=Fage_p,
  # Fage=Fage,
  # Fage_p=Ffull_p,
  # F_flt=F_flt,
  # F_flt_p=F_flt_p,
  # Fage_p=Fage_p,
  # bam_p = bam_p
)

#### Add projected objects to rdat_proj
  rdat_proj <- rdat

  endyr_proj <- max(yp)

  # a.series (no need to change at this time 2025-08-26 NPK)

  # t.series
  t.series.proj <- local({
    a <- t.series[paste(styr:endyr),]
    b <- matrix(NA,nrow=nyp,ncol=ncol(t.series),dimnames=list(yp,colnames(t.series)))
    out <- rbind(a,b)
    out$year <- ybp
    out$recruits <- R
    out$logR.dev <- logRdev
    out$SSB <- S
    out$SSB.SSBmsy <- SSB.SSBmsy
    out$SSB.SSBF30 <- SSB.SSBF30
    out$SSB.msst <- SSB.msst
    out$SSB.msstF30 <- SSB.msst
    out$F.full <- Ffull
    out$Fsum <- Fsum
    out$F.Fmsy <- F.Fmsy
    out$N <- Nsum
    out$B <- Bsum
    out$B.B0 <- out$B/B0
    out
  })

  # Add F values for each flt
  t.series.proj[,key_abb_F[names(Fsum_flt)]] <- Fsum_flt

  ## Add landings and (live) discards for each fleet
  # Landings
  abb_L <- abb_LD[grepl("^L",abb_LD)] # Abbreviations of landings columns in t.series
  # unit_L <- setNames(rep(NA,length(abb_L)),abb_L)
  Cn.L.b <- Cn.L[paste(yb),names(Cn.L)!="year"]
  Cw.L.b <- Cw.L[paste(yb),names(Cw.L)!="year"]
  Cn.L.p <- Cn.L[paste(yp),names(Cn.L)!="year"]
  Cw.L.p <- Cw.L[paste(yp),names(Cw.L)!="year"]

  for(i in seq_along(abb_L)){
    abb_i <- gsub("^L.","",abb_L[i])
    nm.L.ob.ts.i <- paste0("L.",abb_i,".ob") # Name of the observed landings column in t.series
    nm.L.pr.ts.i <- paste0("L.",abb_i,".pr") # Name of the predicted landings column in t.series
    # base years
    Cn.L.b.i <- setNames(Cn.L.b[,abb_i],yb) # predicted landings in numbers
    Cw.L.b.i <- setNames(Cw.L.b[,abb_i],yb) # predicted landings in weight (often in 1000 lb. Should probably be in standard units in all assessments.)
    # projection years
    Cn.L.p.i <- setNames(Cn.L.p[,abb_i],yp) # predicted landings in numbers
    Cw.L.p.i <- setNames(Cw.L.p[,abb_i],yp) # predicted landings in weight
    L.C.p.i <- data.frame("n"=Cn.L.p.i,"w"=Cw.L.p.i)
    L.ob.ts.i <- setNames(t.series[paste(yb),nm.L.ob.ts.i],yb) # observed landings in unknown unit
    L.pr.ts.i <- setNames(t.series[paste(yb),nm.L.pr.ts.i],yb) # predicted landings in unknown unit
    lmfitn.i <- lm(L.pr.ts.i~Cn.L.b.i)
    lmfitw.i <- lm(L.pr.ts.i~Cw.L.b.i)
    r.squared.i <- c("n"=summary(lmfitn.i)$r.squared,"w"=summary(lmfitw.i)$r.squared)
    r.squared.max.i <- max(r.squared.i)
    slope.i <- c("n"=coef(lmfitn.i)[[2]],"w"=coef(lmfitw.i)[[2]])
    # which catch series (numbers or weight) best predicts what's in t.series?
    unit.i <- names(which(r.squared.i==r.squared.max.i))
    # What is the slope of that relationship? (i.e. is a unit conversion necessary)
    slope.unit.i <- signif(slope.i[[unit.i]],6)
    # Get projected landings series in appropriate units and scale by slope of linear model
    L.ob.ts.p.i <- L.C.p.i[,unit.i]*slope.unit.i # Compute projected landings
    if(is.na(L.ob.ts.i[[paste(endyr)]])){L.ob.ts.p.i <- L.ob.ts.p.i*NA} # If observed landings were NA in endyr, set projected values to NA
    ndigi <- round(median(nchar(gsub("([0-9]*.)([0-9]*)","\\2",paste(na.omit(L.ob.ts.i)))))) # Number of digits in observed data
    L.ob.ts.p.i <- round(L.ob.ts.p.i,ndigi)
    t.series.proj[paste(yp),paste0("L.",abb_i,".ob")] <- L.ob.ts.p.i # Add to observed column only
    if(r.squared.max.i<0.95){ # Test to see if the best r.squared is high enough. It should be very close to 1.
      warning(paste(nm.L.pr.ts.i, "is not matching up well with data from Cn or Cw. Check that units are correct in the projection"))
    }else{
      message(paste0(abb_L[i], " appears to be in ",unit.i," units in t.series, scaled by ",slope.unit.i," compared to the analog in CLD.est.mats"))
    }
  }
  nm.L.ob.ts <- paste0(abb_L,".ob") # Names of observed landings, as they appear in t.series

  # Discards (note discards in t.series and CLD.est.mats are all dead, so there is no need to reconsider discard mortality)
  abb_D <- abb_LD[grepl("^D",abb_LD)] # Abbreviations of landings columns in t.series
  # unit_D <- setNames(rep(NA,length(abb_D)),abb_D)
  Cn.D.b <- Cn.D[paste(yb),names(Cn.D)!="year"]
  Cw.D.b <- Cw.D[paste(yb),names(Cw.D)!="year"]
  Cn.D.p <- Cn.D[paste(yp),names(Cn.D)!="year"]
  Cw.D.p <- Cw.D[paste(yp),names(Cw.D)!="year"]

  for(i in seq_along(abb_D)){
    abb_i <- gsub("^D.","",abb_D[i])
    nm.D.ob.ts.i <- paste0("D.",abb_i,".ob") # Name of the observed discards column in t.series
    nm.D.pr.ts.i <- paste0("D.",abb_i,".pr") # Name of the predicted discards column in t.series
    # base years
    Cn.D.b.i <- setNames(Cn.D.b[,abb_i],yb) # predicted discards in numbers
    Cw.D.b.i <- setNames(Cw.D.b[,abb_i],yb) # predicted discards in weight (often in 1000 lb. Should probably be in standard units in all assessments.)
    # projection years
    Cn.D.p.i <- setNames(Cn.D.p[,abb_i],yp) # predicted discards in numbers
    Cw.D.p.i <- setNames(Cw.D.p[,abb_i],yp) # predicted discards in weight
    D.C.p.i <- data.frame("n"=Cn.D.p.i,"w"=Cw.D.p.i)
    D.ob.ts.i <- setNames(t.series[paste(yb),nm.D.ob.ts.i],yb) # observed discards in unknown unit
    D.pr.ts.i <- setNames(t.series[paste(yb),nm.D.pr.ts.i],yb) # predicted discards in unknown unit
    lmfitn.i <- lm(D.pr.ts.i~Cn.D.b.i)
    lmfitw.i <- lm(D.pr.ts.i~Cw.D.b.i)
    r.squared.i <- c("n"=summary(lmfitn.i)$r.squared,"w"=summary(lmfitw.i)$r.squared)
    r.squared.max.i <- max(r.squared.i)
    slope.i <- c("n"=coef(lmfitn.i)[[2]],"w"=coef(lmfitw.i)[[2]])
    # which catch series (numbers or weight) best predicts what's in t.series?
    unit.i <- names(which(r.squared.i==r.squared.max.i))
    # What is the slope of that relationship? (i.e. is a unit conversion necessary)
    slope.unit.i <- signif(slope.i[[unit.i]],6)
    # Get projected discards series in appropriate units and scale by slope of linear model
    D.ob.ts.p.i <- D.C.p.i[,unit.i]*slope.unit.i # Compute projected discards
    if(is.na(D.ob.ts.i[[paste(endyr)]])){D.ob.ts.p.i <- D.ob.ts.p.i*NA} # If observed discards were NA in endyr, set projected values to NA
    ndigi <- round(median(nchar(gsub("([0-9]*.)([0-9]*)","\\2",paste(na.omit(D.ob.ts.i))))))
    D.ob.ts.p.i <- round(D.ob.ts.p.i,ndigi)
    t.series.proj[paste(yp),paste0("D.",abb_i,".ob")] <- D.ob.ts.p.i # Add to observed column only
    if(r.squared.max.i<0.95){ # Test to see if the best r.squared is high enough. It should be very close to 1.
      warning(paste(nm.D.pr.ts.i, "is not matching up well with data from Cn or Cw. Check that units are correct in the projection"))
    }else{
      message(paste0(abb_D[i], " appears to be in ",unit.i," units in t.series, scaled by ",slope.unit.i," compared to the analog in CLD.est.mats"))
    }
  }

  nm.D.ob.ts <- paste0(abb_D,".ob") # Names of observed (dead) discards, as they appear in t.series

  # cv of L and D
  t.series.proj[paste(yp),colnames(cv_C_p)] <- cv_C_p # Add projected cvs

  # Indices of abundance
  t.series.proj[paste(yp),paste0("U.",colnames(U_p),".ob")] <- U_p # Add projected indices

  # cv of indices of abundance
  t.series.proj[paste(yp),colnames(cv_U_p)] <- cv_U_p # Add projected cvs

  # catchability (q) of indices of abundance
  t.series.proj[paste(yp),dimnames(q_p)[[2]]] <- q_mn_p

  # add q.rate.mult values to t.series.proj, when appropriate
  nm_q_ratemult <- paste("q",dimnames(q_ratemult_p)[[2]],"rate.mult",sep=".")
  if(any(nm_q_ratemult%in%names(t.series.proj))){
    nm_q_ratemult_sub <- nm_q_ratemult[nm_q_ratemult%in%names(t.series.proj)]
    fleet_abb_sub <- gsub("(^q.)([a-zA-Z0-9]+)(.rate.mult)","\\2",nm_q_ratemult_sub)
    t.series.proj[paste(yp),nm_q_ratemult_sub] <- q_ratemult_p[paste(yp),fleet_abb_sub]
  }

  # density dependent function as a multiple of q (scaled a la Katsukawa and Matsuda. 2003)
  t.series.proj[paste(yp),"q.DD.mult"] <- q_DD_mult_p

  # Comp nfish and ntrip
  t.series.proj[,names(key_nm_ntrip)] <- ntrip
  t.series.proj[,names(key_nm_neff)] <- neff
  t.series.proj[,names(key_nm_nfish)] <- nfish

  ## Modify rdat_proj
  # parms
  rdat_proj$parms$endyr <- endyr_proj
  rdat_proj$t.series <- t.series.proj

  # CLD.est.mats
  CLD.est.mats[names(key_nm_Cn)] <- Cn
  CLD.est.mats[names(key_nm_Cw)] <- Cw
  CLD.est.mats$Ln.total <- Reduce("+",Cn[grepl("^Cn.L.",names(Cn))])
  CLD.est.mats$Lw.total <- Reduce("+",Cw[grepl("^Cw.L.",names(Cw))])
  CLD.est.mats$Dn.total <- Reduce("+",Cn[grepl("^Cn.D.",names(Cn))])
  CLD.est.mats$Dw.total <- Reduce("+",Cw[grepl("^Cw.D.",names(Cw))])

  rdat_proj$CLD.est.mats <- CLD.est.mats

  # Z.age, N.age, N.mdyr
  rdat_proj$Z.age <- Z
  rdat_proj$N.age <- N
  rdat_proj$N.age.mdyr <- Nmdyr
  rdat_proj$B.age <- B

  # size.age.fishery
  # Fill projection years with values from last year of base (endyr)
  # Note: This projection is done the same way above in code that actually uses
  # the weights and lengths. If any more complicated projection of length or weight
  # at age is desired, it should be coded in both places. At present (2025-09-08)
  # it is far easier to duplicate this simple repeating of endyr values in two
  # places than it is to match the objects above with size.age.fishery in the rdat file.
  nm.size.age.fishery <- names(size.age.fishery)
  size.age.fishery <- lapply(1:length(size.age.fishery),function(j){
    nmj <- names(size.age.fishery)[j]
    xj <- size.age.fishery[[j]]
    xj_p <- matrix(xj[paste(endyr),],nrow=nyp,ncol=ncol(xj),byrow = TRUE,dimnames = list(yp,colnames(xj)))
    rbind(xj,xj_p)
  })
  names(size.age.fishery) <- nm.size.age.fishery

  rdat_proj$size.age.fishery <- size.age.fishery

  # sel.age
  # Fill projection years with values from last year of base (endyr)
  # Note: This projection is done the same way above in code that actually uses
  # the selectivities. If any more complicated projection of length or weight
  # at age is desired, it should be coded in both places. At present (2025-09-09)
  # it is far easier to duplicate this simple repeating of endyr values in two
  # places than it is to match the objects above with sel.age in the rdat file.
  nm.sel.age <- names(sel.age)
  sel.age <- lapply(1:length(sel.age),function(j){
    nmj <- names(sel.age)[j]
    xj <- sel.age[[j]]
    if(is.matrix(xj)){
    xj_p <- matrix(xj[paste(endyr),],nrow=nyp,ncol=ncol(xj),byrow = TRUE,dimnames = list(yp,colnames(xj)))
    rbind(xj,xj_p)
    }else{
      xj
    }
  })
  names(sel.age) <- nm.sel.age

  rdat_proj$sel.age <- sel.age

  # comp.mats
  comp.mats[names(key_nm_acomp_ob)] <- acomp_ob[key_nm_acomp_ob]
  comp.mats[names(key_nm_acomp_pr)] <- acomp_pr[key_nm_acomp_pr]
  comp.mats[names(key_nm_lcomp_ob)] <- lcomp_ob[key_nm_lcomp_ob]
  comp.mats[names(key_nm_lcomp_pr)] <- lcomp_pr[key_nm_lcomp_pr]

  rdat_proj$comp.mats <- comp.mats

# Add to out_proj
  out_proj$rdat_proj <- rdat_proj

#############################################################
###### Extend data inputs and build projected dat file ######
#############################################################
  # Identify the parts of init that need to be changed
  # Identify the new parts and match them with the parts that need to change
  if(project_bam&!is.null(args_bam2r)){
    args_bam2r <- modifyList(x=as.list(formals(bam2r)),val=args_bam2r) # Add user supplied arguments to defaults
    bam <- do.call(bam2r,args_bam2r)
    init_p <- init_b <- bam$init

    ### Identify temporal objects in init that need to change
    ## Check to see if life history inputs are time varying. If so extend them accordingly
    # maturity
    obs_maturity_nm <- names(init_p)[grepl("^obs_maturity",names(init_p))]
    if(length(obs_maturity_nm)>0){
      for(nm_i in obs_maturity_nm){
        xi <- init_p[[nm_i]]
        if(is.matrix(xi)){
          message(paste0(nm_i," is a matrix"))
          # Check if rownames are equal to model years (they should be if it is time varying)
          if(all(dimnames(xi)[[1]]==paste(yb))){
            message(paste0(nm_i," appears to be time varying. endyr values will be projected"))
            xi_p <- matrix(xi[paste(endyr),],nrow=nyp,ncol=ncol(xi),dimnames=list(paste(yp),dimnames(xi)[[2]]),byrow=TRUE)
            init_p[[nm_i]] <- rbind(xi,xi_p)

            }else{
            warning(paste0("rownames of ",nm_i," are not equal to model years"))
          }
        }

      }
    }
    # proportion female or male
    obs_prop_nm <- names(init_p)[grepl("^obs_prop_[fmFM]",names(init_p))]
    if(length(obs_prop_nm)>0){
      for(nm_i in obs_prop_nm){
        xi <- init_p[[nm_i]]
        if(is.matrix(xi)){
          message(paste0(nm_i," is a matrix"))
          # Check if rownames are equal to model years (they should be if it is time varying)
          if(all(dimnames(xi)[[1]]==paste(yb))){
            message(paste0(nm_i," appears to be time varying. endyr values will be projected"))
            xi_p <- matrix(xi[paste(endyr),],nrow=nyp,ncol=ncol(xi),dimnames=list(paste(yp),dimnames(xi)[[2]]),byrow=TRUE)
            init_p[[nm_i]] <- rbind(xi,xi_p)

          }else{
            warning(paste0("rownames of ",nm_i," are not equal to model years"))
          }
        }
      }
    }

    # natural mortality-at-age
    set_M_nm <- names(init_p)[grepl("^set_M",names(init_p))]
    if(length(set_M_nm)>0){
      for(nm_i in set_M_nm){
        xi <- init_p[[nm_i]]
        if(is.matrix(xi)){
          message(paste0(nm_i," is a matrix"))
          # Check if rownames are equal to model years (they should be if it is time varying)
          if(all(dimnames(xi)[[1]]==paste(yb))){
            message(paste0(nm_i," appears to be time varying. endyr values will be projected"))
            xi_p <- matrix(xi[paste(endyr),],nrow=nyp,ncol=ncol(xi),dimnames=list(paste(yp),dimnames(xi)[[2]]),byrow=TRUE)
            init_p[[nm_i]] <- rbind(xi,xi_p)

          }else{
            warning(paste0("rownames of ",nm_i," are not equal to model years"))
          }
        }
      }
    }

    # I'll have to do some more monkeying around to deal with the tv objects in menhaden
    # mostly because they have weird names
    # tv (time varying objects currently in Atlantic Menhaden model)
    # init_tv <- init_b[grepl("_tv$",names(init_b))]

    #### Update temporal values
    ## _endyr
    # Any _endyr for select values and for data sets that extend to the terminal
    # year of the assessment should be extended by nyp
    yr_nm_add_p <- local({
      a <- nm_yr_p
      a[a%in%names(init_b)]
    })
    init_p[yr_nm_add_p] <- lapply(init_p[yr_nm_add_p],function(x){paste(as.numeric(x)+nyp)})

    # Allows fairly flexible naming of rec dev parameters
    nm_styr_dev_rec <- names(init_p)[grepl("^(?=.*styr)(?=.*rec)(?=.*dev).*$",names(init_p),perl=TRUE)]
    nm_endyr_dev_rec <- names(init_p)[grepl("^(?=.*endyr)(?=.*rec)(?=.*dev).*$",names(init_p),perl=TRUE)]
    nm_set_log_dev_vals_rec <- names(init_p)[grepl("^(?=.*set)(?=.*log)(?=.*dev)(?=.*vals)(?=.*rec)(?=.*dev).*$",names(init_p),perl=TRUE)]

    yrs_rec_dev <- init_p[[nm_styr_dev_rec]]:init_p[[nm_endyr_dev_rec]]
    init_p[[nm_set_log_dev_vals_rec]] <- setNames(rep("0.0",length(yrs_rec_dev)),yrs_rec_dev)

    ## landings
    # Note: need to make sure you get the right units (n or w)
    init_obs_L_nm <- names(init_b)[grepl("^obs_L",names(init_b))]
    for(i in seq_along(init_obs_L_nm)){
      nm_i <- init_obs_L_nm[i]
      abb_i <- gsub("^(obs_L)(_)(.*)","\\3",nm_i)
      nm_cv_i <- paste0("obs_cv_L_",abb_i)
      xbi <- setNames(as.numeric(init_b[[nm_i]]),names(init_b[[nm_i]]))
      xcvbi <- setNames(as.numeric(init_b[[nm_cv_i]]),names(init_b[[nm_cv_i]]))
      yrsbi <- names(xbi)

      ndigi <- round(median(nchar(gsub("([0-9]*.)([0-9]*)","\\2",xbi))))

      # Landings for projection years, for fleet i
      # If the last year of xbi was the last year of the base model, then project it forward
      if(paste(endyr)%in%names(xbi)){
        xpi <- round(setNames(rdat_proj$t.series[paste(yp),paste0("L.",abb_i,".ob")],yp),ndigi)
        xi <- c(xbi,xpi)
      }else{
        xi <- xbi
      }

      init_p[[nm_i]] <- setNames(paste(xi),names(xi))
      # reset appropriate _endyr value to make sure it agrees with the projected data (some endyr values might get projected even if the data doesn't)
      endyr_nm_i <- paste0("endyr_L_",abb_i)
      if(endyr_nm_i%in%names(init_p)){init_p[[endyr_nm_i]] <- tail(names(xi),1)}

      # If the last year of xcvbi was the last year of the base model, then project it forward
      if(paste(endyr)%in%names(xcvbi)){
        xcvpi <- round(setNames(rdat_proj$t.series[paste(yp),paste0("cv.L.",abb_i)],yp),ndigi)
        xcvi <- c(xcvbi,xcvpi)
      }else{
        xcvi <- xcvbi
      }

      init_p[[nm_cv_i]] <- setNames(paste(xcvi),names(xcvi))

      ## set_log_dev_vals_F_L
      set_log_dev_vals_F_L_nm_i <- paste0("set_log_dev_vals_F_L_",abb_i)
      if(set_log_dev_vals_F_L_nm_i%in%names(init_b)){
        init_p[[set_log_dev_vals_F_L_nm_i]] <- setNames(rep("0.0",length(xi)),names(xi))
      }
    }

    ## discards
    # Note: Assumed to be in numbers
    init_obs_released_nm <- names(init_b)[grepl("^obs_released",names(init_b))]
    for(i in seq_along(init_obs_released_nm)){
      nm_i <- init_obs_released_nm[i]
      abb_i <- gsub("^(obs_released)(_)(.*)","\\3",nm_i)
      #nm_cv_i <- paste0("obs_cv_D_",abb_i)

      xbi <- setNames(as.numeric(init_b[[nm_i]]),names(init_b[[nm_i]])) # Live discards
      #xcvbi <- setNames(as.numeric(init_b[[nm_cv_i]]),names(init_b[[nm_cv_i]]))
      yrsbi <- names(xbi)

      ndigi <- round(median(nchar(gsub("([0-9]*.)([0-9]*)","\\2",xbi))))

      # Dead discards for fleet i
      xbi_dead <- setNames(rdat_proj$t.series[paste(yb),paste0("D.",abb_i,".ob")],yb) # base years
      xpi_dead <- setNames(rdat_proj$t.series[paste(yp),paste0("D.",abb_i,".ob")],yp) # projection years

      # Discard mortality rate
      nm_Dmort_i <- paste0("set_Dmort_",abb_i)
      # xni_dead <- rowSums(Cn[[paste0("Cn.D.",abb_i)]])
      if(nm_Dmort_i%in%names(init_b)){
        Dmort_i <- as.numeric(init_b[[paste0("set_Dmort_",abb_i)]])
      }else{
        Dmort_i <- mean(xbi_dead[yrsbi]/xbi,na.rm=TRUE) # Compute discard mortality by dividing dead discards by live discards for base years
      }

      # If the last year of xbi was the last year of the base model, then project it forward
      if(paste(endyr)%in%names(xbi)){
        xpi <- round(xpi_dead*(1/Dmort_i),ndigi) # Convert to released (live discards)
        xi <- c(xbi,xpi)
      }else{
        xi <- xbi
      }

      init_p[[nm_i]] <- setNames(paste(xi),names(xi))
      # reset appropriate _endyr value to make sure it agrees with the projected data (some endyr values might get projected even if the data doesn't)
      endyr_nm_i <- paste0("endyr_D_",abb_i)
      if(endyr_nm_i%in%names(init_p)){init_p[[endyr_nm_i]] <- tail(names(xi),1)}

      # xcvpi <- setNames(rdat_proj$t.series[paste(yp),paste0("cv.D.",abb_i)],yp)
      # xcvi <- c(xcvbi,xcvpi)
      #
      # init_p[[nm_cv_i]] <- setNames(paste(xcvi),names(xcvi))
    }

    ## discard cvs
    #  This is done in a separate loop because the discard cvs don't always match discard time series supplied to the dat file
    #  (e.g. Black Sea Bass SEDAR 56)
    init_obs_cv_D_nm <- names(init_b)[grepl("^obs_cv_D",names(init_b))]
    for(i in seq_along(init_obs_cv_D_nm)){
      nm_cv_i <- init_obs_cv_D_nm[i]
      abb_i <- gsub("^(obs_cv_D)(_)(.*)","\\3",nm_cv_i)
      xcvbi <- setNames(as.numeric(init_b[[nm_cv_i]]),names(init_b[[nm_cv_i]]))
      # If the last year of xcvbi was the last year of the base model, then project it forward
      if(paste(endyr)%in%names(xcvbi)){
        xcvpi <- setNames(rdat_proj$t.series[paste(yp),paste0("cv.D.",abb_i)],yp)
        xcvi <- c(xcvbi,xcvpi)
      }else{
        xcvi <- xcvbi
      }
      init_p[[nm_cv_i]] <- xcvi
    }

    ## set_log_dev_vals_F_D
    init_set_log_dev_vals_F_D_nm <- names(init_b)[grepl("^set_log_dev_vals_F_D",names(init_b))]
    for(nm_i in init_set_log_dev_vals_F_D_nm){
      vals_b_i <- init_b[[nm_i]]
      # If the value was provided in the last year of the base model, project it forward
      if(paste(endyr)%in%names(vals_b_i)){
        vals_p_i <- setNames(rep("0.0",nyp),yp)
        vals_i <- c(vals_b_i,vals_p_i)
      }else{
        vals_i <- vals_b_i
      }
      init_p[[nm_i]] <- vals_i
    }

    ## cpue
    init_obs_cpue_nm <- local({
      nm <- names(init_b)[grepl("^obs_cpue",names(init_b))]
      nm_abb <- gsub("^(obs_cpue)(_)(.*)","\\3",nm)
      nm_U_ob_ts <- paste0("U.",nm_abb,".ob")
      nm_yes <- nm[which(nm_U_ob_ts%in%names(rdat_proj$t.series))]
      nm_no <- nm[which(!nm_U_ob_ts%in%names(rdat_proj$t.series))]
      nm_no_ts <- nm_U_ob_ts[match(nm_no,nm)]
      if(length(nm_no)>0){
        message(paste0(paste(nm_no,collapse=", "), " found in the bam tpl but ",paste(nm_no_ts,collapse=", ")," not found in the rdat t.series. May not be included in the likelihood."))
      }
      # Only include names found in U (indices reported in the rdat)
      return(nm_yes)
    })

    for(i in seq_along(init_obs_cpue_nm)){
      nm_i <- init_obs_cpue_nm[i]
      abb_i <- gsub("^(obs_cpue)(_)(.*)","\\3",nm_i)
      nm_cv_i <- paste0("obs_cv_cpue_",abb_i)
      xbi <- setNames(as.numeric(init_b[[nm_i]]),names(init_b[[nm_i]])) # get index values from init

      ndigi <- round(median(nchar(gsub("([0-9]*.)([0-9]*)","\\2",xbi))))

      xpi <- round(setNames(rdat_proj$t.series[paste(yp),paste0("U.",abb_i,".ob")],yp),ndigi)
      xi <- c(xbi,xpi)


      ndigi <- round(median(nchar(gsub("([0-9]*.)([0-9]*)","\\2",xbi))))
      xpi <- round(setNames(rdat_proj$t.series[paste(yp),paste0("U.",abb_i,".ob")],yp),ndigi)
      xi <- c(xbi,xpi[!is.na(xpi)])
      yrsi <- names(xi)

      xcvbi <- init_b[[nm_cv_i]] # get index cv values from init
      ndigcvi <- round(median(nchar(gsub("([0-9]*.)([0-9]*)","\\2",xcvbi))))
      xcvpi <- round(setNames(rdat_proj$t.series[paste(yp),paste0("cv.U.",abb_i)],yp),ndigcvi)
      xcvi <- c(xcvbi,xcvpi[!is.na(xcvpi)])

      init_p[[nm_i]] <- setNames(paste(xi),names(xi))
      # init_p[[nm_cv_i]] <- setNames(paste(round(cv_U[,paste0("cv.U.",abb_i)],ndigcvi)),rownames(cv_U))[names(xi)]
      init_p[[nm_cv_i]] <- setNames(paste(xcvi),names(xcvi))

      # update styr, endyr, yrs, and nyr as necessary
      yrinfoi <- setNames(list(min(yrsi),max(yrsi),yrsi,paste(length(yrsi))),
                          paste0(c("styr","endyr","yrs","nyr"),gsub("obs","",nm_i))
      )
      yrinfoi_nm_is <- names(yrinfoi)[names(yrinfoi)%in%names(init_p)]
      init_p[yrinfoi_nm_is] <- yrinfoi[yrinfoi_nm_is]
    }

    ## agec
    init_obs_agec_nm <- names(init_b)[grepl("^obs_agec",names(init_b))]
    key_nm_init_obs_agec <- setNames(names(acomp_ob),
                                     paste0("obs_agec_",gsub("^(D.)([A-Za-z]+)","\\2_D",gsub("^L.","",names(acomp_ob)))))
    # If the age comp names aren't using the _D suffix, attempt to remove _D suffix from the names in the key
    if(!any(grepl("\\_D$",init_obs_agec_nm))){
      names(key_nm_init_obs_agec) <- gsub("\\_D$","",names(key_nm_init_obs_agec))
    }

    for(i in seq_along(init_obs_agec_nm)){
      nm_i <- init_obs_agec_nm[i] # init object name
      nm2_i <- key_nm_init_obs_agec[[nm_i]] # acomp_ob name
      nm3_i <- names(key_nm_acomp_ob)[match(nm2_i,key_nm_acomp_ob)] # comp.mats name
      xbi <- init_b[[nm_i]]
      #x2i <- acomp_ob[[nm2_i]]
      x3i <- rdat_proj$comp.mats[[nm3_i]]

      ndigi <- round(median(nchar(gsub("([0-9]*.)([0-9]*)","\\2",xbi))))
      obsi <- apply(x3i,2,function(x){sprintf(paste0("%.",ndigi,"f"), x)})
      attributes(obsi) <- attributes(x3i)
      yrsi <- rownames(obsi)
      init_p[[nm_i]] <- obsi

      # update styr, endyr, yrs, and nyr as necessary
      yrinfoi <- setNames(list(min(yrsi),max(yrsi),yrsi,paste(length(yrsi))),
                          paste0(c("styr","endyr","yrs","nyr"),gsub("obs","",nm_i))
      )
      yrinfoi_nm_is <- names(yrinfoi)[names(yrinfoi)%in%names(init_p)]
      init_p[yrinfoi_nm_is] <- yrinfoi[yrinfoi_nm_is]

      # update nfish and nsamp as necessary
      nm_nfish_i <- gsub("^obs","nfish",nm_i) # name in init
      nm_nsamp_i <- gsub("^obs","nsamp",nm_i) # name in init
      nm_nfish_ts_i <- gsub("ob$","nfish",nm3_i) # name in rdat$t.series
      nm_nsamp_ts_i <- gsub("ob$","n",nm3_i) # name in rdat$t.series

      nfishi <- setNames(rdat_proj$t.series[yrsi,nm_nfish_ts_i],yrsi)
      nsampi <- setNames(rdat_proj$t.series[yrsi,nm_nsamp_ts_i],yrsi)


      init_p[[nm_nfish_i]] <- nfishi[yrsi]
      init_p[[nm_nsamp_i]] <- nsampi[yrsi]
    }

    ## lenc
    init_obs_lenc_nm <- names(init_b)[grepl("^obs_lenc",names(init_b))]
    key_nm_init_obs_lenc <- setNames(names(lcomp_ob),
                                     paste0("obs_lenc_",gsub("^(D.)([A-Za-z]+)","\\2_D",gsub("^L.","",names(lcomp_ob)))))
    # If the length comp names aren't using the _D suffix, attempt to remove _D suffix from the names in the key
    if(!any(grepl("\\_D$",init_obs_lenc_nm))){
      names(key_nm_init_obs_lenc) <- gsub("\\_D$","",names(key_nm_init_obs_lenc))
    }

    for(i in seq_along(init_obs_lenc_nm)){
      nm_i <- init_obs_lenc_nm[i] # init object name
      nm2_i <- key_nm_init_obs_lenc[[nm_i]] # lcomp_ob name
      nm3_i <- names(key_nm_lcomp_ob)[match(nm2_i,key_nm_lcomp_ob)] # comp.mats name
      xbi <- init_b[[nm_i]]
      #x2i <- lcomp_ob[[nm2_i]]
      x3i <- rdat_proj$comp.mats[[nm3_i]]

      ndigi <- round(median(nchar(gsub("([0-9]*.)([0-9]*)","\\2",xbi))))
      obsi <- apply(x3i,2,function(x){sprintf(paste0("%.",ndigi,"f"), x)})
      attributes(obsi) <- attributes(x3i)
      yrsi <- rownames(obsi)
      init_p[[nm_i]] <- obsi

      # update styr, endyr, yrs, and nyr as necessary
      yrinfoi <- setNames(list(min(yrsi),max(yrsi),yrsi,paste(length(yrsi))),
                          paste0(c("styr","endyr","yrs","nyr"),gsub("obs","",nm_i))
      )
      yrinfoi_nm_is <- names(yrinfoi)[names(yrinfoi)%in%names(init_p)]
      init_p[yrinfoi_nm_is] <- yrinfoi[yrinfoi_nm_is]

      # update nfish and nsamp as necessary
      nm_nfish_i <- gsub("^obs","nfish",nm_i) # name in init
      nm_nsamp_i <- gsub("^obs","nsamp",nm_i) # name in init
      nm_nfish_ts_i <- gsub("ob$","nfish",nm3_i) # name in rdat$t.series
      nm_nsamp_ts_i <- gsub("ob$","n",nm3_i) # name in rdat$t.series

      nfishi <- setNames(rdat_proj$t.series[yrsi,nm_nfish_ts_i],yrsi)
      nsampi <- setNames(rdat_proj$t.series[yrsi,nm_nsamp_ts_i],yrsi)


      init_p[[nm_nfish_i]] <- nfishi[yrsi]
      init_p[[nm_nsamp_i]] <- nsampi[yrsi]
    }

    bam_p <- bam2r(dat_obj=bam$dat,tpl_obj=bam$tpl,cxx_obj=bam$cxx,init=init_p)
  }else{
    bam_p <- NA
  } # end if(!is.null(args_bam2r))

  out_proj$bam_p <- bam_p

##############################################################
### plot stuff ###############################################
##############################################################
  if(plot){
  par(mfrow=c(2,2),mar=c(3,3,1,1),mgp=c(1,0.2,0),tck=-0.01)
  # N
  plot(as.numeric(names(Nsum)),Nsum,type="o",xlab="year")
  abline(v=endyr,lty=2)
  text(x=endyr,y=par("usr")[4]*0.9,labels = "endyr",srt=90,pos=1)
  # R
  plot(as.numeric(names(R)),R,type="o",xlab="year")
  abline(v=endyr,lty=2)
  text(x=endyr,y=par("usr")[4]*0.9,labels = "endyr",srt=90,pos=1)
  # B
  plot(as.numeric(names(Bsum)),Bsum,type="o",xlab="year")
  abline(v=endyr,lty=2)
  text(x=endyr,y=par("usr")[4]*0.9,labels = "endyr",srt=90,pos=1)
  # Ffull
  plot(as.numeric(names(Ffull)),Ffull,type="o",xlab="year")
  abline(v=endyr,lty=2)
  text(x=endyr,y=par("usr")[4]*0.9,labels = "endyr",srt=90,pos=1)

  # Landings and discards
  # Cn.L.long
  p <- ggplot(Cn.L.long,mapping=aes(x=year,y=Ln))+
    geom_area(aes(fill=fleet))+
    theme_bw()+
    scale_fill_brewer(palette="Spectral")+
    stat_summary(fun = sum, geom = "line", size = 1)+
    stat_summary(fun = sum, geom = "point", size = 2)+
    geom_vline(xintercept = endyr, linetype="dashed", size = 0.3)
  p2 <- p + annotate("text",x=endyr, label="endyr\n",y=Inf, angle=90,hjust=1.5,vjust=1)
  print(p2)


  # Cn.D.long
  if(nrow(Cn.D.long)>0){
  p2 <- p %+% Cn.D.long + aes(y=Dn) +
  annotate("text",x=endyr, label="endyr\n",y=Inf, angle=90,hjust=1.5,vjust=1)
  print(p2)
}

  # Cw.L.long
  p2 <- p %+% Cw.L.long + aes(y=Lw) +
  annotate("text",x=endyr, label="endyr\n",y=Inf, angle=90,hjust=1.5,vjust=1)
  print(p2)


  # Cw.D.long
  if(nrow(Cw.D.long)>0){
  p2 <- p %+% Cw.D.long + aes(y=Dw) +
  annotate("text",x=endyr, label="endyr\n",y=Inf, angle=90,hjust=1.5,vjust=1)
  print(p2)
}

  # cpue
    matplot(as.numeric(dimnames(U)[[1]]),U,type="o",xlab="",xlim=c(styr,endyr+nyp),pch=1,ylim=range(c(0,U),na.rm=TRUE))
    # matpoints(as.numeric(rownames(U_p)),U_p,type="o",pch=1)
    legend("topleft",legend=dimnames(U)[[2]],col=1:ncol(U),lty=1:ncol(U),pch=1)
    abline(v=endyr,lty=2)
    text(x=endyr,y=par("usr")[4]*0.9,labels = "endyr",srt=90,pos=1)
  }

##############################################################
### Return results ###########################################
##############################################################
  dimnames(U)[[2]]   <- paste0("U_", dimnames(U)[[2]])
  dimnames(U_p)[[2]] <- paste0("U_", dimnames(U_p)[[2]])

  invisible(out_proj)
}
