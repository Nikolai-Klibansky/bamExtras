#' run_proj
#'
#' Project BAM forward in time to conduct deterministic or stochastic projections.
#' @param rdat BAM output rdat (list) object read in with dget().
#' @param run_bam_args list of arguments passed internally to run_bam if any are provided
#' @param bam2r_args list of arguments passed internally to bam2r if any are provided
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
#' @param L_cur current level of landings
#' @param ages   ages. numeric vector
#' @param N_styr_proj abundance at age in styr_proj. numeric vector
#' @param S_styr_proj spawning stock size at age in styr_proj. Often biomass in mt, sometimes eggs in n. numeric vector
#' @param M      natural mortality rate, at age
#' @param wgt_mt weight at age of population in metric tons (mt). numeric vector
#' @param wgt_L_klb  weight at age of landings in thousand pounds (klb). numeric vector
#' @param wgt_D_klb  weight at age of discards in thousand pounds (klb). numeric vector
#' @param wgt_F_flt_klb list of weight vectors (klb by age) or matrices (by year,age) for each fleet associated with landings or discards. numeric vector or matrix
#' @param len_F_flt_mm list of length vectors (mm by age) or matrices (by year,age) for each fleet associated with landings or discards. numeric vector or matrix
#' @param reprod reproductive contribution at age to SSB. numeric vector
#' @param sel_L  selectivity at age to compute landings. numeric vector
#' @param sel_D  selectivity at age to compute dead discards. numeric vector
#' @param sel_tot  selectivity at age to compute Z (includes landings and discards). numeric vector
#' @param sel_F_flt list of selectivity vectors (by age) or matrices (by year,age) for each fleet associated with landings or discards. numeric vector or matrix
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
#' @param ia_args list of arguments to pass to interim adjustment function. See Details section below.
#' @param obs_error list of optional arguments for incorporating observation error into projections. See Details section below.
#' @param plot logical. Produce plots of extended base model including projected data
#' @param nm_comp_sfx Character vector. Possible suffix added to fleet abbreviation in comps
#' (e.g. c("a","b") in SEDAR 82 gray triggerfish). Used to match comps with catch when projecting
#' comp data when \code{project_bam=TRUE}
#' @param t_series_na_vals Numeric vector of values to consider as NA in t.series object.
#' @param key_lenprob Optional list for specifying cv values for length-at-age,
#' for each set of length compositions. (e.g. key_lenprob = c("D.rHD"="len.cv.L", "L.rGN"="len.cv.L")). Use keyword 'all', as in default to set all
#' values to a particular value in parm.cons.
#' @details
#' \code{ia_args}:
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
#' \code{obs_error}:
#' \itemize{
#' \item{cv_U_sc}{Number to multiply by cv values of abundance indices in the projection period.}
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
#' # and project bam inputs, by adding the bam files with bam2r_args
#' proj_VeSn <- run_proj(rdat_VermilionSnapper,
#' project_bam=TRUE,
#' bam2r_args = list(
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
                     run_bam_args = NULL,
                     bam2r_args = NULL,
                     stochastic = FALSE,
                     pstar = NULL,
                     styr_proj = NULL,
                     nyb_rcn = list(L=3, U=3, R=5, comp=3),
                     nyp_cur = 3,
                     nyp = 5,
                     nm_yr_p=c("endyr","endyr_rec_dev","endyr_rec_phase2","endyr_rec_spr","styr_regs","endyr_proj"),
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
                     ia_args = list(U_nm = NULL),
                     obs_error = list("cv_U_sc" = 0),
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

ia_args_default = list(
  type = "Uprop", yr_U_lag = -(2:4), yr_ia_by = 5, include_yrlim_p=FALSE
  )
ia_args <- modifyList(ia_args_default,ia_args)

if(!is.null(run_bam_args)){
  if(!"return_obj"%in%names(run_bam_args)){
    run_bam_args <- c(run_bam_args,list(return_obj=c("dat","tpl","cxx","rdat")))
  }
  run_bam_out <- do.call(run_bam,run_bam_args)
  rdat <- run_bam_out$rdat
}

### Get parts of the rdat used later in the function and compute stuff
  # if(!is.null(rdat)){
   parms <- rdat$parms
   a.series <- rdat$a.series
   t.series <- rdat$t.series
   parm.cons <- rdat$parm.cons
   sel_age <- rdat$sel.age
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

    sel_age_1 <- sel_age[names(sel_age)%in%c("sel.v.wgted.L","sel.v.wgted.D","sel.v.wgted.tot")]
    sel_age_2 <- sel_age[!names(sel_age)%in%names(sel_age_1)]

    yrs_L_b <- paste((endyr-(nyb_rcn$L-1)):endyr) # years of recent landings (i.e. from the base model)
    yrs_R_b <- paste((endyr-(nyb_rcn$R-1)):endyr)
    R_b <- t.series[yrs_R_b,"recruits"]
    R_b_gm <- bamExtras::geomean2(R_b)

    # Identify F-at-age for each fleet in endyr
    # L = (F/Z)*N*(1-exp(-Z))
    # L/(N*(1-exp(-Z))) = F/Z
    # L*Z/(N*(1-exp(-Z))) = F
    # F = L*Z/(N*(1-exp(-Z)))
    Cn <- CLD.est.mats[names(CLD.est.mats)[grepl("^(Ln|Dn)((?!total).)*$",names(CLD.est.mats),perl=TRUE)]]
    names(Cn) <- gsub("^([LD])(n)(.*)","Cn.\\1\\3",names(Cn))
    Cw <- CLD.est.mats[names(CLD.est.mats)[grepl("^(Lw|Dw)((?!total).)*$",names(CLD.est.mats),perl=TRUE)]]
    names(Cw) <- gsub("^([LD])(w)(.*)","Cw.\\1\\3",names(Cw))
    C_ts_cv <- t.series[paste(yb),grepl("^cv.[DL]",names(t.series))] # cvs associated with components of the catch (removals) during the base years
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
    # names(F_flt) <- paste("F",gsub("([LD])(n)(.*)","\\1\\3",names(F_flt)),sep=".")
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
      LD_abb <- gsub(".pr$","",names(t.series)[grepl("^[LD].*.pr$",names(t.series),perl=TRUE)]) # landings abbreviations
      LD_abb1_ts <- paste0("F.",gsub("^(D.)(.*)","\\2.D",gsub("^L.","",LD_abb))) # pattern 1 (typical)
      LD_abb2_ts <- paste0("F.",LD_abb) # pattern 2 (less typical but preferred)
      if(all(LD_abb1_ts%in%names(t.series))){
        LD_abb_ts <- LD_abb1_ts
      }else if(all(LD_abb2_ts%in%names(t.series))){
        LD_abb_ts <- LD_abb2_ts
      }else {
        warning(paste0("Can't match Fs. Neither set 1 (",paste(LD_abb1_ts,collapse=", "),
                       ") nor set 2 (",paste(LD_abb2_ts,collapse=", "),
                       ") column names can all be found in t.series"))
      }
      Fsum_flt <- t.series[paste(yb),LD_abb_ts] # F time series for each fleet
      names(Fsum_flt) <- paste0("F.",LD_abb)
      Fsum <- rowSums(Fsum_flt) # Should be equal to t.series$Fsum

      # Compute selectivity for each fleet associated with an F value.
      # NOTE: Computing selectivities is somewhat more reliable than trying to find
      # them in the rdat, in instances when selectivities from one fleet are used
      # for multiple fleets (e.g. when the headboat selectivity is used for the MRIP landings)
      sel_F_flt <- lapply(names(F_flt),function(x){
        a <- pmax(pmin(F_flt[[x]]/Fsum_flt[,x],1),0)
        a[is.na(a)] <- 0 # Replace NA with zero (NaN will occur when Fsum_flt values are zero, as in years when fleets do not have any landings)
        a
      })
      names(sel_F_flt) <- gsub("^F.","sel.",names(F_flt))

      F_ybgm_flt <- apply(tail(Fsum_flt,nyb_rcn$L),2,bamExtras::geomean2) # geomean F by fleet during the last nyb_rcn years of the base model
      F_prop_flt <- F_ybgm_flt/sum(F_ybgm_flt) # Proportion of F attributed to each fleet at the end of the assessment
      # This is just F_ybgm_flt multiplied by the selectivity in the endyr of the base model for each fleet
      # In the bam tpl, these vectors are named F_end although they are summed by landings (F_end_L), discards (F_end_D), or both (F_end)
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
      sel_wgted_F_flt <- F_ybgm_flt_a/F_end_apex     # by extension, scale selectivity at age for each fleet
      # Note that: round((sel_wgted_L+sel_wgted_D) == rowSums(sel_wgted_F_flt)
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
    U_ts_pr <- t.series[paste(yb),grepl("^U.*pr$",names(t.series)),drop=FALSE]
    U_ts_ob <- t.series[paste(yb),grepl("^U.*ob$",names(t.series)),drop=FALSE]
    U_ts_cv <- t.series[paste(yb),grepl("^cv.U",names(t.series)),drop=FALSE] # cvs associated with U during the base years
    # names(U_ts_pr) <- gsub("^(U.)(.*)(.)(pr)$","\\1\\4.\\2",names(U_ts_pr))
    U_abb <- gsub("^U.|.pr$","",names(U_ts_pr))
    names(U_ts_pr) <- names(U_ts_ob) <- U_abb
    # Find q values in rdat
    # Are the q values provided in t.series?
    #q_nm <- paste0("q.",U_abb)
    q_nm_tsY <- names(rdat$t.series)[grepl(paste0("^q..*(",paste(U_abb,collapse="|"),")$"),names(rdat$t.series))]#q_nm[q_nm%in%names(t.series)]
    U_abb_tsY <- U_abb[unlist(lapply(U_abb,function(x){any(grepl(x,q_nm_tsY))}))]
    U_abb_tsN <- U_abb[!U_abb%in%U_abb_tsY]
    # q_nm_tsN <- q_nm[!q_nm%in%names(t.series)]
    #q_nm_tsNpcY <- names(parm.cons)[grepl(paste0("^log.q..*(",paste(U_abb,collapse="|"),")$"),names(parm.cons))]#paste0("log.",q_nm_tsN)[paste0("log.",q_nm_tsN)%in%names(parm.cons)]
    # q_nm_tsNpcN <- q_nm[which((!q_nm%in%names(t.series))&(!paste0("log.",q_nm)%in%names(parm.cons)))]
    if(length(q_nm_tsY)>0){
      q_mn <- t.series[paste(yb),q_nm_tsY,drop=FALSE]  # This is really a kind of mean q, since q_rate and q_DD_mult might scale it to compute U
    }
    # If any of the q values are not in t.series, look in parms.cons
    if(length(U_abb_tsN)>0){
      message(paste0("q values for ",paste(U_abb_tsN,collapse=", "), " index not found in t.series"))
      q_nm_tsNpcY <- names(parm.cons)[grepl(paste0("^log.q..*(",paste(U_abb_tsN,collapse="|"),")$"),names(parm.cons))]
      if(length(q_nm_tsNpcY)>0){
        message(paste0("time invariant values ",paste(q_nm_tsNpcY,collapse=", "), " found in parm.cons will be exp transformed and used instead."))
        U_abb_tsNpcY <- gsub(paste0("^(.*)(",paste(U_abb,collapse="|"),")"),"\\2",q_nm_tsNpcY)
        q_tsNpcY <- setNames(exp(parm.cons[q_nm_tsNpcY][8,]),U_abb_tsNpcY)
        # names(q_tsNpcY) <- gsub("^log.q.","",names(q_tsNpcY))
        for(nm_i in names(q_tsNpcY)){
          a <- U_ts_pr[,nm_i]
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
    q_rate <- q_mn*0+1
    q <- q_mn*NA # Initialize values for final q value (q_mn * q_rate * q_DD_mult) multiplied by N to compute CPUE

    sel_U <- sel_age_2[which(gsub("sel.[vm].","",names(sel_age_2))%in%U_abb)]
    names(sel_U) <- paste0("sel.U.",U_abb)

    for(i in 1:length(sel_U)){
      if(is.vector(sel_U[[i]])){
        sel_U[[i]] <- matrix(sel_U[[i]], nrow=length(yb),ncol=length(sel_U[[i]]),byrow=TRUE,dimnames=list(yb,names(sel_U[[i]])))
      }
      sel_U_i <- sel_U[[i]]
      sel_nm_i <- names(sel_U)[i]
      sel_U_abb_i <- gsub("sel.U.","",sel_nm_i)

      q_nm_i <- paste0("q.",sel_U_abb_i)
      q_mn_i <- q_mn[paste(yb),q_nm_i,drop=FALSE]
      q_rate_nm_i <- paste0("q.",sel_U_abb_i,".rate.mult")
      if(q_rate_nm_i%in%names(t.series)){
        q_rate_i <- t.series[paste(yb),q_rate_nm_i,drop=FALSE]
      }else{
        q_rate_i <- q_mn_i*0+1
        message(paste("message:",q_rate_nm_i,"not found in t.series. Setting q rate multiplier values to 1.\n"))
      }
      q_rate[,q_nm_i] <- q_rate_i

      q_i <- q_mn_i*q_rate_i*q_DD_mult
      q[,q_nm_i] <- q_i

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
      U_ts_ob_i <- U_ts_ob[,sel_U_abb_i]
      U_ts_pr_i <- U_ts_pr[,sel_U_abb_i]

      # Compute the sums of squared deviations to determine whether the indices are
      # in numbers or in weight. Indices in the appropriate units are then added
      # to the the appropriate objects
      SSUi <- unlist(lapply(list(U_n_i=U_n_i,U_w_i=U_w_i),function(x){sum((U_ts_pr_i-x)^2,na.rm=TRUE)}))
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

      ndigi <- round(median(nchar(gsub("([0-9]*.)([0-9]*)","\\2",na.omit(U_ts_ob_i)))))
      xpi <- round(get(names(SSUi)[which.min(SSUi)]),ndigi)


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
    acomp_b_ob <- comp.mats[grepl("^acomp.*.ob$",names(comp.mats))]
    names(acomp_b_ob) <- local({
      a <- gsub("^(acomp.)(.*)(.ob)$","\\2",names(acomp_b_ob))
      b <- gsub("(.*)(.D)$","D.\\1",a) # Convert .D suffix to D. prefix
      c <- gsub("^([cr])(.*)","L.\\1\\2",b) # Add L. prefix to comps starting with c or r (landings not discards. won't modify fishery independent survey comps)
      c
    })
    # predicted
    acomp_b_pr <- comp.mats[grepl("^acomp.*.pr$",names(comp.mats))]
    names(acomp_b_pr) <- local({
      a <- gsub("^(acomp.)(.*)(.pr)$","\\2",names(acomp_b_pr))
      b <- gsub("(.*)(.D)$","D.\\1",a) # Convert .D suffix to D. prefix
      c <- gsub("^([cr])(.*)","L.\\1\\2",b) # Add L. prefix to comps starting with c or r (landings not discards. won't modify fishery independent survey comps)
      c
    })

    ## length compositions from the base years
    # observed
    lcomp_b_ob <- comp.mats[grepl("^lcomp.*.ob$",names(comp.mats))]
    names(lcomp_b_ob) <- local({
      a <- gsub("^(lcomp.)(.*)(.ob)$","\\2",names(lcomp_b_ob))
      b <- gsub("(.*)(\\.D)$","D.\\1",a) # Convert .D suffix to D. prefix
      c <- gsub("^([a-zA-Z]{2})(D)","D.\\1\\2",b) # Add D. prefix if third letter is D
      d <- gsub("^([cr])(.*)","L.\\1\\2",c) # Add L. prefix to comps starting with c or r (landings not discards. won't modify fishery independent survey comps)
      d
    })
    # predicted
    lcomp_b_pr <- comp.mats[grepl("^lcomp.*.pr$",names(comp.mats))]
    names(lcomp_b_pr) <- local({
      a <- gsub("^(lcomp.)(.*)(.pr)$","\\2",names(lcomp_b_pr))
      b <- gsub("(.*)(\\.D)$","D.\\1",a) # Convert .D suffix to D. prefix
      c <- gsub("^([a-zA-Z]{2})(D)","D.\\1\\2",b) # Add D. prefix if third letter is D
      d <- gsub("^([cr])(.*)","L.\\1\\2",c) # Add L. prefix to comps starting with c or r (landings not discards. won't modify fishery independent survey comps)
      d
    })

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

    for(i in 1:length(lcomp_b_pr)){
      nmi <- names(lcomp_b_pr)[i]
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

      lc_i <- lcomp_b_pr[[i]]
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

    names(lenprob) <- names(lcomp_b_pr)


    ## Recompute predicted age and length comps
    # Note: These should be the same or very similar to age and length comps computed by BAM.
    #  But these are recomputed for the base years as a check since the computation
    #  will be used in the projection years
    # acomp_b_pr_2 <-  list()
    # for(nm_i in names(acomp_b_pr)){
    #   yrs_ac_i <- rownames(acomp_b_pr[[nm_i]])
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
    # for(nm_i in names(lcomp_b_pr)){
    #   yrs_lc_i <- rownames(lcomp_b_pr[[nm_i]])
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
    ntrip <- t.series[paste(yb),grepl("^([la]comp).*(.n)$",names(t.series))]
    # ntrip[ntrip==-99999] <- NA
    names(ntrip) <- local({
      a <- gsub("(.n)$","",names(ntrip))
      b <- gsub("^([al]comp.)(.*)(.)(D)$","\\1D.\\2",a) # Convert .D suffix to D. prefix
      c <- gsub("^([al]comp.)([cr])(.*)","\\1L.\\2\\3",b) # Add L. prefix to comps starting with c or r (landings not discards. won't modify fishery independent survey comps)
      c
    })
    # number of fish in compositions
    nfish <- t.series[paste(yb),grepl("^([la]comp).*(.nfish)$",names(t.series))]
    # nfish[nfish==-99999] <- NA
    names(nfish) <- local({
      a <- gsub("(.nfish)$","",names(nfish))
      b <- gsub("^([al]comp.)(.*)(.)(D)$","\\1D.\\2",a) # Convert .D suffix to D. prefix
      c <- gsub("^([al]comp.)([cr])(.*)","\\1L.\\2\\3",b) # Add L. prefix to comps starting with c or r (landings not discards. won't modify fishery independent survey comps)
      c
    })

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
  # sel_F_flt_p <- lapply(1:length(Fsum_flt),function(x){emypa})
  # names(sel_F_flt_p) <- gsub("^F.","sel.",names(Fsum_flt))

  U_a_p <-  emuypa    # predicted cpue at age by fleet
  unit_U_p <- emuypa # weights associated with cpue at age by fleet
  q_p <- emuyp      # catchability (q)

  F_flt_p <- emfypa

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
  C_ts_cv_p <- local({
    a <- matrix(apply(C_ts_cv,2,function(x){rep(bamExtras::geomean2(tail(x,nyb_rcn$L)),nyp)}),
                nrow=nyp,dimnames=list(paste(yp),names(C_ts_cv)))
    # rownames(a) <- paste(yp)
    a
  })

  U_ts_cv_p <- local({
    a <- matrix(apply(U_ts_cv,2,function(x){rep(bamExtras::geomean2(tail(x,nyb_rcn$U)),nyp)}),
                nrow=nyp,dimnames=list(paste(yp),names(U_ts_cv)))
    # rownames(a) <- paste(yp)
    a
  })

  ## sel during projection years (by py, age, fleet)
  # for any source of F (landings or discards)
  # Fill with weighted selectivity from last year of base (endyr)
  sel_F_flt_p <- lapply(1:length(sel_F_flt),function(i){
    nm_i <- names(sel_F_flt)[i]
    sel_endyr_i <- sel_wgted_F_flt[,gsub("^sel.","F.",nm_i)]
    # sel_endyr_i <- sel_F_flt[[nm_i]][paste(endyr),] # Gets unweighted selectivities
    matrix(sel_endyr_i, nrow=nyp, ncol=nages, byrow=TRUE, dimnames = dimnames(emypa))
  })
  names(sel_F_flt_p) <- gsub("^F.","sel.",names(Fsum_flt))

  # sel during projection years (by py, age, fleet) for any source of indices of abundance (cpue)
  # Fill with selectivity from last year of base (endyr)
  sel_U_p <- lapply(1:length(sel_U),function(i){
    nm_i <- names(sel_U)[i]
    sel_endyr_i <- sel_U[[nm_i]][paste(endyr),]
    matrix(sel_endyr_i, nrow=nyp, ncol=nages, byrow=TRUE, dimnames = dimnames(emypa))
  })
  names(sel_U_p) <- names(sel_U)

  # sel during projection years (by py, age, fleet) for any miscellaneous data sources
  # Fill with selectivity from last year of base (endyr)
  if(length(sel_misc)>0){
  sel_misc_p <- lapply(1:length(sel_misc),function(i){
    nm_i <- names(sel_misc)[i]
    sel_endyr_i <- sel_misc[[nm_i]][paste(endyr),]
    matrix(sel_endyr_i, nrow=nyp, ncol=nages, byrow=TRUE, dimnames = dimnames(emypa))
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
  q_p[paste(yp),] <- rep(unlist(q[paste(endyr),paste0("q.",colnames(q_p))]),each=nyp)

  # compositions
  acomp_p <- list()
  for(i in names(acomp_b_pr)){
  ac_i <- acomp_b_pr[[i]]
  acomp_p[[i]] <- matrix(NA,nrow=nyp,ncol=ncol(ac_i),dimnames=list("year"=yp,"age"=colnames(ac_i)))
  }

  lcomp_p <- list()
  for(i in names(lcomp_b_pr)){
    lc_i <- lcomp_b_pr[[i]]
    lcomp_p[[i]] <- matrix(NA,nrow=nyp,ncol=ncol(lc_i),dimnames=list("year"=yp,"lenbin"=colnames(lc_i)))
  }

  ntrip_p <- local({
    a <- ceiling(apply(tail(ntrip,nyb_rcn$comp),2,function(x){bamExtras::geomean2(x)}))
    matrix(a,nrow=nyp,ncol=length(a),dimnames=list(year=paste(yp),fleet=names(a)),byrow=TRUE)
  })

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
      yr_ia <- paste(seq(yrlim_p[1],yrlim_p[2],by=ia_args$yr_ia_by))
      if(!ia_args$include_yrlim_p){
        yr_ia <- yr_ia[!yr_ia%in%yrlim_p]
      }

      if(!is.null(ia_args$U_nm)&yp_i%in%yr_ia){
        F_adj <- F_adjust(data=U,
                          U_nm=ia_args$U_nm,
                          yr_ref = paste(endyr),
                          yr = paste(sort(as.numeric(yp_i)+ia_args$yr_U_lag)),
                          type = ia_args$type
                          )
        message(paste("For year =",yp_i,"full F was adjusted by", F_adj,"based on the",ia_args$U_nm, "index"))

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
        sel_F_flt_p[[j]]*Ffull_p
      })
      names(F_flt_p) <- names(Fsum_flt)

      # 2023-10-16 This isn't quite right. Need to figure out how to compute Fsum
      # for the projections
      # Fage_p <- Reduce("+",F_flt_p)
      # Fsum_p <- rowSums(Fage_p)

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
      B_p[i] <- sum(N_p[i,]*wgt_mt) # biomass of fish at-age at that start of the year
      Nmdyr_p[i,] <- N_p[i,]*(exp(-1.0*Z_p[i,]*0.5))        # number of fish at-age at mid-year
      Nspwn_p[i,] <- N_p[i,]*(exp(-1.0*Z_p[i,]*spawn_time)) # number of fish at-age at spawning time
      S_p[i] <- sum(Nspwn_p[i,]*reprod)                     # spawning stock at-age at spawning time

      ## cpue
      # By fleet
      for(j in 1:length(U_a_p)){
        nm_j <- names(U_a_p)[j]
        NU_p_ji <- Nmdyr_p[yp_i,]*sel_U_p[[paste0("sel.U.",nm_j)]][yp_i,]*unit_U_p[[nm_j]][yp_i,]
        NU_p[[nm_j]][yp_i,] <- NU_p_ji
        U_a_p_ji <- NU_p_ji*q_p[yp_i,nm_j]
        # Add observation error to cpue (has no effect if obs_error$cv_U_sc = 0)
        U_p_ji <- local({
          a <- sum(U_a_p_ji)
          cv_nm_j <- paste0("cv.U.",nm_j)
          U_ts_cv_p_ij <- U_ts_cv_p[yp_i,cv_nm_j]*obs_error$cv_U_sc
          as.numeric(lnorm_vector_boot(a,U_ts_cv_p_ij))
        })
        U_a_p[[nm_j]][yp_i,] <- (U_a_p_ji/sum(U_a_p_ji))*U_p_ji
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
      U <- local({
        a <- rbind(U,U_p[yp_i,,drop=FALSE])
        apply(a,2,function(x){x/mean(x,na.rm=TRUE)}) # standardize indices to mean of 1
      })

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

  #### Post-projection computations
  R_p <- N_p[,1] # Get recruitment vector from N-at-age matrix
  #U_p <- as.data.frame(lapply(U_a_p,rowSums)) # Compute cpue time series (summing across ages)

  ## Compute predicted age and length comps during projection years
  for(nm_i in names(acomp_p)){
    yrs_acomp_i <- rownames(acomp_b_ob[[nm_i]])
    yrs_acomp_p_i <- rownames(acomp_p[[nm_i]])
    # If the last year of comps is within the last nyb_rcn$comp of the assessment, then project them.
    if(max(as.numeric(yrs_acomp_i))%in%tail(as.numeric(yb),nyb_rcn$comp)){
    if(grepl("^[LD]",nm_i)){ # If length comps are associated with catch (landings or discards), use numbers from landings
      # Remove nm_comp_sfx from nm_i, if specified
      nm_Cn_i <- nm_i
      if(nm_comp_sfx[1]!=""){
        pat_i <- paste0(nm_comp_sfx,collapse="|")
        nm_Cn_i <- gsub(paste0("(.*)(",pat_i,")$"),"\\1",nm_i)
      }
      n_i <- Cn_p[[paste0("Cn.",nm_Cn_i)]][yrs_acomp_p_i,,drop=FALSE]
    }else if(nm_i%in%names(NU_p)){
      n_i <- NU_p[[nm_i]][yrs_acomp_p_i,,drop=FALSE] # If not, see if it's associated with an index (e.g. fishery independent) and use numbers associated with the index
    }else if(nm_i%in%names(Nmisc_p)){
      n_i <- Nmisc_p[[nm_i]][yrs_acomp_p_i,,drop=FALSE] # If not, there should be an N-at-age matrix in Nmisc that matches
    }else{
      warning(paste0(nm_i, " from acomp_p doesn't appear to match any sources of F, indices, or other data sources associated with a selectivity."))
    }
    P_n_i <- t(apply(n_i,1,function(x){x/sum(x)}))
    agebins_agec_i <- dimnames(acomp_b_ob[[nm_i]])[[2]]
    # Make sure the agebins match the observed age comps and apply a plus group if necessary
    acomp_p_i <- local({
      a <- t(apply(P_n_i,1,function(x){colSums(x*age_error)}))
      b <- a[,agebins_agec_i]
      bplusgroup <- tail(agebins_agec_i,1)
      b[,bplusgroup] <- b[,bplusgroup]+rowSums(a[,as.numeric(dimnames(a)[[2]])>as.numeric(bplusgroup)])
      b
    })
    acomp_p[[nm_i]] <- acomp_p_i
    }
  }

  for(nm_i in names(lcomp_p)){
    yrs_lcomp_i <- rownames(lcomp_b_ob[[nm_i]])
    yrs_lcomp_p_i <- rownames(lcomp_p[[nm_i]])
    # If the last year of comps is within the last nyb_rcn$comp of the assessment, then project them.
    if(max(as.numeric(yrs_lcomp_i))%in%tail(as.numeric(yb),nyb_rcn$comp)){
      lenprob_i <- lenprob[[nm_i]]
      if(grepl("^[LD]",nm_i)){ # If length comps are associated with catch (landings or discards), used numbers from landings
        n_i <- Cn_p[[paste0("Cn.",nm_i)]][yrs_lcomp_p_i,,drop=FALSE]
      }else if(nm_i%in%names(NU_p)){
        n_i <- NU_p[[nm_i]][yrs_lcomp_p_i,,drop=FALSE] # If not, see if it's associated with an index (e.g. fishery independent) and use numbers associated with the index
      }else if(nm_i%in%names(Nmisc_p)){
        n_i <- Nmisc_p[[nm_i]][yrs_lcomp_p_i,,drop=FALSE] # If not, their should be an N-at-age matrix in Nmisc that matches
      }else{
        warning(paste0(nm_i, "from lcomp_p doesn't appear to match any sources of F, indices, or other data sources associated with a selectivity."))
      }
      P_n_i <- t(apply(n_i,1,function(x){x/sum(x)}))
      lenbins_lenc_i <- dimnames(lcomp_b_ob[[nm_i]])[[2]]
      # Make sure the lenbins match the observed length comps and apply a plus group if necessary
      lcomp_p_i <- local({
        a <- t(apply(P_n_i,1,function(x){colSums(x*lenprob_i)}))
        b <- a[,lenbins_lenc_i]
        bplusgroup <- tail(lenbins_lenc_i,1)
        b[,bplusgroup] <- b[,bplusgroup]+rowSums(a[,as.numeric(dimnames(a)[[2]])>as.numeric(bplusgroup)])
        b
      })
      lcomp_p[[nm_i]] <- lcomp_p_i
    }
  }
  # # Project cvs for catch and cpue
  # C_ts_cv_p <- local({
  #   a <- matrix(apply(C_ts_cv,2,function(x){rep(bamExtras::geomean2(tail(x,nyb_rcn$L)),nyp)}),
  #               nrow=nyp,dimnames=list(paste(yp),names(C_ts_cv)))
  #   # rownames(a) <- paste(yp)
  #   a
  # })
  #
  # U_ts_cv_p <- local({
  #   a <- matrix(apply(U_ts_cv,2,function(x){rep(bamExtras::geomean2(tail(x,nyb_rcn$U)),nyp)}),
  #               nrow=nyp,dimnames=list(paste(yp),names(U_ts_cv)))
  #   # rownames(a) <- paste(yp)
  #   a
  # })

# Add projected values to base values
  ybp <- c(yb,yp)
  nybp <- length(ybp)
  N <- rbind(N,N_p)
  logRdev <- c(logRdev,logRdev_p)
  R <- c(R,R_p)
  S <- c(S,S_p)
  # U <- rbind(U,U_p) # U is updated above
  U_a <- local({
    a <- lapply(seq_along(U_a),function(i){rbind(U_a[[i]],U_a_p[[i]])})
    names(a) <- names(U_a)
    a
  })
  C_ts_cv <- rbind(C_ts_cv,C_ts_cv_p)
  U_ts_cv <- rbind(U_ts_cv,U_ts_cv_p)
  nfish <- rbind(nfish,nfish_p)
  ntrip <- rbind(ntrip,ntrip_p)

  lcomp <- local({
    a <- lapply(seq_along(names(lcomp_b_ob)),function(x){
      nm_x <- names(lcomp_b_ob)[x]
      ax <- rbind(lcomp_b_ob[[x]],lcomp_p[[x]])
      ax[complete.cases(ax),,drop=FALSE]
    })
    names(a) <- names(lcomp_b_ob)
    a
  })

  acomp <- local({
    a <- lapply(seq_along(names(acomp_b_ob)),function(x){
      nm_x <- names(acomp_b_ob)[x]
      ax <- rbind(acomp_b_ob[[x]],acomp_p[[x]])
      ax[complete.cases(ax),,drop=FALSE]
    })
    names(a) <- names(acomp_b_ob)
    a
  })

  ##Fsum <-  c(Fsum, Fsum_p)
  Ffull <- c(Ffull, Ffull_p)

  for(nm_i in names(Cn)){
    Cn[[nm_i]] <- rbind(Cn[[nm_i]],Cn_p[[nm_i]])
  }
  for(nm_i in names(Cw)){
    Cw[[nm_i]] <- rbind(Cw[[nm_i]],Cw_p[[nm_i]])
  }

  } # end if nyp>0

    B <- wgt_mt*N
    Nsum <- rowSums(N)
    Bsum <- rowSums(B)


  # Values by fleet  in long format (for ggplot)
  # Landings (n)
  Cn.L <- local({
    a <- as.data.frame(lapply(Cn[grepl("^Cn.L.",names(Cn))],rowSums))
    b <- cbind("year"=as.numeric(rownames(a)),a)
    colnames(b) <- gsub("Cn.L.","",colnames(b))
    gather(data=b,key=fleet,value=Ln,names(b)[names(b)!="year"],factor_key = TRUE)

  })
  # Discards (n)
  Cn.D <- local({
    a <- as.data.frame(lapply(Cn[grepl("^Cn.D.",names(Cn))],rowSums))
    b <- cbind("year"=as.numeric(rownames(a)),a)
    colnames(b) <- gsub("Cn.D.","",colnames(b))
    gather(data=b,key=fleet,value=Dn,names(b)[names(b)!="year"],factor_key = TRUE)

  })

  # Landings (w)
  Cw.L <- local({
    a <- as.data.frame(lapply(Cw[grepl("^Cw.L.",names(Cw))],rowSums))
    b <- cbind("year"=as.numeric(rownames(a)),a)
    colnames(b) <- gsub("Cw.L.","",colnames(b))
    gather(data=b,key=fleet,value=Lw,names(b)[names(b)!="year"],factor_key = TRUE)

  })
  # Discards (w)
  Cw.D <- local({
    a <- as.data.frame(lapply(Cw[grepl("^Cw.D.",names(Cw))],rowSums))
    b <- cbind("year"=as.numeric(rownames(a)),a)
    colnames(b) <- gsub("Cw.D.","",colnames(b))
    gather(data=b,key=fleet,value=Dw,names(b)[names(b)!="year"],factor_key = TRUE)

  })

  Cn.L.tot <- tapply(Cn.L[,"Ln"],Cn.L[,"year"],sum)
  Cw.L.tot <- tapply(Cw.L[,"Lw"],Cw.L[,"year"],sum)

  if(nrow(Cn.D)>0){
  Cn.D.tot <- tapply(Cn.D[,"Dn"],Cn.D[,"year"],sum)
  Cw.D.tot <- tapply(Cw.D[,"Dw"],Cw.D[,"year"],sum)
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
    Cn.L=Cn.L,
    Cw.L=Cw.L,
    Cn.D=Cn.D,
    Cw.D=Cw.D
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


#############################################################
###### Extend data inputs and build projected dat file ######
#############################################################
  # Identify the parts of init that need to be changed
  # Identify the new parts and match them with the parts that need to change
  if(project_bam&!is.null(bam2r_args)){
    bam2r_args <- modifyList(x=as.list(formals(bam2r)),val=bam2r_args) # Add user supplied arguments to defaults
    bam <- do.call(bam2r,bam2r_args)
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
    for(nm_i in init_obs_L_nm){
      abb_i <- gsub("^(obs_L)(_)(.*)","\\3",nm_i)
      xbi <- setNames(as.numeric(init_b[[nm_i]]),names(init_b[[nm_i]]))
      yrsbi <- names(xbi)
      xni <- rowSums(Cn[[paste0("Cn.L.",abb_i)]])
      xkwi <- rowSums(Cw[[paste0("Cw.L.",abb_i)]]) # weights are usually 1000 lbs in CLD.est.mats
      xkni <- xni/1000
      xwi <- xkwi*1000
      xcvi <- init_b[[paste0("obs_cv_L_",abb_i)]]
      # Compute the sums of squared deviations to determine whether landings are
      # in numbers or in weight. Landings in the appropriate units are then added
      # to the observed landings to add to the new dat file.
      SSxbi <- unlist(lapply(list(xni=xni,xwi=xwi,xkni=xkni,xkwi=xkwi),function(x){sum((x[yrsbi]-xbi)^2)}))
      ndigi <- round(median(nchar(gsub("([0-9]*.)([0-9]*)","\\2",xbi))))
      # ndigcvi <- round(median(nchar(gsub("([0-9]*.)([0-9]*)","\\2",xcvi))))
      xpi <- round(get(names(SSxbi)[which.min(SSxbi)])[paste(yp)],ndigi)
      # If the last year of xbi was the last year of the base model, then project it forward
      if(paste(endyr)%in%names(xbi)){
        xi <- c(xbi,xpi)
      }else{
        xi <- xbi
      }
      init_p[[nm_i]] <- setNames(paste(xi),names(xi))
      # reset appropriate _endyr value to make sure it agrees with the projected data (some endyr values might get projected even if the data doesn't)
      endyr_nm_i <- paste0("endyr_L_",abb_i)
      if(endyr_nm_i%in%names(init_p)){init_p[[endyr_nm_i]] <- tail(names(xi),1)}

      init_p[[paste0("obs_cv_L_",abb_i)]] <- local({
        # setNames(paste(round(C_ts_cv[,paste0("cv.L.",abb_i)],ndigcvi)),rownames(C_ts_cv))[names(xi)]
        a <- setNames(as.numeric(init_p[[nm_i]])*NA,names(init_p[[nm_i]]))
        b <- init_b[[paste0("obs_cv_L_",abb_i)]]
        a[names(b)] <- b
        c <- bamExtras::geomean2(tail(b,nyb_rcn$L))
        d <- paste(yp)[paste(yp)%in%names(a)]
        a[d] <- c
        a
      })

      ## set_log_dev_vals_F_L
      set_log_dev_vals_F_L_nm_i <- paste0("set_log_dev_vals_F_L_",abb_i)
      if(set_log_dev_vals_F_L_nm_i%in%names(init_b)){
        init_p[[set_log_dev_vals_F_L_nm_i]] <- setNames(rep("0.0",length(xi)),names(xi))
      }

    }

    ## discards
    # Note: Assumed to be in numbers
    init_obs_released_nm <- names(init_b)[grepl("^obs_released",names(init_b))]
    for(nm_i in init_obs_released_nm){
      abb_i <- gsub("^(obs_released)(_)(.*)","\\3",nm_i)
      xbi <- setNames(as.numeric(init_b[[nm_i]]),names(init_b[[nm_i]]))
      yrsbi <- names(xbi)
      # If the name of the time series in the init object matches a name reported in the model output (Cn),
      # use the corresponding time series from Cn for the init object in the projection
      # NOTE: These if else statements were included as a workaround to accommodate the Black Sea Bass model
      if(paste0("Cn.D.",abb_i)%in%names(Cn)){
        DMi_nm <- paste0("set_Dmort_",abb_i)
        xni_dead <- rowSums(Cn[[paste0("Cn.D.",abb_i)]])
        if(DMi_nm%in%names(init_b)){
          DMi <- as.numeric(init_b[[paste0("set_Dmort_",abb_i)]])
        }else{
          DMi <- mean(xni_dead[yrsbi]/1000/xbi,na.rm=TRUE)
        }
        xni <- xni_dead*1/DMi # Convert to released (live discards)
        xkni <- xni/1000
      }else{ # ..but if not, compute the init object in the projection as a function init object i
        # from the base and the most similar time series from Cn
        ft_i <- gsub("^([cr]).*","\\1",abb_i) # identify fleet type i should be "c" or "r"
        # sum all discard matrices in Cn with fleet_type_i and sum by year across ages
        Cn_ts_ft_i <- rowSums(Reduce(sum,Cn[names(Cn)[grepl(paste0("^(Cn.D.",ft_i,").*"),names(Cn))]]))
        # identify init time series with ft_i
        init_b_released_ft_i <- local({
          a <- init_b[init_obs_released_nm[grepl(paste0("^(obs_released_",ft_i,").*"),init_obs_released_nm)]]
          c <- lapply(a,function(x){
          b <- setNames(rep(0,nyb),paste(yb))
          b[names(x)] <- as.numeric(x)
          b
        }
        )
          as.data.frame(c)
      })
        # Estimate total released (i.e. live discards) during the base model years
          Ckn_ts_ft_i <- Cn_ts_ft_i[paste(yb)]/1000
          released_tot_obs_kn_i <- rowSums(init_b_released_ft_i)
          Dm_ts_i <- Ckn_ts_ft_i/released_tot_obs_kn_i # Annual Dmort
          Dm_ts_i[!is.finite(Dm_ts_i)] <- NA
          Cn_released_b_ft_i <- 1000*Ckn_ts_ft_i*1/Dm_ts_i

        # Compute observed released (live discards) as a proportion of total released in ft_i
        init_b_released_prop_ft_i <- init_b_released_ft_i/rowSums(init_b_released_ft_i)

        # Estimate xni_b from values that can be used in projections
        xni_b <- Cn_released_b_ft_i*init_b_released_prop_ft_i[,nm_i]
        xni_b[!is.finite(xni_b)] <- 0
        xkni_b <- xni_b/1000

        # Estimate xni for projection years
        Cn_released_p_ft_i <- (Cn_ts_ft_i[paste(yp)])*(1/mean(Dm_ts_i[paste(tail(yb,nyb_rcn$L))]))
        xni_p <- Cn_released_p_ft_i*mean(init_b_released_prop_ft_i[paste(tail(yb,nyb_rcn$L)),nm_i])
        xkni_p <- xni_p/1000
        xni <- c(xni_b,xni_p)
        xkni <- c(xkni_b,xkni_p)
      }

      xcvi <- init_b[[paste0("obs_cv_D_",abb_i)]]
      SSxbi <- unlist(lapply(list(xni=xni,xkni=xkni
                                  #xkwi=xkwi,xwi=xwi
                                  ),function(x){sum((x[yrsbi]-xbi)^2)}))
      ndigi <- round(median(nchar(gsub("([0-9]*.)([0-9]*)","\\2",xbi))))
      # ndigcvi <- round(median(nchar(gsub("([0-9]*.)([0-9]*)","\\2",xcvi))))
      xpi <- round(get(names(SSxbi)[which.min(SSxbi)])[paste(yp)],ndigi)
      # If the last year of xbi was the last year of the base model, then project it forward
      if(paste(endyr)%in%names(xbi)){
        xi <- c(xbi,xpi)
      }else{
        xi <- xbi
      }
      init_p[[nm_i]] <- setNames(paste(xi),names(xi))
      # reset appropriate _endyr value to make sure it agrees with the projected data (some endyr values might get projected even if the data doesn't)
      endyr_nm_i <- paste0("endyr_D_",abb_i)
      if(endyr_nm_i%in%names(init_p)){init_p[[endyr_nm_i]] <- tail(names(xi),1)}

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

    ## discard cvs
    #  This is done in a separate loop because the discard cvs don't always match discard time series supplied to the dat file
    #  (e.g. Black Sea Bass SEDAR 56)
    init_obs_cv_D_nm <- names(init_b)[grepl("^obs_cv_D",names(init_b))]
    for(nm_i in init_obs_cv_D_nm){
      obs_cv_D_b_i <- init_b[[nm_i]]
      # If the cv was provided in the last year of the base model, project it forward
      if(paste(endyr)%in%names(obs_cv_D_b_i)){
        obs_cv_D_p_i <- setNames(rep(bamExtras::geomean2(tail(obs_cv_D_b_i,nyb_rcn$L)),nyp),yp)
        obs_cv_D_i <- c(obs_cv_D_b_i,obs_cv_D_p_i)
      }else{
        obs_cv_D_i <- obs_cv_D_b_i
      }
      init_p[[nm_i]] <- obs_cv_D_i
    }

    ## cpue
    init_obs_cpue_nm <- local({
      nm <- names(init_b)[grepl("^obs_cpue",names(init_b))]
      nm_abb <- gsub("^(obs_cpue)(_)(.*)","\\3",nm)
      nm_yes <- nm[which(nm_abb%in%dimnames(U)[[2]])]
      nm_no <- nm[which(!nm_abb%in%dimnames(U)[[2]])]
      if(length(nm_no)>0){
        message(paste(paste(nm_no,collapse=", "), "were found in the bam tpl but not in the rdat. They may not be included in the likelihood."))
      }
      # Only include names found in U (indices reported in the rdat)
      return(nm_yes)
    })

    for(nm_i in init_obs_cpue_nm){
      abb_i <- gsub("^(obs_cpue)(_)(.*)","\\3",nm_i)
      xbi <- setNames(as.numeric(init_b[[nm_i]]),names(init_b[[nm_i]])) # get index values from init
      # yrsbi <- names(xbi)
      xproji <- setNames(U[,abb_i],dimnames(U)[[1]])
      ndigi <- round(median(nchar(gsub("([0-9]*.)([0-9]*)","\\2",xbi))))
      xpi <- round(xproji[paste(yp)],ndigi) # get projected index values
      xi <- c(xbi,xpi[!is.na(xpi)])
      yrsi <- names(xi)

      xcvbi <- init_b[[paste0("obs_cv_cpue_",abb_i)]] # get index cv values from init
      xprojcvi <- setNames(U_ts_cv[,paste0("cv.U.",abb_i)],dimnames(U_ts_cv)[[1]])
      ndigcvi <- round(median(nchar(gsub("([0-9]*.)([0-9]*)","\\2",xcvbi))))
      xcvpi <- round(xprojcvi[paste(yp)],ndigcvi) # get projected index values
      xcvi <- c(xcvbi,xcvpi[!is.na(xcvpi)])

      init_p[[nm_i]] <- setNames(paste(xi),names(xi))
      # init_p[[paste0("obs_cv_cpue_",abb_i)]] <- setNames(paste(round(U_ts_cv[,paste0("cv.U.",abb_i)],ndigcvi)),rownames(U_ts_cv))[names(xi)]
      init_p[[paste0("obs_cv_cpue_",abb_i)]] <- setNames(paste(xcvi),names(xcvi))

      # update styr, endyr, yrs, and nyr as necessary
      yrinfoi <- setNames(list(min(yrsi),max(yrsi),yrsi,paste(length(yrsi))),
                          paste0(c("styr","endyr","yrs","nyr"),gsub("obs","",nm_i))
      )
      yrinfoi_nm_is <- names(yrinfoi)[names(yrinfoi)%in%names(init_p)]
      init_p[yrinfoi_nm_is] <- yrinfoi[yrinfoi_nm_is]
    }

    ## agec
    init_obs_agec_nm <- names(init_b)[grepl("^obs_agec",names(init_b))]
    init_obs_agec_nm_key <- setNames(names(acomp),
                                     paste0("obs_agec_",gsub("^(D.)([A-Za-z]+)","\\2_D",gsub("^L.","",names(acomp)))))
    # If the age comp names aren't using the _D suffix, attempt to remove _D suffix from the names in the key
    if(!any(grepl("\\_D$",init_obs_agec_nm))){
      names(init_obs_agec_nm_key) <- gsub("\\_D$","",names(init_obs_agec_nm_key))
    }

    for(nm_i in init_obs_agec_nm){
      nm2_i <- init_obs_agec_nm_key[[nm_i]]
      xbi <- init_b[[nm_i]]
      x2bi <- acomp[[nm2_i]]
      ndigi <- round(median(nchar(gsub("([0-9]*.)([0-9]*)","\\2",xbi))))
      obsi <- apply(x2bi,2,function(x){sprintf(paste0("%.",ndigi,"f"), x)})
      attributes(obsi) <- attributes(x2bi)
      yrsi <- rownames(obsi)
      init_p[[nm_i]] <- obsi

      # update styr, endyr, yrs, and nyr as necessary
      yrinfoi <- setNames(list(min(yrsi),max(yrsi),yrsi,paste(length(yrsi))),
                          paste0(c("styr","endyr","yrs","nyr"),gsub("obs","",nm_i))
      )
      yrinfoi_nm_is <- names(yrinfoi)[names(yrinfoi)%in%names(init_p)]
      init_p[yrinfoi_nm_is] <- yrinfoi[yrinfoi_nm_is]

      # update nfish and nsamp as necessary
      nfish_nm_i <- gsub("obs","nfish",nm_i)
      nsamp_nm_i <- gsub("obs","nsamp",nm_i)
      nfishbi <- init_b[[nfish_nm_i]]
      nsampbi <- init_b[[nsamp_nm_i]]
      nfishpi <- setNames(rep(bamExtras::geomean2(nfishbi[paste(tail(yb,nyb_rcn$comp))]),nyp),paste(yp))
      nsamppi <- setNames(rep(bamExtras::geomean2(nsampbi[paste(tail(yb,nyb_rcn$comp))]),nyp),paste(yp))
      nfishdigi <- round(median(nchar(gsub("([0-9]*.)([0-9]*)","\\2",nfishbi))))
      nsampdigi <- round(median(nchar(gsub("([0-9]*.)([0-9]*)","\\2",nsampbi))))
      nfishi <- setNames(sprintf(paste0("%.",nfishdigi,"f"),round(as.numeric(c(nfishbi,nfishpi)))),names(c(nfishbi,nfishpi)))
      nsampi <- setNames(sprintf(paste0("%.",nsampdigi,"f"),round(as.numeric(c(nsampbi,nsamppi)))),names(c(nsampbi,nsamppi)))

      init_p[[nfish_nm_i]] <- nfishi[yrsi]
      init_p[[nsamp_nm_i]] <- nsampi[yrsi]
    }

    ## lenc
    init_obs_lenc_nm <- names(init_b)[grepl("^obs_lenc",names(init_b))]
    init_obs_lenc_nm_key <- setNames(names(lcomp),
                                     paste0("obs_lenc_",gsub("^(D.)([A-Za-z]+)","\\2_D",gsub("^L.","",names(lcomp)))))
    # If the length comp names aren't using the _D suffix, attempt to remove _D suffix from the names in the key
    if(!any(grepl("\\_D$",init_obs_lenc_nm))){
      names(init_obs_lenc_nm_key) <- gsub("\\_D$","",names(init_obs_lenc_nm_key))
    }

    for(nm_i in init_obs_lenc_nm){
      nm2_i <- init_obs_lenc_nm_key[[nm_i]]
      xbi <- init_b[[nm_i]]
      x2bi <- lcomp[[nm2_i]]
      ndigi <- round(median(nchar(gsub("([0-9]*.)([0-9]*)","\\2",xbi))))
      obsi <- apply(x2bi,2,function(x){sprintf(paste0("%.",ndigi,"f"), x)})
      attributes(obsi) <- attributes(x2bi)
      yrsi <- rownames(obsi)
      init_p[[nm_i]] <- obsi

      # update styr, endyr, yrs, and nyr as necessary
      yrinfoi <- setNames(list(min(yrsi),max(yrsi),yrsi,paste(length(yrsi))),
                          paste0(c("styr","endyr","yrs","nyr"),gsub("obs","",nm_i))
      )
      yrinfoi_nm_is <- names(yrinfoi)[names(yrinfoi)%in%names(init_p)]
      init_p[yrinfoi_nm_is] <- yrinfoi[yrinfoi_nm_is]

      # update nfish and nsamp as necessary
      nfish_nm_i <- gsub("obs","nfish",nm_i)
      nsamp_nm_i <- gsub("obs","nsamp",nm_i)
      nfishbi <- init_b[[nfish_nm_i]]
      nsampbi <- init_b[[nsamp_nm_i]]
      nfishpi <- setNames(rep(bamExtras::geomean2(nfishbi[paste(tail(yb,nyb_rcn$comp))]),nyp),paste(yp))
      nsamppi <- setNames(rep(bamExtras::geomean2(nsampbi[paste(tail(yb,nyb_rcn$comp))]),nyp),paste(yp))
      nfishdigi <- round(median(nchar(gsub("([0-9]*.)([0-9]*)","\\2",nfishbi))))
      nsampdigi <- round(median(nchar(gsub("([0-9]*.)([0-9]*)","\\2",nsampbi))))
      nfishi <- setNames(sprintf(paste0("%.",nfishdigi,"f"),round(as.numeric(c(nfishbi,nfishpi)))),names(c(nfishbi,nfishpi)))
      nsampi <- setNames(sprintf(paste0("%.",nsampdigi,"f"),round(as.numeric(c(nsampbi,nsamppi)))),names(c(nsampbi,nsamppi)))

      init_p[[nfish_nm_i]] <- nfishi[yrsi]
      init_p[[nsamp_nm_i]] <- nsampi[yrsi]
    }

    bam_p <- bam2r(dat_obj=bam$dat,tpl_obj=bam$tpl,cxx_obj=bam$cxx,init=init_p)
  }else{
    bam_p <- NULL
  } # end if(!is.null(bam2r_args))

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
  # Cn.L
  p <- ggplot(Cn.L,mapping=aes(x=year,y=Ln))+
    geom_area(aes(fill=fleet))+
    theme_bw()+
    scale_fill_brewer(palette="Spectral")+
    stat_summary(fun = sum, geom = "line", size = 1)+
    stat_summary(fun = sum, geom = "point", size = 2)+
    geom_vline(xintercept = endyr, linetype="dashed", size = 0.3)
  p2 <- p + annotate("text",x=endyr, label="endyr\n",y=Inf, angle=90,hjust=1.5,vjust=1)
  print(p2)


  # Cn.D
  if(nrow(Cn.D)>0){
  p2 <- p %+% Cn.D + aes(y=Dn) +
  annotate("text",x=endyr, label="endyr\n",y=Inf, angle=90,hjust=1.5,vjust=1)
  print(p2)
}

  # Cw.L
  p2 <- p %+% Cw.L + aes(y=Lw) +
  annotate("text",x=endyr, label="endyr\n",y=Inf, angle=90,hjust=1.5,vjust=1)
  print(p2)


  # Cw.D
  if(nrow(Cw.D)>0){
  p2 <- p %+% Cw.D + aes(y=Dw) +
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
