#' Plot results from MCBE uncertainty analysis
#'
#' @param sim_summary Output object from summarize_MCBE
#' @param dir_figs name of directory that will be created to store figures
#' @param fileName Name of BAM base files in dir_bam_base, not including file extensions.
#' @param dir_bam_base Name of directory where the bam base model is located relative to the working directory, for adding base model results to plots (optional).
#' @param filter_info list of ranges used to filter out MCBE results. set to NULL to retain all results.
#' @param tseries_plot_args List of arguments to pass to \code{bamExtras::plot_boot_vec} when plotting time series.
#' @param aseries_plot_args List of arguments to pass to \code{bamExtras::plot_boot_vec} when plotting age series.
#' @param base_tseries_args List of arguments to add to plotting functions (e.g. \code{points()}) that add base results to time series plots
#' @param base_aseries_args List of arguments to add to plotting functions (e.g. \code{points()}) that add base results to age series plots
#' @param nm_rp list specifying names of reference points to include in plots.
#' @param plot2pdf If \code{plot2pdf} = TRUE, each plot will be written to a pdf
#' in \code{dir_figs}. If false, figures will be plotted to the active graphic device.
#' @param sc_pdf Optional scalars for pdf width and height. Multiplied by defaults pdf width and height. Defaults to 1 which has no effect.
#' @keywords bam MCBE stock assessment fisheries
#' @export
#' @examples
#' \dontrun{
#' # Run MCBE, writing files to dir_bam_sim
#' run_MCBE("GrayTriggerfish", dir_bam_sim="sim_GrTr")
#' # Summarize results of MCBE and assign to object
#' MCBE_GrTr <- summarize_MCBE(dir_bam_sim="sim_GrTr")
#' # Plot MCBE results
#' }



plot_MCBE <- function(sim_summary = NULL,
                      dir_figs = "figs",
                      fileName     = "bam",
                      dir_bam_base = NULL,
                      filter_MCBE_args = list("filter_info"=
                                                list("gradient.max"  = c(0,0.01),
                                                     "Fref"          = c(0,5),
                                                     "F.Fref"        = c(0,5),
                                                     "R.sigma.par"   = c(0.2,1.0),
                                                     "R0_prob"       = c(0.005,0.995), # Not limits but probabilities used to compute limits
                                                     "avoid_bounds"  = TRUE # logical. If TRUE, filter out runs with parameters near the bounds
                                                ),
                                              "nm_rp" = list("parms"=c("Fref"="Fmsy","F.Fref"="Fend.Fmsy.mean"))
                                              ),
                      tseries_plot_args = list(),
                      aseries_plot_args = list(),
                      base_tseries_args = list("type"="o","lty"=1,"lwd"=2,"pch"=16, "col"="black"),
                      base_aseries_args = list("type"="o","lty"=1,"lwd"=2,"pch"=16, "col"="black"),
                      nm_rp = list("parms"=c("Fref"="Fmsy","F.Fref"="Fend.Fmsy.mean",
                                             "Sref"="msst","S.Sref"="SSBend.MSST",
                                             "Lref"="msy.klb","Bref"="Bmsy"),
                                    "t.series"=c("F"="F.full","F.Fref"="F.Fmsy","S"="SSB","S.Sref"="SSB.msst")
                      ),
                      plot2pdf=TRUE,
                      sc_pdf=list("w"=1,"h"=1)
                      ){
if(!dir.exists(dir_figs)){
  dir.create(dir_figs,recursive = TRUE)
}

# filter_MCBE_args <- modifyList(filter_MCBE_args,list("nm_rp"=nm_rp))

ss <- sim_summary

parms <- ss$parms
if("spr.brps"%in%names(ss)){
  spr.brps <- ss$spr.brps
  nm_spr.brps <- names(spr.brps)
  nm_parms <- names(parms)
  nm_spr.brps.new <- nm_spr.brps[!nm_spr.brps%in%nm_parms]
  parms <- cbind(parms,spr.brps[,nm_spr.brps.new])
  message(paste("Added the following variables from spr.brps to bootstrapped parms:",paste(nm_spr.brps.new,collapse=", ")))
}else{
  message("spr.brps not found in sim_summary")
}
nm_parms <- names(parms)

# parm.cons <- ss$parm.cons
a.series <- ss$a.series

styr <- unique(parms$styr)
if(length(styr)>1){
  styr <- min(styr)
  warning(paste0("styr not the same among all sims. Will attempt to use earliest styr (",styr,")"))
  }

endyr <- unique(parms$endyr)
if(length(endyr)>1){
  endyr <- max(endyr)
  warning(paste0("endyr not the same among all sims. Will attempt to use latest endyr (",endyr,")"))
  }
yrs <- paste(styr:endyr)
t.series <- lapply(ss$t.series,function(x){x[,yrs]})

is_base <- !is.null(dir_bam_base)

env_rp <- new.env()
env_rp$nm_rp2 <- c() # names of reference points assigned to vectors (e.g. Fref) rather than matrices of time series (e.g. F)
for(i in seq_along(nm_rp)){
  nm_i <- names(nm_rp)[i]
  nm_rp_i <- nm_rp[[nm_i]]
  xi <- ss[[nm_i]]
  # Don't try to assign values from anything other than a matrix or data frame.
  if(is.matrix(xi)|is.data.frame(xi)){
    for(j in seq_along(nm_rp_i)){
      nm_rp_ij <- nm_rp_i[j]
      if(nm_rp_ij%in%names(xi)){
          xij <- xi[,nm_rp_ij]}else{
          message(paste0(nm_rp_ij," not found in ", nm_i,". Set a valid value for ", names(nm_rp_ij) ," in nm_rp."))
          xij <- NA
        }
      if(names(nm_rp_ij)%in%ls(env_rp)){
        warning(paste(names(nm_rp_ij),
                      "was specified twice in nm_rp. Only the first value will be used.",nm_rp_ij,
                      "will be ignored."))}else{
                        assign(names(nm_rp_ij),xij,envir=env_rp)
                        env_rp$nm_rp2 <- c(env_rp$nm_rp2,nm_rp_ij)

        }
    }
  }
}

### Do stuff with base run
if(is_base){
rb <- rdat_base <- dget(file.path(dir_bam_base,paste(fileName,'rdat', sep="."))) # spp

styr_b <- rb$parms$styr
endyr_b <- rb$parms$endyr
yrs_b <- paste(styr_b:endyr_b)
t.series.b <- rb$t.series[yrs_b,]

a.series.b <- rb$a.series

parms.b <- rb$parms

if("spr.brps"%in%names(rb)){
  spr.brps.b <- rb$spr.brps

  nm_spr.brps.b <- names(spr.brps.b)
  nm_parms.b <- names(parms.b)
  nm_spr.brps.b.new <- nm_spr.brps.b[!nm_spr.brps.b%in%nm_parms.b]
  parms.b <- c(parms.b,spr.brps.b[nm_spr.brps.b.new])
  message(paste("Added the following variables from spr.brps to base parms:",paste(nm_spr.brps.b.new,collapse=", ")))

}else{
  message("spr.brps not found in base model rdat.")
}
nm_parms.b <- names(parms.b)
}


# nm_sim <- sprintf(paste("%0",nchar(nsim),".0f",sep=""),1:nsim)
#
#
# ## Get base model stuff
#
# # Identify working directory
#   wd <- getwd()
#   message(paste("working directory:",wd))
#

  nsim <- length(parms$modRunName)
  simID2 <- simID <- 1:nsim

  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ####   plotsConvergence   ####
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Lref_sim <- Fref_sim <- Sref_sim <- rep(NA,nsim)
if(nsim>1){
  nm_Lref <- env_rp$nm_rp2[["Lref"]]
  if(!nm_Lref%in%nm_parms){
    message(paste(nm_Lref,"not found in parms. Will not appear in convergence plots."))
  }else{
    Lref_sim <- parms[,nm_Lref]
  }
  nm_Fref <- env_rp$nm_rp2[["Fref"]]
  if(!nm_Fref%in%nm_parms){
    message(paste(nm_Fref,"not found in parms. Will not appear in convergence plots."))
  }else{
    Fref_sim <- parms[,nm_Fref]
  }
  nm_Sref <- env_rp$nm_rp2[["Sref"]]
  if(!nm_Sref%in%nm_parms){
    message(paste(nm_Sref,"not found in parms. Will not appear in convergence plots."))
  }else{
    Sref_sim <- parms[,nm_Sref]
  }

    if(plot2pdf){pdf(file.path(dir_figs,"Fig-MCBE-converge.pdf"),width=sc_pdf$w*7,height=sc_pdf$h*7)}
      sd.Lref <- sd.Fref <- sd.Sref <- rep(NA,nsim)
      for(i in 1:nsim)
      {
        #nm_msy <- names(parms)[grepl("^msy.(mt|klb)",names(parms))]
        # msy_sim <- parms[[nm_msy]]
        sd.Lref[i]=sd(Lref_sim[1:i])
        sd.Fref[i]=sd(Fref_sim[1:i])
        sd.Sref[i]=sd(Sref_sim[1:i])
      }
      par(mfrow=c(3,1),mar=c(2.5,3,0,0),oma=c(0,0,1,1.5),mgp=c(1.2,.3,0),tck=-.02)

      plot(1:nsim, sd.Lref,xlab="", ylab=paste0("SE(",nm_Lref,")"), type="l", lwd=2)
      plot(1:nsim, sd.Fref,xlab="", ylab=paste0("SE(",nm_Fref,")"), type="l", lwd=2)
      plot(1:nsim, sd.Sref,xlab="Number bootstrap replicates", ylab=paste0("SE(",nm_Sref,")"), type="l", lwd=2)
        if(plot2pdf){dev.off()}

      if(plot2pdf){pdf(file.path(dir_figs,"lk.total.v.gradient.max.pdf"),width=sc_pdf$w*7,height=sc_pdf$h*7)}
      par(mfrow=c(3,1),mar=c(2.5,3,0,0),oma=c(0,0,1,1.5),mgp=c(1.2,.3,0),tck=-.02)
    x <- ss$like$lk.total
    y <- ss$like$gradient.max
    plot(x=x,
         y=y,
         pch=16,col=rgb(0,0,0,0.2),
         xlab="lk.total",ylab="gradient.max",
         xlim=range(x)+c(-1,1),ylim=range(y)+diff(range(y))*c(0,.1))
    # points(x=spp$like[["lk.total"]],
    #        y=spp$like[["gradient.max"]],col="red",pch=3)
    legend("topright",legend=paste("n=",length(x)),bty="n")

    plot_boot_density(ss$like,"lk.total",plot_args = list(xlab = "lk.total", ylab="density",main = "", type = "l"))
    # abline(v=spp$like[["lk.total"]],col="red",lwd=2)
    plot_boot_density(ss$like,"gradient.max",plot_args = list(xlab = "gradient.max", ylab="density",main = "", type = "l"))

    # abline(v=spp$like[["gradient.max"]],col="red",lwd=2)
        if(plot2pdf){dev.off()}
}else{
  message("nsim is not greater than 1. convergence plots will not be drawn")
}
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ####       trimRuns       ####
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    filter_MCBE_res <- do.call(filter_MCBE,c(list(sim_summary = ss), filter_MCBE_args))
    sim_pass <- filter_MCBE_res$sim_pass
    simID2 <- simID[sim_pass]

    if(plot2pdf){pdf(file.path(dir_figs,"gradient.max.density.pdf"),width=sc_pdf$w*7,height=sc_pdf$h*7)}
    par(mfrow=c(1,1))
    if(length(simID2)>1){
    plot_boot_density(data=ss$like[simID2,],x_name="gradient.max")
    }else{
      message("Filtered results do not contain more than 1 sim run. gradient.max.density not plotted.")
    }
        if(plot2pdf){dev.off()}

    # Filter some objects
    parms.2 <- parms[simID2,]
    # parm.cons.2 <-lapply(parm.cons,function(xi){xi[simID2,]})
    a.series.2 <- lapply(a.series,function(xi){xi[simID2,]})
    t.series.2 <- lapply(t.series,function(xi){xi[simID2,]})

  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ####     plotMCBData      ####
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Plot input data
    # Landings
    if(plot2pdf){pdf(file.path(dir_figs,"Fig-MCBE-landings.ob.pdf"),width=sc_pdf$w*7,height=sc_pdf$h*7)}

    names_L_ob <- names(t.series)[grepl("^L.*ob$",names(t.series))]

    par(mfrow=c(length(names_L_ob),1),mar=c(2,5,0.5,0.5),mgp=c(1.5,.3,0),cex.axis=1.5,cex.lab=1.5,tck=0.02)
    for(i in names_L_ob){
      obj.i <- t.series.2[[i]]
      if(is_base){tseries_plot_args$ylim <- range(c(obj.i,t.series.b[,i]),na.rm=TRUE)}
      do.call(plot_boot_vec,c(list(data=obj.i,xlab="",ylab=i),tseries_plot_args))
      if(is_base){
        do.call("points",c(list("x"=t.series.b$year,"y"=t.series.b[,i]),base_tseries_args))
      }
    }
        if(plot2pdf){dev.off()}

    # Discards
    if(plot2pdf){pdf(file.path(dir_figs,"Fig-MCBE-discards.ob.pdf"),width=sc_pdf$w*7,height=sc_pdf$h*7)}

    names_D_ob <- names(t.series)[grepl("^D.*ob$",names(t.series))]
    if(length(names_D_ob)>0){
    par(mfrow=c(length(names_D_ob),1),mar=c(2,5,0.5,0.5),mgp=c(1.5,.3,0),cex.axis=1.5,cex.lab=1.5,tck=0.02)
    for(i in names_D_ob){
      obj.i <- t.series.2[[i]]
      if(is_base){tseries_plot_args$ylim <- range(c(obj.i,t.series.b[,i]),na.rm=TRUE)}
      do.call(plot_boot_vec,c(list(data=obj.i,xlab="",ylab=i),tseries_plot_args))
      if(is_base){
        do.call("points",c(list("x"=t.series.b$year,"y"=t.series.b[,i]),base_tseries_args))
      }
    }
    }else{
      warning("Discard time series not found. No plot produced.")
    }
        if(plot2pdf){dev.off()}

    # Indices of abundance
    if(plot2pdf){pdf(file.path(dir_figs,"Fig-MCBE-indices.ob.pdf"),width=sc_pdf$w*7,height=sc_pdf$h*7)}

    names_U_ob <- names(t.series)[grepl("^U.*ob$",names(t.series))]

    par(mfrow=c(length(names_U_ob),1),mar=c(2,5,0.5,0.5),mgp=c(1.5,.3,0),cex.axis=1.5,cex.lab=1.5,tck=0.02)
    for(i in names_U_ob){
      obj.i <- t.series.2[[i]]
      if(any(which(obj.i<0))){
        message(paste("Negative values found in",i,"will not be plotted."))
      }
      obj.i[obj.i<0] <- NA
      if(is_base){
        t.series.b[,i]
        xi <- t.series.b$year
        yi <- t.series.b[,i]
        yi[yi<0] <- NA
        tseries_plot_args$ylim <- range(c(obj.i,yi),na.rm=TRUE)}
      do.call(plot_boot_vec,c(list(data=obj.i,xlab="",ylab=i),tseries_plot_args))
      if(is_base){
        do.call("points",c(list("x"=xi,"y"=yi),base_tseries_args))
      }
    }
        if(plot2pdf){dev.off()}

    # Natural mortality
    if(plot2pdf){pdf(file.path(dir_figs,"Fig-MCBE-naturalMortality.ob.pdf"),width=sc_pdf$w*7,height=sc_pdf$h*7)}
    par(mfrow=c(2,1),mar=c(3,5,0.5,0.5),mgp=c(1.5,.3,0),cex.axis=1.5,cex.lab=1.5,tck=0.02)
    if("M.constant"%in%names(parms.2)&nsim>1){
    plot_boot_density(data=parms.2,x_name="M.constant")
     abline(v=median(parms.2[,"M.constant"]),lty=2,lwd=2)
     if(is_base){
       abline(v=parms.b$M.constant,lty=1,lwd=2)
     }
    }else{
      warning("M.constant not found in parms or nsim is not greater than 1. M.constant density plot not produced.")
    }

    if("M"%in%names(a.series.2)){
      if(is_base){aseries_plot_args$ylim <- range(c(a.series.2$M,a.series.b[,"M"]),na.rm=TRUE)}
      do.call(plot_boot_vec,c(list(data=a.series.2$M,xlab="age",ylab="natural mortality"),aseries_plot_args))
    # plot_boot_vec(data=a.series.2$M,
    #                #ref_x=spp$a.series$age,ref_y=spp$a.series$M,
    #                xlab="age",ylab="natural mortality")
      if(is_base){
        do.call("points",c(list("x"=a.series.b$age,"y"=a.series.b[,"M"]),base_aseries_args))
      }
    }else{
      warning("M not found in a.series. No plot produced.")
    }
    if(plot2pdf){dev.off()}

  # Discard mortality
    if(plot2pdf){pdf(file.path(dir_figs,"Fig-MCBE-discardMortality.ob.pdf"),width=sc_pdf$w*7,height=sc_pdf$h*7)}
    names_Dmort <- names(parms.2)[grepl("^D.mort",names(parms.2))]
    if(length(names_Dmort)>0&nsim>1){
      par(mfrow=c(length(names_Dmort),1),mar=c(3,5,0.5,0.5),mgp=c(1.5,.3,0),cex.axis=1.5,cex.lab=1.5,tck=0.02)
      for(i in names_Dmort){

        plot_boot_density(data=parms.2,x_name=i)
        abline(v=median(parms.2[,i]),lty=2,lwd=2)
        if(is_base){
          abline(v=parms.b[[i]],lty=1,lwd=2)
        }
      }
    }
         if(plot2pdf){dev.off()}

  # #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # ####     plotTSeries      ####
    nm_ts2 <- names(t.series.2)

    par(mfrow=c(1,1),mar=c(2.5,2.5,0.5,0.5),mgp=c(1,.25,0),cex.axis=1,cex.lab=1,tck=0.02)

    nm_par <- nm_rp$t.series[["F"]]
    if(!nm_par%in%nm_ts2){
      message(paste(nm_par,"not found in sim_summary$t.series. Plot not drawn."))
    }else{
    if(plot2pdf){pdf(file.path(dir_figs,"Fig-MCBE-apicalF.pdf"),width=sc_pdf$w*7,height=sc_pdf$h*7)}
    if(is_base){tseries_plot_args$ylim <- range(c(t.series.2[[nm_par]],t.series.b[,nm_par]),na.rm=TRUE)}
    do.call(plot_boot_vec,c(list(data=t.series.2[[nm_par]],xlab="",ylab=nm_par),tseries_plot_args))
    if(is_base){
      do.call("points",c(list("x"=t.series.b$year,"y"=t.series.b[,nm_par]),base_tseries_args))
    }
    if(plot2pdf){dev.off()}
    }

    nm_par <- nm_rp$t.series[["F.Fref"]]
    if(!nm_par%in%nm_ts2){
      message(paste(nm_par,"not found in sim_summary$t.series. Plot not drawn."))
    }else{
    if(plot2pdf){pdf(file.path(dir_figs,"Fig-MCBE-FdFmsy.pdf"),width=sc_pdf$w*7,height=sc_pdf$h*7)}
    par(mfrow=c(1,1),mar=c(2.5,2.5,0.5,0.5),mgp=c(1,.25,0),cex.axis=1,cex.lab=1,tck=0.02)
    if(is_base){tseries_plot_args$ylim <- range(c(t.series.2[[nm_par]],t.series.b[,nm_par]),na.rm=TRUE)}
    do.call(plot_boot_vec,c(list(data=t.series.2[[nm_par]],xlab="",ylab=nm_par),tseries_plot_args))
    # abline(h=1,lwd=2)
    if(is_base){
      do.call("points",c(list("x"=t.series.b$year,"y"=t.series.b[,nm_par]),base_tseries_args))
    }
    if(plot2pdf){dev.off()}
    }

    nm_par <- "N"
    if(!nm_par%in%nm_ts2){
      message(paste(nm_par,"not found in sim_summary$t.series. Plot not drawn."))
    }else{
    if(plot2pdf){pdf(file.path(dir_figs,"Fig-MCBE-N.pdf"),width=sc_pdf$w*7,height=sc_pdf$h*7)}
      par(mfrow=c(1,1),mar=c(2.5,2.5,0.5,0.5),mgp=c(1,.25,0),cex.axis=1,cex.lab=1,tck=0.02)
    if(is_base){tseries_plot_args$ylim <- range(c(t.series.2[[nm_par]],t.series.b[,nm_par]),na.rm=TRUE)}
    do.call(plot_boot_vec,c(list(data=t.series.2[[nm_par]],xlab="",ylab=nm_par),tseries_plot_args))
    if(is_base){
      do.call("points",c(list("x"=t.series.b$year,"y"=t.series.b[,nm_par]),base_tseries_args))
    }
    if(plot2pdf){dev.off()}
    }

    nm_par <- nm_rp$t.series[["S"]]
    if(!nm_par%in%nm_ts2){
      message(paste(nm_par,"not found in sim_summary$t.series. Plot not drawn."))
    }else{
    if(plot2pdf){pdf(file.path(dir_figs,"Fig-MCBE-SSB.pdf"),width=sc_pdf$w*7,height=sc_pdf$h*7)}
    par(mfrow=c(1,1),mar=c(2.5,2.5,0.5,0.5),mgp=c(1,.25,0),cex.axis=1,cex.lab=1,tck=0.02)
    if(is_base){tseries_plot_args$ylim <- range(c(t.series.2[[nm_par]],t.series.b[,nm_par]),na.rm=TRUE)}
    do.call(plot_boot_vec,c(list(data=t.series.2[[nm_par]],xlab="",ylab=nm_par),tseries_plot_args))
    if(is_base){
      do.call("points",c(list("x"=t.series.b$year,"y"=t.series.b[,nm_par]),base_tseries_args))
    }
    if(plot2pdf){dev.off()}
    }

    nm_par <- nm_rp$t.series[["S.Sref"]]
    if(!nm_par%in%nm_ts2){
      message(paste(nm_par,"not found in sim_summary$t.series. Plot not drawn."))
    }else{
    if(plot2pdf){pdf(file.path(dir_figs,"Fig-MCBE-SSBdMSST.pdf"),width=sc_pdf$w*7,height=sc_pdf$h*7)}
    par(mfrow=c(1,1),mar=c(2.5,2.5,0.5,0.5),mgp=c(1,.25,0),cex.axis=1,cex.lab=1,tck=0.02)
    if(is_base){tseries_plot_args$ylim <- range(c(t.series.2[[nm_par]],t.series.b[,nm_par]),na.rm=TRUE)}
    do.call(plot_boot_vec,c(list(data=t.series.2[[nm_par]],xlab="",ylab=nm_par),tseries_plot_args))
    # abline(h=1,lwd=2)
    if(is_base){
      do.call("points",c(list("x"=t.series.b$year,"y"=t.series.b[,nm_par]),base_tseries_args))
    }
     # abline(h=1,lwd=2)
     if(plot2pdf){dev.off()}
    }

    nm_par <- "SSB.SSBmsy"
    if(!nm_par%in%nm_ts2){
      message(paste(nm_par,"not found in sim_summary$t.series. Plot not drawn."))
    }else{
     if(plot2pdf){pdf(file.path(dir_figs,"Fig-MCBE-SSBdSSBmsy.pdf"),width=sc_pdf$w*7,height=sc_pdf$h*7)}
     par(mfrow=c(1,1),mar=c(2.5,2.5,0.5,0.5),mgp=c(1,.25,0),cex.axis=1,cex.lab=1,tck=0.02)
     if(is_base){tseries_plot_args$ylim <- range(c(t.series.2[[nm_par]],t.series.b[,nm_par]),na.rm=TRUE)}
     if(nm_par%in%names(t.series.2)){
       do.call(plot_boot_vec,c(list(data=t.series.2[[nm_par]],xlab="",ylab=nm_par),tseries_plot_args))
       # abline(h=1,lwd=2)
       if(is_base){
         do.call("points",c(list("x"=t.series.b$year,"y"=t.series.b[,nm_par]),base_tseries_args))
       }
       # abline(h=1,lwd=2)
     }
     if(plot2pdf){dev.off()}
    }

    nm_par <- "recruits"
    if(!nm_par%in%nm_ts2){
      message(paste(nm_par,"not found in sim_summary$t.series. Plot not drawn."))
    }else{
     if(plot2pdf){pdf(file.path(dir_figs,"Fig-MCBE-recruits.pdf"),width=sc_pdf$w*7,height=sc_pdf$h*7)}
     par(mfrow=c(1,1),mar=c(2.5,2.5,0.5,0.5),mgp=c(1,.25,0),cex.axis=1,cex.lab=1,tck=0.02)
     if(is_base){tseries_plot_args$ylim <- range(c(t.series.2[[nm_par]],t.series.b[,nm_par]),na.rm=TRUE)}
     do.call(plot_boot_vec,c(list(data=t.series.2[[nm_par]],xlab="",ylab=nm_par),tseries_plot_args))
     if(is_base){
       do.call("points",c(list("x"=t.series.b$year,"y"=t.series.b[,nm_par]),base_tseries_args))
     }
     if(plot2pdf){dev.off()}
    }

    # Benchmarks time series, two panel plot
     if(plot2pdf){pdf(file.path(dir_figs,"Fig-MCBE-status-ts.pdf"),width=sc_pdf$w*7,height=sc_pdf$h*7)}
    par(mfrow=c(2,1),mar=c(3,3,1,1),mgp=c(1,.2,0),cex.lab=1,cex.axis=1,cex=1,tck=-0.02)
    nm_par <- nm_rp$t.series[["S.Sref"]]
    if(!nm_par%in%nm_ts2){
      message(paste(nm_par,"not found in sim_summary$t.series. Plot not drawn."))
      plot.new()
      legend("center",legend=paste("no plot for",nm_par),bty="n")
    }else{
    if(is_base){tseries_plot_args$ylim <- range(c(t.series.2[[nm_par]],t.series.b[,nm_par]),na.rm=TRUE)}
    do.call(plot_boot_vec,c(list(data=t.series.2[[nm_par]],xlab="",ylab=nm_par),tseries_plot_args))
    # abline(h=1,lwd=2)
    if(is_base){
      do.call("points",c(list("x"=t.series.b$year,"y"=t.series.b[,nm_par]),base_tseries_args))
    }
    # abline(h=1,lwd=2)
    }

    nm_par <- nm_rp$t.series[["F.Fref"]]
    if(!nm_par%in%nm_ts2){
      message(paste(nm_par,"not found in sim_summary$t.series. Plot not drawn."))
      plot.new()
      legend("center",legend=paste("no plot for",nm_par),bty="n")
    }else{
    if(is_base){tseries_plot_args$ylim <- range(c(t.series.2[[nm_par]],t.series.b[,nm_par]),na.rm=TRUE)}
    do.call(plot_boot_vec,c(list(data=t.series.2[[nm_par]],xlab="",ylab=nm_par),tseries_plot_args))
    # abline(h=1,lwd=2)
    if(is_base){
      do.call("points",c(list("x"=t.series.b$year,"y"=t.series.b[,nm_par]),base_tseries_args))
    }
    # abline(h=1,lwd=2)
    }
    if(plot2pdf){dev.off()}

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ####     plotASeries      ####
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(plot2pdf){pdf(file.path(dir_figs,"Fig-MCBE-length-at-age.pdf"),width=sc_pdf$w*7,height=sc_pdf$h*7)}
    nm_par <- "length"
    if(nm_par%in%names(a.series.2)){
      par(mfrow=c(1,1),mar=c(3,2,1,1),mgp=c(1.2,.3,0),cex.lab=1,cex.axis=1,cex=1,tck=-0.02)
      if(is_base){aseries_plot_args$ylim <- range(c(a.series.2[[nm_par]],a.series.b[,nm_par]),na.rm=TRUE)}
      do.call(plot_boot_vec,c(list(data=a.series.2[[nm_par]],xlab="age",ylab=nm_par),aseries_plot_args))
      if(is_base){
        do.call("points",c(list("x"=a.series.b$age,"y"=a.series.b[,nm_par]),base_aseries_args))
      }
    }else{
      warning(paste(nm_par,"not found in a.series. No plot produced."))
    }
    if(plot2pdf){dev.off()}

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ####     plotDensity      ####
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ##### status ####
    nmfile <- "Fig-MCBE-status.pdf"
    if(length(simID2)>1){
      if(plot2pdf){pdf(file.path(dir_figs,nmfile),width=sc_pdf$w*7,height=sc_pdf$h*7)}
      par(mfrow=c(2,1),mar=c(3,3,1,1),mgp=c(1,.2,0),cex.lab=1,cex.axis=1,cex=1,tck=-0.02)
      nm_par <- env_rp$nm_rp2[["S.Sref"]]
      if(!nm_par%in%nm_parms){
        message(paste(nm_par,"not found in parms. Plot not drawn."))
        plot.new()
        legend("center",legend=paste("no plot for",nm_par),bty="n")
      }else{
        dt <- parms[simID2,]
        plot_boot_density(data=dt,x_name=nm_par,quantile_draw = TRUE)
        if(is_base){abline(v=parms.b[[nm_par]],lty=1,lwd=2)}
      }
      nm_par <- env_rp$nm_rp2[["F.Fref"]]
      if(!nm_par%in%nm_parms){
        message(paste(nm_par,"not found in parms. Plot not drawn."))
        plot.new()
        legend("center",legend=paste("no plot for",nm_par),bty="n")
      }else{
        dt <- parms[simID2,]
        plot_boot_density(data=dt,x_name=nm_par,quantile_draw = TRUE)
        if(is_base){abline(v=parms.b[[nm_par]],lty=1,lwd=2)}
      }
      if(plot2pdf){dev.off()}
    }else{
      message(paste("Filtered results do not contain more than 1 sim run.",nmfile,"not produced."))
    }

    ##### benchmarks ####
    nmfile <- "Fig-MCBE-benchmarks.pdf"
    if(length(simID2)>1){
      if(plot2pdf){pdf(file.path(dir_figs,nmfile),width=sc_pdf$w*7,height=sc_pdf$h*7)}
      par(mfrow=c(2,2),mar=c(3,2,1,1),mgp=c(1,.2,0),cex.lab=1,cex.axis=1,cex=1,tck=-0.02)
      nm_par <- env_rp$nm_rp2[["Fref"]]
      if(!nm_par%in%nm_parms){
        message(paste(nm_par,"not found in parms. Plot not drawn."))
        plot.new()
        legend("center",legend=paste("no plot for",nm_par),bty="n")
      }else{
        dt <- parms[simID2,]
        plot_boot_density(data=dt,x_name=nm_par,quantile_draw = TRUE)
        if(is_base){abline(v=parms.b[[nm_par]],lty=1,lwd=2)}
      }
      nm_par <- env_rp$nm_rp2[["Sref"]]
      if(!nm_par%in%nm_parms){
        message(paste(nm_par,"not found in parms. Plot not drawn."))
        plot.new()
        legend("center",legend=paste("no plot for",nm_par),bty="n")
      }else{
        dt <- parms[simID2,]
        plot_boot_density(data=dt,x_name=nm_par,quantile_draw = TRUE)
        if(is_base){abline(v=parms.b[[nm_par]],lty=1,lwd=2)}
      }
      nm_par <- env_rp$nm_rp2[["Lref"]]
      if(!nm_par%in%nm_parms){
        message(paste(nm_par,"not found in parms. Plot not drawn."))
        plot.new()
        legend("center",legend=paste("no plot for",nm_par),bty="n")
      }else{
        dt <- parms[simID2,]
        plot_boot_density(data=dt,x_name=nm_par,quantile_draw = TRUE)
        if(is_base){abline(v=parms.b[[nm_par]],lty=1,lwd=2)}
      }
      nm_par <- env_rp$nm_rp2[["Bref"]]
      if(!nm_par%in%nm_parms){
        message(paste(nm_par,"not found in parms. Plot not drawn."))
        plot.new()
        legend("center",legend=paste("no plot for",nm_par),bty="n")
      }else{
        dt <- parms[simID2,]
        plot_boot_density(data=dt,x_name=nm_par,quantile_draw = TRUE)
        if(is_base){abline(v=parms.b[[nm_par]],lty=1,lwd=2)}
      }
      if(plot2pdf){dev.off()}
    }else{
      message(paste("Filtered results do not contain more than 1 sim run.",nmfile,"not produced."))
    }

    ##### stock-recruit parameters ####
    nmfile <- "Fig-MCBE-SR-parms.pdf"
    if(length(simID2)>1){
      if(plot2pdf){pdf(file.path(dir_figs,nmfile),width=sc_pdf$w*7,height=sc_pdf$h*7)}
    par(mfrow=c(2,2),mar=c(3,2,1,1),mgp=c(1.2,.3,0),tck=-0.01,lend="butt")

    nm_par <- c("BH.R0","R0")
    if(!any(nm_par%in%nm_parms)){
      message(paste(paste(nm_par,collapse=", "),"not found in parms. Plot not drawn."))
      plot.new()
      legend("center",legend=paste("no plot for",nm_par),bty="n")
    }else{
      dt <- parms[simID2,]
      nm_par1 <- nm_par[which(nm_par%in%nm_parms)][1] # Just use the first match in case both versions were added to parms
      plot_boot_density(data=dt,x_name=nm_par1,quantile_draw = TRUE)
      if(is_base){abline(v=parms.b[[nm_par1]],lty=1,lwd=2)}
    }

    nm_par <- c("BH.steep","steep")
    if(!any(nm_par%in%nm_parms)){
      message(paste(paste(nm_par,collapse=", "),"not found in parms. Plot not drawn."))
      plot.new()
      legend("center",legend=paste("no plot for",nm_par),bty="n")
    }else{
      dt <- parms[simID2,]
      nm_par1 <- nm_par[which(nm_par%in%nm_parms)][1] # Just use the first match in case both versions were added to parms
      plot_boot_density(data=dt,x_name=nm_par1,quantile_draw = TRUE)
      if(is_base){abline(v=parms.b[[nm_par1]],lty=1,lwd=2)}
    }
    nm_par <- c("BH.Phi0","Phi0")
    if(!any(nm_par%in%nm_parms)){
      message(paste(paste(nm_par,collapse=", "),"not found in parms. Plot not drawn."))
      plot.new()
      legend("center",legend=paste("no plot for",nm_par),bty="n")
    }else{
      dt <- parms[simID2,]
      nm_par1 <- nm_par[which(nm_par%in%nm_parms)][1] # Just use the first match in case both versions were added to parms
      plot_boot_density(data=dt,x_name=nm_par1,quantile_draw = TRUE)
      if(is_base){abline(v=parms.b[[nm_par1]],lty=1,lwd=2)}
    }
    nm_par <- c("BH.R.sigma.par","R.sigma.par")
    if(!any(nm_par%in%nm_parms)){
      message(paste(paste(nm_par,collapse=", "),"not found in parms. Plot not drawn."))
      plot.new()
      legend("center",legend=paste("no plot for",nm_par),bty="n")
    }else{
      dt <- parms[simID2,]
      nm_par1 <- nm_par[which(nm_par%in%nm_parms)][1] # Just use the first match in case both versions were added to parms
      plot_boot_density(data=dt,x_name=nm_par1,quantile_draw = TRUE)
      if(is_base){abline(v=parms.b[[nm_par1]],lty=1,lwd=2)}
    }
    if(plot2pdf){dev.off()}
    }else{
      message(paste("Filtered results do not contain more than 1 sim run.",nmfile,"not produced."))
    }

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ####      plotPhase      ####
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(plot2pdf){pdf(file.path(dir_figs,"Fig-MCBE-phase.pdf"),width=sc_pdf$w*7,height=sc_pdf$h*7)}
    par(mfrow=c(1,1),mar=c(2.5,2.5,0.5,0.5),mgp=c(1,.25,0),cex.axis=1,cex.lab=1,tck=0.02)
    nm_par_x <- env_rp$nm_rp2[["F.Fref"]]
    nm_par_y <- env_rp$nm_rp2[["S.Sref"]]
    test_tmp <- !c(nm_par_x,nm_par_y)%in%nm_parms
    if(any(test_tmp)){
      paste(c(nm_par_x,nm_par_y)[which(test_tmp)],collapse=" and ")
      message(paste(paste(c(nm_par_x,nm_par_y)[which(test_tmp)],collapse=" and "),"not found in parms. Plot not drawn."))
      plot.new()
      legend("center",legend="no plot",bty="n")
    }else{
      dt <- parms[simID2,]
      if(is_base){
        xref <- parms.b[[nm_par_x]]
        yref <- parms.b[[nm_par_y]]
        }else{
          xref <- NULL
          yref <- NULL
        }
      plot_boot_phase(x=dt[,nm_par_x],y=dt[,nm_par_y],xref=xref,yref=yref,xlab=nm_par_x,ylab=nm_par_y)
    }
    if(plot2pdf){dev.off()}

  }
