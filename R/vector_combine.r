#' combine vectors using Conn 2010 hierarchical method
#'
#' @description This is a version of a function that was written by Paul Conn, Kyle Shertzer, and Erik Williams, using jags to combine indices of abundance. Nikolai Klibansky converted that script into this function.
#' @param x observed vectors as a data frame
#' @param cv observed cvs associated with x, as a data frame
#' @param key_nm If a common column, like "year", is provided in x and cv as a key to relate values between data frames
#' provide the name as key_nm. Should be coercible to a numeric vector. If NULL, then row names will be used.
#' @param jags_args list of arguments to be passed to R2jags::jags() within the function
#' @param cv.na value to to replace missing cv values with
#' @param run_mcmcplot Should mcmcplots::mcmcplot be run to plot jags results?
#' @param mcmcplot_args if run_mcmcplot = TRUE, additional arguements can be passed to mcmcplot
#' @references Conn, P. B. 2010. Hierarchical analysis of multiple noisy abundance indices. Canadian Journal of Fisheries and Aquatic Sciences 67:108-120.
#' @keywords bam stock assessment fisheries Conn method combine indices
#' @author Nikolai Klibansky, Kyle Shertzer, Erik Williams, Paul Conn
#' @details
#' Additional details...
#'   This function was designed to be used to combine multiple indices of abundance into a single index, but it could also be used to combine other types of numeric vectors.
#'   Requires JAGS install and mcmcplots package.
#'
#' @export
#' @examples
#' \dontrun{
#' # Example usage

#' library(ggplot2)
#' library(tidyr)
#'
#' # Combining commercial indices of abundance from the recent blueline tilefish assessment for the US Southeast Atlantic
#' ts <- rdat_BluelineTilefish$t.series
#'
#' x <- ts[,grepl("^U.c.*.ob$",names(ts))]
#' cv <- ts[,grepl("^cv.U.c",names(ts))]
#'
#' names(x) <- gsub(".ob$","",names(x))
#' names(cv) <- gsub(".ob$","",names(cv))
#'
#' # Only include years where at least one index had real values
#' yrs <- rownames(x)[apply(x,1,function(y){any(!is.na(y))})]
#' x <- x[yrs,]
#' cv <- cv[yrs,]
#'
#' # run the function
#' vc_out <- vector_combine(x=x,cv=cv)
#'
#' # plot results with mcmcplot (use dir argument if you want to save the plots somewhere)
#' mcmcplots::mcmcplot(vc_out$mcmcout)
#'
#' # plot original indices and resulting index (with ggplot, to be hip)
#' dfout <- vc_out$dfout
#' x2 <- cbind(x,"combined"=dfout$mu)
#' sd <- x*cv
#' sd2 <- cbind(cv,"combined"=dfout$sd)
#' names(sd2) <- gsub("^cv.","",names(sd2))
#' x2lo <- x2-2*sd2
#' names(x2lo) <- paste0(names(x2lo),".lo")
#' x2up <- x2+2*sd2
#' names(x2up) <- paste0(names(x2up),".up")
#'
#' x2 <- cbind("year"=rownames(x2),x2)
#' sd2 <- cbind("year"=rownames(sd2),sd2)
#'
#' x3 <- tidyr::pivot_longer(x2,cols=!year,names_to="fleet",values_to="index")
#' sd3 <- tidyr::pivot_longer(sd2,cols=!year,names_to="fleet",values_to="sd")
#'
#' x4 <- na.omit(merge(x3,sd3,by=c("year","fleet")))
#'
#' p <- ggplot(x4, aes(x=year, y=index, group=fleet, color=fleet,linetype=fleet,size=fleet)) +
#'   geom_line() +
#'   geom_point()+
#'   geom_errorbar(aes(ymin=index-2*sd, ymax=index+2*sd), width=1, position=position_dodge(0.2)) +
#'   labs(title="indices of abundance", x="year", y = "index of abundance") +
#'   theme_classic() +
#'   scale_color_manual(values=c("black","red","blue")) +
#'   scale_linetype_manual(values=c("solid","dashed","dashed")) +
#'   scale_size_manual(values=c(1,0.5,0.5))
#'
#' p
#' }

vector_combine <- function(x,
                           cv,
                           key_nm = NULL,
                           jags_args = list(
                             parameters.to.save = c("chi","mu","sigma"),
                             n.chains = 3,
                             n.thin = 1,
                             n.burnin = 10000,
                             n.iter = 80000
                           ),
                           cv.na = 0.05,
                           run_mcmcplot = FALSE,
                           mcmcplot_args = list(dir = "mcmcplot")
){
  library(R2jags)
  library(mcmcplots)
  library(R2WinBUGS)

  model.file <- function(){
    for(xi in 1:nx){
      for(ti in 1:nt){
        lnx[xi,ti]~dnorm(temp.mean[xi,ti],temp.tau[xi,ti])
        temp.tau[xi,ti]<-1/(sigma[xi]*sigma[xi]+log(cvsq[xi,ti]+1))  #exact form given by reviewer
      }
      chi[xi]~dnorm(a,2)
    }
    for(xi in 1:nx){
      for(ti in 1:nt){
        temp.mean[xi,ti]<-log(mu[ti])+chi[xi]
      }
      sigma[xi]~dunif(0,5)
    }
    b<-log(100)
    for(ti in 1:nt) {
      nu[ti]~dnorm(b,1)
      mu[ti]<-exp(nu[ti])
    }
  }


  inits <- function(){
    list(chi=rnorm(nx,log(.01),1),nu=rep(0,nt),sigma=rep(0.5,nx))
  }

  if(!is.null(key_nm)){
    xx <- as.numeric(x[,key_nm])
    rownames(x) <- x[,key_nm]
    x <- x[,names(x)!=key_nm]
  }else{
    xx <- as.numeric(dimnames(x)[[1]])
  }

  tx <- t(x)
  nx <- ncol(x)
  nt <- nrow(x)
  lnx <- t(log(x))
  cvsq <- (t(cv))^2
  cvsq[which(is.na(cvsq))] <- cv.na

  data <- list("lnx" = lnx, "nx" = nx, "nt" = nt,"cvsq" = cvsq, "a"=log(0.01))


  #### Run JAGS ####
  # This usually runs pretty quickly
  mcmcout <- do.call(R2jags::jags,c(jags_args,
                           list(data = data,
                                inits = inits,
                                model.file = model.file
                                )))


  #### Run mcmcplot ####
  # This can take a long time to run (>20 minutes)
  if(run_mcmcplot){
  st <- Sys.time()
  if("dir" %in% names(mcmcplot_args)){
    dir <- mcmcplot_args$dir
    if(!dir.exists(dir)){
      dir.create(dir)
    }
  }

  do.call(mcmcplots::mcmcplot,c(list(mcmcout=mcmcout),mcmcplot_args))
  print(Sys.time()-st)
}
  #### Output stuff ###

  rn <- rownames(mcmcout$BUGSoutput$summary)
  Mu_est=mcmcout$BUGSoutput$summary[rn[grepl("mu",rn)],"mean"]
  mean_Mu=mean(Mu_est)

  # standardize index to a mean of 1
  Mu_est=Mu_est/mean_Mu

  cn <- colnames(mcmcout$BUGSoutput$sims.matrix)
  all.iters <- mcmcout$BUGSoutput$sims.matrix[,cn[grepl("mu",cn)]]
  all.std <- all.iters/rowMeans(all.iters)
  sd.std <- apply(all.std,2,sd)
  cv <- sd.std/Mu_est
  dfout <- data.frame(id=xx, mu=Mu_est, sd=sd.std, cv=cv)

  invisible(list(mcmcout=mcmcout,dfout=dfout))
}
