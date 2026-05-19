#' Fit catch curves to age composition data
#'
#' Simple geometric mean function using only is.finite values of log(x).
#' @param M matrix of proportions at age where rows are years (catch years or cohort years) and columns are ages. Rows should sum to 1.
#' @param age_min minimum age to use in catch curves. If a value is not supplied, see methods used in age_min_method
#' @param age_max maximum age to use in catch curves. If the value is not supplied, the default is to use the second oldest age class, since the
#' oldest age class is often a plus group.
#' @param nobs_min the minimum number of usable observations (i.e. catch > 0) required before a catch curve should be calculated
#' @param age_min_method character. allows you to set the method that should be used to calculate the minimum age used in the catch curves.
#                 peak_overall sets the min age at the peak in the age comp, pooling all years while peak_by_year sets age min similarly
#                 but for each year separately
#' @keywords bam stock assessment fisheries
#' @author Nikolai Klibansky
#' @export
#' @examples
#' \dontrun{
#' # read in black ses bass age comps used in the recent stock assessment
#' rdat <- rdat_BlackSeaBass
#' comp <- rdat$comp.mats # List of composition matrices
#' acomp <- comp[grepl("^acomp.*.ob$",names(comp))] # list of observed age composition matrices
#'
#' comp_i <- acomp$acomp.sTV.ob
#' cc_i <- catch_curve(comp_i)
#'
#' # Plot comps and add catch curves (note the log-tranformation of the composition data)
#' comp_plot(log(comp_i),cc=cc_i,fillComp = FALSE,ylab= "log(proportion)",xlab="age",title="black sea bass chevron trap catch curves")
#'
#' # Plot Z estimates over time for multiple sets of age compositions
#' par(mfrow=c(1,1),mgp=c(1,.2,0),mar=c(2,2,1,1),oma=c(0,0,0,0),tck=0.01)
#' M <- rdat$parms$M.msst
#' cc <- lapply(acomp,catch_curve)
#' abb <- gsub("^(acomp.)|(.ob)$","",names(cc)) # abbreviations of age composition fleets
#' cols <- rainbow(length(cc))
#' xlim <- range(unlist(lapply(cc,function(x){as.numeric(x$year)})),na.rm=TRUE)
#' ylim <- range(c(0,unlist(lapply(cc,function(x){-as.numeric(x$slope)}))),na.rm=TRUE)
#' for(i in 1:length(cc)){
#'   xi <- cc[[i]]
#'   year <- xi$year
#'   Z <- -xi$slope
#'   if(i==1){
#'     plot(year,Z,xlim=xlim,ylim=ylim,type="o",col=cols[i])
#'   }else{
#'     points(year,Z,type="o",col=cols[i])
#'   }
#' }
#' # Add reference line for natural mortality
#' abline(h=M,lty=2,lwd=2)
#' text(par("usr")[1]+0.05*diff(par("usr")[1:2]),M,labels="M",pos=3)
#'
#' legend("bottomleft",legend=abb,col=cols,pch=1,lty=1,bty="n")
#'
#' }



catch_curve <- function(M,age_min=NULL,age_max=NULL,nobs_min=3,age_min_method="peak_overall"){
  # Initialize output
  Ve <- rep(NA,nrow(M))
  out <- data.frame("year"=rownames(M),"age_min" = Ve, "age_max" = Ve,
                    "intercept"=Ve,"slope"=Ve,"slope_se"=Ve,"slope_pval"=Ve,"r_sq"=Ve)
  if(age_min_method=="peak_overall"&is.null(age_min)){
    x <- as.numeric(colnames(M))
    y <- as.numeric(colSums(M,na.rm=TRUE))
    age_min <- x[which(y==max(y))[1]]
  }

  for(i in 1:nrow(M)){
    usable_elements <- M[i,which(M[i,]>0)] # Determine which elements from row i or matrix M are usable
    if(length(usable_elements)>=nobs_min){ # Determine if there are enough observations to use
      x <- as.numeric(names(usable_elements))
      y <- as.numeric(usable_elements)

      if(is.null(age_min)){
        if(age_min_method=="peak_by_year"){
          x_min <- x[which(y==max(y))[1]]
        }else{warning("Must specify either age_min or valid age_min_method (peak_overall or peak_by_year)")}
      }else{
        x_min <- age_min
      }

      if(is.null(age_max)){
        x_max <- rev(sort(unique(x)))[2] # Set to second oldest age class to avoid potential plus group
      }else{
        x_max <- age_max
      }
      x_sub <- x[which(x%in%(x_min:x_max))] # Subset x based on x range
      y_sub <- y[which(x%in%x_sub)]       # Subset y to match x

      if(length(x_sub)>=nobs_min){  # Again, determine if there are enough observations to use
        lm1 <- lm(log(y_sub)~x_sub)

        out[i,"age_min"] <- x_min
        out[i,"age_max"] <- x_max
        out[i,"intercept"] <- coef(lm1)[[1]]
        out[i,"slope"]     <- coef(lm1)[[2]]
        out[i,"slope_se"]  <- summary(lm1)$coefficients[2,"Std. Error"]
        out[i,"slope_pval"]  <- summary(lm1)$coefficients[2,"Pr(>|t|)"]
        out[i,"r_sq"] <- summary(lm1)$r.squared
      }
    }
  }
  return(out)
}
