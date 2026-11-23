#' Plot age or length composition data
#'
#' @param data_ob observed composition data. numeric matrix
#' @param data_pr predicted composition data with the same dimensions as data_ob. numeric matrix
#' @param n1 sample size by year associated with comps type 1 (e.g. number of trips)
#' @param n2 sample size by year associated with comps type 2 (e.g. number of fish)
#' @param cc catch curve statistics output from bamExtras::catch_curve() function
#' @param data_type Used in part of y-axis label. Specify type of data: "catch at age", "number at age", "proportion at age". character vector
#' @param year_type Used in part of y-axis label. Specify type of data with respect to year: e.g. "year caught", "cohort year". character vector
#' @param ylab will overwrite the y-axis label pasted together from data_type and year_type
#' @param vlines values for adding vertical lines to plots. numeric vector
#' @param signif_cc Number of significant digits to display when plotting catch curve statistics.
#' @param ... pass to \code{plot}

#' @keywords bam stock assessment fisheries
#' @author Nikolai Klibansky
#' @export
#' @examples
#' \dontrun{
#' # read in commercial handline age comps used in the rex.ob <- rdat_BlackSeaBass$comp.mats$acomp.cHL.ob
#' x.pr <- rdat_BlackSeaBass$comp.mats$acomp.cHL.pr
#' n <- rdat_BlackSeaBass$t.series[rownames(x),c("acomp.cHL.n","acomp.cHL.nfish")]
#' # Plot comps by year
#' comp_plot(data_ob=x.ob,data_pr = x.pr)
#' # Add sample sizes
#' comp_plot(data_ob=x.ob,data_pr = x.pr, n1=n[,1],n2=n[,2])
#' # Plot comps and add catch curves (note the log-transformation of the composition data)
#' comp_plot(log(x),cc=catch_curve(x),fillComp = FALSE,ylab= "log(proportion)",xlab="age",title="black sea bass commercial handline catch curves")
#' }
#'

comp_plot <- function(data_ob,data_pr=NULL,n1=NULL,n2=NULL,cc=NULL,data_type="proportion",year_type="year",
                      title="",xlab="",ylab=NULL,ylim=NULL,
                      addGrid=FALSE,byrow=FALSE,
                      cex_n=1.3,
                      cex_year=1.3,
                      fillComp=TRUE, colFill=rgb(0,0,0,0.25), vlines=NULL,
                      signif_cc=3,
                      ...){
  if(!is.matrix(data_ob)){data_ob <- as.matrix(data_ob)
  message("data_ob coerced to matrix")
  }

  # Calculate dimensions of plot matrix
  nr <- ceiling(sqrt(nrow(data_ob))) # Determine number of rows in figure
  nc <- ceiling(nrow(data_ob)/nr)    # Determine number of columns in figure
  n.panels <- nr*nc # Calculate the number of panels in figure
  n.empty <- n.panels-nrow(data_ob) # Calculate number of empty panels in figure

  # Plot comps by year
  if(byrow){
    par(mfrow=c(nr,nc),mar=c(0,0,0,0),oma=c(3,3,2,0),mgp=c(0.8,.3,0),tck=-0.01)
    xaxt <- c(rep("n",n.panels-nc),rep("s",nc)) # Set up vector so that x axes plot only on the bottom row
    yaxt <- rep(c("s",rep("n",nc-1)),nr) # Set up vector so that y axes plot only on the leftmost column
  }else{
    par(mfcol=c(nr,nc),mar=c(0,0,0,0),oma=c(3,3,2,0),mgp=c(0.8,.3,0),tck=-0.01)
    xaxt <- rep(c(rep("n",nr-1),"s"),nc) # Set up vector so that x axes plot only on the bottom row
    yaxt <- c(rep("s",nr),rep("n",n.panels-nr)) # Set up vector so that y axes plot only on the leftmost column
  }
  mode(data_ob) <- "numeric"
  x <- as.numeric(colnames(data_ob))

  if(is.null(ylim)){
    ylim <- local({ # This stupid thing makes it so that -Inf is not included in ylim when y is log transformed
      V <- as.vector(as.matrix((data_ob[,paste(x)])))
      return(range(V[V>-Inf],na.rm=TRUE))
    })

  }

  for(i in 1:n.panels){
    if(i<=n.empty){ # Plot empty panels
      plot(x,rep(NA,length(x)),
           xaxt=xaxt[i],yaxt=yaxt[i],
           ylim=ylim,...)
    }else{
      r <- i-n.empty
      y <- data_ob[r,paste(x)]
      plot(x,y,
           pch=16, type="p",
           xaxt=xaxt[i],yaxt=yaxt[i],
           ylim=ylim,...)
      if(fillComp){
        polygon(x=c(x,rev(x)),y=c(data_ob[r,paste(x)],rep(par("usr")[3],length(x))),
                border=NA,col=colFill)
      }
      if(!is.null(data_pr)){ # add predicted values from comp matrix of same dimensions
        mode(data_pr) <- "numeric"
        points(x,data_pr[r,paste(x)],type="l",lwd=2)
      }

      # Add fitted catch curve
      if(!is.null(cc)){
        with(cc[r,],{
          if(!any(is.na(c(age_min,age_max)))){ # Basically, if the catch curve was fit to the data, plot it.
            x_pr <<- age_min:age_max
            y_pr <<- intercept + slope*x_pr

            points(x_pr,y_pr,type="l",lwd=2,col="red")

            # Print regression statistics on plot
            legend("right",legend=c(paste("Z =",signif(-slope,signif_cc)),paste("P =",(signif(slope_pval,signif_cc)))),
                   bty="n",text.col="blue")
          }
        })
      }
      if(addGrid){
        grid()
      }

      # Add year to plot
      legend("topright",legend=rownames(data_ob)[r],
             cex=cex_year,bty="n",text.col="blue")

      # Add sample size to plot
      # n1 and n2
      if(!is.null(n1)){
        if(is.null(n2)){
          leg_text <- paste("n =",n1[r])
        }else{
          leg_text <- c(paste("n1 =",n1[r]),paste("n2 =",n2[r]))
        }

        if(any(!is.na(c(n1[r],n2[r])))){ # If the value of n1 is not NA
          legend("right",legend=leg_text,
                 cex=cex_n,bty="n",text.col="blue")
        }
      }
      if(!is.null(vlines)&any(!is.na(y))){
        abline(v=vlines,lty=2)
      }

    }
  }
  if(is.null(ylab)){
    ylab <- paste(data_type," by ",year_type,sep="")}

  mtext(side=2,line=1.5,text=ylab,outer=TRUE) # y-axis label
  mtext(side=1,line=1.5,text=xlab,outer=TRUE) # x-axis label
  mtext(side=3,line=0.2,text=title,cex=1.5,outer=TRUE) # title of figure
}
