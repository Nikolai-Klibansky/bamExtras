#' Time series plot
#'
#' This is a simple plotting fucntion that is useful when plotting time series typically found in stock assessments (e.g. catch, indices)
#' @param x data frame of mean observations by year, where columns are different groups (e.g. fleets)
#' @param err error values for each value in x. Must be same dimensions as x
#' @param errType type of error provided: "cv"= coefficient of variation or "sd"=standard deviation. If "cv" then x is multiplied by err to compute sd.
#' @param errStyle Style of displaying error around the central tendency time series x. Can be"bars" or "bands"
#' @param errUnits Units of err. Used only for displaying text on plot (e.g. "SE")
#' @param nErrUnits Multiplied by err to determine how far above and below x to to plot the errors
#' @param errLim Range for limiting the plotting of error bars or bands.
#' @param plotLabel character string for labeling plot
#' @param matplot_args optional list of additional arguments passed to \code{matplot}, which plots x.
#' @param polygon_args optional list of additional arguments passed to \code{polygon}, which plots error bands when errStyle="bands".
#' @param col_band_alpha alpha (transparency) value between 0 and 1 to adjust opacity of color bands
#' @keywords bam stock assessment fisheries
#' @author Nikolai Klibansky
#' @export
#' @examples
#' \dontrun{
#' rdat <- rdat_BlackSeaBass
#' tser <- rdat$t.series
#' U <- tser[,grepl("^U.*ob$",names(tser))]
#' Ucv <- tser[,grepl("^cv.U.",names(tser))]
#' tseries_plot(U,Ucv,plotLabel=rdat$info$species,errType = "cv"),matplot_args=list(xlab="year",ylab="index of abundance"))
#' }
#'

tseries_plot <- function(x,
                         err,
                         errType="sd",
                         errStyle="bands",
                         errUnits="SE",
                         nErrUnits=2,
                         errLim = c(0,Inf),
                         plotLabel=NULL,
                         matplot_args = list(),
                         polygon_args = list(),
                         col_band_alpha=0.3,
                         pt_buffer=0.1){

  matplot_args_user <- matplot_args
  matplot_args_default <- list(col = rainbow(ncol(x)),type = "o",lty = 1,lwd = 2,pch = 16,ylim=NULL,xlab="x",ylab="y")
  matplot_args <- modifyList(matplot_args_default,matplot_args_user)

  polygon_args_user <- polygon_args
  polygon_args_default = list(border=NA,lwd=polygon_args$lwd)
  polygon_args <- modifyList(polygon_args_default,polygon_args_user)

  if(errType=="sd"){
    err <- err
  }else if(errType=="cv"){
    err <- x*err
  }else{
    errUnits <- errType
    message("errType not a recognized type (sd or cv). errUnits set equal to errType.")
  }

  xLo <- apply(x-nErrUnits*err,2,function(x){pmin(pmax(x,errLim[1]),errLim[2])})
  xUp <- apply(x+nErrUnits*err,2,function(x){pmin(pmax(x,errLim[1]),errLim[2])})

  nSer <- ncol(x)
  col_band <- color_trans(matplot_args$col,col_band_alpha)

  if(is.null(matplot_args$ylim)){
    matplot_args$ylim <- c(range(c(0,x,xUp),na.rm=TRUE))
  }


  # plot
  par(mfrow=c(1,1),mar=c(3,3,1,3),mgp=c(1.5,.3,0),oma=c(0,0,0,0),tck=0.01,cex.lab=1.5,cex.axis=1.5)
  legend_text <- dimnames(x)[[2]]

  do.call(matplot,c(list(x=as.numeric(dimnames(x)[[1]]),y=x),matplot_args))

  for(i in 1:nSer){
    xi <- as.numeric(dimnames(x)[[1]])
    yilo <- xLo[,i]
    yiup <- xUp[,i]

    cc <- complete.cases(xi,yilo,yiup)

    xi <- xi[cc]
    xisect <- sections(xi,1) # Section numbers for non-consecutive years
    yilo <- yilo[cc]
    yiup <- yiup[cc]

    if(errStyle=="bars"){
      arrows(x0=xi,
             y0=yilo,
             y1=yiup,
             length=0,col=col[i]
      )
    }
    if(errStyle=="bands"){
      for(j in unique(xisect)){
      jw <- which(xisect==j)
      xij <- xi[jw]
      yiloj <- yilo[jw]
      yiupj <- yiup[jw]

      # Add two more points to create a thin polygon to show error around isolated points.
      if(length(jw)==1){
      xij <- xij+pt_buffer*c(-1,1)
      yiloj <- rep(yiloj,2)
      yiupj <- rep(yiupj,2)
      }

      do.call(polygon,
              c(list(
                x=c(xij,rev(xij)),
                y=c(yiloj,rev(yiupj)),
                col=col_band[i]
              ),
              polygon_args)
      )
      }
    }
  }

  do.call(legend,c(list("topright",legend=legend_text,
         bty="n"), matplot_args[c("col","lty","lwd","pch")]))

  if(!is.null(plotLabel)){
    text_legend(plotLabel)
  }

  legend("right",bty="n",legend=paste("Error",errStyle,"\nplotted to",nErrUnits,errUnits))

  invisible(list("x"=x,"err"=err,"xLo"=xLo,"xUp"=xUp))
}
