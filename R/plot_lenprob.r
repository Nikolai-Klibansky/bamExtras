#' Plot distribution of size at age (age-length key, lenprob)
#'
#' Plot distribution of size at age (age-length key, lenprob) as an image
#' @param spp bam output rdat object (list). This function is expecting the rdat to contain a list of lenprob matrices \code{lenprob.mats} which will be plotted separately.
#' @param par_args list of arguments to pass to \code{graphics::par}
#' @param plot_args list of arguments to pass to \code{plot.matrix::plot}. This may include standard arguments for \code{base::plot} but there are also special arguments that can be passed to \code{plot.matrix::plot}. See help files for \code{plot.matrix::plot} for details.
#' @keywords bam stock assessment fisheries
#' @author Nikolai Klibansky
#' @note Requires that plot.matrix package be installed to plot properly
#' @export
#' @examples
#' \dontrun{
#' # Most rdat stored in bamExtras contains \code{lenprob.mats} and many have multiple matrices within the list.
#' plot_lenprob(rdat_BlackSeaBass)
#' plot_lenprob(rdat_GagGrouper)
#' plot_lenprob(rdat_GrayTriggerfish)
#' plot_lenprob(rdat_VermilionSnapper)
#' # Change the color palette
#' plot_lenprob(rdat_GrayTriggerfish,plot_args=list(col=heat.colors))
#' # I think it's kind of nice to have the color for values near zero to contrast strongly with the rest of the color palette.
#' plot_lenprob(rdat_GrayTriggerfish,plot_args=list(col=c(rgb(1,1,1),heat.colors(30))))
#' plot_lenprob(rdat_GrayTriggerfish,plot_args=list(col=c(rgb(1,1,1),heat.colors(100))))
#' }
#'
plot_lenprob <- function(spp,
                         par_args=list(),
                         plot_args=list()){
  require(plot.matrix)

  par_args_user <- par_args
  par_args_default <- list(mfrow=c(1,1),mar=c(2,2,1.5,3),oma=c(1,1,1,1),mgp=c(1,.2,0),tck=-0.01)
  par_args <- modifyList(par_args_default,par_args_user)

  plot_args_user <- plot_args

  if("lenprob.mats"%in%names(spp)){
    lpm <- spp$lenprob.mats
    for(i in seq_along(lpm)){
      nm_lpm_i <- names(lpm)[i]
      #pdf(file.path(df1,"misc",paste0(nm_lpm_i,".pdf")),width = 8,height=6)
      lpm_i <- lpm[[nm_lpm_i]]
      lpm_i <- as.matrix(t(lpm_i[,ncol(lpm_i):1])) # flip and transpose to plot correctly

      do.call(par,par_args)

      n_col <- 20
      col_lp <- c(rgb(1,1,1),gray(((n_col*.8):1)/n_col))
      plot_args_default <- list(x=lpm_i,
                                       fmt.cell ="%.2f",
                                       fmt.key ="%.2f",
                                       digits=2,
                                       col=col_lp,
                                       xlab="age bins",
                                       ylab="length bins",
                                       main=nm_lpm_i)
      plot_args <- modifyList(plot_args_default, plot_args_user)
      do.call(plot,plot_args)
      # plot(x=as.matrix(lpm_i),fmt.cell ="%.2f",fmt.key ="%.2f",digits=2,col=col_lp,
      #      xlab="age bins",ylab="length bins",main=nm_lpm_i)
      #dev.off()
    }

  }else{
    message("lenprob.mats not found in spp")
  }
}
