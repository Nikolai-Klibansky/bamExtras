#' Plot standard deviation of normalized residuals (sdnr)
#'
#' Plot standard deviation of normalized residuals (sdnr) values found in bam output (rdat$sdnr) as barplot
#' @param spp bam output rdat object
#' @param par_args list of arguments to pass to \code{graphics::par}
#' @keywords bam stock assessment fisheries
#' @author Nikolai Klibansky
#' @export
#' @examples
#' \dontrun{
#' plot_sdnr(rdat_BlackSeaBass)
#' plot_sdnr(rdat_AtlanticMenhaden)
#' }
#'

plot_sdnr <- function(spp,par_args=list()){
  par_args_user <- par_args
  par_args_default <- list(mar=c(3,4,2,1),mgp=c(1.5,0.3,0),tck=-0.01)
  par_args <- modifyList(par_args_default,par_args_user)
  if("sdnr"%in%names(spp)){
  sdnr <- spp$sdnr
  sdnr <- sdnr[order(names(sdnr))]

  do.call(par,par_args)
  bp <- barplot(sdnr,names.arg=gsub("sdnr.","",names(sdnr)),horiz=TRUE,
                las=1, xlim=c(0,max(sdnr)*1.2),
                xlab="sdnr",main="sdnr")
  abline(v=seq(0,ceiling(max(sdnr)),by=0.5),lty=3,col="gray")
  text(x=sdnr,y=bp,labels = round(sdnr,3),pos=4,col="blue")
  abline(v=1,lty=2,lwd=2)
  }else{
    message("sdnr not found in spp")
  }
}

