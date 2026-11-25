#' Plot likelihood component values
#'
#' Plot likelihood component values values found in bam output (spp$like) as barplot
#' @param spp bam output rdat object
#' @param par_args list of arguments to pass to \code{graphics::par}
#' @param digits number of significant digits to round labels to
#' @keywords bam stock assessment fisheries
#' @author Nikolai Klibansky
#' @export
#' @examples
#' plot_lk(rdat_BlackSeaBass)
#' plot_lk(rdat_RedSnapper)
#' plot_lk(rdat_Cobia)
#' plot_lk(rdat_AtlanticMenhaden)
#'

plot_lk <- function(spp,par_args=list(),digits=6){
  par_args_user <- par_args
  par_args_default <- list(mar=c(3,4,2,1),mgp=c(1.5,0.3,0),tck=-0.01)
  par_args <- modifyList(par_args_default,par_args_user)
  if("like"%in%names(spp)){

      like <- spp$like
      like.type <- unlist(lapply(strsplit(names(like),".",fixed=T),"[[",1))
      lk <- like[which(like.type=="lk")]
      lk.parts <- lk[!names(lk)%in%c("lk.total","lk.unwgt.data")]
      lk.parts <- lk.parts[rev(names(lk.parts))]
      lk.total <- lk["lk.total"]
      lk.unwgt.data <- lk["lk.unwgt.data"]

      do.call(par,par_args)
      bp <- barplot(lk.parts,names.arg=gsub("lk.","",names(lk.parts)),horiz=TRUE,
                    las=1,xlim=c(0,lk.total*1.1),
                    xlab="negative log likelihood",main="likelihood components")

      # Add grid
      abline(v=pretty(lk),lty=3,col="gray")

      # Add lk labels
      text(x=lk.parts,y=bp,labels = signif(lk.parts,digits),pos=4, srt=0,offset=0.2,adj=0,col="blue")

      # Add line at zero
      abline(v=0)

      cols <- c("blue","red")
      ltys <- c(2,3)
      lwds <- c(2,2)

      # Add lk.total reference line
      abline(v=lk.total,lty=ltys[1],lwd=lwds[1],col=cols[1])

      # Add lk.unwgt.data reference line
      abline(v=lk.unwgt.data,lty=ltys[2],lwd=lwds[2],col=cols[2])

      # Add legend
      legend("top",legend=c(paste("lk.total (",signif(lk.total,digits),")",sep=""),
                            paste("lk.unwgt.data (",signif(lk.unwgt.data,digits),")",sep="")),
             bty="n", col=cols, lty=ltys, lwd=lwds)
  }else{
    message("like not found in spp")
  }
}

