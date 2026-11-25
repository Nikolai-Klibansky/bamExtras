#' Plot immature proportion of removals over time
#'
#' Plot immature proportion of removals over time
#' @param spp bam output rdat object
#' @param args_plot list of arguments to pass to \code{graphics::plot}
#' @param nm_L Name of landings matrix in spp$CLD.est.mats (e.g. "Ln.total", "Lw.total")
#' @param nm_D Name of dead discard matrix in spp$CLD.est.mats (e.g. "Dn.total", "Dw.total")
#' @param nm_mat Name of vector of proportion mature-at-age in spp$a.series (e.g. "mat.female")
#' @keywords bam stock assessment fisheries
#' @author Nikolai Klibansky
#' @export
#' @examples
#' plot_imm_prop(rdat_AtlanticMenhaden)
#' plot_imm_prop(rdat_BlackSeaBass)
#' plot_imm_prop(rdat_BluelineTilefish)
#' plot_imm_prop(rdat_Cobia)
#' plot_imm_prop(rdat_GagGrouper)
#' plot_imm_prop(rdat_GrayTriggerfish)
#' plot_imm_prop(rdat_GreaterAmberjack)
#' plot_imm_prop(rdat_RedGrouper)
#' plot_imm_prop(rdat_RedPorgy,nm_mat="mat.fem.endyr")
#' plot_imm_prop(rdat_RedSnapper)
#' plot_imm_prop(rdat_ScampGrouper)
#' plot_imm_prop(rdat_SnowyGrouper)
#' plot_imm_prop(rdat_SpanishMackerel)
#' plot_imm_prop(rdat_Tilefish)
#' plot_imm_prop(rdat_VermilionSnapper)
#'

plot_imm_prop <- function(spp,
                     nm_L="Ln.total",
                     nm_D="Dn.total",
                     nm_mat="mat.female",
                     args_plot=list()
){
  args_plot_user <- args_plot
  args_plot_default <- list(xlab="year",ylab="proportion of removals immature (by numbers)",type="o")
  args_plot <- modifyList(args_plot_default,args_plot_user)

  L_ya <- spp$CLD.est.mats[[nm_L]]
  if(is.null(nm_D)|!nm_D%in%names(spp$CLD.est.mats)){
    message("nm_D not found in names(spp$CLD.est.mats). Discards assumed to be zero.")
    D_ya <- L_ya*0
  }else{
    D_ya <- spp$CLD.est.mats[[nm_D]]
  }

  R_ya <- L_ya + D_ya
  pmat_a <- spp$a.series[,nm_mat]
  pimm_a <- 1-pmat_a
  R_imm_ya <- t(t(R_ya)*pimm_a)
  pimm_y <- rowSums(R_imm_ya)/rowSums(R_ya)

  do.call(plot,c(list(x=as.numeric(names(pimm_y)),y=as.numeric(pimm_y)),args_plot))
  invisible(pimm_y)
  }
