#' sr_eq
#'
#' Spawner-recruit equilibrium function (Beverton-Holt or Ricker). Used inside the get_ref_pts function.
#' This is code from the Beaufort Assessment Model converted from ADMB to R.
#' @param steep Beverton-Holt steepness parameter. numeric vector
#' @param R0 Beverton-Holt R0 parameter (virgin recruitment). Numbers of fish at first age (often age-0 or age-1). numeric vector
#' @param BiasCor bias correction. numeric vector
#' @param spr_F0 spawners per recruit at F=0
#' @param spr_F spawners per recruit at F
#' @param SR_switch Set which stock-recruit relationship to use. 1. "beverton_holt" = Beverton-Holt, 2. "ricker" = Ricker, 3. "none" = none (average recruitment).
#' @keywords bam stock assessment fisheries population dynamics
#' @author Kyle Shertzer, Erik Williams, and Nikolai Klibansky
#' @export
#' @examples
#' \dontrun{
#' }


sr_eq <- function(SR_switch="beverton_holt", steep, R0, BiasCor, spr_F0, spr_F){
  switch(SR_switch,
         # Beverton-Holt
         beverton_holt = (R0/((5.0*steep-1.0)*spr_F))*(BiasCor*4.0*steep*spr_F-spr_F0*(1.0-steep)),
         # Ricker
         ricker = R0/(spr_F/spr_F0)*(1.0+log(BiasCor*spr_F/spr_F0)/steep),
         # None
         none = BiasCor*R0
  )
}
