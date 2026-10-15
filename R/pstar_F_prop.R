#' Function to compute proportion of Fmsy associated with a particular Pstar value, (a proportion between 0 and 1) 
#'
#' This method was developed by Shertzer et al. (2010) to compute the level of F associated with a value Pstar, which is defined as the allowable probability of overfishing.
#' @param Pstar numeric vector. The pstar value the SSC provides to compute Fmsy scaling percentage  (when pstart = 0.5 should match the resulting shift or F.scale value should be 1 and therefore Fmsy isn't really scaled. With a value below 0.5, the resulting scalar will be below 1 and therefore the scaled Fmsy value will be lower and therefore more conservative)
#' @param Fmsy_boot Vector of bootstrap Fmsy values from uncertainty analysis (MCBE). The length of Fmsy_boot is the number of bootstrap runs (possibly filtered to remove nonconvergant or extreme runs)
#' @references Shertzer, K. W., M. H. Prager, and E. H. Williams. 2010. Probabilistic approaches to setting acceptable biological catch and annual catch targets for multiple years: reconciling methodology with national standards guidelines. Marine and Coastal Fisheries 2:451-458.
#' @returns Proportion of FMSY to use to compute Ftarget
#' @keywords bam stock assessment fisheries
#' @export
#' @examples
#' \dontrun{
#' }
#'
pstar_F_prop <- function(Pstar,Fmsy_boot){
  den_Fmsy <- density(Fmsy_boot, n=2^10, from=min(Fmsy_boot)) # Kernal density of bootstrapped Fmsy
  P_Fmsy <- den_Fmsy$y/sum(den_Fmsy$y)      
  Fmsy <- den_Fmsy$x  
  

# Finds the shift (s) in Fmsy distribution that provides Pstar, using double integration method
Pstar_dble_int  <-  function(s, Pstar, den_Fmsy, Fmsy,P_Fmsy){
  Fmsy_s <- s*Fmsy
  dFmsy <- round(den_Fmsy$x[2]-den_Fmsy$x[1],8) # gets width of each interval between Fmsy values
  
  Pstar_est <- 0.0 # Initialize estimate of Pstar
  # Double integral calculation of Pstar
  # I think this is Eq. 3 of Shertzer et al 2010.
  for (i in 1:length(Fmsy)) {
    Pstar_est <- Pstar_est + sum(den_Fmsy$y[Fmsy_s>Fmsy[i]] *dFmsy) * P_Fmsy[i]   
  }  
  val <- (Pstar_est-Pstar)^2 # objective function value
  return(val)
}


fit <- optimize(f=Pstar_dble_int,
                interval=c(0.0, 1.0001),
                Pstar=Pstar, den_Fmsy=den_Fmsy, Fmsy=Fmsy, P_Fmsy=P_Fmsy,
                maximum = FALSE,
                tol = 0.000001
                )
return(fit$minimum)
}