#' Exponential decay
#'
#' Compute numbers of individuals alive at time in a population experiencing exponential decay, allowing age-varying Z. Can be used to compute numbers at age at the beginning of the year.
#' @param age ages at which to compute numbers alive. Can be fractional ages (e.g. daily, monthly).
#' @param Z instantaneous mortality rate
#' @param N1 number of individuals at first age (i.e. \code{age[1]}) at the beginning of the year
#' @param plus_group Should the function include a plus group? logical
#' @keywords bam stock assessment fisheries population dynamics
#' @author Nikolai Klibansky
#' @export
#' @examples
#' \dontrun{
#' exp_decay(age=1:10,Z=0.2,N1=100)
#' exp_decay(age=1:10,Z=0.2,N1=100,plus_group = TRUE)
#' # The age vector in this example triggers a warning, but in this case it's actually okay.
#' exp_decay(age=c(1,5,10),Z=0.2,N1=100,plus_group = TRUE)
#' }
exp_decay <- function(age,Z,N1=1,plus_group=FALSE){
  if(length(unique(signif(diff(age),2)))>1){
    warning("Age intervals not equivalent. Estimates of numbers at age may be inaccurate")
  }
  ac <- seq_along(age) # Vector of all age classes
  # age <- ac-1
  ac_n <- length(ac)
  N_a <-  setNames(rep(0,ac_n),age)    # Numbers at age
  if(length(Z)==1){
    Z <- c(Z*diff(age),Z) # Z scaled by difference between subsequent age
  }
  a_steps <-  ac_n/(max(ac)-min(ac)+1) # Number of age steps per age class
  # Numbers at age (based on Gabriel et al. 1989)
  for (a in ac){
    N_a[a] <- local({
      if(a==1)   {return(N1)}
      if(a>1)    {
        if(plus_group&a==ac_n){return(N_a[ac_n-1]*exp(-Z[ac_n-1]/a_steps)/   # Plus group
                                        (1-exp(-Z[ac_n]/a_steps)))
        }else{return(N_a[a-1]* exp(-Z[a-1]/a_steps))}
      }
    })
  }
  return(N_a)
}
