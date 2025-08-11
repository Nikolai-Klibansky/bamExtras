#' pop_demo
#'
#' Computes population demographics at age and totals across age, for a single value of fishing mortality (F).
#' This function is called by \code{bamExtras::F_calc} to compute values over a range of F and compute reference points.
#' \code{pop_demo} is useful to return and plot the demographic data-at-age for specific values of F (e.g. Fmsy)
#' @param Fi Fishing mortality. numeric vector of length 1.
#' @param age age classes. numeric vector
#' @param h Beverton-Holt steepness parameter
#' @param R0 Beverton-Holt R0 parameter. Numbers of fish at age-a (often age-0 or age-1).
#' @param rec_sigma recruitment standard deviation in log space. used to compute lognormal bias correction -- exp(sigma^2/2)
#' @param M Natural mortality rate
#' @param sel_L Selectivity at age for landings
#' @param sel_D Selectivity at age for discards. If NULL, assumed to be zero for all ages.
#' @param PS Proportion of fish to include in stock calculation, at age (e.g. mature females at age)
#' @param w weight at age in any units, used to convert numbers-at-age to weight-at-age. Results in weight will be in these units.
#' @param ep egg production proxy at age (e.g. weight, fecundity)
#' @param P_st Time of year when spawning occurs, as a proportion.
#' @param plus_group Passed to bamExtras::exp_decay (Should the function include a plus group? logical)
#' @param N_ast_ob Numbers of fish at age, at the time of spawning, observed. Allows the user to provide age structure of the population in the reference year (e.g. terminal year) used to compute observed stock size (i.e. numerator in stock status calculations). If values are provided, additional outputs will be computed in output data, with suffix "ob".
#' @param R_neq Recruitment used to scale non-equilibrium per recruit quantities to compute SPR-based reference points (e.g. stock size at Fi; SSB_40%). An estimate of mean recruitment is often used here.
#' @param plots logical. Draw plots?
#' @param plot_digits number of significant digits to show in plot
#' @keywords stock assessment fisheries population dynamics
#' @author Nikolai Klibansky, Erik Williams, and Kyle Shertzer
#' @export
#' @examples
#' # Example of how to use this function
#'

pop_demo <-  function(Fi,
                      age,
                      h=NULL,
                      R0=NULL,
                      rec_sigma=0,
                      M,
                      sel_L,
                      sel_D=NULL,
                      PS,
                      w,
                      ep,
                      P_st=0.5,
                      plus_group=TRUE,
                      N_ast_ob=NULL,
                      R_neq=NA,
                      plots=FALSE,
                      plot_digits=3

){
  # Check inputs
  if(!is.numeric(age)){warning("age must be a numeric vector")}
  check_vec(M,age,"pop_demo")
  check_vec(sel_L,age,"pop_demo")
  if(is.null(sel_D)){sel_D <- rep(0,length(age))}
  check_vec(sel_D,age,"pop_demo")
  check_vec(PS,age,"pop_demo")
  check_vec(w,age,"pop_demo")
  check_vec(ep,age,"pop_demo")
  if(is.null(N_ast_ob)){N_ast_ob <- rep(NA,length(age))}

  age_n <- length(age) # Number of age classes
  age_steps <- age_n/(max(age)-min(age)+1)  # Number of age steps per age

  # ALL FISH
  F_L <- Fi * sel_L
  F_D <- Fi * sel_D
  Z <- M+F_L+F_D

  # PER RECRUIT (R = 1) CALCULATIONS
  # Numbers at-age at the beginning of the year
  N_a_F0_pr <- exp_decay(age=age,Z=M,N1=1,plus_group=plus_group) # when F=0
  N_a_pr <- exp_decay(age=age,Z=Z,N1=1,plus_group=plus_group)    # when F>0

  # Numbers at-age at the time of spawning
  N_ast_F0_pr <- N_a_F0_pr * exp(-M * P_st) # when F=0
  N_ast_pr    <- N_a_pr * exp(-Z * P_st) # when F>0
  if(plus_group){
    N_ast_F0_pr[age_n] <- (N_ast_F0_pr[age_n-1] * (exp(-(M[age_n-1] * (1.0-P_st) + M[age_n] * P_st) )))/(1.0-exp(-M[age_n]))
    N_ast_pr[age_n] <-    (N_ast_pr[age_n-1]    * (exp(-(Z[age_n-1] * (1.0-P_st) + Z[age_n] * P_st) )))/(1.0-exp(-Z[age_n]))
  }

  # Number of fish in the spawning stock (e.g. number of mature females) at age, at the time of spawning
  N_S_ast_F0_pr <- N_ast_F0_pr * PS # when F=0
  N_S_ast_pr <-    N_ast_pr    * PS # when F>0
  N_S_ast_ob <- N_ast_ob * PS # observed value

  # Spawning stock at age
  S_ast_F0_pr <- N_S_ast_F0_pr * ep
  S_ast_pr <-    N_S_ast_pr    * ep
  S_ast_ob <- N_S_ast_ob * ep

  # Spawning stock summed across ages
  S_F0_pr <- spr_F0 <- sum(S_ast_F0_pr) # when F=0, per recruit
  S_pr <- spr <- sum(S_ast_pr)          # when F>0, per recruit
  S_ob <- sum(S_ast_ob)           # observed value

  # Unfished spawning biomass per recruit
  phi0 <- spr_F0
  # Spawning potential ratio
  SPR <- spr/spr_F0

  # NON-EQUILIBRIUM (R_neq based) CALCULATIONS
  # Numbers at age
  N_a_neq <-   R_neq * N_a_pr
  N_ast_neq <- R_neq * N_ast_pr
  # Number of fish in the spawning stock (e.g. number of mature females) at age
  N_S_ast_neq <- N_ast_neq * PS
  # Spawning stock (biomass, eggs, or possible other units)
  S_ast_neq <- N_S_ast_neq * ep
  S_neq <- sum(S_ast_neq)
  # Total biomass
  B_a_neq <- N_a_neq * w
  B_neq <- sum(B_a_neq)
  # Landings (numbers) at age (Gabriel et al. 1989; BAM code)
  L_a_neq <- N_a_neq * (F_L/Z) * (1-exp(-Z/age_steps))
  Ln_neq <- sum(L_a_neq) # a.k.a. yield
  Lw_neq <- sum(L_a_neq * w) # a.k.a. yield
  # Discards (numbers) at age
  D_a_neq <- N_a_neq * (F_D/Z) * (1-exp(-Z/age_steps))
  Dn_neq <- sum(D_a_neq)
  Dw_neq <- sum(D_a_neq * w)
  # Exploitation rate (total catch/number of fish)
  E_neq <- (Ln_neq+Dn_neq)/sum(N_a_neq)

  # EQUILIBRIUM (R_eq based) CALCULATIONS
  R_eq <- S_eq <- B_eq <- Ln_eq <- Lw_eq <- Dn_eq <- Dw_eq <- E_eq <- NA
  N_a_eq <- N_ast_eq <- N_S_ast_eq <- S_ast_eq <- B_a_eq <- L_a_eq <- D_a_eq <- age*NA
  is_SRR <- !is.null(h)&!is.null(R0) # Are stock recruit relationship parameters provided?
  if(is_SRR){
    # Incorporate the SRR into the following calculations
    # Recruitment (equilibrium recruitment based on SSB/R at F, and R_eq)
    R_eq <- local({
      BC<-exp(rec_sigma^2/2.0)              # multiplicative bias correction
      R_eq <- (R0/((5*h-1.0)*spr))*(BC*4.0*h*spr-phi0 *(1.0-h))
      ifelse(R_eq<1e-7, 1e-7, R_eq)   # Keeps R_eq from getting too close to zero.
    })
    # Numbers at age
    N_a_eq   <- R_eq * N_a_pr
    N_ast_eq <- R_eq * N_ast_pr
    # Number of fish in the spawning stock (e.g. number of mature females) at age
    N_S_ast_eq <- N_ast_eq * PS
    # Spawning stock (biomass, eggs, or possible other units)
    S_ast_eq <- N_S_ast_eq * ep
    S_eq <- sum(S_ast_eq)
    # Total biomass
    B_a_eq <- N_a_eq * w
    B_eq <- sum(B_a_eq)
    # Landings (numbers) at age (Gabriel et al. 1989; BAM code)
    L_a_eq <- N_a_eq * (F_L/Z) * (1-exp(-Z/age_steps))
    Ln_eq <- sum(L_a_eq) # a.k.a. yield
    Lw_eq <- sum(L_a_eq * w) # a.k.a. yield
    # Discards (numbers) at age
    D_a_eq <- N_a_eq * (F_D/Z) * (1-exp(-Z/age_steps))
    Dn_eq <- sum(D_a_eq)
    Dw_eq <- sum(D_a_eq * w)
    # Exploitation rate (total catch/number of fish)
    E_eq <- (Ln_eq+Dn_eq)/sum(N_a_eq)
  }

  # data at-age
  data_age <- data.frame(# inputs
                         age = age,
                         M = M,
                         sel_L = sel_L,
                         sel_D = sel_D,
                         PS = PS,
                         w = w,
                         ep = ep,

                         # computed for SPR or MSY reference points
                         F_L = F_L,
                         F_D = F_D,
                         Z = Z,

                         # computed for SPR quantities at F=0, per recruit
                         N_a_F0_pr = N_a_F0_pr,
                         N_ast_F0_pr = N_ast_F0_pr,
                         N_S_ast_F0_pr = N_S_ast_F0_pr,
                         S_ast_F0_pr = S_ast_F0_pr,

                         # computed for SPR quantities at F>0, per recruit
                         N_a_pr = N_a_pr,
                         N_ast_pr = N_ast_pr,
                         N_S_ast_pr = N_S_ast_pr,
                         S_ast_pr = S_ast_pr,

                         # computed for SPR quantities at F>0, scaled to R_neq
                         N_a_neq = N_a_neq,
                         N_ast_neq = N_ast_neq,
                         N_S_ast_neq = N_S_ast_neq,
                         S_ast_neq = S_ast_neq,

                         # computed for observed age structure
                         N_ast_ob = N_ast_ob,
                         N_S_ast_ob = N_S_ast_ob,
                         S_ast_ob = S_ast_ob,

                         # computed for MSY quantities (equilibrium recruitment based on SRR)
                         N_a_eq = N_a_eq,
                         N_ast_eq = N_ast_eq,
                         N_S_ast_eq = N_S_ast_eq,
                         S_ast_eq = S_ast_eq,
                         B_a_eq = B_a_eq,
                         L_a_eq = L_a_eq,
                         D_a_eq = D_a_eq
                         )

  # data at-F
  data_F <-
     c("F"=Fi,
       spr_F0=spr_F0,
       spr=spr,
       SPR=SPR,

       S_ob=S_ob,

       R_neq=R_neq,
       S_neq=S_neq,
       B_neq=B_neq,
       Ln_neq=Ln_neq,
       Lw_neq=Lw_neq,
       Dn_neq=Dn_neq,
       Dw_neq=Dw_neq,
       E_neq=E_neq,

       R_eq=R_eq,
       S_eq=S_eq,
       B_eq=B_eq,
       Ln_eq=Ln_eq,
       Lw_eq=Lw_eq,
       Dn_eq=Dn_eq,
       Dw_eq=Dw_eq,
       E_eq=E_eq
  )

  invisible(list(data_age=data_age,data_F=data_F))
}
