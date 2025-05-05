#' F_calc
#'
#' Calculations of stock dependent quantities over a range of fishing mortalities (F).
#' This code calculates spawning potential ratio (SPR), and then maximum sustainable
#' yield (MSY) using equilibrium methods incorporating a Beverton-Holt stock-recruit relationship.
#' Many of the calculations, particularly the MSY calculations, are based on code
#' in the Beaufort Assessment Model. This function is essentially the same as the
#' \code{get_msy} function, but is not as specific to BAM.
#' @param ac age classes. numeric vector
#' @param h Beverton-Holt steepness parameter
#' @param R0 Beverton-Holt R0 parameter. Numbers of fish at age-a (often age-0 or age-1).
#' @param rec_sigma recruitment standard deviation in log space. used to compute lognormal bias correction -- exp(sigma^2/2)
#' @param M Natural mortality rate
#' @param Fmax Maximum fishing mortality rate
#' @param F_n Number of fishing mortality rates to try
#' @param sel_L Selectivity at age for landings
#' @param sel_D Selectivity at age for discards
#' @param PS Proportion of fish to include in stock calculation, at age (e.g. mature females at age)
#' @param w weight at age in any units, used to convert numbers-at-age to weight-at-age. Results in weight will be in these units.
#' @param ep egg production proxy at age (e.g. weight, fecundity)
#' @param Px reference proportion of unfished spawning potential ratio (e.g. Px=.40)
#' @param P_st Time of year when spawning occurs, as a proportion.
#' @param plots logical. Draw plots?
#' @param plot_digits number of significant digits to show in plot
#' @keywords bam stock assessment fisheries population dynamics
#' @author Erik Williams, Kyle Shertzer, and Nikolai Klibansky
#' @export
#' @examples
#' rdat <- rdat_VermilionSnapper
#' pr <- rdat$parms
#' as <- rdat$a.series
#' F_calc(ac=as$age, h=pr$BH.steep, R0=pr$R0,rec_sigma=rdat$parms$R.sigma.par, M=as$M,
#'        sel_L=rdat$sel.age$sel.v.wgted.L, sel_D=rdat$sel.age$sel.v.wgted.D,
#'        PS=rep(1,length(as$age)), w=as$wholewgt.wgted.L.klb, ep=as$reprod,
#'        plots=TRUE,P_st=rdat$parms$spawn.time,plot_digits=4)
#'
#' # When there is no SRR, the function will still compute SPR without computing MSY-based equilibrium values
#' rdat <- rdat_RedSnapper
#' pr <- rdat$parms
#' as <- rdat$a.series
#' F_calc(ac=as$age, M=as$M, sel_L=rdat$sel.age$sel.v.wgted.L, sel_D=rdat$sel.age$sel.v.wgted.D,
#'        PS=rep(1,length(as$age)), w=as$wholewgt.wgted.L.klb, ep=as$reprod,
#'        plots=TRUE,P_st=rdat$parms$spawn.time,plot_digits=4)

F_calc <-  function(ac,
                    h=NULL,
                    R0=NULL,
                    rec_sigma=0,
                    M,
                    Fmax=1,
                    F_n=101,
                    sel_L,
                    sel_D,
                    PS,
                    w,
                    ep,
                    Px=.40,
                    P_st=0.5,
                    plots=FALSE,
                    plus_group=TRUE,
                    plot_digits=3){

F <- seq(0,Fmax,length=F_n)
ac_n=length(ac) # Number of age classes
a_steps <- ac_n/(max(ac)-min(ac)+1)  # Number of age steps per age

spr <- SPR <- R_eq <- S_eq <- B_eq <- Ln_eq <- Lw_eq <- Dn_eq <- Dw_eq <- E_eq <- rep(0,length(F)) # Initialize storage vectors

is_SRR <- !is.null(h)&!is.null(R0) # Are stock recruit relationship parameters provided?
if(!is_SRR){
  warning("Either h or R0 is null. Equilibrium calculations will not be computed.")
}

# FOR EACH LEVEL OF F..
for (fi in 1:F_n) {
  # NON-EQUILIBRIUM CALCULATIONS
  # ALL FISH
  F_L <- F[fi]*sel_L
  F_D <- F[fi]*sel_D
  Z <- M+F_L+F_D
  # Numbers at-age at the beginning of the year
  N_a <- exp_decay(age=ac,Z=Z,N1=1,plus_group=plus_group)
  # Numbers at-age at the time of spawning
  N_ast <- N_a*exp(-Z*P_st)
  if(plus_group){
    N_ast[ac_n] <- (N_ast[ac_n-1]*(exp(-(Z[ac_n-1]*(1.0-P_st) + Z[ac_n]*P_st) )))/(1.0-exp(-Z[ac_n]))
  }

  # Numbers at age in the unfished population (when F=0)
  if(F[fi]==0){
    N_a_F0 <- data.frame("age_class"=ac, "N_F0"=N_a)
  }

  # Number of fish in the spawning stock (e.g. number of mature females) at age, at the time of spawning
  N_S_ast <- N_ast*PS
  # Spawners per recruit
  spr[fi] <- sum(N_S_ast*ep)
  # Unfished spawning biomass per recruit
  phi0 <- spr[1]
  # Spawning potential ratio
  SPR[fi] <- spr[fi]/spr[1]

  # EQUILIBRIUM CALCULATIONS
  if(is_SRR){
    # Incorporate the SRR into the following calculations
    # Recruitment (equilibrium recruitment based on SSB/R at F, and R_eq)
    R_eq[fi] <- local({
      BC<-exp(rec_sigma^2/2.0)              # multiplicative bias correction
      #R_eq_i <- (R0*(spr[fi]*4*h*BC-phi0*(1-h)))/(spr[fi]*(5*h-1))
      R_eq_i <- (R0/((5*h-1.0)*spr[fi]))*(BC*4.0*h*spr[fi]-phi0 *(1.0-h))
      ifelse(R_eq_i<1e-7, 1e-7, R_eq_i)   # I think this just keeps R_eq_i from getting too close to zero.
    })
    # Numbers at age
    N_a_eq   <- R_eq[fi]*N_a
    N_ast_eq <- R_eq[fi]*N_ast
    # Number of fish in the spawning stock (e.g. number of mature females) at age
    N_S_ast_eq <- N_ast_eq*PS
    # Spawning stock (biomass, eggs, or possible other units)
    S_eq[fi] <- sum(N_S_ast_eq*ep)
    # Total biomass
    B_eq[fi] <- sum(N_a_eq*w)
    # Landings (numbers) at age (Gabriel et al. 1989; BAM code)
    L_a_eq <- N_a_eq*(F_L/Z)*(1-exp(-Z/a_steps))
    Ln_eq[fi] <- sum(L_a_eq) # a.k.a. yield
    Lw_eq[fi] <- sum(L_a_eq*w) # a.k.a. yield
    # Discards (numbers) at age
    D_a_eq <- N_a_eq*(F_D/Z)*(1-exp(-Z/a_steps))
    Dn_eq[fi] <- sum(D_a_eq)
    Dw_eq[fi] <- sum(D_a_eq*w)
    # Exploitation rate (total catch/number of fish)
    E_eq[fi] <- (Ln_eq[fi]+Dn_eq[fi])/sum(N_a_eq)
  }
}

  # Reference points
    SPR_Px<-Px  # SPR_Px
    #SPR_Px<<-which(abs(SPR-Px)==min(abs(SPR-Px)))  # Index value corresponding to SPR_Px
    # Calculate close approximation of F at SPR_Px, by interpolating between the nearest values
        x1_ix <- which(SPR==max(SPR[SPR<=SPR_Px]))
        x2_ix <- which(SPR==min(SPR[SPR>SPR_Px]))
        x1 <- SPR[x1_ix]  # Value just below SPR_Px
        x2 <- SPR[x2_ix]   # Value just above SPR_Px

        # F_Px
          F_Px <- local({
          y1 <- F[x1_ix]      # F at x1
          y2 <- F[x2_ix]      # F at x2
          B1 <- (y2-y1)/(x2-x1)            # Slope of the line connecting the adjacent points
          B0 <- y1-B1*x1                   # Intercept of the line connecting the adjacent points
          return(B0+B1*SPR_Px)          # F_Px calculated through linear interpolation
          })


    #F_Px <-    F[SPR_Px_ix]    # F at SPR reference value
    #SPR_Px <-  SPR[SPR_Px_ix]  # SPR at SPR reference value
          if(is_SRR){
            msy <-     max(Lw_eq)          # maximum sustainable yield
            msy_ix <-  which(Lw_eq==msy)   # msy F-index
            Fmsy <-    F[msy_ix]           # fishing rate at msy
            sprmsy <-  spr[msy_ix]         # spawners per recruit at msy
            SPRmsy <-  sprmsy/phi0         # spawning potential ratio at msy
            Rmsy <-    R_eq[msy_ix]        # equilibrium recruitment at msy
            Smsy <-    S_eq[msy_ix]        # spawning stock at msy (e.g. SSB, fecundity)
            Bmsy <-    B_eq[msy_ix]        # total biomass (male and female) at msy
            Emsy <-    E_eq[msy_ix]        # exploitation rate at msy
          }else{
            msy <- Fmsy <- sprmsy <- SPRmsy <- Rmsy <- Smsy <- Bmsy <- Emsy <- NA
          }

ref_points <- if(sum(ep)>0){data.frame(F_Px, SPR_Px, msy, Fmsy, sprmsy, SPRmsy, Rmsy, Smsy, Bmsy, Emsy)
           }else{data.frame("F_Px"=NA,"SPR_Px"=SPR_Px,"msy"=NA,"Fmsy"=NA,"sprmsy"=NA,"SPRmsy"=NA,
                            "Rmsy"=NA,"Smsy"=NA,"Bmsy"=NA,"Emsy"=NA)}
data <- data.frame("F"=F, "spr"=spr, "SPR"=SPR,
                "R_eq"=R_eq, "S_eq"=S_eq, "B_eq"=B_eq,
                "Lw_eq"=Lw_eq, "Ln_eq"=Ln_eq,
                "Dw_eq"=Dw_eq, "Dn_eq"=Dn_eq,
                "E_eq"=E_eq)
if(plots){
  par(mar=c(2,2,2,1),mfrow=c(3,3),mgp=c(1.1,0.1,0),tck=0.01)
  plotNK <- function(...,rp=NA){
    plot(type="l",lwd=2,...)
    if(!is.na(rp)){points(Fmsy,rp,type="p",col="blue",pch=16)}
  }

  plotNK(F,spr,  main="Non-Equilibrium", rp=sprmsy)
  plotNK(F,SPR,  main="Non-Equilibrium", rp=SPRmsy)
  usr <- par("usr")
  #text(Fmsy,SPRmsy, labels=bquote(SPR[msy] == .(signif(SPRmsy,plot_digits))),pos=4)
  text(Fmsy,SPRmsy, labels=bquote(list(F[MSY] == .(signif(Fmsy,plot_digits)),SPR[MSY] == .(signif(SPRmsy,plot_digits)))),pos=4)
  points(F_Px,SPR_Px,type="p",col="red",pch=16)
  text(F_Px,SPR_Px, labels=bquote(list(F[.(Px)] == .(signif(F_Px,plot_digits)),SPR[.(Px)] == .(signif(SPR_Px,plot_digits)))),pos=4)
  #text(F_Px,usr[3]+diff(usr[3:4])*.05, labels=bquote(F[.(Px)] == .(signif(F_Px,plot_digits))),pos=4)
  if(is_SRR){
    plotNK(F,Lw_eq, main="Equilibrium",     rp=msy)
    text(Fmsy,msy, labels=bquote(MSY == .(signif(msy,plot_digits))),pos=4)
    plotNK(F,R_eq, main="Equilibrium",     rp=Rmsy)
    text(Fmsy,Rmsy, labels=bquote(R[MSY] == .(signif(Rmsy,plot_digits))),pos=4)
    plotNK(F,S_eq, main="Equilibrium",     rp=Smsy)
    text(Fmsy,Smsy, labels=bquote(S[MSY] == .(signif(Smsy,plot_digits))),pos=4)
    plotNK(F,B_eq, main="Equilibrium",     rp=Bmsy)
    text(Fmsy,Bmsy, labels=bquote(B[MSY] == .(signif(Bmsy,plot_digits))),pos=4)
    plotNK(F,E_eq, main="Equilibrium",     rp=Emsy)
    text(Fmsy,Emsy, labels=bquote(E[MSY] == .(signif(Emsy,plot_digits))),pos=4)
  }
}

invisible(list("N_a_F0"=N_a_F0, "RefPts"=ref_points, "Data"=data))
}
