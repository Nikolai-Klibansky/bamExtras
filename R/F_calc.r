#' F_calc
#'
#' Calculations of stock dependent quantities over a range of fishing mortalities (F).
#' This code calculates spawning potential ratio (SPR), and then maximum sustainable
#' yield (MSY) using equilibrium methods incorporating a Beverton-Holt stock-recruit relationship.
#' Many of the calculations, particularly the MSY calculations, are based on code
#' in the Beaufort Assessment Model. This function is similar to the
#' \code{get_ref_pts} function, but is not as specific to BAM.
#' @param age age classes. numeric vector
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
#' @param P_st Time of year when spawning occurs, as a proportion.
#' @param plus_group Passed to bamExtras::exp_decay (Should the function include a plus group? logical)
#' @param N_ast_ob Numbers of fish at age, at the time of spawning, observed. Allows the user to provide age structure of the population in the reference year (e.g. terminal year) used to compute observed stock size (i.e. numerator in stock status calculations). If values are provided, additional outputs will be computed in output data, with suffix "ob".
#' @param R_neq Recruitment used to scale non-equilibrium per recruit quantities to compute SPR-based reference points (e.g. stock size at Fi; SSB_40%). An estimate of mean recruitment is often used here.
#' @param Px reference proportion of unfished spawning potential ratio (e.g. Px=.40)
#' @param plots logical. Draw plots?
#' @param plot_digits number of significant digits to show in plot
#' @keywords bam stock assessment fisheries population dynamics
#' @author Erik Williams, Kyle Shertzer, and Nikolai Klibansky
#' @export
#' @examples
#' rdat <- rdat_VermilionSnapper
#' pr <- rdat$parms
#' as <- rdat$a.series
#' F_calc(age=as$age, h=pr$BH.steep, R0=pr$R0,rec_sigma=rdat$parms$R.sigma.par, M=as$M,
#'        sel_L=rdat$sel.age$sel.v.wgted.L, sel_D=rdat$sel.age$sel.v.wgted.D,
#'        PS=rep(1,length(as$age)), w=as$wholewgt.wgted.L.klb, ep=as$reprod,
#'        plots=TRUE,P_st=rdat$parms$spawn.time,plot_digits=4)
#'
#' # When there is no SRR, the function will still compute SPR without computing MSY-based equilibrium values
#' rdat <- rdat_RedSnapper
#' pr <- rdat$parms
#' as <- rdat$a.series
#' F_calc(age=as$age, M=as$M, sel_L=rdat$sel.age$sel.v.wgted.L, sel_D=rdat$sel.age$sel.v.wgted.D,
#'        PS=rep(1,length(as$age)), w=as$wholewgt.wgted.L.klb, ep=as$reprod,
#'        plots=TRUE,P_st=rdat$parms$spawn.time,plot_digits=4)

F_calc <-  function(age,
                    h=NULL,
                    R0=NULL,
                    rec_sigma=0,
                    M,
                    Fmax=3,
                    F_n=101,
                    sel_L,
                    sel_D=NULL,
                    PS,
                    w,
                    ep,
                    P_st=0.5,
                    plus_group=TRUE,
                    N_ast_ob=NULL,
                    R_neq=NA,
                    Px=.40,
                    plots=FALSE,
                    plot_digits=3){
# Check inputs
  if(!is.numeric(age)){warning("age must be a numeric vector")}
  check_vec(M,age,"F_calc")
  check_vec(sel_L,age,"F_calc")
  if(is.null(sel_D)){sel_D <- rep(0,length(age))}
  check_vec(sel_D,age,"F_calc")
  check_vec(PS,age,"F_calc")
  check_vec(w,age,"F_calc")
  check_vec(ep,age,"F_calc")
  if(is.null(N_ast_ob)){N_ast_ob <- rep(NA,length(age))}

F <- seq(0,Fmax,length=F_n)
age_n <- length(age) # Number of age classes

# Initialize storage vectors
spr <- SPR <- rep(0,length(F))
R_eq <- S_eq <- B_eq <- Ln_eq <- Lw_eq <- Dn_eq <- Dw_eq <- E_eq <- rep(NA,length(F))

is_SRR <- !is.null(h)&!is.null(R0) # Are stock recruit relationship parameters provided?
if(!is_SRR){
  warning("Either h or R0 is null. Equilibrium calculations will not be computed.")
}

# data inputs at-age

pdout <- list()

# FOR EACH LEVEL OF F..
for (fi in 1:F_n) {
  Fi <- F[fi]
  pdouti <- pop_demo(
    Fi=Fi,
    age=age,
    h=h,
    R0=R0,
    rec_sigma=rec_sigma,
    M=M,
    sel_L=sel_L,
    sel_D=sel_D,
    PS=PS,
    w=w,
    ep=ep,
    P_st=P_st,
    plus_group=plus_group,
    N_ast_ob=N_ast_ob,
    R_neq=R_neq
    )
  pdout[[fi]] <- pdouti$data_F
}

# data  outputs by F
data_F <- as.data.frame(do.call(rbind,pdout)) # Combine data_F vectors into a matrix with F_n rows

spr <- data_F$spr
SPR <- data_F$SPR
phi0 <- data_F$spr_F0[1]

R_eq <- data_F$R_eq
S_eq <- data_F$S_eq
B_eq <- data_F$B_eq
Lw_eq <- data_F$Lw_eq
E_eq <- data_F$E_eq

  ### Reference points
  # SPR-based
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
          B1 <- (y2-y1)/(x2-x1)         # Slope of the line connecting the adjacent points
          B0 <- y1-B1*x1                # Intercept of the line connecting the adjacent points
          return(B0+B1*SPR_Px)          # F_Px calculated through linear interpolation
          })


  # MSY-based
          if(is_SRR){
            msy <-     max(Lw_eq)          # maximum sustainable yield
            msy_ix <-  which(Lw_eq==msy)   # msy F-index
            Fmsy <-    F[msy_ix]           # fishing rate at msy
            if(Fmsy==Fmax){warning("Fmsy estimate is at the upper bound (Fmax). Try increasing Fmax.")}
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

# Get demographic data associated with F-reference values
pdout_SPR <- pop_demo(
  Fi=F_Px,
  age=age,
  h=h,
  R0=R0,
  rec_sigma=rec_sigma,
  M=M,
  sel_L=sel_L,
  sel_D=sel_D,
  PS=PS,
  w=w,
  ep=ep,
  P_st=P_st,
  plus_group=plus_group,
  N_ast_ob=N_ast_ob,
  R_neq=R_neq
  )
data_age_SPR <- pdout_SPR$data_age

data_age_msy <- if(is_SRR){
  pdout_msy <- pop_demo(
    Fi=Fmsy,
    age=age,
    h=h,
    R0=R0,
    rec_sigma=rec_sigma,
    M=M,
    sel_L=sel_L,
    sel_D=sel_D,
    PS=PS,
    w=w,
    ep=ep,
    P_st=P_st,
    plus_group=plus_group,
    N_ast_ob=N_ast_ob,
    R_neq=R_neq
  )
  pdout_msy$data_age
}else{
  NA
}

if(plots){
# # Values by age
#   par(mar=c(2,2,2,1),mfrow=c(3,3),mgp=c(1.1,0.1,0),tck=0.01)
# with(dt_age,
#      {
#      plot(age,M,type="o",xlab="age",ylab="natural mortality (M)")
#      plot(age,sel_L,type="o",xlab="age",ylab="selectivity of landings (sel_L)")
#      plot(age,sel_D,type="o",xlab="age",ylab="selectivity of discards (sel_D)")
#      plot(age,PS,type="o",xlab="age",ylab="proportion spawning (PS)")
#      plot(age,w,type="o",xlab="age",ylab="weight (w)")
#      plot(age,ep,type="o",xlab="age",ylab="egg production proxy (ep)")
#      }
#      )


# Values by F
  mfrow <- if(is_SRR){c(3,3)}else{c(2,1)}
  par(mar=c(2,2,2,1),mfrow=mfrow,mgp=c(1.1,0.1,0),tck=0.01,cex=1)
  plot_F <- function(...,rp=NA){
    plot(type="l",lwd=2,...)
    if(!is.na(rp)){points(Fmsy,rp,type="p",col="blue",pch=16)}
  }

  plot_F(F,spr,  main="Non-Equilibrium", rp=sprmsy)
  plot_F(F,SPR,  main="Non-Equilibrium", rp=SPRmsy,ylim=c(0,1))
  usr <- par("usr")
  #text(Fmsy,SPRmsy, labels=bquote(SPR[msy] == .(signif(SPRmsy,plot_digits))),pos=4)
  text(Fmsy,SPRmsy, labels=bquote(list(F[MSY] == .(signif(Fmsy,plot_digits)),SPR[MSY] == .(signif(SPRmsy,plot_digits)))),pos=4)
  points(F_Px,SPR_Px,type="p",col="red",pch=16)
  text(F_Px,SPR_Px, labels=bquote(list(F[.(Px)] == .(signif(F_Px,plot_digits)),SPR[.(Px)] == .(signif(SPR_Px,plot_digits)))),pos=4)
  #text(F_Px,usr[3]+diff(usr[3:4])*.05, labels=bquote(F[.(Px)] == .(signif(F_Px,plot_digits))),pos=4)
  if(is_SRR){
    plot_F(F,Lw_eq, main="Equilibrium",     rp=msy)
    text(Fmsy,msy, labels=bquote(MSY == .(signif(msy,plot_digits))),pos=4)
    plot_F(F,R_eq, main="Equilibrium",     rp=Rmsy)
    text(Fmsy,Rmsy, labels=bquote(R[MSY] == .(signif(Rmsy,plot_digits))),pos=4)
    plot_F(F,S_eq, main="Equilibrium",     rp=Smsy)
    text(Fmsy,Smsy, labels=bquote(S[MSY] == .(signif(Smsy,plot_digits))),pos=4)
    plot_F(F,B_eq, main="Equilibrium",     rp=Bmsy)
    text(Fmsy,Bmsy, labels=bquote(B[MSY] == .(signif(Bmsy,plot_digits))),pos=4)
    plot_F(F,E_eq, main="Equilibrium",     rp=Emsy)
    text(Fmsy,Emsy, labels=bquote(E[MSY] == .(signif(Emsy,plot_digits))),pos=4)
  }
}

## 2025-06-20 Updating this function dependent on pop_demo function. Need to tweak pop_demo and add it to bamExtras.
## Also need to output dt_age associated with Fmsy and F_Px, and program some relevant plots associated
## with those dt_age data frames.
invisible(list(
               "RefPts"=ref_points, "data_F"=data_F, "data_age_SPR"=data_age_SPR, "data_age_msy"=data_age_msy
               ))
}
