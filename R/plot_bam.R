#' Plot results from Beaufort Assessment Model
#'
#' A generic plotting function that calls plotting functions from FishGraph and bamExtras.
#' @param spp bam output rdat object
#' @param graphics.type from FishGraph help: 'a vector of graphics file types to which graphics are saved. When NULL, no plots are saved.'
#' @param draft from FishGraph help: 'modifies plots for use in a report. When FALSE main titles are omitted.'
#' @param use.color from FishGraph help: 'plots are made in grayscale when FALSE'
#' @param years_plot vector of years to plot in time series. If NULL, all years are plotted
#' @param FGPlot_args set of arguments to pass to all FishGraph plotting functions
#' @keywords bam MCBE stock assessment fisheries
#' @export
#' @examples
#' \dontrun{
#' # Make plots
#' plot_bam(rdat_BlackSeaBass)
#' }

plot_bam <- function(spp=NULL,
                     graphics.type="pdf", # NULL (no quotes) for no plots saved; other options: "pdf", "wmf", "eps"
                     draft=FALSE, # draft type
                     use.color=TRUE,   # color type
                     years_plot = NULL,
                     FGPlot_args=list(),
                     windows_args=list(),
                     Bound.vec.plots.args = list(),
                     BSR.time.plots.args = list(),
                     CLD.total.plots.args = list(),
                     Cohort.plots.args = list(),
                     Comp.plots.args = list(resid.type="pearson"),
                     Comp.yearly.plots.args = list(print.neff=TRUE,print.n=TRUE),
                     Eq.plots.args = list(F.references=list("Fmsy","F40"), user.Eq=list("L.eq.klb", "L.eq.knum", "D.eq.knum")),
                     F.time.plots.args = list(F.references=NULL, F.additional=NULL),
                     Growth.plots.args = list(plot.all = TRUE),
                     Index.plots.args = list(),
                     Landings.plots.args = list(),
                     NFZ.age.plots.args = list(user.plots="N.age.mdyr"),
                     Parm.plots.args = list(),
                     PerRec.plots.args = list(user.PR = list("SPR", "ypr.lb"), user.lab= list("SPR", "ypr.lb") , F.references=list("Fmsy","F40")),
                     Phase.plots.args = list(Xaxis.F=FALSE),
                     Selectivity.plots.args = list(compact=TRUE),
                     StockRec.plots.args = list()
                     )
                     {

library(FishGraph)
  FGPlot_args_user <- FGPlot_args
  FGPlot_args_default <- list(x=quote(spp), DataName="spp",draft=draft, use.color=use.color, graphics.type=graphics.type)
  FGPlot_args <- modifyList(FGPlot_args_default,FGPlot_args_user)

  # set_ts_years()
  # Set years of time series in data frame, when rownames are years
  set_ts_years <- function(x,years){
    # Initialize data frame with desired dimensions
    x <- as.matrix(x)
    x.out <- matrix(NA,nrow=length(years),ncol=ncol(x),
                    byrow=TRUE,dimnames=list(years,colnames(x)))
    x.out[rownames(x)[rownames(x)%in%paste(years)],] <- x[rownames(x)[rownames(x)%in%paste(years)],]
    x.out
  }

# Default settings for setting up plot windows for each FishGraph function
  windows_args_user <- windows_args # Get values specified by the user, if any
  windows_args_default <- list(
    Bound.vec =   list(width =  8, height = 8, record = TRUE),
    BSR.time =    list(width =  8, height = 8, record = TRUE),
    CLD.total =   list(width =  7, height = 5, record = TRUE),
    Cohort =      list(width = 10, height = 8, record = TRUE),
    Comp =        list(width = 7, height = 7, record = TRUE),
    Comp.yearly = list(width = 6, height = 6, record = TRUE),
    Eq =          list(width =  7, height = 5, record = TRUE),
    F.time =      list(width =  7, height = 5, record = TRUE),
    Growth =      list(width =  7, height = 5, record = TRUE),
    Index =       list(width =  7, height = 7, record = TRUE),
    Landings =    list(width =  7, height = 5, record = TRUE),
    NFZ.age =     list(width =  7, height = 5, record = TRUE),
    Parm =        list(width =  8, height = 6, record = TRUE),
    PerRec =      list(width =  7, height = 5, record = TRUE),
    Phase =       list(width =  7, height = 5, record = TRUE),
    Selectivity = list(width =  7, height = 5, record = TRUE, xpos = 10, ypos = 10),
    StockRec =    list(width =  7, height = 5, record = TRUE, xpos = 10, ypos = 10)
  )
  windows_args <- modifyList(windows_args_default,windows_args_user)

styr <- spp$parms$styr
endyr <- spp$parms$endyr

if(is.null(years_plot)){
  years_plot <- as.numeric(styr:endyr)
}


  ### rdat file
  # Modify spp
  # This sorting is helpful because apparently FishGraph gives you errors if the observed comps for a fleet aren't immediately followed by the predicted
  # comps for that same fleet.
  spp$comp.mats <- spp$comp.mats[sort(names(spp$comp.mats))]

  # Trim time series
  spp$N.age <- spp$N.age[rownames(spp$N.age)>=spp$parms$styr,]
  spp$N.age.mdyr <- spp$N.age.mdyr[rownames(spp$N.age.mdyr)>=spp$parms$styr,]
  spp$N.age.spawn <- spp$N.age.spawn[rownames(spp$N.age.spawn)>=spp$parms$styr,]

  # Modify time series to plot a certain set of years
  spp$t.series <- as.data.frame.matrix(set_ts_years(spp$t.series,years_plot))

  spp$CLD.est.mats <- lapply(spp$CLD.est.mats,FUN=set_ts_years,years=years_plot)

  spp$N.age <- set_ts_years(spp$N.age,years_plot)
  spp$N.age.mdyr <- set_ts_years(spp$N.age.mdyr,years_plot)
  spp$N.age.spawn <- set_ts_years(spp$N.age.spawn,years_plot)
  spp$B.age <- set_ts_years(spp$B.age,years_plot)

  # Truncate F-range in eq.series to make nicer plots
  spp$eq.series <- spp$eq.series[which(spp$eq.series$F.eq<=1),]

  # Calculate values to add to (CLD) plots
  spp$parms$Dmsy.num <- spp$parms$Dmsy.knum*1000 # Add value of Dmsy in numbers
  spp$parms$msy.num <- spp$parms$msy.knum*1000 # Add value of msy in numbers

  #Round Effective sample sizes for more enjoyable comp
  spp$t.series[,grepl(pattern="comp",names(spp$t.series))&grepl(pattern="neff",names(spp$t.series))] <-
    round(spp$t.series[,grepl(pattern="comp",names(spp$t.series))&grepl(pattern="neff",names(spp$t.series))],2)

  ###### PLOTAPALOOZA !!! ######
  # Generic plot function for calling various FishGraph plot functions, with a couple of tweaks
  plot_fn <- function(fn){
    cs <-   'if(!is.null(fn.args)){
    cat("\nRunning fn\n")
    FGPlot_args <- FGPlot_args[names(FGPlot_args)%in%names(formals(fn))]
    tryCatch(expr=do.call(FishGraph::fn,c(FGPlot_args,fn.args)),
             error = function(e) {message(paste("Error in fn:",e))}
    )

  }'
    eval(parse(text=gsub("fn",fn,cs)))
  }

  ### Bound.vec.plots ###
  do.call(windows,args=windows_args$Bound.vec)
  plot_fn("Bound.vec.plots")

  ### BSR ###
  do.call(windows,args=windows_args$BSR.time)
  plot_fn("BSR.time.plots")

  ### CLD ###
  do.call(windows,args=windows_args$CLD.total)
  plot_fn("CLD.total.plots")

  ### comp ###
  do.call(windows,args=windows_args$Comp)
  plot_fn("Comp.plots")

  ### compyr ###
  do.call(windows,args=windows_args$Comp.yearly)
  plot_fn("Comp.yearly.plots")

  do.call(windows,args=windows_args$Cohort)
  plot_fn("Cohort.plots")
  # if(!is.null(Cohort.plots.args)){
  #   cat("\nRunning Cohort.plots\n")
  #   # only include args that are in the function formals
  #   FGPlot_args_Cohort <- FGPlot_args[names(FGPlot_args)%in%names(formals(Cohort.plots))]
  #   tryCatch(expr=do.call(FishGraph::Cohort.plots,c(FGPlot_args_Cohort,Cohort.plots.args)),
  #            error = function(e) {message(paste("Error in Cohort.plots:",e))}
  #   )
  # }

  ### EQ ###
  do.call(windows,args=windows_args$Eq)
  plot_fn("Eq.plots")

  ### F ###
  do.call(windows,args=windows_args$F.time)
  plot_fn("F.time.plots")

  ### growth ###
  do.call(windows,args=windows_args$Growth)
  plot_fn("Growth.plots")

  ### index ###
  do.call(windows,args=windows_args$Index)
  plot_fn("Index.plots")

  ### L ###
  do.call(windows,args=windows_args$Landings)
  plot_fn("Landings.plots")

  ### NFZ.age ###
  do.call(windows,args=windows_args$NFZ.age)
  plot_fn("NFZ.age.plots")

  ### parms ###
  do.call(windows,args=windows_args$Parm)
  plot_fn("Parm.plots")

  ### phase ###
  do.call(windows,args=windows_args$Phase)
  plot_fn("Phase.plots")

  ### PR ###
  do.call(windows,args=windows_args$PerRec)
  plot_fn("PerRec.plots")

  ### sel ###
  do.call(windows,args=windows_args$Selectivity)
  plot_fn("Selectivity.plots")

  ### SR ###
  do.call(windows,args=windows_args$StockRec)
  plot_fn("StockRec.plots")

  graphics.off()

}
