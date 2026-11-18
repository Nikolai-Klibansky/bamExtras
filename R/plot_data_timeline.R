#' plot_data_timeline
#'
#' Make a data timeline table from the BAM rdat output t.series object. Note: plot_data_timeline() uses the plot.matrix package
#' @param rdat (list) object read in with dget()
#' @param key_data_type_abb data type abbreviation key
#' @param key_fleet_abb fleet abbreviation key
#' @param data_types  data types to include in plot
#' @param t.series.row.trunc truncate rows of t.series object by this value to avoid plotting projection year or show what available data would look like in retrospective runs
#' @keywords bam stock assessment fisheries population dynamics
#' @author Nikolai Klibansky
#' @export
#' @examples
#' \dontrun{
#' ## Data timeline plot Examples
#' par(mar=c(3,12,1,8),mgp=c(1,0.3,0),tck=0,las=2)
#' rdat <- rdat_GrayTriggerfish
#' plot_data_timeline(rdat,FUN_col = "heat.colors")
#'
#' par(mar=c(3,12,1,8),mgp=c(1,0.3,0),tck=0,las=2)
#' library(viridisLite)
#' plot_data_timeline(rdat,data_types = c("acomp","lcomp"),FUN_col = "viridis")
#' }
#'


# Data input timeline
plot_data_timeline <- function(rdat,
                               key_data_type_abb=c("L"="landings",
                                                "D"="discards",
                                                "U"="index",
                                                "acomp"="age comp",
                                                "lcomp"="length comp"),
                               key_fleet_abb=c(
                                 "cGN"   = "com general",
                                 "cGD"   = "com general disc",                                                                             "cHL"   = "com handline",
                                 "cHD"   = "com handline discards",
                                 "cPT"   = "com pots",
                                 "cPD"   = "com pots discards",
                                 "cOT"   = "com other",
                                 "cOD"   = "com other discards",
                                 "cLL"   = "com longline",
                                 "cLD"   = "com longline discards",
                                 "rGN"   = "rec general",
                                 "rGD"   = "rec general discards",
                                 "rHB"   = "rec headboat",
                                 "rHD"   = "rec headboat discards",
                                 "sTV"   = "SERFS trap video"),
                               data_types = c("L","D","U","acomp","lcomp"),
                               FUN_col="rainbow",
                               t.series.row.trunc = 1
){
  library(plot.matrix)

  t.ser <- rdat$t.series
  t.ser <- t.ser[1:(nrow(t.ser)-t.series.row.trunc),]

  t.ser.nm <- names(t.ser)

  # Landings
  L.nm <- t.ser.nm[grepl("^L..*.ob$",t.ser.nm)]

  # Discards
  D.nm <- t.ser.nm[grepl("^D..*.ob$",t.ser.nm)]

  # Comps
  acomp.nm <- t.ser.nm[grepl("^acomp..*.n$",t.ser.nm)]
  lcomp.nm <- t.ser.nm[grepl("^lcomp..*.n$",t.ser.nm)]

  # Indices
  U.nm <- t.ser.nm[grepl("^U..*.ob$",t.ser.nm)]

  # Combine data
  nm <- as.character(unlist(mget(paste(data_types,"nm",sep="."))))

  t.ser2 <- t.ser[,nm]
  t.ser2[t.ser2==-99999] <- NA
  t.ser2.logical <- !is.na(t.ser2)

  # data set where values are years
  t.ser2.years <- apply(t.ser2.logical,2,function(x){
    out <- rep(NA,length(x))
    out[which(x)] <- as.numeric(rownames(t.ser2.logical)[x])
    out
  })
  rownames(t.ser2.years) <- rownames(t.ser2)

  fleetAbb <-    gsub("(L|D|U|acomp|lcomp)(.)(.*)(.)(ob|n)$","\\3",colnames(t.ser2.years))
  dataTypeAbb <- gsub("(L|D|U|acomp|lcomp)(.)(.*)(.)(ob|n)$","\\1",colnames(t.ser2.years))

  fleetName <- key_fleet_abb[fleetAbb]
  dataTypeName <- key_data_type_abb[dataTypeAbb]

  colnames(t.ser2.years) <- paste(dataTypeName,fleetName)

  # data set where values indicate cell color
  t.ser2.fleet <- t.ser2.years
  for(i in 1:ncol(t.ser2.fleet)){
    x <- t.ser2.fleet[,i]
    x[which(!is.na(x))] <- fleetName[i]
    t.ser2.fleet[,i] <- x
  }

  # Plot using plot.matrix package
  plot(t(t.ser2.fleet),xlab="",ylab="",col=get(FUN_col),main="")
}
