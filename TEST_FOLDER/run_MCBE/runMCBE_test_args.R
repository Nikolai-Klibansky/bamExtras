rm(list=ls())
#library(bamExtras)

# styr <- rdat_RedPorgy$parms$styr
# endyr <- rdat_RedPorgy$parms$endyr

frm <- formals(run_MCBE)

for(i in seq_along(frm)){
  assign(names(frm)[i],eval(frm[[i]]))
}

CommonName="RedPorgy"
dir_bam_base="RePo_base"
dir_bam_sim="RePo_sim"
