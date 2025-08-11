rm(list=ls())

# install.packages("devtools")
# install.packages("rlang")
# install.packages("roxygen2")
# install.packages(Rtools)

# You may need to restart R to avoid errors with help file
# .rs.restartR()

library(devtools)
library(pkgdown)
library(rlang)

library(roxygen2)
library(stringr)
# library(Rtools)


# Set working director
# wd <- getwd()

# Create package3
#create("bamExtras") # I think you only need to do this once

# Create vignette (only run to initially create a vignette)
#usethis::use_vignette("intro")
#usethis::use_vignette(name="RefPtCalc",title="Calculating reference points")

# Add dependencies
#usethis::use_package("stringr","Suggests")

# Install package
devtools::document(file.path(dirname(getwd()),"bamExtras")) # Important for updating package
devtools::install(file.path(dirname(getwd()),"bamExtras"),build_vignettes = TRUE)
#devtools::build_manual(pkg=getwd()) # Didn't work when I tried it
