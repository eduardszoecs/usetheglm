#####--------------------------------------------------------------------------
### Setup project structure
rm(list =ls()) 

# Project Path
## You have to change this!
prj <- "/home/edisz/Documents/Uni/Projects/PHD/6USETHEGLM/"

# Subfolder paths
srcdir <- file.path(prj, "src")     # source code
datadir <- file.path(prj, "data")   # data
cachedir <- file.path(prj, "cache")
figdir <- file.path(prj, "report/fig") # figures for latex

#####--------------------------------------------------------------------------
### install missing if needed!
require(reshape2)
require(ggplot2)
require(plyr)
require(coefplot2)
require(gridExtra)
require(bbmle)
require(lmtest)
require(MASS)
require(gridExtra)

#####--------------------------------------------------------------------------
### Source defined functions
# source(file.path(srcdir, "themes.R"))
source(file.path(srcdir, "functions.R"))

#####--------------------------------------------------------------------------
### Others
# check if load.R already run
ld <- TRUE
# Run simulations? 
sim1 <- FALSE
sim2 <- FALSE
# keep these objects
keep_obj <- c(ls(), 'keep_obj')
# seed
seed <- 1606
# export plots?
exp_plot <- TRUE
# number of simulated datasets
nsims <- 100
