rm(list = ls()[!ls() %in% 'prj'])
### ----------------------------------------------------------------------------
### Setup project structure to run simulations and examples
### Written by Eduard Szöcs
### ----------------------------------------------------------------------------

# Project Paths
## Please uncomment and change accordingly
# eg.
# prj <- "/home/edisz/Documents/Uni/Projects/PHD/6USETHEGLM/" # on Linux
# prj <- 'C:\\Users\\Edi\\Documents\\usetheglm' # on Windows
if(!exists("prj")){
  stop("You need to create a object 'prj' that points to the top folder, 
       e.g. prj <- '/home/edisz/Documents/Uni/Projects/PHD/4BFG/Project'!")
}

# Subfolder paths
srcdir <- file.path(prj, "src")     # source code
datadir <- file.path(prj, "data")   # data
cachedir <- file.path(prj, "cache")
figdir <- file.path(prj, "report", "fig") # figures for latex
markdir <- file.path(prj, "report_markdown", "fig")

#####--------------------------------------------------------------------------
### install missing if needed!
require(reshape2)
require(ggplot2)
require(plyr)
require(gridExtra)
require(MASS)
require(multcomp)
require(exactRankTests)
require(xtable)

#####--------------------------------------------------------------------------
### Source defined functions
# source(file.path(srcdir, "themes.R"))
source(file.path(srcdir, "functions.R"))

#####--------------------------------------------------------------------------
### Other switches
# check if load.R already run
ld <- TRUE
# Run simulations? 
# if TRUE the simulation will run (takes some time...), 
# if FALSE simulation results from cache
sim1 <- FALSE
sim2 <- FALSE
# keep these objects
keep_obj <- c(ls(), 'keep_obj')
# random
# seed <- 1234
# export plots to figure dir?
exp_plot <- TRUE
# number of simulated datasets (may overwritten by simulation scripts)
nsims <- 100
