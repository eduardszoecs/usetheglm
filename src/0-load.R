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
if (!exists("prj")) {
  stop("You need to create a object 'prj' that points to the top folder, 
       e.g. prj <- '/home/edisz/Documents/Uni/Projects/PHD/4BFG/Project'!")
}

# Subfolder paths
srcdir <- file.path(prj, "src")     # source code
datadir <- file.path(prj, "data")   # data
cachedir <- file.path(prj, "cache")
suppdir <- file.path(prj, "supplement")
figdir <- file.path(prj, "manuscript/revision/report") # figures for latex

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

# or get the packages from packrat
source(file.path(prj, 'packrat/init.R'))


#####--------------------------------------------------------------------------
### Source defined functions
# source(file.path(srcdir, "themes.R"))
source(file.path(srcdir, "0-functions.R"))


#####--------------------------------------------------------------------------
### Switches
# check if load.R already run
ld <- TRUE
# Run simulations? 
# if TRUE the simulation will run (takes some time...), 
# if FALSE simulation results from cache
sim1 <- FALSE
sim2 <- FALSE
# export plots to figure dir?
exp_plot <- TRUE
# number of simulated datasets 
nsims <- 1000
# run sims parallel?
parallel = FALSE
if(parallel){
  require(parallel)
  ncores <- detectCores()
  # keep one core for OS
  if(ncores > 1)
    ncores <- ncores -1
}

# run parametric boostrap?
pb <-  TRUE

# keep these objects
keep_obj <- c(ls(), 'keep_obj')
