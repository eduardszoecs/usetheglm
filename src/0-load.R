rm(list = ls()[!ls() %in% 'prj'])
#####--------------------------------------------------------------------------
### Setup project structure

# Project Path
## You have to change this!
# prj <- "/home/edisz/Documents/Uni/Projects/PHD/6USETHEGLM/"
# prj <- 'C:\\Users\\Edi\\Documents\\usetheglm'
if(!exists("prj")){
  stop("You need to create a object 'prj' that points to the top folder, 
       e.g. prj <- '/home/edisz/Documents/Uni/Projects/PHD/4BFG/Project'!")
}


# Subfolder paths
srcdir <- file.path(prj, "src")     # source code
datadir <- file.path(prj, "data")   # data
cachedir <- file.path(prj, "cache")
figdir <- file.path(prj, "report", "fig") # figures for latex

#####--------------------------------------------------------------------------
### install missing if needed!
require(reshape2)
require(ggplot2)
require(plyr)
require(gridExtra)
require(bbmle)
require(lmtest)
require(MASS)
# require(pbkrtest)
# require(multcomp)


#####--------------------------------------------------------------------------
### Source defined functions
# source(file.path(srcdir, "themes.R"))
source(file.path(srcdir, "functions.R"))

#####--------------------------------------------------------------------------
### Others
# check if load.R already run
ld <- TRUE
# Run simulations? 
sim1 <- TRUE
sim2 <- FALSE
sim3 <- FALSE
# keep these objects
keep_obj <- c(ls(), 'keep_obj')
# seed
# seed <- 1234
# export plots?
exp_plot <- TRUE
# number of simulated datasets
nsims <- 100
