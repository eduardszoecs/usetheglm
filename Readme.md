# Ecotoxicology is not normal.
## - How the use of proper statistical models can increase statistical power in ecotoxicological experiments.

## Online repository for the paper submitted to [Environmental Science and Pollution Research](http://www.springer.com/environment/journal/11356)

### Authors: 

[Eduard Szöcs](http://www.uni-koblenz-landau.de/en/campus-landau/faculty7/environmental-sciences/landscape-ecology/Staff/eduardszoecs), [Ralf B. Schäfer](http://www.uni-koblenz-landau.de/en/campus-landau/faculty7/environmental-sciences/landscape-ecology/Staff/ralf-schaefer/ralf-schaefer)


### Structure of this repository

* `~/cache/`	: cached files (simulation results)
* `~/report/` 	: Manuscript files (LaTeX)
* `~/src/`    	: Source code (R)
* `~/supplement/`: Supplemental files (LaTeX)



#### Structure of `~/src/`

R code has been written hierarchically, so script must be run in a specific order.

* `~/0-load.R`   				: Defines the project structure, loads packages, sources functions and sets some switches. This scrip has to be run first!
* `~/1-simulation_counts.R`    	: Code to run count data simulations
* `~/1-simulation_survival.R`	: Code to run binomial data simulations
* `~/2-results.R`				: grabs the results from previous code and makes the graphics
* `~/3-brock.R`					: Code for motivating example (counts)
* `~/3-weber.R`					: Code for motivating example (binomial)
* `functions.R`					: Custom functions needed for the simulations


### How to reproduce the results

* Download and extract this repository as ZIP.
* Extract the repository
* Please uncomment and point the path in `~/0-load.R`, Lines 8-11 to the extracted folder.
* run `~/0-load.R`
* if necessary install missing packages
* then you can run the other R scripts.
* If the switches `sim1` and `sim2` in `~/0-load.R` are set to `FALSE`, the simulations results are grabbed from cache. Otherwise, the actual simulation will run (takes some hours...).


```{r}
> sessionInfo()
R version 3.1.2 (2014-10-31)
Platform: x86_64-pc-linux-gnu (64-bit)

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB.UTF-8       
 [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
[10] LC_TELEPHONE=C             LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] splines   grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] exactRankTests_0.8-27 multcomp_1.3-8        TH.data_1.0-6         survival_2.37-7      
 [5] mvtnorm_1.0-2         MASS_7.3-37           gridExtra_0.9.1       plyr_1.8.1           
 [9] ggplot2_1.0.0         reshape2_1.4.1       

loaded via a namespace (and not attached):
 [1] colorspace_1.2-4 digest_0.6.8     gtable_0.1.2     lattice_0.20-29  munsell_0.4.2    proto_0.3-10    
 [7] Rcpp_0.11.3      sandwich_2.3-2   scales_0.2.4     stringr_0.6.2    tools_3.1.2      zoo_1.7-11 
 ```
