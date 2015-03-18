Ecotoxicology is not normal.
============================
## How the use of proper statistical models can increase statistical power in ecotoxicological experiments.

Online repository for the paper submitted to [Environmental Science and Pollution Research](http://www.springer.com/environment/journal/11356)

### Authors
[Eduard Szöcs](http://www.uni-koblenz-landau.de/en/campus-landau/faculty7/environmental-sciences/landscape-ecology/Staff/eduardszoecs), [Ralf B. Schäfer](http://www.uni-koblenz-landau.de/en/campus-landau/faculty7/environmental-sciences/landscape-ecology/Staff/ralf-schaefer/ralf-schaefer)


### Structure of this repository

* `~/cache/`	: cached files (simulation results)
* `~/manuscript/`  : manuscript files (LaTeX) - including review process
* `~/src/`    	: Source code (R)
* `~/supplement/`: Supplementary material (LaTeX)

#### Structure of `~/src/`

R code has been written hierarchically, that is scripts must be run in a specific order.

* `~/0-load.R`   				: Defines the project structure, loads packages, sources functions and sets some switches. This script has to be run first!
* `~/0-functions.R`					: Custom functions needed for the simulations
* `~/1-simulations.R`   : Run simulations
* `~/2-results.R`				: Compile results
* `~/3-brock.R`					: Script for motivating example (counts)
* `~/3-weber.R`					: Script for motivating example (binomial)


### How to reproduce the results

1. Download ('Download ZIP') this repository.
2. Extract the zip file.
3. Uncomment and point the path in `~/0-load.R`, Lines 8-11 to the extracted folder.
4. source `~/0-load.R` to setup the project
5. source `~/1-simulations.R` to run the simulations
6. source `~/2-results.R`	to compile the results

### Notes

We added a snapshot of all used packages using [packrat](http://rstudio.github.io/packrat/) to this repo.
packrat will try to install these packages from the snapshot on first start of R from this repository.


`~/0-load.R` contains some switches (at the end of the file):

* If the switches `sim1` and `sim2` are set to `FALSE`, the simulations results are grabbed from cache. Otherwise, the actual simulation will run.
* the count data simulations took ~ 12 hours on a amazon EC2 c3.2xlarge instance.
* simulations can be run in parallel, set `parallel = TRUE` to use the `parallel` package for parallelization.



### Information about the R Environment used

```{r}
> sessionInfo()
R version 3.1.3 (2015-03-09)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 14.04.2 LTS

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_GB.UTF-8   
 [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_GB.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] survival_2.38-1 MASS_7.3-39     ggplot2_1.0.0   reshape2_1.4.1 

loaded via a namespace (and not attached):
 [1] codetools_0.2-11 colorspace_1.2-6 digest_0.6.8     gtable_0.1.2     htmltools_0.2.6  lattice_0.20-30  munsell_0.4.2    packrat_0.4.3    plyr_1.8.1
[10] proto_0.3-10     Rcpp_0.11.5      rmarkdown_0.5.1  sandwich_2.3-2   scales_0.2.4     splines_3.1.3    stringr_0.6.2    tools_3.1.3      zoo_1.7-11
 ```
