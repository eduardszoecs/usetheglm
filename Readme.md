# Ecotoxicology is not normal.
## How the use of proper statistical models can increase statistical power in ecotoxicological experiments.

### Authors: 

* Eduard Szöcs
* Ralf B. Schäfer


### Structure of this repository

* `~/cache/`	: cached files (simulation results)
* `~/report/` 	: Manuscript files (LaTeX)
* `~/src/`    	: Source code (R)
* `~/supplement/: Supplement files (LaTeX)



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
* then you can run the other R scripts.
* If the switches `sim1` and `sim2` in `~/0-load.R` are set to `FALSE`, the simulations results are grabbed from cache. Otherwise, the actual simulation will run (takes some hours...).