# Use the GLM, Luke
### How the use of proper statistical models can increase statistical power in ecotoxicological experiments.


#### Structure of this repository

* `~/src`    : R-Code
* `~/report` : Manuscript (LaTeX code)



##### R-Code

* `~/0-load.R`   : Defines the project structure, loads packages, sources functions and sets some switches. This scrip has to be run first!
* `~/1-simulation_counts.R`    :  Code to run count data simulations
* `~/1-simulation_survival.R`	: Code to run binomial data simulations
* `~2-results.R`	: grabs the results from previous code and makes the graphics
* `functions.R`		: Custom functions defined