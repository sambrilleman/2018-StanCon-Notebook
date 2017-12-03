This repository contains the source files required to compile the following RMarkdown notebook, which was submitted to [StanCon 2018](http://mc-stan.org/events/stancon2018/):

&nbsp;&nbsp;&nbsp;&nbsp;Joint longitudinal and time-to-event models via Stan** by Sam Brilleman, Michael Crowther, Margarita Moreno-Betancur, Jacqueline Buros Novik, and Rory Wolfe

A compiled PDF version of the notebook is included in the repository's root directory. To compile the notebook yourself the following R packages will need to be installed (including their dependencies): knitr, ggplot2, rstan, survival, data.table, rstanarm. 

**Note:** At the time of writing this, the version of **rstanarm** that is required is the development version available [here](https://github.com/stan-dev/rstanarm) since this contains the `stan_jm` modelling function, however, this will be uploaded to CRAN in the near future, so it may be worth checking [here](https://cran.r-project.org/web/packages/rstanarm/index.html) whether **rstanarm** version >=2.17.0 has been released on CRAN. If it has, then it would be easier to install the package directly from within R by calling `install.packages("rstanarm")` instead of downloading the development version from GitHub.
