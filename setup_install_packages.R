rm(list = ls())

usethis::create_github_token()

gitcreds::gitcreds_set()
gitcreds::gitcreds_get()

if(!require("devtools"))
  install.packages("devtools")

# after installing devtools, restart RStudio

library(devtools)

options(timeout = 6000)
# download and install spectre 
devtools::install_github("immunedynamics/spectre")

library(Spectre)
Spectre::package.check()

if(!require("tidyverse"))
  install.packages("tidyverse")

library(tidyverse)

install.packages("aws.signature", repos = c(cloudyr = "http://cloudyr.github.io/drat", getOption("repos")))

BiocManager::install("CATALYST", dependencies = TRUE)

library(CATALYST)

BiocManager::install("tidySingleCellExperiment")
library(tidySingleCellExperiment)
