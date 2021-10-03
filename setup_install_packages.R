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

BiocManager::install("CATALYST", dependencies = TRUE)

library(CATALYST)
