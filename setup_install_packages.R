rm(list = ls())

usethis::create_github_token()

gitcreds::gitcreds_set()
gitcreds::gitcreds_get()

if(!require("devtools"))
  install.packages("devtools")

# after installing devtools, restart RStudio

library(devtools)
