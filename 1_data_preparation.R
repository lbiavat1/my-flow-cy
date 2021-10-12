rm(list = ls())

### load libraries

library(Spectre)
Spectre::package.check()
# load all required packages
Spectre::package.load()

library(tidyverse)

### Set PrimaryDirectory
dirname(rstudioapi::getActiveDocumentContext()$path)            # Finds the directory where this script is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Sets the working directory to where the script is located
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory

### Set 'input' directory
setwd(PrimaryDirectory)
if(!dir.exists("data"))
  dir.create("data")
setwd("data/")
InputDirectory <- getwd()
setwd(PrimaryDirectory)

### Set 'metadata' directory
setwd(PrimaryDirectory)
if(!dir.exists("metadata"))
  dir.create("metadata")
setwd("metadata/")
MetaDirectory <- getwd()
setwd(PrimaryDirectory)

### Create output directory
if(!dir.exists("Output"))
  dir.create("Output", showWarnings = FALSE)
setwd("Output")
OutputDirectory <- getwd()
setwd(PrimaryDirectory)

# import data
list.files(InputDirectory)
data.list <- Spectre::read.files(file.loc = InputDirectory, file.type = ".csv", 
                                 do.embed.file.names = TRUE)
check <- do.list.summary(data.list)
check$name.table

# merge data
cell.data <- do.merge.files(data.list)

# prep and read in metadata

Filename <- list.files(InputDirectory)
# write_csv(as.data.frame(Filename), file = paste(MetaDirectory, "sample.details.csv", sep = "/"))
metadata_file <- paste(MetaDirectory, "sample.details.csv", sep = "/")
sample_details <- read_csv(metadata_file)




