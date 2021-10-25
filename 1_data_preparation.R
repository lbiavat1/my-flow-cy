rm(list = ls())

### load libraries

library(Spectre)
Spectre::package.check()
# load all required packages
Spectre::package.load()

library(tidyverse)
library(purrr)

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
if(!dir.exists("output"))
  dir.create("output", showWarnings = FALSE)
setwd("output")
OutputDirectory <- getwd()
setwd(PrimaryDirectory)


# prep and read in metadata

Filename <- list.files(InputDirectory)
# # write_csv(as.data.frame(Filename), file = paste(MetaDirectory, "sample.details.csv", sep = "/"))
# metadata_file <- paste(MetaDirectory, "sample.details.csv", sep = "/")
# sample_details <- read_csv(metadata_file)
# #add cells per sample
# sample_details$`Cells per sample` <- map_dbl(data.list, nrow)
# sample_details
# write_csv(sample_details, file = metadata_file)

metadata_file <- paste(MetaDirectory, "sample.details.csv", sep = "/")
sample_details <- read_csv(metadata_file)
sample_details
# keep only baseline sample
samples_to_keep <- sample_details %>% dplyr::filter(Timepoint == "Baseline")

# import data
list.files(InputDirectory)
data.list <- Spectre::read.files(file.loc = InputDirectory, file.type = ".csv", 
                                 do.embed.file.names = TRUE, header = TRUE)
check <- do.list.summary(data.list)
check$name.table
# subset only baseline samples
# merge data
cell.data <- do.merge.files(data.list)

# note: remove unwanted samples (limit analysis to baseline/exp vs non-exp)
# select only samples from Baseline, exp vs not_exp - remove HD
file.names <- gsub(".csv", "", samples_to_keep$Filename)
data <- cell.data %>% dplyr::filter(FileName %in% file.names)

as.matrix(names(cell.data))

# make sample plot
make.colour.plot(do.subsample(cell.data, 20000), "Comp-PE-A :: TCF1", "Comp-BV711-A :: CD127")
  
# add metadata and set parameters
sample.info <- sample_details %>% select("Filename", "Group")

cell.data <- do.add.cols(data, "FileName", sample.info, "Filename", rmv.ext = TRUE)
class(cell.data)
cell.data

as.matrix(names(cell.data))
# select relevant markers - if needed, subset state vs type markers here
type.markers <- names(cell.data)[c(11, 13:32)]
type.markers

# select markers used for clustering/DR
cluster.markers <- type.markers

# specify sample, group, and batch columns
exp.name <- "J1484 - aMILs expansion"
sample.col <- "Sample"
group.col <- "Group"
batch.col <- "Batch"

data.frame(table(cell.data[[group.col]]))

unique(cell.data[[group.col]])

# Clustering and DR

setwd(OutputDirectory)
if(!dir.exists("output-clustering"))
  dir.create("output-clustering")
clustering.dir <- "output-clustering"
setwd(clustering.dir)

# run flowsom
cell.data <- run.flowsom(cell.data, cluster.markers, meta.k = 12)

# dimrnsionality reduction - DR
unique(cell.data[[group.col]])
subsampling.targets <- c(3000, 3000)
cell.sub <- do.subsample(cell.data, subsampling.targets, group.col)




