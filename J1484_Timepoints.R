rm(list = ls())

### load libraries

library(Spectre)
Spectre::package.check()
# load all required packages
Spectre::package.load()

library(tidyverse)
library(purrr)
library(CATALYST)
library(tidySingleCellExperiment)
library(flowCore)
library(stringr)
library(scater)

library(miloR)
library(patchwork)

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
if(!dir.exists("fcs"))
  dir.create("fcs")
setwd("fcs/")
fcsDir <- getwd()
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
# keep only first infusion pre/post pairs
samples_to_keep <- sample_details %>% dplyr::filter(Group == "exp")
samples_to_keep <- samples_to_keep %>% arrange(MRN) %>% add_count(MRN) %>% 
  dplyr::filter(n > 1) %>% 
  arrange(desc(n)) %>%
  slice(- c(3, 6, 9 , 10, 13, 17)) %>% 
  select(-n) %>% add_count(MRN) %>%
  mutate(pair = (n == 2))

all(samples_to_keep$pair == TRUE)


# prepData for CATALYST - create SCE
CSVfiles <- samples_to_keep$Filename

######################### convert csv files to fcs #############################
csvTofcs <- function(file.names, dest){
  # create an empty list to start
  DataList <- list()
  
  for(file in file.names){
    tmp <- read_csv(file.path(file))
    file <- gsub(".csv", "", file)
    DataList[[file]] <- tmp
  }
  rm(tmp)
  
  filenames <- names(DataList)
  head(DataList)
  
  # convert csv to fcs
  
  for(i in c(1:length(filenames))){
    data_subset <- DataList[i]
    data_subset <- data.table::rbindlist(as.list(data_subset))
    file_name <- names(DataList)[i]
    
    metadata <- data.frame(name = dimnames(data_subset)[[2]], desc = "")
    
    # create FCS file metadata
    # metadata$range <- apply(apply(data_subset, 2, range), 2, diff)
    metadata$minRange <- apply(data_subset, 2, min)
    metadata$maxRange <- apply(data_subset, 2, max)
    
    
    # data as matrix by exprs
    data_subset.ff <- new("flowFrame", exprs = as.matrix(data_subset),
                          parameters = AnnotatedDataFrame(metadata))
    
    head(data_subset.ff)
    write.FCS(data_subset.ff, paste0(dest, "/", file_name, ".fcs"), what = "numeric")
  }
}
setwd(InputDirectory)
# csvTofcs(CSVfiles, fcsDir)

# read fcs files as flowSet and add $CYT keyword

fcsFiles <- list.files(path = fcsDir, pattern = ".fcs")
length(fcsFiles)
fcs_to_keep <- sapply(samples_to_keep$Filename, simplify = TRUE, function(.) gsub(".csv", ".fcs", .))
fcsFiles <- fcsFiles[fcsFiles %in% fcs_to_keep]
fcsFiles
fs <- read.flowSet(files = fcsFiles, path = fcsDir, truncate_max_range = FALSE)
fs
fs[[1]]@description$`$CYT` <- "FACS"
fs[[1]]@description$`$CYT`

############ create sample_md & panel_md - CATALYST constructor ################
# create tibble sample_md - assign real patient_id to pair samples!
sample_metadata <- samples_to_keep %>% select(Filename, Sample, Timepoint, Group, MRN) %>%
  mutate(Filename = gsub(".csv", "", Filename)) %>%
  mutate(file_name = paste0(Filename, ".fcs")) %>%
  mutate(pt_id = match(MRN, unique(MRN))) %>%
  mutate(patient_id = paste("PT", pt_id, sep = "_")) %>%
  mutate(condition = Timepoint) %>%
  mutate(sample_id = paste(patient_id, Timepoint, sep = "_")) %>%
  select(file_name, Sample, Timepoint, Group, MRN, patient_id, condition, sample_id) %>%
  write_csv("sample_metadata.csv")
sample_metadata
deID_key <- sample_metadata %>%
  select(MRN, patient_id) %>%
  distinct() %>%
  write_csv("deID_key.csv")

sample_md <- sample_metadata %>%
  select(file_name, patient_id, condition, sample_id) %>%
  as.data.frame()
sample_md

# create tibble panel_md
fcs_colname <- colnames(fs)
fcs_colname
antigen <- fcs_colname
antigen[10:35] <- sapply(fcs_colname[10:35], function(.) unlist(str_split(., " :: "))[2])
fluorochrome <- fcs_colname
fluorochrome[10:35] <- sapply(fcs_colname[10:35], function(.) unlist(str_split(., " :: "))[1])
marker_class <- fcs_colname
marker_class[ c(1:10, 12, 33:36)] <- "none"
marker_class[c(11, 13:32)] <- "type"
panel_md <- as_tibble(cbind(fcs_colname, antigen, fluorochrome, marker_class))
as.data.frame(panel_md)

######### prepData - create SCE using CATALYST #################################
sce <- prepData(fs, panel_md, sample_md, FACS = TRUE)
assay(sce, "exprs") <- assay(sce, "counts")

# p <- plotExprs(sce, features = NULL, color_by = "condition")
# p$facet$params$ncol <- 9
# p

n_events <- min(n_cells(sce))
n_events
n_cells(sce)
plotCounts(sce, group_by = "sample_id", color_by = "condition")

plotNRS(sce, features = type_markers(sce), color_by = "condition")

std_sce <- sce
std_sce <- cluster(std_sce, features = "type", xdim = 10, ydim = 10, maxK = 20, verbose = TRUE)
delta_area(std_sce)
std_sce <- runDR(std_sce, dr = "UMAP", cells = 2000, features = "type", assay = "exprs")
plotAbundances(std_sce, k = "meta8", by = "cluster_id", group_by = "condition")


######### prep to run miloR - subsample SCE per sample_id ######################

# simple function: x -SCE, n_cells -#cells per sample
subsampleSCE <- function(x, n_cells){
  cs <- split(seq_len(ncol(x)), x$sample_id)
  cs <- unlist(lapply(cs, function(.) sample(., min(n_cells, length(.)))))
  x <- x[, cs]
  return(x)
}

cells <- min(c(2000, n_events))
sub.sce <- subsampleSCE(sce, cells)
logcounts(sub.sce) <- log(counts(sub.sce) + 1)

sub.sce <- runPCA(sub.sce, ncomponents = 15)
sub.sce <- runUMAP(sub.sce, dimred = "PCA", name = "umap")

plotReducedDim(sub.sce, colour_by = "condition", dimred = "umap")

sce_milo <- Milo(sub.sce)
sce_milo

# construct kNN graph
sce_milo <- buildGraph(sce_milo, k = 30, d = 15, reduced.dim = "PCA")
sce_milo

# defining representative nhoods on the kNN graph
sce_milo <- makeNhoods(sce_milo, prop = 0.1, k = 30, d = 15, refined = TRUE, 
                       reduced_dims = "PCA")
plotNhoodSizeHist(sce_milo)

# counting cells in nhoods
sce_milo <- countCells(sce_milo, meta.data = data.frame(colData(sce_milo)), sample = "sample_id")
head(nhoodCounts(sce_milo))

# defining experimental design
sce_design <- data.frame(colData(sce_milo))[,c("sample_id", "condition", "patient_id")]
## Convert batch info from integer to factor
sce_design$condition <- as.factor(sce_design$condition)
sce_design <- distinct(sce_design)
rownames(sce_design) <- sce_design$sample_id
sce_design

# computing nhood connectivity
sce_milo <- calcNhoodDistance(sce_milo, d = 15, reduced.dim = "PCA")
sce_milo

# testing
da_results <- testNhoods(sce_milo, design = ~patient_id + condition, 
                         design.df = sce_design)
da_results %>%
  arrange(SpatialFDR) %>%
  head()

ggplot(da_results, aes(PValue)) +
  geom_histogram(bins = 50)

ggplot(da_results, aes(logFC, -log10(SpatialFDR))) +
  geom_point() +
  geom_hline(yintercept = 1)

sce_milo <- buildNhoodGraph(sce_milo)

umap_plot <- plotReducedDim(sce_milo, dimred = "umap", 
                            colour_by = "condition", text_by = "condition", 
                            text_size = 3) +
  guides(fill = "none")

nh_graph_plot <- plotNhoodGraphDA(sce_milo, da_results, layout = "umap", alpha = 0.05)

umap_plot + nh_graph_plot +
  plot_layout(guides = "collect")
setwd(OutputDirectory)
ggsave("miloPlot.pdf")

da_results <- annotateNhoods(sce_milo, da_results, coldata_col = "condition")
head(da_results)

da_results$condition <- ifelse(da_results$condition_fraction < 0.7, "Mixed", da_results$condition)

plotDAbeeswarm(da_results, group.by = "condition")

sce_milo <- logNormCounts(sce_milo)

# da_results$NhoodGroup <- as.numeric(da_results$SpatialFDR < 0.1 & da_results$logFC < 0)
# as_tibble(da_results)
# da_nhood_markers <- findNhoodGroupMarkers(sce_milo, da_results, 
#                                           subset.row = rownames(sce_milo)[c(11,13:32)],
#                                           aggregate.samples = TRUE, sample_col = "sample_id")

## Run buildNhoodGraph to store nhood adjacency matrix
sce_milo <- buildNhoodGraph(sce_milo)

## Find groups
da_results <- groupNhoods(sce_milo, da_results, max.lfc.delta = 2)
head(da_results)

plotNhoodGroups(sce_milo, da_results, layout="umap")
plotDAbeeswarm(da_results, "NhoodGroup")


