rm(list = ls())

### load libraries

library(Spectre)
Spectre::package.check()
# load all required packages
Spectre::package.load()

library(tidyverse)
library(purrr)
library(CATALYST)
library(diffcyt)
library(tidySingleCellExperiment)
library(flowCore)
library(stringr)
library(scater)
library(scran)


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
# keep only baseline sample
samples_to_keep <- sample_details %>% dplyr::filter(Timepoint == "Baseline") %>%
  dplyr::filter(!grepl("no_harvest", Group))



# prepData for CATALYST - create SCE
CSVfiles <- sample_details$Filename
# keep only baseline sample - AML only!!
samples_to_keep <- sample_details %>% 
  dplyr::filter(!grepl("CL", Disease)) %>%
  dplyr::filter(Timepoint == "Baseline") %>%
  dplyr::filter(!grepl("no_harvest", Group)) %>%
  mutate(pt_id = match(MRN, unique(MRN))) %>%
  mutate(patient_id = paste("PT", pt_id, sep = "_")) %>%
  arrange(pt_id) %>%
  slice(-c(12))
samples_to_keep
as.data.frame(samples_to_keep) %>% select(Filename, Sample, Timepoint, Group, pt_id)
table(samples_to_keep$Group)
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
fs <- read.flowSet(files = fcsFiles, path = fcsDir, truncate_max_range = FALSE)
fs
fs[[1]]@description$`$CYT` <- "FACS"

############ create sample_md & panel_md - CATALYST constructor ################
setwd(MetaDirectory)
sample_metadata <- samples_to_keep %>% select(Filename, Sample, Timepoint, Group, MRN) %>%
  mutate(Filename = gsub(".csv", "", Filename)) %>%
  mutate(file_name = paste0(Filename, ".fcs")) %>%
  mutate(pt_id = match(MRN, unique(MRN))) %>%
  mutate(patient_id = paste("PT", pt_id, sep = "_")) %>%
  mutate(condition = Group) %>%
  mutate(sample_id = paste(patient_id, Group, sep = "_")) %>%
  select(file_name, Sample, Timepoint, Group, MRN, patient_id, condition, sample_id) %>%
  write_csv("sample_metadata_AMLMMExpVsNExp.csv")
sample_metadata
deID_key <- sample_metadata %>%
  select(MRN, patient_id) %>%
  distinct() %>%
  write_csv("deID_key_AMLMMExpVsNExp.csv")

sample_md <- sample_metadata %>%
  select(file_name, patient_id, condition, sample_id) %>%
  as.data.frame()
sample_md




# # create tibble sample_md
# sample_md <- samples_to_keep %>% select(Filename, Sample, Timepoint, Group) %>%
#   mutate(Filename = gsub(".csv", "", Filename)) %>%
#   mutate(file_name = paste0(Filename, ".fcs")) %>%
#   mutate(patient_id = Sample) %>%
#   mutate(condition = Group) %>%
#   mutate(sample_id = paste(Sample, Timepoint, sep = "_")) %>%
#   select(file_name, patient_id, condition, sample_id) %>%
#   as.data.frame()

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
setwd(OutputDirectory)
sce <- prepData(fs, panel_md, sample_md, FACS = TRUE)
assay(sce, "exprs") <- assay(sce, "counts")

n_events <- min(n_cells(sce))
n_events
n_cells(sce)
plotCounts(sce, group_by = "sample_id", color_by = "condition")


plotNRS(sce, features = type_markers(sce), color_by = "condition")
pbMDS(sce, features = type_markers(sce), color_by = "condition")

# TODO: try Spectre, milo, consider CFA 

########### run std analysis workflow ########################################
std_sce <- sce

std_sce <- cluster(std_sce, features = "type", xdim = 10, ydim = 10, maxK = 20, verbose = TRUE,
                   seed = 1)
std_sce <- runDR(std_sce, dr = "UMAP", cells = 2000, features = "type", assay = "exprs")
delta_area(std_sce)
plotAbundances(std_sce, k = "meta8", by = "cluster_id", group_by = "condition")
plotDR(std_sce, dr = "UMAP", color_by = "meta8", facet_by = "condition") +
  geom_density2d(binwidth = 0.006, colour = "black")

annotation_table <- as.data.frame(cbind(c(1:8), c(1:8)))
colnames(annotation_table) <- c("meta8", "Clusters")
annotation_table$Clusters <- factor(annotation_table$Clusters)
std_sce <- mergeClusters(std_sce, k = "meta8", 
                     table = annotation_table, id = "cluster_annotation", overwrite = TRUE)
std_sce$cluster_annotation <- cluster_ids(std_sce, "cluster_annotation")

FDR_cutoff <- 0.05
ei <- std_sce@metadata$experiment_info
plotAbundances(std_sce, k = "cluster_annotation", by = "cluster_id", group_by = "condition")
plotExprHeatmap(std_sce, features = type_markers(std_sce), k = "cluster_annotation", 
                by = "cluster_id", scale = "last", bars = TRUE, perc = TRUE)

# DA using edgeR
design <- createDesignMatrix(ei,
                             cols_design = c("condition"))
contrast <- createContrast(c(0, 1))

nrow(contrast) == ncol(design)

out_DA <- diffcyt(std_sce,
                  experiment_info = ei, design = design, contrast = contrast,
                  analysis_type = "DA", method_DA = "diffcyt-DA-edgeR",
                  clustering_to_use = "cluster_annotation", verbose = TRUE, subsampling = TRUE,
                  transform = FALSE, normalize = FALSE, min_cells = 3, min_samples = 1)

da <- rowData(out_DA$res)
plotDiffHeatmap(std_sce, da, top_n = 8, all = TRUE, fdr = FDR_cutoff)
CATALYST::pbMDS(std_sce, color_by = "condition", features = type_markers(sce), fun = "median")

getwd()
setwd(OutputDirectory)
getwd()
saveRDS(std_sce, file = "std_sceMMAML.rds")

######### prep to run miloR - subsample SCE per sample_id ######################

# simple function: x -SCE, n_cells -#cells per sample
subsampleSCE <- function(x, n_cells){
  cs <- split(seq_len(ncol(x)), x$sample_id)
  cs <- unlist(lapply(cs, function(.) sample(., min(n_cells, length(.)))))
  x <- x[, cs]
  return(x)
}
subsampleSCE_condition <- function(x, n_cells){
  cs <- split(seq_len(ncol(x)), x$condition)
  cs <- unlist(lapply(cs, function(.) sample(., min(n_cells, length(.)))))
  x <- x[, cs]
  return(x)
}

cells <- min(c(3000, n_events))
cells <- 20000
sub.sce <- subsampleSCE_condition(sce, cells)


# std.sce <- sub.sce
# seed <- 123456
# set.seed(seed)
# std.sce <- cluster(std.sce, features = "type", xdim = 10, ydim = 10, maxK = 20, 
#                verbose = TRUE, seed = seed)
# delta_area(std.sce)
# # Run dimensionality reduction
# 
# std.sce <- runDR(std.sce, dr =  "UMAP", features = "type")
# 
# plotAbundances(std.sce, k = "meta8", by = "cluster_id", group_by = "condition")
# plotDR(std.sce, dr = "UMAP", color_by = "meta8", facet_by = "condition") +
#   geom_density2d(binwidth = 0.006, colour = "black")
# 
# plotExprHeatmap(std.sce, features = type_markers(sce), k = "meta8", 
#                 by = "cluster_id", scale = "last", bars = TRUE, perc = TRUE)
# 
# 
logcounts(sub.sce) <- log(counts(sub.sce) + 1)

sub.sce <- runPCA(sub.sce, ncomponents = 15)

percent.var <- attr(reducedDim(sub.sce), "percentVar")
plot(percent.var, log = "y", xlab = "PC", ylab = "Variance explained (%)")
plotReducedDim(sub.sce, colour_by = "condition", dimred = "PCA")

sub.sce <- runUMAP(sub.sce, dimred = "PCA", name = "umap")

plotReducedDim(sub.sce, colour_by = "condition", dimred = "umap")
plotReducedDim(sub.sce, colour_by = "condition", dimred = "PCA", ncomponents = 3)

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
sce_design <- data.frame(colData(sce_milo))[, c("sample_id", "condition")]
## Convert batch info from integer to factor
sce_design$condition <- as.factor(sce_design$condition)
sce_design <- distinct(sce_design)
rownames(sce_design) <- sce_design$sample_id
sce_design

# computing nhood connectivity
sce_milo <- calcNhoodDistance(sce_milo, d = 15, reduced.dim = "PCA")
sce_milo

# testing
da_results <- testNhoods(sce_milo, design = ~condition, 
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
ggsave("miloPlot_AMLMM.pdf")


# assign a "condition" label to each nhood by finding the most abundant "condition" within cells in each neighbourhood. 
da_results <- annotateNhoods(sce_milo, da_results, coldata_col = "condition")
head(da_results)
ggplot(da_results, aes(condition_fraction)) +
  geom_histogram(bins=50)

# define a threshold for condition_fraction to exclude nhoods that are a mix of cell types.
da_results$condition <- ifelse(da_results$condition_fraction < 0.7, "Mixed", da_results$condition)

plotDAbeeswarm(da_results, group.by = "condition")

# add log-normalized counts to milo object
sce_milo <- logNormCounts(sce_milo)

da_results <- groupNhoods(sce_milo, da.res = da_results, da.fdr = 0.1, max.lfc.delta = 2)
as_tibble(da_results)


type.markers <- rownames(sce_milo)[c(11, 13:32)]
cellID <- sce_milo %>% select(cell)
colnames(sce_milo) <- c(1:69000)
da_nhood_markers <- findNhoodGroupMarkers(sce_milo, da_results, subset.row = type.markers)

plotNhoodGroups(sce_milo, da_results, layout = "umap")
plotDAbeeswarm(da_results, group.by = "NhoodGroup")

plotDAbeeswarm(groupNhoods(sce_milo, da_results, max.lfc.delta = 1),
               group.by = "NhoodGroup") + 
  ggtitle("max LFC delta = 1")

plotDAbeeswarm(groupNhoods(sce_milo, da_results, da.fdr = 0.1, max.lfc.delta = 2, overlap = 1), 
               group.by = "NhoodGroup") + ggtitle("overlap = 1")

## Exclude zero counts genes
keep.rows <- rowSums(logcounts(sce_milo)) != 0
sce_milo <- sce_milo[keep.rows, ]

## Find HVGs
dec <- modelGeneVar(sce_milo, subset.row = type.markers)
hvgs <- getTopHVGs(dec, n = NULL)
head(hvgs)

nhood_markers <- findNhoodGroupMarkers(sce_milo, da_results, subset.row = hvgs)

head(nhood_markers)

gr11_markers <- nhood_markers[c("logFC_11", "adj.P.Val_11")] 
colnames(gr11_markers) <- c("logFC", "adj.P.Val")

head(gr11_markers[order(gr11_markers$adj.P.Val), ])

markers <- rownames(nhood_markers)[nhood_markers$adj.P.Val_11 < 0.01 & nhood_markers$logFC_11 > 0]

plotNhoodExpressionGroups(sce_milo, da_results, features = markers,
                          subset.nhoods = da_results$NhoodGroup %in% c('11'), scale=TRUE,
                          grid.space = "fixed")

plotNhoodExpressionGroups(sce_milo, da_results, features = markers, scale=TRUE,
                          subset.nhoods = NULL, cluster_features = TRUE,
                          grid.space = "fixed")



############# Spectre ########################################################
# import data
list.files(InputDirectory)
samples_to_keep$Filename

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
# make.colour.plot(do.subsample(cell.data, 20000), "Comp-PE-A :: TCF1", "Comp-BV711-A :: CD127")

# add metadata and set parameters
sample.info <- samples_to_keep %>% select("Filename", "Group")

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
cell.data

# dimrnsionality reduction - DR
unique(cell.data[[group.col]])
subsampling.targets <- c(20000, 20000)
cell.sub <- do.subsample(cell.data, subsampling.targets, group.col)

cell.sub <- run.umap(cell.sub, cluster.markers)
cell.sub

make.colour.plot(cell.sub, x.axis = "UMAP_X", y.axis = "UMAP_Y", col.axis = "FlowSOM_metacluster", 
                 col.type = "factor", add.label = TRUE)

make.multi.plot(cell.sub, "UMAP_X", "UMAP_Y", "FlowSOM_metacluster", col.type = 'factor',
                divide.by = group.col, add.density = TRUE)

RColorBrewer::display.brewer.all()

ggplot(cell.sub, aes(x = UMAP_X, y = UMAP_Y, col = as.factor(FlowSOM_metacluster))) +
  geom_point(size = 0.1) +
  facet_wrap(~ Group) +
  theme_bw()

exp <- do.aggregate(cell.data, type.markers, by = "FlowSOM_metacluster")
make.pheatmap(exp, "FlowSOM_metacluster", cluster.markers)

setwd(OutputDirectory)
dir.create("Output - annotation")
setwd("Output - annotation")

# annotation 


setwd(OutputDirectory)
dir.create("Output - summary data")
setwd("Output - summary data")

# stats

variance.test <- 'kruskal.test'
pairwise.test <- "wilcox.test"

as.matrix(unique(cell.data[[group.col]]))

comparisons <- list(c("exp", "not_exp"))
comparisons

grp.order <- c("exp", "not_exp")
grp.order

as.matrix(type.markers)

dyn.cols <- type.markers
dyn.cols

meta.data <- sample_details %>% filter(Filename %in% samples_to_keep$Filename)
counts <- meta.data[, c(3,9)]

sum.dat <- create.sumtable(dat = cell.data,
                           sample.col = sample.col,
                           pop.col = "FlowSOM_metacluster",
                           use.cols = dyn.cols,
                           annot.cols = c(group.col),
                           counts = NULL)
sum.dat









