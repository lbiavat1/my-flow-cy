rm(list = ls())

library(miloR)
library(SingleCellExperiment)
library(scater)
library(tidyverse)
library(patchwork)

### Set PrimaryDirectory
dirname(rstudioapi::getActiveDocumentContext()$path)            # Finds the directory where this script is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Sets the working directory to where the script is located
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory

# load data
library(MouseGastrulationData)

select_samples <- c(2,  3,  6, 15,
                    # 4, 19, 
                    10, 14, 20, 30
                    #31, 32
)
sce <- EmbryoAtlasData(samples = select_samples)
sce

sce <- sce[, apply(reducedDim(sce, "pca.corrected"), 1, function(.) !all(is.na(.)))]
sce <- runUMAP(sce, dimred = "pca.corrected", name = "umap")

plotReducedDim(sce, dimred = "umap", colour_by = "stage")

sce_milo <- Milo(sce)
sce_milo

# construct kNN graph
sce_milo <- buildGraph(sce_milo, k = 30, d = 30, reduced.dim = "pca.corrected")
sce_milo

# defining representative nhoods on the kNN graph
sce_milo <- makeNhoods(sce_milo, prop = 0.1, k = 30, d = 30, refined = TRUE, 
                       reduced_dims = "pca.corrected")
plotNhoodSizeHist(sce_milo)

# counting cells in nhoods
sce_milo <- countCells(sce_milo, meta.data = data.frame(colData(sce_milo)), sample = "sample")
head(nhoodCounts(sce_milo))

# defining experimental design
sce_design <- data.frame(colData(sce_milo))[,c("sample", "stage", "sequencing.batch")]
## Convert batch info from integer to factor
sce_design$sequencing.batch <- as.factor(sce_design$sequencing.batch)
sce_design <- distinct(sce_design)
rownames(sce_design) <- sce_design$sample
sce_design

# computing nhood connectivity
sce_milo <- calcNhoodDistance(sce_milo, d = 30, reduced.dim = "pca.corrected")
sce_milo

# testing
da_results <- testNhoods(sce_milo, design = ~ sequencing.batch + stage, 
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
                            colour_by = "celltype", text_by = "celltype", 
                            text_size = 3) +
  guides(fill = "none")

nh_graph_plot <- plotNhoodGraphDA(sce_milo, da_results, layout = "umap", alpha = 0.05)

umap_plot + nh_graph_plot +
  plot_layout(guides = "collect")

da_results <- annotateNhoods(sce_milo, da_results, coldata_col = "celltype")
head(da_results)

da_results$celltype <- ifelse(da_results$celltype_fraction < 0.7, "Mixed", da_results$celltype)

plotDAbeeswarm(da_results, group.by = "celltype")

# identifying signatures of DA populations 



