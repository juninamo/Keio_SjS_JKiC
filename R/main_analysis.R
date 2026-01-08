# ==============================================================================
# scRNA-seq Analysis Pipeline: KEIO SjS Repertoire Study
# Description: This script performs QC, Integration (Harmony), and 
#              Clustering for GEX data.
# Author: Jun Inamo
# Date: 2026/01/08
# ==============================================================================


# 1. Setup and Dependencies ----------------------------------------------------
# setwd("~/Desktop/collaborative_reseaerch/repertoire/SjS_Takeshita/JKiC/R_codes")
source("utils.R") 

data_type = "GEX"
dir = paste0(getwd(),"/../")

sample_ids = list.files(path = paste0(dir,"/data"), pattern = "^LB", full.names = FALSE)
successful_gex_list <- list()

# 2. Data Loading and Preprocessing --------------------------------------------
for (sample_id in sample_ids) {
  # Load the dataset
  gex.data <- Read_CellBender_h5_Mat(paste0(dir,"/data/",sample_id,"/",data_type,"/cellbender_feature_bc_matrix_filtered.h5"))
  # Initialize the Seurat object with the raw (non-normalized data).
  gex <- CreateSeuratObject(counts = gex.data, 
                            project = paste0(data_type,"_",sample_id) 
                            # min.cells = 3, min.features = 200
  )
  successful_gex_list[[sample_id]] <- gex
  
}
print("Success sample_ids; ")
print(names(successful_gex_list))

gex <- Merge_Seurat_List(successful_gex_list, 
                         add.cell.ids = names(successful_gex_list),
                         merge.data = TRUE,
                         project = "SjS_JKiC")

# 3. Quality Control (QC) ------------------------------------------------------
gex[["percent.mt"]] <- PercentageFeatureSet(gex, pattern = "^MT-")

# Filter based on features and mitochondrial content
gex <- subset(gex, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20)

# Remove Doublets
doubletcells = read.table(paste0(dir,"/output/scDblFinderClass.txt"), header=TRUE) %>%
  dplyr::filter(scDblFinder.class == "doublet") %>%
  .$barcode
gex <- subset(gex, cells = doubletcells, invert = TRUE)

# Calculate Ribosomal content
rb.genes <- rownames(gex)[grep("^RP[SL]", rownames(gex))]
gex[["percent.ribo"]] <- colSums(gex[rb.genes, ]) / Matrix::colSums(gex) * 100

# 4. Normalization and Integration ---------------------------------------------
gex <- NormalizeData(gex, normalization.method = "LogNormalize", scale.factor = 10000)

gex <- FindVariableFeatures(gex, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(gex)
gex <- ScaleData(gex, features = all.genes)

set.seed(1234)
gex <- RunPCA(gex, features = VariableFeatures(object = gex))

# Batch Correction using Harmony
set.seed(1234)
gex <- gex %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE)

# Dimensional Reduction and Clustering
DIM_PCA <- 8
gex <- FindNeighbors(gex, reduction = "harmony", 
                     dims = 1:DIM_PCA,
                     k.param = 30,
                     verbose = FALSE
)
gex <- RunUMAP(gex, reduction = "harmony", 
               dims = 1:DIM_PCA,
               n.neighbors = 30L,
               min.dist = min.dist,
               verbose = FALSE)

resolution=0.1
gex <- FindClusters(gex, reduction = "harmony", 
                    dims = 1:DIM_PCA,
                    resolution = resolution,
                    verbose = FALSE)

gex <- JoinLayers(gex, overwrite = TRUE)

gex@meta.data$new_cluster = dplyr::case_when(
  gex@meta.data$seurat_clusters %in% c("0") ~ "T cell",
  gex@meta.data$seurat_clusters %in% c("1","8") ~ "Plasma cell",
  gex@meta.data$seurat_clusters %in% c("2","5") ~ "Aciner",
  gex@meta.data$seurat_clusters %in% c("3") ~ "B cell",
  gex@meta.data$seurat_clusters %in% c("4") ~ "Duct & Stem",
  gex@meta.data$seurat_clusters %in% c("6") ~ "Myeloid",
  gex@meta.data$seurat_clusters %in% c("7") ~ "Fibroblast & Stem",
  gex@meta.data$seurat_clusters %in% c("9") ~ "Plasmablast"
)

ggplot() +
  geom_point(
    data = data.frame(gex[["umap"]]@cell.embeddings,
                      res_cell = Idents(gex)) %>%
      dplyr::sample_frac(1L),
    mapping = aes_string(x = "umap_1", y = "umap_2", color = "res_cell"),
    size = 1, stroke = 0, shape = 16, alpha = 0.5
  ) +
  scale_color_manual(values=manual_colors) +
  scale_fill_manual(values=manual_colors) +
  labs(
    x = "UMAP1",
    y = "UMAP2",
    title = ""
  ) +
  coord_fixed() +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "right",
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(size=15),
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size =15),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    legend.text =  element_text(size = 12),
    legend.key.size = grid::unit(0.5, "lines"),
    legend.title = element_text(size = 0.8, hjust = 0)
  ) + guides(colour = guide_legend(override.aes = list(size=3,alpha=1),
                                   title = "",
                                   ncol = 1),
             fill="none")

ggplot() +
  geom_point(
    data = data.frame(gex[["umap"]]@cell.embeddings,
                      res_cell = gex@meta.data$new_cluster) %>%
      dplyr::sample_frac(1L),
    mapping = aes_string(x = "umap_1", y = "umap_2", color = "res_cell"),
    size = 1, stroke = 0, shape = 16, alpha = 0.5
  ) +
  scale_color_brewer(palette = "Set1") +
  labs(
    x = "UMAP1",
    y = "UMAP2",
    title = ""
  ) +
  coord_fixed() +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "right",
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(size=15),
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size =15),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    legend.text =  element_text(size = 12),
    legend.key.size = grid::unit(0.5, "lines"),
    legend.title = element_text(size = 0.8, hjust = 0)
  ) + guides(colour = guide_legend(override.aes = list(size=3,alpha=1),
                                   title = "",
                                   ncol = 1),
             fill="none")

# 5. T-Cell Sub-analysis -------------------------------------------------------
message("Subsetting T-cells...")

# Broad filtering for T-cell lineage
gex_t = subset(gex, 
               subset = CD19 < 0.5 & MS4A1 < 0.5 & CD79A < 0.5 & MZB1 < 0.5 & IGHG1 < 0.5 & IGHA1 < 0.5 & IGHD < 0.5 & CD3E > 1 & CD14 < 0.5 & CD68 < 0.5 & TYROBP < 0.5 & NCAM1 < 0.5 & MUC7 < 0.5 & MUC5B < 0.5 & KRT19 < 0.5 & KRT5 < 0.5 & KRT14 < 0.5)
# include cells with TCR
gex_t = subset(gex_t, cells = cellranger_rep[cellranger_rep$cell_type=="T",]$sample_barcode)

# Re-process T-cell subset
gex_t <- NormalizeData(gex_t, normalization.method = "LogNormalize", scale.factor = 10000)
gex_t <- FindVariableFeatures(gex_t, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(gex_t)
gex_t <- ScaleData(gex_t, features = all.genes)

set.seed(1234)
gex_t <- RunPCA(gex_t, features = VariableFeatures(object = gex_t))

set.seed(1234)
gex_t <- gex_t %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE)

dim=20
min.dist=0.3
resolution=1.00
gex_t <- FindNeighbors(gex_t, reduction = "harmony", 
                       dims = 1:dim,
                       k.param = 30,
                       verbose = FALSE
)
gex_t <- RunUMAP(gex_t, reduction = "harmony", 
                 dims = 1:dim,
                 n.neighbors = 30L,
                 min.dist = min.dist,
                 verbose = FALSE)
gex_t <- FindClusters(gex_t, reduction = "harmony", 
                      dims = 1:dim , resolution = resolution
)

# 6. Visualization -------------------------------------------------------------
ggplot() +
  geom_point(
    data = data.frame(gex_t[["umap"]]@cell.embeddings,
                      res_cell = Idents(gex_t)),
    mapping = aes_string(x = "umap_1", y = "umap_2", color = "res_cell"),
    size = 2, stroke = 0, shape = 16, alpha = 0.5
  ) +
  scale_color_manual(values=manual_colors) +
  scale_fill_manual(values=manual_colors) +
  labs(
    x = "",
    y = "",
    title = ""
  ) +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "right",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(color="black", size=10)
  ) + guides(colour = guide_legend(override.aes = list(size=3,alpha=1),
                                   title = "",
                                   ncol = 1))

