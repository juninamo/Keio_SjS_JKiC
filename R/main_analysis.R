# ==============================================================================
# scRNA-seq Analysis Pipeline: KEIO SjS Repertoire Study
# Description: This script performs T-Cell sub-analysis for GEX data.
# Author: Jun Inamo
# Date: 2026/01/08
# ==============================================================================


# 1. Setup and Dependencies ----------------------------------------------------
# setwd("~/Desktop/collaborative_reseaerch/repertoire/SjS_Takeshita/JKiC/R_codes")
source("utils.R") 

dir = paste0(getwd(),"/../")

gex = readRDS(file = paste0(dir,"/output/SeuratObj_merged_GEX.rds"))

# filtering for T-cell lineage
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

min_cell_cluster = 50

print(names(table(Idents(gex_t)))[table(Idents(gex_t))>min_cell_cluster])

gex_t = subset(gex_t, 
               idents = names(table(Idents(gex_t)))[table(Idents(gex_t))>min_cell_cluster])


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

Stacked_VlnPlot(gex_t, 
                features = c("CD4", "CD8A", "CCR7", "CXCL13", "PDCD1", "CXCR5", "FOXP3", "PRF1", 
                             "GZMK", "GZMB", "CX3CR1", "GZMH", "GZMA"
                             
                ), 
                x_lab_rotate = TRUE,
                vln_linewidth = 0.2,
                plot_spacing = 0.3) & 
  scale_fill_manual(values = manual_colors)


# Heatmap: Average Expression of T-cell Markers -----------------------------

# Define functional markers for T-cell subsets
t_markers <- c(
  "SELL","CCR7","TCF7","CD69","IL7R","GIMAP5","LEF1", 
  "CD4", "CD8A", "CXCR3","CXCR4","CXCR5","CXCR6","CX3CR1","CCR2","CCR4","CCR6",
  "CCL5","CCL4","CXCL13", "FOXP3","IL2RA", "PDCD1","CTLA4","ICOS","CD40LG", 
  "TIGIT","LAG3", "MAF", "RORC", "GATA3", "TBX21", "HLA-DRA","HLA-DRB1",
  "B3GAT1", "CD27", "CD38", "ITGAE", "KLRB1", "IFI44L","MX1","IFIT3",
  "ZNF683","XCL1", "GZMB","GZMK","GZMA","GZMH","GNLY","PRF1","NKG7","IFNG",
  "TRDV1","TRDV2","TRGV9","SLC4A10","AQP3", "MKI67"
)

# Extract and aggregate data
exp_heat <- as.data.frame(t(GetAssayData(gex_t, slot = "data")[t_markers, ]))
exp_heat$res_cell <- Idents(gex_t)

exp_ave <- aggregate(exp_heat[, 1:(ncol(exp_heat)-1)], list(exp_heat$res_cell), mean)
rownames(exp_ave) <- exp_ave$Group.1
exp_ave <- as.data.frame(t(exp_ave[,-1]))

# Format column names with cluster IDs and names
colnames(exp_ave) <- paste0("T-", cluster_df$cluster[match(colnames(exp_ave), cluster_df$cluster)],
                            ":", cluster_df$clu_name[match(colnames(exp_ave), cluster_df$cluster)])

# Custom function for pheatmap column rotation
draw_colnames_45 <- function (coln, gaps, ...) {
  coord <- pheatmap:::find_coordinates(length(coln), gaps)
  x <- coord$coord - 0.5 * coord$size
  res <- grid::textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"),
                        vjust = 0.75, hjust = 1, rot = 45, gp = grid::gpar(...))
  return(res)
}
assignInNamespace("draw_colnames", "draw_colnames_45", ns = asNamespace("pheatmap"))

pheatmap(
  mat = scale(t(exp_ave)),
  border_color = "white",
  color = colorRampPalette(c("blue", "white", "red"))(100),
  cluster_rows = TRUE, cluster_cols = TRUE,
  fontsize = 9, cellwidth = 9, cellheight = 9,
  main = "Average Marker Expression per Cluster"
)

# TCR Data Processing ----------------------------------------

# Load TCR contig annotations
data_type <- "TCR"
contig_list <- list()

for (sample_id in sample_ids) {
  path <- file.path(dir, "data", sample_id, data_type, "outs", "filtered_contig_annotations.csv")
  if(file.exists(path)) {
    contig_list[[sample_id]] <- read.csv(path)
  }
}

# Combine TCR contigs
combined <- combineTCR(contig_list, samples = sample_ids, ID = sample_ids)

# Strip barcodes for Seurat integration
for (i in seq_along(combined)) {
  combined[[i]] <- stripBarcode(combined[[i]], column = 1, connector = "_", num_connects = 3)
}

combined = lapply(combined, function(df) {
  df$barcode <- paste(df$sample, df$barcode, sep = "_")
  return(df)
})

# Integration of TCR data with Seurat object
seurat_t <- combineExpression(combined, gex_t, cloneCall = "strict", proportion = FALSE)

# Antigen-Specific Clonotype Analysis (Ro60) --------------------------------

tested_t = t_clone_ %>%
  dplyr::filter(orig.ident %in% c("GEX_LB183",
                                  "GEX_LB189",
                                  #"GEX_LB215", # no TCR specificities were identified 
                                  #"GEX_LB216", # no TCR specificities were identified 
                                  "GEX_LB219",
                                  #"GEX_LB221", # no TCR specificities were identified 
                                  "GEX_LB214",
                                  "GEX_LB220")) %>%
  dplyr::filter(clusters_all == "2") %>% # CD4T-enriched ALL-cluster
  merge(.,cluster_df[cluster_df$cell=="T",],by.x="clusters_T",by.y="cluster") %>%
  dplyr::filter(!grepl("CD8$",clu_name) & !grepl("DN$",clu_name) & !grepl("MT+",clu_name) & !grepl("replicating",clu_name)) %>% # CD4T-subclusters in T-cluster
  dplyr::filter(ID!="") %>%
  #dplyr::filter(response=="Ro60") %>%
  .$combination.of.the.nucleotide.and.gene.sequence.CTstrict. %>% unique()
seurat_t <- highlightClonotypes(seurat_t, 
                                cloneCall= "strict", 
                                sequence = tested_t)
seurat_t_full = seurat_t
seurat_t = subset(seurat_t, orig.ident %in% c("GEX_LB183",
                                              "GEX_LB189",
                                              #"GEX_LB215", # no TCR specificities were identified 
                                              #"GEX_LB216", # no TCR specificities were identified 
                                              "GEX_LB219",
                                              #"GEX_LB221", # no TCR specificities were identified 
                                              "GEX_LB214",
                                              "GEX_LB220"))
table(seurat_t@meta.data$orig.ident)
g0 = ggplot() +
  geom_point(
    data = cbind(seurat_t[["umap"]]@cell.embeddings,
                 seurat_t@meta.data) %>% 
      dplyr::mutate(highlight = ifelse(!is.na(highlight),highlight,"hide")) %>%
      dplyr::filter(highlight == "hide"),
    mapping = aes_string(x = "umap_1", y = "umap_2"),
    size = 1, stroke = 0.01, shape = 21, alpha = 0.5, color = "grey75", fill = "grey95"
  ) +
  geom_point(
    data = cbind(seurat_t[["umap"]]@cell.embeddings,
                 seurat_t@meta.data) %>% 
      dplyr::mutate(highlight = ifelse(!is.na(highlight),highlight,"hide")) %>%
      dplyr::filter(highlight != "hide"),
    mapping = aes_string(x = "umap_1", y = "umap_2"),
    size = 2, stroke = 0, shape = 16, alpha = 0.8
  ) +
  labs(
    title = paste0("UMAP colored by tested clones for TCR reporter assay"),
    x = "UMAP1",
    y = "UMAP2"
  ) + 
  theme_classic(base_size = 10) +
  theme(strip.text.x=element_text(size=10, color="black", face="bold"),
        strip.text.y=element_text(size=10, color="black", face="bold"),
        legend.position = "none",
        plot.title = element_text(size=10),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size =10),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust=1),
        legend.text =  element_text(size = 5), 
        legend.key.size = grid::unit(0.4, "lines"),
        legend.title = element_text(size = 7, hjust = 0)) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2), ncol=1),
         fill = guide_legend(override.aes = list(alpha = 1, size = 2)),
         alpha = "none")

ro60_t = t_clone_ %>%
  dplyr::filter(clusters_all == "2") %>% # CD4T-enriched ALL-cluster
  merge(.,cluster_df[cluster_df$cell=="T",],by.x="clusters_T",by.y="cluster") %>%
  dplyr::filter(!grepl("CD8$",clu_name) & !grepl("DN$",clu_name) & !grepl("MT+",clu_name) & !grepl("replicating",clu_name)) %>% # CD4T-subclusters in T-cluster
  dplyr::filter(ID!="") %>%
  dplyr::filter(response=="Ro60") %>%
  .$combination.of.the.nucleotide.and.gene.sequence.CTstrict. %>% unique()
cluster_colors = manual_colors[1:length(ro60_t)]
names(cluster_colors) = ro60_t

nonro60_t = t_clone_ %>%
  dplyr::filter(clusters_all == "2") %>% # CD4T-enriched ALL-cluster
  merge(.,cluster_df[cluster_df$cell=="T",],by.x="clusters_T",by.y="cluster") %>%
  dplyr::filter(!grepl("CD8$",clu_name) & !grepl("DN$",clu_name) & !grepl("MT+",clu_name) & !grepl("replicating",clu_name)) %>% # CD4T-subclusters in T-cluster
  dplyr::filter(ID!="") %>%
  #dplyr::filter(response=="Ro60") %>%
  .$combination.of.the.nucleotide.and.gene.sequence.CTstrict. %>% unique()

seurat_t <- highlightClonotypes(seurat_t, 
                                cloneCall= "strict", 
                                sequence = ro60_t)

g1 = ggplot() +
  geom_point(
    data = cbind(seurat_t[["umap"]]@cell.embeddings,
                 seurat_t@meta.data) %>% 
      dplyr::mutate(highlight = ifelse(!is.na(highlight),highlight,"hide")) %>%
      dplyr::filter(highlight == "hide"),
    mapping = aes_string(x = "umap_1", y = "umap_2"),
    size = 1, stroke = 0.01, shape = 21, alpha = 0.5, color = "grey75", fill = "grey95"
  ) +
  geom_point(
    data = cbind(seurat_t[["umap"]]@cell.embeddings,
                 seurat_t@meta.data) %>% 
      dplyr::mutate(highlight = ifelse(!is.na(highlight),highlight,"hide")) %>%
      dplyr::filter(highlight != "hide"),
    mapping = aes_string(x = "umap_1", y = "umap_2", fill = "highlight"),
    size = 2, stroke = 0.2, shape = 21, alpha = 0.8, color = "black"
  ) +
  labs(
    title = paste0("UMAP colored by Ro60-responded clones"),
    subtitle = "colored by clone IDs",
    x = "UMAP1",
    y = "UMAP2"
  ) + 
  scale_color_manual(values = c(cluster_colors,hide="grey95")) +
  theme_classic(base_size = 10) +
  theme(strip.text.x=element_text(size=10, color="black", face="bold"),
        strip.text.y=element_text(size=10, color="black", face="bold"),
        legend.position = "bottom",
        plot.title = element_text(size=10),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size =10),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust=1),
        legend.text =  element_text(size = 5), 
        legend.key.size = grid::unit(0.4, "lines"),
        legend.title = element_text(size = 7, hjust = 0)) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2), ncol=1),
         fill = guide_legend(override.aes = list(alpha = 1, size = 2)),
         alpha = "none")

table(seurat_t@meta.data$new_cluster)
table(seurat_t@meta.data$seurat_clusters)

ord = cbind(seurat_t[["umap"]]@cell.embeddings,
            seurat_t@meta.data) %>% 
  #merge(.,cluster_df[cluster_df$cell=="T",],by.x="seurat_clusters",by.y="cluster") %>%
  dplyr::mutate(cluster = Idents(seurat_t)) %>%
  #dplyr::filter(!grepl("CD8$",clu_name) & !grepl("DN$",clu_name) & !grepl("MT+",clu_name) & !grepl("replicating",clu_name)) %>% # CD4T-subclusters in T-cluster
  dplyr::filter(!is.na(highlight)) %>%
  .$cluster %>% unique() %>% as.character()
g2 = ggplot() +
  geom_point(
    data = cbind(seurat_t[["umap"]]@cell.embeddings,
                 seurat_t@meta.data) %>%
      #merge(.,cluster_df[cluster_df$cell=="T",],by.x="seurat_clusters",by.y="cluster") %>%
      dplyr::mutate(highlight = ifelse(!is.na(highlight),highlight,"hide")) %>%
      dplyr::mutate(cluster = seurat_t@meta.data$seurat_clusters) %>%
      dplyr::mutate(cluster = factor(ifelse(highlight=="hide",NA,cluster))) %>%
      dplyr::filter(highlight == "hide"),
    mapping = aes_string(x = "umap_1", y = "umap_2"),
    size = 1, stroke = 0.001, shape = 21, alpha = 0.2, color = "grey75", fill = "grey95"
  ) +
  geom_point(
    data = cbind(seurat_t[["umap"]]@cell.embeddings,
                 seurat_t@meta.data) %>% 
      dplyr::mutate(cluster = seurat_t@meta.data$seurat_clusters) %>%
      #merge(.,cluster_df[cluster_df$cell=="T",],by.x="seurat_clusters",by.y="cluster") %>%
      dplyr::filter(!is.na(highlight)),
    mapping = aes_string(x = "umap_1", y = "umap_2", fill = "cluster"),
    size = 2, stroke = 0.2, shape = 21, alpha = 1, color = "black",
  ) +
  labs(
    title = paste0("UMAP colored by Ro60-responded clones"),
    subtitle = "colored by T-cell clusters",
    x = "UMAP1",
    y = "UMAP2"
  ) + 
  scale_fill_manual(values = manual_colors,
                    name ="cluster"
                    #, labels = paste0("T-",cluster_df[cluster_df$cell=="T" & cluster_df$cluster %in% ord,]$cluster,":",cluster_df[cluster_df$cell=="T" & cluster_df$cluster %in% ord,]$clu_name)
  ) +
  theme_classic(base_size = 10) +
  theme(strip.text.x=element_text(size=10, color="black", face="bold"),
        strip.text.y=element_text(size=10, color="black", face="bold"),
        legend.position = "bottom",
        plot.title = element_text(size=10),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size =10),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust=1),
        legend.text =  element_text(size = 5), 
        legend.key.size = grid::unit(0.4, "lines"),
        legend.title = element_text(size = 7, hjust = 0)) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2), ncol=1),
         fill = guide_legend(override.aes = list(alpha = 1, size = 2)),
         alpha = "none")
g0 / g1 / g2
