# ==============================================================================
# utils.R  -  Helper functions and shared objects for the SjD scRNA-seq / TCR
#             repertoire analysis (Takeshita, Inamo et al., Sci. Adv. 2026).
#
# This is a trimmed version that keeps ONLY the libraries, helper functions and
# shared objects required by the analysis notebooks:
#   - 01_preprocessing.ipynb
#   - 02_figure_analysis.ipynb
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Libraries
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)        # single-cell toolkit
  library(harmony)       # batch integration (RunHarmony)
  library(scCustomize)   # DotPlot_scCustom / Stacked_VlnPlot / Plot_Density_Custom / Read_CellBender_h5_Mat
  library(scRepertoire)  # combineTCR / combineExpression / highlightClonotypes
  library(tidyverse)     # dplyr / ggplot2 / tidyr / tibble / stringr / purrr
  library(magrittr)      # %>% and %<>%
  library(ggrepel)       # repelled text labels
  library(patchwork)     # plot composition with / and +
  library(pheatmap)      # heatmaps
  library(ggpubr)        # ggpaired() / stat_cor() used in the TCR comparisons
  library(ggsci)         # scale_fill_jco() palette used in the TCR comparisons
})
# Namespaced (called via pkg::fun, no need to attach): MASS, scales,
# grid, immunarch, ggVennDiagram, ggrastr, Matrix.

# ------------------------------------------------------------------------------
# 2. Project paths
# ------------------------------------------------------------------------------
# `dir` points to the project root (one level above the R working directory).

dir = paste0(getwd(),"/../")
dir.create(paste0(dir,"/output"), showWarnings = FALSE)
dir.create(paste0(dir,"/tmp"), showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 3. Color palette and marker genes
# ------------------------------------------------------------------------------

manual_colors = c(
  "0" = "#E41A1C",
  "1" = "#377EB8",
  "2" = "#4DAF4A",
  "3" = "#984EA3",
  "4" = "#FF7F00",
  "5" = "#FFFF33",
  "6" = "#A65628",
  "7" = "lightgrey",
  "8" = "#999999",
  "9" = "#66C2A5",
  "10" = "#67000D",
  "11" = "#8DA0CB",
  "12" = "#FFD92F",
  "13" = "#A6D854",
  "14" = "#E78AC3",
  "15" = "#FC8D62",
  "16" = "darkgrey",
  "17" = "#FEB24C",
  "18" = "#377EB8",
  "19" = "lightblue",
  "20" = "#FDE0EF",
  "21" = "#B8E186",
  "22" = "#66C2A5",
  "23" = "#A6D855",
  "24" = "#E78AC4",
  "25" = "#FC8D63",
  "26" = "brown",
  "27" = "#FEB25C",
  "28" = "#377EB9",
  "29" = "lightgreen",
  "30" = "#FDE1EF"
)

all_markers = c("PTPRC", # CD45: a glycoprotein expressed on all nucleated hematopoietic cells
                "CD19","MS4A1", "CD79A", # CD79A: BCR-complex-associated protein alpha chain
                "IGHG1","IGHA1","SDC1", # SDC1(CD138): plasma cells
                "RPN2","XBP1","PRDX4", "MZB1", "CD38", # plasmablasts 
                "CD3E", "CD4", "CD8A",
                "CD14", "FCGR3A", "S100A8", "CST3", # (FCGR3A=CD16)  #  LYZ and CST3: monocyte and pDCs # S100A8 and S100A9: neutrophils, monocytes and pDCs
                "CD68", # monocyte and macrophage
                "ITGAX",# "CD11c": monocyte and macrophage
                "NCAM1", # CD56: NK cell
                "PDPN", "THY1", "ICAM1","VCAM1",# THY1(CD90): Fibroblast 
                
                "EPCAM", # EPCAM: epithelial cell adhesion molecule
                "LYZ","PIP","MUC7", # serous gland
                "MUC5B","CEACAM6","BPIFB2", # mucous gland
                "TFCP2L1","KRT19", # duct
                "KRT5","KRT14","TP63","CAV1", # stem cell
                "ACTA2", # myoepithelium
                "CDH5", # VE-cadherin: Endothelial
                "MCAM" # CD146: Pericyte
)

# ------------------------------------------------------------------------------
# 4. Sample list, cluster annotation and 10x cellranger TCR/BCR annotation
# ------------------------------------------------------------------------------
# Sample IDs are the LB* folders under <project>/data.

sample_ids = list.files(path = paste0(dir,"/data"), pattern = "^LB", full.names = FALSE)

sample_colors <- manual_colors[1:length(sample_ids)]
names(sample_colors) = sample_ids

cluster_df = read.table("./cluster_annotation.txt",sep = "\t",header = TRUE)
cluster_df$cluster <- as.character(cluster_df$cluster)

cellranger_rep <- do.call("rbind", lapply(sample_ids, function(sample_id) {
  data <- rbind(
    read.csv(paste0(dir, "/data/", sample_id, "/BCR/outs/filtered_contig_annotations.csv")) %>%
      dplyr::mutate(cell_type = "B/plasma"),
    read.csv(paste0(dir, "/data/", sample_id, "/TCR/outs/filtered_contig_annotations.csv")) %>%
      dplyr::mutate(cell_type = "T")
  )
  data$sample_id <- sample_id
  return(data)
}))

# print(table(cellranger_rep$sample_id))
# print(table(cellranger_rep$chain))
# print(table(cellranger_rep$cell_type,cellranger_rep$sample_id,cellranger_rep$chain))

cellranger_rep = cellranger_rep %>%
  dplyr::distinct(cell_type,barcode,.keep_all = TRUE) %>%
  dplyr::mutate(sample_barcode = paste(sample_id,barcode,sep = "_")) %>%
  dplyr::group_by(sample_barcode) %>%
  dplyr::mutate(count = dplyr::n()) %>%
  dplyr::filter(count == 1) %>% # remove cells with annotation both T and B/plasma
  as.data.frame() %>%
  dplyr::distinct(sample_barcode,cell_type)
rownames(cellranger_rep) = cellranger_rep$sample_barcode


# ------------------------------------------------------------------------------
# 5. Doublet barcodes flagged by scDblFinder (computed in QC_CellBender.Rmd)
# ------------------------------------------------------------------------------

doubletcells = read.table(paste0(dir,"/output/scDblFinderClass.txt"), header=TRUE) %>%
  dplyr::filter(scDblFinder.class == "doublet") %>%
  .$barcode

# ------------------------------------------------------------------------------
# 6. Helper functions
# ------------------------------------------------------------------------------

Merge_Seurat_List <- function(
    list_seurat,
    add.cell.ids = NULL,
    merge.data = TRUE,
    project = "SeuratProject"
) {
  # Check list_seurat is list
  if (!inherits(x = list_seurat, what = "list")) {
    cli::cli_abort(message = "{.code list_seurat} must be environmental variable of class {.val list}")
  }
  
  # Check list_seurat is only composed of Seurat objects
  for (i in 1:length(x = list_seurat)) {
    if (!inherits(x = list_seurat[[i]], what = "Seurat")) {
      cli::cli_abort("One or more of entries in {.code list_seurat} are not objects of class {.val Seurat}")
    }
  }
  
  # Check all barcodes are unique to begin with
  duplicated_barcodes <- list_seurat %>%
    lapply(colnames) %>%
    unlist() %>%
    duplicated() %>%
    any()
  
  if (duplicated_barcodes && is.null(x = add.cell.ids)) {
    cli::cli_abort(message = c("There are overlapping cell barcodes present in the input objects",
                          "i" = "Please rename cells or provide prefixes to {.code add.cell.ids} parameter to make unique.")
    )
  }
  
  # Check right number of suffix/prefix ids are provided
  if (!is.null(x = add.cell.ids) && length(x = add.cell.ids) != length(x = list_seurat)) {
    cli::cli_abort(message = "The number of prefixes in {.code add.cell.ids} must be equal to the number of objects supplied to {.code list_seurat}.")
  }
  
  # Rename cells if provided
  list_seurat <- lapply(1:length(x = list_seurat), function(x) {
    list_seurat[[x]] <- RenameCells(object = list_seurat[[x]], add.cell.id = add.cell.ids[x])
  })
  
  current_index <- 1
  # Merge objects
  merged_object <- purrr::reduce(list_seurat, function(x, y) {
    message("Merging object with index: ", current_index)
    current_index <<- current_index + 1  # advance the index (<<- updates the enclosing variable)
    merge(x = x, y = y, merge.data = merge.data, project = project)
  })
}

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

count_seq_cells <- function(df){
  df %>%
    dplyr::distinct(barcode) %>%
    nrow()
}

sample_data <- function(data, num_samples, seed) {
  set.seed(seed)
  # sample one index from each batch first
  unique_batches <- unique(data$batch)
  
  if(num_samples < length(unique_batches)) {
    stop("num_samples must be greater than or equal to the number of unique batches")
  }
  
  initial_samples <- sapply(unique_batches, function(b) {
    sample(which(data$batch == b), 1)
  })
  
  # then draw the remaining samples at random (excluding the ones already chosen)
  remaining_samples_needed <- num_samples - length(initial_samples)
  remaining_samples <- setdiff(1:nrow(data), initial_samples)
  
  remaining_samples_selected <- sample(remaining_samples, 
                                       min(remaining_samples_needed, length(remaining_samples)),
                                       replace = TRUE)
  
  # combine all sampled indices
  all_samples <- c(initial_samples, remaining_samples_selected)
  
  return(all_samples)
}

exportClones <- function(input.data,
                         format = "paired",
                         group.by = NULL,
                         write.file = TRUE,
                         dir = NULL,
                         file.name = "clones.csv") {
  
  # Validate format parameter
  format <- match.arg(format, c("paired", "airr", "TCRMatch", "tcrpheno", "immunarch"))
  
  # Select the appropriate internal export function
  exportFunc <- switch(format,
                       "paired"     = .pairedExport,
                       "airr"       = .airrExport,
                       "TCRMatch"   = .tcrMatchExport,
                       "tcrpheno"   = .tcrPhenoExport,
                       "immunarch"  = .immunarchExport
  )
  
  # Generate the data matrix/list
  output_data <- exportFunc(input.data, group.by)
  
  # Replace string "NA" with actual NA values
  if (is.data.frame(output_data)) {
    output_data[output_data == "NA"] <- NA
  }
  
  if (!write.file) {
    return(output_data)
  }
  
  # Handle file writing
  if (is.null(dir)) {
    dir <- "."
  }
  filepath <- file.path(dir, file.name)
  
  # Immunarch format returns a list of data frames. To save as a single CSV,
  # we bind them together and add a 'Sample' identifier column.
  if (format == "immunarch") {
    bound_data <- dplyr::bind_rows(output_data$data, .id = "Sample")
    write.csv(bound_data, file = filepath, row.names = FALSE)
  } else {
    write.csv(output_data, file = filepath)
  }
}
environment(exportClones) <- asNamespace("scRepertoire")

.immunarchExport <- function(input.data, group.by) {
  df_list <- .dataWrangle(input.data, group.by, "CTgene", "both")
  meta <- data.frame(Sample = names(df_list))
  
  data_out <- lapply(df_list, function(x) {
    # Summarize by clonotype (CTstrict)
    result <- x %>%
      dplyr::group_by(.data[["CTstrict"]]) %>%
      dplyr::summarise(
        Clones  = dplyr::n(),
        barcode = paste(barcode, collapse = ";"),
        CTaa    = dplyr::first(.data[["CTaa"]]),
        CTnt    = dplyr::first(.data[["CTnt"]]),
        CTgene  = dplyr::first(.data[["CTgene"]]),
        .groups = "drop"
      ) %>%
      dplyr::mutate(Proportion = Clones / sum(Clones))
    
    # Determine chain type from first half
    first_genes <- ._split_and_pad(result$CTgene, "_", 2)[, 1]
    chain_type  <- substr(first_genes, 1, 3)
    
    # Impute Missing Chain info (fix condition)
    if (any(is.na(chain_type)) && length(unique(chain_type)) == 2) {
      tbl.store <- rev(sort(table(chain_type, useNA = "ifany")))
      if (length(tbl.store) == 2) {
        chain_type[is.na(chain_type)] <- names(tbl.store)[1]
      }
    }
    
    # Split sequence/gene info
    aa_split   <- ._split_and_pad(result$CTaa,   "_", 2)
    nt_split   <- ._split_and_pad(result$CTnt,   "_", 2)
    gene_split <- ._split_and_pad(result$CTgene, "_", 2)
    
    .is_layout3_vec <- function(v_gene_first_token) {
      grepl("^(TRAV|TRGV|IGKV|IGLV)", v_gene_first_token %||% "", perl = TRUE)
    }
    
    # Adaptive per-row ordering: alpha/light first
    first.gene.pos <- ifelse(chain_type %in% c("TRA","TRG","IGK","IGL"), 1L, 2L)
    na_pos <- is.na(first.gene.pos)
    if (any(na_pos)) {
      left_is3  <- .is_layout3_vec(gsub("\\..*$", "", gene_split[, 1]))
      right_is3 <- .is_layout3_vec(gsub("\\..*$", "", gene_split[, 2]))
      first.gene.pos[na_pos] <- ifelse(left_is3[na_pos], 1L, 2L)
    }
    second.gene.pos <- 3L - first.gene.pos
    
    n    <- nrow(result)
    ridx <- seq_len(n)
    
    # Per-row pick halves
    aa_first   <- as.character(aa_split[cbind(ridx, first.gene.pos)])
    aa_second  <- as.character(aa_split[cbind(ridx, second.gene.pos)])
    nt_first   <- as.character(nt_split[cbind(ridx, first.gene.pos)])
    nt_second  <- as.character(nt_split[cbind(ridx, second.gene.pos)])
    g_first_s  <- as.character(gene_split[cbind(ridx, first.gene.pos)])
    g_second_s <- as.character(gene_split[cbind(ridx, second.gene.pos)])
    
    # Chain presence flags (entire half missing?)
    present1 <- !(is.na(g_first_s)  | g_first_s  == "" | g_first_s  == "NA")
    present2 <- !(is.na(g_second_s) | g_second_s == "" | g_second_s == "NA")
    
    # Tokenize per selected half
    genes1 <- ._split_and_pad(g_first_s,  "[.]", 4)
    genes2 <- ._split_and_pad(g_second_s, "[.]", 4)
    
    # Remap to V/D/J/C (D=NA for 3-gene chains)
    g1 <- .remap_VDJC(genes1)
    g2 <- .remap_VDJC(genes2)
    
    # Display helpers
    pair_strict <- function(lhs, rhs) {
      lhs <- ifelse(is.na(lhs) | lhs == "", "NA", as.character(lhs))
      rhs <- ifelse(is.na(rhs) | rhs == "", "NA", as.character(rhs))
      paste(lhs, rhs, sep = ";")
    }
    
    # For genes: if chain absent -> "NA"
    show_gene <- function(g, present) {
      out <- ifelse(is.na(g) | g == "", "NA", g)
      out[!present] <- "NA"
      out
    }
    # For D gene: if chain present but D missing (3-gene) -> "None"; if chain absent -> "NA"
    show_D <- function(d, present) {
      out <- d
      out <- ifelse(is.na(out) | out == "", "None", out)  
      out[!present] <- "NA"                               
      out
    }
    
    # Build display vectors
    V1 <- show_gene(g1$V, present1); V2 <- show_gene(g2$V, present2)
    D1 <- show_D(  g1$D, present1); D2 <- show_D(  g2$D, present2)
    J1 <- show_gene(g1$J, present1); J2 <- show_gene(g2$J, present2)
    C1 <- show_gene(g1$C, present1); C2 <- show_gene(g2$C, present2)
    
    # CDR3: if chain absent or blank -> "NA"
    aa1 <- ifelse(present1 & !(is.na(aa_first)  | aa_first  == ""), aa_first,  "NA")
    aa2 <- ifelse(present2 & !(is.na(aa_second) | aa_second == ""), aa_second, "NA")
    nt1 <- ifelse(present1 & !(is.na(nt_first)  | nt_first  == ""), nt_first,  "NA")
    nt2 <- ifelse(present2 & !(is.na(nt_second) | nt_second == ""), nt_second, "NA")
    
    # Assemble; always "lhs;rhs"
    out <- data.frame(
      Clones      = result[["Clones"]],
      Proportion  = result[["Proportion"]],
      CDR3.nt     = pair_strict(nt1, nt2),
      CDR3.aa     = pair_strict(aa1, aa2),
      V.name      = pair_strict(V1, V2),
      D.name      = pair_strict(D1, D2),  # preserves "None" vs "NA"
      J.name      = pair_strict(J1, J2),
      C.name      = pair_strict(C1, C2),
      Barcode     = result[["barcode"]],
      stringsAsFactors = FALSE
    )
    out
  })
  
  names(data_out) <- names(df_list)
  list(data = data_out, meta = meta)
}
environment(.immunarchExport) <- asNamespace('scRepertoire')

.dataWrangle <- function(df, split.by, cloneCall, chain) {
  df <- .list.input.return(df, split.by)
  df <- .checkBlanks(df, cloneCall)
  for (i in seq_along(df)) {
    if (chain != "both") {
      df[[i]] <- .offTheChain(df[[i]], chain, cloneCall)
    }
  }
  return(df)
}
environment(.dataWrangle) <- asNamespace('scRepertoire')

.list.input.return <- function(df, split.by) {
  if (.is.seurat.or.se.object(df)) {
    if(is.null(split.by)){
      split.by <- "ident"
    }
    df <- .expression2List(df, split.by)
  } 
  df
}
environment(.list.input.return) <- asNamespace('scRepertoire')

.is.seurat.or.se.object <- function(obj) {
  .is.seurat.object(obj) || .is.se.object(obj)
}
environment(.is.seurat.or.se.object) <- asNamespace('scRepertoire')

.is.seurat.object <- function(obj) inherits(obj, "Seurat")
environment(.is.seurat.object) <- asNamespace('scRepertoire')

.is.se.object <- function(obj) inherits(obj, "SummarizedExperiment")
environment(.is.se.object) <- asNamespace('scRepertoire')

.checkBlanks <- function(df, cloneCall) {
  count <- NULL
  for (i in seq_along(df)) {
    # First, check if there are no rows
    if (nrow(df[[i]]) == 0) {
      count <- c(i, count)
    } 
    # If there are rows, then proceed with blank checks
    else if (length(df[[i]][,cloneCall]) == length(which(is.na(df[[i]][,cloneCall]))) |
             length(which(!is.na(df[[i]][,cloneCall]))) == 0) {
      count <- c(i, count)
    } else {
      next()
    }
  }
  if (!is.null(count)) {
    df <- df[-count]
  }
  return(df)
}
environment(.checkBlanks) <- asNamespace('scRepertoire')

._split_and_pad <- function(x, split, n_cols) {
  s <- strsplit(x, split)
  # Create a matrix by safely subsetting each list element up to n_cols
  mat <- t(sapply(s, `[`, seq_len(n_cols)))
  mat[mat == "NA"] <- NA
  return(mat)
}
environment(._split_and_pad) <- asNamespace('scRepertoire')

.remap_VDJC <- function(genes_mat) {
  # Coerce to a 4-col character matrix, pad if needed
  if (is.null(genes_mat)) stop(".remap_VDJC: genes_mat is NULL")
  genes_mat <- as.matrix(genes_mat)
  storage.mode(genes_mat) <- "character"
  if (ncol(genes_mat) < 4) {
    pad <- matrix(NA_character_, nrow(genes_mat), 4 - ncol(genes_mat))
    genes_mat <- cbind(genes_mat, pad)
  }
  
  # Clean tokens: drop everything after first ';', trim, map "", "NA" -> NA
  clean_token <- function(x) {
    x <- sub(";.*", "", x, perl = TRUE)
    x <- trimws(x)
    x[x == "" | x == "NA"] <- NA_character_
    x
  }
  genes_mat[] <- clean_token(genes_mat)
  
  vcol <- genes_mat[, 1]
  
  # Detect 3-gene layouts (TRA/TRG/IGK/IGL); handle NAs
  layout3 <- !is.na(vcol) & grepl("^(TRAV|TRGV|IGKV|IGLV)", vcol, perl = TRUE)
  
  # Initialize as 4-gene (V,D,J,C)
  V <- vcol
  D <- genes_mat[, 2]
  J <- genes_mat[, 3]
  C <- genes_mat[, 4]
  
  # Overwrite 3-gene rows: V, (no D), J <- col2, C <- col3
  if (any(layout3, na.rm = TRUE)) {
    idx <- which(layout3)
    D[idx] <- NA_character_
    J[idx] <- genes_mat[idx, 2]
    C[idx] <- genes_mat[idx, 3]
  }
  
  # Ensure character vectors with row length
  stopifnot(length(V) == nrow(genes_mat),
            length(D) == nrow(genes_mat),
            length(J) == nrow(genes_mat),
            length(C) == nrow(genes_mat))
  
  list(V = V, D = D, J = J, C = C)
}
environment(.remap_VDJC) <- asNamespace('scRepertoire')
