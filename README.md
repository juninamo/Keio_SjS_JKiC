# SjS Repertoire ScRNA-seq Analysis Pipeline

This repository contains the R analysis pipeline for the study: **"Single-cell Transcriptomic Profiling of SjS Repertoire"**. The workflow covers the processing of Gene Expression (GEX) data, from raw counts (CellBender output) to integrated clustering and T-cell sub-population analysis.

## 📌 Overview

This pipeline is designed to ensure transparency and reproducibility of the findings. It utilizes the `Seurat` framework and `Harmony` for batch effect correction across multiple samples.

**Key Features:**

* Automated loading of **CellBender** `.h5` filtered matrices.
* Standardized Quality Control (QC) and filtering (Mitochondrial/Ribosomal content).
* Doublet removal integration.
* Batch correction and data integration using **Harmony**.
* Cell type-specific sub-clustering (focused on T-cells).
* Visualization tools (UMAP, DotPlots, Stacked Violin Plots, and Heatmaps).

---

## 📂 Repository Structure

To run the script successfully, please organize your local directory as follows:

```text
.
├── R/
│   ├── main_analysis.R    # Principal analysis script
│   └── utils.R           # Custom utility functions (required)
├── data/                  # Input data folder
│   ├── LB215/             # Sample ID
│   │   └── GEX/           # Data type
│   │       └── cellbender_feature_bc_matrix_filtered.h5
│   └── ...
├── output/                # Generated results (RDS objects)
│   └── figure/            # Exported PDF/PNG plots
└── README.md

```

---

## 🛠 Prerequisites

### Environment

* **R version:** 4.1 or higher is recommended.
* **Required Packages:**
```r
install.packages(c("Seurat", "harmony", "dplyr", "ggplot2", "magrittr", "pheatmap", "tidyr", "ggrepel"))
# For DotPlot_scCustom and other visualization extensions
install.packages("SeuratCustomExtensions") 

```



---

## 🚀 Getting Started

1. **Clone the repository:**
```bash
git clone https://github.com/yourusername/SjS_Repertoire_Analysis.git
cd SjS_Repertoire_Analysis

```


2. **Prepare your data:** Place your `.h5` files and the `scDblFinderClass.txt` (doublet info) into the `data/` and `output/` folders respectively as shown in the structure above.
3. **Run the analysis:** Open `R/main_analysis.R` in RStudio and source the file, or run via terminal:
```bash
Rscript R/main_analysis.R

```



---

## 📊 Methodology Highlights

### 1. Pre-processing & QC

Cells are filtered based on:

* `nFeature_RNA`: > 200 and < 5,000
* `percent.mt`: < 20%
* Doublets identified via `scDblFinder` are excluded.

### 2. Integration

Batch effects between individual samples are corrected using the **Harmony** algorithm on the top 2000 variable features.

### 3. T-cell Sub-analysis

A high-resolution sub-clustering is performed specifically on the T-cell population to identify rare subsets and functional states (e.g., Tph, Tfh, Trm).

---

## 📝 Code Availability Statement

Following the reporting guidelines:

* **Newly generated code:** All custom algorithms and scripts essential for replicating the main findings are included in this repository.
* **Reused code:** Dependencies and cited packages are listed in the script headers and `sessionInfo`.
* **Licensing:** This code is released under the [MIT License](https://www.google.com/search?q=LICENSE).

---

## 📧 Contact

For questions or collaborations, please reach out to:
**[Your Name]** - [Your Email Address]
**Institution:** [Your Department/University]

---

### Would you like me to:

1. **Create a `LICENSE` file** (e.g., MIT or Apache 2.0) to include in your repo?
2. **Add a `requirements.txt` or `renv.lock` section** to help others perfectly replicate your R environment?