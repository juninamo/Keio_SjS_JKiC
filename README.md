# SjD Repertoire scRNA-seq / TCR Analysis Pipeline

R analysis pipeline for the study:

> **Takeshita M, Inamo J, Wakui S, Nagashima R, Nishino T, Tsunoda K, Usuda S, Inokuchi H, Ishigaki K, Sasaki T, Kagoya Y, Suzuki K, Kaneko Y.**
> *Identification of a shared antigen linking CD4⁺ T and B cell pathology in Sjögren's disease.*
> **Science Advances** 12, eaeb2491 (2026). doi:[10.1126/sciadv.aeb2491](https://doi.org/10.1126/sciadv.aeb2491)

The workflow covers single-cell gene expression (GEX) and TCR-repertoire analysis of
salivary-gland tissue from patients with Sjögren's disease (SjD), and reproduces the
figures of the manuscript.

---

## 📂 Repository structure

```text
.
├── notebooks/
│   ├── 01_preprocessing.ipynb     # CellBender + QC + integration; reproduces Fig. 1A UMAP
│   └── 02_figure_analysis.ipynb   # Fig. 1 (cell types, T-cell subclusters, markers), fig. S1, Fig. 6A
├── R/
│   └── utils.R                    # shared libraries, helper functions and objects
├── scripts/
│   └── cellbender_jkic.sh         # CellBender remove-background (SLURM array)
├── annotation/
│   └── cluster_annotation.txt     # cluster → cell-type / subcluster labels
├── LICENSE
└── README.md
```

The notebooks are committed **with their executed outputs and figures embedded**, so the
results can be viewed directly on GitHub without re-running anything.

---

## 🔬 Analysis overview

```
FASTQ  ->  Cell Ranger (count)  ->  CellBender (remove-background)
       ->  QC filtering + scDblFinder (doublet removal)
       ->  LogNormalize -> HVG (2000) -> ScaleData -> PCA
       ->  Harmony (batch = sample)  ->  kNN graph -> Louvain clusters -> UMAP
       ->  Cell-type / T-cell subcluster annotation  (Figure 1)
       ->  TCR repertoire integration & Ro60-reactive clonotype analysis
```

- **`01_preprocessing.ipynb`** — ambient-RNA removal (CellBender), quality control and
  doublet detection, normalization and Harmony batch integration, and reproduction of the
  all-cell UMAP (Figure 1A).
- **`02_figure_analysis.ipynb`** — broad cell-type UMAP and marker dot plot (**Fig. 1A–B**),
  T-cell subclustering with marker violins (**Fig. 1C–D**), the T-cell marker heatmap and
  per-sample distribution (**fig. S1**), and the reporter-assay / Ro60-reactive clonotypes on
  the T-cell UMAP (**Fig. 6A**).

---

## 🧰 Requirements

- R ≥ 4.3 with: `Seurat` (v5), `harmony`, `scCustomize`, `scRepertoire`, `tidyverse`,
  `magrittr`, `ggrepel`, `patchwork`, `pheatmap`, `ggpubr`, `ggsci`, `immunarch`,
  `ggVennDiagram`.
- Jupyter with the [IRkernel](https://irkernel.github.io/) to run the notebooks.

To re-run: place `utils.R` and `annotation/cluster_annotation.txt` in the working
directory, set `dir` to the project root that holds `data/` and `output/`, then execute the
notebooks top to bottom.

> Note: the notebooks do **not** write any `.rds`/`.csv` result files; all figures are
> rendered inline so existing analysis outputs are never overwritten.

---

## 📧 Contact

- Corresponding author: **Masaru Takeshita** — takeshita.a5@keio.jp
- Code / analysis: **Jun Inamo** — juninamo@keio.jp

Department of Microbiology and Immunology / Division of Rheumatology,
Keio University School of Medicine, Tokyo, Japan.
