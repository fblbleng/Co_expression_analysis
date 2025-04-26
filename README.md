# Re-analysis of RNA-Seq Dataset Using UMAP, k-means Clustering, and Co-expression Network Analysis

This project presents a re-analysis of the publicly available **GSE98212** RNA-Seq dataset, combining advanced techniques such as UMAP, k-means clustering, and weighted gene co-expression network analysis (WGCNA).

---

## Table of Contents

- [Background](#background)
- [Objectives](#objectives)
- [Workflow Summary](#workflow-summary)
- [How to Run](#how-to-run)
- [Requirements](#requirements)
- [License](#license)

---

## Background

Whole-transcriptome technologies like microarrays and RNA-Seq provide rich gene expression profiles across tissues and timepoints. Traditional analyses often focus on identifying differentially expressed genes (DEGs) between predefined groups. However, such approaches typically overlook the complex interconnections between genes and the heterogeneity within sample groups.

Identifying genes with correlated expression patterns can reveal shared biological functions or regulatory relationships. Similarly, investigating sample-to-sample heterogeneity is crucial for detecting outliers or hidden subgroups, which, if unaddressed, can bias differential gene expression analyses.

---

## Objectives

- Capture hidden structures within the dataset using **UMAP** for dimensionality reduction and **k-means clustering**.
- Uncover sample heterogeneity and potential subgroups missed by standard group-based comparisons.
- Integrate clustering results into a **weighted gene co-expression network analysis (WGCNA)** to explore gene networks associated with sample groups.

---

## Workflow Summary

### Data Loading
- Import the GSE98212 RNA-Seq dataset.

### Dimensionality Reduction & Clustering
- Apply **UMAP** to project high-dimensional expression data into a low-dimensional space.
- Use **k-means clustering** to identify hidden sample subgroups based on the UMAP projection.

### Co-expression Network Construction
- Perform **WGCNA** to build a weighted gene co-expression network.
- Associate identified gene modules with the sample clusters and highlight key regulatory genes.

---

## How to Run

### Clone the Repository

```bash
git clone https://github.com/yourusername/demo_R.git
cd demo_R
```

### Set Up the Environment

Install the required R packages:

```r
install.packages(c("umap", "WGCNA", "ggplot2", "cluster", "dplyr", "Seurat"))
```

You can also use a `renv` or `conda` environment if preferred.

### Run the Analysis

Execute the main R script:

```r
source("main_analysis.R")
```

This script will perform data loading, clustering, dimensionality reduction, and co-expression network analysis.

### View the Results

Outputs such as plots, clustering assignments, and network graphs will be saved in the `/results` directory.

---

## Requirements

- R version >= 4.2.0

### Required Packages

- `umap`
- `WGCNA`
- `Seurat`
- `ggplot2`
- `dplyr`
- `cluster`

---






