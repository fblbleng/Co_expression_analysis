# hdWGCNA Single-cell Tutorial

## Tutorial Overview

This tutorial provides a comprehensive walkthrough of using **hdWGCNA** (high-dimensional Weighted Gene Co-expression Network Analysis) for single-cell RNA-seq data, specifically demonstrating the analysis pipeline using a preprocessed single-nucleus RNA-seq (snRNA-seq) dataset of human cortical samples.

### [Link to Tutorial](https://smorabit.github.io/hdWGCNA/index.html)

## Key Steps Summarized

### 1. **Dataset Preparation and Setup**

* Load and preprocess single-cell data using Seurat.
* Initialize Seurat object specifically for hdWGCNA, selecting genes based on their expression fraction across cells.

### 2. **Metacell Construction**

* Generate metacells by aggregating similar single cells to reduce sparsity.
* Normalize aggregated expression matrices for improved network analysis performance.

### 3. **Co-expression Network Analysis**

* Focus analysis on specific cell types, like inhibitory neurons (INH).
* Select optimal soft-power threshold to construct robust gene-gene correlation networks.
* Generate and visualize dendrograms depicting co-expression modules.

### 4. **Module Eigengenes and Connectivity**

* Calculate harmonized module eigengenes (hMEs) to summarize expression patterns within modules.
* Compute eigengene-based connectivity (kME) to identify highly connected "hub genes" within each module.

### 5. **Hub Gene Signature Scoring**

* Compute gene signature scores for the most influential genes using UCell or Seurat scoring methods.

### 6. **Visualization Tools**

* Visualize module expression across cells using feature plots and radar plots.
* Utilize custom visualization with Seurat's built-in functions for detailed and customized analysis.

## Final Outputs

* Saved Seurat object containing all analysis results (`hdWGCNA_object.rds`).
* Visual and quantitative summaries facilitating the identification and interpretation of biologically relevant gene networks.
