---
title: "Single-Cell RNA-seq Analysis from Scratch"
date: 2026-01-31
description: "A complete guide to analyzing scRNA-seq data"
featured_image: ""
tags: ["single-cell", "RNA-seq", "tutorial"]
---

A practical guide to analyzing single-cell RNA sequencing data, from raw counts to biological insights.

## What You'll Learn

By the end of this tutorial, you'll be able to:

1. Load and preprocess scRNA-seq data
2. Perform quality control and filtering
3. Normalize and integrate datasets
4. Cluster cells and identify cell types
5. Find marker genes and perform differential expression

## Overview

Single-cell RNA sequencing has revolutionized our understanding of cellular heterogeneity. Unlike bulk RNA-seq, which averages gene expression across millions of cells, scRNA-seq measures expression in individual cells.

```
Bulk RNA-seq:    [Cell1 + Cell2 + Cell3 + ...] → Average expression
scRNA-seq:       Cell1 → Expression₁
                 Cell2 → Expression₂
                 Cell3 → Expression₃
```

## Getting Started

### Install Required Packages

**R (Seurat)**
```r
install.packages("Seurat")
install.packages("tidyverse")
install.packages("patchwork")

library(Seurat)
library(tidyverse)
```

**Python (Scanpy)**
```python
pip install scanpy
pip install leidenalg

import scanpy as sc
import numpy as np
import pandas as pd
```

### Download Example Data

We'll use the PBMC 3k dataset from 10x Genomics — a classic benchmark dataset.

**R**
```r
# Download PBMC dataset
pbmc_data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc_data, 
                           project = "pbmc3k",
                           min.cells = 3, 
                           min.features = 200)
```

**Python**
```python
# Download PBMC dataset
adata = sc.datasets.pbmc3k()
adata
```

## Step 1: Quality Control

The first step is filtering out low-quality cells and genes.

### Key QC Metrics

| Metric | What it measures | Typical threshold |
|--------|-----------------|-------------------|
| nFeature | Genes detected per cell | 200 - 5000 |
| nCount | Total UMI counts | 500 - 50000 |
| percent.mt | Mitochondrial reads | < 10-20% |

High mitochondrial percentage often indicates dying cells.

**R**
```r
# Calculate mitochondrial percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filter cells
pbmc <- subset(pbmc, 
               subset = nFeature_RNA > 200 & 
                        nFeature_RNA < 2500 & 
                        percent.mt < 5)
```

**Python**
```python
# Calculate QC metrics
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

# Visualize
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'])

# Filter
sc.pp.filter_cells(adata, min_genes=200)
adata = adata[adata.obs.pct_counts_mt < 5, :]
```

## Coming Soon

This tutorial is being expanded with:
- Normalization methods
- Dimensionality reduction (PCA, UMAP)
- Clustering algorithms
- Cell type annotation
- Differential expression analysis

**Subscribe** to get notified when the full version is available.

---

*Tutorial by [Boxiang Liu](https://twitter.com/boxiangliu) | Last updated: January 2026*
