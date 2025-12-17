# Aging Reconfigures Extracellular Vesicle Signaling Driving Environmental Inflammatory Susceptibility

## Overview
This repository contains the bioinformatics analysis code used in the study  
“Aging Reconfigures Extracellular Vesicle Signaling Driving Environmental Inflammatory Susceptibility”.

All bioinformatics analyses were performed exclusively using the R programming environment.
This repository contains analysis scripts only and does not include any wet-lab experiments,
sample preparation, or primary data generation.

The provided R scripts reproduce the computational analyses reported in the manuscript,
including UMAP-based cell-type visualization and chemokine expression comparisons
under inflammation and environmentally driven inflammatory conditions.

## Data Sources
This study is based entirely on publicly available datasets.

Transcriptomic data were obtained from the Gene Expression Omnibus (GEO).
The GEO accession numbers used in this study are provided in the manuscript.
GEO accession numbers: GSE145926, GSE58682
All analyses represent secondary bioinformatics analyses of existing datasets.
No new experimental data were generated as part of this repository.

## System Requirements

### Operating system
- Windows 11 x64  
(The analysis was tested on this operating system; other modern operating systems
supporting R ≥ 4.2.0 should be compatible.)

### R version
- R ≥ 4.2.0 (tested on R 4.5.2)

### Required R packages
- Seurat
- dplyr
- ggplot2
- patchwork
- scales
- BiocManager
- oligo
- pd.hugene.1.0.st.v1
- hugene10sttranscriptcluster.db
- AnnotationDbi
- stringr
- writexl

No non-standard hardware (e.g., GPU or high-performance computing resources) is required.

## Installation Guide

1. Install R (version ≥ 4.2.0) from CRAN.
2. Install the required R packages:

```r
install.packages(c(
  "Seurat", "dplyr", "ggplot2", "patchwork", "scales",
  "BiocManager", "stringr", "writexl"
))

BiocManager::install(c(
  "oligo",
  "pd.hugene.1.0.st.v1",
  "hugene10sttranscriptcluster.db",
  "AnnotationDbi"
))


