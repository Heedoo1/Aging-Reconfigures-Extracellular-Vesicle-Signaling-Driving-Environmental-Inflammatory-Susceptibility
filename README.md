# Aging Reconfigures Extracellular Vesicle Signaling Driving Environmental Inflammatory Susceptibility

## Overview
This repository contains the bioinformatics analysis code used in the study  
“Aging Reconfigures Extracellular Vesicle Signaling Driving Environmental Inflammatory Susceptibility.”

All bioinformatics  analyses were performed exclusively using the R programming environment.  
No wet-lab experiments, sample preparation, or data generation were conducted as part of this repository.

The provided R scripts reproduce the computational analyses used to investigate age-associated alterations in extracellular vesicle (EV)–mediated signaling under environmentally driven inflammatory conditions.

---

## Data Sources
This study is based entirely on publicly available datasets.

- Transcriptomic data were obtained from the Gene Expression Omnibus (GEO).
- GEO accession numbers are provided in the manuscript.

All analyses in this repository represent secondary bioinformatics analyses of existing datasets.  
No new experimental data were generated.

---

## Computational Environment

All analyses were conducted exclusively in the R programming environment**.

### R version
- R ≥ 4.2.0

### Required R packages
```r
BiocManager
oligo
pd.hugene.1.0.st.v1
hugene10sttranscriptcluster.db
AnnotationDbi
dplyr
stringr
writexl
