## =========================
## ALI scRNA-seq analysis (Seurat) — reviewer-ready pipeline
## Generates: processed Seurat object + UMAP plot + CCL3/CCL4 stats
## =========================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)

## -------------------------
## 0) User settings (EDIT HERE)
## -------------------------
# Project root: энэ скрипт байрлаж буй фолдер дотроос ажиллуулна гэж үзье.
# Хэрвээ өөр газраас ажиллуулах бол RStudio дээр "Session > Set Working Directory > To Source File Location" сонго.

h5_dir      <- "data_h5"          # H5 файлын фолдер
out_rds_dir <- file.path("results", "rds")
out_fig_dir <- file.path("results", "figs")
out_qc_dir  <- file.path("results", "qc_csv")

dir.create(out_rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_qc_dir,  recursive = TRUE, showWarnings = FALSE)

group_label <- "ALI"

# QC thresholds (EDIT if needed)
min_features <- 200
min_counts   <- 500
max_mt       <- 25

# Downstream settings
n_hvf     <- 3000
npcs      <- 40
dims_use  <- 1:30
resolution <- 0.4

## -------------------------
## 1) Detect input H5 files
## -------------------------
h5_files <- list.files(h5_dir, pattern = "_filtered_feature_bc_matrix\\.h5$", full.names = TRUE)
stopifnot(length(h5_files) > 0)

message("Found H5 files: ", length(h5_files))

## sample_id parse helper
parse_sample_id <- function(f) {
  base <- basename(f)
  sid  <- sub("_filtered.*", "", sub("^GSM\\d+_", "", base))
  sid
}

## -------------------------
## 2) Read each H5 -> QC filter -> save filtered RDS (counts-only)
## -------------------------
qc_one <- function(fpath) {
  sid <- parse_sample_id(fpath)

  mat <- Read10X_h5(fpath)
  so  <- CreateSeuratObject(mat, project = sid, min.features = min_features)
  so$sample_id <- sid
  so$group     <- group_label

  # QC
  so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
  so <- subset(
    so,
    subset = nFeature_RNA >= min_features & nCount_RNA >= min_counts & percent.mt < max_mt
  )

  qc <- data.frame(
    sample_id  = sid,
    cells      = ncol(so),
    features   = nrow(so),
    mt_median  = median(so$percent.mt)
  )
  write.csv(qc, file.path(out_qc_dir, paste0("QC_", sid, ".csv")), row.names = FALSE)

  # save filtered counts-only object
  so <- DietSeurat(so, counts = TRUE, data = FALSE, scale.data = FALSE)
  saveRDS(so, file.path(out_rds_dir, paste0("ALI_filt_", sid, ".rds")))

  rm(mat, so); gc()
  invisible(TRUE)
}

invisible(lapply(h5_files, qc_one))
message("✓ Per-sample QC + filtered RDS saved to: ", out_rds_dir)

## -------------------------
## 3) Merge filtered samples (counts-only)
## -------------------------
filt_rds <- list.files(out_rds_dir, pattern = "^ALI_filt_.*\\.rds$", full.names = TRUE)
stopifnot(length(filt_rds) > 0)

first <- TRUE
for (r in filt_rds) {
  so <- readRDS(r)
  so <- DietSeurat(so, counts = TRUE, data = FALSE, scale.data = FALSE)

  if (first) {
    ali <- so; first <- FALSE
  } else {
    ali <- merge(ali, y = so, merge.data = FALSE)
  }
  rm(so); gc()
}

saveRDS(ali, file.path(out_rds_dir, "ALI_merged_countOnly.rds"))
message("✓ Saved merged object: ", file.path(out_rds_dir, "ALI_merged_countOnly.rds"))

## -------------------------
## 4) Standard Seurat processing (Normalize -> HVF -> PCA -> UMAP)
## -------------------------
ali <- NormalizeData(ali, verbose = FALSE)
ali <- FindVariableFeatures(ali, selection.method = "vst", nfeatures = n_hvf, verbose = FALSE)

vf  <- VariableFeatures(ali)
ali <- ScaleData(ali, features = vf, verbose = FALSE)
ali <- RunPCA(ali, features = vf, npcs = npcs, verbose = FALSE)
ali <- FindNeighbors(ali, dims = dims_use, verbose = FALSE)
ali <- FindClusters(ali, resolution = resolution, verbose = FALSE)
ali <- RunUMAP(ali, dims = dims_use, verbose = FALSE)

saveRDS(ali, file.path(out_rds_dir, "ALI_proc.rds"))
message("✓ Saved processed object: ", file.path(out_rds_dir, "ALI_proc.rds"))

## -------------------------
## 5) UMAP plot: 4 major immune lineages
##    (Requires: ali$celltype column exists OR Idents(ali) already set to cell type)
## -------------------------
# EDIT HERE: if your metadata column name is different
celltype_col <- "celltype"

if (celltype_col %in% colnames(ali@meta.data)) {
  Idents(ali) <- ali@meta.data[[celltype_col]]
} else {
  message("NOTE: meta column '", celltype_col, "' not found. Using current Idents(ali).")
}

keep_types <- c("Alveolar macrophage", "T cell", "B cell", "Neutrophil")
ali_4 <- subset(ali, idents = intersect(levels(ali), keep_types))
ali_4 <- SetIdent(ali_4, value = Idents(ali_4))

# stable order
ali_4@meta.data$celltype_plot <- factor(Idents(ali_4), levels = keep_types)

cols4 <- c(
  "Alveolar macrophage" = "#A7C7E7",
  "T cell"              = "#F4B6C2",
  "B cell"              = "#FFD1A9",
  "Neutrophil"          = "#98DF8A"
)

p_umap <- DimPlot(
  ali_4,
  reduction = "umap",
  group.by  = "celltype_plot",
  cols      = cols4,
  pt.size   = 0.35
) +
  ggtitle("ALI UMAP — Alveolar macrophage / T cell / B cell / Neutrophil") +
  xlab("umap_1") + ylab("umap_2") +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank()
  )

ggsave(file.path(out_fig_dir, "ALI_UMAP_4celltypes.png"), p_umap, width = 8.6, height = 5.4, dpi = 180)

## -------------------------
## 6) CCL3/CCL4 co-expression stats in Alveolar macrophages
## -------------------------
ali_mac <- subset(ali_4, idents = "Alveolar macrophage")
g1 <- grep("(?i)^ccl3$", rownames(ali_mac), value = TRUE)[1]
g2 <- grep("(?i)^ccl4$", rownames(ali_mac), value = TRUE)[1]
stopifnot(!is.na(g1), !is.na(g2))

expr_mac <- FetchData(ali_mac, vars = c(g1, g2))
expr_mac$CCL3_4_coexp <- expr_mac[[g1]] * expr_mac[[g2]]

summary_stats <- data.frame(
  Gene = c("CCL3", "CCL4", "CCL3×CCL4 (co-expression)"),
  Mean = c(mean(expr_mac[[g1]], na.rm = TRUE),
           mean(expr_mac[[g2]], na.rm = TRUE),
           mean(expr_mac$CCL3_4_coexp, na.rm = TRUE)),
  Median = c(median(expr_mac[[g1]], na.rm = TRUE),
             median(expr_mac[[g2]], na.rm = TRUE),
             median(expr_mac$CCL3_4_coexp, na.rm = TRUE)),
  SD = c(sd(expr_mac[[g1]], na.rm = TRUE),
         sd(expr_mac[[g2]], na.rm = TRUE),
         sd(expr_mac$CCL3_4_coexp, na.rm = TRUE))
)

write.csv(summary_stats, file.path("results", "CCL3_CCL4_coexpression_stats_ALI_macrophage.csv"),
          row.names = FALSE)
