## ================= 0) Пакетууд суулгах/ачаалах =================
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")

req_bioc <- c("oligo","pd.hugene.1.0.st.v1","hugene10sttranscriptcluster.db","AnnotationDbi","Biobase")
for (p in req_bioc) if (!requireNamespace(p, quietly=TRUE)) BiocManager::install(p, ask=FALSE, update=FALSE)

req_cran <- c("dplyr","stringr","writexl","ggplot2","ggpubr","tidyr")
for (p in req_cran) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)

library(oligo)
library(Biobase)
library(AnnotationDbi)
library(hugene10sttranscriptcluster.db)
library(dplyr)
library(stringr)
library(writexl)
library(ggplot2)
library(ggpubr)
library(tidyr)

## ================= 1) Файлын зам, CEL хайх, meta хийх =================
base_dir <- "D:/bio-informatic analysis for CCL3_4/ozone inhalation/Ozone Inhalation group"

cel_files <- list.files(base_dir, pattern="\\.CEL$", ignore.case=TRUE,
                        recursive=TRUE, full.names=TRUE)
stopifnot(length(cel_files) > 0)

meta <- tibble(
  file   = cel_files,
  sample = tools::file_path_sans_ext(basename(file)),
  group  = ifelse(str_detect(file, regex("Clean\\s*air", ignore_case=TRUE)), "CleanAir", "Ozone"),
  dose_ppb = case_when(
    str_detect(sample, "_100$") ~ 100L,
    str_detect(sample, "_200$") ~ 200L,
    TRUE ~ NA_integer_
  )
)

## ================= 2) RMA нормчлол (log2) =================
raw  <- read.celfiles(cel_files)
eset <- rma(raw)
expr <- Biobase::exprs(eset)   # probeset x sample (log2)

## ================= 3) Аннотаци (PROBEID -> SYMBOL) =================
map_all <- AnnotationDbi::select(hugene10sttranscriptcluster.db,
                                 keys = rownames(expr),
                                 columns = "SYMBOL",
                                 keytype = "PROBEID") |>
  dplyr::filter(!is.na(SYMBOL)) |>
  dplyr::distinct(PROBEID, SYMBOL)

## ================= 4) CCL3/CCL4 gene-level (probeset mean) =================
get_gene_means <- function(g, map_df, expr_mat){
  pr <- map_df$PROBEID[map_df$SYMBOL == g]
  pr <- intersect(pr, rownames(expr_mat))
  if (length(pr) == 0) return(rep(NA_real_, ncol(expr_mat)))
  colMeans(expr_mat[pr, , drop=FALSE], na.rm=TRUE)
}

mat_cc <- cbind(
  CCL3 = get_gene_means("CCL3", map_all, expr),
  CCL4 = get_gene_means("CCL4", map_all, expr)
)

gene_expr <- as.data.frame(mat_cc)
gene_expr$sample <- tools::file_path_sans_ext(colnames(expr))
gene_expr$CCL3 <- as.numeric(gene_expr$CCL3)
gene_expr$CCL4 <- as.numeric(gene_expr$CCL4)

## ================= 5) Z-score + co-expression =================
zmat <- scale(gene_expr[, c("CCL3","CCL4")]) |> as.data.frame()
colnames(zmat) <- c("zCCL3","zCCL4")
gene_expr$zCCL3 <- zmat$zCCL3
gene_expr$zCCL4 <- zmat$zCCL4
gene_expr$CoExpr_Zmean <- rowMeans(zmat, na.rm=TRUE)

## ================= 6) Source data (meta join) =================
source_df <- meta |>
  dplyr::left_join(gene_expr, by="sample") |>
  dplyr::mutate(
    CCL3_log2 = CCL3,
    CCL4_log2 = CCL4
  ) |>
  dplyr::select(sample, group, dose_ppb,
                CCL3_log2, CCL4_log2, zCCL3, zCCL4, CoExpr_Zmean, file) |>
  dplyr::arrange(group, dose_ppb, sample)

## ================= 7) (Сонголт) AM inflammation signature score =================
am_inflam_genes <- c("CCL3","CCL4","IL1B","IL6","TNF","CXCL8","CXCL2",
                     "CCL2","NFKBIA","PTGS2","S100A8","S100A9","MMP9")
am_inflam_genes <- ifelse(am_inflam_genes == "IL8", "CXCL8", am_inflam_genes)

genes_found   <- intersect(am_inflam_genes, unique(map_all$SYMBOL))
genes_missing <- setdiff(am_inflam_genes, genes_found)
if (length(genes_found) == 0) stop("AM gene set олдсонгүй (annotation дээр байхгүй).")

am_mat <- sapply(genes_found, get_gene_means, map_df=map_all, expr_mat=expr)  # samples x genes
am_df  <- as.data.frame(am_mat)
am_df$sample <- tools::file_path_sans_ext(colnames(expr))

am_z <- scale(am_df[, genes_found, drop=FALSE]) |> as.data.frame()
colnames(am_z) <- paste0("z.", genes_found)
AMInflam_Score <- rowMeans(am_z, na.rm=TRUE)

source_df <- source_df |>
  dplyr::mutate(
    AMInflam_Score  = AMInflam_Score,
    AMInflam_nGenes = length(genes_found)
  )

## ================= 8) Excel + CSV экспорт =================
out_xlsx1 <- file.path(base_dir, "source_CCL3_CCL4_coexpression.xlsx")
writexl::write_xlsx(list(SourceData = source_df, Meta = meta), path = out_xlsx1)

out_csv1 <- file.path(base_dir, "source_CCL3_CCL4_coexpression.csv")
write.csv(source_df, file = out_csv1, row.names=FALSE)

out_xlsx2 <- file.path(base_dir, "source_CCL3_CCL4_AMinflammation.xlsx")
writexl::write_xlsx(
  list(
    All              = source_df,
    CleanAir         = dplyr::filter(source_df, group=="CleanAir"),
    Ozone            = dplyr::filter(source_df, group=="Ozone"),
    Meta             = meta,
    AM_genes_found   = data.frame(gene = genes_found),
    AM_genes_missing = data.frame(gene = genes_missing),
    AMgene_zscores   = cbind(sample = am_df$sample, am_z)
  ),
  path = out_xlsx2
)

message("DONE → ", out_xlsx1)
message("DONE → ", out_xlsx2)

## ================= 9) Богино шалгалт/статистик (сонголт) =================
print(table(source_df$group))
print(utils::head(source_df, 4))

if ("AMInflam_Score" %in% names(source_df)) {
  print(
    source_df %>%
      dplyr::group_by(group) %>%
      dplyr::summarise(mean_AM = mean(AMInflam_Score, na.rm=TRUE),
                       sd_AM   = sd(AMInflam_Score, na.rm=TRUE),
                       n       = dplyr::n())
  )
  print(t.test(AMInflam_Score ~ group, data = source_df))
}

## ================= 10) Figure 1: CCL3 / CCL4 / CCL3+CCL4 (univariate) =================
# ---- 10.1) Figure-д хэрэгтэй дата ----
fig_df <- source_df %>%
  transmute(
    sample,
    group = factor(group, levels = c("CleanAir","Ozone")),
    CCL3 = CCL3_log2,
    CCL4 = CCL4_log2,
    CCL3_CCL4 = CoExpr_Zmean      # таны co-expression (mean z)
  ) %>%
  filter(!is.na(CCL3), !is.na(CCL4), !is.na(CCL3_CCL4))

# ---- 10.2) long формат (3 панель) ----
df_all <- fig_df %>%
  pivot_longer(cols = c(CCL3, CCL4, CCL3_CCL4),
               names_to = "feature", values_to = "value") %>%
  mutate(
    feature = recode(feature,
                     CCL3 = "CCL3",
                     CCL4 = "CCL4",
                     CCL3_CCL4 = "CCL3+CCL4"),
    feature = factor(feature, levels = c("CCL3","CCL4","CCL3+CCL4"))
  )

# ---- 10.3) test method сонгох ----
# method <- "t.test"
method <- "wilcox.test"

# ---- 10.4) plot ----
p1 <- ggboxplot(
    df_all, x="group", y="value",
    fill="group", palette="Pastel1",
    add="jitter", facet.by="feature", nrow=1
  ) +
  stat_compare_means(
    method = method,
    label  = "p.signif",
    comparisons = list(c("CleanAir","Ozone"))
  ) +
  labs(
    title = "Figure 1. CCL3, CCL4, and CCL3+CCL4 (Control vs Ozone)",
    x = NULL, y = "Expression (log2 normalized) / co-expression score"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")

print(p1)

# ---- 10.5) Save figure ----
fig_out <- file.path(base_dir, "Figure1_CCL3_CCL4_CoExpr.png")
ggsave(fig_out, p1, width = 8.5, height = 3.6, dpi = 600)
message("DONE → ", fig_out)

# ---- 10.6) Save figure source data (optional) ----
fig_src <- file.path(base_dir, "Figure1_source_data_CCL3_CCL4_CoExpr.csv")
write.csv(fig_df, fig_src, row.names = FALSE)
message("DONE → ", fig_src)

# ---- 10.7) Print exact p-values (panel-wise) ----
cat("\n=== Figure 1: exact tests ===\n")
cat("\n[CCL3]\n")
print(if (method=="t.test") t.test(CCL3 ~ group, data=fig_df) else wilcox.test(CCL3 ~ group, data=fig_df, exact=FALSE))
cat("\n[CCL4]\n")
print(if (method=="t.test") t.test(CCL4 ~ group, data=fig_df) else wilcox.test(CCL4 ~ group, data=fig_df, exact=FALSE))
cat("\n[CCL3+CCL4]\n")
print(if (method=="t.test") t.test(CCL3_CCL4 ~ group, data=fig_df) else wilcox.test(CCL3_CCL4 ~ group, data=fig_df, exact=FALSE))
}
