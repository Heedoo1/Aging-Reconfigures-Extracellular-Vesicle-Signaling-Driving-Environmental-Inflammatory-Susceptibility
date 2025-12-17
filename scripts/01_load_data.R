## ===== CONTROL (Healthy) UMAP pipeline: display UMAP + CCL3/CCL4 4-panel =====
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(scales)

## ------------------------------------------------------------
## 0) INPUT: obj (Seurat), metadata: celltype_coarse
## ------------------------------------------------------------
stopifnot(inherits(obj, "Seurat"))
stopifnot("celltype_coarse" %in% colnames(obj[[]]))

# ашиглах reduction ("usw" байвал түүнийг, үгүй бол "umap")
red <- if ("usw" %in% names(obj@reductions)) "usw" else "umap"
emb <- Embeddings(obj, red)[, 1:2, drop = FALSE]
colnames(emb) <- c("umap1","umap2")

plot_df <- data.frame(
  cell     = rownames(emb),
  umap1    = emb[,1],
  umap2    = emb[,2],
  celltype = obj$celltype_coarse
)

# (сонголт) зөвхөн 4 coarse бүлэг үлдээнэ
keep_levels <- c("Alveolar macrophage","T cell","B cell","Neutrophil")
plot_df <- plot_df %>% filter(celltype %in% keep_levels)
plot_df$celltype <- factor(plot_df$celltype, levels = keep_levels)

## ------------------------------------------------------------
## 1) (OPTION) “харагдах байдлын” clean: core + buffer (cross-mixing багасгах)
##    - AM: бүх цэг хэвээр
##    - T/B/Neu: core дотор + бусад төвүүдэд хэт ойр бол хасна
## ------------------------------------------------------------
cts <- levels(plot_df$celltype)

# төвүүд (median)
cent <- lapply(cts, function(ct){
  ix <- plot_df$celltype == ct
  c(median(plot_df$umap1[ix], na.rm=TRUE),
    median(plot_df$umap2[ix], na.rm=TRUE))
}); names(cent) <- cts

diag_len <- sqrt(diff(range(plot_df$umap1, na.rm=TRUE))^2 +
                 diff(range(plot_df$umap2, na.rm=TRUE))^2)

q_core <- 0.80                 # T/B/Neu core quantile
buffer <- 0.03 * diag_len      # 3% diagonal

# өөрийн төв хүртэлх зай (d_self), хамгийн ойр өөр төв хүртэлх зай (d_alt)
d_self <- numeric(nrow(plot_df))
d_alt  <- numeric(nrow(plot_df))

for(i in seq_len(nrow(plot_df))){
  ct <- as.character(plot_df$celltype[i])
  p  <- c(plot_df$umap1[i], plot_df$umap2[i])
  d0 <- sqrt(sum((p - cent[[ct]])^2))
  d1 <- sapply(setdiff(cts, ct), function(k){
    sqrt(sum((p - cent[[k]])^2))
  })
  d_self[i] <- d0
  d_alt[i]  <- min(d1)
}

plot_df$d_self <- d_self
plot_df$d_alt  <- d_alt

thr_by_ct <- sapply(cts, function(ct){
  quantile(plot_df$d_self[plot_df$celltype == ct], probs = q_core, na.rm = TRUE)
})

keep_clean <- with(plot_df, {
  isAM <- celltype == "Alveolar macrophage"
  core <- d_self <= thr_by_ct[as.character(celltype)]
  notmix <- (d_self + buffer) <= d_alt
  isAM | (core & notmix)
})

df_clean <- plot_df[keep_clean, c("cell","umap1","umap2","celltype")]

## ------------------------------------------------------------
## 2) (OPTION) “shape” тохируулах: shrink + (OPTION) subsample (visibility ratio)
##    AM хэвээр; T/B бага; Neu хамгийн жижиг/цөөн харагдана
## ------------------------------------------------------------
# shrink factors (cluster-ийн “тархалтыг” жижигрүүлэх)
shrink_by_ct <- c(
  "Alveolar macrophage" = 1.00,
  "T cell"              = 0.90,
  "B cell"              = 0.85,
  "Neutrophil"          = 0.65
)

df2 <- df_clean
for(ct in cts){
  ix <- df2$celltype == ct
  if(!any(ix)) next
  s  <- shrink_by_ct[[ct]]
  c0 <- cent[[ct]]
  df2$umap1[ix] <- c0[1] + s*(df2$umap1[ix] - c0[1])
  df2$umap2[ix] <- c0[2] + s*(df2$umap2[ix] - c0[2])
}

# харагдах цэгийн ratio (дүрслэлд цэг хэт их бол)
vis_ratio <- c(
  "Alveolar macrophage" = 1.00,
  "T cell"              = 0.75,
  "B cell"              = 0.55,
  "Neutrophil"          = 0.35
)

set.seed(123)
df3 <- bind_rows(lapply(cts, function(ct){
  sub <- df2 %>% filter(celltype == ct)
  n_keep <- ceiling(nrow(sub) * vis_ratio[[ct]])
  if(n_keep < nrow(sub)) sub <- sub[sample(nrow(sub), n_keep), , drop=FALSE]
  sub
}))

## ------------------------------------------------------------
## 3) Control UMAP plot (энэ бол “харагдаж буй UMAP”)
## ------------------------------------------------------------
cols_ct <- c(
  "Alveolar macrophage"="#9BBBE6",
  "T cell"              ="#F3A5B5",
  "B cell"              ="#F3B57A",
  "Neutrophil"          ="#9BE49B"
)

p_control_umap <- ggplot(df3, aes(umap1, umap2)) +
  geom_point(aes(color = celltype), size = 0.40, alpha = 0.95) +
  scale_color_manual(values = cols_ct,
                     guide = guide_legend(title = NULL, override.aes = list(size=4))) +
  labs(title = "Healthy (Control) BALF UMAP — cleaned display",
       x = "umap1", y = "umap2") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face="bold"),
        legend.position="right")
print(p_control_umap)

## ------------------------------------------------------------
## 4) CCL3 / CCL4 4-panel (яг df3 дээрх байрлал ашиглана)
## ------------------------------------------------------------
genes <- c("CCL3","CCL4")
stopifnot(all(genes %in% rownames(obj)))
expr <- FetchData(obj, vars = genes)
expr$cell <- rownames(expr)

dfp <- merge(df3, expr, by = "cell", all.x = TRUE)
dfp$CCL3[is.na(dfp$CCL3)] <- 0
dfp$CCL4[is.na(dfp$CCL4)] <- 0

r01 <- function(x){
  rng <- range(x, na.rm=TRUE)
  if(!is.finite(rng[1]) || !is.finite(rng[2]) || diff(rng)==0) return(rep(0, length(x)))
  (x - rng[1]) / diff(rng)
}

dfp$ccl3_n <- r01(dfp$CCL3)
dfp$ccl4_n <- r01(dfp$CCL4)
dfp$sum01  <- r01(dfp$ccl3_n + dfp$ccl4_n)
dfp$prod01 <- r01(dfp$ccl3_n * dfp$ccl4_n)

eps <- 0
is0_ccl3 <- dfp$ccl3_n <= eps
is0_ccl4 <- dfp$ccl4_n <= eps
is0_sum  <- dfp$sum01  <= eps
is0_prod <- dfp$prod01 <= eps

# palettes
col_grey <- "#D9D9D9"
col_ccl3 <- c("#e8f1fb","#cfe0f6","#9fc2ea","#6aa4df","#3182bd")
col_ccl4 <- c("#e6f5e6","#c9ecc9","#9fe09f","#72d272","#31a354")
col_sum  <- c("#f7fbff","#deebf7","#9ecae1","#6baed6","#2171b5")
col_prod <- c("#FAE4F3","#F7CBEA","#F2ACDF","#E58BDD","#D36CD5")  # LOW→HIGH light→dark

base_theme <- theme_classic(base_size=12) +
  theme(legend.position="right", plot.title=element_text(face="bold"))

p1 <- ggplot(dfp, aes(umap1, umap2)) +
  geom_point(data=dfp[is0_ccl3,], color=col_grey, size=0.35, alpha=0.85) +
  geom_point(data=dfp[!is0_ccl3,], aes(color=ccl3_n), size=0.35, alpha=0.95) +
  scale_color_gradientn(colors=col_ccl3, limits=c(0,1), na.value=col_grey) +
  labs(title="CCL3 expression", x="umap1", y="umap2", color=NULL) + base_theme

p2 <- ggplot(dfp, aes(umap1, umap2)) +
  geom_point(data=dfp[is0_ccl4,], color=col_grey, size=0.35, alpha=0.85) +
  geom_point(data=dfp[!is0_ccl4,], aes(color=ccl4_n), size=0.35, alpha=0.95) +
  scale_color_gradientn(colors=col_ccl4, limits=c(0,1), na.value=col_grey) +
  labs(title="CCL4 expression", x="umap1", y="umap2", color=NULL) + base_theme

p3 <- ggplot(dfp, aes(umap1, umap2)) +
  geom_point(data=dfp[is0_sum,], color=col_grey, size=0.35, alpha=0.85) +
  geom_point(data=dfp[!is0_sum,], aes(color=sum01), size=0.35, alpha=0.95) +
  scale_color_gradientn(colors=col_sum, limits=c(0,1), na.value=col_grey) +
  labs(title="CCL3/4 combined (sum, 0–1)", x="umap1", y="umap2", color=NULL) + base_theme

p4 <- ggplot(dfp, aes(umap1, umap2)) +
  geom_point(data=dfp[is0_prod,], color=col_grey, size=0.35, alpha=0.85) +
  geom_point(data=dfp[!is0_prod,], aes(color=prod01, alpha=prod01), size=0.35) +
  scale_color_gradientn(colors=col_prod, limits=c(0,1), na.value=col_grey) +
  scale_alpha(range=c(0.25,1.00), limits=c(0,1), guide="none") +
  guides(color = guide_colorbar(title=NULL, ticks=TRUE)) +
  labs(title="CCL3 & CCL4 co-expression (product, 0–1)", x="umap1", y="umap2") +
  base_theme

p_4panel <- (p1 | p2) / (p3 | p4)
print(p_4panel)

# ggsave("Control_UMAP_clean.png", p_control_umap, width=9, height=5.5, dpi=300)
# ggsave("Control_UMAP_CCL3_CCL4_4panel.png", p_4panel, width=11, height=7.5, dpi=300)
