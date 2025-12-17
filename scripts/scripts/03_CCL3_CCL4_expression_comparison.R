suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

## ---------- 0) Input check ----------
stopifnot(all(c("group","sample","CCL3_expr","CCL4_expr","CCL_sum") %in% colnames(am_cells)))

## group factor (Control -> ALI order)
am_cells <- am_cells %>%
  mutate(group = factor(group, levels = c("Control","ALI")))

## ---------- 1) Cell-level long dataframe ----------
cell_long <- am_cells %>%
  select(group, sample, CCL3_expr, CCL4_expr, CCL_sum) %>%
  pivot_longer(
    cols = c(CCL3_expr, CCL4_expr, CCL_sum),
    names_to = "feature",
    values_to = "expr"
  ) %>%
  mutate(
    feature = factor(feature, levels = c("CCL3_expr","CCL4_expr","CCL_sum"))
  )

## Facet нэр
facet_lab <- c(
  CCL3_expr = "CCL3",
  CCL4_expr = "CCL4",
  CCL_sum   = "CCL3+CCL4"
)

## ---------- 2) Colors ----------
cols_fill  <- c(Control="#F6A5A5", ALI="#A8D7FF")  # pastel fill
cols_point <- c(Control="#E15757", ALI="#2E7CCB")  # stronger points

## ---------- 3) Cell-level stats (⚠️ pseudo-replication risk) ----------
p_cell <- cell_long %>%
  group_by(feature) %>%
  summarise(
    n_control = sum(group=="Control"),
    n_ALI     = sum(group=="ALI"),
    p_value   = wilcox.test(expr[group=="Control"], expr[group=="ALI"],
                            exact = FALSE, correct = TRUE)$p.value,
    .groups="drop"
  ) %>%
  mutate(p_value = format.pval(p_value, digits=6, eps=1e-16))

print(p_cell)

## ---------- 4) Sample-level stats (✅ recommended; avoids pseudo-replication) ----------
sample_long <- am_cells %>%
  group_by(group, sample) %>%
  summarise(
    mean_CCL3   = mean(CCL3_expr, na.rm=TRUE),
    mean_CCL4   = mean(CCL4_expr, na.rm=TRUE),
    mean_CCLsum = mean(CCL_sum,   na.rm=TRUE),
    .groups="drop"
  ) %>%
  pivot_longer(cols = starts_with("mean_"), names_to="feature", values_to="value") %>%
  mutate(feature = factor(feature, levels=c("mean_CCL3","mean_CCL4","mean_CCLsum")))

p_sample <- sample_long %>%
  group_by(feature) %>%
  summarise(
    n_control = sum(group=="Control"),
    n_ALI     = sum(group=="ALI"),
    p_value   = wilcox.test(value[group=="Control"], value[group=="ALI"],
                            exact = FALSE, correct = TRUE)$p.value,
    .groups="drop"
  ) %>%
  mutate(
    feature = recode(feature, mean_CCL3="CCL3", mean_CCL4="CCL4", mean_CCLsum="CCL3+CCL4"),
    p_value = format.pval(p_value, digits=6, eps=1e-16)
  )

print(p_sample)

## ---------- 5) Plot: cell-level violin + box + jitter ----------
set.seed(123)  # jitter reproducible

p_violin_cell <- ggplot(cell_long, aes(x=group, y=expr, fill=group)) +
  geom_violin(scale="width", trim=TRUE, color="grey30") +
  geom_boxplot(width=0.15, outlier.shape=NA, alpha=0.90, color="grey30") +
  geom_point(aes(color=group),
             position=position_jitter(width=0.08, height=0),
             size=0.55, alpha=0.45, show.legend=FALSE) +
  facet_wrap(~feature, nrow=1, labeller = as_labeller(facet_lab)) +
  scale_fill_manual(values=cols_fill) +
  scale_color_manual(values=cols_point) +
  ylab("Normalized expression in AM cells") + xlab("") +
  theme_bw(base_size=11) +
  theme(
    legend.position="none",
    strip.background=element_rect(fill="grey95"),
    panel.grid.minor=element_blank()
  )

print(p_violin_cell)

## ---------- 6) Save (optional) ----------
# ggsave("Fig_CellLevel_violin_pastel_points.png", p_violin_cell, width=12, height=4, dpi=300)
# write.csv(p_cell,   "PVAL_cell_level.csv", row.names=FALSE)
# write.csv(p_sample, "PVAL_sample_level.csv", row.names=FALSE) 
