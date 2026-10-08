#!/usr/bin/env Rscript

# Audited J2 K-resolution selection after the disease-blind K=3-12 scan.
# The predeclared automated score is retained, but the final recommendation
# also considers convergence, stability, fragmentation, patient coverage,
# UMI association, separation, and Reference18 composition interpretability.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: 1015_select_j2_k_resolution.R <resolution_dir>")
}

output_dir <- normalizePath(args[[1]], mustWork = TRUE)
figure_dir <- file.path(output_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

summary <- fread(file.path(output_dir, "03_J2_K3_K12_RESOLUTION_SUMMARY.tsv"))
stopifnot(setequal(summary$k, 3:12), nrow(summary) == 10L)

primary_k <- 5L
lower_k <- 3L
higher_k <- 7L

selected <- summary[k %in% c(primary_k, lower_k, higher_k)]
selected[, role := fifelse(
  k == primary_k,
  "primary",
  fifelse(k == lower_k, "lower-resolution comparator", "higher-resolution sensitivity")
)]
selected[, audited_reason := fifelse(
  k == primary_k,
  paste(
    "K=5 retains 20/20 convergence, median ARI 0.988 and median NMI 0.979,",
    "all niches cover all seven patients, the smallest niche is 15.7%,",
    "within-patient UMI eta-squared is 0.040, and Reference18 composition",
    "separates broad B-enriched, endothelial-associated, alveolar-associated,",
    "airway-associated and stromal-associated gradients more clearly than K=3"
  ),
  fifelse(
    k == lower_k,
    paste(
      "K=3 is the automated-score optimum and highly stable (20/20 converged;",
      "median ARI 0.997), but combines several broad composition gradients and",
      "is therefore retained as an under-resolution comparator"
    ),
    paste(
      "K=7 is the nearest eligible higher-resolution sensitivity with 19/20",
      "convergence, all-patient coverage, no niche below 1%, minimum niche",
      "fraction 9.3%, and lower UMI association than K>=8; its median ARI 0.797",
      "requires sensitivity-only interpretation"
    )
  )
)]
selected[, selection_not_based_on := paste(
  "disease separation; disease P values; similarity to L3; frozen aggregate overlap"
)]
selected[, role_order := match(
  role,
  c("primary", "lower-resolution comparator", "higher-resolution sensitivity")
)]
setorder(selected, role_order)
selected[, role_order := NULL]
fwrite(
  selected,
  file.path(output_dir, "04_J2_AUDITED_K_RECOMMENDATION.tsv"),
  sep = "\t"
)

writeLines(
  c(
    "# J2 K分辨率审核结论",
    "",
    "## 推荐",
    "",
    "- **Primary K = 5**。K=5的20/20次拟合均收敛，重复拟合ARI中位数为0.988、NMI中位数为0.979；五个niche均覆盖7位患者，最小niche占15.7%，患者内UMI关联较低（eta-squared=0.040）。与K=3相比，K=5在保持高度稳定的同时，更清楚地区分B富集、内皮、肺泡、气道和基质相关的连续混合组成。",
    "- **Lower-resolution comparator = K=3**。K=3是预设自动综合分数最高的方案，技术稳定性最好，但把若干可解释的Reference18组成轴合并，因此作为欠分辨对照而非主分辨率。",
    "- **Higher-resolution sensitivity = K=7**。K=7仍覆盖全部患者、无<1%小簇，且UMI关联低于K>=8；但重复拟合稳定性下降（ARI中位数0.797），仅用于高分辨率敏感性。",
    "",
    "## 选择限制",
    "",
    "疾病标签、疾病分离、frozen aggregate以及L3 niche编号均未参与邻域构建、聚类、K选择或命名。K=5的选择不是为了获得疾病差异，也不表示五个niche是离散细胞类型。"
  ),
  file.path(output_dir, "05_J2_K_RECOMMENDATION_ZH.md"),
  useBytes = TRUE
)

plot_data <- melt(
  summary,
  id.vars = "k",
  measure.vars = c(
    "converged_fraction",
    "repeat_ari_median",
    "repeat_nmi_median",
    "min_niche_fraction",
    "within_patient_umi_eta_squared",
    "betweenss_fraction",
    "interpretability_index"
  ),
  variable.name = "metric",
  value.name = "value"
)
plot_data[, selected_role := factor(
  fifelse(
    k == primary_k,
    "Primary",
    fifelse(k == lower_k, "Lower comparator", fifelse(k == higher_k, "Higher sensitivity", "Other"))
  ),
  levels = c("Primary", "Lower comparator", "Higher sensitivity", "Other")
)]

plot <- ggplot(plot_data, aes(k, value, group = metric)) +
  geom_line(color = "grey45", linewidth = 0.6) +
  geom_point(aes(fill = selected_role), shape = 21, color = "black", size = 2.8) +
  facet_wrap(~metric, scales = "free_y", ncol = 2) +
  scale_x_continuous(breaks = 3:12) +
  scale_fill_manual(values = c(
    "Primary" = "#D55E00",
    "Lower comparator" = "#0072B2",
    "Higher sensitivity" = "#009E73",
    "Other" = "white"
  )) +
  labs(
    title = "J2 unbiased niche-resolution audit",
    subtitle = "K=5 primary; K=3 lower-resolution comparator; K=7 higher-resolution sensitivity",
    x = "K",
    y = "Metric value",
    fill = "Audited role"
  ) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())

ggsave(
  file.path(figure_dir, "01_J2_K_RESOLUTION_AUDIT_EN.pdf"),
  plot,
  width = 11,
  height = 10,
  device = cairo_pdf,
  limitsize = FALSE
)
ggsave(
  file.path(figure_dir, "01_J2_K_RESOLUTION_AUDIT_EN.png"),
  plot,
  width = 11,
  height = 10,
  dpi = 320,
  limitsize = FALSE
)

if (requireNamespace("sessioninfo", quietly = TRUE)) {
  capture.output(
    sessioninfo::session_info(),
    file = file.path(output_dir, "K_SELECTION_SESSION_INFO.txt")
  )
} else {
  capture.output(sessionInfo(), file = file.path(output_dir, "K_SELECTION_SESSION_INFO.txt"))
}
