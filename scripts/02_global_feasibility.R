#!/usr/bin/env Rscript
root <- normalizePath(file.path(getwd(), "AD_CD33_reanalysis"))
res_dir <- file.path(root, "results", "feasibility")
fig_dir <- file.path(root, "figures", "feasibility")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
x <- read.csv(file.path(res_dir, "all_probe_contrasts.csv"), check.names = FALSE)

x$significant <- x$adj.P.Val < 0.05
x$moderate <- x$adj.P.Val < 0.05 & abs(x$logFC) >= 0.25
counts <- aggregate(cbind(significant, moderate) ~ series + tissue + contrast, x, sum)
write.csv(counts, file.path(res_dir, "differential_signal_counts.csv"), row.names = FALSE)

comparison_rows <- list()
reversal_rows <- list()
replicated_rows <- list()
k <- 1
for (tissue in unique(x$tissue)) {
  for (contrast in c("Disease", "CD33_specific", "Control_ASO_effect")) {
    a <- x[x$series == "GSE74437" & x$tissue == tissue & x$contrast == contrast, c("probe_id", "logFC", "adj.P.Val")]
    b <- x[x$series == "GSE74438" & x$tissue == tissue & x$contrast == contrast, c("probe_id", "logFC", "adj.P.Val")]
    m <- merge(a, b, by = "probe_id", suffixes = c("_75ug", "_50ug"))
    comparison_rows[[k]] <- data.frame(
      tissue = tissue, contrast = contrast, common_probes = nrow(m),
      spearman_rho = cor(m$logFC_75ug, m$logFC_50ug, method = "spearman"),
      pearson_r = cor(m$logFC_75ug, m$logFC_50ug),
      direction_concordance = mean(sign(m$logFC_75ug) == sign(m$logFC_50ug)),
      replicated_FDR05 = sum(m$adj.P.Val_75ug < 0.05 & m$adj.P.Val_50ug < 0.05 & sign(m$logFC_75ug) == sign(m$logFC_50ug))
    )
    k <- k + 1
  }
  for (series in c("GSE74437", "GSE74438")) {
    d <- x[x$series == series & x$tissue == tissue & x$contrast == "Disease", c("probe_id", "logFC", "adj.P.Val", "Gene symbol")]
    tr <- x[x$series == series & x$tissue == tissue & x$contrast == "CD33_specific", c("probe_id", "logFC", "adj.P.Val", "Gene symbol")]
    m <- merge(d, tr, by = "probe_id", suffixes = c("_disease", "_treatment"))
    rev <- m[m$adj.P.Val_disease < 0.05 & m$adj.P.Val_treatment < 0.05 & sign(m$logFC_disease) == -sign(m$logFC_treatment), ]
    reversal_rows[[paste(series, tissue)]] <- data.frame(series = series, tissue = tissue, reversal_probes = nrow(rev), reversal_genes = length(unique(rev[["Gene symbol_disease"]][rev[["Gene symbol_disease"]] != ""])))
  }
  a <- x[x$series == "GSE74437" & x$tissue == tissue, ]
  b <- x[x$series == "GSE74438" & x$tissue == tissue, ]
  get_reversal <- function(z) {
    d <- z[z$contrast == "Disease", c("probe_id", "logFC", "adj.P.Val", "Gene symbol")]
    tr <- z[z$contrast == "CD33_specific", c("probe_id", "logFC", "adj.P.Val")]
    m <- merge(d, tr, by = "probe_id", suffixes = c("_disease", "_treatment"))
    m[m$adj.P.Val_disease < 0.05 & m$adj.P.Val_treatment < 0.05 & sign(m$logFC_disease) == -sign(m$logFC_treatment), ]
  }
  ra <- get_reversal(a); rb <- get_reversal(b)
  rr <- merge(ra, rb, by = "probe_id", suffixes = c("_75ug", "_50ug"))
  rr <- rr[sign(rr$logFC_disease_75ug) == sign(rr$logFC_disease_50ug) & sign(rr$logFC_treatment_75ug) == sign(rr$logFC_treatment_50ug), ]
  if (nrow(rr)) {
    rr$tissue <- tissue
    replicated_rows[[tissue]] <- rr
  }
}

comparisons <- do.call(rbind, comparison_rows)
reversal <- do.call(rbind, reversal_rows)
replicated <- if (length(replicated_rows)) do.call(rbind, replicated_rows) else data.frame()
write.csv(comparisons, file.path(res_dir, "cross_series_reproducibility.csv"), row.names = FALSE)
write.csv(reversal, file.path(res_dir, "within_stratum_reversal_counts.csv"), row.names = FALSE)
write.csv(replicated, file.path(res_dir, "replicated_reversal_probes.csv"), row.names = FALSE)

pdf(file.path(fig_dir, "cross_series_effect_correlations.pdf"), width = 10, height = 8)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (tissue in c("Hippocampus", "Frontal_cortex")) {
  for (contrast in c("Disease", "CD33_specific")) {
    a <- x[x$series == "GSE74437" & x$tissue == tissue & x$contrast == contrast, c("probe_id", "logFC")]
    b <- x[x$series == "GSE74438" & x$tissue == tissue & x$contrast == contrast, c("probe_id", "logFC")]
    m <- merge(a, b, by = "probe_id", suffixes = c("_75", "_50"))
    smoothScatter(m$logFC_75, m$logFC_50, xlab = "75 ug/day logFC", ylab = "50 ug/day logFC", main = paste(tissue, contrast))
    abline(h = 0, v = 0, col = "grey70", lty = 2); abline(0, 1, col = "#B24C63")
  }
}
dev.off()

lines <- c(
  "# Global feasibility summary", "",
  "## Differential-signal counts", "", paste(capture.output(print(counts)), collapse = "\n"), "",
  "## Cross-series reproducibility", "", paste(capture.output(print(comparisons)), collapse = "\n"), "",
  "## Disease-to-treatment reversal", "", paste(capture.output(print(reversal)), collapse = "\n"), "",
  paste0("Strictly replicated reversal probes across both dose experiments: ", nrow(replicated), ".")
)
writeLines(lines, file.path(res_dir, "GLOBAL_FEASIBILITY.md"))
message("Global feasibility summary completed")
