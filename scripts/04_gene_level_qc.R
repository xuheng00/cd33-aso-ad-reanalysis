#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(limma))
root <- normalizePath(file.path(getwd(), "AD_CD33_reanalysis"))
obj <- readRDS(file.path(root, "data", "processed", "parsed_series.rds"))
out_data <- file.path(root, "data", "processed")
out_res <- file.path(root, "results", "gene_level")
out_fig <- file.path(root, "figures", "qc")
dir.create(out_res, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)

common <- intersect(rownames(obj$GSE74437$expr), rownames(obj$GSE74438$expr))
ann <- obj$annotation[match(common, obj$annotation$ID), ]
valid <- !is.na(ann$`Gene symbol`) & ann$`Gene symbol` != "" & !grepl("///", ann$`Gene symbol`, fixed = TRUE)
common <- common[valid]; ann <- ann[valid, ]
pooled_mean <- (rowMeans(obj$GSE74437$expr[common, , drop = FALSE]) + rowMeans(obj$GSE74438$expr[common, , drop = FALSE])) / 2
probe_map <- data.frame(probe_id = common, gene_symbol = ann$`Gene symbol`, pooled_mean = pooled_mean, stringsAsFactors = FALSE)
probe_map <- probe_map[order(probe_map$gene_symbol, -probe_map$pooled_mean), ]
probe_map <- probe_map[!duplicated(probe_map$gene_symbol), ]
write.csv(probe_map, file.path(out_data, "representative_probe_per_gene.csv"), row.names = FALSE)

collapse <- function(e) {
  z <- e[probe_map$probe_id, , drop = FALSE]
  rownames(z) <- probe_map$gene_symbol
  z
}
g37 <- collapse(obj$GSE74437$expr)
g38 <- collapse(obj$GSE74438$expr)

plot_pca <- function(e, m, series_name, tissue_name) {
  idx <- m$tissue == tissue_name
  z <- e[, idx, drop = FALSE]
  vars <- apply(z, 1, var)
  z <- z[order(vars, decreasing = TRUE)[seq_len(min(1000, nrow(z)))], , drop = FALSE]
  p <- prcomp(t(z), center = TRUE, scale. = FALSE)
  ve <- round(100 * summary(p)$importance[2, 1:2], 1)
  mm <- m[idx, ]
  cols <- c(WT_PBS = "#777777", APPPS1_PBS = "#B24C63", APPPS1_Control_ASO = "#D49A3A", APPPS1_CD33_ASO = "#3B7EA1")
  plot(p$x[,1], p$x[,2], pch = 19, col = cols[mm$group],
       xlab = paste0("PC1 (", ve[1], "%)"), ylab = paste0("PC2 (", ve[2], "%)"),
       main = paste(series_name, tissue_name))
  legend("topright", legend = names(cols), col = cols, pch = 19, cex = .7, bty = "n")
  data.frame(gsm = mm$gsm, series = series_name, tissue = tissue_name, group = mm$group,
             PC1 = p$x[,1], PC2 = p$x[,2])
}

png(file.path(out_fig, "PCA_by_series_region.png"), width = 1800, height = 1400, res = 180)
par(mfrow = c(2,2), mar = c(4,4,3,1))
pca <- rbind(
  plot_pca(g37, obj$GSE74437$meta, "GSE74437", "Hippocampus"),
  plot_pca(g37, obj$GSE74437$meta, "GSE74437", "Frontal_cortex"),
  plot_pca(g38, obj$GSE74438$meta, "GSE74438", "Hippocampus"),
  plot_pca(g38, obj$GSE74438$meta, "GSE74438", "Frontal_cortex")
)
dev.off()
write.csv(pca, file.path(out_res, "pca_coordinates.csv"), row.names = FALSE)

saveRDS(list(probe_map = probe_map, GSE74437 = list(expr = g37, meta = obj$GSE74437$meta), GSE74438 = list(expr = g38, meta = obj$GSE74438$meta)), file.path(out_data, "gene_level_series.rds"), compress = "xz")
writeLines(c(
  "# Gene-level preprocessing and QC", "",
  paste0("- Common platform probes across series: ", length(intersect(rownames(obj$GSE74437$expr), rownames(obj$GSE74438$expr))), "."),
  paste0("- Unique annotated genes retained: ", nrow(probe_map), "."),
  "- For genes with multiple probes, the common probe with the highest pooled mean expression was selected.",
  "- PCA was performed separately by series and brain region using the 1,000 most variable retained genes."
), file.path(out_res, "GENE_LEVEL_QC.md"))
message("Gene-level matrices and PCA completed")
