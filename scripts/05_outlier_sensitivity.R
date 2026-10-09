#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(limma))
root <- normalizePath(file.path(getwd(), "AD_CD33_reanalysis"))
obj <- readRDS(file.path(root, "data", "processed", "gene_level_series.rds"))
out_dir <- file.path(root, "results", "gene_level")

all_outliers <- list(); all_results <- list(); all_counts <- list(); k <- 1
for (series in c("GSE74437", "GSE74438")) {
  for (tissue in c("Hippocampus", "Frontal_cortex")) {
    m <- obj[[series]]$meta[obj[[series]]$meta$tissue == tissue, ]
    e <- obj[[series]]$expr[, obj[[series]]$meta$tissue == tissue, drop = FALSE]
    cmat <- cor(e)
    medcor <- apply(cmat, 2, function(v) median(v[v < 0.999999]))
    rz <- (medcor - median(medcor)) / (1.4826 * mad(medcor, constant = 1))
    exclude <- rz < -3.5
    all_outliers[[k]] <- data.frame(series = series, tissue = tissue, gsm = m$gsm, group = m$group, median_correlation = medcor, robust_z = rz, excluded = exclude)
    k <- k + 1
    m2 <- m[!exclude, ]; e2 <- e[, !exclude, drop = FALSE]
    m2$group <- factor(m2$group, levels = c("WT_PBS", "APPPS1_PBS", "APPPS1_Control_ASO", "APPPS1_CD33_ASO"))
    m2$surgery_date <- factor(m2$surgery_date)
    design <- model.matrix(~0 + group + surgery_date, m2)
    colnames(design) <- make.names(colnames(design))
    fit <- lmFit(e2, design)
    cm <- makeContrasts(Disease = groupAPPPS1_PBS - groupWT_PBS,
                        CD33_specific = groupAPPPS1_CD33_ASO - groupAPPPS1_Control_ASO,
                        levels = design)
    fit2 <- eBayes(contrasts.fit(fit, cm), trend = TRUE, robust = TRUE)
    for (contrast in colnames(cm)) {
      z <- topTable(fit2, coef = contrast, number = Inf, sort.by = "none")
      z$gene_symbol <- rownames(z); z$series <- series; z$tissue <- tissue; z$contrast <- contrast
      all_results[[length(all_results)+1]] <- z
      all_counts[[length(all_counts)+1]] <- data.frame(series=series,tissue=tissue,contrast=contrast,FDR05=sum(z$adj.P.Val<.05),abs_logFC025_FDR05=sum(z$adj.P.Val<.05 & abs(z$logFC)>=.25))
    }
  }
}
outliers <- do.call(rbind, all_outliers)
res <- do.call(rbind, all_results)
counts <- do.call(rbind, all_counts)
write.csv(outliers, file.path(out_dir, "sample_correlation_outlier_audit.csv"), row.names = FALSE)
write.csv(res, file.path(out_dir, "outlier_sensitivity_gene_contrasts.csv"), row.names = FALSE)
write.csv(counts, file.path(out_dir, "outlier_sensitivity_counts.csv"), row.names = FALSE)
cd33 <- res[tolower(res$gene_symbol)=="cd33",c("series","tissue","contrast","logFC","P.Value","adj.P.Val")]
writeLines(c("# Outlier sensitivity analysis","",
             "Samples with within-stratum median-correlation robust z < -3.5 were excluded only for this sensitivity analysis.","",
             "## Flagged samples","",paste(capture.output(print(outliers[outliers$excluded,])),collapse="\n"),"",
             "## CD33 results after exclusion","",paste(capture.output(print(cd33)),collapse="\n"),"",
             "## Signal counts","",paste(capture.output(print(counts)),collapse="\n")),
           file.path(out_dir,"OUTLIER_SENSITIVITY.md"))
message("Outlier sensitivity analysis completed")
