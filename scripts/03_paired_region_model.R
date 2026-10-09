#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(limma))
root <- normalizePath(file.path(getwd(), "AD_CD33_reanalysis"))
obj <- readRDS(file.path(root, "data", "processed", "parsed_series.rds"))
annot <- obj$annotation
out_dir <- file.path(root, "results", "feasibility")

fit_series <- function(dat, series_name) {
  m <- dat$meta
  e <- dat$expr
  m$cell <- factor(paste(m$tissue, m$group, sep = "."))
  m$surgery_date <- factor(m$surgery_date)
  design <- model.matrix(~0 + cell + surgery_date, m)
  colnames(design) <- make.names(colnames(design))
  corfit <- duplicateCorrelation(e, design, block = m$animal_id)
  fit <- lmFit(e, design, block = m$animal_id, correlation = corfit$consensus)
  cm <- makeContrasts(
    HC_Disease = cellHippocampus.APPPS1_PBS - cellHippocampus.WT_PBS,
    HC_CD33 = cellHippocampus.APPPS1_CD33_ASO - cellHippocampus.APPPS1_Control_ASO,
    FC_Disease = cellFrontal_cortex.APPPS1_PBS - cellFrontal_cortex.WT_PBS,
    FC_CD33 = cellFrontal_cortex.APPPS1_CD33_ASO - cellFrontal_cortex.APPPS1_Control_ASO,
    levels = design
  )
  fit2 <- eBayes(contrasts.fit(fit, cm), trend = TRUE, robust = TRUE)
  res <- do.call(rbind, lapply(colnames(cm), function(k) {
    z <- topTable(fit2, coef = k, number = Inf, sort.by = "none")
    z$probe_id <- rownames(z); z$contrast <- k; z$series <- series_name
    merge(z, annot, by.x = "probe_id", by.y = "ID", all.x = TRUE, sort = FALSE)
  }))
  list(results = res, within_animal_correlation = corfit$consensus)
}

a <- fit_series(obj$GSE74437, "GSE74437")
b <- fit_series(obj$GSE74438, "GSE74438")
res <- rbind(a$results, b$results)
write.csv(res, file.path(out_dir, "paired_region_probe_contrasts.csv"), row.names = FALSE)

cd33 <- res[tolower(res[["Gene symbol"]]) == "cd33", ]
counts <- aggregate(I(adj.P.Val < 0.05) ~ series + contrast, res, sum)
names(counts)[3] <- "FDR05_probes"
write.csv(cd33, file.path(out_dir, "paired_region_cd33_results.csv"), row.names = FALSE)
write.csv(counts, file.path(out_dir, "paired_region_signal_counts.csv"), row.names = FALSE)

lines <- c(
  "# Paired-region sensitivity model", "",
  paste0("- GSE74437 estimated within-animal cross-region correlation: ", signif(a$within_animal_correlation, 3)),
  paste0("- GSE74438 estimated within-animal cross-region correlation: ", signif(b$within_animal_correlation, 3)),
  "- Models include surgery date and use animal ID as a repeated-measures block.", "",
  "## CD33 results", "", paste(capture.output(print(cd33[, c("series", "contrast", "logFC", "P.Value", "adj.P.Val", "probe_id")])), collapse = "\n"), "",
  "## FDR-significant probe counts", "", paste(capture.output(print(counts)), collapse = "\n")
)
writeLines(lines, file.path(out_dir, "PAIRED_REGION_SENSITIVITY.md"))
message("Paired-region sensitivity model completed")
