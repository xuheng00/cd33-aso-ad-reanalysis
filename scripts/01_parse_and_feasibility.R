#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(limma))

root <- normalizePath(file.path(dirname(commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))]), ".."), mustWork = FALSE)
root <- normalizePath(file.path(getwd(), "AD_CD33_reanalysis"))
raw_dir <- file.path(root, "data", "raw")
processed_dir <- file.path(root, "data", "processed")
results_dir <- file.path(root, "results", "feasibility")
figures_dir <- file.path(root, "figures", "feasibility")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

parse_matrix <- function(path, series, dose) {
  lines <- readLines(gzfile(path), warn = FALSE)
  meta_lines <- lines[grepl("^!Sample_", lines)]
  parse_meta <- function(prefix, occurrence = 1) {
    hits <- meta_lines[grepl(paste0("^", prefix, "\\t"), meta_lines)]
    if (length(hits) < occurrence) return(NULL)
    x <- strsplit(hits[occurrence], "\\t", fixed = FALSE)[[1]][-1]
    gsub('^"|"$', "", x)
  }
  gsm <- parse_meta("!Sample_geo_accession")
  title <- parse_meta("!Sample_title")
  source <- parse_meta("!Sample_source_name_ch1")
  desc <- parse_meta("!Sample_description")
  char_hits <- meta_lines[grepl("^!Sample_characteristics_ch1\\t", meta_lines)]
  chars <- lapply(char_hits, function(z) gsub('^"|"$', "", strsplit(z, "\\t")[[1]][-1]))
  names(chars) <- vapply(chars, function(v) sub(":.*$", "", v[1]), character(1))
  get_char <- function(pattern) {
    idx <- grep(pattern, names(chars), ignore.case = TRUE)
    if (!length(idx)) return(rep(NA_character_, length(gsm)))
    sub("^[^:]+:\\s*", "", chars[[idx[1]]])
  }
  genotype <- get_char("genotype")
  treatment <- get_char("treatment")
  tissue <- get_char("tissue")
  animal <- get_char("animal")
  sex <- get_char("sex")
  surgery_date <- get_char("surgery date")
  collection_date <- get_char("collection date")
  treatment_std <- ifelse(grepl("CD33", treatment, ignore.case = TRUE), "CD33_ASO",
                    ifelse(grepl("CON|Ctrl", treatment, ignore.case = TRUE), "Control_ASO", "PBS"))
  tissue_std <- ifelse(grepl("hipp", tissue, ignore.case = TRUE), "Hippocampus", "Frontal_cortex")
  genotype_std <- ifelse(grepl("APP", genotype, ignore.case = TRUE), "APPPS1", "WT")
  group <- paste(genotype_std, treatment_std, sep = "_")
  metadata <- data.frame(
    gsm = gsm, series = series, dose_ug_day = dose, title = title,
    source = source, tissue = tissue_std, genotype = genotype_std,
    treatment = treatment_std, group = group, animal_id = animal, sex = sex,
    surgery_date = surgery_date, collection_date = collection_date,
    array_id = desc, stringsAsFactors = FALSE
  )
  begin <- grep("^!series_matrix_table_begin", lines) + 1
  end <- grep("^!series_matrix_table_end", lines) - 1
  tab <- read.delim(text = paste(lines[begin:end], collapse = "\n"), check.names = FALSE, quote = '"')
  rownames(tab) <- tab[[1]]
  tab <- as.matrix(tab[-1])
  storage.mode(tab) <- "numeric"
  colnames(tab) <- gsub('^"|"$', "", colnames(tab))
  stopifnot(identical(colnames(tab), metadata$gsm))
  list(meta = metadata, expr = tab)
}

d37 <- parse_matrix(file.path(raw_dir, "GSE74437_series_matrix.txt.gz"), "GSE74437", 75)
d38 <- parse_matrix(file.path(raw_dir, "GSE74438_series_matrix.txt.gz"), "GSE74438", 50)
meta <- rbind(d37$meta, d38$meta)
write.csv(meta, file.path(processed_dir, "sample_manifest_170.csv"), row.names = FALSE)

counts <- as.data.frame(with(meta, table(series, tissue, group)))
counts <- counts[counts$Freq > 0, ]
write.csv(counts, file.path(results_dir, "group_counts.csv"), row.names = FALSE)

annot <- read.delim(gzfile(file.path(raw_dir, "GPL6885.annot.gz")), skip = 28, check.names = FALSE, quote = "", comment.char = "")
annot <- annot[, c("ID", "Gene symbol", "Gene ID")]
write.csv(annot, file.path(processed_dir, "GPL6885_probe_annotation.csv"), row.names = FALSE)

analyze_stratum <- function(dat, tissue_name) {
  idx <- dat$meta$tissue == tissue_name
  m <- dat$meta[idx, , drop = FALSE]
  e <- dat$expr[, idx, drop = FALSE]
  keep_groups <- c("WT_PBS", "APPPS1_PBS", "APPPS1_Control_ASO", "APPPS1_CD33_ASO")
  stopifnot(all(keep_groups %in% m$group))
  m$group <- factor(m$group, levels = keep_groups)
  design <- model.matrix(~0 + group, m)
  colnames(design) <- levels(m$group)
  fit <- lmFit(e, design)
  cm <- makeContrasts(
    Disease = APPPS1_PBS - WT_PBS,
    CD33_specific = APPPS1_CD33_ASO - APPPS1_Control_ASO,
    Control_ASO_effect = APPPS1_Control_ASO - APPPS1_PBS,
    levels = design
  )
  fit2 <- eBayes(contrasts.fit(fit, cm), trend = TRUE, robust = TRUE)
  out <- lapply(colnames(cm), function(k) {
    tt <- topTable(fit2, coef = k, number = Inf, sort.by = "none")
    tt$probe_id <- rownames(tt)
    tt$contrast <- k
    tt$series <- unique(m$series)
    tt$tissue <- tissue_name
    merge(tt, annot, by.x = "probe_id", by.y = "ID", all.x = TRUE, sort = FALSE)
  })
  do.call(rbind, out)
}

all_res <- rbind(
  analyze_stratum(d37, "Hippocampus"),
  analyze_stratum(d37, "Frontal_cortex"),
  analyze_stratum(d38, "Hippocampus"),
  analyze_stratum(d38, "Frontal_cortex")
)
write.csv(all_res, file.path(results_dir, "all_probe_contrasts.csv"), row.names = FALSE)

cd33 <- all_res[tolower(all_res[["Gene symbol"]]) == "cd33", ]
write.csv(cd33, file.path(results_dir, "cd33_probe_results.csv"), row.names = FALSE)

pdf(file.path(figures_dir, "CD33_effects.pdf"), width = 9, height = 5)
op <- par(mar = c(8, 5, 2, 1))
z <- cd33[cd33$contrast %in% c("Disease", "CD33_specific"), ]
labs <- paste(z$series, z$tissue, z$contrast, sep = "\n")
cols <- ifelse(z$contrast == "Disease", "#B24C63", "#3B7EA1")
bp <- barplot(z$logFC, names.arg = labs, las = 2, col = cols,
              ylab = "CD33 log2 fold change", main = "CD33 disease and treatment effects")
abline(h = 0, lty = 2)
par(op)
dev.off()

summary_lines <- c(
  "# Feasibility results",
  "",
  paste0("- Total samples parsed: ", nrow(meta), "."),
  paste0("- GSE74437 samples: ", sum(meta$series == "GSE74437"), "."),
  paste0("- GSE74438 samples: ", sum(meta$series == "GSE74438"), "."),
  "- All four planned groups are present in both brain regions and both series.",
  "- Primary feasibility gate: CD33 must decrease in CD33-ASO versus control-ASO in both series, preferably in both regions.",
  "",
  "## CD33 probe-level results",
  "",
  paste(capture.output(print(cd33[, c("series", "tissue", "contrast", "logFC", "P.Value", "adj.P.Val", "probe_id")])), collapse = "\n")
)
writeLines(summary_lines, file.path(results_dir, "FEASIBILITY_RESULTS.md"))
saveRDS(list(metadata = meta, GSE74437 = d37, GSE74438 = d38, annotation = annot), file.path(processed_dir, "parsed_series.rds"), compress = "xz")
message("Completed feasibility analysis. Outputs: ", results_dir)
