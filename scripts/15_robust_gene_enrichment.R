#!/usr/bin/env Rscript
# 15_robust_gene_enrichment.R
# Functional characterization of the robust CD33-specific and control-ASO gene sets:
#   (1) GO-BP / KEGG over-representation (hypergeometric, BH-adjusted) for each region;
#   (2) overlap between the CD33-specific and control-ASO robust sets.
suppressPackageStartupMessages(library(limma))

root <- normalizePath(getwd())
robust <- read.csv(file.path(root, "results", "meta_analysis", "robust_meta_signals.csv"), check.names = FALSE)
go_map <- read.csv(file.path(root, "results", "pathway_reversal", "GO_BP_gene_sets.csv"), check.names = FALSE)
kegg_map <- read.csv(file.path(root, "results", "pathway_reversal", "KEGG_gene_sets.csv"), check.names = FALSE)
out_dir <- file.path(root, "results", "robust_gene_enrichment")
fig_dir <- file.path(root, "figures", "robust_gene_enrichment")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

background <- unique(robust$gene_symbol) # not used as gene-space; true background from probe map
probe_map <- read.csv(file.path(root, "data", "processed", "representative_probe_per_gene.csv"), check.names = FALSE)
gene_space <- unique(probe_map$gene_symbol)

make_sets <- function(map) {
  sets <- split(map$gene, map$pathway)
  sets <- lapply(sets, unique)
  sets <- sets[lengths(sets) >= 10 & lengths(sets) <= 500]
  names_lookup <- setNames(map$name, map$pathway)
  list(sets = sets, names = names_lookup)
}
collections <- list(GO_BP = make_sets(go_map), KEGG = make_sets(kegg_map))

ora <- function(query, sets, names_lookup, gene_space, min_genes = 5) {
  query <- unique(query)
  n_back <- length(gene_space)
  n_q <- length(query)
  rows <- list(); k <- 1
  for (p in names(sets)) {
    s <- intersect(sets[[p]], gene_space)
    if (length(s) < min_genes) next
    overlap <- length(intersect(query, s))
    if (overlap < min_genes) next
    # Hypergeometric: Pr(X >= overlap)
    pval <- phyper(overlap - 1, length(s), n_back - length(s), n_q, lower.tail = FALSE)
    rows[[k]] <- data.frame(pathway = p, name = unname(names_lookup[p]),
                            set_size = length(s), query_overlap = overlap,
                            expected = round(n_q * length(s) / n_back, 2),
                            fold_enrichment = round((overlap / n_q) / (length(s) / n_back), 2),
                            P = pval, stringsAsFactors = FALSE)
    k <- k + 1
  }
  z <- do.call(rbind, rows)
  if (!is.null(z) && nrow(z)) z$FDR <- p.adjust(z$P, "BH")
  z
}

results <- list()
for (tissue in c("Frontal_cortex", "Hippocampus")) {
  for (contrast in c("CD33_specific", "Control_ASO_effect")) {
    q <- robust$gene_symbol[robust$tissue == tissue & robust$contrast == contrast]
    for (cn in names(collections)) {
      z <- ora(q, collections[[cn]]$sets, collections[[cn]]$names, gene_space)
      if (!is.null(z) && nrow(z)) {
        z$tissue <- tissue; z$contrast <- contrast; z$collection <- cn
        z <- z[order(z$FDR), ]
      }
      results[[paste(tissue, contrast, cn)]] <- z
    }
  }
}
enrich <- do.call(rbind, results)
write.csv(enrich, file.path(out_dir, "robust_gene_set_enrichment.csv"), row.names = FALSE)

# Overlap between CD33-specific and control-ASO robust sets per region
overlap_rows <- list(); k <- 1
for (tissue in c("Frontal_cortex", "Hippocampus")) {
  cd33 <- unique(robust$gene_symbol[robust$tissue == tissue & robust$contrast == "CD33_specific"])
  ctrl <- unique(robust$gene_symbol[robust$tissue == tissue & robust$contrast == "Control_ASO_effect"])
  ov <- intersect(cd33, ctrl)
  overlap_rows[[k]] <- data.frame(tissue = tissue,
                                  cd33_n = length(cd33), control_n = length(ctrl),
                                  overlap = length(ov),
                                  overlap_frac_cd33 = round(length(ov) / length(cd33), 3),
                                  overlap_frac_control = round(length(ov) / length(ctrl), 3),
                                  jaccard = round(length(ov) / length(union(cd33, ctrl)), 3),
                                  stringsAsFactors = FALSE)
  k <- k + 1
  write.csv(data.frame(gene = sort(ov)), file.path(out_dir, paste0("overlap_", tissue, ".csv")), row.names = FALSE)
}
overlap <- do.call(rbind, overlap_rows)
write.csv(overlap, file.path(out_dir, "cd33_vs_control_overlap.csv"), row.names = FALSE)

# Top enriched GO-BP terms for the CD33-specific sets (for the figure)
top <- enrich[enrich$collection == "GO_BP" & enrich$contrast == "CD33_specific" & enrich$FDR < 0.05, ]
if (nrow(top)) {
  top <- top[order(top$tissue, top$FDR), ]
  plot_data <- top[seq_len(min(nrow(top), 20)), ]
  png(file.path(fig_dir, "top_CD33_robust_GO_BP.png"), width = 2000, height = 1500, res = 200)
  par(mar = c(5, 20, 3, 2))
  pd <- plot_data[order(plot_data$tissue, plot_data$fold_enrichment), ]
  cols <- ifelse(pd$tissue == "Frontal_cortex", "#B24C63", "#3B7EA1")
  barplot(rev(pd$fold_enrichment), names.arg = rev(substr(pd$name, 1, 55)),
          horiz = TRUE, las = 1, col = rev(cols), border = NA,
          xlab = "Fold enrichment (CD33-specific robust genes)",
          main = "Top GO-BP terms: CD33-specific robust gene sets (FDR < 0.05)")
  legend("bottomright", c("Frontal cortex", "Hippocampus"), fill = c("#B24C63", "#3B7EA1"), bty = "n", cex = 0.8)
  dev.off()
}

sink(file.path(out_dir, "ROBUST_ENRICHMENT_SUMMARY.md"))
cat("# Functional characterization of robust gene sets\n\n")
cat("Gene sets: robust response = meta-FDR < 0.05, concordant direction, I2 < 50%.\n")
cat("Over-representation: hypergeometric test against the 12,272-gene measured space; BH FDR.\n\n")
cat("## CD33-specific vs control-ASO overlap\n\n")
print(overlap, row.names = FALSE)
cat("\n## GO-BP enrichment (CD33-specific, FDR < 0.05)\n\n")
if (nrow(top)) print(head(top[, c("tissue", "name", "set_size", "query_overlap", "expected", "fold_enrichment", "P", "FDR")], 20), row.names = FALSE) else cat("None\n")
cat("\n## GO-BP enrichment (control-ASO background, FDR < 0.05, top 20)\n\n")
ctrl_top <- enrich[enrich$collection == "GO_BP" & enrich$contrast == "Control_ASO_effect" & enrich$FDR < 0.05, ]
if (nrow(ctrl_top)) print(head(ctrl_top[order(ctrl_top$tissue, ctrl_top$FDR), c("tissue", "name", "set_size", "query_overlap", "fold_enrichment", "FDR")], 20), row.names = FALSE) else cat("None\n")
cat("\n## KEGG enrichment (CD33-specific, FDR < 0.05)\n\n")
kegg_top <- enrich[enrich$collection == "KEGG" & enrich$contrast == "CD33_specific" & enrich$FDR < 0.05, ]
if (nrow(kegg_top)) print(kegg_top[order(kegg_top$tissue, kegg_top$FDR), c("tissue", "name", "set_size", "query_overlap", "fold_enrichment", "FDR")], row.names = FALSE) else cat("None\n")
sink()

message("Robust gene enrichment completed; results in ", out_dir)
