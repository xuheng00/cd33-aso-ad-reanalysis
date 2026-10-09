#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(limma))
root <- normalizePath(file.path(getwd(), "AD_CD33_reanalysis"))
obj <- readRDS(file.path(root, "data", "processed", "gene_level_series.rds"))
out_dir <- file.path(root, "results", "meta_analysis")
fig_dir <- file.path(root, "figures", "meta_analysis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

fit_stratum <- function(series, tissue) {
  m <- obj[[series]]$meta[obj[[series]]$meta$tissue == tissue, ]
  e <- obj[[series]]$expr[, obj[[series]]$meta$tissue == tissue, drop = FALSE]
  m$group <- factor(m$group, levels = c("WT_PBS", "APPPS1_PBS", "APPPS1_Control_ASO", "APPPS1_CD33_ASO"))
  m$surgery_date <- factor(m$surgery_date)
  design <- model.matrix(~0 + group + surgery_date, m)
  colnames(design) <- make.names(colnames(design))
  fit <- lmFit(e, design)
  cm <- makeContrasts(
    Disease = groupAPPPS1_PBS - groupWT_PBS,
    CD33_specific = groupAPPPS1_CD33_ASO - groupAPPPS1_Control_ASO,
    Control_ASO_effect = groupAPPPS1_Control_ASO - groupAPPPS1_PBS,
    levels = design
  )
  fit2 <- eBayes(contrasts.fit(fit, cm), trend = TRUE, robust = TRUE)
  do.call(rbind, lapply(colnames(cm), function(k) {
    z <- topTable(fit2, coef = k, number = Inf, sort.by = "none")
    z$gene_symbol <- rownames(z)
    z$SE <- abs(z$logFC / z$t)
    z$series <- series; z$tissue <- tissue; z$contrast <- k
    z
  }))
}

res <- do.call(rbind, lapply(c("GSE74437", "GSE74438"), function(s)
  do.call(rbind, lapply(c("Hippocampus", "Frontal_cortex"), function(t) fit_stratum(s, t)))))
write.csv(res, file.path(out_dir, "gene_level_stratum_contrasts.csv"), row.names = FALSE)

meta_pair <- function(a, b) {
  m <- merge(a[, c("gene_symbol", "logFC", "SE", "adj.P.Val")], b[, c("gene_symbol", "logFC", "SE", "adj.P.Val")], by = "gene_symbol", suffixes = c("_75ug", "_50ug"))
  w1 <- 1 / m$SE_75ug^2; w2 <- 1 / m$SE_50ug^2
  m$meta_logFC <- (w1*m$logFC_75ug + w2*m$logFC_50ug) / (w1+w2)
  m$meta_SE <- sqrt(1/(w1+w2))
  m$meta_z <- m$meta_logFC / m$meta_SE
  m$meta_P <- 2*pnorm(-abs(m$meta_z))
  m$Q <- w1*(m$logFC_75ug-m$meta_logFC)^2 + w2*(m$logFC_50ug-m$meta_logFC)^2
  m$Q_P <- pchisq(m$Q, df=1, lower.tail=FALSE)
  m$I2 <- pmax(0, (m$Q-1)/m$Q)*100
  m$direction_concordant <- sign(m$logFC_75ug)==sign(m$logFC_50ug)
  m$meta_FDR <- p.adjust(m$meta_P, "BH")
  m
}

meta_list <- list()
for (tissue in c("Hippocampus", "Frontal_cortex")) for (contrast in c("Disease", "CD33_specific", "Control_ASO_effect")) {
  a <- res[res$series=="GSE74437" & res$tissue==tissue & res$contrast==contrast, ]
  b <- res[res$series=="GSE74438" & res$tissue==tissue & res$contrast==contrast, ]
  m <- meta_pair(a,b); m$tissue <- tissue; m$contrast <- contrast
  meta_list[[paste(tissue,contrast)]] <- m
}
meta <- do.call(rbind, meta_list)
write.csv(meta, file.path(out_dir, "fixed_effect_meta_by_region.csv"), row.names = FALSE)

summary <- aggregate(cbind(meta_FDR05 = meta$meta_FDR < .05,
                           concordant_FDR05 = meta$meta_FDR < .05 & meta$direction_concordant,
                           low_heterogeneity_FDR05 = meta$meta_FDR < .05 & meta$direction_concordant & meta$I2 < 50) ~ tissue + contrast, meta, sum)
write.csv(summary, file.path(out_dir, "meta_signal_counts.csv"), row.names = FALSE)

top <- meta[meta$meta_FDR < .05 & meta$direction_concordant & meta$I2 < 50, ]
top <- top[order(top$tissue, top$contrast, top$meta_FDR), ]
write.csv(top, file.path(out_dir, "robust_meta_signals.csv"), row.names = FALSE)

cd33 <- meta[tolower(meta$gene_symbol)=="cd33" & meta$contrast %in% c("Disease","CD33_specific"), ]
png(file.path(fig_dir,"CD33_meta_forest.png"),width=1700,height=1200,res=180)
par(mfrow=c(2,2),mar=c(4,8,3,1))
for(i in seq_len(nrow(cd33))){
  z<-cd33[i,]; est<-c(z$logFC_75ug,z$logFC_50ug,z$meta_logFC); se<-c(z$SE_75ug,z$SE_50ug,z$meta_SE)
  plot(est,3:1,xlim=range(c(est-1.96*se,est+1.96*se)),ylim=c(.5,3.5),yaxt='n',pch=c(19,19,18),xlab='log2 fold change',ylab='',main=paste(z$tissue,z$contrast)); axis(2,at=3:1,labels=c('75 ug/day','50 ug/day','Fixed-effect meta'),las=1); segments(est-1.96*se,3:1,est+1.96*se,3:1); abline(v=0,lty=2,col='grey60')
}
dev.off()

writeLines(c("# Gene-level cross-cohort meta-analysis","",
             "Fixed-effect inverse-variance meta-analysis was performed separately for hippocampus and frontal cortex. Surgery date was included in each cohort-specific model.","",
             "## Signal counts","",paste(capture.output(print(summary)),collapse='\n'),"",
             "## CD33 meta-analysis","",paste(capture.output(print(cd33[,c('tissue','contrast','logFC_75ug','logFC_50ug','meta_logFC','meta_P','meta_FDR','I2')])),collapse='\n')),
           file.path(out_dir,"META_ANALYSIS_SUMMARY.md"))
message("Gene-level meta-analysis completed")
