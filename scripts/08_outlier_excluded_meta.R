#!/usr/bin/env Rscript
root <- normalizePath(file.path(getwd(), "AD_CD33_reanalysis"))
x <- read.csv(file.path(root,"results","gene_level","outlier_sensitivity_gene_contrasts.csv"),check.names=FALSE)
x$SE <- with(x, ifelse(is.finite(t) & abs(t)>1e-12, abs(logFC/t), NA_real_))

meta_two <- function(z){
  a <- z[z$series=="GSE74437",]; b <- z[z$series=="GSE74438",]
  m <- merge(a,b,by="gene_symbol",suffixes=c("_75ug","_50ug"))
  w1 <- 1/m$SE_75ug^2; w2 <- 1/m$SE_50ug^2
  m$meta_logFC <- (w1*m$logFC_75ug+w2*m$logFC_50ug)/(w1+w2)
  m$meta_SE <- sqrt(1/(w1+w2)); m$meta_z <- m$meta_logFC/m$meta_SE
  m$meta_P <- 2*pnorm(-abs(m$meta_z)); m$meta_FDR <- p.adjust(m$meta_P,"BH")
  m$Q <- w1*(m$logFC_75ug-m$meta_logFC)^2+w2*(m$logFC_50ug-m$meta_logFC)^2
  m$Q_P <- pchisq(m$Q,df=1,lower.tail=FALSE)
  m$I2 <- pmax(0,(m$Q-1)/m$Q)*100
  m$direction_concordant <- sign(m$logFC_75ug)==sign(m$logFC_50ug)
  m
}

out <- list()
for(tissue in unique(x$tissue)) for(contrast in unique(x$contrast)){
  z <- x[x$tissue==tissue & x$contrast==contrast,]
  m <- meta_two(z); m$tissue <- tissue; m$contrast <- contrast
  out[[paste(tissue,contrast)]] <- m[,c("gene_symbol","logFC_75ug","SE_75ug","adj.P.Val_75ug","logFC_50ug","SE_50ug","adj.P.Val_50ug","meta_logFC","meta_SE","meta_z","meta_P","meta_FDR","Q","Q_P","I2","direction_concordant","tissue","contrast")]
}
res <- do.call(rbind,out)
dest <- file.path(root,"results","meta_analysis","fixed_effect_meta_outlier_excluded.csv")
write.csv(res,dest,row.names=FALSE)

counts <- aggregate(cbind(meta_FDR05=res$meta_FDR<0.05,
                          concordant_FDR05=res$meta_FDR<0.05 & res$direction_concordant,
                          low_heterogeneity_FDR05=res$meta_FDR<0.05 & res$direction_concordant & res$I2<50)~tissue+contrast,res,sum)
write.csv(counts,file.path(root,"results","meta_analysis","meta_signal_counts_outlier_excluded.csv"),row.names=FALSE)
message("Outlier-excluded fixed-effect meta-analysis completed")
