#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(limma))
root<-normalizePath(file.path(getwd(),"AD_CD33_reanalysis")); rawdir<-file.path(root,"data","human_raw"); procdir<-file.path(root,"data","human_processed"); outdir<-file.path(root,"results","human_validation"); figdir<-file.path(root,"figures","human_validation")

fit_counts <- function(counts,design,coef,minimum_replicates=2){
  lib<-colSums(counts); cpm<-t(t(counts)/lib*1e6); keep<-rowSums(cpm>1)>=minimum_replicates
  counts<-counts[keep,,drop=FALSE]; v<-voom(counts,design,plot=FALSE); f<-eBayes(lmFit(v,design),trend=FALSE,robust=TRUE)
  tt<-topTable(f,coef=coef,number=Inf,sort.by="none"); tt$human_symbol<-rownames(tt); tt$SE<-abs(tt$logFC/tt$t); tt
}

# GSE138852: pool-level pseudobulk (each library pool contains two donors).
cov<-read.csv(gzfile(file.path(rawdir,"GSE138852_covariates.csv.gz")),check.names=FALSE); names(cov)[1]<-"cell"
mcov<-cov[cov$oupSample.cellType=="mg",]; mcov$pool<-sub("^[^_]+_","",mcov$cell); mcov$status<-ifelse(mcov$oupSample.batchCond=="AD","AD","Control")
count_path<-file.path(rawdir,"GSE138852_counts.csv.gz"); hdr<-names(read.csv(gzfile(count_path),nrows=0,check.names=FALSE)); wanted<-match(mcov$cell,hdr)
cc<-rep("NULL",length(hdr)); cc[1]<-"character"; cc[wanted[!is.na(wanted)]]<-"numeric"
raw<-read.csv(gzfile(count_path),check.names=FALSE,colClasses=cc); genes<-raw[[1]]; mat<-as.matrix(raw[,-1]); rownames(mat)<-genes
mcov<-mcov[match(colnames(mat),mcov$cell),]; pools<-unique(mcov$pool); pb<-sapply(pools,function(p)rowSums(mat[,mcov$pool==p,drop=FALSE])); rownames(pb)<-genes
pmeta<-unique(mcov[,c("pool","status")]); pmeta<-pmeta[match(colnames(pb),pmeta$pool),]; pmeta$status<-relevel(factor(pmeta$status),"Control")
design<-model.matrix(~status,pmeta); gr<-fit_counts(pb,design,"statusAD",minimum_replicates=2); gr$series<-"GSE138852"; gr$tissue<-"Entorhinal_cortex"; gr$comparison<-"AD_vs_control_pool"
write.csv(pmeta,file.path(procdir,"GSE138852_microglia_pool_metadata.csv"),row.names=FALSE)

# GSE243292: donor-level microglia pseudobulk; ordered amyloid/tau pathology trend.
cnt<-read.csv(gzfile(file.path(procdir,"GSE243292_microglia_pseudobulk_counts.csv.gz")),check.names=FALSE); genes<-cnt[[1]]; q<-as.matrix(cnt[,-1]); storage.mode(q)<-"numeric"
valid<-!is.na(genes)&genes!=""; q<-rowsum(q[valid,,drop=FALSE],group=genes[valid],reorder=FALSE)
qm<-read.csv(file.path(procdir,"GSE243292_microglia_sample_metadata.csv"),check.names=FALSE); qm<-qm[match(colnames(q),qm$sampleID),]
qm$pathology_score<-c("A-T-"=0,"A+T-"=1,"A+T+"=2)[qm$atscore]
design2<-model.matrix(~pathology_score,qm); gd<-fit_counts(q,design2,"pathology_score",minimum_replicates=3); gd$series<-"GSE243292"; gd$tissue<-"Frontal_cortex"; gd$comparison<-"A/T_pathology_trend"
allgenes<-rbind(gr,gd); write.csv(allgenes,file.path(outdir,"microglia_gene_contrasts.csv"),row.names=FALSE)

sig<-read.csv(file.path(outdir,"mapped_robust_mouse_signatures.csv"))
test_one<-function(h,series,mouse_tissue){
 s<-sig[sig$tissue==mouse_tissue,]; z<-merge(s,h,by="human_symbol",suffixes=c("_mouse","_human")); rho<-cor(z$meta_logFC,z$logFC,method="spearman"); cp<-cor.test(z$meta_logFC,z$logFC,method="spearman",exact=FALSE)$p.value
 sets<-list(Mouse_ASO_down=s$human_symbol[s$meta_logFC<0],Mouse_ASO_up=s$human_symbol[s$meta_logFC>0]); idx<-ids2indices(sets,h$human_symbol,remove.empty=TRUE); cam<-cameraPR(setNames(h$t,h$human_symbol),idx,use.ranks=TRUE,inter.gene.cor=0.01,sort=FALSE)
 data.frame(series=series,mouse_signature=mouse_tissue,n_mapped=nrow(z),spearman_rho=rho,cor_P=cp,down_N=cam[1,"NGenes"],down_human_direction=cam[1,"Direction"],down_P=cam[1,"PValue"],up_N=cam[2,"NGenes"],up_human_direction=cam[2,"Direction"],up_P=cam[2,"PValue"])
}
tests<-rbind(test_one(gr,"GSE138852","Hippocampus"),test_one(gd,"GSE243292","Frontal_cortex")); adj<-p.adjust(c(tests$down_P,tests$up_P),"BH"); tests$down_FDR<-adj[seq_len(nrow(tests))]; tests$up_FDR<-adj[nrow(tests)+seq_len(nrow(tests))]; tests$expected_complete<-tests$spearman_rho<0 & tests$down_human_direction=="Up" & tests$up_human_direction=="Down"
write.csv(tests,file.path(outdir,"microglia_signature_validation.csv"),row.names=FALSE)

cd33<-allgenes[allgenes$human_symbol=="CD33",c("series","tissue","comparison","logFC","P.Value","adj.P.Val")]; write.csv(cd33,file.path(outdir,"microglia_CD33_results.csv"),row.names=FALSE)

png(file.path(figdir,"microglia_mouse_signature_vs_human_AD.png"),width=1800,height=900,res=180); par(mfrow=c(1,2),mar=c(5,5,3,1))
for(i in seq_len(nrow(tests))){h<-if(i==1)gr else gd; mt<-tests$mouse_signature[i]; s<-sig[sig$tissue==mt,]; z<-merge(s,h,by="human_symbol",suffixes=c("_mouse","_human")); plot(z$meta_logFC,z$logFC,pch=16,col=adjustcolor("#315C8A",.55),xlab="Mouse CD33-ASO meta log2FC",ylab=if(i==1)"Human microglia AD log2FC" else "Human microglia pathology-trend coefficient",main=tests$series[i]); abline(h=0,v=0,col="grey70"); abline(lm(logFC~meta_logFC,z),col="#B24C63",lwd=2); legend("topleft",bty="n",legend=sprintf("rho = %.2f, P = %.3g",tests$spearman_rho[i],tests$cor_P[i]))}; dev.off()

writeLines(c("# Human microglia external validation","","GSE138852 was analysed as six microglial pseudobulk library pools (three AD and three control; two donors per pool). GSE243292 was analysed as a donor-level ordered A/T pathology trend across 15 donors.","","## Signature tests",paste(capture.output(print(tests)),collapse="\n"),"","## CD33",paste(capture.output(print(cd33)),collapse="\n")),file.path(outdir,"MICROGLIA_VALIDATION_SUMMARY.md"))
message("Human microglia external validation completed")
