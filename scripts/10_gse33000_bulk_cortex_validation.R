#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(limma))
root <- normalizePath(file.path(getwd(),"AD_CD33_reanalysis"))
rawdir<-file.path(root,"data","human_raw"); outdir<-file.path(root,"results","human_validation")

read_geo <- function(path){
 z<-readLines(gzfile(path)); s<-grep("!series_matrix_table_begin",z,fixed=TRUE); e<-grep("!series_matrix_table_end",z,fixed=TRUE)
 tab<-read.delim(textConnection(z[(s+1):(e-1)]),check.names=FALSE,quote="\"",stringsAsFactors=FALSE)
 ids<-tab[[1]]; mat<-as.matrix(tab[,-1]); storage.mode(mat)<-"numeric"; rownames(mat)<-ids
 m<-z[seq_len(s-1)]; get<-function(tag){i<-grep(paste0("^",tag,"\\t"),m); lapply(i,function(j)gsub('^"|"$','',strsplit(m[j],"\\t")[[1]][-1]))}
 list(expr=mat,gsm=get("!Sample_geo_accession")[[1]],title=get("!Sample_title")[[1]],chars2=get("!Sample_characteristics_ch2"))
}
field<-function(chars,key,n){for(v in chars){h<-grepl(paste0("^",key,"\\s*:"),v,ignore.case=TRUE); if(sum(h)>n/2)return(trimws(sub("^[^:]+:\\s*","",v)))};rep(NA_character_,n)}
g<-read_geo(file.path(rawdir,"GSE33000_series_matrix.txt.gz")); n<-ncol(g$expr)
meta<-data.frame(series="GSE33000",gsm=g$gsm,title=g$title,age=as.numeric(sub("[^0-9.].*$","",field(g$chars2,"age",n))),sex=tolower(field(g$chars2,"gender",n)),status=field(g$chars2,"disease status",n),stringsAsFactors=FALSE)
meta$status<-ifelse(grepl("Alzheimer",meta$status,ignore.case=TRUE),"AD",ifelse(grepl("non-demented|control",meta$status,ignore.case=TRUE),"Control","Other"))
meta$include<-meta$status%in%c("AD","Control") & is.finite(meta$age) & !is.na(meta$sex)
write.csv(meta,file.path(root,"data","human_processed","GSE33000_sample_manifest.csv"),row.names=FALSE)

az<-readLines(gzfile(file.path(rawdir,"GPL4372.annot.gz"))); s<-grep("!platform_table_begin",az,fixed=TRUE); e<-grep("!platform_table_end",az,fixed=TRUE)
ann<-read.delim(textConnection(az[(s+1):(e-1)]),check.names=FALSE,quote="",comment.char="",stringsAsFactors=FALSE)[,c("ID","Gene symbol")]; names(ann)<-c("probe","human_symbol")
ann<-ann[ann$human_symbol!="" & ann$human_symbol!="---" & !grepl("///",ann$human_symbol,fixed=TRUE),]
a<-ann[match(rownames(g$expr),ann$probe),]; a$variance<-apply(g$expr,1,var,na.rm=TRUE); a<-a[!is.na(a$human_symbol),]; a<-a[order(a$human_symbol,-a$variance),]; a<-a[!duplicated(a$human_symbol),]
y<-g$expr[match(a$probe,rownames(g$expr)),meta$include,drop=FALSE]; rownames(y)<-a$human_symbol
mm<-droplevels(meta[meta$include,]); mm$status<-relevel(factor(mm$status),"Control"); mm$sex<-factor(mm$sex)
design<-model.matrix(~status+age+sex,data=mm); fit<-eBayes(lmFit(y,design),trend=TRUE,robust=TRUE)
tt<-topTable(fit,coef="statusAD",number=Inf,sort.by="none"); tt$human_symbol<-rownames(tt); tt$series<-"GSE33000"; tt$tissue<-"Frontal_cortex"; tt$n_AD<-sum(mm$status=="AD"); tt$n_Control<-sum(mm$status=="Control"); tt$SE<-abs(tt$logFC/tt$t)
write.csv(tt,file.path(outdir,"GSE33000_gene_contrasts.csv"),row.names=FALSE)

sig<-read.csv(file.path(outdir,"mapped_robust_mouse_signatures.csv")); sig<-sig[sig$tissue=="Frontal_cortex",]
z<-merge(sig,tt,by="human_symbol",suffixes=c("_mouse","_human")); rho<-cor(z$meta_logFC,z$logFC,method="spearman"); cp<-cor.test(z$meta_logFC,z$logFC,method="spearman",exact=FALSE)$p.value
sets<-list(Mouse_ASO_down=sig$human_symbol[sig$meta_logFC<0],Mouse_ASO_up=sig$human_symbol[sig$meta_logFC>0]); idx<-ids2indices(sets,tt$human_symbol,remove.empty=TRUE)
cam<-cameraPR(setNames(tt$t,tt$human_symbol),idx,use.ranks=TRUE,inter.gene.cor=0.01,sort=FALSE)
test<-data.frame(series="GSE33000",tissue="Frontal_cortex",n_mapped=nrow(z),spearman_rho=rho,cor_P=cp,down_N=cam[1,"NGenes"],down_human_direction=cam[1,"Direction"],down_P=cam[1,"PValue"],up_N=cam[2,"NGenes"],up_human_direction=cam[2,"Direction"],up_P=cam[2,"PValue"])
write.csv(test,file.path(outdir,"GSE33000_signature_validation.csv"),row.names=FALSE)

# Primary bulk-cortex human meta-analysis: GSE48350 SFG plus GSE33000 DLPFC.
old<-read.csv(file.path(outdir,"human_region_gene_contrasts.csv")); a1<-old[old$series=="GSE48350" & old$tissue=="Frontal_cortex",]; a2<-tt
m<-merge(a1,a2,by="human_symbol",suffixes=c("_48350","_33000")); w1<-1/m$SE_48350^2; w2<-1/m$SE_33000^2
m$meta_logFC<-(w1*m$logFC_48350+w2*m$logFC_33000)/(w1+w2); m$meta_SE<-sqrt(1/(w1+w2)); m$meta_z<-m$meta_logFC/m$meta_SE; m$meta_P<-2*pnorm(-abs(m$meta_z)); m$meta_FDR<-p.adjust(m$meta_P,"BH")
write.csv(m,file.path(outdir,"primary_bulk_frontal_human_meta.csv"),row.names=FALSE)
zz<-merge(sig,m,by="human_symbol",suffixes=c("_mouse","_human")); mrho<-cor(zz$meta_logFC_mouse,zz$meta_logFC_human,method="spearman"); mcp<-cor.test(zz$meta_logFC_mouse,zz$meta_logFC_human,method="spearman",exact=FALSE)$p.value
midx<-ids2indices(sets,m$human_symbol,remove.empty=TRUE); mcam<-cameraPR(setNames(m$meta_z,m$human_symbol),midx,use.ranks=TRUE,inter.gene.cor=0.01,sort=FALSE)
mtest<-data.frame(analysis="Primary_bulk_frontal_meta",n_mapped=nrow(zz),spearman_rho=mrho,cor_P=mcp,down_N=mcam[1,"NGenes"],down_human_direction=mcam[1,"Direction"],down_P=mcam[1,"PValue"],up_N=mcam[2,"NGenes"],up_human_direction=mcam[2,"Direction"],up_P=mcam[2,"PValue"])
write.csv(mtest,file.path(outdir,"primary_bulk_frontal_signature_validation.csv"),row.names=FALSE)
png(file.path(root,"figures","human_validation","primary_bulk_frontal_mouse_vs_human.png"),width=1200,height=1000,res=180)
par(mar=c(5,5,3,1)); plot(zz$meta_logFC_mouse,zz$meta_logFC_human,pch=16,col=adjustcolor("#315C8A",0.55),xlab="Mouse CD33-ASO meta log2FC",ylab="Human AD bulk-cortex meta log2FC",main="Cross-species direction test: frontal cortex"); abline(h=0,v=0,col="grey70"); abline(lm(meta_logFC_human~meta_logFC_mouse,zz),col="#B24C63",lwd=2); legend("topleft",bty="n",legend=sprintf("Spearman rho = %.2f, P = %.3g",mrho,mcp)); dev.off()
writeLines(c("# GSE33000 and primary bulk-cortex validation","",sprintf("Included: %d AD and %d controls.",sum(mm$status=="AD"),sum(mm$status=="Control")),"","## GSE33000",paste(capture.output(print(test)),collapse="\n"),"","## Primary bulk frontal meta-analysis",paste(capture.output(print(mtest)),collapse="\n")),file.path(outdir,"GSE33000_AND_PRIMARY_BULK_SUMMARY.md"))
message("GSE33000 validation completed")
