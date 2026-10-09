#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(limma))
root <- normalizePath(file.path(getwd(), "AD_CD33_reanalysis"))
rawdir <- file.path(root,"data","human_raw")
outdir <- file.path(root,"results","human_validation")
figdir <- file.path(root,"figures","human_validation")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE); dir.create(figdir,recursive=TRUE,showWarnings=FALSE)

read_geo_matrix <- function(path){
  z <- readLines(gzfile(path))
  s <- grep("!series_matrix_table_begin",z,fixed=TRUE)
  e <- grep("!series_matrix_table_end",z,fixed=TRUE)
  tab <- read.delim(textConnection(z[(s+1):(e-1)]),check.names=FALSE,quote="\"",stringsAsFactors=FALSE)
  ids <- tab[[1]]; mat <- as.matrix(tab[,-1]); storage.mode(mat)<-"numeric"; rownames(mat)<-ids
  meta_lines <- z[seq_len(s-1)]
  get_tag <- function(tag){
    i <- grep(paste0("^",tag,"\\t"),meta_lines)
    if(!length(i)) return(NULL)
    lapply(i,function(j){v<-strsplit(meta_lines[j],"\\t",fixed=FALSE)[[1]][-1]; gsub('^"|"$','',v)})
  }
  list(expr=mat,title=get_tag("!Sample_title")[[1]],gsm=get_tag("!Sample_geo_accession")[[1]],chars=get_tag("!Sample_characteristics_ch1"),source=get_tag("!Sample_source_name_ch1")[[1]])
}

field_from_chars <- function(chars,key,n){
  for(v in chars){
    hit <- grepl(paste0("^",key,"\\s*:"),v,ignore.case=TRUE)
    if(sum(hit,na.rm=TRUE)>n/2) return(trimws(sub("^[^:]+:\\s*","",v)))
  }
  rep(NA_character_,n)
}

g483 <- read_geo_matrix(file.path(rawdir,"GSE48350_series_matrix.txt.gz"))
n <- ncol(g483$expr)
m483 <- data.frame(series="GSE48350",gsm=g483$gsm,title=g483$title,stringsAsFactors=FALSE)
m483$region <- ifelse(grepl("hippocampus",m483$title,ignore.case=TRUE),"Hippocampus",ifelse(grepl("superior ?frontal",m483$title,ignore.case=TRUE),"Frontal_cortex",NA))
m483$sex <- ifelse(grepl("female",m483$title,ignore.case=TRUE),"female",ifelse(grepl("male",m483$title,ignore.case=TRUE),"male",NA))
m483$age <- suppressWarnings(as.numeric(sub(".*_(female|male)_([0-9]+)(yrs)?_.*","\\2",m483$title,ignore.case=TRUE)))
m483$status <- ifelse(grepl("_AD_",m483$title,ignore.case=TRUE),"AD","Control")
# Prespecified age restriction: AD range begins at 60 in this series.
m483$include <- !is.na(m483$region) & !is.na(m483$age) & (m483$status=="AD" | m483$age>=60)

g528 <- read_geo_matrix(file.path(rawdir,"GSE5281_series_matrix.txt.gz"))
n <- ncol(g528$expr)
m528 <- data.frame(series="GSE5281",gsm=g528$gsm,title=g528$title,stringsAsFactors=FALSE)
m528$status <- ifelse(grepl("affected",m528$title,ignore.case=TRUE),"AD","Control")
m528$region_code <- sub("[_ ].*","",m528$title)
m528$region <- ifelse(m528$region_code=="HIP","Hippocampus",ifelse(m528$region_code=="SFG","Frontal_cortex",NA))
m528$sex <- tolower(field_from_chars(g528$chars,"sex",n))
age_text <- field_from_chars(g528$chars,"age",n)
m528$age <- suppressWarnings(as.numeric(sub("^>?([0-9.]+).*$","\\1",age_text)))
m528$include <- !is.na(m528$region) & !is.na(m528$age)

write.csv(rbind(m483,m528[,names(m483)]),file.path(root,"data","human_processed","human_sample_manifest.csv"),row.names=FALSE)

# Official GPL570 annotation and one representative probe per human gene.
az <- readLines(gzfile(file.path(rawdir,"GPL570.annot.gz")))
as <- grep("!platform_table_begin",az,fixed=TRUE); ae <- grep("!platform_table_end",az,fixed=TRUE)
ann <- read.delim(textConnection(az[(as+1):(ae-1)]),check.names=FALSE,quote="",comment.char="",stringsAsFactors=FALSE)
ann <- ann[,c("ID","Gene symbol")]; names(ann)<-c("probe","human_symbol")
ann <- ann[ann$human_symbol!="" & ann$human_symbol!="---" & !grepl("///",ann$human_symbol,fixed=TRUE),]

prep_gene <- function(expr){
  if(quantile(expr,0.99,na.rm=TRUE)>20 || quantile(expr,0.01,na.rm=TRUE)>=0) expr <- log2(pmax(expr,0.01))
  a <- ann[match(rownames(expr),ann$probe),]; a$mean <- rowMeans(expr,na.rm=TRUE)
  a <- a[!is.na(a$human_symbol),]; a <- a[order(a$human_symbol,-a$mean),]; a <- a[!duplicated(a$human_symbol),]
  out <- expr[match(a$probe,rownames(expr)),,drop=FALSE]; rownames(out)<-a$human_symbol
  out
}
e483<-prep_gene(g483$expr); e528<-prep_gene(g528$expr)

fit_region <- function(expr,meta,series,tissue){
  keep <- meta$include & meta$region==tissue
  mm <- droplevels(meta[keep,]); y <- expr[,keep,drop=FALSE]
  mm$status <- relevel(factor(mm$status),"Control"); mm$sex<-factor(mm$sex)
  covars <- c("status",if(length(unique(mm$age[is.finite(mm$age)]))>1) "age",if(nlevels(mm$sex)>1) "sex")
  design <- model.matrix(as.formula(paste("~",paste(covars,collapse="+"))),data=mm)
  fit <- eBayes(lmFit(y,design),trend=TRUE,robust=TRUE)
  coefname <- grep("^statusAD$",colnames(design),value=TRUE)
  tt <- topTable(fit,coef=coefname,number=Inf,sort.by="none")
  tt$human_symbol<-rownames(tt); tt$series<-series; tt$tissue<-tissue
  tt$n_AD<-sum(mm$status=="AD"); tt$n_Control<-sum(mm$status=="Control")
  tt$SE<-abs(tt$logFC/tt$t)
  tt
}
human <- do.call(rbind,list(fit_region(e483,m483,"GSE48350","Hippocampus"),fit_region(e483,m483,"GSE48350","Frontal_cortex"),fit_region(e528,m528,"GSE5281","Hippocampus"),fit_region(e528,m528,"GSE5281","Frontal_cortex")))
write.csv(human,file.path(outdir,"human_region_gene_contrasts.csv"),row.names=FALSE)

# MGI one-to-one mouse-human homology map.
hom <- read.delim(file.path(rawdir,"HOM_MouseHumanSequence.rpt"),check.names=FALSE,stringsAsFactors=FALSE)
mouse <- hom[hom$`NCBI Taxon ID`==10090,c("DB Class Key","Symbol")]; names(mouse)<-c("class","mouse_symbol")
humanmap <- hom[hom$`NCBI Taxon ID`==9606,c("DB Class Key","Symbol")]; names(humanmap)<-c("class","human_symbol")
hm <- merge(mouse,humanmap,by="class")
hm <- hm[!duplicated(hm$class) & !duplicated(hm$mouse_symbol) & !duplicated(hm$human_symbol),]
write.csv(hm,file.path(root,"data","human_processed","MGI_one_to_one_mouse_human.csv"),row.names=FALSE)

mousemeta <- read.csv(file.path(root,"results","meta_analysis","fixed_effect_meta_by_region.csv"))
sig <- mousemeta[mousemeta$contrast=="CD33_specific" & mousemeta$meta_FDR<0.05 & mousemeta$direction_concordant & mousemeta$I2<50,]
sig <- merge(sig,hm,by.x="gene_symbol",by.y="mouse_symbol")
write.csv(sig,file.path(outdir,"mapped_robust_mouse_signatures.csv"),row.names=FALSE)

test_signature <- function(h,tissue,series){
  s <- sig[sig$tissue==tissue,]; h <- h[h$tissue==tissue & h$series==series,]
  z <- merge(s,h,by="human_symbol",suffixes=c("_mouse","_human"))
  rho <- cor(z$meta_logFC,z$logFC,method="spearman",use="complete.obs")
  cp <- suppressWarnings(cor.test(z$meta_logFC,z$logFC,method="spearman",exact=FALSE)$p.value)
  universe <- setNames(h$t,h$human_symbol)
  sets <- list(Mouse_ASO_down=s$human_symbol[s$meta_logFC<0],Mouse_ASO_up=s$human_symbol[s$meta_logFC>0])
  idx <- ids2indices(sets,names(universe),remove.empty=TRUE)
  cam <- cameraPR(universe,index=idx,use.ranks=TRUE,inter.gene.cor=0.01,sort=FALSE)
  data.frame(series=series,tissue=tissue,n_mapped=nrow(z),spearman_rho=rho,cor_P=cp,
             down_N=cam["Mouse_ASO_down","NGenes"],down_human_direction=cam["Mouse_ASO_down","Direction"],down_P=cam["Mouse_ASO_down","PValue"],
             up_N=cam["Mouse_ASO_up","NGenes"],up_human_direction=cam["Mouse_ASO_up","Direction"],up_P=cam["Mouse_ASO_up","PValue"],stringsAsFactors=FALSE)
}
tests <- do.call(rbind,lapply(c("GSE48350","GSE5281"),function(series) do.call(rbind,lapply(c("Hippocampus","Frontal_cortex"),function(tissue)test_signature(human,tissue,series)))))
tests$direction_support <- tests$spearman_rho<0
tests$down_expected <- tests$down_human_direction=="Up"
tests$up_expected <- tests$up_human_direction=="Down"
padj <- p.adjust(c(tests$down_P,tests$up_P),"BH")
tests$down_FDR <- padj[seq_len(nrow(tests))]
tests$up_FDR <- padj[nrow(tests)+seq_len(nrow(tests))]
write.csv(tests,file.path(outdir,"cohort_signature_validation.csv"),row.names=FALSE)

# Human fixed-effect meta-analysis across the two cohorts, region-specific.
meta_two <- function(tissue){
  a<-human[human$tissue==tissue & human$series=="GSE48350",]; b<-human[human$tissue==tissue & human$series=="GSE5281",]
  m<-merge(a,b,by="human_symbol",suffixes=c("_48350","_5281")); w1<-1/m$SE_48350^2; w2<-1/m$SE_5281^2
  m$meta_logFC<-(w1*m$logFC_48350+w2*m$logFC_5281)/(w1+w2); m$meta_SE<-sqrt(1/(w1+w2)); m$meta_z<-m$meta_logFC/m$meta_SE
  m$meta_P<-2*pnorm(-abs(m$meta_z)); m$meta_FDR<-p.adjust(m$meta_P,"BH"); m$tissue<-tissue; m
}
hmeta<-do.call(rbind,lapply(c("Hippocampus","Frontal_cortex"),meta_two)); write.csv(hmeta,file.path(outdir,"human_ad_fixed_effect_meta.csv"),row.names=FALSE)

meta_tests<-do.call(rbind,lapply(c("Hippocampus","Frontal_cortex"),function(tissue){
  s<-sig[sig$tissue==tissue,]; h<-hmeta[hmeta$tissue==tissue,]; z<-merge(s,h,by="human_symbol")
  rho<-cor(z$meta_logFC.x,z$meta_logFC.y,method="spearman"); cp<-cor.test(z$meta_logFC.x,z$meta_logFC.y,method="spearman",exact=FALSE)$p.value
  sets<-list(Mouse_ASO_down=s$human_symbol[s$meta_logFC<0],Mouse_ASO_up=s$human_symbol[s$meta_logFC>0]); idx<-ids2indices(sets,h$human_symbol,remove.empty=TRUE)
  cam<-cameraPR(setNames(h$meta_z,h$human_symbol),idx,use.ranks=TRUE,inter.gene.cor=0.01,sort=FALSE)
  data.frame(tissue=tissue,n_mapped=nrow(z),spearman_rho=rho,cor_P=cp,down_N=cam[1,"NGenes"],down_human_direction=cam[1,"Direction"],down_P=cam[1,"PValue"],up_N=cam[2,"NGenes"],up_human_direction=cam[2,"Direction"],up_P=cam[2,"PValue"])
}))
meta_tests$all_expected <- meta_tests$spearman_rho<0 & meta_tests$down_human_direction=="Up" & meta_tests$up_human_direction=="Down"
mpadj <- p.adjust(c(meta_tests$down_P,meta_tests$up_P),"BH")
meta_tests$down_FDR <- mpadj[seq_len(nrow(meta_tests))]
meta_tests$up_FDR <- mpadj[nrow(meta_tests)+seq_len(nrow(meta_tests))]
write.csv(meta_tests,file.path(outdir,"meta_signature_validation.csv"),row.names=FALSE)

# Compact validation figure.
png(file.path(figdir,"mouse_signature_vs_human_AD.png"),width=1800,height=900,res=180)
par(mfrow=c(1,2),mar=c(5,5,3,1))
for(tissue in c("Hippocampus","Frontal_cortex")){
 s<-sig[sig$tissue==tissue,]; h<-hmeta[hmeta$tissue==tissue,]; z<-merge(s,h,by="human_symbol")
 plot(z$meta_logFC.x,z$meta_logFC.y,pch=16,col=adjustcolor("#315C8A",0.55),xlab="Mouse CD33-ASO meta log2FC",ylab="Human AD meta log2FC",main=gsub("_"," ",tissue)); abline(h=0,v=0,col="grey70"); abline(lm(meta_logFC.y~meta_logFC.x,z),col="#B24C63",lwd=2)
 legend("topright",bty="n",legend=sprintf("Spearman rho = %.2f",cor(z$meta_logFC.x,z$meta_logFC.y,method="spearman")))
}
dev.off()

counts <- with(rbind(m483[m483$include,],m528[m528$include,names(m483)]),as.data.frame(table(series,region,status)))
write.csv(counts,file.path(outdir,"human_sample_counts.csv"),row.names=FALSE)
writeLines(c("# Human AD external validation","","## Included sample counts","",paste(capture.output(print(counts)),collapse="\n"),"","## Cohort-level validation","",paste(capture.output(print(tests)),collapse="\n"),"","## Cross-cohort human meta-validation","",paste(capture.output(print(meta_tests)),collapse="\n")),file.path(outdir,"HUMAN_VALIDATION_SUMMARY.md"))
message("Human AD external validation completed")
