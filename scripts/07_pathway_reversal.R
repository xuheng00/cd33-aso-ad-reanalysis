#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(limma))
root <- normalizePath(file.path(getwd(), "AD_CD33_reanalysis"))
meta_input <- Sys.getenv("META_INPUT", unset=file.path(root,"results","meta_analysis","fixed_effect_meta_by_region.csv"))
output_name <- Sys.getenv("PATHWAY_OUTPUT_NAME", unset="pathway_reversal")
meta <- read.csv(meta_input,check.names=FALSE)
probe_map <- read.csv(file.path(root,"data","processed","representative_probe_per_gene.csv"),check.names=FALSE)
out_dir <- file.path(root,"results",output_name)
fig_dir <- file.path(root,"figures",output_name)
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE); dir.create(fig_dir,recursive=TRUE,showWarnings=FALSE)

ann <- read.delim(gzfile(file.path(root,"data","raw","GPL6885.annot.gz")),skip=28,check.names=FALSE,quote="",comment.char="",stringsAsFactors=FALSE)
background <- unique(probe_map$gene_symbol)

go_rows <- list(); kk <- 1
for(i in seq_len(nrow(ann))){
  g <- ann$`Gene symbol`[i]; ids <- ann$`GO:Process ID`[i]; terms <- ann$`GO:Process`[i]
  if(is.na(g)||g==""||grepl("///",g,fixed=TRUE)||is.na(ids)||ids=="") next
  idv <- strsplit(ids,"///",fixed=TRUE)[[1]]; tv <- strsplit(terms,"///",fixed=TRUE)[[1]]
  n <- min(length(idv),length(tv)); if(!n) next
  go_rows[[kk]] <- data.frame(pathway=idv[seq_len(n)],name=tv[seq_len(n)],gene=g,stringsAsFactors=FALSE); kk<-kk+1
}
go_map <- unique(do.call(rbind,go_rows)); go_map <- go_map[go_map$gene %in% background & go_map$pathway!="GO:0008150",]

id_map <- unique(ann[,c("Gene ID","Gene symbol")]); names(id_map)<-c("entrez","gene")
id_map <- id_map[id_map$entrez!="" & id_map$gene!="" & !grepl("///",id_map$entrez,fixed=TRUE) & !grepl("///",id_map$gene,fixed=TRUE),]
kl <- read.delim(file.path(root,"data","pathways","kegg_mmu_links.tsv"),header=FALSE,col.names=c("pathway","entrez"),stringsAsFactors=FALSE)
kn <- read.delim(file.path(root,"data","pathways","kegg_mmu_names.tsv"),header=FALSE,col.names=c("pathway","name"),stringsAsFactors=FALSE)
kl$pathway <- sub("^path:","",kl$pathway); kl$entrez<-sub("^mmu:","",kl$entrez); kn$name<-sub(" - Mus musculus.*$","",kn$name)
kegg_map <- merge(merge(kl,kn,by="pathway"),id_map,by="entrez")
kegg_map <- unique(kegg_map[,c("pathway","name","gene")]); kegg_map<-kegg_map[kegg_map$gene%in%background,]

make_sets <- function(map){
  sets<-split(map$gene,map$pathway); sets<-lapply(sets,unique); sets<-sets[lengths(sets)>=10 & lengths(sets)<=500]
  names_lookup<-setNames(map$name,map$pathway)
  list(sets=sets,names=names_lookup)
}
collections<-list(GO_BP=make_sets(go_map),KEGG=make_sets(kegg_map))
write.csv(go_map,file.path(out_dir,"GO_BP_gene_sets.csv"),row.names=FALSE); write.csv(kegg_map,file.path(out_dir,"KEGG_gene_sets.csv"),row.names=FALSE)

run_collection <- function(tissue,collection){
  d<-meta[meta$tissue==tissue & meta$contrast=="Disease",]; tr<-meta[meta$tissue==tissue & meta$contrast=="CD33_specific",]
  dz<-setNames(d$meta_z,d$gene_symbol); tz<-setNames(tr$meta_z,tr$gene_symbol)
  genes<-intersect(names(dz),names(tz)); dz<-dz[genes]; tz<-tz[genes]
  idx<-ids2indices(collection$sets,genes,remove.empty=TRUE)
  rd<-cameraPR(dz,index=idx,use.ranks=TRUE,inter.gene.cor=0.01,sort=FALSE)
  rt<-cameraPR(tz,index=idx,use.ranks=TRUE,inter.gene.cor=0.01,sort=FALSE)
  rd$pathway<-rownames(rd); rt$pathway<-rownames(rt)
  z<-merge(rd,rt,by="pathway",suffixes=c("_disease","_treatment"))
  z$name<-unname(collection$names[z$pathway]); z$tissue<-tissue
  z$strict_reversal<-z$FDR_disease<0.05 & z$FDR_treatment<0.05 & z$Direction_disease!=z$Direction_treatment
  z$nominal_reversal<-z$PValue_disease<0.05 & z$PValue_treatment<0.05 & z$Direction_disease!=z$Direction_treatment
  z$reversal_correlation<-vapply(z$pathway,function(p){g<-intersect(collection$sets[[p]],genes); if(length(g)<5) NA_real_ else suppressWarnings(cor(d$meta_logFC[match(g,d$gene_symbol)],tr$meta_logFC[match(g,tr$gene_symbol)],method="spearman"))},numeric(1))
  z
}

all<-list()
for(cn in names(collections)) for(tissue in c("Hippocampus","Frontal_cortex")){
  z<-run_collection(tissue,collections[[cn]]); z$collection<-cn; all[[paste(cn,tissue)]]<-z
}
pathres<-do.call(rbind,all)
write.csv(pathres,file.path(out_dir,"pathway_reversal_all.csv"),row.names=FALSE)
strict<-pathres[pathres$strict_reversal,]; strict<-strict[order(strict$collection,strict$tissue,pmax(strict$FDR_disease,strict$FDR_treatment)),]
write.csv(strict,file.path(out_dir,"strict_pathway_reversals.csv"),row.names=FALSE)

key<-interaction(pathres$collection,pathres$pathway,drop=TRUE)
cross<-do.call(rbind,lapply(split(pathres,key),function(z){
  if(nrow(z)!=2 || !all(c("Hippocampus","Frontal_cortex")%in%z$tissue)) return(NULL)
  if(!all(z$strict_reversal)) return(NULL)
  if(length(unique(z$Direction_disease))!=1 || length(unique(z$Direction_treatment))!=1) return(NULL)
  data.frame(collection=z$collection[1],pathway=z$pathway[1],name=z$name[1],disease_direction=z$Direction_disease[1],treatment_direction=z$Direction_treatment[1],max_FDR=max(z$FDR_disease,z$FDR_treatment),mean_reversal_correlation=mean(z$reversal_correlation,na.rm=TRUE))
}))
if(is.null(cross)) cross<-data.frame()
if(nrow(cross)) cross<-cross[order(cross$max_FDR),]
write.csv(cross,file.path(out_dir,"cross_region_strict_reversals.csv"),row.names=FALSE)

counts<-aggregate(cbind(strict=pathres$strict_reversal,nominal=pathres$nominal_reversal)~collection+tissue,pathres,sum)
write.csv(counts,file.path(out_dir,"pathway_reversal_counts.csv"),row.names=FALSE)

plotset<-head(strict,30)
if(nrow(plotset)){
  plotset$score<- -log10(pmax(plotset$FDR_disease,plotset$FDR_treatment)); plotset$label<-paste(plotset$collection,plotset$name)
  png(file.path(fig_dir,"top_strict_pathway_reversals.png"),width=1800,height=1400,res=180)
  par(mar=c(5,18,3,1)); ord<-order(plotset$score); barplot(plotset$score[ord],names.arg=plotset$label[ord],horiz=TRUE,las=1,col=ifelse(plotset$Direction_disease[ord]=="Up","#B24C63","#3B7EA1"),xlab="-log10(max FDR of disease and treatment)",main="Strict pathway reversals"); dev.off()
}

writeLines(c("# Pathway-level reversal test","",
             "A pathway is a strict reversal when rank-based competitive tests show FDR < 0.05 for both the disease and CD33-specific treatment contrasts, with opposite directions.","",
             "## Counts","",paste(capture.output(print(counts)),collapse="\n"),"",
             paste0("Cross-region strict reversals: ",nrow(cross),"."),"",
             "## Cross-region results","",paste(capture.output(print(head(cross,30))),collapse="\n")),
           file.path(out_dir,"PATHWAY_REVERSAL_SUMMARY.md"))
message("Pathway reversal analysis completed")
