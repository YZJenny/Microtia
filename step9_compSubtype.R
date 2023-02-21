.libPaths(c('/home/yzj/R/x86_64-pc-linux-gnu-library/4.0',
            "/home/zy/tools/R-4.0.0/library"))
library(Seurat)
library(pheatmap)
data(c2BroadSets)
library(Biobase)
library(genefilter)
library(limma)
library(RColorBrewer)
library(AnnotationHub)	
library(org.Hs.eg.db)   #人类注释数据库
library(clusterProfiler)
library(DOSE)
library(dplyr)
library(tidyverse)
library(reshape2)

get_enrichr <- function(wp,DEG,TAG){
  gene_erichment_results <- list()
  testgenes<- DEG
  symbol2id=mapIds(org.Hs.eg.db,testgenes,"ENTREZID",'SYMBOL')
  id=symbol2id[which(symbol2id!='')] #提取出非NA的ENTREZID
  
  if(length(id) > 0){
    #GO BP 富集分析#
    ego <- enrichGO(OrgDb="org.Hs.eg.db", gene = id, ont = "BP", pvalueCutoff = 0.05, readable= TRUE) #GO富集分析
    ego <- clusterProfiler::simplify(ego,cutoff=0.7,by="p.adjust",select_fun=min,measure="Wang")
    gene_erichment_results[['BP']] <- ego
    ego_res <- as.data.frame(ego)
    write.csv(ego_res,paste(wp,'/GO/',TAG,'_BP.csv',sep=''),quote=F,row.names = F)
    
    #KEGG分析#
    ekk <- enrichKEGG(gene= id,organism  = 'hsa')	 #KEGG富集分析
    gene_erichment_results[['KEGG']] <- ekk
    ekk_res <- as.data.frame(ekk)
    write.csv(ekk_res,paste(wp,'/KEGG/',TAG,'.csv',sep=''),quote=F,row.names = F)
    
    RES <- list(res=gene_erichment_results)
  }
  return(RES)
}

get_barplot <- function(wp,TAG,topN){
  enrich_res <- readRDS(paste(wp,'/DEG/',TAG,'_enrich_res.RDS',sep=''))
  ego <- enrich_res$res[['BP']]
  ekk <- enrich_res$res[['KEGG']]
  
  NUM=nrow(ego)
  if(0 < NUM & NUM < 20){
    height=8
  } else if (NUM >= 20){
    height=16
  }
  
  if(NUM >=1){
    BP <- barplot(ego, showCategory=topN,title="GO Enrichment")
    pdf(paste(wp,'/GO/',TAG,'_BP.pdf',sep=''),width = 14,height = height)
    print(BP)
    dev.off()
  }
  
  NUM=nrow(ekk)
  if(0 < NUM & NUM < 20){
    height=8
  } else if(20 <= NUM & NUM < 30){
    height=10
  }else if(NUM >= 30){
    height=16
  }
  if(NUM >=1){
    KEGG <- barplot(ekk, showCategory=topN,title="KEGG Enrichment")
    pdf(paste(wp,'/KEGG/',TAG,'.pdf',sep=''),width = 14,height = height)
    print(KEGG)
    dev.off()
  }
  
}


#### Find MK
pbmc <- readRDS('/home/yzj/JingMA/res/Harmony/ALL/RDS/PBMC_harmony_noUNK_ORI.RDS')
resol='0.6_new'
index <- which(colnames(pbmc@meta.data)==paste('RNA_snn_res.',resol,sep=''))

CT1  <- '0'
CT2 <- '1'
CT3 <- '2'
subpbmc <- subset(pbmc,idents = c(CT1,CT2,CT3))
DimPlot(subpbmc)

FC=1.5
MK <- FindAllMarkers(subpbmc,test.use = 'MAST',only.pos = T,logfc.threshold = log(FC))
saveRDS(MK,paste('JingMA/res/Harmony/ALL/RDS/Markers_',resol,'_',CT1,'_',CT2,'_',CT3,'_',FC,'.RDS',sep=''))
####


####. Enrichment: compare RegC/C vs EC in all cells
CT1  <- '0'
CT2 <- '1'
CT3 <- '2'

FC=1.5
Qvalue=0.01
topN=50
resol='0.6_new'

root=paste('JingMA/res/EnrichR_',resol,'/Phylog_',CT1,'_',CT2,'_',CT3,sep='')
wp <- paste(root,'/','FC',FC,'_Qvalue',Qvalue,sep='')
print(wp)
Dir <- c('DEG','GO','KEGG')
for(d in Dir){
  if(!file.exists(file.path(wp, d))){
    dir.create(file.path(wp, d),recursive = TRUE)
  }
}

################### Input MK
DEG <- readRDS(paste('JingMA/res/Harmony/ALL/RDS/Markers_',resol,'_',CT1,'_',CT2,'_',CT3,'_',FC,'.RDS',sep=''))

DEG <- DEG[DEG$avg_logFC > log(FC) & DEG$p_val_adj < Qvalue,]
DEG <- DEG[,c(7,6,2,5,1,3:4)]
DEG <- DEG[order(DEG$avg_logFC,decreasing = T),]
write.csv(DEG,paste(wp,'/DEG/markers.csv',sep=''),row.names = FALSE,quote = FALSE)

################### Enrichment
Type=c(CT1,CT2,CT3)

#DEG <- read.table(paste(wp,'/DEG/markers.csv',sep=''))
for(k in 1:length(Type)){
  type=Type[k]
  print(type)
  markers <- DEG[DEG$cluster==type,]
  if(nrow(markers) > 0){
    enrich_res <- get_enrichr(wp=wp,DEG=markers$gene,TAG=type)
    saveRDS(enrich_res,paste(wp,'/DEG/',type,'_enrich_res.RDS',sep=''))
    
    #enrich_res <- readRDS(paste(wp,'/DEG/',type,'_enrich_res.RDS',sep=''))
    get_barplot(wp=wp,TAG = type,topN = topN)
  }
}



