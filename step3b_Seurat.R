############
## 3. Enrichment
############
.libPaths(c("/home/yzj/R/x86_64-pc-linux-gnu-library/4.0","/home/zy/tools/R-4.0.0/library"))
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
library(xlsx)

get_enrichr <- function(wp,genes,TAG){
  symbol2id=mapIds(org.Hs.eg.db,genes,"ENTREZID",'SYMBOL')
  id=symbol2id[which(symbol2id!='')] #提取出非NA的ENTREZID
  
  if(length(id) > 0){
    #GO BP 富集分析#
    ego <- enrichGO(OrgDb="org.Hs.eg.db", gene = id, ont = "BP", pvalueCutoff = 0.05, readable= TRUE) #GO富集分析
    ego <- clusterProfiler::simplify(ego,cutoff=0.8,by="p.adjust",select_fun=min,measure="Wang")
    ego_res <- as.data.frame(ego)
    if(nrow(ego_res)>0){
      write.xlsx(ego_res,paste(wp,"/GO_BP.xlsx",sep=''),row.names = FALSE,sheetName = TAG,append = TRUE)
    }    
    
    #GO CC 富集分析#
    ego <- enrichGO(OrgDb="org.Hs.eg.db", gene = id, ont = "CC", pvalueCutoff = 0.05, readable= TRUE) #GO富集分析
    ego <- clusterProfiler::simplify(ego,cutoff=0.8,by="p.adjust",select_fun=min,measure="Wang")
    ego_res <- as.data.frame(ego)
    if(nrow(ego_res)>0){
      write.xlsx(ego_res,paste(wp,"/GO_CC.xlsx",sep=''),row.names = FALSE,sheetName = TAG,append = TRUE)
    } 
    
    #GO MF 富集分析#
    ego <- enrichGO(OrgDb="org.Hs.eg.db", gene = id, ont = "MF", pvalueCutoff = 0.05, readable= TRUE) #GO富集分析
    ego <- clusterProfiler::simplify(ego,cutoff=0.8,by="p.adjust",select_fun=min,measure="Wang")
    ego_res <- as.data.frame(ego)
    if(nrow(ego_res)>0){
      write.xlsx(ego_res,paste(wp,"/GO_MF.xlsx",sep=''),row.names = FALSE,sheetName = TAG,append = TRUE)
    } 
    
    #KEGG分析#
    # ekk <- enrichKEGG(gene= id,organism  = 'hsa')	 #KEGG富集分析
    # ekk_res <- as.data.frame(ekk)
    # if(nrow(ekk_res)>0){
    #   write.xlsx(ekk_res,paste(wp,"/KEGG.xlsx",sep=''),row.names = FALSE,sheetName = TAG,append = TRUE)
    # }
  }
}
get_runMK <- function(root,MKlst,FC,adjP){
  #### Create FileDir
  wp <- paste(root,'/','FC',FC,'_Qvalue',adjP,sep='')
  if(!file.exists(file.path(wp))){
    dir.create(file.path(wp),recursive = TRUE)
  }
  CellType <- names(MKlst)
  for(i in 1:length(CellType)){
    CT=CellType[i]
    print(CT)
    TAG=CT
    ################### Find DEG
    MK <- MKlst[[CT]]
    MK <- MK[MK$avg_logFC > log(FC) & MK$p_val_adj < adjP,]
    MK <- tibble::rownames_to_column(MK, "gene")
    MK <- MK[order(MK$avg_logFC,decreasing = T),]
    #write.xlsx(MK,paste(wp,'/MK.xlsx',sep=''),row.names = FALSE, sheetName = TAG,append = TRUE)
    
    ################### Enrichment
    if(nrow(MK) >0){
      get_enrichr(wp=wp,genes = MK$gene,TAG=TAG)
    }
  }
}

###################
## 1.1 Find MK of each cell type in  all cells;
###################
MKlst <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/Markers_celltype.RDS')
root='/home/yzj/JingMA_NEW/res/compCT/All'

FC=1.5
adjP=0.05
get_runMK(root,MKlst,FC,adjP)
print(paste(FC,' DONE!'))


###################
## 1.2 Find MK of each cell type in Chond cells
###################
MKlst <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/Markers_celltype_Chond.RDS')
root='/home/yzj/JingMA_NEW/res/compCT/Chond'

FC=1.5
adjP=0.05
get_runMK(root,MKlst,FC,adjP)
print(paste(FC,' Chond DONE!'))


###################
## 1.3 Find MK of each cell type in Stroma cells
###################
MKlst <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/Markers_celltype_Stroma.RDS')
root='/home/yzj/JingMA_NEW/res/compCT/Stroma'

FC=1.5
adjP=0.05
get_runMK(root,MKlst,FC,adjP)
print(paste(FC,' Stroma DONE!'))

