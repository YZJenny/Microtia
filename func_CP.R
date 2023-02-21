.libPaths(c("/home/yzj/R/x86_64-pc-linux-gnu-library/4.0","/home/zy/tools/R-4.0.0/library"))
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

func_CP <- function(genes){
  genes=unique(RegulonInfo_CSC_up$gene)
  symbol2id=mapIds(org.Hs.eg.db,genes,"ENTREZID",'SYMBOL')
  id=symbol2id[which(symbol2id!='')] #提取出非NA的ENTREZID
  
  #GO BP 富集分析#
  ego <- enrichGO(OrgDb="org.Hs.eg.db", gene = id, ont = "BP", pvalueCutoff = 0.05, readable= TRUE) #GO富集分析
  ego <- clusterProfiler::simplify(ego,cutoff=0.8,by="p.adjust",select_fun=min,measure="Wang")
  ego_res <- as.data.frame(ego)
  return(ego_res)
}
