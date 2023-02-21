#########
## 比较正常儿童和患病儿童(C4/C6 vs M1/M2/M3)
#########
library(Seurat)
library(ggplot2)
library(ggrepel)

###################
## 1. Find DEG of Chond cell type only in children samples
###################
setwd('/local/yzj')
library(future)
plan("multiprocess", workers = 60)


pbmc_chond <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype_Chond.Rdata')
Idents(pbmc_chond) <- pbmc_chond$batch
pbmc_chond <- subset(pbmc_chond,idents = c('C4','C6','M1','M2','M3'))

Idents(pbmc_chond) <- pbmc_chond$celltype
celltypes <- unique(pbmc_chond$celltype)

MK.lst <- list()
for (CT in celltypes){
  print(paste('FindDEG:',CT,'Start!'))
  subpbmc <- subset(pbmc_chond,idents = CT)
  Idents(subpbmc) <- subpbmc$type
  markers <- FindMarkers(subpbmc, ident.1 = 'Microtia',only.pos = FALSE)
  MK.lst[[as.character(CT)]] <- markers
}
saveRDS(MK.lst, file = 'JingMA_NEW/res/Harmony/ALL/RDS/DEGs_inChond_inChildren_NormalMicrotia.RDS')

### 1.2 GSEA: no filtering gene
MK.lst <- list()
for (CT in celltypes){
  print(paste('FindDEG:',CT,'Start!'))
  subpbmc <- subset(pbmc_chond,idents = CT)
  Idents(subpbmc) <- subpbmc$type
  markers <- FindMarkers(subpbmc, ident.1 = 'Microtia',only.pos = FALSE,
                         logfc.threshold = 0,min.pct = 0,min.cells.feature = 0)
  MK.lst[[as.character(CT)]] <- markers
}
saveRDS(MK.lst, file = 'JingMA_NEW/res/Harmony/ALL/RDS/DEGs_inChond_inChildren_NormalMicrotia_GSEA.RDS')


################
## step2. pick DEG
################
library(xlsx)
library(Seurat)
library(ggplot2)

batch1='Microtia'
batch2='Normal'

pbmc_chond <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype_Chond.Rdata')
Idents(pbmc_chond) <- pbmc_chond$batch
subpbmc <- subset(pbmc_chond,idents = c('C4','C6','M1','M2','M3'))
Idents(subpbmc) <- subpbmc$celltype

MK.lst <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/DEGs_inChond_inChildren_NormalMicrotia.RDS')
wp=paste('JingMA_NEW/res/compMicrotia/',batch1,'vs',batch2,'_inChildren',sep='')

FC <- 2
adjP <- 0.01

CellTypes=levels(subpbmc$celltype)
upNum <- c()
dnNum <- c()
for(i in 1:length(CellTypes)){
  CT=CellTypes[i]
  print(CT)
  
  op=paste(wp,'/DEG/FC',FC,'_adjP',adjP,'/',sep='')
  ## mkdir
  if(!file.exists(file.path(op))){
    dir.create(file.path(op),recursive = TRUE)
  }
  ##
  data<- MK.lst[[CT]]
  data <- tibble::rownames_to_column(data,'gene')
  
  ## 相对于患病来说
  sigUP_gene <- data[data$p_val_adj  < adjP & data$avg_logFC > log(FC),]
  sigDN_gene <- data[data$p_val_adj  < adjP & data$avg_logFC < -log(FC),]
  
  sigUP_gene <- sigUP_gene[order(sigUP_gene$avg_logFC,decreasing = T),]
  sigDN_gene <- sigDN_gene[order(sigDN_gene$avg_logFC,decreasing = F),]
  
  upNum <- c(upNum,nrow(sigUP_gene))
  dnNum <- c(dnNum,nrow(sigDN_gene))
  ## 相对于成年人来说,上调的基因
  write.xlsx(sigUP_gene, file = paste(op,CT,'.xlsx',sep=''), row.names = FALSE, sheetName = "UP")
  ## 相对于成年人来说,下调的基因
  write.xlsx(sigDN_gene, file = paste(op,CT,'.xlsx',sep=''), row.names = FALSE, sheetName = "DN",append=TRUE)
}

Num <- data.frame(CellType=CellTypes,upGene=upNum,dnGene=dnNum)
Num$CellType <- factor(Num$CellType,levels = c("CSC","TC","C1","C2"))
library(reshape2)
NUM_meltdf <- melt(Num)

p <- ggplot(data=NUM_meltdf, mapping=aes(x=CellType,y=value,fill=variable))+
  geom_bar(stat="identity",width=0.5,position='dodge')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))+
  scale_fill_manual(values=c('#999999','#E69F00'))+theme_classic()
ggsave(paste(op,'DEGnum.pdf',sep=''),p,height = 4,width = 6)


#####
# step3. Clusterprofiler
#####
.libPaths(c("/home/yzj/R/x86_64-pc-linux-gnu-library/4.0",'/home/zy/tools/R-4.0.0/library'))
library(xlsx)
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

batch1='Microtia'
batch2='Normal'

MK.lst <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/DEGs_inChond_inChildren_NormalMicrotia.RDS')
wp=paste('JingMA_NEW/res/compMicrotia/',batch1,'vs',batch2,'_inChildren',sep='')

FC <- 1.5
adjP <- 0.01

CellTypes=c('CSC','TC','C1','C2')
for(i in 1:length(CellTypes)){
  CT <- CellTypes[i]
  print(CT)
  
  op=paste(wp,'/ClusterPro/FC',FC,'_adjP',adjP,'/',sep='')
  ## mkdir
  if(!file.exists(file.path(op))){
    dir.create(file.path(op),recursive = TRUE)
  }
  
  TYPE <- c('UP','DN')
  for(j in 1:length(TYPE)){
    type=TYPE[j]
    sigG_lst <- read.xlsx(paste(wp,'/DEG/FC',FC,'_adjP',adjP,'/',CT,'.xlsx',sep=''),sheetName = type)[,1]
    testgenes<- sigG_lst
    symbol2id=mapIds(org.Hs.eg.db,testgenes,"ENTREZID",'SYMBOL')
    id=symbol2id[which(symbol2id!='')] #提取出非NA的ENTREZID
    
    if(length(id) > 0){
      #GO BP 富集分析#
      ego <- enrichGO(OrgDb="org.Hs.eg.db", gene = id, ont = "BP", pvalueCutoff = 0.05, readable= TRUE) #GO富集分析
      ego <- clusterProfiler::simplify(ego,cutoff=0.8,by="p.adjust",select_fun=min,measure="Wang")
      ego_res <- as.data.frame(ego)
      if(nrow(ego_res)>0){
        write.xlsx(ego_res,paste(op,CT,"_BP.xlsx",sep=''),row.names = FALSE,sheetName = type,append = TRUE)
      } 
      
      #KEGG分析#
      ekk <- enrichKEGG(gene= id,organism  = 'hsa')	 #KEGG富集分析
      ekk_res <- as.data.frame(ekk)
      if(nrow(ekk_res)>0){
        write.xlsx(ekk_res,paste(op,CT,"_KEGG.xlsx",sep=''),row.names = FALSE,sheetName = type,append = TRUE)
      }
    }
  }
}
