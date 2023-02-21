#########
## 比较儿童和成年人年龄段(C4/C6 vs C1/C2/C3/C5)
#########
.libPaths(c("/home/yzj/R/x86_64-pc-linux-gnu-library/4.0","/home/zy/tools/R-4.0.0/library"))
.libPaths("/lustre/tianlab/yzj/software/R-4.0.0/lib64/R/library")
setwd('/lustre/tianlab/yzj')

library(Seurat)
library(ggplot2)
library(ggrepel)

###############
#  step1. DEG in children and adults
###############
library(future)
plan("multiprocess", workers = 80)

pbmc_chond <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype_Chond.Rdata')
Idents(pbmc_chond) <- pbmc_chond$type
pbmc_C <- subset(pbmc_chond,idents = 'Normal')
pbmc_C@meta.data$Phase <- 'Adults'
pbmc_C$Phase[pbmc_C$batch %in% c('C4','C6')] <- 'Children'
pbmc_C$Phase <- factor(pbmc_C$Phase,levels = c('Children','Adults'))
saveRDS(pbmc_C,'JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype_Control_Chond.Rdata')


## 1.1 filtering gene with cutoff
batch1='Children'
batch2='Adults'
resol='celltype'
MK.lst <- list()
celltypes <- unique(pbmc_C$celltype)

for (CT in celltypes) {
  print(paste('FindMarkers:',CT,'Start!'))
  cells <- colnames(pbmc_C)[pbmc_C$celltype == CT]
  subpbmc <- subset(pbmc_C,cells = cells)
  Idents(subpbmc) <- subpbmc$Phase
  cell_numbers <- as.numeric(table(subpbmc$Phase))
  if(length(cell_numbers) == 2 & all(cell_numbers>3)){
    sub.markers <- FindMarkers(subpbmc, ident.1 = batch1,ident.2 = batch2)
    MK.lst[[as.character(CT)]] <- sub.markers
  }else{
    print(paste(CT,':no correct cell numbers!'))
  }
}
saveRDS(MK.lst, file = paste('JingMA_NEW/res/Harmony/ALL/RDS/DEGs_inChond_in',batch1,batch2,'.RDS',sep=''))


### 1.2 GSEA: no filtering gene
pbmc_C <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype_Control_Chond.Rdata')
MK.lst <- list()
celltypes <- unique(pbmc_C$celltype)

for (CT in celltypes) {
  print(paste('FindMarkers:',CT,'Start!'))
  cells <- colnames(pbmc_C)[pbmc_C$celltype == CT]
  subpbmc <- subset(pbmc_C,cells = cells)
  Idents(subpbmc) <- subpbmc$Phase
  cell_numbers <- as.numeric(table(subpbmc$Phase))
  if(length(cell_numbers) == 2 & all(cell_numbers>3)){
    sub.markers <- FindMarkers(subpbmc, ident.1 = 'Children',ident.2 = 'Adults',
                               logfc.threshold = 0,min.pct = 0,min.cells.feature = 0)
    MK.lst[[as.character(CT)]] <- sub.markers
  }else{
    print(paste(CT,':no correct cell numbers!'))
  }
}
saveRDS(MK.lst, file ='JingMA_NEW/res/Harmony/ALL/RDS/DEGs_inChond_inChildrenAdults_GSEA.RDS')



################
## step2. pick DEG
################
library(xlsx)
library(Seurat)
library(ggplot2)

batch1='Children'
batch2='Adults'
resol='celltype'

pbmc_chond <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype_Chond.Rdata')
Idents(pbmc_chond) <- pbmc_chond$type
pbmc_C <- subset(pbmc_chond,idents = 'Normal')
pbmc_C@meta.data$Phase <- 'Adults'
pbmc_C$Phase[pbmc_C$batch %in% c('C4','C6')] <- 'Children'
pbmc_C$Phase <- factor(pbmc_C$Phase,levels = c('Children','Adults'))
VG <- VariableFeatures(pbmc_C)

MK.lst <- readRDS(paste('JingMA_NEW/res/Harmony/ALL/RDS/DEGs_inChond_in',batch1,batch2,'.RDS',sep=''))
wp=paste('JingMA_NEW/res/compControl/',batch1,'vs',batch2,sep='')

FC <- 1.5
adjP <- 0.05

CellTypes=c('CSC','TC','C1','C2')
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
  ## 相对于成人来说
  data$avg_logFC <- -data$avg_logFC
  
  # sigUP_gene <- data[data$gene %in% VG & data$p_val_adj  < adjP & data$avg_logFC > log(FC),]
  # sigDN_gene <- data[data$gene %in% VG & data$p_val_adj  < adjP & data$avg_logFC < -log(FC),]
  
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

batch1='Children'
batch2='Adults'
resol='celltype'
MK.lst <- readRDS(paste('JingMA_NEW/res/Harmony/ALL/RDS/DEGs_inChond_in',batch1,batch2,'.RDS',sep=''))

wp=paste('/home/yzj/JingMA_NEW/res/compControl/',batch1,'vs',batch2,sep='')

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

# ######
# ## step4. 小提琴图
# ######
# batch1='Children'
# batch2='Adult'
# resol='celltype'
# MK.lst <- readRDS(paste('JingMA_NEW/res/Harmony/ALL/RDS/DEGs_inChond_in',batch1,batch2,'.RDS',sep=''))
# 
# pbmc_chond <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype_Chond.Rdata')
# Idents(pbmc_chond) <- pbmc_chond$type
# pbmc_C <- subset(pbmc_chond,idents = 'Normal')
# pbmc_C@meta.data$Phase <- 'Adults'
# pbmc_C$Phase[pbmc_C$batch %in% c('C4','C6')] <- 'Children'
# pbmc_C$Phase <- factor(pbmc_C$Phase,levels = c('Children','Adults'))
# 
# wp=paste('/home/yzj/JingMA_NEW/res/compControl/',batch1,'vs',batch2,'/DEG/',sep='')
# op_batch=paste('/home/yzj/JingMA_NEW/res/compControl/',batch1,'vs',batch2,'/Vlnplot_batch/',sep='')
# op_phase=paste('/home/yzj/JingMA_NEW/res/compControl/',batch1,'vs',batch2,'/Vlnplot_phase/',sep='')
# 
# FC <- 1.5
# adjP <- 0.01
# 
# CellTypes=c('CSC','TC','C1','C2')
# for(i in 1:length(CellTypes)){
#   CT=CellTypes[i]
#   print(CT)
#   cells <- colnames(pbmc_C)[pbmc_C$celltype==CT]
#   subpbmc <- subset(pbmc_C,cells=cells)
#   Idents(subpbmc) <- subpbmc$batch
# 
#   newCT=paste(unlist(strsplit(CT,' ')),collapse = '_')
#   TAG=paste('FC',FC,'_adjP',adjP,'/',newCT,sep='')
#   ## mkdir
#   if(!file.exists(file.path(op_batch, TAG))){
#     dir.create(file.path(op_batch, TAG),recursive = TRUE)
#   }
#   if(!file.exists(file.path(op_phase, TAG))){
#     dir.create(file.path(op_phase, TAG),recursive = TRUE)
#   }
#   ##
#   sigUP_gene  <- read.csv(paste(wp,TAG,'/sigUP_gene.csv',sep=''))[,1]
#   sigDN_gene <- read.csv(paste(wp,TAG,'/sigDN_gene.csv',sep=''))[,1]
# 
#   NUM <- ceiling(sqrt(length(sigUP_gene)))
#   p_up <- VlnPlot(subpbmc,features = sigUP_gene, pt.size=0, group.by="batch", ncol=NUM)
#   ggsave(paste(op_batch,TAG,'/UPgene.pdf',sep=''), p_up, width=40 ,height=40)
# 
#   p_up <- VlnPlot(subpbmc,features = sigUP_gene, pt.size=0, group.by="Phase", ncol=NUM)
#   ggsave(paste(op_phase,TAG,'/UPgene.pdf',sep=''), p_up, width=40 ,height=40)
# 
#   NUM <- ceiling(sqrt(length(sigDN_gene)))
#   p_dn <- VlnPlot(subpbmc,features = sigDN_gene, pt.size=0, group.by="batch", ncol=NUM)
#   ggsave(paste(op_batch,TAG,'/DNgene.pdf',sep=''), p_dn, width=40 ,height=40)
# 
#   p_dn <- VlnPlot(subpbmc,features = sigDN_gene, pt.size=0, group.by="Phase", ncol=NUM)
#   ggsave(paste(op_phase,TAG,'/DNgene.pdf',sep=''), p_dn, width=40 ,height=40)
# }
# #####
# # step5. 热图
# #####
# library(tibble)
# library(dplyr)
# pbmc <- readRDS('/home/yzj/JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype.Rdata')
# Idents(pbmc) <- pbmc$type
# pbmc_C <- subset(pbmc,idents = 'Control')
# pbmc_C@meta.data$Phase <- 'Adult'
# pbmc_C$Phase[pbmc_C$batch %in% c('C4','C6')] <- 'Children'
# pbmc_C$Phase <- factor(pbmc_C$Phase,levels = c('Children','Adult'))
# 
# batch1='Children'
# batch2='Adult'
# resol='celltype'
# MK.lst <- readRDS(paste('JingMA_NEW/res/Harmony/ALL/RDS/DEGs_inChond_in',batch1,batch2,'.RDS',sep=''))
# 
# wp=paste('/home/yzj/JingMA_NEW/res/compControl/',batch1,'vs',batch2,sep='')
# op=paste('/home/yzj/JingMA_NEW/res/compControl/',batch1,'vs',batch2,'/Heatmap',sep='')
# ## mkdir
# if(!file.exists(file.path(op))){
#   dir.create(file.path(op),recursive = TRUE)
# }
# 
# CellTypes=c('CSC','TC',"Chondrocyte",'C1','C2')
# 
# for(i in 1:length(CellTypes)){
#   CT=CellTypes[i]
#   print(CT)
#   
#   if(CT=='Chondrocyte'){
#     cells <- colnames(pbmc_C)[pbmc_C$celltype %in% c('C1','C2')]
#   }else{
#     cells <- colnames(pbmc_C)[pbmc_C$celltype==CT]
#   }
#   subpbmc <- subset(pbmc_C,cells = cells)
#   Idents(subpbmc) <- subpbmc$Phase
#   
#   data <- MK.lst[[CT]]
#   data$gene <- rownames(data)
#   up <- data %>% as_tibble() %>% arrange(desc(avg_logFC)) %>% filter(p_val_adj < 0.01) %>%
#     head(n=50)
#   dn <- data %>% as_tibble() %>% arrange(desc(avg_logFC)) %>% filter(p_val_adj < 0.01) %>%
#     tail(n= 50)
#   genes <- c(up$gene,dn$gene)
#   pdf(paste(op,'/',CT,'.pdf',sep=''),height = 13,width = 5)
#   p <- DoHeatmap(subpbmc,features = genes,draw.lines = F,angle = 0)
#   print(p)
#   dev.off()
# }
# 
# ######
# ## step6.  GSEA for Children vs Adult
# ######
# library(enrichplot)
# library(clusterProfiler)
# library(org.Hs.eg.db)
# library(ggplot2)
# library(forcats)
# library(msigdbr) #提供MSigdb数据库基因集
# library(fgsea)
# library(Seurat)
# library(tibble)
# 
# run_GSEA <- function(wp,CT,MK.lst){
#   Dir <- c('GSEA')
#   for(d in Dir){
#     if(!file.exists(file.path(wp, d))){
#       dir.create(file.path(wp, d),recursive = TRUE)
#     }
#   }
#   
#   MK <- MK.lst[[as.character(CT)]]
#   geneList <- MK[, 'avg_logFC']
#   names(geneList) <- rownames(MK)
#   geneList <- na.omit(geneList)
#   geneList <- geneList[order(geneList, decreasing = T)]
#   
#   ## GO BP
#   ego <- gseGO(geneList = geneList, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", pvalueCutoff = 0.05)
#   ego <- clusterProfiler::simplify(ego,cutoff=0.8,by="p.adjust",select_fun=min,measure="Wang")
#   ego_res <- as.data.frame(ego)
#   write.csv(ego_res,paste(wp,'/GSEA/',CT,'_BP.csv',sep=''),quote=F,row.names = F)
#   
#   ## KEGG
#   symbol2id=mapIds(org.Hs.eg.db,names(geneList),"ENTREZID",'SYMBOL')
#   geneList_kegg <- geneList
#   names(geneList_kegg) <- symbol2id
#   geneList_kegg=geneList_kegg[which(names(geneList_kegg)!='')] #提取出非NA的ENTREZID
#   
#   ekk <- gseKEGG(geneList  = geneList_kegg,keyType  = 'kegg',organism = 'hsa',pvalueCutoff = 0.05,pAdjustMethod= "BH")
#   ekk_res <- as.data.frame(ekk)
#   write.csv(ekk_res,paste(wp,'/GSEA/',CT,'_KEGG.csv',sep=''),quote=F,row.names = F)
#   
#   enrich_res <- list()
#   enrich_res[['BP']] <- ego
#   enrich_res[['KEGG']] <- ekk
#   saveRDS(enrich_res,paste(wp,'/GSEA/',CT,'_enrich_res.RDS',sep=''))
#   return(enrich_res)
# }
# 
# get_gseaBarplot <- function(wp,TAG,topN){
#   enrich_res <- readRDS(paste(wp,'/GSEA/',TAG,'_enrich_res.RDS',sep=''))
#   ego <- enrich_res[['BP']]
#   ekk <- enrich_res[['KEGG']]
#   
#   NUM=nrow(ego)
#   
#   if(length(NUM) >0 ){
#     if(0 < NUM & NUM < 20){
#       height=8
#     } else if (NUM >= 20){
#       height=16
#     }
#     if(NUM >=1){
#       up <- ego %>% as_tibble() %>% arrange(desc(NES)) %>% filter(p.adjust < 0.05) %>%
#         head(n=topN/2)
#       dn <- ego %>% as_tibble() %>% arrange(desc(NES)) %>% filter(p.adjust < 0.05) %>%
#         tail(n= topN/2)
#       BP <- ggplot(rbind(up,dn) , aes(reorder(Description, NES), NES)) +
#         geom_col(aes(fill= NES < 0)) +
#         theme_classic()+
#         coord_flip() +
#         theme(text = element_text(size = 15))+
#         labs(x="Pathway", y="Normalized Enrichment Score",
#              title="GO BP NES from GSEA")
#       pdf(paste(wp,'/GSEA/',TAG,'_BP.pdf',sep=''),width = 14,height = height)
#       print(BP)
#       dev.off()
#       
#     }
#   }else{
#     print('no GO pdf!')
#   }
#   
#   
#   NUM=nrow(ekk)
#   if(length(NUM) >0 ){
#     if( 0 < NUM & NUM < 20){
#       height=8
#     } else if(20 <= NUM & NUM < 30){
#       height=10
#     }else if(NUM >= 30){
#       height=16
#     }
#     if(NUM >=1){
#       up <- ekk %>% as_tibble() %>% arrange(desc(NES)) %>% filter(p.adjust < 0.05) %>%
#         head(n=topN/2)
#       dn <- ekk %>% as_tibble() %>% arrange(desc(NES)) %>% filter(p.adjust < 0.05) %>%
#         tail(n= topN/2)
#       KEGG <- ggplot(rbind(up,dn) , aes(reorder(Description, NES), NES)) +
#         geom_col(aes(fill= NES < 0)) +
#         theme_classic()+
#         coord_flip() +
#         theme(text = element_text(size = 15))+
#         labs(x="Pathway", y="Normalized Enrichment Score",
#              title="KEGG NES from GSEA")
#       pdf(paste(wp,'/GSEA/',TAG,'_KEGG.pdf',sep=''),width = 14,height = height)
#       print(KEGG)
#       dev.off()
#     }
#   }else{
#     print('no KEGG pdf!')
#   }
# }
# 
# resol='celltype'
# batch1 <- 'Children'
# batch2 <- 'Adult'
# MK.lst <- readRDS(paste('JingMA_NEW/res/Harmony/ALL/RDS/DEGs_inChond_in',batch1,batch2,'.RDS',sep=''))
# wp=paste('/home/yzj/JingMA_NEW/res/compControl/',batch1,'vs',batch2,sep='')
# 
# CellTypes=c('CSC','TC',"Chondrocyte",'C1','C2')
# for(CT in CellTypes){
#   print(CT)
#   res <- run_GSEA(wp,CT,MK.lst)
#   get_gseaBarplot(wp,CT,50)
# }
# 
# #####
# # 7. EnrichR
# #####
# .libPaths("/home/yzj/R/x86_64-pc-linux-gnu-library/4.0")
# library(enrichR)
# library(RColorBrewer)
# library("VennDiagram")
# color <- brewer.pal(9,"Set1");
# mycol <- colorRampPalette(color)(9)
# dbs <- listEnrichrDbs()
# 
# height <- function(df){
#   if(nrow(df) >= 1 & nrow(df) <= 20 ){
#     height=7
#   } else if(nrow(df) > 20 & nrow(df) <= 60 ) {
#     height=15
#   } else if(nrow(df) > 60 & nrow(df) <= 200 ) {
#     height=18
#   } else if (nrow(df) > 200 ) {
#     height=140
#   }
#   return(height)
# }
# 
# batch1='Children'
# batch2='Adult'
# resol='celltype'
# MK.lst <- readRDS(paste('JingMA_NEW/res/Harmony/ALL/RDS/DEGs_inChond_in',batch1,batch2,'.RDS',sep=''))
# 
# wp=paste('/home/yzj/JingMA_NEW/res/compControl/',batch1,'vs',batch2,'/DEG/',sep='')
# op=paste('/home/yzj/JingMA_NEW/res/compControl/',batch1,'vs',batch2,'/EnrichR/',sep='')
# 
# FC <- 2
# adjP <- 0.01
# CUT=paste('FC',FC,'_adjP',adjP,sep='')
# 
# 
# CellTypes=c('CSC','TC',"Chondrocyte",'C1','C2')
# for(i in 1:length(CellTypes)){
#   CT <- CellTypes[i]
#   print(CT)
# 
#   newCT=paste(unlist(strsplit(CT,' ')),collapse = '_')
#   TAG=paste('FC',FC,'_adjP',adjP,'/',newCT,sep='')
# 
#   if(!file.exists(file.path(op, TAG))){
#     dir.create(file.path(op, TAG),recursive = TRUE)
#   }
# 
#   TYPE <- c('UP','DN')
#   for(j in 1:length(TYPE)){
#     type=TYPE[j]
#     sigG_lst <- read.csv(paste(wp,TAG,'/sig',type,'_gene.csv',sep=''))[,1]
# 
#     ## GO BP
#     enrichr_res=enrichr(sigG_lst,databases=c('GO_Biological_Process_2018'))
#     tmp=enrichr_res[[1]][,c(1,4)]
#     print(dim(tmp))
#     if(nrow(tmp) > 0){
#       tmp=tmp[which(tmp$Adjusted.P.value < 0.05),];rownames(tmp)=tmp[,1]
#       tmp=tmp[order(tmp[,2],decreasing=F),]
#       print(dim(tmp))
# 
#       if(nrow(tmp) > 0){
#         write.csv(tmp,paste0(op,TAG,"/BP_",type,".csv"),quote=F, row.names=F )
#         if(nrow(tmp) >50){
#           tmp <-tmp[1:50,]
#         }else{
#           tmp <- tmp
#         }
#         h=height(tmp)
#         mycolor<-rep(mycol[2],nrow(tmp))
#         pdf(paste0(op,TAG,"/BP_",type,".pdf"),width=13,height=h)
#         bb<-barplot(c(-log(tmp[,2],10)),horiz=TRUE,axisnames=F,main=TAG,cex.axis=0.9,col=mycolor,border=NA,xlim=c(0,25),space=0.3,xlab="-log10(qvalue)",ylab='GO BP')
#         text(0.001,bb,rownames(tmp),cex=0.8,xpd=T,srt=0,pos=4)
#         dev.off()
#       }else{
#         print('no sigBP!')
#       }
#     }else{
#       print('no BP!')
#     }
# 
#     ## KEGG
#     enrichr_res=enrichr(sigG_lst,databases=c('KEGG_2019_Human'))
#     tmp=enrichr_res[[1]][,c(1,4)]
#     print(dim(tmp))
#     if(nrow(tmp) > 0){
#       tmp=tmp[which(tmp$Adjusted.P.value < 0.05),];rownames(tmp)=tmp[,1]
#       tmp=tmp[order(tmp[,2],decreasing=F),]
#       print(dim(tmp))
# 
#       if(nrow(tmp) > 0){
#         write.csv(tmp,paste0(op,TAG,"/KEGG_",type,".csv"),quote=F, row.names=F )
#         if(nrow(tmp) >50){
#           tmp <-tmp[1:50,]
#         }else{
#           tmp <- tmp
#         }
#         h=height(tmp)
#         mycolor<-rep(mycol[2],nrow(tmp))
#         pdf(paste0(op,TAG,"/KEGG_",type,".pdf"),width=13,height=h)
#         bb<-barplot(c(-log(tmp[,2],10)),horiz=TRUE,axisnames=F,main=TAG,cex.axis=0.9,col=mycolor,border=NA,xlim=c(0,25),space=0.3,xlab="-log10(qvalue)",ylab='GO BP')
#         text(0.001,bb,rownames(tmp),cex=0.8,xpd=T,srt=0,pos=4)
#         dev.off()
#       }else{
#         print('no sigKEGG!')
#       }
#     }else{
#       print('no KEGG!')
#     }
# 
#   }
# }
# 
# 
# 
# ######
# ## step8. 火山图,并加上感兴趣基因的名字
# ######
# batch1='Children'
# batch2='Adult'
# resol='celltype'
# MK.lst <- readRDS(paste('JingMA_NEW/res/Harmony/ALL/RDS/DEGs_inChond_in',batch1,batch2,'.RDS',sep=''))
# 
# wp=paste('/home/yzj/JingMA_NEW/res/',batch1,'vs',batch2,sep='')
# op=paste('/home/yzj/JingMA_NEW/res/',batch1,'vs',batch2,'/Volcano',sep='')
# 
# FC <- 2
# adjP <- 0.01
# 
# CellTypes=c('CSPC')
# for(i in 1:length(CellTypes)){
#   CT=CellTypes[i]
#   print(CT)
#   CUT=paste('FC',FC,'_adjP',adjP,sep='')
#   TAG=paste('FC',FC,'_adjP',adjP,'/',CT,sep='')
#   ## mkdir
#   if(!file.exists(file.path(op, CUT))){
#     dir.create(file.path(op, CUT),recursive = TRUE)
#   }
#   data <- MK.lst[[CT]]
#   data$gene <- rownames(data)
#   data$significant <- as.factor(ifelse(data$p_val_adj < adjP & abs(data$avg_logFC) >= log(FC),
#                                        ifelse(data$avg_logFC >= log(FC),'UP','DOWN'),'NOT'))
#   print(table(data$significant))
# 
#   gene<-data$gene
#   gene_logFC<-data$avg_logFC
#   gene_logAdjP<--log10(data$p_val_adj)
#   boolean_isInf<-is.infinite(gene_logFC)
#   filter_gene<-as.vector(gene[!boolean_isInf])
#   filter_gene_logFC<-gene_logFC[!boolean_isInf]
#   filter_gene_logAdjP<-gene_logAdjP[!boolean_isInf]
#   significant<-data$significant[!boolean_isInf]
#   filter_data<-data.frame(filter_gene,filter_gene_logFC,filter_gene_logAdjP,significant)
#   boolean_tag=(filter_gene_logFC> log(FC) | filter_gene_logFC< -log(FC)) & filter_gene_logAdjP > -log10(adjP)
#   filter_gene[which(boolean_tag==F)]=""
# 
#   ## CSPC
#   tag_lst=c("SPARC","ELN","TIMP4","MMP3","FST","IL6","IL8","CDKN1A","MMP1")
#   index_tag=which(!filter_gene %in% tag_lst)
#   filter_gene[index_tag]=""
# 
# 
#   pdf(paste(op,'/',CUT,'/',CT,'.pdf',sep=""))
#   figure<-ggplot(filter_data,aes(x=filter_gene_logFC,y=filter_gene_logAdjP))
#   xstt <- floor(min(filter_data$filter_gene_logFC))
#   xend <- ceiling(max(filter_data$filter_gene_logFC))
#   x <- max(abs(xstt),abs(xend))
#   xstt <- -x
#   xend <- x
# 
#   ystt <- floor(min(filter_data$filter_gene_logAdjP))
#   yend <- ceiling(max(filter_data$filter_gene_logAdjP))
# 
#   figure1 <- figure+geom_point(aes(color=significant))+xlim(xstt,xend) + ylim(ystt,yend)+
#     scale_color_manual(values = c("#377EB8","#999999","#E41A1C"))+
#     labs(title=paste(batch1,'vs',batch2,'in',CT),x="logFC",y="-log10(adjP)")+
#     theme(plot.title = element_text(hjust = 0.5))+
#     geom_hline(yintercept = -log10(adjP),linetype=3)+
#     geom_vline(xintercept=c(-log(FC),log(FC)),linetype=3)
#   print(figure1+geom_text_repel(label=filter_gene))
#   dev.off()
# }

