## 1. Find MK of each cell type in  all/Chond cells;
library(Seurat)
library(ggplot2)
library(dplyr)
library(future)
plan("multiprocess", workers = 60)

###################
## 1.1 Find MK of each cell type in  all cells;
###################
setwd('/lustre/tianlab/yzj')

pbmc <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype.Rdata')
Idents(pbmc) <- pbmc$celltype.abbr
celltypes <- levels(pbmc$celltype.abbr)

MK.lst <- list()
for (CT in celltypes) {
  print(paste('FindMarkers:',CT,'Start!'))
  markers <- FindMarkers(pbmc, ident.1 = CT,only.pos = TRUE)
  MK.lst[[as.character(CT)]] <- markers
}
saveRDS(MK.lst, file = 'JingMA_NEW/res/Harmony/ALL/RDS/Markers_celltype.RDS')

### 转成矩阵
MK.lst <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/Markers_celltype.RDS')
MK_df <- c()
for(i in 1:length(MK.lst)){
  ct=names(MK.lst)[i]
  MK <- MK.lst[[ct]]
  MK <-tibble::rownames_to_column(MK,'gene')
  MK$CT <-ct
  #MK <- MK[MK$avg_logFC>log(1.5)&MK$p_val_adj < 0.05,]
  MK_df <- rbind(MK_df,MK)
}
saveRDS(MK_df,'JingMA_NEW/res/Harmony/ALL/RDS/Markers_celltype_mtx.RDS')


###################
## 1.2 Find MK of each cell type in Chond cells
###################
setwd('/lustre/tianlab/yzj')

pbmc_chond <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype_Chond.Rdata')
Idents(pbmc_chond) <- pbmc_chond$celltype
celltypes <- levels(pbmc_chond$celltype)

MK.lst <- list()
for (CT in celltypes) {
  print(paste('FindMarkers:',CT,'Start!'))
  markers <- FindMarkers(pbmc_chond, ident.1 = CT,only.pos = TRUE)
  MK.lst[[as.character(CT)]] <- markers
}
saveRDS(MK.lst, file = 'JingMA_NEW/res/Harmony/ALL/RDS/Markers_celltype_Chond.RDS')

### 转成矩阵
MK.lst <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/Markers_celltype_Chond.RDS')
MK_df <- c()
for(i in 1:length(MK.lst)){
  ct=names(MK.lst)[i]
  MK <- MK.lst[[ct]]
  MK <-tibble::rownames_to_column(MK,'gene')
  MK$CT <-ct
  MK_df <- rbind(MK_df,MK)
}
saveRDS(MK_df,'JingMA_NEW/res/Harmony/ALL/RDS/Markers_celltype_Chond_mtx.RDS')


###################
## 1.3 Find MK of each cell type in Stroma cells
###################
setwd('/lustre/tianlab/yzj')

pbmc_stroma <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype_Stroma.Rdata')
Idents(pbmc_stroma) <- pbmc_stroma$celltype
celltypes <- levels(pbmc_stroma$celltype)

MK.lst <- list()
for (CT in celltypes) {
  print(paste('FindMarkers:',CT,'Start!'))
  markers <- FindMarkers(pbmc_stroma, ident.1 = CT,only.pos = TRUE)
  MK.lst[[as.character(CT)]] <- markers
}
saveRDS(MK.lst, file = 'JingMA_NEW/res/Harmony/ALL/RDS/Markers_celltype_Stroma.RDS')

### 转成矩阵
MK.lst <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/Markers_celltype_Stroma.RDS')
MK_df <- c()
for(i in 1:length(MK.lst)){
  ct=names(MK.lst)[i]
  MK <- MK.lst[[ct]]
  MK <-tibble::rownames_to_column(MK,'gene')
  MK$CT <-ct
  MK_df <- rbind(MK_df,MK)
}
saveRDS(MK_df,'JingMA_NEW/res/Harmony/ALL/RDS/Markers_celltype_Stroma_mtx.RDS')


###################
## 1.4 Find MK of each cell subtype in  all cells;
###################
setwd('/lustre/tianlab/yzj')

pbmc <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype.Rdata')
Idents(pbmc) <- pbmc$celltype
celltypes <- levels(pbmc$celltype)

MK.lst <- list()
for (CT in celltypes) {
  print(paste('FindMarkers:',CT,'Start!'))
  markers <- FindMarkers(pbmc, ident.1 = CT,only.pos = TRUE)
  MK.lst[[as.character(CT)]] <- markers
}
saveRDS(MK.lst, file = 'JingMA_NEW/res/Harmony/ALL/RDS/Markers_subcelltype.RDS')

### 转成矩阵
MK.lst <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/Markers_subcelltype.RDS')
MK_df <- c()
for(i in 1:length(MK.lst)){
  ct=names(MK.lst)[i]
  MK <- MK.lst[[ct]]
  MK <-tibble::rownames_to_column(MK,'gene')
  MK$CT <-ct
  MK_df <- rbind(MK_df,MK)
}
saveRDS(MK_df,'JingMA_NEW/res/Harmony/ALL/RDS/Markers_subcelltype_mtx.RDS')


###################
## 1.5 Find MK of each cell subtype in main cells
###################
setwd('/lustre/tianlab/yzj')

pbmc<- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype.Rdata')
pbmc_main <- subset(pbmc,cells=colnames(pbmc)[pbmc$celltype.abbr %in% c("CSC","C",'SSC','SC')])
pbmc_main$celltype.abbr <- factor(pbmc_main$celltype.abbr,levels = c("CSC","C",'SSC','SC'))
pbmc_main$celltype <- factor(pbmc_main$celltype,levels = c("CSC","C0","C1","C2",'SSC','SC1','SC2'))
celltypes <- levels(pbmc_main$celltype)
Idents(pbmc_main) <- pbmc_main$celltype

MK.lst <- list()
for (CT in celltypes) {
  print(paste('FindMarkers:',CT,'Start!'))
  markers <- FindMarkers(pbmc_main, ident.1 = CT,only.pos = TRUE)
  MK.lst[[as.character(CT)]] <- markers
}
saveRDS(MK.lst, file = 'JingMA_NEW/res/Harmony/ALL/RDS/Markers_celltype_Main.RDS')

### 转成矩阵
MK.lst <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/Markers_celltype_Main.RDS')
MK_df <- c()
for(i in 1:length(MK.lst)){
  ct=names(MK.lst)[i]
  MK <- MK.lst[[ct]]
  MK <-tibble::rownames_to_column(MK,'gene')
  MK$CT <-ct
  MK_df <- rbind(MK_df,MK)
}
saveRDS(MK_df,'JingMA_NEW/res/Harmony/ALL/RDS/Markers_celltype_Main_mtx.RDS')


MK_Chond_df <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/Markers_celltype_Chond_mtx.RDS')
MK_Chond_df$CT[MK_Chond_df$CT=='TC'] <- 'C0'
MK_Stroma_df <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/Markers_celltype_Stroma_mtx.RDS')
library(xlsx)
write.xlsx(MK_Chond_df,'JingMA_NEW/res/Harmony/ALL/RDS/Markers_celltype_Chond_Stroma.xlsx',
           append = TRUE,sheetName = 'Chondral',row.names = FALSE)
write.xlsx(MK_Stroma_df,'JingMA_NEW/res/Harmony/ALL/RDS/Markers_celltype_Chond_Stroma.xlsx',append=TRUE,
           row.names = FALSE,sheetName = 'periChondral')


##############################
#### 补充材料
##############################
#### tableS2: 每种细胞类型的DEG
library(xlsx)
MK.lst <- readRDS( 'JingMA_NEW/res/Harmony/ALL/RDS/Markers_celltype.RDS')
CT=  names(MK.lst)

for(ct in names(MK.lst)){
  DEG <- MK.lst[[ct]]
  DEG <- DEG[DEG$avg_logFC > log(1.5) & DEG$p_val_adj < 0.05,]
  write.xlsx(DEG,'JingMA_NEW/res/Harmony/ALL/FIG/TableS2_MK.xlsx',append = TRUE,sheetName = ct)
}
