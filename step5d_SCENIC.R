library(Seurat)
library(SCENIC)
library(foreach)
library(pheatmap)
library(AUCell)
library(ggplot2)
library(limma)

############
## 根据Binary情况，筛选每种细胞类型的regulon
# ############
pbmc <- readRDS('/local/yzj/JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype_Chond.Rdata')
Meta_df <- pbmc@meta.data
Meta_df <-tibble::rownames_to_column(Meta_df,var = 'cell')

regulonAUC <-  readRDS("~/JingMA_NEW/res/SCENIC_main/int/4.2_binaryRegulonActivity_nonDupl.Rds")

cellInfo <- data.frame(Meta_df[Meta_df$cell %in% colnames(regulonAUC),c('cell','celltype')])
rownames(cellInfo) <- cellInfo$cell
colnames(cellInfo)[2] <- c('celltype')

regulonAUC <- regulonAUC[,colnames(regulonAUC) %in% colnames(pbmc)]
print(all(colnames(pbmc)==colnames(regulonAUC)))

RA_byCellType <- sapply(split(rownames(cellInfo),cellInfo$celltype),
                        function(cells) rowMeans(regulonAUC[,cells]))

RA_byCellType.lst <- list()
for(i in 1:ncol(RA_byCellType)){
  ct <- colnames(RA_byCellType)[i]
  RA_byCellType.lst[[ct]] <- rownames(RA_byCellType)[RA_byCellType[,i]>=0.1] # >10%细胞中活跃即可
}
saveRDS(RA_byCellType.lst,'/local/yzj/JingMA_NEW/res/SCENIC_main/FIG/RA_byCellType.lst.RDS')



############
## regulon在每种细胞类型中的activity
# ############
regulonAUC <- readRDS(file='/local/yzj/JingMA_NEW/res/SCENIC_main/int/3.4_regulonAUC.Rds')
AUC_df <- getAUC(regulonAUC)

cellInfo <- data.frame(Meta_df[Meta_df$cell %in% colnames(regulonAUC),c('cell','celltype')])
rownames(cellInfo) <- cellInfo$cell
colnames(cellInfo)[2] <- c('celltype')

AUC_df <- AUC_df[,cellInfo$cell]
print(all(colnames(AUC_df)==colnames(cellInfo$cell)))

RA_byCellType <- sapply(split(rownames(cellInfo),cellInfo$celltype),
                        function(cells) rowMeans(AUC_df[,cells]))
saveRDS(RA_byCellType,'JingMA_NEW/res/SCENIC_main/FIG/RA_df_noScale_allRegulon.RDS')

############
## regulon在正常和疾病情况下每种细胞类型中的activity
#############
regulonAUC <- readRDS(file='/local/yzj/JingMA_NEW/res/SCENIC_main/int/3.4_regulonAUC.Rds')
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]

regulonAUC <- getAUC(regulonAUC)
regulonAUC <- t(scale(t(regulonAUC), center = T, scale=T))
dim(regulonAUC)

cellInfo <- read.table('JingMA_NEW/res/CellPhoneDB_merge_subcelltype/in/meta.txt',header = T,sep='\t',stringsAsFactors = F)
cellInfo <- cellInfo[cellInfo$Cell %in% colnames(regulonAUC),]
print(all(cellInfo$Cell==colnames(regulonAUC)))

pbmc.AUC=CreateSeuratObject(counts = regulonAUC, min.cells = 0, min.features = 0)

pbmc.AUC@meta.data[1:2,]
pbmc.AUC$type <- 'Normal'
pbmc.AUC$type[pbmc.AUC$orig.ident %in% c('M1','M2','M3')] <- 'Microtia'
pbmc.AUC$celltype <- cellInfo$CellType
Idents(pbmc.AUC) <- pbmc.AUC$celltype

MK.lst <- list()
celltypes <- unique(pbmc.AUC$celltype)

for (CT in celltypes) {
  print(paste('FindMarkers:',CT,'Start!'))
  subpbmc <- subset(pbmc.AUC,idents=CT)
  Idents(subpbmc) <- subpbmc$type
  cell_numbers <- as.numeric(table(subpbmc$type))
  if(length(cell_numbers) == 2 & all(cell_numbers>3)){
    sub.markers <- FindMarkers(subpbmc, ident.1 = 'Microtia',ident.2 = 'Normal',
                               min.pct = 0,logfc.threshold =0)
    MK.lst[[as.character(CT)]] <- sub.markers
  }else{
    print(paste(CT,':no correct cell numbers!'))
  }
}
saveRDS(MK.lst, file = '/local/yzj/JingMA_NEW/res/SCENIC_main/FIG/DARegulon.RDS')




############
## 异常表达的TF
# ############
# pbmc <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype.Rdata')
# Idents(pbmc) <- pbmc$celltype.abbr
# regulonGeneSet <- readRDS('/home/yzj/JingMA_NEW/res/SCENIC_All_2/int/2.6_regulons_asGeneSet.Rds')
# mkLst <- readRDS('/home/yzj/JingMA_NEW/res/Harmony/ALL/RDS/Markers_MicrotiaConChildren_celltype_0_forGSEA.RDS')
# 
# RA_Control <- readRDS('/home/yzj/JingMA_NEW/res/SCENIC_All_2/FIG/heatmap.RDS')
# RA_Control_df <- as.data.frame(RA_Control)
# 
# CSC_TF <- rownames(RA_Control_df)[RA_Control_df$CSC>0]
# Trans_TF <- rownames(RA_Control_df)[RA_Control_df$TransC>0]
# Chond1_TF <- rownames(RA_Control_df)[RA_Control_df$Chond1>0]
# Chond2_TF <- rownames(RA_Control_df)[RA_Control_df$Chond2>0]
# 
# CSC_TF <- apply(as.matrix(CSC_TF),1,function(x) unlist(strsplit(x," "))[1])
# CSC_TF <- apply(as.matrix(CSC_TF),1,function(x) unlist(strsplit(x,"_extended"))[1])
# CSC_TF <- unique(CSC_TF)
# 
# Trans_TF <- apply(as.matrix(Trans_TF),1,function(x) unlist(strsplit(x," "))[1])
# Trans_TF <- apply(as.matrix(Trans_TF),1,function(x) unlist(strsplit(x,"_extended"))[1])
# Trans_TF <- unique(Trans_TF)
# 
# Chond1_TF <- apply(as.matrix(Chond1_TF),1,function(x) unlist(strsplit(x," "))[1])
# Chond1_TF <- apply(as.matrix(Chond1_TF),1,function(x) unlist(strsplit(x,"_extended"))[1])
# Chond1_TF <- unique(Chond1_TF)
# 
# Chond2_TF <- apply(as.matrix(Chond2_TF),1,function(x) unlist(strsplit(x," "))[1])
# Chond2_TF <- apply(as.matrix(Chond2_TF),1,function(x) unlist(strsplit(x,"_extended"))[1])
# Chond2_TF <- unique(Chond2_TF)
# 
# 
# Idents(pbmc) <- pbmc$batch
# subpbmc <- subset(pbmc,idents = c('C4','C6','M1','M2','M3'))
# Idents(subpbmc) <- subpbmc$celltype.abbr
# 
# #### CSC中差异的TF 
# MK <- mkLst$`Chondral stem cell`
# pbmc_CSC <- subset(subpbmc,idents=c("CSC"))
# Idents(pbmc_CSC) <- pbmc_CSC$type
# sigCSC_TF_df <- (MK[rownames(MK)%in% CSC_TF & abs(MK$avg_logFC) > log(1.5) & MK$p_val_adj < 0.05,])
# sigCSC_TF <- rownames(sigCSC_TF_df)
# 
# tmp <- AverageExpression(pbmc_CSC,features =sigCSC_TF,slot = 'data')
# EXP <- as.matrix(tmp$RNA)
# library(pheatmap)
# pheatmap(EXP,scale = 'row',cluster_rows = T,cluster_cols = F,show_rownames = T,treeheight_row = 0)
# VlnPlot(pbmc_CSC,features = sigCSC_TF,pt.size = 0)
# 
# #### Trans中差异的TF 
# MK <- mkLst$`Transitional chondrocyte`
# pbmc_Trans <- subset(subpbmc,idents=c("TC"))
# Idents(pbmc_Trans) <- pbmc_Trans$type
# sigTrans_TF_df <- (MK[rownames(MK)%in% Trans_TF & abs(MK$avg_logFC) > log(1.5) & MK$p_val_adj < 0.05,])
# sigTrans_TF <- rownames(sigTrans_TF_df)
# 
# tmp <- AverageExpression(pbmc_Trans,features =sigTrans_TF,slot = 'data')
# EXP <- as.matrix(tmp$RNA)
# library(pheatmap)
# pheatmap(EXP,scale = 'row',cluster_rows = T,cluster_cols = F,show_rownames = T,treeheight_row = 0)
# VlnPlot(pbmc_Trans,features = sigTrans_TF,pt.size = 0)
# 
# #### Chond1中差异的TF 
# MK <- mkLst$Chondrocyte1
# pbmc_Chond1 <- subset(subpbmc,idents=c("Chondrocyte1"))
# Idents(pbmc_Chond1) <- pbmc_Chond1$type
# sigChond1_TF_df <- (MK[rownames(MK)%in% Chond1_TF & abs(MK$avg_logFC) > log(1.5) & MK$p_val_adj < 0.05,])
# sigChond1_TF <- rownames(sigChond1_TF_df)
# 
# tmp <- AverageExpression(pbmc_Chond1,features =sigChond1_TF,slot = 'data')
# EXP <- as.matrix(tmp$RNA)
# library(pheatmap)
# pheatmap(EXP,scale = 'row',cluster_rows = T,cluster_cols = F,show_rownames = T,treeheight_row = 0)
# VlnPlot(pbmc_Chond1,features = sigChond1_TF,pt.size = 0)
# 
# #### Chond2中差异的TF 
# MK <- mkLst$Chondrocyte2
# pbmc_Chond2 <- subset(subpbmc,idents=c("Chondrocyte2"))
# Idents(pbmc_Chond2) <- pbmc_Chond2$type
# sigChond2_TF_df <- (MK[rownames(MK)%in% Chond2_TF & abs(MK$avg_logFC) > log(1.5) & MK$p_val_adj < 0.05,])
# sigChond2_TF <- rownames(sigChond2_TF_df)
# 
# tmp <- AverageExpression(pbmc_Chond2,features =sigChond2_TF,slot = 'data')
# EXP <- as.matrix(tmp$RNA)
# library(pheatmap)
# pheatmap(EXP,scale = 'row',cluster_rows = T,cluster_cols = F,show_rownames = T,treeheight_row = 0)
# VlnPlot(pbmc_Chond2,features = sigChond2_TF,pt.size = 0,ncol = 2)
# 
# 
# ### pick SOX8画图
# mkLst <- readRDS('/home/yzj/JingMA_NEW/res/Harmony/ALL/RDS/Markers_MicrotiaConChildren_celltype_0_forGSEA.RDS')
# geneLst <- c('SOX8', 'COL11A1', 'COL11A2', 'COL9A1', 'COL9A3', 'FGF1', 'FGFR2', 'FGFR3', 'SCRG1', 'SPARC', 'CHAD', 'CTGF', 'FOS', 'A2M', 'ITGA2', 'SOX5', 'S100B', 'VIT', 'BMP4', 'COL2A1', 'SAMD11')
# Idents(pbmc) <- pbmc$celltype.abbr
# pbmc_CSC <- subset(pbmc,idents='CSC')
# Idents(pbmc_CSC) <- pbmc_CSC$type
# 
# MK <- mkLst$`Chondral stem cell`
# sigGenes <- rownames(MK)[rownames(MK)%in% geneLst & MK$avg_logFC < -log(1.5) & MK$p_val_adj < 0.05]
# tmp <- AverageExpression(pbmc_CSC,features =sigGenes,slot = 'scale.data')
# EXP <- as.matrix(tmp$RNA)
# library(pheatmap)
# p <- pheatmap(EXP,color = colorRampPalette(c("navy", "white", "firebrick3"))(100),cluster_rows = T,cluster_cols = F,show_rownames = T,treeheight_row = 0)
# save_pheatmap_pdf(p,'/home/yzj/JingMA_NEW/res/compMicrotia/SOX8_heatmap.pdf',width = 3,height = 4)
# 
# ### pick ATF3画图
# geneLst <- c('ATF3','BAG3','DEDD2','SKIL','ATF3','DEDD2','SKIL','ATF3','CTH','DDIT3','DEDD2','HSPA1A','NR4A2','SKIL','DDIT3','HSPA1A','PPP1R15A')
# Idents(pbmc) <- pbmc$CellType
# 
# pbmc_CSC <- subset(pbmc,idents='CSC')
# Idents(pbmc_CSC) <- pbmc_CSC$type
# 
# MK <- mkLst$`Chondral stem cell`
# sigGenes <- rownames(MK)[rownames(MK)%in% geneLst & MK$avg_logFC > log(1.5) & MK$p_val_adj < 0.05]
# tmp <- AverageExpression(pbmc_CSC,features =sigGenes,slot = 'scale.data')
# EXP <- as.matrix(tmp$RNA)
# library(pheatmap)
# EXP <- EXP[c(3,2,4,1,5),]
# p <- pheatmap(EXP,color = colorRampPalette(c("navy", "white", "firebrick3"))(100),cluster_rows = F,cluster_cols = F,show_rownames = T,treeheight_row = 0)
# save_pheatmap_pdf(p,'/home/yzj/JingMA_NEW/res/compMicrotia/ATF3_heatmap.pdf',width = 3,height = 3)
