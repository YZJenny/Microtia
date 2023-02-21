library(Seurat)
library(SCENIC)
library(foreach)
library(pheatmap)
library(AUCell)
library(ggplot2)
library(limma)

## SCENIC结果解读
pbmc <- readRDS('/local/yzj/JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype.Rdata')
pbmc_main <- subset(pbmc,cells=colnames(pbmc)[pbmc$celltype %in% c("CSPC","Chond",'SSPC','SC')])
pbmc_main$celltype <- factor(pbmc_main$celltype,levels = c("CSPC","Chond",'SSPC','SC'))
pbmc_main$subcelltype <- factor(pbmc_main$subcelltype,levels = c("CSPC","EC","IC","LC",'SSPC','SC1','SC2'))

regulonAUC <- readRDS(file='/local/yzj/JingMA_NEW/res/SCENIC_main/int/3.4_regulonAUC.Rds')
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]

Meta_df <- pbmc_main@meta.data
Meta_df <-tibble::rownames_to_column(Meta_df,var = 'cell')

############
## 2. TF热图
############

########################
### 只看大类"CSC","C",'SSC','SC'
########################
cellInfo <- data.frame(Meta_df[Meta_df$cell %in% colnames(regulonAUC),c('cell','celltype')])
rownames(cellInfo) <- cellInfo$cell
colnames(cellInfo)[2] <- c('celltype')

RA_byCellType <- sapply(split(rownames(cellInfo),cellInfo$celltype),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
RA_byCellType_Scale <- t(scale(t(RA_byCellType), center = T, scale=T))

RA_df <- RA_byCellType_Scale

p <- pheatmap::pheatmap(RA_df, #fontsize_row=3,
                         color=colorRampPalette(c("navy","white","firebrick3"))(100), breaks=seq(-1.5, 1.5, length.out = 100),
                         treeheight_row=10, treeheight_col=10, border_color=NA,
                         cluster_cols = FALSE,cluster_rows = TRUE)

tmp <- RA_df[rev(p$tree_row$order),]
pickTF <- c(grep('FOS_extend',rownames(RA_df)),grep('DDIT3',rownames(RA_df)),grep('FOXC1',rownames(RA_df)),grep('FOXO1',rownames(RA_df)),
            grep('HES1',rownames(RA_df)),grep('SOX5',rownames(RA_df)),grep('JUN ',rownames(RA_df)),
            grep('ATF3 ',rownames(RA_df)),grep('JUNB',rownames(RA_df)),grep('EGR1',rownames(RA_df)),grep('JUND',rownames(RA_df)),
            grep('SOX8',rownames(RA_df)),grep('NFATC2',rownames(RA_df)),grep('FOXD1',rownames(RA_df)),grep('HMGB2',rownames(RA_df)),
            grep('PITX1',rownames(RA_df)),grep('SOX9',rownames(RA_df)),grep('PPARG',rownames(RA_df)),grep('SMAD3',rownames(RA_df)),grep('HIVEP3',rownames(RA_df)),
            grep('KLF3',rownames(RA_df)),grep('ELF3',rownames(RA_df)),grep('HMGA2',rownames(RA_df)),grep('MSX1',rownames(RA_df)),grep('EZH2',rownames(RA_df)),
            grep('IRX2',rownames(RA_df)),grep('REL',rownames(RA_df)),grep('BLC3',rownames(RA_df)),grep('BHLHE22',rownames(RA_df)),
            grep('FOXP1',rownames(RA_df)),grep('CEBPB',rownames(RA_df)),grep('TWIST',rownames(RA_df))
            )
tmp <- RA_df[pickTF,]

p <- pheatmap::pheatmap(tmp, #fontsize_row=3,
                         color=colorRampPalette(c("navy","white","firebrick3"))(100), breaks=seq(-1.5, 1.5, length.out = 100),
                         treeheight_row=10, treeheight_col=10, border_color=NA,gaps_col = c(2),
                         cluster_cols = FALSE,cluster_rows = FALSE)
system('mkdir -p JingMA_NEW/res/SCENIC_main/FIG')
save_pheatmap_pdf(p,'JingMA_NEW/res/SCENIC_main/FIG/heatmap_celltype.pdf',height = 5.5,width = 4)

saveRDS(RA_df,'JingMA_NEW/res/SCENIC_main/FIG/RA_df.RDS')


########################
### 看小类"CSC","C0","C1","C2",'SSC','SC1',"SC2"
########################
cellInfo <- data.frame(Meta_df[Meta_df$cell %in% colnames(regulonAUC),c('cell','subcelltype')])
rownames(cellInfo) <- cellInfo$cell
colnames(cellInfo)[2] <- c('celltype')
AUC_df <- getAUC(regulonAUC)

RA_byCellType <- sapply(split(rownames(cellInfo),cellInfo$celltype),
                                     function(cells) rowMeans(AUC_df[,cells]))
saveRDS(RA_byCellType,'JingMA_NEW/res/SCENIC_main/FIG/RA_df_noScale.RDS')


RA_byCellType_Scale <- t(scale(t(RA_byCellType), center = T, scale=T))
RA_df <- RA_byCellType_Scale

p <- pheatmap::pheatmap(RA_df, #fontsize_row=3,
                        color=colorRampPalette(c("navy","white","firebrick3"))(100), breaks=seq(-1.5, 1.5, length.out = 100),
                        treeheight_row=10, treeheight_col=10, border_color=NA,
                        cluster_cols = FALSE,cluster_rows = TRUE)
saveRDS(RA_df,'JingMA_NEW/res/SCENIC_main/FIG/RA_df.RDS')


########################################################################
######################## 人工挑选TF #######################
########################################################################
pickTF <- c(grep('FOS_extend',rownames(RA_df)),grep('DDIT3',rownames(RA_df)),grep('FOXC1',rownames(RA_df)),grep('FOXO1',rownames(RA_df)),
            grep('HES1',rownames(RA_df)),grep('SOX5',rownames(RA_df)),grep('JUN ',rownames(RA_df)),grep('ATF3 ',rownames(RA_df)),grep('JUNB',rownames(RA_df)),
            grep('EGR1',rownames(RA_df)),grep('JUND',rownames(RA_df)),grep('SOX8',rownames(RA_df)),grep('NFATC2',rownames(RA_df)),grep('FOXD1',rownames(RA_df)),
            grep('HMGB2',rownames(RA_df)),grep('PITX1',rownames(RA_df)),grep('MYLK',rownames(RA_df)),grep('PPARG',rownames(RA_df)),grep('SMAD3',rownames(RA_df)),
            grep('HIVEP3',rownames(RA_df)),grep('KLF3',rownames(RA_df)),grep('ELF3',rownames(RA_df)),
            grep('SOX9',rownames(RA_df)),grep('HMGA2',rownames(RA_df)),
            grep('MSX1',rownames(RA_df)),grep('IRX2',rownames(RA_df)),grep('BLC3',rownames(RA_df)),grep('BHLHE22',rownames(RA_df)),
            grep('FOXP1',rownames(RA_df)),grep('TWIST',rownames(RA_df)),grep('CEBPB',rownames(RA_df)),grep('ETS2',rownames(RA_df)))
pickRA_df <- RA_df[pickTF,]
saveRDS(pickRA_df,'/home/yzj/JingMA_NEW/res/SCENIC_main/FIG/pickRA_df.RDS')

p <- pheatmap(pickRA_df, fontsize_row=10,
              color=colorRampPalette(c("navy","white","firebrick3"))(100), breaks=seq(-1.5, 1.5, length.out = 100),
              treeheight_row=10, treeheight_col=10, border_color=NA,gaps_col = c(4),
              cluster_cols = FALSE,cluster_rows = FALSE)
system('mkdir -p /home/yzj/JingMA_NEW/res/SCENIC_main/FIG')
save_pheatmap_pdf(p, '/home/yzj/JingMA_NEW/res/SCENIC_main/FIG/heatmap_celltype.pdf',height = 6.5,width = 5)


########################################################################
#################Fig2I limma检验细胞类型之间，挑选TF #######################
########################################################################
regulonAUC <- readRDS(file='/local/yzj/JingMA_NEW/res/SCENIC_main/int/3.4_regulonAUC.Rds')
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]

regulonAUC <- getAUC(regulonAUC)
regulonAUC <- t(scale(t(regulonAUC), center = T, scale=T))

cellInfo <- read.table('JingMA_NEW/res/CellPhoneDB_merge_subcelltype/in/meta.txt',header = T,sep='\t',stringsAsFactors = F)
cellInfo <- cellInfo[cellInfo$Cell %in% colnames(regulonAUC),]
print(all(cellInfo$Cell==colnames(regulonAUC)))

design <- model.matrix(~0+factor(cellInfo$CellType))
colnames(design)=levels(factor(cellInfo$CellType))
rownames(design)=colnames(regulonAUC)

CellType=c('CSC','C0','C1','C2','SSC','SC1','SC2')
FC=1
DEG_mtx <- c()
allG_mtx <- c()
for(i in 1:length(CellType)){
  CT=CellType[i]
  x <- paste0(CT,'-(',paste0(CellType[-which(CellType==CT)],collapse = '+'),')')
  contrast.matrix<-makeContrasts(x,levels=design)
  
  fit <- lmFit(regulonAUC, design)
  fit2 <- contrasts.fit(fit, contrast.matrix) 
  fit2 <- eBayes(fit2) 
  tempOutput=topTable(fit2, adjust="fdr",n=Inf)
  nrDEG = na.omit(tempOutput)
  nrDEG = tibble::rownames_to_column(nrDEG,'Regulon')
  nrDEG$CellType=CT
  allG_mtx = rbind(allG_mtx,nrDEG)
  
  nrDEG = nrDEG[nrDEG$adj.P.Val < 0.05 & nrDEG$logFC > log(FC,2),]
  nrDEG =nrDEG[order(nrDEG$logFC,decreasing = T),]
  DEG_mtx <- rbind(DEG_mtx,nrDEG)
}
saveRDS(allG_mtx,'JingMA_NEW/res/SCENIC_main/FIG/Allregulon_mtx.RDS')
saveRDS(DEG_mtx,'JingMA_NEW/res/SCENIC_main/FIG/DEregulon_mtx.RDS')


DERegulon <- unique(DEG_mtx$Regulon)
## for CSC: delete FOXA3,MSX2,NFATC1
## for SC1: delete CREB5,CREM,REL
DERegulon <- DERegulon[!DERegulon %in% c('FOXA3 (11g)','MSX2_extended (435g)','NFATC1 (43g)','CREB5 (10g)','CREM (61g)','REL (834g)')]

RA_df <- readRDS('/home/yzj/JingMA_NEW/res/SCENIC_main/FIG/RA_df.RDS')
pickRA_df <- RA_df[rownames(RA_df)%in% DERegulon,]
saveRDS(pickRA_df,'/home/yzj/JingMA_NEW/res/SCENIC_main/FIG/DERA_df.RDS')


## 只在chond lineage里计算limma
regulonAUC <- readRDS(file='/local/yzj/JingMA_NEW/res/SCENIC_main/int/3.4_regulonAUC.Rds')
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]

regulonAUC <- getAUC(regulonAUC)
regulonAUC <- t(scale(t(regulonAUC), center = T, scale=T))

chondcells <- colnames(pbmc_main)[pbmc_main$celltype %in% c('CSPC','Chond')]
subregulonAUC <- regulonAUC[,colnames(regulonAUC) %in% chondcells]

cellInfo <- read.table('JingMA_NEW/res/CellPhoneDB_merge_subcelltype/in/meta.txt',header = T,sep='\t',stringsAsFactors = F)
cellInfo <- cellInfo[cellInfo$Cell %in% chondcells,]
print(all(cellInfo$Cell==colnames(subregulonAUC)))

design <- model.matrix(~0+factor(cellInfo$CellType))
colnames(design)=c('C0','C1','C2','CSC')
rownames(design)=colnames(subregulonAUC)

CellType=c('C0','C1','C2','CSC')
FC=1
DEG_mtx <- c()
allG_mtx <- c()
for(i in 1:length(CellType)){
  CT=CellType[i]
  x <- paste0(CT,'-(',paste0(CellType[-which(CellType==CT)],collapse = '+'),')')
  contrast.matrix<-makeContrasts(x,levels=design)
  
  fit <- lmFit(subregulonAUC, design)
  fit2 <- contrasts.fit(fit, contrast.matrix) 
  fit2 <- eBayes(fit2) 
  tempOutput=topTable(fit2, adjust="fdr",n=Inf)
  nrDEG = na.omit(tempOutput)
  nrDEG = tibble::rownames_to_column(nrDEG,'Regulon')
  nrDEG$CellType=CT
  allG_mtx = rbind(allG_mtx,nrDEG)
  
  nrDEG = nrDEG[nrDEG$adj.P.Val < 0.05 & nrDEG$logFC > log(FC,2),]
  nrDEG =nrDEG[order(nrDEG$logFC,decreasing = T),]
  DEG_mtx <- rbind(DEG_mtx,nrDEG)
}
saveRDS(allG_mtx,'JingMA_NEW/res/SCENIC_main/FIG/Allregulon_chond_mtx.RDS')
saveRDS(DEG_mtx,'JingMA_NEW/res/SCENIC_main/FIG/DEregulon_chond_mtx.RDS')


DERegulon <- unique(DEG_mtx$Regulon)
## for CSC: delete FOXA3,MSX2,NFATC1
## for SC1: delete CREB5,CREM,REL
DERegulon <- DERegulon[!DERegulon %in% c('FOXA3 (11g)','MSX2_extended (435g)','NFATC1 (43g)','CREB5 (10g)','CREM (61g)','REL (834g)')]

RA_df <- readRDS('/home/yzj/JingMA_NEW/res/SCENIC_main/FIG/RA_df.RDS')
pickRA_df <- RA_df[rownames(RA_df)%in% DERegulon,]
saveRDS(pickRA_df,'/home/yzj/JingMA_NEW/res/SCENIC_main/FIG/DERA_df.RDS')




############
### 计算差值
############
func_chazhi <- function(ct,df){
  index <- which(colnames(df)==ct)
  order_df <- t(apply(df,1,order))
  aim.auc <- df[,index]
  bg.auc <- df[,-index]
  bg.auc <- apply(bg.auc,1,max)
  
  maxTF <- rownames(order_df)[order_df[,ncol(df)]==index]
  bg.auc <- bg.auc[names(bg.auc)%in%maxTF]
  aim.auc <- df[rownames(df) %in% maxTF,1]
  order_chazhi <- maxTF[order(aim.auc-bg.auc,decreasing = T)]
  return(order_chazhi)
}

diffTF.CSC <- func_chazhi('CSC',pickRA_df);diffTF.C0 <- func_chazhi('C0',pickRA_df);diffTF.C1 <- func_chazhi('C1',pickRA_df);diffTF.C2 <- func_chazhi('C2',pickRA_df)
diffTF.SSC <- func_chazhi('SSC',pickRA_df);diffTF.SC1 <- func_chazhi('SC1',pickRA_df);diffTF.SC2 <- func_chazhi('SC2',pickRA_df)
diffTF <- unique(c(diffTF.CSC,diffTF.C0,diffTF.C1,diffTF.C2,diffTF.SSC,diffTF.SC1,diffTF.SC2))

################
## 差值大+ 是maker gene,两者排序综合
################
Idents(pbmc_main) <- pbmc_main$subcelltype
diffTF.lst <- list(diffTF.CSC,diffTF.C0,diffTF.C1,diffTF.C2,diffTF.SSC,diffTF.SC1,diffTF.SC2)
names(diffTF.lst) <- levels(pbmc_main$subcelltype)
MK <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/Markers_celltype_Main_mtx.RDS')

diffTF_new <- apply(as.matrix(diffTF),1,function(x) unlist(strsplit(x,'[ _extend]'))[1])
averEXP <- AverageExpression(pbmc_main,features = diffTF_new,slot = 'data')
averEXP <- averEXP$RNA

func_mkTF <- function(ct,MK){
  diffTF <- apply(as.matrix(diffTF.lst[[ct]]),1,function(x) unlist(strsplit(x,'[ _extend]'))[1])
  # sigMK <- MK[MK$avg_logFC>log(1.5)&MK$p_val_adj < 0.05&MK$CT==ct,]
  # sigMK <- MK[MK$CT==ct,]
  # sigmkTF <- (sigMK[sigMK$gene%in%diffTF,])
  EXP <- averEXP[rownames(averEXP)%in%diffTF,]
  sigmkTF <- EXP[order(EXP[,ct],decreasing = TRUE),]
  sigmkTF$gene <- rownames(sigmkTF) 
  #return(sigmkTF)
  if(nrow(sigmkTF) > 0){
    MK_rank <- 1:nrow(sigmkTF)
    names(MK_rank) <- sigmkTF$gene
    
    Diff_rank <- 1:length(diffTF)
    names(Diff_rank) <- diffTF
    
    order_df <- data.frame(CellType=ct,TF=sigmkTF$gene,
                           Diff_rank=Diff_rank[sigmkTF$gene],MK_rank=MK_rank[sigmkTF$gene])
    order_df$sum_rank <- order_df$Diff_rank+order_df$MK_rank
    order_df <- order_df[order(order_df$sum_rank),]
    order_df$sum_rank <- 1:nrow(order_df)
    return(order_df)
  }else{
    print('no sigmkTF')
  }
  
}

library(xlsx)
pickTF.CSC <- func_mkTF('CSC',MK);pickTF.C0 <- func_mkTF('C0',MK);pickTF.C1 <- func_mkTF('C1',MK);pickTF.C2 <- func_mkTF('C2',MK)
pickTF.SSC <- func_mkTF('SSC',MK);pickTF.SC1 <- func_mkTF('SC1',MK);pickTF.SC2 <- func_mkTF('SC2',MK)

topTF <- rbind(pickTF.CSC,pickTF.C0,pickTF.C1,pickTF.C2,pickTF.SSC,pickTF.SC1,pickTF.SC2)
write.xlsx(topTF,'JingMA_NEW/res/SCENIC_main/FIG/topTF.xlsx',sheetName = 'TF_rank')

TF2Regulon <- rownames(RA_df)
names(TF2Regulon) <-apply(as.matrix(rownames(RA_df)),1,function(x) unlist(strsplit(x,"[ _]"))[1])

top_TF <- c(pickTF.CSC$TF[1:5],pickTF.C0$TF[1:5],pickTF.C1$TF[1:5],pickTF.C2$TF[1:5],
            pickTF.SSC$TF[1:5],pickTF.SC1$TF[1:5],pickTF.SC2$TF[1:5])
top_TF <- na.omit(top_TF)
pickRA_df <- RA_df[rownames(RA_df) %in% TF2Regulon[top_TF],]

library(pheatmap)
order_TF <- c("FOXO1 (337g)" ,"KLF10 (16g)" ,"HES1 (37g)","ATF3 (83g)","JUN (34g)","RARG_extended (26g)" ,"SOX8 (158g)",
              "NFATC2 (227g)","HMGB2 (574g)","MYLK (73g)","MAF_extended (12g)","PPARG (22g)","NR1D1_extended (21g)","HIVEP3 (303g)",
              "SOX9 (20g)","ZMIZ1_extended (19g)" ,"HMGA2 (213g)","MECOM_extended (10g)", "SIX1 (18g)" ,"KLF11 (13g)" ,"SNAI1_extended (16g)",
              "KLF6_extended (23g)","TBX2_extended (20g)","GLIS2_extended (18g)","FOXP1 (319g)","MSC_extended (36g)" ,"TRPS1 (622g)",
              "CEBPD (199g)", "CEBPB (258g)","ETS2 (16g)","ATF5 (20g)")

pickRA_df  <- dplyr::select(as.data.frame(t(pickRA_df)),order_TF)
pickRA_df <- as.data.frame(t(pickRA_df))

p <- pheatmap(pickRA_df, clustering_method = 'mcquitty',fontsize = 5,fontsize_col = 5,fontsize_row = 5,
         color=colorRampPalette(c("navy","white","firebrick3"))(100), breaks=seq(-1.5, 1.5, length.out = 100),
         treeheight_row=10, treeheight_col=10, border_color=NA,gaps_col = c(4),
         cluster_cols = FALSE,cluster_rows = FALSE)

system('mkdir -p /home/yzj/JingMA_NEW/res/SCENIC_main/FIG')
ggsave('/home/yzj/JingMA_NEW/res/SCENIC_main/FIG/Fig2I_heatmap_subcelltype.pdf',
       p,width = 6,height = 7,units = 'cm')


###########################
### Fig2j 转录因子 feature plot
###########################
library(Seurat)
library(SCENIC)
pbmc <- readRDS('/home/yzj/JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype.Rdata')
Idents(pbmc) <- pbmc$celltype
pbmc_main <- subset(pbmc,idents = c('CSC','C','SSC','SC'))
pickRA_df <- readRDS('/home/yzj/JingMA_NEW/res/SCENIC_main/FIG/pickRA_df.RDS')
TF <- rownames(pickRA_df)


regulonAUC <- readRDS("~/JingMA_NEW/res/SCENIC_main/int/4.1_binaryRegulonActivity.Rds")
regulonAUC <- as.data.frame(regulonAUC)

Dim_df <- pbmc_main@reductions$umap@cell.embeddings

for(gene in TF){
  Gene= gene
  print(Gene)
  Cell2AUC=as.numeric(regulonAUC[rownames(regulonAUC)==Gene,])
  Active_cell <-colnames(regulonAUC)[which(Cell2AUC==1)]
  print(length(Active_cell))
  
  plot.df <- as.data.frame(Dim_df)
  plot.df$AUC <- 0
  plot.df$AUC[rownames(plot.df)%in%Active_cell] <-1
  
  plot.df_0 <- plot.df[plot.df$AUC==0,]
  sample_cells <- sample(rownames(plot.df_0),5000)
  plot.df_0 <- plot.df_0[rownames(plot.df_0)%in%sample_cells,]
  
  plot.df_1 <- plot.df[plot.df$AUC==1,]
  plot.df <- rbind(plot.df_0,plot.df_1)
  
  library(ggplot2)
  library(scales)
  p <- ggplot(data = plot.df,aes(x=UMAP_1,y=UMAP_2,colour=AUC))+
    geom_point(size=0.7)+
    scale_color_gradient(low='grey',high='firebrick')+
    theme_minimal()+
    theme(panel.grid = element_blank(),
          axis.title.x = element_text(size = 0),
          axis.title.y = element_text(size = 0),
          axis.text = element_text(size = 0))+
    labs(title=paste(unlist(strsplit(Gene,'[ _extend]'))[1]," regulon",sep=''))+
    theme(plot.title = element_text(hjust = 0.5))  #也就加上这一行
  ggsave(paste('/home/yzj/JingMA_NEW/res/SCENIC_main/FIG/scenic/',Gene,'_scenic.pdf',sep=''),p,width = 4.5,height = 4)
}

## 调整SOX9
regulonAUC_conti <- readRDS("~/JingMA_NEW/res/SCENIC_main/int/3.4_regulonAUC.Rds")
regulonAUC_conti <- as.data.frame(regulonAUC_conti@assays@data@listData)

Gene= 'SOX9 (20g)'
print(Gene)
Cell2AUC=as.numeric(regulonAUC_conti[rownames(regulonAUC_conti)==Gene,])
Active_cell <-colnames(regulonAUC_conti)[which(Cell2AUC>0.27)]
Active_cell <- apply(as.matrix(Active_cell),1, function(x) sub('\\.','-',unlist(strsplit(x,'AUC.'))[2]))

plot.df <- as.data.frame(Dim_df)
plot.df$AUC <- 0
plot.df$AUC[rownames(plot.df)%in%Active_cell] <-1

plot.df_0 <- plot.df[plot.df$AUC==0,]
sample_cells <- sample(rownames(plot.df_0),5000)
plot.df_0 <- plot.df_0[rownames(plot.df_0)%in%sample_cells,]

plot.df_1 <- plot.df[plot.df$AUC==1,]
plot.df <- rbind(plot.df_0,plot.df_1)

library(ggplot2)
library(scales)
p <- ggplot(data = plot.df,aes(x=UMAP_1,y=UMAP_2,colour=AUC))+
  geom_point(size=0.7)+
  scale_color_gradient(low='grey',high='firebrick')+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 0),
        axis.title.y = element_text(size = 0),
        axis.text = element_text(size = 0))+
  labs(title=paste(Gene," regulon",sep=''))+
  theme(plot.title = element_text(hjust = 0.5))  #也就加上这一行
p

ggsave(paste('/home/yzj/JingMA_NEW/res/SCENIC_main/FIG/',Gene,'_scenic.pdf',sep=''),p,width = 4.5,height = 4)



## 调整FOXP1
regulonAUC_conti <- readRDS("~/JingMA_NEW/res/SCENIC_main/int/3.4_regulonAUC.Rds")
regulonAUC_conti <- as.data.frame(regulonAUC_conti@assays@data@listData)

Gene= 'FOXP1 (319g)'
print(Gene)
Cell2AUC=as.numeric(regulonAUC_conti[rownames(regulonAUC_conti)==Gene,])
Active_cell <-colnames(regulonAUC_conti)[which(Cell2AUC>0.26)]
Active_cell <- apply(as.matrix(Active_cell),1, function(x) sub('\\.','-',unlist(strsplit(x,'AUC.'))[2]))

plot.df <- as.data.frame(Dim_df)
plot.df$AUC <- 0
plot.df$AUC[rownames(plot.df)%in%Active_cell] <-1

plot.df_0 <- plot.df[plot.df$AUC==0,]
sample_cells <- sample(rownames(plot.df_0),5000)
plot.df_0 <- plot.df_0[rownames(plot.df_0)%in%sample_cells,]

plot.df_1 <- plot.df[plot.df$AUC==1,]
plot.df <- rbind(plot.df_0,plot.df_1)

library(ggplot2)
library(scales)
p <- ggplot(data = plot.df,aes(x=UMAP_1,y=UMAP_2,colour=AUC))+
  geom_point(size=0.7)+
  scale_color_gradient(low='grey',high='firebrick')+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 0),
        axis.title.y = element_text(size = 0),
        axis.text = element_text(size = 0))+
  labs(title=paste(Gene," regulon",sep=''))+
  theme(plot.title = element_text(hjust = 0.5))  #也就加上这一行
p

ggsave(paste('/home/yzj/JingMA_NEW/res/SCENIC_main/FIG/',Gene,'_scenic.pdf',sep=''),p,width = 4.5,height = 4)


## 调整SOX5
regulonAUC_conti <- readRDS("~/JingMA_NEW/res/SCENIC_main/int/3.4_regulonAUC.Rds")
regulonAUC_conti <- as.data.frame(regulonAUC_conti@assays@data@listData)

Gene= 'SOX5 (218g)'
print(Gene)
Cell2AUC=as.numeric(regulonAUC_conti[rownames(regulonAUC_conti)==Gene,])
Active_cell <-colnames(regulonAUC_conti)[which(Cell2AUC>0.21)]
Active_cell <- apply(as.matrix(Active_cell),1, function(x) sub('\\.','-',unlist(strsplit(x,'AUC.'))[2]))

plot.df <- as.data.frame(Dim_df)
plot.df$AUC <- 0
plot.df$AUC[rownames(plot.df)%in%Active_cell] <-1

plot.df_0 <- plot.df[plot.df$AUC==0,]
sample_cells <- sample(rownames(plot.df_0),5000)
plot.df_0 <- plot.df_0[rownames(plot.df_0)%in%sample_cells,]

plot.df_1 <- plot.df[plot.df$AUC==1,]
plot.df <- rbind(plot.df_0,plot.df_1)

library(ggplot2)
library(scales)
p <- ggplot(data = plot.df,aes(x=UMAP_1,y=UMAP_2,colour=AUC))+
  geom_point(size=0.7)+
  scale_color_gradient(low='grey',high='firebrick')+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 0),
        axis.title.y = element_text(size = 0),
        axis.text = element_text(size = 0))+
  labs(title=paste(Gene," regulon",sep=''))+
  theme(plot.title = element_text(hjust = 0.5))  #也就加上这一行
p

ggsave(paste('/home/yzj/JingMA_NEW/res/SCENIC_main/FIG/',Gene,'_scenic.pdf',sep=''),p,width = 4.5,height = 4)



## 调整TWIST
regulonAUC_conti <- readRDS("~/JingMA_NEW/res/SCENIC_main/int/3.4_regulonAUC.Rds")
regulonAUC_conti <- as.data.frame(regulonAUC_conti@assays@data@listData)

Gene= 'TWIST1 (124g)'
print(Gene)
Cell2AUC=as.numeric(regulonAUC_conti[rownames(regulonAUC_conti)==Gene,])
Active_cell <-colnames(regulonAUC_conti)[which(Cell2AUC>0.21)]
Active_cell <- apply(as.matrix(Active_cell),1, function(x) sub('\\.','-',unlist(strsplit(x,'AUC.'))[2]))

plot.df <- as.data.frame(Dim_df)
plot.df$AUC <- 0
plot.df$AUC[rownames(plot.df)%in%Active_cell] <-1

plot.df_0 <- plot.df[plot.df$AUC==0,]
sample_cells <- sample(rownames(plot.df_0),5000)
plot.df_0 <- plot.df_0[rownames(plot.df_0)%in%sample_cells,]

plot.df_1 <- plot.df[plot.df$AUC==1,]
plot.df <- rbind(plot.df_0,plot.df_1)

library(ggplot2)
library(scales)
p <- ggplot(data = plot.df,aes(x=UMAP_1,y=UMAP_2,colour=AUC))+
  geom_point(size=0.7)+
  scale_color_gradient(low='grey',high='firebrick')+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 0),
        axis.title.y = element_text(size = 0),
        axis.text = element_text(size = 0))+
  labs(title=paste(Gene," regulon",sep=''))+
  theme(plot.title = element_text(hjust = 0.5))  #也就加上这一行
p

ggsave(paste('/home/yzj/JingMA_NEW/res/SCENIC_main/FIG/',Gene,'_scenic.pdf',sep=''),p,width = 4.5,height = 4)


## 调整ETS2
regulonAUC_conti <- readRDS("~/JingMA_NEW/res/SCENIC_main/int/3.4_regulonAUC.Rds")
regulonAUC_conti <- as.data.frame(regulonAUC_conti@assays@data@listData)

Gene= 'ETS2 (16g)'
print(Gene)
Cell2AUC=as.numeric(regulonAUC_conti[rownames(regulonAUC_conti)==Gene,])
Active_cell <-colnames(regulonAUC_conti)[which(Cell2AUC>0.26)]
Active_cell <- apply(as.matrix(Active_cell),1, function(x) sub('\\.','-',unlist(strsplit(x,'AUC.'))[2]))

plot.df <- as.data.frame(Dim_df)
plot.df$AUC <- 0
plot.df$AUC[rownames(plot.df)%in%Active_cell] <-1

plot.df_0 <- plot.df[plot.df$AUC==0,]
sample_cells <- sample(rownames(plot.df_0),6000)
plot.df_0 <- plot.df_0[rownames(plot.df_0)%in%sample_cells,]

plot.df_1 <- plot.df[plot.df$AUC==1,]
plot.df <- rbind(plot.df_0,plot.df_1)

library(ggplot2)
library(scales)
p <- ggplot(data = plot.df,aes(x=UMAP_1,y=UMAP_2,colour=AUC))+
  geom_point(size=0.7)+
  scale_color_gradient(low='grey',high='firebrick')+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 0),
        axis.title.y = element_text(size = 0),
        axis.text = element_text(size = 0))+
  labs(title=paste(Gene," regulon",sep=''))+
  theme(plot.title = element_text(hjust = 0.5))  #也就加上这一行
p

ggsave(paste('/home/yzj/JingMA_NEW/res/SCENIC_main/FIG/',Gene,'_scenic.pdf',sep=''),p,width = 4.5,height = 4)


## 调整HMGA2
regulonAUC_conti <- readRDS("~/JingMA_NEW/res/SCENIC_main/int/3.4_regulonAUC.Rds")
regulonAUC_conti <- as.data.frame(regulonAUC_conti@assays@data@listData)

Gene= 'HMGA2 (213g)'
print(Gene)
Cell2AUC=as.numeric(regulonAUC_conti[rownames(regulonAUC_conti)==Gene,])
Active_cell <-colnames(regulonAUC_conti)[which(Cell2AUC>0.27)]
Active_cell <- apply(as.matrix(Active_cell),1, function(x) sub('\\.','-',unlist(strsplit(x,'AUC.'))[2]))

plot.df <- as.data.frame(Dim_df)
plot.df$AUC <- 0
plot.df$AUC[rownames(plot.df)%in%Active_cell] <-1

plot.df_0 <- plot.df[plot.df$AUC==0,]
sample_cells <- sample(rownames(plot.df_0),5000)
plot.df_0 <- plot.df_0[rownames(plot.df_0)%in%sample_cells,]

plot.df_1 <- plot.df[plot.df$AUC==1,]
plot.df <- rbind(plot.df_0,plot.df_1)

library(ggplot2)
library(scales)
p <- ggplot(data = plot.df,aes(x=UMAP_1,y=UMAP_2,colour=AUC))+
  geom_point(size=0.7)+
  scale_color_gradient(low='grey',high='firebrick')+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 0),
        axis.title.y = element_text(size = 0),
        axis.text = element_text(size = 0))+
  labs(title=paste(Gene," regulon",sep=''))+
  theme(plot.title = element_text(hjust = 0.5))  #也就加上这一行
p

ggsave(paste('/home/yzj/JingMA_NEW/res/SCENIC_main/FIG/',Gene,'_scenic.pdf',sep=''),p,width = 4.5,height = 4)



#########
#### 补充材料: Vlnplot
#########
library(Seurat)
library(ggplot2)
library(ggpubr)
pbmc <- readRDS('/home/yzj/JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype.Rdata')
Idents(pbmc) <- pbmc$celltype
pbmc_main <- subset(pbmc,idents = c('CSC','C','SSC','SC'))
Color <- c("#EE9572","#B2DF8A" ,"#A6CEE3","#9999FF","#BC80BD", "#80B1D3", "#F4A460")

p <- DimPlot(pbmc_main,group.by = 'subcelltype',cols =Color,label = T,pt.size = 0.2,label.size = 8)
ggsave("/home/yzj/JingMA_NEW/res/Harmony/ALL/FIG/sFig_subcelltypes.pdf",p,width = 9.5,height = 9)


GeneLst <- c('SOX8','NFATC2','SOX9','FOXP1')
subpbmc <- subset(pbmc,feature=GeneLst)
EXP <- t(as.data.frame(subpbmc@assays$RNA@data))
tmp <- apply(EXP,2,function(x) x[x>3]=3)

plot.list <- list()
for(i in 1:length(GeneLst)){
  Gene=GeneLst[i]
  p <- VlnPlot(pbmc_main,Gene,pt.size = 0,group.by = 'celltype',cols = Color,y.max = 3)+
    theme(axis.text = element_text(size=15,face = 'bold',colour = 'black'),
          axis.title = element_text(size=15,face = 'bold',colour = 'black'),
          panel.background=element_rect(fill='transparent', color='black',size = 1.5),
          legend.position = 'none')+
    labs(x="")
  plot.list[[i]] <- p
}

p <- ggarrange(plotlist =plot.list,ncol=1,nrow = 4)
ggsave('JingMA_NEW/res/SCENIC_main/FIG/sFig_TF_EXP.pdf',p,width = 4,height = 14)


############
#### 补充材料,表4
############
library(xlsx)
regulonTargetsInfo <- readRDS("~/JingMA_NEW/res/SCENIC_main/int/2.5_regulonTargetsInfo.Rds")
head(regulonTargetsInfo)
TF <- c('SOX8','NFATC2','SOX9','FOXP1','HMGA2','JUN','ETS2','EGR1','HES1','SOX5','HIVEP3','PPARG',
        'TWIST1','CEBPB')
regulonTargetsInfo <- as.data.frame(regulonTargetsInfo)
pickTF_G <- regulonTargetsInfo[regulonTargetsInfo$TF %in% TF,]
pickTF_G <- pickTF_G[pickTF_G$highConfAnnot==TRUE,]
dim(pickTF_G)
write.xlsx(pickTF_G,'JingMA_NEW/res/SCENIC_main/FIG/regulonTargetsInfo.xlsx',row.names = FALSE)
