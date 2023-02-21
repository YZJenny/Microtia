library(Seurat)
library(ggplot2)
library(SCENIC)
require("RColorBrewer")
library(maSigPro)

pbmc_C <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype_Control_Chond.Rdata')
pbmc_C$celltype <- factor(pbmc_C$celltype,levels = c('CSC','C0','C1','C2'))

sample.cells <- sample(colnames(pbmc_C),5000)
sample.pbmc_C <- subset(pbmc_C,cells = sample.cells)
# plot single gene
Harmony2 <- sample.pbmc_C@reductions$harmony@cell.embeddings[, 'harmony_2']
mat.gene <- sample.pbmc_C@assays$RNA@data
# AUC
regulonAUC <- readRDS(file='/home/yzj/JingMA_NEW/res/SCENIC_main/int/3.4_regulonAUC.Rds')
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)), colnames(mat.gene)]
mat.auc <- as.matrix(regulonAUC@assays@data@listData$AUC)
df.pc.gene <- data.frame(t(rbind(as.matrix(mat.gene), mat.auc)), check.names = F)
df.pc.gene$Harmony2 <- Harmony2
df.pc.gene$celltype <- sample.pbmc_C$celltype
df.pc.gene$status <- sample.pbmc_C$Phase
df.pc.gene <- df.pc.gene[order(Harmony2, decreasing = F),]
df.pc.gene$idx <- 1:nrow(df.pc.gene)
colors <- c("#EE9572","#B2DF8A" ,"#A6CEE3","#9999FF")
names(colors) <- c('CSC', 'C0', 'C1', 'C2')

########################
### Fig 3E: 比较成人和儿童 TF AUC 随发育时间的变化
########################
vec.TF <- c('CEBPB (258g)','CEBPD (199g)')
type='up'
for (i in 1:length(vec.TF)) {
  TF <- vec.TF[i]
  df.plot <- df.pc.gene[, c('idx', 'Harmony2', 'celltype', 'status', TF)]
  names(df.plot) <- c('idx', 'Harmony2', 'celltype', 'status', 'TF')
  p1 <- ggplot(data = df.plot, aes(x = idx, linetype = status,y = TF)) + 
    geom_point(aes(color = celltype), size = 0.0000001) +  theme_classic()+
    scale_color_manual(labels = c('CSC', 'C0', 'C1', 'C2'),values = colors) + 
    geom_smooth(color = '#696969',size=0.5) +
    labs(x = '', y = '') + 
    theme(panel.background=element_rect(fill='transparent', color='black',size = 0.5),
          axis.text = element_blank(),axis.ticks = element_blank(),plot.margin = unit(c(0.2,0.1,-0.5,-0.5), "cm"),
          legend.position = 'none') +
    annotate('text', label = TF, x = 5000, y = max(df.plot$TF), hjust = 1, vjust = 1, size = 2)
  
  if (i == 1) {
    p <- p1
  } else {
    p <- p / p1
  }
}
ggsave(paste('/home/yzj/JingMA_NEW/res/compControl/ChildrenvsAdults/FIG/lineage_AUC/Fig3E_TF_AUC_',type,'.pdf',sep=''),p,
       width = 3,height = 4.5, units = 'cm')


########################
### Fig 3E: 比较成人和儿童的gene expression随发育时间的变化
########################
vec.TF.exp <- c('CEBPB','CEBPD')
type='up'
for (i in 1:length(vec.TF.exp)) {
  TF <- vec.TF.exp[i]
  df.plot <- df.pc.gene[, c('idx', 'Harmony2', 'celltype', 'status', TF)]
  names(df.plot) <- c('idx', 'Harmony2', 'celltype', 'status', 'TF')
  p1 <- ggplot(data = df.plot, aes(x = idx,linetype = status, y = TF)) + 
    geom_point(aes(color = celltype), size = 0.0000001) + theme_classic()+
    scale_color_manual(labels = c('CSC', 'C0', 'C1', 'C2'),values = colors) + 
    geom_smooth(color = '#696969',size=0.5) +
    labs(x = '', y = '') + 
    theme(panel.background=element_rect(fill='transparent', color='black',size = 0.5),plot.margin = unit(c(0.2,0.1,-0.5,-0.5), "cm"),
          axis.text = element_blank(),axis.ticks = element_blank(),legend.position = 'none') +
    annotate('text', label = TF,x = 5000, y = max(df.plot$TF), hjust = 1, vjust = 1, size = 2)
  if (i == 1) {
    p <- p1
  } else {
    p <- p / p1
  }
}
ggsave(paste('/home/yzj/JingMA_NEW/res/compControl/ChildrenvsAdults/FIG/lineage_EXP/Fig3E_TF_EXP_',type,'.pdf',sep=''),p,
       height = 4.5, width = 3, units = 'cm')



####################
## Fig 3E.
####################
library(Seurat)
library(ggplot2)
require("RColorBrewer")

pbmc_C <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype_Control_Chond.Rdata')
pbmc_C$celltype <- factor(pbmc_C$celltype,levels = c('CSC','C0','C1','C2'))

###############
###############
# plot single gene
Harmony2 <- pbmc_C@reductions$harmony@cell.embeddings[, 'harmony_2']
length(Harmony2)

# EXP
mat.auc <- as.data.frame(pbmc_C@assays$RNA@data)

df.pc.gene <- dplyr::select(mat.auc,names(Harmony2))
df.pc.gene <- as.data.frame(t(df.pc.gene))
print(all(rownames(df.pc.gene) == names(Harmony2)))

df.pc.gene$Harmony2 <- Harmony2
df.pc.gene$celltype <- pbmc_C$celltype
df.pc.gene$phase <- pbmc_C$Phase

df.pc.gene <- df.pc.gene[order(Harmony2, decreasing = F),]
df.pc.gene$idx <- 1:nrow(df.pc.gene)
colors <- c("#EE9572","#B2DF8A" ,"#A6CEE3","#9999FF")
names(colors) <- c('CSC', 'C0', 'C1', 'C2')


gene.lst=c('CEBPB','CEBPD')
for(i in 1:length(gene.lst)){
  gene <- gene.lst[i]
  print(gene)
  df.plot <- df.pc.gene[, c('idx', 'Harmony2', 'celltype', 'phase', gene)]
  names(df.plot) <- c('idx', 'Harmony2', 'celltype', 'phase', 'gene')
  p.gene <-
    ggplot(data = df.plot, aes(x = idx,linetype = phase,y = gene)) + 
    geom_point(aes(color = celltype), size = 0.0000001) + 
    scale_color_manual(labels = c('CSC', 'C0', 'C1', 'C2'),values = colors) + 
    geom_smooth(color = '#696969') +
    labs(x = '', y = '') + 
    theme(panel.background=element_rect(fill='transparent', color='black',size = 1),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = 'none',plot.margin = unit(c(0.1,0.1,-0.5,-0.5), "cm")) +
    annotate('text', label = gene, x = 22000, y = max(df.plot$gene), hjust = 1, vjust = 1, size = 2)
  ggsave(paste('/home/yzj/JingMA_NEW/res/compControl/ChildrenvsAdults/FIG/lineage_EXP/',gene,'.pdf',sep=''),
         p.gene,width = 5,height = 5,units = 'cm')
}


####################
## limma比较成人和儿童的AUC Score
####################
library(Seurat)
library(SCENIC)
library(foreach)
library(AUCell)
library(limma)
library(reshape2)
library(tibble)
library(dplyr)

regulonAUC <- readRDS("~/JingMA_NEW/res/SCENIC_main/int/3.4_regulonAUC.Rds")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonAUC <- getAUC(regulonAUC)
regulonAUC <- t(scale(t(regulonAUC), center = T, scale=T))

pbmc_C <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype_Control_Chond.Rdata')
pbmc_C$celltype <- as.character(pbmc_C$celltype)
pbmc_C$celltype[pbmc_C$celltype=='TC'] <- 'C0'
pbmc_C$celltype <- factor(pbmc_C$celltype,levels = c('CSC','C0','C1','C2'))

cellInfo <-pbmc_C@meta.data[,c('Phase','celltype')]

get_DERA <- function(CT,FC){
  sub.cellInfo <- cellInfo[cellInfo$celltype==CT,]
  sub.regulonAUC <- dplyr::select(as.data.frame(regulonAUC),rownames(sub.cellInfo))
  print(all(rownames(sub.cellInfo)==colnames(sub.regulonAUC)))
  
  design <- model.matrix(~0+factor(sub.cellInfo$Phase))
  colnames(design)=levels(factor(sub.cellInfo$Phase))
  rownames(design)=colnames(sub.regulonAUC)
  
  x <- 'Adults-Children'
  contrast.matrix<-makeContrasts(x,levels=design)
  
  fit <- lmFit(sub.regulonAUC, design)
  fit2 <- contrasts.fit(fit, contrast.matrix) 
  fit2 <- eBayes(fit2) 
  tempOutput=topTable(fit2, adjust="fdr",n=Inf)
  nrDEG = na.omit(tempOutput)
  nrDEG = tibble::rownames_to_column(nrDEG,'Regulon')
  nrDEG = nrDEG[nrDEG$adj.P.Val < 0.05 & abs(nrDEG$logFC) > log(FC,2),]
  nrDEG =nrDEG[order(nrDEG$logFC,decreasing = T),]
  nrDEG$CellType=CT
  return(nrDEG)
}

CSC.DERA <- get_DERA('CSC',1.5)
C0.DERA <- get_DERA('C0',1.5)
C1.DERA <- get_DERA('C1',1.5)
C2.DERA <- get_DERA('C2',1.5)

DERA.lst <- list(CSC=CSC.DERA,C0=C0.DERA,C1=C1.DERA,C2=C2.DERA)
saveRDS(DERA.lst,'JingMA_NEW/res/compControl/ChildrenvsAdults/FIG/DERA_lst.RDS')

colnames(CSC.DERA)[2] <- 'CSC'
colnames(C0.DERA)[2] <- 'C0'
colnames(C1.DERA)[2] <- 'C1'
colnames(C2.DERA)[2] <- 'C2'

merge.DERA <- merge(CSC.DERA[,1:2],merge(C0.DERA[,1:2],merge(C1.DERA[,1:2],C2.DERA[,1:2],by='Regulon',all = TRUE),by='Regulon',all = TRUE),by='Regulon',all = TRUE)
merge.DERA[is.na(merge.DERA)] <- 0
rownames(merge.DERA) <- merge.DERA$Regulon
merge.DERA <- merge.DERA[,-1]

##
up_CSC <- CSC.DERA$Regulon[CSC.DERA[,2] > log(FC)]
dn_CSC <- CSC.DERA$Regulon[CSC.DERA[,2] < -log(FC) ]

up_C0 <- C0.DERA$Regulon[C0.DERA[,2] > log(FC)]
dn_C0 <- C0.DERA$Regulon[C0.DERA[,2] < -log(FC) ]

up_C1 <- C1.DERA$Regulon[C1.DERA[,2] > log(FC)]
dn_C1 <- C1.DERA$Regulon[C1.DERA[,2] < -log(FC) ]

up_C2 <- C2.DERA$Regulon[C2.DERA[,2] > log(FC)]
dn_C2 <- C2.DERA$Regulon[C2.DERA[,2] < -log(FC) ]


get_values <- function(sigCSC,sigTC,sigC1,sigC2){
  sigGene <- unique(c(sigCSC,sigTC,sigC1,sigC2))
  values <- matrix(c(rep(0,4*length(sigGene))),ncol = 4,dimnames = list(sigGene,c('CSC','C0','C1','C2')))
  for(i in 1:length(sigGene)){
    g=sigGene[i]
    if(g %in% sigCSC){values[i,1] <-1};
    if(g %in% sigTC){values[i,2] <-1};
    if(g %in% sigC1){values[i,3] <-1};
    if(g %in% sigC2){values[i,4] <-1};
  }
  values_sum <- apply(values, 1, sum)
  values <- values[order(values_sum,decreasing = T),]
  return(values)
}

## 对成人来说，上调矩阵
upValues_mtx <- get_values(up_CSC,up_C0,up_C1,up_C2)
up_sum <- apply(upValues_mtx,1,sum)
up_df <- upValues_mtx[-(which(up_sum>1)),]
annotation_col = data.frame(CellType = factor(c("CSC", "C0","C1","C2")))
rownames(annotation_col) <- colnames(upValues_mtx)
annotation_row = data.frame(GeneClass = factor(rep(c("Common", "CSC", "C0","C1","C2"), 
                                                   c(length(which(up_sum>1)), length(which(up_df[,1]==1)), length(which(up_df[,2]==1)),
                                                     length(which(up_df[,3]==1)),length(which(up_df[,4]==1))))))
rownames(annotation_row) = rownames(upValues_mtx)

ann_colors = list( CellType = c(CSC="#EE9572",C0="#B2DF8A",C1="#A6CEE3",C2="#9999FF"),
                   GeneClass = c(Common='grey',CSC="#EE9572",C0="#B2DF8A",C1="#A6CEE3",C2="#9999FF"))

order_upTF <- rownames(upValues_mtx)
merge.upDERA <- t(select(as.data.frame(t(merge.DERA)),order_upTF))

p_UP <- pheatmap(merge.upDERA,cluster_rows = F,cluster_cols = F,color =  colorRampPalette(c("white","firebrick3"))(100),
                 border_color ='transparent',show_rownames = T,angle_col='45',
                 annotation_row = annotation_row,annotation_colors = ann_colors,legend=T,annotation_legend = FALSE)

save_pheatmap_pdf(p_UP,'JingMA_NEW/res/compControl/ChildrenvsAdults/FIG/DERAHeatmap_UP.pdf',height = 8,width = 4)


## 对成人来说下调矩阵
dnValues_mtx <- get_values(dn_CSC,dn_C0,dn_C1,dn_C2)
dn_sum <- apply(dnValues_mtx,1,sum)
dn_df <- dnValues_mtx[-(which(dn_sum>1)),]

annotation_col = data.frame(CellType = factor(c("CSC", "C0","C1","C2")))
rownames(annotation_col) <- colnames(dnValues_mtx)

annotation_row = data.frame(GeneClass = factor(rep(c("Common", "CSC", "C0","C1","C2"), 
                                                   c(length(which(dn_sum>1)), length(which(dn_df[,1]==1)),length(which(dn_df[,2]==1)),
                                                     length(which(dn_df[,3]==1)),length(which(dn_df[,4]==1))))))
rownames(annotation_row) = rownames(dnValues_mtx)

ann_colors = list( CellType = c(CSC="#EE9572",C0="#B2DF8A",C1="#A6CEE3",C2="#9999FF"),
                   GeneClass = c(Common='grey',CSC="#EE9572",C0="#B2DF8A",C1="#A6CEE3",C2="#9999FF"))

order_dnTF <- rownames(dnValues_mtx)
merge.dnDERA <- t(select(as.data.frame(t(merge.DERA)),order_dnTF))

p_DN <- pheatmap(merge.dnDERA,cluster_rows = F,cluster_cols = F,color = colorRampPalette(c("navy",'white'))(100),
                 border_color ='transparent',show_rownames = T,legend=T,angle_col='45',
                 annotation_row = annotation_row,annotation_colors = ann_colors,annotation_legend = FALSE)
save_pheatmap_pdf(p_DN,'JingMA_NEW/res/compControl/ChildrenvsAdults/FIG/DERAHeatmap_DN.pdf',height = 5,width = 4)



### 单独画每个细胞类型的DERA
DERA.lst <- readRDS('JingMA_NEW/res/compControl/ChildrenvsAdults/FIG/DERA_lst.RDS')
CellType=c('CSC','C0','C1','C2')
for(i in 1:length(CellType)){
  CT=CellType[i]
  sub.cellInfo <- cellInfo[cellInfo$celltype==CT,]
  sub.regulonAUC <- dplyr::select(as.data.frame(regulonAUC),rownames(sub.cellInfo))
  
  print(all(rownames(sub.cellInfo)==colnames(sub.regulonAUC)))
  
  RA_byPhase <- sapply(split(rownames(sub.cellInfo),c(sub.cellInfo$Phase)),
                          function(cells) rowMeans(sub.regulonAUC[,cells]))
  RA_byPhase <- t(scale(t(RA_byPhase), center = T, scale=T))
  RA <- DERA.lst[[CT]][,1]
  
  RA_df <- RA_byPhase[rownames(RA_byPhase)%in%RA,]
  p <- pheatmap::pheatmap(RA_df, #fontsize_row=3,
                          color=colorRampPalette(c("navy","white","firebrick3"))(100), breaks=seq(-0.8, 0.8, length.out = 100),
                          treeheight_row=0, treeheight_col=0, border_color=NA,
                          cluster_cols = FALSE,cluster_rows = TRUE)
  save_pheatmap_pdf(p,paste('/home/yzj/JingMA_NEW/res/compControl/ChildrenvsAdults/FIG/',CT,'_AdultsChildren.pdf',sep=''),height = 5.5,width = 3)
}




####################
## 比较成人和儿童的AUC Score随发育时间的变化
####################
library(Seurat)
library(ggplot2)
library(SCENIC)
library(AUCell)
require("RColorBrewer")

pbmc_C <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype_Control_Chond.Rdata')
pbmc_C$celltype <- as.character(pbmc_C$celltype)
pbmc_C$celltype[pbmc_C$celltype=='TC'] <- 'C0'
pbmc_C$celltype <- factor(pbmc_C$celltype,levels = c('CSC','C0','C1','C2'))

###############
###############
# plot single gene
Harmony2 <- pbmc_C@reductions$harmony@cell.embeddings[, 'harmony_2']
length(Harmony2)

# AUC
regulonAUC <- readRDS(file='/home/yzj/JingMA_NEW/res/SCENIC_main/int/3.4_regulonAUC.Rds')
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
mat.auc <- getAUC(regulonAUC)
dim(mat.auc)

df.pc.gene <- dplyr::select(as.data.frame(mat.auc),names(Harmony2))
df.pc.gene <- as.data.frame(t(df.pc.gene))
print(all(rownames(df.pc.gene) == names(Harmony2)))

df.pc.gene$Harmony2 <- Harmony2
df.pc.gene$celltype <- pbmc_C$celltype
df.pc.gene$phase <- pbmc_C$Phase

df.pc.gene <- df.pc.gene[order(Harmony2, decreasing = F),]
df.pc.gene$idx <- 1:nrow(df.pc.gene)
colors <- c("#EE9572","#B2DF8A" ,"#A6CEE3","#9999FF")
names(colors) <- c('CSC', 'C0', 'C1', 'C2')


DERA.lst <- readRDS('JingMA_NEW/res/compControl/ChildrenvsAdults/FIG/DERA_lst.RDS')
gene.lst <- unique(c(DERA.lst[['CSC']][,1],DERA.lst[['C0']][,1],DERA.lst[['C1']][,1],DERA.lst[['C2']][,1]))

#gene <- colnames(df.pc.gene)[grepl('SOX8',colnames(df.pc.gene))]

gene.lst <- c('CEBPB (258g)','CEBPD (199g)')
for(i in 1:length(gene.lst)){
  gene <- gene.lst[i]
  print(gene)
  df.plot <- df.pc.gene[, c('idx', 'Harmony2', 'celltype', 'phase', gene)]
  names(df.plot) <- c('idx', 'Harmony2', 'celltype', 'phase', 'gene')
  p.gene <-
    ggplot(data = df.plot, aes(x = idx,linetype = phase,y = gene)) + 
    geom_point(aes(color = celltype), size = 0.3) + 
    scale_color_manual(labels = c('CSC', 'TC', 'C1', 'C2'),values = colors) + 
    geom_smooth(color = '#696969') +
    labs(x = '', y = '') + 
    theme(panel.background=element_rect(fill='transparent', color='black',size = 2),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = 'none',plot.margin = unit(c(0.1,0.1,-0.5,-0.5), "cm")) +
    annotate('text', label = gene, x = 22000, y = max(df.plot$gene), hjust = 1, vjust = 1, size = 7)
  ggsave(paste('/home/yzj/JingMA_NEW/res/compControl/ChildrenvsAdults/FIG/lineage_AUC/',gene,'.pdf',sep=''),p.gene,width = 4,height = 4)
}


