library(Seurat)
library(ggplot2)
library(dplyr)

pbmc_chond <- readRDS('/home/yzj/JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype_Chond.Rdata')
subpbmc <- subset(pbmc_chond,cells = colnames(pbmc_chond)[pbmc_chond$batch %in% c('C4','C6','M1','M2','M3')])
subpbmc$celltype <- factor(subpbmc$celltype,levels = rev(c('CSC','TC','C1','C2')))

############
### Fig4A. Proportion.
############
phylog_df <- subpbmc@meta.data[,c('batch',"celltype")]
phylog_df <- table(phylog_df$batch,phylog_df[,"celltype"])
phylog_df <- data.frame(phylog_df)
colnames(phylog_df) <- c('SampleID','CellType','Freq')
phylog_df$CellType <- factor(phylog_df$CellType)

Color <- c("#EE9572","#B2DF8A" ,"#A6CEE3","#9999FF")
p1 <- ggplot(phylog_df,aes(x=SampleID,y=Freq,fill=CellType))+
  geom_col(position = "fill", width = 0.6)+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text = element_text(size = 10,colour = "black"),
        axis.line = element_line(size=0.7, colour = "black"),
        axis.title.y = element_text(size=10),
        plot.title = element_text(hjust = 0.5))+
  labs(x='',y='% of chondrocytes',title = 'Per sample')+guides(fill=FALSE)+
  scale_fill_manual(values = rev(Color))
p1

phylog_df <- subpbmc@meta.data[,c('type',"celltype")]
phylog_df <- table(phylog_df$type,phylog_df[,"celltype"])
phylog_df <- data.frame(phylog_df)
colnames(phylog_df) <- c('Status','CellType','Freq')
phylog_df$CellType <- factor(phylog_df$CellType)

Color <- c("#EE9572","#B2DF8A" ,"#A6CEE3","#9999FF")
p2 <- ggplot(phylog_df,aes(x=Status,y=Freq,fill=CellType))+
  geom_col(position = "fill", width = 0.9)+
  theme_classic()+
  labs(x='',y='',title = 'Mean')+
  theme(legend.position="right",plot.title = element_text(hjust = 0.5))+
  theme(axis.text = element_text(size = 10,colour = "black"),
        axis.line = element_line(size=0.7, colour = "black"),
        axis.title.y = element_text(size=10))+
  labs(x='',y='')+theme(legend.position="right")+
  scale_fill_manual(values = rev(Color))+
  guides(fill = guide_legend(reverse=TRUE))
p2

library(ggpubr)
pdf('/home/yzj/JingMA_NEW/res/Harmony/ALL/FIG/Fig4A_Barplot_PropChond_Microtia.pdf',width = 7,height = 4)
ggarrange(p1, p2,ncol = 2, nrow = 1)
dev.off()

##################
###2. 计算markerDEG的ratio
##################
DEG.lst <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/DEGs_inChond_inChildren_NormalMicrotia.RDS')
print(names(DEG.lst))

MK.lst <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/Markers_celltype_Chond.RDS')
print(names(MK.lst))

up.ratio <- c()
dn.ratio <- c()
ratio <- c()
for(i in 1:length(names(MK.lst))){
  CT=names(MK.lst)[i]
  #print(CT)
  DEG <- DEG.lst[[CT]]
  up.DEG <- rownames(DEG)[DEG$avg_logFC > log(1.5,2) & DEG$p_val_adj < 0.05]
  dn.DEG <- rownames(DEG)[DEG$avg_logFC < -log(1.5,2) & DEG$p_val_adj < 0.05]
  DEG <- c(up.DEG,dn.DEG)
  
  MK <- MK.lst[[CT]]
  MK <- rownames(MK)[MK$avg_logFC > log(1.5,2) & MK$p_val_adj < 0.05]
  
  # up.ratio <- c(round(length(intersect(up.DEG,MK))/length(up.DEG),2),up.ratio)
  # dn.ratio <- c(round(length(intersect(dn.DEG,MK))/length(dn.DEG),2),dn.ratio)
  # ratio <- c(round(length(intersect(DEG,MK))/length(DEG),2),ratio)
  
  up.ratio <- c(round(length(intersect(up.DEG,MK))/length(MK),2),up.ratio)
  dn.ratio <- c(round(length(intersect(dn.DEG,MK))/length(MK),2),dn.ratio)
  ratio <- c(round(length(intersect(DEG,MK))/length(MK),2),ratio)
}

# Load ggplot2
library(ggplot2)

# Create data
data <- data.frame(
  celltype=c("CSPC","EC","IC","LC") ,  
  value=ratio
)

# Barplot
p <- ggplot(data, aes(x=celltype, y=value,fill=celltype)) + 
  geom_bar(stat = "identity")+theme_classic()+
  scale_fill_manual(values =  c("#EE9572","#B2DF8A" ,"#A6CEE3","#9999FF"))+
  labs(x="", y="Ratio of DE-marker gene")+guides(fill=FALSE)+
  theme(axis.text=element_text(size=8,colour="black"),
        axis.title=element_text(size = 8,colour="black"))
ggsave('JingMA_NEW/res/compMicrotia/MicrotiavsNormal_inChildren/FIG/Fig4B_Ratio.pdf',p,width = 6,height = 6,units = 'cm')


############
### Fig4B. Heatmap of DEG
############
MK.lst <- readRDS('/home/disk/drizzle/wgk/microtia_chon_child_M1M2M3/cutoff_0.4/marker_go.Rdata')

data <- MK.lst[['CSC']]
up_CSC <- rownames(data)[data$avg_logFC > log(1.5) & data$p_val_adj < 0.05]
dn_CSC <- rownames(data)[data$avg_logFC < -log(1.5) & data$p_val_adj < 0.05]

data <- MK.lst[["TC"]]
up_TC <- rownames(data)[data$avg_logFC > log(1.5) & data$p_val_adj < 0.05]
dn_TC <- rownames(data)[data$avg_logFC < -log(1.5) & data$p_val_adj < 0.05]

data <- MK.lst[["C1"]]
up_C1 <- rownames(data)[data$avg_logFC > log(1.5) & data$p_val_adj < 0.05]
dn_C1 <- rownames(data)[data$avg_logFC < -log(1.5) & data$p_val_adj < 0.05]

data <- MK.lst[["C2"]]
up_C2 <- rownames(data)[data$avg_logFC > log(1.5) & data$p_val_adj < 0.05]
dn_C2 <- rownames(data)[data$avg_logFC < -log(1.5) & data$p_val_adj < 0.05]

get_values <- function(sigCSC,sigTC,sigC1,sigC2){
  sigGene <- unique(c(sigCSC,sigTC,sigC1,sigC2))
  values <- matrix(c(rep(0,4*length(sigGene))),ncol = 4,dimnames = list(sigGene,c('CSC','TC','C1','C2')))
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

## 对疾病来说，上调矩阵
upValues_mtx <- get_values(up_CSC,up_TC,up_C1,up_C2)
up_sum <- apply(upValues_mtx,1,sum)
up_df <- upValues_mtx[-(which(up_sum>1)),]

annotation_col = data.frame(CellType = factor(c("CSC", "TC","C1","C2")))
rownames(annotation_col) <- colnames(upValues_mtx)

annotation_row = data.frame(GeneClass = factor(rep(c("Common", "CSC", "TC","C1","C2"), 
                                                   c(length(which(up_sum>1)), 
                                                     length(which(up_df[,1]==1)), 
                                                     length(which(up_df[,2]==1)),
                                                     length(which(up_df[,3]==1)),
                                                     length(which(up_df[,4]==1))))))
rownames(annotation_row) = rownames(upValues_mtx)

ann_colors = list( CellType = c(CSC="#EE9572",TC="#B2DF8A",C1="#A6CEE3",C2="#9999FF"),
                   GeneClass = c(Common='grey',CSC="#EE9572",TC="#B2DF8A",C1="#A6CEE3",C2="#9999FF"))

p_UP <- pheatmap(upValues_mtx,cluster_rows = F,cluster_cols = F,color =  colorRampPalette(c("#EFEFEF", "white","#B15E72"))(100),
                 border_color ='transparent',show_rownames = F,angle_col='45',
                 annotation_row = annotation_row,annotation_colors = ann_colors,legend=F,annotation_legend = FALSE)
save_pheatmap_pdf(p_UP,'JingMA_NEW/res/compMicrotia/MicrotiavsNormal_inChildren/FIG/DEGHeatmap_UP.pdf',height = 4,width = 2)
saveRDS(upValues_mtx,'JingMA_NEW/res/compMicrotia/MicrotiavsNormal_inChildren/FIG/DEGHeatmap_UPmtx.RDS')

## 对疾病来说，下调矩阵
dnValues_mtx <- get_values(dn_CSC,dn_TC,dn_C1,dn_C2)
dn_sum <- apply(dnValues_mtx,1,sum)
dn_df <- dnValues_mtx[-(which(dn_sum>1)),]

annotation_col = data.frame(CellType = factor(c("CSC", "TC","C1","C2")))
rownames(annotation_col) <- colnames(dnValues_mtx)

annotation_row = data.frame(GeneClass = factor(rep(c("Common", "CSC", "TC","C1","C2"), 
                                                   c(length(which(dn_sum>1)), 
                                                     length(which(dn_df[,1]==1)), 
                                                     length(which(dn_df[,2]==1)),
                                                     length(which(dn_df[,3]==1)),
                                                     length(which(dn_df[,4]==1))))))
rownames(annotation_row) = rownames(dnValues_mtx)

ann_colors = list( CellType = c(CSC="#EE9572",TC="#B2DF8A",C1="#A6CEE3",C2="#9999FF"),
                   GeneClass = c(Common='grey',CSC="#EE9572",TC="#B2DF8A",C1="#A6CEE3",C2="#9999FF"))

p_DN <- pheatmap(dnValues_mtx,cluster_rows = F,cluster_cols = F,color =  colorRampPalette(c("#EFEFEF", "white","#7F99CE"))(100),
                 border_color ='transparent',show_rownames = F,legend=F,angle_col='45',
                 annotation_row = annotation_row,annotation_colors = ann_colors,annotation_legend = FALSE)
save_pheatmap_pdf(p_DN,'JingMA_NEW/res/compMicrotia/MicrotiavsNormal_inChildren/FIG/DEGHeatmap_DN.pdf',height = 4,width = 2)
saveRDS(dnValues_mtx,'JingMA_NEW/res/compMicrotia/MicrotiavsNormal_inChildren/FIG/DEGHeatmap_DNmtx.RDS')


up_sum <- apply(upValues_mtx,1,sum)
length(which(up_sum==4))
length(which(up_sum>1))
up_df <- upValues_mtx[-(which(up_sum>1)),]
length(which(up_df[,1]==1))
length(which(up_df[,2]==1))
length(which(up_df[,3]==1))
length(which(up_df[,4]==1))

dn_sum <- apply(dnValues_mtx,1,sum)
length(which(dn_sum==4))
length(which(dn_sum>1))
dn_df <- dnValues_mtx[-(which(dn_sum>1)),]
length(which(dn_df[,1]==1))
length(which(dn_df[,2]==1))
length(which(dn_df[,3]==1))
length(which(dn_df[,4]==1))


###############
### Fig4C. Heatmap of GO
###############
fc.cutoff <- 0.5
path.M123 <- '/home/disk/drizzle/wgk/microtia_child_M1M2M3/'
path.cutoff <- paste0(path.M123, 'cutoff_', fc.cutoff, '/')
# path.M12 <- '/home/disk/drizzle/wgk/microtia_child_M1M2/'
# path.cutoff <- paste0(path.M12, 'cutoff_', fc.cutoff, '/')
file.marker.go <- paste0(path.cutoff, 'marker_go.Rdata')
list.marker.go <- readRDS(file.marker.go)

file.go.BP <- paste0(path.cutoff, 'GO_BP.Rdata')
list.go.BP <- readRDS(file.go.BP)
file.go.MF <- paste0(path.cutoff, 'GO_MF.Rdata')
list.go.MF <- readRDS(file.go.MF)


# select GO
df.GO <- data.frame(stringsAsFactors = F)
# Chondral stem cell
GO.BP.CSC.M <- c('response to oxidative stress', 'response to unfolded protein', 
                 'response to tumor necrosis factor', 'RNA splicing',
                 'RNA localization', 'positive regulation of defense response')
GO.BP.CSC.N <- c('ribosome biogenesis', 'response to copper ion', 'oxidative phosphorylation',
                 'extracellular matrix organization', 'skeletal system development',
                 'cell aggregation', 'cellular zinc ion homeostasis')
sel.GO.BP <- c('response to oxidative stress', 
               'response to unfolded protein', 
               'positive regulation of defense response', 
               # 'regulation of inflammatory response',
               'p38MAPK cascade', 'ERK1 and ERK2 cascade',
               # 'regulation of ERK1 and ERK2 cascade', 
               'intrinsic apoptotic signaling pathway', 
               'cell cycle arrest',
               'negative regulation of stem cell differentiation',
               'negative regulation of cell growth',
               'RNA splicing', 'RNA localization', 
               'vascular endothelial growth factor production', 'angiogenesis', 
               'positive regulation of vasculature development',
               'positive regulation of cell migration',
               'negative regulation of cell adhesion',
               'translational initiation', 'ribosome biogenesis', 
               'cartilage condensation', 
               'extracellular matrix organization', 
               # 'skeletal system development',
               'connective tissue development',
               'zinc ion homeostasis')
sel.GO.MF <- c('extracellular matrix structural constituent', 
               'extracellular matrix binding', 'S100 protein binding')

terms <- c("Chondral stem cell_Microtia_increase",
           "Chondral stem cell_Microtia_decrease",
           "Transitional chondrocyte_Microtia_increase",
           "Transitional chondrocyte_Microtia_decrease",
           "Chondrocyte1_Microtia_increase",
           "Chondrocyte1_Microtia_decrease",
           "Chondrocyte2_Microtia_increase",
           "Chondrocyte2_Microtia_decrease")

df.plot <- data.frame(stringsAsFactors = F)
for (term in terms) {
  cell <- strsplit(term, split = '_')[[1]][1]
  status <- strsplit(term, split = '_')[[1]][3]
  sub.BP <- list.go.BP[[term]]
  rownames(sub.BP) <- sub.BP$Description
  sub.BP <- sub.BP[sub.BP$p.adjust < 0.1,]
  sel.BP <- sub.BP[sel.GO.BP, c('Description', 'pvalue', 'geneID')]
  sel.BP$Description <- sel.GO.BP
  sub.MF <- list.go.MF[[term]]
  rownames(sub.MF) <- sub.MF$Description
  sub.MF <- sub.MF[sub.MF$p.adjust < 0.1,]
  sel.MF <- sub.MF[sel.GO.MF, c('Description', 'pvalue', 'geneID')]
  sel.MF$Description <- sel.GO.MF
  sub.plot <- rbind(sel.BP, sel.MF)
  sub.plot$pvalue[is.na(sub.plot$pvalue)] <- 1
  sub.plot$CellType <- rep(cell, nrow(sub.plot))
  if (status == 'increase') {
    sub.plot$Status <- rep('Microtia', nrow(sub.plot))
    sub.plot$Coeff <- rep(1, nrow(sub.plot))
  } else {
    sub.plot$Status <- rep('Normal', nrow(sub.plot))
    sub.plot$Coeff <- rep(-1, nrow(sub.plot))
  }
  df.plot <- rbind(df.plot, sub.plot)
}
df.plot$log10Pval <- -log10(df.plot$pvalue)
df.plot$log10Pval[abs(df.plot$log10Pval) > 10] = 10
df.plot$log10Pval <- df.plot$log10Pval * df.plot$Coeff
df.plot$abs_log10Pval <- abs(df.plot$log10Pval)
df.plot$Description <- factor(df.plot$Description, levels = c(sel.GO.BP, sel.GO.MF))
df.plot$CellType_raw <- df.plot$CellType
df.plot$CellType[df.plot$CellType_raw=='Chondral stem cell']  <- 'CSC'
df.plot$CellType[df.plot$CellType_raw=='Transitional chondrocyte']  <- 'TC'
df.plot$CellType <- factor(df.plot$CellType,levels = c('CSC','TC','Chondrocyte1','Chondrocyte2'))

df.plot$col_name <- paste(df.plot$CellType, df.plot$Status, sep = '_')
mat.plot <- reshape2::dcast(df.plot, Description ~ col_name, value.var = 'log10Pval')
row.names(mat.plot) <- mat.plot$Description
mat.plot$Description <- NULL

ann_colors = list(
  CellType = c(CSC="#33A02C",TC="#B2DF8A",Chondrocyte1="#1F78B4",Chondrocyte2="#A6CEE3"),
  Status = c(Normal = "#637FBF", Microtia = "#6C6C6C")
) 

# col annotation
annotation_col = data.frame(
  CellType = factor(c(rep('CSC', 2),
                      rep('Chondrocyte1', 2),
                      rep('Chondrocyte2', 2),
                      rep('TC', 2)), 
                    levels = c('CSC', 'TC',
                               'Chondrocyte1', 'Chondrocyte2')), 
  Status = factor(rep(c('Microtia', 'Normal'), 4), levels = c('Normal', 'Microtia')),
  row.names = colnames(mat.plot)
)

cols <- c("CSC_Microtia", "TC_Microtia", 
          "Chondrocyte1_Microtia", "Chondrocyte2_Microtia",
          "CSC_Normal", "TC_Normal",
          "Chondrocyte1_Normal", "Chondrocyte2_Normal")
mat.plot <- mat.plot[rev(rownames(mat.plot)), cols]
annotation_col <- annotation_col[cols,]

ann_colors = list(
  CellType = c(CSC="#33A02C",TC="#B2DF8A",Chondrocyte1="#1F78B4",Chondrocyte2="#A6CEE3"),
  Status = c(Normal = "#6C6C6C", Microtia = "#637FBF")
) 

plot.heatmap <- pheatmap::pheatmap(mat.plot,fontsize = 6,fontsize_row = 6,fontsize_col = 6,
                   color = colorRampPalette(c('blue', 'white', 'red'))(100),
                   cluster_rows = F, cluster_cols = F, scale = "none",
                   display_numbers = F,border_color='grey',
                   annotation_col = annotation_col ,annotation_colors = ann_colors,
                   show_rownames = T, show_colnames = F, legend = T, 
                   # fontsize_row = 18, fontsize_col = 15,
                   gaps_col = c(4),lwd=0.5)
ggsave('/home/yzj/JingMA_NEW/res/compMicrotia/Fig4C.pdf',plot.heatmap,width = 13,height = 8,units = 'cm')
