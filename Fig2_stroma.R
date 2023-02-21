library(Seurat)
library(slingshot)
library(ggplot2)

############
## 1. for chond lineage
############
pbmc_chond <- readRDS('/home/yzj/JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype_Chond.Rdata')
Idents(pbmc_chond) <- pbmc_chond$celltype

# sample_cells <- sample(colnames(pbmc_chond),5000)
# pbmc_sample <- subset(pbmc_chond,cells = sample_cells)

rd=Embeddings(pbmc_chond, "umap")
cl1=Idents(pbmc_chond)

crv1 <- slingshot(rd, cl1,start.clus='CSC',end.clus='C2')
Color <- c("#EE9572","#B2DF8A" ,"#A6CEE3","#9999FF")
names(Color) <- c('CSC', 'C0', 'C1','C2')

pdf('/home/yzj/JingMA_NEW/res/slingshot/Chond/Chond_slingshot.pdf')
plot(rd, col = Color[cl1], asp = 1, pch = 16)
lines(crv1, lwd = 4, col = 'black',lty = 4)
dev.off()
saveRDS(crv1,'/home/yzj/JingMA_NEW/res/slingshot/Chond/Chond_slingshot.RDS')


############
## 2. for stromal lineage
############
pbmc_stroma <- readRDS('/home/disk/drizzle/wgk/data/AllSample_2_merge/stroma_lineage/seurat_celltype.Rdata')
Idents(pbmc_stroma) <- pbmc_stroma$celltype
Color <- c("#BC80BD", "#80B1D3", "#F4A460")
names(Color) <- c('SSC', 'SC1', 'SC2')
pbmc_stroma$celltype <- factor(pbmc_stroma$celltype,levels =c('SSC', 'SC1', 'SC2') )

DimPlot(pbmc_stroma)
rd=Embeddings(pbmc_stroma, "umap")

p <- DimPlot(pbmc_stroma, group.by='celltype', label=T ,label.size = 8,pt.size = 1,cols = Color)+
  theme(axis.text = element_text(size=15,face = 'bold',colour = 'black'),
        axis.title = element_text(size=15,face = 'bold',colour = 'black'),
        panel.background=element_rect(fill='transparent', color='black',size = 1),
        legend.key=element_rect(fill='transparent', color='transparent'),
        legend.text = element_text(size=15,face = 'bold',colour = 'black'))

data=as.data.frame(pbmc_stroma@reductions$umap@cell.embeddings)
data$CellType <- pbmc_stroma$celltype

plot.umap <- ggplot(data = data,aes(x=UMAP_1,y=UMAP_2,color=CellType))+
  geom_point(size=0.1)+
  theme_bw()+
  theme(axis.title = element_text(size = 7,colour = 'black'),
        panel.background=element_rect(fill='transparent', color='black',size = 1),
        legend.key=element_rect(fill='transparent', color='transparent'),
        legend.text = element_text(size = 7,colour = 'black'),
        legend.title = element_text(size = 7,colour = 'black'),
        legend.key.height = unit(0.05,'cm'),
        legend.key.size = unit(0.01, "cm"),legend.key.width = unit(0.05,"cm"),
        plot.margin = unit(c(0.1,0,0.1,0.1),'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(), 
        axis.ticks = element_blank()) + 
  guides(colour = guide_legend(override.aes = list(size=5)))+
  scale_color_manual(values = Color)
ggsave('/home/yzj/JingMA_NEW/res/slingshot/Stroma/Fig2E_UMAP_Stroma.pdf',plot.umap,width =6 ,height = 4,units = 'cm')

cl1=Idents(pbmc_stroma)

slingshot <- slingshot(rd, cl1,start.clus='SSC',end.clus=c('SC1','SC2'))
saveRDS(slingshot,'/home/yzj/JingMA_NEW/res/slingshot/Stroma/Stroma_slingshot.RDS')

# plot(rd, col = Color[cl1], asp = 1, pch = 16)
# lines(slingshot, lwd = 3, col = 'black',lty = 6)

###############
### Fig2E. Slingshot 图
###############

library(ggplot2)
slingshot <- readRDS('/home/yzj/JingMA_NEW/res/slingshot/Stroma/Stroma_slingshot.RDS')
data <- data.frame(UMAP_1=rd[,1],UMAP_2=rd[,2],CellType=pbmc_stroma$celltype)
levels(data$CellType)

PS <- slingshot::slingPseudotime(slingshot)
PS <- as.data.frame(PS)

PS_time <- c()
for(i in 1:nrow(PS)){
  if(length(which(is.na(PS[i,])))==0){
    PS_time <-c(PS_time,as.numeric(mean(PS[i,1],PS[i,2])))
  }else{
    PS_time <- c(PS_time,as.numeric(PS[i,][which(!is.na(PS[i,]))]))
  }
}

data$Pseudotime <- PS_time

### arrow location
curve1 <- as.data.frame(slingshot@curves$curve1$s)
curve2 <- as.data.frame(slingshot@curves$curve2$s)
tail(table(curve1$UMAP_2))
index <- order(curve1$UMAP_1)[c(1,239)]
arrow_c1_x1 <- curve1$UMAP_1[index[2]]
arrow_c1_x2 <- curve1$UMAP_1[index[1]]
arrow_c1_y1 <- curve1$UMAP_2[index[2]]
arrow_c1_y2 <- curve1$UMAP_2[index[1]]

index=order(curve2$UMAP_2)[1:2]
arrow_c2_x1 <- curve2$UMAP_1[index[2]]
arrow_c2_x2 <- curve2$UMAP_1[index[1]]
arrow_c2_y1 <- curve2$UMAP_2[index[2]]
arrow_c2_y2 <- curve2$UMAP_2[index[1]]


slingshot.plot <- ggplot(data = data,aes(x=UMAP_1,y=UMAP_2,color=CellType))+
  geom_point(size=0.1)+
  theme_bw()+
  theme(axis.text = element_text(size=5,color='black'),
        axis.title = element_text(size = 5,colour = 'black'),
        panel.background=element_rect(fill='transparent', color='black',size = 0.3),
        legend.key=element_rect(fill='transparent', color='transparent'),
        legend.text = element_text(size = 5,colour = 'black'),
        legend.title = element_text(size = 5,colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_color_manual(values = Color)+
  geom_point(data = curve1, col = 'black',size=0.05)+
  geom_point(data = curve2, col = 'black',size=0.05)+
  geom_segment(x = arrow_c1_x1, y = arrow_c1_y1, xend = arrow_c1_x2, yend = arrow_c1_y2,
               arrow = arrow(length = unit(0.05, "npc")),color = "black")+
  geom_segment(x = arrow_c2_x1, y = arrow_c2_y1, xend = arrow_c2_x2, yend = arrow_c2_y2,
               arrow = arrow(length = unit(0.05, "npc")),color = "black")
slingshot.plot
ggsave('JingMA_NEW/res/slingshot/Stroma/Fig2E_Stroma_slingshot.pdf',slingshot.plot,width = 6,height = 4,units = 'cm')

###############
### Fig2F. Slingshot 图
###############
ps.plot <- ggplot(data = data,aes(x=UMAP_1,y=UMAP_2,colour=Pseudotime))+
  geom_point(size=0.0000001)+
  theme_bw()+
  theme(axis.text = element_text(size=5,color='black'),
        axis.title = element_text(size = 5,colour = 'black'),
        axis.ticks = element_line(size = 0.3,colour = 'black'),
        panel.background=element_rect(fill='transparent', color='black',size = 0.5),
        legend.key=element_rect(fill='transparent', color='transparent'),
        legend.text = element_text(size = 3,colour = 'black'),
        legend.title = element_text(size = 3,colour = 'black'),
        legend.key.size = unit(0.2, "cm"),legend.key.width = unit(0.2,"cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_colour_viridis_c(option = "D")+
  geom_point(data = curve1, col = 'black',size=0.0005)+
  geom_point(data = curve2, col = 'black',size=0.0005)+
  geom_segment(x = arrow_c1_x1, y = arrow_c1_y1, xend = arrow_c1_x2, yend = arrow_c1_y2,
               arrow = arrow(length = unit(0.05, "npc")),color = "black")+
  geom_segment(x = arrow_c2_x1, y = arrow_c2_y1, xend = arrow_c2_x2, yend = arrow_c2_y2,
               arrow = arrow(length = unit(0.05, "npc")),color = "black")
ps.plot
ggsave('JingMA_NEW/res/slingshot/Stroma/Fig2F_Stroma_PS.pdf',ps.plot,width = 6,height = 4,units = 'cm')

saveRDS(data,'JingMA_NEW/res/slingshot/Stroma/Stroma.RDS')

###############
### Fig2G. gene heatmap
###############
library(monocle)
library(pheatmap)
library(colorRamps)
Binner <- function(cds_object,cells_subset){
  df <- data.frame(pData(cds_object[,cells_subset]))
  df <- df[,c("Pseudotime", "State")]
  df <- df[order(df$Pseudotime, decreasing = F),]
  len <- length(df$Pseudotime)
  bin<-round(len/100)
  State <- c()
  value <- c()
  for(i in 0:99){
    if(i < 99){
      start <- 1+(bin*i);stop <- bin+(bin*i)
      value <- median(as.numeric(as.vector(df$State[c(start:stop)])))
      State <- c(State, value)
    }
    else{
      State <- c(State, value)
    }
  }
  return(as.data.frame(State))
}

Color <- c("#BC80BD", "#80B1D3", "#F4A460")
cds <- readRDS('/home/yzj/JingMA_NEW/res/Monocle/Stroma/monocle.RDS')
data <- readRDS('JingMA_NEW/res/slingshot/Stroma/Stroma.RDS')
data <- data[rownames(data) %in% colnames(cds),]
MK.lst <- readRDS('/home/yzj/JingMA_NEW/res/Harmony/ALL/RDS/Markers_celltype_Stroma.RDS')
MK.SSC <- MK.lst$SSC;MK.SC1 <- MK.lst$SC1;MK.SC2 <- MK.lst$SC2

sigMK.SSC <- rownames(MK.SSC)[MK.SSC$avg_logFC > 0.25 & MK.SSC$p_val_adj < 0.05]
sigMK.SC1 <- rownames(MK.SC1)[MK.SC1$avg_logFC > 0.25 & MK.SC1$p_val_adj < 0.05]
sigMK.SC2 <- rownames(MK.SC2)[MK.SC2$avg_logFC > 0.25 & MK.SC2$p_val_adj < 0.05]
print(length(sigMK.SC2))

gene.SC1 <- c(sigMK.SSC[1:350],sigMK.SC1[c(1:161,174,185,189,207,211,237,243,246,248)])
gene.SC2 <- c(sigMK.SSC[1:350],sigMK.SC2[1:170])

SSC.gene <- c('COL1A1','COL1A2','HES1','EGR1')
SC1.gene <- c("VEGRA",'EDNRB',"CXCL5",'FGF2','CCL20','CXCL1','CXCL3','CXCL2',"EGFR","FGF7")
SC2.gene <- c('MMP10',"MMP3",'MMP1','IL6','BMP2','IL23A')
pickSC1.gene <- c(SSC.gene,SC1.gene)
pickSC2.gene <- c(SSC.gene,SC2.gene)

## 修改Pseudotime
Pseudotime <- data$Pseudotime
names(Pseudotime) <- rownames(data)
pData(cds)['Pseudotime'] <- Pseudotime[colnames(cds)]
pData(cds) <- na.omit(pData(cds))

## 修改State
pData <- pData(cds)
pData$State <- NA
pData$Cluster <- as.character(pData$Cluster)
pData$State[pData$Cluster=='SSC'] <- 0
pData$State[pData$Cluster=='SC1'] <- 1
pData$State[pData$Cluster=='SC2'] <- 2
pData(cds)['State'] <- pData$State

SSC.cell <- rownames(data)[data$CellType=='SSC']
SC1.cell <- rownames(data)[data$CellType=='SC1']
SC2.cell <- rownames(data)[data$CellType=='SC2']

cds.SC1 <- cds[gene.SC1,c(SSC.cell,SC1.cell)]
cds.SC2 <- cds[gene.SC2,c(SSC.cell,SC2.cell)]


### SC1 pheatmap: 细胞顺序反过来
newdata.SC1 <- data.frame(Pseudotime = seq(min(pData(cds.SC1)$Pseudotime),
                                       max(pData(cds.SC1)$Pseudotime), length.out = 100))
                      
mtx.SC1 <- genSmoothCurves(cds.SC1, new_data = newdata.SC1)
print(dim(mtx.SC1))
mtx.SC1 <- log10(mtx.SC1+1)

mtx.SC1 = mtx.SC1[!apply(mtx.SC1, 1, sd) == 0, ]
mtx.SC1 = Matrix::t(scale(Matrix::t(mtx.SC1), center = TRUE))
mtx.SC1 = mtx.SC1[is.na(row.names(mtx.SC1)) == FALSE, ]
mtx.SC1[is.nan(mtx.SC1)] = 0
mtx.SC1[mtx.SC1 > 2] = 2
mtx.SC1[mtx.SC1 < -2] = -2

bks <- seq(-2.1, 2.1, by = 0.1)
hmcols <- blue2green2red(length(bks) - 1)
mtx.SC1 <- mtx.SC1[,rev(1:ncol(mtx.SC1))] #细胞顺序反过来!!

bin <- Binner(cds.SC1,colnames(cds.SC1))
anno_col <- data.frame(CellType=rep('SSC',nrow(bin)))
anno_col$CellType[bin=='1'] <- 'SC1'
ann_colors = list( CellType = c(SSC=Color[1],SC1=Color[2]))

labels_row = rep("",nrow(mtx.SC1))
for(g in pickSC1.gene){
  if(g %in% rownames(mtx.SC1)){
    labels_row[which(rownames(mtx.SC1)==g)] <- g
  }
}

plot.SC1 <- pheatmap(mtx.SC1, useRaster = T, cluster_cols = FALSE, fontsize=3,legend=FALSE,annotation_legend=FALSE,
              cluster_rows = FALSE, show_rownames = TRUE, show_colnames = F,  
              silent = TRUE, filename = NA, annotation_col = anno_col,annotation_colors = ann_colors,
              breaks = bks, border_color = NA, color = hmcols,labels_row = labels_row)
ggsave('JingMA_NEW/res/compCT/Stroma/FC1.5_Qvalue0.05/Fig2G_pickSC1_Gene.pdf',plot.SC1,width = 3,height = 7,units = 'cm')

### SC2 pheatmap
#plot.SC2 <- plot_pseudotime_heatmap(cds.SC2,return_heatmap = T,cluster_rows = FALSE,scale_max = 2,scale_min = -2)
newdata.SC2 <- data.frame(Pseudotime = seq(min(pData(cds.SC2)$Pseudotime), 
                                           max(pData(cds.SC2)$Pseudotime), length.out = 100))
mtx.SC2 <- genSmoothCurves(cds.SC2, new_data = newdata.SC2)
mtx.SC2 <- log10(mtx.SC2+1)

mtx.SC2 = mtx.SC2[!apply(mtx.SC2, 1, sd) == 0, ]
mtx.SC2 = Matrix::t(scale(Matrix::t(mtx.SC2), center = TRUE))
mtx.SC2 = mtx.SC2[is.na(row.names(mtx.SC2)) == FALSE, ]
mtx.SC2[is.nan(mtx.SC2)] = 0
mtx.SC2[mtx.SC2 > 2] = 2
mtx.SC2[mtx.SC2 < -2] = -2

bks <- seq(-2.1, 2.1, by = 0.1)
hmcols <- blue2green2red(length(bks) - 1)

bin <- Binner(cds.SC2,colnames(cds.SC2))
anno_col <- data.frame(CellType=rep('SSC',nrow(bin)))
anno_col$CellType[bin=='2'] <- 'SC2'
ann_colors = list( CellType = c(SSC=Color[1],SC2=Color[3]))

labels_row = rep("",nrow(mtx.SC2))
for(g in pickSC2.gene){
  if(g %in% rownames(mtx.SC2)){
    labels_row[which(rownames(mtx.SC2)==g)] <- g
  }
}
plot.SC2 <- pheatmap(mtx.SC2, useRaster = T, cluster_cols = FALSE, fontsize=3,legend=FALSE,annotation_legend=FALSE,
                     cluster_rows = FALSE, show_rownames = TRUE, show_colnames = F,  
                     silent = TRUE, filename = NA, annotation_col = anno_col,annotation_colors = ann_colors,
                     breaks = bks, border_color = NA, color = hmcols,labels_row =labels_row )

#cowplot::plot_grid(plot.SC1$gtable, plot.SC2$gtable, ncol= 2)
ggsave('JingMA_NEW/res/compCT/Stroma/FC1.5_Qvalue0.05/Fig2G_pickSC2_Gene.pdf',plot.SC2,width = 3,height = 7,units = 'cm')



###############
### Fig2H. Enrichment Dot plot
###############
####
## pick term 
####
library(ggplot2)
list.go.BP <- readRDS('/home/disk/drizzle/wgk/data/AllSample_2_merge/stroma_lineage/cutoff_0.4_0.05/GO_BP.Rdata')

list.sel.GO <- list()
list.sel.GO$SSC <- c('transforming growth factor beta receptor signaling pathway', 
                     'extracellular matrix organization', 
                     'connective tissue development',
                     'BMP signaling pathway', 'cell cycle arrest')
list.sel.GO$SC2 <- c('extracellular matrix organization', 
                     'extracellular matrix disassembly',
                     'regulation of receptor signaling pathway via JAK-STAT',
                     'positive regulation of cell migration')
list.sel.GO$SC1 <- c('extracellular matrix organization', 
                     'regulation of inflammatory response', 
                     'negative regulation of cell migration',
                     'regulation of angiogenesis', 
                     'cell chemotaxis',
                     'regulation of smooth muscle cell proliferation')
# barplot
sort.cells <- c('SSC', 'SC1', 'SC2')
colors <- c("#BC80BD", "#80B1D3", "#F4A460")
df.plot <- data.frame()
i = 0
for (cell in sort.cells) {
  i = i + 1
  sub.go <- list.go.BP[[cell]]
  sel.go.term <- list.sel.GO[[cell]]
  sel.go <- sub.go[sub.go$Description %in% sel.go.term, 
                   c('Description', 'pvalue')]
  sel.go$log10Pval <- -log10(sel.go$pvalue)
  sel.go$celltype <- rep(cell, nrow(sel.go))
  # sel.go$Description <- factor(sel.go$Description, levels = rev(sel.go.term))
  df.plot <- rbind(df.plot, sel.go)
}
col_name <- paste(df.plot$celltype, df.plot$Description, sep = '_')
df.plot$col_name <- factor(col_name, levels = rev(col_name))
df.plot$celltype <- factor(df.plot$celltype, levels = sort.cells)

print(summary(df.plot$log10Pval))
p <- ggplot(df.plot, aes(x = celltype, y = col_name, 
                         color = celltype, size = log10Pval)) + 
  geom_point(fill = 'cornsilk') +
  scale_color_manual(breaks = c('SSC', 'SC1', 'SC2'),
                     values = c("#BC80BD", "#80B1D3", "#F4A460")) + 
  scale_size_continuous(range = c(2,8),breaks = c(4,5,6)) +
  scale_y_discrete(breaks = col_name,
                   labels = c(list.sel.GO$SSC, list.sel.GO$SC1, list.sel.GO$SC2)) +
  labs(x = '', title = 'Enriched GO terms', color = 'Cell type',
       size = expression(paste("-log"[10], "(adj", italic("P"), "-value)"))) +
  theme(panel.background = element_rect(color = 'black',fill = 'transparent'), 
        panel.grid.major = element_line(colour = 'gray', size = 0.2, linetype = 5),
        axis.title = element_text(size = 6, color = 'black'), 
        axis.text = element_text(size = 6,  color = 'black'), 
        plot.margin = unit(c(0.2,0.2,-0.3,0),'cm'),
        legend.text = element_text(size = 6, color = 'black'),legend.title = element_text(size = 6, color = 'black'),
        legend.key = element_blank(),legend.key.height = unit(0.01,'cm'),legend.key.width  = unit(0.01,'cm'),legend.key.size = unit(0.01,'cm'),
        legend.position = "bottom",legend.box = "vertical") + 
  guides(colour = guide_legend(override.aes = list(size=2)))
p
ggsave('JingMA_NEW/res/compCT/Stroma/FC1.5_Qvalue0.05/Fig2H_pickStroma_GO.pdf',p, height = 8, width = 8)



#############
## 补充材料
##############
library(Seurat)
library(ggplot2)
pbmc_stroma <- readRDS('/home/disk/drizzle/wgk/data/AllSample_2_merge/stroma_lineage/seurat_celltype.Rdata')

pdf('JingMA_NEW/res/Harmony/ALL/FIG/sFig_strom_lineage_umap.pdf',height = 3, width = 3.5)
DimPlot(pbmc_stroma, group.by = "RNA_snn_res.0.15", label = T, reduction = 'umap', pt.size = 0.35)
dev.off()

pdf('JingMA_NEW/res/Harmony/ALL/FIG/sFig_strom_lineage_type.pdf',height = 3, width = 4)
DimPlot(pbmc_stroma, group.by = "type", label = F, reduction = 'umap', pt.size = 0.01)
dev.off()
