#########
## Fig3 比较儿童和成年人年龄段(C4/C6 vs C1/C2/C3/C5)
#########
library(tibble)
library(dplyr)
library(Seurat)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)

pbmc_C <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype_Control_Chond.Rdata')
pbmc_C$batch= factor(pbmc_C$batch,levels = c('C6','C4','C2','C3','C1','C5'))
pbmc_C$Phase <- factor(pbmc_C$Phase,levels = c('Children','Adults'))

########################
## Fig3A 细胞比例
########################
phylog_batch <- pbmc_C@meta.data[,c('batch',"celltype")]
phylog_batch <- table(phylog_batch$batch,phylog_batch[,"celltype"])
phylog_batch <- data.frame(phylog_batch)
colnames(phylog_batch) <- c('SampleID','CellType','Freq')
print(phylog_batch$CellType)

Color <- rev(c("#EE9572","#B2DF8A" ,"#A6CEE3","#9999FF"))
p1 <- ggplot(phylog_batch,aes(x=SampleID,y=Freq,fill=fct_rev(CellType)))+
  geom_col(position = "fill", width = 0.7)+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5,size = 8,colour = 'black'))+
  theme(axis.text = element_text(size = 8,colour = "black"),
        axis.line = element_line(size=0.5, colour = "black"),
        axis.title.y = element_text(size=8),
        plot.title = element_text(hjust = 0.5))+
  labs(x='',y='% of chondrocytes',title = 'Per sample')+
  scale_fill_manual(values = (Color))+
  guides(fill = guide_legend(reverse=TRUE))+labs(fill = "Cell type")
p1


phylog_phase <- pbmc_C@meta.data[,c('Phase',"celltype")]
phylog_phase <- table(phylog_phase$Phase,phylog_phase[,"celltype"])
phylog_phase <- data.frame(phylog_phase)
colnames(phylog_phase) <- c('Phase','CellType','Freq')
phylog_phase$CellType <- factor(phylog_phase$CellType)
phylog_phase$Phase <- factor(phylog_phase$Phase,levels = c('Children','Adults'))


Color <- c("#EE9572","#B2DF8A" ,"#A6CEE3","#9999FF")
p2 <- ggplot(phylog_phase,aes(x=Phase,y=Freq,fill=fct_rev(CellType)))+
  geom_col(position = "fill", width = 0.9)+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5,size = 8,colour = 'black'))+
  theme(axis.text = element_text(size = 8,colour = "black"),
        axis.line = element_line(size=0.5, colour = "black"),
        axis.title = element_text(size=8),legend.title = element_text(size = 6,color = 'black'),
        legend.text = element_text(size = 6,color = 'black'),legend.key.size = unit(0.3,'cm'),
        legend.key.width  = unit(0.3,'cm'),legend.key.height = unit(0.3,'cm'),)+
  labs(x='',y='',title = 'Mean')+theme(legend.position="right",plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(values = rev(Color))+
  guides(fill = guide_legend(reverse=TRUE))+labs(fill='Cell type')
p2

library(ggpubr)
p <- ggarrange(p1, p2,ncol = 2, nrow = 1)
ggsave('/home/yzj/JingMA_NEW/res/Harmony/ALL/FIG/Fig3A_Barplot_PropChond.pdf',p,
       width = 10,height = 6,units = 'cm')


########################
## Fig3B DEG画heatmap,参照卵巢衰老的文章
########################
save_pheatmap_pdf <- function(x, filename, width=4, height=5) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
DEG.lst <- readRDS('/home/yzj/JingMA_NEW/res/Harmony/ALL/RDS/DEGs_inChond_inChildrenAdults.RDS')
print(names(DEG.lst))

data <- DEG.lst[['CSC']]
up_CSC <- rownames(data)[data$avg_logFC > log(1.5) & data$p_val_adj < 0.05]
dn_CSC <- rownames(data)[data$avg_logFC < -log(1.5) & data$p_val_adj < 0.05]

data <- DEG.lst[["TC"]]
up_C0 <- rownames(data)[data$avg_logFC > log(1.5) & data$p_val_adj < 0.05]
dn_C0 <- rownames(data)[data$avg_logFC < -log(1.5) & data$p_val_adj < 0.05]

data <- DEG.lst[["C1"]]
up_C1 <- rownames(data)[data$avg_logFC > log(1.5) & data$p_val_adj < 0.05]
dn_C1 <- rownames(data)[data$avg_logFC < -log(1.5) & data$p_val_adj < 0.05]

data <- DEG.lst[["C2"]]
up_C2 <- rownames(data)[data$avg_logFC > log(1.5) & data$p_val_adj < 0.05]
dn_C2 <- rownames(data)[data$avg_logFC < -log(1.5) & data$p_val_adj < 0.05]

get_values <- function(sigCSC,sigC0,sigC1,sigC2){
  sigGene <- unique(c(sigCSC,sigC0,sigC1,sigC2))
  values <- matrix(c(rep(0,4*length(sigGene))),ncol = 4,dimnames = list(sigGene,c("CSPC", "EC","IC","LC")))
  for(i in 1:length(sigGene)){
    g=sigGene[i]
    if(g %in% sigCSC){values[i,1] <-1};
    if(g %in% sigC0){values[i,2] <-1};
    if(g %in% sigC1){values[i,3] <-1};
    if(g %in% sigC2){values[i,4] <-1};
  }
  values_sum <- apply(values, 1, sum)
  values <- values[order(values_sum,decreasing = T),]
  return(values)
}


## 对成人来说上调矩阵
dnValues_mtx <- get_values(dn_CSC,dn_C0,dn_C1,dn_C2)
dn_sum <- apply(dnValues_mtx,1,sum)
dn_df <- dnValues_mtx[-(which(dn_sum>1)),]

annotation_col = data.frame(CellType = factor(c("CSC", "C0","C1","C2")))
rownames(annotation_col) <- colnames(dnValues_mtx)

annotation_row = data.frame(GeneClass = factor(rep(c("Common", "CSC", "C0","C1","C2"), c(length(which(dn_sum>1)), length(which(dn_df[,1]==1)), 
                                                     length(which(dn_df[,2]==1)),length(which(dn_df[,3]==1)),length(which(dn_df[,4]==1))))))
rownames(annotation_row) = rownames(dnValues_mtx)

ann_colors = list( CellType = c(CSC="#EE9572",C0="#B2DF8A",C1="#A6CEE3",C2="#9999FF"),
                   GeneClass = c(Common='grey',CSC="#EE9572",C0="#B2DF8A",C1="#A6CEE3",C2="#9999FF"))

p_DN <- pheatmap(dnValues_mtx,cluster_rows = F,cluster_cols = F,color =  colorRampPalette(c("#EFEFEF", "white","#B15E72"))(100),
                 border=FALSE,show_rownames = F,legend=F,angle_col='45',
                 annotation_row = annotation_row,annotation_colors = ann_colors,annotation_legend = FALSE)
save_pheatmap_pdf(p_DN,'JingMA_NEW/res/compControl/ChildrenvsAdults/FIG/Fig3B_DEGHeatmap_DN.pdf',height = 4,width = 2)
saveRDS(dnValues_mtx,'JingMA_NEW/res/compControl/ChildrenvsAdults/FIG/DEGHeatmap_DNmtx.RDS')


## 对成人来说，下调矩阵
upValues_mtx <- get_values(up_CSC,up_C0,up_C1,up_C2)
up_sum <- apply(upValues_mtx,1,sum)
up_df <- upValues_mtx[-(which(up_sum>1)),]

annotation_col = data.frame(CellType = factor(c("CSPC", "EC","IC","LC")))
rownames(annotation_col) <- colnames(upValues_mtx)

annotation_row = data.frame(GeneClass = factor(rep(c("Common", "CSPC", "EC","IC","LC"), c(length(which(up_sum>1)), length(which(up_df[,1]==1)), 
                                                                                          length(which(up_df[,2]==1)),length(which(up_df[,3]==1)),length(which(up_df[,4]==1))))))
rownames(annotation_row) = rownames(upValues_mtx)

ann_colors = list( CellType = c(CSPC="#EE9572",EC="#B2DF8A",IC="#A6CEE3",LC="#9999FF"),
                   GeneClass = c(Common='grey',CSPC="#EE9572",EC="#B2DF8A",IC="#A6CEE3",LC="#9999FF"))

p_UP <- pheatmap(upValues_mtx,cluster_rows = F,cluster_cols = F,color =  colorRampPalette(c("#EFEFEF", "white", "#7F99CE"))(100),
                 border=FALSE,show_rownames = F,angle_col='45',
                 annotation_row = annotation_row,annotation_colors = ann_colors,legend=F,annotation_legend = FALSE)
save_pheatmap_pdf(p_UP,'JingMA_NEW/res/compControl/ChildrenvsAdults/FIG/Fig3B_DEGHeatmap_UP.pdf',height = 4,width = 2)
saveRDS(upValues_mtx,'JingMA_NEW/res/compControl/ChildrenvsAdults/FIG/DEGHeatmap_UPmtx.RDS')

### print数字
up_sum <- apply(upValues_mtx,1,sum)
length(which(up_sum==4))
print(length(which(up_sum>1))/nrow(upValues_mtx))

up_df <- upValues_mtx[-(which(up_sum>1)),]
print(length(which(up_df[,1]==1))/nrow(upValues_mtx))
print(length(which(up_df[,2]==1))/nrow(upValues_mtx))
print(length(which(up_df[,3]==1))/nrow(upValues_mtx))
print(length(which(up_df[,4]==1))/nrow(upValues_mtx))

dn_sum <- apply(dnValues_mtx,1,sum)
length(which(dn_sum==4))
print(length(which(dn_sum>1))/nrow(dnValues_mtx))

dn_df <- dnValues_mtx[-(which(dn_sum>1)),]
print(length(which(dn_df[,1]==1))/nrow(dnValues_mtx))
print(length(which(dn_df[,2]==1))/nrow(dnValues_mtx))
print(length(which(dn_df[,3]==1))/nrow(dnValues_mtx))
print(length(which(dn_df[,4]==1))/nrow(dnValues_mtx))
###



#########
### Fig3C 挑选term组合画图
#########
pickterm_DN <- c("extracellular matrix structural constituent","extracellular matrix binding","collagen fibril organization",
  "cartilage development","chondrocyte differentiation","chondrocyte morphogenesis","NAD metabolic process")

pickterm_UP <- c("aging","response to oxidative stress","reactive oxygen species metabolic process","reactive oxygen species biosynthetic process",
                 "cell death in response to oxidative stress", "response to interleukin-6",'regulation of inflammatory response',"intrinsic apoptotic signaling pathway",
                 'cellular response to tumor necrosis factor',"p38MAPK cascade","ERK1 and ERK2 cascade",
                 "DNA damage response, signal transduction by p53 class mediator","negative regulation of stem cell differentiation","cellular senescence")

library(xlsx)
CSC_DN_BP <- read.xlsx('/home/yzj/JingMA_NEW/res/compControl/ChildrenvsAdults/ClusterPro/FC1.5_adjP0.05/CSC_BP.xlsx',sheetName = 'DN')
CSC_DN_BP <- CSC_DN_BP[CSC_DN_BP$p.adjust < 0.1,]
index <- which(CSC_DN_BP$Description %in% pickterm_DN)
pickCSC_DN_BP <- CSC_DN_BP[index,]

CSC_DN_MF <- read.xlsx('/home/yzj/JingMA_NEW/res/compControl/ChildrenvsAdults/ClusterPro/FC1.5_adjP0.05/CSC_MF_all.xlsx',sheetName = 'DN')
CSC_DN_MF <- CSC_DN_MF[CSC_DN_MF$p.adjust < 0.1,]
index <- which(CSC_DN_MF$Description %in% pickterm_DN)
pickCSC_DN_MF <- CSC_DN_MF[index,]

pickCSC_DN <- rbind(pickCSC_DN_BP,pickCSC_DN_MF)

CSC_UP_BP <- read.xlsx('/home/yzj/JingMA_NEW/res/compControl/ChildrenvsAdults/ClusterPro/FC1.5_adjP0.05/CSC_BP_all.xlsx',sheetName = 'UP')
CSC_UP_BP <- CSC_UP_BP[CSC_UP_BP$p.adjust < 0.1,]
print(CSC_UP_BP$Description)
index <- c(which(CSC_UP_BP$Description %in% pickterm_UP))
pickCSC_UP_BP <- CSC_UP_BP[index,]

pickCSC_UP <- pickCSC_UP_BP

####
Trans_DN_BP <- read.xlsx('/home/yzj/JingMA_NEW/res/compControl/ChildrenvsAdults/ClusterPro/FC1.5_adjP0.05/TC_BP.xlsx',sheetName = 'DN')
Trans_DN_BP <- Trans_DN_BP[Trans_DN_BP$p.adjust < 0.1,]
index <- which(Trans_DN_BP$Description %in% pickterm_DN)
pickTrans_DN_BP <- Trans_DN_BP[index,]

Trans_DN_MF <- read.xlsx('/home/yzj/JingMA_NEW/res/compControl/ChildrenvsAdults/ClusterPro/FC1.5_adjP0.05/TC_MF_all.xlsx',sheetName = 'DN')
Trans_DN_MF <- Trans_DN_MF[Trans_DN_MF$p.adjust < 0.1,]
index <- which(Trans_DN_MF$Description %in% pickterm_DN)
pickTrans_DN_MF <- Trans_DN_MF[index,]

pickTrans_DN <- rbind(pickTrans_DN_BP,pickTrans_DN_MF)

Trans_UP_BP <- read.xlsx('/home/yzj/JingMA_NEW/res/compControl/ChildrenvsAdults/ClusterPro/FC1.5_adjP0.05/TC_BP_all.xlsx',sheetName = 'UP')
Trans_UP_BP <- Trans_UP_BP[Trans_UP_BP$p.adjust < 0.1,]
index <- c(which(Trans_UP_BP$Description %in% pickterm_UP))
pickTrans_UP_BP <- Trans_UP_BP[index,]
pickTrans_UP <- pickTrans_UP_BP

####
Chond1_DN_BP <- read.xlsx('/home/yzj/JingMA_NEW/res/compControl/ChildrenvsAdults/ClusterPro/FC1.5_adjP0.05/C1_BP.xlsx',sheetName = 'DN')
Chond1_DN_BP <- Chond1_DN_BP[Chond1_DN_BP$p.adjust < 0.1,]
index <- which(Chond1_DN_BP$Description %in% pickterm_DN)
pickChond1_DN_BP <- Chond1_DN_BP[index,]

Chond1_DN_MF <- read.xlsx('/home/yzj/JingMA_NEW/res/compControl/ChildrenvsAdults/ClusterPro/FC1.5_adjP0.05/C1_MF_all.xlsx',sheetName = 'DN')
Chond1_DN_MF <- Chond1_DN_MF[Chond1_DN_MF$p.adjust < 0.1,]
index <- which(Chond1_DN_MF$Description %in% pickterm_DN)
pickChond1_DN_MF <- Chond1_DN_MF[index,]

pickChond1_DN <- rbind(pickChond1_DN_BP,pickChond1_DN_MF)

Chond1_UP_BP <- read.xlsx('/home/yzj/JingMA_NEW/res/compControl/ChildrenvsAdults/ClusterPro/FC1.5_adjP0.05/C1_BP_all.xlsx',sheetName = 'UP')
Chond1_UP_BP <- Chond1_UP_BP[Chond1_UP_BP$p.adjust < 0.1,]
index <- which(Chond1_UP_BP$Description %in% pickterm_UP)
pickChond1_UP_BP <- Chond1_UP_BP[index,]

pickChond1_UP <- pickChond1_UP_BP

####
Chond2_DN_BP <- read.xlsx('/home/yzj/JingMA_NEW/res/compControl/ChildrenvsAdults/ClusterPro/FC1.5_adjP0.05/C2_BP.xlsx',sheetName = 'DN')
Chond2_DN_BP <- Chond2_DN_BP[Chond2_DN_BP$p.adjust < 0.1,]
index <-  which(Chond2_DN_BP$Description %in% pickterm_DN)
pickChond2_DN_BP <- Chond2_DN_BP[index,]

Chond2_DN_MF <- read.xlsx('/home/yzj/JingMA_NEW/res/compControl/ChildrenvsAdults/ClusterPro/FC1.5_adjP0.05/C2_MF_all.xlsx',sheetName = 'DN')
Chond2_DN_MF <- Chond2_DN_MF[Chond2_DN_MF$p.adjust < 0.1,]
index <- which(Chond2_DN_MF$Description %in% pickterm_DN)
pickChond2_DN_MF <- Chond2_DN_MF[index,]

pickChond2_DN <- rbind(pickChond2_DN_BP,pickChond2_DN_MF)

Chond2_UP_BP <- read.xlsx('/home/yzj/JingMA_NEW/res/compControl/ChildrenvsAdults/ClusterPro/FC1.5_adjP0.05/C2_BP_all.xlsx',sheetName = 'UP')
Chond2_UP_BP <- Chond2_UP_BP[Chond2_UP_BP$p.adjust < 0.1,]
index <- which(Chond2_UP_BP$Description %in% pickterm_UP)
pickChond2_UP_BP <- Chond2_UP_BP[index,]

pickChond2_UP <- pickChond2_UP_BP

bar_DN <- rbind(pickCSC_DN,pickTrans_DN,pickChond1_DN,pickChond2_DN)
bar_DN$CellType <- c(rep('CSC',nrow(pickCSC_DN)),rep('C0',nrow(pickTrans_DN)),
                     rep('C1',nrow(pickChond1_DN)),rep('C2',nrow(pickChond2_DN)))
bar_DN$CellType <- factor(bar_DN$CellType,levels = c('CSC','C0','C1','C2'))
bar_DN$Group <- 'Children'
bar_DN$log10Pval <- -log(bar_DN$p.adjust,10)

bar_UP <- rbind(pickCSC_UP,pickTrans_UP,pickChond1_UP,pickChond2_UP)
bar_UP$CellType <- c(rep('CSC',nrow(pickCSC_UP)),rep('C0',nrow(pickTrans_UP)),
                     rep('C1',nrow(pickChond1_UP)),rep('C2',nrow(pickChond2_UP)))
bar_UP$CellType <- factor(bar_UP$CellType,levels = c('CSC','C0','C1','C2'))
bar_UP$Group <- 'Adults'
bar_UP$log10Pval <- log(bar_UP$p.adjust,10)
bar_df <- rbind(bar_DN,bar_UP)

levels_DN=rev(pickterm_DN)
setdiff(unique(bar_DN$Description),levels_DN)

levels_UP=rev(pickterm_UP)
setdiff(unique(bar_UP$Description),levels_UP)

bar_df$Description <- factor(bar_df$Description,levels = rev(c(levels_UP,levels_DN)))
bar_df$Count <- as.numeric(bar_df$Count)


library(reshape2)
mat.plot <- bar_df[,c('Description','CellType','Group','log10Pval')]
mat.plot <- dcast(mat.plot,Description~CellType+Group)
mat.plot[is.na(mat.plot)] <- 0
rownames(mat.plot) <- mat.plot$Description
mat.plot <- mat.plot[,-1]
colNames <- c('CSC_Children','C0_Children','C1_Children','C2_Children','CSC_Adults','C0_Adults','C1_Adults','C2_Adults')
mat.plot <- dplyr::select(mat.plot,colNames)

# col annotation
annotation_col = data.frame(
  CellType = factor(c(rep('CSC', 2),rep('C0', 2),rep('C1', 2),rep('C2', 2)), levels = c('CSC', 'C0','C1', 'C2')), 
  Phase = factor(rep(c('Children', 'Adults'), 4), levels = c('Children', 'Adults')),row.names = colNames)

annotation_col = data.frame(
  CellType = factor(rep(c('CSC', 'C0','C1', 'C2'), 2), levels = c('CSC', 'C0','C1', 'C2')), 
  Phase = factor(rep(c('Children', 'Adults'), each=4), levels = c('Children', 'Adults')),row.names = colNames)

ann_colors = list(
  CellType = c(CSC="#EE9572",C0="#B2DF8A",C1="#A6CEE3",C2="#9999FF"),
  Phase = c(Children = "#6C6C6C", Adults = "#DC143C")) 

bk <- c(seq(-8,-0.1,by=0.01),seq(0,8,by=0.01))
plot.heatmap <- pheatmap::pheatmap(mat.plot,fontsize = 6,fontsize_row = 6,fontsize_col = 6,
                   cluster_rows = F, cluster_cols = F, scale = "none",
                   display_numbers = F,border_color='grey',
                   annotation_col = annotation_col ,annotation_colors = ann_colors,
                   show_rownames = T, show_colnames = F, legend = TRUE, 
                   gaps_col = c(4),lwd=0.5,
                   color = c(colorRampPalette(colors = c("red","white"))(length(bk)/2),colorRampPalette(colors = c("white","blue"))(length(bk)/2)),
                   legend_breaks=seq(-8,8,2),
                   breaks=bk)

ggsave('/home/yzj/JingMA_NEW/res/compControl/ChildrenvsAdults/FIG/Fig3C_pickHeatmap.pdf',
       plot.heatmap,width = 13,height = 8,units = 'cm')


################
## Fig3D 与GenAge/AgeAtlas数据库的fisher.test
################
### gene associated aging
GA_mtx <- read.csv('publicData/GenAge/human_genes/genage_human.csv')
Aging.db <- GA_mtx[,2]
length(Aging.db)

GO <- read.gmt(gmtfile = '/home/yzj/publicData/GMT/c5.all.v6.2.symbols.gmt')
Aging.GO <- GO$gene[GO$term=='GO_AGING']
length(Aging.GO)

library(xlsx)
Aging.BIG_df <- read.xlsx('publicData/GeneAtlas/aging_list_V2.0.xlsx',sheetIndex = 1)
Aging.BIG_df <- as.data.frame(Aging.BIG_df)
Aging.BIG <- Aging.BIG_df$Symbol
length(Aging.BIG)

#Aging.genes <- union(union(Aging.db,Aging.GO),Aging.BIG)
Aging.genes <- union(Aging.db,Aging.BIG)
length(Aging.genes)

Aging.lst <- list(GenAge=Aging.db,AgingAtlas=Aging.BIG)

#### 成人上调
dnValues_mtx <- readRDS('JingMA_NEW/res/compControl/ChildrenvsAdults/FIG/DEGHeatmap_DNmtx.RDS')
dnValues_mtx <- as.data.frame(dnValues_mtx)
colnames(dnValues_mtx)[2] <- 'C0'

d <- length(intersect(keys(org.Hs.eg.db, keytype = "SYMBOL"),rownames(pbmc_C)))
upPval <- c()
for(j in 1:length(Aging.lst)){
  Aging.Set <- Aging.lst[[j]]
  upval <- c()
  for(i in 1:length(dnValues_mtx)){
    upGene<- rownames(dnValues_mtx)[which(dnValues_mtx[[i]]==1)]
    #a: DEG in aging.genes b: aging.genes, c: DEG, d: bg gene
    a <- length(intersect(Aging.Set,upGene));b <- length(Aging.Set);c <- length(upGene)
    p <- fisher.test(matrix(c(a,b,c,d), nrow=2), alternative="greater")$p.value
    #upval[i] <- -log(p,10)
    upval[i] <- p
  }
  upPval <- rbind(upPval,upval)
}

colnames(upPval) <- c('CSC','EC','IC','LC');rownames(upPval) <- names(Aging.lst)
print(upPval)
upPval <- -log(upPval,10)

melted_cormat <- reshape2::melt(upPval)
p.up <- ggplot(data = melted_cormat, aes(x=Var2, y=Var1, fill=value)) + geom_tile(color = "white",size=1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.title = element_blank(),axis.text  = element_blank(),axis.ticks = element_blank(),
        legend.position = 'bottom',legend.title  = element_text(size=3),
        legend.text = element_text(size=3),legend.key.width = unit(0.2,'cm'),
        legend.key.size = unit(0.2,'cm'),legend.key.height = unit(0.2,'cm'),
        plot.margin = unit(c(0,0,-0.1,-0.1),units = 'cm'))+
  scale_fill_gradient2(low = '#EFEFEF', high = '#B15E72') 
ggsave('JingMA_NEW/res/compControl/ChildrenvsAdults/FIG/Fig3D_upGenAge_FC.pdf',p.up,width = 3,height =1.7 ,units = 'cm')


#### 成人下调
upValues_mtx <- readRDS('JingMA_NEW/res/compControl/ChildrenvsAdults/FIG/DEGHeatmap_UPmtx.RDS')
upValues_mtx <- as.data.frame(upValues_mtx)
colnames(upValues_mtx)[2] <- 'C0'

dnPval <- c()
for(j in 1:length(Aging.lst)){
  Aging.Set <- Aging.lst[[j]]
  dnpval <- c()
  for(i in 1:length(upValues_mtx)){
    dnGene<- rownames(upValues_mtx)[which(upValues_mtx[[i]]==1)]
    a <- length(intersect(Aging.Set,dnGene));b <- length(Aging.Set);c <- length(dnGene)
    p <- fisher.test(matrix(c(a,b,c,d), nrow=2), alternative="greater")$p.value
    #dnpval[i] <- -log(p,10)
    dnpval[i] <- p
  }
  dnPval <- rbind(dnPval,dnpval)
}

colnames(dnPval) <- c('CSC','EC','IC','LC');rownames(dnPval) <- names(Aging.lst)
print(dnPval)
dnPval <- -log(dnPval,10)

melted_cormat <- reshape2::melt(dnPval)
p.dn <- ggplot(data = melted_cormat, aes(x=Var2, y=Var1, fill=value)) + geom_tile(color = "white",size=1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.title = element_blank(),axis.text  = element_blank(),axis.ticks = element_blank(),
        legend.position = 'bottom',legend.title  = element_text(size=3),
        legend.text = element_text(size=3),legend.key.width = unit(0.2,'cm'),
        legend.key.size = unit(0.2,'cm'),legend.key.height = unit(0.2,'cm'),
        plot.margin = unit(c(0,0,-0.1,-0.1),units = 'cm'))+
  scale_fill_gradient2(low = "white", high = '#7F99CE') 
ggsave('JingMA_NEW/res/compControl/ChildrenvsAdults/FIG/Fig3D_dnGenAge_FC.pdf',p.dn,width = 3,height =2 ,units = 'cm')


###
upOL_GAmtx <- Aging.BIG_df[Aging.BIG_df$Symbol %in% upOL,]
dnOL_GAmtx <- Aging.BIG_df[Aging.BIG_df$Symbol %in% dnOL,]
write.xlsx(upOL_GAmtx,'JingMA_NEW/res/compControl/ChildrenvsAdults/DEG/FC1.5_adjP0.05/OL_GAmtx.xlsx',row.names = FALSE,sheetName = 'UP')
write.xlsx(dnOL_GAmtx,'JingMA_NEW/res/compControl/ChildrenvsAdults/DEG/FC1.5_adjP0.05/OL_GAmtx.xlsx',row.names = FALSE,sheetName = 'DN',append = TRUE)


