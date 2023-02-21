library(Seurat)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)

pbmc_C <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype_Control_Chond.Rdata')
pbmc_C$Phase <- factor(pbmc_C$Phase,levels = c('Children','Adults'))

### gene associated aging
GA_mtx <- read.csv('publicData/GenAge/human_genes/genage_human.csv')
Aging.db <- GA_mtx[,2]
length(Aging.db)

GO <- read.gmt(gmtfile = 'publicData/GMT/c5.all.v6.2.symbols.gmt')
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

################
## Fig3D(上) 与GenAge/AgeAtlas数据库的fisher.test的热图
################
#### 不分上下调
dnValues_mtx <- readRDS('JingMA_NEW/res/compControl/ChildrenvsAdults/FIG/DEGHeatmap_DNmtx.RDS')
dnValues_mtx <- as.data.frame(dnValues_mtx)
colnames(dnValues_mtx)[2] <- 'C0'

upValues_mtx <- readRDS('JingMA_NEW/res/compControl/ChildrenvsAdults/FIG/DEGHeatmap_UPmtx.RDS')
upValues_mtx <- as.data.frame(upValues_mtx)
colnames(upValues_mtx)[2] <- 'C0'

Values_mtx <- rbind(dnValues_mtx,upValues_mtx)

d <- length(intersect(keys(org.Hs.eg.db, keytype = "SYMBOL"),rownames(pbmc_C)))

upPval <- c()
for(j in 1:length(Aging.lst)){
  Aging.Set <- Aging.lst[[j]]
  upval <- c()
  for(i in 1:length(Values_mtx)){
    upGene<- rownames(Values_mtx)[which(Values_mtx[[i]]==1)]
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
colnames(melted_cormat) <- c('database','celltype','logp')
# Load ggplot2
library(ggplot2)
library(RColorBrewer)
col=brewer.pal(n = 12, name ='Set3')

# Barplot
p.up <- ggplot(data=melted_cormat, mapping=aes(x=celltype,y=logp,fill=database))+
  geom_bar(stat="identity",width=0.8,position='dodge',alpha = 1)+theme_classic()+
  theme(axis.text = element_text(size = 6,colour = 'black'),axis.title = element_text(size = 6,colour = 'black'),
        legend.key.height = unit(0.3,'cm'),legend.key.width = unit(0.3,'cm'),
        legend.title=element_text(size=5),legend.text=element_text(size=6),
        plot.margin = unit(c(0.1,0,0,0.1),'cm'))+
  scale_fill_manual(values=col[6:7])+
  labs(x="", y="-log(adjPvalue)")+geom_hline(aes(yintercept=-log(0.05,10)),colour="#990000", linetype="dashed")
ggsave('JingMA_NEW/res/compControl/ChildrenvsAdults/FIG/Fig3D_GenAge_FC.pdf',p.up,width = 7,height = 4 ,units = 'cm')


#### 对成人来说上调
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
colnames(melted_cormat) <- c('database','celltype','logp')
# Load ggplot2
library(ggplot2)
library(RColorBrewer)
col=brewer.pal(n = 8, name ='Dark2')
# Barplot
p.up <- ggplot(data=melted_cormat, mapping=aes(x=celltype,y=logp,fill=database))+
  geom_bar(stat="identity",width=0.8,position='dodge',alpha = 0.7)+theme_classic()+
  theme(axis.text = element_text(size = 6,colour = 'black'),axis.title = element_text(size = 6,colour = 'black'),
        legend.key.height = unit(0.3,'cm'),legend.key.width = unit(0.3,'cm'),
        legend.title=element_text(size=5),legend.text=element_text(size=6),
        plot.margin = unit(c(0.1,0,0,0.1),'cm'))+
  scale_fill_manual(values=col[1:2])+
  labs(x="", y="-log(adjPvalue)")+geom_hline(aes(yintercept=-log(0.05,10)),colour="#990000", linetype="dashed")
ggsave('JingMA_NEW/res/compControl/ChildrenvsAdults/FIG/Fig3D_upGenAge_FC.pdf',p.up,width = 7,height = 5 ,units = 'cm')

# 
# p.up <- ggplot(data = melted_cormat, aes(x=Var2, y=Var1, fill=value)) + geom_tile(color = "white",size=1)+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
#         axis.title = element_blank(),axis.text  = element_blank(),axis.ticks = element_blank(),
#         legend.position = 'bottom',legend.title  = element_text(size=3),
#         legend.text = element_text(size=3),legend.key.width = unit(0.2,'cm'),
#         legend.key.size = unit(0.2,'cm'),legend.key.height = unit(0.2,'cm'),
#         plot.margin = unit(c(0,0,-0.1,-0.1),units = 'cm'))+
#   scale_fill_gradient2(low = '#EFEFEF', high = '#B15E72') 


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
colnames(melted_cormat) <- c('database','celltype','logp')
# Load ggplot2
library(ggplot2)
library(RColorBrewer)
col=brewer.pal(n = 8, name ='Dark2')
# Barplot
p.dn <- ggplot(data=melted_cormat, mapping=aes(x=celltype,y=logp,fill=database))+
  geom_bar(stat="identity",width=0.8,position='dodge',alpha = 0.7)+theme_classic()+
  theme(axis.text = element_text(size = 6,colour = 'black'),axis.title = element_text(size = 6,colour = 'black'),
        legend.key.height = unit(0.3,'cm'),legend.key.width = unit(0.3,'cm'),
        legend.title=element_text(size=5),legend.text=element_text(size=6),
        plot.margin = unit(c(0.1,0,0,0.1),'cm'))+
  scale_fill_manual(values=col[1:2])+
  labs(x="", y="-log(adjPvalue)")+geom_hline(aes(yintercept=-log(0.05,10)),colour="#990000", linetype="dashed")
ggsave('JingMA_NEW/res/compControl/ChildrenvsAdults/FIG/Fig3D_dnGenAge_FC.pdf',p.dn,width = 7,height = 5 ,units = 'cm')



# p.dn <- ggplot(data = melted_cormat, aes(x=Var2, y=Var1, fill=value)) + geom_tile(color = "white",size=1)+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
#         axis.title = element_blank(),axis.text  = element_blank(),axis.ticks = element_blank(),
#         legend.position = 'bottom',legend.title  = element_text(size=3),
#         legend.text = element_text(size=3),legend.key.width = unit(0.2,'cm'),
#         legend.key.size = unit(0.2,'cm'),legend.key.height = unit(0.2,'cm'),
#         plot.margin = unit(c(0,0,-0.1,-0.1),units = 'cm'))+
#   scale_fill_gradient2(low = "white", high = '#7F99CE') 
# ggsave('JingMA_NEW/res/compControl/ChildrenvsAdults/FIG/Fig3D_dnGenAge_FC.pdf',p.dn,width = 3,height =2 ,units = 'cm')


###
upOL_GAmtx <- Aging.BIG_df[Aging.BIG_df$Symbol %in% upOL,]
dnOL_GAmtx <- Aging.BIG_df[Aging.BIG_df$Symbol %in% dnOL,]
write.xlsx(upOL_GAmtx,'JingMA_NEW/res/compControl/ChildrenvsAdults/DEG/FC1.5_adjP0.05/OL_GAmtx.xlsx',row.names = FALSE,sheetName = 'UP')
write.xlsx(dnOL_GAmtx,'JingMA_NEW/res/compControl/ChildrenvsAdults/DEG/FC1.5_adjP0.05/OL_GAmtx.xlsx',row.names = FALSE,sheetName = 'DN',append = TRUE)



################################################################################
##### 补充材料
################################################################################
library(Seurat)
library(ggplot2)
pbmc <- readRDS('/home/yzj/JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype.Rdata')
pbmc_C <- readRDS('/home/yzj/JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype_Control_Chond.Rdata')
GA_mtx <- read.csv('publicData/GenAge/human_genes/genage_human.csv')
Aging.db <- GA_mtx[,2]
length(Aging.db)

library(clusterProfiler)
GO <- read.gmt(gmtfile = '/home/yzj/publicData/GMT/c5.all.v6.2.symbols.gmt')
Aging.GO <- GO$gene[GO$term=='GO_AGING']
length(Aging.GO)

length(intersect(Aging.db,Aging.GO))

Aging.genes <- union(Aging.db,Aging.GO)
length(Aging.genes)

############
### 参考卵巢衰老FigS4A，markergene在每种细胞类型中的children/Adults的差异，发现差异很大
############
pbmc_C$Phase <- factor(pbmc_C$Phase,levels = c('Children','Adults'))
marker.genes <- c('CDH5','CLDN5','PDGFRB','ACTA2','PTPRC','HLA-DRA','COL1A1','LUM','VCAN','ACAN','COL9A2','CYTL1','ELN','COL2A1','EGR1','HES1')
Idents(pbmc_C) <- pbmc_C$Phase


exp <- as.data.frame(t(as.data.frame(pbmc_C@assays$RNA@data[marker.genes,])))
exp.melt <- reshape2::melt(exp)
exp.melt$celltype <- rep(pbmc_C$celltype.abbr,length(marker.genes))
exp.melt$phase <- rep(pbmc_C$Phase,length(marker.genes))
colnames(exp.melt)[1:2] <- c('gene','expvalue')

df.gene <- data.frame(stringsAsFactors = F)
for (gene in marker.genes) {
  df.sub <- data.frame(expvalue = pbmc_C@assays$RNA@data[gene,],
                       gene = rep(gene, ncol(pbmc_C@assays$RNA@data)),
                       celltype = pbmc_C$celltype.abbr,
                       phase=pbmc_C$Phase)
  df.gene <- rbind(df.gene, df.sub)
}
df.plot <- df.gene
df.plot$gene <- factor(df.gene$gene, levels = marker.genes)
df.plot$celltype <- factor(df.gene$celltype, 
                           levels = c('CSC', 'C', 'SSC', 'SC','IC','PVC','EC'))
color.cell <- c("#A6CEE3" ,"#1F78B4","#CAB2D6","#33A02C","#E31A1C","#FF7F00" ,"#6A3D9A")
plot.vln <- 
  ggplot(data = exp.melt, aes(x = gene, y = expvalue, color = phase, fill = phase)) + 
  geom_violin(trim = T, scale = 'width') + 
  facet_grid( ~celltype) + 
  theme_classic() + coord_flip() +
  stat_summary(fun= mean, geom = "point",
               shape = 23, size = 2, color = "black") + 
  labs(x = 'Gene', y = 'Expression Level') + 
  theme(axis.text.y = element_text(
    size = 13, color = "black", face = 'bold.italic'), 
    axis.text.x = element_text(
      size = 10, color = "black", face = "bold"),
    axis.title = element_text(size = 15, face = 'bold'), 
    strip.text.x = element_text(
      size = 12, color = "black", face = "bold"), 
    legend.position = 'none')
plot.vln


############
### 参考卵巢衰老FigS4B,4C，GeneAge是否在所有细胞类型里高表达
############
Idents(pbmc) <- pbmc$celltype
sort.cells <- c('IC', 'EC', 'PVC', 'SC','SSC', 'CSC', 'C')
color.cell <- c("#E31A1C","#6A3D9A","#FF7F00" ,"#33A02C","#CAB2D6" ,"#A6CEE3","#1F78B4")

## 所有细胞类型
subpbmc <- subset(pbmc,features = Aging.db)
tmp <- AverageExpression(subpbmc, return.seurat = TRUE)

exp <- as.matrix(tmp@assays$RNA@scale.data)
p_all <- pheatmap::pheatmap(exp,cluster_rows = T,cluster_cols = F,show_rownames = F)

## 软骨细胞类型
subpbmc <- subset(pbmc_C,features = Aging.db)
tmp <- AverageExpression(subpbmc, return.seurat = TRUE)

exp <- as.matrix(tmp@assays$RNA@scale.data)
p_c <- pheatmap::pheatmap(exp,cluster_rows = T,cluster_cols = F,show_rownames = F)


############
### 参考卵巢衰老FigS4D，GeneAge是否在所有细胞类型里高表达
############
subpbmc <- subset(pbmc_C,features = Aging.db)
subpbmc@meta.data$group <- paste(subpbmc$celltype,subpbmc$Phase,sep='_')
subpbmc$group <- factor(subpbmc$group,levels = c(paste('CSC',c('Children','Adults'),sep='_'),paste('C0',c('Children','Adults'),sep='_'),
                                                 paste('C1',c('Children','Adults'),sep='_'),paste('C2',c('Children','Adults'),sep='_')))
Idents(subpbmc) <- subpbmc$group
tmp <- AverageExpression(subpbmc, return.seurat = TRUE)
exp <- as.matrix(tmp@assays$RNA@scale.data)


p_group <- pheatmap::pheatmap(exp,cluster_rows = T,cluster_cols = F,show_rownames = F)


############
### GA在各种细胞类型里与差异表达基因的overlap
############

## 在所有细胞里
MK.all <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/Markers_celltype.RDS')
print(names(MK.all))

OL.lst <- list()
OL.ratio <- c()
for(i in 1:length(MK.all)){
  CT=names(MK.all)[i]
  mk.df <- MK.all[[CT]]
  mk.gene <- rownames(mk.df)[mk.df$avg_logFC > log(1.5,2) & mk.df$p_val_adj < 0.05]
  ol.gene <- intersect(Aging.genes,mk.gene)
  OL.lst[[CT]] <- ol.gene
  OL.ratio[i] <- round(length(ol.gene)/length(mk.gene),3)
}

names(OL.ratio) <- names(OL.lst)
library(ggplot2)
barplot(OL.ratio)

