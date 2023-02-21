setwd('/home/zy/my_git/bioinformatics/wgk')
.libPaths('/home/yzj/R/x86_64-pc-linux-gnu-library/4.0')
library(Seurat)
library(ggplot2)
library(SCENIC)
require("RColorBrewer")


path.data <- '/home/disk/drizzle/wgk/data/AllSample_2_merge/'
path.lineage <- paste0(path.data, 'chon_lineage/')
file.chon <- paste0(path.lineage, 'seurat_celltype.Rdata')
seurat.chon <- readRDS(file.chon)
seurat.child <- subset(seurat.chon, subset = batch %in% c('C4', 'C6', 'M1', 'M2', 'M3'))
# seurat.child <- NormalizeData(seurat.child)
# seurat.child <- FindVariableFeatures(seurat.child, nfeatures = 3000)
# highvar.genes <- VariableFeatures(seurat.child)
# 
# path.M123 <- '/home/disk/drizzle/wgk/microtia_chon_child_M1M2M3/'
# path.change <- paste0(path.M123, 'chon_lineage_1/')
# if (!file.exists(path.change)) {
#     dir.create(path.change)
# }
# 
# 
# seurat.child.N <- subset(seurat.child, subset = type == 'Normal')
# seurat.child.M <- subset(seurat.child, subset = type == 'Microtia')
# 
# library(foreach)
# library(doParallel)
# registerDoParallel(cores = 10)
# corr.PC.exp <- function(mat.gene, PC2, gene) {
#     vec.exp <- mat.gene[gene,]
#     corr <- cor(PC2, vec.exp, method = 'spearman')
#     return(data.frame(gene = gene, corr = corr))
# }
# 
# # Normal
# Harmony2_N <- seurat.child.N@reductions$harmony@cell.embeddings[, 'harmony_2']
# mat.gene_N <- seurat.child.N@assays$RNA@data
# df.corr_N <- 
#     foreach(gene = highvar.genes, .combine = rbind) %dopar% 
#     corr.PC.exp(mat.gene_N, Harmony2_N, gene)
# rownames(df.corr_N) <- highvar.genes
# 
# # Microtia
# Harmony2_M <- seurat.child.M@reductions$harmony@cell.embeddings[, 'harmony_2']
# mat.gene_M <- seurat.child.M@assays$RNA@data
# df.corr_M <- 
#     foreach(gene = highvar.genes, .combine = rbind) %dopar% 
#     corr.PC.exp(mat.gene_M, Harmony2_M, gene)
# rownames(df.corr_M) <- highvar.genes
# 
# df.corr <- merge(df.corr_N, df.corr_M, by = 'gene')
# df.corr[is.na(df.corr)] <- 0
# df.corr$diff_M_N <- df.corr$corr.y - df.corr$corr.x

###############
###############
# plot single gene
Harmony2 <- seurat.child@reductions$harmony@cell.embeddings[, 'harmony_2']
##mat.gene <- seurat.child@assays$RNA@data
# AUC
regulonAUC <- readRDS(file='/home/yzj/JingMA_NEW/res/SCENIC_main/int/3.4_regulonAUC.Rds')
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)), colnames(mat.gene)]
mat.auc <- as.matrix(regulonAUC@assays@data@listData$AUC)
# mat.gene.TF <- rbind(as.matrix(mat.gene), mat.auc)
##df.pc.gene <- data.frame(t(rbind(as.matrix(mat.gene), mat.auc)), check.names = F)
df.pc.gene <- data.frame(t( mat.auc), check.names = F)
df.pc.gene$Harmony2 <- Harmony2
df.pc.gene$celltype <- seurat.child$celltype
df.pc.gene$status <- seurat.child$type
df.pc.gene <- df.pc.gene[order(Harmony2, decreasing = F),]
df.pc.gene$idx <- 1:nrow(df.pc.gene)
colors <- c("#EE9572","#B2DF8A" ,"#A6CEE3","#9999FF")
names(colors) <- c('CSC', 'TC', 'C1', 'C2')

gene <- 'SOX9'
gene <- 'SOX5 (218g)'
df.plot <- df.pc.gene[, c('idx', 'Harmony2', 'celltype', 'status', gene)]
names(df.plot) <- c('idx', 'Harmony2', 'celltype', 'status', 'gene')
p.gene <-
    ggplot(data = df.plot, aes(x = idx, 
                               #linetype = status, 
                               y = gene)) + 
    geom_point(aes(color = celltype), size = 0.3) + 
    scale_color_manual(labels = c('CSC', 'TC', 'C1', 'C2'),
                       values = colors) + 
    # xlim(-30, 10) + 
    geom_smooth(color = '#696969') +
    # geom_smooth(color = '#696969', formula = y~poly(x, 3), method = lm) + 
    labs(x = '', y = '') + 
    theme(panel.background=element_rect(fill='transparent', color='black',size = 1),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = 'none') +
    theme(panel.background=element_rect(fill='transparent', color='black',size = 1)) +
    annotate('text', label = gene, x = 22000, y = max(df.plot$gene), 
             hjust = 1, vjust = 1, size = 7)

    
# time series
# four time
library(maSigPro)
mat.count <- seurat.child@assays$RNA@counts[highvar.genes,]
Time <- as.character(seurat.child$celltype)
Time[Time == 'CSC'] <- 1
Time[Time == 'TC'] <- 2
Time[Time == 'C1'] <- 3
Time[Time == 'C2'] <- 4

# 20 points
library(maSigPro)
mat.count <- seurat.child@assays$RNA@counts
Time <- as.character(seurat.child$celltype)
Harmony2 <- seurat.child@reductions$harmony@cell.embeddings[, 'harmony_2']
quantile.H2 <- quantile(Harmony2, probs = seq(0, 1, 0.05))
for (i in 1:20) {
    Time[Harmony2 >= quantile.H2[i] & Harmony2 <= quantile.H2[i+1]] <- i
}

# time series analysis
Time <- as.numeric(Time)
Replicates <- as.numeric(as.factor(seurat.child$batch))
Normal <- rep(0, length(seurat.child$type))
Normal[seurat.child$type == 'Normal'] <- 1
Microtia <- rep(0, length(seurat.child$type))
Microtia[seurat.child$type == 'Microtia'] <- 1
edesign <- cbind(Time,Replicates,Normal,Microtia)
# rownames(edesign) = paste("Cell",c(1:nrow(edesign)), sep= "_")
edesign <- data.frame(Time,Replicates,Normal,Microtia,
                      row.names = colnames(mat.count))
d = make.design.matrix(edesign, degree = 1)

# filter genes
path.M123 <- '/home/disk/drizzle/wgk/microtia_chon_child_M1M2M3/'
file.gsea.marker <- paste0(path.M123, 'marker_GSEA.Rdata')
list.marker.gsea <- readRDS(file.gsea.marker)
diff.genes <- c()
for (cell in names(list.marker.gsea)) {
    sub.marker <- list.marker.gsea[[cell]]
    sub.marker.1 <- sub.marker[sub.marker$p_val_adj < 0.01 & 
                                   abs(sub.marker$avg_logFC) > 0.25,]
    diff.genes <- c(diff.genes, rownames(sub.marker.1))
}
diff.genes <- unique(diff.genes)


high.exp <- function(mat.gene, gene) {
    vec.single <- mat.gene[gene,]
    agge.group <- aggregate(vec.single, by = list(cond_times), FUN = mean)
    if (max(agge.group$x) > 0.25) {
        return(data.frame(gene = gene, high.exp = 1))
    } else {
        return(data.frame(gene = gene, high.exp = 0))
    }
}

library(foreach)
library(doParallel)
registerDoParallel(cores = 10)
cond_times <- paste(seurat.child$type, Time, sep = '_')
df.highexp <- 
    foreach(gene = diff.genes, .combine = rbind) %dopar% high.exp(mat.gene, gene)
genes.highexp <- df.highexp[df.highexp$high.exp == 1, 'gene']

# maSigPro
fit = p.vector(mat.count[genes.highexp,], d, Q=0.01, 
               MT.adjust = "bonferroni",counts = T)
file.fit <- paste0(path.change, 'fit.Rdata')
saveRDS(fit, file.fit)

fit <- readRDS(file.fit)
tstep = T.fit(fit,step.method="backward")
file.tstep <- paste0(path.change, 'tstep.Rdata')
saveRDS(tstep, file.tstep)

tstep <- readRDS(file.tstep)

sigs = get.siggenes(tstep, rsq = 0.25, vars = "groups")

sigs.new <- sigs$sig.genes$MicrotiavsNormal
new.genes <- intersect(rownames(sigs.new$sig.profiles), names(df.pc.gene))
length(new.genes)

# XY genes
df.genes.v19 <- read.delim(
    '/home/yzj/publicData/annotation/hg19/gencode_v19_gene_pos.txt',
    header = F)
XY.genes <- df.genes.v19$V1[df.genes.v19$V2 %in% c('chrY', 'chrX')]
new.genes <- setdiff(new.genes, XY.genes)
length(new.genes)


sigs.new$sig.profiles <- sigs.new$sig.profiles[new.genes,]
sigs.new$coefficients <- sigs.new$coefficients[new.genes,]
sigs.new$group.coeffs <- sigs.new$group.coeffs[new.genes,]
sigs.new$sig.pvalues <- sigs.new$sig.pvalues[new.genes,]

# df.rsq <- tstep$sol
df.pval <- sigs.new$sig.pvalues
df.coff <- sigs.new$coefficients

pdf(file = paste(path.change,"cluster_gene_1.pdf",sep=''), width = 20, height = 20);
res.cluster <- see.genes(
    sigs$sig.genes$MicrotiavsNormal,
    show.fit = T,
    dis =d$dis,
    cluster.method="hclust" ,
    cluster.data = 1,
    k = 12,
    # k.mclust = T,
    ylim = c(0, 20),
    newX11 = F)
dev.off()


View(merge(res.cluster$cut, df.pval, by = 'row.names'))
View(merge(res.cluster$cut, df.coff, by = 'row.names'))


# plot by clust
sel.genes <- intersect(names(res.cluster$cut), names(df.pc.gene))
df.pc.gene.sel <- df.pc.gene[, c('idx', 'Harmony2', 'celltype', 'status', sel.genes)]
clusters <- names(table(res.cluster$cut))
for (cluster in clusters) {
    sub.gene <- intersect(names(res.cluster$cut[res.cluster$cut == cluster]), sel.genes)
    if (length(sub.gene) < 1) {next()}
    sub.df.gene <- df.pc.gene.sel[,sub.gene]
    sub.df <- data.frame(idx = df.pc.gene.sel$idx, status = df.pc.gene.sel$status, 
                         celltype = df.pc.gene.sel$celltype, 
                         mean.exp = rowMeans(sub.df.gene))
    p.mean <-
        ggplot(data = sub.df, aes(x = idx, y = mean.exp, linetype = status)) + 
        geom_point(aes(color = celltype), size = 0.3) + 
        scale_color_manual(labels = c('CSC', 'TC', 'C1', 'C2'),
                           values = colors) + 
        # xlim(-30, 10) + 
        geom_smooth(color = '#696969') + 
        labs(x = '', y = '') + 
        # theme(panel.background=element_rect(fill='transparent', color='black',size = 1),
        #       axis.text = element_blank(),
        #       axis.ticks = element_blank(),
        #       legend.position = 'none') +
        theme(panel.background=element_rect(fill='transparent', color='black',size = 1)) +
        annotate('text', label = paste0('Cluster_', cluster), 
                 x = 22000, y = max(sub.df$mean.exp), 
                 hjust = 1, vjust = 1, size = 7)
    ggsave(p.mean, path = path.change, 
           filename = paste0('Cluster_', cluster, '.png'),
           height = 6, width = 9, units = 'cm')
}


# 









