library(Seurat)
library(ggplot2)
library(SCENIC)
require("RColorBrewer")

seurat.chon <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype_Control_Chond.Rdata')
mat.gene <- seurat.chon@assays$RNA@data

# AUC
regulonAUC <- readRDS(file='/local/yzj/JingMA_NEW/res/SCENIC_main/int/3.4_regulonAUC.Rds')
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)), colnames(mat.gene)]
mat.auc <- as.matrix(regulonAUC@assays@data@listData$AUC)
TF <- apply(as.matrix(rownames(mat.auc)),1, function(x) unlist(strsplit(x,'_'))[1] )
TF <- apply(as.matrix(TF),1, function(x) unlist(strsplit(x,' '))[1] )

# # mat.gene.TF <- rbind(as.matrix(mat.gene), mat.auc)
# df.pc.gene <- data.frame(t(rbind(as.matrix(mat.gene), mat.auc)), check.names = F)
# #df.pc.gene <- data.frame(t( mat.auc), check.names = F)
# df.pc.gene$Harmony2 <- Harmony2
# df.pc.gene$celltype <- seurat.chon$celltype
# df.pc.gene$status <- seurat.chon$Phase
# df.pc.gene <- df.pc.gene[order(Harmony2, decreasing = F),]
# df.pc.gene$idx <- 1:nrow(df.pc.gene)
# saveRDS(df.pc.gene,'JingMA_NEW/res/maSigPro/control/df.pc.gene.RDS')
# df.pc.gene <- readRDS('JingMA_NEW/res/maSigPro/RDS/df.pc.gene.RDS')

###############
###############
# time series
library(maSigPro)

# 20 points
mat.count <- seurat.chon@assays$RNA@counts
Time <- as.character(seurat.chon$celltype)
Harmony2 <- seurat.chon@reductions$harmony@cell.embeddings[, 'harmony_2']
quantile.H2 <- quantile(Harmony2, probs = seq(0, 1, 0.05))
for (i in 1:20) {
  Time[Harmony2 >= quantile.H2[i] & Harmony2 <= quantile.H2[i+1]] <- i
}

# time series analysis
Time <- as.numeric(Time)
Replicates <- as.numeric(as.factor(seurat.chon$batch))
Children <- rep(0, length(seurat.chon$Phase))
Children[seurat.chon$Phase == 'Children'] <- 1
Adults <- rep(0, length(seurat.chon$Phase))
Adults[seurat.chon$Phase == 'Adults'] <- 1
edesign <- cbind(Time,Replicates,Children,Adults)
# rownames(edesign) = paste("Cell",c(1:nrow(edesign)), sep= "_")
edesign <- data.frame(Time,Replicates,Children,Adults,
                      row.names = colnames(mat.count))
d = make.design.matrix(edesign, degree = 1)

# filter genes
list.marker.gsea <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/DEGs_inChond_inChildrenAdults.RDS')
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
registerDoParallel(cores = 40)
mat.gene <- seurat.chon@assays$RNA@data
cond_times <- paste(seurat.chon$Phase, Time, sep = '_')
df.highexp <- 
  foreach(gene = diff.genes, .combine = rbind) %dopar% high.exp(mat.gene, gene)
genes.highexp <- df.highexp[df.highexp$high.exp == 1, 'gene']

# maSigPro
path.change='/local/yzj/JingMA_NEW/res/maSigPro/control/'
fit = p.vector(mat.count[genes.highexp,], d, Q=0.01, 
               MT.adjust = "bonferroni",counts = T)
file.fit <- paste0(path.change, 'fit.Rdata')
saveRDS(fit, file.fit)

# fit <- readRDS(file.fit)
tstep = T.fit(fit,step.method="backward")
file.tstep <- paste0(path.change, 'tstep.Rdata')
saveRDS(tstep, file.tstep)

# tstep <- readRDS(file.tstep)
sigs = get.siggenes(tstep, rsq = 0.25, vars = "groups")

sigs.new <- sigs$sig.genes$AdultsvsChildren
new.genes <- intersect(rownames(sigs.new$sig.profiles), TF)
length(new.genes)

# XY genes
df.genes.v19 <- read.delim(
  '/local/yzj/publicData/annotation/hg19/gencode_v19_gene_pos.txt',
  header = F)
XY.genes <- df.genes.v19$V1[df.genes.v19$V2 %in% c('chrY', 'chrX')]
new.genes <- setdiff(new.genes, XY.genes)
length(new.genes)


pdf(file = paste(path.change,"cluster_gene_1.pdf",sep=''), width = 20, height = 20);
res.cluster <- see.genes(
  sigs$sig.genes$AdultsvsChildren,
  show.fit = T,
  dis =d$dis,
  cluster.method="hclust" ,
  cluster.data = 1,
  k = 12,
  # k.mclust = T,
  ylim = c(0, 20),
  newX11 = F)
dev.off()
saveRDS(res.cluster,paste(path.change,'cluster.Rdata',sep=''))

# # df.rsq <- tstep$sol
# df.pval <- sigs.new$sig.pvalues
# df.coff <- sigs.new$coefficients
# 
# View(merge(res.cluster$cut, df.pval, by = 'row.names'))
# View(merge(res.cluster$cut, df.coff, by = 'row.names'))
# 
# 
# # plot by clust
# sel.genes <- intersect(names(res.cluster$cut), names(df.pc.gene))
# df.pc.gene.sel <- df.pc.gene[, c('idx', 'Harmony2', 'celltype', 'status', sel.genes)]
# clusters <- names(table(res.cluster$cut))
# colors <- c("#EE9572","#B2DF8A" ,"#A6CEE3","#9999FF")
# names(colors) <- c('CSC', 'C0', 'C1', 'C2')
# 
# for (cluster in clusters) {
#   sub.gene <- intersect(names(res.cluster$cut[res.cluster$cut == cluster]), sel.genes)
#   if (length(sub.gene) < 1) {next()}
#   sub.df.gene <- df.pc.gene.sel[,sub.gene]
#   sub.df <- data.frame(idx = df.pc.gene.sel$idx, status = df.pc.gene.sel$status, 
#                        celltype = df.pc.gene.sel$celltype, 
#                        mean.exp = rowMeans(sub.df.gene))
#   p.mean <-
#     ggplot(data = sub.df, aes(x = idx, y = mean.exp, linetype = status)) + 
#     geom_point(aes(color = celltype), size = 0.3) + 
#     scale_color_manual(labels = c('CSC', 'C0', 'C1', 'C2'),
#                        values = colors) + 
#     # xlim(-30, 10) + 
#     geom_smooth(color = '#696969') + 
#     labs(x = '', y = '') + 
#     # theme(panel.background=element_rect(fill='transparent', color='black',size = 1),
#     #       axis.text = element_blank(),
#     #       axis.ticks = element_blank(),
#     #       legend.position = 'none') +
#     theme(panel.background=element_rect(fill='transparent', color='black',size = 1)) +
#     annotate('text', label = paste0('Cluster_', cluster), 
#              x = 22000, y = max(sub.df$mean.exp), 
#              hjust = 1, vjust = 1, size = 7)
#   ggsave(p.mean, path = path.change, 
#          filename = paste0('Cluster_', cluster, '.png'),
#          height = 6, width = 9, units = 'cm')
# }
# 
# # 

