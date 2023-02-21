## 在生科院服务器上产生.loom文件(JingMA_ORI/rawData/bam)
library(Seurat)
library(velocyto.R)
library(SeuratWrappers)

pbmc <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype.Rdata')
Idents(pbmc) <- pbmc$celltype.abbr
sample_cells <- sample(colnames(pbmc),10000)
pbmc_sample <- subset(pbmc,cells = sample_cells)
cells <- gsub('_',':',colnames(pbmc_sample))
cells <- gsub('-1','x',cells)
rm(pbmc);rm(pbmc)

print("Step1!")
print(Sys.time())
batch='All'
input <- ReadVelocity(file = paste('/home/yzj/JingMA_ORI/res/RNAvelocity/',batch,'/',batch,'.loom',sep = ''))
RV <- as.Seurat(x = input)
RV <- subset(RV,cells=cells)
RV <- SCTransform(object = RV, assay = "spliced")
saveRDS(RV,paste('/home/yzj/JingMA_NEW/res/RNAvelocity/',batch,'/RV_SCT.RDS',sep=''))

print("Step2!")
print(Sys.time())
RV <- readRDS(paste('/home/yzj/JingMA_NEW/res/RNAvelocity/',batch,'/RV_SCT.RDS',sep=''))
RV <- RunPCA(object = RV, seed.use=123, npcs=50, features = VariableFeatures(object = RV), ndims.print=1,nfeatures.print=1)
RV <- RunUMAP(RV, dims = 1:50, seed.use = 123,n.components=2,n.neighbors = 10,min.dist = 0.5)

RV@reductions$umap@cell.embeddings[,1] <- pbmc_sample@reductions$umap@cell.embeddings[,1]
RV@reductions$umap@cell.embeddings[,2] <- pbmc_sample@reductions$umap@cell.embeddings[,2]
RV@meta.data$Cluster <- pbmc_sample@meta.data$celltype.abbr
Idents(RV) <- RV@meta.data$Cluster

print("Step3!")
print(Sys.time())
RV <- RunVelocity(object = RV, deltaT = 1, kCells = 25, fit.quantile = 0.1,n.cores=15)
saveRDS(RV,paste('/home/yzj/JingMA_NEW/res/RNAvelocity/',batch,'/RV.RDS',sep=''))

batch='All'
print("Step4!")
print(Sys.time())
RV <- readRDS(paste('JingMA_NEW/res/RNAvelocity/',batch,'/RV.RDS',sep=''))

cells <- gsub(':','_',colnames(RV))
cells <- gsub('x','-1',cells)
subpbmc <- subset(pbmc,cells = cells)
print(all(colnames(subpbmc)==cells))
RV@meta.data$Cluster <- subpbmc@meta.data$celltype.abbr
Idents(RV) <- RV@meta.data$Cluster

Color <- c("#A6CEE3","#1F78B4","#CAB2D6" ,"#33A02C","#E31A1C" ,"#FF7F00","#6A3D9A")
ident.colors=Color

names(x = ident.colors) <- levels(x = RV)
cell.colors <- ident.colors[Idents(object = RV)]
names(x = cell.colors) <- colnames(x = RV)

for(n in c(150)){
  N=n
  for(s in c('sqrt')){
    scale=s
    pdf(paste('JingMA_NEW/res/RNAvelocity/',batch,'/RV_',n,'_',scale,'.pdf',sep=''))
    p <- show.velocity.on.embedding.cor(emb = Embeddings(object = RV, reduction = "umap"),
                                        vel = Tool(object = RV, slot = "RunVelocity"),
                                        n = N, scale = scale, cell.colors = ac(x = cell.colors, alpha = 0.5),
                                        cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5,
                                        grid.n = 40, arrow.lwd = 1, n.cores=70,
                                        do.par = FALSE, cell.border.alpha = 0.1)
    dev.off()
  }
}

print('DONE!')
print(Sys.time())

