library(Seurat)
pbmc <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype.Rdata')
batches <- unique(pbmc$batch)
Idents(pbmc) <- pbmc$batch

geneinfo <- read.table('JingMA_NEW/res/Harmony/ALL/RDS/DataAva/features.tsv')
gene2ID <- geneinfo[,1]
names(gene2ID) <- geneinfo[,2]

get_output10X <- function(pbmc,batch){
  subpbmc <- subset(pbmc,idents=batch)
  exp=GetAssayData(object = subpbmc, assay = "RNA", slot = "counts") 
  exp=as.data.frame(exp)
  
  ##
  write.table(data.frame(gene2ID[rownames(exp)],rownames(exp)),
              file = paste('JingMA_NEW/res/Harmony/ALL/RDS/DataAva/',batch,'/features.tsv',sep=''),
              quote = F,sep = '\t',
              col.names = F,row.names = F)
  ##
  write.table(colnames(exp),file = paste('JingMA_NEW/res/Harmony/ALL/RDS/DataAva/',batch,'/barcodes.tsv',sep=''),
              quote = F,
              col.names = F,row.names = F)
  ##
  file=paste('JingMA_NEW/res/Harmony/ALL/RDS/DataAva/',batch,'/matrix.mtx',sep='')
  sink(file)
  cat("%%MatrixMarket matrix coordinate integer general\n")
  cat("%\n")
  cat(paste(nrow(exp),ncol(exp),sum(exp>0),"\n")) 
  sink()
  
  tmp=do.call(rbind,lapply(1:ncol(exp),function(i){
    return(data.frame(row=1:nrow(exp),
                      col=i,
                      exp=exp[,i]))
  }) )
  tmp=tmp[tmp$exp>0,]
  write.table(tmp,file = paste('JingMA_NEW/res/Harmony/ALL/RDS/DataAva/',batch,'/matrix.mtx',sep=''),quote = F,
              col.names = F,row.names = F,append = T )
  
}

get_output10X(pbmc,'M3')

table(pbmc$batch)
data <- Read10X(data.dir='JingMA_NEW/res/Harmony/ALL/RDS/DataAva/C3/')
dim(data)



tmp <- read.table('sample',fill = TRUE)
sample <- c()
for(i in 1:nrow(tmp)){
  sample <- c(sample,unique(as.character(tmp[i,])))
}
sample <- unique(sample)
sample <- sample[sample!='']

a <- data.frame(matrix(sample,ncol = 2,byrow = TRUE))
# a <- data.frame(sample)
write.xlsx(a,'tmp.xlsx',col.names = TRUE,row.names = FALSE)


CT.df <- pbmc@meta.data[,c(2:6,12:13)]
CT.df$subcelltype <- as.character(CT.df$subcelltype)
CT.df$CT <- CT.df$subcelltype
CT.df$subcelltype[CT.df$CT=='IC'] <- 'Immune'
CT.df$subcelltype[CT.df$CT=='EC'] <- 'Endo'
CT.df$subcelltype[CT.df$CT=='CSC'] <- 'CSPC'
CT.df$subcelltype[CT.df$CT=='C0'] <- 'EC'
CT.df$subcelltype[CT.df$CT=='C1'] <- 'IC'
CT.df$subcelltype[CT.df$CT=='C2'] <- 'LC'
CT.df$subcelltype[CT.df$CT=='SSC'] <- 'SSPC'


CT.df$celltype <- as.character(CT.df$celltype)
CT.df$CT <- CT.df$celltype
CT.df$celltype[CT.df$CT=='IC'] <- 'Immune'
CT.df$celltype[CT.df$CT=='EC'] <- 'Endo'
CT.df$celltype[CT.df$CT=='C'] <- 'Chond'
CT.df$celltype[CT.df$CT=='CSC'] <- 'CSPC'
CT.df$celltype[CT.df$CT=='SSC'] <- 'SSPC'


UMAP.df <- pbmc@reductions$umap@cell.embeddings

print(all(rownames(UMAP.df)==rownames(CT.df)))

META.df <- cbind(UMAP.df,CT.df)
META.df <- META.df[,c(1:7,9,8)]
META.df <- rownames_to_column(META.df,'Cell')

#write.xlsx(META.df,'/mdshare/node9/yzj/JingMA_ORI/Data/CellMetaInfo.xlsx',row.names = FALSE)
write.table(META.df,'/mdshare/node9/yzj/JingMA_ORI/Data/CellMetaInfo.txt',sep = '\t',
            row.names = FALSE,col.names = TRUE,quote = F)


pbmc@meta.data$celltype <-CT.df$celltype
pbmc@meta.data$subcelltype <-CT.df$subcelltype

DimPlot(pbmc,group.by = 'celltype',label = T)
DimPlot(pbmc,group.by = 'subcelltype',label = T)
saveRDS(pbmc,'JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype.Rdata')


