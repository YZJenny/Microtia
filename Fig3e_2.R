library(Seurat)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)

pbmc_C <- readRDS('/home/yzj/JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype_Control_Chond.Rdata')
pbmc_C$Phase <- factor(pbmc_C$Phase,levels = c('Children','Adults'))

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

################
## Fig3E 成人来说上调的点图:在aging db里/不在aging db里
################
#### 对成人来说下调的基因
upValues_mtx <- readRDS('JingMA_NEW/res/compControl/ChildrenvsAdults/FIG/DEGHeatmap_UPmtx.RDS')
upValues_mtx <- as.data.frame(upValues_mtx);colnames(upValues_mtx)[2] <- 'C0'
upValues_mtx <- tibble::rownames_to_column(upValues_mtx,var = 'gene')
head(upValues_mtx)

upFreq.sum <- apply(upValues_mtx[,2:5],1,sum)
names(upFreq.sum) <- upValues_mtx$gene

dnGenes <- upValues_mtx$gene

### 计算cell里基因表达的percent和FC
pbmc_dnOL <- subset(pbmc_C,features = dnGenes)
celltypes <- unique(pbmc_dnOL$celltype)
OL_MKlst <- list()
for (CT in celltypes) {
  print(paste('FindMarkers:',CT,'Start!'))
  cells <- colnames(pbmc_dnOL)[pbmc_dnOL$celltype == CT]
  subpbmc <- subset(pbmc_dnOL,cells = cells)
  Idents(subpbmc) <- subpbmc$Phase
  cell_numbers <- as.numeric(table(subpbmc$Phase))
  if(length(cell_numbers) == 2 & all(cell_numbers>3)){
    sub.markers <- FindMarkers(subpbmc, ident.1 = 'Children',ident.2 = 'Adults',min.pct = 0,logfc.threshold = 0)
    OL_MKlst[[as.character(CT)]] <- sub.markers
  }else{
    print(paste(CT,':no correct cell numbers!'))
  }
}

### pct
get_pct <- function(CT,method='diff'){
  df <- OL_MKlst[[CT]]
  if(method=='mean'){
    pct=as.numeric(apply(df[,c('pct.1','pct.2')],1,mean))
  }else if(method=='diff'){
    pct=as.numeric(df[,'pct.2']-df[,'pct.1'])
  }
  names(pct) <- rownames(df)
  pct_df <- data.frame(gene=names(pct),pct=pct)
  colnames(pct_df)[2] <- CT
  return(pct_df)
}
pct_CSC <- get_pct(CT='CSC',method = 'mean');
pct_C0 <- get_pct(CT='C0',method = 'mean')
pct_C1 <- get_pct(CT='C1',method = 'mean')
pct_C2 <- get_pct(CT='C2',method = 'mean')

dnPct_mtx <-merge(merge(merge(pct_CSC,pct_C0,by = 'gene'),pct_C1,by = 'gene'),pct_C2,by = 'gene')
dnPct.sum <- apply(dnPct_mtx[,2:5],1,mean);names(dnPct.sum) <- dnPct_mtx$gene
dnPct.max <- apply(dnPct_mtx[,2:5],1,max);names(dnPct.max) <- dnPct_mtx$gene

### FC
get_fc <- function(CT){
  df <- OL_MKlst[[CT]]
  fc <- df[,'avg_logFC']
  names(fc) <- rownames(df)
  fc_df <- data.frame(gene=names(fc),fc=fc)
  colnames(fc_df)[2] <- CT
  return(fc_df)
}

fc_CSC=get_fc('CSC');fc_C0=get_fc('C0');fc_C1=get_fc('C1');fc_C2=get_fc('C2')
dnFC_mtx <-merge(merge(merge(fc_CSC,fc_C0,by = 'gene'),fc_C1,by = 'gene'),fc_C2,by = 'gene')
dnFC.sum <- apply(dnFC_mtx[,2:5],1,mean);names(dnFC.sum) <- dnFC_mtx$gene

### 合并 Freq, pct, FC
dnSum_mtx <- merge(merge(merge(data.frame(gene=names(upFreq.sum),freq=upFreq.sum),data.frame(gene=names(dnFC.sum),fc=dnFC.sum),by='gene'),
                         data.frame(gene=names(dnPct.sum),pct=dnPct.sum),by='gene'),data.frame(gene=names(dnPct.max),pct.max=dnPct.max),by='gene')
dnSum_mtx$agingdb <- 'no'
dnSum_mtx$agingdb[dnSum_mtx$gene %in% Aging.genes] <- 'yes'
dnSum_mtx <-dnSum_mtx[order(-dnSum_mtx$freq,-dnSum_mtx$pct,dnSum_mtx$fc),]

## pct.max > 0.5
pct.max <- dnSum_mtx[dnSum_mtx$pct.max > 0.5,]

################
## 导入 GO term
GO <- read.gmt('/home/yzj/publicData/GMT/c5.all.v6.2.symbols.gmt')
#ECM_as <- GO[GO$term == 'GO_EXTRACELLULAR_MATRIX_ASSEMBLY','gene'] 
ECM_comp <- GO[GO$term == 'GO_EXTRACELLULAR_MATRIX_COMPONENT','gene']
chond_diff <- GO[GO$term %in% c('GO_CHONDROCYTE_DIFFERENTIATION','GO_CHONDROCYTE_DEVELOPMENT',
                                'GO_REGULATION_OF_CHONDROCYTE_DIFFERENTIATION','GO_POSITIVE_REGULATION_OF_CHONDROCYTE_DIFFERENTIATION'),'gene']
Conn <- GO[GO$term=='GO_CONNECTIVE_TISSUE_DEVELOPMENT','gene']
Ske <- GO[GO$term=='GO_SKELETAL_SYSTEM_DEVELOPMENT','gene']
## 在aging db里
# ECM_comp.genelst.inAging <- intersect(intersect(pct.max$gene,ECM_comp),Aging.genes)
# chond_diff.genelst.inAging <- intersect(intersect(pct.max$gene,chond_diff),Aging.genes)
# show.dnGene.inAging <-unique(c(ECM_comp.genelst.inAging,chond_diff.genelst.inAging))
show.dnGene.inAging <- intersect(pct.max$gene,Aging.genes)
show.dnGene.inAging <- c('CTGF','ELN','PRDX1','TXN','LMNA','JUN','HSPA1A','GAPDH','MYC')


## 不在aging db里
ECM_comp.genelst.noAging <- setdiff(intersect(pct.max$gene,ECM_comp),Aging.genes)
chond_diff.genelst.noAging <- intersect(intersect(pct.max$gene,chond_diff),Aging.genes)
conn.genelst.noAging <- intersect(intersect(pct.max$gene,Conn),Aging.genes)
ske.genelst.noAging <- intersect(intersect(pct.max$gene,Ske),Aging.genes)
show.dnGene.noAging <- unique(c(ECM_comp.genelst.noAging,chond_diff.genelst.noAging))

Cond <- 'no'
if(Cond=='in'){
  ## 在aging db里
  show.gene <- show.dnGene.inAging
}else if(Cond=='no'){
  ## 不在aging db里
  show.gene <- show.dnGene.noAging
}

OL_meltdf <- upValues_mtx[upValues_mtx$gene %in% show.gene,]
OL_meltdf <- reshape2::melt(OL_meltdf)
colnames(OL_meltdf) <- c('gene','CellType','DEG')

pct_meltdf <- dnPct_mtx[dnPct_mtx$gene %in% show.gene,]
pct_meltdf <- reshape2::melt(pct_meltdf)
colnames(pct_meltdf) <- c('gene','CellType','pct')

print(all(OL_meltdf$variable==pct_meltdf$variable))

plot.df <- merge(OL_meltdf,pct_meltdf,by = c('gene','CellType'))
plot.df$pct=plot.df$pct*100
plot.df$gene <- factor(plot.df$gene,levels = rev(intersect(pct.max$gene,show.gene)))
plot.df$group <- 'Unchanged'
plot.df$group[plot.df$DEG==1] <- 'Upregulated'
plot.df$group <- factor(plot.df$group,levels =c('Upregulated','Unchanged'))
library(ggplot2)
DN.plot <- 
  ggplot(data = plot.df, 
         aes(x = CellType, y = gene, size = pct, color = group)) + 
  geom_point(fill = 'cornsilk') + 
  scale_color_manual(values = c('#7F99CE','#EFEFEF'))+
  scale_size_continuous(range = c(1,3),breaks = c(30,60,90)) +
  theme_classic()+
  labs(x='',y='')+
  theme(axis.text.x = element_text(size = 8,colour = "black"),
        axis.text.y = element_text(size = 8,colour = "black",face = 'italic'),
        text = element_text(size = 6,colour = "black"),legend.text = element_text(size = 6, color = 'black'),
        legend.key.size = unit(0.05,'cm'),legend.key.height = unit(0.05,'cm'),
        axis.line = element_line(size=0.2, colour = "black"),
        plot.margin = unit(c(0.1,-0.1,-0.1,-0.1),units = 'cm'))+
  guides(colour = guide_legend(override.aes = list(size=1)))
DN.plot
if(Cond=='in'){
  ggsave('JingMA_NEW/res/compControl/ChildrenvsAdults/FIG/Fig3E_dnGA_dotplot_inAging.pdf',DN.plot,height = 3.5,width = 5.3,units = 'cm')
}else if(Cond=='no'){
  ggsave('JingMA_NEW/res/compControl/ChildrenvsAdults/FIG/Fig3E_dnGA_dotplot_noAging.pdf',DN.plot,height = 3.5,width = 5.5,units = 'cm')
}
