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
## Fig3E. 成人来说上调的点图:在aging db里/不在aging db里
################
#### 对成人来说上调的基因
dnValues_mtx <- readRDS('JingMA_NEW/res/compControl/ChildrenvsAdults/FIG/DEGHeatmap_DNmtx.RDS')
dnValues_mtx <- as.data.frame(dnValues_mtx);colnames(dnValues_mtx)[2] <- 'C0'
dnValues_mtx <- tibble::rownames_to_column(dnValues_mtx,var = 'gene')
head(dnValues_mtx)

upFreq.sum <- apply(dnValues_mtx[,2:5],1,sum)
names(upFreq.sum) <- dnValues_mtx$gene

upGenes <- dnValues_mtx$gene

### 计算cell里基因表达的percent和FC
pbmc_upOL <- subset(pbmc_C,features = upGenes)
celltypes <- unique(pbmc_upOL$celltype)
OL_MKlst <- list()
for (CT in celltypes) {
  print(paste('FindMarkers:',CT,'Start!'))
  cells <- colnames(pbmc_upOL)[pbmc_upOL$celltype == CT]
  subpbmc <- subset(pbmc_upOL,cells = cells)
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

upPct_mtx <-merge(merge(merge(pct_CSC,pct_C0,by = 'gene'),pct_C1,by = 'gene'),pct_C2,by = 'gene')
upPct.sum <- apply(upPct_mtx[,2:5],1,mean);names(upPct.sum) <- upPct_mtx$gene
upPct.max <- apply(upPct_mtx[,2:5],1,max);names(upPct.max) <- upPct_mtx$gene

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
upFC_mtx <-merge(merge(merge(fc_CSC,fc_C0,by = 'gene'),fc_C1,by = 'gene'),fc_C2,by = 'gene')
upFC.sum <- apply(upFC_mtx[,2:5],1,mean);names(upFC.sum) <- upFC_mtx$gene

### 合并 Freq, pct, FC
upSum_mtx <- merge(merge(merge(data.frame(gene=names(upFreq.sum),freq=upFreq.sum),data.frame(gene=names(upFC.sum),fc=upFC.sum),by='gene'),
                   data.frame(gene=names(upPct.sum),pct=upPct.sum),by='gene'),data.frame(gene=names(upPct.max),pct.max=upPct.max),by='gene')
upSum_mtx$agingdb <- 'no'
upSum_mtx$agingdb[upSum_mtx$gene %in% Aging.genes] <- 'yes'
upSum_mtx <-upSum_mtx[order(-upSum_mtx$freq,-upSum_mtx$pct,upSum_mtx$fc),]

## pct.max > 0.8
pct.max <- upSum_mtx[upSum_mtx$pct.max > 0.5,]

################
## 导入 GO term
GO <- read.gmt('/home/yzj/publicData/GMT/c5.all.v6.2.symbols.gmt')
#ECM_as <- GO[GO$term == 'GO_EXTRACELLULAR_MATRIX_ASSEMBLY','gene'] 
ECM_dis <- GO[GO$term=='GO_EXTRACELLULAR_MATRIX_DISASSEMBLY','gene']
Inflam <- GO[GO$term == 'GO_INFLAMMATORY_RESPONSE','gene']
Aging <- GO[GO$term == 'GO_AGING','gene']
Apop <- GO[GO$term %in% c('GO_EXTRINSIC_APOPTOTIC_SIGNALING_PATHWAY','GO_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY'),'gene']
ROS <- GO[GO$term=='GO_RESPONSE_TO_OXIDATIVE_STRESS','gene']

ECM_comp <- GO[GO$term == 'GO_EXTRACELLULAR_MATRIX_COMPONENT','gene']
chond_diff <- GO[GO$term %in% c('GO_CHONDROCYTE_DIFFERENTIATION','GO_CHONDROCYTE_DEVELOPMENT'),'gene']

## 在aging db里
ECM_dis.genelst.inAging <- intersect(intersect(pct.max$gene,ECM_dis),Aging.genes)
Inflam.genelst.inAging <- intersect(intersect(pct.max$gene,Inflam),Aging.genes)
Aging.genelst.inAging <- intersect(intersect(pct.max$gene,Aging),Aging.genes)
Apop.genelst.inAging <- intersect(intersect(pct.max$gene,Apop),Aging.genes)
ROS.genelst.inAging <-intersect(intersect(pct.max$gene,ROS),Aging.genes)
show.upGene.inAging <-unique(c(ECM_dis.genelst.inAging,Inflam.genelst.inAging,Aging.genelst.inAging,Apop.genelst.inAging,ROS.genelst.inAging))

## 不在aging db里
ECM_dis.genelst.noAging <- setdiff(intersect(pct.max$gene,ECM_dis),Aging.genes)
Inflam.genelst.noAging <- setdiff(intersect(pct.max$gene,Inflam),Aging.genes)
Aging.genelst.noAging <- setdiff(intersect(pct.max$gene,Aging),Aging.genes)
Apop.genelst.noAging <- setdiff(intersect(pct.max$gene,Apop),Aging.genes)
ROS.genelst.noAging <-setdiff(intersect(pct.max$gene,ROS),Aging.genes)
show.upGene.noAging <- unique(c(ECM_dis.genelst.noAging,Inflam.genelst.noAging,Aging.genelst.noAging,Apop.genelst.noAging,ROS.genelst.noAging))

Cond='in'
if(Cond=='in'){
  ## 在aging db里
  show.gene <- setdiff(show.upGene.inAging,c('SOD2','ATF4'))
}else if(Cond=='no'){
  ## 不在aging db里
  show.gene <- c('CEBPD',intersect(show.upGene.noAging,c('APOD','ICAM1','SOCS3','GPX3','NFKBIZ','FOSL1','CHI3L1','ID2','BAG3')))
}

OL_meltdf <- dnValues_mtx[dnValues_mtx$gene %in% show.gene,]
OL_meltdf <- reshape2::melt(OL_meltdf)
colnames(OL_meltdf) <- c('gene','CellType','DEG')

pct_meltdf <- upPct_mtx[upPct_mtx$gene %in% show.gene,]
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
UP.plot <- 
  ggplot(data = plot.df, 
         aes(x = CellType, y = gene, size = pct, color = group)) + 
  geom_point(fill = 'cornsilk') + 
  scale_color_manual(values = c('#B15E72','#EFEFEF'))+
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
UP.plot
if(Cond=='in'){
  ggsave('JingMA_NEW/res/compControl/ChildrenvsAdults/FIG/Fig3E_upGA_dotplot_inAging.pdf',UP.plot,height = 6,width = 5.9,units = 'cm')
}else if(Cond=='no'){
  ggsave('JingMA_NEW/res/compControl/ChildrenvsAdults/FIG/Fig3E_upGA_dotplot_noAging.pdf',UP.plot,height = 4,width = 5.2,units = 'cm')
}
