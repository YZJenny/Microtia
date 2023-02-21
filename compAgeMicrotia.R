#########
## 比较老人 和 疾病 (C5 vs M1/M2/M3)
#########
library(Seurat)
library(ggplot2)
library(ggrepel)

###############
#  step1. DEG in children and adults
###############
library(future)
plan("multiprocess", workers = 20)

pbmc_chond <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype_Chond.Rdata')
Idents(pbmc_chond) <- pbmc_chond$batch
pbmc_AM <- subset(pbmc_chond,idents = c('C5','M1','M2','M3') )
pbmc_AM@meta.data$Phase <- 'Aged'
pbmc_AM$Phase[pbmc_AM$batch %in% c('M1','M2','M3')] <- 'Microtia'
pbmc_AM$Phase <- factor(pbmc_AM$Phase,levels = c('Microtia','Aged'))
pbmc_AM <- saveRDS(pbmc_AM,'JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype_AgedMicrotia_Chond.Rdata')

## 1.1 filtering gene with cutoff
batch1='Aged'
batch2='Microtia'
resol='celltype'
MK.lst <- list()
celltypes <- unique(pbmc_AM$celltype)

for (CT in celltypes) {
  print(paste('FindMarkers:',CT,'Start!'))
  cells <- colnames(pbmc_AM)[pbmc_AM$celltype == CT]
  subpbmc <- subset(pbmc_AM,cells = cells)
  Idents(subpbmc) <- subpbmc$Phase
  cell_numbers <- as.numeric(table(subpbmc$Phase))
  if(length(cell_numbers) == 2 & all(cell_numbers>3)){
    sub.markers <- FindMarkers(subpbmc, ident.1 = batch1,ident.2 = batch2)
    MK.lst[[as.character(CT)]] <- sub.markers
  }else{
    print(paste(CT,':no correct cell numbers!'))
  }
}
saveRDS(MK.lst, file = paste('JingMA_NEW/res/Harmony/ALL/RDS/DEGs_inChond_in',batch1,batch2,'.RDS',sep=''))


################
## step2. pick DEG
################
library(xlsx)
library(Seurat)
library(ggplot2)

pbmc_AM <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype_AgedMicrotia_Chond.Rdata')
VG <- VariableFeatures(pbmc_AM)

MK.lst <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/DEGs_inChond_inAgedMicrotia.RDS')
wp= 'JingMA_NEW/res/compMicrotia/MicrotiavsNormal_inC5'

FC <- 1.5
adjP <- 0.05

CellTypes <- unique(pbmc_AM$celltype)
upNum <- c()
dnNum <- c()
for(i in 1:length(CellTypes)){
  CT=CellTypes[i]
  print(CT)
  
  op=paste(wp,'/DEG/FC',FC,'_adjP',adjP,'/',sep='')
  ## mkdir
  if(!file.exists(file.path(op))){
    dir.create(file.path(op),recursive = TRUE)
  }
  ##
  data<- MK.lst[[CT]]
  data <- tibble::rownames_to_column(data,'gene')

  # sigUP_gene <- data[data$gene %in% VG & data$p_val_adj  < adjP & data$avg_logFC > log(FC),]
  # sigDN_gene <- data[data$gene %in% VG & data$p_val_adj  < adjP & data$avg_logFC < -log(FC),]
  
  sigUP_gene <- data[data$p_val_adj  < adjP & data$avg_logFC > log(FC),]
  sigDN_gene <- data[data$p_val_adj  < adjP & data$avg_logFC < -log(FC),]
  
  sigUP_gene <- sigUP_gene[order(sigUP_gene$avg_logFC,decreasing = T),]
  sigDN_gene <- sigDN_gene[order(sigDN_gene$avg_logFC,decreasing = F),]
  
  upNum <- c(upNum,nrow(sigUP_gene))
  dnNum <- c(dnNum,nrow(sigDN_gene))
  ## 相对于疾病来说,上调的基因
  write.xlsx(sigUP_gene, file = paste(op,CT,'.xlsx',sep=''), row.names = FALSE, sheetName = "UP")
  ## 相对于疾病来说,下调的基因
  write.xlsx(sigDN_gene, file = paste(op,CT,'.xlsx',sep=''), row.names = FALSE, sheetName = "DN",append=TRUE)
}

Num <- data.frame(CellType=CellTypes,upGene=upNum,dnGene=dnNum)
Num$CellType <- factor(Num$CellType,levels = c("CSC","C0","C1","C2"))
library(reshape2)
NUM_meltdf <- melt(Num)


#####
# step3. Clusterprofiler
#####
library(xlsx)
data(c2BroadSets)
library(Biobase)
library(genefilter)
library(limma)
library(RColorBrewer)
library(AnnotationHub)
library(org.Hs.eg.db)   #人类注释数据库
library(clusterProfiler)
library(DOSE)
library(dplyr)
library(reshape2)

MK.lst <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/DEGs_inChond_inAgedMicrotia.RDS')
wp= 'JingMA_NEW/res/compMicrotia/MicrotiavsNormal_inC5'

FC <- 1.5
adjP <- 0.05

CellTypes=c('CSC','C0','C1','C2')
for(i in 1:length(CellTypes)){
  CT <- CellTypes[i]
  print(CT)
  
  op=paste(wp,'/ClusterPro/FC',FC,'_adjP',adjP,'/',sep='')
  ## mkdir
  if(!file.exists(file.path(op))){
    dir.create(file.path(op),recursive = TRUE)
  }
  
  TYPE <- c('UP','DN')
  for(j in 1:length(TYPE)){
    type=TYPE[j]
    sigG_lst <- read.xlsx(paste(wp,'/DEG/FC',FC,'_adjP',adjP,'/',CT,'.xlsx',sep=''),sheetName = type)[,1]
    testgenes<- sigG_lst
    symbol2id=mapIds(org.Hs.eg.db,testgenes,"ENTREZID",'SYMBOL')
    id=symbol2id[which(symbol2id!='')] #提取出非NA的ENTREZID
    
    if(length(id) > 0){
      #GO BP 富集分析#
      ego <- enrichGO(OrgDb="org.Hs.eg.db", gene = id, ont = "BP", pvalueCutoff = 0.05, readable= TRUE) #GO富集分析
      ego <- clusterProfiler::simplify(ego,cutoff=0.8,by="p.adjust",select_fun=min,measure="Wang")
      ego_res <- as.data.frame(ego)
      if(nrow(ego_res)>0){
        write.xlsx(ego_res,paste(op,CT,"_BP.xlsx",sep=''),row.names = FALSE,sheetName = type,append = TRUE)
      } 
    }
  }
}



##################
### step4：Geneset score的boxplot
##################
library(xlsx)
library(ggsignif)
library(ggpubr)
pbmc_AM <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype_AgedMicrotia_Chond.Rdata')
geneSet <- read.xlsx('JingMA_NEW/res/compControl/ChildrenvsAdults/geneset_skinAging.xlsx',sheetIndex = 4)
geneSet <- geneSet[,c('DNA.repair.genes','DNA.damage.genes','Inflammatory.response.genes','ROS.genes','NF.κB.pathway.genes.')]
colnames(geneSet)

for(i in 1:length(colnames(geneSet))){
  set.name=colnames(geneSet)[i]
  print(set.name)
  prename <- paste(unlist(strsplit(set.name,'\\.')),collapse ='_')
  pbmc.score <- AddModuleScore(pbmc_AM,features = list(na.omit(geneSet[,set.name])),name = 'score')
  dt <- data.frame(cells = colnames(pbmc.score),score=pbmc.score$score1,celltype=pbmc.score$celltype,phase=pbmc.score$Phase)
  dt$celltype <- as.character(dt$celltype)
  dt$CT=dt$celltype
  dt$celltype[dt$CT=='CSC']='CSPC'
  dt$celltype[dt$CT=='C0']='EC'
  dt$celltype[dt$CT=='C1']='IC'
  dt$celltype[dt$CT=='C2']='LC'
  
  p <- ggplot(dt, aes(x=celltype, y=score,fill=phase)) + 
    geom_boxplot(width=0.5,position=position_dodge(0.9),outlier.size = 0.05)+
    scale_fill_manual(values = c( "#6C6C6C","#DC143C"))+ 
    theme_bw()+
    theme(axis.text.x=element_blank(),axis.ticks.x =element_blank(),
          axis.text.y=element_text(size=6,colour="black"),axis.title.y=element_text(size = 6,colour="black"), 
          axis.ticks.y =element_line(colour="black",size = 0.01),
          panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5,size = 6),
          plot.margin = unit(c(0.1,0,-0.3,-0.2), "cm"))+ 
    ylab("")+xlab("")+ labs(title = set.name)+
    facet_wrap(~celltype,ncol = 4,scales= "free_x")+
    theme(strip.background = element_rect(color = "black", fill = "white",size = 0.5),
          strip.text.x = element_text(size = 6, color = "black",face = 'bold'),
          panel.grid = element_blank(),panel.border = element_rect(color = 'black',size = 0.3))+
    stat_compare_means(method = 'wilcox.test',label = "p.signif",label.x=0.1)
  ggsave(paste('JingMA_NEW/res/compMicrotia/MicrotiavsNormal_inC5/FIG/sFig_',prename,'.pdf',sep=''),p,width = 4,height = 2)
}


SASP <- read.gmt('publicData/GMT/SASP_reactome.gmt')[,'gene']
pbmc.score <- AddModuleScore(pbmc_AM,features = list(SASP),name = 'score')
dt <- data.frame(cells = colnames(pbmc.score),score=pbmc.score$score1,celltype=pbmc.score$celltype,phase=pbmc.score$Phase)
dt$celltype <- as.character(dt$celltype)
dt$CT=dt$celltype
dt$celltype[dt$CT=='CSC']='CSPC'
dt$celltype[dt$CT=='C0']='EC'
dt$celltype[dt$CT=='C1']='IC'
dt$celltype[dt$CT=='C2']='LC'


p <- ggplot(dt, aes(x=celltype, y=score,fill=phase)) + 
  geom_boxplot(width=0.5,position=position_dodge(0.9),outlier.size = 0.05)+
  scale_fill_manual(values = c( "#6C6C6C","#DC143C"))+ 
  theme_bw()+
  theme(axis.text.x=element_blank(),axis.ticks.x =element_blank(),
        axis.text.y=element_text(size=6,colour="black"),axis.title.y=element_text(size = 6,colour="black"), 
        axis.ticks.y =element_line(colour="black",size = 0.01),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 6),
        plot.margin = unit(c(0.1,0,-0.3,-0.2), "cm"))+ 
  ylab("")+xlab("")+ labs(title = 'SASP')+
  facet_wrap(~celltype,ncol = 4,scales= "free_x")+
  theme(strip.background = element_rect(color = "black", fill = "white",size = 0.5),
        strip.text.x = element_text(size = 6, color = "black",face = 'bold'),
        panel.grid = element_blank(),panel.border = element_rect(color = 'black',size = 0.3))+
  stat_compare_means(method = 'wilcox.test',label = "p.signif",label.x=0.1)
ggsave('JingMA_NEW/res/compMicrotia/MicrotiavsNormal_inC5/FIG/sFig_SASP.pdf',p,width = 4,height = 2)


c5_GO <- read.gmt('publicData/GMT/c5.all.v6.2.symbols.gmt')
Aging_GO <- c5_GO$gene[c5_GO$term=='GO_AGING']
pbmc.score <- AddModuleScore(pbmc_AM,features = list(Aging_GO),name = 'score')
dt <- data.frame(cells = colnames(pbmc.score),score=pbmc.score$score1,celltype=pbmc.score$celltype,phase=pbmc.score$Phase)
dt$celltype <- as.character(dt$celltype)
dt$CT=dt$celltype
dt$celltype[dt$CT=='CSC']='CSPC'
dt$celltype[dt$CT=='C0']='EC'
dt$celltype[dt$CT=='C1']='IC'
dt$celltype[dt$CT=='C2']='LC'


p <- ggplot(dt, aes(x=celltype, y=score,fill=phase)) + 
  geom_boxplot(width=0.5,position=position_dodge(0.9),outlier.size = 0.05)+
  scale_fill_manual(values = c( "#6C6C6C","#DC143C"))+ 
  theme_bw()+
  theme(axis.text.x=element_blank(),axis.ticks.x =element_blank(),
        axis.text.y=element_text(size=6,colour="black"),axis.title.y=element_text(size = 6,colour="black"), 
        axis.ticks.y =element_line(colour="black",size = 0.01),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 6),
        plot.margin = unit(c(0.1,0,-0.3,-0.2), "cm"))+ 
  ylab("")+xlab("")+ labs(title = 'Aging')+
  facet_wrap(~celltype,ncol = 4,scales= "free_x")+
  theme(strip.background = element_rect(color = "black", fill = "white",size = 0.5),
        strip.text.x = element_text(size = 6, color = "black",face = 'bold'),
        panel.grid = element_blank(),panel.border = element_rect(color = 'black',size = 0.3))+
  stat_compare_means(method = 'wilcox.test',label = "p.signif",label.x=0.1)
ggsave('JingMA_NEW/res/compMicrotia/MicrotiavsNormal_inC5/FIG/sFig_Aging.pdf',p,width = 4,height = 2)



# #####
# # step5. 热图
# #####
# library(tibble)
# library(dplyr)
# pbmc <- readRDS('/home/yzj/JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype.Rdata')
# Idents(pbmc) <- pbmc$type
# pbmc_AM <- subset(pbmc,idents = 'Control')
# pbmc_AM@meta.data$Phase <- 'Adult'
# pbmc_AM$Phase[pbmc_AM$batch %in% c('C4','C6')] <- 'Children'
# pbmc_AM$Phase <- factor(pbmc_AM$Phase,levels = c('Children','Adult'))
# 
# batch1='Children'
# batch2='Adult'
# resol='celltype'
# MK.lst <- readRDS(paste('JingMA_NEW/res/Harmony/ALL/RDS/DEGs_inChond_in',batch1,batch2,'.RDS',sep=''))
# 
# wp=paste('/home/yzj/JingMA_NEW/res/compControl/',batch1,'vs',batch2,sep='')
# op=paste('/home/yzj/JingMA_NEW/res/compControl/',batch1,'vs',batch2,'/Heatmap',sep='')
# ## mkdir
# if(!file.exists(file.path(op))){
#   dir.create(file.path(op),recursive = TRUE)
# }
# 
# CellTypes=c('CSC','TC',"Chondrocyte",'C1','C2')
# 
# for(i in 1:length(CellTypes)){
#   CT=CellTypes[i]
#   print(CT)
#   
#   if(CT=='Chondrocyte'){
#     cells <- colnames(pbmc_AM)[pbmc_AM$celltype %in% c('C1','C2')]
#   }else{
#     cells <- colnames(pbmc_AM)[pbmc_AM$celltype==CT]
#   }
#   subpbmc <- subset(pbmc_AM,cells = cells)
#   Idents(subpbmc) <- subpbmc$Phase
#   
#   data <- MK.lst[[CT]]
#   data$gene <- rownames(data)
#   up <- data %>% as_tibble() %>% arrange(desc(avg_logFC)) %>% filter(p_val_adj < 0.01) %>%
#     head(n=50)
#   dn <- data %>% as_tibble() %>% arrange(desc(avg_logFC)) %>% filter(p_val_adj < 0.01) %>%
#     tail(n= 50)
#   genes <- c(up$gene,dn$gene)
#   pdf(paste(op,'/',CT,'.pdf',sep=''),height = 13,width = 5)
#   p <- DoHeatmap(subpbmc,features = genes,draw.lines = F,angle = 0)
#   print(p)
#   dev.off()
# }

# 
# ######
# ## step8. 火山图,并加上感兴趣基因的名字
# ######
# batch1='Children'
# batch2='Adult'
# resol='celltype'
# MK.lst <- readRDS(paste('JingMA_NEW/res/Harmony/ALL/RDS/DEGs_inChond_in',batch1,batch2,'.RDS',sep=''))
# 
# wp=paste('/home/yzj/JingMA_NEW/res/',batch1,'vs',batch2,sep='')
# op=paste('/home/yzj/JingMA_NEW/res/',batch1,'vs',batch2,'/Volcano',sep='')
# 
# FC <- 2
# adjP <- 0.01
# 
# CellTypes=c('CSPC')
# for(i in 1:length(CellTypes)){
#   CT=CellTypes[i]
#   print(CT)
#   CUT=paste('FC',FC,'_adjP',adjP,sep='')
#   TAG=paste('FC',FC,'_adjP',adjP,'/',CT,sep='')
#   ## mkdir
#   if(!file.exists(file.path(op, CUT))){
#     dir.create(file.path(op, CUT),recursive = TRUE)
#   }
#   data <- MK.lst[[CT]]
#   data$gene <- rownames(data)
#   data$significant <- as.factor(ifelse(data$p_val_adj < adjP & abs(data$avg_logFC) >= log(FC),
#                                        ifelse(data$avg_logFC >= log(FC),'UP','DOWN'),'NOT'))
#   print(table(data$significant))
# 
#   gene<-data$gene
#   gene_logFC<-data$avg_logFC
#   gene_logAdjP<--log10(data$p_val_adj)
#   boolean_isInf<-is.infinite(gene_logFC)
#   filter_gene<-as.vector(gene[!boolean_isInf])
#   filter_gene_logFC<-gene_logFC[!boolean_isInf]
#   filter_gene_logAdjP<-gene_logAdjP[!boolean_isInf]
#   significant<-data$significant[!boolean_isInf]
#   filter_data<-data.frame(filter_gene,filter_gene_logFC,filter_gene_logAdjP,significant)
#   boolean_tag=(filter_gene_logFC> log(FC) | filter_gene_logFC< -log(FC)) & filter_gene_logAdjP > -log10(adjP)
#   filter_gene[which(boolean_tag==F)]=""
# 
#   ## CSPC
#   tag_lst=c("SPARC","ELN","TIMP4","MMP3","FST","IL6","IL8","CDKN1A","MMP1")
#   index_tag=which(!filter_gene %in% tag_lst)
#   filter_gene[index_tag]=""
# 
# 
#   pdf(paste(op,'/',CUT,'/',CT,'.pdf',sep=""))
#   figure<-ggplot(filter_data,aes(x=filter_gene_logFC,y=filter_gene_logAdjP))
#   xstt <- floor(min(filter_data$filter_gene_logFC))
#   xend <- ceiling(max(filter_data$filter_gene_logFC))
#   x <- max(abs(xstt),abs(xend))
#   xstt <- -x
#   xend <- x
# 
#   ystt <- floor(min(filter_data$filter_gene_logAdjP))
#   yend <- ceiling(max(filter_data$filter_gene_logAdjP))
# 
#   figure1 <- figure+geom_point(aes(color=significant))+xlim(xstt,xend) + ylim(ystt,yend)+
#     scale_color_manual(values = c("#377EB8","#999999","#E41A1C"))+
#     labs(title=paste(batch1,'vs',batch2,'in',CT),x="logFC",y="-log10(adjP)")+
#     theme(plot.title = element_text(hjust = 0.5))+
#     geom_hline(yintercept = -log10(adjP),linetype=3)+
#     geom_vline(xintercept=c(-log(FC),log(FC)),linetype=3)
#   print(figure1+geom_text_repel(label=filter_gene))
#   dev.off()
# }

