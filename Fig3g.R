######
### Fig3E CEBPB and CEBPD 的网络分析
######
library(foreach)
library(reshape2)
library(tibble)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(igraph)
library(ggraph)
library(tidygraph)
data(c2BroadSets)
library(Biobase)
library(genefilter)
library(limma)
library(AnnotationHub)
library(org.Hs.eg.db)   #人类注释数据库
library(clusterProfiler)
library(DOSE)
source('/home/yzj/JingMA_NEW/scr/func_CP.R')

upValues_mtx <- readRDS('JingMA_NEW/res/compControl/ChildrenvsAdults/FIG/DEGHeatmap_UPmtx.RDS')
dnValues_mtx <- readRDS('JingMA_NEW/res/compControl/ChildrenvsAdults/FIG/DEGHeatmap_DNmtx.RDS')

upGene <- rownames(upValues_mtx)
dnGene <- rownames(dnValues_mtx)

up_sum <- apply(upValues_mtx,1,sum)
dn_sum <- apply(dnValues_mtx,1,sum)

commonUP <- names(up_sum)[which(up_sum == 4)]
commonDN <- names(dn_sum)[which(dn_sum == 4)]

regulon2Targets <- readRDS('JingMA_NEW/res/SCENIC_main/int/2.5_regulonTargetsInfo.Rds')

TF <- unique(regulon2Targets$TF)
TF_upGene <- intersect(TF,commonUP)
TF_dnGene <- intersect(TF,commonDN) #"CEBPB","CEBPD"
print(TF_dnGene)


################
## network图
################
sub.regulon2Targets <- regulon2Targets[regulon2Targets$TF %in% TF_dnGene & regulon2Targets$gene %in% dnGene & regulon2Targets$highConfAnnot =='TRUE',]

mat_weight <- reshape2::dcast(sub.regulon2Targets, 
                              TF ~ gene, value.var = 'Genie3Weight')
rownames(mat_weight) <- mat_weight$TF
mat_weight$TF <- NULL
row_add <- setdiff(colnames(mat_weight), rownames(mat_weight))
df.row_add <- data.frame(matrix(rep(NA, length(row_add)*ncol(mat_weight)), 
                                nrow = length(row_add), ncol = ncol(mat_weight)),
                         row.names = row_add)
colnames(df.row_add) <- colnames(mat_weight)
mat_weight <- rbind(mat_weight, df.row_add)
col_add <- setdiff(rownames(mat_weight), colnames(mat_weight))
df.col_add <- data.frame(matrix(rep(NA, length(col_add)*nrow(mat_weight)), 
                                nrow = nrow(mat_weight), ncol = length(col_add)),
                         row.names = rownames(mat_weight))
colnames(df.col_add) <- col_add
mat_weight <- cbind(mat_weight, df.col_add)
mat_weight[is.na(mat_weight)] <- 0
TFgene <- sort(colnames(mat_weight))
mat_weight <- mat_weight[TFgene, TFgene]

igraph_mtx <- 
    graph_from_adjacency_matrix(t(as.matrix(mat_weight)), 
                                mode = 'directed', 
                                weighted = T, diag = T)
V(igraph_mtx)$degree <- degree(igraph_mtx, normalized = T)
V(igraph_mtx)$weight_degree <- strength(igraph_mtx)
V(igraph_mtx)$closeness_centrality <- closeness(igraph_mtx, normalized = T)
V(igraph_mtx)$betweenness_centrality <- betweenness(igraph_mtx, normalized = T)
V(igraph_mtx)$eigenvector_centrality <- evcent(igraph_mtx)$vector
V(igraph_mtx)$page_rank <- page_rank(igraph_mtx)$vector
V(igraph_mtx)$betweenness_centrality_directed <- betweenness(igraph_mtx, directed = T, normalized = T)
V(igraph_mtx)$eigenvector_centrality_directed <- evcent(igraph_mtx, directed = T)$vector

node_list <- data.frame(
    node_id = V(igraph_mtx)$name,
    degree = V(igraph_mtx)$degree,
    weight_degree = V(igraph_mtx)$weight_degree,
    closeness_centrality = V(igraph_mtx)$closeness_centrality,
    betweenness_centrality = V(igraph_mtx)$betweenness_centrality,
    eigenvector_centrality = V(igraph_mtx)$eigenvector_centrality,
    page_rank = V(igraph_mtx)$page_rank,
    betweenness_centrality_directed = V(igraph_mtx)$betweenness_centrality_directed,
    eigenvector_centrality_directed = V(igraph_mtx)$eigenvector_centrality_directed)

igraph_up <- 
    graph_from_adjacency_matrix((as.matrix(mat_weight)), 
                                mode = 'directed', 
                                weighted = T, diag = T)
V(igraph_up)$page_rank <- page_rank(igraph_mtx)$vector

## node group
# library(xlsx)
# GO <- read.xlsx('JingMA_NEW/res/compControl/ChildrenvsAdults/ClusterPro/FC1.5_adjP0.05/C2_BP_all.xlsx',sheetName = 'UP')
# ILgene <- unique(c(unlist(strsplit(GO[GO$Description=='response to interleukin-6','geneID'],'\\/')),
#                    unlist(strsplit(GO[GO$Description=='regulation of inflammatory response','geneID'],'\\/'))))
# 
# Apopgene <- unlist(strsplit(GO[GO$Description=='intrinsic apoptotic signaling pathway','geneID'],'\\/'))
# 
# ROSgene <- unique(c(unlist(strsplit(GO[GO$Description=='response to oxidative stress','geneID'],'\\/')),
#                     unlist(strsplit(GO[GO$Description=='reactive oxygen species biosynthetic process','geneID'],'\\/')),
#                     unlist(strsplit(GO[GO$Description=='reactive oxygen species metabolic process','geneID'],'\\/')),
#                     unlist(strsplit(GO[GO$Description=='cell death in response to oxidative stress','geneID'],'\\/'))
# ))
# ECM <- unlist(strsplit('MMP3/TNFRSF11B/DCN/PTX3/MMP1/CYP1B1/LAMB3/RGCC/FN1/PLA2G2A/ICAM1/SOX9/ELF3/CD44/SDCBP/IL6','\\/'))
# 
# Inflam <- unique(c(unlist(strsplit(GO[GO$Description=='ERK1 and ERK2 cascade','geneID'],'\\/')),
#                  unlist(strsplit(GO[GO$Description=='p38Inflam cascade','geneID'],'\\/'))))
# 
# Senes <- unique(c(unlist(strsplit(GO[GO$Description=='aging','geneID'],'\\/')),
#                         unlist(strsplit(GO[GO$Description=='cellular senescence','geneID'],'\\/')),
#                         c('MMP3','MMP1','IL6','CXCL1','CXCL2','BIRC3','CEBPB','CDKN1A','CCL20','NFKB1','CXCL6','NFKBIA',
# 'HSPA9','HMGB2','C2orf40','MT-CO1','JUND','MMP10')))



## 导入 GO term
GO <- read.gmt('/home/yzj/publicData/GMT/c5.all.v6.2.symbols.gmt')
Apop <- GO[GO$term %in% c('GO_EXTRINSIC_APOPTOTIC_SIGNALING_PATHWAY','GO_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY'),'gene']
ECM <- GO[GO$term=='GO_EXTRACELLULAR_MATRIX_DISASSEMBLY','gene']
Inflam <- GO[GO$term == 'GO_INFLAMMATORY_RESPONSE','gene']
ROS <- GO[GO$term=='GO_RESPONSE_TO_OXIDATIVE_STRESS','gene']
Senes <- GO[GO$term %in% c('GO_AGING','GO_CELLULAR_SENESCENCE'),'gene']

geneset<- unique(sub.regulon2Targets$gene)
Apop.gene <- intersect(Apop,geneset)
ECM.gene <- intersect(ECM,geneset)
Inflam.gene <- intersect(Inflam,geneset)
ROS.gene <- intersect(ROS,geneset)
Senes.gene <- intersect(Senes,geneset)

vec_desc <- c('Apoptotic signaling pathway','Cellular senescence','Extracellular matrix degradation',
              'Inflammatory response','Response to oxidative stress')
group <- c('Apop','Senes','ECM','Inflam','ROS')

library(RColorBrewer)
mycolor=brewer.pal(10,"Set3")
col=mycolor[c(1,3:6)]
names(col) <- group

node_type <- V(igraph_mtx)$name
node_group <- rep('7', length(node_type))
node_group[node_type %in% c(Apopgene)] <- '1'
node_group[node_type %in% c(Senes)] <- '2'
node_group[node_type %in% c(ECM)] <- '3'
node_group[node_type %in% c(Inflam)] <- '4'
node_group[node_type %in% c(ROSgene)] <- '5'
node_group[node_type %in% c('CEBPB', 'CEBPD')] <- '6'
V(igraph_up)$group <- node_group


ggraph <- as_tbl_graph(igraph_up)
s=5
p <- 
    ggraph(ggraph, layout = 'stress') + 
    geom_edge_link(aes(edge_width=weight, alpha = weight),color="gray",
                   arrow = arrow(length = unit(2, 'mm')), 
                   end_cap = circle(1, 'mm'), linejoin = 'bevel', 
                   start_cap = circle(0.3, 'mm')) +
    scale_edge_width(range=c(0.6,1)) + 
    scale_size_continuous(range = c(3,10)) + 
    geom_node_point(aes(size = page_rank, fill = group, color = group, alpha = group),
                    shape=21) + 
    scale_color_manual(values = c(col,"firebrick3", 'gray')) + ## 没有Apop独特基因
    scale_fill_manual(values = c(col,"firebrick3", 'gray')) + ## 没有Apop独特基因
    scale_alpha_manual(values = c( rep(1,length(vec_desc)+2), 0.1)) + 
    geom_node_text(aes(filter = group == 1,label=name),size=s, repel = T) + 
    geom_node_text(aes(filter = group == 2,label=name),size=s, repel = T) + 
    geom_node_text(aes(filter = group == 3,label=name),size=s, repel = T) + 
    geom_node_text(aes(filter = group == 4,label=name),size=s, repel = T) + 
    geom_node_text(aes(filter = group == 5,label=name),size=s, repel = T) + 
    geom_node_text(aes(filter = group == 6,label=name),size=s+1, repel = T) + 
    geom_node_text(aes(filter = group == 7,label=name),size=s-1) +
    theme_void() + theme(legend.position = 'none')
p
ggsave(plot = p, path = '/home/yzj/JingMA_NEW/res/compControl/ChildrenvsAdults/FIG',
       filename = 'Fig3E_TFnet.pdf',height = 8, width = 12, units = 'cm')


show.genes <- union(union(union(union(union(Senes,Apopgene),ECM),ILgene),Inflam),ROSgene)
for(i in 1:length(show.genes)){
    g=show.genes[i]
    print(paste('Gene: ',g))
    if(g %in%  Senes){print('cellular senescence')}
    if(g %in%  Apopgene){print('Apop')}
    if(g %in%  ECM){print('ECM')}
    if(g %in% ILgene){print('IL')}
    if(g %in%  Inflam){print('Inflam')}
    if(g %in%  ROSgene){print('ROS')}
}

### 画饼图
set <- c('Inflam','ROS')
pct <- data.frame(group=set,pct=rep(100/length(set),length(set)),ncol = 1)
p<- ggplot(pct,aes(x="",y=pct,fill=group)) +
    geom_bar(stat = "identity",color="white",size =0.1) + 
    scale_fill_manual(values = col[set]) +
    coord_polar(theta = "y") +
    theme(axis.text.x = element_blank(),axis.title = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),panel.background =  element_blank())+guides(fill=FALSE)
p
ggsave(paste(paste(set,collapse = '_'),'.pdf',sep=''),p,width = 3,height = 3,units = 'cm')


#####################################################################################
#####################################################################################
##################
## SFig5h  fisher for CEBPB/CEBPD with database
##################
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
### gene associated aging
GA_mtx <- read.csv('publicData/GenAge/human_genes/genage_human.csv')
Aging.db <- GA_mtx[,2]
length(Aging.db)

library(xlsx)
Aging.BIG_df <- read.xlsx('publicData/GeneAtlas/aging_list_V2.0.xlsx',sheetIndex = 1)
Aging.BIG_df <- as.data.frame(Aging.BIG_df)
Aging.BIG <- Aging.BIG_df$Symbol
length(Aging.BIG)

#Aging.genes <- union(union(Aging.db,Aging.GO),Aging.BIG)
Aging.genes <- union(Aging.db,Aging.BIG)
length(Aging.genes)
Aging.lst <- list(GenAge=Aging.db,AgingAtlas=Aging.BIG)

TF.lst <- c('CEBPB','CEBPD')
gene.lst <- list()
for(i in 1:2){
    tf=TF.lst[i]
    sub.regulon2Targets <- regulon2Targets[regulon2Targets$TF %in% tf & regulon2Targets$gene %in% dnGene & regulon2Targets$highConfAnnot =='TRUE',]
    testgenes <- sub.regulon2Targets$gene
    gene.lst[[i]] <- testgenes
}


pbmc <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype.Rdata')
d <- length(intersect(keys(org.Hs.eg.db, keytype = "SYMBOL"),rownames(pbmc)))

upPval <- c()
for(j in 1:length(Aging.lst)){
    Aging.Set <- Aging.lst[[j]]
    upval <- c()
    for(i in 1:length(gene.lst)){
        testgenes <- gene.lst[[i]]
        #a: DEG in aging.genes b: aging.genes, c: DEG, d: bg gene
        a <- length(intersect(Aging.Set,testgenes));b <- length(Aging.Set);c <- length(testgenes)
        p <- fisher.test(matrix(c(a,b,c,d), nrow=2), alternative="greater")$p.value
        #upval[i] <- -log(p,10)
        upval[i] <- p
    }
    upPval <- rbind(upPval,upval)
}

colnames(upPval) <- c('CEBPB','CEBPD');rownames(upPval) <- names(Aging.lst)
print(upPval)
upPval <- -log(upPval,10)

pvalue.df <- data.frame(TF=colnames(upPval),GenAge=upPval[1,],AgingAtlas=upPval[2,])

library(reshape2)
library(ggplot2)
library(RColorBrewer)
col=brewer.pal(n = 12, name ='Dark2')

pvalue_meltdf <- reshape2::melt(pvalue.df)

p <- ggplot(data=pvalue_meltdf, mapping=aes(x=TF,y=value,fill=variable))+
    geom_bar(stat="identity",width=0.8,position='dodge')+theme_classic()+
    theme(axis.text = element_text(size = 7,colour = 'black'),axis.title = element_text(size = 7,colour = 'black'),
          legend.key.height = unit(0.3,'cm'),legend.key.width = unit(0.3,'cm'),
          legend.title=element_text(size=5),legend.text=element_text(size=7),
          plot.margin = unit(c(0.1,0,0,0.1),'cm'),legend.position="bottom",legend.direction = "vertical")+
    scale_fill_manual(values=col[6:7])+
    labs(x="", y="-log(adjPvalue)")+geom_hline(aes(yintercept=-log(0.05,10)),colour="#990000", linetype="dashed")
ggsave('JingMA_NEW/res/compControl/ChildrenvsAdults/FIG/sFig5H_CEBP_FC.pdf',
       p,width = 3.5,height = 8,units = 'cm')



##################
## SFig5j GO enrichments for CEBPB/CEBPD,  CEBPB and CEBPD
##################
tf='CEBPD'
sub.regulon2Targets <- regulon2Targets[regulon2Targets$TF %in% tf & regulon2Targets$gene %in% dnGene & regulon2Targets$highConfAnnot =='TRUE',]
upGene <- sub.regulon2Targets$gene

testgenes<- upGene
symbol2id=mapIds(org.Hs.eg.db,testgenes,"ENTREZID",'SYMBOL')
id=symbol2id[which(symbol2id!='')] #提取出非NA的ENTREZID

if(length(id) > 0){
    #GO BP 富集分析#
    ego <- enrichGO(OrgDb="org.Hs.eg.db", gene = id, ont = "BP", pvalueCutoff = 0.05, readable= TRUE) #GO富集分析
    ego <- clusterProfiler::simplify(ego,cutoff=0.8,by="p.adjust",select_fun=min,measure="Wang")
    ego_res <- as.data.frame(ego)
    if(nrow(ego_res)>0){
        write.xlsx(ego_res,paste('JingMA_NEW/res/compControl/ChildrenvsAdults/',tf,"_BP.xlsx",sep=''),row.names = FALSE,sheetName = type,append = TRUE)
    } 
}
print(ego_res$Description)

CEBPB <- read.xlsx('JingMA_NEW/res/compControl/ChildrenvsAdults/CEBPB_BP.xlsx',sheetIndex = 1)
CEBPD <- read.xlsx('JingMA_NEW/res/compControl/ChildrenvsAdults/CEBPD_BP.xlsx',sheetIndex = 1)

print(intersect(CEBPB$Description,CEBPD$Description))

pickCEBPB <- c("cellular senescence","autophagy","extracellular matrix organization","interleukin-8 secretion",
               "response to reactive oxygen species","cell aging","negative regulation of growth")

pickCEBPD <- c("regulation of inflammatory response","extracellular matrix organization",
               "positive regulation of extrinsic apoptotic signaling pathway","cellular response to oxidative stress",
               "interleukin-8 secretion","response to reactive oxygen species","collagen metabolic process",
               "regulation of NIK/NF-kappaB signaling","negative regulation of growth")

pickTerm <- rbind(CEBPB[CEBPB$Description %in% pickCEBPB,],CEBPD[CEBPD$Description %in% pickCEBPD,])
pickTerm$TF <- c(rep('CEBPB',length(pickCEBPB)),rep('CEBPD',length(pickCEBPD)))
pickTerm$TF <- factor(pickTerm$TF,levels = c('CEBPB','CEBPD'))

pickTerm$log10Pval <- -log10(pickTerm$pvalue)
#pickTerm <- pickTerm[order(pickTerm$log10Pval,decreasing = F),]
#pickTerm$Description <- factor(pickTerm$Description, levels = pickTerm$Description)

library(reshape2)
mat.plot <- pickTerm[,c('Description','TF','log10Pval')]
mat.plot <- dcast(mat.plot,Description~TF)
mat.plot[is.na(mat.plot)] <- 0
rownames(mat.plot) <- mat.plot$Description
order <- mat.plot$Description[order(mat.plot$CEBPB,mat.plot$CEBPD)]

mat.plot.melt <- melt(mat.plot)
colnames(mat.plot.melt) <- c('Description','TF','log10Pval')
mat.plot.melt$Description <- factor(mat.plot.melt$Description,levels = order)

p <- ggplot(data = mat.plot.melt, 
            aes(x = TF, y = Description, size = log10Pval,color = TF)) + 
    geom_point(fill = 'cornsilk') + 
    scale_color_manual(breaks = c('CEBPB', 'CEBPD'),values = c("#BC80BD", "#80B1D3"))+
    scale_size_continuous(range = c(0,4),breaks = c(2,4,6)) +
    theme_classic()+
    labs(x='',y='')+
    theme(axis.text.x = element_text(size = 6,colour = "black"),
          axis.text.y = element_text(size = 6,colour = "black"),
          text = element_text(size = 5,colour = "black"),legend.text = element_text(size = 5, color = 'black'),
          legend.key.size = unit(0.05,'cm'),legend.key.height = unit(0.05,'cm'),
          axis.line = element_line(size=0.2, colour = "black"),
          plot.margin = unit(c(0,-0.1,-0.1,-0.1),units = 'cm'))+
    guides(colour = guide_legend(override.aes = list(size=1)))
ggsave('JingMA_NEW/res/compControl/ChildrenvsAdults/sFig_Barplot_CEBPB_CEBPD.pdf',p,
       height = 5, width = 9, units = 'cm')


### CEBPB and CEBPD 靶基因和在一起做富集
tf='CEBPB'
sub.regulon2Targets <- regulon2Targets[regulon2Targets$TF %in% tf & regulon2Targets$gene %in% dnGene & regulon2Targets$highConfAnnot =='TRUE',]
upGene_CEBPB <- sub.regulon2Targets$gene

tf='CEBPD'
sub.regulon2Targets <- regulon2Targets[regulon2Targets$TF %in% tf & regulon2Targets$gene %in% dnGene & regulon2Targets$highConfAnnot =='TRUE',]
upGene_CEBPD <- sub.regulon2Targets$gene

testgenes<- c(upGene_CEBPB,upGene_CEBPD)
symbol2id=mapIds(org.Hs.eg.db,testgenes,"ENTREZID",'SYMBOL')
id=symbol2id[which(symbol2id!='')] #提取出非NA的ENTREZID

if(length(id) > 0){
    #GO BP 富集分析#
    ego <- enrichGO(OrgDb="org.Hs.eg.db", gene = id, ont = "BP", pvalueCutoff = 0.05, readable= TRUE) #GO富集分析
    ego <- clusterProfiler::simplify(ego,cutoff=0.8,by="p.adjust",select_fun=min,measure="Wang")
    ego_res <- as.data.frame(ego)
    if(nrow(ego_res)>0){
        write.xlsx(ego_res,"JingMA_NEW/res/compControl/ChildrenvsAdults/CEBPB_CEBPD_BP.xlsx",row.names = FALSE,sheetName = type,append = TRUE)
    } 
}
print(ego_res$Description)

## barplot
ego_res <- read.xlsx("JingMA_NEW/res/compControl/ChildrenvsAdults/CEBPB_CEBPD_BP.xlsx",sheetIndex = 1)
Common_BP <- ego_res[ego_res$p.adjust < 0.05,]
pickterm <- c('regulation of inflammatory response','extracellular matrix organization','response to interleukin-1',
              'cellular response to oxidative stress',"cell chemotaxis",
              'negative regulation of growth','cellular response to tumor necrosis factor',
              'interleukin-8 secretion',"autophagy","process utilizing autophagic mechanism" ,"positive regulation of extrinsic apoptotic signaling pathway" ,
              "response to interleukin-6","regulation of NIK/NF-kappaB signaling","intrinsic apoptotic signaling pathway by p53 class mediator",
              "cellular senescence")

index <- which(Common_BP$Description %in% pickterm)
pickCommon_BP <- Common_BP[index,]

pickCommon <- pickCommon_BP

pickCommon$log10Pval <- -log10(pickCommon$pvalue)
pickCommon <- pickCommon[order(pickCommon$log10Pval,decreasing = F),]
pickCommon$Description <- factor(pickCommon$Description, levels = pickCommon$Description)

p <- ggplot(pickCommon, aes(x = Description, y = log10Pval)) + 
    geom_bar(stat = 'identity', color = '#B15E72', fill = '#B15E72',width = 0.8) + 
    theme_classic() + coord_flip() +
    labs(x='',y = expression(paste("-log"[10], "(adj", italic("P"), "-value)"))) +
    theme(axis.title = element_text(size = 8, colour = 'black'), 
          axis.text.y = element_text(size = 8, colour = 'black'), 
          axis.text.x = element_text(size = 8, colour = 'black'))
ggsave(paste('JingMA_NEW/res/compControl/ChildrenvsAdults/sFig_Barplot_CEBPB_CEBPD.pdf',sep=''),p,
       height = 7, width = 11, units = 'cm')


