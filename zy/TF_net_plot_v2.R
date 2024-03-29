setwd('/home/zy/my_git/bioinformatics/wgk')
.libPaths('/home/zy/R/x86_64-pc-linux-gnu-library/4.0')
library(Seurat)
library(patchwork)
.libPaths('/home/yzj/R/x86_64-pc-linux-gnu-library/4.0')
library(tidyverse)
library(SCENIC)
library(igraph)
library(ggraph)
library(tidygraph)

path.TF.net <- '/home/disk/drizzle/wgk/TF_net/'

fc.cutoff <- 0.4
path.M123 <- '/home/disk/drizzle/wgk/microtia_chon_child_M1M2M3/'
path.cutoff <- paste0(path.M123, 'cutoff_', fc.cutoff, '/')
file.go.BP <- paste0(path.cutoff, 'GO_BP_all.Rdata')
list.go.BP <- readRDS(file.go.BP)
file.go.MF <- paste0(path.cutoff, 'GO_MF_all.Rdata')
list.go.MF <- readRDS(file.go.MF)

# down
file_CSC_down <- paste0(path.TF.net, 'igraph_CSC_down.RDS')
igraph_down_2 <- readRDS(file_CSC_down)

sel.MF_SCS_down <- c('antioxidant activity',
                     'S100 protein binding',
                     'extracellular matrix structural constituent')
sel.BP_SCS_down <- c('extracellular matrix organization',
                     'cellular response to zinc ion',
                     'cartilage development')

BP_SCS_down <- list.go.BP[['CSC_Microtia_decrease']]
rownames(BP_SCS_down) <- BP_SCS_down$Description
MF_SCS_down <- list.go.MF[['CSC_Microtia_decrease']]
rownames(MF_SCS_down) <- MF_SCS_down$Description
df_GO <- rbind(BP_SCS_down[sel.BP_SCS_down, c('Description', 'geneID')],
               MF_SCS_down[sel.MF_SCS_down, c('Description', 'geneID')])



node_type <- V(igraph_down_2)$name
list_gene_GO <- list()
for (gene in node_type) {
    vec_GO <- c()
    for (term in df_GO$Description) {
        gene.set <- strsplit(df_GO$geneID[df_GO$Description == term], split = '/')[[1]]
        if (gene %in% gene.set) {
            vec_GO <- c(vec_GO, term)
        }
    }
    list_gene_GO[[gene]] <- vec_GO
}

node_group <- rep('3', length(node_type))
node_group[node_type %in% names(list_gene_GO)] <- '2'
node_group[node_type %in% c('SOX5', 'SOX8', 'DBP', 'ARID5B')] <- '1'
V(igraph_down_2)$group <- node_group

ggraph_CSC_down <- as_tbl_graph(igraph_down_2)

plot_CSC_down <-
    # ggraph(ggraph_CSC_down, layout = 'centrality',cent=page_rank) +
    ggraph(ggraph_CSC_down, layout = 'stress') +
    geom_edge_link(aes(edge_width=weight, alpha = weight),color="gray",
                   arrow = arrow(length = unit(2, 'mm')),
                   end_cap = circle(1.5, 'mm'), 
                   start_cap = circle(3, 'mm')) +
    scale_edge_width(range=c(0.5,1)) +
    scale_edge_alpha(range=c(0.2,1)) + 
    scale_size_continuous(range = c(2,15)) +
    geom_node_point(aes(size = page_rank, fill = group, alpha = FC),
                    shape=21, color = 'transparent') +
    # scale_color_manual(values = c('#0000CD', '#DC143C', 'dimgray')) +
    scale_fill_manual(values = c('#4169E1', '#FF4500', 'dimgray')) +
    # scale_alpha_manual(values = c(0.8, 1, 0.6)) +
    geom_node_text(aes(filter = group == 1,label=name),size=5) +
    geom_node_text(aes(filter = group == 2,label=name),size=4, repel = T) +
    theme_void() + theme(legend.position = 'none')
ggsave(plot = plot_CSC_down, path = path.TF.net,
       filename = 'TF_CSC_down.png',
       height = 16, width = 16, units = 'cm')


# up
file_CSC_up <- paste0(path.TF.net, 'igraph_CSC_up.RDS')
igraph_up_2 <- readRDS(file_CSC_up)

sel.BP_SCS_up <- c('cell cycle arrest', 
                   'negative regulation of cell growth',
                   'negative regulation of stem cell differentiation',
                   'response to oxidative stress',
                   'p38MAPK cascade',
                   'I-kappaB kinase/NF-kappaB signaling',
                   'intrinsic apoptotic signaling pathway',
                   'extrinsic apoptotic signaling pathway',
                   'response to unfolded protein',
                   'regulation of RNA stability',
                   'activation of innate immune response',
                   'cellular response to tumor necrosis factor',
                   'regulation of inflammatory response', 
                   'cellular response to interleukin-1')

BP_SCS_up <- list.go.BP[['CSC_Microtia_increase']]
rownames(BP_SCS_up) <- BP_SCS_up$Description
# MF_SCS_up <- list.go.MF[['CSC_Microtia_increase']]
# rownames(MF_SCS_up) <- MF_SCS_up$Description
# df_GO <- rbind(BP_SCS_up[sel.BP_SCS_up, c('Description', 'geneID')],
#                MF_SCS_up[sel.MF_SCS_up, c('Description', 'geneID')])
df_GO_pre <- BP_SCS_up[sel.BP_SCS_up, c('Description', 'geneID')]
vec_desc <- c('Reduction of stem cell ability',
              'Response to oxidative stress',
              'NF-kappaB signaling and p38MAPK cascade',
              'Apoptotic signaling pathway',
              'Stability of protein and RNA',
              'Immune and inflammatory response')
vec_genes <- c(paste(df_GO_pre[1:3, 'geneID'], collapse = '/'),
               paste(df_GO_pre[4, 'geneID'], collapse = '/'),
               paste(df_GO_pre[5:6, 'geneID'], collapse = '/'),
               paste(df_GO_pre[7:8, 'geneID'], collapse = '/'),
               paste(df_GO_pre[9:10, 'geneID'], collapse = '/'),
               paste(df_GO_pre[11:14, 'geneID'], collapse = '/'))
df_GO <- data.frame(Description = vec_desc, geneID = vec_genes)

node_type_up <- V(igraph_up_2)$name
list_gene_GO_up <- list()
for (gene in node_type_up) {
    vec_GO <- c()
    for (term in df_GO$Description) {
        gene.set <- unique(strsplit(df_GO$geneID[df_GO$Description == term], split = '/')[[1]])
        if (gene %in% gene.set) {
            vec_GO <- c(vec_GO, term)
        }
    }
    list_gene_GO_up[[gene]] <- vec_GO
}


ILgene <- c('EGR1','HEXIM1','NFKBIZ','SLC7A2','SOCS3','TAC1','TNFAIP6')
Apop <- c('FNIP2','G0S2','MLLT11','PMAIP1')
ROSgene <- c('AXL','KDM6B','MMP2','NCOA7','NR4A2','NR4A3','RHOB')
Stability <- c('HSPA2','HSPA4L','HSPH1','NAF1')
NF_MAPK <- c('BMP2','MAP2K3','MAP3K8','REL')
Stem <- c('BTG2','IFRD1','RGS2','SOX4','ZFHX3')

# node_group <- rep('3', length(node_type_up))
# node_group[node_type_up %in% names(list_gene_GO_up)] <- '2'
node_group[node_type_up %in% c('EGR1', 'KLF10', 'JUNB', 'REL', 'EGR3',
                            'ATF3', 'HIVEP3', 'IRX2', 'EGR2', 'KLF2', 'BCL3',
                            'CEBPB', 'CEBPD', 'STAT5A')] <- '1'

node_group <- rep('8', length(node_type_up))
node_group[node_type %in% c(ILgene)] <- '1'
node_group[node_type %in% c(Apop)] <- '2'
node_group[node_type %in% c(ROSgene)] <- '3'
node_group[node_type %in% c(Stability)] <- '4'
node_group[node_type %in% c(NF_MAPK)] <- '5'
node_group[node_type %in% c(Stem)] <- '6'
node_group[node_type %in% c('EGR1', 'KLF10', 'JUNB', 'REL', 'EGR3',
                                  'ATF3', 'HIVEP3', 'IRX2', 'EGR2', 'KLF2', 'BCL3',
                                  'CEBPB', 'CEBPD', 'STAT5A')] <- '7'


V(igraph_up_2)$group <- node_group
ggraph_CSC_up <- as_tbl_graph(igraph_up_2)




library(RColorBrewer)
mycolor=brewer.pal(10,"Set3")
col=mycolor[1:6]
names(col) <- c('IL','Apop','ROS','Sta','NF','Stem')

set <- c('ROS','IL')
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




plot_CSC_up <- 
    ggraph(ggraph_CSC_up, layout = 'stress') +
    geom_edge_link(aes(edge_width=weight, alpha = weight),color="gray",
                   arrow = arrow(length = unit(2, 'mm')),
                   end_cap = circle(1.5, 'mm'), 
                   start_cap = circle(1, 'mm')) +
    scale_edge_width(range=c(0.5,1)) +
    scale_edge_alpha(range=c(0.2,1)) + 
    scale_size_continuous(range = c(2,7)) +
    geom_node_point(aes(size = page_rank, fill = group),
                    shape=21, color = 'transparent') +
    scale_fill_manual(values = c(col,'firebrick3', '#BEBEBE')) +
    geom_node_text(aes(filter = group == 1,label=name),size=1.75, repel = T) + 
    geom_node_text(aes(filter = group == 2,label=name),size=1.75, repel = T) + 
    geom_node_text(aes(filter = group == 3,label=name),size=1.75, repel = T) + 
    geom_node_text(aes(filter = group == 4,label=name),size=1.75, repel = T) + 
    geom_node_text(aes(filter = group == 5,label=name),size=1.75, repel = T) + 
    geom_node_text(aes(filter = group == 6,label=name),size=1.75, repel = T) + 
    geom_node_text(aes(filter = group == 7,label=name),size=2) + 
    theme_void() + theme(legend.position = 'none')
ggsave(filename = '/home/yzj/JingMA_NEW/res/compMicrotia/MicrotiavsNormal_inChildren/FIG/Fig5B_TF_CSC_up.pdf',plot_CSC_up,
       height = 10, width = 10, units = 'cm')






