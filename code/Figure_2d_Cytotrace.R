library(Seurat)
library(magrittr)
library(SingleR)
library(ggplot2)
library(dplyr)
library(org.Mm.eg.db)
library(ggrepel)
library(patchwork)
library(stringi)

# --- Load Input Datasets ------------------------------------------------------
# All data files should be placed in the appropriate subfolders.
# Example folder structure:
#   ./processed_data/           - for main input files
#   ./output/         - for saving plots and results


cols <- c("#532C8A","#c19f70","#f9decf","#c9a997","#B51D8D","#9e6762","#3F84AA","#F397C0",
          "#C594BF","#DFCDE4","#eda450","#635547","#C72228","#EF4E22","#f77b59","#989898",
          "#7F6874","#8870ad","#65A83E","#EF5A9D","#647a4f","#FBBE92","#354E23","#139992",
          "#C3C388","#8EC792","#0F4A9C","#8DB5CE","#1A1A1A","#FACB12","#C9EBFB","#DABE99",
          "#ed8f84","#005579","#CDE088","#BBDCA8","#F6BFCB"
)

cols.update <- c('#A4CDE1','#3186BD','#82BF98','#6CBE59','#5EB9A9','#F9B168',
                 '#F68C26','#DD9D82','#F27C7C','#DC432E','#9B7BB7','#927E92',
                 '#EFE584','#C17E77')



## Markers for new subtype
neu.subtype.markers <- c('Cd34','Sox4','Rpl12',
                         'Irf8','Gata1','Gata2','Hexa','Cebpa','Prtn3','Prss57','Cstg','Mpo','Elane','Ctsc','Sod2','P2rx4','Sod1','Spr','Klrk1','C1qbp','Pla2g6',
                         'Runx1','Gfi1','Cebpe','Hexb','Stom','Ltf','Camp','Cybb','Ceacam1','Lcn2','Lyz2','Ltb4r1','Fcnb','Itgb3','Anxa1',
                         'Hp','Mmp25','Mmp9','Mmp8','Cfp','Adam8','Slc11a1','Lyst','Ptk2b','Trem1','Ccl3','Ccl4','Ager','Cxcr3','Syk','Abca1','Cd300a',
                         'Cd81','Pvr','Itgb1','Spn','Cd47','Tfrc','Itga4','Itga5','Cxcr4','Cd24','Ly75','Itgam','Cd88','Cd55','Itgal')
# Cd155: Pvr, Cd29: Itgb1, Cd43: Spn, Cd71: Tfrc, Cd49d: Itga4, Cd49e: Itga5, Cd184: Cxcr4, Cd205: Ly75, Cd11b: Itgam, Cd11a: Itgal

neu.mk1 <- c('Cd34','Sox4','Rpl12','Elane','Cebpa','Prtn3','Mpo','Chil3','Tuba1b','Fcnb','Ltf','Ngp','Camp','Retnlg','Cxcl2','Mmp8','Ccl6','Gm5483','Stfa2l1','Isg15','Rsad2','Ifit3','Fgl2','Gm2a','Gngt2')

Neu.obj <- readRDS('processed_data/Neu_obj.RDS')
Neu.obj.PBS.CIP <- subset(Neu.obj, subset=stim %in% c('PBS','CIP'))

DotPlot(Neu.obj.PBS.CIP, features = neu.mk1, cluster.idents = T) + coord_flip()

Neu.obj.PBS.CIP$subtype_mod <- c('Aged','Mature','Mature','Immature','Mature','Mature','Mature','Immature','Aged','Pro/pre','Mature','Mature',
                     'Mature','Mature','Aged','Aged','Pro/pre','Aged','Mature','Pro/pre','Pro/pre','Pro/pre')[match(Neu.obj.PBS.CIP$seurat_clusters, c(0:21))]

neu.subtype.cnt <- table(Neu.obj.PBS.CIP$subtype_mod, Neu.obj.PBS.CIP$stim, Neu.obj.PBS.CIP$time)
write.csv(neu.subtype.cnt, file = 'output/Neutrophil_subtype_cnt.csv')


# cols.subset <- cols[c(1,2,3,5,6)]
cols.subset <- cols.update[1:5]
types <- c('Pro/pre',"Immature","Mature","Aged")
names(cols.subset) <- types
Neu.obj.PBS.CIP$stim <- factor(Neu.obj.PBS.CIP$stim, levels = c('PBS', 'CIP'), ordered = T)
Neu.obj.PBS.CIP$subtype_mod <- factor(Neu.obj.PBS.CIP$subtype_mod, levels=c('Pro/pre','Immature','Mature','Aged'), ordered = T)

p.neu <- DimPlot(Neu.obj.PBS.CIP, group.by = 'subtype_mod', cols = cols.subset)
ggsave(p.neu, filename = 'output/Neutrophil_dimplot.pdf', width = 5, height = 4)

p.neu.split <- DimPlot(Neu.obj.PBS.CIP, group.by = 'subtype_mod', cols = cols.subset, split.by = 'stim')
ggsave(p.neu.split, filename = 'output/Neutrophil_dimplot_splitby_stim.pdf', width = 10, height = 5)

Idents(Neu.obj.PBS.CIP) <- 'subtype_mod'
neu.markers <- FindAllMarkers(Neu.obj.PBS.CIP)

top20 <- neu.markers %>% group_by(cluster) %>% top_n(n=10,wt=avg_log2FC)
VlnPlot(Neu.obj.PBS.CIP, features = top20$gene, stack = T, flip = T)








## Use Cytotrace
library(CytoTRACE)
exp1 <- as.matrix(Neu.obj.PBS.CIP@assays$RNA@counts)
exp1 <- exp1[apply(exp1 > 0,1,sum) >= 5,]
results <- CytoTRACE(exp1,ncores = 1)
phenot <- Neu.obj.PBS.CIP$subtype_mod
phenot <- as.character(phenot)
names(phenot) <- rownames(Neu.obj.PBS.CIP@meta.data)
emb <- Neu.obj.PBS.CIP@reductions[["umap"]]@cell.embeddings
plotCytoTRACE(results, phenotype = phenot, emb = emb, outputDir = './')
plotCytoGenes(results, numOfGenes = 30, outputDir = './')

plot(Neu.obj.PBS.CIP@reductions$umap@cell.embeddings[,1],
     Neu.obj.PBS.CIP@reductions$umap@cell.embeddings[,2],
     col=ifelse(Neu.obj.PBS.CIP$seurat_clusters %in% c(2), 'red', 'grey'), cex=0.2)





