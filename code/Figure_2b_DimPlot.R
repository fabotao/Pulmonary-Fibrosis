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
cols.subset <- cols.update[1:14]



## Markers for new subtype
neu.subtype.markers <- c('Cd34','Sox4','Rpl12',
                         'Irf8','Gata1','Gata2','Hexa','Cebpa','Prtn3','Prss57','Cstg','Mpo','Elane','Ctsc','Sod2','P2rx4','Sod1','Spr','Klrk1','C1qbp','Pla2g6',
                         'Runx1','Gfi1','Cebpe','Hexb','Stom','Ltf','Camp','Cybb','Ceacam1','Lcn2','Lyz2','Ltb4r1','Fcnb','Itgb3','Anxa1',
                         'Hp','Mmp25','Mmp9','Mmp8','Cfp','Adam8','Slc11a1','Lyst','Ptk2b','Trem1','Ccl3','Ccl4','Ager','Cxcr3','Syk','Abca1','Cd300a',
                         'Cd81','Pvr','Itgb1','Spn','Cd47','Tfrc','Itga4','Itga5','Cxcr4','Cd24','Ly75','Itgam','Cd88','Cd55','Itgal')
# Cd155: Pvr, Cd29: Itgb1, Cd43: Spn, Cd71: Tfrc, Cd49d: Itga4, Cd49e: Itga5, Cd184: Cxcr4, Cd205: Ly75, Cd11b: Itgam, Cd11a: Itgal

neu.mk1 <- c('Cd34','Sox4','Rpl12','Elane','Cebpa','Prtn3','Mpo','Chil3','Tuba1b','Fcnb','Ltf','Ngp','Camp','Retnlg','Cxcl2','Mmp8','Ccl6','Gm5483','Stfa2l1','Isg15','Rsad2','Ifit3','Fgl2','Gm2a','Gngt2')

Neu.obj <- readRDS('processed_data/Neu_obj.RDS')


DotPlot(Neu.obj, features = neu.mk1, cluster.idents = T) + coord_flip()

Neu.obj$subtype_mod <- c('Aged','Mature','Mature','Immature','Mature','Mature','Mature','Immature','Aged','Pro/pre','Mature','Mature',
                         'Mature','Mature','Aged','Aged','Pro/pre','Aged','Mature','Pro/pre','Pro/pre','Pro/pre')[match(Neu.obj$seurat_clusters, c(0:21))]
Neu.obj$subtype_mod <- factor(Neu.obj$subtype_mod, levels = c('Pro/pre', 'Immature', 'Mature', 'Aged'), ordered = T)
Neu.obj$condition <- ifelse(Neu.obj$stim=='PBS', 'PBS', paste0(Neu.obj$stim,'_', Neu.obj$time))

cols <- c('#bea2a6', '#6fcbae', '#f37d80', '#8473b5')
names(cols) <- c('Pro/pre', 'Immature', 'Mature', 'Aged')


p.dim <- DimPlot(Neu.obj, group.by = 'subtype_mod', split.by = 'condition', raster = F, cols = cols) + labs(title=NULL)
ggsave(p.dim, filename = 'output/Neutrophil_Dimplot_by_condition.pdf', height = 5, width = 21)


neu.subtype.cnt <- table(Neu.obj$subtype_mod, Neu.obj$condition)
write.csv(neu.subtype.cnt, file = 'output/Neu_condition_count.csv')




##
celltype.unique <- unique(Neu.obj$subtype_mod)


marker.PBS.D10.BLM <- do.call(rbind,lapply(celltype.unique, function(x){

  sce.obj <- Neu.obj
  sce.obj$group <- ifelse(sce.obj$condition == 'BLM_D10', 'A',
                          ifelse(sce.obj$condition=='PBS' , 'B', 'C'))
  Idents(sce.obj) <- 'group'
  mk <- FindMarkers(sce.obj, ident.1 = 'A', ident.2 = 'B')
  mk$gene <- rownames(mk)
  mk$cluster = x
  return(mk)
}))


marker.PBS.D10.CIP <- do.call(rbind,lapply(celltype.unique, function(x){

  sce.obj <- Neu.obj
  sce.obj$group <- ifelse(sce.obj$condition=='CIP_D10' , 'A',
                          ifelse(sce.obj$condition=='PBS', 'B', 'C'))
  Idents(sce.obj) <- 'group'
  mk <- FindMarkers(sce.obj, ident.1 = 'A', ident.2 = 'B')
  mk$gene <- rownames(mk)
  mk$cluster = x
  return(mk)
}))


marker.PBS.D28.BLM <- do.call(rbind,lapply(celltype.unique, function(x){

  sce.obj <- Neu.obj
  sce.obj$group <- ifelse(sce.obj$condition=='BLM_D28', 'A',
                          ifelse(sce.obj$condition=='PBS', 'B', 'C'))
  Idents(sce.obj) <- 'group'
  mk <- FindMarkers(sce.obj, ident.1 = 'A', ident.2 = 'B')
  mk$gene <- rownames(mk)
  mk$cluster = x
  return(mk)
}))

marker.PBS.D28.CIP <- do.call(rbind,lapply(celltype.unique, function(x){

  sce.obj <- Neu.obj
  sce.obj$group <- ifelse(sce.obj$condition=='CIP_D28', 'A',
                          ifelse(sce.obj$condition=='PBS', 'B', 'C'))
  Idents(sce.obj) <- 'group'
  mk <- FindMarkers(sce.obj, ident.1 = 'A', ident.2 = 'B')
  mk$gene <- rownames(mk)
  mk$cluster = x
  return(mk)
}))


marker.D10.BLM.CIP <- do.call(rbind,lapply(celltype.unique, function(x){

  sce.obj <- Neu.obj
  sce.obj$group <- ifelse(sce.obj$condition=='CIP_D10', 'A',
                          ifelse(sce.obj$condition=='BLM_D10', 'B', 'C'))
  Idents(sce.obj) <- 'group'
  mk <- FindMarkers(sce.obj, ident.1 = 'A', ident.2 = 'B')
  mk$gene <- rownames(mk)
  mk$cluster = x
  return(mk)
}))


marker.D28.BLM.CIP <- do.call(rbind,lapply(celltype.unique, function(x){

  sce.obj <- Neu.obj
  sce.obj$group <- ifelse(sce.obj$condition=='CIP_D28', 'A',
                          ifelse(sce.obj$condition=='BLM_D28', 'B', 'C'))
  Idents(sce.obj) <- 'group'
  mk <- FindMarkers(sce.obj, ident.1 = 'A', ident.2 = 'B')
  mk$gene <- rownames(mk)
  mk$cluster = x
  return(mk)
}))



write.csv(marker.PBS.D10.BLM, file='output/Neu_subtype_marker_PBS_vs_BLM_D10.csv', row.names = F)
write.csv(marker.PBS.D10.CIP, file='output/Neu_subtype_marker_PBS_vs_CIP_D10.csv', row.names = F)
write.csv(marker.PBS.D28.BLM, file='output/Neu_subtype_marker_PBS_vs_BLM_D28.csv', row.names = F)
write.csv(marker.PBS.D28.CIP, file='output/Neu_subtype_marker_PBS_vs_CIP_D28.csv', row.names = F)
write.csv(marker.D10.BLM.CIP, file='output/Neu_subtype_marker_BLM_D10_vs_CIP_D10.csv', row.names = F)
write.csv(marker.D28.BLM.CIP, file='output/Neu_subtype_marker_BLM_D28_vs_CIP_D28.csv', row.names = F)





