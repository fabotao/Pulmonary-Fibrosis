library(Seurat)
library(magrittr)
#library(SingleR)
library(ggplot2)
library(dplyr)
library(org.Mm.eg.db)
library(ggpubr)


# --- Load Input Datasets ------------------------------------------------------
# All data files should be placed in the appropriate subfolders.
# Example folder structure:
#   ./processed_data/           - for main input files
#   ./output/         - for saving plots and results


## Compute cell type fraction in conditions
sce.all <- readRDS('processed_data/sce_all_celltype.RDS')

group.map <- data.frame(sample=c('C56','C58','C60','C62','C63','C64','C66','C67','C69','C70'),
                        stim=c('CIP','BLM','CIP','BLM','BLM','PBS','PBS','CIP','CIP','BLM'),
                        time=c('D28','D28','D10','D28','D10','D10','D10','D10','D28','D10'))
sce.all$stim <- group.map$stim[match(sce.all$orig.ident, group.map$sample)]
sce.all$time <- group.map$time[match(sce.all$orig.ident, group.map$sample)]



cols <- c("#532C8A","#c19f70","#f9decf","#c9a997","#B51D8D","#9e6762","#3F84AA","#F397C0",
          "#C594BF","#DFCDE4","#eda450","#635547","#C72228","#EF4E22","#f77b59","#989898",
          "#7F6874","#8870ad","#65A83E","#EF5A9D","#647a4f","#FBBE92","#354E23","#139992",
          "#C3C388","#8EC792","#0F4A9C","#8DB5CE","#1A1A1A","#FACB12","#C9EBFB","#DABE99",
          "#ed8f84","#005579","#CDE088","#BBDCA8","#F6BFCB"
)

cols.update <- c('#A4CDE1','#3186BD','#82BF98','#6CBE59','#5EB9A9','#F9B168',
                 '#F68C26','#DD9D82','#F27C7C','#DC432E','#9B7BB7','#927E92',
                 '#EFE584','#C17E77')


#cols.subset <- cols[c(1,2,3,5,6,7,8,9,11,13,14,19,20,24,25)]
cols.subset <- cols.update[1:14]
types <- c("Adipocytes","B cells","Basophils","Dendritic cells","Endothelial cells","Epithelial cells","Erythrocytes","Fibroblasts",
           "Macrophages","Monocytes","NK cells","Neutrophils","Platelets","T cells")
names(cols.subset) <- types

sce.all$condition <- ifelse(sce.all$stim=='PBS', 'PBS', paste0(sce.all$time, '_', sce.all$stim))
sce.all$condition <- factor(sce.all$condition, levels=c('PBS', 'D10_BLM', 'D10_CIP', 'D28_BLM', 'D28_CIP'), ordered = T)

p.celltype.umap.split <- DimPlot(sce.all, group.by = 'cell_type', split.by = 'condition', raster = F, cols = cols.subset) + labs(title=NULL)
ggsave(p.celltype.umap.split, filename = 'output/Celltype_umap_splitby_time.pdf', height = 5, width = 25)




get.celltype.prop <- function(time, title){
  cell.time <- sce.all@meta.data
  cell.time.sub <- cell.time[cell.time$time == time,]
  cell.cnt <- as.matrix(table(cell.time.sub$cell_type, cell.time.sub$stim))
  cell.ratio <- t(t(cell.cnt)/colSums(cell.cnt))
  if(time=='D28'){
    pbs.cnt <- table(cell.time$cell_type[cell.time$time=='D10' & cell.time$stim=='PBS'])
    pbs.ratio <- pbs.cnt/sum(pbs.cnt)
    cell.ratio <- cbind(pbs.ratio[rownames(cell.ratio)], cell.ratio)
    colnames(cell.ratio)[1] <- 'PBS'
  }
  plot.df <- data.frame(celltype=rep(rownames(cell.ratio), dim(cell.ratio)[2]),
                        ratio=c(cell.ratio[,'PBS'],
                                cell.ratio[,'BLM'],
                                cell.ratio[,'CIP']),
                        stim=c(rep('PBS', dim(cell.ratio)[1]),
                               rep('BLM', dim(cell.ratio)[1]),
                               rep('CIP', dim(cell.ratio)[1])))
  plot.df$stim <- factor(plot.df$stim, levels=c('PBS','BLM','CIP'),ordered = T)
  cols.time <- cols.subset[names(cols.subset) %in% plot.df$celltype]
  p <- ggplot(plot.df, aes(x=stim, y=ratio, fill=factor(celltype))) + geom_bar(stat="identity", color='white') +
    theme_bw() + scale_fill_manual(NULL,values = cols.subset) + labs(x='Condition', y='Proportion', title = title) +
    theme(axis.text = element_text(colour = 'black'),
          panel.grid = element_blank())
  return(p)
}

p.proportion <- ggarrange(get.celltype.prop('D10', title='Day10'), get.celltype.prop('D28', title='Day28'),
                          nrow=1, common.legend = T, legend = 'right')
ggsave(p.proportion, filename = 'output/Celltype_proportion.pdf', width = 10, height = 5)


cnt.sample.celltype <- t(table(sce.all$condition, sce.all$cell_type))
write.csv(cnt.sample.celltype, file = 'output/Count_condition_celltype.csv')





