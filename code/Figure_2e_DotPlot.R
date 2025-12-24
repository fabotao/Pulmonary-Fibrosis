library(Seurat)
library(ggplot2)

# --- Load Input Datasets ------------------------------------------------------
# All data files should be placed in the appropriate subfolders.
# Example folder structure:
#   ./processed_data/           - for main input files
#   ./output/         - for saving plots and results


Neu.obj <- readRDS('processed_data/Neu_obj.RDS')
Neu.obj.PBS.CIP <- subset(Neu.obj, subset=stim %in% c('PBS','CIP'))


Neu.obj.PBS.CIP$subtype_mod <- c('Aged','Mature','Mature','Immature','Mature','Mature','Mature','Immature','Aged','Pro/pre','Mature','Mature',
                                 'Mature','Mature','Aged','Aged','Pro/pre','Aged','Mature','Pro/pre','Pro/pre','Pro/pre')[match(Neu.obj.PBS.CIP$seurat_clusters, c(0:21))]


Neu.obj.PBS.CIP$stim <- factor(Neu.obj.PBS.CIP$stim, levels = c('PBS', 'CIP'), ordered = T)
Neu.obj.PBS.CIP$subtype_mod <- factor(Neu.obj.PBS.CIP$subtype_mod, levels=c('Pro/pre','Immature','Mature','Aged'), ordered = T)



## N1 vs N2 marker dotplot
N2 <- c('Ly6C','Mmp9','S100a8','S100a9','Bhlhe40') #,'Ccl2','Ccl5', 'Il10', 'Fatp2','Ccl17','Cxcl5','Gpx4','Il6','Tgfb2','Mettl1','Rab27b')
N1 <- c('Cd11b','Cxcr2','Lcn2','Tnfaip8l2','Tnfaip2', # 'Tnfaip8l1','Tnfaip1','Mpo',
        'Tnfaip6','Tnfaip8','Tnfaip3','Itgam','Nos2','Hif1a','Stat3','Prdx1','Vegfa','Tlr4')
N17 <- c('Il21','Il17a','Rorc')

Idents(Neu.obj.PBS.CIP) = 'subtype_mod'
p.N1.N2.dotplot <- DotPlot(Neu.obj.PBS.CIP, features = c(N1, N2, N17), cols=c('#003366','#D32F2F')) + coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ggsave(p.N1.N2.dotplot, filename = 'output/Neu_N1_N2_markers_dotplot.pdf', width = 5, height = 6)

