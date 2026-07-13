# ============================================================
# Human dCCA scRNA-seq — cell-type annotation (exploratory)
#
# Study      : Distinct cancer-associated fibroblast lineages and an
#              Am80-based stromal-reprogramming strategy in extrahepatic
#              cholangiocarcinoma (Meflin/ISLR).
# Repository : eCCA_AM80  (https://github.com/shirakiyukihiro/eCCA_AM80)
#
# Purpose    : Load the human distal cholangiocarcinoma (dCCA) Seurat object and annotate cell types with canonical markers (feature/violin plots).
# Figure(s)  : Human scRNA-seq panels (Fig. 5 and related supplementary figures).
# Input      : data/human_dCCA_scRNAseq/ (Seurat .rds).
# Output     : Annotation feature/violin plots.
#
# Paths have been made relative for public release: place input files
# in ./data and write outputs to ./output (create these folders, or edit
# the paths below).
# ============================================================

library(tidyverse)
path<-"data/human_dCCA_scRNAseq"
list.files(path)
library(Seurat)

sc_hdBDC<-readRDS(path%>%paste0("/2023.11.10_seurat.rds"))#2023.11.10

DimPlot(sc_hdBDC,split.by = "sample")

FeaturePlot(object = sc_hdBDC, features = c('ACTA2', 'CDH11', 'PDGFRB', 'COL1A1', 'COL3A1', 'RGS5', 'IGFBP7',
                                            'PDPN','DCN','MCAM','IL6','APOE','GLI1','GLI2','GLI3','PDGFA'), 
            cols = c("grey", "deeppink"), reduction = "umap", pt.size = .5)

VlnPlot(sc_hdBDC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)#
VlnPlot(sc_hdBDC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)#

#fibro 7 pericyte 14
FeaturePlot(object = sc_hdBDC, features = c('CD3E', 'CD19', 'CD68', 'CD14', 'PECAM1', 'TEK'),label=T,
                                            cols = c("grey", "deeppink"), reduction = "umap", pt.size = .5)
#T:6,1,9,0 myeloid:12,8 endo:3,4 epi:13,2,5,16  
FeaturePlot(object = sc_hdBDC, features = c('SDC1', 'CD79A', 'LYVE1', 'KRT19', 'KRT7', 'KRT14'),label=T,
            cols = c("grey", "deeppink"), reduction = "umap", pt.size = .5)
#B:10
FeaturePlot(object = sc_hdBDC, features = c('MSLN', 'LRRN4', 'UPK3B', 'NEFL', 'NEFH', 'SYP'),label=T,
            cols = c("grey", "deeppink"), reduction = "umap", pt.size = .5)

FeaturePlot(object = sc_hdBDC, features = c('KIT', 'GFAP', 'HLA-DRB1', 'CD34', 'FCGR3A', 'SOX10'),label=T,
            cols = c("grey", "deeppink"), reduction = "umap", pt.size = .5)
#mast:15 schwan 11
FeaturePlot(object = sc_hdBDC, features = c('FCGR2A', 'NGFR', 'GAP43', 'NCAM1', 'CD63', 'ENPP3'),label=T,
            cols = c("grey", "deeppink"), reduction = "umap", pt.size = .5)
#mast:15 schwan 11
FeaturePlot(object = sc_hdBDC, features = c('FCER1A', 'ITGA4', 'ITGB7', 'VCAM1', 'KITLG', 'FCGR2B'),label=T,
            cols = c("grey", "deeppink"), reduction = "umap", pt.size = .5)


sc_hdBDC$cluster<-Idents(sc_hdBDC)#sc_hdBDC$cluster->Idents(sc_hdBDC)
sc_hdBDC <- RenameIdents(sc_hdBDC, 
                      `0` = "T_cell", `1` = "T_cell", 
                      `2` = "epithelial_cell", `3` = "endothelial_cell", 
                      `4` = "endothelial_cell",`5` = "epithelial_cell", 
                      `6` = "T_cell", `7` = "Fibroblast",
                      `8` = "Myeloid_cell", `9` = "T_cell",
                      `10` = "B_cell",`11` = "Schwann_cell",
                      `12` = "Myeloid_cell",`13` = "epithelial_cell",
                      `14` = "Pericyte",`15` = "Mast_cell",
                      `16` = "epithelial_cell")
DimPlot(sc_hdBDC, label = TRUE)
sc_hdBDC$cell_type<-Idents(sc_hdBDC)#sc_hdBDC$cell_type->Idents(sc_hdBDC)

DotPlot(sc_hdBDC,features = c("CD3E","KRT19","PECAM1","COL1A1","ISLR","THY1","ACTA2","CD14","CD68","CD79A",
                   "SOX10","MCAM","KIT"))+RotatedAxis()



