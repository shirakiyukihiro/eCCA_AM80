# ============================================================
# Human dCCA scRNA-seq — post-CellBender import and processing
#
# Study      : Distinct cancer-associated fibroblast lineages and an
#              Am80-based stromal-reprogramming strategy in extrahepatic
#              cholangiocarcinoma (Meflin/ISLR).
# Repository : eCCA_AM80  (https://github.com/shirakiyukihiro/eCCA_AM80)
#
# Purpose    : Import CellBender-corrected .h5 matrices (rhdf5), rebuild Seurat objects, QC, and downstream processing after ambient-RNA removal.
# Figure(s)  : Human scRNA-seq panels (Fig. 5 and related supplementary figures).
# Input      : data/human_dCCA_scRNAseq/*cellbender_filtered.h5.
# Output     : Post-CellBender Seurat objects.
#
# Paths have been made relative for public release: place input files
# in ./data and write outputs to ./output (create these folders, or edit
# the paths below).
# ============================================================

#2023.12.7_reanalyze.rでcellbender前のh5ファイルは作成済
#それを元にcellbender→1回目はパラメーター設定なしで作成したが、あまりきれいではない？

#1回目：list.files("data/human_dCCA_scRNAseq",pattern = "cellbender_filtered.h5",full.names = T)->h5list
#以下1回目の場合、Read10xで直接読み込めないためのスクリプト
library(rhdf5)
h5ls(h5list[1])h5file <- H5File$new(h5list[index], mode = "r")
barcodes <- h5file[["matrix/barcodes"]][]
genes <- h5file[["matrix/features/name"]][]
gene_ids <- h5file[["matrix/features/id"]][]
data <- h5file[["matrix/data"]][]
indices <- h5file[["matrix/indices"]][]
indptr <- h5file[["matrix/indptr"]][]
shape <- h5file[["matrix/shape"]][]

# スパース行列として組み立て
mat <- sparseMatrix(
  i = indices + 1,
  p = indptr,
  x = data,
  dims = shape,
  dimnames = list(make.unique(genes, sep = "-"), barcodes)
)

#2回目 recelbender 直接読み込めるところまで変換済み
list.files("data/human_dCCA_scRNAseq",pattern = "matrix_._recelbender_filtered_seurat.h5",full.names = T)->h5list


library(Seurat)
library(hdf5r)
library(Matrix)
library(tidyverse)
library(scDblFinder)

seurat_obj_cellbender_list<-lapply(1:7, function(index){
  set.seed(1)
  mat<-Read10X_h5(h5list[index])
  # 新しい集約済み行列でSeuratオブジェクトを再作成
  seurat_obj_cellbender <- CreateSeuratObject(counts = mat,min.features = 200)
  # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
  seurat_obj_cellbender[["percent.mt"]] <- PercentageFeatureSet(seurat_obj_cellbender, pattern = "^MT-|^mt-")
  seurat_obj_cellbender[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj_cellbender, pattern = "^RP[SL}|^Rp[sl]")
  seurat_obj_cellbender[["percent.hb"]] <- PercentageFeatureSet(seurat_obj_cellbender, pattern = "^HB[^(P)]|^Hb[(^(p)]")
  
  seurat_obj_cellbender$log1p_ncount<-log1p(seurat_obj_cellbender$nCount_RNA)
  seurat_obj_cellbender$log1p_nfeat<-log1p(seurat_obj_cellbender$nFeature_RNA)
  sce<- scDblFinder(SingleCellExperiment(list(counts=mat)))
  seurat_obj_cellbender<-AddMetaData(seurat_obj_cellbender,metadata = colData(sce)[Cells(seurat_obj_cellbender),]%>%as.data.frame)
  counts <- GetAssayData(seurat_obj_cellbender, slot = "counts") 
  pct_top_20 <- apply(counts, 2, function(cell_counts) { 
    top_20_sum <- sum(sort(cell_counts, decreasing = TRUE)[1:20]) 
    total_sum <- sum(cell_counts) 
    return((top_20_sum / total_sum) * 100) }) 
  seurat_obj_cellbender$pct_counts_in_top_20_genes <- pct_top_20 
  seurat_obj_cellbender<-subset(seurat_obj_cellbender,subset=percent.hb<2&percent.ribo<20&percent.mt<10&scDblFinder.class=="singlet"&pct_counts_in_top_20_genes<50)
  seurat_obj_cellbender <- NormalizeData(seurat_obj_cellbender)
  seurat_obj_cellbender <- FindVariableFeatures(seurat_obj_cellbender)
  seurat_obj_cellbender <- ScaleData(seurat_obj_cellbender, features = rownames(seurat_obj_cellbender))
  seurat_obj_cellbender <- RunPCA(seurat_obj_cellbender)
  seurat_obj_cellbender <- RunUMAP(seurat_obj_cellbender, dims = 1:30)#seq(PCAtools::findElbowPoint(Stdev(seurat_obj_cellbender)^2)))
  seurat_obj_cellbender <- FindNeighbors(seurat_obj_cellbender,dims = 1:30)#seq(PCAtools::findElbowPoint(Stdev(seurat_obj_cellbender)^2)))
  seurat_obj_cellbender <- FindClusters(seurat_obj_cellbender, resolution = c(1,0.5,0.1))
  Idents(seurat_obj_cellbender)<-"RNA_snn_res.1"
  seurat_obj_cellbender$number<-rep(index)
  return(seurat_obj_cellbender)
})

names(seurat_obj_cellbender_list)<-h5list%>%basename%>%str_remove_all("_recelbender_filtered_seurat.h5")

saveRDS(seurat_obj_cellbender_list,"data/human_dCCA_scRNAseq/2025.06.03hBDC_seurat_obj_cellbender_list.rds")#
#seurat_obj_cellbender_list<-readRDS("data/human_dCCA_scRNAseq/2025.06.03hBDC_seurat_obj_cellbender_list.rds")




seurat_obj_cellbender_merge <- merge(seurat_obj_cellbender_list[[1]], y = seurat_obj_cellbender_list[-1])
#Check all objects present
table(seurat_obj_cellbender_merge$orig.ident)
seurat_obj_cellbender_merge <- NormalizeData(object = seurat_obj_cellbender_merge)

#Plot the nFeatures/counts/% Mito to get general idea about the quality of your data
VlnPlot(seurat_obj_cellbender_merge, features = c("nFeature_RNA", "nCount_RNA", "pct_counts_in_top_20_genes"),group.by = "number")
VlnPlot(seurat_obj_cellbender_merge, features = c("percent.hb", "percent.ribo", "percent.mt"),group.by = "number")

seurat_obj_cellbender_merge_subset<-subset(seurat_obj_cellbender_merge,subset = percent.ribo<10&percent.mt<5&nFeature_RNA>1000&nFeature_RNA<5000&pct_counts_in_top_20_genes<25&nCount_RNA<50000)
#FIND VARIABLE GENES
seurat_obj_cellbender_merge_subset<- FindVariableFeatures(object = seurat_obj_cellbender_merge_subset)
#Calculate Cell Cycle Score (S-G2M Difference)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
seurat_obj_cellbender_merge_subset<-JoinLayers(seurat_obj_cellbender_merge_subset)
seurat_obj_cellbender_merge_subset<- CellCycleScoring(seurat_obj_cellbender_merge_subset, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
seurat_obj_cellbender_merge_subset$CC.Difference <- seurat_obj_cellbender_merge_subset$S.Score - seurat_obj_cellbender_merge_subset$G2M.Score

#THIS STEP MAY TAKE A VERY LONG TIME
#Scale Data
seurat_obj_cellbender_merge_subset<- ScaleData(object = seurat_obj_cellbender_merge_subset, features = rownames(seurat_obj_cellbender_merge_subset))
#Run PCA and Determine Dimensions for 90% Variance
seurat_obj_cellbender_merge_subset <- RunPCA(object = seurat_obj_cellbender_merge_subset)
#Find Neighbors + Find CLusters (without harmony batch correction)
seurat_obj_cellbender_merge_subset <- FindNeighbors(object = seurat_obj_cellbender_merge_subset, dims = 1:30)
seurat_obj_cellbender_merge_subset <- FindClusters(object = seurat_obj_cellbender_merge_subset, resolution = 1.2)
#Run UMAP and get unlabelled cluster UMAP and violin plot (without harmony batch correction)
seurat_obj_cellbender_merge_subset <- RunUMAP(object = seurat_obj_cellbender_merge_subset, dims = 1:30)

DimPlot(object = seurat_obj_cellbender_merge_subset, reduction = "umap", label = TRUE, pt.size = 0.5)
library(harmony)
seurat_obj_cellbender_merge.harmony <- seurat_obj_cellbender_merge_subset %>% RunHarmony("number")

seurat_obj_cellbender_merge.harmony <- FindNeighbors(object = seurat_obj_cellbender_merge.harmony, dims = 1:30, reduction ="harmony")
seurat_obj_cellbender_merge.harmony <- FindClusters(object = seurat_obj_cellbender_merge.harmony, resolution = 1.2)

#Run UMAP and get unlabelled cluster UMAP and violin plot
seurat_obj_cellbender_merge.harmony <- RunUMAP(object = seurat_obj_cellbender_merge.harmony, dims = 1:30, reduction = "harmony")
DimPlot(object = seurat_obj_cellbender_merge.harmony, reduction = "umap", label = TRUE, pt.size = 0.5)
DimPlot(object = seurat_obj_cellbender_merge.harmony, reduction = "umap", label = TRUE, pt.size = 0.5,group.by = "Phase")
DimPlot(object = seurat_obj_cellbender_merge.harmony, reduction = "umap", label = TRUE, pt.size = 0.5,group.by = "number")
seurat_obj_cellbender_merge.harmony$number
Idents(object = seurat_obj_cellbender_merge.harmony) <- 'RNA_snn_res.1.2'

#MAKE CLUSTER MARKER TABLE
all_markers<- presto::wilcoxauc(seurat_obj_cellbender_merge.harmony, assay = "data",seurat_assay = "RNA")

all_markers %>% head()
presto::top_markers(all_markers, n = 20, auc_min = 0.5, pct_in_min = 20)%>%dplyr::select(matches("^[0-9]$"))
presto::top_markers(all_markers, n = 20, auc_min = 0.5, pct_in_min = 20)%>%dplyr::select(matches("^[1][0-9]$"))
presto::top_markers(all_markers, n = 20, auc_min = 0.5, pct_in_min = 20)%>%dplyr::select(matches("^[2][0-9]$"))
presto::top_markers(all_markers, n = 20, auc_min = 0.5, pct_in_min = 20)%>%dplyr::select(matches("^20$"))

FeaturePlot(object = seurat_obj_cellbender_merge.harmony, 
            features = c('CD79A', 'MS4A1','SDC1','KIT','CD14','CD68','CD3D','COL1A1','PECAM1',"KRT19","S100A8","ACTA2"), cols = c("grey", "deeppink"),label = T, reduction = "umap", pt.size = .5)

FeaturePlot(object = seurat_obj_cellbender_merge.harmony, 
            features = c('SOX10', 'S100B','MPZ','NCAM1','NGFR',"GAP43"), cols = c("grey", "deeppink"),label = T, reduction = "umap", pt.size = .5)
FeaturePlot(object = seurat_obj_cellbender_merge.harmony, 
            features = c('IL3RA', 'HDC','CD22','ITGAX','ENPP3',"CCR3"), cols = c("grey", "deeppink"),label = T, reduction = "umap", pt.size = .5)
 
#LABEL THE CLUSTERS round 1 for collapsed populations
Idents(object = seurat_obj_cellbender_merge.harmony) <- 'RNA_snn_res.1.2'
levels(seurat_obj_cellbender_merge.harmony)
current.cluster.ids <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25)
new.cluster.ids <- c("T_cell","T_cell","Endothelial_cell","T_cell","Endothelial_cell",#0-4#B_cell
                     "Schwann_cell", "T_reg","Macrophage","NK_cell","Fibroblast",#5-9#Basophil
                     "Endothelial_cell","Macrophage","T_cell","Epithelial_cell", "Schwann_cell",#10-14#Lymphatic_endothelial_cell
                     "Smooth_muscle_cell", "T_cell", "Mast_cell", "B_cell", "Fibroblast",#15-19#Smooth_muscle_cell
                     "Lymphatic_endothelial_cell","Epithelial_cell")#,"","", "",#20-24#Mast_cell
                     #"", "", "Epithelial_cell", "Cancer_cell", "Plasmacyte",#25- 
                     #"Mast_cell", "T_cell", "Schwann_cell", "Fibroblast", "Macrophage", 
                     #"Neutrophil")
names(x = new.cluster.ids) <- levels(x = seurat_obj_cellbender_merge.harmony)
seurat_obj_cellbender_merge.harmony <- RenameIdents(object = seurat_obj_cellbender_merge.harmony, new.cluster.ids)
DimPlot(object =seurat_obj_cellbender_merge.harmony, reduction = "umap", pt.size = 1, label = T)


saveRDS(seurat_obj_cellbender_merge.harmony,"data/human_dCCA_scRNAseq/2025.06.03hBDC_seurat_obj_cellbender_merge.harmony.rds")
#seurat_obj_cellbender_merge.harmony<-readRDS("data/human_dCCA_scRNAseq/2025.06.03hBDC_seurat_obj_cellbender_merge.harmony.rds")


seurat_obj_cellbender_merge.harmony_fibrosubset<-subset(seurat_obj_cellbender_merge.harmony,idents="Fibroblast")

seurat_obj_cellbender_merge.harmony_fibrosubset<- FindVariableFeatures(object = seurat_obj_cellbender_merge.harmony_fibrosubset)
#Scale Data
seurat_obj_cellbender_merge.harmony_fibrosubset<- ScaleData(object = seurat_obj_cellbender_merge.harmony_fibrosubset, features = rownames(seurat_obj_cellbender_merge.harmony_fibrosubset))
#Run PCA and Determine Dimensions for 90% Variance
seurat_obj_cellbender_merge.harmony_fibrosubset <- RunPCA(object = seurat_obj_cellbender_merge.harmony_fibrosubset)
#Find Neighbors + Find CLusters (without harmony batch correction)
seurat_obj_cellbender_merge.harmony_fibrosubset <- FindNeighbors(object = seurat_obj_cellbender_merge.harmony_fibrosubset, dims = 1:30)
seurat_obj_cellbender_merge.harmony_fibrosubset <- FindClusters(object = seurat_obj_cellbender_merge.harmony_fibrosubset, resolution = 1.5)
#Run UMAP and get unlabelled cluster UMAP and violin plot (without harmony batch correction)
seurat_obj_cellbender_merge.harmony_fibrosubset <- RunUMAP(object = seurat_obj_cellbender_merge.harmony_fibrosubset, dims = 1:30)

DimPlot(object = seurat_obj_cellbender_merge.harmony_fibrosubset, reduction = "umap", label = TRUE, pt.size = 0.5)
library(harmony)
seurat_obj_cellbender_merge.harmony_fibrosubset <- seurat_obj_cellbender_merge.harmony_fibrosubset %>% RunHarmony("number")

seurat_obj_cellbender_merge.harmony_fibrosubset <- FindNeighbors(object = seurat_obj_cellbender_merge.harmony_fibrosubset, dims = 1:3, reduction ="harmony")
seurat_obj_cellbender_merge.harmony_fibrosubset <- FindClusters(object = seurat_obj_cellbender_merge.harmony_fibrosubset, resolution = 1.5)

#Run UMAP and get unlabelled cluster UMAP and violin plot
seurat_obj_cellbender_merge.harmony_fibrosubset <- RunUMAP(object = seurat_obj_cellbender_merge.harmony_fibrosubset, dims = 1:3, reduction = "harmony")
DimPlot(object = seurat_obj_cellbender_merge.harmony_fibrosubset, reduction = "umap", label = TRUE, pt.size = 0.5)
DimPlot(object = seurat_obj_cellbender_merge.harmony_fibrosubset, reduction = "umap", label = TRUE, pt.size = 0.5,group.by = "Phase")
DimPlot(object = seurat_obj_cellbender_merge.harmony_fibrosubset, reduction = "umap", label = TRUE, pt.size = 0.5,group.by = "number")
seurat_obj_cellbender_merge.harmony_fibrosubset$number
Idents(object = seurat_obj_cellbender_merge.harmony_fibrosubset) <- 'RNA_snn_res.1.5'

#MAKE CLUSTER MARKER TABLE
all_markers<- presto::wilcoxauc(seurat_obj_cellbender_merge.harmony_fibrosubset, assay = "data",seurat_assay = "RNA")

all_markers %>% head()
presto::top_markers(all_markers, n = 20, auc_min = 0.5, pct_in_min = 60)%>%dplyr::select(matches("^[0-9]$"))
presto::top_markers(all_markers, n = 20, auc_min = 0.5, pct_in_min = 60)%>%dplyr::select(matches("^[1][0-9]$"))
presto::top_markers(all_markers, n = 20, auc_min = 0.5, pct_in_min = 60)%>%dplyr::select(matches("^[2][0-9]$"))
presto::top_markers(all_markers, n = 20, auc_min = 0.5, pct_in_min = 50)%>%dplyr::select(matches("^6$"))

FeaturePlot(object = seurat_obj_cellbender_merge.harmony_fibrosubset, 
            features = c('ISLR', 'IL6','PI16','DPT','ACTA2','CXCL12','PDGFRA','COL1A1','PECAM1',"CD74","THY1","HLA-DRA"), cols = c("grey", "deeppink"),label = T, reduction = "umap", pt.size = .5)
VlnPlot(object = seurat_obj_cellbender_merge.harmony_fibrosubset, 
            features = c('ISLR', 'IL6','PI16','DPT','ACTA2','CXCL12','PDGFRA','COL1A1','PECAM1',"CD74","THY1","HLA-DRA"), 
        pt.size = .5)


FeaturePlot(object = seurat_obj_cellbender_merge.harmony_fibrosubset, 
            c("nCount_RNA","nFeature_RNA","percent.mt","percent.ribo","percent.hb","log1p_ncount","log1p_nfeat",
              "scDblFinder.score","pct_counts_in_top_20_genes"))
FeaturePlot(object = seurat_obj_cellbender_merge.harmony, 
            c("nCount_RNA","nFeature_RNA","percent.mt","percent.ribo","percent.hb","log1p_ncount","log1p_nfeat",
              "scDblFinder.score","pct_counts_in_top_20_genes"))

FeaturePlot(object = seurat_obj_cellbender_merge.harmony_fibrosubset, 
            features = c('ISLR', 'IL6','PI16','DPT','ACTA2','CXCL12','PDGFRA','COL1A1','C3',"CD74","THY1","HLA-DRA"),
            cols = c("grey", "deeppink"),label = T, reduction = "umap", pt.size = .5)

FeaturePlot(object = seurat_obj_cellbender_merge.harmony_fibrosubset, 
            features = c('POSTN', 'MMP11','CXCL1','CCL2','COL15A1','DPP4','COL4A1','HSPG2','HAS1',"LRRC15","CD34","GREM1"),
            cols = c("grey", "deeppink"),label = T, reduction = "umap", pt.size = 2)


VlnPlot(object = seurat_obj_cellbender_merge.harmony_fibrosubset, 
        features = c('ISLR', 'IL6','PI16','DPT','ACTA2','CXCL12','PDGFRA','COL1A1','C3',"CD74","THY1","HLA-DRA"), 
        pt.size = .5)

Idents(object = seurat_obj_cellbender_merge.harmony_fibrosubset) <- 'RNA_snn_res.1.5'

seurat_obj_cellbender_merge.harmony_fibrosubset <- RenameIdents(seurat_obj_cellbender_merge.harmony_fibrosubset, 
                                                    `0` = "myCAF", #
                                                    `1` = "myCAF", #
                                                    `2` = "Thy1_ssl", 
                                                    `3` = "Meflin_ssl",
                                                    `4` = "iCAF", #
                                                    `5` = "Thy1_ssl", 
                                                    `6` = "iCAF",
                                                    `7` = "iCAF")
DimPlot(seurat_obj_cellbender_merge.harmony_fibrosubset, label = TRUE)
seurat_obj_cellbender_merge.harmony_fibrosubset$annotations<-Idents(seurat_obj_cellbender_merge.harmony_fibrosubset)
#seurat_obj_cellbender_merge.harmony_fibrosubset$annotations->Idents(seurat_obj_cellbender_merge.harmony_fibrosubset)

DimPlot(seurat_obj_cellbender_merge.harmony_fibrosubset, label = TRUE,group.by = "annotations")


saveRDS(seurat_obj_cellbender_merge.harmony_fibrosubset,"data/human_dCCA_scRNAseq/2025.06.05_seurat_obj_cellbender_merge.harmony_fibrosubset.rds")
#seurat_obj_cellbender_merge.harmony_fibrosubset<-readRDS("data/human_dCCA_scRNAseq/2025.06.05_seurat_obj_cellbender_merge.harmony_fibrosubset.rds")


#Put these new cluster labels as metadata
seurat_obj_cellbender_merge.harmony$annotations <- Idents(object = seurat_obj_cellbender_merge.harmony)%>%as.character
seurat_obj_cellbender_merge.harmony$annotations[Cells(seurat_obj_cellbender_merge.harmony_fibrosubset)]<-seurat_obj_cellbender_merge.harmony_fibrosubset$annotations%>%as.character
seurat_obj_cellbender_merge.harmony$annotations%>%unique
seurat_obj_cellbender_merge.harmony$annotations<-factor(seurat_obj_cellbender_merge.harmony$annotations,levels = 
                                                          c("Epithelial_cell","Macrophage","T_cell","T_reg","NK_cell","B_cell",
                                                            "Mast_cell","Endothelial_cell","Lymphatic_endothelial_cell",
                                                            "Smooth_muscle_cell",
                                                            "Schwann_cell",
                                                            "iCAF","myCAF","Meflin_ssl","Thy1_ssl"))

                                                            
DimPlot(seurat_obj_cellbender_merge.harmony, label = TRUE,group.by = "annotations")


library(CellChat)
library(patchwork)
cellchat <- createCellChat(object = GetAssayData(seurat_obj_cellbender_merge.harmony, assay = "RNA", layer = "data"), 
                           meta = data.frame(annotations=seurat_obj_cellbender_merge.harmony$annotations),
                            group.by = "annotations")
CellChatDB <-CellChatDB.human  #  CellChatDB.mouse
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
#cellchat@data.signaling

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
netVisual_heatmap(cellchat, measure = "weight")

library(NMF)
library(ggalluvial)
selectK(cellchat, pattern = "outgoing")#グラフで突然低下する前の値
nPatterns = 7
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
netAnalysis_river(cellchat, pattern = "outgoing")
netAnalysis_dot(cellchat, pattern = "outgoing")


selectK(cellchat, pattern = "incoming")
nPatterns = 7
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
netAnalysis_river(cellchat, pattern = "incoming")
netAnalysis_dot(cellchat, pattern = "incoming")

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
cellchat@netP$pathways

subsetCommunication(cellchat)%>%dplyr::filter(source%>%str_detect("ssl"))%>%arrange(-prob)%>%pull(pathway_name)%>%unique
netAnalysis_signalingRole_network(cellchat, signaling = "LAMININ", width = 8, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(cellchat, signaling = "FN1", width = 8, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(cellchat, signaling = "APP", width = 8, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(cellchat, signaling = "ANGPTL", width = 8, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(cellchat, signaling = "COLLAGEN", width = 8, height = 2.5, font.size = 10)

# 1. 通信データを取得
df <- subsetCommunication(cellchat)
# 2. 比較したい2つの送信元（source）→受信先（target）を指定
# 例: Fibroblast と Endothelial が target: Macrophage に信号送っているケース
df_Meflin <- df[df$source == "Meflin_ssl", ]
df_Thy1 <- df[df$source == "Thy1_ssl", ]

# 3. 各pathwayごとの通信強度（prob）を合計
library(dplyr)
df_Meflin_sum <- df_Meflin %>%
  group_by(pathway_name) %>%
  summarise(prob_Meflin = sum(prob))

df_Thy1_sum <- df_Thy1 %>%
  group_by(pathway_name) %>%
  summarise(prob_Thy1 = sum(prob))

# 4. 2つを統合して差を計算
df_diff <- full_join(df_Meflin_sum, df_Thy1_sum, by = "pathway_name") %>%
  mutate(
    prob_Meflin = replace_na(prob_Meflin, 0),
    prob_B = replace_na(prob_Thy1, 0),
    diff = prob_Meflin - prob_Thy1,
    abs_diff = abs(diff)
  ) %>%
  arrange(desc(abs_diff))

# 5. 上位10個の経路を表示
head(df_diff, 20)



subsetCommunication(cellchat)%>%dplyr::filter(source%>%str_detect("ssl"),pathway_name=="COLLAGEN",ligand=="COL4A1" )%>%arrange(-prob)
subsetCommunication(cellchat)%>%dplyr::filter(source%>%str_detect("ssl"),pathway_name=="LAMININ" )%>%arrange(-prob)

netVisual_bubble(cellchat, sources.use = c(14,15), targets.use = 1, signaling = c("COLLAGEN","LAMININ"), remove.isolate = FALSE)

df%>% filter(source == "Meflin_ssl", target == "Epithelial_cell", pathway_name %in% c("COLLAGEN", "LAMININ")) %>%
  select(interaction_name_2, pathway_name, ligand, receptor, prob)%>%left_join(
    df%>% filter(source == "Thy1_ssl", target == "Epithelial_cell", pathway_name %in% c("COLLAGEN", "LAMININ")) %>%
  select(interaction_name_2, pathway_name, ligand, receptor, prob), by = "interaction_name_2")%>%
  mutate(diff = prob.x - prob.y)%>%filter(diff>0)%>%pull(interaction_name_2)

df%>% filter(source == "Meflin_ssl", target == "Epithelial_cell", pathway_name %in% c("COLLAGEN", "LAMININ")) %>%
  select(interaction_name_2, pathway_name, ligand, receptor, prob)%>%left_join(
    df%>% filter(source == "Thy1_ssl", target == "Epithelial_cell", pathway_name %in% c("COLLAGEN", "LAMININ")) %>%
      select(interaction_name_2, pathway_name, ligand, receptor, prob), by = "interaction_name_2")%>%
  mutate(diff = prob.x - prob.y)%>%filter(diff<0)%>%pull(interaction_name_2)


#Meflinは基底膜系、Thy1は通常のECMとのシグナルに関わることがわかった。EMTに基底膜系は関わるのでそれで論じることはできるかも

saveRDS(cellchat,"data/human_dCCA_scRNAseq/2025.06.05_cellchat.rds")
#seurat_obj_cellbender_merge.harmony_fibrosubset<-readRDS("data/human_dCCA_scRNAseq/2025.06.05_cellchat.rds")











cellchat <- computeNetSimilarity(cellchat, type = "functional")
#reticulate::py_install(packages = 'umap-learn')
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)

cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "structural")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "structural", label.size = 3.5)
netVisual_embeddingZoomIn(cellchat, type = "structural", nCol = 2)




