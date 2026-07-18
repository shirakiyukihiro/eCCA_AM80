# ============================================================
# Human dCCA scRNA-seq — figure plotting
#
# Study      : Distinct cancer-associated fibroblast lineages and an
#              Am80-based stromal-reprogramming strategy in extrahepatic
#              cholangiocarcinoma (Meflin/ISLR).
# Repository : eCCA_AM80  (https://github.com/shirakiyukihiro/eCCA_AM80)
#
# Purpose    : Load the Harmony-integrated human dCCA Seurat object and the fibroblast subset, and generate the manuscript UMAP / feature / violin figures.
# Figure(s)  : Human scRNA-seq panels (Figure 5F-H and related supplementary figures).
# Input      : data/human_dCCA_scRNAseq/*.harmony.rds and *fibrosubset.rds.
# Output     : Human scRNA-seq figure panels.
#
# Paths have been made relative for public release: place input files
# in ./data and write outputs to ./output (create these folders, or edit
# the paths below).
# ============================================================

sessionInfo()
library(Seurat)
list.files("data")
basepath = "data"
list.files(paste0(basepath,"/human_dCCA_scRNAseq"))
h_dBDC_path ="/human_dCCA_scRNAseq"
seurat_obj_cellbender_merge.harmony<-readRDS(paste0(basepath,h_dBDC_path,"/2025.06.03hBDC_seurat_obj_cellbender_merge.harmony.rds"))
DimPlot(object = seurat_obj_cellbender_merge.harmony, reduction = "umap", label = TRUE, pt.size = 0.5)

seurat_obj_cellbender_merge.harmony_fibrosubset<-readRDS(paste0(basepath,h_dBDC_path,"/2025.06.05_seurat_obj_cellbender_merge.harmony_fibrosubset.rds"))

DimPlot(object = seurat_obj_cellbender_merge.harmony_fibrosubset, reduction = "umap", label = TRUE, pt.size = 0.5)
DimPlot(seurat_obj_cellbender_merge.harmony_fibrosubset, label = TRUE)

DimPlot(seurat_obj_cellbender_merge.harmony_fibrosubset, label = TRUE,group.by = "annotations")

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

svglite::svglite(paste0(basepath,h_dBDC_path,"/2026.01.24_dimplot.svg"),
                 height = 3,width =3)
DimPlot(seurat_obj_cellbender_merge.harmony, label = TRUE,
        label.size = 2.7, raster = TRUE, raster.dpi = c(5000,5000), 
        pt.size = 20) + NoLegend()+
  theme(
    plot.title = element_text(size = 8, color = "black"),      # タイトル文字サイズ10pt
    axis.title = element_text(size = 8, color = "black"),       # 軸タイトル8pt
    axis.text = element_text(size = 8, color = "black"),        # 軸目盛りテキスト8pt
    legend.title = element_text(size = 8, color = "black"),     # 凡例タイトル8pt
    legend.text = element_text(size = 8, color = "black"))
dev.off()


levels(Idents(seurat_obj_cellbender_merge.harmony))
new_order <- c(
  "Fibroblast",
  "Smooth_muscle_cell",
  "Endothelial_cell",
  "Lymphatic_endothelial_cell",
  "Schwann_cell",
  "Epithelial_cell",
  "Macrophage",
  "Mast_cell",
  "NK_cell",
  "T_cell",
  "T_reg",
  "B_cell")

Idents(seurat_obj_cellbender_merge.harmony) <- factor(
  Idents(seurat_obj_cellbender_merge.harmony),levels = new_order)

# 確認
levels(Idents(seurat_obj_cellbender_merge.harmony))

features_all <- c(
  # Fibroblast
  "COL1A1","COL1A2","DCN","LUM",
  # Smooth muscle
  "ACTA2","TAGLN","MYH11",
  # Endothelial (blood)
  "PECAM1","VWF","KDR","FLT1","RAMP2",
  # Lymphatic endothelial
  "PDPN","PROX1","LYVE1","CCL21",
  # Schwann / Neural
  "SOX10","S100B","NCAM1","NGFR","GAP43",
  # Epithelial
  "EPCAM","KRT19","KRT8","KRT18",
  # Myeloid (Macrophage)
  "LST1","CD14","LILRB1","CD68","S100A8","S100A9",
  # Mast
  "KIT","TPSAB1","HDC",
  # NK
  "NKG7","GNLY",
  # T 
  "CD3D","CD3E",
  # Treg-specific
  "FOXP3","IL2RA","CTLA4","IKZF2",
  # B / Plasma
  "CD79A","MS4A1","CD19")

svglite::svglite(paste0(basepath,h_dBDC_path,"/2026.01.24_dotplot.svg"),
                 height = 3,width =8.2)
DotPlot(
  object = seurat_obj_cellbender_merge.harmony,
  features = features_all,dot.scale = 3,
  cols = c("grey90", "deeppink")
) +
  RotatedAxis() +
  theme(
    plot.title = element_text(size = 10, color = "black"),      # タイトル文字サイズ10pt
    axis.title = element_blank(),       # 軸タイトル8pt
    axis.text = element_text(size = 8, color = "black"),        # 軸目盛りテキスト8pt
    legend.title = element_text(size = 8, color = "black"),     # 凡例タイトル8pt
    legend.text = element_text(size = 8, color = "black"),
    legend.position = "top",
    legend.box = "horizontal")
dev.off()



seurat_obj_cellbender_merge.harmony_fibrosubset$annotations

svglite::svglite(paste0(basepath,h_dBDC_path,"/2026.01.24_dimplot_fibro.svg"),
                 height = 3,width =3)
DimPlot(seurat_obj_cellbender_merge.harmony_fibrosubset, label = TRUE,
        label.size = 2.7, raster = TRUE, raster.dpi = c(5000,5000), 
        pt.size = 70) + NoLegend()+
  theme(
    plot.title = element_text(size = 8, color = "black"),      # タイトル文字サイズ10pt
    axis.title = element_text(size = 8, color = "black"),       # 軸タイトル8pt
    axis.text = element_text(size = 8, color = "black"),        # 軸目盛りテキスト8pt
    legend.title = element_text(size = 8, color = "black"),     # 凡例タイトル8pt
    legend.text = element_text(size = 8, color = "black"))
dev.off()

svglite::svglite(paste0(basepath,h_dBDC_path,"/2026.01.24_dimplot_fibro_ISLR.svg"),
                 height = 3,width =3)
FeaturePlot(seurat_obj_cellbender_merge.harmony_fibrosubset, features = "ISLR",
        raster = TRUE, raster.dpi = c(5000,5000), order = T,
        pt.size = 70,  cols = c("grey90", "deeppink")) + NoLegend()+
  theme(
    plot.title = element_text(size = 8, color = "black"),      # タイトル文字サイズ10pt
    axis.title = element_text(size = 8, color = "black"),       # 軸タイトル8pt
    axis.text = element_text(size = 8, color = "black"),        # 軸目盛りテキスト8pt
    legend.title = element_text(size = 8, color = "black"),     # 凡例タイトル8pt
    legend.text = element_text(size = 8, color = "black"))
dev.off()


svglite::svglite(paste0(basepath,h_dBDC_path,"/2026.01.24_dimplot_fibro_THY1.svg"),
                 height = 3,width =3)
FeaturePlot(seurat_obj_cellbender_merge.harmony_fibrosubset, features = "THY1",
            raster = TRUE, raster.dpi = c(5000,5000),  order = T,
            pt.size = 70,  cols = c("grey90", "deeppink")) + NoLegend()+
  theme(
    plot.title = element_text(size = 8, color = "black"),      # タイトル文字サイズ10pt
    axis.title = element_text(size = 8, color = "black"),       # 軸タイトル8pt
    axis.text = element_text(size = 8, color = "black"),        # 軸目盛りテキスト8pt
    legend.title = element_text(size = 8, color = "black"),     # 凡例タイトル8pt
    legend.text = element_text(size = 8, color = "black"))
dev.off()


svglite::svglite(paste0(basepath,h_dBDC_path,"/2026.01.24_dimplot_fibro_ACTA2.svg"),
                 height = 3,width =3)
FeaturePlot(seurat_obj_cellbender_merge.harmony_fibrosubset, features = "ACTA2",
            raster = TRUE, raster.dpi = c(5000,5000),  order = T,
            pt.size = 70,  cols = c("grey90", "deeppink")) + NoLegend()+
  theme(
    plot.title = element_text(size = 8, color = "black"),      # タイトル文字サイズ10pt
    axis.title = element_text(size = 8, color = "black"),       # 軸タイトル8pt
    axis.text = element_text(size = 8, color = "black"),        # 軸目盛りテキスト8pt
    legend.title = element_text(size = 8, color = "black"),     # 凡例タイトル8pt
    legend.text = element_text(size = 8, color = "black"))
dev.off()



svglite::svglite(paste0(basepath,h_dBDC_path,"/2026.01.24_dimplot_fibro_scale.svg"),
                 height = 8,width =8)
FeaturePlot(seurat_obj_cellbender_merge.harmony_fibrosubset, features = c("ISLR","THY1","ACTA2"),
            cols = c("grey90", "deeppink")) +
  theme(legend.title = element_text(size = 8, color = "black"),     # 凡例タイトル8pt
    legend.text = element_text(size = 8, color = "black"))
dev.off()


svglite::svglite(
  paste0(basepath, h_dBDC_path, "/2026.01.24_dimplot_fibro_PI16.svg"),
  height = 3, width = 3
)

FeaturePlot(
  seurat_obj_cellbender_merge.harmony_fibrosubset,
  features = "PI16",
  raster = TRUE,
  raster.dpi = c(5000, 5000),
  order = TRUE,
  pt.size = 70,
  cols = c("grey90", "deeppink")
) +
  NoLegend() +
  theme(
    plot.title = element_text(size = 8, color = "black"),
    axis.title = element_text(size = 8, color = "black"),
    axis.text  = element_text(size = 8, color = "black")
  )

dev.off()

svglite::svglite(
  paste0(basepath, h_dBDC_path, "/2026.01.24_dimplot_fibro_CXCL12.svg"),
  height = 3, width = 3
)

FeaturePlot(
  seurat_obj_cellbender_merge.harmony_fibrosubset,
  features = "CXCL12",
  raster = TRUE,
  raster.dpi = c(5000, 5000),
  order = TRUE,
  pt.size = 70,
  cols = c("grey90", "deeppink")
) +
  NoLegend() +
  theme(
    plot.title = element_text(size = 8, color = "black"),
    axis.title = element_text(size = 8, color = "black"),
    axis.text  = element_text(size = 8, color = "black")
  )

dev.off()

svglite::svglite(
  paste0(basepath, h_dBDC_path, "/2026.01.24_dimplot_fibro_COL1A1.svg"),
  height = 3, width = 3
)

FeaturePlot(
  seurat_obj_cellbender_merge.harmony_fibrosubset,
  features = "COL1A1",
  raster = TRUE,
  raster.dpi = c(5000, 5000),
  order = TRUE,
  pt.size = 70,
  cols = c("grey90", "deeppink")
) +
  NoLegend() +
  theme(
    plot.title = element_text(size = 8, color = "black"),
    axis.title = element_text(size = 8, color = "black"),
    axis.text  = element_text(size = 8, color = "black")
  )

dev.off()

svglite::svglite(
  paste0(basepath, h_dBDC_path, "/2026.01.24_dimplot_fibro_scale.svg"),
  height = 8, width = 8
)

FeaturePlot(
  seurat_obj_cellbender_merge.harmony_fibrosubset,
  features = c("PI16", "CXCL12", "COL1A1"),
  cols = c("grey90", "deeppink")
) +
  theme(
    legend.title = element_text(size = 8, color = "black"),
    legend.text  = element_text(size = 8, color = "black")
  )

dev.off()

subset_fibro<-seurat_obj_cellbender_merge.harmony_fibrosubset%>%subset(ident = c("myCAF","Thy1_ssl","Meflin_ssl"))
subset_fibro$annotations<-factor(subset_fibro$annotations,levels=c("Thy1_ssl","Meflin_ssl","myCAF"))

svglite::svglite(
  paste0(basepath, h_dBDC_path, "/2026.01.24_vlnplot_fibro_6genes.svg"),
  height = 2.5,
  width  = 8)
VlnPlot(
  subset_fibro,
  features = c(
    "ISLR",
    "THY1",
    "ACTA2",    
    "PI16",
    "COL1A1"
  ),raster = T,group.by = "annotations",
  pt.size = 0.01,      # 
  ncol = 5
) +
  NoLegend() &
  theme(
    axis.title = element_blank(),                 # 軸タイトル不要
    axis.text  = element_text(size = 8, color = "black"),
    strip.text = element_text(size = 8, color = "black"),
    plot.title = element_text(size = 8, color = "black")
  )

dev.off()

pdf(
  paste0(basepath, h_dBDC_path, "/2026.01.31_vlnplot_fibro_6genes.pdf"),
  height = 2.5,
  width  = 8)
VlnPlot(
  subset_fibro,
  features = c(
    "ISLR",
    "THY1",
    "ACTA2",    
    "PI16",
    "COL1A1"
  ),raster = T,group.by = "annotations",
  pt.size = 0.01,      # 
  ncol = 5
) +
  NoLegend() &
  theme(
    axis.title = element_blank(),                 # 軸タイトル不要
    axis.text  = element_text(size = 8, color = "black"),
    strip.text = element_text(size = 8, color = "black"),
    plot.title = element_text(size = 8, color = "black")
  )

dev.off()

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




