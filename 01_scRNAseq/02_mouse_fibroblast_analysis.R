# ============================================================
# Mouse bilio-vascular fibroblast scRNA-seq processing and trajectory analysis
#
# Study      : Distinct cancer-associated fibroblast lineages and an
#              Am80-based stromal-reprogramming strategy in extrahepatic
#              cholangiocarcinoma (Meflin/ISLR).
# Repository : eCCA_AM80  (https://github.com/shirakiyukihiro/eCCA_AM80)
#
# Purpose    : Read the GSE163777 10x matrices, build the Seurat object, extract fibroblasts (excluding Lrat+ HSCs), cluster, and run monocle3/Slingshot-style pseudotime and marker analysis.
# Figure(s)  : Mouse scRNA-seq panels (Figure 5A-E and related supplementary figures).
# Input      : data/GSE163777/ (GSM4987107 matrices).
# Output     : Seurat/monocle objects and plots (mouse fibroblast clusters, pseudotime).
#
# Paths have been made relative for public release: place input files
# in ./data and write outputs to ./output (create these folders, or edit
# the paths below).
# ============================================================

path="data/GSE163777"
library(tidyverse)
list.files(path%>%paste0("/GSE163777"))

library(Matrix)
mat<-readMM(file = path%>%paste0("/GSE163777/GSM4987107_matrix.mtx.gz"))
feature.names = read.delim( path%>%paste0("/GSE163777/GSM4987107_features.tsv.gz"), 
                            header = FALSE,
                            stringsAsFactors = FALSE)
barcode.names = read.delim( path%>%paste0("/GSE163777/GSM4987107_barcodes.tsv.gz"), 
                            header = FALSE,
                            stringsAsFactors = FALSE)
colnames(mat) <- barcode.names$V1
rownames(mat) <- feature.names$V1
rownames(mat)[!duplicated(feature.names$V2)]<-feature.names$V2[!duplicated(feature.names$V2)]
sc<-CreateSeuratObject(mat, project = "portal_fib", min.cells = 3, min.features = 200)

sc[["percent.mt"]] <-PercentageFeatureSet(sc, pattern = "^mt-|^MT-")

VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
sc <- SCTransform(sc, vars.to.regress = "percent.mt", verbose = FALSE)

list.files("Z:/R/data/")
#amigo2からUser filters +	isa_partof_closure: GO:0007165 +	taxon_subset_closure_label: Mus musculus	+	evidence_subset_closure_label: experimental evidence	+	type: protein_coding_gene
sigtra<-read.table("Z:/R/data/Amigo_signal_transduction_mm_experimental_protein_coding_gene.txt",sep = "\t",header = F)
#amigo2からUser filters	+	isa_partof_closure: GO:0003700	+	taxon_subset_closure_label: Mus musculus	+	evidence_subset_closure_label: experimental evidence
DbTFa<-read.table("Z:/R/data/Amigo_DNA_bind_TF_activity_mm_experimental.txt",sep = "\t",header = F)

library(clusterProfiler)#BiocManager::install("clusterProfiler")
library(org.Mm.eg.db)#BiocManager::install("org.Mm.eg.db")
sigtra%>%distinct(V1)%>%unlist%>%bitr(fromType = "MGI",toType ="SYMBOL",OrgDb = org.Mm.eg.db)->st_gene
DbTFa%>%distinct(V1)%>%unlist%>%bitr(fromType = "MGI",toType ="SYMBOL",OrgDb = org.Mm.eg.db)->TF_gene
c(st_gene$SYMBOL,TF_gene$SYMBOL)%>%unique->pcagene

sc <- RunPCA(sc, features =  rownames(sc)%>%intersect(pcagene) )
sc <- FindNeighbors(sc, dims = seq(PCAtools::findElbowPoint(Stdev(sc,"pca")^2)))
sc <- FindClusters(sc, resolution = 0.5) #resolutionを高くするとクラスターの数が増える
head(Idents(sc), 5)
sc <- RunUMAP(sc, dims = seq(PCAtools::findElbowPoint(Stdev(sc,"pca")^2)))
DimPlot(sc, reduction = "umap")
VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)#8,9はpercent.mt高いので除く
sc_<-subset(sc,ident=c(0,1,2,3,4,5,6,7,10,11,12))

DimPlot(sc_, label = TRUE)
FeaturePlot(sc_, features = c("Pecam1", "Col1a1", "Krt18","Cd14"), min.cutoff = "q9", label = TRUE)
FeaturePlot(sc_, features = c("Islr", "Thy1", "Acta2","Pi16"), min.cutoff = "q9", label = TRUE)
FeaturePlot(sc_, features = c("Mcam", "Msln", "Mxra8","Grem1"), min.cutoff = "q9", label = TRUE)
FeaturePlot(sc_, features = c("Lrat", "Gfap", "Ngfr","Synm"), min.cutoff = "q9", label = TRUE)

sc_fib<-subset(sc,ident=c(0,1,2,3,4,5,12))

sc_fib <- RunPCA(sc_fib, features =  rownames(sc_fib)%>%intersect(pcagene))#(ECM_gene$SYMBOL) )
sc_fib <- FindNeighbors(sc_fib, dims = seq(PCAtools::findElbowPoint(Stdev(sc_fib,"pca")^2)))
sc_fib <- FindClusters(sc_fib, resolution = 0.5) #resolutionを高くするとクラスターの数が増える
head(Idents(sc_fib), 5)
sc_fib <- RunUMAP(sc_fib, dims = seq(PCAtools::findElbowPoint(Stdev(sc_fib,"pca")^2)))
DimPlot(sc_fib, reduction = "umap")
VlnPlot(sc_fib, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)#

FeaturePlot(sc_fib, features = c("Islr", "Col1a1", "Thy1","Acta2"), min.cutoff = "q1", label = TRUE)
FeaturePlot(sc_fib, features = c("Pi16", "Mcam", "Grem1","Lrat"), min.cutoff = "q1", label = TRUE)
FeaturePlot(sc_fib, features = c("Ngfr", "Synm", "Rgs5","Cspg4"), min.cutoff = "q1", label = TRUE)


#extracellular matrix Go:0031012
#amigo2からUser filters	+	isa_partof_closure: GO:0031012	+	taxon_subset_closure_label: Mus musculus	+	evidence_subset_closure_label: experimental evidence
ECM<-read.table("Z:/R/data/Amigo_ECM_mm_experimental.txt",sep = "\t",header = F)
ECM%>%distinct(V1)%>%unlist%>%bitr(fromType = "MGI",toType ="SYMBOL",OrgDb = org.Mm.eg.db)->ECM_gene
sc_fib_ECM <- RunPCA(sc_fib, features =  rownames(sc_fib)%>%intersect(ECM_gene$SYMBOL) )
sc_fib_ECM <- FindNeighbors(sc_fib_ECM, dims = seq(PCAtools::findElbowPoint(Stdev(sc_fib_ECM,"pca")^2)))
sc_fib_ECM <- FindClusters(sc_fib_ECM, resolution = 0.5) #resolutionを高くするとクラスターの数が増える
sc_fib_ECM <- RunUMAP(sc_fib_ECM, dims = seq(PCAtools::findElbowPoint(Stdev(sc_fib_ECM,"pca")^2)))
DimPlot(sc_fib_ECM, reduction = "umap")
FeaturePlot(sc_fib_ECM, features = c("Islr", "Col1a1", "Thy1","Acta2"), min.cutoff = "q50", label = TRUE)
FeaturePlot(sc_fib_ECM, features = c("Pi16", "Mcam", "Ngfr","Lrat"), min.cutoff = "q50", label = TRUE)
FeaturePlot(sc_fib_ECM, features = c("Il6", "Cxcl12", "Pdgfra","Cxcl2"), min.cutoff = "q50", label = TRUE)
FeaturePlot(sc_fib_ECM, features = c("Cd74", "H2-Aa", "H2-Ab1","H2-Eb1"), min.cutoff = "q50", label = TRUE)
VlnPlot(sc_fib_ECM, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)#
sc_fib_ECM.markers <- FindAllMarkers(sc_fib_ECM, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
top30_all.markers <- sc_fib_ECM.markers %>% group_by(cluster) %>% top_n(n =30, wt=avg_log2FC)



#devtools::install_github('cole-trapnell-lab/monocle3')
#remotes::install_github('satijalab/seurat-wrappers')
library(monocle3)
library(SeuratWrappers)
cds_ <- as.cell_data_set(sc_fib_ECM)
fData(cds_)$gene_short_name <- rownames(fData(cds_))
head(fData(cds_))

recreate.partitions <- c(rep(1, length(cds_@colData@rownames)))
names(recreate.partitions) <- cds_@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
recreate.partitions

cds_@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions

list.cluster <- sc_fib_ECM@active.ident
cds_@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster
plot_cells(cds_, color_cells_by = "cluster", label_groups_by_cluster = F, 
           group_label_size = 5) + theme(legend.position = "right")
cds_ <- learn_graph(cds_, use_partition = F)
plot_cells(cds_, color_cells_by = "cluster", label_groups_by_cluster = F,
           label_branch_points = T, label_roots = T, label_leaves = F,
           group_label_size = 5,label_principal_points = T)#show node name
plot_cells(cds_,genes="Islr",label_cell_groups=T,show_trajectory_graph=T,cell_size = 0.8,label_principal_points = T)
plot_cells(cds_,genes="Thy1",label_cell_groups=T,show_trajectory_graph=T,cell_size = 0.8,label_principal_points = T)
plot_cells(cds_,genes="Lrat",label_cell_groups=T,show_trajectory_graph=T,cell_size = 0.8,label_principal_points = T)
plot_cells(cds_,genes="Acta2",label_cell_groups=T,show_trajectory_graph=T,cell_size = 0.8,label_principal_points = T)

cds_sub_Islr <- choose_graph_segments(cds_,starting_pr_node = "Y_13",ending_pr_nodes = "Y_257")#Islr_highからActa2_highで選ぶ
cds_sub_Thy1 <- choose_graph_segments(cds_,starting_pr_node = "Y_154",ending_pr_nodes = "Y_257")
cds_sub_Lrat <- choose_graph_segments(cds_,starting_pr_node = "Y_199",ending_pr_nodes = "Y_257")

cds_sub_Islr_order <- order_cells(cds_[,colnames(cds_sub_Islr)], reduction_method = "UMAP",root_pr_nodes = "Y_13")#Islr_highをrootに選ぶ
plot_cells(cds_sub_Islr_order, color_cells_by = "pseudotime", label_groups_by_cluster = T,
           label_branch_points = T, label_roots = F, label_leaves = F)
plot_cells(cds_sub_Islr_order,genes="Islr",label_cell_groups=T,show_trajectory_graph=T,cell_size = 0.8)
plot_cells(cds_sub_Islr_order,genes="Il6",label_cell_groups=T,show_trajectory_graph=T,cell_size = 0.8)
plot_genes_in_pseudotime(cds_sub_Islr_order["Islr",],cell_size = 1)
pseudotime(cds_sub_Islr_order)
sce_ISLR <- as.SingleCellExperiment(sc_fib_ECM[,colnames(cds_sub_Islr)],assay = "RNA")
library(scran)#BiocManager::install("scran")
clust <- quickCluster(sce_ISLR) 
sce_ISLR <- computeSumFactors(sce_ISLR, cluster=clust, min.mean=0.1)
sce_ISLR <- logNormCounts(sce_ISLR)
sce_ISLR$pseudotime<-pseudotime(cds_sub_Islr_order)
library(scater)#BiocManager::install("scater")
plotExpression(sce_ISLR, "Islr", x = "pseudotime", colour_by = "ident", show_violin = FALSE,
               show_smooth = TRUE,exprs_values = "logcounts")
plotExpression(sce_ISLR, "Il6", x = "pseudotime", colour_by = "ident", show_violin = FALSE,
               show_smooth = TRUE,exprs_values = "logcounts")
plotExpression(sce_ISLR, "Cxcl12", x = "pseudotime", colour_by = "ident", show_violin = FALSE,
               show_smooth = TRUE,exprs_values = "logcounts")

cds_sub_Thy1_order <- order_cells(cds_[,colnames(cds_sub_Thy1)], reduction_method = "UMAP",root_pr_nodes = "Y_154")
plot_cells(cds_sub_Thy1_order, color_cells_by = "pseudotime", label_groups_by_cluster = T,
           label_branch_points = T, label_roots = F, label_leaves = F)
plot_cells(cds_sub_Thy1_order,genes="Islr",label_cell_groups=T,show_trajectory_graph=T,cell_size = 0.8)
plot_cells(cds_sub_Thy1_order,genes="Thy1",label_cell_groups=T,show_trajectory_graph=T,cell_size = 0.8)
pseudotime(cds_sub_Thy1_order)
sce_Thy1 <- as.SingleCellExperiment(sc_fib_ECM[,colnames(cds_sub_Thy1)],assay = "RNA")
clust <- quickCluster(sce_Thy1) 
sce_Thy1 <- computeSumFactors(sce_Thy1, cluster=clust, min.mean=0.1)
sce_Thy1 <- logNormCounts(sce_Thy1)
sce_Thy1$pseudotime<-pseudotime(cds_sub_Thy1_order)
plotExpression(sce_Thy1, "Thy1", x = "pseudotime", colour_by = "ident", show_violin = FALSE,
               show_smooth = TRUE,exprs_values = "logcounts")

cds_sub_Lrat_order <- order_cells(cds_[,colnames(cds_sub_Lrat)], reduction_method = "UMAP",root_pr_nodes = "Y_199")
plot_cells(cds_sub_Lrat_order, color_cells_by = "pseudotime", label_groups_by_cluster = T,
           label_branch_points = T, label_roots = F, label_leaves = F)
plot_cells(cds_sub_Lrat_order,genes="Lrat",label_cell_groups=T,show_trajectory_graph=T,cell_size = 0.8)
pseudotime(cds_sub_Lrat_order)
sce_Lrat <- as.SingleCellExperiment(sc_fib_ECM[,colnames(cds_sub_Lrat)],assay = "RNA")
clust <- quickCluster(sce_Lrat) 
sce_Lrat <- computeSumFactors(sce_Lrat, cluster=clust, min.mean=0.1)
sce_Lrat <- logNormCounts(sce_Lrat)
sce_Lrat$pseudotime<-pseudotime(cds_sub_Lrat_order)
plotExpression(sce_Lrat, "Lrat", x = "pseudotime", colour_by = "ident", show_violin = FALSE,
               show_smooth = TRUE,exprs_values = "logcounts")
plotExpression(sce_Lrat, "Ngfr", x = "pseudotime", colour_by = "ident", show_violin = FALSE,
               show_smooth = TRUE,exprs_values = "logcounts")


# Only look at the 1,000 most variable genes when identifying temporally expressesd genes.
# Identify the variable genes by ranking all genes by their variance.
#Y <- log2(counts(sce_ISLR) + 1)
#var1K <- names(sort(apply(Y, 1, var),decreasing = TRUE))[1:1000]
#Y <- Y[var1K, ]  # only counts for variable genes
#デフォルトは変動遺伝子だが、signal transductionのGOの遺伝子にしてみる。
#Y <- Y[rownames(Y)%>%intersect(st_gene$SYMBOL), ]
Y <- logcounts(sce_ISLR)[rownames(sce_ISLR)%>%intersect(st_gene$SYMBOL), ]

library(gam)
# Fit GAM for each gene using pseudotime as independent variable.
t <- sce_ISLR$pseudotime
gam.pval <- apply(Y, 1, function(z){
  d <- data.frame(z=z, t=t)
  tmp <- gam(z ~ lo(t), data=d)
  p <- summary(tmp)[4][[1]][1,5]
  p
})

# Identify genes with the most significant time-dependent model fit.
topgenes <- names(sort(gam.pval, decreasing = FALSE))[1:30]  

# Prepare and plot a heatmap of the top genes that vary their expression over pseudotime.
library(clusterExperiment)#BiocManager::install("clusterExperiment")
library(tidyverse)
heatdata <- sce_ISLR[topgenes, order(sce_ISLR$pseudotime, na.last = NA)]%>%assay("logcounts")
heatclus <- sce_ISLR$ident[order(sce_ISLR$pseudotime, na.last = NA)]
ce <- ClusterExperiment(heatdata%>%as.matrix, heatclus)
clusterExperiment::plotHeatmap(ce, clusterSamplesData = "orderSamplesValue", visualizeData = 'transformed')

t <- sce_ISLR$pseudotime

gam <- apply(Y, 1, function(z){
  d <- data.frame(z=z, t=t)
  tmp <- gam(z ~ lo(t), data=d)
  return(tmp)
})


gam$Ptp4a3[2][[1]][[1]]
(gam$C3%>%summary)
gam[["C3"]][2][[1]][[1]]

coef<-tibble(coef=sapply(names(gam), function(x){gam[[x]][2][[1]][[1]]}),gene=names(gam),
             pval=sapply(names(gam), function(x){summary(gam[[x]])[4][[1]][1,5]}))
coef%>%filter(pval<0.05)%>%arrange(coef)->geneList


#GSEA
library(clusterProfiler)
library(msigdbr)#install.packages("msigdbr")
msigdbr_species()
m_df <- msigdbr(species = "Mus musculus")
head(m_df, 2) %>% as.data.frame
m_df%>%dplyr::select(gs_name, gene_symbol)
em <- enricher(geneList$gene%>%unique(),
               minGSSize    = 20, 
               TERM2GENE=m_df%>%dplyr::select(gs_name, gene_symbol))
head(em)
dotplot(em, showCategory=50) + ggtitle("dotplot for MSigDB")

m_df%>%filter(gs_cat=="H")%>%dplyr::select(gs_name, gene_symbol)

geneList$coef->genelist
GSEA <- GSEA(genelist%>%sort(decreasing = T),
             TERM2GENE = m_df%>%dplyr::select(gs_name, gene_symbol))
head(GSEA)


dotplot(GSEA, showCategory=200,x="Count",font.size=5) + ggtitle("dotplot for GSEA")
dotplot(GSEA, showCategory=200,font.size=5) + ggtitle("dotplot for GSEA")

ridgeplot(GSEA, showCategory=200)

GSEA_H <- GSEA(genelist%>%sort(decreasing = T),
             TERM2GENE = m_df%>%filter(gs_cat=="H")%>%dplyr::select(gs_name, gene_symbol))
dotplot(GSEA_H, showCategory=200,font.size=5) + ggtitle("dotplot for GSEA")
ridgeplot(GSEA_H, showCategory=200)

library(enrichplot)
gseaplot2(GSEA_H, geneSetID = 1, title = GSEA_H$Description[1])
gseaplot2(GSEA_H, geneSetID = 2, title = GSEA_H$Description[2])
gseaplot2(GSEA_H, geneSetID = 3, title = GSEA_H$Description[3])
gseaplot2(GSEA_H, geneSetID = 4, title = GSEA_H$Description[4])
gseaplot2(GSEA_H, geneSetID = 5, title = GSEA_H$Description[5])
gseaplot2(GSEA_H, geneSetID = 6, title = GSEA_H$Description[6])
gseaplot2(GSEA, geneSetID = 1:10)
gseaplot2(GSEA, geneSetID = 11:20)
gseaplot2(GSEA, geneSetID = 21:30)
gseaplot2(GSEA, geneSetID = 31:40)
gseaplot2(GSEA, geneSetID = 41:50)
gseaplot2(GSEA, geneSetID = 51:60)
gseaplot2(GSEA, geneSetID = 61:70)
gseaplot2(GSEA, geneSetID = 71:80)
gseaplot2(GSEA, geneSetID = 81:90)
gseaplot2(GSEA, geneSetID = 91:100)
gseaplot2(GSEA, geneSetID = 101:110)


Y <- logcounts(sce_Thy1)[rownames(sce_Thy1)%>%intersect(st_gene$SYMBOL), ]
t <- sce_Thy1$pseudotime

gam <- apply(Y, 1, function(z){
  d <- data.frame(z=z, t=t)
  tmp <- gam(z ~ lo(t), data=d)
  return(tmp)
})

coef_thy1<-tibble(coef=sapply(names(gam), function(x){gam[[x]][2][[1]][[1]]}),gene=names(gam),
             pval=sapply(names(gam), function(x){summary(gam[[x]])[4][[1]][1,5]}))
coef_thy1%>%filter(pval<0.05)%>%arrange(coef)->geneList_thy1

geneList_thy1$coef->genelist
GSEA_thy1 <- GSEA(genelist%>%sort(decreasing = T),
                  TERM2GENE = m_df%>%filter(gs_cat=="H")%>%dplyr::select(gs_name, gene_symbol))

ridgeplot(GSEA_thy1, showCategory=200)

library(org.Mm.eg.db)
geneList$coef->genelist_Islr
names(genelist_Islr)%>%bitr(fromType = "SYMBOL",toType ="ENTREZID",OrgDb = org.Mm.eg.db)->namesls
names(genelist_Islr)<-namesls$ENTREZID
kk_Islr <- gseKEGG(geneList     = genelist_Islr[!duplicated(names(genelist_Islr))]%>%sort(decreasing = T),
               organism     = 'mmu', #mmu
               pvalueCutoff = 0.05,
               verbose      = FALSE)
dotplot(kk_Islr, showCategory=30) + ggtitle("dotplot for KEGG")

geneList_thy1$coef->genelist_thy1
names(genelist_thy1)%>%bitr(fromType = "SYMBOL",toType ="ENTREZID",OrgDb = org.Mm.eg.db)->namesls
names(genelist_thy1)<-namesls$ENTREZID

kk_Thy1 <- gseKEGG(geneList     = genelist_thy1[!duplicated(names(genelist_thy1))]%>%sort(decreasing = T),
                   organism     = 'mmu', #mmu
                   minGSSize    = 10,
                   pvalueCutoff = 0.05,
                   verbose      = FALSE)
dotplot(kk_Thy1, showCategory=30) + ggtitle("dotplot for KEGG")




path
pdfpath=
pdf(pdfpath)
dotplot(ego_result, showCategory=200,font.size=5) + ggtitle("dotplot for GO")
cnetplot(ego_result,cex_label_gene = 0.5,cex_label_category=0.5, showCategory=10) 
dotplot(order, showCategory=20,font.size=8) + ggtitle("dotplot for GO")
dotplot(kk, showCategory=30,font.size=5) + ggtitle("dotplot for KEGG")
cnetplot(kk,cex_label_gene = 0.5,cex_label_category=0.5, showCategory=10)  
dotplot(WP, showCategory=30,font.size=5) + ggtitle("dotplot for WikiPathways")
dotplot(RP, showCategory=100,font.size=5) + ggtitle("dotplot for REACTOME pathway")

dotplot(GSEA, showCategory=200,font.size=5) + ggtitle("dotplot for MSigDB:GSEA")

dev.off()




















deg_Islr <- graph_test(cds_sub_Islr_order, neighbor_graph = "principal_graph")
deg_Islr %>% arrange(q_value) %>% filter(status == "OK") %>% head()
deg_Thy1 <- graph_test(cds_sub_Thy1_order, neighbor_graph = "principal_graph")
deg_Thy1 %>% arrange(q_value) %>% filter(status == "OK") %>% head()
deg_Lrat <- graph_test(cds_sub_Lrat_order, neighbor_graph = "principal_graph")
deg_Lrat %>% arrange(q_value) %>% filter(status == "OK") %>% head()


deg_Islr %>% arrange(q_value) %>% filter(status == "OK",q_value==0) ->deg_Islr_q0
deg_Thy1 %>% arrange(q_value) %>% filter(status == "OK",q_value==0) ->deg_Thy1_q0
deg_Lrat %>% arrange(q_value) %>% filter(status == "OK",q_value==0) ->deg_Lrat_q0

deg_Islr_q0$gene_short_name%>%setdiff(deg_Thy1_q0$gene_short_name) 
deg_Islr_q0$gene_short_name%>%setdiff(deg_Lrat_q0$gene_short_name) 
deg_Thy1_q0$gene_short_name%>%setdiff(deg_Islr_q0$gene_short_name) 
deg_Thy1_q0$gene_short_name%>%setdiff(deg_Lrat_q0$gene_short_name) 
deg_Lrat_q0$gene_short_name%>%setdiff(deg_Islr_q0$gene_short_name) 
deg_Lrat_q0$gene_short_name%>%setdiff(deg_Thy1_q0$gene_short_name) 





















# Prepare and plot a heatmap of the top genes that vary their expression over pseudotime.
library(clusterExperiment)#BiocManager::install("clusterExperiment")
heatdata <- sce_ISLR[rownames(deg_Islr_q0)%>%setdiff(rownames(deg_Thy1_q0)), order(sce_ISLR$pseudotime, na.last = NA)]%>%assay("logcounts")
heatclus <- sce_ISLR$ident[order(sce_ISLR$pseudotime, na.last = NA)]
ce <- ClusterExperiment(heatdata%>%as.matrix, heatclus, transformation = log1p)
clusterExperiment::plotHeatmap(ce, clusterSamplesData = "orderSamplesValue", visualizeData = 'transformed')

heatdata <- sce_ISLR[rownames(deg_Islr_q0)%>%setdiff(rownames(deg_Lrat_q0)), order(sce_ISLR$pseudotime, na.last = NA)]%>%assay("logcounts")
heatclus <- sce_ISLR$ident[order(sce_ISLR$pseudotime, na.last = NA)]
ce <- ClusterExperiment(heatdata%>%as.matrix, heatclus, transformation = log1p)
clusterExperiment::plotHeatmap(ce, clusterSamplesData = "orderSamplesValue", visualizeData = 'transformed')


table(deg_Islr_q0$morans_I<0)
pdf(file ="./R/data/pdf/2023.10.19_portal_fib_monocle.pdf")
FeaturePlot(sc_fib_ECM, features = "Islr", min.cutoff = "q85", label = TRUE,order = T,cols = c("grey", "magenta"))->plot;print(plot)
FeaturePlot(sc_fib_ECM, features = "Acta2", min.cutoff = "q85", label = TRUE,order = T,cols = c("grey", "magenta"))->plot;print(plot)
FeaturePlot(sc_fib_ECM, features = "Thy1", min.cutoff = "q85", label = TRUE,order = T,cols = c("grey", "magenta"))->plot;print(plot)
FeaturePlot(sc_fib_ECM, features = "Lrat", min.cutoff = "q85", label = TRUE,order = T,cols = c("grey", "magenta"))->plot;print(plot)
FeaturePlot(sc_fib_ECM, features = "Ngfr", min.cutoff = "q85", label = TRUE,order = T,cols = c("grey", "magenta"))->plot;print(plot)
plot_cells(cds_, color_cells_by = "cluster", label_groups_by_cluster = F,
           label_branch_points = T, label_roots = T, label_leaves = F,
           group_label_size = 5,label_principal_points = T)->plot;print(plot)
plot_cells(cds_sub_Islr_order, color_cells_by = "pseudotime", label_groups_by_cluster = T,
           label_branch_points = T, label_roots = F, label_leaves = F)->plot;print(plot)
plotExpression(sce_ISLR, "Islr", x = "pseudotime", colour_by = "ident", show_violin = FALSE,
               show_smooth = TRUE,exprs_values = "logcounts")->plot;print(plot)
plot_cells(cds_sub_Thy1_order, color_cells_by = "pseudotime", label_groups_by_cluster = T,
           label_branch_points = T, label_roots = F, label_leaves = F)->plot;print(plot)
plotExpression(sce_Thy1, "Thy1", x = "pseudotime", colour_by = "ident", show_violin = FALSE,
               show_smooth = TRUE,exprs_values = "logcounts")->plot;print(plot)
plot_cells(cds_sub_Lrat_order, color_cells_by = "pseudotime", label_groups_by_cluster = T,
           label_branch_points = T, label_roots = F, label_leaves = F)->plot;print(plot)
plotExpression(sce_Lrat, "Lrat", x = "pseudotime", colour_by = "ident", show_violin = FALSE,
               show_smooth = TRUE,exprs_values = "logcounts")->plot;print(plot)
plotExpression(sce_Lrat, "Ngfr", x = "pseudotime", colour_by = "ident", show_violin = FALSE,
               show_smooth = TRUE,exprs_values = "logcounts")->plot;print(plot)
dev.off()



pdf(file ="./R/data/pdf/2023.10.19_portal_fib_monocle_deg_heatmap.pdf")
heatdata <- sce_ISLR[rownames(deg_Islr_q0)%>%setdiff(rownames(deg_Thy1_q0)), order(sce_ISLR$pseudotime, na.last = NA)]%>%assay("logcounts")
heatclus <- sce_ISLR$ident[order(sce_ISLR$pseudotime, na.last = NA)]
ce <- ClusterExperiment(heatdata%>%as.matrix, heatclus, transformation = log1p)
clusterExperiment::plotHeatmap(ce, clusterSamplesData = "orderSamplesValue", visualizeData = 'transformed')
heatdata <- sce_ISLR[rownames(deg_Islr_q0)%>%setdiff(rownames(deg_Lrat_q0)), order(sce_ISLR$pseudotime, na.last = NA)]%>%assay("logcounts")
heatclus <- sce_ISLR$ident[order(sce_ISLR$pseudotime, na.last = NA)]
ce <- ClusterExperiment(heatdata%>%as.matrix, heatclus, transformation = log1p)
clusterExperiment::plotHeatmap(ce, clusterSamplesData = "orderSamplesValue", visualizeData = 'transformed')
heatdata <- sce_Thy1[rownames(deg_Thy1_q0)%>%setdiff(rownames(deg_Islr_q0)), order(sce_Thy1$pseudotime, na.last = NA)]%>%assay("logcounts")
heatclus <- sce_Thy1$ident[order(sce_Thy1$pseudotime, na.last = NA)]
ce <- ClusterExperiment(heatdata%>%as.matrix, heatclus, transformation = log1p)
clusterExperiment::plotHeatmap(ce, clusterSamplesData = "orderSamplesValue", visualizeData = 'transformed')
heatdata <- sce_Thy1[rownames(deg_Thy1_q0)%>%setdiff(rownames(deg_Lrat_q0)), order(sce_Thy1$pseudotime, na.last = NA)]%>%assay("logcounts")
heatclus <- sce_Thy1$ident[order(sce_Thy1$pseudotime, na.last = NA)]
ce <- ClusterExperiment(heatdata%>%as.matrix, heatclus, transformation = log1p)
clusterExperiment::plotHeatmap(ce, clusterSamplesData = "orderSamplesValue", visualizeData = 'transformed')
heatdata <- sce_Lrat[rownames(deg_Lrat_q0)%>%setdiff(rownames(deg_Islr_q0)), order(sce_Lrat$pseudotime, na.last = NA)]%>%assay("logcounts")
heatclus <- sce_Lrat$ident[order(sce_Lrat$pseudotime, na.last = NA)]
ce <- ClusterExperiment(heatdata%>%as.matrix, heatclus, transformation = log1p)
clusterExperiment::plotHeatmap(ce, clusterSamplesData = "orderSamplesValue", visualizeData = 'transformed')
heatdata <- sce_Lrat[rownames(deg_Lrat_q0)%>%setdiff(rownames(deg_Thy1_q0)), order(sce_Lrat$pseudotime, na.last = NA)]%>%assay("logcounts")
heatclus <- sce_Lrat$ident[order(sce_Lrat$pseudotime, na.last = NA)]
ce <- ClusterExperiment(heatdata%>%as.matrix, heatclus, transformation = log1p)
clusterExperiment::plotHeatmap(ce, clusterSamplesData = "orderSamplesValue", visualizeData = 'transformed')
dev.off()


library(msigdbr)#install.packages("msigdbr")
msigdbr_show_species()
m_df <- msigdbr(species = "Mus musculus")
head(m_df, 2) %>% as.data.frame
m_df%>%dplyr::select(gs_name, gene_symbol)
em <- enricher(rownames(deg_Islr_q0)%>%setdiff(rownames(deg_Thy1_q0)), TERM2GENE=m_df%>%dplyr::select(gs_name, gene_symbol))
head(em)

library(clusterProfiler)#BiocManager::install("clusterProfiler")
library(org.Mm.eg.db)#BiocManager::install("org.Mm.eg.db")
rownames(deg_Islr_q0)%>%setdiff(rownames(deg_Thy1_q0))%>%bitr(fromType = "SYMBOL",toType ="ENTREZID",OrgDb = org.Mm.eg.db)->Islr_thy1
kk <- enrichKEGG(gene         = Islr_thy1$ENTREZID,
                 organism     = 'mmu', #hsa
                 pvalueCutoff = 0.05)
head(kk)
kk@result$Description#[9] "PI3K-Akt signaling pathway" 
kk@result$ID#[9] "PI3K-Akt signaling pathway" 
dotplot(kk, showCategory=30) + ggtitle("dotplot for KEGG")


#BiocManager::install("tradeSeq")
library(tradeSeq)
library(magrittr)
#install.packages("raster", dependencies = F)
#install.packages("gsl", dependencies = F)#sudo apt-get install libgsl-dev
#install.packages("energy", dependencies = F)
#install.packages("mgcv", dependencies = F)
library(mgcv)
# Get the closest vertice for every cell
y_to_cells <-  principal_graph_aux(cds_sub_Islr_order)$UMAP$pr_graph_cell_proj_closest_vertex %>%as.data.frame
y_to_cells$cells <- rownames(y_to_cells)
y_to_cells$Y <- y_to_cells$V1
# Get the root vertices
# It is the same node as above
root <- cds_sub_Islr_order@principal_graph_aux$UMAP$root_pr_nodes
mst <- principal_graph(cds_sub_Islr_order)$UMAP
# Get the other endpoints
endpoints <- names(which(igraph::degree(mst) == 1))
endpoints <- endpoints[!endpoints %in% root]

# For each endpoint
cellWeights <- lapply(endpoints, function(endpoint) {
  # We find the path between the endpoint and the root
  path <- igraph::shortest_paths(mst, root, endpoint)$vpath[[1]]
  path <- as.character(path)
  # We find the cells that map along that path
  df <- y_to_cells[y_to_cells$Y %in% path, ]
  df <- data.frame(weights = as.numeric(colnames(cds_sub_Islr_order) %in% df$cells))
  colnames(df) <- endpoint
  return(df)
}) %>% do.call(what = 'cbind', args = .) %>%
  as.matrix()
rownames(cellWeights) <- colnames(cds_sub_Islr_order)
pseudotime <- matrix(pseudotime(cds_sub_Islr_order), ncol = ncol(cellWeights),
                     nrow = ncol(cds_sub_Islr_order), byrow = FALSE)

BPPARAM <- BiocParallel::bpparam()
BPPARAM # lists current options
BPPARAM$workers <- 16 # use  cores
icMat <- evaluateK(counts = cds_sub_Islr_order%>%exprs, pseudotime = pseudotime, cellWeights = cellWeights,
                   k=3:7, nGenes = 100, plot = TRUE)
fitGAM <- fitGAM(counts = cds_sub_Islr_order%>%exprs,pseudotime = pseudotime,cellWeights = cellWeights, nknots = 6,  parallel=TRUE, BPPARAM = BPPARAM)
#fitGAM<-sce
assoRes <- associationTest(fitGAM)
head(assoRes)
startRes <- startVsEndTest(fitGAM)
oStart <- order(startRes$waldStat, decreasing = TRUE)
sigGeneStart <- names(fitGAM)[oStart[3]]
plotSmoothers(fitGAM, counts, gene = sigGeneStart)
plotGeneCount(crv, counts, gene = sigGeneStart)


gene_symbol <- as.list(org.Mm.egSYMBOL)
raw_count_data <- GetAssayData(sc_fib, assay = "RNA", slot = "counts")
class(raw_count_data)

cells_info <- sc_fib@meta.data

gene_name <- gene_symbol[rownames(raw_count_data)]  
gene_name <- sapply(gene_name, function(x) x[[1]][1])

#preparing cds
gene_name <- ifelse(is.na(gene_name), names(gene_name), gene_name)
gene_short_name <- gene_name
gene_id <- rownames(raw_count_data)


genes_info <- cbind(gene_id, gene_short_name)
genes_info <- as.data.frame(genes_info)
rownames(genes_info) <- rownames(raw_count_data)

cds <- new_cell_data_set(expression_data = raw_count_data,
                         cell_metadata = cells_info,
                         gene_metadata = genes_info)

cds <- preprocess_cds(cds, num_dim = 100)
plot_pc_variance_explained(cds)

cds <- reduce_dimension(cds)

plot_cells(cds, color_cells_by="seurat_clusters")
plot_cells(cds, genes="Islr")

cds = cluster_cells(cds, resolution=1e-5)
plot_cells(cds)


#replace monocle clusters with seurat
cds@clusters$UMAP$partitions <- sc_fib@meta.data$seurat_clusters
names(cds@clusters$UMAP$partitions) <- rownames(sc_fib@meta.data)
cds@clusters$UMAP$clusters <- sc_fib@meta.data$seurat_clusters
names(cds@clusters$UMAP$clusters) <- rownames(sc_fib@meta.data)

#replace monocle dimensions with seurat
cds@reduce_dim_aux$UMAP <-  sc_fib@reductions$umap@cell.embeddings
#cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- data@reductions$umap@cell.embeddings

plot_cells(cds, color_cells_by="partition", group_cells_by="partition")

marker_test_res <- top_markers(cds, group_cells_by="partition", 
                               reference_cells=1000, cores=8)

top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(1, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))

plot_genes_by_group(cds,
                    top_specific_marker_ids,
                    group_cells_by="partition",
                    ordering_type="maximal_on_diag",
                    max.size=3)











sc_$cluster<-Idents(sc_)
sc_ <- RenameIdents(sc_, `0` = "Fibroblast", `1` = "Fibroblast", `2` = "Fibroblast", `3` = "Fibroblast", 
                    `4` = "Fibroblast", `5` = "Fibroblast",`6` = "Endothelial_cell", `7` = "Endothelial_cell",
                    `10` = "Mesothelial_cells", `11` = "",`11` = "Fibroblast")
DimPlot(sc_, label = TRUE)
Idents(sc_) <- factor(Idents(sc_), levels = c( "Fibroblast","Macrophage","Dendritic_cell","Neutrophil", 
                                               "T_cell","NK_cell","Epithelial_cell","Endothelial_cell"))
markers.to.plot <-  c("Cd3d", "Col1a1", "Cd79a", "Cd74", "Cd8a", "Krt18", "Islr", "Acta2","Pdgfra","Stra6",
                      "Cd14","Cd34","Pecam1","Gzma","Gzmb","Krt8","Cd4","Cd83")
DotPlot(sc_, features = markers.to.plot[1:10], cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()
DotPlot(sc_, features = markers.to.plot[11:18], cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()
sc_$celltype<-Idents(sc_)

