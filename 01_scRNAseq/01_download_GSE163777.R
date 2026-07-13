# ============================================================
# Download mouse bilio-vascular scRNA-seq (GSE163777) from GEO
#
# Study      : Distinct cancer-associated fibroblast lineages and an
#              Am80-based stromal-reprogramming strategy in extrahepatic
#              cholangiocarcinoma (Meflin/ISLR).
# Repository : eCCA_AM80  (https://github.com/shirakiyukihiro/eCCA_AM80)
#
# Purpose    : Download the public mouse portal/bilio-vascular fibroblast scRNA-seq supplementary files from GEO.
# Figure(s)  : Source dataset for the mouse scRNA-seq analysis (Fig. 3 and related supplementary figures).
# Input      : GEO accession GSE163777 (downloaded into data/GSE163777/).
# Output     : data/GSE163777/ (raw 10x matrices).
#
# Paths have been made relative for public release: place input files
# in ./data and write outputs to ./output (create these folders, or edit
# the paths below).
# ============================================================

BiocManager::install("GEOquery")
library(GEOquery)
dir.create("data/GSE163777")

path="data/GSE163777"

setwd(path)
filePaths = getGEOSuppFiles("GSE163777")
filePaths


library(tidyverse)
list.files(path%>%paste0("/GSE163777"))

library(Matrix)
  mat<-readMM(file = "GSE163777/GSM4987107_matrix.mtx.gz")
  feature.names = read.delim("GSE163777/GSM4987107_features.tsv.gz", 
                             header = FALSE,
                             stringsAsFactors = FALSE)
  barcode.names = read.delim("GSE163777/GSM4987107_barcodes.tsv.gz", 
                             header = FALSE,
                             stringsAsFactors = FALSE)
  colnames(mat) <- barcode.names$V1
  rownames(mat) <- feature.names$V1
  rownames(mat)[!duplicated(feature.names$V2)]<-feature.names$V2[!duplicated(feature.names$V2)]
  sc<-CreateSeuratObject(mat, project = "portal_fib", min.cells = 3, min.features = 200)

percent.mt <-PercentageFeatureSet(sc, pattern = "^mt-|^MT-")
percent.mt

sc[["percent.mt"]]<- percent.mt

VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

seurat<-subset(sc,subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

seurat<-NormalizeData(seurat)
seurat<-FindVariableFeatures(seurat,nfeatures = 2000)
seurat<-ScaleData(seurat,features=rownames(seurat))
seurat<-RunPCA(seurat,features=VariableFeatures(seurat))
library(PCAtools)#BiocManager::install("PCAtools")
seurat<-FindNeighbors(seurat, dims = seq(PCAtools::findElbowPoint(Stdev(seurat)^2)))
seurat<-FindClusters(seurat, resolution = 0.1)
seurat<-RunUMAP(seurat,dims = seq(PCAtools::findElbowPoint(Stdev(seurat)^2)))
  
  
saveRDS(seurat,"seurat.rds")




gene<-c("PTPRC","CD68","ITGAM","CD163","CD3D","MS4A1","CD79A","PECAM1","RGS5",
        "ISLR","THY1","ACTA2","PDGFRB","PDGFRA","PDPN","COL1A1","COL3A1","DES",
        "FAP","FN1","GREM1","LRAT","NGFR","MSLN","GFAP","MCAM","SYNM")%>%unique

library(homologene)#install.packages("homologene")

human2mouse(gene)->mousegene


pdfpath <- path%>%paste0("/2023.8.17.pdf")
pdf(pdfpath)
DimPlot(seurat,raster = F,label = T,pt.size = 1)+ NoLegend()+ggtitle("GSE163777")->plot;print(plot)
for (variable in mousegene$mouseGene%>%intersect(rownames(seurat))) {
  FeaturePlot(seurat, features = variable,raster = F,order = T,pt.size = 1, cols = c("lightgrey", "magenta"))->plot;print(plot)
};dev.off()
library(pdftools)
res.info <- pdf_info(pdf = pdfpath)
library(magick)
dpi=100
pdf(pdfpath%>%str_remove(".pdf$")%>%paste0("_.pdf"),width = 14,height = 14)
for(pages in 1:res.info$pages){
  pdf_render_page(pdfpath,page = pages,dpi = dpi)%>%image_read%>%plot
};dev.off()



FeaturePlot(seurat, features = "Islr",raster = F,order = T,pt.size = 1)+
  scale_color_gradient2(low = "lightgrey", mid = "lightgrey", high = "magenta",midpoint = 1.5)

FeaturePlot(seurat, features = variable,raster = F,order = T,pt.size = 1)+
  scale_color_gradient2(low = "lightgrey", mid = "lightgrey", high = "magenta",midpoint = (seurat@assays$RNA@data[variable,]%>%max)/3)

FeaturePlot(seurat, features = "Lrat",raster = F,order = T,pt.size = 1)+
  scale_color_gradient2(low = "lightgrey", mid = "lightgrey", high = "magenta",
                        midpoint = (seurat@assays$RNA@data["Lrat",]%>%max)/2)

FeaturePlot(seurat, features = "Islr",raster = F,order = T,pt.size = 1)+
  scale_color_gradient2(low = "lightgrey", mid = "lightgrey", high = "magenta",
                        midpoint = (seurat@assays$RNA@data["Islr",]%>%max)/2)

pdfpath <- path%>%paste0("/2023.8.17_mod.pdf")
pdf(pdfpath)
DimPlot(seurat,raster = F,label = T,pt.size = 1)+ NoLegend()+ggtitle("GSE163777")->plot;print(plot)
for (variable in mousegene$mouseGene%>%intersect(rownames(seurat))) {
  FeaturePlot(seurat, features = variable,raster = F,order = T,pt.size = 1)+
    scale_color_gradient2(low = "lightgrey", mid = "lightgrey", high = "magenta",midpoint = (seurat@assays$RNA@data[variable,]%>%max)/2)->plot;print(plot)
};dev.off()
library(pdftools)
res.info <- pdf_info(pdf = pdfpath)
library(magick)
dpi=100
pdf(pdfpath%>%str_remove(".pdf$")%>%paste0("_.pdf"),width = 14,height = 14)
for(pages in 1:res.info$pages){
  pdf_render_page(pdfpath,page = pages,dpi = dpi)%>%image_read%>%plot
};dev.off()
