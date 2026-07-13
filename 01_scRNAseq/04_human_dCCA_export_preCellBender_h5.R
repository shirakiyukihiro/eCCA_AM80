# ============================================================
# Human dCCA scRNA-seq — export per-sample 10x matrices as .h5 for CellBender
#
# Study      : Distinct cancer-associated fibroblast lineages and an
#              Am80-based stromal-reprogramming strategy in extrahepatic
#              cholangiocarcinoma (Meflin/ISLR).
# Repository : eCCA_AM80  (https://github.com/shirakiyukihiro/eCCA_AM80)
#
# Purpose    : Read the human dCCA 10x matrices and write per-sample .h5 files (DropletUtils::write10xCounts) as input for CellBender ambient-RNA removal.
# Figure(s)  : Preprocessing for the human scRNA-seq analysis (Fig. 5).
# Input      : data/human_dCCA_scRNAseq/ (raw 10x matrices).
# Output     : data/human_dCCA_scRNAseq/matrix_1..7.h5 (CellBender inputs).
#
# Paths have been made relative for public release: place input files
# in ./data and write outputs to ./output (create these folders, or edit
# the paths below).
# ============================================================


library(tidyverse)
path<-"data/human_dCCA_scRNAseq"

data.table::fread(path%>%paste0("/barcodes.tsv.gz"),header = F,nrows = 1,encoding = "UTF-8")


#
list.files(path)
library(Seurat)
mtx<-Read10X(path)
data.table::fread(path%>%paste0("/features.tsv.gz"),header = F,encoding = "UTF-8")->feature
data.table::fread(path%>%paste0("/barcodes.tsv.gz"),header = F,encoding = "UTF-8")->barcodes

barcodes[33974399:34215268,]

barcodes%>%mutate(SN=V1%>%str_remove("^.+-"))->barcodes

split(barcodes,barcodes$SN)->split_bar

split_bar[[1]]%>%head
mtx[,split_bar[[1]]$V1[1]]
#すべての細胞をひっくるめてcellranger reanalyzeやると3週間たっても終わらなかったので
#バーコードの番号ごとで分けてreanalyzeする
#BiocManager::install("DropletUtils")
DropletUtils::write10xCounts('data/human_dCCA_scRNAseq/matrix_1.h5', 
                             mtx[,split_bar[[1]]$V1], gene.id=feature$V1, 
                             gene.symbol=feature$V2, 
                             barcodes=split_bar[[1]]$V1, version='3')
#保存したものをcellranger reanalyzeへ

DropletUtils::write10xCounts('data/human_dCCA_scRNAseq/matrix_2.h5', 
                             mtx[,split_bar[[2]]$V1], gene.id=feature$V1, 
                             gene.symbol=feature$V2, 
                             barcodes=split_bar[[2]]$V1,version ="3" )


DropletUtils::write10xCounts('data/human_dCCA_scRNAseq/matrix_3.h5', 
                             mtx[,split_bar[[3]]$V1], gene.id=feature$V1, 
                             gene.symbol=feature$V2, 
                             barcodes=split_bar[[3]]$V1,version ="3" )

DropletUtils::write10xCounts('data/human_dCCA_scRNAseq/matrix_4.h5', 
                             mtx[,split_bar[[4]]$V1], gene.id=feature$V1, 
                             gene.symbol=feature$V2, 
                             barcodes=split_bar[[4]]$V1,version ="3" )

DropletUtils::write10xCounts('data/human_dCCA_scRNAseq/matrix_5.h5', 
                             mtx[,split_bar[[5]]$V1], gene.id=feature$V1, 
                             gene.symbol=feature$V2, 
                             barcodes=split_bar[[5]]$V1,version ="3" )

DropletUtils::write10xCounts('data/human_dCCA_scRNAseq/matrix_6.h5', 
                             mtx[,split_bar[[6]]$V1], gene.id=feature$V1, 
                             gene.symbol=feature$V2, 
                             barcodes=split_bar[[6]]$V1,version ="3" )

DropletUtils::write10xCounts('data/human_dCCA_scRNAseq/matrix_7.h5', 
                             mtx[,split_bar[[7]]$V1], gene.id=feature$V1, 
                             gene.symbol=feature$V2, 
                             barcodes=split_bar[[7]]$V1,version ="3" )


warnings()


iconv("<83>h<83><89><83>C<83>u C <82>̃{<83><8a><83><85><81>[<83><80> <83><89><83>x<83><8b><82><cd> OS <82>ł<b7>",
      from = "GB18030", to = "UTF-8")

iconv("<83>h<83><89><83>C<83>u C <82>̃{<83><8a><83><85><81>[<83><80> <83><89><83>x<83><8b><82><cd> OS <82>ł<b7>",
      from = "UTF-16", to = "UTF-8")

iconv("<83>h<83><89><83>C<83>u C <82>̃{<83><8a><83><85><81>[<83><80> <83><89><83>x<83><8b><82><cd> OS <82>ł<b7>",
      from = "GB2312", to = "UTF-8")

iconv("<83>h<83><89><83>C<83>u C <82>̃{<83><8a><83><85><81>[<83><80> <83><89><83>x<83><8b><82><cd> OS <82>ł<b7>",
      from = "UTF-16", to = "GB18030")%>%iconv(from = "GB18030", to = "UTF-8")

iconv("<83>h<83><89><83>C<83>u C <82>̃{<83><8a><83><85><81>[<83><80> <83><89><83>x<83><8b><82><cd> OS <82>ł<b7>",
      from = "UTF-16", to = "Big5")



readLines(path%>%paste0("/matrix.mtx"),n=8)
data.table::fread(path%>%paste0("/matrix.mtx"),header = F,encoding = "UTF-8",skip = 2)->mtx
mtx$V3
