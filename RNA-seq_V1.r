
library(DESeq2)
require(DESeq2)
library(ggplot2)

setwd ("E:/AMIN HDD/STUDY MATERIALS/UNIVERSITY OF TSUKUBA/R Coding/Micro-regional RNA seq data/KK Sensei data/RNAseq with Waseda/RNAseq with Waseda/Ex015")
require("DESeq2")

just.raw.counts = read.delim("count_matrix.txt")

head(just.raw.counts, n = 20)

dim(just.raw.counts)

just.raw.counts = read.delim("count_matrix.txt", row.names = 1)

head(just.raw.counts)

dim(just.raw.counts)

meta.data = read.delim(file="Meta_data.txt", row.names = 1)

head(meta.data, n = 12)

# making sure the row names in colData matches to column name in counts_data
all(colnames(just.raw.counts) %in% rownames(meta.data))

# are they in the same order?
all(colnames(just.raw.counts) == rownames(meta.data))

count.data.set <- DESeqDataSetFromMatrix(countData=just.raw.counts, 
                                         colData=meta.data, design= ~ condition) 


count.data.set.object <- DESeq(count.data.set)

vsd <- vst(count.data.set.object)

norm.data = assay(vsd)

head(norm.data)

write.table(norm.data, sep="\t",file="Norm_data_all_genes_NO_counts_cut_off.txt", row.names=TRUE,col.names=NA,quote=FALSE)

sampleDists <- dist(t(norm.data),  method = "euclidean")

reversed_rows_columns = (t(norm.data))

reversed_rows_columns[1:5,1:5]

sampleDists

clusters=hclust(sampleDists)

plot(clusters)

plotPCA(vsd, intgroup=c("condition")) 

require(ggplot2)
plotPCA(vsd, intgroup=c("condition")) +
  scale_colour_hue(breaks = c("WT","WT","WT","WT","Homo Normal","Homo Normal","Homo Normal","Homo Normal", "Homo Disease", "Homo Disease", "Homo Disease", "Homo Disease" ))




