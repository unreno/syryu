#!/usr/bin/env Rscript




#	The main goal is to “Cluster the cells” and “Finding differentially expressed genes (cluster biomarkers)” as shown in pbmc-tutorial.Rmd.

#	However, additional things like “a heatmap to examine heterogeneity within/between clusters” will be good.

 



#install.packages("devtools")
library(devtools)
#install_github("satijalab/seurat", ref = "3bd092a")
#install.packages("httpuv")
#install.packages("Seurat")
library(Seurat)



#a=read.table("out_cell_readcounts.txt.gz", header=F, stringsAsFactors=F)
#x=cumsum(a$V1)
#x=x/max(x)
#plot(1:length(x), x, type='l', col="blue", xlab="cell barcodes sorted by number of reads [descending]", ylab="cumulative fraction of reads", xlim=c(1,50000))
#plot(1:length(x), x, type='l', col="blue", xlab="cell barcodes sorted by number of reads [descending]", ylab="cumulative fraction of reads", xlim=c(1,5000))


#	from http://satijalab.org/seurat/Seurat_AlignmentTutorial.html


ds.data <- read.table("error_detected.dge.txt.gz",row.names=1,header=T)
ds <- CreateSeuratObject(raw.data = ds.data)
ds <- NormalizeData(object = ds)
ds <- ScaleData(object = ds)
ds <- FindVariableGenes(object = ds, do.plot = FALSE)
#ds@meta.data[, "protocol"] <- "Whatever"




#	from http://satijalab.org/seurat/pbmc3k_tutorial.html

ds <- RunPCA(object = ds, pc.genes = ds@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)

PrintPCA(object = ds, pcs.print = 1:5, genes.print = 5, use.full = FALSE)

VizPCA(object = ds, pcs.use = 1:2)

PCAPlot(object = ds, dim.1 = 1, dim.2 = 2)

ds <- ProjectPCA(object = ds, do.print = FALSE)

PCHeatmap(object = ds, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)

PCHeatmap(object = ds, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)

#ds <- JackStraw(object = ds, num.replicate = 100, do.print = FALSE)
#JackStrawPlot(object = ds, PCs = 1:12)
#PCElbowPlot(object = ds)










ds <- FindClusters(object = ds, reduction.type = "pca", dims.use = 1:10, resolution = 0.6, print.output = 0, save.SNN = TRUE)

PrintFindClustersParams(object = ds)

ds <- RunTSNE(object = ds, dims.use = 1:10, do.fast = TRUE, check_duplicates = FALSE)
#Error in Rtsne.default(X = as.matrix(x = data.use), dims = dim.embed,  : 
#  Remove duplicates before running TSNE.

TSNEPlot(object = ds)

save(ds, file = "~/syryu/meddling.R")





cluster1.markers <- FindMarkers(object = ds, ident.1 = 1, min.pct = 0.25)
print(x = head(x = cluster1.markers, n = 5))

cluster5.markers <- FindMarkers(object = ds, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
print(x = head(x = cluster5.markers, n = 5))

#ds.markers <- FindAllMarkers(object = ds, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
#	Error in rbind(deparse.level, ...) : 
#	  numbers of columns of arguments do not match
#	In addition: Warning messages:
#	1: In data.frame(..., check.names = FALSE) :
#	  row names were found from a short variable and have been discarded
#	2: In data.frame(..., check.names = FALSE) :
#	  row names were found from a short variable and have been discarded
#	3: In data.frame(..., check.names = FALSE) :
#	  row names were found from a short variable and have been discarded


#	This can occur if individual clusters don't have any transcriptomic markers.
#	We'll improve in a future release, but you can either use FindMarkers, or lower the threshold for DE to solve.
#	Trying lower thresh.use
#ds.markers <- FindAllMarkers(object = ds, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.15)
#ds.markers <- FindAllMarkers(object = ds, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.05)
#ds.markers <- FindAllMarkers(object = ds, only.pos = TRUE, min.pct = 0.25 )
#ds.markers <- FindAllMarkers(object = ds, only.pos = TRUE, min.pct = 0.25, thresh.use = 0)
#ds.markers <- FindAllMarkers(object = ds, only.pos = TRUE, min.pct = 0.25, thresh.use = 10)
#
#	Nothing works on this data set






#ds.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)

#cluster1.markers <- FindMarkers(object = ds, ident.1 = 0, thresh.use = 0.25, test.use = "roc", only.pos = TRUE)






#	kinda specific to the example data so can't go any further until understand.
#VlnPlot(object = ds, features.plot = c("MS4A1", "CD79A"))
#VlnPlot(object = ds, features.plot = c("NKG7", "PF4"), use.raw = TRUE, y.log = TRUE)
#FeaturePlot(object = ds, features.plot = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"), cols.use = c("grey", "blue"), reduction.use = "tsne")
#
#top10 <- ds.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
#DoHeatmap(object = ds, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)

