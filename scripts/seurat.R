#!/usr/bin/env Rscript




#	The main goal is to “Cluster the cells” and “Finding differentially expressed genes (cluster biomarkers)” as shown in pbmc-tutorial.Rmd.

#	However, additional things like “a heatmap to examine heterogeneity within/between clusters” will be good.

 

print("Loading libraries")

#install.packages("devtools")
library(devtools)
#install_github("satijalab/seurat", ref = "3bd092a")	#	No. Doesn't include CreateSeuratObject and likely others.
#install.packages("httpuv")
#install.packages("Seurat")
library(Seurat)



#a=read.table("out_cell_readcounts.txt.gz", header=F, stringsAsFactors=F)
#x=cumsum(a$V1)
#x=x/max(x)
#plot(1:length(x), x, type='l', col="blue", xlab="cell barcodes sorted by number of reads [descending]", ylab="cumulative fraction of reads", xlim=c(1,50000))
#plot(1:length(x), x, type='l', col="blue", xlab="cell barcodes sorted by number of reads [descending]", ylab="cumulative fraction of reads", xlim=c(1,5000))


#	from http://satijalab.org/seurat/Seurat_AlignmentTutorial.html

print("Loading data from error_detected.dge.txt.gz")
# load data
ds.data <- read.table("error_detected.dge.txt.gz",row.names=1,header=T)
#ds <- CreateSeuratObject(raw.data = ds.data)
ds <- CreateSeuratObject(raw.data = ds.data, min.cells = 3,  min.genes = 200, is.expr=1)
ds <- NormalizeData(object = ds)
ds <- ScaleData(object = ds)
ds <- FindVariableGenes(object = ds, do.plot = FALSE)
#ds@meta.data[, "protocol"] <- "Whatever"




#	from http://satijalab.org/seurat/pbmc3k_tutorial.html

print("Running PCA")
ds <- RunPCA(object = ds, pc.genes = ds@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)

# Examine and visualize PCA results a few different ways
PrintPCA(object = ds, pcs.print = 1:5, genes.print = 5, use.full = FALSE)

VizPCA(object = ds, pcs.use = 1:2)

PCAPlot(object = ds, dim.1 = 1, dim.2 = 2)

# ProjectPCA scores each gene in the dataset (including genes not included
# in the PCA) based on their correlation with the calculated components.
# Though we don't use this further here, it can be used to identify markers
# that are strongly correlated with cellular heterogeneity, but may not have
# passed through variable gene selection.  The results of the projected PCA
# can be explored by setting use.full=T in the functions above
ds <- ProjectPCA(object = ds, do.print = FALSE)

print("Making a couple Heat Maps")
PCHeatmap(object = ds, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)

PCHeatmap(object = ds, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)

# NOTE: This process can take a long time for big datasets, comment out for
# expediency.  More approximate techniques such as those implemented in
# PCElbowPlot() can be used to reduce computation time
#ds <- JackStraw(object = ds, num.replicate = 100, do.print = FALSE)
#JackStrawPlot(object = ds, PCs = 1:12)
#PCElbowPlot(object = ds)





#	Warning messages:
#	1: In heatmap.2(data.use, Rowv = NA, Colv = NA, trace = "none", col = col.use,  :
#	  Discrepancy: Rowv is FALSE, while dendrogram is `both'. Omitting row dendogram.
#	2: In heatmap.2(data.use, Rowv = NA, Colv = NA, trace = "none", col = col.use,  :
#	  Discrepancy: Colv is FALSE, while dendrogram is `column'. Omitting column dendogram.
#	3: In plot.window(...) : "dimTitle" is not a graphical parameter
#	4: In plot.xy(xy, type, ...) : "dimTitle" is not a graphical parameter
#	5: In title(...) : "dimTitle" is not a graphical parameter
#	Exception in thread "main" java.lang.NumberFormatException: For input string: "1e+05"
#		at java.lang.NumberFormatException.forInputString(NumberFormatException.java:65)
#		at java.lang.Integer.parseInt(Integer.java:580)
#		at java.lang.Integer.parseInt(Integer.java:615)
#		at ModularityOptimizer.readInputFile(ModularityOptimizer.java:184)
#		at ModularityOptimizer.main(ModularityOptimizer.java:74)
#	Error in file(file, "rt") : cannot open the connection
#	Calls: FindClusters -> RunModularityClustering -> read.table -> file
#	In addition: Warning message:
#	In file(file, "rt") :
#	  cannot open file '/Users/jakewendt/.R/output_50527.txt': No such file or directory
#	Execution halted







print("Finding Clusters")

# save.SNN = T saves the SNN so that the clustering algorithm can be rerun
# using the same graph but with a different resolution value (see docs for
# full details)
ds <- FindClusters(object = ds, reduction.type = "pca", dims.use = 1:10, resolution = 0.6, print.output = 0, save.SNN = TRUE)

PrintFindClustersParams(object = ds)

# While we do provide function-specific printing functions, the more general
# function to print calculation parameters is PrintCalcParams().


print("Running TSNE")
ds <- RunTSNE(object = ds, dims.use = 1:10, do.fast = TRUE, check_duplicates = FALSE)
#Error in Rtsne.default(X = as.matrix(x = data.use), dims = dim.embed,  : 
#  Remove duplicates before running TSNE.

# note that you can set do.label=T to help label individual clusters

TSNEPlot(object = ds)



#	No point in this
#save(ds, file = "~/syryu/meddling.R")




#	Error in FindMarkers(object = ds, ident.1 = 1, min.pct = 0.25) : 
#	  No genes pass logfc.threshold threshold
#	Execution halted


print("Finding Markers")

# find all markers of cluster 1
cluster1.markers <- FindMarkers(object = ds, ident.1 = 1, min.pct = 0.25)
print(x = head(x = cluster1.markers, n = 5))

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(object = ds, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
print(x = head(x = cluster5.markers, n = 5))

# find markers for every cluster compared to all remaining cells, report
# only the positive ones
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



print("Ending R script")




#ds.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)

#cluster1.markers <- FindMarkers(object = ds, ident.1 = 0, thresh.use = 0.25, test.use = "roc", only.pos = TRUE)






#	kinda specific to the example data so can't go any further until understand.

# you can plot raw UMI counts as well
#VlnPlot(object = ds, features.plot = c("MS4A1", "CD79A"))
#VlnPlot(object = ds, features.plot = c("NKG7", "PF4"), use.raw = TRUE, y.log = TRUE)
#FeaturePlot(object = ds, features.plot = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"), cols.use = c("grey", "blue"), reduction.use = "tsne")
#
#top10 <- ds.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
# setting slim.col.label to TRUE will print just the cluster IDS instead of
# every cell name
#DoHeatmap(object = ds, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)

