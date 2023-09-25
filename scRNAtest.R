#TESTING SOME FUN STUFF
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#load libraries
library(Seurat)
library(SeuratDisk)
library(tidyverse)

# Load the NSCLC dataset
nsclc.sparse.m <- Read10X_h5(filename = '20k_NSCLC_DTC_3p_nextgem_Multiplex_count_raw_feature_bc_matrix.h5')

#shows us what info is in the data that we uploaded
str(nsclc.sparse.m) 

#here we have 3 modalities, we will take a loot at gene expression
cts <-  nsclc.sparse.m$`Gene Expression`
cts[1:10,10]

# Initialize the Seurat object with the raw (non-normalized data).
nsclc.seurat.obj <- CreateSeuratObject(counts = cts,
                                       project = "NSCLC", #just the name
                                       min.cells = 3,     #keep all features that have show up in at least 3 cells
                                       min.features = 200)#keep all cells with at least 200 features 
str(nsclc.seurat.obj)
nsclc.seurat.obj
# 29552 features across 42081 samples


# 1. QC ------- filtering out low quality cells. 
#want to remove outrageouly low and outrageously high reads. this will remove broken cells, doublets. Before we removed the blank droplets.
View(nsclc.seurat.obj@meta.data)

# % MT reads
#here we are looking at the percentage of mitochondria genes. (mt genes)
nsclc.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(nsclc.seurat.obj, pattern = "^MT-")
View(nsclc.seurat.obj@meta.data)

VlnPlot(nsclc.seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) #these are the cells that need to be sorted out.

#This represents the quality of cells. A good quality datasets should follow the straight line. 
#Why does this happen? 
#If the genes are below the line (bottom right) it means the same genes are being captured and sequenced over and over again.
#If the genes are above the line (top left) it means your experiment is capturing a high number of genes, but they aren't being sequenced enough.
#In both cases, these cells are not of good quality.
FeatureScatter(nsclc.seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')

# 2. Filtering -----------------
nsclc.seurat.obj <- subset(nsclc.seurat.obj,
                           subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
#we can also use DoubleFinder to filter out doublets. 

#Here is the plot once the filtering has been done. 
#still some low quality cells, but a lot better. 
VlnPlot(nsclc.seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(nsclc.seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')


# 3. Normalize data ----------
#nsclc.seurat.obj <- NormalizeData(nsclc.seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
# OR
nsclc.seurat.obj <- NormalizeData(nsclc.seurat.obj) #the values above are default.
str(nsclc.seurat.obj) #see the command slot. you can see what has been done to the data.


# 4. Identify highly variable features -------------- This is DE genes. However, data is still not scaled, gonna be shit!
nsclc.seurat.obj <- FindVariableFeatures(nsclc.seurat.obj, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(nsclc.seurat.obj), 10)# Identify the 10 most highly variable genes

plot1 <- VariableFeaturePlot(nsclc.seurat.obj)
LabelPoints(plot = plot1, points = top10, repel = TRUE)# plot variable features with and without labels


# 5. Scaling ------------- many unwanted sources of variation, want to remove these biases.
all.genes <- rownames(nsclc.seurat.obj)
nsclc.seurat.obj <- ScaleData(nsclc.seurat.obj, features = all.genes)

str(nsclc.seurat.obj) 
#Under the assay slot we have "counts", "data", "scale.data".
#This stores all previous sets that we had. We can just access them again!

# 6. Perform Linear dimensionality reduction --------------
nsclc.seurat.obj <- RunPCA(nsclc.seurat.obj, features = VariableFeatures(object = nsclc.seurat.obj))

# visualize PCA results. 
print(nsclc.seurat.obj[["pca"]], dims = 1:5, nfeatures = 5) #The first 5 features of each.
DimHeatmap(nsclc.seurat.obj, dims = 1, cells = 500, balanced = TRUE) #a heat map of 500 cells for PC_1. 
#dims = will change the PC. 
#cells = changes number of cells.

# determine dimensionality of the data
ElbowPlot(nsclc.seurat.obj) #shows which dimension we should stop at. 
#compare PC1 compared to PC20.
DimHeatmap(nsclc.seurat.obj, dims = 20, cells = 500, balanced = TRUE) #clearly not a lot of separation compared to 1. 


# 7. Clustering ------------
nsclc.seurat.obj <- FindNeighbors(nsclc.seurat.obj, dims = 1:15)

# understanding resolution
nsclc.seurat.obj <- FindClusters(nsclc.seurat.obj, resolution = c(0.1, 0.3, 0.5, 0.7, 1)) 
# the lower the number the fewer clusters. the higher number the more clusters. 
View(nsclc.seurat.obj@meta.data) #shows the names of columns for the next step. 

#looking at the clusters.
DimPlot(nsclc.seurat.obj, group.by = "RNA_snn_res.0.5", label = TRUE)

# setting identity of clusters
Idents(nsclc.seurat.obj)
Idents(nsclc.seurat.obj) <- "RNA_snn_res.0.1"
Idents(nsclc.seurat.obj)

# non-linear dimensionality reduction --------------
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
nsclc.seurat.obj <- RunUMAP(nsclc.seurat.obj, dims = 1:15)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(nsclc.seurat.obj, reduction = "umap") #uses the Indents value. 

