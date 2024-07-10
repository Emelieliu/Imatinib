#----------------------------------------------------------------------------------------------------

# Here's code to load the gene expression data and create SeuratObjects for each dataset
# There's also code to merge the data into one set

# The assumption here is that the files are placed in a folder named "data"

#----------------------------------------------------------------------------------------------------

# Load the gene expression data (Short-read sequenced) for the heat-shock experiment

#SR_hs00.data <- ReadMtx(
 # mtx = "data/Short_read_HS00_counts.mtx", features = "data/Short_read_HS00_features.tsv",
  #cells = "data/Short_read_HS00_barcodes.tsv")

#SR_hs60.data <- ReadMtx(
  #mtx = "data/Short_read_HS60_counts.mtx", features = "data/Short_read_HS60_features.tsv",
  #cells = "data/Short_read_HS60_barcodes.tsv")

#--------------------------------------------------

# Create seuratObjects

#SR_hs00.seurat <- CreateSeuratObject(counts = SR_hs00.data, project = "SR_hs")
#SR_hs00.seurat$stim <- "hs00"

#SR_hs60.seurat <- CreateSeuratObject(counts = SR_hs60.data, project = "SR_hs")
#SR_hs60.seurat$stim <- "hs60"

#----------------------------------------------------------------------------------

library(Seurat)

# Load the gene expression data (Short-read sequenced) for the imatinib experiment

imatinib_ctrl.data <- ReadMtx(
  mtx = "data.tsv/Short_read_imatinib_ctrl_matrix.mtx.gz", features = "data.tsv/Short_read_imatinib_ctrl_features.tsv.gz",
  cells = "data.tsv/Short_read_imatinib_ctrl_barcodes.tsv.gz"
)


imatinib.data <- ReadMtx(
  mtx = "data.tsv/Short_read_imatinib_matrix.mtx.gz", features = "data.tsv/Short_read_imatinib_features.tsv.gz",
  cells = "data.tsv/Short_read_imatinib_barcodes.tsv.gz"
)

#--------------------------------------------------

# Create seuratObjects

ctrl <- CreateSeuratObject(counts = imatinib_ctrl.data, project = "SR_imatinib")
ctrl$stim <- "ctrl"

imat <- CreateSeuratObject(counts = imatinib.data, project = "SR_imatinib")
imat$stim <- "imatinib"

#----------------------------------------------------------------------------------

# If you want to merge your two SeuratObjects into one object, you can use this:

#heat_shock <- merge(SR_hs00.seurat, y = SR_hs60.seurat, add.cell.ids = c("hs00", "hs60"), project = "heat_shock")

imatinib <- merge(ctrl, y = imat, add.cell.ids = c("ctrl", "imatinib"), project = "imatinib")

#----------------------------------------------------------------------------------


# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
imatinib[["percent.mt"]] <- PercentageFeatureSet(imatinib, pattern = "^MT-")


#PRE-PROCESS WORKFLOW
# Visualize QC metrics as a violin plot, imatinib
  #VlnPlot(imatinib, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#dessa reglerar axlarna för violinplot, subset = ta bort, irrelivanta men inte cutta för tidigt 
#imatinib
  #plot1 <- FeatureScatter(imatinib, feature1 = "nCount_RNA", feature2 = "percent.mt")
  #plot2 <- FeatureScatter(imatinib, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  #plot1 + plot2
imatinib <- subset(imatinib, subset = nFeature_RNA > 1000 & nFeature_RNA < 12500 & percent.mt < 15)


#NORMALIZING DATA
#imatinib
imatinib <- NormalizeData(imatinib)


#FEATURE SELECTION
#kolla plot2, samma som 1 + namn
#imatinib
imatinib <- FindVariableFeatures(imatinib, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(imatinib), 10)
# plot variable features with and without labels
  #plot1 <- VariableFeaturePlot(imatinib)
  #plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  #plot1 + plot2
#sticker ut mycket - statistiskt konstig från resten


#SCALING DATA
#imatinib
all.genes <- rownames(imatinib)
imatinib <- ScaleData(imatinib, features = all.genes)
#imatinib <- ScaleData(imatinib, vars.to.regress = "percent.mt")


#PCA (linear reduction)
#första 5 PC, kolla på 10 gener
#imatinib
imatinib <- RunPCA(imatinib, features = VariableFeatures(object = imatinib))
#print(imatinib[["pca"]], dims = 1:5, nfeatures = 10)

  #VizDimLoadings(imatinib, dims = 1:5, reduction = "pca", ncol = 5)

DimPlot(imatinib, reduction = "pca", group.by = "stim", cols = c("blue", "red")) 
#+ NoLegend()

  #DimHeatmap(imatinib, dims = 1:6, cells = 1500, balanced = TRUE)


#DIMENSIONALITY
#Imatinib
  #ElbowPlot(imatinib, ndims = 30)


#CLUSTER CELLS
#imatinib
imatinib<- FindNeighbors(imatinib, dims = 1:8)
imatinib <- FindClusters(imatinib, resolution = 0.4)
  #head(Idents(imatinib), 5)


#UMAP/tSNE (non-linear reduction)
#dessa copy-paste med coden i cluster cell, ändra resolution i cluster cell för 0,5, och glöm inte att först run den en gång och sedan umap)
#imatinib
imatinib <- RunUMAP(imatinib, dims = 1:8)
#få fram cluster
DimPlot(imatinib, reduction = "umap")
#dela upp i blå och röd beroende på vilket dataset datan kommer ifrån
#DimPlot(imatinib, reduction = "umap", group.by = "stim", cols = c("blue", "red"))


#INSTALL PRESTO --> faster
#install.packages("devtools")
#devtools::install_github("immunogenomics/presto")


imatinib <- JoinLayers(imatinib)

#CLUSTER BIOMARKER
#look at p_val_adj value
#imatinib
#markers in imatinib - red (cluster i)
#clusteri.markers <- FindMarkers(imatinib, ident.1 = c(1, 3, 4, 5))
#head(clusteri.markers, n = 5)

#cluster ctrl, blå
#clusterctrl.markers <- FindMarkers(imatinib, ident.1 = c(0, 2))
#head(clusterctrl.markers, n = 5)

#imat vs imat (rosa5 vs resten imat)
#clusterimat.markers <- FindMarkers(imatinib, ident.1 = c(1, 3, 4), ident.2 = 5)
#head(clusterimat.markers, n = 5)

#ctrl vs imatinib
#clusterc.markers <- FindMarkers(imatinib, ident.1 = c(1, 3, 4, 5), ident.2 = c(0, 2))
#head(clusterc.markers, n = 5)



#FeaturePlot(imatinib, features = c("ABCB1", "ABCG2", "CYP1A2", "BCR", "ABL1", "OIP5", "CTSB", "MUL1", "TNFAIP3", "SLC22A1"))
#FeaturePlot(imatinib, features = c("BCR", "ABL1", "DHRS2", "CALB1", "HBZ", "HBA2", "GAPDH", "HSPA1A"))

#ctr vs imatinib
#FeaturePlot(imatinib, features = c("AL034397.3", "KIFC3", "SLC4A1", "ESPN", "AL034397.2"))

#imatinib vs imatinib
#FeaturePlot(imatinib, features = c("GADD45B", "ERRFI1", "PLCG2", "EGR1", "SGK1"))

#random 
#FeaturePlot(imatinib, features = c("HBZ", "HBA2"))
#FeaturePlot(imatinib, features = c("GAPDH", "HSPA1A"))

library(dplyr)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
imatinib.markers <- FindAllMarkers(imatinib, only.pos = TRUE)
imatinib.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

#heatmap top10 markers
#imatinib.markers %>%
  #group_by(cluster) %>%
  #dplyr::filter(avg_log2FC > 1) %>%
  #slice_head(n = 10) %>%
  #ungroup() -> top10
#DoHeatmap(imatinib, features = top10$gene) + NoLegend()

#heatmap top10 markers
imatinib.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(imatinib, features = top10$gene) 
#+ NoLegend()

# Ridge plots - from ggridges. Visualize single cell expression distributions in each cluster
RidgePlot(pbmc3k.final, features = features, ncol = 2)

imatinib.markers <- FindAllMarkers(imatinib, only.pos = TRUE)
imatinib.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

saveRDS(imatinib, file = "imatinib_merge_data_tillcluster")



