#"Medulloblastoma scRNA-seq Data Analysis"
#Luz Ruiz
#February 20, 2025

#Set working directory

setwd("~pathtoworkingdirectory")

#Set up libraries
library(dplyr)
library(Seurat)
library(future)
library(patchwork)
library(ggplot2)
library(cowplot)
library(clustree)
#install.packages("ggpubr")
library(ggpubr)

#Set up options
options(future.globals.maxSize = 2048 * 1024^2)

#Load in metadata from sample series matrix file
series_matrix <- read.delim(file = "GSE119926_series_matrix.txt", skip = 32)
series_matrix

#Optimize the titles
rownames(series_matrix) <- paste0(rownames(series_matrix), "_", series_matrix$X.Sample_title)
series_matrix$X.Sample_title <- NULL
series_mat <- as.data.frame(t(series_matrix))

#Determine the heading of series_matrix
rownames(series_mat)

#Load in Medulloblastoma Patient Tumor Samples

#file locations
file_list <- c(
  "BCH807"  = "GSM3905406_BCH807.txt",
  "BCH825"  = "GSM3905407_BCH825.txt",
  "BCH1031" = "GSM3905408_BCH1031.txt",
  "BCH1205" = "GSM3905409_BCH1205.txt",
  "MUV11"   = "GSM3905410_MUV11.txt",
  "MUV19"   = "GSM3905411_MUV19.txt",
  "MUV27"   = "GSM3905412_MUV27.txt",
  "MUV29"   = "GSM3905413_MUV29.txt",
  "MUV34"   = "GSM3905414_MUV34.txt",
  "MUV37"   = "GSM3905415_MUV37.txt",
  "MUV39"   = "GSM3905416_MUV39.txt",
  "MUV41"   = "GSM3905417_MUV41.txt",
  "MUV44"   = "GSM3905418_MUV44.txt",
  "SJ17"    = "GSM3905419_SJ17.txt",
  "SJ99"    = "GSM3905420_SJ99.txt",
  "SJ129"   = "GSM3905421_SJ129.txt",
  "SJ217"   = "GSM3905422_SJ217.txt",
  "SJ454"   = "GSM3905423_SJ454.txt",
  "SJ516"   = "GSM3905424_SJ516.txt",
  "SJ577"   = "GSM3905425_SJ577.txt",
  "SJ617"   = "GSM3905426_SJ617.txt",
  "SJ625"   = "GSM3905427_SJ625.txt",
  "SJ723"   = "GSM3905428_SJ723.txt",
  "SJ917"   = "GSM3905429_SJ917.txt",
  "SJ970"   = "GSM3905430_SJ970.txt"
)
#custom function called patient_counts creates individual seurat objects for all the patient samples and load their metadata
patient_counts <- function(sample_id, series_mat, file_list) {
  file_name <- file_list[sample_id]  
  counts <- read.table(file = file_name, header = TRUE, sep = "\t", row.names = 1)
  obj <- CreateSeuratObject(counts = counts)
  obj$Methylation_Subgroup <- series_mat[sample_id, ]$`11_!Sample_characteristics_ch1`
  obj$Methylation_Subtype <- series_mat[sample_id, ]$`12_!Sample_characteristics_ch1`
  obj$Metastasis <- series_mat[sample_id, ]$`13_!Sample_characteristics_ch1`
  obj$Histology <- series_mat[sample_id, ]$`14_!Sample_characteristics_ch1`
  obj$ID <- sample_id
  return(obj)
}

#patient names are saved in the patient_id variable 
patient_id <- names(file_list) 

#Store all the objects in MedList
MedList <- lapply(patient_id, patient_counts, series_mat = series_mat, file_list = file_list)
names(MedList) <- patient_id
MedList

# Perform the pre-processing workflow
for (i in 1:length(MedList)){
  MedList[[i]]$percent.mt <- PercentageFeatureSet(MedList[[i]], pattern = "^MT")
}

head(MedList[[1]]$percent.mt)
FeatureScatter(object = MedList[[2]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(object = MedList[[2]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend()

FeatureScatter(object = MedList[[2]], feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(object = MedList[[2]], feature1 = "nCount_RNA", feature2 = "percent.mt") + NoLegend()

# Normalize the data and find variable features
for (i in 1:length(x = MedList)) {
  MedList[[i]] <- FindVariableFeatures(object = MedList[[i]],
                                       selection.method = "vst", nfeatures = 2000, verbose = TRUE)
}

for (i in 1:length(x = MedList)) {
  MedList[[i]] <- NormalizeData(object = MedList[[i]],
  )
}


Med.Anchors <- FindIntegrationAnchors(object.list = MedList, dims = c(1:30),
                                      k.anchor = 5,
                                      k.filter = 15,
                                      k.score = 30)

Med.Int <- IntegrateData(anchorset = Med.Anchors, dims = 1:30, verbose = TRUE, k.weight = 40)

# Run the standard workflow for visualization and clustering
plan(strategy = "multicore", workers = 6)

Med.Int <- ScaleData(Med.Int, features = rownames(Med.Int))

Med.Int <- RunPCA(object = Med.Int, npcs = 30, verbose = FALSE)
ElbowPlot(Med.Int)

# Find cell clusters
Med.Int <- FindNeighbors(Med.Int, dims = 1:15)
#Med.Int <- FindClusters(Med.Int, resolution = 0.475)
Med.Int <- FindClusters(Med.Int, resolution = 0.5)
#Med.Int <- FindClusters(Med.Int, resolution = seq(from = 5.0, to = 0.1, by = -0.025,  algorithm = 2))

# Run non-linear dimensional reduction
Med.Int <- RunUMAP(object = Med.Int, reduction = "pca",
                   dims = 1:15)
unique(Med.Int$Methylation_Subgroup)

DefaultAssay(Med.Int) <- "integrated"
Med.Int <- CellCycleScoring(Med.Int, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes, set.ident = TRUE)
Med.Int$CellCycleIdents <- Idents(Med.Int)
DimPlot(Med.Int, reduction = "umap", label.size = 50)

# Run the standard workflow for visualization and clustering
plan(strategy = "multicore", workers = 6)

#Med.Int <- ScaleData(Med.Int, features = rownames(Med.Int))
Med.Int <- ScaleData(Med.Int, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(Med.Int))

Med.Int <- RunPCA(object = Med.Int, npcs = 30, verbose = FALSE)
ElbowPlot(Med.Int)

# Plot figures of patient data (not anchor integrated)
Med.Int <- FindVariableFeatures(object = Med.Int, selection.method = "vst", nfeatures = 2000, verbose = TRUE)
ElbowPlot(Med.Int)

# Identify all markers
DefaultAssay(Med.Int) <- "integrated"
Idents(Med.Int) <- Med.Int$integrated_snn_res.0.5
Idents(Med.Int)
Med.Int <- ScaleData(Med.Int, features = rownames(Med.Int))
Med.Int.markers <- FindAllMarkers(Med.Int, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")
Med.Int.Topmarkers <- Med.Int.markers %>% group_by(cluster) %>% top_n(1, wt = avg_log2FC)
p1 <- DoHeatmap(object = Med.Int, features = Med.Int.Topmarkers$gene)
p1


Med.Int.markers$gene <- rownames(Med.Int.markers)
Med.Int.markers$gene

# Load in non-tumor brain cells
data <- read.table(file = "Darmanis_et_al_2015_matrix.txt", sep = "\t")
hb <- CreateSeuratObject(counts = data)

hb <- NormalizeData(object = hb)
hb <- FindVariableFeatures(object = hb)
hb <- ScaleData(object = hb)
hb <- RunPCA(object = hb)
ElbowPlot(hb)
hb <- FindNeighbors(object = hb)
hb <- FindClusters(object = hb, resolution = seq(from = 0.00, to = 0.5, by = 0.1))

hb <- RunUMAP(object = hb, dims = 1:20)
DimPlot(object = hb, reduction = "umap")

hb$clusters <- hb$RNA_snn_res.0.5
Idents(hb) <- hb$clusters

# Identify markers in non-tumor brain cells
hb.markers <- FindAllMarkers(hb, test.use = "MAST", only.pos = TRUE)
top.markers <- hb.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top.markers
DoHeatmap(hb, features = top.markers$gene)

# Rename cell identities to the cell types from Darmanis 2017 paper
hb <- RenameIdents(hb, "0" = "Neural", "1" = "Neural", "2" = "Astrocytes", 
                   "3" = "Endothelial", "4" = "Neural", "5" = 
                     "Oligodendrocytes", "6" = "Microglia")
hb$CellTypes <- Idents(hb)

hb <- RunPCA(object = hb, npcs = 30, verbose = FALSE)
hb <- RunUMAP(object = hb, reduction = "pca",
              dims = 1:15)
unique(hb$Methylation_Subgroup)
hb <- CellCycleScoring(hb, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes, set.ident = TRUE)
hb$CellCycleIdents <- Idents(hb)
DimPlot(hb, reduction = "umap")

hb$clusters <- hb$RNA_snn_res.0.5
Idents(hb) <- hb$clusters

DoHeatmap(hb, features = top.markers$gene)

### Integrate this seurat object with your tumor data, and make disease signatures
### by differential expression between tumor cells and these cells (should be comparable
### with other SMARTseq based datasets)

MB.list <- c(Med.Int, hb)

# # And uniformly filter all datasets
for (i in 1:length(MB.list)) {
  MB.list[[i]] <- subset(x = MB.list[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 4500)
}

# Run scTransform on all objects
for (i in 1:length(MB.list)) {
  MB.list[[i]] <- SCTransform(MB.list[[i]], verbose = TRUE, return.only.var.genes	= FALSE, assay = "RNA")
}

# Next, select features for downstream integration, and run PrepSCTIntegration,
# which ensures that all necessary Pearson residuals have been calculated.
options(future.globals.maxSize = 3000 * 1024^2)

MBint.features <- SelectIntegrationFeatures(
  object.list = MB.list,
  nfeatures = 3000
)

MB.list <- PrepSCTIntegration(
  object.list = MB.list,
  anchor.features = MBint.features,
  verbose = TRUE)

# Next, identify anchors and integrate the datasets to compare MB and non-tumor brain cells. 
MBint.anchors <- FindIntegrationAnchors(
  object.list = MB.list,
  normalization.method = "SCT",
  anchor.features = MBint.features,
  verbose = TRUE
)

to_integrate <- Reduce(intersect, lapply(MBint.anchors@object.list, rownames))

MB_hb_int <- IntegrateData(
  anchorset = MBint.anchors,
  features.to.integrate = to_integrate,
  normalization.method = "SCT",
  verbose = TRUE
)

Idents(MB_hb_int)
unique(MB_hb_int$clusters)
MB_hb_int <- RenameIdents(MB_hb_int, "Neural" = "Healthy Cells", "Astrocytes" = "Healthy Cells", 
                          "Endothelial" = "Healthy Cells","Neural" = "Healthy Cells","Oligodendrocytes" = "Healthy Cells","Microglia" = "Healthy Cells")
Idents(MB_hb_int)
MB_hb_int <- RenameIdents(MB_hb_int, "Neural" = "Healthy Cells", "Astrocytes" = "Healthy Cells", 
                          "Endothelial" = "Healthy Cells","Neural" = "Healthy Cells","Oligodendrocytes" = "Healthy Cells","Microglia" = "Healthy Cells")

# To compare with non-tumor brain cells in the integrated data, we need to make sure the NAs are listed as healthy cells
MB_hb_int$isna <- is.na(MB_hb_int$integrated_snn_res.0.5)
MB.HB <- data.frame(MB_hb_int$isna, MB_hb_int$integrated_snn_res.0.5)
MB.HB
MB_hb_int$isna == TRUE
colnames(MB.HB)

MB.HB$NewClusters <- "TEST"
for (i in 1:length(rownames(MB.HB))){
  if (MB.HB$MB_hb_int.isna[i] == TRUE){
    MB.HB$neoplasticvsnon[i] <- "Healthy Cells"
  }
  if (MB.HB$MB_hb_int.isna[i] == FALSE){
    MB.HB$NewClusters[i] <- MB.HB$MB_hb_int.integrated_snn_res.0.5[i]
    MB.HB$neoplasticvsnon[i] <- MB.HB$MB_hb_int.integrated_snn_res.0.5[i]
  }
}

Idents(MB_hb_int) <- MB.HB$neoplasticvsnon
Idents(MB_hb_int)
MB_hb_int <- AddMetaData(MB_hb_int, metadata = MB.HB$neoplasticvsnon, col.name = "SNN_NeoplasticVsNon")

Idents(MB_hb_int) <- MB_hb_int$SNN_NeoplasticVsNon

DefaultAssay(MB_hb_int) <- "RNA"
MB_hb_int <- ScaleData(MB_hb_int, features = rownames(MB_hb_int))
MB_hb_int.markers <- FindAllMarkers(MB_hb_int, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")
MB_hb_int.Topmarkers <- MB_hb_int.markers %>% group_by(cluster) %>% top_n(2, wt = avg_log2FC)
DoHeatmap(object = MB_hb_int, features = unique(MB_hb_int.Topmarkers$gene), label = F)
#mb_hb_int_disease_sig <- MB_hb_int.markers[c("gene", "avg_log2FC")]
#head(mb_hb_int_disease_sig)

DefaultAssay(MB_hb_int) <- "integrated"
MB_hb_int <- RunPCA(object = MB_hb_int, npcs = 30, verbose = FALSE)
MB_hb_int <- RunUMAP(object = MB_hb_int, reduction = "pca",
                     dims = 1:15)
unique(MB_hb_int$Methylation_Subgroup)
DefaultAssay(MB_hb_int) <- "integrated"
MB_hb_int <- CellCycleScoring(MB_hb_int, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes, set.ident = TRUE)
MB_hb_int$CellCycleIdents <- Idents(MB_hb_int)
DimPlot(MB_hb_int, reduction = "umap", order = c(10,9,8,7,6,5,4,3,2,1))
plan(strategy = "multicore", workers = 6)
MB_hb_int <- ScaleData(MB_hb_int, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(MB_hb_int))

MB_hb_int <- RunPCA(object = MB_hb_int, npcs = 30, verbose = FALSE)
ElbowPlot(MB_hb_int)
MB_hb_int <- FindVariableFeatures(object = MB_hb_int, selection.method = "vst", nfeatures = 2000, verbose = TRUE)
ElbowPlot(MB_hb_int)
DefaultAssay(MB_hb_int) <- "integrated"

Idents(MB_hb_int) <- MB_hb_int$SNN_NeoplasticVsNon

MB_hb_int <- AddMetaData(MB_hb_int, metadata = tumorSNNdat, col.name = "OldCluster")
MB_hb_int$OldCluster
MB_hb_int$isna <- is.na(MB_hb_int$OldCluster)

MB.HB <- data.frame(MB_hb_int$isna, MB_hb_int$OldCluster, MB_hb_int$CellTypes)
MB.HB
# MB_hb_int$isna == TRUE

MB.HB$tumorClusters<- "TEST"
for (i in 1:length(rownames(MB.HB))){
  if (MB.HB$MB_hb_int.isna[i] == TRUE){
    MB.HB$tumorClusters[i] <- MB.HB$MB_hb_int.CellTypes[i]
    MB.HB$neoplasticvsnon[i] <- "Healthy Cells"
  }
  if (MB.HB$MB_hb_int.isna[i] == FALSE){
    MB.HB$tumorClusters[i] <- MB.HB$MB_hb_int.OldCluster[i]
    MB.HB$neoplasticvsnon[i] <- MB.HB$MB_hb_int.OldCluster[i]
  }
}
# colnames(tumorSNNdat) <-  MB_hb_int 
#add this to metadata for int obj (w/ healthy cells)

MB_hb_int <- AddMetaData(MB_hb_int, metadata = MB.HB$neoplasticvsnon, col.name = "TumorVsHealthy")
#Cluster_0_diseaseSig <- FindMarkers(MB_hb_int, ident.1 = "0", ident.2 = "healthyCells")
Idents(MB_hb_int) <- MB_hb_int$TumorVsHealthy
head(Idents(MB_hb_int))

# Identify markers of MB clusters compared to non-tumor brain cells
Cluster0 <- FindMarkers(MB_hb_int,  ident.1 = "0", ident.2 = "Healthy Cells", test.use = "MAST", logfc.threshold = 0.25, min.pct = 0.25, only.pos = T)
Cluster1 <- FindMarkers(MB_hb_int,  ident.1 = "1", ident.2 = "Healthy Cells", test.use = "MAST", logfc.threshold = 0.25, min.pct = 0.25, only.pos = T)
Cluster2 <- FindMarkers(MB_hb_int,  ident.1 = "2", ident.2 = "Healthy Cells", test.use = "MAST", logfc.threshold = 0.25, min.pct = 0.25, only.pos = T)
Cluster3 <- FindMarkers(MB_hb_int,  ident.1 = "3", ident.2 = "Healthy Cells", test.use = "MAST", logfc.threshold = 0.25, min.pct = 0.25, only.pos = T)
Cluster4 <- FindMarkers(MB_hb_int,  ident.1 = "4", ident.2 = "Healthy Cells", test.use = "MAST", logfc.threshold = 0.25, min.pct = 0.25, only.pos = T)
Cluster5 <- FindMarkers(MB_hb_int,  ident.1 = "5", ident.2 = "Healthy Cells", test.use = "MAST", logfc.threshold = 0.25, min.pct = 0.25, only.pos = T)
Cluster6 <- FindMarkers(MB_hb_int,  ident.1 = "6", ident.2 = "Healthy Cells", test.use = "MAST", logfc.threshold = 0.25, min.pct = 0.25, only.pos = T)
Cluster7 <- FindMarkers(MB_hb_int,  ident.1 = "7", ident.2 = "Healthy Cells", test.use = "MAST", logfc.threshold = 0.25, min.pct = 0.25, only.pos = T)
Cluster8 <- FindMarkers(MB_hb_int,  ident.1 = "8", ident.2 = "Healthy Cells", test.use = "MAST", logfc.threshold = 0.25, min.pct = 0.25, only.pos = T)
Cluster9 <- FindMarkers(MB_hb_int,  ident.1 = "9", ident.2 = "Healthy Cells", test.use = "MAST", logfc.threshold = 0.25, min.pct = 0.25, only.pos = T)

#need to add data on what cluster each result row is from
Cluster0$cluster <- "0"
Cluster0$gene <- rownames(Cluster0)
Cluster1$cluster <- "1"
Cluster1$gene <- rownames(Cluster1)
Cluster2$cluster <- "2"
Cluster2$gene <- rownames(Cluster2)
Cluster3$cluster <- "3"
Cluster3$gene <- rownames(Cluster3)
Cluster4$cluster <- "4"
Cluster4$gene <- rownames(Cluster4)
Cluster5$cluster <- "5"
Cluster5$gene <- rownames(Cluster5)
Cluster6$cluster <- "6"
Cluster6$gene <- rownames(Cluster6)
Cluster7$cluster <- "7"
Cluster7$gene <- rownames(Cluster7)
Cluster8$cluster <- "8"
Cluster8$gene <- rownames(Cluster8)
Cluster9$cluster <- "9"
Cluster9$gene <- rownames(Cluster9)

diseaseSigDataFrame <- rbind(Cluster0, Cluster1, Cluster2, Cluster3, Cluster4, Cluster5, Cluster6, Cluster7, Cluster8, Cluster9) 

#diseaseSigDataFrame$gene <- rownames(diseaseSigDataFrame)
unique(diseaseSigDataFrame$cluster)
head(diseaseSigDataFrame)

## Analysis for data visualization shown in figure of the paper is shown below 

head(Idents(MB_hb_int))

MB_hb_int_data <- GetAssayData(object = MB_hb_int, assay = "RNA")

MB_hb_int_BAIAP2_CDC42 <- subset(MB_hb_int, features = c("BAIAP2", "CDC42"))
MB_hb_int_BAIAP2_CDC42 <- GetAssayData(object = MB_hb_int, assay = "RNA")

Idents(MB_hb_int) <- MB_hb_int$Methylation_Subgroup
Idents(MB_hb_int)

DimPlot(MB_hb_int,  label = T)
MB_hb_int <- RenameIdents(MB_hb_int, "methylation subgroup: WNT" = "WNT", "methylation subgroup: Group 3" = "Group 3",
                          "methylation subgroup: Group 4" = "Group 4", "methylation subgroup: SHH-infant" = "SHH-infant",
                          "methylation subgroup: SHH-adult" = "SHH-adult")
Idents(MB_hb_int)
Idents(MB_hb_int) -> MB_hb_int$subgroup

Idents(MB_hb_int) <- MB_hb_int$CellTypes
Idents(MB_hb_int)

#Next add clusters so labels are each subgroup then non-neoplastic and use find markers to get markers

Idents(MB_hb_int) <- MB_hb_int$Methylation_Subgroup
Idents(MB_hb_int) 

# Create a new 'subgroupvcelltypes' column that labels cells by either their subgroup or cell type
MB_hb_int$subgroupvcelltypes <- ifelse(is.na(MB_hb_int$subgroup), 
                                       MB_hb_int$CellTypes, 
                                       MB_hb_int$subgroup)

# Add 'subgroupvcelltypes' to the metadata of MB_hb_int object
MB_hb_int <- AddMetaData(MB_hb_int, metadata = MB_hb_int$subgroupvcelltypes, col.name = "subgroupvcelltypes")

# Set the identity of MB_hb_int to the 'subgroupvcelltypes' column
Idents(MB_hb_int) <- MB_hb_int$subgroupvcelltypes

# Check the new identities
head(Idents(MB_hb_int))

# Ensure that 'subgroup' and 'CellTypes' are character vectors
MB_hb_int$subgroup <- as.character(MB_hb_int$subgroup)
MB_hb_int$CellTypes <- as.character(MB_hb_int$CellTypes)

# Use coalesce to prioritize non-NA values from 'subgroup', fallback to 'CellTypes' if NA
MB_hb_int$subgroupvcelltypes <- coalesce(MB_hb_int$subgroup, MB_hb_int$CellTypes)

# Add 'subgroupvcelltypes' to the metadata of MB_hb_int object
MB_hb_int <- AddMetaData(MB_hb_int, metadata = MB_hb_int$subgroupvcelltypes, col.name = "subgroupvcelltypes")

# Set the identity of MB_hb_int to the 'subgroupvcelltypes' column
Idents(MB_hb_int) <- MB_hb_int$subgroupvcelltypes

# Check the new identities
head(Idents(MB_hb_int))
DimPlot(MB_hb_int, label = T)

MB_hb_int <- RunUMAP(MB_hb_int, dims = 1:5, min.dist =	
                       0.1, n.neighbors = 25)

# Visualize BAIAP2 and CDC42 expression
FeaturePlot(MB_hb_int, features = c("BAIAP2", "CDC42"), blend = T)
DimPlot(MB_hb_int)
FeaturePlot(MB_hb_int, features = c("BAIAP2"))
FeaturePlot(MB_hb_int, features = c("CDC42"))
FeaturePlot(MB_hb_int, features = c("BAIAP2"), label = T)
FeaturePlot(MB_hb_int, features = c("CDC42"), label = T)
VlnPlot(MB_hb_int, features = "BAIAP2")
VlnPlot(MB_hb_int, features = "CDC42")
DotPlot(MB_hb_int, features = c("BAIAP2", "CDC42"))

# Group cells by either neoplastic or non-neoplastic

MB_hb_int <- RenameIdents(MB_hb_int, "Group 3" = "Neoplastic", "Group 4" = "Neoplastic", 
                          "SHH-infant" = "Neoplastic","SHH-adult" = "Neoplastic", "WNT" = "Neoplastic",
                          "Oligodendrocytes" = "Non-Neoplastic","Neural" = "Non-Neoplastic",
                          "Endothelial" = "Non-Neoplastic","Microglia" = "Non-Neoplastic",
                          "Astrocytes" = "Non-Neoplastic")

MB_hb_int$NeovsNon <- Idents(MB_hb_int)

table(MB_hb_int$NeovsNon) #3793 neo and 279 non-neoplastic

FeaturePlot(MB_hb_int, features = c("BAIAP2"), split.by = "NeovsNon")
FeaturePlot(MB_hb_int, features = c("CDC42"), split.by = "NeovsNon")
FeatureScatter(MB_hb_int, feature1 = "BAIAP2", feature2 = "CDC42", split.by = "NeovsNon")

MB_hb_int$NeovsNon -> Idents(MB_hb_int)

# Visualize BAIAP2 and CDC42 expression
FeaturePlot(MB_hb_int, features = c("BAIAP2", "CDC42"), blend = T)
DimPlot(MB_hb_int, cols = c("magenta", "black"))
VlnPlot(MB_hb_int, features = "BAIAP2", cols = c("magenta", "black"))
VlnPlot(MB_hb_int, features = "CDC42", cols = c("magenta", "black"))
DotPlot(MB_hb_int, features = c("BAIAP2", "CDC42"))


b <- VlnPlot(MB_hb_int, features = "BAIAP2", cols = c("magenta", "black"))
c <- VlnPlot(MB_hb_int, features = "CDC42", cols = c("magenta", "black"))

# Run MAST test for differential expression in neoplastic vs non-neoplastic
markers <- FindMarkers(MB_hb_int, ident.1 = "Neoplastic", ident.2 = "Non-Neoplastic", test.use = "MAST")
write.csv(markers, file = "markers.20250218.csv")

#BAIAP2 p-val adj 2.06008579481588e-30
#CDC42 p-val adj 8.60336839866911e-06

# To do differential expression analysis by subgroup, group cells by either subgroup or non-neoplastic
Idents(MB_hb_int) <- MB_hb_int$subgroupvcelltypes

MB_hb_int <- RenameIdents(MB_hb_int, 
                          "Oligodendrocytes" = "Non-Neoplastic","Neural" = "Non-Neoplastic",
                          "Endothelial" = "Non-Neoplastic","Microglia" = "Non-Neoplastic",
                          "Astrocytes" = "Non-Neoplastic", "SHH-infant" = "SHH", "SHH-adult" = "SHH")

MB_hb_int$SubgroupvsNonneo <- Idents(MB_hb_int)

SHH.markers <- FindMarkers(MB_hb_int, ident.1 = "SHH", ident.2 = "Non-Neoplastic", test.use = "MAST", assay = "RNA")
WNT.markers <- FindMarkers(MB_hb_int, ident.1 = "WNT", ident.2 = "Non-Neoplastic", test.use = "MAST", assay = "RNA")
G3.markers <- FindMarkers(MB_hb_int, ident.1 = "Group 3", ident.2 = "Non-Neoplastic", test.use = "MAST", assay = "RNA")
G4.markers <- FindMarkers(MB_hb_int, ident.1 = "Group 4", ident.2 = "Non-Neoplastic", test.use = "MAST", assay = "RNA")


write.csv(SHH.markers, file = "SHH.markers.20250219.csv")
write.csv(WNT.markers, file = "WNT.markers.20250219.csv")
write.csv(G3.markers, file = "G3.markers.20250219.csv")
write.csv(G4.markers, file = "G4.markers.20250219.csv")

# save object as an RDS file
saveRDS(MB_hb_int, file = "MB_hb_int.20250219.BAIAP2paper.RDS")
#MB_hb_int <- readRDS("MB_hb_int.20250218.BAIAP2paper.RDS")
