setwd("/Users/jeffery/Desktop/生醫資料處理/single_cell/")
#install.packages('Seurat')
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(grid)
library(cowplot)

#Load data
data <- Read10X_h5("filtered_feature_bc_matrix.h5") #cellranger output
data <- CreateSeuratObject(counts = data, project = "NK",min.features = 200)
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
head(data@meta.data, 5)

#original
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0.0001)
plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt") + 
  theme(legend.position="none")
plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  theme(legend.position="none")
combined_plot <- plot1 + plot2
print(combined_plot)

#filter
data_filter <- subset(data, subset = nFeature_RNA < 6000 & nCount_RNA < 25000 & percent.mt < 10)
VlnPlot(data_filter, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0.0001)
plot1 <- FeatureScatter(data_filter, feature1 = "nCount_RNA", feature2 = "percent.mt") + 
  theme(legend.position="none")
plot2 <- FeatureScatter(data_filter, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  theme(legend.position="none")
combined_plot <- plot1 + plot2
print(combined_plot)

#Log normalization
data_nor <- NormalizeData(data_filter) 

#find HVGs
HVGs <- FindVariableFeatures(data_nor, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
top10 <- head(VariableFeatures(HVGs), 10)
grid.newpage()
plot1 <- VariableFeaturePlot(HVGs) + 
  theme(legend.position="top")
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0) + 
  theme(legend.position="none")
combined_plot <- plot1 + plot2
print(combined_plot)

#Scale data
#HVGs
all.genes <- rownames(HVGs)
HVGs_scale <- ScaleData(HVGs, features = all.genes)

#PCA
HVGs_PCA <- RunPCA(HVGs_scale, features = VariableFeatures(object = HVGs_scale))
HVGs_PCA <- FindNeighbors(HVGs_PCA, dims = 1:10)

#Find cluster by K-means
HVGs_PCA$kmeans_5 <- kmeans(x = HVGs_PCA@reductions[["pca"]]@cell.embeddings, centers = 5)$cluster
HVGs_PCA$kmeans_10 <- kmeans(x = HVGs_PCA@reductions[["pca"]]@cell.embeddings,centers = 10)$cluster
HVGs_PCA$kmeans_15 <- kmeans(x = HVGs_PCA@reductions[["pca"]]@cell.embeddings,centers = 15)$cluster
plot_grid(ncol = 3,
          DimPlot(HVGs_PCA, reduction = "pca", group.by = "kmeans_5")+ggtitle("kmeans_5"),
          DimPlot(HVGs_PCA, reduction = "pca", group.by = "kmeans_10")+ggtitle("kmeans_10"),
          DimPlot(HVGs_PCA, reduction = "pca", group.by = "kmeans_15")+ggtitle("kmeans_15")
)

#Find clusters by Louvain
HVGs_PCA_Louvain <- FindClusters(HVGs_PCA, resolution = 0.5)
DimPlot(HVGs_PCA_Louvain, reduction = "pca")+ggtitle("HVGs_PCA_Louvain")
HVGs_UMAP_Louvain <- RunUMAP(HVGs_PCA_Louvain, dims = 1:10)
DimPlot(HVGs_UMAP_Louvain, reduction = "umap")+ggtitle("HVGs_UMAP_Louvain")
HVGs_TSNE_Louvain <- RunTSNE(HVGs_PCA_Louvain, dims = 1:10)
DimPlot(HVGs_TSNE_Louvain, reduction = "tsne")+ggtitle("HVGs_TSNE_Louvain")

####cell type annotation####
library(SingleR)
library(SingleCellExperiment)
library(celldex)

#Load reference dataset for cellde
ref <- DatabaseImmuneCellExpressionData()
ref2 <- BlueprintEncodeData()
ref3 <- HumanPrimaryCellAtlasData()

####umap####
#Convert Seurat object to SingleCellExperiment object for SingleR
HVGs_sce_umap <- as.SingleCellExperiment(HVGs_UMAP_Louvain)

#SingleR annotation 
pred.HVGs_annotations_umap <- SingleR(test = HVGs_sce_umap, ref = ref, labels = ref$label.main)
pred.HVGs_annotations_umap2 <- SingleR(test = HVGs_sce_umap, ref = ref2, labels = ref2$label.main)
pred.HVGs_annotations_umap3 <- SingleR(test = HVGs_sce_umap, ref = ref3, labels = ref3$label.main)

#Add SingleR annotations to the Seurat object
HVGs_UMAP_Louvain$singleR_celltypes_dice <- pred.HVGs_annotations_umap$pruned.labels
HVGs_UMAP_Louvain$singleR_celltypes_blueprint <- pred.HVGs_annotations_umap2$pruned.labels
HVGs_UMAP_Louvain$singleR_celltypes_hpca <- pred.HVGs_annotations_umap3$pruned.labels

# Plot UMAP with SingleR annotations
DimPlot(HVGs_UMAP_Louvain, reduction = "umap", group.by = "singleR_celltypes_dice") + ggtitle("HVGs UMAP with SingleR dice Annotations")
DimPlot(HVGs_UMAP_Louvain, reduction = "umap", group.by = "singleR_celltypes_blueprint") + ggtitle("HVGs UMAP with SingleR blueprint Annotations")
DimPlot(HVGs_UMAP_Louvain, reduction = "umap", group.by = "singleR_celltypes_hpca") + ggtitle("HVGs UMAP with SingleR hpca Annotations")

####tsne####
#Convert Seurat object to SingleCellExperiment object for SingleR
HVGs_sce_tsne <- as.SingleCellExperiment(HVGs_TSNE_Louvain)

#SingleR annotation 
pred.HVGs_annotations_tsne <- SingleR(test = HVGs_sce_tsne, ref = ref, labels = ref$label.main)
pred.HVGs_annotations_tsne2 <- SingleR(test = HVGs_sce_tsne, ref = ref2, labels = ref2$label.main)
pred.HVGs_annotations_tsne3 <- SingleR(test = HVGs_sce_tsne, ref = ref3, labels = ref3$label.main)

#Add SingleR annotations to the Seurat object
HVGs_TSNE_Louvain$singleR_celltypes_dice <- pred.HVGs_annotations_tsne$pruned.labels
HVGs_TSNE_Louvain$singleR_celltypes_blueprint <- pred.HVGs_annotations_tsne2$pruned.labels
HVGs_TSNE_Louvain$singleR_celltypes_hpca <- pred.HVGs_annotations_tsne3$pruned.labels

# Plot tSNE with SingleR annotations
DimPlot(HVGs_TSNE_Louvain, reduction = "tsne", group.by = "singleR_celltypes_dice") + ggtitle("HVGs tSNE with SingleR dice Annotations")
DimPlot(HVGs_TSNE_Louvain, reduction = "tsne", group.by = "singleR_celltypes_blueprint") + ggtitle("HVGs tSNE with SingleR blueprint Annotations")
DimPlot(HVGs_TSNE_Louvain, reduction = "tsne", group.by = "singleR_celltypes_hpca") + ggtitle("HVGs tSNE with SingleR hpca Annotations")

####scDEA####
library("scDEA")
#Convert Seurat object to SingleCellExperiment object for SingleR
raw_counts <- as.matrix(HVGs_PCA[["RNA"]]$"scale.data")+1#看能不能count+1
raw_counts[raw_counts < 0] <- 0
#raw_counts1 <- data@assays[["RNA"]]@layers[["counts"]]
#raw_counts2 <- as.matrix(HVGs_PCA_Louvain@assays[["RNA"]]@layers[["counts"]])
#sce <- SingleCellExperiment(assays = list(counts = raw_counts))
cell_label <- HVGs_TSNE_Louvain@meta.data[["singleR_celltypes_hpca"]]
cell_label_1_other <- ifelse(!is.na(cell_label) & cell_label == "T_cells", 1, 0)
cell_label_1_other <- as.vector(cell_label_1_other)

raw_counts_subset <- as.matrix(raw_counts[,1:100])
cell_label_subset <- as.vector(cell_label_1_other[1:100])

Pvals <- scDEA_individual_methods(raw.count = raw_counts_subset, cell.label = cell_label_subset,
                                  BPSC.parallel = FALSE, DEsingle.parallel = FALSE, DESeq2 = FALSE,
                                  MAST.parallel = FALSE, scDD = FALSE, zingeR.edgeR = FALSE)
combination.Pvals <- lancaster.combination(Pvals, weight = TRUE, trimmed = 0.2)
adjusted.Pvals <- scDEA.p.adjust(combination.Pvals, adjusted.method = "bonferroni")


library(pheatmap)
library(matrixStats)

pval_matrix <- as.matrix(Pvals[,-1])
rownames(pval_matrix) <- Pvals[,1]
num_cols <- ncol(pval_matrix)
#pearson
library(Hmisc)
pearson_corr <- rcorr(pval_matrix, type = "pearson")
pearson_corr_matrix <- pearson_corr$r
dissimilarity_matrix <- pearson_corr_matrix
# Plot the heatmap
pheatmap(dissimilarity_matrix, 
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         main = "Dissimilarity Heatmap of DE Methods",
         color = colorRampPalette(c("blue", "white", "red"))(50),
         show_rownames = TRUE,  # Show row names
         show_colnames = TRUE)  # Show column names




