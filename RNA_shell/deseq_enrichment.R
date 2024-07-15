#### Import Packages ####
setwd("/home/blue/code/rna-seq")
library(readr)
library(DOSE)
library(org.Hs.eg.db)
library(dplyr)
library(DESeq2)
library(clusterProfiler)
library(ggplot2)


#### Data Loading ####
metadata <- read.csv("data/clean/deseq2_metadata.csv", header = TRUE, row.names = 1)
countData_df <- read.csv("data/clean/deseq2_tumor_count_df.csv", header = TRUE, row.names = 1, check.names=FALSE)#
condition <- "PD_SDPR"
contrast <- c("PD_SDPR", "PD", "SDPR")
design <- formula(~PD_SDPR)

#condition <- "PR_SDPD"
#contrast <- c("PR_SDPD", "PR", "SDPD")
#design <- formula(~PR_SDPD)

#### Find DEgenes ####

# create deseq dataset
dds <- DESeqDataSetFromMatrix(countData=countData_df, colData=metadata, design=design)
#dds <- DESeqDataSetFromMatrix(countData=countData_df, design=design)
dds <- DESeq(dds)

# vst and PCA plot
vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup=condition)

# Dispersion plot
plotDispEsts(dds)

# contrast = c('factorName','numeratorLevel','denominatorLevel') -> PD/SDPR
results <- results(dds, contrast=contrast)

sorted_results <- results[order(results$padj), ]
sorted_results_df <- as.data.frame(sorted_results)
sorted_results_df %>% head()

# get HUGO symbol and ENTREZID
tmp_df <- bitr(row.names(sorted_results_df), 
               fromType = "ENSEMBL", 
               toType = c("SYMBOL", "ENTREZID"), 
               OrgDb = org.Hs.eg.db, 
               drop = FALSE)
tmp_df_unique <- tmp_df[!duplicated(tmp_df$ENSEMBL), ] # drop repeated ENSEMBL
sorted_results_df$SYMBOL <- tmp_df_unique$SYMBOL
sorted_results_df$ENTREZID <- tmp_df_unique$ENTREZID
padj_threshold <- 0.05
fold_change_threshold <- 2

# Vaolcano plot color
sorted_results_df <- sorted_results_df %>%
  mutate(color = ifelse(padj < padj_threshold & log2FoldChange > 2, "#FF6688",
                        ifelse(padj < padj_threshold & log2FoldChange < -2, "#6688FF",
                               ifelse(padj < padj_threshold, "#999999", "#DDDDDD"))))


DEgenes_df <- subset(sorted_results_df, padj < padj_threshold & abs(log2FoldChange) > fold_change_threshold)

# Vacano plot 
par(mfrow=c(1,1))
ggplot(sorted_results_df, aes(x = log2FoldChange, y = -log10(padj), color = color)) +
  geom_point(alpha = 0.5) +
  scale_color_identity() +  # 使用預設的顏色
  theme_minimal() +
  labs(x = "log2(Fold Change)", y = "-log10(Adjusted p-value)", title = "Volcano Plot")

# top 5 degenes
DEgenes <- DEgenes_df %>% row.names()
top_genes_df <- DEgenes_df %>% head(10) 
top_genes <- top_genes_df %>% row.names()

# group expression plot
par(mfrow=c(2, 5))
for (i in 1:10) {
  gene_name <- top_genes[i]
  symbol <- top_genes_df$SYMBOL[i]
  plotCounts(dds, gene_name, intgroup = condition, main = paste(gene_name, "\n", symbol))
}

# heatmap
DEgene_count <- countData_df[DEgenes_df %>% row.names(), ]
DEgene_count

#### Enrichment ####

DEgenesENTREZID <- DEgenes_df$ENTREZID %>% na.omit() %>% as.integer()# %>% sort(decreasing = TRUE)

#### GO ####
run_GO <- function(gene, category) {
  enrichGO(gene = gene,
           OrgDb = org.Hs.eg.db,
           keyType = "ENTREZID",
           ont = category,
           pAdjustMethod = "BH",  
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05,
           readable = TRUE)
}
GO_BP<-run_GO(DEgenesENTREZID, "BP")
dotplot(GO_BP, title = "Enrichment GO_BP_dot")
barplot(GO_BP, showCategory=10, title="Enrichment_GO_BP", label_format=100)
cnetplot(GO_BP, colorEdge=TRUE)

#### KEGG ####
kegg_result <- enrichKEGG(DEgenesENTREZID,
                          pAdjustMethod = "BH",  
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05,)

dotplot(kegg_result, title="KEGG Enrichment")
barplot(kegg_result, showCategory=20, title="KEGG Enrichment", label_format=100)
cnetplot(kegg_result, title="KEGG Enrichment", colorEdge=TRUE)

#### DO ####
do_result = enrichDO(DEgenesENTREZID,
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05,
                     pAdjustMethod = "BH",
                     readable = TRUE)
barplot(do_result, showCategory=20, title="DO Enrichment", label_format=100)
