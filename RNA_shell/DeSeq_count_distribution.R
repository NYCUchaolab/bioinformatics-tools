setwd("/Users/jeffery/Desktop/生醫資料處理/RNAII")
library(readr)
library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)
library(dplyr)
library(DESeq2)
library(tibble)
#count
countData_df<-read_csv("countData_for_DesSeq.csv")
countData_df<-as.data.frame(countData_df)
#change row name
row.names(countData_df) <- countData_df[, 1]
countData_df <- countData_df[, -1]
countData_df<-t(countData_df)
#col
colData_df<-read_csv("colData_for_DesSeq.csv")
colData_df<-as.data.frame(colData_df)
#change row name
rownames(colData_df) <- colData_df$Sample
colData_df <- colData_df[, "Condition", drop = FALSE]
colnames(colData_df)[1]<-"type"

#dds
dds <- DESeqDataSetFromMatrix(countData=countData_df, 
                              colData=colData_df, 
                              design=~factor(type))
dds <- DESeq(dds)
# Take a look at the results table
res <- results(dds)
head(results(dds, tidy=TRUE))
# Summary of differential gene expression
summary(res)
# Sort summary list by p-value
res <- res[order(res$padj),]
head(res)
DEGs <- row.names(subset(res, padj<.01 & abs(log2FoldChange)>2)) 
res_df <- as.data.frame(res)
filtered_res <- res_df[row.names(res_df) %in% DEGs, ]
DEG_LFC_df <- rownames_to_column(filtered_res, var = "gene") %>%
  select(gene, log2FoldChange)

####PCA####
vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="type")
#plotCounts
#par(mfrow=c(2,3))
#top6<-rownames(res)[1:6]
#plotCounts(dds,top6, intgroup="type")
####Volcano Plot####
par(mfrow=c(1,1))
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3), ylim=c(0, 40), col="grey"))
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

####Functional enrichment####
######Transfer ID type
DEG_ENTREZID_df <- bitr(DEG_LFC_df$gene, fromType = "ENSEMBL", toType = "ENTREZID"
                     , OrgDb = org.Hs.eg.db)
colnames(DEG_LFC_df)[1] <- "ENSEMBL"
merged_df <- merge(DEG_ENTREZID_df, DEG_LFC_df, by = "ENSEMBL", all.x = TRUE)
merged_df <- merged_df[!duplicated(merged_df$ENTREZID), ]
DEG_ENTREZID <- merged_df$ENTREZID
rownames(merged_df)<-merged_df[,2]
merged_df <- merged_df %>%
  select(log2FoldChange)
#####enrichGO#######
run_GO <- function(gene,category) {
  enrichGO(gene = gene,
           OrgDb = org.Hs.eg.db,
           keyType = "ENTREZID",
           ont = category,
           pAdjustMethod = "BH", 
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.2,
           readable = TRUE)
}
###
GO_BP<-run_GO(DEG_ENTREZID,"BP")
dotplot(GO_BP, title = "Enrichment GO_BP_dot")
barplot(GO_BP, showCategory=20, title="Enrichment_GO_BP",label_format=100)
goplot(GO_BP)
GO_CC<-run_GO(DEG_ENTREZID,"CC")
dotplot(GO_CC, title = "Enrichment GO_CC_dot")
barplot(GO_CC, showCategory=20, title="EnrichmentGO_CC",label_format=100)
goplot(GO_CC)
GO_MF<-run_GO(DEG_ENTREZID,"MF")
dotplot(GO_MF, title = "Enrichment GO_MF_dot")
barplot(GO_MF, showCategory=20, title="EnrichmentGO_MF",label_format=100)
goplot(GO_MF)

# KEGG pathway
kk<-enrichKEGG(gene=DEG_ENTREZID,organism = 'hsa',pvalueCutoff = 0.05)
dotplot(kk,title="Enrichment KEGG_dot")
merged_df$ENTREZID <- as.numeric(merged_df$ENTREZID)
merged_df$log2FoldChange <- as.numeric(merged_df$log2FoldChange)
print(head(kk@result$ID))
hsa04080 <- pathview(gene.data = merged_df, pathway.id = "04080", species = "hsa")
