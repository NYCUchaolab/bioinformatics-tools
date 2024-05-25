
library(maftools)
maf <- read.maf("/home/blue/code/TKI-multiomics/TKI_metadata - 工作表9 (1).tsv")
#maf <- read.maf("/home/blue/code/TKI-multiomics/NTU_TKI_26TB_samples_2callers_001_vep103_maftools.maf")

summary(maf)# Assuming 'maf' is your loaded MAF data

# Summary plot
plotmafSummary(maf, dashboard = TRUE, showBarcodes=TRUE)

# Oncoplot
oncoplot(maf, top = 20, showTumorSampleBarcodes = TRUE)

# TMB
tmb <- ttmb(maf)

somaticInteractions(maf, top=20)


# 設置圖形的寬度和高度
options(repr.plot.width=80, repr.plot.height=8)
lollipopPlot(
       maf = maf,
       gene = c("EGFR"),
       showMutationRate = TRUE,
       AACol="HGVSp",
       cBioPortal=TRUE,
       collapsePosLabel=FALSE,
       repel = TRUE,
       showDomainLabel=FALSE,
       
       labelPos=c(746, 790, 858, 1047),
       refSeqID = "NM_005228"
       
    )


# Mutation burden plot
mafMutBurden(maf, sampleid = NULL)
