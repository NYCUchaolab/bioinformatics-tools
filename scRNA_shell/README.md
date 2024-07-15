# Single cell RNA-seq pipeline

## Workflow
1. **Quality Control :** FastQC & MultiQC -> TrimMomatic, qualimap
2. **Alignment :** 
   - CellRanger
4. **Scaling :**
   - Seurat
5. **Removal of low-quality cells :** SoupX, Scrublet, Seurat
6. **UMI matrix analysis :**
   - Highly variable genes (HVG)
   - PCA, UMAP, tSNE
   - Clustering
   - Cell type identification
   - DE analysis
   - Trajectory analysis


<img src="https://github.com/Juan-Jeffery/bioinformatics-tools/blob/main/scRNA_shell/img/scRNA_pipeline.png" width="600" height="800">

