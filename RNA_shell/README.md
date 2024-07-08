# DNA variant calling pipeline

## Workflow
1. **Quality Control :** FastQC & MultiQC -> TrimMomatic, qualimap
2. **Alignment :** 
   - Genome Mapping : STAR, HISAT2 + Samtools
   - Pseudo Mapping : Kallisto, salmon
4. **Expression Quantification :**
   - Gene-level:
      - Using GTF/GFFfile : featureCounts
   - Transcriptome-level:
      - genome-mapping : stringTie
      - Transcriptome-mapping : RSEM
      - Alignment-free : salmon, kallisto
   After completing Expression Quantification, draw the express distribution graph.
5. **Data preprocessing :**
   - Remove the Batch effect
   - Outlier detection
   - Within/Between sample normalization
6. **Differential expression analysis :**
   - Methods & tools : limma, edgeR, DESeq
7. **Functional profiling :**
   – map genes onto pathways
   - Methods : ORA, GSEA
   - Database : KEGG, GO

<img src="https://github.com/Juan-Jeffery/bioinformatics-tools/blob/main/RNA_shell/img/RNA_pipeline.png" width="600" height="600">

## Reference file
**GRCh38.d1.vd1.fa & gencode.v36.annotation.gtf**

The ref and ref index file can be downloaded from [GDC Website](https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files).

If you don't want to install the `.fai`, you can also index them yourself.

**gencode.v22.annotation.gtf**

The ref and ref index file can be downloaded from [GDC Website](https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files).

**Homo_sapiens.GRCh38.cdna.all.fa.gz (ref for Kallistal & Salmon)**

This file can be downloaded from [ensemble Website](https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/cdna/).

**grch38_genome.tar.gz (ref for HISAT2)**
This file can be downloaded from [ensemble Website](https://genome-idx.s3.amazonaws.com/hisat/)

## Package version
Make sure to verify and install the required packages by reviewing the `.yml` file.

**RNA_seq**
  ```bash
# Name                    Version                   Build  Channel
fastqc                    0.12.1               hdfd78af_0    bioconda 
multiqc                   1.17                     pypi_0    pypi
qualimap                  2.2.2a                        1    bioconda
------------------------------------------------------------------------
hisat2                    2.2.1                hdbdd923_6    bioconda
samtools                  1.19.2               h50ea8bc_1    bioconda
star                      2.7.11b              h43eeafb_1    bioconda
------------------------------------------------------------------------
salmon                    0.14.2               ha0cc327_0    bioconda
kallisto                  0.50.1               hc877fd6_1    bioconda
------------------------------------------------------------------------
rsem                      1.3.3          pl5321h0033a41_7    bioconda
subread (featureCounts)   2.0.6                he4a0461_0    bioconda
stringtie                 2.2.1                h43eeafb_6    bioconda
  ```
