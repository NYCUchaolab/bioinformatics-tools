# DNA variant calling pipeline

## Workflow
1. **Quality Control :** FastQC & MultiQC -> TrimMomatic, qualimap
2. **Alignment :** BWA -> GATK
3. **Variant Calling :** 
   - SNP Calling : Mutect2, Varscan2, Muse
   - Indel Calling : Mutect2, Varscan2, Pindel
   - CNV Calling : CNVkit -> GISTIC2
5. **Annotation :** VEP
6. **Format Conversion and Visualization :** vcf2maf, IGV-reports, maftools, SigProFiler, MutsigCV

<img src="https://github.com/Juan-Jeffery/DNA_Variant_Calling_pipeline/blob/main/img/DNA_pipeline.png" width="600" height="600">

## Reference file
**GRCh38.d1.vd1.fa & gencode.v36.annotation.gtf**

The ref and ref index file can be downloaded from [GDC Website](https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files).

If you don't want to install the `.fai` and `.dict` files, you can also index them yourself.

**Other**

- Homo_sapiens_assembly38.dbsnp138.vcf
- af-only-gnomad.hg38.vcf.gz
- 1000g_pon.hg38.vcf.gz

The file can be downloaded from [GATK Website](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle) or [genomics-public-data](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0/).

## Package version
Make sure to verify and install the required packages by reviewing the `.yml` file, excluding MutsigCV.

**WES_preprocessing_SNP (QC + Alinment + Variant_Calling(Mutect,varscan,muse))**
  ```bash
  # Name                    Version                   Build  Channel
  bwa                       0.7.17               h5bf99c6_8    bioconda 
  picard                    2.18.29                       0    bioconda
  samtools                  1.6                  hb116620_7    bioconda 
  ------------------------------------------------------------------------
  fastqc                    0.12.1               hdfd78af_0    bioconda 
  multiqc                   1.17                     pypi_0    pypi
  trimmomatic               0.39                 hdfd78af_2    bioconda
  qualimap                  2.2.2a                        1    bioconda
  ------------------------------------------------------------------------
  gatk4                     4.1.0.0                       0    bioconda 
  muse                      1.0.rc               h2e03b76_5    bioconda
  varscan                   2.4.6                hdfd78af_0    bioconda
  ```
**WES_Indel (Variant_Calling(Pindel))**
  ```bash
  # Name                    Version                   Build  Channel
  pindel                    0.2.5b9             h84372a0_10    bioconda
  ```
**WES_CNV (Variant_Calling(CNVkit))**
  ```bash
  # Name                    Version                   Build  Channel
  cnvkit                    0.9.10             pyhdfd78af_0    bioconda
  gistic2                   2.0.23                        0    hcc
  ```
**WES_Annotation (Annotation(vep, vcf2maf, igv-reports, SigProFiler))**
  ```bash
  # Name                    Version                   Build  Channel
  igv-reports               1.12.0             pyh7cba7a3_0    bioconda
  sigmut                    1.0                  hdfd78af_2    bioconda
  samtools                  1.10                 h2e538c0_3    bioconda
  vcf2maf                   1.6.21               hdfd78af_0    bioconda
  ensembl-vep               103.1          pl5262h4a94de4_2    bioconda 
  ```
