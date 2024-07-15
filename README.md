# Bioinformatics Tools
## ENVS
 - create env file
 - WES environment.yml file
 - RNA environment.yml file
 - scRNA(not done) environment.yml file
   
## Download and processing
 - Method to download bam file from TCGA & SRA (not done)
 - bam to fastq & QC

## WES, RNA, scRNA shell
 - Executable files, from bam to the final step (QC included in Download and Processing)

## Other
 - 不知道是啥

---
## TODO
### Download and preprocessing
- [x]  download bam from GDC
- [x]  select samples by sample sheet
- [ ]  download fastq from SRA
- [x]  env file
- [x]  bam → fastq sh
- [x]  fastqc & multiqc
- [ ]  trim fastq (Trimmomatic)

### WES
- [x]  fastq → bam (BWA-mem)
- [x]  variants calling
    - [x] muse
    - [x] mutect2
    - [x] varscan
    - [x] pindel
- [x]  vcf merge python
    - [ ]  preprocessing
    - [ ]  merge
    - [ ]  status table
- [ ]  vep python
- [ ]  vep vcf → maf
- [ ]  maf → maftool format maf

### RNA-seq
- [x]  env file
- [x]  make indices (RSEM, kallisto)
- [x]  fastq → abundance table (kallisto, K-mer base pseudo alignment method)
- [x]  fastq → bam (STAR)
- [x]  bam → count matrix (RSEM)
- [ ]  visualize count matrix & distribution
  
### Methylation
- [ ]  preprocessing
- [ ]  DMP
- [ ]  DMR 

### Others
- [x]  config
- [x]  template
- [x]  send message python script (Telegram)
