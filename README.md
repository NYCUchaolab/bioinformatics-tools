# Bioinformatics Tools

## TODO
### Download and preprocessing
- [x]  download bam from GDC
- [x]  select samples by sample sheet
- [ ]  download fastq from SRA
### WES
- [ ]  env file
- [x]  bam → fastq sh
- [x]  fastqc & multiqc
- [ ]  trim fastq (Trimmomatic)
- [ ]  fastq → bam (BWA-mem)
- [x]  vcf merge python
    - [x]  preprocessing
    - [x]  merge
    - [x]  status table
- [x]  vep python
- [x]  vep vcf → maf
- [x]  maf → maftool format maf
### RNA-seq
- [ ]  env file
- [x]  make indices (RSEM, kallisto)
- [x]  fastq → abundance table (kallisto, K-mer base pseudo alignment method)
- [x]  fastq → bam (STAR)
- [ ]  bam → count matrix (RSEM)
- [ ]  visualize count matrix

### Methylation
- [ ]  preprocessing
- [ ]  DMP
- [ ]  DMR 

### Others
- [ ]  config
- [x]  template
- [x]  send message python script(Telegram)
