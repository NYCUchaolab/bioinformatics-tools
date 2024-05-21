# Bioinformatics Tools

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
- [x]  config
- [x]  template
- [x]  send message python script (Telegram)
