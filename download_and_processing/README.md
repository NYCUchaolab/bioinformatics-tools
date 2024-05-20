# Data Download and Preprocessing

## TCGA (BAM)

1. **Download Sample Sheet**  
   Retrieve the sample sheet from the TCGA database.

2. **Replace Sample Type**  
   Modify the sample sheet to represent sample types with characters T (Tumor), N (Normal), or C (Control).

3. **Download and Rename BAM Files**  
   Download the BAM files as specified in the sample sheet and rename them accordingly.

4. **Convert BAM to FASTQ**  
   Convert the BAM files to FASTQ format using appropriate tools.

## Other Datasets (FASTQ)

1. **Create Sample Sheet**  
   Generate a sample sheet for the dataset.

2. **Rename FASTQ Files**  
   Rename the FASTQ files as per the sample sheet.

## Quality Control

1. **Run FastQC**  
   Perform quality control on the FASTQ files using FastQC.

2. **Aggregate with MultiQC**  
   Aggregate FastQC reports using MultiQC for a comprehensive overview.

3. **Check Base Quality and Length**  
   Review base quality scores and sequence lengths from the QC reports.

4. **Trim Reads if Necessary**  
   Trim the FASTQ files to remove low-quality bases or adapters if needed.
