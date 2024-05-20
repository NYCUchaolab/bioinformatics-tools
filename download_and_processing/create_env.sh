# add channel
conda config --add channels bioconda 

# gdc-client (avoid conda problem)
wget https://gdc.cancer.gov/system/files/public/file/gdc-client_v1.6.1_Ubuntu_x64.zip
unzip gdc-client_v1.6.1_Ubuntu_x64.zip
mv gdc-client ../
rm gdc-client_v1.6.1_Ubuntu_x64.zip

# preprocessing
conda create -n wes-preprocessing python=3.8 -y
conda activate wes-preprocessing
conda install mamba -y
mamba install fastqc multiqc gatk4 samtools bedtools bwa picard qualimap trimmomatic pandas requests -y

# variants calling
conda create -n wes-variants-calling -y
conda activate wes-variants-calling
conda install mamba -y
mamba install pindel muse varscan -y