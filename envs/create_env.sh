# add channel
conda config --add channels bioconda 
conda config --add channels conda-forge 

# download gdc-client (avoid conda problem)
wget https://gdc.cancer.gov/system/files/public/file/gdc-client_v1.6.1_Ubuntu_x64.zip
unzip gdc-client_v1.6.1_Ubuntu_x64.zip
mv gdc-client ../
rm gdc-client_v1.6.1_Ubuntu_x64.zip

# preprocessing
source "/home/${USER}/miniconda3/etc/profile.d/conda.sh" 
conda create -n wes-preprocessing python=3.8 -y
conda activate wes-preprocessing
conda install mamba -y
mamba install fastqc multiqc gatk4 samtools bedtools bwa picard qualimap trimmomatic pandas requests -y

# variants calling
conda create -n wes-variants-calling python=3.6 -y
conda activate wes-variants-calling
conda install mamba -y
mamba install gatk4=4.1 pindel muse varscan bcftools openjdk=8 samtools ensembl-vep=103.1 -y
mamba install pandas matplotlib perl-compress-raw-zlib=2.202 -y
