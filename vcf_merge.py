from vcf_utils import *
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import subprocess
pd.set_option('future.no_silent_downcasting', True)

# config
project_dir = "/home/data/dataset/CHOL_10sample/"
os.chdir(project_dir)
variants_calling_dir = "./variants_calling"
output_merged_vcf_path = f"{variants_calling_dir}/all_samples_2callers.merged.vcf"
output_VEP_vcf_path = f"{variants_calling_dir}/all_samples_2callers.vep.vcf"
output_stat_path = f"{variants_calling_dir}/all_samples_2callers.stat.csv"
output_venn_path = f"{variants_calling_dir}/all_samples_2callers.venn.png"
sample_sheet_path = "./gdc_sample_sheet.2024-02-26.tsv"
caller_list = ["muse", "mutect2", "varscan"]

# merge VCF to dataframe
result_df, all_case_df, caller_stat_df =  merge_vcf(sample_sheet_path, variants_calling_dir)

# save files
result_df.to_csv(output_merged_vcf_path, sep="\t", index=False)
caller_stat_df.to_csv(output_stat_path)

venn_plot(all_case_df.loc[:, caller_list], 3, output_venn_path)

# run VEP
def run_command(command):
    process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True, executable='/bin/bash')
    output, error = process.communicate()
    return output, error

cache_version = "103"
dir_cache="/home/data/database/vep"
assembly="GRCh38"
threads=20
assembly_fasta="/home/data/data_thousand/gatk_index/Homo_sapiens_assembly38.fasta"
env_name = "vep"
vep_command = [
    "vep", 
    "-i", output_merged_vcf_path, 
    "-o", output_VEP_vcf_path, 
    "--vcf", 
    "--cache", "--dir_cache", dir_cache, "--cache_version", cache_version, 
    "--assembly", assembly, "--force_overwrite", 
    "--fork", str(threads), 
    "--offline", "--no_stats", "--everything", "--pick", 
    "--fasta", assembly_fasta
]

command = f"source $HOME/miniconda3/etc/profile.d/conda.sh; conda run -n {env_name} {' '.join(vep_command)}"
output, error = run_command(command)
print(output, error)
