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
