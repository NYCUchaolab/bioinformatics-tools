from vcf_utils import *
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import subprocess
pd.set_option('future.no_silent_downcasting', True)

# config
project_dir = "/home/data/dataset/CHOL_10sample/"
os.chdir(project_dir)
output_merged_vcf_path = "./variants_calling/all_samples_2callers.merged.vcf"
output_stat_path = "./variants_calling/all_samples_2callers.stat.csv"
output_venn_path = "./variants_calling/all_samples_2callers.venn.png"
sample_sheet_path = "./gdc_sample_sheet.2024-02-26.tsv"
sample_sheet = pd.read_csv(sample_sheet_path, sep="\t")
unique_case_IDs = sample_sheet.loc[:, "Case ID"].unique() 

all_case_df = pd.DataFrame()
all_case_stat = pd.DataFrame()
caller_list = ["muse", "mutect2", "varscan"]
result_dict = {}
caller_stat_df = pd.DataFrame(columns=["all"] 
                                    + [f"{caller}_raw" for caller in caller_list]
                                    + [f"{caller}_filtered" for caller in caller_list]
                                    + ["2_caller"],
                              index=pd.Index(unique_case_IDs, name="case_ID"))
caller_stat_df.index.name = "case_ID"

for case_ID in unique_case_IDs: # "TCGA-W5-AA2R" # progressing bar
    
    result_dict["muse"] = preprocessing_vcf(vcf_file=f"./variants_calling/muse/{case_ID}.MuSE.vcf")
    result_dict["mutect2"] = preprocessing_vcf(vcf_file=f"./variants_calling/mutect/{case_ID}.splited.vcf")
    varscan_indel = preprocessing_vcf(vcf_file=f"./variants_calling/varscan/{case_ID}.varscan.indel.Somatic.hc.vcf")
    varscan_snp = preprocessing_vcf(vcf_file=f"./variants_calling/varscan/{case_ID}.varscan.snp.Somatic.hc.vcf")
    result_dict["varscan"] = pd.concat([varscan_indel, varscan_snp])

    # merge dataframes
    case_df = pd.DataFrame()
    for package_name in result_dict:
        
        sub_df = result_dict[package_name]
        caller_stat_df.loc[case_ID, f"{package_name}_raw"] = len(sub_df)
        filtered_df = vcf_df_filter(sub_df)
        caller_stat_df.loc[case_ID, f"{package_name}_filtered"] = len(filtered_df)
        sub_df = pd.DataFrame(True, index=sub_df.index, columns=[package_name])
        case_df = pd.concat([case_df, sub_df], axis=1)
    
    case_df = case_df.fillna(False).infer_objects(copy=False)
    caller_2hit_filter = case_df.loc[:, caller_list].sum(axis=1) >= 2
    caller_stat_df.loc[case_ID, "all"] = len(case_df)
    caller_stat_df.loc[case_ID, "2_caller"] = caller_2hit_filter.sum()
    
    case_df["case"] = case_ID
    all_case_df = pd.concat([all_case_df, case_df])
    print(f"{case_ID} is merged")

# all_case_df
all_case_df_gb = all_case_df.groupby(all_case_df.index)
all_case_df.groupby("case").size()

# merge same ID
all_case_df_dedup = all_case_df_gb.agg({'muse':'max',
                                          'mutect2':'max',
                                          'varscan':'max',
                                          'case': lambda x: ','.join(x)})
callers_2_hit = all_case_df_dedup.loc[:, caller_list].sum(axis=1) >= 2
all_case_df_2_callers = all_case_df_dedup[callers_2_hit] 
splited_info = all_case_df_2_callers.reset_index()['index'].str.split("|", expand=True)
splited_info.columns = ["#CHROM", "POS", "REF",  "ALT"]
splited_info['callers'] = all_case_df_2_callers.apply(lambda row: ';'.join([caller for caller, value in row[caller_list].items() if value]), axis=1).values
splited_info['case'] = all_case_df_2_callers['case'].values

custom_order = ['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY']
sorter = ['chr'+str(i) for i in range(1, 23)] + ['chrX', 'chrY']
splited_info['#CHROM'] = splited_info['#CHROM'].astype('category')
splited_info['#CHROM'] = splited_info['#CHROM'].cat.set_categories(sorter)
splited_info = splited_info.sort_values(by=['#CHROM', 'POS'])
splited_info = splited_info.reset_index(drop=True)   

# save files
splited_info.to_csv(output_merged_vcf_path, sep="\t", index=False)
caller_stat_df.to_csv(output_stat_path)
venn_plot(all_case_df.loc[:, caller_list], 
          3, 
          output_venn_path)
