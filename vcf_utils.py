import os
import pandas as pd
import re
import itertools
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn3

def sample_type(barcode):
    # example: get "01" from TCGA-W5-A22X-01A
    # 00-09 tumor; 10-19 normal; 20-29 control
    match = re.search(r'-([0-2][0-9])[A-Z]*', barcode)
    if match:
        code = int(match.group(1))
        if code < 10:
            return "TUMOR"
        elif code < 20:
            return "NORMAL"
        elif code < 30:
            return "CONTROL"
        else:
            raise SyntaxError("Its not a TCGA barcode format.")

def open_vcf(vcf_file):
    # Find the line number containing string "#CHROM"
    with open(vcf_file, 'r') as f:
        for line_ind, line in enumerate(f):
            if line.startswith('#CHROM'):
                break
    
    # load data and skip headers
    df = pd.read_csv(vcf_file, sep='\t', header=line_ind, dtype=object)
    return df

def parse_info(format, info_str):
    '''
    GT: 基因型 (Genotype)
    DP: 深度 (Depth)
    AD: 等位基因深度 (Allelic Depth)，通常是指各等位基因的深度
    BQ: 鹼基品質 (Base Quality)
    AF = FREQ: Allele Frequency(等位基因頻率) = AD/DP
    RD: Reference Depth(参考序列深度)
    
    varscan, GT:GQ:DP:RD:AD:FREQ:DP4, 0/0:.:68:65:2:2.99%:15,50,1,1
    mutect2, GT:AD:AF:DP:F1R2:F2R1, 0/1:12,4:0.275:16:8,2:4,2
    mutect2, GT:AD:AF:DP:F1R2:F2R1, 0/1/2:13,3,2:0.181,0.145:18:3,1,1:10,2,1
    muse, GT:DP:AD:BQ:SS, 0/1:19:14,5:31,31:2

    * mutect 2 multiallelic問題已經在前處理中使用bcftools norm拆分 *
    '''
    info_list = info_str.split(":")
    # varscan format
    if format == "GT:GQ:DP:RD:AD:FREQ:DP4":
        DP_index, AD_index, AF_index = 2, 4, 5
        DP, AD = int(info_list[DP_index]), int(info_list[AD_index])
        AF = float(info_list[AF_index].replace("%", ""))/100
        
    # mutect2 format
    elif format == "GT:AD:AF:DP:F1R2:F2R1" or format == "GT:AD:AF:DP:F1R2:F2R1:PGT:PID:PS":
        DP_index, AD_index, AF_index = 3, 1, 2
        AD = int(info_list[AD_index].split(",")[1])
        DP, AF = int(info_list[DP_index]), float(info_list[AF_index])
        
    # muse format
    elif format == "GT:DP:AD:BQ:SS":
        DP_index, AD_index = 1, 2
        AD = int(info_list[AD_index].split(",")[1])
        DP = int(info_list[DP_index])
        AF = AD / DP
    else:
        raise ValueError(f"{format} is an Unexpected format")
    return DP, AD, AF

def preprocessing_vcf(vcf_file):
    vcf_df = open_vcf(vcf_file)
    vcf_df = vcf_df.loc[vcf_df.FILTER == "PASS"] # keep "PASS"
    
    # keep only chr1-22 or XY
    vcf_df = vcf_df[vcf_df['#CHROM'].str.contains(r'^chr(?:[1-9]|1\d|2[0-2]|X|Y)$')]
    
    # make sure TUMOR/NORMAL is upper
    vcf_df.columns = ['CHROM'] + list(vcf_df.columns[1:].str.upper()) 
    
    # reindex to CHROM|POS|REF|ALT and drop 7 columns
    vcf_df.index = vcf_df.loc[:, ["CHROM", "POS", "REF", "ALT"]].apply(lambda row: "|".join(row.astype(str)), axis=1)
    vcf_df = vcf_df.drop(columns=["CHROM", "POS", "REF", "ALT", 'ID', 'QUAL', 'FILTER'])
    
    # if no columns TUMOR/NORMAL(mutect2) rename last 2 columns
    if "TUMOR" not in vcf_df.columns or "NORMAL" not in vcf_df.columns :
        vcf_df = vcf_df.rename(columns={vcf_df.columns[-1]: sample_type(vcf_df.columns[-1]), 
                                        vcf_df.columns[-2]: sample_type(vcf_df.columns[-2])})
    # extract DP AD AF columns from TUMOR columns
    vcf_df[['DP', 'AD', 'AF']] = vcf_df.apply(lambda row: parse_info(row['FORMAT'], row['TUMOR']), axis=1, result_type='expand')
    
    return vcf_df

def vcf_df_filter(vcf_df, DP=30, AD=3, AF=0.05):
    filter = (vcf_df.DP >= DP) & (vcf_df.DP >= AD) & (vcf_df.AF > AF)
    return vcf_df[filter]

def venn_plot(caller_df, n_callers, plot_path):
    subset_sizes = {}
    for comb in itertools.product([0, 1], repeat=n_callers):
        filters = np.where(np.array(comb) == 1, True, False)
        subset_sizes["".join(map(str, comb))] = caller_df.loc[:, filters].all(axis=1).sum()
    plt.figure(figsize=(8, 10))
    venn3(subsets=subset_sizes, set_labels=caller_df.columns)
    plt.title(f'Venn Diagram: {", ".join(caller_df.columns)}')
    plt.savefig(plot_path)