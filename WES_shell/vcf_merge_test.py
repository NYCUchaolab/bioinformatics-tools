import os
import pandas as pd
import re
import itertools
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('../others')
import venn
import argparse
pd.set_option('future.no_silent_downcasting', True)

class VCF(pd.DataFrame):
    # Define chromosome order as a class variable
    chromosome_order = ['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY']
    
    def __init__(self, data=None, variant_index_df: pd.DataFrame = None, *args, **kwargs):
        """Initialize VCF object from various input types.
        
        Args:
            data: Input data, can be:
                - str: Path to VCF file
                - pd.DataFrame: DataFrame containing VCF data
                - None: Empty VCF (requires variant_index_df)
            variant_index_df: Optional pre-formatted variant index DataFrame
            *args, **kwargs: Additional arguments passed to pandas DataFrame
            
        Raises:
            ValueError: If input parameters are invalid or missing
        """
        if variant_index_df is not None and data is not None:
            raise ValueError("Cannot specify both data and variant_index_df. Please choose one.")
            
        # Initialize empty DataFrame first
        super().__init__(*args, **kwargs)
        
        if data is None and variant_index_df is None:
            raise ValueError("Either data or variant_index_df must be provided.")
            
        if variant_index_df is not None:
            # Handle variant index DataFrame
            df = self.variant_index_to_vcf(variant_index_df)
            self._update_inplace(df)
            return
            
        if isinstance(data, str):
            # Handle file path input
            self.file_path = data
            df = self.load()
            self._update_inplace(df)
            self.rename_tumor_normal_column()
            self.parse()
            self.sort()
            
        elif isinstance(data, pd.DataFrame):
            # Handle DataFrame input
            self._update_inplace(data)
        else:
            raise ValueError("data must be either a file path (str) or a DataFrame")
            
    def _update_inplace(self, df: pd.DataFrame) -> None:
        """Update the current DataFrame with new data inplace.
        
        Args:
            df: New DataFrame to update with
        """
        # Update the current DataFrame with the new data
        self._mgr = df._mgr
        
    def load(self):
        # Find the line number containing string "#CHROM"
        with open(self.file_path, 'r') as f:
            for line_ind, line in enumerate(f):
                if line.startswith('#CHROM'):
                    break
        
        # load data and skip headers
        data = pd.read_csv(self.file_path, sep='\t', header=line_ind, dtype=object)
        data.columns = ['CHROM'] + list(data.columns[1:].str.upper())  # format column names

        # filter PASS and non CHR
        data = data[data.FILTER == "PASS"]
        data = data[data['CHROM'].str.contains(r'^chr(?:[1-9]|1\d|2[0-2]|X|Y)$')]  # chr1-22, x, y
        
        # Convert 'CHROM' to category type with ordered categories
        data['CHROM'] = data['CHROM'].astype('category')
        data['CHROM'] = data['CHROM'].cat.set_categories(self.chromosome_order)
        
        # Ensure 'POS' is numeric
        data['POS'] = pd.to_numeric(data['POS'], errors='coerce')

        # Set the DataFrame with the loaded data
        #self.update(data)
        return data
    
    def rename_tumor_normal_column(self):
        if "TUMOR" not in self.columns or "NORMAL" not in self.columns:
            try:
                self.rename(columns={
                    self.columns[-2]: self.sample_type(self.columns[-2]),
                    self.columns[-1]: self.sample_type(self.columns[-1])
                }, inplace=True)
            except:
                self.rename(columns={
                    self.columns[-2]: "NORMAL",
                    self.columns[-1]: "TUMOR"
                }, inplace=True)

    def to_variant_index_df(self):
        variant_index_df = self.copy()
        # Merge columns
        variant_index_df.index = variant_index_df.loc[:, ["CHROM", "POS", "REF", "ALT"]].apply(lambda row: "|".join(row.astype(str)), axis=1)
        variant_index_df = variant_index_df.drop(columns=["CHROM", "POS", "REF", "ALT", 'ID', 'QUAL', 'FILTER'])
        return variant_index_df

    def concat(self, vcf: "VCF"):
        #vcf_copy = copy.deepcopy(self)
        concated_df = pd.concat([self, vcf])
        concated_vcf = VCF(concated_df)
        concated_vcf.sort()
        return concated_vcf
        
    @staticmethod
    def sample_type(sample_ID: str):
        if sample_ID[-1] == "N":
            return "NORMAL"
        elif sample_ID[-1] == "T":
            return "TUMOR"
        else:
            raise SyntaxError("Wrong sample ID format.")
    
    def parse(self):
        self[['DP', 'AD', 'AF']] = self.apply(lambda row: self.parse_info(row['FORMAT'], row['TUMOR']), axis=1, result_type='expand')

    @staticmethod
    def parse_info(format: str, info_str:str):
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
        pindel, GT:AD, 0/0:20,0, 0/0:19,1

        * mutect 2 multiallelic問題已經在前處理中使用bcftools norm拆分 *
        '''
        info_list = info_str.split(":")
        # varscan format
        if format == "GT:GQ:DP:RD:AD:FREQ:DP4":
            DP_index, AD_index, AF_index = 2, 4, 5
            DP, AD = int(info_list[DP_index]), int(info_list[AD_index])
            AF = float(info_list[AF_index].replace("%", ""))/100
            
        # mutect2 format
        elif format == "GT:AD:AF:DP:F1R2:F2R1" or format == "GT:AD:AF:DP:F1R2:F2R1:PGT:PID:PS" \
            or format == 'GT:AD:AF:DP:F1R2:F2R1:FAD:SB' or format == 'GT:AD:AF:DP:F1R2:F2R1:FAD:PGT:PID:PS:SB':
            DP_index, AD_index, AF_index = 3, 1, 2
            AD = int(info_list[AD_index].split(",")[1])
            DP, AF = int(info_list[DP_index]), float(info_list[AF_index])
            
        # muse format
        elif format == "GT:DP:AD:BQ:SS":
            DP_index, AD_index = 1, 2
            AD = int(info_list[AD_index].split(",")[1])
            DP = int(info_list[DP_index])
            AF = AD / DP if DP > 0 else 0

        # pindel format
        elif format == "GT:AD":
            AD_index = 1
            D1, D2 = info_list[AD_index].split(",")
            D1, D2 = int(D1), int(D2) 
            DP = D1 + D2
            AD = D2
            AF = AD / DP if DP > 0 else 0

        else:
            raise ValueError(f"{format} is an Unexpected format")
        
        return DP, AD, AF

    def filter(self, DP: int = 30, AD: int = 3, AF: float = 0.05):
        filter_mask = (self['DP'] >= DP) & (self['DP'] >= AD) & (self['AF'] >= AF)
        return VCF(self[filter_mask])
    
    def sort(self):
        '''Sort by CHROM and POS'''
        self.sort_values(by=['CHROM', 'POS'], inplace=True)
        self.reset_index(drop=True, inplace=True)

    def load_from_variant_index(self, variant_index_df):
        split_info = variant_index_df
        split_info["index"] = variant_index_df.index
        split_info[['CHROM', 'POS', 'REF', 'ALT']] = split_info['index'].str.split('|', expand=True)
        split_info = split_info.drop(columns='index')
        split_info = split_info[['CHROM', 'POS', 'REF', 'ALT'] + [col for col in split_info.columns if col not in ['CHROM', 'POS', 'REF', 'ALT']]]
        self._update_inplace(split_info)

class CaseVCFs:
    def __init__(self, case_ID):
        self.case_ID = case_ID
        self.VCFs = {}
        self.filtered_VCFs = {}
        self.caller_list = ["muse", "mutect2", "varscan", "pindel"]
        self.caller_stat_df = pd.DataFrame(
            columns=["all"] 
                    + [f"{caller}_raw" for caller in self.caller_list]
                    + [f"{caller}_filtered" for caller in self.caller_list]
                    + ["2_caller"],
            index=pd.Index([self.case_ID], name="case_ID")  # 將self.case_ID包裝在列表中
            )
        
    def load_VCFs(self, variants_calling_dir=None, muse_suffix=None, mutect_suffix=None, varscan_indel_suffix=None, varscan_snp_suffix=None, pindel_suffix=None):
        # set 
        variants_calling_dir = variants_calling_dir or "variants_calling"
        muse_suffix = muse_suffix or ".MuSE.vcf"
        mutect_suffix = mutect_suffix or ".normed.vcf"
        varscan_indel_suffix = varscan_indel_suffix or ".varscan.snp.Somatic.hc.vcf"
        varscan_snp_suffix = varscan_snp_suffix or ".varscan.indel.Somatic.hc.vcf"
        pindel_suffix = pindel_suffix or ".filtered.vcf"

        # load VCFs
        muse_vcf = VCF(f"{variants_calling_dir}/muse/{self.case_ID}{muse_suffix}")
        mutect_vcf = VCF(f"{variants_calling_dir}/mutect2/{self.case_ID}{mutect_suffix}")
        varscan_snp_vcf = VCF(f"{variants_calling_dir}/varscan/{self.case_ID}{varscan_snp_suffix}")
        varscan_indel_vcf =  VCF(f"{variants_calling_dir}/varscan/{self.case_ID}{varscan_indel_suffix}")
        varscan_vcf = varscan_snp_vcf.concat(varscan_indel_vcf)
        pindel_vcf = VCF(f"{variants_calling_dir}/pindel/{self.case_ID}{pindel_suffix}")

        self.VCFs = {
            "muse": muse_vcf,
            "mutect2": mutect_vcf,
            "varscan": varscan_vcf,
            "pindel": pindel_vcf
        }

    
    def filter_VCFs(self):
        hit_dfs = []
        for package_name, vcf in self.VCFs.items():
            filtered_vcf = vcf.filter()
            filtered_vcf_df = filtered_vcf.to_variant_index_df()
            hit_dfs.append(pd.DataFrame(True, index=filtered_vcf_df.index, columns=[package_name])) # true value mutation dataframes
            self.filtered_VCFs[package_name] = filtered_vcf
        self.is_hit_df = pd.concat(hit_dfs, axis=1).fillna(False)
        return self.is_hit_df
    
    def caculate_caller_info(self):
        self.caller_stat_df.loc[self.case_ID, "all"] = len(self.is_hit_df)
        for package_name in self.VCFs.keys():
            vcf = self.VCFs[package_name]
            filtered_vcf = self.filtered_VCFs[package_name]
            self.caller_stat_df.loc[self.case_ID, [f"{package_name}_raw", f"{package_name}_filtered"]] = len(vcf), len(filtered_vcf)
        caller_2hit_filter = self.is_hit_df.loc[:, self.caller_list].sum(axis=1) >= 2
        self.caller_stat_df.loc[self.case_ID, "2_caller"] = caller_2hit_filter.sum()
        return self.caller_stat_df
    

def split_index_to_columns(df: pd.DataFrame):
    split_info = df.copy()
    split_info["tmp"] = df.index
    split_info[['CHROM', 'POS', 'REF', 'ALT']] = split_info["tmp"].str.split('|', expand=True)
    split_info = split_info.drop(columns='tmp')
    split_info = split_info[['CHROM', 'POS', 'REF', 'ALT'] + [col for col in split_info.columns if col not in ['CHROM', 'POS', 'REF', 'ALT']]]
    split_info.reset_index(drop=True)
    return split_info

def caller_venn_plot(all_case_hit_df: pd.DataFrame, caller_list):
    labels = venn.get_labels(all_case_hit_df.loc[:, caller_list])
    fig, ax = venn.venn4(labels, names=caller_list)
    fig.savefig("variants_calling/venn.png")

def parse_arguments():
    parser = argparse.ArgumentParser(description="Process VCF files.")
    parser.add_argument("project_dir", type=str, help="The project directory")
    parser.add_argument("--muse_suffix", type=str, default=".MuSE.vcf", help="Suffix for MuSE VCF files")
    parser.add_argument("--mutect_suffix", type=str, default=".normed.vcf", help="Suffix for Mutect VCF files")
    parser.add_argument("--varscan_indel_suffix", type=str, default=".varscan.snp.Somatic.hc.vcf", help="Suffix for VarScan indel VCF files")
    parser.add_argument("--varscan_snp_suffix", type=str, default=".varscan.indel.Somatic.hc.vcf", help="Suffix for VarScan SNP VCF files")
    parser.add_argument("--pindel_suffix", type=str, default=".filtered.vcf", help="Suffix for Pindel VCF files")
    
    return parser.parse_args()
def main(args):
    project_dir = args.project_dir
    os.chdir(project_dir)
    os.makedirs("variants_calling/merged_vcf", exist_ok=True)

    print("Reading case IDs from 'case_IDs.tsv'...")
    case_IDs = pd.read_csv("case_IDs.tsv", sep="\t", header=None)[0].tolist()
    print(f"Found {len(case_IDs)} case IDs.")

    all_case_hit_df = pd.DataFrame()
    all_caller_info = pd.DataFrame()

    for case_ID in case_IDs:
        print(f"\nProcessing case ID: {case_ID}")
        case_vcfs = CaseVCFs(case_ID)

        # Load VCFs with specified suffixes
        print("Loading VCFs...")
        case_vcfs.load_VCFs(
            muse_suffix=args.muse_suffix,
            mutect_suffix=args.mutect_suffix,
            varscan_indel_suffix=args.varscan_indel_suffix,
            varscan_snp_suffix=args.varscan_snp_suffix,
            pindel_suffix=args.pindel_suffix
        )

        # Filter VCFs and calculate caller info
        print("Filtering VCFs...")
        is_hit_df = case_vcfs.filter_VCFs()
        caller_info = case_vcfs.caculate_caller_info()
        all_caller_info = pd.concat([all_caller_info, caller_info])
        
        is_hit_df["case"] = case_ID
        all_case_hit_df = pd.concat([all_case_hit_df, is_hit_df])

    caller_list = ["muse", "mutect2", "varscan", "pindel"]
    print("\nIdentifying cases with hits from at least two callers...")
    all_case_2_hit_df = all_case_hit_df[all_case_hit_df.loc[:, caller_list].sum(axis=1) >= 2]
    print(f"Found {len(all_case_2_hit_df)} cases with hits from at least two callers.")

    all_case_2_hit_full_df = split_index_to_columns(all_case_2_hit_df)

    # Save to CSV files
    print("Saving results to CSV files...")
    all_case_2_hit_full_df.to_csv("variants_calling/all_case_2_hit_full.csv", index=False)
    all_caller_info.to_csv("variants_calling/all_caller_info.csv", index=False)

    # Optionally, call the Venn plot function here
    print("Generating Venn plot...")
    caller_venn_plot(all_case_hit_df, caller_list)
    print("Process completed successfully.")


if __name__ == "__main__":
    args = parse_arguments()
    main(args)
