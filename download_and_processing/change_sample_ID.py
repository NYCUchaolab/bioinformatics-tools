import pandas as pd
import re
import os 

def TCGA_sample_type(TCGA_barcode):
    pattern = r"-([0-2][0-9])[A-Z]$"
    match = re.search(pattern, TCGA_barcode)
    sample_type = int(match.group(1))
    if sample_type < 10:
        return "T"
    elif sample_type < 20:
        return "N"
    elif sample_type < 30:
        return "C"
    raise ValueError("wrong sample type code.")

sample_sheet_dir = "/home/blue/code/TCGA_CHOL_2sample"
sample_sheet_name = "gdc_sample_sheet_test.2024-02-26.tsv"
sample_sheet_path = os.path.join(sample_sheet_dir, sample_sheet_name)

sample_sheet = pd.read_csv(sample_sheet_path, sep="\t")
sample_sheet["New sample ID"] = sample_sheet.loc[:, "Case ID"] + "-" + sample_sheet.loc[:, "Sample ID"].apply(TCGA_sample_type)
sample_sheet.to_csv(sample_sheet_path, index=False, sep="\t")

# sample case_IDS, sample_IDs
pd.Series(sample_sheet.loc[:, "Case ID"].unique()).to_csv("/home/blue/code/TCGA_CHOL_2sample/case_IDs.tsv", index=False, sep="\t", header=None)
sample_sheet.loc[:, "New sample ID"].to_csv("/home/blue/code/TCGA_CHOL_2sample/sample_IDs.tsv", index=False, sep="\t", header=None)