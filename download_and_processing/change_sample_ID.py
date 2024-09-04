import pandas as pd
import re
import os
import sys

# python change_sample_ID.py /home/blue/code/test_variant_calling/2sample_sheet.tsv /home/blue/code/test_variant_calling
def TCGA_sample_type(TCGA_barcode):
    pattern = r"-([0-2][0-9])[A-Z]$"
    match = re.search(pattern, TCGA_barcode)
    if not match:
        raise ValueError(f"Barcode {TCGA_barcode} does not match expected format.")
    sample_type = int(match.group(1))
    if sample_type < 10:
        return "T"  # Tumor
    elif sample_type < 20:
        return "N"  # Normal
    elif sample_type < 30:
        return "C"  # Control
    else:
        raise ValueError(f"Sample type code {sample_type} is not recognized.")

# 使用 sys.argv 來獲取命令行參數
sample_sheet_path = sys.argv[1]
output_dir = sys.argv[2]

# 如果 output_dir 不存在，則創建該目錄
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# 使用 os.path.splitext 來分離文件名和擴展名
file_name, file_ext = os.path.splitext(os.path.basename(sample_sheet_path))
new_sample_sheet_name = os.path.join(output_dir, f"{file_name}_new{file_ext}")

# 讀取 sample_sheet
sample_sheet = pd.read_csv(sample_sheet_path, sep="\t")

# 生成新的 Sample ID
sample_sheet["New sample ID"] = sample_sheet["Case ID"] + "-" + sample_sheet["Sample ID"].apply(TCGA_sample_type)

# 輸出新的sample_sheet 和 Case ID 和 New Sample ID
sample_sheet.to_csv(new_sample_sheet_name, index=False, sep="\t")
pd.Series(sample_sheet["Case ID"].unique()).to_csv(os.path.join(output_dir, f"{file_name}_case_IDs.tsv"), index=False, sep="\t", header=None)
sample_sheet["New sample ID"].to_csv(os.path.join(output_dir, f"{file_name}_sample_IDs.tsv"), index=False, sep="\t", header=None)
print("輸出新的sample_sheet, sample_IDs.tsv, case_IDs.tsv")