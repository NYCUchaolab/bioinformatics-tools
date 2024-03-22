# -*- coding: utf-8 -*-
"""
Created on Tue Dec  19 18:45:37 2023

@author: CJ (Thousand)
"""
import argparse
import pandas as pd


parser = argparse.ArgumentParser()
parser.add_argument('-v','--vcf', dest='input_vcf_file', help = '必要參數，例如: ASC_8samples_2callers_vep.vcf', required = True)
parser.add_argument('-m','--maf', dest='output_maf_file', help = '必要参数，例如: ASC_8samples_2callers_vep.maf', required = True)
parser.add_argument('-t','--tools', dest='output_maftools_maf_file', help = '必要参数，例如: ASC_8samples_2callers_vep_maftools.maf', required = True)
args = parser.parse_args()



def vep_vcf2maf(input_vcf_file, output_maf_file, output_maftools_maf_file):
    vcf_file = input_vcf_file
    
    # Open original vcf
    f = open(vcf_file)
    k = f.readlines()
    f.close()
    
    # Find the line include name
    name_line = None
    for i in range(len(k)):
        if '#CHROM' in k[i]:
            name_line = i
    
    # Open vcf as df
    data = pd.read_csv(vcf_file, sep='\t', header=name_line, dtype=object)
    
    # Find the line include INFO VEP annotation
    annotation_line = None
    for i in range(len(k)):
        if 'Consequence annotations from Ensembl VEP' in k[i]:
            annotation_line = i

    # Let INFO split to multiple columns
    data_FORMAT = data['FORMAT'].str.split('|', expand=True)
    data_FORMAT.rename(columns=dict(zip(data_FORMAT.columns, k[annotation_line].split('Format: ')[1].split('">\n')[0].split('|'))), inplace=True)
    
    data_maf = pd.concat([data.drop(columns='FORMAT'), data_FORMAT], axis=1)
    data_maf.to_csv(output_maf_file, index=False, sep='\t')
    
    # Make the maf for maftools used
    # Mandatory fields: Hugo_Symbol, Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2, 
    # Variant_Classification, Variant_Type and Tumor_Sample_Barcode.
    maf_done = data_maf[['SYMBOL', '#CHROM', 'POS', 'REF', 'ALT', 'VARIANT_CLASS', 'Consequence', 'sample']]
    maf_done = maf_done.rename(columns=dict(zip(maf_done.columns, ['Hugo_Symbol', 'Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Seq_Allele2', 
                                                                    'Variant_Type', 'Variant_Classification', 'Tumor_Sample_Barcode'])))
    # Start_Position, End_Position
    start_pos = []
    end_pos = []
    for i in range(len(maf_done)):
        if maf_done['Variant_Type'][i] == 'SNV':
            start_pos.append(int(maf_done['Start_Position'][i]))
            end_pos.append(int(maf_done['Start_Position'][i]))
        elif maf_done['Variant_Type'][i] == 'insertion':
            start_pos.append(int(maf_done['Start_Position'][i]))
            end_pos.append(int(maf_done['Start_Position'][i]) + 1)
        elif maf_done['Variant_Type'][i] == 'deletion':
            start_pos.append(int(maf_done['Start_Position'][i]) + 1)
            end_pos.append(int(maf_done['Start_Position'][i]) + len(maf_done['Reference_Allele'][i]))
    maf_done['Start_Change'] = start_pos
    maf_done['End_Position'] = end_pos
    maf_done = maf_done[['Hugo_Symbol', 'Chromosome', 'Start_Change', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2', 
                          'Variant_Classification', 'Variant_Type', 'Tumor_Sample_Barcode']]
    maf_done = maf_done.rename(columns={'Start_Change':'Start_Position'})
    
    # Variant_Type, Reference_Allele, Tumor_Seq_Allele2
    maf_done['Variant_Type'] = maf_done['Variant_Type'].replace({'SNV':'SNP', 'insertion':'INS', 'deletion':'DEL'})
    change_Reference_Allele = []
    for i in range(len(maf_done)):
        if maf_done['Variant_Type'][i] == 'INS':
            change_Reference_Allele.append('-')
        else:
            change_Reference_Allele.append(maf_done['Reference_Allele'][i])
    maf_done['Reference_Allele'] = change_Reference_Allele
    change_Tumor_Seq_Allele2 = []
    for i in range(len(maf_done)):
        if maf_done['Variant_Type'][i] == 'DEL':
            change_Tumor_Seq_Allele2.append('-')
        else:
            change_Tumor_Seq_Allele2.append(maf_done['Tumor_Seq_Allele2'][i])
    maf_done['Tumor_Seq_Allele2'] = change_Tumor_Seq_Allele2
    
    # Variant_Classification
    # vc_nonSyn: 'Frame_Shift_Del', 'Frame_Shift_Ins', 'Splice_Site', 'Translation_Start_Site', 'Nonsense_Mutation', 
    # 'Nonstop_Mutation', 'In_Frame_Del', 'In_Frame_Ins', 'Missense_Mutation'
    important_dict = {'transcript_ablation':1, 'exon_loss_variant':2, 'splice_donor_variant':3, 'splice_acceptor_variant':4, 'stop_gained':5, 
                      'frameshift_variant':6, 'stop_lost':7, 'start_lost':8, 'initiator_codon_variant':9, 'disruptive_inframe_insertion':10, 
                      'disruptive_inframe_deletion':11, 'conservative_inframe_insertion':12, 'conservative_inframe_deletion':13, 
                      'inframe_insertion':14, 'inframe_deletion':15, 'protein_altering_variant':16, 'missense_variant':17, 
                      'conservative_missense_variant':18, 'rare_amino_acid_variant':19, 'transcript_amplification':20, 'splice_region_variant':21, 
                      'splice_donor_5th_base_variant':22, 'splice_donor_region_variant':23, 'splice_polypyrimidine_tract_variant':24, 
                      'start_retained_variant':25, 'stop_retained_variant':26, 'synonymous_variant':27, 'incomplete_terminal_codon_variant':28, 
                      'coding_sequence_variant':29, 'mature_miRNA_variant':30, 'exon_variant':31, '5_prime_UTR_variant':32, 
                      '5_prime_UTR_premature_start_codon_gain_variant':33, '3_prime_UTR_variant':34, 'non_coding_exon_variant':35, 
                      'non_coding_transcript_exon_variant':36, 'non_coding_transcript_variant':37, 'nc_transcript_variant':38, 
                      'intron_variant':39, 'intragenic_variant':40, 'INTRAGENIC':41, 'NMD_transcript_variant':42, 'upstream_gene_variant':43, 
                      'downstream_gene_variant':44, 'TFBS_ablation':45, 'TFBS_amplification':46, 'TF_binding_site_variant':47, 
                      'regulatory_region_ablation':48, 'regulatory_region_amplification':49, 'regulatory_region_variant':50, 
                      'regulatory_region':51, 'feature_elongation':52, 'feature_truncation':53, 'intergenic_variant':54, 'intergenic_region':55}
    change_classification = []
    for i in range(len(maf_done)):
        if len(maf_done.loc[i, 'Variant_Classification'].split('&')) == 1:
            change_classification.append(maf_done.loc[i, 'Variant_Classification'].split('&')[0])
        elif len(maf_done.loc[i, 'Variant_Classification'].split('&')) != 1:
            largest = maf_done.loc[i, 'Variant_Classification'].split('&')[0]
            for j in range(len(maf_done.loc[i, 'Variant_Classification'].split('&'))):
                if important_dict[maf_done.loc[i, 'Variant_Classification'].split('&')[j]] < important_dict[largest]:
                    largest = maf_done.loc[i, 'Variant_Classification'].split('&')[j]
            change_classification.append(largest)
    maf_done['Variant_Classification'] = change_classification
    
    maf_done['Variant_Classification'] = maf_done['Variant_Classification'].replace({'splice_acceptor_variant':'Splice_Site', 'splice_donor_variant':'Splice_Site', 'transcript_ablation':'Splice_Site', 'exon_loss_variant':'Splice_Site', 
                                                                                     'stop_gained':'Nonsense_Mutation', 'stop_lost':'Nonstop_Mutation', 'splice_region_variant':'Splice_Region', 
                                                                                     'initiator_codon_variant':'Translation_Start_Site', 'start_lost':'Translation_Start_Site', 
                                                                                     'missense_variant':'Missense_Mutation', 'coding_sequence_variant':'Missense_Mutation', 'conservative_missense_variant':'Missense_Mutation', 'rare_amino_acid_variant':'Missense_Mutation', 
                                                                                     'transcript_amplification':'Intron', 'intron_variant':'Intron', 'INTRAGENIC':'Intron', 'intragenic_variant':'Intron', 
                                                                                     'incomplete_terminal_codon_variant':'Silent', 'synonymous_variant':'Silent', 'stop_retained_variant':'Silent', 'NMD_transcript_variant':'Silent', 
                                                                                     'mature_miRNA_variant':'RNA', 'exon_variant':'RNA', 'non_coding_exon_variant':'RNA', 'non_coding_transcript_exon_variant':'RNA', 'non_coding_transcript_variant':'RNA', 'nc_transcript_variant':'RNA', 
                                                                                     '5_prime_UTR_variant':"5'UTR", '5_prime_UTR_premature_start_codon_gain_variant':"5'UTR", '3_prime_UTR_variant':"3'UTR", 
                                                                                     'TF_binding_site_variant':'IGR', 'regulatory_region_variant':'IGR', 'regulatory_region':'IGR', 'intergenic_variant':'IGR', 'intergenic_region':'IGR', 
                                                                                     'upstream_gene_variant':"5'Flank", 'downstream_gene_variant':"3'Flank"})
    change_classification = []
    for i in range(len(maf_done)):
        if maf_done['Variant_Type'][i] == 'SNP':
            change_classification.append(maf_done['Variant_Classification'][i])
        elif maf_done['Variant_Type'][i] == 'DEL':
            if (maf_done['Variant_Classification'][i] == 'frameshift_variant') or (maf_done['Variant_Classification'][i] == 'protein_altering_variant'):
                change_classification.append('Frame_Shift_Del')
            elif maf_done['Variant_Classification'][i] == 'inframe_deletion':
                change_classification.append('In_Frame_Del')
            else:
                change_classification.append(maf_done['Variant_Classification'][i])
        elif maf_done['Variant_Type'][i] == 'INS':
            if (maf_done['Variant_Classification'][i] == 'frameshift_variant') or (maf_done['Variant_Classification'][i] == 'protein_altering_variant'):
                change_classification.append('Frame_Shift_Ins')
            elif maf_done['Variant_Classification'][i] == 'inframe_insertion':
                change_classification.append('In_Frame_Ins')
            else:
                change_classification.append(maf_done['Variant_Classification'][i])
    maf_done['Variant_Classification'] = change_classification

    change_classification = []
    for i in range(len(maf_done)):
        if maf_done['Variant_Classification'][i] in ['Splice_Site', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Splice_Region', 'Translation_Start_Site', 'Missense_Mutation', 'Intron', 
                                                     'Silent', 'RNA', "5'UTR", "3'UTR", 'IGR', "5'Flank", "3'Flank", 'Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins']:
            change_classification.append(maf_done['Variant_Classification'][i])
        else:
            change_classification.append('Targeted_Region')
    maf_done['Variant_Classification'] = change_classification
    
    maf_done.to_csv(output_maftools_maf_file, sep='\t', index=False)


if __name__ == '__main__':
    input_vcf_file = parser.parse_args().input_vcf_file
    output_maf_file = parser.parse_args().output_maf_file
    output_maftools_maf_file = parser.parse_args().output_maftools_maf_file
    vep_vcf2maf(input_vcf_file, output_maf_file, output_maftools_maf_file)





