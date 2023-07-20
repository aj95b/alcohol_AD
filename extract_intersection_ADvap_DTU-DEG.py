#Script to find genes common between differential expression analysis of alcohol subjected mice with AD of P. Sanna et. al, 2023 and differential transcript usage (DTU) and expression analysis (meaning, differentially expressed genes - DEGs) of human AD patients by D. Coelho et. al, 2023

import pandas as pd
import sys

def read_data():
    ######### Parameter 1: #####################
    #Read in the statistically significant differentially expressed genes from AD vpor mice, pseudobulk
    #Stored at: 'AD_data/AD_vapor_2023_geo_sanna_lab/041_ADvap_SIGN.xlsx'
    ad_vap = pd.read_excel(sys.argv[1])
    #Change the gene name all capital letters for better comparison with the other dataset
    ad_vap['Gene'] = ad_vap['Gene'].str.upper()

    ######## Parameter 2: ######################
    #Read in the differentially expressed genes and genes from the DTU and DEG for AD human patients from Temporal Lob Intersection (TLI) by D. Cielho et. al
    #Stored at: 'AD_data/DTU_gene_expr_2021/41514_2020_52_MOESM3_ESM.xlsx'
    deg_tli = pd.read_excel(sys.argv[2],sheet_name=1) 
    dtu_tli = pd.read_excel(sys.argv[2],sheet_name=3)

    ######## Parameter 3: #####################
    #Read in the differential expression details matrix from D. Coelho et. al
    #Stored at: 'AD_data/DTU_gene_expr_2021/41514_2020_52_MOESM2_ESM.xlsx'
    deg_dtu_details = pd.read_excel(sys.argv[3],header=1)

    ######## Parameter 4: ####################
    #Path to output files
    output_path = sys.argv[4]
    
    return ad_vap, deg_tli, dtu_tli, deg_dtu_details, output_path

def get_common_genes_mice_human(ad_vap, deg_tli, dtu_tli, deg_dtu_details, output_path):
    #All four parameters in the same order as returned by the read data function

    #Compute the intersection of genes between AD vapor mice and human AD TLI
    intersect_advap_deg = set(ad_vap.Gene) & set(deg_tli.Genes) #24 in all
    intersect_advap_dtu = set(ad_vap.Gene) & set(dtu_tli.Genes) #28 in all

    #Extract expression details from the intersection of mice and DEG human
    x = deg_dtu_details.loc[deg_dtu_details.GENE_NAME.isin(intersect_advap_deg)]
    #Extract expression details from the intersection of mice and DTU human
    y = deg_dtu_details.loc[deg_dtu_details.GENE_NAME.isin(intersect_advap_dtu)]

    return x, y

def main():
    ad_vap, deg_tli, dtu_tli, deg_dtu_details, output_path = read_data()
    get_x, get_y = get_common_genes_mice_human(ad_vap, deg_tli, dtu_tli, deg_dtu_details, output_path)

    #Write the output frames to files
    get_x.to_csv(output_path+'diff_exp_details_comm_gene_ADvap_deg.tsv', sep='\t')
    get_y.to_csv(output_path+'diff_exp_details_comm_gene_ADvap_dtu.tsv', sep='\t')
  
if __name__ == "__main__":
    main()
