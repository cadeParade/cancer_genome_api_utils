A brief description of each file:


cbio_utils.R: Wrappers for the cBioPortal.org R api. It will coallate RSEM (rna seq), RPPA (protein expression), mutation status, and clinical information (if available) for a list of genes and a cancer type. 

cbio_fetch_all(cancer_type, genes, ab_names)

example usage: cbio_fetch_all(skcm, c('PTEN', 'YAP1'), c('PTEN', 'YAP1_PS127'))

ARGUMENTS:

cancer_type = string abbreviation of cancer type (ex: 'skcm' or 'luad' or 'gbm'). This should also be the name of the directory where your rna_seq_v2 files are stored. See here for standard abbreviations: http://www.cbioportal.org/public-portal/webservice.do?cmd=getTypesOfCancer

genes = vector of strings of hugo gene names. If you're unsure of the correct name, try typing your gene name into the genes box on cbioportal.org. It will tell you if it's a recognized gene symbol or suggest alternatives. They are usally in all capital letters.

ab_names = vector of strings of protein antibody names. These are not always the same as the gene names above so if it doesn't work, check http://bit.ly/1mPBkWA for the correlations.

OUTPUT: a data frame with the following columns:
GENE1_prot  GENE_X_prot tcga_id  GENE1_mut  GENE_X_mut GENE1_rsem GENE_X_rsem clin_col_1  clin_col_2  clin_col_X


Caviats: There are some defaults set for (my personal) ease of use that might not be ideal for all cases.
1. It automatically picks which study to use. See set_study_id_num function.
2. It automatically chooses the 3_way_complete case type.



compile_all_exp_data_for_cancer_type.R: The TCGA firehose_get utility downloads several rna_seq_v2 expression files for each sample. It also stores ids of the samples as UNC barcodes. This file creates one table compiling all expression files for a given cancer type and associates the expression with a correlated TCGA barcode. 

Your rna seq files must be stored in this configuration: 
/base_dir/cancer_type/rna_seq_v2_files


compile_sample_data(cancer_type, out_loc, 
                    base_dir = '/taylorlab/data/tcga/tcga_rna/cancer_type/',  
                    file_type = 'rsem.genes.results', 
                    return_table = F)
example usage: compile_sample_data('skcm', '/path/to/out/dir')


ARGUMENTS: 
cancer_type = string abbreviation of cancer type (ex: 'skcm' or 'luad' or 'gbm'). This should also be the name of the directory where your rna_seq_v2 files are stored. See here for standard abbreviations: http://www.cbioportal.org/public-portal/webservice.do?cmd=getTypesOfCancer

out_location = string of path and filename of where you want the output table to go. 

base_dir = string of where your rna_seq_v2 files are stored without the cancer_type

file_type = can either be "rsem.genes.results" or "rsem.genes.normalized_results" depending on which type of data you want. It defaults to rsem.genes.normalized.

return_table = Boolean. Defaults to FALSE. If this is TRUE, it will return the compiled table as an R variable AND write the table to a file. If it is FALSE it will just write the table to a file.


OUTPUT:
Will always write a file to the location and name specified in out_loc argument.
If return_table is true it will return a data frame of this form: 

gene_id   TCGA.23.2304.123   TCGA.23.2304.123   TCGA.23.2304.112        ....
A1BG      2.59e-5             5.125e-5            4.70e-5               ....
A1CG      4.70e-5             2.59e-5             5.125e-5              ....
A2BP1     0                   4.70e-5             0                     ....
....      ....                ....                ....                  ....




MORE SOON....
