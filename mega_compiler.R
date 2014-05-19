# Input a cancer_type and genes of interest into the function make_table to generate a combined table of 
# expression data and clincial information for that cancer type.

# Inputs: 
# REQUIRED
# cancer_type (eg. 'luad', 'skcm')
# genes (eg. c('BRAF', 'EGFR', 'PTEN'))
# OPTIONAL


# default = cols_to_keep = c('days_to_death', 'days_to_last_followup', 'days_to_new_tumor_event_after_initial_treatment', 'pathologic_stage', 'age_at_initial_pathologic_diagnosis', 'vital_status')) 

# output table 
# barcode   long_barcode    g1_exp      g2_exp     clin_col1_cleaned    clin_col2_cleaned
source('/taylorlab/scripts/R/cancer_genome_api_utils/compile_all_exp_data_for_cancer_type.R')
source('/taylorlab/scripts/R/cancer_genome_api_utils/tcga_utils.R')

load_clinical_table = function(cancer_type){
  clin_dir = paste('/taylorlab/data/tcga/tcga_rna/cancer_type/', cancer_type, sep='')
  clin_loc = list.files(clin_dir, pattern = 'clinical_patient*', full.names = T)
  clin = read.table(clin_loc[1], sep = '\t', header = T, stringsAsFactors = F, fill = T)  
  return(clin)
}

subset_and_clean_clinical_data = function(clin, cols_to_keep, cancer_type){
  clin = subset(clin, select = colnames(clin) %in% c('bcr_patient_barcode', cols_to_keep))
  clin$bcr_patient_barcode = format_tcga_barcodes(clin$bcr_patient_barcode)
  day_cols = cols_to_keep[str_detect(cols_to_keep, 'days')] #find columns counting days --need to be cleaned
  
  for(col_header in day_cols){
    clin[ , col_header ] = tcga_clean_days_to_cols(clin[ , col_header ])
  }

  if(cancer_type == 'prad'){
    clin$pathologic_stage_cleaned = tcga_pathologic_stage_cleaned(clin$pathologic_T)
  }
  else if(cancer_type == 'gbm' | cancer_type == 'lgg'){
    clin$pathologic_stage_cleaned = NA
  }
  else{
    clin$pathologic_stage_cleaned = tcga_pathologic_stage_cleaned(clin$pathologic_stage)
   } 

  clin$days_to_death_or_followup = tcga_combine_death_and_last_followup(clin$days_to_death, clin$days_to_last_followup) #TODO make ths more generic...
  id_col = which(colnames(clin) == 'bcr_patient_barcode')
  colnames(clin)[id_col] = 'tcga_id'
  return(clin)
}

numeric_and_tpm = function(df){
  for(i in 1:ncol(df)){
    df[, i] = as.numeric(as.character(df[ , i])) * 1e6
  }
  return(df)
}

subset_for_genes = function(all_exp_for_cancer, genes){
  gene_sub = as.data.frame(t(all_exp_for_cancer[ all_exp_for_cancer$gene_id %in% genes, ]))
  if(length(genes) == 1){
    gene_sub$temp = 0
  }
  return(gene_sub)
}

gene_sub_colnames = function(gene_sub){
  colnames(gene_sub) = as.character(unlist(gene_sub[1, ]))
  colnames(gene_sub) = paste(colnames(gene_sub), 'tpm', sep='_')
  #deleting first row because it is is characters
  return(gene_sub[-1, ])
}

make_table = function(cancer_type, genes, out_dir = getwd(), all_exp_file = NA, barcodes = 'all', clin_cols_to_keep = c('days_to_death', 'days_to_last_followup', 'pathologic_stage', 'age_at_initial_pathologic_diagnosis', 'vital_status')){

  #prad stage grading is different than ALL THE OTHERS. THANKS PRAD.
  if(cancer_type == 'prad'){
    clin_cols_to_keep  = str_replace(clin_cols_to_keep, 'pathologic_stage', 'pathologic_T')
  }
  else if(cancer_type == 'gbm' | cancer_type == 'lgg'){
    clin_cols_to_keep = c('days_to_death', 'days_to_last_followup', 'age_at_initial_pathologic_diagnosis', 'vital_status')
  }


  if(is.na(all_exp_file)){
    # compile all tables for cancer type
    all_exp_for_cancer = compile_sample_data(cancer_type, paste(out_dir, '/all_data_for_', cancer_type, '.txt', sep=''), return_table = T)  
  }
  else{
    all_exp_for_cancer = read.table(all_exp_file, sep='\t', header = T, stringsAsFactors = F)
  }

  gene_sub = subset_for_genes(all_exp_for_cancer, genes)
  gene_sub = gene_sub_colnames(gene_sub)
  gene_sub = numeric_and_tpm(gene_sub)

  gene_sub$tcga_id = tcga_format_barcodes(rownames(gene_sub))
  gene_sub$tcga_id_long = rownames(gene_sub)
  gene_sub = gene_sub[ !(duplicated(gene_sub$tcga_id)), ] # TODO make this dedup process more specific...
  
  if(length(genes) == 1){
    gene_sub[ , '0_tpm'] = NULL
  }

  sub_filename = paste(out_dir, '/', cancer_type, '_subset_for_', paste(genes, collapse='_'), '.txt', sep='')
  write.table(gene_sub, sub_filename, sep='\t', row.names = F, quote = F)

  clin = load_clinical_table(cancer_type)
  colnames(clin) = normalize_clin_headers(colnames(clin))
  clin = subset_and_clean_clinical_data(clin, clin_cols_to_keep, cancer_type)
  m = merge(gene_sub, clin, by='tcga_id')
  m_filename = paste(out_dir, '/', cancer_type, '_exp_and_clin_for_', paste(genes, collapse='_'), '.txt', sep='')
  write.table(m, m_filename, sep='\t', row.names = F, quote = F)
  return(m)
}
