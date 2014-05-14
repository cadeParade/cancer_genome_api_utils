library(cgdsr)
library(stringr)
library(plyr)


#Wrappers for base functions
#Brief documentation for base functions here: http://www.cbioportal.org/public-portal/cgds_r.jsp

cbio_get_tcga_studies = function(con, cancer_type){
  studies= getCancerStudies(con) #all studies
  tcga_cancer_type = paste(cancer_type, 'tcga', sep='_')
  tcga_studies_for_cancer_type = studies[ str_detect(studies$cancer_study_id, tcga_cancer_type) , 'cancer_study_id']
  return(tcga_studies_for_cancer_type)
}

cbio_get_sample_set_id = function(con, study_id, case_type){
  cases = getCaseLists(con, study_id)
  single_case_id = cases[ str_detect(cases$case_list_id, case_type), 'case_list_id']
  return(single_case_id)
}

cbio_get_sample_sets = function(con, study_id){
  cases = getCaseLists(con, study_id)[1:4]
  return(cases)
}

cbio_get_profile_id = function(con, study_id, profile_type){
  profiles = getGeneticProfiles(con, study_id)
  single_profile_id = profiles[ str_detect(profiles$genetic_profile_id, profile_type), 'genetic_profile_id']
  return(single_profile_id)
}

# cbio_base_fetch_data = function(con, study_id, data_type, genes, case_type = 'tcga_all'){
cbio_base_fetch_data = function(con, study_id, data_type, genes, case_type = '3way_complete'){
  case_id = cbio_get_sample_set_id(con, study_id, case_type)
  profile_id = cbio_get_profile_id(con, study_id, data_type)
  d = getProfileData(con, genes, profile_id, case_id)[-1,]
  return(d)

}


# Wrappers for specific data types from cbio_portal
cbio_clinical = function(con, study_id, case_type = '3way_complete'){
  case_id = cbio_get_sample_set_id(con, study_id, case_type)
  clin = getClinicalData(con, case_id)
  return(clin)
}

cbio_rppa = function(con, study_id, abs){
  #TODO find some way to validate gene/antibody names.
  rppa_df = cbio_base_fetch_data(con, study_id, 'RPPA', abs)
  colnames(rppa_df) = paste(colnames(rppa_df), 'prot', sep='_')
  return(rppa_df)
}

cbio_muts = function(con, study_id, genes){
  mut_df = cbio_base_fetch_data(con, study_id, 'mutations', genes)
  colnames(mut_df) = paste(colnames(mut_df), 'mut', sep='_')
  return(mut_df)
}

cbio_mrna = function(con, study_id, genes){
  mrna_df = cbio_base_fetch_data(con, study_id, 'rna_seq_v2_mrna$', genes)
  colnames(mrna_df) = paste(colnames(mrna_df), 'rsem', sep='_')
  return(mrna_df)
}


# Compiling all datatypes
cbio_fetch_all = function(cancer_type, genes, ab_names){
  library(cgdsr)
  library(stringr)
  library(plyr)
  source('/taylorlab/scripts/R/tcga_utils.R')
  study_id_num = set_study_id_num(cancer_type)
  cbio_con =  make_cbio_connection()
  study_id = cbio_get_tcga_studies(cbio_con, cancer_type)
  study_id = multiple_study_ids_warning(cancer_type, study_id, study_id_num)
  possible_cases = cbio_get_sample_sets(cbio_con, study_id)

  rppa = add_tcga_id_col(cbio_rppa(cbio_con, study_id, ab_names))
  muts = add_tcga_id_col(cbio_muts(cbio_con, study_id, genes))
  mrna = add_tcga_id_col(cbio_mrna(cbio_con, study_id, genes))
  clin = add_tcga_id_col(cbio_clinical(cbio_con, study_id))
  check_clin(clin, study_id)

  a = join_all(list(rppa, muts, mrna, clin), by='tcga_id')
  a$tcga_id = format_tcga_barcodes(a$tcga_id)

  protein_antibody_name_warning()
  return(a)

}

# Helper functions
set_study_id_num = function(cancer_type){
  if(cancer_type == 'luad' | cancer_type == 'lusc' | cancer_type == 'coadread'){
    return(2)
  }
  else{return(1)}
}

make_cbio_connection = function(){
  return(CGDS("http://www.cbioportal.org/public-portal/"))
}

add_tcga_id_col = function(d){
  d$tcga_id = rownames(d)
  return(d)
}

multiple_study_ids_warning = function(cancer_type, study_id, study_id_num){
  if(length(study_id) > 1){
    print(paste('Warning: there is more than one tcga study for', cancer_type, '. [', paste(study_id, collapse=', '), ']. Number', study_id_num, 'was chosen.'))
  }
  return(study_id[study_id_num])
}

protein_antibody_name_warning = function(){
  warning('if your protein data is incomplete, check /tayorlab/data/tcga/RPPA/full_ab_list.txt for possible antibody/gene names.')  
}

check_clin = function(clin, study_id){
  if(length(clin) == 0){
    warning(paste('There was no clinical data from cbio for', study_id, '.'))
  }
}
