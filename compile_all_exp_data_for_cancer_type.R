#########################
# this script finds the TCGA expression data of a certain type (scaled estimate or normalized)
# for all samples of a certain cancer type and compiles them into one table with TCGA barcodes
# and expression data for each gene per each sample

# EXAMPLE USE:
# source('/path/of/this/script/correlate_unc_id_with_tcga_bc.R')
# compile_sample_data('luad', '/path/to/out/file/all_tcga_luad_scaled_estimate_exp.tsv')
# this will take a few minutes to run. 
# it will output a table like so:
# gene_id   TCGA.23.2304.123   TCGA.23.2304.123   TCGA.23.2304.112        ....
# A1BG      2.59e-5             5.125e-5            4.70e-5               ....
# A1CG      4.70e-5             2.59e-5             5.125e-5              ....
# A2BP1     0                   4.70e-5             0                     ....
# ....      ....                ....                ....                  ....
#########################

read.tsv = function(file_loc){
  return(read.table(file_loc, sep='\t', header = T, stringsAsFactors = F))
}

extract_unc_id_from_filename = function(filename){
  library(stringr)
  #regex matches pattern of unc_barcode xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx
  f = str_extract(filename, "[a-z0-9]{8}-[a-z0-9]{4}-[a-z0-9]{4}-[a-z0-9]{4}-[a-z0-9]{12}")
  return(f)
}

find_tcga_bc = function(mage_tab, unc_id){
  tcga_bc_long = unique(mage_tab[mage_tab$Extract.Name == unc_id, 'Comment..TCGA.Barcode.'])
  if(length(tcga_bc_long) == 0){ return('x')}
  tcga_bc_short = paste(strsplit(tcga_bc_long, '-')[[1]][1:4], collapse='.')
  return(tcga_bc_short)
}

add_table = function(compiled_table, table_filename, exp_col_name, mage_tab){
  #reads expression data from each file, adds it to compiled table. 
  df = read.tsv(table_filename)
  df = df[,c('gene_id', exp_col_name)]
  unc_id = extract_unc_id_from_filename(table_filename)
  tcga_bc = find_tcga_bc(mage_tab, as.character(unc_id))
  
  if(tcga_bc == 'x'){
    print(paste("UNC barcode", unc_id, "was not found in mage_tab file. It was left out of the table."))
    return(compiled_table)}
  colnames(df) = c('gene_id', tcga_bc)
  compiled_table = merge(compiled_table, df, by='gene_id')
  return(compiled_table)
}

file_loc_list = function(base_dir, cancer_type, file_type, out_loc){
  file_locs = list( base_dir = paste('/taylorlab/data/tcga/tcga_rna/cancer_type/', cancer_type, sep=''))
  file_locs$gene_exp_files = list.files(file_locs$base_dir, pattern = file_type, full.names=T, recursive=T)
  file_locs$mage_tab_loc = list.files(file_locs$base_dir, pattern = 'sdrf.txt', full.names=T, recursive=T)
  file_locs$out_file = out_loc
  return(file_locs)
}

check_exp_type = function(file_type){
  if(file_type == 'rsem.genes.results'){ return('scaled_estimate') }
  else if(file_type == 'rsem.genes.normalized_results'){ return('normalized_count') }
  else{ stop('file type must be either "rsem.genes.results" or "rsem.genes.normalized_results".') }
}

make_starter_table = function(file_locs, mage_tab, exp_col_name){
 #make starter table to add all the rest onto
  all_gene_exp = read.tsv(file_locs$gene_exp_files[1])
  all_gene_exp = all_gene_exp[,c('gene_id', exp_col_name)] 
  
  first_unc_id = extract_unc_id_from_filename(file_locs$gene_exp_files[1])
  first_tcga_bc = find_tcga_bc(mage_tab, first_unc_id)
  colnames(all_gene_exp) = c('gene_id', first_tcga_bc)
  return(all_gene_exp)
}

remove_num_from_gene_id = function(gene_id_col){
  # removes number ID from gene_id column (they are originally like ABC123|20394 )
  return( unlist(sapply(gene_id_col, function(x){
    return(strsplit(x, split='|', fixed=T)[[1]][1]) }))
  )
}

compile_sample_data = function(cancer_type, out_loc, file_type = 'rsem.genes.results', return_table=F){
  
  file_locs = file_loc_list(base_dir, cancer_type, file_type, out_loc)
  exp_col_name = check_exp_type(file_type)
  mage_tab = read.tsv(file_locs$mage_tab_loc)
  all_gene_exp = make_starter_table(file_locs, mage_tab, exp_col_name)
  # compiles table, adds one column for each sample
  for(f in file_locs$gene_exp_files[2:length(file_locs$gene_exp_files)]){
    all_gene_exp = add_table(all_gene_exp, f, exp_col_name, mage_tab)
  }
  all_gene_exp$gene_id = remove_num_from_gene_id(all_gene_exp$gene_id)
  all_gene_exp = all_gene_exp[!(all_gene_exp$gene_id=='?'),]
  write.table(all_gene_exp, file_locs$out_file , sep='\t', row.names=F, quote=F)

  if(return_table == T){
    return(all_gene_exp)
  }
}
