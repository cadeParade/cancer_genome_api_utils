tcga_format_barcodes = function(barcodes){
  # replaces periods with dashes
  # removes all characters after the patient barcode
  library(stringr)
  barcodes = gsub('.', '-', barcodes, fixed=T)
  barcodes = str_sub(barcodes, 1, 12)
  return(barcodes)
}


tcga_combine_death_and_last_followup = function(d, l){
  # d = column listing days_to_death
  # l = column listing days_to_last_followup
  n = vector()
  for(i in 1:length(d)){
    #if both are na, return NA
    if(is.na(d[i]) & is.na(l[i])){
      n[i] = NA
    }
    #if only last followup is listed, return that one
    else if(is.na(d[i]) & !is.na(l[i])){
      n[i] = l[i]
    }
    # if only death is listed, return that one
    else if(!is.na(d[i]) & is.na(l[i])){
      n[i] = d[i]
    }
    # if both are listed
    else if(!is.na(d[i]) & !is.na(l[i])){
      #if they are both equal, return last followup
      if(d[i] == l[i]){
        n[i] = l[i]
      }
      #otherwise return days to death
      else{
        n[i] = d[i]
      }
    }
  }
  return(n)
}

# removes strings from days_to_* columns and makes it numeric
tcga_clean_days_to_cols = function(column){
  library(stringr)
  column = str_replace(column, '[Not Applicable]', NA)
  column = str_replace(column, '[Not Available]', NA)
  column = as.numeric(column)
  return(column)
}

#consolidateds clinical stage notations
tcga_pathologic_stage_cleaned = function(stage_col){
  cleaned = sapply(stage_col, function(x){
    if(x == 'Stage I' | x == 'Stage IA' | x == 'Stage IB'){
      return('1')
    }
    else if(x == 'Stage II' | x == 'Stage IIA' | x == 'Stage IIB' | x == 'Stage IIC' |
            x == 'T2a' | x == 'T2b' | x == 'T2c'){
      return('2')
    }
    else if(x == 'Stage III' | x == 'Stage IIIA' | x == 'Stage IIIB' | x == 'Stage IIIC' |
            x == 'T3a' | x == 'T3b'){
      return('3')
    }
    else if(x == 'Stage IV' | x == 'T4'){
      return('4')
    }
    else{
      return(NA)
    }
  })
  return(cleaned)
}

#for use in kaplan meier plots
high_or_low = function(cutoff, margin, vec){
  h_o_l = sapply(vec, function(x){
    if(is.na(x)){
      return(NA)
    }
    if(x >= cutoff + margin){
      return('high')
    }
    else if( x < cutoff - margin){
      return('low')
    }
    else{
      return(NA)
    }
  })
  return(h_o_l)
}
