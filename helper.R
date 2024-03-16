# Data loading 
.ct_load_ct = function(ctmat_fn,verbose = F){
  # ctmat_fn - path to matrix_counts/
  require(stringr)
  mat_ct = t(Matrix::as.matrix(Matrix::readMM(file = sprintf("%s/genes.mtx", ctmat_fn))))
  rownames(mat_ct) = read.table(sprintf("%s/genes.genes.txt", ctmat_fn))[,1]
  colnames(mat_ct) = read.table(sprintf("%s/genes.barcodes.txt", ctmat_fn))[,1]
  if(verbose){
    message(sprintf('mat_ct: %s rows; %s columns',nrow(mat_ct),ncol(mat_ct)))
  }
  return(mat_ct)
}
.ct_load_cr = function(cdna_lib_fn){
  #
  #' `cdna_lib_fn` - path to the `raw_feature_bc_matrix/` of `cellranger`.
  #
  return(Seurat::Read10X(
    cdna_lib_fn,
    gene.column = 2,
    cell.column = 1,
    unique.features = TRUE,
    strip.suffix = T
  ))
}
.ct_load_bc = function(ctbcs_fn,verbose = F){
  require(stringr)
  #
  #' `ctbcs_fn` - path to the `_barcodes.pairs.csv` CSV file of the format:
  #' `BC name, barcode sequence, bc batch name`
  #
  bct = read.table(ctbcs_fn, sep = ",", col.names = c("name","barcode","barcode_pair"))
  if(verbose){
    message(sprintf('%s barcodes in %s combinations\n%s BCs per combination',length(unique(bct$name)) ,length(unique(bct$barcode_pair)),length(unique(bct$name))/length(unique(bct$barcode_pair)) ))
  return(bct)
  }
}

table_to_matrix = function(table) {
  mat = matrix(table, nrow = nrow(table))
  rownames(mat) = rownames(table)
  colnames(mat) = colnames(table)
  return(mat)
}

# Plotting 
xavis_green = function(m=NULL,qs = c(.25,.5,.75,.99),reverse = F){
  require("circlize")
  require("ComplexHeatmap")
  # This function mimics the color function from plot_complex_heatmap() function. 
  # It takes the matrix as an argument and outputs the color function to be used 
  # In a ComplexHeatmap call.
  pal = c("white","#d6e72e","#6fb600","#003f4d")
  if(reverse){
    pal = rev(pal)
  }
  return(colorRamp2(c(quantile(m, qs[1]),quantile(m, qs[3]),quantile(m, qs[3]),quantile(m, qs[4])),pal))
}


# Cell classification 
.classify_cells_by_simulation = function(mat_ct = NULL,probs = NULL,bc2lib = NULL,n_sim = 1000,pval_thr = 0.001,min_bc_matches = NULL,pval_mode = 'byCT',labels = labels){
  # labels - a named list containing the names for `label_ambient`, `label_ambiguous` and `label_doublet`
  message(sprintf('Simulating %s events for %s drops = %s M events => grab a coffee ...',n_sim,ncol(mat_ct),n_sim*ncol(mat_ct)/1e6))
  
  .pval_bysum = function(observed,probs,n_sim){
    # given "ambient" proportions, this function simulates an ambient clicktag profile 
    # of a  given depth
    #observed_p = observed/sum(observed)
    sim = rmultinom(n_sim,sum(observed),probs)
    # split by barcode combinations 
    return(sapply(split(names(bc2lib),bc2lib),FUN = function(y) sum(colSums(sim[y,])>=sum(observed[y]))/n_sim))
  }
  
  .pval_bycol = function(observed,probs,n_sim){
    # given "ambient" proportions, this function simulates an ambient clicktag profile 
    # of a  given depth
    observed_p = observed/sum(observed)
    sim = rmultinom(n_sim,sum(observed),probs)
    return(colSums(t(apply(apply(sim,2,FUN = function(z) z/sum(z)),2,FUN = function(x) x>observed_p)))/n_sim)
  }
  ##############################################################################
  
  if(pval_mode == 'byCT'){
    ######### P-value by individual barcodes 
    if(is.null(min_bc_matches)){
      message("Specify min_bc_matches!")
      break()
    }
    pval_mat = t(apply(mat_ct[,],2,FUN = function(x){
      .pval_bycol(x,probs,n_sim)
    }))
    message('Done!')
    
    # Assign the droplets to the libraries:
    pass_mat = sapply(split(names(bc2lib),bc2lib),FUN = function(x){
      rowSums(pval_mat[,x]<=pval_thr)
    })
    pass_mat[is.na(pass_mat)] = 0
    
    status = apply(pass_mat,1,FUN = function(x){
      if(length(x[x>0])>1&all(x[x>0]>=min_bc_matches)){
        return(labels[['label_doublet']])
      }else if(length(x[x>0])==1&all(x[x>0]>=min_bc_matches)){
        #return('singlet')
        return(colnames(pass_mat)[x>0])
      }else if(all(x==0)){
        return(labels[['label_ambient']])
      }else{
        # it matches against
        return(labels[['label_ambiguous']])
      }
    })
  }
  
  else if(pval_mode == "bySUM"){
    ######### P-value by the sum of barcodes 
    pval_mat = t(apply(mat_ct[,],2,FUN = function(x){
      .pval_bysum(x,probs,n_sim)
    }))
    message('Done!')
    
    pass_mat = pval_mat<=pval_thr
    
    # round the values
    pass_mat[pass_mat == 0] = 0
    pass_mat[pass_mat == 1] = 1
    pass_mat[is.na(pass_mat)] = 0
    status = apply(pass_mat,1,FUN = function(x){
      if(length(x[x>0])>1){
        return(labels[['label_doublet']])
      }else if(length(x[x>0])==1){
        #return('singlet')
        return(colnames(pass_mat)[x>0])
      }else if(all(x==0)){
        return(labels[['label_ambient']])
      }else{
        return(labels[['label_ambiguous']])
      }
    })
  }
  levs = c(c(labels[['label_ambient']],labels[['label_ambiguous']],labels[['label_doublet']]),bc2lib[!duplicated(bc2lib)])
  status = factor(status, levels = levs)
  return(status)
}

# Reporting 

.do_ct_report = function(status,labels){
  n_assigned = sum(table(status[!status %in% unlist(labels)]))
  n_doublet = table(status)[labels[['label_doublet']]]
  n_ambig = table(status)[labels[['label_ambiguous']]]
  n_nonambient = sum(table(status)[names(table(status)) != labels[['label_ambient']]])
  # whic rates should you report here?
  n_nonempty = length(putative_cells_ids)
  n_ambient = sum(status == labels[['label_ambient']])
  n_noambient = n_nonempty - n_ambient
  ambig_rate = round(n_ambig/n_nonempty*100,2)
  ambient_rate = round(n_ambient/n_nonempty*100,2) # proportion of cells with "empty" clicktag profile
  assignment_rate = round(n_assigned/n_nonempty*100,2)
  discard_rate = round((n_doublet+n_ambig+n_ambient)/n_nonempty*100,2)
  doublet_rate1 = round((n_doublet)/n_nonempty*100,2)
  doublet_rate2 = round((n_doublet)/n_noambient*100,2)
  message(sprintf('Non-empty droplets (cDNA): %s\n%s (%s%%) droplets assigned unambiguously\n%s (%s%%) droplets with ambient CT profile \n%s (%s%%) droplets with ambiguous CT profile \n%s%% Discard rate ( (doublets + unclassifiable) / all non-ambient droplets)\n%s%% Global doublet rate ( Identified doublets / all non-empty droplets)\n%s%% Clicktag doublet rate ( Identified doublets / all non-empty-CT profiles)',
                  n_nonempty,n_assigned,assignment_rate,n_ambient,ambient_rate,n_ambig,ambig_rate,discard_rate,doublet_rate1,doublet_rate2))
  return(setNames(c(n_nonempty,n_assigned,assignment_rate,n_ambient,ambient_rate,n_ambig,ambig_rate,discard_rate,n_doublet,doublet_rate1,doublet_rate2),c(c("n_nonempty","n_assigned","assignment_rate","n_ambient","ambient_rate","n_ambig","ambig_rate","discard_rate",'n_doublet',"doublet_rate1","doublet_rate2"))))
}

