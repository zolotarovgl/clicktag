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


.ct_import_data = function(ctmat_fn,cdna_lib_fn,ctbcs_fn){
  
  message('1. Loading data')
  message(paste0('Loading CT counts from: ',ctmat_fn))
  mat_ct = .ct_load_ct(ctmat_fn)
  message(paste0('Loading cDNA counts from: ',cdna_lib_fn,'\n'))
  mat_cdna = .ct_load_cr(cdna_lib_fn)
  int = intersect(colnames(mat_ct),colnames(mat_cdna))
  message(sprintf('%s common barcodes.',length(int)))
  mat_ct = mat_ct[,int]
  mat_cdna = mat_cdna[,int]
  # reorder cDNA matrix by the library size:
  mat_cdna = mat_cdna[,order(Matrix::colSums(mat_cdna),decreasing = T)]
  mat_ct = mat_ct[,colnames(mat_cdna)]
  N_CT_UMI=sum(mat_ct)
  message(sprintf('Total number of CT UMIs: %s',N_CT_UMI))
  # Load barcode information:
  bct = .ct_load_bc(ctbcs_fn,verbose = T)
  bc2lib = setNames(bct$barcode_pair,bct$name)
  
  
  assign("mat_ct", mat_ct, envir = .GlobalEnv)
  assign("mat_cdna", mat_cdna, envir = .GlobalEnv)
  assign("bct", bct, envir = .GlobalEnv)
  assign("bc2lib", bc2lib, envir = .GlobalEnv)
}

################################################################################
# Cell classification 
################################################################################

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
    message('Simulation done!')
    
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
  return(list(status = status, pval_mat = pval_mat))
}

.ct_sum_by_sample = function(mat_ct,bct){
  sapply(split(bct$name,bct$barcode_pair),FUN = function(x){
    colSums(mat_ct[x,])
  })
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
  #message(sprintf('Non-empty droplets (cDNA): %s\n%s (%s%%) droplets assigned unambiguously\n%s (%s%%) droplets with ambient CT profile \n%s (%s%%) droplets with ambiguous CT profile \n%s%% Discard rate ( (doublets + unclassifiable) / all non-ambient droplets)\n%s%% Global doublet rate ( Identified doublets / all non-empty droplets)\n%s%% Clicktag doublet rate ( Identified doublets / all non-empty-CT profiles)',
  #                n_nonempty,n_assigned,assignment_rate,n_ambient,ambient_rate,n_ambig,ambig_rate,discard_rate,doublet_rate1,doublet_rate2))
  return(setNames(c(n_nonempty,n_assigned,assignment_rate,n_ambient,ambient_rate,n_ambig,ambig_rate,discard_rate,n_doublet,doublet_rate1,doublet_rate2),c(c("n_nonempty","n_assigned","assignment_rate","n_ambient","ambient_rate","n_ambig","ambig_rate","discard_rate",'n_doublet',"doublet_rate1","doublet_rate2"))))
}

# 19.06.2024 - clearer reporting 
.do_ct_report = function(status,labels,putative_cells_ids){
  o = list()
  o$n_cells = length(putative_cells_ids)

  o$n_singlet = sum(table(status[!status %in% unlist(labels)]))
  o$n_doublet = table(status)[labels[['label_doublet']]]
  
  o$n_nonambient = sum(table(status)[names(table(status)) != labels[['label_ambient']]])
  
  o$n_ambient = sum(status == labels[['label_ambient']])
  o$n_ambig = table(status)[labels[['label_ambiguous']]]

  o$n_stained = o$n_singlet + o$n_doublet
  n_cells = o$n_cells 
  n_singlet = o$n_singlet
  n_doublet = o$n_doublet
  n_ambient = o$n_ambient
  n_ambig = o$n_ambig
  n_stained = n_singlet + n_doublet 

  o$n_nonambient = n_cells - n_ambient
  o$ambig_rate = round(n_ambig/n_cells*100,2)
  o$ambient_rate = round(n_ambient/n_cells*100,2) # proportion of cells with "empty" clicktag profile
  o$singlet_rate = round(n_singlet/n_cells*100,2)
  o$stain_rate = round((n_singlet + n_doublet)/n_cells*100,2)
  o$discard_rate = round((n_doublet+n_ambig+n_ambient)/n_cells*100,2)
  o$doublet_rate_stained = round((n_doublet)/(n_stained)*100,2)
  o$doublet_rate_cells = round((n_doublet)/(n_cells)*100,2)
  
  # print stats 
  #library(dplyr)
  #print(tibble(name = names(o), value = unlist(o)))

  return(o)
}

.estimate_doublet_rate = function(n_cells){
  # Publicly available multiplet rates 
  # https://www.biostars.org/p/9467937/
  exp = data.frame(
    rd = c(0.4, 0.8, 1.6, 2.3, 3.1, 3.9, 4.6, 5.4, 6.1, 6.9, 7.6)/100,
    n_cells_loaded = c(800, 1600, 3200, 4800, 6400, 8000, 9600, 11200, 12800, 14400, 16000),
    n_cells_recovered = c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000)
  )
  colnames(exp) = c('rd','n_cells_loaded','n_cells_recovered')
  lm1 = lm(data = exp,rd ~ n_cells_recovered)
  return(predict(lm1,data.frame(n_cells_recovered = n_cells)))
}

.run_lda = function(mat_ct,status,bl_labels,out_fn,postprob_thr = 0.9, multiplet_rate = NULL, clicktag_efficiency = NULL){
  require(MASS)
  if(is.null(multiplet_rate)){
    message("Multiplet rate hasn't been specified! Estimating using recovered cell counts")
    n_cells = length(status)
    global_multiplet_rate = round(.estimate_doublet_rate(n_cells),2)
    message(sprintf('Estimated multiplet rate for %s cells: %s',n_cells, global_multiplet_rate))
    # account for stained and unstained doublets"
    # Aimed stained doublet rate: Rc^2 * Rd 
    stained_multiplet_rate = global_multiplet_rate * clicktag_efficiency^2
    message(sprintf('Estimated stained multiplet rate: %s', stained_multiplet_rate))
    # ultiplet_rate = stained_multiplet_rate
    # CAVE: force more multiplets than actually expected!
    #multiplet_rate = global_multiplet_rate
    multiplet_rate = stained_multiplet_rate
  }else{
    message(sprintf('Multiplet rate specified: %s',round(multiplet_rate,2)))
  }
   # labels to exclude from training 
  sample_labels = levels(status)[!levels(status) %in% bl_labels]
  
  
  # The training excludes the doublets 
  #############################################
  
  #stop("This feature is not yet implemented.")
  # LDA is trained on the cells labeled by simulation
  message("Training LDA classifier using positive examples")
  message(sprintf('Posterior probability threshold: %s',postprob_thr))
  message(sprintf('Assumed multiplet rate: %s',multiplet_rate))
  
  # use 60% of the cells for training 
  n_train = round(length(status[!status %in% bl_labels])*0.6)
  message(sprintf('%s train examples.',n_train))
  train_ids = names(sample(status[!status %in% bl_labels],n_train,replace = F))
  test_ids = names(status[!status %in% bl_labels])[!names(status[!status %in% bl_labels]) %in% train_ids]
  test_dbl_ids = names(status)[status == 'doublet']
  
  training = as.data.frame(t(mat_ct[,train_ids]))
  training = log(training+1)
  training$label = as.character(status[rownames(training)])
  
  test = as.data.frame(t(mat_ct[,test_ids]))
  test = log(test+1)
  test$label = as.character(status[rownames(test)])
  
  test_dbl = log(as.data.frame(t(mat_ct[,test_dbl_ids]))+1)
  
  linear = lda(label ~ ., data = training)
  pred = predict(linear, training)$class
  tab <- table(Predicted = pred, Actual = training$label)
  train_acc = round(sum(diag(tab))/sum(tab),2)
  
  pred = predict(linear, test)$class
  tab <- table(Predicted = pred, Actual = test$label)
  test_acc = round(sum(diag(tab))/sum(tab),2)
  
  #message(sprintf('Training accuracy: %s',train_acc))
  #message(sprintf('Test accuracy: %s',test_acc))
  
  # CAVE: the ground truth is unknown here 
  train_ids = names(status[!status %in% unlist(labels)])
  training = as.data.frame(t(mat_ct[,train_ids]))
  training = log(training+1)
  training$label = as.character(status[rownames(training)])
  linear = lda(label ~ ., data = training)
  
  fname = sprintf('%s/lda.RData',out_fn)
  saveRDS(linear,fname)
  sprintf('Model saved: %s',fname)
  linear = readRDS(fname)
  
  # ploat loadings 
  m = as.matrix(linear$scaling)
  rownames(m) = paste0(rownames(m),': ',setNames(bct$barcode_pair,bct$name)[rownames(m)])
  cw = ch = unit(10,'mm')
  hm_lda = Heatmap(m,cluster_rows = F,cluster_columns = F,border = T,width = ncol(m)*cw,height = nrow(m)*ch)
  
  # predict all non-ambient profiles:
  topredict_ids = putative_cells_ids[status[putative_cells_ids] != 'ambient']
  topredict = as.data.frame(log(t(mat_ct[,topredict_ids])+1))
  pred = predict(linear, topredict)
  predicted = setNames(pred$class,rownames(topredict))
  ##############################################################################
  # now, identify the doublets by selecting the non-top probability threshold
  # multiplet rate should be 12-20% 
  # can we "slide" the probability threshold to get the correct number of doublets? 
  top_prob = apply(pred$posterior,1,max)
  nontop_prob = 1-apply(pred$posterior,1,FUN = function(x) sort(x,decreasing = T)[1])
  thrs = seq(0.0001,0.1,0.0001)

  multiplet_rates = sapply(thrs,FUN = function(x) sum(nontop_prob>=x)/length(nontop_prob))
  nontop_thr = thrs[which.min(abs(multiplet_rates-multiplet_rate))] # select the threshold that gives the 
  message(sprintf('LDA: Selected non-top posterior prob threshold: %s',nontop_thr))
  
  # closest proportion of multiplets to the aimed proportion 
  # you can also classify the multiplets here.
  # we have to add this on top of the real cells!
  status2 = setNames(as.character(apply(pred$posterior,1,FUN = function(x) colnames(pred$posterior)[which.max(x)])),rownames(pred$posterior))
  status2[top_prob<postprob_thr] = labels['label_ambiguous']
  status2[nontop_prob >= nontop_thr] = 'doublet'
  status2 = c(status2,setNames(as.character(status[status %in% labels['label_ambient']]),names(status[status %in% labels['label_ambient']])))
  status2 = status2[names(status)]
  
  status_res = status2
  status_res[status == 'doublet' | status2 == 'doublet'] = 'doublet'
  status_res = factor(status_res,levels = c(labels,sample_labels))
  # save posterior probabilities
  return(list(status = status_res, pred = pred))
}



.lda_sim_hm = function(status_input,status_lda){
  d = data.frame(input = as.character(status_input),by_lda = as.character(status_lda))
  m = dcast(d,input ~ by_lda,fill = 0)
  rownames(m) = m[,1]
  m = m[,-1]
  m = m[!rownames(m) %in% labels[labels != labels['label_doublet']],]
  m = m[,!colnames(m) %in% labels[labels != labels['label_doublet']]]
  
  mf = m
  m = m/sum(m)
  m[m==0] = NA
  
  ra = rowAnnotation(N = anno_barplot(rowSums(mf)))
  ha = HeatmapAnnotation(N = anno_barplot(colSums(mf)))
  hm = Heatmap(m,cluster_rows = F,cluster_columns = F,
               right_annotation = ra,
               bottom_annotation = ha,
               cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                 if(mf[i, j]>0){
                   grid.text(mf[i, j], x, y)
                 }
               },
               column_title = 'LDA classification',
               row_title = 'Simulation'
  )
  return(hm)
}
