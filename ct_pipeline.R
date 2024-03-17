# Clicktag pipeline with LDA classifier 
# 1. Classify using simulations 
# 2. Train the classifier 
# 3. Classify using LDA 
# 4. Cut off based on posterior probabilities


# Multiplet rate:
# https://kb.10xgenomics.com/hc/en-us/articles/360054599512-What-is-the-cell-multiplet-rate-when-using-the-3-CellPlex-Kit-for-Cell-Multiplexing- 
# ~0.8 per 1000 cells targeted = 12-20% per 15k target cells (usual setup by sebas)


# TODOS:
# rewrite the pipeline so it's built in a more logical way
# if useLDA, classify by simulation, then use this as a training set => run LDA
# output comparison between simulation and LDA
# Implement safety checks
################################################################################
#setwd('~/Documents/projects/mlei_development/ct_pipeline/')
#configfile = 'configs/sample.yaml'
################################################################################
getScriptPath <- function(){
	    cmd.args <- commandArgs()
    m <- regexpr("(?<=^--file=).+", cmd.args, perl=TRUE)
        script.dir <- dirname(regmatches(cmd.args, m))
        if(length(script.dir) == 0) stop("can't determine script dir: please call the script with Rscript")
	    if(length(script.dir) > 1) stop("can't determine script dir: more than one '--file' argument detected")
	    return(script.dir)
}
dir = getScriptPath()

args = commandArgs(trailingOnly=TRUE)
if(length(args)<1){
  message('Please, specify a config.yaml file!')
  quit()
}
options(warn=-1)
options(max.print = 20)
suppressMessages(library('ComplexHeatmap'))
suppressMessages(library("scales"))
suppressMessages(library("metacell"))
suppressMessages(library("yaml"))
suppressMessages(library("stringr"))
suppressMessages(library("reshape2"))
suppressMessages(library('circlize'))
source(sprintf('%s/helper.R',dir))
############ -1. Input & Output ###############
configfile = args[1]
#configfile = 'configs/231219_Mlei09_10x.ct_config.yaml'
message(sprintf('Reading pipeline setttings from %s.',configfile))
config = read_yaml(configfile)
config = unlist(config)
lib = config['lib']
# File paths:
ctmat_fn = config['ctmat_fn'] # path to matrix_counts/ from kallisto
ctbcs_fn = config['ctbcs_fn'] # path to barcodes.pairs.csv file 
cdna_lib_fn = config[['cdna_lib_fn']] # path to raw_feature_bc_matrix/ from cell ranger

out_fn = config['out_fn'] # name of the output directory 
figsdir = sprintf("%s/figs",out_fn)
scdb_path = config['scdb_path'] # a path to where put the metacell object to


########### 0. Pipeline settings ##############
cz_large_thr = as.integer(config[['cz_large_thr']])     # cells smaller or larger than this (in terms of expression) are considered non-cells
cz_small_thr = as.integer(config[['cz_small_thr']])
min_ct = as.integer(config[['min_ct']]) # the droplets with less CT UMIs than this value will be filtered out

# Simulation settings:
n_top_drops = as.integer(config[['n_top_drops']])
n_amb =as.integer( config[['n_amb']])
pval_thr = as.numeric(config[['pval_thr']]) # non-adjusted (yet) p-value threshold for a barcode in droplet. If droplet CT value (or sum, vis classification_mode) is lower 
# than this p-value, it gets assigned to the batch.
n_sim = as.integer(config[['n_sim']]) # number of ambient profile simulations per droplet 
classification_mode = config[['classification_mode']] # the mode by which to classify the cells:
# byCT - assigns the cell to a batch if a p-value is below `pval_thr` in at least `min_bc_matches`. 
#        each barcode is tested individually
# bySUM - instead of performing the tests for every barcode, the sum of the barcodes from each CT 
#         combination is tested.
# 3-4 * significance level (pval_thr)

# number of clicktags per sample: (only used for plotting)
num_ct = as.integer(config['num_ct'])
nth_ix = num_ct + 1
min_bc_matches = as.integer(config['min_bc_matches'])

# set up label names
labels = unlist(str_split(config[['labels']],','))
labels = setNames(labels,c('label_ambient','label_ambiguous','label_doublet'))
if(!dir.exists(out_fn)){
  dir.create(out_fn)
}
figsdir = paste0(out_fn,'/figs')
if(!dir.exists(figsdir)){
  dir.create(figsdir)
}

# check if LDA should be used 
if('useLDA' %in% names(config)){
  useLDA = as.logical(config['useLDA'])
  if(useLDA){
    message('Linear Discriminant Analysis has been specified.')
  }
  if(useLDA){
    #postprob_thr = 0.99
    postprob_thr = as.numeric(config['Posterior_prob_threshold'])
    multiplet_rate = as.numeric(config['Multiplet_rate'])
    #multiplet_rate = 0.16 # assumed multiplet rate
  }
}else{
  useLDA = FALSE
}

# init metacell database to export the filtered matrix
dir.create(scdb_path, showWarnings = TRUE,recursive = T)
metacell::scdb_init(scdb_path,force_reinit=TRUE)
########### 1. Load CT and cDNA data ##############
message('1. Loading the data')
mat_ct = .ct_load_ct(ctmat_fn)
mat_cdna = .ct_load_cr(cdna_lib_fn)
int = intersect(colnames(mat_ct),colnames(mat_cdna))
message(sprintf('%s common barcodes.',length(int)))
mat_ct = mat_ct[,int]
mat_cdna = mat_cdna[,int]
# reorder cDNA matrix by the library size:
mat_cdna = mat_cdna[,order(Matrix::colSums(mat_cdna),decreasing = T)]
mat_ct = mat_ct[,colnames(mat_cdna)]
# Load barcode information:
bct = .ct_load_bc(ctbcs_fn,verbose = T)
bc2lib = setNames(bct$barcode_pair,bct$name)
########### 2. Filter droplets based on cDNA and CT content ####################
message('2. Selecting the droplets by:\n - min cDNA UMI\n - min CT UMI')
cdna_size = Matrix::colSums(mat_cdna)
ct_size = Matrix::colSums(mat_ct)
putative_cells_ids = colnames(mat_cdna)[cdna_size >= cz_small_thr & cdna_size < cz_large_thr & ct_size >= min_ct]
message(sprintf('%s droplets passing the thresholds',length(putative_cells_ids)))
########### 3. Classify droplets by simulation #################################
message(paste0('3. Classifying the droplets by simulation ... mode:',classification_mode))
# This code should take as an input:
# - ambient_ids
# - mat_cdna
# - mat_ct 
# - bct 
# Output:
# `status` vector with cellID to barcode combination assignments 
################################################################################
# Select ambient droplets
message(sprintf('Selecting ambient droplets: lowest %s among top %s',n_amb,n_top_drops))
ambient_ids = tail(head(colnames(mat_cdna),n_top_drops),n_amb)
empty_set = mat_ct[,ambient_ids]
# get ambient proportions 
probs = rowSums(empty_set)/sum(empty_set)
# Simulate the cells:
# - this step is used to tell apart the empty cells vs non-empty cells 
status = .classify_cells_by_simulation(mat_ct[,putative_cells_ids],probs = probs,n_sim = n_sim,pval_thr = pval_thr,bc2lib = bc2lib,labels = labels,pval_mode = classification_mode,min_bc_matches = min_bc_matches)
# LDA should be better here 
if(useLDA){
  library(MASS)
  message("Training LDA classifier using positive examples")
  #postprob_thr = 0.99
  #multiplet_rate = 0.16 # assumed multiplet rate
  
  message(sprintf('Posterior probability threshold: %s',postprob_thr))
  
  bl_labels = unlist(labels) # labels to exclude from training 
  
  n_train = round(length(status[!status %in% bl_labels])*0.6)
  message(sprintf('%s train examples.',n_train))
  train_ids = names(sample(status[!status %in% bl_labels],n_train,replace = F))
  test_ids = names(status[!status %in% bl_labels])[!names(status[!status %in% bl_labels]) %in% train_ids]
  
  training = as.data.frame(t(mat_ct[,train_ids]))
  training = log(training+1)
  training$label = as.character(status[rownames(training)])
  
  test = as.data.frame(t(mat_ct[,test_ids]))
  test = log(test+1)
  test$label = as.character(status[rownames(test)])
  
  
  linear = lda(label ~ ., data = training)
  pred = predict(linear, training)$class
  tab <- table(Predicted = pred, Actual = training$label)
  train_acc = round(sum(diag(tab))/sum(tab),2)
  
  pred = predict(linear, test)$class
  tab <- table(Predicted = pred, Actual = test$label)
  test_acc = round(sum(diag(tab))/sum(tab),2)
  
  message(sprintf('Training accuracy: %s',train_acc))
  message(sprintf('Test accuracy: %s',test_acc))
  
  # CAVE: the ground truth is unknown here 
  train_ids = names(status[!status %in% unlist(labels)])
  training = as.data.frame(t(mat_ct[,train_ids]))
  training = log(training+1)
  training$label = as.character(status[rownames(training)])
  linear = lda(label ~ ., data = training)
  saveRDS(linear,sprintf('%s/lda.RData',out_fn))
  linear = readRDS(sprintf('%s/lda.RData',out_fn))
  
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
  
  # now, identify the doublets
  # multiplet rate should be 12-20% 
  top_prob = apply(pred$posterior,1,max)
  nontop_prob = 1-apply(pred$posterior,1,FUN = function(x) sort(x,decreasing = T)[1])
  thrs = seq(0.0001,0.1,0.0001)
  
  multiplet_rates = sapply(thrs,FUN = function(x) sum(nontop_prob>=x)/length(nontop_prob))
  nontop_thr = thrs[which.min(abs(multiplet_rates-multiplet_rate))] # select the threshold that gives the 
  message(sprintf('LDA: Selected non-top posterior prob threshold: %s',nontop_thr))
  # closest proportion of multiplets to the aimed proportion 
  # you can also classify the multiplets here.
  status2 = setNames(as.character(apply(pred$posterior,1,FUN = function(x) colnames(pred$posterior)[which.max(x)])),rownames(pred$posterior))
  status2[top_prob<postprob_thr] = labels['label_ambiguous']
  status2[nontop_prob >= nontop_thr] = 'doublet'
  status2 = c(status2,setNames(as.character(status[status %in% labels['label_ambient']]),names(status[status %in% labels['label_ambient']])))
  status2 = status2[names(status)]
  # now you can compare 2 classifications 
  d = data.frame(by_sim = as.character(status),by_lda = as.character(status2))
  m = dcast(d,by_sim ~ by_lda,fill = 0)
  rownames(m) = m[,1]
  m = m[,-1]
  m = m[!rownames(m) %in% labels[labels != labels['label_doublet']],]
  m = m[,!colnames(m) %in% labels[labels != labels['label_doublet']]]
  # how do we compare this properly? 
  # gosh, I want to smart as xavi
  # I'm repeating the same analyses over and over - why? 
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

  status_res = status2
  status_res[status == 'doublet' | status2 == 'doublet'] = 'doublet'
  sum(status_res == 'doublet')/length(status_res) # global status res
  sum(status_res == 'doublet')/length(status_res[status_res != 'ambient']) # stained doublet rate
  status = factor(status_res,levels = c(labels,bct$barcode_pair[!duplicated(bct$barcode_pair)]))

  # doublet counts heatmap 
  m = pred$posterior
  nontop_prob = 1-apply(m,1,FUN = function(x) sort(x,decreasing = T)[1])
  m = m[nontop_prob>=nontop_thr,]
  # now, how do we assign them? 
  # types of doublets:
  dbl_thr = 1-postprob_thr
  dbl = apply(m>=dbl_thr,1,FUN = function(x) paste0(sort(colnames(m)[which(x)]),collapse = '.'))
  dbl = dbl[str_count(dbl,'\\.')==1]
  d = as.data.frame(do.call(rbind,str_split(dbl,'\\.')))
  rownames(d) = names(dbl)
  d$V1 = factor(d$V1,levels = bct$barcode_pair[!duplicated(bct$barcode_pair)])
  d$V2 = factor(d$V2,levels = bct$barcode_pair[!duplicated(bct$barcode_pair)])
  m = dcast(d,V1 ~ V2,fill = 0,drop = F)
  rownames(m) = m[,1]
  m = m[,-1]
  m = m/sum(m)
  m[m==0] = NA
  hm_dbl = Heatmap(m,cluster_rows = F,cluster_columns = F,
          height = nrow(m)*ch,
          width = ncol(m)*cw,
          cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
            if(!is.na(m[i, j])){
              grid.text(round(m[i, j]*100,1), x, y)
            }
          },column_title = sprintf('Doublet proportions: %s droplets',nrow(d)))
  
  
  # Plot 
  fname = sprintf('%s/lda.pdf',figsdir)
  pdf(fname,height = 10, width = 10)
  draw(hm_lda)
  draw(hm)
  draw(hm_dbl)
  dev.off()  
  
  # save posterior probabilities
  write.table(pred$posterior,sprintf('%s/lda.posterior.tsv',out_fn),sep = '\t',quote = F)
  
  }
  message(cat(paste0(sapply(seq_along(levels(status)),FUN = function(i) paste0(names(table(status)[i]),' : ',table(status)[i])),collapse = '\n')))
stats = .do_ct_report(status,labels)
write.table(data.frame(stats),sprintf('%s/ct_stats.tsv',out_fn),quote = F,sep = '\t',col.names = F)

# the doublet rate is too low.
# how to become more restrictive 
write.table(status, sprintf('%s/cell_classification.tsv',out_fn),sep = '\t',quote=F,col.names = F)
write.table(t(mat_ct[,names(status)]),sprintf('%s/cell_ct.tsv',out_fn),sep = '\t',quote=F,col.names = T)
write.table(status, sprintf('%s/cell_classification.tsv',out_fn),sep = '\t',quote=F,col.names = F)
########## 4.5 Cell metadata ###################################################
# we need to separate the simulation 
md = data.frame(row.names = names(status),cdna_size = Matrix::colSums(mat_cdna)[names(status)],ct_size = Matrix::colSums(mat_ct)[names(status)])
md$ct_count_1st = apply(mat_ct[,names(status)],2,max)
md$ct_count_2nd = apply(mat_ct[,names(status)],2,FUN = function(x) sort(x,decreasing = T)[2])
md$ct_count_nth = apply(mat_ct[,names(status)],2,FUN = function(x) sort(x,decreasing = T)[nth_ix])
md$ct_name_1st = apply(mat_ct[,names(status)],2,FUN = function(x) rownames(mat_ct)[which.max(x)])
md$ct_name_2nd = apply(mat_ct[,names(status)],2,FUN = function(x) rownames(mat_ct)[which(order(x,decreasing = T)==2)])
md$ct_name_nth = apply(mat_ct[,names(status)],2,FUN = function(x) rownames(mat_ct)[which(order(x,decreasing = T)==nth_ix)])

mat_cdna_cells = mat_cdna[,rownames(md)]
mat_ct_cells = mat_ct[,rownames(md)]

md$ct_sum_topN = apply(mat_ct_cells,2,FUN = function(x) sum(sort(x,decreasing = T)[1:num_ct]))
log2_fist_nth = apply(mat_ct_cells,2,FUN = function(x) log2(sort(x,decreasing = T)[1]/(sort(x,decreasing = T)[nth_ix]+1)))
md$ct_log2_fist_nth = log2_fist_nth
md$ct_log2_sumtopN_other = log2(md$ct_sum_topN/(md$ct_size-md$ct_sum_topN))
md$clicktag_label = status[rownames(md)]
write.table(md, sprintf('%s/cell_ct_metadata.tsv',out_fn),sep = '\t',quote=F,col.names = T)

write.table(md, sprintf('%s/cell_ct_metadata.tsv',out_fn),sep = '\t',quote=F,col.names = T)
write.table(t(mat_ct[,rownames(md)]),sprintf('%s/cell_ct_counts.tsv',out_fn),sep = '\t',quote = F,col.names = T)
########### 4. Create metacell object ##########################################
message('4. Creating a metacell object')
#md = data.frame(row.names = names(status),cdna_size = Matrix::colSums(mat_cdna)[names(status)],ct_size = Matrix::colSums(mat_ct)[names(status)])
#md$clicktag_label = status[rownames(md)]
md = md[!md$clicktag_label %in% c(labels[['label_doublet']],labels[['label_ambiguous']],labels[['label_ambient']]),]
md$lib = lib
mat_f = tgScMat(mat_cdna[,rownames(md)], stat_type = "umi", cell_metadata = md)
# save filtered matrix
matid = sprintf("%s_CTfilt",lib)
metacell::scdb_add_mat(matid, mat_f)
message(sprintf("Added a new mat object to %s:\n%s/mat.%s.Rda",scdb_path,scdb_path,matid))

########### 5. QC and diagnostic plots #########################################
##### 1. Scatterplot of cDNA UMIs vs CT per CT ##########
par(cex.lab=0.5, cex.axis=0.5, cex.main=0.5, cex.sub=0.5)
fname = sprintf('%s/qc_cdna_ct_size.pdf',figsdir)
pdf(fname,height = 10, width = 10)
par(mfrow = c(3,2))
hist(log10(Matrix::colSums(mat_cdna)),breaks = 100,main = 'log10 cDNA UMIs [all droplets]',xlab = '')
abline(v = log10(cz_small_thr),lty = 2,col = 'red')
abline(v = log10(cz_large_thr),lty = 2,col = 'red')
hist(log10(Matrix::colSums(mat_ct)),breaks = 100,main = 'log10 CT UMIs [all droplets]',xlab = '')
abline(v = log10(min_ct),lty = 2,col = 'red')
par(mfrow = c(1,1))
plot(Matrix::colSums(mat_cdna),Matrix::colSums(mat_ct),pch = 16,cex = 0.3,col = scales::alpha('black',0.3),
     log = 'xy',xlab = 'cDNA UMIs',ylab = 'CT UMIs',main = 'cDNA ~ CT [all droplets]')
abline(v = cz_small_thr,lty = 2,col = 'red')
abline(v = cz_large_thr,lty = 2,col = 'red')
abline(h = min_ct,lty = 2,col = 'red')
par(mfrow = c(1,1))
plot(log10(Matrix::colSums(mat_cdna[,putative_cells_ids])),log10(Matrix::colSums(mat_ct[,putative_cells_ids])),
     pch = 16, col = scales::alpha('black',0.3),xlab = 'log10 cDNA UMIs', ylab = 'log10 CT UMIs',main = 'cDNA ~ CT [putative cells]')
lm1 = lm(log10(Matrix::colSums(mat_ct[,putative_cells_ids]))~log10(Matrix::colSums(mat_cdna[,putative_cells_ids])))
abline(lm1,col = 'red')

par(mfrow = c(4,2))
for(label in c(unlist(labels),bc2lib[!duplicated(bc2lib)])){
  if(sum(status == label)>10){
    plot(Matrix::colSums(mat_cdna[,names(status[status == label])]),Matrix::colSums(mat_ct[,names(status[status == label])]),
         pch = 16, col = scales::alpha('black',0.3),xlab = 'cDNA UMIs', ylab = 'CT UMIs',main = label,log = 'xy') 
  }
}
dev.off()

##### 2. CT counts heatmaps ##########
# - for top x cells including the empty 
# - for selected cells
fname = sprintf('%s/qc_ct_heatmap.pdf',figsdir)
pdf(fname,height = 15, width = 10)
par(mfrow = c(3,2))

ids_to_plot = c(putative_cells_ids,ambient_ids)
mr = as.matrix(t(mat_ct[,ids_to_plot]))
m_n = t(t(mr) / Matrix::colSums(mr, na.rm = TRUE)) * 1e4
m_n = apply(m_n, 2, function(c) c + quantile(c, 0.001, na.rm = TRUE))
m_n = m_n [ order(apply(m_n, 1, function(x) which.max(zoo::rollmean(x,2)) )), ]

mr = mr[rownames(m_n),]

# split the rows
rs = status[rownames(m_n)]
rs[is.na(rs)] = labels[['label_ambient']]
rs = factor(as.character(rs),levels = c(levels(rs)[-c(1:3)],levels(rs)[1:3]))
relev = setNames(paste0(names(table(rs)),': ',table(rs)),names(table(rs)))
rs = factor(relev[as.character(rs)],levels = relev[levels(rs)])

# barplot with counts 
barplot(table(rs),las=2)

n_assigned = sum(table(md$clicktag_label))
n_doublet = table(status)[labels[['label_doublet']]]
n_ambiguous = table(status)[labels[['label_ambiguous']]]
n_nonempty = sum(table(status)[names(table(status)) != labels[['label_ambient']]])
l = sprintf('%s droplets assigned unambiguously to one of the batches;\n%s%% doublet rate (doublets / all non-ambient doublets)',
            n_assigned,round((n_doublet+n_ambiguous)/n_nonempty*100,2))
mtext(side=3, line=2, at=-0.07, adj=0, cex=0.7, l)
# heatmap
hm1 = Heatmap(mr,show_row_dend = F,show_column_dend = F,show_row_names = F,
              cluster_columns = F,cluster_rows = F,col = xavis_green(mr),
              border = T,bottom_annotation = NULL,row_title_rot = 0,column_title = 'Raw',
              row_split = rs,row_title_gp = gpar(fontsize=5),show_heatmap_legend = F,use_raster = F)
hm2 = Heatmap(m_n,show_row_dend = F,show_column_dend = F,show_row_names = F,
              cluster_columns = F,cluster_rows = F,col = xavis_green(m_n),
              border = T,bottom_annotation = NULL,row_title_rot = 0,column_title = 'Norm',
              row_split = rs,row_title_gp = gpar(fontsize=5),show_heatmap_legend = F,use_raster = F)
draw(hm1+hm2)
dev.off()

##### 3. CT efficiency ##########
fname = sprintf('%s/qc_ct_efficiency.pdf',figsdir)
pdf(fname)
par(mfrow = c(2,3))
#par(mar = c(3,2,3,2))
par(cex.lab=0.5, cex.axis=0.5, cex.main=0.5, cex.sub=0.5,pch = 16)
# 1. cDNA cell size by clicktag label (cells only)
# 2. CT cell size by clicktag label (cells only)
# 3. CT 1st-nth FC per clicktag label (cells only).
#? 4. cDNA cell size per barcode (cells only)
#? 5. CT cell size per barcode (cells only)
#? 6. CT 1st-nth FC per barcode

mat_cdna_cells = mat_cdna[,rownames(md)]
mat_ct_cells = mat_ct[,rownames(md)]
# 1. cDNA cell size by clicktag label (cells only)
boxplot(split(Matrix::colSums(mat_cdna_cells),status),
        log = 'y',ylab = 'UMI/Cell',
        las=2,
        main = "cDNA cell size by clicktag label\n(cells only)")
# 2. CT cell size by clicktag label (cells only)
boxplot(split(Matrix::colSums(mat_ct_cells),status),
        log = 'y',ylab = 'UMI/Cell',
        las=2,
        main = " CT cell size by clicktag label\n(cells only)")
# 3. CT 1st-nth FC per clicktag label (cells only).
log2_fist_nth = apply(mat_ct_cells[,],2,FUN = function(x) log2(sort(x,decreasing = T)[1]/(sort(x,decreasing = T)[nth_ix]+1)))
boxplot(split(log2_fist_nth,status[names(log2_fist_nth)]),
        ylab = '1st-nth log2(FC)',
        las=2,
        main = "CT 1st-nth FC per clicktag label\n(cells only)")
# 4. cDNA cell size per barcode (cells only)
boxplot(split(Matrix::colSums(mat_cdna_cells),rownames(mat_ct_cells)[apply(mat_ct_cells,2,which.max)]),
        log = 'y',ylab = 'UMI/Cell',
        las=2,
        main = "cDNA cell size by clicktag barcode\n(cells only)")
# 5. CT cell size by clicktag barcode (cells only)
boxplot(split(Matrix::colSums(mat_ct_cells),rownames(mat_ct_cells)[apply(mat_ct_cells,2,which.max)]),
        log = 'y',ylab = 'UMI/Cell',
        las=2,
        main = "CT cell size by clicktag barcode\n(cells only)")
# 6. CT 1st-nth FC per barcode
boxplot(split(log2_fist_nth,rownames(mat_ct_cells)[apply(mat_ct_cells,2,which.max)]),
        ylab = '1st-nth log2(FC)',
        las=2,
        main = "CT 1st-nth FC per clicktag barcode\n(cells only)")

# ambient clicktag proportions
par(mfrow = c(2,3))
barplot(probs,main = sprintf('Background CT proportions estimated from %s droplets',n_amb),ylim = c(0,1),
        sub = 'Clicktag ambient proportions used as an input for classification.\n
        They are estimated by taking X lowest-size droplets (cDNA determined) and averaging across those.')
abline(h = 1/length(probs),lty = 2)
# compute the UMIfrac of the top clicktag combination 
prop_top_ct = setNames(sapply(1:ncol(mat_ct_cells),FUN = function(i){
  x = mat_ct_cells[,i]
  target_ct = names(bc2lib[bc2lib==status[colnames(mat_ct_cells)][i]])
  sum(x[target_ct])/sum(x)
}),colnames(mat_ct_cells))


# the ratio of the off-target clicktags for every cell
boxplot(split(prop_top_ct,status[names(prop_top_ct)]),
        ylab = '[top CT combionation]/[CT library size]',
        las=2,
        outline=FALSE,
        main = "Proportion of UMIs coming from a top combination\n(cells only)",
        sub = 'The proportion of UMIs coming from the clicktag combination the droplet has been assigned to.\nThe points and a line specify the rate of swappers (1-top_umifrac).')
swappers = 1-sapply(split(prop_top_ct,status[names(prop_top_ct)])[-c(1:3)],FUN = function(x) median(x,na.rm = T))
points(c(rep(NA,3),swappers),col = 'red')
abline(h = mean(swappers),col = 'red',lty = 2)
dev.off()
message("Clicktag pipeline done.")

