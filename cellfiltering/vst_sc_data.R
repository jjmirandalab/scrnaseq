# Author: Michael Finlayson
#
# Description:
# This object will create a list of umi_counts and metadata for each library 
# that can be found within or below the current directory or directory 
# provided
#

library(Matrix)
library(sctransform)

######## -----------------------------------------------------------------------
########                 Parameters                                     
######## ---------------------------------------------------------------

VSTDParameters = setClass('VSTDParameters', 
  slots = c(
    batch_var = 'character',
    latent_var_nonreg = 'character'
  )
)

setMethod('initialize', 'VSTDParameters', 
function(.Object, batch_var = NA_character_, latent_var_nonreg = NA_character_) {
  .Object@batch_var = batch_var
  .Object@latent_var_nonreg = latent_var_nonreg
  return(.Object)
})


######## -----------------------------------------------------------------------
########                 Definition                                     
######## ---------------------------------------------------------------

library(Matrix)
library(Seurat)

VSTData = setClass('VSTData', 
  slots=c(
    parameters = 'VSTDParameters',
    corrected_umi_counts = 'Matrix', 
    logp1_umi_counts = 'Matrix',
    pearson_residuals = 'matrix',
    scaled_data = 'matrix',
    gene_stats = 'data.frame'
  )
)


######## -----------------------------------------------------------------------
########                 Initialization                                 
######## ---------------------------------------------------------------

setMethod('initialize', 'VSTData', function(.Object, umi_counts, metadata, parameters=VSTDParameters()) {
  
  stopifnot(all(colnames(umi_counts)==rownames(metadata)))
  
  num_cells = ncol(umi_counts)
  clip_magnitude_residuals = sqrt(num_cells)
  clip_magnitude_scale_data = num_cells^(1/3.7)
  
  batch_var = NULL
  if (!is.na(parameters@batch_var)) { batch_var = parameters@batch_var }
  latent_var_nonreg = NULL
  if (!is.na(parameters@latent_var_nonreg)) { latent_var_nonreg = parameters@batch_var }
  
  # call VST
  vst_out = vst(umi = umi_counts,
                cell_attr = metadata,
                latent_var = 'log_umi',
                batch_var = batch_var,
                latent_var_nonreg = latent_var_nonreg,
                n_genes = 3000,
                n_cells = NULL,
                method = 'poisson',
                do_regularize = T,
                res_clip_range = c(-clip_magnitude_residuals, clip_magnitude_residuals),
                bin_size = 256, 
                min_cells = 5, 
                residual_type = 'pearson',
                return_cell_attr = F, 
                return_gene_attr = T,
                return_corrected_umi = T, 
                bw_adjust = 3, 
                gmean_eps = 1,
                theta_given = NULL, 
                show_progress = TRUE)
  
  
  # scale the residuals
  pearson_residuals = vst_out$y
  
  scaled_data = pearson_residuals
  scaled_data[scaled_data < -clip_magnitude_scale_data] = -clip_magnitude_scale_data
  scaled_data[scaled_data > clip_magnitude_scale_data] = clip_magnitude_scale_data
  scaled_data = t(apply(scaled_data, 1, scale))
  colnames(scaled_data) = colnames(pearson_residuals)
  
  .Object@corrected_umi_counts = vst_out$umi_corrected
  .Object@logp1_umi_counts = log1p(vst_out$umi_corrected)
  .Object@pearson_residuals = pearson_residuals
  .Object@scaled_data = scaled_data
  .Object@gene_stats = vst_out$gene_attr
  
  return(.Object)
})

######## -----------------------------------------------------------------------
########                 Methods                                        
######## ---------------------------------------------------------------

###### ------------------------------
######     vstdata_is_null

setGeneric('vstdata_is_null', def=function(vstdata) 
{standardGeneric('vstdata_is_null')})

setMethod('vstdata_is_null', signature(vstdata='VSTData'), 
function(vstdata) {
  return(all(dim(vstdata@corrected_umi_counts)==c(1,1)) && all(is.na(vstdata@corrected_umi_counts)))
})

###### ------------------------------
######     vstdata_save

setGeneric('vstdata_save', def=function(vstdata, save_dir) 
{standardGeneric('vstdata_save')})

setMethod('vstdata_save', signature(vstdata='VSTData', save_dir='character'), 
function(vstdata, save_dir) {
  
  stopifnot(!vstdata_is_null(vstdata))
  stopifnot(dir.exists(save_dir))
  
  save_dir = paste(strsplit(save_dir, split = '/')[[1]], collapse = '/')
  
  # save parameters
  parameters_filename = paste(save_dir, 'vstdata_parameters.tsv', sep = '/')
  if (file.exists(parameters_filename)) { stop('VSTData parameters file exists already') }
  
  ptable = data.frame(parameter=c('batch_var', 'latent_var_nonreg'), value=rep(NA, times=2))
  for (p in ptable$parameter) {
    if (length(slot(vstdata@parameters, p) > 0)) { 
      ptable[p,'value'] = paste(slot(vstdata@parameters, p), collapse = ' ') 
    }
  }
  write.table(x = ptable,
              file = parameters_filename, 
              row.names = F, 
              sep = '\t')
  
  # save the corrected count matrix
  corrected_umi_counts_filename = paste(save_dir, 'corrected_umi_counts.csv.bz2', sep = '/')
  if (file.exists(corrected_umi_counts_filename)) { stop('corrected umi counts file exists already') }
  
  write.csv(vstdata@corrected_umi_counts, 
            file = bzfile(corrected_umi_counts_filename), 
            row.names = T)

  # save the pearson residuals
  pearson_residuals_filename = paste(save_dir, 'pearson_residuals.csv.bz2', sep = '/')
  if (file.exists(pearson_residuals_filename)) { stop('pearson residuals file exists already') }
  
  write.csv(vstdata@pearson_residuals, 
            file = bzfile(pearson_residuals_filename), 
            row.names = T)

  # save the scale data
  scaled_data_filename = paste(save_dir, 'scaled_data.csv.bz2', sep = '/')
  if (file.exists(scaled_data_filename)) { stop('scaled data file exists already') }
  
  write.csv(vstdata@scaled_data, 
            file = bzfile(scaled_data_filename), 
            row.names = T)
  
  # save the gene stats
  gene_stats_filename = paste(save_dir, 'vst_gene_stats.csv', sep = '/')
  if (file.exists(gene_stats_filename)) { stop('vst gene stats file exists already') }
  
  write.csv(vstdata@gene_stats, 
            file = gene_stats_filename, 
            row.names = T)
})



