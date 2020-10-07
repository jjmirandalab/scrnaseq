# Author: Michael Finlayson
#
# Description:
# This object will perform and encapsulate the standard process of preparing
# single-cell RNA-seq data for analysis
#


######## -----------------------------------------------------------------------
########                 Definition                                     
######## ---------------------------------------------------------------

library(Matrix)

SCDPrep = setClass('SCDPrep', 
  slots=c(
    scdata = 'SCData', 
    metadata = 'data.frame',
    umi_counts = 'Matrix',
    vstdata = 'VSTData',
    cell_filtration = 'CellFiltration',
    info_gene_selection = 'InfoGeneSelection'
  )
)


######## -----------------------------------------------------------------------
########                 Initialization                                 
######## ---------------------------------------------------------------

setMethod('initialize', 'SCDPrep', function(.Object, scdata) {
  stopifnot(class(scdata)=='SCData')
  .Object@scdata = scdata
  .Object@metadata = scdata@metadata
  .Object@umi_counts = scdata@umi_counts
  return(.Object)
})

######## -----------------------------------------------------------------------
########                 Methods                                        
######## ---------------------------------------------------------------

###### ------------------------------
######     scdprep_prep

setGeneric('scdprep_prep', def=function(scdprep) 
{standardGeneric('scdprep_prep')})

setMethod('scdprep_prep', signature(scdprep='SCDPrep'), 
function(scdprep) {
  
  scdprep = scdprep_filter(scdprep)
  scdprep = scdprep_vst(scdprep)
  scdprep = scdprep_select_info_genes(scdprep)
  
  return(scdprep)
})

###### ------------------------------
######     scdprep_filter

setGeneric('scdprep_filter', def=function(scdprep, ...) 
  {standardGeneric('scdprep_filter')})

setMethod('scdprep_filter', signature(scdprep='SCDPrep'), 
function(scdprep, parameters = CFParameters()) {
  
  scdprep@cell_filtration = CellFiltration(scdprep@scdata, parameters = parameters)
  
  scdata_filtered = scdata_limit_to_cells(scdata = scdprep@scdata,
                                          cells_to_keep = scdprep@cell_filtration@cells_to_keep)
  
  scdprep@metadata = scdata_filtered@metadata
  scdprep@umi_counts = scdata_filtered@umi_counts
  
  return(scdprep)
})

###### ------------------------------
######     scdprep_unfilter

setGeneric('scdprep_unfilter', def=function(scdprep, ...) 
{standardGeneric('scdprep_unfilter')})

setMethod('scdprep_unfilter', signature(scdprep='SCDPrep'), 
function(scdprep, parameters = CFParameters()) {
  
  scdprep@metadata = scdprep@scdata@metadata
  scdprep@umi_counts = scdprep@scdata@umi_counts
  scdprep@cell_filtration = cf_nullify(scdprep@cell_filtration)
  
  return(scdprep)
})

###### ------------------------------
######     scdprep_vst

setGeneric('scdprep_vst', def=function(scdprep, ...) 
{standardGeneric('scdprep_vst')})

setMethod('scdprep_vst', signature(scdprep='SCDPrep'), 
function(scdprep, parameters=VSTDParameters()) {
  
  scdprep@vstdata = VSTData(umi_counts=scdprep@umi_counts, 
                            metadata=scdprep@metadata, 
                            parameters=parameters)
  
  return(scdprep)
})

###### ---------------------------------
######     scdprep_select_info_genes

setGeneric('scdprep_select_info_genes', def=function(scdprep, ...) 
  {standardGeneric('scdprep_select_info_genes')})

setMethod('scdprep_select_info_genes', signature(scdprep='SCDPrep'), 
function(scdprep, parameters = IGSParameters()) {
  
  scdprep@info_gene_selection = InfoGeneSelection(scdprep@vstdata@corrected_umi_counts, 
                                                  selection_sets = subset(scdprep@metadata, select = 'sample'), 
                                                  parameters = parameters)
  
  return(scdprep)
})

###### --------------------
######     scdprep_save

setGeneric('scdprep_save', def=function(scdprep, save_dir, ...) 
{standardGeneric('scdprep_save')})

setMethod('scdprep_save', signature(scdprep='SCDPrep', save_dir='character'), 
function(scdprep, save_dir, dev_to_use='pdf') {
  
  stopifnot(dir.exists(save_dir))
  save_dir = paste(strsplit(save_dir, '/')[[1]], collapse = '/')
  
  data_dir = paste(save_dir, 'data', sep = '/')
  dir.create(data_dir)
  
  # save the scdata
  scdata_save(scdprep@scdata, data_dir)
  
  # save the vst data
  if (!vstdata_is_null(scdprep@vstdata)) { 
    vstdata_save(scdprep@vstdata, data_dir)
  }
  
  # save the cell filtration
  if (length(scdprep@cell_filtration@filtration_by_sample)!=0) { 
    cf_save_dir = paste(save_dir, 'cell_filtration', sep = '/')
    cf_save(scdprep@cell_filtration, cf_save_dir) 
  }
  
  # save the info gene selection
  if (length(scdprep@info_gene_selection@selection_data)!=0) { 
    igs_save_dir = paste(save_dir, 'info_gene_selection', sep = '/')
    save_info_gene_selection(scdprep@info_gene_selection, igs_save_dir) 
  }
})



