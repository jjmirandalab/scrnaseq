# Author: Michael Finlayson
#
# Description:
# This object will create a list of umi_counts and metadata for each library 
# that can be found within or below the current directory or directory 
# provided
#


######## -----------------------------------------------------------------------
########                 Definition                                     
######## ---------------------------------------------------------------

library(Matrix)
library(Seurat)

SCData = setClass('SCData', 
  slots=c(
    umi_counts = 'Matrix', 
    metadata = 'data.frame',
    functions = 'list'
    ),
  prototype = list(
    functions = list(
      set_count_derived_metadata = function(umi_counts, metadata) {
        
        colnames_keep = !(colnames(metadata) %in% c('nGene', 'nUMI', 'mt_prop'))
        metadata = metadata[,colnames_keep,drop=F]
        
        # create count derived metadata
        nUMI = colSums(umi_counts)
        nGene = colSums(umi_counts != 0)
        
        mito_genes = grep(pattern="^mt-", x=rownames(umi_counts), ignore.case=T, value=T)
        mt_prop = colSums(umi_counts[mito_genes,,drop=F]) / colSums(umi_counts)
        
        # combine
        if (!is.null(metadata)) {
          
          metadata = data.frame(nUMI=nUMI, nGene=nGene, mt_prop=mt_prop, metadata)
          
        } else {
          
          metadata = data.frame(nUMI=nUMI, nGene=nGene, mt_prop=mt_prop)
        }
        
        metadata = droplevels(metadata)
        
        return(metadata)
      }
    )
  )
)


######## -----------------------------------------------------------------------
########                 Initialization                                 
######## ---------------------------------------------------------------

setMethod('initialize', 'SCData', function(.Object, umi_counts, metadata=NULL) {
  
  if (!is.null(metadata)) { stopifnot(all(colnames(umi_counts)==rownames(metadata))) }
  
  metadata = .Object@functions$set_count_derived_metadata(umi_counts, metadata)
  
  .Object@umi_counts = umi_counts
  .Object@metadata = metadata
  
  return(.Object)
})

######## -----------------------------------------------------------------------
########                 Methods                                        
######## ---------------------------------------------------------------

###### ------------------------------
######     scdata_limit_to_cells

setGeneric('scdata_limit_to_cells', def=function(scdata, cells_to_keep) {standardGeneric('scdata_limit_to_cells')})

setMethod('scdata_limit_to_cells', signature(scdata='SCData', cells_to_keep='character'), 
function(scdata, cells_to_keep) {
  
  stopifnot(all(cells_to_keep %in% colnames(scdata@umi_counts)))
  
  return(SCData(umi_counts=scdata@umi_counts[,cells_to_keep,drop=F], 
                metadata=scdata@metadata[cells_to_keep,,drop=F]))
})

###### ---------------------
######     scdata_limit

setGeneric('scdata_limit', def=function(scdata, cells_to_keep, genes_to_keep) { standardGeneric('scdata_limit') })

setMethod('scdata_limit', signature(scdata='SCData', cells_to_keep='character', genes_to_keep='character'), 
function(scdata, cells_to_keep, genes_to_keep) {
  
  stopifnot(all(cells_to_keep %in% colnames(scdata@umi_counts)))
  stopifnot(all(genes_to_keep %in% rownames(scdata@umi_counts)))
  
  return(SCData(umi_counts=scdata@umi_counts[genes_to_keep,cells_to_keep,drop=F], 
                metadata=data.frame(scdata@metadata[cells_to_keep,,drop=F])))
})

###### --------------------------
######     scdata_downsample

setGeneric('scdata_downsample', def=function(scdata, ...) {standardGeneric('scdata_downsample')})

setMethod('scdata_downsample', signature(scdata='SCData'), 
function(scdata, max_umi=NULL) {
  
  stopifnot(all(rownames(scdata@metadata)==colnames(scdata@umi_counts)))
  
  if (is.null(max_umi)) {
    max_umi = median(scdata@metadata$nUMI)
  }
  
  umi_counts = SampleUMI(scdata@umi_counts, max.umi = floor(max_umi))
  colnames(umi_counts) = colnames(scdata@umi_counts)
  rownames(umi_counts) = rownames(scdata@umi_counts)
  
  # count derived metadata recalculated upon SCData creation
  
  return(SCData(umi_counts=umi_counts, metadata=scdata@metadata))
})

###### --------------------------
######     scdata_split

setGeneric('scdata_split', def=function(scdata, factor, ...) {standardGeneric('scdata_split')})

setMethod('scdata_split', signature(scdata='SCData', factor='character'), 
function(scdata, factor, value=NULL) {
  
  stopifnot(all(rownames(scdata@metadata)==colnames(scdata@umi_counts)))
  stopifnot(factor %in% colnames(scdata@metadata))
  if (!is.null(value)) { value %in% unique(scdata@metadata[,factor]) }

  values = value
  if (is.null(value)) { values = unique(scdata@metadata[,factor]) }
  
  scdata_list = list()
  for (v in values) {
    v_cell_names = rownames(scdata@metadata)[scdata@metadata[,factor]==v]
    scdata_list[[v]] = scdata_limit_to_cells(scdata = scdata, 
                                               cells_to_keep = v_cell_names)
  }
  
  return_value = scdata_list
  if (!is.null(value)) { return_value = scdata_list[[1]] }
  
  return(return_value)
})

###### ----------------------
######     scdata_combine

setGeneric('scdata_combine', def=function(scdata_list) {standardGeneric('scdata_combine')})

setMethod('scdata_combine', signature(scdata_list='list'), 
function(scdata_list) {
  
  # get first library
  umi_counts = scdata_list[[1]]@umi_counts
  metadata = scdata_list[[1]]@metadata
  
  # add the other libraries
  scdata_list = scdata_list[-1]
  for (scdata in scdata_list) {
    umi_counts = cbind(umi_counts, scdata@umi_counts)
    metadata = rbind(metadata, scdata@metadata)
  }
  
  return(SCData(umi_counts=umi_counts, metadata=metadata))
})

###### ------------------------------
######     scdata_save

setGeneric('scdata_save', def=function(scdata, save_dir) 
{standardGeneric('scdata_save')})

setMethod('scdata_save', signature(scdata='SCData', save_dir='character'), 
function(scdata, save_dir) {
  
  stopifnot(dir.exists(save_dir))

  save_dir = paste(strsplit(save_dir, split = '/')[[1]], collapse = '/')
  
  # save the corrected count matrix
  umi_counts_filename = paste(save_dir, 'umi_counts.csv.bz2', sep = '/')
  if (file.exists(umi_counts_filename)) { stop('umi counts file exists already') }
  
  write.csv(scdata@umi_counts, 
            file = bzfile(umi_counts_filename), 
            row.names = T)
  
  # save the metadata
  metadata_filename = paste(save_dir, 'metadata.csv', sep = '/')
  if (file.exists(metadata_filename)) { stop('metadata file exists already') }
  
  write.csv(scdata@metadata, 
            file = metadata_filename, 
            row.names = T)
})



