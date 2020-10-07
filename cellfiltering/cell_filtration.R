# Author: Michael Finlayson
#
# Description: 
# This class models the process of filtering cells in a scRNA sample
#
######## -----------------------------------------------------------------------
########                 Parameters                                     
######## ---------------------------------------------------------------

CFParameters = setClass('CFParameters', 
  slots = c(
   min_nUMI = 'numeric',
   max_nUMI = 'numeric',
   min_nGene = 'numeric',
   max_nGene = 'numeric',
   min_mt_prop = 'numeric',
   max_mt_prop = 'numeric'
  )
)

setMethod('initialize', 'CFParameters', 
function(.Object, 
         min_nUMI = -Inf, max_nUMI = Inf,
         min_nGene = -Inf, max_nGene = Inf,
         min_mt_prop = -Inf, max_mt_prop = Inf) {
  .Object@min_nUMI = min_nUMI
  .Object@max_nUMI = max_nUMI
  .Object@min_nGene = min_nGene
  .Object@max_nGene = max_nGene
  .Object@min_mt_prop = min_mt_prop
  .Object@max_mt_prop = max_mt_prop
  return(.Object)
})

######## -----------------------------------------------------------------------
########                 Definition                                     
######## ---------------------------------------------------------------

library(ggplot2)
library(reshape2)
library(cowplot)
library(MASS)

CellFiltration = setClass('CellFiltration', 
  slots = c(
    parameters = 'CFParameters',
    filtration_by_sample = 'list',
    cells_to_keep = 'character',
    functions = 'list'
  ),
  prototype = list(
    functions = list(
      
      init = function(cf, scdata, parameters) {
        
        # check inputs
        stopifnot(is(scdata, 'SCData'))
        stopifnot(all(c('sample', 'nUMI', 'nGene', 'mt_prop') %in% colnames(scdata@metadata)))
        stopifnot(is(parameters, 'CFParameters'))
        
        metadata = scdata@metadata
        
        filtration_by_sample = list()
        for (sample_name in unique(metadata$sample)) {
          
          filtration_data = metadata[metadata$sample == sample_name, c('nGene', 'nUMI', 'mt_prop'), drop=F]
          
          # remove based on minimums and maximums in the parameters
          filtration_data = cf@functions$apply_parameter_limits(filtration_data, parameters)
          
          # first thing, mt_prop
          filtration_data = cf@functions$get_mt_prop_probablities(filtration_data)
          filtration_data = cf@functions$apply_mt_prop_prob_limits(filtration_data)
          
          # Now, nUMI
          filtration_data = cf@functions$get_nUMI_probablities(filtration_data)
          filtration_data = cf@functions$apply_nUMI_prob_limits(filtration_data)
          
          # Now, nGene
          filtration_data = cf@functions$get_nGene_probabilities(filtration_data)
          filtration_data = cf@functions$apply_nGene_prob_limits(filtration_data)
          
          # get cells to keep by applying the cutoffs
          filtration_data = cf@functions$get_cells_to_keep(filtration_data)
          
          filtration_by_sample[[sample_name]] = filtration_data
        }
        
        # compose cell filtration object
        cf@filtration_by_sample = filtration_by_sample
        
        for(sample_name in names(filtration_by_sample)) { 
          sample_filtration_data = filtration_by_sample[[sample_name]]
          sample_cells_to_keep = rownames(sample_filtration_data)[sample_filtration_data$within_all_limits]
          cf@cells_to_keep = union(cf@cells_to_keep, sample_cells_to_keep) 
        }
        
        cf@parameters = parameters
        
        return(cf)
      },
      
      apply_parameter_limits = function(filtration_data, parameters) {
        
        gt_max_nUMI = filtration_data$nUMI > parameters@max_nUMI
        lt_min_nUMI = filtration_data$nUMI < parameters@min_nUMI
        gt_max_nGene = filtration_data$nGene > parameters@max_nGene
        lt_min_nGene = filtration_data$nGene < parameters@min_nGene
        gt_max_mt_prop = filtration_data$mt_prop > parameters@max_mt_prop
        lt_min_mt_prop = filtration_data$mt_prop < parameters@min_mt_prop
        
        beyond_limits = gt_max_nUMI | lt_min_nUMI | gt_max_nGene 
        beyond_limits = beyond_limits | lt_min_nGene | gt_max_mt_prop 
        beyond_limits = beyond_limits | lt_min_mt_prop
        
        filtration_data = data.frame(filtration_data, 
                                     gt_max_nUMI = gt_max_nUMI,
                                     lt_min_nUMI = lt_min_nUMI,
                                     gt_max_nGene = gt_max_nGene,
                                     lt_min_nGene = lt_min_nGene,
                                     gt_max_mt_prop = gt_max_mt_prop,
                                     lt_min_mt_prop = lt_min_mt_prop,
                                     within_parameter_limits = !beyond_limits)
        
        return(filtration_data)
      },
      
      get_mt_prop_probablities = function(filtration_data) {
        
        # determine which cells to fit on
        log_nUMI = log1p(filtration_data$nUMI)
        log_nUMI_median = median(log_nUMI)
        log_nUMI_sd = sqrt( sum( (log_nUMI - log_nUMI_median)^2 ) / (length(log_nUMI) - 1))
        
        nGene = filtration_data$nGene
        nGene_fit = fitdist(nGene, distr = 'norm')
        nGene_mean = nGene_fit$estimate['mean']
        nGene_sd = nGene_fit$estimate['sd']
        
        fit_data_select = abs(log_nUMI - log_nUMI_median) < log_nUMI_sd & abs(nGene - nGene_mean) < nGene_sd
        
        # log transform the data
        log_mt_prop = log1p(filtration_data$mt_prop)
        names(log_mt_prop) = rownames(filtration_data)
        
        # get data to fit on
        fit_data = log_mt_prop[fit_data_select]
        
        # fit the data
        fit_mean = median(fit_data)
        fit_sd = sqrt( sum( (fit_data - fit_mean)^2 ) / length(fit_data) )
        
        # compute the probabilities
        mt_prop_prob = pnorm(log_mt_prop, mean = fit_mean, sd = fit_sd)
        
        filtration_data = data.frame(filtration_data, mt_prop_prob=mt_prop_prob)
        
        return(filtration_data)
      }, 
      
      apply_mt_prop_prob_limits = function(filtration_data) {
        
        gt_mt_prop_max_prob_limit = filtration_data$mt_prop_prob > 0.96
        within_mt_prop_prob_limits = !(gt_mt_prop_max_prob_limit)
        
        filtration_data = data.frame(filtration_data, 
                                     gt_mt_prop_max_prob_limit=gt_mt_prop_max_prob_limit,
                                     within_mt_prop_prob_limits=within_mt_prop_prob_limits)
        
        return(filtration_data)
      },
      
      get_nUMI_probablities = function(filtration_data) {
        
        # determine which cells to fit on
        fit_data_select = filtration_data$within_mt_prop_prob_limits
        
        # log transform the data
        log_nUMI = log1p(filtration_data$nUMI)
        
        # get data to fit on
        fit_data = log_nUMI[fit_data_select]
        
        # fit the data
        median_log_nUMI = median(fit_data)
        sd_log_nUMI = sqrt( sum( (fit_data - median_log_nUMI)^2) / (length(log_nUMI) - 1) )
        
        # compute the probabilities
        nUMI_prob = pnorm(log_nUMI, mean = median_log_nUMI, sd = sd_log_nUMI)
        
        filtration_data = data.frame(filtration_data, nUMI_prob=nUMI_prob)
        return(filtration_data)
      },
      
      apply_nUMI_prob_limits = function(filtration_data) {
        
        gt_nUMI_max_prob_limit = filtration_data$nUMI_prob > 0.98
        lt_nUMI_min_prob_limit = filtration_data$nUMI_prob < 0.02
        within_nUMI_prob_limits = !(gt_nUMI_max_prob_limit | lt_nUMI_min_prob_limit)
        
        filtration_data = data.frame(filtration_data, 
                                     gt_nUMI_max_prob_limit=gt_nUMI_max_prob_limit,
                                     lt_nUMI_min_prob_limit=lt_nUMI_min_prob_limit,
                                     within_nUMI_prob_limits=within_nUMI_prob_limits)
        
        return(filtration_data)
      },
      
      get_nGene_probabilities = function(filtration_data) {
        
        # determine which cells to fit on
        fit_data_select = filtration_data$within_mt_prop_prob_limits
        
        nGene = filtration_data$nGene
        
        # get data to fit on
        fit_data = nGene[fit_data_select]
        
        # fit the data
        dist_fit = fitdistr(fit_data, densfun = 'normal')
        
        nGene_prob = pnorm(nGene, 
                           mean = dist_fit$estimate['mean'], 
                           sd = dist_fit$estimate['sd'])
        
        filtration_data = data.frame(filtration_data, nGene_prob=nGene_prob)
        return(filtration_data)
      },
      
      apply_nGene_prob_limits = function(filtration_data) {
        
        gt_nGene_max_prob_limit = filtration_data$nUMI_prob > 0.99
        lt_nGene_min_prob_limit = filtration_data$nUMI_prob < 0.01
        within_nGene_prob_limits = !(gt_nGene_max_prob_limit | lt_nGene_min_prob_limit)
        
        filtration_data = data.frame(filtration_data, 
                                     gt_nGene_max_prob_limit=gt_nGene_max_prob_limit,
                                     lt_nGene_min_prob_limit=lt_nGene_min_prob_limit,
                                     within_nGene_prob_limits=within_nGene_prob_limits)
        
        return(filtration_data)
      },
      
      get_cells_to_keep = function(filtration_data) {
        
        within_all_limits = filtration_data$within_parameter_limits 
        within_all_limits = filtration_data$within_mt_prop_prob_limits & within_all_limits
        within_all_limits = filtration_data$within_nUMI_prob_limits & within_all_limits
        within_all_limits = filtration_data$within_nGene_prob_limits & within_all_limits
        
        filtration_data = data.frame(filtration_data, within_all_limits=within_all_limits)
        return(filtration_data)
      }
    )
  )
)


######## -----------------------------------------------------------------------
########                 Initialization                                 
######## ---------------------------------------------------------------

setMethod('initialize', 'CellFiltration', 
function(.Object, scdata, parameters = CFParameters()) { 
  return(.Object@functions$init(.Object, scdata, parameters = parameters)) 
})


######## -----------------------------------------------------------------------
########                 Methods                                        
######## ---------------------------------------------------------------

###### --------------------
######     cf_nullify

setGeneric('cf_nullify', def=function(cf) {standardGeneric('cf_nullify')})

setMethod('cf_nullify', signature(cf='CellFiltration'), function(cf) { 
  cf@parameters = CFParameters()
  cf@filtration_by_sample = list()
  cf@cells_to_keep = c()
  return(cf)
})
  
###### --------------------
######     cf_get_plots

setGeneric('cf_get_plots', def=function(cf) {standardGeneric('cf_get_plots')})

setMethod('cf_get_plots', signature(cf='CellFiltration'), function(cf) { 
  
  plot_list = list()
  
  for (sample_name in names(cf@filtration_by_sample)) {
    
    filtration_data = cf@filtration_by_sample[[sample_name]]
    
    qual_metrics = filtration_data[,c('nGene','nUMI','mt_prop'),drop=F]
    
    qual_metrics_before = qual_metrics
    qual_metrics_before = data.frame(time=rep('before', times=nrow(qual_metrics_before)), 
                                     qual_metrics_before)
    
    qual_metrics_after = qual_metrics[filtration_data$within_all_limits,,drop=F]
    qual_metrics_after = data.frame(time=rep('after', times=nrow(qual_metrics_after)), 
                                    qual_metrics_after)
    
    pdata = rbind(qual_metrics_before, qual_metrics_after)
    pdata[,'mt_prop'] = pdata[,'mt_prop'] + 0.00001
    pdata[,'nUMI'] = pdata[,'nUMI'] + 0.001
    
    stroke_size = 1
    alpha_val = .8
    if (nrow(qual_metrics) > 50) {
      stroke_size = .5
      alpha_val = .5
    }
    if (nrow(qual_metrics) > 200) {
      stroke_size = .2
      alpha_val = .2
    }
    if (nrow(qual_metrics) > 1000) {
      stroke_size = .05
      alpha_val = .01
    }
    if (nrow(qual_metrics) > 2000) {
      stroke_size = .01
      alpha_val = .001
    }
    
    plots = list()
    for (metric in colnames(qual_metrics)) {
      
      ggp = ggplot(melt(pdata[,c('time',metric)], id.vars = 'time'), aes(x=time, y=value)) 
      ggp = ggp + geom_violin(aes(fill='red')) + geom_jitter(aes(alpha = .2, stroke = stroke_size))
      
      if (metric %in% c('nUMI', 'mt_prop')) { ggp = ggp + scale_y_continuous(trans='log10') }
      
      plots[[metric]] = ggp + ggtitle(metric) + scale_fill_discrete(guide=F) + scale_alpha(guide=F)
      
    }
    filtration_plot = plot_grid(plotlist=plots, ncol=3)
    
    title = ggdraw() + draw_label(sample_name, fontface='bold')
    filtration_plot = plot_grid(title, filtration_plot, ncol=1, rel_heights=c(0.1, 1)) 
    
    plot_list[[sample_name]] = filtration_plot
  }
  
  return(plot_list)
})

###### ---------------
######     cf_save

setGeneric('cf_save', def=function(cf, save_dir, ...) 
{standardGeneric('cf_save')})

setMethod('cf_save', 
signature(cf='CellFiltration', save_dir='character'), 
function(cf, save_dir, dev_to_use = 'pdf') { 
  
  # work out that trailing slash ambiguity
  save_dir = paste(strsplit(save_dir, '/')[[1]], collapse='/')
  
  stopifnot(!dir.exists(save_dir) && !(exists(save_dir)))  # don't overwrite with files or directories
  dir.create(save_dir)
  
  # save parameters
  # parameters_filename = paste(save_dir, 'cell_filtration_parameters.csv', sep = '/')
  # ptable = data.frame(parameter=c('min_nUMI', 'max_nUMI', 'min_nGene', 'max_nGene','min_mt_prop','max_mt_prop'),
  #                     value=c(cf@parameters@min_nUMI, 
  #                             cf@parameters@max_nUMI, 
  #                             cf@parameters@min_nGene, 
  #                             cf@parameters@max_nGene, 
  #                             cf@parameters@min_mt_prop, 
  #                             cf@parameters@max_mt_prop))
  # write.csv(ptable, file = parameters_filename)
  
  # save plots and filtration data for each sample
  plot_list = cf_get_plots(cf)
  for (sample in names(plot_list)) {
    
    sample_dir = paste(save_dir, sample, sep = '/')
    dir.create(sample_dir)
    
    # save the filtration data
    fdata_filename = paste(sample, 'filtration_data.csv', sep = '_')
    fdata_filename = paste(sample_dir, fdata_filename, sep = '/')
    
    filtration_data = cf@filtration_by_sample[[sample]]
    filtration_data = data.frame(cell_name=rownames(filtration_data), filtration_data)
    
    write.csv(x = filtration_data, file = fdata_filename, row.names = F)
    
    # save the plot
    plot_filename = paste(sample, paste('_filtration', dev_to_use, sep = '.'), sep = '')
    plot_filename = paste(sample_dir, plot_filename, sep = '/')
    
    cf_plot = plot_list[[sample]]
    
    ggsave(filename = plot_filename, plot = cf_plot, device = dev_to_use, height = 7, width = 9)
  }
})

