setwd('~/cscac_analyses/jj_190513_JJ_NICOLE_1_HUMAN_10X/')
source('~/cscac_analyses/analysis_code/cscac_public_source_awR_shared.R')

scd = readRDS('rds_files/scd.rds')

# combined

scdp = SCDPrep(scdata = scd)
scdp = scdprep_filter(scdprep = scdp)
scdp = scdprep_vst(scdprep = scdp, parameters = VSTDParameters(batch_var = 'sample'))
scdp = scdprep_select_info_genes(scdprep = scdp, parameters = IGSParameters(collect_method = 'intersection'))

saveRDS(scdp, 'rds_files/scdp.rds')

scd_list = scdata_split(scdata = scd, 'sample')

# JN001

scdp_jn001 = SCDPrep(scd_list$JN001)
scdp_jn001 = scdprep_filter(scdprep = scdp_jn001)
scdp_jn001 = scdprep_vst(scdprep = scdp_jn001)
scdp_jn001 = scdprep_select_info_genes(scdprep = scdp_jn001)

saveRDS(scdp_jn001, 'rds_files/scdp_jn001.rds')

# JJ002

scdp_jj002 = SCDPrep(scd_list$JJ002)
scdp_jj002 = scdprep_filter(scdprep = scdp_jj002)
scdp_jj002 = scdprep_vst(scdprep = scdp_jj002)
scdp_jj002 = scdprep_select_info_genes(scdprep = scdp_jj002)

saveRDS(scdp_jj002, 'rds_files/scdp_jj002.rds')
