setwd('~/cscac_analyses/jj_190513_JJ_NICOLE_1_HUMAN_10X/')
source('~/cscac_analyses/analysis_code/cscac_public_source_awR_shared.R')
source('~/cscac_analyses/analysis_code/utility_functions.R')

# read in viral counts
jn1i = Matrix(as.matrix(read.table('experiment_data/JN001/JN001_viral_counts_intersection.tsv', header = T, sep = '\t', row.names = 1)))
colnames(jn1i) = paste(colnames(jn1i), '-1', sep = '')

jj2i = Matrix(as.matrix(read.table('experiment_data/JJ002/JJ002_viral_counts_intersection.tsv', header = T, sep = '\t', row.names = 1)))
colnames(jj2i) = paste(colnames(jj2i), '-2', sep = '')

shared_vgenes = intersect(rownames(jn1i), rownames(jj2i))

# read in regular counts data
JJ002 = readMM(file = 'experiment_data/JJ002/filtered_feature_bc_matrix/matrix.mtx.gz')
colnames(JJ002) = read.table(file = 'experiment_data/JJ002/filtered_feature_bc_matrix/barcodes.tsv.gz', stringsAsFactors = F)[,1]
colnames(JJ002) = sub(pattern = '-1', replacement = '-2', colnames(JJ002))

JN001 = readMM(file = 'experiment_data/JN001/filtered_gene_bc_matrices/GRCh38_EBV_KSHV/matrix.mtx')
colnames(JN001) = read.table(file = 'experiment_data/JN001/filtered_gene_bc_matrices/GRCh38_EBV_KSHV/barcodes.tsv', stringsAsFactors = F)[,1]

gene_names_jj002 = read.table(file = 'experiment_data/JJ002/filtered_feature_bc_matrix/features.tsv.gz', stringsAsFactors = F)
gene_names_jn001 = read.table(file = 'experiment_data/JN001/filtered_gene_bc_matrices/GRCh38_EBV_KSHV/genes.tsv', stringsAsFactors = F)
all(gene_names_jj002[,2]==gene_names_jn001[,2])

gene_names = make.unique(gene_names_jj002[,2])
rownames(JJ002) = gene_names
rownames(JN001) = gene_names

# metadata
md = data.frame(sample=c(rep('JN001', times=ncol(JN001)), rep('JJ002', times=ncol(JJ002))))
rownames(md) = c(colnames(JN001), colnames(JJ002))

# add viral counts
JJ002 = rbind(JJ002, jj2i[shared_vgenes,colnames(JJ002)])
JN001 = rbind(JN001, jn1i[shared_vgenes,colnames(JN001)])
                         
scd_combined = SCData(umi_counts = cbind(JN001, JJ002), metadata = md)

ebv_counts = colSums(scd_combined@umi_counts[grep('EBV_', shared_vgenes, value = T),])
kshv_counts = colSums(scd_combined@umi_counts[grep('KSHV_', shared_vgenes, value = T),])

scd_combined@metadata = data.frame(scd_combined@metadata, ebv_counts=ebv_counts, kshv_counts=kshv_counts)
scd_combined@umi_counts = scd_combined@umi_counts[-which(rownames(scd_combined@umi_counts) %in% c('EBV', 'KSHV')),]

saveRDS(object = scd_combined, file = 'rds_files/scd.rds')
