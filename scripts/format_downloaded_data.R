library(data.table)
library(readxl) 
library(stringr)
library(GEOquery)
library(tabulizer) # needs to install from repo: remotes::install_github(c("ropensci/tabulizerjars", "ropensci/tabulizer"))

args <- commandArgs(trailingOnly = TRUE)
work_dir <- args[1]
annot_dir <- args[2]

# CLIN.txt, CLIN_Sample.txt
out <- extract_tables(file.path(work_dir, "41467_2019_12159_MOESM1_ESM.pdf"), pages = c(12))
clin <- data.frame(out[[1]])
colnames(clin) <- clin[1, ]
clin <- clin[-1, ]
colnames(clin) <- str_replace_all(colnames(clin), '\\W', '.')
clin[ clin == 'na'] <- NA
clin[, colnames(clin)[colnames(clin) != 'Clinical.benefit']] <- sapply(clin[, colnames(clin)[colnames(clin) != 'Clinical.benefit']], as.numeric)
write.table(clin, file.path(work_dir, 'CLIN.txt'), col.names = TRUE, sep='\t')

gunzip(file.path(work_dir, "GSE135222_series_matrix.txt.gz"))
clin_sample <- getGEO(filename=file.path(work_dir, 'GSE135222_series_matrix.txt'), destdir=work_dir)
clin_sample <- pData(clin_sample)
colnames(clin_sample) <- str_replace_all(colnames(clin_sample), '\\W', '.')
clin_sample <- clin_sample[, c('title', 'geo_accession', 'source_name_ch1', 'age.ch1', 'gender.ch1', 'pfs.time.ch1', 'progression.free.survival..pfs..ch1')]
colnames(clin_sample) <- c("patient", "sample", "primary", "age", "sex", "pfs.time", 'pfs')
write.table(clin_sample, file.path(work_dir, 'CLIN_Sample.txt'), col.names = TRUE, row.names=FALSE, sep='\t')

file.remove(file.path(work_dir, 'GPL16791.soft'))
file.remove(file.path(work_dir, 'GSE135222_series_matrix.txt'))

# expr_list.rds
source('https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/process_kallisto_output.R')
load(file.path(annot_dir, "Gencode.v40.annotation.RData"))
process_kallisto_output(work_dir, 'Jung_kallisto.zip', tx2gene)

expr_list <- readRDS(file.path(work_dir, 'expr_list.rds'))
rnaseq_samples <- read.table(file.path(work_dir, 'rnaseq_samples.tsv'), header=TRUE, sep='\t')

for(assay_name in names(expr_list)){
  expr <- data.frame(expr_list[[assay_name]])
  colnames(expr) <- unlist(lapply(colnames(expr), function(col){
    return(
      str_replace(rnaseq_samples$sample_title[rnaseq_samples$run_accession == col], '\\W', '_')
    )
  }))
  expr_list[[assay_name]] <- expr
}

saveRDS(expr_list, file.path(work_dir, 'expr_list.rds'))

# SNV.csv.gz
snv <- read_excel(file.path(work_dir, '41467_2019_12159_MOESM8_ESM.xlsx'), sheet='Sheet1')
colnames(snv) <- snv[3, ]
snv <- snv[-(1:3), ]
snv[, c('Start', 'End', 'Sample ID')] <- sapply(snv[, c('Start', 'End', 'Sample ID')], as.numeric)
gz <- gzfile(file.path(work_dir, 'SNV.csv.gz'), "w")
write.table( snv , file=gz , quote=FALSE , sep=";" , col.names=TRUE, row.names = FALSE )
close(gz)

# CNA_seg.txt.gz
cna <- read_excel(file.path(work_dir, '41467_2019_12159_MOESM9_ESM.xlsx'), sheet='Sheet1')
colnames(cna) <- cna[3, ]
cna <- cna[-(1:3), ]
cna <- sapply(cna, as.numeric)
gz <- gzfile(file.path(work_dir, 'CNA_seg.txt.gz'), "w")
write.table( cna , file=gz , quote=FALSE , sep="\t" , col.names=TRUE, row.names = FALSE )
close(gz)

