library(data.table)
library(biomaRt)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]

expr_list <- readRDS(file.path(input_dir, 'expr_list.rds'))
case = read.csv( file.path(output_dir, "cased_sequenced.csv"), stringsAsFactors=FALSE , sep=";" )

for(assay_name in names(expr_list)){
  expr <- data.frame(expr_list[[assay_name]])
  colnames(expr) <- str_replace(colnames(expr), 'NSCLC_', 'P.')
  expr <- expr[ , colnames(expr) %in% case[ case$expr %in% 1 , ]$patient ]
  write.table( 
    expr , 
    file= file.path(output_dir, paste0('EXPR_', str_replace(assay_name, 'expr_', ''), '.csv')) , 
    quote=FALSE , 
    sep=";" , 
    col.names=TRUE , 
    row.names=TRUE 
  )
}
