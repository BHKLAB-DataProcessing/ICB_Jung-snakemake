library(data.table)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]

cna = as.data.frame( fread( file.path(input_dir, "CNA_seg.txt.gz") , stringsAsFactors=FALSE , sep="\t" , dec="," ))
colnames(cna) = c( "sample" , "chrom" , "loc.start" , "loc.end" , "num.mark" , "seg.mean" )

cna$chrom = sapply( cna$chrom , function(x){ paste( "chr" , x , sep="" ) } )

cna = unique( cna[ !is.na(cna$sample) , ] )
cna$sample = paste( "P." , cna$sample , sep="" )

write.table( cna , file= file.path(output_dir, "CNA_seg.txt") , quote=FALSE , sep="\t" , col.names=TRUE , row.names=FALSE )