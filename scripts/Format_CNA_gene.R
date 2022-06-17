library(data.table)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]

case = read.csv( file.path(output_dir, "cased_sequenced.csv"), stringsAsFactors=FALSE , sep=";" )

cna = as.data.frame( fread( file.path(input_dir, "gistic/all_thresholded.by_genes.txt.gz") , stringsAsFactors=FALSE , sep="\t" , dec="," ))
rownames(cna) = cna[ , 1]
cna = cna[ , -( 1:3 ) ]

cna = unique( cna[ , colnames(cna) %in% case[ case$cna %in% 1 , ]$patient ] )

rownames(cna)[ grep( "|" , rownames(cna) , fixed=TRUE ) ] = sapply( rownames(cna)[ grep( "|" , rownames(cna) , fixed=TRUE ) ] , function(x) unlist( strsplit( x , "|" , fixed=TRUE ) )[1] )

write.table( cna , file=file.path(output_dir, "CNA_gene.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=TRUE )