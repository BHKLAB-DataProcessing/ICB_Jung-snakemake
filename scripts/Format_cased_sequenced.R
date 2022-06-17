library(data.table)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]

clin = read.csv( file.path(input_dir, "CLIN.txt"), stringsAsFactors=FALSE , sep="\t" )
rownames(clin) = clin$Sample.ID

snv = as.data.frame( fread( file.path(input_dir, "SNV.csv.gz") , stringsAsFactors=FALSE  , sep=";"))
snv = sort( unique( snv[ , "Sample ID" ] ))

cna = as.data.frame( fread( file.path(input_dir, "CNA_seg.txt.gz") , stringsAsFactors=FALSE , sep="\t" , dec="," ))
colnames(cna) = c( "sample" , "chrom" , "loc.start" , "loc.end" , "num.mark" , "seg.mean" )
cna = cna[ !is.na(cna$sample) , ]
cna = sort( unique( cna$sample ))

rna = sapply( colnames( as.matrix( fread( file.path(input_dir, "EXPR.txt.gz") , stringsAsFactors=FALSE  , sep="\t" , dec=',') ) )[ -1 ] , 
				function(x){ unlist( strsplit( x , "NSCLC" , fixed=TRUE ) )[2] } )

patient = sort( unique( clin$Sample.ID ) )

case = as.data.frame( cbind( patient , rep( 0 , length(patient) ) , rep( 0 , length(patient) ) , rep( 0 , length(patient) ) ) )
colnames(case) = c( "patient" , "snv" , "cna" , "expr" )
rownames(case) = patient

case$snv = as.numeric( as.character( case$snv ) )
case$cna = as.numeric( as.character( case$cna ) )
case$expr = as.numeric( as.character( case$expr ) )


for( i in 1:nrow(case)){
	if( rownames(case)[i] %in% snv ){
		case$snv[i] = 1
	}
	if( rownames(case)[i] %in% cna ){
		case$cna[i] = 1
	}
	if( rownames(case)[i] %in% rna ){
		case$expr[i] = 1
	}
}

case$patient = paste( "P." , case$patient , sep="" )

write.table( case , file=file.path(output_dir, "cased_sequenced.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=FALSE )


