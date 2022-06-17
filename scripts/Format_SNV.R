library(data.table)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]

case = read.csv( file.path(output_dir, "cased_sequenced.csv"), stringsAsFactors=FALSE , sep=";" )

snv = as.data.frame( fread( file.path(input_dir, "SNV.csv.gz") , stringsAsFactors=FALSE , sep=";" ))

data = cbind( snv[ , c("Start" , "Chromosome" , "Sample ID" , "Ref" , "Alt" , "Function" ) ] , 
				sapply( snv[ , "Gene" ] , function(x){ unlist( strsplit(  unlist( strsplit( x , "(" , fixed=TRUE ))[1] , ":" , fixed=TRUE ))[1] } )
			)
colnames(data) = c( "Pos" , "Chr" , "Sample" , "Ref" , "Alt" , "Effect" , "Gene"  )

data$Ref = ifelse( data$Ref %in% "-" , "" , data$Ref )
data$Alt = ifelse( data$Alt %in% "-" , "" , data$Alt )

data = cbind ( data , apply( data[ , c( "Ref", "Alt" ) ] , 1 , function(x){ ifelse( nchar(x[1]) != nchar(x[2]) , "INDEL", "SNV") } ) )

colnames(data) = c( "Pos" , "Chr" , "Sample" , "Ref" , "Alt" , "Effect" , "Gene" , "MutType"  )

data$Sample = paste( "P." , data$Sample , sep="" ) 
data$Effect[ data$Effect %in% "frameshift deletion" ] = "Frame_Shift_Del"
data$Effect[ data$Effect %in% "nonsynonymous" ] = "Missense_Mutation"
data$Effect[ data$Effect %in% "splicing" ] = "Splice_Site"
data$Effect[ data$Effect %in% "stopgain" ] = "Nonsense_Mutation"
data$Effect[ data$Effect %in% "stoploss" ] = "Nonstop_Mutation"

data = data[ data$Sample %in% case[ case$snv %in% 1 , ]$patient , c( "Sample" , "Gene" , "Chr" , "Pos" , "Ref" , "Alt" , "Effect" , "MutType" ) ]


write.table( data , file=file.path(output_dir, "SNV.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=FALSE )