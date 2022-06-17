args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]

source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/Get_Response.R")

clin = read.csv( file.path(input_dir, "CLIN.txt"), stringsAsFactors=FALSE , sep="\t" )
clin$Sample.ID = paste( "P." , clin$Sample.ID , sep="" )
rownames(clin) = clin$Sample.ID
clin$PFS = as.numeric(as.character(clin$PFS))
clin$Clinical.benefit = toupper(clin$Clinical.benefit)

clin = clin[ , c( "Sample.ID" , "PFS", "PD_Event.1_Censoring.0" , "Clinical.benefit" )]
colnames(clin) = c( "patient" , "t.pfs" , "pfs" , "response.other.info"  )

sample = t( read.csv( file.path(input_dir, "CLIN_Sample.txt"), stringsAsFactors=FALSE , sep="\t" , header=FALSE ) ) 
colnames(sample) = as.character( sample[ 1 , ] )
sample = sample[ -1 , ] 

sample = as.data.frame( sample )
sample[ , "patient"] = paste( "P." , sapply( as.character( sample[ , "patient"] ) , function(x){ unlist( strsplit( x , " " , fixed=TRUE ))[2] } ) ,sep="")
sample[ , "sex"] = sapply( as.character( sample[ , "sex"] ) , function(x){ unlist( strsplit( x , ": " , fixed=TRUE ))[2] } )
sample[ , "age"] = sapply( as.character( sample[ , "age"] ) , function(x){ unlist( strsplit( x , ": " , fixed=TRUE ))[2] } )
sample[ , "sex"] = ifelse( sample[ , "sex"] %in% "male" , "M" , "F" )
rownames(sample) = sample$patient

clin = cbind( clin , sample[ clin$patient , c( "sex" , "age" ) ] , "Lung" , "PD-1/PD-L1" , NA , NA , NA , NA , NA, NA , NA , NA )
colnames(clin) = c( "patient" , "t.pfs" , "pfs", "response.other.info" , "sex" , "age" , "primary" , "drug_type"  , "response" , "recist" , "t.os" ,"os" ,"histo" , "stage" , "dna" , "rna" )

case = read.csv( file.path(output_dir, "cased_sequenced.csv"), stringsAsFactors=FALSE , sep=";" )
clin$rna[ clin$patient %in% case[ case$expr %in% 1 , ]$patient ] = "tpm"
clin$dna[ clin$patient %in% case[ case$cna %in% 1 , ]$patient ] = "wes"

clin$response.other.info = ifelse( clin$response.other.info %in% "DCB" , "R", 
							ifelse( clin$response.other.info %in% "NDB" , "NR" , NA ))

clin$response = Get_Response( data=clin )


clin = clin[ , c("patient" , "sex" , "age" , "primary" , "histo" , "stage" , "response.other.info" , "recist" , "response" , "drug_type" , "dna" , "rna" , "t.pfs" , "pfs" , "t.os" , "os" ) ]

write.table( clin , file=file.path(output_dir, "CLIN.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=FALSE )

