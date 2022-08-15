library(stringr)
library(tibble)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]
annot_dir <- args[3]

source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/Get_Response.R")
source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/format_clin_data.R")
source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/annotate_tissue.R")
source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/annotate_drug.R")

clin_original = read.csv( file.path(input_dir, "CLIN.txt"), stringsAsFactors=FALSE , sep="\t" )
clin_original$Sample.ID = paste( "P." , clin_original$Sample.ID , sep="" )
clin <- clin_original
rownames(clin) = clin$Sample.ID
clin$PFS = as.numeric(as.character(clin$PFS))
clin$Clinical.benefit = toupper(clin$Clinical.benefit)

selected_cols <- c( "Sample.ID" , "PFS", "PD_Event.1_Censoring.0" , "Clinical.benefit" )
clin = clin[ , selected_cols]
colnames(clin) = c( "patient" , "t.pfs" , "pfs" , "response.other.info"  )

sample = read.csv( file.path(input_dir, "CLIN_Sample.txt"), stringsAsFactors=FALSE , sep="\t" )  
sample$patient <- str_replace(sample$patient, 'NSCLC ', 'P.')
rownames(sample) = sample$patient
sample$sex = ifelse( sample$sex %in% "male" , "M" , "F" )
# colnames(sample) = as.character( sample[ 1 , ] )
# sample = sample[ -1 , ] 
# sample = as.data.frame( sample )
# sample[ , "patient"] = paste( "P." , sapply( as.character( sample[ , "patient"] ) , function(x){ unlist( strsplit( x , " " , fixed=TRUE ))[2] } ) ,sep="")
# sample[ , "sex"] = sapply( as.character( sample[ , "sex"] ) , function(x){ unlist( strsplit( x , ": " , fixed=TRUE ))[2] } )
# sample[ , "age"] = sapply( as.character( sample[ , "age"] ) , function(x){ unlist( strsplit( x , ": " , fixed=TRUE ))[2] } )
# sample[ , "sex"] = ifelse( sample[ , "sex"] %in% "male" , "M" , "F" )

clin = cbind( clin , sample[ clin$patient , c( "sex" , "age" ) ] , "Lung" , "PD-1/PD-L1" , NA , NA , NA , NA , NA, NA , NA , NA )
colnames(clin) = c( "patient" , "t.pfs" , "pfs", "response.other.info" , "sex" , "age" , "primary" , "drug_type"  , "response" , "recist" , "t.os" ,"os" ,"histo" , "stage" , "dna" , "rna" )

case = read.csv( file.path(output_dir, "cased_sequenced.csv"), stringsAsFactors=FALSE , sep=";" )
clin$rna[ clin$patient %in% case[ case$expr %in% 1 , ]$patient ] = "tpm"
clin$dna[ clin$patient %in% case[ case$cna %in% 1 , ]$patient ] = "wes"

clin$response.other.info = ifelse( clin$response.other.info %in% "DCB" , "R", 
							ifelse( clin$response.other.info %in% "NDB" , "NR" , NA ))

clin$response = Get_Response( data=clin )

clin = clin[ , c("patient" , "sex" , "age" , "primary" , "histo" , "stage" , "response.other.info" , "recist" , "response" , "drug_type" , "dna" , "rna" , "t.pfs" , "pfs" , "t.os" , "os" ) ]

clin <- format_clin_data(clin_original, 'Sample.ID', selected_cols, clin)

# Tissue and drug annotation
annotation_tissue <- read.csv(file=file.path(annot_dir, 'curation_tissue.csv'))
clin <- annotate_tissue(clin=clin, study='Jung', annotation_tissue=annotation_tissue, check_histo=FALSE)

annotation_drug <- read.csv(file=file.path(annot_dir, 'curation_drug.csv'))
clin <- add_column(clin, unique_drugid='', .after='unique_tissueid')

write.table( clin , file=file.path(output_dir, "CLIN.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=FALSE )

