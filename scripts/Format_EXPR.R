library(data.table)
library(biomaRt)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]

tpm = as.matrix( fread( file.path(input_dir, "EXPR.txt.gz") , stringsAsFactors=FALSE , sep="\t" , dec="," ))

rownames(tpm) = sapply( tpm[,1] , function(x){ unlist( strsplit( x , "." , fixed=TRUE ) )[1] } )
colnames(tpm) = sapply( colnames(tpm) , function(x){ paste( "P." , unlist( strsplit( x , "NSCLC" , fixed=TRUE ) )[2] , sep="") } )
tpm = tpm[,-1]

rid = rownames(tpm)
cid = colnames(tpm)
tpm = apply(apply(tpm,2,as.character),2,as.numeric)
colnames(tpm) = cid
rownames(tpm) = rid

##################################
##################################
genes <- rownames(tpm)
human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
geneID = as.matrix( getBM(attributes=c("hgnc_symbol","ensembl_gene_id"), filters="ensembl_gene_id",values=genes, mart=human) )
rownames(geneID) = geneID[ , 2 ]

geneID = geneID[ !( rownames(geneID) %in% rownames( geneID[ duplicated( rownames(geneID) ) , ] ) ) , ]

tpm = tpm[ rownames(tpm) %in% rownames(geneID) , ]
rownames(tpm) = geneID[ rownames(tpm) , 1 ]

tpm = unique( tpm[ !(rownames(tpm) %in% "") & rowSums(tpm) > 0 , ] )


##################################
##################################
##Remove Duplicated Genes
t_uniq <- tpm[ !( rownames(tpm) %in% rownames( tpm )[ duplicated( rownames( tpm ) ) ] ) , ]
t_dup <- tpm[ ( rownames(tpm) %in% rownames( tpm )[ duplicated( rownames( tpm ) ) ] ) , ]

t_dup <- t_dup[ order( rownames( t_dup ) ) , ]
t.dup.sd <- apply(t_dup,1,sd)
id <- unique(names(t.dup.sd))

t.dup.sd.rm <- NULL
for(i in 1:length(id)){
	tmp <- t.dup.sd[which(names(t.dup.sd)%in%id[i])]
	sd <- max(tmp)		
	t.dup.sd.rm <- c(t.dup.sd.rm,ifelse(tmp==sd,TRUE,FALSE))	
}
tpm <- rbind( t_uniq , t_dup[ t.dup.sd.rm , ] )

#############################################################################
#############################################################################

case = read.csv( file.path(output_dir, "cased_sequenced.csv"), stringsAsFactors=FALSE , sep=";" )
tpm = tpm[ , colnames(tpm) %in% case[ case$expr %in% 1 , ]$patient ]

write.table( tpm , file= file.path(output_dir, "EXPR.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=TRUE )
