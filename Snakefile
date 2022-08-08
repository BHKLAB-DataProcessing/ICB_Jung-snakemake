from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
S3 = S3RemoteProvider(
    access_key_id=config["key"], 
    secret_access_key=config["secret"],
    host=config["host"],
    stay_on_remote=False
)
prefix = config["prefix"]
filename = config["filename"]

rule get_MultiAssayExp:
    output:
        S3.remote(prefix + filename)
    input:
        S3.remote(prefix + "processed/CLIN.csv"),
        # S3.remote(prefix + "processed/CNA_seg.txt"),
        # S3.remote(prefix + "processed/CNA_gene.csv"),
        S3.remote(prefix + "processed/EXPR_gene_tpm.csv"),
        S3.remote(prefix + "processed/EXPR_gene_counts.csv"),
        S3.remote(prefix + "processed/EXPR_isoform_tpm.csv"),
        S3.remote(prefix + "processed/EXPR_isoform_counts.csv"),
        S3.remote(prefix + "processed/SNV.csv"),
        S3.remote(prefix + "processed/cased_sequenced.csv"),
        S3.remote(prefix + "annotation/Gencode.v40.annotation.RData")
    resources:
        mem_mb=3000,
        disk_mb=3000
    shell:
        """
        Rscript -e \
        '
        load(paste0("{prefix}", "annotation/Gencode.v40.annotation.RData"))
        source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/get_MultiAssayExp.R");
        saveRDS(
            get_MultiAssayExp(study = "Jung", input_dir = paste0("{prefix}", "processed"), expr_with_counts_isoforms=TRUE), 
            "{prefix}{filename}"
        );
        '
        """

rule format_snv:
    input:
        S3.remote(prefix + "download/SNV.csv.gz"),
        S3.remote(prefix + "processed/cased_sequenced.csv")
    output:
        S3.remote(prefix + "processed/SNV.csv")
    shell:
        """
        Rscript scripts/Format_SNV.R \
        {prefix}download \
        {prefix}processed \
        """

rule format_expr:
    input:
        S3.remote(prefix + "download/expr_list.rds"),
        S3.remote(prefix + "processed/cased_sequenced.csv")
    output:
        S3.remote(prefix + "processed/EXPR_gene_tpm.csv"),
        S3.remote(prefix + "processed/EXPR_gene_counts.csv"),
        S3.remote(prefix + "processed/EXPR_isoform_tpm.csv"),
        S3.remote(prefix + "processed/EXPR_isoform_counts.csv"),
    shell:
        """
        Rscript scripts/Format_EXPR.R \
        {prefix}download \
        {prefix}processed \
        """

# rule format_cna_seg:
#     input:
#         S3.remote(prefix + "download/CNA_seg.txt.gz")
#     output:
#         S3.remote(prefix + "processed/CNA_seg.txt")
#     shell:
#         """
#         Rscript scripts/Format_CNA_seg.R \
#         {prefix}download \
#         {prefix}processed \
#         """

# rule format_cna_gene:
#     input:
#         S3.remote(prefix + "processed/cased_sequenced.csv"),
#         S3.remote(prefix + "download/all_thresholded.by_genes.txt.gz")
#     output:
#         S3.remote(prefix + "processed/CNA_gene.csv")
#     shell:
#         """
#         Rscript scripts/Format_CNA_gene.R \
#         {prefix}download \
#         {prefix}processed \
#         """

rule format_clin:
    input:
        S3.remote(prefix + "processed/cased_sequenced.csv"),
        S3.remote(prefix + "download/CLIN.txt"),
        S3.remote(prefix + "download/CLIN_Sample.txt")
    output:
        S3.remote(prefix + "processed/CLIN.csv")
    shell:
        """
        Rscript scripts/Format_CLIN.R \
        {prefix}download \
        {prefix}processed \
        """

rule format_cased_sequenced:
    input:
        S3.remote(prefix + "download/CLIN.txt"),
        S3.remote(prefix + "download/expr_list.rds"),
        S3.remote(prefix + "download/SNV.csv.gz"),
        S3.remote(prefix + "download/CNA_seg.txt.gz")
    output:
        S3.remote(prefix + "processed/cased_sequenced.csv")
    shell:
        """
        Rscript scripts/Format_cased_sequenced.R \
        {prefix}download \
        {prefix}processed \
        """

rule format_downloaded_data:
    input:
        S3.remote(prefix + "download/GSE135222_series_matrix.txt.gz"),
        S3.remote(prefix + "download/41467_2019_12159_MOESM1_ESM.pdf"),
        S3.remote(prefix + "download/rnaseq_samples.tsv"),
        S3.remote(prefix + "download/41467_2019_12159_MOESM8_ESM.xlsx"),
        S3.remote(prefix + "download/41467_2019_12159_MOESM9_ESM.xlsx"),
        S3.remote(prefix + "download/Jung_kallisto.zip"),
        S3.remote(prefix + "annotation/Gencode.v40.annotation.RData")
    output:
        S3.remote(prefix + "download/CLIN.txt"),
        S3.remote(prefix + "download/CLIN_Sample.txt"),
        S3.remote(prefix + "download/expr_list.rds"),
        S3.remote(prefix + "download/SNV.csv.gz"),
        S3.remote(prefix + "download/CNA_seg.txt.gz")        
    shell:
        """
        Rscript scripts/format_downloaded_data.R \
        {prefix}download \
        {prefix}annotation 
        """

rule download_annotation:
    output:
        S3.remote(prefix + "annotation/Gencode.v40.annotation.RData")
    shell:
        """
        wget https://github.com/BHKLAB-Pachyderm/Annotations/blob/master/Gencode.v40.annotation.RData?raw=true -O {prefix}annotation/Gencode.v40.annotation.RData 
        """

rule download_data:
    output:
        S3.remote(prefix + "download/GSE135222_series_matrix.txt.gz"),
        S3.remote(prefix + "download/41467_2019_12159_MOESM1_ESM.pdf"),
        S3.remote(prefix + "download/rnaseq_samples.tsv"),
        S3.remote(prefix + "download/41467_2019_12159_MOESM8_ESM.xlsx"),
        S3.remote(prefix + "download/41467_2019_12159_MOESM9_ESM.xlsx"),
        # S3.remote(prefix + "download/all_thresholded.by_genes.txt.gz"),
        S3.remote(prefix + "download/Jung_kallisto.zip")
    shell:
        """
        wget -O {prefix}download/GSE135222_series_matrix.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135222/matrix/GSE135222_series_matrix.txt.gz
        wget -O {prefix}download/41467_2019_12159_MOESM1_ESM.pdf https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-019-12159-9/MediaObjects/41467_2019_12159_MOESM1_ESM.pdf
        wget -O {prefix}download/rnaseq_samples.tsv "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA557841&result=read_run&fields=run_accession,sample_title&format=tsv&download=true&limit=0"
        wget -O {prefix}download/41467_2019_12159_MOESM8_ESM.xlsx https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-019-12159-9/MediaObjects/41467_2019_12159_MOESM8_ESM.xlsx
        wget -O {prefix}download/41467_2019_12159_MOESM9_ESM.xlsx https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-019-12159-9/MediaObjects/41467_2019_12159_MOESM9_ESM.xlsx
        wget -O {prefix}download/Jung_kallisto.zip https://zenodo.org/record/6968389/files/Jung_kallisto.zip?download=1
        """ 