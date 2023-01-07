library("biomaRt")
library(readr)

repository = "candida-infection-pipeline"

# For RStudio only
current_dir = dirname(rstudioapi::getActiveDocumentContext()$path)
dir = strsplit(current_dir, repository)[[1]][1]

matrix = read.csv(file.path(paste0(dir, repository, "/results/deseq_results.csv")), stringsAsFactors = FALSE)
genes = matrix$X

listMarts()
mart = useMart("ensembl")
datasets = listDatasets(mart)
ensembl = useDataset("mmusculus_gene_ensembl", mart = mart)
attributes = listAttributes(ensembl)

annotation = getBM(filters = "ensembl_gene_id",
                   attributes = c("ensembl_gene_id", "description", "external_gene_name",
                                  "mgi_id", "entrezgene_id", "gene_biotype"),
                   values = genes, mart = ensembl)
#table(complement$gene_biotype)

write.csv(as.data.frame(annotation), row.names = FALSE, file = paste0(dir, repository, "/results/gene_annotation.csv"))
