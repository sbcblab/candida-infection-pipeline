library(tximport)
library(GenomicFeatures)

repository = "candida-infection-pipeline"

# For RStudio only
current_dir = dirname(rstudioapi::getActiveDocumentContext()$path)
dir = strsplit(current_dir, repository)[[1]][1]

gtf_path = paste0(dir, "/reference/Mus_musculus.GRCm39.103.gtf")

transcript_db = makeTxDbFromGFF(file = gtf_path, format = "gtf",
                        dataSource = "ftp://ftp.ensembl.org/pub/release-103/gtf/mus_musculus/Mus_musculus.GRCm39.103.gtf.gz",
                        organism = "Mus musculus")

keys = keys(transcript_db, keytype = "TXNAME")
gene_id = select(transcript_db, keys, "GENEID", "TXNAME")

write.csv(gene_id, paste0(dir, repository, "/results/transcript_to_gene.csv"), row.names = FALSE)


