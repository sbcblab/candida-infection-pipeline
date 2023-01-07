library("biomaRt")
library(readr)

dir = "~/Documentos/Mestrado/candida_lungs/candida-infection-pipeline"
mod = "turquoise"
#df = read.csv(file = paste0(dir, "/results/selected_genes_", mod, ".csv"), sep = "\t")
lnames = load(paste0(dir, "/results/save/lncrna_connections_t.RData"))

genes = read.csv(file = paste0(dir, "/results/selected_genes_0.4_turquoise.csv"))
expressao = read.csv(file = paste0(dir, "/results/MMxGS_expression.csv"))
expressao = subset(expressao, expressao$ensembl_gene_id %in% genes$ensembl_gene_id)

listMarts()
mart = useMart("ensembl")
#datasets = listDatasets(mart)
ensembl = useDataset("mmusculus_gene_ensembl", mart = mart)
attributes = listAttributes(ensembl)

annotation = getBM(filters = "ensembl_gene_id",
                   attributes = c("ensembl_gene_id", "ensembl_transcript_id",
                                  "external_gene_name", "gene_exon_intron"), #cdna
                   values = names(df), mart = ensembl)

#annotation = annotation[order(annotation$ensembl_gene_id, decreasing = FALSE), ]

for(lncrna in names(df)){
    df_lncrna = subset(annotation, annotation$ensembl_gene_id == lncrna)
    fasta = NULL
    
    for(x in 1:nrow(df_lncrna)){
        fasta = paste0(fasta, '>', df_lncrna[x, "ensembl_gene_id"], " ",
                       df_lncrna[x, "ensembl_transcript_id"], " ",
                       df_lncrna[x, "external_gene_name"], "\n",
                       df_lncrna[x, "gene_exon_intron"], "\n\n")
    }
    
    arquivo = file(paste0(dir, "/results/intarna/lncrna/", lncrna, ".fasta"))
    write_file(fasta, arquivo)
}

rm(fasta, lncrna, df_lncrna, arquivo)

id_usados = list()

for(lncrna in names(df)[!names(df) %in% id_usados]){
    genes_temp = df[[x]]$ensembl_gene_id
    annotation = getBM(filters = "ensembl_gene_id",
                       attributes = c("ensembl_gene_id", "ensembl_transcript_id",
                                      "external_gene_name", "cdna"),
                       values = genes_temp, mart = ensembl)
    
    fasta = NULL
    
    for(y in 1:nrow(annotation)){
        fasta = paste0(fasta, '>', annotation[y, "ensembl_gene_id"], "_",
                       annotation[y, "ensembl_transcript_id"], " ",
                       annotation[y, "external_gene_name"], "\n", annotation[y, "cdna"], "\n")
    }
    
    print(dim(annotation))
    arquivo = file(paste0(dir, "/results/intarna/protein_coding/", lncrna, "_pc_genes.fasta"))
    write_file(fasta, file = arquivo)
    
    id_usados = c(id_usados, lncrna)
}
