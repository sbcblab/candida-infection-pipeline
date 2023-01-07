library(readr)
library(tibble)

dir = "~/Documentos/Mestrado/candida_lungs/candida-infection-pipeline"
cpc = read.csv(file = paste0(dir, "/results/cpc2_turquoise.tab"), sep = "\t")
exp = read.csv(file = paste0(dir, "/results/selected_genes_turquoise.csv"), sep = "\t")

gene = NULL
transcript = NULL
gene_name = NULL

for(x in 1:nrow(cpc)){
    temp = unlist(strsplit(cpc[x, "X.ID"], "_"))
    gene = append(gene, temp[1])
    transcript = append(transcript, temp[2])
    temp2 = subset(exp, exp$ensembl_gene_id == temp[1])
    gene_name = append(gene_name, temp2$external_gene_name)
}

cpc = add_column(cpc, gene, .after = 1)
cpc = add_column(cpc, transcript, .after = 2)
cpc = add_column(cpc, gene_name, .after = 3)
cpc$X.ID = NULL

coding = subset(cpc, cpc$label == "coding")
result = subset(cpc, cpc$gene %in% coding$gene)

result_exp = subset(exp, exp$ensembl_gene_id %in% result$gene)

