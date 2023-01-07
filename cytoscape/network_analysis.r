dir = "~/Documentos/Mestrado/candida_lungs/candida-infection-pipeline"
setwd(dir)

options(stringsAsFactors = FALSE)

df = read.csv(file = paste0(dir, "/results/MMxGS_expression.csv"))

module = "turquoise"
#rede = read.table(paste0(dir, "/results/", module, "_rede_820genes_filt03.txt"), sep = "\t", header = TRUE)
rede = read.table(paste0(dir, "/results/CytoscapeInput-edges-", module, "-lncrna-filtered.txt"), sep = "\t", header = TRUE)
rede = subset(rede, rede$weight > 0.3)

rede_genes = unique(c(rede$fromNode, rede$toNode))
genes = subset(df, df$ensembl_gene_id %in% rede_genes)

lncrna = subset(genes, genes$gene_biotype == "lncRNA")
n_lncrna = subset(genes, genes$gene_biotype != "lncRNA")

interacao_lncrna = subset(rede, rede$fromNode %in% lncrna$ensembl_gene_id &
                                  rede$toNode %in% lncrna$ensembl_gene_id)

interacao_protein_coding = subset(rede, rede$fromNode %in% n_lncrna$ensembl_gene_id &
                            rede$toNode %in% n_lncrna$ensembl_gene_id)

interacao_lncrna_pc = subset(rede, (rede$fromNode %in% lncrna$ensembl_gene_id & 
                                            rede$toNode %in% n_lncrna$ensembl_gene_id) |
                                          (rede$fromNode %in% n_lncrna$ensembl_gene_id & 
                                             rede$toNode %in% lncrna$ensembl_gene_id))

interacao = data.frame(interacao = c("lncRNA - lncRNA", "lncRNA - Protein coding", "Protein coding - protein coding"),
                       quantidade = c(dim(interacao_lncrna)[1], dim(interacao_lncrna_pc)[1], dim(interacao_protein_coding)[1]),
                       min_weight = c(min(interacao_lncrna$weight), min(interacao_lncrna_pc$weight), min(interacao_protein_coding$weight)),
                       max_weight = c(max(interacao_lncrna$weight), max(interacao_lncrna_pc$weight), max(interacao_protein_coding$weight)),
                       mean_weight = c(mean(interacao_lncrna$weight), mean(interacao_lncrna_pc$weight), mean(interacao_protein_coding$weight)))

# Todos genes estão over expressed
