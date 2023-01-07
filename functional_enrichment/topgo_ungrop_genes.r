library(readr)

dir = "~/Documentos/Mestrado/candida_lungs/candida-infection-pipeline"

mod = "turquoise"
#mod = "blue"
#topgo = read.csv(file = paste0(dir, "/results/topgo_", mod, ".tab"), sep = "\t")
#genes_selecionados = read.csv(file = paste0(dir, "/results/selected_genes_", mod, ".csv"), sep = "\t")
topgo = read.csv(file = paste0(dir, "/results/topgo_selected_genes_final.tab"), sep = "\t")
rede_mod = read.csv(file = paste0(dir, "/results/MMxGS_expression.csv"))

goid = NULL
gene = NULL

for(x in 1:nrow(topgo)){
  gene_temp = unlist(strsplit(topgo[x, "Genes"], ","))
  goid_temp = replicate(length(gene_temp), topgo[x, "ID"])
  goid = append(goid, goid_temp)
  gene = append(gene, gene_temp)
}

df = data.frame(GENE = gene, GOID = goid)

lncrna = subset(rede_mod, rede_mod$gene_biotype == 'lncRNA')
lncrna_topgo = subset(df, df$GENE %in% lncrna$ensembl_gene_id)

topgo_lncrna = subset(topgo, topgo$ID %in% lncrna_topgo$GOID)

teste = subset(rede_mod, rede_mod$ensembl_gene_id %in% lncrna_topgo$GENE)

#write.table(df, paste0(dir, "/rascunho/topgo_selected_genes_", mod, "_ungrouped.tab"), quote = FALSE, sep = "\t", row.names = FALSE)
