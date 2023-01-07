library(tidyr)

#dir = "~/workspace/candida_lungs/candida-infection-pipeline"
dir = "~/Documentos/Mestrado/candida_lungs/candida-infection-pipeline"

mod = "blue"
rede = read.csv(file = paste0(dir, "/results/edges_", mod, ".csv"))
genes = read.csv(file = paste0(dir, "/results/selected_genes_0.4_", mod, ".csv"))

expressao = read.csv(file = paste0(dir, "/results/MMxGS_expression.csv"))
expressao = subset(expressao, expressao$ensembl_gene_id %in% genes$ensembl_gene_id)
lncrna = subset(expressao, expressao$gene_biotype == "lncRNA")

rede_original = read.table(file = paste0(dir, "/results/CytoscapeInput-edges-", mod, ".txt"), 
                         header = TRUE, col.names = c("fromNode", "toNode", "weight", "direction",
                                                      "fromAltName", "toAltName"),
                         sep = "\t")

df = subset(rede_original, rede_original$fromNode %in% genes$ensembl_gene_id &
              rede_original$toNode %in% genes$ensembl_gene_id & rede_original$weight > 0.4)

df = subset(df, df$fromNode %in% lncrna$ensembl_gene_id | df$toNode %in% lncrna$ensembl_gene_id)

table(table(df$weight))

rede_nova = merge(rede, df, by = "weight", all.x = TRUE)

write.table(rede_nova, file = paste0(dir, "/results/network_selected_genes_", mod, ".tab"),
            row.names = FALSE, sep = "\t" , quote = FALSE)

from_biotype = subset(expressao, expressao$ensembl_gene_id %in% rede_nova$fromNode)
to_biotype = subset(expressao, expressao$ensembl_gene_id %in% rede_nova$toNode)

