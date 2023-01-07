library(VennDiagram)
library(tidyverse)
library(stringr)

options(stringsAsFactors = FALSE)
options(digits = 15)

dir = "~/Documentos/Mestrado/candida_lungs/candida-infection-pipeline"

df = read.csv(file = paste0(dir, "/results/MMxGS_expression.csv"))
module = "turquoise"
#module = "blue"

tabela_ml = read.table(file.path(paste0(dir, "/results/27K-ranking.csv")), header = FALSE, sep = ",")
genes_ml = subset(tabela_ml, tabela_ml$V2 == 1)
mod = subset(df, df$moduleColor == module)
mod_ml = subset(mod, mod$ensembl_gene_id %in% genes_ml$V1)

color = c("#D9C230", "#6E764D", "#C63D31", "#3981A7")

genes = list(
  mod %>% filter(abs(mod$log2FoldChange) > 1 & mod$padj < 0.05) %>% select(ensembl_gene_id) %>% unlist(),
  mod %>% filter(gene_biotype == "lncRNA") %>% select(ensembl_gene_id) %>% unlist(),
  mod %>% filter(abs(GS.condition) > 0.6 & MM > 0.6) %>% select(ensembl_gene_id) %>% unlist(),
  mod_ml$ensembl_gene_id
)

#Make the plot
venn = venn.diagram(x = genes,
                    category.names = c("DEGs", "lncRNA", "MM > 0.6 \n GS > 0.6", "ML"),
                    main = "Candidate genes of interest",
                    sub = paste0(str_to_title(module), " module"),
                    fontfamily = "helvica",
                    main.fontfamily = "helvica",
                    sub.fontfamily = "helvica",
                    cat.fontfamily = "helvica",
                    main.cex = 0.45,
                    sub.cex = 0.4,
                    cat.cex = 0.35,
                    cex = 0.35,
                    output = FALSE,
                    lwd = 0.5,
                    lty = "dotted",
                    fill = color,
                    col = "#FFFFFF",
                    alpha = 0.45, 
                    height = 900,
                    width = 1000,
                    filename = paste0(dir, "/results/figures/venn_", module, ".tiff"),
                    imagetype = "tiff"
)
