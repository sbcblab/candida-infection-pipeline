library(readr)

dir = "~/Documentos/Mestrado/candida_lungs/candida-infection-pipeline"

expression = read.csv(file = paste0(dir, "/results/MMxGS_expression.csv"))

svc_results = read.table(file.path(paste0(dir, "/results/27K-ranking.csv")), header = FALSE, sep = ",")
genes_ml = subset(svc_results, svc_results$V2 == 1)

expr_ml = subset(expression, expression$ensembl_gene_id %in% genes_ml$V1)
expr_ml_deg = subset(expr_ml, abs(expr_ml$log2FoldChange) > 1 & expr_ml$padj < 0.05)

expr_ml_deg_mmgs_t = subset(expr_ml_deg, expr_ml_deg$MM > 0.6 & expr_ml_deg$GS < -0.6)
turquoise = subset(expr_ml_deg_mmgs_t, expr_ml_deg_mmgs_t$moduleColor == "turquoise")
turquoise_lncrna = subset(turquoise, turquoise$gene_biotype == "lncRNA")

expr_ml_deg_mmgs_b = subset(expr_ml_deg, expr_ml_deg$MM > 0.6 & expr_ml_deg$GS > 0.6)
blue = subset(expr_ml_deg_mmgs_b, expr_ml_deg_mmgs_b$moduleColor == "blue")
blue_lncrna = subset(blue, blue$gene_biotype == "lncRNA")

write.table(turquoise_lncrna, paste0(dir, "/results/lncRNA_turquoise.csv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

write.table(blue_lncrna, paste0(dir, "/results/lncRNA_blue.csv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
