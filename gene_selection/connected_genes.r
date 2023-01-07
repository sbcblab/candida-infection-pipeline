library(tibble)

options(stringsAsFactors = FALSE)

dir = "~/Documentos/Mestrado/candida_lungs/candida-infection-pipeline"

mod = "turquoise"
network = read.csv(paste0(dir, "/results/network_selected_genes_", mod, ".tab"),
                header = TRUE, sep = "\t")

# Removing genes incorrectly classified as lncRNA
remove = c("ENSMUSG00000108695", "ENSMUSG00000085564", "ENSMUSG00000086480",
            "ENSMUSG00000103041", "ENSMUSG00000105703", "ENSMUSG00000091542",
            "ENSMUSG00000053889")

genes = read.csv(file = paste0(dir, "/results/selected_genes_0.4_", mod, ".csv"))
genes = subset(genes, !genes$ensembl_gene_id %in% remove)

expr_total = read.csv(file = paste0(dir, "/results/MMxGS_expression.csv"))
expression = subset(expr_total, expr_total$ensembl_gene_id %in% genes$ensembl_gene_id)
lncrna = subset(expression, expression$gene_biotype == "lncRNA")
protein_coding = subset(expression, expression$gene_biotype == "protein_coding")

# Filting the network to include only lncnRNA-protein-coding interactions
network_2 = subset(network, (network$fromNode %in% lncrna$ensembl_gene_id & 
                           network$toNode %in% protein_coding$ensembl_gene_id) |
                    (network$fromNode %in% protein_coding$ensembl_gene_id & 
                         network$toNode %in% lncrna$ensembl_gene_id))

genes_total = c(unlist(network_2$fromNode), unlist(network_2$toNode))
freq = data.frame(table(genes_total))
genes_total = unlist(unique(genes_total))

# LncRNAs interacting with protein coding genes
lncrna_total = genes_total[genes_total %in% lncrna$ensembl_gene_id]

freq_lncrna = subset(freq, freq$genes_total %in% lncrna$ensembl_gene_id)
freq_pc = subset(freq, freq$genes_total %in% protein_coding$ensembl_gene_id)


topgo = read.csv(file = paste0(dir, "/results/topgo_selected_genes_final.tab"), sep = "\t")

goid = NULL
gene = NULL

# Ungrouping the genes from the topGO results

for(x in 1:nrow(topgo)){
    gene_temp = unlist(strsplit(topgo[x, "Genes"], ","))
    goid_temp = replicate(length(gene_temp), topgo[x, "ID"])
    goid = append(goid, goid_temp)
    gene = append(gene, gene_temp)
}

topgo_ungroup = data.frame(GENE = gene, GOID = goid)
rm(goid, gene, gene_temp, goid_temp, x)

freq_pc = subset(freq_pc, freq_pc$genes_total %in% topgo_ungroup$GENE)

# Filtering the network only keep interactions from genes included in the topGO results
network_3 = subset(network_2, network_2$fromNode %in% topgo_ungroup$GENE |
                    network_2$toNode %in% topgo_ungroup$GENE)

#write.table(network_3, paste0(dir, "/results/network_selected_genes_topgo_", mod, ".tab"),
#          sep = "\t", quote = FALSE, row.names = FALSE)

genes_total_3 = c(unlist(network_3$fromNode), unlist(network_3$toNode))
freq_3 = data.frame(table(genes_total_3))
genes_total_3 = unlist(unique(genes_total_3))

freq_lncrna_3 = subset(freq_3, freq_3$genes_total_3 %in% lncrna$ensembl_gene_id)
freq_pc_3 = subset(freq_3, freq_3$genes_total_3 %in% protein_coding$ensembl_gene_id)


#write.table(freq_pc, paste0(dir, "/results/freq_protein_coding_genes_t.tab"),
#            row.names = FALSE, sep = "\t" , quote = FALSE)

up = NULL
down = NULL

# Quantifying how many up and down regulated genes each BP was enriched for

for(x in topgo$Genes){
    temp = unlist(strsplit(x, ","))
    up_temp = subset(expr_total, expr_total$ensembl_gene_id %in% temp &
                         expr_total$moduleColor == "turquoise")
    up = c(up, dim(up_temp)[1])
    
    down_temp = subset(expr_total, expr_total$ensembl_gene_id %in% temp &
                           expr_total$moduleColor == "blue")
    down = c(down, dim(down_temp)[1])
}

# Adding the coutings to the topGO data frame
topgo = add_column(topgo, up, .after = 3)
topgo = add_column(topgo, down, .after = 4)

# Lists of data frames - containing all the PC genes each lncRNA interacts with
# in the network_3 (df) and all BPs enriched for each PC gene (bp_lncrna)
df = list()
bp_lncrna = list()

for(x in freq_lncrna_3$genes_total_3){
    count = NULL
    temp = subset(network_3, network_3$fromNode == x | network_3$toNode == x)
    pc_genes = c(unlist(temp$fromNode), unlist(temp$toNode))
    pc_genes = unlist(unique(pc_genes))
    pc_genes = pc_genes[pc_genes != x]
    
    processes = subset(topgo_ungroup, topgo_ungroup$GENE %in% pc_genes)
    temp2 = data.frame(table(processes$GOID))
    df[[x]] = subset(expression, expression$ensembl_gene_id %in% pc_genes)
    bp_lncrna[[x]] = merge(temp2, topgo, by.x = "Var1", by.y = "ID")
}

rm(temp, temp2, pc_genes, processes, x, count)

#save(df, bp_lncrna, topgo, network_3,
#     file = paste0(dir, "/results/save/lncrna_connections_", mod, ".RData"))

#write.csv(lncrna$ensembl_gene_id, paste0(dir, "/results/selected_lncrna_", mod, ".txt"),
#          quote = FALSE, row.names = FALSE)

#write.csv(network_3, paste0(dir, "/results/network_selected_genes_", mod, "_final.tab"),
#          quote = FALSE, row.names = FALSE)
